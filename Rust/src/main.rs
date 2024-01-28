#![allow(warnings, unused)]

use std::{time::Instant, cmp::min};
use collapsed_matrices::CollapsedMatrix;

mod collapsed_matrices{
    use std::ops::{Add, Index, IndexMut};

    use factorial::Factorial;
    use rayon::prelude::*;

    #[derive(Debug)]
    pub struct CollapsedMatrix <T> 
        where T: Send+Clone+Sync+Add<Output = T>
    {
        pub dimension: u8,
        pub n: u16,
        pub length: u32,
        pub indices: Vec<Vec<usize>>,
        pub results: Vec<T>,
        pub sum: Vec<T>
    }
    
    impl <T> CollapsedMatrix <T> 
        where T: Send + Clone + Sync + Add<Output = T>
    {
        pub fn new (dimension: u8, n: u16) -> CollapsedMatrix <T> {
            let length = nCr(n as u32, dimension as u32);
            let indices: Vec<Vec<usize>> = set_indices(dimension, length);
            let results:Vec<T> = Vec::with_capacity(length as usize);
            let sum: Vec<T> = Vec::with_capacity(n as usize);

            CollapsedMatrix {
                dimension,
                n,
                length,
                indices,
                results,
                sum
            }
        }

        fn calculate_sum (&self) -> Vec<T> {
            let sum : Vec<_> = (0..self.n).into_par_iter().map(|x| {
                let mut output:T = self.sum[0].clone();
                let mut i = 0usize;
                for index in &self.indices {
                    for item in index {
                        if *item == x as usize {
                            output = output + self.results[i].clone();
                            break;
                        }
                    }
                    i+=1;
                }
                output
            }).collect();
            sum
        }
    }

    impl <T> Index<Vec<usize>> for CollapsedMatrix <T> 
        where T: Send + Clone + Sync + Add<Output = T>
    {
        type Output = T;
        fn index(&self, mut indices: Vec<usize>) -> &T {
            let mut t_index = 0_usize;
            indices.sort(); indices.reverse();
            for i in 0..indices.len() {t_index += nCr(indices[i] as u32, (indices.len()-i) as u32) as usize;}
            &self.results[t_index]
        }
    }
    impl <T> IndexMut<Vec<usize>> for CollapsedMatrix <T> 
        where T: Send + Clone + Sync + Add<Output = T>
    {
        fn index_mut(&mut self, mut indices: Vec<usize>) -> &mut Self::Output {
            let mut t_index = 0_usize;
            indices.sort(); indices.reverse();
            for i in 0..indices.len() {t_index += nCr(indices[i] as u32, (indices.len()-i) as u32) as usize;}
            &mut self.results[t_index]
        }
    }

    fn set_indices(dimension: u8, length: u32) -> Vec<Vec<usize>> {
        let mut indices: Vec<Vec<usize>> = vec![vec![0;dimension as usize];length as usize];
        let mut item = vec![0usize; dimension as usize];
        for i in 0..dimension as usize {item[i] = dimension as usize - 1 - i;}
        indices[0] = item.clone();
        for i in 1..length as usize{
            for j in (0..dimension as usize).rev(){
                if j == 0usize {
                    item[0] += 1; 
                    for k in (j+1)..dimension as usize {
                        item[k] = dimension as usize - 1 - k;
                    }
                }
                else if item[j]+1 != item[j-1] {
                    item[j] += 1usize; 
                    for k in (j+1)..dimension as usize {
                        item[k] = dimension as usize - 1 - k;
                    }
                    break;
                }
            }
            indices[i] = item.clone();
        }
        indices
    }
    
    fn nCr (n: u32, r: u32) -> u32 {
        if n<r {0u32}
        else{
            ((n+1-r)..(n+1)).product::<u32>()/(r.factorial())
        }
    }
}

use euclid::*;
use rand::{thread_rng, Rng};
use rayon::{prelude::{ParallelIterator, IntoParallelIterator}, slice::ParallelSlice};

#[derive(Debug)]
struct Body {
    position : default::Vector2D<f32>,
    mass : f32,
    velocity : default::Vector2D<f32>,
    acceleration : default::Vector2D<f32>
}

impl Body {
    fn new_rand() -> Body {
        let mut rng = thread_rng();
        Body {
            position : Vector2D::new(rng.gen_range(0.0..100.0), rng.gen_range(0.0..100.0)),
            mass : rng.gen_range(0.0..100.0),
            velocity : Vector2D::new(0.0, 0.0),
            acceleration : Vector2D::new(0.0, 0.0)
        }
    } 
}

fn grav_force (indices : &Vec<usize>, bodies:&Vec<Body>) -> default::Vector2D<f32> {
    let GRAV_CONST: f32 = 5.0;
    let vec12 = bodies[indices[0]].position-bodies[indices[1]].position;
    let force_mag = -GRAV_CONST * (bodies[indices[0]].mass*bodies[indices[1]].mass)/vec12.square_length();
    let force_angle = vec12.angle_from_x_axis();
    Vector2D::from_angle_and_length(force_angle, force_mag)
}

fn par_calculate_grav_froces (grav_forces: &CollapsedMatrix<default::Vector2D<f32>>, objects : &Vec<Body>) -> Vec<default::Vector2D<f32>>{
    // let n_cpus = min(num_cpus::get(), objects.len());
    let n_cpus = min(32_usize, objects.len());
    let chunks = grav_forces.indices.par_chunks(grav_forces.length as usize/(n_cpus) +1);

    let results:Vec<_> = chunks.flat_map(|chunk|{
        let mut ouput: Vec<default::Vector2D<f32>> = Vec::with_capacity(chunk.len());
        for indices in chunk {
            ouput.push(grav_force(indices, objects));
        }
        ouput
    }).collect();/* 
    let results:Vec<_> = grav_forces.indices.clone().into_par_iter().map(|x| {
        grav_force(&x, objects)
    }).collect(); */
    results
}

fn calculate_grav_froces (grav_forces: &CollapsedMatrix<default::Vector2D<f32>>, objects : &Vec<Body>) -> Vec<default::Vector2D<f32>> {
    let mut results: Vec<default::Vector2D<f32>> = Vec::with_capacity(grav_forces.length as usize);
    for i in 0..grav_forces.length as usize {
        results.push(grav_force(&grav_forces.indices[i], objects));
    }
    results
}

fn calculate_grav_sum (grav_forces: &CollapsedMatrix<default::Vector2D<f32>>, signs: Vec<f32>) -> Vec<default::Vector2D<f32>> {
    let sum : Vec<_> = (0..grav_forces.n).into_par_iter().map(|x| {
        let mut output:default::Vector2D<f32> = Vector2D::new(0.0, 0.0);
        let mut i = 0usize;
        for index in &grav_forces.indices {
            let mut j = 0usize;
            for item in index {
                if *item == x as usize {
                    output = output + grav_forces.results[i].clone() * signs[j];
                    break;
                }
                j+=1;
            }
            i+=1;
        }
        output
    }).collect();
    sum
}

fn main() {
    let n_bodies : u16 = 2;
    let n_times : usize = 1;
    let mut objects = Vec::<Body>::with_capacity(n_bodies as usize);
    
    println!("n cores:\t{}", num_cpus::get());
    println!("n_bodies\t-\ttime(us)");
    for i in 1..14_u32{
        let n_bodies:u16 = 2_u16.pow(i);
        let mut objects = Vec::<Body>::with_capacity(n_bodies as usize);
        for _ in 0..n_bodies {
            objects.push(Body::new_rand());
        }
        
        let mut grav_forces = CollapsedMatrix::<default::Vector2D<f32>>::new(2, n_bodies);

        let now = Instant::now();
        for _ in 0..n_times {
            grav_forces.results = par_calculate_grav_froces(&grav_forces, &objects);
            //grav_forces.sum = calculate_grav_sum(&grav_forces, vec![-1.0, 1.0]);
        }
        let time_1 = now.elapsed();
        println!("{}\t\t-\t{}", 2_u16.pow(i), time_1.as_micros())
    }
    //######################################################################################################################################################
    for i in 1..14 {   
        let n_bodies:u16 = 2_u16.pow(i);
        let mut objects = Vec::<Body>::with_capacity(n_bodies as usize);
        for _ in 0..n_bodies {
            objects.push(Body::new_rand());
        }
        let mut grav_forces = CollapsedMatrix::<default::Vector2D<f32>>::new(2, n_bodies);

        let now = Instant::now();
        for _ in 0..n_times {
            grav_forces.results = calculate_grav_froces(&grav_forces, &objects);
            //grav_forces.sum = calculate_grav_sum(&grav_forces, vec![-1.0, 1.0]);
        }
        let time_2 = now.elapsed();
        println!("{}\t\t-\t{}", 2_u16.pow(i), time_2.as_micros())
    }
    //######################################################################################################################################################
    // for i in 1..15 {
    //     let n_bodies:u16 = 2_u16.pow(i);
    //     let mut objects = Vec::<Body>::with_capacity(n_bodies as usize);
    //     for _ in 0..n_bodies {
    //         objects.push(Body::new_rand());
    //     }
    //     let now = Instant::now();
    //     for _ in 0..n_times{
    //         for i in 1..n_bodies as usize {
    //             for j in 0..i{
    //                 objects[i].acceleration = grav_force(&vec![i,j], &objects);
    //             }
    //         }
    //     }
    //     let time_3 = now.elapsed();
    //     println!("{}\t\t-\t{}", 2_u16.pow(i), time_3.as_micros())
    // }

    // println!("{}", time_1.as_millis());
    /* println!("{}", time_2.as_micros());
    println!("{}", time_3.as_micros()); */
    // println!("{:?}", objects);
    // println!("{:?}", grav_forces.indices);
    // println!("{:?}", grav_forces.results);
    // println!("{:?}", grav_forces[vec![3,1]]);
    // println!("{:?}", grav_forces.sum); 
}