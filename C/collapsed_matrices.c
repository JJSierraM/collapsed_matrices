#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include "vector2D.h"
#include "body.h"

typedef unsigned int uint;

inline uint nCr (uint n, uint r){
    uint output = 1;
    if (n<r) {output = 0;}
    else {
        uint i;
        for (i=1; i<r+1; i++) {output=(output*(i+n-r))/i;}
    }
    return output;
}

uint32_t xorshift32(uint32_t* seed)
{
    // Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs"
    // See <https://stackoverflow.com/questions/53886131/how-does-xorshift32-works>
    // https://en.wikipedia.org/wiki/Xorshift
	uint32_t x = *seed;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	return *seed = x;
}

inline float rand_0_to_1(uint32_t* seed){
	return ((float) xorshift32(seed)) / ((float) UINT32_MAX);
}

inline float rand_float(uint32_t* seed, float from, float to)
{
    return rand_0_to_1(seed) * (to-from) + from;
}

typedef struct CollapsedMatrix 
{
    uint dimension;         // Number of items the relations involves   (max 65,535)
    uint n;                 // Number of items to relate                (max 65,535) 
    uint length;            // Number of total relations                (max 65,535)

    uint* indices;          // Vector of indices, contiguous. Length=dimension*length
    Vector2D* results;      // Vector of resulted calculations
    Vector2D* sum;          // Vector of sum of the results by index 
} CollapsedMatrix;

inline void transfer_array(uint* arr1, uint* arr2, uint dim_arr2){
    for (int i=0; i<dim_arr2; i++){
        *(arr1+i)=*(arr2+i);
    }
}

void set_indices(CollapsedMatrix* collapsed_matrix) {
    uint dimension = collapsed_matrix->dimension;
    uint length = collapsed_matrix->length;

    uint *item = malloc(dimension*sizeof(uint));
    uint i,j,k;
    for (i=0; i<dimension; i++) item[i] = dimension -1 -i;
    transfer_array(collapsed_matrix->indices+0, item, dimension);
    for (i=1; i<length; i++){
        for (j=dimension-1; j>=0; j--){
            if (j == 0)
            {
                item[0]++;
                for (k=j+1; k<dimension; k++){
                    item[k] = dimension -1 -k;
                }
                break;
            }
            else if (item[j]+1 != item[j-1])
            {
                item[j]++;
                for (k=j+1; k<dimension; k++){
                    item[k] = dimension -1 -k;
                }
                break;
            }
        }
        transfer_array(collapsed_matrix->indices+(i*dimension), item, dimension);
    }
    free(item);
}

void initialize_results(CollapsedMatrix* collapsed_matrix) {
    for (int i=0; i<collapsed_matrix->length; i++) vector2D_initialize(collapsed_matrix->results+i);
}

void initialize_sum(CollapsedMatrix* collapsed_matrix) {
    for (int i=0; i<collapsed_matrix->n; i++) vector2D_initialize(collapsed_matrix->sum+i);
}

CollapsedMatrix* collapsed_matrix_new(uint dimension_, uint n_) {
    CollapsedMatrix* collapsed_matrix = malloc(sizeof(CollapsedMatrix));
    collapsed_matrix->dimension = dimension_;
    collapsed_matrix->n = n_;

    collapsed_matrix->length = nCr(n_, dimension_);
    collapsed_matrix->indices = malloc(collapsed_matrix->length*dimension_*sizeof(uint));
    collapsed_matrix->results = malloc(collapsed_matrix->length*sizeof(Vector2D));
    collapsed_matrix->sum = malloc(collapsed_matrix->n*sizeof(Vector2D));

    set_indices(collapsed_matrix);

    initialize_results(collapsed_matrix);
    initialize_sum(collapsed_matrix);

    return collapsed_matrix;
}

void collapsed_matrix_destroy(CollapsedMatrix* collapsed_matrix) {
    free(collapsed_matrix->indices);
    free(collapsed_matrix->results);
    free(collapsed_matrix->sum);
    free(collapsed_matrix);
}

static Vector2D grav_force(const Body* bodies) {
    const float GRAV_CONST = 5.0;
    const Vector2D vec12 = vector2D_subtract(bodies[0].position, bodies[1].position);
    const float force_mag = -GRAV_CONST * (bodies[0].mass*bodies[1].mass)/vector2D_sqr_length(vec12);
    const float force_angle = vector2D_angle(vec12);
    return vector2D_from_magnitude_and_angle(force_mag, force_angle);
}

static Vector2D grav_force_2(const Body* body_A, const Body* body_B) {
    const float GRAV_CONST = 5.0;
    const Vector2D vec12 = vector2D_subtract(body_A->position, body_B->position);
    const float force_mag = -GRAV_CONST * (body_A->mass*body_B->mass)/vector2D_sqr_length(vec12);
    const float force_angle = vector2D_angle(vec12);
    return vector2D_from_magnitude_and_angle(force_mag, force_angle);
}

void collapsed_matrix_apply_function(CollapsedMatrix* collapsed_matrix, const Body* items) {
    uint i;
    
    #pragma omp parallel for
    for (i=0; i<collapsed_matrix->length*collapsed_matrix->dimension; i+=collapsed_matrix->dimension){
        collapsed_matrix->results[i/collapsed_matrix->dimension] = grav_force_2(items+collapsed_matrix->indices[i], items+collapsed_matrix->indices[i+1]);
    }
}

void collapsed_matrix_calculate_sum(CollapsedMatrix* collapsed_matrix) {
    uint i;
    
    #pragma omp parallel for
    for (i = 0; i < (collapsed_matrix->length*collapsed_matrix->dimension); i++){
        vector2D_add_equals(&collapsed_matrix->sum[collapsed_matrix->indices[i]], &collapsed_matrix->results[i/collapsed_matrix->dimension]);
    }
}

int main() {
    uint n=1024, d=2;
    srand(1);

    clock_t start_fun, end_fun, start_sum, end_sum;

    Body* bodies = malloc(n*sizeof(Body));
    for (int i=0; i<n; i++) {body_random(bodies+i);}

    CollapsedMatrix* collapsed_matrix = collapsed_matrix_new(d, n);

    start_fun = clock();
    collapsed_matrix_apply_function(collapsed_matrix, bodies);
    // for (int i=0; i<nCr(n,d); i++) {grav_force(bodies);}
    end_fun = clock();

    start_sum = clock();
    collapsed_matrix_calculate_sum(collapsed_matrix);
    end_sum = clock();

    // start_fun = clock();
    // Vector2D force;
    // for (int i = 1; i < n; i++)
    // {
    //     for (int j = 0; j < i; j++)
    //     {
    //         force = grav_force_2(bodies+i, bodies+j);
    //         vector2D_add_equals(collapsed_matrix->results, &force);
    //     }
    // }
    // end_fun = clock();

    printf("\nTime spent:\n");
    printf("Apply function: %3.3f ms\n", ((double) (end_fun-start_fun)) / CLOCKS_PER_SEC * 1000);
    printf("Apply sumation: %3.3f ms\n", ((double) (end_sum-start_sum)) / CLOCKS_PER_SEC * 1000);
    
    return 0;
}