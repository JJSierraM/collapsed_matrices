#include<stdio.h>
#include<stdlib.h>
#include <time.h>
#include <math.h>

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

typedef struct Body {
    float mass;
    float2 position;
    float2 speed;
    float2 accel; 
} Body;

void body_initialize(Body* body) {
    body->mass=0.001;
    body->position = float2();
    body->speed = float2();
    body->accel = float2();
}

//Random mass and position, speed and accel set to 0
void body_random(Body* body) {
    const float rand_mass = ((float)rand()/(float)RAND_MAX)*(1.0-0.001)+0.001; 
    float2 rand_pos;
    vector2D_random(&rand_pos, -1.0, 1.0); 
    body->mass = rand_mass;
    body->position = rand_pos;
    vector2D_initialize(&body->speed);
    vector2D_initialize(&body->accel);
}

Body* body_new() {
    Body* output = (Body*) malloc(sizeof(Body));
    body_initialize(output);
    return output;
}

void body_print(Body* body) {
    printf("Mass:\t\t%06.3f\n", body->mass);
    printf("Position:\tx: %06.3f, y: %06.3f\n", body->position.x, body->position.y);
    printf("Speed:\t\tx: %06.3f, y: %06.3f\n", body->speed.x, body->speed.y);
    printf("Acceleration:\tx: %06.3f, y: %06.3f\n", body->accel.x, body->accel.y);
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

    uint *item = (uint*) malloc(dimension*sizeof(uint));
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
    CollapsedMatrix* collapsed_matrix = (CollapsedMatrix*) malloc(sizeof(CollapsedMatrix));
    collapsed_matrix->dimension = dimension_;
    collapsed_matrix->n = n_;

    collapsed_matrix->length = nCr(n_, dimension_);
    collapsed_matrix->indices = (uint*) malloc(collapsed_matrix->length*dimension_*sizeof(uint));
    collapsed_matrix->results = (Vector2D*) malloc(collapsed_matrix->length*sizeof(Vector2D));
    collapsed_matrix->sum = (Vector2D*) malloc(collapsed_matrix->n*sizeof(Vector2D));

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
    #pragma omp parallel private (i) 
    {
        #pragma omp for
        for (i=0; i<collapsed_matrix->length*collapsed_matrix->dimension; i+=collapsed_matrix->dimension){
            collapsed_matrix->results[i/collapsed_matrix->dimension] = grav_force_2(items+collapsed_matrix->indices[i], items+collapsed_matrix->indices[i+1]);
        }
    }
}

void collapsed_matrix_calculate_sum(CollapsedMatrix* collapsed_matrix) {
    uint i;
    # pragma omp parallel private ( i )
    {
        # pragma omp for
        for (i = 0; i < (collapsed_matrix->length*collapsed_matrix->dimension); i++){
            vector2D_add_equals(&collapsed_matrix->sum[collapsed_matrix->indices[i]], &collapsed_matrix->results[i/collapsed_matrix->dimension]);
        }
    }
}

struct cuBody {
    float mass;
    float2 position;
    float2 speed;
    float2 accel;
};

static __device__ __inline__ float float2_sqr_length(float2 vec) {
    return vec.x*vec.x+vec.y*vec.y;
}

static __device__ __inline__ float float2_angle(float2 vector) {
    return atan2(vector.y, vector.x);
}

static __device__ __inline__ float2 float2_from_magnitude_and_angle(float magnitude, float angle) {
    float2 output;
    output.x = magnitude*cos(angle);
    output.y = magnitude*sin(angle);
    return output;
}

static __device__ __inline__ float2 float2_subtract(float2 A, float2 B) {
    return make_float2(A.x-B.x, A.y-B.y);
}

static __device__ __inline__ void float2_add_equals(float2* a, const float2* b) {
    a->x = a->x + b->x;
    a->y = a->y + b->y;
}

struct cuCollapsedMatrix 
{
    uint dimension;         // Number of items the relations involves   (max 65,535)
    uint n;                 // Number of items to relate                (max 65,535) 
    uint length;            // Number of total relations                (max 65,535)

    uint* indices;          // Vector of indices, contiguous. Length=dimension*length
    float2* results;        // Vector of resulted calculations
    float2* sum;            // Vector of sum of the results by index 
};

__device__ uint cu_nCr(uint n, uint r) {
    uint output = 1;
    if (n<r) {output = 0;}
    else {
        uint i;
        for (i=1; i<r+1; i++) {output=(output*(i+n-r))/i;}
    }
    return output;
}

//grid(1,1,1)
__global__ void cu_collapsed_matrix_new(cuCollapsedMatrix* cu_collapsed_matrix, uint dimension_, uint n_) {
    cu_collapsed_matrix = (cuCollapsedMatrix*) malloc(sizeof(cuCollapsedMatrix));
    cu_collapsed_matrix->dimension = dimension_;
    cu_collapsed_matrix->n = n_;

    cu_collapsed_matrix->length = cu_nCr(n_, dimension_);
    cu_collapsed_matrix->indices = (uint*) malloc(cu_collapsed_matrix->length*dimension_*sizeof(uint));
    cu_collapsed_matrix->results = (float2*) malloc(cu_collapsed_matrix->length*sizeof(float2));
    cu_collapsed_matrix->sum = (float2*) malloc(cu_collapsed_matrix->n*sizeof(float2));
}

__device__ void cu_set_indices(uint* indices, uint dimension, int x) {
    uint prev = 0;
    uint* item = indices+(x*dimension);
    for (uint k=0; k<dimension; k++){
        for (uint j=0; j<999; j++){
            if ((x-prev) < cu_nCr(dimension-k+j, j)) {
                item[k] = j + dimension -k -1;
                prev += cu_nCr(dimension -k +j -1, j-1);
                break;
            }
        }
    }
}

__device__ void cu_initialize_results(uint length, float2* results, int x) {
    if (x<length) {
    results[x] = float2();
    }
}

__device__ void cu_initialize_sum(uint n, float2* sum, int x) {
    if (x<n) {
    sum[x] = float2();
    }
}

__global__ void cu_collapsed_matrix_initialize(uint dimension, uint n, uint* indices, float2* results, float2* sum) {
    int x = blockIdx.x * blockDim.x+ threadIdx.x;
    cu_set_indices(indices, dimension, x);
    // cu_initialize_results(cu_nCr(n,dimension), results, x);
    // cu_initialize_sum(n, sum, x);
}

__device__ float2 dev_grav_force(cuBody bodyA, cuBody bodyB) {
    const float GRAV_CONST = 5.0;
    const float2 vec12 = float2_subtract(bodyA.position, bodyB.position);
    return make_float2(bodyA.position.x, bodyA.position.y);
    const float force_mag = -GRAV_CONST * (bodyA.mass*bodyB.mass)/float2_sqr_length(vec12);
    const float force_angle = float2_angle(vec12);
    return float2_from_magnitude_and_angle(force_mag, force_angle);
}

__global__ void calculate_gravity (const uint* indices, float2* results, const cuBody* bodies) {
    int x = blockIdx.x * blockDim.x+ threadIdx.x;
    if (x==0) {results[x] = bodies[1].position; return;}
    results[x] = dev_grav_force(bodies[indices[2*x]], bodies[indices[2*x+1]]);
}

//grid(length*dimension)
__global__ void calculate_sum (uint length, uint dimension, uint* indices, float2* results, float2* sum) {
    int x = blockIdx.x * blockDim.x+ threadIdx.x;
    for (uint i=0; i< (length*dimension); i++) {
    if (x == indices[i]) float2_add_equals(&sum[indices[i]], &results[i/dimension]);
    }
}

void cu_collapsed_matrix_new (CollapsedMatrix* from, cuCollapsedMatrix* to) {
    to->dimension = from->dimension;
    to->n = from->n;
    to->length = nCr(to->n, to->dimension);
    cudaMalloc((void**) &to->indices, from->length*from->dimension*sizeof(uint));
    cudaMalloc((void**) &to->results, from->length*sizeof(float2));
    cudaMalloc((void**) &to->sum, from->n*sizeof(float2));

    cudaMemset(&to->results, 0, from->length*sizeof(float2));
    cudaMemset(&to->sum, 0, from->n*sizeof(float2));
    // cudaMemcpy(to->indices, from->indices, from->length*from->dimension*sizeof(uint), cudaMemcpyHostToDevice);
    // cudaMemcpy(to->results, from->results, from->length*sizeof(float2), cudaMemcpyHostToDevice);
    // cudaMemcpy(to->sum, from->sum, from->n*sizeof(float2), cudaMemcpyHostToDevice);
}

int main() {
    uint n=3, d=2;
    srand(3);

    clock_t start_fun, end_fun, start_sum, end_sum;

    Body* bodies = (Body*) malloc(n*sizeof(Body));
    for (int i=0; i<n; i++) {body_random(bodies+i);}

    CollapsedMatrix* collapsed_matrix = collapsed_matrix_new(d, n);
    float2* sum = (float2*) malloc(n*sizeof(float2));
    uint* indices = (uint*) malloc(nCr(n,d)*sizeof(uint));
    ////////////////////////////
    cuCollapsedMatrix* dev_CM = (cuCollapsedMatrix*) malloc(sizeof(cuCollapsedMatrix));;
    cu_collapsed_matrix_new(collapsed_matrix, dev_CM);
    cuBody* dev_bodies;
    Body test_body[1];
    cudaMalloc((void**) &dev_bodies, n*sizeof(cuBody));
    cudaMemcpy(dev_bodies, bodies, n*sizeof(cuBody), cudaMemcpyHostToDevice);
    
    collapsed_matrix_apply_function(collapsed_matrix, bodies);
    collapsed_matrix_calculate_sum(collapsed_matrix);
    
    dim3 block(1);
    dim3 grid(1);
    // cu_collapsed_matrix_new<<<grid,1>>>(dev_CM, d, n);
    printf("Matrix created\n");
    block = dim3(dev_CM->n);
    grid = dim3(dev_CM->length/(dev_CM->n-1));
    printf("Initializing\n");
    cu_collapsed_matrix_initialize<<<grid,block>>>(dev_CM->dimension, dev_CM->n, dev_CM->indices, dev_CM->results, dev_CM->sum);
    printf("initialized\n");
    ////////////////////////////
    start_fun = clock();
    calculate_gravity<<<grid,block>>> (dev_CM->indices, dev_CM->results, dev_bodies);
    printf("Gravity calculated\n");
    block = dim3((int)sqrtf(n));
    grid = dim3(n/block.x);
    calculate_sum<<<grid,block>>> (dev_CM->length, dev_CM->dimension, dev_CM->indices, dev_CM->results, dev_CM->sum);
    printf("Gravity added\n");

    // cudaMemcpy(results, dev_CM->sum, collapsed_matrix->n*sizeof(float2), cudaMemcpyDeviceToHost);
    cudaMemcpy(sum, dev_CM->results, n*sizeof(float2), cudaMemcpyDeviceToHost);

    printf("Results transfered\n");
    // for (int i=0; i<nCr(n,d); i++) {grav_force(bodies);}
    end_fun = clock();

    start_sum = clock();
    //collapsed_matrix_calculate_sum(collapsed_matrix);
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
    printf("Apply function: %f\n", ((double) (end_fun-start_fun)) / CLOCKS_PER_SEC * 1000);
    printf("Apply sumation: %f\n", ((double) (end_sum-start_sum)) / CLOCKS_PER_SEC * 1000);
    
    body_print(&bodies[0]);
    body_print(&bodies[1]);
    printf("Sum[0] = \t");
    printf("(%06.3f,%06.3f)\n", collapsed_matrix->results[0].x, collapsed_matrix->results[0].y);
    printf("Dev_Sum[0] = \t");
    printf("(%06.3f,%06.3f)\n", sum[0].x, sum[0].y);

    return 0;
}