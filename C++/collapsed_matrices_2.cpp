#include<stdio.h>
#include<stdlib.h>
#include <time.h>
#include <math.h>

typedef unsigned short int ushort;
typedef unsigned int uint;

template<typename T, typename U>
struct CollapsedMatrix
{
    public:
    int dimension;
    int n;
    int length;
    int *indices;
    T *results;
    T *sum;

    CollapsedMatrix(int dimension_v, int n_v) : 
        dimension(dimension_v),
        n(n_v),
        length(nCr(n, dimension)),
        indices((int*) malloc(length*dimension*sizeof(int))),
        results((T*) malloc(length*sizeof(T))),
        sum((T*) malloc(n*sizeof(T)))
    {
        set_indices_2();
        initialize_results();
        initialize_sum();
    }
    ~CollapsedMatrix()
    {
        free(indices);
        free(results);
        free(sum);
    }

    static int nCr (int n, int r){
        int output = 1;
        int factorial = 1;
        if (n<r) {output = 0;}
        else {
            int i;
            for (i=n+1-r; i<n+1; i++) {output*=i;}
            for (i=1; i<r+1; i++) {factorial*=i;}
        }
        return output/factorial;
    }

    void set_indices() {
        for (int i=0; i<length; i++){
            int prev = 0;
            int* item = indices+(i*dimension);
            for (int k=0; k<dimension; k++){
                for (int j=0; j<999; j++){
                    if ((i-prev) < nCr(dimension-k+j, j)) {
                        item[k] = j + dimension -k -1;
                        prev += nCr(dimension -k +j -1, j-1);
                        break;
                    }
                }
            }
        }
    }

    void initialize_sum() {
        for (int i=0; i<n; i++) sum[i] = T();
    }

    void initialize_results() {
        for (int i=0; i<length; i++) results[i] = T();
    }

    void apply_function(T f(const U*), const U* items){
        U* selected_items = (U*) malloc(dimension*sizeof(U));
        for (int i=0; i<length; i++){
            //printf("%d - Selected: ", i);
            for (int index=0; index<dimension; index++){
                selected_items[index] = items[indices[(int)i*dimension+index]];
            //    printf("%d ", selected_items[index]);
            }
            //printf("\t");

            //results[i] = 
            f(selected_items);
        }
        free(selected_items);
    }

    void apply_function_2(T f(const U*, const int*), const U* items){
        int* selected_items = (int*) malloc(dimension*sizeof(int));
        //# pragma omp parallel private ( selected_items )
        {
            //#pragma omp for
            for (int i=0; i<length; i++){
                for (int index=0; index<dimension; index++){
                    selected_items[index] = indices[i*dimension+index];
                //    printf("%d ", selected_items[index]);
                }
                results[i] = f(items, selected_items);
            }
        }
        free(selected_items);
    }

    void apply_function_3(T f(const U**), const U* items){
        const U** ordered_items = (const U**) malloc(length*dimension*sizeof(U*));

        //# pragma omp parallel private ( ordered_items )
        {
            //#pragma omp for
            for (int i=0; i<length*dimension; i++){
                ordered_items[i] = items+indices[i];
            }
            //#pragma omp for
            for (int i=0; i<length; i++){
                results[i] = f(ordered_items+(i*dimension));
            }
        }
        free(ordered_items);
    }

    void calculate_sum () {
        for (int item=0; item<n; item++){
            int i = 0;
            for (int index=0; index<length; index++){
                for (int k=0; k<dimension; k++){
                    if ((indices+(index*dimension))[k] == item){
                        sum[item] += results[i];
                        break;
                    }
                }
                i++;
            }
        }
    }

    void calculate_sum_2() {
        int i;
        # pragma omp parallel private ( i )
        {
            # pragma omp for
            for (i = 0; i < (length*dimension); i++){
                sum[indices[i]] += results[i/dimension];
            }
        }
    }

    //A ELIMINAR #########################################################
    static void print_arr(int* arr, int dim){
        for (int i=0; i<dim; i++) printf("%d",arr[i]);
        printf("\n");
    }
    // ###################################################################
    void transfer_array(int* arr1, int* arr2, int dim_arr2){
        for (int i=0; i<dim_arr2; i++){
            *(arr1+i)=*(arr2+i);
        }
    }

    //DEPRECATED. USE set_indices() INSTEAD
    void set_indices_2 () {
        int *item = (int*) malloc(dimension*sizeof(int));
        int i,j,k;
        for (i=0; i<dimension; i++) item[i] = dimension -1 -i;
        transfer_array(indices+0, item, dimension);
        for (i=1; i<length; i++){
            for (j=dimension-1; j>=0; j--){
                if (j == 0)
                {
                    item[0]++;
                    for (k=j+1; k<dimension; k++){
                        item[k] = dimension -1 -k;
                    }
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
            transfer_array(indices+(i*dimension), item, dimension);
        }
        free(item);
    }
};

struct Vector2
{
    float x;
    float y;

    Vector2(float x_v=0.0, float y_v=0.0) : 
    x(x_v), y(y_v)
    {}

    Vector2& operator= (const Vector2& other){
        x = other.x;
        y = other.y;
        return *this;
    }

    Vector2 operator+(const Vector2& other) const {
        return Vector2(x+other.x, y+other.y);
    }

    Vector2 operator-(const Vector2& other) const {
        return Vector2(x-other.x, y-other.y);
    }

    bool operator==(const Vector2& other) const {
        return (x == other.x && y == other.y);
    }

    Vector2 operator+=(const Vector2& other)  {
        x = x + other.x;
        y = y + other.y;
        return *this;
    }

    float length() {
        return sqrtf(x*x+y*y);
    }

    float sqr_length() const {
        return (x*x+y*y);
    }

    float angle() const {
        return atan2f(y,x);
    }

    void print() {
        printf("x: %06.3f, y: %06.3f\n", x, y);
    }

    static Vector2 from_magnitude_and_angle(float magnitude, float angle) {
        return Vector2(magnitude*cosf(angle), magnitude*sinf(angle));
    }

    static Vector2 random(const float from, const float to) {
        return Vector2(((float)rand()/(float)RAND_MAX)*(to-from)+from, ((float)rand()/(float)RAND_MAX)*(to-from)+from);
    }
};


struct Body{
    float mass;
    Vector2 position;
    Vector2 speed;
    Vector2 accel; 

    Body(float mass_v = 0.0001, Vector2 position_v = Vector2(), Vector2 speed_v = Vector2(), Vector2 accel_v = Vector2()) : 
    mass(mass_v),
    position(position_v),
    speed(speed_v),
    accel(accel_v)
    {}

    static Body random (float mass_from=0.0001, float mass_to=1, float vec_from=-1, float vec_to=1) {
        const float rand_mass = ((float)rand()/(float)RAND_MAX)*(mass_to-mass_from)+mass_from; 
        const Vector2 rand_pos = Vector2::random(vec_from, vec_to); 
        return Body(rand_mass, rand_pos);
    } 
    void print(){
        printf("Mass:\t\t%06.3f\n", mass);
        printf("Position:\tx: %06.3f, y: %06.3f\n", position.x, position.y);
        printf("Speed:\t\tx: %06.3f, y: %06.3f\n", speed.x, speed.y);
        printf("Acceleration:\tx: %06.3f, y: %06.3f\n", accel.x, accel.y);
    }
};

static Vector2 grav_force(const Body* bodies) {
    const float GRAV_CONST = 5;
    const Vector2 vec12 = bodies[0].position-bodies[1].position;
    const float force_mag = -GRAV_CONST * (bodies[0].mass*bodies[1].mass)/vec12.sqr_length();
    const float force_angle = vec12.angle();
    return Vector2::from_magnitude_and_angle(force_mag, force_angle);
}

static Vector2 grav_force_2(const Body* bodies, const int* indices){
    const float GRAV_CONST = 5;
    const Vector2 vec12 = bodies[indices[0]].position-bodies[indices[1]].position;
    const float force_mag = -GRAV_CONST * (bodies[indices[0]].mass*bodies[indices[1]].mass)/vec12.sqr_length();
    const float force_angle = vec12.angle();
    return Vector2::from_magnitude_and_angle(force_mag, force_angle);
}

static Vector2 grav_force_3(const Body** bodies){
    const float GRAV_CONST = 5;
    const Vector2 vec12 = (*bodies[0]).position-(*bodies[1]).position;
    const float force_mag = -GRAV_CONST * ((*bodies[0]).mass*(*bodies[1]).mass)/vec12.sqr_length();
    const float force_angle = vec12.angle();
    return Vector2::from_magnitude_and_angle(force_mag, force_angle);
}

static void apply_grav_force(CollapsedMatrix<Vector2,Body>* collapsed_matrix, const Body* items){
    int length = collapsed_matrix->length;
    //printf("Length applied: %d\n", length);
    int dimension = collapsed_matrix->dimension;
    //printf("Dimension applied: %d\n", dimension);
    //printf("Result: %d\n", length*dimension);
    Body* ordered_items = (Body*) malloc(length*dimension*sizeof(Body));

    //# pragma omp parallel private ( ordered_items )
    {
        for (int i=0; i<dimension*length; i++){
            ordered_items[i] = items[collapsed_matrix->indices[i]];
        }
        float GRAV_CONST = 5;
        Vector2 vec12;
        float force_mag=0.0;
        float force_angle=0.0;
        //printf("Items ordered\n");
        //#pragma omp for
        for (int i=0; i<length; i++){
            //collapsed_matrix->results[i] = 
            //grav_force(ordered_items+(i*dimension));
            //printf("Calculated item %d\n", i);
            {
                vec12 = (ordered_items+(i*dimension))[0].position-(ordered_items+(i*dimension))[0].position;
                force_mag = -GRAV_CONST * ((ordered_items+(i*dimension))[0].mass*(ordered_items+(i*dimension))[0].mass)/vec12.sqr_length();
                force_angle = vec12.angle();
                //collapsed_matrix->results[i] = Vector2::from_magnitude_and_angle(force_mag, force_angle);
            }
        }
        //printf("Items calculated\n");
    }
    free(ordered_items);
}

int add_ints(const int* ints){
    //printf("Resultado: %d\n", ints[0]+ints[1]);
    return ints[0]+ints[1];
}

int main () {
    int dimension = 2;
    int n = 1024;
    
    CollapsedMatrix<Vector2,Body> collapsed_matrix(dimension, n);

    srand(1);
    //int ints[5] = {2,3,1,4,5};
    Body* bodies = (Body*) malloc(n*sizeof(Body)); for (int i=0; i<n; i++) {bodies[i]=Body::random();}

    clock_t start_fun, end_fun, start_sum, end_sum;

    //
    //Vector2 sample[1024] = {Vector2::random(-1.0, 1.0)};
    start_fun = clock();
    apply_grav_force(&collapsed_matrix, bodies);
    //collapsed_matrix.apply_function(grav_force, bodies);
    //for (int i=0; i<1023; i++) {sample[i] = grav_force(bodies+i);}
    end_fun = clock();
    start_sum = clock();
    collapsed_matrix.calculate_sum_2();
    end_sum = clock();
    

    printf("\n------\n");
    /*
    for (int i=0; i<collapsed_matrix.length; i++){
        printf("(");
        for (int j=0; j<dimension; j++){
            printf("%d,", collapsed_matrix.indices[i*dimension+j]);
        }
        printf(")\n");
    */
    /*
    printf("Cuerpos:\n");
    for (int i=0; i<n; i++){
        printf("Body number %d:\n",i);
        bodies[i].print();
        //printf("%d - %d\n", i, collapsed_matrix.results[i]);
    }
    printf("Resultados:\n");
    for (int i=0; i<collapsed_matrix.length; i++){
        collapsed_matrix.results[i].print();
        //printf("%d - %d\n", i, collapsed_matrix.results[i]);
    }
    printf("Sumas:\n");
    for (int i=0; i<collapsed_matrix.n; i++){
        collapsed_matrix.sum[i].print();
        //printf("%d - %d\n", i, collapsed_matrix.sum[i]);
    }
    */
    ///*
    printf("\nTime spent:\n");
    printf("Apply function: %f\n", ((double) (end_fun-start_fun)) / CLOCKS_PER_SEC * 1000);
    printf("Apply sumation: %f\n", ((double) (end_sum-start_sum)) / CLOCKS_PER_SEC * 1000);
    //*/

    free(bodies);
    return 0;
}