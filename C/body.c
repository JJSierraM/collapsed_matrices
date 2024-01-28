#include "body.h"
#include <stdio.h>

void body_initialize(Body* body) {
    body->mass=0.001;
    vector2D_initialize(&body->position);
    vector2D_initialize(&body->speed);
    vector2D_initialize(&body->accel);
}

//Random mass and position, speed and accel set to 0
void body_random(Body* body) {
    const float rand_mass = ((float)rand()/(float)RAND_MAX)*(1.0-0.001)+0.001; 
    Vector2D rand_pos;
    vector2D_random(&rand_pos, -1.0, 1.0); 
    body->mass = rand_mass;
    body->position = rand_pos;
    vector2D_initialize(&body->speed);
    vector2D_initialize(&body->accel);
}

Body* body_new() {
    Body* output = malloc(sizeof(Body));
    body_initialize(output);
    return output;
}

void body_print(Body* body) {
    printf("Mass:\t\t%06.3f\n", body->mass);
    printf("Position:\tx: %06.3f, y: %06.3f\n", body->position.x, body->position.y);
    printf("Speed:\t\tx: %06.3f, y: %06.3f\n", body->speed.x, body->speed.y);
    printf("Acceleration:\tx: %06.3f, y: %06.3f\n", body->accel.x, body->accel.y);
}
