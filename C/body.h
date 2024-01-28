#include "vector2D.h"

#ifndef BODYH
#define BODYH

typedef struct Body {
    float mass;
    Vector2D position;
    Vector2D speed;
    Vector2D accel; 
} Body;

void body_initialize(Body* body) ;

//Random mass and position, speed and accel set to 0
void body_random(Body* body) ;

Body* body_new() ;

void body_print(Body* body) ;

#endif