#include <stdlib.h>
#include <math.h>

#ifndef VECTOR2DH
#define VECTOR2DH

//Vector 2D ---
typedef struct Vector2D
{
    float x;
    float y;
} Vector2D;

void vector2D_initialize(Vector2D* vector2D) ;

void vector2D_set(Vector2D* vector2D, const float x, const float y) ;

void vector2D_random(Vector2D* vector2D, const float from, const float to) ;

Vector2D vector2D_new(const float x, const float y) ;

Vector2D vector2D_add(Vector2D a, Vector2D b) ;

Vector2D vector2D_add_equals(Vector2D* a, const Vector2D* b) ;

Vector2D vector2D_subtract(Vector2D a, Vector2D b) ;

Vector2D vector2D_subtract_equals(Vector2D* a, const Vector2D* b) ;

float vector2D_sqr_length(Vector2D vector) ;

float vector2D_length(Vector2D vector) ;

float vector2D_angle(Vector2D vector) ;

Vector2D vector2D_from_magnitude_and_angle(float magnitude, float angle) ;

#endif