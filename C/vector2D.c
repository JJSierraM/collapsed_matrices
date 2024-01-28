#include <stdlib.h>
#include "vector2D.h"

//Vector 2D ---

void vector2D_initialize(Vector2D* vector2D) {
    vector2D->x=0;
    vector2D->y=0;
}

void vector2D_set(Vector2D* vector2D, const float x, const float y) {
    vector2D->x=x;
    vector2D->y=y;
}

void vector2D_random(Vector2D* vector2D, const float from, const float to) {
    vector2D->x = ((float)rand()/(float)RAND_MAX)*(to-from)+from;
    vector2D->y = ((float)rand()/(float)RAND_MAX)*(to-from)+from;
}

Vector2D vector2D_new(const float x, const float y) {
    Vector2D output;
    vector2D_set(&output, x, y);
    return output;
}

Vector2D vector2D_add(Vector2D a, Vector2D b) {
    Vector2D output;
    output.x = a.x + b.x;
    output.y = a.y + b.y;
    return output;
}

Vector2D vector2D_add_equals(Vector2D* a, const Vector2D* b) {
    a->x = a->x + b->x;
    a->y = a->y + b->y;
}

Vector2D vector2D_subtract(Vector2D a, Vector2D b) {
    Vector2D output;
    output.x = a.x - b.x;
    output.y = a.y - b.y;
    return output;
}

Vector2D vector2D_subtract_equals(Vector2D* a, const Vector2D* b) {
    a->x = a->x - b->x;
    a->y = a->y - b->y;
}

float vector2D_sqr_length(Vector2D vector) {
    float output;
    output = vector.x*vector.x+vector.y*vector.y;
    return output;
}

float vector2D_length(Vector2D vector) {
    return sqrtf(vector2D_sqr_length(vector));
}

inline float vector2D_angle(Vector2D vector) {
    return atan2f(vector.y, vector.x);
}

Vector2D vector2D_from_magnitude_and_angle(float magnitude, float angle) {
    Vector2D output;
    output.x = magnitude*cosf(angle);
    output.y = magnitude*sinf(angle);
    return output;
}