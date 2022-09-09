#ifndef VECTOR3_CL_H
#define VECTOR3_CL_H

#include "RealConstants.cl"
#include "tyche_i.cl"

void randomUnitVector3(Float* vector, tyche_i_state* state) {
#ifdef USE_DOUBLE
    vector[2] = PLUS_2_0f * tyche_i_double(state[0]) - PLUS_1_0f;
    Float t = PLUS_2_0f * PLUS_PI * tyche_i_double(state[0]);
#else
    vector[2] = PLUS_2_0f * tyche_i_float(state[0]) - PLUS_1_0f;
    Float t = PLUS_2_0f * PLUS_PI * tyche_i_float(state[0]);
#endif
    Float w = sqrt(PLUS_1_0f - vector[2] * vector[2]);
    vector[0] = w * cos(t);
    vector[1] = w * sin(t);
}

void unitVector3(Float* vector, Float* unit, Float length) {
    if(length > PLUS_0_0f) {
        unit[0] = vector[0] / length;
        unit[1] = vector[1] / length;
        unit[2] = vector[2] / length;
    } else {
        unit[0] = vector[0];
        unit[1] = vector[1];
        unit[2] = vector[2];
    }
}

void saveVector3(Float* what, Float* to) {
    to[0] = what[0];
    to[1] = what[1];
    to[2] = what[2];
}

Float lengthNorm3(Float* vector) {
    //xyz.norm() Returns magnitude of vector (or distance from origin)
   return sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);
}

Float lengthNormSq3(Float* vector) {
    //Returns squared magnitude of vector (or distance from origin)
   return vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2];
}

Float dotProduct3(Float* v1, Float* v2) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

void crossProduct3(Float* v1, Float* v2, Float* result) {
    result[0] = v1[1]*v2[2] - v1[2]*v2[1];
    result[1] = v1[2]*v2[0] - v1[0]*v2[2];
    result[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

void multiplyVectorByScalar3(Float* v, Float s, Float* result) {
    result[0] = v[0] * s;
    result[1] = v[1] * s;
    result[2] = v[2] * s;
}

void multiply2Vectors3(Float* v1, Float* v2, Float* result) {
    // Scalar product (coord * coord : component-wise multiplication) RESULT VECTOR!
    result[0] = v1[0] * v2[0];
    result[1] = v1[1] * v2[1];
    result[2] = v1[2] * v2[2];
}

void divideVectorByScalar3(Float* v, Float scalar, Float* result) {
    result[0] = v[0] / scalar;
    result[1] = v[1] / scalar;
    result[2] = v[2] / scalar;
}

void negateVector3(Float* v, Float* result) {
    result[0] = -v[0];
    result[1] = -v[1];
    result[2] = -v[2];
}

void subtract2Vectors3(Float* v1, Float* v2, Float* result) {
    result[0] = v1[0] - v2[0];
    result[1] = v1[1] - v2[1];
    result[2] = v1[2] - v2[2];
}

void add2Vectors3(Float* v1, Float* v2, Float* result) {
    result[0] = v1[0] + v2[0];
    result[1] = v1[1] + v2[1];
    result[2] = v1[2] + v2[2];
}

void zerosNxN(Float* matrix, int n) {
    int i;
    for(i = 0; i < n*n; i++) {
        matrix[i] = PLUS_0_0f;
    }
}

// Returns distance between two coords
Float distance2Points3(Float* v1, Float* v2) {
    Float difference[3];
    subtract2Vectors3(v2, v1, (Float*)difference);
    return lengthNorm3((Float*)difference);
}

// Returns squared distance between two coords
Float distance2PointsSq3(Float* v1, Float* v2) {
    Float difference[3];
    subtract2Vectors3(v2, v1, (Float*)difference);
    return lengthNormSq3((Float*)difference);
}

// 1 = VALID, 0 = NOT VALID
int isPointInsideCuboid3(Float* p, Float* min, Float* max) {
    return ( p[0] >=min[0] && p[0] <= max[0] &&
             p[1] >=min[1] && p[1] <= max[1] &&
             p[2] >=min[2] && p[2] <= max[2]);
}

#endif
