#ifndef QUAT_CL_H
#define QUAT_CL_H

#include "clStructs.h"
#include "RealConstants.cl"
#include "vector3.cl"

void RbtQuat(Float* axis, Float phi, Float* s, Float* v) {
    Float halfPhi = PLUS_0_5f * phi;
    *s = cos(halfPhi);// Scalar component

    Float unit[3];
    unitVector3((Float*)axis, (Float*)unit, lengthNorm3((Float*)axis));

    // Vector component
    v[0] = sin(halfPhi) * unit[0];
    v[1] = sin(halfPhi) * unit[1];
    v[2] = sin(halfPhi) * unit[2];
}

void ToQuat(Float* chromosomeOrientation, Float* s, Float* v) {
  Float c1 = cos(chromosomeOrientation[0] / PLUS_2_0f);
  Float s1 = sin(chromosomeOrientation[0] / PLUS_2_0f);
  Float c2 = cos(chromosomeOrientation[1] / PLUS_2_0f);
  Float s2 = sin(chromosomeOrientation[1] / PLUS_2_0f);
  Float c3 = cos(chromosomeOrientation[2] / PLUS_2_0f);
  Float s3 = sin(chromosomeOrientation[2] / PLUS_2_0f);
  Float c1c2 = c1 * c2;
  Float s1s2 = s1 * s2;

  // Constructor with initial values (vector component passed as 3 RbtDouble's)
  // float s1, float vx, float vy, float vz
  *s = c1c2 * c3 - s1s2 * s3;
  v[0] = c1c2 * s3 + s1s2 * c3;
  v[1] = c1 * s2 * c3 - s1 * c2 * s3;
  v[2] = s1 * c2 * c3 + c1 * s2 * s3;
}

void multiplyQuat(Float* s, Float* v, Float* s_original, Float* v_original) {
    // Constructor with float and vector (float, vector)

    Float s_new = (*s) * (*s_original) - dotProduct3((Float*)v, (Float*)v_original);
    
    Float v_p1[3];
    multiplyVectorByScalar3((Float*)v_original, (*s), (Float*)v_p1);
    
    Float v_p2[3];
    multiplyVectorByScalar3((Float*)v, (*s_original), (Float*)v_p2);
    
    Float v_cross[3];
    crossProduct3((Float*)v, (Float*)v_original, (Float*)v_cross);
    
    v_original[0] = v_p1[0] + v_p2[0] + v_cross[0];
    v_original[1] = v_p1[1] + v_p2[1] + v_cross[1];
    v_original[2] = v_p1[2] + v_p2[2] + v_cross[2];
    *s_original = s_new;

}

void multiplyQuatResult(Float* s, Float* v, Float* s_original, Float* v_original, Float* s_result, Float* v_result) {
    // Constructor with float and vector (float, vector)

    Float s_new = (*s) * (*s_original) - dotProduct3((Float*)v, (Float*)v_original);
    
    Float v_p1[3];
    multiplyVectorByScalar3((Float*)v_original, (*s), (Float*)v_p1);
    
    Float v_p2[3];
    multiplyVectorByScalar3((Float*)v, (*s_original), (Float*)v_p2);
    
    Float v_cross[3];
    crossProduct3((Float*)v, (Float*)v_original, (Float*)v_cross);
    
    v_result[0] = v_p1[0] + v_p2[0] + v_cross[0];
    v_result[1] = v_p1[1] + v_p2[1] + v_cross[1];
    v_result[2] = v_p1[2] + v_p2[2] + v_cross[2];
    *s_result = s_new;

}

void FromQuat(Float* s, Float* v, Float* euler) {
Float test = (v[0] * v[2]) + (v[1] * (*s));
  if (test > PLUS_0_499999f) { // singularity at north pole
    euler[0] = PLUS_2_0f * atan2(v[0], (*s));
    euler[1] = PLUS_PI / PLUS_2_0f;
    euler[2] = PLUS_0_0f;
  } else if (test < MINUS_0_499999f) { // singularity at south pole
    euler[0] = MINUS_2_0f * atan2(v[0], (*s));
    euler[1] = -PLUS_PI / PLUS_2_0f;
    euler[2] = PLUS_0_0f;
  } else {
    Float sq[3];
    multiply2Vectors3((Float*)v, (Float*)v, (Float*)sq);

    euler[0] = atan2(PLUS_2_0f * (v[2] * (*s) - v[0] * v[1]),
                           PLUS_1_0f - PLUS_2_0f * (sq[2] + sq[1]));
    euler[1] = asin(PLUS_2_0f * test);
    euler[2] = atan2(PLUS_2_0f * (v[0] * (*s) - v[2] * v[1]),
                        PLUS_1_0f - PLUS_2_0f * (sq[0] + sq[1]));
  }
}

void ConjQuat(Float* s, Float* v, Float* s_conj, Float* v_conj) {
    // Returns conjugate
    *s_conj = *s;
    v_conj[0] = -v[0];
    v_conj[1] = -v[1];
    v_conj[2] = -v[2];
}

void Rotate(Float* chromosomeOrientation, Float* axis, Float theta) {
    Float s;
    Float v[3];

    Float s_original;
    Float v_original[3];

    // a. get Quat for random axis
    RbtQuat((Float*)axis, theta, &s, (Float*)v);// s and v variable
    // b. get Quat for  initial orientation
    ToQuat((Float*)chromosomeOrientation, &s_original, (Float*)v_original);
    // c. multiply q * ToQuat()  Multiplication (non-commutative)
    multiplyQuat(&s, (Float*)v, &s_original, (Float*)v_original);
    // d. get back euler
    FromQuat(&s_original, (Float*)v_original, (Float*)chromosomeOrientation);
}

void RotateUsingQuat(Float* s, Float* v, Float* w, Float* result) {
  // 1. Part: (s * s - v.Dot(v)) * w
  Float wScalar = ( (*s) * (*s) ) - dotProduct3((Float*)v, (Float*)v);
  Float wScaled[3];
  multiplyVectorByScalar3((Float*)w, wScalar, (Float*)wScaled);
  // 2. Part: 2 * s * v.Cross(w)
  Float ss = PLUS_2_0f * (*s);
  Float vCrossW[3];
  crossProduct3((Float*)v, (Float*)w, (Float*)vCrossW);
  Float vCrossWscaled[3];
  multiplyVectorByScalar3((Float*)vCrossW, ss, (Float*)vCrossWscaled);
  // 3. Part: 2 * v * v.Dot(w);
  Float vv[3];
  multiplyVectorByScalar3((Float*)v, PLUS_2_0f, (Float*)vv);
  Float vDotW = dotProduct3((Float*)v, (Float*)w);
  Float vScaled[3];
  multiplyVectorByScalar3((Float*)vv, vDotW, (Float*)vScaled);

  // Result:
  result[0] = wScaled[0] + vCrossWscaled[0] + vScaled[0];
  result[1] = wScaled[1] + vCrossWscaled[1] + vScaled[1];
  result[2] = wScaled[2] + vCrossWscaled[2] + vScaled[2];
}

void GetQuatFromAlignVectors(Float* v, Float* ref, Float* q_s, Float* q_v) {
  // Set Default:
  *q_s = PLUS_1_0f;
  q_v[0] = PLUS_0_0f;
  q_v[1] = PLUS_0_0f;
  q_v[2] = PLUS_0_0f;
  
  // Unitise the two vectors
  Float len = lengthNorm3(v);
  Float refLen = lengthNorm3(ref);
  if ((len < 0.001f) || (refLen < 0.001f)) {
    // RbtBadArgument
    printf("[Quat] [GetQuatFromAlignVectors] [Zero length vector (v or ref)]\n");
  }
  Float vUnit[3];
  Float refUnit[3];
  divideVectorByScalar3(v, len, (Float*)vUnit);
  divideVectorByScalar3(ref, refLen, (Float*)refUnit);
  // Determine the rotation axis and angle needed to overlay the two vectors
  Float axis[3];
  crossProduct3((Float*)vUnit, (Float*)refUnit, (Float*)axis);
  // DM 15 March 2006: check for zero-length rotation axis
  // This indicates the vectors are already aligned
  Float axisLen = lengthNorm3((Float*)axis);
  if (axisLen > PLUS_0_001f) {
    Float cosPhi = dotProduct3((Float*)vUnit, (Float*)refUnit);
    if (cosPhi < MINUS_1_0f) {
      cosPhi = MINUS_1_0f;
    } else if (cosPhi > PLUS_1_0f) {
      cosPhi = PLUS_1_0f;
    }
    //errno = 0;
    // Convert rotation axis and angle to a quaternion
    Float halfPhi = PLUS_0_5f * acos(cosPhi);
    if ((halfPhi > PLUS_0_001f)) { // && (errno != EDOM)
      Float axisUnit[3];
      divideVectorByScalar3((Float*)axis, axisLen, (Float*)axisUnit);
      // Set return Quat
      *q_s = cos(halfPhi);
      multiplyVectorByScalar3((Float*)axisUnit, sin(halfPhi), q_v);
    }

  }

}

void GetQuatFromAlignAxes(PrincipalAxesSyncGPU* prAxes, PrincipalAxesSyncGPU* refAxes, Float* s, Float* v) {

  // 1) Determine the quaternion needed to align axis1 with reference
  Float q1_v[3];
  Float q1_s;
  GetQuatFromAlignVectors((Float*)prAxes->axis1, (Float*)refAxes->axis1, &q1_s, (Float*)q1_v);
  // 2) Apply the transformation to axis2
  Float axis2[3];
  RotateUsingQuat(&q1_s, (Float*)q1_v, (Float*)prAxes->axis2, (Float*)axis2);
  // 3) Determine the quaternion needed to align axis2 with reference
  Float q2_v[3];
  Float q2_s;
  GetQuatFromAlignVectors((Float*)axis2, (Float*)refAxes->axis2, &q2_s, (Float*)q2_v);
  // Return the quaternion product (equivalent to both transformations combined)
  multiplyQuatResult(&q2_s, (Float*)q2_v, &q1_s, (Float*)q1_v, s, v);// q2 * q1
}

#endif
