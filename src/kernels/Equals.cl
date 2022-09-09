#ifndef EQUALS_CL_H
#define EQUALS_CL_H

#include "clStructs.h"

#include "RealConstants.cl"
#include "constants.cl"

#include "mutate.cl"
#include "vector3.cl"

// TODO: CompareVectorOccupancyElement(...)

Float CompareVectorDihedral(global Float* g1, global Float* g2, constant parametersForGPU* parameters) {

    Float retVal = PLUS_0_0f;
    Float absDiff;
    
    Float stepSize = parameters->dihedral_step;
    if (stepSize > PLUS_0_0f) {
        for(int i = 0; i < parameters->numDihedralElements; i++) {
            absDiff = fabs(StandardisedValue(g1[i] - g2[i]));
            retVal = fmax(retVal, absDiff / stepSize);
        }
    }
    return retVal;
}

Float CompareVectorPosition(global Float* g1, global Float* g2, constant parametersForGPU* parameters) {
    
    Float retVal = PLUS_0_0f;

    // COM:
    if(parameters->transMode != CHROM_m_mode_eMode_FIXED) {

        int startComIndex = (parameters->chromStoreLen) - CHROM_SUBTRACT_FOR_CENTER_OF_MASS_1;
        Float com1[3];
        com1[0] = g1[startComIndex];
        com1[1] = g1[startComIndex + 1];
        com1[2] = g1[startComIndex + 2];
        Float com2[3];
        com2[0] = g2[startComIndex];
        com2[1] = g2[startComIndex + 1];
        com2[2] = g2[startComIndex + 2];

        Float transStepSize = parameters->trans_step;
        // Compare distance between centres of mass
        if (transStepSize > PLUS_0_0f) {
            Float absDiff = distance2Points3((Float*)com1, (Float*)com2);
            Float relDiff = absDiff / transStepSize;
            retVal = fmax(retVal, relDiff);
        }
    }
    // Orientation:
    if(parameters->rotMode != CHROM_m_mode_eMode_FIXED) {

        int startOrientationIndex = (parameters->chromStoreLen) - CHROM_SUBTRACT_FOR_ORIENTATION_1;
        Float orientation1[3];
        orientation1[0] = g1[startOrientationIndex];
        orientation1[1] = g1[startOrientationIndex + 1];
        orientation1[2] = g1[startOrientationIndex + 2];
        Float orientation2[3];
        orientation2[0] = g2[startOrientationIndex];
        orientation2[1] = g2[startOrientationIndex + 1];
        orientation2[2] = g2[startOrientationIndex + 2];

        Float rotStepSize = parameters->rot_step;
        // Compare orientations
        if (rotStepSize > PLUS_0_0f) {
            // Determine the difference between the two orientations
            // in terms of the axis/angle needed to align them
            // q.s = std::cos(phi / 2)
            Float otherOrientation_v[3];
            Float otherOrientation_s;
            ToQuat((Float*)orientation2, &otherOrientation_s, (Float*)otherOrientation_v);

            Float m_orientation_v[3];
            Float m_orientation_s;
            ToQuat((Float*)orientation1, &m_orientation_s, (Float*)m_orientation_v);

            Float m_orientation_Conj_v[3];
            Float m_orientation_Conj_s;
            ConjQuat(&m_orientation_s, (Float*)m_orientation_v, &m_orientation_Conj_s, (Float*)m_orientation_Conj_v);

            Float qAlign_v[3];
            Float qAlign_s;
            multiplyQuatResult(&otherOrientation_s, (Float*)otherOrientation_v, &m_orientation_Conj_s, (Float*)m_orientation_Conj_v, &qAlign_s, (Float*)qAlign_v);

            Float cosHalfTheta = qAlign_s;
            if (cosHalfTheta < MINUS_1_0f) {
                cosHalfTheta = MINUS_1_0f;
            } else if (cosHalfTheta > PLUS_1_0f) {
                cosHalfTheta = PLUS_1_0f;
            }
            Float absDiff = fabs(StandardisedValueRotation(2.0f * acos(cosHalfTheta)));
            Float relDiff = absDiff / rotStepSize;
            retVal = fmax(retVal, relDiff);
        }
    }
    return retVal;
}

// 1 = not equals, 0 = equals
int EqualsGenome(global Float* g1, global Float* g2, constant parametersForGPU* parameters) {
    
    // No check for same length!

    // Checks for maximum difference of any of the underlying chromosome elements;
    // Dihedral element has one float
    // COM and Orientation has 3 each.

    // TODO: CompareVectorOccupancyElement(...)
    Float retValDihedral = CompareVectorDihedral(g1, g2, parameters);
    Float retValPosition = CompareVectorPosition(g1, g2, parameters);
    Float retVal = fmax(retValDihedral, retValPosition);

    return (retVal > parameters->equality_threshold);
}

#endif
