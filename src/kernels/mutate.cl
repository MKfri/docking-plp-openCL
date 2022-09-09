#ifndef MUTATE_CL_H
#define MUTATE_CL_H

#include "clStructs.h"

#include "RealConstants.cl"
#include "constants.cl"
#include "PLPConstants.cl"

#include "tyche_i.cl"
#include "randomExtra.cl"
#include "vector3.cl"
#include "Quat.cl"

// MUTATE DIHEDRAL:

Float StandardisedValue(Float dihedralAngle) {
  while (dihedralAngle >= PLUS_180_0f) {
    dihedralAngle -= PLUS_360_0f;
  }
  while (dihedralAngle < MINUS_180_0f) {
    dihedralAngle += PLUS_360_0f;
  }
  return dihedralAngle;
}

Float CorrectTetheredDihedral(Float initialValue, Float m_value, constant parametersForGPU* parameters) {
  Float maxDelta = parameters->max_dihedral;

  Float delta = StandardisedValue(m_value - initialValue);
  if (delta > maxDelta) {
    return StandardisedValue(initialValue + maxDelta);
  } else if (delta < -maxDelta) {
    return StandardisedValue(initialValue - maxDelta);
  } else {
      return m_value;//DO NOTHING (aka. return current value)
  }
}

void mutateDihedral(global Float* chromosomeDihedral, int numOfRotatableBonds, Float relStepSize, tyche_i_state* state, constant parametersForGPU* parameters) {
    
    Float absStepSize = relStepSize * parameters->dihedral_step;
    Float delta;
    int m_mode;
    Float initialValue;

    Float tempDihedral;

    for(int i = 0; i < numOfRotatableBonds; i++){
 
        tempDihedral = chromosomeDihedral[i];

        if (absStepSize > PLUS_0_0f) {
            m_mode = parameters->dihedralMode;
            
#ifdef USE_DOUBLE
            delta = PLUS_2_0f * absStepSize * tyche_i_double(state[0]) - absStepSize;
#else
            delta = PLUS_2_0f * absStepSize * tyche_i_float(state[0]) - absStepSize;
#endif


            if (m_mode == CHROM_m_mode_eMode_TETHERED) {
                //delta = 2.0f * absStepSize * tyche_i_float(state[0]) - absStepSize;
                initialValue = tempDihedral;
                tempDihedral = StandardisedValue(tempDihedral + delta);
                tempDihedral = CorrectTetheredDihedral(initialValue, tempDihedral, parameters);

            } else if (m_mode == CHROM_m_mode_eMode_FREE) {
                //delta = 2.0f * absStepSize * tyche_i_float(state[0]) - absStepSize;
                tempDihedral = StandardisedValue(tempDihedral + delta);

            }
        }

        chromosomeDihedral[i] = tempDihedral;
    }
}

// MUTATE CENTER OF MASS:

void CorrectTetheredCOM(Float* chromosomeCOM, Float* initCOM, constant parametersForGPU* parameters) {

    Float axis[3];
    Float unit[3];
    Float maxTrans = parameters->max_trans;

    // vector from initial COM to new point
    axis[0] = chromosomeCOM[0] - initCOM[0];
    axis[1] = chromosomeCOM[1] - initCOM[1];
    axis[2] = chromosomeCOM[2] - initCOM[2];

    // avoid square roots until necessary
    Float length = lengthNorm3((Float*)axis);
    Float length2 = length * length;
    unitVector3((Float*)axis, (Float*)unit, length);
    if (length2 > (maxTrans * maxTrans)) {
        // Stay just inside the boundary sphere
       chromosomeCOM[0] = initCOM[0] + (PLUS_0_999f * maxTrans * unit[0]);
       chromosomeCOM[1] = initCOM[1] + (PLUS_0_999f * maxTrans * unit[1]);
       chromosomeCOM[2] = initCOM[2] + (PLUS_0_999f * maxTrans * unit[2]);
    }
}

void mutateCOM(global Float* chromosomeCOM, Float relStepSize, tyche_i_state* state, constant parametersForGPU* parameters) {
    //chromosomeCOM[0] chromosomeCOM[1] chromosomeCOM[2]
    //Float absTransStepSize;
    //Float dist;
    Float axis[3];
    Float initCOM[3];

    int m_transMode = parameters->transMode;

    Float tempCOM[3];
    tempCOM[0] = chromosomeCOM[0];
    tempCOM[1] = chromosomeCOM[1];
    tempCOM[2] = chromosomeCOM[2];
    /*
    if(m_transMode == CHROM_m_mode_eMode_TETHERED) {
        absTransStepSize = relStepSize * parameters->trans_step;
        if (absTransStepSize > PLUS_0_0f) {
            dist = absTransStepSize * tyche_i_float(state[0]);
            //Get random unit vector to axis
            randomUnitVector3((Float*)axis, state);            
            //save init COM
            saveVector3((Float*)tempCOM, (Float*)initCOM);

            tempCOM[0] += dist * axis[0];
            tempCOM[1] += dist * axis[1];
            tempCOM[2] += dist * axis[2];

            CorrectTetheredCOM((Float*)tempCOM, (Float*)initCOM, parameters);
        }
    } else if(m_transMode == CHROM_m_mode_eMode_FREE) {
        absTransStepSize = relStepSize * parameters->trans_step;
        if (absTransStepSize > PLUS_0_0f) {
            dist = absTransStepSize * tyche_i_Float(state[0]);
            //Get random unit vector to axis
            randomUnitVector3((Float*)axis, state);

            tempCOM[0] += dist * axis[0];
            tempCOM[1] += dist * axis[1];
            tempCOM[2] += dist * axis[2];

        }
    }
    // FIXED: Do nothing
*/



    Float dist;
    Float absTransStepSize = relStepSize * parameters->trans_step;
    if (absTransStepSize > PLUS_0_0f) {
        
#ifdef USE_DOUBLE
        dist = absTransStepSize * tyche_i_double(state[0]);
#else
        dist = absTransStepSize * tyche_i_float(state[0]);
#endif
        if (m_transMode == CHROM_m_mode_eMode_TETHERED) {
            //Get random unit vector to axis
            randomUnitVector3((Float*)axis, state);            
            //save init COM
            saveVector3((Float*)tempCOM, (Float*)initCOM);

            tempCOM[0] += dist * axis[0];
            tempCOM[1] += dist * axis[1];
            tempCOM[2] += dist * axis[2];

            CorrectTetheredCOM((Float*)tempCOM, (Float*)initCOM, parameters);
        } else if (m_transMode == CHROM_m_mode_eMode_FREE) {
            //Get random unit vector to axis
            randomUnitVector3((Float*)axis, state);

            tempCOM[0] += dist * axis[0];
            tempCOM[1] += dist * axis[1];
            tempCOM[2] += dist * axis[2];
        }
    }

    chromosomeCOM[0] = tempCOM[0];
    chromosomeCOM[1] = tempCOM[1];
    chromosomeCOM[2] = tempCOM[2];
}

// MUTATE ORIENTATION:

Float StandardisedValueRotation(Float rotationAngle) {
  while (rotationAngle >= PLUS_PI) {
    rotationAngle -= PLUS_2_0f * PLUS_PI;
  }
  while (rotationAngle < -PLUS_PI) {
    rotationAngle += PLUS_2_0f * PLUS_PI;
  }
  return rotationAngle;
}

void CorrectTetheredOrientation(Float* chromosomeOrientation, Float* initOrientation, constant parametersForGPU* parameters) {
    // Check for orientation out of bounds
    Float maxRot = parameters->max_rot;
    
    Float s_qAlign;
    Float v_qAlign[3];

    Float s_initial;
    Float v_initial[3];

    Float s_conj;
    Float v_conj[3];

    ToQuat((Float*)chromosomeOrientation, &s_initial, (Float*)v_initial);

    ConjQuat(&s_initial, (Float*)v_initial, &s_conj, (Float*)v_conj);

    //Copy to
    s_qAlign = s_conj;
    saveVector3((Float*)v_conj, (Float*)v_qAlign);

    multiplyQuat(&s_initial, (Float*)v_initial, &s_qAlign, (Float*)v_qAlign);
    Float cosHalfTheta = s_qAlign;
    if (cosHalfTheta < MINUS_1_0f) {
        cosHalfTheta = MINUS_1_0f;
    } else if (cosHalfTheta > PLUS_1_0f) {
        cosHalfTheta = PLUS_1_0f;
    }
    // Theta is the rotation angle required to realign with reference
    Float theta = StandardisedValueRotation(2.0f * acos(cosHalfTheta));
    // Deal with pos and neg angles independently as we have to
    // invert the rotation axis for negative angles
    Float axis[3];
    Float negate[3];
    if (theta < -maxRot) {

        negateVector3((Float*)v_qAlign, (Float*)negate);

        divideVectorByScalar3((Float*)negate, sin(theta / PLUS_2_0f), (Float*)axis);

        // Adjust theta to bring the orientation just inside the tethered bound
        theta += PLUS_0_999f * maxRot;

        Rotate((Float*)chromosomeOrientation, (Float*)axis, theta);

    } else if (theta > maxRot) {
        
        divideVectorByScalar3((Float*)v_qAlign, sin(theta / PLUS_2_0f), (Float*)axis);

        // Adjust theta to bring the orientation just inside the tethered bound
        theta -= PLUS_0_999f * maxRot;

        Rotate((Float*)chromosomeOrientation, (Float*)axis, theta);
    }
}

void mutateOrientation(global Float* chromosomeOrientation, Float relStepSize, tyche_i_state* state, constant parametersForGPU* parameters) {
    //chromosomeOrientation[0] chromosomeOrientation[1] chromosomeOrientation[2]


    Float axis[3];
    Float initOrientation[3];

    int m_rotMode = parameters->rotMode;

    Float tempOrientation[3];
    tempOrientation[0] = chromosomeOrientation[0];
    tempOrientation[1] = chromosomeOrientation[1];
    tempOrientation[2] = chromosomeOrientation[2];
    
    /*if (m_rotMode == CHROM_m_mode_eMode_TETHERED) {
        absRotStepSize = relStepSize * parameters->rot_step;
        if (absRotStepSize > PLUS_0_0f) {
            theta = absRotStepSize * tyche_i_float(state[0]);

            // Get random unit vector to axis
            randomUnitVector3((Float*)axis, state);
            // save init Orientation
            saveVector3((Float*)tempOrientation, (Float*)initOrientation);
            // Rotate:
            Rotate((Float*)tempOrientation, (Float*)axis, theta);
            
            CorrectTetheredOrientation((Float*)tempOrientation, (Float*)initOrientation, parameters);

        }
    } else if (m_rotMode == CHROM_m_mode_eMode_FREE) {
        absRotStepSize = relStepSize * parameters->rot_step;
        if (absRotStepSize > PLUS_0_0f) {
            theta = absRotStepSize * tyche_i_float(state[0]);

            //Get random unit vector to axis
            randomUnitVector3((Float*)axis, state);
            // Rotate:
            Rotate((Float*)tempOrientation, (Float*)axis, theta);
        }
    }*/
    // FIXED: Do nothing

    Float theta;
    Float absRotStepSize = relStepSize * parameters->rot_step;
    if (absRotStepSize > PLUS_0_0f) {

#ifdef USE_DOUBLE
        theta = absRotStepSize * tyche_i_double(state[0]);
#else
        theta = absRotStepSize * tyche_i_float(state[0]);
#endif
        if (m_rotMode == CHROM_m_mode_eMode_TETHERED) {
            // Get random unit vector to axis
            randomUnitVector3((Float*)axis, state);
            // save init Orientation
            saveVector3((Float*)tempOrientation, (Float*)initOrientation);
            // Rotate:
            Rotate((Float*)tempOrientation, (Float*)axis, theta);
            
            CorrectTetheredOrientation((Float*)tempOrientation, (Float*)initOrientation, parameters);
        } else if (m_rotMode == CHROM_m_mode_eMode_FREE) {
            // Get random unit vector to axis
            randomUnitVector3((Float*)axis, state);
            // Rotate:
            Rotate((Float*)tempOrientation, (Float*)axis, theta);
        }
    }

    chromosomeOrientation[0] = tempOrientation[0];
    chromosomeOrientation[1] = tempOrientation[1];
    chromosomeOrientation[2] = tempOrientation[2];
}

// MUTATE OCCUPANCY ELEMENT:

// mutateOccupancyElement(...) TODO


// MAIN MUTATE:

void mutate(global Float* chromosome, int chromStoreLen, Float relStepSize, tyche_i_state* state, constant parametersForGPU* parameters) {
    mutateDihedral(chromosome, chromStoreLen - CHROM_SUBTRACT_FOR_ROTATABLE_BONDS_LENGTH, relStepSize, state, parameters);
    mutateCOM(&(chromosome[chromStoreLen - CHROM_SUBTRACT_FOR_CENTER_OF_MASS_1]), relStepSize, state, parameters);
    mutateOrientation(&(chromosome[chromStoreLen - CHROM_SUBTRACT_FOR_ORIENTATION_1]), relStepSize, state, parameters);
    // mutateOccupancyElement(...) TODO
}

// CAUCHY MUTATE:

void CauchyMutate(global Float* chromosome, int chromStoreLen, Float mean, Float variance, tyche_i_state* state, constant parametersForGPU* parameters) {
  // Need to convert the Cauchy random variable to a positive number
  // and use this as the relative step size for mutation
  
  Float relStepSize = PositiveCauchyRandomWMeanVariance(mean, variance, state);
  
  mutate(chromosome, chromStoreLen, relStepSize, state, parameters);
}

#endif
