#ifndef SYNC_TO_MODEL_CL_H
#define SYNC_TO_MODEL_CL_H

#include "clStructs.h"

#include "RealConstants.cl"
#include "constants.cl"

#include "vector3.cl"
#include "Quat.cl"
#include "dsyevh3.cl"

// TODO: setModelValueOccupancyElement(...)

// Returns dihedral formed between 3 vectors
Float Dihedral(global AtomGPUsmall* m_atom1, global AtomGPUsmall* m_atom2, global AtomGPUsmall* m_atom3, global AtomGPUsmall* m_atom4) {
    // 4 Atoms Recieved DONE!

    // 4 Coordinates From Atoms DONE!

    // 3 vectors
    Float v1[3];
    Float v2[3];
    Float v3[3];
        // v1
    v1[0] = m_atom1->x - m_atom2->x;
    v1[1] = m_atom1->y - m_atom2->y;
    v1[2] = m_atom1->z - m_atom2->z;
        // v2
    v2[0] = m_atom2->x - m_atom3->x;
    v2[1] = m_atom2->y - m_atom3->y;
    v2[2] = m_atom2->z - m_atom3->z;
        // v3
    v3[0] = m_atom3->x - m_atom4->x;
    v3[1] = m_atom3->y - m_atom4->y;
    v3[2] = m_atom3->z - m_atom4->z;

    // Calculate the cross products
    Float A[3];
    Float B[3];
    Float C[3];
    crossProduct3((Float*)v1, (Float*)v2, (Float*)A);
    crossProduct3((Float*)v2, (Float*)v3, (Float*)B);
    crossProduct3((Float*)v2, (Float*)A, (Float*)C);

    // Calculate the distances
    Float rA = lengthNorm3((Float*)A);
    Float rB = lengthNorm3((Float*)B);
    Float rC = lengthNorm3((Float*)C);

    // Calculate the sin and cos
    Float cos_phi = dotProduct3((Float*)A, (Float*)B) / (rA * rB);
    Float sin_phi = dotProduct3((Float*)C, (Float*)B) / (rC * rB);

    // Get phi and convert to degrees
    Float phi = -atan2(sin_phi, cos_phi);
    return phi * PLUS_180_0f / PLUS_PI;
}

inline global AtomGPUsmall* getAtomGPUsmallFromID(global AtomGPUsmall* ligandAtoms,int id, int popMaxSize) {
    return getAtomGPUsmallFromBase(popMaxSize, id - 1, ligandAtoms);
}

void setModelValueDihedral(global AtomGPUsmall* ligandAtoms, global Float* individual, global DihedralRefDataGPU* dihedralRefData, constant parametersForGPU* parameters) {
    // For Every Dihedral
    for(int i = 0; i < parameters->numDihedralElements; i++) {
        global AtomGPUsmall* m_atom1 = getAtomGPUsmallFromID(ligandAtoms, dihedralRefData[i].atom1ID, parameters->popMaxSize);
        global AtomGPUsmall* m_atom2 = getAtomGPUsmallFromID(ligandAtoms, dihedralRefData[i].atom2ID, parameters->popMaxSize);
        global AtomGPUsmall* m_atom3 = getAtomGPUsmallFromID(ligandAtoms, dihedralRefData[i].atom3ID, parameters->popMaxSize);
        global AtomGPUsmall* m_atom4 = getAtomGPUsmallFromID(ligandAtoms, dihedralRefData[i].atom4ID, parameters->popMaxSize);

        Float delta = individual[i] - Dihedral(m_atom1, m_atom2, m_atom3, m_atom4);

        // Only rotate if delta is non-zero
        if (fabs(delta) > PLUS_0_001f) {

            // Coords of atom 1 ( m_atom2 !!! NOT 1 !! NOT MISTAKE ! )
            Float coord1[3];
            coord1[0] = m_atom2->x;
            coord1[1] = m_atom2->y;
            coord1[2] = m_atom2->z;

            // Vector along the bond between atom 1 and atom 2 (rotation axis)
            // ( m_atom3 !!! NOT 2 !! NOT MISTAKE ! )
            Float coord2[3];
            coord2[0] = m_atom3->x;
            coord2[1] = m_atom3->y;
            coord2[2] = m_atom3->z;
            Float bondVector[3];
            subtract2Vectors3((Float*)coord2, (Float*)coord1, (Float*)bondVector);
            Float toOrigin[3];
            negateVector3((Float*)coord1, (Float*)toOrigin);
            Float quat_v[3];
            Float quat_s;
            RbtQuat((Float*)bondVector, delta * PLUS_PI / PLUS_180_0f, &quat_s, (Float*)quat_v);
            
            Float tempCoord[3];
            Float translated[3];
            Float rotated[3];
            Float translated2[3];
            for(int j=0; j < dihedralRefData[i].numRotAtoms; j++) {
                global AtomGPUsmall* tempAtom = getAtomGPUsmallFromID(ligandAtoms, dihedralRefData[i].rotAtomsIDs[j], parameters->popMaxSize);
                tempCoord[0] = tempAtom->x;
                tempCoord[1] = tempAtom->y;
                tempCoord[2] = tempAtom->z;
                // Translate(toOrigin); (add coordinates together)
                add2Vectors3((Float*)tempCoord, (Float*)toOrigin, (Float*)translated);
                // RotateUsingQuat(quat);
                RotateUsingQuat(&quat_s, (Float*)quat_v, (Float*)translated, (Float*)rotated);
                // Translate(coord1); (add coordinates together)
                add2Vectors3((Float*)rotated, (Float*)coord1, (Float*)translated2);

                tempAtom->x = translated2[0];
                tempAtom->y = translated2[1];
                tempAtom->z = translated2[2];
            }
        }
    }
}

// Returns center of mass of atoms in the list
void GetCenterOfMass(Float* com, global AtomGPUsmall* atoms, global AtomGPU* ligandAtomsData, int numAtoms, int popMaxSize) {
    // Default constructor (initialise to zero)
    com[0] = PLUS_0_0f;
    com[1] = PLUS_0_0f;
    com[2] = PLUS_0_0f;
    
    // Accumulate sum of mass*coord
    Float totalMass = 0.0f;
    for(int i = 0; i < numAtoms; i++) {
        global AtomGPUsmall* tempAtom = getAtomGPUsmallFromBase(popMaxSize, i, atoms);
        Float tempMass = ligandAtomsData[i].atomicMass;
        com[0] += (tempMass * tempAtom->x);
        com[1] += (tempMass * tempAtom->y);
        com[2] += (tempMass * tempAtom->z);

        totalMass += tempMass;// GetTotalAtomicMass
    }
    // Divide by total mass
    com[0] /= totalMass;
    com[1] /= totalMass;
    com[2] /= totalMass;
}

void GetPrincipalAxes(PrincipalAxesSyncGPU* principalAxes, global AtomGPUsmall* atoms, global AtomGPU* ligandAtomsData, global Float* individual, constant parametersForGPU* parameters, int numAtoms) {
    const int N = 3;

    // TODO: No check for atomList.empty() !

    // TODO: Special case for water !

    // Construct default principal axes:
    principalAxes->com[0] = PLUS_0_0f;
    principalAxes->com[1] = PLUS_0_0f;
    principalAxes->com[2] = PLUS_0_0f;

    principalAxes->axis1[0] = PLUS_1_0f;
    principalAxes->axis1[1] = PLUS_0_0f;
    principalAxes->axis1[2] = PLUS_0_0f;

    principalAxes->axis2[0] = PLUS_0_0f;
    principalAxes->axis2[1] = PLUS_1_0f;
    principalAxes->axis2[2] = PLUS_0_0f;

    principalAxes->axis3[0] = PLUS_0_0f;
    principalAxes->axis3[1] = PLUS_0_0f;
    principalAxes->axis3[2] = PLUS_1_0f;

    principalAxes->moment1 = PLUS_1_0f;
    principalAxes->moment2 = PLUS_1_0f;
    principalAxes->moment3 = PLUS_1_0f;

    // Store center of mass of CURRENT MODEL !!!
    GetCenterOfMass((Float*)principalAxes->com, atoms, ligandAtomsData, numAtoms, parameters->popMaxSize);

    // Construct the moment of inertia tensor
    Float inertiaTensor[3*3];// cuda fix (was: float inertiaTensor[N*N])
    // Set to zero
    zerosNxN((Float*)inertiaTensor, N);
    for(int i=0; i < numAtoms; i++) {
        // Vector from center of mass to atom
        Float r[3];
        Float atomCoord[3];
        global AtomGPUsmall* tempAtom = getAtomGPUsmallFromBase(parameters->popMaxSize, i, atoms);
        atomCoord[0] = tempAtom->x;
        atomCoord[1] = tempAtom->y;
        atomCoord[2] = tempAtom->z;
        subtract2Vectors3((Float*)atomCoord, (Float*)principalAxes->com, (Float*)r);
        // Atomic mass
        Float m = ligandAtomsData[i].atomicMass;
        Float rx2 = r[0] * r[0];
        Float ry2 = r[1] * r[1];
        Float rz2 = r[2] * r[2];
        // Diagonal elements (moments of inertia)
        Float dIxx = m * (ry2 + rz2); //=r^2 - x^2
        Float dIyy = m * (rx2 + rz2); //=r^2 - y^2
        Float dIzz = m * (rx2 + ry2); //=r^2 - z^2
        inertiaTensor[0] += dIxx;// [0][0]
        inertiaTensor[4] += dIyy;// [1][1]
        inertiaTensor[8] += dIzz;// [2][2]
        // Off-diagonal elements (products of inertia) - symmetric matrix
        Float dIxy = m * r[0] * r[1];
        Float dIxz = m * r[0] * r[2];
        Float dIyz = m * r[1] * r[2];
        inertiaTensor[1] -= dIxy;// [0][1]
        inertiaTensor[3] -= dIxy;// [1][0]
        inertiaTensor[2] -= dIxz;// [0][2]
        inertiaTensor[6] -= dIxz;// [2][0]
        inertiaTensor[5] -= dIyz;// [1][2]
        inertiaTensor[7] -= dIyz;// [2][1]
    }

    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(inertiaTensor);

    // Eigen::VectorXd eigenValues = eigenSolver.eigenvalues().real();

    // Eigen::MatrixXd eigenVectors = eigenSolver.eigenvectors().real();

    Float eigenVectors[3*3];// cuda fix (was: float eigenVectors[N*N])
    Float eigenValues[3];// cuda fix (was: float eigenValues[N])
    if(dsyevh3((Float*)inertiaTensor, (Float*)eigenVectors, (Float*)eigenValues) != 0) {
        // Failure
        printf("[SyncToModel] [GetPrincipalAxes] [Failed to get eigenVectors and eigenValues]\n");
    } else {
       // Success

       // eigenvectors already normalized (length=1)

       // FIXME: assert that .imag() is Zero(N) ????? (comment from original code)

        // Load the principal axes and moments into the return parameter
        // We need to sort these so that axis1 is the first principal axis, axis2 the
        // second, axis3 the third. With only three elements to sort, this is probably
        // as good a way as any:
        unsigned int idx1 = 0;
        unsigned int idx2 = 1;
        unsigned int idx3 = 2;
        Float swap;
        if(eigenValues[idx1] > eigenValues[idx2]) {
            swap = idx1;
            idx1 = idx2;
            idx2 = swap;
        }
        if(eigenValues[idx1] > eigenValues[idx3]) {
            swap = idx1;
            idx1 = idx3;
            idx3 = swap;
        }
        if(eigenValues[idx2] > eigenValues[idx3]) {
            swap = idx2;
            idx2 = idx3;
            idx3 = swap;
        }
        principalAxes->axis1[0] = eigenVectors[idx1];// (0, idx1)
        principalAxes->axis1[1] = eigenVectors[3+idx1];// (1, idx1)
        principalAxes->axis1[2] = eigenVectors[6+idx1];// (2, idx1)

        principalAxes->axis2[0] = eigenVectors[idx2];// (0, idx2)
        principalAxes->axis2[1] = eigenVectors[3+idx2];// (1, idx2)
        principalAxes->axis2[2] = eigenVectors[6+idx2];// (2, idx2)

        principalAxes->axis3[0] = eigenVectors[idx3];// (0, idx3)
        principalAxes->axis3[1] = eigenVectors[3+idx3];// (1, idx3)
        principalAxes->axis3[2] = eigenVectors[6+idx3];// (2, idx3)

        principalAxes->moment1 = eigenValues[idx1];
        principalAxes->moment2 = eigenValues[idx2];
        principalAxes->moment3 = eigenValues[idx3];

        // DM 28 Jun 2001 - for GA crossovers in particular we need to ensure the
        // principal axes returned are always consistent for a given ligand
        // conformation. Currently the direction of the axes vectors is not controlled
        // as for e.g. (+1,0,0) is degenerate with (-1,0,0). Method is to arbitrarily
        // check the first atom and invert the axes vectors if necessary to ensure the
        // atom coords lie in the positive quadrant of the coordinate space of
        // principal axes 1 and 2 Principal axis 3 is then defined to give a
        // right-handed set of axes
        //
        // LIMITATION: If atom 1 lies exactly on PA#1 or PA#2 this check will fail.
        // Ideally we would like to test an atom on the periphery of the molecule.
        Float c0[3];
        Float firstAtomCoords[3];
        // ! Base atoms addr is atoms[0], so it is OK for atom at index 0. For atoms other than 0, [i] IS NOT OK, use getAtomGPUsmallFromBase !
        firstAtomCoords[0] = atoms[0].x;
        firstAtomCoords[1] = atoms[0].y;
        firstAtomCoords[2] = atoms[0].z;
        subtract2Vectors3((Float*)firstAtomCoords, (Float*)(principalAxes->com), (Float*)c0);
        Float d1 = dotProduct3((Float*)c0, (Float*)(principalAxes->axis1));
        Float d2 = dotProduct3((Float*)c0, (Float*)(principalAxes->axis2));
        if (d1 < PLUS_0_0f) {
            negateVector3((Float*)(principalAxes->axis1), (Float*)(principalAxes->axis1));
        }
        if (d2 < PLUS_0_0f) {
            negateVector3((Float*)(principalAxes->axis2), (Float*)(principalAxes->axis2));
        }
        crossProduct3((Float*)(principalAxes->axis1), (Float*)(principalAxes->axis2), (Float*)(principalAxes->axis3));
        
    }
    
}

// For COM and Orientation, Only Once!
void setModelValuePosition(global AtomGPUsmall* ligandAtoms, global AtomGPU* ligandAtomsData, global Float* individual, constant parametersForGPU* parameters) {
    
    // Determine the principal axes and centre of mass of the reference atoms
    PrincipalAxesSyncGPU prAxes;
        // TODO: if tethered use tetheredAtomList, not whole atomList.
    // m_refAtoms = atomList (all atoms) = ligandAtoms
    GetPrincipalAxes(&prAxes, ligandAtoms, ligandAtomsData, individual, parameters, parameters->ligandNumAtoms);

    // Determine the overall rotation required.
    // 1) Go back to realign with Cartesian axes
    PrincipalAxesSyncGPU CARTESIAN_AXES;
    CARTESIAN_AXES.com[0] = PLUS_0_0f;
    CARTESIAN_AXES.com[1] = PLUS_0_0f;
    CARTESIAN_AXES.com[2] = PLUS_0_0f;

    CARTESIAN_AXES.axis1[0] = PLUS_1_0f;
    CARTESIAN_AXES.axis1[1] = PLUS_0_0f;
    CARTESIAN_AXES.axis1[2] = PLUS_0_0f;

    CARTESIAN_AXES.axis2[0] = PLUS_0_0f;
    CARTESIAN_AXES.axis2[1] = PLUS_1_0f;
    CARTESIAN_AXES.axis2[2] = PLUS_0_0f;

    CARTESIAN_AXES.axis3[0] = PLUS_0_0f;
    CARTESIAN_AXES.axis3[1] = PLUS_0_0f;
    CARTESIAN_AXES.axis3[2] = PLUS_1_0f;

    CARTESIAN_AXES.moment1 = PLUS_1_0f;
    CARTESIAN_AXES.moment2 = PLUS_1_0f;
    CARTESIAN_AXES.moment3 = PLUS_1_0f;

    Float qBack_v[3];
    Float qBack_s;
    GetQuatFromAlignAxes(&prAxes, &CARTESIAN_AXES, &qBack_s, (Float*)qBack_v);

    // 2) Go forward to the desired orientation
    Float qForward_v[3];
    Float qForward_s;
    Float orientation[3];
    int startOrientationIndex = (parameters->chromStoreLen) - CHROM_SUBTRACT_FOR_ORIENTATION_1;
    orientation[0] = individual[startOrientationIndex];
    orientation[1] = individual[startOrientationIndex + 1];
    orientation[2] = individual[startOrientationIndex + 2];
    ToQuat((Float*)orientation, &qForward_s, (Float*)qForward_v);

    // 3 Combine the two rotations
    Float q_v[3];
    Float q_s;
    multiplyQuatResult(&qForward_s, (Float*)qForward_v, &qBack_s, (Float*)qBack_v, &q_s, (Float*)q_v);// RbtQuat q = qForward * qBack;
    
    
    Float tempCoord[3];
    Float negPrAxesCom[3];
    negateVector3((Float*)prAxes.com, (Float*)negPrAxesCom);
    Float translated[3];
    Float rotated[3];
    Float com[3];
    int startComIndex = (parameters->chromStoreLen) - CHROM_SUBTRACT_FOR_CENTER_OF_MASS_1;
    com[0] = individual[startComIndex];
    com[1] = individual[startComIndex + 1];
    com[2] = individual[startComIndex + 2];
    Float translated2[3];
    for(int i = 0; i < parameters->ligandNumAtoms; i++) {
        global AtomGPUsmall* tempAtom = getAtomGPUsmallFromBase(parameters->popMaxSize, i, ligandAtoms);

        tempCoord[0] = tempAtom->x;
        tempCoord[1] = tempAtom->y;
        tempCoord[2] = tempAtom->z;

        // Move to origin
        add2Vectors3((Float*)tempCoord, (Float*)negPrAxesCom, (Float*)translated);
        // Rotate
        RotateUsingQuat(&q_s, (Float*)q_v, (Float*)translated, (Float*)rotated);
        // Move to new centre of mass
        add2Vectors3((Float*)rotated, (Float*)com, (Float*)translated2);

        tempAtom->x = translated2[0];
        tempAtom->y = translated2[1];
        tempAtom->z = translated2[2];

    }
    
}

void syncToModel(global AtomGPUsmall* ligandAtoms, global AtomGPU* ligandAtomsData, global Float* individual, global DihedralRefDataGPU* dihedralRefData, constant parametersForGPU* parameters) {
    // TODO: setModelValueOccupancyElement(...)
    setModelValueDihedral(ligandAtoms, individual, dihedralRefData, parameters);
    setModelValuePosition(ligandAtoms, ligandAtomsData, individual, parameters);
}

#endif
