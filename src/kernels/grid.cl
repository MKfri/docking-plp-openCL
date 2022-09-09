#ifndef GRID_CL_H
#define GRID_CL_H

#include "clStructs.h"

#include "RealConstants.cl"
#include "PLPConstants.cl"
#include "vector3.cl"

#define GET_GRID(gi) grid[(check>0)*gi].point[classification]

inline int index3Dto1D(int* xyz, GridGPU* ownGrid) {

    return xyz[0] * ownGrid->SXYZ[0] + xyz[1] * ownGrid->SXYZ[1] + xyz[2];
}

inline int index3Dto1Dxyz(int x, int y, int z, GridGPU* ownGrid) {

    return x * ownGrid->SXYZ[0] + y * ownGrid->SXYZ[1] + z;
}

inline void index1Dto3D(int i, GridGPU* ownGrid, int* xyz) {

    xyz[0] = i / ownGrid->SXYZ[0];
    i %= ownGrid->SXYZ[0];
    xyz[1] = i / ownGrid->SXYZ[1];
    xyz[2] = i % ownGrid->SXYZ[1];
}

inline void index3DtoCoords(int* xyz, GridGPU* ownGrid, Float* coord) {

    coord[0] = ownGrid->gridMin[0] + xyz[0] * ownGrid->gridStep[0];
    coord[1] = ownGrid->gridMin[1] + xyz[1] * ownGrid->gridStep[1];
    coord[2] = ownGrid->gridMin[2] + xyz[2] * ownGrid->gridStep[2];
}


// DM 20 Jul 2000 - get values smoothed by trilinear interpolation
// D. Oberlin and H.A. Scheraga, J. Comp. Chem. (1998) 19, 71.
Float GetSmoothedValue(global AtomGPUsmall* lAtom, global AtomGPU* atom, Float* min, Float* max, constant parametersForGPU* parameters, GridGPU* ownGrid, global GridPointGPU* grid) {

    // ligand atom coordinates:
    Float c[3];
    c[0] = lAtom->x;
    c[1] = lAtom->y;
    c[2] = lAtom->z;

    // c_site term:
    int check = (c[0] > min[0] && c[0] < max[0] &&
                 c[1] > min[1] && c[1] < max[1] &&
                 c[2] > min[2] && c[2] < max[2]);

    int heavy = (atom->triposType != TRIPOS_TYPE_H && atom->triposType != TRIPOS_TYPE_H_P);

    Float rx = PLUS_1_0f / ownGrid->gridStep[0]; // reciprocal of grid step (x)
    Float ry = PLUS_1_0f / ownGrid->gridStep[1]; // reciprocal of grid step (y)
    Float rz = PLUS_1_0f / ownGrid->gridStep[2]; // reciprocal of grid step (z)

    // Get lower left corner grid point
    //(not necessarily the nearest grid point as returned by GetIX() etc)
    // Need to shift the int(..) argument by half a grid step
    Float diff[3];
    subtract2Vectors3((Float*)c, (Float*)(ownGrid->gridMin), (Float*)diff);
    int iXYZ[3];
    iXYZ[0] = (int)((unsigned int)(rx * diff[0] - PLUS_0_5f));
    iXYZ[1] = (int)((unsigned int)(ry * diff[1] - PLUS_0_5f));
    iXYZ[2] = (int)((unsigned int)(rz * diff[2] - PLUS_0_5f));

    // p is the vector relative to the lower left corner
    Float cXYZ[3];
    Float p[3];
    index3DtoCoords((int*)iXYZ, ownGrid, (Float*)cXYZ);
    subtract2Vectors3((Float*)c, (Float*)(cXYZ), (Float*)p);
    
    Float val = PLUS_0_0f;
    // DM 3/5/2005 - fully unroll this loop
    Float bx1 = rx * p[0];
    Float bx0 = PLUS_1_0f - bx1;
    Float by1 = ry * p[1];
    Float by0 = PLUS_1_0f - by1;
    Float bz1 = rz * p[2];
    Float bz0 = PLUS_1_0f - bz1;
    Float bx0by0 = bx0 * by0;
    Float bx0by1 = bx0 * by1;
    Float bx1by0 = bx1 * by0;
    Float bx1by1 = bx1 * by1;

    int classification = atom->classification; // Used by GET_GRID macro

    val += GET_GRID(index3Dto1Dxyz(iXYZ[0], iXYZ[1], iXYZ[2], ownGrid)) * bx0by0 * bz0;
    val += GET_GRID(index3Dto1Dxyz(iXYZ[0], iXYZ[1], iXYZ[2] + 1, ownGrid)) * bx0by0 * bz1;
    val += GET_GRID(index3Dto1Dxyz(iXYZ[0], iXYZ[1] + 1, iXYZ[2], ownGrid)) * bx0by1 * bz0;
    val += GET_GRID(index3Dto1Dxyz(iXYZ[0], iXYZ[1] + 1, iXYZ[2] + 1, ownGrid)) * bx0by1 * bz1;
    val += GET_GRID(index3Dto1Dxyz(iXYZ[0] + 1, iXYZ[1], iXYZ[2], ownGrid)) * bx1by0 * bz0;
    val += GET_GRID(index3Dto1Dxyz(iXYZ[0] + 1, iXYZ[1], iXYZ[2] + 1, ownGrid)) * bx1by0 * bz1;
    val += GET_GRID(index3Dto1Dxyz(iXYZ[0] + 1, iXYZ[1] + 1, iXYZ[2], ownGrid)) * bx1by1 * bz0;
    val += GET_GRID(index3Dto1Dxyz(iXYZ[0] + 1, iXYZ[1] + 1, iXYZ[2] + 1, ownGrid)) * bx1by1 * bz1;

    // FIXME: does not calculate val for a non heavy atom outside cavity as it should (currently returns 0.0f)
    // c_site:
    //keep val:       inside cavity    |  outside and heavy                          else: 0.0f        
    val = ((Float)(check > 0)) * val + ((Float)(check == 0 && heavy > 0)) * val;

    return val;
}


#endif
