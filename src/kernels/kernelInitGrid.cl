#include "clStructs.h"

#include "RealConstants.cl"
#include "constants.cl"

#include "grid.cl"
#include "PLP.cl"

__kernel void kernelInitGrid(constant parametersForGPU* parameters,
                    global GridPointGPU* grid,
                    global AtomGPU* receptorAtoms,
                    global int* numGoodReceptors,// one int
                    global int* receptorIndex) {

    uint classID=get_global_id(CLASS_ID_2D);
    uint xyzID=get_global_id(XYZ_ID_2D);

    if(xyzID < parameters->grid.N) {

        GridGPU ownGrid;

        ownGrid.gridStep[0] = parameters->grid.gridStep[0];
        ownGrid.gridStep[1] = parameters->grid.gridStep[1];
        ownGrid.gridStep[2] = parameters->grid.gridStep[2];

        ownGrid.gridMin[0] = parameters->grid.gridMin[0];
        ownGrid.gridMin[1] = parameters->grid.gridMin[1];
        ownGrid.gridMin[2] = parameters->grid.gridMin[2];

        ownGrid.SXYZ[0] = parameters->grid.SXYZ[0];
        ownGrid.SXYZ[1] = parameters->grid.SXYZ[1];
        ownGrid.SXYZ[2] = parameters->grid.SXYZ[2];

        int xyz[3];
        index1Dto3D(xyzID, &ownGrid, (int*)xyz);

        Float coord[3];
        index3DtoCoords((int*)xyz, &ownGrid, (Float*)coord);
       
        Float score = PLUS_0_0f;

        for(int i = 0; i < (*numGoodReceptors); i++) {

            score += FplpActual((Float*)coord, classID, &(receptorAtoms[receptorIndex[i]]));
        }

        grid[xyzID].point[classID] = score;

        // Set at index=0 c_site penalty, for coords outside docking site.
        if(xyzID == 0) {
            grid[xyzID].point[classID] = PLUS_50_0f;
        }
    }
}
