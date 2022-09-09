#include "clStructs.h"

#include "RealConstants.cl"
#include "constants.cl"

__kernel void kernelInit2(constant parametersForGPU* parameters,
                    global AtomGPU* ligandAtoms,
                    global AtomGPUsmall* ligandAtomsSmallGlobalAll) {

    uint runID=get_global_id(RUN_ID_2D);
    uint individualID=get_global_id(INDIVIDUAL_ID_2D);

    if(individualID < parameters->popMaxSize) {

        // Copy Own Ligand Atom Instance
        global AtomGPUsmall* ligandAtomsOwn = getAtomGPUsmallBase(parameters->popMaxSize, runID, individualID, parameters->ligandNumAtoms, ligandAtomsSmallGlobalAll);

        for(int i = 0; i < parameters->ligandNumAtoms; i++) {
            global AtomGPUsmall* tempAtom = getAtomGPUsmallFromBase(parameters->popMaxSize, i, ligandAtomsOwn);
            tempAtom->x = ligandAtoms[i].x;
            tempAtom->y = ligandAtoms[i].y;
            tempAtom->z = ligandAtoms[i].z;
        }
    }
}
