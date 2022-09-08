#include <clStructs.h>
#include <constants.cl>

#include <SyncToModel.cl>

__kernel void kernelSyncToModel(constant parametersForGPU* parameters, global float* globalPopulations,
                    global AtomGPUsmall* ligandAtomsSmallGlobalAll,
                    global DihedralRefDataGPU* dihedralRefData,
                    global int* popNewIndex,
                    global AtomGPU* ligandAtoms) {

    uint runID = get_global_id(RUN_ID_2D);
    uint individualID = get_global_id(INDIVIDUAL_ID_2D);

    // What part of population to Sync (existing (only initial) or new pop) (ALL THREADS SAME PATH).
    if (popNewIndex[0] != 0) {
        individualID += parameters->popSize;
    }

    // Sync to Model Individuals
    if (individualID >= popNewIndex[0] && individualID < popNewIndex[1]) {

        global AtomGPUsmall* ligandAtomsOwn = getAtomGPUsmallBase(parameters->popMaxSize, runID, individualID, parameters->ligandNumAtoms, ligandAtomsSmallGlobalAll);
        global float* individual = getIndividual(parameters->popMaxSize, runID, individualID, parameters->chromStoreLen, globalPopulations);

        // Before we evaluate each individual, we set his score to 0.0
        // to allow splitting evaluation of SF into multiple kernels.
        individual[parameters->chromStoreLen - CHROM_SUBTRACT_FOR_SCORE] = 0.0f;

        syncToModel(ligandAtomsOwn, ligandAtoms, individual, dihedralRefData, parameters); 
    }
}
