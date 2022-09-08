#include <clStructs.h>
#include <constants.cl>

#include <PLP.cl>
#include <grid.cl>

__kernel void kernelSfTors(constant parametersForGPU* parameters,
                    global float* globalPopulations,
                    global int* popNewIndex,
                    global DihedralRefDataGPU* dihedralRefData) {

    uint runID = get_global_id(RUN_ID_2D);
    uint individualID = get_global_id(INDIVIDUAL_ID_2D);

    float torsionalScore = 0.0f;

    // What part of population to Score (existing (only initial) or new pop) (ALL THREADS SAME PATH).
    if (popNewIndex[0] != 0) {
        individualID += parameters->popSize;
    }

    // Score Individuals
    if (individualID >= popNewIndex[0] && individualID < popNewIndex[1]) {
        global float* individual = getIndividual(parameters->popMaxSize,
                                                runID, individualID,
                                                parameters->chromStoreLen,
                                                globalPopulations);
        
        // f_tors term:
        for (int i = 0; i < parameters->numDihedralElements; i++) {
            torsionalScore += PLPtorsional(individual[i], dihedralRefData[i].k, dihedralRefData[i].s);
        }

        // Set score:
        individual[parameters->chromStoreLen - CHROM_SUBTRACT_FOR_SCORE] += torsionalScore;
    }
}
