#include "clStructs.h"

#include "RealConstants.cl"
#include "constants.cl"

__kernel void kernelNormalize(constant parametersForGPU* parameters, global Float* globalPopulations,
                    local Float* localScore, global Float* bestScore) {

    uint runID=get_global_id(RUN_ID_2D);

    uint localSize=get_local_size(INDIVIDUAL_ID_2D);
    uint localID=get_local_id(INDIVIDUAL_ID_2D);

    // Normalize

    // Set start and end index
    int size = parameters->popSize;
    int perThread = size / localSize + 1;// +1 to prevent, when localSize=256 and size=500, last thread does almost half, but other problems

    int startIndex = localID * perThread;
    int endIndex = startIndex + perThread;
    
    // Set Individual score:
    int popMaxSize = parameters->popMaxSize;
    int chromStoreLen = parameters->chromStoreLen;

    Float sum = PLUS_0_0f;
    Float sumSq = PLUS_0_0f;
    Float tempScore;

    for(int i=startIndex; i < endIndex && i < size; i++) {
        
        global Float* individual = getIndividual(popMaxSize, runID, i, chromStoreLen, globalPopulations);
        tempScore = individual[chromStoreLen - CHROM_SUBTRACT_FOR_SCORE];

        sum += tempScore;
        sumSq += tempScore * tempScore;
    }

    // length localScore = 2 * localSize    [ sum | sumSq ]
    localScore[localID] = sum;
    localScore[localID + localSize] = sumSq;
    
    // Reduce Local,  Adapted from: https://web.engr.oregonstate.edu/~mjb/cs575/Handouts/opencl.reduction.2pp.pdf
    for (int offset = 1; offset < localSize; offset *= 2) {
        int mask = 2 * offset - 1;
        barrier(CLK_LOCAL_MEM_FENCE);
        if ((localID & mask) == 0) {
            localScore[localID] += localScore[localID + offset];// sum
            localScore[localID + localSize] += localScore[localID + localSize + offset];// sumSq
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    sum = localScore[0];
    sumSq = localScore[localSize];

    barrier(CLK_LOCAL_MEM_FENCE);

    Float m_scoreMean = sum / size;
    Float m_scoreVariance = (sumSq / size) - (m_scoreMean * m_scoreMean);
    Float sigma = sqrt(m_scoreVariance);
    // calculate scaled fitness values using sigma truncation
    // Goldberg page 124
    Float m_c = PLUS_2_0f;// Sigma Truncation Multiplier
    Float offset = m_scoreMean - m_c * sigma;
    Float partialSum = PLUS_0_0f;
    
    // First set the unnormalised fitness values
    Float tempScoreRW;
    // FIXME: this is sequencial code...
    if (localID == 0) {

        for (int i = 0; i < size; i++) {
            
            global Float* individual = getIndividual(popMaxSize, runID, i, chromStoreLen, globalPopulations);

            tempScore = individual[chromStoreLen - CHROM_SUBTRACT_FOR_SCORE];
            tempScoreRW = fmax(PLUS_0_0f, tempScore - offset);

            tempScoreRW += partialSum;

            individual[chromStoreLen - CHROM_SUBTRACT_FOR_RW_FITNESS] = tempScoreRW;

            partialSum = tempScoreRW;
        }
        localScore[0] = partialSum;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    //NormaliseRWFitness
    Float total = localScore[0];

    for(int i=startIndex; i < endIndex && i < size; i++) {
        
        global Float* individual = getIndividual(popMaxSize, runID, i, chromStoreLen, globalPopulations);
        individual[chromStoreLen - CHROM_SUBTRACT_FOR_RW_FITNESS] /= total;

        // set best score
        if(i == size - 1) {
            bestScore[runID] = individual[chromStoreLen - CHROM_SUBTRACT_FOR_SCORE];
        }
    }
}
