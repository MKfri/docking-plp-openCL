#include <clStructs.h>
#include <constants.cl>

#include <PLP.cl>
#include <grid.cl>

__kernel void kernelSfClash(constant parametersForGPU* parameters,
                    global AtomGPUsmall* ligandAtomsSmallGlobalAll,
                    global float* globalPopulations,
                    global int* popNewIndex,
                    global LigandAtomPairsForClash* ligandAtomPairsForClash) {

    uint runID = get_global_id(RUN_ID_2D);
    uint individualID = get_global_id(INDIVIDUAL_ID_2D);

    float clashScore = 0.0f;

    // What part of population to Score (existing (only initial) or new pop) (ALL THREADS SAME PATH).
    if (popNewIndex[0] != 0) {
        individualID += parameters->popSize;
    }

    // Score Individuals
    if (individualID >= popNewIndex[0] && individualID < popNewIndex[1]) {

        global AtomGPUsmall* ligandAtomsOwn = getAtomGPUsmallBase(parameters->popMaxSize, runID, individualID, parameters->ligandNumAtoms, ligandAtomsSmallGlobalAll);

        // f_clash term:
        int currentAtom = 0;
        for (int i = 0; i < parameters->numUniqueAtom1PairsForClash; i++) {

            int numAtom1 = ligandAtomPairsForClash[currentAtom].numAtom1;

            global AtomGPUsmall* atom1 = getAtomGPUsmallFromBase(parameters->popMaxSize, ligandAtomPairsForClash[currentAtom].atom1ID - 1, ligandAtomsOwn);

            float atom1COORD[3];
            atom1COORD[0] = atom1->x;
            atom1COORD[1] = atom1->y;
            atom1COORD[2] = atom1->z;

            for (int a = currentAtom; a < currentAtom + numAtom1; a++) {
                global AtomGPUsmall* atom2 = getAtomGPUsmallFromBase(parameters->popMaxSize, ligandAtomPairsForClash[a].atom2ID - 1, ligandAtomsOwn);
                clashScore += PLPclash((float*)atom1COORD, atom2, ligandAtomPairsForClash[a].rClash);
            }
            currentAtom += numAtom1;
        }
        global float* individual = getIndividual(parameters->popMaxSize,
                                                runID, individualID,
                                                parameters->chromStoreLen,
                                                globalPopulations);

        // Set score:
        individual[parameters->chromStoreLen - CHROM_SUBTRACT_FOR_SCORE] += clashScore;
    }
}
