#ifndef ROULETTE_WHEEL_SELECT_CL_H
#define ROULETTE_WHEEL_SELECT_CL_H

#include "clStructs.h"

#include "RealConstants.cl"
#include "constants.cl"

#include "tyche_i.cl"
#include "randomExtra.cl"

void selectParent(global Float* population, global Float** parent,
                 int popSize, int chromStoreLen, tyche_i_state* state) {

    // Selects one parent
    float cutoff = tyche_i_float(state[0]);
    int lower = 0;
    int upper = popSize - 1;
    int i;
    while (upper >= lower) {
        i = lower + (upper - lower) / 2;
        int currentC = (population[i * chromStoreLen + (chromStoreLen - CHROM_SUBTRACT_FOR_RW_FITNESS)] > cutoff);
        upper = currentC * (i - 1) + (currentC != 1) * upper;
        lower = currentC * lower + (currentC != 1) * (i + 1);
    }
    // make sure lower is a number between 0 and popSize - 1
    lower = fmin((Float)(popSize - 1), (Float)lower);
    lower = fmax(PLUS_0_0f, (Float)lower);

    *parent = &(population[(int)(lower * chromStoreLen)]);
}

 // NOTE: population MUST BE set at the right begining!
void RouletteWheelSelect(global Float* population, global Float** father,
                         global Float** mother, Float popSize,
                         Float chromStoreLen, tyche_i_state* state) {
    
    selectParent(population, father, popSize, chromStoreLen, state);
    selectParent(population, mother, popSize, chromStoreLen, state);
    // Check that mother and father are not the same genome
    // The check is on the pointers, not that the chromosomes are near-equal
    // If we repeatedly get the same genomes selected, this must mean
    // the population lacks diversity
    int j = 0;
    while ( (*father) == (*mother)) {
        selectParent(population, mother, popSize, chromStoreLen, state);
        if (j > 100) {
            printf("Population failure - not enough diversity\n");
        }
        j++;
    }
}

#endif
