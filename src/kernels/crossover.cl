#ifndef CROSSOVER_CL_H
#define CROSSOVER_CL_H

#include "clStructs.h"

#include "RealConstants.cl"
#include "constants.cl"

#include "tyche_i.cl"
#include "randomExtra.cl"

void Crossover(global Float* parent1, global Float* parent2,
               global Float* child1, global Float* child2,
               int chromStoreLen, tyche_i_state* state) {
    // No check for same crossover length!

    // length is numOfRotatableBonds + COM=1 + Orientation=1
    int length = chromStoreLen - CHROM_SUBTRACT_FOR_ROTATABLE_BONDS_LENGTH + 2;
    
    // 2-point crossover
    // In the spirit of STL, ixbegin is the first gene to crossover, ixend is one
    // after the last gene to crossover
    int ixbegin = GetRandomInt(length, state);
    
    // if ixbegin is 0, we need to avoid selecting the whole chromosome
    int ixend = (ixbegin == 0) * (GetRandomInt(length - 1, state) + 1) + (ixbegin != 0) * (GetRandomInt(length - ixbegin, state) + ixbegin + 1);
    
    // Get actual indexes:
    int arr[3];
    arr[0] = ixbegin;// Start at Dihedral
    arr[1] = chromStoreLen - CHROM_SUBTRACT_FOR_ORIENTATION_1;// Start at Orientation
    arr[2] = chromStoreLen - CHROM_SUBTRACT_FOR_CENTER_OF_MASS_1;// Start at COM
  
    int start = arr[(ixbegin == length - 1) + 2 * (ixbegin == length - 2)];
    
    arr[0] = ixend;// End at Dihedral
    arr[1] = chromStoreLen - CHROM_SUBTRACT_FOR_ORIENTATION_3 + 1;// End at Orientation
    arr[2] = chromStoreLen - CHROM_SUBTRACT_FOR_CENTER_OF_MASS_3 + 1;// End at COM

    int end = arr[(ixend == length) + 2 * (ixend == length - 1)];

    // Swaps each element in the range [start,end) or Copy
    for(int i = 0; i < chromStoreLen - CHROM_SUBTRACT_FOR_CHROM_ACTUAL_LENGTH; i++) {
        
        int check = (i >= start && i < end);

        Float p1 = parent1[i];
        Float p2 = parent2[i];
        child1[i] = ((Float)check) * p2 + ((Float)(check != 1)) * p1;
        child2[i] = ((Float)check) * p1 + ((Float)(check != 1)) * p2;
    }
}

void CopyParent2ToChild2(global Float* parent1, global Float* parent2,
                       global Float* child1, global Float* child2,
                       int chromStoreLen) {

    for(int i = 0; i < chromStoreLen - CHROM_SUBTRACT_FOR_CHROM_ACTUAL_LENGTH; i++) {
        // Copy parent to child
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }
}

#endif
