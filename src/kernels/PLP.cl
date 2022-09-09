/*
From:
 J. Chem. Inf. Model. 2009, 49, 84–96
Empirical Scoring Functions for Advanced Protein-Ligand Docking with PLANTS
Oliver Korb,† Thomas Stützle,‡ and Thomas E. Exner*,†
Theoretische Chemische Dynamik, Fachbereich Chemie, Universität Konstanz, 78457 Konstanz, Germany,
and IRIDIA, CoDE, Université Libre de Bruxelles, Brussels, Belgium

Using PLANTS Model5
*/

#ifndef PLP_CL_H
#define PLP_CL_H

#include "clStructs.h"

#include "RealConstants.cl"
#include "PLPConstants.cl"

#include "classification.cl"
#include "vector3.cl"

inline Float plp(Float r, Float A, Float B, Float C, Float D, Float E, Float F) {
    
    return (r < A) * (F * (A - r) / A)
            + (A <= r && r < B) * (E * (r - A) / (B - A))
            + (B <= r && r < C) * (E)
            + (C <= r && r <= D) * (E * (D - r) / (D - C));
}

inline Float rep(Float r, Float A, Float B, Float C, Float D) {

    return (r < A) * (r * (C - D) / A + D)
            + (A <= r && r <= B) * (-C * (r - A) / (B - A) + C);
}

Float FplpActual(Float* v1, int classification, global AtomGPU* rAtom) {

    int interaction = getPLPinteraction(classification, rAtom->classification);
    // receptor atom coordinates:
    Float v2[3];
    v2[0] = rAtom->x;
    v2[1] = rAtom->y;
    v2[2] = rAtom->z;
    // Distance
    Float r = distance2Points3((Float*)v1, (Float*)v2);

    return (interaction == PLP_INTERACTION_hbond ||
            interaction == PLP_INTERACTION_metal ||
            interaction == PLP_INTERACTION_buried ||
            interaction == PLP_INTERACTION_nonpolar)
                *
            plp(r, FPLP_TABLE[interaction * 6], FPLP_TABLE[interaction * 6 + 1],
                   FPLP_TABLE[interaction * 6 + 2], FPLP_TABLE[interaction * 6 + 3],
                   FPLP_TABLE[interaction * 6 + 4], PLUS_20_0f)
            + 
            (interaction == PLP_INTERACTION_repulsive)
                *
            rep(r, PLUS_3_2f, PLUS_5_0f, W_PLP_REP_M5 * PLUS_0_1f, W_PLP_REP_M5 * PLUS_20_0f);
}

Float Fplp(global AtomGPUsmall* lAtom, int classification, global AtomGPU* rAtom) {

    // ligand atom coordinates:
    Float v1[3];
    v1[0] = lAtom->x;
    v1[1] = lAtom->y;
    v1[2] = lAtom->z;

    return FplpActual((Float*)v1, classification, rAtom);
}

Float CsiteAtom(global AtomGPUsmall* lAtom, int triposType, constant parametersForGPU* parameters) {

    Float atom[3];
    atom[0] = lAtom->x;
    atom[1] = lAtom->y;
    atom[2] = lAtom->z;
    Float min[3];
    min[0] = parameters->dockingSiteInfo.minCavity.x;
    min[1] = parameters->dockingSiteInfo.minCavity.y;
    min[2] = parameters->dockingSiteInfo.minCavity.z;
    Float max[3];
    max[0] = parameters->dockingSiteInfo.maxCavity.x;
    max[1] = parameters->dockingSiteInfo.maxCavity.y;
    max[2] = parameters->dockingSiteInfo.maxCavity.z;

    // if heavy atom outside, constant penalty
    return (triposType != TRIPOS_TYPE_H && triposType != TRIPOS_TYPE_H_P && 
            (atom[0] < min[0] || atom[0] > max[0] ||
             atom[1] < min[1] || atom[1] > max[1] ||
             atom[2] < min[2] || atom[2] > max[2])) * PLUS_50_0f;
}

inline int CsiteReceptor(Float x, Float y, Float z, constant parametersForGPU* parameters) {
    Float min[3];
    min[0] = parameters->dockingSiteInfo.minReceptor.x;
    min[1] = parameters->dockingSiteInfo.minReceptor.y;
    min[2] = parameters->dockingSiteInfo.minReceptor.z;
    Float max[3];
    max[0] = parameters->dockingSiteInfo.maxReceptor.x;
    max[1] = parameters->dockingSiteInfo.maxReceptor.y;
    max[2] = parameters->dockingSiteInfo.maxReceptor.z;
    
    return (x > min[0] && x < max[0] && y > min[1] && y < max[1] && z > min[2] && z < max[2]);
}

Float PLPclash(Float* atom1COORD, global AtomGPUsmall* atom2, Float rClash) {

    Float atom2COORD[3];
    atom2COORD[0] = atom2->x;
    atom2COORD[1] = atom2->y;
    atom2COORD[2] = atom2->z;

    Float r = distance2PointsSq3(atom1COORD, (Float*)atom2COORD);
    // r and rClash are squared
    return ( (Float)( r < rClash) * ( ( (rClash - r) / (rClash) ) * PLUS_50_0f ));
}

/*
From:
Clark, Matthew, Richard D. Cramer III, and Nicole Van Opdenbosch.
"Validation of the general purpose tripos 5.2 force field."
Journal of computational chemistry 10.8 (1989): 982-1012.
*/
inline Float PLPtorsional(Float angle, Float k, Float s) {

    return k * (PLUS_1_0f + s/fabs(s) * cos(fabs(s) * angle));
}

#endif
