#ifndef RANDOM_EXTRA_CL_H
#define RANDOM_EXTRA_CL_H

#include "RealConstants.cl"
#include "tyche_i.cl"

int GetRandomInt(int nMax, tyche_i_state* state) {
    //get a random integer between 0 and nMax-1
    int randInt = (int) (tyche_i_float(state[0]) * (Float)nMax);

    int check = (randInt >= nMax);
    return (check != 1) * randInt + check * (nMax - 1);
}

Float PositiveCauchyRandomWMeanVariance(Float mean, Float variance, tyche_i_state* state) {

    Float t = tan(PLUS_PI * (tyche_i_double(state[0]) - PLUS_0_5f));

    t = (t <= PLUS_100_0f) * t + (t > PLUS_100_0f) * PLUS_100_0f;
    t = (t >= MINUS_100_0f) * t + (t < MINUS_100_0f) * MINUS_100_0f;

    return fabs(mean + variance * t);
}

#endif
