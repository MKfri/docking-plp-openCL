#ifndef REAL_CONSTANTS_CL_H
#define REAL_CONSTANTS_CL_H



// Za compilation z Clang:
// clang -v -Xclang -finclude-default-header <kernel.cl>
// -> predvsem za lovljenje bugov s predprocesorjem



#ifdef USE_DOUBLE
    // Use double as floating point precision (8 Bytes)
    #define Float double

    // Set this macro for AMD cards
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable


    #define PLUS_0_0f 0.0
    #define PLUS_SQRT_3f 1.73205080756887729352744634151
    #define PLUS_PI M_PI // pi -> double

    #define PLUS_0_001f 0.001
    #define PLUS_0_04f 0.04
    #define PLUS_0_1f 0.1
    #define PLUS_0_12f 0.12
    #define PLUS_0_126f 0.126
    #define PLUS_0_2f 0.2
    #define PLUS_0_25f 0.25
    #define PLUS_0_273f 0.273
    #define PLUS_0_274f 0.274
    #define PLUS_0_32f 0.32
    #define PLUS_0_4f 0.4
    #define PLUS_0_499999f 0.499999
    #define PLUS_0_5f 0.5
    #define PLUS_0_6f 0.6
    #define PLUS_0_7f 0.7
    #define PLUS_0_999f 0.999

    #define PLUS_1_0f 1.0
    #define PLUS_1_2f 1.2
    #define PLUS_1_4f 1.4
    #define PLUS_1_424f 1.424
    #define PLUS_1_5f 1.5
    #define PLUS_1_6f 1.6
    #define PLUS_2_0f 2.0
    #define PLUS_2_2f 2.2
    #define PLUS_2_3f 2.3
    #define PLUS_2_5f 2.5
    #define PLUS_2_6f 2.6
    #define PLUS_2_75f 2.75
    #define PLUS_2_8f 2.8
    #define PLUS_3_0f 3.0
    #define PLUS_3_1f 3.1
    #define PLUS_3_2f 3.2
    #define PLUS_3_4f 3.4
    #define PLUS_3_6f 3.6
    #define PLUS_4_0f 4.0
    #define PLUS_4_5f 4.5
    #define PLUS_5_0f 5.0
    #define PLUS_5_5f 5.5
    #define PLUS_5_8f 5.8
    #define PLUS_6_46f 6.46

    #define PLUS_12_0f 12.0
    #define PLUS_12_5f 12.5
    #define PLUS_20_0f 20.0
    #define PLUS_27_0f 27.0
    #define PLUS_50_0f 50.0
    #define PLUS_100_0f 100.0
    #define PLUS_180_0f 180.0
    #define PLUS_222_222f 222.222
    #define PLUS_256_0f 256.0
    #define PLUS_333_333f 333.333
    #define PLUS_360_0f 360.0
    #define PLUS_444_444f 444.444
    #define PLUS_555_555f 555.555
    #define PLUS_666_666f 666.666
    #define PLUS_777_777f 777.777
    #define PLUS_888_888f 888.888
    #define PLUS_999_999f 999.999
    #define PLUS_99999_9f 99999.9
    #define PLUS_999999_9f 999999.9


    #define MINUS_0_05f -0.05
    #define MINUS_0_4f -0.4
    #define MINUS_0_499999f -0.499999
    #define MINUS_1_0f -1.0
    #define MINUS_2_0f -2.0
    #define MINUS_3_0f -3.0
    #define MINUS_4_0f -4.0
    #define MINUS_100_0f -100.0
    #define MINUS_180_0f -180.0

#else
    // Use float as floating point precision (4 Bytes)
    #define Float float


    #define PLUS_0_0f 0.0f
    #define PLUS_SQRT_3f 1.73205080756887729352744634151f
    #define PLUS_PI M_PI_F // pi -> float

    #define PLUS_0_001f 0.001f
    #define PLUS_0_04f 0.04f
    #define PLUS_0_1f 0.1f
    #define PLUS_0_12f 0.12f
    #define PLUS_0_126f 0.126f
    #define PLUS_0_2f 0.2f
    #define PLUS_0_25f 0.25f
    #define PLUS_0_273f 0.273f
    #define PLUS_0_274f 0.274f
    #define PLUS_0_32f 0.32f
    #define PLUS_0_4f 0.4f
    #define PLUS_0_499999f 0.499999f
    #define PLUS_0_5f 0.5f
    #define PLUS_0_6f 0.6f
    #define PLUS_0_7f 0.7f
    #define PLUS_0_999f 0.999f

    #define PLUS_1_0f 1.0f
    #define PLUS_1_2f 1.2f
    #define PLUS_1_4f 1.4f
    #define PLUS_1_424f 1.424f
    #define PLUS_1_5f 1.5f
    #define PLUS_1_6f 1.6f
    #define PLUS_2_0f 2.0f
    #define PLUS_2_2f 2.2f
    #define PLUS_2_3f 2.3f
    #define PLUS_2_5f 2.5f
    #define PLUS_2_6f 2.6f
    #define PLUS_2_75f 2.75f
    #define PLUS_2_8f 2.8f
    #define PLUS_3_0f 3.0f
    #define PLUS_3_1f 3.1f
    #define PLUS_3_2f 3.2f
    #define PLUS_3_4f 3.4f
    #define PLUS_3_6f 3.6f
    #define PLUS_4_0f 4.0f
    #define PLUS_5_0f 5.0f
    #define PLUS_4_5f 4.5f
    #define PLUS_5_5f 5.5f
    #define PLUS_5_8f 5.8f
    #define PLUS_6_46f 6.46f

    #define PLUS_12_0f 12.0f
    #define PLUS_12_5f 12.5f
    #define PLUS_20_0f 20.0f
    #define PLUS_27_0f 27.0f
    #define PLUS_50_0f 50.0f
    #define PLUS_100_0f 100.0f
    #define PLUS_180_0f 180.0f
    #define PLUS_222_222f 222.222f
    #define PLUS_256_0f 256.0f
    #define PLUS_333_333f 333.333f
    #define PLUS_360_0f 360.0f
    #define PLUS_444_444f 444.444f
    #define PLUS_555_555f 555.555f
    #define PLUS_666_666f 666.666f
    #define PLUS_777_777f 777.777f
    #define PLUS_888_888f 888.888f
    #define PLUS_999_999f 999.999f
    #define PLUS_99999_9f 99999.9f
    #define PLUS_999999_9f 999999.9f


    #define MINUS_0_05f -0.05f
    #define MINUS_0_4f -0.4f
    #define MINUS_0_499999f -0.499999f
    #define MINUS_1_0f -1.0f
    #define MINUS_2_0f -2.0f
    #define MINUS_3_0f -3.0f
    #define MINUS_4_0f -4.0f
    #define MINUS_100_0f -100.0f
    #define MINUS_180_0f -180.0f
#endif

#endif