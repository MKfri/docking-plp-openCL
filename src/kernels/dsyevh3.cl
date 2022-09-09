/*
Adopted for OpenCL from:

Joachim Kopp
Numerical diagonalization of hermitian 3x3 matrices
arXiv.org preprint: physics/0610206
Int. J. Mod. Phys. C19 (2008) 523-548
*/

#ifndef DSYEVH3_CL_H
#define DSYEVH3_CL_H

#include "RealConstants.cl"

// Constants
//#define M_SQRT3    1.73205080756887729352744634151f   // sqrt(3)

// Macros
#define SQR(x)      ((x)*(x)) // x^2 

int dsyevc3(Float* A, Float* w) {
    //      float A[3][3], float w[3]
    
    Float m, c1, c0;
    // Determine coefficients of characteristic poynomial. We write
    //       | a   d   f  |
    //  A =  | d*  b   e  |
    //       | f*  e*  c  |
    Float de = A[1] * A[5];// d * e    A[0][1] * A[1][2]
    Float dd = SQR(A[1]);// d^2        A[0][1]
    Float ee = SQR(A[5]);// e^2        A[1][2]
    Float ff = SQR(A[2]);// f^2        A[0][2]
    m  = A[0] + A[4] + A[8];// A[0][0] + A[1][1] + A[2][2]

    // a*b + a*c + b*c - d^2 - e^2 - f^2
    c1 = (A[0]*A[4] + A[0]*A[8] + A[4]*A[8]) - (dd + ee + ff);
    //  (A[0][0]*A[1][1] + A[0][0]*A[2][2] + A[1][1]*A[2][2]) - (dd + ee + ff);

    // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)
    c0 = A[8]*dd + A[0]*ee + A[4]*ff - A[0]*A[4]*A[8] - PLUS_2_0f * A[2]*de;
    //  A[2][2]*dd + A[0][0]*ee + A[1][1]*ff - A[0][0]*A[1][1]*A[2][2] - 2.0 * A[0][2]*de;

    Float p, sqrt_p, q, c, s, phi;
    p = SQR(m) - PLUS_3_0f * c1;
    q = m * (p - (PLUS_3_0f / PLUS_2_0f) * c1) - (PLUS_27_0f / PLUS_2_0f) * c0;
    sqrt_p = sqrt(fabs(p));

    phi = PLUS_27_0f * (PLUS_0_25f * SQR(c1) * (p - c1) + c0 * (q + PLUS_27_0f / PLUS_4_0f * c0));
    phi = (PLUS_1_0f / PLUS_3_0f) * atan2(sqrt(fabs(phi)), q);

    c = sqrt_p * cos(phi);
    s = (PLUS_1_0f / PLUS_SQRT_3f) * sqrt_p * sin(phi);

    w[1]  = (PLUS_1_0f / PLUS_3_0f) * (m - c);
    w[2]  = w[1] + s;
    w[0]  = w[1] + c;
    w[1] -= s;

    return 0;
}

void dsytrd3(Float* A, Float* Q, Float* d, Float* e) {
    //      float A[3][3], float Q[3][3], float d[3], float e[2]
    const int n = 3;
    Float u[3], q[3];// cuda fix (was: float u[n], q[n])
    Float omega, f;
    Float K, h, g;
    int i;
    int j;

    for (i=0; i < n; i++) {
        Q[i*3+i] = PLUS_1_0f; // Q[i][i] = 1.0
        for (int j=0; j < i; j++) {
            Q[i*3+j] = Q[j*3+i] = PLUS_0_0f; // Q[i][j] = Q[j][i] = 0.0
        }
    }

    // Bring first row and column to the desired form 
    h = SQR(A[1]) + SQR(A[2]);// SQR(A[0][1]) + SQR(A[0][2])
    if (A[1] > PLUS_0_0f) {  // A[0][1] > 0
        g = -sqrt(h);
    } else {
        g = sqrt(h);
    }
    e[0] = g;
    f    = g * A[1]; // g * A[0][1]
    u[1] = A[1] - g;// A[0][1] - g
    u[2] = A[2];// A[0][2]

    omega = h - f;
    if (omega > PLUS_0_0f) {
        omega = PLUS_1_0f / omega;
        K     = PLUS_0_0f;
        for (i=1; i < n; i++) {
            f    = A[3+i] * u[1] + A[i*3+2] * u[2];// A[1][i] * u[1] + A[i][2] * u[2]
            q[i] = omega * f;// p
            K   += u[i] * f;// u* A u
        }
        K *= PLUS_0_5f * SQR(omega);

        for (i=1; i < n; i++) {
            q[i] = q[i] - K * u[i];
        }

        d[0] = A[0];// A[0][0]
        d[1] = A[4] - PLUS_2_0f*q[1]*u[1];// A[1][1] - 2.0*q[1]*u[1]
        d[2] = A[8] - PLUS_2_0f*q[2]*u[2];// A[2][2] - 2.0*q[2]*u[2]

        // Store inverse Householder transformation in Q
        for (j=1; j < n; j++) {
            f = omega * u[j];
            for (i=1; i < n; i++) {
                Q[i*3+j] = Q[i*3+j] - f*u[i];// Q[i][j] = Q[i][j] - f*u[i]
            }
        }

        // Calculate updated A[1][2] and store it in e[1]
        e[1] = A[5] - q[1]*u[2] - u[1]*q[2];// e[1] = A[1][2] - q[1]*u[2] - u[1]*q[2]
    
    } else {
        for (i=0; i < n; i++) {
            d[i] = A[i*3+i];// A[i][i]
        }
        e[1] = A[5];// A[1][2]
    }

}

int dsyevq3(Float* A, Float* Q, Float* w) {
    //      float A[3][3], float Q[3][3], float w[3]

    const int n = 3;
    Float e[3];                   // The third element is used only as temporary workspace
    Float g, r, p, f, b, s, c, t; // Intermediate storage
    int nIter;
    int m, l, i, k;

    // Transform A to real tridiagonal form by the Householder method
    dsytrd3(A, Q, w, (Float*)e);

    // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
    // with the QL method
    //
    // Loop over all off-diagonal elements
    for (l=0; l < n-1; l++) {
        nIter = 0;
        while (1) {
            // Check for convergence and exit iteration loop if off-diagonal
            // element e(l) is zero
            for (m=l; m <= n-2; m++) {
                g = fabs(w[m])+fabs(w[m+1]);
                if (fabs(e[m]) + g == g) {
                    break;
                }
            }
            if (m == l) {
                break;
            }
            if (nIter++ >= 30) {
                return -1;
            }

            // Calculate g = d_m - k
            g = (w[l+1] - w[l]) / (e[l] + e[l]);
            r = sqrt(SQR(g) + PLUS_1_0f);
            if (g > PLUS_0_0f) {
                g = w[m] - w[l] + e[l]/(g + r);
            } else {
                g = w[m] - w[l] + e[l]/(g - r);
            }

            s = c = PLUS_1_0f;
            p = PLUS_0_0f;
            for (i=m-1; i >= l; i--) {
                f = s * e[i];
                b = c * e[i];
                if (fabs(f) > fabs(g)) {
                    c      = g / f;
                    r      = sqrt(SQR(c) + PLUS_1_0f);
                    e[i+1] = f * r;
                    c     *= (s = PLUS_1_0f/r);
                } else {
                    s      = f / g;
                    r      = sqrt(SQR(s) + PLUS_1_0f);
                    e[i+1] = g * r;
                    s     *= (c = PLUS_1_0f/r);
                }

                g = w[i+1] - p;
                r = (w[i] - g)*s + PLUS_2_0f*c*b;
                p = s * r;
                w[i+1] = g + p;
                g = c*r - b;

                // Form eigenvectors
                for (k=0; k < n; k++) {
                    t = Q[k*3+(i+1)];// Q[k][i+1]
                    Q[k*3+(i+1)] = s*Q[k*3+i] + c*t;// Q[k][i+1] = s*Q[k][i] + c*t
                    Q[k*3+i]   = c*Q[k*3+i] - s*t;// Q[k][i]   = c*Q[k][i] - s*t
                }
            }

            w[l] -= p;
            e[l]  = g;
            e[m]  = PLUS_0_0f;
        }
    }

    return 0;
}


int dsyevh3(Float* A, Float* Q, Float* w) {
    //      float A[3][3], float Q[3][3], float w[3]

    Float norm;// Squared norm or inverse norm of current eigenvector
    Float error;// Estimated maximum roundoff error
    Float t, u;// Intermediate storage

    // Calculate eigenvalues
    dsyevc3(A, w);

    t = fabs(w[0]);
    if ((u=fabs(w[1])) > t) {
        t = u;
    }
    if ((u=fabs(w[2])) > t) {
        t = u;
    }
    if (t < PLUS_1_0f) {
        u = t;
    } else {
        u = SQR(t);
    }
    error = PLUS_256_0f * FLT_EPSILON * SQR(u);

    Q[1] = A[1]*A[5] - A[2]*A[4];// Q[0][1] = A[0][1]*A[1][2] - A[0][2]*A[1][1]
    Q[4] = A[2]*A[1] - A[5]*A[0];// Q[1][1] = A[0][2]*A[0][1] - A[1][2]*A[0][0]
    Q[7] = SQR(A[1]);// Q[2][1] = SQR(A[0][1])

    // Calculate first eigenvector by the formula
    //   v[0] = (A - w[0]).e1 x (A - w[0]).e2
    Q[0] = Q[1] + A[2]*w[0];// Q[0][0] = Q[0][1] + A[0][2]*w[0]
    Q[3] = Q[4] + A[5]*w[0];// Q[1][0] = Q[1][1] + A[1][2]*w[0]
    Q[6] = (A[0] - w[0]) * (A[4] - w[0]) - Q[7];// Q[2][0] = (A[0][0] - w[0]) * (A[1][1] - w[0]) - Q[2][1]
    norm = SQR(Q[0]) + SQR(Q[3]) + SQR(Q[6]);// SQR(Q[0][0]) + SQR(Q[1][0]) + SQR(Q[2][0])

    // If vectors are nearly linearly dependent, or if there might have
    // been large cancellations in the calculation of A[i][i] - w[0], fall
    // back to QL algorithm
    // Note that this simultaneously ensures that multiple eigenvalues do
    // not cause problems: If w[0] = w[1], then A - w[0] * I has rank 1,
    // i.e. all columns of A - w[0] * I are linearly dependent.
    if (norm <= error) {
        return dsyevq3(A, Q, w);
    } else { // This is the standard branch
        norm = sqrt(PLUS_1_0f / norm);

        Q[0] = Q[0] * norm; // Q[j][0] = Q[j][0] * norm     for j=0
        Q[3] = Q[3] * norm; // Q[j][0] = Q[j][0] * norm     for j=1
        Q[6] = Q[6] * norm; // Q[j][0] = Q[j][0] * norm     for j=2
    }
    // Calculate second eigenvector by the formula
    //   v[1] = (A - w[1]).e1 x (A - w[1]).e2
    Q[1] = Q[1] + A[2]*w[1];// Q[0][1]  = Q[0][1] + A[0][2]*w[1]
    Q[4] = Q[4] + A[5]*w[1];// Q[1][1]  = Q[1][1] + A[1][2]*w[1]
    Q[7] = (A[0] - w[1]) * (A[4] - w[1]) - Q[7];// Q[2][1]  = (A[0][0] - w[1]) * (A[1][1] - w[1]) - Q[2][1]
    norm = SQR(Q[1]) + SQR(Q[4]) + SQR(Q[7]);// norm = SQR(Q[0][1]) + SQR(Q[1][1]) + SQR(Q[2][1])
    if (norm <= error) {
        return dsyevq3(A, Q, w);
    } else {
        norm = sqrt(PLUS_1_0f / norm);

        Q[1] = Q[1] * norm; // Q[j][1] = Q[j][1] * norm     for j=0
        Q[4] = Q[4] * norm; // Q[j][1] = Q[j][1] * norm     for j=1
        Q[7] = Q[7] * norm; // Q[j][1] = Q[j][1] * norm     for j=2
    }
    // Calculate third eigenvector according to
    //   v[2] = v[0] x v[1]
    Q[2] = Q[3]*Q[7] - Q[6]*Q[4];// Q[0][2] = Q[1][0]*Q[2][1] - Q[2][0]*Q[1][1]
    Q[5] = Q[6]*Q[1] - Q[0]*Q[7];// Q[1][2] = Q[2][0]*Q[0][1] - Q[0][0]*Q[2][1]
    Q[8] = Q[0]*Q[4] - Q[3]*Q[1];// Q[2][2] = Q[0][0]*Q[1][1] - Q[1][0]*Q[0][1]

    return 0;
}

#endif
