/*  Subroutine to generate a Bezier curve.
    Copyright (c) 2000 David F. Rogers. All rights reserved.
*/

#pragma once

#include "bezier_util.h"

namespace aitn
{

/* Bezier surface subroutine */
template <typename T>
void bezsurf(T b[], int n, int m, int p1, int p2, T q[])
{
    int i;
    int j;
    int k;
    int i1;
    int j1;
    int jbas;
    int uinc;
    int winc;
    int icount;
    int ch;

    T u;
    T w;
    T jin;
    T kjm;
    T stepu;
    T stepw;
    
    icount = 0;
    stepu = 1.0/((T)(p1-1));
    stepw = 1.0/((T)(p2-1));

    u=0.0;
    
    for (uinc = 0; uinc < p1; uinc++){  /* for fixed u calculate various w's */
        if (1.0 - u < 5e-6) u=1.0; /* fix up the u = 1 value because of float */
        w = 0.0;
        for (winc = 0; winc < p2; winc++){
            if (1.0 - w < 5e-6) w=1.0; /* fix up the w = 1 value because of float */
            for (i = 0; i <= n; i++){
                jin = Basis(n,i,u); /* Bernstein basis function in the u direction (see Eq.(5.2)) */
                if (jin != 0.){ /* don't bother no contribution */
                    jbas = 3 * (m + 1) * i; /* column index for lineal array*/
                    for (j = 0; j <= m; j++){
                        kjm = Basis(m,j,w);  /* Bernstein basis function in the w direction (see Eq.(5.2)) */
                        if (kjm !=0.){ /* don't bother no contribution */
                            j1 = jbas + 3 * j;
                            q[icount] = q[icount]+b[j1]*jin*kjm; /* calculate the surface points */
                            q[icount+1] = q[icount+1]+b[j1+1]*jin*kjm;
                            q[icount+2] = q[icount+2]+b[j1+2]*jin*kjm;
                        }
                    }
                }
            }
            icount = icount + 3;
            w = w + stepw;
        }
        u = u + stepu;
    }
}

}