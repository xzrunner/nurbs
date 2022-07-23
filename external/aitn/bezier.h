/*  Subroutine to generate a Bezier curve.
    Copyright (c) 2000 David F. Rogers. All rights reserved.

    b[]        = array containing the defining polygon vertices
                  b[1] contains the x-component of the vertex
                  b[2] contains the y-component of the vertex
                  b[3] contains the z-component of the vertex
    Basis      = function to calculate the Bernstein basis value (see MECG Eq 5-65)
    cpts       = number of points to be calculated on the curve
    Fractrl    = function to calculate the factorial of a number
    j[]        = array containing the basis functions for a single value of t
    npts       = number of defining polygon vertices
    p[]        = array containing the curve points
                 p[1] contains the x-component of the point
                 p[2] contains the y-component of the point
                 p[3] contains the z-component of the point
    t          = parameter value 0 <= t <= 1
*/

#pragma once

#include <math.h>
#include <stdio.h>

#include <assert.h>

namespace nlib
{

/* function to calculate the factorial */

template <typename T>
T factrl(int n)
{
    static int ntop=6;
    static T a[33]={1.0,1.0,2.0,6.0,24.0,120.0,720.0}; /* fill in the first few values */
    int j1;

    if (n < 0) printf("\nNegative factorial in routine FACTRL\n");
    if (n > 32) printf("\nFactorial value too large in routine FACTRL\n");

    while (ntop < n) { /* use the precalulated value for n = 0....6 */
        j1 = ntop++;
        a[n]=a[j1]*ntop;
    }
    return a[n]; /* returns the value n! as a floating point number */
}

/* function to calculate the factorial function for Bernstein basis */

template <typename T>
T Ni(int n,int i)
{
    T ni;
    ni = factrl<T>(n)/(factrl<T>(i)*factrl<T>(n-i));
    return ni;
}

/* function to calculate the Bernstein basis */

template <typename T>
T Basis(int n,int i,T t)
{
    T basis;
    T ti; /* this is t^i */
    T tni; /* this is (1 - t)^i */

    /* handle the special cases to avoid domain problem with pow */

    if (t==0. && i == 0) ti=static_cast<T>(1.0); else ti = static_cast<T>(pow(t,i));
    if (n==i && t==1.) tni= static_cast<T>(1.0); else tni = static_cast<T>(pow((1-t),(n-i)));
    basis = Ni<T>(n,i)*ti*tni; /* calculate Bernstein basis function */
    return basis;
}

/* Bezier curve subroutine */

template <typename T, int N>
void bezier(int npts, const T b[], int cpts, T p[])
{
    int i;
    int j;
    int i1;
    int icount;
    int jcount;

    T step;
    T t;

/*    calculate the points on the Bezier curve */

    icount = 0;
    t = 0;
    step = static_cast<T>(1.0)/((T)(cpts -1));

    for (i1 = 0; i1<cpts; i1++){ /* main loop */

        if ((1.0 - t) < 5e-6) t = 1.0;

        for (j = 0; j < N; j++){ /* generate a point on the curve */
            jcount = j;
            p[icount+j] = 0.;
            for (i = 0; i < npts; i++){ /* Do x,y,z components */
                p[icount + j] = p[icount + j] + Basis(npts-1,i,t)*b[jcount];
                jcount = jcount + N;
            }
        }

        icount = icount + N;
        t = t + step;
    } 
}


/* Bezier curve subroutine */

template <typename T, int N>
void dbezier(int npts, T b[], int cpts, T p[], T d1[], T d2[])
{
    int i;
    int j;
    int i1;
    int icount;
    int jcount;

    T step;
    T t;
    T tempbasis;
    T temp1;
    T temp2;

    /* calculate the points on the Bezier curve */

    icount = 0;
    t = 0;
    step = static_cast<T>(1.0)/((T)(cpts -1));

    for (i1 = 1; i1<=cpts; i1++){ /* main loop */

        if ((1.0 - t) < 5e-6) t = 1.0;

        for (j = 1; j <= N; j++){ /* generate a point on the curve */
            jcount = j;
            p[icount+j] = 0.;
            d1[icount+j] = 0.;
            d2[icount+j] = 0.;

			if (t == 0.){ /* derivatives at t=0, jcount = 1,2,3 for x,y,z */
                d1[icount + j] = (npts-1)*(b[jcount+N]- b[jcount]);
                d2[icount + j] = (npts-1)*(npts-2)*(b[jcount] - 2*b[jcount+N] + b[jcount+6]);
			}
			
            for (i = 1; i <= npts; i++){ /* Do x,y,z components */
                tempbasis = Basis(npts-1,i-1,t);
                p[icount + j] = p[icount + j] + tempbasis*b[jcount];
                if (t != 0. && t != 1.){
                    d1[icount + j] = d1[icount + j] + (((i-1)-(npts-1)*t)/(t*(1-t)))*tempbasis*b[jcount];
                    temp1 = ((i-1)-(npts-1)*t)*((i-1)-(npts-1)*t);
                    temp2 = temp1-(npts-1)*t*t - (i-1)*(1-2*t);
                    d2[icount + j] = d2[icount + j] + ((temp2)/(t*t*(1-t)*(1-t)))*tempbasis*b[jcount];
                }

   			    if (t == 1.){ /* derivatives at t=1, jcount = cpts-2,1,0 for x,y,z */
                    d1[icount + j] = (npts-1)*(b[jcount] - b[jcount-N]);
                    d2[icount + j] = (npts-1)*(npts-2)*(b[jcount] - 2*b[jcount-N] + b[jcount-6]);
                }
                jcount = jcount + N;
            }
        }

        icount = icount + N;
        t = t + step;
    } 
}

}