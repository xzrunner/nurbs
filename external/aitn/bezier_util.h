/*  Subroutine to generate a Bezier curve.
    Copyright (c) 2000 David F. Rogers. All rights reserved.
*/

#pragma once

namespace aitn
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

}