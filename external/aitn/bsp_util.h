/*  Subroutine to generate a B-spline curve.
    Copyright (c) 2000 David F. Rogers. All rights reserved.
*/

#pragma once

namespace aitn
{

/*
    Subroutine to generate a B-spline open knot vector with multiplicity
    equal to the order at the ends.
	
    c            = order of the basis function
    n            = the number of defining polygon vertices
    nplus2       = index of x() for the first occurence of the maximum knot vector value
    nplusc       = maximum value of the knot vector -- $n + c$
    x()          = array containing the knot vector
*/

void knot(int n, int c,  int x[])
{
	int nplusc, nplus2, i;

	nplusc = n + c;
	nplus2 = n + 1;

	x[0] = 0;
	for (i = 1; i < nplusc; i++) {
		if ((i >= c) && (i < nplus2))
			x[i] = x[i - 1] + 1;
		else
			x[i] = x[i - 1];
	}
}

/*  Subroutine to generate a B-spline uniform (periodic) knot vector.

    c            = order of the basis function
    n            = the number of defining polygon vertices
    nplus2       = index of x() for the first occurence of the maximum knot vector value
    nplusc       = maximum value of the knot vector -- $n + c$
    x[]          = array containing the knot vector
*/

void knotu(int n, int c, int x[])
{
    int nplusc, i;

	nplusc = n + c;

	x[0] = 0;
	for (i = 1; i < nplusc; i++){
	    x[i] = i;
	}
}

/*  Subroutine to generate B-spline basis functions for open knot vectors
	
	Name: basis.c
	Language: C
	Subroutines called: none
	Book reference: p. 279

    c        = order of the B-spline basis function
    d        = first term of the basis function recursion relation
    e        = second term of the basis function recursion relation
    npts     = number of defining polygon vertices
    n[]      = array containing the basis functions
               n[1] contains the basis function associated with B1 etc.
    nplusc   = constant -- npts + c -- maximum number of knot values
    t        = parameter value
    temp[]   = temporary array
    x[]      = knot vector
*/	

template <typename T>
void basis(int c, T t, int npts, int x[], T n[])
{
	int nplusc;
	int i,k;
	T d,e;
	T temp[36];

	nplusc = npts + c;

/*		printf("knot vector is \n");
		for (i = 1; i <= nplusc; i++){
			printf(" %d %d \n", i,x[i]);
		}
		printf("t is %f \n", t);
*/

/* calculate the first order basis functions n[i][1]	*/

	for (i = 0; i < nplusc; i++){
    	if (( t >= x[i]) && (t < x[i+1]))
			temp[i] = 1;
	    else
			temp[i] = 0;
	}

/* calculate the higher order basis functions */

	for (k = 2; k <= c; k++){
    	for (i = 0; i < nplusc-k; i++){
        	if (temp[i] != 0)    /* if the lower order basis function is zero skip the calculation */
           		d = ((t-x[i])*temp[i])/(x[i+k-1]-x[i]);
	        else
				d = 0;

    	    if (temp[i+1] != 0)     /* if the lower order basis function is zero skip the calculation */
        		e = ((x[i+k]-t)*temp[i+1])/(x[i+k]-x[i+1]);
	        else
    			e = 0;

    	    temp[i] = d + e;
		}
	}

	if (t == (T)x[nplusc - 1]){		/*    pick up last point	*/
 		temp[npts - 1] = 1;
	}

/* put in n array	*/

	for (i = 0; i < npts; i++) {
    	n[i] = temp[i];
	}
}

}