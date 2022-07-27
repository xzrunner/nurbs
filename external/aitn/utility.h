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

}