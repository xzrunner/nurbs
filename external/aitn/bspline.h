/*  Subroutine to generate a B-spline curve.
    Copyright (c) 2000 David F. Rogers. All rights reserved.
*/

#pragma once

#include 	<stdio.h>

namespace aitn
{

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

/*  Subroutine to generate a B-spline curve using an uniform open knot vector

	Name: bspline.c
	Language: C
	Subroutines called: knot.c, basis.c, fmtmul.c
	Book reference: Section 3.5, Ex. 3.4, Alg. p. 281

    b[]        = array containing the defining polygon vertices
                  b[1] contains the x-component of the vertex
                  b[2] contains the y-component of the vertex
                  b[3] contains the z-component of the vertex
    k           = order of the \bsp basis function
    nbasis      = array containing the basis functions for a single value of t
    nplusc      = number of knot values
    npts        = number of defining polygon vertices
    p[,]        = array containing the curve points
                  p[1] contains the x-component of the point
                  p[2] contains the y-component of the point
                  p[3] contains the z-component of the point
    p1          = number of points to be calculated on the curve
    t           = parameter value 0 <= t <= 1
    x[]         = array containing the knot vector
*/

template <typename T>
void bspline(int npts, int k, int p1, T b[], T p[])
{
	int i,j,icount,jcount;
	int i1;
	int x[30];		/* allows for 20 data points with basis function of order 5 */
	int nplusc;

	T step;
	T t;
	T nbasis[20];
	T temp;

	nplusc = npts + k;

/*  zero and redimension the knot vector and the basis array */

	for(i = 0; i < npts; i++){
		 nbasis[i] = 0.;
	}

	for(i = 0; i < nplusc; i++){
		 x[i] = 0.;
		}

/* generate the uniform open knot vector */

	knot(npts,k,x);

/*
	printf("The knot vector is ");
	for (i = 1; i <= nplusc; i++){
		printf(" %d ", x[i]);
	}
	printf("\n");
*/

	icount = 0;

/*    calculate the points on the bspline curve */

	t = 0;
	step = ((T)x[nplusc - 1])/((T)(p1-1));

	for (i1 = 0; i1< p1; i1++){

		if ((T)x[nplusc - 1] - t < 5e-6){
			t = (T)x[nplusc - 1];
		}

	    basis(k,t,npts,x,nbasis);      /* generate the basis function for this value of t */
/*
		printf("t = %f \n",t);
		printf("nbasis = ");
		for (i = 1; i <= npts; i++){
			printf("%f  ",nbasis[i]);
		}
		printf("\n");
*/
		for (j = 0; j < 3; j++){      /* generate a point on the curve */
			jcount = j;
			p[icount+j] = 0.;

			for (i = 0; i < npts; i++){ /* Do local matrix multiplication */
				temp = nbasis[i]*b[jcount];
			    p[icount + j] = p[icount + j] + temp;
/*
				printf("jcount,nbasis,b,nbasis*b,p = %d %f %f %f %f\n",jcount,nbasis[i],b[jcount],temp,p[icount+j]);
*/
				jcount = jcount + 3;
			}
		}
/*
		printf("icount, p %d %f %f %f \n",icount,p[icount+1],p[icount+2],p[icount+3]);
*/
    	icount = icount + 3;
		t = t + step;
	}
}

/*  Subroutine to generate a B-spline curve using an uniform open knot vector

	C code for An Introduction to NURBS
	by David F. Rogers. Copyright (C) 2000 David F. Rogers,
	All rights reserved.
	
	Name: bsplineu.c
	Language: C
	Subroutines called: knotu.c, basis.c, fmtmul.c
	Book reference: Section 3.8, Ex. 3.7, Alg. p. 282

    b[]        = array containing the defining polygon vertices
                  b[1] contains the x-component of the vertex
                  b[2] contains the y-component of the vertex
                  b[3] contains the z-component of the vertex
    k           = order of the B-spline basis function
    nbasis      = array containing the basis functions for a single value of t
    nplusc      = number of knot values
    npts        = number of defining polygon vertices
    p[,]        = array containing the curve points
                  p[1] contains the x-component of the point
                  p[2] contains the y-component of the point
                  p[3] contains the z-component of the point
    p1          = number of points to be calculated on the curve
    t           = parameter value 0 <= t <= 1
    x[]         = array containing the knot vector
*/

template <typename T>
void bsplineu(int npts, int k, int p1, T b[], T p[])
{
	int i,j,icount,jcount;
	int i1;
	int x[30];		/* allows for 20 data points with basis function of order 5 */
	int nplusc;

	T step;
	T t;
	T nbasis[20];
	T temp;

	nplusc = npts + k;

/*  zero and redimension the knot vector and the basis array */

	for(i = 0; i < npts; i++){
		 nbasis[i] = 0.;
	}

	for(i = 0; i < nplusc; i++){
		 x[i] = 0.;
		}

/* generate the uniform open knot vector */

	knotu(npts,k,x);

/*
	printf("The knot vector is ");
	for (i = 1; i <= nplusc; i++){
		printf(" %d ", x[i]);
	}
	printf("\n");
*/

	icount = 0;

/*    calculate the points on the bspline curve */

	t = k - 1; /* special parameter range for periodic basis functions */
	step = ((T)((npts)-(k - 1))) / ((T)(p1 - 1));

	for (i1 = 0; i1< p1; i1++){

		if ((T)(npts)-t < 5e-6) {
			t = (T)((npts));
		}

	    basis(k,t,npts,x,nbasis);      /* generate the basis function for this value of t */
/*
		printf("t = %f \n",t);
		printf("nbasis = ");
		for (i = 1; i <= npts; i++){
			printf("%f  ",nbasis[i]);
		}
		printf("\n");
*/
		for (j = 0; j < 3; j++){      /* generate a point on the curve */
			jcount = j;
			p[icount+j] = 0.;

			for (i = 0; i < npts; i++){ /* Do local matrix multiplication */
				temp = nbasis[i]*b[jcount];
			    p[icount + j] = p[icount + j] + temp;
/*
				printf("jcount,nbasis,b,nbasis*b,p = %d %f %f %f %f\n",jcount,nbasis[i],b[jcount],temp,p[icount+j]);
*/
				jcount = jcount + 3;
			}
		}
/*
		printf("icount, p %d %f %f %f \n",icount,p[icount+1],p[icount+2],p[icount+3]);
*/
    	icount = icount + 3;
		t = t + step;
	}
}

}