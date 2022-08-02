#pragma once

#include "bsp_util.h"
#include "rbsp_util.h"

namespace aitn
{

/*  Name: rbspline.c
	Language: C
	Subroutines called: knot.c, rbasis.c, fmtmul.c
	Book reference: Chapter 4, Alg. p. 297

    b[]         = array containing the defining polygon vertices
                  b[1] contains the x-component of the vertex
                  b[2] contains the y-component of the vertex
                  b[3] contains the z-component of the vertex
	h[]			= array containing the homogeneous weighting factors 
    k           = order of the B-spline basis function
    nbasis      = array containing the basis functions for a single value of t
    nplusc      = number of knot values
    npts        = number of defining polygon vertices
    p[,]        = array containing the curve points
                  p[1] contains the x-component of the point
                  p[2] contains the y-component of the point
                  p[3] contains the z-component of the point
    p1          = number of points to be calculated on the curve
    t           = parameter value 0 <= t <= npts - k + 1
    x[]         = array containing the knot vector
*/

template <typename T>
void rbspline(int npts, int k, int p1, T b[], T h[], T p[])
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

/*    calculate the points on the rational B-spline curve */

	t = 0;
	step = ((T)x[nplusc - 1])/((T)(p1-1));

	for (i1 = 0; i1< p1; i1++){

		if ((T)x[nplusc - 1] - t < 5e-6){
			t = (T)x[nplusc - 1];
		}

	    rbasis(k,t,npts,x,h,nbasis);      /* generate the basis function for this value of t */
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

/*  Name: rbsplinu.c
	Language: C
	Subroutines called: knotu.c, rbasis.c, fmtmul.c
	Book reference: Chapter 4, Alg. p. 298

    b[]         = array containing the defining polygon vertices
                  b[1] contains the x-component of the vertex
                  b[2] contains the y-component of the vertex
                  b[3] contains the z-component of the vertex
	h[]			= array containing the homogeneous weighting factors 
    k           = order of the B-spline basis function
    nbasis      = array containing the basis functions for a single value of t
    nplusc      = number of knot values
    npts        = number of defining polygon vertices
    p[,]        = array containing the curve points
                  p[1] contains the x-component of the point
                  p[2] contains the y-component of the point
                  p[3] contains the z-component of the point
    p1          = number of points to be calculated on the curve
    t           = parameter value 0 <= t <= npts - k + 1
    x[]         = array containing the knot vector
*/

template <typename T>
void rbsplineu(int npts, int k, int p1, T b[], T h[], T p[])
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

/* generate the uniform periodic knot vector */

	knotu(npts,k,x);

/*
	printf("The knot vector is ");
	for (i = 1; i <= nplusc; i++){
		printf(" %d ", x[i]);
	}
	printf("\n");

	printf("The usable parameter range is ");
	for (i = k; i <= npts+1; i++){
		printf(" %d ", x[i]);
	}
	printf("\n");
*/

	icount = 0;

/*    calculate the points on the rational B-spline curve */

	t = k-1;
	step = ((T)((npts)-(k-1)))/((T)(p1-1));	

	for (i1 = 0; i1< p1; i1++){

		if ((T)x[nplusc - 1] - t < 5e-6){
			t = (T)x[nplusc - 1];
		}

	    rbasis(k,t,npts,x,h,nbasis);      /* generate the basis function for this value of t */
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