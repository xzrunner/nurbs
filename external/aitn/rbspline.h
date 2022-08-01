#pragma once

#include "bspline_util.h"

namespace aitn
{

/*  Name: rbais
	Language: C
	Subroutines called: none
	Book reference: Chapter 4, Sec. 4. , p 296

	c        = order of the B-spline basis function
    d        = first term of the basis function recursion relation
    e        = second term of the basis function recursion relation
	h[]	     = array containing the homogeneous weights
    npts     = number of defining polygon vertices
    nplusc   = constant -- npts + c -- maximum number of knot values
    r[]      = array containing the rationalbasis functions
               r[1] contains the basis function associated with B1 etc.
    t        = parameter value
    temp[]   = temporary array
    x[]      = knot vector
*/	

template <typename T>
void rbasis(int c, T t, int npts, int x[], T h[], T r[])
{
	int nplusc;
	int i,j,k;
	T d,e;
	T sum;
	T temp[36];

	nplusc = npts + c;

/*		printf("knot vector is \n");
		for (i = 1; i <= nplusc; i++){
			printf(" %d %d \n", i,x[i]);
		}
		printf("t is %f \n", t);
*/

/* calculate the first order nonrational basis functions n[i]	*/

	for (i = 0; i< nplusc; i++){
    	if (( t >= x[i]) && (t < x[i+1]))
			temp[i] = 1;
	    else
			temp[i] = 0;
	}

/* calculate the higher order nonrational basis functions */

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
/*
	printf("Nonrational basis functions are \n");
	for (i=1; i<= npts; i++){
		printf("%f ", temp[i]);
	}
	printf("\n");
*/
/* calculate sum for denominator of rational basis functions */

	sum = 0.;
	for (i = 0; i < npts; i++){
		    sum = sum + temp[i]*h[i];
	}

/* form rational basis functions and put in r vector */

	for (i = 0; i < npts; i++){
    	if (sum != 0){
        	r[i] = (temp[i]*h[i])/(sum);}
		else
			r[i] = 0;
	}
}

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