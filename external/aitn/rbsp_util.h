#pragma once

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

}