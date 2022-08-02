/*  Subroutine to generate a Bezier curve.
    Copyright (c) 2000 David F. Rogers. All rights reserved.
*/

#pragma once

#include "bsp_util.h"
#include "rbsp_util.h"

namespace aitn
{

/*  Subroutine to calculate the sum of the nonrational basis functions (see \eq {6--87}).
	
	Name: sumrbas.c
	Language: C
	Subroutines called: none
	Book reference: Chapter 7, Section 7.1, Alg. p.309
	
    b[]         = array containing the polygon net points
                    b[1] = x-component
                    b[2] = y-component
                    b[3] = z-component
                    b[4] = homogeneous coordinate weighting factor        
                    Note: Bi,j = b[] has dimensions of n*m x 4 with j varying fastest
                     The polygon net is n x m
    mbasis[,]  = array containing the nonrational basis functions for one value of w
    mpts       = the number of polygon vertices in w direction
    nbasis[,]  = array containing the nonrational basis functions for one value of u
    npts       = the number of polygon vertices in u direction
    sum        = sum of the basis functions
*/

template <typename T>
T sumrbas(T b[], T nbasis[], T mbasis[], int npts, int mpts)
{

	int i,j,jbas,j1;
	T sum;
	
 /* calculate the sum */

	sum = 0;
/*	printf("npts,mpts %d %d \n", npts,mpts);*/

	for (i = 0; i < npts; i++){
		if (nbasis[i] != 0.){
			jbas = 4*mpts*i;
		    for (j = 0; j < mpts; j++){
/*				printf("i,j,jbas %d %d %d \n",i,j,jbas);*/
		    	if (mbasis[j] != 0.){
		    	    j1 = jbas + 4*j + 3;
/*					printf("in sumrbas j1 = %d \n",j1);*/
	    	    	sum = sum + b[j1]*nbasis[i]*mbasis[j];
				}
			}
		}
	}
	return(sum);
}

/*  Subroutine to calculate a Cartesian product rational B-spline
    surface using open uniform knot vectors (see Eq. (7.1)).

	Name: rbspsurf.c
	Language: C
	Subroutines called: knot.c, basis.c sumrbas.c
	Book reference: Chapter 7, Section 7.1, Alg. p. 308

    b[]         = array containing the polygon net points
                  b[1] = x-component
                  b[2] = y-component
                  b[3] = z-component
                  b[4] = h-component
                  Note: Bi,j = b[] has dimensions of n*m*3 with j varying fastest
                      The polygon net is n x m
    k           = order in the u direction
    l           = order in the w direction
    mbasis[]    = array containing the nonrational basis functions for one value of w (see Eq. (3.2))
    mpts        = the number of polygon vertices in w direction
    nbasis[]    = array containing the nonrational basis functions for one value of u (see Eq. (3.2))
    npts        = the number of polygon vertices in u direction
    p1          = number of parametric lines in the u direction
    p2          = number of parametric lines in the w direction
    q[]         = array containing the resulting surface
                  q[1] = x-component
                  q[2] = y-component
                  q[3] = z-component
                       for a fixed value of u the next m elements contain
                       the values for the curve q[u[sub i],w] q has dimensions
                       of p1*p2*3. The display surface is p1 x p2
*/

template <typename T>
void rbspsurf(T b[], int k, int l, int npts, int mpts, int p1, int p2, T q[])
{

/*   allows for 20 data points with basis function of order 5 */

    int i,j,j1,jbas;
    int icount;
    int uinc,winc;
    int nplusc,mplusc;
    int x[30],y[30];
    int temp;

    T nbasis[30],mbasis[30];
    T pbasis;
    T h;
    T sum;
    T u,w;
    T stepu,stepw;
        
/*  printf("in bsplsurf.c \n"); */

/*    zero and redimension the arrays */

    nplusc = npts + k;
    mplusc = mpts + l;    

    for (i = 0; i < nplusc; i++){
        x[i] = 0;
    }
    for (i = 0; i < mplusc; i++){
        y[i] = 0;
    }
    for (i = 0; i < npts; i++){
        nbasis[i] = 0.;
    }
    for (i = 0; i < mpts; i++){
        mbasis[i] = 0.;
    }

    temp = 3*(p1*p2);

    for (i = 0; i < 3*p1*p2; i++){
        q[i] = 0.;
    }

/*   generate the open uniform knot vectors */

    knot(npts,k,x);       /*  calculate u knot vector */
    knot(mpts,l,y);       /*  calculate w knot vector */

    icount = 0;

/*    calculate the points on the B-spline surface */

    stepu = (T)x[nplusc - 1]/(T)(p1-1);
    stepw = (T)y[mplusc - 1]/(T)(p2-1);
    u = 0.;
    for (uinc = 0; uinc < p1; uinc++)
    {
        if ((T)x[nplusc - 1] - u < 5e-6){
            u = (T)x[nplusc - 1];
        }
        basis(k,u,npts,x,nbasis);    /* basis function for this value of u */
        w = 0.;
        for (winc = 0; winc < p2; winc++)
        {
            if ((T)y[mplusc - 1] - w < 5e-6){
                w = (T)y[mplusc - 1];
            }
            basis(l,w,mpts,y,mbasis);    /* basis function for this value of w */
            sum = sumrbas(b,nbasis,mbasis,npts,mpts);
/*        	printf("in rbspsurf u,w,sum = %f %f %f \n",u,w,sum);*/
            for (i = 0; i < npts; i++)
            {
                if (nbasis[i] != 0.)
                {
                    jbas = 4*(mpts + 1)*i;
                    for (j = 0; j < mpts; j++)
                    {
                        if (mbasis[j] != 0.)
                        {
                            j1 = jbas +4*j;
                            pbasis = b[j1+3]*nbasis[i]*mbasis[j]/sum;
                            q[icount] = q[icount]+b[j1]*pbasis;  /* calculate surface point */
                            q[icount+1] = q[icount+1]+b[j1+1]*pbasis;
                            q[icount+2] = q[icount+2]+b[j1+2]*pbasis;
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