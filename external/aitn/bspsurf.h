/*  Subroutine to generate a Bezier curve.
    Copyright (c) 2000 David F. Rogers. All rights reserved.
*/

#pragma once

#include "bsp_util.h"

namespace aitn
{

/*  Subroutine to calculate a Cartesian product B-spline surface
    using open uniform knot vectors (see Eq. 6.1}).

	Name: bsplsurf.c
	Language: C
	Subroutines called: knot.c, basis.c
	Book reference: Chapter 6, Ex. 6.1, Alg. p. 303

    b[]         = array containing the polygon net points
                  b[1] = x-component
                  b[2] = y-component
                  b[3] = z-component
                  Note: Bi,j = b[] has dimensions of n*m*3 with j varying fastest
                      The polygon net is n x m
    k           = order in the u direction
    l           = order in the w direction
    mbasis[]    = array containing the nonrational basis functions for one value of w (see \eq{5--84})
    mpts        = the number of polygon vertices in w direction
    nbasis[]    = array containing the nonrational basis functions for one value of u (see \eq{5--84})
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
void bsplsurf(T b[], int k, int l, int npts, int mpts, int p1, int p2, T q[])
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
    T u,w;
    T stepu,stepw;
    
/*    printf("in bsplsurf.c \n");*/

/*    zero and redimension the arrays */


/*    printf("k,l,npts,mpts,p1,p2 = %d %d %d %d %d %d \n",k,l,npts,mpts,p1,p2);*/


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

/*  printf("p1*p2*3 = %d \n",temp);*/

    for (i = 0; i < 3*p1*p2; i++){
        q[i] = 0.;
    }

/*   generate the open uniform knot vectors */

    knot(npts,k,x);       /*  calculate u knot vector */
    knot(mpts,l,y);       /*  calculate w knot vector */

    icount = 0;

/*    calculate the points on the \bsp surface */

    stepu = (T)x[nplusc - 1]/(T)(p1-1);
    stepw = (T)y[mplusc - 1]/(T)(p2-1);
    u = 0.;
    for (uinc = 0; uinc < p1; uinc++)
    {
/*      printf("u = %f \n",u); */
        if ((T)x[nplusc - 1] - u < 5e-6){
            u = (T)x[nplusc - 1];
        }
        basis(k,u,npts,x,nbasis);    /* basis function for this value of u */
        w = 0.;
        for (winc = 0; winc < p2; winc++)
        {
/*          printf("w = %f \n",w);*/
            if ((T)y[mplusc - 1] - w < 5e-6){
                w = (T)y[mplusc - 1];
            }
            basis(l,w,mpts,y,mbasis);    /* basis function for this value of w */
            for (i = 0; i < npts; i++)
            {
                if (nbasis[i] != 0. )
                {
                    jbas = 3 * (mpts + 1) * i;
                    for (j = 0; j < mpts; j++)
                    {
                        if (mbasis[j] != 0.)
                        {
                            j1 = jbas +3 * j;
                            pbasis = nbasis[i]*mbasis[j];
                            q[icount] = q[icount]+b[j1]*pbasis;  /* calculate surface point */
                            q[icount+1] = q[icount+1]+b[j1+1]*pbasis;
                            q[icount+2] = q[icount+2]+b[j1+2]*pbasis;
    /*                      printf("j1,i,j = %d %d %d \n",j1,i,j); */
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

/*  Subroutine to calculate a Cartesian product B-spline surface
    using uniform periodic knot vectors (see Eq. 6.1).

	Name: bspsurfu.c
	Language: C
	Subroutines called: knotu.c, basis.c
	Book reference: p. 304
	
    b[]         = array containing the polygon net points
                  b[1] = x-component
                  b[2] = y-component
                  b[3] = z-component
                  Note: Bi,j = b[] has dimensions of n*m*3 with j varying fastest
                      The polygon net is n x m
    k           = order in the u direction
    l           = order in the w direction
    mbasis[]    = array containing the nonrational basis functions for one value of w (see \eq{5--84})
    mpts        = the number of polygon vertices in w direction
    nbasis[]    = array containing the nonrational basis functions for one value of u (see \eq{5--84})
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
void bspsurfu(T b[], int k, int l, int npts, int mpts, int p1, int p2, T q[])
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
    T u,w;
    T stepu,stepw;
    
/*    printf("in bsplsurf.c \n");*/

/*    zero and redimension the arrays */


/*    printf("k,l,npts,mpts,p1,p2 = %d %d %d %d %d %d \n",k,l,npts,mpts,p1,p2);*/


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

/*  printf("p1*p2*3 = %d \n",temp);*/

    for (i = 0; i < 3*p1*p2; i++){
        q[i] = 0.;
    }

/*   generate the open uniform knot vectors */

    knotu(npts,k,x);       /*  calculate u knot vector */
    knotu(mpts,l,y);       /*  calculate w knot vector */

    icount = 0;

/*    calculate the points on the B-spline surface */

    stepu = (T)((npts-1+1)-(k-1))/(T)(p1-1);
    stepw = (T)((mpts-1+1)-(l-1))/(T)(p2-1);
    u = k-1.;
    for (uinc = 0; uinc < p1; uinc++)
    {
/*      printf("u = %f \n",u); */
        if ((T)(npts) - u < 5e-6){
            u = npts;
        }
        basis(k,u,npts,x,nbasis);    /* basis function for this value of u */
        w = l-1.;
        for (winc = 0; winc < p2; winc++)
        {
/*          printf("w = %f \n",w);*/
            if ((T)mpts - w < 5e-6){
                w = mpts;
            }
            basis(l,w,mpts,y,mbasis);    /* basis function for this value of w */
            for (i = 0; i < npts; i++)
            {
                if (nbasis[i] != 0. )
                {
                    jbas = 3 * (mpts + 1) * i;
                    for (j = 0; j < mpts; j++)
                    {
                        if (mbasis[j] != 0.)
                        {
                            j1 = jbas + 3 * j;
                            pbasis = nbasis[i]*mbasis[j];
                            q[icount] = q[icount]+b[j1]*pbasis;  /* calculate surface point */
                            q[icount+1] = q[icount+1]+b[j1+1]*pbasis;
                            q[icount+2] = q[icount+2]+b[j1+2]*pbasis;
    /*                      printf("j1,i,j = %d %d %d \n",j1,i,j); */
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