/*************************************************************************/
/* Copyright (C) 1996-1997 Marlis Hochbruck                              */
/* All rights reserved.                                                  */
/*                                                                       */
/* C code written by Mathias Froehlich.                                  */
/*                                                                       */
/* This code is part of a copyrighted package. For details, see the file */
/* "COPYING" in the top-level directory.                                 */
/*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftypes.h"
#include "fblas.h"
#include "toeplitz.h"

int d_lev_cl(fdouble *x,fdouble *mu,finteger N,fdouble *b)
     /* classical levinson algorithem
	x     pointer to the array for the solution.
        mu    pointer to the diagonal element of the toeplitz matrix.
        N     the dimension of the toeplitz matrix mu.
        b     pointer to the righthandside. */
{
  extern fdouble lau_coef(toeplitz *h,polynom *q,finteger k);

  /* finteger help variables */
  finteger n=1;
  finteger ihelp;
  finteger success=1;

  /* static increments for blas */
  static finteger inc1=1;
  static finteger inc0=0;
  static finteger incm1=-1;

  /* tolerance for singular matrix */
  const fdouble col_tol=0.0;

  /* help variables */
  fdouble help;
  fdouble *temp;

  /* variables for the method */
  polynom q_up,q;
  toeplitz h;
  fdouble p_up,p;
  fdouble v_up,v;
  fdouble e;

  /* initialize q,q_up,h */
  h.data=mu;
  h.N=N;

  q.length=0;
  q.stdinc=1;

  q_up.length=0;
  q_up.stdinc=-1;

  /* memory management */
  q.data=(fdouble *)malloc(3*N*sizeof(fdouble));
  if (q.data==NULL)
  {
    printf("error can't allocate memory for working!\nAbort ...\n");
    return(0);
  }
  q_up.data=&q.data[N];
  temp=&q_up.data[N];

  /* initialize memory */
  help=0.0;
  ihelp=2*N;
  F77CALL (dcopy) (&ihelp,&help,&inc0,q.data,&inc1);
  F77CALL (dcopy) (&N,&help,&inc0,x,&inc1);


  /* LDU decompostion of a scalar */
  p=h.data[0];
  q.data[0]=1.0;
  q.length++;
  q_up.data[0]=1.0;
  q_up.length++;

  /* update the solution */
  x[0]=b[0]/p;

  /* calculating the pi ... */
  e=h.data[1];
  p_up=h.data[-1];

  while (n<N)
  {
    /* calculating the gammas */
    v=-e/p;
    v_up=-p_up/p;

    /* calculating new q and new q_up */
    /* saving q_up into temp */
    F77CALL (dcopy) (&q_up.length,q_up.data,&q_up.stdinc,temp,&inc1);

    /* q_up = [0 q_up] + v_up * [q 0] */
    F77CALL (daxpy) (&q_up.length,&v_up,q.data,&q.stdinc,&q_up.data[1],
		     &q_up.stdinc);
    q_up.length++;

    /* q = [0 q_up] * v + [q 0] */
    F77CALL (daxpy) (&q.length,&v,temp,&inc1,&q.data[1],&q.stdinc);
    q.length++;

    /* updating epsilon */
    p*=(1.0-v_up*v);

    /* if matrix is singular abort */
    if (fabs(p)<=col_tol)
    {
      printf("Matrix is singular\n");
      success=0;
      break;
    }

    n++;


    /* updating the solution x (Here b with incm1 because 
       q=[rho_n rho_n-1 ... rho_0] but the columns of the matrix U in the
       LDU decomposition of the inverse of the toeplitz matrix are:
       [rho_0 rho_1 ... rho_n]. So you have to increment q in the other
       direction as the stdinc says. That is the same as incrementing
       q with stdinc and b with incm1.). */
    help=F77CALL (ddot) (&q.length,q.data,&q.stdinc,b,&incm1)/p;
    F77CALL (daxpy) (&q_up.length,&help,q_up.data,&q_up.stdinc,x,&inc1);


    /* calculating the pi's ... */
    p_up=lau_coef (&h,&q_up,-1);

    e=lau_coef (&h,&q,n);
  }

  free (q.data);
  return (success);
}












