/*************************************************************************/
/* Copyright (C) 1996-1997 Marlis Hochbruck                              */
/* All rights reserved.                                                  */
/*                                                                       */
/* C code written by Mathias Froehlich.                                  */
/*                                                                       */
/* This code is part of a copyrighted package. For details, see the file */
/* "COPYING" in the top-level directory.                                 */
/*************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "ftypes.h"
#include "fblas.h"
#include "toetool.h"

void d_toe_mv(finteger N,finteger first,finteger last,fdouble *alpha,
	      fdouble *a,fdouble *x,fdouble *beta,fdouble *y)
     /* preforms a matrix vector multiplication:
	y = beta * y + alpha * A * x
        where A is represented by the first line and first column
        stored in the array a */
{
  finteger i;
  finteger length;
  finteger index_t,index_x;
  static finteger inc1=1;
  static finteger incm1=-1;

  if (*beta==0.0)
  {
    for(i=0;i<N;i++)
    {
      length=N;
      index_t=max(i-N+1,first);
      length+=i-N+1-index_t;
      index_x=index_t+length-1-min(index_t+length-1,last);
      length-=index_x;
      length=max(length,0);

      y[i]=(*alpha)*F77CALL (ddot) (&length,&a[index_t],&inc1,
	    &x[index_x],&incm1);
    }
  }
  else
  {
    if (*beta!=1.0)
      F77CALL (dscal) (&N,beta,y,&inc1);

    for(i=0;i<N;i++)
    {
      length=N;
      index_t=max(i-N+1,first);
      length+=i-N+1-index_t;
      index_x=index_t+length-1-min(index_t+length-1,last);
      length-=index_x;
      length=max(length,0);

      y[i]+=(*alpha)*F77CALL (ddot) (&length,&a[index_t],&inc1,
	    &x[index_x],&incm1);
    }
  }
}

void d_toe_mm(finteger N,finteger first,finteger last,finteger m,fdouble *alpha,
	      fdouble *a,fdouble *b__,finteger ldb,fdouble *beta,
	      fdouble *c__,finteger ldc)
     /* preforms a matrix matrix multiplication:
	C = beta * C + alpha * A * B
        where A is toeplitz and represented by the first line and first column
        stored in the array a. B,C are full matrices  */
{
  finteger i;
#define b(I,J) b__[(I) + (J)*ldb]
#define c(I,J) c__[(I) + (J)*ldc]

  for (i=0;i<m;i++)
    d_toe_mv(N,first,last,alpha,a,&b(0,i),beta,&c(0,i));

#undef c
#undef b
}

void d_toe2full(finteger N,fdouble *to,fdouble *fu,finteger *ld)
     /* converts a toeplitz matrix a into a full matrix fu */
{
  finteger i;
  static finteger inc1=1;

#define Fu(I,J) fu[(I) + (J) * (*ld)]

  for(i=0;i<N;i++)
    F77CALL (dcopy) (&N,&to[-i],&inc1,&Fu(0,i),&inc1);

#undef Fu
}
