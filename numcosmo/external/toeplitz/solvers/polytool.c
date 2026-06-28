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
#include "polytool.h"

fdouble lau_coef(toeplitz *h,polynom *q,finteger k)
/* returns the k - th coefficient of the laurent row which is the result
   of the product of the laurent row represetned by the toeplitz matrix h
   and the polynomial q.
*/
{
  static finteger incm1=-1;
  finteger h_begin   = max(k - q->length , - h->N ) + 1;
  finteger first     = max(0,k - h->N + 1);
  finteger length    = min(q->length,h->N + k)-first;

  if (q->stdinc==-1)
    first=q->length-length;

  if (length > 0)
    return(F77CALL (ddot) (&length,&h->data[h_begin],&incm1,
	  &q->data[first],&q->stdinc));
  else
    return 0.0;
}

fdouble  poly_coef(polynom *p,finteger k)
/* returns the k - th coefficient of the polynomial p */
{
  if (k<0 || k>p->length)
    return(0.0);
  else
    return((p->stdinc>0)?(p->data[p->stdinc*k]):
	  (p->data[(p->length-k-1)*(-p->stdinc)]));
}

finteger poly_index(polynom *p,finteger k)
/* returns the array index in the array p->data of the k - th coeffitient
   of the polynomial p.
*/
{
  if (k<0 || k>=p->length)
    return -1;
  else
    return (p->stdinc>0) ? (p->stdinc*k) : ((p->length-k-1)*(-p->stdinc));
}

fdouble *data_pointer(polynom *p,finteger k1,finteger k2)
/* returns the pointer to the beginning of the data which containes
   the coefficiants form k1 to k2 of the polynomial p.
   This is for use with blas routines together with length_of_data below. */
{
  if (k1 < 0)
    k1=0;
  if (k2 >= p->length)
    k2=p->length-1;

  return &(p->data[(p->stdinc>0) ? poly_index(p,k1) : poly_index(p,k2)]);
}

finteger length_of_data(polynom *p,finteger k1,finteger *clip1,
		       finteger k2,finteger *clip2)
/* returns the the real length of the data which containes the coefficiants
   form k1 to k2 of the polynomial p. May be that 0 < k1 or length of p < k2 ..
   This is for use with blas routines together with data_pointer from above. */
{
  *clip1 = - min(0,k1);
  *clip2 = - min(p->length-1-k2,0);
  return k2 - k1 + 1 - *clip1 - *clip2;
}

void mult_by_z_up_k(polynom *res,polynom *p,fdouble *alpha,finteger k)
     /* copys the polynomial q into res and multiplies it with alpha * z^k */
{
  static finteger inc1=1;
  static finteger inc0=0;
  static fdouble zero=0.0;

  res->length=k+p->length;
  if (*alpha==0.0)
      F77CALL (dcopy) (&res->length,&zero,&inc0,res->data,&inc1);
  else
  {
    if (p->stdinc==1)
    {
      F77CALL (dcopy) (&k,&zero,&inc0,res->data,&inc1);
      F77CALL (dcopy) (&p->length,p->data,&inc1,&res->data[k],&inc1);
      F77CALL (dscal) (&p->length,alpha,&res->data[k],&inc1);
    }
    else
    {
      F77CALL (dcopy) (&p->length,p->data,&inc1,res->data,&inc1);
      F77CALL (dcopy) (&k,&zero,&inc0,&res->data[p->length],&inc1);
      F77CALL (dscal) (&p->length,alpha,res->data,&inc1);
    }
  }
}

void polynomaxpy(fdouble *alpha,polynom *x,polynom *y)
     /* performs y = alpha * x + y for polynomials similar to daxpy */
{
  if (y->stdinc==1)
    F77CALL (daxpy) (&x->length,alpha,x->data,&x->stdinc,y->data,&y->stdinc);
  else
  {
    if (x->length > y->length)
    {
      printf("polynomaxpy: can not add these polynomials.\n");
      return;
    }
    else
      F77CALL (daxpy) (&x->length,alpha,x->data,&x->stdinc,
	     &y->data[y->length - x->length],&y->stdinc);
  }
  y->length=max(y->length,x->length);
}

void polynommult(polynom *res,polynom *x,polynom *y)
/* multiplies x with y and writes the result into res */
{
  finteger i,length;
  finteger smaller_length = min(x->length,y->length);
  finteger revinc = - y->stdinc;

  res->length = x->length + y->length - 1;

  if (1 == res->stdinc)
  {
    if (1 == x->stdinc)
    {
      if (1 == y->stdinc)
      {
	for (i=0;i<res->length;i++)  
	{
	  length = min(smaller_length,i+1);
	  length = min(length,res->length-i);
	  res->data[i] = F77CALL (ddot) (&length,
               &x->data[max(i-y->length+1,0)],
	       &x->stdinc,&y->data[max(i-x->length+1,0)],&revinc);
	}
      }
      else
      {
	for (i=0;i<res->length;i++)
	{
	  length = min(smaller_length,i+1);
	  length = min(length,res->length-i);
	  res->data[i] = F77CALL (ddot) (&length,
               &x->data[max(i-y->length+1,0)],
	       &x->stdinc,&y->data[max(y->length-1-i,0)],&revinc);
	}
      }
    }
    else
    {
      if (1 == y->stdinc)
      {
	for (i=0;i<res->length;i++)  
	{
	  length = min(smaller_length,i+1);
	  length = min(length,res->length-i);
	  res->data[i] = F77CALL (ddot) (&length,
               &x->data[max(x->length-1-i,0)],
	       &x->stdinc,&y->data[max(i-x->length+1,0)],&revinc);
	}
      }
      else
      {
	for (i=0;i<res->length;i++)
	{
	  length = min(smaller_length,i+1);
	  length = min(length,res->length-i);
	  res->data[i] = F77CALL (ddot) (&length,
               &x->data[max(x->length-1-i,0)],
	       &x->stdinc,&y->data[max(y->length-1-i,0)],&revinc);
	}
      }
    }
  }
  else
  {
    if (1 == x->stdinc)
    {
      if (1 == y->stdinc)
      {
	for (i=0;i<res->length;i++)
	{
	  length = min(smaller_length,i+1);
	  length = min(length,res->length-i);
	  res->data[i] = F77CALL (ddot) (&length,
               &x->data[max(x->length-1-i,0)],
	       &x->stdinc,&y->data[max(y->length-1-i,0)],&revinc);
	}
      }
      else
      {
	for (i=0;i<res->length;i++)
	{
	  length = min(smaller_length,i+1);
	  length = min(length,res->length-i);
	  res->data[i] = F77CALL (ddot) (&length,
               &x->data[max(x->length-1-i,0)],
	       &x->stdinc,&y->data[max(1-x->length+i,0)],&revinc);
	}
      }
    }
    else
    {
      if (1 == y->stdinc)
      {
	for (i=0;i<res->length;i++)
	{
	  length = min(smaller_length,i+1);
	  length = min(length,res->length-i);
	  res->data[i] = F77CALL (ddot) (&length,
               &x->data[max(1-y->length+i,0)],
	       &x->stdinc,&y->data[max(y->length-1-i,0)],&revinc);
	}
      }
      else
      {
	for (i=0;i<res->length;i++)
	{
	  length = min(smaller_length,i+1);
	  length = min(length,res->length-i);
	  res->data[i] = F77CALL (ddot) (&length,
               &x->data[max(1-y->length+i,0)],
	       &x->stdinc,&y->data[max(1-x->length+i,0)],&revinc);
	}
      }
    }
  }
}
