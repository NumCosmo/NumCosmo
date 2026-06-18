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
#include "drefine.h"

extern void d_toe_mv(finteger N,finteger first,finteger last,fdouble *alpha,
       fdouble *a,fdouble *x,fdouble *beta,fdouble *y);

static void invert_reg(fdouble *x,toeplitz *h,fdouble *b,
       polynom *q,polynom *q_up,fdouble *p);

static void invert_antireg(fdouble *x,toeplitz *h,fdouble *b,
       polynom *q_lef,polynom *q_up);

static finteger inc1  =  1;
static finteger incm1 = -1;

static fdouble one  =  1.0;
static fdouble zero =  0.0;
static fdouble mone = -1.0;

/* variables for refinement */
static int        refine_pair_antireg;
static fdouble *refine_work = NULL;
static finteger    refine_work_length;
static polynom    q1,q2;
static fdouble p;
static toeplitz   h;

int malloc_for_refinement(finteger dimension)
{
  dimension++;

  refine_work_length = 4*dimension-1;
  if (refine_work == NULL)
  {
    refine_work = (fdouble *)malloc(sizeof(fdouble)*refine_work_length);
    if (refine_work == NULL)
    {
      fprintf(stderr,"set_refine_flag Error:\nCan't allocate memory for "
	      "working.\nAbort ...\n");
      return 0;
    }
  }
  else
  {
    refine_work = (fdouble *)realloc((void*)refine_work,
					sizeof(fdouble)*refine_work_length);
    if (refine_work == NULL)
    {
      fprintf(stderr,"set_refine_flag Error:\nCan't allocate memory for "
	      "working.\nAbort ...\n");
      return 0;
    }
  }
  q1.data=refine_work;
  q2.data=&refine_work[dimension];
  h.data=&refine_work[3*dimension-1];
  return 1;
}

void free_refinement_workspace()
{
  free(refine_work);
  refine_work = NULL;
}

int set_refinement(toeplitz *h_in,polynom *q,polynom *q_up,polynom *q_lef,
		   fdouble *p_in)
{
  finteger ihelp;

  if (malloc_for_refinement(h_in->N))
  {
    if (p_in != NULL)
      p = *p_in;

    ihelp = 2*h_in->N - 1;
    F77CALL (dcopy) (&ihelp,&h_in->data[1-h_in->N],&inc1,&h.data[1-h_in->N],&inc1);
    h.N=h_in->N;

    if (q != NULL)
    {
      refine_pair_antireg = 0;

      F77CALL (dcopy) (&q->length,q->data,&inc1,q1.data,&inc1);
      q1.stdinc=q->stdinc;
      q1.length=q->length;
    }

    if (q_lef != NULL)
    {
      refine_pair_antireg = 1;

      F77CALL (dcopy) (&q_lef->length,q_lef->data,&inc1,q1.data,&inc1);
      q1.stdinc=q_lef->stdinc;
      q1.length=q_lef->length;
    }

    F77CALL (dcopy) (&q_up->length,q_up->data,&inc1,q2.data,&inc1);
    q2.stdinc=q_up->stdinc;
    q2.length=q_up->length;

    return 1;
  }
  else
    return 0;
}

void refine_solution(fdouble *x,fdouble *b)
{
  fdouble *temp1,*temp2;

  temp1 = (fdouble *)malloc(sizeof(fdouble)*2*h.N);
  if (temp1 == NULL)
  {
    fprintf(stderr,"can't allocate memory for refining.\n"
	    "Abort refining ...\n");
    return;
  }
  temp2=&temp1[h.N];

  d_toe_mv(h.N,-h.N+1,h.N-1,&one,h.data,x,&zero,temp1);
  F77CALL (daxpy) (&h.N,&mone,b,&inc1,temp1,&inc1);

  invert_toeplitz(temp2,temp1);

  F77CALL (daxpy) (&h.N,&mone,temp2,&inc1,x,&inc1);

  free(temp1);
}

void invert_toeplitz(fdouble *x,fdouble *b)
{
  if (refine_pair_antireg)
    invert_antireg(x,&h,b,&q1,&q2);
  else
    invert_reg(x,&h,b,&q1,&q2,&p);
}

static void invert_reg(fdouble *x,toeplitz *h,fdouble *b,
		polynom *q,polynom *q_up,fdouble *p)
{
  fdouble *temp;
  finteger i,length;

  temp=(fdouble *)malloc(sizeof(fdouble)*h->N);
  if (temp==NULL)
  {
    fprintf(stderr,"can't allocate memory for refining.\n"
	    "Abort refining ...\n");
    return;
  }

  for (i=0,length=h->N;i<h->N;i++,length--)
    temp[i] = F77CALL (ddot) (&length,&q_up->data[i],&incm1,&b[i],&inc1);

  for (i=1;i<=h->N;i++)
    x[i-1] = F77CALL (ddot) (&i,&q->data[h->N-i],&inc1,temp,&inc1);

  temp[h->N-1]=0.0;
  for (i=0,length=h->N-1;i<h->N-1;i++,length--)
    temp[i] = F77CALL (ddot) (&length,&q->data[0],&inc1,&b[i+1],&inc1);

  for (i=1;i<h->N;i++)
    x[i] -= F77CALL (ddot) (&i,&q_up->data[0],&incm1,temp,&inc1);

  *temp = 1.0/(*p);
  F77CALL (dscal) (&h->N,temp,x,&inc1);

  free(temp);
}

static void invert_antireg(fdouble *x,toeplitz *h,fdouble *b,
			   polynom *q_lef,polynom *q_up)
{
  fdouble *temp;
  finteger i,length;

  temp=(fdouble *)malloc(sizeof(fdouble)*h->N);
  if (temp==NULL)
  {
    fprintf(stderr,"can't allocate memory for working.\n"
	    "Abort refining ...\n");
    return;
  }

  for (i=0,length=h->N;i<h->N;i++,length--)
    temp[i] = F77CALL (ddot) (&length,&q_up->data[i+1],&incm1,&b[i],&inc1);

  for (i=1;i<=h->N;i++)
    x[i-1] = F77CALL (ddot) (&i,&q_lef->data[h->N-i],&inc1,temp,&inc1);


  temp[h->N-1]=0.0;
  for (i=0,length=h->N-1;i<h->N-1;i++,length--)
    temp[i] = F77CALL (ddot) (&length,&q_lef->data[0],&inc1,&b[i+1],&inc1);

  for (i=1;i<=h->N;i++)
    x[i-1] -= F77CALL (ddot) (&i,&q_up->data[0],&incm1,temp,&inc1);

  free(temp);
}
