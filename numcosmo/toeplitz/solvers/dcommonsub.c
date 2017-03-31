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
#include "polytool.h"
#include "dcommonsub.h"

static finteger inc1  =  1;
static finteger inc0  =  0;
static finteger incm1 = -1;

static finteger info;

static fdouble one  =  1.0;
static fdouble zero =  0.0;
static fdouble mone = -1.0;

#define lapack_err(info,subroutine,laroutine) \
if (info != 0) { \
  fprintf(stderr,"toeplitz solver Error: in subroutine %s %s returns with " \
           "info = %d\n",subroutine,laroutine,(int)info); \
  return 0; \
}

extern void F77CALL (dgesv) (finteger *n, finteger *nrhs, fdouble *a,
      finteger  *lda, finteger *ipiv, fdouble *b, finteger *ldb,
      finteger *info);

extern void F77CALL (dgetrf) (finteger *m, finteger *n, fdouble *a, finteger *
       lda, finteger *ipiv, finteger *info);

extern void F77CALL (dgetrs) (char *trans, finteger *n, finteger *nrhs, 
       fdouble *a, finteger *lda, finteger *ipiv, fdouble *b, finteger *
       ldb, finteger *info);

extern void F77CALL (dgesvd) (char *jobu, char *jobvt, finteger *m,finteger *n,
       fdouble *a, finteger *lda, fdouble *s, fdouble *u, finteger *
       ldu, fdouble *vt, finteger *ldvt, fdouble *work, finteger *lwork, 
       finteger *info);

extern void F77CALL (dorm2r)(char *, char *, finteger *, finteger *, 
       finteger *, fdouble *, finteger *, fdouble *, fdouble *
       , finteger *, fdouble *, finteger *);

extern void F77CALL (dlacpy) (char *uplo, finteger *m, finteger *n, fdouble *a,
       finteger *lda, fdouble *b, finteger *ldb);

extern void F77CALL (dlapmt) (flogical *, finteger *, finteger *, 
        fdouble *, finteger *, finteger *);

extern void F77CALL (dgeqpf) (finteger *, finteger *, fdouble *, finteger *
        , finteger *, fdouble *, fdouble *, finteger *);

int update_x(fdouble *work,finteger lwork,polynom *q,polynom *q_up,
	     finteger *ldq,fdouble *b,fdouble *x,fdouble *d__,
	     finteger *ldd,finteger *d_perm,fdouble *d_tau,finteger n,
	     finteger k)
     /* needs lwork = 4*k workspace */
{
  finteger ihelp=n+k;

#define d(I,J) d__[(I) + (J) * (*ldd)]

  if (k == 1)
  {
    fdouble help = F77CALL (ddot) (&ihelp,q->data,&inc1,b,&inc1)/d(0,0);
    F77CALL (daxpy) (&ihelp,&help,q_up->data,&inc1,x,&inc1);
  }
  else
  {
#ifdef USE_QR_DECOMP
    finteger i;
    fdouble *lawork = &work[k];
    flogical flag=FFALSE;
    lwork -= k;
#ifdef WORK_DEBUG
    if (lwork < 3*k)
      fprintf(stderr,"update_x has got too little workspace\n");
#endif
#else
#ifdef WORK_DEBUG
    lwork -= k;
    if (lwork < 0)
      fprintf(stderr,"update_x has got too little workspace\n");
#endif
#endif

    /* Q' * b */
    F77CALL (dgemv) ("Transpose",&ihelp,&k,&one,q->data,ldq,b,&inc1,&zero,
	  work,&inc1);

    /* D^-1 * Q' * b */
#ifdef USE_QR_DECOMP
    /* get QR decomposition of D */
    for (i=0;i<k;d_perm[i++]=0);   
    F77CALL (dgeqpf) (&k,&k,&d(0,0),ldd,d_perm,d_tau,lawork,&info);
    lapack_err(info,"update_x","dgeqpf");

    /* solve the equation */
    F77CALL (dorm2r)("left","Transpose",&k,&inc1,&k,&d(0,0),ldd,d_tau,
	    work,&k,lawork,&info);
    lapack_err(info,"update_x","dorm2r");

    F77CALL (dtrsv) ("Upper","No transpose","No unit",&k,&d(0,0),ldd,
	  work,&inc1);

    F77CALL (dlapmt) (&flag,&inc1,&k,work,&inc1,d_perm);
#else
    F77CALL (dgesv) (&k,&inc1,&d(0,0),ldd,d_perm,work,&k,&info);
    lapack_err(info,"update_x","dgesv");
#endif

    /* x += Q_up * D^-1 * Q' * b */
    F77CALL (dgemv) ("No transpose",&ihelp,&k,&one,q_up->data,ldq,work,&inc1,
	   &one,x,&inc1);
  }

  return 1;
#undef d
}

fdouble get_singularity(fdouble *work,finteger lwork,fdouble *m,
			   finteger *ldm,finteger k)
/* returns the smallest singular value of m where m is NOT overwritten */
/* requires lwork = k*k + 6*k workspace */
{
  fdouble *temp = work;
  fdouble *sing = &temp[k*k];
  fdouble *lawork = &sing[k];
  lwork -= k*(k+1);
#ifdef WORK_DEBUG
    if (lwork < 5*k)
      fprintf(stderr,"get_singularity has got too little workspace\n");
#endif

  F77CALL (dlacpy) ("All",&k,&k,m,ldm,temp,&k);
  F77CALL (dgesvd) ("N","N",&k,&k,temp,&k,sing,NULL,&inc1,NULL,&inc1,
	lawork,&lwork,&info);
  lapack_err(info,"get_singularity","dgesvd");

  /* take the maximum because dgesvd returns negative singular values
     somtimes */
  return max(sing[k-1],0.0);
}

void update_d(fdouble *d__,finteger *ldd,fdouble *q__,finteger *ldq,
	      fdouble *e_up__,finteger *lde,finteger n,finteger k)
{
  finteger ihelp = k-1;

#define d(I,J)      d__[(I) + (J) * (*ldd)]
#define Q(I,J)      q__[(I) + (J) * (*ldq)]
#define e_up(I,J)   e_up__[(I) + (J) * (*lde)]

  if (ihelp > 0)
    F77CALL (dgemm) ("transpose","no transpose",&inc1,&ihelp,&k,&one,
	  &Q(n,ihelp),ldq,&e_up(0,0),lde,&zero,&d(k-1,0),ldd);

  F77CALL (dcopy) (&k,&e_up(0,k-1),&inc1,&d(0,k-1),&inc1);
  F77CALL (dtrmv) ("Upper","Transpose","No unit",&k,&Q(n,0),ldq,
	&d(0,k-1),&inc1);

#undef d
#undef Q
#undef e_up
}

void update_d_schur(fdouble *d__,finteger *ldd,fdouble *e_up,finteger *lde,
		 fdouble *p,finteger *ldp,finteger k)
{
#define d(I,J)      d__[(I) + (J) * (*ldd)]
  F77CALL (dcopy) (&k,e_up,lde,&d(0,k-1),&incm1);
  if (1 < k)
  {
    finteger ihelp = -(*ldd);
    k--;
    F77CALL (dcopy) (&k,&p[*ldp],ldp,&d(k,0),&ihelp);
  }
#undef d
}

int update_y(fdouble *work,finteger lwork,fdouble *b,fdouble *y,
      fdouble *e_up__,finteger *lde,finteger N,finteger n,finteger k)
     /* k*k+4*k+liwork(k) workspace needed */
{
  finteger ihelp=N-n-k;
  fdouble help;
 
#define e_up(I,J) e_up__[(I) + (J) * (*lde)]
 
  if (k == 1)
  {
    y[0] = b[0]/e_up(0,0);
    help = -y[0];
    F77CALL (daxpy) (&ihelp,&help,&e_up(k,0),&inc1,&b[k],&inc1);
  }
  else
  {
    fdouble *e_up_save = work;
#ifndef USE_QR_DECOMP
    finteger *perm         = (finteger *)&work[k*k];
#ifdef WORK_DEBUG
    lwork -= k*k + liwork(k);
    if (lwork < 0)
      fprintf(stderr,"update_y has got too little workspace\n");
#endif
#else
    flogical flag=FFALSE;
    finteger *perm,i;
    fdouble *tau       = &work[k*k];
    fdouble *lawork    = &tau[k];
    lwork -= k*(k+1)+liwork(k);
    perm = (finteger *)&lawork[lwork];
#ifdef WORK_DEBUG
    if (lwork < 3*k)
      fprintf(stderr,"update_y has got too little workspace\n");
#endif
#endif

    /* save diagonal block of e_up */
    F77CALL (dlacpy) ("All",&k,&k,&e_up(0,0),lde,e_up_save,&k);
    F77CALL (dcopy) (&k,b,&inc1,y,&inc1);

    /* E_up^-1 * b */
#ifdef USE_QR_DECOMP
    /* get QR decomposition of D */

    for (i=0;i<k;perm[i++]=0);

    F77CALL (dgeqpf) (&k,&k,e_up_save,&k,perm,tau,lawork,&info);
    lapack_err(info,"update_y","dgeqpf");
 
    /* solve the equation */
    F77CALL (dorm2r)("left","Transpose",&k,&inc1,&k,e_up_save,&k,tau,
            y,&k,lawork,&info);
    lapack_err(info,"update_y","dorm2r");
 
    F77CALL (dtrsv) ("Upper","No transpose","No unit",&k,e_up_save,&k,y,&inc1);

    F77CALL (dlapmt) (&flag,&inc1,&k,y,&inc1,perm);
#else
    F77CALL (dgesv) (&k,&inc1,e_up_save,&k,perm,y,&k,&info);
    lapack_err(info,"update_y","dgesv");
#endif
 
    /* b -= E_up*y  */
    F77CALL (dgemv) ("No transpose",&ihelp,&k,&mone,&e_up(k,0),lde,y,&inc1,
           &one,&b[k],&inc1);
  }
 
  return 1;
#undef e_up
}


int update_Dy_schur(fdouble *work,finteger lwork,fdouble *b,fdouble *y,
      fdouble *e_up__,fdouble *d,finteger *ldd,finteger N,finteger n,
      finteger k,finteger *ldq)
     /* k*k+4*k+liwork(k) workspace needed */
{

  /* computes  y = D * E_up^{-1} * b */

  finteger ihelp=N-n-k;
  fdouble help;

#define e_up(I,J) e_up__[(I) + (J) * (*ldq)]

  if (k == 1)
  {
    y[0] = b[0]/e_up(0,0);
    help = -y[0];
    /* update right-hand side */
    F77CALL (daxpy) (&ihelp,&help,&e_up(k,0),&inc1,&b[k],&inc1);
    /* multiply by D */
    y[0] = y[0]*d[0];
  }
  else
  {
    fdouble *e_up_save=work;
#ifndef USE_QR_DECOMP
    finteger *perm=(finteger*)&e_up_save[k*k];
#ifdef WORK_DEBUG
    lwork -= k*k + liwork(k);
    if (lwork < 0)
      fprintf(stderr,"update_y_schur has got too little workspace\n");
#endif
#else
    flogical flag=FFALSE;
    finteger i;
    finteger *perm;
    fdouble *tau=&e_up_save[k*k];
    fdouble *lawork=&tau[k];
    lwork-=k*(k+1)+liwork(k);
    perm=(finteger*)&lawork[lwork];
#ifdef WORK_DEBUG
    if (lwork < 3*k)
      fprintf(stderr,"update_y_schur has got too little workspace\n");
#endif
#endif

    /* save diagonal block of E_up */
    F77CALL (dlacpy) ("All",&k,&k,&e_up(0,0),ldq,e_up_save,&k);
    F77CALL (dcopy) (&k,b,&inc1,y,&inc1);
    /* E_up^-1 * b */
#ifdef USE_QR_DECOMP
    /* get QR decomposition of E_up */
    for (i=0;i<k;perm[i++]=0);   
    F77CALL (dgeqpf) (&k,&k,e_up_save,&k,perm,tau,lawork,&info);
    lapack_err(info,"update_y_schur","dgeqpf");

    /* solve the equation */
    F77CALL (dorm2r)("left","Transpose",&k,&inc1,&k,e_up_save,&k,tau,
	    y,&k,lawork,&info);
    lapack_err(info,"update_y_schur","dorm2r");

    F77CALL (dtrsv) ("Upper","No transpose","No unit",&k,e_up_save,&k,y,&inc1);

    F77CALL (dlapmt) (&flag,&inc1,&k,y,&inc1,perm);
#else
    F77CALL (dgesv) (&k,&inc1,e_up_save,&k,perm,y,&k,&info);
    lapack_err(info,"update_y_schur","dgesv");
#endif

    /* b -= E_up*y  */
    F77CALL (dgemv) ("No transpose",&ihelp,&k,&mone,&e_up(k,0),ldq,y,&inc1,
	   &one,&b[k],&inc1);

    /* multiply by D, store y in e_up_save  */
    F77CALL (dcopy) (&k,y,&inc1,e_up_save,&inc1);
    F77CALL (dgemv) ("No transpose",&k,&k,&one,d,ldd,e_up_save,&inc1,
	   &zero,y,&inc1);
  }

  return 1;
#undef e_up
}

int update_y_schur(fdouble *work,finteger lwork,fdouble *b,fdouble *y,
      fdouble *e_up__,finteger N,finteger n,finteger k,finteger *ldq)
     /* k*k+4*k+liwork(k) workspace needed */
{

  /* computes  y = E_up^{-1} * b */

  finteger ihelp=N-n-k;
  fdouble help;

#define e_up(I,J) e_up__[(I) + (J) * (*ldq)]

  if (k == 1)
  {
    y[0] = b[0]/e_up(0,0);
    help = -y[0];
    /* update right-hand side */
    F77CALL (daxpy) (&ihelp,&help,&e_up(k,0),&inc1,&b[k],&inc1);
  }
  else
  {
    fdouble *e_up_save=work;
#ifndef USE_QR_DECOMP
    finteger *perm=(finteger*)&e_up_save[k*k];
#ifdef WORK_DEBUG
    lwork -= k*k + liwork(k);
    if (lwork < 0)
      fprintf(stderr,"update_y_schur has got too little workspace\n");
#endif
#else
    flogical flag=FFALSE;
    finteger i;
    finteger *perm;
    fdouble *tau=&e_up_save[k*k];
    fdouble *lawork=&tau[k];
    lwork-=k*(k+1)+liwork(k);
    perm=(finteger*)&lawork[lwork];
#ifdef WORK_DEBUG
    if (lwork < 3*k)
      fprintf(stderr,"update_y_schur has got too little workspace\n");
#endif
#endif

    /* save diagonal block of E_up */
    F77CALL (dlacpy) ("All",&k,&k,&e_up(0,0),ldq,e_up_save,&k);
    F77CALL (dcopy) (&k,b,&inc1,y,&inc1);
    /* E_up^-1 * b */
#ifdef USE_QR_DECOMP
    /* get QR decomposition of E_up */
    for (i=0;i<k;perm[i++]=0);   
    F77CALL (dgeqpf) (&k,&k,e_up_save,&k,perm,tau,lawork,&info);
    lapack_err(info,"update_y_schur","dgeqpf");

    /* solve the equation */
    F77CALL (dorm2r)("left","Transpose",&k,&inc1,&k,e_up_save,&k,tau,
	    y,&k,lawork,&info);
    lapack_err(info,"update_y_schur","dorm2r");

    F77CALL (dtrsv) ("Upper","No transpose","No unit",&k,e_up_save,&k,y,&inc1);

    F77CALL (dlapmt) (&flag,&inc1,&k,y,&inc1,perm);
#else
    F77CALL (dgesv) (&k,&inc1,e_up_save,&k,perm,y,&k,&info);
    lapack_err(info,"update_y_schur","dgesv");
#endif

    /* b -= E_up*y  */
    F77CALL (dgemv) ("No transpose",&ihelp,&k,&mone,&e_up(k,0),ldq,y,&inc1,
	   &one,&b[k],&inc1);
   }

  return 1;
#undef e_up
}

int update_x_from_y(polynom *q_up,fdouble *y,fdouble *x,
                    finteger n,finteger k,finteger *ldq)
{
  finteger ihelp=n;

  if (k == 1)
  {
    /* x(0:n-1) += Q_up(0:n-1,n) * y(n) */
    F77CALL (daxpy) (&ihelp,y,q_up->data,&inc1,x,&inc1);
    x[n] = y[0];
  }
  else
  {
    /* x(0:n-1) += Q_up(0:n-1,n:n+k) * y(n:n+k) */
    ihelp = n+k;
    F77CALL (dgemv) ("No transpose",&ihelp,&k,&one,q_up->data,ldq,y,&inc1,
           &one,x,&inc1);
  }

  return 1;
}

int update_x_from_Dy_schur(fdouble *work,finteger orglwork,fdouble *L__,
      fdouble *y,fdouble *x,finteger N,finteger *reg_index,
      finteger reg_pos,finteger *ldq)
     /* k+liwork(k) workspace needed */

{
  /* computes  x = L^{-T} * y */

  /* attention: this variant uses the transpose */
#define Lt(I,J) L__[(J) + (I) * (*ldq)]

  finteger n,k;
  fdouble help;

  for (;reg_pos>0;reg_pos--)
  {
    n = reg_index[reg_pos-1];
    k = reg_index[reg_pos]-n;

    if (k == 1)
    {
      x[n] = y[n]/Lt(n,n);
      help = -x[n];
      /* update right-hand side */
      F77CALL (daxpy) (&n,&help,&Lt(0,n),ldq,y,&inc1);
    }
    else
    {
#ifndef USE_QR_DECOMP
      finteger *perm=(finteger*)work;
#ifdef WORK_DEBUG
      finteger lwork = orglwork - liwork(k);
      if (lwork < 0)
	fprintf(stderr,"update_x_from_y has got too little workspace\n");
#endif
#else
      fdouble *tau=work;
      fdouble *lawork = &work[k];
      finteger lwork = orglwork - k - liwork(k);
      finteger *perm=(finteger*)&lawork[lwork];
      finteger i;
      flogical flag=FTRUE;
#ifdef WORK_DEBUG
      if (lwork < 3*k)
	fprintf(stderr,"update_x_from_y has got too little workspace\n");
#endif
#endif
      F77CALL (dcopy) (&k,&y[n],&inc1,&x[n],&inc1);

      /* Lt^-1 * b */
#ifdef USE_QR_DECOMP
      /* get QR decomposition of Lt */
      for (i=0;i<k;perm[i++]=0);
      F77CALL (dgeqpf) (&k,&k,&Lt(n,n),ldq,perm,tau,lawork,&info);
      lapack_err(info,"update_x_from_Dy_schur","dgeqpf");

      F77CALL (dlapmt) (&flag,&inc1,&k,&x[n],&inc1,perm);

      /* solve the equation */
      F77CALL (dtrsv) ("Upper","transpose","No unit",&k,&Lt(n,n),ldq,
	    &x[n],&inc1);

      F77CALL (dorm2r)("left","No transpose",&k,&inc1,&k,&Lt(n,n),ldq,tau,
	    &x[n],&k,lawork,&info);
      lapack_err(info,"update_x_from_Dy_schur","dorm2r");
#else
      F77CALL (dgetrf) (&k,&k,&Lt(n,n),ldq,perm,&info);
      F77CALL (dgetrs) ("Transpose",&k,&inc1,&Lt(n,n),ldq,perm,&x[n],ldq,
	    &info);
      lapack_err(info,"update_x_from_Dy_schur","dgesv");
#endif
    /* y -= Lt*x  */
      F77CALL (dgemv) ("Transpose",&k,&n,&mone,&Lt(0,n),ldq,&x[n],&inc1,
	     &one,y,&inc1);
    }
  }
  return 1;
}

int update_x_from_y_schur(fdouble *work,finteger orglwork,fdouble *schur_par,
      fdouble *y_,fdouble *x_,finteger N,finteger *reg_index,finteger reg_pos,
      finteger *ldq)
{
  /* computes  x = R_up * y from stored Schur parameters */
  finteger i,n,ihelp;
  finteger k=0;
  int switch_to_cl;

  fdouble *xtemp = work;
  fdouble *ytemp = xtemp+(*ldq);
  fdouble *y = ytemp+(*ldq);
  fdouble *x = y+(*ldq);

  F77CALL (dcopy) (&N,y_,&incm1,y,&inc1);
  F77CALL (dcopy) (&N,&zero,&inc0,x,&inc1);

  reg_pos--;
  /* if last index was column-regular, last coefficients not required */
  if (reg_index[reg_pos] == N-1)
    reg_pos--;

  for (;reg_pos>=0;reg_pos--)
  {
    n = reg_index[reg_pos];
    k = reg_index[reg_pos+1]-n;

    /* check if we have a switch from look-ahead to classical step  */
    /* then: parameters alpha_0, beta_0 are stored in h_alpha, h_beta */
    /* have to apply two transformations! */
    if (reg_pos > 0)
      switch_to_cl = (k == 1 && reg_index[reg_pos-1] < n-1);
    else
      switch_to_cl = 0;

#if (DEBUG >= 4)
      printf("    n=%3d, k=%3d, switch_to_cl = %d\n",(int)n,(int)k,
	    switch_to_cl);
#endif

    if (k == 1)
    {
#if (DEBUG >= 2)
      printf(" n=%3d, k=%3d, classical (8.44)\n",(int)n,(int)k);
#endif
#if (DEBUG >= 5)
      printf("   CR_GAMMA=%e, CR_GAMMA_UP=%e\n",CR_GAMMA(n),CR_GAMMA_UP(n));
#endif
      
      ihelp = N-n-1;
      F77CALL (dcopy) (&ihelp,x+n+1,&inc1,xtemp+n+1,&inc1);
      F77CALL (dcopy) (&ihelp,y,&inc1,ytemp,&inc1);
      F77CALL (daxpy) (&ihelp,&CR_GAMMA_UP(n),ytemp,&inc1,x+n+1,&inc1);

      if (CR_GAMMA(n) != 0.0)
	F77CALL (daxpy) (&ihelp,&CR_GAMMA(n),xtemp+n+1,&inc1,y,&inc1);

      if (switch_to_cl)
      {
#if (DEBUG >= 1)
	printf("n=%3d, apply transformation for switch, ",n);
#endif
#if (DEBUG >=5)
	printf("alpha=%e, beta=%e\n",ALPHA(n),BETA(n));
#endif

	ihelp = N-n;
	F77CALL (dcopy) (&ihelp,x+n,&inc1,xtemp+n,&inc1);
	F77CALL (dcopy) (&ihelp,y,&inc1,ytemp,&inc1);

	F77CALL (dscal) (&ihelp,&ALPHA(n),x+n,&inc1);

	F77CALL (daxpy) (&ihelp,&BETA(n),xtemp+n,&inc1,y,&inc1);
      }
    }      
    else
    {
      finteger len = N-n-k;
#if (DEBUG >= 1)
      printf("block of size %3d, n=%3d, take parameter on position %3d\n",
	    (int)k,(int)n,(int)(n+1));
#endif
#if (DEBUG >=5)
      printf("ALPHA_LEF=%e, BETA_LEF=%e\n",ALPHA_LEF(n),BETA_LEF(n));
#endif

      ihelp = N-n+1;
      F77CALL (dcopy) (&ihelp,x+n,&inc1,xtemp+n,&inc1);
      F77CALL (dcopy) (&ihelp,y,&inc1,ytemp,&inc1);

      F77CALL (dscal) (&len,&ALPHA_LEF(n),x+n+k,&inc1);
      for (i=1;i<k-1;i++)
	F77CALL (daxpy) (&len,&ALPHA_LEF(n+i),xtemp+n+k,&inc1,x+n+k-i,&inc1);

      for (i=0;i<k;i++)
	F77CALL (daxpy) (&len,&ALPHA_UP(n+i),ytemp,&inc1,x+n+k-i,&inc1);

/*       ihelp = k-1; */
/*       for (i=0;i<k-1;i++,ihelp--) */
/* 	F77CALL (daxpy)(&ihelp,&GAMMA_UP(n+i),ytemp+len,&inc1,x+N-ihelp,&inc1); */
      for (ihelp=k-1;ihelp>0;ihelp--)
	F77CALL (daxpy)(&ihelp,&GAMMA_UP(n+k-1-ihelp),ytemp+len,&inc1,
	      x+N-ihelp,&inc1);

      for (i=0;i<k;i++)
	F77CALL (daxpy) (&len,&BETA_UP(n+i),ytemp,&inc1,y+k-i,&inc1);

      for (i=0;i<k;i++)
	F77CALL (daxpy) (&len,&BETA_LEF(n+i),xtemp+n+k,&inc1,y+k-i,&inc1);
    }
  }

  if (k == 1)
    F77CALL (daxpy) (&N,&one,x,&inc1,y,&inc1);
  F77CALL (dcopy) (&N,y,&incm1,x_,&inc1);

  return 1;
}

void get_M(fdouble *M__,finteger *ldM,fdouble *e_lef,
	    fdouble *e_up,fdouble *p_lef,
	    fdouble *p_up,finteger k)
{
  finteger ihelp=(*ldM)*(*ldM);
  finteger i;


#define M(I,J)      M__[(I) + (J) * (*ldM)]

  F77CALL (dcopy) (&ihelp,&zero,&inc0,M__,&inc1);

  for (i=1;i<k;i++)
  {
    F77CALL (dcopy) (&i,p_lef,&incm1,&M(0,i-1),&inc1);
    F77CALL (dcopy) (&i,p_up,&incm1,&M(0,k+i-1),&inc1);
    ihelp = k-i+1;
    F77CALL (dcopy) (&ihelp,e_lef,&inc1,&M(k-2+i,i-1),&inc1);
  }
  for (i=k;i>0;i--)
  {
    F77CALL (dcopy) (&i,e_up,&inc1,&M(2*k-i-1,2*k-i-1),&inc1);
  }
#undef M
}

void get_colreg_pair(polynom *res,polynom *res_up,toeplitz *h,polynom *q,
		     polynom *q_up,fdouble *e,fdouble *e_up,
		     fdouble *p,fdouble *p_up,fdouble *v,
		     fdouble *v_up,finteger n)
{
  /* process variables */
  *v=-(*e)/(*p);
  *v_up=-(*p_up)/(*p);


  /* calculating new q and new q_up */
  /* copy q into res and shift the position */
  res->data[0]=0.0;
  F77CALL (dcopy) (&q->length,q->data,&q->stdinc,&(res->data[1]),&res->stdinc);
  res->length=q->length+1;

  /* copy q_up into res_up and shift the position */
  res_up->data[0]=0.0;
  F77CALL (dcopy) (&q_up->length,q_up->data,&q_up->stdinc,&(res_up->data[1]),
	 &res_up->stdinc);
  res_up->length=q_up->length+1;

  /* q_up = [0 q_up] + v_up * [q 0] */
  F77CALL (daxpy) (&q->length,v_up,q->data,&q->stdinc,
	res_up->data,&res_up->stdinc);

  /* q = [0 q_up] * v + [q 0] */
  F77CALL (daxpy) (&q_up->length,v,q_up->data,&q_up->stdinc,
	res->data,&res->stdinc);


  /* updating pi_0 */
  (*p)*=(1.0-(*v_up)*(*v));

  /* calculating the pi_up_0 */
  (*p_up)=lau_coef(h,res_up,-1);

}

void get_colreg_pair_schurlev(polynom *res,polynom *res_up,toeplitz *h,
     polynom *q,polynom *q_up,fdouble *e,fdouble *e_up,fdouble *eres,
     fdouble *eres_up,fdouble *p,fdouble *p_up,fdouble *pres,
     fdouble *pres_up,finteger *ldsmall,fdouble *v, fdouble *v_up,
     finteger n)
{
  finteger ihelp;

  /* process variables */
  *v=-(*e)/(*p);
  *v_up=-(*p_up)/(*p);

  if (res != NULL)
  {
    /* calculating new q and new q_up */
    /* copy q into res and shift the position */
    res->data[0]=0.0;
    F77CALL (dcopy) (&q->length,q->data,&q->stdinc,
	  &(res->data[1]),&res->stdinc);
    res->length=q->length+1;

    /* copy q_up into res_up and shift the position */
    res_up->data[0]=0.0;
    F77CALL (dcopy) (&q_up->length,q_up->data,&q_up->stdinc,&(res_up->data[1]),
	   &res_up->stdinc);
    res_up->length=q_up->length+1;

    /* q_up = [0 q_up] + v_up * [q 0] */
    F77CALL (daxpy) (&q->length,v_up,q->data,&q->stdinc,
	  res_up->data,&res_up->stdinc);

    /* q = [0 q_up] * v + [q 0] */
    F77CALL (daxpy) (&q_up->length,v,q_up->data,&q_up->stdinc,
	  res->data,&res->stdinc);
  }

  /* copy p into pres */
  ihelp = h->N-n;
  F77CALL (dcopy) (&ihelp,p,&inc1,&pres[1],&inc1);
  F77CALL (dcopy) (&ihelp,e_up,&inc1,&eres_up[1],&inc1);
  ihelp = h->N-n-1;
  F77CALL (dcopy) (&ihelp,&p_up[1],&inc1,&pres_up[1],&inc1);
  F77CALL (dcopy) (&ihelp,&e[1],&inc1,&eres[1],&inc1);

  /* updating p,p_up,e,e_up */
  ihelp = h->N-n;
  F77CALL (daxpy) (&ihelp,v,p_up,&inc1,&pres[1],&inc1);
  F77CALL (daxpy) (&ihelp,v_up,e,&inc1,&eres_up[1],&inc1);

  ihelp = h->N-n-1;
  F77CALL (daxpy) (&ihelp,v_up,&p[1],&inc1,&pres_up[1],&inc1);
  F77CALL (daxpy) (&ihelp,v,&e_up[1],&inc1,&eres[1],&inc1);
}

void get_inner_q_up(polynom *res_up,polynom *q_lef,polynom *q_up,
		    fdouble *v_up,finteger k)
{
  /* res_up = q_up *z + v_up * q_lef */
  mult_by_z_up_k(res_up,q_up,&one,1);

  polynomaxpy(v_up,q_lef,res_up);
}

void get_inner_x_up_schur(fdouble *eres_up,fdouble *e_lef,
     fdouble *e_up,fdouble *pres_up,fdouble *p_lef,fdouble *p_up,
       fdouble *v_up,finteger N,finteger n,finteger k)
{
  finteger ihelp = N-n-k;

  if (pres_up != NULL)
  {
    F77CALL (dcopy) (&ihelp,&p_up[1],&inc1,pres_up,&inc1);
    F77CALL (daxpy) (&ihelp,v_up,&p_lef[1],&inc1,pres_up,&inc1);
  }

  eres_up[0] = 0.0;
  ihelp = N-n;

  F77CALL (dcopy) (&ihelp,e_up,&inc1,&eres_up[1],&inc1);
  F77CALL (daxpy) (&ihelp,v_up,e_lef,&inc1,eres_up,&inc1);
}

void get_inner_q_up_schurlev(polynom *res_up,polynom *q_lef,polynom *q_up,
       fdouble *eres_up,fdouble *e_lef,fdouble *e_up,
       fdouble *pres_up,fdouble *p_lef,fdouble *p_up,
       fdouble *v_up,finteger N,finteger n,finteger k)

{
  if (res_up != NULL)
    get_inner_q_up(res_up,q_lef,q_up,v_up,k);
  get_inner_x_up_schur(eres_up,e_lef,e_up,pres_up,p_lef,p_up,v_up,N,n,k);
}

void get_inner_q(polynom *res,polynom *q_lef,polynom *q,fdouble *v,
		 finteger k)
{
  /* inner q's do not have full degree */
  *(res->data)=0.0;
  (res->data)++;

  /* res = q + v * q_lef * z^k */
  mult_by_z_up_k(res,q_lef,v,k);

  polynomaxpy(&one,q,res);

  /* inner q's do not have full degree */
  (res->data)--;
  res->length++;
}

void get_inner_x_schur(fdouble *eres,fdouble *e_lef,fdouble *e,
       fdouble *pres,fdouble *p_lef,fdouble *p,
       fdouble *v,finteger N,finteger n,finteger k)
{
  finteger ihelp = N-n-k;

  if (eres != NULL)
  {
    F77CALL (dcopy) (&ihelp,&zero,&inc0,eres,&inc1);

    ihelp--;
    F77CALL (dcopy) (&ihelp,&e[1],&inc1,eres,&inc1);
    F77CALL (daxpy) (&ihelp,v,&e_lef[1],&inc1,eres,&inc1);
  }

  pres[0] = 0.0;
  ihelp = N-n;
  F77CALL (dcopy) (&ihelp,p,&inc1,&pres[1],&inc1);
  F77CALL (daxpy) (&ihelp,v,p_lef,&inc1,pres,&inc1);
}

void get_inner_q_schurlev(polynom *res,polynom *q_lef,polynom *q,
       fdouble *eres,fdouble *e_lef,fdouble *e,
       fdouble *pres,fdouble *p_lef,fdouble *p,
       fdouble *v,finteger N,finteger n,finteger k)

{
  if (res != NULL)
    get_inner_q(res,q_lef,q,v,k);
  get_inner_x_schur(eres,e_lef,e,pres,p_lef,p,v,N,n,k);
}

void get_new_coefs(fdouble *e__,fdouble *e_up__,fdouble *p__,
		   fdouble *p_up__,fdouble *p_lef,fdouble *e_lef,
		   finteger *ldsmall,toeplitz *h,polynom *q,polynom *q_up,
		   polynom *q_lef,fdouble *v,fdouble *v_up,
		   finteger k,finteger n)
{
  finteger i;

#define p_up(I,J)   p_up__[(I) + (J) * (*ldsmall)]
#define p(I,J)      p__[(I) + (J) * (*ldsmall)]
#define e_up(I,J)   e_up__[(I) + (J) * (*ldsmall)]
#define e(I,J)      e__[(I) + (J) * (*ldsmall)]

  /*     e(k,0)   =lau_coef(h,&q[0],n+k+1);  */
  e(k,0) = (e_lef[k+1]*e_up(0,0) - e_up(k+1,0)*e_lef[0])/p_lef[0];
  /*     p(k,0)   =lau_coef(h,&q[0],-k); */
  p(k,0) = (p_lef[k]*e_up(0,0) - p_up(k-1,0)*e_lef[0])/p_lef[0];

  for (i=1;i<=k;i++)
  {    
    /*     p(k,i)   =lau_coef(h,&q[i],i-k); */
    p(k,i) = p(k-1,i-1) + p_lef[k]* v[i];

    /*     e(k,i)   =lau_coef(h,&q[i],n+k+1); */
    e(k,i) = e(k,i-1) + e_lef[k+1-i]* v[i];

    /*     e_up(k+1,i)=lau_coef(h,&q_up[i],n+k+1); */
    e_up(k+1,i) = e_up(k,i-1) + e_lef[k+1]* v_up[i];

    /*     p_up(k,i)=lau_coef(h,&q_up[i],i-k-1); */
    p_up(k,i) = p_up(k,i-1) + p_lef[k+1-i]* v_up[i];

  }
  for (i=0;i<k;i++)
  {
    /*     p(i,k)=lau_coef(h,&q[k],k-i); */
    /*     e_up(i,k)=lau_coef(h,&q_up[k],n+i); */
    if (i==0)
    {
      p(i,k) = p_lef[i]* v[k];
      e_up(i,k) = e_lef[i]* v_up[k];
    }
    else
    {
      p(i,k) = p(i-1,k-1) + p_lef[i]* v[k];
      e_up(i,k) = e_up(i-1,k-1) + e_lef[i]* v_up[k];
    }
  }
  e_up(k,k) = e_up(k-1,k-1) + e_lef[k]* v_up[k];
      
#undef p_up
#undef p
#undef e_up
#undef e
}

void get_new_e_pup_coefs(fdouble *e__,
		   fdouble *p_up__,fdouble *p_lef,fdouble *e_lef,
		   finteger *ldsmall,fdouble *v,fdouble *v_up,
		   finteger k,finteger n)
{
  finteger i;

#define p_up(I,J)   p_up__[(I) + (J) * (*ldsmall)]
#define p(I,J)      p__[(I) + (J) * (*ldsmall)]
#define e_up(I,J)   e_up__[(I) + (J) * (*ldsmall)]
#define e(I,J)      e__[(I) + (J) * (*ldsmall)]

  for (i=1;i<=k;i++)
  {    
    e(k,i) = e(k,i-1) + e_lef[k+1-i]* v[i];
    p_up(k,i) = p_up(k,i-1) + p_lef[k+1-i]* v_up[i];
  }
      
#undef p_up
#undef p
#undef e_up
#undef e
}

void get_anti_reg_coefs(fdouble *e__,fdouble *e_up__,fdouble *p__,
		   fdouble *p_up__,fdouble *p_lef,fdouble *e_lef,
		   finteger *ldsmall,toeplitz *h,polynom *q,polynom *q_up,
		   polynom *q_lef,fdouble *v,fdouble *v_up,
		   finteger k,finteger n)
{
#define p_up(I,J)   p_up__[(I) + (J) * (*ldsmall)]
#define p(I,J)      p__[(I) + (J) * (*ldsmall)]
#define e_up(I,J)   e_up__[(I) + (J) * (*ldsmall)]
#define e(I,J)      e__[(I) + (J) * (*ldsmall)]


  e_up(k+1,0)=lau_coef(h,&q_up[0],n+k+1);

/*   e_lef[k+1]=lau_coef(h,q_lef,n+k+1); */
  e_lef[k+1] = get_e_lef_k(q_lef,q_up,e_lef,&e_up(0,0),ldsmall,k+1);

  p_up(k,0)=lau_coef(h,&q_up[0],-k-1);

  /*   p_lef[k]=lau_coef(h,q_lef,-k); */
  p_lef[k] = get_p_lef_km1(q_lef,q_up,p_lef,&p_up(0,0),ldsmall,k+1);

#undef p_up
#undef p
#undef e_up
#undef e
}

int get_anti_regular_up(fdouble *work,finteger lwork,polynom *q_lef,
			polynom *q_up,polynom *q,fdouble *e_up__,
			fdouble *e_lef,fdouble *p_up__,fdouble *d__,
			finteger *d_perm,fdouble *d_tau,fdouble *p_lef,
			finteger *ldsmall,finteger *ldq,finteger n,finteger k)
     /* requires ldq + 4*k workspace */

{
#define e_up(I,J)   e_up__[(I) + (J) * (*ldsmall)]
#define p_up(I,J)   p_up__[(I) + (J) * (*ldsmall)]
#define d(I,J)      d__[(I) + (J) * (*ldsmall)]

  polynom q_temp;
  fdouble theta;
  finteger ihelp;
  /* save e_up(k-1,k-1) */
  fdouble help=e_up(k-1,k-1);
  fdouble *e_up_save=&work[*ldq];
#ifdef USE_QR_DECOMP
  fdouble *lawork = &e_up_save[k];
  finteger i;
  flogical flag=FFALSE;
  lwork -=*ldq+k;
#ifdef WORK_DEBUG
  if (lwork < 3*k)
    fprintf(stderr,"get_anti_regular_up has got too little workspace\n");
#endif
#else
#ifdef WORK_DEBUG
  lwork -=*ldq+k;
  if (lwork < 0)
    fprintf(stderr,"get_anti_regular_up has got too little workspace\n");
#endif
#endif

  q_temp.data=work;
  q_temp.stdinc=1;

  /* save k-th row of e_up, since e_up is overwritten by dgesv later */
  F77CALL (dcopy) (&k,&e_up(k,0),ldsmall,e_up_save,&inc1);

  /* calculating anti-regular q_up */
  theta = p_up(k-1,k-1)/p_lef[0];
  F77CALL (dscal) (&k,&theta,e_lef,&inc1);
  
  ihelp=k-1;
  F77CALL (daxpy) (&ihelp,&mone,&e_up(0,k-1),&inc1,&e_lef[1],&inc1);

  /* solve linear system to compute recurrence coefficients */
  /* save solution in e_lef[0:k-1] */
  if (d__ != NULL)
  {
    F77CALL (dtrmv) ("Upper","Transpose","No unit",&k,&q->data[n],ldq,
	  e_lef,&inc1);

#ifdef USE_QR_DECOMP
    /* solve the equation */
    F77CALL (dorm2r)("left","Transpose",&k,&inc1,&k,&d(0,0),ldsmall,d_tau,
	    e_lef,ldsmall,lawork,&info);
    lapack_err(info,"get_anti_regular_up","dorm2r");

    F77CALL (dtrsv) ("Upper","No transpose","No unit",&k,&d(0,0),ldsmall,
	  e_lef,&inc1);

    F77CALL (dlapmt) (&flag,&inc1,&k,e_lef,&inc1,d_perm);
#else
    F77CALL (dgetrs) ("No transpose",&k,&inc1,&d(0,0),ldsmall,d_perm,
	  e_lef,ldsmall,&info);
    lapack_err(info,"get_anti_regular_up","dgetrs");
#endif
  }
  else
  {
#ifdef USE_QR_DECOMP
    /* get QR decomposition of D */
    for (i=0;i<k;d_perm[i++]=0);
    F77CALL (dgeqpf) (&k,&k,&e_up(0,0),ldsmall,d_perm,d_tau,lawork,&info);
    lapack_err(info,"get_anti_regular_up","dgeqpf");

    /* solve the equation */
    F77CALL (dorm2r)("left","Transpose",&k,&inc1,&k,&e_up(0,0),ldsmall,d_tau,e_lef,
	    ldsmall,lawork,&info);
    lapack_err(info,"get_anti_regular_up","dorm2r");

    F77CALL (dtrsv) ("Upper","No transpose","No unit",&k,&e_up(0,0),ldsmall,
	  e_lef,&inc1);

    F77CALL (dlapmt) (&flag,&inc1,&k,e_lef,&inc1,d_perm);
#else
    F77CALL (dgesv) (&k,&inc1,&e_up(0,0),ldsmall,d_perm,e_lef,&k,&info);
    lapack_err(info,"get_anti_regular_up","dgesv");
#endif
  }

  mult_by_z_up_k(&q_temp,&q_up[k-1],&one,1);
      
  ihelp = n+k;
  F77CALL (dgemv) ("No transpose",&ihelp,&k,&one,q_up[0].data,ldq,e_lef,&inc1,
	 &one,q_temp.data,&q_temp.stdinc);

  theta = -theta;
  polynomaxpy(&theta,q_lef,&q_temp);

  q_up[0].length = q_temp.length;
  F77CALL (dcopy) (&q_temp.length,q_temp.data,&q_temp.stdinc,
	 q_up[0].data,&q_up[0].stdinc);

  /* e_lef[k] still contains the correct value of e_lef */
  e_up(0,0) = help + F77CALL (ddot) (&k,e_lef,&inc1,e_up_save,&inc1)
         + theta*e_lef[k];

  return 1;
#undef e_up
#undef p_up
#undef d
}

int get_anti_regular_up_schurlev(fdouble *work,finteger lwork,fdouble *theta_,
      polynom *q_lef,polynom *q_up,polynom *q,fdouble *e_up__,fdouble *e_lef,
      fdouble *p_up__,fdouble *p_lef,fdouble *d__,finteger *d_perm,
      fdouble *d_tau,finteger *ldsmall,finteger *ldq,finteger N,finteger n,
      finteger k)
     /* lwork > ldq+k*k+3*k */
{
#define e_up(I,J)   e_up__[(I) + (J) * (*ldq)]
#define p_up(I,J)   p_up__[(I) + (J) * (*ldq)]
#define d(I,J)      d__[(I) + (J) * (*ldsmall)]

  polynom q_temp;
  fdouble theta;
  fdouble *e_up_save=&work[*ldq];
  finteger ihelp,ldp1;
#ifdef USE_QR_DECOMP
  fdouble *lawork=&e_up_save[k*k];
  flogical flag=FFALSE;
  finteger i;
  lwork -= *ldq+k*k;
#ifdef WORK_DEBUG
  if (lwork < 3*k)
    fprintf(stderr,"get_anti_regular_up_schurlev has got too little workspace"
	    "\n");
#endif
#else
#ifdef WORK_DEBUG
  lwork -= *ldq+k*k;
  if (lwork < 0)
    fprintf(stderr,"get_anti_regular_up_schurlev has got too little workspace"
	    "\n");
#endif
#endif

  q_temp.data=work;
  q_temp.stdinc=1;

  /* calculating anti-regular q_up */
  theta = p_up(k-1,k-1)/p_lef[0];
  F77CALL (dscal) (&k,&theta,e_lef,&inc1);
  
  ihelp=k-1;
  F77CALL (daxpy) (&ihelp,&mone,&e_up(0,k-1),&inc1,&e_lef[1],&inc1);

  /* solve linear system to compute recurrence coefficients */
  /* save solution in e_lef[0:k-1] */
  if (d__ != NULL)
  {
    F77CALL (dtrmv) ("Upper","Transpose","No unit",&k,&q->data[n],ldq,
	  e_lef,&inc1);

#ifdef USE_QR_DECOMP
    /* solve the equation */
    F77CALL (dorm2r)("left","Transpose",&k,&inc1,&k,&d(0,0),ldsmall,d_tau,
	    e_lef,ldsmall,lawork,&info);
    lapack_err(info,"get_anti_regular_up_schurlev","dorm2r");

    F77CALL (dtrsv) ("Upper","No transpose","No unit",&k,&d(0,0),ldsmall,
	  e_lef,&inc1);

    F77CALL (dlapmt) (&flag,&inc1,&k,e_lef,&inc1,d_perm);
#else
    F77CALL (dgetrs) ("No transpose",&k,&inc1,&d(0,0),ldsmall,d_perm,
	  e_lef,ldsmall,&info);
    lapack_err(info,"get_anti_regular_up_schurlev","dgetrs");
#endif
  }
  else
  {
    F77CALL (dlacpy) ("A",&k,&k,&e_up(0,0),ldq,e_up_save,&k);

#ifdef USE_QR_DECOMP
    /* get QR decomposition of D */
    for (i=0;i<k;d_perm[i++]=0);
    F77CALL (dgeqpf) (&k,&k,e_up_save,&k,d_perm,d_tau,lawork,&info);
    lapack_err(info,"get_anti_regular_up_schurlev","dgeqpf");

    /* solve the equation */
    F77CALL (dorm2r)("left","Transpose",&k,&inc1,&k,e_up_save,&k,d_tau,
	    e_lef,ldsmall,lawork,&info);
    lapack_err(info,"get_anti_regular_up_schurlev","dorm2r");

    F77CALL (dtrsv) ("Upper","No transpose","No unit",&k,e_up_save,&k,
	  e_lef,&inc1);

    F77CALL (dlapmt) (&flag,&inc1,&k,e_lef,&inc1,d_perm);
#else
    F77CALL (dgesv) (&k,&inc1,e_up_save,&k,d_perm,e_lef,&k,&info);
    lapack_err(info,"get_anti_regular_up_schurlev","dgesv");
#endif
  }

  theta = -theta;
  /* store vartheta in theta_[k], theta in theta_[0:k-1] */
  theta_[k] = theta;

  if (q != NULL)
  {
    mult_by_z_up_k(&q_temp,&q_up[k-1],&one,1);
      
    ihelp = n+k;

    F77CALL (dgemv) ("No transpose",&ihelp,&k,&one,q_up[0].data,ldq,
	  e_lef,&inc1,&one,q_temp.data,&q_temp.stdinc);

    polynomaxpy(&theta,q_lef,&q_temp);

    q_up[0].length = q_temp.length;
    F77CALL (dcopy) (&q_temp.length,q_temp.data,&q_temp.stdinc,
	   q_up[0].data,&q_up[0].stdinc);
  }

  /* e_lef[k] still contains the correct value of e_lef */
  ihelp = N-n-k;
  F77CALL (dcopy) (&ihelp,&e_up(k-1,k-1),&inc1,work,&inc1);
  F77CALL (daxpy) (&ihelp,&theta,&e_lef[k],&inc1,work,&inc1);
  F77CALL (dgemv) ("No transpose",&ihelp,&k,&one,&e_up(k,0),ldq,e_lef,&inc1,
	 &one,work,&inc1);
  F77CALL (dcopy) (&ihelp,work,&inc1,&e_up(k,0),&inc1);


  ihelp = N-n-k-1;
  F77CALL (dcopy) (&ihelp,&p_up(k,k-1),&inc1,work,&inc1);
  F77CALL (daxpy) (&ihelp,&theta,&p_lef[1],&inc1,work,&inc1);
  ldp1 = *ldq+1;
  F77CALL (dgemv) ("No transpose",&ihelp,&k,&one,&p_up(0,0),&ldp1,e_lef,&inc1,
	 &one,work,&inc1);
  F77CALL (dcopy) (&ihelp,work,&inc1,&p_up(k,0),&inc1);

  /* store theta in theta_[0:k-1] */
  F77CALL (dcopy) (&k,e_lef,&inc1,theta_,&inc1);

  return 1;
#undef e_up
#undef p_up
#undef d
}

int get_anti_regular_lef(fdouble *work,finteger lwork,polynom *q,
			 polynom *q_lef,fdouble *e__,fdouble *e_lef,
			 fdouble *p__,fdouble *d__,finteger *d_perm,
			 fdouble *d_tau,finteger *ldsmall,finteger *ldq,
			 finteger n,finteger k)
     /* requires 4*k workspace */
{
  finteger ihelp = k-1;
  fdouble *temp1=work;
#ifdef USE_QR_DECOMP
  fdouble *lawork=&work[k];
  flogical flag;
  finteger i;
  lwork -=k;
#ifdef WORK_DEBUG
  if (lwork < 3*k)
    fprintf(stderr,"get_anti_regular_lef has got too little workspace\n");
#endif
#else
#ifdef WORK_DEBUG
  lwork -=k;
  if (lwork < 0)
    fprintf(stderr,"get_anti_regular_lef has got too little workspace\n");
#endif
#endif

#define p(I,J)      p__[(I) + (J) * (*ldsmall)]
#define d(I,J)      d__[(I) + (J) * (*ldsmall)]
#define e(I,J)      e__[(I) + (J) * (*ldsmall)]

  /* calculating new q_lef */
  F77CALL (dcopy) (&ihelp,&zero,&inc0,temp1,&inc1);
  temp1[ihelp]=1.0;
  if (d__ !=  NULL)
  {
#ifdef USE_QR_DECOMP
    flag=FTRUE;
    F77CALL (dlapmt) (&flag,&inc1,&k,temp1,&inc1,d_perm);

    /* solve the equation */
    F77CALL (dtrsv) ("Upper","Transpose","No unit",&k,&d(0,0),ldsmall,
	  temp1,&inc1);

    F77CALL (dorm2r)("left","No transpose",&k,&inc1,&k,&d(0,0),ldsmall,d_tau,
	    temp1,ldsmall,lawork,&info);
    lapack_err(info,"get_anti_regular_lef","dorm2r");
#else
    F77CALL (dgetrs) ("Transpose",&k,&inc1,&d(0,0),ldsmall,d_perm,
	  temp1,ldsmall,&info);
    lapack_err(info,"get_anti_regular_lef","dgetrs");
#endif
  }
  else
  {
#ifdef USE_QR_DECOMP
    flag=FFALSE;
    /* get QR decomposition of D */
    for (i=0;i<k;d_perm[i++]=0);
    F77CALL (dgeqpf) (&k,&k,&p(0,0),ldsmall,d_perm,d_tau,lawork,&info);
    lapack_err(info,"get_anti_regular_lef","dgeqpf");

    /* solve the equation */
    F77CALL (dorm2r)("left","Transpose",&k,&inc1,&k,&p(0,0),ldsmall,d_tau,temp1,ldsmall,
	    lawork,&info);
    lapack_err(info,"get_anti_regular_lef","dorm2r");

    F77CALL (dtrsv) ("Upper","No transpose","No unit",&k,&p(0,0),ldsmall,
	  temp1,&inc1);

    F77CALL (dlapmt) (&flag,&inc1,&k,temp1,&inc1,d_perm);
#else
    F77CALL (dgesv) (&k,&inc1,&p(0,0),ldsmall,d_perm,temp1,&k,&info);
    lapack_err(info,"get_anti_regular_lef","dgesv");
#endif
  }

  q_lef->length = n+k;
  F77CALL (dgemv) ("No transpose",&q_lef->length,&k,&one,q[0].data,ldq,
	temp1,&inc1,&zero,q_lef->data,&inc1);

  ihelp=*ldsmall+1;
  (*e_lef) = F77CALL (ddot) (&k,temp1,&inc1,&e(0,0),&ihelp);

  return 1;
#undef p
#undef e
#undef d
}

int get_anti_regular_lef_schurlev(fdouble *work,finteger lwork,fdouble *theta,
      polynom *q,polynom *q_lef,fdouble *e__,fdouble *e_lef,fdouble *p__,
      fdouble *p_lef,fdouble *d__,finteger *d_perm,fdouble *d_tau,
      finteger *ldsmall,finteger *ldq,finteger N,finteger n,finteger k)
     /* requires k*k + 3*k workspace */
{
  finteger ihelp = k-1;
  finteger ldp1;
  fdouble *p_save=work;
#ifdef USE_QR_DECOMP
  finteger i;
  flogical flag;
  fdouble *lawork=work+k*k;
  lwork-=k*k;
#ifdef WORK_DEBUG
  if (lwork < 3*k)
    fprintf(stderr,"get_anti_regular_lef_schurlev has got too little "
	    "workspace\n");
#endif
#else
  lwork-=k*k;
#ifdef WORK_DEBUG
  if (lwork < 0)
    fprintf(stderr,"get_anti_regular_lef_schurlev has got too little "
	    "workspace\n");
#endif
#endif

#define p(I,J)      p__[(I) + (J) * (*ldq)]
#define d(I,J)      d__[(I) + (J) * (*ldsmall)]
#define e(I,J)      e__[(I) + (J) * (*ldq)]

  /* calculating new q_lef */
  F77CALL (dcopy) (&ihelp,&zero,&inc0,theta,&inc1);
  theta[ihelp]=1.0;

  if (d__ != NULL)
  {
#ifdef USE_QR_DECOMP
    /* solve the equation */
    flag=FTRUE;
    F77CALL (dlapmt) (&flag,&inc1,&k,theta,&inc1,d_perm);

    F77CALL (dtrsv) ("Upper","Transpose","No unit",&k,&d(0,0),ldsmall,
	  theta,&inc1);

    F77CALL (dorm2r)("left","No transpose",&k,&inc1,&k,&d(0,0),ldsmall,d_tau,
	    theta,ldsmall,lawork,&info);
    lapack_err(info,"get_anti_regular_lef_schurlev","dorm2r");
#else
    F77CALL (dgetrs) ("Transpose",&k,&inc1,&d(0,0),ldsmall,d_perm,
	  theta,ldsmall,&info);
    lapack_err(info,"get_anti_regular_lef_schurlev","dgetrs");
#endif
  }
  else
  {
    F77CALL (dlacpy) ("A",&k,&k,&p(0,0),ldq,p_save,&k);

#ifdef USE_QR_DECOMP
    /* get QR decomposition of D */
    for (i=0;i<k;d_perm[i++]=0);
    F77CALL (dgeqpf) (&k,&k,p_save,&k,d_perm,d_tau,lawork,&info);
    lapack_err(info,"get_anti_regular_lef_schurlev","dgeqpf");

    /* solve the equation */
    F77CALL (dorm2r)("left","Transpose",&k,&inc1,&k,p_save,&k,d_tau,
	    theta,ldsmall,lawork,&info);
    lapack_err(info,"get_anti_regular_lef_schurlev","dorm2r");

    F77CALL (dtrsv) ("Upper","No transpose","No unit",&k,p_save,&k,theta,&inc1);

    flag=FFALSE;
    F77CALL (dlapmt) (&flag,&inc1,&k,theta,&inc1,d_perm);
#else
    F77CALL (dgesv) (&k,&inc1,p_save,&k,d_perm,theta,&k,&info);
    lapack_err(info,"get_anti_regular_lef_schurlev","dgesv");
#endif
  }

  if (q != NULL)
  {
    q_lef->length = n+k;
    F77CALL (dgemv) ("No transpose",&q_lef->length,&k,&one,q[0].data,ldq,
	  theta,&inc1,&zero,q_lef->data,&inc1);
  }

  ihelp = N-n-k;
  p_lef[k] = 1.0;
  F77CALL (dgemv) ("No transpose",&ihelp,&k,&one,&p(k,0),ldq,theta,&inc1,
	 &zero,&p_lef[k+1],&inc1);

  ihelp = N-n-1;
  ldp1 = *ldq+1;
  F77CALL (dgemv) ("No transpose",&ihelp,&k,&one,&e(0,0),&ldp1,theta,&inc1,
	 &zero,&e_lef[k],&inc1);

  return 1;
#undef p
#undef e
#undef d
}

void get_regular_q(polynom *q,polynom *q_up,polynom *q_lef,fdouble *e_lef,
	       fdouble *e_up)
{
  fdouble help = -e_lef[0];

  /* calculating new q */
  F77CALL (dcopy) (&q_up[0].length,q_up[0].data,&q_up[0].stdinc,
	q[0].data,&q[0].stdinc);
  q[0].length=q_up[0].length;

  F77CALL (dscal) (&q[0].length,&help,q[0].data,&inc1);

  polynomaxpy(e_up,q_lef,q);
}

void get_regular_q_schurlev(polynom *q,polynom *q_up,polynom *q_lef,
		   fdouble *e,fdouble *e_lef,fdouble *e_up,
		   fdouble *p,fdouble *p_lef,fdouble *p_up,
		   finteger N,finteger n,finteger k)

{
  fdouble help = -e_lef[0];
  finteger ihelp;

  if (q != NULL)
  {
    /* calculating new q */
    F77CALL (dcopy) (&q_up[0].length,q_up[0].data,&q_up[0].stdinc,q[0].data,
	   &q[0].stdinc);
    q[0].length=q_up[0].length;
    
    F77CALL (dscal) (&q[0].length,&help,q[0].data,&inc1);

    polynomaxpy(e_up,q_lef,q);
  }

  /* calculating new p,e */
  ihelp = N-n-k-1;
  p[0] = 0.0;
  F77CALL (dcopy) (&ihelp,p_up,&inc1,&p[1],&inc1);
  F77CALL (dscal) (&ihelp,&help,&p[1],&inc1);
  ihelp++;
  F77CALL (daxpy) (&ihelp,e_up,p_lef,&inc1,p,&inc1);

  ihelp = N-n-k;
  F77CALL (dcopy) (&ihelp,&e_up[1],&inc1,e,&inc1);
  F77CALL (dscal) (&ihelp,&help,e,&inc1);
  F77CALL (daxpy) (&ihelp,e_up,&e_lef[1],&inc1,e,&inc1);
}

fdouble get_e_lef_k(polynom* q_lef, polynom* q_up, fdouble* e_lef,
		       fdouble* e_up__, finteger *ldsmall, finteger k)

     /* compute epsilon_lef_{n+k,n} */
{
#define e_up(I,J)   e_up__[(I) + (J) * (*ldsmall)]

  fdouble res;
  finteger clip1,clip2,length;

  if (0 < (length=length_of_data(q_lef,0,&clip1,k,&clip2)))
    res=F77CALL (ddot) (&length,&e_up(clip2,0),&inc1,data_pointer(q_lef,0,k),&inc1);
  else
    res=0.0;

  if (0 < (length=length_of_data(q_up,1,&clip1,k,&clip2)))
    res-=F77CALL (ddot) (&length,&e_lef[clip2],&inc1,data_pointer(q_up,1,k),&incm1);

  res /= poly_coef(q_up,0);

  return res;

#undef e_up
}

fdouble get_p_lef_km1(polynom* q_lef, polynom* q_up, fdouble* p_lef,
			 fdouble* p_up__, finteger *ldsmall, finteger k)
     /* compute pi_lef_{k-1,n} */
{
#define p_up(I,J)   p_up__[(I) + (J) * (*ldsmall)]

  fdouble res;
  finteger clip1,clip2,length;

  finteger n = q_lef->length;

  if (0 < (length=length_of_data(q_lef,n-k+2,&clip1,n-1,&clip2)))
    res=F77CALL (ddot) (&length,&p_up(clip1,0),&inc1,data_pointer(q_lef,n-k+2,n-1),
	      &incm1);
  else
    res=0.0;

  if (0 < (length=length_of_data(q_up,n-k+1,&clip1,n-1,&clip2)))
    res-=F77CALL (ddot) (&length,&p_lef[clip1],&inc1,data_pointer(q_up,n-k+1,n-1),&inc1);

  return(res);

#undef p_up
}

int get_anti_regular_pair(fdouble *work,finteger lwork,fdouble *M,
      finteger *ldM,polynom *q_lef,polynom *q_up,finteger *ldq,fdouble *e_up,
      fdouble *e_lef,fdouble *p_up,fdouble *p_lef,finteger n,finteger k)
     /* requires 2*ldq + 12*k - 3 + liwork(2*k-1) workspace */
{
  finteger i,ihelp=2*k-1;
  fdouble help;
  fdouble *new_q_up=work;
  fdouble *new_q_lef=&new_q_up[*ldq];
  fdouble *temp=&new_q_lef[*ldq];
  fdouble *alpha=&temp[2*k-1];
  fdouble *beta=&alpha[k];
#ifndef USE_QR_DECOMP
  finteger *perm=(finteger *)&beta[k+1];
#ifdef WORK_DEBUG
  lwork -= 2*(*ldq) + 2*k-1 + 2*k+1 + liwork(2*k-1);
  if (lwork < 0)
    fprintf(stderr,"get_anti_regular_pair has got too little workspace\n");
#endif
#else
  flogical flag=FFALSE;
  fdouble *tau=&beta[k+1];
  fdouble *lawork=&tau[2*k-1];
  finteger *perm;
  lwork -= 2*(*ldq) + 2*k-1 + 2*k+1 + 2*k-1 + liwork(2*k-1);
  perm=(finteger*)&lawork[lwork];
#ifdef WORK_DEBUG
  if (lwork < 3*k)
    fprintf(stderr,"get_anti_regular_pair has got too little workspace\n");
#endif
#endif

  /* calculating anti-regular q_up */
  beta[k] = 1.0;
  alpha[k-1] = -p_up[0]/p_lef[0];

  /* compute right-hand side */
  ihelp = k-1;
  F77CALL (dcopy) (&ihelp,&p_up[1],&incm1,temp,&inc1);
  F77CALL (dscal) (&ihelp,&mone,temp,&inc1);  
  F77CALL (dcopy) (&k,&zero,&inc0,&temp[k-1],&inc1);
  help = -alpha[k-1];
  F77CALL (daxpy) (&ihelp,&help,&p_lef[1],&incm1,temp,&inc1);
  temp[2*k-2] -= e_lef[0]*alpha[k-1];

  /* solve linear system to compute recurrence coefficients */
  /* save solution in alpha */

  ihelp = 2*k-1;

#ifdef USE_QR_DECOMP
  /* get QR decomposition of M */
  for (i=0;i<ihelp;perm[i++]=0);   
  F77CALL (dgeqpf) (&ihelp,&ihelp,M,ldM,perm,tau,lawork,&info);
  lapack_err(info,"get_anti_regular_pair","dgeqpf");
  
  /* solve the equation */
  F77CALL (dorm2r)("left","Transpose",&ihelp,&inc1,&ihelp,M,ldM,tau,
	  temp,ldM,lawork,&info);
  lapack_err(info,"get_anti_regular_pair","dorm2r");
  
  F77CALL (dtrsv) ("Upper","No transpose","No unit",&ihelp,M,ldM,temp,&inc1);

  F77CALL (dlapmt) (&flag,&inc1,&ihelp,temp,&inc1,perm);
#else
  F77CALL (dgesv) (&ihelp,&inc1,M,ldM,perm,temp,&ihelp,&info);
  lapack_err(info,"get_anti_regular_pair","dgesv");
#endif

  ihelp = k-1;
  F77CALL (dcopy) (&ihelp,temp,&inc1,alpha,&inc1);
  F77CALL (dcopy) (&k,&temp[k-1],&inc1,beta,&inc1);

  F77CALL (dcopy) (&k,&zero,&inc0,new_q_up,&inc1);
  F77CALL (dcopy) (&q_up->length,q_up->data,&q_up->stdinc,&new_q_up[k],&inc1);

  for (i=0;i<k;i++)
  {
    F77CALL (daxpy) (&q_lef->length,&alpha[i],q_lef->data,&q_lef->stdinc,
	   &new_q_up[i],&inc1);
    F77CALL (daxpy) (&q_up->length,&beta[i],q_up->data,&q_up->stdinc,
	   &new_q_up[i],&inc1);
  }
    
  ihelp = k+1;
  e_up[0] = F77CALL (ddot) (&ihelp,e_up,&incm1,beta,&inc1);
  e_up[0] += F77CALL (ddot) (&k,&e_lef[1],&incm1,alpha,&inc1);

  /* calculating anti-regular q_lef */

  /* compute right-hand side */
  ihelp = 2*k-1;
  F77CALL (dcopy) (&ihelp,&zero,&inc0,temp,&inc1);
  temp[0] = 1.0;

  /* solve linear system to compute recurrence coefficients */
  /* save solution in alpha */

#ifdef USE_QR_DECOMP
  /* solve the equation */
  F77CALL (dorm2r)("left","Transpose",&ihelp,&inc1,&ihelp,M,ldM,tau,
	  temp,&ihelp,lawork,&info);
  lapack_err(info,"get_anti_regular_pair","dorm2r");
  
  F77CALL (dtrsv) ("Upper","No transpose","No unit",&ihelp,M,ldM,temp,&inc1);

  F77CALL (dlapmt) (&flag,&inc1,&ihelp,temp,&inc1,perm);
#else
  F77CALL (dgetrs) ("No transpose",&ihelp,&inc1,M,ldM,perm,temp,&ihelp,&info);
  lapack_err(info,"get_anti_regular_pair","dgetrs");
#endif

  ihelp = k-1;
  F77CALL (dcopy) (&ihelp,temp,&inc1,alpha,&inc1);
  F77CALL (dcopy) (&k,&temp[k-1],&inc1,beta,&inc1);

  F77CALL (dcopy) (&ihelp,&zero,&inc0,new_q_lef,&inc1);
  F77CALL (dcopy) (&q_up->length,q_up->data,&q_up->stdinc,
	&new_q_lef[k-1],&inc1);
  F77CALL (dscal) (&q_up->length,&beta[k-1],&new_q_lef[k-1],&inc1);  

  for (i=0;i<k-1;i++)
  {
    F77CALL (daxpy) (&q_lef->length,&alpha[i],q_lef->data,&q_lef->stdinc,
	   &new_q_lef[i],&inc1);
    F77CALL (daxpy) (&q_up->length,&beta[i],q_up->data,&q_up->stdinc,
	   &new_q_lef[i],&inc1);
  }

  ihelp = k-1;
  e_lef[0] = F77CALL (ddot) (&k,&e_up[1],&incm1,beta,&inc1);
  e_lef[0] += F77CALL (ddot) (&ihelp,&e_lef[2],&incm1,alpha,&inc1);
    
  q_up[0].length = n+k+1;
  q_lef->length   = n+k;
  F77CALL (dcopy) (&q_up->length,new_q_up,&inc1,q_up->data,&q_up->stdinc);
  F77CALL (dcopy) (&q_lef->length,new_q_lef,&inc1,q_lef->data,&q_lef->stdinc);

  return 1;
}

int get_anti_regular_pair_schurlev(fdouble *work,finteger lwork,
      fdouble *M,finteger *ldM,polynom *q_lef,polynom *q_up,finteger *ldq,
      fdouble *e_up,fdouble *e_lef,fdouble *p_up,fdouble *p_lef,
      fdouble *schur_par,finteger N,finteger n,finteger k)
     /* requires 6*ldq + 12*k - 3 + liwork(2*k-1) workspace */
{
  finteger ihelp,ihelp2,i;
  fdouble help;
  fdouble *new_q_up=work;
  fdouble *new_q_lef=&new_q_up[*ldq];
  fdouble *new_p_up=&new_q_lef[*ldq];
  fdouble *new_p_lef=&new_p_up[*ldq];
  fdouble *new_e_up=&new_p_lef[*ldq];
  fdouble *new_e_lef=&new_e_up[*ldq];
  fdouble *alpha=&new_e_lef[*ldq];
  fdouble *beta=&alpha[k];
  fdouble *temp=&beta[k+1];
#ifndef USE_QR_DECOMP  
  finteger *perm=(finteger*)&temp[2*k-1];
  lwork -= 6*(*ldq) + 2*k+1 + 2*k-1 + liwork(2*k-1);
#ifdef WORK_DEBUG
  if (lwork < 0)
    fprintf(stderr,"get_anti_regular_pair_schurlev has got too little "
	    "workspace\n");
#endif
#else
  finteger *perm;
  flogical flag=FFALSE;
  fdouble *tau=&temp[2*k-1];
  fdouble *lawork=&tau[2*k-1];
  lwork -= 6*(*ldq) + 2*k+1 + 2*k-1 + 2*k-1 + liwork(2*k-1);
  perm = (finteger*)&lawork[lwork];
#ifdef WORK_DEBUG
  if (lwork < 3*k)
    fprintf(stderr,"get_anti_regular_pair_schurlev has got too little "
	    "workspace\n");
#endif
#endif

  ihelp = 2*k-1;

  /* calculating anti-regular q_up */
  beta[k] = 1.0;
  alpha[k-1] = -p_up[0]/p_lef[0];

  /* compute right-hand side */
  ihelp = k-1;
  F77CALL (dcopy) (&ihelp,&p_up[1],&incm1,temp,&inc1);
  F77CALL (dscal) (&ihelp,&mone,temp,&inc1);  
  F77CALL (dcopy) (&k,&zero,&inc0,&temp[k-1],&inc1);
  help = -alpha[k-1];
  F77CALL (daxpy) (&ihelp,&help,&p_lef[1],&incm1,temp,&inc1);
  temp[2*k-2] -= e_lef[0]*alpha[k-1];

  /* solve linear system to compute recurrence coefficients */
  /* save solution in alpha */

  ihelp = 2*k-1;

#ifdef USE_QR_DECOMP
  /* get QR decomposition of M */
  for (i=0;i<ihelp;perm[i++]=0);
  F77CALL (dgeqpf) (&ihelp,&ihelp,M,ldM,perm,tau,lawork,&info);
  lapack_err(info,"get_anti_regular_pair_schurlev","dgeqpf");
  
  /* solve the equation */
  F77CALL (dorm2r)("left","Transpose",&ihelp,&inc1,&ihelp,M,ldM,tau,
	  temp,ldM,lawork,&info);
  lapack_err(info,"get_anti_regular_pair_schurlev","dorm2r");
  
  F77CALL (dtrsv) ("Upper","No transpose","No unit",&ihelp,M,ldM,temp,&inc1);

  F77CALL (dlapmt) (&flag,&inc1,&ihelp,temp,&inc1,perm);
#else
  F77CALL (dgesv) (&ihelp,&inc1,M,ldM,perm,temp,&ihelp,&info);
  lapack_err(info,"get_anti_regular_pair_schurlev","dgesv");
#endif

  ihelp = k-1;
  F77CALL (dcopy) (&ihelp,temp,&inc1,alpha,&inc1);
  F77CALL (dcopy) (&k,&temp[k-1],&inc1,beta,&inc1);

  if (schur_par != NULL)
  {
    F77CALL (dcopy) (&k,alpha,&inc1,&ALPHA_UP(n),&inc1);
    F77CALL (dcopy) (&k,beta,&inc1,&BETA_UP(n),&inc1);
  }

  if (q_up != NULL)
  {
    F77CALL (dcopy) (&k,&zero,&inc0,new_q_up,&inc1);
    F77CALL (dcopy) (&q_up->length,q_up->data,&q_up->stdinc,
	  &new_q_up[k],&inc1);
  }

  ihelp = N-n-k;
  F77CALL (dcopy) (&ihelp,&e_up[0],&inc1,new_e_up,&inc1);

  ihelp2 = N-n-k-1;
  F77CALL (dcopy) (&ihelp2,&p_up[k],&inc1,new_p_up,&inc1);

  for (i=0;i<k;i++)
  {
    if (q_up != NULL)
    {
      F77CALL (daxpy) (&q_lef->length,&alpha[i],q_lef->data,&q_lef->stdinc,
	     &new_q_up[i],&inc1);
      F77CALL (daxpy) (&q_up->length,&beta[i],q_up->data,&q_up->stdinc,
	     &new_q_up[i],&inc1);
    }

    F77CALL (daxpy) (&ihelp,&alpha[i],&e_lef[k-i],&inc1,
	   new_e_up,&inc1);
    F77CALL (daxpy) (&ihelp,&beta[i],&e_up[k-i],&inc1,
	   new_e_up,&inc1);

    F77CALL (daxpy) (&ihelp2,&alpha[i],&p_lef[i+1],&inc1,
	   new_p_up,&inc1);
    F77CALL (daxpy) (&ihelp2,&beta[i],&p_up[i],&inc1,
	   new_p_up,&inc1);
  }


  /* calculating anti-regular q_lef */

  /* compute right-hand side */
  ihelp = 2*k-1;
  F77CALL (dcopy) (&ihelp,&zero,&inc0,temp,&inc1);
  temp[0] = 1.0;

  /* solve linear system to compute recurrence coefficients */
  /* save solution in alpha */

#ifdef USE_QR_DECOMP
  /* solve the equation */
  F77CALL (dorm2r)("left","Transpose",&ihelp,&inc1,&ihelp,M,ldM,tau,
	  temp,&ihelp,lawork,&info);
  lapack_err(info,"get_anti_regular_pair_schurlev","dorm2r");
  
  F77CALL (dtrsv) ("Upper","No transpose","No unit",&ihelp,M,ldM,temp,&inc1);

  F77CALL (dlapmt) (&flag,&inc1,&ihelp,temp,&inc1,perm);
#else
  F77CALL (dgetrs) ("No transpose",&ihelp,&inc1,M,ldM,perm,temp,&ihelp,&info);
  lapack_err(info,"get_anti_regular_pair_schurlev","dgetrs");
#endif

  ihelp = k-1;
  F77CALL (dcopy) (&ihelp,temp,&inc1,alpha,&inc1);
  F77CALL (dcopy) (&k,&temp[k-1],&inc1,beta,&inc1);

  if (schur_par != NULL)
  {
    F77CALL (dcopy) (&ihelp,alpha,&inc1,&ALPHA_LEF(n),&inc1);
    F77CALL (dcopy) (&k,beta,&inc1,&BETA_LEF(n),&inc1);
  }

  if (q_up != NULL)
  {
    F77CALL (dcopy) (&ihelp,&zero,&inc0,new_q_lef,&inc1);
    F77CALL (dcopy) (&q_up->length,q_up->data,&q_up->stdinc,
	  &new_q_lef[k-1],&inc1);
    F77CALL (dscal) (&q_up->length,&beta[k-1],&new_q_lef[k-1],&inc1);  
  }

  ihelp = N-n-k;
  F77CALL (dcopy) (&ihelp,&e_up[1],&inc1,new_e_lef,&inc1);
  F77CALL (dscal) (&ihelp,&beta[k-1],new_e_lef,&inc1);

  ihelp2 = N-n-k-1;
  F77CALL (dcopy) (&ihelp2,&p_up[k-1],&inc1,new_p_lef,&inc1);
  F77CALL (dscal) (&ihelp2,&beta[k-1],new_p_lef,&inc1);

  for (i=0;i<k-1;i++)
  {
    if (q_up != NULL)
    {
      F77CALL (daxpy) (&q_lef->length,&alpha[i],q_lef->data,&q_lef->stdinc,
	     &new_q_lef[i],&inc1);
      F77CALL (daxpy) (&q_up->length,&beta[i],q_up->data,&q_up->stdinc,
	     &new_q_lef[i],&inc1);
    }

    F77CALL (daxpy) (&ihelp,&alpha[i],&e_lef[k-i],&inc1,
	   new_e_lef,&inc1);
    F77CALL (daxpy) (&ihelp,&beta[i],&e_up[k-i],&inc1,
	   new_e_lef,&inc1);

    F77CALL (daxpy) (&ihelp2,&alpha[i],&p_lef[i+1],&inc1,
	   new_p_lef,&inc1);
    F77CALL (daxpy) (&ihelp2,&beta[i],&p_up[i],&inc1,
	   new_p_lef,&inc1);
  }

  if (q_up != NULL)
  {
    q_up[0].length = n+k+1;
    q_lef->length   = n+k;
    F77CALL (dcopy) (&q_up->length,new_q_up,&inc1,q_up->data,&q_up->stdinc);
    F77CALL (dcopy) (&q_lef->length,new_q_lef,&inc1,
	  q_lef->data,&q_lef->stdinc);
  }

  ihelp = N-n-k;
  F77CALL (dcopy) (&ihelp,new_e_up,&inc1,&e_up[k],&inc1);
  F77CALL (dcopy) (&ihelp,new_e_lef,&inc1,&e_lef[k],&inc1);

  ihelp = N-n-k-1;
  F77CALL (dcopy) (&ihelp,new_p_up,&inc1,&p_up[k],&inc1);
  F77CALL (dcopy) (&ihelp,new_p_lef,&inc1,&p_lef[k+1],&inc1);
  p_lef[k] = 1.0;

  return 1;
}

void get_alpha_beta_up_from_theta(fdouble *schur_par,fdouble *theta,
      finteger N,finteger n,finteger k)
{
  finteger ihelp = k-1;
/*   flaprintf(stdout,"theta=%ladvte\n",k+1,theta,1); */

  /* Compute coefficients of u_up, v_up of product form algorithm */
  ALPHA_UP(n) = theta[k];
  F77CALL (dcopy) (&ihelp,&GAMMA_UP(n),&incm1,&ALPHA_UP(n+1),&inc1);
/*   flaprintf(stdout,"alpha_up=%ladvte\n",k,&ALPHA_UP(n),1); */

  for (ihelp=k-1;ihelp>0;ihelp--)
    ALPHA_UP(n+k-1-ihelp) +=
      F77CALL (ddot) (&ihelp,&GAMMA_UP(n),&inc1,theta+k-ihelp,&inc1);
/*   flaprintf(stdout,"alpha_up=%ladvte\n",k,&ALPHA_UP(n),1); */
  
  F77CALL (dcopy) (&k,theta,&inc1,&BETA_UP(n),&inc1);
/*   flaprintf(stdout,"beta_up=%ladvte\n",k,&BETA_UP(n),1); */
}

void get_alpha_beta_lef_from_theta(fdouble *schur_par,fdouble *theta,
      fdouble *v,finteger N,finteger n,finteger k)
  /* Compute coefficients of u_lef, v_lef of product form algorithm */
{
  finteger ihelp=k-1;

  F77CALL (dcopy) (&ihelp,&zero,&inc0,&ALPHA_LEF(n),&inc1);

  for (ihelp=1;ihelp<k;ihelp++)
    ALPHA_LEF(n+ihelp-1)=
      F77CALL (ddot) (&ihelp,v,&inc1,theta+k-ihelp,&inc1);
 
  F77CALL (dcopy) (&k,theta,&incm1,&BETA_LEF(n),&inc1);
  if (n>0)
    F77CALL (dscal) (&k,&BETA(n),&BETA_LEF(n),&inc1);
}

int decomp_matrix(fdouble *lawork,finteger *k,fdouble *mat,
		  finteger *ldmat,finteger *perm,fdouble *tau)
{
#ifdef USE_QR_DECOMP
  finteger i;

  for (i=0;i< *k;perm[i++]=0);
  F77CALL (dgeqpf) (k,k,mat,ldmat,perm,tau,lawork,&info);
  lapack_err(info,"decomp_matrix","dgeqpf");
#else
  F77CALL (dgetrf) (k,k,mat,ldmat,perm,&info);
  lapack_err(info,"decomp_matrix","dgetrf");  
#endif
  return 1;
}
