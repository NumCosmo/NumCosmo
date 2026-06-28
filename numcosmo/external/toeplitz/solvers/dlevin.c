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
#include "drefine.h"
#include "dopt_k.h"

#define SING_TOL   1e-12

static finteger inc1  =  1;
static finteger inc0  =  0;

#ifndef CALCULATE_WITH_D
static fdouble one  =  1.0;
static fdouble zero =  0.0;
#endif

int d_lev_inner(fdouble *x,fdouble *mu,finteger N,fdouble *b,
		fdouble col_tol,fdouble tol,int k_max,int refine, fdouble *psol)
     /*
       x     pointer to the array for the solution.
       mu    pointer to the diagonal element of the toeplitz matrix.
       N     the dimension of the toeplitz matrix mu.
       b     pointer to the righthandside. */
{
  /* statistical information */
#if (DEBUG >= 1)
  finteger number_of_steps[100];
#endif

  void *data1,*data2;

  /* finteger help variables */
  finteger n                   = 0;
  finteger k                   = 0;
  finteger success             = 1;
  finteger regular             = 0;
  finteger column_regular      = 0;
  finteger singular            = 0;
  finteger ihelp,i;

  /* variable which indicates in which polynom the current classical 
     polynomials are */
  finteger storage_index=0;

  /* tolerance for singular matrix */
  /*   const fdouble col_tol=COL_TOL; */
  /*   const fdouble tol=TOL; */
  const fdouble sing_tol=SING_TOL;

  /* help variables */
#ifndef CALCULATE_WITH_D
  fdouble help2;
#endif
  fdouble help;
  fdouble *work;
  finteger lwork;

  /* variables for the method */
  polynom *q,*q_up,q_lef;
  finteger ldq=N+1;
  toeplitz h;
  fdouble *p_up__,*p__;
  fdouble *e_up__,*e__;
  fdouble *d__;
  fdouble *p_lef,*e_lef;
  fdouble *v,*v_up;
  finteger ldsmall=k_max+1;
  finteger *d_perm;
  fdouble *d_tau;

#define p_up(I,J)   p_up__[(I) + (J) * ldsmall]
#define p(I,J)      p__[(I) + (J) * ldsmall]
#define e_up(I,J)   e_up__[(I) + (J) * ldsmall]
#define e(I,J)      e__[(I) + (J) * ldsmall]
#define d(I,J)      d__[(I) + (J) * ldsmall]

  /* initialize room for statistical information */
#if (DEBUG >= 1)
  for (i=0;i<k_max;i++)
    number_of_steps[i]=0;
#endif
  /* memory management */

  q=(polynom*)(data1=malloc(2*max(k_max,2)*sizeof(polynom)));
  if (data1==NULL)
  {
    fprintf(stderr,"error can't allocate memory for working!\nAbort ...\n");
    return 0;
  }
  q_up=&q[max(2,k_max)];

  /* initialize q,q_up,h */
  h.data=mu;
  h.N=N;

  for(i=0;i<max(k_max,2);i++)
  {
    q[i].length=0;
    q[i].stdinc=-1;

    q_up[i].length=0;
    q_up[i].stdinc=1;
  }
  q_lef.length=0;
  q_lef.stdinc=-1;

  /* space for q and q_up */
  ihelp=sizeof(fdouble)*ldq*2*max(k_max,2);
  /* space for q_lef */
  ihelp+=ldq*sizeof(fdouble);
  /* space for d,e,e_up,p,p_up */
  ihelp+=5*ldsmall*ldsmall*sizeof(fdouble);
  /* space for e_lef,p_lef,v,v_up */
  ihelp+=4*ldsmall*sizeof(fdouble);
  /* space for d_perm or d_tau */
#ifdef USE_QR_DECOMP
  ihelp+=ldsmall*(sizeof(fdouble)+sizeof(finteger));
#else
  ihelp+=ldsmall*sizeof(finteger);
#endif
  /* space for temporary variables */
  lwork=max(k_max*k_max+6*k_max,k_max+ldq);
#ifdef USE_QR_DECOMP
  lwork+=k_max*3;
#endif
  ihelp+=lwork*sizeof(fdouble);
  q[0].data=(fdouble *)(data2=malloc(ihelp));
  if (data2==NULL)
  {
    fprintf(stderr,"error can't allocate memory for working!\nAbort ...\n");
    free(data1);
    return(0);
  }
  /* distribute memory */
  q_up[0].data=&(q[0].data[max(k_max,2)*ldq]);
  for(i=1;i<max(k_max,2);i++)
  {
    q[i].data=&q[0].data[ldq*i];
    q_up[i].data=&q[0].data[ldq*(i+max(k_max,2))];
  }
  q_lef.data=&q_up[max(k_max,2)-1].data[ldq];
  e__=&q_lef.data[ldq];
  e_up__=&e__[ldsmall*ldsmall];
  p__=&e_up__[ldsmall*ldsmall];
  p_up__=&p__[ldsmall*ldsmall];
  d__=&p_up__[ldsmall*ldsmall];
  p_lef=&d__[ldsmall*ldsmall];
  e_lef=&p_lef[ldsmall];
  v=&e_lef[ldsmall];
  v_up=&v[ldsmall];
  work=&v_up[ldsmall];
#ifdef USE_QR_DECOMP
  d_tau=&work[lwork];
  d_perm=(finteger*)&d_tau[ldsmall];
#else
  d_tau=NULL;
  d_perm=(finteger*)&work[lwork];
#endif

  /* initialize memory */
  help=0.0;
  F77CALL (dcopy) (&N,&help,&inc0,x,&inc1);
  ihelp=2*max(2,k_max)*ldq;
  F77CALL (dcopy) (&ihelp,&help,&inc0,q->data,&inc1);


#if (DEBUG >= 1)
  fprintf(stdout,"******************************************************\n"
	         "*       New look-ahead Levinson algorithm            *\n"
	         "******************************************************\n");
#endif


  /*************************/
  /* initialize iteration */
  /*************************/

  /* LDU decompostion of a scalar */
  e_up(0,0)=p(0,0)=h.data[0];
  q[storage_index].data[0]=1.0;
  q[storage_index].length++;
  q_up[storage_index].data[0]=1.0;
  q_up[storage_index].length++;

  if ((col_tol < fabs(e_up(0,0))) || (N==1))
  {
    /* update the solution */
    x[0]=b[0]/p(0,0);

    /* setting the flag */
    column_regular=1;

#if (DEBUG >= 2)
    fprintf(stdout,"%2sn = %d : Making column regular step with p = %e\n","",
	    (int)n,p(0,0));
#endif
#if (DEBUG >= 1)
    (number_of_steps[0])++;
#endif
  }
  else
  {
    /* preparing q_lef for the first look-ahead step */
    q_lef.data[0]=0.0;
    q_lef.length=0;

    /* preparing p_lef and e_lef for the first look-ahead step */
    p_lef[0]=1.0;
    e_lef[0]=-1.0;

    /* preventing NaN's */
    v[0]=0.0;
    v_up[0]=0.0;

    /* setting the flag */
    column_regular=0;

    /* storing the best singular value up to now */
    set_new_singular_value(1,p(0,0));

#if (DEBUG >= 1)
    fprintf(stdout,"n = %d : Step is not column regular with p = %e\n",
	    (int)n,p(0,0));
#endif
  }

  /* calculating common parameters for the column regular and 
     not column regular first step */
  e(0,0)=h.data[1];
  p_up(0,0)=h.data[-1];

  psol[0] = x[0];
  /*************************/
  /*     the main loop     */
  /*************************/
  while(n<N-1)
  {
    if (column_regular) /* n column regular */
    {
      p_lef[0]=p(0,0);
      e_lef[0]=e(0,0);

      ihelp=(storage_index+1)%2;
      get_colreg_pair(&q[ihelp],&q_up[ihelp],&h,
		      &q[storage_index],&q_up[storage_index],
		      &e(0,0),&e_up(0,0),&p(0,0),&p_up(0,0),v,v_up,n);
      storage_index=ihelp;

      n++;

      if (fabs(p(0,0)) > col_tol)  /* n column regular */
      {
	e(0,0)=lau_coef(&h,&q[storage_index],n+1);

	if (!update_x(work,lwork,&q[storage_index],&q_up[storage_index],&ldq,
		      b,x,&p(0,0),&ldsmall,NULL,NULL,n,1))  break;

        if (refine && (n+1 == N))
          set_refinement(&h,&q[storage_index],&q_up[storage_index],NULL,
			 &p(0,0));
#if (DEBUG >= 2)
        fprintf(stdout,"%2sn = %d : Making column regular step with p = %e\n",
		"",(int)n,p(0,0));
#endif
#if (DEBUG >= 1)
	(number_of_steps[0])++;
#endif
      }
      else
      {
#if (DEBUG >= 1)
        fprintf(stdout,"n = %d : Step is not column regular with p = %e\n",
		(int)n,p(0,0));
#endif
	/* setting the flag */
	column_regular=0;

        /* test if the end of the matrix is reached */
	if (n+1 == N)
	{
	  if (fabs(p(0,0)) < sing_tol)
	  {
	    fprintf(stderr,"Matrix is close to singular. No solution "
		    "computed\n");
            singular=1;
	  }
	  else
	  {
#if (DEBUG >= 1)
	    (number_of_steps[0])++;
	    fprintf(stderr,"Matrix may be ill conditioned, updating x, last"
		    " singular value is %e\n",fabs(p(0,0)));
#endif

	    if (!update_x(work,lwork,&q[storage_index],&q_up[storage_index],
			  &ldq,b,x,&p(0,0),&ldsmall,NULL,NULL,n,1))  break;

	    if (refine)
	      set_refinement(&h,&q[storage_index],&q_up[storage_index],NULL,
			     &p(0,0));
	  }

	  break;
	}

        /* storing the best singular value up to now */
	set_new_singular_value(1,p(0,0));

        /* preparing q_lef, q_up and q for the first look-ahead step */
	if (storage_index == 1)
	{
	  /* copying last q into q_lef */
	  q_lef.length=q[0].length;
	  F77CALL (dcopy) (&q[0].length,q[0].data,&inc1,q_lef.data,&inc1);

	  /* copying q_up and q if they are in the wrong array */
	  q[0].length=q[1].length;
	  F77CALL (dcopy) (&q[0].length,q[1].data,&inc1,q[0].data,&inc1);

	  q_up[0].length=q_up[1].length;
	  F77CALL (dcopy) (&q_up[0].length,q_up[1].data,&inc1,q_up[0].data,&inc1);
	}
        else
	{
	  /* copying last q into q_lef */
	  q_lef.length=q[1].length;
	  F77CALL (dcopy) (&q[1].length,q[1].data,&inc1,q_lef.data,&inc1);
	}
      }
    }
    else /* n column regular */
    {
      regular=0;

      k=1;

      e_up(0,0)=p(0,0);

      e_up(1,0)=lau_coef(&h,&q_up[0],n+k);

/*       e_lef[1]=lau_coef(&h,&q_lef,n+k); */
      e_lef[1] = get_e_lef_k(&q_lef,q_up,e_lef,&e_up(0,0),&ldsmall,k);

/*       e(0,0) = lau_coef(&h,&q[0],n+1); */
      e(0,0) = (e_lef[1]*e_up(0,0) - e_up(1,0)*e_lef[0])/p_lef[0]; 

#ifdef CALCULATE_WITH_D
      update_d(&d(0,0),&ldsmall,q->data,&ldq,&e_up(0,0),&ldsmall,n,k);
#endif

      while ((k<k_max) && (k+n<N) && !regular)
      {
#if (DEBUG >= 1)
        fprintf(stdout,"n = %d, k = %d : Making look-ahead step\n",
		(int)n,(int)k);
#endif
	v[k]    = -e(k-1,k-1)/e_lef[0];
        v_up[k] = -p_up(k-1,k-1)/p_lef[0];

        get_inner_q(&q[k],&q_lef,&q[k-1],&v[k],k);

	get_inner_q_up(&q_up[k],&q_lef,&q_up[k-1],&v_up[k],k);

        get_anti_reg_coefs(&e(0,0),&e_up(0,0),&p(0,0),&p_up(0,0),p_lef,e_lef,
		      &ldsmall,&h,&q[0],&q_up[0],&q_lef,v,v_up,k,n);

        get_new_coefs(&e(0,0),&e_up(0,0),&p(0,0),&p_up(0,0),p_lef,e_lef,
		      &ldsmall,&h,&q[0],&q_up[0],&q_lef,v,v_up,k,n);

	k++;

#ifdef CALCULATE_WITH_D
        update_d(&d(0,0),&ldsmall,q->data,&ldq,&e_up(0,0),&ldsmall,n,k);
        help = get_singularity(work,lwork,&d(0,0),&ldsmall,k);
#else
        help = get_singularity(work,lwork,&p(0,0),&ldsmall,k);
        help2 = get_singularity(work,lwork,&e_up(0,0),&ldsmall,k);
        help = min(help,help2);
#endif

#if (DEBUG>=1)
	fprintf(stdout,"n = %d, k = %d : Found smallest singular value %e\n",
		(int)n,(int)k,help);
#endif

	/* storing the best singular value up to now */
	set_new_singular_value(k,help);

	regular = (help > tol);
      } /* while(k<=k_max) */

      if (!regular) /* n+k not regular */
      {
        /* getting the optimal singular value uo to now */
	get_opt_singular_value(&k,&help);

#if (DEBUG >= 1)
	fprintf(stdout,"n = %d, k = %d : The best singular value up to now was"
		" %e\n",(int)n,(int)k,help);
#endif

	if (n+k == N) /* end of matrix reached */
	{
	  if (help < sing_tol) /* singular */
	  {
	    fprintf(stderr,"Matrix is close to singular, no solution"
		    " computed\n");
	    singular = 1;
	  }
#if (DEBUG >= 1)
	  else
	    fprintf(stderr,"Matrix may be ill conditioned, last"
		    " singular value is %e\n",help);
#endif
	}

        if (k==1)
	{
#if (DEBUG >= 1)
	  (number_of_steps[0])++;
          fprintf(stdout,"n = %d : Accepted classical step with p = %e. "
		  "updating x\n",(int)n,p(0,0));
#endif

	  /* last regular pair now treated as column-regular */
	  column_regular = 1;

/*           e(0,0)=lau_coef(&h,&q[0],n+1); */

	  /* update solution */
	  if (!update_x(work,lwork,&q[0],&q_up[0],&ldq,b,x,&p(0,0),&ldsmall,
			NULL,NULL,n,1)) break;
          storage_index = 0;

	  if ((n+1 == N) && refine)
	    set_refinement(&h,q,q_up,NULL,&p(0,0));

	  continue;
	}
      }  /* n+k not regular */

#if (DEBUG >= 1)
      (number_of_steps[k-1])++;
      fprintf(stdout,"n = %d, k = %d : Updating x from the last look-ahead "
	      "step\n",(int)n,(int)k);
#endif

      /* update solution */
#ifndef CALCULATE_WITH_D
      F77CALL (dgemm) ("Transpose","No transpose",&k,&k,&k,&one,&q->data[n],
             &ldq,&e_up(0,0),&ldsmall,&zero,&d(0,0),&ldsmall);
#endif
      if (!update_x(work,lwork,&q[0],&q_up[0],&ldq,b,x,&d(0,0),&ldsmall,
		    d_perm,d_tau,n,k))
	break;

      if ((n+k < N) || (refine && (n+k == N)))
      {
#if (DEBUG >= 3)
	fprintf(stdout,"%3scalculating new q_up ,q_lef\n","");
#endif
#ifndef CALCULATE_WITH_D
	/* calculating anti-regular q_up */
	if (!get_anti_regular_up(work,lwork,&q_lef,q_up,q,&e_up(0,0),e_lef,
				 &p_up(0,0),NULL,d_perm,d_tau,p_lef,
				 &ldsmall,&ldq,n,k))  break;

	/* calculating new q_lef */
	if (!get_anti_regular_lef(work,lwork,q,&q_lef,&e(0,0),e_lef,&p(0,0),
				  NULL,d_perm,d_tau,&ldsmall,&ldq,n,k))
	  break;
#else
	/* calculating anti-regular q_up */
	if (!get_anti_regular_up(work,lwork,&q_lef,q_up,q,&e_up(0,0),e_lef,
				 &p_up(0,0),&d(0,0),d_perm,d_tau,p_lef,
				 &ldsmall,&ldq,n,k))  break;

	/* calculating new q_lef */
	if (!get_anti_regular_lef(work,lwork,q,&q_lef,&e(0,0),e_lef,&p(0,0),
				  &d(0,0),d_perm,d_tau,&ldsmall,&ldq,n,k))
	  break;
#endif

/* 	e_up(0,0) = lau_coef(&h,&q_up[0],n+k); */
/*      is now computed in get_anti_regular_up */
/* 	e_lef[0]  = lau_coef(&h,&q_lef,n+k); */
/*      is now computed in get_anti_regular_lef */

        p_up(0,0) = lau_coef(&h,&q_up[0],-1);
	p_lef[0]  = 1.0;

        if ((n+k == N) && refine)
	  set_refinement(&h,NULL,q_up,&q_lef,&p(0,0));

	if (n+k < N)
	{
#if (DEBUG >= 3)
	  fprintf(stdout,"%3scalculating new q\n","");
#endif
	  /* calculating new q */
	  get_regular_q(q,q_up,&q_lef,e_lef,&e_up(0,0));

	  p(0,0)    = e_up(0,0); 

	  set_new_singular_value(1,p(0,0));

	  /* check if new pair is column regular */
	  column_regular = (col_tol < fabs(p(0,0)));

	  if (column_regular || ((n+k+1 == N) && (sing_tol < fabs(p(0,0)))) )
	  {
#if (DEBUG >= 1)
	    (number_of_steps[0])++;
            if (column_regular)
	      fprintf(stdout,"n = %d : New column pair is column regular "
		      "with p = %e updating x\n",(int)(n+k),p(0,0));
	    else
	      fprintf(stdout,"n = %d : New column pair is not column regular "
		    "with p = %e updating x anyway\n",(int)(n+k),p(0,0));
#endif
	    /* preparing for the next classical step */
	    storage_index=0;

	    e(0,0) = lau_coef(&h,&q[storage_index],n+k+1);

	    /* updating solution with the new column regular pair */
	    if (!update_x(work,lwork,&q[storage_index],&q_up[storage_index],
			  &ldq,b,x,&p(0,0),&ldsmall,NULL,NULL,n+k,1)) break;

	    if (refine && (n+k+1 == N))
	      set_refinement(&h,&q[storage_index],&q_up[storage_index],NULL,
			     &p(0,0));

	  }  /* if column_regular */
	  if (!column_regular && (n+k+1 == N))
	  {
	    if (fabs(p(0,0)) < sing_tol)
	    {
	      fprintf(stderr,"Matrix is close to singular, no solution"
		      " computed\n");
	      singular = 1;
	    }
#if (DEBUG >= 1)
	    else
	      fprintf(stderr,"Matrix may be ill conditioned, last"
		      " singular value is %e\n",fabs(p(0,0)));
#endif
    }

	} /* */
      } /* if ((n+k < N) || refine) */

      /* updating n */
      n+=k;

    }  /* else   n column regular */
    psol[n] = x[n];
  } /* main loop */

  if (!singular && refine)
  {
#if (DEBUG >= 1)
    fprintf(stdout,"Refining solution\n");
#endif
    refine_solution(x,b);
  }

#if (DEBUG >= 1)
  fprintf(stdout,"The number of steps:\n");
  for (i=0;i<k_max;i++)
    if (number_of_steps[i] != 0)
      fprintf(stdout,"\t%6d steps of size %3d\n",(int)number_of_steps[i],
	      (int)(i+1));
  fflush(stdout);
#endif

  free(data1);
  free(data2);

  return(success && !singular);

#undef p_up
#undef p
#undef e_up
#undef e
}
