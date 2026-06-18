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
static finteger incm1 = -1;

int d_schurlev(fdouble *x,fdouble *mu,finteger N,fdouble *b,
	     fdouble col_tol,fdouble tol,int k_max,int refine)
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
  finteger k_act;
  finteger success             = 1;
  finteger regular             = 0;
  finteger column_regular      = 0;
  finteger singular            = 0;
  finteger ihelp,i;

  /* variable which indicates in which polynom the current classical 
     polynomials are */
  finteger storage_index=0;

  /* tolerance for singular matrix */
  const fdouble sing_tol=SING_TOL;

  /* help variables */
  fdouble help;
  fdouble *work;
  finteger lwork;

  /* variables for the method */
  polynom *q,*q_up,q_lef;
  finteger ld=N+1;
  toeplitz h;
  fdouble *p_up__,*p__;
  fdouble *e_up__,*e__;
  fdouble *d__;
  fdouble *p_lef,*e_lef;
  fdouble *v,*v_up;
  fdouble *M;
  fdouble *b_tmp,*y;
  finteger ldsmall=k_max;
  finteger ldM=2*k_max;

#define p_up(I,J)   p_up__[(I) + (J) * ld]
#define p(I,J)      p__[(I) + (J) * ld]
#define e_up(I,J)   e_up__[(I) + (J) * ld]
#define e(I,J)      e__[(I) + (J) * ld]
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

  /* memory management */

  /* space for q and q_up */
  ihelp=sizeof(fdouble)*ld*2*max(k_max,2);
  /* space for q_lef */
  ihelp+=ld*sizeof(fdouble);
  /* space for e,e_up,p and p_up */
  ihelp+=sizeof(fdouble)*ld*4*max(k_max,2);
  /* space for e_lef,p_lef */
  ihelp+=sizeof(fdouble)*ld*2;
  /* space for updating temporary right-hand side */
  ihelp+=2*ld*sizeof(fdouble);
  /* space for v,v_up */
  ihelp+=sizeof(fdouble)*(k_max+1)*2;
  /* space for M */
  ihelp+=ldM*ldM*sizeof(fdouble);
  /* space for d */
  ihelp+=ldsmall*ldsmall*sizeof(fdouble);
  /* space for temporary variables */
  lwork=max(k_max*k_max+k_max+liwork(k_max),6*ld+9*k_max+liwork(2*k_max));
#ifdef USE_QR_DECOMP
  lwork+=k_max*3;
#endif
  lwork=max(lwork,4*k_max*k_max+12*k_max)+3;
  ihelp+=lwork*sizeof(fdouble);
  q[0].data=(fdouble *)(data2=malloc(ihelp));
  if (data2==NULL)
  {
    fprintf(stderr,"error can't allocate memory for working!\nAbort ...\n");
    free(data1);
    return(0);
  }

  /* distribute memory */
  q_up[0].data=&(q[0].data[max(k_max,2)*ld]);
  for(i=1;i<max(k_max,2);i++)
  {
    q[i].data=&q[0].data[ld*i];
    q_up[i].data=&q[0].data[ld*(i+max(k_max,2))];
  }
  q_lef.data=&q_up[max(k_max,2)-1].data[ld];
  p__=&q_lef.data[ld];
  p_up__=&p__[max(2,k_max)*ld];
  e__=&p_up__[max(2,k_max)*ld];
  e_up__=&e__[max(2,k_max)*ld];
  d__=&e_up__[max(2,k_max)*ld];
  M=&d__[ldsmall*ldsmall];
  e_lef=&M[ldM*ldM];
  p_lef=&e_lef[ld];
  v=&p_lef[ld];
  v_up=&v[k_max+1];
  b_tmp=&v_up[k_max+1];
  y=&b_tmp[ld];
  work=&y[ld];

  /* initialize memory */
  help=0.0;
  F77CALL (dcopy) (&N,&help,&inc0,x,&inc1);
  F77CALL (dcopy) (&ld,&help,&inc0,y,&inc1);
  ihelp=2*max(k_max,2)*ld;
  F77CALL (dcopy) (&ihelp,&help,&inc0,q->data,&inc1);

  e_lef++;
  p_lef++;

#if (DEBUG >= 1)
  fprintf(stdout,"******************************************************\n"
	         "*       New look-ahead Schur algorithm               *\n"
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

  F77CALL (dcopy) (&h.N,&h.data[-h.N+1],&incm1,&p(0,0),&inc1);
  F77CALL (dcopy) (&h.N,h.data,&inc1,&e_up(0,0),&inc1);
  ihelp = h.N-1;
  F77CALL (dcopy) (&ihelp,&h.data[-h.N+1],&incm1,&p_up(0,0),&inc1);
  F77CALL (dcopy) (&ihelp,&h.data[1],&inc1,&e(0,0),&inc1);

  F77CALL (dcopy) (&N,b,&inc1,b_tmp,&inc1);

  if ((col_tol < fabs(e_up(0,0))) || (N==1))
  {
    /* update the solution */
    update_y(work,lwork,b_tmp,y,&e_up(0,0),&ld,N,0,1);
    update_x_from_y(q_up,y,x,0,1,&ld);

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
    help=0.0;
    F77CALL (dcopy) (&N,&help,&inc0,e_lef,&inc1);
    F77CALL (dcopy) (&N,&help,&inc0,p_lef,&inc1);
    p_lef[-1]=1.0;
    e_lef[-1]=-1.0;

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


  /*************************/
  /*     the main loop     */
  /*************************/
  while(n<N-1)
  {
    if (column_regular) /* n column regular */
    {
      ihelp=(storage_index+1)%2;
      get_colreg_pair_schurlev(&q[ihelp],&q_up[ihelp],&h,
		      &q[storage_index],&q_up[storage_index],
		      &e(n,storage_index),&e_up(n,storage_index),
		      &e(n,ihelp),&e_up(n,ihelp),
		      &p(n,storage_index),&p_up(n,storage_index),
		      &p(n,ihelp),&p_up(n,ihelp),
		      &ld,v,v_up,n);
      storage_index=ihelp;

      n++;

      if (fabs(p(n,storage_index)) > col_tol)  /* n column regular */
      {
	if (!update_y(work,lwork,&b_tmp[n],&y[n],&e_up(n,storage_index),
		      &ld,N,n,1))  break;
	update_x_from_y(&q_up[storage_index],&y[n],x,n,1,&ld);

        if (refine && (n+1 == N))
          set_refinement(&h,&q[storage_index],&q_up[storage_index],NULL,
			 &p(n,storage_index));
#if (DEBUG >= 2)
        fprintf(stdout,"%2sn = %d : Making column regular step with p = %e\n",
		"",(int)n,p(n,storage_index));
#endif
#if (DEBUG >= 1)
	(number_of_steps[0])++;
#endif
      }
      else
      {
#if (DEBUG >= 1)
        fprintf(stdout,"n = %d : Step is not column regular with p = %e\n",
		(int)n,p(n,storage_index));
#endif
	/* setting the flag */
	column_regular=0;

        /* test if the end of the matrix is reached */
	if (n+1 == N)
	{
	  if (fabs(p(n,storage_index)) < sing_tol)
	  {
	    fprintf(stderr,"Matrix is close to singular. No solution "
		    "computed\n");
            singular=1;
	  }
	  else
	  {
#if (DEBUG >= 1)
	    (number_of_steps[0])++;
#endif
	    fprintf(stderr,"Matrix may be ill conditioned, updating x, last"
		    " singular value is %e\n",fabs(p(n,storage_index)));

	    if (!update_y(work,lwork,&b_tmp[n],&y[n],&e_up(n,storage_index),
			  &ld,N,n,1))  break;
	    update_x_from_y(&q_up[storage_index],&y[n],x,n,1,&ld);

	    if (refine)
	      set_refinement(&h,&q[storage_index],&q_up[storage_index],NULL,
			     &p(n,storage_index));
	  }

	  break;
	}

        /* storing the best singular value up to now */
	set_new_singular_value(1,p(n,storage_index));

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

	  /* copy old p,e to p_lef,e_lef */
          ihelp = N - n + 2;
	  F77CALL (dcopy) (&ihelp,&p(n-2,0),&inc1,&p_lef[n-2],&inc1);
	  F77CALL (dcopy) (&ihelp,&e(n-2,0),&inc1,&e_lef[n-2],&inc1);

	  /* copying e,p,e_up,p_up if they are in the wrong array */
	  ihelp = h.N-n-1;
	  F77CALL (dcopy) (&ihelp,&p_up(n,1),&inc1,&p_up(n,0),&inc1);
	  F77CALL (dcopy) (&ihelp,&e(n,1),&inc1,&e(n,0),&inc1);
	  ihelp = h.N-n;
	  F77CALL (dcopy) (&ihelp,&p(n,1),&inc1,&p(n,0),&inc1);
	  F77CALL (dcopy) (&ihelp,&e_up(n,1),&inc1,&e_up(n,0),&inc1);
	}
        else
	{
	  /* copying last q into q_lef */
	  q_lef.length=q[1].length;
	  F77CALL (dcopy) (&q[1].length,q[1].data,&inc1,q_lef.data,&inc1);

	  /* copy old p,e to p_lef,e_lef */
          ihelp = N - n + 2;
	  F77CALL (dcopy) (&ihelp,&p(n-2,1),&inc1,&p_lef[n-2],&inc1);
	  F77CALL (dcopy) (&ihelp,&e(n-2,1),&inc1,&e_lef[n-2],&inc1);
	}
      }
    }
    else /* n column regular */
    {
      regular=0;

      k=1;

      while ((k<k_max) && (k+n<N) && !regular)
      {
        k++;
        get_M(M,&ldM,&e_lef[n-1],&e_up(n,0),&p_lef[n-1],&p_up(n,0),k);
        help = get_singularity(work,lwork,M,&ldM,2*k-1);
#if (DEBUG>=1)
        fprintf(stdout,"n = %d, k = %d : Found smallest singular value %e\n",
                (int)n,(int)k,help);
#endif  
        set_new_singular_value(k,help);
        regular = (help > tol);
      }
 
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
	  else
	    fprintf(stderr,"Matrix may be ill conditioned, last"
		    " singular value is %e\n",help);
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

	  /* update solution */
	  if (!update_y(work,lwork,&b_tmp[n],&y[n],&e_up(n,0),&ld,N,n,1))break;
	  update_x_from_y(&q_up[0],&y[n],x,n,1,&ld);

          storage_index = 0;

	  if ((n+1 == N) && refine)
	    set_refinement(&h,q,q_up,NULL,&p(n,0));

	  continue;
	}
	get_M(M,&ldM,&e_lef[n-1],&e_up(n,0),&p_lef[n-1],&p_up(n,0),k);
      }  /* n+k not regular */

      k_act = k;
      k = 1;

      update_d_schur(&d(0,0),&ldsmall,&e_up(n,0),&ld,&p(n,0),&ld,k);

      while (k<k_act)
      {
#if (DEBUG >= 1)
        fprintf(stdout,"n = %d, k = %d : Making look-ahead step\n",
		(int)n,(int)k);
#endif

	v[k]    = -e(n+k-1,k-1)/e_lef[n-1];
        v_up[k] = -p_up(n+k-1,k-1)/p_lef[n-1];

        get_inner_q_schurlev(&q[k],&q_lef,&q[k-1],&e(n+k,k),&e_lef[n-1],
                 &e(n+k-1,k-1),&p(n,k),&p_lef[n-1],&p(n,k-1),&v[k],N,n,k);

	get_inner_q_up_schurlev(&q_up[k],&q_lef,&q_up[k-1],&e_up(n,k),&e_lef[n-1],
		    &e_up(n,k-1),&p_up(n+k,k),&p_lef[n-1],
		    &p_up(n+k-1,k-1),&v_up[k],N,n,k);

	k++;

        update_d_schur(&d(0,0),&ldsmall,&e_up(n,0),&ld,&p(n,0),&ld,k);

#if (DEBUG>=1)
	fprintf(stdout,"n = %d, k = %d : Found smallest singular value %e\n",
		(int)n,(int)k,help);
#endif
      } /* while(k<=k_act) */

#if (DEBUG >= 1)
      (number_of_steps[k-1])++;
      fprintf(stdout,"n = %d, k = %d : Updating x from the last look-ahead "
	      "step\n",(int)n,(int)k);
#endif

      /* update solution */
      if (!update_y(work,lwork,&b_tmp[n],&y[n],&e_up(n,0),&ld,N,n,k))  break;
      update_x_from_y(&q_up[0],&y[n],x,n,k,&ld);


      if ((n+k < N) || (refine && (n+k == N)))
      {
#if (DEBUG >= 3)
	fprintf(stdout,"%3scalculating new q_up ,q_lef\n","");
#endif

	/* calculating anti-regular pair */
	get_anti_regular_pair_schurlev(work,lwork,M,&ldM,&q_lef,q_up,&ld,
	  &e_up(n,0),&e_lef[n-1],&p_up(n,0),&p_lef[n-1],NULL,N,n,k);

        if ((n+k == N) && refine)
	  set_refinement(&h,NULL,q_up,&q_lef,&p(n+k,0));

	if (n+k < N)
	{
#if (DEBUG >= 3)
	  fprintf(stdout,"%3scalculating new q\n","");
#endif
	  /* calculating new q */
	  get_regular_q_schurlev(q,q_up,&q_lef,&e(n+k,0),&e_lef[n+k-1],&e_up(n+k,0),
			&p(n+k,0),&p_lef[n+k-1],&p_up(n+k,0),N,n,k);

	  set_new_singular_value(1,p(n+k,0));

	  /* check if new pair is column regular */
	  column_regular = (col_tol < fabs(p(n+k,0)));

	  if (column_regular || ((n+k+1 == N) && (sing_tol < fabs(p(n+k,0)))) )
	  {
#if (DEBUG >= 1)
	    (number_of_steps[0])++;
            if (column_regular)
	      fprintf(stdout,"n = %d : New column pair is column regular "
		      "with p = %e updating x\n",(int)(n+k),p(n+k,0));
	    else
	      fprintf(stdout,"n = %d : New column pair is not column regular "
		    "with p = %e updating x anyway\n",(int)(n+k),p(n+k,0));
#endif

	    /* preparing for the next classical step */
	    storage_index=0;

	    /* updating solution with the new column regular pair */
	    if (!update_y(work,lwork,&b_tmp[n+k],&y[n+k],&e_up(n+k,0),
			  &ld,N,n+k,1)) break;
	    update_x_from_y(&q_up[storage_index],&y[n+k],x,n+k,1,&ld);

	    if (refine && (n+k+1 == N))
	      set_refinement(&h,&q[storage_index],&q_up[storage_index],NULL,
			     &p(n+k,0));

	  }  /* if column_regular */
	  if (!column_regular && (n+k+1 == N))
	  {
	    if (fabs(p(n+k,0)) < sing_tol)
	    {
	      fprintf(stderr,"Matrix is close to singular, no solution"
		      " computed\n");
	      singular = 1;
	    }
	    else
	      fprintf(stderr,"Matrix may be ill conditioned, last"
		      " singular value is %e\n",fabs(p(n+k,0)));
	  }

	} /* */
      } /* if ((n+k < N) || refine) */

      /* updating n */
      n+=k;

    }  /* else   n column regular */
  } /* main loop */

  if (!singular && refine)
  {
#if (DEBUG >= 1)
    fprintf(stdout,"Refining solution\n");
#endif
    refine_solution(x,b);
  }

  free(data1);
  free(data2);

#if (DEBUG >= 1)
  fprintf(stdout,"The number of steps:\n");
  for (i=0;i<k_max;i++)
    if (number_of_steps[i] != 0)
      fprintf(stdout,"\t%6d steps of size %3d\n",(int)number_of_steps[i],
	      (int)(i+1));
  fflush(stdout);
#endif

  return(success && !singular);

#undef p_up
#undef p
#undef e_up
#undef e
#undef d
}
