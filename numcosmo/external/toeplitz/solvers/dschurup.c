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
#include "dopt_k.h"

#define SING_TOL   1e-12

static finteger inc1  =  1;
static finteger inc0  =  0;
static finteger incm1 = -1;


int d_schurup(fdouble *x,fdouble *mu,finteger N,fdouble *b,
      fdouble col_tol,fdouble tol,int k_max,int dummy)
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

  void *data2;

  /* finteger help variables */
  finteger n                   = 0;
  finteger k                   = 0;
  finteger k_act;
  finteger success             = 1;
  finteger regular             = 0;
  finteger column_regular      = 0;
  finteger singular            = 0;
  finteger ihelp;
#if (DEBUG >= 1)
  finteger i;
#endif

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
  finteger ld=N+1;
  toeplitz h;
  fdouble *p_up__,*p__;
  fdouble *e_up__,*e__;
  fdouble *d__;
  fdouble *p_lef,*e_lef;
  fdouble *v,*v_up;
  fdouble *M;
  fdouble *b_tmp,*y;
  fdouble *schur_par;
  finteger ldsmall=k_max;
  finteger ldM=2*k_max;
  finteger *reg_index;
  finteger reg_pos=0;

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

  /* initialize h */
  h.data=mu;
  h.N=N;

  /* space for e,e_up,p and p_up */
  ihelp=sizeof(fdouble)*ld*4*max(k_max,2);
  /* space for d */
  ihelp+=sizeof(fdouble)*ldsmall*ldsmall;
  /* space for M */
  ihelp+=sizeof(fdouble)*ldM*ldM;
  /* space for e_lef,p_lef */
  ihelp+=sizeof(fdouble)*ld*2;
  /* space for v,v_up */
  ihelp+=sizeof(fdouble)*(k_max+1)*2;
  /* space for updating temporary right-hand side */
  ihelp+=2*ld*sizeof(fdouble);
  /* space for storing the parameters for schur updating formulas */
  ihelp+=5*N*sizeof(fdouble);
  /* space (finteger) for storing regular indices */
  ihelp+=ld*sizeof(finteger);
  /* space for temporary variables */
  lwork=max(k_max*k_max + k_max + liwork(k_max),6*ld+9*k_max+liwork(2*k_max));
#ifdef USE_QR_DECOMP
  lwork+=k_max*3;
#endif
  lwork=max(lwork,4*k_max*k_max+12*k_max);
  ihelp+=lwork*sizeof(fdouble);
  p__=(fdouble *)(data2=malloc(ihelp));
  if (data2==NULL)
  {
    fprintf(stderr,"error can't allocate memory for working!\nAbort ...\n");
    return(0);
  }

  /* distribute memory */
  p_up__=&p__[max(2,k_max)*ld];
  e__=&p_up__[max(2,k_max)*ld];
  e_up__=&e__[max(2,k_max)*ld];
  d__=&e_up__[max(2,k_max)*ld];
  M=&d__[ldsmall*ldsmall];
  v=&M[ldM*ldM];
  v_up=&v[k_max+1];
  e_lef=&v_up[k_max+1];
  p_lef=&e_lef[ld];
  b_tmp=&p_lef[ld];
  y=&b_tmp[ld];
  schur_par=&y[ld];
  work=&schur_par[5*N];
  reg_index=(finteger*)&work[lwork];

  /* initialize memory */
  help=0.0;
/*   F77CALL (dcopy) (&N,&help,&inc0,x,&inc1); */
  F77CALL (dcopy) (&ld,&help,&inc0,y,&inc1);

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
  F77CALL (dcopy) (&h.N,&h.data[-h.N+1],&incm1,&p(0,0),&inc1);
  F77CALL (dcopy) (&h.N,h.data,&inc1,&e_up(0,0),&inc1);
  ihelp = h.N-1;
  F77CALL (dcopy) (&ihelp,&h.data[-h.N+1],&incm1,&p_up(0,0),&inc1);
  F77CALL (dcopy) (&ihelp,&h.data[1],&inc1,&e(0,0),&inc1);

  F77CALL (dcopy) (&N,b,&inc1,b_tmp,&inc1);

  reg_index[reg_pos++] = 0;

  if ((col_tol < fabs(e_up(0,0))) || (N==1))
  {
    /* update the solution */
    update_y_schur(work,lwork,b_tmp,y,&e_up(0,0),N,0,1,&ld);

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
    /* preparing p_lef and e_lef for the first look-ahead step */
    help = 0.0;
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
      get_colreg_pair_schurlev(NULL,NULL,&h,
	    NULL,NULL,
	    &e(n,storage_index),&e_up(n,storage_index),
	    &e(n,ihelp),&e_up(n,ihelp),
	    &p(n,storage_index),&p_up(n,storage_index),
	    &p(n,ihelp),&p_up(n,ihelp),
	    &ld,v,v_up,n);
      storage_index=ihelp;

      CR_GAMMA(n)=*v;
      CR_GAMMA_UP(n)=*v_up;

      n++;

      reg_index[reg_pos++] = n;

      if (fabs(p(n,storage_index)) > col_tol)  /* n column regular */
      {

	if (!update_y_schur(work,lwork,&b_tmp[n],&y[n],&e_up(n,storage_index),
	      N,n,1,&ld))  break;

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

	CR_GAMMA(n-1)=0.0;

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

	    if (!update_y_schur(work,lwork,&b_tmp[n],&y[n],
		  &e_up(n,storage_index),N,n,1,&ld))  break;

	  }

	  break;
	}

        /* storing the best singular value up to now */
	set_new_singular_value(1,p(n,storage_index));

        /* preparing for the first look-ahead step */
	if (storage_index == 1)
	{
	  /* copy old p,e to p_lef,e_lef */
	  ihelp = h.N - n + 2;
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
	  /* copy old p,e to p_lef,e_lef */
	  ihelp = h.N - n + 2;
	  F77CALL (dcopy) (&ihelp,&p(n-2,1),&inc1,&p_lef[n-2],&inc1);
	  F77CALL (dcopy) (&ihelp,&e(n-2,1),&inc1,&e_lef[n-2],&inc1);
	}
	ALPHA(n) = e_up(n,0)/p_lef[n-1];
	BETA(n)  = -e_lef[n-1]/p_lef[n-1];
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
        set_new_singular_value(k,help);
        regular = (help > tol);
#if (DEBUG>=1)
        fprintf(stdout,"n = %d, k = %d : Found smallest singular value %e\n",
	      (int)n,(int)k,help);
#endif  
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

	  /* restore last gamma if previous step was column-regular */
          if (n-reg_index[reg_pos-2] == 1)
            CR_GAMMA(n-1)=*v;

	  /* update solution */
	  if (!update_y_schur(work,lwork,&b_tmp[n],&y[n],&e_up(n,0),
		N,n,1,&ld))  break;

          storage_index = 0;

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

	GAMMA_UP(n+k-1) = v_up[k];

        get_inner_q_schurlev(NULL,NULL,NULL,NULL,NULL,NULL,&p(n,k),&p_lef[n-1],
	      &p(n,k-1),&v[k],N,n,k);

	get_inner_q_up_schurlev(NULL,NULL,NULL,&e_up(n,k),&e_lef[n-1],
	      &e_up(n,k-1),NULL,NULL,NULL,&v_up[k],N,n,k);

	get_new_e_pup_coefs(&e(n,0),&p_up(n,0),&p_lef[n-1],&e_lef[n-1],
	      &ld,v,v_up,k,n);

	k++;

	update_d_schur(&d(0,0),&ldsmall,&e_up(n,0),&ld,&p(n,0),&ld,k);
      } /* while(k<=k_act) */

#if (DEBUG >= 1)
      (number_of_steps[k-1])++;
      fprintf(stdout,"n = %d, k = %d : Updating x from the last look-ahead "
	    "step\n",(int)n,(int)k);
#endif
      if (!update_y_schur(work,lwork,&b_tmp[n],&y[n],&e_up(n,0),
	    N,n,k,&ld))  break;

      if (n+k < N)
      {
#if (DEBUG >= 3)
	fprintf(stdout,"%3scalculating new q_up ,q_lef\n","");
#endif

	/* calculating anti-regular pair */
	get_anti_regular_pair_schurlev(work,lwork,M,&ldM,NULL,NULL,&ld,
	      &e_up(n,0),&e_lef[n-1],&p_up(n,0),&p_lef[n-1],schur_par,N,n,k);

	reg_index[reg_pos++] = n+k;

#if (DEBUG >= 3)
	fprintf(stdout,"%3scalculating new q\n","");
#endif

	ALPHA(n+k) = e_up(n+k,0);

	BETA(n+k)  = -e_lef[n+k-1];

	/* calculating new q */
	get_regular_q_schurlev(NULL,NULL,NULL,&e(n+k,0),&e_lef[n+k-1],
	      &e_up(n+k,0),&p(n+k,0),&p_lef[n+k-1],&p_up(n+k,0),N,n,k);

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
	  if (!update_y_schur(work,lwork,&b_tmp[n+k],&y[n+k],&e_up(n+k,0),
		N,n+k,1,&ld)) break;

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

      } /* if (n+k < N) */

      /* updating n */
      n+=k;
    }  /* else   n column regular */
  } /* main loop */

  reg_index[reg_pos] = N;

  update_x_from_y_schur(work,lwork,schur_par,y,x,N,reg_index,reg_pos,&ld);

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
#undef L
#undef e_up
#undef e
#undef d
}
