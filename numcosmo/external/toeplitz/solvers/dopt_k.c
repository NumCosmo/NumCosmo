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
#include <math.h>
#include "ftypes.h"
#include "dopt_k.h"

typedef struct
{
  finteger k;
  fdouble singular_value;
} optimal_k;

static optimal_k opt_k;

void set_new_singular_value(finteger k, fdouble sing_val)
{
  sing_val = fabs(sing_val);

  if (k==1)
  {

    opt_k.k              = k;
    opt_k.singular_value = sing_val;
#if (DEBUG>=5)
    fprintf(stdout,"%5sopt_k is reset with smallest singular value %e\n","",
	    opt_k.singular_value);
#endif
  }
  else
    if (opt_k.singular_value < sing_val)
    {
      opt_k.k              = k;
      opt_k.singular_value = sing_val;
#if (DEBUG>=5)
      fprintf(stdout,"%5sopt_k is set to the new k = %d with singular value "
	      "%e\n","",(int)opt_k.k,opt_k.singular_value);
#endif
    }
#if (DEBUG>=5)
    else
    {
      fprintf(stdout,"%5sopt_k stays at the old value k = %d at new k = %d\n",
	      "",(int)opt_k.k,(int)k);
    }
#endif
}

void get_opt_singular_value(finteger *k,fdouble *sing_val)
{
  *k        = opt_k.k;
  *sing_val = opt_k.singular_value;
#if (DEBUG>=5)
    fprintf(stdout,"%5sopt_k returning the best singular value %e\n","",
	    opt_k.singular_value);
#endif
}
