/*************************************************************************/
/* Copyright (C) 1996-1997 Marlis Hochbruck                              */
/* All rights reserved.                                                  */
/*                                                                       */
/* C code written by Mathias Froehlich.                                  */
/*                                                                       */
/* This code is part of a copyrighted package. For details, see the file */
/* "COPYING" in the top-level directory.                                 */
/*************************************************************************/
#include "ftypes.h"

#ifdef  __cplusplus
extern "C" { 
#endif

extern void refine_solution(fdouble *x,fdouble *b);

extern void invert_toeplitz(fdouble *x,fdouble *b);

extern int malloc_for_refinement(finteger dimension);

extern void free_refinement_wrokspace();

extern int set_refinement(toeplitz *h_in,polynom *q,polynom *q_up,
       polynom *q_lef,fdouble *p_in);

#ifdef  __cplusplus
}
#endif
