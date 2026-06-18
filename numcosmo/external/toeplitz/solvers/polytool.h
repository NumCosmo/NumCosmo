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
#include "toeplitz.h"

#ifdef  __cplusplus
extern "C" { 
#endif

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

extern fdouble lau_coef(toeplitz *h,polynom *q,finteger k);

extern fdouble  poly_coef(polynom *p,finteger k);

extern finteger poly_index(polynom *p,finteger k);

extern fdouble *data_pointer(polynom *p,finteger k1,finteger k2);

extern finteger length_of_data(polynom *p,finteger k1,finteger *clip1,
			      finteger k2,finteger *clip2);

extern void mult_by_z_up_k(polynom *res,polynom *p,fdouble *alpha,
			   finteger k);

extern void polynomaxpy(fdouble *alpha,polynom *x,polynom *y);

extern void polynommult(polynom *res,polynom *x,polynom *y);

#ifdef  __cplusplus
}
#endif
