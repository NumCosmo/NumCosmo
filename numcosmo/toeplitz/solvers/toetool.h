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

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

extern void d_toe_mv(finteger N,finteger first,finteger last,fdouble *alpha,
       fdouble *a,fdouble *x,fdouble *beta,fdouble *y);

extern void d_toe_mm(finteger N,finteger first,finteger last,finteger m,
       fdouble *alpha,fdouble *a,fdouble *b__,finteger ldb,
       fdouble *beta,fdouble *c__,finteger ldc);

extern void d_toe2full(finteger N,fdouble *to,fdouble *fu,finteger *ld);

#ifdef  __cplusplus
}
#endif
