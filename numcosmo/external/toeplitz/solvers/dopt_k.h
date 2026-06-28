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

extern void set_new_singular_value(finteger k, fdouble sing_val);

extern void get_opt_singular_value(finteger *k,fdouble *sing_val);

#ifdef  __cplusplus
}
#endif
