/*************************************************************************/
/* Copyright (C) 1996-1997 Marlis Hochbruck                              */
/* All rights reserved.                                                  */
/*                                                                       */
/* C code written by Mathias Froehlich.                                  */
/*                                                                       */
/* This code is part of a copyrighted package. For details, see the file */
/* "COPYING" in the top-level directory.                                 */
/*************************************************************************/
#ifndef TOEPLITZ_HEADER
#define TOEPLITZ_HEADER
#include "ftypes.h"

/*
 * This header file containes all Toeplitz solver subroutines
 * which are important for a user.
 * All these solvers share the way how the matrix coeficients are stored 
 * in memory.
 * Given a Toeplitz matrix 
 * 
 *            | a_0     a_1    a_2  ...  a_{n-1} |
 *            | a_{-1}  a_0    a_1  ...  a_{n-2} |
 *            | a_{-2}  a_{-1} a_0  ...  a_{n-3} |
 *        A = |  .                        .      |
 *            |  .                        .      |
 *            |  .                        .      |
 *            | a_{1-n}        ...       a_0     |
 *
 * the solvers expect this matrix in an array, where the 
 * entries are stored in the order of the indizes above.
 * The pointer mu in the declarations below should point to the diagonal
 * element a_0.
 * So the matrix above can be defined as follows
 * 
 * double *A = malloc(sizeof(double)*(2*n-1));
 * double *mu = A+n;
 * for (i=1-n;i<n;i++)
 *   mu[i] = ... set the value a_{i}
 *
 * All the solvers, except the classical Levinson algorithm,
 * share a common function prototype.
 *
 * solver_name (fdouble *x,fdouble *mu,finteger N,fdouble *b,
 *      fdouble col_tol,fdouble tol,int k_max,int refine)
 *
 * Arguments:
 * x          This is the memory where the solution vector is stored,
 *            length should be N.
 * mu         The coefficients of the Toeplitz matrix stored in the
 *            scheme described above.
 * N          The dimension of the system of equations.
 * b          The right-hand side vector of length N.
 * col_tol    Tolerance to accept a column-regular pair.
 * tol        Tolerance to accept a regular pair.
 * k_max      The maximum block size for the look-ahead algorithms
 * refine     A flag which indicates if iterative refinement will be 
 *            done by the solver.
 *
 * The solver_name could be one of (references in README)
 * d_lev_inner          look-ahead Levinson algorithm from [2],
 *                      (using inner FBOPs)
 * d_lev                look-ahead Levinson algorithm from [1],
 *                      (product form)
 * d_schurlev_inner     look-ahead Schur algorithm from [2],
 *                      (using inner FBOPs + Levinson rec.)
 * d_schurlev           look-ahead Schur algorithm from [1],
 *                      (product form + Levinson rec.)
 * d_schur              look-ahead Schur algorithm from [1],
 *                      (product form, requires O(N^2) storage)
 * d_schurup            look-ahead Schur algorithm from [1],
 *                      (product form, only O(N) storage)
 * d_schur_inner        look-ahead Schur algorithm from [2],
 *                      (using inner FBOPs, requires O(N^2) storage)
 * d_schur_innerup      look-ahead Schur algorithm from [2],
 *                      (using inner FBOPs, only O(N) storage)
 *
 * The classical Levinson algorithm
 *
 * d_lev_cl(fdouble *x,fdouble *mu,finteger N,fdouble *b)
 *
 * where the arguments are the same as above, except that all
 * parameters needed for look-ahead are missing.
 *
 * If the refine flag is set to 1, you can apply one of the following 
 * subroutines.
 *
 * refine_solution(fdouble *x,fdouble *b)
 *
 * computes one further refinement step on the right-hand side b, the 
 * vector x is overwritten by this function with the refined solution.
 *
 * invert_toeplitz(fdouble *x,fdouble *b)
 *
 * applies the computed decomposition to another right-hand side x, the
 * vector b is overwritten by this function with the solution.
 *
 * If the refine flag is given to the solver, it allocates some static
 * memory where some information is stored to apply both functions above.
 * To free this storage call the subroutine
 *
 * free_refinement_workspace()
 *
 */

#ifdef  __cplusplus
extern "C" { 
#endif

extern int d_lev_inner(fdouble *x,fdouble *mu,finteger N,fdouble *b,
       fdouble col_tol,fdouble tol,int k_max,int refine, fdouble *psol); 

extern int d_lev(fdouble *x,fdouble *mu,finteger N,fdouble *b,
       fdouble col_tol,fdouble tol,int k_max,int refine); 

extern int d_schurlev_inner(fdouble *x,fdouble *mu,finteger N,
       fdouble *b,fdouble col_tol,fdouble tol,int k_max,int refine); 

extern int d_schurlev(fdouble *x,fdouble *mu,finteger N,
       fdouble *b,fdouble col_tol,fdouble tol,int k_max,int refine); 

extern int d_schur(fdouble *x,fdouble *mu,finteger N,
       fdouble *b,fdouble col_tol,fdouble tol,int k_max,int refine); 

extern int d_schurup(fdouble *x,fdouble *mu,finteger N,
       fdouble *b,fdouble col_tol,fdouble tol,int k_max,int refine); 

extern int d_schur_inner(fdouble *x,fdouble *mu,finteger N,
       fdouble *b,fdouble col_tol,fdouble tol,int k_max,int refine); 

extern int d_schur_innerup(fdouble *x,fdouble *mu,finteger N,
       fdouble *b,fdouble col_tol,fdouble tol,int k_max,int refine); 

extern int d_lev_cl(fdouble *x,fdouble *mu,finteger N,fdouble *b);

extern void refine_solution(fdouble *x,fdouble *b);

extern void invert_toeplitz(fdouble *x,fdouble *b);

extern void free_refinement_workspace();

#ifdef  __cplusplus
}
#endif

typedef struct
{
  fdouble *data;
  finteger length;
  finteger stdinc;
} polynom;

typedef struct
{
  fdouble *data;
  finteger N;
} toeplitz;

/* usefull stuff */
#ifndef max
#define max(a,b) (((a)<(b))?(b):(a))
#endif

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif

#define liwork(much) ((sizeof(finteger)<=sizeof(fdouble)) ? much	\
      : (((sizeof(finteger)%sizeof(fdouble)) == 0) ?			\
	    (sizeof(finteger)/sizeof(fdouble)*much) :			\
	    ((much*sizeof(finteger))/sizeof(fdouble)+1)))

#endif
