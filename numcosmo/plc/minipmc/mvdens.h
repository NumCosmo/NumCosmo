/*
 *  mvdens.h
 *  likely
 *
 *  Created by Karim Benabed on 10/03/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __MVDENS_H
#define __MVDENS_H
#ifdef __PLANCK__

#include "HL2_likely/tools/errorlist.h"
#include "HL2_likely/tools/io.h"
#include "HL2_likely/tools/maths_base.h"

#else

#include "errorlist.h"
#include "io.h"
#include "maths_base.h"

#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_int.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>

/*#ifdef HAS_LAPACK
#ifdef HAS_MKL
#include "mkl_blas.h"
#include "mkl_lapack.h"
#else
#include "blas.h"
#include "lapack.h"
#endif
#endif
*/

#define __MVDENS_PARANOID_DEBUG__

/* errors */
#define mv_base       -700
#define mv_allocate   -1 + mv_base
#define mv_serialize  -2 + mv_base
#define mv_outOfBound -3 + mv_base
#define mv_badComm    -4 + mv_base
#define mv_negWeight  -5 + mv_base
#define mv_cholesky   -6 + mv_base
#define mv_negative   -7 + mv_base
#define mv_undef      -8 + mv_base
#define mv_file       -9 + mv_base
#define mv_io        -10 + mv_base
#define mv_tooManySteps -11 + mv_base
#define mv_dimension -12 + mv_base
#define mv_type      -13 + mv_base
#define mv_negHatCl  -14 + mv_base

#define PI         3.141592653589793
#define LOGSQRT2PI 0.918938533204673
#ifndef DYNMAX
#define DYNMAX 500.0
#endif
#define mv_INF 1e99

#define mv_EXP  2
#define mv_NORM 1
#define mv_UNORM 0

#define MVT_DF 3

typedef struct {
  size_t ndim;
  void *buf;
  int own_buf;
  double *mean;
  double *std;
  double *x_tmp;
  gsl_vector_view  mean_view_container;
  gsl_vector *mean_view;
  gsl_vector_view  x_tmp_view_container;
  gsl_vector *x_tmp_view;
  gsl_matrix_view  std_view_container;
  gsl_matrix  *std_view;
  size_t band_limit;
  int chol;
  int df;
  double detL;   /* Determinant of L, where L^t L = Covariance matrix */
} mvdens;

typedef struct {
  size_t ncomp, ndim;
  int init_cwght;
  double *wght, *cwght;
  gsl_vector_view wght_view_container, cwght_view_container;
  gsl_vector *wght_view, *cwght_view;
  mvdens **comp;
  void *buf_comp;
} mix_mvdens;

// init and stuff
mvdens * mvdens_alloc(size_t ndim, error **err);
mvdens * mvdens_init(mvdens *g, size_t nndim, double* mn, double* std, error **err);
void mvdens_empty(mvdens *g);
void mvdens_free(mvdens ** g);
void mvdens_free_void(void **g);
void mvdens_print(FILE* where, mvdens* what);
mvdens *mvdens_dwnp(FILE* where, error **err);
mvdens *mvdens_read_and_chol(FILE *where, error **err);
void mvdens_dump(FILE* where, mvdens* what);
void mvdens_chdump(const char *name, mvdens* what, error **err);
void mvdens_set_band_limit(mvdens *self, size_t value);
void mvdens_from_meanvar(mvdens *m, const double *pmean, const double *pvar, double correct);

//computations
double determinant(const double *L, size_t ndim);
void mvdens_cholesky_decomp(mvdens* self, error **err);
double* mvdens_ran(double* dest, mvdens * g, gsl_rng * r,error **err);
double scalar_product(mvdens *g, const double *x, error **err);
double mvdens_log_pdf(mvdens *g, const double * x, error ** err);
double mvdens_log_pdf_void(void *g, const double *x, error **err);
double mvdens_inverse(mvdens *m, error **err);

/* mix_mvdens functions */
mix_mvdens *mix_mvdens_alloc(size_t ncomp, size_t size,error **err);
void mix_mvdens_copy(mix_mvdens *target, const mix_mvdens *source, error **err);
void mix_mvdens_free_void(void **m);
void mix_mvdens_free(mix_mvdens **m);
void mix_mvdens_print(FILE* where,mix_mvdens* what);
void mix_mvdens_dump(FILE* where,mix_mvdens* what);
mix_mvdens * mix_mvdens_dwnp(FILE* where,error **err);
void mix_mvdens_set_band_limit(mix_mvdens *self, size_t value);
double effective_number_of_components(const mix_mvdens *self, error **err);

size_t mix_mvdens_size(size_t ncomp,size_t ndim);
void* serialize_mix_mvdens(mix_mvdens* self, size_t *sz,error** err);
mix_mvdens* deserialize_mix_mvdens(void* serialized, size_t sz,error **err);

void mix_mvdens_cholesky_decomp(mix_mvdens* self, error **err);
double* mix_mvdens_ran(double* dest,size_t *index, mix_mvdens *m, gsl_rng * r, error **err);
size_t ranindx(double *cw, size_t nE, gsl_rng * r, error ** err);
double mix_mvdens_log_pdf(mix_mvdens *m, const double *x, error **err);
double mix_mvdens_log_pdf_void(void *m, const double *x, error **err);

#define MALLOC_IF_NEEDED(ret,dest,size,err) { \
  if (dest==NULL) { \
    ret= (double*) malloc_err(size,err); \
  } else { \
    ret=dest; \
  } \
}

#ifdef HAS_HDF5
#include "hdf5.h"
#define hdf5_base mv_io
void mvdens_hdfdump_infile(mvdens* mv,hid_t loc_id,error **err);
void mvdens_hdfdump(mvdens* mv,char *fname,error **err);
void mix_mvdens_hdfdump(mix_mvdens* mmv,char *fname,error **err);
mvdens* mvdens_hdfdwnp_infile(mvdens* mv,hid_t loc_id,error **err);
mvdens* mvdens_hdfdwnp(char *fname,error **err);
mix_mvdens* mix_mvdens_hdfdwnp(char *fname,error **err);
#endif

#endif

