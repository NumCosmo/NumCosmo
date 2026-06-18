/*
 *  lowly_common.h
 *  lowly_project
 *
 *  Created by Karim Benabed on 16/04/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include "pmc.h"
//#include "erfinv.h"


#ifdef HAS_MKL
#include "mkl_lapack.h"
#elif HL2_ACML
#include <acml.h>
#include <acml_mv.h>
#elif LAPACK_CLIK
#include "lapack_clik.h"
#else
#include "clapack.h"
#endif

#ifndef __LOWLY_COM__
#define __LOWLY_COM__

#define PIXEL_SEEN 1
#define PIXEL_UNSEEN 0

#define HP_RING -1
#define HP_NEST  1

#define LOWLY_PARANOID_NAN_TEST

#ifdef HL2_ACML
#define dsymv(...) dsymv_(__VA_ARGS__,1)
#define dpotrf(...) dpotrf_(__VA_ARGS__,1)
#define dpotri(...) dpotri_(__VA_ARGS__,1)
#define dtsrv(...) dtsrv_(__VA_ARGS__,1,1,1)
#define dgemv(...) dgemv_(__VA_ARGS__,1)
#define dsyrk(...) dsyrk_(__VA_ARGS__,1,1)
#define dtrtri(...) dtrtri_(__VA_ARGS__,1,1)
#define dtrmm(...) dtrmm_(__VA_ARGS__,1,1,1,1)
#define dgeqrf(...) dgeqrf_(__VA_ARGS__)
#define dtrmv(...) dtrmv_(__VA_ARGS__,1,1,1)
#define dormqr(...) dormqr_(__VA_ARGS__,1,1)
#define dsyr2k(...) dsyr2k_(__VA_ARGS__,1,1)
#define dgemm(...) dgemm_(__VA_ARGS__,1,1)

#endif
long* lowly_build_pixel_list(unsigned char *Mask,long *pnpix, error **err);

double * lowly_get_posvec(const long nside,
                    const long * pixel_indices,
                    const long npix_seen, 
                    const int ordering,error **err);
void * lowly_ring_reorder2D(void* from,int nside, size_t sz,error **err);
void * lowly_ring_reorder(void* from,int nside, size_t sz,error **err);
double lowly_computeNl(double noisevar, long nside);
int lowly_which_order(char *ordering,error **err);
void lowly_print_stat(FILE*,int doprint, int n_tests,int time_build, int time_chol,int time_tot);
double lowly_XtRX_lkl(double *S,double *X, double* buffer,long npix_tot,error **err);
double lowly_XtRX_lkl_K(double *S,double *X, double* buffer,long npix_tot,double* Ktilde,long npix_other, error **err);

int lowly_get_offset_cl(int *has_cl, int* offset_cl,int nell);
int lowly_nmode(int nell,int *ell);
int *lowly_get_ell(int *pnell, int *inell, int lmax,error **err);

void lowly_dump(char *filename, void* ptr, size_t sz, error **err);


#define TIMER_DEFS struct timeval TIMER_before, TIMER_after;
#define TIMER_IN gettimeofday(&TIMER_before,NULL);
#define TIMER_OUT gettimeofday(&TIMER_after,NULL);
#define TIMER_MSEC (((TIMER_after.tv_sec-TIMER_before.tv_sec)*1000000+(TIMER_after.tv_usec-TIMER_before.tv_usec))/1000)
#define STRUCT_TIMER_DEFS long time_build,time_chol,time_tot,n_tests; int print_stat_flag
#define DO_PRINT_STAT(self) lowly_print_stat(stderr,self->print_stat_flag,self->n_tests,self->time_build,self->time_chol,self->time_tot);
#define SET_PRINT_STAT(self) self->print_stat_flag=1

#define _SZT_(val) ((size_t) val)

#define clowly_base           -9000
#define lowly_chol           -1  + clowly_base
#define lowly_unkorder       -4  + clowly_base
#define lowly_lbig           -12 + clowly_base

#endif
