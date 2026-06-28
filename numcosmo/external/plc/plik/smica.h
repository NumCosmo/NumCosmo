/*
 *  smica.h
 *  lowly_project
 *
 *  Created by Karim Benabed on 30/10/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */


// some includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "pmc.h"
#include "lowly_common.h"

#ifndef _SMICA_
#define _SMICA_

/* typedef struct {
  int nq; // number of bins
  int nell; // lmax+1
  int m;  // number of channels
  double *wq  ; // number of modes
  double *bin ; // bins
  double *rq_hat; // empirical covariance matrices

  double *rq; // if not null contains sum_c R_q^c

  double *A;  // if not null contains mixing matrix 
  double *P;  // if not null contains power spectra
  int    *C;  // if not null contains component dimensions
  int    nc;  // number of component 
  

} smica;

smica* smica_init(int nq,int nell,int m, double *wq,double *bin,double *rq_hat, double *rq, double *A, double *P, int *C,int nc, error **err);

double smica_lkl(void* smica, double *pars, error **err);

void free_smica(smica** psmica);

//void compute_rq(smica* smic, double* cl,double *nuis,double *rq, error **err);
void bin_cl(smica* smic, double* cl , error **err);
void compute_rq(smica* smic, error **err);
double compute_smica_lkl(smica* smic, error **err);

*/
// Newer better interface

typedef void (update_rq)(void* data,double* locpars, double* rq, error **err);
typedef char _smicanames[256];
typedef double (smica_crit)(void* smic, error **err);

typedef struct {
  int ndim,m,nq;
  void* data;
  update_rq *update;
  posterior_log_free *free;
  _smicanames *names;
  _smicanames comp_name;
  int isfg;
  int ismul;
} SmicaComp;

typedef struct {
  int nq; 
  double *wq;
  int m;
  double *rq_hat, *rq_0, *rq;
  double *z_buf; // buffer
  int nc;
  SmicaComp **SC;
  int *offset_nc;
  smica_crit *crit;
  int crit_classic_init;
  double *gvec;
  double *crit_cor;
  double *eig_buf;
  int eig_lwork;
  double* eig_nrm;
  int *quad_mask;
  int quad_sn;
  void * lkl_data;
  posterior_log_free *lkl_data_free;
  int cnt;
} Smica;

Smica* Smica_init(int nq, double *wq, int m, double *rq_hat, double* rq_0, int nc, SmicaComp **SC,error **err);
double smica_crit_gauss(void *vsmic, error **err);

void smica_set_crit_gauss(Smica *smic, double *crit_cor, int *mask, int *ordering,error **err);


double Smica_lkl(void* smic, double* pars, error **err);

void free_Smica(void **smic);


// components
SmicaComp* alloc_SC(int ndim,int nq,int m,void* data, update_rq* update, posterior_log_free* pfree, error **err);
void SC_setnames(SmicaComp *Sc, char** names, error **err);
void SC_set_compname(SmicaComp *SC, char *name);
void SC_isfg(SmicaComp *SC);
void SC_ismul(SmicaComp *SC);


typedef struct {
  int Acst;
  double *AAt;
} SC_1D_data;
void comp_1D_AAt(int m, double *A, double *AAt, error **err);
SmicaComp* comp_1D_init(int q, int m, double *A, error **err);
void free_comp_1D(void** data);
void comp_1D_update(void* data,double* locpars, double* rq, error **err);

typedef struct {
  int nd,Acst;
  double *A,*Ab,*P;
} SC_nD_data;
SmicaComp* comp_nD_init(int q, int m, int nd, double *A, error **err);
void free_comp_nD(void** data);
void comp_nD_update(void* data,double* locpars, double* rq, error **err);

typedef struct {
  double* locpars;
  SmicaComp *SCnD;
  int has_cl[6];
  int jmp_cl[6];
  int trois;
} SC_CMB_data;

SmicaComp * comp_CMB_init(int nbins, int mt,int mp, int *has_cl, double* Acprs, error **err);
void free_comp_CMB(void** data);
void comp_CMB_update(void* data,double* locpars, double* rq, error **err);

void printMat(double* A, int n, int m);


#define smica_uncomp           -102 + clowly_base



typedef struct {
  double *calvec,*w;
  int *im,*other;
  int npar,mT,mP;
  int TEB[3];
} SC_calTP;

typedef struct {
  double *modes,*pars;
  int *im;
  int neigen;
} SC_beamTP;

SmicaComp* comp_beamTP_init(int q, int mT, int mP, int *TEB, int npar, int *im,int neigen, double *modes,error **err );

SmicaComp* comp_calTP_init(int q,int mT, int mP,  int *TEB, int npar, int *im,double*w,int*other, error **err );
SmicaComp* comp_icalTP_init(int q,int mT, int mP,  int *TEB, int npar, int *im,double*w,int*other, error **err );
void comp_calTP_update(void* data,double* locpars, double* rq, error **err);
void comp_icalTP_update(void* data,double* locpars, double* rq, error **err);
void comp_calTP_free(void** data);

void comp_totcal_update(void* data,double* locpars, double* rq, error **err);
SmicaComp* comp_totcal_init(int q, int mT, int mP, int *TEB,error **err );
void comp_totcalP_update(void* data,double* locpars, double* rq, error **err);
SmicaComp* comp_totcalP_init(int q, int mT, int mP, int *TEB,error **err );
SmicaComp* comp_totcalTP_init(int q, int mT, int mP, int *TEB,error **err );
SmicaComp* comp_totcalPP_init(int q, int mT, int mP, int *TEB,error **err );

void Smica_data(void* vsmic, double* fgvec, error **err);
void Smica_fg(void* vsmic, double* pars, double* fgvec, error **err);
void Smica_gcal(void* vsmic, double* pars, double* fgvec, error **err);
int Smica_vecsize(void* vsmic, error **err);

#endif
