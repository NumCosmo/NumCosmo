#include "errorlist.h"
#ifndef CLIK_LENSING_
#define CLIK_LENSING_

typedef char parnam[256];

typedef struct {
  void *plens_payload;
  int lmax[7];
  int type;
  int renorm;
  int ren1;
  double check;
  int has_check;
  double *cl_fid;
} clik_lensing_object;

clik_lensing_object* clik_lensing_init(char *fpath, error **err);
double clik_lensing_compute(clik_lensing_object *lclik, double *pars,error **err);
void clik_lensing_cleanup(clik_lensing_object **plclik);
int clik_try_lensing(char *fpath,error **err);
int clik_lensing_get_lmax(clik_lensing_object *lclik, error **err);
void clik_lensing_get_lmaxs(clik_lensing_object *lclik, int *lmax, error **err);
int clik_lensing_get_extra_parameter_names(clik_lensing_object* lclik, parnam **names, error **err);
double* clik_lensing_clcmb_fid(clik_lensing_object* lclik, error **err);
double* clik_lensing_cltt_fid(clik_lensing_object* lclik, error **err);
double* clik_lensing_clpp_fid(clik_lensing_object* lclik, error **_err);
void clik_lensing_selftest(clik_lensing_object *lclik, char *fpath, error **err);

#endif