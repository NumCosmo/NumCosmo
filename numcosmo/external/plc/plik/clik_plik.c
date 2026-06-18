#include "clik.h"
#include "smica.h"
#include "clik_helper.h"

int clik_is_plik(clik_object* clikid,error **_err) {
  lklbs *lbs;
  int i;
  _dealwitherr;

  lbs = _clik_dig(clikid,err);
  _forwardError(*err,__LINE__,0);
  for(i=0;i<lbs->nlkl;i++) {
    if(lbs->lkls[i]->lkl_func==&Smica_lkl) {
      //_DEBUGHERE_("%p",lbs->lkls[i]->lkl_data)
      return i;
    }  
  }
  return -1;
}

int clik_must_be_plik(clik_object *clikid, error **_err) {
  int r;
  _dealwitherr;
  r = clik_is_plik(clikid,err);
  _forwardError(*err,__LINE__,0);
  _testErrorRetVA(r==-1,-13141,"clik is not plik !!!",*err,__LINE__,-1,NULL);
  return r;
}

void* clik_must_be_plik_ptr(clik_object *clikid, error **_err) {
  int ir;
  lklbs *lbs;
  _dealwitherr;

  ir = clik_must_be_plik(clikid,err);
  _forwardError(*err,__LINE__,0);
  lbs = _clik_dig(clikid,err);
  _forwardError(*err,__LINE__,0);
  return lbs->lkls[ir]->lkl_data;
}

void lklbs_bs_compute(lklbs *self,double* pars, error **err);
double * lklbs_get_cls(lklbs *self,int ilkl, double *pars, error **err);

#ifdef ADD0US
void fortran_plik_get_fg(long* pself, double* cl_and_pars, double* vec) {
#elseif ADD2US
void fortran_plik_get_fg__(long* pself, double* cl_and_pars, double* vec) {
#else
void fortran_plik_get_fg_(long* pself, double* cl_and_pars, double* vec) {
#endif
  clik_object* self;
  self = (clik_object*) *pself;
  void *smic;
  int ir;
  lklbs *lbs;
  double *cls;
  
  ir = clik_must_be_plik(self,NULL);  
  error *_err;
  error **err;
  _err = NULL;
  err = &_err;
  
  lbs = _clik_dig(self,err);
  quitOnError(*err,__LINE__,stderr);
  
  lklbs_bs_compute(lbs,cl_and_pars, err);
  quitOnError(*err,__LINE__,stderr);
  
  cls = lklbs_get_cls(lbs,ir,cl_and_pars,err);
  quitOnError(*err,__LINE__,stderr);
  
  smic = lbs->lkls[ir]->lkl_data;
  
  Smica_fg(smic,cls,vec,err);
  quitOnError(*err,__LINE__,stderr);
  
}

#ifdef ADD0US
void fortran_plik_get_cal_beam(long* pself, double* cl_and_pars, double* vec) {
#elseif ADD2US
void fortran_plik_get_cal_beam__(long* pself, double* cl_and_pars, double* vec) {
#else
void fortran_plik_get_cal_beam_(long* pself, double* cl_and_pars, double* vec) {
#endif
  clik_object* self;
  self = (clik_object*) *pself;
  void *smic;
  int ir;
  lklbs *lbs;
  double *cls;
  
  ir = clik_must_be_plik(self,NULL);  
  error *_err;
  error **err;
  _err = NULL;
  err = &_err;
  
  lbs = _clik_dig(self,err);
  quitOnError(*err,__LINE__,stderr);
  
  lklbs_bs_compute(lbs,cl_and_pars, err);
  quitOnError(*err,__LINE__,stderr);
 
  cls = lklbs_get_cls(lbs,ir,cl_and_pars,err);
  quitOnError(*err,__LINE__,stderr);
  
  smic = lbs->lkls[ir]->lkl_data;
  
  Smica_gcal(smic,cls,vec,err);
  quitOnError(*err,__LINE__,stderr);
  
}

#ifdef ADD0US
void fortran_plik_get_data(long* pself,  double* vec) {
#elseif ADD2US
void fortran_plik_get_data__(long* pself,  double* vec) {
#else
void fortran_plik_get_data_(long* pself,  double* vec) {
#endif
  clik_object* self;
  self = (clik_object*) *pself;
  void *smic;
  
  smic = clik_must_be_plik_ptr(self,NULL);  
  error *_err;
  error **err;
  _err = NULL;
  err = &_err;
  Smica_data(smic,vec,err);
  quitOnError(*err,__LINE__,stderr);
  
}

#ifdef ADD0US
void fortran_plik_get_vecsize(long* pself,  int* vecs) {
#elseif ADD2US
void fortran_plik_get_vecsize__(long* pself,  int* vecs) {
#else
void fortran_plik_get_vecsize_(long* pself,  int* vecs) {
#endif
  clik_object* self;
  self = (clik_object*) *pself;
  void *smic;

  smic = clik_must_be_plik_ptr(self,NULL);  
  error *_err;
  error **err;
  _err = NULL;
  err = &_err;
  *vecs = Smica_vecsize(smic,err);
  quitOnError(*err,__LINE__,stderr);
  
}
