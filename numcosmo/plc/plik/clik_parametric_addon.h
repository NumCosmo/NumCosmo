#ifndef _CPADN_
#define _CPADN_
#include "clik_helper.h"
#include "smica.h"

int base_parametric_cldf_init(cldf *df,int ndet, double** detlist,int *ndef, char ***defkeys, char*** defvalues, int *nvar, char ***varkeys, error **err);
SmicaComp * finalize_parametric_cldf_init(parametric* p_model,cldf *df,int nb, int m, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err);

#define CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(NAME,INIT_FUNC)                                             \
SmicaComp * clik_smica_comp_##NAME##_init(cldf * df,int nb, int mT,int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {  \
  parametric* p_model;                                                                                   \
  SmicaComp *SC;                                                                                         \
  int lmin,lmax;                                                                                         \
  double *detlist;                                                                                       \
  int ndef,nvar;                                                                                         \
  char **defkeys,**defvalues,**varkeys;                                                                  \
  double *template;                                                                                      \
  int dz,ndet;                                                                                           \
  int m;                                                                                                 \
                                                                                                         \
  m = mtot(mT,mP,has_cl);                                                                                \
                                                                                                         \
  lmin = ell[0];                                                                                         \
  lmax = ell[nell-1];                                                                                    \
  testErrorRet(nell!=(lmax-lmin+1),-111,"SAFEGARD",*err,__LINE__,NULL);                                  \
                                                                                                         \
  ndet = base_parametric_cldf_init(df,m, &detlist,&ndef, &defkeys, &defvalues, &nvar, &varkeys, err);    \
  forwardError(*err,__LINE__,NULL);                                                                      \
                                                                                                         \
  dz = -1;                                                                                               \
  template = cldf_readfloatarray(df,"template",&dz, err);                                                \
  forwardError(*err,__LINE__,NULL);                                                                          \
  p_model = INIT_FUNC(ndet, detlist, ndef, defkeys, defvalues, nvar, varkeys, lmin, lmax,template, err); \
  forwardError(*err,__LINE__,NULL);                                                                      \
                                                                                                         \
  free(detlist);                                                                                         \
  if (defkeys[0]!=NULL) {                                                                                \
    free(defkeys[0]);                                                                                    \
    free(defvalues[0]);                                                                                  \
  }                                                                                                      \
  free(defkeys);                                                                                         \
  free(defvalues);                                                                                       \
                                                                                                         \
  if (varkeys[0]!=NULL) {                                                                                \
    free(varkeys[0]);                                                                                    \
  }                                                                                                      \
  free(varkeys);                                                                                         \
                                                                                                         \
  SC = finalize_parametric_cldf_init(p_model,df,nb,m,nell,ell,has_cl,unit,wl,bins,nbins,err);            \
  forwardError(*err,__LINE__,NULL);                                                                      \
                                                                                                         \
  return SC;                                                                                             \
}

#define CREATE_PARAMETRIC_FILE_INIT(NAME,INIT_FUNC)                                                      \
SmicaComp * clik_smica_comp_##NAME##_init(cldf* df,int nb, int mT, int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {  \
  parametric* p_model;                                                                                   \
  SmicaComp *SC;                                                                                         \
  int lmin,lmax;                                                                                         \
  double *detlist;                                                                                       \
  int ndef,nvar;                                                                                         \
  char **defkeys,**defvalues,**varkeys;                                                                  \
  double *template;                                                                                      \
  int dz,ndet;                                                                                           \
  int m;                                                                                                 \
                                                                                                         \
  m = mtot(mT,mP,has_cl);                                                                                \
                                                                                                         \
  lmin = ell[0];                                                                                         \
  lmax = ell[nell-1];                                                                                    \
  testErrorRet(nell!=(lmax-lmin+1),-111,"SAFEGARD",*err,__LINE__,NULL);                                  \
                                                                                                         \
  ndet = base_parametric_cldf_init(df,m, &detlist,&ndef, &defkeys, &defvalues, &nvar, &varkeys, err);    \
  forwardError(*err,__LINE__,NULL);                                                                      \
                                                                                                         \
  p_model = INIT_FUNC(ndet, detlist, ndef, defkeys, defvalues, nvar, varkeys, lmin, lmax, err);          \
  forwardError(*err,__LINE__,NULL);                                                                      \
                                                                                                         \
  free(detlist);                                                                                         \
  if (defkeys[0]!=NULL) {                                                                                \
    free(defkeys[0]);                                                                                    \
    free(defvalues[0]);                                                                                  \
  }                                                                                                      \
  free(defkeys);                                                                                         \
  free(defvalues);                                                                                       \
                                                                                                         \
  if (varkeys[0]!=NULL) {                                                                                \
    free(varkeys[0]);                                                                                    \
  }                                                                                                      \
  free(varkeys);                                                                                         \
                                                                                                         \
  SC = finalize_parametric_cldf_init(p_model,df,nb,m,nell,ell,has_cl,unit,wl,bins,nbins,err);            \
  forwardError(*err,__LINE__,NULL);                                                                      \
                                                                                                         \
  return SC;                                                                                             \
}
#define CREATE_PARAMETRIC_POL_FILE_INIT(NAME,INIT_FUNC)                                                  \
SmicaComp * clik_smica_comp_##NAME##_init(cldf* df,int nb, int mT, int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {  \
  parametric* p_model;                                                                                   \
  SmicaComp *SC;                                                                                         \
  int lmin,lmax;                                                                                         \
  double *detlist;                                                                                       \
  int ndef,nvar;                                                                                         \
  char **defkeys,**defvalues,**varkeys;                                                                  \
  double *template;                                                                                      \
  int dz,ndet,ndet_T,ndet_P;                                                                             \
  int m;                                                                                                 \
  int *hasTEB;                                                                                           \
                                                                                                         \
                                                                                                         \
  m = mtot(mT,mP,has_cl);                                                                                \
                                                                                                         \
  lmin = ell[0];                                                                                         \
  lmax = ell[nell-1];                                                                                    \
  testErrorRet(nell!=(lmax-lmin+1),-111,"SAFEGARD",*err,__LINE__,NULL);                                  \
                                                                                                         \
  ndet = base_parametric_cldf_init(df,m, &detlist,&ndef, &defkeys, &defvalues, &nvar, &varkeys, err);    \
  forwardError(*err,__LINE__,NULL);                                                                      \
                                                                                                         \
  p_model = INIT_FUNC(mT, mP, has_cl, detlist, ndef, defkeys, defvalues, nvar, varkeys, lmin, lmax, err);\
  forwardError(*err,__LINE__,NULL);                                                                      \
                                                                                                         \
  free(detlist);                                                                                         \
  if (defkeys[0]!=NULL) {                                                                                \
    free(defkeys[0]);                                                                                    \
    free(defvalues[0]);                                                                                  \
  }                                                                                                      \
  free(defkeys);                                                                                         \
  free(defvalues);                                                                                       \
                                                                                                         \
  if (varkeys[0]!=NULL) {                                                                                \
    free(varkeys[0]);                                                                                    \
  }                                                                                                      \
  free(varkeys);                                                                                         \
                                                                                                         \
  SC = finalize_parametric_cldf_init(p_model,df,nb,m,nell,ell,has_cl,unit,wl,bins,nbins,err);            \
  forwardError(*err,__LINE__,NULL);                                                                      \
                                                                                                         \
  return SC;                                                                                             \
}

#define CREATE_PARAMETRIC_POL_TEMPLATE_FILE_INIT(NAME,INIT_FUNC)                                             \
SmicaComp * clik_smica_comp_##NAME##_init(cldf * df,int nb, int mT,int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {  \
  parametric* p_model;                                                                                   \
  SmicaComp *SC;                                                                                         \
  int lmin,lmax;                                                                                         \
  double *detlist;                                                                                       \
  int ndef,nvar;                                                                                         \
  char **defkeys,**defvalues,**varkeys;                                                                  \
  double *template;                                                                                      \
  int dz,ndet,ndet_T,ndet_P;                                                                             \
  int m;                                                                                                 \
  int *hasTEB;                                                                                           \
                                                                                                         \
  m = mtot(mT,mP,has_cl);                                                                                \
                                                                                                         \
  lmin = ell[0];                                                                                         \
  lmax = ell[nell-1];                                                                                    \
  testErrorRet(nell!=(lmax-lmin+1),-111,"SAFEGARD",*err,__LINE__,NULL);                                  \
                                                                                                         \
  ndet = base_parametric_cldf_init(df,m, &detlist,&ndef, &defkeys, &defvalues, &nvar, &varkeys, err);    \
  forwardError(*err,__LINE__,NULL);                                                                      \
                                                                                                         \
  dz = -1;                                                                                               \
  template = cldf_readfloatarray(df,"template",&dz, err);                                                \
  forwardError(*err,__LINE__,NULL);                                                                          \
  p_model = INIT_FUNC(mT, mP, has_cl, detlist, ndef, defkeys,defvalues, nvar, varkeys, lmin, lmax,template, err); \
  forwardError(*err,__LINE__,NULL);                                                                      \
                                                                                                         \
  free(detlist);                                                                                         \
  if (defkeys[0]!=NULL) {                                                                                \
    free(defkeys[0]);                                                                                    \
    free(defvalues[0]);                                                                                  \
  }                                                                                                      \
  free(defkeys);                                                                                         \
  free(defvalues);                                                                                       \
                                                                                                         \
  if (varkeys[0]!=NULL) {                                                                                \
    free(varkeys[0]);                                                                                    \
  }                                                                                                      \
  free(varkeys);                                                                                         \
                                                                                                         \
  SC = finalize_parametric_cldf_init(p_model,df,nb,m,nell,ell,has_cl,unit,wl,bins,nbins,err);            \
  forwardError(*err,__LINE__,NULL);                                                                      \
                                                                                                         \
  return SC;                                                                                             \
}

#endif

