#include "clik.h"
#include "clik_helper.h"
#include <errno.h>
#include <string.h>


void gibbs_extra_free_(int*);
void gibbs_extra_lkl_(double*,int*,double*);
void   gibbs_extra_parameter_init_(int*, char*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
void gibbs_gauss_extra_free_(int*);
void gibbs_gauss_extra_lkl_(double*,int*,double*);
void gibbs_gauss_extra_parameter_init_(int*,int*,int*,int*);
void comm_lowl_extra_free_(int*);
void comm_lowl_extra_lkl_(double*,int*,double*);
void comm_lowl_extra_parameter_init_(int*,char*,int*,int*,int*);


typedef struct {
  int handle;
  int handle_transition;
  int ltrans,lmin;
  } gibbs;


void free_gibbs(void **none) {
  gibbs *gb;

  
  gb = *none;
  gibbs_extra_free_(&(gb->handle));
  if (gb->handle_transition!=-1) {
    gibbs_extra_free_(&(gb->handle_transition));
  }
}

double gibbs_lkl(void* none, double* pars, error **err) {
  double lkl;
  gibbs *gb;

  gb = none;
  gibbs_extra_lkl_(&lkl,&gb->handle,pars);
  
  if (gb->handle_transition!=-1) {
    double lkl_trans;
    gibbs_extra_lkl_(&lkl_trans,&gb->handle_transition,&(pars[gb->ltrans-gb->lmin]));
    
    lkl -= lkl_trans;
  }

  return lkl;
}

cmblkl* clik_gibbs_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char directory_name[4096],pwd[4096],pwd2[4096];
  int status;
  int bok;
  cmblkl *cing;
  int mlmax;
  char dir_data[2048];
  int ldd;
  int lmin,lmax;
  int firstchain,lastchain,firstsample,lastsample,step,approx_chi2;
  int hk;
  gibbs *gb;

  lmin = ell[0];
  lmax = ell[nell-1];
  
  // get data and change dir
  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  firstchain = cldf_readint(df,"firstchain",err);
  forwardError(*err,__LINE__,NULL);
  lastchain = cldf_readint(df,"lastchain",err);
  forwardError(*err,__LINE__,NULL);
  firstsample = cldf_readint(df,"firstsample",err);
  forwardError(*err,__LINE__,NULL);
  lastsample = cldf_readint(df,"lastsample",err);
  forwardError(*err,__LINE__,NULL);
  step = cldf_readint(df,"step",err);
  forwardError(*err,__LINE__,NULL);

  approx_chi2 = 0;
  hk = cldf_haskey(df,"approx_chi2",err);
  forwardError(*err,__LINE__,NULL);
  if (hk == 1) {
    approx_chi2 = cldf_readint(df,"approx_chi2",err);
    forwardError(*err,__LINE__,NULL);    
  }
  

  memset(dir_data,' ',sizeof(char)*2048);
  sprintf(dir_data,"data/");
  dir_data[5] = ' ';
  ldd = 5;
  
  gb = malloc_err(sizeof(gibbs),err);
  forwardError(*err,__LINE__,NULL);
  
  gb->handle = 0;
  gb->handle_transition = -1;


  //call
  gibbs_extra_parameter_init_(&(gb->handle),dir_data,&ldd,&lmin,&lmax,&firstchain,&lastchain,&firstsample,&lastsample,&step,&approx_chi2);
  testErrorRetVA(gb->handle<=0,-43255432,"handle return is negative (got %d)",*err,__LINE__,NULL,gb->handle);

  hk = cldf_haskey(df,"ltrans",err);
  forwardError(*err,__LINE__,NULL);
  if (hk == 1) {
    int ltrans;
  
    gb->handle_transition = 1;
    ltrans = cldf_readint(df,"ltrans",err);
    forwardError(*err,__LINE__,NULL);
    gb->lmin = lmin;
    gb->ltrans = ltrans;
    
    gibbs_extra_parameter_init_(&(gb->handle_transition),dir_data,&ldd,&ltrans,&lmax,&firstchain,&lastchain,&firstsample,&lastsample,&step,&approx_chi2);
    testErrorRetVA(gb->handle<=0,-43255432,"handle return is negative (got %d)",*err,__LINE__,NULL,gb->handle);

  }

  cldf_external_cleanup(directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  cing = init_cmblkl(gb, &gibbs_lkl, 
                     &free_gibbs,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);

  return cing;
}



void free_gauss_gibbs(void **none) {
  gibbs *gb;

  
  gb = *none;
  gibbs_gauss_extra_free_(&(gb->handle));
  free(gb);
}

double gibbs_gauss_lkl(void* none, double* pars, error **err) {
  double lkl;
  gibbs *gb;

  gb = none;
  gibbs_gauss_extra_lkl_(&lkl,&gb->handle,pars);
  
  return lkl;
}


cmblkl* clik_gibbs_gauss_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char directory_name[4096],pwd[4096],pwd2[4096];
  int status;
  int bok;
  cmblkl *cing;
  int mlmax;
  char dir_data[2048];
  int ldd;
  int lmin,lmax;
  int delta_l;
  int hk;
  gibbs *gb;

  lmin = ell[0];
  lmax = ell[nell-1];
  
  // get data and change dir
  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  delta_l = cldf_readint(df,"delta_l",err);
  forwardError(*err,__LINE__,NULL);
  
  
  gb = malloc_err(sizeof(gibbs),err);
  forwardError(*err,__LINE__,NULL);
  
  gb->handle = 0;
  
  //call
  gibbs_gauss_extra_parameter_init_(&(gb->handle),&lmin,&lmax,&delta_l);
  testErrorRetVA(gb->handle<=0,-43255432,"handle return is negative (got %d)",*err,__LINE__,NULL,gb->handle);

  cldf_external_cleanup(directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  cing = init_cmblkl(gb, &gibbs_gauss_lkl, 
                     &free_gauss_gibbs,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);

  return cing;
}




void free_comm_lowl(void **phandle) {
  
  comm_lowl_extra_free_(*phandle);
}

double comm_lowl_lkl(void* handle, double* pars, error **err) {
  double lkl;
  
  comm_lowl_extra_lkl_(&lkl,handle,pars);
  
  return lkl;
}

cmblkl* clik_comm_lowl_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char directory_name[4096],pwd[4096],pwd2[4096];
  int status;
  int bok;
  cmblkl *cing;
  int mlmax;
  char parfile[2048];
  int ldd;
  int lmin,lmax;
  int firstchain,lastchain,firstsample,lastsample,step;
  int hk;
  int *gb;

  lmin = ell[0];
  lmax = ell[nell-1];
  
  // get data and change dir
  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);

  memset(parfile,' ',sizeof(char)*2048);
  sprintf(parfile,"comm_lowl.par");
  parfile[13] = ' ';
  ldd = 13;
  
  gb = malloc_err(sizeof(int),err);
  forwardError(*err,__LINE__,NULL);
  
  *((int*) gb) =0;
  //call
  comm_lowl_extra_parameter_init_(gb,parfile,&ldd,&lmin,&lmax);
  testErrorRetVA(*((int*) gb) <=0,-43255432,"handle return is negative (got %d)",*err,__LINE__,NULL,*gb);

  cldf_external_cleanup(directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  cing = init_cmblkl(gb, &comm_lowl_lkl, 
                     &free_comm_lowl,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);

  return cing;
}
