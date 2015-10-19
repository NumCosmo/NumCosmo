#include "clik.h"
#include "clik_helper.h"
#include <errno.h>
#include <string.h>

void bflike_extra_lkl_(double*,double*);
void bflike_extra_free_();
void bflike_extra_parameter_init_(char*,int*,int*,int*);

void bflike_qu_extra_lkl_(double*,double*);
void bflike_qu_extra_free_();
void bflike_qu_extra_parameter_init_(char*,int*,int*,int*);

void bflike_smw_extra_lkl_(double*,double*);
void bflike_smw_extra_free_();
void bflike_smw_extra_parameter_init_(char*,int*,int*,int*);


double bflike_lkl(void* none, double* pars, error **err) {
  double lkl;
  bflike_extra_lkl_(&lkl,pars);
  
  return lkl;
}

void free_bflike(void **none) {
  
  bflike_extra_free_();
}

cmblkl* clik_bflike_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char directory_name[4096],pwd[4096],pwd2[4096];
  int status;
  int bok;
  cmblkl *cing;
  int mlmax;
  char dir_data[2048];
  int ldd;
  int lmin,lmax;
  int firstchain,lastchain,firstsample,lastsample,step;
  int hk;

  lmin = ell[0];
  lmax = ell[nell-1];
  
  // get data and change dir
  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);

  memset(dir_data,' ',sizeof(char)*2048);
  sprintf(dir_data,".");
  dir_data[1] = ' ';
  ldd = 2;
  
  //call

  bflike_extra_parameter_init_(dir_data,&ldd,&lmin,&lmax);
  
  cldf_external_cleanup(directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  cing = init_cmblkl(NULL, &bflike_lkl, 
                     &free_bflike,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);

  return cing;
}

double bflike_QU_lkl(void* none, double* pars, error **err) {
  double lkl;
  bflike_qu_extra_lkl_(&lkl,pars);
  
  return lkl;
}

void free_bflike_QU(void **none) {
  
  bflike_qu_extra_free_();
}

cmblkl* clik_bflike_QU_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char directory_name[4096],pwd[4096],pwd2[4096];
  int status;
  int bok;
  cmblkl *cing;
  int mlmax;
  char dir_data[2048];
  int ldd;
  int lmin,lmax;
  int firstchain,lastchain,firstsample,lastsample,step;
  int hk;

  lmin = ell[0];
  lmax = ell[nell-1];
  
  // get data and change dir
  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);

  memset(dir_data,' ',sizeof(char)*2048);
  sprintf(dir_data,".");
  dir_data[1] = ' ';
  ldd = 2;
  
  //call

  bflike_qu_extra_parameter_init_(dir_data,&ldd,&lmin,&lmax);
  
  cldf_external_cleanup(directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  cing = init_cmblkl(NULL, &bflike_QU_lkl, 
                     &free_bflike_QU,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);

  return cing;
}

double bflike_smw_lkl(void* none, double* pars, error **err) {
  double lkl;

  bflike_smw_extra_lkl_(&lkl,pars);
  return lkl;
}

void free_bflike_smw(void **none) {
  
  bflike_smw_extra_free_();
}

cmblkl* clik_bflike_smw_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char directory_name[4096],pwd[4096],pwd2[4096];
  int status;
  int bok;
  cmblkl *cing;
  int mlmax;
  char dir_data[2048];
  int ldd;
  int lmin,lmax;
  int firstchain,lastchain,firstsample,lastsample,step;
  int hk;

  lmin = ell[0];
  lmax = ell[nell-1];
  
  // get data and change dir
  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);

  memset(dir_data,' ',sizeof(char)*2048);
  sprintf(dir_data,".");
  dir_data[1] = ' ';
  ldd = 2;
  
  //call

  bflike_smw_extra_parameter_init_(dir_data,&ldd,&lmin,&lmax);
  
  cldf_external_cleanup(directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  cing = init_cmblkl(NULL, &bflike_smw_lkl, 
                     &free_bflike_smw,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);

  return cing;
}

