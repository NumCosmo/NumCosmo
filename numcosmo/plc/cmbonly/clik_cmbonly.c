#include "clik.h"
#include "clik_helper.h"
#include <errno.h>
#include <string.h>


void plik_cmbonly_extra_free_();
void plik_cmbonly_extra_lkl_(double*,double*);
void plik_cmbonly_extra_only_one_(int*);
void plik_cmbonly_extra_init_(char*,int*,int*,int*,int*);

typedef struct {
  char tmpdir[800];
  } plik_cmbonly;


void free_plik_cmbonly(void **none) {
  plik_cmbonly_extra_free_();
}

double plik_cmbonly_lkl(void* none, double* pars, error **err) {
  double lkl;
  
  plik_cmbonly_extra_lkl_(&lkl,pars);
  return lkl;
}

cmblkl* clik_plik_cmbonly_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char directory_name[4096],pwd[4096],pwd2[4096];
  int status;
  int bok;
  cmblkl *cing;
  int mlmax;
  char dir_data[2048];
  int ldd,hk;
  char *xnames_def[] = {"A_Planck"};
  int use_tt, use_te, use_ee;
  int version;

  hk = cldf_haskey(df,"cmbonly_version",err);
  forwardError(*err,__LINE__,NULL);
  testErrorRet(hk==0,-132,"cmbonly v1 not supported anymore",*err, __LINE__,NULL);
  version = cldf_readint(df,"cmbonly_version",err);
  
  forwardError(*err,__LINE__,NULL);
  testErrorRet(version!=18,-132,"cmbonly plik <v17 not supported anymore",*err, __LINE__,NULL);
  plik_cmbonly_extra_only_one_(&bok);
  testErrorRet(bok!=0,-100,"plik_cmbonly already initialized",*err,__LINE__,NULL);
  // get data and change dir
  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  memset(dir_data,' ',sizeof(char)*2048);
  //sprintf(dir_data,"");
  //dir_data[5] = ' ';
  ldd = 0;
  
  // call plik_cmbonly_init
  use_tt = has_cl[0];
  use_ee = has_cl[1];
  use_te = has_cl[3];
  
  plik_cmbonly_extra_init_(dir_data,&ldd,&use_tt, &use_ee, &use_te);

  cldf_external_cleanup(directory_name,pwd,err);  
  forwardError(*err,__LINE__,NULL);
  
  cing = init_cmblkl(NULL, &plik_cmbonly_lkl, 
                     &free_plik_cmbonly,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,1,err);
  forwardError(*err,__LINE__,NULL);

  cmblkl_set_names(cing, xnames_def,err);
  forwardError(*err,__LINE__,NULL);


  return cing;
}
