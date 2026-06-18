#include "clik.h"
#include "clik_helper.h"
#include <errno.h>
#include <string.h>

void lowlike_extra_free_();
void lowlike_extra_lkl_(double*,double*);
void lowlike_extra_only_one_(int*);
void lowlike_extra_parameter_init_(int*,int*,int*,int*,int*,int*,int*);


typedef struct {
  char tmpdir[800];
  } lowlike;


void free_lowlike(void **none) {
  lowlike_extra_free_();
}

double lowlike_lkl(void* none, double* pars, error **err) {
  double lkl;
  
  lowlike_extra_lkl_(&lkl,pars);
  return lkl;
}

cmblkl* clik_lowlike_init(cldf * df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char directory_name[4096],pwd[4096];
  int status;
  int bok;
  cmblkl *cing;
  int mlmax;
  int ttmin,ttmax,temin,temax,use_gibbs,use_lowl_pol,use_wmap_pol;
  int hk;
  
  

  lowlike_extra_only_one_(&bok);
  testErrorRet(bok!=0,-100,"lowlike already initialized",*err,__LINE__,NULL);
  
  // get data and change dir
  //clik_external_data_init(directory_name,pwd,group_id,cur_lkl,err);
  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  ttmin = cldf_readint(df,"ttmin",err);
  forwardError(*err,__LINE__,NULL);
  ttmax = cldf_readint(df,"ttmax",err);
  forwardError(*err,__LINE__,NULL);
  temin = cldf_readint(df,"temin",err);
  forwardError(*err,__LINE__,NULL);
  temax = cldf_readint(df,"temax",err);
  forwardError(*err,__LINE__,NULL);
  use_gibbs = cldf_readint(df,"use_gibbs",err);
  forwardError(*err,__LINE__,NULL);
  use_lowl_pol = cldf_readint(df,"use_lowl_pol",err);
  forwardError(*err,__LINE__,NULL);
  
  //hstat = H5LTget_attribute_int( group_id, ".", "ttmin",  &ttmin);
  //testErrorRetVA(hstat<0,hdf5_base,"cannot read ttmin in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  //hstat = H5LTget_attribute_int( group_id, ".", "ttmax",  &ttmax);
  //testErrorRetVA(hstat<0,hdf5_base,"cannot read ttmax in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  //hstat = H5LTget_attribute_int( group_id, ".", "temin",  &temin);
  //testErrorRetVA(hstat<0,hdf5_base,"cannot read temin in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  //hstat = H5LTget_attribute_int( group_id, ".", "temax",  &temax);
  //testErrorRetVA(hstat<0,hdf5_base,"cannot read temax in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);  
  //hstat = H5LTget_attribute_int( group_id, ".", "use_gibbs",  &use_gibbs);
  //testErrorRetVA(hstat<0,hdf5_base,"cannot read use_gibbs in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  //hstat = H5LTget_attribute_int( group_id, ".", "use_lowl_pol",  &use_lowl_pol);
  //testErrorRetVA(hstat<0,hdf5_base,"cannot read use_lowl_pol in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  use_wmap_pol = 0;
  use_wmap_pol = cldf_readint_default(df,"use_wmap_pol",0,err);
  forwardError(*err,__LINE__,NULL);
  
  //hstat = H5LTfind_attribute(group_id, "use_wmap_pol");
  //if (hstat==1) {
  //  hstat = H5LTget_attribute_int( group_id, ".", "use_wmap_pol",  &use_wmap_pol);
  //  testErrorRetVA(hstat<0,hdf5_base,"cannot read use_wmap_pol in %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);
  //} 
  mlmax = ell[nell-1];
  testErrorRet(mlmax<ttmax || mlmax<temax,-101011,"Bad parameters for lowlike likelihood",*err,__LINE__,NULL);

  
  // call lowlike_init
  lowlike_extra_parameter_init_(&ttmin,&ttmax,&temin,&temax,&use_gibbs,&use_lowl_pol,&use_wmap_pol);
  
  //clik_external_data_cleanup(directory_name,pwd,err);  
  cldf_external_cleanup(directory_name,pwd,err);  
  forwardError(*err,__LINE__,NULL);
  
  
  
  cing = init_cmblkl(NULL, &lowlike_lkl, 
                     &free_lowlike,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  forwardError(*err,__LINE__,NULL);
  return cing;
}
