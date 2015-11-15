#include "clik.h"
#include "clik_helper.h"
#include <errno.h>
#include <string.h>


void actspt_extra_free_();
void actspt_extra_lkl_(double*, double*);
void actspt_extra_only_one_(int*);
void actspt_extra_parameter_init_(char*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);

typedef struct {
  char tmpdir[800];
  } actspt;


void free_actspt(void **none) {
  actspt_extra_free_();
}

double actspt_lkl(void* none, double* pars, error **err) {
  double lkl;
  
  //_DEBUGHERE_("%g %g %g %g",pars[0],pars[1],pars[2],pars[3]);
  actspt_extra_lkl_(&lkl,pars);
  return lkl;
}

cmblkl* clik_actspt_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char directory_name[4096],pwd[4096],pwd2[4096];
  int status;
  int bok;
  cmblkl *cing;
  int mlmax;
  int ilmin11,ilmin12,ilmin22,ilmax11,ilmax12,ilmax22,itt_lmax_mc,iuse_act_south  , iuse_act_equa   , iuse_spt_lowell , iuse_spt_highell;
  char dir_data[2048];
  int ldd;
  int xdim;
  char *xnames_def[] = {"A_sz","A_ksz", "xi_sz_cib", "a_ps_act_148","a_ps_act_217","a_ps_spt_95","a_ps_spt_150","a_ps_spt_220","A_cib_143",
        "A_cib_217","n_Dl_cib","r_ps_spt_95x150","r_ps_spt_95x220","r_ps_150x220","r_cib","a_gs","a_ge","cal_acts_148","cal_acts_217","cal_acte_148","cal_acte_217","cal_spt_95","cal_spt_150","cal_spt_220"};

  actspt_extra_only_one_(&bok);
  testErrorRet(bok!=0,-100,"actspt already initialized",*err,__LINE__,NULL);
  
  // get data and change dir
  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  ilmin11 = cldf_readint(df,"lmin11",err);
  forwardError(*err,__LINE__,NULL);
  ilmin12 = cldf_readint(df,"lmin12",err);
  forwardError(*err,__LINE__,NULL);
  ilmin22 = cldf_readint(df,"lmin22",err);
  forwardError(*err,__LINE__,NULL);
  ilmax11 = cldf_readint(df,"lmax11",err);
  forwardError(*err,__LINE__,NULL);
  ilmax12 = cldf_readint(df,"lmax12",err);
  forwardError(*err,__LINE__,NULL);
  ilmax22 = cldf_readint(df,"lmax22",err);
  forwardError(*err,__LINE__,NULL);
  

  itt_lmax_mc = cldf_readint(df,"tt_lmax_mc",err);
  forwardError(*err,__LINE__,NULL);
  
  iuse_act_equa = cldf_readint(df,"use_act_equa",err);
  forwardError(*err,__LINE__,NULL);
  
  iuse_act_south = cldf_readint(df,"use_act_south",err);
  forwardError(*err,__LINE__,NULL);
  
  iuse_spt_highell = cldf_readint(df,"use_spt_highell",err);
  forwardError(*err,__LINE__,NULL);
  

  memset(dir_data,' ',sizeof(char)*2048);
  sprintf(dir_data,"data/");
  dir_data[5] = ' ';
  ldd = 5;
  
  // call actspt_init
  actspt_extra_parameter_init_(dir_data,&ldd,&ilmin11,&ilmin12,&ilmin22,&ilmax11,&ilmax12,&ilmax22,&itt_lmax_mc,&iuse_act_south  , &iuse_act_equa    , &iuse_spt_highell);

  cldf_external_cleanup(directory_name,pwd,err);  
  forwardError(*err,__LINE__,NULL);
  
    
  xdim = 24;

  cing = init_cmblkl(NULL, &actspt_lkl, 
                     &free_actspt,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,xdim,err);
  forwardError(*err,__LINE__,NULL);

  cmblkl_set_names(cing, xnames_def,err);
  forwardError(*err,__LINE__,NULL);

  return cing;
}
