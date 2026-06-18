//#include "clik.h"
#include "clik_helper.h"
#include <errno.h>
#include <string.h>

void camspec_extra_free_();
void camspec_extra_lkl_(double*, double*);
void camspec_extra_only_one_(int*);
void camspec_extra_fg_(double*, double*, int*);


void free_CAMspec(void **none) {
  camspec_extra_free_();
}

double CAMspec_lkl(void* none, double* pars, error **err) {
  double lkl;
  
  camspec_extra_lkl_(&lkl,pars);
  return lkl;
}

void camspec_extra_init_(int*, int*,int*,int*,int*,int*, double*,double*,int*,double*,double*,double*,int*,int*,int*,int*,double*,double*,int*,int*,int*,double*,int*,int*,double*);
//SUBROUTINE CAMSPEC_EXTRA_INIT(iNspec, inX,ilminX,ilmaxX,inp,inpt, ic_inv,iX,ilmax_sz,isz_143_temp,iksz_temp,itszxcib_temp,ibeam_Nspec,inum_modes_per_beam,ibeam_lmax,icov_dim,ibeam_cov_inv,ibeam_modes,ihas_dust,ihas_calib,imarge_flag, imarge_mode,imarge_num, ikeep_num)

int _set_str(char* ptr, char* val, int mlen) {
  int sz;
  memset(ptr,' ',mlen);
  sprintf(ptr,"%s",val);
  sz=strlen(ptr);
  ptr[sz]=' ';
  return sz;


}
#define _flen_ 100

#ifndef CAMSPEC_V1
void camspec_extra_free_v3_();

void free_CAMspec_v3(void **none) {
  camspec_extra_free_v3_();
}


void camspec_extra_init_v3_(int*,char*,int*,char*,int*,char*,int*,char*,int*,
                         char*,int*,char*,int*,int*,char*,int*,char*,
                         int*,char*,int*,char*,int*,char*,
                         char*,int*,char*,int*,int*,int*,
                         int*,int*,int*,int*,int*,int*,double*,int*);
cmblkl* clik_CAMspec_v3_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  cmblkl *cing;
  int hk;
  int pre_marged;
  int *spec_flag,*lmins,*lmaxs;
  char like_file[_flen_], sz143_file[_flen_], tszxcib_file[_flen_], ksz_file[_flen_], beam_file[_flen_],data_vector[_flen_],camspec_fiducial_foregrounds[_flen_],camspec_fiducial_cl[_flen_];
  char cib217_file[_flen_],dust100_file[_flen_],dust143_file[_flen_],dust217_file[_flen_],dust143x217_file[_flen_]; 
  int l_like_file, l_sz143_file, l_tszxcib_file, l_ksz_file, l_beam_file,l_data_vector,l_camspec_fiducial_foregrounds,l_camspec_fiducial_cl;
  int l_cib217_file,l_dust100_file,l_dust143_file,l_dust217_file,l_dust143x217_file ;
  char directory_name[4096],pwd[4096];
  char *xnames_tot[300];
  char *nuisance;
  double bs_factor;
  int camspec_beam_mcmc_num, xdim,i,sz_prior;

  pre_marged = 1;

  hk = cldf_haskey(df,"pre_marged",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    pre_marged = cldf_readint(df,"pre_marged",err);
    forwardError(*err,__LINE__,NULL);
  }
  
  
  hk = cldf_haskey(df,"spec_flag",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    pre_marged = 0;
    lmins = cldf_readintarray(df,"spec_lmins",NULL,err);
    forwardError(*err,__LINE__,NULL);
    lmaxs = cldf_readintarray(df,"spec_lmaxs",NULL,err);
    forwardError(*err,__LINE__,NULL);
    spec_flag = cldf_readintarray(df,"spec_flag",NULL,err);
    forwardError(*err,__LINE__,NULL);
  } else {
    spec_flag = malloc_err(sizeof(int)*6,err);
    forwardError(*err,__LINE__,NULL);
    spec_flag[0]=1;
    spec_flag[1]=1;
    spec_flag[2]=1;
    spec_flag[3]=1;
    spec_flag[4]=1;
    spec_flag[5]=1;
    //memset(spec_flag,1,sizeof(int)*6);
    lmins = malloc_err(sizeof(int)*6,err);
    forwardError(*err,__LINE__,NULL);
    memset(lmins,0,sizeof(int)*6);
    lmaxs = malloc_err(sizeof(int)*6,err);
    forwardError(*err,__LINE__,NULL);
    memset(lmaxs,0,sizeof(int)*6);
  }
  //_DEBUGHERE_("%d %d %d %d %d %d",spec_flag[0],spec_flag[1],spec_flag[2],spec_flag[3],spec_flag[4],spec_flag[5]);
  l_like_file = _set_str(like_file,"like_file",_flen_);
  l_sz143_file = _set_str(sz143_file,"sz143_file",_flen_);
  l_tszxcib_file = _set_str(tszxcib_file,"tszxcib_file",_flen_);
  l_ksz_file = _set_str(ksz_file,"ksz_file",_flen_);
  l_beam_file = _set_str(beam_file,"beam_file",_flen_);

  hk = cldf_haskey(df,"cib_consistency_flag",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    l_cib217_file = _set_str(cib217_file,"cib_file",_flen_);
  } else {
    l_cib217_file = _set_str(cib217_file,"",_flen_);
  }

  hk = cldf_haskey(df,"dust_flag",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    l_dust100_file = _set_str(dust100_file,"dust100_file",_flen_);
    l_dust143_file = _set_str(dust143_file,"dust143_file",_flen_);
    l_dust217_file = _set_str(dust217_file,"dust217_file",_flen_);
    l_dust143x217_file = _set_str(dust143x217_file,"dust143x217_file",_flen_);
  } else {
    l_dust100_file = _set_str(dust100_file,"",_flen_);
    l_dust143_file = _set_str(dust143_file,"",_flen_);
    l_dust217_file = _set_str(dust217_file,"",_flen_);
    l_dust143x217_file = _set_str(dust143x217_file,"",_flen_);
  }

  l_data_vector = _set_str(data_vector," ",_flen_);
  hk = cldf_haskey(df,"data_vector_flag",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    int dv;

    dv = cldf_readint(df,"data_vector_flag",err);
    forwardError(*err,__LINE__,NULL);
    if (dv==1) {
      l_data_vector = _set_str(data_vector,"data_vector",_flen_);  
    }  
  }

  camspec_beam_mcmc_num = cldf_readint(df,"camspec_beam_mcmc_num",err);
  forwardError(*err,__LINE__,NULL);

  xdim = cldf_readint(df,"n_nuisance",err);
  forwardError(*err,__LINE__,NULL);
  

  nuisance = cldf_readstr(df,"nuisance",NULL,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<xdim;i++) {
    xnames_tot[i] = &(nuisance[i*256]);
  }  

  bs_factor = cldf_readfloat_default(df,"bs_factor",2.7,err);
  forwardError(*err,__LINE__,NULL);

  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);

  l_camspec_fiducial_foregrounds = _set_str(camspec_fiducial_foregrounds," ",_flen_);
  l_camspec_fiducial_cl = _set_str(camspec_fiducial_foregrounds," ",_flen_);
  if (pre_marged==0) {
    l_camspec_fiducial_foregrounds = _set_str(camspec_fiducial_foregrounds,"camspec_fiducial_foregrounds",_flen_);
    l_camspec_fiducial_cl = _set_str(camspec_fiducial_cl,"camspec_fiducial_cl",_flen_);
    hk = cldf_haskey(df,"minimum.theory_cl",err);
    forwardError(*err,__LINE__,NULL);
    if (hk==1) {
      int hs;
      hs = cldf_readint(df,"minimum.theory_cl",err);
      forwardError(*err,__LINE__,NULL);
      if (hs==1) {
        l_camspec_fiducial_cl = _set_str(camspec_fiducial_cl,"camspec_fiducial_cl.minimum.theory_cl",_flen_);    
      }
    }    
  }
  
  sz_prior = 0;
  hk = cldf_haskey(df,"sz_prior",err); 
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    sz_prior = cldf_readint(df,"sz_prior",err);
    forwardError(*err,__LINE__,NULL);    
  }
  // call the init fortran code here  
  //_DEBUGHERE_("%d",xdim);
  camspec_extra_init_v3_(&pre_marged,like_file,&l_like_file,sz143_file,&l_sz143_file,tszxcib_file,&l_tszxcib_file,ksz_file,&l_ksz_file,
                         beam_file,&l_beam_file,data_vector,&l_data_vector,&l_cib217_file,cib217_file,&l_dust100_file,dust100_file,
                         &l_dust143_file,dust143_file,&l_dust217_file,dust217_file,&l_dust143x217_file,dust143x217_file,
                         camspec_fiducial_foregrounds,&l_camspec_fiducial_foregrounds,camspec_fiducial_cl,&l_camspec_fiducial_cl,lmins,lmaxs,
                         spec_flag,&camspec_beam_mcmc_num,&xdim,&(ell[0]),&(ell[nell-1]),has_cl,&bs_factor,&sz_prior);

  cldf_external_cleanup(directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  // get some way to count the number of external parameters and pass them...

  cing = init_cmblkl(NULL, &CAMspec_lkl, 
                     &free_CAMspec_v3,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,xdim,err);
  forwardError(*err,__LINE__,NULL);
  
  cmblkl_set_names(cing, xnames_tot,err);
  forwardError(*err,__LINE__,NULL);
  
  free(lmins);
  free(lmaxs);
  free(spec_flag);
  free(nuisance);

  return cing;
}

void camspec_extra_init_v2_(int*,char*,int*,char*,int*,char*,int*,char*,int*,char*,int*,char*,int*,char*,int*,char*,int*,int*,int*,int*,int*,int*,int*,int*,int*,double*);

cmblkl* clik_CAMspec_v2_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  cmblkl *cing;
  int hk;
  int pre_marged;
  int *spec_flag,*lmins,*lmaxs;
  char like_file[_flen_], sz143_file[_flen_], tszxcib_file[_flen_], ksz_file[_flen_], beam_file[_flen_],data_vector[_flen_],camspec_fiducial_foregrounds[_flen_],camspec_fiducial_cl[_flen_]; 
  int l_like_file, l_sz143_file, l_tszxcib_file, l_ksz_file, l_beam_file,l_data_vector,l_camspec_fiducial_foregrounds,l_camspec_fiducial_cl;
  char directory_name[4096],pwd[4096];
  char *xnames_tot[300];
  char *nuisance;
  double bs_factor;
  int camspec_beam_mcmc_num, xdim,i;

  pre_marged = 1;

  hk = cldf_haskey(df,"pre_marged",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    pre_marged = cldf_readint(df,"pre_marged",err);
    forwardError(*err,__LINE__,NULL);
  }
  
  
  hk = cldf_haskey(df,"spec_flag",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    pre_marged = 0;
    lmins = cldf_readintarray(df,"spec_lmins",NULL,err);
    forwardError(*err,__LINE__,NULL);
    lmaxs = cldf_readintarray(df,"spec_lmaxs",NULL,err);
    forwardError(*err,__LINE__,NULL);
    spec_flag = cldf_readintarray(df,"spec_flag",NULL,err);
    forwardError(*err,__LINE__,NULL);
  } else {
    spec_flag = malloc_err(sizeof(int)*6,err);
    forwardError(*err,__LINE__,NULL);
    spec_flag[0]=1;
    spec_flag[1]=1;
    spec_flag[2]=1;
    spec_flag[3]=1;
    spec_flag[4]=1;
    spec_flag[5]=1;
    //memset(spec_flag,1,sizeof(int)*6);
    lmins = malloc_err(sizeof(int)*6,err);
    forwardError(*err,__LINE__,NULL);
    memset(lmins,0,sizeof(int)*6);
    lmaxs = malloc_err(sizeof(int)*6,err);
    forwardError(*err,__LINE__,NULL);
    memset(lmaxs,0,sizeof(int)*6);
  }
  //_DEBUGHERE_("%d %d %d %d %d %d",spec_flag[0],spec_flag[1],spec_flag[2],spec_flag[3],spec_flag[4],spec_flag[5]);
  l_like_file = _set_str(like_file,"like_file",_flen_);
  l_sz143_file = _set_str(sz143_file,"sz143_file",_flen_);
  l_tszxcib_file = _set_str(tszxcib_file,"tszxcib_file",_flen_);
  l_ksz_file = _set_str(ksz_file,"ksz_file",_flen_);
  l_beam_file = _set_str(beam_file,"beam_file",_flen_);

  l_data_vector = _set_str(data_vector," ",_flen_);
  hk = cldf_haskey(df,"data_vector_flag",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    int dv;

    dv = cldf_readint(df,"data_vector_flag",err);
    forwardError(*err,__LINE__,NULL);
    if (dv==1) {
      l_data_vector = _set_str(data_vector,"data_vector",_flen_);  
    }  
  }

  camspec_beam_mcmc_num = cldf_readint(df,"camspec_beam_mcmc_num",err);
  forwardError(*err,__LINE__,NULL);

  xdim = cldf_readint(df,"n_nuisance",err);
  forwardError(*err,__LINE__,NULL);
  

  nuisance = cldf_readstr(df,"nuisance",NULL,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<xdim;i++) {
    xnames_tot[i] = &(nuisance[i*256]);
  }  

  bs_factor = cldf_readfloat_default(df,"bs_factor",2.7,err);
  forwardError(*err,__LINE__,NULL);

  cldf_external(df,directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);

  l_camspec_fiducial_foregrounds = _set_str(camspec_fiducial_foregrounds," ",_flen_);
  l_camspec_fiducial_cl = _set_str(camspec_fiducial_foregrounds," ",_flen_);
  if (pre_marged==0) {
    l_camspec_fiducial_foregrounds = _set_str(camspec_fiducial_foregrounds,"camspec_fiducial_foregrounds",_flen_);
    l_camspec_fiducial_cl = _set_str(camspec_fiducial_cl,"camspec_fiducial_cl",_flen_);
  }
  
  // call the init fortran code here  
  //_DEBUGHERE_("%d",xdim);
  camspec_extra_init_v2_(&pre_marged,like_file,&l_like_file,sz143_file,&l_sz143_file,tszxcib_file,&l_tszxcib_file,ksz_file,&l_ksz_file,beam_file,&l_beam_file,data_vector,&l_data_vector,camspec_fiducial_foregrounds,&l_camspec_fiducial_foregrounds,camspec_fiducial_cl,&l_camspec_fiducial_cl,lmins,lmaxs,spec_flag,&camspec_beam_mcmc_num,&xdim,&(ell[0]),&(ell[nell-1]),has_cl,&bs_factor);

  cldf_external_cleanup(directory_name,pwd,err);
  forwardError(*err,__LINE__,NULL);
  
  // get some way to count the number of external parameters and pass them...

  cing = init_cmblkl(NULL, &CAMspec_lkl, 
                     &free_CAMspec,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,xdim,err);
  forwardError(*err,__LINE__,NULL);
  
  cmblkl_set_names(cing, xnames_tot,err);
  forwardError(*err,__LINE__,NULL);
  
  free(lmins);
  free(lmaxs);
  free(spec_flag);
  free(nuisance);

  return cing;
}
#undef _flen_
#endif


cmblkl* clik_CAMspec_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  int bok;
  cmblkl *cing;
  char likefile[2048],szfile[2048];
  int xcase;
  double xdim;
  char *xnames_def[] = {"A_ps_100","A_ps_143", "A_ps_217", "A_cib_143", "A_cib_217", "A_sz",  
                      "r_ps", "r_cib","n_Dl_cib","cal_100", "cal_143", "cal_217","xi_sz_cib","A_ksz","A_dust"};
  
  char *xnames_1[] = {"A_ps_100","A_ps_143", "A_ps_217", "A_cib_143", "A_cib_217", "A_sz",  
                      "r_ps", "r_cib","cal0", "cal1", "cal2"};
  char *xnames_0[] = {"A_ps_143", "A_cib_143", "A_sz",  
                      "cal1", "cal2"};

  char *xnames_tot[300];

  int Nspec, nX, lmax_sz;
  int *lmaxX, *lminX, *np, *npt;
  double *X,*c_inv,*tsz,*ksz,*tszXcib;
  int X_sz,c_inv_sz,sz_temp_sz;
  int beam_Nspec,num_modes_per_beam,beam_lmax,cov_dim;
  double *beam_cov_inv,*beam_modes;
  int beam_cov_inv_sz,beam_modes_sz;
  int i,j,cnt;
  int has_dust, has_calib_prior;
  int hk;
  int marge_num,keep_num;
  int *marge_flag;
  double *marge_mode;
  int sz;
  double bs_factor;
  int cv;

  //int frq[] = {100,143,217};
  
  camspec_extra_only_one_(&bok);
  testErrorRet(bok!=0,-100,"CAMspec already initialized",*err,__LINE__,NULL);

#ifndef CAMSPEC_V1
  hk = cldf_haskey(df,"camspec_version",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    cv = cldf_readint(df,"camspec_version",err);
    forwardError(*err,__LINE__,NULL);
    if (cv==2) {
      cing = clik_CAMspec_v2_init(df, nell,ell,has_cl,unit,wl,bins,nbins,err);
      forwardError(*err,__LINE__,NULL);
      return cing;
    }
    if (cv==3) {
      cing = clik_CAMspec_v3_init(df, nell,ell,has_cl,unit,wl,bins,nbins,err);
      forwardError(*err,__LINE__,NULL);
      return cing;
    }
  }
#endif

  Nspec = cldf_readint(df,"Nspec",err);
  forwardError(*err,__LINE__,NULL);
  
  lminX = cldf_readintarray(df,"lminX",&Nspec,err);
  forwardError(*err,__LINE__,NULL);
  lmaxX = cldf_readintarray(df,"lmaxX",&Nspec,err);
  forwardError(*err,__LINE__,NULL);

  np = cldf_readintarray(df,"np",&Nspec,err);
  forwardError(*err,__LINE__,NULL);
  npt = cldf_readintarray(df,"npt",&Nspec,err);
  forwardError(*err,__LINE__,NULL);
 
  nX = cldf_readint(df,"nX",err);
  forwardError(*err,__LINE__,NULL);

  lmax_sz = cldf_readint(df,"lmax_sz",err);
  forwardError(*err,__LINE__,NULL);

  X_sz = -1;
  X = cldf_readfloatarray(df,"X",&X_sz, err);
  forwardError(*err,__LINE__,NULL);
  
  c_inv_sz = -1;
  c_inv = cldf_readfloatarray(df,"c_inv",&c_inv_sz, err);
  forwardError(*err,__LINE__,NULL);

  sz_temp_sz = -1;
  tsz = cldf_readfloatarray(df,"tsz",&sz_temp_sz, err);
  forwardError(*err,__LINE__,NULL);
  sz_temp_sz = -1;
  ksz = cldf_readfloatarray(df,"ksz",&sz_temp_sz, err);
  forwardError(*err,__LINE__,NULL);
  sz_temp_sz = -1;
  tszXcib = cldf_readfloatarray(df,"tszXcib",&sz_temp_sz, err);
  forwardError(*err,__LINE__,NULL);


  beam_Nspec = cldf_readint(df,"beam_Nspec",err);
  forwardError(*err,__LINE__,NULL);

  if (beam_Nspec==0) {
    num_modes_per_beam = 0;
    beam_lmax = 0;
    cov_dim   = 0;
    beam_cov_inv = malloc_err(sizeof(double)*1,err);
    forwardError(*err,__LINE__,NULL);
    beam_modes = malloc_err(sizeof(double)*1,err);
    forwardError(*err,__LINE__,NULL);
  } else {
    beam_lmax = cldf_readint(df,"beam_lmax",err);
    forwardError(*err,__LINE__,NULL);
    num_modes_per_beam = cldf_readint(df,"num_modes_per_beam",err);
    forwardError(*err,__LINE__,NULL);
    cov_dim = cldf_readint(df,"cov_dim",err);
    forwardError(*err,__LINE__,NULL);
  
    beam_cov_inv_sz = -1;
    beam_cov_inv = cldf_readfloatarray(df,"beam_cov_inv",&beam_cov_inv_sz, err);
    forwardError(*err,__LINE__,NULL);
    beam_modes_sz = -1;
    beam_modes = cldf_readfloatarray(df,"beam_modes",&beam_modes_sz, err);
    forwardError(*err,__LINE__,NULL);
  }


  has_dust =  cldf_readint_default(df, "has_dust",0,err);
  forwardError(*err,__LINE__,NULL);

  has_calib_prior = cldf_readint_default(df, "has_calib_prior",1,err);

  hk = cldf_haskey(df,"marge_flag",err);
  forwardError(*err,__LINE__,NULL);
  marge_num = 0;
  if (hk==0) {
    marge_flag = malloc_err(sizeof(int)*cov_dim,err);
    forwardError(*err,__LINE__,NULL);
    for(i=0;i<cov_dim;i++) {
      marge_flag[i] = 0;
    }  
  } else {
    sz = cov_dim;
    marge_flag = cldf_readintarray(df,"marge_flag",&sz,err);
    forwardError(*err,__LINE__,NULL);
    for(i=0;i<cov_dim;i++) {
      marge_num += marge_flag[i];  
    }
  }
  keep_num = cov_dim - marge_num;

  if (marge_num>0) {
    sz = keep_num * marge_num;
    marge_mode = cldf_readfloatarray(df,"marge_mode",&sz,err);
    forwardError(*err,__LINE__,NULL);
  } else {
    marge_num = 1;
    marge_mode = malloc_err(sizeof(double)*keep_num,err);
    forwardError(*err,__LINE__,NULL);
  }

  bs_factor = cldf_readfloat_default(df,"bs_factor",1,err);
  forwardError(*err,__LINE__,NULL);
  

  camspec_extra_init_(&Nspec, &nX,lminX,lmaxX,np,npt,c_inv,X,&lmax_sz, tsz,ksz,tszXcib,&beam_Nspec,&num_modes_per_beam,&beam_lmax,&cov_dim,beam_cov_inv,beam_modes,&has_dust,&has_calib_prior,marge_flag,marge_mode,&marge_num,&keep_num,&bs_factor);    
  
  //camspec_extra_getcase_(&xcase);
  xdim = 14 + has_dust + keep_num;
  
  for(i=0;i<14+ has_dust;i++) {
    xnames_tot[i] = xnames_def[i];
  }

  cnt = 14+ has_dust;
  
  for(i=0;i<beam_Nspec;i++) {
    for(j=0;j<num_modes_per_beam;j++) {
      if (marge_flag[cnt-14-has_dust]==0) {
      xnames_tot[cnt] = malloc_err(sizeof(char)*50,err);
      forwardError(*err,__LINE__,NULL);

      sprintf(xnames_tot[cnt],"Bm_%d_%d",i+1,j+1);
      cnt++;        
      }
    }
  }

  cing = init_cmblkl(NULL, &CAMspec_lkl, 
                     &free_CAMspec,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,xdim,err);
  forwardError(*err,__LINE__,NULL);
  /*if (xcase==0) {
    cmblkl_set_names(cing, xnames_0,err);
  forwardError(*err,__LINE__,NULL);  
  } else {*/
  cmblkl_set_names(cing, xnames_tot,err);
  forwardError(*err,__LINE__,NULL);
  //}
  
  free(X);
  free(c_inv);
  free(ksz);
  free(tsz);
  free(tszXcib);
  free(beam_modes);
  free(beam_cov_inv);
  free(marge_flag);
  free(marge_mode);
  for(i=14 + has_dust;i<xdim;i++) {
    free(xnames_tot[i]);
  }

  return cing;
}

double* camspec_get_fg(void* camclik,double *par,int lmax,error **err) {
  double *res;

  res = malloc_err(sizeof(double)*(lmax+1)*4,err);
  forwardError(*err,__LINE__,NULL);

  camspec_extra_fg_(res,par,&lmax);

  return res;
}
