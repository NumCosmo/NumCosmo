#include "clik_helper.h"
#include "smica.h"

typedef SmicaComp * clik_smica_comp_init_func(cldf * df,int nb, int mT,int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err);

# ifndef REL_EXT
cmblkl* clik_smica_init(cldf * df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  Smica *ing;
  cmblkl *cing;  
  double zero;
  int n;
  int *Ml;
  double *wq,*rq_hat,*rq_0;
  int m;
  SmicaComp **SCs;
  int ic;
  int mT,mP,nc;
  double *A_cmb;
  Smica *smic;
  int ncl,icl,nb;
  size_t ddum;
  int cnt,xdim;
  char **xnames;
  parname *xnames_buf;
  int hk;
  int mT_plus_mP;

  zero = 0;

  ncl = 0;
  for(icl=0;icl<6;icl++) {
    if (has_cl[icl]==1) {
      ncl++;
    }
  }
  
  nb = nell;
  if (nbins!=0) {
    nb = nbins/ncl;  
  }  
    
  // try to read the bin weights
  wq = NULL;
  hk = cldf_haskey(df,"wq",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    wq = cldf_readfloatarray(df,"wq",&nb,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  
  // read the number of channels
  mT = cldf_readint(df,"m_channel_T",err);
  forwardError(*err,__LINE__,NULL);
  mP = cldf_readint(df,"m_channel_P",err);
  forwardError(*err,__LINE__,NULL);
  
  m = mtot(mT, mP, has_cl);
  mT_plus_mP = mT+mP;

  // read rq_hat
  int nrq;
  
  nrq = nb*m*m;
  rq_hat = cldf_readfloatarray(df,"Rq_hat",&nrq,err);
  forwardError(*err,__LINE__,NULL);  
  
  // try to read rq_0
  rq_0 = NULL;
  hk = cldf_haskey(df,"Rq_0",err);
  forwardError(*err,__LINE__,NULL);  
  if (hk==1) {
    rq_0 = cldf_readfloatarray(df,"Rq_0",&nrq,err);
    forwardError(*err,__LINE__,NULL);  
  }
  
  // how many components ?
  nc = cldf_readint(df,"n_component",err);
  forwardError(*err,__LINE__,NULL);
  
  SCs = malloc_err(sizeof(SmicaComp*) * nc,err);
  forwardError(*err,__LINE__,NULL);
  
  // now deal with the CMB component
  // read A_cmb
  A_cmb = cldf_readfloatarray(df,"A_cmb",&mT_plus_mP,err);
  forwardError(*err,__LINE__,NULL);    

  // init cmb comp
  SCs[0] = comp_CMB_init(nb, mT,mP, has_cl, A_cmb, err);
  forwardError(*err,__LINE__,NULL);    
  SC_set_compname(SCs[0],"CMB");
  
  free(A_cmb);
  
  // deal with other components
  xdim = 0;
  for(ic=1;ic<nc;ic++) {
    char cur_cmp[256];
    char cur_cmp_tot[256];
    clik_smica_comp_init_func *smica_dl_init;
    void* dlhandle;
    char init_func_name[256];
    parname comp_type;
    cldf *comp_df;

#ifdef HAS_RTLD_DEFAULT 
    dlhandle = RTLD_DEFAULT;
#else
    dlhandle = NULL;
#endif
    
    SCs[ic] = NULL;
    
    sprintf(cur_cmp,"component_%d",ic);
    comp_df = cldf_openchild(df,cur_cmp,err);
    forwardError(*err,__LINE__,NULL);

    // get type
    char * cmt;

    cmt = cldf_readstr(comp_df,"component_type",NULL,err);
    forwardError(*err,__LINE__,NULL);
    sprintf(comp_type,"%s",cmt);
    free(cmt);

    sprintf(init_func_name,"clik_smica_comp_%s_init",comp_type);
    smica_dl_init = dlsym(dlhandle,init_func_name);
    testErrorRetVA(smica_dl_init==NULL,-1111,"Cannot initialize smica component type %s from %s dl error : %s",*err,__LINE__,NULL,comp_type,df->root,dlerror()); 
 
    SCs[ic] = smica_dl_init(comp_df,nb,mT,mP, nell, ell, has_cl, unit, wl, bins,nb,err);
 
    forwardError(*err,__LINE__,NULL);
    sprintf(comp_type,"%s_%d",comp_type,ic);
    SC_set_compname(SCs[ic],comp_type);
      
    cldf_close(&comp_df);
    
    xdim += SCs[ic]->ndim;
 
  
  }
  
  // deal with names and xdims
  if (xdim!=0) {
    xnames = malloc_err(sizeof(char*)*xdim,err);
    forwardError(*err,__LINE__,NULL);
  
    xnames_buf = malloc_err(sizeof(parname)*xdim,err);
    forwardError(*err,__LINE__,NULL);
    cnt = 0;
    for(ic=1;ic<nc;ic++) {
      int ix;
      for(ix=0;ix<SCs[ic]->ndim;ix++) {
        if (SCs[ic]->names!=NULL) {
          xnames[cnt]=(char*)&(SCs[ic]->names[ix]);
        } else {
          sprintf(xnames_buf[cnt],"SMICA_COMP_%d_%d",ic,ix);
          xnames[cnt]=(char*)&(xnames_buf[cnt]);
        }
        cnt++;
      }
    }
  }
    
  smic = Smica_init(nb, wq, m, rq_hat, rq_0, nc, SCs,err);
  forwardError(*err,__LINE__,NULL);
    
  // deal with criterion
  hk = cldf_haskey(df,"criterion",err);
  forwardError(*err,__LINE__,NULL);
  
  if (hk == 1) {
    char *crit_name;
    
    crit_name = cldf_readstr(df,"criterion",NULL,err);
    forwardError(*err,__LINE__,NULL);
    if(strcmp(crit_name,"classic")==0) {
      //nothing to do here, we are fine
    } else if(strcmp(crit_name,"gauss")==0) {
      double *quad_crit;
      int *mask;
      int *ordering;
      int nqu;

      nqu = -1;
      quad_crit = cldf_readfloatarray(df,"criterion_gauss_mat",&nqu,err);
      forwardError(*err,__LINE__,NULL);
      mask = NULL;
      hk = cldf_haskey(df,"criterion_gauss_mask",err);
      forwardError(*err,__LINE__,NULL);
      if (hk == 1) { 
        int i;
        nqu = nb*m*m;
        mask = cldf_readintarray(df,"criterion_gauss_mask",&nqu,err);
        forwardError(*err,__LINE__,NULL);
      }
      ordering = NULL;
      hk = cldf_haskey(df,"criterion_gauss_ordering",err);
      forwardError(*err,__LINE__,NULL);
      if (hk == 1) { 
        int i;
        nqu = m*(m+1);
        ordering = cldf_readintarray(df,"criterion_gauss_ordering",&nqu,err);
        forwardError(*err,__LINE__,NULL);
        //_DEBUGHERE_("","");
      }

      smica_set_crit_gauss(smic, quad_crit,mask,ordering,err);
      forwardError(*err,__LINE__,NULL);
      
      if (mask!=NULL) {
        free(mask);  
      }
      if (ordering!=NULL) {
        free(ordering);  
      }
      free(quad_crit);
    } else {
      testErrorRetVA(1==1,hdf5_base,"does not understand criterion '%s' in %s",*err,__LINE__,NULL,crit_name,df->root);
    }
    free(crit_name);
  }
  
  free(rq_hat);
  
  if (rq_0!=NULL) {
    free(rq_0);
  }    
  
  if (wq!=NULL) {
    free(wq);
  }  
  
  cing = init_cmblkl(smic, &Smica_lkl, 
                     &free_Smica,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,xdim,err);
  forwardError(*err,__LINE__,NULL);
  if (xdim!=0) {
    cmblkl_set_names(cing, xnames,err);
    forwardError(*err,__LINE__,NULL);
  
    free(xnames);
    free(xnames_buf);
  }
  
  
  
  return cing;  
}
#endif


SmicaComp * clik_smica_comp_1d_init(cldf *df,int nb, int mT, int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  double *someA;
  int nA;
  SmicaComp *SC;
  int hk;
  int m;

  m = mtot(mT,mP,has_cl);
  // try to read A

  someA = NULL;
  
  //hstat = H5LTfind_attribute(comp_id, "A");
  hk = cldf_haskey(df,"A",err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    nA = -1;
    //someA = hdf5_double_attarray(comp_id,cur_lkl,"A",&nA,err);
    someA = cldf_readfloatarray(df,"A",&nA,err);
    forwardError(*err,__LINE__,NULL);
  }
  testErrorRetVA(someA!=NULL && nA!=m,-100,"Not enough data in %s (expected %d got %d)",*err,__LINE__,NULL,df->root,nA,m);
  SC = comp_1D_init(nb, m, someA, err);
  forwardError(*err,__LINE__,NULL);    
  if (someA!=NULL) {
    free(someA);
  }
  return SC;
}

SmicaComp * clik_smica_comp_nd_init(cldf *df,int nb, int mT,int mP,int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  double *someA;
  int nA;
  int nd;
  SmicaComp *SC;
  int hk;
  // try to read A
  int m;

  m = mtot(mT,mP,has_cl);
  
  nd = cldf_readint(df,"nd",err);
  forwardError(*err,__LINE__,NULL);    
  
  //hstat = H5LTget_attribute_int(comp_id, ".", "nd",  &nd);
  //testErrorRetVA(hstat<0,hdf5_base,"cannot read nd in component %s (got %d)",*err,__LINE__,NULL,cur_lkl,hstat);

  someA = NULL;
  hk = cldf_haskey(df,"A",err);
  forwardError(*err,__LINE__,NULL);
  //hstat = H5LTfind_attribute(comp_id, "A");
  if (hk==1) {
    nA = -1;
    //someA = hdf5_double_attarray(comp_id,cur_lkl,"A",&nA,err);
    someA = cldf_readfloatarray(df,"A",&nA,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  testErrorRetVA(someA!=NULL && nA!=m*nd,-100,"Not enough data in %s (expected %d got %d)",*err,__LINE__,NULL,df->root,nA,m*nd);
  SC = comp_nD_init(nb, m, nd, someA, err);
  forwardError(*err,__LINE__,NULL);    
  if (someA!=NULL) {
    free(someA);
  }
  return SC;
}



SmicaComp * clik_smica_comp_calTP_init(cldf *df,int nb, int mT,int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {
  SmicaComp *SC; 
  int npar,i;
  int *im,*jm;
  double *tpl;
  int tot_tpl;
  char *bnames, **xnames;
  int m;
  int dsz;
  double *w;
  int *other;
  int hk,im1,im2;

  m = mtot(mT,mP,has_cl);
  
  npar = -1;
  im =  cldf_readintarray(df,"im",&npar,err);
  forwardError(*err,__LINE__,NULL);    


  hk = cldf_haskey(df,"w",err);
  forwardError(*err,__LINE__,NULL);
  //hstat = H5LTfind_attribute(comp_id, "A");
  if (hk==1) {
    int sz;
    sz = m*m*2;
    w = cldf_readfloatarray(df,"w",&sz,err);
    forwardError(*err,__LINE__,NULL);
    other = cldf_readintarray(df,"other",&sz,err);
    forwardError(*err,__LINE__,NULL);
  } else {
    //_DEBUGHERE_("%d",m);
    w = malloc_err(sizeof(double)*m*m*2,err);
    forwardError(*err,__LINE__,NULL);
    other = malloc_err(sizeof(int)*m*m*2,err);
    forwardError(*err,__LINE__,NULL);
    for(im1=0;im1<m;im1++) {
      for(im2=im1;im2<m;im2++) {
        w[(im1*m+im2)*2 + 0] = 1;
        w[(im1*m+im2)*2 + 1] = 0;
        w[(im2*m+im1)*2 + 0] = 1;
        w[(im2*m+im1)*2 + 1] = 0;
        other[(im1*m+im2)*2 + 0] = im1;
        other[(im1*m+im2)*2 + 1] = im2;
        other[(im2*m+im1)*2 + 0] = im2;
        other[(im2*m+im1)*2 + 1] = im1;
        //_DEBUGHERE_("%d w %d,%d->%g,%g, other %d,%d->%d,%d",(im1*m+im2)*2,im1,im2,w[(im1*m+im2)*2 + 0],w[(im1*m+im2)*2 + 1],im1,im2,other[(im1*m+im2)*2 + 0],other[(im1*m+im2)*2 + 1])
        //_DEBUGHERE_("%d %d %d %d | %d",(im1*m+im2)*2 + 0, (im1*m+im2)*2 + 1, (im2*m+im1)*2 + 0, (im2*m+im1)*2 + 1,m*m*2);
      }
    }
  }
  

  SC = comp_calTP_init(nb, mT,mP,has_cl, npar, im,w,other,  err);
  forwardError(*err,__LINE__,NULL);    
  
  free(im);
  free(w);
  free(other);

  dsz = -1;
  bnames = cldf_readstr(df,"names",&dsz, err);
  forwardError(*err,__LINE__,NULL); 
  xnames = malloc_err(sizeof(char*)*npar,err);
  for(i=0;i<npar;i++) {
    xnames[i] =&(bnames[i*256]);
  } 
  SC_setnames(SC, xnames, err);
  forwardError(*err,__LINE__,NULL);
  free(xnames); 
  free(bnames); 
  
  return SC;
}

SmicaComp * clik_smica_comp_icalTP_init(cldf *df,int nb, int mT,int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {
  SmicaComp *SC; 
  int npar,i;
  int *im,*jm;
  double *tpl;
  int tot_tpl;
  char *bnames, **xnames;
  int m;
  int dsz;
  double *w;
  int *other;
  int hk,im1,im2;

  m = mtot(mT,mP,has_cl);
  
  npar = -1;
  im =  cldf_readintarray(df,"im",&npar,err);
  forwardError(*err,__LINE__,NULL);    


  hk = cldf_haskey(df,"w",err);
  forwardError(*err,__LINE__,NULL);
  //hstat = H5LTfind_attribute(comp_id, "A");
  if (hk==1) {
    int sz;
    sz = m*m*2;
    w = cldf_readfloatarray(df,"w",&sz,err);
    forwardError(*err,__LINE__,NULL);
    other = cldf_readintarray(df,"other",&sz,err);
    forwardError(*err,__LINE__,NULL);
  } else {
    //_DEBUGHERE_("%d",m);
    w = malloc_err(sizeof(double)*m*m*2,err);
    forwardError(*err,__LINE__,NULL);
    other = malloc_err(sizeof(int)*m*m*2,err);
    forwardError(*err,__LINE__,NULL);
    for(im1=0;im1<m;im1++) {
      for(im2=im1;im2<m;im2++) {
        w[(im1*m+im2)*2 + 0] = 1;
        w[(im1*m+im2)*2 + 1] = 0;
        w[(im2*m+im1)*2 + 0] = 1;
        w[(im2*m+im1)*2 + 1] = 0;
        other[(im1*m+im2)*2 + 0] = im1;
        other[(im1*m+im2)*2 + 1] = im2;
        other[(im2*m+im1)*2 + 0] = im2;
        other[(im2*m+im1)*2 + 1] = im1;
        //_DEBUGHERE_("%d w %d,%d->%g,%g, other %d,%d->%d,%d",(im1*m+im2)*2,im1,im2,w[(im1*m+im2)*2 + 0],w[(im1*m+im2)*2 + 1],im1,im2,other[(im1*m+im2)*2 + 0],other[(im1*m+im2)*2 + 1])
        //_DEBUGHERE_("%d %d %d %d | %d",(im1*m+im2)*2 + 0, (im1*m+im2)*2 + 1, (im2*m+im1)*2 + 0, (im2*m+im1)*2 + 1,m*m*2);
      }
    }
  }
  

  SC = comp_icalTP_init(nb, mT,mP,has_cl, npar, im,w,other,  err);
  forwardError(*err,__LINE__,NULL);    
  
  free(im);
  free(w);
  free(other);

  dsz = -1;
  bnames = cldf_readstr(df,"names",&dsz, err);
  forwardError(*err,__LINE__,NULL); 
  xnames = malloc_err(sizeof(char*)*npar,err);
  for(i=0;i<npar;i++) {
    xnames[i] =&(bnames[i*256]);
  } 
  SC_setnames(SC, xnames, err);
  forwardError(*err,__LINE__,NULL);
  free(xnames); 
  free(bnames); 
  
  return SC;
}

SmicaComp * clik_smica_comp_beamTP_init(cldf *df,int nb, int mT,int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {
  SmicaComp *SC; 
  int npar,i;
  int *im,*jm;
  double *tpl;
  int tot_tpl;
  char *bnames, **xnames;
  int m;
  int dsz;
  double *w;
  int *other;
  int hk,im1,im2;
  int neigen;
  double *modes;

  m = mtot(mT,mP,has_cl);
  
  neigen = cldf_readint(df,"neigen",err);
  forwardError(*err,__LINE__,NULL);    

  npar = cldf_readint(df,"npar",err);
  forwardError(*err,__LINE__,NULL);

  dsz = m*m*neigen;
  im =  cldf_readintarray(df,"im",&dsz,err);
  forwardError(*err,__LINE__,NULL);    

  dsz = m*m*neigen*nb;
  modes = cldf_readfloatarray(df,"modes",&dsz,err);
  forwardError(*err,__LINE__,NULL);    

  SC = comp_beamTP_init(nb, mT,mP,has_cl,npar, im, neigen, modes ,err);
  forwardError(*err,__LINE__,NULL);    

  free(im);
  free(modes);

  dsz = -1;
  bnames = cldf_readstr(df,"names",&dsz, err);
  forwardError(*err,__LINE__,NULL); 
  xnames = malloc_err(sizeof(char*)*npar,err);
  for(i=0;i<npar;i++) {
    xnames[i] =&(bnames[i*256]);
  } 
  SC_setnames(SC, xnames, err);
  forwardError(*err,__LINE__,NULL);
  free(xnames); 
  free(bnames); 

  return SC;
}

SmicaComp * clik_smica_comp_totcal_init(cldf *df,int nb, int mT,int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {
  SmicaComp *SC; 
  int npar,i;
  int *im,*jm;
  double *tpl;
  int tot_tpl;
  char *calname;
  int m;
  int dsz;
  double *w;
  int *other;
  int hk,im1,im2;
  int neigen;
  double *modes;

  m = mtot(mT,mP,has_cl);
  
  dsz = -1;
  calname = cldf_readstr(df,"calname",&dsz, err);
  forwardError(*err,__LINE__,NULL);  

  SC = comp_totcal_init(nb, mT,mP,has_cl,err);
  forwardError(*err,__LINE__,NULL);    

  SC_setnames(SC, &calname, err);
  forwardError(*err,__LINE__,NULL);
  
  free(calname); 

  return SC;
}

SmicaComp * clik_smica_comp_totcalP_init(cldf *df,int nb, int mT,int mP, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {
  SmicaComp *SC; 
  int npar,i;
  int *im,*jm;
  double *tpl;
  int tot_tpl;
  char *calname;
  int m;
  int dsz;
  double *w;
  int *other;
  int hk,im1,im2;
  int neigen;
  double *modes;

  m = mtot(mT,mP,has_cl);
  
  dsz = -1;
  calname = cldf_readstr(df,"calnameP",&dsz, err);
  forwardError(*err,__LINE__,NULL);  

  SC = comp_totcalP_init(nb, mT,mP,has_cl,err);
  forwardError(*err,__LINE__,NULL);    

  SC_setnames(SC, &calname, err);
  forwardError(*err,__LINE__,NULL);
  
  free(calname); 

  return SC;
}

