#include "clik_parametric.h"
#include "clik_helper.h"
#include "clik_parametric_addon.h"

typedef struct {
  parametric *p_model;
  double *rq,*wrq;
  int nell,nbins,m;
  double unit;
  int *ell;
  double *bins,*wl;
  double *wbins;
  int *bi,*bo;
  int bn;
  double *A;
  } parametric_smica;

void comp_parametric_update(void* data,double* locpars, double* rq, error **err) {
  parametric_smica *p_pay;
  SmicaComp* SC;
  double *wl,*wl0,one;
  int inc,il,im;
  double res;
  int m, ndet;
  //double r10[16];
  char nm[3000];

  SC = data;
  p_pay = SC->data;

  //_DEBUGHERE_("%g",locpars[0]);

  parametric_compute(p_pay->p_model, locpars, p_pay->rq, NULL, err);
  forwardError(*err,__LINE__,);
  
  //sprintf(nm,"rq_%s.dat",SC->comp_name);
  //write_bin_vector(p_pay->rq, nm, sizeof(double)*(p_pay->nell*p_pay->p_model->ndet*p_pay->p_model->ndet), err);   
  //forwardError(*err,__LINE__,);
  
  m = p_pay->m;
  ndet = p_pay->p_model->ndet;
  //_DEBUGHERE_("%d %d",m,ndet);

  // apply wl and binning
  one=1;
  if (p_pay->wl==NULL) {
    wl0 = &one;
    inc = 0;
  } else {
    wl0 = p_pay->wl;
    inc = 1;
  }
  
  //_DEBUGHERE_("","");
  wl=wl0;
  for(il=0;il<p_pay->nell;il++) {
    int ip;
    int im1,im2;
    ip = il * ndet * ndet;

    for(im1=0;im1<ndet;im1++) {
      for(im2=0;im2<ndet;im2++) {
        p_pay->rq[il*ndet*ndet+im1*ndet+im2] = p_pay->rq[il*ndet*ndet+im1*ndet+im2] * *wl * p_pay->unit * p_pay->A[im1]*p_pay->A[im2];  
      }
    }
    wl+=inc;
  }
  //_DEBUGHERE_("%g %g %g %g",egfs_pay->A[0],egfs_pay->A[1],egfs_pay->A[2],egfs_pay->A[3])
  //_DEBUGHERE_("%g %g %g",egfs_pay->rq[0],egfs_pay->rq[2],egfs_pay->unit);
  
  // apply binning if needed
  //sprintf(nm,"rq_before_%s.dat",SC->comp_name);
  //write_bin_vector(rq, nm, sizeof(double)*(p_pay->nbins*p_pay->p_model->ndet*p_pay->p_model->ndet), err);   
  //forwardError(*err,__LINE__,);
  
  if (p_pay->bins!=NULL) {
    char transa,transb;
    int npar;
    double done,dzero;
    int one,nbns,nell;
    int ndim;
  //  int ii;
  //  double rq0,rq2;
  //  double poc;

    transa='N';
    transb='N';
    ndim = m*m;
    one = 1;
    done = 1.;
    dzero = 0;
    nbns = p_pay->nbins;
    nell = p_pay->nell;
    
    /*rq0=rq[0];
    rq2=rq[2];
    _DEBUGHERE_("%g %g",rq[0],rq[2]);
    _DEBUGHERE_("%g %g",rq[0]-rq0,rq[2]-rq2);
    _DEBUGHERE_("m %d n %d k %d",ndim, nbns,nell);*/  
    //_DEBUGHERE_("avant egfs","");
    //memset(r10,0,sizeof(double)*16);

    //printMat(&rq[10*ndim], egfs_pay->m,egfs_pay->m);
    //printMat(r10, egfs_pay->m,egfs_pay->m);
    //{
    //  int il,iq,if1,if2;
    //  for(il=0;il<nell;il++) {
    //    for(iq=0;iq<nbns;iq++) {
    //      for(if1=0;if1<p_pay->m;if1++) {
    //        for(if2=0;if2<p_pay->m;if2++) {
    //          rq[iq*ndim+if1*p_pay->m+if2] += p_pay->bins[iq*nell+il] * p_pay->rq[il*p_pay->m*p_pay->m+if1*p_pay->m+if2];
    //          /*if (iq==10)
    //            r10[if1*egfs_pay->m+if2] += egfs_pay->bins[iq*nell+il] * egfs_pay->rq[il+if1*egfs_pay->m*egfs_pay->nell+if2*egfs_pay->nell];*/
    //        }  
    //      }
    //    }
    //  }
    //}
    {
     int il,iq,if1,if2,bb;
     double w;
      bb=0;
      for(iq=0;iq<nbns;iq++) {
        for(il=p_pay->bi[iq];il<p_pay->bo[iq];il++) {
          w = p_pay->wbins[bb];
          bb++;
          for(if1=0;if1<ndet;if1++) {
            for(if2=0;if2<ndet;if2++) {
              //_DEBUGHERE_("%d %d %g %d,%d %g %g",il,iq,w,if1,if2,rq[iq*ndim+if1*m+if2],p_pay->rq[il*ndet*ndet+if1*ndet+if2]);
              rq[iq*ndim+if1*m+if2] += w * p_pay->rq[il*ndet*ndet+if1*ndet+if2];
              /*if (iq==10)
                  r10[if1*egfs_pay->m+if2] += egfs_pay->bins[iq*nell+il] * egfs_pay->rq[il+if1*egfs_pay->m*egfs_pay->nell+if2*egfs_pay->nell];*/
            }  
          }
        }
      }
    }
    //dgemm(&transa, &transb, &ndim, &nbns,&nell, &done, egfs_pay->rq, &ndim, egfs_pay->bins, &nell, &done, rq, &ndim);
    /*_DEBUGHERE_("apres egfs","");
    printMat(&rq[10*ndim], egfs_pay->m,egfs_pay->m);
    printMat(r10, egfs_pay->m,egfs_pay->m);*/
    /*_DEBUGHERE_("","");
    poc = 0;
    for(ii=0;ii<10;ii++) {
      _DEBUGHERE_("%g %g %g",egfs_pay->rq[ii*ndim],egfs_pay->rq[ii*ndim+2],egfs_pay->bins[ii]);
      poc += egfs_pay->rq[ii*ndim] * egfs_pay->bins[ii];
    }
    _DEBUGHERE_("%g %g %g",rq[0]-rq0,rq[2]-rq2,poc);*/
  } else {
    int if1,if2;
    for(il=0;il<p_pay->nell;il++) {
      for(if1=0;if1<ndet;if1++) {
        for(if2=0;if2<ndet;if2++) {
          rq[il*m*m+if1*m+if2] += p_pay->rq[il*ndet*ndet+if1*ndet+if2];
        }
      }
    }
  }
 //sprintf(nm,"rq_after_%s.dat",SC->comp_name);
 //write_bin_vector(rq, nm, sizeof(double)*(p_pay->nbins*p_pay->p_model->ndet*p_pay->p_model->ndet), err);   
 //forwardError(*err,__LINE__,);
    
}
void free_comp_parametric(void** data) {
  SmicaComp *SC;
  parametric_smica *p_pay;
  
  SC = *data;
  p_pay = SC->data;
  free(p_pay->rq);
  if (p_pay->nbins!=0) {
    free(p_pay->bins);
    free(p_pay->bi);
    free(p_pay->bo);
    free(p_pay->wbins);
  }
  if (p_pay->wl!=NULL) {
    free(p_pay->wl);
  }
  
  parametric_free((void**)&(p_pay->p_model));
  free(p_pay->A);

  free(SC->data);
  free(SC);
  *data = NULL;
}

void apply_rename(int nkey, char* keys, int nrename, char* rename_from, char* rename_to) {
  int i,j;

  if (nrename==0) {
    return;
  }

  for(i=0;i<nkey;i++) {
    for(j=0;j<nrename;j++) {
      if (strcmp(&(keys[i*256]),&(rename_from[j*256]))==0) {
        memcpy(&(keys[i*256]),&(rename_to[j*256]),sizeof(char)*256);
        break;
      }
    }
  }
}

int base_parametric_cldf_init(cldf *df,int m, double** detlist,int *ndef, char ***defkeys, char*** defvalues, int *nvar, char ***varkeys, error **err) {
  int dz,i;
  char *keyvartable,*deftable,*valtable;
  int nrename;
  char *rename_from, *rename_to;
  int hk;
  
  *nvar = cldf_readint(df,"ndim",err);
  forwardError(*err,__LINE__,0);
  
  nrename=0;
  rename_from = NULL;
  rename_to = NULL;

  hk = cldf_haskey(df,"nrename",err);
  forwardError(*err,__LINE__,0);  
  if (hk ==1) {
    nrename = cldf_readint(df,"nrename",err);
    forwardError(*err,__LINE__,0);   
    dz = -1;
    rename_from = cldf_readstr(df,"rename_from",&dz, err);
    forwardError(*err,__LINE__,0);
    dz = -1;
    rename_to = cldf_readstr(df,"rename_to",&dz, err);
    forwardError(*err,__LINE__,0);
  }
  
  dz = -1;
  keyvartable = cldf_readstr(df,"keys",&dz, err);
  forwardError(*err,__LINE__,0);
    
  if (*nvar!=0) {
    *varkeys = malloc_err(sizeof(char*)**nvar,err);
    forwardError(*err,__LINE__,0);
  } else {
    *varkeys = malloc_err(sizeof(char*)*1,err);
    forwardError(*err,__LINE__,0);
    (*varkeys)[0] = NULL;
  }
  
  apply_rename(*nvar,keyvartable,nrename,rename_from,rename_to);
  for(i=0;i<*nvar;i++) {
    (*varkeys)[i] = &(keyvartable[i*256]);
  }
  
  // get defaults
  *ndef = cldf_readint(df,"ndef",err);
  forwardError(*err,__LINE__,0);
  
  dz = -1;
  deftable = cldf_readstr(df,"defaults",&dz, err);
  forwardError(*err,__LINE__,0);
  dz = -1;
  valtable = cldf_readstr(df,"values",&dz, err);
  forwardError(*err,__LINE__,0);
  
  if (*ndef!=0) {
    *defkeys = malloc_err(sizeof(char*)**ndef,err);
    forwardError(*err,__LINE__,0);
    *defvalues = malloc_err(sizeof(char*)**ndef,err);
    forwardError(*err,__LINE__,0);    
  } else {
    *defkeys = malloc_err(sizeof(char*)*1,err);
    forwardError(*err,__LINE__,0);
    *defvalues = malloc_err(sizeof(char*)*1,err);
    forwardError(*err,__LINE__,0);
    (*defkeys)[0] = NULL;
    (*defvalues)[0] = NULL;
  }
  
  apply_rename(*ndef,deftable,nrename,rename_from,rename_to);
  for(i=0;i<*ndef;i++) {
    (*defkeys)[i] = &(deftable[i*256]);
    (*defvalues)[i] = &(valtable[i*256]);
  }  
  
  dz = -1;
  hk = cldf_haskey(df,"dfreq",err);
  forwardError(*err,__LINE__,0);
  if (hk ==1) {
    *detlist = cldf_readfloatarray(df,"dfreq",&dz, err);
    forwardError(*err,__LINE__,0);
  } else {
    int * ietlist;
    ietlist = cldf_readintarray(df,"freq",&dz, err);
    forwardError(*err,__LINE__,0);
    *detlist = malloc_err(sizeof(double)*dz,err);
    forwardError(*err,__LINE__,0);
    
    for(i=0;i<dz;i++) {
      (*detlist)[i]=ietlist[i];
    }
    free(ietlist);
  }
  

  free(rename_from);
  free(rename_to);

  return dz;

}


SmicaComp * finalize_parametric_cldf_init(parametric* p_model,cldf *df,int nb, int m, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins,error **err) {
  parametric_smica *p_pay;
  int i,eb;
  char **xnames;
  SmicaComp *SC;
  int lmin,lmax;
  double *color;
  int dz;
  int nrename;
  char *rename_from, *rename_to;
  int j;
  int nvoid;
  int *voidlist;
  int hk;
  int ncl;

  nrename=0;
  rename_from = NULL;
  rename_to = NULL;

  hk = cldf_haskey(df,"nrename",err);
  forwardError(*err,__LINE__,NULL);  
  if (hk ==1) {
    nrename = cldf_readint(df,"nrename",err);
    forwardError(*err,__LINE__,NULL);   
    dz = -1;
    rename_from = cldf_readstr(df,"rename_from",&dz, err);
    forwardError(*err,__LINE__,NULL);
    dz = -1;
    rename_to = cldf_readstr(df,"rename_to",&dz, err);
    forwardError(*err,__LINE__,NULL);
  }

  hk = cldf_haskey(df,"color",err);
  forwardError(*err,__LINE__,NULL);  
  if (hk ==1) {
    dz  = -1;
    color = cldf_readfloatarray(df,"color",&dz, err);
    forwardError(*err,__LINE__,NULL);
    if (dz == p_model->ndet) {
      double *rolor;
      int ii,jj;
      dz = p_model->ndet * p_model->ndet;
      rolor = malloc_err(sizeof(double)*dz,err);
      forwardError(*err,__LINE__,NULL);
      for (ii=0;ii<p_model->ndet;ii++) {
        for (jj=ii;jj<p_model->ndet;jj++) {
          rolor[ii*p_model->ndet+jj] = color[ii]*color[jj];
          rolor[jj*p_model->ndet+ii] = color[ii]*color[jj];
        }
      }
      free(color);
      color = rolor;
    }
    testErrorRet(dz != p_model->ndet*p_model->ndet,-12443243,"bad size of color array",*err,__LINE__,NULL);
    parametric_set_color(p_model,color,err);
    forwardError(*err,__LINE__,NULL);
    free(color);
  }
  
  hk = cldf_haskey(df,"nvoid",err);
  forwardError(*err,__LINE__,NULL);  
  if (hk ==1) {
    nvoid = cldf_readint(df,"nvoid",err);
    forwardError(*err,__LINE__,NULL);   
    dz = nvoid;
    voidlist = cldf_readintarray(df,"voidlist",&dz, err);
    forwardError(*err,__LINE__,NULL);
    parametric_set_void(p_model,nvoid,voidlist,err);
    forwardError(*err,__LINE__,NULL);
    free(voidlist);
  }
      
  lmin = ell[0];

  lmax = ell[nell-1];
  testErrorRet(nell!=(lmax-lmin+1),-111,"SAFEGARD",*err,__LINE__,NULL);

  /*eb = 0;
  for(i=1;i<6;i++) {
    eb +=has_cl[i];
  }
  testErrorRet(eb!=0,-7693,"parametric does not work with polarized data yet",*err,__LINE__,NULL);
  */

  p_pay = malloc_err(sizeof(parametric_smica),err);
  forwardError(*err,__LINE__,NULL);
    
  p_pay->m = m;
  p_pay->p_model = p_model;  
  p_pay->unit = unit;
  
  p_pay->A = cldf_readfloatarray(df,"A_cmb",&m,err);
  forwardError(*err,__LINE__,NULL);    

  p_pay->nell = nell;

  p_pay->nbins = 0;
  p_pay->bins = NULL;
  if (bins !=NULL) {
    int li,bi,bn;
    //ncl = p_model->has_TEB[0] + p_model->has_TEB[1] + p_model->has_TEB[2] + p_model->has_TEB[0] * p_model->has_TEB[1] + p_model->has_TEB[0] * p_model->has_TEB[2] +p_model->has_TEB[1] * p_model->has_TEB[2];
    ncl = has_cl[0] + has_cl[1] + has_cl[2] + has_cl[3] + has_cl[4] + has_cl[5]; 
    //_DEBUGHERE_("ncl %d",ncl);
    p_pay->nbins = nbins;
    p_pay->bins = malloc_err(sizeof(double)*(ncl*nell*nbins),err);
    forwardError(*err,__LINE__,NULL);
    memcpy(p_pay->bins,bins,sizeof(double)*nbins*ncl*nell);    
    p_pay->bi = malloc_err(sizeof(int)*nbins,err);
    forwardError(*err,__LINE__,NULL);
    p_pay->bo = malloc_err(sizeof(int)*nbins,err);
    forwardError(*err,__LINE__,NULL);
    p_pay->wbins = malloc_err(sizeof(double)*nell,err);
    bn = 0;
    for(bi=0;bi<nbins;bi++) {
      for(li = 0;li<nell;li++) {
        p_pay->bi[bi] = 0;
        if (p_pay->bins[bi*ncl*nell+li]!=0) {
          p_pay->bi[bi] = li;
          break;
        }
      }
      for(li=p_pay->bi[bi];li<nell;li++) {
        p_pay->bo[bi] = nell;
        if (p_pay->bins[bi*ncl*nell+li]==0) {
          p_pay->bo[bi] = li;
          break; 
        }
        p_pay->wbins[bn] = p_pay->bins[bi*ncl*nell+li];
        bn++;
      }
      //memcpy(&(p_pay->wbins[bn]),&(p_pay->bins[bi*nell+ p_pay->bi[bi]]),sizeof(double)*( p_pay->bo[bi]- p_pay->bi[bi])); 
      //bn += p_pay->bo[bi]- p_pay->bi[bi];
    }
    //write_bin_vector(bins, "bins.dat", sizeof(double)*(nell*nbins), err);  
    //forwardError(*err,__LINE__,-1); 
    //write_bin_vector(p_pay->bi, "bi.dat", sizeof(int)*(nbins), err);  
    //forwardError(*err,__LINE__,-1); 
    //write_bin_vector(p_pay->bo, "bo.dat", sizeof(int)*(nbins), err);  
    //forwardError(*err,__LINE__,-1); 
    //_DEBUGHERE_("%d",bn);
    //write_bin_vector(p_pay->wbins, "wb.dat", sizeof(double)*(bn), err);  
    //forwardError(*err,__LINE__,-1); 

  }
  p_pay->wl = NULL;
  if (wl!=NULL) {
    p_pay->wl = malloc_err(sizeof(double)*(nell),err);
    forwardError(*err,__LINE__,NULL);
    memcpy(p_pay->wl,wl,sizeof(double)*nell);    
    
  }
  p_pay->rq = malloc_err(sizeof(double)*(lmax+1-lmin)*m*m,err);
  forwardError(*err,__LINE__,NULL);
  memset(p_pay->rq,0,sizeof(double)*(lmax+1-lmin)*m*m);
  
  SC = alloc_SC(p_model->nvar,nb,m,p_pay,&comp_parametric_update,&free_comp_parametric,err);
  forwardError(*err,__LINE__,NULL);
  
  SC_isfg(SC);  
  

  if (p_model->nvar!=0) {
    xnames = malloc_err(sizeof(char*)*(p_model->nvar),err);
    forwardError(*err,__LINE__,NULL);
  } else{
    xnames = malloc_err(sizeof(char*)*1,err);
    forwardError(*err,__LINE__,NULL);
  }   
  for(i=0;i<p_model->nvar;i++) {
    char *kp;
    kp = p_model->varkey[i];
    for(j=0;j<nrename;j++) {
      if (strcmp(kp,&(rename_to[j*256]))==0) {
        kp = &(rename_from[j*256]);
        break;
      }
    }
    xnames[i] = kp;
  }
  
  SC_setnames(SC, xnames, err);
  forwardError(*err,__LINE__,NULL);
  
  free(xnames);
  free(rename_from);
  free(rename_to);
  
  return SC;    
}

//CREATE_PARAMETRIC_FILE_INIT(powerlaw_free_emissivity,powerlaw_free_emissivity_init);

