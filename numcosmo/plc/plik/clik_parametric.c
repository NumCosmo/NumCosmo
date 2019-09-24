#include "clik_parametric.h"


void get_freq(int ndet, double *detlist, int* pnfreq, double** pfreqlist, int** pdet2freq, error **err) {
  int i,j;
  double *freqlist;
  int *det2freq;
  int nfreq;
  int mdet;

  if (ndet==0) {
    *pnfreq = 0;
    *pfreqlist = NULL;
    *pdet2freq = NULL;
    return;
  }
  mdet = ndet;
  if (ndet<0) {
    mdet=-ndet;
  }

  freqlist = malloc_err(sizeof(double)*mdet,err);
  forwardError(*err,__LINE__,);
  
  det2freq = malloc_err(sizeof(int)*mdet, err);
  forwardError(*err,__LINE__,);
  
  if (ndet<0) {
    nfreq = mdet;
    for(i=0;i<mdet;i++) {
      det2freq[i] = i;
      freqlist[i]=detlist[i];
    }
  }
  else {
    nfreq = 0;
    for(i=0;i<mdet;i++) {
      det2freq[i]=-1;
      for(j=0;j<nfreq;j++) {
        if (fabs(freqlist[j]-detlist[i])<1e-6) {
          det2freq[i] = j;
          break;
        }
      }
      if (det2freq[i]==-1) {
        det2freq[i]=nfreq;
        freqlist[nfreq]=detlist[i];
        //_DEBUGHERE_("f %g",detlist[i]);
        nfreq++;
      }
    }  
  }
  
  *pnfreq = nfreq;
  *pfreqlist = freqlist;
  *pdet2freq = det2freq;
}

parametric *parametric_init(int ndet,double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *epl;
  int has_TEB[3];

  has_TEB[0] = 1;
  has_TEB[1] = 0;
  has_TEB[2] = 0;

  epl = parametric_pol_init(ndet,0,has_TEB,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);
  
  return epl;
}

parametric *parametric_pol_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *epl;
  int i,jj,mnvar;
  char *nop;

  epl = malloc_err(sizeof(parametric),err);
  forwardError(*err,__LINE__,NULL);
  
  epl->pf = pflist_init(err);
  forwardError(*err,__LINE__,NULL);
  
  epl->default_settings = pflist_init(err);
  forwardError(*err,__LINE__,NULL);
  
  nop = NULL;
  for(i=0;i<nvar;i++) {
    pflist_add_item(epl->pf,1,&(varkey[i]),&nop,err);
    forwardError(*err,__LINE__,NULL);
    //_DEBUGHERE_("%s",varkey[i]);
  }

  pflist_add_item(epl->pf,ndef,defkey,defvalue,err);
  forwardError(*err,__LINE__,NULL);

  epl->nvar = nvar;
  epl->ndef = epl->pf->nkey - nvar;

  epl->ndet_T = ndet_T;
  if (ndet_T<0) {
    epl->ndet_T = -ndet_T;
  }
  get_freq(ndet_T, detlist, &(epl->nfreq_T), &(epl->freqlist_T), &(epl->det2freq_T),err);
  forwardError(*err,__LINE__,NULL);
  
  epl->ndet_P = ndet_P;
  if (ndet_P<0) {
    epl->ndet_P = -ndet_P;
  }
  get_freq(ndet_P, detlist + epl->ndet_T, &(epl->nfreq_P), &(epl->freqlist_P), &(epl->det2freq_P),err);
  forwardError(*err,__LINE__,NULL);

  epl->detlist_T = NULL;
  if (epl->ndet_T!=0) {
    epl->detlist_T = malloc_err(sizeof(double)*epl->ndet_T,err);
    forwardError(*err,__LINE__,NULL);
    memcpy(epl->detlist_T,detlist,sizeof(double)*epl->ndet_T);
  }
  
  epl->detlist_P = NULL;
  if (epl->ndet_P!=0) {
    epl->detlist_P = malloc_err(sizeof(double)*epl->ndet_P,err);
    forwardError(*err,__LINE__,NULL);
    memcpy(epl->detlist_P,detlist+epl->ndet_T,sizeof(double)*epl->ndet_P);
  }
  
  epl->payload = NULL;
  epl->eg_compute = NULL;
  epl->eg_free = NULL;

  epl->lmin = lmin;
  epl->lmax = lmax;

  epl->has_TEB[0] =  has_TEB[0];
  epl->has_TEB[1] =  has_TEB[1];
  epl->has_TEB[2] =  has_TEB[2];

  epl->ndet = epl->ndet_T + epl->ndet_P * (epl->has_TEB[1]+epl->has_TEB[2]);
  epl->nfreq = epl->nfreq_T + epl->nfreq_P * (epl->has_TEB[1]+epl->has_TEB[2]);

  epl->detlist = malloc_err(sizeof(double)*epl->ndet,err);
  forwardError(*err,__LINE__,NULL);
  epl->freqlist = malloc_err(sizeof(double)*epl->nfreq,err);
  forwardError(*err,__LINE__,NULL);
  epl->det2freq = malloc_err(sizeof(double)*epl->ndet,err);
  forwardError(*err,__LINE__,NULL);

  for(i=0;i<epl->ndet_T;i++) {
    epl->detlist[i] = epl->detlist_T[i];
    epl->det2freq[i] = epl->det2freq_T[i];
  }
  for(i=0;i<epl->nfreq_T;i++) {
    epl->freqlist[i] = epl->freqlist_T[i];
  }

  for(i=0;i<epl->ndet_P*epl->has_TEB[1];i++) {
    epl->detlist[i+epl->ndet_T] = epl->detlist_P[i];
    epl->det2freq[i+epl->ndet_T] = epl->det2freq_P[i]+epl->ndet_T;
  }
  for(i=0;i<epl->nfreq_P*epl->has_TEB[1];i++) {
    epl->freqlist[i+epl->ndet_T] = epl->freqlist_P[i];
  }

  for(i=0;i<epl->ndet_P*epl->has_TEB[2];i++) {
    epl->detlist[i+epl->ndet_T+epl->ndet_P*epl->has_TEB[1]] = epl->detlist_P[i];
    epl->det2freq[i+epl->ndet_T+epl->ndet_P*epl->has_TEB[1]] = epl->det2freq_P[i]+epl->ndet_T+epl->ndet_P*epl->has_TEB[1];
  }
  for(i=0;i<epl->nfreq_P*epl->has_TEB[2];i++) {
    epl->freqlist[i+epl->ndet_T+epl->ndet_P*epl->has_TEB[1]] = epl->freqlist_P[i];
  }

  epl->sRq = malloc_err(sizeof(double)*(lmax+1-lmin)*epl->nfreq*epl->nfreq,err);
  forwardError(*err,__LINE__,NULL);
  
  if (epl->nvar>0) {
    epl->sdRq = malloc_err(sizeof(double)*(lmax+1-lmin)*epl->nfreq*epl->nfreq*epl->nvar,err);
    forwardError(*err,__LINE__,NULL);
  } else {
      epl->sdRq = malloc_err(sizeof(double)*(lmax+1-lmin)*epl->nfreq*epl->nfreq*1,err);
      forwardError(*err,__LINE__,NULL);
  }
  
  mnvar = 1;
  if (nvar!=0) {
    mnvar = nvar;
  }
  epl->varkey = malloc_err(sizeof(pfchar)*mnvar,err);
  forwardError(*err,__LINE__,NULL);
  epl->ikey = malloc_err(sizeof(int)*mnvar,err);
  forwardError(*err,__LINE__,NULL);
  
  
  for(jj=0;jj<epl->nvar;jj++) {
    int ps;
    sprintf(epl->varkey[jj],"%s",varkey[jj]);
    ps = pflist_key_index(epl->pf,varkey[jj],err);
    forwardError(*err,__LINE__,NULL);
    testErrorRet(ps==-1,-10010,"Fatal",*err,__LINE__,NULL);
    epl->ikey[jj]=ps;
    //dumpHERE_("%d %s",jj,epl->varkey[jj]);
  }


  epl->dnofail=0;

  epl->color = malloc_err(sizeof(double)*epl->ndet*epl->ndet,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<epl->ndet*epl->ndet;i++) {
    epl->color[i]=1;
  }  
  epl->has_color = 0;

  epl->eg_deriv_any = NULL;
  epl->eg_deriv = NULL;
  epl->nderiv = 0;
  epl->eg_compute_and_deriv = NULL;

  epl->nvoid = 0;
  epl->voidlist = NULL;

  return epl;
}


parametric *parametric_bydet_pol_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric * egl;

  egl = parametric_pol_init(-ndet_T, -ndet_P,has_TEB,detlist,ndef,defkey,defvalue, nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}

parametric *parametric_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric * egl;

  egl = parametric_init(-ndet,detlist,ndef,defkey,defvalue, nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}

void parametric_set_color(parametric *egl,double *color, error **err) {
  int i;
  memcpy(egl->color,color,sizeof(double)*egl->ndet*egl->ndet);
  for (i=0;i<egl->ndet*egl->ndet;i++) {
    if (egl->color[i]!=1) {
      egl->has_color = 1;
      break;
    }
  }
}

void parametric_add_derivative_function(parametric *egl, char* vk, exg_deriv* fnc,error **err) {
  if (strcmp(vk,"any")==0) {
    egl->eg_deriv_any = fnc;
    return;
  }
  if (egl->nderiv == 0) {
    egl->eg_deriv = malloc_err(sizeof(exg_compute*)*100,err);
    egl->deriv_key = malloc_err(sizeof(pfchar)*100,err);
  }
  testErrorRet(egl->nderiv+1>100,-1234532,"too many derivative functions",*err,__LINE__,);
  egl->eg_deriv[egl->nderiv] = fnc;
  sprintf(egl->deriv_key[egl->nderiv],"%s",vk);
  egl->nderiv ++;
}


void parametric_dnofail(parametric* egl, int vl) {
  egl->dnofail = vl;
}

void parametric_free(void** pegl) {
  parametric *egl;

  
  egl = *pegl;
  if(egl->eg_free!=NULL) {
    //_DEBUGHERE_("","");
    egl->eg_free(&(egl->payload));
  }
  //_DEBUGHERE_("","");
  free(egl->varkey); 
  free(egl->ikey); 
  //_DEBUGHERE_("","");
  if (egl->ndet_T!=0) {
    free(egl->det2freq_T);  
    free(egl->freqlist_T);
    free(egl->detlist_T);
  }
  //_DEBUGHERE_("","");
  if (egl->ndet_P!=0) {
    free(egl->det2freq_P);
    free(egl->freqlist_P);
    free(egl->detlist_P);
  }
  //_DEBUGHERE_("","");
  free(egl->det2freq);
  free(egl->freqlist);
  free(egl->detlist);
  free(egl->color);
  free(egl->sdRq);
  free(egl->sRq);
  pflist_free((void**)&(egl->pf));
  pflist_free((void**)&(egl->default_settings));
  if (egl->nderiv!=0) {
    //_DEBUGHERE_("","");
    free(egl->eg_deriv);
    free(egl->deriv_key);
  }
  if (egl->nvoid!=0) {
    //_DEBUGHERE_("","");
    free(egl->voidlist);
  }
  //_DEBUGHERE_("","");
  free(egl);
  *pegl = NULL;
}

void parametric_set_void(parametric *egl, int nvoid, int *voidlist,error **err) {
  int i;
  
  egl->nvoid = nvoid;
  egl->voidlist = malloc_err(sizeof(int)*nvoid,err);
  forwardError(*err,__LINE__,);
  for(i=0;i<nvoid;i++) {
    testErrorRet(voidlist[i]<0 || voidlist[i]>=egl->ndet,-3525342,"bad void list",*err,__LINE__,);
    egl->voidlist[i]=voidlist[i];
  }
}

void parametric_compute(parametric *egl, double *pars, double* Rq, double *dRq, error **err) {
  int idet,jdet;
  int ell,dvar;
  int midet, mjdet;
  int m2,f2,mell,fell,ivar,mdvar,fdvar; 
  int dvaro;
  double *sRq, *sdRq;
  double v;

  
  for(ivar=0;ivar<egl->nvar;ivar++) {
    int ps;
    //_DEBUGHERE_("%s %d %g",egl->pf->key[egl->ikey[ivar]],egl->ikey[ivar],pars[ivar]);
    sprintf(egl->pf->value[egl->ikey[ivar]],"%40g",pars[ivar]);
  }
  
  if (egl->ndet!=egl->nfreq) {  
    sRq = NULL;
    if (Rq!=NULL) {
      sRq = egl->sRq;
    }
    sdRq = NULL;
    if (dRq!=NULL) {
      sdRq = egl->sdRq;
    }
    parametric_compute_loop(egl,sRq,sdRq,err);
    forwardError(*err,__LINE__,);
  
    m2 = egl->ndet;
    m2 = m2*m2;
    f2 = egl->nfreq;
    f2 = f2*f2;
    mdvar = (egl->lmax+1-egl->lmin) * m2;
    fdvar = (egl->lmax+1-egl->lmin) * f2;
  
    if (sRq!=NULL) {
      for(ell=0;ell<egl->lmax+1-egl->lmin;ell++) {
        mell = ell*m2;
        fell = ell*f2;
        for(idet=0;idet<egl->ndet;idet++) {
          midet = egl->det2freq[idet];
          for(jdet=idet;jdet<egl->ndet;jdet++) {
            mjdet = egl->det2freq[jdet];
            //_DEBUGHERE_("%d->%d %d->%d %d %d",idet,midet,jdet,mjdet,egl->ndet,egl->nfreq);
            v = sRq[fell + midet*egl->nfreq + mjdet] * egl->color[idet*egl->ndet+jdet];
            Rq[mell + idet*egl->ndet + jdet] = v;
            Rq[mell + jdet*egl->ndet + idet] = v;
          }
        }
      }      
    }
    if (sdRq!=NULL) {
      for(ell=0;ell<egl->lmax+1-egl->lmin;ell++) {
        mell = ell*m2;
        fell = ell*f2; 
        for(idet=0;idet<egl->ndet;idet++) {
          midet = egl->det2freq[idet];
          for(jdet=idet;jdet<egl->ndet;jdet++) {
            mjdet = egl->det2freq[jdet];
            for(dvar=0;dvar<egl->nvar;dvar++) {

              v = sdRq[dvar*fdvar + fell + midet*egl->nfreq + mjdet] * egl->color[idet*egl->ndet+jdet];
              dRq[dvar*mdvar + mell + idet*egl->ndet + jdet] = v;
              dRq[dvar*mdvar + mell + jdet*egl->ndet + idet] = v;
            }
          }
        }
      }      
    }
  } else {
    parametric_compute_loop(egl,Rq,dRq,err);
    forwardError(*err,__LINE__,);
    if (egl->has_color==1) {
      m2 = egl->ndet;
      m2 = m2*m2;
      f2 = egl->nfreq;
      f2 = f2*f2;
      mdvar = (egl->lmax+1-egl->lmin) * m2;
      fdvar = (egl->lmax+1-egl->lmin) * f2;
      
      for(ell=0;ell<egl->lmax+1-egl->lmin;ell++) {
        mell = ell*m2;
        for(idet=0;idet<egl->ndet;idet++) {
          for(jdet=idet;jdet<egl->ndet;jdet++) {
            //_DEBUGHERE_("%d->%d %d->%d %d %d",idet,midet,jdet,mjdet,egl->ndet,egl->nfreq);
            v = Rq[mell + idet*egl->ndet + jdet] * egl->color[idet*egl->ndet+jdet];
            Rq[mell + idet*egl->ndet + jdet] = v;
            Rq[mell + jdet*egl->ndet + idet] = v;
          }
        }
      }
      if (dRq!=NULL) {
        for(ell=0;ell<egl->lmax+1-egl->lmin;ell++) {
          mell = ell*m2;
          fell = ell*f2; 
          for(idet=0;idet<egl->ndet;idet++) {
            for(jdet=idet;jdet<egl->ndet;jdet++) {
              for(dvar=0;dvar<egl->nvar;dvar++) {
                v = dRq[dvar*mdvar + mell + idet*egl->ndet + jdet] * egl->color[idet*egl->ndet+jdet];
                dRq[dvar*mdvar + mell + idet*egl->ndet + jdet] = v;
                dRq[dvar*mdvar + mell + jdet*egl->ndet + idet] = v;
              }
            }
          }
        }
      }
    }
    
  }
  if (egl->nvoid!=0) {
  
    for(ell=0;ell<egl->lmax+1-egl->lmin;ell++) {
      int ii;
      mell = ell*m2;
      for(ii=0;ii<egl->nvoid;ii++) {
        idet = egl->voidlist[ii];
        for(jdet=0;jdet<egl->ndet;jdet++) {
          Rq[mell+idet*egl->ndet + jdet] = 0;
          Rq[mell+jdet*egl->ndet + idet] = 0;
        }
      }
    }
    if (dRq!=NULL) {
      for(dvar=0;dvar<egl->nvar;dvar++) {
        int pvar;
        pvar = dvar*mdvar;
        for(ell=0;ell<egl->lmax+1-egl->lmin;ell++) {
          int ii;
          mell = ell*m2;
          for(ii=0;ii<egl->nvoid;ii++) {
            idet = egl->voidlist[ii];
            for(jdet=0;jdet<egl->ndet;jdet++) {
              Rq[pvar + mell+idet*egl->ndet + jdet] = 0;
              Rq[pvar + mell+jdet*egl->ndet + idet] = 0;
            }
          }
        } 
      }
    }
  }
  
}

void parametric_compute_loop(parametric * egl, double* Rq, double *dRq, error **err) {
  int ok,i,j;

  testErrorRet(egl->eg_compute==NULL && egl->eg_compute_and_deriv==NULL ,-1234,"badly initialized",*err,__LINE__,);
  testErrorRet(egl->eg_compute_and_deriv==NULL && egl->nderiv==0 && egl->eg_deriv_any==NULL && dRq!=NULL,-1234,"Don't know how to compute derivatives",*err,__LINE__,);
  
  
  if (egl->eg_compute_and_deriv!=NULL) {
    egl->eg_compute_and_deriv(egl, Rq, dRq, err);
    forwardError(*err,__LINE__,);
    return;
  }
  
  // first use the compute function
  egl->eg_compute(egl, Rq, err);
  forwardError(*err,__LINE__,);
  
  // then look whether we need derivatives 
  if (dRq==NULL) {
    return;
  } 
  
  
  for(i=0;i<egl->nvar;i++) {
    for(j=0;j<egl->nderiv;j++) {
      if (strcmp(egl->varkey[i],egl->deriv_key[j])==0) {
        egl->eg_deriv[j](egl,i,Rq,PTR_DER(egl,i,dRq),err);
        forwardError(*err,__LINE__,);
        ok = 1;
        break;      
      }
    }
    if (ok==1) {
      continue;
    }
    
    if (egl->eg_deriv_any!=NULL) {
      egl->eg_deriv_any(egl,i,Rq,PTR_DER(egl,i,dRq),err);
      forwardError(*err,__LINE__,);
      continue;
    }
        
    parametric_end_derivative_loop(egl,dRq,egl->varkey[i],err);
    forwardError(*err,__LINE__,);
  }
  return;
}

void parametric_end_derivative_loop(parametric *egl,double* dRq, char* varkey, error **err) {
  testErrorRetVA(egl->dnofail!=1,-1234,"Cannot derive on parameter '%s'",*err,__LINE__,,varkey);
  memset(dRq,0,sizeof(double)*(egl->lmax+1-egl->lmin)*egl->nfreq*egl->nfreq);
}

double parametric_get_default(parametric* egl,char *key, error **err) {
  char *res;
  double vres;
  res = pflist_get_value(egl->default_settings,key,NULL,err);
  forwardError(*err,__LINE__,0);
  testErrorRetVA(res==NULL,-123432,"unknown default setting '%s'",*err,__LINE__,-1,key);
  //_DEBUGHERE_("'%s' -> '%s'",key,res);
  vres = atof(res);
  //_DEBUGHERE_("%g",vres);
  return vres;
}


void parametric_set_default(parametric* egl,char *key, double value,error **err) {
  char calue[200];
  char *nalue;
  int ps;
  char nop;
  char *res;

  //_DEBUGHERE_("%s %g",key,value);
  sprintf(calue,"%30g",value);
  //_DEBUGHERE_("%p",egl->default_settings);
  //_DEBUGHERE_("%s",calue);
  nalue = calue;
  //_DEBUGHERE_("%p %p",&(nalue),calue);
  
  //has it been defined already ?
  nop = '\0';
  res = pflist_get_value(egl->pf,key,&nop,err);
  forwardError(*err,__LINE__,);
  //_DEBUGHERE_("'%s' %d",res,res[0]);
  if (res[0]!='\0') {
    //_DEBUGHERE_("","");
    sprintf(calue,"%s",res);
  }
  //_DEBUGHERE_("value %s",nalue);

  pflist_add_item(egl->default_settings,1, &key, (char**)&(nalue),err);
  forwardError(*err,__LINE__,);
  //_DEBUGHERE_("","");
  
}

double parametric_get_value(parametric *egl, char *key, error **err) {
  int ps;
  double res;

  ps=pflist_key_index(egl->pf,key,err);
  forwardError(*err,__LINE__,0);
  if (ps==-1) { // not in the default/var try in the default_settings
    res = parametric_get_default(egl,key,err);
    forwardError(*err,__LINE__,0);
    return res;
  }
  res = 0;
  res = pflist_get_double_value(egl->pf,key,&res,err);
  forwardError(*err,__LINE__,-1);
  return res;
}

void parametric_declare_mandatory(parametric *egl, char* key, error **err) {
  int ps;

  ps = pflist_key_index(egl->pf,key,err);
  forwardError(*err,__LINE__,);

  testErrorRetVA(ps==-1,-1234332,"Mandatory parameter '%s' is absent from the list of variable or default parameters",*err,__LINE__,,key);
  return;
}

void parametric_norm_derivative(parametric * egl, int iv, double *Rq, double *dRq, error **err) {
  double norm;
  int ell,m1,m2;
  char *key;

  key = egl->varkey[iv];
  norm = parametric_get_value(egl,key,err);
  forwardError(*err,__LINE__,);

  if (norm==0) {
    sprintf(egl->pf->value[iv],"%40g",1.);
    egl->eg_compute(egl,dRq,err);
    forwardError(*err,__LINE__,);
    sprintf(egl->pf->value[iv],"%40g",0.);
    return;
  }

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<egl->nfreq;m1++) {
      for (m2=m1;m2<egl->nfreq;m2++) {
        dRq[IDX_R(egl,ell,m1,m2)] = Rq[IDX_R(egl,ell,m1,m2)]/norm;
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  return;
}


void parametric_triangle_fill(parametric *egl, double *A, error **err) {
  int m1,m2,p;
  pfchar Ac;
  double *Abuf;

  Abuf = A + egl->nfreq*egl->nfreq;
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      sprintf(Ac,"%s%d_%d",egl->tensor_norm_template,m1,m2);    
      Abuf[m1*egl->nfreq+m2] = parametric_get_value(egl,Ac,err);
      forwardError(*err,__LINE__,);  
      //_DEBUGHERE_("%d %d %g",m1,m2,Abuf[m1*egl->nfreq+m2])
    }
  }
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      A[m1*egl->nfreq+m2] = Abuf[m1*egl->nfreq+m2] * Abuf[m2*egl->nfreq+m2];
      //_DEBUGHERE_("%d %d %g %g %g",m1,m2,Abuf[m1*egl->nfreq],Abuf[m2*egl->nfreq],A[m1*egl->nfreq+m2])

      for(p=m2+1;p<egl->nfreq;p++) {
        A[m1*egl->nfreq+m2] += Abuf[m1*egl->nfreq+p] * Abuf[m2*egl->nfreq+p];
        //_DEBUGHERE_("%d %d %d %g %g %g",m1,m2,p,Abuf[m1*egl->nfreq+p],Abuf[m2*egl->nfreq+p],A[m1*egl->nfreq+m2])
      }
      A[m2*egl->nfreq+m1] = A[m1*egl->nfreq+m2];

    }
  }
}


void parametric_triangle_fill_derivative(parametric * egl, int iv, double *A, error **err) {
  double norm;
  int ell,m1,m2;
  char *key,*kp;
  int ic1,ic2,p;
  double *mRq;

  key = egl->varkey[iv];
  
  testErrorRet(strncmp(key,egl->tensor_norm_template,egl->tensor_norm_template_len)!=0,-141241,"Argl",*err,__LINE__,);
  sscanf(&(key[egl->tensor_norm_template_len]),"%d_%d",&ic1,&ic2);

  memset(A,0,sizeof(double)*egl->nfreq*egl->nfreq);

  for(p=0;p<=ic2;p++) {
    pfchar Ac;
    sprintf(Ac,"%s%d_%d",egl->tensor_norm_template,p,ic2);    
    A[ic1*egl->nfreq+p] += parametric_get_value(egl,Ac,err);
    forwardError(*err,__LINE__,); 
    A[p*egl->nfreq+ic1] += A[ic1*egl->nfreq+p];
    //_DEBUGHERE_("%s %d %d -> %d %d %g %d %d %g",key,ic1,ic2,p,ic1, A[p*egl->nfreq+ic1] ,ic1,p, A[ic1*egl->nfreq+p]);
  }
}

void parametric_tanh_fill(parametric *egl, double *A, error **err) {
  int m1,m2,p;
  pfchar Ac;
  double *Abuf;
  double theta;

  Abuf = A + egl->nfreq*egl->nfreq;
  for(m1=0;m1<egl->nfreq;m1++) {
    sprintf(Ac,"%s%d",egl->tensor_norm_template,m1);    
    Abuf[m1] = parametric_get_value(egl,Ac,err);
    forwardError(*err,__LINE__,); 
  }
  for(m1=0;m1<egl->nfreq;m1++) { 
    A[m1*egl->nfreq+m1] = Abuf[m1]*Abuf[m1];
    for(m2=m1+1;m2<egl->nfreq;m2++) {
      sprintf(Ac,"%sT_%d_%d",egl->tensor_norm_template,m1,m2);    
      theta = parametric_get_value(egl,Ac,err);
      forwardError(*err,__LINE__,); 
      A[m1*egl->nfreq+m2] = Abuf[m1]*Abuf[m2]*tanh(theta);
    }
  }
}


void parametric_tanh_fill_derivative(parametric * egl, int iv, double *A, error **err) {
  double norm;
  int ell,m1,m2;
  char *key,*kp;
  int ic1,ic2,p;
  double *mRq;
  double A1,A2,theta,dm;
  pfchar Ac;
    

  key = egl->varkey[iv];
  memset(A,0,sizeof(double)*egl->nfreq*egl->nfreq);
  
  testErrorRet(strncmp(key,egl->tensor_norm_template,egl->tensor_norm_template_len)!=0,-141241,"Argl",*err,__LINE__,);
  if (key[egl->tensor_norm_template_len]=='T') {
    sscanf(&(key[egl->tensor_norm_template_len]),"T_%d_%d",&ic1,&ic2);
    sprintf(Ac,"%s%d",egl->tensor_norm_template,ic1);    
    A1 = parametric_get_value(egl,Ac,err);
    forwardError(*err,__LINE__,); 
    sprintf(Ac,"%s%d",egl->tensor_norm_template,ic2);    
    A2 = parametric_get_value(egl,Ac,err);
    forwardError(*err,__LINE__,); 
    sprintf(Ac,"%sT_%d_%d",egl->tensor_norm_template,ic1,ic2);    
    theta = parametric_get_value(egl,Ac,err);
    forwardError(*err,__LINE__,);
    A[ic1*egl->nfreq+ic2] = A1*A2 * (1-tanh(theta)*tanh(theta));
    A[ic2*egl->nfreq+ic1] = A[ic1*egl->nfreq+ic2];
  } else {
    sscanf(&(key[egl->tensor_norm_template_len]),"%d",&ic1);
    for(p=0;p<ic1;p++) {
      sprintf(Ac,"%s%d",egl->tensor_norm_template,p);    
      A1 = parametric_get_value(egl,Ac,err);
      forwardError(*err,__LINE__,); 
      sprintf(Ac,"%sT_%d_%d",egl->tensor_norm_template,p,ic1);    
      theta = parametric_get_value(egl,Ac,err);
      forwardError(*err,__LINE__,); 
      A[ic1*egl->nfreq+p] = tanh(theta)*A1;
      A[p*egl->nfreq+ic1] = A[ic1*egl->nfreq+p];
    }
    sprintf(Ac,"%s%d",egl->tensor_norm_template,ic1);    
    A1 = parametric_get_value(egl,Ac,err);
    forwardError(*err,__LINE__,); 
    A[ic1*egl->nfreq+ic1] = 2*A1;

    for(p=ic1+1;p<egl->nfreq;p++) {
      sprintf(Ac,"%s%d",egl->tensor_norm_template,p);    
      A1 = parametric_get_value(egl,Ac,err);
      forwardError(*err,__LINE__,); 
      sprintf(Ac,"%sT_%d_%d",egl->tensor_norm_template,ic1,p);    
      theta = parametric_get_value(egl,Ac,err);
      forwardError(*err,__LINE__,); 
      A[ic1*egl->nfreq+p] = tanh(theta)*A1;
      A[p*egl->nfreq+ic1] = A[ic1*egl->nfreq+p];
    }
  }
}


void parametric_tensor_fill(parametric *egl,double *A,error **err) {
  int ic;
  pfchar Ac;

  for(ic=0;ic<egl->nfreq;ic++) {
    sprintf(Ac,"%s%d",egl->tensor_norm_template,ic);    
    A[ic] = parametric_get_value(egl,Ac,err);
    forwardError(*err,__LINE__,);  
  } 
}

void parametric_tensor_norm_derivative(parametric * egl, int iv, double *Rq, double *dRq, error **err) {
  double norm;
  int ell,m1,m2;
  char *key,*kp;
  int ic;
  double *mRq;

  key = egl->varkey[iv];
  
  testErrorRet(strncmp(key,egl->tensor_norm_template,egl->tensor_norm_template_len)!=0,-141241,"Argl",*err,__LINE__,);
  ic = atoi(&(key[egl->tensor_norm_template_len]));

  norm = parametric_get_value(egl,key,err);
  forwardError(*err,__LINE__,);

  mRq = Rq;

  if (norm==0) {
    sprintf(egl->pf->value[iv],"%40g",1.);
    egl->eg_compute(egl,dRq,err);
    forwardError(*err,__LINE__,);
    sprintf(egl->pf->value[iv],"%40g",0.);
    norm = 1;
    mRq = dRq;
  }

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<egl->nfreq;m1++) {
      if (m1==ic) {
        dRq[IDX_R(egl,ell,m1,m1)] = 2*mRq[IDX_R(egl,ell,m1,m1)]/norm;  
      } else {
        dRq[IDX_R(egl,ell,m1,m1)] = 0;
      }
      
      for (m2=m1+1;m2<egl->nfreq;m2++) {
        if (m2==ic || m1==ic) {
          dRq[IDX_R(egl,ell,m1,m2)] = mRq[IDX_R(egl,ell,m1,m2)]/norm;
        } else {
          dRq[IDX_R(egl,ell,m1,m2)] = 0;
        }
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  return;
}


void parametric_index_derivative(parametric * egl, int iv, double* Rq, double *dRq, error **err) {
  double index,l_pivot;
  int ell,m1,m2;

  l_pivot = egl->l_pivot;

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    for (m1=0;m1<egl->nfreq;m1++) {
      for (m2=m1;m2<egl->nfreq;m2++) {
        dRq[IDX_R(egl,ell,m1,m2)] = Rq[IDX_R(egl,ell,m1,m2)]  * log((double)ell/l_pivot);
        dRq[IDX_R(egl,ell,m2,m1)] = dRq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }
  return;
}

// ASTRO PART //

// To get from intensity to \delta_T (CMB)

double dBdT(double nu, double nu0) {

  double x,x0,ex,ex0,res,res0;
  x0 = nu0/56.78;
  ex0 = exp(x0);
  x = nu/56.78;
  ex = exp(x);
  res0 = pow(x0,4.0)*ex0/((ex0-1.)*(ex0-1.));
  res = pow(x,4.0)*ex/((ex-1.)*(ex-1.));
  return(res/res0);

}

double sz_spectrum(double nu, double nu0) {

  // This gives the SZ spectrum in \delta_T (CMB)
  // normalized at nu0
  double x, x0, res, res0;
  x0 = nu0/56.78;
  x  = nu /56.78;
  res = 2.0-x/2.0/tanh(x/2.0);
  res0 = 2.0-x0/2.0/tanh(x0/2.0);
  return (res/res0);
}

// Simple power law in ell, constant emissivity
void powerlaw_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double l_pivot,index,A,v;
  
  l_pivot = parametric_get_value(egl,"pw_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"pw_index",err);
  forwardError(*err,__LINE__,);

  A = parametric_get_value(egl,"pw_A",err);  
  forwardError(*err,__LINE__,);

  
  nfreq = egl->nfreq;
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = A*pow((double) ell/l_pivot,(double) index);
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = v;
        Rq[IDX_R(egl,ell,m2,m1)] = v;
      }  
    }
  }
  return;
}

parametric *powerlaw_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &powerlaw_compute;
  egl->eg_free = NULL;
  
  // set default settings (those will be overide by defkey/value and varkey/value)
  parametric_set_default(egl,"pw_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"pw_index",0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"pw_index",&parametric_index_derivative,err);  
  forwardError(*err,__LINE__,NULL);
   
  parametric_set_default(egl,"pw_A",1,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"pw_A",&parametric_norm_derivative,err);  
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}

// Simple power law in ell, arbitrary emissivity (including arbitrary cross-correlations, i.e. NOT rank=1 a priori)
// This could be used e.g for CIB Poisson

void powerlaw_free_emissivity_compute(parametric* egl, double *Rq,  error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;
  double nrmit;
  int l2norm;

  l_pivot = parametric_get_value(egl,"pwfe_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  index = parametric_get_value(egl,"pwfe_index",err);
  forwardError(*err,__LINE__,);

  l2norm = parametric_get_value(egl,"pwfe_l2_norm",err);
  forwardError(*err,__LINE__,);

  nrmit = 1;
  if (l2norm==1) {
    nrmit = l_pivot*(l_pivot+1.)/2./M_PI;
  }
  

  A = egl->payload;
  nfreq = egl->nfreq;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"pwfe_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      v = 1;
      v = parametric_get_value(egl,name,err);
      A[m1*nfreq+m2] = v/nrmit;
      A[m2*nfreq+m1] = v/nrmit;
    }
  }


  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/l_pivot,index);
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(m1=0;m1<nfreq;m1++) {
      for(m2=m1;m2<nfreq;m2++) {
        lA = A[m1*nfreq+m2];
        Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
        Rq[IDX_R(egl,ell,m2,m1)] = lA*v;
      }  
    }
  }

  return;
}

void powerlaw_free_emissivity_A_derivative(parametric* egl, int iv,double *Rq, double *dRq, error **err) {
  int ell,m1,m2,mell,nfreq,mv;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;
  double nrmit;
  int l2norm;


  l_pivot = parametric_get_value(egl,"pwfe_l_pivot",err);
  forwardError(*err,__LINE__,);
  
  index = parametric_get_value(egl,"pwfe_index",err);
  forwardError(*err,__LINE__,);

  l2norm = parametric_get_value(egl,"pwfe_l2_norm",err);
  forwardError(*err,__LINE__,);

  nrmit = 1;
  if (l2norm==1) {
    nrmit = l_pivot*(l_pivot+1.)/2./M_PI;
  }

  //A = egl->payload;
  nfreq = egl->nfreq;
  //for(m1=0;m1<nfreq;m1++) {
  //  for(m2=m1;m2<nfreq;m2++) {
  //    sprintf(name,"pwfe_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
  //    v = 1;
  //    v = parametric_get_value(egl,name,err);
  //    A[m1*nfreq+m2] = v;
  //    A[m2*nfreq+m1] = v;
  //  }
  //}

  stop = 0;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"pwfe_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      if (strcmp(egl->varkey[iv],name)==0) {
        stop=1;
        memset(dRq,0,sizeof(double)*(egl->lmax+1-egl->lmin)*nfreq*nfreq);
        for(ell=egl->lmin;ell<=egl->lmax;ell++) {
          v = pow(ell/l_pivot,index)/nrmit;
          dRq[IDX_R(egl,ell,m1,m2)] = v;
          dRq[IDX_R(egl,ell,m2,m1)] = v;
        }
        break;
      }
    }
    if (stop==1) {
      return;
    }
  }      
  // error return
  parametric_end_derivative_loop(egl,&(dRq[iv]),egl->varkey[iv],err);
  forwardError(*err,__LINE__,);

  return;
}


parametric *powerlaw_free_emissivity_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int m1,m2;
  pfchar name;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  egl->eg_compute = &powerlaw_free_emissivity_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"pwfe_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"pwfe_index",0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"pwfe_index",&parametric_index_derivative,err);  
  forwardError(*err,__LINE__,NULL);
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      sprintf(name,"pwfe_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  

  parametric_set_default(egl,"pwfe_XX_l2_norm",0,err);
  forwardError(*err,__LINE__,NULL);

  parametric_add_derivative_function(egl,"any",&powerlaw_free_emissivity_A_derivative,err);  
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void parametric_check_freq(parametric *egl, double* freqlist, int nfreq, error **err) {
  int m1,m2, ok;

  for(m1=0;m1<egl->nfreq;m1++) {
    ok = 0;
    for(m2=0;m2<nfreq;m2++) {
      if (fabs(egl->freqlist[m1]-freqlist[m2])<1e-6) {
        ok=1;
        break;
      }
    }
    testErrorRetVA(ok==0,-1234,"Cannot compute prediction for %g Ghz channel",*err,__LINE__,,egl->freqlist[m1]);
  }
}
void parametric_check_freq_T(parametric *egl, double* freqlist, int nfreq, error **err) {
  int m1,m2, ok;

  for(m1=0;m1<egl->nfreq_T;m1++) {
    ok = 0;
    for(m2=0;m2<nfreq;m2++) {
      if (fabs(egl->freqlist_T[m1]-freqlist[m2])<1e-6) {
        ok=1;
        break;
      }
    }
    testErrorRetVA(ok==0,-1234,"Cannot compute prediction for %g Ghz channel",*err,__LINE__,,egl->freqlist[m1]);
  }
}
void parametric_check_freq_P(parametric *egl, double* freqlist, int nfreq, error **err) {
  int m1,m2, ok;

  for(m1=0;m1<egl->nfreq_P;m1++) {
    ok = 0;
    for(m2=0;m2<nfreq;m2++) {
      if (fabs(egl->freqlist_P[m1]-freqlist[m2])<1e-6) {
        ok=1;
        break;
      }
    }
    testErrorRetVA(ok==0,-1234,"Cannot compute prediction for %g Ghz channel",*err,__LINE__,,egl->freqlist[m1]);
  }
}

void parametric_template_payload_init(parametric *egl, double *template, int template_size, double *hfi_freqlist, int nfreqs_hfi, error **err) {
  int m1,m2;
  template_payload *payload;

  parametric_check_freq(egl, hfi_freqlist, nfreqs_hfi,err);
  forwardError(*err,__LINE__,);

  egl->payload = malloc_err(sizeof(template_payload),err);
  forwardError(*err,__LINE__,);

  payload = egl->payload;

  payload->ind_freq = malloc_err(sizeof(int)*egl->nfreq,err);
  forwardError(*err,__LINE__,);

  for (m1=0;m1<egl->nfreq;m1++) {
    for(m2=0;m2<nfreqs_hfi;m2++) {
      if (fabs(egl->freqlist[m1]-hfi_freqlist[m2])<1e-6) {
        payload->ind_freq[m1]=m2;
      }
    }
  }

  payload->template = malloc_err(sizeof(double)*template_size,err);
  forwardError(*err,__LINE__,);

  memcpy(payload->template,template,sizeof(double)*template_size);
}




void parametric_simple_payload_free(void **pp) {
  void *p;

  p = *pp;
  if (p!=NULL) {
    free(p);
  }
  *pp=NULL;
}


void parametric_template_payload_free(void **pp) {
  template_payload *p;
  p = *pp;
  if (p!=NULL) {
    free(p->template);
    free(p->ind_freq);
    free(p);
  }
  *pp=NULL;
}

