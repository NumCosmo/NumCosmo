#include "clik_parametric.h"
#include "clik_parametric_addon.h"


void gal545_compute(parametric *egl, double *Rq, error **err);

parametric *gal545_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  double fac;
  int ell, m1, m2;
  double dell;
  template_payload *payload;
  pfchar name;
  
  egl = parametric_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);

  
  parametric_set_default(egl,"gal545_h",2.3e-11,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gal545_k",5.05,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gal545_t",56,err);
  forwardError(*err,__LINE__,NULL);
  parametric_set_default(egl,"gal545_index",-2.63,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal545_l_pivot",200,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gal545_check_defpo",0,err);
  forwardError(*err,__LINE__,NULL);

  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"gal545_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"gal545_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      parametric_set_default(egl,name,0,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  
  

  // Declare payload, allocate it and fill it  

  egl->eg_compute = &gal545_compute;
  egl->eg_free = &parametric_simple_payload_free;
  
  egl->payload = malloc_err(sizeof(double)*sizeof(double)*egl->nfreq*egl->nfreq,err);

  return egl;
}

void gal545_compute(parametric *egl, double *Rq, error **err) {

  int ell, m1, m2, m2ll, nfreq;
  int *ind_freq;
  int ind1, ind2;
  double *template;
  template_payload *payload;
  double dell;
  double gpe_dust_norm;
  double h,k,t,n,l_pivot;
  double *AA;
  double nrm,v;
  pfchar name;
  int defpo;

  nfreq = egl->nfreq;
  
  l_pivot = parametric_get_value(egl,"gal545_l_pivot",err);
  forwardError(*err,__LINE__,); 

  h = parametric_get_value(egl,"gal545_h",err);
  forwardError(*err,__LINE__,);
  
  k = parametric_get_value(egl,"gal545_k",err);
  forwardError(*err,__LINE__,);

  t = parametric_get_value(egl,"gal545_t",err);
  forwardError(*err,__LINE__,);

  n = parametric_get_value(egl,"gal545_index",err);
  forwardError(*err,__LINE__,);

  AA = (double*) egl->payload;
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"gal545_A_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"gal545_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      AA[m1*nfreq+m2] = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      AA[m2*nfreq+m1] = AA[m1*nfreq+m2];      
    }
  }
  
  defpo = parametric_get_value(egl,"gal545_check_defpo",err);
  forwardError(*err,__LINE__,);

  if (defpo!=0) {
    for(m1=0;m1<egl->nfreq;m1++) {
      if (AA[m1*nfreq+m1]==0) {
        continue;
      }
      for(m2=m1+1;m2<egl->nfreq;m2++) {
      if (AA[m2*nfreq+m2]==0) {
        continue;
      }
      testErrorRetVA(AA[m1*nfreq+m2]>sqrt(AA[m1*nfreq+m1] * AA[m2*nfreq+m2]),-130,"invalid dust amplitude (%d %d)",*err,__LINE__,,m1,m2)
      }
    }
  }

  nrm = (h * pow(l_pivot,k) * exp(-l_pivot/t) + 1) * 200*201/2./M_PI;

  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = (h * pow(ell,k) * exp(-ell/t) + 1) * pow(ell/l_pivot,n);
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
      Rq[IDX_R(egl,ell,m1,m2)] = AA[m1*nfreq+m2] * v/nrm;
      Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;
}

CREATE_PARAMETRIC_FILE_INIT(gal545,gal545_init);

typedef struct {
  int n;
  int *m1,*m2;
  double *nrm;
  char **nrm_names;
  int off1,off2;
} pw_XX_payload;



#define EE_KIND 1
#define BB_KIND 2
#define TE_KIND 3
#define TB_KIND 4
#define EB_KIND 5


pw_XX_payload*  init_pw_XX_payload(int kind,int nT,int nP,int* has_TEB,error **err) {
  pw_XX_payload* payload;
  int mx,i,j;
  int lim1,lim2,off1,off2,mul1;

  payload  = malloc_err(sizeof(pw_XX_payload),err);  
  forwardError(*err,__LINE__,NULL);
  
  if (kind==EE_KIND) {
    lim1 = nP;
    lim2 = nP;
    off1 = nT;
    off2 = nT;
    mul1 = 1;
    payload->n = (nP*(nP+1))/2;
  } else if (kind==BB_KIND) {
    lim1 = nP;
    lim2 = nP;
    off1 = nT+nP*has_TEB[1];
    off2 = nT+nP*has_TEB[1];
    mul1 = 1;
    payload->n = (nP*(nP+1))/2;
  } else if (kind==TE_KIND) {
    lim1 = nT;
    lim2 = nP;
    off1 = 0;
    off2 = nT;
    mul1 = 0;
    payload->n = nT * nP;
   } else if (kind==TB_KIND) {
    lim1 = nT;
    lim2 = nP;
    off1 = 0;
    off2 = nT+nP*has_TEB[1];
    mul1 = 0;
    payload->n = nT * nP;
  } else if (kind==EB_KIND) {
    lim1 = nP;
    lim2 = nP;
    off1 = nT;
    off2 = nT+nP*has_TEB[1];
    mul1 = 0;
    payload->n = nP * nP;
  } else {
    testErrorRetVA(1==1,-130,"invalid kind (%d)",*err,__LINE__,NULL,kind)
  }
    
  payload->off1 = off1;
  payload->off2 = off2;

  payload->m1 = malloc_err(sizeof(int)*payload->n,err);
  forwardError(*err,__LINE__,NULL);
  payload->m2 = malloc_err(sizeof(int)*payload->n,err);
  forwardError(*err,__LINE__,NULL);
  
  mx = 0;

  for(i=0;i<lim1;i++) {
    for(j=0 + i*mul1;j<lim2;j++) {
      payload->m1[mx] = i + off1;
      payload->m2[mx] = j + off2;
      mx++;
    }
  }

  payload->nrm_names = NULL;

  //mx = nT + nP;
  //
  //if (mx < payload->n) {
  //  mx = payload->n;
  //} 

  payload->nrm = malloc_err(sizeof(double)*(nT+nP*2)*(nT+nP*2),err);
  forwardError(*err,__LINE__,NULL);
  
  return payload;
}


void pw_XX_free(void **PP) {
  pw_XX_payload *P;
  P = *PP;
  free(P->nrm);
  free(P->m1);
  free(P->m2);
  if (P->nrm_names!=NULL) {
    free(P->nrm_names[0]);
    free(P->nrm_names);  
  }
  free(P);
  *PP = NULL;
}


void powerlaw_free_emissivity_TE_compute(parametric* egl, double *Rq,  error **err);

parametric *powerlaw_free_emissivity_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int mx,i,j;
  pw_XX_payload *payload;
  int tl, mm, lbf,l,delta,m1,m2;
  pfchar name;
  
  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &powerlaw_free_emissivity_TE_compute;
  egl->eg_free = &pw_XX_free;

  egl->payload  = init_pw_XX_payload(TE_KIND,egl->nfreq_T,egl->nfreq_P,has_TEB,err);
  forwardError(*err,__LINE__,NULL);
  payload = egl->payload;

  parametric_set_default(egl,"pwfe_TE_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"pwfe_TE_l2_norm",1,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"pwfe_TE_index",0,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<payload->n;i++) {
    m1 = payload->m1[i];
    m2 = payload->m2[i];
    if (m1 > m2 - ndet_T) {
      continue;
    }
    if (m1==m2 - ndet_T) {
      sprintf(name,"pwfe_TE_A_%d",(int)egl->freqlist[m1]);
    } else {
      sprintf(name,"pwfe_TE_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
    }
    parametric_set_default(egl,name,0,err);
    forwardError(*err,__LINE__,NULL);
  }

  return egl;
}

void powerlaw_free_emissivity_TE_compute(parametric* egl, double *Rq,  error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,i;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;
  pw_XX_payload *payload;
  double nrmit;
  int l2norm;

  l_pivot = parametric_get_value(egl,"pwfe_TE_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  l2norm = parametric_get_value(egl,"pwfe_TE_l2_norm",err);
  forwardError(*err,__LINE__,);

  nrmit = 1;
  if (l2norm==1) {
    nrmit = l_pivot*(l_pivot+1.)/2./M_PI;
  }
  
  index = parametric_get_value(egl,"pwfe_TE_index",err);
  forwardError(*err,__LINE__,);

  payload = egl->payload;

  A = payload->nrm;
  nfreq = egl->nfreq;

  for(i=0;i<payload->n;i++) {
    m1 = payload->m1[i];
    m2 = payload->m2[i];
    if (m1 > m2 - egl->ndet_T) {
      int mtmp;
      mtmp = m2-egl->ndet_T;
      m2 = m1+egl->ndet_T;
      m1 = mtmp;
    }
    if (m1==m2-egl->ndet_T) {
      sprintf(name,"pwfe_TE_A_%d",(int)egl->freqlist[m1]);
    } else {
      sprintf(name,"pwfe_TE_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
    }
    v = 1;
    v = parametric_get_value(egl,name,err);
    A[i] = v/nrmit;
  }

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/l_pivot,index);
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(i=0;i<payload->n;i++) {
      m1 = payload->m1[i];
      m2 = payload->m2[i];
      lA = A[i];
      Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
      Rq[IDX_R(egl,ell,m2,m1)] = lA*v;
    }  
  }

  return;
}

void powerlaw_free_emissivity_EE_compute(parametric* egl, double *Rq,  error **err);

parametric *powerlaw_free_emissivity_EE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int mx,i,j;
  pw_XX_payload *payload;
  int tl, mm, lbf,l,delta,m1,m2;
  pfchar name;
  
  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &powerlaw_free_emissivity_EE_compute;
  egl->eg_free = &pw_XX_free;

  egl->payload  = init_pw_XX_payload(EE_KIND,egl->nfreq_T,egl->nfreq_P,has_TEB,err);
  forwardError(*err,__LINE__,NULL);
  payload = egl->payload;

  parametric_set_default(egl,"pwfe_EE_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"pwfe_EE_l2_norm",1,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"pwfe_EE_index",0,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"pwfe_EE_check_defpo",0,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<payload->n;i++) {
    m1 = payload->m1[i];
    m2 = payload->m2[i];
    
    if (m1==m2) {
      sprintf(name,"pwfe_EE_A_%d",(int)egl->freqlist[m1]);
    } else {
      sprintf(name,"pwfe_EE_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
    }
    parametric_set_default(egl,name,0,err);
    forwardError(*err,__LINE__,NULL);
  }

  return egl;
}

void powerlaw_free_emissivity_EE_compute(parametric* egl, double *Rq,  error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;
  pw_XX_payload *payload;
  double nrmit;
  int l2norm,i;
  int defpo;

  l_pivot = parametric_get_value(egl,"pwfe_EE_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  l2norm = parametric_get_value(egl,"pwfe_EE_l2_norm",err);
  forwardError(*err,__LINE__,);

  nrmit = 1;
  if (l2norm==1) {
    nrmit = l_pivot*(l_pivot+1.)/2./M_PI;
  }
  
  index = parametric_get_value(egl,"pwfe_EE_index",err);
  forwardError(*err,__LINE__,);

  payload = egl->payload;

  A = payload->nrm;
  nfreq = egl->nfreq;

  for(i=0;i<payload->n;i++) {
    m1 = payload->m1[i];
    m2 = payload->m2[i];
    if (m1==m2) {
      sprintf(name,"pwfe_EE_A_%d",(int)egl->freqlist[m1]);
    } else {
      sprintf(name,"pwfe_EE_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
    }
    v = 1;
    v = parametric_get_value(egl,name,err);
    A[m1*nfreq+m2] = v/nrmit;
    A[m2*nfreq+m1] = A[m1*nfreq+m2];
  }

  defpo = parametric_get_value(egl,"pwfe_EE_check_defpo",err);
  forwardError(*err,__LINE__,);

  if (defpo!=0) {
    for(m1=0;m1<egl->nfreq;m1++) {
      if (A[m1*nfreq+m1]==0) {
        continue;
      }
      for(m2=m1+1;m2<egl->nfreq;m2++) {
      if (A[m2*nfreq+m2]==0) {
        continue;
      }
      testErrorRetVA(A[m1*nfreq+m2]>sqrt(A[m1*nfreq+m1] * A[m2*nfreq+m2]),-130,"invalid dust amplitude (%d %d)",*err,__LINE__,,m1,m2)
      }
    }
  }
  
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/l_pivot,index);
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(i=0;i<payload->n;i++) {
      m1 = payload->m1[i];
      m2 = payload->m2[i];
      lA = A[m1*nfreq+m2];
      Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
      Rq[IDX_R(egl,ell,m2,m1)] = lA*v;
    }  
  }

  return;
}

CREATE_PARAMETRIC_POL_FILE_INIT(powerlaw_free_emissivity_EE,powerlaw_free_emissivity_EE_init);
CREATE_PARAMETRIC_POL_FILE_INIT(powerlaw_free_emissivity_TE,powerlaw_free_emissivity_TE_init);


void pointsource_compute(parametric* egl, double *Rq,  error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv;
  double l_pivot,index,v,lA,nrm;
  double *A;
  pfchar name;
  int stop;

  l_pivot = parametric_get_value(egl,"ps_l_pivot",err);
  forwardError(*err,__LINE__,);

  nrm = l_pivot*(l_pivot+1)/2/M_PI;
  
  A = egl->payload;
  nfreq = egl->nfreq;
  for(m1=0;m1<nfreq;m1++) {
    for(m2=m1;m2<nfreq;m2++) {
      sprintf(name,"ps_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      v = 1;
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      A[m1*nfreq+m2] = v/nrm;
      A[m2*nfreq+m1] = v/nrm;
    }
  }

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = 1;
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


parametric *pointsource_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int m1,m2;
  pfchar name;

  egl = parametric_init( ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  egl->eg_compute = &pointsource_compute;
  egl->eg_free = &parametric_simple_payload_free;
  egl->payload = malloc_err(sizeof(double)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"ps_l_pivot",3000,err);
  forwardError(*err,__LINE__,NULL);
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      sprintf(name,"ps_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      parametric_set_default(egl,name,1,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  

  return egl;
}

CREATE_PARAMETRIC_FILE_INIT(pointsource,pointsource_init);

void fill_offset_freq(int idreq,double *dreq, parametric *egl,int *mv,int def, error **err) {
  int m,i;

  for(m=0;m<egl->nfreq;m++) {
    double f;
    f = egl->freqlist[m];
    mv[m]=def;
    for(i=0;i<idreq;i++) {
      //_DEBUGHERE_("%g %d",f,dreq[i]);
      if (fabs(f-dreq[i])<1e-6) {
        mv[m]=i;
        break;
      }  
    }
    testErrorRetVA(mv[m]==-1,-431432,"Don't know how to compute component for freq %g",*err,__LINE__,,f);
  }    
}

#define gib_lmax_corr_template  10000
#define PRM_NU0 143.
#define lmin_sz_template 2
#define lmax_sz_template  10000                               
#define lmin_corr_template  2
#define lmax_corr_template  9999

void gcib_compute(parametric *egl, double *rq, error **err);

parametric *gcib_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int i,m,*mv,m1,m2;
  double dreq[4];
  pfchar name;
  double *conv,*A;


  egl = parametric_init(ndet, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);

  egl->eg_compute = &gcib_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)* (10001*4*4) + sizeof(double)*(4+egl->nfreq*egl->nfreq)+sizeof(int)*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);


  memcpy(egl->payload,template,sizeof(double)* (10001*4*4));

  mv = egl->payload + sizeof(double)* (10001*4*4) + sizeof(double)*(4+egl->nfreq*egl->nfreq);
  conv = egl->payload + sizeof(double)* (10001*4*4);

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;
  fill_offset_freq(4,dreq, egl,mv,-1,err);
  forwardError(*err,__LINE__,NULL);

  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"A_gib_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"A_gib_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      parametric_set_default(egl,name,70,err);
      forwardError(*err,__LINE__,NULL);
    }
  }  

  parametric_set_default(egl,"gib_index",-1.3,err); 
  forwardError(*err,__LINE__,NULL);

  // conversion factor from table 6
  //4096.68168783,  2690.05218701,  2067.43988919,  3478.86588972
  parametric_set_default(egl,"gib_muK_MJ-2sr_100",4096.68168783/1e6,err); 
  forwardError(*err,__LINE__,NULL);
  
  conv[0] = parametric_get_value(egl,"gib_muK_MJ-2sr_100",err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gib_muK_MJ-2sr_143",2690.05218701/1e6,err); 
  forwardError(*err,__LINE__,NULL);

  conv[1] = parametric_get_value(egl,"gib_muK_MJ-2sr_143",err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gib_muK_MJ-2sr_217",2067.43988919/1e6,err); 
  forwardError(*err,__LINE__,NULL);

  conv[2] = parametric_get_value(egl,"gib_muK_MJ-2sr_217",err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gib_muK_MJ-2sr_353",3478.86588972/1e6,err); 
  forwardError(*err,__LINE__,NULL);

  conv[3] = parametric_get_value(egl,"gib_muK_MJ-2sr_353",err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gib_rigid",217,err); 
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gib_l_pivot",3000,err); 
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void gcib_compute(parametric *egl, double *Rq, error **err) {
  double *template;
  double l_pivot,index,v;
  int m1,m2,ell;
  double nrm;
  int *mv;
  int rigid,irigid;
  double *conv,*A;
  pfchar name;
  

  mv = egl->payload + sizeof(double)* (10001*4*4) + sizeof(double)*(4+egl->nfreq*egl->nfreq);
  conv = egl->payload + sizeof(double)* (10001*4*4);
  A = egl->payload + sizeof(double)* (10001*4*4)+sizeof(double)*4;

  template = egl->payload;
  l_pivot = parametric_get_value(egl,"gib_l_pivot",err);
  forwardError(*err,__LINE__,);;

  rigid = parametric_get_value(egl,"gib_rigid",err);
  forwardError(*err,__LINE__,);
  //_DEBUGHERE_("rigid %d",rigid);

  index = parametric_get_value(egl,"gib_index",err);
  forwardError(*err,__LINE__,);

  if (rigid==0) {
    sprintf(name,"A_gib_%d",217);    
    irigid = 2;
  } else {
    double dreq[4];

    dreq[0] = 100;
    dreq[1] = 143;
    dreq[2] = 217;
    dreq[3] = 353;

    sprintf(name,"A_gib_%d",(int) rigid);  
    irigid=-1;
    for(m1=0;m1<4;m1++) {
      if (fabs(rigid-dreq[m1])<1e-6) {
        irigid=m1;
      }
    }
    testErrorRet(irigid==-1,-55214,"AAAAAAA",*err,__LINE__,);
    rigid = 1;
    
  }


  nrm = parametric_get_value(egl,name,err);
  forwardError(*err,__LINE__,);
  //_DEBUGHERE_("%s %g",name,nrm);
    
  
  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      if (m1==m2) {
        sprintf(name,"A_gib_%d",(int)egl->freqlist[m1]);  
      } else {
        sprintf(name,"A_gib_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
      }
      v = parametric_get_value(egl,name,err);
      forwardError(*err,__LINE__,);
      A[m1*egl->nfreq+m2] = (v/template[((int) l_pivot)*16+mv[m1]*4+mv[m2]]*(1-rigid) +  nrm/template[((int)l_pivot)*16+irigid*4+irigid]*rigid*conv[mv[m1]]*conv[mv[m2]]/conv[irigid]/conv[irigid]) /l_pivot/(l_pivot+1)*2*M_PI ;
      A[m2*egl->nfreq+m1] = A[m1*egl->nfreq+m2];
    }
  }

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow((double) ell/l_pivot,(double) index-(-1.3));
    //_DEBUGHERE_("%d %g",ell,v);
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        //_DEBUGHERE_("%d %d %d %g",ell,mv[m1],mv[m2],template[ell*16+mv[m1]*4+mv[m2]]);
        Rq[IDX_R(egl,ell,m1,m2)] = v*template[ell*16+mv[m1]*4+mv[m2]] * A[m2*egl->nfreq+m1];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
}

void gibXsz_compute(parametric *egl, double *Rq, error **err) {
  double a_cib,xi_sz_cib,a_sz;
  double *conv;
  double *corr_template;
  double l_pivot,index;
  int m1,m2,ell;
  double nrm;
  int *mv;
  double *fnu ,*A;
  double v;

  mv = egl->payload + sizeof(double)* (gib_lmax_corr_template+1) + sizeof(double)*(4+egl->nfreq+egl->nfreq*egl->nfreq);
  conv = egl->payload + sizeof(double)* (gib_lmax_corr_template+1);
  fnu = egl->payload + sizeof(double)* (gib_lmax_corr_template+1) + sizeof(double)*(4);
  A = egl->payload + sizeof(double)* (gib_lmax_corr_template+1) + sizeof(double)*(4+egl->nfreq);
  corr_template = egl->payload;
  
  a_cib = parametric_get_value(egl,"A_cib_217",err);
  forwardError(*err,__LINE__,);
  a_sz = parametric_get_value(egl,"A_sz",err);
  forwardError(*err,__LINE__,);
  xi_sz_cib = parametric_get_value(egl,"xi_sz_cib",err);
  forwardError(*err,__LINE__,);

  for(m1=0;m1<egl->nfreq;m1++) {
    for(m2=m1;m2<egl->nfreq;m2++) {
      A[m1*egl->nfreq+m2] = - xi_sz_cib * sqrt(a_sz) * ( sqrt(fnu[m1]*a_cib*conv[mv[m2]]) + sqrt(fnu[m2]*a_cib*conv[mv[m1]]));
    }
  }

  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = corr_template[ell];
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        Rq[IDX_R(egl,ell,m1,m2)] = A[m1*egl->nfreq+m2] * corr_template[ell] * 2.0*M_PI/(ell*(ell+1.0));
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
        
      }  
    }
  } 
}


parametric *gibXsz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  double fac;
  int l,i;
  double *fnu;
  int m1;
  int *mv;
  double *corr_template;
  double dreq[4];
  double *conv;
  int remove_100;
  double szcolor[4];
  egl = parametric_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);

  egl->payload = malloc_err(sizeof(double)*(gib_lmax_corr_template+1 + egl->nfreq + 4 +egl->nfreq*egl->nfreq)+sizeof(int)*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);
  
  memcpy(egl->payload+sizeof(double)*2,template,sizeof(double)* (gib_lmax_corr_template-1));

  ((double*)egl->payload)[0] = 0;
  ((double*)egl->payload)[1] = 0;
  
  mv = egl->payload + sizeof(double)* (gib_lmax_corr_template+1) + sizeof(double)*(4+egl->nfreq+egl->nfreq*egl->nfreq);
  conv = egl->payload + sizeof(double)* (gib_lmax_corr_template+1);
  fnu = egl->payload + sizeof(double)* (gib_lmax_corr_template+1) + sizeof(double)*(4);
  
  parametric_set_default(egl,"sz_color_143_to_143",0.975,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"sz_color_100_to_143",0.981,err);
  forwardError(*err,__LINE__,NULL);

  szcolor[0] = parametric_get_value(egl,"sz_color_100_to_143",err);
  forwardError(*err,__LINE__,NULL);
  szcolor[1] = parametric_get_value(egl,"sz_color_143_to_143",err);
  forwardError(*err,__LINE__,NULL);

  szcolor[2]=1;
  szcolor[3]=1;

  parametric_set_default(egl,"gibXsz_100_to_217",0.022,err); 
  forwardError(*err,__LINE__,NULL);
  
  conv[0] = parametric_get_value(egl,"gibXsz_100_to_217",err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"gibXsz_143_to_217",0.094,err); 
  forwardError(*err,__LINE__,NULL);
  
  conv[1] = parametric_get_value(egl,"gibXsz_143_to_217",err);
  forwardError(*err,__LINE__,NULL);

  conv[2] = 1;

  parametric_set_default(egl,"gibXsz_353_to_217",46.8,err); 
  forwardError(*err,__LINE__,NULL);
  
  conv[3] = parametric_get_value(egl,"gibXsz_353_to_217",err);
  forwardError(*err,__LINE__,NULL);
  
  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;
  fill_offset_freq(4,dreq, egl,mv,-1,err);
  forwardError(*err,__LINE__,NULL);
  
  
  // Compute SZ spectrum
  for (m1=0;m1<egl->nfreq;m1++) {
    fnu[m1] = sz_spectrum((double)egl->freqlist[m1],PRM_NU0)*szcolor[mv[m1]];
  }
  
  parametric_set_default(egl,"no_szxcib_100",1,err); 
  forwardError(*err,__LINE__,NULL);
  
  remove_100 = parametric_get_value(egl,"no_szxcib_100",err);
  forwardError(*err,__LINE__,NULL);
  
  if (remove_100!=0) {
    conv[0] = 0;
    // this only removes the 100 auto
  }

  egl->eg_compute = &gibXsz_compute;
  egl->eg_free = &parametric_simple_payload_free;
  
  //parametric_declare_mandatory(egl,"A_cib_217",err);
  parametric_set_default(egl,"A_cib_217",70,err); 
  forwardError(*err,__LINE__,NULL);
  //parametric_declare_mandatory(egl,"A_sz",err);
  parametric_set_default(egl,"A_sz",4.0,err);
  forwardError(*err,__LINE__,NULL);
  //parametric_declare_mandatory(egl,"xi_sz_cib",err);
  parametric_set_default(egl,"xi_sz_cib",0.0,err);
  forwardError(*err,__LINE__,NULL);
  
  return egl;
}

CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(gibXsz,gibXsz_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(gcib,gcib_init);

// SZ alone
void sz_compute(parametric* egl, double *Rq, error **err);

parametric *sz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  double fac;
  int l;
  double *fnu;
  int m1;
  double szcolor[4],dreq[4];
  int mv[100];
  
  // make sure l(l+1)/(2pi)*C_l template is normalized to 1 at l=3000 
  fac = template[3000-lmin_sz_template];
  for (l=lmin_sz_template;l<=lmax_sz_template;l++) {
    template[l-lmin_sz_template] /= fac;
  }
  egl = parametric_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);


  // Declare payload, allocate it and fill it

  egl->payload = malloc_err(sizeof(double)*(lmax_sz_template-lmin_sz_template+1 + egl->nfreq),err);
  forwardError(*err,__LINE__,NULL);
  memcpy(egl->payload,template,sizeof(double)*(lmax_sz_template-lmin_sz_template+1));
  
  fnu = egl->payload + (lmax_sz_template-lmin_sz_template+1)*sizeof(double);
  parametric_set_default(egl,"sz_color_143_to_143",1,err);
  forwardError(*err,__LINE__,NULL);
  
  parametric_set_default(egl,"sz_color_100_to_143",1,err);
  forwardError(*err,__LINE__,NULL);

  szcolor[0] = parametric_get_value(egl,"sz_color_100_to_143",err);
  forwardError(*err,__LINE__,NULL);
  szcolor[1] = parametric_get_value(egl,"sz_color_143_to_143",err);
  forwardError(*err,__LINE__,NULL);

  szcolor[2]=1;
  szcolor[3]=1;

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;
  fill_offset_freq(4,dreq, egl,mv,-1,err);
  forwardError(*err,__LINE__,NULL);
  // Compute SZ spectrum
  for (m1=0;m1<egl->nfreq;m1++) {
    fnu[m1] = sz_spectrum((double)egl->freqlist[m1],PRM_NU0)*szcolor[mv[m1]];
  }

  egl->eg_compute = &sz_compute;
  egl->eg_free = &parametric_simple_payload_free;
  
  parametric_set_default(egl,"A_sz",4.0,err);
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"A_sz",&parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}


void sz_compute(parametric* egl, double *Rq, error **err) {

  int ell,m1,m2,mell,nfreq;
  double *cl, *fnu, *A;
  double sz_norm, dell;
  
  nfreq = egl->nfreq;
  A = (double*) egl->payload;
  cl = A;
  fnu = &(A[lmax_sz_template-lmin_sz_template+1]);
  
  sz_norm = parametric_get_value(egl,"A_sz",err);
  forwardError(*err,__LINE__,);

  
  // Create covariance matrix
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    dell = (double)ell;
    mell = (ell-lmin_sz_template);
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
  Rq[IDX_R(egl,ell,m1,m2)] = sz_norm * 2.0*M_PI/(dell*(dell+1.0)) * cl[mell] * fnu[m1] * fnu[m2];
  Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;

}

// kSZ alone

void ksz_compute(parametric* egl, double *Rq, error **err);

parametric *ksz_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int lmin_ksz_template= 0; //CHECK
  int lmax_ksz_template= 5000; // CHECK
  double fac;
  int l;

  // make sure l(l+1)/(2pi)*C_l template is normalized to 1 at l=3000 
  fac = template[3000-lmin_ksz_template];
  for (l=lmin_ksz_template;l<=lmax_ksz_template;l++) {
    template[l-lmin_ksz_template] /= fac;
  }
  egl = parametric_init(ndet,detlist,ndef,defkey,defvalue,nvar,varkey,lmin,lmax,err);
  forwardError(*err,__LINE__,NULL);

  // Declare payload, allocate it and fill it

  egl->payload = malloc_err(sizeof(double)*(lmax_ksz_template-lmin_ksz_template+1),err);
  forwardError(*err,__LINE__,NULL);
  memcpy(egl->payload,template,sizeof(double)*(lmax_ksz_template-lmin_ksz_template+1));
  
  egl->eg_compute = &ksz_compute;
  egl->eg_free = &parametric_simple_payload_free;
  
  parametric_set_default(egl,"ksz_norm",0.0,err); // PICK YOUR FAVORITE
  forwardError(*err,__LINE__,NULL);
  parametric_add_derivative_function(egl,"ksz_norm",&parametric_norm_derivative,err);
  forwardError(*err,__LINE__,NULL);

  return egl;
}


void ksz_compute(parametric* egl, double *Rq, error **err) {

  int ell,m1,m2,mell,nfreq;
  int lmin_ksz_template = 0; // CHECK
  int lmax_ksz_template = 5000; // CHECK
  double *cl,*A;
  double ksz_norm, dell;
  
  nfreq = egl->nfreq;
  A = (double*) egl->payload;
  cl = A;
  
  ksz_norm = parametric_get_value(egl,"ksz_norm",err);
  forwardError(*err,__LINE__,);

  // kSZ spectrum is like CMB

  // Create covariance matrix
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    dell = (double)ell;
    mell = (ell-lmin_ksz_template);
    for (m1=0;m1<nfreq;m1++) {
      for (m2=m1;m2<nfreq;m2++) {
  Rq[IDX_R(egl,ell,m1,m2)] = ksz_norm * 2.0*M_PI/(dell*(dell+1.0)) * cl[mell];
  Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;

}
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(sz,sz_init);
CREATE_PARAMETRIC_TEMPLATE_FILE_INIT(ksz,ksz_init);


void powerlaw_free_emissivity_XX_compute(parametric* egl, double *Rq,  error **err);

parametric *powerlaw_free_emissivity_XX_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err) {
  parametric *egl;
  int mx,i,j,kind;
  pw_XX_payload *payload;
  int tl, mm, lbf,l,delta,m1,m2;
  pfchar name;
  
  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  egl->eg_compute = &powerlaw_free_emissivity_XX_compute;
  egl->eg_free = &pw_XX_free;

  kind = parametric_get_value(egl,"pwfe_XX_kind",err);
  forwardError(*err,__LINE__,NULL);
  
  egl->payload  = init_pw_XX_payload(kind,egl->nfreq_T,egl->nfreq_P,has_TEB,err);
  forwardError(*err,__LINE__,NULL);
  
  payload = egl->payload;

  parametric_set_default(egl,"pwfe_XX_l_pivot",500,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"pwfe_XX_l2_norm",1,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"pwfe_XX_index",0,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<payload->n;i++) {
    m1 = payload->m1[i];
    m2 = payload->m2[i];
    if (m1-payload->off1 > m2-payload->off2) {
      continue;
    }
    if (m1-payload->off1==m2-payload->off2) {
      sprintf(name,"pwfe_XX_A_%d",(int)egl->freqlist[m1]);
    } else {
      sprintf(name,"pwfe_XX_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
    }
    parametric_set_default(egl,name,0,err);
    forwardError(*err,__LINE__,NULL);
  }

  return egl;
}

void powerlaw_free_emissivity_XX_compute(parametric* egl, double *Rq,  error **err) {
  int ell,m1,m2,mell,nfreq,iv,mv,i;
  double l_pivot,index,v,lA;
  double *A;
  pfchar name;
  int stop;
  pw_XX_payload *payload;
  double nrmit;
  int l2norm;

  l_pivot = parametric_get_value(egl,"pwfe_XX_l_pivot",err);
  forwardError(*err,__LINE__,);
  egl->l_pivot = l_pivot;

  l2norm = parametric_get_value(egl,"pwfe_XX_l2_norm",err);
  forwardError(*err,__LINE__,);

  nrmit = 1;
  if (l2norm==1) {
    nrmit = l_pivot*(l_pivot+1.)/2./M_PI;
  }
  
  index = parametric_get_value(egl,"pwfe_XX_index",err);
  forwardError(*err,__LINE__,);

  payload = egl->payload;

  A = payload->nrm;
  nfreq = egl->nfreq;

  for(i=0;i<payload->n;i++) {
    m1 = payload->m1[i];
    m2 = payload->m2[i];
    if (m1-payload->off1 > m2-payload->off2) {
      int mtmp;
      mtmp = m2-payload->off2+payload->off1;
      m2 = m1-payload->off1+payload->off2;
      m1 = mtmp;
    }
    if (m1-payload->off1==m2-payload->off2) {
      sprintf(name,"pwfe_XX_A_%d",(int)egl->freqlist[m1]);
    } else {
      sprintf(name,"pwfe_XX_A_%d_%d",(int)egl->freqlist[m1],(int)egl->freqlist[m2]);
    }
    v = 1;
    v = parametric_get_value(egl,name,err);
    A[i] = v/nrmit;
  }

  #pragma omp parallel for private(ell,v, mell,i,m1,m2,lA)
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    v = pow(ell/l_pivot,index);
    mell = (ell-egl->lmin)*nfreq*nfreq;
    for(i=0;i<payload->n;i++) {
      m1 = payload->m1[i];
      m2 = payload->m2[i];
      lA = A[i];
      Rq[IDX_R(egl,ell,m1,m2)] = lA*v;
      Rq[IDX_R(egl,ell,m2,m1)] = lA*v;
    }  
  }

  return;
}


CREATE_PARAMETRIC_POL_FILE_INIT(powerlaw_free_emissivity_XX,powerlaw_free_emissivity_XX_init);
