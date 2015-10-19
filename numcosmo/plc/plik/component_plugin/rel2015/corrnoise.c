#include "clik_parametric.h"
#include "clik_parametric_addon.h"

void cnoise_compute(parametric *egl, double *rq, error **err);

void fill_offset_freq_TP(int idreq,double *dreq, int nfreq, double* freqlist,int *mv,int off, error **err) {
  int m,i;
  
  for(m=0;m<nfreq;m++) {
    double f;
    f = freqlist[m];
    mv[m]=-1;
    for(i=0;i<idreq;i++) {
      if (fabs(f-dreq[i])<1e-6) {
        mv[m] = i + off;
        break;
      }  
    }
    testErrorRetVA(mv[m]==-1,-431432,"Don't know how to compute component for freq %g",*err,__LINE__,,f);
  }    
}


parametric *cnoise_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int i,m,*mv,m1,m2;
  double dreq[4];
  pfchar name;
  double *conv,*A;
  char teb[3];
  int f1,f2;

  // init
  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  
  egl->eg_compute = &cnoise_compute;
  egl->eg_free = &parametric_simple_payload_free;

  egl->payload = malloc_err(sizeof(double)* (3001*12*12 + 12*12) + sizeof(int)*egl->nfreq*egl->nfreq,err);
  forwardError(*err,__LINE__,NULL);


  memcpy(egl->payload,template,sizeof(double)* (3001*12*12));

  mv = egl->payload + sizeof(double)* (3001*12*12 + 12*12);

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;

  fill_offset_freq_TP(4,dreq, egl->nfreq_T*has_TEB[0],egl->freqlist_T,mv,0,err);
  forwardError(*err,__LINE__,NULL);
  fill_offset_freq_TP(4,dreq, egl->nfreq_P*has_TEB[1],egl->freqlist_P,mv + egl->nfreq_T*has_TEB[0],4,err);
  forwardError(*err,__LINE__,NULL);
  fill_offset_freq_TP(4,dreq, egl->nfreq_P*has_TEB[2],egl->freqlist_P,mv + egl->nfreq_P*has_TEB[1]+egl->nfreq_T*has_TEB[0],8,err);
  forwardError(*err,__LINE__,NULL);

  teb[0] = 'T';
  teb[1] = 'E';
  teb[2] = 'B';

  for(m1=0;m1<4;m1++) {
    for(m2=m1;m2<4;m2++) {
      for(f1=0;f1<3;f1++) {
        for(f2=f1;f2<3;f2++) {
          sprintf(name,"A_cnoise_%d_%d_%c%c",(int)dreq[m1],(int)dreq[m2],teb[f1],teb[f2]);
          parametric_set_default(egl,name,0,err);
          forwardError(*err,__LINE__,NULL);      
        }
      }
    }
  }  

  parametric_set_default(egl,"cnoise_l_pivot",2000,err); 
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"cnoise_abs",0,err); 
  forwardError(*err,__LINE__,NULL);

  return egl;
}

void cnoise_compute(parametric *egl, double *Rq, error **err) {
  double *template;
  double l_pivot,index,v;
  int m1,m2,ell;
  double nrm;
  int *mv;
  int rigid,irigid;
  double *conv,*A;
  pfchar name;
  double dreq[4];
  char teb[3];
  int rm1,rm2;
  int f1,f2;
  int abso;

  dreq[0] = 100;
  dreq[1] = 143;
  dreq[2] = 217;
  dreq[3] = 353;

  teb[0] = 'T';
  teb[1] = 'E';
  teb[2] = 'B';

  mv = egl->payload + sizeof(double)* (3001*12*12 + 12*12);
  A = egl->payload + sizeof(double)* (3001*12*12);

  template = egl->payload;
  l_pivot = parametric_get_value(egl,"cnoise_l_pivot",err);
  forwardError(*err,__LINE__,);;
  abso = parametric_get_value(egl,"cnoise_abs",err);
  forwardError(*err,__LINE__,);;

  if (abso!=0) {
    abso = 1;
  }

  for(m1=0;m1<4;m1++) {
    for(m2=m1;m2<4;m2++) {
      for(f1=0;f1<3;f1++) {
        for(f2=f1;f2<3;f2++) {
          sprintf(name,"A_cnoise_%d_%d_%c%c",(int)dreq[m1],(int)dreq[m2],teb[f1],teb[f2]);
          v = parametric_get_value(egl,name,err);
          forwardError(*err,__LINE__,);
          rm1 = m1 + f1*4;
          rm2 = m2 + f2*4;
          A[rm1*12+rm2] = v/(1-abso+abso*template[((int) l_pivot)*12*12+rm1*12+rm2]/l_pivot/(l_pivot+1)*2*M_PI);  
          A[rm2*12+rm1] = A[rm1*12+rm2];
          //_DEBUGHERE_("%s %g %g %d %d %g",name,v,(1-abso+abso*template[((int) l_pivot)*12*12+rm1*12+rm2]/l_pivot/(l_pivot+1)*2*M_PI),rm1,rm2,A[rm1*12+rm2]);
          rm1 = m1 + f2*4;
          rm2 = m2 + f1*4;
          A[rm1*12+rm2] = v/(1-abso+abso*template[((int) l_pivot)*12*12+rm1*12+rm2]/l_pivot/(l_pivot+1)*2*M_PI);  
          A[rm2*12+rm1] = A[rm1*12+rm2];
          //_DEBUGHERE_("%s %g %g %d %d %g",name,v,(1-abso+abso*template[((int) l_pivot)*12*12+rm1*12+rm2]/l_pivot/(l_pivot+1)*2*M_PI),rm1,rm2,A[rm1*12+rm2]);
                
        }
      }
    }
  }  
  
  for(ell=egl->lmin;ell<=egl->lmax;ell++) {
    for(m1=0;m1<egl->nfreq;m1++) {
      for(m2=m1;m2<egl->nfreq;m2++) {
        if(ell==egl->lmin) {
          //_DEBUGHERE_("%d %d %d %g",ell,mv[m1],mv[m2],template[ell*12*12+mv[m1]*12+mv[m2]]);  
        }
        
        Rq[IDX_R(egl,ell,m1,m2)] = template[ell*12*12+mv[m1]*12+mv[m2]] * A[mv[m1]*12+mv[m2]];
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }  
    }
  }
}

CREATE_PARAMETRIC_POL_TEMPLATE_FILE_INIT(cnoise,cnoise_init);
