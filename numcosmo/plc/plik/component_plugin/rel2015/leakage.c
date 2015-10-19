// SZ alone
#include "clik_parametric.h"
#include "clik_parametric_addon.h"


typedef struct {
  double *templ, *beta_coeff, *epsilon_coeff;
  int *get_cli;
  int np;
} bleak_data;

void bleak_free(void **pdlata) {
  bleak_data *dlata;
  dlata = *pdlata;
  free(dlata->templ);
  free(dlata->beta_coeff);
  free(dlata->epsilon_coeff);
  free(dlata->get_cli);
  free(dlata);
  *pdlata = NULL;
}

void bleak_compute(parametric* egl, double *Rq, error **err);

double bleak_polynomial(long ncoeff, double *coeff, double x) {
  /* return coeff[0] + coeff[1]*x + coeff[2]*x^2 + .... */
  double p;
  long d;
  p = coeff[ncoeff-1];
  for (d = ncoeff-2 ; d >= 0 ; d--){
    p *= x ;  
    p += coeff[d];
  }
  return p;
}

double bleak_apply_t2pmatrix(double e,double b, int nk, double *cl) {
     /* apply matrix coupling TT, EE, ... spectra because of beam mismatch at line nk
       cf Eq 11 of  http://hfilfi.planck.fr/index.php/WG/CorrectingPolClFromBeamMismatch */
  int i,j;
  double blatrix;

  //_DEBUGHERE_("%g %g %d %g",e,b,nk,cl[0])
  if (nk==1) {
    blatrix = cl[0] *e*e + cl[3] * 2*e;
  } else if (nk==2) {
    blatrix = cl[0] *b*b + cl[4] * 2*b;
  } else if (nk==3) {
    blatrix = cl[0] *e;
  } else if (nk==4) {
    blatrix = cl[0] *b;
  } else if (nk==5) {
    blatrix = cl[0] *e*b + cl[3] * b + cl[4] * e;
  }

  return blatrix;

  //double zero=(double)0;
  //double one=(double)1;
  ///* identity matrix */
  //for (i=0; i<nk ; i++){
  //  for (j=0; j<nk ; j++){ matrix[i*nk+j] = zero;}
  //  matrix[i*nk+i] = one;
  //}
  ///* spectra are in order TT, EE, BB, TE, TB, EB */
  //matrix[1*nk + 0] = e*e ; matrix[1*nk + 3] = 2*e ;
  //matrix[2*nk + 0] = b*b ;                        matrix[2*nk + 4] = 2*b ;
  //matrix[3*nk + 0] = e ;
  //matrix[4*nk + 0] = b ;
  //matrix[5*nk + 0] = e*b ; matrix[5*nk + 3] = b ; matrix[5*nk + 4] = e;

}

parametric *bleak_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err) {
  parametric *egl;
  int template_lmax;
  int np;
  int mtot;
  int i,j,ip;
  pfchar name;
  pfchar iN, jN;
  bleak_data *dlata;
  int u1,u2,mf1,mf2;
  int ugly[9];
  double dfreq[4];
  int ell,cli;

  // init
  egl = parametric_pol_init(ndet_T, ndet_P, has_TEB, detlist, ndef, defkey, defvalue, nvar, varkey, lmin, lmax, err);
  forwardError(*err,__LINE__,NULL);
  
  template_lmax = 3000;
  dfreq[0] = 100;
  dfreq[1] = 143;
  dfreq[2] = 217;
  dfreq[3] = 353;

  testErrorRet(lmax>=template_lmax,-24342,"lmax too big",*err,__LINE__,NULL);

  parametric_set_default(egl,"bleak_l_pivot",2e4,err);
  forwardError(*err,__LINE__,NULL);

  parametric_set_default(egl,"bleak_np",1,err);
  forwardError(*err,__LINE__,NULL);        
  
  np = parametric_get_value(egl,"bleak_np",err);
  forwardError(*err,__LINE__,NULL);

  mtot = egl->nfreq_T*has_TEB[0] + egl->nfreq_P*has_TEB[1] + egl->nfreq_P*has_TEB[2];

  dlata = malloc_err(sizeof(bleak_data),err);
  forwardError(*err,__LINE__,NULL);

  egl->payload = dlata;

  dlata -> np = np;
  dlata->beta_coeff = malloc_err(sizeof(double)*mtot*mtot*np,err);
  forwardError(*err,__LINE__,NULL);
  dlata->epsilon_coeff = malloc_err(sizeof(double)*mtot*mtot*np,err);
  forwardError(*err,__LINE__,NULL);
  dlata->get_cli = malloc_err(sizeof(int)*mtot*mtot,err);
  forwardError(*err,__LINE__,NULL);
  
  dlata -> templ = malloc_err(sizeof(double)*(lmax+1)*mtot*mtot*6,err);
  forwardError(*err,__LINE__,NULL);

  ugly[0  +3*  0] = 0;
  ugly[1  +3*  1] = 1;
  ugly[2  +3*  2] = 2;
  ugly[1  +3*  0] = 3;
  ugly[0  +3*  1] = 3;
  ugly[2  +3*  0] = 4;
  ugly[0  +3*  2] = 4;
  ugly[2  +3*  1] = 5;
  ugly[1  +3*  2] = 5;

  for(i=0;i<mtot;i++) {
    if (i<egl->nfreq_T*has_TEB[0]) {
      sprintf(iN,"%dT",i);
      u1=0;
    } else if (i<egl->nfreq_T*has_TEB[0] + egl->nfreq_P*has_TEB[1]) {
      sprintf(iN,"%dE",i-egl->nfreq_T*has_TEB[0]);
      u1=1;
    } else {
      sprintf(iN,"%dB",i-egl->nfreq_T*has_TEB[0]-egl->nfreq_P*has_TEB[1]);
      u1=2;
    }
    for (mf1=0;mf1<4;mf1++) {
      if (abs(dfreq[mf1]-egl->freqlist[i])<1e-6) {
        break;
      }
    }
    for (j=i;j<mtot;j++) {
      if (j<egl->nfreq_T*has_TEB[0]) {
        sprintf(jN,"%dT",j);
        u2=0;
      } else if (j<egl->nfreq_T*has_TEB[0] + egl->nfreq_P*has_TEB[1]) {
        sprintf(jN,"%dE",j-egl->nfreq_T*has_TEB[0]);
        u2=1;
      } else {
        sprintf(jN,"%dB",j-egl->nfreq_T*has_TEB[0]-egl->nfreq_P*has_TEB[1]);
        u2=2;
      }
      for (mf2=0;mf2<4;mf2++) {
        if (abs(dfreq[mf2]-egl->freqlist[j])<1e-6) {
          break;
        }
      }
      
      //_DEBUGHERE_("%d %d - >%d %d %d %d",i,j,u1,u2,mf1,mf2);
      dlata->get_cli[i*mtot+j] = ugly[u1*3+u2];
      dlata->get_cli[j*mtot+i] = ugly[u1*3+u2];
      
      for(ell=0;ell<=lmax;ell++) {
        //_DEBUGHERE_("%d %d %d %d %d %d",i,j,ell, ell*mtot*mtot*6+i*mtot*6+j*6 + 0,(lmax+1)*mtot*mtot*6,ell*12*12+( mf1   )*12+( mf2   ));
        dlata->templ[ell*mtot*mtot*6+i*mtot*6+j*6 + 0] = template[ell*12*12+( mf1   )*12+( mf2   )];
        //_DEBUGHERE_("%d %d %d %d %d %d",i,j,ell, ell*mtot*mtot*6+i*mtot*6+j*6 + 1,(lmax+1)*mtot*mtot*6,ell*12*12+( mf1+4   )*12+( mf2+4   ));
        dlata->templ[ell*mtot*mtot*6+i*mtot*6+j*6 + 1] = template[ell*12*12+( mf1+4 )*12+( mf2+4 )];
        //_DEBUGHERE_("%d %d %d %d %d %d",i,j,ell, ell*mtot*mtot*6+i*mtot*6+j*6 + 2,(lmax+1)*mtot*mtot*6,ell*12*12+( mf1+8   )*12+( mf2+8   ));
        dlata->templ[ell*mtot*mtot*6+i*mtot*6+j*6 + 2] = template[ell*12*12+( mf1+8 )*12+( mf2+8 )];
        //_DEBUGHERE_("%d %d %d %d %d %d",i,j,ell, ell*mtot*mtot*6+i*mtot*6+j*6 + 3,(lmax+1)*mtot*mtot*6,ell*12*12+( mf1   )*12+( mf2+4   ));
        dlata->templ[ell*mtot*mtot*6+i*mtot*6+j*6 + 3] = template[ell*12*12+( mf1   )*12+( mf2+4 )];
        //_DEBUGHERE_("%d %d %d %d %d %d",i,j,ell, ell*mtot*mtot*6+i*mtot*6+j*6 + 4,(lmax+1)*mtot*mtot*6,ell*12*12+( mf1   )*12+( mf2+8   ));
        dlata->templ[ell*mtot*mtot*6+i*mtot*6+j*6 + 4] = template[ell*12*12+( mf1   )*12+( mf2+8 )];
        //_DEBUGHERE_("%d %d %d %d %d %d",i,j,ell, ell*mtot*mtot*6+i*mtot*6+j*6 + 5,(lmax+1)*mtot*mtot*6,ell*12*12+( mf1+4   )*12+( mf2+8   ));
        dlata->templ[ell*mtot*mtot*6+i*mtot*6+j*6 + 5] = template[ell*12*12+( mf1+4 )*12+( mf2+8 )];
        //_DEBUGHERE_("%d %d %d",i,j,ell);
        
      }
      cli = ugly[u1*3+u2];
      if (cli==0) {
        continue;
      }
      if (cli==1 || cli ==3 || cli==5) {
        for(ip=0;ip<np;ip++) {
          sprintf(name,"bleak_epsilon_%d_%s_%s",ip,iN,jN);
          parametric_set_default(egl,name,0,err);
          forwardError(*err,__LINE__,NULL);        
        }  
      }
      if (cli==2 || cli ==4 || cli==5) {
        for(ip=0;ip<np;ip++) {
          sprintf(name,"bleak_beta_%d_%s_%s",ip,iN,jN);
          parametric_set_default(egl,name,0,err);
          forwardError(*err,__LINE__,NULL); 
        }
      }
    }
  }

  // Declare payload, allocate it and fill it
  
  egl->eg_compute = &bleak_compute;
  egl->eg_free = &bleak_free;
  
  return egl;
}

void bleak_compute(parametric* egl, double *Rq, error **err) {
  int ell,m1,m2,mell,nfreq;
  int template_lmax;
  int n1,n2;
  int np;
  int mtot;
  int i,j,ip;
  pfchar name;
  pfchar iN, jN;
  double *templ, *beta_coeff, *epsilon_coeff, *leak_matrix;
  int cli;
  double tmp;
  double epsilon,beta;
  bleak_data *dlata;
  double l_pivot;

  dlata = egl->payload;

  np = dlata->np;
  
  mtot = egl->nfreq_T*egl->has_TEB[0] + egl->nfreq_P*egl->has_TEB[1] + egl->nfreq_P*egl->has_TEB[2];

  templ         = dlata->templ;
  beta_coeff    = dlata->beta_coeff;
  epsilon_coeff = dlata->epsilon_coeff;

  l_pivot = parametric_get_value(egl,"bleak_l_pivot",err);
  forwardError(*err,__LINE__,);

  for(i=0;i<mtot;i++) {
    if (i<egl->nfreq_T*egl->has_TEB[0]) {
      sprintf(iN,"%dT",i);
    } else if (i<egl->nfreq_T*egl->has_TEB[0] + egl->nfreq_P*egl->has_TEB[1]) {
      sprintf(iN,"%dE",i-egl->nfreq_T*egl->has_TEB[0]);
    } else {
      sprintf(iN,"%dB",i-egl->nfreq_T*egl->has_TEB[0]-egl->nfreq_P*egl->has_TEB[1]);
    }
    for (j=i;j<mtot;j++) {
      if (j<egl->nfreq_T*egl->has_TEB[0]) {
        sprintf(jN,"%dT",j);
      } else if (j<egl->nfreq_T*egl->has_TEB[0] + egl->nfreq_P*egl->has_TEB[1]) {
        sprintf(jN,"%dE",j-egl->nfreq_T*egl->has_TEB[0]);
      } else {
        sprintf(jN,"%dB",j-egl->nfreq_T*egl->has_TEB[0]-egl->nfreq_P*egl->has_TEB[1]);
      }
      cli = dlata->get_cli[i*mtot+j];
      if (cli==0) {
        continue;
      }
      if (cli==1 || cli ==3 || cli==5) {
        for(ip=0;ip<np;ip++) {
          sprintf(name,"bleak_epsilon_%d_%s_%s",ip,iN,jN);
          epsilon_coeff[i*mtot*np+j*np+ip] = parametric_get_value(egl,name,err);
          forwardError(*err,__LINE__,);        
        }  
      }
      if (cli==2 || cli ==4 || cli==5) {
        for(ip=0;ip<np;ip++) {
          sprintf(name,"bleak_beta_%d_%s_%s",ip,iN,jN);
          beta_coeff[i*mtot*np+j*np+ip] = parametric_get_value(egl,name,err);
          forwardError(*err,__LINE__,); 
        }
      }
    }
  }


  //Do something clever about the computation here !
  
  // Create covariance matrix
  for (ell=egl->lmin;ell<=egl->lmax;ell++) {
    
    // call the routine that fills the leak matrix
    for (m1=0;m1<mtot;m1++) {      
      for (m2=m1;m2<mtot;m2++) {
        //_DEBUGHERE_("epsilon %d %d %g",m1,m2,epsilon_coeff[m1*mtot*np+m2*np])
        cli = dlata->get_cli[m1*mtot+m2];
        if (cli==0) {
          continue;
        }
        if (cli==1 || cli ==3 || cli==5) {
          epsilon = bleak_polynomial(np,&(epsilon_coeff[m1*mtot*np+m2*np]),ell/l_pivot);
        }
        //_DEBUGHERE_("beta %d %d %g",m1,m2,epsilon_coeff[m1*mtot*np+m2*np])
        if (cli==2 || cli ==4 || cli==5) {
          beta    = bleak_polynomial(np,&(beta_coeff[m1*mtot*np+m2*np]),ell/l_pivot);
        }
        tmp = bleak_apply_t2pmatrix(epsilon,beta, cli, &(templ[ell*mtot*mtot*6+m1*mtot*6+m2*6]));
        Rq[IDX_R(egl,ell,m1,m2)] = tmp;
        Rq[IDX_R(egl,ell,m2,m1)] = Rq[IDX_R(egl,ell,m1,m2)];
      }
    }
  }

  return;

}

CREATE_PARAMETRIC_POL_TEMPLATE_FILE_INIT(bleak,bleak_init);
