/*
 *  smica.c
 *  lowly_project
 *
 *  Created by Karim Benabed on 30/10/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#include "smica.h"
double smica_crit_classic(void *vsmic,error **err);
double kld(int n, double* rq_hat, double* rq, error **err);

void printMat(double* A, int n, int m) {
  int im,in;
  for(in=0;in<n;in++) {
    for(im=0;im<m-1;im++) {
      fprintf(stderr,"%g , ",A[in*m+im]);
    }
    fprintf(stderr,"%g\n",A[in*m+m-1]);
  }
}
// General funcs

Smica* Smica_init(int nq, double *wq, int m, double *rq_hat, double* rq_0, int nc, SmicaComp **SC,error **err) {
  Smica* smic;
  int isc,iq,info;
  int trois;
  char uplo,diag;
  
  //_DEBUGHERE_("","");
  
  smic = malloc_err(sizeof(Smica),err);
  forwardError(*err,__LINE__,NULL);
  
  //_DEBUGHERE_("","");
  smic->nq = nq;
  smic->m = m;
  trois = 3;
  if(rq_0 == NULL) {
    trois = 2;
  }
  //_DEBUGHERE_("","");
  
  smic->rq_hat = malloc_err(sizeof(double)*(trois*nq+1)*m*m, err);
  forwardError(*err,__LINE__,NULL);
  memcpy(smic->rq_hat,rq_hat,sizeof(double)*m*m*nq);
  
  smic->z_buf = smic->rq_hat + m*m*nq;
  smic->rq = smic->z_buf + m*m;
  //_DEBUGHERE_("","");
  if (rq_0!=NULL) {
    smic->rq_0 = smic->rq + m*m*nq;
    memcpy(smic->rq_0,rq_0,sizeof(double)*m*m*nq);
    //_DEBUGHERE_("","");
  } else {
    smic->rq_0 = NULL;
    //_DEBUGHERE_("","");
  }
  //_DEBUGHERE_("","");
  
  smic->nc = nc;
  smic->SC = malloc_err(sizeof(SmicaComp*)*nc, err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");
  smic->offset_nc = malloc_err(sizeof(int)*nc, err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");
  
  smic->offset_nc[0] = 0;
  smic->SC[0] = SC[0];
  //_DEBUGHERE_("","");
  for(isc=1;isc<nc;isc++) {
    //_DEBUGHERE_("","");
    smic->offset_nc[isc] = smic->offset_nc[isc-1] + smic->SC[isc-1]->ndim;
    testErrorRetVA(m!=SC[isc]->m,smica_uncomp,"uncompatible number of band in component %d (got %d expected %d)",*err,__LINE__,NULL,isc,SC[isc]->m,m);
    testErrorRetVA(nq!=SC[isc]->nq,smica_uncomp,"uncompatible number of bins in component %d (got %d expected %d)",*err,__LINE__,NULL,isc,SC[isc]->nq,nq);
    
    //_DEBUGHERE_("","");
    smic->SC[isc] = SC[isc];
  }
  //_DEBUGHERE_("","");
  
  smic->wq = malloc_err(sizeof(double)*nq,err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");
  
  if (wq==NULL) {
    for(iq=0;iq<nq;iq++) {
      //_DEBUGHERE_("","");
      smic->wq[iq] = 1;
    }
  } else {
    //_DEBUGHERE_("","");
    memcpy(smic->wq,wq,sizeof(double)*nq);
  }
  //_DEBUGHERE_("","");
  
  smic->crit = &smica_crit_classic;
  smic->crit_classic_init = 0;
  smic->crit_cor = NULL;
  smic->gvec = NULL;

  smic->eig_buf = NULL;
  smic->eig_lwork = 0;
  smic->eig_nrm = 0;

  smic->lkl_data = NULL;
  smic->lkl_data_free = NULL;

  smic->cnt=0;
  return smic;
  
}

void free_Smica(void **psmic) {
  Smica *smic;
  int isc;
  
  smic = *psmic;
  
  
  free(smic->wq);
  free(smic->rq_hat);
  free(smic->offset_nc);
  for(isc=0;isc<smic->nc;isc++) {
    if (smic->SC[isc]->names!=NULL) {
      free(smic->SC[isc]->names);
    }
    if (smic->SC[isc]->free!=NULL) {
      smic->SC[isc]->free((void**)&(smic->SC[isc]));          
    }
  }
  free(smic->SC);
    
  if (smic->gvec!=NULL) {
    free(smic->gvec);
  }
  if (smic->quad_mask!=NULL) {
    free(smic->quad_mask);
  }
  if (smic->crit_cor!=NULL) {
    free(smic->crit_cor);
  }
  
  if (smic->eig_buf!=NULL) {
    free(smic->eig_buf);
    free(smic->eig_nrm);
  }
  
  if (smic->lkl_data_free!=NULL) {
    smic->lkl_data_free(&(smic->lkl_data));
  }
  free(smic);
  *psmic=NULL;
}

void Smica_fg(void* vsmic, double* pars, double* fgvec, error **err) {
  int isc;
  Smica *smic;
  double res;
  int iq,i,j;
  int iv;
  smic = vsmic;

  // init rq matrix
  memset(smic->rq,0,sizeof(double)*smic->m*smic->m*smic->nq);
  
  // update rq matrix according to each component
  for(isc=0;isc<smic->nc;isc++) {
    char nn[40];
    if (smic->SC[isc]->isfg==0) {
      //_DEBUGHERE_("jump %d",isc);
      continue;
    }
    //_DEBUGHERE_("comp %d update (off %d)",isc,smic->offset_nc[isc]);
    //printMat(smic->rq, smic->m, smic->m);
    //_DEBUGHERE_("%g",*(pars+smic->offset_nc[isc]));
    smic->SC[isc]->update(smic->SC[isc],pars+smic->offset_nc[isc], smic->rq, err);
    forwardError(*err,__LINE__,);
    //sprintf(nn,"pq_%d.la",isc);
    //write_bin_vector(smic->rq, nn, sizeof(double)*(smic->nq*smic->m*smic->m), err);  
    //forwardError(*err,__LINE__,);
  
    //_DEBUGHERE_("comp %d update done",isc);
    //printMat(smic->rq, smic->m, smic->m);
  }
  
  //symetrize matrix
  
  for (iq=0;iq<smic->nq;iq++) {
    double *rq;
    rq = smic->rq+iq*smic->m*smic->m;
    for(i=0;i<smic->m;i++) {
      for (j=0;j<i;j++) {
        rq[j+i*smic->m] = rq[j*smic->m+i];
      }
    }
  }
  //write_bin_vector(smic->rq, "rq.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);  
  //write_bin_vector(smic->rq_hat, "rqhat.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);   
  // ici calculer la vraissemblance a partir de rq et rq_hat

  testErrorRet(smic->crit!=smica_crit_gauss,-24324,"not implemented",*err,__LINE__,);

  for (iv=0;iv<smic->quad_sn;iv++) {
    int civ;
    fgvec[iv] = smic->rq[smic->quad_mask[iv]];
  }
  //_DEBUGHERE_("%g %g %g %g",fgvec[0],fgvec[1],fgvec[smic->quad_sn-2],fgvec[smic->quad_sn-1]);
}

void Smica_data(void* vsmic, double* fgvec, error **err) {
  int isc;
  Smica *smic;
  double res;
  int iq,i,j;
  int iv;
  smic = vsmic;

  testErrorRet(smic->crit!=&smica_crit_gauss,-24324,"not implemented",*err,__LINE__,);

  for (iv=0;iv<smic->quad_sn;iv++) {
    int civ;
    fgvec[iv] = smic->rq_hat[smic->quad_mask[iv]];
  }
  //_DEBUGHERE_("%g %g %g %g",fgvec[0],fgvec[1],fgvec[smic->quad_sn-2],fgvec[smic->quad_sn-1]);
}

int Smica_vecsize(void* vsmic, error **err) {
  int isc;
  Smica *smic;
  double res;
  int iq,i,j;
  int iv;
  smic = vsmic;

  //_DEBUGHERE_("%p %p %p",smic,smic->crit,&smica_crit_gauss);
  
  testErrorRet(smic->crit!=&smica_crit_gauss,-24324,"not implemented",*err,__LINE__,0);

  return smic->quad_sn;
}

void Smica_gcal(void* vsmic, double* pars, double* fgvec, error **err) {
  int isc;
  Smica *smic;
  double res;
  int iq,i,j;
  int iv;
  smic = vsmic;

  // init rq matrix
  for(i=0;i<smic->nq*smic->m*smic->m;i++) {
    smic->rq[i] = 1;
  }
  
  // update rq matrix according to each component
  for(isc=0;isc<smic->nc;isc++) {
    char nn[40];
    if (smic->SC[isc]->ismul==0) {
      //_DEBUGHERE_("jump %d",isc);
      continue;
    }
    //_DEBUGHERE_("comp %d update (off %d)",isc,smic->offset_nc[isc]);
    //printMat(smic->rq, smic->m, smic->m);
    //_DEBUGHERE_("%g",*(pars+smic->offset_nc[isc]));
    smic->SC[isc]->update(smic->SC[isc],pars+smic->offset_nc[isc], smic->rq, err);
    forwardError(*err,__LINE__,);
    //sprintf(nn,"pq_%d.la",isc);
    //write_bin_vector(smic->rq, nn, sizeof(double)*(smic->nq*smic->m*smic->m), err);  
    //forwardError(*err,__LINE__,-1);
  
    //_DEBUGHERE_("comp %d update done",isc);
    //printMat(smic->rq, smic->m, smic->m);
  }
  
  //symetrize matrix
  
  for (iq=0;iq<smic->nq;iq++) {
    double *rq;
    rq = smic->rq+iq*smic->m*smic->m;
    for(i=0;i<smic->m;i++) {
      for (j=0;j<i;j++) {
        rq[j+i*smic->m] = rq[j*smic->m+i];
      }
    }
  }
  //write_bin_vector(smic->rq, "rq.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);  
  //write_bin_vector(smic->rq_hat, "rqhat.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);   
  // ici calculer la vraissemblance a partir de rq et rq_hat

  testErrorRet(smic->crit!=smica_crit_gauss,-24324,"not implemented",*err,__LINE__,);

  for (iv=0;iv<smic->quad_sn;iv++) {
    int civ;
    fgvec[iv] = smic->rq[smic->quad_mask[iv]];
  }
  //_DEBUGHERE_("%g %g %g %g",fgvec[0],fgvec[1],fgvec[smic->quad_sn-2],fgvec[smic->quad_sn-1]);
}

//#define TIMER_MSEC(TIMER_after,TIMER_before) (((TIMER_after.tv_sec-TIMER_before.tv_sec)*1000000+(TIMER_after.tv_usec-TIMER_before.tv_usec))/1000)
double Smica_lkl(void* vsmic, double* pars, error **err) {
  int isc;
  Smica *smic;
  double res;
  int iq,i,j;
  //struct timeval starttime1,starttime2,starttime3,starttime4,endtime;

  smic = vsmic;

  //gettimeofday(&starttime1,NULL);
  // init rq matrix
  if (smic->rq_0 == NULL) {
    memset(smic->rq,0,sizeof(double)*smic->m*smic->m*smic->nq);
  } else {
    memcpy(smic->rq,smic->rq_0,sizeof(double)*smic->m*smic->m*smic->nq);
  }
  
  //gettimeofday(&starttime2,NULL);
 
  // update rq matrix according to each component
  for(isc=0;isc<smic->nc;isc++) {
    char nn[40];
    //_DEBUGHERE_("comp %d update (off %d) %d %d %d",isc,smic->offset_nc[isc],smic->nq,smic->m,smic->nq*smic->m*smic->m);
    //printMat(smic->rq, smic->m, smic->m);
    //_DEBUGHERE_("%g",*(pars+smic->offset_nc[isc]));
    smic->SC[isc]->update(smic->SC[isc],pars+smic->offset_nc[isc], smic->rq, err);
    forwardError(*err,__LINE__,0);
    //sprintf(nn,"pq_%d.la",isc);
    //write_bin_vector(smic->rq, nn, sizeof(double)*(smic->nq*smic->m*smic->m), err);  
    //forwardError(*err,__LINE__,-1);
  
    //_DEBUGHERE_("comp %d update done",isc);
    //printMat(smic->rq, smic->m, smic->m);
  }
  
  //symetrize matrix
  
  //gettimeofday(&starttime3,NULL);
  for (iq=0;iq<smic->nq;iq++) {
    double *rq;
    rq = smic->rq+iq*smic->m*smic->m;
    for(i=0;i<smic->m;i++) {
      for (j=0;j<i;j++) {
        rq[j+i*smic->m] = rq[j*smic->m+i];
      }
    }
  }
  //write_bin_vector(smic->rq, "rq.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);  
  //write_bin_vector(smic->rq_hat, "rqhat.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);   
  // ici calculer la vraissemblance a partir de rq et rq_hat
  
  //gettimeofday(&starttime4,NULL);
  res = smic->crit(smic,err);
  forwardError(*err,__LINE__,0);
  //gettimeofday(&endtime,NULL);
 
  //testErrorRet(smic->cnt==1,-24324,"not implemented",*err,__LINE__,0);
  smic->cnt+=1;
  //_DEBUGHERE_("init %ld, update %ld, sym %ld,crit %ld, total %ld",TIMER_MSEC(starttime2,starttime1),TIMER_MSEC(starttime3,starttime2),TIMER_MSEC(starttime4,starttime3),TIMER_MSEC(endtime,starttime4),TIMER_MSEC(endtime,starttime1));
  return res;
}


//gaussian approx criterion

double smica_crit_gauss(void *vsmic, error **err) {
  int iq,im1,im2,iv,m2,m;
  double res,les;
  Smica *smic;
  char uplo;
  double done,dzero;
  int one,i,j;
  //struct timeval starttime1, starttime2,starttime3, endtime;

  smic = vsmic;
  //_DEBUGHERE_("%d %d",smic->nq,smic->m);
  //write_bin_vector(smic->rq, "rq.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);   
  //forwardError(*err,__LINE__,0);
  //write_bin_vector(smic->rq_hat, "rq_hat.dat", sizeof(double)*(smic->nq*smic->m*smic->m), err);   
  //forwardError(*err,__LINE__,0);
  // reorganize data
  //write_bin_vector(smic->quad_mask, "quad_mask.dat", sizeof(int)*(smic->quad_sn), err);   
  //gettimeofday(&starttime1,NULL);
  m = smic->m;
  m2 = m*m;
  iv = 0;
  for (iv=0;iv<smic->quad_sn;iv++) {
    int civ;
    smic->gvec[iv] = smic->rq[smic->quad_mask[iv]] - smic->rq_hat[smic->quad_mask[iv]];
    civ = smic->quad_mask[iv];
    //_DEBUGHERE_("%d (%d %d %d) %d %g %g %g",iv,civ/m2,(civ-(civ/m2)*m2)/m,(civ-(civ/m2)*m2)%m, smic->quad_mask[iv],smic->rq[smic->quad_mask[iv]], smic->rq_hat[smic->quad_mask[iv]],smic->gvec[iv]);
  }
  //_DEBUGHERE_("%d",iv);
  //gettimeofday(&starttime2,NULL);
 
  //one = 1;
  //done = 1;
  //dzero = 0;
  //uplo = 'L';
  //printMat(smic->crit_cor,smic->quad_sn,smic->quad_sn);
  //_DEBUGHER_("%d",smic->quad_sn);
  //write_bin_vector(smic->gvec, "gvec.dat", sizeof(double)*(smic->quad_sn), err);   
  //write_bin_vector(smic->crit_cor, "crit_cor.dat", sizeof(double)*(smic->quad_sn)*(smic->quad_sn), err);   
  //_DEBUGHERE_("","");
  //for(i=0;i<smic->quad_sn;i++) {
  //  smic->gvec[smic->quad_sn+i] = smic->crit_cor[i*smic->quad_sn]*smic->gvec[0];
  //  for(j=1;j<smic->quad_sn;j++) {
  //    smic->gvec[smic->quad_sn+i] += smic->crit_cor[i*smic->quad_sn+j]*smic->gvec[j];
  //  }
  //}
  //dsymv(&uplo, &smic->quad_sn, &done, smic->crit_cor, &smic->quad_sn, smic->gvec, &one, &dzero, smic->gvec+smic->quad_sn, &one);
  //write_bin_vector(smic->gvec, "gvecCC.dat", sizeof(double)*(smic->quad_sn), err);   
   //gettimeofday(&starttime3,NULL);
 /*
   res = 0;
  for(iq=0;iq<iv;iq++) {
    //_DEBUGHERE_("%g %g %g",res,smic->gvec[iq],smic->gvec[iq+iv])
    res += smic->gvec[iq]*smic->gvec[iq+iv];
  }
  */
  res = 0;
  int qsn = smic->quad_sn;
  double *crit_cor = smic->crit_cor;
  double *gvec = smic->gvec;
#pragma omp parallel for private(les,i,j) reduction(+:res) firstprivate(crit_cor,gvec) schedule(dynamic,512)
  //_DEBUGHERE_("tt","");
  for(i=0;i<qsn;i++) {
    les = crit_cor[i*qsn+i]*gvec[i];
    for(j=i+1;j<qsn;j++) {
      les += 2*gvec[j]*crit_cor[i*qsn+j];
    }
    res += les*gvec[i];
  }

  //gettimeofday(&endtime,NULL);
  //_DEBUGHERE_("select %ld, half product %ld, finish %ld, total %ld",TIMER_MSEC(starttime2,starttime1),TIMER_MSEC(starttime3,starttime2),TIMER_MSEC(endtime,starttime3),TIMER_MSEC(endtime,starttime1));
 return -.5*res;  
}

void smica_set_crit_gauss(Smica *smic, double *crit_cor, int *mask,int *ordering,error **err) {
  int iv,jv,iq;
  int nv;


  nv = (smic->nq * smic->m * (smic->m+1))/2;
  smic->quad_mask = malloc_err(sizeof(int)*smic->nq*smic->m*smic->m,err);
  forwardError(*err,__LINE__,);

  // mask is both in spectra an ell !!!
  nv = 0;
  if (ordering == NULL) {
    for (iv=0;iv<smic->m;iv++) {
      for (jv=iv;jv<smic->m;jv++) {
        for (iq=0;iq<smic->nq;iq++) {
          //_DEBUGHERE_("%d %d %d %d (%d)",iq,iv,jv,mask[iq*smic->m*smic->m+iv*smic->m+jv],nv);
          if (mask == NULL || mask[iq*smic->m*smic->m+iv*smic->m+jv]!=0) {
            smic->quad_mask[nv] = iq*smic->m*smic->m+iv*smic->m+jv;
            nv++;
          }
        }
      }
    }
  } else {
    int ic;
    for (ic=0;ic<(smic->m*(smic->m+1))/2;ic++) {
      iv = ordering[ic*2];
      jv = ordering[ic*2+1];
      for (iq=0;iq<smic->nq;iq++) {
        if (mask == NULL || mask[iq*smic->m*smic->m+iv*smic->m+jv]!=0) {
          //_DEBUGHERE_("%d %d %d %d (%d)",iq,iv,jv,mask[iq*smic->m*smic->m+iv*smic->m+jv],nv);
          smic->quad_mask[nv] = iq*smic->m*smic->m+iv*smic->m+jv;
          nv++;
        }
      }
    }
  }

  smic->gvec = malloc_err(sizeof(double)*nv*2,err);
  forwardError(*err,__LINE__,);
  smic->quad_sn = nv;

  //write_bin_vector(smic->quad_mask, "quad_mask.dat", sizeof(int)*(smic->quad_sn), err);

  smic->crit_cor = malloc_err(sizeof(double)*nv*nv,err);
  forwardError(*err,__LINE__,);
  
  memcpy(smic->crit_cor,crit_cor,sizeof(double)*nv*nv);

  smic->crit = &smica_crit_gauss;
}


// Simple components

// constant

void comp_cst_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  double * rq_0;
  int sz,one;
  double done;

  SC = data;
  rq_0 = SC->data;
  sz = SC->m*SC->m*SC->nq;
  done = 1;
  one = 1;

  daxpy(&sz,&done,rq_0,&one, rq, &one);
}

void free_comp_cst(void** data) {
  SmicaComp *SC;
  
  SC = *data;
  free(SC->data);
  free(SC);
  *data = NULL;
}

SmicaComp* comp_cst_init(int nq, int m, double *rq_0, error **err) {
  SmicaComp *SC;
  double *data;

  data = malloc_err(sizeof(double)*m*m*nq,err);
  forwardError(*err,__LINE__,NULL);

  memcpy(data,rq_0,sizeof(double)*m*m*nq);

  SC = alloc_SC(0,nq,m,data,&comp_cst_update,&free_comp_cst,err);
  forwardError(*err,__LINE__,NULL);

  return SC;
}


// 1D

void comp_1D_AAt(int m, double *A, double *AAt, error **err) {
  int mm,one;
  char uplo;
  double done;
    
  memset(AAt,0,m*m*sizeof(double));
  
  mm = m;
  uplo = 'L';
  done = 1;
  one = 1;
  // -> AAt = 1 * A * A' + AAt
  //printMat(A, m, 1);
  dsyr(&uplo, &mm, &done, A, &one, AAt, &mm);  
  //printMat(AAt, m, m);
}





SmicaComp* alloc_SC(int ndim,int nq,int m,void* data, update_rq* update, posterior_log_free* pfree, error **err) {
  SmicaComp* SC;
  
  SC = malloc_err(sizeof(SmicaComp), err);
  forwardError(*err,__LINE__,NULL);
  SC->m = m;
  SC->nq = nq;
  SC->ndim = ndim;
  SC->update = update;
  SC->data = data;
  SC->free = pfree;
  SC->names=NULL;
  SC_set_compname(SC,"UNK");
  SC->isfg = 0;
  SC->ismul = 0;
  return SC;
}

void SC_set_compname(SmicaComp *SC, char *name) {
  sprintf(SC->comp_name,"%s",name);
}

void SC_isfg(SmicaComp *SC) {
  SC->isfg = 1;
}
void SC_ismul(SmicaComp *SC) {
  SC->ismul = 1;
}

void SC_setnames(SmicaComp *SC, char** names, error **err) {
  int i;
  if (SC->names!=NULL) {
    free(SC->names);
  }
  if (SC->ndim!=0) {
   SC->names = malloc_err(sizeof(_smicanames)*SC->ndim,err);
   forwardError(*err,__LINE__,);  
  } else{
    SC->names = malloc_err(sizeof(_smicanames)*1,err);
   forwardError(*err,__LINE__,);  
  }
  for(i=0;i<SC->ndim;i++) {
    sprintf(SC->names[i],"%s",names[i]);
  }
}


SmicaComp* comp_1D_init(int nq, int m, double *A, error **err) {
  SmicaComp *SC;
  char uplo;
  double done;
  int one;
  SC_1D_data *data;
  int ndim;
  
  data = malloc_err(sizeof(SC_1D_data), err);
  forwardError(*err,__LINE__,NULL);
  
  data->AAt = malloc_err(sizeof(double)*m*m, err);
  forwardError(*err,__LINE__,NULL);
  
  if (A!=NULL) {
    data->Acst=1;
    ndim = nq;

    comp_1D_AAt(m, A,data->AAt, err);    
    forwardError(*err,__LINE__,NULL);
    
  } else {
    data->Acst=0;
    ndim = nq + m;    
  }
  
  SC = alloc_SC(ndim,nq,m,data,&comp_1D_update,&free_comp_1D,err);
  forwardError(*err,__LINE__,NULL);

  return SC;
}




void free_comp_1D(void** data) {
  SmicaComp *SC;
  
  SC = *data;
  //_DEBUGHERE_("","");
  free(((SC_1D_data*) SC->data)->AAt);
  //_DEBUGHERE_("","");
  free(SC->data);
  //_DEBUGHERE_("","");
  free(SC);
  //_DEBUGHERE_("","");
  *data = NULL;
  //_DEBUGHERE_("","");
}




void comp_1D_update(void* data,double* locpars, double* rq, error **err) {
  int iq,one;
  SmicaComp *SC;
  double *AAt;
  int m,m2;
  double *mpars;
  SC_1D_data *SCdat;
  
  SC = data;
  m = SC->m;
  SCdat = SC->data;
  AAt = SCdat->AAt;
  m2 = m*m;
  one = 1;
  
  mpars = locpars;
  if (SCdat->Acst==0) {
    double *A;
    A = locpars;
    mpars = locpars + m;
    comp_1D_AAt(m,A,AAt, err);    
    forwardError(*err,__LINE__,);
  }
  
  //printMat(AAt, m, m);

  for(iq=0;iq<SC->ndim;iq++) {
    //printMat(rq+m2*iq, m, m);
    daxpy(&m2,mpars + iq,AAt,&one, rq+m2*iq, &one);
    //printMat(rq+m2*iq, m, m);
    // -> rq[iq*m2] = locpars[iq] * AAt + rq[iq*m2]
  }
}



// nD

SmicaComp* comp_nD_init(int nq, int m, int nd, double *A, error **err) {
  SmicaComp *SC;
  char uplo;
  double done;
  int one;
  int ndim;
  void* data;
  
  
  //_DEBUGHERE_("","");
  data = malloc_err(sizeof(SC_nD_data), err);
  forwardError(*err,__LINE__,NULL);
  
  //_DEBUGHERE_("","");
  ((SC_nD_data*) data)->nd = nd;

  //_DEBUGHERE_("","");
  ((SC_nD_data*) data)->A = malloc_err(sizeof(double)*(m*nd*2+nd*nd), err);
  forwardError(*err,__LINE__,NULL);
  
  //_DEBUGHERE_("","");
  if (A!=NULL) {
    ndim = nq*(nd*(nd+1))/2;
    ((SC_nD_data*) data)->Acst=1;
    //_DEBUGHERE_("","");
    
    memcpy(((SC_nD_data*) data)->A,A,m*nd*sizeof(double));    
  } else {
    ndim = m*nd+nq*(nd*(nd+1))/2;
    ((SC_nD_data*) data)->Acst=0;
    //_DEBUGHERE_("","");
    
  }
  //_DEBUGHERE_("","");
  ((SC_nD_data*) data)->Ab = ((SC_nD_data*) data)->A + m*nd;
  ((SC_nD_data*) data)->P = ((SC_nD_data*) data)->Ab + m*nd;
  //_DEBUGHERE_("","");
  
  // beware A is C oriented
  //_DEBUGHERE_("","");
  SC = alloc_SC(ndim,nq,m,data,&comp_nD_update,&free_comp_nD,err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");
  
  return SC;
}




void free_comp_nD(void** data) {
  SmicaComp *SC;
  
  SC = *data;
  //_DEBUGHERE_("","");
  free(((SC_nD_data*) SC->data)->A);
  //_DEBUGHERE_("","");
  free(SC->data);
  //_DEBUGHERE_("","");
  free(SC);
  //_DEBUGHERE_("","");
  *data = NULL;
}




void comp_nD_update(void* data,double* locpars, double* rq, error **err) {
  int iq,one;
  SmicaComp *SC;
  double *A,*Ab,*P,*Ppack;
  int m,m2,nd,nq;
  double done,dzero;
  double *mpars;
  char transa,transb,side,uplo;
  
  // locpar is a (C oriented) q*nd*nd U Triangular matrix
  
  SC = data;
  m = SC->m;
  nq = SC->nq;
  
  nd = ((SC_nD_data*) SC->data)->nd;
  Ab = ((SC_nD_data*) SC->data)->Ab;
  P = ((SC_nD_data*) SC->data)->P;
  if (((SC_nD_data*) SC->data)->Acst==1) {
    A = ((SC_nD_data*) SC->data)->A;
    mpars = locpars;
  } else {
    A = locpars;
    mpars = locpars + m*nd;
  }
  
  m2 = m*m;
  one = 1;
  
  done = 1;
  dzero = 0;
  transb = 'N'; // because A is C ordered
  
  for(iq=0;iq<nq;iq++) {
    int ii;
    int ix,iy;
    
    // unpack locpar[q]
    Ppack = locpars + iq*((nd*(nd+1))/2);
    ii = 0;
    for(ix=0;ix<nd;ix++) {
      for(iy=ix;iy<nd;iy++) {
        P[ix*nd+iy] = Ppack[ii];
        //P[iy*nd+ix] = Ppack[ii];
        //P[ix*nd+iy] = 0;
        //P[iy*nd+ix] = 0;
        ii++;
      }
      //P[ix*nd+ix] = 1;
    }
    //printMat(P,2,2);

    //_DEBUGHERE_("P %d",iq);
    //printMat(P, nd, nd);
    //_DEBUGHERE_("A","");
    //printMat(A, m, nd);
    //_DEBUGHERE_("P.A'","");
    //printMat(Ab, m, nd);    
    transa = 'N';
    side = 'L';
    uplo = 'L';
    // Ab = P.A' (Ab is fortran ordered while A is C ordered)
    dsymm(&side, &uplo, &nd, &m, &done, P, &nd, A, &nd, &dzero, Ab, &nd);
    //dgemm(&transa, &transb, &nd, &m, &nd, &done, P, &nd, A, &nd, &dzero, Ab, &nd);

    //_DEBUGHERE_("P %d",iq);
    //printMat(P, nd, nd);
    //_DEBUGHERE_("A","");
    //printMat(A, m, nd);    
    //_DEBUGHERE_("P.A'","");
    //printMat(Ab, m, nd);    
    
    // Rq += A.Ab (=A.P.A')
    transa = 'T'; // because A is c ordered

    //_DEBUGHERE_("rq a","");
    //printMat(rq+m2*iq, m, m);    
    dgemm(&transa, &transb, &m, &m, &nd, &done, A, &nd, Ab, &nd, &done, rq+m2*iq, &m);
    //_DEBUGHERE_("rq b","");
    //printMat(rq+m2*iq, m, m);    
    //printMat(rq+m2*iq,2,2);
  }
}



// CMB

SmicaComp * comp_CMB_init(int nbins, int mt,int mp, int *has_cl, double* Acprs, error **err) {
  double *A;
  SC_CMB_data* data;
  int mtot;
  int trois,im,i,six;
  SmicaComp *SC;
  
  data = malloc_err(sizeof(SC_CMB_data), err);
  forwardError(*err,__LINE__,NULL);
  
  trois = has_cl[0]+has_cl[1]+has_cl[2];
  testErrorRet(trois==0,smica_uncomp,"mismatch",*err,__LINE__,NULL);
  six = trois + has_cl[3]+has_cl[4]+has_cl[5];
  
  mtot = mt*has_cl[0]+(has_cl[1]+has_cl[2])*mp;
  testErrorRet(mtot==0,smica_uncomp,"mismatch",*err,__LINE__,NULL);
  
  testErrorRet(mt!=0 && has_cl[0]==0,smica_uncomp,"mismatch",*err,__LINE__,NULL);
  testErrorRet(mt==0 && has_cl[0]!=0,smica_uncomp,"mismatch",*err,__LINE__,NULL);

  testErrorRet(mp!=0 && has_cl[1]==0 && has_cl[2]==0,smica_uncomp,"mismatch",*err,__LINE__,NULL);
  testErrorRet(mp==0 && (has_cl[1]!=0 || has_cl[2]!=0),smica_uncomp,"mismatch",*err,__LINE__,NULL);

  testErrorRet(has_cl[0]==0 && (has_cl[3]!=0 || has_cl[4]!=0),smica_uncomp,"mismatch",*err,__LINE__,NULL);
  testErrorRet(has_cl[1]==0 && (has_cl[3]!=0 || has_cl[5]!=0),smica_uncomp,"mismatch",*err,__LINE__,NULL);
  testErrorRet(has_cl[2]==0 && (has_cl[4]!=0 || has_cl[5]!=0),smica_uncomp,"mismatch",*err,__LINE__,NULL);

  if (trois==1) {
    // cas particulier 1D
    data->locpars=NULL;

    data->SCnD = comp_1D_init(nbins,mtot,Acprs,err);
    forwardError(*err,__LINE__,NULL);
  
    SC = alloc_SC(nbins,nbins,mtot,data,&comp_CMB_update,&free_comp_CMB,err);
    forwardError(*err,__LINE__,NULL);
    return SC;
    
  }
  
  A = malloc_err(sizeof(double)*mtot*trois, err);
  forwardError(*err,__LINE__,NULL);
  memset(A,0,sizeof(double)*mtot*trois);
  // fill T
  if(has_cl[0]==1) {
    for(im=0;im<mt;im++) {
      A[im*trois] = Acprs[im];
    }
  }
  // fill E and B
  for(im=0;im<mp;im++) {
    if (has_cl[1]==1) {
      A[(im+mt*has_cl[0])*trois+has_cl[0]] = Acprs[im+mt];      
    }
    if (has_cl[2]==1) {
      A[(im+mt*has_cl[0]+has_cl[1]*mp)*trois+has_cl[0]+has_cl[1]] = Acprs[im+mt];      
    }
  }
  
  data->locpars = malloc_err(sizeof(double)*nbins*(trois*(trois+1))/2,err);
  forwardError(*err,__LINE__,NULL);
  
  data->SCnD = comp_nD_init(nbins,mtot,trois,A,err);
  forwardError(*err,__LINE__,NULL);

  free(A);

  for(i=0;i<6;i++) {
    data->has_cl[i] = has_cl[i];
    data->jmp_cl[i] = -1;
  }
  
  if (trois==3) {
    data->jmp_cl[0] = 0;
    data->jmp_cl[1] = 3;
    data->jmp_cl[2] = 5;
    if (has_cl[3]==1) {
      data->jmp_cl[3] = 1;
    }
    if (has_cl[4]==1) {
      data->jmp_cl[4] = 2;
    }
    if (has_cl[5]==1) {
      data->jmp_cl[5] = 4;
    }
  }
  
  if (trois==2) {
    if (has_cl[0]==1) {
      data->jmp_cl[0] = 0;
      if (has_cl[1]==1) {
        data->jmp_cl[1] = 2;
        if (has_cl[3]==1) {
          data->jmp_cl[3] = 1;  
        }
      } else {
        data->jmp_cl[2] = 2;
        if (has_cl[4]==1) {
          data->jmp_cl[4] = 1;
        }
      }
    } else {
      data->jmp_cl[1] = 0;
      data->jmp_cl[2] = 2;
      if (has_cl[5]==1) {
        data->jmp_cl[5] = 1;
      }
    }
  }
  data->trois = trois;
  
  SC = alloc_SC(nbins*six,nbins,mtot,data,&comp_CMB_update,&free_comp_CMB,err);
  forwardError(*err,__LINE__,NULL);
  return SC;
  
}




void comp_CMB_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  SC_CMB_data *SCd;
  int td,ic,im,i,iq;
  
  SC = data;
  SCd = SC->data;
  
  if (SCd->locpars==NULL) {
    // 1d case
    comp_1D_update(SCd->SCnD,locpars,rq,err);
    forwardError(*err,__LINE__,);
    return;
  } 
  
  td  = (SCd->trois*(SCd->trois+1))/2;
  
  // nd il faut que je reoriente mon vecteur
  i=0;
  for(ic=0;ic<6;ic++) {

    if(SCd->jmp_cl[ic]!=-1) {
      for(iq=0;iq<SC->nq;iq++) {
        //_DEBUGHERE_("%d %d %d %g %d",ic,iq,i, locpars[i],iq*td+SCd->jmp_cl[ic]);
        SCd->locpars[iq*td+SCd->jmp_cl[ic]] = locpars[i];
        i+=SCd->has_cl[ic];
      }
    }
  }
  comp_nD_update(SCd->SCnD,SCd->locpars,rq,err);
  forwardError(*err,__LINE__,);
  return;
  
}




void free_comp_CMB(void** data) {
  SmicaComp *SC;
  SC_CMB_data *SCd;
  
  SC = *data;
  SCd = SC->data;
  SCd->SCnD->free((void**)&SCd->SCnD);
  if(SCd->locpars!=NULL) {
    free(SCd->locpars);
  }
  free(SCd);
  free(SC);
  *data = NULL;
}




SmicaComp* comp_calTP_init(int q, int mT, int mP, int *TEB, int npar, int *im,double *w,int *other, error **err ) {
  SC_calTP *gc;
  SmicaComp *SC;
  int m;
  int i;

  gc = malloc_err(sizeof(SC_calTP),err);
  forwardError(*err,__LINE__,NULL);

  gc->npar = npar;


  m = mT*TEB[0] + mP *TEB[1] + mP*TEB[2];
  gc->mT = mT;
  gc->mP = mP;
  gc->TEB[0] = TEB[0];
  gc->TEB[1] = TEB[1];
  gc->TEB[2] = TEB[2];

  gc->calvec = malloc_err(sizeof(double)*m,err);
  forwardError(*err,__LINE__,NULL);
  
  for (i=0;i<m;i++) {
    gc->calvec[i]=1;
  }

  
  gc->im = malloc_err(sizeof(int)*npar,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->im,im,sizeof(int)*npar);

  gc->w = malloc_err(sizeof(double)*m*m*2,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->w,w,sizeof(double)*m*m*2);

  gc->other = malloc_err(sizeof(int)*m*m*2,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->other,other,sizeof(int)*m*m*2);

  SC = alloc_SC(npar,q,m,gc, &comp_calTP_update, &comp_calTP_free,err);
  forwardError(*err,__LINE__,NULL);
  SC_ismul(SC);
  return SC;

}

void comp_calTP_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  SC_calTP *gc;
  int i,iq,im1,im2,m,m2;

  SC = data;
  gc = SC->data;

  for(i=0;i<SC->ndim;i++) {
    int im;
    im = gc->im[i];
    gc->calvec[im] = exp(locpars[i]);
    if (gc->TEB[2]==1 && im>gc->mT) {
      gc->calvec[im+gc->mP*gc->TEB[1]] = gc->calvec[im];
    }
  }
  m = SC->m;
  m2 = m*m;
  
  for(iq=0;iq<SC->nq;iq++) {
    int iqo;
    iqo = iq*m2;
    for(im1=0;im1<SC->m;im1++) {
      int imo;
      double k1;
      k1 = gc->calvec[im1];
      imo = iqo+im1*m;
      for(im2=im1;im2<SC->m;im2++) {
        int mpos,im1_prime,im2_prime;
        double w,w_prime;
        mpos = (im1*m+im2)*2;
        //_DEBUGHERE_("%d %d %d %g %g %g",iq,im1,im2,k1,gc->calvec[im2],rq[imo+im2]);
        w = gc->w[mpos];
        w_prime = gc->w[mpos+1];
        im1_prime = gc->other[mpos];
        im2_prime = gc->other[mpos+1];
        //if (iq==0) {
        //  _DEBUGHERE_("%d | %d %d -> %g %g %g %g, %d %d -> %g %g %g %g,",mpos,im1,im2,w,gc->calvec[im1],gc->calvec[im2],w*gc->calvec[im1]*gc->calvec[im2],im1_prime,im2_prime,w_prime,gc->calvec[im1_prime],gc->calvec[im2_prime],w_prime*gc->calvec[im1_prime]*gc->calvec[im2_prime]);  
        //  _DEBUGHERE_("%g %g",w*gc->calvec[im1]*gc->calvec[im2]+w_prime*gc->calvec[im1_prime]*gc->calvec[im2_prime],gc->calvec[im1]*gc->calvec[im2]);
        // }
        rq[imo+im2] *= w*gc->calvec[im1]*gc->calvec[im2]+w_prime*gc->calvec[im1_prime]*gc->calvec[im2_prime];
        //rq[imo+im2] *= gc->calvec[im1]*gc->calvec[im2];
      } 
    }
  }
}

SmicaComp* comp_icalTP_init(int q, int mT, int mP, int *TEB, int npar, int *im,double *w,int *other, error **err ) {
  SC_calTP *gc;
  SmicaComp *SC;
  int m;
  int i;

  gc = malloc_err(sizeof(SC_calTP),err);
  forwardError(*err,__LINE__,NULL);

  gc->npar = npar;


  m = mT*TEB[0] + mP *TEB[1] + mP*TEB[2];
  gc->mT = mT;
  gc->mP = mP;
  gc->TEB[0] = TEB[0];
  gc->TEB[1] = TEB[1];
  gc->TEB[2] = TEB[2];

  gc->calvec = malloc_err(sizeof(double)*m,err);
  forwardError(*err,__LINE__,NULL);
  
  for (i=0;i<m;i++) {
    gc->calvec[i]=1;
  }

  
  gc->im = malloc_err(sizeof(int)*npar,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->im,im,sizeof(int)*npar);

  gc->w = malloc_err(sizeof(double)*m*m*2,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->w,w,sizeof(double)*m*m*2);

  gc->other = malloc_err(sizeof(int)*m*m*2,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->other,other,sizeof(int)*m*m*2);

  SC = alloc_SC(npar,q,m,gc, &comp_icalTP_update, &comp_calTP_free,err);
  forwardError(*err,__LINE__,NULL);
  SC_ismul(SC);
  return SC;

}

void comp_icalTP_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  SC_calTP *gc;
  int i,iq,im1,im2,m,m2;

  SC = data;
  gc = SC->data;

  for(i=0;i<SC->ndim;i++) {
    int im;
    im = gc->im[i];
    gc->calvec[im] = 1./sqrt(locpars[i]);
    if (gc->TEB[2]==1 && im>gc->mT) {
      gc->calvec[im+gc->mP*gc->TEB[1]] = gc->calvec[im];
    }
  }
  m = SC->m;
  m2 = m*m;
  
  for(iq=0;iq<SC->nq;iq++) {
    int iqo;
    iqo = iq*m2;
    for(im1=0;im1<SC->m;im1++) {
      int imo;
      double k1;
      k1 = gc->calvec[im1];
      imo = iqo+im1*m;
      for(im2=im1;im2<SC->m;im2++) {
        int mpos,im1_prime,im2_prime;
        double w,w_prime;
        mpos = (im1*m+im2)*2;
        //_DEBUGHERE_("%d %d %d %g %g %g",iq,im1,im2,k1,gc->calvec[im2],rq[imo+im2]);
        w = gc->w[mpos];
        w_prime = gc->w[mpos+1];
        im1_prime = gc->other[mpos];
        im2_prime = gc->other[mpos+1];
        //if (iq==0) {
        //  _DEBUGHERE_("%d | %d %d -> %g %g %g %g, %d %d -> %g %g %g %g,",mpos,im1,im2,w,gc->calvec[im1],gc->calvec[im2],w*gc->calvec[im1]*gc->calvec[im2],im1_prime,im2_prime,w_prime,gc->calvec[im1_prime],gc->calvec[im2_prime],w_prime*gc->calvec[im1_prime]*gc->calvec[im2_prime]);  
        //  _DEBUGHERE_("%g %g",w*gc->calvec[im1]*gc->calvec[im2]+w_prime*gc->calvec[im1_prime]*gc->calvec[im2_prime],gc->calvec[im1]*gc->calvec[im2]);
        // }
        if (w+w_prime!=0) {
          rq[imo+im2] *= w*gc->calvec[im1]*gc->calvec[im2]+w_prime*gc->calvec[im1_prime]*gc->calvec[im2_prime];
        }
        //rq[imo+im2] *= gc->calvec[im1]*gc->calvec[im2];
      } 
    }
  }
}

void comp_calTP_free(void** data) {
  SmicaComp *SC;
  SC_calTP *gc;
  
  SC = *data;
  gc = SC->data;

  
  free(gc->im);
  free(gc->calvec);
  free(gc->w);
  free(gc->other);
  free(gc);

  free(SC);
  
  *data = NULL;
}


void comp_beamTP_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  SC_beamTP *gc;
  int t,iq,im1,im2,m,m2,neigen,offm,offq;
  double cal;

  SC = data;
  gc = SC->data;
  memcpy(gc->pars+1,locpars,sizeof(double)*SC->ndim);

  m = SC->m;
  m2 = m*m;
  neigen = gc->neigen;
  
  for(iq=0;iq<SC->nq;iq++) {
    for(im1=0;im1<m;im1++) {
      for(im2=im1;im2<m;im2++) {
        cal = 0;
        offm = im1*m+im2;
        offq = iq*m2+ offm;
        for(t=0;t<neigen;t++) {
          cal += gc->pars[gc->im[offm*neigen + t]] * gc->modes[offq*neigen + t]; 
	        //if(iq==0) _DEBUGHERE_("%d %d %d -> %d %g %d %g | %g %g",im1,im2,t,gc->im[offm*neigen + t], gc->pars[gc->im[offm*neigen + t]],offq*neigen+t,gc->modes[offq*neigen + t],cal,exp(cal));
        }
        rq[offq] *= exp(2*cal);
	      rq[iq*m2 + im2*m + im1] = rq[offq];
      }
    }
  }
}

void comp_beamTP_free(void** data) {
  SmicaComp *SC;
  SC_beamTP *gc;
  
  SC = *data;
  gc = SC->data;

  
  free(gc->im);
  free(gc->pars);
  free(gc->modes);
  free(gc);

  free(SC);
  
  *data = NULL;
}

SmicaComp* comp_beamTP_init(int q, int mT, int mP, int *TEB, int npar, int *im,int neigen, double *modes,error **err ) {
  SC_beamTP *gc;
  SmicaComp *SC;
  int m;
  int i;

  gc = malloc_err(sizeof(SC_beamTP),err);
  forwardError(*err,__LINE__,NULL);
   
  m = mT*TEB[0] + mP *TEB[1] + mP*TEB[2];
  
  gc->pars = malloc_err(sizeof(double)*(npar+1),err);
  forwardError(*err,__LINE__,NULL);
  gc->pars[0] = 0;
   
  gc->neigen = neigen;
  gc->modes = malloc_err(sizeof(double)*neigen*m*m*q,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->modes,modes,sizeof(double)*neigen*m*m*q);  
  
  gc->im = malloc_err(sizeof(int)*neigen*m*m,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(gc->im,im,sizeof(int)*neigen*m*m);

  SC = alloc_SC(npar,q,m,gc, &comp_beamTP_update, &comp_beamTP_free,err);
  forwardError(*err,__LINE__,NULL);
  SC_ismul(SC);
  return SC;

}

void comp_totcal_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  int t,iq,im1,im2,m,m2,neigen,offm,offq;
  double cal;

  SC = data;
  cal = 1./(locpars[0]*locpars[0]);
  
  m = SC->m;
  m2 = m*m;
  
  for(iq=0;iq<SC->nq;iq++) {
    for(im1=0;im1<m;im1++) {
      for(im2=im1;im2<m;im2++) {
        rq[iq*m2 + im1*m + im2] *= cal;
        rq[iq*m2 + im2*m + im1] = rq[iq*m2 + im1*m + im2];
      }
    }
  }
}

void comp_totcal_free(void** data) {
  SmicaComp *SC;

  SC = *data;
  free(SC);
  
  *data = NULL;
}

SmicaComp* comp_totcal_init(int q, int mT, int mP, int *TEB,error **err ) {
  SC_beamTP *gc;
  SmicaComp *SC;
  int m;
  int i;

  m = mT*TEB[0] + mP *TEB[1] + mP*TEB[2];
  
  SC = alloc_SC(1,q,m,NULL, &comp_totcal_update, &comp_totcal_free,err);
  forwardError(*err,__LINE__,NULL);
  SC_ismul(SC);
  return SC;

}

void comp_totcalP_free(void** data) {
  SmicaComp *SC;

  SC = *data;
  free(SC->data);
  free(SC);
  
  *data = NULL;
}

void comp_totcalP_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  int t,iq,im1,im2,m,m2,neigen,offm,offq,im0;
  double cal,scal;
  int *mz;

  SC = data;
  mz = SC->data;

  cal = 1./(locpars[0]*locpars[0]);
  scal = 1./(locpars[0]);

  m = SC->m;
  m2 = m*m;
  
  for(iq=0;iq<SC->nq;iq++) {
    // TP
    for(im1=0;im1<mz[0];im1++) {
      im0 = im1;
      if (im0<mz[0]) {
        im0 = mz[0];
      }
      for(im2=im0;im2<m;im2++) {
        rq[iq*m2 + im1*m + im2] *= scal;
        rq[iq*m2 + im2*m + im1] = rq[iq*m2 + im1*m + im2];        
      }
    }
    //PP
    for(im1=mz[0];im1<m;im1++) {
      for(im2=im1;im2<m;im2++) {
        rq[iq*m2 + im1*m + im2] *= cal;
        rq[iq*m2 + im2*m + im1] = rq[iq*m2 + im1*m + im2];
      }
    }
  }
}


SmicaComp* comp_totcalP_init(int q, int mT, int mP, int *TEB,error **err ) {
  SC_beamTP *gc;
  SmicaComp *SC;
  int m;
  int i;
  int *mz;

  m = mT*TEB[0] + mP *TEB[1] + mP*TEB[2];
  
  mz = malloc_err(sizeof(int)*3,err);
  forwardError(*err,__LINE__,NULL);
  
  mz[0] = mT*TEB[0];
  mz[1] = mT*TEB[1];
  mz[2] = mT*TEB[2];

  SC = alloc_SC(1,q,m,mz, &comp_totcalP_update, &comp_totcalP_free,err);
  forwardError(*err,__LINE__,NULL);
  SC_ismul(SC);
  return SC;

}

void comp_totcalTP_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  int t,iq,im1,im2,m,m2,neigen,offm,offq,im0;
  double cal,scal;
  int *mz;

  SC = data;
  mz = SC->data;

  cal = 1./(locpars[0]*locpars[0]);
  scal = 1./(locpars[0]);

  m = SC->m;
  m2 = m*m;
  
  for(iq=0;iq<SC->nq;iq++) {
    // TP
    for(im1=0;im1<mz[0];im1++) {
      im0 = im1;
      if (im0<mz[0]) {
        im0 = mz[0];
      }
      for(im2=im0;im2<m;im2++) {
        rq[iq*m2 + im1*m + im2] *= scal;
        rq[iq*m2 + im2*m + im1] = rq[iq*m2 + im1*m + im2];        
      }
    }
    //PP
    //for(im1=mz[0];im1<m;im1++) {
    //  for(im2=im1;im2<m;im2++) {
    //    rq[iq*m2 + im1*m + im2] *= cal;
    //    rq[iq*m2 + im2*m + im1] = rq[iq*m2 + im1*m + im2];
    //  }
    //}
  }
}

SmicaComp* comp_totcalTP_init(int q, int mT, int mP, int *TEB,error **err ) {
  SC_beamTP *gc;
  SmicaComp *SC;
  int m;
  int i;
  int *mz;

  m = mT*TEB[0] + mP *TEB[1] + mP*TEB[2];
  
  mz = malloc_err(sizeof(int)*3,err);
  forwardError(*err,__LINE__,NULL);
  
  mz[0] = mT*TEB[0];
  mz[1] = mT*TEB[1];
  mz[2] = mT*TEB[2];

  SC = alloc_SC(1,q,m,mz, &comp_totcalTP_update, &comp_totcalP_free,err);
  forwardError(*err,__LINE__,NULL);
  SC_ismul(SC);
  return SC;

}

void comp_totcalPP_update(void* data,double* locpars, double* rq, error **err) {
  SmicaComp *SC;
  int t,iq,im1,im2,m,m2,neigen,offm,offq,im0;
  double cal,scal;
  int *mz;

  SC = data;
  mz = SC->data;

  cal = 1./(locpars[0]*locpars[0]);
  scal = 1./(locpars[0]);

  m = SC->m;
  m2 = m*m;
  
  for(iq=0;iq<SC->nq;iq++) {
    // TP
    //for(im1=0;im1<mz[0];im1++) {
    //  im0 = im1;
    //  if (im0<mz[0]) {
    //    im0 = mz[0];
    //  }
    //  for(im2=im0;im2<m;im2++) {
    //    rq[iq*m2 + im1*m + im2] *= scal;
    //    rq[iq*m2 + im2*m + im1] = rq[iq*m2 + im1*m + im2];        
    //  }
    //}
    //PP
    for(im1=mz[0];im1<m;im1++) {
      for(im2=im1;im2<m;im2++) {
        rq[iq*m2 + im1*m + im2] *= cal;
        rq[iq*m2 + im2*m + im1] = rq[iq*m2 + im1*m + im2];
      }
    }
  }
}

SmicaComp* comp_totcalPP_init(int q, int mT, int mP, int *TEB,error **err ) {
  SC_beamTP *gc;
  SmicaComp *SC;
  int m;
  int i;
  int *mz;

  m = mT*TEB[0] + mP *TEB[1] + mP*TEB[2];
  
  mz = malloc_err(sizeof(int)*3,err);
  forwardError(*err,__LINE__,NULL);
  
  mz[0] = mT*TEB[0];
  mz[1] = mT*TEB[1];
  mz[2] = mT*TEB[2];

  SC = alloc_SC(1,q,m,mz, &comp_totcalPP_update, &comp_totcalP_free,err);
  forwardError(*err,__LINE__,NULL);
  SC_ismul(SC);
  return SC;

}

/* OLD */
double smica_crit_classic(void *vsmic,error **err) {
  double res;
  int iq,m,nq;
  Smica *smic;


  smic = vsmic;
  m = smic->m;
  nq = smic->nq;
  if (smic->crit_classic_init==0) {
    for(iq = 0; iq < nq;iq++) { //precompute chol decomposed rq_hat
      int mx,my,info;
      double *rql;
      char uplo;
      //_DEBUGHERE_("","");
      
      rql = smic->rq_hat + iq*m*m;
      // chol
      uplo = 'L';
      dpotrf(&uplo,&m,rql,&m,&info);
      testErrorRetVA(info!=0,lowly_chol,"Could not cholesky decompose rq_hat using dpotrf (%d)",*err,__LINE__,0,info);
      //_DEBUGHERE_("","");
      
      // fill the U part with 0 (beware, f90!)
      for(mx=0;mx<m;mx++) {
        for(my=mx+1;my<m;my++) {
          rql[my*m+mx] = 0;
        }
      }
      //_DEBUGHERE_("","");  
    }
    smic->crit_classic_init=1;
  }

  res = 0;
  for(iq=0;iq<nq;iq++) {
    double kdd;
    //_DEBUGHERE_("iq %d -> %g",iq,res);
    memcpy(smic->z_buf,smic->rq_hat+m*m*iq,m*m*sizeof(double));
    kdd = kld(m,smic->z_buf,smic->rq+m*m*iq,err);
    //_DEBUGHERE_("%g",kdd);
    res += smic->wq[iq] * kdd;
    //_DEBUGHERE_(" -> %g (%g)",res,smic->wq[iq]);
    forwardError(*err,__LINE__,0);
  }
  //_DEBUGHERE_(" --> %g (%g)",res);  
  return -res; 
}

double kld(int n, double* rq_hat, double* rq, error **err) {
  char uplo,trans,diag,side;
  int nn, info,i;
  double *z;
  double done;
  double res;
  
  // suppose que rq_hat et la chol rq_hat
  // rq_hat et rq sont detruits en sortie
  
  //_DEBUGHERE_("","");
  //printMat(rq_hat,n,n);
  //_DEBUGHERE_("","");
  //printMat(rq,n,n);
         
  // chol rq
  uplo = 'L';
  nn = n;
  dpotrf(&uplo,&nn,rq,&nn,&info);
  testErrorRetVA(info!=0,lowly_chol,"Could not cholesky decompose rq using dpotrf (%d)",*err,__LINE__,0,info);
  //printMat(rq,nn,nn);
  
  // solve rq z = rq_hat
  side = 'L';
  trans = 'N';
  diag = 'N';
  done = 1;
  dtrsm(&side, &uplo, &trans, &diag, &nn, &nn, &done, rq, &nn, rq_hat, &nn);
  
  z = rq_hat;
  //printMat(z,n,n);
  
  // calcule sum(z^2)
  res = 0;
  //_DEBUGHERE_(": %g",res);
  for(i=0;i<n*n;i++) {
    res += z[i]*z[i];
    //_DEBUGHERE_(": %g",res);
  }
  //_DEBUGHERE_(": %g",res);
  for(i=0;i<n;i++) {
    res -= 2 * log(z[i*n+i]);
    //_DEBUGHERE_(": %g",res);
  }
  //_DEBUGHERE_(": %g",res);
  res -=n;
  //_DEBUGHERE_(": %g",res);
  return 0.5*res;  
}

