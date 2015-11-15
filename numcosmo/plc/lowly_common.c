/*
 *  lowly_common.c
 *  lowly_project
 *
 *  Created by Karim Benabed on 16/04/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */


#ifdef __PLANCK__
#include "HL2_likely/target/lowly_common.h"
#else
#include "lowly_common.h"
#endif


int* lowly_get_ell(int *pnell, int *inell, int lmax,error **err) {
  int *outell;
  int nell,i;
  
  nell = *pnell;
  
  if (nell==0) {
    outell = malloc_err(sizeof(int)*(lmax+1),err);
    forwardError(*err,__LINE__,NULL);
    for(i=0;i<lmax+1;i++) {
      outell[i]=i;
    }
    nell = lmax+1;
  } else {
    outell = malloc_err(sizeof(int)*nell,err);
    forwardError(*err,__LINE__,NULL);
    memcpy(outell,inell,sizeof(int)*nell);
  }
  *pnell = nell;
  if (lmax!=0) {
    for(i=0;i<nell;i++) {
      testErrorRetVA(outell[i]>lmax,lowly_lbig,"ell[%d]=%d>lmax (%d)",*err,__LINE__,NULL,i,outell[i],lmax);
    }
  }
  return outell;
}



int lowly_get_offset_cl(int *has_cl, int* offset_cl,int nell) {
  int offset, i;
  offset = 0;

  for(i=0;i<6;i++) {
    if (has_cl==NULL || has_cl[i]==1) {
      offset_cl[i] = offset;
      offset += nell;
    } else {
      offset_cl[i] = -1;
    }
  }
  return offset;
}

#ifndef NOHEALPIX
double * lowly_get_posvec(const long nside,
                    const long * pixel_indices,
                    const long npix_seen, 
                    const int ordering,error **err) {
  double *posvec;
  int i;
  
  if (npix_seen==0) {
    return NULL;
  }
  posvec = (double*) malloc_err(npix_seen*3*sizeof(double),err);
  forwardError(*err,__LINE__,NULL);
  
  if (ordering==HP_RING) { 
    // match
    // define position vectors for all pixels
    for (i=0;i<npix_seen;i++) {
	    pix2vec_ring(nside,pixel_indices[i],&posvec[3*i]);
    }
  } else {
    // nested
    for (i=0;i<npix_seen;i++) {
	    pix2vec_nest(nside,pixel_indices[i],&posvec[3*i]);
    }
  }
  return posvec;
}

 


int lowly_which_order(char *ordering,error **err) {
  char unknownOrder[1024];
  if (strcasecmp(ordering,"ring")==0) { 
    return HP_RING;
  }
  if (strcasecmp(ordering,"nested")==0) { 
    return HP_NEST;
  }
  sprintf(unknownOrder,"Don't know ordering '%s'",ordering);
  testErrorRet(0,lowly_unkorder,unknownOrder,*err,__LINE__,0);
}

 


double lowly_computeNl(double noisevar, long nside) {
  return noisevar*4*PI/(12.*nside*nside);
}

 



long* lowly_build_pixel_list(unsigned char *Mask,long *pnpix, error **err) {
  long ipix, npix, jpix;
  long  *pixel_indices;
  
  npix = *pnpix;
  // Counting the number of pixels seen and Building the vector of indices of seen pixels
  pixel_indices = (long*) malloc_err(npix*sizeof(long),err);
  forwardError(*err,__LINE__,NULL);
  jpix=0;
  for (ipix =0; ipix < npix; ipix++) { 
    if (Mask[ipix] == PIXEL_SEEN) {
	    pixel_indices[jpix]=ipix;
      jpix++;
    }
  }
  testErrorRet(jpix==0,lowly_unkorder,"Found no pixel !!!",*err,__LINE__,NULL);
  
  *pnpix = jpix;
  return pixel_indices;
}






int lowly_nmode(int nell,int *ell) {
  int l,nmode;

  nmode =0;
  // count the number of modes
  for(l=0;l<nell;l++) {
    nmode += 2*ell[l]+1;
  }
  return nmode;
}



void * lowly_ring_reorder(void* from,int nside, size_t sz,error **err) {
  int i,npix;
  void* dest;
  
  if (from==NULL) {
    return NULL;
  }
  npix=12*nside*nside;
  dest = malloc_err(sz*npix,err);
  forwardError(*err,__LINE__,NULL);
  for(i=0;i<npix;i++) {
    long ip;
    nest2ring(nside,i,&ip);
    memcpy(dest + ip*sz,from + i*sz,sz);
  }
  return dest;
}



void * lowly_ring_reorder2D(void* from,int nside, size_t sz,error **err) {
  int i,j,npix;
  void* dest;
  
  if (from==NULL) {
    return NULL;
  }
  npix=12*nside*nside;
  dest = malloc_err(sz*_SZT_(npix)*_SZT_(npix),err);
  forwardError(*err,__LINE__,NULL);
  for(i=0;i<npix;i++) {
    long ip;
    nest2ring(nside,i,&ip);
    for(j=0;j<npix;j++) {
      long jp;
      nest2ring(nside,j,&jp);      
      memcpy(dest + (ip*npix+jp)*sz,from + (i*npix+j)*sz,sz);
    }
  }
  return dest;
}

double lowly_XtRX_lkl(double *S,double *X, double* buffer,long npix_tot,
                    error **err) {
  double res;
  
  res = lowly_XtRX_lkl_K(S,X,buffer,npix_tot,NULL,0,err);
  forwardError(*err,__LINE__,-1);
  return res;
}

double lowly_XtRX_lkl_K(double *S,double *X, double* buffer,long npix_tot,double *Ktilde,long npix_other,
                    error **err) {
  int npix,info,xinc,npixo;
  char uplo,trans,diag,side;
  double quad,logdet,regul;
  int i;
  int tot_time;
  double fone;
  TIMER_DEFS;
  
  tot_time=0;
  
  npix = npix_tot;
  uplo = 'L';
  
  TIMER_IN;
  dpotrf(&uplo,&npix,S,&npix,&info);
  testErrorRetVA(info!=0,lowly_chol,"Could not cholesky decompose using dpotrf (%d)",*err,__LINE__,0,info);
  TIMER_OUT;
  tot_time+=TIMER_MSEC;
  //fprintf(stderr,"chol : %d msec\n",TIMER_MSEC);
  
  
  TIMER_IN;
  memcpy(buffer,X,npix_tot*sizeof(double));
  // Now solve covmat^1/2.y=X_seen, and place the result in X_seen
  trans='N'; diag='N';xinc=1;
  dtrsv(&uplo,&trans,&diag,&npix,S,&npix,buffer,&xinc);
  TIMER_OUT;
  tot_time+=TIMER_MSEC;
  //fprintf(stderr,"solve : %d msec\n",TIMER_MSEC);
  
  TIMER_IN;
  // compute quad and loget terms
  quad = 0;
  logdet = 0;
  for (i=0;i<npix;i++) {
    double ab,re;
    ab = buffer[i];
    re=S[i*npix_tot+i];
#ifdef LOWLY_PARANOID_NAN_TEST
    testErrorRetVA((isnan(ab)!=0 || isinf(ab)!=0),pmc_infinite,"a_bar is invalid (%g) at indice %d",*err,__LINE__,0,ab,i);
    testErrorRetVA((isnan(re)!=0 || isinf(re)!=0),pmc_infinite,"Re is invalid (%g) at indice %d",*err,__LINE__,0,re,i);
#endif
    logdet += 2 * log(re);
    quad += ab * ab;
  }
#ifdef LOWLY_PARANOID_NAN_TEST
  testErrorRetVA((isnan(quad)!=0 || isinf(quad)!=0 || quad<=0),pmc_infinite,"quad term is invalid (%g)",*err,__LINE__,0,quad);
  testErrorRetVA((isnan(logdet)!=0 || isinf(logdet)!=0),pmc_infinite,"logdet term is invalid (%g)",*err,__LINE__,0,logdet);
#endif
  TIMER_OUT;
  //fprintf(stderr,"final sum : %d msec\n",TIMER_MSEC);
  tot_time+=TIMER_MSEC;
  //fprintf(stderr,"total_time = %d msec\n",tot_time);

  regul = 0;
  if (Ktilde!=NULL) {
    memcpy(buffer,Ktilde,npix_tot*npix_other*sizeof(double));
    side = 'L';
    uplo = 'L';
    trans = 'N';
    diag = 'N';
    npixo = npix_other;
    fone = 1;
    dtrsm(&side, &uplo, &trans, &diag, &npix, &npixo, &fone, S, &npix, buffer, &npixo);
    for(i=0;i<npixo;i++) {
      int j;
      regul += ddot(&npixo,(buffer+i),&xinc,(buffer+i),&xinc);  
    }
  }
  
  return - (quad + logdet + regul) / 2.;  
}




void lowly_print_stat(FILE* where,int doprint, int n_tests,int time_build, int time_chol,int time_tot) {
  if ((doprint==1) && n_tests>0) {
    fprintf(where,"Timers after %d computations\n",n_tests);
    fprintf(where,"-> mean build time %g msec\n",time_build*1./n_tests);
    fprintf(where,"-> mean chol+solve time %g msec\n",time_chol*1./n_tests);
    fprintf(where,"-> mean total time %g msec\n",time_tot*1./n_tests);
  }
}


void lowly_dump(char *filename, void* ptr, size_t sz, error **err) {
  FILE *f;
  
  f = fopen_err(filename,"w",err);
  forwardError(*err,__LINE__,);
  fwrite(ptr,sz,1,f);
  fclose(f);
  return;
}
#endif