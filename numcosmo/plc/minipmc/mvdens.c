/*
 *  mvdens.c
 *  likely
 *
 *  Created by Karim Benabed on 10/03/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifdef __PLANCK__
#include "HL2_likely/tools/mvdens.h"
#else
#include "mvdens.h"
#endif

mvdens * mvdens_alloc(size_t size, error **err) {
  mvdens *g;
  size_t ndim;

  ndim = size;
  g = (mvdens *) malloc_err(sizeof(mvdens),err);
  forwardError(*err,__LINE__,NULL);
  g->buf = calloc_err(ndim*(ndim+1), sizeof(double),err);
  forwardError(*err,__LINE__,NULL);
  
  g = mvdens_init(g,size,g->buf,((double*) g->buf)+ndim,err);
  forwardError(*err,__LINE__,NULL);
  g->own_buf=1;

  return g;
}

mvdens *mvdens_init(mvdens *g, size_t size, double* mn, double *std, error **err) {

  testErrorRet(size<=0, mv_negative, "Invalid size = %d, has to be positive", *err, __LINE__, NULL);

  g->ndim = size;
  g->mean = mn;
  
  /* This is used to store the square root (Cholesky) of the cov.   */
  g->std = std;
  
  gsl_set_error_handler_off();
  g->mean_view_container = gsl_vector_view_array(g->mean,g->ndim);
  g->mean_view = &(g->mean_view_container.vector);
  g->std_view_container = gsl_matrix_view_array(g->std,g->ndim,g->ndim); 
  g->std_view   = &(g->std_view_container.matrix);
  g->band_limit = size;
  g->chol       = 0;
  g->df         = -1;
  g->own_buf    = 0;
  g->detL       = -1.0;

  g->x_tmp = malloc_err(sizeof(double)*g->ndim,err);
  forwardError(*err,__LINE__,NULL);
  g->x_tmp_view_container = gsl_vector_view_array(g->x_tmp,g->ndim);
  g->x_tmp_view = &(g->x_tmp_view_container.vector);
  
  return g;
}

void mvdens_from_meanvar(mvdens *m, const double *pmean, const double *pvar, double correct) {
  int j, k;
  for (j=0; j<m->ndim; j++) {
    for (k=0; k<m->ndim; k++) {
      m->std[j*m->ndim+k] = pvar[j*m->ndim+k]*correct;
    }
    if (pmean!=NULL) 
      m->mean[j] = pmean[j];
  }
}

void mvdens_empty(mvdens *g) {
  if (g->own_buf==1) {
    free(g->buf);
  }
}

void mvdens_free(mvdens ** g) {
  mvdens_empty(*g);
  free((*g)->x_tmp);
  free(*g);
  *g = NULL;
}

void mvdens_free_void(void **g)
{
   mvdens_free((mvdens**)g);
}

void mvdens_print(FILE* where, mvdens* what) {
  FILE *rhere;
  size_t j,k,ndim;
  
  rhere=where;
  if (where==NULL) {
    rhere=stdout;
  }
  ndim=what->ndim;
  if(what->df!=-1)
    fprintf(rhere, "df = %d\n", what->df);
  if(what->chol==1)
    fprintf(rhere, "Beware! Cholesky decomposed std, df=%d\n", what->df);
  for (j=0;j<ndim;j++) {
    fprintf(rhere,"% .3g | ",what->mean[j]);
    for (k=0;k<j;k++) 
      fprintf(rhere,"% .3g ",what->std[j*ndim+k]);
    fprintf(rhere,"% .3g\n",what->std[j*ndim+j]);
    //fprintf(rhere,"\n");
  }
  return;
}

void mvdens_dump(FILE* where, mvdens* what) {
  FILE *rhere;
  size_t j,k,ndim;
  
  rhere=where;
  if (where==NULL) {
    rhere=stdout;
  }
  ndim=what->ndim;
  fprintf(rhere, "%zu %d %zu %d\n", ndim, what->df, what->band_limit, what->chol);
  for (j=0;j<ndim;j++) {
    fprintf(rhere,"% 20.10g ",what->mean[j]);
  }
  fprintf(rhere,"\n");
  for (j=0;j<ndim;j++) {
    for (k=0;k<ndim;k++) 
      fprintf(rhere,"% 20.15g ",what->std[j*ndim+k]);
    fprintf(rhere,"\n");
  }
  return;
}

void mvdens_chdump(const char *name, mvdens* what, error **err)
{
   FILE *F;

   F = fopen_err(name, "w", err);
   forwardError(*err, __LINE__,);
   mvdens_dump(F, what);
   fclose(F);
}

mvdens* mvdens_dwnp(FILE* where, error **err) {
  mvdens *what;
  FILE *rhere;
  size_t j,k,ndim,bdl;
  int df, chol;
  
  rhere=where;
  if (where==NULL) {
    rhere=stdin;
  }
  testErrorRet(fscanf(rhere, "%zu %d %zu %d\n", &ndim, &df, &bdl, &chol)!=4,
	       mv_io, "Error while reading mvdens file: Four integers \n"
	       "(ndim, dof, bandlimit, cholesky flag) expected on first line",
	       *err, __LINE__, NULL);
  
  what=mvdens_alloc(ndim,err);
  forwardError(*err,__LINE__,NULL);
  
  for (j=0; j<ndim; j++) {
    testErrorRetVA(fscanf(rhere, "%lg", what->mean+j)!=1, mv_io,
		   "Cannot read mvdens (mean element #%d) double expected",
		   *err, __LINE__, NULL, j);
  }
  for (j=0; j<ndim; j++) {
    for (k=0; k<ndim; k++) {
      testErrorRetVA(fscanf(rhere, "%lg", what->std+j*ndim+k)!=1, mv_io,
		     "Cannot read mvdens (covariance element #(%d,%d)) double expected",
		     *err, __LINE__, NULL, j, k);
     }
  }
  
  what->band_limit = bdl;
  what->df         = df;
  what->chol       = chol;
  if (what->chol==1) {
     what->detL = determinant(what->std, what->ndim);
  }

  return what;
}

mvdens *mvdens_read_and_chol(FILE *where, error **err)
{
   mvdens *what;

   what = mvdens_dwnp(where, err);
   forwardError(*err, __LINE__, NULL);
   if (what->chol==0) {
      mvdens_cholesky_decomp(what, err);
      forwardError(*err, __LINE__, NULL);
   }

   //fprintf(stderr, "mvdens_read_and_chol: det = %g\n", what->detL);

   return what;
}

/* L has to be a trianglular matrix (e.g. after Cholesky-decomposition) */
double determinant(const double *L, size_t ndim)
{
   size_t i;
   double det;

   for (i=0,det=1.0; i<ndim*ndim; i+=(ndim+1))
     det *= L[i];

   return det;
}

void mvdens_cholesky_decomp(mvdens* self, error **err) {
  gsl_set_error_handler_off();
  if (self->chol==1)
    return;

#ifdef __MVDENS_PARANOID_DEBUG__
  double * _std;
  _std=malloc_err(self->ndim*self->ndim*sizeof(double),err);
  forwardError(*err,__LINE__,);
  memcpy(_std,self->std,sizeof(double)*self->ndim*self->ndim);
#endif

/*#ifdef _WITH_LAPACK_
  char uplo;
  int info; 
  int n;
  uplo = 'L';
  info = 0; 
  n = self->ndim;
  _DEBUGHERE_("%c %d %d %d %d",uplo,n,info,sizeof(char),sizeof(int))
  _DEBUGHERE_("%p %p %p",&uplo,&n,&info)
  
  dpotrf( &uplo, &n, self->std, &n, &info );
  
  if (info!=0) {
#else*/
  if (gsl_linalg_cholesky_decomp(self->std_view) ==  GSL_EDOM) {
//#endif

#ifdef __MVDENS_PARANOID_DEBUG__
    memcpy(self->std, _std, sizeof(double)*self->ndim*self->ndim);
    free(_std);
#endif
/*#ifdef _WITH_LAPACK_
    testErrorRetVA(info!=0,mv_cholesky,"Could not cholesky decompose using dpotrf (%d)",*err,__LINE__,,info);
#else*/
    *err = addError(mv_cholesky,"Cholesky decomposition failed",*err,__LINE__);
//#endif
    return;
  }

  self->detL = determinant(self->std, self->ndim);

#ifdef __MVDENS_PARANOID_DEBUG__
  free(_std);
#endif

  self->chol=1;
}  

double* mvdens_ran(double* dest, mvdens * g, gsl_rng * r,error **err) {
  size_t i;
  double *res;
  gsl_vector * res_view;
  gsl_vector_view res_view_container;
  double val,corr,u;
  
  mvdens_cholesky_decomp(g,err);
  forwardError(*err,__LINE__,NULL);
  
  MALLOC_IF_NEEDED(res,dest,sizeof(double)*g->ndim,err);
  forwardError(*err,__LINE__,NULL);

  gsl_set_error_handler_off();
  /* Generate a N(O,I) distributed vector               */
  for (i = 0; i < g->ndim; i++) {
    val=gsl_ran_gaussian(r, 1.0);
    res[i]=val;
  }
  
  /* Make it N(O,Sigma)                         */
  res_view_container=gsl_vector_view_array(res,g->ndim);
  res_view=&(res_view_container.vector);
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, g->std_view, res_view);

  corr=1;
  if(g->df!=-1) {
    // make it student-t
    u = gsl_ran_chisq(r, g->df);
    corr=sqrt(g->df/u);
  }
 
  // add mean and correct if student-t
  for(i=0;i<g->ndim;i++) {
    res[i]*=corr;
    res[i]+=g->mean[i];
  }
  return res;
}

double scalar_product(mvdens *g, const double *x, error **err) {
  gsl_vector * x_tmp_view;
  double *x_tmp;
  double logl;
  size_t i, ndim;
  
  gsl_set_error_handler_off();
  
  mvdens_cholesky_decomp(g,err);
  forwardError(*err,__LINE__,-1);
  
  ndim = g->ndim;
  /* Copy vector */
  x_tmp = g->x_tmp;
  x_tmp_view = g->x_tmp_view;
  
  memcpy((void*)x_tmp,(void*)x,(size_t)sizeof(double)*ndim);
  
  /* Subtract mean */
  gsl_vector_sub(x_tmp_view, g->mean_view);
  

  /* Compute Sigma^{-1/2}(x-mu) (by substitution) */
  gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, g->std_view, x_tmp_view);

  /* Compute squared norm */
  for (i=0,logl=0.0; i<ndim; i++) {
    logl += x_tmp[i]*x_tmp[i];
  }
  
  return logl;
}

/* ============================================================ *
 * Returns the negative log-likelihood for a mv Gaussian pdf.   *
 * ============================================================ */
double mvdens_log_pdf(mvdens *g, const double *x, error **err) {
  double log_g, logdet;
  size_t ndim, n, p;
  double gamma, norm, r;

  gsl_set_error_handler_off();
  
  mvdens_cholesky_decomp(g, err);
  forwardError(*err,__LINE__,-1);
  
  ndim   = g->ndim;
  logdet = log(g->detL);
  
  /* log_g contains the scalar product xt . var . x */
  log_g = scalar_product(g, x, err);
  forwardError(*err,__LINE__,-1);  
  /* Student-t distribution with df degrees of freedom */
  if (g->df!=-1) {
    n = g->df;
    p = g->ndim;
    
    gamma = gsl_sf_lngamma(0.5*(n+p));
    norm  = gsl_sf_lngamma(0.5*n) + 0.5*p*log(n*PI);
    r     = 0.5*(n+p)*log(1.0+(1.0/n)*log_g);
    return gamma - (norm + logdet + r);
  } else {
     /* logdet = log|Sigma| = 0.5*log|Covariance| ! */
     r = -0.5*(ndim*ln2pi + log_g) - logdet;
     return r;
  }
}

double mvdens_log_pdf_void(void *g, const double *x, error **err) {
   double res;
   res = mvdens_log_pdf((mvdens*)g, x, err);
   forwardError(*err, __LINE__, 0.0);
   return res;
}

/* ============================================================ *
 * Inverts the covariance of m and returns its determinant.     *
 * ============================================================ */
double mvdens_inverse(mvdens *m, error **err)
{
   double det;

   mvdens_cholesky_decomp(m, err);
   forwardError(*err, __LINE__, 0.0);

   /* Note: gsl_linalg_cholesky_invert is defined in gsl version 1.14 and higher */
   gsl_linalg_cholesky_invert(m->std_view);
   // Now, m->std contains the inverse of the original m->std before cholesky decomposition
   m->chol = 0;
   det = determinant(m->std, m->ndim);

   return det;
}

mix_mvdens *mix_mvdens_alloc(size_t ncomp, size_t size,error **err) {
  mix_mvdens * m; 
  size_t i,dcomp,std0,ndim;
  void * data;
  
  m = (mix_mvdens *) malloc_err(sizeof(mix_mvdens),err);
  forwardError(*err,__LINE__,NULL);
  
  m->ncomp = ncomp;
  m->ndim  = size;
  m->wght  = calloc_err(ncomp,sizeof(double),err);
  forwardError(*err,__LINE__,NULL);

  gsl_set_error_handler_off();
  
  m->wght_view_container=gsl_vector_view_array(m->wght,ncomp);
  m->wght_view=&(m->wght_view_container.vector);
  
  m->cwght = calloc_err(ncomp,sizeof(double),err);
  forwardError(*err,__LINE__,NULL);
  m->cwght_view_container=gsl_vector_view_array(m->cwght,ncomp);
  m->cwght_view=&(m->cwght_view_container.vector);

  m->comp = (mvdens **) malloc_err(ncomp * sizeof(mvdens *),err);
  forwardError(*err,__LINE__,NULL);

  ndim = m->ndim;
  m->buf_comp = malloc_err(ncomp*(sizeof(mvdens)+ndim*(ndim+1)*sizeof(double)),err);
  forwardError(*err,__LINE__,NULL);
  memset(((char*)m->buf_comp)+ncomp*sizeof(mvdens),0,ncomp*(ndim*(ndim+1)*sizeof(double)));
  data  = ((char*)m->buf_comp)+ncomp * sizeof(mvdens);
  dcomp = ndim*(ndim+1)*sizeof(double);
  std0  = ndim*sizeof(double);

  for (i=0; i<ncomp; i++) {
    m->comp[i] = mvdens_init((mvdens*)((char*)(m->buf_comp)+i*sizeof(mvdens)), size,
			     (double*)((char*)data+dcomp*i), (double*)((char*)data+dcomp*i+std0), err);
    forwardError(*err,__LINE__,NULL);
  }

  m->init_cwght = 0;

  return m;
}

/* ============================================================ *
 * Copies the content of 'source' to 'target'. Memory for the   *
 * has to be allocated before, e.g. with mix_mvdens_alloc.      *
 * ============================================================ */
void mix_mvdens_copy(mix_mvdens *target, const mix_mvdens *source, error **err)
{
   target->ncomp      = source->ncomp;
   target->ndim       = source->ndim;
   target->init_cwght = source->init_cwght;

   memcpy((void*)target->wght, (void*)source->wght, sizeof(double)*source->ncomp);
   memcpy((void*)target->cwght, (void*)source->cwght, sizeof(double)*source->ncomp);

   memcpy((void*)target->comp, (void*)source->comp, sizeof(mvdens*) * source->ncomp);
}

void mix_mvdens_free_void(void **m)
{
   mix_mvdens_free((mix_mvdens**)m);
}

void mix_mvdens_free(mix_mvdens **m) {
  size_t i;
  
  for (i = 0; i < (*m)->ncomp; i++) {
    mvdens_empty((*m)->comp[i]);
  }
  free((*m)->comp);
  free((*m)->buf_comp);
  free((*m)->wght);
  free((*m)->cwght);
  free(*m);
  *m=NULL;
}

void mix_mvdens_print(FILE* where,mix_mvdens* what) {
  FILE *rhere;
  int i, ncomp;
  
  rhere=where;
  if (where==NULL) {
    rhere=stdout;
  }
  
  ncomp=what->ncomp;
  
  for (i=0; i<ncomp; i++) {
    fprintf(rhere, "*************************************\n");
    fprintf(rhere, "Weight of component #%d is %g\n", i, what->wght[i]);
    if (what->wght[i]!=0) {
      mvdens_print(rhere, what->comp[i]);
    } 
  }
}

void mix_mvdens_dump(FILE* where,mix_mvdens* what) {
  FILE *rhere;
  int i,ncomp;

  rhere=where;
  if (where==NULL) {
    rhere=stdout;
  }

  ncomp=what->ncomp;
  fprintf(rhere, "%d %d\n", (int)what->ncomp, (int)what->ndim);
  for (i=0;i<ncomp;i++) {
    fprintf(rhere,"%.8g\n",what->wght[i]);
    mvdens_dump(rhere,what->comp[i]);
  }
}

mix_mvdens* mix_mvdens_dwnp(FILE* where, error **err) {
  mix_mvdens *what;
  FILE *rhere;
  size_t i,j,k,ndim,ncomp;
  size_t df,bdl,chol;
  mvdens *w;
  
  rhere=where;
  if (where==NULL) {
    rhere=stdin;
  }
  testErrorRet(fscanf(rhere,"%zu %zu\n",&ncomp,&ndim)!=2,mv_io,
	       "Cannot read mix_mvdens first line (two integers expected)",*err,__LINE__,NULL);

  //fprintf(stderr,"mix_mvdens_dwnp -> %ld %ld\n",ncomp,ndim);
  what=mix_mvdens_alloc(ncomp,ndim,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<ncomp;i++) {
    testErrorRetVA(fscanf(rhere,"%lg\n",&(what->wght[i]))!=1, mv_io,
		   "Cannot read mix_mvdens (error for weight of component %d)", *err, __LINE__, NULL, i);
    testErrorRetVA(fscanf(rhere,"%zu %zu %zu %zu\n",&ndim,&df,&bdl, &chol)!=4,mv_io,
		   "Cannot read mix_mvdens (error on first line for component %d)",*err,__LINE__,NULL,i);
    //fprintf(stderr,"mix_mvdens_dwnp(%d) -> %ld %ld %ld %ld\n",i,ndim,df,bdl, chol);
    testErrorRetVA(ndim!=what->ndim,mv_io,
		   "Cannot read mix_mvdens, dimension ndim=%d of component %d invalid (expected %d)",
		   *err, __LINE__, NULL, what->ndim, i, ndim);
    
    w=what->comp[i];
    
    for (j=0; j<ndim; j++) {
      testErrorRetVA(fscanf(rhere, "%lg", w->mean+j)!=1,mv_io,
		     "Cannot read mix_mvdens component %d, (error for mean entry #%d)", *err,__LINE__, NULL, i, j);
      //fprintf(stderr,"mix_mvdens_dwnp(%d) %d->%g\n",i,j,w->mean[j]);
    }
    //fprintf(stderr,"mix_mvdens_dwnp(%d) std\n",i);
    for (j=0; j<ndim; j++) {
      for (k=0; k<ndim; k++) {
        testErrorRetVA(fscanf(rhere, "%lg", w->std+j*ndim+k)!=1,mv_io,
		       "Cannot read mix_mvdens component %d (error for variance entry (%d,%d)",
		       *err,__LINE__,NULL,i,j,k)      //fprintf(stderr,"%g ",w->std+j*ndim+k);
      }
      //fprintf(stderr,"\n");
    }
    
    w->band_limit=bdl;
    w->df=df;
    w->chol=chol;
    if (w->chol==1) 
      w->detL = determinant(w->std, w->ndim);
  }

  return what;
}

size_t mix_mvdens_size(size_t ncomp,size_t ndim) {
  size_t rsz;

  rsz  = 0;
  rsz += sizeof(size_t)*2;                     /* ncomp, ndim */
  rsz += ncomp*sizeof(double);                 /* wght */
  rsz += ncomp*sizeof(int);                    /* df */
  rsz += ncomp*sizeof(size_t);                 /* band_limit */
  rsz += ncomp*ndim*(1+ndim)*sizeof(double);   /* mean, std */

  return rsz;  
}

#define add_size(me,type,n) me=((char*) me)+(sizeof(type)*n)
#define store(me,type,what) (*(type*)me=what); add_size(me,type,1);
#define unstore(where,me,type) where=*(type*)me; add_size(me,type,1);


void* serialize_mix_mvdens(mix_mvdens* self, size_t *sz,error** err)
{
  void *serialized, *crt;
  size_t rsz,i;
  
  //fprintf(stderr,"seri : 1\n");
  
  rsz = mix_mvdens_size(self->ncomp,self->ndim);
  serialized = malloc_err(rsz,err);
  forwardError(*err,__LINE__,NULL);
  
  //fprintf(stderr,"seri : 2\n");
  
  crt=serialized;
  //fprintf(stderr,"%p %p\n",crt,serialized);
  store(crt,size_t,self->ncomp);

  //fprintf(stderr,"%p %p\n",crt,serialized);
  store(crt,size_t,self->ndim);

  //fprintf(stderr,"%p %p\n",crt,serialized);
  memcpy(crt,self->wght,sizeof(double)*self->ncomp);
  add_size(crt,double,self->ncomp);
  //fprintf(stderr,"%p %p\n",crt,serialized);
  for(i=0;i<self->ncomp;i++) {
    //fprintf(stderr,"seri : 2 %d a\n",i);

    store(crt,int,self->comp[i]->df);
    //fprintf(stderr,"%p %p\n",crt,serialized);
    //fprintf(stderr,"seri : 2 %d b\n",i);

    store(crt,size_t,self->comp[i]->band_limit);
    //fprintf(stderr,"%p %p\n",crt,serialized);
    //fprintf(stderr,"seri : 2 %d c\n",i);

    mvdens_cholesky_decomp(self->comp[i],err);
    forwardError(*err,__LINE__,NULL);
    //fprintf(stderr,"seri : 2 %d d\n",i);

    memcpy(crt,self->comp[i]->mean,self->ndim*sizeof(double));
    add_size(crt,double,self->ndim);
    //fprintf(stderr,"%p %p\n",crt,serialized);
    //fprintf(stderr,"seri : 2 %d e\n",i);

    memcpy(crt,self->comp[i]->std,self->ndim*self->ndim*sizeof(double));
    add_size(crt,double,self->ndim*self->ndim);
    //fprintf(stderr,"%p %p\n",crt,serialized);
    //fprintf(stderr,"seri : 2 %d e\n",i);

  }
  //fprintf(stderr,"%p %p\n",crt,serialized);
  testErrorRetVA((char*)crt-(char*)serialized!=rsz, mv_serialize,
		 "Bad size for serialization (expected %d, got %d)",
		 *err, __LINE__, NULL, rsz, (char*)crt-(char*)serialized);
  //fprintf(stderr,"seri : 3\n");
  *sz=rsz;
  return serialized;
}

mix_mvdens* deserialize_mix_mvdens(void* serialized, size_t sz, error **err) {
  size_t ndim,ncomp,rsz,i;
  void *crt;
  mix_mvdens* m;
  
  //fprintf(stderr,"deseri : 1\n");
  
  crt=serialized;
  //fprintf(stderr,"%p %p\n",crt,serialized);
  
  unstore(ncomp,crt,size_t);
  //fprintf(stderr,"%p %p\n",crt,serialized);
  unstore(ndim,crt,size_t);
  //fprintf(stderr,"%p %p\n",crt,serialized);
  rsz = mix_mvdens_size(ncomp,ndim);
  
  testErrorRetVA(rsz!=sz,mv_serialize,"Cannot deserialize, bad size (Expected %d got %d)",*err,__LINE__,NULL,sz,rsz);

  //fprintf(stderr,"deseri : 2\n");
  m = mix_mvdens_alloc(ncomp,ndim,err);
  forwardError(*err,__LINE__,NULL);
  
  memcpy(m->wght,crt,sizeof(double)*ncomp);
  add_size(crt, double,ncomp);
  //fprintf(stderr,"%p %p\n",crt,serialized);
  
  for(i=0;i<ncomp;i++) {
    //fprintf(stderr,"deseri : 3 %d\n",i);

    unstore(m->comp[i]->df,crt,int);
    //fprintf(stderr,"%p %p\n",crt,serialized);

    unstore(m->comp[i]->band_limit,crt,size_t);
    //fprintf(stderr,"%p %p\n",crt,serialized);

    m->comp[i]->chol = 0;
    memcpy(m->comp[i]->mean,crt,sizeof(double)*ndim);
    add_size(crt,double,m->ndim);
    //fprintf(stderr,"%p %p\n",crt,serialized);

    memcpy(m->comp[i]->std,crt,sizeof(double)*ndim*ndim);
    add_size(crt,double,m->ndim*m->ndim);
    //fprintf(stderr,"%p %p\n",crt,serialized);

    m->comp[i]->chol = 1;
    /* comp was Cholesky-decomposed in serialize, therefore we can call determinant */
    m->comp[i]->detL = determinant(m->comp[i]->std, ndim);
  }
  //fprintf(stderr,"deseri : 4\n");
  //fprintf(stderr,"%p %p\n",crt,serialized);
  
  testErrorRetVA((char*)crt-(char*)serialized!=rsz, mv_serialize,
		 "Cannot deserialize, bad size (Expected %d, got %d)",
		 *err, __LINE__, NULL, rsz, (char*)crt-(char*)serialized);
    
  //fprintf(stderr,"deseri : 5\n");

  return m;
}

void mix_mvdens_cholesky_decomp(mix_mvdens* self, error **err) {
  size_t i;
  for (i=0;i<self->ncomp;i++) {
    mvdens_cholesky_decomp(self->comp[i],err);
     forwardError(*err,__LINE__,);
  }
}

double* mix_mvdens_ran(double* dest,size_t *index, mix_mvdens *m, gsl_rng * r,error **err) {
  /* Draw random realization according to mixture */
  size_t i;
  double * x;
  gsl_set_error_handler_off();
  
  if (m->init_cwght == 0) { /* Initialize cumulative weights */
    
    double sum_wght=0.0;
    
    for (i=0;i<m->ncomp;i++) 
      sum_wght += m->wght[i];
    if (fabs(sum_wght-1.0) > 1e-6) {
      testErrorRet(sum_wght==0,mv_negWeight,"All weights to zero !! cannot sample",*err,__LINE__,NULL);
      gsl_vector_scale(m->wght_view,1.0/sum_wght);
    }
    m->cwght[0]=m->wght[0];
    for (i=1;i<m->ncomp;i++) {
      m->cwght[i]=m->cwght[i-1]+m->wght[i];
    }
    m->init_cwght=1;
  }
  
  /* Draw index according to weights */
  *index = ranindx(m->cwght,m->ncomp,r,err);
  forwardError(*err,__LINE__,NULL);
  //fprintf(stderr,"%d -> ",*index);

  /* Draw gaussian random variable of corresponding law,
   the last index corresponding to the initial proposal ? */
  x = mvdens_ran(dest, m->comp[*index], r, err);
  forwardError(*err,__LINE__,NULL);

  return x;
}

size_t ranindx(double *cw, size_t nE, gsl_rng *r, error **err) {
  
  /* Draw index according to probability weights */
  /* Uses cumulative weight vector cw */
  
  size_t index=0;
  double z;
  gsl_set_error_handler_off();
  
  z = gsl_ran_flat(r,0.0,1.0); /* Weights are assumed normalized */
  while (cw[index]< z && index < nE-1) {
    index++;
  }
  //  printf("index is %d \n",index);
  return index;
  
}

double mix_mvdens_log_pdf(mix_mvdens *m, const double *x, error **err) {
  /* Log pdf of MOG */
  
  double *logdens;
  double logdens_max_loc, logdens_tmp, lpos;
  int i,iloc;
  int ncomp = m->ncomp;
  
  //fprintf(stderr,"--->>> %d\n",ncomp);
  logdens = malloc_err(ncomp*sizeof(double),err);
  forwardError(*err,__LINE__,0);
  
  iloc=0;
  while (m->wght[iloc] <= 0)
    iloc++;
  testErrorRetVA(iloc>=ncomp,mv_outOfBound,"Indice out of bonds %d >= %d ",*err,__LINE__,0,iloc,ncomp);

  logdens_max_loc = mvdens_log_pdf(m->comp[iloc],x,err);
  forwardError(*err,__LINE__,0);
  logdens[iloc] = logdens_max_loc;
  
  for (i=iloc+1;i<ncomp;i++) {
    if (m->wght[i] > 0.0) {
      logdens_tmp = mvdens_log_pdf(m->comp[i],x,err);
      forwardError(*err,__LINE__,0);
      logdens[i] = logdens_tmp;
      if (logdens_max_loc < logdens_tmp) 
        logdens_max_loc=logdens_tmp;
    }
  }
  lpos=0.0;
  for (i=0;i<ncomp;i++) {
    if (m->wght[i] > 0.0) {
      logdens[i] = logdens[i]-logdens_max_loc+DYNMAX;
      lpos = lpos + exp(logdens[i])*m->wght[i];
    }
  }
  lpos = logdens_max_loc+log(lpos)-DYNMAX;
  free(logdens);

  return lpos;
}

double mix_mvdens_log_pdf_void(void *m, const double *x, error **err)
{
   double res;
   res = mix_mvdens_log_pdf((mix_mvdens*)m, x, err);
   forwardError(*err, __LINE__, 0);
   return res;
}

void mvdens_set_band_limit(mvdens *self, size_t value) {
  self->band_limit=value;
}

void mix_mvdens_set_band_limit(mix_mvdens *self, size_t value) {
  size_t i;
  for(i=0;i<self->ncomp;i++)
    self->comp[i]->band_limit=value;
}

#define EPS 1.0e-6
double effective_number_of_components(const mix_mvdens *self, error **err)
{
   double enc;
   int i;

   for (i=0,enc=0.0; i<self->ncomp; i++) {
      enc += self->wght[i]*self->wght[i];
   }
   testErrorRet(enc<EPS, mv_negWeight, "All mixture weights are zero", *err, __LINE__, 0.0);

   return 1.0/enc;
}
#undef EPS


#ifdef HAS_HDF5
void mvdens_hdfdump_infile(mvdens* mv,hid_t loc_id,error **err) {
  hid_t       group_id;   
  herr_t      status;
  hsize_t     dims[2];
  int i;
  char fname[200];
  
  H5Iget_name( loc_id, fname, 200 );
  
  // create group 
  group_id = H5Gcreate( loc_id, "mvdens", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  testErrorRetVA(group_id<0,hdf5_base,"cannot create group /mvdens in file %s (got %d)",*err,__LINE__,,fname,group_id);
  
  status = H5LTset_attribute_int( group_id, ".", "ndim", &(mv->ndim), 1);
  testErrorRetVA(status<0,hdf5_base,"cannot save /mvdens/ndim in file %s (got %d)",*err,__LINE__,,fname,status);
  
  status = H5LTset_attribute_int( group_id, ".", "df", &(mv->df), 1);
  testErrorRetVA(status<0,hdf5_base,"cannot save /mvdens/df in file %s (got %d)",*err,__LINE__,,fname,status);
  
  status = H5LTset_attribute_int( group_id, ".", "bdl", &(mv->band_limit), 1);
  testErrorRetVA(status<0,hdf5_base,"cannot save /mvdens/bdl in file %s (got %d)",*err,__LINE__,,fname,status);
  
  status = H5LTset_attribute_int( group_id, ".", "chol", &(mv->chol), 1);
  testErrorRetVA(status<0,hdf5_base,"cannot save /mvdens/chol in file %s (got %d)",*err,__LINE__,,fname,status);
  
  dims[0] = mv->ndim;
  dims[1] = mv->ndim;
  
  status = H5LTmake_dataset(group_id,"mean",1,dims,H5T_NATIVE_DOUBLE,mv->mean);
  testErrorRetVA(status<0,hdf5_base,"cannot save /mvdens/mean in file %s (got %d)",*err,__LINE__,,fname,status);

  status = H5LTmake_dataset(group_id,"var",2,dims,H5T_NATIVE_DOUBLE,mv->std);
  testErrorRetVA(status<0,hdf5_base,"cannot save /mvdens/var in file %s (got %d)",*err,__LINE__,,fname,status);

  status = H5Gclose(group_id);
  testErrorRetVA(status<0,hdf5_base,"cannot save /mvdens in file %s (got %d)",*err,__LINE__,,fname,status);
}

mvdens* mvdens_hdfdwnp_infile(mvdens* mv,hid_t loc_id,error **err) {
  hid_t       group_id;   
  herr_t      status;
  char fname[200];
  int ldim;
  mvdens *rmv;
  
  H5Iget_name( loc_id, fname, 200 );
  
  group_id = H5Gopen(loc_id, "mvdens", H5P_DEFAULT);
  testErrorRetVA(group_id<0,hdf5_base,"cannot open group mvdens in file %s (got %d)",*err,__LINE__,NULL,fname,group_id);
  
  status = H5LTget_attribute_int( group_id, ".", "ndim",  &ldim);
  testErrorRetVA(status<0,hdf5_base,"cannot read /mvdens/ndim in file %s (got %d)",*err,__LINE__,NULL,fname,status);
  
  if (mv == NULL) {
    rmv = mvdens_alloc(ldim,err);
    forwardError(*err,__LINE__,NULL);
  } else {
    rmv = mv;
  }
  testErrorRetVA(ldim!=rmv->ndim,mv_dimension,"bad ndim (got %d expected %d)",*err,__LINE__,NULL,ldim,rmv->ndim);
  
  status = H5LTget_attribute_int( group_id, ".", "df", &(mv->df));
  testErrorRetVA(status<0,hdf5_base,"cannot read /mvdens/df in file %s (got %d)",*err,__LINE__,NULL,fname,status);
  
  status = H5LTget_attribute_int( group_id, ".", "bdl", &(mv->band_limit));
  testErrorRetVA(status<0,hdf5_base,"cannot read /mvdens/bdl in file %s (got %d)",*err,__LINE__,NULL,fname,status);
  
  status = H5LTget_attribute_int( group_id, ".", "chol", &(mv->chol));
  testErrorRetVA(status<0,hdf5_base,"cannot read /mvdens/chol in file %s (got %d)",*err,__LINE__,NULL,fname,status);
  
  status = H5LTread_dataset_double (group_id, "mean",rmv->mean );
  testErrorRetVA(status<0,hdf5_base,"cannot read /mvdens/mean in file %s (got %d)",*err,__LINE__,NULL,fname,status);

  status = H5LTread_dataset_double (group_id, "var",rmv->std );
  testErrorRetVA(status<0,hdf5_base,"cannot read /mvdens/var in file %s (got %d)",*err,__LINE__,NULL,fname,status);
  
  status = H5Gclose(group_id);
  testErrorRetVA(status<0,hdf5_base,"cannot save /mvdens in file %s (got %d)",*err,__LINE__,NULL,fname,status);
  if (rmv->chol ==1) {
    rmv->detL = determinant(rmv->std, rmv->ndim);
  }
  
  return rmv;
}


void mvdens_hdfdump(mvdens* mv,char *fname,error **err) {
  hid_t file_id;
  herr_t      status;
  
  /* Create a new file using default properties. */
  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  testErrorRetVA(file_id<0,hdf5_base,"cannot create hdf5 file %s (got %d)",*err,__LINE__,,fname,file_id);
  
  mvdens_hdfdump_infile(mv,file_id,err);
  forwardError(*err,__LINE__,);
  
  status = H5Fclose(file_id);
  testErrorRetVA(status<0,hdf5_base,"cannot save file %s (got %d)",*err,__LINE__,,fname,status);
}

mvdens* mvdens_hdfdwnp(char *fname,error **err) {
  hid_t file_id;
  herr_t      status;
  mvdens *mv;
  
  file_id = H5Fopen( fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  testErrorRetVA(file_id<0,hdf5_base,"cannot open  file %s (got %d)",*err,__LINE__,NULL,fname,file_id);
  
  mv = mvdens_hdfdwnp_infile(NULL,file_id,err);
  forwardError(*err,__LINE__,NULL);
  
  status = H5Fclose(file_id);
  testErrorRetVA(status<0,hdf5_base,"cannot save file %s (got %d)",*err,__LINE__,NULL,fname,status);
  return mv;
}


void mix_mvdens_hdfdump(mix_mvdens* mmv,char *fname,error **err) {
  hid_t file_id,group_id,comp_id;
  herr_t      status;
  int i;
  
  /* Create a new file using default properties. */
  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  testErrorRetVA(file_id<0,hdf5_base,"cannot create hdf5 file %s (got %d)",*err,__LINE__,,fname,file_id);
  
  group_id = H5Gcreate( file_id, "mix_mvdens", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  testErrorRetVA(group_id<0,hdf5_base,"cannot create group /mix_mvdens in file %s (got %d)",*err,__LINE__,,fname,group_id);
  
  status = H5LTset_attribute_int( group_id, ".", "ndim", &(mmv->ndim), 1);
  testErrorRetVA(status<0,hdf5_base,"cannot save /mvdens/ndim in file %s (got %d)",*err,__LINE__,,fname,status);

  status = H5LTset_attribute_int( group_id, ".", "ncomp", &(mmv->ncomp), 1);
  testErrorRetVA(status<0,hdf5_base,"cannot save /mvdens/ncomp in file %s (got %d)",*err,__LINE__,,fname,status);
  
  status = H5LTset_attribute_double( group_id, ".", "weight", mmv->wght, mmv->ncomp);
  testErrorRetVA(status<0,hdf5_base,"cannot save /mvdens/weight in file %s (got %d)",*err,__LINE__,,fname,status);
  
  for(i=0;i<mmv->ncomp;i++) {
    char scomp[50];
    sprintf(scomp,"component_%d",i);
    comp_id = H5Gcreate( group_id, scomp, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    testErrorRetVA(comp_id<0,hdf5_base,"cannot create group /mix_mvdens/%s in file %s (got %d)",*err,__LINE__,,scomp,fname,comp_id);
    mvdens_hdfdump_infile(mmv->comp[i],comp_id,err);
    forwardError(*err,__LINE__,);
    status = H5Gclose(comp_id);
    testErrorRetVA(status<0,hdf5_base,"cannot save /mix_mvdens/%s in file %s (got %d)",*err,__LINE__,,scomp,fname,status);
  }
  status = H5Gclose(group_id);
  testErrorRetVA(status<0,hdf5_base,"cannot save /mix_mvdens in file %s (got %d)",*err,__LINE__,,fname,status);
  
  status = H5Fclose(file_id);
  testErrorRetVA(status<0,hdf5_base,"cannot save file %s (got %d)",*err,__LINE__,,fname,status);
}

mix_mvdens* mix_mvdens_hdfdwnp(char *fname,error **err) {
  mix_mvdens *mmv;
  hid_t file_id,group_id,comp_id;
  herr_t      status;
  int ncomp,ndim,i;
  
  file_id = H5Fopen( fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  testErrorRetVA(file_id<0,hdf5_base,"cannot open  file %s (got %d)",*err,__LINE__,NULL,fname,file_id);
  
  group_id = H5Gopen(file_id, "mvdens", H5P_DEFAULT);
  testErrorRetVA(group_id<0,hdf5_base,"cannot open group mvdens in file %s (got %d)",*err,__LINE__,NULL,fname,group_id);
  
  status = H5LTget_attribute_int( group_id, ".", "ndim",  &ndim);
  testErrorRetVA(status<0,hdf5_base,"cannot read /mvdens/ndim in file %s (got %d)",*err,__LINE__,NULL,fname,status);

  status = H5LTget_attribute_int( group_id, ".", "ncomp",  &ncomp);
  testErrorRetVA(status<0,hdf5_base,"cannot read /mvdens/ncomp in file %s (got %d)",*err,__LINE__,NULL,fname,status);
  
  mmv = mix_mvdens_alloc(ncomp,ndim,err);
  forwardError(*err,__LINE__,NULL);
  
  status = H5LTget_attribute_double( group_id, ".", "weight", mmv->wght);
  testErrorRetVA(status<0,hdf5_base,"cannot read /mvdens/weight in file %s (got %d)",*err,__LINE__,NULL,fname,status);
  
  for(i=0;i<ncomp;i++) {
    char scomp[50];
    sprintf(scomp,"component_%d",i);
    comp_id = H5Gopen( group_id, scomp, H5P_DEFAULT);
    testErrorRetVA(comp_id<0,hdf5_base,"cannot open group /mix_mvdens/%s in file %s (got %d)",*err,__LINE__,NULL,scomp,fname,comp_id);
    mvdens_hdfdwnp_infile(mmv->comp[i],comp_id,err);
    forwardError(*err,__LINE__,NULL);
    status = H5Gclose(comp_id);
    testErrorRetVA(status<0,hdf5_base,"cannot read /mix_mvdens/%s in file %s (got %d)",*err,__LINE__,NULL,scomp,fname,status);
  }
  status = H5Fclose(file_id);
  testErrorRetVA(status<0,hdf5_base,"cannot read file %s (got %d)",*err,__LINE__,NULL,fname,status);
  

}

#endif
