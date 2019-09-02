#include "pmc.h"
#include "cldf/cldf.h"
#include "clik_helper.h"

typedef struct {
  int nstepsEE,nstepsBB,nstepsTE; 
  float stepEE,stepBB,stepTE;
  double *probEE,*probBB,*probTE;
  int has_cl[6];
  int lmin,nell;
  int clmin,clmax;
} simall_data;

void free_simall(void **pft) {
  simall_data *ft;
  
  ft = *pft;
  if (ft->has_cl[1]) {
    free(ft->probEE);
  }
  if (ft->has_cl[2]) {
    free(ft->probBB);
  }
  if (ft->has_cl[3]) {
    free(ft->probTE);
  } 
  free(ft);
}

double simall_lkl(void* ot, double *pars, error **err) {
  simall_data *ft;
  double res,dl;
  int ell,il,position,offset;

//  printf("Calc SimAll\n"); 
  
  ft = ot;
  res = 0;

  offset=0;
  
  if (ft->has_cl[1]) {
    //EE
    for(il=0;il<ft->nell;il++) {
      ell=il+2;
      dl=pars[il+offset]*ell*(ell+1)/2/M_PI;
      position = (int) (dl/ft->stepEE);
      testErrorRetVA(position>ft->nstepsEE,-1233,"multipole EE %d too large (got %g expected <%g)",*err,__LINE__,-1e10,ell,dl,ft->stepEE*ft->nstepsEE);
      res+=ft->probEE[position+il*ft->nstepsEE];
//      printf("%i %e %i %e\n",ell,pars[il],position,ft->probEE[position+ft->nell*ft->nstepsEE]);
    }
    offset+=ft->nell;
  }

  if (ft->has_cl[2]) {
    //BB
    for(il=0;il<ft->nell;il++) {
      ell=il+2;
      dl=pars[il+offset]*ell*(ell+1)/2/M_PI;
      position = (int) (dl/ft->stepBB);
      testErrorRetVA(position>ft->nstepsBB,-1233,"multipole BB %d too large (got %g expected <%g)",*err,__LINE__,-1e10,ell,dl,ft->stepBB*ft->nstepsBB);
      res+=ft->probBB[position+il*ft->nstepsBB];
//      printf("%i %e %e %e\n",ell,pars[il],like,res);
    }
    offset+=ft->nell;
  }
  
  if (ft->has_cl[3]) {
    //TE
    for(il=0;il<ft->nell;il++) {
      ell=il+2;
      dl=pars[il+offset]*ell*(ell+1)/2/M_PI;
      position = (int) (dl/ft->stepTE);
      testErrorRetVA(position>ft->nstepsTE,-1233,"multipole TE %d too large (got %g expected <%g)",*err,__LINE__,-1e10,ell,dl,ft->stepTE*ft->nstepsTE);
      res+=ft->probTE[position+il*ft->nstepsTE];
//      printf("%i %e %e %e\n",ell,pars[il],like,res);
    }
  }
  return res;
}


cmblkl* clik_simall_init(cldf * df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  cmblkl *cing;  
  simall_data *ft;
  int idata,ndat;

  printf("Initializing SimAll\n");
  
  ft = malloc_err(sizeof(simall_data),err);
  forwardError(*err,__LINE__,NULL);
  
  for(idata=0;idata<6;idata++) {
    ft->has_cl[idata] = has_cl[idata];
  }
  
  if (ft->has_cl[1]) {
    ft->nstepsEE = cldf_readint(df,"nstepsEE",err);
    forwardError(*err,__LINE__,NULL);
    ft->stepEE = cldf_readfloat(df,"stepEE",err);
    forwardError(*err,__LINE__,NULL);
  }	  
  
  if (ft->has_cl[2]) {
    ft->nstepsBB = cldf_readint(df,"nstepsBB",err);
    forwardError(*err,__LINE__,NULL);
    ft->stepBB = cldf_readfloat(df,"stepBB",err);
    forwardError(*err,__LINE__,NULL);
  }
  
  if (ft->has_cl[3]) {
    ft->nstepsTE = cldf_readint(df,"nstepsTE",err);
    forwardError(*err,__LINE__,NULL);
    ft->stepTE = cldf_readfloat(df,"stepTE",err);
    forwardError(*err,__LINE__,NULL);
  }
  
  ft->lmin = cldf_readint(df,"lmin",err);
  forwardError(*err,__LINE__,NULL);
  
  ft->nell = cldf_readint(df,"nell",err);
  forwardError(*err,__LINE__,NULL);
  
  ft->clmin = ell[0];
  ft->clmax = ell[nell-1];
  testErrorRetVA(ft->clmax>ft->lmin+ft->nell+1,-1233,"lmax too large (got %d expected %d at most)",*err,__LINE__,NULL,ft->clmax,ft->lmin+ft->nell); 
  
  if (ft->has_cl[1]) {
    ndat = ft->nell*ft->nstepsEE;
    ft->probEE = cldf_readfloatarray(df,"probEE",&(ndat),err);
    forwardError(*err,__LINE__,NULL);
  }
  
  if (ft->has_cl[2]) {
    ndat = ft->nell*ft->nstepsBB;
    ft->probBB = cldf_readfloatarray(df,"probBB",&(ndat),err);
    forwardError(*err,__LINE__,NULL);
  }
  
  if (ft->has_cl[3]) {
    ndat = ft->nell*ft->nstepsTE;
    ft->probTE = cldf_readfloatarray(df,"probTE",&(ndat),err);
    forwardError(*err,__LINE__,NULL);
  }
//  printf("SimAll initialized\n");
  
  cing = init_cmblkl(ft, &simall_lkl, 
                     &free_simall,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,0,err);
  
  forwardError(*err,__LINE__,NULL);  
  return cing;  
  
  printf("SimAll initialized\n");
}

