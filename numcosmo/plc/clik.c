/*
 *  clik.c
 *  lowly_project
 *
 *  Created by Karim Benabed on 16/03/11.
 *  Copyright 2011 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */


#include "clik.h"
#include "pmc.h"

#include "clik_helper.h"

// ARE YOU STILL READING ?

// YOU HAVE BEEN WARNED !
char* clik_get_version(clik_object *clikid,error **_err) {
  char *version_str;
  lklbs *lbs;
  int ilkl;
  _dealwitherr;

  version_str = malloc_err(sizeof(char)*500,err);
  _forwardError(*err,__LINE__,NULL);

  sprintf(version_str,"clik version %s",CLIKSVNVERSION);

  if (clikid!=NULL) {
    lbs = _clik_dig(clikid,err);
    _forwardError(*err,__LINE__,NULL);
    for(ilkl=0;ilkl<lbs->nlkl;ilkl++) {
      sprintf(version_str,"%s\n  %s",version_str,lbs->lkls[ilkl]->version);
    }

  }
  return version_str;
}

clik_object* clik_init(char* hdffilepath, error **_err) {
  int n_lkl,i_lkl;
  int *lmax;
  int sz;
  char cur_lkl[100];
  cmblkl **clkl;
  int cli,n_cl;
  zero_bs* zbs;
  distribution *target;
  parname lkl_type;
  cldf *df,*cdf;
  int hk;

  _dealwitherr;

  df  = cldf_open(hdffilepath,err);
  _forwardError(*err,__LINE__,NULL);
  
  n_lkl = cldf_readint(df,"clik/n_lkl_object",err);
  _forwardError(*err,__LINE__,NULL);
  
  sz = 6;
  lmax = cldf_readintarray(df,"clik/lmax",&sz,err);
  _forwardError(*err,__LINE__,NULL);
  

  clkl = malloc_err(sizeof(cmblkl*)*n_lkl,err);
  _forwardError(*err,__LINE__,NULL);
  
  for (i_lkl=0;i_lkl<n_lkl;i_lkl++) {
    sprintf(cur_lkl,"clik/lkl_%d",i_lkl);
  
    cdf = cldf_openchild(df,cur_lkl,err);
    _forwardError(*err,__LINE__,NULL);

    clkl[i_lkl] = clik_lklobject_init(cdf,err);
    _forwardError(*err,__LINE__,NULL);
    
    cmblkl_check_lmax(clkl[i_lkl],lmax,err);
    _forwardError(*err,__LINE__,NULL);
    
    cldf_close(&cdf);
  }
  
    
  n_cl = 0;
  for(cli=0;cli<6;cli++) {
    n_cl += lmax[cli]+1;
  }
  
  zbs = init_zero_bs(lmax, err);
  _forwardError(*err,__LINE__,NULL);
  
  target = init_multilklbs_distribution(n_cl , clkl,n_lkl,
                                        zbs, &zero_bs_compute, &free_zero_bs, lmax, err);
  _forwardError(*err,__LINE__,NULL);
  
  cdf = cldf_openchild(df,"clik",err);
  _forwardError(*err,__LINE__,NULL);
  
  hk = cldf_haskey(cdf,"default",err);
  _forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    int nepar,ndef,idef,j;
    char *defname;
    int *ldef;
    double *loc;
    parname *pn;
    cldf *def_df;

    def_df = cldf_openchild(cdf,"default",err);
    _forwardError(*err,__LINE__,NULL);
    
    nepar = clik_get_extra_parameter_names(target,&pn,err);
    _forwardError(*err,__LINE__,NULL);
    _testErrorRetVA(nepar==0,hdf5_base,"cannot add defaults without extra parameters",*err,__LINE__,NULL,"");    
    ndef = -1;

    defname = cldf_readstr(def_df,"name",&ndef, err);
    _forwardError(*err,__LINE__,NULL);
    ndef = ndef/256;
    _testErrorRetVA(nepar<ndef,hdf5_base,"too many defaults ! Expected less than %d got %d",*err,__LINE__,NULL,nepar,ndef);

    ldef = malloc_err(sizeof(int)*ndef,err);
    _forwardError(*err,__LINE__,NULL);
    
    for(idef=0;idef<ndef;idef++) {
      ldef[idef] = -1;
      for (j=0;j<nepar;j++) {
        if (strcmp(pn[j],&(defname[256*idef]))==0) {
          ldef[idef] = n_cl+j;
          break;   
        }
      }
      _testErrorRetVA(ldef[idef]==-1,hdf5_base,"Unknown extra parameter %s",*err,__LINE__,NULL,defname[256*idef]);
    }
    free(defname);
    free(pn);

    loc = cldf_readfloatarray(def_df,"loc",&ndef,err);
    _forwardError(*err,__LINE__,NULL);  
    distribution_set_default(target, ndef, ldef, loc,err);
    _forwardError(*err,__LINE__,NULL);  
    cldf_close(&def_df);
  }
  
  hk = cldf_haskey(cdf,"prior",err);
  _forwardError(*err,__LINE__,NULL); 
  if (hk==1) {
    int nepar;
    parname *pn;
    int nprior,iprior,j;
    int *lprior;
    char *priorname;
    int nvar;
    double *loc,*var;
    cldf * prior_df;

    prior_df = cldf_openchild(cdf,"prior",err);
    _forwardError(*err,__LINE__,NULL);
    
    nepar = clik_get_extra_parameter_names(target,&pn,err);
    _forwardError(*err,__LINE__,NULL);
    _testErrorRetVA(nepar==0,hdf5_base,"cannot add a prior without extra parameters",*err,__LINE__,NULL,"");    
    
    nprior = -1;
    priorname = cldf_readstr(prior_df,"name",&nprior, err);
    _forwardError(*err,__LINE__,NULL);
    nprior = nprior/256;
    _testErrorRetVA(nepar<nprior,hdf5_base,"too many priors ! Expected less than %d got %d",*err,__LINE__,NULL,nepar,nprior);
      
    lprior = malloc_err(sizeof(int)*nprior,err);
    _forwardError(*err,__LINE__,NULL);
        
    for(iprior=0;iprior<nprior;iprior++) {
      lprior[iprior] = -1;
      for (j=0;j<nepar;j++) {
        if (strcmp(pn[j],&(priorname[256*iprior]))==0) {
          lprior[iprior] = n_cl+j;
          break;   
        }
      }
      _testErrorRetVA(lprior[iprior]==-1,hdf5_base,"Unknown extra parameter %s",*err,__LINE__,NULL,priorname[256*iprior]);
    }   
      
    free(priorname);
    free(pn);
      
    loc = cldf_readfloatarray(prior_df,"loc",&nprior,err);
    _forwardError(*err,__LINE__,NULL);  
        
    nvar=-1;
    var = cldf_readfloatarray(prior_df,"var",&nvar,err);
    _forwardError(*err,__LINE__,NULL);  
    
    if (nvar==nprior) {
      target = add_gaussian_prior(target, nprior, lprior, loc, var, err);
      _forwardError(*err,__LINE__,NULL);  
    } else if (nvar==nprior*nprior) {
      target = add_gaussian_prior_2(target, nprior, lprior, loc, var, err);
      _forwardError(*err,__LINE__,NULL);        
    } else {
      _testErrorRetVA(1==1,hdf5_base,"I don't feel well",*err,__LINE__,NULL,"");
    }
    free(loc);
    free(var);
    free(lprior);
  
    cldf_close(&prior_df);
  }
  

  {
    char *version;
    version = clik_get_version(target,err);
    _forwardError(*err,__LINE__,NULL);
    
    printf("----\n%s\n",version);
    free(version);
  }

  hk = cldf_haskey(cdf,"check_param",err);
  _forwardError(*err,__LINE__,NULL);
    
  if (hk==1) {
    int npar;
    double *chkp;
    double res,res2;
    npar = clik_get_extra_parameter_names(target,NULL,err) + n_cl;
    _forwardError(*err,__LINE__,NULL);
    
    chkp = cldf_readfloatarray(cdf, "check_param",&npar,err);
    _forwardError(*err,__LINE__,NULL);
    
    res = cldf_readfloat(cdf,"check_value",err);
    _forwardError(*err,__LINE__,NULL);
    
    res2 = clik_compute(target,chkp,err);
    _forwardError(*err,__LINE__,NULL);
    
    printf("Checking likelihood '%s' on test data. got %g expected %g (diff %g)\n",hdffilepath,res2,res,res-res2);
    free(chkp);
  }
  printf("----\n");
  
  cldf_close(&cdf);

  cldf_close(&df);
  free(lmax);
  

  return target;
}

void clik_get_has_cl(clik_object *clikid, int has_cl[6],error **_err) {
  distribution *target;
  lklbs *lbs;
  int cli;
  _dealwitherr;

  lbs = _clik_dig(clikid,err);
  _forwardError(*err,__LINE__,);
  for(cli=0;cli<6;cli++) {
    if (lbs->offset_lmax[cli]!=-1) {
      has_cl[cli]=1;
    } else {
      has_cl[cli]=0;
    }
  }
}

void clik_get_lmax(clik_object *clikid, int lmax[6],error **_err) {
  distribution *target;
  lklbs *lbs;
  zero_bs* zbs;
  int cli;
  _dealwitherr;
  
  lbs = _clik_dig(clikid,err);
  _forwardError(*err,__LINE__,);
  zbs = lbs->rbs->bs;
  
  for(cli=0;cli<6;cli++) {
    lmax[cli] = zbs->lmax[cli];
  }
}

int clik_get_extra_parameter_names(clik_object* clikid, parname **names, error **_err) {
  parname *pn;
  distribution *target;
  lklbs *lbs;
  int i;
  _dealwitherr;
  int n_cl=0;
  int lmax[6];
  int cli;
  int ii;

  target = _clik_dig2(clikid,err);
  _forwardError(*err,__LINE__,-1);
  lbs = _clik_dig(clikid,err);
  _forwardError(*err,__LINE__,-1);
  
  clik_get_lmax(clikid,lmax,err);
  _forwardError(*err,__LINE__,-1);
  
  n_cl = 0;
  for(cli=0;cli<6;cli++) {
    n_cl += lmax[cli]+1;
  }

  ii = 0;
  if (names!=NULL) {
    if (lbs->xdim==0) {
      //for now, no extr parameters
      pn = malloc_err(1*sizeof(parname),err);
      _forwardError(*err,__LINE__,-1);
    } else {
      
      pn = malloc_err(lbs->xdim*sizeof(parname),err);
      //_DEBUGHERE_("%d %d",lbs->xdim,sizeof(parname));
      _forwardError(*err,__LINE__,-1);
    }

    for(i=0;i<lbs->xdim;i++) {
      if  (target->ndef==0 || target->def[i+n_cl]==0) {
        sprintf(pn[ii],"%s",lbs->xnames[i]);
        //_DEBUGHERE_("%d %s %d %s",ii,pn[ii],i,lbs->xnames[i])
        ii++;
      }
    }
    *names = pn;  
  } else {
    for(i=0;i<lbs->xdim;i++) {
      if  (target->ndef==0 || target->def[i+n_cl]==0) {
        ii++;
      }
    }
  }
  return ii;
}

int clik_get_extra_parameter_names_by_lkl(clik_object* clikid, int ilkl,parname **names, error **_err) {
  parname *pn;
  distribution *target;
  lklbs *lbs;
  int i;
  _dealwitherr;

  lbs = _clik_dig(clikid,err);
  _forwardError(*err,__LINE__,-1);
  _testErrorRetVA(ilkl>lbs->nlkl,-11010,"Asked for lkl %d, while there are only %d objects",*err,__LINE__,-1,ilkl,lbs->nlkl);
  
  if (lbs->lkls[ilkl]->xdim==0) {
    //for now, no extr parameters
    pn = malloc_err(1*sizeof(parname),err);
    _forwardError(*err,__LINE__,-1);
  } else {
    pn = malloc_err(lbs->lkls[ilkl]->xdim*sizeof(parname),err);
    _forwardError(*err,__LINE__,-1);
  }
  for(i=0;i<lbs->lkls[ilkl]->xdim;i++) {
    sprintf(pn[i],"%s",lbs->lkls[ilkl]->xnames[i]);
  }
  *names = pn;
  return lbs->lkls[ilkl]->xdim;
}

void clik_cleanup(clik_object** pclikid) {
  free_distribution((distribution **)pclikid);
}

double clik_compute(clik_object* clikid, double* cl_and_pars,error **_err) {
  double res;
  //int i,n,lmax[6],lm;
  //parname *names;
  _dealwitherr;

  //_DEBUGHERE_("","");
  //clik_get_lmax(clikid,  lmax,err);
  //_forwardError(*err,__LINE__,-1);
  //lm = lmax[0]+lmax[1]+lmax[2]+lmax[3]+lmax[4]+lmax[5]+6;
  //n = clik_get_extra_parameter_names(clikid, &names, err);
  //_forwardError(*err,__LINE__,-1);
  //for(i=0;i<n;i++) {
  //  _DEBUGHERE_("%s = %g",names[i],cl_and_pars[lm+i]);
  //}
  //free(names);

  res = distribution_lkl(clikid, cl_and_pars,err);
  _forwardError(*err,__LINE__,-1);
  return res;
}

void* _clik_dig(clik_object* clikid, error **err) {
  distribution *target;
  target = clikid;
  if ((void*)target->log_pdf == (void*)&combine_lkl) { 
    // return the first clik likelihood
    int i;
    comb_dist_data* cbd;
    cbd = target->data;
    for (i=0;i<cbd->ndist;i++) {
      if (cbd->dist[i]->log_pdf == &lklbs_lkl) {
        return cbd->dist[i]->data;
      }
    }
  }
  if (target->log_pdf==&lklbs_lkl) {
    return target->data;
  }
  testErrorRet(1==1,-111,"No clik likelihood found",*err,__LINE__,NULL);
  return NULL;
}

void* _clik_dig2(clik_object* clikid, error **err) {
  distribution *target;
  target = clikid;
  if ((void*)target->log_pdf == (void*)&combine_lkl) { 
    // return the first clik likelihood
    int i;
    comb_dist_data* cbd;
    cbd = target->data;
    for (i=0;i<cbd->ndist;i++) {
      if (cbd->dist[i]->log_pdf == &lklbs_lkl) {
        return cbd->dist[i];
      }
    }
  }
  if (target->log_pdf==&lklbs_lkl) {
    return target;
  }
  testErrorRet(1==1,-111,"No clik likelihood found",*err,__LINE__,NULL);
  return NULL;
}

