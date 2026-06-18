#include "clik_dic.h"

char* cdic_newchar(char* fromchar,error **err) {
  char * rchar;
  int len;
  if (fromchar==NULL || fromchar[0]=='\0') {
    rchar = malloc_err(sizeof(char)*256,err);
    forwardError(*err,__LINE__,NULL);
    rchar[0]='\0';
    return rchar;
  } 
  len = strlen(fromchar);
  if (len<256) {
    len=256;
  }
  rchar = malloc_err(sizeof(char)*(len+1),err);
  forwardError(*err,__LINE__,NULL);
  sprintf(rchar,"%s",fromchar);
  return rchar;
}

long cdic_hash(char* key) {
  long res;
  int i;
  res = 0;
  //_DEBUGHERE_("%s %d",key,strlen(key));
  for(i=0;key[i]!='\0';i++) {
    //_DEBUGHERE_("%d",i);
    res = res * 31 + key[i];
  } 
  //_DEBUGHERE_("%s -> %d",key,res);
  return res;
}

void cdic_dump(cdic *pf, FILE* strm, error **err) {
  int i;

  for(i=0;i<pf->nkey;i++) {
    fprintf(strm, "%s = %s\n", pf->key[i],pf->value[i]);
  }
}

cdic* cdic_init(error **err) {
  cdic *pf;
 
  pf = malloc_err(sizeof(cdic),err);
  forwardError(*err,__LINE__,NULL);
 
  pf->nmax = 20;
  pf->nkey = 0;

  pf->key = malloc_err(sizeof(char*)*pf->nmax,err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("%p %d %p",pf,pf->nmax,pf->key);
  pf->value = malloc_err(sizeof(char*)*pf->nmax,err);
  forwardError(*err,__LINE__,NULL);
  pf->dvalue = malloc_err(sizeof(double)*pf->nmax,err);
  forwardError(*err,__LINE__,NULL);    
  pf->hash = malloc_err(sizeof(long)*pf->nmax,err);
  forwardError(*err,__LINE__,NULL);    

  return pf;
}

char** cdic_get_keys(cdic *hbf, int *nkeys, error **err) {
  *nkeys = hbf->nkey;
  return hbf->key;
}

void cdic_free(void **ppf) {
  cdic *pf;
  int i;

  pf = *ppf;
  //_DEBUGHERE_("%p %d %p",pf,pf->nmax,pf->key);
  if (pf->nmax !=0) {
    for(i=0;i<pf->nkey;i++) {
      //_DEBUGHERE_("%p",pf->key[i]);
      //_DEBUGHERE_("%s",pf->key[i]);
      free(pf->key[i]);
      free(pf->value[i]);
    }
    free(pf->key);
    free(pf->value);
    free(pf->hash);
  }
  //_DEBUGHERE_("","");
  free(pf);
  *ppf = NULL;
}

void cdic_add_item(cdic* pf,int nit, char** key, char **value,error **err) {
  int i;
  int cur;

  if (pf->nkey+nit>= pf->nmax) {
     //grow list;
    //_DEBUGHERE_("-> %p %d %d %d",pf,pf->nmax, pf->nkey, nit);

    int nnmax =  (pf->nkey+nit)*2;
    pf->key=resize_err(pf->key, sizeof(char*)*pf->nmax, sizeof(char*)*nnmax, 1, err);
    forwardError(*err,__LINE__,);
    pf->value=resize_err(pf->value, sizeof(char*)*pf->nmax, sizeof(char*)*nnmax, 1, err);
    forwardError(*err,__LINE__,);
    pf->hash=resize_err(pf->hash, sizeof(long)*pf->nmax, sizeof(long)*nnmax, 1, err);
    forwardError(*err,__LINE__,);
    //pf->dvalue=resize_err(pf->dvalue, sizeof(double)*pf->nmax, sizeof(hbchar)*nnmax, 1, err);
    //forwardError(*err,__LINE__,);
    pf->nmax =  nnmax;
  }
  //_DEBUGHERE_("-> %p %d %d %d",pf,pf->nmax, pf->nkey, nit);
  cur = pf->nkey;
  for(i=0;i<nit;i++) {
    int idx;
    //_DEBUGHERE_("%s",key[i]);
    //_DEBUGHERE_("%d",strlen(key[i]));
    idx = cdic_key_index(pf,key[i],err);
    forwardError(*err,__LINE__,);
    //_DEBUGHERE_("cur %d idx %d",cur,idx)
    if (idx==-1) {
      idx = cur;
      cur++;
    }
    //_DEBUGHERE_("cur %d idx %d",cur,idx)
    //_DEBUGHERE_("%p",value);
    //_DEBUGHERE_("%p",value[i]);
    //_DEBUGHERE_("%s",value[i]);
    //_DEBUGHERE_("","");
    //_DEBUGHERE_("%s (%d) -> %s (%d)",key[i],strlen(key[i]),value[i],strlen(value[i]));
    pf->hash[idx] = cdic_hash(key[i]);
    pf->key[idx] = cdic_newchar(key[i],err);
    forwardError(*err,__LINE__,);
    pf->value[idx] = cdic_newchar(value[i],err);
    forwardError(*err,__LINE__,);
    
  }
  pf->nkey = cur;  

}

void cdic_remove_item(cdic* pf, int index,error **err) {
  int i;

  free(pf->key[index]);
  free(pf->value[index]);
  for(i=index;i<pf->nkey-1;i++) {
    pf->key[i] = pf->key[i+1];
    pf->value[i] = pf->value[i+1];
    pf->hash[i] = pf->hash[i+1];
  }
  pf->nkey--;
}

int cdic_key_index(cdic *pf, char *key, error **err) {
  int i;
  long hash;

  //_DEBUGHERE_("","");
  hash = cdic_hash(key);
  //_DEBUGHERE_("","");
  
  for(i=0;i<pf->nkey;i++) {
    if (hash==pf->hash[i]) {
      if (strcmp(key,pf->key[i])==0) {
        return i;
      }  
    }
  }
  return -1;
}

char* cdic_get(cdic* pf, char* key, char* safeguard,error **err) {
  int ps;
  
  ps=cdic_key_index(pf,key,err);
  forwardError(*err,__LINE__,NULL);
  if (ps==-1) {
    testErrorRetVA(safeguard==NULL,-1234,"key '%s' absent",*err,__LINE__,NULL,key);
    return safeguard;
  }
  return pf->value[ps];
}

char* cdic_sget(cdic* pf, char* key, char* safeguard,error **err) {
  int ps;
  
  ps=cdic_key_index(pf,key,err);
  forwardError(*err,__LINE__,NULL);
  if (ps==-1) {
    testErrorRetVA(safeguard==NULL,-1234,"key '%s' absent",*err,__LINE__,NULL,key);
    cdic_add_item(pf,1,&key,&safeguard,err);
    forwardError(*err,__LINE__,NULL);
    return safeguard;
  }
  return pf->value[ps];
}

char* cdic_set(cdic* pf, char* key, char* safeguard,error **err) {
  int ps;
  cdic_add_item(pf,1,&key,&safeguard,err);
  forwardError(*err,__LINE__,NULL);
  return safeguard;
}

long cdic_get_int(cdic *pf, char *key,long* safeguard, error **err) {
  char *res;

  res = cdic_get(pf,key,(char*) safeguard,err);
  forwardError(*err,__LINE__,0);
  
  if (res==(char*)safeguard) {
    return *safeguard;
  }
  
  return atol(res);
}

long cdic_sget_int(cdic *pf, char *key,long* safeguard, error **err) {
  char *res;

  res = cdic_get(pf,key,(char*) safeguard,err);
  forwardError(*err,__LINE__,0);
  
  if (res==(char*)safeguard) {
    cdic_tempchar value;
    char *dvalue;
    dvalue = value;
    sprintf(value,"%ld",*safeguard);
    cdic_add_item(pf,1,&key,&dvalue,err);
    forwardError(*err,__LINE__,0); 
    return *safeguard;
  }
  
  return atol(res);
}

long cdic_set_int(cdic *pf, char *key,long safeguard, error **err) {
  cdic_tempchar value;
  char *dvalue;
  dvalue = value;
  sprintf(value,"%ld",safeguard);
  cdic_add_item(pf,1,&key,&dvalue,err);
  forwardError(*err,__LINE__,0); 
  return safeguard;
}


double cdic_get_double(cdic *pf, char *key,double *safeguard, error **err) {
  int ps;

  ps=cdic_key_index(pf,key,err);
  forwardError(*err,__LINE__,M_NAN);

  if (ps==-1) {
    return *safeguard;
  }

  return atof(pf->value[ps]);
}

double cdic_sget_double(cdic *pf, char *key,double *safeguard, error **err) {
  int ps;

  ps=cdic_key_index(pf,key,err);
  forwardError(*err,__LINE__,M_NAN);

  if (ps==-1) {
    cdic_tempchar value;
    char *dvalue;
    dvalue = value;
    sprintf(value,"%g",*safeguard);
    cdic_add_item(pf,1,&key,&dvalue,err);
    forwardError(*err,__LINE__,M_NAN); 
    return *safeguard;
  }

  return atof(pf->value[ps]);
}

double cdic_set_double(cdic *pf, char *key,double safeguard, error **err) {
  cdic_tempchar value;
  char *dvalue;
  dvalue = value;
  sprintf(value,"%g",safeguard);
  cdic_add_item(pf,1,&key,&dvalue,err);
  forwardError(*err,__LINE__,M_NAN); 
  return safeguard;
}


#ifdef ADD0US
void f90_cdic_get_int(int* res,long *lpf, char *fkey, int *lfkey) {
#elif ADD2US
void f90_cdic_get_int__(int* res,long *lpf, char *fkey, int *lfkey) {
#else
void f90_cdic_get_int_(int* res,long *lpf, char *fkey, int *lfkey) {
#endif
  cdic *pf;
  cdic_tempchar key;
  long safe;
  int ires;
  error *_err,**err;

  _err = NULL;
  err = &_err;

  pf = (cdic*) *lpf;
  memcpy(key,fkey,(*lfkey)*sizeof(char));
  key[*lfkey] = '\0';

  safe = *res;
  ires = cdic_sget_int(pf,key,&safe,err);
  quitOnError(*err,__LINE__,stderr);

  *res = ires;
}

#ifdef ADD0US
void f90_cdic_set_int(int* res,long *lpf, char *fkey, int *lfkey) {
#elif ADD2US
void f90_cdic_set_int__(int* res,long *lpf, char *fkey, int *lfkey) {
#else
void f90_cdic_set_int_(int* res,long *lpf, char *fkey, int *lfkey) {
#endif
  cdic *pf;
  cdic_tempchar key;
  long safe;
  int ires;
  error *_err,**err;

  _err = NULL;
  err = &_err;

  pf = (cdic*) *lpf;
  memcpy(key,fkey,(*lfkey)*sizeof(char));
  key[*lfkey] = '\0';

  safe = *res;
  ires = cdic_set_int(pf,key,safe,err);
  quitOnError(*err,__LINE__,stderr);

  *res = ires;
}

#ifdef ADD0US
void f90_cdic_get_double(double* res,long *lpf, char *fkey, int *lfkey) {
#elif ADD2US
void f90_cdic_get_double__(double* res,long *lpf, char *fkey, int *lfkey) {
#else
void f90_cdic_get_double_(double* res,long *lpf, char *fkey, int *lfkey) {
#endif
  cdic *pf;
  cdic_tempchar key;
  double safe;
  double ires;
  error *_err,**err;

  _err = NULL;
  err = &_err;

  pf = (cdic*) *lpf;
  memcpy(key,fkey,(*lfkey)*sizeof(char));
  key[*lfkey] = '\0';

  safe = *res;
  ires = cdic_sget_double(pf,key,&safe,err);
  quitOnError(*err,__LINE__,stderr);
  *res = ires;
}

#ifdef ADD0US
void f90_cdic_set_double(double* res,long *lpf, char *fkey, int *lfkey) {
#elif ADD2US
void f90_cdic_set_double__(double* res,long *lpf, char *fkey, int *lfkey) {
#else
void f90_cdic_set_double_(double* res,long *lpf, char *fkey, int *lfkey) {
#endif
  cdic *pf;
  cdic_tempchar key;
  double safe;
  double ires;
  error *_err,**err;

  _err = NULL;
  err = &_err;

  pf = (cdic*) *lpf;
  memcpy(key,fkey,(*lfkey)*sizeof(char));
  key[*lfkey] = '\0';

  safe = *res;
  ires = cdic_set_double(pf,key,safe,err);
  quitOnError(*err,__LINE__,stderr);
  *res = ires;
}

#ifdef ADD0US
void f90_cdic_get_str(char* res,int *lfres, long *lpf, char *fkey, int *lfkey) {
#elif ADD2US
void f90_cdic_get_str__(char* res,int *lfres, long *lpf, char *fkey, int *lfkey) {
#else
void f90_cdic_get_str_(char* res,int *lfres, long *lpf, char *fkey, int *lfkey) {
#endif
  cdic *pf;
  cdic_tempchar key;
  cdic_tempchar safe;
  char* ires;
  error *_err,**err;
  int len;

  _err = NULL;
  err = &_err;

  pf = (cdic*) *lpf;
  memcpy(key,fkey,(*lfkey)*sizeof(char));
  key[*lfkey] = '\0';

  memcpy(safe,res,(*lfres)*sizeof(char));
  safe[*lfres] = '\0';

  ires = cdic_sget(pf,key,safe,err);
  quitOnError(*err,__LINE__,stderr);

  memset(res,' ',*lfres);
  len = strlen(ires);
  memcpy(res,ires,len*sizeof(char));
  *lfres = len;

}

#ifdef ADD0US
void f90_cdic_set_str(char* res,int *lfres, long *lpf, char *fkey, int *lfkey) {
#elif ADD2US
void f90_cdic_set_str__(char* res,int *lfres, long *lpf, char *fkey, int *lfkey) {
#else
void f90_cdic_set_str_(char* res,int *lfres, long *lpf, char *fkey, int *lfkey) {
#endif
  cdic *pf;
  cdic_tempchar key;
  cdic_tempchar safe;
  char* ires;
  error *_err,**err;
  int len;

  _err = NULL;
  err = &_err;

  pf = (cdic*) *lpf;
  memcpy(key,fkey,(*lfkey)*sizeof(char));
  key[*lfkey] = '\0';

  memcpy(safe,res,(*lfres)*sizeof(char));
  safe[*lfres] = '\0';

  ires = cdic_set(pf,key,safe,err);
  quitOnError(*err,__LINE__,stderr);

  memset(res,' ',*lfres);
  len = strlen(ires);
  memcpy(res,ires,len*sizeof(char));
  *lfres = len;

}

#ifdef ADD0US
void f90_cdic_test(long *lpf) {
#elif ADD2US
void f90_cdic_test__(long *lpf) {
#else
void f90_cdic_test_(long *lpf) {
#endif
  cdic *pf;
  error *_err,**err;
  cdic_tempchar key,value;
  char *dkey,*dvalue;

  _err = NULL;
  err = &_err;

  dkey = (char*)key;
  dvalue = (char*)value;

  pf = cdic_init(err);
  quitOnError(*err,__LINE__,stderr);

  sprintf(key,"chr");
  sprintf(value,"toto");

  cdic_add_item(pf,1, &(dkey), &(dvalue), err);
  quitOnError(*err,__LINE__,stderr);

  sprintf(key,"int");
  sprintf(value,"1");

  cdic_add_item(pf,1, &dkey, &dvalue, err);
  quitOnError(*err,__LINE__,stderr);

  sprintf(key,"flt");
  sprintf(value,"1.3455");

  cdic_add_item(pf,1, &dkey, &dvalue, err);
  quitOnError(*err,__LINE__,stderr);

  *lpf = (long) pf;
}


#ifdef ADD0US
void f90_cdic_free(long *lpf) {
#elif ADD2US
void f90_cdic_free__(long *lpf) {
#else
void f90_cdic_free_(long *lpf) {
#endif
  cdic *pf;
  pf = (cdic*) *lpf;
  cdic_free((void**) &pf);
}

#ifdef ADD0US
void f90_cdic_dump(long *lpf) {
#elif ADD2US
void f90_cdic_dump__(long *lpf) {
#else
void f90_cdic_dump_(long *lpf) {
#endif
  cdic *pf;
  error *_err,**err;

  _err = NULL;
  err = &_err;

  pf = (cdic*) *lpf;
  cdic_dump(pf,stdout,err);
  quitOnError(*err,__LINE__,stderr);
}