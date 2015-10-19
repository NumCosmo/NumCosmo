#include "cldf.h"


void _read_meta(cldf *df,error **err) {
  char cuf[4096];
  char li[8192];
  char *li2,*li3;
  int ls;
  int nmax;
  FILE *ff;


  if (df->root[0]=='\0') {
    sprintf(cuf,"%s/%s",df->name,MDB);  
  } else {
    sprintf(cuf,"%s/%s/%s",df->root,df->name,MDB);  
  }
  
  ff = fopen_err(cuf,"r",err);
  forwardError(*err,__LINE__,);
  nmax = 1;

  df->metakey = malloc_err(sizeof(kv)*nmax,err);
  forwardError(*err,__LINE__,);
  df->metatype = malloc_err(sizeof(kv)*nmax,err);
  forwardError(*err,__LINE__,);
  df->metavalue = malloc_err(sizeof(kv)*nmax,err);
  forwardError(*err,__LINE__,);
  df->nmeta = 0;
  //_DEBUGHERE_("","");  
  while(1) {
    int wh;
    ls = read_line(ff, li, 8192, err);
    testErrorRet(ls==-1,-1234,"bad line",*err,__LINE__,);  
    if (ls==1) {
      break;
    }
    
    if (df->nmeta==nmax) {
      //_DEBUGHERE_("","");  
      df->metakey=resize_err(df->metakey, sizeof(kv)*nmax, sizeof(kv)*nmax*2, 1, err);
      forwardError(*err,__LINE__,);
      df->metavalue = resize_err(df->metavalue, sizeof(kv)*nmax, sizeof(kv)*nmax*2, 1, err);
      forwardError(*err,__LINE__,);
      df->metatype = resize_err(df->metatype, sizeof(kv)*nmax, sizeof(kv)*nmax*2, 1, err);
      forwardError(*err,__LINE__,);
      nmax*=2;
    }

    for(wh=0;wh<ls;wh++) {
      if (li[wh]==' ') break;
    }
    li[wh]='\0';
    strcpy(df->metakey[df->nmeta],li);

    li2 = &(li[wh+1]);
    ls -= wh+1;
    for(wh=0;wh<ls;wh++) {
      if (li2[wh]==' ') break;
    }
    li2[wh]='\0';
    strcpy(df->metatype[df->nmeta],li2);

    li3 = &(li2[wh+1]);
    ls -= wh+1;
/*    for(wh=0;wh<ls;wh++) {
      if (li3[wh]==' ') break;
    }
    li3[wh]='\0';*/
    strcpy(df->metavalue[df->nmeta],li3);
    df->nmeta++;

  }
  fclose(ff);
}



cldf * cldf_open_sub(char *path, char* sub,error **err) {
  cldf * df;
  struct stat buf;
  int erri;
  DIR *dd;
  struct dirent *dc;
  char *mpath;
  char mmpath[8192];
  char lpath[8192];
  int maxchild,maxdata;
  struct stat m;

  df = malloc_err(sizeof(cldf),err);
  forwardError(*err,__LINE__,NULL);

#ifdef HDF5_COMPAT_MODE
  df->hid = 0;
#endif

  df->root[0] = '\0';

  mpath = path;
  if (sub!=NULL && sub[0]!='\0') {
    strcpy(df->root,sub);  
    mpath = mmpath;
    sprintf(mpath,"%s/%s",sub,path);
  }

  erri = stat(mpath, &buf);
  testErrorRetVA(erri!=0,-1234,"cannot stat file (error %d)", *err,__LINE__,NULL,erri);
  testErrorRetVA(S_ISDIR(buf.st_mode)==0,-1234,"%s in not a directory",*err,__LINE__,NULL,path);

  strcpy(df->name,path);

  _read_meta(df,err);
  forwardError(*err,__LINE__,NULL);

  dd = opendir(mpath);
  testErrorRet(dd==NULL,-1234,"bad bad bad",*err,__LINE__,NULL);
  maxchild = 20;
  maxdata = 20;
  
  df->datakey = malloc_err(sizeof(kv)*maxdata,err);
  forwardError(*err,__LINE__,NULL);
  df->ndata=0;

  df->child = malloc_err(sizeof(void*)*maxchild,err);
  forwardError(*err,__LINE__,NULL);
  df->nchild=0;

  while(1) {
    dc = readdir(dd);
    if (dc==NULL) {
      break;
    }
    
    sprintf(lpath,"%s/%s",mpath,dc->d_name);
    
    if (dc->d_name[0]=='.' || dc->d_name[0]=='_') {
      continue;
    }
    testErrorRetVA(stat(lpath,&m)==-1,-1234,"Can't stat %s (%s)",*err,__LINE__,NULL,strerror(errno));
    if (S_ISDIR(m.st_mode)) {
      if (df->nchild==maxchild) {
        resize_err(df->child, sizeof(void*)*maxchild, sizeof(void*)*maxchild*2, 1, err);
        forwardError(*err,__LINE__,NULL);
        maxchild*=2;
      }
      df->child[df->nchild] = cldf_open_sub(dc->d_name,mpath,err);
      forwardError(*err,__LINE__,NULL);      
      df->nchild++;
      continue;
    }
    if (df->ndata==maxdata) {
      resize_err(df->datakey, sizeof(void*)*maxdata, sizeof(void*)*maxdata*2, 1, err);
      forwardError(*err,__LINE__,NULL);
      maxdata*=2;
    }
    strcpy(df->datakey[df->ndata],dc->d_name);
    df->ndata++;
  }
  return df;
}

cldf * cldf_open(char *path, error **err) {
  cldf *df;

#ifdef HDF5_COMPAT_MODE
  struct stat buf;
  int erri;
  hid_t file_id;
   
  erri = stat(path, &buf);
  testErrorRetVA(erri!=0,-1234,"cannot stat file (error %d)", *err,__LINE__,NULL,erri);

  if (S_ISDIR(buf.st_mode)==0) {
    df = malloc_err(sizeof(cldf),err);
    forwardError(*err,__LINE__,NULL);

    file_id = H5Fopen( path, H5F_ACC_RDONLY, H5P_DEFAULT);
    testErrorRetVA(file_id<0,-331234,"cannot open  file %s (got %d)",*err,__LINE__,NULL,path,file_id);
  
  
    df->hid = file_id;
    df->isgroup = 0;
  
    sprintf(df->root,"%s",path);
    return df;
  }
#endif
  
  df = cldf_open_sub(path,NULL,err);
  forwardError(*err,__LINE__,NULL);
  return df;
}

cldf* cldf_tolast(cldf* df, char* key, char* kp, error **err) {
  int ls;
  int l0;
  int lk,i;
  cldf *cdf;

  l0 = 0;
  ls = 0;
  lk = strlen(key);
  cdf = df;

  while(1) {                
    if (key[lk]=='/') {     
      lk--;                 
    } else {                
      break;                
    }                       
  }             
  
  while(1) {
      
    if (key[ls]=='\0' || key[ls]=='/') {
      if (ls==l0) {
        l0++;
        ls=l0;
        continue;
      }
      memcpy(kp,&(key[l0]),sizeof(char)*(ls-l0));
      kp[ls-l0] = '\0';
      l0 = ls+1;
      ls = l0;
      
      // tester si c'est le dernier element
      if (ls>=lk) { //dernier element !
        return cdf;
      }
      // pas le dernier, deplace !
      for(i=0;i<cdf->nchild;i++) {
        if (strcmp(kp,((cldf*)cdf->child[i])->name)==0) {
            cdf = cdf->child[i];
            break;
        }
      }
    }
    ls++;
  }
}


#ifdef HDF5_COMPAT_MODE
int cldf_cut_key(char* key, char *pk) {
  herr_t hstat;
  int lk;
  int value;

  lk = strlen(key);
  while(lk>0) {                
    if (key[lk]!='/') {     
      lk--;                 
    } else {                
      break;                
    }                       
  }
  if (lk==0) {
    sprintf(pk,".");
    return -1;
  }             
  memcpy(pk,key,sizeof(char)*lk);
  pk[lk]='\0';
  return lk;      
}
int cldf_has_attribute(cldf* df, char* pk, char* key, int lk, error **err) {
  hid_t group_id;
  herr_t hstat;

  group_id = H5Oopen( df->hid, pk, H5P_DEFAULT);  
  testErrorRetVA(group_id<0,-331234,"cannot read object %s in %s (got %d)",*err,__LINE__,0,pk,df->root,group_id);
  hstat = H5LTfind_attribute(group_id, &(key[lk+1]));
  H5Oclose(group_id);  
  return hstat;
}

#endif

int cldf_haskey(cldf *df, char *key, error **err) {
  char kp[256];
  int i;
  cldf * cdf;
  
#ifdef HDF5_COMPAT_MODE
  if (df->hid!=0) {
      herr_t hstat;
      int lk;
      char pk[1000];
      int value;
      int hk;

    hstat = H5Lexists(df->hid, key, H5P_DEFAULT);
    if (hstat==1) {
      return 1;
    }
    lk = cldf_cut_key(key,pk);
    hk = cldf_has_attribute(df,pk,key,lk,err);
    forwardError(*err, __LINE__,0);
    return hk;
  }  
  
#endif

  cdf = cldf_tolast(df,key,kp,err);
  forwardError(*err,__LINE__,0);

  for(i=0;i<cdf->nmeta;i++) {
    if (strcmp(kp,cdf->metakey[i])==0) {
      return 1;
    }
  }
  for(i=0;i<cdf->ndata;i++) {
    if (strcmp(kp,cdf->datakey[i])==0) {
      return 1;
    }
  }
  for(i=0;i<cdf->nchild;i++) {
    if (strcmp(kp,((cldf*) cdf->child[i])->name)==0) {
      return 1;
    }
  }
  return 0;
}

long cldf_readint_default(cldf *df, char *key,long def, error **err) {
  int hk;
  long res;

  res = def;
  hk = cldf_haskey(df,key,err);
  forwardError(*err,__LINE__,def);
  if (hk == 1) {
    res = cldf_readint(df,key,err);
    forwardError(*err,__LINE__,def);
  }
  return res;
}

long cldf_readint(cldf *df, char *key, error **err) {
  char kp[256];
  int i;
  cldf * cdf;

#ifdef HDF5_COMPAT_MODE
  if (df->hid!=0) {
      herr_t hstat;
      int lk;
      char pk[1000];
      int value;
      int hk;

      lk = cldf_cut_key(key,pk);
      hk = cldf_has_attribute(df,pk,key,lk,err);
      forwardError(*err, __LINE__,0);
      if (hk == 1) {
        hstat = H5LTget_attribute_int( df->hid, pk, &(key[lk+1]),  &value);
        testErrorRetVA(hstat<0,-331234,"cannot read %s int in file %s (got %d)",*err,__LINE__,0,key,df->root,hstat);
        return value;  
      } else {
        hstat = H5LTread_dataset_int( df->hid, key,  &value);
        testErrorRetVA(hstat<0,-331234,"cannot read %s int in file %s (got %d)",*err,__LINE__,0,key,df->root,hstat);
        return value;  
      }
  }
#endif
  
  cdf = cldf_tolast(df,key,kp,err);
  forwardError(*err,__LINE__,0);

  for(i=0;i<cdf->nmeta;i++) {
    if (strcmp(kp,cdf->metakey[i])==0) {
      testErrorRetVA(strcmp(cdf->metatype[i],"int")!=0,-1234,"bad type for element '%s'",*err,__LINE__,0,key);
      return atol(cdf->metavalue[i]);          
    }
  }
  testErrorRetVA(1==1,-1234,"unknown element '%s'",*err,__LINE__,0,key);
  return 0;
}

double cldf_readfloat_default(cldf *df, char *key,double def, error **err) {
  int hk;
  double res;

  res = def;
  hk = cldf_haskey(df,key,err);
  forwardError(*err,__LINE__,def);
  if (hk == 1) {
    res = cldf_readfloat(df,key,err);
    forwardError(*err,__LINE__,def);
  }
  return res;
}

double cldf_readfloat(cldf *df, char *key, error **err) {
  char kp[256];
  int i;
  cldf * cdf;
  double res;

#ifdef HDF5_COMPAT_MODE
  if (df->hid!=0) {
      herr_t hstat;
      int lk;
      char pk[1000];
      double value;
      int hk;

      lk = cldf_cut_key(key,pk);
      hk = cldf_has_attribute(df,pk,key,lk,err);
      forwardError(*err, __LINE__,0);
      if (hk == 1) {
        hstat = H5LTget_attribute_double( df->hid, pk, &(key[lk+1]),  &value);
        testErrorRetVA(hstat<0,-331234,"cannot read %s int in file %s (got %d)",*err,__LINE__,0,key,df->root,hstat);
        return value;  
      } else {
        hstat = H5LTread_dataset_double( df->hid, key,  &value);
        testErrorRetVA(hstat<0,-331234,"cannot read %s int in file %s (got %d)",*err,__LINE__,0,key,df->root,hstat);
        return value;  
      }
  }
#endif

  cdf = cldf_tolast(df,key,kp,err);
  forwardError(*err,__LINE__,0);

  for(i=0;i<cdf->nmeta;i++) {
  if (strcmp(kp,cdf->metakey[i])==0) {
      testErrorRetVA(strcmp(cdf->metatype[i],"str")==0,-1234,"bad type for element '%s'",*err,__LINE__,0,key);
      testErrorRet(sscanf(cdf->metavalue[i],"%lg",&res)!=1,-1234,"gloups",*err,__LINE__,0);      
      return res;
    }
  }
  testErrorRetVA(1==1,-1234,"unknown element '%s'",*err,__LINE__,0,key);
  return 0;
}

/*char *cldf_readstrarray(cldf *df, char *key, int *sz, error **err) {
  int hk;
  char kt[1000];
  char *typ;
  int tsz;
  char **res;
  char *ves;
  int vsz;
  int i;

  sprintf(kt,"%s__type__",key);
  hk =  cldf_haskey(df, kt, err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    typ = cldf_readstr(df,kt,NULL,err);
    forwardError(*err,__LINE__,NULL);
    testErrorRetVA(strcmp(typ,"str_array"),-24235,"%s not a str array",*err,__LINE__,NULL,key);
    ves = cldf_readstr(df,key,&vsz,err);
    forwardError(*err,__LINE__,NULL);
    for(i=0;i<vsz;i++) {
      if (ves[i]=='\n') {
        ves[i] = '\0';
        sscanf(ves,"%d",sz);
        i++;
        break;
      }
    }
    testErrorRetVA(i==vsz,-22311,"%d not a str array",*err,__LINE__,NULL,key);
    res = malloc_err(sizeof(char)*vsz+sizeof(char*)**sz,err);
    cn = 0;
    ci = 0;
    res[cn] = (((void*) res) + *sz*sizeof(char*) + ci*sizeof(char));
    ip = i;
    for(;i<vsz;i++) {
      if (ves[i]=='\n') {
        ves[i] = '\0';
        sscanf(&(ves[ip]),"%d",&gn);
        i++;
        ves[i+gn]='\0'
        memcpy(&(res(cn)),&(ves[i]),sizeof(char)*(gn+1));
        i+=gn;
        cn++;
        ci+=gn+1;
        ip=i+1;
        res[cn] = (((void*) res) + *sz*sizeof(char*) + ci*sizeof(char));    
      }
    }
    free(ves);
    return res;    
  }


  ves = cldf_readstr(df,key,&vsz,err);
  forwardError(*err,__LINE__,NULL);
  testErrorRetVA(vsz%256!=0,-242352,"%s is not a str array (size %d)",*err,__LINE__,NULL,key,vsz);

  *sz = vsz/256;
  res = malloc_err(sizeof(char)*vsz+sizeof(char*)**sz,err);
  forwardError(*err,__LINE__,NULL);

  for(i=0;i<*sz;i++) {
    res[i] = (((void*) res) + *sz*sizeof(char*) + i*256*sizeof(char));
    sprintf(res[i],"%s",&(ves[i*256]));
  }    
  free(ves);
  return res;

  
}
*/
char* cldf_readstr(cldf *df, char *key, int *sz,error **err) {
  char kp[256];
  int i;
  cldf * cdf;
  char* res;
  int nsz;

  if (sz==NULL) {
    nsz = -1;
  } else {
    nsz = *sz;
  }

  #ifdef HDF5_COMPAT_MODE
  if (df->hid!=0) {
    herr_t hstat;
    int lk;
    char pk[1000];
    hsize_t ndum;
    H5T_class_t dum;
    size_t ddum;
    char *res;
    int j;
    int hk;

    lk = cldf_cut_key(key,pk);
    hk = cldf_has_attribute(df,pk,key,lk,err);
    forwardError(*err,__LINE__,NULL);
      
    if (hk == 1) {
      hstat = H5LTget_attribute_info( df->hid, pk,&(key[lk+1]), &ddum, &dum, &ndum);
      testErrorRetVA(hstat<0,-331234,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,key,df->root,hstat);
      testErrorRetVA((ndum!=nsz && nsz>0),-331234,"Bad size for %s in %s (got %d expected %d)",*err,__LINE__,NULL,key,df->root,ndum,nsz);
      res = malloc_err(sizeof(char)*(ndum+1),err);
      forwardError(*err,__LINE__,NULL);
      hstat = H5LTget_attribute_string(df->hid, pk,&(key[lk+1]),res);
      testErrorRetVA(hstat<0,-331234,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,key,df->root,hstat);
      if (nsz<0 && sz!=NULL) {
        *sz = ndum;
      }
      res[ndum] = '\0';
      return res;
    } else {
      hstat = H5LTget_dataset_info(df->hid, key, &ddum, &dum, &ndum);
      testErrorRetVA(hstat<0,-331234,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,key,df->root,hstat);
      testErrorRetVA((ndum!=nsz && nsz>0),-331234,"Bad size for %s in %s (got %d expected %d)",*err,__LINE__,NULL,key,df->root,ndum,nsz);
      res = malloc_err(sizeof(char)*(ndum+1),err);
      forwardError(*err,__LINE__,NULL);
      hstat = H5LTread_dataset_string(df->hid, key,res);
      testErrorRetVA(hstat<0,-331234,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,key,df->root,hstat);
      if (nsz<0 && sz !=NULL) {
        *sz = ndum;
      }
      res[ndum] = '\0';
      return res;
    }

  }
#endif

  cdf = cldf_tolast(df,key,kp,err);
  forwardError(*err,__LINE__,0);


  for(i=0;i<cdf->nmeta;i++) {
    if (strcmp(kp,cdf->metakey[i])==0) {
      testErrorRetVA(strcmp(cdf->metatype[i],"str")!=0,-1234,"bad type for element '%s'",*err,__LINE__,0,key);
      res = malloc_err(sizeof(char)*(strlen(cdf->metavalue[i])+1),err);
      forwardError(*err,__LINE__,NULL);
      strcpy(res,cdf->metavalue[i]);
      return res;
    }
  }
  for(i=0;i<cdf->ndata;i++) {
    if (strcmp(kp,cdf->datakey[i])==0) {
      FILE *f;
      struct stat st;
      char pth[8192];
      size_t size;
      if (cdf->root[0]=='\0') {
        sprintf(pth,"%s/%s",cdf->name,kp);  
      } else {
        sprintf(pth,"%s/%s/%s",cdf->root,cdf->name,kp);  
      }

      stat(pth, &st);
      size = st.st_size;
      testErrorRetVA((size!=nsz && nsz>0),-331234,"Bad size for %s in %s (got %d expected %d)",*err,__LINE__,NULL,key,df->root,size,nsz);
      res = malloc_err(size+1,err);
      forwardError(*err,__LINE__,NULL);
      f = fopen_err(pth,"r",err);
      fread(res, 1, size, f);
      res[size] = '\0';
      fclose(f);
      if (nsz<0 && sz !=NULL) {
        *sz = size;
      }
      //_DEBUGHERE_("%s %d",key,size);
      return res; 
    }
  }
  testErrorRetVA(1==1,-1234,"unknown element '%s'",*err,__LINE__,0,key);
  return 0;
}


cldf* cldf_openchild(cldf *df, char* key, error **err) {
  char kp[256];
  int i;
  cldf *cdf;

#ifdef HDF5_COMPAT_MODE
  if (df->hid!=0) {
    cdf = malloc_err(sizeof(cldf),err);
    forwardError(*err,__LINE__,NULL);
    
    cdf->hid = H5Gopen(df->hid, key, H5P_DEFAULT );
    testErrorRetVA(cdf->hid<0,-331234,"cannot open  group %s in %s (got %d)",*err,__LINE__,NULL,key,df->root,cdf->hid);
    
    sprintf(cdf->root,"%s/%s",df->root,key);    
    cdf->isgroup=1;
    return cdf;
  }
#endif
  cdf = cldf_tolast(df,key,kp,err);
  forwardError(*err,__LINE__,0);

  for(i=0;i<cdf->nchild;i++) {
    if (strcmp(kp,((cldf*)cdf->child[i])->name)==0) {
      return cdf->child[i];
    }
  }
  testErrorRetVA(1==1,-1234,"unknown element '%s'",*err,__LINE__,0,key);
  return 0;
}

void cldf_close(cldf** pdf) {
  cldf *df,*cdf;
  int i;

  df = *pdf;
  //_DEBUGHERE_("","");
#ifdef HDF5_COMPAT_MODE
  if(df->hid!=0) {
    if (df->isgroup==1) {
      H5Gclose(df->hid);
    } else {
      H5Fclose(df->hid);
    }
    free(df);
    *pdf = NULL;
    return;
  }
#endif
  //_DEBUGHERE_("","");
  df = *pdf;
  if (df->root[0] == '\0') {
    //_DEBUGHERE_("","");
    free(df->metavalue);
    free(df->metakey);
    free(df->metatype);
    free(df->datakey);
    //_DEBUGHERE_("","");
    for(i=0;i<df->nchild;i++) {
      //_DEBUGHERE_("","");
      cdf = df->child[i];
      cdf->root[0] = '\0';
      cldf_close(&cdf);
      //_DEBUGHERE_("","");
    } 
    //_DEBUGHERE_("","");
  }
  *pdf = NULL;
}

void* cldf_readanyfits(char* path, int *sz, int typ, error **err) {
  void * res;
  int fitserr;
  fitsfile *fitsptr;
  long fsz;
  int bitpix;
  double dzero;
  int ploc;
  char ferrchar[80];
  int datatype;

  //_DEBUGHERE_("","");
  fitserr = 0;
  //_DEBUGHERE_("%s %d",path,strlen(path));

  fits_open_data(&fitsptr, path, READONLY, &fitserr);
  //_DEBUGHERE_("","");
  fits_get_errstatus(fitserr,ferrchar);
  //_DEBUGHERE_("","");
  testErrorRetVA(fitserr!=0,-1234,"Fits error (%d : %s) while opening %s",*err,__LINE__,NULL,fitserr,ferrchar,path);
  //_DEBUGHERE_("","");

  //_DEBUGHERE_("","");
  fitserr = 0;
  fits_get_img_type(fitsptr, &bitpix, &fitserr);
  fits_get_errstatus(fitserr,ferrchar);
  testErrorRetVA(fitserr!=0,-1234,"Fits error (%d : %s)  while reading %s",*err,__LINE__,NULL,fitserr,ferrchar,path);
  testErrorRetVA(bitpix!=typ ,-1234,"bad type for file  %s",*err,__LINE__,NULL,path);

  //_DEBUGHERE_("","");
  fitserr = 0;
  fits_get_img_size(fitsptr, 1, &fsz, &fitserr);
  fits_get_errstatus(fitserr,ferrchar);
  testErrorRetVA(fitserr!=0,-1234,"Fits error (%d : %s)  while reading %s",*err,__LINE__,NULL,fitserr,ferrchar,path);
  testErrorRetVA(fsz!=*sz && *sz>0 ,-1234,"bad size for file  %s (expected %d got %d)",*err,__LINE__,NULL,path,*sz,fsz);

  //_DEBUGHERE_("","");
  *sz = fsz;

  //_DEBUGHERE_("","");
  res = malloc_err((abs(typ)/8)*fsz,err);
  forwardError(*err,__LINE__,NULL);

  //_DEBUGHERE_("","");
  dzero = 0;
  fitserr = 0;
  datatype= TDOUBLE;
  if (typ>0) {
    datatype = TLONG;
  }
  //_DEBUGHERE_("","");
  fits_read_img(fitsptr, datatype, 1, fsz, &dzero, res, &ploc, &fitserr);
  fits_get_errstatus(fitserr,ferrchar);
  testErrorRetVA(fitserr!=0,-1234,"Fits error (%d : %s)  while reading %s",*err,__LINE__,NULL,fitserr,ferrchar,path);
  //_DEBUGHERE_("","");

  fitserr = 0;
  fits_close_file(fitsptr, &fitserr);
  fits_get_errstatus(fitserr,ferrchar);
  testErrorRetVA(fitserr!=0,-1234,"Fits error (%d : %s)  while closing %s",*err,__LINE__,NULL,fitserr,ferrchar,path);
  //_DEBUGHERE_("","");

  return res;
 
}

double* cldf_readfloatarray(cldf *df, char *key, int* sz, error **err) {
  double *res;
  int fitserr;
  fitsfile *fitsptr;
  long fsz;
  int bitpix;
  double dzero; 
  int ploc;
  char kp[256];
  int i;
  cldf *cdf;
  int nsz;

  nsz = -1;
  if (sz !=NULL) {
    nsz = *sz;
  }

#ifdef HDF5_COMPAT_MODE
  if (df->hid!=0) {
    herr_t hstat;
    int lk;
    char pk[1000];
    hsize_t ndum;
    H5T_class_t dum;
    size_t ddum;
    double *res;
    int j;
    int hk;

    lk = cldf_cut_key(key,pk);
    hk = cldf_has_attribute(df,pk,key,lk,err);
    forwardError(*err,__LINE__,NULL);
    
    if (hk == 1) {
      hstat = H5LTget_attribute_info( df->hid, pk,&(key[lk+1]), &ndum, &dum, &ddum);
      testErrorRetVA(hstat<0,-331234,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,key,df->root,hstat);
      testErrorRetVA((ndum!=nsz && nsz>0),-331234,"Bad size for %s in %s (got %d expected %d)",*err,__LINE__,NULL,key,df->root,ndum,nsz);
      res = malloc_err(sizeof(double)*ndum,err);
      forwardError(*err,__LINE__,NULL);
      hstat = H5LTget_attribute_double(df->hid, pk,&(key[lk+1]),res);
      testErrorRetVA(hstat<0,-331234,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,key,df->root,hstat);
      if (nsz<0 && sz !=NULL) {
        *sz = ndum;
      }
      return res;
    } else {
      hstat = H5LTget_dataset_info(df->hid, key, &ndum, &dum, &ddum);
      testErrorRetVA(hstat<0,-331234,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,key,df->root,hstat);
      testErrorRetVA((ndum!=nsz && nsz>0),-331234,"Bad size for %s in %s (got %d expected %d)",*err,__LINE__,NULL,key,df->root,ndum,nsz);
      res = malloc_err(sizeof(double)*ndum,err);
      forwardError(*err,__LINE__,NULL);
      hstat = H5LTread_dataset_double(df->hid, key,res);
      testErrorRetVA(hstat<0,-331234,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,key,df->root,hstat);
      if (nsz<0 && sz !=NULL) {
        *sz = ndum;
      }
      return res;
    }

  }
#endif

  cdf = cldf_tolast(df,key,kp,err);
  forwardError(*err,__LINE__,0);

  for(i=0;i<cdf->ndata;i++) {
    if (strcmp(kp,cdf->datakey[i])==0) {
      char pth[8192];
      
      if (cdf->root[0]=='\0') {
        sprintf(pth,"%s/%s",cdf->name,kp);  
      } else {
        sprintf(pth,"%s/%s/%s",cdf->root,cdf->name,kp);  
      }

      res = cldf_readanyfits(pth,&nsz,-64,err);
      forwardError(*err,__LINE__,NULL);

      if (sz!=NULL) {
        *sz = nsz;
      }
      return res;
    }
  }
  testErrorRetVA(1==1,-1234,"unknown element '%s'",*err,__LINE__,0,key);
  return 0;
}

int* cldf_readintarray(cldf *df, char *key, int* sz, error **err) {
  int *ires;
  long *res;
  int j;
  int nsz;

  nsz = -1;
  if (sz !=NULL) {
    nsz = *sz;
  }
  //_DEBUGHERE_("","");

  res = cldf_readlongarray(df,key,&nsz,err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");

  ires = malloc_err(sizeof(int)*(nsz),err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");
  
  for(j=0;j<nsz;j++) {
    ires[j] = res[j];
  }
  free(res);
  if (sz !=NULL) {
    *sz = nsz;
  }
  //_DEBUGHERE_("","");

  return ires;
}

long* cldf_readlongarray(cldf *df, char *key, int* sz, error **err) {
  long *res;
  int fitserr;
  fitsfile *fitsptr;
  long fsz;
  int bitpix;
  double dzero;
  int ploc;
  char kp[256];
  int i;
  cldf *cdf;
  int nsz;

  nsz = -1;
  if (sz!=NULL) {
    nsz = *sz;
  }
  //_DEBUGHERE_("","");

  #ifdef HDF5_COMPAT_MODE
  if (df->hid!=0) {
    herr_t hstat;
    int lk;
    char pk[1000];
    hsize_t ndum;
    H5T_class_t dum;
    size_t ddum;
    int *ires;
    long *res;
    int j;
    int hk;

    lk = cldf_cut_key(key,pk);
    hk = cldf_has_attribute(df,pk,key,lk,err);
    forwardError(*err,__LINE__,NULL);
    
    if (hk == 1) {
      hstat = H5LTget_attribute_info( df->hid, pk,&(key[lk+1]), &ndum, &dum, &ddum);
      testErrorRetVA(hstat<0,-331234,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,key,df->root,hstat);
      testErrorRetVA((ndum!=nsz && nsz>0),-331234,"Bad size for %s in %s (got %d expected %d)",*err,__LINE__,NULL,key,df->root,ndum,nsz);
      ires = malloc_err(sizeof(int)*ndum,err);
      res = malloc_err(sizeof(long)*ndum,err);
      forwardError(*err,__LINE__,NULL);
      hstat = H5LTget_attribute_int(df->hid, pk,&(key[lk+1]),ires);
      testErrorRetVA(hstat<0,-331234,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,key,df->root,hstat);
      if (nsz<0 && sz!=NULL) {
        *sz = ndum;
      }
      for(j=0;j<ndum;j++) {
        res[j] = ires[j];
      }
      free(ires);
      return res;
    } else {
      hstat = H5LTget_dataset_info(df->hid, key, &ndum, &dum, &ddum);
      testErrorRetVA(hstat<0,-331234,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,key,df->root,hstat);
      testErrorRetVA((ndum!=nsz && nsz>0),-331234,"Bad size for %s in %s (got %d expected %d)",*err,__LINE__,NULL,key,df->root,ndum,nsz);
      ires = malloc_err(sizeof(int)*ndum,err);
      res = malloc_err(sizeof(long)*ndum,err);
      forwardError(*err,__LINE__,NULL);
      hstat = H5LTread_dataset_int(df->hid, key,ires);
      testErrorRetVA(hstat<0,-331234,"cannot read %s in %s (got %d)",*err,__LINE__,NULL,key,df->root,hstat);
      if (nsz<0 & sz!=NULL) {
        *sz = ndum;
      }
      for(j=0;j<ndum;j++) {
        res[j] = ires[j];
      }
      free(ires);
      return res;
    }

  }
  #endif
  //_DEBUGHERE_("","");

  cdf = cldf_tolast(df,key,kp,err);
  forwardError(*err,__LINE__,0);
  //_DEBUGHERE_("","");

  for(i=0;i<cdf->ndata;i++) {
    if (strcmp(kp,cdf->datakey[i])==0) {
      char pth[8192];
      
      if (cdf->root[0]=='\0') {
        sprintf(pth,"%s/%s",cdf->name,kp);  
      } else {
        sprintf(pth,"%s/%s/%s",cdf->root,cdf->name,kp);  
      }

      res = cldf_readanyfits(pth,&nsz,64,err);
      forwardError(*err,__LINE__,NULL);

      if (sz !=NULL) {
        *sz = nsz;
      }

    return res;
    }
  }
  testErrorRetVA(1==1,-1234,"unknown element '%s'",*err,__LINE__,0,key);
  return 0;
}

void cldf_external(cldf* df, char *dirname, char* pwd, error **err) {
  char dirtmpl[2048*4];
  char *drn;
  char *fpix_data_name, fpix_data_n[10000];
  FILE *fpix_data;
  char command[4096*4];
  int status;

  testErrorRetVA(getcwd(pwd,4096)==NULL,-101010,"can't get cwd name (cause = '%s')",*err,__LINE__,,strerror(errno));
  
#ifdef HDF5_COMPAT_MODE
  herr_t hstat;
  hsize_t ndum;
  hid_t group_id;
  char *cur_lkl;
  if (df->hid!=0) {
    group_id = df->hid;
    cur_lkl = df->root;
    hstat = H5LTfind_dataset (group_id, "external_data");
    if (hstat==1) {
      char *data;
      
      // yes !
      sprintf(dirtmpl,"/tmp/clik_XXXXXX");
      drn = mkdtemp(dirtmpl);
      testErrorRetVA(drn==NULL,-100,"cannot create temporary dir (cause = '%s')",*err,__LINE__,,strerror(errno));

      // read tarfile from hdffile
      hstat = H5LTget_dataset_info( group_id, "external_data", &ndum, NULL,NULL);
      testErrorRetVA(hstat<0,-331234,"cannot read %s in %s (got %d)",*err,__LINE__,,"tardata",cur_lkl,hstat);
      data = malloc_err(sizeof(char)*ndum,err);
      forwardError(*err,__LINE__,);
      hstat = H5LTread_dataset(group_id,"external_data",H5T_NATIVE_UINT8,data);
      testErrorRetVA(hstat<0,-331234,"cannot read %s in %s (got %d)",*err,__LINE__,,"tardata",cur_lkl,hstat);

      fpix_data_name = malloc_err(sizeof(char)*4096,err);
      forwardError(*err,__LINE__,);

      // save to file !
      sprintf(fpix_data_name,"%s/data.tar",drn);
      fpix_data = fopen_err(fpix_data_name,"w",err);
      forwardError(*err,__LINE__,);
      testErrorRetVA(fwrite(data,1,ndum,fpix_data)<ndum,-100,"Cannot write to file %s",*err,__LINE__,,fpix_data_name);
      fclose(fpix_data);
      free(data);

      // change dir
      testErrorRetVA(chdir(drn)!=0,-100,"Cannot change dir to %s (cause = '%s')",*err,__LINE__,,drn,strerror(errno));

      // call tar to recreate the files  
      sprintf(command,"tar xf %s",fpix_data_name);
      status = system(command);
      testErrorRetVA(status!=0,-100,"cannot untar, command '%s' got status %d",*err,__LINE__,,command,status);
      sprintf(dirname,"%s",drn);
      free(fpix_data_name);
      return;    
    }  
  }
#endif
  fpix_data_name = cldf_readstr(df,"external_dir",NULL,err);
  forwardError(*err,__LINE__,);
  if (strlen(fpix_data_name)==1 && fpix_data_name[0]=='.') {
    free(fpix_data_name);
    fpix_data_name = malloc_err(sizeof(char)*(strlen(df->root)+strlen(df->name)+2+9+1),err);
    forwardError(*err,__LINE__,);
    sprintf(fpix_data_name,"%s/%s/_external",df->root,df->name);
  }
  testErrorRetVA(chdir(fpix_data_name)!=0,-100,"Cannot change dir to %s (cause = '%s')",*err,__LINE__,,fpix_data_name,strerror(errno));
  dirname[0]='\0';
  free(fpix_data_name);
}

void cldf_external_cleanup(char *dirname,char *pwd,error **err) {
  char command[4096*4];
  int status;
  
   // delete all files (like a macho !)
  testErrorRetVA(chdir(pwd)!=0,-100,"Cannot change dir to %s (cause = '%s')",*err,__LINE__,,pwd,strerror(errno));
  
  if (dirname[0]!='\0') {
    // remove files
    sprintf(command,"rm -rf %s",dirname); 
    //_DEBUGHERE_("%s",dirname);
    status = system(command);
    testErrorRetVA(status!=0,-100,"cannot delete files, command '%s' got status %d",*err,__LINE__,,command,status);    
  }
}


