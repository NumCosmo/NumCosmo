// this fakes hdf5 API with a bunch of directories along with fits and text files

#include "errorlist.h"
#include "io.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include "fitsio.h"
#include <dlfcn.h>
#include <string.h>
#include <errno.h>

#define MDB "_mdb"

#ifdef HDF5_COMPAT_MODE
#include "hdf5.h"
#endif

#ifndef CLDF_
#define CLDF_
typedef char kv[256];

typedef struct {
  char name[2048];
  char root[2048];
  int nmeta,ndata,nchild;
  kv* metakey;
  kv* metavalue;
  kv* metatype;
  kv* datakey;
  void **child;
#ifdef HDF5_COMPAT_MODE
  hid_t hid;
  int isgroup;
#endif
} cldf;

cldf * cldf_open(char *path, error **err);
cldf* cldf_openchild(cldf *df, char* ipath, error **err);
int cldf_haskey(cldf *df, char *key, error **err);
long cldf_readint(cldf *df, char *key, error **err);
double cldf_readfloat(cldf *df, char *key, error **err);
long cldf_readint_default(cldf *df, char *key, long def,error **err);
double cldf_readfloat_default(cldf *df, char *key,double def, error **err);
char* cldf_readstr(cldf *df, char *key, int* sz, error **err);
long* cldf_readlongarray(cldf *df, char *key, int* sz, error **err);
int* cldf_readintarray(cldf *df, char *key, int* sz, error **err);
double* cldf_readfloatarray(cldf *df, char *key, int* sz, error **err);
void cldf_external(cldf* df, char *dirname, char* pwd, error **err);
void cldf_external_cleanup(char* pwd,char *dirname,error **err);
void cldf_close(cldf** pdf);

#endif