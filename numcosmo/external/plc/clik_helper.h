#include "clik.h"
#include "lklbs.h"
#include <dlfcn.h>
#include <string.h>
#include <errno.h>
#include "cldf/cldf.h"
#if 0
#include "aplowly.h"
#include "fowly.h"
#endif

#ifndef _CLIK_HLP_
#define _CLIK_HLP_

#ifndef hdf5_base
#define hdf5_base  -1000
#endif

#define _dealwitherr error *lerr,**err; if(_err==NULL) {lerr=NULL;err=&lerr;} else {err=_err;}

#define _forwardError(A,B,C) if(_err!=NULL) {forwardError(A,B,C);} else {quitOnError(A,B,stderr);}
#define _testErrorRetVA(A,B,C,D,E,F,...) if(_err!=NULL) {testErrorRetVA(A,B,C,D,E,F,__VA_ARGS__);} else {testErrorExitVA(A,B,C,D,E,__VA_ARGS__);}


// get environ parameter
int clik_getenviron_integer(char* name, int sfg, error **err);
double clik_getenviron_real(char* name, double sfg, error **err);
char* clik_getenviron_string(char* name, char* sfg, error **err);
int clik_getenviron_numthread(char* name, int sfg, error **err);

// init lkls
cmblkl * clik_lklobject_init(cldf *df,error **err);
typedef cmblkl* clik_lkl_init_func(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err);
//typedef cmblkl* clik_addon_init_func(cmblkl* base, hid_t group_id, char* cur_lkl, error **err);

//void clik_external_data_init(char *pwd,char * dirname,hid_t group_id, char* cur_lkl,error **err);
//void clik_external_data_cleanup(char *pwd,char* dirname,error **err);

/*typedef {
  int nc;
  cmblkl** clkl;
  int* w;
  int cw; 
  int ic;
} wcmblkl;

double* wcmblkl_lkl(void* data, double *pars, error **err);
*/

int mtot(int mT,int mP,int *has_cl);

#endif
