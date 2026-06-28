#include "errorlist.h"
#include "io.h"
#include <math.h>

#ifndef CLIK_DIC
#define CLIK_DIC

#define hbcharsize 256

#ifndef NAN
#define M_NAN -1.2345e32
#define m_isnan(v) (v == M_NAN)
#else
#define M_NAN NAN
#define m_isnan(v) (isnan(v))
#endif

typedef struct {
  char **key;
  char **value;
  long *hash;
  double *dvalue;
  int nkey;
  int nmax;
} cdic;


typedef char cdic_tempchar[256];

char* cdic_newchar(char* fromchar,error **err);

long cdic_hash(char* key);

cdic* cdic_init(error **err);
void cdic_add_item(cdic* pf,int nit, char** key, char **value,error **err);
void cdic_free(void **ppf);
int cdic_key_index(cdic *pf, char *key, error **err);
void cdic_remove_item(cdic* pf, int index,error **err);
char* cdic_get(cdic* pf, char* key,char* safeguard, error **err);
long cdic_get_int(cdic *pf, char *key,long* safeguard, error **err);
double cdic_get_double(cdic *pf, char *key,double *safeguard, error **err);
char* cdic_sget(cdic* pf, char* key,char* safeguard, error **err);
long cdic_sget_int(cdic *pf, char *key,long* safeguard, error **err);
double cdic_sget_double(cdic *pf, char *key,double *safeguard, error **err);
char* cdic_set(cdic* pf, char* key,char* safeguard, error **err);
long cdic_set_int(cdic *pf, char *key,long safeguard, error **err);
double cdic_set_double(cdic *pf, char *key,double safeguard, error **err);
void cdic_dump(cdic *pf, FILE* strm, error **err);
char** cdic_get_keys(cdic *hbf, int *nkeys, error **err);
#endif