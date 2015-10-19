/*
 *  clik.h
 *  lowly_project
 *
 *  Created by Karim Benabed on 16/03/11.
 *  Copyright 2011 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#ifndef _CLIK_
#define _CLIK_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef CLIKSVNVERSION
#define CLIKSVNVERSION "NOPE"
#endif

#include "errorlist.h"

#define _pn_size 256
typedef char parname[_pn_size];

typedef void clik_object;
 
// in each of the following functions, if err is set to NULL the code exit as soon as it encounters an error (with a call to exit, bad for you if you're using mpi...)

// initialize the planck likelihood from an cldf file
clik_object* clik_init(char* hdffilepath, error **err);

// retrieve the list of cls as a list of flags 
// order is  TT EE BB TE TB EB
// for example for a likelihood acting on TT 
// has_cl = 1 0 0 0 0 0
void clik_get_has_cl(clik_object *clikid, int has_cl[6],error **err);


// retrieve the number of extra parameters and their names
// names is allocated by the function. It has to be cleaned after.
int clik_get_extra_parameter_names(clik_object* clikid, parname **names, error **err);

// retrieve the lmax for each power spectrum
// -1 --> no cl
// order is TT EE BB TE TB EB
// for example for a likelihood acting ont TT EE TE with same lmax 2000 
// lmax = 2000 2000 -1 2000 -1 -1 -1
void clik_get_lmax(clik_object *clikid, int lmax[6],error **err);

// compute a log likelyhood value
// cl_and_pars is order this way
// first the powerspectra from l=0 to l=lmax[cli] (included) in the order 
// TT EE BB TE TB EB. Only the one where lmax[cli]!=-1 have to be included
// then the extra parameters in the order given by clik_get_extra_parameters.
// The power spectra are in microK^2
// for example, for a likelihood acting on TT, EE and TE with 3 extra parameters 
// will expect an array ordered this way
// C_0^TT ... C_lmax[0]^TT C_0^EE ... C_lmax[1]^EE C_0^TE ... C_lmax[3]^T3 extrapar1 extrapar2 extrapar3
double clik_compute(clik_object* clikid, double* cl_and_pars,error **err);

// cleanup
void clik_cleanup(clik_object** pclikid);

//internal
void* _clik_dig(clik_object* clikid, error **err);
void* _clik_dig2(clik_object* clikid, error **err);

char* clik_get_version(clik_object *clikid,error **err);

// lensing
#include "clik_lensing.h"

#ifdef __cplusplus
}
#endif

#endif