/* Copyright (C) 2017 Matthew Pitkin
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#ifndef _LINTEGRATE_H
#define _LINTEGRATE_H

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#ifdef HAVE_PYTHON_LINTEGRATE
typedef double (*pylintfunc)(double x, void *funcdata, void *args);
#endif

#ifdef HAVE_PYTHON_LINTEGRATE
typedef void py_integration_rule (pylintfunc f, void *funcdata, void *args,
                                  double a, double b,
                                  double *result, double *abserr,
                                  double *defabs, double *resabs);

int lintegration_qag (pylintfunc f, void *funcdata, void *args,
                      double a, double b,
                      double epsabs, double epsrel, size_t limit,
                      int key,
                      gsl_integration_workspace * workspace,
                      double * result, double * abserr);

void lintegration_qk (const int n, const double xgk[],
                      const double wg[], const double wgk[],
                      double fv1[], double fv2[],
                      pylintfunc f, void *funcdata, void *args, double a, double b,
                      double * result, double * abserr,
                      double * resabs, double * resasc);

void lintegration_qk15 (pylintfunc f, void *funcdata, void *args, double a, double b,
                        double *result, double *abserr,
                        double *resabs, double *resasc);

void lintegration_qk21 (pylintfunc f, void *funcdata, void *args, double a, double b,
                        double *result, double *abserr,
                        double *resabs, double *resasc);

void lintegration_qk31 (pylintfunc f, void *funcdata, void *args, double a, double b,
                        double *result, double *abserr,
                        double *resabs, double *resasc);

void lintegration_qk41 (pylintfunc f, void *funcdata, void *args, double a, double b,
                        double *result, double *abserr,
                        double *resabs, double *resasc);

void lintegration_qk51 (pylintfunc f, void *funcdata, void *args, double a, double b,
                        double *result, double *abserr,
                        double *resabs, double *resasc);

void lintegration_qk61 (pylintfunc f, void *funcdata, void *args, double a, double b,
                        double *result, double *abserr,
                        double *resabs, double *resasc);
#else
int lintegration_qag (const gsl_function *f,
                      double a, double b,
                      double epsabs, double epsrel, size_t limit,
                      int key,
                      gsl_integration_workspace * workspace,
                      double * result, double * abserr);

int lintegration_qag_split (const gsl_function *f, double *splitpts, size_t npts,
                            double epsabs, double epsrel, size_t limit,
                            int key, double * result, double * abserr);

void lintegration_qk (const int n, const double xgk[],
                      const double wg[], const double wgk[],
                      double fv1[], double fv2[],
                      const gsl_function *f, double a, double b,
                      double * result, double * abserr,
                      double * resabs, double * resasc);

void lintegration_qk15 (const gsl_function * f, double a, double b,
                        double *result, double *abserr,
                        double *resabs, double *resasc);

void lintegration_qk21 (const gsl_function * f, double a, double b,
                        double *result, double *abserr,
                        double *resabs, double *resasc);

void lintegration_qk31 (const gsl_function * f, double a, double b,
                        double *result, double *abserr,
                        double *resabs, double *resasc);

void lintegration_qk41 (const gsl_function * f, double a, double b,
                        double *result, double *abserr,
                        double *resabs, double *resasc);

void lintegration_qk51 (const gsl_function * f, double a, double b,
                        double *result, double *abserr,
                        double *resabs, double *resasc);

void lintegration_qk61 (const gsl_function * f, double a, double b,
                        double *result, double *abserr,
                        double *resabs, double *resasc);
#endif

#ifdef HAVE_PYTHON_LINTEGRATE
int lintegration_qng (pylintfunc f, void *funcdata, void *args,
                      double a, double b,
                      double epsabs, double epsrel,
                      double * result, double * abserr, size_t * neval);
#else
int lintegration_qng (const gsl_function *f,
                      double a, double b,
                      double epsabs, double epsrel,
                      double * result, double * abserr, size_t * neval);

int lintegration_qng_split (const gsl_function *f, double *splitpts, size_t npts,
                            double epsabs, double epsrel,
                            double * result, double * abserr, size_t * neval);
#endif

#ifdef HAVE_PYTHON_LINTEGRATE
int lintegration_cquad (pylintfunc f, void *funcdata, void *args, double a, double b,
                        double epsabs, double epsrel,
                        gsl_integration_cquad_workspace * ws,
                        double *result, double *abserr, size_t * nevals);
#else
int lintegration_cquad (const gsl_function * f, double a, double b,
                        double epsabs, double epsrel,
                        gsl_integration_cquad_workspace * ws,
                        double *result, double *abserr, size_t * nevals);

int lintegration_cquad_split (const gsl_function * f, double *splitpts, size_t npts,
                              double epsabs, double epsrel, size_t wsints,
                              double *result, double *abserr, size_t * nevals);
#endif

#if HAVE_EXTENDED_PRECISION_REGISTERS
#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#else
#define GSL_COERCE_DBL(x) (x)
#endif

#endif /* _LINTEGRATE_H */

