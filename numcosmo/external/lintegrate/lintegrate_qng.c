/* Copyright (C) 2017 Matthew Pitkin
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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

#include "lintegrate.h"
#include "qng.h"
#include "err.c"
#include "logadd.c"
#include "logsub.c"

#ifdef HAVE_PYTHON_LINTEGRATE
int lintegration_qng (pylintfunc f, void *funcdata, void *args,
                      double a, double b,
                      double epsabs, double epsrel,
                      double * result, double * abserr, size_t * neval) {
#else
int lintegration_qng (const gsl_function *f,
                      double a, double b,
                      double epsabs, double epsrel,
                      double * result, double * abserr, size_t * neval) {
#endif
  double fv1[5], fv2[5], fv3[5], fv4[5];
  double savfun[21];  /* array of function values which have been computed */
  double res10, res21, res43, res87;    /* 10, 21, 43 and 87 point results */
  double result_kronrod, err ;
  double resabs; /* approximation to the integral of abs(f) */
  double resasc; /* approximation to the integral of abs(f-i/(b-a)) */

  const double half_length =  0.5 * (b - a);
  const double abs_half_length = fabs (half_length);
  const double center = 0.5 * (b + a);
#ifdef HAVE_PYTHON_LINTEGRATE
  const double f_center = f(center, funcdata, args);
#else
  const double f_center = GSL_FN_EVAL(f, center);
#endif

  const double lepsabs = log(epsabs), lepsrel = log(epsrel);

  int k ;

  if (epsabs <= 0 && (epsrel < 50 * GSL_DBL_EPSILON || epsrel < 0.5e-28)){
    * result = 0;
    * abserr = 0;
    * neval = 0;
    GSL_ERROR ("tolerance cannot be achieved with given epsabs and epsrel", GSL_EBADTOL);
  }

  /* Compute the integral using the 10- and 21-point formula. */

  res10 = -INFINITY;
  res21 = log(w21b[5]) + f_center;
  resabs = log(w21b[5]) + f_center;

  for (k = 0; k < 5; k++) {
    const double abscissa = half_length * x1[k];
#ifdef HAVE_PYTHON_LINTEGRATE
    const double fval1 = f(center + abscissa, funcdata, args);
    const double fval2 = f(center - abscissa, funcdata, args);
#else
    const double fval1 = GSL_FN_EVAL(f, center + abscissa);
    const double fval2 = GSL_FN_EVAL(f, center - abscissa);
#endif
    const double fval = logaddexp(fval1, fval2);
    res10 = logaddexp(res10, log(w10[k]) + fval);
    res21 = logaddexp(res21, log(w21a[k]) + fval);
    savfun[k] = fval;
    fv1[k] = fval1;
    fv2[k] = fval2;
  }

  for (k = 0; k < 5; k++) {
    const double abscissa = half_length * x2[k];
#ifdef HAVE_PYTHON_LINTEGRATE
    const double fval1 = f(center + abscissa, funcdata, args);
    const double fval2 = f(center - abscissa, funcdata, args);
#else
    const double fval1 = GSL_FN_EVAL(f, center + abscissa);
    const double fval2 = GSL_FN_EVAL(f, center - abscissa);
#endif
    const double fval = logaddexp(fval1, fval2);
    res21 = logaddexp(res21, log(w21b[k]) + fval);
    savfun[k + 5] = fval;
    fv3[k] = fval1;
    fv4[k] = fval2;
  }

  resabs = res21 + log(abs_half_length);

  const double mean = -M_LN2 + res21;

  resasc = log(w21b[5]) + LOGDIFF(f_center, mean);

  double resasctmp = 0.;
  for (k = 0; k < 5; k++) {
    resasctmp = logaddexp(log(w21a[k]) + logaddexp(LOGDIFF(fv1[k], mean), LOGDIFF(fv2[k], mean)),  log(w21b[k]) + logaddexp(LOGDIFF(fv3[k], mean), LOGDIFF(fv4[k], mean)));
    resasc = logaddexp(resasc, resasctmp);
  }
  resasc += log(abs_half_length);

  result_kronrod = res21 + log(half_length);

  err = rescale_error (LOGDIFF(res21, res10) + log(half_length), resabs, resasc) ;

  /*   test for convergence. */

  if (err < lepsabs || err < lepsrel + fabs (result_kronrod)){
    *result = result_kronrod;
    *abserr = err;
    *neval = 21;
    return GSL_SUCCESS;
  }

  /* compute the integral using the 43-point formula. */

  res43 = log(w43b[11]) + f_center;

  for (k = 0; k < 10; k++){
    res43 = logaddexp(res43, savfun[k] + log(w43a[k]));
  }

  for (k = 0; k < 11; k++){
    const double abscissa = half_length * x3[k];
#ifdef HAVE_PYTHON_LINTEGRATE
    const double fval = logaddexp(f(center + abscissa, funcdata, args), f(center - abscissa, funcdata, args));
#else
    const double fval = logaddexp(GSL_FN_EVAL(f, center + abscissa), GSL_FN_EVAL(f, center - abscissa));
#endif
    res43 = logaddexp(res43, fval + log(w43b[k]));
    savfun[k + 10] = fval;
  }

  /*  test for convergence */

  result_kronrod = res43 + log(half_length);
  err = rescale_error(LOGDIFF(res43, res21) + log(half_length), resabs, resasc);

  if (err < lepsabs || err < lepsrel + fabs (result_kronrod)) {
    *result = result_kronrod;
    *abserr = err;
    *neval = 43;
    return GSL_SUCCESS;
  }

  /* compute the integral using the 87-point formula. */

  res87 = log(w87b[22]) + f_center;

  for (k = 0; k < 21; k++){
    res87 = logaddexp(res87, savfun[k] + log(w87a[k]));
  }

  for (k = 0; k < 22; k++){
    const double abscissa = half_length * x4[k];
#ifdef HAVE_PYTHON_LINTEGRATE
    res87 = logaddexp(res87, log(w87b[k]) + logaddexp(f(center + abscissa, funcdata, args), f(center - abscissa, funcdata, args)));
#else
    res87 = logaddexp(res87, log(w87b[k]) + logaddexp(GSL_FN_EVAL(f, center + abscissa), GSL_FN_EVAL(f, center - abscissa)));
#endif
  }

  /*  test for convergence */

  result_kronrod = res87 + log(half_length);

  err = rescale_error (LOGDIFF(res87, res43) + log(half_length), resabs, resasc);

  if (err < lepsabs || err < lepsrel + fabs (result_kronrod)){
    *result = result_kronrod;
    *abserr = err;
    *neval = 87;
    return GSL_SUCCESS;
  }

  /* failed to converge */

  *result = result_kronrod;
  *abserr = err;
  *neval = 87;

  GSL_ERROR("failed to reach tolerance with highest-order rule", GSL_ETOL) ;
}

