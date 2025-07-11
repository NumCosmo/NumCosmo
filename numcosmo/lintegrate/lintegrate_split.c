/* lintegrate_split.c
 *
 * Copyright (C) 2017 Matthew Pitkin
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

/* Here we provide wrappers to the lintegrate functions that can split the integrals
 * at a number of user defined points. This can be useful for, e.g., integrating a
 * sharply peaked function.
 */

#include "lintegrate.h"
#include "logadd.c"

/* wrapper for lintegrate_qag */
int lintegration_qag_split (const gsl_function *f, double *splitpts, size_t npts,
                            double epsabs, double epsrel, size_t limit,
                            int key, double * result, double * abserr) {
  const size_t nint = npts - 1; /* number of intervals */
  size_t i = 0;

  double sumres = -INFINITY, sumabserr = -INFINITY;

  for (i = 0; i < nint; i++){
    if (splitpts[i + 1] < splitpts[i]){
      GSL_ERROR ("Points are not in ascending order", GSL_EINVAL);
    }
  }

  /* perform integral on each interval and sum */
  for (i = 0; i < nint; i++){
    double tmpresult, tmperror;
    int retval;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (limit);

    retval = lintegration_qag(f, splitpts[i], splitpts[i+1], epsabs, epsrel, limit, key, w, &tmpresult, &tmperror);

    if ( retval != GSL_SUCCESS ){ return retval; }

    /* combined integrals from each interval */
    sumres = logaddexp(sumres, tmpresult);
    sumabserr = logaddexp(sumabserr, tmpresult);

    gsl_integration_workspace_free(w);
  }

  *result = sumres;
  *abserr = sumabserr;

  return GSL_SUCCESS;
}


/* wrapper for lintegrate_qng */
int lintegration_qng_split (const gsl_function *f, double *splitpts, size_t npts,
                            double epsabs, double epsrel, double * result, double * abserr,
                            size_t * nevals) {
  const size_t nint = npts - 1; /* number of intervals */
  size_t i = 0;

  double sumres = -INFINITY, sumabserr = -INFINITY;

  for (i = 0; i < nint; i++){
    if (splitpts[i + 1] < splitpts[i]){
      GSL_ERROR ("Points are not in ascending order", GSL_EINVAL);
    }
  }

  /* perform integral on each interval and sum */
  for (i = 0; i < nint; i++){
    double tmpresult, tmperror;
    int retval;
    size_t tmpnevals;

    retval = lintegration_qng(f, splitpts[i], splitpts[i+1], epsabs, epsrel, &tmpresult, &tmperror, &tmpnevals);

    if ( retval != GSL_SUCCESS ){ return retval; }

    /* combined integrals from each interval */
    sumres = logaddexp(sumres, tmpresult);
    sumabserr = logaddexp(sumabserr, tmpresult);
    *nevals += tmpnevals;
  }

  *result = sumres;
  *abserr = sumabserr;

  return GSL_SUCCESS;
}


/* wrapper for lintegrate_cquad */
int lintegration_cquad_split (const gsl_function * f, double *splitpts, size_t npts,
                              double epsabs, double epsrel, size_t wsints,
                              double *result, double *abserr, size_t * nevals) {
  const size_t nint = npts - 1; /* number of intervals */
  size_t i = 0;

  double sumres = -INFINITY, sumabserr = -INFINITY;

  for (i = 0; i < nint; i++){
    if (splitpts[i + 1] < splitpts[i]){
      GSL_ERROR ("Points are not in ascending order", GSL_EINVAL);
    }
  }

  /* perform integral on each interval and sum */
  for (i = 0; i < nint; i++){
    double tmpresult, tmperror;
    int retval;
    size_t tmpnevals;
    gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(wsints);

    retval = lintegration_cquad(f, splitpts[i], splitpts[i+1], epsabs, epsrel, w, &tmpresult, &tmperror, &tmpnevals);

    if ( retval != GSL_SUCCESS ){ return retval; }

    /* combined integrals from each interval */
    sumres = logaddexp(sumres, tmpresult);
    sumabserr = logaddexp(sumabserr, tmpresult);
    *nevals += tmpnevals;

    gsl_integration_cquad_workspace_free(w);
  }

  *result = sumres;
  *abserr = sumabserr;

  return GSL_SUCCESS;
}
