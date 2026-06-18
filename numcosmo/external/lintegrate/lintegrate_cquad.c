/* lintegrate_cquad.c
 *
 * Copyright (C) 2017 Matthew Pitkin
 * Copyright (C) 2010 Pedro Gonnet
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

#include "cquad_const.c"
#include "logadd.c"
#include "logsub.c"

/* Compute the product of the fx with one of the inverse
    Vandermonde-like matrices. */

static void Vinvfx (const double *fx, double *c, const int d) {
  int i, j;

  switch (d)
    {
    case 0:
      for (i = 0; i <= 4; i++)
	{
	  double pluspart = -INFINITY, minuspart = -INFINITY;
          for (j = 0; j <= 4; j++){
	    if ( V1inv[i * 5 + j] >= 0. ) { pluspart = logaddexp(pluspart, log(V1inv[i * 5 + j]) + fx[j * 8]); }
	    else { minuspart = logaddexp(minuspart, log(-V1inv[i * 5 + j]) + fx[j * 8]); }
          }
          if ( minuspart > pluspart ) { c[i] = -INFINITY; }
          else { c[i] = logsubexp(pluspart, minuspart); }
	}
      break;
    case 1:
      for (i = 0; i <= 8; i++)
	{
	  double pluspart = -INFINITY, minuspart = -INFINITY;
	  for (j = 0; j <= 8; j++){
            if ( V2inv[i * 9 + j] >= 0. ) { pluspart = logaddexp(pluspart, log(V2inv[i * 9 + j]) + fx[j * 4]); }
	    else { minuspart = logaddexp(minuspart, log(-V2inv[i * 9 + j]) + fx[j * 4]); }
          }
          if ( minuspart > pluspart ) { c[i] = -INFINITY; }
          else { c[i] = logsubexp(pluspart, minuspart); }
	}
      break;
    case 2:
      for (i = 0; i <= 16; i++)
	{
	  double pluspart = -INFINITY, minuspart = -INFINITY;
	  for (j = 0; j <= 16; j++){
	    if ( V3inv[i * 17 + j] >= 0. ) { pluspart = logaddexp(pluspart, log(V3inv[i * 17 + j]) + fx[j * 2]); }
	    else { minuspart = logaddexp(minuspart, log(-V3inv[i * 17 + j]) + fx[j * 2]); }
          }
          if ( minuspart > pluspart ) { c[i] = -INFINITY; }
          else { c[i] = logsubexp(pluspart, minuspart); }
	}
      break;
    case 3:
      for (i = 0; i <= 32; i++)
	{
	  double pluspart = -INFINITY, minuspart = -INFINITY;
	  for (j = 0; j <= 32; j++){
            if ( V4inv[i * 33 + j] >= 0. ) { pluspart = logaddexp(pluspart, log(V4inv[i * 33 + j]) + fx[j]); }
	    else { minuspart = logaddexp(minuspart, log(-V4inv[i * 33 + j]) + fx[j]); }
          }
          if ( minuspart > pluspart ) { c[i] = -INFINITY; }
          else { c[i] = logsubexp(pluspart, minuspart); }
	}
      break;
    }

}


/* Downdate the interpolation given by the n coefficients c
    by removing the nodes with indices in nans. */

static void downdate (double *c, int n, int d, int *nans, int nnans) {

  static const int bidx[4] = { 0, 6, 16, 34 };
  double b_new[34], alpha;
  int i, j;

  for (i = 0; i <= n + 1; i++)
    b_new[i] = bee[bidx[d] + i];
  for (i = 0; i < nnans; i++)
    {
      b_new[n + 1] = b_new[n + 1] / Lalpha[n];
      b_new[n] = (b_new[n] + xi[nans[i]] * b_new[n + 1]) / Lalpha[n - 1];
      for (j = n - 1; j > 0; j--)
	b_new[j] =
	  (b_new[j] + xi[nans[i]] * b_new[j + 1] -
	   Lgamma[j + 1] * b_new[j + 2]) / Lalpha[j - 1];
      for (j = 0; j <= n; j++)
	b_new[j] = b_new[j + 1];
      //alpha = c[n] / b_new[n];
      alpha = c[n] - log(b_new[n]);
      for (j = 0; j < n; j++)
	c[j] = logsubexp(c[j], alpha + log(b_new[j]));
      c[n] = -INFINITY;
      n--;
    }
}


/* The actual integration routine.
    */

#ifdef HAVE_PYTHON_LINTEGRATE
int lintegration_cquad (pylintfunc f, void *funcdata, void *args, double a, double b,
                        double epsabs, double epsrel,
                        gsl_integration_cquad_workspace * ws,
                        double *result, double *abserr, size_t * nevals){
#else
int lintegration_cquad (const gsl_function * f, double a, double b,
                        double epsabs, double epsrel,
                        gsl_integration_cquad_workspace * ws,
                        double *result, double *abserr, size_t * nevals){
#endif
  /* Some constants that we will need. */
  static const int n[4] = { 4, 8, 16, 32 };
  static const int skip[4] = { 8, 4, 2, 1 };
  static const int idx[4] = { 0, 5, 14, 31 };
  static const double w = M_SQRT2 / 2;
  static const int ndiv_max = 20;

  /* Actual variables (as opposed to constants above). */
  double m, h, temp;
  double igral, err, igral_final, err_final;
  int nivals, neval = 0;
  int i, j, d, split, t;
  int nnans, nans[32];
  gsl_integration_cquad_ival *iv, *ivl, *ivr;
  double nc, ncdiff;

  /* Check the input arguments. */
  if (f == NULL)
    GSL_ERROR ("function pointer shouldn't be NULL", GSL_EINVAL);
  if (result == NULL)
    GSL_ERROR ("result pointer shouldn't be NULL", GSL_EINVAL);
  if (ws == NULL)
    GSL_ERROR ("workspace pointer shouldn't be NULL", GSL_EINVAL);

  /* Check for unreasonable accuracy demands */
  if (epsabs < 0.0 || epsrel < 0.0)
    GSL_ERROR ("tolerances may not be negative", GSL_EBADTOL);
  if (epsabs <= 0 && epsrel < GSL_DBL_EPSILON)
    GSL_ERROR ("unreasonable accuracy requirement", GSL_EBADTOL);

  /* Create the first interval. */
  iv = &(ws->ivals[0]);
  m = (a + b) / 2.;
  h = (b - a) / 2.;
  nnans = 0;
  for (i = 0; i <= n[3]; i++){
#ifdef HAVE_PYTHON_LINTEGRATE
    iv->fx[i] = f(m + xi[i] * h, funcdata, args);
#else
    iv->fx[i] = GSL_FN_EVAL (f, m + xi[i] * h);
#endif
    neval++;
    if (gsl_isnan(iv->fx[i]) || gsl_isinf(iv->fx[i]) == 1 ){ /* only count NaNs for actual NaNs and +infinity (-infinity is fine) */
      nans[nnans++] = i;
      iv->fx[i] = -INFINITY;
    }
  }
  Vinvfx (iv->fx, &(iv->c[idx[0]]), 0);
  Vinvfx (iv->fx, &(iv->c[idx[3]]), 3);
  Vinvfx (iv->fx, &(iv->c[idx[2]]), 2);
  for (i = 0; i < nnans; i++)
    iv->fx[nans[i]] = GSL_NAN;
  iv->a = a;
  iv->b = b;
  iv->depth = 3;
  iv->rdepth = 1;
  iv->ndiv = 0;
  iv->igral = iv->c[idx[3]] + log(2. * h * w);

  nc = -INFINITY;
  for (i = n[2] + 1; i <= n[3]; i++)
    {
      temp = iv->c[idx[3] + i];
      //nc += temp * temp;
      nc = logaddexp(nc, logaddexp(temp, temp));
    }
  ncdiff = nc;
  for (i = 0; i <= n[2]; i++){
    temp = LOGDIFF(iv->c[idx[2] + i], iv->c[idx[3] + i]);
    ncdiff = logaddexp(ncdiff, logaddexp(temp, temp));
    nc = logaddexp(nc, logaddexp(iv->c[idx[3] + i], iv->c[idx[3] + i]));
  }
  ncdiff *= 0.5;
  nc *= 0.5;

  iv->err = ncdiff + log(2. * h);
  if ( exp(ncdiff - nc) > 0.1 && iv->err < log(2. * h) + nc ) { iv->err = log(2. * h) + nc; }

  /* Initialize the heaps. */
  for (i = 0; i < (int)ws->size; i++) { ws->heap[i] = i; }

  /* Initialize some global values. */
  igral = iv->igral;
  err = iv->err;
  nivals = 1;
  igral_final = -INFINITY;
  err_final = -INFINITY;

  double lepsabs = log(epsabs);
  double lepsrel = log(epsrel);

  /* Main loop. */
  while (nivals > 0 &&
	 !(err <= igral + lepsrel || err <= lepsabs)
	 && !(err_final > igral + lepsrel && LOGDIFF(err, err_final) < igral + lepsrel)
	 && !(err_final > lepsabs && LOGDIFF(err, err_final) < lepsabs))
    {

      /* Put our finger on the interval with the largest error. */
      iv = &(ws->ivals[ws->heap[0]]);
      m = (iv->a + iv->b) / 2.;
      h = (iv->b - iv->a) / 2.;

      /* Should we try to increase the degree? */
      if (iv->depth < 3)
	{

	  /* Keep tabs on some variables. */
	  d = ++iv->depth;

	  /* Get the new (missing) function values */
	  for (i = skip[d]; i <= 32; i += 2 * skip[d])
	    {
#ifdef HAVE_PYTHON_LINTEGRATE
              iv->fx[i] = f(m + xi[i] * h, funcdata, args);
#else
              iv->fx[i] = GSL_FN_EVAL (f, m + xi[i] * h);
#endif
	      neval++;
	    }
	  nnans = 0;
	  for (i = 0; i <= 32; i += skip[d])
	    {
              if (gsl_isnan(iv->fx[i]) || gsl_isinf(iv->fx[i]) == 1 ){ /* only count NaNs for actual NaNs and +infinity (-infinity is fine) */
                nans[nnans++] = i;
                iv->fx[i] = -INFINITY;
              }
	    }

	  /* Compute the new coefficients. */
	  Vinvfx (iv->fx, &(iv->c[idx[d]]), d);

	  /* Downdate any NaNs. */
	  if (nnans > 0)
	    {
	      downdate (&(iv->c[idx[d]]), n[d], d, nans, nnans);
	      for (i = 0; i < nnans; i++)
		iv->fx[nans[i]] = GSL_NAN;
	    }

	  /* Compute the error estimate. */
	  nc = -INFINITY;
	  for (i = n[d - 1] + 1; i <= n[d]; i++) {
	      temp = iv->c[idx[d] + i];
              nc = logaddexp(nc, logaddexp(temp, temp));
	  }
	  ncdiff = nc;
	  for (i = 0; i <= n[d - 1]; i++) {
	      temp = LOGDIFF(iv->c[idx[d - 1] + i], iv->c[idx[d] + i]);
	      ncdiff = logaddexp(ncdiff, logaddexp(temp, temp));
              nc = logaddexp(nc, logaddexp(iv->c[idx[d] + i], iv->c[idx[d] + i]));
          }
	  ncdiff *= 0.5;
          nc *= 0.5;

          iv->err = ncdiff + log(2. * h);

	  /* Compute the local integral. */
          iv->igral = log(2. * h * w) + iv->c[idx[d]];

	  /* Split the interval prematurely? */
	  split = (exp(ncdiff - nc) > 0.1);

	}

      /* Maximum degree reached, just split. */
      else {
	  split = 1;
      }


      /* Should we drop this interval? */
      if ((m + h * xi[0]) >= (m + h * xi[1])
	  || (m + h * xi[31]) >= (m + h * xi[32])
	  || iv->err < iv->igral + log(GSL_DBL_EPSILON * 10))
	{

	  /* Keep this interval's contribution */
          err_final = logaddexp(err_final, iv->err);
          igral_final = logaddexp(igral_final, iv->igral);

	  /* Swap with the last element on the heap */
	  t = ws->heap[nivals - 1];
	  ws->heap[nivals - 1] = ws->heap[0];
	  ws->heap[0] = t;
	  nivals--;

	  /* Fix up the heap */
	  i = 0;
	  while (2 * i + 1 < nivals)
	    {

	      /* Get the kids */
	      j = 2 * i + 1;

	      /* If the j+1st entry exists and is larger than the jth,
	         use it instead. */
	      if (j + 1 < nivals
		  && ws->ivals[ws->heap[j + 1]].err >=
		  ws->ivals[ws->heap[j]].err)
		j++;

	      /* Do we need to move the ith entry up? */
	      if (ws->ivals[ws->heap[j]].err <= ws->ivals[ws->heap[i]].err)
		break;
	      else
		{
		  t = ws->heap[j];
		  ws->heap[j] = ws->heap[i];
		  ws->heap[i] = t;
		  i = j;
		}
	    }

	}

      /* Do we need to split this interval? */
      else if (split) {

	  /* Some values we will need often... */
	  d = iv->depth;

	  /* Generate the interval on the left */
	  ivl = &(ws->ivals[ws->heap[nivals++]]);
	  ivl->a = iv->a;
	  ivl->b = m;
	  ivl->depth = 0;
	  ivl->rdepth = iv->rdepth + 1;
	  ivl->fx[0] = iv->fx[0];
	  ivl->fx[32] = iv->fx[16];
	  for (i = skip[0]; i < 32; i += skip[0]) {
#ifdef HAVE_PYTHON_LINTEGRATE
              ivl->fx[i] = f((ivl->a + ivl->b) / 2. + xi[i] * h / 2., funcdata, args);
#else
              ivl->fx[i] = GSL_FN_EVAL (f, (ivl->a + ivl->b) / 2. + xi[i] * h / 2.);
#endif
	      neval++;
	  }
	  nnans = 0;
	  for (i = 0; i <= 32; i += skip[0]) {
	      if (gsl_isnan(ivl->fx[i]) || gsl_isinf(ivl->fx[i]) == 1 ){ /* only count NaNs for actual NaNs and +infinity (-infinity is fine) */
                nans[nnans++] = i;
                ivl->fx[i] = -INFINITY;
              }
          }
	  Vinvfx (ivl->fx, ivl->c, 0);
	  if (nnans > 0) {
	      downdate (ivl->c, n[0], 0, nans, nnans);
	      for (i = 0; i < nnans; i++)
		ivl->fx[nans[i]] = GSL_NAN;
	  }
	  for (i = 0; i <= n[d]; i++)  {
	      double pluspart = -INFINITY, minuspart = -INFINITY;
	      for (j = i; j <= n[d]; j++){
                if ( Tleft[i * 33 + j] >= 0. ){ pluspart = logaddexp(pluspart, log(Tleft[i * 33 + j]) + iv->c[idx[d] + j]); }
                else { minuspart = logaddexp(pluspart, log(-Tleft[i * 33 + j]) + iv->c[idx[d] + j]); }
		//ivl->c[idx[d] + i] += Tleft[i * 33 + j] * iv->c[idx[d] + j];
              }
              if ( minuspart > pluspart ) { ivl->c[idx[d] + i] = -INFINITY; }
              else { ivl->c[idx[d] + i]  = logsubexp(pluspart, minuspart); }
          }
	  ncdiff = -INFINITY;
	  for (i = 0; i <= n[0]; i++) {
              temp = LOGDIFF(ivl->c[i], ivl->c[idx[d] + i]);
              ncdiff = logaddexp(ncdiff, logaddexp(temp, temp));
	  }
	  for (i = n[0] + 1; i <= n[d]; i++) {
	      temp = ivl->c[idx[d] + i];
              ncdiff = logaddexp(ncdiff, logaddexp(temp, temp));
	    }
	  ncdiff *= 0.5;
          ivl->err = ncdiff + log(h);

	  /* Check for divergence. */
	  ivl->ndiv = iv->ndiv + (ivl->c[0] - iv->c[0] > M_LN2);
	  if (ivl->ndiv > ndiv_max && 2 * ivl->ndiv > ivl->rdepth) {
              /* need copysign(INFINITY, igral) */
	      *result = (igral >= 0) ? GSL_POSINF : GSL_NEGINF;
	      if (nevals != NULL)
		*nevals = neval;
	      return GSL_EDIVERGE;
          }

	  /* Compute the local integral. */
          ivl->igral = log(h * w) + ivl->c[0];

	  /* Generate the interval on the right */
	  ivr = &(ws->ivals[ws->heap[nivals++]]);
	  ivr->a = m;
	  ivr->b = iv->b;
	  ivr->depth = 0;
	  ivr->rdepth = iv->rdepth + 1;
	  ivr->fx[0] = iv->fx[16];
	  ivr->fx[32] = iv->fx[32];
	  for (i = skip[0]; i < 32; i += skip[0]) {
#ifdef HAVE_PYTHON_LINTEGRATE
              ivr->fx[i] = f((ivr->a + ivr->b) / 2. + xi[i] * h / 2., funcdata, args);
#else
              ivr->fx[i] = GSL_FN_EVAL (f, (ivr->a + ivr->b) / 2. + xi[i] * h / 2.);
#endif
	      neval++;
          }
	  nnans = 0;
	  for (i = 0; i <= 32; i += skip[0]) {
              if (gsl_isnan(ivr->fx[i]) || gsl_isinf(ivr->fx[i]) == 1 ){ /* only count NaNs for actual NaNs and +infinity (-infinity is fine) */
                nans[nnans++] = i;
                ivr->fx[i] = -INFINITY;
              }
          }
	  Vinvfx (ivr->fx, ivr->c, 0);
	  if (nnans > 0) {
	      downdate (ivr->c, n[0], 0, nans, nnans);
	      for (i = 0; i < nnans; i++)
		ivr->fx[nans[i]] = GSL_NAN;
          }
	  for (i = 0; i <= n[d]; i++) {
	      double pluspart = -INFINITY, minuspart = -INFINITY;
              for (j = i; j <= n[d]; j++){
                if ( Tright[i * 33 + j] >= 0. ){ pluspart = logaddexp(pluspart, log(Tright[i * 33 + j]) + iv->c[idx[d] + j]); }
                else { minuspart = logaddexp(pluspart, log(-Tright[i * 33 + j]) + iv->c[idx[d] + j]); }
              }
              if ( minuspart > pluspart ) { ivr->c[idx[d] + i] = -INFINITY; }
              else { ivr->c[idx[d] + i] = logsubexp(pluspart, minuspart); }
          }
	  ncdiff = -INFINITY;
	  for (i = 0; i <= n[0]; i++) {
	      temp = LOGDIFF(ivr->c[i], ivr->c[idx[d] + i]);
              ncdiff = logaddexp(ncdiff, logaddexp(temp, temp));
          }
	  for (i = n[0] + 1; i <= n[d]; i++) {
	      temp = ivr->c[idx[d] + i];
              ncdiff = logaddexp(ncdiff, logaddexp(temp, temp));
          }
	  ncdiff *= 0.5;
          ivr->err = ncdiff + log(h);

	  /* Check for divergence. */
	  ivr->ndiv = iv->ndiv + (ivr->c[0] - iv->c[0] > M_LN2);
	  if (ivr->ndiv > ndiv_max && 2 * ivr->ndiv > ivr->rdepth) {
              /* need copysign(INFINITY, igral) */
	      *result = (igral >= 0) ? GSL_POSINF : GSL_NEGINF;
	      if (nevals != NULL)
		*nevals = neval;
	      return GSL_EDIVERGE;
          }

	  /* Compute the local integral. */
          ivr->igral = log(h * w) + ivr->c[0];

	  /* Fix-up the heap: we now have one interval on top
	     that we don't need any more and two new, unsorted
	     ones at the bottom. */

	  /* Flip the last interval to the top of the heap and
	     sift down. */
	  t = ws->heap[nivals - 1];
	  ws->heap[nivals - 1] = ws->heap[0];
	  ws->heap[0] = t;
	  nivals--;

	  /* Sift this interval back down the heap. */
	  i = 0;
	  while (2 * i + 1 < nivals - 1) {
	      j = 2 * i + 1;
	      if (j + 1 < nivals - 1
		  && ws->ivals[ws->heap[j + 1]].err >=
		  ws->ivals[ws->heap[j]].err)
		j++;
	      if (ws->ivals[ws->heap[j]].err <= ws->ivals[ws->heap[i]].err)
		break;
	      else
		{
		  t = ws->heap[j];
		  ws->heap[j] = ws->heap[i];
		  ws->heap[i] = t;
		  i = j;
		}
          }

	  /* Now grab the last interval and sift it up the heap. */
	  i = nivals - 1;
	  while (i > 0) {
	      j = (i - 1) / 2;
	      if (ws->ivals[ws->heap[j]].err < ws->ivals[ws->heap[i]].err)
		{
		  t = ws->heap[j];
		  ws->heap[j] = ws->heap[i];
		  ws->heap[i] = t;
		  i = j;
		}
	      else
		break;
          }
	}

      /* Otherwise, just fix-up the heap. */
      else {
	  i = 0;
	  while (2 * i + 1 < nivals) {
	      j = 2 * i + 1;
	      if (j + 1 < nivals
		  && ws->ivals[ws->heap[j + 1]].err >=
		  ws->ivals[ws->heap[j]].err)
		j++;
	      if (ws->ivals[ws->heap[j]].err <= ws->ivals[ws->heap[i]].err)
		break;
	      else {
		  t = ws->heap[j];
		  ws->heap[j] = ws->heap[i];
		  ws->heap[i] = t;
		  i = j;
              }
          }
      }

      /* If the heap is about to overflow, remove the last two
         intervals. */
      while (nivals > (int)ws->size - 2) {
	  iv = &(ws->ivals[ws->heap[nivals - 1]]);

	  err_final = logaddexp(err_final, iv->err);
          igral_final = logaddexp(igral_final, iv->igral);
	  nivals--;
      }

      /* Collect the value of the integral and error. */
      igral = igral_final;
      err = err_final;
      for (i = 0; i < nivals; i++) {
          igral = logaddexp(igral, ws->ivals[ws->heap[i]].igral);
          err = logaddexp(err, ws->ivals[ws->heap[i]].err);
      }
    }

  /* Clean up and present the results. */
  *result = igral;
  if (abserr != NULL)
    *abserr = err;
  if (nevals != NULL)
    *nevals = neval;

  /* All is well that ends well. */
  return GSL_SUCCESS;

}
