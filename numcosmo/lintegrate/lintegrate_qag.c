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
#include "err.c"
#include "logadd.c"
#include "logsub.c"
#include "qkrules.c"

/* copied from GSL integration initialise.c http://git.savannah.gnu.org/cgit/gsl.git/tree/integration/initialise.c */
static inline
void initialise (gsl_integration_workspace * workspace, double a, double b);

static inline
void initialise (gsl_integration_workspace * workspace, double a, double b)
{
  workspace->size = 0;
  workspace->nrmax = 0;
  workspace->i = 0;
  workspace->alist[0] = a;
  workspace->blist[0] = b;
  workspace->rlist[0] = 0.0;
  workspace->elist[0] = 0.0;
  workspace->order[0] = 0;
  workspace->level[0] = 0;

  workspace->maximum_level = 0;
}


/* copied from GSL integration set_initial.c http://git.savannah.gnu.org/cgit/gsl.git/tree/integration/set_initial.c */
static inline
void set_initial_result (gsl_integration_workspace * workspace,
                         double result, double error);

static inline
void set_initial_result (gsl_integration_workspace * workspace,
                         double result, double error)
{
  workspace->size = 1;
  workspace->rlist[0] = result;
  workspace->elist[0] = error;
}


/* copied from GSL integration qpsrt.c http://git.savannah.gnu.org/cgit/gsl.git/tree/integration/qpsrt.c */
static inline void
qpsrt (gsl_integration_workspace * workspace);

static inline
void qpsrt (gsl_integration_workspace * workspace)
{
  const size_t last = workspace->size - 1;
  const size_t limit = workspace->limit;

  double * elist = workspace->elist;
  size_t * order = workspace->order;

  double errmax ;
  double errmin ;
  int i, k, top;

  size_t i_nrmax = workspace->nrmax;
  size_t i_maxerr = order[i_nrmax] ;

  /* Check whether the list contains more than two error estimates */

  if (last < 2)
    {
      order[0] = 0 ;
      order[1] = 1 ;
      workspace->i = i_maxerr ;
      return ;
    }

  errmax = elist[i_maxerr] ;

  /* This part of the routine is only executed if, due to a difficult
     integrand, subdivision increased the error estimate. In the normal
     case the insert procedure should start after the nrmax-th largest
     error estimate. */

  while (i_nrmax > 0 && errmax > elist[order[i_nrmax - 1]])
    {
      order[i_nrmax] = order[i_nrmax - 1] ;
      i_nrmax-- ;
    }

  /* Compute the number of elements in the list to be maintained in
     descending order. This number depends on the number of
     subdivisions still allowed. */

  if(last < (limit/2 + 2))
    {
      top = last ;
    }
  else
    {
      top = limit - last + 1;
    }

  /* Insert errmax by traversing the list top-down, starting
     comparison from the element elist(order(i_nrmax+1)). */

  i = i_nrmax + 1 ;

  /* The order of the tests in the following line is important to
     prevent a segmentation fault */

  while (i < top && errmax < elist[order[i]])
    {
      order[i-1] = order[i] ;
      i++ ;
    }

  order[i-1] = i_maxerr ;

  /* Insert errmin by traversing the list bottom-up */

  errmin = elist[last] ;

  k = top - 1 ;

  while (k > i - 2 && errmin >= elist[order[k]])
    {
      order[k+1] = order[k] ;
      k-- ;
    }

  order[k+1] = last ;

  /* Set i_max and e_max */

  i_maxerr = order[i_nrmax] ;

  workspace->i = i_maxerr ;
  workspace->nrmax = i_nrmax ;
}


/* update and retrieve functions copied from GSL integration util.c http://git.savannah.gnu.org/cgit/gsl.git/tree/integration/util.c */
static inline
void update (gsl_integration_workspace * workspace,
                 double a1, double b1, double area1, double error1,
                 double a2, double b2, double area2, double error2);

static inline void
retrieve (const gsl_integration_workspace * workspace,
          double * a, double * b, double * r, double * e);


static inline
void update (gsl_integration_workspace * workspace,
             double a1, double b1, double area1, double error1,
             double a2, double b2, double area2, double error2)
{
  double * alist = workspace->alist ;
  double * blist = workspace->blist ;
  double * rlist = workspace->rlist ;
  double * elist = workspace->elist ;
  size_t * level = workspace->level ;

  const size_t i_max = workspace->i ;
  const size_t i_new = workspace->size ;

  const size_t new_level = workspace->level[i_max] + 1;

  /* append the newly-created intervals to the list */

  if (error2 > error1)
    {
      alist[i_max] = a2;        /* blist[maxerr] is already == b2 */
      rlist[i_max] = area2;
      elist[i_max] = error2;
      level[i_max] = new_level;

      alist[i_new] = a1;
      blist[i_new] = b1;
      rlist[i_new] = area1;
      elist[i_new] = error1;
      level[i_new] = new_level;
    }
  else
    {
      blist[i_max] = b1;        /* alist[maxerr] is already == a1 */
      rlist[i_max] = area1;
      elist[i_max] = error1;
      level[i_max] = new_level;

      alist[i_new] = a2;
      blist[i_new] = b2;
      rlist[i_new] = area2;
      elist[i_new] = error2;
      level[i_new] = new_level;
    }

  workspace->size++;

  if (new_level > workspace->maximum_level)
    {
      workspace->maximum_level = new_level;
    }

  qpsrt (workspace) ;
}


static inline void
retrieve (const gsl_integration_workspace * workspace,
          double * a, double * b, double * r, double * e)
{
  const size_t i = workspace->i;
  double * alist = workspace->alist;
  double * blist = workspace->blist;
  double * rlist = workspace->rlist;
  double * elist = workspace->elist;

  *a = alist[i] ;
  *b = blist[i] ;
  *r = rlist[i] ;
  *e = elist[i] ;
}


/* modified sum_results function from GSL integration util.c http://git.savannah.gnu.org/cgit/gsl.git/tree/integration/util.c */
static inline double
log_sum_results (const gsl_integration_workspace * workspace);

static inline double
log_sum_results (const gsl_integration_workspace * workspace)
{
  const double * const rlist = workspace->rlist;
  const size_t n = workspace->size;

  size_t k;
  double result_sum = -INFINITY;

  for (k = 0; k < n; k++){
    result_sum = logaddexp(result_sum, rlist[k]);
  }

  return result_sum;
}


static inline int
subinterval_too_small (double a1, double a2, double b2);

static inline int
subinterval_too_small (double a1, double a2, double b2)
{
  const double e = GSL_DBL_EPSILON;
  const double u = GSL_DBL_MIN;

  double tmp = (1 + 100 * e) * (fabs (a2) + 1000 * u);

  int status = fabs (a1) <= tmp && fabs (b2) <= tmp;

  return status;
}


/* function to perform part of integral - heavily based on the GSL gsl_integration_qk function */
#ifdef HAVE_PYTHON_LINTEGRATE
void lintegration_qk (const int n, const double xgk[],
                      const double wg[], const double wgk[],
                      double fv1[], double fv2[],
                      pylintfunc f, void *funcdata, void *args, double a, double b,
                      double * result, double * abserr,
                      double * resabs, double * resasc){
#else
void lintegration_qk (const int n,
                      const double xgk[], const double wg[], const double wgk[],
                      double fv1[], double fv2[],
                      const gsl_function * f, double a, double b,
                      double *result, double *abserr,
                      double *resabs, double *resasc){
#endif
  const double center = 0.5 * (a + b);
  const double half_length = 0.5 * (b - a);
  const double abs_half_length = fabs (half_length);
#ifdef HAVE_PYTHON_LINTEGRATE
  const double f_center =  f(center, funcdata, args);
#else
  const double f_center = GSL_FN_EVAL (f, center);
#endif

  double result_gauss = -INFINITY;
  double result_kronrod = f_center + log(wgk[n - 1]);

  double result_abs = 0.;
  double result_asc = -INFINITY;
  double mean = 0, err = 0;

  int j;

  if (n % 2 == 0) {
    result_gauss = f_center + log(wg[n / 2 - 1]);
  }

  for (j = 0; j < (n - 1) / 2; j++){
    const int jtw = j * 2 + 1;  /* in original fortran j=1,2,3 jtw=2,4,6 */
    const double abscissa = half_length * xgk[jtw];
#ifdef HAVE_PYTHON_LINTEGRATE
    const double fval1 = f(center - abscissa, funcdata, args);
    const double fval2 = f(center + abscissa, funcdata, args);
#else
    const double fval1 = GSL_FN_EVAL (f, center - abscissa);
    const double fval2 = GSL_FN_EVAL (f, center + abscissa);
#endif
    const double fsum = logaddexp(fval1, fval2);
    fv1[jtw] = fval1;
    fv2[jtw] = fval2;
    result_gauss = logaddexp(result_gauss, log(wg[j]) + fsum);
    result_kronrod = logaddexp(result_kronrod, log(wgk[jtw]) + fsum);
  }

  for (j = 0; j < n / 2; j++){
    int jtwm1 = j * 2;
    const double abscissa = half_length * xgk[jtwm1];
#ifdef HAVE_PYTHON_LINTEGRATE
    const double fval1 = f(center - abscissa, funcdata, args);
    const double fval2 = f(center + abscissa, funcdata, args);
#else
    const double fval1 = GSL_FN_EVAL (f, center - abscissa);
    const double fval2 = GSL_FN_EVAL (f, center + abscissa);
#endif
    fv1[jtwm1] = fval1;
    fv2[jtwm1] = fval2;
    result_kronrod = logaddexp(result_kronrod, log(wgk[jtwm1]) + logaddexp(fval1, fval2));
  }

  result_abs = result_kronrod;

  mean = result_kronrod - M_LN2;

  result_asc = log(wgk[n - 1]) + LOGDIFF(f_center, mean);

  for (j = 0; j < n - 1; j++){
    result_asc = logaddexp(result_asc, log(wgk[j]) + logaddexp(LOGDIFF(fv1[j], mean), LOGDIFF(fv2[j], mean)));
  }

  /* scale by the width of the integration region */
  err = LOGDIFF(result_kronrod, result_gauss) + log(half_length);

  result_kronrod += log(half_length);
  result_abs += log(abs_half_length);
  result_asc += log(abs_half_length);

  *result = result_kronrod;
  *resabs = result_abs;
  *resasc = result_asc;
  *abserr = rescale_error (err, result_abs, result_asc);
}

#ifdef HAVE_PYTHON_LINTEGRATE
static int lqag (pylintfunc f, void *funcdata, void *args,
                 const double a, const double b,
                 const double epsabs, const double epsrel,
                 const size_t limit,
                 gsl_integration_workspace * workspace,
                 double * result, double * abserr,
                 py_integration_rule * q);
#else
static int lqag (const gsl_function *f,
                 const double a, const double b,
                 const double epsabs, const double epsrel,
                 const size_t limit,
                 gsl_integration_workspace * workspace,
                 double * result, double * abserr,
                 gsl_integration_rule * q);
#endif

#ifdef HAVE_PYTHON_LINTEGRATE
int lintegration_qag (pylintfunc f, void *funcdata, void *args,
                      double a, double b,
                      double epsabs, double epsrel, size_t limit,
                      int key,
                      gsl_integration_workspace * workspace,
                      double * result, double * abserr){
#else
int lintegration_qag (const gsl_function *f,
                      double a, double b,
                      double epsabs, double epsrel, size_t limit,
                      int key,
                      gsl_integration_workspace * workspace,
                      double * result, double * abserr){
#endif
  int status;
#ifdef HAVE_PYTHON_LINTEGRATE
  py_integration_rule * integration_rule = lintegration_qk15 ;
#else
  gsl_integration_rule * integration_rule = lintegration_qk15 ;
#endif

  if (key < GSL_INTEG_GAUSS15) {
    key = GSL_INTEG_GAUSS15;
  }
  else if (key > GSL_INTEG_GAUSS61) {
    key = GSL_INTEG_GAUSS61;
  }

  switch (key) {
    case GSL_INTEG_GAUSS15:
      integration_rule = lintegration_qk15;
      break;
    case GSL_INTEG_GAUSS21:
      integration_rule = lintegration_qk21;
      break;
    case GSL_INTEG_GAUSS31:
      integration_rule = lintegration_qk31;
      break;
    case GSL_INTEG_GAUSS41:
      integration_rule = lintegration_qk41;
      break;
    case GSL_INTEG_GAUSS51:
      integration_rule = lintegration_qk51;
      break;
    case GSL_INTEG_GAUSS61:
      integration_rule = lintegration_qk61;
      break;
    default:
      GSL_ERROR("value of key does specify a known integration rule", GSL_EINVAL) ;
    }

#ifdef HAVE_PYTHON_LINTEGRATE
  status = lqag(f, funcdata, args, a, b, epsabs, epsrel, limit,
                workspace,
                result, abserr,
                integration_rule);
#else
  status = lqag(f, a, b, epsabs, epsrel, limit,
                workspace,
                result, abserr,
                integration_rule);
#endif

  return status;
}

#ifdef HAVE_PYTHON_LINTEGRATE
static int lqag (pylintfunc f, void *funcdata, void *args,
                 const double a, const double b,
                 const double epsabs, const double epsrel,
                 const size_t limit,
                 gsl_integration_workspace * workspace,
                 double *result, double *abserr,
                 py_integration_rule * q) {
#else
static int lqag (const gsl_function * f,
                 const double a, const double b,
                 const double epsabs, const double epsrel,
                 const size_t limit,
                 gsl_integration_workspace * workspace,
                 double *result, double *abserr,
                 gsl_integration_rule * q) {
#endif
  double area, errsum;
  double result0, abserr0, resabs0, resasc0;
  double tolerance;
  size_t iteration = 0;
  int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;

  double round_off;
  const double lepsabs = log(epsabs), lepsrel = log(epsrel);

  /* Initialize results */

  initialise (workspace, a, b);

  *result = 0;
  *abserr = 0;

  if (limit > workspace->limit){
    GSL_ERROR ("iteration limit exceeds available workspace", GSL_EINVAL) ;
  }

  if (epsabs <= 0 && (epsrel < 50 * GSL_DBL_EPSILON || epsrel < 0.5e-28)){
      GSL_ERROR ("tolerance cannot be achieved with given epsabs and epsrel", GSL_EBADTOL);
  }

  /* perform the first integration */
#ifdef HAVE_PYTHON_LINTEGRATE
  q (f, funcdata, args, a, b, &result0, &abserr0, &resabs0, &resasc0);
#else
  q (f, a, b, &result0, &abserr0, &resabs0, &resasc0);
#endif

  set_initial_result (workspace, result0, abserr0);

  /* Test on accuracy */

  tolerance = GSL_MAX_DBL (lepsabs, lepsrel + fabs (result0));

  /* need IEEE rounding here to match original quadpack behavior */

  round_off = GSL_COERCE_DBL (log(50 * GSL_DBL_EPSILON) + resabs0);

  if (abserr0 <= round_off && abserr0 > tolerance){
    *result = result0;
    *abserr = abserr0;

    GSL_ERROR ("cannot reach tolerance because of roundoff error on first attempt", GSL_EROUND);
  }
  else if ((abserr0 <= tolerance && abserr0 != resasc0) || gsl_isinf(abserr0) == -1) {
    *result = result0;
    *abserr = abserr0;
    return GSL_SUCCESS;
  }
  else if (limit == 1) {
    *result = result0;
    *abserr = abserr0;

    GSL_ERROR ("a maximum of one iteration was insufficient", GSL_EMAXITER);
  }

  area = result0;
  errsum = abserr0;

  iteration = 1;

  do{
    double a1, b1, a2, b2;
    double a_i, b_i, r_i, e_i;
    double area1 = 0, area2 = 0, area12 = 0;
    double error1 = 0, error2 = 0, error12 = 0;
    double resasc1, resasc2;
    double resabs1, resabs2;

    /* Bisect the subinterval with the largest error estimate */

    retrieve (workspace, &a_i, &b_i, &r_i, &e_i);

    a1 = a_i;
    b1 = 0.5 * (a_i + b_i);
    a2 = b1;
    b2 = b_i;

#ifdef HAVE_PYTHON_LINTEGRATE
    q (f, funcdata, args, a1, b1, &area1, &error1, &resabs1, &resasc1);
    q (f, funcdata, args, a2, b2, &area2, &error2, &resabs2, &resasc2);
#else
    q (f, a1, b1, &area1, &error1, &resabs1, &resasc1);
    q (f, a2, b2, &area2, &error2, &resabs2, &resasc2);
#endif

    area12 = logaddexp(area1, area2);
    error12 = logaddexp(error1, error2);

    errsum = logsubexp(logaddexp(errsum, error12), e_i);
    area = logsubexp(logaddexp(area, area12), r_i);

    if (resasc1 != error1 && resasc2 != error2) {
      double delta = LOGDIFF(r_i, area12); /* stay in log-space when checking round off error */
      if ( delta <= log(1.0e-5) - fabs(area12) && error12 >= log(0.99) + e_i) {
        roundoff_type1++;
      }
      if (iteration >= 10 && error12 > e_i) {
        roundoff_type2++;
      }
    }

    tolerance = GSL_MAX_DBL (lepsabs, lepsrel + fabs (area));

    if (errsum > tolerance) {
      if (roundoff_type1 >= 6 || roundoff_type2 >= 20) {
        error_type = 2;   /* round off error */
      }

      /* set error flag in the case of bad integrand behaviour at
         a point of the integration range */

      if (subinterval_too_small (a1, a2, b2)) {
        error_type = 3;
      }
    }

    update (workspace, a1, b1, area1, error1, a2, b2, area2, error2);

    retrieve (workspace, &a_i, &b_i, &r_i, &e_i);

    iteration++;
  }
  while (iteration < limit && !error_type && errsum > tolerance);

  *result = log_sum_results (workspace);
  *abserr = errsum;

  if (errsum <= tolerance) {
    return GSL_SUCCESS;
  }
  else if (error_type == 2) {
   GSL_ERROR ("roundoff error prevents tolerance from being achieved", GSL_EROUND);
  }
  else if (error_type == 3) {
    GSL_ERROR ("bad integrand behavior found in the integration interval", GSL_ESING);
  }
  else if (iteration == limit) {
    GSL_ERROR ("maximum number of subdivisions reached", GSL_EMAXITER);
  }
  else {
    GSL_ERROR ("could not integrate function", GSL_EFAILED);
  }
}

