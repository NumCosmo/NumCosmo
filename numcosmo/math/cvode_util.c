/***************************************************************************
 *            cvode_util.c
 *
 *  Wed Nov 12 15:37:31 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION:cvode_util
 * @title: Sundials CVODE interface
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/cvode_util.h"

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

static N_Vector    _ncm_nvector_gsl_clone (N_Vector);
static N_Vector    _ncm_nvector_gsl_cloneempty (N_Vector);
static void        _ncm_nvector_gsl_destroy (N_Vector);
static void        _ncm_nvector_gsl_space (N_Vector, glong *, glong *);
static realtype *  _ncm_nvector_gsl_getarraypointer (N_Vector);
static void        _ncm_nvector_gsl_setarraypointer (realtype *, N_Vector);
static void        _ncm_nvector_gsl_linearsum (realtype, N_Vector, realtype, N_Vector, N_Vector); 
static void        _ncm_nvector_gsl_const (realtype, N_Vector);
static void        _ncm_nvector_gsl_prod (N_Vector, N_Vector, N_Vector);
static void        _ncm_nvector_gsl_div (N_Vector, N_Vector, N_Vector);
static void        _ncm_nvector_gsl_scale (realtype, N_Vector, N_Vector);
static void        _ncm_nvector_gsl_abs (N_Vector, N_Vector);
static void        _ncm_nvector_gsl_inv (N_Vector, N_Vector);
static void        _ncm_nvector_gsl_addconst (N_Vector, realtype, N_Vector);
static realtype    _ncm_nvector_gsl_dotprod (N_Vector, N_Vector);
static realtype    _ncm_nvector_gsl_maxnorm (N_Vector);
static realtype    _ncm_nvector_gsl_wrmsnorm (N_Vector, N_Vector);
static realtype    _ncm_nvector_gsl_wrmsnormmask (N_Vector, N_Vector, N_Vector);
static realtype    _ncm_nvector_gsl_min (N_Vector);
static realtype    _ncm_nvector_gsl_wl2norm (N_Vector, N_Vector);
static realtype    _ncm_nvector_gsl_l1norm (N_Vector);
static void        _ncm_nvector_gsl_compare (realtype, N_Vector, N_Vector);
static booleantype _ncm_nvector_gsl_invtest (N_Vector, N_Vector);
static booleantype _ncm_nvector_gsl_constrmask (N_Vector, N_Vector, N_Vector);
static realtype    _ncm_nvector_gsl_minquotient (N_Vector, N_Vector);

struct _generic_N_Vector_Ops _ncm_nvector_gsl[1] = {{
  &_ncm_nvector_gsl_clone,
  &_ncm_nvector_gsl_cloneempty,
  &_ncm_nvector_gsl_destroy,
  &_ncm_nvector_gsl_space,
  &_ncm_nvector_gsl_getarraypointer,
  &_ncm_nvector_gsl_setarraypointer,
  &_ncm_nvector_gsl_linearsum,
  &_ncm_nvector_gsl_const,
  &_ncm_nvector_gsl_prod,
  &_ncm_nvector_gsl_div,
  &_ncm_nvector_gsl_scale,
  &_ncm_nvector_gsl_abs,
  &_ncm_nvector_gsl_inv,
  &_ncm_nvector_gsl_addconst,
  &_ncm_nvector_gsl_dotprod,
  &_ncm_nvector_gsl_maxnorm,
  &_ncm_nvector_gsl_wrmsnorm,
  &_ncm_nvector_gsl_wrmsnormmask,
  &_ncm_nvector_gsl_min,
  &_ncm_nvector_gsl_wl2norm,
  &_ncm_nvector_gsl_l1norm,
  &_ncm_nvector_gsl_compare,
  &_ncm_nvector_gsl_invtest,
  &_ncm_nvector_gsl_constrmask,
  &_ncm_nvector_gsl_minquotient
}};

/**
 * ncm_cvode_util_nvector_new: (skip)
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
N_Vector 
ncm_cvode_util_nvector_new (guint n)
{
  N_Vector v = g_slice_new (struct _generic_N_Vector);
  v->ops = _ncm_nvector_gsl;
  v->content = gsl_vector_alloc (n);
  return v;
}

/**
 * ncm_cvode_util_check_flag:
 * @flagvalue: FIXME
 * @funcname: FIXME
 * @opt: FIXME 
 *
 * FIXME
 *
 * Returns: FIXME
*/
gboolean
ncm_cvode_util_check_flag (gpointer flagvalue, gchar *funcname, gint opt)
{
  gint *errflag;

  switch (opt)
  {
    case 0:
    {
      if (flagvalue == NULL) 
      {
        g_message ("\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return FALSE; 
      }
      break;
    }
    case 1:
    {   
      errflag = (int *) flagvalue;
      if (*errflag < 0) 
      {
        g_message ("\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
        return FALSE; 
      }
      break;
    }
    case 2:
    {
      if (flagvalue == NULL) 
      {
        g_message ("\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return FALSE; 
      }
      break;
    }
    default:
      g_assert_not_reached ();
  }  
  return TRUE;
}

/**
 * ncm_cvode_util_print_stats:
 * @cvode: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gboolean 
ncm_cvode_util_print_stats (gpointer cvode)
{
  glong nsteps, nfunceval, nlinsetups, njaceval, ndiffjaceval, nnonliniter, 
		nconvfail, nerrortests, nrooteval;
  gint flag, qcurorder, qlastorder;
	gdouble hinused, hlast, hcur, tcur;

	flag = CVodeGetIntegratorStats(cvode, &nsteps, &nfunceval,
	                               &nlinsetups, &nerrortests, &qlastorder, &qcurorder,
	                               &hinused, &hlast, &hcur, &tcur);
	ncm_cvode_util_check_flag(&flag, "CVodeGetIntegratorStats", 1);
	
	
  flag = CVodeGetNumNonlinSolvIters(cvode, &nnonliniter);
  ncm_cvode_util_check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode, &nconvfail);
  ncm_cvode_util_check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1); 

  flag = CVDlsGetNumJacEvals (cvode, &njaceval);
  ncm_cvode_util_check_flag (&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals (cvode, &ndiffjaceval);
  ncm_cvode_util_check_flag (&flag, "CVDlsGetNumRhsEvals", 1);  
    
  flag = CVodeGetNumGEvals(cvode, &nrooteval);
  ncm_cvode_util_check_flag(&flag, "CVodeGetNumGEvals", 1);

  flag = CVodeGetLastOrder(cvode, &qlastorder);
  ncm_cvode_util_check_flag(&flag, "CVodeGetLastOrder", 1);

  flag = CVodeGetCurrentOrder(cvode, &qcurorder);
  ncm_cvode_util_check_flag(&flag, "CVodeGetCurrentOrder", 1);

  g_message ("# Final Statistics:\n");
	g_message ("# nerrortests = %-6ld | %.5e %.5e %.5e %.5e\n", nerrortests, hinused, hlast, hcur, tcur);
  g_message ("# nsteps = %-6ld nfunceval  = %-6ld nlinsetups = %-6ld njaceval = %-6ld ndiffjaceval = %ld\n",
         nsteps, nfunceval, nlinsetups, njaceval, ndiffjaceval);
  g_message ("# nnonliniter = %-6ld nconvfail = %-6ld nerrortests = %-6ld nrooteval = %-6ld qcurorder = %-6d qlastorder = %-6d\n",
         nnonliniter, nconvfail, nerrortests, nrooteval, qcurorder, qlastorder);
  
  return TRUE;
}

#define NV2GSLV(v) ((gsl_vector *)(v)->content)

static N_Vector
_ncm_nvector_gsl_cloneempty (N_Vector src)
{
  N_Vector dest = g_slice_new (struct _generic_N_Vector);
  dest->content = NULL;
  dest->ops = src->ops;
  return dest;
}

static N_Vector
_ncm_nvector_gsl_clone (N_Vector src)
{
  N_Vector dest = _ncm_nvector_gsl_cloneempty (src);
  dest->content = gsl_vector_alloc (NV2GSLV(src)->size);
  return dest;
}

static void
_ncm_nvector_gsl_destroy (N_Vector v)
{
  gsl_vector_free (NV2GSLV(v));
  g_slice_free (struct _generic_N_Vector, v);
}

static void
_ncm_nvector_gsl_space (N_Vector v, glong *lrw, glong *lri)
{
  *lrw = NV2GSLV(v)->size;
  *lri = 2;
}

static realtype *
_ncm_nvector_gsl_getarraypointer (N_Vector v)
{
  return NV2GSLV(v)->data;
}

static void 
_ncm_nvector_gsl_setarraypointer (realtype *data, N_Vector v)
{
   NV2GSLV(v)->data = data;
}

static void
_ncm_nvector_gsl_linearsum (realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  register guint i;
//g_message ("%.15g %.15g %p %p %p\n", a, b, x, y ,z);
  if (b == 1 && y == z)
    gsl_blas_daxpy (a, NV2GSLV(x), NV2GSLV(z));
  else if (a == 1 && x == z)
    gsl_blas_daxpy (b, NV2GSLV(y), NV2GSLV(z));
  else
    for (i = 0; i < NV2GSLV(z)->size; i++)
      NV2GSLV(z)->data[i] = a * NV2GSLV(x)->data[i] + b * NV2GSLV(y)->data[i];  
}

static void
_ncm_nvector_gsl_const (realtype a, N_Vector v)
{
  gsl_vector_set_all (NV2GSLV(v), a);
}

static void
_ncm_nvector_gsl_prod (N_Vector x, N_Vector y, N_Vector z)
{
  register guint i;
  for (i = 0; i < NV2GSLV(z)->size; i++)
    NV2GSLV(z)->data[i] = NV2GSLV(x)->data[i] * NV2GSLV(y)->data[i];    
}

static void
_ncm_nvector_gsl_div (N_Vector x, N_Vector y, N_Vector z)
{
  register guint i;
  for (i = 0; i < NV2GSLV(z)->size; i++)
    NV2GSLV(z)->data[i] = NV2GSLV(x)->data[i] / NV2GSLV(y)->data[i];      
}

static void
_ncm_nvector_gsl_scale (realtype a, N_Vector x, N_Vector z)
{
  register guint i;
  if (x == z)
    gsl_blas_dscal (a, NV2GSLV(x));
  else
    for (i = 0; i < NV2GSLV(z)->size; i++)
      NV2GSLV(z)->data[i] = a * NV2GSLV(x)->data[i];        
}

static void
_ncm_nvector_gsl_abs (N_Vector x, N_Vector z)
{
  register guint i;
  for (i = 0; i < NV2GSLV(z)->size; i++)
    NV2GSLV(z)->data[i] = fabs(NV2GSLV(x)->data[i]);
}

static void
_ncm_nvector_gsl_inv (N_Vector x, N_Vector z)
{
  register guint i;
  for (i = 0; i < NV2GSLV(z)->size; i++)
    NV2GSLV(z)->data[i] = 1.0 / (NV2GSLV(x)->data[i]);
}

static void
_ncm_nvector_gsl_addconst (N_Vector x, realtype a, N_Vector z)
{
  register guint i;
  for (i = 0; i < NV2GSLV(z)->size; i++)
    NV2GSLV(z)->data[i] = a + (NV2GSLV(x)->data[i]);  
}

static realtype
_ncm_nvector_gsl_dotprod (N_Vector x, N_Vector y)
{
  gdouble result;
  gsl_blas_ddot (NV2GSLV(x), NV2GSLV(y), &result);
  return result;
}

static realtype
_ncm_nvector_gsl_maxnorm (N_Vector x)
{
  return NV2GSLV(x)->data[gsl_blas_idamax (NV2GSLV(x))];
}

static realtype
_ncm_nvector_gsl_wrmsnorm (N_Vector x, N_Vector w)
{
  gdouble result = 0.0;
  register guint i;
  for (i = 0; i < NV2GSLV(x)->size; i++)
  {
    const gdouble xw = (NV2GSLV(x)->data[i] * NV2GSLV(w)->data[i]);
    result += xw * xw;
  }
  return sqrt(result / (1.0 * NV2GSLV(x)->size));
}

static realtype
_ncm_nvector_gsl_wrmsnormmask (N_Vector x, N_Vector w, N_Vector id)
{
  gdouble result = 0.0;
  register guint i;
  for (i = 0; i < NV2GSLV(x)->size; i++)
  {
    const gdouble xw = (NV2GSLV(x)->data[i] * NV2GSLV(w)->data[i]);
    if (NV2GSLV(id)->data[i] > 0)
      result += xw * xw;
  }
  return sqrt(result / (1.0 * NV2GSLV(x)->size));
}

static realtype
_ncm_nvector_gsl_min (N_Vector x)
{
  return gsl_vector_min (NV2GSLV(x));
}

static realtype
_ncm_nvector_gsl_wl2norm (N_Vector x, N_Vector w)
{
  gdouble result = 0.0;
  register guint i;
  for (i = 0; i < NV2GSLV(x)->size; i++)
  {
    const gdouble xw = (NV2GSLV(x)->data[i] * NV2GSLV(w)->data[i]);
    result += xw * xw;
  }
  return sqrt(result);  
}

static realtype
_ncm_nvector_gsl_l1norm (N_Vector x)
{
  return gsl_blas_dasum (NV2GSLV(x));
}

static void
_ncm_nvector_gsl_compare (realtype c, N_Vector x, N_Vector z)
{
  register guint i;
  for (i = 0; i < NV2GSLV(z)->size; i++)
    NV2GSLV(z)->data[i] = fabs(NV2GSLV(x)->data[i]) > c ? 1.0 : 0.0;
}

static booleantype
_ncm_nvector_gsl_invtest (N_Vector x, N_Vector z)
{
  register guint i;
  for (i = 0; i < NV2GSLV(z)->size; i++)
  {
    NV2GSLV(z)->data[i] = 1.0 / (NV2GSLV(x)->data[i]);
    if (NV2GSLV(x)->data[i] == 0.0) return FALSE;
  }
  return TRUE;
}

static booleantype
_ncm_nvector_gsl_constrmask (N_Vector c, N_Vector x, N_Vector m)
{
  register guint i;
  gboolean ok = TRUE;
  for (i = 0; i < NV2GSLV(x)->size; i++)
  {
    gint t = NV2GSLV(c)->data[i];
    switch (t)
    {
      case 0:
        NV2GSLV(m)->data[i] = 0.0;
        break;
      case 2:
        NV2GSLV(m)->data[i] = NV2GSLV(x)->data[i] > 0.0 ? 0.0 : 1.0;
        break;
      case 1:
        NV2GSLV(m)->data[i] = NV2GSLV(x)->data[i] >= 0.0 ? 0.0 : 1.0;
        break;
      case -1:
        NV2GSLV(m)->data[i] = NV2GSLV(x)->data[i] <= 0.0 ? 0.0 : 1.0;
        break;
      case -2:
        NV2GSLV(m)->data[i] = NV2GSLV(x)->data[i] < 0.0 ? 0.0 : 1.0;
        break;
    }
    ok = ok && !NV2GSLV(m)->data[i];
  }
  return ok;
}

static realtype
_ncm_nvector_gsl_minquotient (N_Vector x, N_Vector y)
{
  register guint i;
  gdouble min = BIG_REAL;
  for (i = 0; i < NV2GSLV(x)->size; i++)
  {
    if (NV2GSLV(y)->data[i] == 0.0)
      continue;
    min = GSL_MIN (min, NV2GSLV(x)->data[i] / (NV2GSLV(y)->data[i]));
  }
  return min;
}
