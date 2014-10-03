/***************************************************************************
 *            scalefactor.c
 *
 *  Wed Nov 12 14:46:09 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:scalefactor
 * @title: Scale Factor
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "scalefactor.h"
#include "math/ncm_spline_cubic_notaknot.h"

#include <nvector/nvector_serial.h>

G_DEFINE_BOXED_TYPE (NcScaleFactor, nc_scale_factor, nc_scale_factor_copy, nc_scale_factor_free);

static gint dz_dt_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint dz_dt_J (_NCM_SUNDIALS_INT_TYPE N, realtype lambda, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static gint dz_dt_conformal_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data);
static gint dz_dt_conformal_J (_NCM_SUNDIALS_INT_TYPE N, realtype lambda, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/**
 * nc_scale_factor_new:
 * @ttype: a #NcScaleFactorTimeType
 * @zf: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcScaleFactor *
nc_scale_factor_new (NcScaleFactorTimeType ttype, gdouble zf)
{
  NcScaleFactor *a = g_slice_new (NcScaleFactor);
  a->ttype = ttype;
	a->a_t = ncm_spline_cubic_notaknot_new ();
	a->t_a = ncm_spline_cubic_notaknot_new ();
  a->spline_init = FALSE;
	a->dist = nc_distance_new (zf);
	a->ctrl = ncm_model_ctrl_new (NULL);
	a->zf = zf;

  a->cvode = CVodeCreate (CV_BDF, CV_NEWTON);
  NCM_CVODE_CHECK ((void *)a->cvode, "CVodeCreate", 0, NULL);
  a->cvode_malloc = FALSE;

  a->reltol = 1e-13;
  a->abstol = 1e-20;

  switch (ttype)
  {
    case NC_TIME_COSMIC:
      a->dz_dt_f = &dz_dt_f;
      a->dz_dt_J = &dz_dt_J;
      break;
    case NC_TIME_CONFORMAL:
      a->dz_dt_f = &dz_dt_conformal_f;
      a->dz_dt_J = &dz_dt_conformal_J;
      break;
    default:
      g_assert_not_reached ();
  }

  a->y = N_VNew_Serial(1);

  return a;
}

/**
 * nc_scale_factor_copy:
 * @a: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcScaleFactor *
nc_scale_factor_copy (NcScaleFactor *a)
{
  NcScaleFactor *new_a = nc_scale_factor_new (a->ttype, a->zf);
  return new_a;
}

/**
 * nc_scale_factor_free:
 * @a: FIXME
 *
 * FIXME
 */
void
nc_scale_factor_free (NcScaleFactor *a)
{
  if (a->a_t != NULL)
    ncm_spline_free (a->a_t);
  if (a->t_a != NULL)
    ncm_spline_free (a->t_a);
  CVodeFree (a->cvode);
  N_VDestroy (a->y);

	ncm_model_ctrl_free (a->ctrl);
	nc_distance_free (a->dist);

  g_slice_free (NcScaleFactor, a);
}

static void nc_scale_factor_init_cvode (NcScaleFactor *a, NcHICosmo *model);

/**
 * nc_scale_factor_init:
 * @a: FIXME
 * @zf: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
static void
_nc_scale_factor_init (NcScaleFactor *a, NcHICosmo *model)
{
  const gdouble Omega_k = nc_hicosmo_Omega_k (model);

  NV_Ith_S(a->y, 0) = a->zf;

  switch (a->ttype)
  {
    case NC_TIME_COSMIC:
      a->ti = nc_distance_cosmic_time (a->dist, model, a->zf);
      a->tf = nc_distance_cosmic_time (a->dist, model, 0.0);
      break;
    case NC_TIME_CONFORMAL:
      a->ti = nc_distance_conformal_time (a->dist, model, a->zf);
      a->tf = nc_distance_conformal_time (a->dist, model, 0.0);
      break;
    default:
      g_assert_not_reached ();
  }

  if (fabs (Omega_k) < 1e-14)
  {
    a->a0 = 1.0;
  }
  else
  {
    const gdouble RH_Mpc = (ncm_c_c () / (1.0e3 * nc_hicosmo_H0 (model)));
    a->a0 = RH_Mpc / sqrt (fabs (Omega_k));
  }

  nc_scale_factor_init_cvode (a, model);

  return;
}

static void
nc_scale_factor_init_cvode (NcScaleFactor *a, NcHICosmo *model)
{
  gint flag;

  if (a->cvode_malloc)
  {
    flag = CVodeReInit (a->cvode, a->ti, a->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }
  else
  {
    flag = CVodeInit (a->cvode, a->dz_dt_f, a->ti, a->y);
    NCM_CVODE_CHECK (&flag, "CVodeMalloc", 1, );
    a->cvode_malloc = TRUE;
  }

  flag = CVodeSStolerances (a->cvode, a->reltol, a->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

  flag = CVodeSetStopTime(a->cvode, a->tf);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  flag = CVodeSetUserData (a->cvode, model);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

  flag = CVodeSetMaxNumSteps(a->cvode, 500);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVDense(a->cvode, 1);
  NCM_CVODE_CHECK (&flag, "CVDense", 1, );

  flag = CVDlsSetDenseJacFn (a->cvode, a->dz_dt_J);
  NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );

  return;
}

/**
 * nc_scale_factor_z_t:
 * @a: FIXME
 * @t: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_scale_factor_z_t (NcScaleFactor *a, gdouble t)
{
  return -ncm_spline_eval (a->a_t, t);
}

/**
 * nc_scale_factor_t_z:
 * @a: FIXME
 * @z: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_scale_factor_t_z (NcScaleFactor *a, gdouble z)
{
  return ncm_spline_eval (a->t_a, -z);
}

/**
 * nc_scale_factor_t_x:
 * @a: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_scale_factor_t_x (NcScaleFactor *a, gdouble x)
{
  return ncm_spline_eval (a->t_a, -(x - 1.0));
}

/**
 * nc_scale_factor_a_t:
 * @a: FIXME
 * @t: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_scale_factor_a_t (NcScaleFactor *a, gdouble t)
{
  gdouble mz;
  mz = ncm_spline_eval (a->a_t, t);
  return a->a0 / (1.0 - mz);
}

/**
 * nc_scale_factor_calc_spline:
 * @a: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
static void
_nc_scale_factor_calc_spline (NcScaleFactor *a)
{
  GArray *x, *y;
  gdouble t, mzi;

	if (!ncm_spline_is_empty (a->a_t))
	{
		x = ncm_vector_get_array (a->a_t->xv);
		y = ncm_vector_get_array (a->a_t->yv);
		g_array_set_size (x, 0);
		g_array_set_size (y, 0);
	}
	else
	{
		x = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 1000);
		y = g_array_sized_new (FALSE, FALSE, sizeof(gdouble), 1000);
	}

  mzi = -a->zf;
  g_array_append_val (x, a->ti);
  g_array_append_val (y, mzi);

  while (1)
  {
    gint flag;

    flag = CVode(a->cvode, a->tf, a->y, &t, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode", 1, );
    mzi = -NV_Ith_S (a->y, 0);
    g_array_append_val (x, t);
    g_array_append_val (y, mzi);
    if (a->tf == t)
      break;
  }

	if (fabs (g_array_index (y, gdouble, y->len - 1)) < 1e-10)
	{
		g_array_index (y, gdouble, y->len - 1) = 0.0;
	}
	else
		g_error ("_nc_scale_factor_calc_spline today redshift must be zero not % 20.15g\n",
		         fabs (g_array_index (y, gdouble, y->len - 1)));

	ncm_spline_set_array (a->a_t, x, y, TRUE);
	ncm_spline_set_array (a->t_a, y, x, TRUE);

  g_array_unref (x);
  g_array_unref (y);

  return;
}

/**
 * nc_scale_factor_prepare:
 * @a: FIXME
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_scale_factor_prepare (NcScaleFactor *a, NcHICosmo *model)
{
	_nc_scale_factor_init (a, model);
	_nc_scale_factor_calc_spline (a);
	ncm_model_ctrl_update (a->ctrl, NCM_MODEL (model));
}

/**
 * nc_scale_factor_prepare_if_needed:
 * @a: FIXME
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
void
nc_scale_factor_prepare_if_needed (NcScaleFactor *a, NcHICosmo *model)
{
	if (ncm_model_ctrl_update (a->ctrl, NCM_MODEL (model)))
		nc_scale_factor_prepare (a, model);
}

static gint
dz_dt_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHICosmo *model = NC_HICOSMO (f_data);
  const gdouble z = NV_Ith_S (y, 0);
  const gdouble x = 1.0 + z;
  const gdouble E2 = nc_hicosmo_E2 (model, z);
  const gdouble E = sqrt(E2);

  NCM_UNUSED (t);
  
  NV_Ith_S (ydot, 0) = -x * E;
  return 0;
}

static gint
dz_dt_J (_NCM_SUNDIALS_INT_TYPE N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHICosmo *model = NC_HICOSMO (jac_data);
  const gdouble z = NV_Ith_S (y, 0);
  const gdouble x = 1.0 + z;
  const gdouble E2 = nc_hicosmo_E2 (model, z);
  const gdouble E = sqrt(E2);
  const gdouble dE2_dz = nc_hicosmo_dE2_dz (model, z);

  NCM_UNUSED (N);
  NCM_UNUSED (t);
  NCM_UNUSED (fy);
  NCM_UNUSED (tmp1);
  NCM_UNUSED (tmp2);
  NCM_UNUSED (tmp3);

  
  DENSE_ELEM (J, 0, 0) = - x * dE2_dz / (2.0 * E);

  return 0;
}

static gint
dz_dt_conformal_f (realtype t, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHICosmo *model = NC_HICOSMO (f_data);
  const gdouble z = NV_Ith_S (y, 0);
  const gdouble E2 = nc_hicosmo_E2 (model, z);
  const gdouble E = sqrt(E2);

  NCM_UNUSED (t);
  
  NV_Ith_S(ydot, 0) = -E;
  return 0;
}

static gint
dz_dt_conformal_J (_NCM_SUNDIALS_INT_TYPE N, realtype t, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHICosmo *model = NC_HICOSMO (jac_data);
  const gdouble z = NV_Ith_S(y,0);
  const gdouble E2 = nc_hicosmo_E2 (model, z);
  const gdouble E = sqrt(E2);
  const gdouble dE2_dz = nc_hicosmo_dE2_dz (model, z);

  NCM_UNUSED (N);
  NCM_UNUSED (t);
  NCM_UNUSED (fy);
  NCM_UNUSED (tmp1);
  NCM_UNUSED (tmp2);
  NCM_UNUSED (tmp3);
  
  DENSE_ELEM (J,0,0) = - dE2_dz / (2.0 * E);

  return 0;
}
