/***************************************************************************
 *            nc_growth_func.c
 *
 *  Tue Apr  6 01:12:58 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
 *
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
 * SECTION:nc_growth_func
 * @title: Perturbations Growth Function
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <cvodes/cvodes_dense.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

G_DEFINE_TYPE (NcGrowthFunc, nc_growth_func, G_TYPE_OBJECT);

/**
 * nc_growth_func_new:
 *
 * This function allocates memory for a new #NcGrowthFunc object.
 *
 * Returns: A new #NcGrowthFunc.
*/
NcGrowthFunc *
nc_growth_func_new (void)
{
  return g_object_new (NC_TYPE_GROWTH_FUNC, NULL);
}

/**
 * nc_growth_func_copy:
 * @gf: a #NcGrowthFunc.
 *
 * This function duplicates @gf.
 *
 * Returns: (transfer full): A #NcGrowthFunc.
*/
NcGrowthFunc *
nc_growth_func_copy (NcGrowthFunc *gf)
{
  return nc_growth_func_new ();
}


/**
 * nc_growth_func_free:
 * @gf: a #NcGrowthFunc.
 *
 * Atomically decrements the reference count of @gf by one. If the reference count drops to 0,
 * all memory allocated by @gf is released.
 *
*/
void
nc_growth_func_free (NcGrowthFunc *gf)
{
  g_clear_object (&gf);
}

static gint
growth_f (realtype z, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcmModel *model = NCM_MODEL (f_data);
  gdouble E2 = nc_hicosmo_E2 (NC_HICOSMO (model), z);
  gdouble dE2dz = nc_hicosmo_dE2_dz (NC_HICOSMO (model), z);
  gdouble x = 1.0 + z;
  gdouble Omega_m = nc_hicosmo_Omega_m (NC_HICOSMO (model));

  NV_Ith_S (ydot, 0) = NV_Ith_S (y, 1);
  NV_Ith_S (ydot, 1) =
    (1.0 / x - dE2dz / (2.0 * E2)) * NV_Ith_S (y, 1) + 3.0 * Omega_m * x * NV_Ith_S (y, 0) / (2.0 * E2);

  return 0;
}

static int
growth_J (int N, realtype z, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2,
	  N_Vector tmp3)
{
  NcmModel *model = NCM_MODEL (jac_data);
  gdouble E2 = nc_hicosmo_E2 (NC_HICOSMO (model), z);
  gdouble dE2dz = nc_hicosmo_dE2_dz (NC_HICOSMO (model), z);
  gdouble x = 1.0 + z;
  gdouble Omega_m = nc_hicosmo_Omega_m (NC_HICOSMO (model));

  DENSE_ELEM (J, 0, 0) = 0.0;
  DENSE_ELEM (J, 0, 1) = 1.0;

  DENSE_ELEM (J, 1, 0) = 3.0 * Omega_m * x / (2.0 * E2);
  DENSE_ELEM (J, 1, 1) = (1.0 / x - dE2dz / (2.0 * E2));

  return 0;
}

#define _NC_START_Z (1000.0)
#define _NC_MAX_SPLINE_POINTS 50000

/**
 * nc_growth_func_prepare:
 * @gf: a #NcGrowthFunc.
 * @model: a #NcHICosmo.
 *
 * FIXME
 *
*/
void
nc_growth_func_prepare (NcGrowthFunc * gf, NcHICosmo * model)
{
  gpointer cvode;
  N_Vector yv;
  gdouble zf;
  gint i = _NC_MAX_SPLINE_POINTS - 1;
  GArray *x_array, *y_array;

  if (gf->s != NULL)
  {
    x_array = gf->s->xv->a;
    y_array = gf->s->yv->a;
    /* FIXME LEAK!! */
    ncm_spline_free (gf->s);
  }
  else
  {
    x_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NC_MAX_SPLINE_POINTS);
    y_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NC_MAX_SPLINE_POINTS);
  }

  g_array_set_size (x_array, _NC_MAX_SPLINE_POINTS);
  g_array_set_size (y_array, _NC_MAX_SPLINE_POINTS);

  yv = N_VNew_Serial (2);

  NV_Ith_S (yv, 0) = 1.0 / _NC_START_Z;
  NV_Ith_S (yv, 1) = -1.0 / gsl_pow_2 (_NC_START_Z);

//  cvode = CVodeCreate(CV_BDF, CV_NEWTON);
  cvode = CVodeCreate (CV_ADAMS, CV_FUNCTIONAL);
  CVodeInit (cvode, &growth_f, _NC_START_Z, yv);
  CVodeSStolerances (cvode, 1e-13, 0.0);
  CVodeSetUserData (cvode, model);
  CVodeSetMaxNumSteps (cvode, 50000);
  CVDense (cvode, 2);
  CVDlsSetDenseJacFn (cvode, &growth_J);

  CVodeSetStopTime (cvode, 0.0);

  g_array_index (y_array, gdouble, i) = 1.0 / _NC_START_Z;
  g_array_index (x_array, gdouble, i) = _NC_START_Z;
  i--;

  while (1)
  {
    CVode (cvode, 0.0, yv, &zf, CV_ONE_STEP);
    g_array_index (y_array, gdouble, i) = NV_Ith_S (yv, 0);
    g_array_index (x_array, gdouble, i) = zf;
    if (zf == 0.0)
      break;

    i--;
    if (i < 0)
      g_error ("Error: More than %d points to compute the growth spline.\n", _NC_MAX_SPLINE_POINTS);
  }
  gf->one_gf0 = 1.0 / g_array_index (y_array, gdouble, i);

  CVodeFree (&cvode);
  N_VDestroy (yv);
  g_array_remove_range (y_array, 0, i);
  g_array_remove_range (x_array, 0, i);

  gf->s = ncm_spline_cubic_notaknot_new ();
  ncm_spline_set_array (gf->s, x_array, y_array, TRUE);

  return;
}

/**
 * nc_growth_func_eval:
 * @gf: a #NcGrowthFunc.
 * @model: a #NcHICosmo.
 * @z: redshift.
 *
 * FIXME
 *
 * Returns: The normalized groth function at @z.
*/
gdouble
nc_growth_func_eval (NcGrowthFunc *gf, NcHICosmo *model, gdouble z)
{
  if (z > _NC_START_Z || z < 0.0)
    g_error ("Growth spline prepared for (%g, %g), %g is not included.\n", 0.0, _NC_START_Z, z);
  if (ncm_model_ctrl_update (gf->ctrl, NCM_MODEL(model)))
    nc_growth_func_prepare (gf, model);
  return ncm_spline_eval (gf->s, z) * gf->one_gf0;	/* Spline normalized */
}

/**
 * nc_growth_func_eval_deriv:
 * @gf: a #NcGrowthFunc.
 * @model: a #NcHICosmo.
 * @z: redshift.
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_growth_func_eval_deriv (NcGrowthFunc *gf, NcHICosmo *model, gdouble z)
{
  if (z > _NC_START_Z || z < 0.0)
    g_error ("Growth spline prepared for (%g, %g), %g is not included.\n", 0.0, _NC_START_Z, z);
  if (ncm_model_ctrl_update (gf->ctrl, NCM_MODEL(model)))
    nc_growth_func_prepare (gf, model);
  return ncm_spline_eval_deriv (gf->s, z) * gf->one_gf0;	/* Spline normalized */
}

static void
nc_growth_func_init (NcGrowthFunc *gf)
{
  gf->s = NULL;
  gf->one_gf0 = 0.0;
  gf->zf = 0.0;
  gf->ctrl = ncm_model_ctrl_new (NULL);
}

static void
_nc_growth_func_dispose (GObject * object)
{
  NcGrowthFunc *gf = NC_GROWTH_FUNC (object);
  ncm_model_ctrl_free (gf->ctrl);
  G_OBJECT_CLASS (nc_growth_func_parent_class)->dispose (object);
}

static void
_nc_growth_func_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_growth_func_parent_class)->finalize (object);
}

static void
nc_growth_func_class_init (NcGrowthFuncClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->dispose = _nc_growth_func_dispose;
  object_class->finalize = _nc_growth_func_finalize;
}

