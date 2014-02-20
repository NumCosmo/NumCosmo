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
#include "build_cfg.h"

#include "lss/nc_growth_func.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/cvode_util.h"
#include "math/ncm_cfg.h"

#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <cvode/cvode_band.h>
#include <nvector/nvector_serial.h>

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
  NCM_UNUSED (gf);
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
  g_object_unref (gf);
}

/**
 * nc_growth_func_clear:
 * @gf: a #NcGrowthFunc.
 *
 * Atomically decrements the reference count of @gf by one. If the reference count drops to 0,
 * all memory allocated by @gf is released. Set pointer to NULL.
 *
 */
void
nc_growth_func_clear (NcGrowthFunc **gf)
{
  g_clear_object (gf);
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

  //printf ("res = % 20.15g y1 = % 20.15g y0 = % 20.15g E2 = % 20.15g dE2dz = % 20.15g xOm = % 20.15g 1/x = % 20.15g\n", NV_Ith_S (ydot, 1), NV_Ith_S (y, 1), NV_Ith_S (y, 0), E2, dE2dz, x * Omega_m, 1/x);
  return 0;
}

static gint
growth_J (_NCM_SUNDIALS_INT_TYPE N, realtype z, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2,
          N_Vector tmp3)
{
  NcmModel *model = NCM_MODEL (jac_data);
  gdouble E2 = nc_hicosmo_E2 (NC_HICOSMO (model), z);
  gdouble dE2dz = nc_hicosmo_dE2_dz (NC_HICOSMO (model), z);
  gdouble x = 1.0 + z;
  gdouble Omega_m = nc_hicosmo_Omega_m (NC_HICOSMO (model));

  NCM_UNUSED (N);
  NCM_UNUSED (y);
  NCM_UNUSED (fy);
  NCM_UNUSED (tmp1);
  NCM_UNUSED (tmp2);
  NCM_UNUSED (tmp3);
  
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
nc_growth_func_prepare (NcGrowthFunc *gf, NcHICosmo *model)
{
  gdouble zf, z;
  gint i = _NC_MAX_SPLINE_POINTS - 1;
  gint ii = 0;
  GArray *x_array, *y_array;

  if (gf->s != NULL)
  {
    x_array = ncm_vector_get_array (gf->s->xv);
    y_array = ncm_vector_get_array (gf->s->yv);
  }
  else
  {
    x_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NC_MAX_SPLINE_POINTS);
    y_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NC_MAX_SPLINE_POINTS);
  }

  g_array_set_size (x_array, _NC_MAX_SPLINE_POINTS);
  g_array_set_size (y_array, _NC_MAX_SPLINE_POINTS);

  NV_Ith_S (gf->yv, 0) = 1.0 / _NC_START_Z;
  NV_Ith_S (gf->yv, 1) = -1.0 / gsl_pow_2 (_NC_START_Z);

  if (gf->cvode == NULL)
  {
    gf->cvode = CVodeCreate (CV_ADAMS, CV_FUNCTIONAL);
    CVodeInit (gf->cvode, &growth_f, _NC_START_Z, gf->yv);
    if (FALSE)
    {
      CVDense (gf->cvode, 2);
      CVDlsSetDenseJacFn (gf->cvode, &growth_J);
    }
  }
  else
  {
    CVodeReInit (gf->cvode, _NC_START_Z, gf->yv);
  }

  CVodeSStolerances (gf->cvode, 1e-13, 0.0);
  CVodeSetMaxNumSteps (gf->cvode, 50000);
  CVodeSetUserData (gf->cvode, model);
  CVodeSetStopTime (gf->cvode, 0.0);

  g_array_index (y_array, gdouble, i) = 1.0 / _NC_START_Z;
  g_array_index (x_array, gdouble, i) = _NC_START_Z;
  i--;

  zf = _NC_START_Z;
  while (TRUE)
  {
    gint flag = CVode (gf->cvode, 0.0, gf->yv, &z, CV_ONE_STEP);
    if (!ncm_cvode_util_check_flag (&flag, "CVode", 1))
      g_error ("nc_growth_func_prepare: cannot integrate function.");

    if (zf == z && ii++ > 10)
      g_error ("nc_growth_func_prepare: cannot integrate function.");
    zf = z;
    
    g_array_index (y_array, gdouble, i) = NV_Ith_S (gf->yv, 0);
    g_array_index (x_array, gdouble, i) = zf;

    if (zf == 0.0)
      break;

    i--;
    if (i < 0)
    {
      /* i = i + 101; */
      for (i = 100; i >= 0; i--)
      {
        printf ("zf = %20.16g D = % 20.16g\n", g_array_index (x_array, gdouble, i), g_array_index (y_array, gdouble, i));
      }
      g_error ("Error: More than %d points to compute the growth spline.\n", _NC_MAX_SPLINE_POINTS);
    }

  }

  g_array_remove_range (y_array, 0, i);
  g_array_remove_range (x_array, 0, i);

  if (gf->s == NULL)
    gf->s = ncm_spline_cubic_notaknot_new ();

  {
    NcmVector *xv = ncm_vector_new_array (x_array);
    NcmVector *yv = ncm_vector_new_array (y_array);
    ncm_vector_scale (yv, 1.0 / ncm_vector_get (yv, 0));
    ncm_spline_set (gf->s, xv, yv, TRUE);
    ncm_vector_free (xv);
    ncm_vector_free (yv);
  }

  g_array_unref (x_array);
  g_array_unref (y_array);

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
/**
 * nc_growth_func_eval_both:
 * @gf: a #NcGrowthFunc.
 * @model: a #NcHICosmo.
 * @z: redshift.
 * @d: Growth function.
 * @f: Growth function derivative.
 *
 * FIXME
 *
 * Returns: FIXME
 */

static void
nc_growth_func_init (NcGrowthFunc *gf)
{
  gf->s       = NULL;
  gf->cvode   = NULL;
  gf->yv      = N_VNew_Serial (2);
  gf->zf      = 0.0;
  gf->ctrl    = ncm_model_ctrl_new (NULL);
}

static void
_nc_growth_func_dispose (GObject * object)
{
  NcGrowthFunc *gf = NC_GROWTH_FUNC (object);

  ncm_spline_clear (&gf->s);
  ncm_model_ctrl_clear (&gf->ctrl);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_growth_func_parent_class)->dispose (object);
}

static void
_nc_growth_func_finalize (GObject *object)
{
  NcGrowthFunc *gf = NC_GROWTH_FUNC (object);

  CVodeFree (&gf->cvode);
  N_VDestroy (gf->yv);

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

