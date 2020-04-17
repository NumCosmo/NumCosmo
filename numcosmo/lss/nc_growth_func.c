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
 * @title: NcGrowthFunc
 * @short_description: Growth function of linear perturbations.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_growth_func.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_cfg.h"


#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */

#ifndef NUMCOSMO_GIR_SCAN
#include <cvode/cvode.h>

#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#define SUN_DENSE_ACCESS SM_ELEMENT_D

#include <nvector/nvector_serial.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGrowthFuncPrivate
{
  gpointer cvode;
  N_Vector yv;
  SUNMatrix A;
  SUNLinearSolver LS;
  gdouble zf;
  NcmModelCtrl *ctrl_cosmo;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGrowthFunc, nc_growth_func, G_TYPE_OBJECT);

static void
nc_growth_func_init (NcGrowthFunc *gf)
{
  NcGrowthFuncPrivate * const self = gf->priv = nc_growth_func_get_instance_private (gf);
  self->cvode      = NULL;
  self->yv         = N_VNew_Serial (3);
  self->A          = SUNDenseMatrix (3, 3);
  self->LS         = SUNDenseLinearSolver (self->yv, self->A);
  self->zf         = 0.0;
  self->ctrl_cosmo = ncm_model_ctrl_new (NULL);
  
  gf->s            = NULL;
  gf->Da0          = 0.0;

  NCM_CVODE_CHECK ((gpointer)self->A, "SUNDenseMatrix", 0, );
  NCM_CVODE_CHECK ((gpointer)self->LS, "SUNDenseLinearSolver", 0, );
}

static void
_nc_growth_func_dispose (GObject * object)
{
  NcGrowthFunc *gf = NC_GROWTH_FUNC (object);
  NcGrowthFuncPrivate * const self = gf->priv;

  ncm_spline_clear (&gf->s);
  ncm_model_ctrl_clear (&self->ctrl_cosmo);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_growth_func_parent_class)->dispose (object);
}

static void
_nc_growth_func_finalize (GObject *object)
{
  NcGrowthFunc *gf = NC_GROWTH_FUNC (object);
  NcGrowthFuncPrivate * const self = gf->priv;

  CVodeFree (&self->cvode);
  N_VDestroy (self->yv);

  if (self->A != NULL)
  {
    SUNMatDestroy (self->A);
    self->A = NULL;
  }
  if (self->LS != NULL)
  {
    SUNLinSolFree (self->LS);
    self->LS = NULL;
  }
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_growth_func_parent_class)->finalize (object);
}

static void
nc_growth_func_class_init (NcGrowthFuncClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->dispose  = &_nc_growth_func_dispose;
  object_class->finalize = &_nc_growth_func_finalize;
}

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
  NcGrowthFunc *gf = g_object_new (NC_TYPE_GROWTH_FUNC, 
                                   NULL);
  return gf;
}

/**
 * nc_growth_func_ref:
 * @gf: a #NcGrowthFunc
 *
 * Increases the reference count of @gf atomically.
 *
 * Returns: (transfer full): @gf.
   */
NcGrowthFunc *
nc_growth_func_ref (NcGrowthFunc *gf)
{
  return g_object_ref (gf);
}

/**
 * nc_growth_func_free:
 * @gf: a #NcGrowthFunc
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
 * @gf: a #NcGrowthFunc
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
growth_f (realtype a, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHICosmo *cosmo       = NC_HICOSMO (f_data);
  const gdouble a2       = a * a;
  const gdouble a5       = a2 * gsl_pow_3 (a);
  const gdouble z        = 1.0 / a - 1.0;
  const gdouble E2       = nc_hicosmo_E2 (cosmo, z);
  const gdouble E        = sqrt (E2);
  const gdouble dE2dz    = nc_hicosmo_dE2_dz (cosmo, z);
  
  const gdouble Omega_m0 = nc_hicosmo_Omega_m0 (cosmo);
  const gdouble D        = NV_Ith_S (y, 0);
  const gdouble B        = NV_Ith_S (y, 1);

  NV_Ith_S (ydot, 0) = B;
  NV_Ith_S (ydot, 1) = (dE2dz / (2.0 * a2 * E2) - 3.0 / a) * B + 3.0 * Omega_m0 * D / (2.0 * a5 * E2);
  NV_Ith_S (ydot, 2) = 1.0 / gsl_pow_3 (a * E);
  
  return 0;
}

static gint
growth_J (realtype a, N_Vector y, N_Vector fy, SUNMatrix J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHICosmo *cosmo       = NC_HICOSMO (jac_data);
  const gdouble a2       = a * a;
  const gdouble a5       = a2 * gsl_pow_3 (a);
  const gdouble z        = 1.0 / a - 1.0;
  const gdouble E2       = nc_hicosmo_E2 (cosmo, z);
  const gdouble dE2dz    = nc_hicosmo_dE2_dz (cosmo, z);
  const gdouble Omega_m0 = nc_hicosmo_Omega_m0 (cosmo);

  NCM_UNUSED (y);
  NCM_UNUSED (fy);
  NCM_UNUSED (tmp1);
  NCM_UNUSED (tmp2);
  NCM_UNUSED (tmp3);
  
  SUN_DENSE_ACCESS (J, 0, 0) = 0.0;
  SUN_DENSE_ACCESS (J, 0, 1) = 1.0;
  SUN_DENSE_ACCESS (J, 0, 2) = 0.0;

  SUN_DENSE_ACCESS (J, 1, 0) = 3.0 * Omega_m0 / (2.0 * a5 * E2);
  SUN_DENSE_ACCESS (J, 1, 1) = (dE2dz / (2.0 * a2 * E2) - 3.0 / a);
  SUN_DENSE_ACCESS (J, 1, 2) = 0.0;

  SUN_DENSE_ACCESS (J, 2, 0) = 0.0;
  SUN_DENSE_ACCESS (J, 2, 1) = 0.0;
  SUN_DENSE_ACCESS (J, 2, 2) = 0.0;
  
  return 0;
}

#define _NC_GROWTH_FUNC_START_A (1.0e-12)

/**
 * nc_growth_func_prepare:
 * @gf: a #NcGrowthFunc
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 */
void
nc_growth_func_prepare (NcGrowthFunc *gf, NcHICosmo *cosmo)
{
  NcGrowthFuncPrivate * const self = gf->priv;
  GArray *x_array, *y_array;
  gdouble ai, a;
  const gdouble Omega_m0 = nc_hicosmo_Omega_m0 (cosmo);
  const gdouble Omega_r0 = nc_hicosmo_Omega_r0 (cosmo);
  gdouble dDa0;
  gint flag;

  ai = _NC_GROWTH_FUNC_START_A;

  if (gf->s != NULL)
  {
    x_array = ncm_vector_get_array (gf->s->xv);
    y_array = ncm_vector_get_array (gf->s->yv);
  }
  else
  {
    x_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    y_array = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  }

  g_array_set_size (x_array, 0);
  g_array_set_size (y_array, 0);

  dDa0 = 1.0e-20 / gsl_pow_3 (ai * nc_hicosmo_E (cosmo, 1.0 / ai - 1.0));
  
  NV_Ith_S (self->yv, 0) = 1.0;
  NV_Ith_S (self->yv, 1) = (Omega_r0 > 0.0) ? ((3.0 / 2.0) * Omega_m0 / Omega_r0) : 1.0;
  NV_Ith_S (self->yv, 2) = dDa0;

  if (self->cvode == NULL)
  {
    self->cvode = CVodeCreate (CV_BDF);

    flag = CVodeInit (self->cvode, &growth_f, ai, self->yv);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

    flag = CVodeSetLinearSolver (self->cvode, self->LS, self->A);
    NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

    flag = CVodeSetJacFn (self->cvode, &growth_J);
    NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );
  }
  else
  {
    flag = CVodeReInit (self->cvode, ai, self->yv);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );

    flag = CVodeSetLinearSolver (self->cvode, self->LS, self->A);
    NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

    flag = CVodeSetJacFn (self->cvode, &growth_J);
    NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );
  }

  flag = CVodeSStolerances (self->cvode, 1e-13, 0.0);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );
  
  flag = CVodeSetMaxNumSteps (self->cvode, 500000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );
  
  flag = CVodeSetUserData (self->cvode, cosmo);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );
  
  flag = CVodeSetStopTime (self->cvode, 1.0);
  NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

  g_array_append_val (x_array, ai);
  g_array_append_val (y_array, NV_Ith_S (self->yv, 0));

  while (TRUE)
  {
    gint flag = CVode (self->cvode, 1.0, self->yv, &a, CV_ONE_STEP);
    NCM_CVODE_CHECK (&flag, "CVode", 1, );

    g_array_append_val (x_array, a);
    g_array_append_val (y_array, NV_Ith_S (self->yv, 0));
    
    if (a == 1.0)
      break;
  }

  gf->Da0 = 2.5 * Omega_m0 * (NV_Ith_S (self->yv, 2) - dDa0);

  if (gf->s == NULL)
  {
    gf->s = ncm_spline_cubic_notaknot_new ();
  }

  {
    NcmVector *xv = ncm_vector_new_array (x_array);
    NcmVector *yv = ncm_vector_new_array (y_array);

    ncm_vector_scale (yv, 1.0 / ncm_vector_get (yv, y_array->len - 1));
    ncm_spline_set (gf->s, xv, yv, TRUE);

    ncm_vector_free (xv);
    ncm_vector_free (yv);
  }

  g_array_unref (x_array);
  g_array_unref (y_array);

  ncm_model_ctrl_update (self->ctrl_cosmo, NCM_MODEL (cosmo));
  return;
}

/**
 * nc_growth_func_prepare_if_needed:
 * @gf: a #NcGrowthFunc
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 */
void
nc_growth_func_prepare_if_needed (NcGrowthFunc *gf, NcHICosmo *cosmo)
{
  NcGrowthFuncPrivate * const self = gf->priv;
  gboolean cosmo_up = ncm_model_ctrl_update (self->ctrl_cosmo, NCM_MODEL (cosmo));

  if (cosmo_up)
    nc_growth_func_prepare (gf, cosmo);
}

/**
 * nc_growth_func_eval:
 * @gf: a #NcGrowthFunc
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: The normalized growth function at @z.
 */
/**
 * nc_growth_func_eval_deriv:
 * @gf: a #NcGrowthFunc
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_growth_func_eval_both:
 * @gf: a #NcGrowthFunc
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 * @d: (out): Growth function
 * @f: (out): Growth function derivative
 *
 * FIXME
 *
 */
/**
 * nc_growth_func_get_dust_norma_Da0:
 * @gf: a #NcGrowthFunc
 *
 * FIXME
 *
 */
