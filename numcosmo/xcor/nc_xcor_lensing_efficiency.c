/***************************************************************************
 *            nc_xcor_lensing_efficiency.c
 *
 *  Thu January 02 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_xcor_lensing_efficiency.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcXcorLensingEfficiency:
 *
 * Abstract class for computing lensing efficiency.
 *
 * This class provides infrastructure for computing the lensing efficiency function:
 * \begin{equation}
 * g(z) = \int_z^{z_{\max}} dz' \left(1 - \frac{\chi(z)}{\chi(z')}\right) W_{\mathrm{src}}(z')
 * \end{equation}
 * where $W_{\mathrm{src}}(z')$ is a source weight function that must be implemented by
 * subclasses.
 * The integration is performed using CVODE (BDF method) to solve the equivalent ODE
 * system:
 * \begin{align}
 * \frac{df}{d(-z)} &= -\frac{g}{E(z)} \\
 * \frac{dg}{d(-z)} &= -\frac{W_{\mathrm{src}}(z)}{d_t(z)} - \Omega_{k0} \frac{f}{E(z)}
 * \end{align}
 * where $f(z) = g(z)$, $E(z) = H(z)/H_0$, $d_t(z)$ is the transverse comoving
 * distance, and the integration is performed backwards from $z_{\max}$ to $z \approx
 * 0$.
 *
 * Subclasses must implement two virtual methods:
 *
 * - eval_source(): returns $W_{\mathrm{src}}(z)$
 * - get_z_range(): returns the source redshift range
 *
 * Example subclasses:
 *
 * - Weak lensing: $W_{\mathrm{src}}(z) = \frac{dn}{dz}$ (source galaxy distribution)
 * - Galaxy magnification bias: $W_{\mathrm{src}}(z) = (5s-2) \frac{dn}{dz}$ (weighted
 *   by magnification bias parameter)
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_cfg.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "xcor/nc_xcor_lensing_efficiency.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <cvode/cvode.h>

#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#define SUN_DENSE_ACCESS SM_ELEMENT_D

#include <nvector/nvector_serial.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcXcorLensingEfficiencyPrivate
{
  /* ODE solver components */
  gpointer cvode;
  SUNContext sunctx;
  N_Vector yv;
  SUNMatrix A;
  SUNLinearSolver LS;
  gdouble reltol;
  gdouble abstol;

  /* Distance object (needed for comoving distances) */
  NcDistance *dist;

  /* Output: the g(z) function stored as g(-z) */
  NcmSpline *g_mz;

  /* Control */
  NcmModelCtrl *ctrl_cosmo;
} NcXcorLensingEfficiencyPrivate;

enum
{
  PROP_0,
  PROP_DISTANCE,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcXcorLensingEfficiency, nc_xcor_lensing_efficiency, G_TYPE_OBJECT)

static void
nc_xcor_lensing_efficiency_init (NcXcorLensingEfficiency *lens_eff)
{
  NcXcorLensingEfficiencyPrivate *self = nc_xcor_lensing_efficiency_get_instance_private (lens_eff);

  if (SUNContext_Create (SUN_COMM_NULL, &self->sunctx))
    g_error ("ERROR: SUNContext_Create failed\n");

  self->cvode      = NULL;
  self->yv         = N_VNew_Serial (2, self->sunctx);
  self->A          = SUNDenseMatrix (2, 2, self->sunctx);
  self->LS         = SUNLinSol_Dense (self->yv, self->A, self->sunctx);
  self->reltol     = NC_XCOR_LENSING_EFFICIENCY_DEFAULT_RELTOL;
  self->abstol     = NC_XCOR_LENSING_EFFICIENCY_DEFAULT_ABSTOL;
  self->dist       = NULL;
  self->g_mz       = NULL;
  self->ctrl_cosmo = ncm_model_ctrl_new (NULL);

  NCM_CVODE_CHECK ((gpointer) self->A, "SUNDenseMatrix", 0, );
  NCM_CVODE_CHECK ((gpointer) self->LS, "SUNDenseLinearSolver", 0, );
}

static void
_nc_xcor_lensing_efficiency_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorLensingEfficiency *lens_eff = NC_XCOR_LENSING_EFFICIENCY (object);

  g_return_if_fail (NC_IS_XCOR_LENSING_EFFICIENCY (object));

  switch (prop_id)
  {
    case PROP_DISTANCE:
      nc_xcor_lensing_efficiency_set_distance (lens_eff, g_value_get_object (value));
      break;
    case PROP_RELTOL:
      nc_xcor_lensing_efficiency_set_reltol (lens_eff, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      nc_xcor_lensing_efficiency_set_abstol (lens_eff, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_lensing_efficiency_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorLensingEfficiency *lens_eff = NC_XCOR_LENSING_EFFICIENCY (object);

  g_return_if_fail (NC_IS_XCOR_LENSING_EFFICIENCY (object));

  switch (prop_id)
  {
    case PROP_DISTANCE:
      g_value_set_object (value, nc_xcor_lensing_efficiency_peek_distance (lens_eff));
      break;
    case PROP_RELTOL:
      g_value_set_double (value, nc_xcor_lensing_efficiency_get_reltol (lens_eff));
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, nc_xcor_lensing_efficiency_get_abstol (lens_eff));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_lensing_efficiency_dispose (GObject *object)
{
  NcXcorLensingEfficiency *lens_eff    = NC_XCOR_LENSING_EFFICIENCY (object);
  NcXcorLensingEfficiencyPrivate *self = nc_xcor_lensing_efficiency_get_instance_private (lens_eff);

  nc_distance_clear (&self->dist);
  ncm_spline_clear (&self->g_mz);
  ncm_model_ctrl_clear (&self->ctrl_cosmo);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_lensing_efficiency_parent_class)->dispose (object);
}

static void
_nc_xcor_lensing_efficiency_finalize (GObject *object)
{
  NcXcorLensingEfficiency *lens_eff    = NC_XCOR_LENSING_EFFICIENCY (object);
  NcXcorLensingEfficiencyPrivate *self = nc_xcor_lensing_efficiency_get_instance_private (lens_eff);

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

  SUNContext_Free (&self->sunctx);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_lensing_efficiency_parent_class)->finalize (object);
}

static void
nc_xcor_lensing_efficiency_class_init (NcXcorLensingEfficiencyClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_xcor_lensing_efficiency_set_property;
  object_class->get_property = &_nc_xcor_lensing_efficiency_get_property;
  object_class->dispose      = &_nc_xcor_lensing_efficiency_dispose;
  object_class->finalize     = &_nc_xcor_lensing_efficiency_finalize;

  /**
   * NcXcorLensingEfficiency:distance:
   *
   * The #NcDistance object used for comoving distance computations.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DISTANCE,
                                   g_param_spec_object ("distance",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorLensingEfficiency:reltol:
   *
   * Relative tolerance used when integrating the ODE.
   * Default value: $10^{-13}$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance",
                                                        GSL_DBL_EPSILON, 1.0e-1, NC_XCOR_LENSING_EFFICIENCY_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorLensingEfficiency:abstol:
   *
   * Absolute tolerance used when integrating the ODE.
   * Default value: $10^{-50}$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Absolute tolerance",
                                                        0.0, G_MAXDOUBLE, NC_XCOR_LENSING_EFFICIENCY_DEFAULT_ABSTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Virtual methods must be implemented by subclasses */
  klass->eval_source = NULL;
  klass->get_z_range = NULL;
}

/**
 * nc_xcor_lensing_efficiency_ref:
 * @lens_eff: a #NcXcorLensingEfficiency
 *
 * Increases the reference count of @lens_eff by one.
 *
 * Returns: (transfer full): @lens_eff
 */
NcXcorLensingEfficiency *
nc_xcor_lensing_efficiency_ref (NcXcorLensingEfficiency *lens_eff)
{
  return g_object_ref (lens_eff);
}

/**
 * nc_xcor_lensing_efficiency_free:
 * @lens_eff: a #NcXcorLensingEfficiency
 *
 * Decreases the reference count of @lens_eff by one.
 *
 */
void
nc_xcor_lensing_efficiency_free (NcXcorLensingEfficiency *lens_eff)
{
  g_object_unref (lens_eff);
}

/**
 * nc_xcor_lensing_efficiency_clear:
 * @lens_eff: a #NcXcorLensingEfficiency
 *
 * If @lens_eff is different from NULL, decreases the reference count of
 * @lens_eff by one and sets @lens_eff to NULL.
 *
 */
void
nc_xcor_lensing_efficiency_clear (NcXcorLensingEfficiency **lens_eff)
{
  g_clear_object (lens_eff);
}

/**
 * nc_xcor_lensing_efficiency_set_distance:
 * @lens_eff: a #NcXcorLensingEfficiency
 * @dist: a #NcDistance
 *
 * Sets the distance object used for comoving distance computations.
 *
 */
void
nc_xcor_lensing_efficiency_set_distance (NcXcorLensingEfficiency *lens_eff, NcDistance *dist)
{
  NcXcorLensingEfficiencyPrivate *self = nc_xcor_lensing_efficiency_get_instance_private (lens_eff);

  g_assert (dist != NULL);

  nc_distance_clear (&self->dist);
  self->dist = nc_distance_ref (dist);
}

/**
 * nc_xcor_lensing_efficiency_peek_distance:
 * @lens_eff: a #NcXcorLensingEfficiency
 *
 * Gets the distance object.
 *
 * Returns: (transfer none): the #NcDistance object
 */
NcDistance *
nc_xcor_lensing_efficiency_peek_distance (NcXcorLensingEfficiency *lens_eff)
{
  NcXcorLensingEfficiencyPrivate *self = nc_xcor_lensing_efficiency_get_instance_private (lens_eff);

  return self->dist;
}

/**
 * nc_xcor_lensing_efficiency_set_reltol:
 * @lens_eff: a #NcXcorLensingEfficiency
 * @reltol: relative tolerance
 *
 * Sets the relative tolerance for ODE integration.
 *
 */
void
nc_xcor_lensing_efficiency_set_reltol (NcXcorLensingEfficiency *lens_eff, gdouble reltol)
{
  NcXcorLensingEfficiencyPrivate *self = nc_xcor_lensing_efficiency_get_instance_private (lens_eff);

  self->reltol = reltol;
}

/**
 * nc_xcor_lensing_efficiency_set_abstol:
 * @lens_eff: a #NcXcorLensingEfficiency
 * @abstol: absolute tolerance
 *
 * Sets the absolute tolerance for ODE integration.
 *
 */
void
nc_xcor_lensing_efficiency_set_abstol (NcXcorLensingEfficiency *lens_eff, gdouble abstol)
{
  NcXcorLensingEfficiencyPrivate *self = nc_xcor_lensing_efficiency_get_instance_private (lens_eff);

  self->abstol = abstol;
}

/**
 * nc_xcor_lensing_efficiency_get_reltol:
 * @lens_eff: a #NcXcorLensingEfficiency
 *
 * Gets the relative tolerance for ODE integration.
 *
 * Returns: relative tolerance
 */
gdouble
nc_xcor_lensing_efficiency_get_reltol (NcXcorLensingEfficiency *lens_eff)
{
  NcXcorLensingEfficiencyPrivate *self = nc_xcor_lensing_efficiency_get_instance_private (lens_eff);

  return self->reltol;
}

/**
 * nc_xcor_lensing_efficiency_get_abstol:
 * @lens_eff: a #NcXcorLensingEfficiency
 *
 * Gets the absolute tolerance for ODE integration.
 *
 * Returns: absolute tolerance
 */
gdouble
nc_xcor_lensing_efficiency_get_abstol (NcXcorLensingEfficiency *lens_eff)
{
  NcXcorLensingEfficiencyPrivate *self = nc_xcor_lensing_efficiency_get_instance_private (lens_eff);

  return self->abstol;
}

typedef struct _lens_eff_int_params
{
  NcXcorLensingEfficiency *lens_eff;
  NcHICosmo *cosmo;
  const gdouble Omega_k0;
} lens_eff_int_params;

static gint
_lens_eff_ode_f (sunrealtype mz, N_Vector y, N_Vector ydot, gpointer f_data)
{
  lens_eff_int_params *params          = (lens_eff_int_params *) f_data;
  NcXcorLensingEfficiency *lens_eff    = params->lens_eff;
  NcXcorLensingEfficiencyPrivate *self = nc_xcor_lensing_efficiency_get_instance_private (lens_eff);
  NcXcorLensingEfficiencyClass *klass  = NC_XCOR_LENSING_EFFICIENCY_GET_CLASS (lens_eff);
  NcHICosmo *cosmo                     = params->cosmo;
  NcDistance *dist                     = self->dist;
  const gdouble z                      = -mz;
  const gdouble E                      = nc_hicosmo_E (cosmo, z);
  const gdouble W_src                  = klass->eval_source (lens_eff, z);
  const gdouble dt                     = nc_distance_transverse (dist, cosmo, z);

  NV_Ith_S (ydot, 0) = -NV_Ith_S (y, 1) / E;
  NV_Ith_S (ydot, 1) = -W_src / dt - params->Omega_k0 * NV_Ith_S (y, 0) / E;

  return 0;
}

static gint
_lens_eff_ode_J (sunrealtype mz, N_Vector y, N_Vector fy, SUNMatrix J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  lens_eff_int_params *params = (lens_eff_int_params *) jac_data;
  NcHICosmo *cosmo            = params->cosmo;
  const gdouble z             = -mz;
  const gdouble E             = nc_hicosmo_E (cosmo, z);

  SUN_DENSE_ACCESS (J, 0, 0) = 0.0;
  SUN_DENSE_ACCESS (J, 0, 1) = -1.0 / E;

  SUN_DENSE_ACCESS (J, 1, 0) = -params->Omega_k0 / E;
  SUN_DENSE_ACCESS (J, 1, 1) = 0.0;

  return 0;
}

/**
 * nc_xcor_lensing_efficiency_prepare:
 * @lens_eff: a #NcXcorLensingEfficiency
 * @cosmo: a #NcHICosmo
 *
 * Prepares the lensing efficiency object by computing the $g(z)$ function.
 * This involves solving the ODE system backwards from $z_{\max}$ to $z \approx 0$
 * using CVODE.
 *
 */
void
nc_xcor_lensing_efficiency_prepare (NcXcorLensingEfficiency *lens_eff, NcHICosmo *cosmo)
{
  NcXcorLensingEfficiencyPrivate *self = nc_xcor_lensing_efficiency_get_instance_private (lens_eff);
  NcXcorLensingEfficiencyClass *klass  = NC_XCOR_LENSING_EFFICIENCY_GET_CLASS (lens_eff);
  gdouble zmin, zmax;

  g_assert (self->dist != NULL);
  g_assert (klass->eval_source != NULL);
  g_assert (klass->get_z_range != NULL);

  nc_distance_prepare_if_needed (self->dist, cosmo);

  klass->get_z_range (lens_eff, &zmin, &zmax);

  {
    lens_eff_int_params params = {lens_eff, cosmo, nc_hicosmo_Omega_k0 (cosmo)};
    gdouble mz_ini             = -zmax;
    const gdouble mz_end       = -1.0e-10;
    GArray *x_array, *y_array;
    gint flag;

    if (self->g_mz != NULL)
    {
      NcmVector *xv = ncm_spline_peek_xv (self->g_mz);
      NcmVector *yv = ncm_spline_peek_yv (self->g_mz);

      x_array = ncm_vector_get_array (xv);
      y_array = ncm_vector_get_array (yv);

      g_array_set_size (x_array, 0);
      g_array_set_size (y_array, 0);
    }
    else
    {
      self->g_mz = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
      x_array    = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
      y_array    = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    }

    {
      /*
       * Setting up initial conditions based on second order approximation of the
       * original integral.
       */
      const gdouble dz    = zmax * sqrt (GSL_DBL_EPSILON);
      const gdouble W_src = klass->eval_source (lens_eff, zmax);
      const gdouble f     = 0.5 * gsl_pow_2 (dz) * W_src / nc_distance_transverse (self->dist, cosmo, zmax) / nc_hicosmo_E (cosmo, zmax);
      const gdouble g     = -dz * W_src / nc_distance_transverse (self->dist, cosmo, zmax);

      NV_Ith_S (self->yv, 0) = f;
      NV_Ith_S (self->yv, 1) = g;

      mz_ini += dz;
    }

    if (self->cvode == NULL)
    {
      self->cvode = CVodeCreate (CV_BDF, self->sunctx);

      flag = CVodeInit (self->cvode, &_lens_eff_ode_f, mz_ini, self->yv);
      NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

      flag = CVodeSetLinearSolver (self->cvode, self->LS, self->A);
      NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

      flag = CVodeSetJacFn (self->cvode, &_lens_eff_ode_J);
      NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );
    }
    else
    {
      flag = CVodeReInit (self->cvode, mz_ini, self->yv);
      NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );

      flag = CVodeSetLinearSolver (self->cvode, self->LS, self->A);
      NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

      flag = CVodeSetJacFn (self->cvode, &_lens_eff_ode_J);
      NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );
    }

    flag = CVodeSStolerances (self->cvode, self->reltol, self->abstol);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetMaxNumSteps (self->cvode, 500000);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

    flag = CVodeSetUserData (self->cvode, &params);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVodeSetStopTime (self->cvode, mz_end);
    NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

    flag = CVodeSetMaxStep (self->cvode, 1.0e-1);

    g_array_append_val (x_array, mz_ini);
    g_array_append_val (y_array, NV_Ith_S (self->yv, 0));

    while (TRUE)
    {
      flag = CVode (self->cvode, mz_end, self->yv, &mz_ini, CV_ONE_STEP);

      NCM_CVODE_CHECK (&flag, "CVode", 1, );

      g_array_append_val (x_array, mz_ini);
      g_array_append_val (y_array, NV_Ith_S (self->yv, 0));

      if (mz_ini == mz_end)
        break;
    }

    {
      NcmVector *xv = ncm_vector_new_array (x_array);
      NcmVector *yv = ncm_vector_new_array (y_array);

      ncm_spline_set (self->g_mz, xv, yv, TRUE);
      ncm_vector_free (xv);
      ncm_vector_free (yv);

      g_array_unref (x_array);
      g_array_unref (y_array);
    }
  }
}

/**
 * nc_xcor_lensing_efficiency_eval:
 * @lens_eff: a #NcXcorLensingEfficiency
 * @z: redshift
 *
 * Evaluates the lensing efficiency function $g(z)$ at redshift @z.
 * The function must have been prepared first using nc_xcor_lensing_efficiency_prepare().
 *
 * Returns: $g(z)$
 */
gdouble
nc_xcor_lensing_efficiency_eval (NcXcorLensingEfficiency *lens_eff, gdouble z)
{
  NcXcorLensingEfficiencyPrivate *self = nc_xcor_lensing_efficiency_get_instance_private (lens_eff);

  return ncm_spline_eval (self->g_mz, -z);
}

