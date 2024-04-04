/***************************************************************************
 *            ncm_powspec.c
 *
 *  Tue February 16 17:00:52 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_powspec.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:ncm_powspec
 * @title: NcmPowspec
 * @short_description: Abstrac class for power spectrum implementation.
 * @stability: Stable
 * @include: numcosmo/math/ncm_powspec.h
 *
 * This module comprises the set of functions to compute a power spectrum and
 * derived quantities.
 *
 * Given a field $\delta(\vec{x})$ at position $\vec{x}$, the power spectrum is
 * defined as the Fourier transform of the two-point correlation point, i.e.,
 * $$\xi(\vec{x} - \vec{x}^\prime) = \int \frac{d^3 k}{(2 \pi)^3} e^{i \vec{k}.(\vec{x} - \vec{x}^\prime)} P(k),$$
 * where $\langle \delta(\vec{k} - \vec{k}^\prime)\rangle = (2\pi)^3 \delta_D(\vec{k} - \vec{k}^\prime) P(k)$
 * and $\delta_D$ is the Dirac's delta function.
 * The standard output is the dimensional power spectrum, not the dimensionless one $\Delta(k)^2$,
 * $$P(k) \equiv \frac{2 \pi^2 \Delta(k)^2}{k^3}.$$
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_powspec.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_integral1d_ptr.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_sf_sbessel.h"
#include "math/ncm_c.h"
#include "math/ncm_spline2d_bicubic.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_bessel.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmPowspecPrivate
{
  /*< private > */
  GObject parent_instance;
  gdouble zi;
  gdouble zf;
  gdouble kmin;
  gdouble kmax;
  NcmIntegral1dPtr *var_tophat_R;
  NcmIntegral1dPtr *corr3D;
  NcmIntegral1dPtr *sproj;
  gdouble reltol_spline;
  NcmModelCtrl *ctrl;
} NcmPowspecPrivate;

enum
{
  PROP_0,
  PROP_ZI,
  PROP_ZF,
  PROP_KMIN,
  PROP_KMAX,
  PROP_RELTOL_SPLINE,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmPowspec, ncm_powspec, G_TYPE_OBJECT)

static gdouble _ncm_powspec_var_tophat_R_integ (gpointer user_data, gdouble lnk, gdouble weight);
static gdouble _ncm_powspec_corr3D_integ (gpointer user_data, gdouble lnk, gdouble weight);
static gdouble _ncm_powspec_sproj_integ (gpointer user_data, gdouble lnk, gdouble weight);

static void
ncm_powspec_init (NcmPowspec *powspec)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  self->zi            = 0.0;
  self->zf            = 0.0;
  self->kmin          = 0.0;
  self->kmax          = 0.0;
  self->reltol_spline = 0.0;

  self->var_tophat_R = ncm_integral1d_ptr_new (&_ncm_powspec_var_tophat_R_integ, NULL);
  self->corr3D       = ncm_integral1d_ptr_new (&_ncm_powspec_corr3D_integ, NULL);
  self->sproj        = ncm_integral1d_ptr_new (&_ncm_powspec_sproj_integ, NULL);

  self->ctrl = ncm_model_ctrl_new (NULL);
}

static void
_ncm_powspec_dispose (GObject *object)
{
  NcmPowspec *powspec            = NCM_POWSPEC (object);
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  ncm_model_ctrl_clear (&self->ctrl);
  ncm_integral1d_ptr_clear (&self->var_tophat_R);
  ncm_integral1d_ptr_clear (&self->corr3D);
  ncm_integral1d_ptr_clear (&self->sproj);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_parent_class)->dispose (object);
}

static void
_ncm_powspec_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_parent_class)->finalize (object);
}

static void
_ncm_powspec_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPowspec *powspec = NCM_POWSPEC (object);

  g_return_if_fail (NCM_IS_POWSPEC (object));

  switch (prop_id)
  {
    case PROP_ZI:
      ncm_powspec_set_zi (powspec, g_value_get_double (value));
      break;
    case PROP_ZF:
      ncm_powspec_set_zf (powspec, g_value_get_double (value));
      break;
    case PROP_KMIN:
      ncm_powspec_set_kmin (powspec, g_value_get_double (value));
      break;
    case PROP_KMAX:
      ncm_powspec_set_kmax (powspec, g_value_get_double (value));
      break;
    case PROP_RELTOL_SPLINE:
      ncm_powspec_set_reltol_spline (powspec, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_powspec_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPowspec *powspec = NCM_POWSPEC (object);

  g_return_if_fail (NCM_IS_POWSPEC (object));

  switch (prop_id)
  {
    case PROP_ZI:
      g_value_set_double (value, ncm_powspec_get_zi (powspec));
      break;
    case PROP_ZF:
      g_value_set_double (value, ncm_powspec_get_zf (powspec));
      break;
    case PROP_KMIN:
      g_value_set_double (value, ncm_powspec_get_kmin (powspec));
      break;
    case PROP_KMAX:
      g_value_set_double (value, ncm_powspec_get_kmax (powspec));
      break;
    case PROP_RELTOL_SPLINE:
      g_value_set_double (value, ncm_powspec_get_reltol_spline (powspec));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_powspec_prepare (NcmPowspec *powspec, NcmModel *model)
{
  g_error ("_ncm_powspec_prepare: no default implementation, all children must implement it.");
}

static gdouble
_ncm_powspec_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  g_error ("_ncm_powspec_eval: no default implementation, all children must implement it.");

  return 0.0;
}

static void _ncm_powspec_eval_vec (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk);

typedef struct __NcmPowspecSplineData
{
  NcmPowspec *powspec;
  NcmModel *model;
  gdouble z_m;
  gdouble k_m;
} _NcmPowspecSplineData;

static gdouble
__P_z (gdouble z, gpointer p)
{
  _NcmPowspecSplineData *data = (_NcmPowspecSplineData *) p;

  return ncm_powspec_eval (data->powspec, data->model, z, data->k_m);
}

static gdouble
__P_k (gdouble k, gpointer p)
{
  _NcmPowspecSplineData *data = (_NcmPowspecSplineData *) p;

  return ncm_powspec_eval (data->powspec, data->model, data->z_m, k);
}

static NcmSpline2d *
_ncm_powspec_get_spline_2d (NcmPowspec *powspec, NcmModel *model)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);
  NcmSpline2d *pk_s2d            = ncm_spline2d_bicubic_notaknot_new ();
  _NcmPowspecSplineData data     = {
    powspec, model,
    0.5 * (self->zi + self->zf),
    exp (0.5 * (log (self->kmin) + log (self->kmax)))
  };
  gsl_function Fx, Fy;
  guint i, j;

  Fx.function = __P_z;
  Fx.params   = &data;

  Fy.function = __P_k;
  Fy.params   = &data;

  ncm_spline2d_set_function (pk_s2d,
                             NCM_SPLINE_FUNCTION_SPLINE,
                             &Fx, &Fy,
                             self->zi, self->zf,
                             self->kmin, self->kmax,
                             self->reltol_spline);
  {
    NcmVector *xv = ncm_spline2d_peek_xv (pk_s2d);
    NcmVector *yv = ncm_spline2d_peek_yv (pk_s2d);
    NcmMatrix *zm = ncm_spline2d_peek_zm (pk_s2d);

    const guint nz = ncm_vector_len (xv);
    const guint nk = ncm_vector_len (yv);

    for (i = 0; i < nk; i++)
    {
      const gdouble k = ncm_vector_get (yv, i);

      for (j = 0; j < nz; j++)
      {
        const gdouble z   = ncm_vector_get (xv, j);
        const gdouble Pkz = ncm_powspec_eval (powspec, model, z, k);

        ncm_matrix_set (zm, i, j, Pkz);
      }
    }
  }
  ncm_spline2d_prepare (pk_s2d);

  return pk_s2d;
}

static void
ncm_powspec_class_init (NcmPowspecClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_powspec_set_property;
  object_class->get_property = &_ncm_powspec_get_property;

  object_class->dispose  = &_ncm_powspec_dispose;
  object_class->finalize = &_ncm_powspec_finalize;

  /**
   * NcmPowspec:zi:
   *
   * The initial time (redshift) to compute $P(k,z)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_ZI,
                                   g_param_spec_double ("zi",
                                                        NULL,
                                                        "Initial time",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspec:zf:
   *
   * The final time (redshift) to compute $P(k,z)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_ZF,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "Final time",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspec:kmin:
   *
   * The minimum mode (wave-number) value to compute $P(k,z)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_KMIN,
                                   g_param_spec_double ("kmin",
                                                        NULL,
                                                        "Minimum mode value",
                                                        0.0, G_MAXDOUBLE, 1.0e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspec:kmax:
   *
   * The maximum mode (wave-number) value to compute $P(k,z)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_KMAX,
                                   g_param_spec_double ("kmax",
                                                        NULL,
                                                        "Maximum mode value",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmPowspec:reltol:
   *
   * The relative tolerance on the interpolation error.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL_SPLINE,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance on the interpolation error",
                                                        GSL_DBL_EPSILON, 1.0, sqrt (GSL_DBL_EPSILON),
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->prepare       = &_ncm_powspec_prepare;
  klass->eval          = &_ncm_powspec_eval;
  klass->eval_vec      = &_ncm_powspec_eval_vec;
  klass->get_spline_2d = &_ncm_powspec_get_spline_2d;
}

static void
_ncm_powspec_eval_vec (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk)
{
  const guint len = ncm_vector_len (k);
  guint i;

  for (i = 0; i < len; i++)
  {
    const gdouble ki  = ncm_vector_get (k, i);
    const gdouble Pki = ncm_powspec_eval (powspec, model, z, ki);

    ncm_vector_set (Pk, i, Pki);
  }

  return;
}

/**
 * ncm_powspec_ref:
 * @powspec: a #NcmPowspec
 *
 * Increases the reference count of @powspec by one atomically.
 *
 * Returns: (transfer full): @powspec.
 */
NcmPowspec *
ncm_powspec_ref (NcmPowspec *powspec)
{
  return g_object_ref (powspec);
}

/**
 * ncm_powspec_free:
 * @powspec: a #NcmPowspec
 *
 * Atomically decrements the reference count of @powspec by one.
 * If the reference count drops to 0,
 * all memory allocated by @powspec is released.
 *
 */
void
ncm_powspec_free (NcmPowspec *powspec)
{
  g_object_unref (powspec);
}

/**
 * ncm_powspec_clear:
 * @powspec: a #NcmPowspec
 *
 * If @powspec is different from NULL,
 * atomically decrements the reference count of @powspec by one.
 * If the reference count drops to 0,
 * all memory allocated by @powspec is released and @powspec is set to NULL.
 *
 */
void
ncm_powspec_clear (NcmPowspec **powspec)
{
  g_clear_object (powspec);
}

/**
 * ncm_powspec_set_zi:
 * @powspec: a #NcmPowspec
 * @zi: initial time $z_i$
 *
 * Sets the initial time $z_i$.
 *
 */
void
ncm_powspec_set_zi (NcmPowspec *powspec, const gdouble zi)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  if (self->zi != zi)
  {
    self->zi = zi;
    ncm_model_ctrl_force_update (self->ctrl);
  }
}

/**
 * ncm_powspec_set_zf:
 * @powspec: a #NcmPowspec
 * @zf: final time $z_f$
 *
 * Sets the final time $z_f$.
 *
 */
void
ncm_powspec_set_zf (NcmPowspec *powspec, const gdouble zf)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  if (self->zf != zf)
  {
    self->zf = zf;
    ncm_model_ctrl_force_update (self->ctrl);
  }
}

/**
 * ncm_powspec_set_kmin:
 * @powspec: a #NcmPowspec
 * @kmin: minimum mode $k_\mathrm{min}$
 *
 * Sets the minimum mode value $k_\mathrm{min}$.
 *
 */
void
ncm_powspec_set_kmin (NcmPowspec *powspec, const gdouble kmin)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  if (self->kmin != kmin)
  {
    self->kmin = kmin;
    ncm_model_ctrl_force_update (self->ctrl);
  }
}

/**
 * ncm_powspec_set_kmax:
 * @powspec: a #NcmPowspec
 * @kmax: maxmimum mode $k_\mathrm{max}$
 *
 * Sets the maximum mode value $k_\mathrm{max}$.
 *
 */
void
ncm_powspec_set_kmax (NcmPowspec *powspec, const gdouble kmax)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  if (self->kmax != kmax)
  {
    self->kmax = kmax;
    ncm_model_ctrl_force_update (self->ctrl);
  }
}

/**
 * ncm_powspec_set_reltol_spline:
 * @powspec: a #NcmPowspec
 * @reltol: relative tolerance for interpolation errors
 *
 * Sets the relative tolerance for interpolation errors to @reltol.
 *
 */
void
ncm_powspec_set_reltol_spline (NcmPowspec *powspec, const gdouble reltol)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  self->reltol_spline = reltol;
}

/**
 * ncm_powspec_require_zi:
 * @powspec: a #NcmPowspec
 * @zi: initial time $z_i$
 *
 * Requires the initial time to be less or equal to $z_i$.
 *
 */
void
ncm_powspec_require_zi (NcmPowspec *powspec, const gdouble zi)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  if (zi < self->zi)
    ncm_powspec_set_zi (powspec, zi);
}

/**
 * ncm_powspec_require_zf:
 * @powspec: a #NcmPowspec
 * @zf: final time $z_f$
 *
 * Requires the final time to be greater or equal to $z_f$.
 *
 */
void
ncm_powspec_require_zf (NcmPowspec *powspec, const gdouble zf)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  if (zf > self->zf)
    ncm_powspec_set_zf (powspec, zf);
}

/**
 * ncm_powspec_require_kmin:
 * @powspec: a #NcmPowspec
 * @kmin: minimum mode $k_\mathrm{min}$
 *
 * Requires the minimum mode value to be less or equal to $k_\mathrm{min}$.
 *
 */
void
ncm_powspec_require_kmin (NcmPowspec *powspec, const gdouble kmin)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  if (kmin < self->kmin)
    ncm_powspec_set_kmin (powspec, kmin);
}

/**
 * ncm_powspec_require_kmax:
 * @powspec: a #NcmPowspec
 * @kmax: maxmimum mode $k_\mathrm{max}$
 *
 * Sets the maximum mode value $k_\mathrm{max}$.
 *
 */
void
ncm_powspec_require_kmax (NcmPowspec *powspec, const gdouble kmax)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  if (kmax > self->kmax)
    ncm_powspec_set_kmax (powspec, kmax);
}

/**
 * ncm_powspec_get_zi:
 * @powspec: a #NcmPowspec
 *
 * Gets the initial value $z_i$.
 *
 */
gdouble
ncm_powspec_get_zi (NcmPowspec *powspec)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  return self->zi;
}

/**
 * ncm_powspec_get_zf:
 * @powspec: a #NcmPowspec
 *
 * Gets the final value $z_f$.
 *
 */
gdouble
ncm_powspec_get_zf (NcmPowspec *powspec)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  return self->zf;
}

/**
 * ncm_powspec_get_kmin:
 * @powspec: a #NcmPowspec
 *
 * Gets the minimum mode value $k_\mathrm{min}$.
 *
 */
gdouble
ncm_powspec_get_kmin (NcmPowspec *powspec)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  return self->kmin;
}

/**
 * ncm_powspec_get_kmax:
 * @powspec: a #NcmPowspec
 *
 * Gets the maximum mode value $k_\mathrm{max}$.
 *
 */
gdouble
ncm_powspec_get_kmax (NcmPowspec *powspec)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  return self->kmax;
}

/**
 * ncm_powspec_get_reltol_spline:
 * @powspec: a #NcmPowspec
 *
 * Gets the relative tolerance for interpolation errors.
 *
 * Returns: the relative tolerance for interpolation errors.
 */
gdouble
ncm_powspec_get_reltol_spline (NcmPowspec *powspec)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  return self->reltol_spline;
}

/**
 * ncm_powspec_get_nknots: (virtual get_nknots)
 * @powspec: a #NcmPowspec
 * @Nz: (out): number of knots in $z$
 * @Nk: (out): number of knots in $k$
 *
 * Gets the number of knots used to calculate the power spectrum.
 *
 */
void
ncm_powspec_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk)
{
  NCM_POWSPEC_GET_CLASS (powspec)->get_nknots (powspec, Nz, Nk);
}

/**
 * ncm_powspec_prepare:
 * @powspec: a #NcmPowspec
 * @model: a #NcmModel
 *
 * Prepares the power spectrum @powspec using the model @model.
 *
 */

void
ncm_powspec_prepare (NcmPowspec *powspec, NcmModel *model)
{
  NCM_POWSPEC_GET_CLASS (powspec)->prepare (powspec, model);
}

/**
 * ncm_powspec_prepare_if_needed:
 * @powspec: a #NcmPowspec
 * @model: a #NcmModel
 *
 * Prepares the object @powspec using the model @model if it was changed
 * since last preparation.
 *
 */
void
ncm_powspec_prepare_if_needed (NcmPowspec *powspec, NcmModel *model)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);
  gboolean model_up              = ncm_model_ctrl_update (self->ctrl, NCM_MODEL (model));

  if (model_up)
    ncm_powspec_prepare (powspec, model);
}

/**
 * ncm_powspec_eval:
 * @powspec: a #NcmPowspec
 * @model: a #NcmModel
 * @z: time $z$
 * @k: mode $k$
 *
 * Evaluates the power spectrum @powspec at $(z, k)$.
 *
 * Returns: $P(z, k)$.
 */
gdouble
ncm_powspec_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  return NCM_POWSPEC_GET_CLASS (powspec)->eval (powspec, model, z, k);
}

/**
 * ncm_powspec_eval_vec:
 * @powspec: a #NcmPowspec
 * @model: a #NcmModel
 * @z: time $z$
 * @k: a #NcmVector
 * @Pk: a #NcmVector
 *
 * Evaluates the power spectrum @powspec at $z$ and in the knots
 * contained in @k and puts the result in @Pk.
 *
 */
void
ncm_powspec_eval_vec (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk)
{
  NCM_POWSPEC_GET_CLASS (powspec)->eval_vec (powspec, model, z, k, Pk);
}

/**
 * ncm_powspec_get_spline_2d:
 * @powspec: a #NcmPowspec
 * @model: a #NcmModel compatible with @powspec
 *
 * Compute a 2D spline for the power spectrum.
 *
 * Returns: (transfer full): a #NcmSpline2d interpolating spline as a function of $(z, k)$.
 */
NcmSpline2d *
ncm_powspec_get_spline_2d (NcmPowspec *powspec, NcmModel *model)
{
  return NCM_POWSPEC_GET_CLASS (powspec)->get_spline_2d (powspec, model);
}

/**
 * ncm_powspec_peek_model_ctrl:
 * @powspec: a #NcmPowspec
 *
 * Gets the #NcmModelCtrl used by @powspec.
 *
 * Returns: (transfer none): the #NcmModelCtrl used by @powspec.
 */
NcmModelCtrl *
ncm_powspec_peek_model_ctrl (NcmPowspec *powspec)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  return self->ctrl;
}

typedef struct _NcmPowspecInt
{
  const gdouble z;
  const gdouble R;
  NcmPowspec *ps;
  NcmModel *model;
  const gdouble z2;
  const gdouble xi1;
  const gdouble xi2;
  const gint ell;
} NcmPowspecInt;

static gdouble
_ncm_powspec_var_tophat_R_integ (gpointer user_data, gdouble lnk, gdouble weight)
{
  NcmPowspecInt *data = (NcmPowspecInt *) user_data;
  const gdouble k     = exp (lnk);
  const gdouble x     = k * data->R;
  const gdouble Pk    = ncm_powspec_eval (data->ps, data->model, data->z, k);
  const gdouble W     = 3.0 * gsl_sf_bessel_j1 (x) / x;
  const gdouble W2    = W * W;

  return gsl_pow_3 (k) * Pk * W2;
}

/**
 * ncm_powspec_var_tophat_R:
 * @powspec: a #NcmPowspec
 * @model: a #NcmModel
 * @reltol: relative tolerance for integration
 * @z: the value of $z$
 * @R: the value of $R$
 *
 * This function computes the value of the linearly extrapolated
 * rms fluctuations of mass in a sphere of radius $R$ applying a top-hat filter at redshift $z$,
 * $$\sigma_{R}^{2}(z) = \frac{1}{2\pi^2} \int_{k_{\mathrm{min}}}^{k_\mathrm{max}}  W^{2}_{TH}(kR) \, P(k,z) \, k^2 \, \mathrm{d}k \, .$$
 * Where, $W_{TH}(t)$ is the top-hat filter in Fourier space,
 * $$W_{TH}(t) = \frac{3}{t^3} \left( \sin t - t \cos t  \right) = \frac{3}{t} j_{1}(t),$$
 * and $j_1(t)$ is the first order spherical Bessel function of the first kind.
 * This function is recommended for a small set of $z$ and $R$ values.
 * For a wide range of values it is best to apply #NcmPowspecFilter, instead.
 *
 * Returns: $\sigma_{R}^{2}(z)$.
 */
gdouble
ncm_powspec_var_tophat_R (NcmPowspec *powspec, NcmModel *model, const gdouble reltol, const gdouble z, const gdouble R)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);
  NcmPowspecInt data             = {z, R, powspec, model, 0.0, 0.0, 0.0, 0};
  const gdouble kmin             = ncm_powspec_get_kmin (powspec);
  const gdouble kmax             = ncm_powspec_get_kmax (powspec);
  const gdouble lnkmin           = log (kmin);
  const gdouble lnkmax           = log (kmax);
  const gdouble one_2pi2         = 1.0 / ncm_c_two_pi_2 ();
  gdouble error, sigma2_2pi2;

  ncm_powspec_prepare_if_needed (powspec, model);

  ncm_integral1d_ptr_set_userdata (self->var_tophat_R, &data);
  ncm_integral1d_set_reltol (NCM_INTEGRAL1D (self->var_tophat_R), reltol);

  sigma2_2pi2 = ncm_integral1d_eval (NCM_INTEGRAL1D (self->var_tophat_R), lnkmin, lnkmax, &error);

  return sigma2_2pi2 * one_2pi2;
}

/**
 * ncm_powspec_sigma_tophat_R:
 * @powspec: a #NcmPowspec
 * @model: a #NcmModel
 * @reltol: relative tolerance for integration
 * @z: the value of $z$
 * @R: the value of $R$
 *
 * Computes $\sigma_R(z) = \sqrt{\sigma_{R}^{2}(z)}$. See ncm_powspec_var_tophat_R().
 *
 * Returns: $\sigma_R(z)$.
 */
gdouble
ncm_powspec_sigma_tophat_R (NcmPowspec *powspec, NcmModel *model, const gdouble reltol, const gdouble z, const gdouble R)
{
  return sqrt (ncm_powspec_var_tophat_R (powspec, model, reltol, z, R));
}

static gdouble
_ncm_powspec_corr3D_integ (gpointer user_data, gdouble lnk, gdouble weight)
{
  NcmPowspecInt *data = (NcmPowspecInt *) user_data;
  const gdouble k     = exp (lnk);
  const gdouble x     = k * data->R;
  const gdouble Pk    = ncm_powspec_eval (data->ps, data->model, data->z, k);
  const gdouble W     = gsl_sf_bessel_j0 (x);

  return gsl_pow_3 (k) * Pk * W;
}

/**
 * ncm_powspec_corr3d:
 * @powspec: a #NcmPowspec
 * @model: a #NcmModel
 * @reltol: relative tolerance for integration
 * @z: the value of $z$
 * @r: the value of $r$
 *
 * Computes the spatial correlation function in configuration space at redshift $z$ and position $r$,
 * $$\xi(r,z) = \frac{1}{2\pi^2} \int_{k_{\mathrm{min}}}^{k_\mathrm{max}} P(k,z) \, j_{0}(kr) \, k^2 \, \mathrm{d}k \, ,$$
 * where, $j_0(t)$ is the zero order spherical Bessel function of the first kind.
 * This function is recommended for a small set of $r$ and $z$ values.
 * For a wide range of values it is best to apply #NcmPowspecCorr3d, instead.
 *
 * Returns: $\xi(r,z)$.
 */
gdouble
ncm_powspec_corr3d (NcmPowspec *powspec, NcmModel *model, const gdouble reltol, const gdouble z, const gdouble r)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);
  NcmPowspecInt data             = {z, r, powspec, model, 0.0, 0.0, 0.0, 0};
  const gdouble kmin             = ncm_powspec_get_kmin (powspec);
  const gdouble kmax             = ncm_powspec_get_kmax (powspec);
  const gdouble lnkmin           = log (kmin);
  const gdouble lnkmax           = log (kmax);
  const gdouble one_2pi2         = 1.0 / ncm_c_two_pi_2 ();
  gdouble error, xi_2pi2;

  ncm_powspec_prepare_if_needed (powspec, model);

  ncm_integral1d_ptr_set_userdata (self->corr3D, &data);
  ncm_integral1d_set_reltol (NCM_INTEGRAL1D (self->corr3D), reltol);

  xi_2pi2 = ncm_integral1d_eval (NCM_INTEGRAL1D (self->corr3D), lnkmin, lnkmax, &error);

  return xi_2pi2 * one_2pi2;
}

static gdouble
_ncm_powspec_sproj_integ (gpointer user_data, gdouble lnk, gdouble weight)
{
  NcmPowspecInt *data = (NcmPowspecInt *) user_data;
  const gdouble k     = exp (lnk);
  const gdouble x1    = k * data->xi1;
  const gdouble x2    = k * data->xi2;
  const gdouble Pk    = sqrt (ncm_powspec_eval (data->ps, data->model, data->z, k) * ncm_powspec_eval (data->ps, data->model, data->z2, k));
  const gdouble W     = ncm_sf_sbessel (data->ell, x1) * ncm_sf_sbessel (data->ell, x2);

  return gsl_pow_3 (k) * Pk * W;
}

/**
 * ncm_powspec_sproj:
 * @powspec: a #NcmPowspec
 * @model: a #NcmModel
 * @reltol: relative tolerance for integration
 * @ell: the value of $\ell$
 * @z1: the value of $z_1$
 * @z2: the value of $z_2$
 * @xi1: the value of $\xi_1$
 * @xi2: the value of $\xi_2$
 *
 * Computes \(C_\ell (z_1, z_2) = \int\dots\). This method calculates the angular power
 * spectrum directly from the power spectrum by integrating over the wave-numbers. It
 * is slow and intended for testing purposes only.
 *
 */
gdouble
ncm_powspec_sproj (NcmPowspec *powspec, NcmModel *model, const gdouble reltol, const gint ell, const gdouble z1, const gdouble z2, const gdouble xi1, const gdouble xi2)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);
  NcmPowspecInt data             = {z1, 0.0, powspec, model, z2, xi1, xi2, ell};
  const gdouble kmin             = ncm_powspec_get_kmin (powspec);
  const gdouble kmax             = ncm_powspec_get_kmax (powspec);
  const gdouble lnkmin           = log (kmin);
  const gdouble lnkmax           = log (kmax);
  const gdouble two_pi           = 2.0 / ncm_c_pi ();
  gdouble error, xi_two_pi;

  ncm_powspec_prepare_if_needed (powspec, model);

  ncm_integral1d_ptr_set_userdata (self->sproj, &data);
  ncm_integral1d_set_reltol (NCM_INTEGRAL1D (self->sproj), reltol);

  xi_two_pi = ncm_integral1d_eval (NCM_INTEGRAL1D (self->sproj), lnkmin, lnkmax, &error);

  return xi_two_pi * two_pi;
}

