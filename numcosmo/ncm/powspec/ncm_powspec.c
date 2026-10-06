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
 * NcmPowspec:
 *
 * Abstract class for a power spectrum $P(k, z)$.
 *
 * For a field $\delta(\vec{x})$ with Fourier transform $\tilde{\delta}(\vec{k})$,
 * $$\langle \tilde{\delta}(\vec{k}) \tilde{\delta}^*(\vec{k}^\prime) \rangle = (2\pi)^3 \delta_D(\vec{k} - \vec{k}^\prime) P(k),$$
 * and the two-point correlation function is
 * $$\xi(|\vec{x} - \vec{x}^\prime|) = \int \frac{\mathrm{d}^3 k}{(2\pi)^3} e^{i \vec{k} \cdot (\vec{x} - \vec{x}^\prime)} P(k).$$
 * $P$ is the dimensional spectrum, $P(k) = 2\pi^2 \Delta^2(k) / k^3$, with $k$ in
 * $\mathrm{Mpc}^{-1}$ and $P$ in $\mathrm{Mpc}^3$.
 *
 * A subclass implements prepare(), eval() and get_nknots(); eval_vec(), deriv_z(),
 * deriv_k() and get_spline_2d() have defaults built on eval(). Each of
 * ncm_powspec_var_tophat_R(), ncm_powspec_corr3d() and ncm_powspec_sproj() uses one
 * integrator stored in the object, so concurrent calls on the same object are not safe.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/powspec/ncm_powspec.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/integration/ncm_integral1d_ptr.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/specfunc/ncm_sf_sbessel.h"
#include "ncm/core/ncm_c.h"
#include "ncm/spline/ncm_spline2d_bicubic.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_bessel.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmPowspecPrivate
{
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

/* LCOV_EXCL_START */
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

/* LCOV_EXCL_STOP */

/*
 * Default derivatives: fourth-order finite differences of eval(), error O(h^4),
 * with h = 1e-3 (1 + z) in z and h = 1e-3 in ln k. Below z = 2h the z stencil is
 * one-sided (z to z + 4h), so it never evaluates z < 0; either stencil can reach
 * outside [zi, zf].
 */
#define NCM_POWSPEC_DERIV_REL_STEP (1.0e-3)

static gdouble
_ncm_powspec_deriv_z (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  const gdouble h = NCM_POWSPEC_DERIV_REL_STEP * (1.0 + z);

  if (z >= 2.0 * h)
  {
    const gdouble Pm2 = ncm_powspec_eval (powspec, model, z - 2.0 * h, k);
    const gdouble Pm1 = ncm_powspec_eval (powspec, model, z - h, k);
    const gdouble Pp1 = ncm_powspec_eval (powspec, model, z + h, k);
    const gdouble Pp2 = ncm_powspec_eval (powspec, model, z + 2.0 * h, k);

    return (Pm2 - 8.0 * Pm1 + 8.0 * Pp1 - Pp2) / (12.0 * h);
  }
  else
  {
    const gdouble P0 = ncm_powspec_eval (powspec, model, z, k);
    const gdouble P1 = ncm_powspec_eval (powspec, model, z + h, k);
    const gdouble P2 = ncm_powspec_eval (powspec, model, z + 2.0 * h, k);
    const gdouble P3 = ncm_powspec_eval (powspec, model, z + 3.0 * h, k);
    const gdouble P4 = ncm_powspec_eval (powspec, model, z + 4.0 * h, k);

    return (-25.0 * P0 + 48.0 * P1 - 36.0 * P2 + 16.0 * P3 - 3.0 * P4) / (12.0 * h);
  }
}

static gdouble
_ncm_powspec_deriv_k (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  const gdouble h   = NCM_POWSPEC_DERIV_REL_STEP;
  const gdouble Pm2 = ncm_powspec_eval (powspec, model, z, k * exp (-2.0 * h));
  const gdouble Pm1 = ncm_powspec_eval (powspec, model, z, k * exp (-h));
  const gdouble Pp1 = ncm_powspec_eval (powspec, model, z, k * exp (h));
  const gdouble Pp2 = ncm_powspec_eval (powspec, model, z, k * exp (2.0 * h));

  /* d/dk = (d/dln k) / k */
  return (Pm2 - 8.0 * Pm1 + 8.0 * Pp1 - Pp2) / (12.0 * h * k);
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
   * Lower end of the redshift range of $P(k, z)$.
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
   * Upper end of the redshift range of $P(k, z)$.
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
   * Lower end of the wavenumber range of $P(k, z)$, in $\mathrm{Mpc}^{-1}$.
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
   * Upper end of the wavenumber range of $P(k, z)$, in $\mathrm{Mpc}^{-1}$.
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
   * Relative tolerance of the spline built by the default get_spline_2d().
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
  klass->deriv_z       = &_ncm_powspec_deriv_z;
  klass->deriv_k       = &_ncm_powspec_deriv_k;
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
 * Returns: (transfer full): @powspec
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
 * If *@powspec is not %NULL, decrements its reference count and sets
 * *@powspec to %NULL.
 */
void
ncm_powspec_clear (NcmPowspec **powspec)
{
  g_clear_object (powspec);
}

/**
 * ncm_powspec_set_zi:
 * @powspec: a #NcmPowspec
 * @zi: lowest redshift $z_i$
 *
 * Sets NcmPowspec:zi to @zi.
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
 * @zf: highest redshift $z_f$
 *
 * Sets NcmPowspec:zf to @zf.
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
 * @kmin: lowest wavenumber $k_\mathrm{min}$ in $\mathrm{Mpc}^{-1}$
 *
 * Sets NcmPowspec:kmin to @kmin.
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
 * @kmax: highest wavenumber $k_\mathrm{max}$ in $\mathrm{Mpc}^{-1}$
 *
 * Sets NcmPowspec:kmax to @kmax.
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
 * @reltol: relative tolerance
 *
 * Sets NcmPowspec:reltol to @reltol.
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
 * @zi: lowest redshift $z_i$
 *
 * Lowers NcmPowspec:zi to @zi if @zi is below it.
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
 * @zf: highest redshift $z_f$
 *
 * Raises NcmPowspec:zf to @zf if @zf is above it.
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
 * @kmin: lowest wavenumber $k_\mathrm{min}$ in $\mathrm{Mpc}^{-1}$
 *
 * Lowers NcmPowspec:kmin to @kmin if @kmin is below it.
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
 * @kmax: highest wavenumber $k_\mathrm{max}$ in $\mathrm{Mpc}^{-1}$
 *
 * Raises NcmPowspec:kmax to @kmax if @kmax is above it.
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
 * Returns: NcmPowspec:zi
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
 * Returns: NcmPowspec:zf
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
 * Returns: NcmPowspec:kmin
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
 * Returns: NcmPowspec:kmax
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
 * Returns: NcmPowspec:reltol
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
 * Gets the number of knots of the table behind @powspec.
 */
void
ncm_powspec_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk)
{
  NCM_POWSPEC_GET_CLASS (powspec)->get_nknots (powspec, Nz, Nk);
}

/**
 * ncm_powspec_prepare:
 * @powspec: a #NcmPowspec
 * @model: (allow-none): a #NcmModel
 *
 * Prepares @powspec for @model and records the state of @model, so
 * ncm_powspec_prepare_if_needed() does nothing until @model or the range of
 * @powspec changes. With @model %NULL nothing is recorded.
 */
void
ncm_powspec_prepare (NcmPowspec *powspec, NcmModel *model)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);

  NCM_POWSPEC_GET_CLASS (powspec)->prepare (powspec, model);

  /* ncm_model_ctrl_update() dereferences its model. */
  if (model != NULL)
    ncm_model_ctrl_update (self->ctrl, model);
}

/**
 * ncm_powspec_prepare_if_needed:
 * @powspec: a #NcmPowspec
 * @model: (allow-none): a #NcmModel
 *
 * Calls ncm_powspec_prepare() if @model or the range of @powspec changed since
 * the last preparation, and always when @model is %NULL.
 */
void
ncm_powspec_prepare_if_needed (NcmPowspec *powspec, NcmModel *model)
{
  NcmPowspecPrivate * const self = ncm_powspec_get_instance_private (powspec);
  gboolean model_up;

  /* No model to compare against: prepare unconditionally. */
  if (model == NULL)
  {
    ncm_powspec_prepare (powspec, NULL);

    return;
  }

  model_up = ncm_model_ctrl_update (self->ctrl, NCM_MODEL (model));

  if (model_up)
    ncm_powspec_prepare (powspec, model);
}

/**
 * ncm_powspec_eval:
 * @powspec: a #NcmPowspec
 * @model: (allow-none): a #NcmModel
 * @z: redshift
 * @k: wavenumber in $\mathrm{Mpc}^{-1}$
 *
 * Returns: $P(k, z)$ in $\mathrm{Mpc}^3$
 */
gdouble
ncm_powspec_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  return NCM_POWSPEC_GET_CLASS (powspec)->eval (powspec, model, z, k);
}

/**
 * ncm_powspec_eval_vec:
 * @powspec: a #NcmPowspec
 * @model: (allow-none): a #NcmModel
 * @z: redshift
 * @k: wavenumbers in $\mathrm{Mpc}^{-1}$
 * @Pk: output, same length as @k
 *
 * Sets @Pk to $P(k_i, z)$ for each $k_i$ in @k.
 */
void
ncm_powspec_eval_vec (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk)
{
  NCM_POWSPEC_GET_CLASS (powspec)->eval_vec (powspec, model, z, k, Pk);
}

/**
 * ncm_powspec_deriv_z:
 * @powspec: a #NcmPowspec
 * @model: (allow-none): a #NcmModel
 * @z: redshift
 * @k: wavenumber in $\mathrm{Mpc}^{-1}$
 *
 * The default is a fourth-order finite difference of ncm_powspec_eval() with
 * step $h = 10^{-3}(1+z)$; it evaluates $P$ in $[\max(0, z - 2h), z + 4h]$,
 * which can reach outside [NcmPowspec:zi, NcmPowspec:zf].
 *
 * Returns: $\partial P(k, z) / \partial z$
 */
gdouble
ncm_powspec_deriv_z (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  return NCM_POWSPEC_GET_CLASS (powspec)->deriv_z (powspec, model, z, k);
}

/**
 * ncm_powspec_deriv_k:
 * @powspec: a #NcmPowspec
 * @model: (allow-none): a #NcmModel
 * @z: redshift
 * @k: wavenumber in $\mathrm{Mpc}^{-1}$
 *
 * The default is a fourth-order finite difference of ncm_powspec_eval() with
 * step $10^{-3}$ in $\ln k$.
 *
 * Returns: $\partial P(k, z) / \partial k$
 */
gdouble
ncm_powspec_deriv_k (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  return NCM_POWSPEC_GET_CLASS (powspec)->deriv_k (powspec, model, z, k);
}

/**
 * ncm_powspec_get_spline_2d:
 * @powspec: a #NcmPowspec
 * @model: (allow-none): a #NcmModel
 *
 * Builds a spline of $P$ in $(z, k)$. The default is a bicubic not-a-knot
 * spline with knots in $z$ and in $k$ (not $\ln k$), placed to meet
 * NcmPowspec:reltol along $z$ at the geometric middle of the $k$ range and
 * along $k$ at the middle of the $z$ range.
 *
 * Returns: (transfer full): a new #NcmSpline2d
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
 * Returns: (transfer none): the #NcmModelCtrl that records the state @powspec
 * was prepared for
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
 * @model: (allow-none): a #NcmModel
 * @reltol: relative tolerance of the quadrature
 * @z: redshift
 * @R: radius in Mpc
 *
 * Computes the variance of the field smoothed by a top-hat of radius @R,
 * $$\sigma_R^2(z) = \frac{1}{2\pi^2} \int_{k_\mathrm{min}}^{k_\mathrm{max}} W^2(kR) \, P(k, z) \, k^2 \, \mathrm{d}k, \qquad W(x) = \frac{3 j_1(x)}{x},$$
 * by quadrature in $\ln k$, after ncm_powspec_prepare_if_needed(). The integral
 * covers $[k_\mathrm{min}, k_\mathrm{max}]$ only, so for $R$ within a few
 * e-foldings of $1/k_\mathrm{max}$ it differs from #NcmPowspecFilter, which
 * continues the table into its padding. For many values of $z$ and $R$ use
 * #NcmPowspecFilter.
 *
 * Returns: $\sigma_R^2(z)$
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
 * @model: (allow-none): a #NcmModel
 * @reltol: relative tolerance of the quadrature
 * @z: redshift
 * @R: radius in Mpc
 *
 * Returns: $\sigma_R(z)$, the square root of ncm_powspec_var_tophat_R()
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
 * @model: (allow-none): a #NcmModel
 * @reltol: relative tolerance of the quadrature
 * @z: redshift
 * @r: separation in Mpc
 *
 * Computes the correlation function
 * $$\xi(r, z) = \frac{1}{2\pi^2} \int_{k_\mathrm{min}}^{k_\mathrm{max}} P(k, z) \, j_0(kr) \, k^2 \, \mathrm{d}k$$
 * by quadrature in $\ln k$, after ncm_powspec_prepare_if_needed(). For many
 * values of $r$ and $z$ use #NcmPowspecCorr3d.
 *
 * Returns: $\xi(r, z)$
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
 * @model: (allow-none): a #NcmModel
 * @reltol: relative tolerance of the quadrature
 * @ell: multipole $\ell$
 * @z1: redshift of the first sphere
 * @z2: redshift of the second sphere
 * @xi1: comoving radius of the first sphere in Mpc
 * @xi2: comoving radius of the second sphere in Mpc
 *
 * Computes the angular power spectrum of the field on two spheres,
 * $$C_\ell = \frac{2}{\pi} \int_{k_\mathrm{min}}^{k_\mathrm{max}} \sqrt{P(k, z_1) P(k, z_2)} \, j_\ell(k \xi_1) \, j_\ell(k \xi_2) \, k^2 \, \mathrm{d}k,$$
 * by quadrature in $\ln k$, after ncm_powspec_prepare_if_needed(). The
 * unequal-time spectrum is the geometric mean. Slow; for tests.
 *
 * Returns: $C_\ell$
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

