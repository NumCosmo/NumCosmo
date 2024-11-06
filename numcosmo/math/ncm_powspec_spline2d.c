/***************************************************************************
 *            ncm_powspec_spline2d.c
 *
 *  Tue February 16 17:00:52 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_powspec_spline2d.c
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
 * SECTION:ncm_powspec_spline2d
 * @title: NcmPowspecSpline2d
 * @short_description: Power spectrum implementation using a 2D spline
 * @stability: Stable
 * @include: numcosmo/math/ncm_powspec_spline2d.h
 *
 * #NcmPowspecSpline2d is a power spectrum implementation using a 2D spline.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_powspec_spline2d.h"
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

typedef struct _NcmPowspecSpline2dPrivate
{
  /*< private > */
  NcmSpline2d *spline2d;
  gdouble intern_lnkmin;
  gdouble intern_lnkmax;
  gdouble intern_lnkmax_m1;
} NcmPowspecSpline2dPrivate;

struct _NcmPowspecSpline2d
{
  NcmPowspec parent_instance;
};

enum
{
  PROP_0,
  PROP_SPLINE2D,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmPowspecSpline2d, ncm_powspec_spline2d, NCM_TYPE_POWSPEC)

static void
ncm_powspec_spline2d_init (NcmPowspecSpline2d *ps_s2d)
{
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);

  self->spline2d         = NULL;
  self->intern_lnkmin    = 0.0;
  self->intern_lnkmax    = 0.0;
  self->intern_lnkmax_m1 = 0.0;
}

static void
_ncm_powspec_spline2d_dispose (GObject *object)
{
  NcmPowspecSpline2d *ps_s2d             = NCM_POWSPEC_SPLINE2D (object);
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);

  ncm_spline2d_clear (&self->spline2d);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_spline2d_parent_class)->dispose (object);
}

static void
_ncm_powspec_spline2d_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_spline2d_parent_class)->finalize (object);
}

static void
_ncm_powspec_spline2d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPowspecSpline2d *ps_s2d = NCM_POWSPEC_SPLINE2D (object);

  g_return_if_fail (NCM_IS_POWSPEC_SPLINE2D (object));

  switch (prop_id)
  {
    case PROP_SPLINE2D:
      ncm_powspec_spline2d_set_spline2d (ps_s2d, g_value_get_object (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_powspec_spline2d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPowspecSpline2d *ps_s2d = NCM_POWSPEC_SPLINE2D (object);

  g_return_if_fail (NCM_IS_POWSPEC_SPLINE2D (object));

  switch (prop_id)
  {
    case PROP_SPLINE2D:
      g_value_set_object (value, ncm_powspec_spline2d_peek_spline2d (ps_s2d));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_powspec_spline2d_prepare (NcmPowspec *powspec, NcmModel *model)
{
  NcmPowspecSpline2d *ps_s2d             = NCM_POWSPEC_SPLINE2D (powspec);
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);

  if (!ncm_spline2d_is_init (self->spline2d))
    ncm_spline2d_prepare (self->spline2d);
}

static gdouble
_ncm_powspec_spline2d_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  NcmPowspecSpline2d *ps_s2d             = NCM_POWSPEC_SPLINE2D (powspec);
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);
  const gdouble lnk                      = log (k);

  if (lnk < self->intern_lnkmin)
  {
    const gdouble lnPkmin = ncm_spline2d_eval (self->spline2d, z, self->intern_lnkmin);

    return exp (lnPkmin + 3.0 * (lnk - self->intern_lnkmin));
  }
  else if (lnk > self->intern_lnkmax)
  {
    const gdouble lnkmax     = self->intern_lnkmax;
    const gdouble lnPkmax    = ncm_spline2d_eval (self->spline2d, z, self->intern_lnkmax);
    const gdouble lnPkmax_m1 = ncm_spline2d_eval (self->spline2d, z, self->intern_lnkmax_m1);
    const gdouble delta_lnk  = lnk - lnkmax;
    const gdouble lambda     = (lnPkmax - lnPkmax_m1) / (self->intern_lnkmax - self->intern_lnkmax_m1);

    return exp (lnPkmax - 0.5 * 10.0 * gsl_pow_2 (delta_lnk) + lambda * delta_lnk);
  }
  else
  {
    return exp (ncm_spline2d_eval (self->spline2d, z, lnk));
  }
}

static void
_ncm_powspec_spline2d_eval_vec (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk)
{
  guint n = ncm_vector_len (k);
  guint i;

  g_assert_cmpuint (n, ==, ncm_vector_len (Pk));

  for (i = 0; i < n; i++)
  {
    const gdouble k_i = ncm_vector_get (k, i);

    ncm_vector_set (Pk, i, _ncm_powspec_spline2d_eval (powspec, model, z, k_i));
  }
}

void
_ncm_powspec_spline2d_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk)
{
  NcmPowspecSpline2d *ps_s2d             = NCM_POWSPEC_SPLINE2D (powspec);
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);

  Nz[0] = ncm_vector_len (ncm_spline2d_peek_xv (self->spline2d));
  Nk[0] = ncm_vector_len (ncm_spline2d_peek_yv (self->spline2d));
}

static NcmSpline2d *
_ncm_powspec_spline2d_get_spline_2d (NcmPowspec *powspec, NcmModel *model)
{
  NcmPowspecSpline2d *ps_s2d             = NCM_POWSPEC_SPLINE2D (powspec);
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);

  return ncm_spline2d_ref (self->spline2d);
}

static void
ncm_powspec_spline2d_class_init (NcmPowspecSpline2dClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmPowspecClass *powspec_class = NCM_POWSPEC_CLASS (klass);

  object_class->set_property = &_ncm_powspec_spline2d_set_property;
  object_class->get_property = &_ncm_powspec_spline2d_get_property;
  object_class->dispose      = &_ncm_powspec_spline2d_dispose;
  object_class->finalize     = &_ncm_powspec_spline2d_finalize;

  /**
   * NcmPowspecSpline2d:reltol:
   *
   * The relative tolerance on the interpolation error.
   */
  g_object_class_install_property (object_class,
                                   PROP_SPLINE2D,
                                   g_param_spec_object ("spline2d",
                                                        NULL,
                                                        "Spline2d representing the values of the power-spectrum",
                                                        NCM_TYPE_SPLINE2D,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  powspec_class->prepare       = &_ncm_powspec_spline2d_prepare;
  powspec_class->eval          = &_ncm_powspec_spline2d_eval;
  powspec_class->eval_vec      = &_ncm_powspec_spline2d_eval_vec;
  powspec_class->get_nknots    = &_ncm_powspec_spline2d_get_nknots;
  powspec_class->get_spline_2d = &_ncm_powspec_spline2d_get_spline_2d;
}

NcmPowspecSpline2d *
ncm_powspec_spline2d_new (NcmSpline2d *spline2d)
{
  NcmPowspecSpline2d *ps_s2d = g_object_new (NCM_TYPE_POWSPEC_SPLINE2D,
                                             "spline2d", spline2d,
                                             NULL);

  return ps_s2d;
}

/**
 * ncm_powspec_spline2d_ref:
 * @ps_s2d: a #NcmPowspecSpline2d
 *
 * Increases the reference count of @ps_s2d by one atomically.
 *
 * Returns: (transfer full): @ps_s2d.
 */
NcmPowspecSpline2d *
ncm_powspec_spline2d_ref (NcmPowspecSpline2d *ps_s2d)
{
  return g_object_ref (ps_s2d);
}

/**
 * ncm_powspec_spline2d_free:
 * @ps_s2d: a #NcmPowspecSpline2d
 *
 * Atomically decrements the reference count of @ps_s2d by one.
 * If the reference count drops to 0,
 * all memory allocated by @ps_s2d is released.
 *
 */
void
ncm_powspec_spline2d_free (NcmPowspecSpline2d *ps_s2d)
{
  g_object_unref (ps_s2d);
}

/**
 * ncm_powspec_spline2d_clear:
 * @ps_s2d: a #NcmPowspecSpline2d
 *
 * If @ps_s2d is different from NULL,
 * atomically decrements the reference count of @powspec by one.
 * If the reference count drops to 0,
 * all memory allocated by @powspec is released and @powspec is set to NULL.
 *
 */
void
ncm_powspec_spline2d_clear (NcmPowspecSpline2d **ps_s2d)
{
  g_clear_object (ps_s2d);
}

/**
 * ncm_powspec_spline2d_set_spline2d:
 * @ps_s2d: a #NcmPowspecSpline2d
 * @spline2d: a NcmSpline2d
 *
 * Sets the #NcmSpline2d to @spline2d.
 *
 */
void
ncm_powspec_spline2d_set_spline2d (NcmPowspecSpline2d *ps_s2d, NcmSpline2d *spline2d)
{
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);

  g_assert_nonnull (spline2d);

  ncm_spline2d_clear (&self->spline2d);

  self->spline2d = ncm_spline2d_ref (spline2d);

  {
    NcmVector *z_vec        = ncm_spline2d_peek_xv (self->spline2d);
    NcmVector *lnk_vec      = ncm_spline2d_peek_yv (self->spline2d);
    const gdouble lnkmin    = ncm_vector_get (lnk_vec, 0);
    const gdouble lnkmax    = ncm_vector_get (lnk_vec, ncm_vector_len (lnk_vec) - 1);
    const gdouble lnkmax_m1 = lnkmax - 1.0;

    self->intern_lnkmin    = lnkmin;
    self->intern_lnkmax    = lnkmax;
    self->intern_lnkmax_m1 = lnkmax_m1;

    ncm_powspec_set_kmin (NCM_POWSPEC (ps_s2d), exp (lnkmin));
    ncm_powspec_set_kmax (NCM_POWSPEC (ps_s2d), exp (lnkmax));

    ncm_powspec_set_zi (NCM_POWSPEC (ps_s2d), ncm_vector_get (z_vec, 0));
    ncm_powspec_set_zf (NCM_POWSPEC (ps_s2d), ncm_vector_get (z_vec, ncm_vector_len (z_vec) - 1));
  }
}

/**
 * ncm_powspec_spline2d_peek_spline2d:
 * @ps_s2d: a #NcmPowspecSpline2d
 *
 * Peeks the current #NcmSpline2d.
 *
 * Returns: (transfer none): the current #NcmSpline2d.
 */
NcmSpline2d *
ncm_powspec_spline2d_peek_spline2d (NcmPowspecSpline2d *ps_s2d)
{
  NcmPowspecSpline2dPrivate * const self = ncm_powspec_spline2d_get_instance_private (ps_s2d);

  return self->spline2d;
}

