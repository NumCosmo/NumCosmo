/***************************************************************************
 *            ncm_sbessel_integrator.c
 *
 *  Thu January 09 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sbessel_integrator.c
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
 * NcmSBesselIntegrator:
 *
 * Base class for spherical Bessel function integrators.
 *
 * This class provides a framework for integrating functions multiplied by
 * spherical Bessel functions $j_\ell(x)$.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sbessel_integrator.h"

typedef struct _NcmSBesselIntegratorPrivate
{
  guint lmin;
  guint lmax;
} NcmSBesselIntegratorPrivate;

enum
{
  PROP_0,
  PROP_LMIN,
  PROP_LMAX,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmSBesselIntegrator, ncm_sbessel_integrator, G_TYPE_OBJECT)

static void
ncm_sbessel_integrator_init (NcmSBesselIntegrator *sbi)
{
  NcmSBesselIntegratorPrivate *self = ncm_sbessel_integrator_get_instance_private (sbi);

  self->lmin = 0;
  self->lmax = 0;
}

static void
_ncm_sbessel_integrator_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSBesselIntegrator *sbi         = NCM_SBESSEL_INTEGRATOR (object);
  NcmSBesselIntegratorPrivate *self = ncm_sbessel_integrator_get_instance_private (sbi);

  g_return_if_fail (NCM_IS_SBESSEL_INTEGRATOR (object));

  switch (prop_id)
  {
    case PROP_LMIN:
      self->lmin = g_value_get_uint (value);
      break;
    case PROP_LMAX:
      self->lmax = g_value_get_uint (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sbessel_integrator_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSBesselIntegrator *sbi         = NCM_SBESSEL_INTEGRATOR (object);
  NcmSBesselIntegratorPrivate *self = ncm_sbessel_integrator_get_instance_private (sbi);

  g_return_if_fail (NCM_IS_SBESSEL_INTEGRATOR (object));

  switch (prop_id)
  {
    case PROP_LMIN:
      g_value_set_uint (value, self->lmin);
      break;
    case PROP_LMAX:
      g_value_set_uint (value, self->lmax);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sbessel_integrator_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sbessel_integrator_parent_class)->dispose (object);
}

static void
_ncm_sbessel_integrator_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sbessel_integrator_parent_class)->finalize (object);
}

static void _ncm_sbessel_integrator_prepare_not_implemented (NcmSBesselIntegrator *sbi);
static gdouble _ncm_sbessel_integrator_integrate_ell_not_implemented (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gpointer user_data, gdouble a, gdouble b, gint ell);
static void _ncm_sbessel_integrator_integrate_default (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gpointer user_data, gdouble a, gdouble b, gdouble *result);

static void
ncm_sbessel_integrator_class_init (NcmSBesselIntegratorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_sbessel_integrator_set_property;
  object_class->get_property = &_ncm_sbessel_integrator_get_property;
  object_class->dispose      = &_ncm_sbessel_integrator_dispose;
  object_class->finalize     = &_ncm_sbessel_integrator_finalize;

  /**
   * NcmSBesselIntegrator:lmin:
   *
   * Minimum multipole.
   */
  g_object_class_install_property (object_class,
                                   PROP_LMIN,
                                   g_param_spec_uint ("lmin",
                                                      NULL,
                                                      "Minimum multipole",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSBesselIntegrator:lmax:
   *
   * Maximum multipole.
   */
  g_object_class_install_property (object_class,
                                   PROP_LMAX,
                                   g_param_spec_uint ("lmax",
                                                      NULL,
                                                      "Maximum multipole",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->prepare       = &_ncm_sbessel_integrator_prepare_not_implemented;
  klass->integrate_ell = &_ncm_sbessel_integrator_integrate_ell_not_implemented;
  klass->integrate     = &_ncm_sbessel_integrator_integrate_default;
}

static void
_ncm_sbessel_integrator_prepare_not_implemented (NcmSBesselIntegrator *sbi)
{
  g_error ("ncm_sbessel_integrator_prepare: method not implemented for `%s'",
           G_OBJECT_TYPE_NAME (sbi));
}

static gdouble
_ncm_sbessel_integrator_integrate_ell_not_implemented (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gpointer user_data, gdouble a, gdouble b, gint ell)
{
  g_error ("ncm_sbessel_integrator_integrate_ell: method not implemented for `%s'",
           G_OBJECT_TYPE_NAME (sbi));

  return 0.0;
}

static void
_ncm_sbessel_integrator_integrate_default (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gpointer user_data, gdouble a, gdouble b, gdouble *result)
{
  NcmSBesselIntegratorPrivate *self = ncm_sbessel_integrator_get_instance_private (sbi);
  gint ell;

  for (ell = self->lmin; ell <= (gint) self->lmax; ell++)
  {
    result[ell - self->lmin] = NCM_SBESSEL_INTEGRATOR_GET_CLASS (sbi)->integrate_ell (sbi, F, user_data, a, b, ell);
  }
}

/**
 * ncm_sbessel_integrator_ref:
 * @sbi: a #NcmSBesselIntegrator
 *
 * Increases the reference count of @sbi by one.
 *
 * Returns: (transfer full): @sbi
 */
NcmSBesselIntegrator *
ncm_sbessel_integrator_ref (NcmSBesselIntegrator *sbi)
{
  return g_object_ref (sbi);
}

/**
 * ncm_sbessel_integrator_free:
 * @sbi: a #NcmSBesselIntegrator
 *
 * Decreases the reference count of @sbi by one.
 *
 */
void
ncm_sbessel_integrator_free (NcmSBesselIntegrator *sbi)
{
  g_object_unref (sbi);
}

/**
 * ncm_sbessel_integrator_clear:
 * @sbi: a #NcmSBesselIntegrator
 *
 * If @sbi is different from NULL, decreases the reference count of
 * @sbi by one and sets @sbi to NULL.
 *
 */
void
ncm_sbessel_integrator_clear (NcmSBesselIntegrator **sbi)
{
  g_clear_object (sbi);
}

/**
 * ncm_sbessel_integrator_get_lmin:
 * @sbi: a #NcmSBesselIntegrator
 *
 * Gets the minimum multipole.
 *
 * Returns: the minimum multipole
 */
guint
ncm_sbessel_integrator_get_lmin (NcmSBesselIntegrator *sbi)
{
  NcmSBesselIntegratorPrivate *self = ncm_sbessel_integrator_get_instance_private (sbi);

  return self->lmin;
}

/**
 * ncm_sbessel_integrator_get_lmax:
 * @sbi: a #NcmSBesselIntegrator
 *
 * Gets the maximum multipole.
 *
 * Returns: the maximum multipole
 */
guint
ncm_sbessel_integrator_get_lmax (NcmSBesselIntegrator *sbi)
{
  NcmSBesselIntegratorPrivate *self = ncm_sbessel_integrator_get_instance_private (sbi);

  return self->lmax;
}

/**
 * ncm_sbessel_integrator_set_lmin:
 * @sbi: a #NcmSBesselIntegrator
 * @lmin: minimum multipole
 *
 * Sets the minimum multipole.
 *
 */
void
ncm_sbessel_integrator_set_lmin (NcmSBesselIntegrator *sbi, guint lmin)
{
  NcmSBesselIntegratorPrivate *self = ncm_sbessel_integrator_get_instance_private (sbi);

  self->lmin = lmin;
}

/**
 * ncm_sbessel_integrator_set_lmax:
 * @sbi: a #NcmSBesselIntegrator
 * @lmax: maximum multipole
 *
 * Sets the maximum multipole.
 *
 */
void
ncm_sbessel_integrator_set_lmax (NcmSBesselIntegrator *sbi, guint lmax)
{
  NcmSBesselIntegratorPrivate *self = ncm_sbessel_integrator_get_instance_private (sbi);

  self->lmax = lmax;
}

/**
 * ncm_sbessel_integrator_prepare: (virtual prepare)
 * @sbi: a #NcmSBesselIntegrator
 *
 * Prepares the integrator for integration.
 *
 */
void
ncm_sbessel_integrator_prepare (NcmSBesselIntegrator *sbi)
{
  NCM_SBESSEL_INTEGRATOR_GET_CLASS (sbi)->prepare (sbi);
}

/**
 * ncm_sbessel_integrator_integrate_ell: (virtual integrate_ell)
 * @sbi: a #NcmSBesselIntegrator
 * @F: (scope call): function to integrate
 * @user_data: user data passed to @F
 * @a: lower integration limit
 * @b: upper integration limit
 * @ell: multipole
 *
 * Integrates the function @F multiplied by the spherical Bessel function
 * $j_\ell(x)$ from @a to @b for a single multipole.
 *
 * Returns: the integral value
 */
gdouble
ncm_sbessel_integrator_integrate_ell (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gpointer user_data, gdouble a, gdouble b, gint ell)
{
  return NCM_SBESSEL_INTEGRATOR_GET_CLASS (sbi)->integrate_ell (sbi, F, user_data, a, b, ell);
}

/**
 * ncm_sbessel_integrator_integrate: (virtual integrate)
 * @sbi: a #NcmSBesselIntegrator
 * @F: (scope call): function to integrate
 * @user_data: user data passed to @F
 * @a: lower integration limit
 * @b: upper integration limit
 * @result: (array) (out caller-allocates): array to store results
 *
 * Integrates the function @F multiplied by the spherical Bessel function
 * $j_\ell(x)$ from @a to @b for all multipoles from lmin to lmax.
 * The results are stored in @result, which must be pre-allocated with
 * size at least (lmax - lmin + 1).
 *
 */
void
ncm_sbessel_integrator_integrate (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gpointer user_data, gdouble a, gdouble b, gdouble *result)
{
  NCM_SBESSEL_INTEGRATOR_GET_CLASS (sbi)->integrate (sbi, F, user_data, a, b, result);
}

