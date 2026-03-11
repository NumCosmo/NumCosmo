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
#include "math/ncm_dtuple.h"

typedef struct _NcmSBesselIntegratorPrivate
{
  guint ell_min;
  guint ell_max;
} NcmSBesselIntegratorPrivate;

enum
{
  PROP_0,
  PROP_ELL_RANGE,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmSBesselIntegrator, ncm_sbessel_integrator, G_TYPE_OBJECT)

static void
ncm_sbessel_integrator_init (NcmSBesselIntegrator *sbi)
{
  NcmSBesselIntegratorPrivate *self = ncm_sbessel_integrator_get_instance_private (sbi);

  self->ell_min = 0;
  self->ell_max = 0;
}

static void
_ncm_sbessel_integrator_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (object);

  g_return_if_fail (NCM_IS_SBESSEL_INTEGRATOR (object));

  switch (prop_id)
  {
    case PROP_ELL_RANGE:
    {
      NcmDTuple2 *ell_range = g_value_get_boxed (value);

      if (ell_range == NULL)
        g_error ("_ncm_sbessel_integrator_set_property: ell_range is NULL.");

      /* Convert from double to uint with validation */
      if ((ell_range->elements[0] < 0.0) || (ell_range->elements[1] < 0.0))
        g_error ("_ncm_sbessel_integrator_set_property: ell values must be non-negative.");

      if ((ell_range->elements[0] != floor (ell_range->elements[0])) ||
          (ell_range->elements[1] != floor (ell_range->elements[1])))
        g_error ("_ncm_sbessel_integrator_set_property: ell values must be integers.");

      if ((ell_range->elements[0] > (gdouble) G_MAXUINT) ||
          (ell_range->elements[1] > (gdouble) G_MAXUINT))
        g_error ("_ncm_sbessel_integrator_set_property: ell values out of range.");

      {
        const guint ell_min = (guint) ell_range->elements[0];
        const guint ell_max = (guint) ell_range->elements[1];

        if (ell_min > ell_max)
          g_error ("_ncm_sbessel_integrator_set_property: ell_min (%u) must be <= ell_max (%u).",
                   ell_min, ell_max);

        ncm_sbessel_integrator_set_ell_range (sbi, ell_min, ell_max);
      }
      break;
    }
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
    case PROP_ELL_RANGE:
    {
      g_value_take_boxed (value, ncm_dtuple2_new ((gdouble) self->ell_min, (gdouble) self->ell_max));
      break;
    }
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

static void _ncm_sbessel_integrator_set_ell_range_default (NcmSBesselIntegrator *sbi, guint ell_min, guint ell_max);
static gdouble _ncm_sbessel_integrator_integrate_ell_default (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, gint ell, gpointer user_data);
static void _ncm_sbessel_integrator_integrate_not_implemented (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, NcmVector *result, gpointer user_data);

static void
ncm_sbessel_integrator_class_init (NcmSBesselIntegratorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_sbessel_integrator_set_property;
  object_class->get_property = &_ncm_sbessel_integrator_get_property;
  object_class->dispose      = &_ncm_sbessel_integrator_dispose;
  object_class->finalize     = &_ncm_sbessel_integrator_finalize;

  /**
   * NcmSBesselIntegrator:ell-range:
   *
   * Multipole range [ell_min, ell_max]. Both values must be non-negative integers
   * with ell_min <= ell_max.
   */
  g_object_class_install_property (object_class,
                                   PROP_ELL_RANGE,
                                   g_param_spec_boxed ("ell-range",
                                                       NULL,
                                                       "Multipole range [ell_min, ell_max]",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->set_ell_range = &_ncm_sbessel_integrator_set_ell_range_default;
  klass->integrate_ell = &_ncm_sbessel_integrator_integrate_ell_default;
  klass->integrate     = &_ncm_sbessel_integrator_integrate_not_implemented;
}

static void
_ncm_sbessel_integrator_set_ell_range_default (NcmSBesselIntegrator *sbi, guint ell_min, guint ell_max)
{
  NcmSBesselIntegratorPrivate *self = ncm_sbessel_integrator_get_instance_private (sbi);

  /* Default implementation simply updates ell_min and ell_max */
  self->ell_min = ell_min;
  self->ell_max = ell_max;
}

static gdouble
_ncm_sbessel_integrator_integrate_ell_default (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, gint ell, gpointer user_data)
{
  NcmSBesselIntegratorPrivate *self = ncm_sbessel_integrator_get_instance_private (sbi);
  const guint old_ell_min           = self->ell_min;
  const guint old_ell_max           = self->ell_max;
  NcmVector *result                 = ncm_vector_new (1);
  gdouble val;

  /* Temporarily set range to single ell */
  ncm_sbessel_integrator_set_ell_range (sbi, ell, ell);

  /* Call vectorized integrate */
  NCM_SBESSEL_INTEGRATOR_GET_CLASS (sbi)->integrate (sbi, F, a, b, k, result, user_data);
  val = ncm_vector_get (result, 0);

  /* Restore original range */
  ncm_sbessel_integrator_set_ell_range (sbi, old_ell_min, old_ell_max);

  ncm_vector_free (result);

  return val;
}

static void
_ncm_sbessel_integrator_integrate_not_implemented (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, NcmVector *result, gpointer user_data)
{
  g_error ("ncm_sbessel_integrator_integrate: method not implemented for `%s'",
           G_OBJECT_TYPE_NAME (sbi));
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
 * ncm_sbessel_integrator_get_ell_range:
 * @sbi: a #NcmSBesselIntegrator
 * @ell_min: (out): location to store minimum multipole
 * @ell_max: (out): location to store maximum multipole
 *
 * Gets the multipole range.
 *
 */
void
ncm_sbessel_integrator_get_ell_range (NcmSBesselIntegrator *sbi, guint *ell_min, guint *ell_max)
{
  NcmSBesselIntegratorPrivate *self = ncm_sbessel_integrator_get_instance_private (sbi);

  *ell_min = self->ell_min;
  *ell_max = self->ell_max;
}

/**
 * ncm_sbessel_integrator_set_ell_range: (virtual set_ell_range)
 * @sbi: a #NcmSBesselIntegrator
 * @ell_min: minimum multipole
 * @ell_max: maximum multipole
 *
 * Sets the multipole range for integration. If the range has changed from
 * the previous call, subclasses may perform preparation work (e.g., allocating
 * operators for the new range). The default implementation simply updates
 * ell_min and ell_max properties.
 *
 */
void
ncm_sbessel_integrator_set_ell_range (NcmSBesselIntegrator *sbi, guint ell_min, guint ell_max)
{
  NCM_SBESSEL_INTEGRATOR_GET_CLASS (sbi)->set_ell_range (sbi, ell_min, ell_max);
}

/**
 * ncm_sbessel_integrator_integrate_ell: (virtual integrate_ell)
 * @sbi: a #NcmSBesselIntegrator
 * @F: (scope call) (closure user_data): function to integrate
 * @a: lower integration limit
 * @b: upper integration limit
 * @k: wave number parameter
 * @ell: multipole
 * @user_data: (nullable): user data passed to @F
 *
 * Integrates the function @F(x, k) multiplied by the spherical Bessel function
 * $j_\ell(kx)$ from @a to @b for a single multipole.
 * Computes: $\int_a^b K(x,k) j_\ell(kx) dx$
 *
 * Returns: the integral value
 */
gdouble
ncm_sbessel_integrator_integrate_ell (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, gint ell, gpointer user_data)
{
  return NCM_SBESSEL_INTEGRATOR_GET_CLASS (sbi)->integrate_ell (sbi, F, a, b, k, ell, user_data);
}

/**
 * ncm_sbessel_integrator_integrate: (virtual integrate)
 * @sbi: a #NcmSBesselIntegrator
 * @F: (scope call) (closure user_data): function to integrate
 * @a: lower integration limit
 * @b: upper integration limit
 * @k: wave number parameter
 * @result: a #NcmVector to store results
 * @user_data: (nullable): user data passed to @F
 *
 * Integrates the function @F(x, k) multiplied by the spherical Bessel function
 * $j_\ell(kx)$ from @a to @b for all multipoles from ell_min to ell_max.
 * Computes: $\int_a^b K(x,k) j_\ell(kx) dx$ for each $\ell$.
 * The results are stored in @result, which must have length (ell_max - ell_min + 1).
 *
 */
void
ncm_sbessel_integrator_integrate (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, NcmVector *result, gpointer user_data)
{
  NCM_SBESSEL_INTEGRATOR_GET_CLASS (sbi)->integrate (sbi, F, a, b, k, result, user_data);
}

typedef struct _NcmSBesselIntegratorGaussianData
{
  gdouble center;
  gdouble std;
  gdouble k;
} NcmSBesselIntegratorGaussianData;

static gdouble
_ncm_sbessel_integrator_gaussian_func (gpointer user_data, gdouble x, gdouble k)
{
  NcmSBesselIntegratorGaussianData *data = (NcmSBesselIntegratorGaussianData *) user_data;
  const gdouble z                        = (x - data->center) / data->std;

  return exp (-0.5 * z * z);
}

/**
 * ncm_sbessel_integrator_integrate_gaussian_ell:
 * @sbi: a #NcmSBesselIntegrator
 * @center: center of the Gaussian
 * @std: standard deviation of the Gaussian
 * @k: wave number parameter
 * @a: lower integration limit
 * @b: upper integration limit
 * @ell: multipole
 *
 * Integrates a Gaussian function $\exp(-\frac{1}{2}(\frac{x - center}{std})^2)$
 * multiplied by the spherical Bessel function $j_\ell(kx)$ from @a to @b
 * for a single multipole.
 *
 * This is a convenience function optimized for testing against truth tables,
 * avoiding the overhead of Python callbacks.
 *
 * Returns: the integral value
 */
gdouble
ncm_sbessel_integrator_integrate_gaussian_ell (NcmSBesselIntegrator *sbi, gdouble center, gdouble std, gdouble k, gdouble a, gdouble b, gint ell)
{
  NcmSBesselIntegratorGaussianData data = {center, std, k};

  return ncm_sbessel_integrator_integrate_ell (sbi, &_ncm_sbessel_integrator_gaussian_func, a, b, k, ell, &data);
}

/**
 * ncm_sbessel_integrator_integrate_gaussian:
 * @sbi: a #NcmSBesselIntegrator
 * @center: center of the Gaussian
 * @std: standard deviation of the Gaussian
 * @k: wave number parameter
 * @a: lower integration limit
 * @b: upper integration limit
 * @result: a #NcmVector to store results
 *
 * Integrates a Gaussian function $\exp(-\frac{1}{2}(\frac{x - center}{std})^2)$
 * multiplied by the spherical Bessel function $j_\ell(kx)$ from @a to @b
 * for all multipoles from ell_min to ell_max.
 * The results are stored in @result, which must have length (ell_max - ell_min + 1).
 *
 * This is a convenience function optimized for testing against truth tables,
 * avoiding the overhead of Python callbacks.
 *
 */
void
ncm_sbessel_integrator_integrate_gaussian (NcmSBesselIntegrator *sbi, gdouble center, gdouble std, gdouble k, gdouble a, gdouble b, NcmVector *result)
{
  NcmSBesselIntegratorGaussianData data = {center, std, k};

  ncm_sbessel_integrator_integrate (sbi, &_ncm_sbessel_integrator_gaussian_func, a, b, k, result, &data);
}

typedef struct _NcmSBesselIntegratorRationalData
{
  gdouble center;
  gdouble std;
  gdouble k;
} NcmSBesselIntegratorRationalData;

static gdouble
_ncm_sbessel_integrator_rational_func (gpointer user_data, gdouble x, gdouble k)
{
  NcmSBesselIntegratorRationalData *data = (NcmSBesselIntegratorRationalData *) user_data;
  const gdouble z                        = (x - data->center) / data->std;
  const gdouble denom                    = 1.0 + z * z;
  const gdouble denom_cubed              = denom * denom * denom;

  return x * x / denom_cubed;
}

/**
 * ncm_sbessel_integrator_integrate_rational_ell:
 * @sbi: a #NcmSBesselIntegrator
 * @center: center of the rational function
 * @std: standard deviation parameter
 * @k: wave number parameter
 * @a: lower integration limit
 * @b: upper integration limit
 * @ell: multipole
 *
 * Integrates a rational function $\frac{x^2}{(1+((x - center)/std)^2)^3}$
 * multiplied by the spherical Bessel function $j_\ell(kx)$ from @a to @b for a single
 * multipole.
 *
 * This is a convenience function optimized for testing against truth tables, avoiding
 * the overhead of Python callbacks.
 *
 * Returns: the integral value
 */
gdouble
ncm_sbessel_integrator_integrate_rational_ell (NcmSBesselIntegrator *sbi, gdouble center, gdouble std, gdouble k, gdouble a, gdouble b, gint ell)
{
  NcmSBesselIntegratorRationalData data = {center, std, k};

  return ncm_sbessel_integrator_integrate_ell (sbi, &_ncm_sbessel_integrator_rational_func, a, b, k, ell, &data);
}

/**
 * ncm_sbessel_integrator_integrate_rational:
 * @sbi: a #NcmSBesselIntegrator
 * @center: center of the rational function
 * @std: standard deviation parameter
 * @k: wave number parameter
 * @a: lower integration limit
 * @b: upper integration limit
 * @result: a #NcmVector to store results
 *
 * Integrates a rational function $\frac{x^2}{(1+((x - center)/std)^2)^3}$
 * multiplied by the spherical Bessel function $j_\ell(kx)$ from @a to @b for all
 * multipoles from ell_min to ell_max. The results are stored in @result, which must have
 * length (ell_max - ell_min + 1).
 *
 * This is a convenience function optimized for testing against truth tables, avoiding
 * the overhead of Python callbacks.
 *
 */
void
ncm_sbessel_integrator_integrate_rational (NcmSBesselIntegrator *sbi, gdouble center, gdouble std, gdouble k, gdouble a, gdouble b, NcmVector *result)
{
  NcmSBesselIntegratorRationalData data = {center, std, k};

  ncm_sbessel_integrator_integrate (sbi, &_ncm_sbessel_integrator_rational_func, a, b, k, result, &data);
}

