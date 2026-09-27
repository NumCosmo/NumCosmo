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
 * Base class for spherical Bessel integrators.
 *
 * Computes
 * $$
 * I_\ell(k) = \int_a^b F(\chi, k)\, j_\ell(k\chi)\, \mathrm{d}\chi
 * $$
 * for every multipole in #NcmSBesselIntegrator:ell-range at once, with
 * ncm_sbessel_integrator_integrate(), or for a single one with
 * ncm_sbessel_integrator_integrate_ell(). ncm_sbessel_integrator_integrate_deriv()
 * replaces $j_\ell$ by one of its first two derivatives.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/specfunc/ncm_sbessel_integrator.h"
#include "ncm/core/ncm_dtuple.h"

typedef struct _NcmSBesselIntegratorPrivate
{
  guint ell_min;
  guint ell_max;
} NcmSBesselIntegratorPrivate;

enum
{
  PROP_0,
  PROP_ELL_RANGE,
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

      /* Not given at construction: the range [0, 0] */
      if (ell_range == NULL)
      {
        ncm_sbessel_integrator_set_ell_range (sbi, 0, 0);
        break;
      }

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
static void _ncm_sbessel_integrator_integrate_deriv_not_implemented (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, guint deriv, NcmVector *result, gpointer user_data);

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
   * Multipole range $[\ell_\mathrm{min}, \ell_\mathrm{max}]$, two non-negative integers
   * in order. Defaults to $[0, 0]$.
   */
  g_object_class_install_property (object_class,
                                   PROP_ELL_RANGE,
                                   g_param_spec_boxed ("ell-range",
                                                       NULL,
                                                       "Multipole range [ell_min, ell_max]",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->set_ell_range   = &_ncm_sbessel_integrator_set_ell_range_default;
  klass->integrate_ell   = &_ncm_sbessel_integrator_integrate_ell_default;
  klass->integrate       = &_ncm_sbessel_integrator_integrate_not_implemented;
  klass->integrate_deriv = &_ncm_sbessel_integrator_integrate_deriv_not_implemented;
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

/* LCOV_EXCL_START */

static void
_ncm_sbessel_integrator_integrate_not_implemented (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, NcmVector *result, gpointer user_data)
{
  g_error ("ncm_sbessel_integrator_integrate: method not implemented for `%s'",
           G_OBJECT_TYPE_NAME (sbi));
}

static void
_ncm_sbessel_integrator_integrate_deriv_not_implemented (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, guint deriv, NcmVector *result, gpointer user_data)
{
  g_error ("ncm_sbessel_integrator_integrate_deriv: method not implemented for `%s'",
           G_OBJECT_TYPE_NAME (sbi));
}

/* LCOV_EXCL_STOP */

static void
_ncm_sbessel_integrator_check_result (NcmSBesselIntegrator *sbi, NcmVector *result)
{
  NcmSBesselIntegratorPrivate *self = ncm_sbessel_integrator_get_instance_private (sbi);
  const guint n_ell                 = self->ell_max - self->ell_min + 1;

  if (ncm_vector_len (result) != n_ell)
    g_error ("ncm_sbessel_integrator: result has length %u, the multipole range [%u, %u] needs %u.",
             ncm_vector_len (result), self->ell_min, self->ell_max, n_ell);
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
 * If *@sbi is not NULL, decreases its reference count by one and sets *@sbi to NULL.
 */
void
ncm_sbessel_integrator_clear (NcmSBesselIntegrator **sbi)
{
  g_clear_object (sbi);
}

/**
 * ncm_sbessel_integrator_get_ell_range:
 * @sbi: a #NcmSBesselIntegrator
 * @ell_min: (out): lowest multipole
 * @ell_max: (out): highest multipole
 *
 * Gets #NcmSBesselIntegrator:ell-range.
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
 * @ell_min: lowest multipole
 * @ell_max: highest multipole
 *
 * Sets #NcmSBesselIntegrator:ell-range. Subclasses may prepare work for the new range
 * here, such as allocating operators.
 */
void
ncm_sbessel_integrator_set_ell_range (NcmSBesselIntegrator *sbi, guint ell_min, guint ell_max)
{
  NCM_SBESSEL_INTEGRATOR_GET_CLASS (sbi)->set_ell_range (sbi, ell_min, ell_max);
}

/**
 * ncm_sbessel_integrator_integrate_ell: (virtual integrate_ell)
 * @sbi: a #NcmSBesselIntegrator
 * @F: (scope call) (closure user_data): the function $F(\chi, k)$
 * @a: lower limit
 * @b: upper limit
 * @k: wavenumber $k$
 * @ell: multipole $\ell \geq 0$
 * @user_data: (nullable): user data passed to @F
 *
 * Computes $I_\ell(k)$ for a single multipole. The default implementation sets the
 * range to $[\ell, \ell]$, calls ncm_sbessel_integrator_integrate() and restores the
 * range, so it is not reentrant and may trigger the preparation work of
 * ncm_sbessel_integrator_set_ell_range() twice per call.
 *
 * Returns: $I_\ell(k)$
 */
gdouble
ncm_sbessel_integrator_integrate_ell (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, gint ell, gpointer user_data)
{
  if (ell < 0)
    g_error ("ncm_sbessel_integrator_integrate_ell: negative multipole %d.", ell);

  return NCM_SBESSEL_INTEGRATOR_GET_CLASS (sbi)->integrate_ell (sbi, F, a, b, k, ell, user_data);
}

/**
 * ncm_sbessel_integrator_integrate: (virtual integrate)
 * @sbi: a #NcmSBesselIntegrator
 * @F: (scope call) (closure user_data): the function $F(\chi, k)$
 * @a: lower limit
 * @b: upper limit
 * @k: wavenumber $k$
 * @result: output, one value per multipole
 * @user_data: (nullable): user data passed to @F
 *
 * Computes $I_\ell(k)$ for every multipole in #NcmSBesselIntegrator:ell-range, storing
 * $I_{\ell_\mathrm{min}+i}(k)$ in element $i$ of @result, whose length must be
 * $\ell_\mathrm{max} - \ell_\mathrm{min} + 1$.
 */
void
ncm_sbessel_integrator_integrate (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, NcmVector *result, gpointer user_data)
{
  _ncm_sbessel_integrator_check_result (sbi, result);
  NCM_SBESSEL_INTEGRATOR_GET_CLASS (sbi)->integrate (sbi, F, a, b, k, result, user_data);
}

/**
 * ncm_sbessel_integrator_integrate_deriv: (virtual integrate_deriv)
 * @sbi: a #NcmSBesselIntegrator
 * @F: (scope call) (closure user_data): the function $F(\chi, k)$
 * @a: lower limit
 * @b: upper limit
 * @k: wavenumber $k$
 * @deriv: derivative order $d \leq 2$
 * @result: output, one value per multipole
 * @user_data: (nullable): user data passed to @F
 *
 * Same as ncm_sbessel_integrator_integrate() with $j_\ell(k\chi)$ replaced by
 * $j_\ell^{(d)}(k\chi)$, the derivative with respect to the argument. Order zero is
 * ncm_sbessel_integrator_integrate().
 */
void
ncm_sbessel_integrator_integrate_deriv (NcmSBesselIntegrator *sbi, NcmSBesselIntegratorF F, gdouble a, gdouble b, gdouble k, guint deriv, NcmVector *result, gpointer user_data)
{
  _ncm_sbessel_integrator_check_result (sbi, result);

  if (deriv > 2)
    g_error ("ncm_sbessel_integrator_integrate_deriv: derivative order %u not supported (up to 2)", deriv);

  if (deriv == 0)
    NCM_SBESSEL_INTEGRATOR_GET_CLASS (sbi)->integrate (sbi, F, a, b, k, result, user_data);
  else
    NCM_SBESSEL_INTEGRATOR_GET_CLASS (sbi)->integrate_deriv (sbi, F, a, b, k, deriv, result, user_data);
}

typedef struct _NcmSBesselIntegratorGaussianData
{
  gdouble center;
  gdouble std;
} NcmSBesselIntegratorGaussianData;

static gdouble
_ncm_sbessel_integrator_gaussian_func (gpointer user_data, gdouble chi, gdouble k)
{
  NcmSBesselIntegratorGaussianData *data = (NcmSBesselIntegratorGaussianData *) user_data;
  const gdouble z                        = (chi - data->center) / data->std;

  return exp (-0.5 * z * z);
}

/**
 * ncm_sbessel_integrator_integrate_gaussian_ell:
 * @sbi: a #NcmSBesselIntegrator
 * @center: center $\chi_c$
 * @std: width $\sigma$
 * @a: lower limit
 * @b: upper limit
 * @k: wavenumber $k$
 * @ell: multipole $\ell \geq 0$
 *
 * ncm_sbessel_integrator_integrate_ell() with $F = e^{-(\chi - \chi_c)^2 / (2\sigma^2)}$,
 * the shape of the Gaussian truth tables, evaluated in C to avoid callback overhead.
 *
 * Returns: $I_\ell(k)$
 */
gdouble
ncm_sbessel_integrator_integrate_gaussian_ell (NcmSBesselIntegrator *sbi, gdouble center, gdouble std, gdouble a, gdouble b, gdouble k, gint ell)
{
  NcmSBesselIntegratorGaussianData data = {center, std};

  return ncm_sbessel_integrator_integrate_ell (sbi, &_ncm_sbessel_integrator_gaussian_func, a, b, k, ell, &data);
}

/**
 * ncm_sbessel_integrator_integrate_gaussian:
 * @sbi: a #NcmSBesselIntegrator
 * @center: center $\chi_c$
 * @std: width $\sigma$
 * @a: lower limit
 * @b: upper limit
 * @k: wavenumber $k$
 * @result: output, one value per multipole
 *
 * ncm_sbessel_integrator_integrate() with the Gaussian of
 * ncm_sbessel_integrator_integrate_gaussian_ell().
 */
void
ncm_sbessel_integrator_integrate_gaussian (NcmSBesselIntegrator *sbi, gdouble center, gdouble std, gdouble a, gdouble b, gdouble k, NcmVector *result)
{
  NcmSBesselIntegratorGaussianData data = {center, std};

  ncm_sbessel_integrator_integrate (sbi, &_ncm_sbessel_integrator_gaussian_func, a, b, k, result, &data);
}

typedef struct _NcmSBesselIntegratorRationalData
{
  gdouble center;
  gdouble std;
} NcmSBesselIntegratorRationalData;

static gdouble
_ncm_sbessel_integrator_rational_func (gpointer user_data, gdouble chi, gdouble k)
{
  NcmSBesselIntegratorRationalData *data = (NcmSBesselIntegratorRationalData *) user_data;
  const gdouble z                        = (chi - data->center) / data->std;
  const gdouble denom                    = 1.0 + z * z;
  const gdouble denom_cubed              = denom * denom * denom;

  return chi * chi / denom_cubed;
}

/**
 * ncm_sbessel_integrator_integrate_rational_ell:
 * @sbi: a #NcmSBesselIntegrator
 * @center: center $\chi_c$
 * @std: width $\sigma$
 * @a: lower limit
 * @b: upper limit
 * @k: wavenumber $k$
 * @ell: multipole $\ell \geq 0$
 *
 * ncm_sbessel_integrator_integrate_ell() with
 * $F = \chi^2 / [1 + (\chi - \chi_c)^2 / \sigma^2]^3$, the shape of the rational truth
 * tables, evaluated in C to avoid callback overhead.
 *
 * Returns: $I_\ell(k)$
 */
gdouble
ncm_sbessel_integrator_integrate_rational_ell (NcmSBesselIntegrator *sbi, gdouble center, gdouble std, gdouble a, gdouble b, gdouble k, gint ell)
{
  NcmSBesselIntegratorRationalData data = {center, std};

  return ncm_sbessel_integrator_integrate_ell (sbi, &_ncm_sbessel_integrator_rational_func, a, b, k, ell, &data);
}

/**
 * ncm_sbessel_integrator_integrate_rational:
 * @sbi: a #NcmSBesselIntegrator
 * @center: center $\chi_c$
 * @std: width $\sigma$
 * @a: lower limit
 * @b: upper limit
 * @k: wavenumber $k$
 * @result: output, one value per multipole
 *
 * ncm_sbessel_integrator_integrate() with the rational function of
 * ncm_sbessel_integrator_integrate_rational_ell().
 */
void
ncm_sbessel_integrator_integrate_rational (NcmSBesselIntegrator *sbi, gdouble center, gdouble std, gdouble a, gdouble b, gdouble k, NcmVector *result)
{
  NcmSBesselIntegratorRationalData data = {center, std};

  ncm_sbessel_integrator_integrate (sbi, &_ncm_sbessel_integrator_rational_func, a, b, k, result, &data);
}

