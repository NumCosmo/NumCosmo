/***************************************************************************
 *            nc_multiplicity_func_bhattacharya.c
 *
 *  Mon May 18 14:26:00 2026
 *  Copyright  2026  Cinthia Nunes de Lima / Henrique C. N. Lettieri
 *  <cinthia.nlima@gmail.com> <henrique.cnl@hotmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026  Cinthia Nunes de Lima / Henrique C. N. Lettieri
 * <cinthia.nlima@gmail.com> <henrique.cnl@hotmail.com>
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
 * NcMultiplicityFuncBhattacharya:
 *
 * Dark matter halo -- Bhattacharya multiplicity function.
 *
 * Dark matter halo multiplicity function fitted to LCDM simulations using the
 * friends of friends algorithm, FoF(0.2). See reference arXiv:1005.2239.
 *
 * $$f(\sigma, z) = A(z) \sqrt{\frac{2}{\pi}} e^{-\frac{a(z) \delta_c^2}{2 \sigma^2}}
 * \left[ 1 + \left( \frac{\sigma^2}{a(z) \delta_c^2} \right)^{p} \right]
 * \left( \frac{\delta_c \sqrt{a(z)}}{\sigma} \right)^{q},$$
 * where $A(z) = A_0 (1+z)^{-0.11}$, $a(z) = a_0 (1+z)^{-0.01}$, and the default values
 * fitted in the range $0 \leq z \leq 2$ are $A_0 = 0.333$, $a_0 = 0.788$, $p = 0.807$
 * and $q = 1.795$.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_bhattacharya.h"

struct _NcMultiplicityFuncBhattacharyaPrivate
{
  NcMultiplicityFuncMassDef mdef;
  gdouble A;
  gdouble a;
  gdouble p;
  gdouble q;
  gdouble delta_c;
  gdouble Delta;
};

enum
{
  PROP_0,
  PROP_A,
  PROP_a,
  PROP_p,
  PROP_q,
  PROP_DELTA_C,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncBhattacharya, nc_multiplicity_func_bhattacharya, NC_TYPE_MULTIPLICITY_FUNC)

static void
nc_multiplicity_func_bhattacharya_init (NcMultiplicityFuncBhattacharya *mbt)
{
  NcMultiplicityFuncBhattacharyaPrivate * const self = mbt->priv = nc_multiplicity_func_bhattacharya_get_instance_private (mbt);

  self->mdef    = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;
  self->A       = 0.0;
  self->a       = 0.0;
  self->p       = 0.0;
  self->q       = 0.0;
  self->delta_c = 0.0;
  self->Delta   = 0.0;
}

static void
_nc_multiplicity_func_bhattacharya_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncBhattacharya *mbt = NC_MULTIPLICITY_FUNC_BHATTACHARYA (object);

  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_BHATTACHARYA (object));

  switch (prop_id)
  {
    case PROP_A:
      nc_multiplicity_func_bhattacharya_set_A (mbt, g_value_get_double (value));
      break;
    case PROP_a:
      nc_multiplicity_func_bhattacharya_set_a (mbt, g_value_get_double (value));
      break;
    case PROP_p:
      nc_multiplicity_func_bhattacharya_set_p (mbt, g_value_get_double (value));
      break;
    case PROP_q:
      nc_multiplicity_func_bhattacharya_set_q (mbt, g_value_get_double (value));
      break;
    case PROP_DELTA_C:
      nc_multiplicity_func_bhattacharya_set_delta_c (mbt, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_multiplicity_func_bhattacharya_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncBhattacharya *mbt = NC_MULTIPLICITY_FUNC_BHATTACHARYA (object);

  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_BHATTACHARYA (object));

  switch (prop_id)
  {
    case PROP_A:
      g_value_set_double (value, nc_multiplicity_func_bhattacharya_get_A (mbt));
      break;
    case PROP_a:
      g_value_set_double (value, nc_multiplicity_func_bhattacharya_get_a (mbt));
      break;
    case PROP_p:
      g_value_set_double (value, nc_multiplicity_func_bhattacharya_get_p (mbt));
      break;
    case PROP_q:
      g_value_set_double (value, nc_multiplicity_func_bhattacharya_get_q (mbt));
      break;
    case PROP_DELTA_C:
      g_value_set_double (value, nc_multiplicity_func_bhattacharya_get_delta_c (mbt));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_multiplicity_func_bhattacharya_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_bhattacharya_parent_class)->finalize (object);
}

static void _nc_multiplicity_func_bhattacharya_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef);
static void _nc_multiplicity_func_bhattacharya_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta);
static NcMultiplicityFuncMassDef _nc_multiplicity_func_bhattacharya_get_mdef (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_bhattacharya_get_Delta (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_bhattacharya_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

static void
nc_multiplicity_func_bhattacharya_class_init (NcMultiplicityFuncBhattacharyaClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass *parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = _nc_multiplicity_func_bhattacharya_set_property;
  object_class->get_property = _nc_multiplicity_func_bhattacharya_get_property;
  object_class->finalize     = _nc_multiplicity_func_bhattacharya_finalize;

  /**
   * NcMultiplicityFuncBhattacharya:A:
   *
   * Normalization parameter in the Bhattacharya multiplicity function.
   */
  g_object_class_install_property (object_class,
                                   PROP_A,
                                   g_param_spec_double ("A",
                                                        NULL,
                                                        "A",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.333,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcMultiplicityFuncBhattacharya:a:
   *
   * Shape parameter in the Bhattacharya multiplicity function.
   */
  g_object_class_install_property (object_class,
                                   PROP_a,
                                   g_param_spec_double ("a",
                                                        NULL,
                                                        "a",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.788,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcMultiplicityFuncBhattacharya:p:
   *
   * Shape parameter in the Bhattacharya multiplicity function.
   */
  g_object_class_install_property (object_class,
                                   PROP_p,
                                   g_param_spec_double ("p",
                                                        NULL,
                                                        "p",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.807,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcMultiplicityFuncBhattacharya:q:
   *
   * Power-law index parameter in the Bhattacharya multiplicity function.
   */
  g_object_class_install_property (object_class,
                                   PROP_q,
                                   g_param_spec_double ("q",
                                                        NULL,
                                                        "q",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.795,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcMultiplicityFuncBhattacharya:critical-delta:
   *
   * Critical overdensity for spherical collapse.
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA_C,
                                   g_param_spec_double ("critical-delta",
                                                        NULL,
                                                        "Critical delta",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, NC_MULTIPLICITY_FUNC_DELTA_C0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  parent_class->set_mdef  = &_nc_multiplicity_func_bhattacharya_set_mdef;
  parent_class->get_mdef  = &_nc_multiplicity_func_bhattacharya_get_mdef;
  parent_class->set_Delta = &_nc_multiplicity_func_bhattacharya_set_Delta;
  parent_class->get_Delta = &_nc_multiplicity_func_bhattacharya_get_Delta;
  parent_class->eval      = &_nc_multiplicity_func_bhattacharya_eval;
}

static void
_nc_multiplicity_func_bhattacharya_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncBhattacharya *mbt                = NC_MULTIPLICITY_FUNC_BHATTACHARYA (mulf);
  NcMultiplicityFuncBhattacharyaPrivate * const self = mbt->priv;

  switch (mdef)
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_FOF:
      /* nothing to do */
      break;
    /* LCOV_EXCL_START */
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
      g_error ("NcMultiplicityFuncBhattacharya does not support mean mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
      g_error ("NcMultiplicityFuncBhattacharya does not support critical mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
      g_error ("NcMultiplicityFuncBhattacharya does not support virial mass def");
      break;
    default:
      g_assert_not_reached ();
      break;
      /* LCOV_EXCL_STOP */
  }

  self->mdef = mdef;
}

static NcMultiplicityFuncMassDef
_nc_multiplicity_func_bhattacharya_get_mdef (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncBhattacharya *mba                = NC_MULTIPLICITY_FUNC_BHATTACHARYA (mulf);
  NcMultiplicityFuncBhattacharyaPrivate * const self = mba->priv;

  return self->mdef;
}

static void
_nc_multiplicity_func_bhattacharya_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta)
{
  NcMultiplicityFuncBhattacharya *mba                = NC_MULTIPLICITY_FUNC_BHATTACHARYA (mulf);
  NcMultiplicityFuncBhattacharyaPrivate * const self = mba->priv;

  self->Delta = Delta;
}

static gdouble
_nc_multiplicity_func_bhattacharya_get_Delta (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncBhattacharya *mba                = NC_MULTIPLICITY_FUNC_BHATTACHARYA (mulf);
  NcMultiplicityFuncBhattacharyaPrivate * const self = mba->priv;

  return self->Delta;
}

static gdouble
_nc_multiplicity_func_bhattacharya_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z) /* f(\sigma) - Bhattacharya */
{
  NcMultiplicityFuncBhattacharya *mba                = NC_MULTIPLICITY_FUNC_BHATTACHARYA (mulf);
  NcMultiplicityFuncBhattacharyaPrivate * const self = mba->priv;

  const gdouble A   = self->A / pow (1.0 + z, 0.11);
  const gdouble a   = self->a / pow (1.0 + z, 0.01);
  const gdouble bc1 = sqrt (2.0 / M_PI);
  const gdouble p   = self->p;
  const gdouble q   = self->q;
  const gdouble x   = self->delta_c / sigma;
  const gdouble x2  = x * x;

  const gdouble f_Bhattacharya = A * bc1 * exp (-0.5 * a * x2) * (1.0 + pow (a * x2, -p)) * pow (x * sqrt (a), q);

  NCM_UNUSED (cosmo);

  return f_Bhattacharya;
}

/**
 * nc_multiplicity_func_bhattacharya_new:
 *
 * Creates a new #NcMultiplicityFuncBhattacharya with default parameters.
 *
 * Returns: A new #NcMultiplicityFuncBhattacharya.
 */
NcMultiplicityFuncBhattacharya *
nc_multiplicity_func_bhattacharya_new (void)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_BHATTACHARYA,
                       "mass-def", NC_MULTIPLICITY_FUNC_MASS_DEF_FOF,
                       NULL);
}

/**
 * nc_multiplicity_func_bhattacharya_ref:
 * @mbt: a #NcMultiplicityFuncBhattacharya
 *
 * Increases the reference count of @mbt by one.
 *
 * Returns: (transfer full): @mbt
 */
NcMultiplicityFuncBhattacharya *
nc_multiplicity_func_bhattacharya_ref (NcMultiplicityFuncBhattacharya *mbt)
{
  return g_object_ref (mbt);
}

/**
 * nc_multiplicity_func_bhattacharya_free:
 * @mbt: a #NcMultiplicityFuncBhattacharya
 *
 * Atomically decrements the reference count of @mbt by one. If the reference count drops to 0,
 * all memory allocated by @mbt is released.
 *
 */
void
nc_multiplicity_func_bhattacharya_free (NcMultiplicityFuncBhattacharya *mbt)
{
  g_object_unref (mbt);
}

/**
 * nc_multiplicity_func_bhattacharya_clear:
 * @mbt: a #NcMultiplicityFuncBhattacharya
 *
 * Atomically decrements the reference count of @mbt by one. If the reference count drops to 0,
 * all memory allocated by @mbt is released. Set the pointer to NULL;
 *
 */
void
nc_multiplicity_func_bhattacharya_clear (NcMultiplicityFuncBhattacharya **mbt)
{
  g_clear_object (mbt);
}

/**
 * nc_multiplicity_func_bhattacharya_set_A:
 * @mbt: a #NcMultiplicityFuncBhattacharya.
 * @A: value of #NcMultiplicityFuncBhattacharya:A.
 *
 * Sets the value @A to the #NcMultiplicityFuncBhattacharya:A property.
 *
 */
void
nc_multiplicity_func_bhattacharya_set_A (NcMultiplicityFuncBhattacharya *mbt, gdouble A)
{
  NcMultiplicityFuncBhattacharyaPrivate * const self = mbt->priv;

  g_assert (A >= 0);

  self->A = A;
}

/**
 * nc_multiplicity_func_bhattacharya_get_A:
 * @mbt: a #NcMultiplicityFuncBhattacharya.
 *
 * Returns: the value of #NcMultiplicityFuncBhattacharya:A property.
 */
gdouble
nc_multiplicity_func_bhattacharya_get_A (const NcMultiplicityFuncBhattacharya *mbt)
{
  NcMultiplicityFuncBhattacharyaPrivate * const self = mbt->priv;

  return self->A;
}

/**
 * nc_multiplicity_func_bhattacharya_set_a:
 * @mbt: a #NcMultiplicityFuncBhattacharya.
 * @a: value of #NcMultiplicityFuncBhattacharya:a.
 *
 * Sets the value @a to the #NcMultiplicityFuncBhattacharya:a property.
 *
 */
void
nc_multiplicity_func_bhattacharya_set_a (NcMultiplicityFuncBhattacharya *mbt, gdouble a)
{
  NcMultiplicityFuncBhattacharyaPrivate * const self = mbt->priv;

  g_assert (a >= 0);

  self->a = a;
}

/**
 * nc_multiplicity_func_bhattacharya_get_a:
 * @mbt: a #NcMultiplicityFuncBhattacharya.
 *
 * Returns: the value of #NcMultiplicityFuncBhattacharya:a property.
 */
gdouble
nc_multiplicity_func_bhattacharya_get_a (const NcMultiplicityFuncBhattacharya *mbt)
{
  NcMultiplicityFuncBhattacharyaPrivate * const self = mbt->priv;

  return self->a;
}

/**
 * nc_multiplicity_func_bhattacharya_set_p:
 * @mbt: a #NcMultiplicityFuncBhattacharya.
 * @p: value of #NcMultiplicityFuncBhattacharya:p.
 *
 * Sets the value @p to the #NcMultiplicityFuncBhattacharya:p property.
 *
 */
void
nc_multiplicity_func_bhattacharya_set_p (NcMultiplicityFuncBhattacharya *mbt, gdouble p)
{
  NcMultiplicityFuncBhattacharyaPrivate * const self = mbt->priv;

  g_assert (p >= 0);

  self->p = p;
}

/**
 * nc_multiplicity_func_bhattacharya_get_p:
 * @mbt: a #NcMultiplicityFuncBhattacharya.
 *
 * Returns: the value of #NcMultiplicityFuncBhattacharya:p property.
 */
gdouble
nc_multiplicity_func_bhattacharya_get_p (const NcMultiplicityFuncBhattacharya *mbt)
{
  NcMultiplicityFuncBhattacharyaPrivate * const self = mbt->priv;

  return self->p;
}

/**
 * nc_multiplicity_func_bhattacharya_set_q:
 * @mbt: a #NcMultiplicityFuncBhattacharya.
 * @q: value of #NcMultiplicityFuncBhattacharya:q.
 *
 * Sets the value @q to the #NcMultiplicityFuncBhattacharya:q property.
 *
 */
void
nc_multiplicity_func_bhattacharya_set_q (NcMultiplicityFuncBhattacharya *mbt, gdouble q)
{
  NcMultiplicityFuncBhattacharyaPrivate * const self = mbt->priv;

  g_assert (q >= 0);

  self->q = q;
}

/**
 * nc_multiplicity_func_bhattacharya_get_q:
 * @mbt: a #NcMultiplicityFuncBhattacharya.
 *
 * Returns: the value of #NcMultiplicityFuncBhattacharya:q property.
 */
gdouble
nc_multiplicity_func_bhattacharya_get_q (const NcMultiplicityFuncBhattacharya *mbt)
{
  NcMultiplicityFuncBhattacharyaPrivate * const self = mbt->priv;

  return self->q;
}

/**
 * nc_multiplicity_func_bhattacharya_set_delta_c:
 * @mbt: a #NcMultiplicityFuncBhattacharya.
 * @delta_c: value of #NcMultiplicityFuncBhattacharya:critical-delta.
 *
 * Sets the value @delta_c to the #NcMultiplicityFuncBhattacharya:critical-delta property.
 *
 */
void
nc_multiplicity_func_bhattacharya_set_delta_c (NcMultiplicityFuncBhattacharya *mbt, gdouble delta_c)
{
  NcMultiplicityFuncBhattacharyaPrivate * const self = mbt->priv;

  g_assert (delta_c >= 0);

  self->delta_c = delta_c;
}

/**
 * nc_multiplicity_func_bhattacharya_get_delta_c:
 * @mbt: a #NcMultiplicityFuncBhattacharya.
 *
 * Returns: the value of #NcMultiplicityFuncBhattacharya:critical-delta property.
 */
gdouble
nc_multiplicity_func_bhattacharya_get_delta_c (const NcMultiplicityFuncBhattacharya *mbt)
{
  NcMultiplicityFuncBhattacharyaPrivate * const self = mbt->priv;

  return self->delta_c;
}

