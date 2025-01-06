/***************************************************************************
 *            nc_powspec_ml.c
 *
 *  Thu February 18 12:32:13 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_powspec_ml.c
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
 * NcPowspecML:
 *
 * Abstract class for linear matter power spectrum implementation.
 *
 * This module comprises the set of functions to compute the linear matter power
 * spectrum and derived quantities.
 *
 * Following the description presented in #NcmPowspec, in this case we have that the
 * field $\delta(\vec{x})$ represents the matter density fluctuations, i.e.,
 * $$\delta(\vec{x}) = \frac{\rho(\vec{x}) - \bar{\rho}}{\bar{\rho}},$$ where $\rho$ is
 * the cold matter density field and $\bar{\rho}$ its mean.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_powspec_ml.h"

enum
{
  PROP_0,
  PROP_ZI,
  PROP_ZF,
  PROP_KMIN,
  PROP_KMAX,
  PROP_SIZE
};

G_DEFINE_ABSTRACT_TYPE (NcPowspecML, nc_powspec_ml, NCM_TYPE_POWSPEC)

static void
nc_powspec_ml_init (NcPowspecML *nc_powspec_ml)
{
}

static void
_nc_powspec_ml_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcPowspecML *ps_mlt = NC_POWSPEC_ML (object);
  NcmPowspec *ps      = NCM_POWSPEC (ps_mlt);

  g_return_if_fail (NC_IS_POWSPEC_ML (object));

  switch (prop_id)
  {
    case PROP_ZI:
      ncm_powspec_set_zi (ps, g_value_get_double (value));
      break;
    case PROP_ZF:
      ncm_powspec_set_zf (ps, g_value_get_double (value));
      break;
    case PROP_KMIN:
      ncm_powspec_set_kmin (ps, g_value_get_double (value));
      break;
    case PROP_KMAX:
      ncm_powspec_set_kmax (ps, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_powspec_ml_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcPowspecML *ps_mlt = NC_POWSPEC_ML (object);
  NcmPowspec *ps      = NCM_POWSPEC (ps_mlt);

  g_return_if_fail (NC_IS_POWSPEC_ML (object));

  switch (prop_id)
  {
    case PROP_ZI:
      g_value_set_double (value, ncm_powspec_get_zi (ps));
      break;
    case PROP_ZF:
      g_value_set_double (value, ncm_powspec_get_zf (ps));
      break;
    case PROP_KMIN:
      g_value_set_double (value, ncm_powspec_get_kmin (ps));
      break;
    case PROP_KMAX:
      g_value_set_double (value, ncm_powspec_get_kmax (ps));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_powspec_ml_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_powspec_ml_parent_class)->finalize (object);
}

static void
nc_powspec_ml_class_init (NcPowspecMLClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_powspec_ml_set_property;
  object_class->get_property = &_nc_powspec_ml_get_property;
  object_class->finalize     = &_nc_powspec_ml_finalize;

  /**
   * NcPowspecML:zi:
   *
   * The initial time (redshift) to compute $P(k,z)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_ZI,
                                   g_param_spec_double ("zi",
                                                        NULL,
                                                        "Initial redshift",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcPowspecML:zf:
   *
   * The final time (redshift) to compute $P(k,z)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_ZF,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "Final redshift",
                                                        0.0, G_MAXDOUBLE, 5.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcPowspecML:kmin:
   *
   * The minimum mode (wave-number) value to compute $P(k,z)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_KMIN,
                                   g_param_spec_double ("kmin",
                                                        NULL,
                                                        "Minimum mode value",
                                                        0.0, G_MAXDOUBLE, 1.0e-6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcPowspecML:kmax:
   *
   * The maximum mode (wave-number) value to compute $P(k,z)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_KMAX,
                                   g_param_spec_double ("kmax",
                                                        NULL,
                                                        "Maximum mode value",
                                                        0.0, G_MAXDOUBLE, 1.0e3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_powspec_ml_ref:
 * @ps_ml: a #NcPowspecML
 *
 * Increases the reference count of @ps_ml by one atomically.
 *
 * Returns: (transfer full): @ps_ml.
 */
NcPowspecML *
nc_powspec_ml_ref (NcPowspecML *ps_ml)
{
  return g_object_ref (ps_ml);
}

/**
 * nc_powspec_ml_free:
 * @ps_ml: a #NcPowspecML
 *
 * Atomically decrements the reference count of @ps_ml by one.
 * If the reference count drops to 0, all memory allocated by @ps_ml is released.
 *
 */
void
nc_powspec_ml_free (NcPowspecML *ps_ml)
{
  g_object_unref (ps_ml);
}

/**
 * nc_powspec_ml_clear:
 * @ps_ml: a #NcPowspecML
 *
 *  If @ps_ml is different from NULL,
 *  atomically decrements the reference count of @ps_ml by one.
 *  If the reference count drops to 0,
 *  all memory allocated by @ps_ml is released and @ps_ml is set to NULL.
 *
 */
void
nc_powspec_ml_clear (NcPowspecML **ps_ml)
{
  g_clear_object (ps_ml);
}

