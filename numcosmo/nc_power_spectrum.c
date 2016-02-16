/***************************************************************************
 *            nc_power_spectrum.c
 *
 *  Tue February 16 17:00:52 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_power_spectrum.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_power_spectrum
 * @title: NcPowerSpectrum
 * @short_description: Abstrac class for power spectrum implementation.
 *
 * This module comprises the set of functions to compute a power spectrum and
 * derived quantities.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_power_spectrum.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_ZI,
  PROP_ZF,
  PROP_KMIN,
  PROP_KMAX
};

G_DEFINE_ABSTRACT_TYPE (NcPowerSpectrum, nc_power_spectrum, G_TYPE_OBJECT);

static void
nc_power_spectrum_init (NcPowerSpectrum *powspec)
{
  powspec->zi   = 0.0;
  powspec->zf   = 0.0;
  powspec->kmin = 0.0;
  powspec->kmax = 0.0;

  powspec->ctrl_cosmo = ncm_model_ctrl_new (NULL);
}

static void
nc_power_spectrum_dispose (GObject *object)
{
  NcPowerSpectrum *powspec = NC_POWER_SPECTRUM (object);

  ncm_model_ctrl_clear (&powspec->ctrl_cosmo);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_power_spectrum_parent_class)->dispose (object);
}

static void
nc_power_spectrum_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_power_spectrum_parent_class)->finalize (object);
}

static void
nc_power_spectrum_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcPowerSpectrum *powspec = NC_POWER_SPECTRUM (object);
  g_return_if_fail (NC_IS_POWER_SPECTRUM (object));

  switch (prop_id)
  {
    case PROP_ZI:
      nc_power_spectrum_set_zi (powspec, g_value_get_double (value));
      break;
    case PROP_ZF:
      nc_power_spectrum_set_zf (powspec, g_value_get_double (value));
      break;
    case PROP_KMIN:
      nc_power_spectrum_set_kmin (powspec, g_value_get_double (value));
      break;
    case PROP_KMAX:
      nc_power_spectrum_set_kmax (powspec, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_power_spectrum_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcPowerSpectrum *powspec = NC_POWER_SPECTRUM (object);
  g_return_if_fail (NC_IS_POWER_SPECTRUM (object));

  switch (prop_id)
  {
    case PROP_ZI:
      g_value_set_double (value, powspec->zi);
      break;
    case PROP_ZF:
      g_value_set_double (value, powspec->zf);
      break;
    case PROP_KMIN:
      g_value_set_double (value, powspec->kmin);
      break;
    case PROP_KMAX:
      g_value_set_double (value, powspec->kmax);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void _nc_power_spectrum_prepare (NcPowerSpectrum *powspec, NcHICosmo *cosmo) { g_error ("_nc_power_spectrum_prepare: no default implementation, all children must implement it."); } 
static gdouble _nc_power_spectrum_eval (NcPowerSpectrum *powspec, NcHICosmo *cosmo, const gdouble z, const gdouble k) { g_error ("_nc_power_spectrum_eval: no default implementation, all children must implement it."); return 0.0; }

static void
nc_power_spectrum_class_init (NcPowerSpectrumClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = nc_power_spectrum_set_property;
  object_class->get_property = nc_power_spectrum_get_property;

  object_class->dispose      = nc_power_spectrum_dispose;
  object_class->finalize     = nc_power_spectrum_finalize;

  g_object_class_install_property (object_class,
                                   PROP_ZI,
                                   g_param_spec_double ("zi",
                                                        NULL,
                                                        "Initial redshift",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_ZF, 
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "Final redshift",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_KMIN,
                                   g_param_spec_double ("kmin",
                                                        NULL,
                                                       "Minimum mode value",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_KMAX,
                                   g_param_spec_double ("kmax",
                                                        NULL,
                                                        "Maximum mode value",
                                                        0.0, G_MAXDOUBLE, 1.0e5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->prepare = &_nc_power_spectrum_prepare;
  klass->eval    = &_nc_power_spectrum_eval;
}

/**
 * nc_power_spectrum_ref:
 * @powspec: a #NcPowerSpectrum
 *
 * Increases the reference count of @powspec by one.
 *
 * Returns: (transfer full): @powspec.
 */
NcPowerSpectrum *
nc_power_spectrum_ref (NcPowerSpectrum *powspec)
{
  return g_object_ref (powspec);
}

/**
 * nc_power_spectrum_free:
 * @powspec: a #NcPowerSpectrum
 *
 * Decreases the reference count of @powspec by one.
 *
 */
void 
nc_power_spectrum_free (NcPowerSpectrum *powspec)
{
  g_object_unref (powspec);
}

/**
 * nc_power_spectrum_clear:
 * @powspec: a #NcPowerSpectrum
 *
 * If *@powspec is different from NULL, decreases the reference count of 
 * *@powspec by one and sets *@powspec to NULL.
 *
 */
void 
nc_power_spectrum_clear (NcPowerSpectrum **powspec)
{
  g_clear_object (powspec);
}

/**
 * nc_power_spectrum_set_zi:
 * @powspec: a #NcPowerSpectrum
 * @zi: initial redshift $z_i$
 * 
 * Sets the initial redshift $z_i$.
 *
 */
void 
nc_power_spectrum_set_zi (NcPowerSpectrum *powspec, const gdouble zi)
{
  powspec->zi = zi;
}

/**
 * nc_power_spectrum_set_zf:
 * @powspec: a #NcPowerSpectrum
 * @zf: minimum redshift $z_f$
 * 
 * Sets the final redshift $z_i$.
 *
 */
void 
nc_power_spectrum_set_zf (NcPowerSpectrum *powspec, const gdouble zf)
{
  powspec->zf = zf;
}

/**
 * nc_power_spectrum_set_kmin:
 * @powspec: a #NcPowerSpectrum
 * @kmin: minimum mode $k_\mathrm{min}$
 * 
 * Sets the minimum mode value $k_\mathrm{min}$.
 *
 */
void 
nc_power_spectrum_set_kmin (NcPowerSpectrum *powspec, const gdouble kmin)
{
  powspec->kmin = kmin;
}

/**
 * nc_power_spectrum_set_kmax:
 * @powspec: a #NcPowerSpectrum
 * @kmax: maxmimum mode $k_\mathrm{max}$
 * 
 * Sets the maximum mode value $k_\mathrm{max}$.
 *
 */
void 
nc_power_spectrum_set_kmax (NcPowerSpectrum *powspec, const gdouble kmax)
{
  powspec->kmax = kmax;
}

/**
 * nc_power_spectrum_prepare:
 * @powspec: a #NcPowerSpectrum
 * @cosmo: a #NcHICosmo
 * 
 * Prepares the power spectrum @powspec using the cosmology @cosmo.
 * 
 */
/**
 * nc_power_spectrum_prepare_if_needed:
 * @powspec: a #NcPowerSpectrum
 * @cosmo: a #NcHICosmo
 *
 * Prepare the object using the model @cosmo if it was changed
 * since last preparation.
 *
 */
/**
 * nc_power_spectrum_eval:
 * @powspec: a #NcPowerSpectrum
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 * @k: mode $k$
 * 
 * Evaluate the power spectrum @powspec at $(z, k)$.
 * 
 * Returns: $P(z, k)$.
 */
