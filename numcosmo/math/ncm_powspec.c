/***************************************************************************
 *            ncm_powspec.c
 *
 *  Tue February 16 17:00:52 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_powspec.c
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
 * SECTION:ncm_powspec
 * @title: NcmPowspec
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

#include "math/ncm_powspec.h"
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

G_DEFINE_ABSTRACT_TYPE (NcmPowspec, ncm_powspec, G_TYPE_OBJECT);

static void
ncm_powspec_init (NcmPowspec *powspec)
{
  powspec->zi   = 0.0;
  powspec->zf   = 0.0;
  powspec->kmin = 0.0;
  powspec->kmax = 0.0;

  powspec->ctrl = ncm_model_ctrl_new (NULL);
}

static void
ncm_powspec_dispose (GObject *object)
{
  NcmPowspec *powspec = NCM_POWSPEC (object);

  ncm_model_ctrl_clear (&powspec->ctrl);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_parent_class)->dispose (object);
}

static void
ncm_powspec_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_parent_class)->finalize (object);
}

static void
ncm_powspec_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_powspec_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPowspec *powspec = NCM_POWSPEC (object);
  g_return_if_fail (NCM_IS_POWSPEC (object));

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

static void _ncm_powspec_prepare (NcmPowspec *powspec, NcmModel *model) { g_error ("_ncm_powspec_prepare: no default implementation, all children must implement it."); } 
static gdouble _ncm_powspec_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k) { g_error ("_ncm_powspec_eval: no default implementation, all children must implement it."); return 0.0; }

static void
ncm_powspec_class_init (NcmPowspecClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = ncm_powspec_set_property;
  object_class->get_property = ncm_powspec_get_property;

  object_class->dispose      = ncm_powspec_dispose;
  object_class->finalize     = ncm_powspec_finalize;

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

  klass->prepare = &_ncm_powspec_prepare;
  klass->eval    = &_ncm_powspec_eval;
}

/**
 * ncm_powspec_ref:
 * @powspec: a #NcmPowspec
 *
 * Increases the reference count of @powspec by one.
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
 * Decreases the reference count of @powspec by one.
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
 * If *@powspec is different from NULL, decreases the reference count of 
 * *@powspec by one and sets *@powspec to NULL.
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
 * @zi: initial redshift $z_i$
 * 
 * Sets the initial redshift $z_i$.
 *
 */
void 
ncm_powspec_set_zi (NcmPowspec *powspec, const gdouble zi)
{
  powspec->zi = zi;
}

/**
 * ncm_powspec_set_zf:
 * @powspec: a #NcmPowspec
 * @zf: minimum redshift $z_f$
 * 
 * Sets the final redshift $z_i$.
 *
 */
void 
ncm_powspec_set_zf (NcmPowspec *powspec, const gdouble zf)
{
  powspec->zf = zf;
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
  powspec->kmin = kmin;
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
  powspec->kmax = kmax;
}

/**
 * ncm_powspec_prepare:
 * @powspec: a #NcmPowspec
 * @cosmo: a #NcHICosmo
 * 
 * Prepares the power spectrum @powspec using the cosmology @cosmo.
 * 
 */
/**
 * ncm_powspec_prepare_if_needed:
 * @powspec: a #NcmPowspec
 * @cosmo: a #NcHICosmo
 *
 * Prepare the object using the model @cosmo if it was changed
 * since last preparation.
 *
 */
/**
 * ncm_powspec_eval:
 * @powspec: a #NcmPowspec
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 * @k: mode $k$
 * 
 * Evaluate the power spectrum @powspec at $(z, k)$.
 * 
 * Returns: $P(z, k)$.
 */
