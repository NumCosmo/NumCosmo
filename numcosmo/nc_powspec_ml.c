/***************************************************************************
 *            nc_powspec_ml.c
 *
 *  Thu February 18 12:32:13 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_powspec_ml.c
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
 * SECTION:nc_powspec_ml
 * @title: NcPowspecML
 * @short_description: Abstrac class for linear matter power spectrum implementation.
 *
 * This module comprises the set of functions to compute the linear matter power spectrum and
 * derived quantities.
 * 
 * Following the description presented in #NcmPowspec, in this case we have that the field $\delta(\vec{x})$ 
 * represents the matter density fluctuations, i.e.,
 * $$\delta(\vec{x}) = \frac{\rho(\vec{x}) - \bar{\rho}}{\bar{\rho}},$$ 
 * where $\rho$ is the cold matter density field and $\bar{\rho}$ its mean.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_powspec_ml.h"
#include "math/ncm_c.h"
#include "math/integral.h"
#include "math/ncm_memory_pool.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_bessel.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_ABSTRACT_TYPE (NcPowspecML, nc_powspec_ml, NCM_TYPE_POWSPEC);

static void
nc_powspec_ml_init (NcPowspecML *nc_powspec_ml)
{
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
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = &_nc_powspec_ml_finalize;
}

/**
 * nc_powspec_ml_new_from_name:
 * @ps_ml_name: string which specifies the linear matter power spectrum object to be used
 *
 * This function returns a new #NcPowspecML whose type is defined by @ps_ml_name.
 *
 * Returns: A new #NcPowspecML.
 */
NcPowspecML *
nc_powspec_ml_new_from_name (const gchar *ps_ml_name)
{
  GObject *obj = ncm_serialize_global_from_string (ps_ml_name);

  if (!NC_IS_POWSPEC_ML (obj))
    g_error ("nc_powspec_ml_new_from_name: NcPowspecML %s do not descend from %s.", ps_ml_name, g_type_name (NC_TYPE_POWSPEC_ML));

  return NC_POWSPEC_ML (obj);
}

/**
 * nc_powspec_ml_ref:
 * @ps_ml: a #NcPowspecML
 *
 * Increases the reference count of @ps_ml atomically.
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
 * Decreases the reference count of @ps_ml atomically.
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
 * Decreses the reference count of *@ps_ml atomically and sets the pointer *@ps_ml to null.
 *
 */
void 
nc_powspec_ml_clear (NcPowspecML **ps_ml)
{
  g_clear_object (ps_ml);
}

typedef struct _NcPowspecMLSigmaInt
{
  const gdouble z;
  const gdouble R;
  NcPowspecML *ps_ml;
  NcmModel *model;
} NcPowspecMLSigmaInt;

static gdouble
_nc_powspec_ml_sigma_R_integ (gdouble k, gpointer user_data)
{
  NcPowspecMLSigmaInt *data = (NcPowspecMLSigmaInt *) user_data;
  const gdouble x  = k * data->R;
  const gdouble Pk = ncm_powspec_eval (NCM_POWSPEC (data->ps_ml), data->model, data->z, k);
  const gdouble W  = 3.0 * gsl_sf_bessel_j1 (x) / x;
  const gdouble W2 = W * W;
  
  return k * k * Pk * W2;
}

/**
 * nc_powspec_ml_sigma_R:
 * @ps_ml: a #NcPowspecML
 * @model: a #NcmModel
 * @reltol: relative tolerance for integration
 * @z: the value of $z$
 * @R: the value of $R$
 * 
 * Computes $\sigma_R = \sqrt{}$. FIXME
 * 
 */
gdouble 
nc_powspec_ml_sigma_R (NcPowspecML *ps_ml, NcmModel *model, const gdouble reltol, const gdouble z, const gdouble R)
{
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  NcPowspecMLSigmaInt data      = {z, R, ps_ml, model};
  const gdouble kmin            = ncm_powspec_get_kmin (NCM_POWSPEC (ps_ml));
  const gdouble kmax            = ncm_powspec_get_kmax (NCM_POWSPEC (ps_ml));
  const gdouble one_2pi2        = 1.0 / ncm_c_2_pi_2 ();
  gdouble error, sigma2_2pi2;
  gsl_function F;
  
  ncm_powspec_prepare_if_needed (NCM_POWSPEC (ps_ml), model);

  F.function = &_nc_powspec_ml_sigma_R_integ;
  F.params   = &data;
  
  gsl_integration_qag (&F, kmin, kmax, 0.0, reltol, NCM_INTEGRAL_PARTITION, 6, *w, &sigma2_2pi2, &error);

  ncm_memory_pool_return (w);  

  return sqrt (sigma2_2pi2 * one_2pi2);
}
