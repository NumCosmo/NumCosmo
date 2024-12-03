/***************************************************************************
 *            nc_halo_density_profile_einasto.c
 *
 *  Wed July 17 12:33:27 2019
 *  Copyright  2014
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_halo_density_profile_einasto.c
 * Copyright (C) 2014 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_halo_density_profile_einasto
 * @title: NcHaloDensityProfileEinasto
 * @short_description: Density profile of Einasto type.
 *
 * This object implements the #NcHaloDensityProfile class for the Einasto density profile.
 *
 * As described #NcHaloDensityProfile, we just need to implement the dimensionless 3D density $\hat{\rho}(x)$
 * [which refers to the virtual function nc_halo_density_profile_eval_dl_density()].
 * In particular, the Einasto profile is given by
 * \begin{equation}
 * \hat{\rho}(x) = \exp \left[ - \frac{2}{\alpha} \left( x^\alpha - 1 \right) \right],
 * \end{equation}
 * where $x = r/r_s$ and $r_s$ is the scale radius and $\alpha$ is a parameter that defines how the
 * profile steepens with slope.
 *
 * Both the mass $M_\Delta$ and the scale profile $\rho_s$ are written in terms of the integral
 * $I_{x^2\hat\rho}(c_\Delta)$ [virtual function nc_halo_density_profile_eval_dl_spher_mass()].
 * The respective Einasto implementation provides
 * \begin{equation}
 * I_{x^2\hat\rho}(x) = \left(\frac{\alpha}{2}\right)^{3/\alpha} \frac{e^{2/\alpha}}{\alpha}
 * \Gamma \left( \frac{3}{\alpha} \right) \left[1 - \frac{\Gamma \left( 3/\alpha, 2c_{\Delta}^\alpha / \alpha\right)}{\Gamma \left( 3/\alpha \right)} \right],
 * \end{equation}
 * where $c_{\Delta}$ is the concentration parameter.
 *
 * References: Einasto (1965), arXiv:1712.04512. $\alpha$ parameter range values: Gao et al. (2008) https://ui.adsabs.harvard.edu/abs/2008MNRAS.387..536G/abstract,
 * Dutton \& Maccio (2014) https://academic.oup.com/mnras/article/441/4/3359/1209689.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_halo_density_profile_einasto.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_integrate.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_gamma.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcHaloDensityProfileEinastoPrivate
{
  NcHaloMassSummary *hms;
} NcHaloDensityProfileEinastoPrivate;

struct _NcHaloDensityProfileEinasto
{
  NcHaloDensityProfile parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloDensityProfileEinasto, nc_halo_density_profile_einasto, NC_TYPE_HALO_DENSITY_PROFILE);

#define VECTOR  (NCM_MODEL (dp))
/*#define M_DELTA (ncm_model_orig_param_get (VECTOR, NC_HALO_DENSITY_PROFILE_M_DELTA)) */
/*#define C_DELTA (ncm_model_orig_param_get (VECTOR, NC_HALO_DENSITY_PROFILE_C_DELTA)) */
#define ALPHA   (ncm_model_orig_param_get (VECTOR, NC_HALO_DENSITY_PROFILE_EINASTO_ALPHA))

enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_halo_density_profile_einasto_init (NcHaloDensityProfileEinasto *dpe)
{
  NcHaloDensityProfileEinastoPrivate * const self = nc_halo_density_profile_einasto_get_instance_private (dpe);

  self->hms = NULL;
}

static void
_nc_halo_density_profile_einasto_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /* NcHaloDensityProfileEinasto *dpe = NC_HALO_DENSITY_PROFILE_EINASTO (object); */
  g_return_if_fail (NC_IS_HALO_DENSITY_PROFILE_EINASTO (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_density_profile_einasto_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcHaloDensityProfileEinasto *dpe = NC_HALO_DENSITY_PROFILE_EINASTO (object); */
  g_return_if_fail (NC_IS_HALO_DENSITY_PROFILE_EINASTO (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_halo_density_profile_einasto_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_density_profile_einasto_parent_class)->finalize (object);
}

static void _nc_halo_density_profile_einasto_add_submodel (NcmModel *model, NcmModel *submodel);

static gdouble _nc_halo_density_profile_einasto_eval_dl_density (NcHaloDensityProfile *dp, const gdouble X);
static gdouble _nc_halo_density_profile_einasto_eval_dl_spher_mass (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble X);

static void
nc_halo_density_profile_einasto_class_init (NcHaloDensityProfileEinastoClass *klass)
{
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);
  NcHaloDensityProfileClass *dp_class = NC_HALO_DENSITY_PROFILE_CLASS (klass);
  NcmModelClass *model_class          = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_halo_density_profile_einasto_set_property;
  model_class->get_property = &_nc_halo_density_profile_einasto_get_property;
  object_class->finalize    = &nc_halo_density_profile_einasto_finalize;

  ncm_model_class_set_name_nick (model_class, "Einasto Density Profile", "Einasto");
  ncm_model_class_add_params (model_class, NC_HALO_DENSITY_PROFILE_EINASTO_LOCAL_SPARAM_LEN, 0, PROP_SIZE);

  /**
   * NcHaloDensityProfileEinasto:alpha:
   *
   * Defines how the profile steepens with slope.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_DENSITY_PROFILE_EINASTO_ALPHA, "\\alpha", "alpha",
                              0.12,  0.35, 0.01,
                              NC_HALO_DENSITY_PROFILE_EINASTO_DEFAULT_PARAMS_ABSTOL, NC_HALO_DENSITY_PROFILE_EINASTO_DEFAULT_ALPHA,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  model_class->add_submodel = &_nc_halo_density_profile_einasto_add_submodel;


  dp_class->eval_dl_density    = &_nc_halo_density_profile_einasto_eval_dl_density;
  dp_class->eval_dl_spher_mass = &_nc_halo_density_profile_einasto_eval_dl_spher_mass;
}

static void
_nc_halo_density_profile_einasto_add_submodel (NcmModel *model, NcmModel *submodel)
{
  /* Chain up : start */
  NCM_MODEL_CLASS (nc_halo_density_profile_einasto_parent_class)->add_submodel (model, submodel);
  {
    NcHaloDensityProfileEinasto *dpe                = NC_HALO_DENSITY_PROFILE_EINASTO (model);
    NcHaloDensityProfileEinastoPrivate * const self = nc_halo_density_profile_einasto_get_instance_private (dpe);

    if (ncm_model_id (submodel) == nc_halo_mass_summary_id ())
    {
      if (self->hms != NULL)
        g_error ("Halo mass summary was already defined in `%s' class.", G_OBJECT_CLASS_NAME (dpe));

      self->hms = NC_HALO_MASS_SUMMARY (submodel);
    }
  }
}

static gdouble
_nc_halo_density_profile_einasto_eval_dl_density (NcHaloDensityProfile *dp, const gdouble x)
{
  return exp (-2.0 / ALPHA * (pow (x, ALPHA) - 1.0));
}

static gdouble
_nc_halo_density_profile_einasto_eval_dl_spher_mass (NcHaloDensityProfile *dp, NcHICosmo *cosmo, const gdouble x)
{
  NcHaloDensityProfileEinasto *dpe                = NC_HALO_DENSITY_PROFILE_EINASTO (dp);
  NcHaloDensityProfileEinastoPrivate * const self = nc_halo_density_profile_einasto_get_instance_private (dpe);
  const gdouble cDelta                            = nc_halo_mass_summary_concentration (self->hms, cosmo);
  const gdouble gamma_3_alpha                     = gsl_sf_gamma (3.0 / ALPHA);
  const gdouble arg_2                             = 2.0 * pow (cDelta, ALPHA) / ALPHA;
  const gdouble gamma_inc_P                       = gsl_sf_gamma_inc_P (3.0 / ALPHA, arg_2);

  return (pow (ALPHA / 2.0, 3.0 / ALPHA) * exp (2.0 / ALPHA) / ALPHA * gamma_3_alpha * gamma_inc_P);
}

/**
 * nc_halo_density_profile_einasto_new:
 * @hms: a #NcHaloMassSummary
 *
 * This function returns the #NcHaloDensityProfileEinasto implementation of
 * #NcHaloDensityProfile given #NcHaloMassSummary.
 *
 * Returns: a new instance of #NcHaloDensityProfileEinasto.
 */
NcHaloDensityProfileEinasto *
nc_halo_density_profile_einasto_new (NcHaloMassSummary *hms)
{
  NcmObjArray *submodels = ncm_obj_array_new ();

  ncm_obj_array_add (submodels, G_OBJECT (hms));
  {
    NcHaloDensityProfileEinasto *dp_einasto = g_object_new (NC_TYPE_HALO_DENSITY_PROFILE_EINASTO,
                                                            "submodel-array", submodels,
                                                            NULL);

    ncm_obj_array_unref (submodels);

    return dp_einasto;
  }
}

