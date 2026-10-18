/***************************************************************************
 *            nc_dist_mu_nuisance.c
 *
 *  Sun Oct 18 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_dist_mu_nuisance.c
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
 * NcDistMuNuisance:
 *
 * Nuisance parameters of a distance-modulus likelihood.
 *
 * Two systematics of a catalog of distance moduli: the calibration offset
 * $\delta\mu_s$ of each sample $s$ of the catalog, added to the model modulus of every
 * object of that sample, and the peculiar velocity dispersion $\sigma_v$, which adds
 * $$
 * \sigma_{\mu,v}(z) = \frac{5}{\ln 10}\,\frac{\sigma_v}{c}\,\frac{1 + z}{z\,(1 + z/2)}
 * $$
 * in quadrature to the uncertainty of the modulus of an object at redshift $z$, the
 * same form as #NcSNIADistCov uses for its peculiar velocity term. The
 * number of samples is the construction property #NcDistMuNuisance:dmu-length; with no
 * samples the model carries only $\sigma_v$. The defaults, $\sigma_v = 0$ and
 * $\delta\mu_s = 0$, leave a likelihood unchanged.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_dist_mu_nuisance.h"
#include "ncm/core/ncm_c.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcDistMuNuisance
{
  /*< private >*/
  NcmModel parent_instance;
};

enum
{
  PROP_0,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDistMuNuisance, nc_dist_mu_nuisance, NCM_TYPE_MODEL)

#define VECTOR  (NCM_MODEL (nuisance))
#define SIGMA_V (ncm_model_orig_param_get (VECTOR, NC_DIST_MU_NUISANCE_SIGMA_V))

static void
nc_dist_mu_nuisance_init (NcDistMuNuisance *nuisance)
{
}

static void
_nc_dist_mu_nuisance_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_DIST_MU_NUISANCE (object));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_dist_mu_nuisance_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_DIST_MU_NUISANCE (object));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_dist_mu_nuisance_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_dist_mu_nuisance_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_dist_mu_nuisance, NC_TYPE_DIST_MU_NUISANCE);

static void
nc_dist_mu_nuisance_class_init (NcDistMuNuisanceClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_dist_mu_nuisance_set_property;
  model_class->get_property = &_nc_dist_mu_nuisance_get_property;
  object_class->finalize    = &_nc_dist_mu_nuisance_finalize;

  ncm_model_class_set_name_nick (model_class, "Distance modulus nuisance parameters", "DistMuNuisance");
  ncm_model_class_add_params (model_class, NC_DIST_MU_NUISANCE_SPARAM_LEN, NC_DIST_MU_NUISANCE_VPARAM_LEN, PROP_SIZE);

  ncm_mset_model_register_id (model_class,
                              "NcDistMuNuisance",
                              "Nuisance parameters of a distance-modulus likelihood.",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_set_sparam (model_class, NC_DIST_MU_NUISANCE_SIGMA_V, "\\sigma_v", "sigmav",
                              0.0, 2000.0, 10.0,
                              NC_DIST_MU_NUISANCE_DEFAULT_PARAMS_ABSTOL, NC_DIST_MU_NUISANCE_DEFAULT_SIGMA_V,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_vparam (model_class, NC_DIST_MU_NUISANCE_DMU, 0, "\\delta\\mu", "dmu",
                              -2.0, 2.0, 0.01,
                              NC_DIST_MU_NUISANCE_DEFAULT_PARAMS_ABSTOL, NC_DIST_MU_NUISANCE_DEFAULT_DMU,
                              NCM_PARAM_TYPE_FREE);

  ncm_model_class_check_params_info (model_class);
}

/**
 * nc_dist_mu_nuisance_new:
 * @nsamples: number of samples with an offset
 *
 * Creates a new #NcDistMuNuisance with @nsamples offsets and the default parameters.
 *
 * Returns: (transfer full): a new #NcDistMuNuisance
 */
NcDistMuNuisance *
nc_dist_mu_nuisance_new (guint nsamples)
{
  NcDistMuNuisance *nuisance = g_object_new (NC_TYPE_DIST_MU_NUISANCE,
                                             "dmu-length", nsamples,
                                             NULL);

  return nuisance;
}

/**
 * nc_dist_mu_nuisance_ref:
 * @nuisance: a #NcDistMuNuisance
 *
 * Increases the reference count of @nuisance by one.
 *
 * Returns: (transfer full): @nuisance
 */
NcDistMuNuisance *
nc_dist_mu_nuisance_ref (NcDistMuNuisance *nuisance)
{
  return g_object_ref (nuisance);
}

/**
 * nc_dist_mu_nuisance_free:
 * @nuisance: a #NcDistMuNuisance
 *
 * Decreases the reference count of @nuisance by one.
 *
 */
void
nc_dist_mu_nuisance_free (NcDistMuNuisance *nuisance)
{
  g_object_unref (nuisance);
}

/**
 * nc_dist_mu_nuisance_clear:
 * @nuisance: a #NcDistMuNuisance
 *
 * If *@nuisance is not %NULL, decreases its reference count and sets *@nuisance to %NULL.
 *
 */
void
nc_dist_mu_nuisance_clear (NcDistMuNuisance **nuisance)
{
  g_clear_object (nuisance);
}

/**
 * nc_dist_mu_nuisance_get_nsamples:
 * @nuisance: a #NcDistMuNuisance
 *
 * Returns: the number of samples with an offset
 */
guint
nc_dist_mu_nuisance_get_nsamples (NcDistMuNuisance *nuisance)
{
  return ncm_model_vparam_len (VECTOR, NC_DIST_MU_NUISANCE_DMU);
}

/**
 * nc_dist_mu_nuisance_get_sigma_v:
 * @nuisance: a #NcDistMuNuisance
 *
 * Returns: $\sigma_v$ [km/s]
 */
gdouble
nc_dist_mu_nuisance_get_sigma_v (NcDistMuNuisance *nuisance)
{
  return SIGMA_V;
}

/**
 * nc_dist_mu_nuisance_get_dmu:
 * @nuisance: a #NcDistMuNuisance
 * @sample: sample index, below nc_dist_mu_nuisance_get_nsamples()
 *
 * Returns: the offset $\delta\mu_s$ of @sample [mag]
 */
gdouble
nc_dist_mu_nuisance_get_dmu (NcDistMuNuisance *nuisance, const guint sample)
{
  return ncm_model_orig_vparam_get (VECTOR, NC_DIST_MU_NUISANCE_DMU, sample);
}

/**
 * nc_dist_mu_nuisance_sigma_mu_v:
 * @nuisance: a #NcDistMuNuisance
 * @z: redshift
 *
 * Returns: $\sigma_{\mu,v}(z)$ [mag], see #NcDistMuNuisance
 */
gdouble
nc_dist_mu_nuisance_sigma_mu_v (NcDistMuNuisance *nuisance, const gdouble z)
{
  const gdouble sigma_z = SIGMA_V / (1.0e-3 * ncm_c_c ());

  return 5.0 / M_LN10 * sigma_z * (1.0 + z) / (z * (1.0 + 0.5 * z));
}

