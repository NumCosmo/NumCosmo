/***************************************************************************
 *            nc_galaxy_shape_pop_shear_scatter_gauss.c
 *
 *  Thu Nov 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop_shear_scatter_gauss.c
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyShapePopShearScatterGauss:
 *
 * Gaussian scatter of the tangential shear estimate.
 *
 * The residual $\delta = \hat g_t - \mu_t$ is Gaussian with variance
 * $\sigma_g^2 + \sigma_{\mathrm{noise},g}^2$, where
 * $\sigma_{\mathrm{noise},g}^2$ is the measurement-noise variance in shear
 * units (dropped when #NcGalaxyShapePopShearScatter:include-noise is
 * %FALSE). The width $\sigma_g$ is either the model parameter
 * #NcGalaxyShapePopShearScatterGauss:sigma, in reduced-shear units and
 * common to every galaxy, or, when
 * #NcGalaxyShapePopShearScatterGauss:sigma-from-catalog is %TRUE, the
 * per-galaxy catalog column "e_rms" (the per-component intrinsic
 * ellipticity dispersion, as in the HSC shape catalogs) converted to shear
 * units, $\sigma_g = e_\mathrm{rms}/(sR)$. In the latter mode the parameter
 * sigma is not used and should be kept fixed.
 *
 * With sigma-from-catalog, #NcGalaxyShapeFactorWtG reproduces the
 * per-galaxy Gaussian likelihood commonly used as the Gaussian variant of
 * the Applegate et al. (2014) method, with variance
 * $(e_\mathrm{rms}^2 + \sigma_\mathrm{noise}^2)/(sR)^2$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_pop_shear_scatter_gauss.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyShapePopShearScatterGauss
{
  NcGalaxyShapePopShearScatter parent_instance;
  gboolean sigma_from_catalog;
};

enum
{
  PROP_0,
  PROP_SIGMA_FROM_CATALOG,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxyShapePopShearScatterGauss, nc_galaxy_shape_pop_shear_scatter_gauss, NC_TYPE_GALAXY_SHAPE_POP_SHEAR_SCATTER)

static void
nc_galaxy_shape_pop_shear_scatter_gauss_init (NcGalaxyShapePopShearScatterGauss *gspssg)
{
  gspssg->sigma_from_catalog = FALSE;
}

static void
_nc_galaxy_shape_pop_shear_scatter_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapePopShearScatterGauss *gspssg = NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS (gspssg));

  switch (prop_id)
  {
    case PROP_SIGMA_FROM_CATALOG:
      gspssg->sigma_from_catalog = g_value_get_boolean (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_pop_shear_scatter_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapePopShearScatterGauss *gspssg = NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS (object);

  g_return_if_fail (NC_IS_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS (gspssg));

  switch (prop_id)
  {
    case PROP_SIGMA_FROM_CATALOG:
      g_value_set_boolean (value, gspssg->sigma_from_catalog);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_pop_shear_scatter_gauss_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_shape_pop_shear_scatter_gauss_parent_class)->finalize (object);
}

#define VECTOR (NCM_MODEL (gspss))
#define SIGMA  (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_SIGMA))

static void
_nc_galaxy_shape_pop_shear_scatter_gauss_ldata_noop (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_nc_galaxy_shape_pop_shear_scatter_gauss_ldata_no_columns (NcGalaxyShapePopData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_pop_shear_scatter_gauss_ldata_read_e_rms (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i)
{
  data->e_rms = ncm_catalog_get (NCM_CATALOG (obs), NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_COL_E_RMS, i, NULL);
}

static void
_nc_galaxy_shape_pop_shear_scatter_gauss_ldata_write_e_rms (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i)
{
  ncm_catalog_set (NCM_CATALOG (obs), NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_COL_E_RMS, i, data->e_rms, NULL);
}

static void
_nc_galaxy_shape_pop_shear_scatter_gauss_ldata_e_rms_column (NcGalaxyShapePopData *data, GList **columns)
{
  *columns = g_list_append (*columns, g_strdup (NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_COL_E_RMS));
}

static void
_nc_galaxy_shape_pop_shear_scatter_gauss_data_init (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  NcGalaxyShapePopShearScatterGauss *gspssg = NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS (gsp);

  data->ldata         = NULL;
  data->ldata_destroy = &g_free;

  if (gspssg->sigma_from_catalog)
  {
    data->ldata_read_row         = &_nc_galaxy_shape_pop_shear_scatter_gauss_ldata_read_e_rms;
    data->ldata_write_row        = &_nc_galaxy_shape_pop_shear_scatter_gauss_ldata_write_e_rms;
    data->ldata_required_columns = &_nc_galaxy_shape_pop_shear_scatter_gauss_ldata_e_rms_column;
  }
  else
  {
    data->ldata_read_row         = &_nc_galaxy_shape_pop_shear_scatter_gauss_ldata_noop;
    data->ldata_write_row        = &_nc_galaxy_shape_pop_shear_scatter_gauss_ldata_noop;
    data->ldata_required_columns = &_nc_galaxy_shape_pop_shear_scatter_gauss_ldata_no_columns;
  }
}

static void
_nc_galaxy_shape_pop_shear_scatter_gauss_prepare (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  /* Nothing to resolve: the width is read at evaluation. */
}

static inline gdouble
_nc_galaxy_shape_pop_shear_scatter_gauss_var (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, const gdouble var_noise_g, const gdouble eps_to_g)
{
  NcGalaxyShapePopShearScatterGauss *gspssg = NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS (gspss);
  const gdouble sigma_g                     = gspssg->sigma_from_catalog ? data->e_rms * eps_to_g : SIGMA;

  return sigma_g * sigma_g + var_noise_g;
}

static gdouble
_nc_galaxy_shape_pop_shear_scatter_gauss_eval_ln_p (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, const gdouble delta_g, const gdouble var_noise_g, const gdouble eps_to_g)
{
  const gdouble var = _nc_galaxy_shape_pop_shear_scatter_gauss_var (gspss, data, var_noise_g, eps_to_g);

  return -0.5 * delta_g * delta_g / var - 0.5 * log (2.0 * M_PI * var);
}

static gdouble
_nc_galaxy_shape_pop_shear_scatter_gauss_gen_delta (NcGalaxyShapePopShearScatter *gspss, NcGalaxyShapePopData *data, NcmRNG *rng, const gdouble var_noise_g, const gdouble eps_to_g)
{
  const gdouble var = _nc_galaxy_shape_pop_shear_scatter_gauss_var (gspss, data, var_noise_g, eps_to_g);

  return ncm_rng_gaussian_gen (rng, 0.0, sqrt (var));
}

static void
nc_galaxy_shape_pop_shear_scatter_gauss_class_init (NcGalaxyShapePopShearScatterGaussClass *klass)
{
  NcGalaxyShapePopShearScatterClass *gspss_class = NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_CLASS (klass);
  NcGalaxyShapePopClass *gsp_class               = NC_GALAXY_SHAPE_POP_CLASS (klass);
  GObjectClass *object_class                     = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class                     = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_shape_pop_shear_scatter_gauss_set_property;
  model_class->get_property = &_nc_galaxy_shape_pop_shear_scatter_gauss_get_property;
  object_class->finalize    = &_nc_galaxy_shape_pop_shear_scatter_gauss_finalize;

  ncm_model_class_set_name_nick (model_class, "Gaussian shear-estimator scatter", "ShearScatterGauss");
  ncm_model_class_add_params (model_class, NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxyShapePopShearScatterGauss:sigma-from-catalog:
   *
   * Whether the Gaussian width is taken per galaxy from the catalog column
   * "e_rms" (in ellipticity units, converted to shear units by the shape
   * factor) instead of the parameter sigma.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SIGMA_FROM_CATALOG,
                                   g_param_spec_boolean ("sigma-from-catalog",
                                                         "Sigma from catalog",
                                                         "Whether the width is the per-galaxy e_rms column",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyShapePopShearScatterGauss:sigma:
   *
   * The width $\sigma_g$ of the Gaussian scatter of the tangential shear
   * estimate, in reduced-shear units. Unused when
   * #NcGalaxyShapePopShearScatterGauss:sigma-from-catalog is %TRUE.
   *
   */
  ncm_model_class_set_sparam (model_class,
                              NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_SIGMA,
                              "\\sigma_g",
                              "sigma", 1.0e-3, 1.0, 1.0e-2,
                              NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_DEFAULT_PARAMS_ABSTOL,
                              NC_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS_DEFAULT_SIGMA,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  gsp_class->data_init = &_nc_galaxy_shape_pop_shear_scatter_gauss_data_init;
  gsp_class->prepare   = &_nc_galaxy_shape_pop_shear_scatter_gauss_prepare;

  gspss_class->eval_ln_p = &_nc_galaxy_shape_pop_shear_scatter_gauss_eval_ln_p;
  gspss_class->gen_delta = &_nc_galaxy_shape_pop_shear_scatter_gauss_gen_delta;
}

/**
 * nc_galaxy_shape_pop_shear_scatter_gauss_new:
 * @sigma_from_catalog: whether the width is the per-galaxy e_rms column
 *
 * Creates a new #NcGalaxyShapePopShearScatterGauss.
 *
 * Returns: (transfer full): a new #NcGalaxyShapePopShearScatterGauss.
 */
NcGalaxyShapePopShearScatterGauss *
nc_galaxy_shape_pop_shear_scatter_gauss_new (gboolean sigma_from_catalog)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_POP_SHEAR_SCATTER_GAUSS,
                       "sigma-from-catalog", sigma_from_catalog,
                       NULL);
}

/**
 * nc_galaxy_shape_pop_shear_scatter_gauss_ref:
 * @gspssg: a #NcGalaxyShapePopShearScatterGauss
 *
 * Increases the reference count of @gspssg by one.
 *
 * Returns: (transfer full): @gspssg.
 */
NcGalaxyShapePopShearScatterGauss *
nc_galaxy_shape_pop_shear_scatter_gauss_ref (NcGalaxyShapePopShearScatterGauss *gspssg)
{
  return g_object_ref (gspssg);
}

/**
 * nc_galaxy_shape_pop_shear_scatter_gauss_free:
 * @gspssg: a #NcGalaxyShapePopShearScatterGauss
 *
 * Decreases the reference count of @gspssg by one.
 *
 */
void
nc_galaxy_shape_pop_shear_scatter_gauss_free (NcGalaxyShapePopShearScatterGauss *gspssg)
{
  g_object_unref (gspssg);
}

/**
 * nc_galaxy_shape_pop_shear_scatter_gauss_clear:
 * @gspssg: a #NcGalaxyShapePopShearScatterGauss
 *
 * Decreases the reference count of *@gspssg by one, and sets the pointer
 * *@gspssg to NULL.
 *
 */
void
nc_galaxy_shape_pop_shear_scatter_gauss_clear (NcGalaxyShapePopShearScatterGauss **gspssg)
{
  g_clear_object (gspssg);
}

/**
 * nc_galaxy_shape_pop_shear_scatter_gauss_get_sigma_from_catalog:
 * @gspssg: a #NcGalaxyShapePopShearScatterGauss
 *
 * Returns: whether the width is the per-galaxy e_rms column.
 */
gboolean
nc_galaxy_shape_pop_shear_scatter_gauss_get_sigma_from_catalog (NcGalaxyShapePopShearScatterGauss *gspssg)
{
  return gspssg->sigma_from_catalog;
}

