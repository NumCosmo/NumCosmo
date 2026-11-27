/***************************************************************************
 *            nc_galaxy_wl_binned_estimator.c
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_binned_estimator.c
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyWLBinnedEstimator:
 *
 * Abstract estimator of a binned cluster weak-lensing profile.
 *
 * An estimator does two things, which must agree with each other:
 *
 * - it builds a binned profile (#NcGalaxyWLObsBinned) from a per-galaxy catalog
 *   (#NcGalaxyWLObs), nc_galaxy_wl_binned_estimator_build();
 * - it predicts the binned profile for the models of a #NcmMSet,
 *   nc_galaxy_wl_binned_estimator_eval_model().
 *
 * Each concrete estimator implements a classic method of the literature, with its own
 * weights, calibration and approximations to the model: #NcGalaxyWLBinnedEstimatorDeltaSigma
 * (the excess surface mass density of Umetsu et al. 2020) and #NcGalaxyWLBinnedEstimatorGt
 * (the reduced tangential shear, with the model of Schrabback et al. 2018).
 *
 * The build uses the models of the #NcmMSet only for the geometry: the projected radius and
 * the tangential rotation of each galaxy (#NcHaloPosition) and its critical densities
 * (#NcWLSurfaceMassDensity), averaged over its redshift density. That density is the one
 * of #NcGalaxyWLBinnedEstimator:redshift-factor when it is set, the same density
 * #NcDataClusterWLFactor integrates; otherwise the #NcGalaxyPz the catalog holds for the
 * row (nc_galaxy_wl_obs_peek_pz()), or its column "z" when it has none. The catalog must have the columns "ra", "dec",
 * "epsilon_obs_1", "epsilon_obs_2" and "std_noise"; the columns "m", "c1", "c2" are taken as
 * zero when absent, and the column "e_rms" is replaced by #NcGalaxyWLBinnedEstimator:e-rms
 * when that property is non-negative. Radii are proper Mpc (NumCosmo's convention; Umetsu
 * et al. 2020 use comoving units). Bins without galaxies are dropped.
 *
 * Besides its classic model, every estimator supports the model
 * nc_galaxy_wl_binned_estimator_eval_model_exact(): the expected value of the estimator,
 * the weighted mean of the reduced shear of each galaxy at its own radius, averaged over its
 * own redshift density, with the weights of the estimator (see #NcGalaxyWLObsBinned). It
 * needs the galaxies of the profile, so it is not available for a profile built from
 * external arrays. It validates the approximations of the classic model.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_wl_binned_estimator.h"
#include "nc/lss/galaxy/nc_galaxy_wl_binned_estimator_private.h"
#include "nc/lss/galaxy/nc_galaxy_position_factor.h"
#include "nc/lss/galaxy/nc_galaxy_redshift_factor.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor.h"
#include "nc/lss/halo/nc_halo_position.h"
#include "nc/lss/halo/nc_halo_density_profile.h"
#include "nc/lss/wl/nc_wl_surface_mass_density.h"
#include "nc/lss/wl/nc_wl_ellipticity.h"
#include "nc/background/nc_hicosmo.h"
#include "ncm/integration/ncm_integrate.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

#define _NC_WLBE_GL_ORDER (5)

typedef struct _NcGalaxyWLBinnedEstimatorPrivate
{
  gdouble e_rms;
  NcGalaxyRedshiftFactor *redshift_factor;
  guint pz_intervals;
  guint exact_intervals;

  /* Cache of eval_model_exact(): the redshift nodes of every member */
  NcmCatalog *exact_members;
  guint64 exact_members_key;
  guint64 exact_crit_hash;
  guint64 exact_z_hash;
  NcGalaxyRedshiftFactor *exact_zf; /* compared by address only, never dereferenced */
  guint exact_n_int;
  GArray *exact_offset; /* guint, members + 1 */
  GArray *exact_s;      /* Sigma_cr^-1 at the nodes */
  GArray *exact_w;      /* weight of the nodes, P(z) / norm included */
} NcGalaxyWLBinnedEstimatorPrivate;

enum
{
  PROP_0,
  PROP_E_RMS,
  PROP_PZ_INTERVALS,
  PROP_EXACT_INTERVALS,
  PROP_REDSHIFT_FACTOR,
  PROP_LEN,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcGalaxyWLBinnedEstimator, nc_galaxy_wl_binned_estimator, G_TYPE_OBJECT)

static void
nc_galaxy_wl_binned_estimator_init (NcGalaxyWLBinnedEstimator *est)
{
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);

  self->e_rms             = -1.0;
  self->redshift_factor   = NULL;
  self->exact_z_hash      = 0;
  self->exact_zf          = NULL;
  self->pz_intervals      = NC_GALAXY_WL_BINNED_ESTIMATOR_DEFAULT_PZ_INTERVALS;
  self->exact_intervals   = NC_GALAXY_WL_BINNED_ESTIMATOR_DEFAULT_EXACT_INTERVALS;
  self->exact_members     = NULL;
  self->exact_members_key = NCM_CATALOG_KEY_PENDING;
  self->exact_crit_hash   = 0;
  self->exact_n_int       = 0;
  self->exact_offset      = g_array_new (FALSE, FALSE, sizeof (guint));
  self->exact_s           = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->exact_w           = g_array_new (FALSE, FALSE, sizeof (gdouble));
}

static void
_nc_galaxy_wl_binned_estimator_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLBinnedEstimator *est = NC_GALAXY_WL_BINNED_ESTIMATOR (object);

  switch (prop_id)
  {
    case PROP_E_RMS:
      nc_galaxy_wl_binned_estimator_set_e_rms (est, g_value_get_double (value));
      break;
    case PROP_PZ_INTERVALS:
      nc_galaxy_wl_binned_estimator_set_pz_intervals (est, g_value_get_uint (value));
      break;
    case PROP_EXACT_INTERVALS:
      nc_galaxy_wl_binned_estimator_set_exact_intervals (est, g_value_get_uint (value));
      break;
    case PROP_REDSHIFT_FACTOR:
      nc_galaxy_wl_binned_estimator_set_redshift_factor (est, g_value_get_object (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_wl_binned_estimator_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLBinnedEstimator *est = NC_GALAXY_WL_BINNED_ESTIMATOR (object);

  switch (prop_id)
  {
    case PROP_E_RMS:
      g_value_set_double (value, nc_galaxy_wl_binned_estimator_get_e_rms (est));
      break;
    case PROP_PZ_INTERVALS:
      g_value_set_uint (value, nc_galaxy_wl_binned_estimator_get_pz_intervals (est));
      break;
    case PROP_EXACT_INTERVALS:
      g_value_set_uint (value, nc_galaxy_wl_binned_estimator_get_exact_intervals (est));
      break;
    case PROP_REDSHIFT_FACTOR:
      g_value_set_object (value, nc_galaxy_wl_binned_estimator_peek_redshift_factor (est));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_wl_binned_estimator_dispose (GObject *object)
{
  NcGalaxyWLBinnedEstimator *est                = NC_GALAXY_WL_BINNED_ESTIMATOR (object);
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);

  ncm_catalog_clear (&self->exact_members);
  nc_galaxy_redshift_factor_clear (&self->redshift_factor);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_binned_estimator_parent_class)->dispose (object);
}

static void
_nc_galaxy_wl_binned_estimator_finalize (GObject *object)
{
  NcGalaxyWLBinnedEstimator *est                = NC_GALAXY_WL_BINNED_ESTIMATOR (object);
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);

  g_clear_pointer (&self->exact_offset, g_array_unref);
  g_clear_pointer (&self->exact_s, g_array_unref);
  g_clear_pointer (&self->exact_w, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_binned_estimator_parent_class)->finalize (object);
}

/* LCOV_EXCL_START */
static GStrv
_nc_galaxy_wl_binned_estimator_extra_columns (NcGalaxyWLBinnedEstimator *est)
{
  return g_new0 (gchar *, 1);
}

static void
_nc_galaxy_wl_binned_estimator_build_bins (NcGalaxyWLBinnedEstimator *est, const NcGalaxyWLBinnedEstimatorSources *src, NcGalaxyWLObsBinned *binned, GArray *factor)
{
  g_error ("_nc_galaxy_wl_binned_estimator_build_bins: method not implemented.");
}

static void
_nc_galaxy_wl_binned_estimator_eval_model (NcGalaxyWLBinnedEstimator *est, NcmMSet *mset, NcGalaxyWLObsBinned *binned, NcmVector *model)
{
  g_error ("_nc_galaxy_wl_binned_estimator_eval_model: method not implemented.");
}

/* LCOV_EXCL_STOP */

static void
nc_galaxy_wl_binned_estimator_class_init (NcGalaxyWLBinnedEstimatorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_wl_binned_estimator_set_property;
  object_class->get_property = &_nc_galaxy_wl_binned_estimator_get_property;
  object_class->dispose      = &_nc_galaxy_wl_binned_estimator_dispose;
  object_class->finalize     = &_nc_galaxy_wl_binned_estimator_finalize;

  /**
   * NcGalaxyWLBinnedEstimator:e-rms:
   *
   * The per-component intrinsic ellipticity rms of every galaxy, in the convention of the
   * catalog. A negative value, the default, reads it from the column "e_rms" of the catalog.
   */
  g_object_class_install_property (object_class,
                                   PROP_E_RMS,
                                   g_param_spec_double ("e-rms",
                                                        NULL,
                                                        "Intrinsic ellipticity rms",
                                                        -1.0, 1.0, -1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLBinnedEstimator:pz-intervals:
   *
   * The number of intervals of the Gauss-Legendre rule (5 nodes each) of the averages over
   * the redshift density of a galaxy, from the cluster redshift to the end of its support,
   * used to build the profile.
   */
  g_object_class_install_property (object_class,
                                   PROP_PZ_INTERVALS,
                                   g_param_spec_uint ("pz-intervals",
                                                      NULL,
                                                      "Intervals of the p(z) averages of the build",
                                                      1, G_MAXUINT, NC_GALAXY_WL_BINNED_ESTIMATOR_DEFAULT_PZ_INTERVALS,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyWLBinnedEstimator:exact-intervals:
   *
   * The number of intervals of the Gauss-Legendre rule (5 nodes each) of the averages over
   * the redshift density of nc_galaxy_wl_binned_estimator_eval_model_exact(). The nodes of
   * every galaxy are kept in memory.
   */
  g_object_class_install_property (object_class,
                                   PROP_EXACT_INTERVALS,
                                   g_param_spec_uint ("exact-intervals",
                                                      NULL,
                                                      "Intervals of the p(z) averages of the exact model",
                                                      1, G_MAXUINT, NC_GALAXY_WL_BINNED_ESTIMATOR_DEFAULT_EXACT_INTERVALS,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  /**
   * NcGalaxyWLBinnedEstimator:redshift-factor:
   *
   * The redshift factor whose per-galaxy density $p(z_\mathrm{phot}, z \mid I)$, normalized
   * by nc_galaxy_redshift_factor_norm(), the averages over redshift use: the same density
   * #NcDataClusterWLFactor integrates. When %NULL, the default, the averages use the
   * #NcGalaxyPz of each row (nc_galaxy_wl_obs_peek_pz()) or its column "z".
   */
  g_object_class_install_property (object_class,
                                   PROP_REDSHIFT_FACTOR,
                                   g_param_spec_object ("redshift-factor",
                                                        NULL,
                                                        "Redshift factor of the averages over redshift",
                                                        NC_TYPE_GALAXY_REDSHIFT_FACTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->name          = "abstract";
  klass->extra_columns = &_nc_galaxy_wl_binned_estimator_extra_columns;
  klass->build_bins    = &_nc_galaxy_wl_binned_estimator_build_bins;
  klass->eval_model    = &_nc_galaxy_wl_binned_estimator_eval_model;
}

/**
 * nc_galaxy_wl_binned_estimator_ref:
 * @est: a #NcGalaxyWLBinnedEstimator
 *
 * Increases the reference count of @est by one.
 *
 * Returns: (transfer full): @est.
 */
NcGalaxyWLBinnedEstimator *
nc_galaxy_wl_binned_estimator_ref (NcGalaxyWLBinnedEstimator *est)
{
  return g_object_ref (est);
}

/**
 * nc_galaxy_wl_binned_estimator_free:
 * @est: a #NcGalaxyWLBinnedEstimator
 *
 * Decreases the reference count of @est by one.
 */
void
nc_galaxy_wl_binned_estimator_free (NcGalaxyWLBinnedEstimator *est)
{
  g_object_unref (est);
}

/**
 * nc_galaxy_wl_binned_estimator_clear:
 * @est: a #NcGalaxyWLBinnedEstimator
 *
 * Decreases the reference count of *@est by one and sets *@est to %NULL.
 */
void
nc_galaxy_wl_binned_estimator_clear (NcGalaxyWLBinnedEstimator **est)
{
  g_clear_object (est);
}

/**
 * nc_galaxy_wl_binned_estimator_get_name:
 * @est: a #NcGalaxyWLBinnedEstimator
 *
 * Returns: (transfer none): the name of the estimator, stored in the profiles it builds.
 */
const gchar *
nc_galaxy_wl_binned_estimator_get_name (NcGalaxyWLBinnedEstimator *est)
{
  return NC_GALAXY_WL_BINNED_ESTIMATOR_GET_CLASS (est)->name;
}

/**
 * nc_galaxy_wl_binned_estimator_set_e_rms:
 * @est: a #NcGalaxyWLBinnedEstimator
 * @e_rms: the intrinsic ellipticity rms, or a negative value
 *
 * Sets #NcGalaxyWLBinnedEstimator:e-rms.
 */
void
nc_galaxy_wl_binned_estimator_set_e_rms (NcGalaxyWLBinnedEstimator *est, const gdouble e_rms)
{
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);

  self->e_rms = e_rms;
}

/**
 * nc_galaxy_wl_binned_estimator_get_e_rms:
 * @est: a #NcGalaxyWLBinnedEstimator
 *
 * Returns: #NcGalaxyWLBinnedEstimator:e-rms.
 */
gdouble
nc_galaxy_wl_binned_estimator_get_e_rms (NcGalaxyWLBinnedEstimator *est)
{
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);

  return self->e_rms;
}

/**
 * nc_galaxy_wl_binned_estimator_set_pz_intervals:
 * @est: a #NcGalaxyWLBinnedEstimator
 * @pz_intervals: the number of intervals
 *
 * Sets #NcGalaxyWLBinnedEstimator:pz-intervals.
 */
void
nc_galaxy_wl_binned_estimator_set_pz_intervals (NcGalaxyWLBinnedEstimator *est, const guint pz_intervals)
{
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);

  g_assert_cmpuint (pz_intervals, >, 0);
  self->pz_intervals = pz_intervals;
}

/**
 * nc_galaxy_wl_binned_estimator_get_pz_intervals:
 * @est: a #NcGalaxyWLBinnedEstimator
 *
 * Returns: #NcGalaxyWLBinnedEstimator:pz-intervals.
 */
guint
nc_galaxy_wl_binned_estimator_get_pz_intervals (NcGalaxyWLBinnedEstimator *est)
{
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);

  return self->pz_intervals;
}

/**
 * nc_galaxy_wl_binned_estimator_set_exact_intervals:
 * @est: a #NcGalaxyWLBinnedEstimator
 * @exact_intervals: the number of intervals
 *
 * Sets #NcGalaxyWLBinnedEstimator:exact-intervals.
 */
void
nc_galaxy_wl_binned_estimator_set_exact_intervals (NcGalaxyWLBinnedEstimator *est, const guint exact_intervals)
{
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);

  g_assert_cmpuint (exact_intervals, >, 0);
  self->exact_intervals = exact_intervals;
}

/**
 * nc_galaxy_wl_binned_estimator_get_exact_intervals:
 * @est: a #NcGalaxyWLBinnedEstimator
 *
 * Returns: #NcGalaxyWLBinnedEstimator:exact-intervals.
 */
guint
nc_galaxy_wl_binned_estimator_get_exact_intervals (NcGalaxyWLBinnedEstimator *est)
{
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);

  return self->exact_intervals;
}

/**
 * nc_galaxy_wl_binned_estimator_set_redshift_factor:
 * @est: a #NcGalaxyWLBinnedEstimator
 * @redshift_factor: (nullable): a #NcGalaxyRedshiftFactor
 *
 * Sets #NcGalaxyWLBinnedEstimator:redshift-factor, referencing @redshift_factor.
 */
void
nc_galaxy_wl_binned_estimator_set_redshift_factor (NcGalaxyWLBinnedEstimator *est, NcGalaxyRedshiftFactor *redshift_factor)
{
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);

  nc_galaxy_redshift_factor_clear (&self->redshift_factor);

  if (redshift_factor != NULL)
    self->redshift_factor = nc_galaxy_redshift_factor_ref (redshift_factor);
}

/**
 * nc_galaxy_wl_binned_estimator_peek_redshift_factor:
 * @est: a #NcGalaxyWLBinnedEstimator
 *
 * Returns: (transfer none) (nullable): #NcGalaxyWLBinnedEstimator:redshift-factor.
 */
NcGalaxyRedshiftFactor *
nc_galaxy_wl_binned_estimator_peek_redshift_factor (NcGalaxyWLBinnedEstimator *est)
{
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);

  return self->redshift_factor;
}

/**
 * nc_galaxy_wl_binned_estimator_required_columns:
 * @est: a #NcGalaxyWLBinnedEstimator
 *
 * Gets the columns the per-galaxy catalog must have. The optional columns ("m", "c1", "c2",
 * "e_rms" and "z" for the rows without a redshift density) are not listed.
 *
 * Returns: (element-type utf8) (transfer full): the required columns.
 */
GList *
nc_galaxy_wl_binned_estimator_required_columns (NcGalaxyWLBinnedEstimator *est)
{
  GList *columns = NULL;

  columns = g_list_append (columns, g_strdup (NC_GALAXY_POSITION_FACTOR_COL_RA));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_POSITION_FACTOR_COL_DEC));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_1));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_2));
  columns = g_list_append (columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_COL_STD_NOISE));

  return columns;
}

typedef struct _NcWLBEModels
{
  NcHICosmo *cosmo;
  NcHaloPosition *hp;
  NcWLSurfaceMassDensity *smd;
  NcHaloDensityProfile *dp;
  gdouble z_cl;

  /* Set only when the estimator has a redshift factor */
  NcmMSet *mset;
  NcGalaxyRedshiftFactor *zf;
  NcGalaxyRedshiftFactorData *z_data;
  NcGalaxyRedshiftFactorIntegrand *z_integ;
} NcWLBEModels;

static void
_nc_wlbe_models (NcmMSet *mset, NcWLBEModels *models)
{
  models->cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  models->hp    = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  models->smd   = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  models->dp    = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));

  if ((models->cosmo == NULL) || (models->hp == NULL) || (models->smd == NULL))
    g_error ("nc_galaxy_wl_binned_estimator: the mset must hold a NcHICosmo, a NcHaloPosition and a NcWLSurfaceMassDensity.");

  nc_halo_position_prepare_if_needed (models->hp, models->cosmo);
  nc_wl_surface_mass_density_prepare_if_needed (models->smd, models->cosmo);

  models->z_cl = nc_halo_position_get_redshift (models->hp);

  models->mset    = mset;
  models->zf      = NULL;
  models->z_data  = NULL;
  models->z_integ = NULL;
}

/* Prepares the redshift factor of the estimator, if any, for _nc_wlbe_append_nodes() */
static void
_nc_wlbe_models_zf_begin (NcGalaxyWLBinnedEstimator *est, NcWLBEModels *models)
{
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);

  if (self->redshift_factor == NULL)
    return;

  models->zf = self->redshift_factor;
  nc_galaxy_redshift_factor_prepare (models->zf, models->mset);
  models->z_data  = nc_galaxy_redshift_factor_data_new (models->zf, models->mset);
  models->z_integ = nc_galaxy_redshift_factor_integ (models->zf, models->mset, FALSE);
  nc_galaxy_redshift_factor_integrand_prepare (models->z_integ, models->mset);
}

static void
_nc_wlbe_models_zf_end (NcWLBEModels *models)
{
  g_clear_pointer (&models->z_data, nc_galaxy_redshift_factor_data_unref);
  g_clear_pointer (&models->z_integ, nc_galaxy_redshift_factor_integrand_free);
  models->zf = NULL;
}

/* Appends the Gauss-Legendre nodes of [z_min, z_max] with weight p (z) dz / norm */
static void
_nc_wlbe_append_gl (NcWLBEModels *models, gdouble z_min, const gdouble z_max, const guint n_int,
                    gdouble (*p_eval) (gpointer arg, const gdouble z), gpointer arg, const gdouble norm,
                    GArray *s, GArray *w)
{
  z_min = GSL_MAX (z_min, models->z_cl);

  if ((z_max > z_min) && (norm > 0.0))
  {
    const gdouble dz = (z_max - z_min) / n_int;
    const gdouble *gl_x, *gl_w;
    guint k, l;

    ncm_integral_gauss_legendre_rule (_NC_WLBE_GL_ORDER, &gl_x, &gl_w);

    for (k = 0; k < n_int; k++)
    {
      const gdouble za = z_min + k * dz;

      for (l = 0; l < _NC_WLBE_GL_ORDER; l++)
      {
        const gdouble z  = za + 0.5 * dz * (1.0 + gl_x[l]);
        const gdouble wz = 0.5 * dz * gl_w[l] * p_eval (arg, z) / norm;

        if (wz == 0.0)
          continue;

        {
          const gdouble inv_sigma = 1.0 / nc_wl_surface_mass_density_sigma_critical (models->smd, models->cosmo, z, models->z_cl, models->z_cl);

          g_array_append_val (w, wz);
          g_array_append_val (s, inv_sigma);
        }
      }
    }
  }
}

static gdouble
_nc_wlbe_pz_eval (gpointer arg, const gdouble z)
{
  return nc_galaxy_pz_eval_p (NC_GALAXY_PZ (arg), z);
}

static gdouble
_nc_wlbe_zf_eval (gpointer arg, const gdouble z)
{
  NcWLBEModels *models = (NcWLBEModels *) arg;

  return nc_galaxy_redshift_factor_integrand_eval (models->z_integ, z, models->z_data);
}

/*
 * Appends to s and w the nodes of the average over the redshift density of galaxy i above
 * the cluster: s = Sigma_cr^-1 (z) and w = P (z) dz / norm, or one node with w = 1 for a
 * galaxy with a point redshift behind the cluster. The density is the integrand of the
 * redshift factor of the estimator when it has one (models->zf, prepared by
 * _nc_wlbe_models_zf_begin()), otherwise the NcGalaxyPz of the row or its "z" column.
 */
static void
_nc_wlbe_append_nodes (NcWLBEModels *models, NcGalaxyWLObs *obs, const guint i, const guint n_int, GArray *s, GArray *w)
{
  if (models->zf != NULL)
  {
    gdouble z_min, z_max, norm;

    nc_galaxy_redshift_factor_data_read_row (models->z_data, obs, i);
    nc_galaxy_redshift_factor_update_data (models->zf, models->z_data);
    nc_galaxy_redshift_factor_get_integ_lim (models->zf, models->mset, models->z_data, &z_min, &z_max);
    norm = nc_galaxy_redshift_factor_norm (models->zf, models->mset, models->z_data);

    _nc_wlbe_append_gl (models, z_min, z_max, n_int, &_nc_wlbe_zf_eval, models, norm, s, w);
  }
  else
  {
    NcGalaxyPz *pz = nc_galaxy_wl_obs_peek_pz (obs, i);

    if (pz != NULL)
    {
      gdouble z_min, z_max;

      nc_galaxy_pz_get_range (pz, &z_min, &z_max);
      _nc_wlbe_append_gl (models, z_min, z_max, n_int, &_nc_wlbe_pz_eval, pz, nc_galaxy_pz_get_norm (pz), s, w);
    }
    else
    {
      const gdouble z = nc_galaxy_wl_obs_get (obs, NC_GALAXY_REDSHIFT_FACTOR_COL_Z, i, NULL);

      if (z > models->z_cl)
      {
        const gdouble inv_sigma = 1.0 / nc_wl_surface_mass_density_sigma_critical (models->smd, models->cosmo, z, models->z_cl, models->z_cl);
        const gdouble one       = 1.0;

        g_array_append_val (s, inv_sigma);
        g_array_append_val (w, one);
      }
    }
  }
}

static gdouble
_nc_wlbe_mean_inv_sigma_crit (NcWLBEModels *models, NcGalaxyWLObs *obs, const guint i, const guint n_int, GArray *s, GArray *w, gdouble *mean2)
{
  gdouble mean1 = 0.0;
  guint k;

  g_array_set_size (s, 0);
  g_array_set_size (w, 0);
  _nc_wlbe_append_nodes (models, obs, i, n_int, s, w);

  *mean2 = 0.0;

  for (k = 0; k < s->len; k++)
  {
    const gdouble s_k = g_array_index (s, gdouble, k);
    const gdouble w_k = g_array_index (w, gdouble, k);

    mean1  += w_k * s_k;
    *mean2 += w_k * s_k * s_k;
  }

  return mean1;
}

/**
 * nc_galaxy_wl_binned_estimator_mean_inv_sigma_crit:
 * @est: a #NcGalaxyWLBinnedEstimator
 * @mset: a #NcmMSet
 * @obs: a #NcGalaxyWLObs
 * @i: the row of the galaxy
 * @mean_inv_sigma_crit2: (out): the mean of $\Sigma_\mathrm{cr}^{-2}$
 *
 * Computes the mean of the inverse critical density of galaxy @i over its redshift density
 * (Umetsu et al. 2020, Eq. 6), $\langle\Sigma_\mathrm{cr}^{-1}\rangle = \int P(z)
 * \Sigma_\mathrm{cr}^{-1}(z)\,\mathrm{d}z / \int P(z)\,\mathrm{d}z$, with
 * $\Sigma_\mathrm{cr}^{-1} = 0$ in front of the cluster, by the rule of
 * #NcGalaxyWLBinnedEstimator:pz-intervals.
 *
 * Returns: $\langle\Sigma_\mathrm{cr}^{-1}\rangle$ in $\mathrm{Mpc}^2 / M_\odot$.
 */
gdouble
nc_galaxy_wl_binned_estimator_mean_inv_sigma_crit (NcGalaxyWLBinnedEstimator *est, NcmMSet *mset, NcGalaxyWLObs *obs, const guint i, gdouble *mean_inv_sigma_crit2)
{
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);
  GArray *s                                     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *w                                     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcWLBEModels models;
  gdouble mean1;

  _nc_wlbe_models (mset, &models);
  _nc_wlbe_models_zf_begin (est, &models);
  mean1 = _nc_wlbe_mean_inv_sigma_crit (&models, obs, i, self->pz_intervals, s, w, mean_inv_sigma_crit2);
  _nc_wlbe_models_zf_end (&models);

  g_array_unref (s);
  g_array_unref (w);

  return mean1;
}

static gdouble
_nc_wlbe_get_opt (NcGalaxyWLObs *obs, const gchar *col, const guint i)
{
  if (ncm_catalog_has_column (NCM_CATALOG (obs), col))
    return nc_galaxy_wl_obs_get (obs, col, i, NULL);
  else
    return 0.0;
}

static NcGalaxyWLBinnedEstimatorSources *
_nc_wlbe_sources_new (void)
{
  NcGalaxyWLBinnedEstimatorSources *src = g_new0 (NcGalaxyWLBinnedEstimatorSources, 1);

  src->row             = g_array_new (FALSE, FALSE, sizeof (guint));
  src->bin             = g_array_new (FALSE, FALSE, sizeof (guint));
  src->radius          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  src->e_t             = g_array_new (FALSE, FALSE, sizeof (gdouble));
  src->e_x             = g_array_new (FALSE, FALSE, sizeof (gdouble));
  src->c_t             = g_array_new (FALSE, FALSE, sizeof (gdouble));
  src->c_x             = g_array_new (FALSE, FALSE, sizeof (gdouble));
  src->c1              = g_array_new (FALSE, FALSE, sizeof (gdouble));
  src->c2              = g_array_new (FALSE, FALSE, sizeof (gdouble));
  src->rot_re          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  src->rot_im          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  src->m               = g_array_new (FALSE, FALSE, sizeof (gdouble));
  src->std_noise       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  src->e_rms           = g_array_new (FALSE, FALSE, sizeof (gdouble));
  src->inv_sigma_crit  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  src->inv_sigma_crit2 = g_array_new (FALSE, FALSE, sizeof (gdouble));

  return src;
}

static void
_nc_wlbe_sources_free (NcGalaxyWLBinnedEstimatorSources *src)
{
  g_array_unref (src->row);
  g_array_unref (src->bin);
  g_array_unref (src->radius);
  g_array_unref (src->e_t);
  g_array_unref (src->e_x);
  g_array_unref (src->c_t);
  g_array_unref (src->c_x);
  g_array_unref (src->c1);
  g_array_unref (src->c2);
  g_array_unref (src->rot_re);
  g_array_unref (src->rot_im);
  g_array_unref (src->m);
  g_array_unref (src->std_noise);
  g_array_unref (src->e_rms);
  g_array_unref (src->inv_sigma_crit);
  g_array_unref (src->inv_sigma_crit2);
  g_free (src);
}

/*
 * Drops the bins whose value or variance is not finite (for instance a bin of galaxies in
 * front of the cluster only, which have zero weight). On return compact[b] is the index of
 * bin b of src in the returned profile, or G_MAXUINT for a dropped bin. Takes binned.
 */
static NcGalaxyWLObsBinned *
_nc_wlbe_drop_invalid (NcGalaxyWLObsBinned *binned, const NcGalaxyWLBinnedEstimatorSources *src, GArray *compact)
{
  const guint n_bins = nc_galaxy_wl_obs_binned_len (binned);
  NcmMatrix *cov     = nc_galaxy_wl_obs_binned_peek_cov (binned);
  guint n_valid      = 0;
  guint b, b2;

  g_array_set_size (compact, n_bins);

  for (b = 0; b < n_bins; b++)
  {
    const gdouble value = nc_galaxy_wl_obs_binned_get (binned, NC_GALAXY_WL_OBS_BINNED_COL_VALUE, b);
    const gdouble var   = ncm_matrix_get (cov, b, b);

    if (gsl_finite (value) && gsl_finite (var) && (var > 0.0))
      g_array_index (compact, guint, b) = n_valid++;
    else
      g_array_index (compact, guint, b) = G_MAXUINT;
  }

  if (n_valid == n_bins)
    return binned;

  if (n_valid == 0)
    g_error ("nc_galaxy_wl_binned_estimator_build: no bin has a finite value and variance.");

  {
    GStrv cols                 = ncm_catalog_peek_columns (NCM_CATALOG (binned));
    NcGalaxyWLObsBinned *valid = g_object_new (NC_TYPE_GALAXY_WL_OBS_BINNED, "columns", cols, "len", n_valid, NULL);
    NcmMatrix *cov_valid       = ncm_matrix_new (n_valid, n_valid);
    guint c;

    for (b = 0; b < n_bins; b++)
    {
      const guint bc = g_array_index (compact, guint, b);

      if (bc == G_MAXUINT)
        continue;

      for (c = 0; cols[c] != NULL; c++)
        nc_galaxy_wl_obs_binned_set (valid, cols[c], bc, nc_galaxy_wl_obs_binned_get (binned, cols[c], b));

      for (b2 = 0; b2 < n_bins; b2++)
      {
        const guint bc2 = g_array_index (compact, guint, b2);

        if (bc2 != G_MAXUINT)
          ncm_matrix_set (cov_valid, bc, bc2, ncm_matrix_get (cov, b, b2));
      }
    }

    nc_galaxy_wl_obs_binned_set_cov (valid, cov_valid);
    nc_galaxy_wl_obs_binned_set_z_cluster (valid, nc_galaxy_wl_obs_binned_get_z_cluster (binned));
    nc_galaxy_wl_obs_binned_set_estimator (valid, nc_galaxy_wl_obs_binned_peek_estimator (binned));

    ncm_matrix_free (cov_valid);
    nc_galaxy_wl_obs_binned_free (binned);

    return valid;
  }
}

/**
 * nc_galaxy_wl_binned_estimator_build:
 * @est: a #NcGalaxyWLBinnedEstimator
 * @mset: the #NcmMSet of the geometry
 * @obs: the per-galaxy catalog
 * @edges: (element-type gdouble): the bin edges in projected radius (proper Mpc), increasing
 *
 * Builds the binned profile of @obs: the galaxies with projected radius in
 * $[e_0, e_n)$ enter the bin $[e_b, e_{b+1})$ that holds them; bins without galaxies are
 * dropped. The models of @mset fix the geometry only (see the class description); the
 * profile keeps the cluster redshift it was built with.
 *
 * Returns: (transfer full): the binned profile.
 */
NcGalaxyWLObsBinned *
nc_galaxy_wl_binned_estimator_build (NcGalaxyWLBinnedEstimator *est, NcmMSet *mset, NcGalaxyWLObs *obs, GArray *edges)
{
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);
  NcGalaxyWLBinnedEstimatorClass *klass         = NC_GALAXY_WL_BINNED_ESTIMATOR_GET_CLASS (est);
  const guint n_edges                           = edges->len;
  const guint n_rows                            = nc_galaxy_wl_obs_len (obs);
  const NcWLEllipticityFrame coord              = nc_galaxy_wl_obs_get_coord (obs);
  const gboolean has_e_rms                      = ncm_catalog_has_column (NCM_CATALOG (obs), NC_GALAXY_WL_BINNED_ESTIMATOR_COL_E_RMS);
  NcGalaxyWLBinnedEstimatorSources *src         = _nc_wlbe_sources_new ();
  GArray *s                                     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *w                                     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *count                                 = g_array_new (FALSE, TRUE, sizeof (guint));
  GArray *compact                               = g_array_new (FALSE, TRUE, sizeof (guint));
  NcGalaxyWLObsBinned *binned;
  NcWLBEModels models;
  guint i, b;

  if (n_edges < 2)
    g_error ("nc_galaxy_wl_binned_estimator_build: at least two bin edges are needed.");

  for (b = 0; b + 1 < n_edges; b++)
  {
    if (!(g_array_index (edges, gdouble, b + 1) > g_array_index (edges, gdouble, b)))
      g_error ("nc_galaxy_wl_binned_estimator_build: the bin edges must increase.");
  }

  if ((self->e_rms < 0.0) && !has_e_rms)
    g_error ("nc_galaxy_wl_binned_estimator_build: the catalog has no column `%s' and the property e-rms is not set.",
             NC_GALAXY_WL_BINNED_ESTIMATOR_COL_E_RMS);

  _nc_wlbe_models (mset, &models);
  _nc_wlbe_models_zf_begin (est, &models);

  switch (nc_galaxy_wl_obs_get_ellip_conv (obs))
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      src->shear_scale      = 2.0;
      src->use_responsivity = TRUE;
      break;
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      src->shear_scale      = 1.0;
      src->use_responsivity = FALSE;
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
      break;                   /* LCOV_EXCL_LINE */
  }

  src->sigma_crit_inf = nc_wl_surface_mass_density_sigma_critical_infinity (models.smd, models.cosmo, models.z_cl, models.z_cl);
  g_array_set_size (count, n_edges - 1);

  for (i = 0; i < n_rows; i++)
  {
    const gdouble ra  = nc_galaxy_wl_obs_get (obs, NC_GALAXY_POSITION_FACTOR_COL_RA, i, NULL);
    const gdouble dec = nc_galaxy_wl_obs_get (obs, NC_GALAXY_POSITION_FACTOR_COL_DEC, i, NULL);
    gdouble theta, phi, radius;

    nc_halo_position_polar_angles (models.hp, ra, dec, &theta, &phi);
    radius = nc_halo_position_projected_radius (models.hp, models.cosmo, theta);

    if ((radius < g_array_index (edges, gdouble, 0)) || (radius >= g_array_index (edges, gdouble, n_edges - 1)))
      continue;

    for (b = 0; b + 2 < n_edges; b++)
    {
      if (radius < g_array_index (edges, gdouble, b + 1))
        break;
    }

    {
      const gdouble e1        = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_1, i, NULL);
      const gdouble e2        = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_2, i, NULL);
      const gdouble std_noise = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SHAPE_FACTOR_COL_STD_NOISE, i, NULL);
      const gdouble m         = _nc_wlbe_get_opt (obs, NC_GALAXY_SHAPE_FACTOR_COL_M, i);
      const gdouble c1        = _nc_wlbe_get_opt (obs, NC_GALAXY_SHAPE_FACTOR_COL_C1, i);
      const gdouble c2        = _nc_wlbe_get_opt (obs, NC_GALAXY_SHAPE_FACTOR_COL_C2, i);
      const gdouble e_rms     = (self->e_rms >= 0.0) ? self->e_rms : nc_galaxy_wl_obs_get (obs, NC_GALAXY_WL_BINNED_ESTIMATOR_COL_E_RMS, i, NULL);
      const gdouble phi_f     = nc_wl_ellipticity_celestial_to_frame_angle (coord, phi);
      const gdouble rot_re    = cos (2.0 * phi_f);
      const gdouble rot_im    = sin (2.0 * phi_f);
      /* (e1 + i e2) exp (-2 i phi) */
      const gdouble e_t = e1 * rot_re + e2 * rot_im;
      const gdouble e_x = e2 * rot_re - e1 * rot_im;
      const gdouble c_t = c1 * rot_re + c2 * rot_im;
      const gdouble c_x = c2 * rot_re - c1 * rot_im;
      gdouble mean2;
      const gdouble mean1 = _nc_wlbe_mean_inv_sigma_crit (&models, obs, i, self->pz_intervals, s, w, &mean2);

      g_array_append_val (src->row, i);
      g_array_append_val (src->bin, b);
      g_array_append_val (src->radius, radius);
      g_array_append_val (src->e_t, e_t);
      g_array_append_val (src->e_x, e_x);
      g_array_append_val (src->c_t, c_t);
      g_array_append_val (src->c_x, c_x);
      g_array_append_val (src->c1, c1);
      g_array_append_val (src->c2, c2);
      g_array_append_val (src->rot_re, rot_re);
      g_array_append_val (src->rot_im, rot_im);
      g_array_append_val (src->m, m);
      g_array_append_val (src->std_noise, std_noise);
      g_array_append_val (src->e_rms, e_rms);
      g_array_append_val (src->inv_sigma_crit, mean1);
      g_array_append_val (src->inv_sigma_crit2, mean2);

      g_array_index (count, guint, b)++;
    }
  }

  src->len = src->row->len;

  if (src->len == 0)
    g_error ("nc_galaxy_wl_binned_estimator_build: no galaxy falls in the bins.");

  /* Drop the empty bins */
  g_array_set_size (compact, n_edges - 1);
  src->n_bins = 0;

  for (b = 0; b + 1 < n_edges; b++)
  {
    g_array_index (compact, guint, b) = src->n_bins;

    if (g_array_index (count, guint, b) > 0)
      src->n_bins++;
  }

  for (i = 0; i < src->len; i++)
    NC_WLBE_SRC_BIN (src, i) = g_array_index (compact, guint, NC_WLBE_SRC_BIN (src, i));

  {
    GStrv extra_cols = klass->extra_columns (est);

    binned = nc_galaxy_wl_obs_binned_new (src->n_bins, extra_cols);
    g_strfreev (extra_cols);
  }

  for (b = 0; b + 1 < n_edges; b++)
  {
    const guint n_gal = g_array_index (count, guint, b);

    if (n_gal > 0)
    {
      const guint bc = g_array_index (compact, guint, b);

      nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_OBS_BINNED_COL_R_LO, bc, g_array_index (edges, gdouble, b));
      nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_OBS_BINNED_COL_R_HI, bc, g_array_index (edges, gdouble, b + 1));
      nc_galaxy_wl_obs_binned_set (binned, NC_GALAXY_WL_OBS_BINNED_COL_N_GAL, bc, n_gal);
    }
  }

  nc_galaxy_wl_obs_binned_set_z_cluster (binned, models.z_cl);
  nc_galaxy_wl_obs_binned_set_estimator (binned, klass->name);

  {
    GArray *factor = g_array_new (FALSE, TRUE, sizeof (gdouble));

    g_array_set_size (factor, src->len);
    klass->build_bins (est, src, binned, factor);
    binned = _nc_wlbe_drop_invalid (binned, src, compact);

    {
      const gchar *cols[] = {
        NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_ROW,
        NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_BIN,
        NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_RADIUS,
        NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_FACTOR,
        NULL
      };
      const NcmCatalogColType col_types[] = {
        NCM_CATALOG_COL_TYPE_INT,
        NCM_CATALOG_COL_TYPE_INT,
        NCM_CATALOG_COL_TYPE_DOUBLE,
        NCM_CATALOG_COL_TYPE_DOUBLE,
      };
      guint n_members = 0;
      NcmCatalog *members;

      /* compact now maps the bins of src to the valid bins, G_MAXUINT for a dropped one */
      for (i = 0; i < src->len; i++)
      {
        if ((g_array_index (compact, guint, NC_WLBE_SRC_BIN (src, i)) != G_MAXUINT) && (g_array_index (factor, gdouble, i) != 0.0))
          n_members++;
      }

      if (n_members == 0)
        g_error ("nc_galaxy_wl_binned_estimator_build: no galaxy has a non-zero weight.");

      members   = ncm_catalog_new_full (n_members, (GStrv) cols, col_types, 4);
      n_members = 0;

      for (i = 0; i < src->len; i++)
      {
        const guint bc  = g_array_index (compact, guint, NC_WLBE_SRC_BIN (src, i));
        const gdouble f = g_array_index (factor, gdouble, i);

        if ((bc == G_MAXUINT) || (f == 0.0))
          continue;

        ncm_catalog_set_int (members, NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_ROW, n_members, g_array_index (src->row, guint, i), NULL);
        ncm_catalog_set_int (members, NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_BIN, n_members, bc, NULL);
        ncm_catalog_set (members, NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_RADIUS, n_members, NC_WLBE_SRC (src, radius, i), NULL);
        ncm_catalog_set (members, NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_FACTOR, n_members, f, NULL);
        n_members++;
      }

      nc_galaxy_wl_obs_binned_set_members (binned, obs, members);
      ncm_catalog_free (members);
    }

    g_array_unref (factor);
  }

  _nc_wlbe_models_zf_end (&models);
  _nc_wlbe_sources_free (src);
  g_array_unref (s);
  g_array_unref (w);
  g_array_unref (count);
  g_array_unref (compact);

  return binned;
}

/**
 * nc_galaxy_wl_binned_estimator_eval_model:
 * @est: a #NcGalaxyWLBinnedEstimator
 * @mset: a #NcmMSet
 * @binned: a binned profile built by an estimator of the same type
 * @model: the model, one value per bin, set on return
 *
 * Evaluates the classic model of the estimator for the models of @mset at the bins of
 * @binned.
 */
void
nc_galaxy_wl_binned_estimator_eval_model (NcGalaxyWLBinnedEstimator *est, NcmMSet *mset, NcGalaxyWLObsBinned *binned, NcmVector *model)
{
  g_assert_cmpuint (ncm_vector_len (model), ==, nc_galaxy_wl_obs_binned_len (binned));

  NC_GALAXY_WL_BINNED_ESTIMATOR_GET_CLASS (est)->eval_model (est, mset, binned, model);
}

static void
_nc_wlbe_exact_prepare (NcGalaxyWLBinnedEstimator *est, NcWLBEModels *models, NcGalaxyWLObsBinned *binned)
{
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);
  NcmCatalog *members                           = nc_galaxy_wl_obs_binned_peek_members (binned);
  NcGalaxyWLObs *obs                            = nc_galaxy_wl_obs_binned_peek_source (binned);
  guint64 crit_hash, z_cl_bits, z_hash = 0;

  memcpy (&z_cl_bits, &models->z_cl, sizeof (z_cl_bits));
  crit_hash = ncm_model_state_get_pkey (NCM_MODEL (models->cosmo)) ^
              (ncm_model_state_get_pkey (NCM_MODEL (models->smd)) * 37U) ^
              (z_cl_bits * 43U);

  /* The density of a redshift factor changes with the parameters of its models */
  if (self->redshift_factor != NULL)
  {
    nc_galaxy_redshift_factor_prepare (self->redshift_factor, models->mset);
    z_hash = nc_galaxy_redshift_factor_get_hash (self->redshift_factor);
  }

  if ((members == self->exact_members) &&
      ncm_catalog_key_is_current (members, &self->exact_members_key) &&
      (crit_hash == self->exact_crit_hash) &&
      (self->redshift_factor == self->exact_zf) &&
      (z_hash == self->exact_z_hash) &&
      (self->exact_n_int == self->exact_intervals))
    return;

  {
    const guint n_members = ncm_catalog_len (members);
    guint j;

    ncm_catalog_clear (&self->exact_members);
    self->exact_members     = ncm_catalog_ref (members);
    self->exact_members_key = ncm_catalog_get_key (members);
    self->exact_crit_hash   = crit_hash;
    self->exact_z_hash      = z_hash;
    self->exact_zf          = self->redshift_factor;
    self->exact_n_int       = self->exact_intervals;

    g_array_set_size (self->exact_offset, 0);
    g_array_set_size (self->exact_s, 0);
    g_array_set_size (self->exact_w, 0);

    _nc_wlbe_models_zf_begin (est, models);

    for (j = 0; j < n_members; j++)
    {
      const guint row = (guint) ncm_catalog_get_int (members, NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_ROW, j, NULL);
      const guint off = self->exact_s->len;

      g_array_append_val (self->exact_offset, off);
      _nc_wlbe_append_nodes (models, obs, row, self->exact_intervals, self->exact_s, self->exact_w);
    }

    {
      const guint off = self->exact_s->len;

      g_array_append_val (self->exact_offset, off);
    }

    _nc_wlbe_models_zf_end (models);
  }
}

/**
 * nc_galaxy_wl_binned_estimator_eval_model_exact:
 * @est: a #NcGalaxyWLBinnedEstimator
 * @mset: a #NcmMSet
 * @binned: a binned profile built by @est from a per-galaxy catalog
 * @model: the model, one value per bin, set on return
 *
 * Evaluates the expected value of the estimator for the models of @mset,
 * $\langle\hat d_b\rangle = o_b + \sum_{i \in b} f_i\,\bar g_i$ (see #NcGalaxyWLObsBinned),
 * with $\bar g_i = \int P_i(z)\, g(R_i, z)\,\mathrm{d}z / \int P_i(z)\,\mathrm{d}z$ the
 * reduced shear $g = \Sigma_\mathrm{cr}^{-1}\Delta\Sigma / (1 - \Sigma_\mathrm{cr}^{-1}\Sigma)$
 * of galaxy $i$ at its radius $R_i$, averaged over its redshift density. A profile without
 * members aborts.
 */
void
nc_galaxy_wl_binned_estimator_eval_model_exact (NcGalaxyWLBinnedEstimator *est, NcmMSet *mset, NcGalaxyWLObsBinned *binned, NcmVector *model)
{
  NcGalaxyWLBinnedEstimatorPrivate * const self = nc_galaxy_wl_binned_estimator_get_instance_private (est);
  const guint n_bins                            = nc_galaxy_wl_obs_binned_len (binned);
  NcmCatalog *members;
  NcWLBEModels models;
  guint j, b;

  g_assert_cmpuint (ncm_vector_len (model), ==, n_bins);

  if (!nc_galaxy_wl_obs_binned_has_members (binned))
    g_error ("nc_galaxy_wl_binned_estimator_eval_model_exact: the profile has no members "
             "(built from external arrays); only the classic model is available.");

  _nc_wlbe_models (mset, &models);

  if (models.dp == NULL)
    g_error ("nc_galaxy_wl_binned_estimator_eval_model_exact: the mset must hold a NcHaloDensityProfile.");

  _nc_wlbe_exact_prepare (est, &models, binned);
  members = nc_galaxy_wl_obs_binned_peek_members (binned);

  for (b = 0; b < n_bins; b++)
    ncm_vector_set (model, b, nc_galaxy_wl_obs_binned_get (binned, NC_GALAXY_WL_OBS_BINNED_COL_OFFSET, b));

  for (j = 0; j < ncm_catalog_len (members); j++)
  {
    const guint bin      = (guint) ncm_catalog_get_int (members, NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_BIN, j, NULL);
    const gdouble radius = ncm_catalog_get (members, NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_RADIUS, j, NULL);
    const gdouble factor = ncm_catalog_get (members, NC_GALAXY_WL_OBS_BINNED_MEMBER_COL_FACTOR, j, NULL);
    const guint k0       = g_array_index (self->exact_offset, guint, j);
    const guint k1       = g_array_index (self->exact_offset, guint, j + 1);
    gdouble g_mean       = 0.0;
    guint k;

    if ((factor == 0.0) || (k0 == k1))
      continue;

    {
      const gdouble sigma   = nc_wl_surface_mass_density_sigma (models.smd, models.dp, models.cosmo, radius, models.z_cl);
      const gdouble d_sigma = nc_wl_surface_mass_density_sigma_excess (models.smd, models.dp, models.cosmo, radius, models.z_cl);

      for (k = k0; k < k1; k++)
      {
        const gdouble s_k = g_array_index (self->exact_s, gdouble, k);
        const gdouble w_k = g_array_index (self->exact_w, gdouble, k);

        g_mean += w_k * s_k * d_sigma / (1.0 - s_k * sigma);
      }
    }

    ncm_vector_addto (model, bin, factor * g_mean);
  }
}

