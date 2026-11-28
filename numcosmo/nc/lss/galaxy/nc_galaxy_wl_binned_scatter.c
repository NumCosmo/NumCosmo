/***************************************************************************
 *            nc_galaxy_wl_binned_scatter.c
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_binned_scatter.c
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
 * NcGalaxyWLBinnedScatter:
 *
 * Abstract model of the scatter of a binned weak-lensing profile around its model.
 *
 * The distribution of the measured binned profile (#NcGalaxyWLObsBinned) given the model
 * of its estimator (nc_galaxy_wl_binned_estimator_eval_model()). It is a #NcmModel of the
 * #NcmMSet, so a concrete scatter may have parameters to fit, as an intrinsic scatter of
 * the profile. #NcGalaxyWLBinnedScatterGauss is the Gaussian of the covariance of the
 * profile.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_wl_binned_scatter.h"

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_ABSTRACT_TYPE (NcGalaxyWLBinnedScatter, nc_galaxy_wl_binned_scatter, NCM_TYPE_MODEL)

static void
nc_galaxy_wl_binned_scatter_init (NcGalaxyWLBinnedScatter *scatter)
{
}

static void
_nc_galaxy_wl_binned_scatter_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_wl_binned_scatter_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_wl_binned_scatter, NC_TYPE_GALAXY_WL_BINNED_SCATTER);

/* LCOV_EXCL_START */
static gdouble
_nc_galaxy_wl_binned_scatter_m2lnL (NcGalaxyWLBinnedScatter *scatter, NcGalaxyWLObsBinned *binned, NcmVector *model)
{
  g_error ("_nc_galaxy_wl_binned_scatter_m2lnL: method not implemented.");

  return 0.0;
}

static void
_nc_galaxy_wl_binned_scatter_gen (NcGalaxyWLBinnedScatter *scatter, NcGalaxyWLObsBinned *binned, NcmVector *model, NcmRNG *rng, NcmVector *data)
{
  g_error ("_nc_galaxy_wl_binned_scatter_gen: method not implemented.");
}

/* LCOV_EXCL_STOP */

static void
nc_galaxy_wl_binned_scatter_class_init (NcGalaxyWLBinnedScatterClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_wl_binned_scatter_finalize;

  ncm_model_class_set_name_nick (model_class, "Binned weak-lensing profile scatter", "WLBinnedScatter");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);
  ncm_mset_model_register_id (model_class, "NcGalaxyWLBinnedScatter", "Binned weak-lensing profile scatter", NULL, FALSE, NCM_MSET_MODEL_MAIN);
  ncm_model_class_check_params_info (model_class);

  klass->m2lnL = &_nc_galaxy_wl_binned_scatter_m2lnL;
  klass->gen   = &_nc_galaxy_wl_binned_scatter_gen;
}

/**
 * nc_galaxy_wl_binned_scatter_ref:
 * @scatter: a #NcGalaxyWLBinnedScatter
 *
 * Increases the reference count of @scatter by one.
 *
 * Returns: (transfer full): @scatter.
 */
NcGalaxyWLBinnedScatter *
nc_galaxy_wl_binned_scatter_ref (NcGalaxyWLBinnedScatter *scatter)
{
  return g_object_ref (scatter);
}

/**
 * nc_galaxy_wl_binned_scatter_free:
 * @scatter: a #NcGalaxyWLBinnedScatter
 *
 * Decreases the reference count of @scatter by one.
 */
void
nc_galaxy_wl_binned_scatter_free (NcGalaxyWLBinnedScatter *scatter)
{
  g_object_unref (scatter);
}

/**
 * nc_galaxy_wl_binned_scatter_clear:
 * @scatter: a #NcGalaxyWLBinnedScatter
 *
 * Decreases the reference count of *@scatter by one and sets *@scatter to %NULL.
 */
void
nc_galaxy_wl_binned_scatter_clear (NcGalaxyWLBinnedScatter **scatter)
{
  g_clear_object (scatter);
}

/**
 * nc_galaxy_wl_binned_scatter_m2lnL: (virtual m2lnL)
 * @scatter: a #NcGalaxyWLBinnedScatter
 * @binned: a #NcGalaxyWLObsBinned
 * @model: the model of the profile, one value per bin
 *
 * Computes $-2\ln L$ of the measured profile of @binned given @model.
 *
 * Returns: $-2\ln L$.
 */
gdouble
nc_galaxy_wl_binned_scatter_m2lnL (NcGalaxyWLBinnedScatter *scatter, NcGalaxyWLObsBinned *binned, NcmVector *model)
{
  return NC_GALAXY_WL_BINNED_SCATTER_GET_CLASS (scatter)->m2lnL (scatter, binned, model);
}

/**
 * nc_galaxy_wl_binned_scatter_gen: (virtual gen)
 * @scatter: a #NcGalaxyWLBinnedScatter
 * @binned: a #NcGalaxyWLObsBinned
 * @model: the model of the profile, one value per bin
 * @rng: a #NcmRNG
 * @data: the generated profile, one value per bin, set on return
 *
 * Draws a profile from the scatter around @model.
 */
void
nc_galaxy_wl_binned_scatter_gen (NcGalaxyWLBinnedScatter *scatter, NcGalaxyWLObsBinned *binned, NcmVector *model, NcmRNG *rng, NcmVector *data)
{
  NC_GALAXY_WL_BINNED_SCATTER_GET_CLASS (scatter)->gen (scatter, binned, model, rng, data);
}

