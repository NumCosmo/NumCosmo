/***************************************************************************
 *            nc_galaxy_redshift_obs.c
 *
 *  Tue Jul 1 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_obs.c
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyRedshiftObs:
 *
 * Abstract model for the photometric-redshift observable model $P(\mathrm{data}|z)$.
 *
 * This small #NcmModel describes the conditional density of a per-galaxy
 * photometric-redshift observation given the true redshift $z$. The whole
 * observation (e.g. the point estimate $z_\mathrm{phot}$ and its scatter
 * $\sigma_0$) is carried per-galaxy in a #NcGalaxyRedshiftObsData; the
 * observation's structure is defined by the observable model itself, so the redshift
 * calculator that convolves this kernel with the true-redshift distribution
 * $P(z|I)$ stays observation-agnostic and only supplies the integration
 * variable $z$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_redshift_obs.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_ABSTRACT_TYPE (NcGalaxyRedshiftObs, nc_galaxy_redshift_obs, NCM_TYPE_MODEL);
G_DEFINE_BOXED_TYPE (NcGalaxyRedshiftObsData, nc_galaxy_redshift_obs_data, nc_galaxy_redshift_obs_data_ref, nc_galaxy_redshift_obs_data_unref); /* LCOV_EXCL_LINE */

static void
nc_galaxy_redshift_obs_init (NcGalaxyRedshiftObs *gsdre)
{
}

static void
_nc_galaxy_redshift_obs_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_redshift_obs_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_redshift_obs, NC_TYPE_GALAXY_REDSHIFT_OBS);

/* LCOV_EXCL_START */
static void
_nc_galaxy_redshift_obs_data_init (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data)
{
  g_error ("_nc_galaxy_redshift_obs_data_init: method not implemented.");
}

static gdouble
_nc_galaxy_redshift_obs_eval (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z)
{
  g_error ("_nc_galaxy_redshift_obs_eval: method not implemented.");

  return 0.0;
}

static gdouble
_nc_galaxy_redshift_obs_gen (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z, NcmRNG *rng)
{
  g_error ("_nc_galaxy_redshift_obs_gen: method not implemented.");

  return 0.0;
}

static gdouble
_nc_galaxy_redshift_obs_window_mass (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z, const gdouble obs_lo, const gdouble obs_hi)
{
  g_error ("_nc_galaxy_redshift_obs_window_mass: method not implemented.");

  return 0.0;
}

static void
_nc_galaxy_redshift_obs_get_true_z_lim (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, gdouble *z_min, gdouble *z_max)
{
  g_error ("_nc_galaxy_redshift_obs_get_true_z_lim: method not implemented.");
}

/* LCOV_EXCL_STOP */

static void
nc_galaxy_redshift_obs_class_init (NcGalaxyRedshiftObsClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_redshift_obs_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy photometric-redshift observable model", "GalaxyRedshiftObs");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);
  ncm_mset_model_register_id (model_class, "NcGalaxyRedshiftObs", "Galaxy photometric-redshift observable model", NULL, FALSE, NCM_MSET_MODEL_MAIN);
  ncm_model_class_check_params_info (model_class);

  klass->data_init      = &_nc_galaxy_redshift_obs_data_init;
  klass->eval           = &_nc_galaxy_redshift_obs_eval;
  klass->gen            = &_nc_galaxy_redshift_obs_gen;
  klass->window_mass    = &_nc_galaxy_redshift_obs_window_mass;
  klass->get_true_z_lim = &_nc_galaxy_redshift_obs_get_true_z_lim;
}

/**
 * nc_galaxy_redshift_obs_ref:
 * @gsdre: a #NcGalaxyRedshiftObs
 *
 * Increases the reference count of @gsdre by one.
 *
 * Returns: (transfer full): @gsdre.
 */
NcGalaxyRedshiftObs *
nc_galaxy_redshift_obs_ref (NcGalaxyRedshiftObs *gsdre)
{
  return g_object_ref (gsdre);
}

/**
 * nc_galaxy_redshift_obs_free:
 * @gsdre: a #NcGalaxyRedshiftObs
 *
 * Decreases the reference count of @gsdre by one.
 *
 */
void
nc_galaxy_redshift_obs_free (NcGalaxyRedshiftObs *gsdre)
{
  g_object_unref (gsdre);
}

/**
 * nc_galaxy_redshift_obs_clear:
 * @gsdre: a #NcGalaxyRedshiftObs
 *
 * Decreases the reference count of *@gsdre by one, and sets the pointer *@gsdre
 * to NULL.
 *
 */
void
nc_galaxy_redshift_obs_clear (NcGalaxyRedshiftObs **gsdre)
{
  g_clear_object (gsdre);
}

/**
 * nc_galaxy_redshift_obs_data_new:
 * @gsdre: a #NcGalaxyRedshiftObs
 *
 * Creates a new per-galaxy #NcGalaxyRedshiftObsData for @gsdre.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftObsData.
 */
NcGalaxyRedshiftObsData *
nc_galaxy_redshift_obs_data_new (NcGalaxyRedshiftObs *gsdre)
{
  NcGalaxyRedshiftObsData *data = g_new0 (NcGalaxyRedshiftObsData, 1);

  data->ldata                  = NULL;
  data->ldata_destroy          = NULL;
  data->ldata_read_row         = NULL;
  data->ldata_write_row        = NULL;
  data->ldata_required_columns = NULL;

  g_atomic_ref_count_init (&data->ref_count);
  NC_GALAXY_REDSHIFT_OBS_GET_CLASS (gsdre)->data_init (gsdre, data);

  g_assert_nonnull (data->ldata_destroy);
  g_assert_nonnull (data->ldata_read_row);
  g_assert_nonnull (data->ldata_write_row);
  g_assert_nonnull (data->ldata_required_columns);

  return data;
}

/**
 * nc_galaxy_redshift_obs_data_ref:
 * @data: a #NcGalaxyRedshiftObsData
 *
 * Increases the reference count of @data by one.
 *
 * Returns: (transfer full): @data.
 */
NcGalaxyRedshiftObsData *
nc_galaxy_redshift_obs_data_ref (NcGalaxyRedshiftObsData *data)
{
  g_atomic_ref_count_inc (&data->ref_count);

  return data;
}

/**
 * nc_galaxy_redshift_obs_data_unref:
 * @data: a #NcGalaxyRedshiftObsData
 *
 * Decreases the reference count of @data by one. If it reaches 0, @data is freed.
 *
 */
void
nc_galaxy_redshift_obs_data_unref (NcGalaxyRedshiftObsData *data)
{
  if (g_atomic_ref_count_dec (&data->ref_count))
  {
    g_assert_nonnull (data->ldata_destroy);
    data->ldata_destroy (data->ldata);
    g_free (data);
  }
}

/**
 * nc_galaxy_redshift_obs_data_read_row:
 * @data: a #NcGalaxyRedshiftObsData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Reads the photometric-redshift observation of row @i.
 *
 */
void
nc_galaxy_redshift_obs_data_read_row (NcGalaxyRedshiftObsData *data, NcGalaxyWLObs *obs, const guint i)
{
  data->ldata_read_row (data, obs, i);
}

/**
 * nc_galaxy_redshift_obs_data_write_row:
 * @data: a #NcGalaxyRedshiftObsData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Writes the photometric-redshift observation of row @i.
 *
 */
void
nc_galaxy_redshift_obs_data_write_row (NcGalaxyRedshiftObsData *data, NcGalaxyWLObs *obs, const guint i)
{
  data->ldata_write_row (data, obs, i);
}

/**
 * nc_galaxy_redshift_obs_data_required_columns:
 * @data: a #NcGalaxyRedshiftObsData
 *
 * Returns: (element-type utf8) (transfer full): the required columns.
 */
GList *
nc_galaxy_redshift_obs_data_required_columns (NcGalaxyRedshiftObsData *data)
{
  GList *columns = NULL;

  data->ldata_required_columns (data, &columns);

  return columns;
}

/**
 * nc_galaxy_redshift_obs_eval: (virtual eval)
 * @gsdre: a #NcGalaxyRedshiftObs
 * @data: a #NcGalaxyRedshiftObsData carrying the per-galaxy observation
 * @z: the true redshift $z$
 *
 * Evaluates the conditional density $P(\mathrm{data}|z)$ of the per-galaxy
 * observation in @data given the true redshift @z.
 *
 * Returns: the conditional density $P(\mathrm{data}|z)$.
 */
gdouble
nc_galaxy_redshift_obs_eval (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z)
{
  return NC_GALAXY_REDSHIFT_OBS_GET_CLASS (gsdre)->eval (gsdre, data, z);
}

/**
 * nc_galaxy_redshift_obs_gen: (virtual gen)
 * @gsdre: a #NcGalaxyRedshiftObs
 * @data: a #NcGalaxyRedshiftObsData whose observation inputs are set
 * @z: the true redshift $z$
 * @rng: a #NcmRNG
 *
 * Samples a photometric-redshift observation given the true redshift @z, storing
 * the sampled observable(s) into @data and returning the sampled point estimate.
 *
 * Returns: the sampled photometric redshift.
 */
gdouble
nc_galaxy_redshift_obs_gen (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z, NcmRNG *rng)
{
  return NC_GALAXY_REDSHIFT_OBS_GET_CLASS (gsdre)->gen (gsdre, data, z, rng);
}

/**
 * nc_galaxy_redshift_obs_window_mass: (virtual window_mass)
 * @gsdre: a #NcGalaxyRedshiftObs
 * @data: a #NcGalaxyRedshiftObsData whose observation inputs are set
 * @z: the true redshift $z$
 * @obs_lo: the lower edge of the observable window
 * @obs_hi: the upper edge of the observable window
 *
 * Computes the probability mass of the observable conditional in the window
 * $[\mathtt{obs\_lo}, \mathtt{obs\_hi}]$ at fixed true redshift @z, i.e.
 * $\int_{\mathtt{obs\_lo}}^{\mathtt{obs\_hi}} P(\mathrm{obs} \mid z)\,
 * \mathrm{d}\,\mathrm{obs}$. This is the selection/normalization factor used by
 * the redshift calculator to condition on a photometric selection window and by
 * the binned $\mathrm{d}n/\mathrm{d}z$; it depends only on @z and the per-galaxy
 * scatter, not on the sampled point estimate.
 *
 * Returns: the probability mass in the window.
 */
gdouble
nc_galaxy_redshift_obs_window_mass (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, const gdouble z, const gdouble obs_lo, const gdouble obs_hi)
{
  return NC_GALAXY_REDSHIFT_OBS_GET_CLASS (gsdre)->window_mass (gsdre, data, z, obs_lo, obs_hi);
}

/**
 * nc_galaxy_redshift_obs_get_true_z_lim: (virtual get_true_z_lim)
 * @gsdre: a #NcGalaxyRedshiftObs
 * @data: a #NcGalaxyRedshiftObsData whose observation inputs are set
 * @z_min: (out): the lower edge of the effective true-redshift support
 * @z_max: (out): the upper edge of the effective true-redshift support
 *
 * Returns the range of true redshift $z$ over which this galaxy's observable
 * conditional $P(\mathrm{obs} \mid z)$ is non-negligible (e.g. the Gaussian
 * kernel restricted to a few sigma around the point estimate). The redshift
 * calculator intersects this with the population support to build the effective
 * $z$-integration limits, keeping the adaptive quadrature away from resolving a
 * needle in a haystack. Unclamped: the caller applies the population bounds.
 *
 */
void
nc_galaxy_redshift_obs_get_true_z_lim (NcGalaxyRedshiftObs *gsdre, NcGalaxyRedshiftObsData *data, gdouble *z_min, gdouble *z_max)
{
  NC_GALAXY_REDSHIFT_OBS_GET_CLASS (gsdre)->get_true_z_lim (gsdre, data, z_min, z_max);
}

