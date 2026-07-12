/***************************************************************************
 *            nc_galaxy_redshift_obs_sel.c
 *
 *  Tue Jul 1 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_obs_sel.c
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
 * NcGalaxyRedshiftObsSel:
 *
 * Abstract population-level photo-z observable distribution.
 *
 * Models the distribution of a photo-z observable across the galaxy population at
 * a given true redshift, `P(obs | z)` at the population level. Distinct from the
 * per-galaxy conditional #NcGalaxyRedshiftObs (whose scatter is
 * per-galaxy data): here the scatter is a population-level model parameter, and
 * the distribution is conceptually a mixture over the per-galaxy scatter. It
 * supplies the selection mass `window_mass(z, obs_lo, obs_hi)` that the binning
 * calculator uses to build the binned true-redshift distribution `P(z | I in W)`.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_redshift_obs_sel.h"

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_ABSTRACT_TYPE (NcGalaxyRedshiftObsSel, nc_galaxy_redshift_obs_sel, NCM_TYPE_MODEL);

static void
nc_galaxy_redshift_obs_sel_init (NcGalaxyRedshiftObsSel *gsdrop)
{
}

static void
_nc_galaxy_redshift_obs_sel_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_redshift_obs_sel_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_redshift_obs_sel, NC_TYPE_GALAXY_REDSHIFT_OBS_SEL);

/* LCOV_EXCL_START */
static gdouble
_nc_galaxy_redshift_obs_sel_eval (NcGalaxyRedshiftObsSel *gsdrop, const gdouble z, const gdouble obs)
{
  g_error ("_nc_galaxy_redshift_obs_sel_eval: method not implemented.");

  return 0.0;
}

static gdouble
_nc_galaxy_redshift_obs_sel_window_mass (NcGalaxyRedshiftObsSel *gsdrop, const gdouble z, const gdouble obs_lo, const gdouble obs_hi)
{
  g_error ("_nc_galaxy_redshift_obs_sel_window_mass: method not implemented.");

  return 0.0;
}

/* LCOV_EXCL_STOP */

static void
nc_galaxy_redshift_obs_sel_class_init (NcGalaxyRedshiftObsSelClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_redshift_obs_sel_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy population photo-z observable distribution", "GalaxyRedshiftObsSel");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);
  ncm_mset_model_register_id (model_class, "NcGalaxyRedshiftObsSel", "Galaxy population photo-z observable distribution", NULL, FALSE, NCM_MSET_MODEL_MAIN);
  ncm_model_class_check_params_info (model_class);

  klass->eval        = &_nc_galaxy_redshift_obs_sel_eval;
  klass->window_mass = &_nc_galaxy_redshift_obs_sel_window_mass;
}

/**
 * nc_galaxy_redshift_obs_sel_ref:
 * @gsdrop: a #NcGalaxyRedshiftObsSel
 *
 * Increases the reference count of @gsdrop by one.
 *
 * Returns: (transfer full): @gsdrop
 */
NcGalaxyRedshiftObsSel *
nc_galaxy_redshift_obs_sel_ref (NcGalaxyRedshiftObsSel *gsdrop)
{
  return g_object_ref (gsdrop);
}

/**
 * nc_galaxy_redshift_obs_sel_free:
 * @gsdrop: a #NcGalaxyRedshiftObsSel
 *
 * Decreases the reference count of @gsdrop by one.
 *
 */
void
nc_galaxy_redshift_obs_sel_free (NcGalaxyRedshiftObsSel *gsdrop)
{
  g_object_unref (gsdrop);
}

/**
 * nc_galaxy_redshift_obs_sel_clear:
 * @gsdrop: a #NcGalaxyRedshiftObsSel
 *
 * Decreases the reference count of @gsdrop by one, and sets the pointer *@gsdrop
 * to NULL.
 *
 */
void
nc_galaxy_redshift_obs_sel_clear (NcGalaxyRedshiftObsSel **gsdrop)
{
  g_clear_object (gsdrop);
}

/**
 * nc_galaxy_redshift_obs_sel_eval: (virtual eval)
 * @gsdrop: a #NcGalaxyRedshiftObsSel
 * @z: the true redshift $z$
 * @obs: the observable value
 *
 * Computes the population-level probability density of the observable at value
 * @obs and fixed true redshift @z, i.e. the density whose integral over a window
 * is nc_galaxy_redshift_obs_sel_window_mass(). Used by the binning
 * calculator to build the marginal observable distribution $P(z_p)$.
 *
 * Returns: the population observable density at @obs.
 */
gdouble
nc_galaxy_redshift_obs_sel_eval (NcGalaxyRedshiftObsSel *gsdrop, const gdouble z, const gdouble obs)
{
  return NC_GALAXY_REDSHIFT_OBS_SEL_GET_CLASS (gsdrop)->eval (gsdrop, z, obs);
}

/**
 * nc_galaxy_redshift_obs_sel_window_mass: (virtual window_mass)
 * @gsdrop: a #NcGalaxyRedshiftObsSel
 * @z: the true redshift $z$
 * @obs_lo: the lower edge of the observable window
 * @obs_hi: the upper edge of the observable window
 *
 * Computes the population-level probability mass of the observable in the window
 * $[\mathtt{obs\_lo}, \mathtt{obs\_hi}]$ at fixed true redshift @z, i.e. the
 * fraction of the population at @z whose observable falls in the window. This is
 * the selection/normalization factor used by the binning calculator.
 *
 * Returns: the population probability mass in the window.
 */
gdouble
nc_galaxy_redshift_obs_sel_window_mass (NcGalaxyRedshiftObsSel *gsdrop, const gdouble z, const gdouble obs_lo, const gdouble obs_hi)
{
  return NC_GALAXY_REDSHIFT_OBS_SEL_GET_CLASS (gsdrop)->window_mass (gsdrop, z, obs_lo, obs_hi);
}

