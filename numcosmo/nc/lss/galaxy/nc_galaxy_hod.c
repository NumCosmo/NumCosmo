/***************************************************************************
 *            nc_galaxy_hod.c
 *
 *  Sun Jun 14 12:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_hod.c
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
 * NcGalaxyHOD:
 *
 * Abstract halo occupation distribution.
 *
 * Base class for halo occupation distribution (HOD) models, describing how many
 * central and satellite galaxies populate a halo of a given mass. Concrete
 * subclasses provide the mean central and satellite occupations as functions of
 * the halo mass (nc_galaxy_hod_mean_n_central() and
 * nc_galaxy_hod_mean_n_satellite()); this base class turns those means into an
 * integer realization in nc_galaxy_hod_gen().
 *
 * The central galaxy is drawn as a Bernoulli trial with probability equal to the
 * mean central occupation (clamped to [0, 1]) when
 * #NcGalaxyHOD:stochastic-central is %TRUE, or placed deterministically when the
 * mean is at least one half otherwise. Satellites are drawn from a Poisson
 * distribution and only when a central is present.
 *
 * As an #NcmModel its occupation parameters live in the model parameter vector
 * and can be fit; it is meant to be composed by the mock member generators.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include "nc/lss/galaxy/nc_galaxy_hod.h"
#include "math/ncm_rng.h"

typedef struct _NcGalaxyHODPrivate
{
  gboolean stochastic_central;
} NcGalaxyHODPrivate;

enum
{
  PROP_0,
  PROP_STOCHASTIC_CENTRAL,
  PROP_LEN,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcGalaxyHOD, nc_galaxy_hod, NCM_TYPE_MODEL)

static void
nc_galaxy_hod_init (NcGalaxyHOD *hod)
{
  NcGalaxyHODPrivate * const self = nc_galaxy_hod_get_instance_private (hod);

  self->stochastic_central = TRUE;
}

static void
_nc_galaxy_hod_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyHOD *hod = NC_GALAXY_HOD (object);

  g_return_if_fail (NC_IS_GALAXY_HOD (hod));

  switch (prop_id)
  {
    case PROP_STOCHASTIC_CENTRAL:
      nc_galaxy_hod_set_stochastic_central (hod, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_hod_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyHOD *hod = NC_GALAXY_HOD (object);

  g_return_if_fail (NC_IS_GALAXY_HOD (hod));

  switch (prop_id)
  {
    case PROP_STOCHASTIC_CENTRAL:
      g_value_set_boolean (value, nc_galaxy_hod_get_stochastic_central (hod));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_hod_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_hod_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_hod, NC_TYPE_GALAXY_HOD);

/* LCOV_EXCL_START */
static gdouble
_nc_galaxy_hod_mean_n_central (NcGalaxyHOD *hod, const gdouble lnM)
{
  g_error ("_nc_galaxy_hod_mean_n_central: not implemented by `%s'.", G_OBJECT_TYPE_NAME (hod));

  return 0.0;
}

static gdouble
_nc_galaxy_hod_mean_n_satellite (NcGalaxyHOD *hod, const gdouble lnM)
{
  g_error ("_nc_galaxy_hod_mean_n_satellite: not implemented by `%s'.", G_OBJECT_TYPE_NAME (hod));

  return 0.0;
}

/* LCOV_EXCL_STOP */

static void
nc_galaxy_hod_class_init (NcGalaxyHODClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_hod_set_property;
  model_class->get_property = &_nc_galaxy_hod_get_property;
  object_class->finalize    = &_nc_galaxy_hod_finalize;

  ncm_model_class_set_name_nick (model_class, "Halo occupation distribution", "GalaxyHOD");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  /**
   * NcGalaxyHOD:stochastic-central:
   *
   * Whether the central galaxy is a Bernoulli draw with probability equal to the
   * mean central occupation (%TRUE) or placed deterministically when the mean is
   * at least one half (%FALSE).
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_STOCHASTIC_CENTRAL,
                                   g_param_spec_boolean ("stochastic-central",
                                                         "Stochastic central",
                                                         "Whether the central occupation is a Bernoulli draw",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  ncm_mset_model_register_id (model_class, "NcGalaxyHOD", "Halo occupation distribution", NULL, TRUE, NCM_MSET_MODEL_MAIN);
  ncm_model_class_check_params_info (model_class);

  klass->mean_n_central   = &_nc_galaxy_hod_mean_n_central;
  klass->mean_n_satellite = &_nc_galaxy_hod_mean_n_satellite;
}

/**
 * nc_galaxy_hod_ref:
 * @hod: a #NcGalaxyHOD
 *
 * Increases the reference count of @hod by one.
 *
 * Returns: (transfer full): @hod.
 */
NcGalaxyHOD *
nc_galaxy_hod_ref (NcGalaxyHOD *hod)
{
  return g_object_ref (hod);
}

/**
 * nc_galaxy_hod_free:
 * @hod: a #NcGalaxyHOD
 *
 * Decreases the reference count of @hod by one.
 *
 */
void
nc_galaxy_hod_free (NcGalaxyHOD *hod)
{
  g_object_unref (hod);
}

/**
 * nc_galaxy_hod_clear:
 * @hod: a #NcGalaxyHOD
 *
 * If *@hod is different from %NULL, decreases its reference count and sets
 * *@hod to %NULL.
 *
 */
void
nc_galaxy_hod_clear (NcGalaxyHOD **hod)
{
  g_clear_object (hod);
}

/**
 * nc_galaxy_hod_set_stochastic_central:
 * @hod: a #NcGalaxyHOD
 * @stochastic_central: whether the central occupation is a Bernoulli draw
 *
 * Sets whether the central galaxy is drawn stochastically.
 *
 */
void
nc_galaxy_hod_set_stochastic_central (NcGalaxyHOD *hod, gboolean stochastic_central)
{
  NcGalaxyHODPrivate * const self = nc_galaxy_hod_get_instance_private (hod);

  self->stochastic_central = stochastic_central;
}

/**
 * nc_galaxy_hod_get_stochastic_central:
 * @hod: a #NcGalaxyHOD
 *
 * Gets whether the central galaxy is drawn stochastically.
 *
 * Returns: %TRUE if the central occupation is a Bernoulli draw, %FALSE otherwise.
 */
gboolean
nc_galaxy_hod_get_stochastic_central (NcGalaxyHOD *hod)
{
  NcGalaxyHODPrivate * const self = nc_galaxy_hod_get_instance_private (hod);

  return self->stochastic_central;
}

/**
 * nc_galaxy_hod_mean_n_central:
 * @hod: a #NcGalaxyHOD
 * @lnM: the natural logarithm of the halo mass (solar masses)
 *
 * Computes the mean central galaxy occupation of a halo of mass $e^{\ln M}$.
 *
 * Returns: the mean central occupation $\langle N_\mathrm{cen} \rangle$.
 */
gdouble
nc_galaxy_hod_mean_n_central (NcGalaxyHOD *hod, const gdouble lnM)
{
  return NC_GALAXY_HOD_GET_CLASS (hod)->mean_n_central (hod, lnM);
}

/**
 * nc_galaxy_hod_mean_n_satellite:
 * @hod: a #NcGalaxyHOD
 * @lnM: the natural logarithm of the halo mass (solar masses)
 *
 * Computes the mean satellite galaxy occupation of a halo of mass $e^{\ln M}$.
 *
 * Returns: the mean satellite occupation $\langle N_\mathrm{sat} \rangle$.
 */
gdouble
nc_galaxy_hod_mean_n_satellite (NcGalaxyHOD *hod, const gdouble lnM)
{
  return NC_GALAXY_HOD_GET_CLASS (hod)->mean_n_satellite (hod, lnM);
}

/**
 * nc_galaxy_hod_gen:
 * @hod: a #NcGalaxyHOD
 * @lnM: the natural logarithm of the halo mass (solar masses)
 * @rng: a #NcmRNG
 * @n_central: (out): the realized number of central galaxies (0 or 1)
 * @n_satellite: (out): the realized number of satellite galaxies
 *
 * Draws an integer realization of the central and satellite occupations of a
 * halo of mass $e^{\ln M}$. Satellites are drawn only when a central is present.
 *
 */
void
nc_galaxy_hod_gen (NcGalaxyHOD *hod, const gdouble lnM, NcmRNG *rng, gint *n_central, gint *n_satellite)
{
  NcGalaxyHODPrivate * const self = nc_galaxy_hod_get_instance_private (hod);
  const gdouble mean_n_central    = nc_galaxy_hod_mean_n_central (hod, lnM);
  gint n_cen;
  gint n_sat = 0;

  ncm_rng_lock (rng);

  if (self->stochastic_central)
  {
    const gdouble p = CLAMP (mean_n_central, 0.0, 1.0);

    n_cen = (ncm_rng_uniform01_gen (rng) < p) ? 1 : 0;
  }
  else
  {
    n_cen = (mean_n_central >= 0.5) ? 1 : 0;
  }

  if (n_cen >= 1)
  {
    const gdouble mean_n_satellite = nc_galaxy_hod_mean_n_satellite (hod, lnM);

    if (mean_n_satellite > 0.0)
      n_sat = ncm_rng_poisson_gen (rng, mean_n_satellite);
  }

  ncm_rng_unlock (rng);

  *n_central   = n_cen;
  *n_satellite = n_sat;
}
