/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_true_redshift.c
 *
 *  Wed Jul 31 20:52:43 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_true_redshift.c
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * SECTION: nc_galaxy_sd_true_redshift
 * @title: NcGalaxySDTrueRedshift
 * @short_description: Class describing galaxy sample redshift distributions.
 * @stability: Unstable
 *
 *
 * This class describes a galaxy sample redshift distributions.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "galaxy/nc_galaxy_sd_true_redshift.h"
#include "galaxy/nc_galaxy_sd_obs_redshift.h"
#include "math/ncm_dtuple.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_stats_dist1d_spline.h"

typedef struct _NcGalaxySDTrueRedshiftPrivate
{
  gint placeholder;
} NcGalaxySDTrueRedshiftPrivate;

enum
{
  PROP_0,
  PROP_LIM,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDTrueRedshift, nc_galaxy_sd_true_redshift, NCM_TYPE_MODEL)

static void
nc_galaxy_sd_true_redshift_init (NcGalaxySDTrueRedshift *gsdtr)
{
}

static void
_nc_galaxy_sd_true_redshift_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDTrueRedshift *gsdtr = NC_GALAXY_SD_TRUE_REDSHIFT (object);

  g_return_if_fail (NC_IS_GALAXY_SD_TRUE_REDSHIFT (gsdtr));

  switch (prop_id)
  {
    case PROP_LIM:
    {
      NcmDTuple2 *lim = g_value_get_boxed (value);

      if (lim == NULL)
        g_error ("_nc_galaxy_sd_true_redshift_set_property: lim is NULL");

      nc_galaxy_sd_true_redshift_set_lim (gsdtr, lim->elements[0], lim->elements[1]);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_true_redshift_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDTrueRedshift *gsdtr = NC_GALAXY_SD_TRUE_REDSHIFT (object);

  g_return_if_fail (NC_IS_GALAXY_SD_TRUE_REDSHIFT (gsdtr));

  switch (prop_id)
  {
    case PROP_LIM:
    {
      gdouble z_min, z_max;

      nc_galaxy_sd_true_redshift_get_lim (gsdtr, &z_min, &z_max);

      g_value_take_boxed (value, ncm_dtuple2_new (z_min, z_max));
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_sd_true_redshift_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_true_redshift_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_sd_true_redshift, NC_TYPE_GALAXY_SD_TRUE_REDSHIFT);

/* LCOV_EXCL_START */
static gdouble
_nc_galaxy_sd_true_redshift_gen (NcGalaxySDTrueRedshift *gsdtr, NcmRNG *rng)
{
  g_error ("_nc_galaxy_sd_true_redshift_gen_z: not implemented");

  return 0.0;
}

static gdouble
_nc_galaxy_sd_true_redshift_integ (NcGalaxySDTrueRedshift *gsdtr, gdouble z)
{
  g_error ("_nc_galaxy_sd_true_redshift_integ: not implemented");

  return 0.0;
}

static gboolean
_nc_galaxy_sd_true_redshift_set_lim (NcGalaxySDTrueRedshift *gsdtr, const gdouble z_min, const gdouble z_max)
{
  g_error ("_nc_galaxy_sd_true_redshift_set_lim: not implemented");

  return FALSE;
}

static gboolean
_nc_galaxy_sd_true_redshift_get_lim (NcGalaxySDTrueRedshift *gsdtr, gdouble *z_min, gdouble *z_max)
{
  g_error ("_nc_galaxy_sd_true_redshift_get_lim: not implemented");

  return FALSE;
}

/* LCOV_EXCL_STOP */

static void
nc_galaxy_sd_true_redshift_class_init (NcGalaxySDTrueRedshiftClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_sd_true_redshift_set_property;
  model_class->get_property = &_nc_galaxy_sd_true_redshift_get_property;
  object_class->finalize    = &_nc_galaxy_sd_true_redshift_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy sample redshift distribution", "GalaxySDTrueRedshift");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  /**
   * NcGalaxySDTrueRedshift:lim:
   *
   * Galaxy sample redshift distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LIM,
                                   g_param_spec_boxed ("lim",
                                                       NULL,
                                                       "Galaxy sample redshift distribution limits",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_mset_model_register_id (model_class, "NcGalaxySDTrueRedshift", "Galaxy sample redshift distribution", NULL, FALSE, nc_galaxy_sd_obs_redshift_id ());
  ncm_model_class_check_params_info (model_class);

  klass->gen     = &_nc_galaxy_sd_true_redshift_gen;
  klass->integ   = &_nc_galaxy_sd_true_redshift_integ;
  klass->set_lim = &_nc_galaxy_sd_true_redshift_set_lim;
  klass->get_lim = &_nc_galaxy_sd_true_redshift_get_lim;
}

/**
 * nc_galaxy_sd_true_redshift_ref:
 * @gsdtr: a #NcGalaxySDTrueRedshift
 *
 * Increases the reference count of @gsdtr by one.
 *
 * Returns: (transfer full): @gsdtr.
 */
NcGalaxySDTrueRedshift *
nc_galaxy_sd_true_redshift_ref (NcGalaxySDTrueRedshift *gsdtr)
{
  return g_object_ref (gsdtr);
}

/**
 * nc_galaxy_sd_true_redshift_free:
 * @gsdtr: a #NcGalaxySDTrueRedshift
 *
 * Decreases the reference count of @gsdtr by one.
 *
 */
void
nc_galaxy_sd_true_redshift_free (NcGalaxySDTrueRedshift *gsdtr)
{
  g_object_unref (gsdtr);
}

/**
 * nc_galaxy_sd_true_redshift_clear:
 * @gsdtr: a #NcGalaxySDTrueRedshift
 *
 * Decreases the reference count of @gsdtr by one, and sets the
 * pointer @gsdtr to NULL.
 *
 */
void
nc_galaxy_sd_true_redshift_clear (NcGalaxySDTrueRedshift **gsdtr)
{
  g_clear_object (gsdtr);
}

/**
 * nc_galaxy_sd_true_redshift_set_lim:
 * @gsdtr: a #NcGalaxySDTrueRedshift
 * @z_min: a #gdouble representing minimum redshift
 * @z_max: a #gdouble representing maximum redshift
 *
 * Sets the redshift limits of the galaxy sample redshift distribution.
 *
 */
gboolean
nc_galaxy_sd_true_redshift_set_lim (NcGalaxySDTrueRedshift *gsdtr, const gdouble z_min, const gdouble z_max)
{
  return NC_GALAXY_SD_TRUE_REDSHIFT_GET_CLASS (gsdtr)->set_lim (gsdtr, z_min, z_max);
}

/**
 * nc_galaxy_sd_true_redshift_get_lim:
 * @gsdtr: a #NcGalaxySDTrueRedshift
 * @z_min: a #gdouble pointer representing minimum redshift
 * @z_max: a #gdouble pointer representing maximum redshift
 *
 * Gets the redshift limits of the galaxy sample redshift distribution.
 *
 */
gboolean
nc_galaxy_sd_true_redshift_get_lim (NcGalaxySDTrueRedshift *gsdtr, gdouble *z_min, gdouble *z_max)
{
  return NC_GALAXY_SD_TRUE_REDSHIFT_GET_CLASS (gsdtr)->get_lim (gsdtr, z_min, z_max);
}

/**
 * nc_galaxy_sd_true_redshift_gen:
 * @gsdtr: a #NcGalaxySDTrueRedshift
 * @rng: a #NcmRNG
 *
 * Generates a redshift value from the galaxy sample redshift distribution.
 *
 * Returns: the generated redshift.
 */
gdouble
nc_galaxy_sd_true_redshift_gen (NcGalaxySDTrueRedshift *gsdtr, NcmRNG *rng)
{
  return NC_GALAXY_SD_TRUE_REDSHIFT_GET_CLASS (gsdtr)->gen (gsdtr, rng);
}

/**
 * nc_galaxy_sd_true_redshift_integ:
 * @gsdtr: a #NcGalaxySDTrueRedshift
 * @z: the redshift
 *
 * Evaluates the galaxy sample redshift distribution at redshift @z.
 *
 * Returns: the probability density at $z$, $P(z)$.
 */
gdouble
nc_galaxy_sd_true_redshift_integ (NcGalaxySDTrueRedshift *gsdtr, gdouble z)
{
  return NC_GALAXY_SD_TRUE_REDSHIFT_GET_CLASS (gsdtr)->integ (gsdtr, z);
}

typedef struct _NcGalaxySDTrueRedshiftIntegData
{
  NcGalaxySDTrueRedshift *gsdtr;
  gdouble abstol;
} NcGalaxySDTrueRedshiftIntegData;

static gdouble
_nc_galaxy_sd_true_redshift_integ_f (gdouble z, gpointer user_data)
{
  NcGalaxySDTrueRedshiftIntegData *data = (NcGalaxySDTrueRedshiftIntegData *) user_data;
  NcGalaxySDTrueRedshift *gsdtr         = data->gsdtr;

  return -2.0 * log (nc_galaxy_sd_true_redshift_integ (gsdtr, z) + data->abstol);
}

/**
 * nc_galaxy_sd_true_redshift_dist:
 * @gsdtr: a #NcGalaxySDTrueRedshift
 * @reltol: relative tolerance
 * @abstol: absolute tolerance
 *
 * Computes the redshift distribution of the galaxy sample redshift distribution.
 *
 * Returns: (transfer full): the redshift distribution object #NcmStatsDist1d.
 */
NcmStatsDist1d *
nc_galaxy_sd_true_redshift_dist (NcGalaxySDTrueRedshift *gsdtr, const gdouble reltol, const gdouble abstol)
{
  NcmSpline *s                               = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  const gdouble z_min                        = 0.0;
  const gdouble z_max                        = 20.0;
  NcGalaxySDTrueRedshiftIntegData integ_data = {gsdtr, abstol};
  gsl_function F;

  F.function = _nc_galaxy_sd_true_redshift_integ_f;
  F.params   = &integ_data;

  ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, z_min, z_max, 0, reltol);

  {
    NcmStatsDist1dSpline *spline = ncm_stats_dist1d_spline_new (s);

    g_object_set (G_OBJECT (spline), "reltol", reltol, "abstol", abstol, NULL);

    ncm_spline_free (s);

    ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (spline));

    return NCM_STATS_DIST1D (spline);
  }
}

