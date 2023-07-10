/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_sd_position_srd_y1.c
 *
 *  Tue June 22 14:59:27 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_sd_position_srd_y1.c
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * SECTION:nc_galaxy_sd_position_srd_y1
 * @title: NcGalaxySDPositionSRDY1
 * @short_description: Class describing galaxy sample position distributions with SRD year 1 distribution
 * @stability: Unstable
 *
 *
 * Class defining a galaxy sample position distribution with SRD year 1
 * probability distribution $P(z)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

/* #include "nc_enum_types.h" */
#include "galaxy/nc_galaxy_sd_position_srd_y1.h"
#include "galaxy/nc_galaxy_sd_position.h"
#include "math/ncm_stats_dist1d.h"
#include "math/ncm_stats_dist1d_spline.h"
#include "math/ncm_spline.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_spline_cubic.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include <math.h>
#include <gsl/gsl_math.h>

struct _NcGalaxySDPositionSRDY1Private
{
  NcmVector *z_lim;
  NcmVector *r_lim;
  NcmStatsDist1d *z_dist;
};

enum
{
  PROP_0,
  PROP_Z_LIM,
  PROP_R_LIM,
  PROP_Z_DIST,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySDPositionSRDY1, nc_galaxy_sd_position_srd_y1, NC_TYPE_GALAXY_SD_POSITION);

static void
nc_galaxy_sd_position_srd_y1_init (NcGalaxySDPositionSRDY1 *gsdpsrdy1)
{
  NcGalaxySDPositionSRDY1Private * const self = gsdpsrdy1->priv = nc_galaxy_sd_position_srd_y1_get_instance_private (gsdpsrdy1);

  self->z_lim  = NULL;
  self->r_lim  = NULL;
  self->z_dist = NULL;
}

static void
_nc_galaxy_sd_position_srd_y1_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySDPositionSRDY1 *gsdpsrdy1 = NC_GALAXY_SD_POSITION_SRD_Y1 (object);

  g_return_if_fail (NC_IS_GALAXY_SD_POSITION_SRD_Y1 (object));

  switch (prop_id)
  {
    case PROP_Z_LIM:
      nc_galaxy_sd_position_srd_y1_set_z_lim (gsdpsrdy1, g_value_get_object (value));
      break;
    case PROP_R_LIM:
      nc_galaxy_sd_position_srd_y1_set_r_lim (gsdpsrdy1, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_sd_position_srd_y1_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySDPositionSRDY1 *gsdpsrdy1 = NC_GALAXY_SD_POSITION_SRD_Y1 (object);

  g_return_if_fail (NC_IS_GALAXY_SD_POSITION_SRD_Y1 (object));

  switch (prop_id)
  {
    case PROP_Z_LIM:
      g_value_set_object (value, nc_galaxy_sd_position_srd_y1_peek_z_lim (gsdpsrdy1));
      break;
    case PROP_R_LIM:
      g_value_set_object (value, nc_galaxy_sd_position_srd_y1_peek_r_lim (gsdpsrdy1));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}


static void
_nc_galaxy_sd_position_srd_y1_dispose (GObject *object)
{
  NcGalaxySDPositionSRDY1 *gsdpsrdy1          = NC_GALAXY_SD_POSITION_SRD_Y1 (object);
  NcGalaxySDPositionSRDY1Private * const self = gsdpsrdy1->priv;
  
  ncm_vector_clear (&self->z_lim);
  ncm_vector_clear (&self->r_lim);

  G_OBJECT_CLASS (nc_galaxy_sd_position_srd_y1_parent_class)->dispose (object);
}

static void
_nc_galaxy_sd_position_srd_y1_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_sd_position_srd_y1_parent_class)->finalize (object);
}

static void _nc_galaxy_sd_position_srd_y1_gen_r (NcGalaxySDPosition *gsdp, NcmRNG *rng, gdouble *gen_r);
static void _nc_galaxy_sd_position_srd_y1_gen_z (NcGalaxySDPosition *gsdp, NcmRNG *rng, gdouble *gen_z);
static gdouble _nc_galaxy_sd_position_srd_y1_integ (NcGalaxySDPosition *gsdp, NcmVector *pos);

static void
nc_galaxy_sd_position_srd_y1_class_init (NcGalaxySDPositionSRDY1Class *klass)
{
  NcGalaxySDPositionClass *sd_position_class = NC_GALAXY_SD_POSITION_CLASS (klass);
  GObjectClass *object_class                 = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_sd_position_srd_y1_set_property;
  object_class->get_property = &_nc_galaxy_sd_position_srd_y1_get_property;
  object_class->dispose      = &_nc_galaxy_sd_position_srd_y1_dispose;
  object_class->finalize     = &_nc_galaxy_sd_position_srd_y1_finalize;

  /**
   * NcGalaxySDPositionSRDY1:z-lim:
   *
   * Galaxy sample redshift distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_LIM,
                                   g_param_spec_object ("z-lim",
                                                        NULL,
                                                        "Galaxy sample redshift distribution limits",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxySDPositionSRDY1:R-lim:
   *
   * Galaxy sample radius distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_R_LIM,
                                   g_param_spec_object ("r-lim",
                                                        NULL,
                                                        "Galaxy sample radius distribution limits",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  sd_position_class->gen_r = &_nc_galaxy_sd_position_srd_y1_gen_r;
  sd_position_class->gen_z = &_nc_galaxy_sd_position_srd_y1_gen_z;
  sd_position_class->integ = &_nc_galaxy_sd_position_srd_y1_integ;
}


static void
_nc_galaxy_sd_position_srd_y1_gen_r (NcGalaxySDPosition *gsdp, NcmRNG *rng, gdouble *gen_r)
{
  NcGalaxySDPositionSRDY1 *gsdpsrdy1          = NC_GALAXY_SD_POSITION_SRD_Y1 (gsdp);
  NcGalaxySDPositionSRDY1Private * const self = gsdpsrdy1->priv;
  const gdouble r_lb                          = ncm_vector_get (self->r_lim, 0);
  const gdouble r_ub                          = ncm_vector_get (self->r_lim, 1);
  const gdouble r_lb2                         = r_lb * r_lb;
  const gdouble r_ub2                         = r_ub * r_ub;
  gdouble cumul_gen                           = ncm_rng_uniform_gen (rng, 0.0, 1.0);
  gen_r[0]                                    = sqrt (cumul_gen * (r_ub2 - r_lb2) + r_lb2);
}

static void
_nc_galaxy_sd_position_srd_y1_gen_z (NcGalaxySDPosition *gsdp, NcmRNG *rng, gdouble *gen_z)
{
  NcGalaxySDPositionSRDY1 *gsdpsrdy1          = NC_GALAXY_SD_POSITION_SRD_Y1 (gsdp);
  NcGalaxySDPositionSRDY1Private * const self = gsdpsrdy1->priv;
  gen_z[0]                                    = ncm_stats_dist1d_gen (self->z_dist, rng);
}

static gdouble
_nc_galaxy_sd_position_srd_y1_integ (NcGalaxySDPosition *gsdp, NcmVector *pos)
{
  return 0.0;
}

/**
 * nc_galaxy_sd_position_srd_y1_new:
 *
 * Creates a new #NcGalaxySDPositionSRDY1
 *
 * Returns: (transfer full): a new NcGalaxySDPositionSRDY1.
 */
NcGalaxySDPositionSRDY1 *
nc_galaxy_sd_position_srd_y1_new ()
{
  NcGalaxySDPositionSRDY1 *gsdpsrdy1 = g_object_new (NC_TYPE_GALAXY_SD_POSITION_SRD_Y1,
                                                     NULL);

  return gsdpsrdy1;
}

/**
 * nc_galaxy_sd_position_srd_y1_ref:
 * @gsdpsrdy1: a #NcGalaxySDPositionSRDY1
 *
 * Increase the reference of @gsdpsrdy1 by one.
 *
 * Returns: (transfer full): @gsdpsrdy1.
 */
NcGalaxySDPositionSRDY1 *
nc_galaxy_sd_position_srd_y1_ref (NcGalaxySDPositionSRDY1 *gsdpsrdy1)
{
  return g_object_ref (gsdpsrdy1);
}

/**
 * nc_galaxy_sd_position_srd_y1_free:
 * @gsdpsrdy1: a #NcGalaxySDPositionSRDY1
 *
 * Decrease the reference count of @gsdpsrdy1 by one.
 *
 */
void
nc_galaxy_sd_position_srd_y1_free (NcGalaxySDPositionSRDY1 *gsdpsrdy1)
{
  g_object_unref (gsdpsrdy1);
}

/**
 * nc_galaxy_sd_position_srd_y1_clear:
 * @gsdpsrdy1: a #NcGalaxySDPositionSRDY1
 *
 * Decrease the reference count of @gsdpsrdy1 by one, and sets the pointer *@gsdpsrdy1 to
 * NULL.
 *
 */
void
nc_galaxy_sd_position_srd_y1_clear (NcGalaxySDPositionSRDY1 **gsdpsrdy1)
{
  g_clear_object (gsdpsrdy1);
}

/**
 * nc_galaxy_sd_position_srd_y1_set_z_lim:
 * @gsdpsrdy1: a #NcGalaxySDPositionSRDY1
 * @lim: a #NcmVector
 *
 * Sets the redshift limits @lim.
 */
void
nc_galaxy_sd_position_srd_y1_set_z_lim (NcGalaxySDPositionSRDY1 *gsdpsrdy1, NcmVector *lim)
{
  NcGalaxySDPositionSRDY1Private * const self = gsdpsrdy1->priv;

  g_assert_cmpuint (ncm_vector_len (lim), ==, 2);

  ncm_vector_clear (&self->z_lim);
  ncm_stats_dist1d_clear (&self->z_dist);

  gdouble z_ll = ncm_vector_get (lim, 0);
  gdouble z_ul = ncm_vector_get (lim, 1);

  gdouble
  m2lnp (gdouble z, void * p)
  {
    (void)(p);
    gdouble alpha = 0.78;
    gdouble beta = 2.0;
    gdouble z0 = 0.13;

    return - 2 * log (pow (z, beta) * exp (- pow (z/z0, alpha)));
  }

  gsl_function F;

  F.function = &m2lnp;
  F.params = 0;

  NcmSpline *spline = ncm_spline_cubic_notaknot_new ();
  ncm_spline_set_func (spline, NCM_SPLINE_FUNCTION_SPLINE, &F, z_ll, z_ul, 10000, 0.01);

  NcmStatsDist1dSpline *dist = ncm_stats_dist1d_spline_new (spline);

  ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (dist));

  self->z_lim = ncm_vector_ref (lim);
  self->z_dist = NCM_STATS_DIST1D (dist);
}


/**
 * nc_galaxy_sd_position_srd_y1_peek_z_lim:
 * @gsdpsrdy1: a #NcGalaxySDPositionSRDY1
 *
 * Gets the redshift limits.
 *
 * Returns: (transfer none): the redshift limits.
 */
NcmVector *
nc_galaxy_sd_position_srd_y1_peek_z_lim (NcGalaxySDPositionSRDY1 *gsdpsrdy1)
{
  NcGalaxySDPositionSRDY1Private * const self = gsdpsrdy1->priv;

  return self->z_lim;
}

/**
 * nc_galaxy_sd_position_srd_y1_set_r_lim:
 * @gsdpsrdy1: a #NcGalaxySDPositionSRDY1
 * @lim: a #NcmVector
 *
 * Sets the radius limits @lim.
 */
void
nc_galaxy_sd_position_srd_y1_set_r_lim (NcGalaxySDPositionSRDY1 *gsdpsrdy1, NcmVector *lim)
{
  NcGalaxySDPositionSRDY1Private * const self = gsdpsrdy1->priv;

  g_assert_cmpuint (ncm_vector_len (lim), ==, 2);

  ncm_vector_clear (&self->r_lim);

  self->r_lim = ncm_vector_ref (lim);
}


/**
 * nc_galaxy_sd_position_srd_y1_peek_r_lim:
 * @gsdpsrdy1: a #NcGalaxySDPositionSRDY1
 *
 * Gets the radius limits.
 *
 * Returns: (transfer none): the radius limits.
 */
NcmVector *
nc_galaxy_sd_position_srd_y1_peek_r_lim (NcGalaxySDPositionSRDY1 *gsdpsrdy1)
{
  NcGalaxySDPositionSRDY1Private * const self = gsdpsrdy1->priv;

  return self->r_lim;
}