/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_galaxy_redshift.c
 *
 *  Tue April 17 11:12:53 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_redshift.c
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
 * Copyright (C) 2018 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_galaxy_redshift
 * @title: NcGalaxyRedshift
 * @short_description: Abstract class describing galaxy redshifts.
 *
 * Abstract class used to define a generic galaxy redshift probability
 * distribution $P_g(z)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_galaxy_redshift.h"

struct _NcGalaxyRedshiftPrivate
{
	gint placeholder;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyRedshift, nc_galaxy_redshift, G_TYPE_OBJECT);

static void
nc_galaxy_redshift_init (NcGalaxyRedshift *gz)
{
	gz->priv = nc_galaxy_redshift_get_instance_private (gz);
	
}

static void
_nc_galaxy_redshift_finalize (GObject *object)
{

	/* Chain up : end */
	G_OBJECT_CLASS (nc_galaxy_redshift_parent_class)->finalize (object);
}

static gboolean _nc_galaxy_redshift_has_dist (NcGalaxyRedshift *gz) { g_error ("_nc_galaxy_redshift_has_dist: method not implemented."); return FALSE; }
static gdouble _nc_galaxy_redshift_mode (NcGalaxyRedshift *gz) { g_error ("_nc_galaxy_redshift_mode: method not implemented."); return 0.0; }
static guint _nc_galaxy_redshift_nintervals (NcGalaxyRedshift *gz) { g_error ("_nc_galaxy_redshift_nintervals: method not implemented."); return 0; }
static gdouble _nc_galaxy_redshift_interval_weight (NcGalaxyRedshift *gz, const guint di) { g_error ("_nc_galaxy_redshift_interval_weight: method not implemented."); return 0; }
static void _nc_galaxy_redshift_pdf_limits (NcGalaxyRedshift *gz, const guint di, gdouble *zmin, gdouble *zmax) { g_error ("_nc_galaxy_redshift_pdf_limits: method not implemented."); }
static gdouble _nc_galaxy_redshift_pdf (NcGalaxyRedshift *gz, const guint di, const gdouble z) { g_error ("_nc_galaxy_redshift_pdf: method not implemented."); return 0.0; }
static gdouble _nc_galaxy_redshift_gen (NcGalaxyRedshift *gz, NcmRNG *rng) { g_error ("_nc_galaxy_redshift_gen: method not implemented."); return 0.0; }
static gdouble _nc_galaxy_redshift_quantile (NcGalaxyRedshift *gz, const gdouble q) { g_error ("_nc_galaxy_redshift_quantile: method not implemented."); return 0.0; }

static void
nc_galaxy_redshift_class_init (NcGalaxyRedshiftClass *klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);

	object_class->finalize = &_nc_galaxy_redshift_finalize;

	klass->has_dist        = &_nc_galaxy_redshift_has_dist;
	klass->mode            = &_nc_galaxy_redshift_mode;
	klass->nintervals      = &_nc_galaxy_redshift_nintervals;
	klass->interval_weight = &_nc_galaxy_redshift_interval_weight;
	klass->pdf_limits      = &_nc_galaxy_redshift_pdf_limits;
	klass->pdf             = &_nc_galaxy_redshift_pdf;
	klass->gen             = &_nc_galaxy_redshift_gen;
	klass->quantile        = &_nc_galaxy_redshift_quantile;
}

/**
 * nc_galaxy_redshift_ref:
 * @gz: a #NcGalaxyRedshift
 *
 * Increase the reference of @gz by one.
 *
 * Returns: (transfer full): @gz.
 */
NcGalaxyRedshift *
nc_galaxy_redshift_ref (NcGalaxyRedshift *gz)
{
  return g_object_ref (gz);
}

/**
 * nc_galaxy_redshift_free:
 * @gz: a #NcGalaxyRedshift
 *
 * Decrease the reference count of @gz by one.
 *
 */
void
nc_galaxy_redshift_free (NcGalaxyRedshift *gz)
{
  g_object_unref (gz);
}

/**
 * nc_galaxy_redshift_clear:
 * @gz: a #NcGalaxyRedshift
 *
 * Decrease the reference count of @gz by one, and sets the pointer *@gz to
 * NULL.
 *
 */
void
nc_galaxy_redshift_clear (NcGalaxyRedshift **gz)
{
  g_clear_object (gz);
}

/**
 * nc_galaxy_redshift_has_dist: (virtual has_dist)
 * @gz: a #NcGalaxyRedshift
 *
 * Returns: whether the galaxy has a redshift distribution.
 */
/**
 * nc_galaxy_redshift_mode: (virtual mode)
 * @gz: a #NcGalaxyRedshift
 *
 * Returns: the value of the distribution mode or the spectroscopic redshift when available.
 */
/**
 * nc_galaxy_redshift_nintervals: (virtual nintervals)
 * @gz: a #NcGalaxyRedshift
 *
 * Returns: the number of disjointed intervals in the distribution.
 */
/**
 * nc_galaxy_redshift_interval_weight: (virtual interval_weight)
 * @gz: a #NcGalaxyRedshift
 * @di: distribution index $i$
 *
 * Returns: the weight associated to interval $i$.
 */
/**
 * nc_galaxy_redshift_pdf_limits: (virtual pdf_limits)
 * @gz: a #NcGalaxyRedshift
 * @di: distribution index $i$
 * @zmin: (out caller-allocates): the redshift lower limit
 * @zmax: (out caller-allocates): the redshift upper limit
 *
 * This method provides the limits of the redshift probability density
 * $p_i(z)$.
 * 
 */
/**
 * nc_galaxy_redshift_pdf: (virtual pdf)
 * @gz: a #NcGalaxyRedshift
 * @di: distribution index $i$
 * @z: the redshift $z$
 *
 * Returns: the probability density at @z, $p_{i}(z)$.
 */
/**
 * nc_galaxy_redshift_gen: (virtual gen)
 * @gz: a #NcGalaxyRedshift
 * @rng: a #NcmRNG
 * 
 * Generates a redshift from the distribution using @rng.
 * 
 * Returns: the generated value $z$.
 */
/**
 * nc_galaxy_redshift_quantile: (virtual quantile)
 * @gz: a #NcGalaxyRedshift
 * @q: the quantile $q \in [0, 1]$
 * 
 * Computes the $q$ quantile.
 * 
 * Returns: the $q$ quantile.
 */
