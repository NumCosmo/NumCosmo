/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_galaxy_redshift_spec.c
 *
 *  Tue April 17 14:26:17 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_spec.c
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
 * SECTION:nc_galaxy_redshift_spec
 * @title: NcGalaxyRedshiftSpec
 * @short_description: Class describing spectroscopic galaxy redshifts.
 *
 * Class used to define a generic galaxy redshift probability distribution 
 * $P_g(z) = \delta(z - z_\mathrm{spec})$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_galaxy_redshift_spec.h"

struct _NcGalaxyRedshiftSpecPrivate
{
	gdouble z_spec;
};

enum
{
	PROP_0,
	PROP_Z_SPEC
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyRedshiftSpec, nc_galaxy_redshift_spec, NC_TYPE_GALAXY_REDSHIFT);

static void
nc_galaxy_redshift_spec_init (NcGalaxyRedshiftSpec *gzs)
{
	NcGalaxyRedshiftSpecPrivate * const self = gzs->priv = nc_galaxy_redshift_spec_get_instance_private (gzs);

	self->z_spec = -1.0;
}

static void
_nc_galaxy_redshift_spec_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
	NcGalaxyRedshiftSpec *gzs = NC_GALAXY_REDSHIFT_SPEC (object);
	g_return_if_fail (NC_IS_GALAXY_REDSHIFT_SPEC (object));

	switch (prop_id)
	{
		case PROP_Z_SPEC:
			nc_galaxy_redshift_spec_set_z (gzs, g_value_get_double (value));
			break;
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
_nc_galaxy_redshift_spec_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
	NcGalaxyRedshiftSpec *gzs = NC_GALAXY_REDSHIFT_SPEC (object);
	g_return_if_fail (NC_IS_GALAXY_REDSHIFT_SPEC (object));

	switch (prop_id)
	{
		case PROP_Z_SPEC:
			g_value_set_double (value, nc_galaxy_redshift_spec_get_z (gzs));
			break;
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
_nc_galaxy_redshift_spec_finalize (GObject *object)
{

	/* Chain up : end */
	G_OBJECT_CLASS (nc_galaxy_redshift_spec_parent_class)->finalize (object);
}

static gboolean _nc_galaxy_redshift_spec_has_dist (NcGalaxyRedshift *gz);
static gdouble _nc_galaxy_redshift_spec_mode (NcGalaxyRedshift *gz);
static guint _nc_galaxy_redshift_spec_nintervals (NcGalaxyRedshift *gz);
static gdouble _nc_galaxy_redshift_spec_interval_weight (NcGalaxyRedshift *gz, const guint di);
static void _nc_galaxy_redshift_spec_pdf_limits (NcGalaxyRedshift *gz, const guint di, gdouble *zmin, gdouble *zmax);
static gdouble _nc_galaxy_redshift_spec_pdf (NcGalaxyRedshift *gz, const guint di, const gdouble z);
static gdouble _nc_galaxy_redshift_spec_gen (NcGalaxyRedshift *gz, NcmRNG *rng);
static gdouble _nc_galaxy_redshift_spec_quantile (NcGalaxyRedshift *gz, const gdouble q);

static void
nc_galaxy_redshift_spec_class_init (NcGalaxyRedshiftSpecClass *klass)
{
	GObjectClass* object_class      = G_OBJECT_CLASS (klass);
	NcGalaxyRedshiftClass *gz_class = NC_GALAXY_REDSHIFT_CLASS (klass);

	object_class->set_property = &_nc_galaxy_redshift_spec_set_property;
	object_class->get_property = &_nc_galaxy_redshift_spec_get_property;
	object_class->finalize     = &_nc_galaxy_redshift_spec_finalize;

	g_object_class_install_property (object_class,
	                                 PROP_Z_SPEC,
	                                 g_param_spec_double ("z-spec",
	                                                      NULL,
	                                                      "Spectroscopic redshift",
	                                                      0.0, G_MAXDOUBLE, 0.0,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	gz_class->has_dist        = &_nc_galaxy_redshift_spec_has_dist;
	gz_class->mode            = &_nc_galaxy_redshift_spec_mode;
	gz_class->nintervals      = &_nc_galaxy_redshift_spec_nintervals;
	gz_class->interval_weight = &_nc_galaxy_redshift_spec_interval_weight;
	gz_class->pdf_limits      = &_nc_galaxy_redshift_spec_pdf_limits;
  gz_class->pdf             = &_nc_galaxy_redshift_spec_pdf;
  gz_class->gen             = &_nc_galaxy_redshift_spec_gen;
  gz_class->quantile        = &_nc_galaxy_redshift_spec_quantile;
	
}

static gboolean 
_nc_galaxy_redshift_spec_has_dist (NcGalaxyRedshift *gz)
{
	return FALSE;
}

static gdouble 
_nc_galaxy_redshift_spec_mode (NcGalaxyRedshift *gz)
{
	NcGalaxyRedshiftSpec *gzs = NC_GALAXY_REDSHIFT_SPEC (gz);
	NcGalaxyRedshiftSpecPrivate * const self = gzs->priv;

	return self->z_spec;
}

static guint 
_nc_galaxy_redshift_spec_nintervals (NcGalaxyRedshift *gz)
{
	return 1;
}

static gdouble 
_nc_galaxy_redshift_spec_interval_weight (NcGalaxyRedshift *gz, const guint di)
{
	g_assert_cmpuint (di, ==, 0);
	return 1.0;
}

static void 
_nc_galaxy_redshift_spec_pdf_limits (NcGalaxyRedshift *gz, const guint di, gdouble *zmin, gdouble *zmax)
{
	NcGalaxyRedshiftSpec *gzs = NC_GALAXY_REDSHIFT_SPEC (gz);
	NcGalaxyRedshiftSpecPrivate * const self = gzs->priv;

	g_assert_cmpuint (di, ==, 0);
	
	zmin[0] = zmax[0] = self->z_spec;
}

static gdouble 
_nc_galaxy_redshift_spec_pdf (NcGalaxyRedshift *gz, const guint di, const gdouble z)
{
	g_error ("_nc_galaxy_redshift_spec_pdf: delta distribution not implementable!");
	return 0.0;
}

static gdouble 
_nc_galaxy_redshift_spec_gen (NcGalaxyRedshift *gz, NcmRNG *rng)
{
	NcGalaxyRedshiftSpec *gzs = NC_GALAXY_REDSHIFT_SPEC (gz);
	NcGalaxyRedshiftSpecPrivate * const self = gzs->priv;

	return self->z_spec;
}

static gdouble 
_nc_galaxy_redshift_spec_quantile (NcGalaxyRedshift *gz, const gdouble q)
{
	NcGalaxyRedshiftSpec *gzs = NC_GALAXY_REDSHIFT_SPEC (gz);
	NcGalaxyRedshiftSpecPrivate * const self = gzs->priv;

	return self->z_spec;
}

/**
 * nc_galaxy_redshift_spec_new:
 * @z_spec: the spectroscopic redshift $z_\mathrm{spec}$
 * 
 * Creates a new #NcGalaxyRedshiftSpec with $z_\mathrm{spec}$ = @z_spec.
 *
 * Returns: (transfer full): The newly created #NcGalaxyRedshiftSpec.
 */
NcGalaxyRedshiftSpec *
nc_galaxy_redshift_spec_new (const gdouble z_spec)
{
	NcGalaxyRedshiftSpec *gzs = g_object_new (NC_TYPE_GALAXY_REDSHIFT_SPEC, 
	                                          "z-spec", z_spec,
	                                          NULL);

	return gzs;
}

/**
 * nc_galaxy_redshift_spec_ref:
 * @gzs: a #NcGalaxyRedshiftSpec
 *
 * Increase the reference of @gzs by one.
 *
 * Returns: (transfer full): @gzs.
 */
NcGalaxyRedshiftSpec *
nc_galaxy_redshift_spec_ref (NcGalaxyRedshiftSpec *gzs)
{
  return g_object_ref (gzs);
}

/**
 * nc_galaxy_redshift_spec_free:
 * @gzs: a #NcGalaxyRedshiftSpec
 *
 * Decrease the reference count of @gzs by one.
 *
 */
void
nc_galaxy_redshift_spec_free (NcGalaxyRedshiftSpec *gzs)
{
  g_object_unref (gzs);
}

/**
 * nc_galaxy_redshift_spec_clear:
 * @gzs: a #NcGalaxyRedshiftSpec
 *
 * Decrease the reference count of @gzs by one, and sets the pointer *@gzs to
 * NULL.
 *
 */
void
nc_galaxy_redshift_spec_clear (NcGalaxyRedshiftSpec **gzs)
{
  g_clear_object (gzs);
}

/**
 * nc_galaxy_redshift_spec_set_z:
 * @gzs: a #NcGalaxyRedshiftSpec
 * @z_spec: the spectroscopic redshift $z_\mathrm{spec}$
 * 
 * Sets $z_\mathrm{spec}$ = @z_spec.
 *
 */
void
nc_galaxy_redshift_spec_set_z (NcGalaxyRedshiftSpec *gzs, const gdouble z_spec)
{
	NcGalaxyRedshiftSpecPrivate * const self = gzs->priv;
	
	self->z_spec = z_spec;
}

/**
 * nc_galaxy_redshift_spec_get_z:
 * @gzs: a #NcGalaxyRedshiftSpec
 * 
 * Gets $z_\mathrm{spec}$.
 *
 * Returns: $z_\mathrm{spec}$.
 */
gdouble
nc_galaxy_redshift_spec_get_z (NcGalaxyRedshiftSpec *gzs)
{
	NcGalaxyRedshiftSpecPrivate * const self = gzs->priv;
	
	return self->z_spec;
}
