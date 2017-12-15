/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_sf_spherical_harmonics.c
 *
 *  Thu December 14 11:18:00 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_sf_spherical_harmonics.c
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
 *
 * NumCosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * NumCosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION:ncm_sf_spherical_harmonics
 * @title: NcmSFSphericalHarmonics
 * @short_description: Spherical Harmonics object
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sf_spherical_harmonics.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
	PROP_0,
	PROP_LMAX,
	PROP_X,
};

G_DEFINE_TYPE (NcmSFSphericalHarmonics, ncm_sf_spherical_harmonics, G_TYPE_OBJECT);

static void
ncm_sf_spherical_harmonics_init (NcmSFSphericalHarmonics *spha)
{
	spha->lmax   = 0;
	spha->l      = 0;
	spha->m      = 0;
	spha->sqrt_n = g_array_new (FALSE, FALSE, sizeof (gdouble));
	spha->x      = 0.0;
}

static void
_ncm_sf_spherical_harmonics_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
	NcmSFSphericalHarmonics *spha = NCM_SF_SPHERICAL_HARMONICS (object);
	g_return_if_fail (NCM_IS_SF_SPHERICAL_HARMONICS (object));

	switch (prop_id)
	{
		case PROP_LMAX:
			ncm_sf_spherical_harmonics_set_lmax (spha, g_value_get_uint (value));
			break;
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
_ncm_sf_spherical_harmonics_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
	NcmSFSphericalHarmonics *spha = NCM_SF_SPHERICAL_HARMONICS (object);
	g_return_if_fail (NCM_IS_SF_SPHERICAL_HARMONICS (object));

	switch (prop_id)
	{
		case PROP_LMAX:
			g_value_set_uint (value, ncm_sf_spherical_harmonics_get_lmax (spha));
			break;
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
_ncm_sf_spherical_harmonics_dispose (GObject *object)
{
	NcmSFSphericalHarmonics *spha = NCM_SF_SPHERICAL_HARMONICS (object);

	g_clear_pointer (&spha->sqrt_n, g_array_unref);
	
	/* Chain up : end */
	G_OBJECT_CLASS (ncm_sf_spherical_harmonics_parent_class)->dispose (object);
}

static void
_ncm_sf_spherical_harmonics_finalize (GObject *object)
{

	/* Chain up : end */
	G_OBJECT_CLASS (ncm_sf_spherical_harmonics_parent_class)->finalize (object);
}

static void
ncm_sf_spherical_harmonics_class_init (NcmSFSphericalHarmonicsClass *klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);

	object_class->set_property = &_ncm_sf_spherical_harmonics_set_property;
	object_class->get_property = &_ncm_sf_spherical_harmonics_get_property;
	object_class->dispose      = &_ncm_sf_spherical_harmonics_dispose;
	object_class->finalize     = &_ncm_sf_spherical_harmonics_finalize;

	g_object_class_install_property (object_class,
	                                 PROP_LMAX,
	                                 g_param_spec_uint ("lmax",
	                                                    NULL,
	                                                    "max l",
	                                                    0, G_MAXUINT, 1024,
	                                                    G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_sf_spherical_harmonics_new:
 * 
 * Creates a new #NcmSFSphericalHarmonics object.
 * 
 * Returns: a new #NcmSFSphericalHarmonics.
 */
NcmSFSphericalHarmonics *
ncm_sf_spherical_harmonics_new (const guint lmax)
{
	NcmSFSphericalHarmonics *spha = g_object_new (NCM_TYPE_SF_SPHERICAL_HARMONICS,
	                                              "lmax", lmax,
	                                              NULL);
	return spha;
}


/**
 * ncm_sf_spherical_harmonics_ref:
 * @spha: a #NcmSFSphericalHarmonics
 *
 * Increase the reference of @spha by one.
 *
 * Returns: (transfer full): @spha.
 */
NcmSFSphericalHarmonics *
ncm_sf_spherical_harmonics_ref (NcmSFSphericalHarmonics *spha)
{
  return g_object_ref (spha);
}

/**
 * ncm_sf_spherical_harmonics_free:
 * @spha: a #NcmSFSphericalHarmonics
 *
 * Decrease the reference count of @spha by one.
 *
 */
void
ncm_sf_spherical_harmonics_free (NcmSFSphericalHarmonics *spha)
{
  g_object_unref (spha);
}

/**
 * ncm_sf_spherical_harmonics_clear:
 * @spha: a #NcmSFSphericalHarmonics
 *
 * Decrease the reference count of @spha by one, and sets the pointer *spha to
 * NULL.
 *
 */
void
ncm_sf_spherical_harmonics_clear (NcmSFSphericalHarmonics **spha)
{
  g_clear_object (spha);
}

/**
 * ncm_sf_spherical_harmonics_set_lmax:
 * @spha: a #NcmSFSphericalHarmonics
 * @lmax: maximum l
 *
 * Sets the maximum l to @lmax.
 * 
 */
void 
ncm_sf_spherical_harmonics_set_lmax (NcmSFSphericalHarmonics *spha, const guint lmax)
{
	if (lmax != spha->lmax)
	{
		const guint nmax = 2 * lmax + 1;
		guint n;

		g_array_set_size (spha->sqrt_n, nmax + 1);

		for (n = 0; n <= nmax; n++)
			g_array_index (spha->sqrt_n, gdouble, n) = sqrt (n);

		spha->lmax = lmax;
	}
}

/**
 * ncm_sf_spherical_harmonics_get_lmax:
 * @spha: a #NcmSFSphericalHarmonics
 *
 * Gets the maximum l to @lmax.
 * 
 * Returns: the currently used lmax.
 */
guint 
ncm_sf_spherical_harmonics_get_lmax (NcmSFSphericalHarmonics *spha)
{
	return spha->lmax;
}

/**
 * ncm_sf_spherical_harmonics_start_rec:
 * @spha: a #NcmSFSphericalHarmonics
 * @x: $x$
 * @sqrt1mx2: $\sqrt{1-x^2}$
 * 
 * Start recursion for $x$ at $l = 0,\; m = 0$.
 *
 */
/**
 * ncm_sf_spherical_harmonics_next_l:
 * @spha: a #NcmSFSphericalHarmonics
 * 
 * Move the recursion for $x$ to $l = l + 1$.
 *
 */
/**
 * ncm_sf_spherical_harmonics_next_m:
 * @spha: a #NcmSFSphericalHarmonics
 * 
 * Restart the recursion for $x$ at $l = m + 1,\; m = m + 1$.
 *
 */

/**
 * ncm_sf_spherical_harmonics_get_Yblm:
 * @spha: a #NcmSFSphericalHarmonics
 *
 * Returns: the current value of $\bar{Y}_l^m (x)$.
 */
/**
 * ncm_sf_spherical_harmonics_get_Yblp1m:
 * @spha: a #NcmSFSphericalHarmonics
 *
 * Returns: the current value of $\bar{Y}_{l+1}^m (x)$.
 */
/**
 * ncm_sf_spherical_harmonics_get_x:
 * @spha: a #NcmSFSphericalHarmonics
 *
 * Returns: the current value of $x$.
 */
/**
 * ncm_sf_spherical_harmonics_get_l:
 * @spha: a #NcmSFSphericalHarmonics
 *
 * Returns: the current value of $l$.
 */
/**
 * ncm_sf_spherical_harmonics_get_m:
 * @spha: a #NcmSFSphericalHarmonics
 *
 * Returns: the current value of $m$.
 */
