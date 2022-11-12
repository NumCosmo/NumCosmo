/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_sf_spherical_harmonics.c
 *
 *  Thu December 14 11:18:00 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sf_spherical_harmonics.c
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
};

G_DEFINE_BOXED_TYPE (NcmSFSphericalHarmonicsY, ncm_sf_spherical_harmonics_Y, ncm_sf_spherical_harmonics_Y_dup, ncm_sf_spherical_harmonics_Y_free);
G_DEFINE_BOXED_TYPE (NcmSFSphericalHarmonicsYArray, ncm_sf_spherical_harmonics_Y_array, ncm_sf_spherical_harmonics_Y_array_dup, ncm_sf_spherical_harmonics_Y_array_free);
G_DEFINE_TYPE (NcmSFSphericalHarmonics, ncm_sf_spherical_harmonics, G_TYPE_OBJECT);

static void
ncm_sf_spherical_harmonics_init (NcmSFSphericalHarmonics *spha)
{
	spha->lmax     = 0;
	spha->sqrt_n   = g_array_new (FALSE, FALSE, sizeof (gdouble));
	spha->sqrtm1_n = g_array_new (FALSE, FALSE, sizeof (gdouble));
	spha->K_array  = g_ptr_array_new ();
	spha->Klm_m    = -1;

	g_ptr_array_set_free_func (spha->K_array, (GDestroyNotify) g_array_unref);
}

static void
_ncm_sf_spherical_harmonics_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
	NcmSFSphericalHarmonics *spha = NCM_SF_SPHERICAL_HARMONICS (object);
	g_return_if_fail (NCM_IS_SF_SPHERICAL_HARMONICS (object));

	switch (prop_id)
	{
		case PROP_LMAX:
			ncm_sf_spherical_harmonics_set_lmax (spha, g_value_get_int (value));
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
			g_value_set_int (value, ncm_sf_spherical_harmonics_get_lmax (spha));
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
	g_clear_pointer (&spha->sqrtm1_n, g_array_unref);
	g_clear_pointer (&spha->K_array, g_ptr_array_unref);
	
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
	                                 g_param_spec_int ("lmax",
	                                                   NULL,
	                                                   "max l",
	                                                   0, G_MAXINT, 1024,
	                                                   G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_sf_spherical_harmonics_Y_new:
 * @spha: a #NcmSFSphericalHarmonics
 * @abstol: absolute tolerance
 *
 * FIXME 	
 *
 * Returns: FIXME
 */
NcmSFSphericalHarmonicsY *
ncm_sf_spherical_harmonics_Y_new (NcmSFSphericalHarmonics *spha, const gdouble abstol)
{
  NcmSFSphericalHarmonicsY *sphaY = g_slice_new0 (NcmSFSphericalHarmonicsY);
	
  sphaY->x        = 0.0;
  sphaY->sqrt1mx2 = 0.0;

	sphaY->l        = 0;
	sphaY->l0       = 0;
	sphaY->m        = 0;
	sphaY->Klm      = NULL;
	sphaY->Pl0m     = 0.0;
	sphaY->Pl0p1m   = 0.0;
	sphaY->Plm      = 0.0;
	sphaY->Plp1m    = 0.0;

	sphaY->spha     = ncm_sf_spherical_harmonics_ref (spha);
	sphaY->abstol   = abstol;

  return sphaY;
}

/**
 * ncm_sf_spherical_harmonics_Y_dup:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmSFSphericalHarmonicsY *
ncm_sf_spherical_harmonics_Y_dup (NcmSFSphericalHarmonicsY *sphaY)
{
  NcmSFSphericalHarmonicsY *sphaY_dup = ncm_sf_spherical_harmonics_Y_new (sphaY->spha, sphaY->abstol);
	sphaY_dup[0] = sphaY[0];
  return sphaY_dup;
}

/**
 * ncm_sf_spherical_harmonics_Y_free:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * FIXME
 *
 */
void
ncm_sf_spherical_harmonics_Y_free (NcmSFSphericalHarmonicsY *sphaY)
{
	ncm_sf_spherical_harmonics_clear (&sphaY->spha);
  g_slice_free (NcmSFSphericalHarmonicsY, sphaY);
}

/**
 * ncm_sf_spherical_harmonics_Y_array_new:
 * @spha: a #NcmSFSphericalHarmonics
 * @len: array length
 * @abstol: absolute tolerance
 * 
 * FIXME 	
 *
 * Returns: FIXME
 */
NcmSFSphericalHarmonicsYArray *
ncm_sf_spherical_harmonics_Y_array_new (NcmSFSphericalHarmonics *spha, const gint len, const gdouble abstol)
{
  NcmSFSphericalHarmonicsYArray *sphaYa = g_slice_new0 (NcmSFSphericalHarmonicsYArray);

	sphaYa->l        = 0;
	sphaYa->l0       = 0;
	sphaYa->m        = 0;
	sphaYa->Klm      = NULL;
	sphaYa->len      = len;

	g_assert_cmpuint (len, <=, NCM_SF_SPHERICAL_HARMONICS_MAX_LEN);

	sphaYa->spha     = ncm_sf_spherical_harmonics_ref (spha);
	sphaYa->abstol   = abstol;

  return sphaYa;
}

/**
 * ncm_sf_spherical_harmonics_Y_array_dup:
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmSFSphericalHarmonicsYArray *
ncm_sf_spherical_harmonics_Y_array_dup (NcmSFSphericalHarmonicsYArray *sphaYa)
{
  NcmSFSphericalHarmonicsYArray *sphaYa_dup = ncm_sf_spherical_harmonics_Y_array_new (sphaYa->spha, sphaYa->len, sphaYa->abstol);

	sphaYa_dup[0] = sphaYa[0];

  return sphaYa_dup;
}

/**
 * ncm_sf_spherical_harmonics_Y_array_free:
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 *
 * FIXME
 *
 */
void
ncm_sf_spherical_harmonics_Y_array_free (NcmSFSphericalHarmonicsYArray *sphaYa)
{
	ncm_sf_spherical_harmonics_clear (&sphaYa->spha);

  g_slice_free (NcmSFSphericalHarmonicsYArray, sphaYa);
}

/**
 * ncm_sf_spherical_harmonics_Y_next_l:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 * 
 * Move the recursion for $x$ to $l = l + 1$.
 *
 */
/**
 * ncm_sf_spherical_harmonics_Y_next_l2:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 * @Yblm: (array fixed-size=2) (element-type gdouble) (out caller-allocates): the four Yblm from l to l + 2  
 * 
 * Move the recursion for $x$ to $l = l + 2$.
 *
 */
/**
 * ncm_sf_spherical_harmonics_Y_next_l4:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 * @Yblm: (array fixed-size=4) (element-type gdouble) (out caller-allocates): the four Yblm from l to l + 3  
 * 
 * Move the recursion for $x$ to $l = l + 4$.
 *
 */
/**
 * ncm_sf_spherical_harmonics_Y_next_lnp2:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 * @Yblm: (array length=n) (element-type gdouble) (out caller-allocates): the four Yblm from l to l + n + 2
 * @n: an integer 
 * 
 * Move the recursion for $x$ to $l = l + n + 2$.
 *
 */
/**
 * ncm_sf_spherical_harmonics_Y_next_m:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 * 
 * Restart the recursion for $x$ at $l = m + 1,\; m = m + 1$.
 * If the value of $Ybll < a$ where $a$ is the absolute tolerance,
 * advance $l$ until the tolerance is reached.
 *
 */

/**
 * ncm_sf_spherical_harmonics_Y_get_lm:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Returns: the current value of $\bar{Y}_l^m (x)$.
 */
/**
 * ncm_sf_spherical_harmonics_Y_get_lp1m:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Returns: the current value of $\bar{Y}_{l+1}^m (x)$.
 */
/**
 * ncm_sf_spherical_harmonics_Y_get_x:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Returns: the current value of $x$.
 */
/**
 * ncm_sf_spherical_harmonics_Y_get_l:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Returns: the current value of $l$.
 */
/**
 * ncm_sf_spherical_harmonics_Y_get_m:
 * @sphaY: a #NcmSFSphericalHarmonicsY
 *
 * Returns: the current value of $m$.
 */

/**
 * ncm_sf_spherical_harmonics_new:
 * @lmax: $\ell_\mathrm{max}$
 * 
 * Creates a new #NcmSFSphericalHarmonics object.
 * 
 * Returns: a new #NcmSFSphericalHarmonics.
 */
NcmSFSphericalHarmonics *
ncm_sf_spherical_harmonics_new (const gint lmax)
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
 * Decrease the reference count of @spha by one, and sets the pointer *@spha to
 * NULL.
 *
 */
void
ncm_sf_spherical_harmonics_clear (NcmSFSphericalHarmonics **spha)
{
  g_clear_object (spha);
}

#define SN(n)   g_array_index (spha->sqrt_n, gdouble, (n))
#define SNM1(n) g_array_index (spha->sqrtm1_n, gdouble, (n))

/**
 * ncm_sf_spherical_harmonics_set_lmax:
 * @spha: a #NcmSFSphericalHarmonics
 * @lmax: maximum l
 *
 * Sets the maximum l to @lmax.
 * 
 */
void 
ncm_sf_spherical_harmonics_set_lmax (NcmSFSphericalHarmonics *spha, const gint lmax)
{
	if (lmax != spha->lmax)
	{
		const gint nmax_old = spha->sqrt_n->len;
		const gint nmax     = 2 * lmax + 3;
		gint n;

		g_array_set_size (spha->sqrt_n,   nmax + 1);
		g_array_set_size (spha->sqrtm1_n, nmax + 1);

		for (n = nmax_old; n <= nmax; n++)
		{
			const gdouble sqrt_n = sqrt (n);
			g_array_index (spha->sqrt_n, gdouble, n)   = sqrt_n;
			g_array_index (spha->sqrtm1_n, gdouble, n) = 1.0 / sqrt_n;
		}

		for (n = 0; n <= lmax; n++)
		{
			GArray *Km_array = NULL;
			gint l;

			if (n < spha->K_array->len)
			{
				Km_array = g_ptr_array_index (spha->K_array, n);
			}
			else if (n == spha->K_array->len)
			{
				Km_array = g_array_new (TRUE, TRUE, sizeof (NcmSFSphericalHarmonicsK));
				g_ptr_array_add (spha->K_array, Km_array);
			}
			else
			{
				g_assert_not_reached ();
			}

			for (l = Km_array->len + n; l < lmax; l++)
			{
				NcmSFSphericalHarmonicsK Km;
				const gint m       = n;
				const gint twol    = 2 * l;
				const gint lmm     = l - m;
				const gint lpm     = l + m;
				const gdouble pref = SN (twol + 5) * SNM1 (lpm + 2) * SNM1 (lmm + 2);

				Km.lp1  = pref * SN (twol + 3);
				Km.l    = pref * SN (lmm + 1) * SN (lpm + 1) * SNM1 (twol + 1);

				g_array_append_val (Km_array, Km);
			}

			g_array_set_size (Km_array, lmax - n);
		}

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
gint 
ncm_sf_spherical_harmonics_get_lmax (NcmSFSphericalHarmonics *spha)
{
	return spha->lmax;
}

/**
 * ncm_sf_spherical_harmonics_start_rec:
 * @spha: a #NcmSFSphericalHarmonics
 * @sphaY: a #NcmSFSphericalHarmonicsY
 * @theta: $\theta \in [0, \pi]$
 * 
 * Start recursion for $\theta$ at $l = 0,\; m = 0$.
 *
 */
/**
 * ncm_sf_spherical_harmonics_start_rec_array:
 * @spha: a #NcmSFSphericalHarmonics
 * @sphaYa: a #NcmSFSphericalHarmonicsYArray
 * @len: array length
 * @theta: (array length=len) (element-type gdouble): array of angles $\theta_i \in [0, \pi]$
 * 
 * Start recursion for the array $\theta_i$ at $l = 0,\; m = 0$. The array @theta must have
 * length @len.
 *
 */
/**
 * ncm_sf_spherical_harmonics_get_Klm: (skip)
 * @spha: a #NcmSFSphericalHarmonics
 * @l0: First $l_0$
 * @m: $m$ index
 *
 * Gets an array of #NcmSFSphericalHarmonicsK with the coeficients 
 * necessary to move the recurence from $(l_0, m)\; \to\; (l_\mathrm{max}, m)$. 
 * 
 * Returns: (array) (element-type NcmSFSphericalHarmonicsK): a pointer to an array of #NcmSFSphericalHarmonicsK.
 */

