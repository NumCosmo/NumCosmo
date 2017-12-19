/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_sf_spherical_harmonics.h
 *
 *  Thu December 14 11:18:00 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_sf_spherical_harmonics.h
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

#ifndef _NCM_SF_SPHERICAL_HARMONICS_H_
#define _NCM_SF_SPHERICAL_HARMONICS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_c.h>

G_BEGIN_DECLS

#define NCM_TYPE_SF_SPHERICAL_HARMONICS             (ncm_sf_spherical_harmonics_get_type ())
#define NCM_SF_SPHERICAL_HARMONICS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SF_SPHERICAL_HARMONICS, NcmSFSphericalHarmonics))
#define NCM_SF_SPHERICAL_HARMONICS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SF_SPHERICAL_HARMONICS, NcmSFSphericalHarmonicsClass))
#define NCM_IS_SF_SPHERICAL_HARMONICS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SF_SPHERICAL_HARMONICS))
#define NCM_IS_SF_SPHERICAL_HARMONICS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SF_SPHERICAL_HARMONICS))
#define NCM_SF_SPHERICAL_HARMONICS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SF_SPHERICAL_HARMONICS, NcmSFSphericalHarmonicsClass))

typedef struct _NcmSFSphericalHarmonicsClass NcmSFSphericalHarmonicsClass;
typedef struct _NcmSFSphericalHarmonics NcmSFSphericalHarmonics;

struct _NcmSFSphericalHarmonicsClass
{
	/*< private >*/
	GObjectClass parent_class;
};

typedef struct _NcmSFSphericalHarmonicsK
{
	gdouble l;
	gdouble lp1;
} NcmSFSphericalHarmonicsK;

struct _NcmSFSphericalHarmonics
{
	/*< private >*/
	GObject parent_instance;
	guint lmax;
	guint l;
	guint l0;
	guint m;
	NcmSFSphericalHarmonicsK * restrict Klm;
	GArray *sqrt_n;
	GArray *sqrtm1_n;
	GPtrArray *K_array;
	gdouble x;
	gdouble sqrt1mx2;
	gdouble Pl0m;
	gdouble Pl0p1m;
	gdouble Plm;
	gdouble Plp1m;
	gdouble abstol;
};

GType ncm_sf_spherical_harmonics_get_type (void) G_GNUC_CONST;

NcmSFSphericalHarmonics *ncm_sf_spherical_harmonics_new (const guint lmax);
NcmSFSphericalHarmonics *ncm_sf_spherical_harmonics_ref (NcmSFSphericalHarmonics *spha);
void ncm_sf_spherical_harmonics_free (NcmSFSphericalHarmonics *spha);
void ncm_sf_spherical_harmonics_clear (NcmSFSphericalHarmonics **spha);

void ncm_sf_spherical_harmonics_set_lmax (NcmSFSphericalHarmonics *spha, const guint lmax);
guint ncm_sf_spherical_harmonics_get_lmax (NcmSFSphericalHarmonics *spha);

void ncm_sf_spherical_harmonics_set_abstol (NcmSFSphericalHarmonics *spha, const gdouble abstol);
gdouble ncm_sf_spherical_harmonics_get_abstol (NcmSFSphericalHarmonics *spha);

G_INLINE_FUNC void ncm_sf_spherical_harmonics_start_rec (NcmSFSphericalHarmonics *spha, const gdouble x, const gdouble sqrt1mx2);
G_INLINE_FUNC void ncm_sf_spherical_harmonics_next_l (NcmSFSphericalHarmonics *spha);
G_INLINE_FUNC void ncm_sf_spherical_harmonics_next_m (NcmSFSphericalHarmonics *spha);

G_INLINE_FUNC gdouble ncm_sf_spherical_harmonics_get_Yblm (NcmSFSphericalHarmonics *spha);
G_INLINE_FUNC gdouble ncm_sf_spherical_harmonics_get_Yblp1m (NcmSFSphericalHarmonics *spha);
G_INLINE_FUNC gdouble ncm_sf_spherical_harmonics_get_x (NcmSFSphericalHarmonics *spha);
G_INLINE_FUNC gint ncm_sf_spherical_harmonics_get_l (NcmSFSphericalHarmonics *spha);
G_INLINE_FUNC gint ncm_sf_spherical_harmonics_get_m (NcmSFSphericalHarmonics *spha);

G_END_DECLS

#endif /* _NCM_SF_SPHERICAL_HARMONICS_H_ */

#ifndef _NCM_SF_SPHERICAL_HARMONICS_INLINE_H_
#define _NCM_SF_SPHERICAL_HARMONICS_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

#define SN(n)   g_array_index (spha->sqrt_n, gdouble, (n))
#define SNM1(n) g_array_index (spha->sqrtm1_n, gdouble, (n))

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_start_rec (NcmSFSphericalHarmonics *spha, const gdouble x, const gdouble sqrt1mx2)
{
	GArray *Km_array = g_ptr_array_index (spha->K_array, 0);
	
	spha->x        = x;
	spha->sqrt1mx2 = sqrt1mx2;
	spha->l     	 = 0;
	spha->l0    	 = 0;
	spha->m     	 = 0;
	spha->Klm      = &g_array_index (Km_array, NcmSFSphericalHarmonicsK, 0);
	spha->Pl0m     = ncm_c_sqrt_1_4pi ();
	spha->Pl0p1m   = x * SN (2 * spha->l + 3) * spha->Pl0m;
	spha->Plm      = spha->Pl0m;
	spha->Plp1m    = spha->Pl0p1m;
}

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_next_l (NcmSFSphericalHarmonics *spha)
{
	const gdouble x     = spha->x;
	const gdouble Plp2m = spha->Klm->lp1 * x * spha->Plp1m - spha->Klm->l * spha->Plm;

	spha->Plm           = spha->Plp1m;
	spha->Plp1m         = Plp2m;
	spha->l++;
	spha->Klm++;
}

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_next_m (NcmSFSphericalHarmonics *spha)
{
	const gdouble sqrt1mx2 = spha->sqrt1mx2;
	const gdouble x        = spha->x;

	if (spha->m == spha->l0)
	{
		const gint l0 = spha->m;

		spha->Pl0m    = - sqrt1mx2 * SN (2 * l0 + 3) * SNM1 (2 * l0 + 2) * spha->Pl0m;
		spha->Pl0p1m  = x * SN (2 * l0 + 5) * spha->Pl0m;

		spha->m++;
		spha->l0      = spha->m;
		
		spha->l       = spha->l0;
		
		spha->Plm     = spha->Pl0m;
		spha->Plp1m   = spha->Pl0p1m;
	}
	else
	{
		const gdouble sqrt1mx2 = spha->sqrt1mx2;
		const gdouble x        = spha->x;
		const gint l0          = spha->l0;
		const gint twol0       = 2 * l0;
		const gint m           = spha->m;
		const gint l0mm        = l0 - m;
		const gint l0pm        = l0 + m;
		const gdouble Pl0p1m   = spha->Pl0p1m;
		const gdouble Pl0m     = spha->Pl0m;
		
		spha->Pl0m  = (SN (twol0 + 1) * SN (l0mm + 1)   * SNM1 (twol0 + 3) * SNM1 (l0mm)     * Pl0p1m - SN (l0pm + 1)  * SNM1 (l0mm)                                        * x * Pl0m) / sqrt1mx2;
		spha->Plp1m = (SN (l0mm + 1)  * SNM1 (l0pm + 2)                                  * x * Pl0p1m - SN (twol0 + 3) * SN (l0pm + 1) * SNM1 (twol0 + 1) * SNM1 (l0pm + 2)     * Pl0m) / sqrt1mx2;
		spha->Plm   = spha->Pl0m;

		spha->m++;

		spha->l     = spha->l0;
	}

	{
		GArray *Km_array = g_ptr_array_index (spha->K_array, spha->m);
		spha->Klm = &g_array_index (Km_array, NcmSFSphericalHarmonicsK, spha->l - spha->m);
	}
	
	/*printf ("#(%6d, %6d)[% 22.15g]:", spha->l0, spha->m, fabs (spha->Plm));*/
	if (fabs (spha->Plm) < spha->abstol)
	{
		do {
		  ncm_sf_spherical_harmonics_next_l (spha);
		  /*printf (".");*/
	  } while (fabs (spha->Plm) < spha->abstol);

		spha->l0     = ncm_sf_spherical_harmonics_get_l (spha);
		spha->Pl0m   = spha->Plm;
		spha->Pl0p1m = spha->Plp1m;
	}
	/*printf ("\n");*/
}

G_INLINE_FUNC gdouble 
ncm_sf_spherical_harmonics_get_Yblm (NcmSFSphericalHarmonics *spha)
{
	return spha->Plm;
}

G_INLINE_FUNC gdouble 
ncm_sf_spherical_harmonics_get_Yblp1m (NcmSFSphericalHarmonics *spha)
{
	return spha->Plp1m;
}

G_INLINE_FUNC gdouble 
ncm_sf_spherical_harmonics_get_x (NcmSFSphericalHarmonics *spha)
{
	return spha->x;
}

G_INLINE_FUNC gint 
ncm_sf_spherical_harmonics_get_l (NcmSFSphericalHarmonics *spha)
{
	return spha->l;
}

G_INLINE_FUNC gint 
ncm_sf_spherical_harmonics_get_m (NcmSFSphericalHarmonics *spha)
{
	return spha->m;
}

#undef SN

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_SF_SPHERICAL_HARMONICS_INLINE_H_ */
