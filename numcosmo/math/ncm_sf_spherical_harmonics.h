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

struct _NcmSFSphericalHarmonics
{
	/*< private >*/
	GObject parent_instance;
	guint lmax;
	guint l;
	guint m;
	GArray *sqrt_n;
	gdouble x;
	gdouble sqrt1mx2;
	gdouble Pll;
	gdouble Plm;
	gdouble Plp1m;
};

GType ncm_sf_spherical_harmonics_get_type (void) G_GNUC_CONST;

NcmSFSphericalHarmonics *ncm_sf_spherical_harmonics_new (const guint lmax);
NcmSFSphericalHarmonics *ncm_sf_spherical_harmonics_ref (NcmSFSphericalHarmonics *spha);
void ncm_sf_spherical_harmonics_free (NcmSFSphericalHarmonics *spha);
void ncm_sf_spherical_harmonics_clear (NcmSFSphericalHarmonics **spha);

void ncm_sf_spherical_harmonics_set_lmax (NcmSFSphericalHarmonics *spha, const guint lmax);
guint ncm_sf_spherical_harmonics_get_lmax (NcmSFSphericalHarmonics *spha);

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

#define SN(n) g_array_index (spha->sqrt_n, gdouble, (n))

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_start_rec (NcmSFSphericalHarmonics *spha, const gdouble x, const gdouble sqrt1mx2)
{
	spha->x        = x;
	spha->sqrt1mx2 = sqrt1mx2;
	spha->l     	 = 0;
	spha->m     	 = 0;
	spha->Pll      = ncm_c_sqrt_1_4pi ();
	spha->Plm      = spha->Pll;
	spha->Plp1m    = x * SN (2 * spha->l + 3) * spha->Pll;
}

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_next_l (NcmSFSphericalHarmonics *spha)
{
	const gdouble x     = spha->x;
	const gint l        = spha->l;
	const gint twol     = 2 * l;
	const gint m        = spha->m;
	const gint lmm      = l - m;
	const gint lpm      = l + m;
	const gdouble Plp2m = SN (twol + 5) / (SN (lpm + 2) * SN (lmm + 2)) * ( x * SN (twol + 3) * spha->Plp1m - SN (lmm + 1) * SN (lpm + 1) * spha->Plm / SN (twol + 1));

	spha->Plm           = spha->Plp1m;
	spha->Plp1m         = Plp2m;
	spha->l++;
}

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_next_m (NcmSFSphericalHarmonics *spha)
{
	const gdouble sqrt1mx2 = spha->sqrt1mx2;
	const gdouble x        = spha->x;
	const gint l           = spha->m;

	spha->m++;

	spha->Pll   = - sqrt1mx2 * SN (2 * l + 3) / SN (2 * l + 2) * spha->Pll;
	spha->l     = spha->m;
	spha->Plm   = spha->Pll;

	spha->Plp1m = x * SN (2 * spha->l + 3) * spha->Pll;
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
