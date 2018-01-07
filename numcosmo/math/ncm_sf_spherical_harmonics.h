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
	GArray *sqrt_n;
	GArray *sqrtm1_n;
	GPtrArray *K_array;
	gdouble abstol;
};

typedef struct _NcmSFSphericalHarmonicsY
{
	gdouble x;
	gdouble sqrt1mx2;
	guint l;
	guint l0;
	guint m;
	NcmSFSphericalHarmonicsK * restrict Klm;
	gdouble Pl0m;
	gdouble Pl0p1m;
	gdouble Plm;
	gdouble Plp1m;
	NcmSFSphericalHarmonics *spha;
} NcmSFSphericalHarmonicsY;

typedef	struct _NcmSFSphericalHarmonicsP 
{
	gdouble x;
	gdouble sqrt1mx2;
	gdouble l0m;
	gdouble l0p1m;
	gdouble lm;
	gdouble lp1m;
} NcmSFSphericalHarmonicsP;

typedef struct _NcmSFSphericalHarmonicsYArray
{
	guint l;
	guint l0;
	guint m;
	guint len;
	NcmSFSphericalHarmonicsP * restrict P;
	NcmSFSphericalHarmonicsK * restrict Klm;
	NcmSFSphericalHarmonics *spha;
} NcmSFSphericalHarmonicsYArray;


GType ncm_sf_spherical_harmonics_Y_get_type (void) G_GNUC_CONST;
GType ncm_sf_spherical_harmonics_Y_array_get_type (void) G_GNUC_CONST;
GType ncm_sf_spherical_harmonics_get_type (void) G_GNUC_CONST;

NcmSFSphericalHarmonicsY *ncm_sf_spherical_harmonics_Y_new (NcmSFSphericalHarmonics *spha);
NcmSFSphericalHarmonicsY *ncm_sf_spherical_harmonics_Y_dup (NcmSFSphericalHarmonicsY *sphaY);
void ncm_sf_spherical_harmonics_Y_free (NcmSFSphericalHarmonicsY *sphaY);

G_INLINE_FUNC gdouble ncm_sf_spherical_harmonics_Y_get_lm (NcmSFSphericalHarmonicsY *sphaY);
G_INLINE_FUNC gdouble ncm_sf_spherical_harmonics_Y_get_lp1m (NcmSFSphericalHarmonicsY *sphaY);
G_INLINE_FUNC gdouble ncm_sf_spherical_harmonics_Y_get_x (NcmSFSphericalHarmonicsY *sphaY);
G_INLINE_FUNC gint ncm_sf_spherical_harmonics_Y_get_l (NcmSFSphericalHarmonicsY *sphaY);
G_INLINE_FUNC gint ncm_sf_spherical_harmonics_Y_get_m (NcmSFSphericalHarmonicsY *sphaY);

G_INLINE_FUNC void ncm_sf_spherical_harmonics_Y_next_l (NcmSFSphericalHarmonicsY *sphaY);
G_INLINE_FUNC void ncm_sf_spherical_harmonics_Y_next_l2 (NcmSFSphericalHarmonicsY *sphaY, gdouble Yblm[2]);
G_INLINE_FUNC void ncm_sf_spherical_harmonics_Y_next_l4 (NcmSFSphericalHarmonicsY *sphaY, gdouble Yblm[4]);
G_INLINE_FUNC void ncm_sf_spherical_harmonics_Y_next_l2pn (NcmSFSphericalHarmonicsY *sphaY, gdouble *Yblm, const gint n);
G_INLINE_FUNC void ncm_sf_spherical_harmonics_Y_next_m (NcmSFSphericalHarmonicsY *sphaY);

NcmSFSphericalHarmonicsYArray *ncm_sf_spherical_harmonics_Y_array_new (NcmSFSphericalHarmonics *spha, const guint len);
NcmSFSphericalHarmonicsYArray *ncm_sf_spherical_harmonics_Y_array_dup (NcmSFSphericalHarmonicsYArray *sphaYa);
void ncm_sf_spherical_harmonics_Y_array_free (NcmSFSphericalHarmonicsYArray *sphaYa);

G_INLINE_FUNC gdouble ncm_sf_spherical_harmonics_Y_array_get_lm (NcmSFSphericalHarmonicsYArray *sphaYa, const guint i);
G_INLINE_FUNC gdouble ncm_sf_spherical_harmonics_Y_array_get_lp1m (NcmSFSphericalHarmonicsYArray *sphaYa, const guint i);
G_INLINE_FUNC gdouble ncm_sf_spherical_harmonics_Y_array_get_x (NcmSFSphericalHarmonicsYArray *sphaYa, const guint i);
G_INLINE_FUNC gint ncm_sf_spherical_harmonics_Y_array_get_l (NcmSFSphericalHarmonicsYArray *sphaYa);
G_INLINE_FUNC gint ncm_sf_spherical_harmonics_Y_array_get_m (NcmSFSphericalHarmonicsYArray *sphaYa);

G_INLINE_FUNC void ncm_sf_spherical_harmonics_Y_array_next_l (NcmSFSphericalHarmonicsYArray *sphaYa, const guint len);
G_INLINE_FUNC void ncm_sf_spherical_harmonics_Y_array_next_l2 (NcmSFSphericalHarmonicsYArray *sphaYa, const guint len, gdouble *Yblm);
G_INLINE_FUNC void ncm_sf_spherical_harmonics_Y_array_next_l4 (NcmSFSphericalHarmonicsYArray *sphaYa, const guint len, gdouble *Yblm);
G_INLINE_FUNC void ncm_sf_spherical_harmonics_Y_array_next_l2pn (NcmSFSphericalHarmonicsYArray *sphaYa, const guint len, gdouble *Yblm, const gint n);
G_INLINE_FUNC void ncm_sf_spherical_harmonics_Y_array_next_m (NcmSFSphericalHarmonicsYArray *sphaYa, const guint len);

NcmSFSphericalHarmonics *ncm_sf_spherical_harmonics_new (const guint lmax);
NcmSFSphericalHarmonics *ncm_sf_spherical_harmonics_ref (NcmSFSphericalHarmonics *spha);
void ncm_sf_spherical_harmonics_free (NcmSFSphericalHarmonics *spha);
void ncm_sf_spherical_harmonics_clear (NcmSFSphericalHarmonics **spha);

void ncm_sf_spherical_harmonics_set_lmax (NcmSFSphericalHarmonics *spha, const guint lmax);
guint ncm_sf_spherical_harmonics_get_lmax (NcmSFSphericalHarmonics *spha);

void ncm_sf_spherical_harmonics_set_abstol (NcmSFSphericalHarmonics *spha, const gdouble abstol);
gdouble ncm_sf_spherical_harmonics_get_abstol (NcmSFSphericalHarmonics *spha);

G_INLINE_FUNC void ncm_sf_spherical_harmonics_start_rec (NcmSFSphericalHarmonics *spha, NcmSFSphericalHarmonicsY *sphaY, const gdouble theta);
G_INLINE_FUNC void ncm_sf_spherical_harmonics_start_rec_array (NcmSFSphericalHarmonics *spha, NcmSFSphericalHarmonicsYArray *sphaYa, const guint len, const gdouble *theta);

G_END_DECLS

#endif /* _NCM_SF_SPHERICAL_HARMONICS_H_ */

#ifndef _NCM_SF_SPHERICAL_HARMONICS_INLINE_H_
#define _NCM_SF_SPHERICAL_HARMONICS_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

#define SN(n)   g_array_index (spha->sqrt_n, gdouble, (n))
#define SNM1(n) g_array_index (spha->sqrtm1_n, gdouble, (n))

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_Y_next_l (NcmSFSphericalHarmonicsY *sphaY)
{
	const gdouble x     = sphaY->x;
	const gdouble Plp2m = sphaY->Klm->lp1 * x * sphaY->Plp1m - sphaY->Klm->l * sphaY->Plm;

	sphaY->Plm           = sphaY->Plp1m;
	sphaY->Plp1m         = Plp2m;
	sphaY->l++;
	sphaY->Klm++;
}

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_Y_next_l2 (NcmSFSphericalHarmonicsY *sphaY, gdouble Yblm[2])
{
	const gdouble x     = sphaY->x;

	Yblm[0]     = sphaY->Plm;
	Yblm[1]     = sphaY->Plp1m;
	sphaY->Plm   = sphaY->Klm[0].lp1 * x * Yblm[1]    - sphaY->Klm[0].l * Yblm[0];
	sphaY->Plp1m = sphaY->Klm[1].lp1 * x * sphaY->Plm - sphaY->Klm[1].l * Yblm[1];

	sphaY->l            += 2;
	sphaY->Klm          += 2;
}

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_Y_next_l4 (NcmSFSphericalHarmonicsY *sphaY, gdouble Yblm[4])
{
	const gdouble x = sphaY->x;

	Yblm[0]      = sphaY->Plm;
	Yblm[1]      = sphaY->Plp1m;
	Yblm[2]      = sphaY->Klm[0].lp1 * x * Yblm[1]    - sphaY->Klm[0].l * Yblm[0];
	Yblm[3]      = sphaY->Klm[1].lp1 * x * Yblm[2]    - sphaY->Klm[1].l * Yblm[1];
	sphaY->Plm   = sphaY->Klm[2].lp1 * x * Yblm[3]    - sphaY->Klm[2].l * Yblm[2];
	sphaY->Plp1m = sphaY->Klm[3].lp1 * x * sphaY->Plm - sphaY->Klm[3].l * Yblm[3];

	sphaY->l    += 4;
	sphaY->Klm  += 4;
}

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_Y_next_l2pn (NcmSFSphericalHarmonicsY *sphaY, gdouble *Yblm, const gint n)
{
	const gdouble x = sphaY->x;
	const gint np2  = n + 2;
	gint i;

	Yblm[0]     = sphaY->Plm;
	Yblm[1]     = sphaY->Plp1m;

	for (i = 0; i < n; i++)
	{
		const NcmSFSphericalHarmonicsK *Klm = &sphaY->Klm[i];
		Yblm[i + 2] = Klm->lp1 * x * Yblm[i + 1]   - Klm->l * Yblm[i + 0];
	}
	
	sphaY->Plm   = sphaY->Klm[i].lp1     * x * Yblm[i + 1]  - sphaY->Klm[i].l     * Yblm[i + 0];
	sphaY->Plp1m = sphaY->Klm[i + 1].lp1 * x * sphaY->Plm   - sphaY->Klm[i + 1].l * Yblm[i + 1];

	sphaY->l    += np2;
	sphaY->Klm  += np2;
}

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_Y_next_m (NcmSFSphericalHarmonicsY *sphaY)
{
	NcmSFSphericalHarmonics *spha = sphaY->spha;
	const gdouble sqrt1mx2        = sphaY->sqrt1mx2;
	const gdouble x               = sphaY->x;

	if (sphaY->m == sphaY->l0)
	{
		const gint l0  = sphaY->m;

		sphaY->Pl0m    = - sqrt1mx2 * SN (2 * l0 + 3) * SNM1 (2 * l0 + 2) * sphaY->Pl0m;
		sphaY->Pl0p1m  = x * SN (2 * l0 + 5) * sphaY->Pl0m;

		sphaY->m++;
		sphaY->l0      = sphaY->m;
		
		sphaY->l       = sphaY->l0;
		
		sphaY->Plm     = sphaY->Pl0m;
		sphaY->Plp1m   = sphaY->Pl0p1m;
	}
	else
	{
    const gdouble sqrt1mx2 = sphaY->sqrt1mx2;
	  const gdouble x        = sphaY->x;
		const gint l0          = sphaY->l0;
		const gint twol0       = 2 * l0;
		const gint m           = sphaY->m;
		const gint l0mm        = l0 - m;
		const gint l0pm        = l0 + m;
		const gdouble Pl0p1m   = sphaY->Pl0p1m;
		const gdouble Pl0m     = sphaY->Pl0m;
		
		sphaY->Pl0m  = (SN (twol0 + 1) * SN (l0mm + 1)   * SNM1 (twol0 + 3) * SNM1 (l0mm)     * Pl0p1m - SN (l0pm + 1)  * SNM1 (l0mm)                                        * x * Pl0m) / sqrt1mx2;
		sphaY->Plp1m = (SN (l0mm + 1)  * SNM1 (l0pm + 2)                                  * x * Pl0p1m - SN (twol0 + 3) * SN (l0pm + 1) * SNM1 (twol0 + 1) * SNM1 (l0pm + 2)     * Pl0m) / sqrt1mx2;
		sphaY->Plm   = sphaY->Pl0m;

		sphaY->m++;

		sphaY->l     = sphaY->l0;
	}

	{
		GArray *Km_array = g_ptr_array_index (sphaY->spha->K_array, sphaY->m);
		sphaY->Klm       = &g_array_index (Km_array, NcmSFSphericalHarmonicsK, sphaY->l - sphaY->m);
	}
	
	/*printf ("#(%6d, %6d)[% 22.15g]:", spha->l0, spha->m, fabs (spha->Plm));*/
	if (fabs (sphaY->Plm) < sphaY->spha->abstol)
	{
		const guint lmax = ncm_sf_spherical_harmonics_get_lmax (spha);
		guint l;
		
		do {
		  ncm_sf_spherical_harmonics_Y_next_l (sphaY);
			l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
		  /*printf (".");*/
	  } while ((fabs (sphaY->Plm) < sphaY->spha->abstol) && (l <= lmax));

		sphaY->l0     = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
		sphaY->Pl0m   = sphaY->Plm;
		sphaY->Pl0p1m = sphaY->Plp1m;
	}
	/*printf ("\n");*/
}

G_INLINE_FUNC gdouble 
ncm_sf_spherical_harmonics_Y_get_lm (NcmSFSphericalHarmonicsY *sphaY)
{
	return sphaY->Plm;
}

G_INLINE_FUNC gdouble 
ncm_sf_spherical_harmonics_Y_get_lp1m (NcmSFSphericalHarmonicsY *sphaY)
{
	return sphaY->Plp1m;
}

G_INLINE_FUNC gdouble 
ncm_sf_spherical_harmonics_Y_get_x (NcmSFSphericalHarmonicsY *sphaY)
{
	return sphaY->x;
}

G_INLINE_FUNC gint 
ncm_sf_spherical_harmonics_Y_get_l (NcmSFSphericalHarmonicsY *sphaY)
{
	return sphaY->l;
}

G_INLINE_FUNC gint 
ncm_sf_spherical_harmonics_Y_get_m (NcmSFSphericalHarmonicsY *sphaY)
{
	return sphaY->m;
}

/* Array methods */

G_INLINE_FUNC gdouble 
ncm_sf_spherical_harmonics_Y_array_get_lm (NcmSFSphericalHarmonicsYArray *sphaYa, const guint i)
{
	return sphaYa->P[i].lm;
}

G_INLINE_FUNC gdouble 
ncm_sf_spherical_harmonics_Y_array_get_lp1m (NcmSFSphericalHarmonicsYArray *sphaYa, const guint i)
{
	return sphaYa->P[i].lp1m;
}

G_INLINE_FUNC gdouble 
ncm_sf_spherical_harmonics_Y_array_get_x (NcmSFSphericalHarmonicsYArray *sphaYa, const guint i)
{
	return sphaYa->P[i].x;
}

G_INLINE_FUNC gint 
ncm_sf_spherical_harmonics_Y_array_get_l (NcmSFSphericalHarmonicsYArray *sphaYa)
{
	return sphaYa->l;
}

G_INLINE_FUNC gint 
ncm_sf_spherical_harmonics_Y_array_get_m (NcmSFSphericalHarmonicsYArray *sphaYa)
{
	return sphaYa->m;
}

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_Y_array_next_l (NcmSFSphericalHarmonicsYArray *sphaYa, const guint len)
{
	guint i;
	for (i = 0; i < len; i++)
	{
		const gdouble x     = sphaYa->P[i].x;
		const gdouble Plp2m = sphaYa->Klm->lp1 * x * sphaYa->P[i].lp1m - sphaYa->Klm->l * sphaYa->P[i].lm;

		sphaYa->P[i].lm           = sphaYa->P[i].lp1m;
		sphaYa->P[i].lp1m         = Plp2m;
	}

	sphaYa->l++;
	sphaYa->Klm++;
}

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_Y_array_next_l2 (NcmSFSphericalHarmonicsYArray *sphaYa, const guint len, gdouble *Yblm)
{
	guint i;
	for (i = 0; i < len; i++)
	{
		const gdouble x = sphaYa->P[i].x;

		Yblm[len * 0 + i]   = sphaYa->P[i].lm;
		Yblm[len * 1 + i]   = sphaYa->P[i].lp1m;

		sphaYa->P[i].lm   = sphaYa->Klm[0].lp1 * x * Yblm[len * 1 + i] - sphaYa->Klm[0].l * Yblm[len * 0 + i];
		sphaYa->P[i].lp1m = sphaYa->Klm[1].lp1 * x * sphaYa->P[i].lm   - sphaYa->Klm[1].l * Yblm[len * 1 + i];
	}

	sphaYa->l            += 2;
	sphaYa->Klm          += 2;
}

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_Y_array_next_l4 (NcmSFSphericalHarmonicsYArray *sphaYa, const guint len, gdouble *Yblm)
{
	guint i;
	for (i = 0; i < len; i++)
	{	
		const gdouble x = sphaYa->P[i].x;

		Yblm[len * 0 + i] = sphaYa->P[i].lm;
		Yblm[len * 1 + i] = sphaYa->P[i].lp1m;
		Yblm[len * 2 + i] = sphaYa->Klm[0].lp1 * x * Yblm[len * 1 + i] - sphaYa->Klm[0].l * Yblm[len * 0 + i];
		Yblm[len * 3 + i] = sphaYa->Klm[1].lp1 * x * Yblm[len * 2 + i] - sphaYa->Klm[1].l * Yblm[len * 1 + i];
		sphaYa->P[i].lm   = sphaYa->Klm[2].lp1 * x * Yblm[len * 3 + i] - sphaYa->Klm[2].l * Yblm[len * 2 + i];
		sphaYa->P[i].lp1m = sphaYa->Klm[3].lp1 * x * sphaYa->P[i].lm   - sphaYa->Klm[3].l * Yblm[len * 3 + i];
	}
	
	sphaYa->l    += 4;
	sphaYa->Klm  += 4;
}

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_Y_array_next_l2pn (NcmSFSphericalHarmonicsYArray *sphaYa, const guint len, gdouble *Yblm, const gint n)
{
	const gint np2  = n + 2;
	guint i, j;

	for (i = 0; i < len; i++)
	{	
		Yblm[len * 0 + i] = sphaYa->P[i].lm;
		Yblm[len * 1 + i] = sphaYa->P[i].lp1m;
	}
	
	for (j = 0; j < n; j++)
	{
		const NcmSFSphericalHarmonicsK *Klm = &sphaYa->Klm[j];
		for (i = 0; i < len; i++)
		{	
			const gdouble x = sphaYa->P[i].x;
			Yblm[len * (j + 2) + i] = Klm->lp1 * x * Yblm[len * (j + 1) + i] - Klm->l * Yblm[len * (j + 0) + i];
		}
	}

	for (i = 0; i < len; i++)
	{	
		const gdouble x = sphaYa->P[i].x;
		sphaYa->P[i].lm   = sphaYa->Klm[j].lp1     * x * Yblm[len * (j + 1) + i] - sphaYa->Klm[j].l     * Yblm[len * (j + 0) + i];
		sphaYa->P[i].lp1m = sphaYa->Klm[j + 1].lp1 * x * sphaYa->P[i].lm         - sphaYa->Klm[j + 1].l * Yblm[len * (j + 1) + i];
	}

	sphaYa->l    += np2;
	sphaYa->Klm  += np2;
}

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_Y_array_next_m (NcmSFSphericalHarmonicsYArray *sphaYa, const guint len)
{
	NcmSFSphericalHarmonics *spha = sphaYa->spha;
	guint i;

	if (sphaYa->m == sphaYa->l0)
	{
		const gint l0            = sphaYa->m;
		const gdouble sn_2l0_3   = SN (2 * l0 + 3);
		const gdouble sn_2l0_5   = SN (2 * l0 + 5);
		const gdouble snm1_2l0_2 = SNM1 (2 * l0 + 2);
		
		for (i = 0; i < len; i++)
		{
			const gdouble sqrt1mx2 = sphaYa->P[i].sqrt1mx2;
			const gdouble x        = sphaYa->P[i].x;

			sphaYa->P[i].l0m   = - sqrt1mx2 * sn_2l0_3 * snm1_2l0_2 * sphaYa->P[i].l0m;
			sphaYa->P[i].l0p1m = x * sn_2l0_5 * sphaYa->P[i].l0m;

			sphaYa->P[i].lm     = sphaYa->P[i].l0m;
			sphaYa->P[i].lp1m   = sphaYa->P[i].l0p1m;
		}

		sphaYa->m++;
		sphaYa->l0 = sphaYa->m;
		sphaYa->l  = sphaYa->l0;
	}
	else
	{
		const gint l0    = sphaYa->l0;
		const gint m     = sphaYa->m;
		const gint twol0 = 2 * l0;
		const gint l0mm  = l0 - m;
		const gint l0pm  = l0 + m;
		
		for (i = 0; i < len; i++)
		{
			const gdouble sqrt1mx2 = sphaYa->P[i].sqrt1mx2;
			const gdouble x        = sphaYa->P[i].x;
			const gdouble Pl0p1m   = sphaYa->P[i].l0p1m;
			const gdouble Pl0m     = sphaYa->P[i].l0m;

			sphaYa->P[i].l0m  = (SN (twol0 + 1) * SN (l0mm + 1)   * SNM1 (twol0 + 3) * SNM1 (l0mm)     * Pl0p1m - SN (l0pm + 1)  * SNM1 (l0mm)                                        * x * Pl0m) / sqrt1mx2;
			sphaYa->P[i].lp1m = (SN (l0mm + 1)  * SNM1 (l0pm + 2)                                  * x * Pl0p1m - SN (twol0 + 3) * SN (l0pm + 1) * SNM1 (twol0 + 1) * SNM1 (l0pm + 2)     * Pl0m) / sqrt1mx2;
			sphaYa->P[i].lm   = sphaYa->P[i].l0m;
		}

		sphaYa->m++;
		sphaYa->l = sphaYa->l0;
	}

	{
		GArray *Km_array = g_ptr_array_index (sphaYa->spha->K_array, sphaYa->m);
		sphaYa->Klm      = &g_array_index (Km_array, NcmSFSphericalHarmonicsK, sphaYa->l - sphaYa->m);
	}

	{
		gdouble max_Plm = 0.0;

		for (i = 0; i < len; i++)
		{
			const gdouble abs_Plm = fabs (sphaYa->P[i].lm);
			max_Plm = MAX (max_Plm, abs_Plm);
		}
		
		if (max_Plm < sphaYa->spha->abstol)
		{
			const guint lmax = ncm_sf_spherical_harmonics_get_lmax (spha);
			guint l;
						
			do {
			  ncm_sf_spherical_harmonics_Y_array_next_l (sphaYa, len);
				l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
				
			  max_Plm = 0.0;
			  
			  for (i = 0; i < len; i++)
			  {
				  const gdouble abs_Plm = fabs (sphaYa->P[i].lm);
				  max_Plm = MAX (max_Plm, abs_Plm);
			  }
		  } while ((max_Plm < sphaYa->spha->abstol) && (l <= lmax));
			
			sphaYa->l0     = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);

			for (i = 0; i < len; i++)
			{
				sphaYa->P[i].l0m   = sphaYa->P[i].lm;
				sphaYa->P[i].l0p1m = sphaYa->P[i].lp1m;
			}
		}
	}
}

/* Start methods */

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_start_rec (NcmSFSphericalHarmonics *spha, NcmSFSphericalHarmonicsY *sphaY, const gdouble theta)
{
	GArray *Km_array = g_ptr_array_index (spha->K_array, 0);

	sincos (theta, &sphaY->sqrt1mx2, &sphaY->x);
	
	sphaY->l        = 0;
	sphaY->l0       = 0;
	sphaY->m        = 0;
	sphaY->Klm      = &g_array_index (Km_array, NcmSFSphericalHarmonicsK, 0);
	sphaY->Pl0m     = ncm_c_sqrt_1_4pi ();
	sphaY->Pl0p1m   = sphaY->x * SN (2 * sphaY->l + 3) * sphaY->Pl0m;
	sphaY->Plm      = sphaY->Pl0m;
	sphaY->Plp1m    = sphaY->Pl0p1m;
}

G_INLINE_FUNC void 
ncm_sf_spherical_harmonics_start_rec_array (NcmSFSphericalHarmonics *spha, NcmSFSphericalHarmonicsYArray *sphaYa, const guint len, const gdouble *theta)
{
	GArray *Km_array = g_ptr_array_index (spha->K_array, 0);
	guint i;

	g_assert_cmpuint (len, ==, sphaYa->len);
	
	sphaYa->Klm = &g_array_index (Km_array, NcmSFSphericalHarmonicsK, 0);
	sphaYa->l   = 0;
	sphaYa->l0  = 0;
	sphaYa->m   = 0;

	for (i = 0; i < len; i++)
	{
		sincos (theta[i], &sphaYa->P[i].sqrt1mx2, &sphaYa->P[i].x);

		sphaYa->P[i].l0m     = ncm_c_sqrt_1_4pi ();
		sphaYa->P[i].l0p1m   = sphaYa->P[i].x * SN (2 * sphaYa->l + 3) * sphaYa->P[i].l0m;
		sphaYa->P[i].lm      = sphaYa->P[i].l0m;
		sphaYa->P[i].lp1m    = sphaYa->P[i].l0p1m;
	}
}

#undef SN
#undef SNM1

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_SF_SPHERICAL_HARMONICS_INLINE_H_ */
