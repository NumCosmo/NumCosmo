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

#include <gsl/gsl_sf_legendre.h>

G_BEGIN_DECLS

#define NCM_TYPE_SF_SPHERICAL_HARMONICS             (ncm_sf_spherical_harmonics_get_type ())
#define NCM_SF_SPHERICAL_HARMONICS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SF_SPHERICAL_HARMONICS, NcmSFSphericalHarmonics))
#define NCM_SF_SPHERICAL_HARMONICS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SF_SPHERICAL_HARMONICS, NcmSFSphericalHarmonicsClass))
#define NCM_IS_SF_SPHERICAL_HARMONICS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SF_SPHERICAL_HARMONICS))
#define NCM_IS_SF_SPHERICAL_HARMONICS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SF_SPHERICAL_HARMONICS))
#define NCM_SF_SPHERICAL_HARMONICS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SF_SPHERICAL_HARMONICS, NcmSFSphericalHarmonicsClass))

typedef struct _NcmSFSphericalHarmonicsClass NcmSFSphericalHarmonicsClass;
typedef struct _NcmSFSphericalHarmonics NcmSFSphericalHarmonics;
typedef struct _NcmSFSphericalHarmonicsK NcmSFSphericalHarmonicsK;
typedef struct _NcmSFSphericalHarmonicsY NcmSFSphericalHarmonicsY;
typedef	struct _NcmSFSphericalHarmonicsP NcmSFSphericalHarmonicsP;
typedef struct _NcmSFSphericalHarmonicsYArray NcmSFSphericalHarmonicsYArray;

struct _NcmSFSphericalHarmonicsClass
{
	/*< private >*/
	GObjectClass parent_class;
};

/**
 * NcmSFSphericalHarmonicsK:
 * @l: $K_l$
 * @lp1: $K_{l+1}$
 * 
 * Recurrence coefficients.
 * 
 */
struct _NcmSFSphericalHarmonicsK
{
	gdouble l;
	gdouble lp1;
};

struct _NcmSFSphericalHarmonics
{
	/*< private >*/
	GObject parent_instance;
	gint lmax;
	GArray *sqrt_n;
	GArray *sqrtm1_n;
	GPtrArray *K_array;
	gint Klm_m;
};

/**
 * NcmSFSphericalHarmonicsY:
 * @x: $x$
 * @sqrt1mx2: $\sqrt{1-x^2}$ 
 * @l: $l$
 * @l0: $l_0$
 * @m: $m$
 * @Klm: #NcmSFSphericalHarmonicsK pointer
 * @Pl0m: $P_{l_0}^m$
 * @Pl0p1m: $P_{l_0+1}^m$
 * @Plm: $P_{l}^m$
 * @Plp1m: $P_{l+1}^m$
 * @spha: pointer to parent #NcmSFSphericalHarmonics
 * @abstol: absolute tolerance
 * 
 * Recurrence boxed object.
 * 
 */
struct _NcmSFSphericalHarmonicsY
{
	gdouble x;
	gdouble sqrt1mx2;
	gint l;
	gint l0;
	gint m;
	NcmSFSphericalHarmonicsK * restrict Klm;
	gdouble Pl0m;
	gdouble Pl0p1m;
	gdouble Plm;
	gdouble Plp1m;
	NcmSFSphericalHarmonics *spha;
	gdouble abstol;
};

/**
 * NcmSFSphericalHarmonicsP:
 * @x: $x$
 * @sqrt1mx2: $\sqrt{1-x^2}$ 
 * @l0m: $P_{l_0}^m$
 * @l0p1m: $P_{l_0+1}^m$
 * @lm: $P_{l}^m$
 * @lp1m: $P_{l+1}^m$
 * 
 * Boxed P values.
 * 
 */
struct _NcmSFSphericalHarmonicsP 
{
	gdouble x;
	gdouble sqrt1mx2;
	gdouble l0m;
	gdouble l0p1m;
	gdouble lm;
	gdouble lp1m;
};

#define NCM_SF_SPHERICAL_HARMONICS_MAX_LEN 6

/**
 * NcmSFSphericalHarmonicsYArray:
 * @l: $l$
 * @l0: $l_0$
 * @m: $m$
 * @x: array of $x$
 * @sqrt1mx2: array of $\sqrt{1-x^2}$
 * @Yl0m: array of $Y_{l_0}^m$
 * @Ylm: array of $Y_{l}^m$
 * @Klm: #NcmSFSphericalHarmonicsK pointer
 * @spha: pointer to parent #NcmSFSphericalHarmonics
 * @abstol: absolute tolerance
 * 
 * Recurrence array boxed object.
 * 
 */
struct _NcmSFSphericalHarmonicsYArray
{
	gint l;
	gint l0;
	gint m;
	gint len;
	gdouble x[NCM_SF_SPHERICAL_HARMONICS_MAX_LEN];
	gdouble sqrt1mx2[NCM_SF_SPHERICAL_HARMONICS_MAX_LEN];
	gdouble Yl0m[NCM_SF_SPHERICAL_HARMONICS_MAX_LEN * 2];
	gdouble Ylm[NCM_SF_SPHERICAL_HARMONICS_MAX_LEN * 2];
	NcmSFSphericalHarmonicsK * restrict Klm;
	NcmSFSphericalHarmonics *spha;
	gdouble abstol;
};

GType ncm_sf_spherical_harmonics_Y_get_type (void) G_GNUC_CONST;
GType ncm_sf_spherical_harmonics_Y_array_get_type (void) G_GNUC_CONST;
GType ncm_sf_spherical_harmonics_get_type (void) G_GNUC_CONST;

NcmSFSphericalHarmonicsY *ncm_sf_spherical_harmonics_Y_new (NcmSFSphericalHarmonics *spha, const gdouble abstol);
NcmSFSphericalHarmonicsY *ncm_sf_spherical_harmonics_Y_dup (NcmSFSphericalHarmonicsY *sphaY);
void ncm_sf_spherical_harmonics_Y_free (NcmSFSphericalHarmonicsY *sphaY);

NCM_INLINE gdouble ncm_sf_spherical_harmonics_Y_get_lm (NcmSFSphericalHarmonicsY *sphaY);
NCM_INLINE gdouble ncm_sf_spherical_harmonics_Y_get_lp1m (NcmSFSphericalHarmonicsY *sphaY);
NCM_INLINE gdouble ncm_sf_spherical_harmonics_Y_get_x (NcmSFSphericalHarmonicsY *sphaY);
NCM_INLINE gint ncm_sf_spherical_harmonics_Y_get_l (NcmSFSphericalHarmonicsY *sphaY);
NCM_INLINE gint ncm_sf_spherical_harmonics_Y_get_m (NcmSFSphericalHarmonicsY *sphaY);

NCM_INLINE void ncm_sf_spherical_harmonics_Y_next_l (NcmSFSphericalHarmonicsY *sphaY);
NCM_INLINE void ncm_sf_spherical_harmonics_Y_next_l2 (NcmSFSphericalHarmonicsY *sphaY, gdouble * restrict Yblm);
NCM_INLINE void ncm_sf_spherical_harmonics_Y_next_l4 (NcmSFSphericalHarmonicsY *sphaY, gdouble * restrict Yblm);
NCM_INLINE void ncm_sf_spherical_harmonics_Y_next_l2pn (NcmSFSphericalHarmonicsY *sphaY, gdouble * restrict Yblm, const gint n);
NCM_INLINE void ncm_sf_spherical_harmonics_Y_next_m (NcmSFSphericalHarmonicsY *sphaY);

NCM_INLINE void ncm_sf_spherical_harmonics_Y_reset (NcmSFSphericalHarmonicsY *sphaY);

NcmSFSphericalHarmonicsYArray *ncm_sf_spherical_harmonics_Y_array_new (NcmSFSphericalHarmonics *spha, const gint len, const gdouble abstol);
NcmSFSphericalHarmonicsYArray *ncm_sf_spherical_harmonics_Y_array_dup (NcmSFSphericalHarmonicsYArray *sphaYa);
void ncm_sf_spherical_harmonics_Y_array_free (NcmSFSphericalHarmonicsYArray *sphaYa);

NCM_INLINE gdouble ncm_sf_spherical_harmonics_Y_array_get_lm (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len, const gint i);
NCM_INLINE gdouble ncm_sf_spherical_harmonics_Y_array_get_lp1m (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len, const gint i);
NCM_INLINE gdouble ncm_sf_spherical_harmonics_Y_array_get_x (NcmSFSphericalHarmonicsYArray *sphaYa, const gint i);
NCM_INLINE gint ncm_sf_spherical_harmonics_Y_array_get_l (NcmSFSphericalHarmonicsYArray *sphaYa);
NCM_INLINE gint ncm_sf_spherical_harmonics_Y_array_get_m (NcmSFSphericalHarmonicsYArray *sphaYa);

NCM_INLINE void ncm_sf_spherical_harmonics_Y_array_next_l (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len);
NCM_INLINE void ncm_sf_spherical_harmonics_Y_array_next_l2 (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len, gdouble * restrict Yblm);
NCM_INLINE void ncm_sf_spherical_harmonics_Y_array_next_l4 (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len, gdouble * restrict Yblm);
NCM_INLINE void ncm_sf_spherical_harmonics_Y_array_next_l2pn (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len, gdouble * restrict Yblm, const gint n);
NCM_INLINE void ncm_sf_spherical_harmonics_Y_array_next_m (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len);

NCM_INLINE void ncm_sf_spherical_harmonics_Y_array_reset (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len);

NcmSFSphericalHarmonics *ncm_sf_spherical_harmonics_new (const gint lmax);
NcmSFSphericalHarmonics *ncm_sf_spherical_harmonics_ref (NcmSFSphericalHarmonics *spha);
void ncm_sf_spherical_harmonics_free (NcmSFSphericalHarmonics *spha);
void ncm_sf_spherical_harmonics_clear (NcmSFSphericalHarmonics **spha);

void ncm_sf_spherical_harmonics_set_lmax (NcmSFSphericalHarmonics *spha, const gint lmax);
gint ncm_sf_spherical_harmonics_get_lmax (NcmSFSphericalHarmonics *spha);

NCM_INLINE void ncm_sf_spherical_harmonics_start_rec (NcmSFSphericalHarmonics *spha, NcmSFSphericalHarmonicsY *sphaY, const gdouble theta);
NCM_INLINE void ncm_sf_spherical_harmonics_start_rec_array (NcmSFSphericalHarmonics *spha, NcmSFSphericalHarmonicsYArray *sphaYa, const gint len, const gdouble *theta);
NCM_INLINE NcmSFSphericalHarmonicsK *ncm_sf_spherical_harmonics_get_Klm (NcmSFSphericalHarmonics *spha, const gint l0, const gint m);

#define NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL (1.0e-20)
#define NCM_SF_SPHERICAL_HARMONICS_ARRAY_DEFAULT_ABSTOL (1.0e-40)
#define NCM_SF_SPHERICAL_HARMONICS_EPS (1.0e-280)
#define NCM_SF_SPHERICAL_HARMONICS_LATERAL_MOVE 1
#define NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX(ai,li,len) ((li) * (len) + (ai))

G_END_DECLS

#endif /* _NCM_SF_SPHERICAL_HARMONICS_H_ */

#ifndef _NCM_SF_SPHERICAL_HARMONICS_INLINE_H_
#define _NCM_SF_SPHERICAL_HARMONICS_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

#define _SN(n)   g_array_index (spha->sqrt_n, gdouble, (n))
#define _SNM1(n) g_array_index (spha->sqrtm1_n, gdouble, (n))

NCM_INLINE void 
ncm_sf_spherical_harmonics_Y_next_l (NcmSFSphericalHarmonicsY *sphaY)
{
	const gdouble x     = sphaY->x;
	const gdouble Plp2m = sphaY->Klm->lp1 * x * sphaY->Plp1m - sphaY->Klm->l * sphaY->Plm;

	sphaY->Plm           = sphaY->Plp1m;
	sphaY->Plp1m         = Plp2m;
	sphaY->l++;
	sphaY->Klm++;
}

NCM_INLINE void 
ncm_sf_spherical_harmonics_Y_next_l2 (NcmSFSphericalHarmonicsY *sphaY, gdouble * restrict Yblm)
{
	const gdouble x     = sphaY->x;

	Yblm[0]      = sphaY->Plm;
	Yblm[1]      = sphaY->Plp1m;
	sphaY->Plm   = sphaY->Klm[0].lp1 * x * Yblm[1]    - sphaY->Klm[0].l * Yblm[0];
	sphaY->Plp1m = sphaY->Klm[1].lp1 * x * sphaY->Plm - sphaY->Klm[1].l * Yblm[1];

	sphaY->l            += 2;
	sphaY->Klm          += 2;
}

NCM_INLINE void 
ncm_sf_spherical_harmonics_Y_next_l4 (NcmSFSphericalHarmonicsY *sphaY, gdouble * restrict Yblm)
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

NCM_INLINE void 
ncm_sf_spherical_harmonics_Y_next_l2pn (NcmSFSphericalHarmonicsY *sphaY, gdouble * restrict Yblm, const gint n)
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

NCM_INLINE void 
ncm_sf_spherical_harmonics_Y_next_m (NcmSFSphericalHarmonicsY *sphaY)
{
	NcmSFSphericalHarmonics *spha = sphaY->spha;
	const gdouble sqrt1mx2        = sphaY->sqrt1mx2;
	const gdouble x               = sphaY->x;

	if (sphaY->m == sphaY->l0)
	{
		const gint l0   = sphaY->m;

		sphaY->Pl0m     = - sqrt1mx2 * _SN (2 * l0 + 3) * _SNM1 (2 * l0 + 2) * sphaY->Pl0m;
		sphaY->Pl0p1m   = x * _SN (2 * l0 + 5) * sphaY->Pl0m;

		sphaY->m++;
		sphaY->l0       = sphaY->m;
		
		sphaY->l        = sphaY->l0;
		
		sphaY->Plm      = sphaY->Pl0m   * NCM_SF_SPHERICAL_HARMONICS_EPS;
		sphaY->Plp1m    = sphaY->Pl0p1m * NCM_SF_SPHERICAL_HARMONICS_EPS;
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
		const gdouble Pl0m     = sphaY->Pl0m;
		const gdouble Pl0p1m   = sphaY->Pl0p1m;

		const gdouble Llp1     = _SN (twol0 + 1) * _SN (l0mm + 1)   * _SNM1 (twol0 + 3) * _SNM1 (l0mm);
		const gdouble Ll       = _SN (l0pm + 1)  * _SNM1 (l0mm);
		const gdouble Mlp1     = _SN (l0mm + 1)  * _SNM1 (l0pm + 2);
		const gdouble Ml       = _SN (twol0 + 3) * _SN (l0pm + 1) * _SNM1 (twol0 + 1) * _SNM1 (l0pm + 2);

		sphaY->Pl0m   = (Llp1     * Pl0p1m - Ll * x * Pl0m) / sqrt1mx2;
		sphaY->Pl0p1m = (Mlp1 * x * Pl0p1m - Ml     * Pl0m) / sqrt1mx2;

		sphaY->Plm    = sphaY->Pl0m   * NCM_SF_SPHERICAL_HARMONICS_EPS;
		sphaY->Plp1m  = sphaY->Pl0p1m * NCM_SF_SPHERICAL_HARMONICS_EPS;

		sphaY->m++;

		sphaY->l     = sphaY->l0;
	}

	sphaY->Klm = ncm_sf_spherical_harmonics_get_Klm (spha, sphaY->l, sphaY->m);
	
	/*printf ("#(%6d, %6d)[% 22.15g]:", spha->l0, spha->m, fabs (spha->Plm));*/
	if (fabs (sphaY->Plm) < sphaY->abstol)
	{
		const gint lmax = ncm_sf_spherical_harmonics_get_lmax (spha);
		gdouble Pl0p1m = sphaY->Pl0p1m;
		gdouble Pl0m   = sphaY->Pl0m;
		
		do {
			const gdouble Pl0p2m = sphaY->Klm->lp1 * x * Pl0p1m - sphaY->Klm->l * Pl0m;

			Pl0m   = Pl0p1m;
			Pl0p1m = Pl0p2m;

			sphaY->l++;
			sphaY->Klm++;

			/*printf (".");*/
	  } while ((fabs (Pl0m * NCM_SF_SPHERICAL_HARMONICS_EPS) < sphaY->abstol) && (sphaY->l <= lmax));

		sphaY->Plm   = Pl0m * NCM_SF_SPHERICAL_HARMONICS_EPS;
		sphaY->Plp1m = Pl0p1m * NCM_SF_SPHERICAL_HARMONICS_EPS;
#ifdef NCM_SF_SPHERICAL_HARMONICS_LATERAL_MOVE
		sphaY->l0     = sphaY->l;
		sphaY->Pl0m   = Pl0m;
		sphaY->Pl0p1m = Pl0p1m;
#endif /* NCM_SF_SPHERICAL_HARMONICS_LATERAL_MOVE */
	}
	/*printf ("\n");*/
}

NCM_INLINE gdouble 
ncm_sf_spherical_harmonics_Y_get_lm (NcmSFSphericalHarmonicsY *sphaY)
{
	return sphaY->Plm;
}

NCM_INLINE gdouble 
ncm_sf_spherical_harmonics_Y_get_lp1m (NcmSFSphericalHarmonicsY *sphaY)
{
	return sphaY->Plp1m;
}

NCM_INLINE gdouble 
ncm_sf_spherical_harmonics_Y_get_x (NcmSFSphericalHarmonicsY *sphaY)
{
	return sphaY->x;
}

NCM_INLINE gint 
ncm_sf_spherical_harmonics_Y_get_l (NcmSFSphericalHarmonicsY *sphaY)
{
	return sphaY->l;
}

NCM_INLINE gint 
ncm_sf_spherical_harmonics_Y_get_m (NcmSFSphericalHarmonicsY *sphaY)
{
	return sphaY->m;
}

/* Array methods */

NCM_INLINE gdouble 
ncm_sf_spherical_harmonics_Y_array_get_lm (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len, const gint i)
{
	return sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)];
}

NCM_INLINE gdouble 
ncm_sf_spherical_harmonics_Y_array_get_lp1m (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len, const gint i)
{
	return sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)];
}

NCM_INLINE gdouble 
ncm_sf_spherical_harmonics_Y_array_get_x (NcmSFSphericalHarmonicsYArray *sphaYa, const gint i)
{
	return sphaYa->x[i];
}

NCM_INLINE gint 
ncm_sf_spherical_harmonics_Y_array_get_l (NcmSFSphericalHarmonicsYArray *sphaYa)
{
	return sphaYa->l;
}

NCM_INLINE gint 
ncm_sf_spherical_harmonics_Y_array_get_m (NcmSFSphericalHarmonicsYArray *sphaYa)
{
	return sphaYa->m;
}

NCM_INLINE void 
ncm_sf_spherical_harmonics_Y_array_next_l (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len)
{
	gint i;
	for (i = 0; i < len; i++)
	{
		const gdouble x     = sphaYa->x[i];
		const gdouble Plp2m = sphaYa->Klm->lp1 * x * sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] - sphaYa->Klm->l * sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)];

		sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] = sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)];
		sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] = Plp2m;
	}

	sphaYa->l++;
	sphaYa->Klm++;
}

NCM_INLINE void 
ncm_sf_spherical_harmonics_Y_array_next_l2 (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len, gdouble * restrict Yblm)
{
	gint i;
	for (i = 0; i < len; i++)
	{
		const gdouble x = sphaYa->x[i];
		
		Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)]        = sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)];
		Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)]        = sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)];

		sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] = sphaYa->Klm[0].lp1 * x * Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)]        - sphaYa->Klm[0].l * Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)];
		sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] = sphaYa->Klm[1].lp1 * x * sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] - sphaYa->Klm[1].l * Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)];
	}

	sphaYa->l            += 2;
	sphaYa->Klm          += 2;
}

NCM_INLINE void 
ncm_sf_spherical_harmonics_Y_array_next_l4 (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len, gdouble * restrict Yblm)
{
	NcmSFSphericalHarmonicsK * restrict Klm = sphaYa->Klm;
	gint i;

	for (i = 0; i < len; i++)
	{	
		const gdouble x = sphaYa->x[i];

		Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)]        = sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)];
		Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)]        = sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)];
		Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 2, len)]        = Klm[0].lp1 * x * Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)]        - Klm[0].l * Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)];
		Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 3, len)]        = Klm[1].lp1 * x * Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 2, len)]        - Klm[1].l * Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)];
		sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] = Klm[2].lp1 * x * Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 3, len)]        - Klm[2].l * Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 2, len)];
		sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] = Klm[3].lp1 * x * sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] - Klm[3].l * Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 3, len)];
	}
	
	sphaYa->l    += 4;
	sphaYa->Klm  += 4;
}

NCM_INLINE void 
ncm_sf_spherical_harmonics_Y_array_next_l2pn (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len, gdouble * restrict Yblm, const gint n)
{
	NcmSFSphericalHarmonicsK * restrict Klm = sphaYa->Klm;
	const gint np2  = n + 2;
	gint i, j;
	
	for (i = 0; i < len; i++)
	{	
		Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] = sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)];
		Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] = sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)];
	}
	
	for (j = 0; j < n; j++)
	{
		const gint jp0        = (j + 0) * len;
		const gint jp1        = (j + 1) * len;
		const gint jp2        = (j + 2) * len;
		
		for (i = 0; i < len; i++)
		{	
			const gdouble x = sphaYa->x[i];			
			Yblm[jp2 + i] = Klm[j].lp1 * x * Yblm[jp1 + i] - Klm[j].l * Yblm[jp0 + i];
		}
	}

	for (i = 0; i < len; i++)
	{	
		const gdouble x = sphaYa->x[i];
		sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] = Klm[n + 0].lp1 * x * Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, j + 1, len)]    - Klm[n + 0].l * Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, j + 0, len)];
		sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] = Klm[n + 1].lp1 * x * sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] - Klm[n + 1].l * Yblm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, j + 1, len)];
	}

	sphaYa->l    += np2;
	sphaYa->Klm  += np2;
}

NCM_INLINE void 
ncm_sf_spherical_harmonics_Y_array_next_m (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len)
{
	NcmSFSphericalHarmonics *spha = sphaYa->spha;
	gint i;

	if (sphaYa->m == sphaYa->l0)
	{
		const gint l0            = sphaYa->m;
		const gdouble sn_2l0_3   = _SN (2 * l0 + 3);
		const gdouble sn_2l0_5   = _SN (2 * l0 + 5);
		const gdouble snm1_2l0_2 = _SNM1 (2 * l0 + 2);
		
		for (i = 0; i < len; i++)
		{
			const gdouble sqrt1mx2 = sphaYa->sqrt1mx2[i];
			const gdouble x        = sphaYa->x[i];

			sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] = - sqrt1mx2 * sn_2l0_3 * snm1_2l0_2 * sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)];
			sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] = x * sn_2l0_5 * sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)];

			sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)]  = sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)]   * NCM_SF_SPHERICAL_HARMONICS_EPS;
			sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)]  = sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] * NCM_SF_SPHERICAL_HARMONICS_EPS;
		}

		sphaYa->m++;
		sphaYa->l0 = sphaYa->m;
		sphaYa->l  = sphaYa->l0;
	}
	else
	{
		const gint l0      = sphaYa->l0;
		const gint m       = sphaYa->m;
		const gint twol0   = 2 * l0;
		const gint l0mm    = l0 - m;
		const gint l0pm    = l0 + m;
		const gdouble Llp1 = _SN (twol0 + 1) * _SN (l0mm + 1)   * _SNM1 (twol0 + 3) * _SNM1 (l0mm);
		const gdouble Ll   = _SN (l0pm + 1)  * _SNM1 (l0mm);
		const gdouble Mlp1 = _SN (l0mm + 1)  * _SNM1 (l0pm + 2);
		const gdouble Ml   = _SN (twol0 + 3) * _SN (l0pm + 1) * _SNM1 (twol0 + 1) * _SNM1 (l0pm + 2);
		
		for (i = 0; i < len; i++)
		{
			const gdouble sqrt1mx2 = sphaYa->sqrt1mx2[i];
			const gdouble x        = sphaYa->x[i];
			const gdouble Pl0m     = sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)];
			const gdouble Pl0p1m   = sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)];

			sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] = (Llp1     * Pl0p1m - Ll * x * Pl0m) / sqrt1mx2;
			sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] = (Mlp1 * x * Pl0p1m - Ml     * Pl0m) / sqrt1mx2;
/*
			printf ("[%d] <%6d %6d> % 22.15g % 22.15g % e == % 22.15g | % 22.15g % 22.15g % e == % 22.15g \n", i, l0, m, 
			        Llp1     * Pl0p1m, Ll * x * Pl0m, Llp1     * Pl0p1m /( Ll * x * Pl0m) - 1.0, sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)],
			        Mlp1 * x * Pl0p1m, Ml     * Pl0m, Mlp1 * x * Pl0p1m / (Ml     * Pl0m) - 1.0, sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)]);
*/			
			sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] = sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)]   * NCM_SF_SPHERICAL_HARMONICS_EPS;
			sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] = sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] * NCM_SF_SPHERICAL_HARMONICS_EPS;
		}

		sphaYa->m++;
		sphaYa->l = sphaYa->l0;
	}

	sphaYa->Klm = ncm_sf_spherical_harmonics_get_Klm (spha, sphaYa->l, sphaYa->m);

	{
		gdouble min_Plm = GSL_POSINF;

		for (i = 0; i < len; i++)
		{
			const gdouble abs_Plm = fabs (sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)]);
			min_Plm = MIN (min_Plm, abs_Plm);
		}
		
		if (min_Plm < sphaYa->abstol)
		{
			const gint lmax = ncm_sf_spherical_harmonics_get_lmax (spha);
			gint l;
						
			do {
			  ncm_sf_spherical_harmonics_Y_array_next_l (sphaYa, len);
				l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
				
			  min_Plm = GSL_POSINF;
			  
			  for (i = 0; i < len; i++)
			  {
				  const gdouble abs_Plm = fabs (sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)]);
				  min_Plm = MIN (min_Plm, abs_Plm);
			  }
		  } while ((min_Plm < sphaYa->abstol) && (l <= lmax));

#ifdef NCM_SF_SPHERICAL_HARMONICS_LATERAL_MOVE
			sphaYa->l0 = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
			for (i = 0; i < len; i++)
			{
				sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] = sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] / NCM_SF_SPHERICAL_HARMONICS_EPS;
				sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] = sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] / NCM_SF_SPHERICAL_HARMONICS_EPS;
			}
#endif /* NCM_SF_SPHERICAL_HARMONICS_LATERAL_MOVE */
		}
	}
}

/* Reset methods */

NCM_INLINE void 
ncm_sf_spherical_harmonics_Y_reset (NcmSFSphericalHarmonicsY *sphaY)
{
	NcmSFSphericalHarmonics *spha = sphaY->spha;

	sphaY->l      = 0;
	sphaY->l0     = 0;
	sphaY->m      = 0;
	sphaY->Klm    = ncm_sf_spherical_harmonics_get_Klm (spha, sphaY->l, sphaY->m);
	sphaY->Pl0m   = ncm_c_sqrt_1_4pi () / NCM_SF_SPHERICAL_HARMONICS_EPS;
	sphaY->Pl0p1m = sphaY->x * _SN (2 * sphaY->l + 3) * sphaY->Pl0m;
	sphaY->Plm    = sphaY->Pl0m * NCM_SF_SPHERICAL_HARMONICS_EPS;
	sphaY->Plp1m  = sphaY->Pl0p1m * NCM_SF_SPHERICAL_HARMONICS_EPS;
}

NCM_INLINE void 
ncm_sf_spherical_harmonics_Y_array_reset (NcmSFSphericalHarmonicsYArray *sphaYa, const gint len)
{
	NcmSFSphericalHarmonics *spha = sphaYa->spha;
	gint i;

	g_assert_cmpuint (len, ==, sphaYa->len);
	
	sphaYa->l   = 0;
	sphaYa->l0  = 0;
	sphaYa->m   = 0;
	sphaYa->Klm = ncm_sf_spherical_harmonics_get_Klm (spha, sphaYa->l, sphaYa->m);

	for (i = 0; i < len; i++)
	{
		sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] = ncm_c_sqrt_1_4pi () / NCM_SF_SPHERICAL_HARMONICS_EPS;
		sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] = sphaYa->x[i] * _SN (2 * sphaYa->l + 3) * sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)];
		sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)]  = sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] * NCM_SF_SPHERICAL_HARMONICS_EPS;
		sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)]  = sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] * NCM_SF_SPHERICAL_HARMONICS_EPS;
	}
}

/* Start methods */

NCM_INLINE void 
ncm_sf_spherical_harmonics_start_rec (NcmSFSphericalHarmonics *spha, NcmSFSphericalHarmonicsY *sphaY, const gdouble theta)
{
	sincos (theta, &sphaY->sqrt1mx2, &sphaY->x);
	
	sphaY->l      = 0;
	sphaY->l0     = 0;
	sphaY->m      = 0;
	sphaY->Klm    = ncm_sf_spherical_harmonics_get_Klm (spha, sphaY->l, sphaY->m);
	sphaY->Pl0m   = ncm_c_sqrt_1_4pi () / NCM_SF_SPHERICAL_HARMONICS_EPS;
	sphaY->Pl0p1m = sphaY->x * _SN (2 * sphaY->l + 3) * sphaY->Pl0m;
	sphaY->Plm    = sphaY->Pl0m * NCM_SF_SPHERICAL_HARMONICS_EPS;
	sphaY->Plp1m  = sphaY->Pl0p1m * NCM_SF_SPHERICAL_HARMONICS_EPS;
}

NCM_INLINE void 
ncm_sf_spherical_harmonics_start_rec_array (NcmSFSphericalHarmonics *spha, NcmSFSphericalHarmonicsYArray *sphaYa, const gint len, const gdouble *theta)
{
	gint i;

	g_assert_cmpuint (len, ==, sphaYa->len);
	
	sphaYa->Klm = ncm_sf_spherical_harmonics_get_Klm (spha, 0, 0);
	sphaYa->l   = 0;
	sphaYa->l0  = 0;
	sphaYa->m   = 0;

	for (i = 0; i < len; i++)
	{
		sincos (theta[i], &sphaYa->sqrt1mx2[i], &sphaYa->x[i]);

		sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)] = ncm_c_sqrt_1_4pi () / NCM_SF_SPHERICAL_HARMONICS_EPS;
		sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] = sphaYa->x[i] * _SN (2 * sphaYa->l + 3) * sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)];
		sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)]  = sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 0, len)]   * NCM_SF_SPHERICAL_HARMONICS_EPS;
		sphaYa->Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)]  = sphaYa->Yl0m[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (i, 1, len)] * NCM_SF_SPHERICAL_HARMONICS_EPS;
	}
}

NCM_INLINE NcmSFSphericalHarmonicsK *
ncm_sf_spherical_harmonics_get_Klm (NcmSFSphericalHarmonics *spha, const gint l0, const gint m)
{
	GArray *Km_array = g_ptr_array_index (spha->K_array, m);
	return &g_array_index (Km_array, NcmSFSphericalHarmonicsK, l0 - m);
}

#undef NCM_SF_SPHERICAL_HARMONICS_EPS
#undef _SN
#undef _SNM1

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_SF_SPHERICAL_HARMONICS_INLINE_H_ */
