/***************************************************************************
 *            ncm_laurent_series.h
 *
 *  Tue Jul 8 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_laurent_series.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
#ifndef _NCM_LAURENT_SERIES_H_
#define _NCM_LAURENT_SERIES_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/core/ncm_util.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_LAURENT_SERIES (ncm_laurent_series_get_type ())

GType ncm_laurent_series_get_type (void) G_GNUC_CONST;

typedef struct _NcmLaurentSeries NcmLaurentSeries;

#ifndef NUMCOSMO_GIR_SCAN
struct _NcmLaurentSeries
{
  gint hmin;
  gint hmax;
  gint c_cap; /* allocated length of c, >= hmax-hmin+1; see ncm_laurent_series_reset() */
  complex double *c;
};
#endif /* NUMCOSMO_GIR_SCAN */

/* Native, `complex double`-based secondary interface -- what the hot loop in
 * nc_wl_ellipticity_series.c actually calls; direct struct-field access
 * (a->c[h-a->hmin], bounds-checked via ncm_laurent_series_get_c only) is
 * also fine internally within that file. */
#ifndef NUMCOSMO_GIR_SCAN
NcmLaurentSeries *ncm_laurent_series_new_single (gint h, complex double val);

#endif /* NUMCOSMO_GIR_SCAN */
NcmLaurentSeries *ncm_laurent_series_new (gint hmin, gint hmax);
NcmLaurentSeries *ncm_laurent_series_copy (const NcmLaurentSeries *a);
void ncm_laurent_series_free (NcmLaurentSeries *a);
void ncm_laurent_series_clear (NcmLaurentSeries **a);

#ifndef NUMCOSMO_GIR_SCAN
void ncm_laurent_series_reset (NcmLaurentSeries *a, gint hmin, gint hmax);

#endif /* NUMCOSMO_GIR_SCAN */

gint ncm_laurent_series_hmin (const NcmLaurentSeries *a);
gint ncm_laurent_series_hmax (const NcmLaurentSeries *a);

/*
 * Introspectable (#NcmComplex-based) API -- usable and testable from Python.
 * Every arithmetic operation returns a brand new #NcmLaurentSeries
 * (transfer full); see the native, `complex double`-based secondary
 * interface below (#ifndef NUMCOSMO_GIR_SCAN) for the same operations at
 * native-arithmetic speed, matching nc_wl_ellipticity.h's own
 * introspectable/`_c`-suffixed dual-interface convention.
 */
void ncm_laurent_series_get (const NcmLaurentSeries *a, gint h, NcmComplex *out);
void ncm_laurent_series_set (NcmLaurentSeries *a, gint h, const NcmComplex *val);

/* a + sb*b */
NcmLaurentSeries *ncm_laurent_series_add (const NcmLaurentSeries *a, const NcmLaurentSeries *b, gdouble sb);
NcmLaurentSeries *ncm_laurent_series_scale (const NcmLaurentSeries *a, const NcmComplex *s);
NcmLaurentSeries *ncm_laurent_series_conv (const NcmLaurentSeries *a, const NcmLaurentSeries *b);
NcmLaurentSeries *ncm_laurent_series_conj (const NcmLaurentSeries *a);

/* Evaluates $a(w)=\sum_h c_h w^h$ at a concrete point @w. */
void ncm_laurent_series_eval (const NcmLaurentSeries *a, const NcmComplex *w, NcmComplex *out);

/* Jacobi-Anger reduction of Int_0^2pi w^h * exp(z*(cos(theta-phi)-1)) dtheta,
 * summed against @cm's own harmonic content:
 * cm_0*I0(z) + 2*sum_{h>=1} Ih(z)*Re(cm_h*e^{i h phi}). @Ik must hold scaled
 * Bessel values (exp(-z)*I_h(z), e.g. gsl_sf_bessel_In_scaled) for h=0..n_Ik-1.
 * Already fully introspectable (no complex-valued parameters), so there is
 * no separate native variant of this one. */
gdouble ncm_laurent_series_jacobi_anger_reduce (const NcmLaurentSeries *cm, gdouble phi, const gdouble *Ik, gint n_Ik);

#ifndef NUMCOSMO_GIR_SCAN
complex double ncm_laurent_series_get_c (const NcmLaurentSeries *a, gint h);

void ncm_laurent_series_set_c (NcmLaurentSeries *a, gint h, complex double val);
NcmLaurentSeries *ncm_laurent_series_scale_c (const NcmLaurentSeries *a, complex double s);

complex double ncm_laurent_series_eval_c (const NcmLaurentSeries *a, complex double w);

/* Reuse (no-allocation-in-the-common-case) counterparts of add/conv/
 * scale_c/conj/new_single: write into a caller-supplied @out (which
 * ncm_laurent_series_reset()'s to the correct range internally) instead of
 * allocating a fresh #NcmLaurentSeries. @out must not alias any input.
 * Meant for hot loops that keep their own scratch #NcmLaurentSeries around
 * (see #NcmLaurentSeriesTPS below and nc_wl_ellipticity_series.c) -- the
 * "returns new" siblings above are unaffected and remain the
 * introspectable/general-purpose API. */
void ncm_laurent_series_set_single_into (NcmLaurentSeries *out, gint h, complex double val);
void ncm_laurent_series_add_into (NcmLaurentSeries *out, const NcmLaurentSeries *a, const NcmLaurentSeries *b, gdouble sb);
void ncm_laurent_series_conv_into (NcmLaurentSeries *out, const NcmLaurentSeries *a, const NcmLaurentSeries *b);
void ncm_laurent_series_scale_c_into (NcmLaurentSeries *out, const NcmLaurentSeries *a, complex double s);
void ncm_laurent_series_conj_into (NcmLaurentSeries *out, const NcmLaurentSeries *a);

#endif /* NUMCOSMO_GIR_SCAN */

/**
 * NcmLaurentSeriesTPS:
 *
 * A truncated power series of order $N$,
 * $\sum_{n=0}^N L_n\,g^n \mod g^{N+1}$, whose coefficients $L_n$ are
 * themselves #NcmLaurentSeries -- the self-describing-length container that
 * used to be a bare `NcmLaurentSeries **` plus a separately-threaded `guint
 * order` at every call site. Order is fixed at construction
 * (ncm_laurent_series_tps_new()) and the $N+1$ coefficients are owned
 * directly (no external pool): a #NcmLaurentSeriesTPS is meant to be a
 * long-lived, repeatedly-refilled object (see nc_wl_ellipticity_series.c),
 * not a short-lived per-call temporary.
 *
 * Boxed with reference-count "copy" (ncm_laurent_series_tps_ref()), not a
 * deep copy, matching #NcGalaxySDShapeData and siblings: a "copy" shares
 * the same underlying storage, so later mutations via the owning object's
 * own eval()/compute step are visible through every outstanding reference.
 *
 * Because this is a *truncated* power series, not a polynomial-ring
 * element, ncm_laurent_series_tps_conv() truncates its product at the
 * shared order of its operands rather than extending to the naive
 * deg(a)+deg(b) -- every op below requires @a, @b (if present) and @out to
 * already share the same order.
 */
typedef struct _NcmLaurentSeriesTPS NcmLaurentSeriesTPS;

#define NCM_TYPE_LAURENT_SERIES_TPS (ncm_laurent_series_tps_get_type ())

GType ncm_laurent_series_tps_get_type (void) G_GNUC_CONST;

NcmLaurentSeriesTPS *ncm_laurent_series_tps_new (guint order);
NcmLaurentSeriesTPS *ncm_laurent_series_tps_ref (NcmLaurentSeriesTPS *tps);
void ncm_laurent_series_tps_unref (NcmLaurentSeriesTPS *tps);
void ncm_laurent_series_tps_clear (NcmLaurentSeriesTPS **tps);

guint ncm_laurent_series_tps_order (const NcmLaurentSeriesTPS *tps);
NcmLaurentSeries *ncm_laurent_series_tps_get (const NcmLaurentSeriesTPS *tps, guint n);

/* Evaluates $\mathrm{tps}(w,g)=\sum_{n=0}^N L_n(w)\,g^n$ at concrete points
 * @w (the coefficients' own formal variable) and @g (the truncation
 * variable) -- e.g. a caller mapping @w back to some angle-parameterized
 * quantity (as nc_wl_ellipticity_series.c's own consumers do via
 * $w=e^{i\theta}$) does that mapping itself; this function only knows
 * about the two formal variables, not what they represent. */
void ncm_laurent_series_tps_eval (const NcmLaurentSeriesTPS *tps, const NcmComplex *w, const NcmComplex *g, NcmComplex *out);

#ifndef NUMCOSMO_GIR_SCAN

/* Native, hot-loop tier (not introspectable: @s below is a native complex
 * double) -- what nc_wl_ellipticity_series.c's evaluators actually call.
 * Only conv() needs any scratch beyond @out's own slots (the fold's
 * ping-pong accumulator, see ncm_laurent_series_tps_new()'s own comment on
 * @conv_acc/@conv_term); it draws that scratch from @out's own private
 * fields, never from an external pool, so no pool/workspace parameter
 * appears anywhere in this API. conj/add/scale write straight into @out's
 * existing slots via the plain `_into` primitives, no scratch at all. */
void ncm_laurent_series_tps_conv (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, const NcmLaurentSeriesTPS *b);
void ncm_laurent_series_tps_conj (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a);
void ncm_laurent_series_tps_add (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, const NcmLaurentSeriesTPS *b, gdouble sb);
void ncm_laurent_series_tps_scale (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, complex double s);

complex double ncm_laurent_series_tps_eval_c (const NcmLaurentSeriesTPS *tps, complex double w, complex double g);

#endif /* NUMCOSMO_GIR_SCAN */

G_END_DECLS

#endif /* _NCM_LAURENT_SERIES_H_ */

