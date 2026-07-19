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

G_BEGIN_DECLS

#define NCM_TYPE_LAURENT_SERIES (ncm_laurent_series_get_type ())

GType ncm_laurent_series_get_type (void) G_GNUC_CONST;

typedef struct _NcmLaurentSeries NcmLaurentSeries;

struct _NcmLaurentSeries
{
  gint hmin;
  gint hmax;
  gint c_cap; /* allocated length of c, >= hmax-hmin+1; see ncm_laurent_series_reset() */
  NcmComplex *c;
  gatomicrefcount ref_count;
};

NcmLaurentSeries *ncm_laurent_series_new_single (gint h, NcmComplex val);
NcmLaurentSeries *ncm_laurent_series_new (gint hmin, gint hmax);
NcmLaurentSeries *ncm_laurent_series_copy (const NcmLaurentSeries *a);
NcmLaurentSeries *ncm_laurent_series_ref (NcmLaurentSeries *a);
void ncm_laurent_series_free (NcmLaurentSeries *a);
void ncm_laurent_series_clear (NcmLaurentSeries **a);

void ncm_laurent_series_reset (NcmLaurentSeries *a, gint hmin, gint hmax);
gint ncm_laurent_series_get_hmin (const NcmLaurentSeries *a);
gint ncm_laurent_series_get_hmax (const NcmLaurentSeries *a);
void ncm_laurent_series_get_ptr (const NcmLaurentSeries *a, gint h, NcmComplex *out);
void ncm_laurent_series_set_ptr (NcmLaurentSeries *a, gint h, const NcmComplex *val);

/* a + sb*b */
NcmLaurentSeries *ncm_laurent_series_add (const NcmLaurentSeries *a, const NcmLaurentSeries *b, gdouble sb);
NcmLaurentSeries *ncm_laurent_series_scale_ptr (const NcmLaurentSeries *a, const NcmComplex *s);
NcmLaurentSeries *ncm_laurent_series_conv (const NcmLaurentSeries *a, const NcmLaurentSeries *b);
NcmLaurentSeries *ncm_laurent_series_conj (const NcmLaurentSeries *a);

void ncm_laurent_series_eval_ptr (const NcmLaurentSeries *a, const NcmComplex *w, NcmComplex *out);
gdouble ncm_laurent_series_jacobi_anger_reduce (const NcmLaurentSeries *cm, gdouble phi, const gdouble *Ik, gint n_Ik);
NcmComplex ncm_laurent_series_get (const NcmLaurentSeries *a, gint h);

void ncm_laurent_series_set (NcmLaurentSeries *a, gint h, NcmComplex val);
NcmLaurentSeries *ncm_laurent_series_scale (const NcmLaurentSeries *a, NcmComplex s);
NcmComplex ncm_laurent_series_eval (const NcmLaurentSeries *a, NcmComplex w);

void ncm_laurent_series_set_single_into (NcmLaurentSeries *out, gint h, NcmComplex val);
void ncm_laurent_series_add_into (NcmLaurentSeries *out, const NcmLaurentSeries *a, const NcmLaurentSeries *b, gdouble sb);
void ncm_laurent_series_conv_into (NcmLaurentSeries *out, const NcmLaurentSeries *a, const NcmLaurentSeries *b);
void ncm_laurent_series_scale_into (NcmLaurentSeries *out, const NcmLaurentSeries *a, NcmComplex s);
void ncm_laurent_series_conj_into (NcmLaurentSeries *out, const NcmLaurentSeries *a);

/**
 * NcmLaurentSeriesTPS:
 *
 * A truncated power series of order $N$,
 * $\sum_{n=0}^N L_n\,g^n \mod g^{N+1}$, whose coefficients $L_n$ are
 * themselves #NcmLaurentSeries. Order is fixed at construction
 * (ncm_laurent_series_tps_new()) and the $N+1$ coefficients are owned
 * directly (no external pool): a #NcmLaurentSeriesTPS is meant to be a
 * long-lived, repeatedly-refilled object (see nc_wl_ellipticity_series.c),
 * not a short-lived per-call temporary.
 *
 * Boxed with reference-count "copy" (ncm_laurent_series_tps_ref()), not a
 * deep copy, matching #NcGalaxyShapeFactorData and siblings: a "copy"
 * shares the same underlying storage, so later mutations via the owning
 * object's own eval()/compute step are visible through every outstanding
 * reference.
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
void ncm_laurent_series_tps_eval_ptr (const NcmLaurentSeriesTPS *tps, const NcmComplex *w, const NcmComplex *g, NcmComplex *out);

/**
 * ncm_laurent_series_tps_pow:
 * @out: a #NcmLaurentSeriesTPS, same order as @a, must not alias @a
 * @a: a #NcmLaurentSeriesTPS whose order-0 coefficient $a_0=L_0(w=1)$
 * (a single harmonic-0 term) is nonzero
 * @p: the (real) exponent
 *
 * Raises the truncated power series $a(g)$ to the real power @p:
 * $out(g)=a(g)^p \mod g^{N+1}$. Only a plain `gdouble` exponent is
 * supported (not `complex double`): every current use (e.g.
 * #NcGalaxyShapePopBeta's own $\rho^{2(\alpha-1)}$ composition) only ever
 * needs a real exponent, and restricting to real keeps this fully
 * introspectable (no native-only guard needed at all).
 *
 * $a$ must have $a_0\ne0$: factor $a(g)=a_0(1+u(g))$ with $u(0)=0$, then
 * $(1+u)^p=\sum_n c_n g^n$ ($c_0=1$) follows the generalized-binomial
 * recursion $n c_n=\sum_{k=1}^n[kp-(n-k)]u_k c_{n-k}$ (from differentiating
 * $F=(1+u)^p$: $(1+u)F'=p u'F$), and $out=a_0^p\,\sum_n c_n g^n$.
 */
void ncm_laurent_series_tps_pow (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, gdouble p);

/* In-place tier: writes into a caller-supplied @out instead of allocating,
 * for the hot loop in nc_wl_ellipticity_series.c. Only conv() needs any
 * scratch beyond @out's own slots (the fold's ping-pong accumulator, see
 * ncm_laurent_series_tps_new()'s own comment on @conv_acc/@conv_term); it
 * draws that scratch from @out's own private fields, never from an external
 * pool. conj/add write straight into @out's existing slots via the plain
 * `_into` primitives, no scratch at all. */
void ncm_laurent_series_tps_conv (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, const NcmLaurentSeriesTPS *b);
void ncm_laurent_series_tps_conj (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a);
void ncm_laurent_series_tps_add (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, const NcmLaurentSeriesTPS *b, gdouble sb);

/* @s by value: fast, `(skip)`-ed from introspection (see the file's own top
 * doc comment on the bare/`_ptr` convention). No `_ptr` sibling exists yet
 * (nothing calls one). */
void ncm_laurent_series_tps_scale (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, NcmComplex s);
NcmComplex ncm_laurent_series_tps_eval (const NcmLaurentSeriesTPS *tps, NcmComplex w, NcmComplex g);

G_END_DECLS

#endif /* _NCM_LAURENT_SERIES_H_ */

