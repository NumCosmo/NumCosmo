/***************************************************************************
 *            nc_galaxy_shape_factor_moment_series.c
 *
 *  Thu Sep 3 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_moment_series.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyShapeFactorMomentSeries:
 *
 * Truncated-$g$-power-series evaluation of the intrinsic-ellipticity marginal.
 *
 * Approximates
 * $$
 * P(\epsilon_\mathrm{obs} \mid g) = \int_{|\chi_I|<1} \mathrm{d}^2\chi_I\,
 *   P_\mathrm{pop}(\chi_I)\, N_2\big(\epsilon_\mathrm{obs} - S(g,\chi_I);
 *   \sigma_\mathrm{noise}^2\big)
 * $$
 * by the Gaussian carrying the model's own mean and covariance,
 * $\mu(g)=\langle S(g,\chi_I)\rangle$,
 * $C(g)=\mathrm{Cov}[S(g,\chi_I)]+\sigma_\mathrm{noise}^2 I_2$, with both
 * expectations truncated as power series in the reduced shear $g$ (order
 * `trunc-order`) rather than in the intrinsic width. The forward map
 * $S(g,\cdot)$ is a ratio of low-degree polynomials in $g$ for both
 * ellipticity conventions, so its coefficients follow a three-term
 * recursion in the tangential gauge ($g$ real; applied to complex $g$ by
 * rotating $(g,\epsilon_\mathrm{obs})$ together, exact, not an
 * approximation) -- no symbolic differentiation, arbitrary order is as
 * cheap as any low order. Under isotropy every population enters only
 * through its own radial moments
 * $M_{2k}=\langle|\chi_I|^{2k}\rangle$ (nc_galaxy_shape_pop_moment_2k()),
 * up to $M_{n+1}$ at order $n$.
 *
 * Unlike #NcGalaxyShapeFactorCGF and #NcGalaxyShapeFactorVarAdd, which
 * truncate in the intrinsic width and therefore require
 * nc_galaxy_shape_pop_get_sigma() (a Gaussian-only capability), this class
 * needs only radial moments -- every #NcGalaxyShapePop provides them,
 * through a closed form where one exists or the base class' own
 * fixed-quadrature default otherwise -- so it is the only closed-form
 * shape factor compatible with every population, including
 * #NcGalaxyShapePopBeta.
 *
 * A truncated-moment Gaussian carrying the population's own mean and
 * covariance, evaluated as functions of $g$, has exactly zero score bias
 * wherever those two moments are exact: non-Gaussianity of the true
 * marginal is not itself a source of shear bias here, only the moment
 * truncation is, and the truncation order is the accuracy knob (see the
 * design note this class was built from,
 * `MOMENT_SERIES.md` at the repository root, for the derivation and the
 * measured accuracy envelope, cross-checked at
 * tests/python/nc/lss/galaxy/test_galaxy_shape_factor_moment_series.py).
 *
 * The truncated covariance is not guaranteed positive at every order and
 * $g$: it is asserted, not silently clamped, so a breakdown reads as an
 * abort naming `trunc-order` and $|g|$ rather than a confidently wrong
 * number. At the default order (5) this does not occur anywhere in the
 * physical shear range; a caller lowering `trunc-order` should expect the
 * guard to become reachable well inside that range.
 *
 * The score (log-likelihood derivative in $g$) of this model can have more
 * than one root in a wide bracket, of which only one is the maximum: an
 * optimizer using this class must maximise the likelihood (e.g. a coarse
 * scan then a local refinement), never root-find the score directly.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_moment_series.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

/*
 * ---- Population-independent shear-map series algebra ----
 *
 * Builds, once per (ellip-conv, trunc-order) at construction time, the
 * rational coefficient tables that turn a population's radial moments
 * M_2k into this order's mean/covariance g-series coefficients. c[j] is
 * the coefficient array of chi_I^a * chibar_I^b for g^j in the
 * tangential-gauge (real g) forward map S(g,chi_I); the three-term
 * recursion below (and its single-term TRACE_DET analogue) needs no
 * symbolic differentiation, so building it to any order is O(n^3) in a
 * small dense array, run once per class instance, never per galaxy.
 *
 * The TRACE_DET recursion here is c_j = n_j - chi_I c_{j-1} (n_0=chi_I,
 * n_1=1), matching nc_wl_ellipticity_apply_shear_trace_det()'s real map
 * (e+g)/(1+conj(g)*e) -- holomorphic in the intrinsic ellipticity, with
 * the conjugate on g. This differs from a c_j = n_j - chibar_I c_{j-1}
 * recursion, which would instead correspond to a D=1+g*chibar_I
 * denominator not used anywhere in this codebase; that alternative gives
 * <epsilon_obs> = (1-M_2) g, at odds with nc_galaxy_shape_factor_direct_estimate()'s
 * TRACE_DET branch (<epsilon>=g, no responsivity factor) and with
 * NcGalaxyShapeFactorCGF's isotropic TRACE_DET response moments -- the
 * map used here reproduces both.
 */

typedef struct _Poly2
{
  gdouble *c; /* c[a*sz+b]: coefficient of chi_I^a * chibar_I^b */
  guint sz;
} Poly2;

static Poly2
poly2_new (guint sz)
{
  Poly2 p;

  p.sz = sz;
  p.c  = g_new0 (gdouble, sz * sz);

  return p;
}

static void
poly2_clear (Poly2 *p)
{
  g_clear_pointer (&p->c, g_free);
}

static inline gdouble
poly2_get (const Poly2 *p, guint a, guint b)
{
  return p->c[a * p->sz + b];
}

static inline void
poly2_addto (Poly2 *p, guint a, guint b, gdouble v)
{
  p->c[a * p->sz + b] += v;
}

/* out += scale * in, shifted by (+1,0): multiplication by chi_I. */
static void
poly2_axpy_mul_a (Poly2 *out, const Poly2 *in, gdouble scale)
{
  guint a, b;

  for (a = 0; a + 1 < out->sz; a++)
    for (b = 0; b < out->sz; b++)
      poly2_addto (out, a + 1, b, scale * poly2_get (in, a, b));
}

/* out += scale * in, shifted by (0,+1): multiplication by chibar_I. */
static void
poly2_axpy_mul_b (Poly2 *out, const Poly2 *in, gdouble scale)
{
  guint a, b;

  for (a = 0; a < out->sz; a++)
    for (b = 0; b + 1 < out->sz; b++)
      poly2_addto (out, a, b + 1, scale * poly2_get (in, a, b));
}

/* out += scale * in, unshifted. */
static void
poly2_axpy (Poly2 *out, const Poly2 *in, gdouble scale)
{
  guint i;

  for (i = 0; i < out->sz * out->sz; i++)
    out->c[i] += scale * in->c[i];
}

/* out += scale * (A * B), as a formal 2D polynomial product; terms landing
 * outside @out's array are dropped (never happens here: @out is always
 * sized for the full product). */
static void
poly2_axpy_conv (Poly2 *out, const Poly2 *A, const Poly2 *B, gdouble scale)
{
  guint a1, b1;

  for (a1 = 0; a1 < A->sz; a1++)
  {
    for (b1 = 0; b1 < A->sz; b1++)
    {
      const gdouble va = poly2_get (A, a1, b1);
      guint a2, b2;

      if (va == 0.0)
        continue;

      for (a2 = 0; a2 + a1 < out->sz; a2++)
      {
        for (b2 = 0; b2 + b1 < out->sz; b2++)
        {
          const gdouble vb = poly2_get (B, a2, b2);

          if (vb != 0.0)
            poly2_addto (out, a1 + a2, b1 + b2, scale * va * vb);
        }
      }
    }
  }
}

/* out[a][b] += in[b][a]: formal complex conjugation (chi_I <-> chibar_I),
 * valid because every coefficient here is real. */
static void
poly2_axpy_transpose (Poly2 *out, const Poly2 *in)
{
  guint a, b;

  for (a = 0; a < out->sz; a++)
    for (b = 0; b < out->sz; b++)
      poly2_addto (out, a, b, poly2_get (in, b, a));
}

/* c_j = n_j - (chi_I+chibar_I) c_{j-1} - c_{j-2}, n_0=chi_I, n_1=2,
 * n_2=chibar_I, n_{j>2}=0, c_{j<0}=0 (MOMENT_SERIES.md eq. 5). */
static void
_build_c_trace (Poly2 *c, guint n, guint sz)
{
  guint j;

  for (j = 0; j <= n; j++)
  {
    c[j] = poly2_new (sz);

    if (j == 0)
      poly2_addto (&c[j], 1, 0, 1.0);
    else if (j == 1)
      poly2_addto (&c[j], 0, 0, 2.0);
    else if (j == 2)
      poly2_addto (&c[j], 0, 1, 1.0);

    if (j >= 1)
    {
      poly2_axpy_mul_a (&c[j], &c[j - 1], -1.0);
      poly2_axpy_mul_b (&c[j], &c[j - 1], -1.0);
    }

    if (j >= 2)
      poly2_axpy (&c[j], &c[j - 2], -1.0);
  }
}

/* c_j = n_j - chi_I c_{j-1}, n_0=chi_I, n_1=1, n_{j>1}=0 -- see this file's
 * top-of-file note on why the shift is by chi_I, not chibar_I. */
static void
_build_c_trace_det (Poly2 *c, guint n, guint sz)
{
  guint j;

  for (j = 0; j <= n; j++)
  {
    c[j] = poly2_new (sz);

    if (j == 0)
      poly2_addto (&c[j], 1, 0, 1.0);
    else if (j == 1)
      poly2_addto (&c[j], 0, 0, 1.0);

    if (j >= 1)
      poly2_axpy_mul_a (&c[j], &c[j - 1], -1.0);
  }
}

struct _NcGalaxyShapeFactorMomentSeries
{
  NcGalaxyShapeFactor parent_instance;
};

typedef struct _NcGalaxyShapeFactorMomentSeriesPrivate
{
  guint trunc_order;

  /* Table shapes: n_m mean coefficients (odd powers of g, m_{2l+1}), n_v
   * covariance coefficients (even powers, v_{2l}/w_{2l}), each a linear
   * functional of n_moments radial moments M_0..M_{2(n_moments-1)}. */
  guint n_m;
  guint n_v;
  guint n_moments;

  gdouble *tab_m; /* n_m * n_moments, row-major */
  gdouble *tab_v; /* n_v * n_moments */
  gdouble *tab_w; /* n_v * n_moments */

  /* Population generation, refreshed in prepare(); see _peek_coeffs()'s
   * own comment for why this is read here rather than through
   * nc_galaxy_shape_factor_get_pop_hash() per evaluation (same reasoning
   * as NcGalaxyShapeFactorCGF's own pop_hash field). */
  guint64 pop_hash;
} NcGalaxyShapeFactorMomentSeriesPrivate;

/*
 * Per-galaxy scratch: the coefficient triple {m,v,w}, constant across a
 * galaxy's z-nodes. Invalidated on the same two axes as CGF's own
 * per-galaxy V cache: @pop_hash_seen catches a population MODEL change,
 * while ldata_read_row() clears @valid because a per-galaxy population
 * (NcGalaxyShapePopGaussLocal) takes its moments from the catalog row --
 * data, which moves with no model pkey bump at all.
 */
typedef struct _NcGalaxyShapeFactorMomentSeriesLData
{
  gdouble *m; /* n_m */
  gdouble *v; /* n_v */
  gdouble *w; /* n_v */
  guint64 pop_hash_seen;
  gboolean valid;
} NcGalaxyShapeFactorMomentSeriesLData;

enum
{
  PROP_0,
  PROP_TRUNC_ORDER,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorMomentSeries, nc_galaxy_shape_factor_moment_series, NC_TYPE_GALAXY_SHAPE_FACTOR)

/* Builds tab_m/tab_v/tab_w for the given convention and self->trunc_order.
 * The raw covariance tables (tab_v/tab_w) hold <(Re S)^2>_j and
 * <(Im S)^2>_j themselves, linear in the moments; the mu*mu subtraction
 * that turns <(Re S)^2>_j into Var(Re S)_j is population-dependent (it
 * uses the concrete m_l values), so it happens per galaxy in
 * _peek_coeffs(), not here -- see that function's own comment. */
static void
_nc_galaxy_shape_factor_moment_series_build_tables (NcGalaxyShapeFactorMomentSeriesPrivate * const self,
                                                    NcGalaxyWLObsEllipConv                         ellip_conv)
{
  const guint n   = self->trunc_order;
  const guint sz  = n + 3;
  Poly2 *c        = g_new (Poly2, n + 1);
  const guint n_m = (n + 1) / 2;
  const guint n_v = n / 2 + 1;
  gdouble *tab_m  = g_new0 (gdouble, n_m * sz);
  gdouble *tab_v  = g_new0 (gdouble, n_v * sz);
  gdouble *tab_w  = g_new0 (gdouble, n_v * sz);
  guint max_a     = 0;
  guint n_moments;
  guint j, l, a;

  switch (ellip_conv)
  {
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
      _build_c_trace (c, n, sz);
      break;
    case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
      _build_c_trace_det (c, n, sz);
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }

  self->n_m = n_m;
  self->n_v = n_v;

  /* Mean series: m_j = diag(c_j), j odd, l = (j-1)/2 -- reading the
   * diagonal is exactly the isotropic average <chi_I^a chibar_I^a> = M_2a
   * (MOMENT_SERIES.md sec. 3.4). */
  for (l = 0; l < n_m; l++)
  {
    j = 2 * l + 1;

    for (a = 0; a < sz; a++)
    {
      const gdouble val = poly2_get (&c[j], a, a);

      tab_m[l * sz + a] = val;

      if (val != 0.0)
        max_a = MAX (max_a, a);
    }
  }

  /* Covariance series: raw_v_j = diag((c*c + c*conj(c))/2) = <(Re S)^2>_j,
   * raw_w_j = diag((c*conj(c) - c*c)/2) = <(Im S)^2>_j, j even, l = j/2. */
  for (l = 0; l < n_v; l++)
  {
    Poly2 p2 = poly2_new (sz); /* sum_i c_i * c_{j-i}:       <S^2>_j    */
    Poly2 ps = poly2_new (sz); /* sum_i c_i * conj(c_{j-i}): <|S|^2>_j  */
    guint i;

    j = 2 * l;

    for (i = 0; i <= j; i++)
    {
      Poly2 ct = poly2_new (sz);

      poly2_axpy_transpose (&ct, &c[j - i]);
      poly2_axpy_conv (&p2, &c[i], &c[j - i], 1.0);
      poly2_axpy_conv (&ps, &c[i], &ct, 1.0);
      poly2_clear (&ct);
    }

    for (a = 0; a < sz; a++)
    {
      const gdouble p2aa = poly2_get (&p2, a, a);
      const gdouble psaa = poly2_get (&ps, a, a);
      const gdouble vval = 0.5 * (p2aa + psaa);
      const gdouble wval = 0.5 * (psaa - p2aa);

      tab_v[l * sz + a] = vval;
      tab_w[l * sz + a] = wval;

      if ((vval != 0.0) || (wval != 0.0))
        max_a = MAX (max_a, a);
    }

    poly2_clear (&p2);
    poly2_clear (&ps);
  }

  n_moments = max_a + 1;

  self->n_moments = n_moments;
  self->tab_m     = g_new (gdouble, n_m * n_moments);
  self->tab_v     = g_new (gdouble, n_v * n_moments);
  self->tab_w     = g_new (gdouble, n_v * n_moments);

  for (l = 0; l < n_m; l++)
    for (a = 0; a < n_moments; a++)
      self->tab_m[l * n_moments + a] = tab_m[l * sz + a];

  for (l = 0; l < n_v; l++)
  {
    for (a = 0; a < n_moments; a++)
    {
      self->tab_v[l * n_moments + a] = tab_v[l * sz + a];
      self->tab_w[l * n_moments + a] = tab_w[l * sz + a];
    }
  }

  g_free (tab_m);
  g_free (tab_v);
  g_free (tab_w);

  for (j = 0; j <= n; j++)
    poly2_clear (&c[j]);

  g_free (c);
}

static void
nc_galaxy_shape_factor_moment_series_init (NcGalaxyShapeFactorMomentSeries *gsfms)
{
  NcGalaxyShapeFactorMomentSeriesPrivate * const self = nc_galaxy_shape_factor_moment_series_get_instance_private (gsfms);

  self->trunc_order = 5;
  self->n_m         = 0;
  self->n_v         = 0;
  self->n_moments   = 0;
  self->tab_m       = NULL;
  self->tab_v       = NULL;
  self->tab_w       = NULL;
  self->pop_hash    = 0;
}

static void
_nc_galaxy_shape_factor_moment_series_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorMomentSeries *gsfms              = NC_GALAXY_SHAPE_FACTOR_MOMENT_SERIES (object);
  NcGalaxyShapeFactorMomentSeriesPrivate * const self = nc_galaxy_shape_factor_moment_series_get_instance_private (gsfms);

  switch (prop_id)
  {
    case PROP_TRUNC_ORDER:
      self->trunc_order = g_value_get_uint (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_moment_series_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorMomentSeries *gsfms              = NC_GALAXY_SHAPE_FACTOR_MOMENT_SERIES (object);
  NcGalaxyShapeFactorMomentSeriesPrivate * const self = nc_galaxy_shape_factor_moment_series_get_instance_private (gsfms);

  switch (prop_id)
  {
    case PROP_TRUNC_ORDER:
      g_value_set_uint (value, self->trunc_order);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_moment_series_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_moment_series_parent_class)->constructed (object);
  {
    NcGalaxyShapeFactor *gsf                            = NC_GALAXY_SHAPE_FACTOR (object);
    NcGalaxyShapeFactorMomentSeriesPrivate * const self = nc_galaxy_shape_factor_moment_series_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENT_SERIES (object));
    const NcGalaxyWLObsEllipConv ellip_conv             = nc_galaxy_shape_factor_get_ellip_conv (gsf);

    /* Departure from a per-population build: the coefficient tables built
     * here depend only on (ellip-conv, trunc-order), never on the
     * population, so they are built once per class instance instead of
     * once per prepare() call. Under a per-galaxy population
     * (NcGalaxyShapePopGaussLocal) prepare() runs per galaxy, so building
     * these O(n^3) tables there would repeat that cost per galaxy per MCMC
     * step instead of once at construction. */
    _nc_galaxy_shape_factor_moment_series_build_tables (self, ellip_conv);
  }
}

static void
_nc_galaxy_shape_factor_moment_series_finalize (GObject *object)
{
  NcGalaxyShapeFactorMomentSeriesPrivate * const self = nc_galaxy_shape_factor_moment_series_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENT_SERIES (object));

  g_clear_pointer (&self->tab_m, g_free);
  g_clear_pointer (&self->tab_v, g_free);
  g_clear_pointer (&self->tab_w, g_free);

  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_moment_series_parent_class)->finalize (object);
}

static void
_nc_galaxy_shape_factor_moment_series_ldata_destroy (gpointer p)
{
  NcGalaxyShapeFactorMomentSeriesLData *ldata = (NcGalaxyShapeFactorMomentSeriesLData *) p;

  g_free (ldata->m);
  g_free (ldata->v);
  g_free (ldata->w);
  g_free (ldata);
}

static void
_nc_galaxy_shape_factor_moment_series_ldata_noop (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
}

/* A per-galaxy population reads its moments straight from the catalog row
 * (through nc_galaxy_shape_pop_moment_2k()'s own @data), so a new row
 * invalidates the cached coefficients without any model pkey moving --
 * mirrors NcGalaxyShapeFactorCGF's own ldata_read_row(). */
static void
_nc_galaxy_shape_factor_moment_series_ldata_read_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyShapeFactorMomentSeriesLData *ldata = (NcGalaxyShapeFactorMomentSeriesLData *) data->ldata;

  ldata->valid = FALSE;
}

static void
_nc_galaxy_shape_factor_moment_series_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_factor_moment_series_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorMomentSeriesPrivate * const self = nc_galaxy_shape_factor_moment_series_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENT_SERIES (gsf));
  NcGalaxyShapeFactorMomentSeriesLData *ldata         = g_new0 (NcGalaxyShapeFactorMomentSeriesLData, 1);

  /* g_new0 leaves @valid FALSE, so the first evaluation populates it. */
  ldata->m = g_new0 (gdouble, self->n_m);
  ldata->v = g_new0 (gdouble, self->n_v);
  ldata->w = g_new0 (gdouble, self->n_v);

  data->ldata                  = ldata;
  data->ldata_destroy          = &_nc_galaxy_shape_factor_moment_series_ldata_destroy;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_moment_series_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_moment_series_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_moment_series_ldata_required_columns;
}

static void
_nc_galaxy_shape_factor_moment_series_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
  NcGalaxyShapeFactorMomentSeriesPrivate * const self = nc_galaxy_shape_factor_moment_series_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENT_SERIES (gsf));

  /* No capability gate here, unlike CGF/VarAdd's prepare(): this class
   * needs only nc_galaxy_shape_pop_moment_2k(), which every
   * NcGalaxyShapePop provides (a closed form where one exists, the base
   * class' own fixed-quadrature default otherwise), so it is compatible
   * with every population, including NcGalaxyShapePopBeta.
   *
   * The base sets its own pop_hash before calling this vfunc, so this
   * snapshot is current; the per-galaxy coefficient caches are keyed on
   * it. */
  self->pop_hash = nc_galaxy_shape_factor_get_pop_hash (gsf);
}

/*
 * Refreshes the cached coefficient triple {m,v,w} from this galaxy's
 * current radial moments, when the population model generation moved or a
 * new catalog row was read into @data (see NcGalaxyShapeFactorMomentSeriesLData's
 * own comment). The raw covariance tables give <(Re S)^2>_j / <(Im S)^2>_j;
 * turning the former into Var(Re S)_j needs the mu*mu subtraction, a
 * self-convolution of the just-computed m array (MOMENT_SERIES.md sec. 4,
 * and this file's own build_tables() comment) -- Var(Im S)_j needs none,
 * since Im(mu(g))=0 identically for real g.
 */
static inline void
_nc_galaxy_shape_factor_moment_series_peek_coeffs (NcGalaxyShapeFactorMomentSeriesPrivate * const self,
                                                   NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                                   const gdouble **m_out, const gdouble **v_out, const gdouble **w_out)
{
  NcGalaxyShapeFactorMomentSeriesLData *ldata = (NcGalaxyShapeFactorMomentSeriesLData *) data->ldata;

  if (G_UNLIKELY (!ldata->valid || (ldata->pop_hash_seen != self->pop_hash)))
  {
    gdouble *M = g_new (gdouble, self->n_moments);
    guint l, a;

    for (a = 0; a < self->n_moments; a++)
      M[a] = nc_galaxy_shape_pop_moment_2k (pop, data->pop_data, a);

    for (l = 0; l < self->n_m; l++)
    {
      gdouble s = 0.0;

      for (a = 0; a < self->n_moments; a++)
        s += self->tab_m[l * self->n_moments + a] * M[a];

      ldata->m[l] = s;
    }

    for (l = 0; l < self->n_v; l++)
    {
      gdouble raw_v = 0.0;
      gdouble raw_w = 0.0;
      gdouble mu2   = 0.0;
      guint p;

      for (a = 0; a < self->n_moments; a++)
      {
        raw_v += self->tab_v[l * self->n_moments + a] * M[a];
        raw_w += self->tab_w[l * self->n_moments + a] * M[a];
      }

      for (p = 0; p < l; p++)
        mu2 += ldata->m[p] * ldata->m[l - 1 - p];

      ldata->v[l] = raw_v - mu2;
      ldata->w[l] = raw_w;
    }

    g_free (M);

    ldata->pop_hash_seen = self->pop_hash;
    ldata->valid         = TRUE;
  }

  *m_out = ldata->m;
  *v_out = ldata->v;
  *w_out = ldata->w;
}

/*
 * Gauge-fixes (g,eps_obs) together by -arg(g) (a rotation, exact -- see
 * this file's class doc comment), then Horner-evaluates the mean and
 * covariance series in u=|g|^2 and closes the diagonal bivariate Gaussian
 * (diagonal in the gauge-fixed frame: the tangential/cross cross term is
 * identically zero, unlike NcGalaxyShapeFactorCGF's TRACE covariance).
 * The want_log=FALSE branch avoids log(det) the way CGF's own closer
 * does, since eval_at_nodes() hits this once per galaxy per z-node.
 */
static gdouble
_nc_galaxy_shape_factor_moment_series_eval (NcGalaxyShapeFactorMomentSeriesPrivate * const self,
                                            const gdouble *m, const gdouble *v, const gdouble *w,
                                            const gdouble g_1, const gdouble g_2,
                                            const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                            const gdouble std_noise, const gboolean want_log)
{
  const gdouble g_mag  = hypot (g_1, g_2);
  const gdouble cos_pg = (g_mag > 0.0) ? g_1 / g_mag : 1.0;
  const gdouble sin_pg = (g_mag > 0.0) ? g_2 / g_mag : 0.0;
  const gdouble dx0    = epsilon_obs_1 * cos_pg + epsilon_obs_2 * sin_pg;
  const gdouble dy0    = -epsilon_obs_1 * sin_pg + epsilon_obs_2 * cos_pg;
  const gdouble u      = g_mag * g_mag;
  const gdouble sn2    = gsl_pow_2 (std_noise);
  gdouble mu, Ct, Cx, dx, dy, chi2;
  gint l;

  mu = m[self->n_m - 1];

  for (l = (gint) self->n_m - 2; l >= 0; l--)
    mu = mu * u + m[l];

  mu *= g_mag;

  Ct = v[self->n_v - 1];
  Cx = w[self->n_v - 1];

  for (l = (gint) self->n_v - 2; l >= 0; l--)
  {
    Ct = Ct * u + v[l];
    Cx = Cx * u + w[l];
  }

  Ct += sn2;
  Cx += sn2;

  if (G_UNLIKELY ((Ct <= 0.0) || (Cx <= 0.0)))
    g_error ("NcGalaxyShapeFactorMomentSeries: the truncated covariance is not "
             "positive (C_t=%g, C_x=%g) at trunc-order=%u, |g|=%g. This marks a "
             "truncation breakdown outside this order's valid range, not a "
             "physical marginal -- raise trunc-order or restrict the shear range.",
             Ct, Cx, self->trunc_order, g_mag);

  dx = dx0 - mu;
  dy = dy0;

  chi2 = dx * dx / Ct + dy * dy / Cx;

  if (want_log)
    return -0.5 * chi2 - log (2.0 * M_PI) - 0.5 * log (Ct) - 0.5 * log (Cx);
  else
    return exp (-0.5 * chi2) / (2.0 * M_PI * sqrt (Ct * Cx));
}

static gdouble
_nc_galaxy_shape_factor_moment_series_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorMomentSeriesPrivate * const self = nc_galaxy_shape_factor_moment_series_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENT_SERIES (gsf));
  const gdouble *m, *v, *w;

  _nc_galaxy_shape_factor_moment_series_peek_coeffs (self, pop, data, &m, &v, &w);

  return _nc_galaxy_shape_factor_moment_series_eval (self, m, v, w, g_1, g_2, epsilon_obs_1, epsilon_obs_2, data->std_noise, FALSE);
}

static gdouble
_nc_galaxy_shape_factor_moment_series_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorMomentSeriesPrivate * const self = nc_galaxy_shape_factor_moment_series_get_instance_private (NC_GALAXY_SHAPE_FACTOR_MOMENT_SERIES (gsf));
  const gdouble *m, *v, *w;

  _nc_galaxy_shape_factor_moment_series_peek_coeffs (self, pop, data, &m, &v, &w);

  return _nc_galaxy_shape_factor_moment_series_eval (self, m, v, w, g_1, g_2, epsilon_obs_1, epsilon_obs_2, data->std_noise, TRUE);
}

static void
nc_galaxy_shape_factor_moment_series_class_init (NcGalaxyShapeFactorMomentSeriesClass *klass)
{
  NcGalaxyShapeFactorClass *gsf_class = NC_GALAXY_SHAPE_FACTOR_CLASS (klass);
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_shape_factor_moment_series_set_property;
  object_class->get_property = &_nc_galaxy_shape_factor_moment_series_get_property;
  object_class->constructed  = &_nc_galaxy_shape_factor_moment_series_constructed;
  object_class->finalize     = &_nc_galaxy_shape_factor_moment_series_finalize;

  /**
   * NcGalaxyShapeFactorMomentSeries:trunc-order:
   *
   * Truncation order $n$ of the $g$-power series (both conventions'
   * closed forms generate any order equally easily). The covariance
   * converges an order of magnitude slower than the mean, and is
   * therefore the binding constraint on the useful order; there is little
   * accuracy to gain past $n\approx7$. Default 5.
   */
  g_object_class_install_property (object_class,
                                   PROP_TRUNC_ORDER,
                                   g_param_spec_uint ("trunc-order",
                                                      "Truncation order",
                                                      "Truncation order n of the g-power series",
                                                      1, G_MAXUINT, 5,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  gsf_class->data_init        = &_nc_galaxy_shape_factor_moment_series_data_init;
  gsf_class->prepare          = &_nc_galaxy_shape_factor_moment_series_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_moment_series_eval_marginal;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_moment_series_eval_ln_marginal;
}

/**
 * nc_galaxy_shape_factor_moment_series_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 * @trunc_order: truncation order $n$ of the $g$-power series, $n\ge1$
 *
 * Creates a new #NcGalaxyShapeFactorMomentSeries.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorMomentSeries.
 */
NcGalaxyShapeFactorMomentSeries *
nc_galaxy_shape_factor_moment_series_new (NcGalaxyWLObsEllipConv ellip_conv, guint trunc_order)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_MOMENT_SERIES,
                       "ellip-conv", ellip_conv,
                       "trunc-order", trunc_order,
                       NULL);
}

/**
 * nc_galaxy_shape_factor_moment_series_ref:
 * @gsfms: a #NcGalaxyShapeFactorMomentSeries
 *
 * Increases the reference count of @gsfms by one.
 *
 * Returns: (transfer full): @gsfms.
 */
NcGalaxyShapeFactorMomentSeries *
nc_galaxy_shape_factor_moment_series_ref (NcGalaxyShapeFactorMomentSeries *gsfms)
{
  return g_object_ref (gsfms);
}

/**
 * nc_galaxy_shape_factor_moment_series_free:
 * @gsfms: a #NcGalaxyShapeFactorMomentSeries
 *
 * Decreases the reference count of @gsfms by one.
 *
 */
void
nc_galaxy_shape_factor_moment_series_free (NcGalaxyShapeFactorMomentSeries *gsfms)
{
  g_object_unref (gsfms);
}

/**
 * nc_galaxy_shape_factor_moment_series_clear:
 * @gsfms: a #NcGalaxyShapeFactorMomentSeries
 *
 * Decreases the reference count of *@gsfms by one, and sets the pointer
 * *@gsfms to NULL.
 *
 */
void
nc_galaxy_shape_factor_moment_series_clear (NcGalaxyShapeFactorMomentSeries **gsfms)
{
  g_clear_object (gsfms);
}

