/***************************************************************************
 *            nc_galaxy_shape_pop_beta.c
 *
 *  Thu Jun 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop_beta.c
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyShapePopBeta:
 *
 * Beta intrinsic ellipticity distribution (BetaGlobal).
 *
 * The intrinsic ellipticity modulus $r = |\chi_I|$ follows
 * $r \sim \mathrm{Beta}(\alpha,\beta)$:
 * $$P_\mathrm{pop}(r) = \frac{r^{\alpha-1}(1-r)^{\beta-1}}{B(\alpha,\beta)}, \qquad r\in[0,1).$$
 *
 * The modulus is the coordinate the weak-lensing literature reports fits
 * against, so $\alpha$ and $\beta$ here are directly comparable to published
 * values. eval_p() returns the expression above with no change of variables
 * and no Jacobian.
 *
 * $\alpha$ and $\beta$ are both bounded to $\ge 1$, and $P_\mathrm{pop}(r)$
 * never diverges: it vanishes at $r=0$ for $\alpha>1$, is a nonzero constant
 * at $\alpha=1$, and behaves symmetrically at $r=1$ in $\beta$. What does
 * diverge for $\alpha<2$ is the derived 2D area density
 * $P_\mathrm{pop}(r)/(2\pi r)$. This class does not compute that; its
 * consumers do (#NcGalaxyShapeFactorFixedQuad, #NcGalaxyShapeFactorQuad), and
 * only as a pointwise artifact of their own quadrature, which is not polar in
 * $\chi_I$. See #NcGalaxyShapePop for why the exact 2D integral stays finite
 * through $r=0$. FixedQuad's cache is guarded against this by
 * `g_spline_pop_safe`.
 *
 * #NcGalaxyShapeFactorSeriesLensed's $g$-Taylor composition needs
 * $\sqrt{x(g)}$, whose branch point at $x(g)=0$ shrinks the series' radius of
 * convergence to an unusable value for any $\alpha<2$, including this class's
 * $\alpha=1.4$ default. SeriesLensed callers should keep $\alpha\ge2$; see
 * `docs/theory/wl_shape_factor_history.md`.
 *
 * nc_galaxy_shape_pop_beta_get_e_rms() and
 * nc_galaxy_shape_pop_beta_get_mode() report
 * $e_\mathrm{rms}=\sqrt{\langle x\rangle/2}$, with $\langle x\rangle=\langle
 * r^2\rangle=\alpha(\alpha+1)/[(\alpha+\beta)(\alpha+\beta+1)]$, and
 * $\mathrm{mode}(r) = (\alpha-1)/(\alpha+\beta-2)$, the argmax of
 * $P_\mathrm{pop}(r)$.
 * The shape parameters satisfy $\alpha,\beta \ge 1$. The class also exposes
 * $e_\mathrm{rms}=\sqrt{\langle r^2\rangle/2}$ and
 * $\mathrm{mode}(r)=(\alpha-1)/(\alpha+\beta-2)$ as derived quantities.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_pop_beta.h"
#include "ncm/model/ncm_mset_func_list.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <complex.h>
#include <gsl/gsl_sf_gamma.h>
#endif /* NUMCOSMO_GIR_SCAN */

#include "ncm/algebra/ncm_laurent_series.h"

struct _NcGalaxyShapePopBeta
{
  NcGalaxyShapePop parent_instance;
};

typedef struct _NcGalaxyShapePopBetaLData
{
  gdouble lnnorm_x; /* -ln(2 B(alpha,beta)): x-space normalization, used only by eval_p_rho2_g_series */
  gdouble lnnorm_r; /* -ln B(alpha,beta): eval_p()/eval_p_array()'s own P_pop(r) normalization */
  gdouble alpha;    /* cached copy of the alpha model param */
  gdouble beta;     /* cached copy of the beta model param */
} NcGalaxyShapePopBetaLData;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxyShapePopBeta, nc_galaxy_shape_pop_beta, NC_TYPE_GALAXY_SHAPE_POP);

static void
nc_galaxy_shape_pop_beta_init (NcGalaxyShapePopBeta *gspb)
{
}

static void
_nc_galaxy_shape_pop_beta_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_shape_pop_beta_parent_class)->finalize (object);
}

static void _nc_galaxy_shape_pop_beta_data_init (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
static void _nc_galaxy_shape_pop_beta_prepare (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
static gdouble _nc_galaxy_shape_pop_beta_ldata_get_mode_r (NcGalaxyShapePopData *data);
static gdouble _nc_galaxy_shape_pop_beta_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble r);
static void _nc_galaxy_shape_pop_beta_eval_p_array (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const GArray *r, GArray **p);
static void _nc_galaxy_shape_pop_beta_gen (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2);
static gdouble _nc_galaxy_shape_pop_beta_e_rms (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
static gdouble _nc_galaxy_shape_pop_beta_exponent_at_origin (NcGalaxyShapePop *gsp);
static gdouble _nc_galaxy_shape_pop_beta_moment_2k (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const guint k);
static void _nc_galaxy_shape_pop_beta_eval_p_rho2_g_series (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data,
                                                            const NcmLaurentSeriesTPS *x_series, NcmLaurentSeriesTPS *out);

static void
nc_galaxy_shape_pop_beta_class_init (NcGalaxyShapePopBetaClass *klass)
{
  NcGalaxyShapePopClass *gsp_class = NC_GALAXY_SHAPE_POP_CLASS (klass);
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_shape_pop_beta_finalize;

  ncm_model_class_set_name_nick (model_class, "Beta intrinsic ellipticity distribution", "BetaIntrinsic");
  ncm_model_class_add_params (model_class, NC_GALAXY_SHAPE_POP_BETA_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxyShapePopBeta:alpha:
   *
   * The shape parameter $\alpha$ of the Beta distribution of $r = |\chi_I|$.
   * Bounded to $\ge 1$ (see the class documentation).
   *
   */
  ncm_model_class_set_sparam (model_class,
                              NC_GALAXY_SHAPE_POP_BETA_ALPHA,
                              "\\alpha",
                              "alpha", 1.0, 1.0e2, 1.0e-1,
                              NC_GALAXY_SHAPE_POP_BETA_DEFAULT_PARAMS_ABSTOL,
                              NC_GALAXY_SHAPE_POP_BETA_DEFAULT_ALPHA,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxyShapePopBeta:beta:
   *
   * The shape parameter $\beta$ of the Beta distribution of $r = |\chi_I|$.
   * Bounded to $\ge 1$ (see the class documentation).
   *
   */
  ncm_model_class_set_sparam (model_class,
                              NC_GALAXY_SHAPE_POP_BETA_BETA,
                              "\\beta",
                              "beta", 1.0, 1.0e2, 1.0e-1,
                              NC_GALAXY_SHAPE_POP_BETA_DEFAULT_PARAMS_ABSTOL,
                              NC_GALAXY_SHAPE_POP_BETA_DEFAULT_BETA,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  gsp_class->data_init            = &_nc_galaxy_shape_pop_beta_data_init;
  gsp_class->prepare              = &_nc_galaxy_shape_pop_beta_prepare;
  gsp_class->eval_p               = &_nc_galaxy_shape_pop_beta_eval_p;
  gsp_class->eval_p_array         = &_nc_galaxy_shape_pop_beta_eval_p_array;
  gsp_class->gen                  = &_nc_galaxy_shape_pop_beta_gen;
  gsp_class->e_rms                = &_nc_galaxy_shape_pop_beta_e_rms;
  gsp_class->exponent_at_origin   = &_nc_galaxy_shape_pop_beta_exponent_at_origin;
  gsp_class->eval_p_rho2_g_series = &_nc_galaxy_shape_pop_beta_eval_p_rho2_g_series;
  gsp_class->moment_2k            = &_nc_galaxy_shape_pop_beta_moment_2k;
}

#define VECTOR (NCM_MODEL (gsp))
#define ALPHA  (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SHAPE_POP_BETA_ALPHA))
#define BETA   (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SHAPE_POP_BETA_BETA))

static void
_nc_galaxy_shape_pop_beta_ldata_noop (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_nc_galaxy_shape_pop_beta_ldata_required_columns (NcGalaxyShapePopData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_pop_beta_data_init (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  NcGalaxyShapePopBetaLData *ldata = g_new0 (NcGalaxyShapePopBetaLData, 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &g_free;
  data->ldata_read_row         = &_nc_galaxy_shape_pop_beta_ldata_noop;
  data->ldata_write_row        = &_nc_galaxy_shape_pop_beta_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_pop_beta_ldata_required_columns;
  data->ldata_get_mode_r       = &_nc_galaxy_shape_pop_beta_ldata_get_mode_r;
}

static gdouble
_nc_galaxy_shape_pop_beta_ldata_get_mode_r (NcGalaxyShapePopData *data)
{
  NcGalaxyShapePopBetaLData *ldata = (NcGalaxyShapePopBetaLData *) data->ldata;

  /* The argmax of P_pop(r) itself. */
  if ((ldata->alpha > 1.0) && (ldata->beta > 1.0))
    return (ldata->alpha - 1.0) / (ldata->alpha + ldata->beta - 2.0);
  else
    return 0.0;
}

static void
_nc_galaxy_shape_pop_beta_prepare (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  NcGalaxyShapePopBetaLData *ldata = (NcGalaxyShapePopBetaLData *) data->ldata;
  const gdouble alpha              = ALPHA;
  const gdouble beta               = BETA;

  /* lnnorm_x: x-space normalization (-ln(2 B(alpha,beta)), the -ln(2)
   * folding in the r=sqrt(x) change-of-variables Jacobian), used only by
   * eval_p_rho2_g_series below. lnnorm_r: -ln B(alpha,beta), the class
   * doc's own P_pop(r) formula's normalization, no Jacobian at all. */
  ldata->lnnorm_x = -gsl_sf_lnbeta (alpha, beta) - M_LN2;
  ldata->lnnorm_r = -gsl_sf_lnbeta (alpha, beta);
  ldata->alpha    = alpha;
  ldata->beta     = beta;
}

static gdouble
_nc_galaxy_shape_pop_beta_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble r)
{
  NcGalaxyShapePopBetaLData *ldata = (NcGalaxyShapePopBetaLData *) data->ldata;

  /* P_pop(r) = r^(alpha-1) * (1-r)^(beta-1) / B(alpha,beta). Log-space eval
   * avoids overflow/NaN for concentrated populations (large alpha, beta). */
  return exp ((ldata->alpha - 1.0) * log (r) + (ldata->beta - 1.0) * log1p (-r) + ldata->lnnorm_r);
}

/* Batched form of eval_p() above: alpha/beta/lnnorm are invariant across
 * the whole call, so this is a straight-line loop with no per-element
 * vfunc dispatch (#NcGalaxyShapeFactorFixedQuad's performance-critical path).
 * A three-pass
 * split (separate log/log1p/exp loops, to let the compiler vectorize via
 * libmvec) was tried and measured no benefit on this toolchain -- GCC's
 * vectorizer found no vectype for these calls regardless -- so this stays
 * the simpler single-pass form. */
static void
_nc_galaxy_shape_pop_beta_eval_p_array (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const GArray *r, GArray **p)
{
  NcGalaxyShapePopBetaLData *ldata = (NcGalaxyShapePopBetaLData *) data->ldata;
  const gdouble alpha_m1           = ldata->alpha - 1.0;
  const gdouble beta_m1            = ldata->beta - 1.0;
  const gdouble lnnorm_r           = ldata->lnnorm_r;
  const guint n                    = r->len;
  const gdouble * restrict r_data;
  gdouble * restrict p_data;
  guint i;

  if (*p == NULL)
    *p = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n);

  g_array_set_size (*p, n);

  r_data = (const gdouble *) r->data;
  p_data = (gdouble *) (*p)->data;

  #pragma omp simd

  for (i = 0; i < n; i++)
    p_data[i] = exp (alpha_m1 * log (r_data[i]) + beta_m1 * log1p (-r_data[i]) + lnnorm_r);
}

/* Composes the g-series of x=|chi_I(chi_L,g)|^2 with this population's x-space
 * density P(x) ~ x^(alpha/2-1)*(1-sqrt(x))^(beta-1). Note that this is the
 * x-space density (see ldata->lnnorm_x's own comment), NOT what eval_p(r)
 * returns.
 *
 * sqrt_x is @x_series^0.5 via ncm_laurent_series_tps_pow(), the same
 * generalized-binomial recursion the fractional exponents below rely on.
 * 1-sqrt(x) is built by scaling sqrt_x by -1 and adding one to its order-0
 * term, a single harmonic-0 entry, since x(g) and therefore sqrt_x are real
 * at every order. See
 * tests/c/nc/lss/galaxy/test_nc_galaxy_shape_pop_series.c for the cross-check
 * against eval_p(sqrt(x))/(2*sqrt(x)) at various g.
 *
 * sqrt(x) has a branch point at x=0 even at integer alpha, which limits this
 * series' radius of convergence in g to wherever x(g) first reaches 0 in the
 * complex g-plane. That radius becomes unusably small for any alpha<2, where
 * P(x) has a pole at x=0 rather than only a branch point.
 * #NcGalaxyShapeFactorSeriesLensed callers should keep alpha>=2; see
 * docs/theory/wl_shape_factor_history.md. */
static void
_nc_galaxy_shape_pop_beta_eval_p_rho2_g_series (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data,
                                                const NcmLaurentSeriesTPS *x_series, NcmLaurentSeriesTPS *out)
{
  NcGalaxyShapePopBetaLData *ldata      = (NcGalaxyShapePopBetaLData *) data->ldata;
  const guint order                     = ncm_laurent_series_tps_order (x_series);
  NcmLaurentSeriesTPS *sqrt_x           = ncm_laurent_series_tps_new (order);
  NcmLaurentSeriesTPS *one_minus_sqrt_x = ncm_laurent_series_tps_new (order);
  NcmLaurentSeriesTPS *num_pow          = ncm_laurent_series_tps_new (order);
  NcmLaurentSeriesTPS *den_pow          = ncm_laurent_series_tps_new (order);
  NcmLaurentSeriesTPS *unnormalized     = ncm_laurent_series_tps_new (order);
  NcmLaurentSeries *slot0;

  ncm_laurent_series_tps_pow (sqrt_x, x_series, 0.5);

  ncm_laurent_series_tps_scale (one_minus_sqrt_x, sqrt_x, -1.0);
  slot0 = ncm_laurent_series_tps_get (one_minus_sqrt_x, 0);
  g_assert_cmpint (ncm_laurent_series_get_hmin (slot0), ==, 0);
  g_assert_cmpint (ncm_laurent_series_get_hmax (slot0), ==, 0);
  ncm_laurent_series_set (slot0, 0, ncm_laurent_series_get (slot0, 0) + 1.0);

  ncm_laurent_series_tps_pow (num_pow, x_series, 0.5 * ldata->alpha - 1.0);
  ncm_laurent_series_tps_pow (den_pow, one_minus_sqrt_x, ldata->beta - 1.0);
  ncm_laurent_series_tps_conv (unnormalized, num_pow, den_pow);
  ncm_laurent_series_tps_scale (out, unnormalized, exp (ldata->lnnorm_x));

  ncm_laurent_series_tps_unref (sqrt_x);
  ncm_laurent_series_tps_unref (one_minus_sqrt_x);
  ncm_laurent_series_tps_unref (num_pow);
  ncm_laurent_series_tps_unref (den_pow);
  ncm_laurent_series_tps_unref (unnormalized);
}

static void
_nc_galaxy_shape_pop_beta_gen (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2)
{
  NcGalaxyShapePopBetaLData *ldata = (NcGalaxyShapePopBetaLData *) data->ldata;

  /* r=|chi_I| is the Beta-distributed variable, sampled directly -- no
   * sqrt() needed (see the class doc). */
  const gdouble r     = ncm_rng_beta_gen (rng, ldata->alpha, ldata->beta);
  const gdouble theta = ncm_rng_uniform_gen (rng, 0.0, 2.0 * M_PI);

  *e_int_1 = r * cos (theta);
  *e_int_2 = r * sin (theta);
}

static gdouble
_nc_galaxy_shape_pop_beta_e_rms (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  return nc_galaxy_shape_pop_beta_get_e_rms (NC_GALAXY_SHAPE_POP_BETA (gsp));
}

static gdouble
_nc_galaxy_shape_pop_beta_exponent_at_origin (NcGalaxyShapePop *gsp)
{
  const gdouble alpha = ALPHA;

  return (alpha - 1.0);
}

/* M_2k = <r^2k> = B(alpha+2k,beta)/B(alpha,beta) = prod_{j=0}^{2k-1}
 * (alpha+j)/(alpha+beta+j), the r=|chi_I| reparametrization's own moment
 * (note alpha+2k, not alpha+k -- see the class docs' x=r^2 vs r discussion).
 * The k=1 case reduces to nc_galaxy_shape_pop_beta_get_e_rms()'s own
 * alpha(alpha+1)/[(alpha+beta)(alpha+beta+1)]. Elementary product, no
 * gsl_sf_lnbeta() call needed for this integer-order case -- verified to
 * machine precision. Uses the already-resolved @data->ldata, matching
 * eval_p()'s own read contract. */
static gdouble
_nc_galaxy_shape_pop_beta_moment_2k (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const guint k)
{
  if (k == 0)
  {
    return 1.0;
  }
  else
  {
    NcGalaxyShapePopBetaLData *ldata = (NcGalaxyShapePopBetaLData *) data->ldata;
    const gdouble alpha              = ldata->alpha;
    const gdouble beta               = ldata->beta;
    gdouble moment                   = 1.0;
    guint j;

    for (j = 0; j < 2 * k; j++)
      moment *= (alpha + j) / (alpha + beta + j);

    return moment;
  }
}

/**
 * nc_galaxy_shape_pop_beta_new:
 *
 * Creates a new #NcGalaxyShapePopBeta.
 *
 * Returns: (transfer full): a new #NcGalaxyShapePopBeta.
 */
NcGalaxyShapePopBeta *
nc_galaxy_shape_pop_beta_new (void)
{
  NcGalaxyShapePopBeta *gspb = g_object_new (NC_TYPE_GALAXY_SHAPE_POP_BETA,
                                             NULL);

  return gspb;
}

/**
 * nc_galaxy_shape_pop_beta_ref:
 * @gspb: a #NcGalaxyShapePopBeta
 *
 * Increases the reference count of @gspb by one.
 *
 * Returns: (transfer full): @gspb.
 */
NcGalaxyShapePopBeta *
nc_galaxy_shape_pop_beta_ref (NcGalaxyShapePopBeta *gspb)
{
  return g_object_ref (gspb);
}

/**
 * nc_galaxy_shape_pop_beta_free:
 * @gspb: a #NcGalaxyShapePopBeta
 *
 * Decreases the reference count of @gspb by one.
 *
 */
void
nc_galaxy_shape_pop_beta_free (NcGalaxyShapePopBeta *gspb)
{
  g_object_unref (gspb);
}

/**
 * nc_galaxy_shape_pop_beta_clear:
 * @gspb: a #NcGalaxyShapePopBeta
 *
 * Decreases the reference count of *@gspb by one, and sets the pointer *@gspb to
 * NULL.
 *
 */
void
nc_galaxy_shape_pop_beta_clear (NcGalaxyShapePopBeta **gspb)
{
  g_clear_object (gspb);
}

/**
 * nc_galaxy_shape_pop_beta_get_e_rms:
 * @gspb: a #NcGalaxyShapePopBeta
 *
 * The induced $e_\mathrm{rms}=\sqrt{\langle x\rangle/2}$, with $\langle
 * x\rangle=\langle r^2\rangle=\alpha(\alpha+1)/[(\alpha+\beta)(\alpha+\beta+1)]$
 * for $r=|\chi_I|\sim\mathrm{Beta}(\alpha,\beta)$ -- a reporting quantity,
 * not itself a model parameter (see the class documentation).
 *
 * Returns: $e_\mathrm{rms}$.
 */
gdouble
nc_galaxy_shape_pop_beta_get_e_rms (NcGalaxyShapePopBeta *gspb)
{
  NcGalaxyShapePop *gsp = NC_GALAXY_SHAPE_POP (gspb);
  const gdouble alpha   = ALPHA;
  const gdouble beta    = BETA;
  const gdouble mean_r2 = alpha * (alpha + 1.0) / ((alpha + beta) * (alpha + beta + 1.0));

  return sqrt (0.5 * mean_r2);
}

/**
 * nc_galaxy_shape_pop_beta_get_mode:
 * @gspb: a #NcGalaxyShapePopBeta
 *
 * The induced $\mathrm{mode}(r) = (\alpha-1)/(\alpha+\beta-2)$, $r=|\chi_I|$
 * -- the "peak/typical ellipticity" number directly comparable to
 * literature fits (see the class documentation), literally the argmax of
 * $P_\mathrm{pop}(r)$ itself (see
 * _nc_galaxy_shape_pop_beta_ldata_get_mode_r()) -- a reporting quantity,
 * not itself a model parameter.
 *
 * Returns: $\mathrm{mode}(r)$.
 */
gdouble
nc_galaxy_shape_pop_beta_get_mode (NcGalaxyShapePopBeta *gspb)
{
  NcGalaxyShapePop *gsp = NC_GALAXY_SHAPE_POP (gspb);
  const gdouble alpha   = ALPHA;
  const gdouble beta    = BETA;

  if ((alpha > 1.0) && (beta > 1.0))
    return (alpha - 1.0) / (alpha + beta - 2.0);
  else
    return 0.0;
}

static void
_nc_galaxy_shape_pop_beta_flist_e_rms (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  NcGalaxyShapePopBeta *gspb = NC_GALAXY_SHAPE_POP_BETA (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));

  g_assert (NC_IS_GALAXY_SHAPE_POP_BETA (gspb));

  res[0] = nc_galaxy_shape_pop_beta_get_e_rms (gspb);
}

static void
_nc_galaxy_shape_pop_beta_flist_mode (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  NcGalaxyShapePopBeta *gspb = NC_GALAXY_SHAPE_POP_BETA (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));

  g_assert (NC_IS_GALAXY_SHAPE_POP_BETA (gspb));

  res[0] = nc_galaxy_shape_pop_beta_get_mode (gspb);
}

/* Called once from ncm_cfg_register_functions(), same convention as other
 * models' derived NcmMSetFuncList entries. */
void
_nc_galaxy_shape_pop_beta_register_functions (void)
{
  ncm_mset_func_list_register ("e_rms", "e_\\mathrm{rms}", "NcGalaxyShapePopBeta",
                               "RMS intrinsic ellipticity sqrt(<x>/2)", G_TYPE_NONE,
                               _nc_galaxy_shape_pop_beta_flist_e_rms, 0, 1);
  ncm_mset_func_list_register ("mode", "\\mathrm{mode}(r)", "NcGalaxyShapePopBeta",
                               "Mode of r = |chi_I|", G_TYPE_NONE,
                               _nc_galaxy_shape_pop_beta_flist_mode, 0, 1);
}

