/***************************************************************************
 *            nc_galaxy_shape_factor_cgf.c
 *
 *  Wed Jul 15 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_cgf.c
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
 * NcGalaxyShapeFactorCGF:
 *
 * CGF-expansion evaluation of the intrinsic-ellipticity marginal.
 *
 * Approximates
 * $$P(\epsilon_\mathrm{obs} \mid g) = \int_{|\chi_I|<1} \mathrm{d}^2\chi_I\,
 *   P_\mathrm{pop}(\chi_I)\, N_2\big(\epsilon_\mathrm{obs} - f_g(\chi_I);
 *   \sigma_\mathrm{noise}^2\big)$$
 * by a single Gaussian in $\epsilon_\mathrm{obs}$ obtained from a
 * cumulant-generating-function (moment) expansion of the pushforward of the
 * intrinsic distribution through the forward shear map $S(g,\cdot) \equiv f_g$.
 * Writing the analytic low-order response moments of $S$ around $\chi_I=0$ —
 * the value $S_0 = S(g,0)$, the (real $2\times2$) Jacobian $A = \partial S/
 * \partial\chi_I|_0$, and the Laplacian $\nabla^2 S|_0$ — the observed
 * ellipticity is modelled as $\mathcal{N}(\mu, C)$ with
 * $$\mu = S_0 + \tfrac12 V\,\nabla^2 S|_0, \qquad
 *   C = V\,A A^{\mathsf T} + \sigma_\mathrm{noise}^2\, I_2,$$
 * where $V$ is the per-component variance of the intrinsic ellipticity.
 *
 * Unlike #NcGalaxyShapeFactorVarAdd (which pulls $\epsilon_\mathrm{obs}$ back
 * through the inverse map and adds scalar variances), CGF keeps the map's
 * curvature: the mean picks up the Laplacian correction and the covariance is
 * the full $A A^{\mathsf T}$ pushforward of the intrinsic variance, not an
 * isotropic sum. It is therefore more accurate than VarAdd, while remaining a
 * single closed-form arithmetic evaluation — cheaper than
 * #NcGalaxyShapeFactorLaplace, which needs a per-galaxy Newton mode search and
 * a Hessian.
 *
 * The variance $V$ used here is the *truncated* per-component second moment of
 * the intrinsic ellipticity, i.e. $V = e_\mathrm{rms}^2$ with $e_\mathrm{rms}$
 * from nc_galaxy_shape_pop_e_rms(), NOT the square of the untruncated Gaussian
 * width nc_galaxy_shape_pop_get_sigma() (which always exceeds it and would bias
 * the recovered shear high).
 *
 * The moment expansion keeps only the intrinsic second moment, so it is a
 * Gaussian-population method: like #NcGalaxyShapeFactorVarAdd it requires the
 * population resolved from the #NcmMSet to support
 * nc_galaxy_shape_pop_get_sigma() (currently #NcGalaxyShapePopGauss or
 * #NcGalaxyShapePopGaussLocal, Global or per-galaxy).
 *
 * The approximation is exact in the doubly-linear ($g\to0$ or $V\to0$) limit
 * and degrades as the reduced shear or the intrinsic width grows (the dropped
 * higher moments and higher map derivatives then matter). See the accuracy
 * envelope measured in tests/python/nc/lss/galaxy/test_galaxy_shape_factor_cgf.py.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_cgf.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

/*
 * Analytic response moments of the forward shear map S(g,chi_I) at chi_I=0.
 * S(g,.) is exactly nc_wl_ellipticity_apply_shear_*(). Derived in closed form
 * per ellipticity convention:
 *
 * TRACE_DET (epsilon convention), |g|<=1: S is holomorphic in chi_I (no
 * conj(chi_I) dependence), so S(g,0)=g, the Jacobian is the real scalar
 * (1-|g|^2) times the identity, and S is harmonic, so its Laplacian vanishes.
 * |g|>1: S depends on chi_I only through conj(chi_I) (anti-holomorphic),
 * S(g,0)=g/|g|^2, Jacobian is a reflection built from beta=(|g|^2-1)/conj(g)^2,
 * and the Laplacian still vanishes.
 *
 * TRACE (chi/distortion convention): S depends on both chi_I and conj(chi_I),
 * so neither holomorphicity shortcut applies; Jacobian and Laplacian were
 * obtained by direct Wirtinger-calculus differentiation of the Möbius-style
 * map.
 *
 * The |g|>1 (strong-lensing) branches are outside this project's shear range.
 */
typedef struct _NcGalaxyShapeFactorCGFMoments
{
  complex double S0;
  gdouble A11, A12, A21, A22;
  complex double lap_S;
} NcGalaxyShapeFactorCGFMoments;

struct _NcGalaxyShapeFactorCGF
{
  NcGalaxyShapeFactor parent_instance;
};

typedef struct _NcGalaxyShapeFactorCGFPrivate
{
  /* Convention specialization resolved once at construction (ellip-conv is
   * CONSTRUCT_ONLY), keeping the switch out of the per-evaluation path. */
  void (*response_moments) (complex double g, NcGalaxyShapeFactorCGFMoments *out);

  /* Population generation, refreshed in prepare(). The per-galaxy V cache
   * (see NcGalaxyShapeFactorCGFLData) is keyed on it. Held here rather than
   * re-read through nc_galaxy_shape_factor_get_pop_hash() per evaluation
   * because the eval hooks already load this private struct, making the
   * cache check two derefs and a compare -- cheaper than the vfunc dispatch
   * it replaces. The base sets its own pop_hash before invoking
   * klass->prepare(), so reading it here always sees the current value. */
  guint64 pop_hash;
} NcGalaxyShapeFactorCGFPrivate;

/*
 * Per-galaxy scratch: the intrinsic per-component variance V = e_rms^2.
 * V is constant across a galaxy's z-nodes, but eval_marginal() is invoked once
 * per node by nc_galaxy_shape_factor_eval_at_nodes(), so fetching it through
 * nc_galaxy_shape_pop_e_rms()'s vfunc every time cost ~1% of an MCMC run.
 * Invalidated on two independent axes: @pop_hash_seen catches a change in the
 * population MODEL (e.g. a fitted NcGalaxyShapePopGauss:sigma), while
 * ldata_read_row() clears @valid because a per-galaxy population
 * (NcGalaxyShapePopGaussLocal) takes e_rms from the catalog row -- data, which
 * moves with no model pkey bump at all.
 */
typedef struct _NcGalaxyShapeFactorCGFLData
{
  gdouble V;
  guint64 pop_hash_seen;
  gboolean valid;
} NcGalaxyShapeFactorCGFLData;

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorCGF, nc_galaxy_shape_factor_cgf, NC_TYPE_GALAXY_SHAPE_FACTOR)

static void
_nc_galaxy_shape_factor_cgf_response_moments_trace (complex double g, NcGalaxyShapeFactorCGFMoments *out)
{
  const gdouble abs_g2 = creal (g * conj (g));
  const gdouble D0     = 1.0 + abs_g2;

  /* One reciprocal, then multiplies: D0 appears as 1/D0, 1/D0^2 and 1/D0^3,
   * and each spelled-out division is a ~15-20 cycle divsd that does not
   * pipeline. This function is ~8.8% of an MCMC run's self time (it runs once
   * per galaxy per z-node), so the four divisions it used to issue were worth
   * removing. */
  const gdouble inv_D0    = 1.0 / D0;
  const gdouble inv_D0_2  = inv_D0 * inv_D0;
  const gdouble one_m_g2  = 1.0 - abs_g2;
  const gdouble kappa     = one_m_g2 * inv_D0_2;
  const complex double g2 = g * g;

  out->S0  = 2.0 * g * inv_D0;
  out->A11 = kappa * (1.0 - creal (g2));
  out->A12 = -kappa *cimag (g2);

  out->A21   = out->A12;
  out->A22   = kappa * (1.0 + creal (g2));
  out->lap_S = -4.0 * g * (one_m_g2 * one_m_g2) * (inv_D0_2 * inv_D0);
}

static void
_nc_galaxy_shape_factor_cgf_response_moments_trace_det (complex double g, NcGalaxyShapeFactorCGFMoments *out)
{
  const gdouble abs_g2 = creal (g * conj (g));

  if (abs_g2 <= 1.0)
  {
    out->S0    = g;
    out->A11   = 1.0 - abs_g2;
    out->A12   = 0.0;
    out->A21   = 0.0;
    out->A22   = 1.0 - abs_g2;
    out->lap_S = 0.0;
  }
  else
  {
    const complex double g_conj = conj (g);
    const complex double beta   = (abs_g2 - 1.0) / (g_conj * g_conj);

    out->S0    = g / abs_g2;
    out->A11   = creal (beta);
    out->A12   = cimag (beta);
    out->A21   = cimag (beta);
    out->A22   = -creal (beta);
    out->lap_S = 0.0;
  }
}

/*
 * Evaluate the 2D Gaussian N(mu, C) at e_o, where mu and C are built from the
 * response moments and the intrinsic variance V:
 *   mu = S0 + 0.5 * V * lap_S,   C = V * A A^T + sigma_noise^2 * I_2.
 */
static gdouble
_nc_galaxy_shape_factor_cgf_eval (const NcGalaxyShapeFactorCGFMoments *mom, complex double e_o, const gdouble V, const gdouble sigma_noise, const gboolean want_log)
{
  const gdouble sigma_noise2 = gsl_pow_2 (sigma_noise);
  const complex double mu    = mom->S0 + 0.5 * V * mom->lap_S;
  gdouble Cxx, Cxy, Cyy, det, dx, dy, chi2;

  /* C = V * A A^T + sigma_noise^2 * I_2 */
  Cxx = V * (gsl_pow_2 (mom->A11) + gsl_pow_2 (mom->A12)) + sigma_noise2;
  Cxy = V * (mom->A11 * mom->A21 + mom->A12 * mom->A22);
  Cyy = V * (gsl_pow_2 (mom->A21) + gsl_pow_2 (mom->A22)) + sigma_noise2;

  det = Cxx * Cyy - gsl_pow_2 (Cxy);
  dx  = creal (e_o) - creal (mu);
  dy  = cimag (e_o) - cimag (mu);

  /* C^-1 = (1/det) [[Cyy,-Cxy],[-Cxy,Cxx]] */
  chi2 = (Cyy * dx * dx - 2.0 * Cxy * dx * dy + Cxx * dy * dy) / det;

  /* The linear branch is the hot one -- FIXED_NODES evaluates through
   * eval_at_nodes(), i.e. want_log = FALSE, once per galaxy per z-node -- so it
   * avoids log(det) entirely: exp(-chi2/2) / (2 pi sqrt(det)) is the same
   * quantity as exp(-chi2/2 - log(2 pi) - log(det)/2) but trades a ~40-60 cycle
   * log for a ~15-20 cycle sqrt. (log(2 pi) is constant-folded by the compiler
   * in either spelling, so only log(det) was ever actually at stake.) */
  if (want_log)
    return -0.5 * chi2 - log (2.0 * M_PI) - 0.5 * log (det);
  else
    return exp (-0.5 * chi2) / (2.0 * M_PI * sqrt (det));
}

/* Cached V = e_rms^2; recomputed only when the population model generation
 * moves or a new catalog row was read into @data. */
static inline gdouble
_nc_galaxy_shape_factor_cgf_peek_V (NcGalaxyShapeFactorCGFPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorCGFLData *ldata = (NcGalaxyShapeFactorCGFLData *) data->ldata;

  if (G_UNLIKELY (!ldata->valid || (ldata->pop_hash_seen != self->pop_hash)))
  {
    ldata->V             = gsl_pow_2 (nc_galaxy_shape_pop_e_rms (pop, data->pop_data));
    ldata->pop_hash_seen = self->pop_hash;
    ldata->valid         = TRUE;
  }

  return ldata->V;
}

static void
nc_galaxy_shape_factor_cgf_init (NcGalaxyShapeFactorCGF *gsfcgf)
{
  NcGalaxyShapeFactorCGFPrivate * const self = nc_galaxy_shape_factor_cgf_get_instance_private (gsfcgf);

  self->response_moments = NULL;
  self->pop_hash         = 0;
}

static void
_nc_galaxy_shape_factor_cgf_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_cgf_parent_class)->constructed (object);
  {
    NcGalaxyShapeFactorCGF *gsfcgf             = NC_GALAXY_SHAPE_FACTOR_CGF (object);
    NcGalaxyShapeFactorCGFPrivate * const self = nc_galaxy_shape_factor_cgf_get_instance_private (gsfcgf);
    const NcGalaxyWLObsEllipConv ellip_conv    = nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (gsfcgf));

    switch (ellip_conv)
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
        self->response_moments = &_nc_galaxy_shape_factor_cgf_response_moments_trace;
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        self->response_moments = &_nc_galaxy_shape_factor_cgf_response_moments_trace_det;
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
    }
  }
}

static void
_nc_galaxy_shape_factor_cgf_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_cgf_parent_class)->finalize (object);
}

static void
_nc_galaxy_shape_factor_cgf_ldata_noop (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
}

/* A per-galaxy population reads e_rms straight from the catalog row, so a new
 * row invalidates the cached V without any model pkey moving. */
static void
_nc_galaxy_shape_factor_cgf_ldata_read_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyShapeFactorCGFLData *ldata = (NcGalaxyShapeFactorCGFLData *) data->ldata;

  ldata->valid = FALSE;
}

static void
_nc_galaxy_shape_factor_cgf_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_factor_cgf_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  /* Per-galaxy scratch holding the cached intrinsic variance V (see
   * NcGalaxyShapeFactorCGFLData); g_new0 leaves it invalid, so the first
   * evaluation populates it. */
  data->ldata                  = g_new0 (NcGalaxyShapeFactorCGFLData, 1);
  data->ldata_destroy          = &g_free;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_cgf_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_cgf_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_cgf_ldata_required_columns;
}

static void
_nc_galaxy_shape_factor_cgf_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
  NcGalaxyShapeFactorCGFPrivate * const self = nc_galaxy_shape_factor_cgf_get_instance_private (NC_GALAXY_SHAPE_FACTOR_CGF (gsf));

  /* Sigma-support is a property of the pop MODEL, identical for every galaxy
   * sharing it -- checked once here via a throwaway pop_data, not per-galaxy.
   * The CGF moment expansion is a Gaussian-population method, so it gates on
   * the same capability as VarAdd (even though V itself comes from e_rms). */
  NcGalaxyShapePop *pop              = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));
  NcGalaxyShapePopData *tmp_pop_data = nc_galaxy_shape_pop_data_new (pop);
  const gboolean has_sigma           = tmp_pop_data->ldata_get_sigma != NULL;

  nc_galaxy_shape_pop_data_unref (tmp_pop_data);

  if (!has_sigma)
    g_error ("NcGalaxyShapeFactorCGF: the CGF moment expansion requires a Gaussian population "
             "supporting nc_galaxy_shape_pop_get_sigma(), got %s.", G_OBJECT_TYPE_NAME (pop));

  /* The base sets its own pop_hash before calling this vfunc, so this snapshot
   * is current; the per-galaxy V caches are keyed on it. */
  self->pop_hash = nc_galaxy_shape_factor_get_pop_hash (gsf);
}

static gdouble
_nc_galaxy_shape_factor_cgf_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorCGFPrivate * const self = nc_galaxy_shape_factor_cgf_get_instance_private (NC_GALAXY_SHAPE_FACTOR_CGF (gsf));
  const complex double g                     = g_1 + I * g_2;
  const complex double e_o                   = epsilon_obs_1 + I * epsilon_obs_2;
  const gdouble V                            = _nc_galaxy_shape_factor_cgf_peek_V (self, pop, data);
  NcGalaxyShapeFactorCGFMoments mom;

  self->response_moments (g, &mom);

  return _nc_galaxy_shape_factor_cgf_eval (&mom, e_o, V, data->std_noise, FALSE);
}

static gdouble
_nc_galaxy_shape_factor_cgf_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorCGFPrivate * const self = nc_galaxy_shape_factor_cgf_get_instance_private (NC_GALAXY_SHAPE_FACTOR_CGF (gsf));
  const complex double g                     = g_1 + I * g_2;
  const complex double e_o                   = epsilon_obs_1 + I * epsilon_obs_2;
  const gdouble V                            = _nc_galaxy_shape_factor_cgf_peek_V (self, pop, data);
  NcGalaxyShapeFactorCGFMoments mom;

  self->response_moments (g, &mom);

  return _nc_galaxy_shape_factor_cgf_eval (&mom, e_o, V, data->std_noise, TRUE);
}

static void
nc_galaxy_shape_factor_cgf_class_init (NcGalaxyShapeFactorCGFClass *klass)
{
  NcGalaxyShapeFactorClass *gsf_class = NC_GALAXY_SHAPE_FACTOR_CLASS (klass);
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);

  object_class->constructed = &_nc_galaxy_shape_factor_cgf_constructed;
  object_class->finalize    = &_nc_galaxy_shape_factor_cgf_finalize;

  gsf_class->data_init        = &_nc_galaxy_shape_factor_cgf_data_init;
  gsf_class->prepare          = &_nc_galaxy_shape_factor_cgf_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_cgf_eval_marginal;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_cgf_eval_ln_marginal;
}

/**
 * nc_galaxy_shape_factor_cgf_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 *
 * Creates a new #NcGalaxyShapeFactorCGF.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorCGF.
 */
NcGalaxyShapeFactorCGF *
nc_galaxy_shape_factor_cgf_new (NcGalaxyWLObsEllipConv ellip_conv)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_CGF,
                       "ellip-conv", ellip_conv,
                       NULL);
}

/**
 * nc_galaxy_shape_factor_cgf_ref:
 * @gsfcgf: a #NcGalaxyShapeFactorCGF
 *
 * Increases the reference count of @gsfcgf by one.
 *
 * Returns: (transfer full): @gsfcgf.
 */
NcGalaxyShapeFactorCGF *
nc_galaxy_shape_factor_cgf_ref (NcGalaxyShapeFactorCGF *gsfcgf)
{
  return g_object_ref (gsfcgf);
}

/**
 * nc_galaxy_shape_factor_cgf_free:
 * @gsfcgf: a #NcGalaxyShapeFactorCGF
 *
 * Decreases the reference count of @gsfcgf by one.
 *
 */
void
nc_galaxy_shape_factor_cgf_free (NcGalaxyShapeFactorCGF *gsfcgf)
{
  g_object_unref (gsfcgf);
}

/**
 * nc_galaxy_shape_factor_cgf_clear:
 * @gsfcgf: a #NcGalaxyShapeFactorCGF
 *
 * Decreases the reference count of *@gsfcgf by one, and sets the pointer
 * *@gsfcgf to NULL.
 *
 */
void
nc_galaxy_shape_factor_cgf_clear (NcGalaxyShapeFactorCGF **gsfcgf)
{
  g_clear_object (gsfcgf);
}

