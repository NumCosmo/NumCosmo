/***************************************************************************
 *            nc_galaxy_shape_factor_quad.c
 *
 *  Fri Jul 3 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_quad.c
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
 * NcGalaxyShapeFactorQuad:
 *
 * Exact quadrature evaluation of the intrinsic-ellipticity marginal.
 *
 * Computes
 * $$P(\epsilon_\mathrm{obs} \mid g) = \int_{|\chi_I|<1} \mathrm{d}^2\chi_I\,
 *   P_\mathrm{pop}(\chi_I)\, N_2\big(\epsilon_\mathrm{obs} - f_g(\chi_I);
 *   \sigma_\mathrm{noise}^2\big)$$
 * exactly, with no linearization of the shear map and no truncation of the
 * intrinsic population to an untruncated Gaussian (contrast
 * #NcGalaxyShapeFactorVarAdd, which approximates both).
 *
 * The substitution variable is the LENSED (noiseless, pre-noise) ellipticity
 * $\chi_L=f_g(\chi_I)$, not the intrinsic $\chi_I$: change variables in the
 * integral above to $\chi_L$,
 * $$P(\epsilon_\mathrm{obs} \mid g) = \int_{|\chi_L|<1} \mathrm{d}^2\chi_L\,
 *   P_\mathrm{pop}\big(f_g^{-1}(\chi_L)\big)\,
 *   \left|\det J_{f_g^{-1}}(\chi_L)\right|\,
 *   N_2\big(\epsilon_\mathrm{obs} - \chi_L; \sigma_\mathrm{noise}^2\big),$$
 * where $\det J_{f_g^{-1}}$ is supplied by the existing
 * nc_galaxy_shape_factor_lndet_jac() machinery (already used by `VarAdd`), so
 * no new map algebra is needed. $\chi_L$ is turned into a plane integral by
 * the usual bijection,
 * $$\chi_L(u,v) = (u,v)/\sqrt{1+u^2+v^2}, \qquad
 *   \left|\det\frac{\partial\chi_L}{\partial(u,v)}\right| = \frac{1}{(1+u^2+v^2)^2}.$$
 * The point of using $\chi_L$ rather than $\chi_I$ is that the noise kernel
 * $N_2(\epsilon_\mathrm{obs}-\chi_L;\sigma_\mathrm{noise}^2)$ becomes an
 * ordinary Gaussian of EXACTLY known width $\sigma_\mathrm{noise}$ centered at
 * EXACTLY $\epsilon_\mathrm{obs}$ -- no shear map, no approximation. Because
 * the population is evaluated at $\chi_I=f_g^{-1}(\chi_L)$ instead of
 * directly at the substitution point, nc_galaxy_shape_pop_eval_p() is used
 * (not nc_galaxy_shape_pop_eval_p_rho2(), whose $\rho^2$ contract is
 * specifically $\rho^2=u^2+v^2$ for the point being substituted for, which is
 * $\chi_L$ here, not $\chi_I$); $x_I=\lvert\chi_I\rvert^2$ is computed
 * directly from the complex division (a sum of two squares, not itself
 * subtraction-sensitive), only losing the extra conditioning
 * nc_galaxy_shape_pop_eval_p_rho2() offers for the $(1-x)^a$ factor right at
 * $x\to1^-$ (relevant mainly for #NcGalaxyShapePopBeta at the disc boundary).
 *
 * The plane integral is evaluated by the Divonne cubature algorithm from the
 * Cuba library (ncm_integrate_2dim_divonne()), over a FIXED box of
 * half-width @bound centered on the exact $(u,v)$ preimage of
 * $\epsilon_\mathrm{obs}$ (no shear map needed: $\chi_L$'s noise kernel peaks
 * exactly there), seeded with two explicit peak hints: that same point (the
 * noise peak, exact) and the preimage of $f_g(0)$ (where $\chi_L$ sits when
 * $\chi_I=0$, i.e. the population's peak location under the common
 * assumption that $P_\mathrm{pop}$ is radially symmetric about $\chi_I=0$,
 * e.g. #NcGalaxyShapePopGauss / #NcGalaxyShapePopGaussLocal -- for a
 * population peaked elsewhere, e.g. #NcGalaxyShapePopBeta with $\mu$ away
 * from 0, this hint is only approximate, but Divonne's own stratified search
 * still reliably finds the true peak from there in every case tested, see
 * below). $f_g(0)$ is convention-dependent -- exactly $g$ in the TRACE_DET
 * (ellipticity) convention, but $2g/(1+\lvert g\rvert^2)$ in the TRACE
 * (distortion) convention -- so it is computed via the actual forward map
 * rather than hardcoded as $g$, since hardcoding it would miss the peak for
 * narrow, off-center TRACE-convention populations.
 *
 * Divonne's explicit peak hints matter: a fixed-degree base rule with no
 * explicit hints (e.g. Cuba's Cuhre) has no mechanism to notice an isolated
 * feature it never samples near, and can return a confidently wrong result
 * (off by orders of magnitude, at tight reltol) whenever the population is
 * narrow relative to the integration box and off-center. Because Divonne is
 * told where to look, @bound can be a single generous fixed constant:
 * shrinking it does not help an already-hint-found narrow feature, but does
 * progressively truncate a broad population's disc-spanning support (a few
 * units already covers the disc for any population). Verified against an
 * independent scipy reference across populations from $\sigma=0.30$ down to
 * $\sigma=0.001$, five different $g$/$\epsilon_\mathrm{obs}$/
 * $\sigma_\mathrm{noise}$ geometries (including near the disc boundary), and
 * #NcGalaxyShapePopBeta concentrations up to $\nu=10^5$: every case matches
 * to $\sim\!10^{-6}$ relative accuracy or better. See
 * docs/theory/wl_shape_factor_history.md for earlier implementations of
 * this class that were tried and rejected.
 *
 * The integrand clamps any non-finite evaluation to zero before returning:
 * Cuba can segfault outright on a NaN/Inf sample (reproduced directly), and
 * some populations have a genuine, mathematically correct divergence at a
 * disc point (e.g. #NcGalaxyShapePopBeta with $\alpha<1$ diverges at $x=0$,
 * which can coincide with a peak hint).
 *
 * This scheme is exact but substantially more expensive per evaluation than
 * `VarAdd` (a full 2D cubature vs. one closed-form expression): Divonne with
 * two hints typically costs $\sim\!100\,\mathrm{ms}$ per evaluation, rising to several
 * seconds for extremely concentrated populations. It is meant as an accuracy
 * reference / fallback for regimes where the variance-add approximation is
 * not trusted, not as a routine replacement in large-catalog likelihood
 * evaluations. Cuba's own fork()-based internal parallelism is disabled
 * globally by ncm_cfg_init() (`cubacores(0, 0)`), so this is safe to call
 * concurrently from multiple OpenMP threads.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_quad.h"
#include "nc/lss/wl/nc_wl_ellipticity.h"
#include "ncm/integration/ncm_integrate.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

#define NC_GALAXY_SHAPE_FACTOR_QUAD_DEFAULT_BOUND  8.0
#define NC_GALAXY_SHAPE_FACTOR_QUAD_DEFAULT_RELTOL 1.0e-7

struct _NcGalaxyShapeFactorQuad
{
  NcGalaxyShapeFactor parent_instance;
};

typedef struct _NcGalaxyShapeFactorQuadPrivate
{
  /* Convention specialization resolved once at construction (ellip-conv is
   * CONSTRUCT_ONLY), keeping the switch out of the per-evaluation path.
   * apply_shear_inv and lndet_jac are evaluated at chi_L (the substitution
   * point) inside the integrand; no shear-map call is needed for
   * re-centering the box (see the class doc). apply_shear (forward) is used
   * only once per evaluation, to place the population-peak hint at the
   * convention-correct f_g(0) -- this is g exactly for TRACE_DET, but
   * 2g/(1+|g|^2) for TRACE, so it cannot be hardcoded as g. */
  complex double (*apply_shear) (complex double g, complex double chi_i);
  complex double (*apply_shear_inv) (complex double g, complex double chi_L);

  gdouble (*lndet_jac) (complex double g, complex double chi_L);

  gdouble bound;
  gdouble reltol;
} NcGalaxyShapeFactorQuadPrivate;

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorQuad, nc_galaxy_shape_factor_quad, NC_TYPE_GALAXY_SHAPE_FACTOR)

enum
{
  PROP_0,
  PROP_BOUND,
  PROP_RELTOL,
  PROP_LEN,
};

/* Divonne integrand ------------------------------------------------------ */

typedef struct _QuadIntArg
{
  NcGalaxyShapePop *pop;
  NcGalaxyShapePopData *pop_data;
  complex double (*apply_shear_inv) (complex double g, complex double chi_L);

  gdouble (*lndet_jac) (complex double g, complex double chi_L);

  complex double g;
  complex double eps_obs;
  gdouble std_noise;
} QuadIntArg;

static gdouble
_nc_galaxy_shape_factor_quad_integrand (gdouble u, gdouble v, gpointer userdata)
{
  QuadIntArg * const arg     = (QuadIntArg *) userdata;
  const gdouble rho2         = u * u + v * v;
  const gdouble s2           = 1.0 + rho2;
  const gdouble jac          = 1.0 / gsl_pow_2 (s2);
  const complex double chi_L = (u + I * v) / sqrt (s2);
  const complex double chi_i = arg->apply_shear_inv (arg->g, chi_L);
  const gdouble x_i          = gsl_pow_2 (creal (chi_i)) + gsl_pow_2 (cimag (chi_i));
  const gdouble p_pop        = nc_galaxy_shape_pop_eval_p (arg->pop, arg->pop_data, x_i) / M_PI;
  const gdouble jac_inv      = exp (arg->lndet_jac (arg->g, chi_L));
  const complex double delta = arg->eps_obs - chi_L;
  const gdouble noise_var    = gsl_pow_2 (arg->std_noise);
  const gdouble noise        = exp (-(gsl_pow_2 (creal (delta)) + gsl_pow_2 (cimag (delta))) / (2.0 * noise_var)) / (2.0 * M_PI * noise_var);
  const gdouble ret          = p_pop * jac_inv * noise * jac;

  /* Some populations have a genuine divergence at a disc point (e.g.
   * #NcGalaxyShapePopBeta with alpha<1 at x_i=0, which can coincide with a
   * peak hint); Cuba can segfault outright on a non-finite sample
   * (reproduced directly), so clamp rather than propagate. */
  return isfinite (ret) ? ret : 0.0;
}

static void
_nc_galaxy_shape_factor_quad_ldata_noop (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_nc_galaxy_shape_factor_quad_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_factor_quad_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  /* Divonne needs no persistent per-galaxy state: each evaluation is a
   * single, self-contained call to ncm_integrate_2dim_divonne(). */
  data->ldata                  = NULL;
  data->ldata_destroy          = &g_free;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_quad_ldata_noop;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_quad_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_quad_ldata_required_columns;
}

static void
_nc_galaxy_shape_factor_quad_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
}

/* (u,v) preimage of a disc point chi under the (u,v)->chi bijection, clamped
 * away from |chi|=1 (the disc boundary, where the preimage diverges). */
static void
_nc_galaxy_shape_factor_quad_preimage (complex double chi, gdouble *u, gdouble *v)
{
  const gdouble abs_chi         = cabs (chi);
  const gdouble abs_chi_clamped = MIN (abs_chi, 0.999);
  const complex double centered = (abs_chi > 0.0) ? chi * (abs_chi_clamped / abs_chi) : chi;
  const gdouble scale           = 1.0 / sqrt (1.0 - gsl_pow_2 (abs_chi_clamped));

  *u = creal (centered) * scale;
  *v = cimag (centered) * scale;
}

/* Refines theta, the angle of a candidate chi_I=rho_mode*exp(I*theta) on the
 * ring where a non-radially-symmetric population peaks (rho_mode>0, e.g.
 * NcGalaxyShapePopBeta with mu away from 0), so that its image f_g(chi_I)
 * lands as close as possible to eps_obs. The population density is exactly
 * flat along this ring by construction (it depends only on x=|chi_I|^2), so
 * this is the whole optimization: no population term to weigh in. A fixed
 * handful of Newton steps on d/dtheta |eps_obs-f_g(rho_mode e^{i theta})|^2,
 * with central-difference derivatives through the already-dispatched forward
 * map (convention-agnostic, no per-convention algebra duplicated), seeded at
 * the phase of the naive noiseless intrinsic preimage -- verified numerically
 * to already sit within a few hundredths of a radian of the true optimum,
 * with 3 Newton steps then closing the rest of the gap to floating-point
 * precision. This runs once per _eval() call (the hint, not the integrand),
 * so a handful of extra apply_shear() calls is negligible next to Divonne's
 * own ~100ms/evaluation cost. */
static gdouble
_nc_galaxy_shape_factor_quad_refine_theta (complex double ( *apply_shear ) (complex double, complex double),
                                           const complex double g, const gdouble rho_mode,
                                           const complex double eps_obs, const gdouble theta0)
{
  const gdouble h = 1.0e-4;
  gdouble theta   = theta0;
  guint i;

  for (i = 0; i < 3; i++)
  {
    const complex double delta_p = eps_obs - apply_shear (g, rho_mode * cexp (I * (theta + h)));
    const complex double delta_0 = eps_obs - apply_shear (g, rho_mode * cexp (I * theta));
    const complex double delta_m = eps_obs - apply_shear (g, rho_mode * cexp (I * (theta - h)));
    const gdouble fp             = gsl_pow_2 (creal (delta_p)) + gsl_pow_2 (cimag (delta_p));
    const gdouble f0             = gsl_pow_2 (creal (delta_0)) + gsl_pow_2 (cimag (delta_0));
    const gdouble fm             = gsl_pow_2 (creal (delta_m)) + gsl_pow_2 (cimag (delta_m));
    const gdouble grad           = (fp - fm) / (2.0 * h);
    const gdouble hess           = (fp - 2.0 * f0 + fm) / gsl_pow_2 (h);

    if (hess <= 0.0)
      break;  /* not a local min along this direction; keep the current estimate */

    theta -= grad / hess;
  }

  return theta;
}

static gdouble
_nc_galaxy_shape_factor_quad_eval (NcGalaxyShapeFactorQuad *gsfq, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorQuadPrivate * const self = nc_galaxy_shape_factor_quad_get_instance_private (gsfq);
  const complex double g                      = g_1 + I * g_2;
  const complex double eps_obs                = epsilon_obs_1 + I * epsilon_obs_2;
  const gdouble mode_x                        = nc_galaxy_shape_pop_get_mode_x (pop, data->pop_data);
  QuadIntArg arg;
  NcmIntegrand2dim integ = {&arg, &_nc_galaxy_shape_factor_quad_integrand};
  gdouble u0, v0, ug, vg, result, error;
  gdouble xgiven[6];
  gint ngiven;

  /* Two peak hints for Divonne, always: the exact preimage of eps_obs (the
   * noise kernel's peak, exact -- also the box center) and the preimage of
   * f_g(0) = self->apply_shear(g, 0) (the population's peak location under
   * the common radially-symmetric-about-0 assumption; only approximate
   * otherwise, but sufficient in every case tested -- see the class doc).
   * f_g(0) is convention-dependent (g exactly for TRACE_DET, 2g/(1+|g|^2)
   * for TRACE), so it must go through the actual forward map, not be
   * hardcoded as g. ncm_integrate_2dim_divonne() mutates xgiven in place
   * (normalizes it to the box), so it must be freshly populated on every
   * call; it is a local, per-call array here, so that is automatic. */
  _nc_galaxy_shape_factor_quad_preimage (eps_obs, &u0, &v0);
  _nc_galaxy_shape_factor_quad_preimage (self->apply_shear (g, 0.0), &ug, &vg);
  xgiven[0] = ug;
  xgiven[1] = vg;
  xgiven[2] = u0;
  xgiven[3] = v0;
  ngiven    = 2;

  /* Third hint, only for a population whose peak is not at chi_I=0 (e.g. a
   * concentrated NcGalaxyShapePopBeta with mu away from 0): the point on the
   * peak ring |chi_I|=rho_mode closest, after the forward map, to eps_obs
   * (see _nc_galaxy_shape_factor_quad_refine_theta() doc). */
  if (mode_x > 0.0)
  {
    const gdouble rho_mode           = sqrt (mode_x);
    const complex double chi_i_naive = self->apply_shear_inv (g, eps_obs);
    const gdouble theta0             = carg (chi_i_naive);
    const gdouble theta              = _nc_galaxy_shape_factor_quad_refine_theta (self->apply_shear, g, rho_mode, eps_obs, theta0);
    const complex double chi_i_hint  = rho_mode * cexp (I * theta);
    gdouble u2, v2;

    _nc_galaxy_shape_factor_quad_preimage (self->apply_shear (g, chi_i_hint), &u2, &v2);
    xgiven[4] = u2;
    xgiven[5] = v2;
    ngiven    = 3;
  }

  arg.pop             = pop;
  arg.pop_data        = data->pop_data;
  arg.apply_shear_inv = self->apply_shear_inv;
  arg.lndet_jac       = self->lndet_jac;
  arg.g               = g;
  arg.eps_obs         = eps_obs;
  arg.std_noise       = data->std_noise;

  ncm_integrate_2dim_divonne (&integ, u0 - self->bound, v0 - self->bound, u0 + self->bound, v0 + self->bound,
                              self->reltol, 0.0, ngiven, 2, xgiven, &result, &error);

  return result;
}

static gdouble
_nc_galaxy_shape_factor_quad_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return _nc_galaxy_shape_factor_quad_eval (NC_GALAXY_SHAPE_FACTOR_QUAD (gsf), pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2);
}

static gdouble
_nc_galaxy_shape_factor_quad_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return log (_nc_galaxy_shape_factor_quad_eval (NC_GALAXY_SHAPE_FACTOR_QUAD (gsf), pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2));
}

/* GObject boilerplate --------------------------------------------------- */

static void
nc_galaxy_shape_factor_quad_init (NcGalaxyShapeFactorQuad *gsfq)
{
  NcGalaxyShapeFactorQuadPrivate * const self = nc_galaxy_shape_factor_quad_get_instance_private (gsfq);

  self->apply_shear     = NULL;
  self->apply_shear_inv = NULL;
  self->lndet_jac       = NULL;
  self->bound           = NC_GALAXY_SHAPE_FACTOR_QUAD_DEFAULT_BOUND;
  self->reltol          = NC_GALAXY_SHAPE_FACTOR_QUAD_DEFAULT_RELTOL;
}

static void
_nc_galaxy_shape_factor_quad_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_quad_parent_class)->constructed (object);
  {
    NcGalaxyShapeFactorQuad *gsfq               = NC_GALAXY_SHAPE_FACTOR_QUAD (object);
    NcGalaxyShapeFactorQuadPrivate * const self = nc_galaxy_shape_factor_quad_get_instance_private (gsfq);
    const NcGalaxyWLObsEllipConv ellip_conv     = nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (gsfq));

    switch (ellip_conv)
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
        self->apply_shear     = &nc_wl_ellipticity_apply_shear_trace;
        self->apply_shear_inv = &nc_wl_ellipticity_apply_shear_inv_trace;
        self->lndet_jac       = &nc_wl_ellipticity_lndet_jac_trace;
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        self->apply_shear     = &nc_wl_ellipticity_apply_shear_trace_det;
        self->apply_shear_inv = &nc_wl_ellipticity_apply_shear_inv_trace_det;
        self->lndet_jac       = &nc_wl_ellipticity_lndet_jac_trace_det;
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
    }
  }
}

static void
_nc_galaxy_shape_factor_quad_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorQuad *gsfq = NC_GALAXY_SHAPE_FACTOR_QUAD (object);

  g_return_if_fail (NC_IS_GALAXY_SHAPE_FACTOR_QUAD (gsfq));

  switch (prop_id)
  {
    case PROP_BOUND:
      nc_galaxy_shape_factor_quad_set_bound (gsfq, g_value_get_double (value));
      break;
    case PROP_RELTOL:
      nc_galaxy_shape_factor_quad_set_reltol (gsfq, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_quad_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyShapeFactorQuad *gsfq = NC_GALAXY_SHAPE_FACTOR_QUAD (object);

  g_return_if_fail (NC_IS_GALAXY_SHAPE_FACTOR_QUAD (gsfq));

  switch (prop_id)
  {
    case PROP_BOUND:
      g_value_set_double (value, nc_galaxy_shape_factor_quad_get_bound (gsfq));
      break;
    case PROP_RELTOL:
      g_value_set_double (value, nc_galaxy_shape_factor_quad_get_reltol (gsfq));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_shape_factor_quad_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_quad_parent_class)->finalize (object);
}

static void
nc_galaxy_shape_factor_quad_class_init (NcGalaxyShapeFactorQuadClass *klass)
{
  NcGalaxyShapeFactorClass *gsf_class = NC_GALAXY_SHAPE_FACTOR_CLASS (klass);
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_nc_galaxy_shape_factor_quad_constructed;
  object_class->set_property = &_nc_galaxy_shape_factor_quad_set_property;
  object_class->get_property = &_nc_galaxy_shape_factor_quad_get_property;
  object_class->finalize     = &_nc_galaxy_shape_factor_quad_finalize;

  /**
   * NcGalaxyShapeFactorQuad:bound:
   *
   * Half-width $B$ of the box $[u_0-B,u_0+B]\times[v_0-B,v_0+B]$ (centered
   * on the exact preimage $(u_0,v_0)$ of $\epsilon_\mathrm{obs}$, see the
   * class documentation) over which the plane-substituted integral is
   * evaluated. Keep this generous: shrinking it does not help Divonne find
   * an already-hinted narrow feature, but does progressively truncate a
   * broad population's disc-spanning support; a few units already covers
   * the disc for any population.
   */
  g_object_class_install_property (object_class,
                                   PROP_BOUND,
                                   g_param_spec_double ("bound",
                                                        NULL,
                                                        "Plane-integration box half-width",
                                                        0.0, G_MAXDOUBLE, NC_GALAXY_SHAPE_FACTOR_QUAD_DEFAULT_BOUND,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcGalaxyShapeFactorQuad:reltol:
   *
   * Relative tolerance passed to the underlying Divonne cubature.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Cubature relative tolerance",
                                                        0.0, 1.0, NC_GALAXY_SHAPE_FACTOR_QUAD_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  gsf_class->data_init        = &_nc_galaxy_shape_factor_quad_data_init;
  gsf_class->prepare          = &_nc_galaxy_shape_factor_quad_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_quad_eval_marginal;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_quad_eval_ln_marginal;
}

/**
 * nc_galaxy_shape_factor_quad_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 *
 * Creates a new #NcGalaxyShapeFactorQuad.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorQuad.
 */
NcGalaxyShapeFactorQuad *
nc_galaxy_shape_factor_quad_new (NcGalaxyWLObsEllipConv ellip_conv)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_QUAD,
                       "ellip-conv", ellip_conv,
                       NULL);
}

/**
 * nc_galaxy_shape_factor_quad_ref:
 * @gsfq: a #NcGalaxyShapeFactorQuad
 *
 * Increases the reference count of @gsfq by one.
 *
 * Returns: (transfer full): @gsfq.
 */
NcGalaxyShapeFactorQuad *
nc_galaxy_shape_factor_quad_ref (NcGalaxyShapeFactorQuad *gsfq)
{
  return g_object_ref (gsfq);
}

/**
 * nc_galaxy_shape_factor_quad_free:
 * @gsfq: a #NcGalaxyShapeFactorQuad
 *
 * Decreases the reference count of @gsfq by one.
 *
 */
void
nc_galaxy_shape_factor_quad_free (NcGalaxyShapeFactorQuad *gsfq)
{
  g_object_unref (gsfq);
}

/**
 * nc_galaxy_shape_factor_quad_clear:
 * @gsfq: a #NcGalaxyShapeFactorQuad
 *
 * Decreases the reference count of *@gsfq by one, and sets the pointer
 * *@gsfq to NULL.
 *
 */
void
nc_galaxy_shape_factor_quad_clear (NcGalaxyShapeFactorQuad **gsfq)
{
  g_clear_object (gsfq);
}

/**
 * nc_galaxy_shape_factor_quad_set_bound:
 * @gsfq: a #NcGalaxyShapeFactorQuad
 * @bound: the new box half-width $B$
 *
 * Sets the half-width of the $[-B,B]^2$ plane-integration box.
 *
 */
void
nc_galaxy_shape_factor_quad_set_bound (NcGalaxyShapeFactorQuad *gsfq, const gdouble bound)
{
  NcGalaxyShapeFactorQuadPrivate * const self = nc_galaxy_shape_factor_quad_get_instance_private (gsfq);

  g_assert_cmpfloat (bound, >, 0.0);

  self->bound = bound;
}

/**
 * nc_galaxy_shape_factor_quad_get_bound:
 * @gsfq: a #NcGalaxyShapeFactorQuad
 *
 * Returns: the box half-width $B$.
 */
gdouble
nc_galaxy_shape_factor_quad_get_bound (NcGalaxyShapeFactorQuad *gsfq)
{
  NcGalaxyShapeFactorQuadPrivate * const self = nc_galaxy_shape_factor_quad_get_instance_private (gsfq);

  return self->bound;
}

/**
 * nc_galaxy_shape_factor_quad_set_reltol:
 * @gsfq: a #NcGalaxyShapeFactorQuad
 * @reltol: the new cubature relative tolerance
 *
 * Sets the relative tolerance passed to the underlying Divonne cubature.
 *
 */
void
nc_galaxy_shape_factor_quad_set_reltol (NcGalaxyShapeFactorQuad *gsfq, const gdouble reltol)
{
  NcGalaxyShapeFactorQuadPrivate * const self = nc_galaxy_shape_factor_quad_get_instance_private (gsfq);

  g_assert_cmpfloat (reltol, >, 0.0);

  self->reltol = reltol;
}

/**
 * nc_galaxy_shape_factor_quad_get_reltol:
 * @gsfq: a #NcGalaxyShapeFactorQuad
 *
 * Returns: the cubature relative tolerance.
 */
gdouble
nc_galaxy_shape_factor_quad_get_reltol (NcGalaxyShapeFactorQuad *gsfq)
{
  NcGalaxyShapeFactorQuadPrivate * const self = nc_galaxy_shape_factor_quad_get_instance_private (gsfq);

  return self->reltol;
}

