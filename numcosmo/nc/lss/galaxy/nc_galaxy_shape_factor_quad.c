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
 * **Substitution.** The integration variable is $\psi$, not $\chi_I$ or the
 * lensed $\chi_L=f_g(\chi_I)$ directly: $\chi_L=f_h(\psi)$, where $h$ is the
 * shear with $f_h(0)=\epsilon_\mathrm{obs}$ exactly
 * (nc_wl_ellipticity_shear_at_origin_trace() /
 * `_trace_det()`). $\psi=0$ therefore always maps to $\chi_L=\epsilon_\mathrm{obs}$,
 * so the noise kernel's peak sits at the box center by construction, with
 * no per-call box re-centering. $\chi_I=f_g^{-1}(\chi_L)$ as usual, with
 * $\det J_{f_g^{-1}}$ from nc_galaxy_shape_factor_lndet_jac()'s sibling
 * kernels (already used by `VarAdd`). $\psi$ is turned into a plane
 * integral by the usual bijection,
 * $$\psi(u,v) = (u,v)/\sqrt{1+u^2+v^2}, \qquad
 *   \left|\det\frac{\partial\psi}{\partial(u,v)}\right| = \frac{1}{(1+u^2+v^2)^2},$$
 * so $(u,v)$ ranges over the whole plane and $\chi_I$ correspondingly over
 * the whole disc, with no truncation. $x_I=\lvert\chi_I\rvert^2$ is
 * computed directly from the complex division (a sum of two squares, not
 * itself subtraction-sensitive), then $r_I=\sqrt{x_I}$ before the
 * nc_galaxy_shape_pop_eval_p() call (see that vfunc's r-native contract in
 * nc_galaxy_shape_pop.h), and the returned r-marginal density is converted
 * to the 2D area density this integrand needs via $P_\mathrm{pop}(r_I)/(2\pi r_I)$.
 *
 * **Puncture correction.** The integrand evaluates only the bracket
 * $N_2(\epsilon_\mathrm{obs}-\chi_L)-N_0$, $N_0=N_2(\epsilon_\mathrm{obs}-\chi_{L0})$,
 * $\chi_{L0}=f_g(0)$ (the population's own peak location under $g$,
 * unrelated to $h$): $N_0$ is added back once after the cubature, using
 * $\int P_\mathrm{pop}=1$ exactly (valid because $\chi_I$ ranges over the
 * whole disc, never a truncated subset). This makes the integrand vanish
 * smoothly at $\chi_I=0$ for $\alpha\geq1$ populations instead of the raw
 * integrand's sharp near-divergent peak there, so Divonne needs far fewer
 * evaluations. The bracket itself is $N_0\,\mathrm{expm1}(\Delta)$,
 * $\Delta=(d_0-d_i)/(2\sigma_\mathrm{noise}^2)$, whenever $\Delta$ stays
 * below a generous safety bound (expm1 avoids the cancellation direct
 * subtraction would suffer near $\Delta=0$); above that bound expm1 would
 * overflow exactly like $\exp$ would, with no cancellation left to protect
 * against either, so the bracket falls back to the plain difference of two
 * individually-bounded noise evaluations instead.
 *
 * **Cubature.** Divonne (ncm_integrate_2dim_divonne()) integrates $\psi$
 * directly in its own polar coordinates $(\rho,\theta)\in[0,1)\times[-\pi,\pi)$
 * -- an exact, finite domain covering the whole disc at $\rho=1$, with no
 * asymptotic tail to truncate -- seeded with two peak hints: $\psi\approx0$
 * (trivial, always the exact noise peak) and the $\psi$-preimage of
 * $\chi_{L0}=f_g(0)$ (the population's peak location under the common
 * assumption that $P_\mathrm{pop}$ is radially symmetric about $\chi_I=0$;
 * only approximate otherwise, e.g. #NcGalaxyShapePopBeta with $\mu$ away
 * from 0, but Divonne's stratified search still finds the true peak from
 * there). A third hint (the peak-ring point closest to $\epsilon_\mathrm{obs}$)
 * is added when the population's mode is not at $\chi_I=0$.
 *
 * The integrand clamps any non-finite evaluation to zero before returning:
 * Cuba can segfault outright on a NaN/Inf sample, and some populations have
 * a genuine, mathematically correct divergence at a disc point (e.g.
 * #NcGalaxyShapePopBeta with $\alpha<1$ diverges at $x=0$, which can
 * coincide with a peak hint).
 *
 * This scheme is exact but substantially more expensive per evaluation than
 * `VarAdd` (a full 2D cubature vs. one closed-form expression). It is meant
 * as an accuracy reference / fallback for regimes where the variance-add
 * approximation is not trusted, not as a routine replacement in
 * large-catalog likelihood evaluations. Cuba's own fork()-based internal
 * parallelism is disabled globally by ncm_cfg_init() (`cubacores(0, 0)`),
 * so this is safe to call concurrently from multiple OpenMP threads.
 *
 * See docs/theory/wl_shape_factor_history.md for the evidence behind this
 * design and earlier approaches tried and rejected.
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

#define NC_GALAXY_SHAPE_FACTOR_QUAD_DEFAULT_RELTOL 1.0e-7

/* exp(50) ~ 5e21, nowhere near double overflow (~exp(709)) -- generous
 * margin, not a tight cutoff (see the integrand's own doc). */
#define NC_GALAXY_SHAPE_FACTOR_QUAD_EXPM1_SAFE_BOUND 50.0

/* rho=0 is a coordinate singularity of the polar (rho,theta) domain
 * (jac_polar=rho vanishes there); nudged off it so the psi~0 hint lands
 * where the integrand is actually large instead of exactly zero. */
#define NC_GALAXY_SHAPE_FACTOR_QUAD_RHO_HINT_EPS 1.0e-6

struct _NcGalaxyShapeFactorQuad
{
  NcGalaxyShapeFactor parent_instance;
};

typedef struct _NcGalaxyShapeFactorQuadPrivate
{
  /* Convention specialization resolved once at construction (ellip-conv is
   * CONSTRUCT_ONLY), keeping the switch out of the per-evaluation path. */
  complex double (*apply_shear) (complex double g, complex double chi_i);
  complex double (*apply_shear_inv) (complex double g, complex double chi_L);
  complex double (*shear_at_origin) (complex double target);

  gdouble (*det_jac) (complex double g, complex double chi_L);

  gdouble reltol;
} NcGalaxyShapeFactorQuadPrivate;

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorQuad, nc_galaxy_shape_factor_quad, NC_TYPE_GALAXY_SHAPE_FACTOR)

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_LEN,
};

/* Divonne integrand ------------------------------------------------------ */

typedef struct _QuadIntArg
{
  NcGalaxyShapePop *pop;
  NcGalaxyShapePopData *pop_data;
  complex double (*apply_shear) (complex double g, complex double chi_i);
  complex double (*apply_shear_inv) (complex double g, complex double chi_L);

  gdouble (*det_jac) (complex double g, complex double chi_L);

  complex double g;
  complex double h;
  complex double eps_obs;
  gdouble noise_var;
  gdouble d0;
  gdouble N0;
} QuadIntArg;

/* See this file's class doc for the psi-substitution and puncture
 * correction this implements. */
static gdouble
_nc_galaxy_shape_factor_quad_integrand (gdouble rho, gdouble theta, gpointer userdata)
{
  QuadIntArg * const arg     = (QuadIntArg *) userdata;
  const gdouble jac_polar    = rho;
  const complex double psi   = rho * cexp (I * theta);
  const complex double chi_L = arg->apply_shear (arg->h, psi);
  const complex double chi_i = arg->apply_shear_inv (arg->g, chi_L);
  const gdouble x_i          = gsl_pow_2 (creal (chi_i)) + gsl_pow_2 (cimag (chi_i));
  const gdouble r_i          = sqrt (x_i);
  const gdouble p_pop        = nc_galaxy_shape_pop_eval_p (arg->pop, arg->pop_data, r_i) / (2.0 * M_PI * r_i);
  const gdouble jac_inv      = arg->det_jac (arg->g, chi_L);
  const gdouble jac_psi      = arg->det_jac (-arg->h, psi);
  const complex double delta = arg->eps_obs - chi_L;
  const gdouble d_i          = gsl_pow_2 (creal (delta)) + gsl_pow_2 (cimag (delta));
  const gdouble delta_ratio  = (arg->d0 - d_i) / (2.0 * arg->noise_var);
  gdouble bracket, ret;

  if (delta_ratio <= NC_GALAXY_SHAPE_FACTOR_QUAD_EXPM1_SAFE_BOUND)
  {
    bracket = arg->N0 * expm1 (delta_ratio);
  }
  else
  {
    const gdouble noise_i = exp (-d_i / (2.0 * arg->noise_var)) / (2.0 * M_PI * arg->noise_var);

    bracket = noise_i - arg->N0;
  }

  ret = p_pop * jac_inv * jac_psi * jac_polar * bracket;

  /* Cuba can segfault on a non-finite sample, so clamp rather than
   * propagate (a genuine population divergence, e.g. Beta with alpha<1 at
   * r_i=0, can coincide with a peak hint). */
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

/* Refines theta so chi_I=rho_mode*exp(I*theta), on the ring where a
 * non-radially-symmetric population peaks, maps as close as possible to
 * eps_obs under f_g. The population density is flat along this ring (it
 * depends only on x=|chi_I|^2), so this is a pure geometric optimization:
 * 3 Newton steps on d/dtheta |eps_obs-f_g(rho_mode e^{i theta})|^2 with
 * central-difference derivatives, seeded at the naive noiseless preimage's
 * phase. */
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
_nc_galaxy_shape_factor_quad_eval (NcGalaxyShapeFactorQuad *gsfq, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, const gdouble std_noise, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorQuadPrivate * const self = nc_galaxy_shape_factor_quad_get_instance_private (gsfq);
  const complex double g                      = g_1 + I * g_2;
  const complex double eps_obs                = epsilon_obs_1 + I * epsilon_obs_2;

  /* shear_at_origin() requires |eps_obs|<1 strictly; unlike FixedQuad, this
   * class has no fallback for an out-of-disc eps_obs (a real possibility
   * once noise is added), so fail loudly instead of returning a wrong
   * answer. */
  if (cabs (eps_obs) >= 1.0)
    g_error ("nc_galaxy_shape_factor_quad: |eps_obs|=% .6g >= 1 at eps_obs=(% .6g,% .6g) -- "
             "outside this class' valid domain; use NcGalaxyShapeFactorFixedQuad instead.",
             cabs (eps_obs), epsilon_obs_1, epsilon_obs_2);

  const complex double h      = self->shear_at_origin (eps_obs);
  const complex double chi_L0 = self->apply_shear (g, 0.0);
  const complex double delta0 = eps_obs - chi_L0;
  const gdouble noise_var     = gsl_pow_2 (std_noise);
  const gdouble d0            = gsl_pow_2 (creal (delta0)) + gsl_pow_2 (cimag (delta0));
  const gdouble N0            = exp (-d0 / (2.0 * noise_var)) / (2.0 * M_PI * noise_var);
  const gdouble mode_r        = nc_galaxy_shape_pop_get_mode_r (pop, pop_data);
  QuadIntArg arg;
  NcmIntegrand2dim integ = {&arg, &_nc_galaxy_shape_factor_quad_integrand};
  gdouble result, error;
  gdouble xgiven[6];
  gint ngiven;

  /* Two peak hints, as (rho,theta) pairs: psi~0 (the noise peak, nudged off
   * the rho=0 coordinate singularity) and the population's peak chi_L0,
   * converted to psi-space via apply_shear(-h,.). */
  xgiven[0] = NC_GALAXY_SHAPE_FACTOR_QUAD_RHO_HINT_EPS;
  xgiven[1] = 0.0;
  {
    const complex double psi_pop_peak = self->apply_shear (-h, chi_L0);

    xgiven[2] = cabs (psi_pop_peak);
    xgiven[3] = carg (psi_pop_peak);
  }
  ngiven = 2;

  /* Third hint, only when the population's peak ring is not at chi_I=0: the
   * point on |chi_I|=rho_mode closest, after the forward map, to eps_obs
   * (_nc_galaxy_shape_factor_quad_refine_theta()), converted to psi-space
   * the same way. */
  if (mode_r > 0.0)
  {
    const gdouble rho_mode           = mode_r;
    const complex double chi_i_naive = self->apply_shear_inv (g, eps_obs);
    const gdouble theta0             = carg (chi_i_naive);
    const gdouble theta              = _nc_galaxy_shape_factor_quad_refine_theta (self->apply_shear, g, rho_mode, eps_obs, theta0);
    const complex double chi_i_hint  = rho_mode * cexp (I * theta);
    const complex double chi_L_hint  = self->apply_shear (g, chi_i_hint);
    const complex double psi_hint    = self->apply_shear (-h, chi_L_hint);

    xgiven[4] = cabs (psi_hint);
    xgiven[5] = carg (psi_hint);
    ngiven    = 3;
  }

  arg.pop             = pop;
  arg.pop_data        = pop_data;
  arg.apply_shear     = self->apply_shear;
  arg.apply_shear_inv = self->apply_shear_inv;
  arg.det_jac         = self->det_jac;
  arg.g               = g;
  arg.h               = h;
  arg.eps_obs         = eps_obs;
  arg.noise_var       = noise_var;
  arg.d0              = d0;
  arg.N0              = N0;

  ncm_integrate_2dim_divonne (&integ, 0.0, -M_PI, 1.0, M_PI,
                              self->reltol, 0.0, ngiven, 2, xgiven, &result, &error);

  return N0 + result;
}

static gdouble
_nc_galaxy_shape_factor_quad_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return _nc_galaxy_shape_factor_quad_eval (NC_GALAXY_SHAPE_FACTOR_QUAD (gsf), pop, data->pop_data, data->std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2);
}

static gdouble
_nc_galaxy_shape_factor_quad_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return log (_nc_galaxy_shape_factor_quad_eval (NC_GALAXY_SHAPE_FACTOR_QUAD (gsf), pop, data->pop_data, data->std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2));
}

/* GObject boilerplate --------------------------------------------------- */

static void
nc_galaxy_shape_factor_quad_init (NcGalaxyShapeFactorQuad *gsfq)
{
  NcGalaxyShapeFactorQuadPrivate * const self = nc_galaxy_shape_factor_quad_get_instance_private (gsfq);

  self->apply_shear     = NULL;
  self->apply_shear_inv = NULL;
  self->shear_at_origin = NULL;
  self->det_jac         = NULL;
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
        self->shear_at_origin = &nc_wl_ellipticity_shear_at_origin_trace;
        self->det_jac         = &nc_wl_ellipticity_det_jac_trace;
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        self->apply_shear     = &nc_wl_ellipticity_apply_shear_trace_det;
        self->apply_shear_inv = &nc_wl_ellipticity_apply_shear_inv_trace_det;
        self->shear_at_origin = &nc_wl_ellipticity_shear_at_origin_trace_det;
        self->det_jac         = &nc_wl_ellipticity_det_jac_trace_det;
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

/**
 * nc_galaxy_shape_factor_quad_eval_direct:
 * @gsfq: a #NcGalaxyShapeFactorQuad
 * @pop: a #NcGalaxyShapePop
 * @g_1: reduced shear, real component
 * @g_2: reduced shear, imaginary component
 * @epsilon_obs_1: observed ellipticity/distortion, real component
 * @epsilon_obs_2: observed ellipticity/distortion, imaginary component
 * @std_noise: measurement noise standard deviation
 *
 * Evaluates the marginal directly from the raw physical parameters,
 * building and preparing a fresh #NcGalaxyShapePopData internally so a
 * caller never needs the #NcGalaxyShapeFactorData/#NcmMSet machinery (and
 * cannot forget to prepare it -- an unprepared #NcGalaxyShapePopData
 * silently reads back a wrong placeholder density, see
 * nc_galaxy_shape_pop_prepare()). Intended for tests and standalone
 * evaluation; nc_galaxy_shape_factor_eval_marginal() is the entry point for
 * catalog-scale use, which reuses one #NcGalaxyShapePopData per galaxy
 * across many @g_1/@g_2 evaluations instead of rebuilding it every call.
 *
 * Returns: $P(\epsilon_\mathrm{obs} \mid g)$.
 */
gdouble
nc_galaxy_shape_factor_quad_eval_direct (NcGalaxyShapeFactorQuad *gsfq, NcGalaxyShapePop *pop, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble std_noise)
{
  NcGalaxyShapePopData *pop_data = nc_galaxy_shape_pop_data_new (pop);
  gdouble result;

  nc_galaxy_shape_pop_prepare (pop, pop_data);
  result = _nc_galaxy_shape_factor_quad_eval (gsfq, pop, pop_data, std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2);
  nc_galaxy_shape_pop_data_unref (pop_data);

  return result;
}

