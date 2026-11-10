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
 * At $|g|=1$ every point of the open disc maps to $g$, and the class returns the
 * closed form $P=N_2(\epsilon_\mathrm{obs}-g;\sigma_\mathrm{noise}^2)$.
 *
 * Otherwise the default path integrates by Divonne over the unit disc in
 * $\psi=\rho\,e^{i\theta}$, $\rho\in[0,1]$, with $\chi_L=f_h(\psi)$ and $h$ chosen
 * so that $f_h(0)$ is a center point $c$. For $|\epsilon_\mathrm{obs}|<1$,
 * $c=\epsilon_\mathrm{obs}$. For $|\epsilon_\mathrm{obs}|\ge1$ no point of the disc
 * maps to $\epsilon_\mathrm{obs}$, and $c$ is the point in the direction of
 * $\epsilon_\mathrm{obs}$ at depth
 * $\sigma_\mathrm{noise}^2/(|\epsilon_\mathrm{obs}|-1+\sigma_\mathrm{noise})$ inside
 * the circle (clamped to $[10^{-8},0.5]$), where the noise mass lies. The
 * measure is the population density at $\chi_I=f_g^{-1}(\chi_L)$ times the
 * Jacobians of $\psi\to\chi_L\to\chi_I$ and the polar factor $\rho$.
 * nc_galaxy_shape_pop_eval_p() is called on $r_I=|\chi_I|$, following the
 * r-native contract in nc_galaxy_shape_pop.h, and the r-marginal density it
 * returns is converted to the 2D area density as $P_\mathrm{pop}(r_I)/(2\pi r_I)$.
 *
 * The integrand is divided by its peak value, located on a 24 x 48 grid in
 * $(\rho,\theta)$ and at the peak hints, so that neither $P$ nor any factor of
 * the integrand leaves the double range; the class computes $\ln P$ and
 * nc_galaxy_shape_factor_eval_marginal() returns its exponential. The peak
 * hints passed to Divonne are $\psi\approx0$, the preimage of $f_g(0)$, the
 * point of the population's mode ring closest to $\epsilon_\mathrm{obs}$ after
 * the map (when the mode is not at $\chi_I=0$), and, for
 * $|\epsilon_\mathrm{obs}|\ge1$, the circle point that maps to
 * $\epsilon_\mathrm{obs}/|\epsilon_\mathrm{obs}|$.
 *
 * Divonne integrates the scaled integrand itself. Subtracting the noise kernel
 * at $f_g(0)$ from the scaled integrand made Divonne return wrong results, not
 * return, or crash at ordinary points. Divonne's convergence flag is not
 * checked: it is false at points where the result is accurate, and the
 * class's accuracy is established against the Arb reference instead.
 *
 * A second path is selected by the temporary property
 * #NcGalaxyShapeFactorQuad:use-lp, for the comparison against the
 * arbitrary-precision reference tests/tools/nc_galaxy_shape_marginal_arb.c: it
 * integrates in log-polar coordinates around the noise peak. One of the two
 * paths will be kept and the property removed.
 *
 * Non-finite integrand values are returned as zero. Cuba can segfault on a
 * NaN or infinite sample, and some populations diverge at a disc point as a
 * matter of mathematics rather than error: #NcGalaxyShapePopBeta with
 * $\alpha<1$ diverges at $x=0$, which can coincide with a peak hint.
 *
 * Uses 2D cubature, so it costs substantially more per evaluation than
 * #NcGalaxyShapeFactorVarAdd's closed-form expression. It is intended as an
 * accuracy reference and fallback where the variance-add approximation is not
 * trusted, not as a routine replacement in large-catalog likelihood
 * evaluations.
 *
 * Cuba's fork()-based internal parallelism is disabled globally by
 * ncm_cfg_init() (`cubacores(0, 0)`), so this is safe to call concurrently
 * from multiple OpenMP threads.
 *
 * See docs/theory/nc/lss/galaxy/wl_shape_factor_history.md for the design history.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_quad.h"
#include "nc/lss/galaxy/nc_galaxy_shape_pop_beta.h"
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

/* Coarse grid in psi (rho x theta) on which the integrand's peak is located. */
#define NC_GALAXY_SHAPE_FACTOR_QUAD_PEAK_NR 24
#define NC_GALAXY_SHAPE_FACTOR_QUAD_PEAK_NT 48

/* Lower bound on r_I in the log-polar integrand: eval_p (r) / (2 pi r) is 0/0 at r = 0. */
#define NC_GALAXY_SHAPE_FACTOR_QUAD_LP_R_MIN 1.0e-150

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
  gboolean use_lp;
} NcGalaxyShapeFactorQuadPrivate;

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorQuad, nc_galaxy_shape_factor_quad, NC_TYPE_GALAXY_SHAPE_FACTOR)

enum
{
  PROP_0,
  PROP_RELTOL,
  PROP_USE_LP,
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
  gdouble s_p;
  gdouble s_N;
} QuadIntArg;

/*
 * At psi = rho e^{i theta}: writes to @pj the population's 2D density at chi_I times the
 * Jacobians of psi -> chi_L -> chi_I and of the polar coordinates, and to @d the squared
 * distance |eps_obs - chi_L|^2.
 */
static void
_nc_galaxy_shape_factor_quad_point (const QuadIntArg *arg, const gdouble rho, const gdouble theta, gdouble *pj, gdouble *d)
{
  const complex double psi   = rho * cexp (I * theta);
  const complex double chi_L = arg->apply_shear (arg->h, psi);
  const complex double chi_i = arg->apply_shear_inv (arg->g, chi_L);
  const gdouble r_i          = sqrt (gsl_pow_2 (creal (chi_i)) + gsl_pow_2 (cimag (chi_i)));
  const gdouble p_pop        = nc_galaxy_shape_pop_eval_p (arg->pop, arg->pop_data, r_i) / (2.0 * M_PI * r_i);
  const complex double delta = arg->eps_obs - chi_L;

  *pj = p_pop * arg->det_jac (arg->g, chi_L) * arg->det_jac (-arg->h, psi) * rho;
  *d  = gsl_pow_2 (creal (delta)) + gsl_pow_2 (cimag (delta));
}

/*
 * The integrand scaled to be of order one at its peak: pj e^{-s_p} e^{-d/2s^2 + s_N}, where s_p
 * and s_N are the logs of pj and of the noise factor at the peak found before the integral.
 * See this file's class doc for the psi-substitution.
 */
static gdouble
_nc_galaxy_shape_factor_quad_integrand (gdouble rho, gdouble theta, gpointer userdata)
{
  QuadIntArg * const arg = (QuadIntArg *) userdata;
  gdouble pj, d_i, ret;

  _nc_galaxy_shape_factor_quad_point (arg, rho, theta, &pj, &d_i);

  ret = pj * exp (-arg->s_p - d_i / (2.0 * arg->noise_var) + arg->s_N);

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

/* Log-polar coordinates around the recentered noise peak ---------------- */

typedef struct _QuadLPArg
{
  NcGalaxyShapePop *pop;
  NcGalaxyShapePopData *pop_data;
  complex double (*apply_shear_inv) (complex double g, complex double chi_L);

  gdouble (*det_jac) (complex double g, complex double chi_L);

  complex double g;
  complex double c;
  complex double eps_obs;
  gdouble noise_var;
  gdouble s_p;
  gdouble s_N;
  gdouble xi_p;
  gdouble th_p;
  gdouble ea;
  gdouble eb;
  gdouble ed;
  gint tri;
} QuadLPArg;

/*
 * At psi = e^{xi + i theta}, chi_L = (psi + c) / (1 + c^* psi): writes to @pj the
 * population's 2D density at chi_I = f_g^{-1} (chi_L) times the measure
 * e^{2 xi} |dchi_L/dpsi|^2 |det dchi_I/dchi_L|, and to @d the squared distance
 * |eps_obs - chi_L|^2.
 */
static void
_nc_galaxy_shape_factor_quad_lp_point (const QuadLPArg *arg, const gdouble xi, const gdouble th, gdouble *pj, gdouble *d)
{
  const complex double psi   = exp (xi) * cexp (I * th);
  const complex double den   = 1.0 + conj (arg->c) * psi;
  const complex double chi_L = (psi + arg->c) / den;
  const complex double chi_i = arg->apply_shear_inv (arg->g, chi_L);
  const gdouble r_i          = GSL_MAX (cabs (chi_i), NC_GALAXY_SHAPE_FACTOR_QUAD_LP_R_MIN);
  const gdouble p_pop        = nc_galaxy_shape_pop_eval_p (arg->pop, arg->pop_data, r_i) / (2.0 * M_PI * r_i);
  const gdouble mjac         = gsl_pow_2 (1.0 - gsl_pow_2 (cabs (arg->c))) / gsl_pow_4 (cabs (den));
  const complex double delta = arg->eps_obs - chi_L;

  *pj = p_pop * arg->det_jac (arg->g, chi_L) * mjac * exp (2.0 * xi);
  *d  = gsl_pow_2 (creal (delta)) + gsl_pow_2 (cimag (delta));
}

static gdouble
_nc_galaxy_shape_factor_quad_lp_integrand (gdouble xi, gdouble th, gpointer userdata)
{
  QuadLPArg * const arg = (QuadLPArg *) userdata;
  gdouble pj, d, ret;

  _nc_galaxy_shape_factor_quad_lp_point (arg, xi, th, &pj, &d);
  ret = pj * exp (-arg->s_p - d / (2.0 * arg->noise_var) + arg->s_N);

  return isfinite (ret) ? ret : 0.0;
}

/*
 * Triangle arg->tri of the rectangle [xi_p - ea, xi_p + eb] x [th_p - ed, th_p + ed] around
 * the population's center (xi_p, th_p), with apex at the center and base on one edge:
 * (xi, th) = (xi_p, th_p) + t (P0 + s P1), s in [-1, 1], Jacobian t |P0 x P1|, which cancels
 * the |chi_I|^(alpha - 2) singularity at t = 0. Edges: xi = xi_p + eb, xi = xi_p - ea,
 * th = th_p + ed, th = th_p - ed.
 */
static gdouble
_nc_galaxy_shape_factor_quad_lp_tri_integrand (gdouble t, gdouble s, gpointer userdata)
{
  QuadLPArg * const arg = (QuadLPArg *) userdata;
  const gdouble mid     = 0.5 * (arg->eb - arg->ea);
  const gdouble half    = 0.5 * (arg->ea + arg->eb);
  gdouble dx, dy, jac;

  switch (arg->tri)
  {
    case 0:
      dx  = arg->eb;
      dy  = s * arg->ed;
      jac = arg->eb * arg->ed;
      break;
    case 1:
      dx  = -arg->ea;
      dy  = s * arg->ed;
      jac = arg->ea * arg->ed;
      break;
    case 2:
      dx  = mid + s * half;
      dy  = arg->ed;
      jac = half * arg->ed;
      break;
    default:
      dx  = mid + s * half;
      dy  = -arg->ed;
      jac = half * arg->ed;
      break;
  }

  return _nc_galaxy_shape_factor_quad_lp_integrand (arg->xi_p + t * dx, arg->th_p + t * dy, userdata) * jac * t;
}

/* Adds the Divonne integral of @f over [x0, x1] x [y0, y1] to @sum, with the hints inside it. */
static void
_nc_galaxy_shape_factor_quad_lp_box (gdouble ( *f ) (gdouble, gdouble, gpointer), QuadLPArg *arg, const gdouble reltol,
                                     const gdouble x0, const gdouble y0, const gdouble x1, const gdouble y1,
                                     const gdouble *hints, const gint nhints, gdouble *sum)
{
  NcmIntegrand2dim integ = {arg, f};
  gdouble xg[16];
  gint ng = 0, i;
  gdouble result, error;

  if ((x1 <= x0) || (y1 <= y0))
    return;

  for (i = 0; i < nhints; i++)
  {
    if ((hints[2 * i] > x0) && (hints[2 * i] < x1) && (hints[2 * i + 1] > y0) && (hints[2 * i + 1] < y1))
    {
      xg[2 * ng]     = hints[2 * i];
      xg[2 * ng + 1] = hints[2 * i + 1];
      ng++;
    }
  }

  ncm_integrate_2dim_divonne (&integ, x0, y0, x1, y1, reltol, 0.0, ng, 2, xg, &result, &error);
  *sum += result;
}

/*
 * ln P in log-polar coordinates (xi, theta) around the noise peak, recentered by the disc
 * automorphism chi_L = (psi + c) / (1 + c^* psi), c = eps_obs or, for |eps_obs| >= 1, the
 * point at the depth of the noise layer inside the circle. The noise band sits at
 * xi ~ ln sigma_psi, sigma_psi = std_noise / (1 - |c|^2), with width of order one for any
 * std_noise. For a Beta population the point where chi_I = 0 is cut out by a square of
 * half-side delta and integrated on four triangles; when that point is closer to psi = 0
 * than e^{xi_min}, the log-polar measure already cancels the power law.
 */
static gdouble
_nc_galaxy_shape_factor_quad_eval_lp (NcGalaxyShapeFactorQuadPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, const gdouble std_noise, const complex double g, const complex double eps_obs)
{
  const gdouble eps_abs       = cabs (eps_obs);
  const gboolean outside      = (eps_abs >= 1.0);
  const gdouble layer         = outside ? GSL_MIN (GSL_MAX (gsl_pow_2 (std_noise) / (eps_abs - 1.0 + std_noise), 1.0e-8), 0.5) : 0.0;
  const complex double c      = outside ? (1.0 - layer) * eps_obs / eps_abs : eps_obs;
  const gdouble noise_var     = gsl_pow_2 (std_noise);
  const gdouble sig_psi       = GSL_MIN (std_noise / (1.0 - gsl_pow_2 (cabs (c))), 1.0);
  const gdouble xi_min        = log (sig_psi) - 15.0;
  const complex double chi_L0 = self->apply_shear (g, 0.0);
  const complex double psi_p  = (chi_L0 - c) / (1.0 - conj (c) * chi_L0);
  const gdouble th_p          = carg (psi_p);
  const gdouble th0           = th_p - M_PI;
  const gdouble th1           = th_p + M_PI;
  const gdouble mode_r        = nc_galaxy_shape_pop_get_mode_r (pop, pop_data);
  gdouble hints[8];
  gint nh = 0;
  QuadLPArg arg;
  gdouble J = 0.0;

  arg.pop             = pop;
  arg.pop_data        = pop_data;
  arg.apply_shear_inv = self->apply_shear_inv;
  arg.det_jac         = self->det_jac;
  arg.g               = g;
  arg.c               = c;
  arg.eps_obs         = eps_obs;
  arg.noise_var       = noise_var;
  arg.s_p             = 0.0;
  arg.s_N             = 0.0;
  arg.xi_p            = log (GSL_MAX (cabs (psi_p), 1.0e-300));
  arg.th_p            = th_p;
  arg.ea              = 0.0;
  arg.eb              = 0.0;
  arg.ed              = 0.0;
  arg.tri             = 0;

  /* Hints: the noise band, the population's center and ring, the boundary point. */
  hints[2 * nh]     = GSL_MIN (log (sig_psi), -1.0e-3);
  hints[2 * nh + 1] = th_p + 0.5 * M_PI;
  nh++;

  if ((arg.xi_p > xi_min) && (arg.xi_p < 0.0))
  {
    hints[2 * nh]     = arg.xi_p;
    hints[2 * nh + 1] = th_p;
    nh++;
  }

  if (mode_r > 0.0)
  {
    const complex double chi_i_naive = self->apply_shear_inv (g, eps_obs);
    const gdouble theta              = _nc_galaxy_shape_factor_quad_refine_theta (self->apply_shear, g, mode_r, eps_obs, carg (chi_i_naive));
    const complex double chi_L_hint  = self->apply_shear (g, mode_r * cexp (I * theta));
    const complex double psi_h       = (chi_L_hint - c) / (1.0 - conj (c) * chi_L_hint);
    gdouble a                        = carg (psi_h);

    while (a < th0)
      a += 2.0 * M_PI;

    while (a > th1)
      a -= 2.0 * M_PI;

    hints[2 * nh]     = log (GSL_MAX (cabs (psi_h), 1.0e-300));
    hints[2 * nh + 1] = a;
    nh++;
  }

  if (outside)
  {
    const complex double b   = eps_obs / eps_abs;
    const complex double psb = (b - c) / (1.0 - conj (c) * b);
    gdouble a                = carg (psb);

    while (a < th0)
      a += 2.0 * M_PI;

    while (a > th1)
      a -= 2.0 * M_PI;

    hints[2 * nh]     = -1.0e-6;
    hints[2 * nh + 1] = a;
    nh++;
  }

  /* Peak of the integrand on a coarse grid and the hints: the shifts s_p and s_N. */
  {
    gdouble best = -GSL_POSINF, pj, d;
    gint i;

    for (i = -nh; i < NC_GALAXY_SHAPE_FACTOR_QUAD_PEAK_NR * NC_GALAXY_SHAPE_FACTOR_QUAD_PEAK_NT; i++)
    {
      const gdouble xi = (i < 0) ? hints[2 * (i + nh)] : xi_min + (0.0 - xi_min) * ((i / NC_GALAXY_SHAPE_FACTOR_QUAD_PEAK_NT) + 0.5) / NC_GALAXY_SHAPE_FACTOR_QUAD_PEAK_NR;
      const gdouble th = (i < 0) ? hints[2 * (i + nh) + 1] : th0 + 2.0 * M_PI * ((i % NC_GALAXY_SHAPE_FACTOR_QUAD_PEAK_NT) + 0.5) / NC_GALAXY_SHAPE_FACTOR_QUAD_PEAK_NT;
      gdouble score;

      _nc_galaxy_shape_factor_quad_lp_point (&arg, xi, th, &pj, &d);
      score = log (pj) - d / (2.0 * noise_var);

      if (isfinite (score) && (score > best))
      {
        best    = score;
        arg.s_p = log (pj);
        arg.s_N = d / (2.0 * noise_var);
      }
    }
  }

  /*
   * The population's center gets its own square of half-side delta in (xi, theta): three
   * times the population's width e_rms mapped to log-polar coordinates by the local scales
   * of chi_I -> chi_L (|det dchi_I/dchi_L|^{-1/2}) and chi_L -> psi, over |psi_p|; at most
   * 0.25 and inside the box. Beta integrates it on four triangles (the power law at the
   * center), Gauss as one box. When the center is inside e^{xi_min}, the log-polar measure
   * already cancels the power law and nothing is cut out.
   */
  if (arg.xi_p > xi_min)
  {
    /*
     * Largest stretch of chi_I -> psi at chi_I = 0, the largest singular value of its
     * Jacobian by central differences: the TRACE map is not conformal, so the population's
     * image is an ellipse, and the square must cover its long axis.
     */
    const gdouble e_rms = nc_galaxy_shape_pop_e_rms (pop, pop_data);
    const gdouble hd    = 1.0e-7;
    complex double dpx, dpy;
    gdouble s1, w_log, dl;

    {
      complex double fp, fm;

      fp  = self->apply_shear (g, hd);
      fm  = self->apply_shear (g, -hd);
      dpx = ((fp - c) / (1.0 - conj (c) * fp) - (fm - c) / (1.0 - conj (c) * fm)) / (2.0 * hd);
      fp  = self->apply_shear (g, I * hd);
      fm  = self->apply_shear (g, -I * hd);
      dpy = ((fp - c) / (1.0 - conj (c) * fp) - (fm - c) / (1.0 - conj (c) * fm)) / (2.0 * hd);
    }
    {
      /* singular values of [[Re dpx, Re dpy], [Im dpx, Im dpy]] */
      const gdouble a11 = gsl_pow_2 (cabs (dpx));
      const gdouble a22 = gsl_pow_2 (cabs (dpy));
      const gdouble a12 = creal (dpx) * creal (dpy) + cimag (dpx) * cimag (dpy);
      const gdouble tr  = a11 + a22;
      const gdouble dt  = a11 * a22 - a12 * a12;

      s1 = sqrt (0.5 * (tr + sqrt (GSL_MAX (tr * tr - 4.0 * dt, 0.0))));
    }

    w_log = 3.0 * e_rms * s1 / exp (arg.xi_p);
    dl    = GSL_MIN (0.25, w_log);

    if (isfinite (dl) && (dl > 1.0e-12))
    {
      /* The rectangle may reach the circle xi = 0 in front, not the lower edge xi_min. */
      arg.ea = GSL_MIN (dl, 0.5 * (arg.xi_p - xi_min));
      arg.eb = GSL_MIN (dl, -arg.xi_p);
      arg.ed = GSL_MIN (dl, 0.5 * M_PI);
    }
  }

  if ((arg.ea > 0.0) && (arg.ed > 0.0))
  {
    const gdouble xp = arg.xi_p, ea = arg.ea, eb = arg.eb, ed = arg.ed;

    _nc_galaxy_shape_factor_quad_lp_box (&_nc_galaxy_shape_factor_quad_lp_integrand, &arg, self->reltol, xi_min, th0, xp - ea, th1, hints, nh, &J);
    _nc_galaxy_shape_factor_quad_lp_box (&_nc_galaxy_shape_factor_quad_lp_integrand, &arg, self->reltol, xp + eb, th0, 0.0, th1, hints, nh, &J);
    _nc_galaxy_shape_factor_quad_lp_box (&_nc_galaxy_shape_factor_quad_lp_integrand, &arg, self->reltol, xp - ea, th0, xp + eb, th_p - ed, hints, nh, &J);
    _nc_galaxy_shape_factor_quad_lp_box (&_nc_galaxy_shape_factor_quad_lp_integrand, &arg, self->reltol, xp - ea, th_p + ed, xp + eb, th1, hints, nh, &J);

    if (NC_IS_GALAXY_SHAPE_POP_BETA (pop))
    {
      for (arg.tri = 0; arg.tri < 4; arg.tri++)
        _nc_galaxy_shape_factor_quad_lp_box (&_nc_galaxy_shape_factor_quad_lp_tri_integrand, &arg, self->reltol, 0.0, -1.0, 1.0, 1.0, hints, 0, &J);
    }
    else
    {
      const gdouble center[2] = {xp, th_p};

      _nc_galaxy_shape_factor_quad_lp_box (&_nc_galaxy_shape_factor_quad_lp_integrand, &arg, self->reltol, xp - ea, th_p - ed, xp + eb, th_p + ed, center, 1, &J);
    }
  }
  else
  {
    _nc_galaxy_shape_factor_quad_lp_box (&_nc_galaxy_shape_factor_quad_lp_integrand, &arg, self->reltol, xi_min, th0, 0.0, th1, hints, nh, &J);
  }

  /* P = e^{s_p - s_N} J / (2 pi s^2) */
  return arg.s_p - arg.s_N + log (J) - log (2.0 * M_PI * noise_var);
}

/*
 * ln P. At |g| = 1 every point of the open disc maps to g, so P = N_2 (eps_obs - g) exactly.
 * Otherwise the integral is computed relative to the peak of its integrand, located by a
 * coarse search over psi and the hints, so that neither P nor any factor of the integrand
 * leaves the double range.
 */
static gdouble
_nc_galaxy_shape_factor_quad_eval_ln (NcGalaxyShapeFactorQuad *gsfq, NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data, const gdouble std_noise, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorQuadPrivate * const self = nc_galaxy_shape_factor_quad_get_instance_private (gsfq);
  const complex double g                      = g_1 + I * g_2;
  const complex double eps_obs                = epsilon_obs_1 + I * epsilon_obs_2;

  if (g_1 * g_1 + g_2 * g_2 == 1.0)
  {
    const complex double delta = eps_obs - g;
    const gdouble noise_var    = gsl_pow_2 (std_noise);

    return -(gsl_pow_2 (creal (delta)) + gsl_pow_2 (cimag (delta))) / (2.0 * noise_var) - log (2.0 * M_PI * noise_var);
  }

  if (self->use_lp)
    return _nc_galaxy_shape_factor_quad_eval_lp (self, pop, pop_data, std_noise, g, eps_obs);

  /*
   * The substitution is centered at psi = 0, the point f_h (0) = center. For |eps_obs| < 1
   * the center is eps_obs. Otherwise no point of the disc maps to eps_obs, and the noise
   * mass lies in a layer inside the circle next to eps_obs / |eps_obs|, of depth about
   * std_noise^2 / (|eps_obs| - 1 + std_noise); the center is placed at that depth.
   */
  const gdouble eps_abs       = cabs (eps_obs);
  const gboolean outside      = (eps_abs >= 1.0);
  const gdouble layer         = outside ? GSL_MIN (GSL_MAX (gsl_pow_2 (std_noise) / (eps_abs - 1.0 + std_noise), 1.0e-8), 0.5) : 0.0;
  const complex double center = outside ? (1.0 - layer) * eps_obs / eps_abs : eps_obs;
  const complex double h      = self->shear_at_origin (center);
  const complex double chi_L0 = self->apply_shear (g, 0.0);
  const gdouble noise_var     = gsl_pow_2 (std_noise);
  const gdouble mode_r        = nc_galaxy_shape_pop_get_mode_r (pop, pop_data);
  QuadIntArg arg;
  NcmIntegrand2dim integ = {&arg, &_nc_galaxy_shape_factor_quad_integrand};
  gdouble result, error;
  gdouble xgiven[8];
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

  /* For |eps_obs| >= 1, the psi-disc boundary point that maps to eps_obs / |eps_obs|. */
  if (outside)
  {
    const complex double psi_b = self->apply_shear (-h, eps_obs / eps_abs);

    xgiven[2 * ngiven]     = 1.0 - NC_GALAXY_SHAPE_FACTOR_QUAD_RHO_HINT_EPS;
    xgiven[2 * ngiven + 1] = carg (psi_b);
    ngiven++;
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
  arg.s_p             = 0.0;
  arg.s_N             = 0.0;

  /*
   * Peak of the integrand: the largest ln (pj) - d / 2s^2 over a coarse grid in psi and the
   * hints. Its two parts become the shifts s_p and s_N.
   */
  {
    gdouble best = -GSL_POSINF, pj, d;
    gint i;

    for (i = -ngiven; i < NC_GALAXY_SHAPE_FACTOR_QUAD_PEAK_NR * NC_GALAXY_SHAPE_FACTOR_QUAD_PEAK_NT; i++)
    {
      const gdouble rho   = (i < 0) ? xgiven[2 * (i + ngiven)] : (i / NC_GALAXY_SHAPE_FACTOR_QUAD_PEAK_NT + 0.5) / NC_GALAXY_SHAPE_FACTOR_QUAD_PEAK_NR;
      const gdouble theta = (i < 0) ? xgiven[2 * (i + ngiven) + 1] : -M_PI + 2.0 * M_PI * ((i % NC_GALAXY_SHAPE_FACTOR_QUAD_PEAK_NT) + 0.5) / NC_GALAXY_SHAPE_FACTOR_QUAD_PEAK_NT;
      gdouble score;

      _nc_galaxy_shape_factor_quad_point (&arg, rho, theta, &pj, &d);
      score = log (pj) - d / (2.0 * noise_var);

      if (isfinite (score) && (score > best))
      {
        best    = score;
        arg.s_p = log (pj);
        arg.s_N = d / (2.0 * noise_var);
      }
    }
  }

  ncm_integrate_2dim_divonne (&integ, 0.0, -M_PI, 1.0, M_PI,
                              self->reltol, 0.0, ngiven, 2, xgiven, &result, &error);

  /* P = e^{s_p - s_N} J / (2 pi s^2) */
  return arg.s_p - arg.s_N + log (result) - log (2.0 * M_PI * noise_var);
}

static gdouble
_nc_galaxy_shape_factor_quad_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return exp (_nc_galaxy_shape_factor_quad_eval_ln (NC_GALAXY_SHAPE_FACTOR_QUAD (gsf), pop, data->pop_data, data->std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2));
}

static gdouble
_nc_galaxy_shape_factor_quad_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return _nc_galaxy_shape_factor_quad_eval_ln (NC_GALAXY_SHAPE_FACTOR_QUAD (gsf), pop, data->pop_data, data->std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2);
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
  self->use_lp          = FALSE;
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
    case PROP_USE_LP:
    {
      NcGalaxyShapeFactorQuadPrivate * const self = nc_galaxy_shape_factor_quad_get_instance_private (gsfq);

      self->use_lp = g_value_get_boolean (value);
      break;
    }
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
    case PROP_USE_LP:
    {
      NcGalaxyShapeFactorQuadPrivate * const self = nc_galaxy_shape_factor_quad_get_instance_private (gsfq);

      g_value_set_boolean (value, self->use_lp);
      break;
    }
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

  /**
   * NcGalaxyShapeFactorQuad:use-lp:
   *
   * Integrates in log-polar coordinates around the noise peak recentered by a disc
   * automorphism, with the population's center excised for Beta. Temporary, for the
   * comparison against the Arb reference.
   */
  g_object_class_install_property (object_class,
                                   PROP_USE_LP,
                                   g_param_spec_boolean ("use-lp",
                                                         NULL,
                                                         "Integrate in log-polar coordinates around the noise peak",
                                                         FALSE,
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
  result = exp (_nc_galaxy_shape_factor_quad_eval_ln (gsfq, pop, pop_data, std_noise, g_1, g_2, epsilon_obs_1, epsilon_obs_2));
  nc_galaxy_shape_pop_data_unref (pop_data);

  return result;
}

