/***************************************************************************
 *            nc_xcor_kernel_radial_kdep.c
 *
 *  Thu August 21 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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

/**
 * NcXcorKernelRadialKDep:
 *
 * Closed-form scale dependence that a #NcXcorKernelRadial can carry.
 *
 * Attached to a kernel, it multiplies the radial integrand,
 * \begin{equation}
 *   K(\chi,k) = W(\chi)\,g(\chi,k)\,\sqrt{P(k, z=0)} ,
 * \end{equation}
 *
 * which is what makes the integrand genuinely a function of both variables. A
 * kernel without one has $g \equiv 1$, and its radial integral is the same
 * shape at every $k$ up to an overall factor.
 *
 * That separability is not a detail: a method that reduces the radial integral
 * to a single transform of $W$ can only produce a $k$-independent shape, so it
 * cannot represent this at all. Any shape can carry any scale dependence, so
 * the two vary independently in a benchmark rather than being baked into
 * separate kernel classes.
 *
 * Note that a code compared against this must be able to evaluate a
 * $k$-dependent kernel at all; several cannot, and that is the claim the
 * construction exists to exercise rather than a limitation of it.
 *
 */

/**
 * NcXcorKernelRadialKDepGrowth:
 *
 * Scale-dependent growth across a smooth transition in $k$.
 *
 * \begin{equation}
 *   g(\chi,k) = \exp\left[-\alpha\,\sigma(k)\,
 *   \frac{\chi_\mathrm{ref} - \chi}{\chi_\mathrm{ref}}\right] ,
 *   \qquad
 *   \sigma(k) = \frac{(k/k_t)^2}{1 + (k/k_t)^2} .
 * \end{equation}
 *
 * Below $k_t$ the growth is scale-free, $g \to 1$; above it the effect
 * saturates at $\alpha$ and accumulates along the line of sight, in the manner
 * of free-streaming or a modified growth rate. It is an emulation with the
 * right structure, closed form and exactly evaluable, not a fit to any
 * particular model: $\alpha = 0$ recovers the scale-independent case exactly.
 *
 * $g = 1$ at $\chi_\mathrm{ref}$ and departs from it on both sides --
 * suppression at smaller $\chi$, enhancement at larger. Setting
 * $\chi_\mathrm{ref}$ at or beyond the outer edge of the window therefore makes
 * it a pure suppression across the whole support. The factor is deliberately
 * left monotone and unclamped: bounding it would introduce a kink, and a kink
 * inside an integration domain is exactly what these kernels exist to avoid.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "nc/xcor/nc_xcor_kernel_radial_kdep.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcXcorKernelRadialKDep, nc_xcor_kernel_radial_kdep, G_TYPE_OBJECT)

static void
nc_xcor_kernel_radial_kdep_init (NcXcorKernelRadialKDep *kdep)
{
}

static void
nc_xcor_kernel_radial_kdep_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_radial_kdep_parent_class)->finalize (object);
}

static void
nc_xcor_kernel_radial_kdep_class_init (NcXcorKernelRadialKDepClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = &nc_xcor_kernel_radial_kdep_finalize;
}

/**
 * nc_xcor_kernel_radial_kdep_ref:
 * @kdep: a #NcXcorKernelRadialKDep
 *
 * Increases the reference count of @kdep by one.
 *
 * Returns: (transfer full): @kdep
 */
NcXcorKernelRadialKDep *
nc_xcor_kernel_radial_kdep_ref (NcXcorKernelRadialKDep *kdep)
{
  return g_object_ref (kdep);
}

/**
 * nc_xcor_kernel_radial_kdep_free:
 * @kdep: a #NcXcorKernelRadialKDep
 *
 * Decreases the reference count of @kdep by one.
 */
void
nc_xcor_kernel_radial_kdep_free (NcXcorKernelRadialKDep *kdep)
{
  g_object_unref (kdep);
}

/**
 * nc_xcor_kernel_radial_kdep_clear:
 * @kdep: a #NcXcorKernelRadialKDep
 *
 * If *@kdep is not %NULL, decreases its reference count by one and sets
 * *@kdep to %NULL.
 */
void
nc_xcor_kernel_radial_kdep_clear (NcXcorKernelRadialKDep **kdep)
{
  g_clear_object (kdep);
}

/**
 * nc_xcor_kernel_radial_kdep_eval: (virtual eval)
 * @kdep: a #NcXcorKernelRadialKDep
 * @chi: comoving distance $\chi$ in Mpc
 * @k: wavenumber $k$ in $\mathrm{Mpc}^{-1}$
 *
 * Evaluates the scale-dependent factor $g(\chi,k)$ multiplying the radial
 * integrand.
 *
 * Returns: the value of $g(\chi,k)$, dimensionless.
 */
gdouble
nc_xcor_kernel_radial_kdep_eval (NcXcorKernelRadialKDep *kdep, gdouble chi, gdouble k)
{
  return NC_XCOR_KERNEL_RADIAL_KDEP_GET_CLASS (kdep)->eval (kdep, chi, k);
}

/*
 * Growth transition
 */

struct _NcXcorKernelRadialKDepGrowth
{
  /*< private >*/
  NcXcorKernelRadialKDep parent_instance;

  gdouble amplitude;
  gdouble k_transition;
  gdouble chi_ref;
};

enum
{
  PROP_GROWTH_0,
  PROP_AMPLITUDE,
  PROP_K_TRANSITION,
  PROP_CHI_REF,
  PROP_GROWTH_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelRadialKDepGrowth, nc_xcor_kernel_radial_kdep_growth, NC_TYPE_XCOR_KERNEL_RADIAL_KDEP)

static void
nc_xcor_kernel_radial_kdep_growth_init (NcXcorKernelRadialKDepGrowth *kdepg)
{
  kdepg->amplitude    = 0.0;
  kdepg->k_transition = 0.0;
  kdepg->chi_ref      = 0.0;
}

static void
nc_xcor_kernel_radial_kdep_growth_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelRadialKDepGrowth *kdepg = NC_XCOR_KERNEL_RADIAL_KDEP_GROWTH (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_RADIAL_KDEP_GROWTH (object));

  switch (prop_id)
  {
    case PROP_AMPLITUDE:
      kdepg->amplitude = g_value_get_double (value);
      break;
    case PROP_K_TRANSITION:
      kdepg->k_transition = g_value_get_double (value);
      break;
    case PROP_CHI_REF:
      kdepg->chi_ref = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_xcor_kernel_radial_kdep_growth_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelRadialKDepGrowth *kdepg = NC_XCOR_KERNEL_RADIAL_KDEP_GROWTH (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_RADIAL_KDEP_GROWTH (object));

  switch (prop_id)
  {
    case PROP_AMPLITUDE:
      g_value_set_double (value, kdepg->amplitude);
      break;
    case PROP_K_TRANSITION:
      g_value_set_double (value, kdepg->k_transition);
      break;
    case PROP_CHI_REF:
      g_value_set_double (value, kdepg->chi_ref);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static gdouble _nc_xcor_kernel_radial_kdep_growth_eval (NcXcorKernelRadialKDep *kdep, gdouble chi, gdouble k);

static void
nc_xcor_kernel_radial_kdep_growth_class_init (NcXcorKernelRadialKDepGrowthClass *klass)
{
  GObjectClass *object_class                = G_OBJECT_CLASS (klass);
  NcXcorKernelRadialKDepClass *parent_class = NC_XCOR_KERNEL_RADIAL_KDEP_CLASS (klass);

  object_class->set_property = &nc_xcor_kernel_radial_kdep_growth_set_property;
  object_class->get_property = &nc_xcor_kernel_radial_kdep_growth_get_property;

  /**
   * NcXcorKernelRadialKDepGrowth:amplitude:
   *
   * Saturated suppression $\alpha$ reached well above the transition. Zero
   * recovers the scale-independent case exactly.
   */
  g_object_class_install_property (object_class,
                                   PROP_AMPLITUDE,
                                   g_param_spec_double ("amplitude",
                                                        NULL,
                                                        "Saturated suppression",
                                                        0.0, G_MAXDOUBLE, 0.1,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelRadialKDepGrowth:k-transition:
   *
   * Wavenumber $k_t$ of the transition, in $\mathrm{Mpc}^{-1}$.
   */
  g_object_class_install_property (object_class,
                                   PROP_K_TRANSITION,
                                   g_param_spec_double ("k-transition",
                                                        NULL,
                                                        "Transition wavenumber in 1/Mpc",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 0.05,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelRadialKDepGrowth:chi-ref:
   *
   * Comoving distance $\chi_\mathrm{ref}$ at which the suppression vanishes,
   * in Mpc. It sets where along the line of sight the factor is normalized to
   * one, and how fast it accumulates towards the observer.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHI_REF,
                                   g_param_spec_double ("chi-ref",
                                                        NULL,
                                                        "Reference comoving distance in Mpc",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 1500.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  parent_class->eval = &_nc_xcor_kernel_radial_kdep_growth_eval;
}

/**
 * nc_xcor_kernel_radial_kdep_growth_new:
 * @amplitude: saturated suppression $\alpha$
 * @k_transition: transition wavenumber $k_t$, in $\mathrm{Mpc}^{-1}$
 * @chi_ref: reference comoving distance $\chi_\mathrm{ref}$, in Mpc
 *
 * Creates a new #NcXcorKernelRadialKDepGrowth.
 *
 * Returns: (transfer full): a new #NcXcorKernelRadialKDepGrowth
 */
NcXcorKernelRadialKDepGrowth *
nc_xcor_kernel_radial_kdep_growth_new (gdouble amplitude, gdouble k_transition, gdouble chi_ref)
{
  NcXcorKernelRadialKDepGrowth *kdepg = g_object_new (NC_TYPE_XCOR_KERNEL_RADIAL_KDEP_GROWTH,
                                                      "amplitude", amplitude,
                                                      "k-transition", k_transition,
                                                      "chi-ref", chi_ref,
                                                      NULL);

  return kdepg;
}

/**
 * nc_xcor_kernel_radial_kdep_growth_get_amplitude:
 * @kdepg: a #NcXcorKernelRadialKDepGrowth
 *
 * Returns: the saturated suppression $\alpha$.
 */
gdouble
nc_xcor_kernel_radial_kdep_growth_get_amplitude (NcXcorKernelRadialKDepGrowth *kdepg)
{
  return kdepg->amplitude;
}

/**
 * nc_xcor_kernel_radial_kdep_growth_get_k_transition:
 * @kdepg: a #NcXcorKernelRadialKDepGrowth
 *
 * Returns: the transition wavenumber $k_t$, in $\mathrm{Mpc}^{-1}$.
 */
gdouble
nc_xcor_kernel_radial_kdep_growth_get_k_transition (NcXcorKernelRadialKDepGrowth *kdepg)
{
  return kdepg->k_transition;
}

/**
 * nc_xcor_kernel_radial_kdep_growth_get_chi_ref:
 * @kdepg: a #NcXcorKernelRadialKDepGrowth
 *
 * Returns: the reference comoving distance $\chi_\mathrm{ref}$, in Mpc.
 */
gdouble
nc_xcor_kernel_radial_kdep_growth_get_chi_ref (NcXcorKernelRadialKDepGrowth *kdepg)
{
  return kdepg->chi_ref;
}

static gdouble
_nc_xcor_kernel_radial_kdep_growth_eval (NcXcorKernelRadialKDep *kdep, gdouble chi, gdouble k)
{
  NcXcorKernelRadialKDepGrowth *kdepg = NC_XCOR_KERNEL_RADIAL_KDEP_GROWTH (kdep);
  const gdouble x                     = k / kdepg->k_transition;
  const gdouble x2                    = x * x;
  const gdouble sigma                 = x2 / (1.0 + x2);

  return exp (-kdepg->amplitude * sigma * (kdepg->chi_ref - chi) / kdepg->chi_ref);
}

