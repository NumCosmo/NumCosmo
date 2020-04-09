/***************************************************************************
 *            nc_halo_density_profile_nfw.c
 *
 *  Tue June 10 16:40:38 2014
 *  Copyright  2014
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_halo_density_profile_nfw.c
 * Copyright (C) 2014 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_halo_density_profile_nfw
 * @title: NcHaloDensityProfileNFW
 * @short_description: Density profile of Navarro-Frenk-White type.
 *
 * This object implements the #NcHaloDensityProfile class for a Navarro-Frenk-White (NFW) density profile.
 *
 * The NFW profile is defined as
 * \begin{equation}
 * \rho(r) = \frac{\delta_c \rho_{crit}}{(r/r_s)(1 + r/r_s)^2},
 * \end{equation}
 * where $\rho_{crit} (z) = \frac{3 H^2(z)}{8\pi G} [M_\odot / Mpc^3]$,
 * \begin{equation}
 * \delta_c = \frac{\Delta}{3} \frac{c^3}{\ln (1 + c) - \frac{c}{1 + c}},
 * \end{equation}
 * $c$ is the concentration parameter and $r_s$ is the scale radius,
 * \begin{equation}
 * r_s [Mpc] \equiv \frac{r_{\Delta}}{c} = \left(\frac{3}{4\pi} \frac{M}{\Delta \rho_{crit}(z) c^3}\right)^{1/3},
 * \end{equation}
 * where $M$ is the halo mass $[M_\odot]$, $\Delta$ is the overdensity parameter (as defined in #NcMultiplicityFunc).
 *
 * FIXME
 * The normalized NFW density profile ($u_M(r) = \rho(r) / M$) in the Fourier space is given by
 * \begin{equation}
 * \tilde{u}_M(k) = \frac{1}{m_{nfw}(c)} \left[ \sin(x) \left[\text{Si}((1+c)x) - \text{Si}(x) \right] + \cos(x) \left[\text{Ci}((1+c)x) - \text{Ci}(x) \right] - \frac{\sin(cx)}{(1+c)x} \right],
 * \end{equation}
 * where $x \equiv (1+z)kr_s$, and $\text{Si}(x)$ and $\text{Ci}(x)$ are the sine and cosine integrals, namely.
 * \begin{equation}
 * \text{Si}(x) = \int_0^x \frac{\sin(t)}{t} dt \quad \text{and} \quad \text{Ci}(x) = - \int_x^\infty \frac{\cos(t)}{t} dt.
 * \end{equation}
 *
 * The concentration parameter is (change this!)
 * \begin{equation}
 * c(M, z) = A_{vir} \left( \frac{M}{2 \times 10^{12} \text{h}^{-1}M_{\odot}}\right)^{B_{vir}} (1+z)^{C_{vir}}.
 * \end{equation}
 *
 * References: [Navarro (1996)][XNavarro1996], [Wright (2000)][XWright2000], astro-ph/0206508 and arxiv:1010.0744.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_halo_density_profile_nfw.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/integral.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_expint.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcHaloDensityProfileNFW, nc_halo_density_profile_nfw, NC_TYPE_HALO_DENSITY_PROFILE);

#define VECTOR  (NCM_MODEL (dpnfw)->params)
#define M_DELTA (ncm_vector_get (VECTOR, NC_HALO_DENSITY_PROFILE_M_DELTA))
#define C_DELTA (ncm_vector_get (VECTOR, NC_HALO_DENSITY_PROFILE_C_DELTA))

enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_halo_density_profile_nfw_init (NcHaloDensityProfileNFW *dpnfw)
{
}

static void
_nc_halo_density_profile_nfw_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /* NcHaloDensityProfileNFW *dpnfw = NC_HALO_DENSITY_PROFILE_NFW (object); */
  g_return_if_fail (NC_IS_HALO_DENSITY_PROFILE_NFW (object));
  
  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_density_profile_nfw_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcHaloDensityProfileNFW *dpnfw = NC_HALO_DENSITY_PROFILE_NFW (object); */
  g_return_if_fail (NC_IS_HALO_DENSITY_PROFILE_NFW (object));
  
  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_density_profile_nfw_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_density_profile_nfw_parent_class)->finalize (object);
}

static gdouble _nc_halo_density_profile_nfw_eval_dl_density (NcHaloDensityProfile *dp, const gdouble x);
static gdouble _nc_halo_density_profile_nfw_eval_dl_spher_mass (NcHaloDensityProfile *dp, const gdouble x);
static gdouble _nc_halo_density_profile_nfw_eval_dl_2d_density (NcHaloDensityProfile *dp, const gdouble X);
static gdouble _nc_halo_density_profile_nfw_eval_dl_cyl_mass (NcHaloDensityProfile *dp, const gdouble X);

static void
nc_halo_density_profile_nfw_class_init (NcHaloDensityProfileNFWClass *klass)
{
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);
  NcHaloDensityProfileClass *dp_class = NC_HALO_DENSITY_PROFILE_CLASS (klass);
  NcmModelClass *model_class          = NCM_MODEL_CLASS (klass);
  
  model_class->set_property = &_nc_halo_density_profile_nfw_set_property;
  model_class->get_property = &_nc_halo_density_profile_nfw_get_property;
  object_class->finalize    = &_nc_halo_density_profile_nfw_finalize;
  
  ncm_model_class_set_name_nick (model_class, "NFW Density Profile", "NFW");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  dp_class->eval_dl_density    = &_nc_halo_density_profile_nfw_eval_dl_density;
  dp_class->eval_dl_spher_mass = &_nc_halo_density_profile_nfw_eval_dl_spher_mass;
  dp_class->eval_dl_2d_density = &_nc_halo_density_profile_nfw_eval_dl_2d_density;
  dp_class->eval_dl_cyl_mass   = &_nc_halo_density_profile_nfw_eval_dl_cyl_mass;
}

static gdouble 
_nc_halo_density_profile_nfw_eval_dl_density (NcHaloDensityProfile *dp, const gdouble x)
{
  return 1.0 / (x * gsl_pow_2 (1.0 + x));
}

static gdouble 
_nc_halo_density_profile_nfw_eval_dl_spher_mass (NcHaloDensityProfile *dp, const gdouble x)
{
  return log1p (x) - x / (1.0 + x);
}

static gdouble 
_nc_halo_density_profile_nfw_eval_dl_2d_density (NcHaloDensityProfile *dp, const gdouble X)
{
  const gdouble Xm1           = X - 1.0;
  const gdouble Xp1           = X + 1.0;
  const gdouble X2m1          = Xm1 * Xp1;
  const gdouble sqrt_abs_X2m1 = sqrt (fabs (X2m1));

  if (fabs (Xm1) < 1.0e-6)
  {
    return 2.0 * (1.0 / 3.0 + (-2.0 / 5.0 + (13.0 / 35.0 - 20.0 / 63.0 * Xm1) * Xm1) * Xm1);
  }
  else if (Xm1 < 0.0)
  {
    if (X < 1.0e-6)
    {
      const gdouble X_2     = 0.5 * X;
      const gdouble log_X_2 = log (X_2);
      return 2.0 / X2m1 * ((1.0 + log_X_2) + 0.25 * (1.0 + 2.0 * log_X_2) * X_2 * X_2 + (7.0 / 2.0 + 6.0 * log_X_2) * X_2 * X_2 * X_2 * X_2);
    }
    else
      return 2.0 / X2m1 * (1.0 - atanh (sqrt_abs_X2m1) / sqrt_abs_X2m1);
  }
  else
  {
    return 2.0 / X2m1 * (1.0 - atan (sqrt_abs_X2m1) / sqrt_abs_X2m1);
  }
}

static gdouble 
_nc_halo_density_profile_nfw_eval_dl_cyl_mass (NcHaloDensityProfile *dp, const gdouble X)
{
  const gdouble Xm1           = X - 1.0;
  const gdouble Xp1           = X + 1.0;
  const gdouble X2m1          = Xm1 * Xp1;
  const gdouble sqrt_abs_X2m1 = sqrt (fabs (X2m1));

  if (fabs (Xm1) < 1.0e-6)
  {
    return 2.0 * (1.0 - M_LN2) + (2.0 / 3.0 + (-1.0 / 15.0 - 2.0 / 105.0 * Xm1) * Xm1) * Xm1;
  }
  else if (Xm1 < 0.0)
  {
    return 2.0 * (atanh (sqrt_abs_X2m1) / sqrt_abs_X2m1 + log (0.5 * X));
  }
  else
  {
    return 2.0 * (atan  (sqrt_abs_X2m1) / sqrt_abs_X2m1 + log (0.5 * X));
  }

}

/**
 * nc_halo_density_profile_nfw_new:
 * @mdef: a #NcHaloDensityProfileMassDef
 * @Delta: cluster threshold mass definition $\Delta$
 *
 * This function returns a #NcHaloDensityProfile with a #NcHaloDensityProfileNFW implementation.
 *
 * Returns: A new #NcHaloDensityProfileNFW.
 */
NcHaloDensityProfileNFW *
nc_halo_density_profile_nfw_new (const NcHaloDensityProfileMassDef mdef, const gdouble Delta)
{
  NcHaloDensityProfileNFW *dp_nfw = g_object_new (NC_TYPE_HALO_DENSITY_PROFILE_NFW,
                                                  "mass-def", mdef,
                                                  "Delta",    Delta,
                                                  NULL);

  return dp_nfw;
}
