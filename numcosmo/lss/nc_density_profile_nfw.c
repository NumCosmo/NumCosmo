/***************************************************************************
 *            nc_density_profile_nfw.c
 *
 *  Tue June 10 16:40:38 2014
 *  Copyright  2014  
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_density_profile_nfw.c
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
 * SECTION:nc_density_profile_nfw
 * @title: NcDensityProfileNFW
 * @short_description: Density profile of Navarro-Frenk-White type.
 * 
 * This object implements the #NcDensityProfile class for a Navarro-Frenk-White (NFW) density profile.
 * 
 * The NFW profile is defined as
 * \begin{equation}
 * \rho(r) = \frac{\rho_s}{(r/r_s)(1 + r/r_s)^2}, 
 * \end{equation}
 * where $\rho_s$ is ... and $r_s$ is the scale radius, 
 * \begin{equation}
 * r_s \equiv \frac{r_{vir}}{c} = \left(\frac{3}{4\pi} \frac{M}{\Delta \overline{\rho}_m(z) c^3}\right)^{1/3},
 * \end{equation}
 * where $M$ is the virial mass, $\overline{\rho}_m(z)$ is the mean matter density, $\Delta$ is the overdensity 
 * parameter (as defined in #NcMultiplicityFunc). 
 *
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
 * References: astro-ph/0206508 and arxiv:1010.0744.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_density_profile_nfw.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include <gsl/gsl_sf_expint.h>
#include <math.h>

G_DEFINE_TYPE (NcDensityProfileNFW, nc_density_profile_nfw, NC_TYPE_DENSITY_PROFILE);

/**
 * nc_density_profile_nfw_new:
 *  
 * This function returns a #NcDensityProfile with a #NcDensityProfileNFW implementation.
 *
 * Returns: A new #NcDensityProfile.
 */
NcDensityProfile *
nc_density_profile_nfw_new ()
{
  return g_object_new (NC_TYPE_DENSITY_PROFILE_NFW, NULL);
} 

static gdouble
_nc_density_profile_nfw_scale_radius (NcHICosmo *model, const gdouble M, const gdouble z, const gdouble Delta)
{
  const gdouble Omega_m = nc_hicosmo_Omega_m (model);
  const gdouble v = 4.0 * M_PI / 3.0;
  gdouble rs_mnfw = cbrt (M / (Omega_m * v * Delta * ncm_c_crit_mass_density_solar_Mpc ())) / (1.0 + z);
  //printf("Delta = %.5g den_crit = %.5g omegam = %.5g M = %.5g\n", Delta, ncm_c_crit_mass_density_solar_Mpc(), Omega_m, M);
  //printf ("M = %.5g\n", 969.6 * Omega_m * Delta * ncm_c_crit_mass_density_solar_Mpc ());
  return rs_mnfw;
}

// trying to reproduce fig. 9 page 26 Cooray 2002 

static gdouble
_nc_density_profile_nfw_concentration_parameter (const gdouble M, const gdouble z)
{
  const gdouble M_star = 7.6 * 1.0e13; /* M value where sigma(M, z=0) = 1 */
  const gdouble c = 9.0 / (1.0 + z) * pow(M / M_star, -0.2);

  return c;
}

//end Cooray comparison

static gdouble
_nc_density_profile_nfw_eval_fourier (NcDensityProfile *dp, NcHICosmo *model, const gdouble k, const gdouble M, const gdouble z)
{
  NcDensityProfileNFW *dpnfw = NC_DENSITY_PROFILE_NFW (dp);
  gdouble c = _nc_density_profile_nfw_concentration_parameter (M, z);
  gdouble onepc = (1.0 + c);
  gdouble m_nfw = log(onepc) - c / onepc;
  gdouble factor_rs = 1.0 / c;
  gdouble rs = _nc_density_profile_nfw_scale_radius (model, M, z, dpnfw->Delta) * factor_rs;
  gdouble x = (1.0 + z) * k * rs;
  gdouble onepcx = onepc * x;
  gdouble u = 1.0 / m_nfw * (sin(x) * (gsl_sf_Si (onepcx) -  gsl_sf_Si (x)) + 
                             cos(x) * (gsl_sf_Ci (onepcx) -  gsl_sf_Ci (x)) -
                             sin(c * x) / onepcx);
  
  if (x == 0.0)
    return u = 1.0;
  
  return u;
} 

 
static void
nc_density_profile_nfw_init (NcDensityProfileNFW *dpnfw)
{
  dpnfw->Delta = 200.0;
}

static void
nc_density_profile_nfw_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_density_profile_nfw_parent_class)->finalize (object);
}

static void
nc_density_profile_nfw_class_init (NcDensityProfileNFWClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcDensityProfileClass *parent_class = NC_DENSITY_PROFILE_CLASS (klass);

  parent_class->eval_fourier = &_nc_density_profile_nfw_eval_fourier;
  object_class->finalize     = nc_density_profile_nfw_finalize;
}


