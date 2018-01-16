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

#include "lss/nc_density_profile_nfw.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/integral.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_expint.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcDensityProfileNFW, nc_density_profile_nfw, NC_TYPE_DENSITY_PROFILE);

#define VECTOR (NCM_MODEL (dpnfw)->params)
#define C_DELTA   (ncm_vector_get (VECTOR, NC_DENSITY_PROFILE_NFW_C_DELTA))
#define M_DELTA   (ncm_vector_get (VECTOR, NC_DENSITY_PROFILE_NFW_M_DELTA))

enum
{
  PROP_0,
  PROP_DELTA,
  PROP_SIZE,
};


static void
nc_density_profile_nfw_init (NcDensityProfileNFW *dpnfw)
{
  dpnfw->Delta = 200.0;
  dpnfw->r_Delta = 0.0;
}

static void
_nc_density_profile_nfw_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcDensityProfileNFW *dpnfw = NC_DENSITY_PROFILE_NFW (object);
  g_return_if_fail (NC_IS_DENSITY_PROFILE_NFW (object));

  switch (prop_id)
  {
    case PROP_DELTA:
      dpnfw->Delta = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_density_profile_nfw_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDensityProfileNFW *dpnfw = NC_DENSITY_PROFILE_NFW (object);
  g_return_if_fail (NC_IS_DENSITY_PROFILE_NFW (object));

  switch (prop_id)
  {
    case PROP_DELTA:
      g_value_set_double (value, dpnfw->Delta);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_density_profile_nfw_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_density_profile_nfw_parent_class)->finalize (object);
}

static gdouble _nc_density_profile_nfw_eval_density (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z);
static gdouble _nc_density_profile_nfw_integral_density_los (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);
static gdouble _nc_density_profile_nfw_integral_density_2d (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z);
static gdouble _nc_density_profile_nfw_eval_fourier (NcDensityProfile *dp, NcHICosmo *model, const gdouble k, const gdouble M, const gdouble z);

static void
nc_density_profile_nfw_class_init (NcDensityProfileNFWClass *klass)
{

  GObjectClass* object_class       = G_OBJECT_CLASS (klass);
  NcDensityProfileClass *parent_class = NC_DENSITY_PROFILE_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_density_profile_nfw_set_property;
  model_class->get_property = &_nc_density_profile_nfw_get_property;
  object_class->finalize    = &nc_density_profile_nfw_finalize;

  ncm_model_class_set_name_nick (model_class, "NFW Density Profile", "NFW");
  ncm_model_class_add_params (model_class, NC_DENSITY_PROFILE_NFW_SPARAM_LEN, 0, PROP_SIZE);

  /**
   * NcDensityProfileNFW:Delta:
   *
   * Consant that indicates the overdensity with respect to the critical density. 
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA,
                                   g_param_spec_double ("Delta",
                                                        NULL,
                                                        "Overdensity constant",
                                                        200, 1500, 200,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDensityProfileNFW:cDelta:
   * 
   * Concentration parameter.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_DENSITY_PROFILE_NFW_C_DELTA, "c_{\\Delta}", "cDelta",
                              0.5,  10.0, 1.0e-1,
                              NC_DENSITY_PROFILE_NFW_DEFAULT_PARAMS_ABSTOL, NC_DENSITY_PROFILE_NFW_DEFAULT_C_DELTA,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcDensityProfileNFW:MDelta:
   * 
   * Cluster mass within $R_\Delta$, where $\Delta$ is the overdensity.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_DENSITY_PROFILE_NFW_M_DELTA, "M_{\\Delta}", "MDelta",
                              1.0e13,  1.0e16, 1.0e13,
                              NC_DENSITY_PROFILE_NFW_DEFAULT_PARAMS_ABSTOL, NC_DENSITY_PROFILE_NFW_DEFAULT_M_DELTA,
                              NCM_PARAM_TYPE_FIXED);
  
  parent_class->eval_density         = &_nc_density_profile_nfw_eval_density;
  parent_class->integral_density_los = &_nc_density_profile_nfw_integral_density_los;
  parent_class->integral_density_2d  = &_nc_density_profile_nfw_integral_density_2d;
  parent_class->eval_fourier         = &_nc_density_profile_nfw_eval_fourier;
  
}

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

/// Old code: review it! //////////////////////////////

static gdouble
_nc_density_profile_nfw_scale_radius_matter (NcHICosmo *cosmo, const gdouble M, const gdouble z, const gdouble Delta)
{
  const gdouble rho_mz = nc_hicosmo_E2Omega_m (cosmo, z) * ncm_c_crit_mass_density_h2_solar_mass_Mpc3 ();
  const gdouble v = 4.0 * M_PI / 3.0;
  gdouble rs_mnfw = cbrt (M / (rho_mz * v * Delta ));
  //printf("Delta = %.5g den_crit = %.5g omegam = %.5g M = %.5g\n", Delta, ncm_c_crit_mass_density_h2_solar_mass_Mpc3(), Omega_m0, M);
  //printf ("M = %.5g\n", 969.6 * Omega_m0 * Delta * ncm_c_crit_mass_density_h2_solar_mass_Mpc3 ());
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
_nc_density_profile_nfw_eval_fourier (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble k, const gdouble M, const gdouble z)
{
  NcDensityProfileNFW *dpnfw = NC_DENSITY_PROFILE_NFW (dp);
  gdouble c = _nc_density_profile_nfw_concentration_parameter (M, z);
  gdouble onepc = (1.0 + c);
  gdouble m_nfw = log(onepc) - c / onepc;
  gdouble factor_rs = 1.0 / c;
  gdouble rs = _nc_density_profile_nfw_scale_radius_matter (cosmo, M, z, dpnfw->Delta) * factor_rs;
  gdouble x = (1.0 + z) * k * rs;
  gdouble onepcx = onepc * x;
  gdouble u = 1.0 / m_nfw * (sin(x) * (gsl_sf_Si (onepcx) -  gsl_sf_Si (x)) +
                             cos(x) * (gsl_sf_Ci (onepcx) -  gsl_sf_Ci (x)) -
                             sin(c * x) / onepcx);

  if (x == 0.0)
    return u = 1.0;

  return u;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////

static gdouble
_nc_density_profile_nfw_r_delta (NcDensityProfileNFW *dpnfw, NcHICosmo *cosmo, gdouble z)
{
  gdouble rho_c = nc_hicosmo_crit_density (cosmo) * nc_hicosmo_E2 (cosmo, z);
  gdouble rD3 = 3.0 * M_DELTA / (4.0 * M_PI * dpnfw->Delta * rho_c);
  gdouble r_Delta = cbrt (rD3);
    
  return r_Delta; 
}

static gdouble
_nc_density_profile_nfw_scale_radius (NcDensityProfileNFW *dpnfw, NcHICosmo *cosmo, gdouble z)
{
  gdouble r_Delta = _nc_density_profile_nfw_r_delta (dpnfw, cosmo, z);
    
  return r_Delta / C_DELTA;
}

static gdouble 
_nc_density_profile_nfw_deltac (NcDensityProfileNFW *dpnfw)
{
  gdouble onepc = 1.0 + C_DELTA;
  gdouble c2    = C_DELTA * C_DELTA;
  gdouble c3    = c2 * C_DELTA;

  gdouble delta_c = (dpnfw->Delta / 3.0) * c3 / (log(onepc) - C_DELTA /onepc); 

  return delta_c;
}

static gdouble
_nc_density_profile_nfw_eval_density (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z)
{
  NcDensityProfileNFW *dpnfw = NC_DENSITY_PROFILE_NFW (dp);

  gdouble rs      = _nc_density_profile_nfw_scale_radius (dpnfw, cosmo, z);
  gdouble x       = r / rs;
  gdouble delta_c = _nc_density_profile_nfw_deltac (dpnfw);
  gdouble rho_c   = nc_hicosmo_crit_density (cosmo) * nc_hicosmo_E2 (cosmo, z);
  gdouble onepx   = 1.0 + x;
  gdouble onepx2  = onepx * onepx; 
    
  return delta_c * rho_c/(x * onepx2);
}

/* los = line of sight */
static gdouble 
_nc_density_profile_nfw_integral_density_los (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z)
{
  NcDensityProfileNFW *dpnfw = NC_DENSITY_PROFILE_NFW (dp);

  gdouble rho_c   = nc_hicosmo_crit_density (cosmo) * nc_hicosmo_E2 (cosmo, z);
  gdouble delta_c = _nc_density_profile_nfw_deltac (dpnfw);
  gdouble rs      = _nc_density_profile_nfw_scale_radius (dpnfw, cosmo, z); 
  gdouble A       = rs * delta_c * rho_c;
  gdouble x       = R / rs;
  gdouble x2      = x * x;
  
  if (x < 1.0)
    return A / (x2 - 1.0) * (1.0 - 2.0/sqrt(1.0 - x2) * gsl_atanh(sqrt((1.0 - x)/(1.0 + x))));

  else if (x == 1.0)
    return A / 3.0;

  else
    return A / (x2 - 1.0) * (1.0 - 2.0/sqrt(x2 - 1.0) * gsl_atanh(sqrt((x - 1.0)/(1.0 + x))));
  
}

static gdouble 
_nc_density_profile_nfw_integral_density_2d (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z)
{
  NcDensityProfileNFW *dpnfw = NC_DENSITY_PROFILE_NFW (dp);

  gdouble rho_c   = nc_hicosmo_crit_density (cosmo) * nc_hicosmo_E2 (cosmo, z);
  gdouble delta_c = _nc_density_profile_nfw_deltac (dpnfw);
  gdouble rs      = _nc_density_profile_nfw_scale_radius (dpnfw, cosmo, z); 
  gdouble A       = rs * delta_c * rho_c;
  gdouble x       = R / rs;
  gdouble x2      = x * x;
    
  if (x < 1.0)
    return 2.0 * A * (2.0/sqrt(1.0 - x2) * gsl_atanh(sqrt((1.0 - x)/(1.0 + x))) + log(x) - M_LN2);

  else if (x == 1.0)
    return 2.0 * A * x2 * (1.0 - M_LN2);

  else
    return 2.0 * A * (2.0/sqrt(x2 - 1.0) * gsl_atanh(sqrt((x - 1.0)/(1.0 + x))) + log(x) - M_LN2);
  
}
