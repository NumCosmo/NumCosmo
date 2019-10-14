/***************************************************************************
 *            nc_density_profile_einasto.c
 *
 *  Wed July 17 12:33:27 2019
 *  Copyright  2014
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_density_profile_einasto.c
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
 * SECTION:nc_density_profile_einasto
 * @title: NcDensityProfileEinasto
 * @short_description: Density profile of Einasto type.
 *
 * This object implements the #NcDensityProfile class for the Einasto density profile.
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

#include "lss/nc_density_profile_einasto.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/integral.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_expint.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcDensityProfileEinasto, nc_density_profile_einasto, NC_TYPE_DENSITY_PROFILE);

#define VECTOR (NCM_MODEL (dpe)->params)
#define RHOS  (ncm_vector_get (VECTOR, NC_DENSITY_PROFILE_EINASTO_RHOS))
#define RS    (ncm_vector_get (VECTOR, NC_DENSITY_PROFILE_EINASTO_RS))
#define ALPHA (ncm_vector_get (VECTOR, NC_DENSITY_PROFILE_EINASTO_ALPHA))

enum
{
  PROP_0,
  PROP_DELTA,
  PROP_SIZE,
};


static void
nc_density_profile_einasto_init (NcDensityProfileEinasto *dpe)
{
  dpe->Delta = 200.0;
  dpe->r_Delta = 0.0;
}

static void
_nc_density_profile_einasto_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcDensityProfileEinasto *dpnfw = NC_DENSITY_PROFILE_EINASTO (object);
  g_return_if_fail (NC_IS_DENSITY_PROFILE_EINASTO (object));

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
_nc_density_profile_einasto_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDensityProfileEinasto *dpnfw = NC_DENSITY_PROFILE_EINASTO (object);
  g_return_if_fail (NC_IS_DENSITY_PROFILE_EINASTO (object));

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
nc_density_profile_einasto_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_density_profile_einasto_parent_class)->finalize (object);
}

static gdouble _nc_density_profile_einasto_eval_density (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z);
static gdouble _nc_density_profile_einasto_integral_density_los (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z);
static gdouble _nc_density_profile_einasto_integral_density_2d (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z);
static gdouble _nc_density_profile_einasto_eval_fourier (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble k, const gdouble M, const gdouble z);
static gdouble _nc_density_profile_einasto_scale_radius (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble z);

static void
nc_density_profile_einasto_class_init (NcDensityProfileEinastoClass *klass)
{
  GObjectClass* object_class          = G_OBJECT_CLASS (klass);
  NcDensityProfileClass *parent_class = NC_DENSITY_PROFILE_CLASS (klass);
  NcmModelClass *model_class          = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_density_profile_einasto_set_property;
  model_class->get_property = &_nc_density_profile_einasto_get_property;
  object_class->finalize    = &nc_density_profile_einasto_finalize;

  ncm_model_class_set_name_nick (model_class, "Einasto Density Profile", "Einasto");
  ncm_model_class_add_params (model_class, NC_DENSITY_PROFILE_EINASTO_SPARAM_LEN, 0, PROP_SIZE);

  /**
   * NcDensityProfileEinasto:Delta:
   *
   * Constant that indicates the overdensity with respect to the critical density. 
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
   * NcDensityProfileEinasto:rhos:
   * 
   * Central density.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_DENSITY_PROFILE_EINASTO_RHOS, "\rho_{s}", "rhos",
                              0.5,  10.0, 1.0e-1,
                              NC_DENSITY_PROFILE_EINASTO_DEFAULT_PARAMS_ABSTOL, NC_DENSITY_PROFILE_EINASTO_DEFAULT_RHOS,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcDensityProfileEinasto:rs:
   * 
   * Scale radius $r_s$.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_DENSITY_PROFILE_EINASTO_RS, "r_{s}", "rs",
                              1.0e13,  1.0e16, 1.0e13,
                              NC_DENSITY_PROFILE_EINASTO_DEFAULT_PARAMS_ABSTOL, NC_DENSITY_PROFILE_EINASTO_DEFAULT_RS,
                              NCM_PARAM_TYPE_FIXED);
  
  parent_class->eval_density         = &_nc_density_profile_einasto_eval_density;
  parent_class->integral_density_los = &_nc_density_profile_einasto_integral_density_los;
  parent_class->integral_density_2d  = &_nc_density_profile_einasto_integral_density_2d;
  parent_class->eval_fourier         = &_nc_density_profile_einasto_eval_fourier;
	parent_class->scale_radius         = &_nc_density_profile_einasto_scale_radius;
  
}

/**
 * nc_density_profile_einasto_new:
 *
 * This function returns a #NcDensityProfile with a #NcDensityProfileEinasto implementation.
 *
 * Returns: A new #NcDensityProfile.
 */
NcDensityProfile *
nc_density_profile_einasto_new ()
{
  return g_object_new (NC_TYPE_DENSITY_PROFILE_EINASTO, NULL);
}

/// Old code: review it! //////////////////////////////

static gdouble
_nc_density_profile_einasto_scale_radius_matter (NcHICosmo *cosmo, const gdouble M, const gdouble z, const gdouble Delta)
{
  const gdouble rho_mz = nc_hicosmo_E2Omega_m (cosmo, z) * ncm_c_crit_mass_density_h2_solar_mass_Mpc3 ();
  const gdouble v = 4.0 * M_PI / 3.0;
  gdouble rs_mnfw = cbrt (M / (rho_mz * v * Delta ));
  //printf("z = % 22.15g Delta = %.5g den_crit = %.5g omegam = %.5g M = %.5g\n", z, Delta, ncm_c_crit_mass_density_h2_solar_mass_Mpc3(), nc_hicosmo_E2Omega_m (cosmo, z), M);
  //printf ("M = %.5g\n", 969.6 * Omega_m0 * Delta * ncm_c_crit_mass_density_h2_solar_mass_Mpc3 ());
  return rs_mnfw;
}

// trying to reproduce fig. 9 page 26 Cooray 2002

static gdouble
_nc_density_profile_einasto_concentration_parameter (const gdouble M, const gdouble z)
{
  const gdouble M_star = 7.6 * 1.0e13; /* M value where sigma(M, z=0) = 1 */
  const gdouble c = 9.0 / (1.0 + z) * pow(M / M_star, -0.2);

  return c;
}

//end Cooray comparison

static gdouble
_nc_density_profile_einasto_eval_fourier (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble k, const gdouble M, const gdouble z)
{
  NcDensityProfileEinasto *dpnfw = NC_DENSITY_PROFILE_EINASTO (dp);
  gdouble c = _nc_density_profile_einasto_concentration_parameter (M, z);
  gdouble onepc = (1.0 + c);
  gdouble m_nfw = log(onepc) - c / onepc;
  gdouble factor_rs = 1.0 / c;
  gdouble rs = _nc_density_profile_einasto_scale_radius_matter (cosmo, M, z, dpnfw->Delta) * factor_rs;
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
_rho_crit_solar_mass_Mpc3 (NcHICosmo *cosmo, gdouble z)
{
	/* Critical density in M_solar / Mpc^3 units */
  return ncm_c_crit_mass_density_h2_solar_mass_Mpc3 () * nc_hicosmo_h2 (cosmo) * nc_hicosmo_E2 (cosmo, z);
}

static gdouble
_nc_density_profile_einasto_r_delta (NcDensityProfileEinasto *dpnfw, NcHICosmo *cosmo, gdouble z)
{
  const gdouble rho_c   = _rho_crit_solar_mass_Mpc3 (cosmo, z);
  const gdouble rD3     = 3.0 / (4.0 * M_PI * dpnfw->Delta * rho_c);
  const gdouble r_Delta = cbrt (rD3);

  return r_Delta; 
}

static gdouble
_nc_density_profile_einasto_scale_radius (NcDensityProfile *dp, NcHICosmo *cosmo, gdouble z)
{
	NcDensityProfileEinasto *dpnfw = NC_DENSITY_PROFILE_EINASTO (dp);
  gdouble r_Delta = _nc_density_profile_einasto_r_delta (dpnfw, cosmo, z);

	/* IMPLEMENT if it is the case!!! */
	return r_Delta;
}

static gdouble 
_nc_density_profile_einasto_deltac (NcDensityProfileEinasto *dpnfw)
{
  /*gdouble onepc = 1.0;*/
  gdouble delta_c = (dpnfw->Delta / 3.0); 
	/* IMPLEMENT if it is the case! */

  return delta_c;
}

////////////////////////////////// Cluster toolkit ///////////////////////////////
#define rhomconst 2.77533742639e+11

int calc_xi_einasto(double*r, int Nr, double Mass, double conc, int delta, double om, double*xi_nfw){
  int i;
  double rhom = om * ncm_c_crit_mass_density_h2_solar_mass_Mpc3() * 0.7 * 0.7; //om*rhomconst;//SM h^2/Mpc^3
  //double rho0_rhom = delta/(3.*(log(1.+conc)-conc/(1.+conc)));
  double rdelta = pow(Mass/(1.33333333333*M_PI*rhom*delta), 0.33333333333);
  double rscale = rdelta/conc;
  double fc = log(1.+conc)-conc/(1.+conc);
  double r_rs = 0.0;
  for(i = 0; i < Nr; i++){
    r_rs = r[i]/rscale;
    //xi_nfw[i] = rho0_rhom/(r_rs*(1+r_rs)*(1+r_rs)) - 1.;
    xi_nfw[i] = Mass/(4.*M_PI*rscale*rscale*rscale*fc)/(r_rs*(1+r_rs)*(1+r_rs))/rhom - 1.0;

  }
	printf ("factor dele % 22.15g\n", Mass/(4.*M_PI*rscale*rscale*rscale*fc)/(r_rs*(1+r_rs)*(1+r_rs)));
	//printf ("rscale dele % 22.15g | % 22.15g % 22.15g % 22.15g %d % 22.15g\n", rscale, rdelta, conc, rhom, delta, Mass);
  return 0;
}

int calc_rho_einasto(double*r, int Nr, double Mass, double conc, int delta, double Omega_m, double*rho_nfw){
  int i;
  double rhom = Omega_m * ncm_c_crit_mass_density_h2_solar_mass_Mpc3() * 0.7 * 0.7; //Omega_m*rhomconst;//Msun h^2/Mpc^3
  calc_xi_einasto(r, Nr, Mass, conc, delta, Omega_m, rho_nfw); //rho_nfw actually holds xi_nfw here
  for(i = 0; i < Nr; i++){
    rho_nfw[i] = rhom*(1+rho_nfw[i]);
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

static gdouble
_nc_density_profile_einasto_eval_density (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z)
{
  NcDensityProfileEinasto *dpe = NC_DENSITY_PROFILE_EINASTO (dp);

  gdouble rs      = _nc_density_profile_einasto_scale_radius (dp, cosmo, z);
  gdouble x       = r / rs;
  gdouble delta_c = _nc_density_profile_einasto_deltac (dpe);
  gdouble rho_c   = _rho_crit_solar_mass_Mpc3 (cosmo, z);
  gdouble onepx   = 1.0 + x;
  gdouble onepx2  = onepx * onepx; 

	{
		const gdouble Omega_m0 = 1.0;//nc_hicosmo_Omega_m0 (cosmo);
		gdouble ret = 0.0;
		gdouble rrr = r;

		/* FIX this function!!! Copy from NFW */
		calc_rho_einasto (&rrr, 1, RHOS, RS, dpe->Delta, Omega_m0, &ret);
		//printf ("rho_c % 22.15g rho % 22.15g\n", rho_c / nc_hicosmo_h2 (cosmo), ncm_c_crit_mass_density_h2_solar_mass_Mpc3());
		printf ("z = % 22.15g rscale meu % 22.15g\n", z, rs);
		{
			const gdouble ret_nc = delta_c * rho_c/(x * onepx2);
			printf ("factor NC % 22.15g\n", delta_c * rho_c / (x * onepx2));
	    printf ("MEU % 22.15g DELE % 22.15g diff %22.15g\n", ret_nc, ret, ret / ret_nc - 1.0);	
		}
	}
	
	return delta_c * rho_c/(x * onepx2);
}

/* los = line of sight */
static gdouble 
_nc_density_profile_einasto_integral_density_los (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z)
{
  NcDensityProfileEinasto *dpe = NC_DENSITY_PROFILE_EINASTO (dp);

  gdouble rho_c   = _rho_crit_solar_mass_Mpc3 (cosmo, z); 
  gdouble delta_c = _nc_density_profile_einasto_deltac (dpe);
  gdouble rs      = _nc_density_profile_einasto_scale_radius (dp, cosmo, z); 
  gdouble A       = rs * delta_c * rho_c;
  gdouble x       = R / rs;
	gdouble xm1     = x - 1.0;
	gdouble xp1     = x + 1.0;
	gdouble x2m1    = xm1 * xp1;

	if (fabs (xm1) < 1.0e-6)
    return A * (1.0 / 3.0 + (- 2.0 / 5.0 + (13.0 / 35.0 - 20.0 / 63.0 * xm1) * xm1) * xm1);	
  else if (xm1 < 0.0)
    return A / x2m1 * (1.0 - 2.0 / sqrt (-x2m1) * atanh (sqrt (-xm1 / xp1)));
  else
    return A / x2m1 * (1.0 - 2.0 / sqrt (x2m1) * atan (sqrt (xm1 / xp1)));
  
}

static gdouble 
_nc_density_profile_einasto_integral_density_2d (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z)
{
  NcDensityProfileEinasto *dpe = NC_DENSITY_PROFILE_EINASTO (dp);

  gdouble rho_c   = _rho_crit_solar_mass_Mpc3 (cosmo, z); 
  gdouble delta_c = _nc_density_profile_einasto_deltac (dpe);
  gdouble rs      = _nc_density_profile_einasto_scale_radius (dp, cosmo, z); 
  gdouble A       = rs * delta_c * rho_c;
  gdouble x       = R / rs;
	gdouble xm1     = x - 1.0;
	gdouble xp1     = x + 1.0;
	gdouble x2m1    = xm1 * xp1;

  if (fabs (xm1) < 1.0e-6)
    return A * (2.0 * (1.0 - M_LN2) + (2.0 / 3.0 + (- 1.0 / 15.0 - 2.0 / 105.0 * xm1) * xm1) * xm1);	
	
  else if (xm1 < 0.0)
    return 2.0 * A * (2.0/sqrt(- x2m1) * atanh(sqrt(-xm1 / xp1)) + log(x) - M_LN2);

  else
    return 2.0 * A * (2.0/sqrt(x2m1) * atan(sqrt(xm1 / xp1)) + log(x) - M_LN2);
  
}
