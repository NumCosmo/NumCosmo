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
 * The Einasto profile is defined as
 * \begin{equation}
 * \rho(r) = \rho_s \exp{\left{-\frac{2}{\alpha}\left[\left(\frac{r}{r_s}\right)^\alpha} - 1 \right]\right}},
 * \end{equation}
 * FIXME 
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
 * References: arxiv:1712.04512 FIXME
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
#define C_DELTA (ncm_vector_get (VECTOR, NC_DENSITY_PROFILE_C_DELTA))
#define M_DELTA (ncm_vector_get (VECTOR, NC_DENSITY_PROFILE_M_DELTA))
#define ALPHA (ncm_vector_get (VECTOR, NC_DENSITY_PROFILE_EINASTO_ALPHA))

enum
{
  PROP_0,
  PROP_SIZE,
};


static void
nc_density_profile_einasto_init (NcDensityProfileEinasto *dpe)
{
}

static void
_nc_density_profile_einasto_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  /* NcDensityProfileEinasto *dpe = NC_DENSITY_PROFILE_EINASTO (object); */
  g_return_if_fail (NC_IS_DENSITY_PROFILE_EINASTO (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_density_profile_einasto_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcDensityProfileEinasto *dpe = NC_DENSITY_PROFILE_EINASTO (object); */
  g_return_if_fail (NC_IS_DENSITY_PROFILE_EINASTO (object));

  switch (prop_id)
  {
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
  ncm_model_class_add_params (model_class, NC_DENSITY_PROFILE_EINASTO_SPARAM_LEN - NC_DENSITY_PROFILE_SPARAM_LEN, 0, PROP_SIZE);

  /**
   * NcDensityProfileEinasto:alpha:
   * 
   * Parameter that characterizes the profile slope.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_DENSITY_PROFILE_EINASTO_ALPHA, "\alpha", "alpha",
                              0.01,  10.0, 1.0e-1,
                              NC_DENSITY_PROFILE_EINASTO_DEFAULT_PARAMS_ABSTOL, NC_DENSITY_PROFILE_EINASTO_DEFAULT_ALPHA,
                              NCM_PARAM_TYPE_FIXED);
 
  parent_class->eval_density         = &_nc_density_profile_einasto_eval_density;
  parent_class->integral_density_los = &_nc_density_profile_einasto_integral_density_los;
  parent_class->integral_density_2d  = &_nc_density_profile_einasto_integral_density_2d;
  parent_class->eval_fourier         = &_nc_density_profile_einasto_eval_fourier;
  
}

/**
 * nc_density_profile_einasto_new:
 * @mdef: a #NcDensityProfileMassDef
 * @Delta: cluster threshold mass definition $\Delta$ 
 *
 * This function returns a #NcDensityProfile with a #NcDensityProfileEinasto implementation.
 *
 * Returns: A new #NcDensityProfileEinasto.
 */
NcDensityProfile *
nc_density_profile_einasto_new (const NcDensityProfileMassDef mdef, const gdouble Delta)
{
  return g_object_new (NC_TYPE_DENSITY_PROFILE_EINASTO, 
                       "mass-def", mdef,
                       "Delta",    Delta,
                       NULL);
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

static gdouble
_nc_density_profile_einasto_eval_fourier (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble k, const gdouble M, const gdouble z)
{
  /* FIXME */
	return 0.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////


static gdouble
_nc_density_profile_einasto_eval_density (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble r, const gdouble z)
{
  NcDensityProfileEinasto *dpe = NC_DENSITY_PROFILE_EINASTO (dp);
  const gdouble rs      = nc_density_profile_scale_radius (dp, cosmo, z);
  const gdouble x       = r / rs;
  const gdouble delta_c = nc_density_profile_deltac (dp, cosmo, z);
  const gdouble rho     = nc_density_profile_mass_density (dp, cosmo, z); 

	return delta_c * rho * exp (- (2.0 / ALPHA) * (pow(x, ALPHA) - 1.0));
}

/* los = line of sight */
static gdouble 
_nc_density_profile_einasto_integral_density_los (NcDensityProfile *dp, NcHICosmo *cosmo, const gdouble R, const gdouble z)
{
  //NcDensityProfileEinasto *dpe = NC_DENSITY_PROFILE_EINASTO (dp);

	const gdouble rho     = nc_density_profile_mass_density (dp, cosmo, z); 
  const gdouble delta_c = nc_density_profile_deltac (dp, cosmo, z);
  const gdouble rs      = nc_density_profile_scale_radius (dp, cosmo, z); 
  gdouble A       = rs * delta_c * rho;
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
  //NcDensityProfileEinasto *dpe = NC_DENSITY_PROFILE_EINASTO (dp);

	const gdouble rho     = nc_density_profile_mass_density (dp, cosmo, z); 
  const gdouble delta_c = nc_density_profile_deltac (dp, cosmo, z);
  const gdouble rs      = nc_density_profile_scale_radius (dp, cosmo, z);
  gdouble A       = rs * delta_c * rho;
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
