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
 * As described #NcHaloDensityProfile, we just need to implement the dimensionless 3D density $\hat{\rho}(x)$
 * [which refers to the virtual function nc_halo_density_profile_eval_dl_density()].
 * In particular, the NFW profile is given by
 * \begin{equation}
 * \hat{\rho}(x) = \frac{1}{x(1 + x)^2},
 * \end{equation}
 * where $x = r/r_s$ and $r_s$ is the scale radius.
 *
 * Both the mass $M_\Delta$ and the scale profile $\rho_s$ are written in terms of the integral
 * $I_{x^2\hat\rho}(c_\Delta)$ [virtual function nc_halo_density_profile_eval_dl_spher_mass()].
 * The respective NFW implementation provides
 * \begin{equation}
 * I_{x^2\hat\rho}(x) = \ln(1 + x) - \frac{x}{1 + x}.
 * \end{equation}
 *
 * The NFW dimensionless surface mass density [virtual function nc_halo_density_profile_eval_dl_2d_density()] is
 * \begin{equation}
 * \hat{\Sigma}(X) = \frac{2}{X^2 -1} \left[1 - \frac{\arctan (\sqrt{X^2 -1})}{\sqrt{X^2 - 1}}\right].
 * \end{equation}
 * For $X^2 - 1 < 0$ the equation above can be written in terms of $\mathrm{arctanh}(\sqrt{1 - X^2})/\sqrt{1 - X^2}$.
 * If $\vert X - 1 \vert < 10^{-6}$ or $X < 10^{-6}$, $\hat{\Sigma} (X)$ is computed using
 * the Taylor series expansion at $1$ or $0$ respectively (with sufficient terms in order to obtain double precision).
 *
 * The NFW enclosed mass is [virtual function nc_halo_density_profile_eval_dl_cyl_mass()]
 * \begin{equation}
 * \hat{\overline{\Sigma}} (< X) = 2 \left[\frac{\arctan  (\sqrt{X^2 - 1})}{\sqrt{X^2 - 1}} + \ln (X / 2)\right],
 * \end{equation}
 * Similar expressions in terms of $\mathrm{arctanh}$ and approximations, as described above, are used here.
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
#include "math/ncm_integrate.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_expint.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcHaloDensityProfileNFW, nc_halo_density_profile_nfw, NC_TYPE_HALO_DENSITY_PROFILE)

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
  
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
  
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
  const gdouble Xm1              = X - 1.0;
  const gdouble Xp1              = X + 1.0;
  const gdouble abs_Xm1          = fabs (Xm1);
  const gdouble X2m1             = Xm1 * Xp1;
  const gdouble sqrt_abs_X2m1    = sqrt (fabs (X2m1));
  const gdouble sqrt_abs_Xm1_Xp1 = sqrt (abs_Xm1 / Xp1);
  
  if (abs_Xm1 < pow (GSL_DBL_EPSILON, 1.0 / 4.0))
  {
    return 2.0 * (1.0 / 3.0 + (-2.0 / 5.0 + (13.0 / 35.0 - 20.0 / 63.0 * Xm1) * Xm1) * Xm1);
  }
  else if (Xm1 < 0.0)
  {
    if (X < pow (GSL_DBL_EPSILON, 1.0 / 10.0))
    {
      const gdouble X_2     = 0.5 * X;
      const gdouble log_X_2 = log (X_2);
      
      return 2.0 / X2m1 * (
        (1.0        +  1.0 * log_X_2) +
        (1.0        +  2.0 * log_X_2) * X_2 * X_2 +
        (7.0 / 2.0  +  6.0 * log_X_2) * X_2 * X_2 * X_2 * X_2 +
        (37.0 / 3.0  + 20.0 * log_X_2) * X_2 * X_2 * X_2 * X_2 * X_2 * X_2 +
        (533.0 / 12.0 + 70.0 * log_X_2) * X_2 * X_2 * X_2 * X_2 * X_2 * X_2 * X_2 * X_2
                          );
    }
    else
    {
      return 2.0 / X2m1 * (1.0 - 2.0 * atanh (sqrt_abs_Xm1_Xp1) / sqrt_abs_X2m1);
    }
  }
  else
  {
    return 2.0 / X2m1 * (1.0 - 2.0 * atan (sqrt_abs_Xm1_Xp1) / sqrt_abs_X2m1);
  }
}

static gdouble
_nc_halo_density_profile_nfw_eval_dl_cyl_mass (NcHaloDensityProfile *dp, const gdouble X)
{
  const gdouble Xm1              = X - 1.0;
  const gdouble Xp1              = X + 1.0;
  const gdouble abs_Xm1          = fabs (Xm1);
  const gdouble sqrt_abs_X2m1    = sqrt (abs_Xm1 * Xp1);
  const gdouble sqrt_abs_Xm1_Xp1 = sqrt (abs_Xm1 / Xp1);
  
  if (abs_Xm1 < pow (GSL_DBL_EPSILON, 1.0 / 4.0))
  {
    return 2.0 * (1.0 - M_LN2) + (2.0 / 3.0 + (-1.0 / 15.0 - 2.0 / 105.0 * Xm1) * Xm1) * Xm1;
  }
  else if (Xm1 < 0.0)
  {
    if (X < pow (GSL_DBL_EPSILON, 1.0 / 12.0))
    {
      const gdouble X_2    = 0.5 * X;
      const gdouble X_22   = X_2 * X_2;
      const gdouble ln_X_2 = log (X_2);
      
      return -2.0 * (
        (
          (1.0 + 2.0 * ln_X_2) +
          (
            (7.0 / 2.0 + 6.0 * ln_X_2) +
            (
              (37.0 / 3.0  + 20.0 * ln_X_2) +
              (
                (533.0 / 12.0 + 70.0 * ln_X_2) +
                (
                  (1627.0 / 10.0 + 252.0 * ln_X_2)
                ) * X_22
              ) * X_22
            ) * X_22
          ) * X_22
        ) * X_22
                    );
    }
    else
    {
      return 2.0 * (2.0 * atanh (sqrt_abs_Xm1_Xp1) / sqrt_abs_X2m1 + log (0.5 * X));
    }
  }
  else
  {
    return 2.0 * (2.0 * atan (sqrt_abs_Xm1_Xp1) / sqrt_abs_X2m1 + log (0.5 * X));
  }
}

/**
 * nc_halo_density_profile_nfw_class_set_ni:
 * @num: boolean (true - numeric; false - analytic)  
 *
 * This function substitutes the child methods by the parent ones if num is TRUE. 
 * It enforces the computations to be performed by numerical integration.
 *
 * WARNING: this function modifies the behavior object. It should be used for testing only!   
 *
 */
void 
nc_halo_density_profile_nfw_class_set_ni (gboolean num)
{
  NcHaloDensityProfileClass *dp_class        = g_type_class_ref (NC_TYPE_HALO_DENSITY_PROFILE_NFW);
  NcHaloDensityProfileClass *dp_parent_class = g_type_class_peek_parent (dp_class);
  
  if (num)
  {
    dp_class->eval_dl_spher_mass = dp_parent_class->eval_dl_spher_mass;
    dp_class->eval_dl_2d_density = dp_parent_class->eval_dl_2d_density;
    dp_class->eval_dl_cyl_mass   = dp_parent_class->eval_dl_cyl_mass;
  }
  else
  {
    dp_class->eval_dl_spher_mass = &_nc_halo_density_profile_nfw_eval_dl_spher_mass;
    dp_class->eval_dl_2d_density = &_nc_halo_density_profile_nfw_eval_dl_2d_density;
    dp_class->eval_dl_cyl_mass   = &_nc_halo_density_profile_nfw_eval_dl_cyl_mass;
  }

  g_type_class_unref (dp_class);  
}


/**
 * nc_halo_density_profile_nfw_new:
 * @mdef: a #NcHaloDensityProfileMassDef
 * @Delta: cluster threshold mass definition $\Delta$
 *
 * This function returns the #NcHaloDensityProfileNFW implementation of
 * #NcHaloDensityProfile setting #NcHaloDensityProfile:mass-def to @mdef
 * and #NcHaloDensityProfile:Delta to @Delta.
 *
 * Returns: a new instance of #NcHaloDensityProfileNFW.
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

