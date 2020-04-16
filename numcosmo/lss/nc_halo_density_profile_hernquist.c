/***************************************************************************
 *            nc_halo_density_profile_hernquist.c
 *
 *  Sat April 11 17:23:06 2020
 *  Copyright  2014
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_halo_density_profile_hernquist.c
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_halo_density_profile_hernquist
 * @title: NcHaloDensityProfileHernquist
 * @short_description: Density profile of Hernquist type.
 *
 * This object implements the #NcHaloDensityProfile class for a Hernquist density profile.
 * 
 * As described #NcHaloDensityProfile, we just need to implement the dimensionless 3D density $\hat{\rho}(x)$ 
 * [which refers to the virtual function nc_halo_density_profile_eval_dl_density()]. 
 * In particular, the Hernquist profile is given by
 * \begin{equation}
 * \hat{\rho}(x) = \frac{1}{x(1 + x)^3},
 * \end{equation}
 * where $x = r/r_s$ and $r_s$ is the scale radius.
 *
 * Both the mass $M_\Delta$ and the scale profile $\rho_s$ are written in terms of the integral 
 * $I_{x^2\hat\rho}(c_\Delta)$ [virtual function nc_halo_density_profile_eval_dl_spher_mass()]. 
 * The respective Hernquist implementation provides  
 * \begin{equation}
 * I_{x^2\hat\rho}(x) = \frac{x^2}{2(1 + x)^2}.
 * \end{equation}
 *
 * The Hernquist dimensionless surface mass density [virtual function nc_halo_density_profile_eval_dl_2d_density()] is 
 * \begin{equation}
 * \hat{\Sigma}(X) = \frac{1}{\left( X^2 - 1 \right)^2} \left[\frac{\left(2 + X^2 \right) \arctan \left(\sqrt{X^2 - 1} \right)}{\sqrt{\left( X^2 - 1 \right)}} - 3\right].
 * \end{equation} 
 * For $X^2 - 1 < 0$ the equation above can be written in terms of $\mathrm{arctanh}(\sqrt{\vert X^2 - 1 \vert})$. 
 * If $\vert X - 1 \vert < 10^{-6}$ or $X < 10^{-6}$, $\hat{\Sigma} (X)$ is computed using 
 * the Taylor series expansion at $1$ or $0$ respectively (with sufficient terms in order to obtain double precision).
 *  
 * The Hernquist enclosed mass is [virtual function nc_halo_density_profile_eval_dl_cyl_mass()]
 * \begin{equation}
 * \hat{\overline{\Sigma}} (< X) = X^2 / (X^2 - 1) * \eft[1 - 2 * \arctan \left(\sqrt{\vert X - 1 \vert / (X + 1)\right)\right];,
 * \end{equation}
 * Similar expressions in terms of $\mathrm{arctanh}$ and approximations, as described above, are used here.
 * 
 * References: , arxiv:1712.04512.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_halo_density_profile_hernquist.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/integral.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_expint.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcHaloDensityProfileHernquist, nc_halo_density_profile_hernquist, NC_TYPE_HALO_DENSITY_PROFILE);

#define VECTOR  (NCM_MODEL (dph)->params)
#define M_DELTA (ncm_vector_get (VECTOR, NC_HALO_DENSITY_PROFILE_M_DELTA))
#define C_DELTA (ncm_vector_get (VECTOR, NC_HALO_DENSITY_PROFILE_C_DELTA))

enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_halo_density_profile_hernquist_init (NcHaloDensityProfileHernquist *dph)
{
}

static void
_nc_halo_density_profile_hernquist_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /* NcHaloDensityProfileHernquist *dph = NC_HALO_DENSITY_PROFILE_HERNQUIST (object); */
  g_return_if_fail (NC_IS_HALO_DENSITY_PROFILE_HERNQUIST (object));
  
  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_density_profile_hernquist_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcHaloDensityProfileHernquist *dph = NC_HALO_DENSITY_PROFILE_HERNQUIST (object); */
  g_return_if_fail (NC_IS_HALO_DENSITY_PROFILE_HERNQUIST (object));
  
  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_density_profile_hernquist_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_density_profile_hernquist_parent_class)->finalize (object);
}

static gdouble _nc_halo_density_profile_hernquist_eval_dl_density (NcHaloDensityProfile *dp, const gdouble x);
static gdouble _nc_halo_density_profile_hernquist_eval_dl_spher_mass (NcHaloDensityProfile *dp, const gdouble x);
static gdouble _nc_halo_density_profile_hernquist_eval_dl_2d_density (NcHaloDensityProfile *dp, const gdouble X);
static gdouble _nc_halo_density_profile_hernquist_eval_dl_cyl_mass (NcHaloDensityProfile *dp, const gdouble X);

static void
nc_halo_density_profile_hernquist_class_init (NcHaloDensityProfileHernquistClass *klass)
{
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);
  NcHaloDensityProfileClass *dp_class = NC_HALO_DENSITY_PROFILE_CLASS (klass);
  NcmModelClass *model_class          = NCM_MODEL_CLASS (klass);
  
  model_class->set_property = &_nc_halo_density_profile_hernquist_set_property;
  model_class->get_property = &_nc_halo_density_profile_hernquist_get_property;
  object_class->finalize    = &_nc_halo_density_profile_hernquist_finalize;
  
  ncm_model_class_set_name_nick (model_class, "Hernquist Density Profile", "Hernquist");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
  
  dp_class->eval_dl_density    = &_nc_halo_density_profile_hernquist_eval_dl_density;
  dp_class->eval_dl_spher_mass = &_nc_halo_density_profile_hernquist_eval_dl_spher_mass;
  dp_class->eval_dl_2d_density = &_nc_halo_density_profile_hernquist_eval_dl_2d_density;
  dp_class->eval_dl_cyl_mass   = &_nc_halo_density_profile_hernquist_eval_dl_cyl_mass;
}

static gdouble 
_nc_halo_density_profile_hernquist_eval_dl_density (NcHaloDensityProfile *dp, const gdouble x)
{
  return 1.0 / (x * gsl_pow_3 (1.0 + x));
}

static gdouble 
_nc_halo_density_profile_hernquist_eval_dl_spher_mass (NcHaloDensityProfile *dp, const gdouble x)
{
  return gsl_pow_2 (x) / ( 2.0 * gsl_pow_2 (1.0 + x));
}

static gdouble 
_nc_halo_density_profile_hernquist_eval_dl_2d_density (NcHaloDensityProfile *dp, const gdouble X)
{
  const gdouble Xm1              = X - 1.0;
  const gdouble Xp1              = X + 1.0;
  const gdouble abs_Xm1          = fabs (Xm1);
  const gdouble X2m1             = Xm1 * Xp1;
  const gdouble X2m1_2           = X2m1 * X2m1;
  const gdouble sqrt_abs_X2m1    = sqrt (fabs (X2m1));
  const gdouble sqrt_abs_Xm1_Xp1 = sqrt (abs_Xm1 / Xp1);

  if (abs_Xm1 < pow (GSL_DBL_EPSILON, 1.0 / 4.0))
  {
    return 4.0 * (1.0 / 15.0 + (-4.0 / 35.0 + (2.0 / 15.0 - 92.0 / 693.0 * Xm1) * Xm1) * Xm1);
  }
  else if (X2m1 < 0.0)
  {
    return 1.0 / X2m1_2 * (- 3.0 + (2.0 + X * X) * 2.0 * atanh (sqrt_abs_Xm1_Xp1) / sqrt_abs_X2m1);
  }
  else
  {
    return 1.0 / X2m1_2 * (- 3.0 + (2.0 + X * X) * 2.0 * atan (sqrt_abs_Xm1_Xp1) / sqrt_abs_X2m1);
  }
}

static gdouble 
_nc_halo_density_profile_hernquist_eval_dl_cyl_mass (NcHaloDensityProfile *dp, const gdouble X)
{
  const gdouble Xm1              = X - 1.0;
  const gdouble Xp1              = X + 1.0;
  const gdouble X2m1             = Xm1 * Xp1;
  const gdouble abs_Xm1          = fabs (Xm1);
  const gdouble sqrt_abs_X2m1    = sqrt (fabs (X2m1));
  const gdouble sqrt_abs_Xm1_Xp1 = sqrt (abs_Xm1 / Xp1);

  if (abs_Xm1 < pow (GSL_DBL_EPSILON, 1.0 / 4.0))
  {
    return 1.0 / 3.0 + (4.0 / 15.0 + (- 2.0 / 21.0 + 8.0 / 315.0 * Xm1) * Xm1) * Xm1;
  }
  else if (Xm1 < 0.0)
  {
    if (X < pow (GSL_DBL_EPSILON, 1.0 / 10.0))
    {
      const gdouble X_2    = 0.5 * X;
      const gdouble X_22   = X_2 * X_2;
      const gdouble ln_X_2 = log (X_2);
      
      return -4.0 * (
                     (
                      (1.0 + ln_X_2) + 
                      (
                       (5.0 + 6.0 * ln_X_2) + 
                       (
                        (47.0 / 2.0  + 30.0 * ln_X_2) +
                        (
                         (319.0 / 3.0 + 140.0 * ln_X_2)
                        ) * X_22
                       ) * X_22
                      ) * X_22
                     ) * X_22
                    );
    }
    else
    {
      return X * X / X2m1 * (1.0 - 2.0 * atanh (sqrt_abs_Xm1_Xp1) / sqrt_abs_X2m1);
    }
  }
  else
  {
    return X * X / X2m1 * (1.0 - 2.0 * atan (sqrt_abs_Xm1_Xp1) / sqrt_abs_X2m1);
  }
}

/**
 * nc_halo_density_profile_hernquist_new:
 * @mdef: a #NcHaloDensityProfileMassDef
 * @Delta: cluster threshold mass definition $\Delta$
 *
 * This function returns the #NcHaloDensityProfileHernquist implementation of 
 * #NcHaloDensityProfile setting #NcHaloDensityProfile:mass-def to @mdef
 * and #NcHaloDensityProfile:Delta to @Delta.
 *
 * Returns: a new instance of #NcHaloDensityProfileHernquist.
 */
NcHaloDensityProfileHernquist *
nc_halo_density_profile_hernquist_new (const NcHaloDensityProfileMassDef mdef, const gdouble Delta)
{
  NcHaloDensityProfileHernquist *dph = g_object_new (NC_TYPE_HALO_DENSITY_PROFILE_HERNQUIST,
                                                  "mass-def", mdef,
                                                  "Delta",    Delta,
                                                  NULL);

  return dph;
}
