/***************************************************************************
 *            nc_cluster_mass_benson_xray.c
 *
 *  Tue July 9 17:01:24 2012
 *  Copyright  2012  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_cluster_mass_benson_xray
 * @title: NcClusterMassBensonXRay
 * @short_description: Sunyaev-Zel'dovich and x-ray cluster abundance mass distribution.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass_benson_xray.h"
#include "math/integral.h"
#include "math/memory_pool.h"
#include "math/ncm_cfg.h"

#include <gsl/gsl_randist.h>

G_DEFINE_TYPE (NcClusterMassBensonXRay, nc_cluster_mass_benson_xray, NC_TYPE_CLUSTER_MASS_BENSON);

#define VECTOR (NCM_MODEL (mx)->params)
#define A_X    (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_XRAY_A_X))
#define B_X    (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_XRAY_B_X))
#define C_X    (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_XRAY_C_X))
#define D_X    (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_BENSON_XRAY_D_X))

enum
{
  PROP_0,
  PROP_YX_OBS_MIN,
  PROP_YX_OBS_MAX,
  PROP_M0X,
  PROP_Y0,
  PROP_SIZE,
};

typedef struct _integrand_data
{
  NcClusterMassBensonXRay *mx;
  NcHICosmo *model;
  gdouble *xi_params;
  gdouble z;
  gdouble lnM;
  gdouble *xi;
  gdouble lnAxh;
  gdouble lnM0;
  gdouble lnE;
  gdouble mu;
  gdouble Dx2_2;
} integrand_data;

static gdouble
_nc_cluster_mass_benson_xray_m_p (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *xi, gdouble *xi_params)
{
  gdouble Psz = NC_CLUSTER_MASS_CLASS (nc_cluster_mass_benson_xray_parent_class)->P (clusterm, model, lnM, z, xi, xi_params);
  gdouble Px = 1.0;
  if (xi[1] != 0.0)
  {
    NcClusterMassBensonXRay *mx = NC_CLUSTER_MASS_BENSON_XRAY (clusterm);
    const gdouble E = nc_hicosmo_E (model, z);
    const gdouble h3_2 = sqrt (gsl_pow_3 (nc_hicosmo_h (model)));
    const gdouble lnAxh3_2 = log (A_X * h3_2);
    const gdouble lnM0x = log (mx->M0x);
    const gdouble lnE = log (E);
    const gdouble lnYx_th_Y0 = ((lnM - lnM0x) - C_X * lnE - lnAxh3_2) / B_X;
    const gdouble lnYx_th = lnYx_th_Y0 + log (mx->Y0);
    const gdouble Yx = xi[1];
    const gdouble x = log (Yx) - lnYx_th;
    const gdouble plnYx = exp (-(x * x / (2.0 * D_X * D_X))) / (Yx * D_X * M_SQRT2 * M_SQRTPI);
    Px = plnYx;
  }

  return Psz * Px;
}

static gdouble
_nc_cluster_mass_benson_xray_intp (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z)
{
  NcClusterMassBensonXRay *mx = NC_CLUSTER_MASS_BENSON_XRAY (clusterm);
  gdouble t_min, norma;
  const gdouble E = nc_hicosmo_E (model, z);
  const gdouble h3_2 = sqrt (gsl_pow_3 (nc_hicosmo_h (model)));
  const gdouble lnAxh3_2 = log (A_X * h3_2);
  const gdouble lnM0x = log (mx->M0x);
  const gdouble lnE = log (E);
  const gdouble lnYx_th_Y0 = ((lnM - lnM0x) - C_X * lnE - lnAxh3_2) / B_X;
  const gdouble lnYx_th = lnYx_th_Y0 + log (mx->Y0);
  
  t_min = (log (mx->Yx_obs_min) - lnYx_th) / (M_SQRT2 * D_X);

  norma = erfc (t_min) / 2.0; 

  return norma * NC_CLUSTER_MASS_CLASS (nc_cluster_mass_benson_xray_parent_class)->intP (clusterm, model, lnM, z);
}

static gdouble
_Yx_to_mass (NcClusterMass *clusterm, NcHICosmo *model, gdouble z, gdouble Yx)
{
  NcClusterMassBensonXRay *mx = NC_CLUSTER_MASS_BENSON_XRAY (clusterm);
  const gdouble E = nc_hicosmo_E (model, z);
  const gdouble h3_2 = sqrt (gsl_pow_3 (nc_hicosmo_h (model)));
  
  const gdouble lnAxh3_2 = log (A_X * h3_2);
  const gdouble lnE = log (E);
  const gdouble lnM = log (mx->M0x) + lnAxh3_2 + B_X * log (Yx / mx->Y0) + C_X * lnE;
  
   //printf("z= %.10g Yx = %.10g lnM = %.10g\n", z, xi, lnM);
  
  return lnM;
}

static void
_nc_cluster_mass_benson_xray_p_limits (NcClusterMass *clusterm, NcHICosmo *model, gdouble *xi, gdouble *xi_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NC_CLUSTER_MASS_CLASS (nc_cluster_mass_benson_xray_parent_class)->P_limits (clusterm, model, xi, xi_params, lnM_lower, lnM_upper);
  if (xi[1] != 0.0)
  {
    NcClusterMassBensonXRay *mx = NC_CLUSTER_MASS_BENSON_XRAY (clusterm);
    const gdouble lnMl = GSL_MAX (_Yx_to_mass (clusterm, model, 2.0, xi[1] - 7.0 * D_X), log (NC_CLUSTER_MASS_BENSON_M_LOWER_BOUND));
    const gdouble lnMu = _Yx_to_mass (clusterm, model, 0.0, xi[1] + 7.0 * D_X);

    *lnM_lower = GSL_MIN (lnMl, *lnM_lower);
    *lnM_upper = GSL_MAX (lnMu, *lnM_upper);
  }
  return;
}

static void
_nc_cluster_mass_benson_xray_n_limits (NcClusterMass *clusterm, NcHICosmo *model, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassBensonXRay *mx = NC_CLUSTER_MASS_BENSON_XRAY (clusterm);
  const gdouble lnMl = GSL_MAX (_Yx_to_mass (clusterm, model, 2.0, mx->Yx_obs_min - 7.0 * D_X), log (NC_CLUSTER_MASS_BENSON_M_LOWER_BOUND));
  const gdouble lnMu = _Yx_to_mass (clusterm, model, 0.0, mx->Yx_obs_max + 7.0 * D_X);

  NC_CLUSTER_MASS_CLASS (nc_cluster_mass_benson_xray_parent_class)->N_limits (clusterm, model, lnM_lower, lnM_upper);
  
  *lnM_lower = GSL_MIN (lnMl, *lnM_lower);
  *lnM_upper = GSL_MAX (lnMu, *lnM_upper);
  
  return;  
}

static gboolean
_nc_cluster_mass_benson_xray_resample (NcClusterMass *clusterm, NcHICosmo *model, gdouble lnM, gdouble z, gdouble *xi, gdouble *xi_params, NcmRNG *rng)
{
  NcClusterMassBensonXRay *mx = NC_CLUSTER_MASS_BENSON_XRAY (clusterm);
  gboolean xi_obs_return;
  gdouble lnYx, lnYx_obs, mu3;
  const gdouble E = nc_hicosmo_E (model, z);
  const gdouble h3_2 = sqrt( gsl_pow_3(nc_hicosmo_h (model)));
  const gdouble lnAxh = log (A_X * h3_2);
  const gdouble lnM0 = log (mx->M0x);
  const gdouble lnE = log (E);

  mu3 = ((lnM - lnM0) - C_X * lnE - lnAxh) / B_X;
  lnYx = mu3 / log (mx->Y0);

  ncm_rng_lock (rng);
  lnYx_obs = lnYx + gsl_ran_gaussian (rng->r, D_X);
  ncm_rng_unlock (rng);
  
  xi[1] = exp (lnYx_obs);

  //printf("M = %e z = %.5g xi = %.5g xiobs = %.5g\n", exp(lnM), z, xi_mean, xi[0]);

  xi_obs_return = NC_CLUSTER_MASS_CLASS (nc_cluster_mass_benson_xray_parent_class)->resample (clusterm, model, lnM, z, xi, xi_params, rng);
  
  return (xi_obs_return && xi[1] >= mx->Yx_obs_min);
  
}

guint _nc_cluster_mass_benson_xray_obs_len (NcClusterMass *clusterm) { NCM_UNUSED (clusterm); return 2; }
guint _nc_cluster_mass_benson_xray_obs_params_len (NcClusterMass *clusterm) { NCM_UNUSED (clusterm); return 0; }

static void
_nc_cluster_mass_benson_xray_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcClusterMassBensonXRay *mx = NC_CLUSTER_MASS_BENSON_XRAY (object);
  g_return_if_fail (NC_IS_CLUSTER_MASS_BENSON_XRAY (object));

  switch (prop_id)
  {
    case PROP_YX_OBS_MIN:
      mx->Yx_obs_min = g_value_get_double (value);
      break;
    case PROP_YX_OBS_MAX:
      mx->Yx_obs_max = g_value_get_double (value);
      break;
    case PROP_M0X:
      mx->M0x = g_value_get_double (value);
      break;
    case PROP_Y0:
      mx->Y0 = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_benson_xray_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassBensonXRay *mx = NC_CLUSTER_MASS_BENSON_XRAY (object);
  g_return_if_fail (NC_IS_CLUSTER_MASS_BENSON_XRAY (object));

  switch (prop_id)
  {
    case PROP_YX_OBS_MIN:
      g_value_set_double (value, mx->Yx_obs_min);
      break;
    case PROP_YX_OBS_MAX:
      g_value_set_double (value, mx->Yx_obs_max);
      break;
    case PROP_M0X:
      g_value_set_double (value, mx->M0x);
      break;
    case PROP_Y0:
      g_value_set_double (value, mx->Y0);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cluster_mass_benson_xray_init (NcClusterMassBensonXRay *mx)
{
  mx->Yx_obs_min = 0.0;
  mx->Yx_obs_max = 0.0;
  mx->M0x = 0.0;
  mx->Y0 = 0.0;
}

static void
_nc_cluster_mass_benson_xray_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_benson_xray_parent_class)->finalize (object);
}

static void
nc_cluster_mass_benson_xray_class_init (NcClusterMassBensonXRayClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcClusterMassClass* grand_parent_class = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  grand_parent_class->P = &_nc_cluster_mass_benson_xray_m_p;
  grand_parent_class->intP = &_nc_cluster_mass_benson_xray_intp;
  grand_parent_class->P_limits = &_nc_cluster_mass_benson_xray_p_limits;
  grand_parent_class->N_limits = &_nc_cluster_mass_benson_xray_n_limits;
  grand_parent_class->resample = &_nc_cluster_mass_benson_xray_resample;
  
  grand_parent_class->obs_len = &_nc_cluster_mass_benson_xray_obs_len;
  grand_parent_class->obs_params_len = &_nc_cluster_mass_benson_xray_obs_params_len;

  object_class->finalize = _nc_cluster_mass_benson_xray_finalize;

  model_class->set_property = &_nc_cluster_mass_benson_xray_set_property;
  model_class->get_property = &_nc_cluster_mass_benson_xray_get_property;

  ncm_model_class_set_name_nick (model_class, "Benson- SZ and XRay", "Benson_SZ_XRay");
  ncm_model_class_add_params (model_class, 4, 0, PROP_SIZE);

  /**
   * NcClusterMassBensonXRay:Yx_obs_min:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_YX_OBS_MIN,
                                   g_param_spec_double ("Yx-obs-min",
                                                        NULL,
                                                        "Minimum obsevational Yx",
                                                        1.0e-20, G_MAXDOUBLE, 1.0e-1,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterMassBensonXray:Yx_obs_max:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_YX_OBS_MAX,
                                   g_param_spec_double ("Yx-obs-max",
                                                        NULL,
                                                        "Maximum obsevational Yx",
                                                        2.0, G_MAXDOUBLE, 50.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterMassBensonXRay:M0:
   *
   * Reference mass (in h^(-1) * M_sun unit) in the X-Ray proxy-mass scaling relation.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_M0X,
                                   g_param_spec_double ("M0x",
                                                        NULL,
                                                        "X Ray Reference mass",
                                                        1.0e13, G_MAXDOUBLE, 1.0e14,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterMassBensonXRay:Y0:
   *
   * Yx reference (in 10^{14} * M_sun * keV unit) in the X-Ray proxy-mass scaling relation.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_Y0,
                                   g_param_spec_double ("Y0",
                                                        NULL,
                                                        "Reference Yx",
                                                        1.0, G_MAXDOUBLE, 3.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /*
   * X-ray signal-mass scaling parameter: Ax.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_XRAY_A_X, "A_{X}", "Ax",
                               1e-8,  10.0, 1.0e-2,
                               NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_A_X,
                               NCM_PARAM_TYPE_FREE);

  /*
   * X-ray signal-mass scaling parameter: Bx.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_XRAY_B_X, "B_{X}", "Bx",
                               1e-8,  10.0, 1.0e-2,
                               NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_B_X,
                               NCM_PARAM_TYPE_FIXED);

  /*
   * X-ray signal-mass scaling parameter: Cx.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_XRAY_C_X, "C_{X}", "Cx",
                               -2.0,  8.0, 1.0e-2,
                               NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_C_X,
                               NCM_PARAM_TYPE_FIXED);
  
 /*
   * X-ray signal-mass scaling parameter: Dx.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_BENSON_XRAY_D_X, "D_{X}", "Dx",
                               1e-8,  10.0, 1.0e-2,
                               NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_BENSON_XRAY_DEFAULT_D_X,
                               NCM_PARAM_TYPE_FIXED);
  
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

}

