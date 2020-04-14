/***************************************************************************
 *            nc_cluster_mass_plcl.c
 *
 *  Sun Mar 1 22:00:23 2015
 *  Copyright  2015  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2015 <pennalima@gmail.com>
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
 * SECTION:nc_cluster_mass_plcl
 * @title: NcClusterMassPlCL
 * @short_description: Planck-CLASH Cluster Mass Distribution
 *
 * FIXME Planck-CLASH Cluster Mass Distribution (SZ - Lensing).
 *
 * Do not use this object to perform cluster abundance analyses. For now, it is suitable just for cluster pseudo counts.  
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass_plcl.h"
#include "math/integral.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcClusterMassPlCL, nc_cluster_mass_plcl, NC_TYPE_CLUSTER_MASS);

#define VECTOR (NCM_MODEL (mszl)->params)
#define A_SZ   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_PLCL_A_SZ))
#define B_SZ   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_PLCL_B_SZ))
#define SD_SZ  (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_PLCL_SD_SZ))
#define A_L    (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_PLCL_A_L))
#define B_L    (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_PLCL_B_L))
#define SD_L   (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_PLCL_SD_L))
#define COR    (ncm_vector_get (VECTOR, NC_CLUSTER_MASS_PLCL_COR)) 

enum
{
  PROP_0,
  PROP_M0,
  PROP_SIZE,
};

static void
nc_cluster_mass_plcl_init (NcClusterMassPlCL *mszl)
{
  mszl->M0 = 0.0;
  mszl->T = gsl_multifit_fdfsolver_lmsder;
  mszl->s = gsl_multifit_fdfsolver_alloc (mszl->T, 4, 2);
}

static void
_nc_cluster_mass_plcl_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (object);
  g_return_if_fail (NC_IS_CLUSTER_MASS_PLCL (object));

  switch (prop_id)
  {
    case PROP_M0:
      mszl->M0 = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_plcl_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (object);
  g_return_if_fail (NC_IS_CLUSTER_MASS_PLCL (object));

  switch (prop_id)
  {
    case PROP_M0:
      g_value_set_double (value, mszl->M0);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_mass_plcl_finalize (GObject *object)
{
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (object);
  gsl_multifit_fdfsolver_free (mszl->s);
    
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_plcl_parent_class)->finalize (object);
}

guint _nc_cluster_mass_plcl_obs_len (NcClusterMass *clusterm) { NCM_UNUSED (clusterm); return 2; }
guint _nc_cluster_mass_plcl_obs_params_len (NcClusterMass *clusterm) { NCM_UNUSED (clusterm); return 2; }
static gdouble _nc_cluster_mass_plcl_Msz_Ml_M500_p (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *Mobs, const gdouble *Mobs_params);
static gdouble _nc_cluster_mass_plcl_intp (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z);
static gboolean _nc_cluster_mass_plcl_resample (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnMobs, const gdouble *lnMobs_params, NcmRNG *rng);
static void _nc_cluster_mass_plcl_n_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper);

static void
nc_cluster_mass_plcl_class_init (NcClusterMassPlCLClass *klass)
{
  GObjectClass* object_class       = G_OBJECT_CLASS (klass);
  NcClusterMassClass* parent_class = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_cluster_mass_plcl_set_property;
  model_class->get_property = &_nc_cluster_mass_plcl_get_property;
  object_class->finalize    = &_nc_cluster_mass_plcl_finalize;

  ncm_model_class_set_name_nick (model_class, "Planck - CLASH cluster mass distribution", "Planck_CLASH");
  ncm_model_class_add_params (model_class, NC_CLUSTER_MASS_PLCL_SPARAM_LEN, 0, PROP_SIZE);

  /**
   * NcClusterMassPlCL:M0:
   *
   * Reference mass (in h^(-1) * M_sun unit) in the SZ signal-mass scaling relation.
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_M0,
                                   g_param_spec_double ("M0",
                                                        NULL,
                                                        "Reference mass",
                                                        1.0e13, G_MAXDOUBLE, 5.7e14,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterMassPlCL:Asz:
   * 
   * SZ signal-mass scaling parameter.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PLCL_A_SZ, "\\alpha_{SZ}", "Asz",
                              0.5, 1.5, 0.1,
                              NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PLCL_DEFAULT_A_SZ,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassPlCL:Bsz:
   * 
   * SZ signal-mass scaling parameter.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PLCL_B_SZ, "b_{SZ}", "Bsz",
                              -1.0,  0.9999, 2.0e-2,
                              NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PLCL_DEFAULT_B_SZ,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassPlCL:sigma_sz:
   * 
   * Standard deviation of the SZ signal-mass scaling relation.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PLCL_SD_SZ, "\\sigma_{SZ}", "sigma_sz",
                              1e-2,  1.0, 5.0e-2,
                              NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PLCL_DEFAULT_SD_SZ,
                              NCM_PARAM_TYPE_FIXED);
  /**
   * NcClusterMassPlCL:Al:
   * 
   * Lensing signal-mass scaling parameter.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PLCL_A_L, "\\alpha_{L}", "Al",
                              0.5,  1.5, 1.0e-1,
                              NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PLCL_DEFAULT_A_L,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassPlCL:Bl:
   * 
   * Lensing signal-mass scaling parameter.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PLCL_B_L, "b_{L}", "Bl",
                              -1.0,  0.9999, 2.0e-2,
                              NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PLCL_DEFAULT_B_L,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcClusterMassPlCL:sigma_l:
   * 
   * Standard deviation of the lensing signal-mass scaling relation.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PLCL_SD_L, "\\sigma_{L}", "sigma_l",
                              1e-2,  1.0, 5.0e-2,
                              NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PLCL_DEFAULT_SD_L,
                              NCM_PARAM_TYPE_FIXED);
  /**
   * NcClusterMassPlCL:cor:
   * 
   * SZ-Lensing signal-mass correlation, $0.0 \leq \rho \leq 1.0$.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PLCL_COR, "\\rho", "cor",
                              -0.99999,  0.99999, 2.0e-1,
                              NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PLCL_DEFAULT_COR,
                              NCM_PARAM_TYPE_FIXED);
    
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);  

  parent_class->P              = &_nc_cluster_mass_plcl_Msz_Ml_M500_p;
  parent_class->intP           = &_nc_cluster_mass_plcl_intp;
  //parent_class->P_limits       = &_nc_cluster_mass_plcl_p_limits;
  parent_class->N_limits       = &_nc_cluster_mass_plcl_n_limits;
  parent_class->resample       = &_nc_cluster_mass_plcl_resample;
  parent_class->obs_len        = &_nc_cluster_mass_plcl_obs_len;
  parent_class->obs_params_len = &_nc_cluster_mass_plcl_obs_params_len;

  ncm_model_class_add_impl_flag (model_class, NC_CLUSTER_MASS_IMPL_ALL);
}

typedef struct _integrand_data
{
  NcClusterMassPlCL *mszl;
  NcHICosmo *cosmo;
  const gdouble *mobs_params;
  const gdouble *mobs; /*observables: Msz, Ml*/
  gdouble lnM;
  gdouble mu_sz; /* SZ mean mass given the MO relation */
  gdouble mu_l; /* Lensing mean mass given the MO relation */
  gdouble lnnorma_p; /* normalization of the combined density prob. distributions of P(M_Pl|Msz), P(M_CL|Ml), P(lnMsz, lnMl|lnM500)*/
  gdouble sd_sz2;
  gdouble sd_l2;
  gdouble twocor_sdsz_sdl;
  gdouble Mcut;
  gdouble peak[2];
  gdouble func_peak; /* function computed at the peak*/
  gdouble c1;
  gdouble c2;
  gdouble erf_Mlow;
  gdouble erf_Mup;
  gdouble erf_const_Msz;
  gdouble erf_const_Ml;
} integrand_data;

static gdouble
_SZ_lnmass_mean (NcClusterMassPlCL *mszl,  gdouble lnM)
{
  const gdouble lnM0 = log (mszl->M0);
  const gdouble lnMsz_M0_mean = log1p (- B_SZ) + A_SZ * (lnM - lnM0);
  return lnMsz_M0_mean;
}

static gdouble
_Lens_lnmass_mean (NcClusterMassPlCL *mszl, gdouble lnM)
{
  const gdouble lnM0 = log (mszl->M0);
  const gdouble lnMlens_M0_mean = log1p (- B_L) + A_L * (lnM - lnM0);
  return lnMlens_M0_mean;
}

/**
 * nc_cluster_mass_plcl_pdf_only_lognormal:
 * @clusterm: a #NcClusterMass
 * @lnM:logarithm base e of the true mass
 * @lnMsz_M0: logarithm base e of the SZ mass minus lnM0 (pivot mass)
 * @lnMl_M0: logarithm base e of the lensing mass minus lnM0 (pivot mass)
 *
 *
 * Returns: $P(\ln M_{SZ}, \ln M_L | M_{500})$
*/
gdouble 
nc_cluster_mass_plcl_pdf_only_lognormal (NcClusterMass *clusterm, gdouble lnM, gdouble lnMsz_M0, gdouble lnMl_M0)
{
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (clusterm);

  const gdouble cor2     = COR * COR;
  const gdouble norma    = 2.0 * M_PI * SD_SZ * SD_L * sqrt (1 - cor2);
  const gdouble Msz_mean = _SZ_lnmass_mean (mszl, lnM);
  const gdouble Ml_mean  = _Lens_lnmass_mean (mszl, lnM);
  const gdouble ysz      = (lnMsz_M0 - Msz_mean) / SD_SZ;
  const gdouble arg_ysz  = ysz * ysz;
  const gdouble yl       = (lnMl_M0 - Ml_mean) / SD_L;
  const gdouble arg_yl   = yl * yl;
  const gdouble arg_sz_l = 2.0 * COR * ysz * yl;
  
  const gdouble exp_arg = - (arg_ysz + arg_yl - arg_sz_l) / (2.0 * (1.0 - cor2));

  if (exp_arg < GSL_LOG_DBL_MIN)
    return exp (-200.0);
  else
  {
    const gdouble result = exp (exp_arg) + exp (-200.0);
    //printf ("Msz = %.8g Ml = %.8g m_dist = %.8g\n", Msz, Ml, result);
    //printf ("M0 = %.2e lnM_M0 = %.5g\n", mszl->M0, lnM_M0);
    //printf ("Mpl = %.8g sd_pl = %.8g Mcl = %.8g sd_cl = %.8g\n", M_Pl, sd_Pl, M_CL, sd_CL);
    //printf ("mean1 = %.8g mean2 = %.8g\n", M1_mean/(M_SQRT2 * sd_Pl), M2_mean/(M_SQRT2 * sd_CL));
    //printf ("norma_partial = %.10g\n", norma_partial);
    return result / norma;
  }
}

/**
 * nc_cluster_mass_plcl_pdf:
 * @clusterm: a #NcClusterMass
 * @lnM_M0:logarithm base e of the true mass minus lnM0 (pivot mass)
 * @w1: new variable 1 
 * @w2: new variable 2
 * @Mobs: (array) (element-type double): observed masses
 * @Mobs_params: (array) (element-type double): observed mass paramaters
 *
 * FIXME
 * Integrals in $M_{sz}$ and $M_l$ performed in the dimensionless quantities $\ln (M_{sz} / M_0)$ 
 * and $\ln (M_l / M_0)$, respectively. The Gaussian distributions between $M_{Pl}$ and $M_{CL}$ 
 * are written in terms of the dimensionless quantities $M_{Pl}/M_0$, $M_{CL}/M_0$, $\sigma_{Pl}/M_0$ 
 * and $\sigma_{CL}/M_0$. 
 * 
 * This distribution is "partially" normalized. The constant normalization factor is included 
 * only in nc_cluster_pseudo_counts_posterior_numerator_plcl().
 *
 * Returns: FIXME
*/
gdouble 
nc_cluster_mass_plcl_pdf (NcClusterMass *clusterm, gdouble lnM_M0, gdouble w1, gdouble w2, const gdouble *Mobs, const gdouble *Mobs_params)
{
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (clusterm);
  
  /* These four variables are dimensionless */
  const gdouble M_Pl  = Mobs[NC_CLUSTER_MASS_PLCL_MPL];
  const gdouble sd_Pl = Mobs_params[NC_CLUSTER_MASS_PLCL_SD_PL];
  const gdouble M_CL  = Mobs[NC_CLUSTER_MASS_PLCL_MCL];
  const gdouble sd_CL = Mobs_params[NC_CLUSTER_MASS_PLCL_SD_CL];

  const gdouble M1_mean   = (1.0 - B_SZ) * exp (A_SZ * lnM_M0 + (sqrt (1.0 - COR * COR) * w1 + COR * w2) * SD_SZ);
  const gdouble M2_mean   = (1.0 - B_L) * exp (A_L * lnM_M0 + w2 * SD_L);
  const gdouble ysz       = (M_Pl - M1_mean) / sd_Pl;
  const gdouble arg_ysz   = ysz * ysz / 2.0;
  const gdouble yl        = (M_CL - M2_mean) / sd_CL;
  const gdouble arg_yl    = yl * yl / 2.0;
  const gdouble arg_gauss = (w1 * w1 + w2 * w2) / 2.0;
  
  const gdouble exp_arg = - arg_ysz - arg_yl - arg_gauss;

  if (exp_arg < GSL_LOG_DBL_MIN)
    return exp (-200.0);
  else
  {
    const gdouble norma_partial = (1.0 + erf(M1_mean/(M_SQRT2 * sd_Pl))) * (1.0 + erf(M2_mean/(M_SQRT2 * sd_CL)));   
    const gdouble result = exp (exp_arg) + exp (-200.0);
    //printf ("Msz = %.8g Ml = %.8g m_dist = %.8g\n", Msz, Ml, result);
    //printf ("M0 = %.2e lnM_M0 = %.5g\n", mszl->M0, lnM_M0);
    //printf ("Mpl = %.8g sd_pl = %.8g Mcl = %.8g sd_cl = %.8g\n", M_Pl, sd_Pl, M_CL, sd_CL);
    //printf ("mean1 = %.8g mean2 = %.8g\n", M1_mean/(M_SQRT2 * sd_Pl), M2_mean/(M_SQRT2 * sd_CL));
    //printf ("norma_partial = %.10g\n", norma_partial);
    return result / norma_partial;
  }
}

static gdouble
_nc_cluster_mass_plcl_Msz_Ml_M500_p_integrand (gdouble lnMsz, gdouble lnMl, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterMassPlCL *mszl = data->mszl;
  const gdouble Msz = exp (lnMsz) * mszl->M0;
  const gdouble Ml  = exp (lnMl) * mszl->M0;
  const gdouble diff_Msz = lnMsz - data->mu_sz;
  const gdouble diff_Ml  = lnMl - data->mu_l;
  
  const gdouble M_Pl  = data->mobs[NC_CLUSTER_MASS_PLCL_MPL];
  const gdouble sd_Pl = data->mobs_params[NC_CLUSTER_MASS_PLCL_SD_PL];
  const gdouble M_CL  = data->mobs[NC_CLUSTER_MASS_PLCL_MCL];
  const gdouble sd_CL = data->mobs_params[NC_CLUSTER_MASS_PLCL_SD_CL];
  
  const gdouble ysz     = (M_Pl - Msz) / sd_Pl;
  const gdouble arg_ysz = ysz * ysz / 2.0;
  const gdouble yl      = (M_CL - Ml) / sd_CL;
  const gdouble arg_yl  = yl * yl / 2.0;
  
  const gdouble xsz       = diff_Msz / SD_SZ;
  const gdouble arg_xsz   = xsz * xsz;
  const gdouble xl        = diff_Ml / SD_L;
  const gdouble arg_xl    = xl * xl;
  const gdouble arg_x_szl = data->twocor_sdsz_sdl * diff_Msz * diff_Ml;

  const gdouble exp_arg = - arg_ysz - arg_yl - (arg_xsz + arg_xl - arg_x_szl) / (2.0 * (1.0 - COR * COR));
  //const gdouble exp_arg = - arg_ysz - arg_yl;
  //const gdouble exp_arg = - (arg_xsz + arg_xl - arg_x_szl) / (2.0 * (1.0 - COR * COR));
  //printf ("%20.15g %20.15g %20.15g %20.15g\n", lnMsz, lnMl, exp_arg, exp (exp_arg - data->lnnorma_p));
  //printf ("%20.15g %20.15g %20.15g %20.15g %20.15g\n", - arg_ysz, - arg_yl, - arg_xsz / (2.0 * (1.0 - COR * COR)),
  //        - arg_xl / (2.0 * (1.0 - COR * COR)), + arg_x_szl / (2.0 * (1.0 - COR * COR)));
  //printf ("Massas: %8.5g %8.5g %8.5e %8.5e %8.5g %8.5g %8.5e %8.5e\n", lnMsz, data->mu_sz, M_Pl, Msz, lnMl, data->mu_l, M_CL, Ml);
  if (exp_arg < GSL_LOG_DBL_MIN)
    return exp (-200.0);
  else
  {
    const gdouble result = exp (exp_arg - data->lnnorma_p) + exp (-200.0);
    //printf ("===> %12.8g %12.8g %12.8g %12.8g\n", data->lnM, lnMsz, lnMl, result);
    return result;
  }
}

/**
 * nc_cluster_mass_plcl_peak_new_variables:
 * @N: FIXME
 * @lb: lower bounds
 * @ub: upper bounds 
 * @mszl: a #NcClusterMassPlCL
 * @lnM: logarithm base e of the mass
 * @Mobs: (array) (element-type double): observed mass
 * @Mobs_params: (array) (element-type double): observed mass paramaters
 *
 * FIXME
 *
*/
void
nc_cluster_mass_plcl_peak_new_variables (gdouble N, gdouble *lb, gdouble *ub, NcClusterMassPlCL *mszl, gdouble lnM, const gdouble *Mobs, const gdouble *Mobs_params)
{
  integrand_data data;
  const gdouble onemcor2  = sqrt (1.0 - COR * COR);
  gdouble lnM500, M_PL, M_CL, sd_PL, sd_CL, w1_m, w1_p, w2_m, w2_p;
  
  data.mszl =            mszl;
  data.lnM =             lnM;
  data.mobs =            Mobs;
  data.mobs_params =     Mobs_params;
   
  lnM500 = lnM - log (mszl->M0);
  M_PL = data.mobs[NC_CLUSTER_MASS_PLCL_MPL] / mszl->M0;
  M_CL = data.mobs[NC_CLUSTER_MASS_PLCL_MCL] / mszl->M0;
  sd_PL = data.mobs_params[NC_CLUSTER_MASS_PLCL_SD_PL] / mszl->M0;
  sd_CL = data.mobs_params[NC_CLUSTER_MASS_PLCL_SD_CL] / mszl->M0;
  
  if (N == 0.0)
  {
    w1_m = (log (M_PL / (1.0 - B_SZ)) / SD_SZ + log ((1.0 - B_L) / M_CL) * COR / SD_L 
        + (A_L * COR / SD_L - A_SZ / SD_SZ) * lnM500) / onemcor2; 
    w2_m = (log (M_CL / (1.0 - B_L)) - A_L * lnM500) / SD_L;
    
    //printf ("M_PL = %.5g M_CL = %.5g ASZ = %.5g BSZ = %.5g AL = %.5g BL = %.5g COR = %.5g\n", M_PL, M_CL, A_SZ, B_SZ, A_L, B_L, COR);
    //printf ("lnM500 = %.5g SDSZ = %.5g SDL = %.5g\n", lnM500, SD_SZ, SD_L);
    //printf ("w1 = %.5g w2 = %.5g\n", w1_m, w2_m);
    
    lb[0] = w1_m;
    lb[1] = w2_m;
    ub[0] = lb[0];
    ub[1] = lb[1];
    //printf ("w1 = %.5g w2 = %.5g\n", lb[0], lb[1]);
  }
  else 
  {
    gdouble w1_min, w1_max, w2_min, w2_max;
    gdouble dif_Mpl = M_PL - N * sd_PL;
    gdouble dif_Mcl = M_CL - N * sd_CL;
    if (dif_Mpl < 0.0)
      dif_Mpl = 1.0e-2; /* 1.0e12 / M0 */
    if (dif_Mcl < 0.0)
      dif_Mpl = 1.0e-2;
    
    w1_m = (log ( dif_Mpl / (1.0 - B_SZ)) / SD_SZ + log ((1.0 - B_L) / M_CL) * COR / SD_L 
        + (A_L * COR / SD_L - A_SZ / SD_SZ) * lnM500) / onemcor2;
    w1_p = (log ((M_PL + N * sd_PL) / (1.0 - B_SZ)) / SD_SZ + log ((1.0 - B_L) / M_CL) * COR / SD_L 
        + (A_L * COR / SD_L - A_SZ / SD_SZ) * lnM500) / onemcor2;
    w2_m = (log (dif_Mcl / (1.0 - B_L)) - A_L * lnM500) / SD_L;
    w2_p = (log ((M_CL + N * sd_CL) / (1.0 - B_L)) - A_L * lnM500) / SD_L;

    w1_min = GSL_MIN (w1_m, w1_p);
    w1_max = GSL_MAX (w1_m, w1_p);
    w2_min = GSL_MIN (w2_m, w2_p);
    w2_max = GSL_MAX (w2_m, w2_p);
        
    lb[0] = w1_min;
    lb[1] = w2_min;
    ub[0] = w1_max;
    ub[1] = w2_max;
  }
}

/**
 * nc_cluster_mass_plcl_gsl_f_new_variables: (skip)
 * @p: FIXME
 * @hx: FIXME
 * @mszl: a #NcClusterMassPlCL
 * @lnM_M0: logarithm base e of the mass divided by the pivot mass
 * @Mobs: (array) (element-type double): observed mass
 * @Mobs_params: (array) (element-type double): observed mass paramaters
 *
 * FIXME
 *
*/
void
nc_cluster_mass_plcl_gsl_f_new_variables (const gsl_vector *p, gsl_vector *hx, NcClusterMassPlCL *mszl, gdouble lnM_M0, const gdouble *Mobs, const gdouble *Mobs_params)
{
  const gdouble onemcor2  = sqrt (1.0 - COR * COR);
  const gdouble w1 = gsl_vector_get (p, 0); //p[0];
  const gdouble w2 = gsl_vector_get (p, 1); //p[1];
  gdouble M_PL, M_CL, sd_PL, sd_CL, dw1, dw2;

  /* Both masses and errors are given in units of the pivot mass, i.e., M_PL -> M_PL / M0 */
  M_PL  = Mobs[NC_CLUSTER_MASS_PLCL_MPL]; 
  M_CL  = Mobs[NC_CLUSTER_MASS_PLCL_MCL];
  sd_PL = Mobs_params[NC_CLUSTER_MASS_PLCL_SD_PL];
  sd_CL = Mobs_params[NC_CLUSTER_MASS_PLCL_SD_CL];

  dw1 = (M_PL - (1.0 - B_SZ) * exp (A_SZ * lnM_M0 + (onemcor2 * w1 + COR * w2) * SD_SZ)) / sd_PL;
  dw2 = (M_CL - (1.0 - B_L) * exp (A_L * lnM_M0 + w2 * SD_L)) / sd_CL;
  
  //printf ("p[0] = %.5g p[1] = %.5g\n", p[0], p[1]); 
  //printf ("Mpl = %.5g Mcl = %.5g sd_pl = %.5g sd_cl = %.5g\n", M_PL, M_CL, sd_PL, sd_CL);
  //printf ("lnMsz/M0 = %.5g lnMl/M0 = %.5g\n", lnMsz_M0, lnMl_M0);
  //printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", w1, w2, lnM_M0, (1.0 - B_SZ), exp ((onemcor2 * w1 + COR * w2) * SD_SZ));
  gsl_vector_set (hx, 0, w1); //hx[0] = w1;
  gsl_vector_set (hx, 1, w2); //hx[1] = w2;
  gsl_vector_set (hx, 2, dw1); //hx[2] = dw1;
  gsl_vector_set (hx, 3, dw2); //hx[3] = dw2;

  /*printf ("Calling! % 20.15g % 20.15g % 20.15g % 20.15g\n", w1, w2, dw1, dw2);*/
  
  //printf ("gsl_f h0 = %.5g h1 = %.5g h2 = %.5g h3 = %.5g\n", hx[0], hx[1], hx[2], hx[3]);
}

/**
 * nc_cluster_mass_plcl_gsl_J_new_variables: (skip)
 * @p: FIXME
 * @j: FIXME 
 * @mszl: a #NcClusterMassPlCL
 * @lnM_M0: logarithm base e of the mass divided by the pivot mass
 * @Mobs: (array) (element-type double): observed mass
 * @Mobs_params: (array) (element-type double): observed mass paramaters
 *
 * FIXME
 *
*/
void
nc_cluster_mass_plcl_gsl_J_new_variables (const gsl_vector *p, gsl_matrix *j, NcClusterMassPlCL *mszl, gdouble lnM_M0, const gdouble *Mobs, const gdouble *Mobs_params)
{
  const gdouble onemcor2  = sqrt (1.0 - COR * COR);
  const gdouble w1 = gsl_vector_get (p, 0); //p[0];
  const gdouble w2 = gsl_vector_get (p, 1); //p[1];
  gdouble sd_PL, sd_CL, kernel_w1, kernel_w2;
   
  sd_PL = Mobs_params[NC_CLUSTER_MASS_PLCL_SD_PL];
  sd_CL = Mobs_params[NC_CLUSTER_MASS_PLCL_SD_CL];

  kernel_w1 = (1.0 - B_SZ) * exp (A_SZ * lnM_M0 + (onemcor2 * w1 + COR * w2) * SD_SZ) / sd_PL;
  kernel_w2 = (1.0 - B_L) * exp (A_L * lnM_M0 + w2 * SD_L) / sd_CL;
  
  gsl_matrix_set (j, 0, 0, 1.0); // j[0 * 2 + 0] = 1.0;
  gsl_matrix_set (j, 0, 1, 0.0); //j[0 * 2 + 1] = 0.0;
  gsl_matrix_set (j, 1, 0, 0.0); //j[1 * 2 + 0] = 0.0; 
  gsl_matrix_set (j, 1, 1, 1.0); //j[1 * 2 + 1] = 1.0;
  gsl_matrix_set (j, 2, 0, - kernel_w1 * onemcor2 * SD_SZ); //j[2 * 2 + 0] = - kernel_w1 * onemcor2 * SD_SZ;
  gsl_matrix_set (j, 2, 1, - kernel_w1 * COR * SD_SZ); //j[2 * 2 + 1] = - kernel_w1 * COR * SD_SZ;
  gsl_matrix_set (j, 3, 0, 0.0); //j[3 * 2 + 0] = 0.0;
  gsl_matrix_set (j, 3, 1, - kernel_w2 * SD_L); //j[3 * 2 + 1] = - kernel_w2 * SD_L;
  
  //printf ("j0 = %.5g j1 = %.5g j2 = %.5g j3 = %.5g j4 = %.5g j5 = %.5g j6 = %.5g j7 = %.5g\n", j[0*2 + 0], j[0*2 + 1])
  
}

/**
 * nc_cluster_mass_plcl_gsl_f:
 * @p: number of parameters 
 * @hx: components of the $\chi^2$ 
 * @n: number of elements of the $\chi^2$ 
 * @mszl: a #NcClusterMassPlCL
 * @lnM: logarithm base e of the mass divided by the pivot mass
 * @Mobs: (array) (element-type double): observed mass
 * @Mobs_params: (array) (element-type double): observed mass paramaters
 *
 * The $\chi^2$ is minimized with respect to the parameters $\ln\left(M_{SZ}/M_0\right)$ and 
 *  $\ln\left(M_{L}/M_0\right)$, therefore @p = 2. FIXME
 *
*/
void
nc_cluster_mass_plcl_gsl_f (const gdouble *p, gdouble *hx, gint n, NcClusterMassPlCL *mszl, gdouble lnM, const gdouble *Mobs, const gdouble *Mobs_params)
{
  integrand_data data;
  const gdouble onemcor2  = sqrt (1.0 - COR * COR);
  const gdouble lnMsz_M0  = p[0];
  const gdouble lnMl_M0   = p[1];
  gdouble Msz, Ml, dMsz, dMl, dlnMsz, dlnMl;
  
  data.mszl =            mszl;
  data.lnM =             lnM;
  data.mobs =            Mobs;
  data.mobs_params =     Mobs_params;
  data.mu_sz =           _SZ_lnmass_mean (data.mszl, data.lnM);
  data.mu_l =            _Lens_lnmass_mean (data.mszl, data.lnM); 
   
  Msz    = exp (lnMsz_M0) * mszl->M0;
  Ml     = exp (lnMl_M0) * mszl->M0;
  dMsz   = Msz - data.mobs[NC_CLUSTER_MASS_PLCL_MPL];
  dMl    = Ml - data.mobs[NC_CLUSTER_MASS_PLCL_MCL];
  dlnMsz = lnMsz_M0 - data.mu_sz;
  dlnMl  = lnMl_M0 - data.mu_l;

  //printf ("m = %d n = %d M0 = %.5g\n", m, n, mszl->M0);
  //printf ("p[0] = %.5g p[1] = %.5g lnM = %.5g\n", p[0], p[1], lnM); 
  //printf ("lnMsz/M0 = %.5g lnMl/M0 = %.5g\n", lnMsz_M0, lnMl_M0); 
  hx[0] = (COR * dlnMsz / SD_SZ - dlnMl / SD_L ) / onemcor2;
  hx[1] = dlnMsz / SD_SZ;
  hx[2] = dMsz / data.mobs_params[NC_CLUSTER_MASS_PLCL_SD_PL];
  hx[3] = dMl / data.mobs_params[NC_CLUSTER_MASS_PLCL_SD_CL];
  
  //printf ("gsl_f h0 = %.5g h1 = %.5g h2 = %.5g h3 = %.5g\n", hx[0], hx[1], hx[2], hx[3]);
}

static gint
_internal_nc_cluster_mass_plcl_gsl_f (const gsl_vector *p, gpointer adata, gsl_vector *hx) 
{
  integrand_data *data = (integrand_data *) adata;
  nc_cluster_mass_plcl_gsl_f (gsl_vector_const_ptr (p, 0), gsl_vector_ptr (hx, 0), 4, data->mszl, data->lnM, data->mobs, data->mobs_params);

  return GSL_SUCCESS;
}

static gint
_nc_cluster_mass_plcl_gsl_J (const gsl_vector *p, gpointer adata, gsl_matrix *J)
{
  integrand_data *data = (integrand_data *) adata;
  NcClusterMassPlCL *mszl = data->mszl;
  const gdouble onemcor2  = sqrt (1.0 - COR * COR);
  const gdouble lnMsz_M0  = gsl_vector_get (p, 0); 
  const gdouble lnMl_M0   = gsl_vector_get (p, 1); 
  const gdouble Msz       = exp (lnMsz_M0) * mszl->M0;
  const gdouble Ml        = exp (lnMl_M0) * mszl->M0;
  gsl_matrix_set (J, 0, 0, COR / (SD_SZ * onemcor2)); 
  gsl_matrix_set (J, 0, 1, - 1.0 / (SD_L * onemcor2)); 
  gsl_matrix_set (J, 1, 0, 1.0 / SD_SZ);               
  gsl_matrix_set (J, 1, 1, 0.0);                       
  gsl_matrix_set (J, 2, 0, Msz / data->mobs_params[NC_CLUSTER_MASS_PLCL_SD_PL]); 
  gsl_matrix_set (J, 2, 1, 0.0);                       
  gsl_matrix_set (J, 3, 0, 0.0);                      
  gsl_matrix_set (J, 3, 1, Ml / data->mobs_params[NC_CLUSTER_MASS_PLCL_SD_CL]); 
  
  //printf ("j0 = %.5g j1 = %.5g j2 = %.5g j3 = %.5g j4 = %.5g j5 = %.5g j6 = %.5g j7 = %.5g\n", j[0*2 + 0], j[0*2 + 1])

  return GSL_SUCCESS;
}

static void
_peakfinder (const gint *ndim, const gdouble bounds[], gint *n, gdouble x[], void *userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterMassPlCL *mszl = data->mszl;
  gsl_multifit_function_fdf f;

  gdouble p0[] = {log (data->mobs[NC_CLUSTER_MASS_PLCL_MPL]/mszl->M0), log (data->mobs[NC_CLUSTER_MASS_PLCL_MCL]/mszl->M0)};
  gdouble lb[] = {bounds[0], bounds[2]};
  gdouble ub[] = {bounds[1], bounds[3]};

#ifdef HAVE_GSL_2_2
  gint status;
  gint info;
  const gdouble xtol = 1e-11;
  const gdouble gtol = 1e-11;
  const gdouble ftol = 0.0;
#endif /* HAVE_GSL_2_2 */

  p0[0] = GSL_MAX (p0[0], lb[0]);
  p0[1] = GSL_MAX (p0[1], lb[1]);
  p0[0] = GSL_MIN (p0[0], ub[0]);
  p0[1] = GSL_MIN (p0[1], ub[1]);

  gsl_vector_view p0_vec = gsl_vector_view_array (p0, 2);
  
  f.f  = &_internal_nc_cluster_mass_plcl_gsl_f;
  f.df = &_nc_cluster_mass_plcl_gsl_J; 
  f.n = 4;
  f.p = 2;
  f.params = data;

  /* initialize solver with starting point */
  gsl_multifit_fdfsolver_set (mszl->s, &f, &p0_vec.vector);

  /* solve the system with a maximum of 20 iterations */
#ifdef HAVE_GSL_2_2
  status = gsl_multifit_fdfsolver_driver (mszl->s, 20000, xtol, gtol, ftol, &info);
  if (status != GSL_SUCCESS)
    g_error ("_peakfinder: NcClusterMassPlCL peakfinder function.\n");
#else /* HAVE_GSL_2_2 */
  g_error ("_peakfinder: this model requires gsl >= 2.2.");
#endif /* HAVE_GSL_2_2 */

  //printf ("Inicial: p0 = %.5g p1 = %.5g\n", p0[0], p0[1]);
  //printf ("2d Minimo : p  = %.5g p  = %.5g\n", gsl_vector_get (mszl->s->x, 0), gsl_vector_get (mszl->s->x, 1));
  x[0] = gsl_vector_get (mszl->s->x, 0);
  x[1] = gsl_vector_get (mszl->s->x, 1);

  *n = 1;
}

static gdouble
_function_at (gdouble *p, integrand_data *data)
{
  gdouble fp[4];
  gsl_vector_view fp_vec = gsl_vector_view_array (fp, 4);
  gsl_vector_view p_vec  = gsl_vector_view_array (p, 2);
  _internal_nc_cluster_mass_plcl_gsl_f (&p_vec.vector, data, &fp_vec.vector);
  return (fp[0] * fp[0] + fp[1] * fp[1] + fp[2] * fp[2] + fp[3] * fp[3]) / 2.0;
}

#define PEAK_DEC (100.0 * M_LN10)

static gdouble
_function_bounds_Msz (gdouble x, void *params)
{
  integrand_data *data = (integrand_data *) params;
  gdouble p[2] = {x, data->peak[1]};
  return _function_at (p, data) - (data->func_peak + PEAK_DEC);
}

static gdouble
_function_bounds_Ml (gdouble x, void *params)
{
  integrand_data *data = (integrand_data *) params;
  gdouble p[2] = {data->peak[0], x};
  return _function_at (p, data) - (data->func_peak + PEAK_DEC);
}

static gdouble
_border_finder (gsl_function *F, gdouble x_lo, gdouble x_hi, gdouble prec, guint max_iter)
{
  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
  gdouble r = 0.0;
  guint iter = 0;
  gint status;

  //printf ("x_lo = %.5g F em x_lo = %.5g x_hi= %.5g F em x_hi = %.5g\n", x_lo, F->function (x_lo, F->params), x_hi, F->function (x_hi, F->params));
  gsl_root_fsolver_set (s, F, x_lo, x_hi);

  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    r = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, 0, prec);
    if (status == GSL_SUCCESS)
      break;
  } while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  return r;
}

static void
_p_bounds (gdouble *lb, gdouble *ub, integrand_data *data)
{
  gint max_iter = 1000000;
  gdouble prec = 1e-1;
  gsl_function F;

  F.params = data;

  F.function = &_function_bounds_Msz;
  lb[0] = _border_finder (&F, -10.0, data->peak[0], prec, max_iter);
  ub[0] = _border_finder (&F, data->peak[0], 60.0, prec, max_iter);

  F.function = &_function_bounds_Ml;
  lb[1] = _border_finder (&F, -10.0, data->peak[1], prec, max_iter);
  ub[1] = _border_finder (&F, data->peak[1], 60.0, prec, max_iter);

  if (FALSE)
  {
    gdouble p[2];
    p[0] = lb[0];
    p[1] = lb[1];
    printf ("% 20.15g % 20.15g % 20.15g\n", p[0], p[1], _function_at (p, data) - (data->func_peak + PEAK_DEC));
    p[0] = lb[0];
    p[1] = ub[1];
    printf ("% 20.15g % 20.15g % 20.15g\n", p[0], p[1], _function_at (p, data) - (data->func_peak + PEAK_DEC));
    p[0] = ub[0];
    p[1] = lb[1];
    printf ("% 20.15g % 20.15g % 20.15g\n", p[0], p[1], _function_at (p, data) - (data->func_peak + PEAK_DEC));
    p[0] = ub[0];
    p[1] = ub[1];
    printf ("% 20.15g % 20.15g % 20.15g\n", p[0], p[1], _function_at (p, data) - (data->func_peak + PEAK_DEC));
  }

}

static gdouble
_nc_cluster_mass_plcl_Msz_Ml_M500_p (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *Mobs, const gdouble *Mobs_params)
{
  integrand_data data;
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (clusterm);
  gdouble sd_Pl, sd_CL; 
  const gdouble four_pi2 = 4.0 * M_PI * M_PI;
  const gdouble sdsz_sdl = SD_SZ * SD_L;
  const gdouble norma_factor =  sdsz_sdl * sqrt (1.0 - COR * COR);
  NcmIntegrand2dim integ;
  gdouble P, err;

  data.mszl =            mszl;
  data.lnM =             lnM;
  data.mobs =            Mobs;
  data.mobs_params =     Mobs_params;
  data.mu_sz =           _SZ_lnmass_mean (data.mszl, lnM);
  data.mu_l =            _Lens_lnmass_mean (data.mszl, lnM);
  data.twocor_sdsz_sdl = 2.0 * COR / sdsz_sdl;

  data.sd_sz2 =          SD_SZ * SD_SZ;
  data.sd_l2 =           SD_L * SD_L;

  sd_Pl = data.mobs_params[NC_CLUSTER_MASS_PLCL_SD_PL];
  sd_CL = data.mobs_params[NC_CLUSTER_MASS_PLCL_SD_CL];
  
  // Dividing norma by mszl->M0 in order to have dimensionless 
  // sd_PL and sd_CL, since they are given in units of 10^{14} h^{-1} M_solar, and given that 
  // P(M_PL, M_CL; M500) dM_PL dM_CL and, therefore, the units are canceled.
  data.lnnorma_p = log (four_pi2 * norma_factor * sd_Pl * sd_CL / (mszl->M0 * mszl->M0));
  //data.lnnorma_p = log (2.0 * M_PI * norma_factor);
  //data.lnnorma_p = log (2.0 * M_PI * sd_Pl * sd_CL / (mszl->M0 * mszl->M0));
  
  integ.f = _nc_cluster_mass_plcl_Msz_Ml_M500_p_integrand;
  integ.userdata = &data;

  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);

  {
    gint n = 1;
    gint ndim = 2;
    gdouble a_sz, a_l, b_sz, b_l, bounds[ndim * 2], x[ndim];
    gdouble lb[2], ub[2]; 
    const gint ngiven   = 1;
    const gint ldxgiven = 2;
    //gdouble xgiven[ldxgiven * ngiven];
    
    a_sz = -10.0;  
    b_sz = 60.0; 
    a_l  = -10.0; 
    b_l  = 60.0; 
    bounds[0] = a_sz;
    bounds[2] = a_l; 
    bounds[1] = b_sz; 
    bounds[3] = b_l;

    //printf ("[%.5g %.5g] [%.8g %.8g]\n", a_sz, b_sz, a_l, b_l);
    //printf ("exp: asz = %.8g bsz = %.8g al = %.8g bl = %.8g\n", exp(a_sz), exp(b_sz), exp(a_l), exp(b_l));
    _peakfinder (&ndim, bounds, &n, x, &data);
    //printf ("start %20.15g %20.15g [%20.15g %20.15g, %20.15g %20.15g]\n", x[0], x[1], a_sz, a_l, b_sz, b_l);
    
    data.peak[0] = x[0];
    data.peak[1] = x[1];
    data.func_peak = _function_at (x, &data);
    _p_bounds (lb, ub, &data);
    //printf ("[%.5g, %.5g] [%.5g, %.5g]\n", lb[0], ub[0], lb[1], ub[1]);

    //lb[0] = 1.8716; 
    //ub[0] = 2.2833;
    //lb[1] = -3.0626;  
    //ub[1] = 3.0843;
    
    //ncm_integrate_2dim_divonne (&integ, a_sz, a_l, b_sz, b_l, 1.0e-5, 0.0, ngiven, ldxgiven, x, &P, &err);
    ncm_integrate_2dim_divonne (&integ, lb[0], lb[1], ub[0], ub[1], 1.0e-5, 0.0, ngiven, ldxgiven, x, &P, &err);
 
    //printf ("P1 % 20.15g P2 % 20.15g P3 % 20.15g | % 20.15g % 20.15g <<% 20.15g, % 20.15g, % 20.15g>>\n", P, P2, P3, P/P3, P2/P3, err/P, err2/P2, err3/P3);
 
  }

  //printf("M500 = %.10g P = %.10g err = %.10g\n", exp(lnM), P, err);
  return P;
}

/* Functions to compute the integrand and the 2D integral over ln(M_SZ/M0) and ln(M_L/M0) when ndet = 1.  See paper! */

static gdouble
_nc_cluster_mass_plcl_Msz_Ml_p_ndetone_integrand (gdouble lnMsz_M0, gdouble lnMl_M0, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterMassPlCL *mszl = data->mszl;  
  const gdouble Msz_M0 = exp (lnMsz_M0);
  const gdouble Ml_M0  = exp (lnMl_M0);
  const gdouble lnMsz_mlnBsz = lnMsz_M0 - log1p (- B_SZ);
  const gdouble lnMl_mlnBl   = lnMl_M0 - log1p (- B_L);
  
  const gdouble M_Pl_M0  = data->mobs[NC_CLUSTER_MASS_PLCL_MPL];
  const gdouble sd_Pl_M0 = data->mobs_params[NC_CLUSTER_MASS_PLCL_SD_PL];
  const gdouble M_CL_M0  = data->mobs[NC_CLUSTER_MASS_PLCL_MCL];
  const gdouble sd_CL_M0 = data->mobs_params[NC_CLUSTER_MASS_PLCL_SD_CL];
  
  const gdouble ysz     = (M_Pl_M0 - Msz_M0) / sd_Pl_M0;
  const gdouble arg_ysz = ysz * ysz / 2.0;
  const gdouble yl      = (M_CL_M0 - Ml_M0) / sd_CL_M0;
  const gdouble arg_yl  = yl * yl / 2.0;
  
  const gdouble xsz       = A_L * lnMsz_mlnBsz;
  const gdouble xl        =  A_SZ * lnMl_mlnBl;
  const gdouble diff_xl_xsz = xl - xsz;
  const gdouble diff_szl2    = diff_xl_xsz * diff_xl_xsz;
  const gdouble arg_xszl = diff_szl2 / data->c1;

  const gdouble exp_arg = - arg_ysz - arg_yl - arg_xszl;

  const gdouble erf_arg     = data->erf_const_Msz * lnMsz_mlnBsz + data->erf_const_Ml * lnMl_mlnBl;
  const gdouble erf_arg_low = (erf_arg + data->erf_Mlow) / data->c2;
  const gdouble erf_arg_up  = (erf_arg + data->erf_Mup) / data->c2; 

  //printf ("%20.15g %20.15g %20.15g %20.15g\n", lnMsz, lnMl, exp_arg, exp (exp_arg - data->lnnorma_p));
  //printf ("Massas: %8.5g %8.5g %8.5e %8.5e %8.5g %8.5g %8.5e %8.5e\n", lnMsz, data->mu_sz, M_Pl, Msz, lnMl, data->mu_l, M_CL, Ml);
  if (exp_arg < GSL_LOG_DBL_MIN)
    return exp (-200.0);
  else
  {
    const gdouble result = exp (exp_arg - data->lnnorma_p) * (erf(erf_arg_up) - erf (erf_arg_low)) + exp (-200.0);
    //printf ("===> %12.8g %12.8g %12.8g %12.8g\n", data->lnM, lnMsz, lnMl, result);
    return result;
  }
}

/**
 * nc_cluster_mass_plcl_Msz_Ml_p_ndetone:
 * @clusterm: a #NcClusterMass
 * @lnMcut: lower threshold of the true mass
 * @z: redshift
 * @Mpl: Planck cluster mass
 * @Mcl: CLASH cluster mass
 * @sigma_pl: Planck mass error
 * @sigma_cl: CLASH mass error
 *
 * This function computes the i-th term of the posterior given flat priors for 
 * the selection function and mass function. See function nc_cluster_pseudo_counts_posterior_ndetone().
 *
 * Returns: FIXME
*/
gdouble
nc_cluster_mass_plcl_Msz_Ml_p_ndetone (NcClusterMass *clusterm, gdouble lnMcut, const gdouble z, const gdouble Mpl, const gdouble Mcl, const gdouble sigma_pl, const gdouble sigma_cl)
{
  integrand_data data;
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (clusterm);
  const gdouble AszSDl = A_SZ * SD_L;
  const gdouble AlSDsz = A_L * SD_SZ;
  const gdouble term1 = AszSDl * AszSDl + AlSDsz * AlSDsz - 2.0 * AszSDl * AlSDsz * COR;
  const gdouble norma_factor = 2.0 * M_SQRT2 * M_SQRTPI * sqrt (term1);
  const gdouble c1 = 2.0 * term1;
  const gdouble c2 = SD_SZ * SD_L * sqrt(1.0 - COR * COR) * sqrt (c1);
  const gdouble erf_Mlow = term1 * (lnMcut - log(mszl->M0));
  const gdouble erf_Mup = term1 * log(1.0e16 / mszl->M0);
  const gdouble erf_const_Msz = (AlSDsz * COR - AszSDl) * SD_L;
  const gdouble erf_const_Ml = (AszSDl * COR - AlSDsz) * SD_SZ;
  const gdouble Mobs[] = {Mpl / mszl->M0, Mcl / mszl->M0};
  const gdouble Mobs_params[] = {sigma_pl / mszl->M0, sigma_cl / mszl->M0};
  /* Normalization of the prior of the selection function: true mass top-hat function between Mcut and 10^{16} h^{-1} Msun. */
  const gdouble norm_Mtrue = 16.0 * M_LN10 - lnMcut; 
  gdouble sd_Pl, sd_CL;
  gdouble P, err;
  NcmIntegrand2dim integ;

  data.mszl          = mszl;
  data.mobs          = Mobs;
  data.mobs_params   = Mobs_params;
  data.c1            = c1;
  data.c2            = c2;
  data.erf_Mlow      = erf_Mlow;
  data.erf_Mup       = erf_Mup;
  data.erf_const_Msz = erf_const_Msz;
  data.erf_const_Ml  = erf_const_Ml;

  sd_Pl = data.mobs_params[NC_CLUSTER_MASS_PLCL_SD_PL];
  sd_CL = data.mobs_params[NC_CLUSTER_MASS_PLCL_SD_CL];
  
  data.lnnorma_p = log (2.0 * M_PI * norma_factor * sd_Pl * sd_CL);
    
  integ.f = _nc_cluster_mass_plcl_Msz_Ml_p_ndetone_integrand;
  integ.userdata = &data;

  NCM_UNUSED (z);

  gdouble a_sz, a_l, b_sz, b_l;
  a_sz = a_l = log(1.0e12 / mszl->M0);
  b_sz = b_l = log(1.0e16 / mszl->M0); 
  ncm_integrate_2dim (&integ, a_sz, a_l, b_sz, b_l, NCM_DEFAULT_PRECISION, 0.0, &P, &err);
    
  return P / norm_Mtrue;
  
}

// integrand to compute 
static gdouble
_nc_cluster_mass_plcl_int_Mobs_cut_inf (gdouble Msz_l, gdouble Mcut, gdouble sigma_Mobs)
{
  const gdouble a = (Mcut - Msz_l) / (M_SQRT2 * sigma_Mobs);

  if (a < 0.0)
    return (1.0 - erf (a)) * 0.5;
  else
    return erfc (a) * 0.5;
}

static gdouble
_nc_cluster_mass_plcl_Msz_Ml_M500_intp_integrand (gdouble lnMsz, gdouble lnMl, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterMassPlCL *mszl = data->mszl;
  const gdouble Msz = exp (lnMsz);
  const gdouble Ml = exp (lnMl);
  const gdouble diff_Msz = lnMsz - data->mu_sz;
  const gdouble diff_Ml = lnMl - data->mu_l;
  const gdouble Mcut = data->Mcut;
  
  const gdouble sd_Pl = 0.2; //data->mobs_params[NC_CLUSTER_MASS_PLCL_SD_PL];
  const gdouble sd_CL = 0.2; //data->mobs_params[NC_CLUSTER_MASS_PLCL_SD_CL];
  
  const gdouble intM_Pl = _nc_cluster_mass_plcl_int_Mobs_cut_inf (Msz, Mcut, sd_Pl);
  const gdouble intM_CL = _nc_cluster_mass_plcl_int_Mobs_cut_inf (Ml, Mcut, sd_CL);
   
  const gdouble xsz = diff_Msz / SD_SZ;
  const gdouble arg_xsz = xsz * xsz;
  const gdouble xl =  diff_Ml / SD_L;
  const gdouble arg_xl = xl * xl;
  const gdouble arg_x_szl = data->twocor_sdsz_sdl * diff_Msz * diff_Ml;

  const gdouble exp_arg = - (arg_xsz + arg_xl - arg_x_szl) / (2.0 * (1.0 - COR * COR));
  if (exp_arg < GSL_LOG_DBL_MIN)
    return 0.0;
  else
  {
    const gdouble result = intM_Pl * intM_CL * exp (exp_arg - data->lnnorma_p);
    return result;
  } 
}


static gdouble
_nc_cluster_mass_plcl_intp (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  integrand_data data;
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (clusterm);
  //gdouble sd_Pl, sd_CL; 
  const gdouble sdsz_sdl = SD_SZ * SD_L;
  const gdouble norma_factor =  sdsz_sdl * sqrt(1.0 - COR * COR);
  gdouble P, err;
  NcmIntegrand2dim integ;

  data.mszl =            mszl;
  data.lnM =             lnM;
  data.mobs_params =     NULL;
  data.mu_sz =           _SZ_lnmass_mean (data.mszl, lnM);
  data.mu_l =            _Lens_lnmass_mean (data.mszl, lnM);
  data.twocor_sdsz_sdl = 2.0 * COR / sdsz_sdl;
  data.sd_sz2 =          SD_SZ * SD_SZ;
  data.sd_l2 =           SD_L * SD_L;
  data.lnnorma_p =       log ((2.0 * M_PI *  norma_factor));
  
  //sd_Pl = data.mobs_params[NC_CLUSTER_MASS_PLCL_SD_PL];
  //sd_CL = data.mobs_params[NC_CLUSTER_MASS_PLCL_SD_CL];
  
  integ.f = _nc_cluster_mass_plcl_Msz_Ml_M500_intp_integrand;
  integ.userdata = &data;

  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);

  gdouble a_sz, a_l, b_sz, b_l;
  a_sz = a_l = log (1.0e10);
  b_sz = b_l = log (1.0e17); 
  ncm_integrate_2dim (&integ, a_sz, a_l, b_sz, b_l, NCM_DEFAULT_PRECISION, 0.0, &P, &err);
    
  return P;
  
}

static gboolean
_nc_cluster_mass_plcl_resample (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnMobs, const gdouble *lnMobs_params, NcmRNG *rng)
{
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (clusterm);
  gdouble r_SZ, r_L, sz_ran, l_ran, M_pl, M_cl;
  gdouble lnM_SZ_ran, lnM_L_ran, M_SZ_ran, M_L_ran;
  const gdouble lnM_SZ = _SZ_lnmass_mean (mszl, lnM);
  const gdouble lnM_L = _Lens_lnmass_mean (mszl, lnM);

  NCM_UNUSED (clusterm);
  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);
  
  ncm_rng_lock (rng);
  gsl_ran_bivariate_gaussian (rng->r, SD_SZ, SD_L, COR, &r_SZ, &r_L);
  lnM_SZ_ran = lnM_SZ + r_SZ;
  lnM_L_ran  = lnM_L + r_L;
  ncm_rng_unlock (rng);

  M_SZ_ran = exp (lnM_SZ_ran) * mszl->M0; /* M is in units of M0, i.e., M -> M/M0 */
  M_L_ran = exp (lnM_L_ran) * mszl->M0;
  
  ncm_rng_lock (rng);
  do {
    sz_ran  = gsl_ran_gaussian (rng->r, lnMobs_params[0]);
    M_pl = M_SZ_ran + sz_ran;
  } while (M_pl < 0.0);
  /*printf ("%8.5f % 20.15e % 20.15e % 20.15f\n", z, M_SZ_ran, sz_ran, 100.0 * fabs (sz_ran / M_SZ_ran));*/
  do {
    l_ran   = gsl_ran_gaussian (rng->r, lnMobs_params[1]);
    M_cl = M_L_ran + l_ran;
  } while (M_cl < 0.0);
  /*printf ("%8.5f % 20.15e % 20.15e % 20.15f\n", z, M_L_ran, l_ran, 100.0 * fabs (l_ran / M_L_ran));*/
  ncm_rng_unlock (rng);
  
  lnMobs[NC_CLUSTER_MASS_PLCL_MPL] = M_pl;  
  lnMobs[NC_CLUSTER_MASS_PLCL_MCL] = M_cl;

  return TRUE; 
  /* Information to select these masses are not implemented here. 
   * The selection function is in nc_cluster_pseudo_counts
   * return (lnMobs[NC_CLUSTER_MASS_PLCL_MPL] >= mszl->lnMobs_min) && (lnMobs[NC_CLUSTER_MASS_PLCL_MCL] >= mszl->lnMobs_min);
   */
}

/*
static void
_nc_cluster_mass_plcl_p_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble *xi, gdouble *xi_params, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (clusterm);
  const gdouble xil = GSL_MAX (xi[0] - 7.0, msz->signif_obs_min);
  const gdouble zetal = _significance_to_zeta (clusterm, model, 2.0, xil) - 7.0 * D_SZ;
  const gdouble lnMl = GSL_MAX (_zeta_to_mass (clusterm, model, 2.0, zetal), log (NC_CLUSTER_MASS_BENSON_M_LOWER_BOUND));

  const gdouble xiu = xi[0] + 7.0;
  const gdouble zetau = _significance_to_zeta (clusterm, model, 0.0, xiu) + 7.0 * D_SZ;
  const gdouble lnMu = _zeta_to_mass (clusterm, model, 0.0, zetau);

  NCM_UNUSED (xi_params);
  
  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;
}
*/


static void
_nc_cluster_mass_plcl_n_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (clusterm);

  NCM_UNUSED (cosmo);
  NCM_UNUSED (mszl);
  
  *lnM_lower = 12.0 * M_LN10; /* This value must be compatible with M_cut and sigma_cut from cluster pseudo counts */ 
  *lnM_upper = 16.0 * M_LN10; /* Largest mass compatible with Tinker multiplicity function */

  return;  
}

