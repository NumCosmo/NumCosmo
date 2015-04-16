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
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass_plcl.h"
#include "math/integral.h"
#include "math/memory_pool.h"
#include "math/ncm_cfg.h"

#include <gsl/gsl_randist.h>

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
  PROP_SIZE,
};

static void
nc_cluster_mass_plcl_init (NcClusterMassPlCL *mszl)
{
  NCM_UNUSED (mszl);
}

static void
_nc_cluster_mass_plcl_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (object);
  g_return_if_fail (NC_IS_CLUSTER_MASS_PLCL (object));

  NCM_UNUSED (mszl);
/*
  switch (prop_id)
  {
    case PROP_:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }*/
}

static void
_nc_cluster_mass_plcl_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (object);
  g_return_if_fail (NC_IS_CLUSTER_MASS_PLCL (object));

  NCM_UNUSED (mszl);
/*
  switch (prop_id)
  {
    case PROP_:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
  */
}

static void
_nc_cluster_mass_plcl_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_mass_plcl_parent_class)->finalize (object);
}

guint _nc_cluster_mass_plcl_obs_len (NcClusterMass *clusterm) { NCM_UNUSED (clusterm); return 2; }
guint _nc_cluster_mass_plcl_obs_params_len (NcClusterMass *clusterm) { NCM_UNUSED (clusterm); return 2; }
static gdouble _nc_cluster_mass_plcl_Msz_Ml_M500_p (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *Mobs, const gdouble *Mobs_params);
static gdouble _nc_cluster_mass_plcl_intp (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z);
static gboolean _nc_cluster_mass_plcl_resample (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnMobs, const gdouble *lnMobs_params, NcmRNG *rng);

static void
nc_cluster_mass_plcl_class_init (NcClusterMassPlCLClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcClusterMassClass* parent_class = NC_CLUSTER_MASS_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_cluster_mass_plcl_set_property;
  model_class->get_property = &_nc_cluster_mass_plcl_get_property;
  object_class->finalize    = &_nc_cluster_mass_plcl_finalize;

  ncm_model_class_set_name_nick (model_class, "Planck - CLASH cluster mass distribution", "Planck_CLASH");
  ncm_model_class_add_params (model_class, NC_CLUSTER_MASS_PLCL_SPARAM_LEN, 0, PROP_SIZE);

  /*
   * SZ signal-mass scaling parameter: Asz.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PLCL_A_SZ, "\\alpha_{SZ}", "Asz",
                              1e-5,  10.0, 1.0e-1,
                              NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PLCL_DEFAULT_A_SZ,
                              NCM_PARAM_TYPE_FIXED);

  /*
   * SZ signal-mass scaling parameter: Bsz.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PLCL_B_SZ, "b_{SZ}", "Bsz",
                              -2.0,  0.9999, 1.0e-1,
                              NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PLCL_DEFAULT_B_SZ,
                              NCM_PARAM_TYPE_FIXED);

  /*
   * SZ signal-mass scaling parameter: SDsz.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PLCL_SD_SZ, "\\sigma_{SZ}", "sigma_sz",
                              1e-5,  2.0, 1.0e-1,
                              NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PLCL_DEFAULT_SD_SZ,
                              NCM_PARAM_TYPE_FIXED);
  /*
   * Lensing signal-mass scaling parameter: Al.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PLCL_A_L, "\\alpha_{L}", "Al",
                              1e-5,  10.0, 1.0e-1,
                              NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PLCL_DEFAULT_A_L,
                              NCM_PARAM_TYPE_FIXED);

  /*
   * Lensing signal-mass scaling parameter: Bl.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PLCL_B_L, "b_{L}", "Bl",
                              0.0,  10.0, 1.0e-1,
                              NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PLCL_DEFAULT_B_L,
                              NCM_PARAM_TYPE_FIXED);

  /*
   * Lensing signal-mass scaling parameter: SDl.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PLCL_SD_L, "\\sigma_{L}", "sigma_l",
                              1e-5,  2.0, 1.0e-1,
                              NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PLCL_DEFAULT_SD_L,
                              NCM_PARAM_TYPE_FIXED);
  /*
   * SZ-Lensing signal-mass correlation: COR.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_PLCL_COR, "\\rho", "cor",
                              1e-4,  1.0, 1.0e-1,
                              NC_CLUSTER_MASS_PLCL_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_PLCL_DEFAULT_COR,
                              NCM_PARAM_TYPE_FIXED);
    
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);  

  parent_class->P = &_nc_cluster_mass_plcl_Msz_Ml_M500_p;
  parent_class->intP = &_nc_cluster_mass_plcl_intp;
  //parent_class->P_limits = &_nc_cluster_mass_plcl_p_limits;
  //parent_class->N_limits = &_nc_cluster_mass_plcl_n_limits;
  parent_class->resample = &_nc_cluster_mass_plcl_resample;
  parent_class->obs_len = &_nc_cluster_mass_plcl_obs_len;
  parent_class->obs_params_len = &_nc_cluster_mass_plcl_obs_params_len;

  parent_class->impl = NC_CLUSTER_MASS_IMPL_ALL;
}

typedef struct _integrand_data
{
  NcClusterMassPlCL *mszl;
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
} integrand_data;

static gdouble
_SZ_lnmass_mean (NcClusterMassPlCL *mszl,  gdouble lnM)
{
  const gdouble lnMsz_mean = log1p (- B_SZ) + A_SZ * lnM + (1.0 - A_SZ) * 14.0 * M_LN10;
  return lnMsz_mean;
}

static gdouble
_Lens_lnmass_mean (NcClusterMassPlCL *mszl, gdouble lnM)
{
  const gdouble lnMlens_mean = B_L + A_L * lnM + (1.0 - A_L) * 14.0 * M_LN10;
  return lnMlens_mean;
}

/**
 * nc_cluster_mass_plcl_pdf:
 * @clusterm: a #NcClusterMass
 * @lnM500:logarithm base e of the true mass
 * @lnMsz: logarithm base e of the SZ mass
 * @lnMl: logarithm base e of the lensing mass
 * @Mobs: (array) (element-type double): observed masses
 * @Mobs_params: (array) (element-type double): observed mass paramaters
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble 
nc_cluster_mass_plcl_pdf (NcClusterMass *clusterm, gdouble lnM500, gdouble lnMsz, gdouble lnMl, const gdouble *Mobs, const gdouble *Mobs_params)
{
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (clusterm);
  const gdouble Msz = exp (lnMsz);
  const gdouble Ml  = exp (lnMl);
  const gdouble diff_Msz = lnMsz - _SZ_lnmass_mean (mszl, lnM500);
  const gdouble diff_Ml  = lnMl - _Lens_lnmass_mean (mszl, lnM500);
  
  const gdouble M_Pl  = Mobs[NC_CLUSTER_MASS_PLCL_MPL];
  const gdouble sd_Pl = Mobs_params[NC_CLUSTER_MASS_PLCL_SD_PL];
  const gdouble M_CL  = Mobs[NC_CLUSTER_MASS_PLCL_MCL];
  const gdouble sd_CL = Mobs_params[NC_CLUSTER_MASS_PLCL_SD_CL];
  
  const gdouble ysz     = (M_Pl - Msz) / sd_Pl;
  const gdouble arg_ysz = ysz * ysz / 2.0;
  const gdouble yl      = (M_CL - Ml) / sd_CL;
  const gdouble arg_yl  = yl * yl / 2.0;
  
  const gdouble xsz             = diff_Msz / SD_SZ;
  const gdouble arg_xsz         = xsz * xsz;
  const gdouble xl              =  diff_Ml / SD_L;
  const gdouble arg_xl          = xl * xl;
  const gdouble sdsz_sdl        = SD_SZ * SD_L;
  const gdouble twocor_sdsz_sdl = 2.0 * COR / sdsz_sdl;
  const gdouble arg_x_szl       = twocor_sdsz_sdl * diff_Msz * diff_Ml;

  const gdouble norma_factor =  sdsz_sdl * sqrt(1.0 - COR * COR);
  const gdouble lnnorma_p    = log (4.0 * M_PI * M_PI * (sd_Pl * sd_CL * 1.0e-28) *  norma_factor);
  
  const gdouble exp_arg = - arg_ysz - arg_yl - (arg_xsz + arg_xl - arg_x_szl) / (2.0 * (1.0 - COR * COR));
  if (exp_arg < GSL_LOG_DBL_MIN)
    return 0.0;
  else
  {
    const gdouble result = exp (exp_arg - lnnorma_p);
    return result;
  }
}

static gdouble
_nc_cluster_mass_plcl_Msz_Ml_M500_p_integrand (gdouble lnMsz, gdouble lnMl, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterMassPlCL *mszl = data->mszl;
  const gdouble Msz = exp (lnMsz);
  const gdouble Ml  = exp (lnMl);
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
  const gdouble xl        =  diff_Ml / SD_L;
  const gdouble arg_xl    = xl * xl;
  const gdouble arg_x_szl = data->twocor_sdsz_sdl * diff_Msz * diff_Ml;

  const gdouble exp_arg = - arg_ysz - arg_yl - (arg_xsz + arg_xl - arg_x_szl) / (2.0 * (1.0 - COR * COR));
  if (exp_arg < GSL_LOG_DBL_MIN)
    return 0.0;
  else
  {
    const gdouble result = exp (exp_arg - data->lnnorma_p);
    return result;
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
  const gdouble norma_factor =  sdsz_sdl * sqrt(1.0 - COR * COR);
  gdouble P, err;
  NcmIntegrand2dim integ;

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
  
  // Multiplying norma by 1.0e^-14 * 1.0e^-14 in order to have adimensional 
  // sd_PL and sd_CL, since they are given in units of 10^{14} h^{-1} M_solar, and given that 
  // P(M_PL, M_CL; M500) dM_PL dM_CL and, therefore, the units are canceled.
  data.lnnorma_p = log ((four_pi2 * (sd_Pl * sd_CL * 1.0e-28) *  norma_factor));
  
  integ.f = _nc_cluster_mass_plcl_Msz_Ml_M500_p_integrand;
  integ.userdata = &data;

  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);

  {
    gdouble a_sz, a_l, b_sz, b_l;
    a_sz = log (fabs (Mobs[NC_CLUSTER_MASS_PLCL_MPL] - 8.0 * sd_Pl));
    b_sz = log (Mobs[NC_CLUSTER_MASS_PLCL_MPL] + 8.0 * sd_Pl);
    a_l = log (fabs (Mobs[NC_CLUSTER_MASS_PLCL_MCL] - 8.0 * sd_CL));
    b_l = log (Mobs[NC_CLUSTER_MASS_PLCL_MCL] + 8.0 * sd_CL);

    ncm_integrate_2dim (&integ, a_sz, a_l, b_sz, b_l, NCM_DEFAULT_PRECISION, 0.0, &P, &err);
  }
  /*
  {
    // Testing the precision of the integral in 4 subregions: in the cases I've tested the improvement is irrelevant
    gdouble a_sz, a_l, b_sz, b_l, P1, P2, P3, P4, err1, err2, err3, err4;

    {
      a_sz = log(Mobs[NC_CLUSTER_MASS_PLCL_MPL]);
      b_sz = log(Mobs[NC_CLUSTER_MASS_PLCL_MPL] + 10.0 * sd_Pl);
      a_l = log(Mobs[NC_CLUSTER_MASS_PLCL_MCL]);
      b_l = log(Mobs[NC_CLUSTER_MASS_PLCL_MCL] + 10.0 * sd_CL);

      printf("asz = %.5e bsz = %.5e al = %.5e bl = %.5e\n", a_sz, b_sz, a_l, b_l);
      ncm_integrate_2dim (&integ, a_sz, a_l, b_sz, b_l, NCM_DEFAULT_PRECISION, 0.0, &P1, &err1);
    }
    {
      a_sz = log(fabs(Mobs[NC_CLUSTER_MASS_PLCL_MPL] - 10.0 * sd_Pl));
      b_sz = log(Mobs[NC_CLUSTER_MASS_PLCL_MPL]);
      a_l = log(Mobs[NC_CLUSTER_MASS_PLCL_MCL]);
      b_l = log(Mobs[NC_CLUSTER_MASS_PLCL_MCL] + 10.0 * sd_CL);

      printf("asz = %.5e bsz = %.5e al = %.5e bl = %.5e\n", a_sz, b_sz, a_l, b_l);
      ncm_integrate_2dim (&integ, a_sz, a_l, b_sz, b_l, NCM_DEFAULT_PRECISION, 0.0, &P2, &err2);
    }
    {
      a_sz = log(fabs(Mobs[NC_CLUSTER_MASS_PLCL_MPL] - 10.0 * sd_Pl));
      b_sz = log(Mobs[NC_CLUSTER_MASS_PLCL_MPL]);
      a_l = log(fabs(Mobs[NC_CLUSTER_MASS_PLCL_MCL] - 10.0 * sd_CL));
      b_l = log(Mobs[NC_CLUSTER_MASS_PLCL_MCL]);

      printf("asz = %.5e bsz = %.5e al = %.5e bl = %.5e\n", a_sz, b_sz, a_l, b_l);
      ncm_integrate_2dim (&integ, a_sz, a_l, b_sz, b_l, NCM_DEFAULT_PRECISION, 0.0, &P3, &err3);
    }
    {
      a_sz = log(Mobs[NC_CLUSTER_MASS_PLCL_MPL]);
      b_sz = log(Mobs[NC_CLUSTER_MASS_PLCL_MPL] + 10.0 * sd_Pl);
      a_l = log(fabs(Mobs[NC_CLUSTER_MASS_PLCL_MCL] - 10.0 * sd_CL));
      b_l = log(Mobs[NC_CLUSTER_MASS_PLCL_MCL]);

      //printf("M_Pl = %.5e 3sd_Pl = %.5e\n", Mobs[NC_CLUSTER_MASS_PLCL_MPL], 3.0 * sd_Pl);
      //printf("M_CL = %.5e 3sd_CL = %.5e\n", Mobs[NC_CLUSTER_MASS_PLCL_MCL], 3.0 * sd_CL);
      printf("asz = %.5e bsz = %.5e al = %.5e bl = %.5e\n", a_sz, b_sz, a_l, b_l);
      ncm_integrate_2dim (&integ, a_sz, a_l, b_sz, b_l, NCM_DEFAULT_PRECISION, 0.0, &P4, &err4);
    }
    printf("P1 = %.10g err1 = %.10g\n", P1, err1);
    printf("P2 = %.10g err2 = %.10g\n", P2, err2);
    printf("P3 = %.10g err3 = %.10g\n", P3, err3);
    printf("P4 = %.10g err4 = %.10g\n", P4, err4);
    P = P1 + P2 + P3 + P4;
    err = err1 + err2 +err3 + err4;
  }
*/
 
  return P;
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
  data.Mcut =            5.0e14; /* implement Mcut as a property of this object */

  //sd_Pl = data.mobs_params[NC_CLUSTER_MASS_PLCL_SD_PL];
  //sd_CL = data.mobs_params[NC_CLUSTER_MASS_PLCL_SD_CL];
  
  integ.f = _nc_cluster_mass_plcl_Msz_Ml_M500_intp_integrand;
  integ.userdata = &data;

  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);

  gdouble a_sz, a_l, b_sz, b_l;
  a_sz = a_l = log(1.0e10);
  b_sz = b_l = log(1.0e17); 
  ncm_integrate_2dim (&integ, a_sz, a_l, b_sz, b_l, NCM_DEFAULT_PRECISION, 0.0, &P, &err);
    
  return P;
  
}

static gboolean
_nc_cluster_mass_plcl_resample (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnMobs, const gdouble *lnMobs_params, NcmRNG *rng)
{
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (clusterm);
  gdouble r_SZ, r_L;

  const gdouble lnM_SZ = _SZ_lnmass_mean (mszl, lnM);
  const gdouble lnM_L = _Lens_lnmass_mean (mszl, lnM);
  
  NCM_UNUSED (lnMobs_params);
  NCM_UNUSED (clusterm);
  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);
  
  ncm_rng_lock (rng);
  gsl_ran_bivariate_gaussian (rng->r, SD_SZ, SD_L, COR, &r_SZ, &r_L);
  lnMobs[NC_CLUSTER_MASS_PLCL_MPL] = lnM_SZ + r_SZ;
  lnMobs[NC_CLUSTER_MASS_PLCL_MCL] = lnM_L + r_L;
  ncm_rng_unlock (rng);

  return FALSE;
  //return (lnMobs[NC_CLUSTER_MASS_PLCL_MPL] >= mszl->lnMobs_min) && (lnMobs[NC_CLUSTER_MASS_PLCL_MCL] >= mszl->lnMobs_min); 
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

/*
static void
_nc_cluster_mass_plcl_n_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper)
{
  NcClusterMassBenson *msz = NC_CLUSTER_MASS_BENSON (clusterm);

  const gdouble xil = msz->signif_obs_min;
  const gdouble zetal = _significance_to_zeta (clusterm, model, 2.0, xil) - 7.0 * D_SZ;
  const gdouble lnMl = GSL_MAX (_zeta_to_mass (clusterm, model, 2.0, zetal), log (NC_CLUSTER_MASS_BENSON_M_LOWER_BOUND));

  const gdouble xiu = msz->signif_obs_max;
  const gdouble zetau = _significance_to_zeta (clusterm, model, 0.0, xiu) + 7.0 * D_SZ;
  const gdouble lnMu = _zeta_to_mass (clusterm, model, 0.0, zetau);

  *lnM_lower = lnMl;
  *lnM_upper = lnMu;

  return;  
}
*/
