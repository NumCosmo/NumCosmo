/***************************************************************************
 *            nc_cor_cluster_cmb_lens_limber.c
 *
 *  Wed June 11 17:18:58 2014
 *  Copyright  2014  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_cor_cluster_cmb_lens_limber.c
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
 * SECTION:nc_cor_cluster_cmb_lens_limber
 * @title: NcCorClusterCmbLensLimber
 * @short_description: Cluster and CMB lensing correlation using halo model and Limber approximation.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cor_cluster_cmb_lens_limber.h"
#include "math/integral.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_cfg.h"

G_DEFINE_TYPE (NcCorClusterCmbLensLimber, nc_cor_cluster_cmb_lens_limber, G_TYPE_OBJECT);

#define _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY 6

/**
 * nc_cor_cluster_cmb_lens_limber_new:
 *
 * This function allocates memory for a new #NcCorClusterCmbLensLimber object.
 *
 * Returns: A new #NcCorClusterCmbLensLimber.
 */
NcCorClusterCmbLensLimber *
nc_cor_cluster_cmb_lens_limber_new ()
{
  NcCorClusterCmbLensLimber *cccll = g_object_new (NC_TYPE_COR_CLUSTER_CMB_LENS_LIMBER,
                                          NULL);
  return cccll;
}

static void
nc_cor_cluster_cmb_lens_limber_init (NcCorClusterCmbLensLimber *cccll)
{
  cccll->oneh_int_mass_spline = NULL;
}

static void
nc_cor_cluster_cmb_lens_limber_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_cor_cluster_cmb_lens_limber_parent_class)->finalize (object);
}

static void
nc_cor_cluster_cmb_lens_limber_class_init (NcCorClusterCmbLensLimberClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = nc_cor_cluster_cmb_lens_limber_finalize;
}

typedef struct _integrand_data_1h_m
{
  NcCorClusterCmbLensLimber *cccll;
  NcClusterAbundance *cad;
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcHaloBiasFunc *hbf;
  NcHaloDensityProfile *dp;
  NcClusterMass *clusterm;
  gdouble z;
  gdouble lnM;
  gdouble *lnM_obs;
  gdouble *lnM_obs_params;
  gdouble k;
} integrand_data_1h_m;

/*
static gdouble
_integrand_mass_1h (gdouble lnM, gpointer userdata)
{
  integrand_data_1h_m *int_data = (integrand_data_1h_m *) userdata;
  gdouble M = exp(lnM);
  gdouble dn_dlnM = nc_halo_mass_function_dn_dlnM (int_data->cad->mfp, int_data->cosmo, lnM, int_data->z);
  //printf("integrando mass, k = %.5g M = %.5e dndM = %.5g \n", int_data->k, M, dn_dlnM);
  gdouble u = nc_halo_density_profile_eval_fourier (int_data->dp, int_data->cosmo, int_data->k, M, int_data->z);
  //printf("integrando mass, u = %.15g\n", u);
  gdouble integrand_mass_1h = M * dn_dlnM * u * nc_cluster_mass_intp (int_data->cad->m, int_data->cosmo, lnM, int_data->z);
  //printf("integrando mass, lognormal = %.15g\n", nc_cluster_mass_intp (int_data->cad->m, int_data->cosmo, lnM, int_data->z));

  return integrand_mass_1h;
}
*/

static gdouble
_integrand_powspec_1h (gdouble lnM, gpointer userdata)
{
  integrand_data_1h_m *int_data = (integrand_data_1h_m *) userdata;
  gdouble M = exp(lnM);
  gdouble dn_dlnM = nc_halo_mass_function_dn_dlnM (int_data->cad->mfp, int_data->cosmo, lnM, int_data->z);
  gdouble u = 1.0;//nc_halo_density_profile_eval_fourier (int_data->dp, int_data->cosmo, int_data->k, M, int_data->z);
  gdouble rho_mz = nc_hicosmo_E2Omega_m (int_data->cosmo, int_data->z) * ncm_c_crit_mass_density_h2_solar_mass_Mpc3 ();
  //gdouble integrand_powspec_1h = M * dn_dlnM * u; //* nc_cluster_mass_intp (int_data->cad->m, int_data->cosmo, lnM, int_data->z);
  gdouble integrand_powspec_1h = (M / rho_mz) * (M / rho_mz) * dn_dlnM * u * u;
  //printf("M = %.5e, rho = %.5g, n = %.5g, u = %.5g\n", M, rho_mz, dn_dlnM, u);
  //printf("M = %.5e\n", M);
  //printf("integrando mass, lognormal = %.15g\n", nc_cluster_mass_intp (int_data->cad->m, int_data->cosmo, lnM, int_data->z));

  g_assert_not_reached (); /* FIXME: nc_halo_density_profile_eval_fourier */
  
  return integrand_powspec_1h;
}

/**
 * nc_cor_cluster_cmb_lens_limber_oneh_int_mass:
 * @cccll: a #NcCorClusterCmbLensLimber
 * @cad: a #NcClusterAbundance
 * @clusterm: a #NcClusterMass
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @k: mode
 * @z: redshift
 * @lnM_obs: (in) (array) (element-type double): FIXME
 * @lnM_obs_params: (in) (array) (element-type double): FIXME
 *
 * This function computes the 1-halo integral on mass of the cluster and CMB lensing potential $\psi$, using Limber approximation.
 *
 * Returns: FIXME
 */
gdouble
nc_cor_cluster_cmb_lens_limber_oneh_int_mass (NcCorClusterCmbLensLimber *cccll, NcClusterAbundance *cad, NcClusterMass *clusterm, NcHICosmo *cosmo, NcHaloDensityProfile *dp, gdouble k, gdouble z, gdouble *lnM_obs, gdouble *lnM_obs_params)
{
  integrand_data_1h_m int_data;
  gdouble int_mass_1h, err, lnMl, lnMu;
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  int_data.cad            = cad;
  int_data.clusterm       = clusterm;
  int_data.cosmo          = cosmo;
  int_data.z              = z;
  int_data.lnM_obs        = lnM_obs;
  int_data.lnM_obs_params = lnM_obs_params;
  int_data.dp             = dp;
  int_data.k              = k;

  //F.function = &_integrand_mass_1h;
  F.function = &_integrand_powspec_1h;
  F.params = &int_data;

  nc_cluster_mass_n_limits (clusterm, cosmo, &lnMl, &lnMu);

  //printf("z = %.5g lnMl = %.5g lnMu = %.5g lnMobs0 = %.5g lnMobs1 = %.5g\n", z, lnMl, lnMu, lnM_obs_params[0], lnM_obs_params[1]);
  //gsl_integration_qag (&F, lnMl, lnMu, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &int_mass_1h, &err);
  gsl_integration_qag (&F, -50.0, 50.0, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &int_mass_1h, &err);

  ncm_memory_pool_return (w);

  return int_mass_1h;
}

typedef struct _integrand_data_1h_z
{
  NcCorClusterCmbLensLimber *cccll;
  NcClusterAbundance *cad;
  NcClusterMass *clusterm;
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcHaloBiasFunc *hbf;
  NcHaloDensityProfile *dp;
  gint l;
  gdouble z;
  gdouble lnM;
  gdouble *lnM_obs;
  gdouble *lnM_obs_params;
  gdouble k;
  gdouble dc_z;
  gdouble dc_zdec; /* comoving distance at z decoupling */
} integrand_data_1h_z;

static gdouble
_integrand_redshift_1h (gdouble z, gpointer userdata)
{
  integrand_data_1h_z *int_data = (integrand_data_1h_z *) userdata;
  int_data->dc_z = nc_distance_comoving (int_data->dist, int_data->cosmo, z);
  //int_data->k = int_data->l / int_data->dc_z; /* dimensionless */
  int_data->k = int_data->l / (int_data->dc_z * ncm_c_hubble_radius_hm1_Mpc ()); /* in units of h Mpc^-1 */
  //printf("l = %4.d k = %.5g kdim = %.5g\n", int_data->l, int_data->k, int_data->k * 1.0e5 / ncm_c_c());
  gdouble dcdec_m_dc = int_data->dc_zdec - int_data->dc_z;
  gdouble ds = int_data->dc_z * dcdec_m_dc / int_data->dc_zdec;
  gdouble rho_mz = nc_hicosmo_E2Omega_m (int_data->cosmo, int_data->z) * ncm_c_crit_mass_density_h2_solar_mass_Mpc3 ();
  gdouble integral_mass = nc_cor_cluster_cmb_lens_limber_oneh_int_mass (int_data->cccll, int_data->cad, int_data->clusterm, int_data->cosmo, int_data->dp, int_data->k, z, int_data->lnM_obs, int_data->lnM_obs_params);
  //printf("int_mass= %.5g\n", integral_mass);
  gdouble integrand_z_1h = (1.0 + z) * ds / (nc_hicosmo_E (int_data->cosmo, z) * rho_mz) * integral_mass;

  return integrand_z_1h;
}

/**
 * nc_cor_cluster_cmb_lens_limber_oneh_term:
 * @cccll: a #NcCorClusterCmbLensLimber
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @dist: a #NcDistance
 * @dp: a #NcHaloDensityProfile
 * @l: spherical harmonis index
 * @lnM_obs: (in) (array) (element-type double): FIXME
 * @lnM_obs_params: (in) (array) (element-type double): FIXME
 * @z_obs: (in) (array) (element-type double): FIXME
 * @z_obs_params: (in) (array) (element-type double): FIXME
 *
 * This function computes the 1-halo term of the cluster and CMB lensing potential $\psi$, using Limber approximation:
 * \begin{equation}
 * \left(C_l^{cl \psi} \right)_{1h} = -3 \frac{\Omega_{m0} H_0^2}{l^2} \int dz \frac{c}{H(z)} (1+z) \chi(z)^2
 *  \frac{\left( \chi(z_\ast) - \chi(z) \right)}{\chi(z_\ast) \chi(z)} \int d\ln M S(\ln M, z) \frac{M}{\overline{\rho}(z)} \frac{dn(M, z)}{d\ln M}
 * \tilde{u}^\ast (k = l/\chi(z) \vert M),
 * \end{equation}
 * where ...
 *
 * Returns: $\left(C_l^{cl \psi} \right)_{1h}$
 */
gdouble
nc_cor_cluster_cmb_lens_limber_oneh_term (NcCorClusterCmbLensLimber *cccll, NcClusterAbundance *cad, NcHICosmo *cosmo, NcDistance *dist, NcHaloDensityProfile *dp, gint l, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble *z_obs, gdouble *z_obs_params)
{
  gdouble cor_1h, ll, cons_factor, zl, zu, err;
  integrand_data_1h_z int_data;
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  int_data.cccll = cccll;
  int_data.cad   = cad;
  int_data.cosmo = cosmo;
  int_data.dist  = dist;
  int_data.dp    = dp;
  int_data.l     = l;
  int_data.lnM_obs = lnM_obs;
  int_data.lnM_obs_params = lnM_obs_params;

  F.function = &_integrand_redshift_1h;
  F.params = &int_data;

  zl = z_obs_params[0];
  zu = z_obs_params[1];
  //printf("zl = %.5g zu = %.5g\n", zl, zu);

  ll = l * l;
  cons_factor      = - 3.0 * nc_hicosmo_Omega_m0 (cosmo) / ll;
  int_data.dc_zdec = nc_distance_comoving_lss (dist, cosmo);

  gsl_integration_qag (&F, zl, zu, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &cor_1h, &err);

//  ncm_spline_clear (&cccll->oneh_int_mass);
//  cccll->oneh_int_mass = ncm_spline_cubic_notaknot_new ();
//  ncm_spline_set_func (cccll->oneh_int_mass, NCM_SPLINE_FUNCTION_SPLINE, &F, zl, zu, 1000000, 1e-7);

  ncm_memory_pool_return (w);

//  return ncm_spline_eval_integ (cccll->oneh_int_mass, zl, zu) * cons_factor;

  return cor_1h * cons_factor;
}

typedef struct _integrand_data_2h_mass1
{
  NcCorClusterCmbLensLimber *cccll;
  NcClusterAbundance *cad;
  NcClusterMass *clusterm;
  NcHICosmo *cosmo;
  NcHaloBiasFunc *hbf;
  gdouble *lnMobs_params;
  gdouble z;
} integrand_data_2h_mass1;

static gdouble
_integrand_mass_2h_first (gdouble lnM, gpointer userdata)
{
  integrand_data_2h_mass1 *int_data = (integrand_data_2h_mass1 *) userdata;
  gdouble dn_dlnM_times_b = nc_halo_bias_func_integrand (int_data->hbf, int_data->cosmo, lnM, int_data->z);
  gdouble integrand_mass_2h_first = dn_dlnM_times_b * nc_cluster_mass_intp (int_data->clusterm, int_data->cosmo, lnM, int_data->z);

  return integrand_mass_2h_first;
}

/**
 * nc_cor_cluster_cmb_lens_limber_twoh_int_mass1:
 * @cccll: a #NcCorClusterCmbLensLimber
 * @cad: a #NcClusterAbundance
 * @clusterm: a #NcClusterMass
 * @cosmo: a #NcHICosmo
 * @z: redshift
 *
 * This function computes the first integral on mass of the 2-halo term of the cluster and CMB lensing potential
 * $\psi$, using Limber approximation.
 *
 * Returns: FIXME
 */
gdouble
nc_cor_cluster_cmb_lens_limber_twoh_int_mass1 (NcCorClusterCmbLensLimber *cccll, NcClusterAbundance *cad, NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble z)
{
  integrand_data_2h_mass1 int_data;
  gdouble int_mass_2h_first, err, lnMl, lnMu;
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  int_data.cccll          = cccll;
  int_data.cad            = cad;
  int_data.clusterm       = clusterm;
  int_data.hbf            = cad->mbiasf;
  int_data.cosmo          = cosmo;
  int_data.z              = z;

  F.function = &_integrand_mass_2h_first;
  F.params = &int_data;

  nc_cluster_mass_n_limits (clusterm, cosmo, &lnMl, &lnMu);

  //printf("z = %.5g lnMl = %.5g lnMu = %.5g\n", z, lnMl, lnMu);

  gsl_integration_qag (&F, lnMl, lnMu, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &int_mass_2h_first, &err);

  ncm_memory_pool_return (w);

  return int_mass_2h_first;
}

typedef struct _integrand_data_2h_mass2
{
  NcCorClusterCmbLensLimber *cccll;
  NcClusterAbundance *cad;
  NcClusterMass *clusterm;
  NcHICosmo *cosmo;
  NcHaloBiasFunc *hbf;
  NcHaloDensityProfile *dp;
  gdouble z;
  gdouble k;
} integrand_data_2h_mass2;


static gdouble
_integrand_powspec_2h (gdouble lnM, gpointer userdata)
{
  integrand_data_2h_mass2 *int_data = (integrand_data_2h_mass2 *) userdata;

  const gdouble lnR                   = nc_halo_mass_function_lnM_to_lnR (int_data->hbf->mfp, int_data->cosmo, lnM);
  const gdouble V                     = ncm_powspec_filter_volume_rm3 (nc_halo_mass_function_peek_psf (int_data->hbf->mfp)) * exp (3.0 * lnR);
  const gdouble dn_dlnM_times_b       = nc_halo_bias_func_integrand (int_data->hbf, int_data->cosmo, lnM, int_data->z);
  const  gdouble integrand_powspec_2h = V * dn_dlnM_times_b;

  return integrand_powspec_2h;
}

/**
 * nc_cor_cluster_cmb_lens_limber_twoh_int_mm:
 * @cccll: a #NcCorClusterCmbLensLimber
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @k: mode
 * @z: redshift
 *
 * This function computes the 2-halo term of the matter-matter power spectrum.
 *
 * Returns: FIXME
 */
gdouble
nc_cor_cluster_cmb_lens_limber_twoh_int_mm (NcCorClusterCmbLensLimber *cccll, NcClusterAbundance *cad, NcHICosmo *cosmo, NcHaloDensityProfile *dp, gdouble k, gdouble z)
{
  integrand_data_2h_mass2 int_data;
  gdouble ps_2h_mm, err;
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  gdouble int_powspec_mm_2h = 0.0;
  gdouble a = 0.0, b = 0.0;
  gboolean conv1 = FALSE;
  gboolean conv2 = FALSE;
  const gdouble step = 2.0;

  nc_halo_mass_function_prepare_if_needed (int_data.hbf->mfp, cosmo);

  int_data.cccll          = cccll;
  int_data.cad            = cad;
  int_data.hbf            = cad->mbiasf;
  int_data.cosmo          = cosmo;
  int_data.dp             = dp;
  int_data.k              = k;
  int_data.z              = z;

  F.function = &_integrand_powspec_2h;
  F.params = &int_data;

  do {
    gdouble res1 = 0.0;
    gdouble res2 = 0.0;

    if (!conv1)
      gsl_integration_qag (&F, a, a + step, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &res1, &err);
    if (!conv2)
      gsl_integration_qag (&F, b - step, b, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &res2, &err);

    int_powspec_mm_2h += res1 + res2;

    if (!conv1 && fabs (res1 / int_powspec_mm_2h) < NCM_DEFAULT_PRECISION)
      conv1 = TRUE;
    if (!conv2 && fabs (res2 / int_powspec_mm_2h) < NCM_DEFAULT_PRECISION)
      conv2 = TRUE;

    a += step;
    b -= step;

  } while (!conv1 || !conv2);

  ncm_memory_pool_return (w);

  ps_2h_mm = int_powspec_mm_2h * int_powspec_mm_2h * ncm_powspec_eval (nc_halo_mass_function_peek_psf (cad->mbiasf->mfp)->ps, NCM_MODEL (cosmo), z, k);

  return ps_2h_mm;
}

static gdouble
_integrand_mass_2h_second (gdouble lnM, gpointer userdata)
{
  integrand_data_2h_mass2 *int_data = (integrand_data_2h_mass2 *) userdata;
  gdouble M = exp(lnM);
  gdouble dn_dlnM_times_b = nc_halo_bias_func_integrand (int_data->hbf, int_data->cosmo, lnM, int_data->z);
  gdouble u = 1.0; //nc_halo_density_profile_eval_fourier (int_data->dp, int_data->cosmo, int_data->k, M, int_data->z);
  gdouble integrand_mass_2h_second = M * dn_dlnM_times_b * u;

  g_assert_not_reached (); /* FIXME: nc_halo_density_profile_eval_fourier */
  
  return integrand_mass_2h_second;
}

/**
 * nc_cor_cluster_cmb_lens_limber_twoh_int_mass2:
 * @cccll: a #NcCorClusterCmbLensLimber
 * @cad: a #NcClusterAbundance
 * @clusterm: a #NcClusterMass
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @k: mode
 * @z: redshift
 *
 * This function computes the second integral on mass of the 2-halo term of the cluster and CMB lensing potential
 * $\psi$, using Limber approximation.
 *
 * Returns: FIXME
 */
gdouble
nc_cor_cluster_cmb_lens_limber_twoh_int_mass2 (NcCorClusterCmbLensLimber *cccll, NcClusterAbundance *cad, NcClusterMass *clusterm, NcHICosmo *cosmo, NcHaloDensityProfile *dp, gdouble k, gdouble z)
{
  integrand_data_2h_mass2 int_data;
  gdouble int_mass_2h_second, err, lnMl, lnMu;
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  int_data.cccll          = cccll;
  int_data.cad            = cad;
  int_data.clusterm       = clusterm;
  int_data.hbf            = cad->mbiasf;
  int_data.cosmo          = cosmo;
  int_data.dp             = dp;
  int_data.k              = k;
  int_data.z              = z;

  F.function = &_integrand_mass_2h_second;
  F.params = &int_data;

  nc_cluster_mass_n_limits (clusterm, cosmo, &lnMl, &lnMu);

  //gsl_integration_qag (&F, lnMl, lnMu, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &int_mass_2h_second, &err);
  gsl_integration_qag (&F, -50.0, 50.0, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &int_mass_2h_second, &err);

  ncm_memory_pool_return (w);

  return int_mass_2h_second;
}

typedef struct _integrand_data_2hz
{
  NcCorClusterCmbLensLimber *cccll;
  NcClusterAbundance *cad;
  NcClusterMass *clusterm;
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcHaloBiasFunc *hbf;
  NcHaloDensityProfile *dp;
  gint l;
  gdouble dc_zdec; /* comoving distance at z decoupling */
} integrand_data_2hz;

static gdouble
_integrand_redshift_2h (gdouble z, gpointer userdata)
{
  integrand_data_2hz *int_data = (integrand_data_2hz *) userdata;
  const gdouble dc_z       = nc_distance_comoving (int_data->dist, int_data->cosmo, z);
  const gdouble dcdec_m_dc = int_data->dc_zdec - dc_z;
  const gdouble ds         = dc_z * dcdec_m_dc / int_data->dc_zdec;

  const gdouble k         = int_data->l / (dc_z * ncm_c_hubble_radius_hm1_Mpc ()); /* in units of h Mpc^-1 */
  const gdouble rho_mz    = nc_hicosmo_E2Omega_m (int_data->cosmo, z) * ncm_c_crit_mass_density_h2_solar_mass_Mpc3 ();
  const gdouble ps_linear = ncm_powspec_eval (nc_halo_mass_function_peek_psf (int_data->cad->mbiasf->mfp)->ps, NCM_MODEL (int_data->cosmo), z, k);

  const gdouble integral_mass1 = nc_cor_cluster_cmb_lens_limber_twoh_int_mass1 (int_data->cccll, int_data->cad, int_data->clusterm, int_data->cosmo, z);
  const gdouble integral_mass2 = nc_cor_cluster_cmb_lens_limber_twoh_int_mass2 (int_data->cccll, int_data->cad, int_data->clusterm, int_data->cosmo, int_data->dp, k, z);
  const gdouble integrand_z_2h = (1.0 + z) * ds / (nc_hicosmo_E (int_data->cosmo, z) * rho_mz) * ps_linear
    * integral_mass1 * integral_mass2;

  return integrand_z_2h;
}

/**
 * nc_cor_cluster_cmb_lens_limber_twoh_term:
 * @cccll: a #NcCorClusterCmbLensLimber
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @dist: a #NcDistance
 * @dp: a #NcHaloDensityProfile
 * @l: spherical harmonis index
 * @z_obs: (in) (array) (element-type double): FIXME
 * @z_obs_params: (in) (array) (element-type double): FIXME
 *
 * This function computes the 2-halo term of the cluster and CMB lensing potential $\psi$, using Limber approximation:
 * \begin{equation}
 * \left(C_l^{cl \psi} \right)_{1h} = -3 \frac{\Omega_{m0} H_0^2}{l^2 c^2} \int dz \frac{c}{H(z)} (1+z) \chi(z)^2
 *  \frac{\left( \chi(z_\ast) - \chi(z) \right)}{\chi(z_\ast) \chi(z)} P_{lin}(k = l/\chi, z) \int d\ln M S(\ln M, z) \frac{dn(M, z)}{d\ln M} b(M, z)
 * \int dln M^{\prime} \frac{M^\prime}{\overline{\rho}(z)} \frac{dn(M^\prime, z)}{d\ln M} b(M^\prime, z) \tilde{u}^\ast (k = l/\chi(z) \vert M^\prime),
 * \end{equation}
 * where ...
 *
 * Returns: $\left(C_l^{cl \psi} \right)_{2h}$
 */
gdouble
nc_cor_cluster_cmb_lens_limber_twoh_term (NcCorClusterCmbLensLimber *cccll, NcClusterAbundance *cad, NcHICosmo *cosmo, NcDistance *dist, NcHaloDensityProfile *dp, gint l, gdouble *z_obs, gdouble *z_obs_params)
{
  gdouble cor_2h, ll, cons_factor, zl, zu, err;
  integrand_data_2hz int_data;
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  int_data.cccll = cccll;
  int_data.cad   = cad;
  int_data.cosmo = cosmo;
  int_data.dist  = dist;
  int_data.hbf   = cad->mbiasf;
  int_data.dp    = dp;
  int_data.l     = l;

  nc_halo_mass_function_prepare_if_needed (int_data.hbf->mfp, cosmo);

  F.function = &_integrand_redshift_2h;
  F.params   = &int_data;

  zl = z_obs_params[0];
  zu = z_obs_params[1];

  ll               = l * l;
  cons_factor      = - 3.0 * nc_hicosmo_Omega_m0 (cosmo) / ll;
  int_data.dc_zdec = nc_distance_comoving_lss (dist, cosmo);

  gsl_integration_qag (&F, zl, zu, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &cor_2h, &err);

  ncm_memory_pool_return (w);

  return cor_2h * cons_factor;
}
