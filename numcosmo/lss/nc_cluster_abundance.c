/***************************************************************************
 *            nc_cluster_abundance.c
 *
 *  Tue Apr 20 10:59:01 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com> & <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_cluster_abundance
 * @title: NcClusterAbundance
 * @short_description: Cluster abundance distribution
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_abundance.h"
#include "lss/nc_cluster_mass_benson.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline2d_bicubic.h"
#include "math/integral.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_histogram.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_MASS_FUNCTION,
  PROP_MEANBIAS,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcClusterAbundance, nc_cluster_abundance, G_TYPE_OBJECT);

#define INTEG_D2NDZDLNM_NNODES (200)
#define LNM_MIN (10.0 * M_LN10)
#define _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY 6

static gdouble
_intp_d2N (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble lnM, gdouble z)
{
  NCM_UNUSED (cad);
  NCM_UNUSED (cosmo);
  NCM_UNUSED (lnM);
  NCM_UNUSED (z);
  g_error ("Function d2NdzdlnM_val not implemented or cad not prepared.");
  
  return 0.0;
};

static gdouble
_N (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm)
{
  NCM_UNUSED (cad);
  NCM_UNUSED (cosmo);
  g_error ("Function N_val not implemented or cad not prepared.");
  
  return 0.0;
};

static void
nc_cluster_abundance_init (NcClusterAbundance *cad)
{
  NcmVector *u1_knots          = ncm_vector_new (INTEG_D2NDZDLNM_NNODES);
  NcmVector *u2_knots          = ncm_vector_new (INTEG_D2NDZDLNM_NNODES);
  NcmVector *u3_knots          = ncm_vector_new (INTEG_D2NDZDLNM_NNODES);
  NcmVector *integ_lnM_knots   = ncm_vector_new (INTEG_D2NDZDLNM_NNODES);
  NcmVector *integ_z_knots     = ncm_vector_new (INTEG_D2NDZDLNM_NNODES);
  NcmMatrix *integ_lnM_z_knots = ncm_matrix_new (INTEG_D2NDZDLNM_NNODES, INTEG_D2NDZDLNM_NNODES);
  
  cad->N        = &_N;
  cad->intp_d2N = &_intp_d2N;
  
  cad->purity    = NULL;
  cad->sd_lnM    = NULL;
  cad->norma     = 0.0;
  cad->log_norma = 0.0;
  
#define D2NDZDLNM_Z(cad) ((cad)->d2NdzdlnM->yv)
#define D2NDZDLNM_LNM(cad) ((cad)->d2NdzdlnM->xv)
#define D2NDZDLNM_VAL(cad) ((cad)->d2NdzdlnM->zm)
  
#define DBDLNM_Z(cad) ((cad)->dbdlnM->yv)
#define DBDLNM_LNM(cad) ((cad)->dbdlnM->xv)
#define DBDLNM_VAL(cad) ((cad)->dbdlnM->zm)
  
  cad->dbdlnM = NULL;
  
  cad->inv_z   = ncm_spline_cubic_notaknot_new ();
  cad->inv_lnM = ncm_spline_cubic_notaknot_new ();
  
  ncm_spline_set (cad->inv_z, u1_knots, integ_z_knots, FALSE);
  ncm_spline_set (cad->inv_lnM, u2_knots, integ_lnM_knots, FALSE);
  
  cad->inv_lnM_z = ncm_spline2d_bicubic_notaknot_new ();
  ncm_spline2d_set (cad->inv_lnM_z, u3_knots, integ_z_knots, integ_lnM_z_knots, FALSE);
  
  ncm_vector_free (u1_knots);
  ncm_vector_free (u2_knots);
  ncm_vector_free (u3_knots);
  ncm_vector_free (integ_lnM_knots);
  ncm_vector_free (integ_z_knots);
  ncm_matrix_free (integ_lnM_z_knots);
  
  cad->ctrl_cosmo = ncm_model_ctrl_new (NULL);
  cad->ctrl_reion = ncm_model_ctrl_new (NULL);
  cad->ctrl_z     = ncm_model_ctrl_new (NULL);
  cad->ctrl_m     = ncm_model_ctrl_new (NULL);
}

static void
_nc_cluster_abundance_dispose (GObject *object)
{
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (object);
  
  nc_halo_mass_function_clear (&cad->mfp);
  nc_halo_bias_func_clear (&cad->mbiasf);
  
  ncm_spline_clear (&cad->inv_z);
  ncm_spline_clear (&cad->inv_lnM);
  
  ncm_spline2d_clear (&cad->inv_lnM_z);
  
  ncm_model_ctrl_clear (&cad->ctrl_cosmo);
  ncm_model_ctrl_clear (&cad->ctrl_reion);
  ncm_model_ctrl_clear (&cad->ctrl_z);
  ncm_model_ctrl_clear (&cad->ctrl_m);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_abundance_parent_class)->dispose (object);
}

static void
_nc_cluster_abundance_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_abundance_parent_class)->finalize (object);
}

static void
_nc_cluster_abundance_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (object);
  
  g_return_if_fail (NC_IS_CLUSTER_ABUNDANCE (object));
  
  switch (prop_id)
  {
    case PROP_MASS_FUNCTION:
      cad->mfp = g_value_dup_object (value);
      break;
    case PROP_MEANBIAS:
      cad->mbiasf = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_abundance_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (object);
  
  g_return_if_fail (NC_IS_CLUSTER_ABUNDANCE (object));
  
  switch (prop_id)
  {
    case PROP_MASS_FUNCTION:
      g_value_set_object (value, cad->mfp);
      break;
    case PROP_MEANBIAS:
      g_value_set_object (value, cad->mbiasf);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cluster_abundance_class_init (NcClusterAbundanceClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  
  object_class->set_property = &_nc_cluster_abundance_set_property;
  object_class->get_property = &_nc_cluster_abundance_get_property;
  object_class->dispose      = &_nc_cluster_abundance_dispose;
  object_class->finalize     = &_nc_cluster_abundance_finalize;
  
  /**
   * NcClusterAbundance:halo-mass-function:
   *
   * Number  density of halos with mass [M, M + dM].
   */
  g_object_class_install_property (object_class,
                                   PROP_MASS_FUNCTION,
                                   g_param_spec_object ("halo-mass-function",
                                                        NULL,
                                                        "Mass Function",
                                                        NC_TYPE_HALO_MASS_FUNCTION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcClusterAbundance:mean-bias:
   *
   * FIXME
   */
  g_object_class_install_property (object_class,
                                   PROP_MEANBIAS,
                                   g_param_spec_object ("mean-bias",
                                                        NULL,
                                                        "Mean Halo Bias Function",
                                                        NC_TYPE_HALO_BIAS_FUNC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_cluster_abundance_new:
 * @mfp: a #NcHaloMassFunction
 * @mbiasf: (allow-none): a #NcHaloBiasFunc
 *
 * This function allocates memory for a new #NcClusterAbundance object and sets its properties to the values from
 * the input arguments.
 *
 * Returns: A new #NcClusterAbundance.
 */
NcClusterAbundance *
nc_cluster_abundance_new (NcHaloMassFunction *mfp, NcHaloBiasFunc *mbiasf)
{
  NcClusterAbundance *cad = g_object_new (NC_TYPE_CLUSTER_ABUNDANCE,
                                          "halo-mass-function", mfp,
                                          "mean-bias",          mbiasf,
                                          NULL);
  
  return cad;
}

/**
 * nc_cluster_abundance_nodist_new:
 * @mfp: a #NcHaloMassFunction.
 * @mbiasf: (allow-none): a #NcHaloBiasFunc.
 *
 * This function allocates memory for a new #NcClusterAbundance object and sets its properties to the values from
 * the input arguments.
 *
 * Returns: A new #NcClusterAbundance.
 */
NcClusterAbundance *
nc_cluster_abundance_nodist_new (NcHaloMassFunction *mfp, NcHaloBiasFunc *mbiasf)
{
  NcClusterAbundance *cad = g_object_new (NC_TYPE_CLUSTER_ABUNDANCE,
                                          "halo-mass-function", mfp,
                                          "mean-bias",          mbiasf,
                                          NULL);
  
  return cad;
}

/**
 * nc_cluster_abundance_ref:
 * @cad: a #NcClusterAbundance
 *
 * Increases the reference count of @cad by one.
 *
 * Returns: (transfer full): @cad
 */
NcClusterAbundance *
nc_cluster_abundance_ref (NcClusterAbundance *cad)
{
  return g_object_ref (cad);
}

/**
 * nc_cluster_abundance_free:
 * @cad: a #NcClusterAbundance
 *
 * Atomically decrements the reference count of @cad by one. If the reference count drops to 0,
 * all memory allocated by @cad is released.
 *
 */
void
nc_cluster_abundance_free (NcClusterAbundance *cad)
{
  g_object_unref (cad);
}

/**
 * nc_cluster_abundance_clear:
 * @cad: a #NcClusterAbundance
 *
 * Atomically decrements the reference count of @cad by one. If the reference count drops to 0,
 * all memory allocated by @cad is released.
 *
 */
void
nc_cluster_abundance_clear (NcClusterAbundance **cad)
{
  g_clear_object (cad);
}

typedef struct _NcClusterAbundanceInt
{
  NcClusterAbundance *cad;
  NcHICosmo *cosmo;
  NcClusterRedshift *clusterz;
  NcClusterMass *clusterm;
  gdouble z;
  gdouble lnM;
  const gdouble *z_obs;
  const gdouble *z_obs_lower;
  const gdouble *z_obs_upper;
  const gdouble *z_obs_params;
  const gdouble *lnM_obs;
  const gdouble *lnM_obs_lower;
  const gdouble *lnM_obs_upper;
  const gdouble *lnM_obs_params;
} NcClusterAbundanceInt;

static gdouble
_nc_cluster_abundance_z_p_lnM_p_d2n_integrand (gdouble lnM, gdouble z, gpointer userdata)
{
  NcClusterAbundanceInt *obs_data = (NcClusterAbundanceInt *) userdata;
  NcClusterAbundance *cad         = obs_data->cad;
  const gdouble p_z_zr            = nc_cluster_redshift_p (obs_data->clusterz, obs_data->cosmo, lnM, z, obs_data->z_obs, obs_data->z_obs_params);
  const gdouble p_M_Mobs          = nc_cluster_mass_p (obs_data->clusterm, obs_data->cosmo, lnM, z, obs_data->lnM_obs, obs_data->lnM_obs_params);
  const gdouble d2NdzdlnM         = nc_halo_mass_function_d2n_dzdlnM (cad->mfp, obs_data->cosmo, lnM, z);
  
  return p_z_zr * p_M_Mobs * d2NdzdlnM;
}

/**
 * nc_cluster_abundance_z_p_lnM_p_d2n:
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @clusterz: a #NcClusterRedshift
 * @clusterm: a #NcClusterMass
 * @lnM_obs: (array) (element-type gdouble): logarithm base e of the observed mass
 * @lnM_obs_params: (array) (element-type gdouble): FIXME
 * @z_obs: (array) (element-type gdouble): observed redshift
 * @z_obs_params: (array) (element-type gdouble): FIXME
 *
 * This function computes $ \int_0^\infty dz \int_{-\infty}^\infty d\ln M \frac{d^2N(\ln M, z)}{dzd\ln M} * P(z^{phot}|z) *
 * P(\ln M^{obs}|\ln M, z) $. We studied the convergence of this integral to optimize this function. We verified
 * that it converges to 5 decimal places at the redshift interval $ [z^{phot} - 10\sigma^{phot}, z^{phot} +
 * 10\sigma^{phot}] $ and the mass interval $ [\ln M^{obs} - 7\sigma_{\ln M}, \ln M^{obs} + 7\sigma_{\ln M}] $.
 *
 * Returns: a gdouble which represents $ \frac{d^2N(\ln M^{obs}, z^{phot})}{dzd\ln M} $.
 */
gdouble
nc_cluster_abundance_z_p_lnM_p_d2n (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble *z_obs, gdouble *z_obs_params)
{
  gdouble d2N, zl, zu, lnMl, lnMu, err;
  NcClusterAbundanceInt obs_data;
  NcmIntegrand2dim integ;
  
  obs_data.cad            = cad;
  obs_data.cosmo          = cosmo;
  obs_data.clusterz       = clusterz;
  obs_data.clusterm       = clusterm;
  obs_data.lnM_obs        = lnM_obs;
  obs_data.lnM_obs_params = lnM_obs_params;
  obs_data.z_obs          = z_obs;
  obs_data.z_obs_params   = z_obs_params;
  
  integ.f        = &_nc_cluster_abundance_z_p_lnM_p_d2n_integrand;
  integ.userdata = &obs_data;
  
  nc_cluster_redshift_p_limits (clusterz, cosmo, z_obs, z_obs_params, &zl, &zu);
  nc_cluster_mass_p_limits (clusterm, cosmo, lnM_obs, lnM_obs_params, &lnMl, &lnMu);
  
  ncm_integrate_2dim (&integ, lnMl, zl, lnMu, zu, NCM_DEFAULT_PRECISION, 0.0, &d2N, &err);
  
  return d2N;
}

static gdouble
_nc_cluster_abundance_z_p_d2n_integrand (gdouble z, gpointer params)
{
  NcClusterAbundanceInt *obs_data = (NcClusterAbundanceInt *) params;
  NcClusterAbundance *cad         = obs_data->cad;
  const gdouble p_z_zr            = nc_cluster_redshift_p (obs_data->clusterz, obs_data->cosmo, obs_data->lnM, z, obs_data->z_obs, obs_data->z_obs_params);
  const gdouble d2NdzdlnM         = nc_halo_mass_function_d2n_dzdlnM (cad->mfp, obs_data->cosmo, obs_data->lnM, z);
  
  return p_z_zr * d2NdzdlnM;
}

/**
 * nc_cluster_abundance_z_p_d2n:
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @clusterz: a #NcClusterRedshift
 * @clusterm: a #NcClusterMass
 * @lnM: the logarithm base e of the mass (gdouble)
 * @z_obs: (array) (element-type gdouble): observed redshift
 * @z_obs_params: (array) (element-type gdouble): FIXME
 *
 * This function computes $ \int_{z_{phot} - 10\sigma_{phot}}^{z_{phot} + 10\sigma_{phot}} dz \,
 * \frac{d^2N}{dzdlnM} * P(z^{photo}|z) $. The integral limits were determined requiring a precision
 * to five decimal places.
 *
 * Returns: a gdouble which corresponds to $ \int_{z_{phot} - 10\sigma_{phot}}^{z_{phot} + 10\sigma_{phot}} dz \,
 * \frac{d^2N}{dzdlnM} * P(z^{photo}|z) $.
 */
gdouble
nc_cluster_abundance_z_p_d2n (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble lnM, gdouble *z_obs, gdouble *z_obs_params)
{
  NcClusterAbundanceInt obs_data;
  gdouble d2N, zl, zu, err;
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  
  obs_data.cad          = cad;
  obs_data.cosmo        = cosmo;
  obs_data.clusterz     = clusterz;
  obs_data.clusterm     = clusterm;
  obs_data.lnM          = lnM;
  obs_data.z_obs        = z_obs;
  obs_data.z_obs_params = z_obs_params;
  
  F.function = &_nc_cluster_abundance_z_p_d2n_integrand;
  F.params   = &obs_data;
  
  nc_cluster_redshift_p_limits (clusterz, cosmo, z_obs, z_obs_params, &zl, &zu);
  
  gsl_integration_qag (&F, zl, zu, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &d2N, &err);
  
  ncm_memory_pool_return (w);
  
  return d2N;
}

/**
 * _nc_cluster_abundance_lnM_p_d2n_integrand:
 * @lnM: logarithm base e of the mass
 * @params: a gpointer
 *
 * This function computes \f$ \frac{d^2N}{dzdlnM} \times P(\ln M^{obs}| \ln M) \f$
 * where \f$ P(M^{obs}| M) \f$ is the
 * probability distribution of the observed mass \f$ M^{obs} \f$ (mass-observable relation).
 *
 * Returns: a gdouble which is the integrand value at lnM.
 */
static gdouble
_nc_cluster_abundance_lnM_p_d2n_integrand (gdouble lnM, gpointer params)
{
  NcClusterAbundanceInt *obs_data = (NcClusterAbundanceInt *) params;
  NcClusterAbundance *cad         = obs_data->cad;
  const gdouble p_M_Mobs          = nc_cluster_mass_p (obs_data->clusterm, obs_data->cosmo, lnM, obs_data->z, obs_data->lnM_obs, obs_data->lnM_obs_params);
  const gdouble d2NdzdlnM         = nc_halo_mass_function_d2n_dzdlnM (cad->mfp, obs_data->cosmo, lnM, obs_data->z);
  
  /*printf ("POINT: % 22.15g % 22.15g % 22.15g\n", lnM, p_M_Mobs, d2NdzdlnM);*/

  return p_M_Mobs * d2NdzdlnM;
}

/**
 * nc_cluster_abundance_lnM_p_d2n:
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @clusterz: a #NcClusterRedshift
 * @clusterm: a #NcClusterMass
 * @lnM_obs: (array) (element-type gdouble): logarithm base e of the observed mass
 * @lnM_obs_params: (array) (element-type gdouble) (allow-none): other information of the observed mass, such as its error
 * @z: redshift
 *
 * This function computes $ \int_{\ln M^{obs} - 7\sigma_{\ln M}}^{\ln M^{obs} + 7\sigma_{\ln M}} d\ln M \,
 * \frac{d^2N}{dzdlnM} * P(\ln M^{obs}|\ln M) $. The integral limits were determined requiring a precision
 * to five decimal places.
 *
 * Returns: a gdouble which corresponds to $ \int_{\ln M^{obs} - 7\sigma_{\ln M}}^{\ln M^{obs} + 7\sigma_{\ln M}} d\ln M \,
 * \frac{d^2N}{dzdlnM} * P(\ln M^{obs}|\ln M) $.
 */
gdouble
nc_cluster_abundance_lnM_p_d2n (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble z)
{
  NcClusterAbundanceInt obs_data;
  gdouble d2N, lnMl, lnMu, err;
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  
  obs_data.cad            = cad;
  obs_data.cosmo          = cosmo;
  obs_data.clusterz       = clusterz;
  obs_data.clusterm       = clusterm;
  obs_data.lnM_obs        = lnM_obs;
  obs_data.lnM_obs_params = lnM_obs_params;
  obs_data.z              = z;
  
  F.function = &_nc_cluster_abundance_lnM_p_d2n_integrand;
  F.params   = &obs_data;
  
  nc_cluster_mass_p_limits (clusterm, cosmo, lnM_obs, lnM_obs_params, &lnMl, &lnMu);
  
  gsl_integration_qag (&F, lnMl, lnMu, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &d2N, &err);
  ncm_memory_pool_return (w);
  
  return d2N;
}

/**
 * nc_cluster_abundance_d2n:
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @clusterz: a #NcClusterRedshift
 * @clusterm: a #NcClusterMass
 * @lnM: true mass
 * @z: true redshift
 *
 * This function computes $ \int_{\ln M^{obs} - 7\sigma_{\ln M}}^{\ln M^{obs} + 7\sigma_{\ln M}} d\ln M \,
 * \frac{d^2N}{dzdlnM} * P(\ln M^{obs}|\ln M) $. The integral limits were determined requiring a precision
 * to five decimal places.
 *
 * Returns: a gdouble which corresponds to $ \int_{\ln M^{obs} - 7\sigma_{\ln M}}^{\ln M^{obs} + 7\sigma_{\ln M}} d\ln M \,
 * \frac{d^2N}{dzdlnM} * P(\ln M^{obs}|\ln M) $.
 */
gdouble
nc_cluster_abundance_d2n (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble lnM, gdouble z)
{
  NCM_UNUSED (clusterz);
  NCM_UNUSED (clusterm);
  
  return nc_halo_mass_function_d2n_dzdlnM (cad->mfp, cosmo, lnM, z);
}

static gdouble
_nc_cluster_abundance_z_intp_lnM_intp_d2N (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble lnM, gdouble z)
{
  const gdouble z_intp   = nc_cluster_redshift_intp (clusterz, cosmo, lnM, z);
  const gdouble lnM_intp = nc_cluster_mass_intp (clusterm, cosmo, lnM, z);
  gdouble d2NdzdlnM      = nc_halo_mass_function_d2n_dzdlnM (cad->mfp, cosmo, lnM, z);
  
  return z_intp * lnM_intp * d2NdzdlnM;
}

static gdouble
_nc_cluster_abundance_z_intp_lnM_intp_N_integrand (gdouble lnM, gdouble z, gpointer userdata)
{
  NcClusterAbundanceInt *obs_data = (NcClusterAbundanceInt *) userdata;
  
  return _nc_cluster_abundance_z_intp_lnM_intp_d2N (obs_data->cad, obs_data->cosmo, obs_data->clusterz, obs_data->clusterm, lnM, z);
}

static gdouble
_nc_cluster_abundance_z_intp_lnM_intp_bin_N_integrand (gdouble lnM, gdouble z, gpointer userdata)
{
  NcClusterAbundanceInt *obs_data = (NcClusterAbundanceInt *) userdata;
  const gdouble z_intp            = nc_cluster_redshift_intp_bin (obs_data->clusterz, obs_data->cosmo, lnM, z, obs_data->z_obs_lower, obs_data->z_obs_upper, obs_data->z_obs_params);
  const gdouble lnM_intp          = nc_cluster_mass_intp_bin (obs_data->clusterm, obs_data->cosmo, lnM, z, obs_data->lnM_obs_lower, obs_data->lnM_obs_upper, obs_data->lnM_obs_params);
  const gdouble d2NdzdlnM         = nc_halo_mass_function_d2n_dzdlnM (obs_data->cad->mfp, obs_data->cosmo, lnM, z);
  
  return z_intp * lnM_intp * d2NdzdlnM;
}

static gdouble
_nc_cluster_abundance_z_intp_lnM_intp_N (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm)
{
  gdouble N, zl, zu, lnMl, lnMu, err;
  NcClusterAbundanceInt obs_data;
  NcmIntegrand2dim integ;
  
  obs_data.cad      = cad;
  obs_data.cosmo    = cosmo;
  obs_data.clusterz = clusterz;
  obs_data.clusterm = clusterm;
  
  integ.f        = &_nc_cluster_abundance_z_intp_lnM_intp_N_integrand;
  integ.userdata = &obs_data;
  
  nc_cluster_redshift_n_limits (clusterz, cosmo, &zl, &zu);
  nc_cluster_mass_n_limits (clusterm, cosmo, &lnMl, &lnMu);
  
  ncm_integrate_2dim (&integ, lnMl, zl, lnMu, zu, NCM_DEFAULT_PRECISION, 0.0, &N, &err);
  
  return N;
}

static gdouble
_nc_cluster_abundance_z_intp_d2N (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble lnM, gdouble z)
{
  const gdouble z_intp = nc_cluster_redshift_intp (clusterz, cosmo, lnM, z);
  gdouble d2NdzdlnM    = nc_halo_mass_function_d2n_dzdlnM (cad->mfp, cosmo, lnM, z);
  
  NCM_UNUSED (clusterm);
  
  return z_intp * d2NdzdlnM;
}

static gdouble
_nc_cluster_abundance_z_intp_N_integrand (gdouble lnM, gdouble z, gpointer userdata)
{
  NcClusterAbundanceInt *obs_data = (NcClusterAbundanceInt *) userdata;
  
  return _nc_cluster_abundance_z_intp_d2N (obs_data->cad, obs_data->cosmo, obs_data->clusterz, obs_data->clusterm, lnM, z);
}

static gdouble
_nc_cluster_abundance_z_intp_N (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm)
{
  gdouble N, zl, zu, lnMl, lnMu, err;
  NcClusterAbundanceInt obs_data;
  NcmIntegrand2dim integ;
  
  obs_data.cad      = cad;
  obs_data.cosmo    = cosmo;
  obs_data.clusterz = clusterz;
  obs_data.clusterm = clusterm;
  
  integ.f        = &_nc_cluster_abundance_z_intp_N_integrand;
  integ.userdata = &obs_data;
  
  nc_cluster_redshift_n_limits (clusterz, cosmo, &zl, &zu);
  nc_cluster_mass_n_limits (clusterm, cosmo, &lnMl, &lnMu);
  
  ncm_integrate_2dim (&integ, lnMl, zl, lnMu, zu, NCM_DEFAULT_PRECISION, 0.0, &N, &err);
  
  return N;
}

static gdouble
_nc_cluster_abundance_lnM_intp_d2N (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble lnM, gdouble z)
{
  const gdouble lnM_intp = nc_cluster_mass_intp (clusterm, cosmo, lnM, z);
  gdouble d2NdzdlnM      = nc_halo_mass_function_d2n_dzdlnM (cad->mfp, cosmo, lnM, z);
  
  return lnM_intp * d2NdzdlnM;
}

static gdouble
_nc_cluster_abundance_lnM_intp_N_integrand (gdouble lnM, gdouble z, gpointer userdata)
{
  NcClusterAbundanceInt *obs_data = (NcClusterAbundanceInt *) userdata;
  
  return _nc_cluster_abundance_lnM_intp_d2N (obs_data->cad, obs_data->cosmo, obs_data->clusterz, obs_data->clusterm, lnM, z);
}

static gdouble
_nc_cluster_abundance_lnM_intp_N (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm)
{
  gdouble N, zl, zu, lnMl, lnMu, err;
  NcClusterAbundanceInt obs_data;
  NcmIntegrand2dim integ;
  
  obs_data.cad      = cad;
  obs_data.cosmo    = cosmo;
  obs_data.clusterz = clusterz;
  obs_data.clusterm = clusterm;
  
  integ.f        = &_nc_cluster_abundance_lnM_intp_N_integrand;
  integ.userdata = &obs_data;
  
  nc_cluster_redshift_n_limits (clusterz, cosmo, &zl, &zu);
  nc_cluster_mass_n_limits (clusterm, cosmo, &lnMl, &lnMu);
  
  ncm_integrate_2dim (&integ, lnMl, zl, lnMu, zu, NCM_DEFAULT_PRECISION, 0.0, &N, &err);
  
  return N;
}

static void
_nc_cluster_abundance_funcs (NcClusterAbundance *cad, NcClusterRedshift *clusterz, NcClusterMass *clusterm)
{
  gboolean z_intp   = ncm_model_check_impl_opt (NCM_MODEL (clusterz), NC_CLUSTER_REDSHIFT_INTP);
  gboolean lnM_intp = ncm_model_check_impl_opt (NCM_MODEL (clusterm), NC_CLUSTER_MASS_INTP);
  
  if (z_intp && lnM_intp)
  {
    cad->N        = &_nc_cluster_abundance_z_intp_lnM_intp_N;
    cad->intp_d2N = &_nc_cluster_abundance_z_intp_lnM_intp_d2N;
  }
  else if (z_intp && !lnM_intp)
  {
    cad->N        = &_nc_cluster_abundance_z_intp_N;
    cad->intp_d2N = &_nc_cluster_abundance_z_intp_d2N;
  }
  else if (!z_intp && lnM_intp)
  {
    cad->N        = &_nc_cluster_abundance_lnM_intp_N;
    cad->intp_d2N = &_nc_cluster_abundance_lnM_intp_d2N;
  }
  else
  {
    cad->N        = &nc_cluster_abundance_true_n;
    cad->intp_d2N = &nc_cluster_abundance_d2n;
  }
}

/**
 * nc_cluster_abundance_prepare:
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @clusterz: a #NcClusterRedshift
 * @clusterm: a #NcClusterMass
 *
 * This function prepares ...
 *
 */
void
nc_cluster_abundance_prepare (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm)
{
  _nc_cluster_abundance_funcs (cad, clusterz, clusterm);
  nc_cluster_redshift_n_limits (clusterz, cosmo, &cad->zi, &cad->zf);
  nc_cluster_mass_n_limits (clusterm, cosmo, &cad->lnMi, &cad->lnMf);

  nc_halo_mass_function_set_eval_limits (cad->mfp, cosmo, cad->lnMi, cad->lnMf, cad->zi, cad->zf);
  nc_halo_mass_function_prepare_if_needed (cad->mfp, cosmo);

  if (cad->zi == 0.0)
    cad->zi = 1.0e-6;
  
  cad->norma     = nc_cluster_abundance_true_n (cad, cosmo, clusterz, clusterm);
  cad->log_norma = log (cad->norma);
  
  ncm_model_ctrl_update (cad->ctrl_cosmo, NCM_MODEL (cosmo));
  ncm_model_ctrl_update (cad->ctrl_z, NCM_MODEL (clusterz));
  ncm_model_ctrl_update (cad->ctrl_m, NCM_MODEL (clusterm));
}

/**
 * nc_cluster_abundance_set_area:
 * @cad: a #NcClusterAbundance
 * @area: FIXME
 *
 * This function prepares ...
 *
 */
void
nc_cluster_abundance_set_area (NcClusterAbundance *cad, const gdouble area)
{
  nc_halo_mass_function_set_area (cad->mfp, area);
  ncm_model_ctrl_force_update (cad->ctrl_cosmo);
  ncm_model_ctrl_force_update (cad->ctrl_reion);
}

gdouble
_nc_cad_inv_dNdz_convergence_f (gdouble n, gdouble epsilon)
{
  return -log1p (epsilon - n);
}

static gdouble
_nc_cad_inv_dNdz_convergence_f_onemn (gdouble onemn, gdouble epsilon)
{
  return -log (epsilon + onemn);
}

/**
 * nc_cluster_abundance_prepare_inv_dNdz:
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @lnMi: logarithm base e of the minimum mass $\ln(M_i)$
 *
 * This function prepares a bidimensional spline...
 *
 */
void
nc_cluster_abundance_prepare_inv_dNdz (NcClusterAbundance *cad, NcHICosmo *cosmo, const gdouble lnMi)
{
  NcHaloMassFunctionSplineOptimize sp_optimize = NC_HALO_MASS_FUNCTION_SPLINE_Z;
  const gdouble delta_z = (cad->zf - cad->zi) / (cad->inv_z->len - 1.0);
  const gdouble delta_lnM = (cad->lnMf - lnMi) / (cad->inv_lnM->len - 1.0);
  const gdouble norma = nc_halo_mass_function_n (cad->mfp, cosmo, lnMi, cad->lnMf, cad->zi, cad->zf, NC_HALO_MASS_FUNCTION_SPLINE_LNM);
  gboolean use_spline = FALSE;
  guint middle = cad->inv_z->len / 2;
  gdouble z0 = cad->zi;
  guint i, j;
  
  g_assert (cad->zi != 0);
  
  {
    const gdouble zfm1   = cad->zf - delta_z;
    const gdouble lnMfm1 = cad->lnMf - delta_lnM;
    
    cad->z_epsilon   = fabs (nc_halo_mass_function_n (cad->mfp, cosmo, lnMi, cad->lnMf, zfm1, cad->zf, sp_optimize) / norma);
    cad->lnM_epsilon = fabs (nc_halo_mass_function_dn_dz (cad->mfp, cosmo, lnMfm1, cad->lnMf, cad->zf, use_spline) /
                             nc_halo_mass_function_dn_dz (cad->mfp, cosmo, lnMi, cad->lnMf, cad->zf, use_spline));
  }
  
  {
    gdouble zm = cad->zi + delta_z * middle;
    
    nc_cluster_abundance_prepare_inv_dNdlnM_z (cad, cosmo, lnMi, zm);
    
    ncm_vector_set (cad->inv_lnM_z->xv, 0, _nc_cad_inv_dNdz_convergence_f (0.0, cad->lnM_epsilon));
    ncm_matrix_set (cad->inv_lnM_z->zm, middle, 0, lnMi);
    
    for (j = 1; j < ncm_vector_len (cad->inv_lnM_z->xv) - 1; j++)
    {
      gdouble u2 = ncm_vector_get (cad->inv_lnM->xv, j);
      
      ncm_vector_set (cad->inv_lnM_z->xv, j, u2);
      ncm_matrix_set (cad->inv_lnM_z->zm, middle, j, ncm_spline_eval (cad->inv_lnM, u2));
    }
    
    ncm_vector_set (cad->inv_lnM_z->xv, j, _nc_cad_inv_dNdz_convergence_f_onemn (0.0, cad->lnM_epsilon));
    ncm_matrix_set (cad->inv_lnM_z->zm, middle, j, cad->lnMf);
  }
  
  nc_cluster_abundance_prepare_inv_dNdlnM_z (cad, cosmo, lnMi, z0);
  ncm_matrix_set (cad->inv_lnM_z->zm, 0, 0, lnMi);
  
  for (j = 1; j < ncm_vector_len (cad->inv_lnM_z->xv) - 1; j++)
  {
    gdouble u2 = ncm_vector_get (cad->inv_lnM_z->xv, j);
    
    ncm_matrix_set (cad->inv_lnM_z->zm, 0, j, ncm_spline_eval (cad->inv_lnM, u2));
  }
  
  ncm_matrix_set (cad->inv_lnM_z->zm, 0, j, cad->lnMf);
  
  {
    gdouble nztot = 0.0;
    gdouble f     = _nc_cad_inv_dNdz_convergence_f (0.0, cad->z_epsilon);
    
    ncm_vector_set (cad->inv_z->xv, 0, f);
    ncm_vector_set (cad->inv_z->yv, 0, z0);
    
    for (i = 1; i < cad->inv_z->len; i++)
    {
      gdouble z1 = cad->zi + delta_z * i;
      
      if (nztot < 0.99)
      {
        gdouble delta = nc_halo_mass_function_n (cad->mfp, cosmo, lnMi, cad->lnMf, z0, z1, sp_optimize) / norma;
        
        nztot += fabs (delta);
        f      = _nc_cad_inv_dNdz_convergence_f (nztot, cad->z_epsilon);
      }
      else
      {
        gdouble onemn = nc_halo_mass_function_n (cad->mfp, cosmo, lnMi, cad->lnMf, z1, cad->zf, sp_optimize) / norma;
        gdouble f_try = _nc_cad_inv_dNdz_convergence_f_onemn (onemn, cad->z_epsilon);
        
        if (f_try < f)
          f = f * (1.0 + GSL_SIGN (f) * 0.01);
        else
          f = f_try;
      }
      
      ncm_vector_set (cad->inv_z->xv, i, f);
      ncm_vector_set (cad->inv_z->yv, i, z1);
      
      z0 = z1;
      
      if (i == middle)
        continue;
      
      nc_cluster_abundance_prepare_inv_dNdlnM_z (cad, cosmo, lnMi, z1);
      
      ncm_matrix_set (cad->inv_lnM_z->zm, i, 0, lnMi);
      
      for (j = 1; j < ncm_vector_len (cad->inv_lnM_z->xv) - 1; j++)
      {
        gdouble u2 = ncm_vector_get (cad->inv_lnM_z->xv, j);
        
        ncm_matrix_set (cad->inv_lnM_z->zm, i, j, ncm_spline_eval (cad->inv_lnM, u2));
      }
      
      ncm_matrix_set (cad->inv_lnM_z->zm, i, j, cad->lnMf);
    }
  }
  
  ncm_spline2d_prepare (cad->inv_lnM_z);
  ncm_spline_prepare (cad->inv_z);
}

/**
 * nc_cluster_abundance_prepare_inv_dNdlnM_z:
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @lnMi: logarithm base e of the minimum mass $\ln(M_i)$
 * @z: redshift $z$
 *
 * This function prepares a spline where the x array corresponds to the value
 * of $\int_{\ln M_0} ^{\ln M_1} d^2N/dzd\ln M dM/ \int_lnMi^lnMf dN/dz dM$ given a redshift $z$
 * and the y array contains the values of logarithms base e of the mass.
 * It is used to generate a sample of $\ln M$ values.
 *
 */
void
nc_cluster_abundance_prepare_inv_dNdlnM_z (NcClusterAbundance *cad, NcHICosmo *cosmo, const gdouble lnMi, gdouble z)
{
  gboolean use_spline = FALSE;
  gdouble dNdz        = nc_halo_mass_function_dn_dz (cad->mfp, cosmo, lnMi, cad->lnMf, z, use_spline);
  gdouble lnM0        = lnMi;
  gdouble ntot        = 0.0;
  gdouble f           = _nc_cad_inv_dNdz_convergence_f (0.0, cad->lnM_epsilon);
  const gdouble dlnM  = (cad->lnMf - lnMi) / (cad->inv_lnM->len - 1.0);
  guint i;
  
  g_assert (z > 0.0);
  
  ncm_vector_set (cad->inv_lnM->xv, 0, f);
  ncm_vector_set (cad->inv_lnM->yv, 0, lnM0);
  
  for (i = 1; i < cad->inv_lnM->len; i++)
  {
    const gdouble lnM1 = lnMi + dlnM * i;
    
    if (ntot < 0.99)
    {
      const gdouble Delta = nc_halo_mass_function_dn_dz (cad->mfp, cosmo, lnM0, lnM1, z, use_spline) / dNdz;
      
      ntot += fabs (Delta);
      f     = _nc_cad_inv_dNdz_convergence_f (ntot, cad->lnM_epsilon);
    }
    else
    {
      const gdouble onemn = nc_halo_mass_function_dn_dz (cad->mfp, cosmo, lnM1, cad->lnMf, z, use_spline) / dNdz;
      const gdouble f_try = _nc_cad_inv_dNdz_convergence_f_onemn (fabs (onemn), cad->lnM_epsilon);
      
      if (f_try < f)
        f = f * (1.0 + GSL_SIGN (f) * 0.01);
      else
        f = f_try;
    }
    
    ncm_vector_set (cad->inv_lnM->xv, i, f);
    ncm_vector_set (cad->inv_lnM->yv, i, lnM1);
    lnM0 = lnM1;
  }
  
  ncm_spline_prepare (cad->inv_lnM);
}

/**
 * nc_cluster_abundance_true_n:
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @clusterz: a #NcClusterRedshift
 * @clusterm: a #NcClusterMass
 *
 * This function computes the total number of "true" clusters, i.e., halos
 * within redshift and mass intervals, and over a sky area.
 *
 * Returns: total number of halos (normalization)
 */
gdouble
nc_cluster_abundance_true_n (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm)
{
  NCM_UNUSED (clusterz);
  NCM_UNUSED (clusterm);
  
  return nc_halo_mass_function_n (cad->mfp, cosmo, cad->lnMi, cad->lnMf, cad->zi, cad->zf, NC_HALO_MASS_FUNCTION_SPLINE_LNM);
}

/**
 * nc_cluster_abundance_n:
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @clusterz: a #NcClusterRedshift
 * @clusterm: a #NcClusterMass
 *
 * This function computes the total number of clusters within specific redshift
 * and mass intervals, which are defined in @clusterz and @clusterm, respectively,
 * and over a sky area.
 *
 * Returns: total number of clusters.
 */
gdouble
nc_cluster_abundance_n (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm)
{
  return cad->N (cad, cosmo, clusterz, clusterm);
}

/**
 * nc_cluster_abundance_intp_d2n:
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @clusterz: a #NcClusterRedshift
 * @clusterm: a #NcClusterMass
 * @lnM: logarithm base e of the cluster mass
 * @z: redshift
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_cluster_abundance_intp_d2n (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble lnM, gdouble z)
{
  return cad->intp_d2N (cad, cosmo, clusterz, clusterm, lnM, z);
}

/**
 * nc_cluster_abundance_intp_bin_d2n:
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @clusterz: a #NcClusterRedshift
 * @clusterm: a #NcClusterMass
 * @lnM_obs_lower: (array) (element-type gdouble): FIXME
 * @lnM_obs_upper: (array) (element-type gdouble): FIXME
 * @lnM_obs_params: (array) (element-type gdouble) (allow-none): FIXME
 * @z_obs_lower: (array) (element-type gdouble): FIXME
 * @z_obs_upper: (array) (element-type gdouble): FIXME
 * @z_obs_params: (array) (element-type gdouble) (allow-none): FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_cluster_abundance_intp_bin_d2n (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, const gdouble *z_obs_lower, const gdouble *z_obs_upper, const gdouble *z_obs_params)
{
  gdouble N, zl, zu, lnMl, lnMu, err;
  NcClusterAbundanceInt obs_data;
  NcmIntegrand2dim integ;
  
  obs_data.cad            = cad;
  obs_data.cosmo          = cosmo;
  obs_data.clusterz       = clusterz;
  obs_data.clusterm       = clusterm;
  obs_data.lnM_obs_lower  = lnM_obs_lower;
  obs_data.lnM_obs_upper  = lnM_obs_upper;
  obs_data.lnM_obs_params = lnM_obs_params;
  obs_data.z_obs_lower    = z_obs_lower;
  obs_data.z_obs_upper    = z_obs_upper;
  obs_data.z_obs_params   = z_obs_params;
  
  integ.f        = &_nc_cluster_abundance_z_intp_lnM_intp_bin_N_integrand;
  integ.userdata = &obs_data;
  
  nc_cluster_redshift_p_bin_limits (clusterz, cosmo, z_obs_lower, z_obs_upper, z_obs_params, &zl, &zu);
  nc_cluster_mass_p_bin_limits (clusterm, cosmo, lnM_obs_lower, lnM_obs_upper, lnM_obs_params, &lnMl, &lnMu);
  
  ncm_integrate_2dim (&integ, lnMl, zl, lnMu, zu, NCM_DEFAULT_PRECISION, 0.0, &N, &err);

  return N;
}

/**
 * nc_cluster_abundance_bin_realization: (skip)
 * @zr: FIXME
 * @h: FIXME
 *
 * FIXME
 */
void
nc_cluster_abundance_bin_realization (GArray *zr, gsl_histogram **h)
{
  guint i;
  
  for (i = 0; i < zr->len; i++)
  {
    int j;
    const gdouble z = g_array_index (zr, gdouble, i);
    
    for (j = 0; h[j] != NULL; j++)
      gsl_histogram_increment (h[j], z);
  }
}

void
nc_cluster_abundance_realizations_save_to_file (GPtrArray *realizations, gchar *filename)
{
  FILE *out;
  guint i, j;
  
  out = fopen (filename, "w");
  
  for (j = 0; j < realizations->len; j++)
  {
    GArray *z_real = g_ptr_array_index (realizations, j);
    
    for (i = 0; i < z_real->len; i++)
    {
      gdouble zi = g_array_index (z_real, gdouble, i);
      
      fprintf (out, "%.15g\n", zi);
      fflush (out);
    }
    
    fprintf (out, "\n\n");
  }
  
  fclose (out);
}

/**
 * nc_cluster_abundance_realizations_read_from_file:
 * @file_realization: (array length=n_realizations): is the file's name which contains the redshift values of the clusters obtained with random Poisson generator.
 * @n_realizations: is the number of realizations generated.
 *
 * GPtrArray *realizations is an array of array with the z values of all realizations. To complete...
 *
 * Returns: (transfer full): FIXME
 */
GPtrArray *
nc_cluster_abundance_realizations_read_from_file (gchar *file_realization, guint n_realizations)
{
  FILE *frealization = fopen (file_realization, "r");
  GPtrArray *realizations = g_ptr_array_sized_new (n_realizations);
  guint i, j;
  glong file_position, goby;
  
  if (frealization == NULL)
    g_error ("abundance_random_generator_read_from_file: file %s, do not exist.", file_realization);
  
  file_position = ftell (frealization);
  
  for (j = 0; j < n_realizations; j++)
  {
    guint z_total;
    gint n_enter = 0;
    gchar line[5000];
    
    z_total = 0; /* Counting the number of redshifts, which corresponds to the total number of clusters, for each realization. */
    
    while (fgets (line, 5000, frealization) != NULL)
    {
      if (strlen (line) == 1)
      {
        n_enter++;
        
        if (n_enter == 2)
          break;
        
        continue;
      }
      
      if (n_enter == 1)
        g_error ("Error, it must not exist just one 'enter'.");
      
      n_enter = 0;
      z_total++;
    }
    
    {
      GArray *realization_z = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), z_total); /* Array which contains the z values of each realization. */
      
      goby = ftell (frealization) - file_position;
      fseek (frealization, -goby, SEEK_CUR);
      
      for (i = 0; i < z_total; i++)
      {
        gdouble z;
        
        if (fscanf (frealization, "%lg\n", &z) < 1)
          g_error ("nc_cluster_abundance_realizations_read_from_file: io error");
        
        /*        printf ("%lg %u\n", z, z_total); */
        g_array_append_val (realization_z, z);
      }
      
      g_ptr_array_add (realizations, realization_z);
    }
    
    if (fscanf (frealization, " ") < 1)
      g_error ("nc_cluster_abundance_realizations_read_from_file: io error");
    
    if (fscanf (frealization, " ") < 1)
      g_error ("nc_cluster_abundance_realizations_read_from_file: io error");
    
    file_position = ftell (frealization);
  }
  
  fclose (frealization);
  
  return realizations;
}

/********* Function to compute the mean bias ************************/

static gdouble
_nc_ca_mean_bias_numerator_integrand (gdouble lnM, gpointer params)
{
  NcClusterAbundanceInt *obs_data = (NcClusterAbundanceInt *) params;
  NcClusterAbundance *cad         = obs_data->cad;
  
  gdouble dbdlnM = nc_halo_bias_func_integrand (cad->mbiasf, obs_data->cosmo, lnM, obs_data->z);
  
  return dbdlnM;
}

gdouble
nc_ca_mean_bias_numerator (NcClusterAbundance *cad, NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  NcClusterAbundanceInt obs_data;
  gdouble mean_bias_numerator;
  gsl_function F;
  
  obs_data.cad   = cad;
  obs_data.cosmo = cosmo;
  
  F.function = &_nc_ca_mean_bias_numerator_integrand;
  F.params   = &obs_data;
  
  {
    gdouble res, err;
    gdouble lnMf = cad->lnMf;
    
    obs_data.z   = z;
    obs_data.lnM = lnM;
    /*printf ("%5.5e, %5.5e\n", exp(lnM), exp(lnMf)); */
    
    {
      gsl_integration_workspace **w = ncm_integral_get_workspace ();
      
      gsl_integration_qag (&F, lnM, lnMf, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &res, &err);
      ncm_memory_pool_return (w);
    }
    mean_bias_numerator = res;
  }
  
  return mean_bias_numerator;
}

static gdouble
_nc_ca_mean_bias_denominator_integrand (gdouble lnM, gpointer params)
{
  NcClusterAbundanceInt *obs_data = (NcClusterAbundanceInt *) params;
  NcClusterAbundance *cad         = obs_data->cad;
  
  gdouble n = nc_halo_mass_function_dn_dlnM (cad->mfp, obs_data->cosmo, lnM, obs_data->z);
  
  /*printf ("% 20.8e % 20.8g % 20.15g\n", obs_data->lnMobs, obs_data->zp, dbdlnM); */
  return n;
}

gdouble
nc_ca_mean_bias_denominator (NcClusterAbundance *cad, NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  NcClusterAbundanceInt obs_data;
  gdouble mean_bias_denominator;
  gsl_function F;
  
  obs_data.cad   = cad;
  obs_data.cosmo = cosmo;
  
  F.function = &_nc_ca_mean_bias_denominator_integrand;
  F.params   = &obs_data;
  
  {
    gdouble res, err;
    gdouble lnMf = cad->lnMf;
    
    obs_data.z   = z;
    obs_data.lnM = lnM;
    {
      gsl_integration_workspace **w = ncm_integral_get_workspace ();
      
      gsl_integration_qag (&F, lnM, lnMf, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &res, &err);
      ncm_memory_pool_return (w);
    }
    mean_bias_denominator = res;
  }
  
  return mean_bias_denominator;
}

gdouble
nc_ca_mean_bias (NcClusterAbundance *cad, NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  gdouble numerator   = nc_ca_mean_bias_numerator (cad, cosmo, lnM, z);
  gdouble denominator = nc_ca_mean_bias_denominator (cad, cosmo, lnM, z);
  gdouble mean_bias   = numerator / denominator;
  
  return mean_bias;
}

static gdouble
_nc_ca_mean_bias_Mobs_numerator_integrand (gdouble lnMobs, gpointer params)
{
  NcClusterAbundanceInt *obs_data = (NcClusterAbundanceInt *) params;
  NcClusterAbundance *cad         = obs_data->cad;
  
  /* In this case zp is the true redshift, i.e., without uncertainty. */
  gdouble p_M_Mobs = 0.0; /*_nc_cluster_abundance_lognormal_mass_dist (obs_data, lnMobs); */
  gdouble dbdlnM   = nc_halo_bias_func_integrand (cad->mbiasf, obs_data->cosmo, lnMobs, obs_data->z);
  
  g_assert_not_reached ();
  
  return dbdlnM * p_M_Mobs;
}

gdouble
nc_ca_mean_bias_Mobs_numerator (NcClusterAbundance *cad, NcHICosmo *cosmo, gdouble lnMobs, gdouble z)
{
  NcClusterAbundanceInt obs_data;
  gdouble mean_bias_Mobs_numerator;
  gsl_function F;
  
  obs_data.cad   = cad;
  obs_data.cosmo = cosmo;
  
  F.function = &_nc_ca_mean_bias_Mobs_numerator_integrand;
  F.params   = &obs_data;
  
  {
    obs_data.z   = z;
    obs_data.lnM = lnMobs;
    {
      gdouble res, err;
      gdouble lnMl = 0.0, lnMu = 0.0;
      
      /*lnMl = GSL_MAX (lnMobs - 7.0 * obs_data.cad->lnM_sigma0, LNM_MIN); */
      /*lnMu = GSL_MIN (lnMobs + 7.0 * obs_data.cad->lnM_sigma0, 16.3 * M_LN10); */
      /*lnMl = lnMobs - 7.0 * obs_data.cad->lnM_sigma0; */
      /*lnMu = lnMobs + 7.0 * obs_data.cad->lnM_sigma0; */
      g_assert_not_reached ();
      {
        gsl_integration_workspace **w = ncm_integral_get_workspace ();
        
        gsl_integration_qag (&F, lnMl, lnMu, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &res, &err);
        ncm_memory_pool_return (w);
      }
      mean_bias_Mobs_numerator = res;
    }
  }
  
  return mean_bias_Mobs_numerator;
}

static gdouble
_nc_ca_mean_bias_Mobs_denominator_integrand (gdouble lnMobs, gpointer params)
{
  NcClusterAbundanceInt *obs_data = (NcClusterAbundanceInt *) params;
  NcClusterAbundance *cad         = obs_data->cad;
  
  /* In this case zp is the true redshift, i.e., without uncertainty. */
  gdouble p_M_Mobs = 0.0; /*_nc_cluster_abundance_lognormal_mass_dist (obs_data, lnMobs); */
  gdouble n        = nc_halo_mass_function_dn_dlnM (cad->mfp, obs_data->cosmo, lnMobs, obs_data->z);
  
  g_assert_not_reached ();
  
  return n * p_M_Mobs;
}

gdouble
nc_ca_mean_bias_Mobs_denominator (NcClusterAbundance *cad, NcHICosmo *cosmo, gdouble lnMobs, gdouble z)
{
  NcClusterAbundanceInt obs_data;
  gdouble mean_bias_Mobs_denominator;
  gsl_function F;
  
  obs_data.cad   = cad;
  obs_data.cosmo = cosmo;
  
  F.function = &_nc_ca_mean_bias_Mobs_denominator_integrand;
  F.params   = &obs_data;
  
  {
    obs_data.z   = z;
    obs_data.lnM = lnMobs;
    {
      gdouble res, err;
      gdouble lnMl = 0.0, lnMu = 0.0;
      
      /*lnMl = GSL_MAX (lnMobs - 7.0 * obs_data.cad->lnM_sigma0, LNM_MIN); */
      /*lnMu = GSL_MIN (lnMobs + 7.0 * obs_data.cad->lnM_sigma0, 16.3 * M_LN10); */
      /*lnMl = lnMobs - 7.0 * obs_data.cad->lnM_sigma0; */
      /*lnMu = lnMobs + 7.0 * obs_data.cad->lnM_sigma0; */
      g_assert_not_reached ();
      {
        gsl_integration_workspace **w = ncm_integral_get_workspace ();
        
        gsl_integration_qag (&F, lnMl, lnMu, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_ABUNDANCE_DEFAULT_INT_KEY, *w, &res, &err);
        ncm_memory_pool_return (w);
      }
      mean_bias_Mobs_denominator = res;
    }
  }
  
  return mean_bias_Mobs_denominator;
}

/**
 * nc_cluster_abundance_prepare_if_needed:
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @clusterz: a #NcClusterRedshift
 * @clusterm: a #NcClusterMass
 *
 * FIXME
 *
 */


static gdouble
nc_cluster_abundance_intp_d2n_bias_integrand (gdouble lnM, gdouble z, gpointer userdata)
{
  NcClusterAbundanceInt *obs_data = (NcClusterAbundanceInt *) userdata;
  NcClusterAbundance *cad         = obs_data->cad;


  const gdouble p_z_zr            = nc_cluster_redshift_p (obs_data->clusterz, obs_data->cosmo, lnM, z, obs_data->z_obs, obs_data->z_obs_params);
  const gdouble p_M_Mobs          = nc_cluster_mass_p (obs_data->clusterm, obs_data->cosmo, lnM, z, obs_data->lnM_obs, obs_data->lnM_obs_params);
  const gdouble dbdlnM            = nc_halo_bias_func_integrand (cad->mbiasf, obs_data->cosmo, lnM, z);
  
  return p_z_zr * p_M_Mobs * dbdlnM;
}




/**
 * nc_cluster_abundance_intp_d2n_bias:
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @clusterz: a #NcClusterRedshift
 * @clusterm: a #NcClusterMass
 * @lnM_obs: (array) (element-type gdouble): FIXME
 * @lnM_obs_params: (array) (element-type gdouble) (allow-none): FIXME
 * @z_obs: (array) (element-type gdouble): FIXME
 * @z_obs_params: (array) (element-type gdouble) (allow-none): FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */




gdouble
nc_cluster_abundance_intp_d2n_bias (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm, gdouble *lnM_obs, gdouble *lnM_obs_params, gdouble *z_obs, gdouble *z_obs_params)
{
  gdouble d2N_bias, zl, zu, lnMl, lnMu, err;
  NcClusterAbundanceInt obs_data;
  NcmIntegrand2dim integ;
  
  obs_data.cad            = cad;
  obs_data.cosmo          = cosmo;
  obs_data.clusterz       = clusterz;
  obs_data.clusterm       = clusterm;
  obs_data.lnM_obs        = lnM_obs;
  obs_data.lnM_obs_params = lnM_obs_params;
  obs_data.z_obs          = z_obs;
  obs_data.z_obs_params   = z_obs_params;
  
  integ.f        = &nc_cluster_abundance_intp_d2n_bias_integrand;
  integ.userdata = &obs_data;
  
  nc_cluster_redshift_p_limits (clusterz, cosmo, z_obs, z_obs_params, &zl, &zu);
  nc_cluster_mass_p_limits (clusterm, cosmo, lnM_obs, lnM_obs_params, &lnMl, &lnMu);
  
  ncm_integrate_2dim (&integ, lnMl, zl, lnMu, zu, NCM_DEFAULT_PRECISION, 0.0, &d2N_bias, &err);
  
  return d2N_bias;
}


static gdouble
nc_cluster_abundance_intp_bin_d2n_bias_integrand (gdouble lnM, gdouble z, gpointer userdata)
{
  NcClusterAbundanceInt *obs_data = (NcClusterAbundanceInt *) userdata;
  NcClusterAbundance *cad         = obs_data->cad;


const gdouble z_intp            = nc_cluster_redshift_intp_bin (obs_data->clusterz, obs_data->cosmo, lnM, z, obs_data->z_obs_lower, obs_data->z_obs_upper, obs_data->z_obs_params);
const gdouble lnM_intp          = nc_cluster_mass_intp_bin (obs_data->clusterm, obs_data->cosmo, lnM, z, obs_data->lnM_obs_lower, obs_data->lnM_obs_upper, obs_data->lnM_obs_params);
const gdouble dbdlnM            = nc_halo_bias_func_integrand (cad->mbiasf, obs_data->cosmo, lnM, z);
  
  return z_intp * lnM_intp * dbdlnM;
}



/**
 * nc_cluster_abundance_intp_bin_d2n_bias:
 * @cad: a #NcClusterAbundance
 * @cosmo: a #NcHICosmo
 * @clusterz: a #NcClusterRedshift
 * @clusterm: a #NcClusterMass
 * @lnM_obs_lower: (array) (element-type gdouble): FIXME
 * @lnM_obs_upper: (array) (element-type gdouble): FIXME
 * @lnM_obs_params: (array) (element-type gdouble) (allow-none): FIXME
 * @z_obs_lower: (array) (element-type gdouble): FIXME
 * @z_obs_upper: (array) (element-type gdouble): FIXME
 * @z_obs_params: (array) (element-type gdouble) (allow-none): FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */

gdouble
nc_cluster_abundance_intp_bin_d2n_bias (NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, const gdouble *z_obs_lower, const gdouble *z_obs_upper, const gdouble *z_obs_params)
{
  gdouble N_bias, zl, zu, lnMl, lnMu, err;
  NcClusterAbundanceInt obs_data;
  NcmIntegrand2dim integ;
  
  obs_data.cad            = cad;
  obs_data.cosmo          = cosmo;
  obs_data.clusterz       = clusterz;
  obs_data.clusterm       = clusterm;
  obs_data.lnM_obs_lower  = lnM_obs_lower;
  obs_data.lnM_obs_upper  = lnM_obs_upper;
  obs_data.lnM_obs_params = lnM_obs_params;
  obs_data.z_obs_lower    = z_obs_lower;
  obs_data.z_obs_upper    = z_obs_upper;
  obs_data.z_obs_params   = z_obs_params;
  
  integ.f        = &nc_cluster_abundance_intp_bin_d2n_bias_integrand;
  integ.userdata = &obs_data;
  
  nc_cluster_redshift_p_bin_limits (clusterz, cosmo, z_obs_lower, z_obs_upper, z_obs_params, &zl, &zu);
  nc_cluster_mass_p_bin_limits (clusterm, cosmo, lnM_obs_lower, lnM_obs_upper, lnM_obs_params, &lnMl, &lnMu);
  
  ncm_integrate_2dim (&integ, lnMl, zl, lnMu, zu, NCM_DEFAULT_PRECISION, 0.0, &N_bias, &err);

  return N_bias;
}


static gdouble
_nc_cluster_abundance_mean_bias_integrand (gdouble lnM, gdouble z, gpointer userdata)
{
  NcClusterAbundanceInt *obs_data = (NcClusterAbundanceInt *) userdata;
  NcClusterAbundance *cad         = obs_data->cad;
  const gdouble z_intp   = nc_cluster_redshift_intp (obs_data->clusterz, obs_data->cosmo, lnM, z);
  const gdouble lnM_intp = nc_cluster_mass_intp (obs_data->clusterm, obs_data->cosmo, lnM, z);
  const gdouble dbdlnM   = nc_halo_bias_func_integrand (cad->mbiasf, obs_data->cosmo, lnM, z);
  
  return z_intp * lnM_intp * dbdlnM;
}



gdouble
nc_cluster_abundance_mean_bias(NcClusterAbundance *cad, NcHICosmo *cosmo, NcClusterRedshift *clusterz, NcClusterMass *clusterm)
{
gdouble N_bias, zl, zu, lnMl, lnMu, err;
  NcClusterAbundanceInt obs_data;
  NcmIntegrand2dim integ;
  
  obs_data.cad      = cad;
  obs_data.cosmo    = cosmo;
  obs_data.clusterz = clusterz;
  obs_data.clusterm = clusterm;
  
  integ.f        = &_nc_cluster_abundance_mean_bias_integrand;
  integ.userdata = &obs_data;
  
  nc_cluster_redshift_n_limits (clusterz, cosmo, &zl, &zu);
  nc_cluster_mass_n_limits (clusterm, cosmo, &lnMl, &lnMu);
  
  ncm_integrate_2dim (&integ, lnMl, zl, lnMu, zu, NCM_DEFAULT_PRECISION, 0.0, &N_bias, &err);
  
  return N_bias;

}