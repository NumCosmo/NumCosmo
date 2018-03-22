/***************************************************************************
 *            nc_reduced_shear_cluster_mass.c
 *
 *  Mon Mar 19 15:42:23 2018
 *  Copyright  2018  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2018 <pennalima@gmail.com> 
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
 * SECTION:nc_reduced_shear_cluster_mass
 * @title: NcReducedShearClusterMass
 * @short_description: FIXME cluster mass estimation via reduced shear
 *
 * FIXME 
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_reduced_shear_cluster_mass.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/integral.h"
#include "math/memory_pool.h"
#include "math/Faddeeva.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_roots.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcReducedShearClusterMass, nc_reduced_shear_cluster_mass, NCM_TYPE_MODEL);

#define VECTOR  (NCM_MODEL (rscm)->params)
#define A  (ncm_vector_fast_get (VECTOR, NC_REDUCED_SHEAR_CLUSTER_MASS_A))
#define B (ncm_vector_fast_get (VECTOR, NC_REDUCED_SHEAR_CLUSTER_MASS_B))
#define C    (ncm_vector_fast_get (VECTOR, NC_REDUCED_SHEAR_CLUSTER_MASS_C))
#define XP    (ncm_vector_fast_get (VECTOR, NC_REDUCED_SHEAR_CLUSTER_MASS_XP))
#define VSIGMA  (ncm_vector_fast_get (VECTOR, NC_REDUCED_SHEAR_CLUSTER_MASS_VSIGMA))
#define VGAMMA  (ncm_vector_fast_get (VECTOR, NC_REDUCED_SHEAR_CLUSTER_MASS_VGAMMA))

enum
{
  PROP_0,
  PROP_R,
  PROP_NZBINS,
  PROP_SIZE,
};

static void
nc_reduced_shear_cluster_mass_init (NcReducedShearClusterMass *rscm)
{
  rscm->R_Mpc  = 0.0;
  rscm->nzbins = 1;
  rscm->T      = gsl_multifit_fdfsolver_lmsder;
  rscm->s      = gsl_multifit_fdfsolver_alloc (rscm->T, 4, 2);
}

static void
_nc_reduced_shear_cluster_mass_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcReducedShearClusterMass *rscm = NC_REDUCED_SHEAR_CLUSTER_MASS (object);
  g_return_if_fail (NC_IS_REDUCED_SHEAR_CLUSTER_MASS (object));

  switch (prop_id)
  {
    case PROP_R:
      rscm->R_Mpc = g_value_get_double (value);
      break;
    case PROP_NZBINS:
      rscm->nzbins = g_value_get_uint (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_reduced_shear_cluster_mass_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcReducedShearClusterMass *rscm = NC_REDUCED_SHEAR_CLUSTER_MASS (object);
  g_return_if_fail (NC_IS_REDUCED_SHEAR_CLUSTER_MASS (object));

  switch (prop_id)
  {
    case PROP_R:
      g_value_set_double (value, rscm->R_Mpc);
      break;
    case PROP_NZBINS:
      g_value_set_uint (value, rscm->nzbins);
      break;  
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_reduced_shear_cluster_mass_dispose (GObject *object)
{
  /*NcReducedShearClusterMass *rscm = NC_REDUCED_SHEAR_CLUSTER_MASS (object);*/
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_reduced_shear_cluster_mass_parent_class)->dispose (object);
}

static void
_nc_reduced_shear_cluster_mass_finalize (GObject *object)
{
  NcReducedShearClusterMass *rscm = NC_REDUCED_SHEAR_CLUSTER_MASS (object);
  gsl_multifit_fdfsolver_free (rscm->s);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_reduced_shear_cluster_mass_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_reduced_shear_cluster_mass, NC_TYPE_REDUCED_SHEAR_CLUSTER_MASS);

static void
nc_reduced_shear_cluster_mass_class_init (NcReducedShearClusterMassClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_reduced_shear_cluster_mass_set_property;
  model_class->get_property = &_nc_reduced_shear_cluster_mass_get_property;
  object_class->dispose     = &_nc_reduced_shear_cluster_mass_dispose;
  object_class->finalize    = &_nc_reduced_shear_cluster_mass_finalize;

  ncm_model_class_set_name_nick (model_class, "Lensing observable for cluster mass estimation: reduced shear", "ReducedShearClusterMass");
  ncm_model_class_add_params (model_class, NC_REDUCED_SHEAR_CLUSTER_MASS_SPARAM_LEN, 0, PROP_SIZE);
  

   /**
   * NcClusterPseudoCounts:R_Mpc:
   *
   * Scale/distance in Mpc from the center of the lens.
   */
  g_object_class_install_property (object_class,
                                   PROP_R,
                                   g_param_spec_double ("R",
                                                        NULL,
                                                        "Distance from the center of the lens",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
   /**
   * NcClusterPseudoCounts:nzbins:
   *
   * Number of redshift bins.
   */
  g_object_class_install_property (object_class,
                                   PROP_NZBINS,
                                   g_param_spec_uint ("number-z-bins",
                                                      NULL,
                                                      "Number of redshift bins",
                                                      1, G_MAXUINT, 10,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcReducedShearClusterMass:a:
   * 
   * FIXME 
   *  
   */
  /**
   * NcReducedShearClusterMass:a-fit:
   * 
   * FIXME
   *  
   */
  ncm_model_class_set_sparam (model_class, NC_REDUCED_SHEAR_CLUSTER_MASS_A, "a", "a",
                              0.0, 16.0, 2.0,
                              NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_PARAMS_ABSTOL, NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_A,
                              NCM_PARAM_TYPE_FIXED);
  /**
   * NcReducedShearClusterMass:b:
   * 
   * FIXME
   */
  ncm_model_class_set_sparam (model_class, NC_REDUCED_SHEAR_CLUSTER_MASS_B, "b", "b",
                              0.0,  0.9, 2.0e-2,
                              NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_PARAMS_ABSTOL, NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_B,
                              NCM_PARAM_TYPE_FIXED);
  /**
   * NcReducedShearClusterMass:c:
   * 
   * FIXME
   */
  ncm_model_class_set_sparam (model_class, NC_REDUCED_SHEAR_CLUSTER_MASS_C, "c", "c",
                              0.0,  2.0, 1.0e-2,
                              NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_PARAMS_ABSTOL, NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_C,
                              NCM_PARAM_TYPE_FIXED);
  /**
   * NcReducedShearClusterMass:xp:
   * 
   * FIXME
   */
  ncm_model_class_set_sparam (model_class, NC_REDUCED_SHEAR_CLUSTER_MASS_XP, "xp", "xp",
                              0.0,  2.0, 1.0e-2,
                              NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_PARAMS_ABSTOL, NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_XP,
                              NCM_PARAM_TYPE_FIXED);
  /**
   * NcReducedShearClusterMass:Vsigma:
   * 
   * Voigt profile parameter, $\sigma$ is the standard deviation of the Gaussian distribution. Range: $\sigma \in [0.15, 0.5]$.
   */
  ncm_model_class_set_sparam (model_class, NC_REDUCED_SHEAR_CLUSTER_MASS_VSIGMA, "\\sigma", "sigma",
                              0.15,  0.5, 1.0e-2,
                              NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_PARAMS_ABSTOL, NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_VSIGMA,
                              NCM_PARAM_TYPE_FIXED);
  /**
   * NcReducedShearClusterMass:VGamma:
   * 
   * Voigt profile parameter, $\Gamma$ is the width of the  Lorentzian profile. Range: $\Gamma \in [0.003, 0.1]$.
   */
  ncm_model_class_set_sparam (model_class, NC_REDUCED_SHEAR_CLUSTER_MASS_VGAMMA, "\\Gamma", "Gamma",
                              3e-3,  0.1, 1.0e-2,
                              NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_PARAMS_ABSTOL, NC_REDUCED_SHEAR_CLUSTER_MASS_DEFAULT_VGAMMA,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  ncm_mset_model_register_id (model_class,
                              "NcReducedShearClusterMass",
                              "Lensing observable - cluster mass estimation: reduced shear",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);
}

/**
 * nc_reduced_shear_cluster_mass_new: 
 * @nclusters: total number of clusters (resample)
 *
 * This function instantiates a new object of type #NcReducedShearClusterMass.
 *
 * Returns: A new #NcReducedShearClusterMass.
 */
NcReducedShearClusterMass *
nc_reduced_shear_cluster_mass_new (guint nclusters)
{
  NcReducedShearClusterMass *rscm = g_object_new (NC_TYPE_REDUCED_SHEAR_CLUSTER_MASS,
                                          NULL);
  return rscm;
}

/**
 * nc_reduced_shear_cluster_mass_ref:
 * @rscm: a #NcReducedShearClusterMass
 *
 * Increases the reference count of @rscm by one.
 *
 * Returns: (transfer full): @rscm
 */
NcReducedShearClusterMass *
nc_reduced_shear_cluster_mass_ref (NcReducedShearClusterMass *rscm)
{
  return g_object_ref (rscm);
}

/**
 * nc_reduced_shear_cluster_mass_free:
 * @rscm: a #NcReducedShearClusterMass
 *
 * Atomically decreases the reference count of @rscm by one. If the reference count drops to 0,
 * all memory allocated by @rscm is released.
 *  
 */
void
nc_reduced_shear_cluster_mass_free (NcReducedShearClusterMass *rscm)
{
  g_object_unref (rscm);
}

/**
 * nc_reduced_shear_cluster_mass_clear:
 * @rscm: a #NcReducedShearClusterMass
 *
 * The reference count of @rscm is decreased and the pointer is set to NULL.
 *
 */
void
nc_reduced_shear_cluster_mass_clear (NcReducedShearClusterMass **rscm)
{
  g_clear_object (rscm);
}

//////////////////////////////////////////////////////////////////////////////

typedef struct _integrand_data
{
  NcReducedShearClusterMass *rscm;
  NcWLSurfaceMassDensity *smd;
  NcDensityProfile *dp;
  NcHICosmo *cosmo;
  gdouble R;
  gdouble z;
  gdouble g_obs;
} integrand_data;

static gdouble
_posterior_integrand_no_shear_calibration (gdouble z, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcReducedShearClusterMass *rscm = data->rscm;
  gdouble g  = nc_wl_surface_mass_density_reduced_shear_infinity (data->smd, data->dp, data->cosmo, data->R);
  gdouble mu = data->g_obs - g;
  double complex Z = (mu + I * VGAMMA) / (M_SQRT2 * VSIGMA);
  gdouble voigt = creal (Faddeeva_w (Z, 1e-6)) / (M_SQRT2 * M_SQRTPI * VSIGMA) ;
  gdouble pz = 0.0; /* FIXME histogram, spline???? */
  const gdouble small = exp (-200.0);
  gdouble result;
  
  if (pz == 0.0)
    result = small;
  else
  {
    result =  voigt *pz + small;
  }  
  return result;
}

/**
 * nc_reduced_shear_cluster_mass_posterior_no_shear_calibration:
 * @rscm: a #NcReducedShearClusterMass
 * @cosmo: a #NcHICosmo 
 * @z: spectroscopic redshift
 * @g_obs: observed reduced shear
 * 
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_reduced_shear_cluster_mass_posterior_no_shear_calibration (NcReducedShearClusterMass *rscm, NcHICosmo *cosmo, const gdouble z, const gdouble g_obs)
{
  integrand_data data;
  gdouble P, err;

  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  data.rscm  = rscm;
  data.cosmo = cosmo;
  data.g_obs = g_obs;
  data.z     = z;
    
  F.function = &_posterior_integrand_no_shear_calibration;
  F.params = &data;

  {
    gdouble lnM_min, lnM_max;
    lnM_min = log (1.0e12); //32.292; //    
    lnM_max = log (1.0e16); //34.561; //log (1.0e16);     
    gsl_integration_qag (&F, lnM_min, lnM_max, 0.0, 1.0e2 * NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &P, &err);
  }

  ncm_memory_pool_return (w);
  /*printf ("numerator = %.8g err = %.8g\n", P, err / P);*/
  return P;
  
}
