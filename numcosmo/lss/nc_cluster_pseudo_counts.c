/***************************************************************************
 *            nc_cluster_pseudo_counts.c
 *
 *  Mon Mar 30 02:04:07 2015
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
 * SECTION:nc_cluster_pseudo_counts
 * @title: NcClusterPseudoCounts
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_pseudo_counts.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/integral.h"
#include "math/memory_pool.h"
#include "levmar/levmar.h"

#include <gsl/gsl_roots.h>

G_DEFINE_TYPE (NcClusterPseudoCounts, nc_cluster_pseudo_counts, NCM_TYPE_MODEL);

#define VECTOR  (NCM_MODEL (cpc)->params)
#define LNMCUT  (ncm_vector_get (VECTOR, NC_CLUSTER_PSEUDO_COUNTS_LNMCUT))
#define SD_MCUT (ncm_vector_get (VECTOR, NC_CLUSTER_PSEUDO_COUNTS_SD_MCUT))
#define ZMIN    (ncm_vector_get (VECTOR, NC_CLUSTER_PSEUDO_COUNTS_ZMIN))
#define DELTAZ  (ncm_vector_get (VECTOR, NC_CLUSTER_PSEUDO_COUNTS_DELTAZ))

enum
{
  PROP_0,
  PROP_MASS_FUNCTION,
  PROP_SIZE,
};

/**
 * nc_cluster_pseudo_counts_new: 
 * @mfp: a #NcMassFunction
 *
 * This function allocates memory for a new #NcClusterPseudoCounts object and sets its properties to the values from
 * the input argument.
 *
 * Returns: A new #NcClusterPseudoCounts.
 */
NcClusterPseudoCounts *
nc_cluster_pseudo_counts_new (NcMassFunction *mfp)
{
  NcClusterPseudoCounts *cpc = g_object_new (NC_TYPE_CLUSTER_PSEUDO_COUNTS,
                                          "mass-function", mfp,                                        
                                          NULL);
  return cpc;
}

/**
 * nc_cluster_pseudo_counts_copy:
 * @cpc: a #NcClusterPseudoCounts
 *
 * Duplicates the #NcClusterPseudoCounts object setting the same values of the original propertities.
 *
 * Returns: (transfer full): A new #NcClusterPseudoCounts.
   */
NcClusterPseudoCounts *
nc_cluster_pseudo_counts_copy (NcClusterPseudoCounts *cpc)
{
  return nc_cluster_pseudo_counts_new (cpc->mfp);
}

/**
 * nc_cluster_pseudo_counts_ref:
 * @cpc: a #NcClusterPseudoCounts
 *
 * Increases the reference count of @cpc by one.
 *
 * Returns: (transfer full): @cpc
 */
NcClusterPseudoCounts *
nc_cluster_pseudo_counts_ref (NcClusterPseudoCounts *cpc)
{
  return g_object_ref (cpc);
}

/**
 * nc_cluster_pseudo_counts_free:
 * @cpc: a #NcClusterPseudoCounts
 *
 * Atomically decreases the reference count of @cpc by one. If the reference count drops to 0,
 * all memory allocated by @cpc is released.
 *  
 */
void
nc_cluster_pseudo_counts_free (NcClusterPseudoCounts *cpc)
{
  g_object_unref (cpc);
}

/**
 * nc_cluster_pseudo_counts_clear:
 * @cpc: a #NcClusterPseudoCounts
 *
 * The reference count of @cpc is decreased and the pointer is set to NULL.
 *
 */
void
nc_cluster_pseudo_counts_clear (NcClusterPseudoCounts **cpc)
{
  g_clear_object (cpc);
}

static void
_nc_cluster_pseudo_counts_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterPseudoCounts *cpc = NC_CLUSTER_PSEUDO_COUNTS (object);
  g_return_if_fail (NC_IS_CLUSTER_PSEUDO_COUNTS (object));

  switch (prop_id)
  {
    case PROP_MASS_FUNCTION:
      cpc->mfp = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_cluster_pseudo_counts_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterPseudoCounts *cpc = NC_CLUSTER_PSEUDO_COUNTS (object);
  g_return_if_fail (NC_IS_CLUSTER_PSEUDO_COUNTS (object));

  switch (prop_id)
  {
    case PROP_MASS_FUNCTION:
      g_value_set_object (value, cpc->mfp);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_cluster_pseudo_counts_init (NcClusterPseudoCounts *cpc)
{
  cpc->mfp = NULL;
  cpc->workz = g_new0 (gdouble, LM_BC_DER_WORKSZ (3, 5));
}

static void
_nc_cluster_pseudo_counts_dispose (GObject *object)
{
  NcClusterPseudoCounts *cpc = NC_CLUSTER_PSEUDO_COUNTS (object);
  
  nc_mass_function_clear (&cpc->mfp);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_pseudo_counts_parent_class)->dispose (object);
}

static void
_nc_cluster_pseudo_counts_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_cluster_pseudo_counts_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_cluster_pseudo_counts, NC_TYPE_CLUSTER_PSEUDO_COUNTS);

static void
nc_cluster_pseudo_counts_class_init (NcClusterPseudoCountsClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_cluster_pseudo_counts_set_property;
  model_class->get_property = &_nc_cluster_pseudo_counts_get_property;
  object_class->dispose     = &_nc_cluster_pseudo_counts_dispose;
  object_class->finalize    = &_nc_cluster_pseudo_counts_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy Cluster observable: pseudo number counts", "PseudoClusterCounts");
  ncm_model_class_add_params (model_class, NC_CLUSTER_PSEUDO_COUNTS_SPARAM_LEN, 0, PROP_SIZE);
  
  /**
   * NcClusterPseudoCounts:mass-function:
   *
   * FIXME
   * 
   */
  g_object_class_install_property (object_class,
                                   PROP_MASS_FUNCTION,
                                   g_param_spec_object ("mass-function",
                                                        NULL,
                                                        "Mass Function",
                                                        NC_TYPE_MASS_FUNCTION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterPseudoCounts:lnMcut:
   * 
   * Logarithm base e of the lower mass cut-off, $\ln (M_{CUT}) \in [12.0 \ln(10), 16.0 \ln(10)]$.
   *  
   */
  /**
   * NcClusterPseudoCounts:lnMcut-fit:
   * 
   * FIXME
   *  
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_PSEUDO_COUNTS_LNMCUT, "\\ln(M_{CUT})", "lnMcut",
                              12.0 * M_LN10, 16.0 * M_LN10, 2.0,
                              NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_LNMCUT,
                              NCM_PARAM_TYPE_FIXED);
  /**
   * NcClusterPseudoCounts:sigma_Mcut:
   * 
   * Standard deviation of the selection function, $\sigma_{CUT} \in [0.1, 0.9]$.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_PSEUDO_COUNTS_SD_MCUT, "\\sigma_{MCUT}", "sigma_Mcut",
                              1.0e-2,  0.9, 2.0e-2,
                              NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_SD_MCUT,
                              NCM_PARAM_TYPE_FIXED);
  /**
   * NcClusterPseudoCounts:zmin:
   * 
   * Range: $z_{min} \in [10^{-5}, 2.0]$
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_PSEUDO_COUNTS_ZMIN, "z_{min}", "zmin",
                              1e-5,  2.0, 1.0e-2,
                              NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_ZMIN,
                              NCM_PARAM_TYPE_FIXED);
  /**
   * NcClusterPseudoCounts:Deltaz:
   * 
   * Maximum redsift is $z_{max} = z_{min} + \Delta z$. Range: $\Delta z \in [0.1, 2.0]$.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_PSEUDO_COUNTS_DELTAZ, "\\delta{}z", "Deltaz",
                              1e-1,  2.0, 1.0e-2,
                              NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_DELTAZ,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  ncm_mset_model_register_id (model_class,
                              "NcClusterPseudoCounts",
                              "Galaxy cluster observable: pseudo number counts.",
                              NULL);
}

//////////////////////////////////////////////////////////////////////////////

typedef struct _integrand_data
{
  NcClusterPseudoCounts *cpc;
  NcClusterMass *clusterm;
  NcHICosmo *cosmo;
  gdouble z;
  const gdouble *Mobs;
  const gdouble *Mobs_params;
  gdouble peak[3];
  gdouble func_peak;
  gdouble lnM_M0;
  gdouble lnM0;
} integrand_data;

static gdouble
_selection_function (NcClusterPseudoCounts *cpc,  gdouble lnM500)
{
  const gdouble sqrt2_sdcut = sqrt(2) * SD_MCUT;
  const gdouble difM = lnM500 - LNMCUT;

  if (difM < 0.0)
    return 0.5 * erfc (fabs(difM) / sqrt2_sdcut);
  else
    return 0.5 * (1.0 + erf (difM / sqrt2_sdcut));
}

static gdouble
_Ndet_wout_volume_integrand (gdouble lnM500, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterPseudoCounts *cpc = data->cpc;
  gdouble sf = _selection_function (cpc,  lnM500);
  gdouble mf = nc_mass_function_dn_dlnm (cpc->mfp, data->cosmo, lnM500, data->z);

  gdouble result = sf * mf;
  //printf("M = %.5g sf = %.5g mfp = %.5g result = %.5g\n", exp(lnM500), sf, mf, result);
  
  return result;
}

/**
 * nc_cluster_pseudo_counts_ndet_no_z_integral:
 * @cpc: a #NcClusterPseudoCounts
 * @cosmo: a #NcHICosmo 
 * @z: redshift
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble 
nc_cluster_pseudo_counts_ndet_no_z_integral (NcClusterPseudoCounts *cpc, NcHICosmo *cosmo, gdouble z)
{
  integrand_data data;
  gdouble P, err;
  gdouble lnM_min, lnM_max;
  gsl_function F;
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  //gdouble dV_dzdOmega; 
  
  data.cpc =   cpc;
  data.cosmo = cosmo;
  data.z = z;

  //dV_dzdOmega = nc_mass_function_dv_dzdomega (data.cpc->mfp, data.cosmo, z);
  
  F.function = _Ndet_wout_volume_integrand;
  F.params = &data;

  lnM_min= log (1.0e12);
  lnM_max = log (1.0e16);
  gsl_integration_qag (&F, lnM_min, lnM_max, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &P, &err);

  //printf ("Ndet wout z integral = %.8g\n", P);
  return P; // * dV_dzdOmega;
}

static gdouble
_Ndet_integrand (gdouble lnM500, gdouble z, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterPseudoCounts *cpc = data->cpc;
  gdouble sf = _selection_function (cpc,  lnM500);
  gdouble mf = nc_mass_function_dn_dlnm (cpc->mfp, data->cosmo, lnM500, z);
  gdouble dV_dzdOmega = nc_mass_function_dv_dzdomega (cpc->mfp, data->cosmo, z);

  gdouble result = sf * mf * dV_dzdOmega;
  //printf("M = %.5g sf = %.5g mfp = %.5g dV = %.5g result = %.5g\n", exp(lnM500), sf, mf, dV_dzdOmega, result);
  
  return result;
}

/**
 * nc_cluster_pseudo_counts_ndet:
 * @cpc: a #NcClusterPseudoCounts
 * @cosmo: a #NcHICosmo 
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble 
nc_cluster_pseudo_counts_ndet (NcClusterPseudoCounts *cpc, NcHICosmo *cosmo)
{
  integrand_data data;
  gdouble P, err;
  gdouble lnM_min, lnM_max, z_min, z_max;
  NcmIntegrand2dim integ;

  data.cpc =   cpc;
  data.cosmo = cosmo;
  
  integ.f = _Ndet_integrand;
  integ.userdata = &data;

  lnM_min= log (1.0e12);
  lnM_max = log (1.0e16);
  z_min = ZMIN;
  z_max = z_min + fabs (DELTAZ);
  ncm_integrate_2dim (&integ, lnM_min, z_min, lnM_max, z_max, NCM_DEFAULT_PRECISION, 0.0, &P, &err);

  return P;
}

static gdouble
_posterior_numerator_integrand (gdouble lnM, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterPseudoCounts *cpc = data->cpc;
  gdouble sf = _selection_function (cpc,  lnM);
  if (sf == 0.0)
    return exp (-200.0);
  else
  {
    gdouble mf = nc_mass_function_dn_dlnm (cpc->mfp, data->cosmo, lnM, data->z);
    //gdouble dV_dzdOmega = nc_mass_function_dv_dzdomega (cpc->mfp, data->cosmo, data->z);
    gdouble pdf_Mobs_Mtrue = nc_cluster_mass_p (data->clusterm, data->cosmo, lnM, data->z, data->Mobs, data->Mobs_params);

    gdouble result = sf * mf * pdf_Mobs_Mtrue + exp (-200.0); //* dV_dzdOmega;

    //printf ("===> %12.8g %12.8g %12.8g %12.8g %12.8g\n", lnM, sf, mf, pdf_Mobs_Mtrue, result);

    return result;    
  }  
}

/**
 * nc_cluster_pseudo_counts_posterior_numerator:
 * @cpc: a #NcClusterPseudoCounts
 * @clusterm: a #NcClusterMass
 * @cosmo: a #NcHICosmo 
 * @z: spectroscopic redshift
 * @Mobs: (array) (element-type double): logarithm base e of the observed mass
 * @Mobs_params: (array) (element-type double): observed mass paramaters
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_cluster_pseudo_counts_posterior_numerator (NcClusterPseudoCounts *cpc, NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble z, const gdouble *Mobs, const gdouble *Mobs_params)
{
  integrand_data data;
  gdouble P, err;
  gdouble dV_dzdOmega;

  if (z < ZMIN || z > (ZMIN + DELTAZ))
    return 0.0;
  else
  {
    gsl_function F;
    gsl_integration_workspace **w = ncm_integral_get_workspace ();

    data.cpc = cpc;
    data.clusterm = clusterm;
    data.cosmo = cosmo;
    data.z = z;
    data.Mobs = Mobs;
    data.Mobs_params = Mobs_params;

    dV_dzdOmega = nc_mass_function_dv_dzdomega (data.cpc->mfp, data.cosmo, data.z);

    F.function = &_posterior_numerator_integrand;
    F.params = &data;

    {
      gdouble lnM_min, lnM_max;
      lnM_min = log (1.0e12); //32.292; //    
      lnM_max = log (1.0e16); //34.561; //log (1.0e16);     
      gsl_integration_qag (&F, lnM_min, lnM_max, 0.0, 1.0e2 * NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &P, &err);
    }

    //printf ("numerator = %.8g err = %.8g dV = %.8g result = %.8g\n", P, err, dV_dzdOmega, P * dV_dzdOmega);
    ncm_memory_pool_return (w);
    return P * dV_dzdOmega;
  }
}

/* 3-dimensional integration: real, SZ and Lensing masses*/

static void
_nc_cluster_pseudo_counts_levmar_f (gdouble *p, gdouble *hx, gint m, gint n, gpointer adata)
{
  integrand_data *data = (integrand_data *) adata;
  NcClusterMass *clusterm = data->clusterm;
  g_assert (NC_IS_CLUSTER_MASS_PLCL (clusterm));
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (clusterm);
  const gint m1 = m;
  const gint n1 = n;
  
  //printf ("m = %d n = %d pa = %.5g p1 = %.5g\n", m1, n1, p[0], p[1]);
  //printf ("lnMtrue = %.5g mf = %.5g sf = %.5g\n", lnM_true, mf, sf); 
  //printf ("lnMsz/M0 = %.5g lnMl/M0 = %.5g\n", lnMsz_M0, lnMl_M0); 

  nc_cluster_mass_plcl_levmar_f_new_variables (p, hx, m1, n1, mszl, data->lnM_M0, data->Mobs, data->Mobs_params);
  //printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n", p[0], p[1], hx[0], hx[1], hx[2], hx[3]);
  //printf ("levmarf h0 = %.5g h1 = %.5g h2 = %.5g h3 = %.5g\n", hx[0], hx[1], hx[2], hx[3]);
}

static void
_nc_cluster_pseudo_counts_levmar_J (gdouble *p, gdouble *J, gint m, gint n, gpointer adata)
{
  integrand_data *data = (integrand_data *) adata;
  NcClusterMass *clusterm = data->clusterm;
  g_assert (NC_IS_CLUSTER_MASS_PLCL (clusterm));
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (clusterm); 
  const gint m1 = m;
  const gint n1 = n;
  
  //printf ("m = %d n = %d pa = %.5g p1 = %.5g\n", m1, n1, p[0], p[1]);
  //printf ("lnMtrue = %.5g mf = %.5g sf = %.5g\n", lnM_true, mf, sf); 
  //printf ("lnMsz/M0 = %.5g lnMl/M0 = %.5g\n", lnMsz_M0, lnMl_M0); 

  nc_cluster_mass_plcl_levmar_J_new_variables (p, J, m1, n1, mszl, data->lnM_M0, data->Mobs, data->Mobs_params);

  //printf ("levmarJ J0 = %.5g J1 = %.5g J2 = %.5g J3 = %.5g J4 = %.5g J5 = %.5g J6 = %.5g J7 = %.5g\n", J[0], J[1], J[2], J[3], J[4], J[5], J[6], J[7]);
}

static void
peakfinder (gdouble lnM_M0, gdouble p0[], const gint *ndim, const gdouble bounds[], gint *n, gdouble x[], void *userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterPseudoCounts *cpc = data->cpc;
  NcClusterMass *clusterm = data->clusterm;
  gdouble lb[] = {bounds[0], bounds[2]};
  gdouble ub[] = {bounds[1], bounds[3]};
  gdouble info[LM_INFO_SZ];
  gdouble opts[LM_OPTS_SZ];
  gint ret;

  p0[0] = log (data->Mobs[NC_CLUSTER_MASS_PLCL_MPL]);
  p0[1] = log (data->Mobs[NC_CLUSTER_MASS_PLCL_MCL]);
  
  data->lnM_M0 = lnM_M0;
  
  g_assert (NC_IS_CLUSTER_MASS_PLCL (clusterm));
  opts[0] = LM_INIT_MU; 
  opts[1] = 1.0e-7; 
  opts[2] = 1.0e-7;
  opts[3] = 1.0e-7;
  
  p0[0] = GSL_MAX (p0[0], lb[0]);
  p0[1] = GSL_MAX (p0[1], lb[1]);
  p0[0] = GSL_MIN (p0[0], ub[0]);
  p0[1] = GSL_MIN (p0[1], ub[1]);

  //printf ("p0[0] = %.5g p0[1] = %.5g lb[0] = %.5g lb[1] = %.5g ub[0] = %.5g ub[1] = %.5g\n", p0[0], p0[1], lb[0], lb[1], ub[0], ub[1]);
  ret = dlevmar_der (
                     &_nc_cluster_pseudo_counts_levmar_f, 
                     &_nc_cluster_pseudo_counts_levmar_J,
                     p0, NULL, 2, 4, 1.0e5, opts, info, cpc->workz, NULL, data
                     );
  
  if (ret < 0)
    g_error ("error: NcClusterPseudoCounts peakfinder function.\n");

  //printf ("Min %g %g Max %g %g\n", lb[0], lb[1], ub[0], ub[1]);
  //printf ("%g %g %g %g %g %g %g %g %g %g\n", 
  //        info[0], info[1], info[2], info[3], info[4], info[5],
   //       info[6], info[7], info[8], info[9]);

  x[0] = p0[0];
  x[1] = p0[1];
  x[2] = lnM_M0; 
  //printf ("Minimo: x0 = %.5g x1 = %.5g x2 = %.5g\n", x[0], x[1], x[2]);
  *n = 1;
}

static gdouble
_posterior_numerator_integrand_plcl (gdouble w1, gdouble w2, gdouble lnM_M0, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterPseudoCounts *cpc = data->cpc;   
  const gdouble lnM  = lnM_M0 + data->lnM0;
  const gdouble sf = _selection_function (cpc, lnM);
  gdouble res;
    
  if (sf == 0.0)
    res = exp (-200.0);
  else
  {
    const gdouble mf = nc_mass_function_dn_dlnm (cpc->mfp, data->cosmo, lnM, data->z);
    const gdouble pdf_Mobs_Mtrue = nc_cluster_mass_plcl_pdf (data->clusterm, lnM_M0, w1, w2, data->Mobs, data->Mobs_params);
    res = sf * mf * pdf_Mobs_Mtrue + exp (-200.0); 
  }
  //printf ("% 20.15g % 20.15g % 20.15g % 20.15g\n", w1, w2, lnM_M0, res);
  return res;
}

/**
 * nc_cluster_pseudo_counts_posterior_numerator_plcl:
 * @cpc: a #NcClusterPseudoCounts
 * @clusterm: a #NcClusterMass
 * @cosmo: a #NcHICosmo 
 * @z: spectroscopic redshift
 * @Mpl: Planck mass
 * @Mcl: CLASH mass
 * @sigma_pl: standard deviation of Planck mass
 * @sigma_cl: standard deviation of CLASH mass
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_cluster_pseudo_counts_posterior_numerator_plcl (NcClusterPseudoCounts *cpc, NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble z, const gdouble Mpl, const gdouble Mcl, const gdouble sigma_pl, const gdouble sigma_cl)
{
  integrand_data data;
  gdouble P, err;
  gdouble dV_dzdOmega;
  NcmIntegrand3dim integ;
  gdouble norma_p;
  const gdouble M0 = 1.0e14;
  const gdouble Mobs[] = {Mpl / M0, Mcl / M0};
  const gdouble Mobs_params[] = {sigma_pl / M0, sigma_cl / M0};
  
  if (z < ZMIN || z > (ZMIN + DELTAZ))
    return 0.0;
  else
  {
    g_assert (NC_IS_CLUSTER_MASS_PLCL (clusterm)); 

    data.cpc = cpc;
    data.clusterm = clusterm;
    data.cosmo = cosmo;
    data.z = z;
    data.Mobs = Mobs;
    data.Mobs_params = Mobs_params;
    data.lnM0 = log (M0);
    
    dV_dzdOmega = nc_mass_function_dv_dzdomega (data.cpc->mfp, data.cosmo, data.z);

    integ.f = _posterior_numerator_integrand_plcl;
    integ.userdata = &data;

    {
      gint n = 1;
      gint ndim = 3;
      gdouble a_w1, a_w2, b_w1, b_w2, a_true, b_true;
      gdouble lb[3], ub[3]; 
      gdouble bounds[ndim * 2], x1[ndim], p0[ndim];
      const gint ngiven   = 10;
      const gint ldxgiven = 3;
      gdouble x[ngiven * ldxgiven];
      gint i;
      
      a_w1 = -10.0; 
      b_w1 = 10.0;
      a_w2 = -10.0;
      b_w2 = 10.0;
      lb[2] = log (1.0e12 / M0);
      ub[2] = log (1.0e16 / M0);

      a_true = log (1.0e12 / M0);
      b_true = log (1.0e16 / M0);
      bounds[0] = a_w1;
      bounds[2] = a_w2;
      bounds[4] = a_true;
      bounds[1] = b_w1; 
      bounds[3] = b_w2;
      bounds[5] = b_true;
      
      lb[0] = -20.0;
      lb[1] = -20.0;
      ub[0] = 20.0;
      ub[1] = 20.0;

      /* initial point to find the peaks */
      p0[0] = log (data.Mobs[NC_CLUSTER_MASS_PLCL_MPL]);
      p0[1] = log (data.Mobs[NC_CLUSTER_MASS_PLCL_MCL]);

      //printf ("mpl = %.5g ml = %.5g\n", data.Mobs[NC_CLUSTER_MASS_PLCL_MPL], data.Mobs[NC_CLUSTER_MASS_PLCL_MCL]);
      
      for (i = 0; i < ngiven; i++)
      {
        gdouble lnM_M0 = log (2.0e12) + (log (7.0e15) - log (2.0e12)) * i * 1.0 / (ngiven - 1.0) - data.lnM0;
        peakfinder (lnM_M0, p0, &ndim, bounds, &n, x1, &data); 
        x[i * ldxgiven + 0] = x1[0]; /* w1 */
        x[i * ldxgiven + 1] = x1[1]; /* w2 */
        x[i * ldxgiven + 2] = x1[2]; /* True mass */
        //printf ("Picos: %d [%.5g, %.5g, %.5g] lnMcut % 20.15e\n", i, x[i * ldxgiven + 0], x[i * ldxgiven + 1], x[i * ldxgiven + 2], exp (LNMCUT));
          
      }    

      ncm_integrate_3dim_divonne (&integ, lb[0], lb[1], lb[2], ub[0], ub[1], ub[2], 1e-5, 0.0, ngiven, ldxgiven, x, &P, &err);
      
      norma_p = 4.0 * M_PI * M_PI * (Mobs_params[NC_CLUSTER_MASS_PLCL_MPL] * Mobs_params[NC_CLUSTER_MASS_PLCL_MCL]);
    }

    //printf ("Mpl = %.5g Mcl = %.5g\n", Mobs[0], Mobs[1]);
    //printf ("P = %.8g err = %.8e dV = %.8g result = %.8g\n", P / norma_p, err / P, dV_dzdOmega, P * dV_dzdOmega / norma_p);
    return P * dV_dzdOmega / norma_p;
  }
}

