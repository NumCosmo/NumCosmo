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
 * Atomically decrements the reference count of @cpc by one. If the reference count drops to 0,
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
 * Atomically decrements the reference count of @cpc by one. If the reference count drops to 0,
 * all memory allocated by @cpc is released.
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
  cpc->mfp      = NULL;
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
nc_cluster_pseudo_counts_finalize (GObject *object)
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

  object_class->dispose = _nc_cluster_pseudo_counts_dispose;
  object_class->finalize = nc_cluster_pseudo_counts_finalize;
  model_class->set_property = _nc_cluster_pseudo_counts_set_property;
  model_class->get_property = _nc_cluster_pseudo_counts_get_property;

  ncm_model_class_set_name_nick (model_class, "Galaxy Cluster observable: pseudo number counts", "PseudoClusterCounts");
  ncm_model_class_add_params (model_class, NC_CLUSTER_PSEUDO_COUNTS_SPARAM_LEN, 0, PROP_SIZE);
  
  /**
   * NcClusterPseudoCounts:mass-function:
   *
   * FIXME
   */
  g_object_class_install_property (object_class,
                                   PROP_MASS_FUNCTION,
                                   g_param_spec_object ("mass-function",
                                                        NULL,
                                                        "Mass Function",
                                                        NC_TYPE_MASS_FUNCTION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /*
   * Lower mass cut-off: LNMCUT.
   * Logarithm base e of the lower mass cut-off, $\ln (M_{CUT}) \in [12.0 \ln(10), 16.0 \ln(10)]$.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_PSEUDO_COUNTS_LNMCUT, "\\ln(M_{CUT})", "LNMCUT",
                              12.0 * M_LN10, 16.0 * M_LN10, 1.0,
                              NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_LNMCUT,
                              NCM_PARAM_TYPE_FIXED);
  /*
   * Standard deviation of the selection function: SD_MCUT.
   * Standard deviation of the selection function, $\sigma_{CUT} \in [0.1, 0.9]$.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_PSEUDO_COUNTS_SD_MCUT, "\\sigma_{MCUT}", "SD_MCUT",
                              0.1,  0.9, 1.0e-2,
                              NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_SD_MCUT,
                              NCM_PARAM_TYPE_FIXED);
  /*
   * Minimum redshift: ZMIN.
   * Range: $z_{min} \in [10^{-5}, 2.0]$
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_PSEUDO_COUNTS_ZMIN, "\\z_{min}", "ZMIN",
                              1e-5,  2.0, 1.0e-1,
                              NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_PSEUDO_COUNTS_DEFAULT_ZMIN,
                              NCM_PARAM_TYPE_FIXED);
  /*
   * Redshift interval size: DELTAZ.
   * Maximum redsift is $z_{max} = z_{min} + \Delta z$. Range: $z_{max} \in [0.1, 2.0]$.
   */
  ncm_model_class_set_sparam (model_class, NC_CLUSTER_PSEUDO_COUNTS_DELTAZ, "\\delta{}z", "DELTAZ",
                              1e-1,  2.0, 1.0e-1,
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
  gdouble mf = nc_mass_function_dn_dlnm (cpc->mfp, data->cosmo, lnM, data->z);
  gdouble dV_dzdOmega = nc_mass_function_dv_dzdomega (cpc->mfp, data->cosmo, data->z);
  gdouble pdf_Mobs_Mtrue = nc_cluster_mass_p (data->clusterm, data->cosmo, lnM, data->z, data->Mobs, data->Mobs_params);

  gdouble result = sf * mf * dV_dzdOmega * pdf_Mobs_Mtrue;
  
  return result;  
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

    F.function = &_posterior_numerator_integrand;
    F.params = &data;

    {
      gdouble lnM_min, lnM_max;
      lnM_min = log (1.0e12);
      lnM_max = log (1.0e16);
      gsl_integration_qag (&F, lnM_min, lnM_max, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &P, &err);
    }

    ncm_memory_pool_return (w);
    return P;
  }
}

static gdouble
_posterior_numerator_integrand_plcl (gdouble lnM, gdouble lnMsz, gdouble lnMl, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterPseudoCounts *cpc = data->cpc;
  gdouble sf = _selection_function (cpc,  lnM);
  gdouble mf = nc_mass_function_dn_dlnm (cpc->mfp, data->cosmo, lnM, data->z);
  gdouble dV_dzdOmega = nc_mass_function_dv_dzdomega (cpc->mfp, data->cosmo, data->z);
  gdouble pdf_Mobs_Mtrue = nc_cluster_mass_plcl_pdf (data->clusterm, lnM, lnMsz, lnMl, data->Mobs, data->Mobs_params);
    
  gdouble result = sf * mf * dV_dzdOmega * pdf_Mobs_Mtrue;
  
  return result;  
}

/**
 * nc_cluster_pseudo_counts_posterior_numerator_plcl:
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
nc_cluster_pseudo_counts_posterior_numerator_plcl (NcClusterPseudoCounts *cpc, NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble z, const gdouble *Mobs, const gdouble *Mobs_params)
{
  integrand_data data;
  gdouble P, err;
  NcmIntegrand3dim integ;

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

    integ.f = _posterior_numerator_integrand_plcl;
    integ.userdata = &data;

    {
      gdouble a_true, a_sz, a_l, b_true, b_sz, b_l;
      a_true = log (1.0e12);
      b_true = log (1.0e16);
      a_sz = log (fabs (Mobs[0] - 8.0 * Mobs_params[0]));
      b_sz = log (Mobs[0] + 8.0 *  Mobs_params[0]);
      a_l = log (fabs (Mobs[1] - 8.0 * Mobs_params[1]));
      b_l = log (Mobs[1] + 8.0 * Mobs_params[1]);

      ncm_integrate_3dim (&integ, a_true, a_sz, a_l, b_true, b_sz, b_l, 1e-30, 0.0, &P, &err);
    }

    /*
    {
      gdouble a_true, a_sz, a_l, b_true, b_sz, b_l;
      gdouble P1, P2, P3, P4, P5, P6, P7, P8;
      gdouble err1, err2, err3, err4, err5, err6, err7, err8;
      {
        a_sz = log (fabs (Mobs[0] - 8.0 * Mobs_params[0]));
        b_sz = log (Mobs[0]);
        a_l = log (fabs (Mobs[1] - 8.0 * Mobs_params[1]));
        b_l = log (Mobs[1]);
        {
          a_true = log (1.0e12);
          b_true = log (1.0e14);

          ncm_integrate_3dim (&integ, a_true, a_sz, a_l, b_true, b_sz, b_l, 1e-20, 0.0, &P1, &err1);
        }
        {
          a_true = log (1.0e14);
          b_true = log (1.0e16);

          ncm_integrate_3dim (&integ, a_true, a_sz, a_l, b_true, b_sz, b_l, 1e-20, 0.0, &P2, &err2);
        }
      }
      {
        a_sz = log (fabs (Mobs[0] - 8.0 * Mobs_params[0]));
        b_sz = log (Mobs[0]);
        a_l  = log (Mobs[1]);
        b_l  = log (Mobs[1] + 8.0 * Mobs_params[1]);
        {
          a_true = log (1.0e12);
          b_true = log (1.0e14);

          ncm_integrate_3dim (&integ, a_true, a_sz, a_l, b_true, b_sz, b_l, 1e-20, 0.0, &P3, &err3);
        }
        {
          a_true = log (1.0e14);
          b_true = log (1.0e16);

          ncm_integrate_3dim (&integ, a_true, a_sz, a_l, b_true, b_sz, b_l, 1e-20, 0.0, &P4, &err4);
        }
      }
      {
        a_sz = log (Mobs[0]);
        b_sz = log (Mobs[0] + 8.0 * Mobs_params[0]);
        a_l  = log (fabs (Mobs[1] - 8.0 * Mobs_params[1]));
        b_l  = log (Mobs[1]);
        {
          a_true = log (1.0e12);
          b_true = log (1.0e14);

          ncm_integrate_3dim (&integ, a_true, a_sz, a_l, b_true, b_sz, b_l, 1e-20, 0.0, &P5, &err5);
        }
        {
          a_true = log (1.0e14);
          b_true = log (1.0e16);

          ncm_integrate_3dim (&integ, a_true, a_sz, a_l, b_true, b_sz, b_l, 1e-20, 0.0, &P6, &err6);
        }
      }
      {
        a_sz = log (Mobs[0]);
        b_sz = log (Mobs[0] + 8.0 * Mobs_params[0]);
        a_l  = log (Mobs[1]);
        b_l  = log (Mobs[1] + 8.0 * Mobs_params[1]);
        {
          a_true = log (1.0e12);
          b_true = log (1.0e14);

          ncm_integrate_3dim (&integ, a_true, a_sz, a_l, b_true, b_sz, b_l, 1e-20, 0.0, &P7, &err7);
        }
        {
          a_true = log (1.0e14);
          b_true = log (1.0e16);

          ncm_integrate_3dim (&integ, a_true, a_sz, a_l, b_true, b_sz, b_l, 1e-20, 0.0, &P8, &err8);
        }
      }
      P = P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8;
      err = err1 + err2 + err3 + err4 + err5 + err6 + err7 + err8; 
      printf ("P1 = %.8g P2 = %.8g P3 = %.8g P4 = %.8g P5 = %.8g P6 = %.8g P7 = %.8g P8 = %.8g\n", P1, P2, P3, P4, P5, P6, P7, P8);
      printf ("err1 = %.8g err2 = %.8g err3 = %.8g err4 = %.8g err5 = %.8g err6 = %.8g err7 = %.8g err8 = %.8g\n", err1, err2, err3, err4, err5, err6, err7, err8);
    }
   */
    printf ("P = %.8g err = %.8g\n", P, err);
    return P;
  }
}
