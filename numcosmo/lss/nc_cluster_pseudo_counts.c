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
                              1.0e-3,  0.9, 2.0e-2,
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

  printf ("Ndet wout z integral = %.8g\n", P);
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
  gdouble mf = nc_mass_function_dn_dlnm (cpc->mfp, data->cosmo, lnM, data->z);
  //gdouble dV_dzdOmega = nc_mass_function_dv_dzdomega (cpc->mfp, data->cosmo, data->z);
  gdouble pdf_Mobs_Mtrue = nc_cluster_mass_p (data->clusterm, data->cosmo, lnM, data->z, data->Mobs, data->Mobs_params);

  gdouble result = sf * mf * pdf_Mobs_Mtrue; //* dV_dzdOmega;
  //gdouble result = mf;
  //printf ("===> %12.8g %12.8g %12.8g %12.8g %12.8g\n", exp(lnM), sf, mf, pdf_Mobs_Mtrue, result);
  
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
      lnM_min = log (1.0e12);
      lnM_max = log (1.0e16);
      gsl_integration_qag (&F, lnM_min, lnM_max, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &P, &err);
    }

    printf ("numerator = %.8g err = %.8g dV = %.8g result = %.8g\n", P, err, dV_dzdOmega, P * dV_dzdOmega);
    ncm_memory_pool_return (w);
    return P * dV_dzdOmega;
  }
}

/* 3-dimensional integration: real, SZ and Lensing masses*/

static void
_nc_cluster_pseudo_counts_levmar_f (gdouble *p, gdouble *hx, gint m, gint n, gpointer adata)
{
  integrand_data *data = (integrand_data *) adata;
  NcClusterPseudoCounts *cpc = data->cpc;
  NcClusterMass *clusterm = data->clusterm;
  g_assert (NC_IS_CLUSTER_MASS_PLCL (clusterm));
  NcClusterMassPlCL *mszl = NC_CLUSTER_MASS_PLCL (clusterm);
  const gdouble lnM_true  = p[2]; 
  const gint m1 = m - 1;
  const gint n1 = n - 1;
  gdouble sf = _selection_function (cpc,  lnM_true);
  gdouble mf = nc_mass_function_dn_dlnm (cpc->mfp, data->cosmo, lnM_true, data->z);

  //ncm_model_params_log_all (NCM_MODEL (clusterm));
  
  //printf ("levmar f 1\n");
  //printf ("m = %d n = %d pa = %.5g p1 = %.5g\n", m1, n1, p[0], p[1]);
  //printf ("mf = %.5g sf = %.5g\n", mf, sf);
  //printf ("p[0] = %.5g p[1] = %.5g rho = %.5g\n", p[0], p[1], COR); 
  //printf ("lnMsz/M0 = %.5g lnMl/M0 = %.5g\n", lnMsz_M0, lnMl_M0); 
  //printf ("%.5g %.5g\n", Msz, Ml);

  nc_cluster_mass_plcl_levmar_f (p, hx, m1, n1, mszl, lnM_true, data->Mobs, data->Mobs_params); 
  hx[4] = sqrt (-2.0 * log (sf));
  hx[5] = sqrt (-2.0 * log (mf));  
  //printf ("levmarf h0 = %.5g h1 = %.5g h2 = %.5g h3 = %.5g h4 = %.5g h5 = %.5g\n", hx[0], hx[1], hx[2], hx[3], hx[4], hx[5]);
  NCM_UNUSED (m);
}

static void
peakfinder (const gint *ndim, const gdouble bounds[], gint *n, gdouble x[], void *userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterPseudoCounts *cpc = data->cpc;
  NcClusterMass *clusterm = data->clusterm;
  gdouble p0[] = {log (data->Mobs[NC_CLUSTER_MASS_PLCL_MPL]/1.0e14), log (data->Mobs[NC_CLUSTER_MASS_PLCL_MCL]/1.0e14), LNMCUT};
  gdouble lb[] = {bounds[0], bounds[2], bounds[4]};
  gdouble ub[] = {bounds[1], bounds[3], bounds[5]};
  gdouble info[LM_INFO_SZ];
  gdouble opts[LM_OPTS_SZ];
  gint ret;

  g_assert (NC_IS_CLUSTER_MASS_PLCL (clusterm));
  opts[0] = LM_INIT_MU; 
  opts[1] = 1.0e-15; 
  opts[2] = 1.0e-15;
  opts[3] = 1.0e-20;
  
  p0[0] = GSL_MAX (p0[0], lb[0]);
  p0[1] = GSL_MAX (p0[1], lb[1]);
  p0[2] = GSL_MAX (p0[2], lb[2]);
  p0[0] = GSL_MIN (p0[0], ub[0]);
  p0[1] = GSL_MIN (p0[1], ub[1]);
  p0[2] = GSL_MIN (p0[2], ub[2]);

  //printf ("p0[0] = %.5g p0[1] = %.5g p0[2] = %.5g lb[0] = %.5g lb[1] = %.5g lb[2] = %.5g ub[0] = %.5g ub[1] = %.5g ub[2] = %.5g\n", p0[0], p0[1], p0[2], lb[0], lb[1], lb[2], ub[0], ub[1], ub[2]);
  ret = dlevmar_bc_dif (
                        &_nc_cluster_pseudo_counts_levmar_f, p0, NULL, 3, 5, lb, ub, NULL, 1.0e5, 
                        opts, info, cpc->workz, NULL, data
                        );
  if (ret < 0)
    g_error ("error: NcClusterPseudoCounts peakfinder function.\n");

//  printf ("Min %g %g Max %g %g\n", lb[0], lb[1], ub[0], ub[1]);
//  printf ("%g %g %g %g %g %g %g %g %g %g\n", 
//          info[0], info[1], info[2], info[3], info[4], info[5],
//          info[6], info[7], info[8], info[9]);

  x[0] = p0[0];
  x[1] = p0[1];
  x[2] = p0[2]; 
  //printf ("Minimo: x0 = %.5g x1 = %.5g x2 = %.5g\n", x[0], x[1], x[2]);
  *n = 1;
}

static gdouble
_function_at (gdouble *p, integrand_data *data)
{
  gdouble fp[6];
  _nc_cluster_pseudo_counts_levmar_f (p, fp, 3, 5, data);
  return (fp[0] * fp[0] + fp[1] * fp[1] + fp[2] * fp[2] + fp[3] * fp[3] + fp[4] * fp[4] + fp[5] * fp[5]) / 2.0;
}

#define PEAK_DEC (100.0 * M_LN10)

static gdouble
_function_bounds_Msz (gdouble x, void *params)
{
  integrand_data *data = (integrand_data *) params;
  gdouble p[3] = {x, data->peak[1], data->peak[2]};
  return _function_at (p, data) - (data->func_peak + PEAK_DEC);
}

static gdouble
_function_bounds_Ml (gdouble x, void *params)
{
  integrand_data *data = (integrand_data *) params;
  gdouble p[3] = {data->peak[0], x, data->peak[2]};
  return _function_at (p, data) - (data->func_peak + PEAK_DEC);
}

static gdouble
_function_bounds_Mtrue (gdouble x, void *params)
{
  integrand_data *data = (integrand_data *) params;
  gdouble p[3] = {data->peak[0], data->peak[1], x};
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

  F.function = &_function_bounds_Mtrue;
  lb[2] = _border_finder (&F, log (1.0e12), data->peak[2], prec, max_iter);
  ub[2] = _border_finder (&F, data->peak[2], log (1.0e16), prec, max_iter);

  if (FALSE)
  {
    gdouble p[3];
    p[0] = lb[0];
    p[1] = lb[1];
    p[2] = lb[2];
    printf ("% 20.15g % 20.15g % 20.15g % 20.15g\n", p[0], p[1], p[2], _function_at (p, data) - (data->func_peak + PEAK_DEC));
    p[0] = lb[0];
    p[1] = lb[1];
    p[2] = ub[2];
    printf ("% 20.15g % 20.15g % 20.15g % 20.15g\n", p[0], p[1], p[2], _function_at (p, data) - (data->func_peak + PEAK_DEC));
    p[0] = lb[0];
    p[1] = ub[1];
    p[2] = lb[2];
    printf ("% 20.15g % 20.15g % 20.15g % 20.15g\n", p[0], p[1], p[2], _function_at (p, data) - (data->func_peak + PEAK_DEC));
    p[0] = ub[0];
    p[1] = lb[1];
    p[2] = lb[2];
    printf ("% 20.15g % 20.15g % 20.15g % 20.15g\n", p[0], p[1], p[2], _function_at (p, data) - (data->func_peak + PEAK_DEC));
    p[0] = lb[0];
    p[1] = ub[1];
    p[2] = ub[2];
    printf ("% 20.15g % 20.15g % 20.15g % 20.15g\n", p[0], p[1], p[2], _function_at (p, data) - (data->func_peak + PEAK_DEC));
    p[0] = ub[0];
    p[1] = lb[1];
    p[2] = ub[2];
    printf ("% 20.15g % 20.15g % 20.15g % 20.15g\n", p[0], p[1], p[2], _function_at (p, data) - (data->func_peak + PEAK_DEC));
    p[0] = ub[0];
    p[1] = ub[1];
    p[2] = lb[2];
    printf ("% 20.15g % 20.15g % 20.15g % 20.15g\n", p[0], p[1], p[2], _function_at (p, data) - (data->func_peak + PEAK_DEC));
    p[0] = ub[0];
    p[1] = ub[1];
    p[2] = ub[2];
    printf ("% 20.15g % 20.15g % 20.15g % 20.15g\n", p[0], p[1], p[2], _function_at (p, data) - (data->func_peak + PEAK_DEC));
  }

}

static gdouble
_posterior_numerator_integrand_plcl (gdouble lnMsz, gdouble lnMl, gdouble lnM, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterPseudoCounts *cpc = data->cpc;
  gdouble sf = _selection_function (cpc,  lnM);
  gdouble mf = nc_mass_function_dn_dlnm (cpc->mfp, data->cosmo, lnM, data->z);
  //gdouble dV_dzdOmega = nc_mass_function_dv_dzdomega (cpc->mfp, data->cosmo, data->z);
  gdouble pdf_Mobs_Mtrue = nc_cluster_mass_plcl_pdf (data->clusterm, lnM, lnMsz, lnMl, data->Mobs, data->Mobs_params);
      
  gdouble result = sf * mf * pdf_Mobs_Mtrue; // * dV_dzdOmega;

  //printf ("%.8g %.8g | sf = %.8g mf = %.8g pdf = %.8g\n", exp(lnM), result, sf, mf, pdf_Mobs_Mtrue);
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
  gdouble dV_dzdOmega;
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
    
    dV_dzdOmega = nc_mass_function_dv_dzdomega (data.cpc->mfp, data.cosmo, data.z);

    integ.f = _posterior_numerator_integrand_plcl;
    integ.userdata = &data;

    {
      gint n = 1;
      gint ndim = 3;
      gdouble a_true, a_sz, a_l, b_true, b_sz, b_l, bounds[ndim * 2], x[ndim];
      gdouble lb[3], ub[3]; 
      const gint ngiven   = 1;
      const gint ldxgiven = 3;
      gdouble P2, err2;

      a_sz = -10.0; 
      b_sz = 60.0;
      a_l = -10.0;
      b_l = 60.0;
      a_true = log (1.0e12);
      b_true = log (1.0e16);
      bounds[0] = a_sz;
      bounds[2] = a_l;
      bounds[4] = a_true;
      bounds[1] = b_sz; 
      bounds[3] = b_l;
      bounds[5] = b_true;

      peakfinder (&ndim, bounds, &n, x, &data);
      data.peak[0] = x[0];
      data.peak[1] = x[1];
      data.peak[2] = x[2];
      data.func_peak = _function_at (x, &data);

      _p_bounds (lb, ub, &data);

      printf ("[%.5g, %.5g] [%.5g, %.5g] [%.5g, %.5g]\n", lb[0], ub[0], lb[1], ub[1], lb[2], ub[2]);
      //ncm_integrate_3dim (&integ, lb[0], lb[1], lb[2], ub[0], ub[1], ub[2], 1e-7, 0.0, &P3, &err3);
      //ncm_integrate_3dim (&integ, bounds[0], bounds[2], bounds[4], bounds[1], bounds[3], bounds[5], 1e-5, 0.0, &P, &err);
      
      //ncm_integrate_3dim_divonne (&integ, bounds[0], bounds[2], bounds[4], bounds[1], bounds[3], bounds[5], 1e-7, 0.0, ngiven, ldxgiven, x, &P2, &err2);
      //ncm_integrate_3dim_divonne (&integ, lb[0], lb[1], bounds[4], ub[0], ub[1], bounds[5], 1e-5, 0.0, ngiven, ldxgiven, x, &P, &err);
      ncm_integrate_3dim_divonne (&integ, lb[0], lb[1], lb[2], ub[0], ub[1], ub[2], 1e-7, 0.0, ngiven, ldxgiven, x, &P, &err);

      printf ("%10.8g %10.8g\n", P, err/P);
      //printf ("%10.8g %10.8g P/P2 = %10.8g %10.8g %10.8g \n", P, P2, P/P2, err/P, err2/P2);
    }

    printf ("Mpl = %.5g Mcl = %.5g\n", Mobs[0], Mobs[1]);
    printf ("P = %.8g err = %.8g dV = %.8g result = %.8g\n", P, err, dV_dzdOmega, P * dV_dzdOmega);
    return P * dV_dzdOmega;
  }
}

/* Integral over Plack and CLASH masses */

static gdouble
_posterior_numerator_integrand_normalization (gdouble lnMpl, gdouble lnMcl, gpointer userdata)
{
  integrand_data *data = (integrand_data *) userdata;
  NcClusterPseudoCounts *cpc = data->cpc;
  gdouble Mobs[2], result;;
  Mobs[0] = exp (lnMpl);
  Mobs[1] = exp (lnMcl);
  data->Mobs = Mobs;

  result = nc_cluster_pseudo_counts_posterior_numerator_plcl (cpc, data->clusterm, data->cosmo, data->z, data->Mobs, data->Mobs_params);

  return result;  
}

/**
 * nc_cluster_pseudo_counts_posterior_numerator_normalization:
 * @cpc: a #NcClusterPseudoCounts
 * @clusterm: a #NcClusterMass
 * @cosmo: a #NcHICosmo 
 * @z: spectroscopic redshift
 * @Mobs_params: (array) (element-type double): observed mass paramaters
 *
 * Integral over Planck and Clash masses. FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_cluster_pseudo_counts_posterior_numerator_normalization (NcClusterPseudoCounts *cpc, NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble z, const gdouble *Mobs_params)
{
  integrand_data data;
  gdouble P, err;
  NcmIntegrand2dim integ;

  if (z < ZMIN || z > (ZMIN + DELTAZ))
    return 0.0;
  else
  {
    data.cpc = cpc;
    data.clusterm = clusterm;
    data.cosmo = cosmo;
    data.z = z;
    data.Mobs_params = Mobs_params;

    integ.f = _posterior_numerator_integrand_normalization;
    integ.userdata = &data;
    
    {
      gdouble lnM_min, lnM_max;
      lnM_min = log (1.0e12);
      lnM_max = log (1.0e16);
      ncm_integrate_2dim (&integ, lnM_min, lnM_min, lnM_max, lnM_max, 1.0e-5, 0.0, &P, &err);
    }

    printf ("numerator = %.8g err = %.8g\n", P, err);
    return P;
  }
}