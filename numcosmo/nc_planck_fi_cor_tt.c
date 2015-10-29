/***************************************************************************
 *            nc_planck_fi_cor_tt.c
 *
 *  Thu October 22 16:21:36 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_planck_fi_cor_tt.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_planck_fi_cor_tt
 * @title: NcPlanckFICorTT
 * @short_description: Planck Foreground and Instrument model for TT correlation maps
 *
 * FIXME
 *
 * If you use this object, cite [Planck 2015 results XI (2015)][XPlanckCollaboration2015a]
 * and related papers.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_planck_fi_cor_tt.h"
#include "math/ncm_priors.h"

enum
{
  PROP_0,
  PROP_CMP,
  PROP_SIZE
};

G_DEFINE_TYPE (NcPlanckFICorTT, nc_planck_fi_cor_tt, NC_TYPE_PLANCK_FI);

static void
nc_planck_fi_cor_tt_init (NcPlanckFICorTT *nc_planck_fi_cor_tt)
{
}

static void
nc_planck_fi_cor_tt_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_PLANCK_FI_COR_TT (object));

  switch (prop_id)
  {
    case PROP_CMP:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_planck_fi_cor_tt_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_PLANCK_FI_COR_TT (object));

  switch (prop_id)
  {
    case PROP_CMP:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_planck_fi_cor_tt_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_planck_fi_cor_tt_parent_class)->finalize (object);
}

static void
nc_planck_fi_cor_tt_class_init (NcPlanckFICorTTClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = nc_planck_fi_cor_tt_set_property;
  model_class->get_property = nc_planck_fi_cor_tt_get_property;
  object_class->finalize    = nc_planck_fi_cor_tt_finalize;

  ncm_model_class_set_name_nick (model_class, "Planck Foreground and Instument Model -- TT", "PlanckFI_TT");
  ncm_model_class_add_params (model_class, NC_PLANCK_FI_COR_TT_SPARAM_LEN, 0, PROP_SIZE);

  g_object_class_install_property (object_class,
                                   PROP_CMP,
                                   g_param_spec_uint ("cmp",
                                                      NULL,
                                                      "mmm",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_A_cib_217, "A^{\\mathrm{CIB}}_{217}", "A_cib_217",
                              0.0, 400.0, 5.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_A_cib_217,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_cib_index, "n^{\\mathrm{CIB}}", "cib_index",
                              -2.0, 0.0, 0.01,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_cib_index,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_xi_sz_cib, "\\xi^{\\mathrm{tSZ}\\times \\mathrm{CIB}}", "xi_sz_cib",
                              0.0, 1.0, 0.1,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_xi_sz_cib,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_A_sz, "A^{\\mathrm{tSZ}}", "A_sz",
                              0.0, 10.0, 0.5,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_A_sz,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_ps_A_100_100, "A^{\\mathrm{PS}}_{100}", "ps_A_100_100",
                              0.0, 400.0, 5.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_ps_A_100_100,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_ps_A_143_143, "A^{\\mathrm{PS}}_{143}", "ps_A_143_143",
                              0.0, 400.0, 5.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_ps_A_143_143,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_ps_A_143_217, "A^{\\mathrm{PS}}_{143\\times 217}", "ps_A_143_217",
                              0.0, 400.0, 5.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_ps_A_143_217,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_ps_A_217_217, "A^{\\mathrm{PS}}_{217}", "ps_A_217_217",
                              0.0, 400.0, 5.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_ps_A_217_217,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_ksz_norm, "A^{\\mathrm{kSZ}}", "ksz_norm",
                              0.0, 10.0, 0.5,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_ksz_norm,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_gal545_A_100, "A^{\\mathrm{dust}TT}_{100}", "gal545_A_100",
                              0.0, 50.0, 2.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_gal545_A_100,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_gal545_A_143, "A^{\\mathrm{dust}TT}_{143}", "gal545_A_143",
                              0.0, 50.0, 2.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_gal545_A_143,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_gal545_A_143_217, "A^{\\mathrm{dust}TT}_{143 \\times 217}", "gal545_A_143_217",
                              0.0, 100.0, 8.5,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_gal545_A_143_217,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_gal545_A_217, "A^{\\mathrm{dust}TT}_{217}", "gal545_A_217",
                              0.0, 400.0, 20.0,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_gal545_A_217,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_calib_100T, "c_{100}", "calib_100T",
                              0.0, 3.0, 0.001,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_calib_100T,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_calib_217T, "c_{217}", "calib_217T",
                              0.0, 3.0, 0.002,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_calib_217T,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_PLANCK_FI_COR_TT_A_planck, "y_{\\mathrm{cal}}", "A_planck",
                              0.9, 1.1, 0.0025,
                              NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL, NC_PLANCK_FI_COR_TT_DEFAULT_A_planck,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
}

/**
 * nc_planck_fi_cor_tt_add_gal_priors:
 * @lh: a #NcmLikelihood
 * @mean: a vector containing the means
 * @sigma: a vector containing the standard deviations
 *
 * Add the galaxy dust priors as described in [Planck 2015 results XI (2015)][XPlanckCollaboration2015a].
 *
 */
void
nc_planck_fi_cor_tt_add_gal_priors (NcmLikelihood *lh, NcmVector *mean, NcmVector *sigma)
{
  g_assert_cmpuint (ncm_vector_len (mean), ==, 4);
  g_assert_cmpuint (ncm_vector_len (sigma), ==, 4);

  ncm_prior_add_gaussian_data (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TT_gal545_A_100,     ncm_vector_get (mean, 0), ncm_vector_get (sigma, 0));
  ncm_prior_add_gaussian_data (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TT_gal545_A_143,     ncm_vector_get (mean, 1), ncm_vector_get (sigma, 1));
  ncm_prior_add_gaussian_data (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TT_gal545_A_143_217, ncm_vector_get (mean, 2), ncm_vector_get (sigma, 2));
  ncm_prior_add_gaussian_data (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TT_gal545_A_217,     ncm_vector_get (mean, 3), ncm_vector_get (sigma, 3));
}

/**
 * nc_planck_fi_cor_tt_add_default_gal_priors:
 * @lh: a #NcmLikelihood
 *
 * Add the galaxy dust priors as described in [Planck 2015 results XI (2015)][XPlanckCollaboration2015a].
 * It uses the default values.
 *
 */
void
nc_planck_fi_cor_tt_add_default_gal_priors (NcmLikelihood *lh)
{
  gdouble mean[4]  = {7.0, 9.0, 21.0, 80.0};
  gdouble sigma[4] = {2.0, 2.0,  8.5, 20.0};
  NcmVector *mean_vec = ncm_vector_new_data_static (mean, 4, 1);
  NcmVector *sigma_vec = ncm_vector_new_data_static (sigma, 4, 1);

  nc_planck_fi_cor_tt_add_gal_priors (lh, mean_vec, sigma_vec);

  ncm_vector_free (mean_vec);
  ncm_vector_free (sigma_vec);
}

/**
 * nc_planck_fi_cor_tt_add_calib_priors:
 * @lh: a #NcmLikelihood
 * @mean: a vector containing the means
 * @sigma: a vector containing the standard deviations
 *
 * Add the calibration priors as described in [Planck 2015 results XI (2015)][XPlanckCollaboration2015a].
 *
 */
void
nc_planck_fi_cor_tt_add_calib_priors (NcmLikelihood *lh, NcmVector *mean, NcmVector *sigma)
{
  g_assert_cmpuint (ncm_vector_len (mean), ==, 3);
  g_assert_cmpuint (ncm_vector_len (sigma), ==, 3);

  ncm_prior_add_gaussian_data (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TT_calib_100T, ncm_vector_get (mean, 0), ncm_vector_get (sigma, 0));
  ncm_prior_add_gaussian_data (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TT_calib_217T, ncm_vector_get (mean, 1), ncm_vector_get (sigma, 1));
  ncm_prior_add_gaussian_data (lh, nc_planck_fi_id (), NC_PLANCK_FI_COR_TT_A_planck, ncm_vector_get (mean, 2), ncm_vector_get (sigma, 2));
}

/**
 * nc_planck_fi_cor_tt_add_default_calib_priors:
 * @lh: a #NcmLikelihood
 *
 * Add the calibration priors as described in [Planck 2015 results XI (2015)][XPlanckCollaboration2015a].
 * It uses the default values.
 *
 */
void
nc_planck_fi_cor_tt_add_default_calib_priors (NcmLikelihood *lh)
{
  gdouble mean[3]  = {0.999, 0.995, 1.0000};
  gdouble sigma[3] = {0.001, 0.002, 0.0025};
  NcmVector *mean_vec = ncm_vector_new_data_static (mean, 3, 1);
  NcmVector *sigma_vec = ncm_vector_new_data_static (sigma, 3, 1);

  nc_planck_fi_cor_tt_add_calib_priors (lh, mean_vec, sigma_vec);

  ncm_vector_free (mean_vec);
  ncm_vector_free (sigma_vec);
}


typedef struct _NcPlanckFICorSZPrior
{
  gdouble f_tSZ;
  gdouble mean;
  gdouble sigma;
} NcPlanckFICorSZPrior;

static void
_nc_planck_fi_cor_tt_sz_prior_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcPlanckFICorSZPrior *tp = (NcPlanckFICorSZPrior *)obj;
  NcPlanckFICorTT *fi_cor_tt = NC_PLANCK_FI_COR_TT (ncm_mset_peek (mset, nc_planck_fi_id ()));
  const gdouble A_tSZ = ncm_model_orig_param_get (NCM_MODEL (fi_cor_tt), NC_PLANCK_FI_COR_TT_A_sz);
  const gdouble A_kSZ = ncm_model_orig_param_get (NCM_MODEL (fi_cor_tt), NC_PLANCK_FI_COR_TT_ksz_norm);

  NCM_UNUSED (x);

  g_assert (NC_IS_PLANCK_FI_COR_TT (fi_cor_tt));

  f[0] = (A_kSZ + tp->f_tSZ * A_tSZ - tp->mean) / tp->mean;
}

static void
_nc_planck_fi_cor_tt_sz_prior_free (NcPlanckFICorSZPrior *tp)
{
  g_slice_free (NcPlanckFICorSZPrior, tp);
}

/**
 * nc_planck_fi_cor_tt_add_sz_prior:
 * @lh: a #NcmLikelihood
 * @f_tSZ: the tSZ factor
 * @mean: the mean $\mu$
 * @sigma: the standard deviation
 *
 * Add the SZ prior as described in [Planck 2015 results XI (2015)][XPlanckCollaboration2015a].
 * The prior is given by a $\chi^2$ factor in the form $$\frac{(A^{\\mathrm{kSZ}} + f_\\mathrm{tSZ} A^{\\mathrm{tSZ}} - \\mu)^2}{\sigma^2}.$$
 *
 */
void
nc_planck_fi_cor_tt_add_sz_prior (NcmLikelihood *lh, gdouble f_tSZ, gdouble mean, gdouble sigma)
{
  NcPlanckFICorSZPrior *tp = g_slice_new (NcPlanckFICorSZPrior);
  NcmMSetFunc *prior = ncm_mset_func_new (&_nc_planck_fi_cor_tt_sz_prior_f, 0, 1, tp, (GDestroyNotify) _nc_planck_fi_cor_tt_sz_prior_free);

  tp->f_tSZ = f_tSZ;
  tp->mean  = mean;
  tp->sigma = sigma;

  ncm_likelihood_priors_add (lh, prior, FALSE);
  ncm_mset_func_free (prior);
}

/**
 * nc_planck_fi_cor_tt_add_default_sz_prior:
 * @lh: a #NcmLikelihood
 *
 * Add the SZ prior as described in [Planck 2015 results XI (2015)][XPlanckCollaboration2015a].
 * The prior is given by a $\chi^2$ factor in the form $$\frac{(A^{\\mathrm{kSZ}} + f_\\mathrm{tSZ} A^{\\mathrm{tSZ}} - \\mu)^2}{\sigma^2}.$$
 * It uses the default values.
 *
 */
void
nc_planck_fi_cor_tt_add_default_sz_prior (NcmLikelihood *lh)
{
  nc_planck_fi_cor_tt_add_sz_prior (lh, 1.6, 9.5, 3.0);
}


/**
 * nc_planck_fi_cor_tt_add_all_default_priors:
 * @lh: a #NcmLikelihood
 *
 * Adds all default priors above:
 * - nc_planck_fi_cor_tt_add_default_gal_priors()
 * - nc_planck_fi_cor_tt_add_default_calib_priors()
 * - nc_planck_fi_cor_tt_add_default_sz_prior()
 *
 */
void
nc_planck_fi_cor_tt_add_all_default_priors (NcmLikelihood *lh)
{
  nc_planck_fi_cor_tt_add_default_gal_priors (lh);
  nc_planck_fi_cor_tt_add_default_calib_priors (lh);
  nc_planck_fi_cor_tt_add_default_sz_prior (lh);
}
