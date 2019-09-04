/***************************************************************************
 *            ncm_fit_levmar.c
 *
 *  Wed Feb 24 21:20:09 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_fit_levmar
 * @title: NcmFitLevmar
 * @short_description: Best-fit finder -- Levenberg-Marquardt nonlinear least squares algorithm library.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_levmar.h"
#include "math/ncm_cfg.h"
#include "ncm_enum_types.h"
#include "levmar/levmar.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_ALGO,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcmFitLevmar, ncm_fit_levmar, NCM_TYPE_FIT);

static void
ncm_fit_levmar_init (NcmFitLevmar *fit_levmar)
{
  fit_levmar->workz = NULL;
  fit_levmar->ub    = NULL;
  fit_levmar->lb    = NULL;
}

static void
_ncm_fit_levmar_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fit_levmar_parent_class)->constructed (object);
  {
    NcmFitLevmar *fit_levmar = NCM_FIT_LEVMAR (object);
    NcmFit *fit = NCM_FIT (object);

    fit_levmar->fparam_len = fit->fstate->fparam_len;

    fit_levmar->lb = ncm_vector_new (fit_levmar->fparam_len);
    fit_levmar->ub = ncm_vector_new (fit_levmar->fparam_len);
    
    fit_levmar->data_len   = fit->fstate->data_len;
    ncm_fit_levmar_set_algo (fit_levmar, fit_levmar->algo);
  }
}

static void
_ncm_fit_levmar_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitLevmar *fit_levmar = NCM_FIT_LEVMAR (object);
  g_return_if_fail (NCM_IS_FIT_LEVMAR (object));

  switch (prop_id)
  {
    case PROP_ALGO:
    {
      if (fit_levmar->workz == NULL)
        fit_levmar->algo = g_value_get_enum (value);
      else
        ncm_fit_levmar_set_algo (fit_levmar, g_value_get_enum (value));
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_levmar_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitLevmar *fit_levmar = NCM_FIT_LEVMAR (object);
  g_return_if_fail (NCM_IS_FIT_LEVMAR (object));

  switch (prop_id)
  {
    case PROP_ALGO:
    g_value_set_enum (value, fit_levmar->algo);
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_levmar_dispose (GObject *object)
{
  NcmFitLevmar *fit_levmar = NCM_FIT_LEVMAR (object);

  ncm_vector_clear (&fit_levmar->lb);
  ncm_vector_clear (&fit_levmar->ub);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_levmar_parent_class)->dispose (object);
}

static void
_ncm_fit_levmar_finalize (GObject *object)
{
  NcmFitLevmar *fit_levmar = NCM_FIT_LEVMAR (object);

  g_clear_pointer (&fit_levmar->workz, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_levmar_parent_class)->finalize (object);
}

static NcmFit *_ncm_fit_levmar_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);
static void _ncm_fit_levmar_reset (NcmFit *fit);
static gboolean _ncm_fit_levmar_run (NcmFit *fit, NcmFitRunMsgs mtype);
static const gchar *_ncm_fit_levmar_get_desc (NcmFit *fit);

static void
ncm_fit_levmar_class_init (NcmFitLevmarClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmFitClass* fit_class     = NCM_FIT_CLASS (klass);

  object_class->constructed  = &_ncm_fit_levmar_constructed;
  object_class->set_property = &_ncm_fit_levmar_set_property;
  object_class->get_property = &_ncm_fit_levmar_get_property;
  object_class->dispose      = &_ncm_fit_levmar_dispose;
  object_class->finalize     = &_ncm_fit_levmar_finalize;

  g_object_class_install_property (object_class,
                                   PROP_ALGO,
                                   g_param_spec_enum ("algorithm",
                                                      NULL,
                                                      "Levmar least squares library algorithm",
                                                      NCM_TYPE_FIT_LEVMAR_ALGOS, NCM_FIT_LEVMAR_DIF,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  fit_class->copy_new = &_ncm_fit_levmar_copy_new;
  fit_class->reset    = &_ncm_fit_levmar_reset;
  fit_class->run      = &_ncm_fit_levmar_run;
  fit_class->get_desc = &_ncm_fit_levmar_get_desc;

  fit_class->is_least_squares = TRUE;
}

static NcmFit *
_ncm_fit_levmar_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype)
{
  NcmFitLevmar *fit_levmar = NCM_FIT_LEVMAR (fit);
  return ncm_fit_levmar_new (lh, mset, gtype, fit_levmar->algo);
}

static void 
_ncm_fit_levmar_reset (NcmFit *fit)
{
  /* Chain up : start */
  NCM_FIT_CLASS (ncm_fit_levmar_parent_class)->reset (fit);
  {
    NcmFitLevmar *fit_levmar = NCM_FIT_LEVMAR (fit);
    if (fit_levmar->fparam_len != fit->fstate->fparam_len || fit_levmar->data_len != fit->fstate->data_len)
    {
      fit_levmar->fparam_len = fit->fstate->fparam_len;
      
      ncm_vector_clear (&fit_levmar->lb);
      ncm_vector_clear (&fit_levmar->ub);

      fit_levmar->lb = ncm_vector_new (fit_levmar->fparam_len);
      fit_levmar->ub = ncm_vector_new (fit_levmar->fparam_len);
      
      fit_levmar->data_len   = fit->fstate->data_len;
      ncm_fit_levmar_set_algo (fit_levmar, fit_levmar->algo);
    }
  }
}

static gboolean ncm_fit_levmar_der_run (NcmFit *fit, NcmFitRunMsgs mtype);
static gboolean ncm_fit_levmar_dif_run (NcmFit *fit, NcmFitRunMsgs mtype);
static gboolean ncm_fit_levmar_bc_der_run (NcmFit *fit, NcmFitRunMsgs mtype);
static gboolean ncm_fit_levmar_bc_dif_run (NcmFit *fit, NcmFitRunMsgs mtype);

static gboolean
_ncm_fit_levmar_run (NcmFit *fit, NcmFitRunMsgs mtype)
{
  NcmFitLevmar *fit_levmar = NCM_FIT_LEVMAR (fit);
  guint i;
  
  g_assert (fit->fstate->fparam_len != 0);

  if (ncm_fit_has_equality_constraints (fit) || ncm_fit_has_inequality_constraints (fit))
    g_error ("_ncm_fit_levmar_run: GSL algorithms do not support constraints.");

  for (i = 0; i < fit_levmar->fparam_len; i++)
  {
    ncm_vector_set (fit_levmar->lb, i, ncm_mset_fparam_get_lower_bound (fit->mset, i));
    ncm_vector_set (fit_levmar->ub, i, ncm_mset_fparam_get_upper_bound (fit->mset, i));
  }
  
  switch (fit_levmar->algo)
  {
    case NCM_FIT_LEVMAR_DER:
      return ncm_fit_levmar_der_run (fit, mtype);
      break;
    case NCM_FIT_LEVMAR_DIF:
      return ncm_fit_levmar_dif_run (fit, mtype);
      break;
    case NCM_FIT_LEVMAR_BC_DER:
      return ncm_fit_levmar_bc_der_run (fit, mtype);
      break;
    case NCM_FIT_LEVMAR_BC_DIF:
      return ncm_fit_levmar_bc_dif_run (fit, mtype);
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}

static void nc_residual_levmar_f (gdouble *p, gdouble *hx, gint m, gint n, gpointer adata);
static void nc_residual_levmar_J (gdouble *p, gdouble *j, gint m, gint n, gpointer adata);

static gboolean
ncm_fit_levmar_der_run (NcmFit *fit, NcmFitRunMsgs mtype)
{
  NcmFitLevmar *fit_levmar = NCM_FIT_LEVMAR (fit);
  gdouble info[LM_INFO_SZ];
  gdouble opts[LM_OPTS_SZ];
  gint ret;

  NCM_UNUSED (mtype);

  opts[0] = LM_INIT_MU; 
  opts[1] = 1.0e-15; 
  opts[2] = 1.0e-15;
  opts[3] = 1.0e-20;
  opts[4] = LM_DIFF_DELTA;
  
  ncm_mset_fparams_get_vector (fit->mset, fit->fstate->fparams);
  
  g_assert (ncm_vector_stride (fit->fstate->ls_f) == 1 &&
            ncm_vector_stride (fit->fstate->fparams) == 1 &&
            ncm_matrix_tda (fit->fstate->covar) == ncm_matrix_ncols (fit->fstate->covar));

  ret = dlevmar_der (
                     &nc_residual_levmar_f, &nc_residual_levmar_J,
                     ncm_vector_data (fit->fstate->fparams), NULL, fit->fstate->fparam_len, fit->fstate->data_len,
                     fit->maxiter, opts, info, fit_levmar->workz, ncm_matrix_data (fit->fstate->covar), fit
                     );
  
  if (ret < 0)
    ncm_fit_log_step_error (fit, "(%d)", ret);

  ncm_fit_ls_f (fit, fit->fstate->ls_f);
  ncm_fit_state_set_m2lnL_curval (fit->fstate, info[1]);
  ncm_fit_state_set_m2lnL_prec (fit->fstate, info[2] / info[1]);
  ncm_fit_state_set_niter (fit->fstate, info[5]);

  ncm_fit_params_set_vector (fit, fit->fstate->fparams);

  return TRUE;
}

static gboolean
ncm_fit_levmar_dif_run (NcmFit *fit, NcmFitRunMsgs mtype)
{
  NcmFitLevmar *fit_levmar = NCM_FIT_LEVMAR (fit);
  gdouble info[LM_INFO_SZ];
  gdouble opts[LM_OPTS_SZ];
  gint ret;

  NCM_UNUSED (mtype);

  opts[0] = LM_INIT_MU; 
  opts[1] = 1.0e-15; 
  opts[2] = 1.0e-15;
  opts[3] = 1.0e-20;
  opts[4] = LM_DIFF_DELTA;
  
  ncm_mset_fparams_get_vector (fit->mset, fit->fstate->fparams);

  g_assert (ncm_vector_stride (fit->fstate->ls_f) == 1 &&
            ncm_vector_stride (fit->fstate->fparams) == 1 &&
            ncm_matrix_tda (fit->fstate->covar) == ncm_matrix_ncols (fit->fstate->covar));

  ret = dlevmar_dif (
                        &nc_residual_levmar_f,
                        ncm_vector_data (fit->fstate->fparams), NULL, fit->fstate->fparam_len, fit->fstate->data_len,  
                        fit->maxiter, opts, info, fit_levmar->workz, ncm_matrix_data (fit->fstate->covar), fit
                        );

  if (ret < 0)
    ncm_fit_log_step_error (fit, "(%d)", ret);

  ncm_fit_ls_f (fit, fit->fstate->ls_f);
  ncm_fit_state_set_m2lnL_curval (fit->fstate, info[1]);
  ncm_fit_state_set_m2lnL_prec (fit->fstate, info[2] / info[1]);
  ncm_fit_state_set_niter (fit->fstate, info[5]);

  ncm_fit_params_set_vector (fit, fit->fstate->fparams);

  return TRUE;
}

static gboolean
ncm_fit_levmar_bc_der_run (NcmFit *fit, NcmFitRunMsgs mtype)
{
  NcmFitLevmar *fit_levmar = NCM_FIT_LEVMAR (fit);
  gdouble info[LM_INFO_SZ];
  gdouble opts[LM_OPTS_SZ];
  gint ret;

  NCM_UNUSED (mtype);

  opts[0] = LM_INIT_MU; 
  opts[1] = 1.0e-15; 
  opts[2] = 1.0e-15;
  opts[3] = 1.0e-20;
  opts[4] = LM_DIFF_DELTA;
  
  ncm_mset_fparams_get_vector (fit->mset, fit->fstate->fparams);
  
  g_assert (ncm_vector_stride (fit->fstate->ls_f) == 1 &&
            ncm_vector_stride (fit->fstate->fparams) == 1 &&
            ncm_matrix_tda (fit->fstate->covar) == ncm_matrix_ncols (fit->fstate->covar));

  ret = dlevmar_bc_der (
                        &nc_residual_levmar_f, &nc_residual_levmar_J,
                        ncm_vector_data (fit->fstate->fparams), NULL, fit->fstate->fparam_len, fit->fstate->data_len,
                        ncm_vector_data (fit_levmar->lb), ncm_vector_data (fit_levmar->ub), NULL,
                        fit->maxiter, opts, info, fit_levmar->workz, ncm_matrix_data (fit->fstate->covar), fit
                        );
  
  if (ret < 0)
    ncm_fit_log_step_error (fit, "(%d)", ret);

  ncm_fit_ls_f (fit, fit->fstate->ls_f);
  ncm_fit_state_set_m2lnL_curval (fit->fstate, info[1]);
  ncm_fit_state_set_m2lnL_prec (fit->fstate, info[2] / info[1]);
  ncm_fit_state_set_niter (fit->fstate, info[5]);

  ncm_fit_params_set_vector (fit, fit->fstate->fparams);

  return TRUE;
}

static gboolean
ncm_fit_levmar_bc_dif_run (NcmFit *fit, NcmFitRunMsgs mtype)
{
  NcmFitLevmar *fit_levmar = NCM_FIT_LEVMAR (fit);
  gdouble info[LM_INFO_SZ];
  gdouble opts[LM_OPTS_SZ];
  gint ret;

  NCM_UNUSED (mtype);

  opts[0] = LM_INIT_MU; 
  opts[1] = 1.0e-15; 
  opts[2] = 1.0e-15;
  opts[3] = 1.0e-20;
  opts[4] = LM_DIFF_DELTA;
  
  ncm_mset_fparams_get_vector (fit->mset, fit->fstate->fparams);

  g_assert (ncm_vector_stride (fit->fstate->ls_f) == 1 &&
            ncm_vector_stride (fit->fstate->fparams) == 1 &&
            ncm_matrix_tda (fit->fstate->covar) == ncm_matrix_ncols (fit->fstate->covar));

  ret = dlevmar_bc_dif (
                        &nc_residual_levmar_f,
                        ncm_vector_data (fit->fstate->fparams), NULL, fit->fstate->fparam_len, fit->fstate->data_len,
                        ncm_vector_data (fit_levmar->lb), ncm_vector_data (fit_levmar->ub), NULL,  
                        fit->maxiter, opts, info, fit_levmar->workz, ncm_matrix_data (fit->fstate->covar), fit
                        );

  if (ret < 0)
    ncm_fit_log_step_error (fit, "(%d)", ret);

  ncm_fit_ls_f (fit, fit->fstate->ls_f);
  ncm_fit_state_set_m2lnL_curval (fit->fstate, info[1]);
  ncm_fit_state_set_m2lnL_prec (fit->fstate, info[2] / info[1]);
  ncm_fit_state_set_niter (fit->fstate, info[5]);

  ncm_fit_params_set_vector (fit, fit->fstate->fparams);

  return TRUE;
}

static void
nc_residual_levmar_f (gdouble *p, gdouble *hx, gint m, gint n, gpointer adata)
{
  NcmFit *fit = NCM_FIT (adata);
  NcmVector *f = ncm_vector_new_data_static (hx, n, 1);

  NCM_UNUSED (m);
  
  ncm_fit_params_set_array (fit, p);
  if (!ncm_mset_params_valid (fit->mset))
    g_warning ("nc_residual_levmar_f: stepping in a invalid parameter point, continuing anyway.");
  
  ncm_fit_ls_f (fit, f);

  ncm_fit_state_set_m2lnL_curval (fit->fstate, gsl_pow_2 (ncm_vector_dnrm2 (f)));
  ncm_fit_log_step (fit);
  ncm_vector_free (f);
}

static void
nc_residual_levmar_J (gdouble *p, gdouble *j, gint m, gint n, gpointer adata)
{
  NcmFit *fit = NCM_FIT (adata);
  NcmMatrix *J = ncm_matrix_new_data_static (j, n, m);

  ncm_fit_params_set_array (fit, p);
  if (!ncm_mset_params_valid (fit->mset))
    g_warning ("nc_residual_levmar_J: stepping in a invalid parameter point, continuing anyway.");

  ncm_fit_ls_J (fit, J);

  ncm_matrix_free (J);
}

static const gchar *
_ncm_fit_levmar_get_desc (NcmFit *fit)
{
  NcmFitLevmar *fit_levmar = NCM_FIT_LEVMAR (fit);
  switch (fit_levmar->algo)
  {
    case NCM_FIT_LEVMAR_DER:
      return "Levmar:external derivatives";
      break;
    case NCM_FIT_LEVMAR_DIF:
      return "Levmar:internal derivatives";
      break;
    case NCM_FIT_LEVMAR_BC_DER:
      return "Levmar:external derivatives and box constraints";
      break;
    case NCM_FIT_LEVMAR_BC_DIF:
      return "Levmar:internal derivatives and box constraints";
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}

/**
 * ncm_fit_levmar_new:
 * @lh: FIXME
 * @mset: FIXME
 * @gtype: FIXME
 * @algo: FIXME
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcmFit *
ncm_fit_levmar_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, NcmFitLevmarAlgos algo)
{
  return g_object_new (NCM_TYPE_FIT_LEVMAR, 
                       "likelihood", lh,
                       "mset", mset,
                       "grad-type", gtype,
                       "algorithm", algo,
                       NULL
                       );
}

/**
 * ncm_fit_levmar_new_default:
 * @lh: FIXME
 * @mset: FIXME
 * @gtype: FIXME
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcmFit *
ncm_fit_levmar_new_default (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype)
{
  return g_object_new (NCM_TYPE_FIT_LEVMAR, 
                       "likelihood", lh,
                       "mset", mset,
                       "grad-type", gtype,
                       NULL
                       );
}

/**
 * ncm_fit_levmar_new_by_name:
 * @lh: FIXME
 * @mset: FIXME
 * @gtype: FIXME
 * @algo_name: FIXME
 *
 * FIXME
 * 
 * Returns: FIXME
 */
NcmFit *
ncm_fit_levmar_new_by_name (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, gchar *algo_name)
{
  if (algo_name != NULL)
  {
    const GEnumValue *algo = ncm_cfg_get_enum_by_id_name_nick (NCM_TYPE_FIT_LEVMAR_ALGOS,
                                                               algo_name);
    if (algo == NULL)
      g_error ("ncm_fit_levmar_new_by_name: algorithm %s not found.", algo_name);
    return ncm_fit_levmar_new (lh, mset, gtype, algo->value);
  }
  else
    return ncm_fit_levmar_new_default (lh, mset, gtype);
}

/**
 * ncm_fit_levmar_set_algo:
 * @fit_levmar: a #NcmFitLevmar.
 * @algo: a #levmar_algorithm.
 *
 * FIXME
 *
 */
void
ncm_fit_levmar_set_algo (NcmFitLevmar *fit_levmar, NcmFitLevmarAlgos algo)
{
  if (fit_levmar->algo != algo)
    g_clear_pointer (&fit_levmar->workz, g_free);

  if (fit_levmar->workz == NULL)
  {
    switch (fit_levmar->algo)
    {
      case NCM_FIT_LEVMAR_DER: 
        fit_levmar->workz = g_new0 (gdouble, LM_DER_WORKSZ (fit_levmar->fparam_len, fit_levmar->data_len));
        break;
      case NCM_FIT_LEVMAR_DIF:
        fit_levmar->workz = g_new0 (gdouble, LM_DIF_WORKSZ (fit_levmar->fparam_len, fit_levmar->data_len));
        break;
      case NCM_FIT_LEVMAR_BC_DER: 
        fit_levmar->workz = g_new0 (gdouble, LM_BC_DER_WORKSZ (fit_levmar->fparam_len, fit_levmar->data_len));
        break;
      case NCM_FIT_LEVMAR_BC_DIF:
        fit_levmar->workz = g_new0 (gdouble, LM_BC_DIF_WORKSZ (fit_levmar->fparam_len, fit_levmar->data_len));
        break;
      default:
        g_assert_not_reached ();
        break;
    }
  }
}
