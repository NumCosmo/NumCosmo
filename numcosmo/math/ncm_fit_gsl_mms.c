/***************************************************************************
 *            ncm_fit_gsl_mms.c
 *
 *  Tue Jun 19 10:55:06 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:ncm_fit_gsl_mms
 * @title: NcmFitGSLMMS
 * @short_description: Best-fit finder -- GSL non-linear minimization (simplex) algorithms.
 *
 * This object implements a best-fit finder using the GSL non-linear
 * minimization (simplex) algorithms. It is a subclass of #NcmFit.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_gsl_mms.h"
#include "math/ncm_cfg.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_ALGO,
  PROP_SIZE,
};

struct _NcmFitGSLMMS
{
  /*< private >*/
  NcmFit parent_instance;
  gsl_multimin_fminimizer *mms;
  gsl_multimin_function f;
  NcmFitGSLMMSAlgos algo;
  gchar *desc;
  NcmVector *ss;
};


G_DEFINE_TYPE (NcmFitGSLMMS, ncm_fit_gsl_mms, NCM_TYPE_FIT);

static void
ncm_fit_gsl_mms_init (NcmFitGSLMMS *fit_gsl_mms)
{
  fit_gsl_mms->mms  = NULL;
  fit_gsl_mms->algo = 0;
  fit_gsl_mms->desc = NULL;
  fit_gsl_mms->ss   = NULL;
}

static void
_ncm_fit_gsl_mms_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitGSLMMS *fit_gsl_mms = NCM_FIT_GSL_MMS (object);

  g_return_if_fail (NCM_IS_FIT_GSL_MMS (object));

  switch (prop_id)
  {
    case PROP_ALGO:
    {
      if (fit_gsl_mms->mms == NULL)
        fit_gsl_mms->algo = g_value_get_enum (value);
      else if (fit_gsl_mms->f.n > 0)
        ncm_fit_gsl_mms_set_algo (fit_gsl_mms, g_value_get_enum (value));

      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_gsl_mms_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitGSLMMS *fit_gsl_mms = NCM_FIT_GSL_MMS (object);

  g_return_if_fail (NCM_IS_FIT_GSL_MMS (object));

  switch (prop_id)
  {
    case PROP_ALGO:
      g_value_set_enum (value, fit_gsl_mms->algo);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static gdouble nc_residual_multimin_f (const gsl_vector *x, gpointer p);

static void
_ncm_fit_gsl_mms_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fit_gsl_mms_parent_class)->constructed (object);
  {
    NcmFitGSLMMS *fit_gsl_mms = NCM_FIT_GSL_MMS (object);
    NcmFit *fit               = NCM_FIT (fit_gsl_mms);
    NcmFitState *fstate       = ncm_fit_peek_state (fit);
    NcmMSet *mset             = ncm_fit_peek_mset (fit);
    const guint fparam_len    = ncm_fit_state_get_fparam_len (fstate);
    guint i;

    if (fparam_len > 0)
    {
      fit_gsl_mms->f.f      = &nc_residual_multimin_f;
      fit_gsl_mms->f.n      = fparam_len;
      fit_gsl_mms->f.params = fit_gsl_mms;

      fit_gsl_mms->ss = ncm_vector_new (fparam_len);

      for (i = 0; i < ncm_mset_fparams_len (mset); i++)
      {
        gdouble pscale = ncm_mset_fparam_get_scale (mset, i);

        ncm_vector_set (fit_gsl_mms->ss, i, pscale * 1e-3);
      }

      ncm_fit_gsl_mms_set_algo (fit_gsl_mms, fit_gsl_mms->algo);
    }
  }
}

static void
_ncm_fit_gsl_mms_dispose (GObject *object)
{
  NcmFitGSLMMS *fit_gsl_mms = NCM_FIT_GSL_MMS (object);

  ncm_vector_clear (&fit_gsl_mms->ss);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_gsl_mms_parent_class)->dispose (object);
}

static void
_ncm_fit_gsl_mms_finalize (GObject *object)
{
  NcmFitGSLMMS *fit_gsl_mms = NCM_FIT_GSL_MMS (object);

  g_clear_pointer (&fit_gsl_mms->mms, gsl_multimin_fminimizer_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_gsl_mms_parent_class)->finalize (object);
}

static NcmFit *_ncm_fit_gsl_mms_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);
static void _ncm_fit_gsl_mms_reset (NcmFit *fit);
static gboolean _ncm_fit_gsl_mms_run (NcmFit *fit, NcmFitRunMsgs mtype);
static const gchar *_ncm_fit_gsl_mms_get_desc (NcmFit *fit);

static void
ncm_fit_gsl_mms_class_init (NcmFitGSLMMSClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmFitClass *fit_class     = NCM_FIT_CLASS (klass);

  object_class->set_property = &_ncm_fit_gsl_mms_set_property;
  object_class->get_property = &_ncm_fit_gsl_mms_get_property;
  object_class->constructed  = &_ncm_fit_gsl_mms_constructed;
  object_class->dispose      = &_ncm_fit_gsl_mms_dispose;
  object_class->finalize     = &_ncm_fit_gsl_mms_finalize;

  g_object_class_install_property (object_class,
                                   PROP_ALGO,
                                   g_param_spec_enum ("algorithm",
                                                      NULL,
                                                      "GSL multidimensional minimization algorithm [simplex]",
                                                      NCM_TYPE_FIT_GSLMMS_ALGOS, NCM_FIT_GSL_MMS_NMSIMPLEX2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  fit_class->copy_new = &_ncm_fit_gsl_mms_copy_new;
  fit_class->reset    = &_ncm_fit_gsl_mms_reset;
  fit_class->run      = &_ncm_fit_gsl_mms_run;
  fit_class->get_desc = &_ncm_fit_gsl_mms_get_desc;
}

static NcmFit *
_ncm_fit_gsl_mms_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype)
{
  NcmFitGSLMMS *fit_gsl_mms = NCM_FIT_GSL_MMS (fit);

  return ncm_fit_gsl_mms_new (lh, mset, gtype, fit_gsl_mms->algo);
}

static void
_ncm_fit_gsl_mms_reset (NcmFit *fit)
{
  /* Chain up : start */
  NCM_FIT_CLASS (ncm_fit_gsl_mms_parent_class)->reset (fit);
  {
    NcmFitGSLMMS *fit_gsl_mms = NCM_FIT_GSL_MMS (fit);
    NcmFitState *fstate       = ncm_fit_peek_state (fit);
    const guint fparam_len    = ncm_fit_state_get_fparam_len (fstate);

    if (fit_gsl_mms->f.n != fparam_len)
    {
      g_clear_pointer (&fit_gsl_mms->mms, gsl_multimin_fminimizer_free);
      fit_gsl_mms->f.n = fparam_len;
      fit_gsl_mms->mms = NULL;

      if (fparam_len > 0)
        ncm_fit_gsl_mms_set_algo (fit_gsl_mms, fit_gsl_mms->algo);
    }
  }
}

static gboolean
_ncm_fit_gsl_mms_run (NcmFit *fit, NcmFitRunMsgs mtype)
{
  NcmFitGSLMMS *fit_gsl_mms = NCM_FIT_GSL_MMS (fit);
  NcmFitState *fstate       = ncm_fit_peek_state (fit);
  NcmMSet *mset             = ncm_fit_peek_mset (fit);
  NcmVector *fparams        = ncm_fit_state_peek_fparams (fstate);
  const gdouble prec        = ncm_fit_get_params_reltol (fit);
  const guint fparam_len    = ncm_fit_state_get_fparam_len (fstate);
  gdouble last_size         = 1e300;
  gdouble last_min          = GSL_POSINF;
  gulong still_count        = 0;
  guint restart             = 10;
  const gdouble rfac        = 0.99;
  gint status;

  if (ncm_fit_equality_constraints_len (fit) || ncm_fit_inequality_constraints_len (fit))
    g_error ("_ncm_fit_gsl_mms_run: GSL algorithms do not support constraints.");

  g_assert (fparam_len != 0);

  ncm_mset_fparams_get_vector (mset, fparams);
  gsl_multimin_fminimizer_set (fit_gsl_mms->mms, &fit_gsl_mms->f, ncm_vector_gsl (fparams), ncm_vector_gsl (fit_gsl_mms->ss));

  do {
    gdouble size;

    ncm_fit_state_add_iter (fstate, 1);
    status = gsl_multimin_fminimizer_iterate (fit_gsl_mms->mms);

    if ((ncm_fit_state_get_niter (fstate) == 1) && !gsl_finite (fit_gsl_mms->mms->fval))
    {
      ncm_fit_params_set_vector (fit, fparams);

      return FALSE;
    }

    if (status)
      if (mtype > NCM_FIT_RUN_MSGS_NONE)
        ncm_fit_log_step_error (fit, gsl_strerror (status));


    size = gsl_multimin_fminimizer_size (fit_gsl_mms->mms);

    if (size == last_size)
    {
      if (++still_count == 3)
      {
        ncm_fit_log_step_error (fit, "size do not improve [prec: %8.5e size: %8.5e]", prec, size);
        status = GSL_SUCCESS;
      }
      else
      {
        status = GSL_CONTINUE;
      }
    }
    else
    {
      still_count = 0;
      last_size   = size;
      status      = gsl_multimin_test_size (size, prec);

      if (status == GSL_ETOL)
        status = GSL_CONTINUE;
    }

    if ((restart > 0) && (status == GSL_SUCCESS))
    {
      if (fit_gsl_mms->mms->fval < (last_min * rfac))
      {
        ncm_mset_fparams_get_vector (mset, fparams);
        gsl_multimin_fminimizer_set (fit_gsl_mms->mms, &fit_gsl_mms->f, ncm_vector_gsl (fparams), ncm_vector_gsl (fit_gsl_mms->ss));
        status = GSL_CONTINUE;
        restart--;
        last_min = fit_gsl_mms->mms->fval;
      }
    }

    ncm_fit_state_set_m2lnL_curval (fstate, fit_gsl_mms->mms->fval);
    ncm_fit_log_step (fit);
  } while ((status == GSL_CONTINUE) && (ncm_fit_state_get_niter (fstate) < ncm_fit_get_maxiter (fit)));

  ncm_fit_params_set_gsl_vector (fit, fit_gsl_mms->mms->x);

  ncm_mset_fparams_get_vector (mset, fparams);
  ncm_fit_state_set_m2lnL_curval (fstate, fit_gsl_mms->mms->fval);
  ncm_fit_state_set_m2lnL_prec (fstate, last_size);

  return TRUE;
}

static gdouble
nc_residual_multimin_f (const gsl_vector *x, gpointer p)
{
  NcmFit *fit   = NCM_FIT (p);
  NcmMSet *mset = ncm_fit_peek_mset (fit);
  gdouble result;

  ncm_fit_params_set_gsl_vector (fit, x);

  if (!ncm_mset_params_valid (mset))
    return GSL_POSINF;

  ncm_fit_m2lnL_val (fit, &result);

  return result;
}

static const gchar *
_ncm_fit_gsl_mms_get_desc (NcmFit *fit)
{
  NcmFitGSLMMS *fit_gsl_mms = NCM_FIT_GSL_MMS (fit);

  if (fit_gsl_mms->desc == NULL)
    fit_gsl_mms->desc = g_strdup_printf ("GSL Multidimensional Minimization:%s",
                                         fit_gsl_mms->mms != NULL ? gsl_multimin_fminimizer_name (fit_gsl_mms->mms) :
                                         "not-set");


  return fit_gsl_mms->desc;
}

/**
 * ncm_fit_gsl_mms_new:
 * @lh: a #NcmLikelihood
 * @mset: a #NcmMSet
 * @gtype: a #NcmFitGradType
 * @algo: a #NcmFitGSLMMSAlgos
 *
 * Creates a new #NcmFitGSLMMS object from the given likelihood, model set,
 * gradient type and algorithm.
 *
 *
 * Returns: (transfer full): a new #NcmFitGSLMMS object.
 */
NcmFit *
ncm_fit_gsl_mms_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, NcmFitGSLMMSAlgos algo)
{
  return g_object_new (NCM_TYPE_FIT_GSL_MMS,
                       "likelihood", lh,
                       "mset", mset,
                       "grad-type", gtype,
                       "algorithm", algo,
                       NULL
                      );
}

/**
 * ncm_fit_gsl_mms_new_default:
 * @lh: a #NcmLikelihood
 * @mset: a #NcmMSet
 * @gtype: a #NcmFitGradType
 *
 * Creates a new #NcmFitGSLMMS object from the given likelihood, model set and
 * gradient type. The algorithm is set to the default one (#NCM_FIT_GSL_MMS_NMSIMPLEX2).
 *
 * Returns: (transfer full): a new #NcmFitGSLMMS object.
 */
NcmFit *
ncm_fit_gsl_mms_new_default (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype)
{
  return g_object_new (NCM_TYPE_FIT_GSL_MMS,
                       "likelihood", lh,
                       "mset", mset,
                       "grad-type", gtype,
                       NULL
                      );
}

/**
 * ncm_fit_gsl_mms_new_by_name:
 * @lh: a #NcmLikelihood
 * @mset: a #NcmMSet
 * @gtype: a #NcmFitGradType
 * @algo_name: a string with the name of the algorithm to be used
 *
 * Creates a new #NcmFitGSLMMS object from the given likelihood, model set,
 * gradient type and algorithm name. If @algo_name is NULL, the default
 * algorithm (#NCM_FIT_GSL_MMS_NMSIMPLEX2) is used.
 *
 * Returns: (transfer full): a new #NcmFitGSLMMS object.
 */
NcmFit *
ncm_fit_gsl_mms_new_by_name (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, gchar *algo_name)
{
  if (algo_name != NULL)
  {
    const GEnumValue *algo = ncm_cfg_get_enum_by_id_name_nick (NCM_TYPE_FIT_GSLMMS_ALGOS,
                                                               algo_name);

    if (algo == NULL)
      g_error ("ncm_fit_gsl_mms_new_by_name: algorithm %s not found.", algo_name);

    return ncm_fit_gsl_mms_new (lh, mset, gtype, algo->value);
  }
  else
  {
    return ncm_fit_gsl_mms_new_default (lh, mset, gtype);
  }
}

/**
 * ncm_fit_gsl_mms_set_algo:
 * @fit_gsl_mms: a #NcmFitGSLMMS.
 * @algo: a #gsl_mms_algorithm.
 *
 * Sets the algorithm to be used by the #NcmFitGSLMMS object.
 *
 */
void
ncm_fit_gsl_mms_set_algo (NcmFitGSLMMS *fit_gsl_mms, NcmFitGSLMMSAlgos algo)
{
  const gsl_multimin_fminimizer_type *ncm_fit_gsl_mms_algos[] = {
    gsl_multimin_fminimizer_nmsimplex2,
    gsl_multimin_fminimizer_nmsimplex,
    gsl_multimin_fminimizer_nmsimplex2rand,
  };
  NcmFit *fit         = NCM_FIT (fit_gsl_mms);
  NcmFitState *fstate = ncm_fit_peek_state (fit);

  g_assert (fit_gsl_mms->algo < NCM_FIT_GSL_MMS_NUM_ALGOS);
  g_assert (fit_gsl_mms->f.n  > 0);

  if (fit_gsl_mms->algo != algo)
  {
    fit_gsl_mms->algo = algo;

    if (fit_gsl_mms->mms != NULL)
      gsl_multimin_fminimizer_free (fit_gsl_mms->mms);

    fit_gsl_mms->mms = NULL;

    if (fit_gsl_mms->desc != NULL)
      g_free (fit_gsl_mms->desc);

    fit_gsl_mms->desc = NULL;
  }

  if (fit_gsl_mms->mms == NULL)
    fit_gsl_mms->mms = gsl_multimin_fminimizer_alloc (ncm_fit_gsl_mms_algos[fit_gsl_mms->algo],
                                                      ncm_fit_state_get_fparam_len (fstate));
}

