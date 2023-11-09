/***************************************************************************
 *            ncm_fit_nlopt.c
 *
 *  Sat Apr  3 16:07:02 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
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
 * SECTION:ncm_fit_nlopt
 * @title: NLopt Interface Object
 * @short_description: Interface for NLopt optmization library
 *
 * A subclass of #NcmFit that uses the NLopt library to perform the
 * optimization.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#ifdef NUMCOSMO_HAVE_NLOPT

#include "math/ncm_fit_nlopt.h"
#include "math/ncm_cfg.h"
#include "ncm_fit_nlopt_enum.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <nlopt.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_ALGO,
  PROP_LOCAL_ALGO,
  PROP_SIZE,
};

struct _NcmFitNLOpt
{
  /*< private >*/
  NcmFit parent_instance;
#ifdef NUMCOSMO_HAVE_NLOPT
  nlopt_opt nlopt;
  nlopt_opt local_nlopt;
  NcmFitNloptAlgorithm nlopt_algo;
  NcmFitNloptAlgorithm local_nlopt_algo;
#endif /* NUMCOSMO_HAVE_NLOPT */
  NcmVector *lb;
  NcmVector *ub;
  NcmVector *pabs;
  NcmVector *pscale;
  gchar *desc;
  guint fparam_len;
};

G_DEFINE_TYPE (NcmFitNLOpt, ncm_fit_nlopt, NCM_TYPE_FIT);

static void
ncm_fit_nlopt_init (NcmFitNLOpt *fit_nlopt)
{
  fit_nlopt->nlopt            = NULL;
  fit_nlopt->local_nlopt      = NULL;
  fit_nlopt->nlopt_algo       = 0;
  fit_nlopt->local_nlopt_algo = 0;
  fit_nlopt->desc             = NULL;
}

static void
_ncm_fit_nlopt_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fit_nlopt_parent_class)->constructed (object);
  {
    NcmFitNLOpt *fit_nlopt = NCM_FIT_NLOPT (object);
    NcmFit *fit            = NCM_FIT (fit_nlopt);
    NcmFitState *fstate    = ncm_fit_peek_state (fit);

    fit_nlopt->fparam_len = ncm_fit_state_get_fparam_len (fstate);

    if (fit_nlopt->fparam_len > 0)
    {
      fit_nlopt->lb     = ncm_vector_new (fit_nlopt->fparam_len);
      fit_nlopt->ub     = ncm_vector_new (fit_nlopt->fparam_len);
      fit_nlopt->pabs   = ncm_vector_new (fit_nlopt->fparam_len);
      fit_nlopt->pscale = ncm_vector_new (fit_nlopt->fparam_len);
    }

    ncm_fit_nlopt_set_algo (fit_nlopt, fit_nlopt->nlopt_algo);
  }
}

static void
_ncm_fit_nlopt_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitNLOpt *fit_nlopt = NCM_FIT_NLOPT (object);

  g_return_if_fail (NCM_IS_FIT_NLOPT (object));

  switch (prop_id)
  {
    case PROP_ALGO:

      if (fit_nlopt->nlopt == NULL)
        fit_nlopt->nlopt_algo = g_value_get_enum (value);
      else
        ncm_fit_nlopt_set_algo (fit_nlopt, g_value_get_enum (value));

      break;
    case PROP_LOCAL_ALGO:
      ncm_fit_nlopt_set_local_algo (fit_nlopt, g_value_get_enum (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_fit_nlopt_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitNLOpt *fit_nlopt = NCM_FIT_NLOPT (object);

  g_return_if_fail (NCM_IS_FIT_NLOPT (object));

  switch (prop_id)
  {
    case PROP_ALGO:
      g_value_set_enum (value, fit_nlopt->nlopt_algo);
      break;
    case PROP_LOCAL_ALGO:
      g_value_set_enum (value, fit_nlopt->local_nlopt_algo);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
ncm_fit_nlopt_dispose (GObject *object)
{
  NcmFitNLOpt *fit_nlopt = NCM_FIT_NLOPT (object);

  ncm_vector_clear (&fit_nlopt->lb);
  ncm_vector_clear (&fit_nlopt->ub);
  ncm_vector_clear (&fit_nlopt->pabs);
  ncm_vector_clear (&fit_nlopt->pscale);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_nlopt_parent_class)->dispose (object);
}

static void
ncm_fit_nlopt_finalize (GObject *object)
{
  NcmFitNLOpt *fit_nlopt = NCM_FIT_NLOPT (object);

  if (fit_nlopt->nlopt != NULL)
    nlopt_destroy (fit_nlopt->nlopt);

  if (fit_nlopt->local_nlopt != NULL)
    nlopt_destroy (fit_nlopt->local_nlopt);

  if (fit_nlopt->desc != NULL)
    g_free (fit_nlopt->desc);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_nlopt_parent_class)->finalize (object);
}

static NcmFit *_ncm_fit_nlopt_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);
static void _ncm_fit_nlopt_reset (NcmFit *fit);
static gboolean _ncm_fit_nlopt_run (NcmFit *fit, NcmFitRunMsgs mtype);
static const gchar *_ncm_fit_nlopt_get_desc (NcmFit *fit);

static void
ncm_fit_nlopt_class_init (NcmFitNLOptClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmFitClass *fit_class     = NCM_FIT_CLASS (klass);

  object_class->constructed  = &_ncm_fit_nlopt_constructed;
  object_class->set_property = &_ncm_fit_nlopt_set_property;
  object_class->get_property = &_ncm_fit_nlopt_get_property;
  object_class->dispose      = &ncm_fit_nlopt_dispose;
  object_class->finalize     = &ncm_fit_nlopt_finalize;

  g_object_class_install_property (object_class,
                                   PROP_ALGO,
                                   g_param_spec_enum ("algorithm",
                                                      NULL,
                                                      "NLOpt algorithm",
                                                      NCM_TYPE_FIT_NLOPT_ALGORITHM, NLOPT_LN_NELDERMEAD,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_LOCAL_ALGO,
                                   g_param_spec_enum ("local-algorithm",
                                                      NULL,
                                                      "NLOpt local algorithm",
                                                      NCM_TYPE_FIT_NLOPT_ALGORITHM, NLOPT_LN_NELDERMEAD,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  fit_class->copy_new = &_ncm_fit_nlopt_copy_new;
  fit_class->reset    = &_ncm_fit_nlopt_reset;
  fit_class->run      = &_ncm_fit_nlopt_run;
  fit_class->get_desc = &_ncm_fit_nlopt_get_desc;
}

static NcmFit *
_ncm_fit_nlopt_copy_new (NcmFit *fit, NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype)
{
  NcmFitNLOpt *fit_nlopt = NCM_FIT_NLOPT (fit);

  if (fit_nlopt->local_nlopt_algo == 0)
    return ncm_fit_nlopt_new (lh, mset, gtype, fit_nlopt->nlopt_algo);
  else
    return ncm_fit_nlopt_local_new (lh, mset, gtype, fit_nlopt->nlopt_algo, fit_nlopt->local_nlopt_algo);
}

static void
_ncm_fit_nlopt_reset (NcmFit *fit)
{
  /* Chain up : start */
  NCM_FIT_CLASS (ncm_fit_nlopt_parent_class)->reset (fit);
  {
    NcmFitNLOpt *fit_nlopt = NCM_FIT_NLOPT (fit);
    NcmFitState *fstate    = ncm_fit_peek_state (fit);
    const guint fparam_len = ncm_fit_state_get_fparam_len (fstate);

    if (fit_nlopt->fparam_len != fparam_len)
    {
      fit_nlopt->fparam_len = fparam_len;

      ncm_vector_clear (&fit_nlopt->lb);
      ncm_vector_clear (&fit_nlopt->ub);
      ncm_vector_clear (&fit_nlopt->pabs);
      ncm_vector_clear (&fit_nlopt->pscale);

      if (fit_nlopt->fparam_len > 0)
      {
        fit_nlopt->lb     = ncm_vector_new (fit_nlopt->fparam_len);
        fit_nlopt->ub     = ncm_vector_new (fit_nlopt->fparam_len);
        fit_nlopt->pabs   = ncm_vector_new (fit_nlopt->fparam_len);
        fit_nlopt->pscale = ncm_vector_new (fit_nlopt->fparam_len);
      }

      g_clear_pointer (&fit_nlopt->nlopt, nlopt_destroy);
      g_clear_pointer (&fit_nlopt->local_nlopt, nlopt_destroy);

      ncm_fit_nlopt_set_algo (fit_nlopt, fit_nlopt->nlopt_algo);

      if (fit_nlopt->local_nlopt_algo != 0)
        ncm_fit_nlopt_set_algo (fit_nlopt, fit_nlopt->local_nlopt_algo);
    }
  }
}

typedef gdouble (*_NcmFitNLOptOldFunc) (gint n, const gdouble *x, gdouble *grad, gpointer userdata);
static gdouble _ncm_fit_nlopt_func (guint n, const gdouble *x, gdouble *grad, gpointer userdata);
static gdouble _ncm_fit_nlopt_func_constraint (guint n, const gdouble *x, gdouble *grad, gpointer userdata);

typedef struct _NcmFitNLOptConst
{
  NcmMSetFunc *func;
  NcmFit *fit;
} NcmFitNLOptConst;

static gboolean
_ncm_fit_nlopt_run (NcmFit *fit, NcmFitRunMsgs mtype)
{
  NcmFitNLOpt *fit_nlopt = NCM_FIT_NLOPT (fit);
  NcmFitState *fstate    = ncm_fit_peek_state (fit);
  NcmVector *fparams     = ncm_fit_state_peek_fparams (fstate);
  NcmMSet *mset          = ncm_fit_peek_mset (fit);
  const guint fparam_len = ncm_fit_state_get_fparam_len (fstate);
  gdouble minf           = 0.0;

  NCM_UNUSED (mtype);

  g_assert (fparam_len != 0);

  ncm_mset_fparams_get_vector (mset, fparams);

  {
    GArray *ca = g_array_new (FALSE, FALSE, sizeof (NcmFitNLOptConst));
    nlopt_result ret;
    guint i;

    for (i = 0; i < fit_nlopt->fparam_len; i++)
    {
      ncm_vector_set (fit_nlopt->lb,     i, ncm_mset_fparam_get_lower_bound (mset, i));
      ncm_vector_set (fit_nlopt->ub,     i, ncm_mset_fparam_get_upper_bound (mset, i));
      ncm_vector_set (fit_nlopt->pabs,   i, ncm_mset_fparam_get_abstol (mset, i));
      ncm_vector_set (fit_nlopt->pscale, i, ncm_mset_fparam_get_scale (mset, i));
    }

    nlopt_remove_inequality_constraints (fit_nlopt->nlopt);
    nlopt_remove_equality_constraints (fit_nlopt->nlopt);

    for (i = 0; i < ncm_fit_inequality_constraints_len (fit); i++)
    {
      NcmFitNLOptConst fc = {NULL, fit};
      gdouble tot;

      ncm_fit_get_inequality_constraint (fit, i, &fc.func, &tot);

      g_array_append_val (ca, fc);

      ret = nlopt_add_inequality_constraint (fit_nlopt->nlopt, &_ncm_fit_nlopt_func_constraint, &g_array_index (ca, NcmFitNLOptConst, ca->len - 1), tot);

      if (ret < 0)
        g_error ("_ncm_fit_nlopt_run: cannot add inequality constrain: (%d)", ret);
    }

    for (i = 0; i < ncm_fit_equality_constraints_len (fit); i++)
    {
      NcmFitNLOptConst fc = {NULL, fit};
      gdouble tot;

      ncm_fit_get_equality_constraint (fit, i, &fc.func, &tot);

      g_array_append_val (ca, fc);

      ret = nlopt_add_equality_constraint (fit_nlopt->nlopt, &_ncm_fit_nlopt_func_constraint, &g_array_index (ca, NcmFitNLOptConst, ca->len - 1), tot);

      if (ret < 0)
        g_error ("_ncm_fit_nlopt_run: cannot add equality constrain: (%d)", ret);
    }

    nlopt_set_min_objective (fit_nlopt->nlopt, &_ncm_fit_nlopt_func, fit);
    nlopt_set_lower_bounds (fit_nlopt->nlopt, ncm_vector_data (fit_nlopt->lb));
    nlopt_set_upper_bounds (fit_nlopt->nlopt, ncm_vector_data (fit_nlopt->ub));

    nlopt_set_ftol_rel (fit_nlopt->nlopt, ncm_fit_get_m2lnL_reltol (fit));
    nlopt_set_ftol_abs (fit_nlopt->nlopt, ncm_fit_get_m2lnL_abstol (fit));
    nlopt_set_xtol_rel (fit_nlopt->nlopt, ncm_fit_get_params_reltol (fit));
    nlopt_set_xtol_abs (fit_nlopt->nlopt, ncm_vector_data (fit_nlopt->pabs));
    nlopt_set_maxeval (fit_nlopt->nlopt, ncm_fit_get_maxiter (fit));
    nlopt_set_initial_step (fit_nlopt->nlopt, ncm_vector_data (fit_nlopt->pscale));

    if (fit_nlopt->local_nlopt != NULL)
    {
      nlopt_set_ftol_rel (fit_nlopt->local_nlopt, ncm_fit_get_m2lnL_reltol (fit));
      nlopt_set_ftol_abs (fit_nlopt->local_nlopt, ncm_fit_get_m2lnL_abstol (fit));
      nlopt_set_xtol_rel (fit_nlopt->local_nlopt, ncm_fit_get_params_reltol (fit));
      nlopt_set_xtol_abs (fit_nlopt->local_nlopt, ncm_vector_data (fit_nlopt->pabs));
      nlopt_set_maxeval (fit_nlopt->local_nlopt, ncm_fit_get_maxiter (fit));
      ret = nlopt_set_local_optimizer (fit_nlopt->nlopt, fit_nlopt->local_nlopt);

      if (ret < 0)
        g_error ("_ncm_fit_nlopt_run[local_nlopt]: (%d)", ret);
    }

    ret = nlopt_optimize (fit_nlopt->nlopt, ncm_vector_data (fparams), &minf);

    ncm_fit_state_set_m2lnL_prec (fstate,
                                  GSL_MAX (nlopt_get_ftol_rel (fit_nlopt->nlopt),
                                           nlopt_get_ftol_abs (fit_nlopt->nlopt) / minf)
                                 );
    ncm_fit_state_set_params_prec (fstate, nlopt_get_xtol_rel (fit_nlopt->nlopt));

    if (ret < 0)
    {
      switch (ret)
      {
        case -2:
          ncm_fit_log_step_error (fit, "Algorithm not supported or inconsistent parameters", ret);
          break;
        default:
          ncm_fit_log_step_error (fit, "(%d)", ret);
          break;
      }
    }

    g_array_unref (ca);
  }

  {
    const gdouble m2lnL_prec = ncm_fit_state_get_m2lnL_prec (fstate);
    gdouble m2lnL            = 0.0;

    ncm_fit_params_set_vector (fit, fparams);
    ncm_fit_m2lnL_val (fit, &m2lnL);

    if (ncm_cmp (m2lnL, minf, m2lnL_prec, 0.0) != 0)
      g_warning ("_ncm_fit_nlopt_run: algorithm minimum differs from evaluated m2lnL % 22.15g != % 22.15g (prec = %e)\n",
                 m2lnL, minf, m2lnL_prec);

    ncm_fit_state_set_m2lnL_curval (fstate, minf);
  }

  return TRUE;
}

static gdouble
_ncm_fit_nlopt_func (guint n, const gdouble *x, gdouble *grad, gpointer userdata)
{
  NcmFit *fit         = NCM_FIT (userdata);
  NcmFitState *fstate = ncm_fit_peek_state (fit);
  NcmMSet *mset       = ncm_fit_peek_mset (fit);
  gdouble m2lnL;
  guint i;

  ncm_fit_state_add_iter (fstate, 1);

  for (i = 0; i < n; i++)
  {
    if (!gsl_finite (x[i]))
      return GSL_POSINF;
  }

  ncm_fit_params_set_array (fit, x);

  if (!ncm_mset_params_valid (mset))
    return GSL_POSINF;

  if (grad != NULL)
  {
    NcmVector *gradv = ncm_vector_new_data_static (grad, n, 1);

    ncm_fit_m2lnL_val_grad (fit, &m2lnL, gradv);
    ncm_vector_free (gradv);
  }
  else
  {
    ncm_fit_m2lnL_val (fit, &m2lnL);
  }

  ncm_fit_state_set_m2lnL_curval (fstate, m2lnL);
  ncm_fit_log_step (fit);

  return m2lnL;
}

static gdouble
_ncm_fit_nlopt_func_constraint (guint n, const gdouble *x, gdouble *grad, gpointer userdata)
{
  NcmFitNLOptConst *fc = (NcmFitNLOptConst *) userdata;
  NcmFit *fit          = fc->fit;
  NcmMSet *mset        = ncm_fit_peek_mset (fit);
  NcmFitState *fstate  = ncm_fit_peek_state (fit);
  gdouble constraint;

  ncm_fit_state_add_func_eval (fstate, 1);
  ncm_fit_params_set_array (fit, x);

  if (!ncm_mset_params_valid (mset))
    return GSL_NAN;

  constraint = ncm_mset_func_eval0 (fc->func, mset);

  if (grad != NULL)
  {
    NcmVector *gradv = ncm_vector_new_data_static (grad, n, 1);

    ncm_mset_func_numdiff_fparams (fc->func, mset, NULL, &gradv);
    ncm_vector_scale (gradv, 2.0 * constraint);
    ncm_vector_free (gradv);
  }

  return constraint * constraint;
}

static const gchar *
_ncm_fit_nlopt_get_desc (NcmFit *fit)
{
  NcmFitNLOpt *fit_nlopt = NCM_FIT_NLOPT (fit);

  if (fit_nlopt->desc == NULL)
  {
    GEnumClass *enum_class = g_type_class_ref (NCM_TYPE_FIT_NLOPT_ALGORITHM);
    GEnumValue *res        = g_enum_get_value (enum_class, fit_nlopt->nlopt_algo);

    if (fit_nlopt->local_nlopt_algo != 0)
    {
      GEnumValue *local_res = g_enum_get_value (enum_class, fit_nlopt->local_nlopt_algo);

      fit_nlopt->desc = g_strdup_printf ("NLOpt:%s:%s", res->value_nick, local_res->value_nick);
    }
    else
    {
      fit_nlopt->desc = g_strdup_printf ("NLOpt:%s", res->value_nick);
    }

    g_type_class_unref (enum_class);
  }

  return fit_nlopt->desc;
}

/**
 * ncm_fit_nlopt_new:
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
ncm_fit_nlopt_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, NcmFitNloptAlgorithm algo)
{
  return g_object_new (NCM_TYPE_FIT_NLOPT,
                       "likelihood", lh,
                       "mset", mset,
                       "grad-type", gtype,
                       "algorithm", algo,
                       NULL
                      );
}

/**
 * ncm_fit_nlopt_local_new:
 * @lh: FIXME
 * @mset: FIXME
 * @gtype: FIXME
 * @algo: FIXME
 * @local_algo: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmFit *
ncm_fit_nlopt_local_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, NcmFitNloptAlgorithm algo, NcmFitNloptAlgorithm local_algo)
{
  return g_object_new (NCM_TYPE_FIT_NLOPT,
                       "likelihood", lh,
                       "mset", mset,
                       "grad-type", gtype,
                       "algorithm", algo,
                       "local-algorithm", local_algo,
                       NULL
                      );
}

/**
 * ncm_fit_nlopt_new_default:
 * @lh: FIXME
 * @mset: FIXME
 * @gtype: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmFit *
ncm_fit_nlopt_new_default (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype)
{
  return g_object_new (NCM_TYPE_FIT_NLOPT,
                       "likelihood", lh,
                       "mset", mset,
                       "grad-type", gtype,
                       NULL
                      );
}

/**
 * ncm_fit_nlopt_new_by_name:
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
ncm_fit_nlopt_new_by_name (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, gchar *algo_name)
{
  if (algo_name != NULL)
  {
    gchar **algo_names = g_strsplit (algo_name, ":", 2);
    guint nnames       = g_strv_length (algo_names);

    if (nnames == 1)
    {
      const GEnumValue *algo = ncm_cfg_get_enum_by_id_name_nick (NCM_TYPE_FIT_NLOPT_ALGORITHM,
                                                                 algo_name);

      if (algo == NULL)
        g_error ("ncm_fit_nlopt_new_by_name: algorithm %s not found.", algo_name);

      g_strfreev (algo_names);

      return ncm_fit_nlopt_new (lh, mset, gtype, algo->value);
    }
    else if (nnames == 2)
    {
      const GEnumValue *algo = ncm_cfg_get_enum_by_id_name_nick (NCM_TYPE_FIT_NLOPT_ALGORITHM,
                                                                 algo_names[0]);
      const GEnumValue *local_algo = ncm_cfg_get_enum_by_id_name_nick (NCM_TYPE_FIT_NLOPT_ALGORITHM,
                                                                       algo_names[1]);

      if (algo == NULL)
        g_error ("ncm_fit_nlopt_new_by_name: algorithm %s not found.", algo_names[0]);

      if (local_algo == NULL)
        g_error ("ncm_fit_nlopt_new_by_name: algorithm %s not found.", algo_names[1]);

      g_strfreev (algo_names);

      return ncm_fit_nlopt_local_new (lh, mset, gtype, algo->value, local_algo->value);
    }
    else
    {
      g_error ("ncm_fit_nlopt_new_by_name: cannot parse algorithm name ``%s''.", algo_name);

      return NULL;
    }
  }
  else
  {
    return ncm_fit_nlopt_new_default (lh, mset, gtype);
  }
}

/**
 * ncm_fit_nlopt_set_algo: (skip)
 * @fit_nlopt: a #NcmFitNLOpt.
 * @algo: a #NcmFitNloptAlgorithm.
 *
 * FIXME
 *
 */
void
ncm_fit_nlopt_set_algo (NcmFitNLOpt *fit_nlopt, NcmFitNloptAlgorithm algo)
{
  NcmFit *fit            = NCM_FIT (fit_nlopt);
  NcmFitState *fstate    = ncm_fit_peek_state (fit);
  const guint fparam_len = ncm_fit_state_get_fparam_len (fstate);

  if (fit_nlopt->nlopt_algo != algo)
    g_clear_pointer (&fit_nlopt->nlopt, nlopt_destroy);

  if (fit_nlopt->nlopt == NULL)
  {
    fit_nlopt->nlopt      = nlopt_create (algo, fparam_len);
    fit_nlopt->nlopt_algo = algo;
  }

  if (fit_nlopt->nlopt_algo != algo)
  {
    fit_nlopt->nlopt_algo = algo;
    g_clear_pointer (&fit_nlopt->desc, g_free);
  }
}

/**
 * ncm_fit_nlopt_set_local_algo:
 * @fit_nlopt: a #NcmFitNLOpt.
 * @algo: a #NcmFitNloptAlgorithm.
 *
 * FIXME
 *
 */
void
ncm_fit_nlopt_set_local_algo (NcmFitNLOpt *fit_nlopt, NcmFitNloptAlgorithm algo)
{
  NcmFit *fit            = NCM_FIT (fit_nlopt);
  NcmFitState *fstate    = ncm_fit_peek_state (fit);
  const guint fparam_len = ncm_fit_state_get_fparam_len (fstate);

  if (fit_nlopt->local_nlopt_algo != algo)
    g_clear_pointer (&fit_nlopt->local_nlopt, nlopt_destroy);

  if (fit_nlopt->local_nlopt == NULL)
  {
    fit_nlopt->local_nlopt      = nlopt_create (algo, fparam_len);
    fit_nlopt->local_nlopt_algo = algo;
  }

  if (fit_nlopt->local_nlopt_algo != algo)
  {
    fit_nlopt->local_nlopt_algo = algo;
    g_clear_pointer (&fit_nlopt->desc, g_free);
  }
}

#endif /* NUMCOSMO_HAVE_NLOPT */

