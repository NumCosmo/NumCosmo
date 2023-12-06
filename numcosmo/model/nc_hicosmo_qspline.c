/***************************************************************************
 *            nc_hicosmo_qspline.c
 *
 *  Wed February 15 11:31:28 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
 * SECTION:nc_hicosmo_qspline
 * @title: NcHICosmoQSpline
 * @short_description: Kinetic model -- Spline deceleration function
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_qspline.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline_gsl.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_fit.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcHICosmoQSpline, nc_hicosmo_qspline, NC_TYPE_HICOSMO)

enum
{
  PROP_0,
  PROP_SPLINE,
  PROP_NKNOTS,
  PROP_Z_F,
  PROP_SIZE,
};

static void
nc_hicosmo_qspline_init (NcHICosmoQSpline *qspline)
{
  qspline->nknots = 0;
  qspline->size   = 0;
  qspline->z_f    = 0.0;
  qspline->q_z    = NULL;
  qspline->E2_z   = NULL;
}

static gdouble _nc_hicosmo_qspline_dE2dz (gdouble E2, gdouble z, gpointer userdata);

static void
_nc_hicosmo_qspline_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hicosmo_qspline_parent_class)->constructed (object);
  {
    NcHICosmoQSpline *qspline  = NC_HICOSMO_QSPLINE (object);
    NcmModel *model            = NCM_MODEL (qspline);
    NcmVector *orig_vec        = ncm_model_orig_params_peek_vector (model);
    NcmModelClass *model_class = NCM_MODEL_GET_CLASS (model);
    NcmVector *zv, *qv;
    guint i, qvi;
    guint qz_size = ncm_model_vparam_len (model, NC_HICOSMO_QSPLINE_Q);

    qspline->nknots = qz_size;
    qspline->size   = model_class->sparam_len + qspline->nknots;

    qvi = ncm_model_vparam_index (model, NC_HICOSMO_QSPLINE_Q, 0);

    zv = ncm_vector_new (qz_size);
    qv = ncm_vector_get_subvector (orig_vec, qvi, qz_size);

    for (i = 0; i < qz_size; i++)
    {
      gdouble zi  = qspline->z_f / (qz_size - 1.0) * i;
      gdouble xi  = zi + 1.0;
      gdouble xi2 = xi * xi;
      gdouble xi3 = xi2 * xi;
      gdouble xi4 = xi2 * xi2;
      gdouble Ei2 = 1e-5 * xi4 + 0.25 * xi3 + (0.75 - 1e-5);
      gdouble qi  = (1e-5 * xi4 + 0.25 * xi3 * 0.5 - (0.75 - 1e-5)) / Ei2;

      ncm_vector_set (zv, i, zi);
      ncm_vector_set (qv, i, qi);
    }

    {
      NcmSpline *s = ncm_spline_cubic_notaknot_new ();

      if (qspline->q_z == NULL)
        qspline->q_z = ncm_spline_new (s, zv, qv, FALSE);
      else
        ncm_spline_set (qspline->q_z, zv, qv, FALSE);

      qspline->E2_z = ncm_ode_spline_new_full (s,
                                               _nc_hicosmo_qspline_dE2dz,
                                               1.0, 0.0, qspline->z_f);
      ncm_ode_spline_set_reltol (qspline->E2_z, 1.0e-9);
      ncm_ode_spline_set_abstol (qspline->E2_z, 0.0);

      ncm_spline_free (s);
      ncm_vector_free (zv);
      ncm_vector_free (qv);
    }

    return;
  }
}

static void
_nc_hicosmo_qspline_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHICosmoQSpline *qspline = NC_HICOSMO_QSPLINE (object);

  g_return_if_fail (NC_IS_HICOSMO_QSPLINE (object));

  switch (prop_id)
  {
    case PROP_SPLINE:
      g_value_set_object (value, qspline->q_z);
      break;
    case PROP_NKNOTS:
      g_value_set_uint (value, qspline->nknots);
      break;
    case PROP_Z_F:
      g_value_set_double (value, qspline->z_f);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hicosmo_qspline_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHICosmoQSpline *qspline = NC_HICOSMO_QSPLINE (object);

  g_return_if_fail (NC_IS_HICOSMO_QSPLINE (object));

  switch (prop_id)
  {
    case PROP_SPLINE:
      qspline->q_z = g_value_dup_object (value);
      break;
    case PROP_NKNOTS:
      qspline->nknots = g_value_get_uint (value);
      break;
    case PROP_Z_F:
      qspline->z_f = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hicosmo_qspline_dispose (GObject *object)
{
  NcHICosmoQSpline *qspline = NC_HICOSMO_QSPLINE (object);

  ncm_spline_clear (&qspline->q_z);
  ncm_ode_spline_clear (&qspline->E2_z);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qspline_parent_class)->dispose (object);
}

static void
_nc_hicosmo_qspline_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qspline_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_qspline_H0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_qspline_E2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_qspline_dE2_dz (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_qspline_d2E2_dz2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_qspline_Omega_t0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_qspline_as_drag (NcHICosmo *cosmo);

static void
nc_hicosmo_qspline_class_init (NcHICosmoQSplineClass *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass *parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_hicosmo_qspline_set_property;
  model_class->get_property = &_nc_hicosmo_qspline_get_property;
  object_class->constructed = &_nc_hicosmo_qspline_constructed;
  object_class->dispose     = &_nc_hicosmo_qspline_dispose;
  object_class->finalize    = &_nc_hicosmo_qspline_finalize;

  ncm_model_class_set_name_nick (model_class, "Q Spline", "qspline");
  ncm_model_class_add_params (model_class, NC_HICOSMO_QSPLINE_SPARAM_LEN, NC_HICOSMO_QSPLINE_VPARAM_LEN, PROP_SIZE);

  /**
   * NcHICosmoQSpline:H0:
   *
   * Hubble parameter today (z = 0).
   */
  /**
   * NcHICosmoQSpline:H0-fit:
   *
   * FIXME
   */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QSPLINE_H0, "H_0", "H0",
                              10.0, 500.0, 1.0, NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QSPLINE_DEFAULT_H0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcHICosmoQSpline:Omegat:
   *
   * FIXME
   */
  /**
   * NcHICosmoQSpline:Omegat-fit:
   *
   * FIXME
   */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QSPLINE_OMEGA_T, "\\Omega_{t0}", "Omegat",
                              0.05, 2.0, 1.0e-1,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QSPLINE_DEFAULT_OMEGA_T,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcHICosmoQSpline:asdrag:
   *
   * FIXME
   */
  /**
   * NcHICosmoQSpline:asdrag-fit:
   *
   * FIXME
   */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QSPLINE_AS_DRAG, "A_s", "asdrag",
                              1.0e-4, 1.0, 1.0e-3,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QSPLINE_DEFAULT_AS_DRAG,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcHICosmoQSpline:qparam:
   *
   * FIXME
   */
  /**
   * NcHICosmoQSpline:qparam-fit:
   *
   * FIXME
   */
  /**
   * NcHICosmoQSpline:qparam-length:
   *
   * FIXME
   */
  ncm_model_class_set_vparam (model_class, NC_HICOSMO_QSPLINE_Q, NC_HICOSMO_QSPLINE_DEFAULT_Q_LEN, "q", "qparam",
                              -50.0, 50.0, 2.0e-1, NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QSPLINE_DEFAULT_Q,
                              NCM_PARAM_TYPE_FREE);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  g_object_class_install_property (object_class,
                                   PROP_SPLINE,
                                   g_param_spec_object ("spline",
                                                        NULL,
                                                        "Spline object",
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Z_F,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "final redshift",
                                                        0.0, 100.0, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  nc_hicosmo_set_H0_impl       (parent_class, &_nc_hicosmo_qspline_H0);
  nc_hicosmo_set_E2_impl       (parent_class, &_nc_hicosmo_qspline_E2);
  nc_hicosmo_set_dE2_dz_impl   (parent_class, &_nc_hicosmo_qspline_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl (parent_class, &_nc_hicosmo_qspline_d2E2_dz2);
  nc_hicosmo_set_Omega_t0_impl (parent_class, &_nc_hicosmo_qspline_Omega_t0);
  nc_hicosmo_set_as_drag_impl  (parent_class, &_nc_hicosmo_qspline_as_drag);
}

#define VECTOR (ncm_model_orig_params_peek_vector (NCM_MODEL (cosmo)))
#define QSPLINE_H0 (ncm_vector_get (VECTOR, NC_HICOSMO_QSPLINE_H0))
#define OMEGA_T    (ncm_vector_get (VECTOR, NC_HICOSMO_QSPLINE_OMEGA_T))
#define AS_DRAG    (ncm_vector_get (VECTOR, NC_HICOSMO_QSPLINE_AS_DRAG))

static void
_nc_hicosmo_qspline_prepare (NcHICosmoQSpline *qs)
{
  if (!ncm_model_state_is_update (NCM_MODEL (qs)))
  {
    /*static guint count = 0;*/
    ncm_spline_prepare (qs->q_z);
    ncm_model_state_set_update (NCM_MODEL (qs));

    /*printf ("Preparing! ==>(%20u)<==\n", count++);*/
    ncm_ode_spline_prepare (qs->E2_z, qs);
  }
  else
  {
    return;
  }
}

static gdouble
_nc_hicosmo_qspline_E2 (NcHICosmo *cosmo, gdouble z)
{
  NcHICosmoQSpline *qs = NC_HICOSMO_QSPLINE (cosmo);

  _nc_hicosmo_qspline_prepare (qs);

  if (z > qs->z_f)
  {
    gdouble q = ncm_spline_eval (qs->q_z, qs->z_f);

    return ncm_spline_eval (ncm_ode_spline_peek_spline (qs->E2_z), qs->z_f) * pow ((1.0 + z) / (1.0 + qs->z_f), 2.0 * (q + 1.0));
  }
  else
  {
    return ncm_spline_eval (ncm_ode_spline_peek_spline (qs->E2_z), z);
  }
}

static gdouble
_nc_hicosmo_qspline_dE2dz (gdouble E2, gdouble z, gpointer userdata)
{
  NcHICosmoQSpline *qs = NC_HICOSMO_QSPLINE (userdata);
  gdouble q;

  _nc_hicosmo_qspline_prepare (qs);

  if (z > qs->z_f)
    q = ncm_spline_eval (qs->q_z, qs->z_f);
  else
    q = ncm_spline_eval (qs->q_z, z);

  return 2.0 * E2 * (q + 1.0) / (1.0 + z);
}

static gdouble
_nc_hicosmo_qspline_dE2_dz (NcHICosmo *cosmo, gdouble z)
{
  NcHICosmoQSpline *qs = NC_HICOSMO_QSPLINE (cosmo);

  _nc_hicosmo_qspline_prepare (qs);

  if (z > qs->z_f)
  {
    const gdouble q = ncm_spline_eval (qs->q_z, qs->z_f);
    const gdouble x = 1.0 + z;

    return ncm_spline_eval (ncm_ode_spline_peek_spline (qs->E2_z), qs->z_f) *
           pow (x / (1.0 + qs->z_f), 2.0 * (q + 1.0)) *
           2.0 * (q + 1.0) / x;
  }
  else
  {
    return ncm_spline_eval_deriv (ncm_ode_spline_peek_spline (qs->E2_z), z);
  }
}

static gdouble
_nc_hicosmo_qspline_d2E2_dz2 (NcHICosmo *cosmo, gdouble z)
{
  NcHICosmoQSpline *qs = NC_HICOSMO_QSPLINE (cosmo);

  _nc_hicosmo_qspline_prepare (qs);

  if (z > qs->z_f)
  {
    const gdouble q  = ncm_spline_eval (qs->q_z, qs->z_f);
    const gdouble x  = 1.0 + z;
    const gdouble x2 = x * x;

    return ncm_spline_eval (ncm_ode_spline_peek_spline (qs->E2_z), qs->z_f) *
           pow (x / (1.0 + qs->z_f), 2.0 * (q + 1.0)) *
           2.0 * (q + 1.0) * (2.0 * q + 1.0) / x2;
  }
  else
  {
    return ncm_spline_eval_deriv2 (ncm_ode_spline_peek_spline (qs->E2_z), z);
  }
}

static gdouble
_nc_hicosmo_qspline_H0 (NcHICosmo *cosmo)
{
  return QSPLINE_H0;
}

static gdouble
_nc_hicosmo_qspline_Omega_t0 (NcHICosmo *cosmo)
{
  return OMEGA_T;
}

static gdouble
_nc_hicosmo_qspline_as_drag (NcHICosmo *cosmo)
{
  return AS_DRAG;
}

/**
 * nc_hicosmo_qspline_new:
 * @s: a #NcmSpline
 * @np: number of knots
 * @z_f: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcHICosmoQSpline *
nc_hicosmo_qspline_new (NcmSpline *s, gsize np, gdouble z_f)
{
  NcHICosmoQSpline *qspline = g_object_new (NC_TYPE_HICOSMO_QSPLINE,
                                            "spline", s,
                                            "zf", z_f,
                                            "qparam-length", np,
                                            NULL);

  return qspline;
}

/**
 * nc_hicosmo_qspline_add_continuity_priors:
 * @qspline: a #NcHICosmoQSpline
 * @lh: a #NcmLikelihood
 * @sigma: FIXME
 * @abstol: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcHICosmoQSplineContPrior *
nc_hicosmo_qspline_add_continuity_priors (NcHICosmoQSpline *qspline, NcmLikelihood *lh, gdouble sigma, gdouble abstol)
{
  g_assert (sigma > 0 && qspline->nknots > 2);
  {
    NcHICosmoQSplineContPrior *qspline_cp = nc_hicosmo_qspline_cont_prior_new (qspline->nknots - 2);
    NcPriorQSplineCont *prior_qz_cont     = nc_prior_qspline_cont_new ();

    nc_hicosmo_qspline_cont_prior_set_all_lnsigma (qspline_cp, log (sigma));
    nc_hicosmo_qspline_cont_prior_set_abstol (qspline_cp, abstol);

    ncm_likelihood_priors_add (lh, NCM_PRIOR (prior_qz_cont));

    return qspline_cp;
  }
}

/****************************** Continuity Prior Model ******************************/

enum
{
  PROP_CP_0,
  PROP_CP_SIZE,
};

G_DEFINE_TYPE (NcHICosmoQSplineContPrior, nc_hicosmo_qspline_cont_prior, NCM_TYPE_MODEL)

static void
nc_hicosmo_qspline_cont_prior_init (NcHICosmoQSplineContPrior *qspline_cp)
{
  NCM_UNUSED (qspline_cp);
}

static void
_nc_hicosmo_qspline_cont_prior_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hicosmo_qspline_cont_prior_parent_class)->constructed (object);
  {
  }
}

static void
_nc_hicosmo_qspline_cont_prior_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qspline_cont_prior_parent_class)->dispose (object);
}

static void
_nc_hicosmo_qspline_cont_prior_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qspline_cont_prior_parent_class)->finalize (object);
}

static void
_nc_hicosmo_qspline_cont_prior_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
/*  NcHICosmoQSplineContPrior *qspline_cp = NC_HICOSMO_QSPLINE_CONT_PRIOR (object); */
  g_return_if_fail (NC_IS_HICOSMO_QSPLINE_CONT_PRIOR (object));
  NCM_UNUSED (value);

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hicosmo_qspline_cont_prior_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
/*  NcHICosmoQSplineContPrior *qspline_cp = NC_HICOSMO_QSPLINE_CONT_PRIOR (object); */
  g_return_if_fail (NC_IS_HICOSMO_QSPLINE_CONT_PRIOR (object));

  NCM_UNUSED (value);

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static gboolean _nc_hicosmo_qspline_cont_prior_valid (NcmModel *model);

NCM_MSET_MODEL_REGISTER_ID (nc_hicosmo_qspline_cont_prior, NC_TYPE_HICOSMO_QSPLINE_CONT_PRIOR);

static void
nc_hicosmo_qspline_cont_prior_class_init (NcHICosmoQSplineContPriorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->constructed = &_nc_hicosmo_qspline_cont_prior_constructed;
  object_class->dispose     = &_nc_hicosmo_qspline_cont_prior_dispose;
  object_class->finalize    = &_nc_hicosmo_qspline_cont_prior_finalize;

  model_class->set_property = &_nc_hicosmo_qspline_cont_prior_set_property;
  model_class->get_property = &_nc_hicosmo_qspline_cont_prior_get_property;

  ncm_model_class_set_name_nick (model_class, "Q Spline Cont Prior", "qspline_cp");
  ncm_model_class_add_params (model_class, 1, 1, PROP_CP_SIZE);

  /**
   * NcHICosmoQSplineContPrior:abstol:
   *
   * FIXME
   */
  /**
   * NcHICosmoQSplineContPrior:abstol-fit:
   *
   * FIXME
   */
  /**
   * NcHICosmoQSplineContPrior:abstol-length:
   *
   * FIXME
   */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QSPLINE_CONT_PRIOR_ABSTOL, "abstol", "abstol",
                              0.0, G_MAXDOUBLE, 1.0e-1, 0.0, 1.0e-1,
                              NCM_PARAM_TYPE_FREE);

  /**
   * NcHICosmoQSplineContPrior:lnsigma:
   *
   * FIXME
   */
  /**
   * NcHICosmoQSplineContPrior:lnsigma-fit:
   *
   * FIXME
   */
  /**
   * NcHICosmoQSplineContPrior:lnsigma-length:
   *
   * FIXME
   */
  ncm_model_class_set_vparam (model_class, NC_HICOSMO_QSPLINE_CONT_PRIOR_LNSIGMA, 3, "lnsigma", "lnsigma",
                              log (1e-200), log (1.0e200), fabs (log (1.0e-1)), 0.0, log (1.0e-1),
                              NCM_PARAM_TYPE_FREE);

  ncm_model_class_check_params_info (model_class);

  ncm_mset_model_register_id (model_class,
                              "NcHICosmoQSplineContPrior",
                              "NcHICosmoQSplineContPrior",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);

  model_class->valid = &_nc_hicosmo_qspline_cont_prior_valid;
}

static gboolean
_nc_hicosmo_qspline_cont_prior_valid (NcmModel *model)
{
  if (!NCM_MODEL_CLASS (nc_hicosmo_qspline_cont_prior_parent_class)->valid (model))
    return FALSE;

  /* Chain up : start */

  return TRUE;
}

/**
 * nc_hicosmo_qspline_cont_prior_new:
 * @npriors: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcHICosmoQSplineContPrior *
nc_hicosmo_qspline_cont_prior_new (guint npriors)
{
  NcHICosmoQSplineContPrior *qspline_cp = g_object_new (NC_TYPE_HICOSMO_QSPLINE_CONT_PRIOR,
                                                        "lnsigma-length", npriors,
                                                        NULL);

  return qspline_cp;
}

/**
 * nc_hicosmo_qspline_cont_prior_ref:
 * @qspline_cp: a #NcHICosmoQSplineContPrior
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcHICosmoQSplineContPrior *
nc_hicosmo_qspline_cont_prior_ref (NcHICosmoQSplineContPrior *qspline_cp)
{
  return g_object_ref (qspline_cp);
}

/**
 * nc_hicosmo_qspline_cont_prior_free:
 * @qspline_cp: a #NcHICosmoQSplineContPrior
 *
 * FIXME
 *
 */
void
nc_hicosmo_qspline_cont_prior_free (NcHICosmoQSplineContPrior *qspline_cp)
{
  g_object_unref (qspline_cp);
}

/**
 * nc_hicosmo_qspline_cont_prior_set_lnsigma:
 * @qspline_cp: a #NcHICosmoQSplineContPrior
 * @i: FIXME
 * @ln_sigma: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_qspline_cont_prior_set_lnsigma (NcHICosmoQSplineContPrior *qspline_cp, guint i, gdouble ln_sigma)
{
  NcmModel *model = NCM_MODEL (qspline_cp);
  guint npriors   = ncm_model_vparam_len (model, NC_HICOSMO_QSPLINE_CONT_PRIOR_LNSIGMA);

  g_assert_cmpint (i, <, npriors);

  ncm_model_orig_vparam_set (model, NC_HICOSMO_QSPLINE_CONT_PRIOR_LNSIGMA, i, ln_sigma);
}

/**
 * nc_hicosmo_qspline_cont_prior_set_all_lnsigma:
 * @qspline_cp: a #NcHICosmoQSplineContPrior
 * @ln_sigma: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_qspline_cont_prior_set_all_lnsigma (NcHICosmoQSplineContPrior *qspline_cp, gdouble ln_sigma)
{
  NcmModel *model = NCM_MODEL (qspline_cp);
  guint npriors   = ncm_model_vparam_len (model, NC_HICOSMO_QSPLINE_CONT_PRIOR_LNSIGMA);
  guint i;

  for (i = 0; i < npriors; i++)
    ncm_model_orig_vparam_set (model, NC_HICOSMO_QSPLINE_CONT_PRIOR_LNSIGMA, i, ln_sigma);
}

/**
 * nc_hicosmo_qspline_cont_prior_get_lnsigma:
 * @qspline_cp: a #NcHICosmoQSplineContPrior
 * @i: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hicosmo_qspline_cont_prior_get_lnsigma (NcHICosmoQSplineContPrior *qspline_cp, guint i)
{
  NcmModel *model = NCM_MODEL (qspline_cp);
  guint npriors   = ncm_model_vparam_len (model, NC_HICOSMO_QSPLINE_CONT_PRIOR_LNSIGMA);

  g_assert_cmpint (i, <, npriors);

  return ncm_model_orig_vparam_get (model, NC_HICOSMO_QSPLINE_CONT_PRIOR_LNSIGMA, i);
}

/**
 * nc_hicosmo_qspline_cont_prior_set_abstol:
 * @qspline_cp: a #NcHICosmoQSplineContPrior
 * @abstol: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_qspline_cont_prior_set_abstol (NcHICosmoQSplineContPrior *qspline_cp, gdouble abstol)
{
  NcmModel *model = NCM_MODEL (qspline_cp);

  g_assert_cmpfloat (abstol, >, 0.0);
  ncm_model_orig_param_set (model, NC_HICOSMO_QSPLINE_CONT_PRIOR_ABSTOL, abstol);
}

/**
 * nc_hicosmo_qspline_cont_prior_get_abstol:
 * @qspline_cp: a #NcHICosmoQSplineContPrior
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hicosmo_qspline_cont_prior_get_abstol (NcHICosmoQSplineContPrior *qspline_cp)
{
  NcmModel *model = NCM_MODEL (qspline_cp);

  return ncm_model_orig_param_get (model, NC_HICOSMO_QSPLINE_CONT_PRIOR_ABSTOL);
}

/****************************** Continuity Prior ******************************/

G_DEFINE_TYPE (NcPriorQSplineCont, nc_prior_qspline_cont, NCM_TYPE_PRIOR)

static void
nc_prior_qspline_cont_init (NcPriorQSplineCont *nc_prior_qspline_cont)
{
}

static void _nc_prior_qspline_cont_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res);

static void
nc_prior_qspline_cont_class_init (NcPriorQSplineContClass *klass)
{
  NcmPriorClass *prior_class   = NCM_PRIOR_CLASS (klass);
  NcmMSetFuncClass *func_class = NCM_MSET_FUNC_CLASS (klass);

  prior_class->is_m2lnL = TRUE;
  func_class->eval      = &_nc_prior_qspline_cont_eval;
}

static void
_nc_prior_qspline_cont_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  NcHICosmoQSpline *qspline             = NC_HICOSMO_QSPLINE (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHICosmoQSplineContPrior *qspline_cp = NC_HICOSMO_QSPLINE_CONT_PRIOR (ncm_mset_peek (mset, nc_hicosmo_qspline_cont_prior_id ()));
  guint knot;

  NCM_UNUSED (x);

  res[0] = 0.0;

  for (knot = 0; knot < qspline->nknots - 2; knot++)
  {
    const gdouble sigma  = exp (nc_hicosmo_qspline_cont_prior_get_lnsigma (qspline_cp, knot));
    const gdouble abstol = nc_hicosmo_qspline_cont_prior_get_abstol (qspline_cp);
    const gdouble zi     = ncm_vector_get (qspline->q_z->xv, knot + 0);
    const gdouble zip1   = ncm_vector_get (qspline->q_z->xv, knot + 1);
    const gdouble zip2   = ncm_vector_get (qspline->q_z->xv, knot + 2);
    const gdouble q1     = ncm_spline_eval (qspline->q_z, zi);
    const gdouble q2     = ncm_spline_eval (qspline->q_z, zip1);
    const gdouble q3     = ncm_spline_eval (qspline->q_z, zip2);
    const gdouble qmean  = q1 + (q3 - q1) * (zip1 - zi) / (zip2 - zi);

    res[0] += gsl_pow_2 ((q2 - qmean) / (abstol + (fabs (q2) + fabs (qmean)) * 0.5 * sigma));
    /*res[0] += gsl_pow_2 ((q2 - qmean) / (sigma + fabs (qmean) * abstol));*/
  }
}

/**
 * nc_prior_qspline_cont_new:
 *
 * FIXME
 *
 * Returns: (transfer full): the newly created #NcPriorQSplineCont.
 */
NcPriorQSplineCont *
nc_prior_qspline_cont_new (void)
{
  NcPriorQSplineCont *prior_qz_cp = g_object_new (NC_TYPE_PRIOR_QSPLINE_CONT,
                                                  NULL);

  return prior_qz_cp;
}

