/***************************************************************************
 *            nc_hicosmo_qspline.c
 *
 *  Wed February 15 11:31:28 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * @title: Spline Desceleration Parameter Model
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_qspline.h"
#include "math/ncm_spline_cubic_notaknot.h"

G_DEFINE_TYPE (NcHICosmoQSpline, nc_hicosmo_qspline, NC_TYPE_HICOSMO);

#define VECTOR     (model->params)
#define QSPLINE_H0 (ncm_vector_get (VECTOR, NC_HICOSMO_QSPLINE_H0))
#define OMEGA_T    (ncm_vector_get (VECTOR, NC_HICOSMO_QSPLINE_OMEGA_T))
#define AS_DRAG    (ncm_vector_get (VECTOR, NC_HICOSMO_QSPLINE_AS_DRAG))

static gdouble
_nc_hicosmo_qspline_dE2dz (gdouble E2, gdouble z, gpointer userdata)
{
  NcHICosmoQSpline *qs = NC_HICOSMO_QSPLINE (userdata);
  gdouble q;
  if (z > qs->z_f)
    q = ncm_spline_eval (qs->q_z, qs->z_f);
  else
    q = ncm_spline_eval (qs->q_z, z);
  return 2.0 * E2 * (q + 1.0) / (1.0 + z);
}

static void
_nc_hicosmo_qspline_prepare (NcHICosmoQSpline *qs)
{
  if (NCM_MODEL (qs)->pkey > qs->pkey)
  {
    ncm_spline_prepare (qs->q_z);
    ncm_ode_spline_prepare (qs->E2_z, qs);
    qs->pkey = NCM_MODEL (qs)->pkey;
  }
  else
    return;
}

static gdouble
_nc_hicosmo_qspline_E2 (NcmModel *model, gdouble z)
{
  NcHICosmoQSpline *qs = NC_HICOSMO_QSPLINE (model);
  _nc_hicosmo_qspline_prepare (qs);
  if (z > qs->z_f)
  {
    gdouble q = ncm_spline_eval (qs->q_z, qs->z_f);
    return ncm_spline_eval (qs->E2_z->s, qs->z_f) * pow ((1.0 + z) / (1.0 + qs->z_f), 2.0 * (q + 1.0));
  }
  else
    return ncm_spline_eval (qs->E2_z->s, z);
}

static gdouble
_nc_hicosmo_qspline_dE2_dz (NcmModel *model, gdouble z)
{
  NcHICosmoQSpline *qs = NC_HICOSMO_QSPLINE (model);
  _nc_hicosmo_qspline_prepare (qs);
  if (z > qs->z_f)
  {
    gdouble q = ncm_spline_eval (qs->q_z, qs->z_f);
    gdouble x = 1.0 + z;
    return ncm_spline_eval (qs->E2_z->s, qs->z_f) *
      pow (x / (1.0 + qs->z_f), 2.0 * (q + 1.0)) *
      2.0 * (q + 1.0) / x;
  }
  else
    return ncm_spline_eval_deriv (qs->E2_z->s, z);
}

static gdouble
_nc_hicosmo_qspline_d2E2_dz2 (NcmModel *model, gdouble z)
{
  NcHICosmoQSpline *qs = NC_HICOSMO_QSPLINE (model);
  _nc_hicosmo_qspline_prepare (qs);
  if (z > qs->z_f)
  {
    gdouble q = ncm_spline_eval (qs->q_z, qs->z_f);
    gdouble x = 1.0 + z;
    gdouble x2 = x * x;
    return ncm_spline_eval (qs->E2_z->s, qs->z_f) *
      pow (x / (1.0 + qs->z_f), 2.0 * (q + 1.0)) *
      2.0 * (q + 1.0) * (2.0 * q + 1.0) / x2;
  }
  else
    return ncm_spline_eval_deriv2 (qs->E2_z->s, z);
}

/****************************************************************************
 * Hubble constant
 ****************************************************************************/
static gdouble _nc_hicosmo_qspline_H0 (NcmModel *model) { return QSPLINE_H0; }
static gdouble _nc_hicosmo_qspline_Omega_t (NcmModel *model) { return OMEGA_T; }
static gdouble _nc_hicosmo_qspline_as_drag (NcmModel *model) { return AS_DRAG; }

/**
 * nc_hicosmo_qspline_new:
 * @s: FIXME
 * @np: FIXME
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

typedef struct _NcHICosmoSplineContPrior
{
  /*< private >*/
  gint knot;
  gdouble sigma;
} NcHICosmoSplineContPrior;

static void
continuity_prior_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcHICosmoQSpline *qspline = NC_HICOSMO_QSPLINE (ncm_mset_peek (mset, NC_HICOSMO_ID));
  NcHICosmoSplineContPrior *acp = (NcHICosmoSplineContPrior *) obj;
  const gdouble x_i = ncm_vector_get (qspline->q_z->xv, acp->knot);
  const gdouble x_ip1 = ncm_vector_get (qspline->q_z->xv, acp->knot + 1);
  const gdouble x_ip2 = ncm_vector_get (qspline->q_z->xv, acp->knot + 2);
  const gdouble mx1 = (x_ip1 + x_i) * 0.5;
  const gdouble mx2 = (x_ip2 + x_ip1) * 0.5;
  const gdouble d1 = ncm_spline_eval_deriv_nmax (qspline->q_z, mx1);
  const gdouble d2 = ncm_spline_eval_deriv_nmax (qspline->q_z, mx2);
  const gdouble mean_d12 = 1.0;//(d1 + d2) * 0.5;
  const gdouble mu = (d2 - d1) / mean_d12;
  //printf ("# [%u] meio [% 10.7g % 10.7g] derivs [% 10.7g % 10.7g] (% 10.7g)\n", acp->knot, mx1, mx2, d1, d2, mean_d12);
  f[0] = (mu / acp->sigma);
}

static void
_nc_hicosmo_spline_continuity_prior_free (gpointer obj)
{
  g_slice_free (NcHICosmoSplineContPrior, obj);
}

/**
 * nc_hicosmo_qspline_add_continuity_prior:
 * @qspline: FIXME
 * @lh: FIXME
 * @knot: FIXME
 * @sigma: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_qspline_add_continuity_prior (NcHICosmoQSpline *qspline, NcmLikelihood *lh, gint knot, gdouble sigma)
{
  NcHICosmoSplineContPrior *cp = g_slice_new (NcHICosmoSplineContPrior);
  NcmMSetFunc *func = ncm_mset_func_new (continuity_prior_f, 0, 1, cp, _nc_hicosmo_spline_continuity_prior_free);
  g_assert (knot < qspline->nknots - 1);
  g_assert (sigma > 0.0);
  cp->knot = knot;
  cp->sigma = sigma;
  ncm_likelihood_priors_add (lh, func);
  return;
}

/**
 * nc_hicosmo_qspline_add_continuity_priors:
 * @qspline: FIXME
 * @lh: FIXME
 * @sigma: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_qspline_add_continuity_priors (NcHICosmoQSpline *qspline, NcmLikelihood *lh, gdouble sigma)
{
  guint i;
  for (i = 0; i < qspline->nknots - 2; i++)
    nc_hicosmo_qspline_add_continuity_prior (qspline, lh, i, sigma);
  return;
}

enum {
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
  qspline->size = 0;
  qspline->z_f = 0.0;
  qspline->q_z = NULL;
  qspline->E2_z = NULL;
  qspline->pkey = NCM_MODEL (qspline)->pkey;
}

static void
_nc_hicosmo_qspline_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hicosmo_qspline_parent_class)->constructed (object);
  {
    NcHICosmoQSpline *qspline = NC_HICOSMO_QSPLINE (object);
    NcmModel *model = NCM_MODEL (qspline);
    NcmModelClass *model_class = NCM_MODEL_GET_CLASS (model);
    NcmVector *zv, *qv;
    guint i, qvi;

    qspline->nknots = ncm_model_vparam_len (model, NC_HICOSMO_QSPLINE_Q);
    qspline->size = model_class->sparam_len + qspline->nknots;

    qvi = ncm_model_vparam_index (model, NC_HICOSMO_QSPLINE_Q, 0);
    
    zv = ncm_vector_new (qspline->nknots);
    qv = ncm_vector_get_subvector (model->params, qvi, qspline->nknots);
    
    for (i = 0; i < qspline->nknots; i++)
    {
      gdouble zi = qspline->z_f / (qspline->nknots - 1.0) * i;
      gdouble xi = zi + 1.0;
      gdouble xi2 = xi * xi;
      gdouble xi3 = xi2 * xi;
      gdouble xi4 = xi2 * xi2;
      gdouble Ei2 = 1e-5 * xi4 + 0.25 * xi3 + (0.75 - 1e-5);
      gdouble qi = (1e-5 * xi4 + 0.25 * xi3 * 0.5 - (0.75 - 1e-5)) / Ei2;
      ncm_vector_set (zv, i, zi);
      ncm_vector_set (qv, i, qi);
    }

    {
      NcmSpline *s = ncm_spline_cubic_notaknot_new ();
      if (qspline->q_z == NULL)
        qspline->q_z = ncm_spline_new (s, zv, qv, FALSE);
      else
        ncm_spline_set (qspline->q_z, zv, qv, FALSE);
      
      qspline->E2_z = ncm_ode_spline_new (s,
                                          _nc_hicosmo_qspline_dE2dz, qspline,
                                          1.0, 0.0, qspline->z_f);
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
nc_hicosmo_qspline_dispose (GObject *object)
{
  NcHICosmoQSpline *qspline = NC_HICOSMO_QSPLINE (object);

  ncm_spline_clear (&qspline->q_z);
  ncm_ode_spline_clear (&qspline->E2_z);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qspline_parent_class)->dispose (object);
}

static void
nc_hicosmo_qspline_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qspline_parent_class)->finalize (object);
}

static void
nc_hicosmo_qspline_class_init (NcHICosmoQSplineClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcHICosmoClass* parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->set_property = &ncm_model_class_set_property;
  object_class->get_property = &ncm_model_class_get_property;
  object_class->constructed  = &_nc_hicosmo_qspline_constructed;
  object_class->dispose      = &nc_hicosmo_qspline_dispose;
  object_class->finalize     = &nc_hicosmo_qspline_finalize;

  model_class->set_property = &_nc_hicosmo_qspline_set_property;
  model_class->get_property = &_nc_hicosmo_qspline_get_property;

  ncm_model_class_add_params (model_class, NC_HICOSMO_QSPLINE_SPARAM_LEN, NC_HICOSMO_QSPLINE_VPARAM_LEN, PROP_SIZE);
  ncm_model_class_set_name_nick (model_class, "Q Spline", "qspline");

  /**
   * NcHICosmoQSpline:H0:
   *
   * FIXME
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
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QSPLINE_OMEGA_T, "Omega_t", "Omegat",
                              -5.0, 5.0, 1.0e-1,
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
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QSPLINE_AS_DRAG, "A_s drag", "asdrag",
                              0.0, 5.0, 1.0e-3,
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
  ncm_model_class_set_vparam (model_class, NC_HICOSMO_QSPLINE_Q, NC_HICOSMO_QSPLINE_DEFAULT_Q_LEN, "qparam", "qparam",
                              -50.0, 50.0, 1.0e-1, NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QSPLINE_DEFAULT_Q,
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
  nc_hicosmo_set_Omega_t_impl  (parent_class, &_nc_hicosmo_qspline_Omega_t);
  nc_hicosmo_set_as_drag_impl  (parent_class, &_nc_hicosmo_qspline_as_drag);
}
