/***************************************************************************
 *            nc_hicosmo_qpw.c
 *
 *  Mon Jun 18 18:16:35 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
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
 * SECTION:nc_hicosmo_qpw
 * @title: Linear Spline Desceleration Parameter Model
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_qpw.h"
#include "math/ncm_reparam_linear.h"

G_DEFINE_TYPE (NcHICosmoQPW, nc_hicosmo_qpw, NC_TYPE_HICOSMO);

#define VECTOR    (model->params)
#define QPW_H0    (ncm_vector_get (VECTOR, NC_HICOSMO_QPW_H0))
#define OMEGA_T   (ncm_vector_get (VECTOR, NC_HICOSMO_QPW_OMEGA_T))
#define QPW_Q0    (ncm_vector_get (VECTOR, NC_HICOSMO_QPW_Q0))
#define QPW_QP(i) (ncm_vector_get (VECTOR, NC_HICOSMO_QPW_QP + i))

#define LOG(x) log(x)
#define EXP(x) exp(x)

static gdouble
_nc_hicosmo_qpw_E2 (NcmModel *model, gdouble z)
{
  NcHICosmoQPW *qpw = NC_HICOSMO_QPW (model);

  gint i, n = nc_hicosmo_qpw_index (qpw, z);
  gdouble result = 0.0;
  gdouble q_i = QPW_Q0;
  gdouble x_n = (1.0 + n * qpw->piece);
  gdouble x = 1.0 + z;
  for (i = 0; i < n; i++)
  {
    gdouble qp_i = QPW_QP(i);
    gdouble x_i = 1.0 + i * qpw->piece;
    gdouble x_i_1 = x_i + qpw->piece;
    result += (1.0 + q_i - qp_i * x_i) * LOG (x_i_1 / x_i) + qp_i * qpw->piece;
    q_i += qpw->piece * qp_i;
  }
  result += (1.0 + q_i - QPW_QP(n) * x_n) * LOG (x / x_n) + QPW_QP(n) * (x - x_n);
  return EXP ( 2.0 * result );
}

static gdouble
_nc_hicosmo_qpw_dE2_dz (NcmModel *model, gdouble z)
{
  NcHICosmoQPW *qpw = NC_HICOSMO_QPW (model);
  gint i, n = nc_hicosmo_qpw_index (qpw, z);
  gdouble q_i = QPW_Q0;
  gdouble z_n = n * qpw->piece;
  for (i = 0; i < n; i++)
    q_i += qpw->piece * QPW_QP(i);

  return 2.0 * (1.0 + q_i + QPW_QP(i) * (z - z_n)) / (1.0 + z) * _nc_hicosmo_qpw_E2 (model, z);
}

static gdouble
_nc_hicosmo_qpw_d2E2_dz2 (NcmModel *model, gdouble z)
{
  NcHICosmoQPW *qpw = NC_HICOSMO_QPW (model);
  gint i, n = nc_hicosmo_qpw_index (qpw, z);
  gdouble q_i = QPW_Q0;
  gdouble z_n = n * qpw->piece;
  gdouble inte, dinte;
  for (i = 0; i < n; i++)
    q_i += qpw->piece * QPW_QP(i);
  inte = (1.0 + q_i + QPW_QP(i) * (z - z_n)) / (1.0 + z);
  dinte = QPW_QP(i) / (1.0 + z) - inte / (1.0 + z);

  return  (4.0 * gsl_pow_2 (inte) + 2.0 * dinte) * _nc_hicosmo_qpw_E2 (model, z);
}

/**
 * nc_hicosmo_qpw_index:
 * @qpw: FIXME
 * @z: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint
nc_hicosmo_qpw_index (NcHICosmoQPW *qpw, gdouble z)
{
  gint wpiece = floor(z / qpw->piece);
  if (fabs(floor(z / qpw->piece) - (z / qpw->piece)) < NCM_ZERO_LIMIT)
    wpiece = wpiece == 0 ? 0 : wpiece - 1;
  return GSL_MIN(wpiece, qpw->npieces - 1);
}

/****************************************************************************
 * Hubble constant
 ****************************************************************************/
static gdouble _nc_hicosmo_qpw_H0 (NcmModel *model) { return QPW_H0; }
static gdouble _nc_hicosmo_qpw_Omega_t (NcmModel *model) { return OMEGA_T; }

static void
AsymptoticCDM_prior_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcmModel *model = ncm_mset_peek (mset, nc_hicosmo_id ());
  NcHICosmoQPWAsymCDMPrior *acp = (NcHICosmoQPWAsymCDMPrior *) obj;
  f[0] = (nc_hicosmo_q (NC_HICOSMO (model), acp->z) - acp->q) / acp->sigma;
}

static void
continuity_prior_f (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcmModel *model = ncm_mset_peek (mset, nc_hicosmo_id ());
  NcHICosmoQPWContPrior *acp = (NcHICosmoQPWContPrior *) obj;
  f[0] = ((QPW_QP (acp->knot) - QPW_QP (acp->knot + 1)) / acp->sigma);
}

static void
_nc_hicosmo_qpw_continuity_prior_free (gpointer obj)
{
  g_slice_free (NcHICosmoQPWContPrior, obj);
}

/**
 * nc_hicosmo_qpw_add_continuity_prior:
 * @qpw: FIXME
 * @lh: FIXME
 * @knot: FIXME
 * @sigma: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_qpw_add_continuity_prior (NcHICosmoQPW *qpw, NcmLikelihood *lh, gint knot, gdouble sigma)
{
  NcHICosmoQPWContPrior *cp = g_slice_new (NcHICosmoQPWContPrior);
  NcmMSetFunc *func = ncm_mset_func_new (continuity_prior_f, 0, 1, cp, _nc_hicosmo_qpw_continuity_prior_free);
  g_assert (knot < ncm_model_len (NCM_MODEL (qpw)) - 4);
  g_assert (sigma > 0.0);
  cp->knot = knot;
  cp->sigma = sigma;
  ncm_likelihood_priors_add (lh, func, FALSE);
  return;
}

/**
 * nc_hicosmo_qpw_add_continuity_priors:
 * @qpw: FIXME
 * @lh: FIXME
 * @sigma: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_qpw_add_continuity_priors (NcHICosmoQPW *qpw, NcmLikelihood *lh, gdouble sigma)
{
  guint i;
  for (i = 0; i < ncm_model_len (NCM_MODEL (qpw)) - 4; i++)
    nc_hicosmo_qpw_add_continuity_prior (qpw, lh, i, sigma);
  return;
}

static void
_nc_hicosmo_qpw_asymptotic_cdm_free (gpointer obj)
{
  g_slice_free (NcHICosmoQPWAsymCDMPrior, obj);
}

/**
 * nc_hicosmo_qpw_add_asymptotic_cdm_prior:
 * @qpw: FIXME
 * @lh: FIXME
 * @z: FIXME
 * @q: FIXME
 * @sigma: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_qpw_add_asymptotic_cdm_prior (NcHICosmoQPW *qpw, NcmLikelihood *lh, gdouble z, gdouble q, gdouble sigma)
{
  NcHICosmoQPWAsymCDMPrior *cp = g_slice_new (NcHICosmoQPWAsymCDMPrior);
  NcmMSetFunc *func = ncm_mset_func_new (AsymptoticCDM_prior_f, 0, 1, cp, _nc_hicosmo_qpw_asymptotic_cdm_free);
  cp->z = z;
  cp->q = q;
  cp->sigma = sigma;
  ncm_likelihood_priors_add (lh, func, FALSE);
  return;
}

/**
 * nc_hicosmo_qpw_change_params:
 * @qpw: FIXME
 * @z: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_qpw_change_params (NcHICosmoQPW *qpw, gdouble z)
{
  guint wpiece = nc_hicosmo_qpw_index (qpw, z);
  NcmMatrix *full_T = ncm_matrix_new (ncm_model_len (NCM_MODEL (qpw)), ncm_model_len (NCM_MODEL (qpw)));
  NcmVector *full_v = ncm_vector_new (ncm_model_len (NCM_MODEL (qpw)));
  NcmMatrix *T = ncm_matrix_get_submatrix (full_T, 2, 2, ncm_model_len (NCM_MODEL (qpw)) - 2, ncm_model_len (NCM_MODEL (qpw)) - 2);
  NcmVector *v = ncm_vector_get_subvector (full_v, 2, ncm_model_len (NCM_MODEL (qpw)) - 2);
  NcmReparamLinear *relin;

  ncm_matrix_set_identity (full_T);
  ncm_vector_set_zero (full_v);

  if (wpiece == 0)
  {
    gdouble x = 1.0 + z;
    gdouble ln_x = LOG (x);
    gdouble denom = x * ln_x - z;

    ncm_matrix_set (T, 0, 0, z / denom);
    ncm_matrix_set (T, 0, 1, (ln_x - z) / denom);
    ncm_vector_set (v, 0, -z * ln_x / denom);

    ncm_matrix_set (T, 1, 0, -1.0 / denom);
    ncm_matrix_set (T, 1, 1, ln_x / denom);
    ncm_vector_set (v, 1, ln_x / denom);
  }
  else
  {
    guint i;
    gdouble x = 1.0 + z;
    gdouble ln_x = LOG (x);
    gdouble x_ln_x = x * ln_x;
    gdouble x_n = 1.0 + wpiece * qpw->piece;
    gdouble x_n_ln_x_n = x_n * LOG (x_n);
    gdouble denom = x - x_n - x_ln_x + x_n_ln_x_n;

    ncm_matrix_set (T, 0, 0, (x_n - x) / denom);

    for (i = 0; i < wpiece; i++)
    {
      gdouble x_i = 1.0 + qpw->piece * i;
      gdouble x_i_1 = x_i + qpw->piece;
      ncm_matrix_set (T, 0, i + 1,
                      -qpw->piece +
                      (x - x_n) *
                      (qpw->piece + x_i * LOG (x_i) - x_i_1 * LOG (x_i_1)) / denom
                      );
    }
    ncm_matrix_set (T, 0, wpiece + 1, 1.0 + (x - x_n) * ln_x / denom);
    ncm_vector_set (v, 0, (x - x_n) * ln_x / denom);

    ncm_matrix_set (T, wpiece + 1, 0, 1.0 / denom);
    for (i = 0; i < wpiece; i++)
    {
      gdouble x_i = 1.0 + qpw->piece * i;
      gdouble x_i_1 = x_i + qpw->piece;
      ncm_matrix_set (T, wpiece + 1, i + 1, -(qpw->piece + x_i * LOG (x_i) - x_i_1 * LOG (x_i_1)) / denom );
    }
    ncm_matrix_set (T, wpiece + 1, wpiece + 1, -ln_x / denom);
    ncm_vector_set (v, wpiece + 1, -ln_x/denom);
  }

  if (FALSE)
  {
    gint i,j;
    printf ("wpiece: %d, npieces: %d\n", wpiece, qpw->npieces);
    for (i = 0; i < ncm_model_len (NCM_MODEL (qpw)); i++)
    {
      for (j = 0; j < ncm_model_len (NCM_MODEL (qpw)); j++)
        printf ("% -12.2f ", ncm_matrix_get (full_T, i, j));
      printf ("\n");
    }
    printf("########################\n");
    for (i = 0; i < ncm_model_len (NCM_MODEL (qpw)); i++)
      printf ("% -12.2f ", ncm_vector_get (full_v, i));
    printf("\n########################\n");
  }

  relin = ncm_reparam_linear_new (ncm_model_len (NCM_MODEL (qpw)), full_T, full_v);

  ncm_model_set_reparam (NCM_MODEL (qpw), NCM_REPARAM (relin));

  ncm_vector_free (v);
  ncm_matrix_free (T);
  ncm_vector_free (full_v);
  ncm_matrix_free (full_T);

  return;
}

/**
 * nc_hicosmo_qpw_change_params_qpp:
 * @qpw: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_qpw_change_params_qpp (NcHICosmoQPW *qpw)
{
  NcmMatrix *T = ncm_matrix_new (ncm_model_len (NCM_MODEL (qpw)), ncm_model_len (NCM_MODEL (qpw)));
  NcmVector *v = ncm_vector_new (ncm_model_len (NCM_MODEL (qpw)));
  NcmReparamLinear *relin;
  gint i;
  ncm_matrix_set_zero (T);
  ncm_vector_set_zero (v);

  ncm_matrix_set (T, 0, 0, 1.0);

  for (i = 0; i < qpw->npieces; i++)
  {
    ncm_matrix_set (T, 1 + i, 1, 1.0);
    ncm_matrix_set (T, 1 + i, 2, i*i);
  }

  if (FALSE)
  {
    gint j;
    printf ("npieces: %d\n", qpw->npieces);
    for (i = 0; i < ncm_model_len (NCM_MODEL (qpw)); i++)
    {
      for (j = 0; j < ncm_model_len (NCM_MODEL (qpw)); j++)
        printf ("% -12.2f ",ncm_matrix_get (T, i, j));
      printf ("\n");
    }
    printf("########################\n");
    for (i = 0; i < ncm_model_len (NCM_MODEL (qpw)); i++)
      printf ("% -12.2f ",ncm_vector_get (v, i));
    printf("\n########################\n");
  }

  relin = ncm_reparam_linear_new (ncm_model_len (NCM_MODEL (qpw)), T, v);

  ncm_model_set_reparam (NCM_MODEL (qpw), NCM_REPARAM (relin));

  ncm_vector_free (v);
  ncm_matrix_free (T);

  return;
}

NcHICosmoQPW *
nc_hicosmo_qpw_new (guint npieces, gdouble z_f, gboolean flat)
{
  return g_object_new (NC_TYPE_HICOSMO_QPW, "qp-length", npieces, "zf", z_f, "flat", flat, NULL);
}

enum {
  PROP_0,
  PROP_QP_LENGTH,
  PROP_Z_F,
  PROP_FLAT,
  PROP_SIZE,
};

static void
nc_hicosmo_qpw_init (NcHICosmoQPW *qpw)
{
  qpw->z_f = 0.0;
  qpw->piece = 0.0;
  qpw->npieces = 0;
}

static void
_nc_hicosmo_qpw_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hicosmo_qpw_parent_class)->constructed (object);
  {
    NcHICosmoQPW *qpw = NC_HICOSMO_QPW (object);
    NcmModel *model = NCM_MODEL (qpw);
    NcmModelClass *model_class = NCM_MODEL_GET_CLASS (model);

    qpw->piece = qpw->z_f / qpw->npieces;
    qpw->size = model_class->sparam_len + qpw->npieces;

    return;
  }
}

static void
_nc_hicosmo_qpw_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHICosmoQPW *qpw = NC_HICOSMO_QPW (object);

  g_return_if_fail (NC_IS_HICOSMO_QPW (object));

  switch (prop_id)
  {
    case PROP_QP_LENGTH:
      g_value_set_uint (value, qpw->npieces);
      break;
    case PROP_Z_F:
      g_value_set_double (value, qpw->z_f);
      break;
    case PROP_FLAT:
      g_value_set_boolean (value, qpw->flat);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hicosmo_qpw_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHICosmoQPW *qpw = NC_HICOSMO_QPW (object);
  g_return_if_fail (NC_IS_HICOSMO_QPW (object));

  switch (prop_id)
  {
    case PROP_QP_LENGTH:
      qpw->npieces = g_value_get_double (value);
      break;
    case PROP_Z_F:
      qpw->z_f = g_value_get_double (value);
      break;
    case PROP_FLAT:
      qpw->flat = g_value_get_boolean (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hicosmo_qpw_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qpw_parent_class)->finalize (object);
}

static void
nc_hicosmo_qpw_class_init (NcHICosmoQPWClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcHICosmoClass* parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->constructed  = &_nc_hicosmo_qpw_constructed;
  object_class->finalize     = &nc_hicosmo_qpw_finalize;

  model_class->set_property = &_nc_hicosmo_qpw_set_property;
  model_class->get_property = &_nc_hicosmo_qpw_get_property;

  ncm_model_class_add_params (model_class, NC_HICOSMO_QPW_SPARAM_LEN, NC_HICOSMO_QPW_VPARAM_LEN, PROP_SIZE);
  ncm_model_class_set_name_nick (model_class, "Q piecewise", "qpw");

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QPW_H0, "H_0", "H0",
                              10.0, 500.0, 1.0,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QPW_DEFAULT_H0,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QPW_OMEGA_T, "\\Omega_t", "Omegat",
                              -5.0, 5.0, 1.0e-1,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QPW_DEFAULT_OMEGA_T,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QPW_Q0, "q_0", "q0",
                              -50.0, 50.0, 1.0e-1,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QPW_DEFAULT_Q0,
                              NCM_PARAM_TYPE_FREE);

  ncm_model_class_set_vparam (model_class, NC_HICOSMO_QPW_QP, NC_HICOSMO_QPW_DEFAULT_QP_LEN, "q^\\prime", "qp",
                              -50.0, 50.0, 1.0e-1,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL,
                              NC_HICOSMO_QPW_DEFAULT_QP,
                              NCM_PARAM_TYPE_FREE);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  g_object_class_install_property (object_class,
                                   PROP_Z_F,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "maximum redshift",
                                                        0.0, 100.0, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_FLAT,
                                   g_param_spec_boolean ("flat",
                                                         NULL,
                                                         "set model flat",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  nc_hicosmo_set_H0_impl       (parent_class, &_nc_hicosmo_qpw_H0);
  nc_hicosmo_set_E2_impl       (parent_class, &_nc_hicosmo_qpw_E2);
  nc_hicosmo_set_dE2_dz_impl   (parent_class, &_nc_hicosmo_qpw_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl (parent_class, &_nc_hicosmo_qpw_d2E2_dz2);
  nc_hicosmo_set_Omega_t_impl  (parent_class, &_nc_hicosmo_qpw_Omega_t);
}
