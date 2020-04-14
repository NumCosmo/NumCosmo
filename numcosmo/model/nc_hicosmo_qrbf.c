/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hicosmo_qrbf.c
 *
 *  Fri November 01 14:18:09 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hicosmo_qrbf.c
 *
 * Copyright (C) 2019 - Sandro Dias Pinto Vitenti
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * SECTION:nc_hicosmo_qrbf
 * @title: NcHICosmoQRBF
 * @short_description: Kinetic model -- Spline deceleration function
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_qrbf.h"
#include "math/integral.h"
#include "math/ncm_memory_pool.h"


#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_fit.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcHICosmoQRBFPrivate 
{
  NcmVector *centers;
  NcmVector *coeffs;
  NcmVector *int_z0;
  guint rbf_size;
  gdouble h;
  gdouble h2;
  gdouble z_f;
  gboolean constructed;
};

enum {
  PROP_0,
  PROP_Z_F,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHICosmoQRBF, nc_hicosmo_qrbf, NC_TYPE_HICOSMO);

static void
nc_hicosmo_qrbf_init (NcHICosmoQRBF *qrbf)
{
  NcHICosmoQRBFPrivate * const self = qrbf->priv = nc_hicosmo_qrbf_get_instance_private (qrbf);
  self->centers     = NULL;
  self->coeffs      = NULL;
  self->int_z0      = NULL;
  self->rbf_size    = 0;
  self->h           = 0.0;
  self->h2          = 0.0;
  self->z_f         = 0.0;
  self->constructed = FALSE;
}

static void
_nc_hicosmo_qrbf_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hicosmo_qrbf_parent_class)->constructed (object);
  {
    NcHICosmoQRBF *qrbf = NC_HICOSMO_QRBF (object);
    NcHICosmoQRBFPrivate * const self = qrbf->priv;
    NcmModel *model      = NCM_MODEL (qrbf);
    const guint rbf_size = ncm_model_vparam_len (model, NC_HICOSMO_QRBF_RBF_CENTERS);
    guint centers_i, coeffs_i;
    gint i;

    g_assert_cmpuint (rbf_size, ==, ncm_model_vparam_len (model, NC_HICOSMO_QRBF_RBF_COEFFS));

    centers_i = ncm_model_vparam_index (model, NC_HICOSMO_QRBF_RBF_CENTERS, 0);
    coeffs_i  = ncm_model_vparam_index (model, NC_HICOSMO_QRBF_RBF_COEFFS, 0);    

    self->centers  = ncm_vector_get_subvector (model->params, centers_i, rbf_size);
    self->coeffs   = ncm_vector_get_subvector (model->params, coeffs_i, rbf_size);
    self->int_z0   = ncm_vector_dup (self->centers);
    self->rbf_size = rbf_size;

    for (i = 0; i < rbf_size; i++)
    {
      ncm_vector_set (self->centers, i, self->z_f * (i + 0.0) / (0.0 + rbf_size));
      ncm_model_param_set_lower_bound (model, centers_i + i, 0.0);
      ncm_model_param_set_upper_bound (model, centers_i + i, self->z_f);
    }

    ncm_model_state_mark_outdated (model);
    self->constructed = TRUE;
  }
}

static void
_nc_hicosmo_qrbf_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHICosmoQRBF *qrbf = NC_HICOSMO_QRBF (object);
  NcHICosmoQRBFPrivate * const self = qrbf->priv;

  g_return_if_fail (NC_IS_HICOSMO_QRBF (object));

  switch (prop_id)
  {
    case PROP_Z_F:
      g_value_set_double (value, self->z_f);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hicosmo_qrbf_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHICosmoQRBF *qrbf = NC_HICOSMO_QRBF (object);
  /*NcHICosmoQRBFPrivate * const self = qrbf->priv;*/
  g_return_if_fail (NC_IS_HICOSMO_QRBF (object));

  switch (prop_id)
  {
    case PROP_Z_F:
      nc_hicosmo_qrbf_set_z_f (qrbf, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hicosmo_qrbf_dispose (GObject *object)
{
  NcHICosmoQRBF *qrbf = NC_HICOSMO_QRBF (object);
  NcHICosmoQRBFPrivate * const self = qrbf->priv;

  ncm_vector_clear (&self->coeffs);
  ncm_vector_clear (&self->centers);
  ncm_vector_clear (&self->int_z0);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qrbf_parent_class)->dispose (object);
}

static void
_nc_hicosmo_qrbf_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qrbf_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_qrbf_H0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_qrbf_E2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_qrbf_dE2_dz (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_qrbf_d2E2_dz2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_qrbf_Omega_t0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_qrbf_as_drag (NcHICosmo *cosmo);

static void
nc_hicosmo_qrbf_class_init (NcHICosmoQRBFClass *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass* parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_hicosmo_qrbf_set_property;
  model_class->get_property = &_nc_hicosmo_qrbf_get_property;
  object_class->constructed = &_nc_hicosmo_qrbf_constructed;
  object_class->dispose     = &_nc_hicosmo_qrbf_dispose;
  object_class->finalize    = &_nc_hicosmo_qrbf_finalize;

  ncm_model_class_set_name_nick (model_class, "Q RBF", "qrbf");
  ncm_model_class_add_params (model_class, NC_HICOSMO_QRBF_SPARAM_LEN, NC_HICOSMO_QRBF_VPARAM_LEN, PROP_SIZE);

  /**
   * NcHICosmoQRBF:H0:
   *
   * Hubble parameter today (z = 0).
   */
  /**
   * NcHICosmoQRBF:H0-fit:
   *
   * FIXME
   */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QRBF_H0, "H_0", "H0",
                              10.0, 500.0, 1.0, NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QRBF_DEFAULT_H0,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcHICosmoQRBF:Omegat:
   *
   * FIXME
   */
  /**
   * NcHICosmoQRBF:Omegat-fit:
   *
   * FIXME
   */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QRBF_OMEGA_T, "Omega_t0", "Omegat",
                              0.05, 2.0, 1.0e-1,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QRBF_DEFAULT_OMEGA_T,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QRBF_AS_DRAG, "A_s", "asdrag",
                              1.0e-4, 1.0, 1.0e-3,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QRBF_DEFAULT_AS_DRAG,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QRBF_RBF_H, "h_r", "hr",
                              1.0e-1, 1.0e1, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QRBF_DEFAULT_RBF_H,
                              NCM_PARAM_TYPE_FREE);
  
  ncm_model_class_set_vparam (model_class, NC_HICOSMO_QRBF_RBF_CENTERS, NC_HICOSMO_QRBF_DEFAULT_RBF_CENTERS_LEN, "x_i", "xi",
                              0.0, 1.0, 1.0e-1, NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QRBF_DEFAULT_RBF_CENTERS,
                              NCM_PARAM_TYPE_FREE);

  ncm_model_class_set_vparam (model_class, NC_HICOSMO_QRBF_RBF_COEFFS, NC_HICOSMO_QRBF_DEFAULT_RBF_COEFFS_LEN, "c_i", "ci",
                              -1.0e2, 1.0e2, 1.0e-1, NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QRBF_DEFAULT_RBF_COEFFS,
                              NCM_PARAM_TYPE_FREE);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  g_object_class_install_property (object_class,
                                   PROP_Z_F,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "final redshift",
                                                        0.0, 100.0, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  nc_hicosmo_set_H0_impl       (parent_class, &_nc_hicosmo_qrbf_H0);
  nc_hicosmo_set_E2_impl       (parent_class, &_nc_hicosmo_qrbf_E2);
  nc_hicosmo_set_dE2_dz_impl   (parent_class, &_nc_hicosmo_qrbf_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl (parent_class, &_nc_hicosmo_qrbf_d2E2_dz2);
  nc_hicosmo_set_Omega_t0_impl (parent_class, &_nc_hicosmo_qrbf_Omega_t0);
  nc_hicosmo_set_as_drag_impl  (parent_class, &_nc_hicosmo_qrbf_as_drag);
}

#define VECTOR     (NCM_MODEL (cosmo)->params)
#define QRBF_H0 (ncm_vector_get (VECTOR, NC_HICOSMO_QRBF_H0))
#define OMEGA_T    (ncm_vector_get (VECTOR, NC_HICOSMO_QRBF_OMEGA_T))
#define AS_DRAG    (ncm_vector_get (VECTOR, NC_HICOSMO_QRBF_AS_DRAG))
#define HRBF       (ncm_vector_get (VECTOR, NC_HICOSMO_QRBF_RBF_H))

static gdouble 
_nc_hicosmo_qrbf_f (const gdouble z, const gdouble zi, const gdouble h)
{
  const gdouble u = h * (z - zi);

  return 1.0 / (1.0 + u * u);
}

static gdouble 
_nc_hicosmo_qrbf_df (const gdouble z, const gdouble zi, const gdouble h)
{
  const gdouble u = h * (z - zi);
  
  return - 2.0 * h * u / gsl_pow_2 (1.0 + u * u);
}

static gdouble 
_nc_hicosmo_qrbf_d2f (const gdouble z, const gdouble zi, const gdouble h)
{
  const gdouble u  = h * (z - zi);
  const gdouble u2 = u * u;
  
  return 2.0 * h * h * (3.0 * u2 - 1.0) / gsl_pow_3 (1.0 + u2);
}

static gdouble 
_nc_hicosmo_qrbf_int_2f_1pz (const gdouble z, const gdouble zi, const gdouble h)
{
  const gdouble w  = h * (1.0 + zi);
  const gdouble u  = h * (z - zi);
  const gdouble x  = 1.0 + z;
  
  return (2.0 * w * atan (u) + log (x * x / (1.0 + u * u))) / (1.0 + w * w);
}

static void
_nc_hicosmo_qrbf_prepare (NcHICosmoQRBF *qrbf)
{
  if (!ncm_model_state_is_update (NCM_MODEL (qrbf)))
  {
    NcHICosmoQRBFPrivate * const self = qrbf->priv;
    NcHICosmoQRBF *cosmo = qrbf;
    gint i;
    
    /*printf ("# Preparing!\n");*/

    self->h = HRBF;

    for (i = 0; i < self->rbf_size; i++)
    {
      const gdouble z_i      = ncm_vector_get (self->centers, i);
      const gdouble int_z0_i = _nc_hicosmo_qrbf_int_2f_1pz (0.0, z_i, self->h);
      
      ncm_vector_set (self->int_z0, i, int_z0_i);
    }

    ncm_model_state_set_update (NCM_MODEL (qrbf));
  }
  else
    return;
}

static gdouble
_nc_hicosmo_qrbf_q (NcHICosmo *cosmo, gdouble z)
{
  NcHICosmoQRBF *qrbf = NC_HICOSMO_QRBF (cosmo);
  NcHICosmoQRBFPrivate * const self = qrbf->priv;
  gdouble qp1 = 0.0;
  gint i;

  _nc_hicosmo_qrbf_prepare (qrbf);
  
  for (i = 0; i < self->rbf_size; i++)
  {
    const gdouble z_i     = ncm_vector_get (self->centers, i);
    const gdouble coeff_i = ncm_vector_get (self->coeffs, i);
    
    qp1 += coeff_i * _nc_hicosmo_qrbf_f (z, z_i, self->h);
  }

  /*printf ("# q % 22.15g % 22.15g\n", z, qp1 - 1.0);*/
  return qp1 - 1.0;
}

static gdouble
_nc_hicosmo_qrbf_dq_dz (NcHICosmo *cosmo, gdouble z)
{
  NcHICosmoQRBF *qrbf = NC_HICOSMO_QRBF (cosmo);
  NcHICosmoQRBFPrivate * const self = qrbf->priv;
  gdouble dq = 0.0;
  gint i;

  _nc_hicosmo_qrbf_prepare (qrbf);
  
  for (i = 0; i < self->rbf_size; i++)
  {
    const gdouble z_i     = ncm_vector_get (self->centers, i);
    const gdouble coeff_i = ncm_vector_get (self->coeffs, i);
    
    dq += coeff_i * _nc_hicosmo_qrbf_df (z, z_i, self->h);
  }
  
  return dq;
}

static gdouble
_nc_hicosmo_qrbf_d2q_dz2 (NcHICosmo *cosmo, gdouble z)
{
  NcHICosmoQRBF *qrbf = NC_HICOSMO_QRBF (cosmo);
  NcHICosmoQRBFPrivate * const self = qrbf->priv;
  gdouble d2q = 0.0;
  gint i;

  _nc_hicosmo_qrbf_prepare (qrbf);
  
  for (i = 0; i < self->rbf_size; i++)
  {
    const gdouble z_i     = ncm_vector_get (self->centers, i);
    const gdouble coeff_i = ncm_vector_get (self->coeffs, i);
    
    d2q += coeff_i * _nc_hicosmo_qrbf_d2f (z, z_i, self->h);
  }
  
  return d2q;
}

static gdouble 
_nc_hicosmo_qrbf_d2q2 (gdouble z, gpointer user_data)
{
  NcHICosmo *cosmo = NC_HICOSMO (user_data);
  return gsl_pow_2 (_nc_hicosmo_qrbf_d2q_dz2 (cosmo, z));
}

static gdouble
_nc_hicosmo_qrbf_E2 (NcHICosmo *cosmo, gdouble z)
{
  NcHICosmoQRBF *qrbf = NC_HICOSMO_QRBF (cosmo);
  NcHICosmoQRBFPrivate * const self = qrbf->priv;
  gdouble lnE2 = 0.0;
  gint i;

  _nc_hicosmo_qrbf_prepare (qrbf);

  for (i = 0; i < self->rbf_size; i++)
  {
    const gdouble z_i     = ncm_vector_get (self->centers, i);
    const gdouble int_z_i = _nc_hicosmo_qrbf_int_2f_1pz (z, z_i, self->h);
    const gdouble coeff_i = ncm_vector_get (self->coeffs, i);

    lnE2 += coeff_i * (int_z_i - ncm_vector_get (self->int_z0, i));
    /*printf ("% 22.15g % 22.15g % 22.15g % 22.15g % 22.15g\n", z, z_i, int_z_i, ncm_vector_get (self->int_z0, i), coeff_i);*/
  }

  /*printf ("# lnE2 % 22.15g % 22.15g\n", z, lnE2);*/
  return exp (lnE2);
}

static gdouble
_nc_hicosmo_qrbf_dE2_dz (NcHICosmo *cosmo, gdouble z)
{
  const gdouble q      = _nc_hicosmo_qrbf_q (cosmo, z);
  const gdouble E2     = _nc_hicosmo_qrbf_E2 (cosmo, z);
  const gdouble dE2_dz = 2.0 * E2 * (q + 1.0) / (1.0 + z);

  /*printf ("# dE2_dz % 22.15g % 22.15g % 22.15g % 22.15g\n", z, E2, q, (dE2_dz * (1.0 + z) / (2.0 * E2) - 1.0));*/
  
  return dE2_dz;
}

static gdouble
_nc_hicosmo_qrbf_d2E2_dz2 (NcHICosmo *cosmo, gdouble z)
{  
  const gdouble q   = _nc_hicosmo_qrbf_q (cosmo, z);
  const gdouble dq  = _nc_hicosmo_qrbf_dq_dz (cosmo, z);
  const gdouble E2  = _nc_hicosmo_qrbf_E2 (cosmo, z);
  const gdouble dE2 = _nc_hicosmo_qrbf_dE2_dz (cosmo, z);
  const gdouble x   = 1.0 + z;

  return dE2 * dE2 / E2 + 2.0 * E2 * (dq - q / x) / x;
}

static gdouble _nc_hicosmo_qrbf_H0 (NcHICosmo *cosmo) { return QRBF_H0; }
static gdouble _nc_hicosmo_qrbf_Omega_t0 (NcHICosmo *cosmo) { return OMEGA_T; }
static gdouble _nc_hicosmo_qrbf_as_drag (NcHICosmo *cosmo) { return AS_DRAG; }

/**
 * nc_hicosmo_qrbf_new:
 * @np: number of knots
 * @z_f: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcHICosmoQRBF *
nc_hicosmo_qrbf_new (gsize np, gdouble z_f)
{
  NcHICosmoQRBF *qrbf = g_object_new (NC_TYPE_HICOSMO_QRBF,
                                      "zf", z_f,
                                      "xi-length", np,
                                      "ci-length", np,
                                      NULL);
  return qrbf;
}

/**
 * nc_hicosmo_qrbf_set_z_f:
 * @qrbf: a #NcHICosmoQRBF
 * @z_f: FIXME
 *
 * FIXME
 * 
 */
void 
nc_hicosmo_qrbf_set_z_f (NcHICosmoQRBF *qrbf, const gdouble z_f)
{
  NcHICosmoQRBFPrivate * const self = qrbf->priv;
  
  self->z_f = z_f;

  if (self->constructed)
  {
    NcmModel *model       = NCM_MODEL (qrbf);
    const guint rbf_size  = ncm_model_vparam_len (model, NC_HICOSMO_QRBF_RBF_CENTERS);
    const guint centers_i = ncm_model_vparam_index (model, NC_HICOSMO_QRBF_RBF_CENTERS, 0);

    gint i;
    
    for (i = 0; i < rbf_size; i++)
    {
      ncm_model_param_set_lower_bound (model, centers_i + i, 0.0);
      ncm_model_param_set_upper_bound (model, centers_i + i, z_f);
    }
  }
}

/**
 * nc_hicosmo_qrbf_q_roughness:
 * @qrbf: a #NcHICosmoQRBF
 * 
 * 
 * Returns: the roughness of $q(z)$.
 */
gdouble
nc_hicosmo_qrbf_q_roughness (NcHICosmoQRBF *qrbf)
{
  NcHICosmoQRBFPrivate * const self = qrbf->priv;
  gsl_integration_workspace **w = ncm_integral_get_workspace();
  gsl_function F;
  gdouble rn, err;

  F.function = &_nc_hicosmo_qrbf_d2q2;
  F.params   = qrbf;
  
  gsl_integration_qag (&F, 0.0, self->z_f, 0.0, 1.0e-4, NCM_INTEGRAL_PARTITION, 6, *w, &rn, &err);

  ncm_memory_pool_return (w);

  return rn;
}

/*----------------------------------------------------------------------------*/

struct _NcHICosmoQRBFRpriorPrivate
{
  gdouble lambda;
};

enum
{
  PROP_01,
  PROP_LAMBDA,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHICosmoQRBFRprior, nc_hicosmo_qrbf_rprior, NCM_TYPE_PRIOR);

static void
nc_hicosmo_qrbf_rprior_init (NcHICosmoQRBFRprior *qrprior)
{
  NcHICosmoQRBFRpriorPrivate * const self = qrprior->priv = nc_hicosmo_qrbf_rprior_get_instance_private (qrprior);
  self->lambda = 0.0;
}

static void
_nc_hicosmo_qrbf_rprior_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHICosmoQRBFRprior *qrprior = NC_HICOSMO_QRBF_RPRIOR (object);
  NcHICosmoQRBFRpriorPrivate * const self = qrprior->priv;
  g_return_if_fail (NC_IS_HICOSMO_QRBF_RPRIOR (object));

  switch (prop_id)
  {
    case PROP_LAMBDA:
      self->lambda = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hicosmo_qrbf_rprior_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHICosmoQRBFRprior *qrprior = NC_HICOSMO_QRBF_RPRIOR (object);
  NcHICosmoQRBFRpriorPrivate * const self = qrprior->priv;
  g_return_if_fail (NC_IS_HICOSMO_QRBF_RPRIOR (object));

  switch (prop_id)
  {
    case PROP_LAMBDA:
      g_value_set_double (value, self->lambda);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hicosmo_qrbf_rprior_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qrbf_rprior_parent_class)->finalize (object);
}

static void 
_nc_hicosmo_qrbf_rprior_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  NcHICosmoQRBFRprior *qrprior = NC_HICOSMO_QRBF_RPRIOR (func);
  NcHICosmoQRBFRpriorPrivate * const self = qrprior->priv;
  NcHICosmoQRBF *qrbf = NC_HICOSMO_QRBF (ncm_mset_peek (mset, nc_hicosmo_id ()));

  /*g_assert (qrbf != NULL);*/
  /*g_assert (NC_IS_HICOSMO_QRBF (qrbf));*/

  res[0] = self->lambda * nc_hicosmo_qrbf_q_roughness (qrbf);
}

static void
nc_hicosmo_qrbf_rprior_class_init (NcHICosmoQRBFRpriorClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmMSetFuncClass *mset_func_class = NCM_MSET_FUNC_CLASS (klass);

  object_class->set_property = &_nc_hicosmo_qrbf_rprior_set_property;
  object_class->get_property = &_nc_hicosmo_qrbf_rprior_get_property;
  object_class->finalize     = &_nc_hicosmo_qrbf_rprior_finalize;

  g_object_class_install_property (object_class,
                                   PROP_LAMBDA,
                                   g_param_spec_double ("lambda",
                                                        NULL,
                                                        "\\lambda",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  NCM_PRIOR_CLASS (klass)->is_m2lnL = TRUE;
  mset_func_class->eval             = &_nc_hicosmo_qrbf_rprior_eval;
}

/**
 * nc_hicosmo_qrbf_rprior_new:
 * @lambda: penalization height
 *
 * Returns: FIXME
 */
NcHICosmoQRBFRprior *
nc_hicosmo_qrbf_rprior_new (gdouble lambda)
{
  NcHICosmoQRBFRprior *qrprior = g_object_new (NC_TYPE_HICOSMO_QRBF_RPRIOR,
                                               "lambda", lambda,
                                               NULL);
  return qrprior;
}
