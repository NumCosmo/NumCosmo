/***************************************************************************
 *            nc_hicosmo_qlinear.c
 *
 *  Tue Jun  5 14:25:00 2007
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
 * SECTION:nc_hicosmo_qlinear
 * @title: Linear Desceleration Parameter Model
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_qlinear.h"
#include "math/integral.h"
#include "math/memory_pool.h"
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_exp.h>

G_DEFINE_TYPE (NcHICosmoQLinear, nc_hicosmo_qlinear, NC_TYPE_HICOSMO);

#define VECTOR   (model->params)
#define MACRO_H0 (ncm_vector_get (VECTOR, NC_HICOSMO_QLINEAR_H0))
#define OMEGA_T  (ncm_vector_get (VECTOR, NC_HICOSMO_QLINEAR_OMEGA_T))
#define QLIN_CD  (ncm_vector_get (VECTOR, NC_HICOSMO_QLINEAR_CD))
#define QLIN_E   (ncm_vector_get (VECTOR, NC_HICOSMO_QLINEAR_E))
#define QLIN_Q   (ncm_vector_get (VECTOR, NC_HICOSMO_QLINEAR_Q))
#define QLIN_QP  (ncm_vector_get (VECTOR, NC_HICOSMO_QLINEAR_QP))
#define QLIN_Z1  (ncm_vector_get (VECTOR, NC_HICOSMO_QLINEAR_Z1))

static gdouble
cd_th_int_integrand (gdouble z, gpointer p)
{
  gdouble *params = (gdouble *)p;
  gdouble z1 = params[0];
  gdouble q  = params[1];
  gdouble qp  = params[2];
  gdouble qq = qq = (1.0 + z1)*qp - q;
  gsl_sf_result result;

  gsl_sf_exp_mult_e (qp * (z1 - z), pow ((1.0 + z)/(1.0 + z1), (qq - 1.0) ), &result);

  return result.val;
}

static gdouble
cd_th_int (gdouble z2, gdouble z1, gdouble E, gdouble q, gdouble qp)
{
  gsl_integration_workspace **w;
  gsl_function F;
  gdouble res, err;
  gdouble p[] = {z1, q, qp};

  F.function = &cd_th_int_integrand;
  F.params = p;

  w = ncm_integral_get_workspace ();
  gsl_integration_qag (&F, z1, z2, NCM_INTEGRAL_ERROR, 0.0, NCM_INTEGRAL_PARTITION, NCM_INTEGRAL_ALG, *w, &res, &err);
  ncm_memory_pool_return (w);

  return res/E;
}

/****************************************************************************
 * Comoving Distance
 ****************************************************************************/
static gdouble
_nc_hicosmo_qlinear_cd (NcmModel *model, gdouble z)
{
  if (QLIN_Z1 == z)
	return QLIN_CD;
  return QLIN_CD + cd_th_int (z, QLIN_Z1, QLIN_E, QLIN_Q, QLIN_QP);
}

/**
 * nc_hicosmo_qlinear_dE:
 * @z2: FIXME
 * @z1: FIXME
 * @q: FIXME
 * @qp: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hicosmo_qlinear_dE (gdouble z2, gdouble z1, gdouble q, gdouble qp)
{
  gdouble x1 = 1.0 + z1;
  gdouble x2 = 1.0 + z2;
  gdouble delta = z2 - z1;
  gdouble qq = x1 * qp - q;
  return gsl_sf_exp_mult (qp * delta, pow (x2 / x1, 1.0 - qq));
}

/****************************************************************************
 * Hubble constant
 ****************************************************************************/
static gdouble _nc_hicosmo_qlinear_H0 (NcmModel *model) { return MACRO_H0; }
static gdouble _nc_hicosmo_qlinear_Omega_t (NcmModel *model) { return OMEGA_T; }

/**
 * nc_hicosmo_qlinear_new:
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcHICosmoQLinear *
nc_hicosmo_qlinear_new (void)
{
  NcHICosmoQLinear *qlinear = g_object_new (NC_TYPE_HICOSMO_QLINEAR, NULL);
  return qlinear;
}

enum {
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_qlinear_init (NcHICosmoQLinear *qlinear)
{
}

static void
nc_hicosmo_qlinear_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qlinear_parent_class)->finalize (object);
}

static void
nc_hicosmo_qlinear_class_init (NcHICosmoQLinearClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcHICosmoClass* parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize     = &nc_hicosmo_qlinear_finalize;

  ncm_model_class_add_params (model_class, 7, 0, PROP_SIZE);
  ncm_model_class_set_name_nick (model_class, "Q Linear", "qlinear");

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QLINEAR_H0, "H_0", "H0",
                               10.0, 500.0, 1.0,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QLINEAR_DEFAULT_H0,
                               NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QLINEAR_OMEGA_T, "\\Omega_t", "Omegat",
                               -5.0, 5.0, 1.0e-1,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QLINEAR_DEFAULT_OMEGA_T,
                               NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QLINEAR_CD, "d_c", "cd",
                               -50.0, 50.0, 1.0e-1,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QLINEAR_DEFAULT_CD,
                               NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QLINEAR_E, "E", "E",
                               0.0, 50.0, 1.0e-1,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QLINEAR_DEFAULT_E,
                               NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QLINEAR_Q, "q", "q",
                               -50.0, 50.0, 1.0e-1,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QLINEAR_DEFAULT_Q,
                               NCM_PARAM_TYPE_FREE);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QLINEAR_QP, "q^\\prime", "qp",
                               -50.0, 50.0, 1.0e-1,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QLINEAR_DEFAULT_QP,
                               NCM_PARAM_TYPE_FREE);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QLINEAR_Z1, "z_\\star", "zs",
                               0.0, 5.0, 1.0e-1,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QLINEAR_DEFAULT_Z1,
                               NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hicosmo_set_H0_impl (parent_class, &_nc_hicosmo_qlinear_H0);
  nc_hicosmo_set_cd_impl (parent_class, &_nc_hicosmo_qlinear_cd);
  nc_hicosmo_set_Omega_t_impl (parent_class, &_nc_hicosmo_qlinear_Omega_t);
}
