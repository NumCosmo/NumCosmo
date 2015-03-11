/***************************************************************************
 *            nc_hicosmo_qconst.c
 *
 *  Wed Jul 11 14:31:13 2007
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
 * SECTION:nc_hicosmo_qconst
 * @title: Constant Deceleration Parameter Model
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_qconst.h"
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>

G_DEFINE_TYPE (NcHICosmoQConst, nc_hicosmo_qconst, NC_TYPE_HICOSMO);

#define VECTOR   (model->params)
#define MACRO_H0 (ncm_vector_get (VECTOR, NC_HICOSMO_QCONST_H0))
#define OMEGA_T  (ncm_vector_get (VECTOR, NC_HICOSMO_QCONST_OMEGA_T))
#define CD       (ncm_vector_get (VECTOR, NC_HICOSMO_QCONST_CD))
#define E        (ncm_vector_get (VECTOR, NC_HICOSMO_QCONST_E))
#define Q        (ncm_vector_get (VECTOR, NC_HICOSMO_QCONST_Q))
#define Z1       (ncm_vector_get (VECTOR, NC_HICOSMO_QCONST_Z1))

static gdouble
_nc_hicosmo_qconst_cd (NcmModel *model, gdouble z)
{
  gdouble x1, x, ln_x_x1;
  x1 = 1.0 + Z1;
  x = 1.0 + z;
  ln_x_x1 = gsl_sf_log (x/x1);

  if (Z1 == z)
	return CD;

  return CD + x1 * ln_x_x1 / E * gsl_sf_exprel (-Q * ln_x_x1);
}

/****************************************************************************
 * Hubble constant
 ****************************************************************************/
static gdouble _nc_hicosmo_qconst_H0 (NcmModel *model) { return MACRO_H0; }
static gdouble _nc_hicosmo_qconst_Omega_t (NcmModel *model) { return OMEGA_T; }

/**
 * nc_hicosmo_qconst_new:
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcHICosmoQConst *
nc_hicosmo_qconst_new (void)
{
  NcHICosmoQConst *qconst = g_object_new (NC_TYPE_HICOSMO_QCONST, NULL);
  return qconst;
}

enum {
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_qconst_init (NcHICosmoQConst *qconst)
{
  NCM_UNUSED (qconst);
}

static void
nc_hicosmo_qconst_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_qconst_parent_class)->finalize (object);
}

static void
nc_hicosmo_qconst_class_init (NcHICosmoQConstClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcHICosmoClass* parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize     = &nc_hicosmo_qconst_finalize;

  ncm_model_class_set_name_nick (model_class, "Q Constant", "qconst");
  ncm_model_class_add_params (model_class, 6, 0, PROP_SIZE);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QCONST_H0, "H_0", "H0",
                               10.0, 500.0, 1.0,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QCONST_DEFAULT_H0,
                               NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QCONST_OMEGA_T, "\\Omega_t", "Omegat",
                               -5.0, 5.0, 1.0e-1,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QCONST_DEFAULT_OMEGA_T,
                               NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QCONST_CD, "d_c", "cd",
                               -50.0, 50.0, 1.0e-1,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QCONST_DEFAULT_CD,
                               NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QCONST_E, "E", "E",
                               0.0, 50.0, 1.0e-1,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QCONST_DEFAULT_E,
                               NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QCONST_Q, "q", "q",
                               -50.0, 50.0, 1.0e-1,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QCONST_DEFAULT_Q,
                               NCM_PARAM_TYPE_FREE);

  ncm_model_class_set_sparam (model_class, NC_HICOSMO_QCONST_Z1, "z_\\star", "zs",
                               0.0, 5.0, 1.0e-1,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_QCONST_DEFAULT_Z1,
                               NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hicosmo_set_H0_impl (parent_class, &_nc_hicosmo_qconst_H0);
  nc_hicosmo_set_cd_impl (parent_class, &_nc_hicosmo_qconst_cd);
  nc_hicosmo_set_Omega_t_impl (parent_class, &_nc_hicosmo_qconst_Omega_t);
}
