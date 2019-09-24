/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_reduced_shear_calib_wtg.c
 *
 *  Tue December 04 23:47:36 2018
 *  Copyright  2018  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_reduced_shear_calib_wtg.c
 * Copyright (C) 2018 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_reduced_shear_calib_wtg
 * @title: NcReducedShearCalibWtg
 * @short_description: Reduced Shear Calibration 
 *
 * Abstract model for the reduced shear calibration.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_reduced_shear_calib_wtg.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_SIZE,
};

struct _NcReducedShearCalibWtgPrivate
{
  gdouble m_slope;
  gdouble m_b;
  gdouble c;
  gdouble size_ratio;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcReducedShearCalibWtg, nc_reduced_shear_calib_wtg, NC_TYPE_REDUCED_SHEAR_CALIB);

static void
nc_reduced_shear_calib_wtg_init (NcReducedShearCalibWtg *rs_wtg)
{
  NcReducedShearCalibWtgPrivate * const self = rs_wtg->priv = nc_reduced_shear_calib_wtg_get_instance_private (rs_wtg);

  self->m_slope    = 0.0;
  self->m_b        = 0.0;
  self->c          = 0.0;
  self->size_ratio = 0.0;
}

static void
_nc_reduced_shear_calib_wtg_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /* NcReducedShearCalibWtg *rs_wtg = NC_REDUCED_SHEAR_CALIB_WTG (object); */
  g_return_if_fail (NC_IS_REDUCED_SHEAR_CALIB_WTG (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_reduced_shear_calib_wtg_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcReducedShearCalibWtg *rs_wtg = NC_REDUCED_SHEAR_CALIB_WTG (object); */
  g_return_if_fail (NC_IS_REDUCED_SHEAR_CALIB_WTG (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_reduced_shear_calib_wtg_dispose (GObject *object)
{
  /* NcReducedShearCalibWtg *rs_wtg = NC_REDUCED_SHEAR_CALIB_WTG (object); */
    
  /* Chain up : end */
  G_OBJECT_CLASS (nc_reduced_shear_calib_wtg_parent_class)->dispose (object);
}

static void
_nc_reduced_shear_calib_wtg_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_reduced_shear_calib_wtg_parent_class)->finalize (object);
}

static gdouble _nc_reduced_shear_calib_wtg_eval (NcReducedShearCalib *rs_calib, const gdouble g_th, const gdouble psf_size, const gdouble gal_size);

static void
nc_reduced_shear_calib_wtg_class_init (NcReducedShearCalibWtgClass *klass)
{
  GObjectClass* object_class               = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class               = NCM_MODEL_CLASS (klass);
  NcReducedShearCalibClass *rs_calib_class = NC_REDUCED_SHEAR_CALIB_CLASS (klass);

  model_class->set_property = &_nc_reduced_shear_calib_wtg_set_property;
  model_class->get_property = &_nc_reduced_shear_calib_wtg_get_property;
  
  object_class->dispose      = &_nc_reduced_shear_calib_wtg_dispose;
  object_class->finalize     = &_nc_reduced_shear_calib_wtg_finalize;

  ncm_model_class_set_name_nick (model_class, "RedShearWtG", "NcReducedShearCalibWtg");
  ncm_model_class_add_params (model_class, NNC_REDUCED_SHEAR_CALIB_WTG_SPARAM_LEN, 0, PROP_SIZE);

  ncm_model_class_set_sparam (model_class, NC_REDUCED_SHEAR_CALIB_WTG_MSLOPE, "m_s", "mslope",
                              0.01, 2.0, 0.05, 0.0, 0.2, NCM_PARAM_TYPE_FREE);
  ncm_model_class_set_sparam (model_class, NC_REDUCED_SHEAR_CALIB_WTG_MB, "m_b", "mb",
                              -0.5, 0.5, 0.05, 0.0, -0.028, NCM_PARAM_TYPE_FREE);
  ncm_model_class_set_sparam (model_class, NC_REDUCED_SHEAR_CALIB_WTG_C, "c", "c",
                              -1.0e-3, 1.0e-3, 1.0e-4, 0.0, -1.0e-5, NCM_PARAM_TYPE_FREE);
  ncm_model_class_set_sparam (model_class, NC_REDUCED_SHEAR_CALIB_WTG_SIZE_RATIO, "x_p", "xp",
                              1.0, 3.0, 0.1, 0.0, 1.97, NCM_PARAM_TYPE_FREE);
  
  ncm_model_class_check_params_info (model_class);

  rs_calib_class->eval = &_nc_reduced_shear_calib_wtg_eval;
}

#define VECTOR (NCM_MODEL (rs_wtg)->params)
#define MSLOPE     (ncm_vector_get (VECTOR, NC_REDUCED_SHEAR_CALIB_WTG_MSLOPE))
#define MB         (ncm_vector_get (VECTOR, NC_REDUCED_SHEAR_CALIB_WTG_MB))
#define C          (ncm_vector_get (VECTOR, NC_REDUCED_SHEAR_CALIB_WTG_C))
#define SIZE_RATIO (ncm_vector_get (VECTOR, NC_REDUCED_SHEAR_CALIB_WTG_SIZE_RATIO))


static gdouble 
_nc_reduced_shear_calib_wtg_eval (NcReducedShearCalib *rs_calib, const gdouble g_th, const gdouble psf_size, const gdouble gal_size)
{
  NcReducedShearCalibWtg *rs_wtg = NC_REDUCED_SHEAR_CALIB_WTG (rs_calib);
  const gdouble rh_psf_size = gal_size / psf_size;
  gdouble m;
  
  if (rh_psf_size >= SIZE_RATIO)
	{
		m = MB;
	}
	else
	{
		m = MSLOPE * rh_psf_size + MB;
	}
  
	return (1.0 + m) * g_th + C;
}

/**
 * nc_reduced_shear_calib_wtg_new:
 * 
 * Creates a new MVND mean model of @dim dimensions.
 * 
 * Returns: (transfer full): the newly created #NcReducedShearCalibWtg
 */
NcReducedShearCalibWtg *
nc_reduced_shear_calib_wtg_new (void)
{
  NcReducedShearCalibWtg *rs_wtg = g_object_new (NC_TYPE_REDUCED_SHEAR_CALIB_WTG,
                                                 NULL);
  return rs_wtg;
}

/**
 * nc_reduced_shear_calib_wtg_ref:
 * @rs_wtg: a #NcReducedShearCalibWtg
 * 
 * Increases the reference count of @rs_wtg by one.
 * 
 * Returns: (transfer full): @rs_wtg
 */
NcReducedShearCalibWtg *
nc_reduced_shear_calib_wtg_ref (NcReducedShearCalibWtg *rs_wtg)
{
  return g_object_ref (rs_wtg);
}

/**
 * nc_reduced_shear_calib_wtg_free:
 * @rs_wtg: a #NcReducedShearCalibWtg
 * 
 * Decreases the reference count of @rs_wtg by one.
 * 
 */
void 
nc_reduced_shear_calib_wtg_free (NcReducedShearCalibWtg *rs_wtg)
{
  g_object_unref (rs_wtg);
}

/**
 * nc_reduced_shear_calib_wtg_clear:
 * @rs_wtg: a #NcReducedShearCalibWtg
 * 
 * If @rs_wtg is different from NULL, decreases the reference count of 
 * @rs_wtg by one and sets @rs_wtg to NULL.
 * 
 */
void 
nc_reduced_shear_calib_wtg_clear (NcReducedShearCalibWtg **rs_wtg)
{
  g_clear_object (rs_wtg);
}
