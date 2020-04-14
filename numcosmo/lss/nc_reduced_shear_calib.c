/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_reduced_shear_calib.c
 *
 *  Tue December 04 23:40:36 2018
 *  Copyright  2018  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_reduced_shear_calib.c
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
 * SECTION:nc_reduced_shear_calib
 * @title: NcReducedShearCalib
 * @short_description: Reduced Shear Calibration 
 *
 * Abstract model for the reduced shear calibration.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_reduced_shear_calib.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_SIZE,
};

struct _NcReducedShearCalibPrivate
{
  gint a;
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcReducedShearCalib, nc_reduced_shear_calib, NCM_TYPE_MODEL);

static void
nc_reduced_shear_calib_init (NcReducedShearCalib *rs_calib)
{
  rs_calib->priv = nc_reduced_shear_calib_get_instance_private (rs_calib);
}

static void
_nc_reduced_shear_calib_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /* NcReducedShearCalib *rs_calib = NC_REDUCED_SHEAR_CALIB (object); */
  g_return_if_fail (NC_IS_REDUCED_SHEAR_CALIB (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_reduced_shear_calib_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcReducedShearCalib *rs_calib = NC_REDUCED_SHEAR_CALIB (object); */
  g_return_if_fail (NC_IS_REDUCED_SHEAR_CALIB (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_reduced_shear_calib_dispose (GObject *object)
{
  /* NcReducedShearCalib *rs_calib = NC_REDUCED_SHEAR_CALIB (object); */
    
  /* Chain up : end */
  G_OBJECT_CLASS (nc_reduced_shear_calib_parent_class)->dispose (object);
}

static void
_nc_reduced_shear_calib_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_reduced_shear_calib_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_reduced_shear_calib, NC_TYPE_REDUCED_SHEAR_CALIB);

static gdouble _nc_reduced_shear_calib_eval (NcReducedShearCalib *rs_calib, const gdouble g_th, const gdouble psf_size, const gdouble gal_size) { g_error ("_nc_reduced_shear_calib_eval: not implemented"); return 0.0; }

static void
nc_reduced_shear_calib_class_init (NcReducedShearCalibClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_reduced_shear_calib_set_property;
  model_class->get_property = &_nc_reduced_shear_calib_get_property;
  
  object_class->dispose      = &_nc_reduced_shear_calib_dispose;
  object_class->finalize     = &_nc_reduced_shear_calib_finalize;

  ncm_model_class_set_name_nick (model_class, "ReducShearCalib", "NcReducedShearCalib");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  ncm_mset_model_register_id (model_class,
                              "NcReducedShearCalib",
                              "ReducShearCalib",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_check_params_info (model_class);

  klass->eval = &_nc_reduced_shear_calib_eval;
}

/**
 * nc_reduced_shear_calib_ref:
 * @rs_calib: a #NcReducedShearCalib
 * 
 * Increases the reference count of @rs_calib by one.
 * 
 * Returns: (transfer full): @rs_calib
 */
NcReducedShearCalib *
nc_reduced_shear_calib_ref (NcReducedShearCalib *rs_calib)
{
  return g_object_ref (rs_calib);
}

/**
 * nc_reduced_shear_calib_free:
 * @rs_calib: a #NcReducedShearCalib
 * 
 * Decreases the reference count of @rs_calib by one.
 * 
 */
void 
nc_reduced_shear_calib_free (NcReducedShearCalib *rs_calib)
{
  g_object_unref (rs_calib);
}

/**
 * nc_reduced_shear_calib_clear:
 * @rs_calib: a #NcReducedShearCalib
 * 
 * If @rs_calib is different from NULL, decreases the reference count of 
 * @rs_calib by one and sets @rs_calib to NULL.
 * 
 */
void 
nc_reduced_shear_calib_clear (NcReducedShearCalib **rs_calib)
{
  g_clear_object (rs_calib);
}

/**
 * nc_reduced_shear_calib_eval: (virtual eval)
 * @rs_calib: a #NcReducedShearCalib
 * @g_th: FIXME
 * @psf_size: FIXME
 * @gal_size: FIXME
 * 
 * 
 * Returns: FIXME
 */
gdouble 
nc_reduced_shear_calib_eval (NcReducedShearCalib *rs_calib, const gdouble g_th, const gdouble psf_size, const gdouble gal_size)
{
  return NC_REDUCED_SHEAR_CALIB_GET_CLASS (rs_calib)->eval (rs_calib, g_th, psf_size, gal_size);
}
