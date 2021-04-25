/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_data_rosenbrock.c
 *
 *  Sat April 17 11:11:28 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_data_rosenbrock.c
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_data_rosenbrock
 * @title: NcmDataRosenbrock
 * @short_description: Rosenbrock distribution.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_data_rosenbrock.h"
#include "math/ncm_model_rosenbrock.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmDataRosenbrockPrivate
{
  gint unused;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmDataRosenbrock, ncm_data_rosenbrock, NCM_TYPE_DATA);

static void
ncm_data_rosenbrock_init (NcmDataRosenbrock *drb)
{
  drb->priv = ncm_data_rosenbrock_get_instance_private (drb);
}

static void
ncm_data_rosenbrock_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_data_rosenbrock_parent_class)->constructed (object);

  ncm_data_set_init (NCM_DATA (object), TRUE);
}

static void
ncm_data_rosenbrock_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_rosenbrock_parent_class)->finalize (object);
}

static guint _ncm_data_rosenbrock_get_length (NcmData *data);
static guint _ncm_data_rosenbrock_get_dof (NcmData *data);
static void _ncm_data_rosenbrock_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);

static void
ncm_data_rosenbrock_class_init (NcmDataRosenbrockClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->constructed = ncm_data_rosenbrock_constructed;
  object_class->finalize    = ncm_data_rosenbrock_finalize;

  data_class->get_length = &_ncm_data_rosenbrock_get_length;
  data_class->get_dof    = &_ncm_data_rosenbrock_get_dof;
  data_class->m2lnL_val  = &_ncm_data_rosenbrock_m2lnL_val;
}

static guint _ncm_data_rosenbrock_get_length (NcmData *data) { return 10; }
static guint _ncm_data_rosenbrock_get_dof (NcmData *data) { return 10; }

static void
_ncm_data_rosenbrock_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcmModelRosenbrock *mrb = NCM_MODEL_ROSENBROCK (ncm_mset_peek (mset, ncm_model_rosenbrock_id ()));
  const gdouble x1 = ncm_model_param_get (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X1);
  const gdouble x2 = ncm_model_param_get (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X2);

  m2lnL[0] = (/*100.0 * */gsl_pow_2 (x2 - x1 * x1) + gsl_pow_2 (1.0 - x1)) * 1.0e-1;
}

/**
 * ncm_data_rosenbrock_new:
 * 
 * Creates a new @dim-dimensional MVND.
 * 
 * Returns: the newly created object.
 */ 
NcmDataRosenbrock *
ncm_data_rosenbrock_new (void)
{
  NcmDataRosenbrock *drb = g_object_new (NCM_TYPE_DATA_ROSENBROCK,
                                         NULL);
  return drb;
}

/**
 * ncm_data_rosenbrock_ref:
 * @drb: a #NcmDataRosenbrock
 * 
 * Increases the reference count of @drb by one.
 * 
 * Returns: (transfer full): @drb
 */
NcmDataRosenbrock *
ncm_data_rosenbrock_ref (NcmDataRosenbrock *drb)
{
  return g_object_ref (drb);
}

/**
 * ncm_data_rosenbrock_free:
 * @drb: a #NcmDataRosenbrock
 * 
 * Decreases the reference count of @drb by one.
 * 
 */
void 
ncm_data_rosenbrock_free (NcmDataRosenbrock *drb)
{
  g_object_unref (drb);
}

/**
 * ncm_data_rosenbrock_clear:
 * @drb: a #NcmDataRosenbrock
 * 
 * If @drb is different from NULL, decreases the reference count of
 * @drb by one and sets @drb to NULL.
 * 
 */
void 
ncm_data_rosenbrock_clear (NcmDataRosenbrock **drb)
{
  g_clear_object (drb);
}
