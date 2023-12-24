/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_model_rosenbrock.c
 *
 *  Sat April 17 10:57:36 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_model_rosenbrock.c
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:ncm_model_rosenbrock
 * @title: NcmModelRosenbrock
 * @short_description: Multivariate Normal Distribution mean model.
 *
 * Multivariate Normal distribution model of the mean.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_model_rosenbrock.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_SIZE,
};

typedef struct _NcmModelRosenbrockPrivate
{
  gint place_holder;
} NcmModelRosenbrockPrivate;

struct _NcmModelRosenbrock
{
  NcmModel parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmModelRosenbrock, ncm_model_rosenbrock, NCM_TYPE_MODEL)

static void
ncm_model_rosenbrock_init (NcmModelRosenbrock *model_rosenbrock)
{
  NcmModelRosenbrockPrivate * const self = ncm_model_rosenbrock_get_instance_private (model_rosenbrock);

  self->place_holder = 0;
}

static void
_ncm_model_rosenbrock_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /*NcmModelRosenbrock *model_rosenbrock = NCM_MODEL_ROSENBROCK (object);*/
  g_return_if_fail (NCM_IS_MODEL_ROSENBROCK (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_model_rosenbrock_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /*NcmModelRosenbrock *model_rosenbrock = NCM_MODEL_ROSENBROCK (object);*/
  g_return_if_fail (NCM_IS_MODEL_ROSENBROCK (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_model_rosenbrock_dispose (GObject *object)
{
  /*NcmModelRosenbrock *model_rosenbrock = NCM_MODEL_ROSENBROCK (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_model_rosenbrock_parent_class)->dispose (object);
}

static void
_ncm_model_rosenbrock_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_model_rosenbrock_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (ncm_model_rosenbrock, NCM_TYPE_MODEL_ROSENBROCK);

static void
ncm_model_rosenbrock_class_init (NcmModelRosenbrockClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_ncm_model_rosenbrock_set_property;
  model_class->get_property = &_ncm_model_rosenbrock_get_property;

  object_class->dispose  = &_ncm_model_rosenbrock_dispose;
  object_class->finalize = &_ncm_model_rosenbrock_finalize;

  ncm_model_class_set_name_nick (model_class, "MRB", "NcmModelRosenbrock");
  ncm_model_class_add_params (model_class, NNCM_MODEL_ROSENBROCK_SPARAM_LEN, 0, PROP_SIZE);

  ncm_mset_model_register_id (model_class,
                              "NcmModelRosenbrock",
                              "MRB",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_set_sparam (model_class, NCM_MODEL_ROSENBROCK_X1, "x_1", "x1",
                              -200.0, 200.0, 1.0e-1, 0.0, 0.0, NCM_PARAM_TYPE_FREE);
  ncm_model_class_set_sparam (model_class, NCM_MODEL_ROSENBROCK_X2, "x_2", "x2",
                              -400.0, 800.0, 1.0e-1, 0.0, 0.0, NCM_PARAM_TYPE_FREE);

  ncm_model_class_check_params_info (model_class);
}

/**
 * ncm_model_rosenbrock_new:
 *
 * Creates a new Rosenbrock model.
 *
 * Returns: (transfer full): the newly created #NcmModelRosenbrock
 */
NcmModelRosenbrock *
ncm_model_rosenbrock_new (void)
{
  NcmModelRosenbrock *mrb = g_object_new (NCM_TYPE_MODEL_ROSENBROCK,
                                          NULL);

  return mrb;
}

/**
 * ncm_model_rosenbrock_ref:
 * @mrb: a #NcmModelRosenbrock
 *
 * Increases the reference count of @mrb by one.
 *
 * Returns: (transfer full): @mrb
 */
NcmModelRosenbrock *
ncm_model_rosenbrock_ref (NcmModelRosenbrock *mrb)
{
  return g_object_ref (mrb);
}

/**
 * ncm_model_rosenbrock_free:
 * @mrb: a #NcmModelRosenbrock
 *
 * Decreases the reference count of @mrb by one.
 *
 */
void
ncm_model_rosenbrock_free (NcmModelRosenbrock *mrb)
{
  g_object_unref (mrb);
}

/**
 * ncm_model_rosenbrock_clear:
 * @mrb: a #NcmModelRosenbrock
 *
 * If @mrb is different from NULL, decreases the reference count of
 * @mrb by one and sets @mrb to NULL.
 *
 */
void
ncm_model_rosenbrock_clear (NcmModelRosenbrock **mrb)
{
  g_clear_object (mrb);
}

