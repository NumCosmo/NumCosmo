/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_grav.c
 *
 *  Thu October 12 14:32:22 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_grav.c
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_hipert_grav
 * @title: NcHIPertGrav
 * @short_description: Abstract class describing a general first order gravitation theory.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/nc_hipert_grav.h"
#include "nc_enum_types.h"

struct _NcHIPertGravPrivate
{
  NcHIPertCompGauge gauge;
};

enum
{
  PROP_0,
  PROP_GAUGE,
  PROP_LEN
};

G_DEFINE_BOXED_TYPE (NcHIPertGravScalar, nc_hipert_grav_scalar, nc_hipert_grav_scalar_dup, nc_hipert_grav_scalar_free);
G_DEFINE_BOXED_TYPE (NcHIPertGravVector, nc_hipert_grav_vector, nc_hipert_grav_vector_dup, nc_hipert_grav_vector_free);
G_DEFINE_BOXED_TYPE (NcHIPertGravTensor, nc_hipert_grav_tensor, nc_hipert_grav_tensor_dup, nc_hipert_grav_tensor_free);
G_DEFINE_TYPE (NcHIPertGrav, nc_hipert_grav, G_TYPE_OBJECT);

static void
nc_hipert_grav_init (NcHIPertGrav *grav)
{
  grav->priv = G_TYPE_INSTANCE_GET_PRIVATE (grav, NC_TYPE_HIPERT_GRAV, NcHIPertGravPrivate);

  grav->priv->gauge = NC_HIPERT_GRAV_GAUGE_LEN;
}

static void
_nc_hipert_grav_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPertGrav *grav = NC_HIPERT_GRAV (object);
  g_return_if_fail (NC_IS_HIPERT_GRAV (object));

  switch (prop_id)
  {
    case PROP_GAUGE:
      nc_hipert_grav_set_gauge (grav, g_value_get_enum (value));
      break;
   default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_grav_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIPertGrav *grav = NC_HIPERT_GRAV (object);
  g_return_if_fail (NC_IS_HIPERT_GRAV (object));

  switch (prop_id)
  {
    case PROP_GAUGE:
      g_value_set_enum (value, nc_hipert_grav_get_gauge (grav));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_grav_dispose (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_grav_parent_class)->dispose (object);
}

static void
_nc_hipert_grav_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_grav_parent_class)->finalize (object);
}

static void _nc_hipert_grav_set_gauge (NcHIPertGrav *grav, NcHIPertCompGauge gauge);
static NcHIPertCompGauge _nc_hipert_grav_get_gauge (NcHIPertGrav *grav);

static void
nc_hipert_grav_class_init (NcHIPertGravClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  g_type_class_add_private (klass, sizeof (NcHIPertGravPrivate));

  object_class->set_property = &_nc_hipert_grav_set_property;
  object_class->get_property = &_nc_hipert_grav_get_property;
  object_class->dispose      = &_nc_hipert_grav_dispose;
  object_class->finalize     = &_nc_hipert_grav_finalize;

  g_object_class_install_property (object_class,
                                   PROP_GAUGE,
                                   g_param_spec_enum ("gauge",
                                                      NULL,
                                                      "gauge",
                                                      NC_TYPE_HIPERT_GRAV_GAUGE,
                                                      NC_HIPERT_GRAV_GAUGE_SYNCHRONOUS,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->set_gauge = &_nc_hipert_grav_set_gauge;
  klass->get_gauge = &_nc_hipert_grav_get_gauge;
}

static void 
_nc_hipert_grav_set_gauge (NcHIPertGrav *grav, NcHIPertCompGauge gauge)
{
  grav->priv->gauge = gauge;
}

static NcHIPertCompGauge 
_nc_hipert_grav_get_gauge (NcHIPertGrav *grav)
{
  return grav->priv->gauge;
}

/**
 * nc_hipert_grav_scalar_new:
 * 
 * Creates a new #NcHIPertGravScalar with all
 * entries set to zero.
 * 
 * Returns: (transfer full): a new #NcHIPertGravScalar.
 */
/**
 * nc_hipert_grav_scalar_dup:
 * @gs: a #NcHIPertGravScalar
 * 
 * Duplicates @gs.
 * 
 * Returns: (transfer full): a copy of @gs.
 */
/**
 * nc_hipert_grav_scalar_free:
 * @gs: a #NcHIPertGravScalar
 * 
 * Frees @gs.
 * 
 */
/**
 * nc_hipert_grav_scalar_set_zero:
 * @gs: a #NcHIPertGravScalar
 * 
 * Sets all @gs entries to zero.
 * 
 */

/**
 * nc_hipert_grav_vector_new:
 * 
 * Creates a new #NcHIPertGravVector with all
 * entries set to zero.
 * 
 * Returns: (transfer full): a new #NcHIPertGravVector.
 */
/**
 * nc_hipert_grav_vector_dup:
 * @gv: a #NcHIPertGravVector
 * 
 * Duplicates @gv.
 * 
 * Returns: (transfer full): a copy of @gv.
 */
/**
 * nc_hipert_grav_vector_free:
 * @gv: a #NcHIPertGravVector
 * 
 * Frees @gv.
 * 
 */
/**
 * nc_hipert_grav_vector_set_zero:
 * @gv: a #NcHIPertGravVector
 * 
 * Sets all @gv entries to zero.
 * 
 */

/**
 * nc_hipert_grav_tensor_new:
 * 
 * Creates a new #NcHIPertGravTensor with all
 * entries set to zero.
 * 
 * Returns: (transfer full): a new #NcHIPertGravTensor.
 */
/**
 * nc_hipert_grav_tensor_dup:
 * @gt: a #NcHIPertGravTensor
 * 
 * Duplicates @gt.
 * 
 * Returns: (transfer full): a copy of @gt.
 */
/**
 * nc_hipert_grav_tensor_free:
 * @gt: a #NcHIPertGravTensor
 * 
 * Frees @gt.
 * 
 */
/**
 * nc_hipert_grav_tensor_set_zero:
 * @gt: a #NcHIPertGravTensor
 * 
 * Sets all @gt entries to zero.
 * 
 */

/**
 * nc_hipert_grav_ref:
 * @grav: a #NcHIPertGrav
 *
 * Increases the reference count of @grav.
 *
 * Returns: (transfer full): @grav.
 */
NcHIPertGrav *
nc_hipert_grav_ref (NcHIPertGrav *grav)
{
  return g_object_ref (grav);
}

/**
 * nc_hipert_grav_free:
 * @grav: a #NcHIPertGrav
 *
 * Decreases the reference count of @grav.
 *
 */
void 
nc_hipert_grav_free (NcHIPertGrav *grav)
{
  g_object_unref (grav);
}

/**
 * nc_hipert_grav_clear:
 * @grav: a #NcHIPertGrav
 *
 * Decreases the reference count of *@grav and sets the pointer *@grav to NULL.
 *
 */
void 
nc_hipert_grav_clear (NcHIPertGrav **grav)
{
  g_clear_object (grav);
}

/**
 * nc_hipert_grav_get_id:
 * @grav: a #NcHIPertGrav
 *
 * Returns: the #NcHIPertBGVar id tied to this gravitation object.
 */

/**
 * nc_hipert_grav_set_gauge: (virtual set_gauge)
 * @grav: a #NcHIPertGrav
 * @gauge: a #NcHIPertCompGauge
 * 
 * Sets the gauge #NcHIPertCompGauge that @grav should use.
 * 
 */
void 
nc_hipert_grav_set_gauge (NcHIPertGrav *grav, NcHIPertCompGauge gauge)
{
  NC_HIPERT_GRAV_GET_CLASS (grav)->set_gauge (grav, gauge);
}

/**
 * nc_hipert_grav_get_gauge: (virtual get_gauge)
 * @grav: a #NcHIPertGrav
 * 
 * Gets the gauge #NcHIPertCompGauge used by the gravonent @grav.
 * 
 * Returns: current gauge of @grav.
 */
NcHIPertCompGauge 
nc_hipert_grav_get_gauge (NcHIPertGrav *grav)
{
  return NC_HIPERT_GRAV_GET_CLASS (grav)->get_gauge (grav);
}

