/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_comp.c
 *
 *  Wed October 11 15:54:13 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_comp.c
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
 * SECTION:nc_hipert_comp
 * @title: NcHIPertComp
 * @short_description: Abstract class describing a general perturbation compoment.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/nc_hipert_comp.h"

struct _NcHIPertCompPrivate
{
  gint a;
};

enum
{
  PROP_0,
  PROP_DIM
};

G_DEFINE_BOXED_TYPE (NcHIPertCompTScalar, nc_hipert_comp_T_scalar, nc_hipert_comp_T_scalar_dup, nc_hipert_comp_T_scalar_free);
G_DEFINE_BOXED_TYPE (NcHIPertCompTVector, nc_hipert_comp_T_vector, nc_hipert_comp_T_vector_dup, nc_hipert_comp_T_vector_free);
G_DEFINE_BOXED_TYPE (NcHIPertCompTTensor, nc_hipert_comp_T_tensor, nc_hipert_comp_T_tensor_dup, nc_hipert_comp_T_tensor_free);
G_DEFINE_TYPE (NcHIPertComp, nc_hipert_comp, G_TYPE_OBJECT);

static void
nc_hipert_comp_init (NcHIPertComp *comp)
{
  comp->priv = G_TYPE_INSTANCE_GET_PRIVATE (comp, NC_TYPE_HIPERT_COMP, NcHIPertCompPrivate);
}

static void
_nc_hipert_comp_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_HIPERT_COMP (object));

  switch (prop_id)
  {
    case PROP_DIM:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_comp_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_HIPERT_COMP (object));

  switch (prop_id)
  {
    case PROP_DIM:
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_comp_dispose (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_comp_parent_class)->dispose (object);
}

static void
_nc_hipert_comp_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_comp_parent_class)->finalize (object);
}

static void
nc_hipert_comp_class_init (NcHIPertCompClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  g_type_class_add_private (klass, sizeof (NcHIPertCompPrivate));

  object_class->set_property = &_nc_hipert_comp_set_property;
  object_class->get_property = &_nc_hipert_comp_get_property;
  object_class->dispose      = &_nc_hipert_comp_dispose;
  object_class->finalize     = &_nc_hipert_comp_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DIM,
                                   g_param_spec_uint ("dim",
                                                      NULL,
                                                      "dimension",
                                                      0, 
                                                      G_MAXUINT,
                                                      2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_hipert_comp_T_scalar_new:
 * 
 * Creates a new #NcHIPertCompTScalar with all
 * entries set to zero.
 * 
 * Returns: (transfer full): a new #NcHIPertCompTScalar.
 */
/**
 * nc_hipert_comp_T_scalar_dup:
 * @Ts: a #NcHIPertCompTScalar
 * 
 * Duplicates @Ts.
 * 
 * Returns: (transfer full): a copy of @Ts.
 */
/**
 * nc_hipert_comp_T_scalar_free:
 * @Ts: a #NcHIPertCompTScalar
 * 
 * Frees @Ts.
 * 
 */
/**
 * nc_hipert_comp_T_scalar_add:
 * @Ts: a #NcHIPertCompTScalar
 * @Ts1: a #NcHIPertCompTScalar
 * @Ts2: a #NcHIPertCompTScalar
 * 
 * Sums @Ts1 and @Ts2 and attribute the result to @Ts.
 * 
 */
/**
 * nc_hipert_comp_T_scalar_set_zero:
 * @Ts: a #NcHIPertCompTScalar
 * 
 * Sets all @Ts entries to zero.
 * 
 */

/**
 * nc_hipert_comp_T_vector_new:
 * 
 * Creates a new #NcHIPertCompTVector with all
 * entries set to zero.
 * 
 * Returns: (transfer full): a new #NcHIPertCompTVector.
 */
/**
 * nc_hipert_comp_T_vector_dup:
 * @Tv: a #NcHIPertCompTVector
 * 
 * Duplicates @Tv.
 * 
 * Returns: (transfer full): a copy of @Tv.
 */
/**
 * nc_hipert_comp_T_vector_free:
 * @Tv: a #NcHIPertCompTVector
 * 
 * Frees @Tv.
 * 
 */
/**
 * nc_hipert_comp_T_vector_add:
 * @Tv: a #NcHIPertCompTVector
 * @Tv1: a #NcHIPertCompTVector
 * @Tv2: a #NcHIPertCompTVector
 * 
 * Sums @Tv1 and @Tv2 and attribute the result to @Tv.
 * 
 */
/**
 * nc_hipert_comp_T_vector_set_zero:
 * @Tv: a #NcHIPertCompTVector
 * 
 * Sets all @Tv entries to zero.
 * 
 */

/**
 * nc_hipert_comp_T_tensor_new:
 * 
 * Creates a new #NcHIPertCompTTensor with all
 * entries set to zero.
 * 
 * Returns: (transfer full): a new #NcHIPertCompTTensor.
 */
/**
 * nc_hipert_comp_T_tensor_dup:
 * @Tt: a #NcHIPertCompTTensor
 * 
 * Duplicates @Tt.
 * 
 * Returns: (transfer full): a copy of @Tt.
 */
/**
 * nc_hipert_comp_T_tensor_free:
 * @Tt: a #NcHIPertCompTTensor
 * 
 * Frees @Tt.
 * 
 */
/**
 * nc_hipert_comp_T_tensor_add:
 * @Tt: a #NcHIPertCompTTensor
 * @Tt1: a #NcHIPertCompTTensor
 * @Tt2: a #NcHIPertCompTTensor
 * 
 * Sums @Tt1 and @Tt2 and attribute the result to @Tt.
 * 
 */
/**
 * nc_hipert_comp_T_tensor_set_zero:
 * @Tt: a #NcHIPertCompTTensor
 * 
 * Sets all @Tt entries to zero.
 * 
 */
