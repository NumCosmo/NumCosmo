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

G_DEFINE_BOXED_TYPE (NcHIPertCompTScalar,  nc_hipert_comp_T_scalar, nc_hipert_comp_T_scalar_dup, nc_hipert_comp_T_scalar_free);
G_DEFINE_BOXED_TYPE (NcHIPertCompTVector,  nc_hipert_comp_T_vector, nc_hipert_comp_T_vector_dup, nc_hipert_comp_T_vector_free);
G_DEFINE_BOXED_TYPE (NcHIPertCompTTensor,  nc_hipert_comp_T_tensor, nc_hipert_comp_T_tensor_dup, nc_hipert_comp_T_tensor_free);
G_DEFINE_BOXED_TYPE (NcHIPertCompCoupling, nc_hipert_comp_coupling, nc_hipert_comp_coupling_dup, nc_hipert_comp_coupling_free);
G_DEFINE_ABSTRACT_TYPE (NcHIPertComp, nc_hipert_comp, G_TYPE_OBJECT);

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

  klass->bg_var_id = -1;
}

G_LOCK_DEFINE_STATIC (last_bg_var_id);

/**
 * nc_hipert_comp_register_bg_var_id: (skip)
 * @comp_class: a #NcHIPertCompClass
 * @cstruct_size: component struct size
 * @ns: component namespace
 * @desc: short description
 * @long_desc: long description
 *
 * FIXME
 *
 */
void
nc_hipert_comp_register_bg_var_id (NcHIPertCompClass *comp_class, guint cstruct_size, const gchar *ns, const gchar *desc, const gchar *long_desc)
{
  static NcHIPertBGVarID last_bg_var_id = 0;
  NcHIPertBGVarClass *bg_var_class     = g_type_class_ref (NC_TYPE_HIPERT_BG_VAR);
  NcHIPertBGVarDesc *bg_var_desc       = NULL;
  guint id;

  G_LOCK (last_bg_var_id);

  comp_class->bg_var_id = last_bg_var_id;
  id                    = last_bg_var_id;

  last_bg_var_id++;
  
  g_array_set_size (bg_var_class->bg_var_desc_array, last_bg_var_id);

  bg_var_desc       = &g_array_index (bg_var_class->bg_var_desc_array, NcHIPertBGVarDesc, id);
  bg_var_desc->init = TRUE;

  if (ns == NULL)
    g_error ("Cannot register background variables without a namespace.");
  if (desc == NULL)
    g_error ("Cannot register background variables without a description.");

  bg_var_desc->ns           = g_strdup (ns);
  bg_var_desc->desc         = g_strdup (desc);
  bg_var_desc->cstruct_size = cstruct_size;

  if (long_desc != NULL)
    bg_var_desc->long_desc = g_strdup (long_desc);
  else
    bg_var_desc->long_desc = NULL;

  if (g_hash_table_lookup (bg_var_class->ns_table, ns) != NULL)
    g_error ("Background variable namespace <%s> already registered.", ns);

  g_hash_table_insert (bg_var_class->ns_table, bg_var_desc->ns, GINT_TO_POINTER (comp_class->bg_var_id));

  G_UNLOCK (last_bg_var_id);
  g_type_class_unref (bg_var_class);
  
  return;
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

/**
 * nc_hipert_comp_coupling_new:
 * 
 * Creates a new #NcHIPertCompCoupling with all
 * entries set to zero.
 * 
 * Returns: (transfer full): a new #NcHIPertCompCoupling.
 */
/**
 * nc_hipert_comp_coupling_dup:
 * @c: a #NcHIPertCompCoupling
 * 
 * Duplicates @c.
 * 
 * Returns: (transfer full): a copy of @c.
 */
/**
 * nc_hipert_comp_coupling_free:
 * @c: a #NcHIPertCompCoupling
 * 
 * Frees @c.
 * 
 */

/**
 * nc_hipert_comp_ref:
 * @comp: a #NcHIPertComp
 *
 * Increases the reference count of @comp.
 *
 * Returns: (transfer full): @comp.
 */
NcHIPertComp *
nc_hipert_comp_ref (NcHIPertComp *comp)
{
  return g_object_ref (comp);
}

/**
 * nc_hipert_comp_free:
 * @comp: a #NcHIPertComp
 *
 * Decreases the reference count of @comp.
 *
 */
void 
nc_hipert_comp_free (NcHIPertComp *comp)
{
  g_object_unref (comp);
}

/**
 * nc_hipert_comp_clear:
 * @comp: a #NcHIPertComp
 *
 * Decreases the reference count of *@comp and sets the pointer *@comp to NULL.
 *
 */
void 
nc_hipert_comp_clear (NcHIPertComp **comp)
{
  g_clear_object (comp);
}

/**
 * nc_hipert_comp_get_id:
 * @comp: a #NcHIPertComp
 *
 * Returns: the #NcHIPertBGVar id tied to this component.
 */

/**
 * nc_hipert_comp_ndyn_var: (virtual ndyn_var)
 * @comp: a #NcHIPertComp
 *
 * Returns: the number of dynamical components in @comp.
 */
guint 
nc_hipert_comp_ndyn_var (NcHIPertComp *comp)
{
  return NC_HIPERT_COMP_GET_CLASS (comp)->ndyn_var (comp);
}

/**
 * nc_hipert_comp_coupling_graph: (virtual cgraph)
 * @comp: a #NcHIPertComp
 *
 * Returns: (transfer full) (array) (element-type NcHIPertCompCoupling): the number of dynamical components in @comp.
 */
GArray *
nc_hipert_comp_coupling_graph (NcHIPertComp *comp)
{
  return NC_HIPERT_COMP_GET_CLASS (comp)->cgraph (comp);
}
