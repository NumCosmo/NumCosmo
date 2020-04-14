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
#include "nc_enum_types.h"

struct _NcHIPertCompPrivate
{
  NcHIPertGravGauge gauge;
};

enum
{
  PROP_0,
  PROP_GAUGE,
  PROP_LEN
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcHIPertComp, nc_hipert_comp, G_TYPE_OBJECT);

static void
nc_hipert_comp_init (NcHIPertComp *comp)
{
  comp->priv = nc_hipert_comp_get_instance_private (comp);

  comp->priv->gauge = NC_HIPERT_GRAV_GAUGE_LEN;
}

static void
_nc_hipert_comp_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPertComp *comp = NC_HIPERT_COMP (object);
  g_return_if_fail (NC_IS_HIPERT_COMP (object));

  switch (prop_id)
  {
    case PROP_GAUGE:
      nc_hipert_comp_set_gauge (comp, g_value_get_enum (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_comp_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIPertComp *comp = NC_HIPERT_COMP (object);
  g_return_if_fail (NC_IS_HIPERT_COMP (object));

  switch (prop_id)
  {
    case PROP_GAUGE:
      g_value_set_enum (value, nc_hipert_comp_get_gauge (comp));
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

static guint _nc_hipert_comp_ndyn_var (NcHIPertComp *comp);
static GArray *_nc_hipert_comp_get_deps (NcHIPertComp *comp, guint vindex);

static void _nc_hipert_comp_set_gauge (NcHIPertComp *comp, NcHIPertGravGauge gauge);
static NcHIPertGravGauge _nc_hipert_comp_get_gauge (NcHIPertComp *comp);

static NcHIPertGravTScalarInfo *_nc_hipert_comp_get_T_scalar_info (NcHIPertComp *comp);

static void _nc_hipert_comp_get_T_scalar (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar);
static void _nc_hipert_comp_get_T_vector (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTVector *T_vector);
static void _nc_hipert_comp_get_T_tensor (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTTensor *T_tensor);

static void _nc_hipert_comp_get_dy_scalar (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar);

static void
nc_hipert_comp_class_init (NcHIPertCompClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_hipert_comp_set_property;
  object_class->get_property = &_nc_hipert_comp_get_property;
  object_class->dispose      = &_nc_hipert_comp_dispose;
  object_class->finalize     = &_nc_hipert_comp_finalize;

  g_object_class_install_property (object_class,
                                   PROP_GAUGE,
                                   g_param_spec_enum ("gauge",
                                                      NULL,
                                                      "gauge",
                                                      NC_TYPE_HIPERT_GRAV_GAUGE,
                                                      NC_HIPERT_GRAV_GAUGE_SYNCHRONOUS,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  klass->ndyn_var          = &_nc_hipert_comp_ndyn_var;
  klass->get_deps          = &_nc_hipert_comp_get_deps;
  klass->set_gauge         = &_nc_hipert_comp_set_gauge;
  klass->get_gauge         = &_nc_hipert_comp_get_gauge;
  klass->get_T_scalar_info = &_nc_hipert_comp_get_T_scalar_info;
  klass->get_T_scalar      = &_nc_hipert_comp_get_T_scalar;
  klass->get_T_vector      = &_nc_hipert_comp_get_T_vector;
  klass->get_T_tensor      = &_nc_hipert_comp_get_T_tensor;
  klass->get_dy_scalar     = &_nc_hipert_comp_get_dy_scalar;
}

static guint _nc_hipert_comp_ndyn_var (NcHIPertComp *comp) { g_error ("_nc_hipert_comp_ndyn_var: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); return 0; }
static GArray *_nc_hipert_comp_get_deps (NcHIPertComp *comp, guint vindex) { g_error ("_nc_hipert_comp_get_deps: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); return NULL; }

static void 
_nc_hipert_comp_set_gauge (NcHIPertComp *comp, NcHIPertGravGauge gauge)
{
  comp->priv->gauge = gauge;
}

static NcHIPertGravGauge 
_nc_hipert_comp_get_gauge (NcHIPertComp *comp)
{
  return comp->priv->gauge;
}

static NcHIPertGravTScalarInfo *_nc_hipert_comp_get_T_scalar_info (NcHIPertComp *comp)                                                                                    { g_error ("_nc_hipert_comp_get_T_scalar_info: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); return NULL; }

static void _nc_hipert_comp_get_T_scalar (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar)                                { g_error ("_nc_hipert_comp_get_T_scalar: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); }
static void _nc_hipert_comp_get_T_vector (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTVector *T_vector)                                { g_error ("_nc_hipert_comp_get_T_vector: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); }
static void _nc_hipert_comp_get_T_tensor (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTTensor *T_tensor)                                { g_error ("_nc_hipert_comp_get_T_tensor: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); }

static void _nc_hipert_comp_get_dy_scalar (NcHIPertComp *comp, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar) { g_error ("_nc_hipert_comp_get_dy_scalar: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); }

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
 * nc_hipert_comp_get_deps: (virtual get_deps)
 * @comp: a #NcHIPertComp
 * @vindex: a variable index
 *
 * Returns: (transfer full) (array) (element-type gint): the array of dependencies of variable @vindex in @comp.
 */
GArray *
nc_hipert_comp_get_deps (NcHIPertComp *comp, guint vindex)
{
  return NC_HIPERT_COMP_GET_CLASS (comp)->get_deps (comp, vindex);
}

/**
 * nc_hipert_comp_set_gauge: (virtual set_gauge)
 * @comp: a #NcHIPertComp
 * @gauge: a #NcHIPertGravGauge
 * 
 * Sets the gauge #NcHIPertGravGauge that the component @comp
 * should use.
 * 
 */
void 
nc_hipert_comp_set_gauge (NcHIPertComp *comp, NcHIPertGravGauge gauge)
{
  NC_HIPERT_COMP_GET_CLASS (comp)->set_gauge (comp, gauge);
}

/**
 * nc_hipert_comp_get_gauge: (virtual get_gauge)
 * @comp: a #NcHIPertComp
 * 
 * Gets the gauge #NcHIPertGravGauge used by the component @comp.
 * 
 * Returns: current gauge of @comp.
 */
NcHIPertGravGauge 
nc_hipert_comp_get_gauge (NcHIPertComp *comp)
{
  return NC_HIPERT_COMP_GET_CLASS (comp)->get_gauge (comp);
}

/**
 * nc_hipert_comp_get_Tscalar_coupling: (virtual get_T_scalar_info)
 * @comp: a #NcHIPertComp
 * 
 * Provides all information about the scalar energy momentum tensor.
 * 
 * Returns: a #NcHIPertGravTScalarInfo.
 */
NcHIPertGravTScalarInfo *
nc_hipert_comp_get_T_scalar_info (NcHIPertComp *comp)
{
  return NC_HIPERT_COMP_GET_CLASS (comp)->get_T_scalar_info (comp);
}

/**
 * nc_hipert_comp_get_T_scalar: (virtual get_T_scalar)
 * @comp: a #NcHIPertComp
 * @bg_var: a #NcHIPertBGVar
 * @ydy: a #NcHIPertBGVarYDY
 * @T_scalar: a #NcHIPertGravTScalar
 * 
 * Calculates the current value of the energy momentum tensor 
 * and stores it in @T_scalar.
 * 
 */
/**
 * nc_hipert_comp_get_dy_scalar: (virtual get_dy_scalar)
 * @comp: a #NcHIPertComp
 * @bg_var: a #NcHIPertBGVar
 * @ydy: a #NcHIPertBGVarYDY
 * @T_scalar: a #NcHIPertGravTScalar
 * @G_scalar: a #NcHIPertGravScalar
 * 
 * Calculates the time derivative of the dynamical variables.
 * 
 */
