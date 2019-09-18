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
  NcHIPertGravGauge gauge;
};

enum
{
  PROP_0,
  PROP_GAUGE,
  PROP_LEN
};

G_DEFINE_BOXED_TYPE (NcHIPertGravScalar,  nc_hipert_grav_scalar,   nc_hipert_grav_scalar_dup,   nc_hipert_grav_scalar_free);
G_DEFINE_BOXED_TYPE (NcHIPertGravVector,  nc_hipert_grav_vector,   nc_hipert_grav_vector_dup,   nc_hipert_grav_vector_free);
G_DEFINE_BOXED_TYPE (NcHIPertGravTensor,  nc_hipert_grav_tensor,   nc_hipert_grav_tensor_dup,   nc_hipert_grav_tensor_free);

G_DEFINE_BOXED_TYPE (NcHIPertGravTScalar, nc_hipert_grav_T_scalar, nc_hipert_grav_T_scalar_dup, nc_hipert_grav_T_scalar_free);
G_DEFINE_BOXED_TYPE (NcHIPertGravTVector, nc_hipert_grav_T_vector, nc_hipert_grav_T_vector_dup, nc_hipert_grav_T_vector_free);
G_DEFINE_BOXED_TYPE (NcHIPertGravTTensor, nc_hipert_grav_T_tensor, nc_hipert_grav_T_tensor_dup, nc_hipert_grav_T_tensor_free);

G_DEFINE_BOXED_TYPE (NcHIPertGravInfo,        nc_hipert_grav_info,          nc_hipert_grav_info_dup,          nc_hipert_grav_info_free);
G_DEFINE_BOXED_TYPE (NcHIPertGravTScalarInfo, nc_hipert_grav_T_scalar_info, nc_hipert_grav_T_scalar_info_dup, nc_hipert_grav_T_scalar_info_free);

G_DEFINE_TYPE_WITH_PRIVATE (NcHIPertGrav, nc_hipert_grav, G_TYPE_OBJECT);

static void
nc_hipert_grav_init (NcHIPertGrav *grav)
{
  grav->priv = nc_hipert_grav_get_instance_private (grav);

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

static guint _nc_hipert_grav_ndyn_var (NcHIPertGrav *grav);
static GArray *_nc_hipert_grav_get_deps (NcHIPertGrav *grav, guint vindex);

static void _nc_hipert_grav_set_gauge (NcHIPertGrav *grav, NcHIPertGravGauge gauge);
static NcHIPertGravGauge _nc_hipert_grav_get_gauge (NcHIPertGrav *grav);
static NcHIPertGravInfo *_nc_hipert_grav_get_G_scalar_info (NcHIPertGrav *grav);
static void _nc_hipert_grav_get_G_scalar (NcHIPertGrav *grav, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar);
static void _nc_hipert_grav_get_dy_scalar (NcHIPertGrav *grav, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar);

static void
nc_hipert_grav_class_init (NcHIPertGravClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

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

  klass->ndyn_var          = &_nc_hipert_grav_ndyn_var;
  klass->get_deps          = &_nc_hipert_grav_get_deps;

  klass->set_gauge         = &_nc_hipert_grav_set_gauge;
  klass->get_gauge         = &_nc_hipert_grav_get_gauge;
  klass->get_G_scalar_info = &_nc_hipert_grav_get_G_scalar_info;
  klass->get_G_scalar      = &_nc_hipert_grav_get_G_scalar;
  klass->get_dy_scalar     = &_nc_hipert_grav_get_dy_scalar;
}

static guint _nc_hipert_grav_ndyn_var (NcHIPertGrav *grav) { g_error ("_nc_hipert_grav_ndyn_var: not implemented by `%s'.", G_OBJECT_TYPE_NAME (grav)); return 0; }
static GArray *_nc_hipert_grav_get_deps (NcHIPertGrav *grav, guint vindex) { g_error ("_nc_hipert_grav_get_deps: not implemented by `%s'.", G_OBJECT_TYPE_NAME (grav)); return NULL; }

static void 
_nc_hipert_grav_set_gauge (NcHIPertGrav *grav, NcHIPertGravGauge gauge)
{
  grav->priv->gauge = gauge;
}

static NcHIPertGravGauge 
_nc_hipert_grav_get_gauge (NcHIPertGrav *grav)
{
  return grav->priv->gauge;
}

static NcHIPertGravInfo *_nc_hipert_grav_get_G_scalar_info (NcHIPertGrav *grav) { g_error ("_nc_hipert_grav_get_G_scalar_info: not implemented by `%s'.", G_OBJECT_TYPE_NAME (grav)); return NULL; }
static void _nc_hipert_grav_get_G_scalar (NcHIPertGrav *grav, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar) { g_error ("_nc_hipert_grav_get_G_scalar: not implemented by `%s'.", G_OBJECT_TYPE_NAME (grav)); }
static void _nc_hipert_grav_get_dy_scalar (NcHIPertGrav *grav, NcHIPertBGVar *bg_var, NcHIPertBGVarYDY *ydy, NcHIPertGravTScalar *T_scalar, NcHIPertGravScalar *G_scalar) { g_error ("_nc_hipert_grav_get_dy_scalar: not implemented by `%s'.", G_OBJECT_TYPE_NAME (grav)); }

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
 * nc_hipert_grav_T_scalar_new:
 * 
 * Creates a new #NcHIPertGravTScalar with all
 * entries set to zero.
 * 
 * Returns: (transfer full): a new #NcHIPertGravTScalar.
 */
/**
 * nc_hipert_grav_T_scalar_dup:
 * @Ts: a #NcHIPertGravTScalar
 * 
 * Duplicates @Ts.
 * 
 * Returns: (transfer full): a copy of @Ts.
 */
/**
 * nc_hipert_grav_T_scalar_free:
 * @Ts: a #NcHIPertGravTScalar
 * 
 * Frees @Ts.
 * 
 */
/**
 * nc_hipert_grav_T_scalar_add:
 * @Ts: a #NcHIPertGravTScalar
 * @Ts1: a #NcHIPertGravTScalar
 * @Ts2: a #NcHIPertGravTScalar
 * 
 * Sums @Ts1 and @Ts2 and attribute the result to @Ts.
 * 
 */
/**
 * nc_hipert_grav_T_scalar_set_zero:
 * @Ts: a #NcHIPertGravTScalar
 * 
 * Sets all @Ts entries to zero.
 * 
 */

/**
 * nc_hipert_grav_T_vector_new:
 * 
 * Creates a new #NcHIPertGravTVector with all
 * entries set to zero.
 * 
 * Returns: (transfer full): a new #NcHIPertGravTVector.
 */
/**
 * nc_hipert_grav_T_vector_dup:
 * @Tv: a #NcHIPertGravTVector
 * 
 * Duplicates @Tv.
 * 
 * Returns: (transfer full): a copy of @Tv.
 */
/**
 * nc_hipert_grav_T_vector_free:
 * @Tv: a #NcHIPertGravTVector
 * 
 * Frees @Tv.
 * 
 */
/**
 * nc_hipert_grav_T_vector_add:
 * @Tv: a #NcHIPertGravTVector
 * @Tv1: a #NcHIPertGravTVector
 * @Tv2: a #NcHIPertGravTVector
 * 
 * Sums @Tv1 and @Tv2 and attribute the result to @Tv.
 * 
 */
/**
 * nc_hipert_grav_T_vector_set_zero:
 * @Tv: a #NcHIPertGravTVector
 * 
 * Sets all @Tv entries to zero.
 * 
 */

/**
 * nc_hipert_grav_T_tensor_new:
 * 
 * Creates a new #NcHIPertGravTTensor with all
 * entries set to zero.
 * 
 * Returns: (transfer full): a new #NcHIPertGravTTensor.
 */
/**
 * nc_hipert_grav_T_tensor_dup:
 * @Tt: a #NcHIPertGravTTensor
 * 
 * Duplicates @Tt.
 * 
 * Returns: (transfer full): a copy of @Tt.
 */
/**
 * nc_hipert_grav_T_tensor_free:
 * @Tt: a #NcHIPertGravTTensor
 * 
 * Frees @Tt.
 * 
 */
/**
 * nc_hipert_grav_T_tensor_add:
 * @Tt: a #NcHIPertGravTTensor
 * @Tt1: a #NcHIPertGravTTensor
 * @Tt2: a #NcHIPertGravTTensor
 * 
 * Sums @Tt1 and @Tt2 and attribute the result to @Tt.
 * 
 */
/**
 * nc_hipert_grav_T_tensor_set_zero:
 * @Tt: a #NcHIPertGravTTensor
 * 
 * Sets all @Tt entries to zero.
 * 
 */

/**
 * nc_hipert_grav_info_new:
 * 
 * Creates a new #NcHIPertGravInfo with all
 * entries set to zero.
 * 
 * Returns: (transfer full): a new #NcHIPertGravInfo.
 */
/**
 * nc_hipert_grav_info_dup:
 * @ginfo: a #NcHIPertGravInfo
 * 
 * Duplicates @ginfo.
 * 
 * Returns: (transfer full): a copy of @ginfo.
 */
/**
 * nc_hipert_grav_info_free:
 * @ginfo: a #NcHIPertGravInfo
 * 
 * Frees @ginfo.
 * 
 */
/**
 * nc_hipert_grav_info_set_zero:
 * @ginfo: a #NcHIPertGravInfo
 * 
 * Sets all @ginfo entries to zero.
 * 
 */
/**
 * nc_hipert_grav_info_set_phi_deps:
 * @ginfo: a #NcHIPertGravInfo
 * @phi_deps: (array) (element-type gint): the $\phi$ dependencies
 * 
 * Sets $\phi$ dependencies to @phi_deps.
 */
/**
 * nc_hipert_grav_info_set_dsigma_deps:
 * @ginfo: a #NcHIPertGravInfo
 * @dsigma_deps: (array) (element-type gint): the $\dsigma$ dependencies
 * 
 * Sets $\delta\sigma$ dependencies to @dsigma_deps.
 */
/**
 * nc_hipert_grav_info_set_psi_deps:
 * @ginfo: a #NcHIPertGravInfo
 * @psi_deps: (array) (element-type gint): the $\psi$ dependencies
 * 
 * Sets $\psi$ dependencies to @psi_deps.
 */
/**
 * nc_hipert_grav_info_set_dotpsi_deps:
 * @ginfo: a #NcHIPertGravInfo
 * @dotpsi_deps: (array) (element-type gint): the $\dotpsi$ dependencies
 * 
 * Sets $\dot\psi$ dependencies to @dotpsi_deps.
 */
/**
 * nc_hipert_grav_info_get_phi_deps:
 * @ginfo: a #NcHIPertGravInfo
 * 
 * Returns: (array) (element-type gint) (transfer full): the $\phi$ dependencies
 */
/**
 * nc_hipert_grav_info_get_dsigma_deps:
 * @ginfo: a #NcHIPertGravInfo
 * 
 * Returns: (array) (element-type gint) (transfer full): the $\dsigma$ dependencies
 */
/**
 * nc_hipert_grav_info_get_psi_deps:
 * @ginfo: a #NcHIPertGravInfo
 * 
 * Returns: (array) (element-type gint) (transfer full): the $\psi$ dependencies
 */
/**
 * nc_hipert_grav_info_get_dotpsi_deps:
 * @ginfo: a #NcHIPertGravInfo
 * 
 * Returns: (array) (element-type gint) (transfer full): the $\dotpsi$ dependencies
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
 * nc_hipert_grav_ndyn_var: (virtual ndyn_var)
 * @grav: a #NcHIPertGrav
 *
 * Returns: the number of dynamical components in @grav.
 */
/**
 * nc_hipert_grav_get_deps: (virtual get_deps)
 * @grav: a #NcHIPertGrav
 * @vindex: a variable index
 *
 * Returns: (transfer full) (array) (element-type gint): the array of dependencies of variable @vindex in @grav.
 */
/**
 * nc_hipert_grav_set_gauge: (virtual set_gauge)
 * @grav: a #NcHIPertGrav
 * @gauge: a #NcHIPertGravGauge
 * 
 * Sets the gauge #NcHIPertGravGauge that @grav should use.
 * 
 */
/**
 * nc_hipert_grav_get_gauge: (virtual get_gauge)
 * @grav: a #NcHIPertGrav
 * 
 * Gets the gauge #NcHIPertGravGauge used by the gravonent @grav.
 * 
 * Returns: current gauge of @grav.
 */
/**
 * nc_hipert_grav_get_G_scalar_info: (virtual get_G_scalar_info)
 * @grav: a #NcHIPertGrav
 * 
 * Gets the determination type and dependencies for each gravitation pontential.
 * 
 * Returns: a #NcHIPertGravInfo describing the scalar sector.
 */
/**
 * nc_hipert_grav_get_G_scalar: (virtual get_G_scalar)
 * @grav: a #NcHIPertGrav
 * @bg_var: a #NcHIPertBGVar
 * @ydy: a #NcHIPertBGVarYDY
 * @T_scalar: a #NcHIPertGravTScalar
 * @G_scalar: a #NcHIPertGravScalar
 * 
 * Gets the scalar part of the gravitation potentials.
 */
/**
 * nc_hipert_grav_get_dy_scalar: (virtual get_dy_scalar)
 * @grav: a #NcHIPertGrav
 * @bg_var: a #NcHIPertBGVar
 * @ydy: a #NcHIPertBGVarYDY
 * @T_scalar: a #NcHIPertGravTScalar
 * @G_scalar: a #NcHIPertGravScalar
 * 
 * Gets the scalar part of the gravitation potentials.
 * 
 */
