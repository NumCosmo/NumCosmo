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
  NcHIPertCompGauge gauge;
};

enum
{
  PROP_0,
  PROP_GAUGE,
  PROP_LEN
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
static GArray *_nc_hipert_comp_cgraph (NcHIPertComp *comp);

static void _nc_hipert_comp_set_gauge (NcHIPertComp *comp, NcHIPertCompGauge gauge);
static NcHIPertCompGauge _nc_hipert_comp_get_gauge (NcHIPertComp *comp);

static void _nc_hipert_comp_get_Tscalar_coupling (NcHIPertComp *comp, GArray **drho, GArray **v, GArray **dp, GArray **Pi);

static void _nc_hipert_comp_dy (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcmVector *dy);
static void _nc_hipert_comp_J (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcmMatrix *J);
static void _nc_hipert_comp_dy_J (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcmVector *dy, NcmMatrix *J);
static void _nc_hipert_comp_Tscalar (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcHIPertCompTScalar *TScalar);
static void _nc_hipert_comp_Tvector (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcHIPertCompTScalar *TVector);
static void _nc_hipert_comp_Ttensor (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcHIPertCompTScalar *TTensor);

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
                                   PROP_GAUGE,
                                   g_param_spec_enum ("gauge",
                                                      NULL,
                                                      "gauge",
                                                      NC_TYPE_HIPERT_GRAV_GAUGE,
                                                      NC_HIPERT_GRAV_GAUGE_SYNCHRONOUS,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  klass->ndyn_var             = &_nc_hipert_comp_ndyn_var;
  klass->cgraph               = &_nc_hipert_comp_cgraph;
  klass->set_gauge            = &_nc_hipert_comp_set_gauge;
  klass->get_gauge            = &_nc_hipert_comp_get_gauge;
  klass->get_Tscalar_coupling = &_nc_hipert_comp_get_Tscalar_coupling;
  klass->dy                   = &_nc_hipert_comp_dy;
  klass->J                    = &_nc_hipert_comp_J;
  klass->dy_J                 = &_nc_hipert_comp_dy_J;
  klass->Tscalar              = &_nc_hipert_comp_Tscalar;
  klass->Tvector              = &_nc_hipert_comp_Tvector;
  klass->Ttensor              = &_nc_hipert_comp_Ttensor;
}

static guint _nc_hipert_comp_ndyn_var (NcHIPertComp *comp) { g_error ("_nc_hipert_comp_ndyn_var: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); return 0; }
static GArray *_nc_hipert_comp_cgraph (NcHIPertComp *comp) { g_error ("_nc_hipert_comp_cgraph: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); return NULL; }

static void 
_nc_hipert_comp_set_gauge (NcHIPertComp *comp, NcHIPertCompGauge gauge)
{
  comp->priv->gauge = gauge;
}

static NcHIPertCompGauge 
_nc_hipert_comp_get_gauge (NcHIPertComp *comp)
{
  return comp->priv->gauge;
}

static void _nc_hipert_comp_get_Tscalar_coupling (NcHIPertComp *comp, GArray **drho, GArray **v, GArray **dp, GArray **Pi)        { g_error ("_nc_hipert_comp_get_Tscalar_coupling: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); }

static void _nc_hipert_comp_dy (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcmVector *dy)                     { g_error ("_nc_hipert_comp_dy: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); }
static void _nc_hipert_comp_J (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcmMatrix *J)                       { g_error ("_nc_hipert_comp_J: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); }
static void _nc_hipert_comp_dy_J (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcmVector *dy, NcmMatrix *J)     { g_error ("_nc_hipert_comp_dy_J: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); }
static void _nc_hipert_comp_Tscalar (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcHIPertCompTScalar *TScalar) { g_error ("_nc_hipert_comp_Tscalar: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); }
static void _nc_hipert_comp_Tvector (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcHIPertCompTScalar *TVector) { g_error ("_nc_hipert_comp_Tvector: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); }
static void _nc_hipert_comp_Ttensor (NcHIPertComp *comp, NcHIPertBGVar *bg_var, const NcmVector *y, NcHIPertCompTScalar *TTensor) { g_error ("_nc_hipert_comp_Ttensor: not implemented by `%s'.", G_OBJECT_TYPE_NAME (comp)); }

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
 * nc_hipert_comp_set_gauge: (virtual set_gauge)
 * @comp: a #NcHIPertComp
 * @gauge: a #NcHIPertCompGauge
 * 
 * Sets the gauge #NcHIPertCompGauge that the component @comp
 * should use.
 * 
 */
void 
nc_hipert_comp_set_gauge (NcHIPertComp *comp, NcHIPertCompGauge gauge)
{
  NC_HIPERT_COMP_GET_CLASS (comp)->set_gauge (comp, gauge);
}

/**
 * nc_hipert_comp_get_gauge: (virtual get_gauge)
 * @comp: a #NcHIPertComp
 * 
 * Gets the gauge #NcHIPertCompGauge used by the component @comp.
 * 
 * Returns: current gauge of @comp.
 */
NcHIPertCompGauge 
nc_hipert_comp_get_gauge (NcHIPertComp *comp)
{
  return NC_HIPERT_COMP_GET_CLASS (comp)->get_gauge (comp);
}

/**
 * nc_hipert_comp_get_Tscalar_coupling: (virtual get_Tscalar_coupling)
 * @comp: a #NcHIPertComp
 * @drho: (out) (array) (element-type gint) (transfer full): The list of components used in $\delta\rho$
 * @v: (out) (array) (element-type gint) (transfer full): The list of components used in $\mathcal{V}$
 * @dp: (out) (array) (element-type gint) (transfer full): The list of components used in $\delta{}p$
 * @Pi: (out) (array) (element-type gint) (transfer full): The list of components used in $\Pi$
 * 
 * Provides the lists of components present in each components of the 
 * @comp energy momentum tensor.
 * 
 */
void 
nc_hipert_comp_get_Tscalar_coupling (NcHIPertComp *comp, GArray **drho, GArray **v, GArray **dp, GArray **Pi)
{
  NC_HIPERT_COMP_GET_CLASS (comp)->get_Tscalar_coupling (comp, drho, v, dp, Pi);
}

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
