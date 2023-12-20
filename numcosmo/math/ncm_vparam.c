/***************************************************************************
 *            ncm_vparam.c
 *
 *  Thu May 10 15:50:20 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
 * SECTION:ncm_vparam
 * @title: NcmVParam
 * @short_description: Properties of a vector-like parameter.
 *
 * This object comprises the necessary properties to define a vector parameter.
 * It is used by #NcmModel to store the description of the vector model parameters.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_vparam.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

enum
{
  PROP_0,
  PROP_DEFAULT_SPARAM,
  PROP_LEN,
};

struct _NcmVParam
{
  /*< private >*/
  GObject parent_instance;
  guint len;
  NcmSParam *default_sparam;
  GPtrArray *sparam;
};


G_DEFINE_TYPE (NcmVParam, ncm_vparam, G_TYPE_OBJECT)

static void
ncm_vparam_init (NcmVParam *vp)
{
  vp->len            = 0;
  vp->default_sparam = NULL;
  vp->sparam         = g_ptr_array_new_with_free_func ((GDestroyNotify) & ncm_sparam_free);
}

static void
_ncm_vparam_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmVParam *vparam = NCM_VPARAM (object);

  g_return_if_fail (NCM_IS_VPARAM (object));

  switch (prop_id)
  {
    case PROP_DEFAULT_SPARAM:
      vparam->default_sparam = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_vparam_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmVParam *vparam = NCM_VPARAM (object);

  g_return_if_fail (NCM_IS_VPARAM (vparam));

  switch (prop_id)
  {
    case PROP_DEFAULT_SPARAM:
      g_value_set_object (value, vparam->default_sparam);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_vparam_dispose (GObject *object)
{
  NcmVParam *vp = NCM_VPARAM (object);

  vp->len = 0;
  ncm_sparam_clear (&vp->default_sparam);

  if (vp->sparam != NULL)
  {
    g_ptr_array_unref (vp->sparam);
    vp->sparam = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_vparam_parent_class)->dispose (object);
}

static void
_ncm_vparam_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_vparam_parent_class)->finalize (object);
}

static void
ncm_vparam_class_init (NcmVParamClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = _ncm_vparam_set_property;
  object_class->get_property = _ncm_vparam_get_property;
  object_class->dispose      = _ncm_vparam_dispose;
  object_class->finalize     = _ncm_vparam_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DEFAULT_SPARAM,
                                   g_param_spec_object  ("default-sparam",
                                                         NULL,
                                                         "Default sparam for the vector components",
                                                         NCM_TYPE_SPARAM,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_vparam_new:
 * @len: vector length.
 * @default_param: a #NcmSParam.
 *
 * This function allocates memory for a new #NcmVParam object and sets its properties to the values from
 * the input arguments. @len provides the number of components.
 *
 * Returns: A new #NcmVParam.
 */
NcmVParam *
ncm_vparam_new (guint len, NcmSParam *default_param)
{
  NcmVParam *vp = g_object_new (NCM_TYPE_VPARAM,
                                "default-sparam", default_param,
                                NULL);

  ncm_vparam_set_len (vp, len);

  return vp;
}

/**
 * ncm_vparam_full_new:
 * @len: vector length.
 * @name: #NcmSParam:name.
 * @symbol: #NcmSParam:symbol.
 * @lower_bound: value of #NcmSParam:lower-bound.
 * @upper_bound: value of #NcmSParam:upper-bound.
 * @scale: value of #NcmSParam:scale.
 * @abstol: value of #NcmSParam:absolute-tolerance.
 * @default_val: value of #NcmSParam:default-value.
 * @ftype: a #NcmParamType.
 *
 * This function allocates memory for a new #NcmVParam object and sets its properties to the values from
 * the input arguments.
 *
 * The @name parameter is restricted to the interval [@lower_bound, @upper_bound].
 * @scale is an initial step for the statistical algorithms.
 * @abstol is the absolute error tolerance of the parameter.
 * @ftype indicates if the parameter will be fitted or not.
 *
 * Returns: A new #NcmVParam.
 */
NcmVParam *
ncm_vparam_full_new (guint len, const gchar *name, const gchar *symbol, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_val, NcmParamType ftype)
{
  NcmSParam *default_param = ncm_sparam_new (name, symbol, lower_bound, upper_bound, scale, abstol, default_val, ftype);
  NcmVParam *vparam        = ncm_vparam_new (len, default_param);

  ncm_sparam_free (default_param);

  return vparam;
}

/**
 * ncm_vparam_ref:
 * @vparam: a #NcmVParam.
 *
 * Increases the reference count of @vparam by one.
 *
 * Returns: (transfer full): @vparam
 */
NcmVParam *
ncm_vparam_ref (NcmVParam *vparam)
{
  return g_object_ref (vparam);
}

/**
 * ncm_vparam_copy:
 * @vparam: a #NcmVParam.
 *
 * Duplicates the #NcmVParam object setting the same values of the original propertities.
 *
 * Returns: (transfer full): A new #NcmVParam.
 */
NcmVParam *
ncm_vparam_copy (NcmVParam *vparam)
{
  NcmSParam *default_param = ncm_sparam_copy (vparam->default_sparam);
  NcmVParam *vparam_new    = ncm_vparam_new (0, default_param);
  guint i;

  ncm_sparam_free (default_param);

  g_ptr_array_set_size (vparam_new->sparam, vparam->len);
  vparam_new->len = vparam->len;

  for (i = 0; i < vparam->len; i++)
  {
    NcmSParam *sp = ncm_sparam_copy (ncm_vparam_peek_sparam (vparam, i));

    g_ptr_array_index (vparam_new->sparam, i) = sp;
  }

  return vparam_new;
}

/**
 * ncm_vparam_free:
 * @vparam: a #NcmVParam.
 *
 * Atomically decrements the reference count of @vparam by one. If the reference count drops to 0,
 * all memory allocated by @vparam is released.
 *
 */
void
ncm_vparam_free (NcmVParam *vparam)
{
  g_object_unref (vparam);
}

/**
 * ncm_vparam_clear:
 * @vparam: a #NcmVParam.
 *
 * Atomically decrements the reference count of @vparam by one. If the reference count drops to 0,
 * all memory allocated by @vparam is released.
 *
 */
void
ncm_vparam_clear (NcmVParam **vparam)
{
  g_clear_object (vparam);
}

/**
 * ncm_vparam_set_len:
 * @vparam: a #NcmVParam.
 * @len: lenght of the #NcmVParam.
 *
 * Sets the length of @vparam to @len.
 *
 */
void
ncm_vparam_set_len (NcmVParam *vparam, guint len)
{
  guint i;

  g_ptr_array_set_size (vparam->sparam, len);

  for (i = vparam->len; i < len; i++)
  {
    g_ptr_array_index (vparam->sparam, i) = ncm_sparam_copy (vparam->default_sparam);
    ncm_sparam_take_name (g_ptr_array_index (vparam->sparam, i),
                          g_strdup_printf ("%s_%u",
                                           ncm_sparam_name (vparam->default_sparam), i)
                         );
    ncm_sparam_take_symbol (g_ptr_array_index (vparam->sparam, i),
                            g_strdup_printf ("{%s}_%u",
                                             ncm_sparam_symbol (vparam->default_sparam), i)
                           );
  }

  vparam->len = len;
}

/**
 * ncm_vparam_get_len:
 * @vparam: a #NcmVParam.
 *
 * Returns: The length of @vparam.
 */
guint
ncm_vparam_get_len (NcmVParam *vparam)
{
  return vparam->len;
}

/**
 * ncm_vparam_set_sparam:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 * @spn: a #NcmSParam.
 *
 * Sets the #NcmSParam associated with the @n-th component of #NcmVParam.
 *
 */
void
ncm_vparam_set_sparam (NcmVParam *vparam, guint n, NcmSParam *spn)
{
  g_assert (n < vparam->len);
  ncm_sparam_free (g_ptr_array_index (vparam->sparam, n));
  g_ptr_array_index (vparam->sparam, n) = spn;
  ncm_sparam_ref (spn);
}

/**
 * ncm_vparam_set_sparam_full:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 * @name: #NcmSParam:name.
 * @symbol: #NcmSParam:symbol.
 * @lower_bound: value of #NcmSParam:lower-bound.
 * @upper_bound: value of #NcmSParam:upper-bound.
 * @scale: value of #NcmSParam:scale.
 * @abstol: value of #NcmSParam:absolute-tolerance.
 * @default_val: value of #NcmSParam:default-value.
 * @ftype: a #NcmParamType.
 *
 * This function sets the properties of the @n-th @vparam component.
 *
 */
void
ncm_vparam_set_sparam_full (NcmVParam *vparam, guint n, gchar *name, gchar *symbol, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_val, NcmParamType ftype)
{
  NcmSParam *spn = ncm_sparam_new (name, symbol, lower_bound, upper_bound, scale, abstol, default_val, ftype);

  g_assert (n < vparam->len);
  ncm_sparam_free (g_ptr_array_index (vparam->sparam, n));
  g_ptr_array_index (vparam->sparam, n) = spn;
}

/**
 * ncm_vparam_peek_sparam:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 *
 * This function does not increment the reference count of #NcmSParam.
 *
 * Returns: (transfer none): A #NcmSParam, which is the @n-th component of @vparam.
 */
NcmSParam *
ncm_vparam_peek_sparam (const NcmVParam *vparam, guint n)
{
  return g_ptr_array_index (vparam->sparam, n);
}

/**
 * ncm_vparam_get_sparam:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 *
 * This function returns the @n-th component of @vparam increasing its reference count.
 *
 * Returns: (transfer full): A #NcmSParam.
 */
NcmSParam *
ncm_vparam_get_sparam (NcmVParam *vparam, guint n)
{
  NcmSParam *sp = ncm_vparam_peek_sparam (vparam, n);

  ncm_sparam_ref (sp);

  return sp;
}

/**
 * ncm_vparam_set_lower_bound:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 * @lb: value of #NcmSParam:lower-bound.
 *
 * Sets the value @lb to the #NcmSParam:lower-bound property of the @n-th component of @vparam.
 */
void
ncm_vparam_set_lower_bound (NcmVParam *vparam, guint n, const gdouble lb)
{
  NcmSParam *sp = ncm_vparam_peek_sparam (vparam, n);

  ncm_sparam_set_lower_bound (sp, lb);
}

/**
 * ncm_vparam_set_upper_bound:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 * @ub: value of #NcmSParam:upper-bound.
 *
 * Sets the value @ub to the #NcmSParam:upper-bound property of the @n-th component of @vparam.
 */
void
ncm_vparam_set_upper_bound (NcmVParam *vparam, guint n, const gdouble ub)
{
  NcmSParam *sp = ncm_vparam_peek_sparam (vparam, n);

  ncm_sparam_set_upper_bound (sp, ub);
}

/**
 * ncm_vparam_set_scale:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 * @scale: value of #NcmSParam:scale.
 *
 * Sets the value @scale to the #NcmSParam:scale property of the @n-th component of @vparam.
 */
void
ncm_vparam_set_scale (NcmVParam *vparam, guint n, const gdouble scale)
{
  NcmSParam *sp = ncm_vparam_peek_sparam (vparam, n);

  ncm_sparam_set_scale (sp, scale);
}

/**
 * ncm_vparam_set_absolute_tolerance:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 * @abstol: value of #NcmSParam:absolute-tolerance.
 *
 * Sets the value @abstol to the #NcmSParam:absolute-tolerance property of the @n-th component of @vparam.
 */
void
ncm_vparam_set_absolute_tolerance (NcmVParam *vparam, guint n, const gdouble abstol)
{
  NcmSParam *sp = ncm_vparam_peek_sparam (vparam, n);

  ncm_sparam_set_absolute_tolerance (sp, abstol);
}

/**
 * ncm_vparam_set_default_value:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 * @default_val: value of #NcmSParam:default-value.
 *
 * Sets the value @default_val to the #NcmSParam:default-value property of the @n-th component of @vparam.
 */
void
ncm_vparam_set_default_value (NcmVParam *vparam, guint n, const gdouble default_val)
{
  NcmSParam *sp = ncm_vparam_peek_sparam (vparam, n);

  ncm_sparam_set_default_value (sp, default_val);
}

/**
 * ncm_vparam_set_fit_type:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 * @ftype: a #NcmParamType.
 *
 * Sets @ftype to the #NcmSParam:fit-type property of the @n-th component of @vparam.
 */
void
ncm_vparam_set_fit_type (NcmVParam *vparam, guint n, const NcmParamType ftype)
{
  NcmSParam *sp = ncm_vparam_peek_sparam (vparam, n);

  ncm_sparam_set_default_value (sp, ftype);
}

/**
 * ncm_vparam_name:
 * @vparam: a #NcmVParam.
 *
 * Gets the @vparam base name.
 *
 * Returns: (transfer none): @vparam base name
 */
const gchar *
ncm_vparam_name (const NcmVParam *vparam)
{
  return ncm_sparam_name (vparam->default_sparam);
}

/**
 * ncm_vparam_symbol:
 * @vparam: a #NcmVParam.
 *
 * Gets the @vparam base symbol.
 *
 * Returns: (transfer none): @vparam base symbol
 */
const gchar *
ncm_vparam_symbol (const NcmVParam *vparam)
{
  return ncm_sparam_symbol (vparam->default_sparam);
}

/**
 * ncm_vparam_get_lower_bound:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 *
 * Returns: The value of #NcmSParam:lower-bound property of the @n-th component of @vparam.
 */
gdouble
ncm_vparam_get_lower_bound (const NcmVParam *vparam, guint n)
{
  NcmSParam *sp = ncm_vparam_peek_sparam (vparam, n);

  return ncm_sparam_get_lower_bound (sp);
}

/**
 * ncm_vparam_get_upper_bound:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 *
 * Returns: The value of #NcmSParam:upper-bound property of the @n-th component of @vparam.
 */
gdouble
ncm_vparam_get_upper_bound (const NcmVParam *vparam, guint n)
{
  NcmSParam *sp = ncm_vparam_peek_sparam (vparam, n);

  return ncm_sparam_get_upper_bound (sp);
}

/**
 * ncm_vparam_get_scale:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 *
 * Returns: The value of #NcmSParam:scale property of the @n-th component of @vparam.
 */
gdouble
ncm_vparam_get_scale (const NcmVParam *vparam, guint n)
{
  NcmSParam *sp = ncm_vparam_peek_sparam (vparam, n);

  return ncm_sparam_get_scale (sp);
}

/**
 * ncm_vparam_get_absolute_tolerance:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 *
 * Returns: The value of #NcmSParam:absolute-tolerance property of the @n-th component of @vparam.
 */
gdouble
ncm_vparam_get_absolute_tolerance (const NcmVParam *vparam, guint n)
{
  NcmSParam *sp = ncm_vparam_peek_sparam (vparam, n);

  return ncm_sparam_get_absolute_tolerance (sp);
}

/**
 * ncm_vparam_get_default_value:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 *
 * Returns: The value of #NcmSParam:default-value property of the @n-th component of @vparam.
 */
gdouble
ncm_vparam_get_default_value (const NcmVParam *vparam, guint n)
{
  NcmSParam *sp = ncm_vparam_peek_sparam (vparam, n);

  return ncm_sparam_get_default_value (sp);
}

/**
 * ncm_vparam_get_fit_type:
 * @vparam: a #NcmVParam.
 * @n: vector index.
 *
 * Returns: The value of #NcmSParam:fit-type property of the @n-th component of @vparam.
 */
NcmParamType
ncm_vparam_get_fit_type (const NcmVParam *vparam, guint n)
{
  NcmSParam *sp = ncm_vparam_peek_sparam (vparam, n);

  return ncm_sparam_get_fit_type (sp);
}

/**
 * ncm_vparam_len:
 * @vparam: a #NcmVParam.
 *
 * Returns: The length of @vparam.
 */
guint
ncm_vparam_len (const NcmVParam *vparam)
{
  return vparam->len;
}

