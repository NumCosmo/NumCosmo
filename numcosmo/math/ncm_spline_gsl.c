/***************************************************************************
 *            ncm_spline_gsl.c
 *
 *  Wed Nov 21 19:09:20 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:ncm_spline_gsl
 * @title: GSL Spline
 * @short_description: GSL spline object wrapper
 *
 * This object comprises the proper functions to use the GNU Scientific 
 * Library (GSL) spline functions and interpolation methods.
   * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_gsl.h"
#include "math/ncm_cfg.h"
#include "ncm_enum_types.h"

G_DEFINE_TYPE (NcmSplineGsl, ncm_spline_gsl, NCM_TYPE_SPLINE);

enum
{
  PROP_0,
  PROP_TYPE_ID,
  PROP_TYPE_NAME,
  PROP_SIZE,
};

static void _ncm_spline_gsl_reset (NcmSpline *s);

/**
 * ncm_spline_gsl_new: (skip)
 * @type: gsl interpolation method.
 * 
 * This function returns a new gsl #NcmSpline which will use @type 
 * interpolation method.
 * 
 * Returns: a new #NcmSpline.
 */
NcmSpline *
ncm_spline_gsl_new (const gsl_interp_type *type)
{
  NcmSplineGsl *sg = g_object_new (NCM_TYPE_SPLINE_GSL, NULL);
  ncm_spline_gsl_set_type (sg, type);
  return NCM_SPLINE (sg);
}

/**
 * ncm_spline_gsl_new_full: (skip)
 * @type: gsl interpolation method.
 * @xv: #NcmVector of knots.
 * @yv: #NcmVector of the values of the function, to be interpolated, computed at @xv.
 * @init: TRUE to prepare the new #NcmSpline or FALSE to not prepare it.
 * 
 * This function returns a new gsl #NcmSpline setting all its members. 
 * 
 * Returns: a new #NcmSpline.
 */
NcmSpline *
ncm_spline_gsl_new_full (const gsl_interp_type *type, NcmVector *xv, NcmVector *yv, gboolean init)
{
  NcmSpline *s = ncm_spline_gsl_new (type);
  ncm_spline_set (s, xv, yv, init);
  return s;
}

/**
 * ncm_spline_gsl_set_type: (skip)
 * @sg: a #NcmSplineGsl.
 * @type: gsl interpolation method.
 * 
 * This function sets the interpolation method @type to @sg.
 * 
 */
void
ncm_spline_gsl_set_type (NcmSplineGsl *sg, const gsl_interp_type *type)
{
  if (sg->interp != NULL)
  {
    if (sg->type != type)
    {
      sg->type = type;
      gsl_interp_free (sg->interp);
      sg->interp = NULL;
      _ncm_spline_gsl_reset (NCM_SPLINE (sg));
    }
  }
  else
  {
    sg->type = type;
    g_free (sg->inst_name);
    sg->inst_name = g_strdup_printf ("NcmSplineGsl[%s]", type->name);
  }
}

/**
 * ncm_spline_gsl_set_type_by_id:
 * @sg: a #NcmSplineGsl.
 * @type_id: gsl interpolation method id.
 * 
 * This function sets the interpolation method @type_id to @sg.
 * 
 */
void 
ncm_spline_gsl_set_type_by_id (NcmSplineGsl *sg, NcmSplineGslType type_id)
{
  sg->type_id = type_id;
  switch (type_id)
  {
    case NCM_SPLINE_GSL_LINEAR:
      ncm_spline_gsl_set_type (sg, gsl_interp_linear);
      break;
    case NCM_SPLINE_GSL_POLYNOMIAL:
      ncm_spline_gsl_set_type (sg, gsl_interp_polynomial);
      break;
    case NCM_SPLINE_GSL_CSPLINE:
      ncm_spline_gsl_set_type (sg, gsl_interp_cspline);
      break;
    case NCM_SPLINE_GSL_CSPLINE_PERIODIC:
      ncm_spline_gsl_set_type (sg, gsl_interp_cspline_periodic);
      break;
    case NCM_SPLINE_GSL_AKIMA:
      ncm_spline_gsl_set_type (sg, gsl_interp_akima);
      break;
    case NCM_SPLINE_GSL_AKIMA_PERIODIC:
      ncm_spline_gsl_set_type (sg, gsl_interp_akima_periodic);
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}

/**
 * ncm_spline_gsl_set_type_by_name:
 * @sg: a #NcmSplineGsl.
 * @type_name: gsl interpolation method name.
 * 
 * This function sets the interpolation method @type_name to @sg.
 * 
 */
void 
ncm_spline_gsl_set_type_by_name (NcmSplineGsl *sg, const gchar *type_name)
{
  const GEnumValue *type_id = ncm_cfg_get_enum_by_id_name_nick (NCM_TYPE_SPLINE_GSL_TYPE, type_name);
  if (type_id == NULL)
  {
    ncm_cfg_enum_print_all (NCM_TYPE_SPLINE_GSL_TYPE, "Error");
    g_error ("NcmSplineGsl type '%s' not found. Availables types above.", type_name);
  }
  ncm_spline_gsl_set_type_by_id (sg, type_id->value);
}

static NcmSpline *
_ncm_spline_gsl_copy_empty (const NcmSpline *s)
{
  NcmSplineGsl *sg = NCM_SPLINE_GSL (s);
  return ncm_spline_gsl_new (sg->type);
}

static const gchar *
_ncm_spline_gsl_name (NcmSpline *s)
{
  NcmSplineGsl *sg = NCM_SPLINE_GSL (s);
  return sg->inst_name;
}

static void 
_ncm_spline_gsl_reset (NcmSpline *s)
{ 
  NcmSplineGsl *sg = NCM_SPLINE_GSL (s);

  if (sg->interp != NULL)
  {
    if (sg->interp->size != s->len)
    {
      gsl_interp_free (sg->interp);
      sg->interp = gsl_interp_alloc (sg->type, s->len);
    }
  }
  else
  {
    sg->interp = gsl_interp_alloc (sg->type, s->len);
    g_free (sg->inst_name);
    sg->inst_name = g_strdup_printf ("NcmSplineGsl[%s]", gsl_interp_name (sg->interp));
  }
}

static void 
_ncm_spline_gsl_prepare (NcmSpline *s) 
{ 
  NcmSplineGsl *sg = NCM_SPLINE_GSL (s);
  gsl_interp_init (sg->interp, ncm_vector_ptr (s->xv, 0), ncm_vector_ptr (s->yv, 0), s->len); 
}

static gsize 
_ncm_spline_gsl_min_size (const NcmSpline *s)
{ 
  NcmSplineGsl *sg = NCM_SPLINE_GSL (s);
  return sg->type->min_size;
}

static gdouble 
_ncm_spline_gsl_eval (const NcmSpline *s, const gdouble x)
{ 
  NcmSplineGsl *sg = NCM_SPLINE_GSL (s);
  return gsl_interp_eval (sg->interp, ncm_vector_ptr (s->xv, 0), ncm_vector_ptr (s->yv, 0), x, s->acc); 
}

static gdouble 
_ncm_spline_gsl_deriv (const NcmSpline *s, const gdouble x) 
{ 
  NcmSplineGsl *sg = NCM_SPLINE_GSL (s);
  return gsl_interp_eval_deriv (sg->interp, ncm_vector_ptr (s->xv, 0), ncm_vector_ptr (s->yv, 0), x, s->acc);
}

static gdouble 
_ncm_spline_gsl_deriv2 (const NcmSpline *s, const gdouble x)
{ 
  NcmSplineGsl *sg = NCM_SPLINE_GSL (s);
  return gsl_interp_eval_deriv2 (sg->interp, ncm_vector_ptr (s->xv, 0), ncm_vector_ptr (s->yv, 0), x, s->acc); 
}

static gdouble 
_ncm_spline_gsl_integ (const NcmSpline *s, const gdouble x0, const gdouble x1)
{ 
  NcmSplineGsl *sg = NCM_SPLINE_GSL (s);
  return gsl_interp_eval_integ (sg->interp, ncm_vector_ptr (s->xv, 0), ncm_vector_ptr (s->yv, 0), x0, x1, s->acc); 
}

static void
ncm_spline_gsl_init (NcmSplineGsl *sg)
{
  sg->interp = NULL;
  sg->type = NULL;
  sg->type_id = NCM_SPLINE_GSL_TYPES_LEN;
  sg->inst_name = NULL;
}

static void
_ncm_spline_gsl_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSplineGsl *fit = NCM_SPLINE_GSL (object);
  g_return_if_fail (NCM_IS_SPLINE_GSL (object));

  switch (prop_id)
  {
    case PROP_TYPE_ID:
      ncm_spline_gsl_set_type_by_id (fit, g_value_get_enum (value));
      break;
    case PROP_TYPE_NAME:
      ncm_spline_gsl_set_type_by_name (fit, g_value_get_string (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_spline_gsl_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSplineGsl *sgsl = NCM_SPLINE_GSL (object);
  g_return_if_fail (NCM_IS_SPLINE_GSL (object));

  switch (prop_id)
  {
    case PROP_TYPE_ID:
      g_value_set_enum (value, sgsl->type_id);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_spline_gsl_finalize (GObject *object)
{
  NcmSplineGsl *sg = NCM_SPLINE_GSL (object);

  gsl_interp_free (sg->interp);
  sg->type = gsl_interp_linear;
  g_free (sg->inst_name);

  G_OBJECT_CLASS (ncm_spline_gsl_parent_class)->finalize (object);
}

static void
ncm_spline_gsl_class_init (NcmSplineGslClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmSplineClass *s_class = NCM_SPLINE_CLASS (klass);

  object_class->set_property = &_ncm_spline_gsl_set_property;
  object_class->get_property = &_ncm_spline_gsl_get_property;
  object_class->finalize     = &ncm_spline_gsl_finalize;
  
  g_object_class_install_property (object_class,
                                   PROP_TYPE_ID,
                                   g_param_spec_enum ("type",
                                                      NULL,
                                                      "GSL Interpolation method",
                                                      NCM_TYPE_SPLINE_GSL_TYPE, NCM_SPLINE_GSL_CSPLINE,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  

  g_object_class_install_property (object_class,
                                   PROP_TYPE_NAME,
                                   g_param_spec_string ("type-name",
                                                        NULL,
                                                        "GSL Interpolation method name",
                                                        NULL,
                                                        G_PARAM_WRITABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  
  
  s_class->name         = &_ncm_spline_gsl_name;
  s_class->reset        = &_ncm_spline_gsl_reset;
  s_class->prepare      = &_ncm_spline_gsl_prepare;
  s_class->prepare_base = NULL;
  s_class->min_size     = &_ncm_spline_gsl_min_size;
  s_class->eval         = &_ncm_spline_gsl_eval;
  s_class->deriv        = &_ncm_spline_gsl_deriv;
  s_class->deriv2       = &_ncm_spline_gsl_deriv2;
  s_class->integ        = &_ncm_spline_gsl_integ;
  s_class->copy_empty   = &_ncm_spline_gsl_copy_empty;
  
}
