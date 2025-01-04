/***************************************************************************
 *            ncm_spline_gsl.c
 *
 *  Wed Nov 21 19:09:20 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
 * NcmSplineGsl:
 *
 * GSL spline object wrapper.
 *
 * This object comprises the proper functions to use the [GNU Scientific Library
 * (GSL)](https://www.gnu.org/software/gsl/) spline functions and interpolation methods.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_gsl.h"
#include "math/ncm_cfg.h"
#include "ncm_enum_types.h"

struct _NcmSplineGsl
{
  /*< private >*/
  NcmSpline parent_instance;
  gsl_interp *interp;
  NcmSplineGslType type_id;
  gchar *inst_name;
  const gsl_interp_type *type;
};

G_DEFINE_TYPE (NcmSplineGsl, ncm_spline_gsl, NCM_TYPE_SPLINE)

enum
{
  PROP_0,
  PROP_TYPE_ID,
  PROP_TYPE_NAME,
  PROP_SIZE,
};

static void
ncm_spline_gsl_init (NcmSplineGsl *sg)
{
  sg->interp    = NULL;
  sg->type      = NULL;
  sg->type_id   = NCM_SPLINE_GSL_TYPES_LEN;
  sg->inst_name = NULL;
}

static void
_ncm_spline_gsl_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSplineGsl *fit = NCM_SPLINE_GSL (object);

  g_return_if_fail (NCM_IS_SPLINE_GSL (object));

  switch (prop_id)
  {
    case PROP_TYPE_NAME:
      ncm_spline_gsl_set_type_by_name (fit, g_value_get_string (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_spline_gsl_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSplineGsl *sgsl = NCM_SPLINE_GSL (object);

  g_return_if_fail (NCM_IS_SPLINE_GSL (object));

  switch (prop_id)
  {
    case PROP_TYPE_NAME:
    {
      const GEnumValue *e = ncm_cfg_enum_get_value (NCM_TYPE_SPLINE_GSL_TYPE, sgsl->type_id);

      g_value_set_string (value, e->value_name);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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

static void _ncm_spline_gsl_reset (NcmSpline *s);
static const gchar *_ncm_spline_gsl_name (NcmSpline *s);
static void _ncm_spline_gsl_prepare (NcmSpline *s);
static gsize _ncm_spline_gsl_min_size (const NcmSpline *s);
static gdouble _ncm_spline_gsl_eval (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_gsl_deriv (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_gsl_deriv2 (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_gsl_deriv_nmax (const NcmSpline *s, const gdouble x);
static gdouble _ncm_spline_gsl_integ (const NcmSpline *s, const gdouble x0, const gdouble x1);
static NcmSpline *_ncm_spline_gsl_copy_empty (const NcmSpline *s);

static void
ncm_spline_gsl_class_init (NcmSplineGslClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmSplineClass *s_class    = NCM_SPLINE_CLASS (klass);

  object_class->set_property = &_ncm_spline_gsl_set_property;
  object_class->get_property = &_ncm_spline_gsl_get_property;
  object_class->finalize     = &ncm_spline_gsl_finalize;

  /**
   * NcmSplineGsl:type-name:
   *
   * The name of the interpolation method from [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/).
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_TYPE_NAME,
                                   g_param_spec_string ("type-name",
                                                        NULL,
                                                        "GSL Interpolation method name",
                                                        "NCM_SPLINE_GSL_CSPLINE",
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  s_class->name         = &_ncm_spline_gsl_name;
  s_class->reset        = &_ncm_spline_gsl_reset;
  s_class->prepare      = &_ncm_spline_gsl_prepare;
  s_class->prepare_base = NULL;
  s_class->min_size     = &_ncm_spline_gsl_min_size;
  s_class->eval         = &_ncm_spline_gsl_eval;
  s_class->deriv        = &_ncm_spline_gsl_deriv;
  s_class->deriv2       = &_ncm_spline_gsl_deriv2;
  s_class->deriv_nmax   = &_ncm_spline_gsl_deriv_nmax;
  s_class->integ        = &_ncm_spline_gsl_integ;
  s_class->copy_empty   = &_ncm_spline_gsl_copy_empty;
}

static void
_ncm_spline_gsl_reset (NcmSpline *s)
{
  NcmSplineGsl *sg  = NCM_SPLINE_GSL (s);
  const guint s_len = ncm_spline_get_len (s);

  if (sg->interp != NULL)
  {
    if (sg->interp->size != s_len)
    {
      gsl_interp_free (sg->interp);
      sg->interp = gsl_interp_alloc (sg->type, s_len);
    }
  }
  else
  {
    sg->interp = gsl_interp_alloc (sg->type, s_len);
    g_free (sg->inst_name);
    sg->inst_name = g_strdup_printf ("NcmSplineGsl[%s]", gsl_interp_name (sg->interp));
  }
}

static const gchar *
_ncm_spline_gsl_name (NcmSpline *s)
{
  NcmSplineGsl *sg = NCM_SPLINE_GSL (s);

  return sg->inst_name;
}

static void
_ncm_spline_gsl_prepare (NcmSpline *s)
{
  NcmSplineGsl *sg  = NCM_SPLINE_GSL (s);
  const guint s_len = ncm_spline_get_len (s);
  NcmVector *s_xv   = ncm_spline_peek_xv (s);
  NcmVector *s_yv   = ncm_spline_peek_yv (s);

  g_assert_cmpint (ncm_vector_stride (s_xv), ==, 1);
  g_assert_cmpint (ncm_vector_stride (s_yv), ==, 1);

  gsl_interp_init (sg->interp, ncm_vector_ptr (s_xv, 0), ncm_vector_ptr (s_yv, 0), s_len);
}

static gsize
_ncm_spline_gsl_min_size (const NcmSpline *s)
{
  NcmSplineGsl *sg = NCM_SPLINE_GSL ((NcmSpline *) s);

  return sg->type->min_size;
}

static gdouble
_ncm_spline_gsl_eval (const NcmSpline *s, const gdouble x)
{
  NcmSplineGsl *sg        = NCM_SPLINE_GSL ((NcmSpline *) s);
  NcmVector *s_xv         = ncm_spline_peek_xv ((NcmSpline *) s);
  NcmVector *s_yv         = ncm_spline_peek_yv ((NcmSpline *) s);
  gsl_interp_accel *s_acc = ncm_spline_peek_acc ((NcmSpline *) s);

  return gsl_interp_eval (sg->interp, ncm_vector_ptr (s_xv, 0), ncm_vector_ptr (s_yv, 0), x, s_acc);
}

static gdouble
_ncm_spline_gsl_deriv (const NcmSpline *s, const gdouble x)
{
  NcmSplineGsl *sg        = NCM_SPLINE_GSL ((NcmSpline *) s);
  NcmVector *s_xv         = ncm_spline_peek_xv ((NcmSpline *) s);
  NcmVector *s_yv         = ncm_spline_peek_yv ((NcmSpline *) s);
  gsl_interp_accel *s_acc = ncm_spline_peek_acc ((NcmSpline *) s);

  return gsl_interp_eval_deriv (sg->interp, ncm_vector_ptr (s_xv, 0), ncm_vector_ptr (s_yv, 0), x, s_acc);
}

static gdouble
_ncm_spline_gsl_deriv2 (const NcmSpline *s, const gdouble x)
{
  NcmSplineGsl *sg        = NCM_SPLINE_GSL ((NcmSpline *) s);
  NcmVector *s_xv         = ncm_spline_peek_xv ((NcmSpline *) s);
  NcmVector *s_yv         = ncm_spline_peek_yv ((NcmSpline *) s);
  gsl_interp_accel *s_acc = ncm_spline_peek_acc ((NcmSpline *) s);

  return gsl_interp_eval_deriv2 (sg->interp, ncm_vector_ptr (s_xv, 0), ncm_vector_ptr (s_yv, 0), x, s_acc);
}

static gdouble
_ncm_spline_gsl_deriv_nmax (const NcmSpline *s, const gdouble x)
{
  NcmSplineGsl *sg        = NCM_SPLINE_GSL ((NcmSpline *) s);
  NcmVector *s_xv         = ncm_spline_peek_xv ((NcmSpline *) s);
  NcmVector *s_yv         = ncm_spline_peek_yv ((NcmSpline *) s);
  gsl_interp_accel *s_acc = ncm_spline_peek_acc ((NcmSpline *) s);

  if (sg->type == gsl_interp_linear)
  {
    return gsl_interp_eval_deriv (sg->interp, ncm_vector_ptr (s_xv, 0), ncm_vector_ptr (s_yv, 0), x, s_acc);
  }
  else if ((sg->type == gsl_interp_cspline) || (sg->type == gsl_interp_cspline_periodic) ||
           (sg->type == gsl_interp_akima) || (sg->type == gsl_interp_akima_periodic))
  {
    const guint knot_i        = ncm_spline_get_index (s, x);
    const gdouble x_i         = ncm_vector_get (s_xv, knot_i);
    const gdouble x_ip1       = ncm_vector_get (s_xv, knot_i + 1);
    const gdouble dx          = x_ip1 - x_i;
    gdouble two_c_i           = gsl_interp_eval_deriv2 (sg->interp, ncm_vector_ptr (s_xv, 0), ncm_vector_ptr (s_yv, 0), x_i, s_acc);
    gdouble two_c_i_p_6d_i_dx = gsl_interp_eval_deriv2 (sg->interp, ncm_vector_ptr (s_xv, 0), ncm_vector_ptr (s_yv, 0), x_ip1, s_acc);

    return (two_c_i_p_6d_i_dx - two_c_i) / dx;
  }
  else
  {
    g_error ("ncm_spline_gsl_deriv_nmax: Calculation of the nmax derivative not supported.");

    return 0.0;
  }
}

static gdouble
_ncm_spline_gsl_integ (const NcmSpline *s, const gdouble x0, const gdouble x1)
{
  NcmSplineGsl *sg        = NCM_SPLINE_GSL ((NcmSpline *) s);
  NcmVector *s_xv         = ncm_spline_peek_xv ((NcmSpline *) s);
  NcmVector *s_yv         = ncm_spline_peek_yv ((NcmSpline *) s);
  gsl_interp_accel *s_acc = ncm_spline_peek_acc ((NcmSpline *) s);

  return gsl_interp_eval_integ (sg->interp, ncm_vector_ptr (s_xv, 0), ncm_vector_ptr (s_yv, 0), x0, x1, s_acc);
}

static NcmSpline *
_ncm_spline_gsl_copy_empty (const NcmSpline *s)
{
  NcmSplineGsl *sg = NCM_SPLINE_GSL ((NcmSpline *) s);

  return NCM_SPLINE (ncm_spline_gsl_new (sg->type));
}

/**
 * ncm_spline_gsl_new:
 * @type: gsl interpolation method
 *
 * This function returns a new gsl #NcmSpline which will use @type
 * interpolation method.
 *
 * Returns: a new #NcmSpline.
 */
NcmSplineGsl *
ncm_spline_gsl_new (const gsl_interp_type *type)
{
  NcmSplineGsl *sg = g_object_new (NCM_TYPE_SPLINE_GSL, NULL);

  ncm_spline_gsl_set_type (sg, type);

  return sg;
}

/**
 * ncm_spline_gsl_new_by_id:
 * @type_id: gsl interpolation method id
 *
 * This function returns a new gsl #NcmSpline which will use @type
 * interpolation method.
 *
 * Returns: a new #NcmSpline.
 */
NcmSplineGsl *
ncm_spline_gsl_new_by_id (NcmSplineGslType type_id)
{
  NcmSplineGsl *sg = g_object_new (NCM_TYPE_SPLINE_GSL, NULL);

  ncm_spline_gsl_set_type_by_id (sg, type_id);

  return sg;
}

/**
 * ncm_spline_gsl_new_full:
 * @type: gsl interpolation method
 * @xv: #NcmVector of knots
 * @yv: #NcmVector of the values of the function, to be interpolated, computed at @xv
 * @init: TRUE to prepare the new #NcmSpline or FALSE to not prepare it
 *
 * This function returns a new gsl #NcmSpline setting all its members.
 *
 * Returns: a new #NcmSpline.
 */
NcmSplineGsl *
ncm_spline_gsl_new_full (const gsl_interp_type *type, NcmVector *xv, NcmVector *yv, gboolean init)
{
  NcmSplineGsl *s = ncm_spline_gsl_new (type);

  ncm_spline_set (NCM_SPLINE (s), xv, yv, init);

  return s;
}

/**
 * ncm_spline_gsl_new_full_by_id:
 * @type_id: gsl interpolation method id
 * @xv: #NcmVector of knots
 * @yv: #NcmVector of the values of the function, to be interpolated, computed at @xv
 * @init: TRUE to prepare the new #NcmSpline or FALSE to not prepare it
 *
 * This function returns a new gsl #NcmSplineGsl setting all its members.
 *
 * Returns: a new #NcmSplineGsl.
 */
NcmSplineGsl *
ncm_spline_gsl_new_full_by_id (NcmSplineGslType type_id, NcmVector *xv, NcmVector *yv, gboolean init)
{
  NcmSplineGsl *s = ncm_spline_gsl_new_by_id (type_id);

  ncm_spline_set (NCM_SPLINE (s), xv, yv, init);

  return s;
}

/**
 * ncm_spline_gsl_set_type:
 * @sg: a #NcmSplineGsl
 * @type: gsl interpolation method
 *
 * This function sets the interpolation method @type to @sg.
 *
 */
void
ncm_spline_gsl_set_type (NcmSplineGsl *sg, const gsl_interp_type *type)
{
  const GEnumValue *type_id = ncm_cfg_get_enum_by_id_name_nick (NCM_TYPE_SPLINE_GSL_TYPE, type->name);

  if (sg->interp != NULL)
  {
    if (sg->type != type)
    {
      sg->type    = type;
      sg->type_id = type_id->value;
      g_free (sg->inst_name);
      sg->inst_name = g_strdup_printf ("NcmSplineGsl[%s]", type->name);

      gsl_interp_free (sg->interp);
      sg->interp = NULL;
      _ncm_spline_gsl_reset (NCM_SPLINE (sg));
    }
  }
  else
  {
    sg->type    = type;
    sg->type_id = type_id->value;
    g_free (sg->inst_name);
    sg->inst_name = g_strdup_printf ("NcmSplineGsl[%s]", type->name);
  }
}

/**
 * ncm_spline_gsl_set_type_by_id:
 * @sg: a #NcmSplineGsl
 * @type_id: gsl interpolation method id
 *
 * This function sets the interpolation method @type_id to @sg.
 *
 */
void
ncm_spline_gsl_set_type_by_id (NcmSplineGsl *sg, NcmSplineGslType type_id)
{
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
 * @sg: a #NcmSplineGsl
 * @type_name: gsl interpolation method name
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

/**
 * ncm_spline_gsl_get_type_id:
 * @sg: a #NcmSplineGsl
 *
 * This function returns the interpolation method id of @sg.
 *
 * Returns: the interpolation method id.
 */
NcmSplineGslType
ncm_spline_gsl_get_type_id (NcmSplineGsl *sg)
{
  return sg->type_id;
}

/**
 * ncm_spline_gsl_get_gsl_type:
 * @sg: a #NcmSplineGsl
 *
 * This function returns the interpolation method of @sg.
 *
 * Returns: the gsl interpolation method.
 */
const gsl_interp_type *
ncm_spline_gsl_get_gsl_type (NcmSplineGsl *sg)
{
  return sg->type;
}

