/***************************************************************************
 *            ncm_integral1d_ptr.c
 *
 *  Mon December 05 17:23:42 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_integral1d_ptr.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmIntegral1dPtr:
 *
 * One-dimensional integral of an integrand given by a function pointer.
 *
 * The integrand is a #NcmIntegral1dPtrF called with #NcmIntegral1dPtr:userdata. The
 * function pointer is not serialized, so this type cannot be serialized.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/integration/ncm_integral1d_ptr.h"

typedef struct _NcmIntegral1dPtrPrivate
{
  NcmIntegral1dF F;
  gpointer userdata;
  GDestroyNotify userfree;
} NcmIntegral1dPtrPrivate;

struct _NcmIntegral1dPtr
{
  NcmIntegral1d parent_instance;
};

enum
{
  PROP_0,
  PROP_INTEGRAND,
  PROP_USERDATA,
  PROP_USERFREE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmIntegral1dPtr, ncm_integral1d_ptr, NCM_TYPE_INTEGRAL1D)

static void
ncm_integral1d_ptr_init (NcmIntegral1dPtr *int1d_ptr)
{
  NcmIntegral1dPtrPrivate * const self = ncm_integral1d_ptr_get_instance_private (int1d_ptr);

  self->F        = NULL;
  self->userdata = NULL;
  self->userfree = NULL;
}

static void
_ncm_integral1d_ptr_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmIntegral1dPtr *int1d_ptr          = NCM_INTEGRAL1D_PTR (object);
  NcmIntegral1dPtrPrivate * const self = ncm_integral1d_ptr_get_instance_private (int1d_ptr);

  g_return_if_fail (NCM_IS_INTEGRAL1D_PTR (object));

  switch (prop_id)
  {
    case PROP_INTEGRAND:
      self->F = g_value_get_pointer (value);
      break;
    case PROP_USERDATA:
      ncm_integral1d_ptr_set_userdata (int1d_ptr, g_value_get_pointer (value));
      break;
    case PROP_USERFREE:
      self->userfree = g_value_get_pointer (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_integral1d_ptr_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmIntegral1dPtr *int1d_ptr          = NCM_INTEGRAL1D_PTR (object);
  NcmIntegral1dPtrPrivate * const self = ncm_integral1d_ptr_get_instance_private (int1d_ptr);

  g_return_if_fail (NCM_IS_INTEGRAL1D_PTR (object));

  switch (prop_id)
  {
    case PROP_INTEGRAND:
      g_value_set_pointer (value, self->F);
      break;
    case PROP_USERDATA:
      g_value_set_pointer (value, self->userdata);
      break;
    case PROP_USERFREE:
      g_value_set_pointer (value, self->userfree);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_integral1d_ptr_finalize (GObject *object)
{
  NcmIntegral1dPtr *int1d_ptr          = NCM_INTEGRAL1D_PTR (object);
  NcmIntegral1dPtrPrivate * const self = ncm_integral1d_ptr_get_instance_private (int1d_ptr);

  if ((self->userdata != NULL) && (self->userfree != NULL))
    self->userfree (self->userdata);

  self->userdata = NULL;

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_integral1d_ptr_parent_class)->finalize (object);
}

static gdouble _ncm_integral1d_ptr_integrand (NcmIntegral1d *int1d, const gdouble x, const gdouble w);

static void
ncm_integral1d_ptr_class_init (NcmIntegral1dPtrClass *klass)
{
  GObjectClass *object_class      = G_OBJECT_CLASS (klass);
  NcmIntegral1dClass *int1d_class = NCM_INTEGRAL1D_CLASS (klass);

  object_class->set_property = &_ncm_integral1d_ptr_set_property;
  object_class->get_property = &_ncm_integral1d_ptr_get_property;
  object_class->finalize     = &_ncm_integral1d_ptr_finalize;

  /**
   * NcmIntegral1dPtr:integrand:
   *
   * The integrand, a #NcmIntegral1dPtrF.
   */
  g_object_class_install_property (object_class,
                                   PROP_INTEGRAND,
                                   g_param_spec_pointer ("integrand",
                                                         NULL,
                                                         "Integrand function pointer",
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmIntegral1dPtr:userdata:
   *
   * The user data passed to the integrand.
   */
  g_object_class_install_property (object_class,
                                   PROP_USERDATA,
                                   g_param_spec_pointer ("userdata",
                                                         NULL,
                                                         "Integrand function user data",
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmIntegral1dPtr:userfree:
   *
   * The function that frees #NcmIntegral1dPtr:userdata when it is replaced and on
   * finalization, or %NULL.
   */
  g_object_class_install_property (object_class,
                                   PROP_USERFREE,
                                   g_param_spec_pointer ("userfree",
                                                         NULL,
                                                         "Integrand function user data free function",
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  int1d_class->integrand = &_ncm_integral1d_ptr_integrand;
}

static gdouble
_ncm_integral1d_ptr_integrand (NcmIntegral1d *int1d, const gdouble x, const gdouble w)
{
  NcmIntegral1dPtr *int1d_ptr          = NCM_INTEGRAL1D_PTR (int1d);
  NcmIntegral1dPtrPrivate * const self = ncm_integral1d_ptr_get_instance_private (int1d_ptr);

  return self->F (self->userdata, x, w);
}

/**
 * ncm_integral1d_ptr_new:
 * @F: (scope notified): the integrand
 * @userfree: (scope notified): the function that frees the user data
 *
 * Creates a #NcmIntegral1dPtr with the default #NcmIntegral1d configuration.
 *
 * Returns: (transfer full): a new #NcmIntegral1dPtr.
 */
NcmIntegral1dPtr *
ncm_integral1d_ptr_new (NcmIntegral1dPtrF F, GDestroyNotify userfree)
{
  NcmIntegral1dPtr *int1d_ptr = g_object_new (NCM_TYPE_INTEGRAL1D_PTR,
                                              "integrand", F,
                                              "userfree",  userfree,
                                              NULL);

  return int1d_ptr;
}

/**
 * ncm_integral1d_ptr_new_full:
 * @F: (scope notified): the integrand
 * @userfree: (scope notified): the function that frees the user data
 * @reltol: the relative tolerance
 * @abstol: the absolute tolerance
 * @partition: the maximum number of subintervals
 * @rule: the Gauss-Kronrod rule
 *
 * Creates a #NcmIntegral1dPtr, see #NcmIntegral1d for the configuration.
 *
 * Returns: (transfer full): a new #NcmIntegral1dPtr.
 */
NcmIntegral1dPtr *
ncm_integral1d_ptr_new_full (NcmIntegral1dPtrF F, GDestroyNotify userfree, gdouble reltol, gdouble abstol, guint partition, guint rule)
{
  NcmIntegral1dPtr *int1d_ptr = g_object_new (NCM_TYPE_INTEGRAL1D_PTR,
                                              "integrand", F,
                                              "userfree", userfree,
                                              "reltol",    reltol,
                                              "abstol",    abstol,
                                              "partition", partition,
                                              "rule",      rule,
                                              NULL);

  return int1d_ptr;
}

/**
 * ncm_integral1d_ptr_ref:
 * @int1d_ptr: a #NcmIntegral1dPtr
 *
 * Increases the reference count of @int1d_ptr by one.
 *
 * Returns: (transfer full): @int1d_ptr.
 */
NcmIntegral1dPtr *
ncm_integral1d_ptr_ref (NcmIntegral1dPtr *int1d_ptr)
{
  return g_object_ref (int1d_ptr);
}

/**
 * ncm_integral1d_ptr_free:
 * @int1d_ptr: a #NcmIntegral1dPtr
 *
 * Decreases the reference count of @int1d_ptr by one.
 */
void
ncm_integral1d_ptr_free (NcmIntegral1dPtr *int1d_ptr)
{
  g_object_unref (int1d_ptr);
}

/**
 * ncm_integral1d_ptr_clear:
 * @int1d_ptr: a #NcmIntegral1dPtr
 *
 * If *@int1d_ptr is not %NULL, decreases its reference count by one and sets *@int1d_ptr
 * to %NULL.
 */
void
ncm_integral1d_ptr_clear (NcmIntegral1dPtr **int1d_ptr)
{
  g_clear_object (int1d_ptr);
}

/**
 * ncm_integral1d_ptr_set_userdata:
 * @int1d_ptr: a #NcmIntegral1dPtr
 * @userdata: the user data
 *
 * Sets #NcmIntegral1dPtr:userdata, freeing the previous one with
 * #NcmIntegral1dPtr:userfree.
 */
void
ncm_integral1d_ptr_set_userdata (NcmIntegral1dPtr *int1d_ptr, gpointer userdata)
{
  NcmIntegral1dPtrPrivate * const self = ncm_integral1d_ptr_get_instance_private (int1d_ptr);

  if ((self->userdata != NULL) && (self->userfree != NULL))
  {
    self->userfree (self->userdata);
    self->userdata = NULL;
  }

  self->userdata = userdata;
}

