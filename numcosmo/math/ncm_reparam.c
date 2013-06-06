/***************************************************************************
 *            ncm_reparam.c
 *
 *  Thu March 08 00:36:24 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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
 * SECTION:ncm_reparam
 * @title: Reparametrization Abstract Class
 * @short_description: Base class for model reparametrization
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_reparam.h"
#include "math/ncm_model.h"
#include "math/ncm_cfg.h"

#include <gsl/gsl_blas.h>

G_DEFINE_ABSTRACT_TYPE (NcmReparam, ncm_reparam, G_TYPE_OBJECT);

/**
 * ncm_reparam_copyto:
 * @reparam: a #NcmReparam
 * @reparam_dest: a #NcmReparam
 *
 * FIXME
 */
void
ncm_reparam_copyto (NcmReparam *reparam, NcmReparam *reparam_dest)
{
  NCM_REPARAM_GET_CLASS (reparam)->copyto (reparam, reparam_dest);
}

/**
 * ncm_reparam_copy: (skip)
 * @reparam: a #NcmReparam
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmReparam *
ncm_reparam_copy (NcmReparam *reparam)
{
  if (!NCM_REPARAM_GET_CLASS (reparam)->copy)
    g_error ("NcmReparam[%s] base class do not implement copy.", G_OBJECT_TYPE_NAME (reparam));
  else
    return NCM_REPARAM_GET_CLASS (reparam)->copy (reparam);
}

/**
 * ncm_reparam_ref:
 * @reparam: a #NcmReparam
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
   */
NcmReparam *
ncm_reparam_ref (NcmReparam *reparam)
{
  return NCM_REPARAM (g_object_ref (reparam));
}

enum
{
  PROP_0,
  PROP_LEN,
};

/**
 * ncm_reparam_old2new:
 * @reparam: a #NcmReparam
 * @model: a #NcmModel
 * @src: a #NcmVector
 * @dest: a #NcmVector
 *
 * FIXME
 */
void
ncm_reparam_old2new (NcmReparam *reparam, NcmModel *model, NcmVector *src, NcmVector *dest)
{
  NCM_REPARAM_GET_CLASS (reparam)->old2new (reparam, model, src, dest);
}

/**
 * ncm_reparam_new2old:
 * @reparam: a #NcmReparam
 * @model: a #NcmModel
 * @src: a #NcmVector
 * @dest: a #NcmVector
 *
 * FIXME
 */
void
ncm_reparam_new2old (NcmReparam *reparam, NcmModel *model, NcmVector *src, NcmVector *dest)
{
  NCM_REPARAM_GET_CLASS (reparam)->new2old (reparam, model, src, dest);
}

/**
 * ncm_reparam_jac:
 * @reparam: a #NcmReparam
 * @model: a #NcmModel
 * @jac: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_reparam_jac (NcmReparam *reparam, NcmModel *model, NcmMatrix *jac)
{
  NCM_REPARAM_GET_CLASS (reparam)->jac (reparam, model, jac);
}

/**
 * ncm_reparam_grad_old2new:
 * @reparam: a #NcmReparam
 * @model: FIXME
 * @jac: a #NcmMatrix
 * @old_grad: a #NcmVector
 * @new_grad: a #NcmVector
 *
 * FIXME
 */
void
ncm_reparam_grad_old2new (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac, NcmVector *old_grad, NcmVector *new_grad)
{
  gint ret;
  NCM_REPARAM_GET_CLASS (reparam)->jac (reparam, model, jac);
  ret = gsl_blas_dgemv (CblasTrans, 1.0, NCM_MATRIX_GSL (jac), ncm_vector_gsl (old_grad), 0.0, ncm_vector_gsl (new_grad));
  NCM_TEST_GSL_RESULT("ncm_reparam_grad_old2new", ret);
}

/**
 * ncm_reparam_M_old2new:
 * @reparam: a #NcmReparam
 * @model: FIXME
 * @jac: a #NcmMatrix
 * @old_M: a #NcmMatrix
 * @new_M: a #NcmMatrix
 *
 * FIXME
 */
void
ncm_reparam_M_old2new (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac, NcmMatrix *old_M, NcmMatrix *new_M)
{
  gint ret;
  NCM_REPARAM_GET_CLASS (reparam)->jac (reparam, model, jac);
  ret = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, NCM_MATRIX_GSL (old_M), NCM_MATRIX_GSL (jac), 0.0, NCM_MATRIX_GSL (new_M));
  NCM_TEST_GSL_RESULT("ncm_reparam_jac_old2new", ret);
}

/**
 * ncm_reparam_free:
 * @reparam: a #NcmReparam
 *
 * FIXME
 */
void
ncm_reparam_free (NcmReparam *reparam)
{
  g_object_unref (reparam);
}

/**
 * ncm_reparam_clear:
 * @reparam: a #NcmReparam
 *
 * FIXME
 */
void
ncm_reparam_clear (NcmReparam **reparam)
{
  g_clear_object (reparam);
}

static void
_ncm_reparam_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_reparam_parent_class)->constructed (object);
  {
    NcmReparam *reparam = NCM_REPARAM (object);

    reparam->new_params = ncm_vector_new (reparam->length);
    reparam->sparams = g_ptr_array_sized_new (reparam->length);
    g_ptr_array_set_size (reparam->sparams, reparam->length);
  }
}

static void
_ncm_reparam_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmReparam *reparam = NCM_REPARAM (object);

  g_return_if_fail (NCM_IS_REPARAM (object));

  switch (prop_id)
  {
    case PROP_LEN:
      g_value_set_uint (value, reparam->length);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_reparam_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmReparam *reparam = NCM_REPARAM (object);

  g_return_if_fail (NCM_IS_REPARAM (object));

  switch (prop_id)
  {
    case PROP_LEN:
      reparam->length = g_value_get_uint (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_reparam_init (NcmReparam *reparam)
{
  reparam->length = 0;
  reparam->sparams = NULL;
  reparam->new_params = NULL;
}

static void
ncm_reparam_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_reparam_parent_class)->finalize (object);
}

static void
_ncm_reparam_copyto (NcmReparam *reparam, NcmReparam *reparam_dest)
{
  guint i;
  g_assert (G_OBJECT_TYPE (reparam) == G_OBJECT_TYPE (reparam_dest));
  if (reparam->length != reparam_dest->length)
  {
    ncm_vector_free (reparam_dest->new_params);
    reparam_dest->new_params = ncm_vector_new (reparam->length);
    g_ptr_array_set_size (reparam_dest->sparams, 0);
    g_ptr_array_set_size (reparam_dest->sparams, reparam->length);
  }

  ncm_vector_memcpy (reparam_dest->new_params, reparam->new_params);
  if (reparam->sparams)
  {
    for (i = 0; i < reparam->length; i++)
    {
      NcmSParam *pinfo = g_ptr_array_index (reparam->sparams, i);
      if (pinfo)
        g_ptr_array_index (reparam_dest->sparams, i) = ncm_sparam_copy (pinfo);
    }
  }
}

static void
ncm_reparam_class_init (NcmReparamClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  klass->copyto = &_ncm_reparam_copyto;
  klass->copy = NULL;
  klass->old2new = NULL;
  klass->new2old = NULL;
  klass->jac = NULL;

  object_class->set_property = &_ncm_reparam_set_property;
  object_class->get_property = &_ncm_reparam_get_property;
  object_class->constructed  = &_ncm_reparam_constructed;

  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("length",
                                                      NULL,
                                                      "System's length",
                                                      0.0, G_MAXUINT, 0.0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  object_class->finalize = ncm_reparam_finalize;
}
