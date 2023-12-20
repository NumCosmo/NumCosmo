/***************************************************************************
 *            ncm_reparam_linear.c
 *
 *  Thu March 08 11:05:07 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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
 * SECTION:ncm_reparam_linear
 * @title: NcmReparamLinear
 * @short_description: Linear reparametrization object.
 *
 * Object implementing a linear reparametrization of the model's parameters.
 * It uses as imput a matrix $M$ (#NcmReparamLinear:matrix) and a vector $v$
 * (#NcmReparamLinear:vector), such that the new parameters
 * $\vec{w}$ are given by $$\vec{w} = M\cdot\vec{y} + \vec{v},$$ where $\vec{y}$
 * represents the original model's parameters.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_reparam_linear.h"
#include "math/ncm_vector.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmReparamLinear
{
  /*< private >*/
  NcmReparam parent_instance;
  NcmMatrix *T;
  NcmVector *v;
  NcmVector *vp;
  NcmMatrix *T_LU;
  gsl_permutation *p;
  gint signum;
};

enum
{
  PROP_0,
  PROP_V,
  PROP_T
};

G_DEFINE_TYPE (NcmReparamLinear, ncm_reparam_linear, NCM_TYPE_REPARAM)

static void
ncm_reparam_linear_init (NcmReparamLinear *relin)
{
  relin->T      = NULL;
  relin->v      = NULL;
  relin->T_LU   = NULL;
  relin->p      = NULL;
  relin->signum = 0;
}

static void
_ncm_reparam_linear_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmReparamLinear *relin = NCM_REPARAM_LINEAR (object);

  g_return_if_fail (NCM_IS_REPARAM_LINEAR (object));

  switch (prop_id)
  {
    case PROP_V:
      g_value_set_object (value, relin->v);
      break;
    case PROP_T:
      g_value_set_object (value, relin->T);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_reparam_linear_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmReparamLinear *relin = NCM_REPARAM_LINEAR (object);

  g_return_if_fail (NCM_IS_REPARAM_LINEAR (object));

  switch (prop_id)
  {
    case PROP_V:
      ncm_vector_substitute (&relin->v, g_value_get_object (value), TRUE);
      break;
    case PROP_T:
      ncm_matrix_substitute (&relin->T, g_value_get_object (value), TRUE);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_reparam_linear_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_reparam_linear_parent_class)->constructed (object);
  {
    NcmReparamLinear *relin = NCM_REPARAM_LINEAR (object);
    NcmReparam *reparam     = NCM_REPARAM (relin);
    const guint length      = ncm_reparam_get_length (reparam);

    g_assert ((ncm_matrix_nrows (relin->T) == ncm_matrix_ncols (relin->T)) &&
              (ncm_matrix_nrows (relin->T) == length) &&
              (ncm_vector_len (relin->v) == length));

    relin->vp   = ncm_vector_new (length);
    relin->T_LU = ncm_matrix_new (length, length);
    relin->p    = gsl_permutation_alloc (length);

    ncm_matrix_memcpy (relin->T_LU, relin->T);
    gsl_linalg_LU_decomp (ncm_matrix_gsl (relin->T_LU), relin->p, &relin->signum);
    gsl_linalg_LU_solve (ncm_matrix_gsl (relin->T_LU), relin->p, ncm_vector_gsl (relin->v), ncm_vector_gsl (relin->vp));
  }
}

static void
_ncm_reparam_linear_dispose (GObject *object)
{
  NcmReparamLinear *relin = NCM_REPARAM_LINEAR (object);

  ncm_matrix_clear (&relin->T_LU);
  ncm_matrix_clear (&relin->T);
  ncm_vector_clear (&relin->v);
  ncm_vector_clear (&relin->vp);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_reparam_linear_parent_class)->dispose (object);
}

static void
_ncm_reparam_linear_finalize (GObject *object)
{
  NcmReparamLinear *relin = NCM_REPARAM_LINEAR (object);

  gsl_permutation_free (relin->p);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_reparam_linear_parent_class)->finalize (object);
}

static gboolean _ncm_reparam_linear_old2new (NcmReparam *reparam, NcmModel *model);
static gboolean _ncm_reparam_linear_new2old (NcmReparam *reparam, NcmModel *model);

static void
ncm_reparam_linear_class_init (NcmReparamLinearClass *klass)
{
  GObjectClass *object_class    = G_OBJECT_CLASS (klass);
  NcmReparamClass *parent_class = NCM_REPARAM_CLASS (klass);

  object_class->set_property = &_ncm_reparam_linear_set_property;
  object_class->get_property = &_ncm_reparam_linear_get_property;
  object_class->constructed  = &_ncm_reparam_linear_constructed;
  object_class->dispose      = &_ncm_reparam_linear_dispose;
  object_class->finalize     = &_ncm_reparam_linear_finalize;

  /**
   * NcmReparamLinear:vector:
   *
   * The vector $\vec{v}$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_V,
                                   g_param_spec_object ("vector",
                                                        NULL,
                                                        "Vector of parameters shift",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmReparamLinear:matrix:
   *
   * The matrix $M$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_T,
                                   g_param_spec_object ("matrix",
                                                        NULL,
                                                        "Matrix of parameters mixing",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  parent_class->old2new = &_ncm_reparam_linear_old2new;
  parent_class->new2old = &_ncm_reparam_linear_new2old;
}

static gboolean
_ncm_reparam_linear_old2new (NcmReparam *reparam, NcmModel *model)
{
  NcmReparamLinear *relin = NCM_REPARAM_LINEAR (reparam);
  NcmVector *params       = ncm_model_orig_params_peek_vector (model);
  NcmVector *new_params   = ncm_reparam_peek_params (reparam);

  ncm_vector_memcpy (new_params, relin->v);
  gsl_blas_dgemv (CblasNoTrans, 1.0, ncm_matrix_gsl (relin->T), ncm_vector_gsl (params), 1.0, ncm_vector_gsl (new_params));

  return TRUE;
}

static gboolean
_ncm_reparam_linear_new2old (NcmReparam *reparam, NcmModel *model)
{
  NcmReparamLinear *relin = NCM_REPARAM_LINEAR (reparam);
  NcmVector *params       = ncm_model_orig_params_peek_vector (model);
  NcmVector *new_params   = ncm_reparam_peek_params (reparam);

  gsl_linalg_LU_solve (ncm_matrix_gsl (relin->T_LU), relin->p, ncm_vector_gsl (new_params), ncm_vector_gsl (params));
  ncm_vector_sub (params, relin->vp);

  return TRUE;
}

/**
 * ncm_reparam_linear_new:
 * @size: model's length.
 * @T: a #NcmMatrix
 * @v: a #NcmVector
 *
 * Creates a new reparametrization using the parameters transformation matrix
 * @T and the shift vector @v, i.e., the new parameters vector $\vec{p}_n$ is
 * given by $\vec{p}_n = T\cdot{}\vec{p} + \vec{v}$, where $p$ are the old
 * parameters vector.
 *
 * Returns: (transfer full): a new #NcmReparamLinear.
 */
NcmReparamLinear *
ncm_reparam_linear_new (guint size, NcmMatrix *T, NcmVector *v)
{
  NcmReparamLinear *relin = g_object_new (NCM_TYPE_REPARAM_LINEAR,
                                          "length", size,
                                          "matrix", T,
                                          "vector", v,
                                          NULL);

  return relin;
}

/**
 * ncm_reparam_linear_set_compat_type:
 * @lin: a #NcmReparamLinear
 * @compat_type: a #GType
 *
 * Sets the object's type compatible with this reparametrization.
 *
 */
void
ncm_reparam_linear_set_compat_type (NcmReparamLinear *lin, GType compat_type)
{
  ncm_reparam_set_compat_type (NCM_REPARAM (lin), compat_type);
}

