/***************************************************************************
 *            ncm_reparam_linear.c
 *
 *  Thu March 08 11:05:07 2012
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
 * SECTION:ncm_reparam_linear
 * @title: Linear Reparametrization
 * @short_description: Linear reparametrization object
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

#include "ncm_reparam_linear.h"

G_DEFINE_TYPE (NcmReparamLinear, ncm_reparam_linear, NCM_TYPE_REPARAM);

NcmReparamLinear *
ncm_reparam_linear_new (guint size, NcmMatrix *T, NcmVector *v)
{
  GVariant *vvar = ncm_vector_peek_variant (v);
  GVariant *Tvar = ncm_matrix_peek_variant (T);
  
  NcmReparamLinear *relin = g_object_new (NCM_TYPE_REPARAM_LINEAR, 
                                          "length", size, 
                                          "matrix", Tvar, 
                                          "vector", vvar, 
                                          NULL);
  g_variant_unref (vvar);
  g_variant_unref (Tvar);

  return relin;
}

enum
{
  PROP_0,
  PROP_V,
  PROP_T
};

static gboolean
linear_old2new (NcmReparam *reparam, NcmModel *model, NcmVector *src, NcmVector *dest)
{
  NcmReparamLinear *relin = NCM_REPARAM_LINEAR (reparam);
  NCM_UNUSED (model);
  ncm_vector_memcpy (dest, relin->v);
  gsl_blas_dgemv (CblasNoTrans, 1.0, ncm_matrix_gsl (relin->T), ncm_vector_gsl(src), 1.0, ncm_vector_gsl (dest));
  return TRUE;
}

static gboolean
linear_new2old (NcmReparam *reparam, NcmModel *model, NcmVector *src, NcmVector *dest)
{
  NcmReparamLinear *relin = NCM_REPARAM_LINEAR (reparam);
  NCM_UNUSED (model);
  gsl_linalg_LU_solve (ncm_matrix_gsl (relin->T_LU), relin->p, ncm_vector_gsl(src), ncm_vector_gsl (dest));
  ncm_vector_sub (dest, relin->vp);
  return TRUE;
}

static gboolean
linear_jac (NcmReparam *reparam, NcmModel *model, NcmMatrix *jac)
{
  NCM_UNUSED (reparam);
  NCM_UNUSED (model);
  NCM_UNUSED (jac);
  g_assert_not_reached ();
  /*
   NcmReparamLinear *relin = NCM_REPARAM_LINEAR (reparam);
   gint i, j = 0;
   for (i = 0; i < ncm_model_len (model); i++)
   {
     if (NCM_FIT_PARAMS_IS_FREE(pt, i))
     {
       NcmVector *col = ncm_matrix_get_col (relin->T, i);
       ncm_matrix_set_col (jac, j++, col);
       ncm_vector_free (col);
}
}
*/
  return TRUE;
}

static void
_ncm_reparam_linear_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_reparam_linear_parent_class)->constructed (object);
  {
    NcmReparamLinear *relin = NCM_REPARAM_LINEAR (object);
    NcmReparam *reparam = NCM_REPARAM (relin);

    g_assert ((ncm_matrix_nrows (relin->T) == ncm_matrix_ncols (relin->T)) &&
              (ncm_matrix_nrows(relin->T) == reparam->length) &&
              (ncm_vector_len(relin->v) == reparam->length));

    relin->vp = ncm_vector_new (reparam->length);
    relin->T_LU = ncm_matrix_new (reparam->length, reparam->length);
    relin->p = gsl_permutation_alloc (reparam->length);

    ncm_matrix_memcpy (relin->T_LU, relin->T);
    gsl_linalg_LU_decomp (ncm_matrix_gsl (relin->T_LU), relin->p, &relin->signum);
    gsl_linalg_LU_solve (ncm_matrix_gsl (relin->T_LU), relin->p, ncm_vector_gsl(relin->v), ncm_vector_gsl(relin->vp));
  }
}

static void
_ncm_reparam_linear_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmReparamLinear *relin = NCM_REPARAM_LINEAR (object);

  g_return_if_fail (NCM_IS_REPARAM_LINEAR (object));

  switch (prop_id)
  {
    case PROP_V:
      g_value_take_variant (value, ncm_vector_get_variant (relin->v));
      break;
    case PROP_T:
      g_value_take_variant (value, ncm_matrix_get_variant (relin->T));
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
      relin->v = ncm_vector_new_variant (g_value_get_variant (value));
      break;
    case PROP_T:
      relin->T = ncm_matrix_new_variant (g_value_get_variant (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

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
ncm_reparam_linear_dispose (GObject *object)
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
ncm_reparam_linear_finalize (GObject *object)
{
  NcmReparamLinear *relin = NCM_REPARAM_LINEAR (object);

  gsl_permutation_free (relin->p);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_reparam_linear_parent_class)->finalize (object);
}

static void
_ncm_reparam_linear_copyto (NcmReparam *reparam, NcmReparam *reparam_dest)
{
  /* save old length since chain up copy will change it */
  guint old_length = reparam_dest->length;

  /* Chain up : start */
  NCM_REPARAM_CLASS (ncm_reparam_linear_parent_class)->copyto (reparam, reparam_dest);

  {
    NcmReparamLinear *relin = NCM_REPARAM_LINEAR (reparam);
    NcmReparamLinear *relin_dest = NCM_REPARAM_LINEAR (reparam_dest);
    if (old_length != reparam->length)
    {
      ncm_matrix_free (relin->T_LU);
      ncm_matrix_free (relin->T);
      ncm_vector_free (relin->v);
      ncm_vector_free (relin->vp);
      gsl_permutation_free (relin->p);

      relin->T_LU = ncm_matrix_new (reparam->length, reparam->length);
      relin->T = ncm_matrix_new (reparam->length, reparam->length);
      relin->v = ncm_vector_new (reparam->length);
      relin->vp = ncm_vector_new (reparam->length);
      relin->p = gsl_permutation_alloc (reparam->length);
    }

    ncm_matrix_memcpy (relin_dest->T, relin->T);
    ncm_matrix_memcpy (relin_dest->T_LU, relin->T_LU);
    ncm_vector_memcpy (relin_dest->v, relin->v);
    ncm_vector_memcpy (relin_dest->vp, relin->vp);
    gsl_permutation_memcpy (relin_dest->p, relin->p);
    relin_dest->signum = relin->signum;
  }
}

static NcmReparam *
_ncm_reparam_linear_copy (NcmReparam *reparam)
{
  NcmReparamLinear *relin = NCM_REPARAM_LINEAR (reparam);
  NcmReparamLinear *relin_new = ncm_reparam_linear_new (reparam->length, ncm_matrix_dup (relin->T), ncm_vector_dup (relin->v));
  NCM_REPARAM_CLASS (ncm_reparam_linear_parent_class)->copyto (reparam, NCM_REPARAM (relin_new));
  return NCM_REPARAM (relin_new);
}

static void
ncm_reparam_linear_class_init (NcmReparamLinearClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmReparamClass *parent_class = NCM_REPARAM_CLASS (klass);
  parent_class->copyto = &_ncm_reparam_linear_copyto;
  parent_class->copy = &_ncm_reparam_linear_copy;
  parent_class->old2new = &linear_old2new;
  parent_class->new2old = &linear_new2old;
  parent_class->jac = &linear_jac;

  object_class->constructed  = &_ncm_reparam_linear_constructed;
  object_class->set_property = &_ncm_reparam_linear_set_property;
  object_class->get_property = &_ncm_reparam_linear_get_property;

  /**
   * NcmReparamLinear:vector:
   *
   * The vector $\vec{v}$.
   * 
   */
  g_object_class_install_property (object_class,
                                   PROP_V,
                                   g_param_spec_variant ("vector",
                                                        NULL,
                                                        "Vector of parameters shift", 
                                                         G_VARIANT_TYPE ("ad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmReparamLinear:matrix:
   *
   * The matrix $M$.
   * 
   */
  g_object_class_install_property (object_class,
                                   PROP_T,
                                   g_param_spec_variant ("matrix",
                                                         NULL,
                                                         "Matrix of parameters mixing",
                                                         G_VARIANT_TYPE ("aad"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  object_class->dispose = ncm_reparam_linear_dispose;
  object_class->finalize = ncm_reparam_linear_finalize;
}
