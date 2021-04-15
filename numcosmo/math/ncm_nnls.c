/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_nnls.c
 *
 *  Sun April 10 10:06:36 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_nnls.c
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_nnls
 * @title: NcmNNLS
 * @short_description: Non-negative linear least-squares
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_nnls.h"
#include "math/ncm_iset.h"
#include "math/ncm_lapack.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include "misc/libqp.h"
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmNNLSPrivate
{
  NcmNNLSUMethod umethod;
  gdouble reltol;
  guint nrows;
  guint ncols;
  guint uncols;
  NcmMatrix *A_QR;
  NcmMatrix *sub_A_QR;
  NcmMatrix *M;
  NcmMatrix *M_U;
  NcmVector *b;
  NcmVector *x_tmp;
  NcmVector *x_try;
  NcmVector *residuals;
  NcmVector *mgrad;
  NcmMatrix *sub_M_U;
  NcmVector *sub_x_tmp;
  NcmISet *Pset;
  NcmISet *Pset_try;
  NcmISet *invalid;
  gboolean LU_alloc;
  gboolean QR_alloc;
  GArray *ipiv;
  GArray *work;
};

enum
{
  PROP_0,
  PROP_UMETHOD,
  PROP_NROWS,
  PROP_NCOLS,
  PROP_RELTOL,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmNNLS, ncm_nnls, G_TYPE_OBJECT);

static void
ncm_nnls_init (NcmNNLS *nnls)
{
  NcmNNLSPrivate *const self = nnls->priv = ncm_nnls_get_instance_private (nnls);
  self->umethod   = NCM_NNLS_UMETHOD_LEN;
  self->reltol    = 0.0;
  self->nrows     = 0;
  self->ncols     = 0;
  self->uncols    = 0;
  self->A_QR      = NULL;
  self->sub_A_QR  = NULL;
  self->M         = NULL;
  self->M_U       = NULL;
  self->b         = NULL;
  self->x_tmp     = NULL;
  self->x_try     = NULL;
  self->residuals = NULL;
  self->mgrad     = NULL;
  self->sub_M_U   = NULL;
  self->sub_x_tmp = NULL;
  self->Pset      = NULL;
  self->Pset_try  = NULL;
  self->invalid   = NULL;
  self->LU_alloc  = FALSE;
  self->QR_alloc  = FALSE;
  self->ipiv      = g_array_new (FALSE, FALSE, sizeof (gint));
  self->work      = g_array_new (FALSE, FALSE, sizeof (gdouble));
}

static void _ncm_nnls_set_nrows (NcmNNLS *nnls, const guint nrows);
static void _ncm_nnls_set_ncols (NcmNNLS *nnls, const guint ncols);

static void
_ncm_nnls_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmNNLS *nnls = NCM_NNLS (object);
  g_return_if_fail (NCM_IS_NNLS (object));

  switch (prop_id)
  {
    case PROP_UMETHOD:
      ncm_nnls_set_umethod (nnls, g_value_get_enum (value));
      break;
    case PROP_RELTOL:
      ncm_nnls_set_reltol (nnls, g_value_get_double (value));
      break;
    case PROP_NROWS:
      _ncm_nnls_set_nrows (nnls, g_value_get_uint (value));
      break;
    case PROP_NCOLS:
      _ncm_nnls_set_ncols (nnls, g_value_get_uint (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_nnls_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmNNLS *nnls = NCM_NNLS (object);
  g_return_if_fail (NCM_IS_NNLS (object));

  switch (prop_id)
  {
    case PROP_UMETHOD:
      g_value_set_enum (value, ncm_nnls_get_umethod (nnls));
      break;
    case PROP_RELTOL:
      g_value_set_double (value, ncm_nnls_get_reltol (nnls));
      break;
    case PROP_NROWS:
      g_value_set_uint (value, ncm_nnls_get_nrows (nnls));
      break;
    case PROP_NCOLS:
      g_value_set_uint (value, ncm_nnls_get_ncols (nnls));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_nnls_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_nnls_parent_class)->constructed (object);
  {
    NcmNNLS *nnls = NCM_NNLS (object);
    NcmNNLSPrivate *const self = nnls->priv;

    self->b         = ncm_vector_new (self->ncols);
    self->x_tmp     = ncm_vector_new (self->ncols);
    self->x_try     = ncm_vector_new (self->ncols);
    self->residuals = ncm_vector_new (self->nrows);
    self->mgrad     = ncm_vector_new (self->ncols);
    self->Pset      = ncm_iset_new (self->ncols);
    self->Pset_try  = ncm_iset_new (self->ncols);
    self->invalid   = ncm_iset_new (self->ncols);

  }
}

static void
_ncm_nnls_dispose (GObject *object)
{
  NcmNNLS *nnls = NCM_NNLS (object);
  NcmNNLSPrivate * const self = nnls->priv;
  
  ncm_matrix_clear (&self->A_QR);
  ncm_matrix_clear (&self->M);
  ncm_matrix_clear (&self->M_U);
  ncm_vector_clear (&self->b);
  ncm_vector_clear (&self->x_tmp);
  ncm_vector_clear (&self->x_try);
  ncm_vector_clear (&self->residuals);
  ncm_vector_clear (&self->mgrad);

  ncm_matrix_clear (&self->sub_A_QR);
  ncm_matrix_clear (&self->sub_M_U);
  ncm_vector_clear (&self->sub_x_tmp);

  ncm_iset_clear (&self->Pset);
  ncm_iset_clear (&self->Pset_try);
  ncm_iset_clear (&self->invalid);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_nnls_parent_class)->dispose (object);
}

static void
_ncm_nnls_finalize (GObject *object)
{
  NcmNNLS *nnls = NCM_NNLS (object);
  NcmNNLSPrivate * const self = nnls->priv;

  g_array_unref (self->ipiv);
  g_array_unref (self->work);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_nnls_parent_class)->finalize (object);
}

static void
ncm_nnls_class_init (NcmNNLSClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_nnls_set_property;
  object_class->get_property = &_ncm_nnls_get_property;
  object_class->constructed  = &_ncm_nnls_constructed;
  object_class->dispose      = &_ncm_nnls_dispose;
  object_class->finalize     = &_ncm_nnls_finalize;

  g_object_class_install_property (object_class,
                                   PROP_UMETHOD,
                                   g_param_spec_enum ("umethod",
                                                      NULL,
                                                      "Unconstrained method",
                                                      NCM_TYPE_NNLS_UMETHOD, NCM_NNLS_UMETHOD_NORMAL,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                      NULL,
                                                      "Relative tolerance",
                                                      GSL_DBL_MIN, 1.0e-1, GSL_DBL_EPSILON,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NROWS,
                                   g_param_spec_uint ("nrows",
                                                      NULL,
                                                      "Number of rows",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NCOLS,
                                   g_param_spec_uint ("ncols",
                                                      NULL,
                                                      "Number of cols",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static void
_ncm_nnls_set_nrows (NcmNNLS *nnls, const guint nrows)
{
  NcmNNLSPrivate *const self = nnls->priv;
  self->nrows = nrows;
}

static void
_ncm_nnls_set_ncols (NcmNNLS *nnls, const guint ncols)
{
  NcmNNLSPrivate *const self = nnls->priv;
  self->ncols = ncols;
}

/**
 * ncm_nnls_new:
 * @nrows: number of rows
 * @ncols: number of columns
 * 
 * Creates a new #NcmNNLS object.
 * 
 * Returns: a new #NcmNNLS.
 */
NcmNNLS *
ncm_nnls_new (guint nrows, guint ncols)
{
  NcmNNLS *nnls = g_object_new (NCM_TYPE_NNLS,
                                "nrows", nrows,
                                "ncols", ncols,
                                NULL);
  return nnls;
}

/**
 * ncm_nnls_ref:
 * @nnls: a #NcmNNLS
 *
 * Increase the reference of @nnls by one.
 *
 * Returns: (transfer full): @nnls.
 */
NcmNNLS *
ncm_nnls_ref (NcmNNLS *nnls)
{
  return g_object_ref (nnls);
}

/**
 * ncm_nnls_free:
 * @nnls: a #NcmNNLS
 *
 * Decrease the reference count of @nnls by one.
 *
 */
void
ncm_nnls_free (NcmNNLS *nnls)
{
  g_object_unref (nnls);
}

/**
 * ncm_nnls_clear:
 * @nnls: a #NcmNNLS
 *
 * Decrease the reference count of @nnls by one, and sets the pointer *@nnls to
 * NULL.
 *
 */
void
ncm_nnls_clear (NcmNNLS **nnls)
{
  g_clear_object (nnls);
}

/**
 * ncm_nnls_set_umethod:
 * @nnls: a #NcmNNLS
 *
 * Sets which unconstrained least-squares method to use.
 *
 */
void
ncm_nnls_set_umethod (NcmNNLS *nnls, NcmNNLSUMethod umethod)
{
  NcmNNLSPrivate *const self = nnls->priv;

  self->umethod = umethod;
}

/**
 * ncm_nnls_get_umethod:
 * @nnls: a #NcmNNLS
 *
 * Gets the unconstrained least-squares method being used.
 *
 * Returns: a #NcmNNLSUMethod.
 */
NcmNNLSUMethod
ncm_nnls_get_umethod (NcmNNLS *nnls)
{
  NcmNNLSPrivate *const self = nnls->priv;

  return self->umethod;
}

/**
 * ncm_nnls_set_reltol:
 * @nnls: a #NcmNNLS
 * @reltol: a double
 *
 * Sets relative tolerance to @reltol.
 *
 */
void
ncm_nnls_set_reltol (NcmNNLS *nnls, const gdouble reltol)
{
  NcmNNLSPrivate *const self = nnls->priv;

  g_assert_cmpfloat (reltol, >=, GSL_DBL_MIN);
  g_assert_cmpfloat (reltol, <, 1.0);

  self->reltol = reltol;
}

/**
 * ncm_nnls_get_reltol:
 * @nnls: a #NcmNNLS
 *
 * Gets the relative tolerance being used.
 *
 * Returns: the current relative tolerance.
 */
gdouble
ncm_nnls_get_reltol (NcmNNLS *nnls)
{
  NcmNNLSPrivate *const self = nnls->priv;

  return self->reltol;
}

/**
 * ncm_nnls_get_nrows:
 * @nnls: a #NcmNNLS
 *
 * Returns: number of rows.
 */
guint
ncm_nnls_get_nrows (NcmNNLS *nnls)
{
  NcmNNLSPrivate *const self = nnls->priv;
  return self->nrows;
}

/**
 * ncm_nnls_get_ncols:
 * @nnls: a #NcmNNLS
 *
 * Returns: number of rows.
 */
guint
ncm_nnls_get_ncols (NcmNNLS *nnls)
{
  NcmNNLSPrivate *const self = nnls->priv;
  return self->ncols;
}

static void
_ncm_nnls_prepare_usys_QR (NcmNNLSPrivate *const self, NcmISet *Pset, NcmMatrix *A, NcmVector *x, NcmVector *f)
{
  self->uncols = ncm_iset_get_len (Pset);

  if (self->uncols == self->ncols)
  {
    ncm_matrix_memcpy (self->A_QR, A);

    ncm_matrix_clear (&self->sub_A_QR);
    ncm_vector_clear (&self->sub_x_tmp);

    self->sub_A_QR = ncm_matrix_ref (self->A_QR);
    self->sub_x_tmp = ncm_vector_get_subvector (self->residuals, 0, self->uncols);
  }
  else
  {
    ncm_matrix_clear (&self->sub_A_QR);
    ncm_vector_clear (&self->sub_x_tmp);

    self->sub_A_QR = ncm_iset_get_submatrix_cols (Pset, A, self->A_QR);
    self->sub_x_tmp = ncm_vector_get_subvector (self->residuals, 0, self->uncols);
  }

  ncm_vector_memcpy (self->residuals, f);
}


static void
_ncm_nnls_prepare_usys_normal (NcmNNLSPrivate *const self, NcmISet *Pset, NcmMatrix *A, NcmVector *x, NcmVector *f)
{
  self->uncols = ncm_iset_get_len (Pset);

  if (self->uncols == self->ncols)
  {
    ncm_matrix_memcpy (self->M_U, self->M);
    ncm_vector_memcpy (self->x_tmp, self->b);

    ncm_matrix_clear (&self->sub_M_U);
    ncm_vector_clear (&self->sub_x_tmp);

    self->sub_M_U   = ncm_matrix_ref (self->M_U);
    self->sub_x_tmp = ncm_vector_ref (self->x_tmp);
  }
  else
  {
    ncm_matrix_clear (&self->sub_M_U);
    ncm_vector_clear (&self->sub_x_tmp);

    self->sub_M_U   = ncm_iset_get_submatrix (Pset, self->M, self->M_U);
    self->sub_x_tmp = ncm_iset_get_subvector (Pset, self->b, self->x_tmp);
  }
}

static void
_ncm_nnls_solve_normal_LU (NcmNNLSPrivate *const self, NcmISet *Pset, NcmMatrix *A, NcmVector *x, NcmVector *f)
{
  gint ret, lwork;

  _ncm_nnls_prepare_usys_normal (self, Pset, A, x, f);

  if (self->work->len == 0)
    g_array_set_size (self->work, self->ncols);
  g_array_set_size (self->ipiv, self->ncols);

  lwork = -1;
  ret = ncm_lapack_dsysv ('U', self->uncols, 1,
                          ncm_matrix_data (self->sub_M_U), ncm_matrix_tda (self->sub_M_U),
                          &g_array_index (self->ipiv, gint, 0),
                          ncm_vector_data (self->sub_x_tmp), self->uncols,
                          &g_array_index (self->work, gdouble, 0), lwork);
  g_assert_cmpint (ret, ==, 0);

  lwork = g_array_index (self->work, gdouble, 0);
  if (lwork > self->work->len)
    g_array_set_size (self->work, lwork);

  ret = ncm_lapack_dsysv ('U', self->uncols, 1,
                          ncm_matrix_data (self->M_U), ncm_matrix_tda (self->M_U),
                          &g_array_index (self->ipiv, gint, 0),
                          ncm_vector_data (self->sub_x_tmp), self->uncols,
                          &g_array_index (self->work, gdouble, 0), lwork);
  g_assert_cmpint (ret, ==, 0);
}

static void
_ncm_nnls_solve_normal_QR (NcmNNLSPrivate *const self, NcmISet *Pset, NcmMatrix *A, NcmVector *x, NcmVector *f)
{
  gint ret, lwork;

  _ncm_nnls_prepare_usys_QR (self, Pset, A, x, f);

  if (self->work->len == 0)
    g_array_set_size (self->work, self->ncols);

  lwork = -1;
  ret = ncm_lapack_dgels ('N', self->nrows, self->uncols, 1,
                          ncm_matrix_data (self->sub_A_QR), ncm_matrix_tda (self->sub_A_QR),
                          ncm_vector_data (self->residuals), self->nrows,
                          &g_array_index (self->work, gdouble, 0), lwork);
  g_assert_cmpint (ret, ==, 0);

  lwork = g_array_index (self->work, gdouble, 0);
  if (lwork > self->work->len)
    g_array_set_size (self->work, lwork);

  ret = ncm_lapack_dgels ('N', self->nrows, self->uncols, 1,
                          ncm_matrix_data (self->sub_A_QR), ncm_matrix_tda (self->sub_A_QR),
                          ncm_vector_data (self->residuals), self->nrows,
                          &g_array_index (self->work, gdouble, 0), lwork);
  g_assert_cmpint (ret, ==, 0);
}

static void
_ncm_nnls_solve_normal_cholesky (NcmNNLSPrivate *const self, NcmISet *Pset, NcmMatrix *A, NcmVector *x, NcmVector *f)
{
  gint ret;

  _ncm_nnls_prepare_usys_normal (self, Pset, A, x, f);

  ret = ncm_matrix_cholesky_solve (self->sub_M_U, self->sub_x_tmp, 'U');
  if (ret > 0)
    _ncm_nnls_solve_normal_LU (self, Pset, A, x, f);
}

static void
_ncm_nnls_solve_unconstrained (NcmNNLSPrivate *const self, NcmISet *Pset, NcmMatrix *A, NcmVector *x, NcmVector *f)
{
  switch (self->umethod)
  {
    case NCM_NNLS_UMETHOD_NORMAL:
      _ncm_nnls_solve_normal_cholesky (self, Pset, A, x, f);
      break;
    case NCM_NNLS_UMETHOD_NORMAL_LU:
      _ncm_nnls_solve_normal_LU (self, Pset, A, x, f);
      break;
    case NCM_NNLS_UMETHOD_QR:
      _ncm_nnls_solve_normal_QR (self, Pset, A, x, f);
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}

static gdouble
_ncm_nnls_compute_residuals (NcmNNLSPrivate *const self, NcmMatrix *A, NcmVector *x, NcmVector *f, NcmVector *residuals)
{
  gdouble rnorm;

  ncm_vector_memcpy (residuals, f);
  ncm_matrix_update_vector (A, 'N', -1.0, x, 1.0, residuals);
  rnorm = ncm_vector_dnrm2 (residuals);

  return rnorm;
}

static void
_ncm_nnls_compute_mgrad (NcmNNLSPrivate *const self, NcmMatrix *A, NcmVector *x, NcmVector *f, NcmVector *residuals, NcmVector *mgrad)
{
  ncm_matrix_update_vector (A, 'T', 1.0, residuals, 0.0, mgrad);
}

static void
_ncm_nnls_solve_feasible (NcmNNLSPrivate *const self, NcmISet *Pset, NcmMatrix *A, NcmVector *x, NcmVector *f, guint max_remove)
{
  _ncm_nnls_solve_unconstrained (self, Pset, A, x, f);

  ncm_vector_set_zero (x);
  ncm_iset_set_subvector (Pset, x, self->sub_x_tmp);

  ncm_iset_get_subset_vec_lt (Pset, self->invalid, x, 0.0);

  while (ncm_iset_get_len (self->invalid))
  {
    ncm_iset_remove_smallest_subset (self->invalid, Pset, x, max_remove);

    _ncm_nnls_solve_unconstrained (self, Pset, A, x, f);

    ncm_vector_set_zero (x);
    ncm_iset_set_subvector (Pset, x, self->sub_x_tmp);

    ncm_iset_get_subset_vec_lt (Pset, self->invalid, x, 0.0);
  }
}

/**
 * ncm_nnls_solve:
 * @nnls: a #NcmNNLS
 * @A: a #NcmMatrix $A$
 * @x: a #NcmVector $\vec{x}$
 * @f: a #NcmVector $\vec{f}$
 *
 * Solves the system $A\vec{x} = \vec{f}$ for $\vec{x}$
 * imposing the non negativity constraint on $\vec{x}$,
 * i.e., $\vec{x} > 0$.
 *
 *
 * Returns: the Euclidean norm of the residuals.
 */
gdouble
ncm_nnls_solve (NcmNNLS *nnls, NcmMatrix *A, NcmVector *x, NcmVector *f)
{
  NcmNNLSPrivate *const self = nnls->priv;
  gdouble rnorm;

  g_assert_cmpuint (ncm_matrix_nrows (A), ==, self->nrows);
  g_assert_cmpuint (ncm_matrix_ncols (A), ==, self->ncols);

  switch (self->umethod)
  {
    case NCM_NNLS_UMETHOD_NORMAL:
    case NCM_NNLS_UMETHOD_NORMAL_LU:
    {
      if (!self->LU_alloc)
      {
        self->M   = ncm_matrix_new (self->ncols, self->ncols);
        self->M_U = ncm_matrix_new (self->ncols, self->ncols);
        self->LU_alloc = TRUE;
      }

      ncm_matrix_square_to_sym (A, 'T', 'U', self->M);
      ncm_matrix_update_vector (A, 'T', 1.0, f, 0.0, self->b);
      break;
    }
    case NCM_NNLS_UMETHOD_QR:
    {
      if (!self->QR_alloc)
      {
        self->A_QR     = ncm_matrix_new (self->nrows, self->ncols);
        self->QR_alloc = TRUE;
      }
      break;
    }
    default:
      g_assert_not_reached ();
      break;
  }

  ncm_iset_reset (self->Pset);

  /*ncm_iset_add (self->Pset, ncm_vector_get_max_index (self->b));*/
  ncm_iset_add_range (self->Pset, 0, self->ncols);

  _ncm_nnls_solve_feasible (self, self->Pset, A, x, f, self->ncols);
  rnorm = _ncm_nnls_compute_residuals (self, A, x, f, self->residuals);
  _ncm_nnls_compute_mgrad (self, A, x, f, self->residuals, self->mgrad);

  while (TRUE)
  {
    gdouble add_frac = 1.0;
    gboolean finish  = FALSE;
    gdouble lrnorm;
    guint added;

    while (TRUE)
    {
      add_frac *= 0.2;

      ncm_iset_copy (self->Pset, self->Pset_try);
      added = ncm_iset_add_largest_subset (self->Pset_try, self->mgrad, 0.0, add_frac);
      if (added == 0)
      {
        finish = TRUE;
        break;
      }

      _ncm_nnls_solve_feasible (self, self->Pset_try, A, self->x_try, f, added);
      lrnorm = _ncm_nnls_compute_residuals (self, A, self->x_try, f, self->residuals);

      if (lrnorm < rnorm * (1.0 - self->reltol))
      {
        ncm_iset_copy (self->Pset_try, self->Pset);
        ncm_vector_memcpy (x, self->x_try);
        rnorm = lrnorm;
        break;
      }

      if (added == 1)
      {
        finish = TRUE;
        break;
      }
    }

    if (finish)
      break;

    _ncm_nnls_compute_mgrad (self, A, x, f, self->residuals, self->mgrad);
  }

  return rnorm;
}

int nnls_c (double *a, const int *mda, const int *m, const int *n, double *b,
            double *x, double* rnorm, double* w, double* zz, int *index,
            int *mode);

/**
 * ncm_nnls_solve_LH:
 * @nnls: a #NcmNNLS
 * @A: a #NcmMatrix $A$
 * @x: a #NcmVector $\vec{x}$
 * @f: a #NcmVector $\vec{f}$
 *
 * Solves the system $A\vec{x} = \vec{f}$ for $\vec{x}$
 * imposing the non negativity constraint on $\vec{x}$,
 * i.e., $\vec{x} > 0$. This method solves the system
 * using the original code by Charles L. Lawson and
 * Richard J. Hanson translated to C using f2c.
 *
 * Returns: the Euclidean norm of the residuals.
 */
gdouble
ncm_nnls_solve_LH (NcmNNLS *nnls, NcmMatrix *A, NcmVector *x, NcmVector *f)
{
  NcmNNLSPrivate *const self = nnls->priv;
  gint nrows = self->nrows;
  gint ncols = self->ncols;
  gint mode = -31;
  gdouble rnorm;

  g_assert_cmpuint (ncm_matrix_nrows (A), ==, self->nrows);
  g_assert_cmpuint (ncm_matrix_ncols (A), ==, self->ncols);

  if (!self->QR_alloc)
  {
    self->A_QR     = ncm_matrix_new (self->nrows, self->ncols);
    self->QR_alloc = TRUE;
  }

  ncm_matrix_memcpy_to_colmajor (self->A_QR, A);
  ncm_vector_memcpy (self->residuals, f);

  if (self->work->len < self->nrows)
    g_array_set_size (self->work, self->nrows);
  if (self->ipiv->len < self->ncols)
    g_array_set_size (self->ipiv, self->ncols);

  nnls_c (ncm_matrix_data (A), &nrows, &nrows, &ncols, ncm_vector_data (self->residuals),
      ncm_vector_data (x), &rnorm, ncm_vector_data (self->x_tmp), &g_array_index (self->work, gdouble, 0),
      &g_array_index (self->ipiv, gint, 0), &mode);

  return rnorm;

}

void LowRankQP (gint *n, gint *m, gint *p, gint *method, gint *verbose, gint *niter,
                gdouble *Q, gdouble *c, gdouble *A, gdouble *b, gdouble *u, gdouble *alpha,
                gdouble *beta, gdouble *xi, gdouble *zeta);

/**
 * ncm_nnls_solve_lowrankqp:
 * @nnls: a #NcmNNLS
 * @A: a #NcmMatrix $A$
 * @x: a #NcmVector $\vec{x}$
 * @f: a #NcmVector $\vec{f}$
 *
 * Solves the system $A\vec{x} = \vec{f}$ for $\vec{x}$
 * imposing the non negativity constraint on $\vec{x}$,
 * i.e., $\vec{x} > 0$. This method solves the system
 * using the LowRankQP quadratic programming code.
 *
 * Returns: the Euclidean norm of the residuals.
 */
gdouble
ncm_nnls_solve_lowrankqp (NcmNNLS *nnls, NcmMatrix *A, NcmVector *x, NcmVector *f)
{
  NcmNNLSPrivate *const self = nnls->priv;
  gint nrows   = self->nrows;
  gint ncols   = self->ncols;
  gint nc      = 0;
  gint verbose = 0;
  gint maxiter = 4000;
  gint method;

  /* Setting upper-bound */
  ncm_vector_set_all (self->x_tmp, 1.0e20);

  g_assert_cmpuint (self->ncols, !=, self->nrows);
  g_assert_cmpuint (ncm_matrix_nrows (A), ==, self->nrows);
  g_assert_cmpuint (ncm_matrix_ncols (A), ==, self->ncols);

  switch (self->umethod)
  {
    case NCM_NNLS_UMETHOD_NORMAL:
      method = 2;
      break;
    case NCM_NNLS_UMETHOD_NORMAL_LU:
      method = 1;
      break;
    default:
      g_error ("ncm_nnls_solve_lowrankqp: LowRankQP only support UMETHOD_NORMAL and UMETHOD_NORMAL_QP.");
      method = 2;
      break;
  }

  if (!self->QR_alloc)
  {
    self->A_QR     = ncm_matrix_new (self->nrows, self->ncols);
    self->QR_alloc = TRUE;
  }
  if (self->work->len < self->nrows)
    g_array_set_size (self->work, self->nrows);

  ncm_matrix_memcpy (self->A_QR, A);
  ncm_matrix_update_vector (A, 'T', -1.0, f, 0.0, self->b);

  LowRankQP (&ncols, &nrows, &nc, &method, &verbose, &maxiter,
             ncm_matrix_data (self->A_QR), ncm_vector_data (self->b),
             NULL, NULL, ncm_vector_data (self->x_tmp),
             ncm_vector_data (x),
             ncm_vector_data (self->x_try),
             ncm_vector_data (self->residuals),
             &g_array_index (self->work, gdouble, 0));

  return _ncm_nnls_compute_residuals (self, A, x, f, self->residuals);
}

static void
print_state (libqp_state_T state)
{
  ncm_message ("niter %d QP % 22.15g QD % 22.15g %d\n", state.nIter, state.QP, state.QD, state.exitflag);
}

/**
 * ncm_nnls_solve_splx:
 * @nnls: a #NcmNNLS
 * @A: a #NcmMatrix $A$
 * @x: a #NcmVector $\vec{x}$
 * @f: a #NcmVector $\vec{f}$
 *
 * Solves the system $A\vec{x} = \vec{f}$ for $\vec{x}$
 * imposing the non negativity constraint on $\vec{x}$,
 * i.e., $\vec{x} > 0$. This method solves the system
 * using function libqp_splx_solver from [libqp](https://cmp.felk.cvut.cz/~xfrancv/libqp/html/).
 *
 * Returns: the Euclidean norm of the residuals.
 */
gdouble
ncm_nnls_solve_splx (NcmNNLS *nnls, NcmMatrix *A, NcmVector *x, NcmVector *f)
{
  NcmNNLSPrivate *const self = nnls->priv;
  GPtrArray *col   = g_ptr_array_new ();
  const gint ncols = self->ncols;
  uint32_t *II = g_new (uint32_t, ncols);
  gdouble bbb = 1.0;
  uint8_t S = 1;
  libqp_state_T res;
  gint i;

  if (!self->LU_alloc)
  {
    self->M   = ncm_matrix_new (ncols, ncols);
    self->M_U = ncm_matrix_new (ncols, ncols);
    self->LU_alloc = TRUE;
  }

  ncm_matrix_square_to_sym (A, 'T', 'U', self->M);
  ncm_matrix_update_vector (A, 'T', -1.0, f, 0.0, self->b);

  for (i = 0; i < ncols; i++)
    g_ptr_array_add (col, ncm_matrix_ptr (self->M, i, 0));

  ncm_vector_set_all (self->x_tmp, 1.0);
  ncm_vector_set_all (x, 0.0);

  for (i = 0; i < ncols; i++)
    II[i] = 1;

  res = libqp_splx_solver ((gdouble **)col->pdata,
                           ncm_vector_data (self->x_tmp),
                           ncm_vector_data (self->b),
                           &bbb,
                           II,
                           &S,
                           ncm_vector_data (x),
                           ncols, 1000,
                           0.0, self->reltol, GSL_NEGINF,
                           /*&print_state*/ NULL);

  g_ptr_array_unref (col);
  g_free (II);

  g_assert_cmpint (res.exitflag, >=, 0);
  if (FALSE)
    print_state (res);

  return _ncm_nnls_compute_residuals (self, A, x, f, self->residuals);
}

/**
 * ncm_nnls_solve_gmso:
 * @nnls: a #NcmNNLS
 * @A: a #NcmMatrix $A$
 * @x: a #NcmVector $\vec{x}$
 * @f: a #NcmVector $\vec{f}$
 *
 * Solves the system $A\vec{x} = \vec{f}$ for $\vec{x}$
 * imposing the non negativity constraint on $\vec{x}$,
 * i.e., $\vec{x} > 0$. This method solves the system
 * using function libqp_gmso_solver from [libqp](https://cmp.felk.cvut.cz/~xfrancv/libqp/html/).
 *
 * Returns: the Euclidean norm of the residuals.
 */
gdouble
ncm_nnls_solve_gsmo (NcmNNLS *nnls, NcmMatrix *A, NcmVector *x, NcmVector *f)
{
  NcmNNLSPrivate *const self = nnls->priv;
  GPtrArray *col   = g_ptr_array_new ();
  const gint ncols = self->ncols;
  uint32_t *II = g_new (uint32_t, ncols);
  libqp_state_T res;
  gint i;

  if (!self->LU_alloc)
  {
    self->M   = ncm_matrix_new (ncols, ncols);
    self->M_U = ncm_matrix_new (ncols, ncols);
    self->LU_alloc = TRUE;
  }

  ncm_matrix_square_to_sym (A, 'T', 'U', self->M);
  ncm_matrix_update_vector (A, 'T', -1.0, f, 0.0, self->b);

  for (i = 0; i < ncols; i++)
    g_ptr_array_add (col, ncm_matrix_ptr (self->M, i, 0));

  ncm_vector_set_all (self->x_tmp, 1.0);
  ncm_vector_set_all (self->x_try, 1.0);
  ncm_vector_set_all (self->mgrad, 0.0);
  ncm_vector_set_all (self->residuals, GSL_POSINF);
  ncm_vector_set_all (x, 1.0 / ncols);

  for (i = 0; i < ncols; i++)
    II[i] = 1;

  res = libqp_gsmo_solver ((gdouble **)col->pdata,
                           ncm_vector_data (self->x_tmp),
                           ncm_vector_data (self->b),
                           ncm_vector_data (self->x_try),
                           1.0,
                           ncm_vector_data (self->mgrad),
                           ncm_vector_data (self->residuals),
                           ncm_vector_data (x),
                           ncols, 1000,
                           self->reltol,
                           /*&print_state*/ NULL);

  g_ptr_array_unref (col);
  g_free (II);

  g_assert_cmpint (res.exitflag, >=, 0);
  if (FALSE)
    print_state (res);

  return _ncm_nnls_compute_residuals (self, A, x, f, self->residuals);
}
