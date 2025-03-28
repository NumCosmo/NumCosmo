/***************************************************************************
 *            ncm_mset_func1.c
 *
 *  Sun May 20 21:32:30 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_func1.c
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmMSetFunc1:
 *
 * Abstract class for arbitrary MSet functions - bindable version
 *
 * This class is an abstract class for arbitrary MSet functions it behaves exactly like
 * #NcmMSetFunc but its virtual function is bindable.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mset_func1.h"

typedef struct _NcmMSetFunc1Private
{
  gint placeholder;
} NcmMSetFunc1Private;

G_DEFINE_TYPE_WITH_PRIVATE (NcmMSetFunc1, ncm_mset_func1, NCM_TYPE_MSET_FUNC)

static void
ncm_mset_func1_init (NcmMSetFunc1 *f1)
{
  NcmMSetFunc1Private *self = ncm_mset_func1_get_instance_private (f1);

  self->placeholder = 0;
}

static void
_ncm_mset_func1_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_func1_parent_class)->finalize (object);
}

void _ncm_mset_func1_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res);

static void
ncm_mset_func1_class_init (NcmMSetFunc1Class *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcmMSetFuncClass *func_class = NCM_MSET_FUNC_CLASS (klass);

  object_class->finalize = &_ncm_mset_func1_finalize;

  func_class->eval = _ncm_mset_func1_eval;
}

void
_ncm_mset_func1_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  GArray *x_a      = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const guint nvar = ncm_mset_func_get_nvar (func);
  const guint dim  = ncm_mset_func_get_dim (func);
  GArray *res_a;

  g_array_append_vals (x_a, x, nvar);

  res_a = ncm_mset_func1_eval1 (NCM_MSET_FUNC1 (func), mset, x_a);

  g_assert_cmpint (res_a->len, ==, dim);
  memcpy (res, res_a->data, res_a->len * sizeof (gdouble));

  g_array_unref (x_a);
  g_array_unref (res_a);
}

/**
 * ncm_mset_func1_ref:
 * @f1: a #NcmMSetFunc1
 *
 * Increments the reference count of @f1 by one.
 *
 * Returns: (transfer full): @f1.
 */
NcmMSetFunc1 *
ncm_mset_func1_ref (NcmMSetFunc1 *f1)
{
  return g_object_ref (f1);
}

/**
 * ncm_mset_func1_free:
 * @f1: a #NcmMSetFunc1
 *
 * Decrements the reference count of @f1 by one. If the reference count
 * reaches zero, @f1 is freed.
 *
 */
void
ncm_mset_func1_free (NcmMSetFunc1 *f1)
{
  g_object_unref (f1);
}

/**
 * ncm_mset_func1_clear:
 * @f1: a #NcmMSetFunc1
 *
 * If *@f1 is non-%NULL, unrefs it and sets *@f1 to %NULL.
 *
 */
void
ncm_mset_func1_clear (NcmMSetFunc1 **f1)
{
  g_clear_object (f1);
}

/**
 * ncm_mset_func1_eval1: (virtual eval1)
 * @mset: a #NcmMSet
 * @f1: a #NcmMSetFunc1
 * @x: (array) (element-type double): function argument
 *
 * Evaluates the function at @x.
 *
 * Returns: (array) (element-type double) (transfer full): function result
 */
GArray *
ncm_mset_func1_eval1 (NcmMSetFunc1 *f1, NcmMSet *mset, GArray *x)
{
  return NCM_MSET_FUNC1_GET_CLASS (f1)->eval1 (f1, mset, x);
}

