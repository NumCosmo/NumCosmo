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
 * Abstract class for #NcmMSetFunc functions written in a language binding.
 *
 * The eval virtual function of #NcmMSetFunc takes and fills plain C arrays whose
 * lengths language bindings cannot know, so it cannot be overridden from them.
 * This class implements it through the eval1 virtual function, which takes the
 * arguments and returns the values as #GArray. In Python a subclass overrides
 * `do_eval1` and sets the number of variables and values with the "nvariables"
 * and "dimension" properties.
 *
 * Evaluate the function with the #NcmMSetFunc methods, for instance
 * ncm_mset_func_eval_array() or ncm_mset_func_eval0().
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/model/ncm_mset_func1.h"

G_DEFINE_ABSTRACT_TYPE (NcmMSetFunc1, ncm_mset_func1, NCM_TYPE_MSET_FUNC)

static void
ncm_mset_func1_init (NcmMSetFunc1 *f1)
{
}

static void _ncm_mset_func1_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res);
static GArray *_ncm_mset_func1_eval1 (NcmMSetFunc1 *f1, NcmMSet *mset, GArray *x);

static void
ncm_mset_func1_class_init (NcmMSetFunc1Class *klass)
{
  NcmMSetFuncClass *func_class = NCM_MSET_FUNC_CLASS (klass);

  func_class->eval = &_ncm_mset_func1_eval;
  klass->eval1     = &_ncm_mset_func1_eval1;
}

static void
_ncm_mset_func1_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  NcmMSetFunc1 *f1 = NCM_MSET_FUNC1 (func);
  const guint nvar = ncm_mset_func_get_nvar (func);
  const guint dim  = ncm_mset_func_get_dim (func);
  GArray *x_a      = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), nvar);
  GArray *res_a;

  g_array_append_vals (x_a, x, nvar);

  res_a = NCM_MSET_FUNC1_GET_CLASS (f1)->eval1 (f1, mset, x_a);

  if ((res_a == NULL) || (res_a->len != dim))
    g_error ("_ncm_mset_func1_eval: function `%s' has dimension %u, but eval1 returned %u value(s).",
             ncm_mset_func_peek_name (func), dim, (res_a != NULL) ? res_a->len : 0);

  memcpy (res, res_a->data, dim * sizeof (gdouble));

  g_array_unref (x_a);
  g_array_unref (res_a);
}

static GArray *
_ncm_mset_func1_eval1 (NcmMSetFunc1 *f1, NcmMSet *mset, GArray *x)
{
  g_error ("_ncm_mset_func1_eval1: no eval1 function implemented.");

  return NULL;
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

