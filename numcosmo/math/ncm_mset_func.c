/***************************************************************************
 *            ncm_mset_func.c
 *
 *  Wed June 06 15:32:21 2012
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
 * NcmMSetFunc:
 *
 * Abstract class for arbitrary MSet functions.
 *
 * This abstract class provides a framework for functions that operate on any model
 * in a #NcmMSet and additional extra variables. It establishes methods to
 * inquire about the function's expectations and characteristics.
 *
 * The functions implemented by subclasses may depend on any model specified by
 * #NcmMSet and may incorporate extra variables. The method
 * `ncm_mset_func_get_nvar()` retrieves the count of extra variables expected
 * by the function, and `ncm_mset_func_get_dim()` returns the number of values
 * returned by the function.
 *
 * Functions can be categorized as scalar or vectorial. A scalar function returns
 * a single value, while a vectorial function returns an array of values. The
 * method `ncm_mset_func_is_scalar()` returns %TRUE if the function is scalar.
 *
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mset_func.h"
#include "math/ncm_util.h"

enum
{
  PROP_0,
  PROP_NVAR,
  PROP_DIM,
  PROP_EVAL_X,
};

typedef struct _NcmMSetFuncPrivate
{
  /*< private >*/
  GObject parent_instance;
  guint nvar;
  guint dim;
  NcmVector *eval_x;
  gchar *name;
  gchar *symbol;
  gchar *ns;
  gchar *desc;
  gchar *uname;
  gchar *usymbol;
  NcmDiff *diff;
} NcmMSetFuncPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmMSetFunc, ncm_mset_func, G_TYPE_OBJECT)

static void
ncm_mset_func_init (NcmMSetFunc *func)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  self->nvar    = 0;
  self->dim     = 0;
  self->eval_x  = NULL;
  self->name    = NULL;
  self->symbol  = NULL;
  self->ns      = NULL;
  self->desc    = NULL;
  self->uname   = NULL;
  self->usymbol = NULL;
  self->diff    = ncm_diff_new ();
}

static void
_ncm_mset_func_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetFunc *func               = NCM_MSET_FUNC (object);
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  g_return_if_fail (NCM_IS_MSET_FUNC (object));

  switch (prop_id)
  {
    case PROP_NVAR:
      self->nvar = g_value_get_uint (value);
      break;
    case PROP_DIM:
      self->dim = g_value_get_uint (value);
      break;
    case PROP_EVAL_X:
      ncm_vector_clear (&self->eval_x);
      self->eval_x = g_value_dup_object (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_func_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetFunc *func               = NCM_MSET_FUNC (object);
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  g_return_if_fail (NCM_IS_MSET_FUNC (object));

  switch (prop_id)
  {
    case PROP_NVAR:
      g_value_set_uint (value, ncm_mset_func_get_nvar (func));
      break;
    case PROP_DIM:
      g_value_set_uint (value, ncm_mset_func_get_dim (func));
      break;
    case PROP_EVAL_X:
      g_value_set_object (value, self->eval_x);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_func_dispose (GObject *object)
{
  NcmMSetFunc *func               = NCM_MSET_FUNC (object);
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  ncm_vector_clear (&self->eval_x);
  ncm_diff_clear (&self->diff);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_func_parent_class)->dispose (object);
}

static void
_ncm_mset_func_finalize (GObject *object)
{
  NcmMSetFunc *func               = NCM_MSET_FUNC (object);
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  g_clear_pointer (&self->name,    g_free);
  g_clear_pointer (&self->symbol,  g_free);
  g_clear_pointer (&self->ns,      g_free);
  g_clear_pointer (&self->desc,    g_free);

  g_clear_pointer (&self->uname,   g_free);
  g_clear_pointer (&self->usymbol, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_func_parent_class)->finalize (object);
}

static void _ncm_mset_func_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res);

static void
ncm_mset_func_class_init (NcmMSetFuncClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_mset_func_set_property;
  object_class->get_property = &_ncm_mset_func_get_property;
  object_class->dispose      = &_ncm_mset_func_dispose;
  object_class->finalize     = &_ncm_mset_func_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NVAR,
                                   g_param_spec_uint ("nvariables",
                                                      NULL,
                                                      "Number of variables",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DIM,
                                   g_param_spec_uint ("dimension",
                                                      NULL,
                                                      "Function dimension",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_EVAL_X,
                                   g_param_spec_object ("eval-x",
                                                        NULL,
                                                        "Evaluation point x",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->eval = &_ncm_mset_func_eval;
}

static void
_ncm_mset_func_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  g_error ("_ncm_mset_func_eval: no eval function implemented.");
}

/**
 * ncm_mset_func_ref:
 * @func: a #NcmMSetFunc.
 *
 * Increases the reference count of @func by one.
 *
 * Returns: (transfer full): @func.
 */
NcmMSetFunc *
ncm_mset_func_ref (NcmMSetFunc *func)
{
  return g_object_ref (func);
}

/**
 * ncm_mset_func_free:
 * @func: a #NcmMSetFunc.
 *
 * Decreases the reference count of @func by one. If the reference count
 * reaches zero, @func is freed.
 *
 */
void
ncm_mset_func_free (NcmMSetFunc *func)
{
  g_object_unref (func);
}

/**
 * ncm_mset_func_clear:
 * @func: a #NcmMSetFunc.
 *
 * If *@func is not %NULL, decreases the reference count of @func by one
 * and sets *@func to %NULL.
 *
 */
void
ncm_mset_func_clear (NcmMSetFunc **func)
{
  g_clear_object (func);
}

/**
 * ncm_mset_func_array_new:
 *
 * Creates a new #GPtrArray to hold #NcmMSetFunc pointers.
 *
 * Returns: (element-type NcmMSetFunc) (transfer full): the new #GPtrArray.
 */
GPtrArray *
ncm_mset_func_array_new (void)
{
  return g_ptr_array_new_with_free_func ((GDestroyNotify) & ncm_mset_func_free);
}

/**
 * ncm_mset_func_eval: (virtual eval)
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x: (array) (element-type double): function arguments
 * @res: (array) (element-type double): function values
 *
 * Evaluate the function @func at @x and store the result in @res.
 *
 */
void
ncm_mset_func_eval (NcmMSetFunc *func, NcmMSet *mset, gdouble *x, gdouble *res)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  if (self->eval_x != NULL)
  {
    if (x != NULL)
      g_warning ("ncm_mset_func_eval: function called with arguments while an eval x was already used, ignoring argument.");

    NCM_MSET_FUNC_GET_CLASS (func)->eval (func, mset, ncm_vector_data (self->eval_x), res);
  }
  else
  {
    NCM_MSET_FUNC_GET_CLASS (func)->eval (func, mset, x, res);
  }
}

/**
 * ncm_mset_func_eval_nvar:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x: function arguments
 *
 * Evaluate the function @func at @x and return the result. This
 * function is only valid if @func is a scalar function.
 *
 * Returns: function value.
 */
gdouble
ncm_mset_func_eval_nvar (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x)
{
  gdouble res;

  NCM_MSET_FUNC_GET_CLASS (func)->eval (func, mset, x, &res);

  return res;
}

/**
 * ncm_mset_func_eval0:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 *
 * Evaluate the function @func and return the result. This
 * function is only valid if @func is a scalar function.
 * The function's arguments are either none, if the function is constant,
 * or the arguments passed to ncm_mset_func_set_eval_x().
 *
 * Returns: function value.
 */
gdouble
ncm_mset_func_eval0 (NcmMSetFunc *func, NcmMSet *mset)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);
  gdouble res;

  NCM_MSET_FUNC_GET_CLASS (func)->eval (func, mset, (self->eval_x != NULL) ? ncm_vector_data (self->eval_x) : NULL, &res);

  return res;
}

/**
 * ncm_mset_func_eval1:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x: function argument
 *
 * Evaluate the function @func at @x and return the result. This
 * function is only valid if @func is a scalar function.
 *
 * Returns: function value.
 */
gdouble
ncm_mset_func_eval1 (NcmMSetFunc *func, NcmMSet *mset, const gdouble x)
{
  gdouble res;

  NCM_MSET_FUNC_GET_CLASS (func)->eval (func, mset, &x, &res);

  return res;
}

/**
 * ncm_mset_func_eval_vector:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x_v: function arguments in a #NcmVector
 * @res_v: a #NcmVector to store the function values
 *
 * Compute the function @func at @x_v and store the result in @res_v. This function is
 * only valid if @func is a vectorial function.
 *
 */
void
ncm_mset_func_eval_vector (NcmMSetFunc *func, NcmMSet *mset, NcmVector *x_v, NcmVector *res_v)
{
  guint i;

  for (i = 0; i < ncm_vector_len (x_v); i++)
  {
    NCM_MSET_FUNC_GET_CLASS (func)->eval (func, mset, ncm_vector_ptr (x_v, i), ncm_vector_ptr (res_v, i));
  }
}

/**
 * ncm_mset_func_set_eval_x:
 * @func: a #NcmMSetFunc
 * @x: (in) (array length=len): function arguments
 * @len: length of @x
 *
 * Sets the function's arguments to @x. Once this function is called, the
 * function's arguments are fixed and @func becomes a constant function.
 *
 */
void
ncm_mset_func_set_eval_x (NcmMSetFunc *func, gdouble *x, guint len)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  ncm_vector_clear (&self->eval_x);

  g_assert_cmpuint (self->nvar, ==, len);

  self->eval_x = ncm_vector_new_data_dup (x, self->nvar, 1);

  ncm_mset_func_peek_name (func);
  ncm_mset_func_peek_desc (func);
  ncm_mset_func_peek_symbol (func);

  {
    GString *args_s = g_string_new ("(");
    gchar *args;
    guint i;

    for (i = 0; i < len; i++)
      g_string_append_printf (args_s, "%.15g", x[i]);

    g_string_append (args_s, ")");

    args = g_string_free (args_s, FALSE);

    g_clear_pointer (&self->usymbol, g_free);
    g_clear_pointer (&self->uname,   g_free);

    self->usymbol = g_strjoin (NULL, self->symbol, args, NULL);

    {
      GRegex *reg  = g_regex_new ("[^0-9]+", 0, 0, NULL);
      gchar *nargs = g_regex_replace_literal (reg, args, -1, 0, "_", 0, NULL);

      self->uname = g_strjoin (NULL, self->name, nargs, NULL);

      g_regex_unref (reg);
    }

    g_free (args);
  }
}

/**
 * ncm_mset_func_is_scalar:
 * @func: a #NcmMSetFunc
 *
 * Checks if @func is a scalar function.
 *
 * Returns: %TRUE if @func is scalar.
 */
gboolean
ncm_mset_func_is_scalar (NcmMSetFunc *func)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  return (self->dim == 1);
}

/**
 * ncm_mset_func_is_vector:
 * @func: a #NcmMSetFunc
 * @dim: function dimension
 *
 * Checks if @func is a vectorial function with dimension @dim.
 *
 * Returns: %TRUE if @func is vectorial with dimension @dim.
 */
gboolean
ncm_mset_func_is_vector (NcmMSetFunc *func, guint dim)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  return (self->dim == dim);
}

/**
 * ncm_mset_func_is_const:
 * @func: a #NcmMSetFunc
 *
 * Checks if @func is a constant function.
 *
 * Returns: %TRUE if @func is constant.
 */
gboolean
ncm_mset_func_is_const (NcmMSetFunc *func)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  return ((self->nvar == 0) || (self->eval_x != NULL));
}

/**
 * ncm_mset_func_has_nvar:
 * @func: a #NcmMSetFunc
 * @nvar: number of variables
 *
 * Checks if @func expects @nvar extra variables.
 *
 * Returns: %TRUE if @func expects @nvar extra variables.
 */
gboolean
ncm_mset_func_has_nvar (NcmMSetFunc *func, guint nvar)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  return (self->nvar == nvar);
}

typedef struct __ncm_mset_func_numdiff_fparams_1
{
  NcmMSetFunc *func;
  NcmMSet *mset;
  const gdouble *x;
} _ncm_mset_func_numdiff_fparams_1;

static gdouble
_mset_func_numdiff_fparams_1_val (NcmVector *x_v, gpointer userdata)
{
  _ncm_mset_func_numdiff_fparams_1 *nd = (_ncm_mset_func_numdiff_fparams_1 *) userdata;

  ncm_mset_fparams_set_vector (nd->mset, x_v);

  return ncm_mset_func_eval_nvar (nd->func, nd->mset, nd->x);
}

/**
 * ncm_mset_func_numdiff_fparams:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x: (array): function arguments
 * @out: (out) (transfer full): function gradient
 *
 * Computes the gradient of @func at @x and stores the result in @out.
 * This function is only valid if @func is a scalar function.
 *
 */
void
ncm_mset_func_numdiff_fparams (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, NcmVector **out)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);
  const guint fparam_len          = ncm_mset_fparam_len (mset);
  GArray *x_a                     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *x_v                  = NULL;
  GArray *grad_a                  = NULL;
  _ncm_mset_func_numdiff_fparams_1 nd;

  g_array_set_size (x_a, fparam_len);
  x_v = ncm_vector_new_array (x_a);
  ncm_mset_fparams_get_vector (mset, x_v);

  nd.mset = mset;
  nd.func = func;
  nd.x    = x;

  grad_a = ncm_diff_rf_d1_N_to_1 (self->diff, x_a, _mset_func_numdiff_fparams_1_val, &nd, NULL);

  ncm_mset_fparams_set_vector (mset, x_v);

  if (*out == NULL)
  {
    *out = ncm_vector_new_array (grad_a);
  }
  else
  {
    g_assert_cmpuint (fparam_len, ==, ncm_vector_len (*out));
    ncm_vector_set_array (*out, grad_a);
  }

  g_array_unref (x_a);
  g_array_unref (grad_a);
}

/**
 * ncm_mset_func_get_nvar:
 * @func: a #NcmMSetFunc
 *
 * Gets the number of variables of @func.
 *
 * Returns: number of variables expected by @func.
 */
guint
ncm_mset_func_get_nvar (NcmMSetFunc *func)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  return self->nvar;
}

/**
 * ncm_mset_func_get_dim:
 * @func: a #NcmMSetFunc
 *
 * Gets the dimension of @func.
 *
 * Returns: number values returned by @func.
 */
guint
ncm_mset_func_get_dim (NcmMSetFunc *func)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  return self->dim;
}

/**
 * ncm_mset_func_set_meta:
 * @func: a #NcmMSetFunc
 * @name: function name
 * @symbol: function symbol
 * @ns: function namespace
 * @desc: function description
 * @nvar: number of variables
 * @dim: function dimension
 *
 * Sets the function's metadata. This function is called by subclasses'
 * to set the function's metadata. It should not be called by users.
 *
 */
void
ncm_mset_func_set_meta (NcmMSetFunc *func, const gchar *name, const gchar *symbol, const gchar *ns, const gchar *desc, const guint nvar, const guint dim)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  g_clear_pointer (&self->name,   g_free);
  g_clear_pointer (&self->symbol, g_free);
  g_clear_pointer (&self->ns,     g_free);
  g_clear_pointer (&self->desc,   g_free);

  if (name != NULL)
    self->name = g_strdup (name);

  if (symbol != NULL)
    self->symbol = g_strdup (symbol);

  if (ns != NULL)
    self->ns = g_strdup (ns);

  if (desc != NULL)
    self->desc = g_strdup (desc);

  self->nvar = nvar;
  self->dim  = dim;
}

/**
 * ncm_mset_func_peek_name:
 * @func: a #NcmMSetFunc
 *
 * Returns: (transfer none): @func name.
 */
const gchar *
ncm_mset_func_peek_name (NcmMSetFunc *func)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  if (self->ns == NULL)
    self->ns = g_strdup (g_type_name (G_OBJECT_TYPE (func)));

  if (self->name == NULL)
    self->name = g_strdup_printf ("%s:no-name", self->ns);

  return self->name;
}

/**
 * ncm_mset_func_peek_symbol:
 * @func: a #NcmMSetFunc
 *
 * Returns: (transfer none): @func symbol.
 */
const gchar *
ncm_mset_func_peek_symbol (NcmMSetFunc *func)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  if (self->ns == NULL)
    self->ns = g_strdup (g_type_name (G_OBJECT_TYPE (func)));

  if (self->symbol == NULL)
    self->symbol = g_strdup_printf ("%s:no-symbol", self->ns);

  return self->symbol;
}

/**
 * ncm_mset_func_peek_ns:
 * @func: a #NcmMSetFunc
 *
 * Returns: (transfer none): @func ns.
 */
const gchar *
ncm_mset_func_peek_ns (NcmMSetFunc *func)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  if (self->ns == NULL)
    self->ns = g_strdup (g_type_name (G_OBJECT_TYPE (func)));

  return self->ns;
}

/**
 * ncm_mset_func_peek_desc:
 * @func: a #NcmMSetFunc
 *
 * Returns: (transfer none): @func desc.
 */
const gchar *
ncm_mset_func_peek_desc (NcmMSetFunc *func)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  if (self->ns == NULL)
    self->ns = g_strdup (g_type_name (G_OBJECT_TYPE (func)));

  if (self->desc == NULL)
    self->desc = g_strdup_printf ("%s:no-desc", self->ns);

  return self->desc;
}

/**
 * ncm_mset_func_peek_uname:
 * @func: a #NcmMSetFunc
 *
 * Peeks unique name.
 *
 * Returns: (transfer none): @func unique name.
 */
const gchar *
ncm_mset_func_peek_uname (NcmMSetFunc *func)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  if (self->uname != NULL)
    return self->uname;
  else
    return ncm_mset_func_peek_name (func);
}

/**
 * ncm_mset_func_peek_usymbol:
 * @func: a #NcmMSetFunc
 *
 * Peeks unique symbol.
 *
 * Returns: (transfer none): @func unique name.
 */
const gchar *
ncm_mset_func_peek_usymbol (NcmMSetFunc *func)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  if (self->usymbol != NULL)
    return self->usymbol;
  else
    return ncm_mset_func_peek_symbol (func);
}

