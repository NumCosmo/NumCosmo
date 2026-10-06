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
 * Abstract class for functions of the models in a #NcmMSet.
 *
 * A function may also take variables besides the models, for instance a redshift,
 * and may return several values. ncm_mset_func_get_nvar() gives the number of
 * variables and ncm_mset_func_get_dim() the number of values; a function with one
 * value is scalar (see ncm_mset_func_is_scalar()).
 *
 * ncm_mset_func_set_eval_x() binds the variables to an evaluation point. The
 * function then becomes constant and its unique name and symbol, which label its
 * column in a #NcmMSetCatalog, encode that point. Every evaluation follows one
 * rule: explicit arguments always win, and %NULL arguments mean the evaluation
 * point. ncm_mset_func_eval1() and ncm_mset_func_eval_vector() always take their
 * arguments.
 *
 * ncm_mset_func_eval() and ncm_mset_func_eval_array() return every value; the
 * scalar evaluations ncm_mset_func_eval0(), ncm_mset_func_eval_nvar() and
 * ncm_mset_func_eval1() abort for a function with more than one value.
 * Subclasses implement the eval virtual function; in language bindings, subclass
 * #NcmMSetFunc1 instead.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/model/ncm_mset_func.h"
#include "ncm/core/ncm_util.h"

enum
{
  PROP_0,
  PROP_NVAR,
  PROP_DIM,
  PROP_EVAL_X,
};

typedef struct _NcmMSetFuncPrivate
{
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

static gchar *_ncm_mset_func_shortest_double (const gdouble x);

/*
 * _ncm_mset_func_update_unames:
 * @func: a #NcmMSetFunc
 *
 * Rebuilds the unique name and symbol from the evaluation point @eval_x. The unique
 * name is the column name of the function in a #NcmMSetCatalog, so two evaluation
 * points must never give the same name. Otherwise a function evaluated over a
 * redshift grid, all named "wDE_z", would produce repeated column names.
 *
 * Each component is written in the shortest form that reads back to the same double.
 * The symbol lists the components separated by commas, as in f(1.5,-2e-05). The
 * name maps each character one to one, '-' to 'm' and '.' to 'p', drops the '+' of
 * the exponent and separates the components by '_', as in f_1p5_m2em05.
 *
 * Both ncm_mset_func_set_eval_x() and the "eval-x" property, restored on
 * deserialization, call this function.
 */
static void
_ncm_mset_func_update_unames (NcmMSetFunc *func)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  g_clear_pointer (&self->usymbol, g_free);
  g_clear_pointer (&self->uname,   g_free);

  if (self->eval_x == NULL)
    return;

  {
    const guint len    = ncm_vector_len (self->eval_x);
    GString *usymbol_s = g_string_new (ncm_mset_func_peek_symbol (func));
    GString *uname_s   = g_string_new (ncm_mset_func_peek_name (func));
    guint i;

    g_string_append_c (usymbol_s, '(');

    for (i = 0; i < len; i++)
    {
      gchar *x_s = _ncm_mset_func_shortest_double (ncm_vector_get (self->eval_x, i));
      gchar *c;

      if (i > 0)
        g_string_append_c (usymbol_s, ',');

      g_string_append (usymbol_s, x_s);
      g_string_append_c (uname_s, '_');

      for (c = x_s; *c != '\0'; c++)
      {
        switch (*c)
        {
          case '-':
            g_string_append_c (uname_s, 'm');
            break;
          case '.':
            g_string_append_c (uname_s, 'p');
            break;
          case '+':
            break;
          default:
            g_string_append_c (uname_s, *c);
            break;
        }
      }

      g_free (x_s);
    }

    g_string_append_c (usymbol_s, ')');

    self->usymbol = g_string_free (usymbol_s, FALSE);
    self->uname   = g_string_free (uname_s, FALSE);
  }
}

/*
 * _ncm_mset_func_shortest_double:
 * @x: a double
 *
 * Writes @x with the fewest significant digits, from 15 to 17, that read back to
 * the same double. The output does not depend on the locale.
 *
 * Returns: (transfer full): the string.
 */
static gchar *
_ncm_mset_func_shortest_double (const gdouble x)
{
  static const gchar *formats[] = {"%.15g", "%.16g"};
  gchar buf[G_ASCII_DTOSTR_BUF_SIZE];
  guint i;

  for (i = 0; i < G_N_ELEMENTS (formats); i++)
  {
    g_ascii_formatd (buf, G_ASCII_DTOSTR_BUF_SIZE, formats[i], x);

    if (g_ascii_strtod (buf, NULL) == x)
      return g_strdup (buf);
  }

  g_ascii_formatd (buf, G_ASCII_DTOSTR_BUF_SIZE, "%.17g", x);

  return g_strdup (buf);
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
    {
      NcmVector *eval_x = g_value_get_object (value);

      if (eval_x == NULL)
      {
        ncm_vector_clear (&self->eval_x);
        _ncm_mset_func_update_unames (func);
      }
      else
      {
        /* The copy is contiguous, @eval_x may have a stride. */
        NcmVector *eval_x_c = ncm_vector_dup (eval_x);

        ncm_mset_func_set_eval_x (func, ncm_vector_data (eval_x_c), ncm_vector_len (eval_x_c));
        ncm_vector_free (eval_x_c);
      }

      break;
    }
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

  /**
   * NcmMSetFunc:nvariables:
   *
   * The number of variables the function takes besides the models.
   */
  g_object_class_install_property (object_class,
                                   PROP_NVAR,
                                   g_param_spec_uint ("nvariables",
                                                      NULL,
                                                      "Number of variables",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetFunc:dimension:
   *
   * The number of values the function returns.
   */
  g_object_class_install_property (object_class,
                                   PROP_DIM,
                                   g_param_spec_uint ("dimension",
                                                      NULL,
                                                      "Function dimension",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetFunc:eval-x:
   *
   * The evaluation point set by ncm_mset_func_set_eval_x(), or %NULL. Setting it
   * copies the vector, which must have #NcmMSetFunc:nvariables components.
   */
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

static const gdouble *_ncm_mset_func_get_x (NcmMSetFunc *func, const gdouble *x);

/**
 * ncm_mset_func_eval: (virtual eval)
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x: (array) (element-type double) (allow-none): function arguments
 * @res: (array) (element-type double): function values
 *
 * Evaluates @func at @x and stores its values in @res. If @x is %NULL, @func is
 * evaluated at the point set by ncm_mset_func_set_eval_x(); a function without
 * variables takes no arguments.
 *
 * @res must hold ncm_mset_func_get_dim() values. From language bindings @res is an
 * input array, so the values do not reach the caller; there, evaluate scalar
 * functions with ncm_mset_func_eval0(), ncm_mset_func_eval_nvar() or
 * ncm_mset_func_eval1().
 *
 */
void
ncm_mset_func_eval (NcmMSetFunc *func, NcmMSet *mset, gdouble *x, gdouble *res)
{
  NCM_MSET_FUNC_GET_CLASS (func)->eval (func, mset, _ncm_mset_func_get_x (func, x), res);
}

/**
 * ncm_mset_func_eval_array:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x: (array) (element-type double) (allow-none): function arguments
 *
 * Evaluates @func at @x and returns its values. If @x is %NULL, @func is evaluated
 * at the point set by ncm_mset_func_set_eval_x(); a function without variables
 * takes no arguments. This is ncm_mset_func_eval() for language bindings, where
 * the values of ncm_mset_func_eval() do not reach the caller.
 *
 * Returns: (array) (element-type double) (transfer full): the
 * ncm_mset_func_get_dim() values of @func.
 */
GArray *
ncm_mset_func_eval_array (NcmMSetFunc *func, NcmMSet *mset, GArray *x)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);
  GArray *res                     = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), self->dim);

  if ((x != NULL) && (x->len != self->nvar))
    g_error ("ncm_mset_func_eval_array: function `%s' takes %u variable(s), but %u argument(s) were given.",
             ncm_mset_func_peek_name (func), self->nvar, x->len);

  g_array_set_size (res, self->dim);
  NCM_MSET_FUNC_GET_CLASS (func)->eval (func, mset,
                                        _ncm_mset_func_get_x (func, (x != NULL) ? &g_array_index (x, gdouble, 0) : NULL),
                                        &g_array_index (res, gdouble, 0));

  return res;
}

/**
 * ncm_mset_func_eval_nvar:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x: (array) (element-type double) (allow-none): function arguments
 *
 * Evaluates the scalar function @func at @x and returns its value. If @x is %NULL,
 * @func is evaluated at the point set by ncm_mset_func_set_eval_x(); a function
 * without variables takes no arguments.
 *
 * Returns: function value.
 */
gdouble
ncm_mset_func_eval_nvar (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);
  gdouble res;

  if (self->dim != 1)
    g_error ("ncm_mset_func_eval_nvar: function `%s' has dimension %u, but only scalar functions return a single value.",
             ncm_mset_func_peek_name (func), self->dim);

  NCM_MSET_FUNC_GET_CLASS (func)->eval (func, mset, _ncm_mset_func_get_x (func, x), &res);

  return res;
}

/**
 * ncm_mset_func_eval0:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 *
 * Evaluates the scalar function @func and returns its value. The arguments are the
 * point set by ncm_mset_func_set_eval_x(), or none if @func has no variables.
 *
 * Returns: function value.
 */
gdouble
ncm_mset_func_eval0 (NcmMSetFunc *func, NcmMSet *mset)
{
  return ncm_mset_func_eval_nvar (func, mset, NULL);
}

/**
 * ncm_mset_func_eval1:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x: function argument
 *
 * Evaluates the scalar function @func at @x and returns its value. @func takes at
 * most one variable; a function without variables ignores @x. The point set by
 * ncm_mset_func_set_eval_x(), if any, is not used.
 *
 * Returns: function value.
 */
gdouble
ncm_mset_func_eval1 (NcmMSetFunc *func, NcmMSet *mset, const gdouble x)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);
  gdouble res;

  if ((self->dim != 1) || (self->nvar > 1))
    g_error ("ncm_mset_func_eval1: function `%s' takes %u variable(s) and has dimension %u, "
             "but only scalar functions of at most one variable are evaluated here.",
             ncm_mset_func_peek_name (func), self->nvar, self->dim);

  NCM_MSET_FUNC_GET_CLASS (func)->eval (func, mset, &x, &res);

  return res;
}

/*
 * _ncm_mset_func_get_x:
 * @func: a #NcmMSetFunc
 * @x: (nullable): function arguments
 *
 * Chooses the arguments of an evaluation. An explicit @x always wins; a %NULL
 * @x means the evaluation point set by ncm_mset_func_set_eval_x(). A function
 * without variables takes no arguments.
 *
 * Returns: the arguments to pass to the eval virtual function.
 */
static const gdouble *
_ncm_mset_func_get_x (NcmMSetFunc *func, const gdouble *x)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  if (x != NULL)
    return x;
  else if (self->eval_x != NULL)
    return ncm_vector_data (self->eval_x);
  else if (self->nvar == 0)
    return NULL;

  g_error ("ncm_mset_func: function `%s' takes %u variable(s), "
           "but it was called without arguments and no evaluation point is set.",
           ncm_mset_func_peek_name (func), self->nvar);

  return NULL;
}

/**
 * ncm_mset_func_eval_vector:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x_v: function arguments in a #NcmVector
 * @res_v: a #NcmVector to store the function values
 *
 * Evaluates the scalar function @func at each component of @x_v and stores the
 * values in the matching components of @res_v. @func takes at most one variable; a
 * function without variables ignores @x_v. The point set by
 * ncm_mset_func_set_eval_x(), if any, is not used.
 *
 */
void
ncm_mset_func_eval_vector (NcmMSetFunc *func, NcmMSet *mset, NcmVector *x_v, NcmVector *res_v)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);
  guint i;

  if ((self->dim != 1) || (self->nvar > 1))
    g_error ("ncm_mset_func_eval_vector: function `%s' takes %u variable(s) and has dimension %u, "
             "but only scalar functions of at most one variable are evaluated here.",
             ncm_mset_func_peek_name (func), self->nvar, self->dim);

  if (ncm_vector_len (res_v) != ncm_vector_len (x_v))
    g_error ("ncm_mset_func_eval_vector: %u argument(s) but room for %u value(s).",
             ncm_vector_len (x_v), ncm_vector_len (res_v));

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
 * Sets the evaluation point of @func to a copy of @x. Evaluations without
 * arguments use it, @func becomes constant (see ncm_mset_func_is_const()) and its
 * unique name and symbol encode it.
 *
 */
void
ncm_mset_func_set_eval_x (NcmMSetFunc *func, const gdouble *x, guint len)
{
  NcmMSetFuncPrivate * const self = ncm_mset_func_get_instance_private (func);

  if (len != self->nvar)
    g_error ("ncm_mset_func_set_eval_x: function `%s' takes %u variable(s), but the evaluation point has %u.",
             ncm_mset_func_peek_name (func), self->nvar, len);

  ncm_vector_clear (&self->eval_x);

  if (len > 0)
  {
    self->eval_x = ncm_vector_new (len);
    memcpy (ncm_vector_data (self->eval_x), x, len * sizeof (gdouble));
  }

  _ncm_mset_func_update_unames (func);
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
 * @x: (array) (element-type double) (allow-none): function arguments
 * @out: (inout) (allow-none) (transfer full): function gradient
 *
 * Computes the gradient of the scalar function @func at @x with respect to the free
 * parameters of @mset and stores it in @out. If *@out is %NULL, a new #NcmVector
 * is allocated; otherwise *@out must have one component per free parameter and is
 * overwritten.
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

  ncm_vector_free (x_v);
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

  _ncm_mset_func_update_unames (func);
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

  if ((self->uname == NULL) && (self->eval_x != NULL))
    _ncm_mset_func_update_unames (func);

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

  if ((self->usymbol == NULL) && (self->eval_x != NULL))
    _ncm_mset_func_update_unames (func);

  if (self->usymbol != NULL)
    return self->usymbol;
  else
    return ncm_mset_func_peek_symbol (func);
}

