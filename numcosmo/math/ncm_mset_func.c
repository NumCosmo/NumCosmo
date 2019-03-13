/***************************************************************************
 *            ncm_mset_func.c
 *
 *  Wed June 06 15:32:21 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:ncm_mset_func
 * @title: NcmMSetFunc
 * @short_description: Abstract class for arbitrary MSet functions.
 *
 * FIXME
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

G_DEFINE_ABSTRACT_TYPE (NcmMSetFunc, ncm_mset_func, G_TYPE_OBJECT);

static void
ncm_mset_func_init (NcmMSetFunc *func)
{
  func->nvar    = -1;
  func->dim     = -1;
  func->eval_x  = NULL;
  func->name    = NULL;
  func->symbol  = NULL;
  func->ns      = NULL;
  func->desc    = NULL;
  func->uname   = NULL;
  func->usymbol = NULL;
  func->diff    = ncm_diff_new ();
}

static void
_ncm_mset_func_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetFunc *func = NCM_MSET_FUNC (object);
  g_return_if_fail (NCM_IS_MSET_FUNC (object));

  switch (prop_id)
  {
    case PROP_NVAR:
		{
		  guint nvar = g_value_get_uint (value);
		  if (func->nvar == -1)
			{
				func->nvar = nvar;
			}
			else
			{
				g_assert_cmpuint (func->nvar, ==, nvar);
			}
      break;
		}
    case PROP_DIM:
		{
			guint dim = g_value_get_uint (value);
		  if (func->dim == -1)
			{			
				func->dim = dim;
			}
			else
			{
				g_assert_cmpuint (func->dim, ==, dim);
			}
      break;
		}
    case PROP_EVAL_X:
      ncm_vector_clear (&func->eval_x);
      func->eval_x = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_mset_func_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetFunc *func = NCM_MSET_FUNC (object);
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
      g_value_set_object (value, func->eval_x);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_mset_func_dispose (GObject *object)
{
  NcmMSetFunc *func = NCM_MSET_FUNC (object);

  ncm_vector_clear (&func->eval_x);
  ncm_diff_clear (&func->diff);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_func_parent_class)->dispose (object);
}

static void
_ncm_mset_func_finalize (GObject *object)
{
  NcmMSetFunc *func = NCM_MSET_FUNC (object);

  g_clear_pointer (&func->name,    g_free);
  g_clear_pointer (&func->symbol,  g_free);
  g_clear_pointer (&func->ns,      g_free);
  g_clear_pointer (&func->desc,    g_free);

  g_clear_pointer (&func->uname,   g_free);
  g_clear_pointer (&func->usymbol, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_func_parent_class)->finalize (object);
}

static void _ncm_mset_func_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res);

static void
ncm_mset_func_class_init (NcmMSetFuncClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

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
  
  klass->eval        = &_ncm_mset_func_eval;
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
 * FIXME
 *
 * Returns: (transfer full): FIXME
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
 * FIXME
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
 * FIXME
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
 * FIXME
 *
 * Returns: (element-type NcmMSetFunc) (transfer full): FIXME
 */
GPtrArray *
ncm_mset_func_array_new (void)
{
  return g_ptr_array_new_with_free_func ((GDestroyNotify) &ncm_mset_func_free);
}

/**
 * ncm_mset_func_eval: (virtual eval)
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x: (array) (element-type double): FIXME
 * @res: (array) (element-type double): FIXME
 * 
 * Evaluate the function.
 * 
 */
/**
 * ncm_mset_func_eval_nvar:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_mset_func_eval0:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_mset_func_eval1:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * ncm_mset_func_eval_vector:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x_v: FIXME
 * @res_v: FIXME
 *
 * FIXME
 *
 */

/**
 * ncm_mset_func_set_eval_x:
 * @func: a #NcmMSetFunc
 * @x: (in) (array length=len): FIXME
 * @len: FIXME
 *
 * FIXME
 *
 */
void 
ncm_mset_func_set_eval_x (NcmMSetFunc *func, gdouble *x, guint len)
{
  ncm_vector_clear (&func->eval_x);
  
  g_assert_cmpuint (func->nvar, ==, len);

  func->eval_x = ncm_vector_new_data_dup (x, func->nvar, 1);

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

    g_clear_pointer (&func->usymbol, g_free);
    g_clear_pointer (&func->uname,   g_free);

    func->usymbol = g_strjoin (NULL, func->symbol, args, NULL);

    {
      GRegex *reg  = g_regex_new ("[^0-9]+", 0, 0, NULL);
      gchar *nargs = g_regex_replace_literal (reg, args, -1, 0, "_", 0, NULL);

      func->uname  = g_strjoin (NULL, func->name, nargs, NULL);
      
      g_regex_unref (reg);      
    }
    
    g_free (args);
  }
}

/**
 * ncm_mset_func_is_scalar:
 * @func: a #NcmMSetFunc
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean 
ncm_mset_func_is_scalar (NcmMSetFunc *func)
{
  return (func->dim == 1);
}

/**
 * ncm_mset_func_is_vector:
 * @func: a #NcmMSetFunc
 * @dim: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean 
ncm_mset_func_is_vector (NcmMSetFunc *func, guint dim)
{
  return (func->dim == dim);
}

/**
 * ncm_mset_func_is_const:
 * @func: a #NcmMSetFunc
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean 
ncm_mset_func_is_const (NcmMSetFunc *func)
{
  return ((func->nvar == 0) || (func->eval_x != NULL));
}

/**
 * ncm_mset_func_has_nvar:
 * @func: a #NcmMSetFunc
 * @nvar: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean 
ncm_mset_func_has_nvar (NcmMSetFunc *func, guint nvar)
{
  return (func->nvar == nvar);
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
  _ncm_mset_func_numdiff_fparams_1 *nd = (_ncm_mset_func_numdiff_fparams_1 *)userdata;

  ncm_mset_fparams_set_vector (nd->mset, x_v);
  
  return ncm_mset_func_eval_nvar (nd->func, nd->mset, nd->x);
}

/**
 * ncm_mset_func_numdiff_fparams:
 * @func: a #NcmMSetFunc
 * @mset: a #NcmMSet
 * @x: FIXME
 * @out: (out) (transfer full): FIXME
 *
 * FIXME
 *
 */
void
ncm_mset_func_numdiff_fparams (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, NcmVector **out)
{
  const guint fparam_len = ncm_mset_fparam_len (mset);
  GArray *x_a            = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *x_v         = NULL;
  GArray *grad_a         = NULL;
  _ncm_mset_func_numdiff_fparams_1 nd;

  g_array_set_size (x_a, fparam_len);
  x_v = ncm_vector_new_array (x_a);
  ncm_mset_fparams_get_vector (mset, x_v);
  
  nd.mset  = mset;
  nd.func  = func;
  nd.x     = x;

  grad_a = ncm_diff_rf_d1_N_to_1 (func->diff, x_a, _mset_func_numdiff_fparams_1_val, &nd, NULL);

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
  return func->nvar;
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
  return func->dim;
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
  if (func->ns == NULL)
    func->ns = g_strdup (g_type_name (G_OBJECT_TYPE (func)));
  
  if (func->name == NULL)
    func->name = g_strdup_printf ("%s:no-name", func->ns);
  
  return func->name;
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
  if (func->ns == NULL)
    func->ns = g_strdup (g_type_name (G_OBJECT_TYPE (func)));
  
  if (func->symbol == NULL)
    func->symbol = g_strdup_printf ("%s:no-symbol", func->ns);

  return func->symbol;
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
  if (func->ns == NULL)
    func->ns = g_strdup (g_type_name (G_OBJECT_TYPE (func)));
  return func->ns;
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
  if (func->ns == NULL)
    func->ns = g_strdup (g_type_name (G_OBJECT_TYPE (func)));

  if (func->desc == NULL)
    func->desc = g_strdup_printf ("%s:no-desc", func->ns);

  return func->desc;
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
  if (func->uname != NULL)
    return func->uname;
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
  if (func->usymbol != NULL)
    return func->usymbol;
  else
    return ncm_mset_func_peek_symbol (func);
}
