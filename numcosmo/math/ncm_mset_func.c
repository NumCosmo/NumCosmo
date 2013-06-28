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
 * @title: A Function of NcmMSet
 * @short_description: Object representing a function on NcmMSet
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mset_func.h"
#include "math/util.h"

G_DEFINE_TYPE (NcmMSetFunc, ncm_mset_func, G_TYPE_OBJECT);

/**
 * ncm_mset_func_new:
 * @func: FIXME
 * @np: FIXME
 * @dim: FIXME
 * @obj: FIXME
 * @free: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmMSetFunc *
ncm_mset_func_new (NcmMSetFuncN func, guint np, guint dim, gpointer obj, GDestroyNotify free)
{
  NcmMSetFunc *gfunc = g_object_new (NCM_TYPE_MSET_FUNC, NULL);
  gfunc->func = func;
  gfunc->np  = np;
  gfunc->dim = dim;
  gfunc->obj = obj;
  gfunc->free = free;
  if (free)
	g_assert (obj != NULL);

  return gfunc;
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
 * ncm_mset_func_eval:
 * @func: FIXME
 * @mset: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_func_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x)
{
  gdouble res;
  func->func (mset, func->obj, x, &res);
  return res;
}

/**
 * ncm_mset_func_eval0:
 * @func: FIXME
 * @mset: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_func_eval0 (NcmMSetFunc *func, NcmMSet *mset)
{
  gdouble res;
#ifdef NCM_MSET_FUNC_CHECK_TYPE
  g_assert (func->dim == 1 && func->np == 0);
#endif
  func->func (mset, func->obj, NULL, &res);
  return res;
}

/**
 * ncm_mset_func_eval1:
 * @func: FIXME
 * @mset: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_mset_func_eval1 (NcmMSetFunc *func, NcmMSet *mset, const gdouble x)
{
  gdouble res;
#ifdef NCM_MSET_FUNC_CHECK_TYPE
  g_assert (func->dim == 1 && func->np == 1);
#endif
  func->func (mset, func->obj, &x, &res);
  return res;
}

/**
 * ncm_mset_func_is_scalar:
 * @func: FIXME
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
 * @func: FIXME
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
 * @func: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean 
ncm_mset_func_is_const (NcmMSetFunc *func)
{
  return (func->np == 0);
}

/**
 * ncm_mset_func_has_params:
 * @func: FIXME
 * @np: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean 
ncm_mset_func_has_params (NcmMSetFunc *func, guint np)
{
  return (func->np == np);
}

typedef struct __ncm_mset_func_numdiff_fparams_1
{
  NcmMSetFunc *func;
  NcmMSet *mset;
  const gdouble *x;
  guint n;
} __ncm_mset_func_numdiff_fparams_1;

static gdouble
_mset_func_numdiff_fparams_1_val (gdouble x, gpointer userdata)
{
  __ncm_mset_func_numdiff_fparams_1 *nd = (__ncm_mset_func_numdiff_fparams_1 *)userdata;
  ncm_mset_fparam_set (nd->mset, nd->n, x);
  return ncm_mset_func_eval (nd->func, nd->mset, nd->x);
}

/**
 * ncm_mset_func_numdiff_fparams:
 * @func: FIXME
 * @mset: FIXME
 * @x: FIXME
 * @out: (out) (transfer full): FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmVector *
ncm_mset_func_numdiff_fparams (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, NcmVector *out)
{
  gsl_function F;
  __ncm_mset_func_numdiff_fparams_1 nd;
  guint fparam_len = ncm_mset_fparam_len (mset);
  gint i;

  nd.mset  = mset;
  nd.func  = func;
  nd.x     = x;
  F.params = &nd;
  F.function = &_mset_func_numdiff_fparams_1_val;

  if (out == NULL)
    out = ncm_vector_new (fparam_len);
  else if (ncm_vector_len (out) != fparam_len)
  {
    ncm_vector_free (out);
    out = ncm_vector_new (fparam_len);
  }

  for (i = 0; i < fparam_len; i++)
  {
    const gdouble p = ncm_mset_fparam_get (mset, i);
    const gdouble p_scale = ncm_mset_fparam_get_scale (mset, i);
    gdouble err, diff;
    nd.n = i;
    diff = ncm_numdiff_1 (&F, p, p_scale, &err);
    ncm_vector_set (out, i, diff);
    ncm_mset_fparam_set (mset, i, p);
  }

  return out;
}

static void
ncm_mset_func_init (NcmMSetFunc *func)
{
  func->func = NULL;
  func->obj = NULL;
  func->free = NULL;
  func->np = 0;
  func->dim = 0;
}

static void
ncm_mset_func_dispose (GObject *object)
{
  NcmMSetFunc *func = NCM_MSET_FUNC (object);

  if (func->free && func->obj != NULL)
  {
    func->free (func->obj);
    func->obj = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_func_parent_class)->dispose (object);
}

static void
ncm_mset_func_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_func_parent_class)->finalize (object);
}

static void
ncm_mset_func_class_init (NcmMSetFuncClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->dispose  = ncm_mset_func_dispose;
  object_class->finalize = ncm_mset_func_finalize;
}
