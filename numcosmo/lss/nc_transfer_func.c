/***************************************************************************
 *            nc_transfer_func.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_transfer_func
 * @title: Transfer Function Abstract Class
 * @short_description: Defines the prototype of the #NcTransferFunc object.
 *
 * This module comprises the set of functions to compute the transfer function and
 * derived quantities.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_transfer_func.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"

G_DEFINE_ABSTRACT_TYPE (NcTransferFunc, nc_transfer_func, G_TYPE_OBJECT);

/**
 * nc_transfer_func_new_from_name:
 * @transfer_name: string which specifies the transfer function type.
 *
 * This function returns a new #NcTransferFunc whose type is defined by @transfer_name.
 *
 * Returns: A new #NcTransferFunc.
 */
NcTransferFunc *
nc_transfer_func_new_from_name (gchar *transfer_name)
{
  GObject *obj = ncm_serialize_global_from_string (transfer_name);
  GType transfer_type = G_OBJECT_TYPE (obj);

  if (!g_type_is_a (transfer_type, NC_TYPE_TRANSFER_FUNC))
	g_error ("nc_transfer_func_new_from_name: NcTransferFunc %s do not descend from %s\n", transfer_name, g_type_name (NC_TYPE_TRANSFER_FUNC));
  return NC_TRANSFER_FUNC (obj);
}

/**
 * nc_transfer_func_prepare:
 * @tf: a #NcTransferFunc.
 * @model: a #NcHICosmo.
 *
 * FIXME
 *
*/
void
nc_transfer_func_prepare (NcTransferFunc *tf, NcHICosmo *model)
{
  if (tf->ctrl == NULL)
    g_error ("You should allocate a new tf using transfer_func_new_from_name(...)");

  if (ncm_model_ctrl_update (tf->ctrl, NCM_MODEL(model)))
    NC_TRANSFER_FUNC_GET_CLASS (tf)->prepare (tf, model);
}

/**
 * nc_transfer_func_eval:
 * @tf: a #NcTransferFunc.
 * @model: a #NcHICosmo.
 * @kh: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_transfer_func_eval (NcTransferFunc *tf, NcHICosmo *model, gdouble kh)
{
  if (tf->ctrl == NULL)
    g_error ("You should allocate a new tf using transfer_func_new_from_name(...)");

  if (ncm_model_ctrl_update (tf->ctrl, NCM_MODEL(model)))
    NC_TRANSFER_FUNC_GET_CLASS (tf)->prepare (tf, model);

  return NC_TRANSFER_FUNC_GET_CLASS (tf)->calc (tf, kh);
}

/**
 * nc_transfer_func_matter_powerspectrum:
 * @tf: a #NcTransferFunc.
 * @model: a #NcHICosmo.
 * @kh: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_transfer_func_matter_powerspectrum (NcTransferFunc *tf, NcHICosmo *model, gdouble kh)
{
  if (tf->ctrl == NULL)
    g_error ("You should allocate a new tf using transfer_func_new_from_name(...)");

  if (ncm_model_ctrl_update (tf->ctrl, NCM_MODEL(model)))
    NC_TRANSFER_FUNC_GET_CLASS (tf)->prepare (tf, model);

  return NC_TRANSFER_FUNC_GET_CLASS (tf)->calc_matter_P (tf, model, kh);
}

/**
 * nc_transfer_func_free:
 * @tf: a #NcTransferFunc.
 *
 * Atomically decrements the reference count of @tf by one. If the reference count drops to 0,
 * all memory allocated by @tf is released.
 *
*/
void
nc_transfer_func_free (NcTransferFunc *tf)
{
  g_object_unref (tf);
}

/**
 * nc_transfer_func_clear:
 * @tf: a #NcTransferFunc.
 *
 * Atomically decrements the reference count of @tf by one. If the reference count drops to 0,
 * all memory allocated by @tf is released. Set the pointer to NULL.
 *
*/
void
nc_transfer_func_clear (NcTransferFunc **tf)
{
  g_clear_object (tf);
}

static void
nc_transfer_func_init (NcTransferFunc *tf)
{
  /* TODO: Add initialization code here */
  tf->ctrl = ncm_model_ctrl_new (NULL);
}

static void
_nc_transfer_func_dispose (GObject * object)
{
  NcTransferFunc *tf = NC_TRANSFER_FUNC (object);

  ncm_model_ctrl_clear (&tf->ctrl);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_transfer_func_parent_class)->dispose (object);
}

static void
_nc_transfer_func_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_transfer_func_parent_class)->finalize (object);
}

static void
nc_transfer_func_class_init (NcTransferFuncClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->dispose = _nc_transfer_func_dispose;
  object_class->finalize = _nc_transfer_func_finalize;
}

