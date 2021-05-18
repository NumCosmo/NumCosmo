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
 * @title: NcTransferFunc
 * @short_description: Abstrac class for perturbation transfer function.
 * @stability: Stable
 * @include: numcosmo/lss/nc_transfer_func.h
 *
 * This module comprises the set of functions to compute the transfer function and
 * derived quantities. The applied $k$ unit is $[\mathrm{Mpc}^{-1}]$. 
 *
 * The transfer function, $T(k)$, is defined as,
 * \begin{equation*}
 *    T(k) \equiv \frac{\hat{\delta}(k, z=0)}{\hat{\delta}(k, z=\infty)} \frac{\hat{\delta}(k=0, z=\infty)}{\hat{\delta}(k=0, z=0)} \, ,
 * \end{equation*}
 * where $\hat{\delta}(k, z)$ is the density perturbation, in Fourier space,
 * for mode (wavenumber) $k$ at redshift $z$. By definition, we have
 * $$ \lim_{k \rightarrow 0} T(k) \rightarrow  1 \, .$$
 *
 * See [Eisenstein and Hu (1998)][XEisenstein1998] [[arXiv](https://arxiv.org/abs/astro-ph/9709112)] for more details.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_transfer_func.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"

G_DEFINE_ABSTRACT_TYPE (NcTransferFunc, nc_transfer_func, G_TYPE_OBJECT);

static void
nc_transfer_func_init (NcTransferFunc *tf)
{
  tf->ctrl_cosmo = ncm_model_ctrl_new (NULL);
  tf->ctrl_reion = ncm_model_ctrl_new (NULL);
}

static void
_nc_transfer_func_dispose (GObject *object)
{
  NcTransferFunc *tf = NC_TRANSFER_FUNC (object);
  
  ncm_model_ctrl_clear (&tf->ctrl_cosmo);
  ncm_model_ctrl_clear (&tf->ctrl_reion);
  
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
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  
  object_class->dispose  = _nc_transfer_func_dispose;
  object_class->finalize = _nc_transfer_func_finalize;
}

/**
 * nc_transfer_func_new_from_name:
 * @transfer_name: string which specifies the transfer function type
 *
 * This function returns a new #NcTransferFunc whose type is defined by @transfer_name.
 *
 * Returns: (transfer full): A new #NcTransferFunc.
 */
NcTransferFunc *
nc_transfer_func_new_from_name (gchar *transfer_name)
{
  GObject *obj        = ncm_serialize_global_from_string (transfer_name);
  GType transfer_type = G_OBJECT_TYPE (obj);
  
  if (!g_type_is_a (transfer_type, NC_TYPE_TRANSFER_FUNC))
    g_error ("nc_transfer_func_new_from_name: NcTransferFunc %s do not descend from %s.",
             transfer_name,
             g_type_name (NC_TYPE_TRANSFER_FUNC));
  
  return NC_TRANSFER_FUNC (obj);
}

/**
 * nc_transfer_func_ref:
 * @tf: a #NcTransferFunc
 *
 * Increases the reference count of @tf atomically.
 *
 * Returns: (transfer full): @tf.
 */
NcTransferFunc *
nc_transfer_func_ref (NcTransferFunc *tf)
{
  return g_object_ref (tf);
}

/**
 * nc_transfer_func_free:
 * @tf: a #NcTransferFunc
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
 * @tf: a #NcTransferFunc
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

/**
 * nc_transfer_func_prepare:
 * @tf: a #NcTransferFunc
 * @cosmo: a #NcHICosmo
 *
 * Prepares the transfer function @tf with model @cosmo,
 * such that one can evaluate it (#nc_transfer_func_eval).
 *
 */
void
nc_transfer_func_prepare (NcTransferFunc *tf, NcHICosmo *cosmo)
{
  NC_TRANSFER_FUNC_GET_CLASS (tf)->prepare (tf, cosmo);
  
  ncm_model_ctrl_update (tf->ctrl_cosmo, NCM_MODEL (cosmo));
}

/**
 * nc_transfer_func_prepare_if_needed:
 * @tf: a #NcTransferFunc
 * @cosmo: a #NcHICosmo
 *
 * Prepares (if necessary) the transfer function @tf with model @cosmo.
 *
 */
void
nc_transfer_func_prepare_if_needed (NcTransferFunc *tf, NcHICosmo *cosmo)
{
  gboolean cosmo_up = ncm_model_ctrl_update (tf->ctrl_cosmo, NCM_MODEL (cosmo));
  
  if (cosmo_up)
    NC_TRANSFER_FUNC_GET_CLASS (tf)->prepare (tf, cosmo);
}

/**
 * nc_transfer_func_eval:
 * @tf: a #NcTransferFunc $T(k)$
 * @cosmo: a #NcHICosmo
 * @kh: mode (wavenumber)
 *
 * The transfer function @tf value at mode (wavenumber) 
 * @kh (in $Mpc^{-1}$ units) with model @cosmo.
 *
 * Returns: $T(k)$.
 */
gdouble
nc_transfer_func_eval (NcTransferFunc *tf, NcHICosmo *cosmo, gdouble kh)
{
  NCM_CHECK_PREPARED (tf, nc_transfer_func_eval);
  
  return NC_TRANSFER_FUNC_GET_CLASS (tf)->calc (tf, kh);
}

