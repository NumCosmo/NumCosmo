/***************************************************************************
 *            nc_hicosmo_de_reparam_ok.c
 *
 *  Mon October 26 10:50:22 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hicosmo_de_reparam_ok.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_hicosmo_de_reparam_ok
 * @title: NcHICosmoDEReparamOk
 * @short_description: Darkenergy model reparametrization $\Omega_{x0} \to \Omega_{k0}$.
 *
 * Object implementing a reparametrization for darkenergy models. It changes
 * $\Omega_{x0} \to \Omega_{k0}$.
 *
 */

#include "model/nc_hicosmo_de_reparam_ok.h"
#include "model/nc_hicosmo_de.h"

G_DEFINE_TYPE (NcHICosmoDEReparamOk, nc_hicosmo_de_reparam_ok, NCM_TYPE_REPARAM);

static void
nc_hicosmo_de_reparam_ok_init (NcHICosmoDEReparamOk *de_reparam_ok)
{
}

static void
nc_hicosmo_de_reparam_ok_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hicosmo_de_reparam_ok_parent_class)->constructed (object);
  {
    NcHICosmoDEReparamOk *reparam_Ok = NC_HICOSMO_DE_REPARAM_OK (object);
    ncm_reparam_set_param_desc_full (NCM_REPARAM (reparam_Ok), NC_HICOSMO_DE_OMEGA_X,
                                     "Omegak","\\Omega_{k0}", -5.0e-1, 5.0e-1, 1.0e-2,
                                     NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, 0.0, NCM_PARAM_TYPE_FIXED);
    NCM_REPARAM (reparam_Ok)->compat_type = NC_TYPE_HICOSMO_DE;
  }
}

static void
nc_hicosmo_de_reparam_ok_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_de_reparam_ok_parent_class)->finalize (object);
}

static gboolean _nc_hicosmo_de_reparam_ok_old2new (NcmReparam *reparam, NcmModel *model);
static gboolean _nc_hicosmo_de_reparam_ok_new2old (NcmReparam *reparam, NcmModel *model);
static gboolean _nc_hicosmo_de_reparam_ok_jac (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac);

static void
nc_hicosmo_de_reparam_ok_class_init (NcHICosmoDEReparamOkClass *klass)
{
  GObjectClass* object_class     = G_OBJECT_CLASS (klass);
  NcmReparamClass *reparam_class = NCM_REPARAM_CLASS (klass);

  object_class->constructed = nc_hicosmo_de_reparam_ok_constructed;
  object_class->finalize    = nc_hicosmo_de_reparam_ok_finalize;

  reparam_class->old2new = &_nc_hicosmo_de_reparam_ok_old2new;
  reparam_class->new2old = &_nc_hicosmo_de_reparam_ok_new2old;
  reparam_class->jac     = &_nc_hicosmo_de_reparam_ok_jac;
}

static gboolean
_nc_hicosmo_de_reparam_ok_old2new (NcmReparam *reparam, NcmModel *model)
{
  NcmVector *params = ncm_model_orig_params_peek_vector (model);
  const gdouble Omega_k0 = nc_hicosmo_Omega_k0 (NC_HICOSMO (model));
  ncm_vector_memcpy (reparam->new_params, params);
  ncm_vector_set (reparam->new_params, NC_HICOSMO_DE_OMEGA_X, Omega_k0);
  return TRUE;
}

static gboolean
_nc_hicosmo_de_reparam_ok_new2old (NcmReparam *reparam, NcmModel *model)
{
  NcmVector *params = ncm_model_orig_params_peek_vector (model);
  ncm_vector_memcpy (params, reparam->new_params);
  {
    NcHICosmo *cosmo = NC_HICOSMO (model);
    const gdouble Omega_x0 = 1.0 - (nc_hicosmo_Omega_m0 (cosmo) + nc_hicosmo_Omega_r0 (cosmo) + ncm_vector_get (reparam->new_params, NC_HICOSMO_DE_OMEGA_X));
    ncm_vector_set (params, NC_HICOSMO_DE_OMEGA_X, Omega_x0);
  }
  return TRUE;
}

static gboolean
_nc_hicosmo_de_reparam_ok_jac (NcmReparam *reparam, NcmModel *model, NcmMatrix *jac)
{
  g_assert_not_reached ();
}

/**
 * nc_hicosmo_de_reparam_ok_new: (constructor)
 * @length: number of parameters
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcHICosmoDEReparamOk
 */
NcHICosmoDEReparamOk *
nc_hicosmo_de_reparam_ok_new (guint length)
{
  NcHICosmoDEReparamOk *de_reparam_ok = g_object_new (NC_TYPE_HICOSMO_DE_REPARAM_OK,
                                                      "length", length,
                                                      NULL);
  return de_reparam_ok;
}
