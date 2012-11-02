/***************************************************************************
 *            nc_hicosmo.c
 *
 *  Tue May 29 19:23:52 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:nc_hicosmo
 * @title: Cosmological Model Pure Virtual Class
 * @short_description: FIXME
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hicosmo.h"
#include "math/ncm_cfg.h"

G_DEFINE_ABSTRACT_TYPE (NcHICosmo, nc_hicosmo, NCM_TYPE_MODEL);

static void
nc_hicosmo_init (NcHICosmo *object)
{
}

static void
nc_hicosmo_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_parent_class)->finalize (object);
}

gint32 NC_HICOSMO_ID = -1;

static void
nc_hicosmo_class_init (NcHICosmoClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  object_class->finalize = nc_hicosmo_finalize;
  ncm_model_class_register_id (NCM_MODEL_CLASS (klass));
  NC_HICOSMO_ID = NCM_MODEL_CLASS (klass)->model_id;
}

static void
_nc_hicosmo_log_all_models_go (GType model_type, guint n)
{
  guint nc, i, j;
  GType *models = g_type_children (model_type, &nc);
  for (i = 0; i < nc; i++)
  {
	guint ncc;
	GType *modelsc = g_type_children (models[i], &ncc);

	g_message ("#  ");
	for (j = 0; j < n; j++) g_message (" ");
	g_message ("%s\n", g_type_name (models[i]));
	if (ncc)
	  _nc_hicosmo_log_all_models_go (models[i], n + 2);

	g_free (modelsc);
  }
  g_free (models);
}

/**
 * nc_hicosmo_log_all_models:
 * @parent: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_log_all_models (GType parent)
{
  g_message ("# Registred NcHICosmos:%s are:\n", g_type_name (parent));
  _nc_hicosmo_log_all_models_go (parent, 0);
}

/**
 * nc_hicosmo_new_from_name:
 * @parent_type: FIXME
 * @model_name: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcHICosmo *
nc_hicosmo_new_from_name (GType parent_type, gchar *model_name)
{
  GObject *obj = ncm_cfg_create_from_string (model_name);
  GType model_type = G_OBJECT_TYPE (obj);

  if (!g_type_is_a (model_type, parent_type))
	g_error ("nc_hicosmo_new_from_name: NcHICosmo %s do not descend from %s\n", model_name, g_type_name (parent_type));
  return NC_HICOSMO (obj);
}

/**
 * nc_hicosmo_free:
 * @hic: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_free (NcHICosmo *hic)
{
  g_object_unref (hic);
}


/**
 * nc_hicosmo_set_H0_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc0,H0)
/**
 * nc_hicosmo_set_Omega_b_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc0,Omega_b)
/**
 * nc_hicosmo_set_Omega_r_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc0,Omega_r)
/**
 * nc_hicosmo_set_Omega_c_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc0,Omega_c)
/**
 * nc_hicosmo_set_Omega_t_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc0,Omega_t)
/**
 * nc_hicosmo_set_sigma_8_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc0,sigma_8)
/**
 * nc_hicosmo_set_T_gamma0_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc0,T_gamma0)

/**
 * nc_hicosmo_set_E2_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc1,E2)
/**
 * nc_hicosmo_set_dE2_dz_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc1,dE2_dz)
/**
 * nc_hicosmo_set_d2E2_dz2_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc1,d2E2_dz2)
/**
 * nc_hicosmo_set_cd_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc1,cd)
/**
 * nc_hicosmo_set_powspec_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc1,powspec)

/*
 * Inlined functions
 */

/**
 * nc_hicosmo_H0:
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_Omega_b:
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_Omega_r:
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_Omega_c:
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_Omega_t:
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_T_gamma0:
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_sigma_8:
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_z_lss:
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_E2:
 * @model: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_dE2_dz:
 * @model: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_d2E2_dz2:
 * @model: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_powspec:
 * @model: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_c_H0:
 * @model: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_Omega_k:
 * @model: a #NcHICosmo.
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_Omega_m:
 * @model: a #NcHICosmo.
 *
 * The matter density parameter is given by the baryonic plus
 * the cold dark matter density parameters.
 *
 * Returns: The matter density parameter at redshift zero.
 */
/**
 * nc_hicosmo_H:
 * @model: FIXME
 * @z: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_h:
 * @model: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_grad_h:
 * @model: a #NcHICosmo
 * @pt: a #NcmFitParams
 * @grad: a #NcmVector
 *
 * FIXME
 */
/**
 * nc_hicosmo_h2:
 * @model: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_grad_h2:
 * @model: a #NcHICosmo
 * @pt: a #NcmFitParams
 * @grad: a #NcmVector
 *
 * FIXME
 */
/**
 * nc_hicosmo_E:
 * @model: a #NcHICosmo
 * @z: redshift
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_grad_E:
 * @model: a #NcHICosmo
 * @pt: a #NcmFitParams
 * @z: redshift
 * @grad: a #NcmVector
 *
 * FIXME
 */
/**
 * nc_hicosmo_dH_dz:
 * @model: a #NcHICosmo
 * @z: redshift
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_q:
 * @model: a #NcHICosmo
 * @z: redshift
 *
 * FIXME
 *
 * Returns: FIXME
*/
/**
 * nc_hicosmo_qp:
 * @model: a #NcHICosmo
 * @z: redshift
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_j:
 * @model: a #NcHICosmo
 * @z: redshift
 *
 * FIXME
 *
 * Returns: FIXME
 */

static void
_nc_hicosmo_func0 (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcHICosmoFunc0 f0 = (NcHICosmoFunc0) obj;
  f[0] = f0 (NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID)));
}

static void
_nc_hicosmo_func1 (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcHICosmoFunc1 f1 = (NcHICosmoFunc1) obj;
  f[0] = f1 (NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID)), x[0]);
}

/**
 * ncm_mset_func_new_hicosmo_func0:
 * @f0: (scope notified): FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmMSetFunc *
ncm_mset_func_new_hicosmo_func0 (NcHICosmoFunc0 f0)
{
  return ncm_mset_func_new (&_nc_hicosmo_func0, 0, 1, f0, NULL);
}

/**
 * ncm_mset_func_new_hicosmo_func1:
 * @f1: (scope notified): FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmMSetFunc *
ncm_mset_func_new_hicosmo_func1 (NcHICosmoFunc1 f1)
{
  return ncm_mset_func_new (&_nc_hicosmo_func1, 1, 1, f1, NULL);
}
