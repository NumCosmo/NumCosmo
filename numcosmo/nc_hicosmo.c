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
 * @title: Cosmological Model Abstract Class
 * @short_description: Class for implementing homogeneous and isotropic cosmological models 
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hicosmo.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"

G_DEFINE_ABSTRACT_TYPE (NcHICosmo, nc_hicosmo, NCM_TYPE_MODEL);
G_DEFINE_BOXED_TYPE (NcHICosmoEOMAdiabZeta, nc_hicosmo_eom_adiab_zeta, nc_hicosmo_eom_adiab_zeta_dup, nc_hicosmo_eom_adiab_zeta_free);

static void
nc_hicosmo_init (NcHICosmo *object)
{
  NCM_UNUSED (object);
}

static void
nc_hicosmo_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_parent_class)->finalize (object);
}

static gboolean _nc_hicosmo_valid (NcmModel *model);
NCM_MSET_MODEL_REGISTER_ID (nc_hicosmo, NC_TYPE_HICOSMO);

static void
nc_hicosmo_class_init (NcHICosmoClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);
  
  object_class->finalize = nc_hicosmo_finalize;

  ncm_model_class_add_params (model_class, 0, 0, 1);
  
  ncm_mset_model_register_id (model_class, 
                              "NcHICosmo",
                              "Homogeneous and isotropic cosmological models.",
                              NULL);

  ncm_model_class_check_params_info (model_class);
  
  model_class->valid = &_nc_hicosmo_valid;
}

static gboolean 
_nc_hicosmo_valid (NcmModel *model)
{
  if (!NCM_MODEL_CLASS (nc_hicosmo_parent_class)->valid (model))
    return FALSE;
  /* Chain up : start */
  
  return (nc_hicosmo_E2 (NC_HICOSMO (model), 0) >= 0.0);
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
 * @cosmo_name: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcHICosmo *
nc_hicosmo_new_from_name (GType parent_type, gchar *cosmo_name)
{
  GObject *obj = ncm_serialize_global_from_string (cosmo_name);
  GType model_type = G_OBJECT_TYPE (obj);

  if (!g_type_is_a (model_type, parent_type))
	g_error ("nc_hicosmo_new_from_name: NcHICosmo %s do not descend from %s\n", cosmo_name, g_type_name (parent_type));
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
 * nc_hicosmo_set_z_lss_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc0,z_lss)

/**
 * nc_hicosmo_set_as_drag_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc0,as_drag)

/**
 * nc_hicosmo_set_xb_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc0,xb)

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
 * nc_hicosmo_set_cs2_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc1,cs2)

/**
 * nc_hicosmo_set_rhopp_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc1,rhopp)

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

/**
 * nc_hicosmo_set_wkb_adiab_theta_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc2,wkb_adiab_theta)

/**
 * nc_hicosmo_set_wkb_adiab_dmtheta_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc2,wkb_adiab_dmtheta)

/**
 * nc_hicosmo_set_eom_adiab_zeta_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFuncEOMAdiabZeta,eom_adiab_zeta)

/*
 * Inlined functions
 */

/**
 * nc_hicosmo_H0:
 * @cosmo: a #NcHICosmo.
 *
 * The value of the Hubble constant in unity of $ms^{-1}kpc^{-1}$.
 *
 * Returns: $H_0$
 */
/**
 * nc_hicosmo_Omega_b:
 * @cosmo: a #NcHICosmo.
 *
 * FIXME
 *
 * Returns: $\Omega_b$
 */
/**
 * nc_hicosmo_Omega_r:
 * @cosmo: a #NcHICosmo.
 *
 * FIXME
 *
 * Returns: $\Omega_r$
 */
/**
 * nc_hicosmo_Omega_c:
 * @cosmo: a #NcHICosmo.
 *
 * FIXME
 *
 * Returns: $\Omega_c$
 */
/**
 * nc_hicosmo_Omega_t:
 * @cosmo: a #NcHICosmo.
 *
 * FIXME
 *
 * Returns: $\Omega_t$
 */
/**
 * nc_hicosmo_T_gamma0:
 * @cosmo: a #NcHICosmo.
 *
 * FIXME
 *
 * Returns: $T_{\gamma0}$
 */
/**
 * nc_hicosmo_sigma_8:
 * @cosmo: a #NcHICosmo.
 *
 * FIXME
 *
 * Returns: $\sigma_8$
 */
/**
 * nc_hicosmo_z_lss:
 * @cosmo: a #NcHICosmo.
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_E2:
 * @cosmo: a #NcHICosmo.
 * @x: redshift
 *
 * Normalized Hubble function squared.
 *
 * Returns: $H^2/H_0^2$
 */
/**
 * nc_hicosmo_dE2_dz:
 * @cosmo: a #NcHICosmo.
 * @x: redshift
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_d2E2_dz2:
 * @cosmo: a #NcHICosmo.
 * @x: redshift
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_cs2:
 * @cosmo: a #NcHICosmo.
 * @x: redshift
 *
 * Speed of sound squared.
 *
 * Returns: $c_s^2$.
 */
/**
 * nc_hicosmo_rhopp:
 * @cosmo: a #NcHICosmo.
 * @x: redshift
 *
 * Energy density plus pressure.
 *
 * Returns: $\rho + p$.
 */
/**
 * nc_hicosmo_powspec:
 * @cosmo: a #NcHICosmo.
 * @x: redshift
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_adiabatic_zeta:
 * @cosmo: a #NcHICosmo.
 * @alpha: log-redshift time.
 * @k: mode.
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
/**
 * nc_hicosmo_c_H0:
 * @cosmo: a #NcHICosmo.
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_Omega_k:
 * @cosmo: a #NcHICosmo.
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_Omega_m:
 * @cosmo: a #NcHICosmo.
 *
 * The matter density parameter is given by the baryonic plus
 * the cold dark matter density parameters.
 *
 * Returns: The matter density parameter at redshift zero.
 */
/**
 * nc_hicosmo_H:
 * @cosmo: a #NcHICosmo.
 * @z: FIXME
 *
 * The value of the Hubble function in unity of $ms^{-1}kpc^{-1}$.
 *
 * Returns: $H(z)$
 */
/**
 * nc_hicosmo_h:
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_grad_h:
 * @cosmo: a #NcHICosmo
 * @pt: a #NcmFitParams
 * @grad: a #NcmVector
 *
 * FIXME
 */
/**
 * nc_hicosmo_h2:
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_Omega_bh2:
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_Omega_ch2:
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_Omega_rh2:
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_grad_h2:
 * @cosmo: a #NcHICosmo
 * @pt: a #NcmFitParams
 * @grad: a #NcmVector
 *
 * FIXME
 */
/**
 * nc_hicosmo_E:
 * @cosmo: a #NcHICosmo
 * @z: redshift
 *
 * This function computes the normalized Hubble function $E(z)$.
 *
 * Returns: $E(z)$.
 */
/**
 * nc_hicosmo_Em2:
 * @cosmo: a #NcHICosmo
 * @z: redshift
 *
 * This function computes the inverse of the square normalized Hubble function.
 *
 * Returns: $E(z)^{-2}$.
 */
/**
 * nc_hicosmo_grad_E:
 * @cosmo: a #NcHICosmo
 * @pt: a #NcmFitParams
 * @z: redshift
 * @grad: a #NcmVector
 *
 * FIXME
 */
/**
 * nc_hicosmo_dH_dz:
 * @cosmo: a #NcHICosmo
 * @z: redshift
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_q:
 * @cosmo: a #NcHICosmo
 * @z: redshift
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_dec:
 * @cosmo: a #NcHICosmo
 * @z: redshift
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_wec:
 * @cosmo: a #NcHICosmo
 * @z: redshift
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_qp:
 * @cosmo: a #NcHICosmo
 * @z: redshift
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_j:
 * @cosmo: a #NcHICosmo
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
  NCM_UNUSED (x);
  f[0] = f0 (NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ())));
}

static void
_nc_hicosmo_func1 (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcHICosmoFunc1 f1 = (NcHICosmoFunc1) obj;
  f[0] = f1 (NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ())), x[0]);
}

/**
 * nc_hicosmo_create_mset_func0:
 * @f0: (scope notified): FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmMSetFunc *
nc_hicosmo_create_mset_func0 (NcHICosmoFunc0 f0)
{
  return ncm_mset_func_new (&_nc_hicosmo_func0, 0, 1, f0, NULL);
}

/**
 * nc_hicosmo_create_mset_func1:
 * @f1: (scope notified): FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmMSetFunc *
nc_hicosmo_create_mset_func1 (NcHICosmoFunc1 f1)
{
  return ncm_mset_func_new (&_nc_hicosmo_func1, 1, 1, f1, NULL);
}


/**
 * nc_hicosmo_eom_adiab_zeta_dup:
 * @adiab_zeta: a #NcHICosmoEOMAdiabZeta.
 *
 * Duplicates @adiab_zeta.
 * 
 * Returns: (transfer full): a copy of @adiab_zeta.
 */
NcHICosmoEOMAdiabZeta *
nc_hicosmo_eom_adiab_zeta_dup (NcHICosmoEOMAdiabZeta *adiab_zeta)
{
  NcHICosmoEOMAdiabZeta *adiab_zeta_dup = g_new (NcHICosmoEOMAdiabZeta, 1);
  *adiab_zeta_dup = *adiab_zeta;
  return adiab_zeta_dup;
}

/**
 * nc_hicosmo_eom_adiab_zeta_free:
 * @adiab_zeta: a #NcHICosmoEOMAdiabZeta.
 *
 * Frees @adiab_zeta.
 * 
 */
void
nc_hicosmo_eom_adiab_zeta_free (NcHICosmoEOMAdiabZeta *adiab_zeta)
{
  g_free (adiab_zeta);
}
