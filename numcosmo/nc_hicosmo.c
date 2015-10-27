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
 * @title: NcHICosmo
 * @short_description: Abstract class for implementing homogeneous and isotropic cosmological models.
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

  klass->func_table   = g_array_new (FALSE, FALSE, sizeof (NcHICosmoFunc));
  klass->func_z_table = g_array_new (FALSE, FALSE, sizeof (NcHICosmoFuncZ));
  klass->func_hash    = g_hash_table_new (g_str_hash, g_str_equal);
  klass->func_z_hash  = g_hash_table_new (g_str_hash, g_str_equal);

  {
    NcHICosmoFunc func_table[] =
    {
      {"H0",         "Hubble constant.",                           &nc_hicosmo_H0,         NC_HICOSMO_IMPL_H0},
      {"Omega_b",    "Baryons density today.",                     &nc_hicosmo_Omega_b,    NC_HICOSMO_IMPL_Omega_b},
      {"Omega_g",    "Photons density today.",                     &nc_hicosmo_Omega_g,    NC_HICOSMO_IMPL_Omega_g},
      {"Omega_nu",   "Ultra-relativistic neutrino density today.", &nc_hicosmo_Omega_nu,   NC_HICOSMO_IMPL_Omega_nu},
      {"Omega_c",    "CDM density today.",                         &nc_hicosmo_Omega_c,    NC_HICOSMO_IMPL_Omega_c},
      {"Omega_r",    "Radiation density today.",                   &nc_hicosmo_Omega_r,    NC_HICOSMO_IMPL_Omega_r},
      {"Omega_t",    "Total energy density today.",                &nc_hicosmo_Omega_t,    NC_HICOSMO_IMPL_Omega_t},
      {"T_gamma0",   "Photons temperature today.",                 &nc_hicosmo_T_gamma0,   NC_HICOSMO_IMPL_T_gamma0},
      {"sigma_8",    "sigma_8.",                                   &nc_hicosmo_sigma_8,    NC_HICOSMO_IMPL_sigma_8},
      {"z_lss",      "redshift at lss.",                           &nc_hicosmo_z_lss,      NC_HICOSMO_IMPL_z_lss},
      {"as_drag",    "as_drag.",                                   &nc_hicosmo_as_drag,    NC_HICOSMO_IMPL_as_drag},
      {"xb",         "Bounce scale.",                              &nc_hicosmo_xb,         NC_HICOSMO_IMPL_xb},
      {"c_H0",       "Hubble radius.",                             &nc_hicosmo_c_H0,       NC_HICOSMO_IMPL_c_H0},
      {"Omega_k",    "Curvature scale.",                           &nc_hicosmo_Omega_k,    NC_HICOSMO_IMPL_Omega_k},
      {"Omega_m",    "Total cold matter density today.",           &nc_hicosmo_Omega_m,    NC_HICOSMO_IMPL_Omega_m},
      {"h",          "Adimensional Hubble constant.",              &nc_hicosmo_h,          NC_HICOSMO_IMPL_h},
      {"h2",         "Adimensional Hubble constant square.",       &nc_hicosmo_h2,         NC_HICOSMO_IMPL_h2},
      {"Omega_bh2",  "Baryons density today times h^2.",           &nc_hicosmo_Omega_bh2,  NC_HICOSMO_IMPL_Omega_bh2},
      {"Omega_gh2",  "Photons density today times h^2.",           &nc_hicosmo_Omega_gh2,  NC_HICOSMO_IMPL_Omega_gh2},
      {"Omega_nuh2", "UR Neutrinos density today times h^2.",      &nc_hicosmo_Omega_nuh2, NC_HICOSMO_IMPL_Omega_nuh2},
      {"Omega_ch2",  "CDM density today times h^2.",               &nc_hicosmo_Omega_ch2,  NC_HICOSMO_IMPL_Omega_ch2},
      {"Omega_rh2",  "Total radiation today times h^2.",           &nc_hicosmo_Omega_rh2,  NC_HICOSMO_IMPL_Omega_rh2},
      {"Omega_mh2",  "Total cold matter density today times h^2.", &nc_hicosmo_Omega_mh2,  NC_HICOSMO_IMPL_Omega_mh2},
    };
    const guint nfuncs = sizeof (func_table) / sizeof (NcHICosmoFunc);
    guint i;
    for (i = 0; i < nfuncs; i++)
    {
      g_array_append_val (klass->func_table, func_table[i]);
      g_hash_table_insert (klass->func_hash, (gchar *)func_table[i].name, GINT_TO_POINTER (i + 1));
    }
  }

  {
    NcHICosmoFuncZ func_z_table[] =
    {
      {"H",        "Hubble function.",                            &nc_hicosmo_H,               NC_HICOSMO_IMPL_H},
      {"dH_dz",    "Derivative of the Hubble function.",          &nc_hicosmo_dH_dz,           NC_HICOSMO_IMPL_dH_dz},
      {"E",        "Hubble function over H_0.",                   &nc_hicosmo_E,               NC_HICOSMO_IMPL_E},
      {"E2",       "Hubble function over H_0 squared.",           &nc_hicosmo_E2,              NC_HICOSMO_IMPL_E2},
      {"Em2",      "One over Hubble function over H_0 squared.",  &nc_hicosmo_Em2,             NC_HICOSMO_IMPL_Em2},
      {"dE2_dz",   "Derivative of the E2 function.",              &nc_hicosmo_dE2_dz,          NC_HICOSMO_IMPL_dE2_dz},
      {"d2E2_dz2", "Second derivative of the E2 function.",       &nc_hicosmo_d2E2_dz2,        NC_HICOSMO_IMPL_d2E2_dz2},
      {"q",        "Deceleration function.",                      &nc_hicosmo_q,               NC_HICOSMO_IMPL_q},
      {"qp",       "Derivative of the deceleration function.",    &nc_hicosmo_qp,              NC_HICOSMO_IMPL_j},
      {"j",        "Jerk function.",                              &nc_hicosmo_j,               NC_HICOSMO_IMPL_j},
      {"wec",      "WEC violation function.",                     &nc_hicosmo_wec,             NC_HICOSMO_IMPL_wec},
      {"dec",      "DEC violation function.",                     &nc_hicosmo_dec,             NC_HICOSMO_IMPL_dec},
    };
    const guint nfuncs = sizeof (func_z_table) / sizeof (NcHICosmoFuncZ);
    guint i;
    for (i = 0; i < nfuncs; i++)
    {
      g_array_append_val (klass->func_z_table, func_z_table[i]);
      g_hash_table_insert (klass->func_z_hash, (gchar *)func_z_table[i].name, GINT_TO_POINTER (i + 1));
    }
  }

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
 * nc_hicosmo_class_func_table:
 *
 * Returns: (transfer container) (element-type NcHICosmoFunc): the function table.
 */
GArray *
nc_hicosmo_class_func_table (void)
{
  NcHICosmoClass *hicosmo_class = g_type_class_ref (NC_TYPE_HICOSMO);
  GArray *func_table = g_array_ref (hicosmo_class->func_table);
  g_type_class_unref (hicosmo_class);
  return func_table;
}

/**
 * nc_hicosmo_class_func_z_table:
 *
 * Returns: (transfer container) (element-type NcHICosmoFuncZ): the function table.
 */
GArray *
nc_hicosmo_class_func_z_table (void)
{
  NcHICosmoClass *hicosmo_class = g_type_class_ref (NC_TYPE_HICOSMO);
  GArray *func_z_table = g_array_ref (hicosmo_class->func_z_table);
  g_type_class_unref (hicosmo_class);
  return func_z_table;
}

/**
 * nc_hicosmo_class_get_func:
 * @name: function name.
 *
 * Returns: (transfer none): the function @name or null if not found.
 */
NcHICosmoFunc *
nc_hicosmo_class_get_func (const gchar *name)
{
  NcHICosmoClass *hicosmo_class = g_type_class_ref (NC_TYPE_HICOSMO);
  gpointer f_i_ptr = g_hash_table_lookup (hicosmo_class->func_hash, name);
  guint f_i = GPOINTER_TO_INT (f_i_ptr) - 1;
  NcHICosmoFunc *f = (f_i_ptr != NULL) ? &g_array_index (hicosmo_class->func_table, NcHICosmoFunc, f_i) : NULL;
  g_type_class_unref (hicosmo_class);
  return f;
}

/**
 * nc_hicosmo_class_get_func_z:
 * @name: function name.
 *
 * Returns: (transfer none): the function @name or null if not found.
 */
NcHICosmoFuncZ *
nc_hicosmo_class_get_func_z (const gchar *name)
{
  NcHICosmoClass *hicosmo_class = g_type_class_ref (NC_TYPE_HICOSMO);
  gpointer f_i_ptr = g_hash_table_lookup (hicosmo_class->func_z_hash, name);
  guint f_i = GPOINTER_TO_INT (f_i_ptr) - 1;
  NcHICosmoFuncZ *fz = (f_i_ptr != NULL) ? &g_array_index (hicosmo_class->func_z_table, NcHICosmoFuncZ, f_i) : NULL;
  g_type_class_unref (hicosmo_class);
  return fz;
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
    g_error ("nc_hicosmo_new_from_name: NcHICosmo %s do not descend from %s.", cosmo_name, g_type_name (parent_type));
  return NC_HICOSMO (obj);
}

/**
 * nc_hicosmo_ref:
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcHICosmo *
nc_hicosmo_ref (NcHICosmo *cosmo)
{
  return g_object_ref (cosmo);
}

/**
 * nc_hicosmo_free:
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 */
void
nc_hicosmo_free (NcHICosmo *cosmo)
{
  g_object_unref (cosmo);
}

/**
 * nc_hicosmo_clear:
 * @cosmo: a #NcHICosmo
 *
 * The reference count of @cosmo is decreased and the pointer is set to NULL.
 *
 */
void
nc_hicosmo_clear (NcHICosmo **cosmo)
{
  g_clear_object (cosmo);
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
 * nc_hicosmo_set_Omega_g_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc0,Omega_g)
/**
 * nc_hicosmo_set_Omega_nu_impl: (skip)
 * @model_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcmModelFunc0,Omega_nu)
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
 * nc_hicosmo_Omega_g:
 * @cosmo: a #NcHICosmo.
 *
 * FIXME
 *
 * Returns: $\Omega_\gamma$
 */
/**
 * nc_hicosmo_Omega_nu:
 * @cosmo: a #NcHICosmo.
 *
 * FIXME
 *
 * Returns: $\Omega_\nu$
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
 * nc_hicosmo_Omega_gh2:
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_Omega_nuh2:
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

typedef struct __NcHICosmoArrayFunc
{
  guint size;
  NcHICosmoFunc1 f1;
} _NcHICosmoArrayFunc;

static void
_nc_hicosmo_arrayfunc1 (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  _NcHICosmoArrayFunc *af1 = (_NcHICosmoArrayFunc *) obj;
  guint i;

  for (i = 0; i < af1->size; i++)
  {
    f[i] = af1->f1 (NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ())), x[i]);
  }
}

/**
 * nc_hicosmo_create_mset_arrayfunc1:
 * @f1: (scope notified): FIXME
 * @size: function dimension
 *
 * Returns: (transfer full): FIXME
 */
NcmMSetFunc *
nc_hicosmo_create_mset_arrayfunc1 (NcHICosmoFunc1 f1, guint size)
{
  g_assert_cmpuint (size, !=, 0);
  _NcHICosmoArrayFunc *af1 = g_new (_NcHICosmoArrayFunc, 1);

  af1->size = size;
  af1->f1 = f1;

  return ncm_mset_func_new (&_nc_hicosmo_arrayfunc1, size, size, af1, g_free);
}
