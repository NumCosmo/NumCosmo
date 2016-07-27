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
#include "nc_hiprim.h"
#include "nc_hireion.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"

G_DEFINE_ABSTRACT_TYPE (NcHICosmo, nc_hicosmo, NCM_TYPE_MODEL);

static void
nc_hicosmo_init (NcHICosmo *cosmo)
{
  cosmo->prim  = NULL;
  cosmo->reion = NULL;
}

static void
nc_hicosmo_dispose (GObject *object)
{
  NcHICosmo *cosmo = NC_HICOSMO (object);

  nc_hiprim_clear (&cosmo->prim);
  nc_hireion_clear (&cosmo->reion);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_parent_class)->dispose (object);
}

static void
nc_hicosmo_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_parent_class)->finalize (object);
}

static gboolean _nc_hicosmo_valid (NcmModel *model);
NCM_MSET_MODEL_REGISTER_ID (nc_hicosmo, NC_TYPE_HICOSMO);
static void _nc_hicosmo_add_submodel (NcmModel *model, NcmModel *submodel);

static void
nc_hicosmo_class_init (NcHICosmoClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->dispose  = &nc_hicosmo_dispose;
  object_class->finalize = &nc_hicosmo_finalize;

  ncm_model_class_set_name_nick (model_class, "Abstract class for HI cosmological models.", "NcHICosmo");
  ncm_model_class_add_params (model_class, 0, 0, 1);

  ncm_mset_model_register_id (model_class,
                              "NcHICosmo",
                              "Homogeneous and isotropic cosmological models.",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_check_params_info (model_class);

  model_class->valid        = &_nc_hicosmo_valid;
  model_class->add_submodel = &_nc_hicosmo_add_submodel;
  
  klass->func_table   = g_array_new (FALSE, FALSE, sizeof (NcHICosmoFunc));
  klass->func_z_table = g_array_new (FALSE, FALSE, sizeof (NcHICosmoFuncZ));
  klass->func_hash    = g_hash_table_new (g_str_hash, g_str_equal);
  klass->func_z_hash  = g_hash_table_new (g_str_hash, g_str_equal);

  
  {
    NcHICosmoFunc func_table[] =
    {
      {"H0",          "Hubble constant.",                           &nc_hicosmo_H0,          NC_HICOSMO_IMPL_H0},
      {"RH_Mpc",      "Hubble radius today in Mpc.",                &nc_hicosmo_RH_Mpc,      NC_HICOSMO_IMPL_RH_Mpc},
      {"h",           "Dimensionless Hubble constant.",             &nc_hicosmo_h,           NC_HICOSMO_IMPL_h},
      {"h2",          "Dimensionless Hubble constant square.",      &nc_hicosmo_h2,          NC_HICOSMO_IMPL_h2},
      {"Omega_b0",    "Baryons density today.",                     &nc_hicosmo_Omega_b0,    NC_HICOSMO_IMPL_Omega_b0},
      {"Omega_c0",    "CDM density today.",                         &nc_hicosmo_Omega_c0,    NC_HICOSMO_IMPL_Omega_c0},
      {"Omega_g0",    "Photons density today.",                     &nc_hicosmo_Omega_g0,    NC_HICOSMO_IMPL_Omega_g0},
      {"Omega_nu0",   "Ultra-relativistic neutrino density today.", &nc_hicosmo_Omega_nu0,   NC_HICOSMO_IMPL_Omega_nu0},
      {"Omega_m0",    "Total dust matter density today.",           &nc_hicosmo_Omega_m0,    NC_HICOSMO_IMPL_Omega_m0},
      {"Omega_r0",    "Radiation density today.",                   &nc_hicosmo_Omega_r0,    NC_HICOSMO_IMPL_Omega_r0},
      {"Omega_t0",    "Total energy density today.",                &nc_hicosmo_Omega_t0,    NC_HICOSMO_IMPL_Omega_t0},
      {"Omega_k0",    "Curvature scale.",                           &nc_hicosmo_Omega_k0,    NC_HICOSMO_IMPL_Omega_k0},
      {"Omega_b0h2",  "Baryons density today times h^2.",           &nc_hicosmo_Omega_b0h2,  NC_HICOSMO_IMPL_Omega_b0h2},
      {"Omega_c0h2",  "CDM density today times h^2.",               &nc_hicosmo_Omega_c0h2,  NC_HICOSMO_IMPL_Omega_c0h2},
      {"Omega_g0h2",  "Photons density today times h^2.",           &nc_hicosmo_Omega_g0h2,  NC_HICOSMO_IMPL_Omega_g0h2},
      {"Omega_nu0h2", "UR Neutrinos density today times h^2.",      &nc_hicosmo_Omega_nu0h2, NC_HICOSMO_IMPL_Omega_nu0h2},
      {"Omega_m0h2",  "Total dust matter density today times h^2.", &nc_hicosmo_Omega_m0h2,  NC_HICOSMO_IMPL_Omega_m0h2},
      {"Omega_r0h2",  "Total radiation today times h^2.",           &nc_hicosmo_Omega_r0h2,  NC_HICOSMO_IMPL_Omega_r0h2},
      {"T_gamma0",    "Photons temperature today.",                 &nc_hicosmo_T_gamma0,    NC_HICOSMO_IMPL_T_gamma0},
      {"Yp_4He",      "Primordial Helium mass fraction.",           &nc_hicosmo_Yp_4He,      NC_HICOSMO_IMPL_Yp_4He},
      {"H_Yp",        "Primordial Hydrogen mass fraction.",         &nc_hicosmo_Yp_1H,       NC_HICOSMO_IMPL_H_Yp},
      {"XHe",         "Primordial Helium abundance.",               &nc_hicosmo_XHe,         NC_HICOSMO_IMPL_XHe},
      {"z_lss",       "redshift at lss.",                           &nc_hicosmo_z_lss,       NC_HICOSMO_IMPL_z_lss},
      {"as_drag",     "as_drag.",                                   &nc_hicosmo_as_drag,     NC_HICOSMO_IMPL_as_drag},
      {"xb",          "Bounce scale.",                              &nc_hicosmo_xb,          NC_HICOSMO_IMPL_xb},
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
      {"dec",      "DEC violation function.",                     &nc_hicosmo_dec,             NC_HICOSMO_IMPL_dec},
      {"wec",      "WEC violation function.",                     &nc_hicosmo_wec,             NC_HICOSMO_IMPL_wec},
      {"qp",       "Derivative of the deceleration function.",    &nc_hicosmo_qp,              NC_HICOSMO_IMPL_j},
      {"j",        "Jerk function.",                              &nc_hicosmo_j,               NC_HICOSMO_IMPL_j},
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

  return (nc_hicosmo_E2 (NC_HICOSMO (model), 0.0) >= 0.0);
}

static void 
_nc_hicosmo_add_submodel (NcmModel *model, NcmModel *submodel)
{
  /* Chain up : start */
  NCM_MODEL_CLASS (nc_hicosmo_parent_class)->add_submodel (model, submodel);
  {
    NcHICosmo *cosmo = NC_HICOSMO (model);

    if (ncm_model_id (submodel) == nc_hiprim_id ())
    {
      nc_hiprim_clear (&cosmo->prim);
      cosmo->prim = nc_hiprim_ref (NC_HIPRIM (submodel));
    }
    else if (ncm_model_id (submodel) == nc_hireion_id ())
    {
      nc_hireion_clear (&cosmo->reion);
      cosmo->reion = nc_hireion_ref (NC_HIREION (submodel));
    }
  }
}

/**
 * nc_hicosmo_set_H0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: an implementation of H0.
 *
 * Sets the implementation of H0 to @f.
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,H0)
/**
 * nc_hicosmo_set_Omega_b0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,Omega_b0)
/**
 * nc_hicosmo_set_Omega_g0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,Omega_g0)
/**
 * nc_hicosmo_set_Omega_nu0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,Omega_nu0)
/**
 * nc_hicosmo_set_Omega_r0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,Omega_r0)
/**
 * nc_hicosmo_set_Omega_c0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,Omega_c0)
/**
 * nc_hicosmo_set_Omega_t0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,Omega_t0)
/**
 * nc_hicosmo_set_T_gamma0_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,T_gamma0)
/**
 * nc_hicosmo_set_Yp_4He_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,Yp_4He)
/**
 * nc_hicosmo_set_z_lss_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,z_lss)

/**
 * nc_hicosmo_set_as_drag_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,as_drag)

/**
 * nc_hicosmo_set_xb_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc0,xb)

/**
 * nc_hicosmo_set_E2_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc1Z,E2)

/**
 * nc_hicosmo_set_dE2_dz_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc1Z,dE2_dz)

/**
 * nc_hicosmo_set_d2E2_dz2_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc1Z,d2E2_dz2)

/**
 * nc_hicosmo_set_bgp_cs2_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc1Z,bgp_cs2)

/**
 * nc_hicosmo_set_Dc_impl: (skip)
 * @model_class: a #NcmModelClass
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO,NcHICosmo,nc_hicosmo,NcHICosmoFunc1Z,Dc)

/**
 * nc_hicosmo_new_from_name:
 * @parent_type: parent's #GType
 * @cosmo_name: Cosmological model's name
 *
 * Creates a new instance of @cosmo_name,
 * asserting that it descends from @parent_type.
 *
 * Returns: newly created @cosmo_name object.
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
 * Increases the reference count of @cosmo by one.
 *
 * Returns: (transfer full): @cosmo.
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
 * Decreases the reference count of @cosmo by one.
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
 * @parent: #GType of the parent model
 *
 * Logs all models descending from @parent.
 *
 */
void
nc_hicosmo_log_all_models (GType parent)
{
  g_message ("# Registred NcHICosmo:%s are:\n", g_type_name (parent));
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

/*
 * Inlined functions
 */
/**
 * nc_hicosmo_H0: (virtual H0)
 * @cosmo: a #NcHICosmo
 *
 * The value of the Hubble constant in unit of $\mathrm{m}\,\mathrm{s}^{-1}\,\mathrm{kpc}^{-1}$,
 * see ncm_c_kpc().
 *
 * Returns: $H_0 \left[\mathrm{m}\,\mathrm{s}^{-1}\,\mathrm{kpc}^{-1}\right]$
 */
/**
 * nc_hicosmo_RH_Mpc:
 * @cosmo: a #NcHICosmo
 *
 * Calculates the Hubble radius in unit of
 * Mpc, i.e., $R_H = (c / (H_0 \times 1\,\mathrm{Mpc}))$. 
 *
 * Returns: $R_H \left[\mathrm{Mpc}\right]$.
 */
/**
 * nc_hicosmo_h:
 * @cosmo: a #NcHICosmo
 *
 * Reduced Hubble constant, $h \equiv H_0 / (1\times\mathrm{m}\mathrm{s}^{-1}\mathrm{kpc}^{-1})$.
 *
 * Returns: $h$.
 */
/**
 * nc_hicosmo_h2:
 * @cosmo: a #NcHICosmo
 *
 * Reduced Hubble constant [nc_hicosmo_h()] squared $h^2$.
 *
 * Returns: $h^2$.
 */
/**
 * nc_hicosmo_Omega_b0: (virtual Omega_b0)
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless baryon density today $\Omega_{b0} = \rho_{b0} / \rho_{\mathrm{crit}0}$,
 * see nc_hicosmo_crit_density().
 *
 * Returns: $\Omega_{b0}$
 */
/**
 * nc_hicosmo_Omega_c0: (virtual Omega_c0)
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless cold dark matter density today $\Omega_{c0} = \rho_{c0} / \rho_{\mathrm{crit}0}$,
 * see nc_hicosmo_crit_density().
 *
 * Returns: $\Omega_{c0}$
 */
/**
 * nc_hicosmo_Omega_g0: (virtual Omega_g0)
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless photon density today $\Omega_{\gamma0} = \rho_{\gamma0} / \rho_{\mathrm{crit}0}$,
 * see nc_hicosmo_crit_density().
 *
 * Returns: $\Omega_{\gamma0}$
 */
/**
 * nc_hicosmo_Omega_nu0: (virtual Omega_nu0)
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless relativistic neutrinos density today $\Omega_{\nu0} = \rho_{\nu0} / \rho_{\mathrm{crit}0}$,
 * see nc_hicosmo_crit_density().
 *
 * Returns: $\Omega_{\nu0}$
 */
/**
 * nc_hicosmo_Omega_m0:
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless total dust density today $\Omega_{m0} = \rho_{m0} / \rho_{\mathrm{crit}0}$,
 * see nc_hicosmo_crit_density().
 *
 * Returns: $\Omega_{m0}$.
 */
/**
 * nc_hicosmo_Omega_r0: (virtual Omega_r0)
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless total radiation density today $\Omega_{r0} = \rho_{r0} / \rho_{\mathrm{crit}0}$,
 * see nc_hicosmo_crit_density().
 *
 * Returns: $\Omega_{r0}$
 */
/**
 * nc_hicosmo_Omega_b0h2:
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless baryon density today [nc_hicosmo_Omega_b0()] times $h^2$.
 *
 * Returns: $\Omega_{b0}h^2$.
 */
/**
 * nc_hicosmo_Omega_c0h2:
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless cold dark matter density today [nc_hicosmo_Omega_c0()] times $h^2$.
 *
 * Returns: $\Omega_{c0}h^2$.
 */
/**
 * nc_hicosmo_Omega_g0h2:
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless photon density today [nc_hicosmo_Omega_g0()] times $h^2$.
 *
 * Returns: $\Omega_{\gamma0}h^2$.
 */
/**
 * nc_hicosmo_Omega_nu0h2:
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless relativistic neutrinos density today [nc_hicosmo_Omega_nu0()] times $h^2$.
 *
 * Returns: $\Omega_{\nu0}h^2$.
 */
/**
 * nc_hicosmo_Omega_m0h2:
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless total dust density today [nc_hicosmo_Omega_m0()] times $h^2$.
 *
 * Returns: $\Omega_{m0}h^2$.
 */
/**
 * nc_hicosmo_Omega_r0h2:
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless total radiation density today [nc_hicosmo_Omega_r0()] times $h^2$.
 *
 * Returns: $\Omega_{r0}h^2$.
 */
/**
 * nc_hicosmo_Omega_t0: (virtual Omega_t0)
 * @cosmo: a #NcHICosmo
 *
 * Dimensionless total matter density today $\Omega_{t0} = \rho_{t0} / \rho_{\mathrm{crit}0}$,
 * see nc_hicosmo_crit_density().
 *
 * Returns: $\Omega_{t0}$
 */
/**
 * nc_hicosmo_Omega_k0:
 * @cosmo: a #NcHICosmo
 *
 * The curvature parameter today, $\Omega_{k0}$.
 *
 * Returns: $\Omega_{k0}$.
 */

/**
 * nc_hicosmo_T_gamma0: (virtual T_gamma0)
 * @cosmo: a #NcHICosmo
 *
 * Gets the cosmic microwave background radiation temperature today.
 *
 * Returns: $T_{\gamma0} \left[\mathrm{K}\right]$.
 */
/**
 * nc_hicosmo_Yp_4He: (virtual Yp_4He)
 * @cosmo: a #NcHICosmo
 *
 * Gets the primordial Helium mass fraction, i.e., 
 * $$Y_p = \frac{m_\mathrm{He}n_\mathrm{He}}
 * {m_\mathrm{He}n_\mathrm{He} + m_\mathrm{H}n_\mathrm{H}},$$ where $m_\mathrm{He}$, 
 * $n_\mathrm{He}$, $m_\mathrm{H}$ and $m_\mathrm{H}$ are respectively Helium-4 mass and 
 * number density and Hydrogen-1 mass and number density.
 *
 * Returns: $Y_p$.
 */
/**
 * nc_hicosmo_Yp_1H:
 * @cosmo: a #NcHICosmo
 *
 * The primordial hydrogen mass fraction $$Y_{\text{1H}p} = 1 - Y_p,$$
 * where $Y_p$ is the helium mass fraction, see nc_hicosmo_Yp_4He().
 *
 * Returns: $Y_{\text{1H}p}$.
 */
/**
 * nc_hicosmo_XHe:
 * @cosmo: a #NcHICosmo
 *
 * The primordial Helium to Hydrogen ratio $$X_\text{He} = 
 * \frac{n_\text{He}}{n_\text{H}} = \frac{m_\text{1H}}{m_\text{4He}}
 * \frac{Y_p}{Y_{\text{1H}p}},$$ see nc_hicosmo_Yp_1H() and nc_hicosmo_Yp_4He().
 * 
 * Returns: The primordial Helium to Hydrogen ratio $X_\text{He}$.
 */
/**
 * nc_hicosmo_crit_density:
 * @cosmo: a #NcHICosmo
 *
 * Calculares the critical density $\rho_\mathrm{crit}$ using 
 * ncm_c_crit_density_h2() $\times$ nc_hicosmo_h2().
 * 
 * Returns: The critical density $\rho_{\mathrm{crit}0}$.
 */
/**
 * nc_hicosmo_baryon_density:
 * @cosmo: a #NcHICosmo
 * 
 * Calculares the baryon density $\rho_{b0} = \rho_{\mathrm{crit}0} \Omega_{b0}$ 
 * using nc_hicosmo_crit_density() $\times$ nc_hicosmo_Omega_b0().
 * 
 * Returns: The baryon density $\rho_{b0}$.
 */
/**
 * nc_hicosmo_He_number_density:
 * @cosmo: a #NcHICosmo
 *
 * Calculares the Helium-4 number density $n_\mathrm{4He} = Y_p n_{b0} / m_\mathrm{4He}$
 * using nc_hicosmo_Yp_4He() $\times$ nc_hicosmo_baryon_density() / ncm_c_rest_energy_4He().
 * 
 * Returns: The baryon density $n_\mathrm{4He}$.
 */
/**
 * nc_hicosmo_H_number_density:
 * @cosmo: a #NcHICosmo
 *
 * Calculares the Hydrogen-1 number density $n_\mathrm{1H} = Y_{\mathrm{1H}p} n_{b0} / m_\mathrm{1H}$
 * using nc_hicosmo_Yp_1H() $\times$ nc_hicosmo_baryon_density() / ncm_c_rest_energy_1H().
 * 
 * Returns: The baryon density $n_\mathrm{1H}$.
 */

/**
 * nc_hicosmo_z_lss: (virtual z_lss)
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_as_drag: (virtual as_drag)
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_xb: (virtual xb)
 * @cosmo: a #NcHICosmo
 * 
 * FIXME
 * 
 * Returns: FIXME
 */

/**
 * nc_hicosmo_E2Omega_t:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * The value of the normalized Hubble function square times
 * the total dimensionless density $\Omega_t$.
 *
 * Returns: $E^2\Omega_t$.
 */
/**
 * nc_hicosmo_H:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * The value of the Hubble function in unity of $\mathrm{m}\,\mathrm{s}^{-1}\,\mathrm{kpc}^{-1}$,
 * see ncm_c_kpc().
 *
 * Returns: $H(z) \left[\mathrm{m}\,\mathrm{s}^{-1}\,\mathrm{kpc}^{-1}\right]$
 */
/**
 * nc_hicosmo_dH_dz:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_E:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the normalized Hubble function $E(z)$.
 *
 * Returns: $E(z)$.
 */
/**
 * nc_hicosmo_E2: (virtual E2)
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Normalized Hubble function squared.
 *
 * Returns: $H^2 / H_0^2$.
 */
/**
 * nc_hicosmo_Em2:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the inverse of the square normalized Hubble function.
 *
 * Returns: $E(z)^{-2}$.
 */
/**
 * nc_hicosmo_dE2_dz: (virtual dE2_dz)
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_d2E2_dz2: (virtual d2E2_dz2)
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */

/**
 * nc_hicosmo_bgp_cs2: (virtual bgp_cs2)
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * Baryon-photon plasma speed of sound squared, 
 * $$c_s^{b\gamma2} = (\dot{\rho}_b + \dot{\rho}_\gamma) / (p_b + p_\gamma).$$
 *
 * Returns: $c_s^{b\gamma2}$.
 */
/**
 * nc_hicosmo_Dc: (virtual Dc)
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_powspec: (virtual powspec)
 * @cosmo: a #NcHICosmo
 * @k: wavenumber $k$
 *
 * FIXME
 *
 * Returns: FIXME
 */

/**
 * nc_hicosmo_q:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_dec:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_wec:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_qp:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_j:
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */

/**
 * nc_hicosmo_abs_alpha:
 * @cosmo: a #NcHICosmo
 * @x: redshift variable $x = 1 + z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_x_alpha:
 * @cosmo: a #NcHICosmo
 * @alpha: redshift $\alpha$
 *
 * FIXME
 *
 * Returns: FIXME
 */

/**
 * nc_hicosmo_peek_prim:
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: (transfer none): the #NcHIPrim submodel.
 */
/**
 * nc_hicosmo_peek_reion:
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: (transfer none): the #NcHIReion submodel.
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
  NcHICosmoFunc1Z f1 = (NcHICosmoFunc1Z) obj;
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
nc_hicosmo_create_mset_func1 (NcHICosmoFunc1Z f1)
{
  return ncm_mset_func_new (&_nc_hicosmo_func1, 1, 1, f1, NULL);
}

typedef struct __NcHICosmoArrayFunc
{
  guint size;
  NcHICosmoFunc1Z f1;
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
nc_hicosmo_create_mset_arrayfunc1 (NcHICosmoFunc1Z f1, guint size)
{
  g_assert_cmpuint (size, !=, 0);
  _NcHICosmoArrayFunc *af1 = g_new (_NcHICosmoArrayFunc, 1);

  af1->size = size;
  af1->f1 = f1;

  return ncm_mset_func_new (&_nc_hicosmo_arrayfunc1, size, size, af1, g_free);
}
