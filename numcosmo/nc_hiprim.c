/***************************************************************************
 *            nc_hiprim.c
 *
 *  Tue October 27 12:12:41 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hiprim.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:nc_hiprim
 * @title: NcHIPrim
 * @short_description: Abstract class for implementing homogeneous and isotropic primordial cosmological models.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_serialize.h"
#include "nc_hicosmo.h"
#include "nc_hiprim.h"
#include "math/ncm_cfg.h"
#include "math/ncm_mset_func_list.h"

enum
{
  PROP_0,
  PROP_K_PIVOT,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcHIPrim, nc_hiprim, NCM_TYPE_MODEL);

static void
nc_hiprim_init (NcHIPrim *prim)
{
  prim->k_pivot   = 0.0;
  prim->lnk_pivot = 0.0;
}

static void
nc_hiprim_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hiprim_parent_class)->finalize (object);
}

static void
nc_hiprim_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPrim *prim = NC_HIPRIM (object);
  g_return_if_fail (NC_IS_HIPRIM (object));

  switch (prop_id)
  {
    case PROP_K_PIVOT:
      nc_hiprim_set_k_pivot (prim, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hiprim_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIPrim *prim = NC_HIPRIM (object);
  g_return_if_fail (NC_IS_HIPRIM (object));

  switch (prop_id)
  {
    case PROP_K_PIVOT:
      g_value_set_double (value, nc_hiprim_get_k_pivot (prim));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

NCM_MSET_MODEL_REGISTER_ID (nc_hiprim, NC_TYPE_HIPRIM);

static void
nc_hiprim_class_init (NcHIPrimClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property  = nc_hiprim_set_property;
  model_class->get_property  = nc_hiprim_get_property;

  object_class->finalize     = nc_hiprim_finalize;

  ncm_model_class_set_name_nick (model_class, "Homogeneous and isotropic primordial cosmological models.", "NcHIPrim");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  ncm_mset_model_register_id (model_class,
                              "NcHIPrim",
                              "Homogeneous and isotropic primordial cosmological models.",
                              NULL,
                              FALSE,
                              nc_hicosmo_id ());

  ncm_model_class_check_params_info (model_class);


  g_object_class_install_property (object_class,
                                   PROP_K_PIVOT,
                                   g_param_spec_double ("k-pivot",
                                                        NULL,
                                                        "Pivotal value of k",
                                                        G_MINDOUBLE, G_MAXDOUBLE, NC_HIPRIM_DEFAULT_K_PIVOT,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_hiprim_new_from_name:
 * @parent_type: #GType of the parent model
 * @prim_name: name of the primordial spcetrum model
 *
 * This function instantiates a new object of type #NcHIPrim.
 *
 * Returns: A new #NcHIPrim
 */
NcHIPrim *
nc_hiprim_new_from_name (GType parent_type, gchar *prim_name)
{
  GObject *obj = ncm_serialize_global_from_string (prim_name);
  GType model_type = G_OBJECT_TYPE (obj);

  if (!g_type_is_a (model_type, parent_type))
    g_error ("nc_hiprim_new_from_name: NcHIPrim %s do not descend from %s.", prim_name, g_type_name (parent_type));
  return NC_HIPRIM (obj);
}

/**
 * nc_hiprim_ref:
 * @prim: a #NcHIPrim
 *
 * Increases the reference count of @prim by one.
 *
 * Returns: (transfer full): @prim
 */
NcHIPrim *
nc_hiprim_ref (NcHIPrim *prim)
{
  return g_object_ref (prim);
}

/**
 * nc_hiprim_free:
 * @prim: a #NcHIPrim
 *
 * Atomically decreases the reference count of @prim by one. If the reference count drops to 0,
 * all memory allocated by @prim is released.
 *
 */
void
nc_hiprim_free (NcHIPrim *prim)
{
  g_object_unref (prim);
}

/**
 * nc_hiprim_clear:
 * @prim: a #NcHIPrim
 *
 * The reference count of @prim is decreased and the pointer is set to NULL.
 *
 */
void
nc_hiprim_clear (NcHIPrim **prim)
{
  g_clear_object (prim);
}

static void
_nc_hiprim_log_all_models_go (GType model_type, guint n)
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
      _nc_hiprim_log_all_models_go (models[i], n + 2);

    g_free (modelsc);
  }
  g_free (models);
}

/**
 * nc_hiprim_log_all_models:
 * @parent: #GType of the parent model
 *
 * Logs all models descending from @parent.
 *
 */
void
nc_hiprim_log_all_models (GType parent)
{
  g_message ("# Registred NcHIPrim:%s are:\n", g_type_name (parent));
  _nc_hiprim_log_all_models_go (parent, 0);
}

/**
 * nc_hiprim_set_k_pivot:
 * @prim: a #NcHIPrim
 * @k_pivot: pivotal $k$ in units of $1/\mathrm{Mpc}$
 *
 * Sets @k_pivot to the respective property. 
 *
 */
void
nc_hiprim_set_k_pivot (NcHIPrim *prim, gdouble k_pivot)
{
  prim->k_pivot   = k_pivot;
  prim->lnk_pivot = log (k_pivot);
  g_assert (gsl_finite (prim->lnk_pivot));
}

/**
 * nc_hiprim_get_k_pivot:
 * @prim: a #NcHIPrim
 *
 * Gets the value of the pivotal $k$.
 *
 * Returns: pivotal $k$ in units of $1/\mathrm{Mpc}$
 */
gdouble
nc_hiprim_get_k_pivot (NcHIPrim *prim)
{
  return prim->k_pivot;
}

/**
 * nc_hiprim_get_lnk_pivot:
 * @prim: a #NcHIPrim
 *
 * Gets the value of the pivotal $k$.
 *
 * Returns: $\ln(k_\mathrm{pivot}\mathrm{Mpc})$
 */
gdouble
nc_hiprim_get_lnk_pivot (NcHIPrim *prim)
{
  return log (prim->k_pivot);
}

/**
 * nc_hiprim_lnSA_powspec_lnk: (virtual lnSA_powspec_lnk)
 * @prim: a #NcHIPrim
 * @lnk: $\ln(k\mathrm{Mpc})$
 *
 * Gets the natural logarithm of the scalar adiabatic power spectrum as a
 * function of $\ln(k\mathrm{Mpc})$
 *
 * Return: $\log(P_{SA})$
 */
/**
 * nc_hiprim_lnT_powspec_lnk: (virtual lnT_powspec_lnk)
 * @prim: a #NcHIPrim
 * @lnk: $\ln(k\mathrm{Mpc})$
 *
 * Gets the natural logarithm of the tensor power spectrum as a
 * function of $\ln(k\mathrm{Mpc})$
 *
 * Return: $\log(P_{T})$
 */
/**
 * nc_hiprim_SA_powspec_k:
 * @prim: a #NcHIPrim
 * @k: $k$ in units of $1/\mathrm{Mpc}$
 *
 * Gets the scalar adiabatic power spectrum as a function of $k$.
 *
 * Return: $P_{SA}$
 */
/**
 * nc_hiprim_T_powspec_k:
 * @prim: a #NcHIPrim
 * @k: $k$ in units of $1/\mathrm{Mpc}$
 *
 * Gets the tensor power spectrum as a function of $k$.
 *
 * Return: $P_{T}$
 */
/**
 * nc_hiprim_SA_Ampl:
 * @prim: a #NcHIPrim
 *
 * Gets the scalar adiabatic power spectrum amplitude,
 * i.e., $P_{SA}(k_\mathrm{pivot})$.
 *
 * Return: $P_{SA}(k_\mathrm{pivot})$
 */
/**
 * nc_hiprim_T_Ampl:
 * @prim: a #NcHIPrim
 *
 * Gets the tensor power spectrum amplitude,
 * i.e., $P_{T}(k_\mathrm{pivot})$.
 *
 * Return: $P_{T}(k_\mathrm{pivot})$
 */
/**
 * nc_hiprim_T_SA_ratio:
 * @prim: a #NcHIPrim
 *
 * Gets the tensor-to-scalar ratio.
 *
 * Return: $P_{T}(k_\mathrm{pivot})/P_{SA}(k_\mathrm{pivot})$
 */

/**
 * nc_hiprim_set_lnSA_powspec_lnk_impl: (skip)
 * @model_class: FIXME
 * @f: (scope notified): FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HIPRIM,NcHIPrim,nc_hiprim,NcHIPrimFunc1,lnSA_powspec_lnk)
/**
 * nc_hiprim_set_lnT_powspec_lnk_impl: (skip)
 * @model_class: FIXME
 * @f: (scope notified): FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HIPRIM,NcHIPrim,nc_hiprim,NcHIPrimFunc1,lnT_powspec_lnk)

#define _NC_HIPRIM_FUNC1_TO_FLIST(fname) \
static void _nc_hiprim_flist_##fname (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *res) \
{ \
 NcHIPrim *prim = NC_HIPRIM (ncm_mset_peek (mset, nc_hiprim_id ())); \
 g_assert (prim != NULL); \
 res[0] = nc_hiprim_##fname (prim, x[0]); \
}

_NC_HIPRIM_FUNC1_TO_FLIST (lnSA_powspec_lnk)
_NC_HIPRIM_FUNC1_TO_FLIST (lnT_powspec_lnk)
_NC_HIPRIM_FUNC1_TO_FLIST (SA_powspec_k)
_NC_HIPRIM_FUNC1_TO_FLIST (T_powspec_k)

void
_nc_hiprim_register_functions (void)
{
  ncm_mset_func_list_register ("lnSA_powspec_lnk", "\\ln(P_\\mathrm{SA})", "NcHIPrim", "Logarithm of the SA power spectrum", G_TYPE_NONE, _nc_hiprim_flist_lnSA_powspec_lnk, 1, 1);
  ncm_mset_func_list_register ("lnT_powspec_lnk",  "\\ln(P_\\mathrm{T})",  "NcHIPrim", "Logarithm of the T power spectrum",  G_TYPE_NONE, _nc_hiprim_flist_lnT_powspec_lnk,  1, 1);
  ncm_mset_func_list_register ("SA_powspec_k",     "P_\\mathrm{SA}",       "NcHIPrim", "SA power spectrum",                  G_TYPE_NONE, _nc_hiprim_flist_SA_powspec_k,     1, 1);
  ncm_mset_func_list_register ("T_powspec_k",      "P_\\mathrm{T}",        "NcHIPrim", "T power spectrum",                   G_TYPE_NONE, _nc_hiprim_flist_T_powspec_k,      1, 1);
}
