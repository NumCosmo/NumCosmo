/***************************************************************************
 *            nc_hireion.c
 *
 *  Tue December 08 14:51:37 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hireion.c
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
 * SECTION:nc_hireion
 * @title: NcHIReion
 * @short_description: Abstract class for implementing homogeneous and isotropic reionization models.
 *
 * $
 *  \newcommand{\He}{\text{He}}
 *  \newcommand{\HeI}{\text{HeI}}
 *  \newcommand{\HeII}{\text{HeII}}
 *  \newcommand{\HeIII}{\text{HeIII}}
 *  \newcommand{\Hy}{\text{H}}
 *  \newcommand{\HyI}{\text{HI}}
 *  \newcommand{\HyII}{\text{HII}}
 *  \newcommand{\e}{{\text{e}^-}}
 * $
 *
 * See [NcRecomb][NcRecomb.description] for symbol definitions.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hireion.h"
#include "nc_recomb.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_integrate.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_mset_func_list.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_integration.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_ABSTRACT_TYPE (NcHIReion, nc_hireion, NCM_TYPE_MODEL)

enum
{
  PROP_0,
  PROP_PREC,
  PROP_SIZE,
};

static void
nc_hireion_init (NcHIReion *reion)
{
  reion->prec = 0.0;
}

static void
nc_hireion_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hireion_parent_class)->finalize (object);
}

static void
nc_hireion_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIReion *reion = NC_HIREION (object);

  g_return_if_fail (NC_IS_HIREION (object));

  switch (prop_id)
  {
    case PROP_PREC:
      reion->prec = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hireion_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIReion *reion = NC_HIREION (object);

  g_return_if_fail (NC_IS_HIREION (object));

  switch (prop_id)
  {
    case PROP_PREC:
      g_value_set_double (value, reion->prec);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

NCM_MSET_MODEL_REGISTER_ID (nc_hireion, NC_TYPE_HIREION);

static gdouble _nc_hireion_get_init_x (NcHIReion *reion, NcHICosmo *cosmo);
static gdouble _nc_hireion_get_tau (NcHIReion *reion, NcHICosmo *cosmo);
static gdouble _nc_hireion_get_Xe (NcHIReion *reion, NcHICosmo *cosmo, const gdouble lambda, const gdouble Xe_recomb);

static void
nc_hireion_class_init (NcHIReionClass *klass)
{
  GObjectClass *object_class  = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class  = NCM_MODEL_CLASS (klass);
  NcHIReionClass *reion_class = NC_HIREION_CLASS (klass);

  model_class->set_property = nc_hireion_set_property;
  model_class->get_property = nc_hireion_get_property;
  object_class->finalize    = nc_hireion_finalize;

  ncm_model_class_set_name_nick (model_class, "Abstract class for H and He reionization models.", "NcHIReion");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  ncm_mset_model_register_id (model_class,
                              "NcHIReion",
                              "Homogeneous and isotropic reionization models.",
                              NULL,
                              FALSE,
                              nc_hicosmo_id ());

  g_object_class_install_property (object_class,
                                   PROP_PREC,
                                   g_param_spec_double ("prec",
                                                        NULL,
                                                        "Precision for reionization calculations",
                                                        0.0, 1.0, 1.0e-7,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_model_class_check_params_info (model_class);

  reion_class->get_init_x = &_nc_hireion_get_init_x;
  reion_class->get_tau    = &_nc_hireion_get_tau;
  reion_class->get_Xe     = &_nc_hireion_get_Xe;
}

static gdouble
_nc_hireion_get_init_x (NcHIReion *reion, NcHICosmo *cosmo)
{
  g_error ("_nc_hireion_get_init_x: error object `%s' do not implement this virtual function.",
           g_type_name (G_OBJECT_TYPE (reion)));

  return 0.0;
}

static gdouble
_nc_hireion_get_Xe (NcHIReion *reion, NcHICosmo *cosmo, const gdouble lambda, const gdouble Xe_recomb)
{
  g_error ("_nc_hireion_get_Xe: error object `%s' do not implement this virtual function.",
           g_type_name (G_OBJECT_TYPE (reion)));

  return 0.0;
}

/**
 * nc_hireion_ref:
 * @reion: a #NcHIReion
 *
 * Increses the reference count of @reion by one.
 *
 * Returns: (transfer full): @reion.
 */
NcHIReion *
nc_hireion_ref (NcHIReion *reion)
{
  return g_object_ref (reion);
}

/**
 * nc_hireion_free:
 * @reion: a #NcHIReion
 *
 * Decreases the reference count of @reion by one.
 *
 */
void
nc_hireion_free (NcHIReion *reion)
{
  g_object_unref (reion);
}

/**
 * nc_hireion_clear:
 * @reion: a #NcHIReion
 *
 * If @reion is different from NULL, decreses the reference
 * count of *@reion by one and sets *reion to NULL.
 *
 */
void
nc_hireion_clear (NcHIReion **reion)
{
  g_clear_object (reion);
}

/**
 * nc_hireion_get_init_x: (virtual get_init_x)
 * @reion: a #NcHIReion
 * @cosmo: a #NcHICosmo
 *
 * Gets the redshift ($x = 1 + z$) where the reionization begins.
 *
 * Returns: $x_\mathrm{init}$.
 */
gdouble
nc_hireion_get_init_x (NcHIReion *reion, NcHICosmo *cosmo)
{
  return NC_HIREION_GET_CLASS (reion)->get_init_x (reion, cosmo);
}

/**
 * nc_hireion_get_Xe: (virtual get_Xe)
 * @reion: a #NcHIReion
 * @cosmo: a #NcHICosmo
 * @lambda: redshift time
 * @Xe_recomb: recombination value for $X_\e$
 *
 * Gets the electron fraction from reionization $X_\e$.
 *
 * Returns: $X_\e$.
 */
gdouble
nc_hireion_get_Xe (NcHIReion *reion, NcHICosmo *cosmo, const gdouble lambda, const gdouble Xe_recomb)
{
  return NC_HIREION_GET_CLASS (reion)->get_Xe (reion, cosmo, lambda, Xe_recomb);
}

typedef struct _NcHIReionTauInt
{
  NcHIReion *reion;
  NcHICosmo *cosmo;
} NcHIReionTauInt;

static gdouble
_nc_hireion_dtau_dlambda (gdouble lambda, gpointer data)
{
  NcHIReionTauInt *tau_int = (NcHIReionTauInt *) data;

  return nc_hireion_get_Xe (tau_int->reion, tau_int->cosmo, lambda, 0.0) * nc_recomb_dtau_dlambda_Xe (tau_int->cosmo, lambda);
}

static gdouble
_nc_hireion_get_tau (NcHIReion *reion, NcHICosmo *cosmo)
{
  gdouble result                = 0.0;
  gdouble error                 = 0.0;
  const gdouble lambda_i        = -log (nc_hireion_get_init_x (reion, cosmo));
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  NcHIReionTauInt tau_int       = {reion, cosmo};
  gsl_function F;

  F.function = &_nc_hireion_dtau_dlambda;
  F.params   = &tau_int;

  gsl_integration_qag (&F, lambda_i, 0.0, 1.0, reion->prec, NCM_INTEGRAL_PARTITION, 6, *w, &result, &error);

  ncm_memory_pool_return (w);

  return -result;
}

/**
 * nc_hireion_get_tau: (virtual get_tau)
 * @reion: a #NcHIReion
 * @cosmo: a #NcHICosmo
 *
 * Calculates the reionization optical depth $\tau_\mathrm{reion}$.
 *
 * Returns: $\tau_\mathrm{reion}$.
 */
gdouble
nc_hireion_get_tau (NcHIReion *reion, NcHICosmo *cosmo)
{
  return NC_HIREION_GET_CLASS (reion)->get_tau (reion, cosmo);
}

static void
_nc_hireion_flist_tau (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHIReion *reion = NC_HIREION (ncm_mset_peek (mset, nc_hireion_id ()));

  res[0] = nc_hireion_get_tau (reion, cosmo);
}

void
_nc_hireion_register_functions (void)
{
  ncm_mset_func_list_register ("tau", "\\tau", "NcHIReion", "Reionization tau", G_TYPE_NONE, _nc_hireion_flist_tau, 0, 1);
}

