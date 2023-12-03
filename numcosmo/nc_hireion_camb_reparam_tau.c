/***************************************************************************
 *            nc_hireion_camb_reparam_tau.c
 *
 *  Tue December 15 04:31:34 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hireion_camb_reparam_tau.c
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
 * SECTION:nc_hireion_camb_reparam_tau
 * @title: NcHIReionCambReparamTau
 * @short_description: CAMB reionization reparametrization $z_\mathrm{reion} \to \tau_\mathrm{reion}$.
 *
 * Object implementing a reparametrization for CAMB reionization. It changes
 * $z_\mathrm{reion} \to \tau_\mathrm{reion}$.
 *
 */

#include "nc_hireion_camb_reparam_tau.h"
#include "nc_hireion_camb.h"

G_DEFINE_TYPE (NcHIReionCambReparamTau, nc_hireion_camb_reparam_tau, NCM_TYPE_REPARAM)

enum
{
  PROP_0,
  PROP_COSMO,
};


static void
nc_hireion_camb_reparam_tau_init (NcHIReionCambReparamTau *reparam_tau)
{
  reparam_tau->ctrl = ncm_model_ctrl_new (NULL);
}

static void
nc_hireion_camb_reparam_tau_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hireion_camb_reparam_tau_parent_class)->constructed (object);
  {
    NcHIReionCambReparamTau *reparam_tau = NC_HIREION_CAMB_REPARAM_TAU (object);
    ncm_reparam_set_param_desc_full (NCM_REPARAM (reparam_tau), NC_HIREION_CAMB_HII_HEII_Z,
                                     "tau_reion","\\tau_\\mathrm{reion}", 0.0, 1.0, 1.0e-2,
                                     NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, 0.0, NCM_PARAM_TYPE_FIXED);
    NCM_REPARAM (reparam_tau)->compat_type = NC_TYPE_HIREION_CAMB;
  }
}

static void
nc_hireion_camb_reparam_tau_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIReionCambReparamTau *reparam_tau = NC_HIREION_CAMB_REPARAM_TAU (object);
  g_return_if_fail (NC_IS_HIREION_CAMB_REPARAM_TAU (object));

  switch (prop_id)
  {
    case PROP_COSMO:
      ncm_model_ctrl_update (reparam_tau->ctrl, NCM_MODEL (g_value_get_object (value)));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hireion_camb_reparam_tau_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIReionCambReparamTau *reparam_tau = NC_HIREION_CAMB_REPARAM_TAU (object);
  g_return_if_fail (NC_IS_HIREION_CAMB_REPARAM_TAU (object));

  switch (prop_id)
  {
    case PROP_COSMO:
      g_value_take_object (value, ncm_model_ctrl_get_model (reparam_tau->ctrl));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hireion_camb_reparam_tau_dispose (GObject *object)
{
  NcHIReionCambReparamTau *reparam_tau = NC_HIREION_CAMB_REPARAM_TAU (object);

  ncm_model_ctrl_clear (&reparam_tau->ctrl);  

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hireion_camb_reparam_tau_parent_class)->dispose (object);
}

static void
nc_hireion_camb_reparam_tau_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hireion_camb_reparam_tau_parent_class)->finalize (object);
}

static gboolean _nc_hireion_camb_reparam_tau_old2new (NcmReparam *reparam, NcmModel *model);
static gboolean _nc_hireion_camb_reparam_tau_new2old (NcmReparam *reparam, NcmModel *model);
static gboolean _nc_hireion_camb_reparam_tau_jac (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac);

static void
nc_hireion_camb_reparam_tau_class_init (NcHIReionCambReparamTauClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmReparamClass *reparam_class = NCM_REPARAM_CLASS (klass);

  object_class->constructed  = &nc_hireion_camb_reparam_tau_constructed;
  object_class->set_property = &nc_hireion_camb_reparam_tau_set_property;
  object_class->get_property = &nc_hireion_camb_reparam_tau_get_property;
  object_class->dispose      = &nc_hireion_camb_reparam_tau_dispose;
  object_class->finalize     = &nc_hireion_camb_reparam_tau_finalize;

  g_object_class_install_property (object_class,
                                   PROP_COSMO,
                                   g_param_spec_object ("cosmo",
                                                        NULL,
                                                        "Cosmological model used to transform tau <=> z",
                                                        NC_TYPE_HICOSMO,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  
  reparam_class->old2new = &_nc_hireion_camb_reparam_tau_old2new;
  reparam_class->new2old = &_nc_hireion_camb_reparam_tau_new2old;
  reparam_class->jac     = &_nc_hireion_camb_reparam_tau_jac;
}

static gboolean
_nc_hireion_camb_reparam_tau_old2new (NcmReparam *reparam, NcmModel *model)
{
  NcHICosmo *cosmo = NC_HICOSMO (ncm_model_ctrl_get_model (NC_HIREION_CAMB_REPARAM_TAU (reparam)->ctrl));
  NcmVector *params = ncm_model_orig_params_peek_vector (model);

  ncm_vector_memcpy (reparam->new_params, params);
  {
    const gdouble tau = nc_hireion_get_tau (NC_HIREION (model), cosmo);
    ncm_vector_set (reparam->new_params, NC_HIREION_CAMB_HII_HEII_Z, tau);
  }

  nc_hicosmo_free (cosmo);
  return TRUE;
}

static gboolean
_nc_hireion_camb_reparam_tau_new2old (NcmReparam *reparam, NcmModel *model)
{
  NcHICosmo *cosmo = NC_HICOSMO (ncm_model_ctrl_get_model (NC_HIREION_CAMB_REPARAM_TAU (reparam)->ctrl));
  NcmVector *params = ncm_model_orig_params_peek_vector (model);

  ncm_vector_memcpy (params, reparam->new_params);

  {
    const gdouble tau     = ncm_vector_get (reparam->new_params, NC_HIREION_CAMB_HII_HEII_Z);
    const gdouble z_reion = nc_hireion_camb_calc_z_from_tau (NC_HIREION_CAMB (model), cosmo, tau);
    
    ncm_vector_set (params, NC_HIREION_CAMB_HII_HEII_Z, z_reion);
  }

  nc_hicosmo_free (cosmo);
  return TRUE;
}

static gboolean
_nc_hireion_camb_reparam_tau_jac (NcmReparam *reparam, NcmModel *model, NcmMatrix *jac)
{
  g_assert_not_reached ();
}

/**
 * nc_hireion_camb_reparam_tau_new: (constructor)
 * @length: number of parameters
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcHIReionCambReparamTau
 */
NcHIReionCambReparamTau *
nc_hireion_camb_reparam_tau_new (guint length, NcHICosmo *cosmo)
{
  NcHIReionCambReparamTau *reparam_tau = g_object_new (NC_TYPE_HIREION_CAMB_REPARAM_TAU,
                                                       "length", length,
                                                       "cosmo", cosmo,
                                                       NULL);
  return reparam_tau;
}
