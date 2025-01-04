/***************************************************************************
 *            nc_hireion_camb.c
 *
 *  Thu December 10 11:56:38 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>>
 ****************************************************************************/
/*
 * nc_hireion_camb.c
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
 * NcHIReionCamb:
 *
 * CAMB-like reionization object.
 *
 * This object implements the reionization as done in CAMB. It is a simple
 * model that allows to set the reionization redshifts and the width of the
 * reionization windows.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hireion_camb.h"
#include "nc_hireion_camb_reparam_tau.h"
#include "math/ncm_util.h"

enum
{
  PROP_0,
  PROP_HII_HEII_REION_DELTA,
  PROP_HEIII_REION_DELTA,
  PROP_HII_HEII_REION_EXPO,
  PROP_HEIII_REIONIZED,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcHIReionCamb, nc_hireion_camb, NC_TYPE_HIREION)

#define VECTOR     (NCM_MODEL (reion))
#define HII_HEII_Z (ncm_model_orig_param_get (VECTOR, NC_HIREION_CAMB_HII_HEII_Z))
#define HEIII_Z    (ncm_model_orig_param_get (VECTOR, NC_HIREION_CAMB_HEIII_Z))

static void
nc_hireion_camb_init (NcHIReionCamb *reion_camb)
{
  reion_camb->HII_HeII_reion_delta      = 0.0;
  reion_camb->HeIII_reion_delta         = 0.0;
  reion_camb->HII_HeII_reion_expo       = 0.0;
  reion_camb->HII_HeII_reion_delta_eff  = 0.0;
  reion_camb->HII_HeII_reion_x_pow_expo = 0.0;
  reion_camb->HEII_reionized            = FALSE;
  reion_camb->fsol                      = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
  reion_camb->tau_ctrl                  = ncm_model_ctrl_new (NULL);
}

static void
nc_hireion_camb_dispose (GObject *object)
{
  NcHIReionCamb *reion_camb = NC_HIREION_CAMB (object);

  ncm_model_ctrl_clear (&reion_camb->tau_ctrl);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hireion_camb_parent_class)->dispose (object);
}

static void
nc_hireion_camb_finalize (GObject *object)
{
  NcHIReionCamb *reion_camb = NC_HIREION_CAMB (object);

  gsl_root_fsolver_free (reion_camb->fsol);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hireion_camb_parent_class)->finalize (object);
}

static void
nc_hireion_camb_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIReionCamb *reion_camb = NC_HIREION_CAMB (object);

  g_return_if_fail (NC_IS_HIREION_CAMB (object));

  switch (prop_id)
  {
    case PROP_HII_HEII_REION_DELTA:
      reion_camb->HII_HeII_reion_delta = g_value_get_double (value);
      break;
    case PROP_HEIII_REION_DELTA:
      reion_camb->HeIII_reion_delta = g_value_get_double (value);
      break;
    case PROP_HII_HEII_REION_EXPO:
      reion_camb->HII_HeII_reion_expo = g_value_get_double (value);
      break;
    case PROP_HEIII_REIONIZED:
      reion_camb->HEII_reionized = g_value_get_boolean (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hireion_camb_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIReionCamb *reion_camb = NC_HIREION_CAMB (object);

  g_return_if_fail (NC_IS_HIREION_CAMB (object));

  switch (prop_id)
  {
    case PROP_HII_HEII_REION_DELTA:
      g_value_set_double (value, reion_camb->HII_HeII_reion_delta);
      break;
    case PROP_HEIII_REION_DELTA:
      g_value_set_double (value, reion_camb->HeIII_reion_delta);
      break;
    case PROP_HII_HEII_REION_EXPO:
      g_value_set_double (value, reion_camb->HII_HeII_reion_expo);
      break;
    case PROP_HEIII_REIONIZED:
      g_value_set_boolean (value, reion_camb->HEII_reionized);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static gdouble _nc_hireion_camb_get_init_x (NcHIReion *reion, NcHICosmo *cosmo);
static gdouble _nc_hireion_camb_get_Xe (NcHIReion *reion, NcHICosmo *cosmo, const gdouble lambda, const gdouble Xe_recomb);

static void
nc_hireion_camb_class_init (NcHIReionCambClass *klass)
{
  GObjectClass *object_class  = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class  = NCM_MODEL_CLASS (klass);
  NcHIReionClass *reion_class = NC_HIREION_CLASS (klass);

  object_class->dispose  = nc_hireion_camb_dispose;
  object_class->finalize = nc_hireion_camb_finalize;

  model_class->set_property = &nc_hireion_camb_set_property;
  model_class->get_property = &nc_hireion_camb_get_property;

  ncm_model_class_set_name_nick (model_class, "Reion-CAMB", "REION_CAMB");
  ncm_model_class_add_params (model_class, NC_HIREION_CAMB_SPARAM_LEN, 0, PROP_SIZE);

  /* Set HII_HEII_Z param info */
  ncm_model_class_set_sparam (model_class, NC_HIREION_CAMB_HII_HEII_Z, "z_\\mathrm{re}", "z_re",
                              0.0, 50.0, 1.0,
                              NC_HIREION_DEFAULT_PARAMS_ABSTOL, NC_HIREION_CAMB_DEFAULT_HII_HEII_Z,
                              NCM_PARAM_TYPE_FIXED);
  /* Set HEIII_Z param info */
  ncm_model_class_set_sparam (model_class, NC_HIREION_CAMB_HEIII_Z, "z^\\mathrm{He}_\\mathrm{re}", "z_He_re",
                              0.0,  10.0, 1.0e-1,
                              NC_HIREION_DEFAULT_PARAMS_ABSTOL, NC_HIREION_CAMB_DEFAULT_HEIII_Z,
                              NCM_PARAM_TYPE_FIXED);

  g_object_class_install_property (object_class,
                                   PROP_HII_HEII_REION_DELTA,
                                   g_param_spec_double ("HII-HeII-reion-delta",
                                                        NULL,
                                                        "Window size for HII and HeII reionization",
                                                        0.0, G_MAXDOUBLE, NC_HIREION_CAMB_DEFAULT_HII_HEII_REION_DELTA,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_HEIII_REION_DELTA,
                                   g_param_spec_double ("HeIII-reion-delta",
                                                        NULL,
                                                        "Window size for HeIII reionization",
                                                        0.0, G_MAXDOUBLE, NC_HIREION_CAMB_DEFAULT_HEIII_REION_DELTA,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_HII_HEII_REION_EXPO,
                                   g_param_spec_double ("HII-HeII-reion-exponent",
                                                        NULL,
                                                        "Exponent for HII and HeII reionization transition",
                                                        0.0, G_MAXDOUBLE, NC_HIREION_CAMB_DEFAULT_HII_HEII_REION_EXPO,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HEIII_REIONIZED,
                                   g_param_spec_boolean ("HeII-reionized",
                                                         NULL,
                                                         "Whether HeIII is reionized",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  reion_class->get_init_x = &_nc_hireion_camb_get_init_x;
  reion_class->get_Xe     = &_nc_hireion_camb_get_Xe;
}

static void
_nc_hireion_camb_prepare_if_needed (NcHIReion *reion, NcHICosmo *cosmo)
{
  NcHIReionCamb *reion_camb = NC_HIREION_CAMB (reion);
  gboolean need_update      = !ncm_model_state_is_update (NCM_MODEL (reion));

  if (need_update)
  {
    const gdouble xre = 1.0 + HII_HEII_Z;

    reion_camb->HII_HeII_reion_x_pow_expo = pow (xre, reion_camb->HII_HeII_reion_expo);
    reion_camb->HII_HeII_reion_delta_eff  =
      reion_camb->HII_HeII_reion_expo *
      reion_camb->HII_HeII_reion_x_pow_expo *
      reion_camb->HII_HeII_reion_delta / xre;
    ncm_model_state_set_update (NCM_MODEL (reion));
  }
}

static gdouble
_nc_hireion_camb_get_init_x (NcHIReion *reion, NcHICosmo *cosmo)
{
  NcHIReionCamb *reion_camb = NC_HIREION_CAMB (reion);

  _nc_hireion_camb_prepare_if_needed (reion, cosmo);

  return pow (reion_camb->HII_HeII_reion_x_pow_expo - reion_camb->HII_HeII_reion_delta_eff * GSL_LOG_DBL_EPSILON, 1.0 / reion_camb->HII_HeII_reion_expo);
}

static gdouble
_nc_hireion_camb_get_Xe (NcHIReion *reion, NcHICosmo *cosmo, const gdouble lambda, const gdouble Xe_recomb)
{
  NcHIReionCamb *reion_camb = NC_HIREION_CAMB (reion);
  const gdouble XHe         = nc_hicosmo_XHe (cosmo);
  const gdouble XH_XHe      = 1.0 + XHe;
  gdouble Xe                = 0.0;

  _nc_hireion_camb_prepare_if_needed (reion, cosmo);

  if (reion_camb->HEII_reionized)
  {
    const gdouble x       = exp (-lambda);
    const gdouble x_HEIII = 1.0 + HEIII_Z;
    const gdouble arg     = 2.0 * (x - x_HEIII) / reion_camb->HeIII_reion_delta;

    if (arg < GSL_LOG_DBL_EPSILON)
      Xe = XHe;
    else
      Xe = XHe / (1.0 + exp (arg));
  }

  {
    const gdouble x_pow_expo = exp (-reion_camb->HII_HeII_reion_expo * lambda);
    const gdouble arg        = 2.0 * (x_pow_expo - reion_camb->HII_HeII_reion_x_pow_expo) / reion_camb->HII_HeII_reion_delta_eff;
    const gdouble exp_arg    = exp (arg);

    if (arg < GSL_LOG_DBL_EPSILON)
      Xe += XH_XHe;
    else
      Xe += XH_XHe / (1.0 + exp_arg);

    if (Xe_recomb != 0.0)
      Xe += Xe_recomb / (1.0 + 1.0 / exp_arg);
  }

  return Xe;
}

/**
 * nc_hireion_camb_new:
 *
 * FIXME
 *
 * Returns: a newly created #NcHIReionCamb.
 */
NcHIReionCamb *
nc_hireion_camb_new (void)
{
  NcHIReionCamb *reion_camb = g_object_new (NC_TYPE_HIREION_CAMB,
                                            NULL);

  return reion_camb;
}

typedef struct _NcHIReionCambTauToZ
{
  NcHIReionCamb *reion_camb;
  NcHICosmo *cosmo;
  gdouble tau;
} NcHIReionCambTauToZ;

static gdouble
_nc_hireion_camb_set_tau_m_tau (gdouble z_re, gpointer data)
{
  NcHIReionCambTauToZ *params = (NcHIReionCambTauToZ *) data;

  ncm_model_orig_param_set (NCM_MODEL (params->reion_camb), NC_HIREION_CAMB_HII_HEII_Z, z_re);

  return nc_hireion_get_tau (NC_HIREION (params->reion_camb), params->cosmo) / params->tau - 1.0;
}

static gdouble
_nc_hireion_camb_with_reparam_set_tau_m_tau (gdouble z_re, gpointer data)
{
  NcHIReionCambTauToZ *params = (NcHIReionCambTauToZ *) data;

  ncm_model_orig_param_set (NCM_MODEL (params->reion_camb), NC_HIREION_CAMB_HII_HEII_Z, z_re);

  return ncm_model_param_get (NCM_MODEL (params->reion_camb), NC_HIREION_CAMB_HII_HEII_Z) / params->tau - 1.0;
}

/**
 * nc_hireion_camb_calc_z_from_tau:
 * @reion_camb: a #NcHIReionCamb
 * @cosmo: a #NcHICosmo
 * @tau: reionization optical depth
 *
 * Calculates the reionization redshift from the value of the reionization
 * optical depth and the cosmological model @cosmo.
 *
 * Returns: $z_\mathrm{reion}$.
 */
gdouble
nc_hireion_camb_calc_z_from_tau (NcHIReionCamb *reion_camb, NcHICosmo *cosmo, const gdouble tau)
{
  gint status;
  gint iter = 0, max_iter = 1000000;
  gdouble z_reion = 0.0;
  gdouble z_reion_l = 0.0, z_reion_u = 0.0;
  const gdouble prec = GSL_MIN (NC_HIREION (reion_camb)->prec, 1e-1);
  gsl_function F;
  NcHIReionCambTauToZ params = {reion_camb, cosmo, tau};

  if (NC_IS_HIREION_CAMB_REPARAM_TAU (ncm_model_peek_reparam (NCM_MODEL (reion_camb))))
    F.function = &_nc_hireion_camb_with_reparam_set_tau_m_tau;
  else
    F.function = &_nc_hireion_camb_set_tau_m_tau;

  F.params = &params;

  if (F.function (0.0, &params) < 0.0)
  {
    do {
      z_reion_u += 20.0;
    } while (F.function (z_reion_u, &params) < 0.0);

    gsl_root_fsolver_set (reion_camb->fsol, &F, z_reion_l, z_reion_u);

    do {
      iter++;
      status = gsl_root_fsolver_iterate (reion_camb->fsol);

      if (status)
        g_error ("nc_hireion_camb_set_z_from_tau: Cannot find root (%s)", gsl_strerror (status));

      z_reion   = gsl_root_fsolver_root (reion_camb->fsol);
      z_reion_l = gsl_root_fsolver_x_lower (reion_camb->fsol);
      z_reion_u = gsl_root_fsolver_x_upper (reion_camb->fsol);

      status = gsl_root_test_interval (z_reion_l, z_reion_u, 0, prec);
    } while (status == GSL_CONTINUE && iter < max_iter);
  }

  return z_reion;
}

/**
 * nc_hireion_camb_set_z_from_tau:
 * @reion_camb: a #NcHIReionCamb
 * @cosmo: a #NcHICosmo
 * @tau: reionization optical depth
 *
 * Sets the reionization redshift from the value of the reionization
 * optical depth and the cosmological model @cosmo.
 *
 */
void
nc_hireion_camb_set_z_from_tau (NcHIReionCamb *reion_camb, NcHICosmo *cosmo, const gdouble tau)
{
  NcmModel *model = NCM_MODEL (reion_camb);

  if (NC_IS_HIREION_CAMB_REPARAM_TAU (ncm_model_peek_reparam (model)))
  {
    ncm_model_param_set (model, NC_HIREION_CAMB_HII_HEII_Z, tau);
  }
  else
  {
    const gdouble z_reion = nc_hireion_camb_calc_z_from_tau (reion_camb, cosmo, tau);

    ncm_model_orig_param_set (NCM_MODEL (reion_camb), NC_HIREION_CAMB_HII_HEII_Z, z_reion);

    return;
  }
}

/**
 * nc_hireion_camb_z_to_tau:
 * @reion_camb: a #NcHIReionCamb
 * @cosmo: a #NcHICosmo
 * @error: a #GError
 *
 * Changes the parametrization to use $\tau_\mathrm{reion}$ instead of $z_\mathrm{reion}$.
 *
 */
void
nc_hireion_camb_z_to_tau (NcHIReionCamb *reion_camb, NcHICosmo *cosmo, GError **error)
{
  g_return_if_fail (error == NULL || *error == NULL);
  {
    NcHIReionCambReparamTau *reparam_tau = nc_hireion_camb_reparam_tau_new (ncm_model_len (NCM_MODEL (reion_camb)), cosmo);
    NcmReparam *reparam                  = NCM_REPARAM (reparam_tau);

    ncm_model_set_reparam (NCM_MODEL (reion_camb), reparam, error);
    NCM_UTIL_ON_ERROR_RETURN (error, ncm_reparam_clear (&reparam), );

    ncm_reparam_clear (&reparam);

    return;
  }
}

