/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_halo_cm_diemer15.c
 *
 *  Wed Jun 11 16:13:15 2025
 *  Copyright  2025  Mariana Penna-Lima <pennalima@unb.br>, Thais Mikami Ornellas <thais.ornellas@uel.br>
 ****************************************************************************/
/*
 * nc_halo_cm_diemer15.c
 * Copyright (C) 2025 Mariana Penna-Lima <pennalima@unb.br>, Thais Mikami Ornellas <thais.ornellas@uel.br>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcHaloCMDiemer15:
 *
 * Class defining the Diemer & Kravtsov 2015 concentration-mass relation.
 * FIXME include reference and equation
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

/* #include "nc_enum_types.h" */
#include "lss/nc_halo_cm_diemer15.h"
#include "nc_powspec_ml_transfer.h"
#include "lss/nc_transfer_func_eh_no_baryon.h"
#include <math.h>
#include <gsl/gsl_math.h>

typedef struct _NcHaloCMDiemer15Private
{
  gdouble Delta;
  NcHaloMassSummaryMassDef mdef;
  NcHaloMassFunction *mfp;
  NcmPowspec *powspec;

  gdouble (*concentration) (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z);
} NcHaloCMDiemer15Private;

struct _NcHaloCMDiemer15
{
  NcHaloMassSummary parent_instance;
};

enum
{
  PROP_0,
  PROP_MFP,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloCMDiemer15, nc_halo_cm_diemer15, NC_TYPE_HALO_MASS_SUMMARY);

#define VECTOR       (NCM_MODEL (hcmdk))
#define LOG10M_DELTA (ncm_model_orig_param_get (VECTOR, NC_HALO_CM_DIEMER15_LOG10M_DELTA))

static void
nc_halo_cm_diemer15_init (NcHaloCMDiemer15 *hcmdk)
{
  NcHaloCMDiemer15Private * const self = nc_halo_cm_diemer15_get_instance_private (hcmdk);
  NcTransferFunc *tf                   = nc_transfer_func_eh_no_baryon_new ();
  NcPowspecMLTransfer *ps_mlt          = nc_powspec_ml_transfer_new (tf);

  self->Delta   = 0.0;
  self->mdef    = NC_HALO_MASS_SUMMARY_MASS_DEF_LEN;
  self->mfp     = NULL;
  self->powspec = NCM_POWSPEC (ps_mlt);

  self->concentration = NULL;

  nc_transfer_func_free (tf);
}

static void
_nc_halo_cm_diemer15_dispose (GObject *object)
{
  NcHaloCMDiemer15 *hcmdk              = NC_HALO_CM_DIEMER15 (object);
  NcHaloCMDiemer15Private * const self = nc_halo_cm_diemer15_get_instance_private (hcmdk);

  nc_halo_mass_function_clear (&self->mfp);
  ncm_powspec_clear (&self->powspec);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_cm_diemer15_parent_class)->dispose (object);
}

static void
_nc_halo_cm_diemer15_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_halo_cm_diemer15_parent_class)->finalize (object);
}

static gdouble _nc_halo_cm_diemer15_mass (NcHaloMassSummary *hms);
static gdouble _nc_halo_cm_diemer15_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z);
static void _nc_halo_cm_diemer15_set_Delta (NcHaloMassSummary *hms, gdouble Delta);
static void _nc_halo_cm_diemer15_set_mdef (NcHaloMassSummary *hms, NcHaloMassSummaryMassDef mdef);

static void
_nc_halo_cm_diemer15_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloCMDiemer15 *hcmdk              = NC_HALO_CM_DIEMER15 (object);
  NcHaloCMDiemer15Private * const self = nc_halo_cm_diemer15_get_instance_private (hcmdk);

  g_return_if_fail (NC_IS_HALO_CM_DIEMER15 (object));

  switch (prop_id)
  {
    case PROP_MFP:
      nc_halo_mass_function_clear (&self->mfp);
      self->mfp = g_value_dup_object (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_halo_cm_diemer15_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloCMDiemer15 *hcmdk              = NC_HALO_CM_DIEMER15 (object);
  NcHaloCMDiemer15Private * const self = nc_halo_cm_diemer15_get_instance_private (hcmdk);

  g_return_if_fail (NC_IS_HALO_CM_DIEMER15 (object));

  switch (prop_id)
  {
    case PROP_MFP:
      g_value_set_object (value, self->mfp);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_halo_cm_diemer15_class_init (NcHaloCMDiemer15Class *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcHaloMassSummaryClass *hms_class = NC_HALO_MASS_SUMMARY_CLASS (klass);
  NcmModelClass *model_class        = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_halo_cm_diemer15_set_property;
  model_class->get_property = &_nc_halo_cm_diemer15_get_property;
  object_class->dispose     = &_nc_halo_cm_diemer15_dispose;
  object_class->finalize    = &_nc_halo_cm_diemer15_finalize;

  ncm_model_class_set_name_nick (model_class, "Diemer & Kravtsov (2015) concentration-mass relation", "CM_DIEMER15");
  ncm_model_class_add_params (model_class, NC_HALO_CM_DIEMER15_LOCAL_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcHaloCMDiemer15:mass-function:
   *
   * The halo mass function object.
   */
  g_object_class_install_property (object_class,
                                   PROP_MFP,
                                   g_param_spec_object ("mass-function",
                                                        NULL,
                                                        "Halo mass function",
                                                        NC_TYPE_HALO_MASS_FUNCTION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloCMDiemer15:log10MDelta:
   *
   * Logarithm base 10 of the cluster mass $M_\Delta$ in units of solar masses $M_\odot$
   * (ncm_c_mass_solar()) within $r_\Delta$, where $\Delta$ is
   * the over-density, see Eq. \eqref{eq:mrr}.
   *
   */
  /**
   * NcHaloCMDiemer15:log10MDelta-fit:
   *
   * Boolean property that controls whether the parameter
   * #NcHaloMCParam:log10MDelta should be included in
   * a statistical analysis.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_CM_DIEMER15_LOG10M_DELTA, "\\log_{10}(M_{\\Delta})", "log10MDelta",
                              10.0,  17.0, 0.5,
                              NC_HALO_CM_DIEMER15_DEFAULT_PARAMS_ABSTOL, NC_HALO_CM_DIEMER15_DEFAULT_LOG10M_DELTA,
                              NCM_PARAM_TYPE_FIXED);


  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  hms_class->mass          = &_nc_halo_cm_diemer15_mass;
  hms_class->concentration = &_nc_halo_cm_diemer15_concentration;
  hms_class->set_Delta     = &_nc_halo_cm_diemer15_set_Delta;
  hms_class->set_mdef      = &_nc_halo_cm_diemer15_set_mdef;
}

static gdouble
_nc_halo_cm_diemer15_mass (NcHaloMassSummary *hms)
{
  NcHaloCMDiemer15 *hcmdk = NC_HALO_CM_DIEMER15 (hms);

  return exp10 (LOG10M_DELTA);
}

static gdouble
_nc_halo_cm_diemer15_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z)
{
  NcHaloCMDiemer15 *hcmdk              = NC_HALO_CM_DIEMER15 (hms);
  NcHaloCMDiemer15Private * const self = nc_halo_cm_diemer15_get_instance_private (hcmdk);

  return self->concentration (hms, cosmo, z);
}

static void
_nc_halo_cm_diemer15_set_Delta (NcHaloMassSummary *hms, gdouble Delta)
{
  NcHaloCMDiemer15 *hcmdk              = NC_HALO_CM_DIEMER15 (hms);
  NcHaloCMDiemer15Private * const self = nc_halo_cm_diemer15_get_instance_private (hcmdk);

  if (self->mdef < NC_HALO_MASS_SUMMARY_MASS_DEF_LEN)
    if ((self->mdef != NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL) || (Delta != 200.0))
      g_error ("Diemer15 concentration: mdef must be NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL and Delta must be 200.0");

  self->Delta = Delta;
}

static gdouble
_dlnpk_dlnk (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z, gdouble k_R)
{
  NcHaloCMDiemer15 *hcmdk              = NC_HALO_CM_DIEMER15 (hms);
  NcHaloCMDiemer15Private * const self = nc_halo_cm_diemer15_get_instance_private (hcmdk);
  const gdouble Pk                     = ncm_powspec_eval (self->powspec, NCM_MODEL (cosmo), z, k_R);
  const gdouble dPk_dk                 = ncm_powspec_deriv_k (self->powspec, NCM_MODEL (cosmo), z, k_R);

  return k_R * dPk_dk / Pk;
}

static gdouble
_nc_halo_cm_diemer15_concentration_critical (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z)
{
  NcHaloCMDiemer15 *hcmdk              = NC_HALO_CM_DIEMER15 (hms);
  NcHaloCMDiemer15Private * const self = nc_halo_cm_diemer15_get_instance_private (hcmdk);

  nc_halo_mass_function_prepare_if_needed (self->mfp, cosmo);
  ncm_powspec_prepare_if_needed (self->powspec, NCM_MODEL (cosmo));

  const gdouble mass      = _nc_halo_cm_diemer15_mass (hms);
  const gdouble lnM       = log (mass);
  const gdouble sigma     = nc_halo_mass_function_sigma_lnM (self->mfp, cosmo, lnM, z);
  const gdouble R         = exp (nc_halo_mass_function_lnM_to_lnR (self->mfp, cosmo, lnM));
  const gdouble kappa     = 1.0; /* TODO: colossus uses kappa = 1.0 but paper uses kappa = 0.69 (original parameters) */
  const gdouble k_R       = 2.0 * M_PI * kappa / R;
  const gdouble dlnk_dlnk = _dlnpk_dlnk (hms, cosmo, z, k_R);
  const gdouble cmin      = 6.58 + 1.27 * dlnk_dlnk;
  const gdouble nmin      = 7.28 + 1.56 * dlnk_dlnk;
  const gdouble nu        = 1.686 / sigma;

  return cmin * 0.5 * (pow ((nmin / nu), 1.08) + pow (nu / nmin, 1.77));
}

static void
_nc_halo_cm_diemer15_set_mdef (NcHaloMassSummary *hms, NcHaloMassSummaryMassDef mdef)
{
  NcHaloCMDiemer15 *hcmdk              = NC_HALO_CM_DIEMER15 (hms);
  NcHaloCMDiemer15Private * const self = nc_halo_cm_diemer15_get_instance_private (hcmdk);

  if (self->Delta != 0.0)
    if ((mdef != NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL) || (self->Delta != 200.0))
      g_error ("Diemer15 concentration: mdef must be NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL and Delta must be 200.0");

  self->mdef = mdef;

  switch (mdef)
  {
    case NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL:
      self->concentration = _nc_halo_cm_diemer15_concentration_critical;
      break;

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }
}

/**
 * nc_halo_cm_diemer15_new:
 * @mdef: a #NcHaloMassSummaryMassDef
 * @Delta: cluster threshold mass definition $\Delta$
 * @mfp: a #NcHaloMassFunction
 *
 * This function returns the #NcHaloCMDiemer15 implementation of
 * #NcHaloMassSummary setting #NcHaloMassSummary:mass-def to @mdef
 * and #NcHaloMassSummary:Delta to @Delta.
 *
 * Returns: a new instance of #NcHaloCMDiemer15.
 */
NcHaloCMDiemer15 *
nc_halo_cm_diemer15_new (const NcHaloMassSummaryMassDef mdef, const gdouble Delta, NcHaloMassFunction *mfp)
{
  NcHaloCMDiemer15 *hcmdk = g_object_new (NC_TYPE_HALO_CM_DIEMER15,
                                          "mass-def", mdef,
                                          "Delta",    Delta,
                                          "mass-function", mfp,
                                          NULL);

  return hcmdk;
}

/**
 * nc_halo_cm_diemer15_ref:
 * @hcmdk: a #NcHaloCMDiemer15
 *
 * Increase the reference of @hcmdk by one.
 *
 * Returns: (transfer full): @hcmdk.
 */
NcHaloCMDiemer15 *
nc_halo_cm_diemer15_ref (NcHaloCMDiemer15 *hcmdk)
{
  return g_object_ref (hcmdk);
}

/**
 * nc_halo_cm_diemer15_free:
 * @hcmdk: a #NcHaloCMDiemer15
 *
 * Decrease the reference count of @hcmdk by one.
 *
 */
void
nc_halo_cm_diemer15_free (NcHaloCMDiemer15 *hcmdk)
{
  g_object_unref (hcmdk);
}

/**
 * nc_halo_cm_diemer15_clear:
 * @hcmdk: a #NcHaloCMDiemer15
 *
 * Decrease the reference count of @hcmdk by one, and sets the pointer *@hcmdk to
 * NULL.
 *
 */
void
nc_halo_cm_diemer15_clear (NcHaloCMDiemer15 **hcmdk)
{
  g_clear_object (hcmdk);
}

