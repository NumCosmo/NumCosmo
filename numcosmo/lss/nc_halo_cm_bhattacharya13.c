/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_halo_cm_bhattacharya13.c
 *
 *  Mon Dec 09 14:58:42 2024
 *  Copyright  2024  Mariana Penna-Lima <pennalima@unb.br>, Thais Mikami Ornellas <thais.ornellas@uel.br> 
 ****************************************************************************/
/*
 * nc_halo_cm_bhattacharya13.c
 * Copyright (C) 2024 Mariana Penna-Lima <pennalima@unb.br>, Thais Mikami Ornellas <thais.ornellas@uel.br>
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
 * NcHaloCMBhattacharya13
 *
 *  Class defining the Bhattacharya et al. 2013 concentration-mass relation
 *  
 * This class implements the Bhattacharya et al. 2013 concentration-mass relation.
 * FIXME include reference, equation, ranges of mass and redshift.
 * TODO: include the following conditions for the mass (parameter): 
 * Range -  M_min = 2E12 up to M_max = 2E15 
 * But, if z > 0.5, then M_max = 2E14 and, if z > 1.5, then M_max = 1E14.
 * These values are in units of Msun/h.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

/* #include "nc_enum_types.h" */
#include "lss/nc_halo_cm_bhattacharya13.h"
#include <math.h>
#include <gsl/gsl_math.h>

typedef struct _NcHaloCMBhattacharya13Private
{
  gdouble Delta;
  NcHaloMassSummaryMassDef mdef;
  NcHaloMassFunction *mfp;
  NcGrowthFunc *gf;
  gdouble (*concentration) (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z);
} NcHaloCMBhattacharya13Private;

struct _NcHaloCMBhattacharya13
{
  NcHaloMassSummary parent_instance;
};

enum
{
  PROP_0,
  PROP_GF,
  PROP_MFP,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloCMBhattacharya13, nc_halo_cm_bhattacharya13, NC_TYPE_HALO_MASS_SUMMARY);

#define VECTOR       (NCM_MODEL (hcmb))
#define LOG10M_DELTA (ncm_model_orig_param_get (VECTOR, NC_HALO_CM_BHATTACHARYA13_LOG10M_DELTA))

static void
nc_halo_cm_bhattacharya13_init (NcHaloCMBhattacharya13 *hcmb)
{
  NcHaloCMBhattacharya13Private * const self = nc_halo_cm_bhattacharya13_get_instance_private (hcmb);
 
  self->Delta = 0.0;
  self->mdef = NC_HALO_MASS_SUMMARY_MASS_DEF_LEN;
  self->mfp = NULL;
  self->gf = NULL;
  self->concentration = NULL;
}

static void
_nc_halo_cm_battacharya13_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloCMBhattacharya13 *hcmb               = NC_HALO_CM_BHATTACHARYA13 (object);
  NcHaloCMBhattacharya13Private * const self =  nc_halo_cm_bhattacharya13_get_instance_private (hcmb);

  g_return_if_fail (NC_IS_HALO_CM_BHATTACHARYA13 (object));

  switch (prop_id)
  {
    case PROP_MFP:
      self->mfp = g_value_dup_object (value);
      break;  
    case PROP_GF:
      self->gf = g_value_dup_object (value);
      break;  
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_cm_bhattacharya13_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloCMBhattacharya13 *hcmb               = NC_HALO_CM_BHATTACHARYA13 (object);
  NcHaloCMBhattacharya13Private * const self =  nc_halo_cm_bhattacharya13_get_instance_private (hcmb);

  g_return_if_fail (NC_IS_HALO_CM_BHATTACHARYA13 (object));

  switch (prop_id)
  {
    case PROP_MFP:
      g_value_set_object (value, self->mfp);
      break;
    case PROP_GF:
      g_value_set_object (value, self->gf);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_halo_cm_bhattacharya13_dispose (GObject *object)
{
  NcHaloCMBhattacharya13 *hcmb               = NC_HALO_CM_BHATTACHARYA13 (object);
  NcHaloCMBhattacharya13Private * const self =  nc_halo_cm_bhattacharya13_get_instance_private (hcmb);

  nc_halo_mass_function_clear (&self->mfp);
  nc_growth_func_clear (&self->gf);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_cm_bhattacharya13_parent_class)->dispose (object);
}

static void
_nc_halo_cm_bhattacharya13_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_halo_cm_bhattacharya13_parent_class)->finalize (object);
}

static gdouble _nc_halo_cm_bhattacharya13_mass (NcHaloMassSummary *hms);
static gdouble _nc_halo_cm_bhattacharya13_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z);
static void _nc_halo_cm_bhattacharya13_set_Delta (NcHaloMassSummary *hms, gdouble Delta);
static void _nc_halo_cm_bhattacharya13_set_mdef (NcHaloMassSummary *hms, NcHaloMassSummaryMassDef mdef);

static void
nc_halo_cm_bhattacharya13_class_init (NcHaloCMBhattacharya13Class *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcHaloMassSummaryClass *hms_class = NC_HALO_MASS_SUMMARY_CLASS (klass);
  NcmModelClass *model_class        = NCM_MODEL_CLASS (klass);

  object_class->dispose  = &_nc_halo_cm_bhattacharya13_dispose;
  object_class->finalize = &_nc_halo_cm_bhattacharya13_finalize;
  object_class->set_property = &_nc_halo_cm_battacharya13_set_property;
  object_class->get_property = &_nc_halo_cm_bhattacharya13_get_property;

  ncm_model_class_set_name_nick (model_class, "Bhattacharya et al. (2013) concentration-mass relation", "CM_BHATTACHARYA13");
  ncm_model_class_add_params (model_class, NC_HALO_CM_BHATTACHARYA13_LOCAL_SPARAM_LEN, 0, PROP_LEN);


  /**
   * NcHaloCMBhattacharya13:mass-function:
   *
   * This property keeps the halo mass function object.
   */
  g_object_class_install_property (object_class,
                                   PROP_MFP,
                                   g_param_spec_object ("mass-function",
                                                        NULL,
                                                        "Halo mass function",
                                                        NC_TYPE_HALO_MASS_FUNCTION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME
                                                        | G_PARAM_STATIC_BLURB));


  /**
   * NcHaloCMBhattacharya13:growth-function:
   *
   * This property keeps the growth function object.
   */
  g_object_class_install_property (object_class,
                                   PROP_GF,
                                   g_param_spec_object ("growth-function",
                                                        NULL,
                                                        "Growth function",
                                                        NC_TYPE_GROWTH_FUNC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME
                                                        | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloCMBhattacharya13:log10MDelta:
   *
   * Logarithm base 10 of the cluster mass $M_\Delta$ in units of solar masses $M_\odot$
   * (ncm_c_mass_solar()) within $r_\Delta$, where $\Delta$ is
   * the over-density, see Eq. \eqref{eq:mrr}.
   *
   */
  /**
   * NcHaloCMBhattacharya13:log10MDelta-fit:
   *
   * Boolean property that controls whether the parameter
   * #NcHaloMCParam:log10MDelta should be included in
   * a statistical analysis.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_CM_BHATTACHARYA13_LOG10M_DELTA, "\\log_{10}(M_{\\Delta})", "log10MDelta",
                              10.0,  17.0, 0.5,
                              NC_HALO_CM_BHATTACHARYA13_DEFAULT_PARAMS_ABSTOL, NC_HALO_CM_BHATTACHARYA13_DEFAULT_LOG10M_DELTA,
                              NCM_PARAM_TYPE_FIXED);


  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  hms_class->mass          = &_nc_halo_cm_bhattacharya13_mass;
  hms_class->concentration = &_nc_halo_cm_bhattacharya13_concentration;
  hms_class->set_Delta     = &_nc_halo_cm_bhattacharya13_set_Delta;
  hms_class->set_mdef      = &_nc_halo_cm_bhattacharya13_set_mdef;
}

static gdouble
_nc_halo_cm_bhattacharya13_mass (NcHaloMassSummary *hms)
{
  NcHaloCMBhattacharya13 *hcmb = NC_HALO_CM_BHATTACHARYA13 (hms);

  return exp10 (LOG10M_DELTA);
}

static gdouble _nc_halo_cm_bhattacharya13_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z)
{
  NcHaloCMBhattacharya13 *hcmb               = NC_HALO_CM_BHATTACHARYA13 (hms);
  NcHaloCMBhattacharya13Private * const self = nc_halo_cm_bhattacharya13_get_instance_private (hcmb);

  return self->concentration (hms, cosmo, z);
}

static void
_nc_halo_cm_bhattacharya13_set_Delta (NcHaloMassSummary *hms, gdouble Delta)
{
  NcHaloCMBhattacharya13 *hcmb               = NC_HALO_CM_BHATTACHARYA13 (hms);
  NcHaloCMBhattacharya13Private * const self = nc_halo_cm_bhattacharya13_get_instance_private (hcmb);

  if (self->mdef < NC_HALO_MASS_SUMMARY_MASS_DEF_LEN)
    if ((self->mdef != NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL) && (Delta != 200.0))
      g_error ("Bhattacharya13 concentration: mdef must be NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL or Delta must be 200.0 (mean and critical)");

  self->Delta = Delta;
}

static gdouble
_nc_halo_cm_bhattacharya13_concentration_mean (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z)
{
  NcHaloCMBhattacharya13 *hcmb               = NC_HALO_CM_BHATTACHARYA13 (hms);
  NcHaloCMBhattacharya13Private * const self = nc_halo_cm_bhattacharya13_get_instance_private (hcmb);

  nc_halo_mass_function_prepare_if_needed (self->mfp, cosmo);
  nc_growth_func_prepare_if_needed (self->gf, cosmo);

  gdouble mass  = _nc_halo_cm_bhattacharya13_mass (hms);
  gdouble lnM   = log(mass);
  gdouble sigma = nc_halo_mass_function_sigma_lnM (self->mfp, cosmo, lnM, z);
  gdouble nu    = 1.686 / sigma;
  gdouble D     = nc_growth_func_eval (self->gf, cosmo, z);

  return 9.0 * pow (D, 1.15) * pow (nu, -0.29);
}

static gdouble
_nc_halo_cm_bhattacharya13_concentration_critical (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z)
{
  NcHaloCMBhattacharya13 *hcmb               = NC_HALO_CM_BHATTACHARYA13 (hms);
  NcHaloCMBhattacharya13Private * const self = nc_halo_cm_bhattacharya13_get_instance_private (hcmb);

  nc_halo_mass_function_prepare_if_needed (self->mfp, cosmo);
  nc_growth_func_prepare_if_needed (self->gf, cosmo);

  gdouble mass  = _nc_halo_cm_bhattacharya13_mass (hms);
  gdouble lnM   = log(mass);
  gdouble sigma = nc_halo_mass_function_sigma_lnM (self->mfp, cosmo, lnM, z);
  gdouble nu    = 1.686 / sigma;
  gdouble D     = nc_growth_func_eval (self->gf, cosmo, z);

  return 5.9 * pow (D, 0.54) * pow (nu, -0.35);
}

static gdouble
_nc_halo_cm_bhattacharya13_concentration_virial (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z)
{
  NcHaloCMBhattacharya13 *hcmb               = NC_HALO_CM_BHATTACHARYA13 (hms);
  NcHaloCMBhattacharya13Private * const self = nc_halo_cm_bhattacharya13_get_instance_private (hcmb);

  nc_halo_mass_function_prepare_if_needed (self->mfp, cosmo);
  nc_growth_func_prepare_if_needed (self->gf, cosmo);

  gdouble mass  = _nc_halo_cm_bhattacharya13_mass (hms);
  gdouble lnM   = log(mass);
  gdouble sigma = nc_halo_mass_function_sigma_lnM (self->mfp, cosmo, lnM, z);
  gdouble nu    = 1.686 / sigma;
  gdouble D     = nc_growth_func_eval (self->gf, cosmo, z);
  
  return 7.7 * pow (D, 0.90) * pow (nu, -0.29);
}

static void
_nc_halo_cm_bhattacharya13_set_mdef (NcHaloMassSummary *hms, NcHaloMassSummaryMassDef mdef)
{
  NcHaloCMBhattacharya13 *hcmd               = NC_HALO_CM_BHATTACHARYA13 (hms);
  NcHaloCMBhattacharya13Private * const self = nc_halo_cm_bhattacharya13_get_instance_private (hcmd);

  if (self->Delta != 0.0)
    if ((mdef != NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL) && (self->Delta != 200.0))
      g_error ("Bhattacharya13 concentration: mdef must be NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL or Delta must be 200.0 (mean and critical)");

  self->mdef = mdef;

  switch (mdef)
  {
    case NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN:
      self->concentration = _nc_halo_cm_bhattacharya13_concentration_mean;
      break;

    case NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL:
      self->concentration = _nc_halo_cm_bhattacharya13_concentration_critical;
      break;

    case NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL:
      self->concentration = _nc_halo_cm_bhattacharya13_concentration_virial;
      break;

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }
}

/**
 * nc_halo_cm_bhattacharya13_new:
 * @mdef: a #NcHaloMassSummaryMassDef
 * @Delta: cluster threshold mass definition $\Delta$
 *
 * This function returns the #NcHaloCMBhattacharya13 implementation of
 * #NcHaloMassSummary setting #NcHaloMassSummary:mass-def to @mdef
 * and #NcHaloMassSummary:Delta to @Delta.
 *
 * Returns: a new instance of #NcHaloCMBhattacharya13.
 */
NcHaloCMBhattacharya13 *
nc_halo_cm_bhattacharya13_new (const NcHaloMassSummaryMassDef mdef, const gdouble Delta)
{
  NcHaloCMBhattacharya13 *hcmb = g_object_new (NC_TYPE_HALO_CM_BHATTACHARYA13,
                                         "mass-def", mdef,
                                         "Delta",    Delta,
                                         NULL);

  return hcmb;
}

/**
 * nc_halo_cm_bhattacharya13_ref:
 * @hcmb: a #NcHaloCMBhattacharya13
 *
 * Increase the reference of @hcmb by one.
 *
 * Returns: (transfer full): @hcmb.
 */
NcHaloCMBhattacharya13 *
nc_halo_cm_bhattacharya13_ref (NcHaloCMBhattacharya13 *hcmb)
{
  return g_object_ref (hcmb);
}

/**
 * nc_halo_cm_bhattacharya13_free:
 * @hcmb: a #NcHaloCMBhattacharya13
 *
 * Decrease the reference count of @hcmb by one.
 *
 */
void
nc_halo_cm_bhattacharya13_free (NcHaloCMBhattacharya13 *hcmb)
{
  g_object_unref (hcmb);
}

/**
 * nc_halo_cm_bhattacharya13_clear:
 * @hcmb: a #NcHaloCMBhattacharya13
 *
 * Decrease the reference count of @hcmb by one, and sets the pointer *@hcmb to
 * NULL.
 *
 */
void
nc_halo_cm_bhattacharya13_clear (NcHaloCMBhattacharya13 **hcmb)
{
  g_clear_object (hcmb);
}

