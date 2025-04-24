/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_halo_cm_prada12.c
 *
 *  Mon Dec 16 15:01:55 2024
 *  Copyright  2024  Mariana Penna-Lima <pennalima@unb.br>, Thais Mikami Ornellas <thais.ornellas@uel.br> 
 ****************************************************************************/
/*
 * nc_halo_cm_prada12.c
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
 * SECTION:nc_halo_cm_prada12
 * @title: NcHaloCMPrada12
 * @short_description: Class defining the Prada et al. 2012 concentration-mass relation
 * @stability: Unstable
 *
 *
 * Class defining the Prada et al. 2012 concentration-mass relation.
 * FIXME include reference and equation
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

/* #include "nc_enum_types.h" */
#include "lss/nc_halo_cm_prada12.h"
#include <math.h>
#include <gsl/gsl_math.h>

typedef struct _NcHaloCMPrada12Private
{
  gdouble Delta;
  NcHaloMassSummaryMassDef mdef;

  gdouble (*concentration) (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z);
} NcHaloCMPrada12Private;

struct _NcHaloCMPrada12
{
  NcHaloMassSummary parent_instance;
  NcHaloMassFunction *mfp;
  NcmPowspecFilter *psf;
  NcHICosmoDE *cosmo_de;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloCMPrada12, nc_halo_cm_prada12, NC_TYPE_HALO_MASS_SUMMARY);

#define VECTOR       (NCM_MODEL (hcmp))
#define LOG10M_DELTA (ncm_model_orig_param_get (VECTOR, NC_HALO_CM_PRADA12_LOG10M_DELTA))

static void
nc_halo_cm_prada12_init (NcHaloCMPrada12 *hcmp)
{
  NcHaloCMPrada12Private * const self = nc_halo_cm_prada12_get_instance_private (hcmp);

  self->Delta = 0.0;
  self->mdef  = NC_HALO_MASS_SUMMARY_MASS_DEF_LEN;

  self->concentration = NULL;
}

static void
_nc_halo_cm_prada12_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_cm_prada12_parent_class)->dispose (object);
}

static void
_nc_halo_cm_prada12_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_halo_cm_prada12_parent_class)->finalize (object);
}

static gdouble _nc_halo_cm_prada12_mass (NcHaloMassSummary *hms);
static gdouble _nc_halo_cm_prada12_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z);
static void _nc_halo_cm_prada12_set_Delta (NcHaloMassSummary *hms, gdouble Delta);
static void _nc_halo_cm_prada12_set_mdef (NcHaloMassSummary *hms, NcHaloMassSummaryMassDef mdef);

static void
nc_halo_cm_prada12_class_init (NcHaloCMPrada12Class *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcHaloMassSummaryClass *hms_class = NC_HALO_MASS_SUMMARY_CLASS (klass);
  NcmModelClass *model_class        = NCM_MODEL_CLASS (klass);

  object_class->dispose  = &_nc_halo_cm_prada12_dispose;
  object_class->finalize = &_nc_halo_cm_prada12_finalize;

  ncm_model_class_set_name_nick (model_class, "Prada et al. (2012) concentration-mass relation", "CM_PRADA12");
  ncm_model_class_add_params (model_class, NC_HALO_CM_PRADA12_LOCAL_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcHaloCMPrada12:log10MDelta:
   *
   * Logarithm base 10 of the cluster mass $M_\Delta$ in units of solar masses $M_\odot$
   * (ncm_c_mass_solar()) within $r_\Delta$, where $\Delta$ is
   * the over-density, see Eq. \eqref{eq:mrr}.
   *
   */
  /**
   * NcHaloCMPrada12:log10MDelta-fit:
   *
   * Boolean property that controls whether the parameter
   * #NcHaloMCParam:log10MDelta should be included in
   * a statistical analysis.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_CM_PRADA12_LOG10M_DELTA, "\\log_{10}(M_{\\Delta})", "log10MDelta",
                              10.0,  17.0, 0.5,
                              NC_HALO_CM_PRADA12_DEFAULT_PARAMS_ABSTOL, NC_HALO_CM_PRADA12_DEFAULT_LOG10M_DELTA,
                              NCM_PARAM_TYPE_FIXED);


  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  hms_class->mass          = &_nc_halo_cm_prada12_mass;
  hms_class->concentration = &_nc_halo_cm_prada12_concentration;
  hms_class->set_Delta     = &_nc_halo_cm_prada12_set_Delta;
  hms_class->set_mdef      = &_nc_halo_cm_prada12_set_mdef;
}

static gdouble
_nc_halo_cm_prada12_mass (NcHaloMassSummary *hms)
{
  NcHaloCMPrada12 *hcmp = NC_HALO_CM_PRADA12 (hms);

  return exp10 (LOG10M_DELTA);
}

static gdouble
_nc_halo_cm_prada12_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z)
{
  NcHaloCMPrada12 *hcmp               = NC_HALO_CM_PRADA12 (hms);
  NcHaloCMPrada12Private * const self = nc_halo_cm_prada12_get_instance_private (hcmp);

  return self->concentration (hms, cosmo, z);
}

static void
_nc_halo_cm_prada12_set_Delta (NcHaloMassSummary *hms, gdouble Delta)
{
  NcHaloCMPrada12 *hcmp               = NC_HALO_CM_PRADA12 (hms);
  NcHaloCMPrada12Private * const self = nc_halo_cm_prada12_get_instance_private (hcmp);

  if (self->mdef < NC_HALO_MASS_SUMMARY_MASS_DEF_LEN)
    if ((self->mdef != NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL) || (Delta != 200.0))
      g_error ("Prada12 concentration: mdef must be NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL AND Delta must be 200.0");

  self->Delta = Delta;
}

static gdouble
_cmin (const gdouble x)
{
  return 3.681 + (5.033 - 3.681) * (1 / M_PI * atan(6.948 * (x - 0.424)) + 0.5);
}

static gdouble
_sigmin (const gdouble x)
{
  return 1.047 + (1.646 - 1.047) * (1 / M_PI * atan(7.386 * (x - 0.526)) + 0.5);
}

static gdouble
_nc_halo_cm_prada12_concentration_critical (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z)
{
  
  NcHaloCMPrada12 *hcmp = NC_HALO_CM_PRADA12 (hms);
  gdouble mass          = _nc_halo_cm_prada12_mass (hms);
  gdouble lnM           = log(mass);
  gdouble R             = exp(nc_halo_mass_function_lnM_to_lnR (hcmp->mfp, cosmo, lnM));
  gdouble sigma         = ncm_powspec_filter_eval_sigma (hcmp->psf, z, R);
  gdouble a             = 1.0 / (1.0 + z);
  gdouble Ode           = nc_hicosmo_de_E2Omega_de (hcmp->cosmo_de, z);
  gdouble Om0           = nc_hicosmo_Omega_m0 (cosmo);
  
  gdouble x             = pow (Ode / Om0, (1.0 / 3.0) * a);
  gdouble B0, B1, sig, C;

  B0 = _cmin(x) / _cmin(1.393);
  B1 = _sigmin(x) / _sigmin(1.393);
  sig = B1 * sigma;
  C = 2.881 * (pow(sig / 1.257, 1.022 + 1.0) * exp(pow(0.060 / sig, 2.00)));

  return B0 * C;
}

static void
_nc_halo_cm_prada12_set_mdef (NcHaloMassSummary *hms, NcHaloMassSummaryMassDef mdef)
{
  NcHaloCMPrada12 *hcmp               = NC_HALO_CM_PRADA12 (hms);
  NcHaloCMPrada12Private * const self = nc_halo_cm_prada12_get_instance_private (hcmp);

  if (self->Delta != 0.0)
    if ((mdef != NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL) || (self->Delta != 200.0))
      g_error ("Prada12 concentration: mdef must be NC_HALO_MASS_SUMMARY_MASS_DEF_CRITIAL and Delta must be 200.0");

  self->mdef = mdef;

  switch (mdef)
  {
    case NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL:
      self->concentration = _nc_halo_cm_prada12_concentration_critical;
      break;

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }
}

/**
 * nc_halo_cm_prada12_new:
 * @mdef: a #NcHaloMassSummaryMassDef
 * @Delta: cluster threshold mass definition $\Delta$
 *
 * This function returns the #NcHaloCMPrada12 implementation of
 * #NcHaloMassSummary setting #NcHaloMassSummary:mass-def to @mdef
 * and #NcHaloMassSummary:Delta to @Delta.
 *
 * Returns: a new instance of #NcHaloCMPrada12.
 */
NcHaloCMPrada12 *
nc_halo_cm_prada12_new (const NcHaloMassSummaryMassDef mdef, const gdouble Delta)
{
  NcHaloCMPrada12 *hcmp = g_object_new (NC_TYPE_HALO_CM_PRADA12,
                                         "mass-def", mdef,
                                         "Delta",    Delta,
                                         NULL);

  return hcmp;
}

/**
 * nc_halo_cm_prada12_ref:
 * @hcmp: a #NcHaloCMPrada12
 *
 * Increase the reference of @hcmp by one.
 *
 * Returns: (transfer full): @hcmp.
 */
NcHaloCMPrada12 *
nc_halo_cm_prada12_ref (NcHaloCMPrada12 *hcmp)
{
  return g_object_ref (hcmp);
}

/**
 * nc_halo_cm_prada12_free:
 * @hcmp: a #NcHaloCMPrada12
 *
 * Decrease the reference count of @hcmp by one.
 *
 */
void
nc_halo_cm_prada12_free (NcHaloCMPrada12 *hcmp)
{
  g_object_unref (hcmp);
}

/**
 * nc_halo_cm_prada12_clear:
 * @hcmp: a #NcHaloCMPrada12
 *
 * Decrease the reference count of @hcmp by one, and sets the pointer *@hcmp to
 * NULL.
 *
 */
void
nc_halo_cm_prada12_clear (NcHaloCMPrada12 **hcmp)
{
  g_clear_object (hcmp);
}

