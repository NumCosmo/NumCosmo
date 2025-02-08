/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_halo_cm_duffy08.c
 *
 *  Thu Dec 05 09:42:15 2024
 *  Copyright  2024  Mariana Penna-Lima <pennalima@unb.br>, Thais Mikami Ornellas <thais.ornellas@uel.br>
 ****************************************************************************/
/*
 * nc_halo_cm_duffy08.c
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
 * SECTION:nc_halo_cm_duffy08
 * @title: NcHaloCMDuffy08
 * @short_description: Class defining the Duffy et al. 2008 concentration-mass relation
 * @stability: Unstable
 *
 *
 * Class defining the Duffy et al. 2008 concentration-mass relation.
 * FIXME include reference and equation
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

/* #include "nc_enum_types.h" */
#include "lss/nc_halo_cm_duffy08.h"
#include <math.h>
#include <gsl/gsl_math.h>

typedef struct _NcHaloCMDuffy08Private
{
  gdouble Delta;
  NcHaloMassSummaryMassDef mdef;

  gdouble (*concentration) (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z);
} NcHaloCMDuffy08Private;

struct _NcHaloCMDuffy08
{
  NcHaloMassSummary parent_instance;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloCMDuffy08, nc_halo_cm_duffy08, NC_TYPE_HALO_MASS_SUMMARY);

#define VECTOR       (NCM_MODEL (hcmd))
#define LOG10M_DELTA (ncm_model_orig_param_get (VECTOR, NC_HALO_CM_DUFFY08_LOG10M_DELTA))

static void
nc_halo_cm_duffy08_init (NcHaloCMDuffy08 *hcmd)
{
  NcHaloCMDuffy08Private * const self = nc_halo_cm_duffy08_get_instance_private (hcmd);

  self->Delta = 0.0;
  self->mdef  = NC_HALO_MASS_SUMMARY_MASS_DEF_LEN;

  self->concentration = NULL;
}

static void
_nc_halo_cm_duffy08_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_cm_duffy08_parent_class)->dispose (object);
}

static void
_nc_halo_cm_duffy08_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_halo_cm_duffy08_parent_class)->finalize (object);
}

static gdouble _nc_halo_cm_duffy08_mass (NcHaloMassSummary *hms);
static gdouble _nc_halo_cm_duffy08_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z);
static void _nc_halo_cm_duffy08_set_Delta (NcHaloMassSummary *hms, gdouble Delta);
static void _nc_halo_cm_duffy08_set_mdef (NcHaloMassSummary *hms, NcHaloMassSummaryMassDef mdef);

static void
nc_halo_cm_duffy08_class_init (NcHaloCMDuffy08Class *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcHaloMassSummaryClass *hms_class = NC_HALO_MASS_SUMMARY_CLASS (klass);
  NcmModelClass *model_class        = NCM_MODEL_CLASS (klass);

  object_class->dispose  = &_nc_halo_cm_duffy08_dispose;
  object_class->finalize = &_nc_halo_cm_duffy08_finalize;

  ncm_model_class_set_name_nick (model_class, "Duffy et al. (2008) concentration-mass relation", "CM_DUFFY08");
  ncm_model_class_add_params (model_class, NC_HALO_CM_DUFFY08_LOCAL_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcHaloCMDuffy08:log10MDelta:
   *
   * Logarithm base 10 of the cluster mass $M_\Delta$ in units of solar masses $M_\odot$
   * (ncm_c_mass_solar()) within $r_\Delta$, where $\Delta$ is
   * the over-density, see Eq. \eqref{eq:mrr}.
   *
   */
  /**
   * NcHaloCMDuffy08:log10MDelta-fit:
   *
   * Boolean property that controls whether the parameter
   * #NcHaloMCParam:log10MDelta should be included in
   * a statistical analysis.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_CM_DUFFY08_LOG10M_DELTA, "\\log_{10}(M_{\\Delta})", "log10MDelta",
                              10.0,  17.0, 0.5,
                              NC_HALO_CM_DUFFY08_DEFAULT_PARAMS_ABSTOL, NC_HALO_CM_DUFFY08_DEFAULT_LOG10M_DELTA,
                              NCM_PARAM_TYPE_FIXED);


  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  hms_class->mass          = &_nc_halo_cm_duffy08_mass;
  hms_class->concentration = &_nc_halo_cm_duffy08_concentration;
  hms_class->set_Delta     = &_nc_halo_cm_duffy08_set_Delta;
  hms_class->set_mdef      = &_nc_halo_cm_duffy08_set_mdef;
}

static gdouble
_nc_halo_cm_duffy08_mass (NcHaloMassSummary *hms)
{
  NcHaloCMDuffy08 *hcmd = NC_HALO_CM_DUFFY08 (hms);

  return exp10 (LOG10M_DELTA);
}

static gdouble
_nc_halo_cm_duffy08_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z)
{
  NcHaloCMDuffy08 *hcmd               = NC_HALO_CM_DUFFY08 (hms);
  NcHaloCMDuffy08Private * const self = nc_halo_cm_duffy08_get_instance_private (hcmd);

  return self->concentration (hms, cosmo, z);
}

static void
_nc_halo_cm_duffy08_set_Delta (NcHaloMassSummary *hms, gdouble Delta)
{
  NcHaloCMDuffy08 *hcmd               = NC_HALO_CM_DUFFY08 (hms);
  NcHaloCMDuffy08Private * const self = nc_halo_cm_duffy08_get_instance_private (hcmd);

  if (self->mdef < NC_HALO_MASS_SUMMARY_MASS_DEF_LEN)
    if ((self->mdef != NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL) && (Delta != 200.0))
      g_error ("Duffy08 concentration: mdef must be NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL or Delta must be 200.0 (mean and critical)");

  self->Delta = Delta;
}

static gdouble
_nc_halo_cm_duffy08_concentration_mean (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z)
{
  gdouble mass = _nc_halo_cm_duffy08_mass (hms);
  gdouble h    = nc_hicosmo_h (cosmo);

  return 10.14 * pow (mass * h / 2.0e12, -0.081) * pow (1.0 + z, -1.01);
}

static gdouble
_nc_halo_cm_duffy08_concentration_critical (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z)
{
  gdouble mass = _nc_halo_cm_duffy08_mass (hms);
  gdouble h    = nc_hicosmo_h (cosmo);

  return 5.71 * pow (mass * h / 2.0e12, -0.084) * pow (1.0 + z, -0.47);
}

static gdouble
_nc_halo_cm_duffy08_concentration_virial (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z)
{
  gdouble mass = _nc_halo_cm_duffy08_mass (hms);
  gdouble h    = nc_hicosmo_h (cosmo);

  return 7.85 * pow (mass * h / 2.0e12, -0.081) * pow (1.0 + z, -0.71);
}

static void
_nc_halo_cm_duffy08_set_mdef (NcHaloMassSummary *hms, NcHaloMassSummaryMassDef mdef)
{
  NcHaloCMDuffy08 *hcmd               = NC_HALO_CM_DUFFY08 (hms);
  NcHaloCMDuffy08Private * const self = nc_halo_cm_duffy08_get_instance_private (hcmd);

  if (self->Delta != 0.0)
    if ((mdef != NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL) && (self->Delta != 200.0))
      g_error ("Duffy08 concentration: mdef must be NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL or Delta must be 200.0 (mean and critical)");

  self->mdef = mdef;

  switch (mdef)
  {
    case NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN:
      self->concentration = _nc_halo_cm_duffy08_concentration_mean;
      break;

    case NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL:
      self->concentration = _nc_halo_cm_duffy08_concentration_critical;
      break;

    case NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL:
      self->concentration = _nc_halo_cm_duffy08_concentration_virial;
      break;

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }
}

/**
 * nc_halo_cm_duffy08_new:
 * @mdef: a #NcHaloMassSummaryMassDef
 * @Delta: cluster threshold mass definition $\Delta$
 *
 * This function returns the #NcHaloCMDuffy08 implementation of
 * #NcHaloMassSummary setting #NcHaloMassSummary:mass-def to @mdef
 * and #NcHaloMassSummary:Delta to @Delta.
 *
 * Returns: a new instance of #NcHaloCMDuffy08.
 */
NcHaloCMDuffy08 *
nc_halo_cm_duffy08_new (const NcHaloMassSummaryMassDef mdef, const gdouble Delta)
{
  NcHaloCMDuffy08 *hcmd = g_object_new (NC_TYPE_HALO_CM_DUFFY08,
                                        "mass-def", mdef,
                                        "Delta",    Delta,
                                        NULL);

  return hcmd;
}

/**
 * nc_halo_cm_duffy08_ref:
 * @hcmd: a #NcHaloCMDuffy08
 *
 * Increase the reference of @hcmd by one.
 *
 * Returns: (transfer full): @hcmd.
 */
NcHaloCMDuffy08 *
nc_halo_cm_duffy08_ref (NcHaloCMDuffy08 *hcmd)
{
  return g_object_ref (hcmd);
}

/**
 * nc_halo_cm_duffy08_free:
 * @hcmd: a #NcHaloCMDuffy08
 *
 * Decrease the reference count of @hcmd by one.
 *
 */
void
nc_halo_cm_duffy08_free (NcHaloCMDuffy08 *hcmd)
{
  g_object_unref (hcmd);
}

/**
 * nc_halo_cm_duffy08_clear:
 * @hcmd: a #NcHaloCMDuffy08
 *
 * Decrease the reference count of @hcmd by one, and sets the pointer *@hcmd to
 * NULL.
 *
 */
void
nc_halo_cm_duffy08_clear (NcHaloCMDuffy08 **hcmd)
{
  g_clear_object (hcmd);
}

