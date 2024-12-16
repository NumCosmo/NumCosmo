/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_halo_cm_dutton14.c
 *
 *  Mon Dec 16 13:06:37 2024
 *  Copyright  2024  Mariana Penna-Lima <pennalima@unb.br>, Thais Mikami Ornellas <thais.ornellas@uel.br> 
 ****************************************************************************/
/*
 * nc_halo_cm_dutton14.c
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
 * SECTION:nc_halo_cm_dutton14
 * @title: NcHaloCMDutton14
 * @short_description: Class defining the Dutton & Macciò (2014) concentration-mass relation
 * @stability: Unstable
 *
 *
 * Class defining the Dutton & Macciò  2014 concentration-mass relation.
 * FIXME include reference and equation
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

/* #include "nc_enum_types.h" */
#include "lss/nc_halo_cm_dutton14.h"
#include <math.h>
#include <gsl/gsl_math.h>

typedef struct _NcHaloCMDutton14Private
{
  gint placeholder;
} NcHaloCMDutton14Private;

struct _NcHaloCMDutton14
{
  NcHaloMassSummary parent_instance;
  NcHaloMassSummaryMassDef mdef;
  gdouble z;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloCMDutton14, nc_halo_cm_dutton14, NC_TYPE_HALO_MASS_SUMMARY);

#define VECTOR       (NCM_MODEL (hcmdm))
#define LOG10M_DELTA (ncm_model_orig_param_get (VECTOR, NC_HALO_CM_DUTTON14_LOG10M_DELTA))

static void
nc_halo_cm_dutton14_init (NcHaloCMDutton14 *hcmdm)
{
}

static void
_nc_halo_cm_dutton14_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_cm_dutton14_parent_class)->dispose (object);
}

static void
_nc_halo_cm_dutton14_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_halo_cm_dutton14_parent_class)->finalize (object);
}

static gdouble _nc_halo_cm_dutton14_mass (NcHaloMassSummary *hms);
static gdouble _nc_halo_cm_dutton14_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z);

static void
nc_halo_cm_dutton14_class_init (NcHaloCMDutton14Class *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcHaloMassSummaryClass *hms_class = NC_HALO_MASS_SUMMARY_CLASS (klass);
  NcmModelClass *model_class        = NCM_MODEL_CLASS (klass);

  object_class->dispose  = &_nc_halo_cm_dutton14_dispose;
  object_class->finalize = &_nc_halo_cm_dutton14_finalize;

  ncm_model_class_set_name_nick (model_class, "Dutton & Macciò (2014) concentration-mass relation", "CM_DUTTON14");
  ncm_model_class_add_params (model_class, NC_HALO_CM_DUTTON14_LOCAL_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcHaloCMDutton14:log10MDelta:
   *
   * Logarithm base 10 of the cluster mass $M_\Delta$ in units of solar masses $M_\odot$
   * (ncm_c_mass_solar()) within $r_\Delta$, where $\Delta$ is
   * the over-density, see Eq. \eqref{eq:mrr}.
   *
   */
  /**
   * NcHaloCMDutton14:log10MDelta-fit:
   *
   * Boolean property that controls whether the parameter
   * #NcHaloMCParam:log10MDelta should be included in
   * a statistical analysis.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_CM_DUTTON14_LOG10M_DELTA, "\\log_{10}(M_{\\Delta})", "log10MDelta",
                              10.0,  17.0, 0.5,
                              NC_HALO_CM_DUTTON14_DEFAULT_PARAMS_ABSTOL, NC_HALO_CM_DUTTON14_DEFAULT_LOG10M_DELTA,
                              NCM_PARAM_TYPE_FIXED);


  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  hms_class->mass          = &_nc_halo_cm_dutton14_mass;
  hms_class->concentration = &_nc_halo_cm_dutton14_concentration;
}

static gdouble
_nc_halo_cm_dutton14_mass (NcHaloMassSummary *hms)
{
  NcHaloCMDutton14 *hcmdm = NC_HALO_CM_DUTTON14 (hms);

  return exp10 (LOG10M_DELTA);
}

static gdouble
_nc_halo_cm_dutton14_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo, gdouble z)
{
  NcHaloCMDutton14 *hcmdm = NC_HALO_CM_DUTTON14 (hms);
  gdouble mass            = _nc_halo_cm_dutton14_mass (hms);
  gdouble h               = nc_hicosmo_h (cosmo);
  gdouble Delta           = nc_halo_mass_summary_Delta (hms, cosmo, z);
  gdouble a, b;

  hcmdm->mdef             = NC_HALO_MASS_SUMMARY_MASS_DEF_LEN;

  switch (hcmdm->mdef)
  {
    case NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL:

      if (Delta == 200.0){
        a = 0.520 + (0.905 - 0.520) * exp(-0.617 * pow (z, 1.21));
        b = -0.101 + 0.026 * z;
        return pow (10, a + b * log10 (mass / 1.0e12 * h));
      }

      else{
        g_assert_not_reached();
        return 0.0;    
      } 

    case NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL:

        a = 0.537 + (1.025 - 0.537) * exp(-0.718 * pow (z, 1.08));
        b = -0.097 + 0.024 * z;
        return pow (10, a + b * log10 (mass / 1.0e12 * h));

    default:                   
      g_assert_not_reached ();

      return 0.0;
  }

}

/**
 * nc_halo_cm_dutton14_new:
 * @mdef: a #NcHaloMassSummaryMassDef
 * @Delta: cluster threshold mass definition $\Delta$
 *
 * This function returns the #NcHaloCMDutton14 implementation of
 * #NcHaloMassSummary setting #NcHaloMassSummary:mass-def to @mdef
 * and #NcHaloMassSummary:Delta to @Delta.
 *
 * Returns: a new instance of #NcHaloCMDutton14.
 */
NcHaloCMDutton14 *
nc_halo_cm_dutton14_new (const NcHaloMassSummaryMassDef mdef, const gdouble Delta)
{
  NcHaloCMDutton14 *hcmdm = g_object_new (NC_TYPE_HALO_CM_DUTTON14,
                                         "mass-def", mdef,
                                         "Delta",    Delta,
                                         NULL);

  return hcmdm;
}

/**
 * nc_halo_cm_dutton14_ref:
 * @hcmdm: a #NcHaloCMDutton14
 *
 * Increase the reference of @hcmdm by one.
 *
 * Returns: (transfer full): @hcmdm.
 */
NcHaloCMDutton14 *
nc_halo_cm_dutton14_ref (NcHaloCMDutton14 *hcmdm)
{
  return g_object_ref (hcmdm);
}

/**
 * nc_halo_cm_dutton14_free:
 * @hcmdm: a #NcHaloCMDutton14
 *
 * Decrease the reference count of @hcmdm by one.
 *
 */
void
nc_halo_cm_dutton14_free (NcHaloCMDutton14 *hcmdm)
{
  g_object_unref (hcmdm);
}

/**
 * nc_halo_cm_dutton14_clear:
 * @hcmdm: a #NcHaloCMDutton14
 *
 * Decrease the reference count of @hcmdm by one, and sets the pointer *@hcmdm to
 * NULL.
 *
 */
void
nc_halo_cm_dutton14_clear (NcHaloCMDutton14 **hcmdm)
{
  g_clear_object (hcmdm);
}

