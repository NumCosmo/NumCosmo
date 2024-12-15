/***************************************************************************
 *            nc_halo_cm_klypin11.c
 *
 *  Thu Nov 28 20:52:31 2024
 *  Copyright  2024  Mariana Penna-Lima <pennalima@unb.br>,
 *  Copyright  2024  Thais Mikami Ornellas <thais.ornellas@uel.br>
 ****************************************************************************/
/*
 * nc_halo_cm_klypin11.c
 * Copyright (C) 2024 Mariana Penna-Lima <pennalima@unb.br>,
 * Copyright (C) 2024 Thais Mikami Ornellas <thais.ornellas@uel.br>
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
 * SECTION:nc_halo_cm_klypin11
 * @title: NcHaloCMKlypin11
 * @short_description: Class defining the Klypin et al. 2011 concentration-mass relation
 * @stability: Unstable
 *
 *
 * Class defining the Klypin et al. 2011 concentration-mass relation.
 * FIXME include reference and equation
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

/* #include "nc_enum_types.h" */
#include "lss/nc_halo_cm_klypin11.h"
#include <math.h>
#include <gsl/gsl_math.h>

typedef struct _NcHaloCMKlypin11Private
{
  gint placeholder;
} NcHaloCMKlypin11Private;

struct _NcHaloCMKlypin11
{
  NcHaloMassSummary parent_instance;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloCMKlypin11, nc_halo_cm_klypin11, NC_TYPE_HALO_MASS_SUMMARY);

#define VECTOR       (NCM_MODEL (hcmk))
#define LOG10M_DELTA (ncm_model_orig_param_get (VECTOR, NC_HALO_CM_KLYPIN11_LOG10M_DELTA))

static void
nc_halo_cm_klypin11_init (NcHaloCMKlypin11 *hcmk)
{
}

static void
_nc_halo_cm_klypin11_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_cm_klypin11_parent_class)->dispose (object);
}

static void
_nc_halo_cm_klypin11_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_halo_cm_klypin11_parent_class)->finalize (object);
}

static gdouble _nc_halo_cm_klypin11_mass (NcHaloMassSummary *hms);
static gdouble _nc_halo_cm_klypin11_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo);

static void
nc_halo_cm_klypin11_class_init (NcHaloCMKlypin11Class *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcHaloMassSummaryClass *hms_class = NC_HALO_MASS_SUMMARY_CLASS (klass);
  NcmModelClass *model_class        = NCM_MODEL_CLASS (klass);

  object_class->dispose  = &_nc_halo_cm_klypin11_dispose;
  object_class->finalize = &_nc_halo_cm_klypin11_finalize;

  ncm_model_class_set_name_nick (model_class, "Klypin et al. (2011) concentration-mass relation", "CM_KLYPIN11");
  ncm_model_class_add_params (model_class, NC_HALO_CM_KLYPIN11_LOCAL_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcHaloCMKlypin11:log10MDelta:
   *
   * Logarithm base 10 of the cluster mass $M_\Delta$ in units of solar masses $M_\odot$
   * (ncm_c_mass_solar()) within $r_\Delta$, where $\Delta$ is
   * the over-density, see Eq. \eqref{eq:mrr}.
   *
   */
  /**
   * NcHaloCMKlypin11:log10MDelta-fit:
   *
   * Boolean property that controls whether the parameter
   * #NcHaloMCParam:log10MDelta should be included in
   * a statistical analysis.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_CM_KLYPIN11_LOG10M_DELTA, "\\log_{10}(M_{\\Delta})", "log10MDelta",
                              10.0,  17.0, 0.5,
                              NC_HALO_CM_KLYPIN11_DEFAULT_PARAMS_ABSTOL, NC_HALO_CM_KLYPIN11_DEFAULT_LOG10M_DELTA,
                              NCM_PARAM_TYPE_FIXED);


  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  hms_class->mass          = &_nc_halo_cm_klypin11_mass;
  hms_class->concentration = &_nc_halo_cm_klypin11_concentration;
}

static gdouble
_nc_halo_cm_klypin11_mass (NcHaloMassSummary *hms)
{
  NcHaloCMKlypin11 *hcmk = NC_HALO_CM_KLYPIN11 (hms);

  return exp10 (LOG10M_DELTA);
}

static gdouble
_nc_halo_cm_klypin11_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo)
{
  gdouble mass = _nc_halo_cm_klypin11_mass (hms);
  gdouble h    = nc_hicosmo_h (cosmo);

  return 9.6 * pow (mass * h / 1.0e12, -0.075);
}

/**
 * nc_halo_cm_klypin11_new:
 * @mdef: a #NcHaloMassSummaryMassDef
 * @Delta: cluster threshold mass definition $\Delta$
 *
 * This function returns the #NcHaloCMKlypin11 implementation of
 * #NcHaloMassSummary setting #NcHaloMassSummary:mass-def to @mdef
 * and #NcHaloMassSummary:Delta to @Delta.
 *
 * Returns: a new instance of #NcHaloCMKlypin11.
 */
NcHaloCMKlypin11 *
nc_halo_cm_klypin11_new (const NcHaloMassSummaryMassDef mdef, const gdouble Delta)
{
  NcHaloCMKlypin11 *hcmk = g_object_new (NC_TYPE_HALO_CM_KLYPIN11,
                                         "mass-def", mdef,
                                         "Delta",    Delta,
                                         NULL);

  return hcmk;
}

/**
 * nc_halo_cm_klypin11_ref:
 * @hcmk: a #NcHaloCMKlypin11
 *
 * Increase the reference of @hcmk by one.
 *
 * Returns: (transfer full): @hcmk.
 */
NcHaloCMKlypin11 *
nc_halo_cm_klypin11_ref (NcHaloCMKlypin11 *hcmk)
{
  return g_object_ref (hcmk);
}

/**
 * nc_halo_cm_klypin11_free:
 * @hcmk: a #NcHaloCMKlypin11
 *
 * Decrease the reference count of @hcmk by one.
 *
 */
void
nc_halo_cm_klypin11_free (NcHaloCMKlypin11 *hcmk)
{
  g_object_unref (hcmk);
}

/**
 * nc_halo_cm_klypin11_clear:
 * @hcmk: a #NcHaloCMKlypin11
 *
 * Decrease the reference count of @hcmk by one, and sets the pointer *@hcmk to
 * NULL.
 *
 */
void
nc_halo_cm_klypin11_clear (NcHaloCMKlypin11 **hcmk)
{
  g_clear_object (hcmk);
}

