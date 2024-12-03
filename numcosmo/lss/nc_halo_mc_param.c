/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_halo_mc_param.c
 *
 *  Sat Oct 12 11:51:53 2024
 *  Copyright  2024  Mariana Penna-Lima
 *  <pennalima@unb.br>
 ****************************************************************************/
/*
 * nc_halo_mc_param.c
 * Copyright (C) 2024 Mariana Penna-Lima <pennalima@unb.br>
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
 * SECTION:nc_halo_mc_param
 * @title: NcHaloMCParam
 * @short_description: Class defining mass and concentration as parameters
 * @stability: Unstable
 *
 *
 * Class defining the halo mass and concentration as parameters for the halo mass
 * density profile.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

/* #include "nc_enum_types.h" */
#include "lss/nc_halo_mc_param.h"
#include <math.h>
#include <gsl/gsl_math.h>

typedef struct _NcHaloMCParamPrivate
{
  gint placeholder;
} NcHaloMCParamPrivate;

struct _NcHaloMCParam
{
  NcHaloMassSummary parent_instance;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloMCParam, nc_halo_mc_param, NC_TYPE_HALO_MASS_SUMMARY);

#define VECTOR       (NCM_MODEL (hmcp))
#define LOG10M_DELTA (ncm_model_orig_param_get (VECTOR, NC_HALO_MC_PARAM_LOG10M_DELTA))
#define C_DELTA      (ncm_model_orig_param_get (VECTOR, NC_HALO_MC_PARAM_C_DELTA))

static void
nc_halo_mc_param_init (NcHaloMCParam *hmcp)
{
}

static void
_nc_halo_mc_param_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_mc_param_parent_class)->dispose (object);
}

static void
_nc_halo_mc_param_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_halo_mc_param_parent_class)->finalize (object);
}

static gdouble _nc_halo_mc_param_mass (NcHaloMassSummary *hms);
static gdouble _nc_halo_mc_param_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo);

static void
nc_halo_mc_param_class_init (NcHaloMCParamClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcHaloMassSummaryClass *hms_class = NC_HALO_MASS_SUMMARY_CLASS (klass);
  NcmModelClass *model_class        = NCM_MODEL_CLASS (klass);

  object_class->dispose  = &_nc_halo_mc_param_dispose;
  object_class->finalize = &_nc_halo_mc_param_finalize;

  ncm_model_class_set_name_nick (model_class, "Mass and concentration as parameters", "MC_PARAM");
  ncm_model_class_add_params (model_class, NC_HALO_MC_PARAM_LOCAL_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcHaloMCParam:log10MDelta:
   *
   * Logarithm base 10 of the cluster mass $M_\Delta$ in units of solar masses $M_\odot$
   * (ncm_c_mass_solar()) within $r_\Delta$, where $\Delta$ is
   * the over-density, see Eq. \eqref{eq:mrr}.
   *
   */
  /**
   * NcHaloMCParam:log10MDelta-fit:
   *
   * Boolean property that controls whether the parameter
   * #NcHaloMCParam:log10MDelta should be included in
   * a statistical analysis.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_MC_PARAM_LOG10M_DELTA, "\\log_{10}(M_{\\Delta})", "log10MDelta",
                              10.0,  17.0, 0.5,
                              NC_HALO_MC_PARAM_DEFAULT_PARAMS_ABSTOL, NC_HALO_MC_PARAM_DEFAULT_LOG10M_DELTA,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcHaloMCParam:cDelta:
   *
   * Concentration parameter, $c_\Delta$, see Eq \eqref{def:cDelta}.
   *
   */
  /**
   * NcHaloMCParam:cDelta-fit:
   *
   * Boolean property that controls whether the parameter
   * #NcHaloMCParam:cDelta should be included in
   * a statistical analysis.
   *
   */
  ncm_model_class_set_sparam (model_class, NC_HALO_MC_PARAM_C_DELTA, "c_{\\Delta}", "cDelta",
                              1.0e-1,  30.0, 1.0e-1,
                              NC_HALO_MC_PARAM_DEFAULT_PARAMS_ABSTOL, NC_HALO_MC_PARAM_DEFAULT_C_DELTA,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  hms_class->mass          = &_nc_halo_mc_param_mass;
  hms_class->concentration = &_nc_halo_mc_param_concentration;
}

static gdouble
_nc_halo_mc_param_mass (NcHaloMassSummary *hms)
{
  NcHaloMCParam *hmcp = NC_HALO_MC_PARAM (hms);

  return exp10 (LOG10M_DELTA);
}

static gdouble
_nc_halo_mc_param_concentration (NcHaloMassSummary *hms, NcHICosmo *cosmo)
{
  NcHaloMCParam *hmcp = NC_HALO_MC_PARAM (hms);

  return C_DELTA;
}

/**
 * nc_halo_mc_param_new:
 * @mdef: a #NcHaloMassSummaryMassDef
 * @Delta: cluster threshold mass definition $\Delta$
 *
 * This function returns the #NcHaloMCParam implementation of
 * #NcHaloMassSummary setting #NcHaloMassSummary:mass-def to @mdef
 * and #NcHaloMassSummary:Delta to @Delta.
 *
 * Returns: a new instance of #NcHaloMCParam.
 */
NcHaloMCParam *
nc_halo_mc_param_new (const NcHaloMassSummaryMassDef mdef, const gdouble Delta)
{
  NcHaloMCParam *hmcp = g_object_new (NC_TYPE_HALO_MC_PARAM,
                                      "mass-def", mdef,
                                      "Delta",    Delta,
                                      NULL);

  return hmcp;
}

/**
 * nc_halo_mc_param_ref:
 * @hmcp: a #NcHaloMCParam
 *
 * Increase the reference of @hmcp by one.
 *
 * Returns: (transfer full): @hmcp.
 */
NcHaloMCParam *
nc_halo_mc_param_ref (NcHaloMCParam *hmcp)
{
  return g_object_ref (hmcp);
}

/**
 * nc_halo_mc_param_free:
 * @hmcp: a #NcHaloMCParam
 *
 * Decrease the reference count of @hmcp by one.
 *
 */
void
nc_halo_mc_param_free (NcHaloMCParam *hmcp)
{
  g_object_unref (hmcp);
}

/**
 * nc_halo_mc_param_clear:
 * @hmcp: a #NcHaloMCParam
 *
 * Decrease the reference count of @hmcp by one, and sets the pointer *@hmcp to
 * NULL.
 *
 */
void
nc_halo_mc_param_clear (NcHaloMCParam **hmcp)
{
  g_clear_object (hmcp);
}

