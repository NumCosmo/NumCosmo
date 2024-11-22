/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_halo_mass_summary.c
 *
 *  Thu Oct 10 14:12:27 2024
 *  Copyright  2024  Mariana Penna-Lima
 *  <pennalima@unb.br>
 ****************************************************************************/
/*
 * nc_halo_mass_summary.c
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
 * SECTION: nc_halo_mass_summary
 * @title: NcHaloMassSummary
 * @short_description: Class describing halo mass summary
 * @stability:
 *
 * This class describes a halo mass summary, i.e. the mass definition
 * of a halo, mass-concentration relationship, etc.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include <math.h>
#include <gsl/gsl_math.h>
#include "nc_enum_types.h"
#include "lss/nc_halo_density_profile.h"
#include "lss/nc_halo_mass_summary.h"

typedef struct _NcHaloMassSummaryPrivate
{
  NcHaloMassSummaryMassDef mdef;
  gdouble z;
  gdouble Delta;
} NcHaloMassSummaryPrivate;

enum
{
  PROP_0,
  PROP_MDEF,
  PROP_DELTA,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloMassSummary, nc_halo_mass_summary, NCM_TYPE_MODEL);

static void
nc_halo_mass_summary_init (NcHaloMassSummary *hms)
{
  NcHaloMassSummaryPrivate * const self = nc_halo_mass_summary_get_instance_private (hms);

  self->mdef  = NC_HALO_MASS_SUMMARY_MASS_DEF_LEN;
  self->z     = 0.0;
  self->Delta = 0.0;
}

static void
_nc_halo_mass_summary_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloMassSummary *hms                = NC_HALO_MASS_SUMMARY (object);
  NcHaloMassSummaryPrivate * const self = nc_halo_mass_summary_get_instance_private (hms);

  g_return_if_fail (NC_IS_HALO_MASS_SUMMARY (object));

  switch (prop_id)
  {
    case PROP_MDEF:
      self->mdef = g_value_get_enum (value);
      break;
    case PROP_DELTA:
      self->Delta = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_halo_mass_summary_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloMassSummary *hms                = NC_HALO_MASS_SUMMARY (object);
  NcHaloMassSummaryPrivate * const self = nc_halo_mass_summary_get_instance_private (hms);

  g_return_if_fail (NC_IS_HALO_MASS_SUMMARY (object));

  switch (prop_id)
  {
    case PROP_MDEF:
      g_value_set_enum (value, self->mdef);
      break;
    case PROP_DELTA:
      g_value_set_double (value, self->Delta);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_halo_mass_summary_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_mass_summary_parent_class)->dispose (object);
}

static void
_nc_halo_mass_summary_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_halo_mass_summary_parent_class)->finalize (object);
}

/* LCOV_EXCL_START */
static gdouble
_nc_halo_mass_summary_mass (NcHaloMassSummary *hms)
{
  g_error ("_nc_halo_mass_summary_mass: method not implemented.");

  return 0.0;
}

static gdouble
_nc_halo_mass_summary_concentration (NcHaloMassSummary *hms)
{
  g_error ("_nc_halo_mass_summary_concentration: method not implemented.");

  return 0.0;
}

/* LCOV_EXCL_STOP */

NCM_MSET_MODEL_REGISTER_ID (nc_halo_mass_summary, NC_TYPE_HALO_MASS_SUMMARY)

static void
nc_halo_mass_summary_class_init (NcHaloMassSummaryClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);


  model_class->set_property = &_nc_halo_mass_summary_set_property;
  model_class->get_property = &_nc_halo_mass_summary_get_property;
  object_class->dispose     = &_nc_halo_mass_summary_dispose;
  object_class->finalize    = &_nc_halo_mass_summary_finalize;

  ncm_model_class_set_name_nick (model_class, "Halo mass summary", "NcHaloMassSummary");
  ncm_model_class_add_params (model_class, 0, 0, 1);

  ncm_mset_model_register_id (model_class, "NcHaloMassSummary", "Halo mass summary.", NULL, TRUE, nc_halo_density_profile_id ());
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  ncm_model_class_check_params_info (NCM_MODEL_CLASS (klass));

  /**
   * NcHaloMassSummary:mass-def:
   *
   * Background density $\rho_\mathrm{bg}$ used in the mass definition \eqref{eq:mrr}.
   * See the enumerator #NcHaloMassSummaryMassDef for more details about the
   * background density definition.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_MDEF,
                                   g_param_spec_enum ("mass-def",
                                                      NULL,
                                                      "Mass definition",
                                                      NC_TYPE_HALO_MASS_SUMMARY_MASS_DEF, NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloMassSummary:Delta:
   *
   * Constant that indicates the overdensity with respect to the background density $\rho_\mathrm{bg}$.
   * See #NcHaloMassSummary:mass-def.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA,
                                   g_param_spec_double ("Delta",
                                                        NULL,
                                                        "Overdensity constant",
                                                        100.0, 3200.0, 200.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->mass          = &_nc_halo_mass_summary_mass;
  klass->concentration = &_nc_halo_mass_summary_concentration;
}

/**
 * nc_halo_mass_summary_ref:
 * @hms: a #NcHaloMassSummary
 *
 * Increase the reference of @hms by one.
 *
 * Returns: (transfer full): @hms.
 */
NcHaloMassSummary *
nc_halo_mass_summary_ref (NcHaloMassSummary *hms)
{
  return g_object_ref (hms);
}

/**
 * nc_halo_mass_summary_free:
 * @hms: a #NcHaloMassSummary
 *
 * Decrease the reference count of @hms by one.
 *
 */
void
nc_halo_mass_summary_free (NcHaloMassSummary *hms)
{
  g_object_unref (hms);
}

/**
 * nc_halo_mass_summary_clear:
 * @hms: a #NcHaloMassSummary
 *
 * Decrease the reference count of @hms by one, and sets the pointer *@hms to
 * NULL.
 *
 */
void
nc_halo_mass_summary_clear (NcHaloMassSummary **hms)
{
  g_clear_object (hms);
}

/**
 * nc_halo_mass_summary_mass: (virtual mass)
 * @hms: a #NcHaloMassSummary
 *
 * Computes the halo mass.
 * The specific implementation is provided by the child classes.
 * In general, mass will be a parameter of the model.
 *
 * Returns: the halo mass.
 */

gdouble
nc_halo_mass_summary_mass (NcHaloMassSummary *hms)
{
  return NC_HALO_MASS_SUMMARY_GET_CLASS (hms)->mass (hms);
}

/**
 * nc_halo_mass_summary_concentration: (virtual concentration)
 * @hms: a #NcHaloMassSummary
 *
 * Computes the concentration.
 * The specific implementation is provided by the child classes.
 * Concentration can be a parameter or defined by a mass-concentration relation.
 *
 * Returns: the concentration.
 */
gdouble
nc_halo_mass_summary_concentration (NcHaloMassSummary *hms)
{
  return NC_HALO_MASS_SUMMARY_GET_CLASS (hms)->concentration (hms);
}

#define _VIRIAL_DELTA(x) (18.0 * M_PI * M_PI + 82.0 * (x) - 39.0 * (x) * (x))

/**
 * nc_halo_mass_summary_Delta:
 * @hms: a #NcHaloMassSummary
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the overdensity with respect to the mass density $\Delta$.
 *
 * The virial overdensity in units of the critical density.
 * Following Colossus code (Diemer 2018) INCLUIR REF!
 * This function uses the fitting formula of Bryan & Norman 1998 INCLUIR REF!
 *
 * Returns: the value of $\Delta$.
 */
gdouble
nc_halo_mass_summary_Delta (NcHaloMassSummary *hms, NcHICosmo *cosmo, const gdouble z)
{
  NcHaloMassSummaryPrivate * const self = nc_halo_mass_summary_get_instance_private (hms);

  switch (self->mdef)
  {
    case NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN:
    case NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL:

      return self->Delta;

    case NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL:
    {
      const gdouble x = nc_hicosmo_E2Omega_m (cosmo, z) / nc_hicosmo_E2 (cosmo, z) - 1.0;

      return _VIRIAL_DELTA (x);
    }
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return 0.0; /* LCOV_EXCL_LINE */
  }
}

/**
 * nc_halo_mass_summary_rho_bg:
 * @hms: a #NcHaloMassSummary
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the background mass density $\rho_\mathrm{bg}$ in $M_\odot\mathrm{Mpc}^{-3}$.
 *
 * Returns: the value of $\rho_\mathrm{bg}\;\left[M_\odot\mathrm{Mpc}^{-3}\right]$.
 */
gdouble
nc_halo_mass_summary_rho_bg (NcHaloMassSummary *hms, NcHICosmo *cosmo, const gdouble z)
{
  NcHaloMassSummaryPrivate * const self = nc_halo_mass_summary_get_instance_private (hms);

  switch (self->mdef)
  {
    case NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN:

      return ncm_c_crit_mass_density_h2_solar_mass_Mpc3 () * nc_hicosmo_h2 (cosmo) * nc_hicosmo_E2Omega_m (cosmo, z);

    case NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL:
    case NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL:

      return ncm_c_crit_mass_density_h2_solar_mass_Mpc3 () * nc_hicosmo_h2 (cosmo) * nc_hicosmo_E2 (cosmo, z);

    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return 0.0; /* LCOV_EXCL_LINE */
  }
}

/**
 * nc_halo_mass_summary_Delta_rho_bg:
 * @hms: a #NcHaloMassSummary
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the mass density threshold $\Delta\,\rho_bg$ in $M_\odot\mathrm{Mpc}^{-3}$.
 *
 * Returns: the value of $\Delta\,\rho_bg\;\left[M_\odot\mathrm{Mpc}^{-3}\right]$.
 */
gdouble
nc_halo_mass_summary_Delta_rho_bg (NcHaloMassSummary *hms, NcHICosmo *cosmo, const gdouble z)
{
  NcHaloMassSummaryPrivate * const self = nc_halo_mass_summary_get_instance_private (hms);

  switch (self->mdef)
  {
    case NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN:

      return self->Delta * ncm_c_crit_mass_density_h2_solar_mass_Mpc3 () * nc_hicosmo_h2 (cosmo) * nc_hicosmo_E2Omega_m (cosmo, z);

    case NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL:

      return self->Delta * ncm_c_crit_mass_density_h2_solar_mass_Mpc3 () * nc_hicosmo_h2 (cosmo) * nc_hicosmo_E2 (cosmo, z);

    case NC_HALO_MASS_SUMMARY_MASS_DEF_VIRIAL:
    {
      const gdouble E2 = nc_hicosmo_E2 (cosmo, z);
      const gdouble x  = nc_hicosmo_E2Omega_m (cosmo, z) / E2 - 1.0;

      return _VIRIAL_DELTA (x) * ncm_c_crit_mass_density_h2_solar_mass_Mpc3 () * nc_hicosmo_h2 (cosmo) * E2;
    }
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */

      return 0.0; /* LCOV_EXCL_LINE */
  }
}

