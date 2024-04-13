/***************************************************************************
 *            nc_halo_mass_function.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_halo_mass_function
 * @title: NcHaloMassFunction
 * @short_description: Clusters mass function.
 *
 * Class that implements the mass function of clusters
 * dark matter halos.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_halo_mass_function.h"
#include "math/ncm_integrate.h"
#include "math/ncm_spline2d_bicubic.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_DISTANCE,
  PROP_PSF,
  PROP_MULTIPLICITY,
  PROP_AREA,
  PROP_PREC,
  PROP_LNMI,
  PROP_LNMF,
  PROP_ZI,
  PROP_ZF,
  PROP_MF_LB,
  PROP_SIZE,
};

struct _NcHaloMassFunctionPrivate
{
  NcDistance *dist;
  NcMultiplicityFunc *mulf;
  NcmPowspecFilter *psf;
  gdouble area_survey;
  gdouble lnMi;
  gdouble lnMf;
  gdouble zi;
  gdouble zf;
  gdouble prec;
  NcmModelCtrl *ctrl_cosmo;
  NcmModelCtrl *ctrl_reion;
  gboolean constructed;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHaloMassFunction, nc_halo_mass_function, G_TYPE_OBJECT)

static void
nc_halo_mass_function_init (NcHaloMassFunction *mfp)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv = nc_halo_mass_function_get_instance_private (mfp);

  mfp->d2NdzdlnM    = NULL;
  mfp->mf_lb        = 0.0;
  self->lnMi        = 0.0;
  self->lnMf        = 0.0;
  self->zi          = 0.0;
  self->zf          = 0.0;
  self->area_survey = 0.0;
  self->dist        = NULL;
  self->mulf        = NULL;
  self->psf         = NULL;
  self->prec        = 0.0;
  self->ctrl_cosmo  = ncm_model_ctrl_new (NULL);
  self->ctrl_reion  = ncm_model_ctrl_new (NULL);
  self->constructed = FALSE;
}

static void
_nc_halo_mass_function_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_halo_mass_function_parent_class)->constructed (object);

  {
    NcHaloMassFunction *mfp                = NC_HALO_MASS_FUNCTION (object);
    NcHaloMassFunctionPrivate * const self = mfp->priv;

    g_assert_cmpfloat (self->lnMi, <, self->lnMf);
    g_assert_cmpfloat (self->zi, <, self->zf);

    ncm_powspec_filter_require_zi (self->psf, self->zi);
    ncm_powspec_filter_require_zf (self->psf, self->zf);

    self->constructed = TRUE;
  }
}

static void
_nc_halo_mass_function_dispose (GObject *object)
{
  NcHaloMassFunction *mfp                = NC_HALO_MASS_FUNCTION (object);
  NcHaloMassFunctionPrivate * const self = mfp->priv;

  nc_distance_clear (&self->dist);
  nc_multiplicity_func_clear (&self->mulf);

  ncm_powspec_filter_clear (&self->psf);

  ncm_model_ctrl_clear (&self->ctrl_cosmo);
  ncm_model_ctrl_clear (&self->ctrl_reion);

  ncm_spline2d_clear (&mfp->d2NdzdlnM);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_mass_function_parent_class)->dispose (object);
}

static void
_nc_halo_mass_function_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_mass_function_parent_class)->finalize (object);
}

static void
_nc_halo_mass_function_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloMassFunction *mfp                = NC_HALO_MASS_FUNCTION (object);
  NcHaloMassFunctionPrivate * const self = mfp->priv;

  g_return_if_fail (NC_IS_HALO_MASS_FUNCTION (object));

  switch (prop_id)
  {
    case PROP_DISTANCE:
      self->dist = g_value_dup_object (value);
      break;
    case PROP_PSF:
      self->psf = g_value_dup_object (value);
      break;
    case PROP_MULTIPLICITY:
      self->mulf = g_value_dup_object (value);
      break;
    case PROP_AREA:
      nc_halo_mass_function_set_area (mfp, g_value_get_double (value));
      break;
    case PROP_PREC:
      nc_halo_mass_function_set_prec (mfp, g_value_get_double (value));
      break;
    case PROP_LNMI:
      self->lnMi = g_value_get_double (value);
      break;
    case PROP_LNMF:
      self->lnMf = g_value_get_double (value);
      break;
    case PROP_ZI:
      self->zi = g_value_get_double (value);

      if (self->constructed)
        ncm_powspec_filter_require_zi (self->psf, self->zi);

      break;
    case PROP_ZF:
      self->zf = g_value_get_double (value);

      if (self->constructed)
        ncm_powspec_filter_require_zf (self->psf, self->zf);

      break;
    case PROP_MF_LB:
      mfp->mf_lb = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_mass_function_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloMassFunction *mfp                = NC_HALO_MASS_FUNCTION (object);
  NcHaloMassFunctionPrivate * const self = mfp->priv;

  g_return_if_fail (NC_IS_HALO_MASS_FUNCTION (object));

  switch (prop_id)
  {
    case PROP_DISTANCE:
      g_value_set_object (value, self->dist);
      break;
    case PROP_PSF:
      g_value_set_object (value, self->psf);
      break;
    case PROP_MULTIPLICITY:
      g_value_set_object (value, self->mulf);
      break;
    case PROP_AREA:
      g_value_set_double (value, self->area_survey);
      break;
    case PROP_PREC:
      g_value_set_double (value, self->prec);
      break;
    case PROP_LNMI:
      g_value_set_double (value, self->lnMi);
      break;
    case PROP_LNMF:
      g_value_set_double (value, self->lnMf);
      break;
    case PROP_ZI:
      g_value_set_double (value, self->zi);
      break;
    case PROP_ZF:
      g_value_set_double (value, self->zf);
      break;
    case PROP_MF_LB:
      g_value_set_double (value, mfp->mf_lb);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_halo_mass_function_class_init (NcHaloMassFunctionClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  /*GObjectClass* parent_class = G_OBJECT_CLASS (klass); */

  object_class->constructed  = &_nc_halo_mass_function_constructed;
  object_class->dispose      = &_nc_halo_mass_function_dispose;
  object_class->finalize     = &_nc_halo_mass_function_finalize;
  object_class->set_property = &_nc_halo_mass_function_set_property;
  object_class->get_property = &_nc_halo_mass_function_get_property;

  /**
   * NcHaloMassFunction:distance:
   *
   * This property keeps the distance object.
   */
  g_object_class_install_property (object_class,
                                   PROP_DISTANCE,
                                   g_param_spec_object ("distance",
                                                        NULL,
                                                        "Distance",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME
                                                        | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloMassFunction:powerspectrum-filtered:
   *
   * This property keeps the filtered power spectrum.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_PSF,
                                   g_param_spec_object ("powerspectrum-filtered",
                                                        NULL,
                                                        "Filtered power-spectrum",
                                                        NCM_TYPE_POWSPEC_FILTER,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloMassFunction:multiplicity:
   *
   * This property keeps the multiplicity function object.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_MULTIPLICITY,
                                   g_param_spec_object ("multiplicity",
                                                        NULL,
                                                        "Multiplicity function",
                                                        NC_TYPE_MULTIPLICITY_FUNC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloMassFunction:area:
   *
   * This property sets the angular area in steradian.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_AREA,
                                   g_param_spec_double ("area",
                                                        NULL,
                                                        "Angular area in steradian",
                                                        0.0, 4.0 * M_PI, 200.0 * M_PI * M_PI / (180.0 * 180.0),
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloMassFunction:prec:
   *
   * This property sets the precision.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_PREC,
                                   g_param_spec_double ("prec",
                                                        NULL,
                                                        "Precision",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloMassFunction:lnMi:
   *
   * This property sets the minimum halo mass (logarithm base e).
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LNMI,
                                   g_param_spec_double ("lnMi",
                                                        NULL,
                                                        "Lower mass",
                                                        log (1.0e11), log (1.0e17), log (1.0e13),
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloMassFunction:lnMf:
   *
   * This property sets the maximum halo mass (logarithm base e) $\ln M_f$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LNMF,
                                   g_param_spec_double ("lnMf",
                                                        NULL,
                                                        "Upper mass",
                                                        log (1.0e11), log (1.0e17), log (1.0e16),
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloMassFunction:zi:
   *
   * This property sets the initial redshift $z_i$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ZI,
                                   g_param_spec_double ("zi",
                                                        NULL,
                                                        "Lower redshift",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloMassFunction:zf:
   *
   * This property sets the final redshift $z_f$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ZF,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "Upper redshift",
                                                        0.0, G_MAXDOUBLE, 1.4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcHaloMassFunction:mf-lb:
   *
   * Mass function lower bound.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_MF_LB,
                                   g_param_spec_double ("mf-lb",
                                                        NULL,
                                                        "Upper redshift",
                                                        0.0, G_MAXDOUBLE, 1.0e-30,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_halo_mass_function_new:
 * @dist: a #NcDistance sets to #NcHaloMassFunction:distance
 * @psf: a #NcmPowspecFilter sets to #NcHaloMassFunction:powerspectrum-filtered
 * @mulf: a #NcMultiplicityFunc sets to #NcHaloMassFunction:multiplicity
 *
 * This function allocates memory for a new #NcHaloMassFunction object and sets its properties to the values from
 * the input arguments.
 *
 * Returns: A new #NcHaloMassFunction.
 */
NcHaloMassFunction *
nc_halo_mass_function_new (NcDistance *dist, NcmPowspecFilter *psf, NcMultiplicityFunc *mulf)
{
  NcHaloMassFunction *mfp = g_object_new (NC_TYPE_HALO_MASS_FUNCTION,
                                          "distance", dist,
                                          "powerspectrum-filtered", psf,
                                          "multiplicity", mulf,
                                          NULL);

  return mfp;
}

/**
 * nc_halo_mass_function_free:
 * @mfp: a #NcHaloMassFunction
 *
 * Atomically decrements the reference count of @mfp by one. If the reference count drops to 0,
 * all memory allocated by @mfp is released.
 *
 */
void
nc_halo_mass_function_free (NcHaloMassFunction *mfp)
{
  g_object_unref (mfp);
}

/**
 * nc_halo_mass_function_clear:
 * @mfp: a #NcHaloMassFunction
 *
 * Atomically decrements the reference count of @mfp by one. If the reference count drops to 0,
 * all memory allocated by @mfp is released. Set pointer to NULL.
 *
 */
void
nc_halo_mass_function_clear (NcHaloMassFunction **mfp)
{
  g_clear_object (mfp);
}

/**
 * nc_halo_mass_function_lnM_to_lnR:
 * @mfp: a #NcHaloMassFunction
 * @cosmo: a #NcHICosmo
 * @lnM: logarithm base $e$ of the mass enclosed in the filter volume $\ln M$
 *
 * This function computes the ln-radius of related to the ln-mass $\ln(M / M_\odot)$.
 *
 * Returns: $\ln(R / 1 \text{Mpc})$
 */
gdouble
nc_halo_mass_function_lnM_to_lnR (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnM)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;
  const gdouble omega_m0                 = nc_hicosmo_Omega_m0h2 (cosmo);

  return (lnM - log (omega_m0 * ncm_powspec_filter_volume_rm3 (self->psf) * ncm_c_crit_mass_density_h2_solar_mass_Mpc3 ())) / 3.0;
}

/**
 * nc_halo_mass_function_lnR_to_lnM:
 * @mfp: a #NcHaloMassFunction
 * @cosmo: a #NcHICosmo
 * @lnR: ln-radius of the related volume $\ln R$
 *
 * This function computes the ln-mass of the mass enclosed in the filter volume.
 *
 * Returns: $\ln(M / 1 M_\odot)$
 */
gdouble
nc_halo_mass_function_lnR_to_lnM (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnR)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;
  const gdouble omega_m0                 = nc_hicosmo_Omega_m0h2 (cosmo);

  return (3.0 * lnR + log (omega_m0 * ncm_powspec_filter_volume_rm3 (self->psf) * ncm_c_crit_mass_density_h2_solar_mass_Mpc3 ()));
}

/**
 * nc_halo_mass_function_sigma_lnR:
 * @mfp: a #NcHaloMassFunction
 * @cosmo: a #NcHICosmo
 * @lnR: logarithm base e of the radius $\ln R$
 * @z: redshift $z$
 *
 * This function computes the matter variance for $\ln R$ = @lnR
 * at redshift @z.
 *
 */
gdouble
nc_halo_mass_function_sigma_lnR (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnR, gdouble z)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;
  const gdouble sigma                    = ncm_powspec_filter_eval_sigma_lnr (self->psf, z, lnR);

  return sigma;
}

/**
 * nc_halo_mass_function_sigma_lnM:
 * @mfp: a #NcHaloMassFunction
 * @cosmo: a #NcHICosmo
 * @lnM: logarithm base e of mass enclosed in the filter volume $\ln M$
 * @z: redshift $z$
 *
 * This function computes the matter variance for $\ln R$ = @lnR
 * at redshift @z.
 *
 */
gdouble
nc_halo_mass_function_sigma_lnM (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;
  const gdouble lnR                      = nc_halo_mass_function_lnM_to_lnR (mfp, cosmo, lnM);
  const gdouble sigma                    = ncm_powspec_filter_eval_sigma_lnr (self->psf, z, lnR);

  return sigma;
}

/**
 * nc_halo_mass_function_dn_dlnR:
 * @mfp: a #NcHaloMassFunction
 * @cosmo: a #NcHICosmo
 * @lnR: logarithm base e of the radius $\ln R$
 * @z: redshift $z$
 *
 * This function computes the comoving number density of dark matter halos per redshift @z and
 * volume with ln-radius @lnR.
 *
 * Returns: $\frac{\mathrm{d}n}{\mathrm{d}\ln(R)}$.
 */
gdouble
nc_halo_mass_function_dn_dlnR (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnR, gdouble z)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;
  const gdouble V                        = ncm_powspec_filter_volume_rm3 (self->psf) * exp (3.0 * lnR);
  const gdouble sigma                    = ncm_powspec_filter_eval_sigma_lnr (self->psf, z, lnR);
  const gdouble dlnvar_dlnR              = ncm_powspec_filter_eval_dlnvar_dlnr (self->psf, z, lnR);
  const gdouble f                        = nc_multiplicity_func_eval (self->mulf, cosmo, sigma, z);
  const gdouble dn_dlnR                  = -(1.0 / V) * f * 0.5 * dlnvar_dlnR;

  if (nc_multiplicity_func_has_correction_factor (self->mulf))
  {
    const gdouble lnM = nc_halo_mass_function_lnR_to_lnM (mfp, cosmo, lnR);

    return dn_dlnR * nc_multiplicity_func_correction_factor (self->mulf, cosmo, sigma, z, lnM);
  }
  else
  {
    return dn_dlnR;
  }
}

/**
 * nc_halo_mass_function_dn_dlnM:
 * @mfp: a #NcHaloMassFunction
 * @cosmo: a #NcHICosmo
 * @lnM: logarithm base e of mass enclosed in the filter volume $\ln M$
 * @z: redshift $z$
 *
 * This function computes the comoving number density of dark matter halos at redshift @z and
 * mass $M$.
 *
 * Returns: $\frac{\mathrm{d}n}{\mathrm{d}\ln(M)}$.
 */
gdouble
nc_halo_mass_function_dn_dlnM (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  const gdouble lnR     = nc_halo_mass_function_lnM_to_lnR (mfp, cosmo, lnM);
  const gdouble dn_dlnR = nc_halo_mass_function_dn_dlnR (mfp, cosmo, lnR, z);
  const gdouble dn_dlnM = dn_dlnR / 3.0;

  return dn_dlnM;
}

typedef struct _nc_ca_integ
{
  NcHaloMassFunction *mfp;
  NcHICosmo *cosmo;
  gdouble z;
} nc_ca_integ;

/**
 * nc_halo_mass_function_dv_dzdomega:
 * @mfp: a #NcHaloMassFunction
 * @cosmo: a #NcHICosmo
 * @z: redshift $z$
 *
 * This function computes the comoving volume (flat universe) element per unit solid angle $d\Omega$
 * given @z, namely, $$\frac{\mathrm{d}^2V}{\mathrm{d}z\mathrm{d}\Omega} = \frac{c}{H(z)} D_c^2(z),$$
 * where $H(z)$ is the Hubble function and $D_c$ is the comoving distance.
 *
 * Returns: comoving volume element $\frac{\mathrm{d}^2V}{\mathrm{d}z\mathrm{d}\Omega} \,\left[\mathrm{Mpc}^3\right]$.
 */
gdouble
nc_halo_mass_function_dv_dzdomega (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble z)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;
  const gdouble RH                       = nc_hicosmo_RH_Mpc (cosmo);
  const gdouble VH                       = gsl_pow_3 (RH);
  const gdouble E                        = sqrt (nc_hicosmo_E2 (cosmo, z));
  gdouble dc                             = nc_distance_comoving (self->dist, cosmo, z);
  gdouble dV_dzdOmega                    = VH * gsl_pow_2 (dc) / E;

  return dV_dzdOmega;
}

static void _nc_halo_mass_function_generate_2Dspline_knots (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble rel_error);

/**
 * nc_halo_mass_function_set_area:
 * @mfp: a #NcHaloMassFunction
 * @area: area in steradian
 *
 * Sets the area of the survey in steradian.
 *
 */
void
nc_halo_mass_function_set_area (NcHaloMassFunction *mfp, gdouble area)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;

  if (self->area_survey != area)
  {
    self->area_survey = area;
    ncm_spline2d_clear (&mfp->d2NdzdlnM);

    ncm_model_ctrl_force_update (self->ctrl_cosmo);
    ncm_model_ctrl_force_update (self->ctrl_reion);
  }
}

/**
 * nc_halo_mass_function_set_prec:
 * @mfp: a #NcHaloMassFunction
 * @prec: precision
 *
 * Sets the precision of the integration.
 *
 */
void
nc_halo_mass_function_set_prec (NcHaloMassFunction *mfp, gdouble prec)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;

  if (self->prec != prec)
  {
    self->prec = prec;
    ncm_spline2d_clear (&mfp->d2NdzdlnM);

    ncm_model_ctrl_force_update (self->ctrl_cosmo);
    ncm_model_ctrl_force_update (self->ctrl_reion);
  }
}

/**
 * nc_halo_mass_function_set_area_sd:
 * @mfp: a #NcHaloMassFunction
 * @area_sd: area in square degree
 *
 * Sets the area of the survey in square degree.
 *
 */
void
nc_halo_mass_function_set_area_sd (NcHaloMassFunction *mfp, gdouble area_sd)
{
  const gdouble conversion_factor = M_PI * M_PI / (180.0 * 180.0);

  nc_halo_mass_function_set_area (mfp, area_sd * conversion_factor);
}

/**
 * nc_halo_mass_function_set_eval_limits:
 * @mfp: a #NcHaloMassFunction
 * @cosmo: a #NcHICosmo
 * @lnMi: minimum logarithm base e of mass
 * @lnMf: maximum logarithm base e of mass
 * @zi: minimum redshift
 * @zf: maximum redshift
 *
 * Sets the limits of the integration.
 *
 */
void
nc_halo_mass_function_set_eval_limits (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnMi, gdouble lnMf, gdouble zi, gdouble zf)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;

  g_assert_cmpfloat (lnMi, <, lnMf);
  g_assert_cmpfloat (zi, <, zf);

  if ((lnMi != self->lnMi) || (self->lnMf != lnMf) || (self->zi != zi) || (self->zf != zf))
  {
    g_object_set (G_OBJECT (mfp),
                  "lnMi", lnMi,
                  "lnMf", lnMf,
                  "zi", zi,
                  "zf", zf,
                  NULL);
    ncm_powspec_filter_require_zi (self->psf, self->zi);
    ncm_powspec_filter_require_zf (self->psf, self->zf);

    ncm_spline2d_clear (&mfp->d2NdzdlnM);

    ncm_model_ctrl_force_update (self->ctrl_cosmo);
    ncm_model_ctrl_force_update (self->ctrl_reion);
  }
}

/**
 * nc_halo_mass_function_prepare:
 * @mfp: a #NcHaloMassFunction
 * @cosmo: a #NcHICosmo
 *
 * Prepares the halo mass function object.
 *
 */
void
nc_halo_mass_function_prepare (NcHaloMassFunction *mfp, NcHICosmo *cosmo)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;
  guint i, j;

  nc_distance_prepare_if_needed (self->dist, cosmo);
  ncm_powspec_filter_prepare_if_needed (self->psf, NCM_MODEL (cosmo));

  if (mfp->d2NdzdlnM == NULL)
    _nc_halo_mass_function_generate_2Dspline_knots (mfp, cosmo, self->prec);

  {
    NcmVector *z_vec   = ncm_spline2d_peek_yv (mfp->d2NdzdlnM);
    NcmVector *lnM_vec = ncm_spline2d_peek_xv (mfp->d2NdzdlnM);
    NcmMatrix *val_mat = ncm_spline2d_peek_zm (mfp->d2NdzdlnM);

    for (i = 0; i < ncm_vector_len (z_vec); i++)
    {
      const gdouble z    = ncm_vector_get (z_vec, i);
      const gdouble dVdz = self->area_survey * nc_halo_mass_function_dv_dzdomega (mfp, cosmo, z);

      for (j = 0; j < ncm_vector_len (lnM_vec); j++)
      {
        const gdouble lnM          = ncm_vector_get (lnM_vec, j);
        const gdouble d2NdzdlnM_ij = (dVdz * nc_halo_mass_function_dn_dlnM (mfp, cosmo, lnM, z));

        ncm_matrix_set (val_mat, i, j, GSL_MAX (d2NdzdlnM_ij, mfp->mf_lb));
      }
    }
  }
  ncm_spline2d_prepare (mfp->d2NdzdlnM);

  ncm_model_ctrl_update (self->ctrl_cosmo, NCM_MODEL (cosmo));
}

/**
 * nc_halo_mass_function_dn_dz:
 * @mfp: a #NcHaloMassFunction
 * @cosmo: a #NcHICosmo
 * @lnMl: logarithm base e of mass, lower threshold $\ln M_l$
 * @lnMu: logarithm base e of mass, upper threshold $\ln M_u$
 * @z: redshift $z$
 * @spline: whenever to create an intermediary spline of the mass integration
 *
 * Computes the number density of halos in the mass range $[M_l, M_u]$ at redshift $z$.
 *
 * Returns: the number density of halos in the mass range $[M_l, M_u]$
 */
gdouble
nc_halo_mass_function_dn_dz (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnMl, gdouble lnMu, gdouble z, gboolean spline)
{
  gdouble dN_dz;

  if (spline)
    dN_dz = ncm_spline2d_integ_dx_spline_val (mfp->d2NdzdlnM, lnMl, lnMu, z);
  else
    dN_dz = ncm_spline2d_integ_dx (mfp->d2NdzdlnM, lnMl, lnMu, z);

  return dN_dz;
}

/**
 * nc_halo_mass_function_n:
 * @mfp: a #NcHaloMassFunction
 * @cosmo: a #NcHICosmo
 * @lnMl: logarithm base e of mass, lower threshold $\ln M_l$
 * @lnMu: logarithm base e of mass, upper threshold $\ln M_u$
 * @zl: minimum redshift
 * @zu: maximum redshift
 * @spline: whenever to create an intermediary spline of the integration
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_halo_mass_function_n (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble lnMl, gdouble lnMu, gdouble zl, gdouble zu, NcHaloMassFunctionSplineOptimize spline)
{
  gdouble N;

  switch (spline)
  {
    case NC_HALO_MASS_FUNCTION_SPLINE_NONE:
      N = ncm_spline2d_integ_dxdy (mfp->d2NdzdlnM, lnMl, lnMu, zl, zu);
      break;
    case NC_HALO_MASS_FUNCTION_SPLINE_LNM:
      N = ncm_spline2d_integ_dxdy_spline_x (mfp->d2NdzdlnM, lnMl, lnMu, zl, zu);
      break;
    case NC_HALO_MASS_FUNCTION_SPLINE_Z:
      N = ncm_spline2d_integ_dxdy_spline_y (mfp->d2NdzdlnM, lnMl, lnMu, zl, zu);
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  return N;
}

typedef struct __encapsulated_function_args
{
  NcHaloMassFunction *mfp;
  NcHaloMassFunctionPrivate *self;
  NcHICosmo *cosmo;
  gdouble z;
  gdouble lnM;
  gdouble dVdz;
} _encapsulated_function_args;

static gdouble
_encapsulated_z (gdouble z, gpointer p)
{
  _encapsulated_function_args *args = (_encapsulated_function_args *) p;

  gdouble A = args->self->area_survey *
              nc_halo_mass_function_dv_dzdomega (args->mfp, args->cosmo, z) *
              nc_halo_mass_function_dn_dlnM (args->mfp, args->cosmo, args->lnM, z);

  /*printf ("z   % 22.15g % 22.15g\n", z, A);fflush(stdout);*/
  return A;
}

static gdouble
_encapsulated_lnM (gdouble lnM, gpointer p)
{
  _encapsulated_function_args *args = (_encapsulated_function_args *) p;

  gdouble A = args->dVdz * nc_halo_mass_function_dn_dlnM (args->mfp, args->cosmo, lnM, args->z);

  /*printf ("lnM % 22.15g % 22.15g\n", lnM, A);fflush (stdout);*/
  return A;
}

static void
_nc_halo_mass_function_generate_2Dspline_knots (NcHaloMassFunction *mfp, NcHICosmo *cosmo, gdouble rel_error)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;
  gsl_function Fx, Fy;
  _encapsulated_function_args args;

  g_assert (mfp->d2NdzdlnM == NULL);
  g_assert_cmpfloat (self->lnMi, <, self->lnMf);
  g_assert_cmpfloat (self->zi, <, self->zf);

  args.mfp   = mfp;
  args.self  = self;
  args.cosmo = cosmo;
  args.z     = (self->zf + self->zi) / 2.0;
  args.lnM   = (self->lnMf + self->lnMi) / 2.0;
  args.dVdz  = self->area_survey * nc_halo_mass_function_dv_dzdomega (mfp, cosmo, args.z);

  Fx.function = _encapsulated_lnM;
  Fx.params   = &args;

  Fy.function = _encapsulated_z;
  Fy.params   = &args;
  /*ncm_model_orig_params_log_all (NCM_MODEL (cosmo));*/
  mfp->d2NdzdlnM = ncm_spline2d_bicubic_notaknot_new ();
  ncm_spline2d_set_function (mfp->d2NdzdlnM,
                             NCM_SPLINE_FUNCTION_SPLINE, &Fx, &Fy, self->lnMi, self->lnMf, self->zi, self->zf, rel_error);
}

/**
 * nc_halo_mass_function_prepare_if_needed:
 * @mfp: a #NcHaloMassFunction
 * @cosmo: a #NcHICosmo
 *
 * Prepare the object if necessary.
 *
 */
void
nc_halo_mass_function_prepare_if_needed (NcHaloMassFunction *mfp, NcHICosmo *cosmo)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;
  gboolean cosmo_up                      = ncm_model_ctrl_update (self->ctrl_cosmo, NCM_MODEL (cosmo));

  if (cosmo_up)
    nc_halo_mass_function_prepare (mfp, cosmo);
}

/**
 * nc_halo_mass_function_peek_psf:
 * @mfp: a #NcHaloMassFunction
 *
 * Peeks the power spectrum filter.
 *
 * Returns: (transfer none): the power spectrum filter.
 */
NcmPowspecFilter *
nc_halo_mass_function_peek_psf (NcHaloMassFunction *mfp)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;

  return self->psf;
}

/**
 * nc_halo_mass_function_peek_multiplicity_function:
 * @mfp: a #NcHaloMassFunction
 *
 * Peeks the multiplicity function.
 *
 * Returns: (transfer none): the multiplicity function.
 */
NcMultiplicityFunc *
nc_halo_mass_function_peek_multiplicity_function (NcHaloMassFunction *mfp)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;

  return self->mulf;
}

/**
 * nc_halo_mass_function_d2n_dzdlnM:
 * @mfp: a #NcHaloMassFunction
 * @cosmo: a #NcHICosmo
 * @lnM: logarithm base e of mass $\ln M$
 * @z: redshift
 *
 * Computes the halo mass function $d^2 N / d\ln M d z$.
 *
 * Returns: $d^2 N / d\ln M d z$.
 */

/**
 * nc_halo_mass_function_peek_survey_area:
 * @mfp: a #NcHaloMassFunction
 *
 * Peeks the survey area.
 *
 * Returns: the survey area.
 */
gdouble
nc_halo_mass_function_peek_survey_area (NcHaloMassFunction *mfp)
{
  NcHaloMassFunctionPrivate * const self = mfp->priv;

  return self->area_survey;
}

