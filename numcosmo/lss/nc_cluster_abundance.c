/***************************************************************************
 *            nc_cluster_abundance.c
 *
 *  Tue Apr 20 10:59:01 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com> & <sandro@isoftware.com.br>
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
 * SECTION:nc_cluster_abundance
 * @title: Cluster Abundance Distribution
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>
#include <gsl/gsl_sf_erf.h>

enum
{
  PROP_0,
  PROP_OPT,
  PROP_MASS_FUNCTION,
  PROP_MEANBIAS,
  PROP_ZI,
  PROP_ZF,
  PROP_LNMI,
  PROP_LNMF,
  PROP_PHOTOZ,
  PROP_LNMS0
};

G_DEFINE_TYPE (NcClusterAbundance, nc_cluster_abundance, G_TYPE_OBJECT);

#define INTEG_D2NDZDLNM_NNODES (200)
#define LNM_MIN (10.0 * M_LN10)

static gdouble _d2NdzdlnM_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z) {g_error ("Function d2NdzdlnM_val not implemented or cad not prepared."); return 0.0;};
static gdouble _dNdz_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble z) {g_error ("Function dNdz_val not implemented or cad not prepared."); return 0.0;};
static gdouble _dNdlnM_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble zl, gdouble zu) {g_error ("Function dNdlnM_val not implemented or cad not prepared."); return 0.0;};
static gdouble _N_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble zl, gdouble zu) {g_error ("Function N_val not implemented or cad not prepared."); return 0.0;};

/**
 * nc_cluster_abundance_new:
 * @opt: a #NcClusterAbundanceOpt.
 * @mfp: a #NcMassFunction.
 * @mbiasf: (allow-none): a #NcHaloBiasFunc.
 * @zi: minimum redshift.
 * @zf: maximum redshift.
 * @lnMi: minimum logarithm base e of mass.
 * @lnMf: maximum logarithm base e of mass.
 * @photoz: a #NcClusterPhotoz.
 * @lnM_sigma0: FIXME
 *
 * This function allocates memory for a new #NcClusterAbundance object and sets its properties to the values from
 * the input arguments.
 *
 * @Returns: A new #NcClusterAbundance.
 */
NcClusterAbundance *
nc_cluster_abundance_new (NcClusterAbundanceOpt opt, NcMassFunction *mfp, NcHaloBiasFunc *mbiasf, gdouble zi, gdouble zf, gdouble lnMi, gdouble lnMf, NcClusterPhotoz *photoz, gdouble lnM_sigma0)
{
  NcClusterAbundance *cad = g_object_new (NC_TYPE_CLUSTER_ABUNDANCE,
                                          "options",                opt,
                                          "mass-function",          mfp,
                                          "mean-bias-function",     mbiasf,
                                          "minimum-redshift",       zi,
                                          "maximum-redshift",       zf,
                                          "minimum-mass",           lnMi,
                                          "maximum-mass",           lnMf,
                                          "photoz",                 photoz,
                                          "mass-observable-sigma0", lnM_sigma0,
                                          NULL);
  return cad;
}

static void
nc_cluster_abundance_init (NcClusterAbundance *cad)
{
  /* TODO: Add initialization code here */
  NcmVector *u1_knots = ncm_vector_new (INTEG_D2NDZDLNM_NNODES);
  NcmVector *u2_knots = ncm_vector_new (INTEG_D2NDZDLNM_NNODES);
  NcmVector *u3_knots = ncm_vector_new (INTEG_D2NDZDLNM_NNODES);
  NcmVector *integ_lnM_knots = ncm_vector_new (INTEG_D2NDZDLNM_NNODES);
  NcmVector *integ_z_knots = ncm_vector_new (INTEG_D2NDZDLNM_NNODES);
  NcmMatrix *integ_lnM_z_knots = ncm_matrix_new (INTEG_D2NDZDLNM_NNODES, INTEG_D2NDZDLNM_NNODES);

  cad->d2NdzdlnM_val = &_d2NdzdlnM_val;
  cad->dNdz_val      = &_dNdz_val;
  cad->dNdlnM_val    = &_dNdlnM_val;
  cad->N_val         = &_N_val;

  cad->completeness = NULL;
  cad->purity = NULL;
  cad->sd_lnM = NULL;
  cad->norma = 0.0;
  cad->log_norma = 0.0;
  cad->completeness_factor = 1.0; //0.874737945; (work with Camila and Alex)
  cad->optimize = TRUE;

#define D2NDZDLNM_Z(cad) ((cad)->d2NdzdlnM->yv)
#define D2NDZDLNM_LNM(cad) ((cad)->d2NdzdlnM->xv)
#define D2NDZDLNM_VAL(cad) ((cad)->d2NdzdlnM->zm)

#define DBDLNM_Z(cad) ((cad)->dbdlnM->yv)
#define DBDLNM_LNM(cad) ((cad)->dbdlnM->xv)
#define DBDLNM_VAL(cad) ((cad)->dbdlnM->zm)

  cad->dNdz = NULL;
  cad->d2NdzdlnM = NULL;
  cad->dbdlnM = NULL;

  cad->inv_z     = ncm_spline_cubic_notaknot_new ();
  cad->inv_lnM   = ncm_spline_cubic_notaknot_new ();

  ncm_spline_set (cad->inv_z, u1_knots, integ_z_knots, FALSE);
  ncm_spline_set (cad->inv_lnM, u2_knots, integ_lnM_knots, FALSE);

  cad->inv_lnM_z = ncm_spline2d_bicubic_notaknot_new ();
  ncm_spline2d_set (cad->inv_lnM_z, u3_knots, integ_z_knots, integ_lnM_z_knots, FALSE);

  cad->ctrl = ncm_model_ctrl_new (NULL);
}

/**
 * nc_cluster_abundance_copy:
 * @cad: a #NcClusterAbundance.
 *
 * Duplicates the #NcClusterAbundance object setting the same values of the original propertities.
 *
 * Returns: (transfer full): A new #NcClusterAbundance.
*/
NcClusterAbundance *
nc_cluster_abundance_copy (NcClusterAbundance *cad)
{
  return nc_cluster_abundance_new (cad->opt, cad->mfp, cad->mbiasf, cad->zi, cad->zf, cad->lnMi, cad->lnMf, cad->photoz, cad->lnM_sigma0);
}

/**
 * nc_cluster_abundance_free:
 * @cad: a #NcClusterAbundance.
 *
 * Atomically decrements the reference count of @cad by one. If the reference count drops to 0,
 * all memory allocated by @cad is released.
 *
 */
void
nc_cluster_abundance_free (NcClusterAbundance *cad)
{
  g_clear_object (&cad);
}

/**
 * nc_cluster_abundance_set_zi:
 * @cad: a #NcClusterAbundance.
 * @zi: value of #NcClusterAbundance:minimum-redshift.
 *
 * Sets the value @zi to the #NcClusterAbundance:minimum-redshift property.
 *
 */
void
nc_cluster_abundance_set_zi (NcClusterAbundance *cad, const gdouble zi)
{
  cad->zi = zi;
}

/**
 * nc_cluster_abundance_get_zi:
 * @cad: a #NcClusterAbundance.
 *
 * Gets the value of the #NcClusterAbundance:minimum-redshift property.
 *
 * Returns: the value of #NcClusterAbundance:minimum-redshift property.
 */
gdouble
nc_cluster_abundance_get_zi (NcClusterAbundance *cad)
{
  return cad->zi;
}

/**
 * nc_cluster_abundance_set_zf:
 * @cad: a #NcClusterAbundance.
 * @zf: value of #NcClusterAbundance:maximum-redshift.
 *
 * Sets the value @zf to the #NcClusterAbundance:maximum-redshift property.
 *
 */
void
nc_cluster_abundance_set_zf (NcClusterAbundance *cad, const gdouble zf)
{
  g_assert (zf > cad->zi);
  cad->zf = zf;
}

/**
 * nc_cluster_abundance_get_zf:
 * @cad: a #NcClusterAbundance.
 *
 * Gets the value of the #NcClusterAbundance:maximum-redshift property.
 *
 * Returns: the value of #NcClusterAbundance:maximum-redshift property.
 */
gdouble
nc_cluster_abundance_get_zf (NcClusterAbundance *cad)
{
  return cad->zf;
}

/**
 * nc_cluster_abundance_set_lnMi:
 * @cad: a #NcClusterAbundance.
 * @lnMi: value of #NcClusterAbundance:minimum-mass.
 *
 * Sets the value @lnMi to the #NcClusterAbundance:minimum-mass property.
 *
 */
void
nc_cluster_abundance_set_lnMi (NcClusterAbundance *cad, const gdouble lnMi)
{
  cad->lnMi = lnMi;
}

/**
 * nc_cluster_abundance_get_lnMi:
 * @cad: a #NcClusterAbundance.
 *
 * Gets the value of the #NcClusterAbundance:minimum-mass property.
 *
 * Returns: the value of #NcClusterAbundance:minimum-mass property.
 */
gdouble
nc_cluster_abundance_get_lnMi (NcClusterAbundance *cad)
{
  return cad->lnMi;
}

/**
 * nc_cluster_abundance_set_lnMf:
 * @cad: a #NcClusterAbundance.
 * @lnMf: value of #NcClusterAbundance:maximum-mass.
 *
 * Sets the value @lnMf to the #NcClusterAbundance:maximum-mass property.
 *
 */
void
nc_cluster_abundance_set_lnMf (NcClusterAbundance *cad, const gdouble lnMf)
{
  g_assert (lnMf > cad->lnMi);
  cad->lnMf = lnMf;
}

/**
 * nc_cluster_abundance_get_lnMf:
 * @cad: a #NcClusterAbundance.
 *
 * Gets the value of the #NcClusterAbundance:maximum-mass property.
 *
 * Returns: the value of #NcClusterAbundance:maximum-mass property.
 */
gdouble
nc_cluster_abundance_get_lnMf (NcClusterAbundance *cad)
{
  return cad->lnMf;
}

/**
 * nc_cluster_abundance_set_photoz:
 * @cad: a #NcClusterAbundance.
 * @photoz: value of #NcClusterAbundance:photoz.
 *
 * Sets the value @photoz to the #NcClusterAbundance:photoz property.
 *
 */
void
nc_cluster_abundance_set_photoz (NcClusterAbundance *cad, NcClusterPhotoz *photoz)
{
  cad->photoz = photoz;
}

/**
 * nc_cluster_abundance_get_photoz:
 * @cad: a #NcClusterAbundance.
 *
 * Gets the value of the #NcClusterAbundance:photoz property.
 *
 * Returns: (transfer none): the value of #NcClusterAbundance:photoz property.
 */
NcClusterPhotoz *
nc_cluster_abundance_get_photoz (NcClusterAbundance *cad)
{
  return cad->photoz;
}

/**
 * nc_cluster_abundance_set_lnM_sigma0:
 * @cad: a #NcClusterAbundance.
 * @lnM_sigma0: value of #NcClusterAbundance:lnM-sigma0.
 *
 * Sets the value @lnM_sigma0 to the #NcClusterAbundance:lnM-sigma0 property.
 *
 */
void
nc_cluster_abundance_set_lnM_sigma0 (NcClusterAbundance *cad, const gdouble lnM_sigma0)
{
  g_assert (lnM_sigma0 > cad->lnM_sigma0);
  cad->lnM_sigma0 = lnM_sigma0;
}

/**
 * nc_cluster_abundance_get_lnM_sigma0:
 * @cad: a #NcClusterAbundance.
 *
 * Gets the value of the #NcClusterAbundance:lnM-sigma0 property.
 *
 * Returns: the value of #NcClusterAbundance:lnM-sigma0 property.
 */
gdouble
nc_cluster_abundance_get_lnM_sigma0 (NcClusterAbundance *cad)
{
  return cad->lnM_sigma0;
}

static void
_nc_cluster_abundance_dispose (GObject *object)
{
  /* TODO: Add deinitalization code here */
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (object);
  nc_mass_function_free (cad->mfp);
  nc_cluster_photoz_free (cad->photoz);
  nc_halo_bias_func_free (cad->mbiasf);
  G_OBJECT_CLASS (nc_cluster_abundance_parent_class)->dispose (object);
}

static void
_nc_cluster_abundance_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_cluster_abundance_parent_class)->finalize (object);
}

static void
_nc_cluster_abundance_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (object);
  g_return_if_fail (NC_IS_CLUSTER_ABUNDANCE (object));

  switch (prop_id)
  {
	case PROP_OPT:
	  cad->opt = g_value_get_flags (value);
	  break;
	case PROP_MASS_FUNCTION:
	  cad->mfp = g_value_get_object (value);
	  break;
	case PROP_MEANBIAS:
	  cad->mbiasf = g_value_get_object (value);
	  break;
	case PROP_ZI:
	  nc_cluster_abundance_set_zi (cad, g_value_get_double (value));
	  break;
	case PROP_ZF:
	  nc_cluster_abundance_set_zf (cad, g_value_get_double (value));
	  break;
	case PROP_LNMI:
	  nc_cluster_abundance_set_lnMi (cad, g_value_get_double (value));
	  break;
	case PROP_LNMF:
	  nc_cluster_abundance_set_lnMf (cad, g_value_get_double (value));
	  break;
	case PROP_PHOTOZ:
	  cad->photoz = g_value_get_object (value);
	  break;
	case PROP_LNMS0:
	  nc_cluster_abundance_set_lnM_sigma0 (cad, g_value_get_double (value));
	  break;
	default:
	  G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
	  break;
  }
}

static void
_nc_cluster_abundance_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcClusterAbundance *cad = NC_CLUSTER_ABUNDANCE (object);
  g_return_if_fail (NC_IS_CLUSTER_ABUNDANCE (object));

  switch (prop_id)
  {
	case PROP_OPT:
	  g_value_set_flags (value, cad->opt);
	  break;
	case PROP_MASS_FUNCTION:
	  g_value_set_object (value, cad->mfp);
	  break;
	case PROP_MEANBIAS:
	  g_value_set_object (value, cad->mbiasf);
	  break;
	case PROP_ZI:
	  g_value_set_double (value, nc_cluster_abundance_get_zi(cad));
	  break;
	case PROP_ZF:
	  g_value_set_double (value, nc_cluster_abundance_get_zf(cad));
	  break;
	case PROP_LNMI:
	  g_value_set_double (value, nc_cluster_abundance_get_lnMi(cad));
	  break;
	case PROP_LNMF:
	  g_value_set_double (value, nc_cluster_abundance_get_lnMf(cad));
	  break;
	case PROP_PHOTOZ:
	  g_value_set_object (value, cad->photoz);
	  break;
	case PROP_LNMS0:
	  g_value_set_double (value, nc_cluster_abundance_get_lnM_sigma0(cad));
	  break;
	default:
	  G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
	  break;
  }
}

static void
nc_cluster_abundance_class_init (NcClusterAbundanceClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->dispose = _nc_cluster_abundance_dispose;
  object_class->finalize = _nc_cluster_abundance_finalize;
  object_class->set_property = _nc_cluster_abundance_set_property;
  object_class->get_property = _nc_cluster_abundance_get_property;

  /**
   * NcClusterAbundance:options:
   *
   * FIXME
   */
  g_object_class_install_property (object_class,
                                   PROP_OPT,
                                   g_param_spec_flags ("options",
                                                       NULL,
                                                       "Options",
                                                       NC_TYPE_CLUSTER_ABUNDANCE_OPT,
                                                       NC_CLUSTER_ABUNDANCE_NONE,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterAbundance:mass-function:
   *
   * FIXME
   */
  g_object_class_install_property (object_class,
                                   PROP_MASS_FUNCTION,
                                   g_param_spec_object ("mass-function",
                                                        NULL,
                                                        "Mass Function",
                                                        NC_TYPE_MASS_FUNCTION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterAbundance:mean-bias:
   *
   * FIXME
   */
  g_object_class_install_property (object_class,
                                   PROP_MEANBIAS,
                                   g_param_spec_object ("mean-bias",
                                                        NULL,
                                                        "Mean Halo Bias Function",
                                                        NC_TYPE_HALO_BIAS_FUNC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcClusterAbundance:minimum-redshift:
   *
   * Minimum redshift value is restricted to [_NC_CLUSTER_ABUNDANCE_MIN_Z, G_MAXDOUBLE].
	 */
  g_object_class_install_property (object_class,
                                   PROP_ZI,
                                   g_param_spec_double ("minimum-redshift",
                                                        NULL,
                                                        "Minimum redshift",
                                                        _NC_CLUSTER_ABUNDANCE_MIN_Z, G_MAXDOUBLE, _NC_CLUSTER_ABUNDANCE_MIN_Z,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterAbundance:maximum-redshift:
   *
   * Maximum redshift value is restricted to [_NC_CLUSTER_ABUNDANCE_MIN_Z, G_MAXDOUBLE].
	 */
  g_object_class_install_property (object_class,
                                   PROP_ZF,
                                   g_param_spec_double ("maximum-redshift",
                                                        NULL,
                                                        "Maximum redshift",
                                                        _NC_CLUSTER_ABUNDANCE_MIN_Z, G_MAXDOUBLE, _NC_CLUSTER_ABUNDANCE_MIN_Z,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterAbundance:minimum-mass:
   *
   * Minimum mass value is restricted to [11.0 * ln10, 16.0 * ln10] in units of h^{-1} M_sun.
*/
  g_object_class_install_property (object_class,
                                   PROP_LNMI,
                                   g_param_spec_double ("minimum-mass",
                                                        NULL,
                                                        "Minimum mass",
                                                        (11.0 * M_LN10), (16.0 * M_LN10), (11.0 * M_LN10),
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterAbundance:maximum-mass:
   *
   * Maximum mass value is restricted to [11.0 * ln10, 16.0 * ln10] in units of h^{-1} M_sun.
*/
  g_object_class_install_property (object_class,
                                   PROP_LNMF,
                                   g_param_spec_double ("maximum-mass",
                                                        NULL,
                                                        "Maximum mass",
                                                        (11.0 * M_LN10), (16.0 * M_LN10), (11.0 * M_LN10),
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterAbundance:photoz:
   *
   * FIXME
	 */
  g_object_class_install_property (object_class,
                                   PROP_PHOTOZ,
                                   g_param_spec_object ("photoz",
                                                        NULL,
                                                        "Photoz",
                                                        NC_TYPE_CLUSTER_PHOTOZ,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcClusterAbundance:lnM-sigma0:
   *
   * Satandard deviation factor of a log-normal mass-observable distribution: sd = sigma0.
   */
  g_object_class_install_property (object_class,
                                   PROP_LNMS0,
                                   g_param_spec_double ("lnM-sigma0",
                                                        NULL,
                                                        "Mass-observable sigma0",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

void
nc_cluster_abundance_set_options (NcClusterAbundance *cad, NcClusterAbundanceOpt opt)
{
  if (cad->opt != opt)
  {
	cad->opt = opt;
	ncm_model_ctrl_force_update (cad->ctrl);
  }
}

void
nc_cluster_abundance_and_options (NcClusterAbundance *cad, NcClusterAbundanceOpt opt)
{
  if ((cad->opt | opt) != cad->opt)
  {
	cad->opt |= opt;
	ncm_model_ctrl_force_update (cad->ctrl);
  }
}

typedef struct _observables_integrand_data
{
  NcClusterAbundance *cad;
  NcHICosmo *model;
  gdouble zp;
  gdouble lnMobs;
  gpointer data;
} observables_integrand_data;

/* Mass observable Distributions */

static gdouble
_nc_cluster_abundance_lognormal_mass_dist (observables_integrand_data *md_data, gdouble lnM)
{
  const gdouble lnMobs = md_data->lnMobs;
  const gdouble sigma_lnM = md_data->cad->lnM_sigma0;
  const gdouble sqrt2_sigma = M_SQRT2 * sigma_lnM;
  const gdouble x = (lnMobs - lnM) / sqrt2_sigma;

  return M_2_SQRTPI / (2.0 * M_SQRT2) * exp (- x * x) / (sigma_lnM);
}

/* Discussion: This function is apropriated when the survey has a known physical limitation and you know the lower mass
 threshold that can be detected.
 static gdouble
 _nc_cluster_abundance_lognormal_mass_dist_Mtoinfinity (observables_integrand_data *md_data, gdouble lnM)
 {
   const gdouble lnMobs = md_data->lnMobs;
   const gdouble sigma_lnM = md_data->cad->lnM_sigma0;
   const gdouble sqrt2_sigma = M_SQRT2 * sigma_lnM;
   const gdouble x = (lnMobs - lnM) / sqrt2_sigma;

   return M_2_SQRTPI * M_SQRT1_2 * exp (- x * x) / (sigma_lnM * (1.0 + gsl_sf_erf ((lnM - LNM_MIN)/sqrt2_sigma)));
   }*/

/* Distribution with integration interval from 0 to infinity. Local = Standard deviation of each cluster. */
/* If the assumption of the above function is right, I will have to implement _nc_cluster_abundance_lognormal_mass_dist_Mtoinfinity_local. */
static gdouble
_nc_cluster_abundance_lognormal_mass_dist_local (observables_integrand_data *md_data, gdouble lnM, gdouble z)
{
  gsize i, j;
  const gdouble lnMobs = md_data->lnMobs;
  gdouble sigma_lnM, sqrt2_sigma, x;
  gsl_histogram2d_find (md_data->cad->sd_lnM, z, lnM, &i, &j);
  sigma_lnM = gsl_histogram2d_get (md_data->cad->sd_lnM, i, j);
  sqrt2_sigma = M_SQRT2 * sigma_lnM;
  x = (lnMobs - lnM) / sqrt2_sigma;
  if (sigma_lnM <= 0)
	return 0.0;

  //printf ("sigma_lnM = %.5g\n", sigma_lnM);
  return M_2_SQRTPI / (2.0 * M_SQRT2) * exp (- x * x) / (sigma_lnM);
}

/* Selection function: completeness and purity*/

/**
 * _nc_cluster_abundance_completeness:
 * @cad: pointer to #NcClusterAbundance
 * @lnM: is the logarithm (base e) of the mass.
   * @z: is the redshift.
 *
 * This function returns the completeness /f$ c /f$ of a catalog for a cluster with mass /f$ lnM /f$ and redshift /f$ z /f$.
 *
 * Returns: a gdouble which represents /f$ c /f$.
 */
static gdouble
_nc_cluster_abundance_completeness (NcClusterAbundance *cad, gdouble lnM, gdouble z)
{
  gsize i, j;
  gdouble c;
  gsl_histogram2d_find (cad->completeness, z, lnM, &i, &j);
  c = gsl_histogram2d_get (cad->completeness, i, j);
  if (c < 0)
	return 0.0;

  return c;
}

/**
 * _nc_cluster_abundance_one_over_purity:
 * @md_data: FIXME
 * @lnM: is the logarithm (base e) of the mass.
   * @z: is the redshift.
 *
 * This function returns one over purity /f$ \frac{1}{p} /f$ of a catalog for a cluster with mass /f$ lnM /f$ and redshift /f$ z /f$..
 *
 * Returns: a gdouble which represents /f$ \frac{1}{p} /f$.
 */
static gdouble
_nc_cluster_abundance_one_over_purity (NcClusterAbundance *cad, gdouble lnM, gdouble z)
{
  gsize i, j;
  gdouble p;
  gsl_histogram2d_find (cad->purity, z, lnM, &i, &j);
  p = gsl_histogram2d_get (cad->purity, i, j);
  if (p <= 0)
	return 0.0;

  return 1.0 / p;
}

#define SMALL_FACTOR (1.0 - GSL_DBL_EPSILON)

static void
_check_histogram_limits (gsl_histogram2d *hist2dim, gdouble lnM, gdouble z, gsize *i_lnM, gsize *i_z)
{
  gint ret;
  if (z == gsl_histogram2d_xmax (hist2dim)) z *= SMALL_FACTOR;
  if (lnM == gsl_histogram2d_ymax (hist2dim)) lnM *= SMALL_FACTOR;
  if ((ret = gsl_histogram2d_find (hist2dim, z, lnM, i_z, i_lnM)) != GSL_SUCCESS)
	g_error ("gsl_histogram2d_find: %s\n\tfind % 20.16g in (% 20.16g, % 20.16g) and % 20.16g in (% 20.16g, % 20.16g).", gsl_strerror (ret),
	         z,   gsl_histogram2d_xmin (hist2dim), gsl_histogram2d_xmax (hist2dim),
	         lnM, gsl_histogram2d_ymin (hist2dim), gsl_histogram2d_ymax (hist2dim));
}

static void
_check_histogram_limits_include_sup (gsl_histogram2d *hist2dim, gdouble lnM, gdouble z, gsize *i_lnM, gsize *i_z)
{
  gint ret;
  if (z != gsl_histogram2d_xmin (hist2dim)) z *= SMALL_FACTOR;
  if (lnM != gsl_histogram2d_ymin (hist2dim)) lnM *= SMALL_FACTOR;
  if ((ret = gsl_histogram2d_find (hist2dim, z, lnM, i_z, i_lnM)) != GSL_SUCCESS)
	g_error ("gsl_histogram2d_find: %s\n\tfind % 20.16g in (% 20.16g, % 20.16g) and % 20.16g in (% 20.16g, % 20.16g).", gsl_strerror (ret),
	         z,   gsl_histogram2d_xmin (hist2dim), gsl_histogram2d_xmax (hist2dim),
	         lnM, gsl_histogram2d_ymin (hist2dim), gsl_histogram2d_ymax (hist2dim));
}

/* Functions with respect to d2N/dzdlogM. */

typedef struct __encapsulated_function_args
{
  NcClusterAbundance *cad;
  NcHICosmo *model;
  gdouble z;
  gdouble lnM;
  gdouble dVdz;
  gpointer data;
} _encapsulated_function_args;

/*
 static void
 _generate_2Dspline_knots (NcClusterAbundance *cad, NcmModel *model, gdouble (*encapsulated_z)(gdouble, gpointer), gdouble (*encapsulated_lnM)(gdouble, gpointer), gdouble rel_error)
 {
   NcmSpline *s_z, *s_lnM;
   GArray *w;
   gsl_function F, G;
   _encapsulated_function_args args;

   args.cad = cad;
   args.model = model;
   args.z = (cad->zf + cad->zi) / 2.0;
   args.lnM = (cad->lnMf + cad->lnMi) / 2.0;

   F.function = encapsulated_lnM;
   F.params = &args;
   s_lnM = ncm_spline_new_function (_NC_CAD_SPLINE, NCM_SPLINE_FUNCTION_SPLINE, &F, cad->lnMi, cad->lnMf, 0, rel_error);

   G.function = encapsulated_z;
   G.params = &args;
   s_z = ncm_spline_new_function (_NC_CAD_SPLINE, NCM_SPLINE_FUNCTION_SPLINE, &G, cad->zi, cad->zf, 0, rel_error);

   w = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), s_lnM->x_array->len * s_z->x_array->len);
   g_array_set_size (w, s_lnM->x_array->len * s_z->x_array->len);

   //    for (i = 0; i < s_lnM->x_array->len; i++)
   //	  printf("no lnM = %.5g\n", g_array_index(s_lnM->x_array, gdouble, i));

   cad->d2NdzdlnM = ncm_spline2d_new (_NC_CAD_SPLINE, s_lnM->x_array, s_z->x_array, w, FALSE);

   g_array_free (s_lnM->y_array, TRUE); // Ignoring y_array, keeping only x_array with the knots.
   g_array_free (s_z->y_array, TRUE);

   ncm_spline_free (s_lnM, FALSE); // Free related to the spline not to the array.
   ncm_spline_free (s_z, FALSE);
   printf ("d2N/dzdlnM len z = %u len M = %u\n", cad->d2NdzdlnM->y_array->len, cad->d2NdzdlnM->x_array->len);
   }
   */

static void
_nc_cluster_abundance_generate_2Dspline_knots (NcClusterAbundance *cad, NcHICosmo *model, gdouble (*encapsulated_z)(gdouble, gpointer), gdouble (*encapsulated_lnM)(gdouble, gpointer), gdouble rel_error)
{
  gsl_function Fx, Fy;
  _encapsulated_function_args args;

  args.cad = cad;
  args.model = model;
  args.z = (cad->zf + cad->zi) / 2.0;
  args.lnM = (cad->lnMf + cad->lnMi) / 2.0;

  Fx.function = encapsulated_lnM;
  Fx.params = &args;

  Fy.function = encapsulated_z;
  Fy.params = &args;

  cad->d2NdzdlnM = ncm_spline2d_bicubic_notaknot_new ();
  ncm_spline2d_set_function (cad->d2NdzdlnM,
                             NCM_SPLINE_FUNCTION_SPLINE,
                             &Fx, &Fy, cad->lnMi, cad->lnMf, cad->zi, cad->zf, 1e-5);

}

static gdouble
_nc_cluster_abundance_photoz_integrand (gdouble z, gpointer params)
{
  observables_integrand_data *photoz_data = (observables_integrand_data *) params;
  NcClusterAbundance *cad = photoz_data->cad;

  /* In this case Mobs is the true mass, i.e., without uncertainty. */
  gdouble p_z_zr = nc_cluster_photoz_dist_eval (cad->photoz, NULL, photoz_data->zp, z);
  gdouble d2NdzdlnM = nc_mass_function_d2NdzdlnM (cad->mfp, photoz_data->model, photoz_data->lnMobs, z);
g_assert_not_reached ();
  //printf ("% 20.15g % 20.15g\n", p_z_zr, d2NdzdlnM);
  return p_z_zr * d2NdzdlnM;
}

static gdouble
_encapsulated_zfunction_photoz (gdouble z, gpointer p)
{
  _encapsulated_function_args *args = (_encapsulated_function_args *) p;

  return nc_cluster_abundance_d2NdzdlnM_photoz (args->cad, args->model, args->lnM, z);
}

static gdouble
_encapsulated_lnMfunction_photoz (gdouble lnM, gpointer p)
{
  _encapsulated_function_args *args = (_encapsulated_function_args *) p;

  return nc_cluster_abundance_d2NdzdlnM_photoz (args->cad, args->model, lnM, args->z);
}

/**
 * _nc_cluster_abundance_Mobs_integrand:
 * @lnM: a gdouble which is the logarithm base e of the mass.
 * @params: a gpointer.
 *
 * This function computes \f$ \frac{d^2N}{dzdlnM} \times P(\ln M^{obs}| \ln M) \f$
   * where \f$ P(M^{obs}| M) \f$ is the log-normal
	 * probability distribution of the observable mass \f$ M^{obs} \f$ with a fixed standard deviation.
	   *
 * Returns: a gdouble which is the integrand value at lnM.
 */
static gdouble
_nc_cluster_abundance_Mobs_integrand (gdouble lnM, gpointer params)
{
  observables_integrand_data *md_data = (observables_integrand_data *) params;
  NcClusterAbundance *cad = md_data->cad;

  /* In this case zp is the true redshift, i.e., without uncertainty. */
  /* We have to define which distribution is more apropriate to be used (take into account or not the physical limitation of a survey). */
  gdouble p_M_Mobs = _nc_cluster_abundance_lognormal_mass_dist (md_data, lnM);
  //gdouble p_M_Mobs = _nc_cluster_abundance_lognormal_mass_dist_Mtoinfinity (md_data, lnM);
  gdouble d2NdzdlnM = nc_mass_function_d2NdzdlnM (cad->mfp, md_data->model, lnM, md_data->zp);

  //printf ("% 20.8e % 20.8g % 20.15g % 20.15g\n", md_data->lnMobs, md_data->zp, p_M_Mobs, d2NdzdlnM);
  return p_M_Mobs * d2NdzdlnM;
}

static gdouble
_encapsulated_zfunction_Mobs (gdouble z, gpointer p)
{
  _encapsulated_function_args *args = (_encapsulated_function_args *) p;

  return nc_cluster_abundance_d2NdzdlnM_Mobs (args->cad, args->model, args->lnM, z);
}

static gdouble
_encapsulated_lnMfunction_Mobs (gdouble lnM, gpointer p)
{
  _encapsulated_function_args *args = (_encapsulated_function_args *) p;

  return nc_cluster_abundance_d2NdzdlnM_Mobs (args->cad, args->model, lnM, args->z);
}

static gdouble
_nc_cluster_abundance_photoz_Mobs_integrand (gdouble lnM, gdouble z, gpointer userdata)
{
  observables_integrand_data *obs_data = (observables_integrand_data *) userdata;
  NcClusterAbundance *cad = obs_data->cad;
  gdouble p_z_zr = nc_cluster_photoz_dist_eval (cad->photoz, NULL, obs_data->zp, z);
  gdouble p_M_Mobs = _nc_cluster_abundance_lognormal_mass_dist (obs_data, lnM);
  gdouble d2NdzdlnM = nc_mass_function_d2NdzdlnM (cad->mfp, obs_data->model, lnM, z);
g_assert_not_reached ();
  //printf ("% 20.8e % 20.8g % 20.15g % 20.15g\n", md_data->lnMobs, md_data->zp, p_M_Mobs, d2NdzdlnM);
  return p_z_zr * p_M_Mobs * d2NdzdlnM;
}

static gdouble
_encapsulated_zfunction_photoz_Mobs (gdouble z, gpointer p)
{
  _encapsulated_function_args *args = (_encapsulated_function_args *) p;

  return nc_cluster_abundance_d2NdzdlnM_photoz_Mobs (args->cad, args->model, args->lnM, z);
}

static gdouble
_encapsulated_lnMfunction_photoz_Mobs (gdouble lnM, gpointer p)
{
  _encapsulated_function_args *args = (_encapsulated_function_args *) p;

  return nc_cluster_abundance_d2NdzdlnM_photoz_Mobs (args->cad, args->model, lnM, args->z);
}

/**
 * _nc_cluster_abundance_Mobs_local_selection_integrand:
 * @lnM: is the logarithm (base e) of the mass.
   * @params: FIXME
 *
 * This function computes /f$ \frac{d^2N}{dzd\lnM} * P(\ln M^{obs}|\ln M, z) /f$.
 * Completeness and purity are included later in the integral.
 *
 * Returns: A double which correspond to the function above.
 */
static gdouble
_nc_cluster_abundance_Mobs_local_selection_integrand (gdouble lnM, gpointer params)
{
  observables_integrand_data *obs_data = (observables_integrand_data *) params;
  NcClusterAbundance *cad = obs_data->cad;
  gdouble p_Mobs_local = _nc_cluster_abundance_lognormal_mass_dist_local (obs_data, lnM, obs_data->zp);
  gdouble d2NdzdlnM = nc_mass_function_dcomoving_volume_dzdomega (cad->mfp, obs_data->model, obs_data->zp) * nc_mass_function (cad->mfp, obs_data->model, lnM, obs_data->zp);

  //printf ("lnMobs =% 10.5e p_M = % 10.5g d2N = % 10.5g\n", obs_data->lnMobs, p_Mobs_local, d2NdzdlnM);
  return p_Mobs_local * d2NdzdlnM;
}

static gdouble
_nc_cluster_abundance_dNdz_Mobs_local_selection_integrand (gdouble lnMobs, gpointer params)
{
  observables_integrand_data *obs_data = (observables_integrand_data *) params;
  NcClusterAbundance *cad = obs_data->cad;

  gdouble d2NdzdlnMobs = nc_cluster_abundance_d2NdzdlnM_Mobs_local_selection (cad, obs_data->model, lnMobs, obs_data->zp);

  //printf ("lnMobs =% 10.5e d2N = % 10.5g\n", lnMobs, d2NdzdlnM);

  return d2NdzdlnMobs;
}

static gdouble
_nc_cluster_abundance_dNdz_photoz_Mobs_local_selection_integrand (gdouble lnMobs, gpointer params)
{
  observables_integrand_data *obs_data = (observables_integrand_data *) params;
  NcClusterAbundance *cad = obs_data->cad;

  gdouble d2NdzdlnMobs = nc_cluster_abundance_d2NdzdlnM_photoz_Mobs_local_selection (cad, obs_data->model, lnMobs, obs_data->zp);

  //printf ("lnMobs =% 10.5e d2N = % 10.5g\n", lnMobs, d2NdzdlnM);

  return d2NdzdlnMobs;
}

static gdouble
_nc_cluster_abundance_N_Mobs_local_selection_integrand (gdouble lnMobs, gdouble z, gpointer params)
{
  observables_integrand_data *obs_data = (observables_integrand_data *) params;
  NcClusterAbundance *cad = obs_data->cad;

  gdouble d2NdzdlnMobs = nc_cluster_abundance_d2NdzdlnM_Mobs_local_selection (cad, obs_data->model, lnMobs, z);

  //printf ("z = %.5g lnM = %.5g d2N = %.5g\n", z, lnMobs, d2NdzdlnMobs);

  return d2NdzdlnMobs;
}

static gdouble
_nc_cluster_abundance_N_photoz_Mobs_local_selection_integrand (gdouble lnMobs, gdouble zp, gpointer params)
{
  observables_integrand_data *obs_data = (observables_integrand_data *) params;
  NcClusterAbundance *cad = obs_data->cad;

  gdouble d2NdzdlnMobs = nc_cluster_abundance_d2NdzdlnM_photoz_Mobs_local_selection (cad, obs_data->model, lnMobs, zp);

  //printf ("lnMobs =% 10.5g  zp = %.5g d2N = % 10.5g\n", lnMobs, zp, d2NdzdlnMobs);
  printf(".");fflush(stdout);

  return d2NdzdlnMobs;
}

/*
 static gdouble
 _encapsulated_zfunction_Mobs_local_selection (gdouble z, gpointer p)
 {
   _encapsulated_function_args *args = (_encapsulated_function_args *) p;

   return nc_cluster_abundance_d2NdzdlnM_Mobs_local_selection (args->cad, args->model, args->lnM, z);
}
*/
/*
 static gdouble
 _encapsulated_lnMfunction_Mobs_local_selection (gdouble lnM, gpointer p)
 {
   _encapsulated_function_args *args = (_encapsulated_function_args *) p;

   return nc_cluster_abundance_d2NdzdlnM_Mobs_local_selection (args->cad, args->model, lnM, args->z);
   }
   */

/*
static void
_nc_cluster_abundance_Mobs_local_selection_prepare (NcClusterAbundance *cad, NcHICosmo *model)
{
  gint i, j;
  printf("Aqui 1\n");

  if (cad->d2NdzdlnM == NULL)
  {
	_nc_cluster_abundance_generate_2Dspline_knots (cad, model, &_encapsulated_zfunction_Mobs_local_selection, &_encapsulated_lnMfunction_Mobs_local_selection, 1.0e-5);
  }

  printf ("n de nos z: %u lnM: %u\n", D2NDZDLNM_Z(cad)->len, D2NDZDLNM_LNM(cad)->len);
  for (i = 0; i < D2NDZDLNM_Z(cad)->len; i++)
  {
	const gdouble z = g_array_index (D2NDZDLNM_Z(cad), gdouble, i);
	for (j = 0; j < D2NDZDLNM_LNM(cad)->len; j++)
	{
	  const gdouble lnMobs = g_array_index (D2NDZDLNM_LNM(cad), gdouble, j);
	  const gdouble d2NdzdlnM_ij = nc_cluster_abundance_d2NdzdlnM_Mobs_local_selection (cad, model, lnMobs, z);
	  g_array_index (D2NDZDLNM_VAL(cad), gdouble, i * D2NDZDLNM_LNM(cad)->len + j) = d2NdzdlnM_ij;
	}
  }
  ncm_spline2d_prepare (cad->d2NdzdlnM);
}
*/

/**
 * _nc_cluster_abundance_photoz_Mobs_local_integrand:
 * @lnM: logarithm base e of mass.
 * @z: redshift.
 * @params: a pointer.
 *
 * This function returns the integrand d2N/dzdlnM * P(lnMobs|lnM) *
 * P(zphot|z) / (Survey area), when we use a matching catalog. Therefore
 * the mass probability distribution has a standard deviation which is function
 * of z and lnM. As the completeness and purity are constant functions (in which
 * in mass and redshift bins), they are included after the integral is performed.
 * See #nc_cluster_abundance_d2NdzdlnM_Mobs_local_selection and
 * #nc_cluster_abundance_d2NdzdlnM_purity_val for example.
 *
 * Returns: a gdouble which represents d2N/dzdlnM * P(lnMobs|lnM) *
 * P(zphot|z) / (Survey area).
 */
static gdouble
_nc_cluster_abundance_photoz_Mobs_local_integrand (gdouble lnM, gdouble z, gpointer params)
{
  observables_integrand_data *obs_data = (observables_integrand_data *) params;
  NcClusterAbundance *cad = obs_data->cad;
  gdouble p_zp = nc_cluster_photoz_dist_eval (cad->photoz, NULL, obs_data->zp, z);
  gdouble p_Mobs_local = _nc_cluster_abundance_lognormal_mass_dist_local (obs_data, lnM, z);
  gdouble d2NdzdlnM = nc_mass_function_dcomoving_volume_dzdomega (cad->mfp, obs_data->model, z) * nc_mass_function (cad->mfp, obs_data->model, lnM, z);

  //printf ("% 20.8e % 20.15g % 20.15g % 20.15g\n", z, p_zp, p_Mobs_local, d2NdzdlnM);
  return p_zp * p_Mobs_local * d2NdzdlnM;
}

/*
static gdouble
_encapsulated_zfunction_photoz_Mobs_local_selection (gdouble z, gpointer p)
{
	_encapsulated_function_args *args = (_encapsulated_function_args *) p;

  return nc_cluster_abundance_d2NdzdlnM_photoz_Mobs_local_selection (args->cad, args->model, args->lnM, z);
}
*/
/*
static gdouble
_encapsulated_lnMfunction_photoz_Mobs_local_selection (gdouble lnM, gpointer p)
{
	_encapsulated_function_args *args = (_encapsulated_function_args *) p;

  return nc_cluster_abundance_d2NdzdlnM_photoz_Mobs_local_selection (args->cad, args->model, lnM, args->z);
}
*/

/**
 * _nc_cluster_abundance_photoz_Mobs_local_selection_unbinned_prepare:
 * @cad: a pointer to #NcClusterAbundance.
 * @model: a pointer to #NcHICosmo.
 *
 * This function prepares a bidimensional spline of d2N/dz_phot dlnMobs
 * taking into account the completeness. It will be used when one performs an
 * unbinned analysis in both redshift and mass and will be called by
 * #nc_cluster_abundance_d2NdzdlnM_purity_val, where the purity is included.
 *
 *
static void
_nc_cluster_abundance_photoz_Mobs_local_selection_unbinned_prepare (NcClusterAbundance *cad, NcHICosmo *model)
{
  gint i, j;

  if (cad->d2NdzdlnM == NULL)
  {
	_nc_cluster_abundance_generate_2Dspline_knots (cad, model, &_encapsulated_zfunction_photoz_Mobs_local_selection, &_encapsulated_lnMfunction_photoz_Mobs_local_selection, 1.0e-5);
  }

  //printf ("n de nos z: %u lnM: %u\n", D2NDZDLNM_Z(cad)->len, D2NDZDLNM_LNM(cad)->len);
  for (i = 0; i < D2NDZDLNM_Z(cad)->len; i++)
  {
	const gdouble zp = g_array_index (D2NDZDLNM_Z(cad), gdouble, i);
	for (j = 0; j < D2NDZDLNM_LNM(cad)->len; j++)
	{
	  const gdouble lnMobs = g_array_index(D2NDZDLNM_LNM(cad), gdouble, j);
	  const gdouble d2NdzdlnM_ij = nc_cluster_abundance_d2NdzdlnM_photoz_Mobs_local_selection (cad, model, lnMobs, zp);
	  g_array_index (D2NDZDLNM_VAL(cad), gdouble, i * D2NDZDLNM_LNM(cad)->len + j) = d2NdzdlnM_ij;
	}
  }
  ncm_spline2d_prepare (cad->d2NdzdlnM);
}
*/

static gdouble _nc_cluster_abundance_d2NdzdlnM_spline_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z);
static gdouble _nc_cluster_abundance_dNdz_spline_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble z);
static gdouble _nc_cluster_abundance_N_spline_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble zl, gdouble zu);
static gdouble _nc_cluster_abundance_d2NdzdlnM_obs_spline_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z);
static gdouble _nc_cluster_abundance_dNdz_obs_spline_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble z);
static gdouble _nc_cluster_abundance_N_obs_spline_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble zl, gdouble zu);
static gdouble _nc_cluster_abundance_d2NdzdlnM_Mobs_local_selection_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs, gdouble z);
static gdouble _nc_cluster_abundance_dNdz_Mobs_local_selection_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble z);
static gdouble _nc_cluster_abundance_N_Mobs_local_selection_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble zl, gdouble zu);
static gdouble _nc_cluster_abundance_d2NdzdlnM_photoz_Mobs_local_selection_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs, gdouble zp);
static gdouble _nc_cluster_abundance_dNdz_photoz_Mobs_local_selection_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble z);
static gdouble _nc_cluster_abundance_N_photoz_Mobs_local_selection_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble zl, gdouble zu);

void
nc_cluster_abundance_prepare (NcClusterAbundance *cad, NcHICosmo *model)
{
  const gint test_mask = NC_CLUSTER_ABUNDANCE_PHOTOZ |
	NC_CLUSTER_ABUNDANCE_MOBS |
	NC_CLUSTER_ABUNDANCE_MOBS_LOCAL |
	NC_CLUSTER_ABUNDANCE_COMPLETENESS |
	NC_CLUSTER_ABUNDANCE_PURITY;
  const gint opt_mask = cad->opt & test_mask;
  GTimer *bench = g_timer_new ();

  if (cad->optimize)
	nc_mass_function_d2NdzdlnM_optimize (cad->mfp, model, cad->lnMi, cad->lnMf, cad->zi, cad->zf);

  switch (opt_mask)
  {
	case NC_CLUSTER_ABUNDANCE_PHOTOZ:
	  cad->d2NdzdlnM_val = &_nc_cluster_abundance_d2NdzdlnM_obs_spline_val;
	  cad->dNdz_val = &_nc_cluster_abundance_dNdz_obs_spline_val;
	  cad->N_val = &_nc_cluster_abundance_N_obs_spline_val;
	  break;
	case NC_CLUSTER_ABUNDANCE_MOBS:
	  cad->d2NdzdlnM_val = &_nc_cluster_abundance_d2NdzdlnM_obs_spline_val;
	  cad->dNdz_val = &_nc_cluster_abundance_dNdz_obs_spline_val;
	  cad->N_val = &_nc_cluster_abundance_N_obs_spline_val;
	  break;
	case NC_CLUSTER_ABUNDANCE_PHOTOZ | NC_CLUSTER_ABUNDANCE_MOBS:
	  cad->d2NdzdlnM_val = &_nc_cluster_abundance_d2NdzdlnM_obs_spline_val;
	  cad->dNdz_val = &_nc_cluster_abundance_dNdz_obs_spline_val;
	  cad->N_val = &_nc_cluster_abundance_N_obs_spline_val;
	  break;
	case NC_CLUSTER_ABUNDANCE_MOBS_LOCAL | NC_CLUSTER_ABUNDANCE_COMPLETENESS | NC_CLUSTER_ABUNDANCE_PURITY:
	  cad->d2NdzdlnM_val = &_nc_cluster_abundance_d2NdzdlnM_Mobs_local_selection_val;
	  cad->dNdz_val = &_nc_cluster_abundance_dNdz_Mobs_local_selection_val;
	  cad->N_val = &_nc_cluster_abundance_N_Mobs_local_selection_val;
	  break;
	case NC_CLUSTER_ABUNDANCE_PHOTOZ | NC_CLUSTER_ABUNDANCE_MOBS_LOCAL | NC_CLUSTER_ABUNDANCE_COMPLETENESS | NC_CLUSTER_ABUNDANCE_PURITY:
	  cad->d2NdzdlnM_val = &_nc_cluster_abundance_d2NdzdlnM_photoz_Mobs_local_selection_val;
	  cad->dNdz_val = &_nc_cluster_abundance_dNdz_photoz_Mobs_local_selection_val;
	  cad->N_val = &_nc_cluster_abundance_N_photoz_Mobs_local_selection_val;
	  break;
	case 0:
	  cad->d2NdzdlnM_val = &_nc_cluster_abundance_d2NdzdlnM_spline_val;
	  cad->dNdz_val = &_nc_cluster_abundance_dNdz_spline_val;
	  cad->N_val = &_nc_cluster_abundance_N_spline_val;
	  break;
	default:
	  g_assert_not_reached ();
	  break;
  }
  cad->norma = nc_cluster_abundance_N_val (cad, model, cad->lnMi, cad->lnMf, cad->zi, cad->zf);
  cad->log_norma = log (cad->norma);

  ncm_model_ctrl_update (cad->ctrl, NCM_MODEL(model));
  //printf ("# preparado %.15f\n", g_timer_elapsed (bench, NULL));
  g_timer_destroy (bench);
  //printf ("HERE! %.15g [(% 20.15g % 20.15g) (% 20.15g % 20.15g)] \n", cad->norma, cad->lnMi, cad->lnMf, cad->zi, cad->zf);
}

gdouble
_nc_cad_inv_dNdz_convergence_f (gdouble n, gdouble epsilon)
{
  return -log1p (epsilon - n);
}

static gdouble
_nc_cad_inv_dNdz_convergence_f_onemn (gdouble onemn, gdouble epsilon)
{
  return -log (epsilon + onemn);
}

/**
 * nc_cluster_abundance_prepare_inv_dNdz:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 *
 * This function prepares a bidimensional spline...
 *
 */
void
nc_cluster_abundance_prepare_inv_dNdz (NcClusterAbundance *cad, NcHICosmo *model)
{
  gint i, j;
  gdouble z0 = cad->zi;
  gint middle = cad->inv_z->len / 2;
  g_assert (cad->zi != 0);

  cad->z_epsilon = nc_cluster_abundance_N_val (cad, model, cad->lnMi, cad->lnMf,
                                               cad->zi + (cad->zf - cad->zi) / (cad->inv_z->len - 1.0) * (cad->inv_z->len - 2.0),
                                               cad->zf) / cad->norma;
  cad->lnM_epsilon = nc_cluster_abundance_dNdz_val (cad, model,
                                                    cad->lnMi + (cad->lnMf - cad->lnMi) / (cad->inv_lnM->len - 1.0) * (cad->inv_lnM->len - 2.0),
                                                    cad->lnMf, cad->zf) /
	nc_cluster_abundance_dNdz_val (cad, model, cad->lnMi, cad->lnMf, cad->zf);
  //printf ("z epsilon % 20.15g lnM epsilon % 20.15g\n", cad->z_epsilon, cad->lnM_epsilon);
  {
	gdouble zm = cad->zi + (cad->zf - cad->zi) / (cad->inv_z->len - 1.0) * middle;
	nc_cluster_abundance_prepare_inv_dNdlnM_z (cad, model, zm);

	ncm_vector_set (cad->inv_lnM_z->xv, 0, _nc_cad_inv_dNdz_convergence_f (0.0, cad->lnM_epsilon)); // 0.0;
	ncm_matrix_set (cad->inv_lnM_z->zm, middle, 0, cad->lnMi);

	//printf ("%zu % 20.15g  %zd % 20.15g\n", cad->inv_z->len, zm, ncm_vector_len(cad->inv_lnM_z->xv), cad->lnMi);
	//printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", 0, ncm_vector_get (cad->inv_lnM->xv, 0), ncm_spline_eval (cad->inv_lnM, ncm_vector_get (cad->inv_lnM->xv, 0)));
	for (j = 1; j < ncm_vector_len(cad->inv_lnM_z->xv) - 1; j++)
	{
	  gdouble u2 = ncm_vector_get (cad->inv_lnM->xv, j);
	  ncm_vector_set (cad->inv_lnM_z->xv, j, u2);
	  ncm_matrix_set (cad->inv_lnM_z->zm, middle, j, ncm_spline_eval (cad->inv_lnM, u2));
	  //printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", j, u2, ncm_spline_eval (cad->inv_lnM, u2));
	}
	ncm_vector_set (cad->inv_lnM_z->xv, j, _nc_cad_inv_dNdz_convergence_f_onemn (0.0, cad->lnM_epsilon));
	ncm_matrix_set (cad->inv_lnM_z->zm, middle, j, cad->lnMf);
	//printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", j, ncm_vector_get (cad->inv_lnM_z->xv, j), ncm_spline_eval (cad->inv_lnM, ncm_vector_get (cad->inv_lnM_z->xv, j)));
  }

  nc_cluster_abundance_prepare_inv_dNdlnM_z (cad, model, z0);
  ncm_matrix_set (cad->inv_lnM_z->zm, 0, 0, cad->lnMi);

  //printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", 0, ncm_vector_get (cad->inv_lnM_z->xv, 0), ncm_spline_eval (cad->inv_lnM, ncm_vector_get (cad->inv_lnM_z->xv, 0)));
  for (j = 1; j < ncm_vector_len(cad->inv_lnM_z->xv) - 1; j++)
  {
	gdouble u2 = ncm_vector_get (cad->inv_lnM_z->xv, j);
	ncm_matrix_set (cad->inv_lnM_z->zm, 0, j, ncm_spline_eval (cad->inv_lnM, u2));
	//printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", j, u2, ncm_spline_eval (cad->inv_lnM, u2));
  }
  ncm_matrix_set (cad->inv_lnM_z->zm, 0, j, cad->lnMf);
  //printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", j, ncm_vector_get (cad->inv_lnM_z->xv, j), ncm_spline_eval (cad->inv_lnM, ncm_vector_get (cad->inv_lnM_z->xv, j)));
  //printf ("\n\n");

  {
	gdouble nztot = 0.0;
	gdouble f = _nc_cad_inv_dNdz_convergence_f (0.0, cad->z_epsilon);
	ncm_vector_set (cad->inv_z->xv, 0, f);
	ncm_vector_set (cad->inv_z->yv, 0, z0);
	//printf ("# f % 20.15g z % 20.15g\n", f, z0);

	for (i = 1; i < cad->inv_z->len; i++)
	{
	  gdouble z1 = cad->zi + (cad->zf - cad->zi) / (cad->inv_z->len - 1.0) * i;
	  if (nztot < 0.99)
	  {
		gdouble delta = nc_cluster_abundance_N_val (cad, model, cad->lnMi, cad->lnMf, z0, z1) / cad->norma;
		nztot += delta;
		f = _nc_cad_inv_dNdz_convergence_f (nztot, cad->z_epsilon);
	  }
	  else
	  {
		gdouble onemn = nc_cluster_abundance_N_val (cad, model, cad->lnMi, cad->lnMf, z1, cad->zf) / cad->norma;
		f = _nc_cad_inv_dNdz_convergence_f_onemn (onemn, cad->z_epsilon);
	  }
	  ncm_vector_set (cad->inv_z->xv, i, f);
	  ncm_vector_set (cad->inv_z->yv, i, z1);
	  //printf ("# f % 20.15g z % 20.15g\n", f, z1);
	  z0 = z1;
	  if (i == middle)
		continue;

	  nc_cluster_abundance_prepare_inv_dNdlnM_z (cad, model, z1);

	  ncm_matrix_set (cad->inv_lnM_z->zm, i, 0, cad->lnMi);
	  //printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", 0, ncm_vector_get (cad->inv_lnM_z->xv, 0), ncm_spline_eval (cad->inv_lnM, ncm_vector_get (cad->inv_lnM_z->xv, 0)));
	  for (j = 1; j < ncm_vector_len(cad->inv_lnM_z->xv) - 1; j++)
	  {
		gdouble u2 = ncm_vector_get (cad->inv_lnM_z->xv, j);
		ncm_matrix_set (cad->inv_lnM_z->zm, i, j, ncm_spline_eval (cad->inv_lnM, u2));
		//printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", j, u2, ncm_spline_eval (cad->inv_lnM, u2));
	  }
	  ncm_matrix_set (cad->inv_lnM_z->zm, i, j, cad->lnMf);
	  //printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", j, ncm_vector_get (cad->inv_lnM_z->xv, j), ncm_spline_eval (cad->inv_lnM, ncm_vector_get (cad->inv_lnM_z->xv, j)));
	  //printf ("\n\n");
	}
  }
  ncm_spline2d_prepare (cad->inv_lnM_z);
  ncm_spline_prepare (cad->inv_z);
}

/**
 * nc_cluster_abundance_prepare_inv_dNdz_no_obs:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 *
 * This function prepares a bidimensional spline...
 *
 */
void
nc_cluster_abundance_prepare_inv_dNdz_no_obs (NcClusterAbundance *cad, NcHICosmo *model)
{
  gint i, j;
  gdouble z0 = cad->zi;
  gint middle = cad->inv_z->len / 2;
  gboolean use_spline = TRUE;
  NcMassFunctionSplineOptimize sp_optimize = NC_MASS_FUNCTION_SPLINE_Z;
  g_assert (cad->zi != 0);

  cad->z_epsilon = nc_mass_function_N (cad->mfp, model, cad->lnMi, cad->lnMf, cad->zi + (cad->zf - cad->zi) / (cad->inv_z->len - 1.0) *
                                       (cad->inv_z->len - 2.0), cad->zf, sp_optimize) / cad->norma;

  cad->lnM_epsilon = nc_mass_function_dNdz (cad->mfp, model, cad->lnMi + (cad->lnMf - cad->lnMi) / (cad->inv_lnM->len - 1.0) *
                                            (cad->inv_lnM->len - 2.0), cad->lnMf, cad->zf, use_spline) /
	nc_mass_function_dNdz (cad->mfp, model, cad->lnMi, cad->lnMf, cad->zf, use_spline);

  printf ("Aqui!!!\n");
  //printf ("z epsilon % 20.15g lnM epsilon % 20.15g\n", cad->z_epsilon, cad->lnM_epsilon);
  {
	gdouble zm = cad->zi + (cad->zf - cad->zi) / (cad->inv_z->len - 1.0) * middle;
	nc_cluster_abundance_prepare_inv_dNdlnM_z_no_obs (cad, model, zm);

	ncm_vector_set (cad->inv_lnM_z->xv, 0, _nc_cad_inv_dNdz_convergence_f (0.0, cad->lnM_epsilon)); // 0.0;
	ncm_matrix_set (cad->inv_lnM_z->zm, middle, 0, cad->lnMi);

	//printf ("%zu % 20.15g  %zd % 20.15g\n", cad->inv_z->len, zm, ncm_vector_len(cad->inv_lnM_z->xv), cad->lnMi);
	//printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", 0, ncm_vector_get (cad->inv_lnM->xv, 0), ncm_spline_eval (cad->inv_lnM, ncm_vector_get (cad->inv_lnM->xv, 0)));
	for (j = 1; j < ncm_vector_len(cad->inv_lnM_z->xv) - 1; j++)
	{
	  gdouble u2 = ncm_vector_get (cad->inv_lnM->xv, j);
	  ncm_vector_set (cad->inv_lnM_z->xv, j, u2);
	  ncm_matrix_set (cad->inv_lnM_z->zm, middle, j, ncm_spline_eval (cad->inv_lnM, u2));
	  //printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", j, u2, ncm_spline_eval (cad->inv_lnM, u2));
	}
	ncm_vector_set (cad->inv_lnM_z->xv, j, _nc_cad_inv_dNdz_convergence_f_onemn (0.0, cad->lnM_epsilon));
	ncm_matrix_set (cad->inv_lnM_z->zm, middle, j, cad->lnMf);
	//printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", j, ncm_vector_get (cad->inv_lnM_z->xv, j), ncm_spline_eval (cad->inv_lnM, ncm_vector_get (cad->inv_lnM_z->xv, j)));
  }

  nc_cluster_abundance_prepare_inv_dNdlnM_z_no_obs (cad, model, z0);
  ncm_matrix_set (cad->inv_lnM_z->zm, 0, 0, cad->lnMi);

  //printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", 0, ncm_vector_get (cad->inv_lnM_z->xv, 0), ncm_spline_eval (cad->inv_lnM, ncm_vector_get (cad->inv_lnM_z->xv, 0)));
  for (j = 1; j < ncm_vector_len(cad->inv_lnM_z->xv) - 1; j++)
  {
	gdouble u2 = ncm_vector_get (cad->inv_lnM_z->xv, j);
	ncm_matrix_set (cad->inv_lnM_z->zm, 0, j, ncm_spline_eval (cad->inv_lnM, u2));
	//printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", j, u2, ncm_spline_eval (cad->inv_lnM, u2));
  }
  ncm_matrix_set (cad->inv_lnM_z->zm, 0, j, cad->lnMf);
  //printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", j, ncm_vector_get (cad->inv_lnM_z->xv, j), ncm_spline_eval (cad->inv_lnM, ncm_vector_get (cad->inv_lnM_z->xv, j)));
  //printf ("\n\n");

  {
	gdouble nztot = 0.0;
	gdouble f = _nc_cad_inv_dNdz_convergence_f (0.0, cad->z_epsilon);
	ncm_vector_set (cad->inv_z->xv, 0, f);
	ncm_vector_set (cad->inv_z->yv, 0, z0);
	//printf ("# f % 20.15g z % 20.15g\n", f, z0);

	for (i = 1; i < cad->inv_z->len; i++)
	{
	  gdouble z1 = cad->zi + (cad->zf - cad->zi) / (cad->inv_z->len - 1.0) * i;
	  if (nztot < 0.99)
	  {
		gdouble delta = nc_mass_function_N (cad->mfp, model, cad->lnMi, cad->lnMf, z0, z1, sp_optimize) / cad->norma;
		nztot += delta;
		f = _nc_cad_inv_dNdz_convergence_f (nztot, cad->z_epsilon);
	  }
	  else
	  {
		gdouble onemn = nc_mass_function_N (cad->mfp, model, cad->lnMi, cad->lnMf, z1, cad->zf, sp_optimize) / cad->norma;
		f = _nc_cad_inv_dNdz_convergence_f_onemn (onemn, cad->z_epsilon);
	  }
	  ncm_vector_set (cad->inv_z->xv, i, f);
	  ncm_vector_set (cad->inv_z->yv, i, z1);
	  //printf ("# f % 20.15g z % 20.15g\n", f, z1);
	  z0 = z1;
	  if (i == middle)
		continue;

	  nc_cluster_abundance_prepare_inv_dNdlnM_z_no_obs (cad, model, z1);

	  ncm_matrix_set (cad->inv_lnM_z->zm, i, 0, cad->lnMi);
	  //printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", 0, ncm_vector_get (cad->inv_lnM_z->xv, 0), ncm_spline_eval (cad->inv_lnM, ncm_vector_get (cad->inv_lnM_z->xv, 0)));
	  for (j = 1; j < ncm_vector_len(cad->inv_lnM_z->xv) - 1; j++)
	  {
		gdouble u2 = ncm_vector_get (cad->inv_lnM_z->xv, j);
		ncm_matrix_set (cad->inv_lnM_z->zm, i, j, ncm_spline_eval (cad->inv_lnM, u2));
		//printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", j, u2, ncm_spline_eval (cad->inv_lnM, u2));
	  }
	  ncm_matrix_set (cad->inv_lnM_z->zm, i, j, cad->lnMf);
	  //printf ("# eval[%d] u2 lnM % 20.15g % 20.15g\n", j, ncm_vector_get (cad->inv_lnM_z->xv, j), ncm_spline_eval (cad->inv_lnM, ncm_vector_get (cad->inv_lnM_z->xv, j)));
	  //printf ("\n\n");
	}
  }
  ncm_spline2d_prepare (cad->inv_lnM_z);
  ncm_spline_prepare (cad->inv_z);
}

/**
 * nc_cluster_abundance_prepare_inv_dNdlnM_z:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @z: redshift.
 *
 * This function prepares a spline where the x array corresponds to the value
 * of \int_lnM0 ^lnM1 d2N/dzdlnM dM/ \int_lnMi^lnMf dN/dz dM given a redshift z
 * and the y array contains the values of logarithms base e of the mass.
 * It is used to generate a sample of lnM values.
 *
 */
void
nc_cluster_abundance_prepare_inv_dNdlnM_z (NcClusterAbundance *cad, NcHICosmo *model, gdouble z)
{
  gdouble dNdz = nc_cluster_abundance_dNdz_val (cad, model, cad->lnMi, cad->lnMf, z);
  gdouble lnM0 = cad->lnMi;
  gdouble ntot = 0.0;
  gdouble f = _nc_cad_inv_dNdz_convergence_f (0.0, cad->lnM_epsilon);
  gdouble Delta;
  gint i;
  g_assert (z > 0.0);

  ncm_vector_set (cad->inv_lnM->xv, 0, f); //dNdz_u2 / dNdz;
  ncm_vector_set (cad->inv_lnM->yv, 0, lnM0);

  for (i = 1; i < cad->inv_lnM->len; i++)
  {
	gdouble lnM1 = cad->lnMi + (cad->lnMf - cad->lnMi) / (cad->inv_lnM->len - 1.0) * i;
	Delta = nc_cluster_abundance_dNdz_val (cad, model, lnM0, lnM1, z) / dNdz;
	ntot += Delta;
	if (ntot > 0.99)
	  break;
	else
	  f = _nc_cad_inv_dNdz_convergence_f (ntot, cad->lnM_epsilon);
	ncm_vector_set (cad->inv_lnM->xv, i, f); //dNdz_u2 / dNdz;
	ncm_vector_set (cad->inv_lnM->yv, i, lnM1);
	//printf ("prep[%d] % 20.15g % 20.15g % 20.15g % 20.15g\n", i, ntot, f, lnM1, z);
	lnM0 = lnM1;
  }

  for (; i < cad->inv_lnM->len; i++)
  {
	gdouble lnM1 = cad->lnMi + (cad->lnMf - cad->lnMi) / (cad->inv_lnM->len - 1.0) * i;
	gdouble onemn = nc_cluster_abundance_dNdz_val (cad, model, lnM1, cad->lnMf, z) / dNdz;
	f = _nc_cad_inv_dNdz_convergence_f_onemn (onemn, cad->lnM_epsilon);
	ncm_vector_set (cad->inv_lnM->xv, i, f); //dNdz_u2 / dNdz;
	ncm_vector_set (cad->inv_lnM->yv, i, lnM1);
	//printf ("prep[%d] % 20.15g % 20.15g % 20.15g % 20.15g\n", i, onemn, f, lnM1, z);
  }

  ncm_spline_prepare (cad->inv_lnM);
}

/**
 * nc_cluster_abundance_prepare_inv_dNdlnM_z_no_obs:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @z: redshift.
 *
 * This function prepares a spline where the x array corresponds to the value
 * of \int_lnM0 ^lnM1 d2N/dzdlnM dM/ \int_lnMi^lnMf dN/dz dM given a redshift z
 * and the y array contains the values of logarithms base e of the mass.
 * It is used to generate a sample of lnM values.
 *
 */
void
nc_cluster_abundance_prepare_inv_dNdlnM_z_no_obs (NcClusterAbundance *cad, NcHICosmo *model, gdouble z)
{
  gboolean use_spline = TRUE;
  gdouble dNdz = nc_mass_function_dNdz (cad->mfp, model, cad->lnMi, cad->lnMf, z, use_spline);
  gdouble lnM0 = cad->lnMi;
  gdouble ntot = 0.0;
  gdouble f = _nc_cad_inv_dNdz_convergence_f (0.0, cad->lnM_epsilon);
  gdouble Delta;
  gint i;
  g_assert (z > 0.0);

  ncm_vector_set (cad->inv_lnM->xv, 0, f); //dNdz_u2 / dNdz;
  ncm_vector_set (cad->inv_lnM->yv, 0, lnM0);

  for (i = 1; i < cad->inv_lnM->len; i++)
  {
	gdouble lnM1 = cad->lnMi + (cad->lnMf - cad->lnMi) / (cad->inv_lnM->len - 1.0) * i;
	Delta = nc_mass_function_dNdz (cad->mfp, model, lnM0, lnM1, z, use_spline) / dNdz;
	ntot += Delta;
	if (ntot > 0.99)
	  break;
	else
	  f = _nc_cad_inv_dNdz_convergence_f (ntot, cad->lnM_epsilon);
	ncm_vector_set (cad->inv_lnM->xv, i, f); //dNdz_u2 / dNdz;
	ncm_vector_set (cad->inv_lnM->yv, i, lnM1);
	//printf ("prep[%d] % 20.15g % 20.15g % 20.15g % 20.15g\n", i, ntot, f, lnM1, z);
	lnM0 = lnM1;
  }

  for (; i < cad->inv_lnM->len; i++)
  {
	gdouble lnM1 = cad->lnMi + (cad->lnMf - cad->lnMi) / (cad->inv_lnM->len - 1.0) * i;
	gdouble onemn = nc_mass_function_dNdz (cad->mfp, model, lnM1, cad->lnMf, z, use_spline) / dNdz;
	f = _nc_cad_inv_dNdz_convergence_f_onemn (onemn, cad->lnM_epsilon);
	ncm_vector_set (cad->inv_lnM->xv, i, f); //dNdz_u2 / dNdz;
	ncm_vector_set (cad->inv_lnM->yv, i, lnM1);
	//printf ("prep[%d] % 20.15g % 20.15g % 20.15g % 20.15g\n", i, onemn, f, lnM1, z);
  }

  ncm_spline_prepare (cad->inv_lnM);
}

/**
 * nc_cluster_abundance_d2NdzdlnM_spline_val:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @lnM: logarithm base e of mass.
 * @z: redshift.
 *
 * This function computes the value of d2N/dzdlnM at redshift z and logarithm base e of mass lnM.
 *
 * Returns: a gdouble wich represents d2N(lnM, z)/dzdlnM.
 */
static gdouble
_nc_cluster_abundance_d2NdzdlnM_spline_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z)
{
  //printf ("lnM = % 20.5g\n", lnM);
  //return (ncm_spline2d_eval (cad->d2NdzdlnM, lnM, z));
  return nc_mass_function_d2NdzdlnM (cad->mfp, model, lnM, z);
}

/**
 * nc_cluster_abundance_d2NdzdlnM_obs_spline_val:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @lnM: logarithm base e of mass.
 * @z: redshift.
 *
 * This function computes the value of d2N/dzdlnM at redshift z and logarithm base e of mass lnM.
 *
 * Returns: a gdouble wich represents d2N(lnM, z)/dzdlnM.
 */
static gdouble
_nc_cluster_abundance_d2NdzdlnM_obs_spline_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z)
{
  if (ncm_model_ctrl_update (cad->ctrl, NCM_MODEL(model)))
	nc_cluster_abundance_prepare (cad, model);

  return (ncm_spline2d_eval (cad->d2NdzdlnM, lnM, z));
}

/**
 * _nc_cluster_abundance_d2NdzdlnM_Mobs_local_selection_val:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @lnMobs: logarithm base e of observable mass.
 * @z: redshift.
 *
 * This function computes the value of d2N/dzdlnM / p(lnMobs, z) at redshift z and logarithm
 * base e of mass lnMobs, where p(lnMobs, z) is the purity of the respective redshift and mass bins.
 *
 * Returns: a gdouble wich represents d2N(lnMobs, z)/dzdlnM / p(lnMobs, z).
 */
static gdouble
_nc_cluster_abundance_d2NdzdlnM_Mobs_local_selection_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs, gdouble z)
{
  gdouble p, res;

  printf("Estou na purity d2N\n");
  p = _nc_cluster_abundance_one_over_purity (cad, lnMobs, z);
  res = nc_cluster_abundance_d2NdzdlnM_Mobs_local_selection (cad, model, lnMobs, z);
  printf ("p = %.5g res = %.5g\n", p, res);

  return res * p;
}

static gdouble
_nc_cluster_abundance_d2NdzdlnM_photoz_Mobs_local_selection_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs, gdouble zp)
{
  gdouble p, res;

  p = _nc_cluster_abundance_one_over_purity (cad, lnMobs, zp);
  res = nc_cluster_abundance_d2NdzdlnM_photoz_Mobs_local_selection (cad, model, lnMobs, zp);
  //printf ("p = %.5g res = %.5g\n", p, res);

  return res * p;
}

/**
 * nc_cluster_abundance_dNdz_spline_val:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @lnMl: logarithm base e of mass, lower threshold.
 * @lnMu: logarithm base e of mass, upper threshold.
 * @z: redshift.
 *
 * This function computes the value of dN/dz at redshift z and logarithm base e of mass lnM.
 *
 * Returns: a gdouble wich represents dN/dz.
 */
static gdouble
_nc_cluster_abundance_dNdz_spline_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble z)
{
  if (lnMl == cad->lnMi && lnMu == cad->lnMf)
	return (nc_mass_function_dNdz (cad->mfp, model, lnMl, lnMu, z, TRUE));
  else
	return (nc_mass_function_dNdz (cad->mfp, model, lnMl, lnMu, z, FALSE));
}

/**
 * nc_cluster_abundance_dNdz_obs_spline_val:
 * @cad: a pointer to #NcClusterAbundance.
 * @model: a pointer to #NcHICosmo.
 * @lnMl: logarithm base e of mass, lower threshold.
 * @lnMu: logarithm base e of mass, upper threshold.
 * @z: redshift.
 *
 * This function computes the value of dN/dz at redshift z and logarithm base e of mass lnM
 * taking into account photometric redshift and/or mass-observable uncertainties.
 *
 * Returns: a gdouble wich represents dN/dz.
 */
static gdouble
_nc_cluster_abundance_dNdz_obs_spline_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble z)
{
  if (ncm_model_ctrl_update (cad->ctrl, NCM_MODEL(model)))
	nc_cluster_abundance_prepare (cad, model);

  if (lnMl == cad->lnMi && lnMu == cad->lnMf)
	return (ncm_spline2d_integ_dx_spline_val (cad->d2NdzdlnM, lnMl, lnMu, z));
  else
	return (ncm_spline2d_integ_dx (cad->d2NdzdlnM, lnMl, lnMu, z));
}

/**
 * _nc_cluster_abundance_dNdz_Mobs_local_selection_val:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @lnMl: logarithm base e of mass, lower threshold.
 * @lnMu: logarithm base e of mass, upper threshold.
 * @z: redshift.
 *
 * This function computes the value of dN/dz = \int_lnMl^lnMu dlnMobs p(lnMobs, z) * d2N/dzdlnM
 * at redshift z, where p(lnMobs, z) is the purity of the respective redshift and mass bins.
 *
 * Returns: a gdouble wich represents dN/dz taking into account purity and completeness.
 */
static gdouble
_nc_cluster_abundance_dNdz_Mobs_local_selection_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble z)
{
  return nc_cluster_abundance_dNdz_Mobs_local_selection (cad, model, lnMl, lnMu, z);
}

/**
 * _nc_cluster_abundance_dNdz_photoz_Mobs_local_selection_val:
 * @cad: a pointer to #NcClusterAbundance.
 * @model: a pointer to #NcHICosmo.
 * @lnMl: logarithm base e of mass, lower threshold.
 * @lnMu: logarithm base e of mass, upper threshold.
 * @z: redshift.
 *
 * This function computes the value of dN/dz = \int_lnMl^lnMu dlnMobs p(lnMobs, z) * d2N/dzdlnM
 * at redshift z, where p(lnMobs, z) is the purity of the respective redshift and mass bins.
 *
 * Returns: a gdouble wich represents dN/dz taking into account purity and completeness.
 */
static gdouble
_nc_cluster_abundance_dNdz_photoz_Mobs_local_selection_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble z)
{
  return nc_cluster_abundance_dNdz_photoz_Mobs_local_selection (cad, model, lnMl, lnMu, z);
}

static gdouble
_nc_cluster_abundance_N_spline_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble zl, gdouble zu)
{
  if (lnMl == cad->lnMi && lnMu == cad->lnMf)
	return (nc_mass_function_N (cad->mfp, model, lnMl, lnMu, zl, zu, NC_MASS_FUNCTION_SPLINE_LNM));
  else if (zl == cad->zi && zu == cad->zf)
	return (nc_mass_function_N (cad->mfp, model, lnMl, lnMu, zl, zu, NC_MASS_FUNCTION_SPLINE_Z));
  else
	return (nc_mass_function_N (cad->mfp, model, lnMl, lnMu, zl, zu, NC_MASS_FUNCTION_SPLINE_NONE));
}

static gdouble
_nc_cluster_abundance_N_obs_spline_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble zl, gdouble zu)
{
  if (ncm_model_ctrl_update (cad->ctrl, NCM_MODEL(model)))
	nc_cluster_abundance_prepare (cad, model);

  if (lnMl == cad->lnMi && lnMu == cad->lnMf)
	return (ncm_spline2d_integ_dxdy_spline_x (cad->d2NdzdlnM, lnMl, lnMu, zl, zu));
  else
	return (ncm_spline2d_integ_dxdy (cad->d2NdzdlnM, lnMl, lnMu, zl, zu));
}

static gdouble
_nc_cluster_abundance_N_Mobs_local_selection_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble zl, gdouble zu)
{
  return nc_cluster_abundance_N_Mobs_local_selection (cad, model, lnMl, lnMu, zl, zu);
}

static gdouble
_nc_cluster_abundance_N_photoz_Mobs_local_selection_val (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMl, gdouble lnMu, gdouble zl, gdouble zu)
{
  return nc_cluster_abundance_N_photoz_Mobs_local_selection (cad, model, lnMl, lnMu, zl, zu);
}

/* Functions below are with respect to dN/dz */

/**
 * nc_cluster_abundance_bin_realization: (skip)
 * @zr: FIXME
 * @h: FIXME
 *
 * FIXME
 */
void
nc_cluster_abundance_bin_realization (GArray *zr, gsl_histogram **h)
{
  int i;
  for (i = 0; i < zr->len; i++)
  {
	int j;
	const gdouble z = g_array_index (zr, gdouble, i);
	for (j = 0; h[j] != NULL; j++)
	  gsl_histogram_increment (h[j], z);
  }
}

void
nc_cluster_abundance_realizations_save_to_file (GPtrArray *realizations, gchar *filename)
{
  FILE *out;
  gint i, j;
  out = fopen (filename,"w");

  for (j = 0; j < realizations->len; j++)
  {
	GArray *z_real = g_ptr_array_index(realizations, j);
	for (i = 0; i < z_real->len; i++)
	{
	  gdouble zi = g_array_index(z_real, gdouble, i);
	  fprintf (out, "%.15g\n", zi);
	  fflush (out);
	}
	fprintf(out, "\n\n");
  }
  fclose (out);
}

/**
 * nc_cluster_abundance_realizations_read_from_file:
 * @file_realization: (array length=n_realizations): is the file's name which contains the redshift values of the clusters obtained with random Poisson generator.
 * @n_realizations: is the number of realizations generated.
 *
 * GPtrArray *realizations is an array of array with the z values of all realizations. To complete...
 *
 * Returns: (transfer full): FIXME
 */
GPtrArray *
nc_cluster_abundance_realizations_read_from_file (gchar *file_realization, gint n_realizations)
{
  FILE *frealization = fopen (file_realization, "r");
  GPtrArray *realizations = g_ptr_array_sized_new (n_realizations);
  int i, j;
  long int file_position, goby;

  if (frealization == NULL)
  {
	fprintf (stderr, "abundance_random_generator_read_from_file: file %s, do not exist.\n", file_realization);
	exit (0);
  }

  file_position = ftell(frealization);

  for (j = 0; j < n_realizations; j++)
  {
	guint z_total;
	gint n_enter = 0;
	gchar line[5000];
	z_total = 0;    /* Counting the number of redshifts, which corresponds to the total number of clusters, for each realization. */
	while (fgets (line, 5000, frealization) != NULL)
	{
	  if (strlen(line) == 1)
	  {
		n_enter++;
		if (n_enter == 2)
		  break;
		continue;
	  }
	  if (n_enter ==1)
	  {
		fprintf (stderr, "Error, it must not exist just one 'enter'.\n");
		exit(-1);
	  }
	  n_enter = 0;
	  z_total++;
	}

	{
	  GArray *realization_z = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), z_total); /* Array which contains the z values of each realization. */
	  goby = ftell(frealization) - file_position;
	  fseek(frealization, -goby, SEEK_CUR);
	  for (i = 0; i < z_total; i++)
	  {
		gdouble z;
		fscanf(frealization, "%lg\n", &z);
		//        printf ("%lg %u\n", z, z_total);
		g_array_append_val(realization_z, z);
	  }
	  g_ptr_array_add (realizations, realization_z);
	}

	fscanf(frealization, " ");
	fscanf(frealization, " ");
	file_position = ftell(frealization);
  }

  fclose (frealization);
  return realizations;
}

/**
 * nc_cluster_abundance_d2NdzdlnM_photoz:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @lnM: the logarithm base e of the mass (gdouble).
   * @z_phot: the photometric redshift (gdouble).
	 *
 * This function computes /f$ \int_{z_{phot} - 10\sigma_{phot}}^{z_{phot} + 10\sigma_{phot}} dz \,
 * \frac{d^2N}{dzdlnM} * P(z^{photo}|z) /f$. The integral limits were determined requiring a precision
 * to five decimal places.
 *
 * Returns: a gdouble which corresponds to /f$ \int_{z_{phot} - 10\sigma_{phot}}^{z_{phot} + 10\sigma_{phot}} dz \,
 * \frac{d^2N}{dzdlnM} * P(z^{photo}|z) /f$.
 */
gdouble
nc_cluster_abundance_d2NdzdlnM_photoz (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z_phot)
{
  observables_integrand_data photoz_data;
  gdouble d2NdzdlnM;
  gsl_function F;
  photoz_data.cad = cad;
  photoz_data.model = model;

  F.function = &_nc_cluster_abundance_photoz_integrand;
  F.params = &photoz_data;

  {
	gdouble res, err;
	gdouble zl, zu;
	photoz_data.zp = z_phot;
	nc_cluster_photoz_integ_limits (cad->photoz, NULL, z_phot, &zl, &zu); /* FIXME HERE */
	photoz_data.lnMobs = lnM;

	{
	  gsl_integration_workspace **w = nc_integral_get_workspace ();
	  gsl_integration_qag (&F, zl, zu, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
	  ncm_memory_pool_return (w);
	}
	d2NdzdlnM = res;
  }

  return d2NdzdlnM;
}

/**
 * nc_cluster_abundance_d2NdzdlnM_Mobs:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @lnMobs: the logarithm base e of the "observable" mass (gdouble).
 * @z: redshift (gdouble).
 *
 * This function computes /f$ \int_{\ln M^{obs} - 7\sigma_{\ln M}}^{\ln M^{obs} + 7\sigma_{\ln M}} d\ln M \,
 * \frac{d^2N}{dzdlnM} * P(\ln M^{obs}|\ln M) /f$. The integral limits were determined requiring a precision
 * to five decimal places.
 *
 * Returns: a gdouble which corresponds to /f$ \int_{\ln M^{obs} - 7\sigma_{\ln M}}^{\ln M^{obs} + 7\sigma_{\ln M}} d\ln M \,
 * \frac{d^2N}{dzdlnM} * P(\ln M^{obs}|\ln M) /f$.
 */
gdouble
nc_cluster_abundance_d2NdzdlnM_Mobs (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs, gdouble z)
{
  observables_integrand_data md_data;
  gdouble d2NdzdlnM;
  gsl_function F;
  md_data.cad = cad;
  md_data.model = model;

  F.function = &_nc_cluster_abundance_Mobs_integrand;
  F.params = &md_data;

  {
	md_data.zp = z;
	md_data.lnMobs = lnMobs;
	{
	  gdouble res, err;
	  gdouble lnMl, lnMu;
	  lnMl = lnMobs - 7.0 * md_data.cad->lnM_sigma0;
	  lnMu = lnMobs + 7.0 * md_data.cad->lnM_sigma0;
	  {
		gsl_integration_workspace **w = nc_integral_get_workspace ();

		gsl_integration_qag (&F, lnMl, lnMu, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
		ncm_memory_pool_return (w);
	  }
	  d2NdzdlnM = res;
	}
  }
  //printf ("d2N = %.5g", d2NdzdlnM);
  return d2NdzdlnM;
}

/**
 * nc_cluster_abundance_d2NdzdlnM_photoz_Mobs:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @lnMobs: logarithm (base e) of the observable mass.
 * @zp: photometric redshift.
 *
 * This function computes /f$ \int_0^\infty dz \int_0^\infty d\ln M \frac{d^2N(\ln M, z)}{dzd\ln M} * P(z^{phot}|z) *
 * P(\ln M^{obs}|\ln M, z) /f$. We studied the convergence of this integral to optimize this function. We verified
 * that it converges to 5 decimal places at the redshift interval /f$ [z^{phot} - 10\sigma^{phot}, z^{phot} +
 * 10\sigma^{phot}] /f$ and the mass interval /f$ [\ln M^{obs} - 7\sigma_{\ln M}, \ln M^{obs} + 7\sigma_{\ln M}] /f$.
 *
 * Returns: a gdouble which represents /f$ \frac{d^2N(\ln M^{obs}, z^{phot})}{dzd\lnM} /f$.
 */
gdouble
nc_cluster_abundance_d2NdzdlnM_photoz_Mobs (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs, gdouble zp)
{
  gdouble d2NdzdlnM;
  observables_integrand_data obs_data;
  NcIntegrand2dim integ;

  obs_data.cad = cad;
  obs_data.model = model;

  integ.f = _nc_cluster_abundance_photoz_Mobs_integrand;
  integ.userdata = &obs_data;

  {
	gdouble zl, zu;
	obs_data.zp = zp;
	obs_data.lnMobs = lnMobs;
	nc_cluster_photoz_integ_limits (cad->photoz, NULL, zp, &zl, &zu); /* FIXME HERE */
	{
	  gdouble res, err;
	  gdouble lnMl, lnMu;
	  lnMl = lnMobs - 7.0 * obs_data.cad->lnM_sigma0;
	  lnMu = lnMobs + 7.0 * obs_data.cad->lnM_sigma0;
	  ncm_integrate_2dim (&integ, lnMl, zl, lnMu, zu, NC_DEFAULT_PRECISION, 0.0, &res, &err);

	  d2NdzdlnM = res;
	}
  }

  return d2NdzdlnM;
}

/**
 * nc_cluster_abundance_d2NdzdlnM_Mobs_local_selection:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @lnMobs: logarithm (base e) of the observable mass
 * @z: redshift
 *
 * This function computes \int_0^\infty d\ln M \frac{d^2N(\ln M, z)}{dzd\ln M} * c(\lnM, z) *
 * P(\ln M^{obs}|\ln M, z). The sample/data has specific redshift and mass ranges, therefore we have to test if
 * zp and lnMobs are inside these ranges. We have also to do the integral in pieces with limits corresponding to the lnM range of
 * the histogram setting the completeness /f$ c(\ln M, z) /f$ and /f$ \sigma_{\lnM} /f$ spcecific of each bin.
 * We determined the integral limits /f$ [\ln M^{obs} - 7\sigma_{\ln M}, \ln M^{obs} + 7\sigma_{\ln M}] /f$ to obtain
 * a precision to five decimal places.
 *
 * It is worth emphasizing that purity is not included, since it is a function of observable mass and photometric redshift.
 * See #nc_cluster_abundance_d2NdzdlnM_purity_val, #nc_cluster_abundance_dNdz_purity_val and
 * #nc_cluster_abundance_N_purity_val.
 *
 * Returns: a gdouble which represents /f$ \frac{d^2N(\ln M^{obs}, z)}{dzd\lnM} /f$ taking into account the
 * completeness.
 */
gdouble
nc_cluster_abundance_d2NdzdlnM_Mobs_local_selection (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs, gdouble z)
{
  gsize i, j;
  gdouble d2NdzdlnM;
  observables_integrand_data obs_data;
  gsl_function F;

  obs_data.cad = cad;
  obs_data.model = model;
  obs_data.zp = z;

  F.function = &_nc_cluster_abundance_Mobs_local_selection_integrand;
  F.params = &obs_data;

  //printf("entrando d2N z = %.5g\n", z);
  {
	gdouble lnM_max = gsl_histogram2d_ymax (obs_data.cad->completeness);
	gdouble lnM_min = gsl_histogram2d_ymin (obs_data.cad->completeness);

	_check_histogram_limits (obs_data.cad->completeness, lnMobs, z, &i, &j);

	obs_data.lnMobs = lnMobs;
	obs_data.cad->lnM_sigma0 = gsl_histogram2d_get (obs_data.cad->sd_lnM, i, j);
	//obs_data.cad->lnM_sigma0 = gsl_histogram2d_max_val (obs_data.cad->sd_lnM);
	//printf ("sd_lnM max = %.4g\n", obs_data.cad->lnM_sigma0);

	if (obs_data.cad->lnM_sigma0 < 0)
	  return 0.0;

	{
	  gsl_integration_workspace **w = nc_integral_get_workspace ();
	  gdouble res, err;
	  gdouble lnMl, lnMu;
	  gsize i_lnMmin, i_lnMmax, k;
	  gint l, n_integrals;
	  lnMl = obs_data.lnMobs - 3.0 * obs_data.cad->lnM_sigma0;
	  lnMu = obs_data.lnMobs + 3.0 * obs_data.cad->lnM_sigma0;
	  if (lnMl < lnM_min)
		lnMl = lnM_min;
	  if (lnMu > lnM_max)
		lnMu = lnM_max;

      lnMu -= 1e-5;
	  gsl_histogram2d_find (obs_data.cad->sd_lnM, z, lnMl, &k, &i_lnMmin);
	  gsl_histogram2d_find (obs_data.cad->sd_lnM, z, lnMu, &k, &i_lnMmax);
	  n_integrals = i_lnMmax - i_lnMmin + 1;
	  //printf ("lnMl = %g lnMu = %g ymin= %g ymax= %g\n", lnMl, lnMu, obs_data.cad->sd_lnM->yrange[0], obs_data.cad->sd_lnM->yrange[obs_data.cad->sd_lnM->ny]);
	  //printf ("i_max = %zu i_min = %zu n integrals = %d\n", i_lnMmax, i_lnMmin, n_integrals);
	  if (n_integrals == 1)
	  {
		gsl_integration_qag (&F, lnMl, lnMu, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
		d2NdzdlnM = _nc_cluster_abundance_completeness (obs_data.cad, (lnMu + lnMl) / 2.0, obs_data.zp) * res * cad->mfp->area_survey;
		//printf ("Ml = %.5g Mu = %.5g res = %.8g\n", lnMl, lnMu, d2NdzdlnM);
	  }
	  else
	  {
		gdouble lnM_lower, lnM_upper, d2NdzdlnM_over_area;
		gsl_histogram2d_get_yrange (obs_data.cad->sd_lnM, i_lnMmin, &lnM_lower, &lnM_upper);

		gsl_integration_qag (&F, lnMl, lnM_upper, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
		d2NdzdlnM_over_area = res * _nc_cluster_abundance_completeness (obs_data.cad, (lnM_upper + lnMl) / 2.0, obs_data.zp);
		//printf ("bin = %zu Ml = %.5g Mu = %.5g res = %.8g\n", i_lnMmin, lnM_lower, lnM_upper, res);
		for (l = 1; l < (n_integrals - 1); l++)
		{
		  i_lnMmin++;
		  gsl_histogram2d_get_yrange (obs_data.cad->sd_lnM, i_lnMmin, &lnM_lower, &lnM_upper);

		  gsl_integration_qag (&F, lnM_lower, lnM_upper, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
		  d2NdzdlnM_over_area += res * _nc_cluster_abundance_completeness (obs_data.cad, (lnM_lower + lnM_upper) / 2.0, obs_data.zp);
		  //printf ("bin = %zu Ml = %.5g Mu = %.5g res = %.8g d2N = %.8g\n", i_lnMmin, lnM_lower, lnM_upper, res, d2NdzdlnM_over_area);
		}
		gsl_integration_qag (&F, lnM_upper, lnMu, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
		d2NdzdlnM_over_area += res * _nc_cluster_abundance_completeness (obs_data.cad, (lnM_upper + lnMu) / 2.0, obs_data.zp);
		d2NdzdlnM = d2NdzdlnM_over_area * cad->mfp->area_survey;
	  }
	  ncm_memory_pool_return (w);
	}
  }

  //printf("d2NdzdlnM = %.5g\n", d2NdzdlnM);
  return d2NdzdlnM;
}

/**
 * nc_cluster_abundance_dNdz_Mobs_local_selection:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @lnMobs_i: logarithm base e of the observable minimum mass threshold.
 * @lnMobs_f: logarithm base e of the observable maximum mass threshold.
 * @z: redshift.
 *
 * This function computes \int_lnMi^lnMf d lnMobs \frac{d^2N(\ln M, z)}{dzd\ln Mobs} / p(lnMobs, z).
 * The sample/data has specific redshift and mass ranges, therefore we have to test if
 * z and lnMobs are inside these ranges. We have also to do the integral in pieces with limits corresponding to the lnM range of
 * the histogram setting the specific purity /f$ p(\ln Mobs, z) /f$ of each bin.
 *
 * Returns: a gdouble which represents /f$ \frac{dN(\ln M^{obs}, z)}{dz} /f$ taking into account the
 * completeness, purity and P(lnMobs | lnM) with varying standard deviation.
 */
gdouble
nc_cluster_abundance_dNdz_Mobs_local_selection (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs_i, gdouble lnMobs_f, gdouble z)
{
  gdouble dNdz;
  observables_integrand_data obs_data;
  gsl_function F;

  obs_data.cad = cad;
  obs_data.model = model;
  obs_data.zp = z;

  F.function = &_nc_cluster_abundance_dNdz_Mobs_local_selection_integrand;
  F.params = &obs_data;

  {
	gsl_integration_workspace **w = nc_integral_get_workspace ();
	gdouble res, err;
	gsize i_lnMmin, i_lnMmax, k;
	gint l, n_integrals;

	_check_histogram_limits (obs_data.cad->purity, lnMobs_i, z, &i_lnMmin, &k);
	_check_histogram_limits_include_sup (obs_data.cad->purity, lnMobs_f, z, &i_lnMmax, &k);

	n_integrals = i_lnMmax - i_lnMmin + 1;
	//printf ("lnMl = %g lnMu = %g ymin= %g ymax= %g\n", lnMl, lnMu, obs_data.cad->sd_lnM->yrange[0], obs_data.cad->sd_lnM->yrange[obs_data.cad->sd_lnM->ny]);
	//printf ("i_max = %zu i_min = %zu n integrals = %d\n", i_lnMmax, i_lnMmin, n_integrals);
	if (n_integrals == 1)
	{
	  gsl_integration_qag (&F, lnMobs_i, lnMobs_f, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
	  dNdz = _nc_cluster_abundance_one_over_purity (obs_data.cad, (lnMobs_i + lnMobs_f) / 2.0, z) * res;
	  //printf ("Ml = %.5g Mu = %.5g res = %.8g\n", lnMl, lnMu, d2NdzdlnM);
	}
	else
	{
	  gdouble lnM_lower, lnM_upper;
	  gsl_histogram2d_get_yrange (obs_data.cad->purity, i_lnMmin, &lnM_lower, &lnM_upper);

	  gsl_integration_qag (&F, lnMobs_i, lnM_upper, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
	  dNdz = res * _nc_cluster_abundance_one_over_purity (obs_data.cad, (lnM_upper + lnMobs_i) / 2.0, z);
	  //printf ("bin = %zu Ml = %.5g Mu = %.5g res = %.8g\n", i_lnMmin, lnM_lower, lnM_upper, res);
	  for (l = 1; l < (n_integrals - 1); l++)
	  {
		i_lnMmin++;
		gsl_histogram2d_get_yrange (obs_data.cad->purity, i_lnMmin, &lnM_lower, &lnM_upper);

		gsl_integration_qag (&F, lnM_lower, lnM_upper, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
		dNdz += res * _nc_cluster_abundance_one_over_purity (obs_data.cad, (lnM_lower + lnM_upper) / 2.0, z);
		//printf ("bin = %zu Ml = %.5g Mu = %.5g res = %.8g d2N = %.8g\n", i_lnMmin, lnM_lower, lnM_upper, res, dNdz);
	  }
	  gsl_integration_qag (&F, lnM_upper, lnMobs_f, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
	  dNdz += res * _nc_cluster_abundance_one_over_purity (obs_data.cad, (lnM_upper + lnMobs_f) / 2.0, z);

	}
	ncm_memory_pool_return (w);
  }

  //printf("dNdz = %.5g\n", dNdz);
  return dNdz;
}

/**
 * nc_cluster_abundance_N_Mobs_local_selection:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @lnMobs_i: logarithm base e of the observable minimum mass threshold.
 * @lnMobs_f: logarithm base e of the observable maximum mass threshold.
 * @z_i: redshift, minimum threshold of a ith bin.
 * @z_f: redshift, maximum threshold of a ith bin.
 *
 * This function computes the number of clusters with mass between [lnMobs_i, lnMobs_f]
 * in redshift bins, i.e., N_i = /f$ \int_z_i^\z_f dz \int_{lnMobs_i}^{lnMobs_f} d\lnMobs 1/p(lnMobs, z)
 * \int_0^\infty d\ln M \frac{d^2N(\ln M, z)}{dzd\ln M} * c(\lnM, z) * P(\ln M^{obs}|\ln M, z) /f$.
 *
 * Returns: a gdouble which represents the number of clusters per redshift bin taking into account the
 * completeness and purity.
 */
gdouble
nc_cluster_abundance_N_Mobs_local_selection (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs_i, gdouble lnMobs_f, gdouble z_i, gdouble z_f)
{
  gsize i, j;
  gdouble N_per_bin = 0.0;
  observables_integrand_data obs_data;
  NcIntegrand2dim integ;

  obs_data.cad = cad;
  obs_data.model = model;

  integ.f = _nc_cluster_abundance_N_Mobs_local_selection_integrand;
  integ.userdata = &obs_data;

  //printf("uhuuuuu, entrando!\n");
  //printf ("zi = %.5g zf = %.5g\n", z_i, z_f);

  {
	gdouble res, err;
	gsize i_lnMmin, i_lnMmax, i_zmin, i_zmax;
	gsize n_lnM, n_z;

	_check_histogram_limits (obs_data.cad->purity, lnMobs_i, z_i, &i_lnMmin, &i_zmin);
	_check_histogram_limits_include_sup (obs_data.cad->purity, lnMobs_f, z_f, &i_lnMmax, &i_zmax);

	n_lnM = i_lnMmax - i_lnMmin + 2;
	n_z = i_zmax - i_zmin + 2;
	//printf ("lnMl = %g lnMu = %g ymin= %g ymax= %g\n", lnMl, lnMu, obs_data.cad->purity->yrange[0], obs_data.cad->sd_lnM->yrange[obs_data.cad->purity->ny]);
	//printf ("%zu %zu n_z = %zu %zu %zu n_lnM = %zu\n", i_zmax, i_zmin, n_z, i_lnMmax, i_lnMmin, n_lnM);
	if ((n_lnM - 1) * (n_z - 1) == 1)
	{
	  ncm_integrate_2dim (&integ, lnMobs_i, z_i, lnMobs_f, z_f, NC_DEFAULT_PRECISION, 0.0, &res, &err);
	  N_per_bin = _nc_cluster_abundance_one_over_purity (obs_data.cad, (lnMobs_i + lnMobs_f) / 2.0, (z_i + z_f) / 2.0) * res;
	}
	else
	{
	  gdouble z_nodes[1000], lnM_nodes[1000];
	  if (n_z > 2)
		memcpy (&z_nodes[1], &obs_data.cad->purity->xrange[i_zmin + 1], sizeof (gdouble) * (n_z - 2));
	  if (n_lnM > 2)
		memcpy (&lnM_nodes[1], &obs_data.cad->purity->yrange[i_lnMmin + 1], sizeof (gdouble) * (n_lnM - 2));
	  z_nodes[0] = z_i;
	  z_nodes[n_z - 1] = z_f;
	  lnM_nodes[0] = lnMobs_i;
	  lnM_nodes[n_lnM - 1] = lnMobs_f;

	  for (i = 0; i < n_z - 1; i++)
	  {
		for (j = 0; j < n_lnM - 1; j++)
		{
		  //printf ("Vai integrar [% 20.15g % 20.15g] [% 20.15g % 20.15g]\n", lnM_nodes[j], lnM_nodes[j+1], z_nodes[i], z_nodes[i+1]);
		  ncm_integrate_2dim (&integ, lnM_nodes[j], z_nodes[i], lnM_nodes[j+1], z_nodes[i+1], NC_DEFAULT_PRECISION, 0.0, &res, &err);
		  N_per_bin += res * _nc_cluster_abundance_one_over_purity (obs_data.cad, (lnM_nodes[j+1] + lnM_nodes[j]) / 2.0, (z_nodes[i+1] + z_nodes[i]) / 2.0);
		  //printf ("integrou\n");
		}
	  }
	  //printf ("bin = %zu Ml = %.5g Mu = %.5g res = %.8g\n", i_lnMmin, lnM_lower, lnM_upper, res);
	}
	return N_per_bin;
  }
}

/**
 * nc_cluster_abundance_d2NdzdlnM_photoz_Mobs_local_selection:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @lnMobs: the logarithm (base e) of the observable mass (gdouble).
 * @zp: the photometric redshift (gdouble).
 *
 * This function computes /f$ \int_0^\infty dz \int_0^\infty d\ln M \frac{d^2N(\ln M, z)}{dzd\ln M} *
 * c(\lnM, z) * P(\ln M^{obs}|\ln M, z) * P(z^{photo}, z) /f$.
 * Analogously to our previous work, we determine the integration limits considering the
 * maximum value of /f$ \sigma_{\ln M} /f$. Therefore, the mass interval is
 * /f$ [\ln M^{obs} - 7\sigma_{\ln M}, \ln M^{obs} + 7\sigma_{\ln M}] /f$ and the redshift interval is
 * /f$ [z^{photo} - 10\sigma, z^{photo} + 10\sigma] /f$. Purity is included later because it only depends
 * on observables mass and redshift.
 *
 * Returns: a gdouble which represents /f$ \frac{d^2N(\ln M^{obs}, z^{photo})}{dzd\lnM} /f$ taking into account the
 * completeness.
 */
gdouble
nc_cluster_abundance_d2NdzdlnM_photoz_Mobs_local_selection (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs, gdouble zp)
{
  gsize i, j;
  gdouble d2NdzdlnM;
  observables_integrand_data obs_data;
  NcIntegrand2dim integ;

  obs_data.cad = cad;
  obs_data.model = model;

  integ.f = _nc_cluster_abundance_photoz_Mobs_local_integrand;
  integ.userdata = &obs_data;

  {
	gdouble z_max = gsl_histogram2d_xmax (obs_data.cad->completeness);
	gdouble z_min = gsl_histogram2d_xmin (obs_data.cad->completeness);
	gdouble lnM_max = gsl_histogram2d_ymax (obs_data.cad->completeness);
	gdouble lnM_min = gsl_histogram2d_ymin (obs_data.cad->completeness);

	_check_histogram_limits (obs_data.cad->completeness, lnMobs, zp, &i, &j);

	obs_data.zp = zp;
	obs_data.lnMobs = lnMobs;
	obs_data.cad->lnM_sigma0 = gsl_histogram2d_get (obs_data.cad->sd_lnM, i, j);
	if (obs_data.cad->lnM_sigma0 < 0)
	  return 0.0;

	{
	  gdouble res, err;
	  gdouble zl, zu, lnMl, lnMu;
	  gsize i_lnMmin, i_lnMmax, i_zmin, i_zmax;
	  gsize n_lnM, n_z;
	  nc_cluster_photoz_integ_limits (cad->photoz, NULL, zp, &zl, &zu); /* FIXME HERE */
	  if (zl < z_min)
		zl = z_min;
	  if (zu > z_max)
		zu = z_max;
	  lnMl = lnMobs - 7.0 * obs_data.cad->lnM_sigma0;
	  if (lnMl < lnM_min)
		lnMl = lnM_min;
	  lnMu = lnMobs + 7.0 * obs_data.cad->lnM_sigma0;
	  if (lnMu > lnM_max)
		lnMu = lnM_max;

	  _check_histogram_limits (obs_data.cad->sd_lnM, lnMl, zl, &i_lnMmin, &i_zmin);
	  _check_histogram_limits_include_sup (obs_data.cad->sd_lnM, lnMu, zu, &i_lnMmax, &i_zmax);

	  n_lnM = i_lnMmax - i_lnMmin + 2;
	  n_z = i_zmax - i_zmin + 2;
	  //printf ("zl = %.5g zu = %.5g ymin= %g ymax= %g\n", zl, zu, obs_data.cad->sd_lnM->yrange[0], obs_data.cad->sd_lnM->yrange[obs_data.cad->sd_lnM->ny]);
	  //printf ("i_max = %zu i_min = %zu n_lnM = %zu n_z = %zu\n", i_lnMmax, i_lnMmin, n_lnM, n_z);
	  if ((n_lnM - 1) * (n_z - 1) == 1)
	  {
		ncm_integrate_2dim (&integ, lnMl, zl, lnMu, zu, NC_DEFAULT_PRECISION, 0.0, &res, &err);
		d2NdzdlnM = _nc_cluster_abundance_completeness (obs_data.cad, obs_data.lnMobs, obs_data.zp) * res * cad->mfp->area_survey;
	  }
	  else
	  {
		gdouble d2NdzdlnM_over_area = 0.0;
		gdouble z_nodes[1000], lnM_nodes[1000];
		if (n_z > 2)
		  memcpy (&z_nodes[1], &obs_data.cad->sd_lnM->xrange[i_zmin + 1], sizeof (gdouble) * (n_z - 2));
		if (n_lnM > 2)
		  memcpy (&lnM_nodes[1], &obs_data.cad->sd_lnM->yrange[i_lnMmin + 1], sizeof (gdouble) * (n_lnM - 2));
		z_nodes[0] = zl;
		z_nodes[n_z - 1] = zu;
		lnM_nodes[0] = lnMl;
		lnM_nodes[n_lnM - 1] = lnMu;

		for (i = 0; i < n_z - 1; i++)
		{
		  //printf ("integral [%.5g %.5g]\n", z_nodes[i],  z_nodes[i+1]);
		  for (j = 0; j < n_lnM - 1; j++)
		  {
			//printf("Vai entrar! [%.5g %.5g] [%.5g %.5g]\n", z_nodes[i],  z_nodes[i+1], lnM_nodes[j], lnM_nodes[j+1]);
			ncm_integrate_2dim (&integ, lnM_nodes[j], z_nodes[i], lnM_nodes[j+1], z_nodes[i+1], NC_DEFAULT_PRECISION, 0.0, &res, &err);
			d2NdzdlnM_over_area += res * _nc_cluster_abundance_completeness (obs_data.cad, (lnM_nodes[j+1] + lnM_nodes[j]) / 2.0, (z_nodes[i+1] + z_nodes[i]) / 2.0);
			//printf("Vai sair!\n");
		  }
		}
		//printf ("bin = %zu Ml = %.5g Mu = %.5g res = %.8g\n", i_lnMmin, lnM_lower, lnM_upper, res);
		d2NdzdlnM = d2NdzdlnM_over_area * cad->mfp->area_survey;
	  }

	}

	return d2NdzdlnM;
  }
}

/**
 * nc_cluster_abundance_dNdz_photoz_Mobs_local_selection:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcHICosmo.
 * @lnMobs_i: logarithm base e of the observable minimum mass threshold.
 * @lnMobs_f: logarithm base e of the observable maximum mass threshold.
 * @z: redshift.
 *
 * This function computes \int_lnMi^lnMf d lnMobs \frac{d^2N(\ln M, z)}{dzd\ln Mobs} / p(lnMobs, z).
 * The sample/data has specific redshift and mass ranges, therefore we have to test if
 * z and lnMobs are inside these ranges. We have also to do the integral in pieces with limits corresponding to the lnM range of
 * the histogram setting the specific purity /f$ p(\ln Mobs, z) /f$ of each bin.
 *
 * Returns: a gdouble which represents /f$ \frac{dN(\ln M^{obs}, z)}{dz} /f$ taking into account the
 * completeness, purity and P(lnMobs | lnM) with varying standard deviation.
 */
gdouble
nc_cluster_abundance_dNdz_photoz_Mobs_local_selection (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs_i, gdouble lnMobs_f, gdouble z)
{
  gdouble dNdz;
  observables_integrand_data obs_data;
  gsl_function F;

  obs_data.cad = cad;
  obs_data.model = model;

  F.function = &_nc_cluster_abundance_dNdz_photoz_Mobs_local_selection_integrand;
  F.params = &obs_data;

  {
	gsize i_lnMmin, i_lnMmax, k;
	gsl_integration_workspace **w = nc_integral_get_workspace ();
	gdouble res, err;
	gint l, n_integrals;

	_check_histogram_limits (obs_data.cad->purity, lnMobs_i, z, &i_lnMmin, &k);
	_check_histogram_limits_include_sup (obs_data.cad->purity, lnMobs_f, z, &i_lnMmax, &k);

	n_integrals = i_lnMmax - i_lnMmin + 1;
	//printf ("lnMl = %g lnMu = %g ymin= %g ymax= %g\n", lnMl, lnMu, obs_data.cad->sd_lnM->yrange[0], obs_data.cad->sd_lnM->yrange[obs_data.cad->sd_lnM->ny]);
	//printf ("i_max = %zu i_min = %zu n integrals = %d\n", i_lnMmax, i_lnMmin, n_integrals);
	if (n_integrals == 1)
	{
	  gsl_integration_qag (&F, lnMobs_i, lnMobs_f, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
	  dNdz = _nc_cluster_abundance_one_over_purity (obs_data.cad, (lnMobs_i + lnMobs_f) / 2.0, z) * res;
	  //printf ("Ml = %.5g Mu = %.5g res = %.8g\n", lnMl, lnMu, d2NdzdlnM);
	}
	else
	{
	  gdouble lnM_lower, lnM_upper;
	  gsl_histogram2d_get_yrange (obs_data.cad->purity, i_lnMmin, &lnM_lower, &lnM_upper);

	  gsl_integration_qag (&F, lnMobs_i, lnM_upper, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
	  dNdz = res * _nc_cluster_abundance_one_over_purity (obs_data.cad, (lnM_upper + lnMobs_i) / 2.0, z);
	  //printf ("bin = %zu Ml = %.5g Mu = %.5g res = %.8g\n", i_lnMmin, lnM_lower, lnM_upper, res);
	  for (l = 1; l < (n_integrals - 1); l++)
	  {
		i_lnMmin++;
		gsl_histogram2d_get_yrange (obs_data.cad->purity, i_lnMmin, &lnM_lower, &lnM_upper);

		gsl_integration_qag (&F, lnM_lower, lnM_upper, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
		dNdz += res * _nc_cluster_abundance_one_over_purity (obs_data.cad, (lnM_lower + lnM_upper) / 2.0, z);
		//printf ("bin = %zu Ml = %.5g Mu = %.5g res = %.8g d2N = %.8g\n", i_lnMmin, lnM_lower, lnM_upper, res, dNdz);
	  }
	  gsl_integration_qag (&F, lnM_upper, lnMobs_f, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
	  dNdz += res * _nc_cluster_abundance_one_over_purity (obs_data.cad, (lnM_upper + lnMobs_f) / 2.0, z);

	}
	ncm_memory_pool_return (w);

  }

  //printf("dNdz = %.5g\n", dNdz);
  return dNdz;
}

/**
 * nc_cluster_abundance_N_photoz_Mobs_local_selection:
 * @cad: a #NcClusterAbundance.
 * @model: a #NcmModel.
 * @lnMobs_i: logarithm base e of the observable minimum mass threshold.
 * @lnMobs_f: logarithm base e of the observable maximum mass threshold.
 * @zp_i: photometric redshift, minimum threshold of a ith bin.
 * @zp_f: photometric redshift, maximum threshold of a ith bin.
 *
 * This function computes the number of clusters with mass between [lnMobs_i, lnMobs_f]
 * in redshift bins, i.e., N_i = /f$ \int_z_i^\z_f dz \int_{lnMobs_i}^{lnMobs_f} d\lnMobs 1/p(lnMobs, z)
 * \int_0^\infty d\ln M \frac{d^2N(\ln M, z)}{dzd\ln M} * c(\lnM, z) * P(\ln M^{obs}|\ln M, z) *
 * P(zphot | z)/f$.
 *
 * Returns: a gdouble which represents the number of clusters per redshift bin taking into account the
 * completeness and purity.
 */
gdouble
nc_cluster_abundance_N_photoz_Mobs_local_selection (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs_i, gdouble lnMobs_f, gdouble zp_i, gdouble zp_f)
{
  gsize i, j;
  gdouble N_per_bin = 0.0;
  observables_integrand_data obs_data;
  NcIntegrand2dim integ;

  obs_data.cad = cad;
  obs_data.model = model;

  integ.f = _nc_cluster_abundance_N_photoz_Mobs_local_selection_integrand;
  integ.userdata = &obs_data;

  {
	gdouble res, err;
	gsize i_lnMmin, i_lnMmax, i_zmin, i_zmax;
	gsize n_lnM, n_z;

	_check_histogram_limits (obs_data.cad->purity, lnMobs_i, zp_i, &i_lnMmin, &i_zmin);
	_check_histogram_limits_include_sup (obs_data.cad->purity, lnMobs_f, zp_f, &i_lnMmax, &i_zmax);

    n_lnM = i_lnMmax - i_lnMmin + 2;
	n_z = i_zmax - i_zmin + 2;
	//printf ("lnMl = %g lnMu = %g ymin= %g ymax= %g\n", lnMl, lnMu, obs_data.cad->sd_lnM->yrange[0], obs_data.cad->sd_lnM->yrange[obs_data.cad->sd_lnM->ny]);
	//printf ("i_max = %zu i_min = %zu nM = %zu nz = %zu\n", i_lnMmax, i_lnMmin, n_lnM, n_z);
	if ((n_lnM - 1) * (n_z - 1) == 1)
	{
	  ncm_integrate_2dim (&integ, lnMobs_i, zp_i, lnMobs_f, zp_f, NC_DEFAULT_PRECISION, 0.0, &res, &err);
	  N_per_bin = _nc_cluster_abundance_one_over_purity (obs_data.cad, (lnMobs_i + lnMobs_f) / 2.0, (zp_i + zp_f) / 2.0) * res;
	}
	else
	{
	  gdouble z_nodes[1000], lnM_nodes[1000];
	  if (n_z > 2)
		memcpy (&z_nodes[1], &obs_data.cad->purity->xrange[i_zmin + 1], sizeof (gdouble) * (n_z - 2));
	  if (n_lnM > 2)
		memcpy (&lnM_nodes[1], &obs_data.cad->purity->yrange[i_lnMmin + 1], sizeof (gdouble) * (n_lnM - 2));
	  z_nodes[0] = zp_i;
	  z_nodes[n_z - 1] = zp_f;
	  lnM_nodes[0] = lnMobs_i;
	  lnM_nodes[n_lnM - 1] = lnMobs_f;

	  for (i = 0; i < n_z - 1; i++)
	  {
		for (j = 0; j < n_lnM - 1; j++)
		{
		  //printf ("[%.5g %.5g] [%.5g %.5g]\n", z_nodes[i], z_nodes[i+1], lnM_nodes[j], lnM_nodes[j+1]);
		  ncm_integrate_2dim (&integ, lnM_nodes[j], z_nodes[i], lnM_nodes[j+1], z_nodes[i+1], NC_DEFAULT_PRECISION, 0.0, &res, &err);
		  N_per_bin += res * _nc_cluster_abundance_one_over_purity (obs_data.cad, (lnM_nodes[j+1] + lnM_nodes[j]) / 2.0, (z_nodes[i+1] + z_nodes[i]) / 2.0);
		}
	  }
	  //printf ("bin = %zu Ml = %.5g Mu = %.5g res = %.8g\n", i_lnMmin, lnM_lower, lnM_upper, res);

	}

	return N_per_bin;
	//	return 0.0;
  }
}

/********* Function to compute the mean bias ************************/

static gdouble
_nc_ca_mean_bias_numerator_integrand (gdouble lnM, gpointer params)
{
  observables_integrand_data *md_data = (observables_integrand_data *) params;
  NcClusterAbundance *cad = md_data->cad;

  /* In this case zp is the true redshift, i.e., without uncertainty. */
  gdouble dbdlnM = nc_halo_bias_func_integrand (cad->mbiasf, md_data->model, lnM, md_data->zp);

  //printf ("% 20.8e % 20.8g % 20.15g\n", md_data->lnMobs, md_data->zp, dbdlnM);
  return dbdlnM;
}

gdouble
nc_ca_mean_bias_numerator (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z)
{
  observables_integrand_data md_data;
  gdouble mean_bias_numerator;
  gsl_function F;
  md_data.cad = cad;
  md_data.model = model;

  F.function = &_nc_ca_mean_bias_numerator_integrand;
  F.params = &md_data;

  {
	gdouble res, err;
	gdouble lnMf = cad->lnMf;
	md_data.zp = z;
	md_data.lnMobs = lnM;
	//printf ("%5.5e, %5.5e\n", exp(lnM), exp(lnMf));

	{
	  gsl_integration_workspace **w = nc_integral_get_workspace ();
	  gsl_integration_qag (&F, lnM, lnMf, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
	  ncm_memory_pool_return (w);
	}
	mean_bias_numerator = res;
  }

  return mean_bias_numerator;
}

static gdouble
_nc_ca_mean_bias_denominator_integrand (gdouble lnM, gpointer params)
{
  observables_integrand_data *md_data = (observables_integrand_data *) params;
  NcClusterAbundance *cad = md_data->cad;

  /* In this case zp is the true redshift, i.e., without uncertainty. */
  gdouble n = nc_mass_function (cad->mfp, md_data->model, lnM, md_data->zp);

  //printf ("% 20.8e % 20.8g % 20.15g\n", md_data->lnMobs, md_data->zp, dbdlnM);
  return n;
}

gdouble
nc_ca_mean_bias_denominator (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z)
{
  observables_integrand_data md_data;
  gdouble mean_bias_denominator;
  gsl_function F;
  md_data.cad = cad;
  md_data.model = model;

  F.function = &_nc_ca_mean_bias_denominator_integrand;
  F.params = &md_data;

  {
	gdouble res, err;
	gdouble lnMf = cad->lnMf;
	md_data.zp = z;
	md_data.lnMobs = lnM;
	//printf ("%5.5e, %5.5e\n", exp(lnM), exp(lnMf));

	{
	  gsl_integration_workspace **w = nc_integral_get_workspace ();
	  gsl_integration_qag (&F, lnM, lnMf, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
	  ncm_memory_pool_return (w);
	}
	mean_bias_denominator= res;
  }

  return mean_bias_denominator;
}

gdouble
nc_ca_mean_bias (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnM, gdouble z)
{
  //gdouble M1 = exp(lnM);
  gdouble numerator = nc_ca_mean_bias_numerator (cad, model, lnM, z);
  gdouble denominator = nc_ca_mean_bias_denominator (cad, model, lnM, z);
  gdouble mean_bias = numerator / denominator;

  //printf ("M1 = %5.5e num = %10.5g den = %10.5g\n", M1, numerator, denominator);
  return mean_bias;
}

static gdouble
_nc_ca_mean_bias_Mobs_numerator_integrand (gdouble lnMobs, gpointer params)
{
  observables_integrand_data *md_data = (observables_integrand_data *) params;
  NcClusterAbundance *cad = md_data->cad;

  /* In this case zp is the true redshift, i.e., without uncertainty. */
  gdouble p_M_Mobs = _nc_cluster_abundance_lognormal_mass_dist (md_data, lnMobs);
  gdouble dbdlnM = nc_halo_bias_func_integrand (cad->mbiasf, md_data->model, lnMobs, md_data->zp);

  //printf ("% 20.8e % 20.8g % 20.15g % 20.15g % 20.15g\n", md_data->lnMobs, md_data->zp, p_M_Mobs, dbdlnM, p_M_Mobs * dbdlnM);
  return dbdlnM * p_M_Mobs;
}

gdouble
nc_ca_mean_bias_Mobs_numerator (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs, gdouble z)
{
  observables_integrand_data md_data;
  gdouble mean_bias_Mobs_numerator;
  gsl_function F;
  md_data.cad = cad;
  md_data.model = model;

  F.function = &_nc_ca_mean_bias_Mobs_numerator_integrand;
  F.params = &md_data;

  {
	md_data.zp = z;
	md_data.lnMobs = lnMobs;
	{
	  gdouble res, err;
	  gdouble lnMl, lnMu;
	  lnMl = GSL_MAX (lnMobs - 7.0 * md_data.cad->lnM_sigma0, LNM_MIN);
	  //lnMu = GSL_MIN (lnMobs + 7.0 * md_data.cad->lnM_sigma0, 16.3 * M_LN10);
	  //lnMl = lnMobs - 7.0 * md_data.cad->lnM_sigma0;
	  lnMu = lnMobs + 7.0 * md_data.cad->lnM_sigma0;
	  {
		gsl_integration_workspace **w = nc_integral_get_workspace ();

		gsl_integration_qag (&F, lnMl, lnMu, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
		ncm_memory_pool_return (w);
	  }
	  mean_bias_Mobs_numerator = res;
	}
  }

  return mean_bias_Mobs_numerator;
}

static gdouble
_nc_ca_mean_bias_Mobs_denominator_integrand (gdouble lnMobs, gpointer params)
{
  observables_integrand_data *md_data = (observables_integrand_data *) params;
  NcClusterAbundance *cad = md_data->cad;

  /* In this case zp is the true redshift, i.e., without uncertainty. */
  gdouble p_M_Mobs = _nc_cluster_abundance_lognormal_mass_dist (md_data, lnMobs);
  gdouble n = nc_mass_function (cad->mfp, md_data->model, lnMobs, md_data->zp);

  //printf ("% 20.8e % 20.8g % 20.15g\n", md_data->lnMobs, md_data->zp, n * p_M_Mobs);
  return n * p_M_Mobs;
}

gdouble
nc_ca_mean_bias_Mobs_denominator (NcClusterAbundance *cad, NcHICosmo *model, gdouble lnMobs, gdouble z)
{
  observables_integrand_data md_data;
  gdouble mean_bias_Mobs_denominator;
  gsl_function F;
  md_data.cad = cad;
  md_data.model = model;

  F.function = &_nc_ca_mean_bias_Mobs_denominator_integrand;
  F.params = &md_data;

  {
	md_data.zp = z;
	md_data.lnMobs = lnMobs;
	{
	  gdouble res, err;
	  gdouble lnMl, lnMu;
	  lnMl = GSL_MAX (lnMobs - 7.0 * md_data.cad->lnM_sigma0, LNM_MIN);
	  //lnMu = GSL_MIN (lnMobs + 7.0 * md_data.cad->lnM_sigma0, 16.3 * M_LN10);
	  //lnMl = lnMobs - 7.0 * md_data.cad->lnM_sigma0;
	  lnMu = lnMobs + 7.0 * md_data.cad->lnM_sigma0;
	  {
		gsl_integration_workspace **w = nc_integral_get_workspace ();

		gsl_integration_qag (&F, lnMl, lnMu, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, *w, &res, &err);
		ncm_memory_pool_return (w);
	  }
	  mean_bias_Mobs_denominator = res;
	}
  }

  return mean_bias_Mobs_denominator;
}
