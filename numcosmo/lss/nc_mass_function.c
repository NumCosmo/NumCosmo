/***************************************************************************
 *            nc_mass_function.c
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
 * SECTION:nc_mass_function
 * @title: Mass Function
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glib.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_const_mksa.h>
#include <glib.h>

enum
{
  PROP_0,
  PROP_DISTANCE,
  PROP_MATTER_VAR,
  PROP_GROWTH,
  PROP_MULTIPLICITY
};

G_DEFINE_TYPE (NcMassFunction, nc_mass_function, G_TYPE_OBJECT);

/**
 * nc_mass_function_new:
 * @dist: a #NcDistance sets to #NcMassFunction:distance.
 * @vp: a #NcMatterVar sets to #NcMassFunction:variance.
 * @gf: a #NcGrowthFunc sets to #NcMassFunction:growth.
 * @mulf: a #NcMultiplicityFunc sets to #NcMassFunction:multiplicity.
 *
 * This function allocates memory for a new #NcMassFunction object and sets its properties to the values from
 * the input arguments.
 *
 * Returns: A new #NcMassFunction.
*/
NcMassFunction *
nc_mass_function_new (NcDistance *dist, NcMatterVar *vp, NcGrowthFunc *gf, NcMultiplicityFunc *mulf)
{
  NcMassFunction *mfp = g_object_new (NC_TYPE_MASS_FUNCTION,
				      "distance", dist,
				      "variance", vp,
				      "growth", gf,
				      "multiplicity", mulf,
				      NULL);
  return mfp;
}

/**
 * nc_mass_function_copy:
 * @mfp: a #NcMassFunction.
 *
 * This function duplicates the #NcMassFunction object setting the same values of the original propertities.
 *
 * Returns: (transfer full): A new #NcMassFunction.
*/
NcMassFunction *
nc_mass_function_copy (NcMassFunction * mfp)
{
  return nc_mass_function_new (mfp->dist, mfp->vp, mfp->gf, mfp->mulf);
}

/**
 * nc_mass_function_free:
 * @mfp: a #NcMassFunction.
 *
 * Atomically decrements the reference count of @mfp by one. If the reference count drops to 0,
 * all memory allocated by @mfp is released.
 *
*/
void
nc_mass_function_free (NcMassFunction * mfp)
{
  g_clear_object (&mfp);
}

static void
nc_mass_function_init (NcMassFunction * mfp)
{
  mfp->lnMi = 0.0;
  mfp->lnMf = 0.0;
  mfp->zi = 0.0;
  mfp->zf = 0.0;
  mfp->area_survey = 0.0;
  mfp->N_sigma = 0.0;
  mfp->growth = 0.0;
  mfp->d2NdzdlnM = NULL;
  mfp->ctrl = ncm_model_ctrl_new (NULL);
}

static void
_nc_mass_function_dispose (GObject * object)
{
  /* TODO: Add deinitalization code here */
  NcMassFunction *mfp = NC_MASS_FUNCTION (object);
  nc_distance_free (mfp->dist);
  nc_matter_var_free (mfp->vp);
  nc_growth_func_free (mfp->gf);
  nc_multiplicity_func_free (mfp->mulf);
  ncm_model_ctrl_free (mfp->ctrl);
  if (mfp->d2NdzdlnM != NULL)
	ncm_spline2d_free (mfp->d2NdzdlnM);
  G_OBJECT_CLASS (nc_mass_function_parent_class)->dispose (object);
}

static void
_nc_mass_function_finalize (GObject * object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_mass_function_parent_class)->finalize (object);
}

static void
_nc_mass_function_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcMassFunction *mfp = NC_MASS_FUNCTION (object);
  g_return_if_fail (NC_IS_MASS_FUNCTION (object));

  switch (prop_id)
  {
    case PROP_DISTANCE:
      mfp->dist = g_value_get_object (value);
      break;
    case PROP_MATTER_VAR:
      mfp->vp = g_value_get_object (value);
      break;
    case PROP_GROWTH:
      mfp->gf = g_value_get_object (value);
      break;
    case PROP_MULTIPLICITY:
      mfp->mulf = g_value_get_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_mass_function_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMassFunction *mfp = NC_MASS_FUNCTION (object);
  g_return_if_fail (NC_IS_MASS_FUNCTION (object));

  switch (prop_id)
  {
    case PROP_DISTANCE:
      g_value_set_object (value, mfp->dist);
      break;
    case PROP_MATTER_VAR:
      g_value_set_object (value, mfp->vp);
      break;
    case PROP_GROWTH:
      g_value_set_object (value, mfp->gf);
      break;
    case PROP_MULTIPLICITY:
      g_value_set_object (value, mfp->mulf);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_mass_function_class_init (NcMassFunctionClass * klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  //GObjectClass* parent_class = G_OBJECT_CLASS (klass);

  object_class->dispose = _nc_mass_function_dispose;
  object_class->finalize = _nc_mass_function_finalize;
  object_class->set_property = _nc_mass_function_set_property;
  object_class->get_property = _nc_mass_function_get_property;

  /**
   * NcMassFunction:distance:
   *
   * This property keeps the distance object.
   */
  g_object_class_install_property (object_class,
				   PROP_DISTANCE,
				   g_param_spec_object ("distance",
							NULL,
							"Distance.",
							NC_TYPE_DISTANCE,
							G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME
							| G_PARAM_STATIC_BLURB));

  /**
   * NcMassFunction:variance:
   *
   * This property keeps the matter variance object.
   */
  g_object_class_install_property (object_class,
				   PROP_MATTER_VAR,
				   g_param_spec_object ("variance",
						       NULL,
						       "Variance.",
						       NC_TYPE_MATTER_VAR,
						       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME
						       | G_PARAM_STATIC_BLURB));
  /**
   * NcMassFunction:growth:
   *
   * This property keeps the growth function object.
   */
  g_object_class_install_property (object_class,
				   PROP_GROWTH,
				   g_param_spec_object ("growth",
						       NULL,
						       "Growth function.",
						       NC_TYPE_GROWTH_FUNC,
						       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME
						       | G_PARAM_STATIC_BLURB));
	/**
   * NcMassFunction:multiplicity:
   *
   * This property keeps the multiplicity function object.
   */
  g_object_class_install_property (object_class,
				   PROP_MULTIPLICITY,
				   g_param_spec_object ("multiplicity",
						       NULL,
						       "Multiplicity function.",
						       NC_TYPE_MULTIPLICITY_FUNC,
						       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME
						       | G_PARAM_STATIC_BLURB));

}

/**
 * nc_mass_function:
 * @mfp: a #NcMassFunction.
 * @model: a #NcHICosmo.
 * @lnM: logarithm base e of mass.
 * @z: redshift.
 *
 * This function computes the comoving number density of dark matter halos at redshift z and
 * mass M.
 *
 * Returns: dn/dlnM.
 */
gdouble
nc_mass_function (NcMassFunction *mfp, NcHICosmo *model, gdouble lnM, gdouble z)
{
  gdouble lnR = nc_matter_var_lnM_to_lnR (mfp->vp, model, lnM);
  gdouble R = exp (lnR);
  gdouble M_rho = nc_window_volume (mfp->vp->wp) * gsl_pow_3 (R);	/* M/\rho comoving, rho independe de z */
  gdouble var0 = nc_matter_var_var0 (mfp->vp, model, lnR);
  gdouble dlnvar0_dlnR = nc_matter_var_dlnvar0_dlnR (mfp->vp, model, lnR);
  gdouble growth = nc_growth_func_eval (mfp->gf, model, z);
  gdouble sigma = nc_matter_var_sigma8_sqrtvar0 (mfp->vp, model) * sqrt (var0) * growth;
  //gdouble sigma = sqrt(var0) * growth;
  gdouble f = nc_multiplicity_func_eval (mfp->mulf, model, sigma, z);
  gdouble dn_dlnM = -(1.0 / (3.0 * M_rho)) * f * 0.5 * dlnvar0_dlnR;	/* dn/dlnM = - (\rho/3M) * f * (R/\sigma)* dsigma_dR */
  //printf ("multip % .5g\n", f);

  return dn_dlnM;
}

/* Com a funcao abaixo calculamos \sigma somente uma vez e podemos acessar os valores de \sigma e dn/dM para fazer o grafico */
/* coloca os valores de dn/dm e sigma nos enderecos apontados pelos ponteiros */
/**
 * nc_mass_function_sigma:
 * @mfp: a #NcMassFunction.
 * @model: a #NcHICosmo.
 * @lnM: logarithm base e of mass.
 * @z: redshift.
 * @dn_dlnM_ptr: pointer to comoving number density of halos.
 * @sigma_ptr: pointer to the standard deviation of the density contrast.
 *
 * This function computes the standard deviation of density contrast of the matter fluctuations and
 * the the comoving number density of dark matter halos at redshift z and mass M.
 * These values are stored in @sigma_ptr and @dn_dlnM_ptr, respectively.
 *
*/
void
nc_mass_function_sigma (NcMassFunction * mfp, NcHICosmo * model, gdouble lnM, gdouble z, gdouble * dn_dlnM_ptr,
			gdouble * sigma_ptr)
{
  const gdouble lnR = nc_matter_var_lnM_to_lnR (mfp->vp, model, lnM);
  const gdouble M_rho = nc_window_volume (mfp->vp->wp) * gsl_pow_3 (exp (lnR));	/* M/\rho comoving, rho independe de z */
  const gdouble growth = nc_growth_func_eval (mfp->gf, model, z);
  const gdouble sigma8_sqrtvar0 = nc_matter_var_sigma8_sqrtvar0 (mfp->vp, model);
  const gdouble var0 = nc_matter_var_var0 (mfp->vp, model, lnR);
  const gdouble dlnvar0_dlnR = nc_matter_var_dlnvar0_dlnR (mfp->vp, model, lnR);
  const gdouble sigma = sigma8_sqrtvar0 * sqrt (var0) * growth;
  const gdouble f = nc_multiplicity_func_eval (mfp->mulf, model, sigma, z);

  *sigma_ptr = sigma;
  *dn_dlnM_ptr = -(1.0 / (3.0 * M_rho)) * f * (1.0 / 2.0) * dlnvar0_dlnR;	/* dn/dlnM = - (\rho/3M) * f * (R/\sigma)* dsigma_dR */

  return;
}

/* A funcao abaixo foi feita para compararmos o nosso calculo com o resultado (equacao 9) do artigo astro-ph/0110246.*/
/**
 * nc_mass_function_alpha_eff:
 * @vp: a #NcMatterVar.
 * @model: a #NcHICosmo.
 * @lnM: logarithm base e of mass.
 * @a_eff_ptr: FIXME
 *
 * FIXME
 *
 */
void
nc_mass_function_alpha_eff (NcMatterVar *vp, NcHICosmo * model, gdouble lnM, gdouble * a_eff_ptr)
{
  gdouble lnR = nc_matter_var_lnM_to_lnR (vp, model, lnM);
  gdouble dlnvar0_dlnR = nc_matter_var_dlnvar0_dlnR (vp, model, lnR);
  *a_eff_ptr = -(1.0 / 3.0) * (1.0 / 2.0) * dlnvar0_dlnR;	/* (M/rho)*(1/f)*dn/dlnM = - (1/3) * (1/\sigma)* dsigma_dlnR */

  return;
}

typedef struct _nc_ca_integ
{
  NcMassFunction *mfp;
  NcHICosmo *model;
  gdouble z;
} nc_ca_integ;

/**
 * nc_mass_function_cluster_abundance_integrand:
 * @R: FIXME
 * @params: FIXME
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_mass_function_cluster_abundance_integrand (gdouble R, gpointer params)
{
  nc_ca_integ *ca_integ = (nc_ca_integ *) params;
  NcHICosmo *model = ca_integ->model;
  NcMassFunction *mfp = ca_integ->mfp;
  const gdouble lnR = log (R);
  {
    const gdouble M_rho = nc_window_volume (mfp->vp->wp) * gsl_pow_3 (R);
    const gdouble var0 = nc_matter_var_var0 (mfp->vp, model, lnR);
    const gdouble dlnvar0_dR = nc_matter_var_dlnvar0_dR (mfp->vp, model, lnR);
    gdouble sigma, f, d2cadvdR;

    sigma = mfp->N_sigma * sqrt (var0) * mfp->growth;
    f = nc_multiplicity_func_eval (mfp->mulf, model, sigma, ca_integ->z);
    d2cadvdR = -(1.0 / M_rho) * f * (1.0 / 2.0) * dlnvar0_dR;

    return d2cadvdR;
  }
}

/**
 * nc_mass_function_dcluster_abundance_dv:
 * @mfp: a #NcMassFunction.
 * @model: a #NcHICosmo.
 * @M: mass in units of h^{-1} * M_sun.
 * @z: redshift.
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_mass_function_dn_M_to_inf_dv (NcMassFunction * mfp, NcHICosmo * model, gdouble M, gdouble z)	/* integracao em M. */
{
  static gsl_integration_workspace *w = NULL;
  gsl_function F;
  gdouble n, n_part, error, R;
  nc_ca_integ ca_integ;

  ca_integ.mfp = mfp;
  ca_integ.model = model;
  ca_integ.z = z;
  F.function = &nc_mass_function_cluster_abundance_integrand;
  F.params = &ca_integ;

  if (w == NULL)
    w = gsl_integration_workspace_alloc (NC_INT_PARTITION);

  n = 0.0;
  R = nc_matter_var_mass_to_R (mfp->vp, model, M);
  mfp->growth = nc_growth_func_eval (mfp->gf, model, z);

#define _NC_STEP 20.0

  while (1)
  {
    gsl_integration_qag (&F, R, R + _NC_STEP, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, w, &n_part, &error);
    R += _NC_STEP;
    n += n_part;

    if (gsl_fcmp (n, n + n_part, NC_DEFAULT_PRECISION * 1e-1) == 0)
      break;
  }

  return n;
}

/**
 * nc_mass_function_dn_M1_to_M2_dv:
 * @mfp: a #NcMassFunction.
 * @model: a #NcHICosmo.
 * @M1: mass in units of h^{-1} M_sun - lower limit of integration.
 * @M2: mass in units of h^{-1} M_sun- upper limit of integration.
 * @z: redshift.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_mass_function_dn_M1_to_M2_dv (NcMassFunction * mfp, NcHICosmo * model, gdouble M1, gdouble M2, gdouble z)	/* integracao em M. */
{
  static gsl_integration_workspace *w = NULL;
  gsl_function F;
  gdouble n, error, R1, R2;
  nc_ca_integ ca_integ;

  ca_integ.mfp = mfp;
  ca_integ.model = model;
  ca_integ.z = z;
  F.function = &nc_mass_function_cluster_abundance_integrand;
  F.params = &ca_integ;

  if (w == NULL)
    w = gsl_integration_workspace_alloc (NC_INT_PARTITION);

  n = 0.0;
  R1 = nc_matter_var_mass_to_R (mfp->vp, model, M1);
  R2 = nc_matter_var_mass_to_R (mfp->vp, model, M2);
  mfp->growth = nc_growth_func_eval (mfp->gf, model, z);

  gsl_integration_qag (&F, R1, R2, 0.0, NC_DEFAULT_PRECISION, NC_INT_PARTITION, 1, w, &n, &error);

  printf ("R1 = %5.5g R2 = %5.5g z = %5.5g n = %5.5g\n", R1, R2, z, n);
  return n;
}

/**
 * nc_mass_function_dcomoving_volume_dzdomega:
 * @mfp: a #NcMassFunction.
 * @model: a #NcHICosmo.
 * @z: redshift.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_mass_function_dcomoving_volume_dzdomega (NcMassFunction * mfp, NcHICosmo * model, gdouble z)
{
  const gdouble E = sqrt (nc_hicosmo_E2 (model, z));
  gdouble dc = NC_C_HUBBLE_RADIUS * nc_distance_comoving (mfp->dist, model, z);
  gdouble dV_dzdOmega = gsl_pow_2 (dc) * NC_C_HUBBLE_RADIUS / E;
  return dV_dzdOmega;
}

static void _nc_mass_function_generate_2Dspline_knots (NcMassFunction * mfp, NcHICosmo * model, gdouble rel_error);

/**
 * nc_mass_function_d2NdzdlnM_optimize:
 * @mfp: a #NcMassFunction.
 * @model: a #NcHICosmo.
 * @lnMi: minimum logarithm base e of mass.
 * @lnMf: maximum logarithm base e of mass.
 * @zi: minimum redshift.
 * @zf: maximum redshift.
 *
 * FIXME
 *
 */
void
nc_mass_function_d2NdzdlnM_optimize (NcMassFunction *mfp, NcHICosmo *model, gdouble lnMi, gdouble lnMf, gdouble zi, gdouble zf)
{
  if (lnMi != mfp->lnMi || mfp->lnMf != lnMf || mfp->zi != zi || mfp->zf != zf)
  {
    mfp->lnMi = lnMi;
    mfp->lnMf = lnMf;
    mfp->zi = zi;
    mfp->zf = zf;
    if (mfp->d2NdzdlnM != NULL)
    {
      ncm_spline2d_free (mfp->d2NdzdlnM);
      mfp->d2NdzdlnM = NULL;
    }
  }

  if (mfp->d2NdzdlnM == NULL)
    _nc_mass_function_generate_2Dspline_knots (mfp, model, 1e-5);
}

/**
 * nc_mass_function_d2NdzdlnM_prepare:
 * @mfp: a #NcMassFunction.
 * @model: a #NcHICosmo.
 *
 * FIXME
 *
 */
void
nc_mass_function_d2NdzdlnM_prepare (NcMassFunction *mfp, NcHICosmo *model)
{
  gint i, j;
  if (mfp->d2NdzdlnM == NULL)
  {
    g_error
      ("nc_mass_function_d2NdzdlnM_prepare: called without a previous call of nc_mass_function_d2NdzdlnM_optimize");
    return;
  }
#define D2NDZDLNM_Z(cad) ((cad)->d2NdzdlnM->yv)
#define D2NDZDLNM_LNM(cad) ((cad)->d2NdzdlnM->xv)
#define D2NDZDLNM_VAL(cad) ((cad)->d2NdzdlnM->zm)

  for (i = 0; i < ncm_vector_len (D2NDZDLNM_Z (mfp)); i++)
  {
    const gdouble z = ncm_vector_get (D2NDZDLNM_Z (mfp), i);
    const gdouble dVdz = mfp->area_survey * nc_mass_function_dcomoving_volume_dzdomega (mfp, model, z);

    for (j = 0; j < ncm_vector_len (D2NDZDLNM_LNM (mfp)); j++)
    {
      const gdouble lnM = ncm_vector_get (D2NDZDLNM_LNM (mfp), j);
      const gdouble d2NdzdlnM_ij = (dVdz * nc_mass_function (mfp, model, lnM, z));
      ncm_matrix_set (D2NDZDLNM_VAL (mfp), i, j, d2NdzdlnM_ij);
    }
  }
  ncm_spline2d_prepare (mfp->d2NdzdlnM);

  ncm_model_ctrl_update (mfp->ctrl, NCM_MODEL(model));
}

/**
 * nc_mass_function_d2NdzdlnM:
 * @mfp: a #NcMassFunction.
 * @model: a #NcHICosmo.
 * @lnM: logarithm base e of mass.
 * @z: redshift.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_mass_function_d2NdzdlnM (NcMassFunction * mfp, NcHICosmo * model, gdouble lnM, gdouble z)
{
  const gboolean has_d2NdzdlnM_spline = (mfp->d2NdzdlnM != NULL);
  const gboolean in_spline = has_d2NdzdlnM_spline && (lnM >= mfp->lnMi) && (lnM <= mfp->lnMf) && (z >= mfp->zi)
    && (z <= mfp->zf);
  gdouble d2N_dzdlnM;
  if (in_spline)
  {
    if (ncm_model_ctrl_update (mfp->ctrl, NCM_MODEL(model)))
      nc_mass_function_d2NdzdlnM_prepare (mfp, model);
    d2N_dzdlnM = ncm_spline2d_eval (mfp->d2NdzdlnM, lnM, z);
  }
  else
  {
	d2N_dzdlnM = mfp->area_survey *
	  nc_mass_function (mfp, model, lnM, z) *
	  nc_mass_function_dcomoving_volume_dzdomega (mfp, model, z);
  }
  return d2N_dzdlnM;
}

/**
 * nc_mass_function_dNdz:
 * @mfp: a #NcMassFunction.
 * @model: a #NcHICosmo.
 * @lnMl: logarithm base e of mass, lower threshold.
 * @lnMu: logarithm base e of mass, upper threshold.
 * @z: redshift.
 * @spline: Whenever to create an intermediary spline of the mass integration.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_mass_function_dNdz (NcMassFunction * mfp, NcHICosmo * model, gdouble lnMl, gdouble lnMu, gdouble z, gboolean spline)
{
  const gboolean has_d2NdzdlnM_spline = (mfp->d2NdzdlnM != NULL);
  const gboolean in_spline = has_d2NdzdlnM_spline && (lnMl >= mfp->lnMi) && (lnMu <= mfp->lnMf) && (z >= mfp->zi)
    && (z <= mfp->zf);
  gdouble dN_dz;
  if (in_spline)
  {
    if (ncm_model_ctrl_update (mfp->ctrl, NCM_MODEL(model)))
      nc_mass_function_d2NdzdlnM_prepare (mfp, model);
    if (spline)
      dN_dz = ncm_spline2d_integ_dx_spline_val (mfp->d2NdzdlnM, lnMl, lnMu, z);
    else
      dN_dz = ncm_spline2d_integ_dx (mfp->d2NdzdlnM, lnMl, lnMu, z);
  }
  else
    g_error ("This option is not implemented. It must be used the optimize option.\n");
  return dN_dz;
}

/**
 * nc_mass_function_N:
 * @mfp: a #NcMassFunction.
 * @model: a #NcHICosmo.
 * @lnMl: logarithm base e of mass, lower threshold.
 * @lnMu: logarithm base e of mass, upper threshold.
 * @zl: minimum redshift.
 * @zu: maximum redshift.
 * @spline: Whenever to create an intermediary spline of the integration.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_mass_function_N (NcMassFunction * mfp, NcHICosmo * model, gdouble lnMl, gdouble lnMu, gdouble zl, gdouble zu,
		    NcMassFunctionSplineOptimize spline)
{
  const gboolean has_d2NdzdlnM_spline = (mfp->d2NdzdlnM != NULL);
  const gboolean in_spline = has_d2NdzdlnM_spline && (lnMl >= mfp->lnMi) && (lnMu <= mfp->lnMf) && (zl >= mfp->zi)
    && (zu <= mfp->zf);
  gdouble N;
  if (in_spline)
  {
    if (ncm_model_ctrl_update (mfp->ctrl, NCM_MODEL(model)))
    {
      nc_mass_function_d2NdzdlnM_prepare (mfp, model);
    }
    switch (spline)
    {
      case NC_MASS_FUNCTION_SPLINE_NONE:
	N = ncm_spline2d_integ_dxdy (mfp->d2NdzdlnM, lnMl, lnMu, zl, zu);
	break;
      case NC_MASS_FUNCTION_SPLINE_LNM:
	N = ncm_spline2d_integ_dxdy_spline_x (mfp->d2NdzdlnM, lnMl, lnMu, zl, zu);
	break;
      case NC_MASS_FUNCTION_SPLINE_Z:
	N = ncm_spline2d_integ_dxdy_spline_y (mfp->d2NdzdlnM, lnMl, lnMu, zl, zu);
	break;
      default:
	g_assert_not_reached ();
	break;
    }
  }
  else
    g_error ("This option is not implemented. It must be used the optimize option.\n");
  return N;
}

typedef struct __encapsulated_function_args
{
  NcMassFunction *mfp;
  NcHICosmo *model;
  gdouble z;
  gdouble lnM;
  gdouble dVdz;
} _encapsulated_function_args;

static gdouble
_encapsulated_z (gdouble z, gpointer p)
{
  _encapsulated_function_args *args = (_encapsulated_function_args *) p;

  gdouble A = args->mfp->area_survey *
    nc_mass_function_dcomoving_volume_dzdomega (args->mfp, args->model, z) *
    nc_mass_function (args->mfp, args->model, args->lnM, z);

  //printf ("1 z %.5g lnM %.5g d2N % 20.15g mf % 20.15g\n", z, args->lnM, A, nc_mass_function (args->mfp, args->model, args->lnM, z));
  return A;

  //return args->mfp->area_survey *
  //nc_dcomoving_volume_dzdomega (args->mfp, args->model, z) *
  //nc_mass_function (args->mfp, args->model, args->lnM, z);
}

static gdouble
_encapsulated_lnM (gdouble lnM, gpointer p)
{
  _encapsulated_function_args *args = (_encapsulated_function_args *) p;

  gdouble A = args->dVdz * nc_mass_function (args->mfp, args->model, lnM, args->z);

  //printf ("2 z %.5g lnM %.5g d2N % 20.15g\n", args->z, lnM, A);
  return A;

  //return args->dVdz * nc_mass_function (args->mfp, args->model, lnM, args->z);
}

static void
_nc_mass_function_generate_2Dspline_knots (NcMassFunction * mfp, NcHICosmo * model, gdouble rel_error)
{
  gsl_function Fx, Fy;
  _encapsulated_function_args args;
  g_assert (mfp->d2NdzdlnM == NULL);

  args.mfp = mfp;
  args.model = model;
  args.z = (mfp->zf + mfp->zi) / 2.0;
  args.lnM = (mfp->lnMf + mfp->lnMi) / 2.0;
  args.dVdz = mfp->area_survey * nc_mass_function_dcomoving_volume_dzdomega (mfp, model, args.z);

  Fx.function = _encapsulated_lnM;
  Fx.params = &args;

  Fy.function = _encapsulated_z;
  Fy.params = &args;


  mfp->d2NdzdlnM = ncm_spline2d_bicubic_notaknot_new ();
  ncm_spline2d_set_function (mfp->d2NdzdlnM,
			     NCM_SPLINE_FUNCTION_SPLINE, &Fx, &Fy, mfp->lnMi, mfp->lnMf, mfp->zi, mfp->zf, 1e-5);
}
