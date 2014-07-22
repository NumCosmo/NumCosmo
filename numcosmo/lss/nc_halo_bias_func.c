/***************************************************************************
 *            nc_halo_bias_func.c
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
 * SECTION:nc_halo_bias_func
 * @title: Mean Halo Bias Function
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_halo_bias_func.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_MASS_FUNCTION,
  PROP_BIAS_TYPE
};

G_DEFINE_TYPE (NcHaloBiasFunc, nc_halo_bias_func, G_TYPE_OBJECT);

/**
 * nc_halo_bias_func_new:
 * @mfp: a #NcMassFunction.
 * @biasf: (allow-none): a #NcHaloBiasType.
 *
 * This function allocates memory for a new #NcHaloBiasFunc object and sets its properties to the values from
 * the input arguments.
 *
 * Returns: A new #NcHaloBiasFunc.
 */
NcHaloBiasFunc *
nc_halo_bias_func_new (NcMassFunction *mfp, NcHaloBiasType *biasf)
{
  NcHaloBiasFunc *mbiasf = g_object_new (NC_TYPE_HALO_BIAS_FUNC,
                                     "mass-function", mfp,
                                     "bias-type", biasf,
                                     NULL);
  return mbiasf;
}

/**
 * nc_halo_bias_func_copy:
 * @mbiasf: a #NcHaloBiasFunc.
 *
 * Duplicates the #NcHaloBiasFunc object setting the same values of the original propertities.
 *
 * Returns: (transfer full): A new #NcHaloBiasFunc.
 */
NcHaloBiasFunc *
nc_halo_bias_func_copy (NcHaloBiasFunc *mbiasf)
{
  return nc_halo_bias_func_new (mbiasf->mfp, mbiasf->biasf);
}

/**
 * nc_halo_bias_func_free:
 * @mbiasf: a #NcHaloBiasFunc.
 *
 * Atomically decrements the reference count of @mbiasf by one. If the reference count drops to 0,
 * all memory allocated by @mbiasf is released.
 *
*/
void
nc_halo_bias_func_free (NcHaloBiasFunc *mbiasf)
{
  g_object_unref (mbiasf);
}

/**
 * nc_halo_bias_func_clear:
 * @mbiasf: a #NcHaloBiasFunc.
 *
 * Atomically decrements the reference count of @mbiasf by one. If the reference count drops to 0,
 * all memory allocated by @mbiasf is released. Set pointer to NULL.
 *
*/
void
nc_halo_bias_func_clear (NcHaloBiasFunc **mbiasf)
{
  g_clear_object (mbiasf);
}

static void
nc_halo_bias_func_init (NcHaloBiasFunc *mbiasf)
{
  mbiasf->mfp = NULL;
  mbiasf->biasf = NULL;
}

static void
_nc_halo_bias_func_dispose (GObject *object)
{
  NcHaloBiasFunc *mbiasf = NC_HALO_BIAS_FUNC (object);
  
  nc_mass_function_clear (&mbiasf->mfp);
  nc_halo_bias_type_clear (&mbiasf->biasf);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_bias_func_parent_class)->dispose (object);
}

static void
_nc_halo_bias_func_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_bias_func_parent_class)->finalize (object);
}

static void
_nc_halo_bias_func_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloBiasFunc *mbiasf = NC_HALO_BIAS_FUNC (object);
  g_return_if_fail (NC_IS_HALO_BIAS_FUNC (object));

  switch (prop_id)
  {
	case PROP_MASS_FUNCTION:
	  mbiasf->mfp = g_value_dup_object (value);
	  break;
	case PROP_BIAS_TYPE:
	  mbiasf->biasf = g_value_dup_object (value);
	  break;
	default:
	  G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
	  break;
  }
}

static void
_nc_halo_bias_func_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloBiasFunc *mbiasf = NC_HALO_BIAS_FUNC (object);
  g_return_if_fail (NC_IS_HALO_BIAS_FUNC (object));

  switch (prop_id)
  {
	case PROP_MASS_FUNCTION:
	  g_value_set_object (value, mbiasf->mfp);
	  break;
	case PROP_BIAS_TYPE:
	  g_value_set_object (value, mbiasf->biasf);
	  break;
	default:
	  G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
	  break;
  }
}

static void
nc_halo_bias_func_class_init (NcHaloBiasFuncClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  //GObjectClass* parent_class = G_OBJECT_CLASS (klass);

  object_class->dispose = _nc_halo_bias_func_dispose;
  object_class->finalize = _nc_halo_bias_func_finalize;
  object_class->set_property = _nc_halo_bias_func_set_property;
  object_class->get_property = _nc_halo_bias_func_get_property;

  /**
   * NcHaloBiasFunc:mass-function:
   *
   * This property keeps the mass function object.
   */
  g_object_class_install_property (object_class,
                                   PROP_MASS_FUNCTION,
                                   g_param_spec_object ("mass-function",
                                                        NULL,
                                                        "Mass Function.",
                                                        NC_TYPE_MASS_FUNCTION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /*
   * NcHaloBiasFunc:bias-type:
   *
   * FIXME
   */
   g_object_class_install_property (object_class,
                                    PROP_BIAS_TYPE,
                                    g_param_spec_object ("bias-type",
                                                         NULL,
                                                         "Bias Function Type.",
                                                         NC_TYPE_HALO_BIAS_TYPE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

}

/**
 * nc_halo_bias_func_integrand:
 * @mbiasf: a #NcHaloBiasFunc.
 * @model: a #NcHICosmo.
 * @lnM: logarithm base e of the mass.
 * @z: redshift.
 *
 * This function is the integrand of the mean bias, i.e., the product of the mass function with the bias function.
 * As both functions depend on the standard deviation of the matter density contrast, we implement this function to
 * compute \f$ \sigma (M, z) \f$ just once.
   *
 * It is worth noting that the multiplicity function must be compatible with the bias function.
 *
 * Returns: a double which corresponds to the mean bias integrand for lnM and at redshift z.
 */
gdouble
nc_halo_bias_func_integrand (NcHaloBiasFunc *mbiasf, NcHICosmo *model, gdouble lnM, gdouble z)
{
  gdouble lnR = nc_matter_var_lnM_to_lnR (mbiasf->mfp->vp, model, lnM);
  gdouble R = exp (lnR);
  gdouble M_rho = nc_window_volume (mbiasf->mfp->vp->wp) * gsl_pow_3(R);   /* M/\rho comoving, rho independe de z */
  gdouble var0 = nc_matter_var_var0 (mbiasf->mfp->vp, model, lnR);
  //printf("var0 = %.5g\n", var0);
  gdouble dlnvar0_dlnR = nc_matter_var_dlnvar0_dlnR (mbiasf->mfp->vp, model, lnR);
  //printf("dlnvar0_dlnR = %.5g\n", dlnvar0_dlnR);
  gdouble growth = nc_growth_func_eval (mbiasf->mfp->gf, model, z);
  //printf("growth = %.5g\n", growth);
  gdouble sigma = nc_matter_var_sigma8_sqrtvar0 (mbiasf->mfp->vp, model) * sqrt(var0) * growth;
  gdouble f = nc_multiplicity_func_eval (mbiasf->mfp->mulf, model, sigma, z);
  gdouble dn_dlnM = -(1.0 / (3.0 * M_rho)) * f * 0.5 * dlnvar0_dlnR;  /* dn/dlnM = - (\rho/3M) * f * (R/\sigma)* dsigma_dR */
  gdouble bias = nc_halo_bias_type_eval (mbiasf->biasf, sigma, z);

  /*printf ("%.5g %.5g %.5g %.5g %.5g\n", lnM, lnR / M_LN10, -(1.0 / (3.0 * M_rho)), f, 0.5 * dlnvar0_dlnR);*/
  /*printf("log_nu = %.5g sigma = %.5g f = %.5g n = %.5g log_b = %.5g\n", log10(1.686/sigma), sigma, f, dn_dlnM, log10(bias));*/

  /*printf("nu = %.5g logsigma = %.5g sigma = %.5g f = %.5g b = %.5g\n", log10(1.686/sigma), log10(sigma), sigma, f, bias);*/

  return dn_dlnM * bias;
}
