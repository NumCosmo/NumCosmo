/***************************************************************************
 *            nc_halo_bias.c
 *
 *  Tue June 28 15:41:57 2011
 *  Copyright  2011  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION:nc_halo_bias
 * @title: NcHaloBias
 * @short_description: Abstract class for halo bias function type.
 *
 * Abstract object to compute the halo bias function given
 * a mass function. The bias function is given by the
 * eval function, which should be implemented by the child
 * objects of the #NcHaloBias class. Below there is a quick
 * review about dark matter halo bias.
 *
 * Dark Matter Halos are tracers of the Large Scale Structure (LSS)
 * of the universe. On sufficient large scales, the correlations
 * are weak and these objects are expected to be well described by
 * perturbation theory. Given the background matter distribution, 
 * $\rho$, one can assume that the mass distribution in the universe
 * is given by a small perturbation around this value. In this context, 
 * we define the density contrast $\delta$ as 
 * 
 * \begin{align}
 * \delta_b &= \frac{\delta \rho}{\rho}
 * ,\end{align}  
 * which represents this perturbation around the background quantity. 
 *
 * The halos are formed when the density contrast passes a certain
 * critical limit, that is, $\delta_b > \delta_c$. Even though
 * $\delta_b$ is a gaussian random field, the real cluster formation
 * is a rare event that involves non-linear mechanisms that cannot
 * be explained perturbatively, and thus we must consider Peak-Background split
 * effects. To connect the theory to observations, we introduce "bias parameters"
 * in the observed halos distributions.
 *
 * These bias parameters are linear corrections to the density contrast
 * field, such that 
 * \begin{align}
 * \delta &= [1 + b(z, M)] \delta_b
 * ,\end{align}
 * where $\delta$ is the new density contrast that accounts for the bias and
 * $b(z, M)$ is the bias. There different biases for each mass function
 * used in the cosmological computations. Check the child objects
 * for their definitions.
 *
 * Since this is an abstract class implementation, there is no initialization
 * function of the object.  
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_halo_bias.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

G_DEFINE_ABSTRACT_TYPE (NcHaloBias, nc_halo_bias, G_TYPE_OBJECT)

enum
{
  PROP_0,
  PROP_MASS_FUNCTION,
  PROP_SIZE,
};

static void
nc_halo_bias_init (NcHaloBias *bias)
{
  bias->mfp = NULL;
}

static void
_nc_halo_bias_dispose (GObject *object)
{
  NcHaloBias *bias = NC_HALO_BIAS   (object);

  nc_halo_mass_function_clear (&bias->mfp);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_bias_parent_class)->dispose (object);
}

static void
_nc_halo_bias_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_bias_parent_class)->finalize (object);
}

static void
_nc_halo_bias_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloBias *bias = NC_HALO_BIAS (object);

  g_return_if_fail (NC_IS_HALO_BIAS (object));

  switch (prop_id)
  {
    case PROP_MASS_FUNCTION:
      bias->mfp = g_value_dup_object (value);
      g_assert (bias->mfp != NULL);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_halo_bias_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloBias *bias = NC_HALO_BIAS (object);

  g_return_if_fail (NC_IS_HALO_BIAS (object));

  switch (prop_id)
  {
    case PROP_MASS_FUNCTION:
      g_value_set_object (value, bias->mfp);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_halo_bias_class_init (NcHaloBiasClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  /*GObjectClass* parent_class = G_OBJECT_CLASS (klass); */

  object_class->dispose      = &_nc_halo_bias_dispose;
  object_class->finalize     = &_nc_halo_bias_finalize;
  object_class->set_property = &_nc_halo_bias_set_property;
  object_class->get_property = &_nc_halo_bias_get_property;

  /**
   * NcHaloBias:mass-function:
   *
   * This property keeps the mass function object.
   */
  g_object_class_install_property (object_class,
                                   PROP_MASS_FUNCTION,
                                   g_param_spec_object ("mass-function",
                                                        NULL,
                                                        "Mass Function.",
                                                        NC_TYPE_HALO_MASS_FUNCTION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_halo_bias_eval:
 * @bias: a #NcHaloBias
 * @cosmo: a #NcHICosmo
 * @sigma: density constrast variance
 * @z: redshift
 *
 * Computes the Halo Bias at a given redshift.
 *
 * Returns: a double, the halo bias.
 */
gdouble
nc_halo_bias_eval (NcHaloBias *bias, NcHICosmo *cosmo, gdouble sigma, gdouble z, gdouble lnM)
{
  return NC_HALO_BIAS_GET_CLASS (bias)->eval (bias, cosmo, sigma, z, lnM);
}

/**
 * nc_halo_bias_free:
 * @bias: a #NcHaloBias
 *
 * Atomically decrements the reference count of @bias by one. If the reference count drops to 0,
 * all memory allocated by @bias is released.
 *
 */
void
nc_halo_bias_free (NcHaloBias *bias)
{
  g_object_unref (bias);
}

/**
 * nc_halo_bias_clear:
 * @bias: a #NcHaloBias
 *
 * Atomically decrements the reference count of @bias by one. If the reference count drops to 0,
 * all memory allocated by @bias is released. Set pointer to NULL.
 *
 */
void
nc_halo_bias_clear (NcHaloBias **bias)
{
  g_clear_object (bias);
}

/**
 * nc_halo_bias_integrand:
 * @mbiasf: a #NcHaloBias
 * @cosmo: a #NcHICosmo
 * @lnM: logarithm base e of the mass
 * @z: redshift
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
nc_halo_bias_integrand (NcHaloBias *mbiasf, NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
  const gdouble d2n_dzdlnM = nc_halo_mass_function_d2n_dzdlnM (mbiasf->mfp, cosmo, lnM, z);
  const gdouble sigma      = nc_halo_mass_function_sigma_lnM (mbiasf->mfp, cosmo, lnM, z);
  const gdouble bias       = nc_halo_bias_eval (mbiasf, cosmo, sigma, z, lnM);

  return d2n_dzdlnM * bias;
}
