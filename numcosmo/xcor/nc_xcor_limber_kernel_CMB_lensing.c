/***************************************************************************
 *            nc_xcor_limber_kernel_cmb_lensing.c
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
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
 * SECTION:nc_xcor_limber_kernel_CMB_lensing
 * @title: NcXcorLimberKernelCMBLensing
 * @short_description: implementation of #NcXcorLimberKernel for CMB lensing
 *
 * The kernel is given by
 * \begin{equation}
 *    W^{\kappa_\mathrm{CMB}} (z) = \frac{3}{2} \frac{\Omega_m H_0^2}{c} \frac{(1+z)}{H(z)} \chi(z) \frac{\chi(z_*) - \chi(z)}{\chi(z_*)}.
 * \end{equation}
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_cfg.h"
#include "xcor/nc_xcor_limber_kernel_CMB_lensing.h"
#include "xcor/nc_xcor.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_TYPE (NcXcorLimberKernelCMBLensing, nc_xcor_limber_kernel_cmb_lensing, NC_TYPE_XCOR_LIMBER_KERNEL);

#define VECTOR (NCM_MODEL (xclkl)->params)

enum
{
  PROP_0,
  PROP_DIST,
  PROP_RECOMB,
  PROP_NL,
  PROP_SIZE,
};

static void
nc_xcor_limber_kernel_cmb_lensing_init (NcXcorLimberKernelCMBLensing* xclkl)
{
  xclkl->dist   = NULL;
  xclkl->recomb = NULL;

  xclkl->Nl     = NULL;
  xclkl->Nlmax  = 0;
}

static void
_nc_xcor_limber_kernel_cmb_lensing_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
  NcXcorLimberKernelCMBLensing* xclkl = NC_XCOR_LIMBER_KERNEL_CMB_LENSING (object);
  g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL_CMB_LENSING (object));

  switch (prop_id)
  {
    case PROP_DIST:
      xclkl->dist = g_value_dup_object (value);
      break;
    case PROP_RECOMB:
      xclkl->recomb = g_value_dup_object (value);
      break;
    case PROP_NL:
      xclkl->Nl = g_value_dup_object (value);
      xclkl->Nlmax = ncm_vector_len (xclkl->Nl) - 1;
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_xcor_limber_kernel_cmb_lensing_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorLimberKernelCMBLensing *xclkl = NC_XCOR_LIMBER_KERNEL_CMB_LENSING (object);
  g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL_CMB_LENSING (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, xclkl->dist);
      break;
    case PROP_RECOMB:
      g_value_set_object (value, xclkl->recomb);
      break;
    case PROP_NL:
      g_value_set_object (value, xclkl->Nl);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_xcor_limber_kernel_cmb_lensing_dispose (GObject *object)
{
  NcXcorLimberKernelCMBLensing *xclkl = NC_XCOR_LIMBER_KERNEL_CMB_LENSING (object);

  nc_distance_clear (&xclkl->dist);
  nc_recomb_clear (&xclkl->recomb);
  ncm_vector_clear (&xclkl->Nl);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_limber_kernel_cmb_lensing_parent_class)->dispose (object);
}

static void
_nc_xcor_limber_kernel_cmb_lensing_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_limber_kernel_cmb_lensing_parent_class)->finalize (object);
}

static gdouble _nc_xcor_limber_kernel_cmb_lensing_eval (NcXcorLimberKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
static void _nc_xcor_limber_kernel_cmb_lensing_prepare (NcXcorLimberKernel *xclk, NcHICosmo *cosmo);
static void _nc_xcor_limber_kernel_cmb_lensing_add_noise (NcXcorLimberKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
static guint _nc_xcor_limber_kernel_cmb_lensing_obs_len (NcXcorLimberKernel *xclk);
static guint _nc_xcor_limber_kernel_cmb_lensing_obs_params_len (NcXcorLimberKernel *xclk);

static void
nc_xcor_limber_kernel_cmb_lensing_class_init (NcXcorLimberKernelCMBLensingClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcXcorLimberKernelClass* parent_class = NC_XCOR_LIMBER_KERNEL_CLASS (klass);
  NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize    = &_nc_xcor_limber_kernel_cmb_lensing_finalize;
  object_class->dispose     = &_nc_xcor_limber_kernel_cmb_lensing_dispose;
  model_class->set_property = &_nc_xcor_limber_kernel_cmb_lensing_set_property;
  model_class->get_property = &_nc_xcor_limber_kernel_cmb_lensing_get_property;

  ncm_model_class_set_name_nick (model_class, "Xcor lensing distribution", "Xcor-lensing");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorLimberKernelCMBLensing:dist:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorLimberKernelCMBLensing:recomb:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_RECOMB,
                                   g_param_spec_object ("recomb",
                                                        NULL,
                                                        "Recombination object",
                                                        NC_TYPE_RECOMB,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorLimberKernelCMBLensing:Nl:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_NL,
                                   g_param_spec_object ("Nl",
                                                        NULL,
                                                        "Noise spectrum",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->eval           = &_nc_xcor_limber_kernel_cmb_lensing_eval;
  parent_class->prepare        = &_nc_xcor_limber_kernel_cmb_lensing_prepare;
  /*parent_class->noise_spec = &_nc_xcor_limber_kernel_cmb_lensing_noise_spec;*/
  parent_class->add_noise      = &_nc_xcor_limber_kernel_cmb_lensing_add_noise;

  parent_class->obs_len        = &_nc_xcor_limber_kernel_cmb_lensing_obs_len;
  parent_class->obs_params_len = &_nc_xcor_limber_kernel_cmb_lensing_obs_params_len;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_LIMBER_KERNEL_IMPL_ALL);
}

/**
 * nc_xcor_limber_kernel_cmb_lensing_new:
 * @dist: a #NcDistance
 * @recomb: a #NcRecomb
 * @Nl: a #NcmVector
 *
 * FIXME
 *
 * Returns: FIXME
 *
 */
NcXcorLimberKernelCMBLensing *
nc_xcor_limber_kernel_cmb_lensing_new (NcDistance* dist, NcRecomb* recomb, NcmVector* Nl) //, gdouble zl, gdouble zu)
{
  NcXcorLimberKernelCMBLensing *xclkl = g_object_new (NC_TYPE_XCOR_LIMBER_KERNEL_CMB_LENSING,
                                                      "dist", dist,
                                                      "recomb", recomb,
                                                      "Nl", Nl,
                                                      NULL);
  return xclkl;
}

static gdouble
_nc_xcor_limber_kernel_cmb_lensing_eval (NcXcorLimberKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)//, gdouble geo_z[])
{
  NcXcorLimberKernelCMBLensing* xclkl = NC_XCOR_LIMBER_KERNEL_CMB_LENSING (xclk);

  NCM_UNUSED (l);

  return ((1.0 + z) * xck->xi_z * (xclkl->xi_lss - xck->xi_z)) / (xck->E_z * xclkl->xi_lss);
}

static void
_nc_xcor_limber_kernel_cmb_lensing_prepare (NcXcorLimberKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorLimberKernelCMBLensing* xclkl = NC_XCOR_LIMBER_KERNEL_CMB_LENSING (xclk);

  nc_distance_prepare_if_needed (xclkl->dist, cosmo);

  xclkl->xi_lss     = nc_distance_comoving_lss (xclkl->dist, cosmo);
  xclk->cons_factor = (3.0 * nc_hicosmo_Omega_m0 (cosmo)) / 2.0;

  // nc_recomb_prepare (xclkl->recomb, cosmo);
  // gdouble lamb = nc_recomb_tau_zstar (xclkl->recomb, cosmo);

  xclk->zmax = 1090.0; //exp (-lamb) - 1.0;
  xclk->zmin = 0.0;
}

static void
_nc_xcor_limber_kernel_cmb_lensing_add_noise (NcXcorLimberKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  NcXcorLimberKernelCMBLensing* xclkl = NC_XCOR_LIMBER_KERNEL_CMB_LENSING (xclk);

  if (xclkl->Nl == NULL)
    g_error ("nc_xcor_limber_kernel_cmb_lensing_noise_spec : noise spectrum empty");

  if (lmin + ncm_vector_len(vp1) > xclkl->Nlmax)
    g_error ("nc_xcor_limber_kernel_cmb_lensing_noise_spec : too high multipole");

  ncm_vector_memcpy(vp2, vp1);

  {
    NcmVector *Nl_sub = ncm_vector_get_subvector(xclkl->Nl, lmin, ncm_vector_len(vp1));

    ncm_vector_add (vp2, Nl_sub);
    ncm_vector_free (Nl_sub);
  }
  // return ncm_vector_get (xclkl->Nl, l);
}

static guint
_nc_xcor_limber_kernel_cmb_lensing_obs_len (NcXcorLimberKernel *xclk)
{
  NCM_UNUSED (xclk);
  return 1;
}

static guint
_nc_xcor_limber_kernel_cmb_lensing_obs_params_len (NcXcorLimberKernel *xclk)
{
  NCM_UNUSED (xclk);
  return 0;
}
