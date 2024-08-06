/***************************************************************************
 *            nc_xcor_limber_kernel_cmb_isw.c
 *
 *  Tue Sept 28 17:17:26 2021
 *  Copyright  2021  Mariana Penna-Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2021 Mariana Penna-Lima  <pennalima@gmail.com>
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
 * SECTION:nc_xcor_limber_kernel_cmb_isw
 * @title: NcXcorLimberKernelCMBISW
 * @short_description: implementation of #NcXcorLimberKernel for integrated Sachs-Wolfe (ISW)
 *
 * The kernel is given by
 * \begin{equation}
 *    W^{T_\mathrm{ISW}} (z) = \frac{3 \Omega_m H_0^2}{c^2 (l + 1/2)^2} \chi^2(z) \frac{\rm{d}}{\rm{d}z}\left((1 + z)D(z)\right).
 * \end{equation}
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_cfg.h"
#include "xcor/nc_xcor_limber_kernel_cmb_isw.h"
#include "xcor/nc_xcor.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorLimberKernelCMBISW
{
  /*< private >*/
  NcXcorLimberKernel parent_instance;
};

typedef struct _NcXcorLimberKernelCMBISWPrivate
{
  NcDistance *dist;
  NcmPowspec *ps;
  NcRecomb *recomb;
  NcmVector *Nl;
  guint Nlmax;
  gdouble xi_lss;
  gdouble cons_factor;
} NcXcorLimberKernelCMBISWPrivate;

enum
{
  PROP_0,
  PROP_DIST,
  PROP_PS,
  PROP_RECOMB,
  PROP_NL,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcXcorLimberKernelCMBISW, nc_xcor_limber_kernel_cmb_isw, NC_TYPE_XCOR_LIMBER_KERNEL);

static void
nc_xcor_limber_kernel_cmb_isw_init (NcXcorLimberKernelCMBISW *xcisw)
{
  NcXcorLimberKernelCMBISWPrivate * const self = nc_xcor_limber_kernel_cmb_isw_get_instance_private (xcisw);

  self->dist   = NULL;
  self->ps     = NULL;
  self->recomb = NULL;
  self->Nl     = NULL;
  self->Nlmax  = 0;
  self->xi_lss = 0.0;
}

static void
_nc_xcor_limber_kernel_cmb_isw_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorLimberKernelCMBISW *xcisw              = NC_XCOR_LIMBER_KERNEL_CMB_ISW (object);
  NcXcorLimberKernelCMBISWPrivate * const self = nc_xcor_limber_kernel_cmb_isw_get_instance_private (xcisw);

  g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL_CMB_ISW (object));

  switch (prop_id)
  {
    case PROP_DIST:
      self->dist = g_value_dup_object (value);
      break;
    case PROP_PS:
      self->ps = g_value_dup_object (value);
      break;
    case PROP_RECOMB:
      self->recomb = g_value_dup_object (value);
      break;
    case PROP_NL:
      self->Nl    = g_value_dup_object (value);
      self->Nlmax = ncm_vector_len (self->Nl) - 1;
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_xcor_limber_kernel_cmb_isw_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorLimberKernelCMBISW *xcisw              = NC_XCOR_LIMBER_KERNEL_CMB_ISW (object);
  NcXcorLimberKernelCMBISWPrivate * const self = nc_xcor_limber_kernel_cmb_isw_get_instance_private (xcisw);

  g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL_CMB_ISW (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, self->dist);
      break;
    case PROP_PS:
      g_value_set_object (value, self->ps);
      break;
    case PROP_RECOMB:
      g_value_set_object (value, self->recomb);
      break;
    case PROP_NL:
      g_value_set_object (value, self->Nl);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_xcor_limber_kernel_cmb_isw_dispose (GObject *object)
{
  NcXcorLimberKernelCMBISW *xcisw              = NC_XCOR_LIMBER_KERNEL_CMB_ISW (object);
  NcXcorLimberKernelCMBISWPrivate * const self = nc_xcor_limber_kernel_cmb_isw_get_instance_private (xcisw);

  nc_distance_clear (&self->dist);
  ncm_powspec_clear (&self->ps);
  nc_recomb_clear (&self->recomb);
  ncm_vector_clear (&self->Nl);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_limber_kernel_cmb_isw_parent_class)->dispose (object);
}

static void
_nc_xcor_limber_kernel_cmb_isw_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_limber_kernel_cmb_isw_parent_class)->finalize (object);
}

static gdouble _nc_xcor_limber_kernel_cmb_isw_eval (NcXcorLimberKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
static void _nc_xcor_limber_kernel_cmb_isw_prepare (NcXcorLimberKernel *xclk, NcHICosmo *cosmo);
static void _nc_xcor_limber_kernel_cmb_isw_add_noise (NcXcorLimberKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
static guint _nc_xcor_limber_kernel_cmb_isw_obs_len (NcXcorLimberKernel *xclk);
static guint _nc_xcor_limber_kernel_cmb_isw_obs_params_len (NcXcorLimberKernel *xclk);

static void
nc_xcor_limber_kernel_cmb_isw_class_init (NcXcorLimberKernelCMBISWClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcXcorLimberKernelClass *parent_class = NC_XCOR_LIMBER_KERNEL_CLASS (klass);
  NcmModelClass *model_class            = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_xcor_limber_kernel_cmb_isw_set_property;
  model_class->get_property = &_nc_xcor_limber_kernel_cmb_isw_get_property;
  object_class->finalize    = &_nc_xcor_limber_kernel_cmb_isw_finalize;
  object_class->dispose     = &_nc_xcor_limber_kernel_cmb_isw_dispose;

  ncm_model_class_set_name_nick (model_class, "Xcor ISW effect", "Xcor-ISW");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorLimberKernelCMBISW:dist:
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
   * NcXcorLimberKernelCMBISW:ps:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_PS,
                                   g_param_spec_object ("ps",
                                                        NULL,
                                                        "Power Spectrum object",
                                                        NCM_TYPE_POWSPEC,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  /**
   * NcXcorLimberKernelCMBISW:recomb:
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
   * NcXcorLimberKernelCMBISW:Nl:
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

  parent_class->eval      = &_nc_xcor_limber_kernel_cmb_isw_eval;
  parent_class->prepare   = &_nc_xcor_limber_kernel_cmb_isw_prepare;
  parent_class->add_noise = &_nc_xcor_limber_kernel_cmb_isw_add_noise;

  parent_class->obs_len        = &_nc_xcor_limber_kernel_cmb_isw_obs_len;
  parent_class->obs_params_len = &_nc_xcor_limber_kernel_cmb_isw_obs_params_len;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_LIMBER_KERNEL_IMPL_ALL);
}

static gdouble
_nc_xcor_limber_kernel_cmb_isw_eval (NcXcorLimberKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  NcXcorLimberKernelCMBISW *xcisw              = NC_XCOR_LIMBER_KERNEL_CMB_ISW (xclk);
  NcXcorLimberKernelCMBISWPrivate * const self = nc_xcor_limber_kernel_cmb_isw_get_instance_private (xcisw);
  const gdouble k_pivot                        = 1.0;
  const gdouble powspec                        = ncm_powspec_eval (NCM_POWSPEC (self->ps), NCM_MODEL (cosmo), z, k_pivot);
  const gdouble dpowspec_dz                    = ncm_powspec_deriv_z (NCM_POWSPEC (self->ps), NCM_MODEL (cosmo), z, k_pivot);
  const gdouble d1pz_growth_dz                 = 1.0 + (1.0 + z) * dpowspec_dz / (2.0 * powspec);
  const gdouble nu                             = l + 0.5;
  const gdouble cor_factor                     = 1.0 / (nu * nu);

  return cor_factor * gsl_pow_2 (xck->xi_z) * d1pz_growth_dz;
}

static void
_nc_xcor_limber_kernel_cmb_isw_prepare (NcXcorLimberKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorLimberKernelCMBISW *xcisw              = NC_XCOR_LIMBER_KERNEL_CMB_ISW (xclk);
  NcXcorLimberKernelCMBISWPrivate * const self = nc_xcor_limber_kernel_cmb_isw_get_instance_private (xcisw);

  nc_distance_prepare_if_needed (self->dist, cosmo);
  ncm_powspec_prepare_if_needed (self->ps, NCM_MODEL (cosmo));
  {
    const gdouble z_lss    = nc_distance_decoupling_redshift (self->dist, cosmo);
    const gdouble Omega_c0 = nc_hicosmo_Omega_c0 (cosmo);
    const gdouble Omega_b0 = nc_hicosmo_Omega_b0 (cosmo);
    const gdouble Omega_m0 = Omega_c0 + Omega_b0;
    const gdouble T_gamma0 = nc_hicosmo_T_gamma0 (cosmo);

    self->xi_lss      = nc_distance_comoving_lss (self->dist, cosmo);
    self->cons_factor = (3.0 * nc_hicosmo_Omega_m0 (cosmo));

    nc_xcor_limber_kernel_set_const_factor (xclk, 3.0 * T_gamma0 * Omega_m0);
    nc_xcor_limber_kernel_set_z_range (xclk, 0.0, z_lss, 2.0);
  }
}

static void
_nc_xcor_limber_kernel_cmb_isw_add_noise (NcXcorLimberKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  NcXcorLimberKernelCMBISW *xcisw              = NC_XCOR_LIMBER_KERNEL_CMB_ISW (xclk);
  NcXcorLimberKernelCMBISWPrivate * const self = nc_xcor_limber_kernel_cmb_isw_get_instance_private (xcisw);

  if (self->Nl == NULL)
    g_error ("nc_xcor_limber_kernel_cmb_isw_noise_spec : noise spectrum empty");

  if (lmin + ncm_vector_len (vp1) > self->Nlmax)
    g_error ("nc_xcor_limber_kernel_cmb_isw_noise_spec : too high multipole");

  ncm_vector_memcpy (vp2, vp1);

  {
    NcmVector *Nl_sub = ncm_vector_get_subvector (self->Nl, lmin, ncm_vector_len (vp1));

    ncm_vector_add (vp2, Nl_sub);
    ncm_vector_free (Nl_sub);
  }
  /* return ncm_vector_get (xcisw->Nl, l); */
}

static guint
_nc_xcor_limber_kernel_cmb_isw_obs_len (NcXcorLimberKernel *xclk)
{
  NCM_UNUSED (xclk);

  return 1;
}

static guint
_nc_xcor_limber_kernel_cmb_isw_obs_params_len (NcXcorLimberKernel *xclk)
{
  NCM_UNUSED (xclk);

  return 0;
}

/**
 * nc_xcor_limber_kernel_cmb_isw_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @recomb: a #NcRecomb
 * @Nl: a #NcmVector
 *
 * FIXME
 *
 * Returns: FIXME
 *
 */
NcXcorLimberKernelCMBISW *
nc_xcor_limber_kernel_cmb_isw_new (NcDistance *dist, NcmPowspec *ps, NcRecomb *recomb, NcmVector *Nl)
{
  NcXcorLimberKernelCMBISW *xcisw = g_object_new (NC_TYPE_XCOR_LIMBER_KERNEL_CMB_ISW,
                                                  "dist", dist,
                                                  "ps", ps,
                                                  "recomb", recomb,
                                                  "Nl", Nl,
                                                  NULL);

  return xcisw;
}

