/***************************************************************************
 *            nc_xcor_kernel_cmb_isw.c
 *
 *  Tue Sept 28 17:17:26 2021
 *  Copyright  2021  Mariana Penna-Lima
 *  <pennalima@gmail.com>
 *  Sat December 27 20:21:01 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2021 Mariana Penna-Lima  <pennalima@gmail.com>
 * Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcXcorKernelCMBISW:
 *
 * Implementation of #NcXcorKernel for integrated Sachs-Wolfe (ISW).
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
#include "xcor/nc_xcor_kernel_cmb_isw.h"
#include "xcor/nc_xcor.h"


#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorKernelCMBISW
{
  /*< private >*/
  NcXcorKernel parent_instance;
};

typedef struct _NcXcorKernelCMBISWPrivate
{
  NcRecomb *recomb;
  NcmVector *Nl;
  guint Nlmax;
  gdouble xi_lss;
  NcDistance *dist;
  NcmPowspec *ps;
} NcXcorKernelCMBISWPrivate;

enum
{
  PROP_0,
  PROP_RECOMB,
  PROP_NL,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcXcorKernelCMBISW, nc_xcor_kernel_cmb_isw, NC_TYPE_XCOR_KERNEL);

static void
nc_xcor_kernel_cmb_isw_init (NcXcorKernelCMBISW *xcisw)
{
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  self->recomb = NULL;
  self->Nl     = NULL;
  self->Nlmax  = 0;
  self->xi_lss = 0.0;
  self->dist   = NULL;
}

static void
_nc_xcor_kernel_cmb_isw_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelCMBISW *xcisw              = NC_XCOR_KERNEL_CMB_ISW (object);
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  g_return_if_fail (NC_IS_XCOR_KERNEL_CMB_ISW (object));

  switch (prop_id)
  {
    case PROP_RECOMB:
      self->recomb = g_value_dup_object (value);
      break;
    case PROP_NL:
      self->Nl    = g_value_dup_object (value);
      self->Nlmax = ncm_vector_len (self->Nl) - 1;
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_kernel_cmb_isw_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelCMBISW *xcisw              = NC_XCOR_KERNEL_CMB_ISW (object);
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  g_return_if_fail (NC_IS_XCOR_KERNEL_CMB_ISW (object));

  switch (prop_id)
  {
    case PROP_RECOMB:
      g_value_set_object (value, self->recomb);
      break;
    case PROP_NL:
      g_value_set_object (value, self->Nl);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_kernel_cmb_isw_dispose (GObject *object)
{
  NcXcorKernelCMBISW *xcisw              = NC_XCOR_KERNEL_CMB_ISW (object);
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  nc_recomb_clear (&self->recomb);
  ncm_vector_clear (&self->Nl);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_cmb_isw_parent_class)->dispose (object);
}

static void
_nc_xcor_kernel_cmb_isw_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_cmb_isw_parent_class)->finalize (object);
}

static gdouble _nc_xcor_kernel_cmb_isw_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
static gdouble _nc_xcor_kernel_cmb_isw_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
static gdouble _nc_xcor_kernel_cmb_isw_eval_kernel_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble k, gint l);
static gdouble _nc_xcor_kernel_cmb_isw_eval_kernel_prefactor_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
static void _nc_xcor_kernel_cmb_isw_get_k_range_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l, gdouble *kmin, gdouble *kmax);
static void _nc_xcor_kernel_cmb_isw_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo);
static void _nc_xcor_kernel_cmb_isw_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
static guint _nc_xcor_kernel_cmb_isw_obs_len (NcXcorKernel *xclk);
static guint _nc_xcor_kernel_cmb_isw_obs_params_len (NcXcorKernel *xclk);

static void
_nc_xcor_kernel_cmb_isw_constructed (GObject *object)
{
  NcXcorKernelCMBISW *xcisw            = NC_XCOR_KERNEL_CMB_ISW (object);
  NcXcorKernel *xclk                   = NC_XCOR_KERNEL (xcisw);
  NcXcorKernelIntegMethod integ_method = nc_xcor_kernel_get_integ_method (xclk);

  switch (integ_method)
  {
    case NC_XCOR_KERNEL_INTEG_METHOD_LIMBER:
      nc_xcor_kernel_set_eval_kernel_func (xclk,
                                           _nc_xcor_kernel_cmb_isw_eval_kernel_limber,
                                           _nc_xcor_kernel_cmb_isw_eval_kernel_prefactor_limber
                                          );
      nc_xcor_kernel_set_get_k_range_func (xclk, _nc_xcor_kernel_cmb_isw_get_k_range_limber);
      break;
    default:
      g_error ("Unknown integration method %d", integ_method);
      break;
  }

  /* Chain up : middle */
  G_OBJECT_CLASS (nc_xcor_kernel_cmb_isw_parent_class)->constructed (object);
}

static void
nc_xcor_kernel_cmb_isw_class_init (NcXcorKernelCMBISWClass *klass)
{
  GObjectClass *object_class      = G_OBJECT_CLASS (klass);
  NcXcorKernelClass *parent_class = NC_XCOR_KERNEL_CLASS (klass);
  NcmModelClass *model_class      = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_xcor_kernel_cmb_isw_set_property;
  model_class->get_property = &_nc_xcor_kernel_cmb_isw_get_property;
  object_class->constructed = &_nc_xcor_kernel_cmb_isw_constructed;
  object_class->finalize    = &_nc_xcor_kernel_cmb_isw_finalize;
  object_class->dispose     = &_nc_xcor_kernel_cmb_isw_dispose;

  ncm_model_class_set_name_nick (model_class, "Xcor ISW effect", "Xcor-ISW");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorKernelCMBISW:recomb:
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
   * NcXcorKernelCMBISW:Nl:
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

  parent_class->eval_limber_z           = &_nc_xcor_kernel_cmb_isw_eval_limber_z;
  parent_class->eval_limber_z_prefactor = &_nc_xcor_kernel_cmb_isw_eval_limber_z_prefactor;
  parent_class->prepare                 = &_nc_xcor_kernel_cmb_isw_prepare;
  parent_class->add_noise               = &_nc_xcor_kernel_cmb_isw_add_noise;

  parent_class->obs_len        = &_nc_xcor_kernel_cmb_isw_obs_len;
  parent_class->obs_params_len = &_nc_xcor_kernel_cmb_isw_obs_params_len;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

static gdouble
_nc_xcor_kernel_cmb_isw_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  NcmPowspec *ps               = nc_xcor_kernel_peek_powspec (xclk);
  const gdouble k_pivot        = 1.0;
  const gdouble powspec        = ncm_powspec_eval (ps, NCM_MODEL (cosmo), z, k_pivot);
  const gdouble dpowspec_dz    = ncm_powspec_deriv_z (ps, NCM_MODEL (cosmo), z, k_pivot);
  const gdouble d1pz_growth_dz = 1.0 + (1.0 + z) * dpowspec_dz / (2.0 * powspec);
  const gdouble nu             = l + 0.5;
  const gdouble cor_factor     = 1.0 / (nu * nu);

  return cor_factor * xck->E_z * gsl_pow_2 (xck->xi_z) * d1pz_growth_dz;
}

static gdouble
_nc_xcor_kernel_cmb_isw_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  const gdouble Omega_c0 = nc_hicosmo_Omega_c0 (cosmo);
  const gdouble Omega_b0 = nc_hicosmo_Omega_b0 (cosmo);
  const gdouble Omega_m0 = Omega_c0 + Omega_b0;
  const gdouble T_gamma0 = nc_hicosmo_T_gamma0 (cosmo);

  return 3.0 * T_gamma0 * Omega_m0;
}

static gdouble
_nc_xcor_kernel_cmb_isw_eval_kernel_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble k, gint l)
{
  NcXcorKernelCMBISW *xcisw              = NC_XCOR_KERNEL_CMB_ISW (xclk);
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);
  const gdouble nu                       = l + 0.5;
  const gdouble xi_nu                    = nu / k;
  const gdouble z                        = nc_distance_inv_comoving (self->dist, cosmo, xi_nu);
  const gdouble powspec                  = ncm_powspec_eval (self->ps, NCM_MODEL (cosmo), z, k);
  const gdouble dpowspec_dz              = ncm_powspec_deriv_z (self->ps, NCM_MODEL (cosmo), z, k);
  const gdouble d1pz_growth_dz           = 1.0 + (1.0 + z) * dpowspec_dz / (2.0 * powspec);
  const gdouble operator                 = 1.0 / (k * k);

  return sqrt (M_PI / 2.0 / nu) / k * operator * d1pz_growth_dz * sqrt (powspec);
}

static gdouble
_nc_xcor_kernel_cmb_isw_eval_kernel_prefactor_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  const gdouble Omega_c0 = nc_hicosmo_Omega_c0 (cosmo);
  const gdouble Omega_b0 = nc_hicosmo_Omega_b0 (cosmo);
  const gdouble Omega_m0 = Omega_c0 + Omega_b0;
  const gdouble T_gamma0 = nc_hicosmo_T_gamma0 (cosmo);

  return 3.0 * T_gamma0 * Omega_m0;
}

static void
_nc_xcor_kernel_cmb_isw_get_k_range_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l, gdouble *kmin, gdouble *kmax)
{
  NcXcorKernelCMBISW *xcisw              = NC_XCOR_KERNEL_CMB_ISW (xclk);
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);
  NcmPowspec *ps                         = nc_xcor_kernel_peek_powspec (xclk);
  const gdouble ps_kmin                  = ncm_powspec_get_kmin (ps);
  const gdouble ps_kmax                  = ncm_powspec_get_kmax (ps);
  const gdouble nu                       = l + 0.5;
  const gdouble kmin_limber              = nu / self->xi_lss;

  *kmin = GSL_MAX (ps_kmin, kmin_limber);
  *kmax = ps_kmax;
}

static void
_nc_xcor_kernel_cmb_isw_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorKernelCMBISW *xcisw              = NC_XCOR_KERNEL_CMB_ISW (xclk);
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);
  NcDistance *dist                       = nc_xcor_kernel_peek_dist (xclk);
  NcmPowspec *ps                         = nc_xcor_kernel_peek_powspec (xclk);
  const gdouble z_lss                    = nc_distance_decoupling_redshift (dist, cosmo);

  self->dist = dist;
  self->ps   = ps;

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  self->xi_lss = nc_distance_comoving_lss (dist, cosmo);
  nc_xcor_kernel_set_z_range (xclk, 0.0, z_lss, 2.0);
}

static void
_nc_xcor_kernel_cmb_isw_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  NcXcorKernelCMBISW *xcisw              = NC_XCOR_KERNEL_CMB_ISW (xclk);
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  if (self->Nl == NULL)
    g_error ("nc_xcor_kernel_cmb_isw_noise_spec : noise spectrum empty");

  if (lmin + ncm_vector_len (vp1) > self->Nlmax)
    g_error ("nc_xcor_kernel_cmb_isw_noise_spec : too high multipole");

  ncm_vector_memcpy (vp2, vp1);

  {
    NcmVector *Nl_sub = ncm_vector_get_subvector (self->Nl, lmin, ncm_vector_len (vp1));

    ncm_vector_add (vp2, Nl_sub);
    ncm_vector_free (Nl_sub);
  }
  /* return ncm_vector_get (xcisw->Nl, l); */
}

static guint
_nc_xcor_kernel_cmb_isw_obs_len (NcXcorKernel *xclk)
{
  return 1;
}

static guint
_nc_xcor_kernel_cmb_isw_obs_params_len (NcXcorKernel *xclk)
{
  return 0;
}

/**
 * nc_xcor_kernel_cmb_isw_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @recomb: a #NcRecomb
 * @Nl: a #NcmVector
 *
 * Creates a new #NcXcorLimberKernelCMBISW for computing the CMB integrated
 * Sachs-Wolfe (ISW) effect kernel. This kernel describes the correlation between
 * the CMB temperature anisotropies from the ISW effect and large-scale structure.
 *
 * Returns: a new #NcXcorLimberKernelCMBISW
 *
 */
NcXcorKernelCMBISW *
nc_xcor_kernel_cmb_isw_new (NcDistance *dist, NcmPowspec *ps, NcRecomb *recomb, NcmVector *Nl)
{
  NcXcorKernelCMBISW *xcisw = g_object_new (NC_TYPE_XCOR_KERNEL_CMB_ISW,
                                            "dist", dist,
                                            "powspec", ps,
                                            "recomb", recomb,
                                            "Nl", Nl,
                                            NULL);

  return xcisw;
}

