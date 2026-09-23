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
 *    W^{T_\mathrm{ISW}} (z) = \frac{3 \Omega_m H_0^2}{c^2 (\ell + 1/2)^2} \chi^2(z) \frac{\rm{d}}{\rm{d}z}\left((1 + z)D(z)\right).
 * \end{equation}
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "ncm/integration/ncm_integrate.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/spline/ncm_spline2d_bicubic.h"
#include "nc/xcor/nc_xcor_kernel_component.h"
#include "nc/xcor/nc_xcor_kernel_cmb_isw.h"
#include "nc/xcor/nc_xcor.h"
#include "nc_enum_types.h"


#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_bessel.h>
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
  gdouble chi_lss;
  gdouble z_lss;
  NcXcorKernelCMBISWSource source;
  gdouble z_src_min;
  gdouble z_src_max;
  gdouble exp_mtau_min; /* e^{-tau} at the near and far edges of the source range */
  gdouble exp_mtau_max;
  NcXcorKernelComponent *isw_comp;
} NcXcorKernelCMBISWPrivate;

enum
{
  PROP_0,
  PROP_RECOMB,
  PROP_NL,
  PROP_SOURCE,
  PROP_SIZE,
};


G_DEFINE_TYPE_WITH_PRIVATE (NcXcorKernelCMBISW, nc_xcor_kernel_cmb_isw, NC_TYPE_XCOR_KERNEL);

/*
 * ISW Component Definition
 * Uses the component macro to define NcXcorKernelComponentISW
 */

typedef struct _ISWComponentData
{
  NcDistance *dist;
  NcmPowspec *ps;
  NcRecomb *recomb;
  NcXcorKernelCMBISWSource source;
  gdouble z_src_min;
  gdouble z_src_max;
  gdouble exp_mtau_min;
  gdouble exp_mtau_max;
} ISWComponentData;

/*
 * The survival fraction F(z) of the visibility sources, the fraction of photons
 * that last scatter beyond z. Since v_tau = e^{-tau} dtau/dlambda, the
 * cumulative visibility is e^{-tau} itself and F needs no table:
 * F = (e^{-tau(z)} - e^{-tau_max}) / (e^{-tau_min} - e^{-tau_max}), one at the
 * near edge of the source range and zero at the far one. The thin screen is
 * the step at the decoupling redshift, handled through the support.
 */
static gdouble
_nc_xcor_kernel_cmb_isw_survival (NcXcorKernelCMBISWSource source, NcRecomb *recomb, NcHICosmo *cosmo,
                                  gdouble z_src_min, gdouble z_src_max, gdouble exp_mtau_min, gdouble exp_mtau_max, gdouble z)
{
  if (source == NC_XCOR_KERNEL_CMB_ISW_SOURCE_THIN_SCREEN)
    return 1.0;

  if (z <= z_src_min)
    return 1.0;

  if (z >= z_src_max)
    return 0.0;

  return (exp (-nc_recomb_tau (recomb, cosmo, -log1p (z))) - exp_mtau_max) / (exp_mtau_min - exp_mtau_max);
}

/* Helper to get data from component - uses pointer arithmetic to access
 * the data member that comes after the parent_instance in the struct
 * defined by the macro below */
#define _NC_XCOR_KERNEL_COMPONENT_ISW_GET_DATA(comp) \
        ((ISWComponentData *) ((guint8 *) (comp) + sizeof (NcXcorKernelComponent)))


static gdouble _isw_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble chi, gdouble k);
static gdouble _isw_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l);
static void _isw_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *chi_min, gdouble *chi_max, gdouble *k_min, gdouble *k_max);
static void _isw_component_data_clear (ISWComponentData *data);
static NcXcorKernelComponent *_nc_xcor_kernel_component_isw_new (NcDistance *dist, NcmPowspec *ps);

NC_XCOR_KERNEL_COMPONENT_DEFINE_TYPE (NC, XCOR_KERNEL_COMPONENT_ISW,
                                      NcXcorKernelComponentISW,
                                      nc_xcor_kernel_component_isw,
                                      _isw_component_eval_kernel,
                                      _isw_component_eval_prefactor,
                                      _isw_component_get_limits,
                                      ISWComponentData,
                                      _isw_component_data_clear)

static void
nc_xcor_kernel_cmb_isw_init (NcXcorKernelCMBISW *xcisw)
{
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  self->recomb       = NULL;
  self->Nl           = NULL;
  self->Nlmax        = 0;
  self->chi_lss      = 0.0;
  self->z_lss        = 0.0;
  self->source       = NC_XCOR_KERNEL_CMB_ISW_SOURCE_THIN_SCREEN;
  self->z_src_min    = 0.0;
  self->z_src_max    = 0.0;
  self->exp_mtau_min = 1.0;
  self->exp_mtau_max = 0.0;
  self->isw_comp     = NULL;
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
    case PROP_SOURCE:
      nc_xcor_kernel_cmb_isw_set_source (xcisw, g_value_get_enum (value));
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
    case PROP_SOURCE:
      g_value_set_enum (value, self->source);
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
  nc_xcor_kernel_component_clear (&self->isw_comp);

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
static void _nc_xcor_kernel_cmb_isw_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo);
static void _nc_xcor_kernel_cmb_isw_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
static guint _nc_xcor_kernel_cmb_isw_obs_len (NcXcorKernel *xclk);
static guint _nc_xcor_kernel_cmb_isw_obs_params_len (NcXcorKernel *xclk);
static void _nc_xcor_kernel_cmb_isw_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid);
static GPtrArray *_nc_xcor_kernel_cmb_isw_get_component_list (NcXcorKernel *xclk);

static void
_nc_xcor_kernel_cmb_isw_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_kernel_cmb_isw_parent_class)->constructed (object);
  {
    NcXcorKernelCMBISW *xcisw              = NC_XCOR_KERNEL_CMB_ISW (object);
    NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);
    NcDistance *dist                       = nc_xcor_kernel_peek_dist (NC_XCOR_KERNEL (object));
    NcmPowspec *ps                         = nc_xcor_kernel_peek_powspec (NC_XCOR_KERNEL (object));

    g_assert_null (self->isw_comp);
    self->isw_comp = NC_XCOR_KERNEL_COMPONENT (_nc_xcor_kernel_component_isw_new (dist, ps));
  }
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
   * Recombination object used for computing the visibility function.
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
   * Noise power spectrum $N_\ell$.
   */
  g_object_class_install_property (object_class,
                                   PROP_NL,
                                   g_param_spec_object ("Nl",
                                                        NULL,
                                                        "Noise spectrum",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKernelCMBISW:source:
   *
   * Where the CMB photons are placed along the line of sight: a single plane at
   * the decoupling redshift, or the visibility function of
   * #NcXcorKernelCMBISW:recomb, with or without the reionization bump.
   */
  g_object_class_install_property (object_class,
                                   PROP_SOURCE,
                                   g_param_spec_enum ("source",
                                                      NULL,
                                                      "Placement of the CMB sources along the line of sight",
                                                      NC_TYPE_XCOR_KERNEL_CMBISW_SOURCE,
                                                      NC_XCOR_KERNEL_CMB_ISW_SOURCE_THIN_SCREEN,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->get_z_range             = &_nc_xcor_kernel_cmb_isw_get_z_range;
  parent_class->eval_limber_z           = &_nc_xcor_kernel_cmb_isw_eval_limber_z;
  parent_class->eval_limber_z_prefactor = &_nc_xcor_kernel_cmb_isw_eval_limber_z_prefactor;
  parent_class->prepare                 = &_nc_xcor_kernel_cmb_isw_prepare;
  parent_class->add_noise               = &_nc_xcor_kernel_cmb_isw_add_noise;
  parent_class->obs_len                 = &_nc_xcor_kernel_cmb_isw_obs_len;
  parent_class->obs_params_len          = &_nc_xcor_kernel_cmb_isw_obs_params_len;
  parent_class->get_component_list      = &_nc_xcor_kernel_cmb_isw_get_component_list;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

/*
 * Implementation of the old limber_z interface.
 */

static void
_nc_xcor_kernel_cmb_isw_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid)
{
  NcXcorKernelCMBISW *xcisw              = NC_XCOR_KERNEL_CMB_ISW (xclk);
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  *zmin = 0.0;
  *zmax = (self->source == NC_XCOR_KERNEL_CMB_ISW_SOURCE_THIN_SCREEN) ? self->z_lss : self->z_src_max;
  *zmid = 2.0;
}

static gdouble
_nc_xcor_kernel_cmb_isw_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  NcXcorKernelCMBISW *xcisw              = NC_XCOR_KERNEL_CMB_ISW (xclk);
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);
  NcmPowspec *ps                         = nc_xcor_kernel_peek_powspec (xclk);
  const gdouble k_pivot                  = 1.0;
  const gdouble powspec                  = ncm_powspec_eval (ps, NCM_MODEL (cosmo), z, k_pivot);
  const gdouble dpowspec_dz              = ncm_powspec_deriv_z (ps, NCM_MODEL (cosmo), z, k_pivot);
  const gdouble d1pz_growth_dz           = 1.0 + (1.0 + z) * dpowspec_dz / (2.0 * powspec);
  const gdouble F_z                      = _nc_xcor_kernel_cmb_isw_survival (self->source, self->recomb, cosmo, self->z_src_min, self->z_src_max,
                                                                             self->exp_mtau_min, self->exp_mtau_max, z);

  return xck->E_z * gsl_pow_2 (xck->chi_z) * d1pz_growth_dz * F_z;
}

static gdouble
_nc_xcor_kernel_cmb_isw_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  const gdouble nu       = l + 0.5;
  const gdouble Omega_c0 = nc_hicosmo_Omega_c0 (cosmo);
  const gdouble Omega_b0 = nc_hicosmo_Omega_b0 (cosmo);
  const gdouble Omega_m0 = Omega_c0 + Omega_b0;
  const gdouble T_gamma0 = nc_hicosmo_T_gamma0 (cosmo);

  return 3.0 * T_gamma0 * Omega_m0 / (nu * nu);
}

/*
 * Implementation of the new Component interface.
 *
 * The ISW component corresponds to the full kernel in this case, since there are no separate contributions.
 */

static void
_isw_component_data_clear (ISWComponentData *data)
{
  /* No need to clear, these are weak references from parent kernel */
}

static gdouble
_isw_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble chi, gdouble k)
{
  /* Access data using offset from base class */
  ISWComponentData *data       = _NC_XCOR_KERNEL_COMPONENT_ISW_GET_DATA (comp);
  const gdouble z              = nc_distance_inv_comoving (data->dist, cosmo, chi);
  const gdouble E_z            = nc_hicosmo_E (cosmo, z);
  const gdouble powspec        = ncm_powspec_eval (data->ps, NCM_MODEL (cosmo), z, k / nc_hicosmo_RH_Mpc (cosmo));
  const gdouble dpowspec_dz    = ncm_powspec_deriv_z (data->ps, NCM_MODEL (cosmo), z, k / nc_hicosmo_RH_Mpc (cosmo));
  const gdouble d1pz_growth_dz = 1.0 + (1.0 + z) * dpowspec_dz / (2.0 * powspec);
  const gdouble operator       = 1.0 / gsl_pow_2 (k);
  const gdouble F_z            = _nc_xcor_kernel_cmb_isw_survival (data->source, data->recomb, cosmo, data->z_src_min, data->z_src_max,
                                                                   data->exp_mtau_min, data->exp_mtau_max, z);

  return operator * E_z * d1pz_growth_dz * F_z * sqrt (powspec);
}

static gdouble
_isw_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l)
{
  const gdouble Omega_c0 = nc_hicosmo_Omega_c0 (cosmo);
  const gdouble Omega_b0 = nc_hicosmo_Omega_b0 (cosmo);
  const gdouble Omega_m0 = Omega_c0 + Omega_b0;
  const gdouble T_gamma0 = nc_hicosmo_T_gamma0 (cosmo);

  return 3.0 * T_gamma0 * Omega_m0;
}

static void
_isw_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *chi_min, gdouble *chi_max, gdouble *k_min, gdouble *k_max)
{
  ISWComponentData *data = _NC_XCOR_KERNEL_COMPONENT_ISW_GET_DATA (comp);
  NcDistance *dist       = data->dist;
  NcmPowspec *ps         = data->ps;

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  {
    const gdouble chi_lss = (data->source == NC_XCOR_KERNEL_CMB_ISW_SOURCE_THIN_SCREEN) ?
                            nc_distance_comoving_lss (dist, cosmo) :
                            nc_distance_comoving (dist, cosmo, data->z_src_max);

    *chi_min = nc_distance_comoving (dist, cosmo, 1.0e-6);
    *chi_max = chi_lss;
    *k_min   = ncm_powspec_get_kmin (ps) * nc_hicosmo_RH_Mpc (cosmo);
    *k_max   = ncm_powspec_get_kmax (ps) * nc_hicosmo_RH_Mpc (cosmo);
  }
}

static NcXcorKernelComponent *
_nc_xcor_kernel_component_isw_new (NcDistance *dist, NcmPowspec *ps)
{
  NcXcorKernelComponent *comp = g_object_new (nc_xcor_kernel_component_isw_get_type (), NULL);
  ISWComponentData *data      = _NC_XCOR_KERNEL_COMPONENT_ISW_GET_DATA (comp);

  data->dist         = dist;
  data->ps           = ps;
  data->recomb       = NULL;
  data->source       = NC_XCOR_KERNEL_CMB_ISW_SOURCE_THIN_SCREEN;
  data->z_src_min    = 0.0;
  data->z_src_max    = 0.0;
  data->exp_mtau_min = 1.0;
  data->exp_mtau_max = 0.0;

  return comp;
}

/*
 * The visibility sources. The far edge of the source range is where the
 * visibility has dropped to 1e-4 of its recombination peak. The near edge is
 * z = 0 with reionization and, without it, the redshift between the
 * reionization bump and the shell where the visibility is smallest: starting
 * the source where nothing scatters keeps the survival fraction smooth to the
 * tolerances the forcing is fitted to.
 */
#define NC_XCOR_KERNEL_CMB_ISW_VISIBILITY_LOGREF (4.0 * M_LN10)

static void
_nc_xcor_kernel_cmb_isw_prepare_visibility (NcXcorKernelCMBISW *xcisw, NcHICosmo *cosmo)
{
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);
  gdouble lambda_max, lambda_l, lambda_u;

  if (self->recomb == NULL)
    g_error ("nc_xcor_kernel_cmb_isw_prepare: the visibility source needs the recomb property set.");

  nc_recomb_prepare_if_needed (self->recomb, cosmo);
  nc_recomb_v_tau_lambda_features (self->recomb, cosmo, NC_XCOR_KERNEL_CMB_ISW_VISIBILITY_LOGREF,
                                   &lambda_max, &lambda_l, &lambda_u);

  /*
   * lambda decreases with z: lambda_u is the near edge of the shell, lambda_l
   * the far one. With reionization the sources extend to z = 0; without it they
   * stop at the visibility minimum between the shell and the reionization bump,
   * where nothing scatters, so that no photon is cut in the middle of a source.
   */
  if (self->source == NC_XCOR_KERNEL_CMB_ISW_SOURCE_VISIBILITY_REIONIZATION)
    lambda_u = 0.0;
  else
    lambda_u = nc_recomb_get_v_tau_reion_min_lambda (self->recomb, cosmo);

  self->z_src_min    = expm1 (-lambda_u);
  self->z_src_max    = expm1 (-lambda_l);
  self->exp_mtau_min = exp (-nc_recomb_tau (self->recomb, cosmo, lambda_u));
  self->exp_mtau_max = exp (-nc_recomb_tau (self->recomb, cosmo, lambda_l));
}

/*
 * Kernel implementation.
 */

static void
_nc_xcor_kernel_cmb_isw_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorKernelCMBISW *xcisw              = NC_XCOR_KERNEL_CMB_ISW (xclk);
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);
  NcDistance *dist                       = nc_xcor_kernel_peek_dist (xclk);
  NcmPowspec *ps                         = nc_xcor_kernel_peek_powspec (xclk);
  const gdouble z_lss                    = nc_distance_decoupling_redshift (dist, cosmo);

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  self->chi_lss = nc_distance_comoving_lss (dist, cosmo);
  self->z_lss   = z_lss;

  if (self->source != NC_XCOR_KERNEL_CMB_ISW_SOURCE_THIN_SCREEN)
    _nc_xcor_kernel_cmb_isw_prepare_visibility (xcisw, cosmo);

  {
    ISWComponentData *data = _NC_XCOR_KERNEL_COMPONENT_ISW_GET_DATA (self->isw_comp);

    data->recomb       = self->recomb;
    data->source       = self->source;
    data->z_src_min    = self->z_src_min;
    data->z_src_max    = self->z_src_max;
    data->exp_mtau_min = self->exp_mtau_min;
    data->exp_mtau_max = self->exp_mtau_max;
  }

  g_assert_nonnull (self->isw_comp);
  nc_xcor_kernel_component_prepare (self->isw_comp, cosmo);
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

static GPtrArray *
_nc_xcor_kernel_cmb_isw_get_component_list (NcXcorKernel *xclk)
{
  NcXcorKernelCMBISW *xcisw              = NC_XCOR_KERNEL_CMB_ISW (xclk);
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);
  GPtrArray *comp_list                   = g_ptr_array_new_with_free_func (g_object_unref);

  if (self->isw_comp != NULL)
    g_ptr_array_add (comp_list, g_object_ref (self->isw_comp));

  return comp_list;
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

/**
 * nc_xcor_kernel_cmb_isw_eval_k_max:
 * @xcisw: a #NcXcorKernelCMBISW
 * @x: the x value (x = k * chi)
 *
 * Evaluates k_max at the given x value from kernel analysis.
 *
 * Returns: the k_max value at x
 */
gdouble
nc_xcor_kernel_cmb_isw_eval_k_max (NcXcorKernelCMBISW *xcisw, gdouble x)
{
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  g_assert (self->isw_comp != NULL);

  return nc_xcor_kernel_component_eval_k_max (self->isw_comp, x);
}

/**
 * nc_xcor_kernel_cmb_isw_eval_KL_max:
 * @xcisw: a #NcXcorKernelCMBISW
 * @x: the x value (x = k * chi)
 *
 * Evaluates KL_max at the given x value from kernel analysis using the Limber approximation.
 *
 * Returns: the KL_max value at x
 */
gdouble
nc_xcor_kernel_cmb_isw_eval_KL_max (NcXcorKernelCMBISW *xcisw, gdouble x)
{
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  g_assert (self->isw_comp != NULL);

  return nc_xcor_kernel_component_eval_KL_max (self->isw_comp, x);
}

/**
 * nc_xcor_kernel_cmb_isw_eval_k_epsilon:
 * @xcisw: a #NcXcorKernelCMBISW
 * @x: the x value (x = k * chi)
 *
 * Evaluates k_epsilon at the given x value from kernel analysis.
 *
 * Returns: the k_epsilon value at x
 */
gdouble
nc_xcor_kernel_cmb_isw_eval_k_epsilon (NcXcorKernelCMBISW *xcisw, gdouble x)
{
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  g_assert (self->isw_comp != NULL);

  return nc_xcor_kernel_component_eval_k_epsilon (self->isw_comp, x);
}

/**
 * nc_xcor_kernel_cmb_isw_set_epsilon:
 * @xcisw: a #NcXcorKernelCMBISW
 * @epsilon: the epsilon value for kernel analysis
 *
 * Sets the epsilon value used in kernel analysis to determine where the kernel
 * drops to epsilon * K_max. Default value is 1.0e-3.
 */
void
nc_xcor_kernel_cmb_isw_set_epsilon (NcXcorKernelCMBISW *xcisw, gdouble epsilon)
{
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  if (self->isw_comp != NULL)
    nc_xcor_kernel_component_set_epsilon (self->isw_comp, epsilon);
}

/**
 * nc_xcor_kernel_cmb_isw_get_epsilon:
 * @xcisw: a #NcXcorKernelCMBISW
 *
 * Gets the epsilon value used in kernel analysis.
 *
 * Returns: the epsilon value
 */
gdouble
nc_xcor_kernel_cmb_isw_get_epsilon (NcXcorKernelCMBISW *xcisw)
{
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  if (self->isw_comp == NULL)
    return NC_XCOR_KERNEL_COMPONENT_DEFAULT_EPSILON;

  return nc_xcor_kernel_component_get_epsilon (self->isw_comp);
}

/**
 * nc_xcor_kernel_cmb_isw_set_source:
 * @xcisw: a #NcXcorKernelCMBISW
 * @source: a #NcXcorKernelCMBISWSource
 *
 * Sets where the CMB photons are placed along the line of sight. The visibility
 * sources require the #NcXcorKernelCMBISW:recomb object. The kernel is marked
 * outdated and is prepared again on the next use.
 */
void
nc_xcor_kernel_cmb_isw_set_source (NcXcorKernelCMBISW *xcisw, NcXcorKernelCMBISWSource source)
{
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  if (self->source != source)
  {
    self->source = source;
    nc_xcor_kernel_mark_outdated (NC_XCOR_KERNEL (xcisw));
  }
}

/**
 * nc_xcor_kernel_cmb_isw_get_source:
 * @xcisw: a #NcXcorKernelCMBISW
 *
 * Returns: the #NcXcorKernelCMBISWSource in use.
 */
NcXcorKernelCMBISWSource
nc_xcor_kernel_cmb_isw_get_source (NcXcorKernelCMBISW *xcisw)
{
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  return self->source;
}

