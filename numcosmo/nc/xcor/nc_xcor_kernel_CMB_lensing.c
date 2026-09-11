/***************************************************************************
 *            nc_xcor_kernel_cmb_lensing.c
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 *  Sat December 27 20:21:01 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
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
 * NcXcorKernelCMBLensing:
 *
 * Implementation of #NcXcorKernel for CMB lensing
 *
 * With the thin-screen source, #NC_XCOR_KERNEL_CMB_LENSING_SOURCE_THIN_SCREEN,
 * the kernel is
 * \begin{equation}
 *    W^{\kappa_\mathrm{CMB}} (z) = \frac{3}{2} \frac{\Omega_m H_0^2}{c} \frac{(1+z)}{H(z)} \chi(z) \frac{\chi(z_*) - \chi(z)}{\chi(z_*)},
 * \end{equation}
 * every photon having last scattered at the decoupling redshift $z_*$ of the
 * #NcDistance. With #NC_XCOR_KERNEL_CMB_LENSING_SOURCE_VISIBILITY the single
 * plane is replaced by the recombination visibility function $g(z)$ of the
 * #NcXcorKernelCMBLensing:recomb object, normalized to unit integral over the
 * last-scattering shell, and the geometric factor becomes the lensing efficiency
 * \begin{equation}
 *    q(z) = \int_z^{z_{\max}} \mathrm{d}z'\, g(z') \frac{\chi(z') - \chi(z)}{\chi(z')},
 * \end{equation}
 * computed by #NcXcorLensingEfficiency. The two agree to the width of the shell
 * over $\chi_*$; the difference is at the far end of the support, where the
 * thin screen ends with a kink and the visibility source goes to zero smoothly.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "nc/xcor/nc_xcor_kernel_CMB_lensing.h"
#include "nc/xcor/nc_xcor_kernel_component.h"
#include "nc/xcor/nc_xcor_lensing_efficiency.h"
#include "nc/xcor/nc_xcor.h"
#include "nc_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#endif /* NUMCOSMO_GIR_SCAN */


struct _NcXcorKernelCMBLensing
{
  /*< private >*/
  NcXcorKernel parent_instance;

  NcRecomb *recomb;

  NcmVector *Nl;
  guint Nlmax;

  gdouble z_lss;
  gdouble chi_lss;
  gdouble dt_lss;
  gdouble dt;

  NcXcorKernelCMBLensingSource source;
  NcXcorLensingEfficiency *lens_eff;
  NcHICosmo *cosmo_prep; /* valid while prepare() runs; read by the visibility source */
  gdouble z_src_min;
  gdouble z_src_max;
  gdouble src_norm;

  NcDistance *dist;
  NcmPowspec *ps;
  NcXcorKernelComponent *cmb_lens_comp;
};

enum
{
  PROP_0,
  PROP_RECOMB,
  PROP_NL,
  PROP_SOURCE,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelCMBLensing, nc_xcor_kernel_cmb_lensing, NC_TYPE_XCOR_KERNEL)

/*
 * The visibility source: W_src(z) = g(z) / (1 + z) / norm, with g the visibility
 * of the recomb object in lambda = -ln(1 + z), restricted to the last-scattering
 * shell found by nc_recomb_v_tau_lambda_features() and normalized to unit
 * integral over it. The lensing efficiency it feeds is the same object the weak
 * lensing kernel uses for its dn/dz.
 */
#define NC_XCOR_KERNEL_CMB_LENSING_VISIBILITY_LOGREF (4.0 * M_LN10)

static gdouble _nc_xcor_kernel_cmb_lensing_lens_eff_eval_source (NcXcorLensingEfficiency *lens_eff, gdouble z);
static void _nc_xcor_kernel_cmb_lensing_lens_eff_get_z_range (NcXcorLensingEfficiency *lens_eff, gdouble *zmin, gdouble *zmax);

NC_XCOR_LENSING_EFFICIENCY_DEFINE_TYPE (NC, XCOR_KERNEL_CMB_LENSING_LENS_EFF,
                                        NcXcorKernelCMBLensingLensEff,
                                        nc_xcor_kernel_cmb_lensing_lens_eff,
                                        _nc_xcor_kernel_cmb_lensing_lens_eff_eval_source,
                                        _nc_xcor_kernel_cmb_lensing_lens_eff_get_z_range,
                                        NcXcorKernelCMBLensing *)

/*
 * CMB Lensing Component Definition
 * Handles the CMB lensing convergence kernel
 */

typedef struct _CMBLensingComponentData
{
  NcDistance *dist;
  NcmPowspec *ps;
  gdouble z_lss;
  gdouble dt_lss;
  NcXcorKernelCMBLensingSource source;
  NcXcorLensingEfficiency *lens_eff;
  gdouble z_src_max;
} CMBLensingComponentData;

#define _NC_XCOR_KERNEL_COMPONENT_CMB_LENSING_GET_DATA(comp) \
        ((CMBLensingComponentData *) ((guint8 *) (comp) + sizeof (NcXcorKernelComponent)))

static gdouble _cmb_lensing_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble chi, gdouble k);
static gdouble _cmb_lensing_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l);
static void _cmb_lensing_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *chi_min, gdouble *chi_max, gdouble *k_min, gdouble *k_max);
static void _cmb_lensing_component_data_clear (CMBLensingComponentData *data);
static NcXcorKernelComponent *_nc_xcor_kernel_component_cmb_lensing_new (NcDistance *dist, NcmPowspec *ps);

NC_XCOR_KERNEL_COMPONENT_DEFINE_TYPE (NC, XCOR_KERNEL_COMPONENT_CMB_LENSING,
                                      NcXcorKernelComponentCMBLensing,
                                      nc_xcor_kernel_component_cmb_lensing,
                                      _cmb_lensing_component_eval_kernel,
                                      _cmb_lensing_component_eval_prefactor,
                                      _cmb_lensing_component_get_limits,
                                      CMBLensingComponentData,
                                      _cmb_lensing_component_data_clear)

static void
nc_xcor_kernel_cmb_lensing_init (NcXcorKernelCMBLensing *xclkl)
{
  xclkl->recomb = NULL;

  xclkl->Nl    = NULL;
  xclkl->Nlmax = 0;

  xclkl->z_lss         = 0.0;
  xclkl->chi_lss       = 0.0;
  xclkl->dt_lss        = 0.0;
  xclkl->source        = NC_XCOR_KERNEL_CMB_LENSING_SOURCE_THIN_SCREEN;
  xclkl->lens_eff      = NULL;
  xclkl->cosmo_prep    = NULL;
  xclkl->z_src_min     = 0.0;
  xclkl->z_src_max     = 0.0;
  xclkl->src_norm      = 1.0;
  xclkl->dt            = 0.0;
  xclkl->dist          = NULL;
  xclkl->ps            = NULL;
  xclkl->cmb_lens_comp = NULL;
}

static void
_nc_xcor_kernel_cmb_lensing_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelCMBLensing *xclkl = NC_XCOR_KERNEL_CMB_LENSING (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_CMB_LENSING (object));

  switch (prop_id)
  {
    case PROP_RECOMB:
      xclkl->recomb = g_value_dup_object (value);
      break;
    case PROP_NL:
      xclkl->Nl    = g_value_dup_object (value);
      xclkl->Nlmax = ncm_vector_len (xclkl->Nl) - 1;
      break;
    case PROP_SOURCE:
      nc_xcor_kernel_cmb_lensing_set_source (xclkl, g_value_get_enum (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_kernel_cmb_lensing_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelCMBLensing *xclkl = NC_XCOR_KERNEL_CMB_LENSING (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_CMB_LENSING (object));

  switch (prop_id)
  {
    case PROP_RECOMB:
      g_value_set_object (value, xclkl->recomb);
      break;
    case PROP_NL:
      g_value_set_object (value, xclkl->Nl);
      break;
    case PROP_SOURCE:
      g_value_set_enum (value, xclkl->source);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_kernel_cmb_lensing_dispose (GObject *object)
{
  NcXcorKernelCMBLensing *xclkl = NC_XCOR_KERNEL_CMB_LENSING (object);

  nc_recomb_clear (&xclkl->recomb);
  ncm_vector_clear (&xclkl->Nl);
  nc_xcor_kernel_component_clear (&xclkl->cmb_lens_comp);
  nc_xcor_lensing_efficiency_clear (&xclkl->lens_eff);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_cmb_lensing_parent_class)->dispose (object);
}

static void
_nc_xcor_kernel_cmb_lensing_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_cmb_lensing_parent_class)->finalize (object);
}

static gdouble _nc_xcor_kernel_cmb_lensing_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
static gdouble _nc_xcor_kernel_cmb_lensing_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
static void _nc_xcor_kernel_cmb_lensing_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo);
static void _nc_xcor_kernel_cmb_lensing_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
static guint _nc_xcor_kernel_cmb_lensing_obs_len (NcXcorKernel *xclk);
static guint _nc_xcor_kernel_cmb_lensing_obs_params_len (NcXcorKernel *xclk);
static void _nc_xcor_kernel_cmb_lensing_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid);
static GPtrArray *_nc_xcor_kernel_cmb_lensing_get_component_list (NcXcorKernel *xclk);

static void
_nc_xcor_kernel_cmb_lensing_constructed (GObject *object)
{
  /* Chain up to parent constructed */
  G_OBJECT_CLASS (nc_xcor_kernel_cmb_lensing_parent_class)->constructed (object);
  {
    NcXcorKernelCMBLensing *xclkl = NC_XCOR_KERNEL_CMB_LENSING (object);
    NcDistance *dist              = nc_xcor_kernel_peek_dist (NC_XCOR_KERNEL (object));
    NcmPowspec *ps                = nc_xcor_kernel_peek_powspec (NC_XCOR_KERNEL (object));

    g_assert_null (xclkl->cmb_lens_comp);
    xclkl->cmb_lens_comp = _nc_xcor_kernel_component_cmb_lensing_new (dist, ps);

    {
      NcXcorKernelCMBLensingLensEff *lens_eff_obj = g_object_new (
        nc_xcor_kernel_cmb_lensing_lens_eff_get_type (),
        "distance", dist,
        NULL);

      lens_eff_obj->data = xclkl;
      xclkl->lens_eff    = NC_XCOR_LENSING_EFFICIENCY (lens_eff_obj);
    }
  }
}

static void
nc_xcor_kernel_cmb_lensing_class_init (NcXcorKernelCMBLensingClass *klass)
{
  GObjectClass *object_class      = G_OBJECT_CLASS (klass);
  NcXcorKernelClass *parent_class = NC_XCOR_KERNEL_CLASS (klass);
  NcmModelClass *model_class      = NCM_MODEL_CLASS (klass);

  object_class->constructed = &_nc_xcor_kernel_cmb_lensing_constructed;
  object_class->finalize    = &_nc_xcor_kernel_cmb_lensing_finalize;
  object_class->dispose     = &_nc_xcor_kernel_cmb_lensing_dispose;
  model_class->set_property = &_nc_xcor_kernel_cmb_lensing_set_property;
  model_class->get_property = &_nc_xcor_kernel_cmb_lensing_get_property;

  ncm_model_class_set_name_nick (model_class, "Xcor lensing distribution", "Xcor-lensing");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorKernelCMBLensing:recomb:
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
   * NcXcorKernelCMBLensing:Nl:
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
   * NcXcorKernelCMBLensing:source:
   *
   * Where the CMB photons are placed along the line of sight: a single plane at
   * the decoupling redshift, or the recombination visibility function of
   * #NcXcorKernelCMBLensing:recomb.
   */
  g_object_class_install_property (object_class,
                                   PROP_SOURCE,
                                   g_param_spec_enum ("source",
                                                      NULL,
                                                      "Placement of the CMB sources along the line of sight",
                                                      NC_TYPE_XCOR_KERNEL_CMB_LENSING_SOURCE,
                                                      NC_XCOR_KERNEL_CMB_LENSING_SOURCE_THIN_SCREEN,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->eval_limber_z           = &_nc_xcor_kernel_cmb_lensing_eval_limber_z;
  parent_class->eval_limber_z_prefactor = &_nc_xcor_kernel_cmb_lensing_eval_limber_z_prefactor;
  parent_class->prepare                 = &_nc_xcor_kernel_cmb_lensing_prepare;
  parent_class->add_noise               = &_nc_xcor_kernel_cmb_lensing_add_noise;

  parent_class->obs_len            = &_nc_xcor_kernel_cmb_lensing_obs_len;
  parent_class->obs_params_len     = &_nc_xcor_kernel_cmb_lensing_obs_params_len;
  parent_class->get_z_range        = &_nc_xcor_kernel_cmb_lensing_get_z_range;
  parent_class->get_component_list = &_nc_xcor_kernel_cmb_lensing_get_component_list;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

static gdouble
_nc_xcor_kernel_cmb_lensing_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l) /*, gdouble geo_z[]) */
{
  NcXcorKernelCMBLensing *xclkl = NC_XCOR_KERNEL_CMB_LENSING (xclk);
  NcDistance *dist              = nc_xcor_kernel_peek_dist (xclk);
  const gdouble dt              = nc_distance_transverse (dist, cosmo, z);

  if (xclkl->source != NC_XCOR_KERNEL_CMB_LENSING_SOURCE_THIN_SCREEN)
  {
    return xck->chi_z * xck->chi_z * (1.0 + z) * nc_xcor_lensing_efficiency_eval (xclkl->lens_eff, z) / dt;
  }
  else
  {
    const gdouble dt_z_zlss = nc_distance_transverse_z1_z2 (dist, cosmo, z, xclkl->z_lss);

    return xck->chi_z * xck->chi_z * (1.0 + z) * dt_z_zlss / (xclkl->dt_lss * dt);
  }
}

static gdouble
_nc_xcor_kernel_cmb_lensing_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  const gdouble nu      = l + 0.5;
  const gdouble lfactor = l * (l + 1.0);

  return 1.5 * nc_hicosmo_Omega_m0 (cosmo) * lfactor / (nu * nu);
}

/*
 * Implementation of the new Component interface.
 */

static void
_cmb_lensing_component_data_clear (CMBLensingComponentData *data)
{
  /* No need to clear, these are weak references from parent kernel */
}

static gdouble
_cmb_lensing_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble chi, gdouble k)
{
  CMBLensingComponentData *data = _NC_XCOR_KERNEL_COMPONENT_CMB_LENSING_GET_DATA (comp);
  const gdouble z               = nc_distance_inv_comoving (data->dist, cosmo, chi);
  const gdouble powspec         = ncm_powspec_eval (data->ps, NCM_MODEL (cosmo), z, k / nc_hicosmo_RH_Mpc (cosmo));
  const gdouble dt              = nc_distance_transverse (data->dist, cosmo, z);
  const gdouble operator_k      = 1.0 / gsl_pow_2 (k);

  if (data->source != NC_XCOR_KERNEL_CMB_LENSING_SOURCE_THIN_SCREEN)
  {
    return operator_k * (1.0 + z) * nc_xcor_lensing_efficiency_eval (data->lens_eff, z) / dt * sqrt (powspec);
  }
  else
  {
    const gdouble dt_z_zlss = nc_distance_transverse_z1_z2 (data->dist, cosmo, z, data->z_lss);

    return operator_k * (1.0 + z) * dt_z_zlss / (data->dt_lss * dt) * sqrt (powspec);
  }
}

static gdouble
_cmb_lensing_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l)
{
  const gdouble lfactor = l * (l + 1.0);

  return 1.5 * nc_hicosmo_Omega_m0 (cosmo) * lfactor;
}

static void
_cmb_lensing_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *chi_min, gdouble *chi_max, gdouble *k_min, gdouble *k_max)
{
  CMBLensingComponentData *data = _NC_XCOR_KERNEL_COMPONENT_CMB_LENSING_GET_DATA (comp);
  NcDistance *dist              = data->dist;
  NcmPowspec *ps                = data->ps;

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  {
    const gdouble chi_src = (data->source != NC_XCOR_KERNEL_CMB_LENSING_SOURCE_THIN_SCREEN) ?
                            nc_distance_comoving (dist, cosmo, data->z_src_max) :
                            nc_distance_comoving_lss (dist, cosmo);

    *chi_min = nc_distance_comoving (dist, cosmo, 1.0e-6);
    *chi_max = chi_src * (1.0 - 1.0e-6);
    *k_min   = ncm_powspec_get_kmin (ps) * nc_hicosmo_RH_Mpc (cosmo);
    *k_max   = ncm_powspec_get_kmax (ps) * nc_hicosmo_RH_Mpc (cosmo);
  }
}

static NcXcorKernelComponent *
_nc_xcor_kernel_component_cmb_lensing_new (NcDistance *dist, NcmPowspec *ps)
{
  NcXcorKernelComponent *comp   = g_object_new (nc_xcor_kernel_component_cmb_lensing_get_type (), NULL);
  CMBLensingComponentData *data = _NC_XCOR_KERNEL_COMPONENT_CMB_LENSING_GET_DATA (comp);

  data->dist      = dist;
  data->ps        = ps;
  data->z_lss     = 0.0;
  data->dt_lss    = 0.0;
  data->source    = NC_XCOR_KERNEL_CMB_LENSING_SOURCE_THIN_SCREEN;
  data->lens_eff  = NULL;
  data->z_src_max = 0.0;

  return comp;
}

/*
 * The visibility source and its normalization.
 */

static gdouble
_nc_xcor_kernel_cmb_lensing_lens_eff_eval_source (NcXcorLensingEfficiency *lens_eff, gdouble z)
{
  NcXcorKernelCMBLensingLensEff *data_obj = NC_XCOR_KERNEL_CMB_LENSING_LENS_EFF (lens_eff);
  NcXcorKernelCMBLensing *xclkl           = data_obj->data;

  if ((z < xclkl->z_src_min) || (z > xclkl->z_src_max))
    return 0.0;

  /* g(z) dz = v_tau dlambda with lambda = -ln(1 + z), so |dlambda/dz| = 1/(1 + z). */
  return nc_recomb_v_tau (xclkl->recomb, xclkl->cosmo_prep, -log1p (z)) / ((1.0 + z) * xclkl->src_norm);
}

static void
_nc_xcor_kernel_cmb_lensing_lens_eff_get_z_range (NcXcorLensingEfficiency *lens_eff, gdouble *zmin, gdouble *zmax)
{
  NcXcorKernelCMBLensingLensEff *data_obj = NC_XCOR_KERNEL_CMB_LENSING_LENS_EFF (lens_eff);
  NcXcorKernelCMBLensing *xclkl           = data_obj->data;

  *zmin = xclkl->z_src_min;
  *zmax = xclkl->z_src_max;
}

static gdouble
_nc_xcor_kernel_cmb_lensing_src_integrand (gdouble z, gpointer user_data)
{
  NcXcorKernelCMBLensing *xclkl = user_data;

  return nc_recomb_v_tau (xclkl->recomb, xclkl->cosmo_prep, -log1p (z)) / (1.0 + z);
}

static void
_nc_xcor_kernel_cmb_lensing_prepare_visibility (NcXcorKernelCMBLensing *xclkl, NcHICosmo *cosmo)
{
  gdouble lambda_max, lambda_l, lambda_u;

  if (xclkl->recomb == NULL)
    g_error ("nc_xcor_kernel_cmb_lensing_prepare: the visibility source needs the recomb property set.");

  nc_recomb_prepare_if_needed (xclkl->recomb, cosmo);
  nc_recomb_v_tau_lambda_features (xclkl->recomb, cosmo, NC_XCOR_KERNEL_CMB_LENSING_VISIBILITY_LOGREF,
                                   &lambda_max, &lambda_l, &lambda_u);

  /*
   * lambda decreases with z: lambda_u is the near edge of the shell, lambda_l
   * the far one. With reionization the source starts at z = 0 and the low-z
   * bump of the visibility is part of it.
   */
  xclkl->cosmo_prep = cosmo;
  xclkl->z_src_min  = (xclkl->source == NC_XCOR_KERNEL_CMB_LENSING_SOURCE_VISIBILITY_REIONIZATION) ? 0.0 : expm1 (-lambda_u);
  xclkl->z_src_max  = expm1 (-lambda_l);
  xclkl->src_norm   = 1.0;

  {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
    gsl_function F;
    gdouble norm, err;

    F.function = &_nc_xcor_kernel_cmb_lensing_src_integrand;
    F.params   = xclkl;

    gsl_integration_qag (&F, xclkl->z_src_min, xclkl->z_src_max, 0.0, 1.0e-11, 1000, GSL_INTEG_GAUSS61, w, &norm, &err);
    gsl_integration_workspace_free (w);

    xclkl->src_norm = norm;
  }

  nc_xcor_lensing_efficiency_prepare (xclkl->lens_eff, cosmo);
}

/*
 * Kernel implementation.
 */

static void
_nc_xcor_kernel_cmb_lensing_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorKernelCMBLensing *xclkl = NC_XCOR_KERNEL_CMB_LENSING (xclk);
  NcDistance *dist              = nc_xcor_kernel_peek_dist (xclk);
  NcmPowspec *ps                = nc_xcor_kernel_peek_powspec (xclk);

  xclkl->dist = dist;
  xclkl->ps   = ps;

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  xclkl->z_lss   = nc_distance_decoupling_redshift (dist, cosmo);
  xclkl->chi_lss = nc_distance_comoving_lss (dist, cosmo);
  xclkl->dt_lss  = nc_distance_transverse (dist, cosmo, xclkl->z_lss);

  if (xclkl->source != NC_XCOR_KERNEL_CMB_LENSING_SOURCE_THIN_SCREEN)
    _nc_xcor_kernel_cmb_lensing_prepare_visibility (xclkl, cosmo);

  /* Update component data with computed values */
  {
    CMBLensingComponentData *data = _NC_XCOR_KERNEL_COMPONENT_CMB_LENSING_GET_DATA (xclkl->cmb_lens_comp);

    data->z_lss     = xclkl->z_lss;
    data->dt_lss    = xclkl->dt_lss;
    data->source    = xclkl->source;
    data->lens_eff  = xclkl->lens_eff;
    data->z_src_max = xclkl->z_src_max;
  }

  g_assert_nonnull (xclkl->cmb_lens_comp);
  nc_xcor_kernel_component_prepare (xclkl->cmb_lens_comp, cosmo);
}

static void
_nc_xcor_kernel_cmb_lensing_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  NcXcorKernelCMBLensing *xclkl = NC_XCOR_KERNEL_CMB_LENSING (xclk);

  if (xclkl->Nl == NULL)
    g_error ("nc_xcor_kernel_cmb_lensing_noise_spec : noise spectrum empty");

  if (lmin + ncm_vector_len (vp1) > xclkl->Nlmax)
    g_error ("nc_xcor_kernel_cmb_lensing_noise_spec : too high multipole");

  ncm_vector_memcpy (vp2, vp1);

  {
    NcmVector *Nl_sub = ncm_vector_get_subvector (xclkl->Nl, lmin, ncm_vector_len (vp1));

    ncm_vector_add (vp2, Nl_sub);
    ncm_vector_free (Nl_sub);
  }
  /* return ncm_vector_get (xclkl->Nl, l); */
}

static guint
_nc_xcor_kernel_cmb_lensing_obs_len (NcXcorKernel *xclk)
{
  return 1;
}

static guint
_nc_xcor_kernel_cmb_lensing_obs_params_len (NcXcorKernel *xclk)
{
  return 0;
}

static void
_nc_xcor_kernel_cmb_lensing_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid)
{
  NcXcorKernelCMBLensing *xclkl = NC_XCOR_KERNEL_CMB_LENSING (xclk);

  *zmin = 0.0;
  *zmax = (xclkl->source != NC_XCOR_KERNEL_CMB_LENSING_SOURCE_THIN_SCREEN) ? xclkl->z_src_max : xclkl->z_lss;
  *zmid = 2.0;
}

static GPtrArray *
_nc_xcor_kernel_cmb_lensing_get_component_list (NcXcorKernel *xclk)
{
  NcXcorKernelCMBLensing *xclkl = NC_XCOR_KERNEL_CMB_LENSING (xclk);
  GPtrArray *comp_list          = g_ptr_array_new_with_free_func (g_object_unref);

  if (xclkl->cmb_lens_comp != NULL)
    g_ptr_array_add (comp_list, g_object_ref (xclkl->cmb_lens_comp));

  return comp_list;
}

/**
 * nc_xcor_kernel_cmb_lensing_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @recomb: a #NcRecomb
 * @Nl: a #NcmVector
 *
 * Creates a new #NcXcorLimberKernelCMBLensing for computing the CMB lensing
 * convergence kernel. This kernel describes the lensing of CMB photons by
 * intervening large-scale structure.
 *
 * Returns: a new #NcXcorLimberKernelCMBLensing
 *
 */
NcXcorKernelCMBLensing *
nc_xcor_kernel_cmb_lensing_new (NcDistance *dist, NcmPowspec *ps, NcRecomb *recomb, NcmVector *Nl) /*, gdouble zl, gdouble zu) */
{
  NcXcorKernelCMBLensing *xclkl = g_object_new (NC_TYPE_XCOR_KERNEL_CMB_LENSING,
                                                "dist", dist,
                                                "powspec", ps,
                                                "recomb", recomb,
                                                "Nl", Nl,
                                                NULL);

  return xclkl;
}

/**
 * nc_xcor_kernel_cmb_lensing_set_source:
 * @xclkl: a #NcXcorKernelCMBLensing
 * @source: a #NcXcorKernelCMBLensingSource
 *
 * Sets where the CMB photons are placed along the line of sight. The visibility
 * source requires the #NcXcorKernelCMBLensing:recomb object. The kernel is
 * marked outdated and is prepared again on the next use.
 */
void
nc_xcor_kernel_cmb_lensing_set_source (NcXcorKernelCMBLensing *xclkl, NcXcorKernelCMBLensingSource source)
{
  if (xclkl->source != source)
  {
    xclkl->source = source;
    nc_xcor_kernel_mark_outdated (NC_XCOR_KERNEL (xclkl));
  }
}

/**
 * nc_xcor_kernel_cmb_lensing_get_source:
 * @xclkl: a #NcXcorKernelCMBLensing
 *
 * Returns: the #NcXcorKernelCMBLensingSource in use.
 */
NcXcorKernelCMBLensingSource
nc_xcor_kernel_cmb_lensing_get_source (NcXcorKernelCMBLensing *xclkl)
{
  return xclkl->source;
}

