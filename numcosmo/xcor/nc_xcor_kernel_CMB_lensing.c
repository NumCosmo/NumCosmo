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
#include "xcor/nc_xcor_kernel_CMB_lensing.h"
#include "xcor/nc_xcor_kernel_component.h"
#include "xcor/nc_xcor.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */


struct _NcXcorKernelCMBLensing
{
  /*< private >*/
  NcXcorKernel parent_instance;

  NcRecomb *recomb;

  NcmVector *Nl;
  guint Nlmax;

  gdouble z_lss;
  gdouble xi_lss;
  gdouble dt_lss;
  gdouble dt;

  NcDistance *dist;
  NcmPowspec *ps;
  NcXcorKernelComponent *cmb_lens_comp;
};

enum
{
  PROP_0,
  PROP_RECOMB,
  PROP_NL,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelCMBLensing, nc_xcor_kernel_cmb_lensing, NC_TYPE_XCOR_KERNEL)

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
} CMBLensingComponentData;

#define _NC_XCOR_KERNEL_COMPONENT_CMB_LENSING_GET_DATA(comp) \
        ((CMBLensingComponentData *) ((guint8 *) (comp) + sizeof (NcXcorKernelComponent)))

static gdouble _cmb_lensing_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble xi, gdouble k);
static gdouble _cmb_lensing_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l);
static void _cmb_lensing_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *xi_min, gdouble *xi_max, gdouble *k_min, gdouble *k_max);
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
  xclkl->xi_lss        = 0.0;
  xclkl->dt_lss        = 0.0;
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
   * NcXcorKernelCMBLensing:Nl:
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
  const gdouble dt_z_zlss       = nc_distance_transverse_z1_z2 (dist, cosmo, z, xclkl->z_lss);

  return xck->xi_z * xck->xi_z * (1.0 + z) * dt_z_zlss / (xclkl->dt_lss * dt);
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
_cmb_lensing_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble xi, gdouble k)
{
  CMBLensingComponentData *data = _NC_XCOR_KERNEL_COMPONENT_CMB_LENSING_GET_DATA (comp);
  const gdouble z               = nc_distance_inv_comoving (data->dist, cosmo, xi);
  const gdouble powspec         = ncm_powspec_eval (data->ps, NCM_MODEL (cosmo), z, k / nc_hicosmo_RH_Mpc (cosmo));
  const gdouble dt              = nc_distance_transverse (data->dist, cosmo, z);
  const gdouble dt_z_zlss       = nc_distance_transverse_z1_z2 (data->dist, cosmo, z, data->z_lss);
  const gdouble operator_k      = 1.0 / gsl_pow_3 (k);

  return operator_k * (1.0 + z) * dt_z_zlss / (data->dt_lss * dt) * sqrt (powspec);
}

static gdouble
_cmb_lensing_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l)
{
  const gdouble lfactor = l * (l + 1.0);

  return 1.5 * nc_hicosmo_Omega_m0 (cosmo) * lfactor;
}

static void
_cmb_lensing_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *xi_min, gdouble *xi_max, gdouble *k_min, gdouble *k_max)
{
  CMBLensingComponentData *data = _NC_XCOR_KERNEL_COMPONENT_CMB_LENSING_GET_DATA (comp);
  NcDistance *dist              = data->dist;
  NcmPowspec *ps                = data->ps;

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  {
    const gdouble xi_lss = nc_distance_comoving_lss (dist, cosmo);

    *xi_min = nc_distance_comoving (dist, cosmo, 1.0e-2);
    *xi_max = xi_lss;
    *k_min  = ncm_powspec_get_kmin (ps) * nc_hicosmo_RH_Mpc (cosmo);
    *k_max  = ncm_powspec_get_kmax (ps) * nc_hicosmo_RH_Mpc (cosmo);
  }
}

static NcXcorKernelComponent *
_nc_xcor_kernel_component_cmb_lensing_new (NcDistance *dist, NcmPowspec *ps)
{
  NcXcorKernelComponent *comp   = g_object_new (nc_xcor_kernel_component_cmb_lensing_get_type (), NULL);
  CMBLensingComponentData *data = _NC_XCOR_KERNEL_COMPONENT_CMB_LENSING_GET_DATA (comp);

  data->dist   = dist;
  data->ps     = ps;
  data->z_lss  = 0.0;
  data->dt_lss = 0.0;

  return comp;
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

  xclkl->z_lss  = nc_distance_decoupling_redshift (dist, cosmo);
  xclkl->xi_lss = nc_distance_comoving_lss (dist, cosmo);
  xclkl->dt_lss = nc_distance_transverse (dist, cosmo, xclkl->z_lss);

  /* Update component data with computed values */
  {
    CMBLensingComponentData *data = _NC_XCOR_KERNEL_COMPONENT_CMB_LENSING_GET_DATA (xclkl->cmb_lens_comp);

    data->z_lss  = xclkl->z_lss;
    data->dt_lss = xclkl->dt_lss;
  }

  g_assert_nonnull (xclkl->cmb_lens_comp);
  nc_xcor_kernel_component_prepare (xclkl->cmb_lens_comp, cosmo);

  /* nc_recomb_prepare (xclkl->recomb, cosmo); */
  /* gdouble lamb = nc_recomb_tau_zstar (xclkl->recomb, cosmo); */
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
  *zmax = xclkl->z_lss;
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

