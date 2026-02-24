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

#include "math/ncm_cfg.h"
#include "math/ncm_integrate.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_spline2d_bicubic.h"
#include "xcor/nc_xcor_kernel_component.h"
#include "xcor/nc_xcor_kernel_cmb_isw.h"
#include "xcor/nc_xcor.h"


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
  gdouble xi_lss;
  gdouble z_lss;
  NcXcorKernelComponent *isw_comp;
} NcXcorKernelCMBISWPrivate;

enum
{
  PROP_0,
  PROP_RECOMB,
  PROP_NL,
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
} ISWComponentData;

/* Helper to get data from component - uses pointer arithmetic to access
 * the data member that comes after the parent_instance in the struct
 * defined by the macro below */
#define _NC_XCOR_KERNEL_COMPONENT_ISW_GET_DATA(comp) \
        ((ISWComponentData *) ((guint8 *) (comp) + sizeof (NcXcorKernelComponent)))


static gdouble _isw_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble xi, gdouble k);
static gdouble _isw_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l);
static void _isw_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *xi_min, gdouble *xi_max, gdouble *k_min, gdouble *k_max);
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

  self->recomb   = NULL;
  self->Nl       = NULL;
  self->Nlmax    = 0;
  self->xi_lss   = 0.0;
  self->z_lss    = 0.0;
  self->isw_comp = NULL;
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
  *zmax = self->z_lss;
  *zmid = 2.0;
}

static gdouble
_nc_xcor_kernel_cmb_isw_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  NcmPowspec *ps               = nc_xcor_kernel_peek_powspec (xclk);
  const gdouble k_pivot        = 1.0;
  const gdouble powspec        = ncm_powspec_eval (ps, NCM_MODEL (cosmo), z, k_pivot);
  const gdouble dpowspec_dz    = ncm_powspec_deriv_z (ps, NCM_MODEL (cosmo), z, k_pivot);
  const gdouble d1pz_growth_dz = 1.0 + (1.0 + z) * dpowspec_dz / (2.0 * powspec);

  return xck->E_z * gsl_pow_2 (xck->xi_z) * d1pz_growth_dz;
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
_isw_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble xi, gdouble k)
{
  /* Access data using offset from base class */
  ISWComponentData *data       = _NC_XCOR_KERNEL_COMPONENT_ISW_GET_DATA (comp);
  const gdouble z              = nc_distance_inv_comoving (data->dist, cosmo, xi);
  const gdouble E_z            = nc_hicosmo_E (cosmo, z);
  const gdouble powspec        = ncm_powspec_eval (data->ps, NCM_MODEL (cosmo), z, k / nc_hicosmo_RH_Mpc (cosmo));
  const gdouble dpowspec_dz    = ncm_powspec_deriv_z (data->ps, NCM_MODEL (cosmo), z, k / nc_hicosmo_RH_Mpc (cosmo));
  const gdouble d1pz_growth_dz = 1.0 + (1.0 + z) * dpowspec_dz / (2.0 * powspec);
  const gdouble operator       = 1.0 / gsl_pow_2 (k);

  return operator * E_z * d1pz_growth_dz * sqrt (powspec);
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
_isw_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *xi_min, gdouble *xi_max, gdouble *k_min, gdouble *k_max)
{
  ISWComponentData *data = _NC_XCOR_KERNEL_COMPONENT_ISW_GET_DATA (comp);
  NcDistance *dist       = data->dist;
  NcmPowspec *ps         = data->ps;

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  {
    const gdouble xi_lss = nc_distance_comoving_lss (dist, cosmo);

    *xi_min = nc_distance_comoving (dist, cosmo, 1.0e-6);
    *xi_max = xi_lss;
    *k_min  = ncm_powspec_get_kmin (ps) * nc_hicosmo_RH_Mpc (cosmo);
    *k_max  = ncm_powspec_get_kmax (ps) * nc_hicosmo_RH_Mpc (cosmo);
  }
}

static NcXcorKernelComponent *
_nc_xcor_kernel_component_isw_new (NcDistance *dist, NcmPowspec *ps)
{
  NcXcorKernelComponent *comp = g_object_new (nc_xcor_kernel_component_isw_get_type (), NULL);
  ISWComponentData *data      = _NC_XCOR_KERNEL_COMPONENT_ISW_GET_DATA (comp);

  data->dist = dist;
  data->ps   = ps;

  return comp;
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

  self->xi_lss = nc_distance_comoving_lss (dist, cosmo);
  self->z_lss  = z_lss;

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
 * @y: the y value (y = k * xi)
 *
 * Evaluates k_max at the given y value from kernel analysis.
 *
 * Returns: the k_max value at y
 */
gdouble
nc_xcor_kernel_cmb_isw_eval_k_max (NcXcorKernelCMBISW *xcisw, gdouble y)
{
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  g_assert (self->isw_comp != NULL);

  return nc_xcor_kernel_component_eval_k_max (self->isw_comp, y);
}

/**
 * nc_xcor_kernel_cmb_isw_eval_K_max:
 * @xcisw: a #NcXcorKernelCMBISW
 * @y: the y value (y = k * xi)
 *
 * Evaluates K_max at the given y value from kernel analysis.
 *
 * Returns: the K_max value at y
 */
gdouble
nc_xcor_kernel_cmb_isw_eval_K_max (NcXcorKernelCMBISW *xcisw, gdouble y)
{
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  g_assert (self->isw_comp != NULL);

  return nc_xcor_kernel_component_eval_K_max (self->isw_comp, y);
}

/**
 * nc_xcor_kernel_cmb_isw_eval_k_epsilon:
 * @xcisw: a #NcXcorKernelCMBISW
 * @y: the y value (y = k * xi)
 *
 * Evaluates k_epsilon at the given y value from kernel analysis.
 *
 * Returns: the k_epsilon value at y
 */
gdouble
nc_xcor_kernel_cmb_isw_eval_k_epsilon (NcXcorKernelCMBISW *xcisw, gdouble y)
{
  NcXcorKernelCMBISWPrivate * const self = nc_xcor_kernel_cmb_isw_get_instance_private (xcisw);

  g_assert (self->isw_comp != NULL);

  return nc_xcor_kernel_component_eval_k_epsilon (self->isw_comp, y);
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

