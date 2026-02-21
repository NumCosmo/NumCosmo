/***************************************************************************
 *            nc_xcor_kernel_tSZ.c
 *
 *  Tue January 10 12:00:00 2022
 *  Copyright  2022  Arthur de Souza Molina
 *  <arthur.souza.molina@uel.br>
 *  Sat December 27 20:21:01 2025
 *  Copyright  2025  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2022 Arthur de Souza Molina <arthur.souza.molina@uel.br>
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
 * NcXcorKerneltSZ:
 *
 * Thermal Sunyaev Zel'dovich effect kernel.
 *
 * The thermal Sunyaev Zel'dovich (tSZ) effect is a modification in the observed
 * temperature of the cosmic microwave background (CMB) due to the inverse Compton
 * scattering of CMB photons with high-energy electrons along the line-of-sight. These
 * electrons are present in the intra-cluster medium (ICM) of galaxy clusters, for
 * example.
 *
 * ## Compton-y parameter
 * The kernel is given by
 * \begin{equation}
 *    W(\\chi) = \\frac{\\sigma_T}{m_ec^2} \\frac{1}{1+z}.
 * \end{equation}
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include "math/ncm_cfg.h"
#include "xcor/nc_xcor_kernel_tSZ.h"
#include "xcor/nc_xcor_kernel_component.h"
#include "xcor/nc_xcor.h"

#ifndef NUMCOSMO_GIR_SCAN

#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorKerneltSZ
{
  /*< private >*/
  NcXcorKernel parent_instance;
  gdouble noise;
  gdouble zmax;
  NcXcorKernelComponent *tsz_comp;
};

G_DEFINE_TYPE (NcXcorKerneltSZ, nc_xcor_kernel_tsz, NC_TYPE_XCOR_KERNEL);

/*
 * tSZ Component Definition
 * Handles the thermal Sunyaev Zel'dovich effect kernel
 */

typedef struct _tSZComponentData
{
  NcDistance *dist;
  NcmPowspec *ps;
  gdouble zmax;
} tSZComponentData;

#define _NC_XCOR_KERNEL_COMPONENT_TSZ_GET_DATA(comp) \
        ((tSZComponentData *) ((guint8 *) (comp) + sizeof (NcXcorKernelComponent)))

static gdouble _tsz_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble xi, gdouble k);
static gdouble _tsz_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l);
static void _tsz_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *xi_min, gdouble *xi_max, gdouble *k_min, gdouble *k_max);
static void _tsz_component_data_clear (tSZComponentData *data);
static NcXcorKernelComponent *_nc_xcor_kernel_component_tsz_new (NcDistance *dist, NcmPowspec *ps, gdouble zmax);

NC_XCOR_KERNEL_COMPONENT_DEFINE_TYPE (NC, XCOR_KERNEL_COMPONENT_TSZ,
                                      NcXcorKernelComponenttSZ,
                                      nc_xcor_kernel_component_tsz,
                                      _tsz_component_eval_kernel,
                                      _tsz_component_eval_prefactor,
                                      _tsz_component_get_limits,
                                      tSZComponentData,
                                      _tsz_component_data_clear)

enum
{
  PROP_0,
  PROP_NOISE,
  PROP_ZMAX,
  PROP_SIZE,
};

static void
nc_xcor_kernel_tsz_init (NcXcorKerneltSZ *xclkl)
{
  xclkl->noise    = 0.0;
  xclkl->zmax     = 0.0;
  xclkl->tsz_comp = NULL;
}

static void
_nc_xcor_kernel_tsz_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKerneltSZ *xclkl = NC_XCOR_KERNEL_TSZ (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_TSZ (object));

  switch (prop_id)
  {
    case PROP_NOISE:
      xclkl->noise = g_value_get_double (value);
      break;
    case PROP_ZMAX:
      xclkl->zmax = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_kernel_tsz_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKerneltSZ *xclkl = NC_XCOR_KERNEL_TSZ (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_TSZ (object));

  switch (prop_id)
  {
    case PROP_NOISE:
      g_value_set_double (value, xclkl->noise);
      break;
    case PROP_ZMAX:
      g_value_set_double (value, xclkl->zmax);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_kernel_tsz_dispose (GObject *object)
{
  NcXcorKerneltSZ *xclkl = NC_XCOR_KERNEL_TSZ (object);

  nc_xcor_kernel_component_clear (&xclkl->tsz_comp);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_tsz_parent_class)->dispose (object);
}

static void
_nc_xcor_kernel_tsz_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_tsz_parent_class)->finalize (object);
}

static gdouble _nc_xcor_kernel_tsz_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
static gdouble _nc_xcor_kernel_tsz_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
static void _nc_xcor_kernel_tsz_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo);
static void _nc_xcor_kernel_tsz_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
static guint _nc_xcor_kernel_tsz_obs_len (NcXcorKernel *xclk);
static guint _nc_xcor_kernel_tsz_obs_params_len (NcXcorKernel *xclk);
static void _nc_xcor_kernel_tsz_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid);
static GPtrArray *_nc_xcor_kernel_tsz_get_component_list (NcXcorKernel *xclk);

static void
_nc_xcor_kernel_tsz_constructed (GObject *object)
{
  G_OBJECT_CLASS (nc_xcor_kernel_tsz_parent_class)->constructed (object);
  {
    NcXcorKerneltSZ *xclkl = NC_XCOR_KERNEL_TSZ (object);
    NcDistance *dist       = nc_xcor_kernel_peek_dist (NC_XCOR_KERNEL (object));
    NcmPowspec *ps         = nc_xcor_kernel_peek_powspec (NC_XCOR_KERNEL (object));

    g_assert_null (xclkl->tsz_comp);
    xclkl->tsz_comp = _nc_xcor_kernel_component_tsz_new (dist, ps, xclkl->zmax);
  }
}

static void
nc_xcor_kernel_tsz_class_init (NcXcorKerneltSZClass *klass)
{
  GObjectClass *object_class      = G_OBJECT_CLASS (klass);
  NcXcorKernelClass *parent_class = NC_XCOR_KERNEL_CLASS (klass);
  NcmModelClass *model_class      = NCM_MODEL_CLASS (klass);

  object_class->constructed = &_nc_xcor_kernel_tsz_constructed;
  object_class->finalize    = &_nc_xcor_kernel_tsz_finalize;
  object_class->dispose     = &_nc_xcor_kernel_tsz_dispose;
  model_class->set_property = &_nc_xcor_kernel_tsz_set_property;
  model_class->get_property = &_nc_xcor_kernel_tsz_get_property;

  ncm_model_class_set_name_nick (model_class, "Xcor with the thermal Sunyaev Zel'dovich Compton-y parameter", "Xcor-tSZ");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorKerneltSZ:noise:
   *
   * Constant noise level.
   */
  g_object_class_install_property (object_class,
                                   PROP_NOISE,
                                   g_param_spec_double ("noise",
                                                        NULL,
                                                        "Constant noise level",
                                                        -10.0, 10.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorKerneltSZ:zmax:
   *
   * Maximum redshift for the kernel.
   */
  g_object_class_install_property (object_class,
                                   PROP_ZMAX,
                                   g_param_spec_double ("zmax",
                                                        NULL,
                                                        "Maximum redshift",
                                                        0.0, 10.0, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->eval_limber_z           = &_nc_xcor_kernel_tsz_eval_limber_z;
  parent_class->eval_limber_z_prefactor = &_nc_xcor_kernel_tsz_eval_limber_z_prefactor;
  parent_class->prepare                 = &_nc_xcor_kernel_tsz_prepare;
  parent_class->add_noise               = &_nc_xcor_kernel_tsz_add_noise;

  parent_class->obs_len            = &_nc_xcor_kernel_tsz_obs_len;
  parent_class->obs_params_len     = &_nc_xcor_kernel_tsz_obs_params_len;
  parent_class->get_z_range        = &_nc_xcor_kernel_tsz_get_z_range;
  parent_class->get_component_list = &_nc_xcor_kernel_tsz_get_component_list;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

static gdouble
_nc_xcor_kernel_tsz_eval_radial_weight (NcXcorKernel *xclk, NcHICosmo *cosmo, const gdouble z, const gdouble xi, const gdouble E)
{
  return 1.0 / (1.0 + z);
}

static gdouble
_nc_xcor_kernel_tsz_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  const gdouble kernel = _nc_xcor_kernel_tsz_eval_radial_weight (xclk, cosmo, z, xck->xi_z, xck->E_z);

  return kernel;
}

static gdouble
_nc_xcor_kernel_tsz_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  const gdouble nc_pre_fac = ncm_c_thomson_cs () / (ncm_c_mass_e () * ncm_c_c () * ncm_c_c ());
  const gdouble units      = nc_hicosmo_RH_Mpc (cosmo) * ncm_c_Mpc () * ncm_c_eV () / (1.0e-6);

  /* unit: eV / cm^3 */
  return nc_pre_fac * units;
}

/*
 * Implementation of the new Component interface.
 */

static void
_tsz_component_data_clear (tSZComponentData *data)
{
  /* No need to clear, these are weak references from parent kernel */
}

static gdouble
_tsz_component_eval_kernel (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble xi, gdouble k)
{
  tSZComponentData *data = _NC_XCOR_KERNEL_COMPONENT_TSZ_GET_DATA (comp);
  const gdouble z        = nc_distance_inv_comoving (data->dist, cosmo, xi);
  const gdouble powspec  = ncm_powspec_eval (data->ps, NCM_MODEL (cosmo), z, k / nc_hicosmo_RH_Mpc (cosmo));
  const gdouble kernel   = 1.0 / (1.0 + z);

  return kernel * sqrt (powspec);
}

static gdouble
_tsz_component_eval_prefactor (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble k, gint l)
{
  const gdouble nc_pre_fac = ncm_c_thomson_cs () / (ncm_c_mass_e () * ncm_c_c () * ncm_c_c ());
  const gdouble units      = nc_hicosmo_RH_Mpc (cosmo) * ncm_c_Mpc () * ncm_c_eV () / (1.0e-6);

  return nc_pre_fac * units;
}

static void
_tsz_component_get_limits (NcXcorKernelComponent *comp, NcHICosmo *cosmo, gdouble *xi_min, gdouble *xi_max, gdouble *k_min, gdouble *k_max)
{
  tSZComponentData *data = _NC_XCOR_KERNEL_COMPONENT_TSZ_GET_DATA (comp);
  NcDistance *dist       = data->dist;
  NcmPowspec *ps         = data->ps;

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  *xi_min = nc_distance_comoving (dist, cosmo, 1.0e-2);
  *xi_max = nc_distance_comoving (dist, cosmo, data->zmax);
  *k_min  = ncm_powspec_get_kmin (ps) * nc_hicosmo_RH_Mpc (cosmo);
  *k_max  = ncm_powspec_get_kmax (ps) * nc_hicosmo_RH_Mpc (cosmo);
}

static NcXcorKernelComponent *
_nc_xcor_kernel_component_tsz_new (NcDistance *dist, NcmPowspec *ps, gdouble zmax)
{
  NcXcorKernelComponent *comp = g_object_new (nc_xcor_kernel_component_tsz_get_type (), NULL);
  tSZComponentData *data      = _NC_XCOR_KERNEL_COMPONENT_TSZ_GET_DATA (comp);

  data->dist = dist;
  data->ps   = ps;
  data->zmax = zmax;

  return comp;
}

/*
 * Kernel implementation.
 */

static void
_nc_xcor_kernel_tsz_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorKerneltSZ *xclkl = NC_XCOR_KERNEL_TSZ (xclk);
  NcDistance *dist       = nc_xcor_kernel_peek_dist (xclk);
  NcmPowspec *ps         = nc_xcor_kernel_peek_powspec (xclk);

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  /* Update component data with zmax if it changed */
  {
    tSZComponentData *data = _NC_XCOR_KERNEL_COMPONENT_TSZ_GET_DATA (xclkl->tsz_comp);

    data->zmax = xclkl->zmax;
  }

  g_assert_nonnull (xclkl->tsz_comp);
  nc_xcor_kernel_component_prepare (xclkl->tsz_comp, cosmo);
}

static void
_nc_xcor_kernel_tsz_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  NcXcorKerneltSZ *xclkl = NC_XCOR_KERNEL_TSZ (xclk);

  ncm_vector_memcpy (vp2, vp1);
  ncm_vector_add_constant (vp2, xclkl->noise);
}

static guint
_nc_xcor_kernel_tsz_obs_len (NcXcorKernel *xclk)
{
  return 1;
}

static guint
_nc_xcor_kernel_tsz_obs_params_len (NcXcorKernel *xclk)
{
  return 0;
}

static void
_nc_xcor_kernel_tsz_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid)
{
  NcXcorKerneltSZ *xclkl = NC_XCOR_KERNEL_TSZ (xclk);

  *zmin = 0.0;
  *zmax = xclkl->zmax;
  *zmid = xclkl->zmax / 2.0;
}

static GPtrArray *
_nc_xcor_kernel_tsz_get_component_list (NcXcorKernel *xclk)
{
  NcXcorKerneltSZ *xclkl = NC_XCOR_KERNEL_TSZ (xclk);
  GPtrArray *comp_list   = g_ptr_array_new_with_free_func (g_object_unref);

  if (xclkl->tsz_comp != NULL)
    g_ptr_array_add (comp_list, g_object_ref (xclkl->tsz_comp));

  return comp_list;
}

/**
 * nc_xcor_kernel_tsz_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @zmax: a gdouble
 *
 * Creates a new instance of the tSZ kernel.
 *
 * Returns: (transfer full): a new #NcXcorKerneltSZ
 */
NcXcorKerneltSZ *
nc_xcor_kernel_tsz_new (NcDistance *dist, NcmPowspec *ps, gdouble zmax)
{
  NcXcorKerneltSZ *xclkl = g_object_new (NC_TYPE_XCOR_KERNEL_TSZ,
                                         "dist", dist,
                                         "powspec", ps,
                                         "zmax", zmax,
                                         NULL);

  return xclkl;
}

