/***************************************************************************
 *            nc_xcor_limber_kernel_tSZ.c
 *
 *  Tue January 10 12:00:00 2022
 *  Copyright  2022  Arthur de Souza Molina
 *  <arthur.souza.molina@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2022 Arthur de Souza Molina <arthur.souza.molina@uel.br>
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
 * NcXcorLimberKerneltSZ:
 *
 * Thermal Sunyaev Zel'dovich effect kernel.
 *
 * The thermal Sunyaev Zel'dovich (tSZ) effect is a modification in the observed
 * temperature of the cosmic microwave background (CMB) due to the inverse Compton
 * scattering of CMB photons with high-energy electrons along the line-of-sight. These
 * electrons are present in the intracluster medium (ICM) of galaxy clusters, for
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
#include "xcor/nc_xcor_limber_kernel_tSZ.h"
#include "xcor/nc_xcor.h"

#ifndef NUMCOSMO_GIR_SCAN

#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorLimberKerneltSZ
{
  /*< private >*/
  NcXcorLimberKernel parent_instance;
  gdouble noise;
};


G_DEFINE_TYPE (NcXcorLimberKerneltSZ, nc_xcor_limber_kernel_tsz, NC_TYPE_XCOR_LIMBER_KERNEL);

#define VECTOR (NCM_MODEL (xclkl)->params)

enum
{
  PROP_0,
  PROP_NOISE,
  PROP_SIZE,
};

static void
nc_xcor_limber_kernel_tsz_init (NcXcorLimberKerneltSZ *xclkl)
{
  xclkl->noise = 0.0;
}

static void
_nc_xcor_limber_kernel_tsz_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorLimberKerneltSZ *xclkl = NC_XCOR_LIMBER_KERNEL_TSZ (object);

  g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL_TSZ (object));

  switch (prop_id)
  {
    case PROP_NOISE:
      xclkl->noise = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_xcor_limber_kernel_tsz_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorLimberKerneltSZ *xclkl = NC_XCOR_LIMBER_KERNEL_TSZ (object);

  g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL_TSZ (object));

  switch (prop_id)
  {
    case PROP_NOISE:
      g_value_set_double (value, xclkl->noise);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_xcor_limber_kernel_tsz_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_limber_kernel_tsz_parent_class)->dispose (object);
}

static void
_nc_xcor_limber_kernel_tsz_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_limber_kernel_tsz_parent_class)->finalize (object);
}

static gdouble _nc_xcor_limber_kernel_tsz_eval (NcXcorLimberKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
static void _nc_xcor_limber_kernel_tsz_prepare (NcXcorLimberKernel *xclk, NcHICosmo *cosmo);
static void _nc_xcor_limber_kernel_tsz_add_noise (NcXcorLimberKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
static guint _nc_xcor_limber_kernel_tsz_obs_len (NcXcorLimberKernel *xclk);
static guint _nc_xcor_limber_kernel_tsz_obs_params_len (NcXcorLimberKernel *xclk);

static void
nc_xcor_limber_kernel_tsz_class_init (NcXcorLimberKerneltSZClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcXcorLimberKernelClass *parent_class = NC_XCOR_LIMBER_KERNEL_CLASS (klass);
  NcmModelClass *model_class            = NCM_MODEL_CLASS (klass);

  object_class->finalize    = &_nc_xcor_limber_kernel_tsz_finalize;
  object_class->dispose     = &_nc_xcor_limber_kernel_tsz_dispose;
  model_class->set_property = &_nc_xcor_limber_kernel_tsz_set_property;
  model_class->get_property = &_nc_xcor_limber_kernel_tsz_get_property;

  ncm_model_class_set_name_nick (model_class, "Xcor with the thermal Sunyaev Zel'dovich Compton-y parameter", "Xcor-tSZ");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcXcorLimberKerneltSZ:zmax:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_NOISE,
                                   g_param_spec_double ("noise",
                                                        NULL,
                                                        "Constant noise level",
                                                        -10.0, 10.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->eval      = &_nc_xcor_limber_kernel_tsz_eval;
  parent_class->prepare   = &_nc_xcor_limber_kernel_tsz_prepare;
  parent_class->add_noise = &_nc_xcor_limber_kernel_tsz_add_noise;

  parent_class->obs_len        = &_nc_xcor_limber_kernel_tsz_obs_len;
  parent_class->obs_params_len = &_nc_xcor_limber_kernel_tsz_obs_params_len;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_LIMBER_KERNEL_IMPL_ALL);
}

/**
 * nc_xcor_limber_kernel_tsz_new:
 * @zmax: a gdouble
 *
 * Creates a new instance of the tSZ kernel.
 *
 * Returns: (transfer full): a new #NcXcorLimberKerneltSZ
 */
NcXcorLimberKerneltSZ *
nc_xcor_limber_kernel_tsz_new (gdouble zmax)
{
  NcXcorLimberKerneltSZ *xclkl = g_object_new (NC_TYPE_XCOR_LIMBER_KERNEL_TSZ,
                                               "zmin", 0.0,
                                               "zmax", zmax,
                                               NULL);

  return xclkl;
}

static gdouble
_nc_xcor_limber_kernel_tsz_eval (NcXcorLimberKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  const gdouble nc_prefac = ncm_c_thomson_cs () / (ncm_c_mass_e () * ncm_c_c () * ncm_c_c ());
  const gdouble units     = nc_hicosmo_RH_Mpc (cosmo) * ncm_c_Mpc () * ncm_c_eV () / (1.0e-6);

  /* unit: eV / cm^3 */

  return nc_prefac * units / (1.0 + z) / xck->E_z;
}

static void
_nc_xcor_limber_kernel_tsz_prepare (NcXcorLimberKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorLimberKerneltSZ *xclkl = NC_XCOR_LIMBER_KERNEL_TSZ (xclk);

  g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL_TSZ (xclkl));
  nc_xcor_limber_kernel_set_const_factor (xclk, 1.0);
}

static void
_nc_xcor_limber_kernel_tsz_add_noise (NcXcorLimberKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  NcXcorLimberKerneltSZ *xclkl = NC_XCOR_LIMBER_KERNEL_TSZ (xclk);

  ncm_vector_memcpy (vp2, vp1);
  ncm_vector_add_constant (vp2, xclkl->noise);
}

static guint
_nc_xcor_limber_kernel_tsz_obs_len (NcXcorLimberKernel *xclk)
{
  return 1;
}

static guint
_nc_xcor_limber_kernel_tsz_obs_params_len (NcXcorLimberKernel *xclk)
{
  return 0;
}

