/***************************************************************************
 *            nc_xcor_limber_kernel_tSZ.c
 *
 *  Tue January 10 12:00:00 2021
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
 * SECTION:nc_xcor_limber_kernel_tSZ
 * @title: NcXcorLimberKerneltSZ
 * @short_description: implementation of #NcXcorLimberKernel with the thermal Sunyaev Zel'dovich
    Compton-y parameter
 *
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

G_DEFINE_TYPE (NcXcorLimberKerneltSZ, nc_xcor_limber_kernel_tsz, NC_TYPE_XCOR_LIMBER_KERNEL);

#define VECTOR (NCM_MODEL (xclkl)->params)

enum
{
  PROP_0,
  PROP_ZMAX,
  PROP_NCHI,
  PROP_SIZE,
};

static void
nc_xcor_limber_kernel_tsz_init (NcXcorLimberKerneltSZ *xclkl)
{
  xclkl->Zmax = 6.0;
  xclkl->Nchi = 1024;
}

static void
_nc_xcor_limber_kernel_tsz_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorLimberKerneltSZ *xclkl = NC_XCOR_LIMBER_KERNEL_TSZ (object);
  
  g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL_TSZ (object));
  
  switch (prop_id)
  {
    case PROP_ZMAX:
      xclkl->Zmax = g_value_get_double(value);
      break;

    case PROP_NCHI:
      xclkl->Nchi = g_value_get_double (value);
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
    case PROP_ZMAX:
      g_value_set_double (value, xclkl->Zmax);
      break;
    case PROP_NCHI:
      g_value_set_double(value, xclkl->Nchi);
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
  ncm_model_class_add_params (model_class, NC_XCOR_LIMBER_KERNEL_TSZ_SPARAM_LEN, 0, PROP_SIZE);
  
  /**
   * NcXcorLimberKerneltSZ:Zmax:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_ZMAX,
                                   g_param_spec_object ("Zmax",
                                                        NULL,
                                                        "Maximum redshift up to which we define the kernel",
                                                        NC_TYPE_XCOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));



  /**
   * NcXcorLimberKerneltSZ:Nchi:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_NCHI,
                                   g_param_spec_object ("Nchi",
                                                        NULL,
                                                        "Number of intervals in the radial comoving distance on which we sample the kernel",
                                                        NC_TYPE_XCOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);
  
  parent_class->eval    = &_nc_xcor_limber_kernel_tsz_eval;
  parent_class->prepare = &_nc_xcor_limber_kernel_tsz_prepare;
  parent_class->add_noise = &_nc_xcor_limber_kernel_tsz_add_noise;
  
  parent_class->obs_len        = &_nc_xcor_limber_kernel_tsz_obs_len;
  parent_class->obs_params_len = &_nc_xcor_limber_kernel_tsz_obs_params_len;
  
  ncm_model_class_add_impl_flag (model_class, NC_XCOR_LIMBER_KERNEL_IMPL_ALL);
}

/**
 * nc_xcor_limber_kernel_tsz_new:
 * @Zmax: a #Gdouble
 * @Nchi: a #Gdouble
 *
 * FIXME
 *
 * Returns: FIXME
 *
 */
NcXcorLimberKerneltSZ *
nc_xcor_limber_kernel_tsz_new (gdouble Zmax, gdouble Nchi)
{
  NcXcorLimberKerneltSZ *xclkl = g_object_new (NC_TYPE_XCOR_LIMBER_KERNEL_TSZ,
                                                      "Zmax", Zmax,
													  "Nchi", Nchi,
                                                      NULL);
  
  return xclkl;
}

static gdouble
_nc_xcor_limber_kernel_tsz_eval (NcXcorLimberKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{

  NcXcorLimberKerneltSZ *xclkl = NC_XCOR_LIMBER_KERNEL_TSZ (xclk);

  NCM_UNUSED (l);
  NCM_UNUSED (cosmo);

  const gdouble nc_prefac = ncm_c_thomson_cs() / ( ncm_c_mass_e() * ncm_c_c() * ncm_c_c() );

  return nc_prefac * 1.0 / (1+xclkl->Zmax);
}

static void
_nc_xcor_limber_kernel_tsz_prepare (NcXcorLimberKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorLimberKerneltSZ *xclkl = NC_XCOR_LIMBER_KERNEL_TSZ (xclk);
  g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL_TSZ (xclkl));

}

static void
_nc_xcor_limber_kernel_tsz_add_noise (NcXcorLimberKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  NcXcorLimberKerneltSZ *xclkl = NC_XCOR_LIMBER_KERNEL_TSZ (xclk);

  if (xclkl->Zmax <= 0){
	  g_error ("nc_xcor_limber_kernel_tsz_noise_spec : Zmax Error.");
  }
  if (xclkl->Nchi == 0){
  	  g_error ("N_chi Error: N_chi cannot be zero.");
    }
  NCM_UNUSED (xclk);
}

static guint
_nc_xcor_limber_kernel_tsz_obs_len (NcXcorLimberKernel *xclk)
{
  NCM_UNUSED (xclk);
  
  return 1;
}

static guint
_nc_xcor_limber_kernel_tsz_obs_params_len (NcXcorLimberKernel *xclk)
{
  NCM_UNUSED (xclk);
  
  return 0;
}

