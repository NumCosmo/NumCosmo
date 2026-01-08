/***************************************************************************
 *            nc_xcor_kernel_weak_lensing.c
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
 * NcXcorKernelWeakLensing:
 *
 * Implementation of #NcXcorKernel for galaxy weak lensing.
 *
 * The kernel is given by:
 * \begin{equation}
 *    W^{\kappa_{\mathrm{gal}}} = \frac{3}{2} \frac{\Omega_m H_0^2}{c} \frac{(1+z)}{H(z)} \chi(z) \int dz^\prime \frac{dn}{dz^\prime} \frac{\chi(z^\prime) - \chi(z)}{\chi(z^\prime)}
 * \end{equation}
 * where $\frac{dn}{dz}$ is the redshift distribution of galaxies.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_cfg.h"
#include "xcor/nc_xcor.h"
#include "xcor/nc_xcor_kernel_weak_lensing.h"
#include "xcor/nc_xcor_lensing_efficiency.h"

#include "math/ncm_integrate.h"
#include "math/ncm_spline_gsl.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorKernelWeakLensing
{
  /*< private >*/
  NcXcorKernel parent_instance;

  /* Lensing efficiency */
  NcXcorLensingEfficiency *lens_eff;

  /* Redshift */
  NcmSpline *dn_dz;
  gdouble dn_dz_zmin;
  gdouble dn_dz_zmax;
  gdouble dn_dz_min;
  gdouble dn_dz_max;

  gdouble nbar;
  gdouble intr_shear;

  NcDistance *dist;
  NcmPowspec *ps;

  gdouble noise;
  NcXcorKernelIntegMethod integ_method;
};

enum
{
  PROP_0,
  PROP_DN_DZ,
  PROP_NBAR,
  PROP_INTR_SHEAR,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelWeakLensing, nc_xcor_kernel_weak_lensing, NC_TYPE_XCOR_KERNEL)

#define VECTOR     (NCM_MODEL (xclkg))
#define MAG_BIAS   (ncm_model_orig_param_get (VECTOR, NC_XCOR_KERNEL_WEAK_LENSING_MAG_BIAS))
#define NOISE_BIAS (ncm_model_orig_param_get (VECTOR, NC_XCOR_KERNEL_WEAK_LENSING_NOISE_BIAS))

static void
nc_xcor_kernel_weak_lensing_init (NcXcorKernelWeakLensing *xclkg)
{
  xclkg->lens_eff     = NULL;
  xclkg->dn_dz        = NULL;
  xclkg->dn_dz_zmin   = 0.0;
  xclkg->dn_dz_zmax   = 0.0;
  xclkg->dn_dz_min    = 0.0;
  xclkg->dn_dz_max    = 0.0;
  xclkg->nbar         = 0.0;
  xclkg->intr_shear   = 0.0;
  xclkg->dist         = NULL;
  xclkg->ps           = NULL;
  xclkg->noise        = 0.0;
  xclkg->integ_method = NC_XCOR_KERNEL_INTEG_METHOD_LEN;
}

static void
_nc_xcor_kernel_weak_lensing_take_dndz (NcXcorKernelWeakLensing *xclkg, NcmSpline *dn_dz)
{
  NcmVector *z_vec      = ncm_spline_peek_xv (dn_dz);
  const guint z_vec_len = ncm_vector_len (z_vec);

  xclkg->dn_dz      = dn_dz;
  xclkg->dn_dz_zmin = ncm_vector_get (z_vec, 0);
  xclkg->dn_dz_zmax = ncm_vector_get (z_vec, z_vec_len - 1);
}

static void
_nc_xcor_kernel_weak_lensing_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_WEAK_LENSING (object));

  switch (prop_id)
  {
    case PROP_DN_DZ:
      _nc_xcor_kernel_weak_lensing_take_dndz (xclkg, g_value_dup_object (value));
      break;
    case PROP_NBAR:
      xclkg->nbar = g_value_get_double (value);
      break;
    case PROP_INTR_SHEAR:
      xclkg->intr_shear = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_kernel_weak_lensing_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_WEAK_LENSING (object));

  switch (prop_id)
  {
    case PROP_DN_DZ:
      g_value_set_object (value, xclkg->dn_dz);
      break;
    case PROP_NBAR:
      g_value_set_double (value, xclkg->nbar);
      break;
    case PROP_INTR_SHEAR:
      g_value_set_double (value, xclkg->intr_shear);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static gdouble _nc_xcor_kernel_weak_lensing_lens_eff_eval_source (NcXcorLensingEfficiency *lens_eff, gdouble z);
static void _nc_xcor_kernel_weak_lensing_lens_eff_get_z_range (NcXcorLensingEfficiency *lens_eff, gdouble *zmin, gdouble *zmax);

NC_XCOR_LENSING_EFFICIENCY_DEFINE_TYPE (NC, XCOR_KERNEL_WEAK_LENSING_LENS_EFF,
                                        NcXcorKernelWeakLensingLensEff,
                                        nc_xcor_kernel_weak_lensing_lens_eff,
                                        _nc_xcor_kernel_weak_lensing_lens_eff_eval_source,
                                        _nc_xcor_kernel_weak_lensing_lens_eff_get_z_range,
                                        NcXcorKernelWeakLensing *)

static gdouble _nc_xcor_kernel_weak_lensing_eval_kernel_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble k, gint l);
static gdouble _nc_xcor_kernel_weak_lensing_eval_kernel_limber_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
static void _nc_xcor_kernel_weak_lensing_get_k_range_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l, gdouble *kmin, gdouble *kmax);

static void
_nc_xcor_kernel_weak_lensing_constructed (GObject *object)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (object);
  NcXcorKernel *xclk             = NC_XCOR_KERNEL (xclkg);
  NcXcorKernelIntegMethod integ  = nc_xcor_kernel_get_integ_method (xclk);

  switch (integ)
  {
    case NC_XCOR_KERNEL_INTEG_METHOD_GSL_QAG:
      g_error ("_nc_xcor_kernel_weak_lensing_constructed: GSL_QAG integration not implemented yet");
      break;
    case NC_XCOR_KERNEL_INTEG_METHOD_LIMBER:
      nc_xcor_kernel_set_eval_kernel_func (xclk,
                                           _nc_xcor_kernel_weak_lensing_eval_kernel_limber,
                                           _nc_xcor_kernel_weak_lensing_eval_kernel_limber_prefactor
                                          );
      nc_xcor_kernel_set_get_k_range_func (xclk, _nc_xcor_kernel_weak_lensing_get_k_range_limber);
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  /* Chain up : middle */
  G_OBJECT_CLASS (nc_xcor_kernel_weak_lensing_parent_class)->constructed (object);
  {
    NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (object);
    NcXcorKernel *xclk             = NC_XCOR_KERNEL (xclkg);
    gdouble zmin, zmax, zmid;

    nc_xcor_kernel_set_z_range (xclk, 0.0, xclkg->dn_dz_zmax, 0.5 * xclkg->dn_dz_zmax);
    nc_xcor_kernel_get_z_range (xclk, &zmin, &zmax, &zmid);

    ncm_spline_prepare (xclkg->dn_dz);
    /* Normalize the redshift distribution */
    {
      gdouble ngal  = ncm_spline_eval_integ (xclkg->dn_dz, xclkg->dn_dz_zmin, xclkg->dn_dz_zmax);
      NcmVector *yv = ncm_spline_peek_yv (xclkg->dn_dz);

      ncm_vector_scale (yv, 1.0 / ngal);
      ncm_spline_prepare (xclkg->dn_dz);

      xclkg->dn_dz_min = ncm_spline_eval (xclkg->dn_dz, xclkg->dn_dz_zmin);
      xclkg->dn_dz_max = ncm_spline_eval (xclkg->dn_dz, xclkg->dn_dz_zmax);
    }

    {
      NcDistance *dist                             = nc_xcor_kernel_peek_dist (xclk);
      NcXcorKernelWeakLensingLensEff *lens_eff_obj = g_object_new (
        nc_xcor_kernel_weak_lensing_lens_eff_get_type (),
        "distance", dist,
        NULL);

      lens_eff_obj->data = xclkg;
      xclkg->lens_eff    = NC_XCOR_LENSING_EFFICIENCY (lens_eff_obj);
    }

    /* Noise level */
    xclkg->noise = gsl_pow_2 (xclkg->intr_shear) / xclkg->nbar;
  }
}

static void
_nc_xcor_kernel_weak_lensing_dispose (GObject *object)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (object);

  ncm_spline_clear (&xclkg->dn_dz);
  nc_xcor_lensing_efficiency_clear (&xclkg->lens_eff);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_weak_lensing_parent_class)->dispose (object);
}

static void
_nc_xcor_kernel_weak_lensing_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_weak_lensing_parent_class)->finalize (object);
}

static gdouble _nc_xcor_kernel_weak_lensing_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
static gdouble _nc_xcor_kernel_weak_lensing_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
static void _nc_xcor_kernel_weak_lensing_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo);

/*static gdouble _nc_xcor_kernel_weak_lensing_noise_spec (NcXcorKernel* xclk, guint l);*/
static void _nc_xcor_kernel_weak_lensing_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
guint _nc_xcor_kernel_weak_lensing_obs_len (NcXcorKernel *xclk);
guint _nc_xcor_kernel_weak_lensing_obs_params_len (NcXcorKernel *xclk);

static void
nc_xcor_kernel_weak_lensing_class_init (NcXcorKernelWeakLensingClass *klass)
{
  GObjectClass *object_class      = G_OBJECT_CLASS (klass);
  NcXcorKernelClass *parent_class = NC_XCOR_KERNEL_CLASS (klass);
  NcmModelClass *model_class      = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_xcor_kernel_weak_lensing_set_property;
  model_class->get_property = &_nc_xcor_kernel_weak_lensing_get_property;
  object_class->constructed = &_nc_xcor_kernel_weak_lensing_constructed;
  object_class->dispose     = &_nc_xcor_kernel_weak_lensing_dispose;
  object_class->finalize    = &_nc_xcor_kernel_weak_lensing_finalize;

  ncm_model_class_set_name_nick (model_class, "Xcor weak lensing", "Xcor-WL");
  ncm_model_class_add_params (model_class, NC_XCOR_KERNEL_WEAK_LENSING_SPARAM_LEN, NC_XCOR_KERNEL_WEAK_LENSING_VPARAM_LEN, PROP_SIZE);

  g_object_class_install_property (object_class,
                                   PROP_DN_DZ,
                                   g_param_spec_object ("dndz",
                                                        NULL,
                                                        "Source redshift distribution",
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_CONSTRUCT | G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NBAR,
                                   g_param_spec_double ("nbar",
                                                        NULL,
                                                        "nbar (galaxy angular density)",
                                                        0.0, GSL_POSINF, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_INTR_SHEAR,
                                   g_param_spec_double ("intr-shear",
                                                        NULL,
                                                        "Intrinsic galaxy shear",
                                                        0.0, 20.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->eval_limber_z           = &_nc_xcor_kernel_weak_lensing_eval_limber_z;
  parent_class->eval_limber_z_prefactor = &_nc_xcor_kernel_weak_lensing_eval_limber_z_prefactor;
  parent_class->prepare                 = &_nc_xcor_kernel_weak_lensing_prepare;
  parent_class->add_noise               = &_nc_xcor_kernel_weak_lensing_add_noise;

  parent_class->obs_len        = &_nc_xcor_kernel_weak_lensing_obs_len;
  parent_class->obs_params_len = &_nc_xcor_kernel_weak_lensing_obs_params_len;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

static gdouble
_nc_xcor_kernel_weak_lensing_lens_eff_eval_source (NcXcorLensingEfficiency *lens_eff, gdouble z)
{
  NcXcorKernelWeakLensingLensEff *data_obj = NC_XCOR_KERNEL_WEAK_LENSING_LENS_EFF (lens_eff);
  NcXcorKernelWeakLensing *xclkg           = data_obj->data;
  const gdouble alpha                      = 1.0e-2;

  if (z < xclkg->dn_dz_zmin)
    return xclkg->dn_dz_min * exp (-gsl_pow_2 ((z - xclkg->dn_dz_zmin) / alpha));

  if (z > xclkg->dn_dz_zmax)
    return xclkg->dn_dz_max * exp (-gsl_pow_2 ((z - xclkg->dn_dz_zmax) / alpha));

  return ncm_spline_eval (xclkg->dn_dz, z);
}

static void
_nc_xcor_kernel_weak_lensing_lens_eff_get_z_range (NcXcorLensingEfficiency *lens_eff, gdouble *zmin, gdouble *zmax)
{
  NcXcorKernelWeakLensingLensEff *data_obj = NC_XCOR_KERNEL_WEAK_LENSING_LENS_EFF (lens_eff);
  NcXcorKernelWeakLensing *xclkg           = data_obj->data;

  *zmin = xclkg->dn_dz_zmin;
  *zmax = xclkg->dn_dz_zmax;
}

static void
_nc_xcor_kernel_weak_lensing_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (xclk);
  NcDistance *dist               = nc_xcor_kernel_peek_dist (xclk);
  NcmPowspec *ps                 = nc_xcor_kernel_peek_powspec (xclk);
  gdouble zmin, zmax, zmid;

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  xclkg->dist = dist;
  xclkg->ps   = ps;

  nc_xcor_kernel_get_z_range (xclk, &zmin, &zmax, &zmid);

  g_assert_cmpfloat (zmin, ==, 0.0);
  g_assert (gsl_finite (ncm_spline_eval (xclkg->dn_dz, zmin)));
  g_assert (gsl_finite (ncm_spline_eval (xclkg->dn_dz, zmax)));

  zmid = ncm_vector_get (ncm_spline_get_xv (xclkg->dn_dz), ncm_vector_get_max_index (ncm_spline_get_yv (xclkg->dn_dz))) / 2.0;
  nc_xcor_kernel_set_z_range (xclk, zmin, zmax, zmid);

  nc_xcor_lensing_efficiency_prepare (xclkg->lens_eff, cosmo);
}

static gdouble
_nc_xcor_kernel_weak_lensing_eval_radial_weight (NcXcorKernel *xclk, NcHICosmo *cosmo, const gdouble z, const gdouble xi, const gdouble E)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (xclk);

  return (1.0 + z) / xi * nc_xcor_lensing_efficiency_eval (xclkg->lens_eff, z);
}

static gdouble
_nc_xcor_kernel_weak_lensing_eval_kernel_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble k, gint l)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (xclk);
  const gdouble nu               = l + 0.5;
  const gdouble xi_nu            = nu / k;
  const gdouble z                = nc_distance_inv_comoving (xclkg->dist, cosmo, xi_nu);
  const gdouble E_z              = nc_hicosmo_E (cosmo, z);
  const gdouble powspec          = ncm_powspec_eval (xclkg->ps, NCM_MODEL (cosmo), z, k / nc_hicosmo_RH_Mpc (cosmo));
  const gdouble kernel           = _nc_xcor_kernel_weak_lensing_eval_radial_weight (xclk, cosmo, z, xi_nu, E_z);
  const gdouble operator_limber  = 1.0 / gsl_pow_3 (k);

  return operator_limber * kernel * sqrt (powspec);
}

static gdouble
_nc_xcor_kernel_weak_lensing_eval_kernel_limber_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  const gdouble cosmo_factor = 1.5 * nc_hicosmo_Omega_m0 (cosmo);
  const gdouble nu           = l + 0.5;
  const gdouble lfactor      = sqrt ((l + 2.0) * (l + 1.0) * l * (l - 1.0));

  return sqrt (M_PI / 2.0 / nu) * cosmo_factor * lfactor;
}

static void
_nc_xcor_kernel_weak_lensing_get_k_range_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l, gdouble *kmin, gdouble *kmax)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (xclk);
  NcDistance *dist               = nc_xcor_kernel_peek_dist (xclk);
  NcmPowspec *ps                 = nc_xcor_kernel_peek_powspec (xclk);
  const gdouble ps_kmin          = ncm_powspec_get_kmin (ps) * nc_hicosmo_RH_Mpc (cosmo);
  const gdouble ps_kmax          = ncm_powspec_get_kmax (ps) * nc_hicosmo_RH_Mpc (cosmo);
  const gdouble nu               = l + 0.5;
  const gdouble xi_max           = nc_distance_comoving (dist, cosmo, xclkg->dn_dz_zmax);
  const gdouble kmin_limber      = nu / xi_max;

  *kmin = GSL_MAX (ps_kmin, kmin_limber);
  *kmax = ps_kmax;
}

static gdouble
_nc_xcor_kernel_weak_lensing_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  return gsl_pow_2 (xck->xi_z) * _nc_xcor_kernel_weak_lensing_eval_radial_weight (xclk, cosmo, z, xck->xi_z, xck->E_z);
}

static gdouble
_nc_xcor_kernel_weak_lensing_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  const gdouble lfactor = sqrt ((l + 2.0) * (l + 1.0) * l * (l - 1.0));
  const gdouble nu      = l + 0.5;

  return 1.5 * nc_hicosmo_Omega_m0 (cosmo) * lfactor / (nu * nu);
}

static void
_nc_xcor_kernel_weak_lensing_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  NcXcorKernelWeakLensing *xclkg = NC_XCOR_KERNEL_WEAK_LENSING (xclk);

  /* ncm_vector_add_constant (vp1, NOISE_BIAS); / * take noise_bias into account * / */
  ncm_vector_memcpy (vp2, vp1);
  ncm_vector_add_constant (vp2, xclkg->noise);
}

guint
_nc_xcor_kernel_weak_lensing_obs_len (NcXcorKernel *xclk)
{
  return 2;
}

guint
_nc_xcor_kernel_weak_lensing_obs_params_len (NcXcorKernel *xclk)
{
  return 1;
}

/**
 * nc_xcor_kernel_weak_lensing_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @zmin: a gdouble
 * @zmax: a gdouble
 * @dn_dz: a #NcmSpline
 * @nbar: a gdouble, gal density
 * @intr_shear: a gdouble, intrinsic galaxy shear
 *
 * Returns: a #NcXcorKernelWeakLensing
 */
NcXcorKernelWeakLensing *
nc_xcor_kernel_weak_lensing_new (NcDistance *dist, NcmPowspec *ps, gdouble zmin, gdouble zmax, NcmSpline *dn_dz, gdouble nbar, gdouble intr_shear)
{
  NcXcorKernelWeakLensing *xclkg = g_object_new (NC_TYPE_XCOR_KERNEL_WEAK_LENSING,
                                                 "dist", dist,
                                                 "powspec", ps,
                                                 "zmin", zmin,
                                                 "zmax", zmax,
                                                 "dndz", dn_dz,
                                                 "nbar", nbar,
                                                 "intr-shear", intr_shear,
                                                 NULL);

  return xclkg;
}

