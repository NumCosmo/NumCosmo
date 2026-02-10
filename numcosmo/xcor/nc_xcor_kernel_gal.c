/***************************************************************************
 *            nc_xcor_kernel_gal.c
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
 * NcXcorKernelGal:
 *
 * Implementation of #NcXcorKernel for galaxy density.
 *
 * The kernel is given by
 * \begin{equation}
 *     W^g (z) = b(z) \frac{dn}{dz} + \frac{3\Omega_m}{2c} \frac{H^2_0}{H(z)} (1+ z) \, \chi(z) \, (5s-2) \, g(z)
 * \end{equation}
 * where
 * \begin{equation}
 * g(z) = \int_z^{z_{*}} dz^\prime \left( 1 - \frac{\chi(z)}{\chi(z^\prime)}\right) \frac{dn}{dz^\prime}$.
 * \end{equation}
 *
 * $\frac{dn}{dz}$ is the (automatically normalized) redshift distribution of galaxies (implemented as a #NcmSpline), $ b(z)$ is the large-scale clustering bias and $s$ is the magnification bias. $b(z)$ can be a single number or a #NcmSpline.
 *
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_cfg.h"
#include "xcor/nc_xcor_kernel_gal.h"
#include "xcor/nc_xcor_lensing_efficiency.h"
#include "xcor/nc_xcor.h"

#include "math/ncm_integrate.h"
#include "math/ncm_spline_gsl.h"
#include "math/ncm_spline_func.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorKernelGal
{
  /*< private >*/
  NcXcorKernel parent_instance;

  NcmSpline *dn_dz;
  gdouble dn_dz_zmin;
  gdouble dn_dz_zmax;
  gdouble dn_dz_min;
  gdouble dn_dz_max;

  NcmSpline *bias_spline;
  guint nknots;
  gdouble *bias;

  NcXcorLensingEfficiency *lens_eff;
  gboolean domagbias;

  gboolean fast_update;
  gdouble bias_old;
  gdouble noise_bias_old;

  gdouble nbarm1;

  NcDistance *dist;
  NcmPowspec *ps;
};

enum
{
  PROP_0,
  PROP_DN_DZ,
  PROP_BIAS,
  PROP_DOMAGBIAS,
  PROP_NBARM1,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorKernelGal, nc_xcor_kernel_gal, NC_TYPE_XCOR_KERNEL)

#define VECTOR     (NCM_MODEL (xclkg))
#define MAG_BIAS   (ncm_model_orig_param_get (VECTOR, NC_XCOR_KERNEL_GAL_MAG_BIAS))
#define NOISE_BIAS (ncm_model_orig_param_get (VECTOR, NC_XCOR_KERNEL_GAL_NOISE_BIAS))

static void
nc_xcor_kernel_gal_init (NcXcorKernelGal *xclkg)
{
  xclkg->dn_dz      = NULL;
  xclkg->dn_dz_zmin = 0.0;
  xclkg->dn_dz_zmax = 0.0;
  xclkg->dn_dz_min  = 0.0;
  xclkg->dn_dz_max  = 0.0;

  xclkg->bias_spline = NULL;
  xclkg->nknots      = 0;
  xclkg->bias        = NULL;

  xclkg->lens_eff  = NULL;
  xclkg->domagbias = FALSE;

  xclkg->fast_update    = FALSE;
  xclkg->bias_old       = 0.0;
  xclkg->noise_bias_old = 0.0;

  xclkg->nbarm1 = 0.0;
  xclkg->dist   = NULL;
  xclkg->ps     = NULL;
}

static void
_nc_xcor_kernel_gal_take_dndz (NcXcorKernelGal *xclkg, NcmSpline *dn_dz)
{
  NcmVector *z_vec      = ncm_spline_peek_xv (dn_dz);
  const guint z_vec_len = ncm_vector_len (z_vec);

  xclkg->dn_dz      = dn_dz;
  xclkg->dn_dz_zmin = ncm_vector_get (z_vec, 0);
  xclkg->dn_dz_zmax = ncm_vector_get (z_vec, z_vec_len - 1);
}

static void
_nc_xcor_kernel_gal_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorKernelGal *xclkg = NC_XCOR_KERNEL_GAL (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_GAL (object));

  switch (prop_id)
  {
    case PROP_BIAS:
      xclkg->bias_spline = g_value_dup_object (value);
      break;
    case PROP_DOMAGBIAS:
      xclkg->domagbias = g_value_get_boolean (value);
      break;
    case PROP_NBARM1:
      xclkg->nbarm1 = g_value_get_double (value);
      break;
    case PROP_DN_DZ:
      _nc_xcor_kernel_gal_take_dndz (xclkg, g_value_dup_object (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_kernel_gal_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorKernelGal *xclkg = NC_XCOR_KERNEL_GAL (object);

  g_return_if_fail (NC_IS_XCOR_KERNEL_GAL (object));

  switch (prop_id)
  {
    case PROP_BIAS:
      g_value_set_object (value, xclkg->bias_spline);
      break;
    case PROP_DOMAGBIAS:
      g_value_set_boolean (value, xclkg->domagbias);
      break;
    case PROP_NBARM1:
      g_value_set_double (value, xclkg->nbarm1);
      break;
    case PROP_DN_DZ:
      g_value_set_object (value, xclkg->dn_dz);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static gdouble _nc_xcor_kernel_gal_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
static gdouble _nc_xcor_kernel_gal_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
static gdouble _nc_xcor_kernel_gal_eval_kernel_prefactor_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
static gdouble _nc_xcor_kernel_gal_eval_kernel_magbias_prefactor_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);
static void _nc_xcor_kernel_gal_get_k_range_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l, gdouble *kmin, gdouble *kmax);
static void _nc_xcor_kernel_gal_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo);
static void _nc_xcor_kernel_gal_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
static guint _nc_xcor_kernel_gal_obs_len (NcXcorKernel *xclk);
static guint _nc_xcor_kernel_gal_obs_params_len (NcXcorKernel *xclk);
static void _nc_xcor_kernel_gal_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid);
static gdouble _nc_xcor_kernel_gal_dndz (NcXcorKernelGal *xclkg, gdouble z);
static NcXcorKernelIntegrand *_nc_xcor_kernel_gal_get_eval (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l);

static gdouble _nc_xcor_kernel_gal_lens_eff_eval_source (NcXcorLensingEfficiency *lens_eff, gdouble z);
static void _nc_xcor_kernel_gal_lens_eff_get_z_range (NcXcorLensingEfficiency *lens_eff, gdouble *zmin, gdouble *zmax);

NC_XCOR_LENSING_EFFICIENCY_DEFINE_TYPE (NC, XCOR_KERNEL_GAL_LENS_EFF,
                                        NcXcorKernelGalLensEff,
                                        nc_xcor_kernel_gal_lens_eff,
                                        _nc_xcor_kernel_gal_lens_eff_eval_source,
                                        _nc_xcor_kernel_gal_lens_eff_get_z_range,
                                        NcXcorKernelGal *)

static void
_nc_xcor_kernel_gal_constructed (GObject *object)
{
  NcXcorKernelGal *xclkg = NC_XCOR_KERNEL_GAL (object);
  NcXcorKernel *xclk     = NC_XCOR_KERNEL (xclkg);

  /* Chain up : middle */
  G_OBJECT_CLASS (nc_xcor_kernel_gal_parent_class)->constructed (object);
  {
    NcmModel *model = NCM_MODEL (xclkg);
    guint i;

    /* Initialize g function for magnification bias */
    if (xclkg->domagbias)
    {
      NcDistance *dist = nc_xcor_kernel_peek_dist (xclk);
      NcXcorKernelGalLensEff *lens_eff_obj;

      g_clear_object (&xclkg->lens_eff);

      lens_eff_obj = g_object_new (nc_xcor_kernel_gal_lens_eff_get_type (),
                                   "distance", dist,
                                   NULL);
      lens_eff_obj->data = xclkg;
      xclkg->lens_eff    = NC_XCOR_LENSING_EFFICIENCY (lens_eff_obj);
    }

    /* Normalize the redshift distribution */
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

    /* Prepare the bias spline and link it to the bias parameter vector */
    NcmVector *zv, *bv;
    guint bvi;
    guint bz_size = ncm_model_vparam_len (model, NC_XCOR_KERNEL_GAL_BIAS);

    xclkg->nknots = bz_size;

    {
      NcmVector *orig_vec = ncm_model_orig_params_peek_vector (model);

      bvi = ncm_model_vparam_index (model, NC_XCOR_KERNEL_GAL_BIAS, 0);
      zv  = ncm_vector_new (bz_size);
      bv  = ncm_vector_get_subvector (orig_vec, bvi, bz_size);
    }

    {
      const gdouble zmin = 0.0;
      const gdouble zmax = xclkg->dn_dz_zmax;

      switch (bz_size)
      {
        case 1:
          xclkg->bias_spline = NULL;
          xclkg->bias        = ncm_vector_ptr (bv, 0);
          break;
        case 2:
          ncm_vector_set (zv, 0, zmin);
          ncm_vector_set (zv, 1, zmax);
          xclkg->bias_spline = NCM_SPLINE (ncm_spline_gsl_new_full (gsl_interp_linear, zv, bv, FALSE));
          break;
        default:
        {
          for (i = 0; i < bz_size; i++)
          {
            gdouble zi = zmin + (zmax - zmin) / (bz_size - 1) * i;

            ncm_vector_set (zv, i, zi);
          }

          xclkg->bias_spline = NCM_SPLINE (ncm_spline_gsl_new_full (gsl_interp_polynomial, zv, bv, FALSE));
          break;
        }
      }
    }

    ncm_vector_free (bv);
    ncm_vector_free (zv);

    /* If bias is constant, it's just a multiplicative factor, so possible to accelerate recomputation of C_l's */
    xclkg->fast_update = (bz_size == 1) && !(xclkg->domagbias);
  }
}

static void
_nc_xcor_kernel_gal_dispose (GObject *object)
{
  NcXcorKernelGal *xclkg = NC_XCOR_KERNEL_GAL (object);

  ncm_spline_clear (&xclkg->bias_spline);
  ncm_spline_clear (&xclkg->dn_dz);
  nc_xcor_lensing_efficiency_clear (&xclkg->lens_eff);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_gal_parent_class)->dispose (object);
}

static void
_nc_xcor_kernel_gal_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_kernel_gal_parent_class)->finalize (object);
}

static void
nc_xcor_kernel_gal_class_init (NcXcorKernelGalClass *klass)
{
  GObjectClass *object_class      = G_OBJECT_CLASS (klass);
  NcXcorKernelClass *parent_class = NC_XCOR_KERNEL_CLASS (klass);
  NcmModelClass *model_class      = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_xcor_kernel_gal_set_property;
  model_class->get_property = &_nc_xcor_kernel_gal_get_property;
  object_class->constructed = &_nc_xcor_kernel_gal_constructed;
  object_class->dispose     = &_nc_xcor_kernel_gal_dispose;
  object_class->finalize    = &_nc_xcor_kernel_gal_finalize;

  ncm_model_class_set_name_nick (model_class, "Xcor galaxy distribution", "Xcor-gal");
  ncm_model_class_add_params (model_class, NC_XCOR_KERNEL_GAL_SPARAM_LEN, NC_XCOR_KERNEL_GAL_VPARAM_LEN, PROP_SIZE);

  g_object_class_install_property (object_class,
                                   PROP_DN_DZ,
                                   g_param_spec_object ("dndz",
                                                        NULL,
                                                        "Galaxy redshift distribution",
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_CONSTRUCT | G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BIAS,
                                   g_param_spec_object ("bias",
                                                        NULL,
                                                        "Bias spline object",
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_DOMAGBIAS,
                                   g_param_spec_boolean ("domagbias",
                                                         NULL,
                                                         "Do magnification bias",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_NBARM1,
                                   g_param_spec_double ("nbarm1",
                                                        NULL,
                                                        "One over nbar (galaxy angular density)",
                                                        0.0, 20.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /*
   * Distribution's magnification bias: mag_bias.
   * FIXME Set correct values (limits)
   */
  ncm_model_class_set_sparam (model_class, NC_XCOR_KERNEL_GAL_MAG_BIAS,
                              "mag_bias", "mag_bias",
                              -10., 10., 0.1,
                              NC_XCOR_KERNEL_GAL_DEFAULT_PARAMS_ABSTOL, NC_XCOR_KERNEL_GAL_DEFAULT_MAG_BIAS, NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_XCOR_KERNEL_GAL_NOISE_BIAS,
                              "noise_bias", "noise_bias",
                              -1., 1., 1e-6,
                              NC_XCOR_KERNEL_GAL_DEFAULT_PARAMS_ABSTOL, NC_XCOR_KERNEL_GAL_DEFAULT_NOISE_BIAS, NCM_PARAM_TYPE_FIXED);


  ncm_model_class_set_vparam (model_class, NC_XCOR_KERNEL_GAL_BIAS, NC_XCOR_KERNEL_GAL_BIAS_DEFAULT_LEN,
                              "bparam", "bparam",
                              0.0, 10.0, 0.1,
                              NC_XCOR_KERNEL_GAL_DEFAULT_PARAMS_ABSTOL, NC_XCOR_KERNEL_GAL_DEFAULT_BIAS,
                              NCM_PARAM_TYPE_FREE);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->eval_limber_z           = &_nc_xcor_kernel_gal_eval_limber_z;
  parent_class->eval_limber_z_prefactor = &_nc_xcor_kernel_gal_eval_limber_z_prefactor;
  parent_class->prepare                 = &_nc_xcor_kernel_gal_prepare;
  parent_class->add_noise               = &_nc_xcor_kernel_gal_add_noise;

  parent_class->obs_len        = &_nc_xcor_kernel_gal_obs_len;
  parent_class->obs_params_len = &_nc_xcor_kernel_gal_obs_params_len;
  parent_class->get_z_range    = &_nc_xcor_kernel_gal_get_z_range;
  parent_class->get_k_range    = &_nc_xcor_kernel_gal_get_k_range_limber;
  parent_class->get_eval       = &_nc_xcor_kernel_gal_get_eval;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_KERNEL_IMPL_ALL);
}

static gdouble
_nc_xcor_kernel_gal_lens_eff_eval_source (NcXcorLensingEfficiency *lens_eff, gdouble z)
{
  NcXcorKernelGalLensEff *data_obj = NC_XCOR_KERNEL_GAL_LENS_EFF (lens_eff);
  NcXcorKernelGal *xclkg           = data_obj->data;
  const gdouble mag_bias           = MAG_BIAS;

  return (5.0 * mag_bias - 2.0) * _nc_xcor_kernel_gal_dndz (xclkg, z);
}

static void
_nc_xcor_kernel_gal_lens_eff_get_z_range (NcXcorLensingEfficiency *lens_eff, gdouble *zmin, gdouble *zmax)
{
  NcXcorKernelGalLensEff *data_obj = NC_XCOR_KERNEL_GAL_LENS_EFF (lens_eff);
  NcXcorKernelGal *xclkg           = data_obj->data;

  *zmin = xclkg->dn_dz_zmin;
  *zmax = xclkg->dn_dz_zmax;
}

/**
 * nc_xcor_kernel_gal_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @np: number of points in the interpolation
 * @nbarm1: a gdouble, noise spectrum
 * @dn_dz: a #NcmSpline
 * @domagbias: whether to do magnification bias
 *
 * Returns: a #NcXcorKernelGal
 */
NcXcorKernelGal *
nc_xcor_kernel_gal_new (NcDistance *dist, NcmPowspec *ps, gsize np, gdouble nbarm1, NcmSpline *dn_dz, gboolean domagbias)
{
  NcXcorKernelGal *xclkg = g_object_new (NC_TYPE_XCOR_KERNEL_GAL,
                                         "dist", dist,
                                         "powspec", ps,
                                         "bparam-length", np,
                                         "nbarm1", nbarm1,
                                         "dndz", dn_dz,
                                         "domagbias", domagbias,
                                         NULL);

  return xclkg;
}

static gdouble
_nc_xcor_kernel_gal_dndz (NcXcorKernelGal *xclkg, gdouble z)
{
  const gdouble alpha = 1.0e-2;

  if (z < xclkg->dn_dz_zmin)
    return xclkg->dn_dz_min * exp (-gsl_pow_2 ((z - xclkg->dn_dz_zmin) / alpha));

  if (z > xclkg->dn_dz_zmax)
    return xclkg->dn_dz_max * exp (-gsl_pow_2 ((z - xclkg->dn_dz_zmax) / alpha));

  return ncm_spline_eval (xclkg->dn_dz, z);
}

static void
_nc_xcor_kernel_gal_prepare (NcXcorKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorKernelGal *xclkg = NC_XCOR_KERNEL_GAL (xclk);
  NcmModel *model        = NCM_MODEL (xclk);
  NcDistance *dist       = nc_xcor_kernel_peek_dist (xclk);
  NcmPowspec *ps         = nc_xcor_kernel_peek_powspec (xclk);

  xclkg->dist = dist;
  xclkg->ps   = ps;

  nc_distance_prepare_if_needed (dist, cosmo);
  ncm_powspec_prepare_if_needed (ps, NCM_MODEL (cosmo));

  if (xclkg->fast_update)
  {
    xclkg->bias_old       = *(xclkg->bias);
    xclkg->noise_bias_old = NOISE_BIAS;
  }

  if (!ncm_model_state_is_update (model))
  {
    if (xclkg->nknots > 1)
      ncm_spline_prepare (xclkg->bias_spline);

    ncm_model_state_set_update (model);
  }

  if (xclkg->domagbias)
    nc_xcor_lensing_efficiency_prepare (xclkg->lens_eff, cosmo);
}

static gdouble
_nc_xcor_kernel_gal_bias (NcXcorKernelGal *xclkg, gdouble z)
{
  gdouble res;

  switch (xclkg->nknots)
  {
    case 1:
      res = *(xclkg->bias);
      break;
    default:
      res = ncm_spline_eval (xclkg->bias_spline, z);
      break;
  }

  return gsl_finite (res) ? res : 0.0;
}

static gdouble
_nc_xcor_kernel_gal_eval_limber_z (NcXcorKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  NcXcorKernelGal *xclkg = NC_XCOR_KERNEL_GAL (xclk);
  const gdouble dn_dz_z  = _nc_xcor_kernel_gal_dndz (xclkg, z);
  const gdouble bias_z   = _nc_xcor_kernel_gal_bias (xclkg, z);
  gdouble res            = bias_z * dn_dz_z * xck->E_z;

  if (xclkg->domagbias)
  {
    const gdouble nu   = l + 0.5;
    const gdouble llp1 = l * (l + 1.0) / nu / nu;
    const gdouble g_z  = nc_xcor_lensing_efficiency_eval (xclkg->lens_eff, z) *
                         (1.0 + z) *
                         xck->xi_z *
                         1.5 * nc_hicosmo_Omega_m0 (cosmo);

    res += llp1 * g_z;
  }

  return res;
}

static gdouble
_nc_xcor_kernel_gal_eval_limber_z_prefactor (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  return 1.0;
}

/*
 * Limber integrand callback.
 */

typedef struct _IntegData
{
  NcXcorKernelGal *xclkg;
  NcHICosmo *cosmo;
  gdouble RH_Mpc;
  gdouble l;
  gdouble nu;
  gdouble prefactor;
  gboolean domagbias;
} IntegData;

static gdouble
_nc_xcor_kernel_gal_eval_kernel_prefactor_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  const gdouble nu = l + 0.5;

  return sqrt (M_PI / 2.0 / nu);
}

static void
_integ_data_free (gpointer data)
{
  IntegData *int_data = (IntegData *) data;

  nc_hicosmo_clear (&int_data->cosmo);
  g_free (data);
}

static gpointer
_integ_data_copy (gpointer data)
{
  IntegData *src = (IntegData *) data;
  IntegData *dst = g_new0 (IntegData, 1);

  dst->xclkg     = src->xclkg;
  dst->cosmo     = nc_hicosmo_ref (src->cosmo);
  dst->l         = src->l;
  dst->nu        = src->nu;
  dst->RH_Mpc    = src->RH_Mpc;
  dst->prefactor = src->prefactor;
  dst->domagbias = src->domagbias;

  return dst;
}

static void
_integ_data_prepare (gpointer data, NcmMSet *mset)
{
  /* Nothing to prepare */
}

static gdouble
_nc_xcor_kernel_gal_eval_limber (gpointer callback_data, const gdouble k)
{
  IntegData *int_data    = (IntegData *) callback_data;
  NcXcorKernelGal *xclkg = int_data->xclkg;
  const gdouble xi_nu    = int_data->nu / k;
  const gdouble k_Mpc    = k / int_data->RH_Mpc;
  const gdouble z        = nc_distance_inv_comoving (xclkg->dist, int_data->cosmo, xi_nu);
  const gdouble E_z      = nc_hicosmo_E (int_data->cosmo, z);
  const gdouble powspec  = ncm_powspec_eval (xclkg->ps, NCM_MODEL (int_data->cosmo), z, k_Mpc);
  const gdouble dn_dz_z  = _nc_xcor_kernel_gal_dndz (xclkg, z);
  const gdouble bias_z   = _nc_xcor_kernel_gal_bias (xclkg, z);
  const gdouble limber_k = 1.0 / k;

  if (int_data->domagbias)
  {
    const gdouble llp1       = int_data->l * (int_data->l + 1.0);
    const gdouble g_z        = nc_xcor_lensing_efficiency_eval (xclkg->lens_eff, z) * (1.0 + z) / xi_nu;
    const gdouble Omega_m0   = nc_hicosmo_Omega_m0 (int_data->cosmo);
    const gdouble operator_k = 1.0 / (k * k);

    return int_data->prefactor * limber_k * (bias_z * dn_dz_z * E_z + 1.5 * Omega_m0 * llp1 * operator_k * g_z) * sqrt (powspec);
  }
  else
  {
    return int_data->prefactor * limber_k * bias_z * dn_dz_z * E_z * sqrt (powspec);
  }
}

static NcXcorKernelIntegrand *
_nc_xcor_kernel_gal_get_eval_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  NcXcorKernelGal *xclkg       = NC_XCOR_KERNEL_GAL (xclk);
  IntegData *int_data          = g_new0 (IntegData, 1);
  NcXcorKernelIntegrand *integ = nc_xcor_kernel_integrand_new (_nc_xcor_kernel_gal_eval_limber,
                                                               _integ_data_free,
                                                               _integ_data_copy,
                                                               _integ_data_prepare,
                                                               int_data);

  int_data->xclkg     = xclkg;
  int_data->cosmo     = cosmo;
  int_data->l         = l;
  int_data->nu        = l + 0.5;
  int_data->RH_Mpc    = nc_hicosmo_RH_Mpc (cosmo);
  int_data->domagbias = xclkg->domagbias;

  if (xclkg->domagbias)
    int_data->prefactor = _nc_xcor_kernel_gal_eval_kernel_magbias_prefactor_limber (xclk, cosmo, l);
  else
    int_data->prefactor = _nc_xcor_kernel_gal_eval_kernel_prefactor_limber (xclk, cosmo, l);

  return integ;
}

/*
 * End Limber integrand callback.
 */

static NcXcorKernelIntegrand *
_nc_xcor_kernel_gal_get_eval (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  return _nc_xcor_kernel_gal_get_eval_limber (xclk, cosmo, l);
}

static gdouble
_nc_xcor_kernel_gal_eval_kernel_magbias_prefactor_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l)
{
  const gdouble nu = l + 0.5;

  return sqrt (M_PI / 2.0 / nu);
}

static void
_nc_xcor_kernel_gal_get_k_range_limber (NcXcorKernel *xclk, NcHICosmo *cosmo, gint l, gdouble *kmin, gdouble *kmax)
{
  NcXcorKernelGal *xclkg    = NC_XCOR_KERNEL_GAL (xclk);
  NcDistance *dist          = nc_xcor_kernel_peek_dist (xclk);
  NcmPowspec *ps            = nc_xcor_kernel_peek_powspec (xclk);
  const gdouble ps_kmin     = ncm_powspec_get_kmin (ps) * nc_hicosmo_RH_Mpc (cosmo);
  const gdouble ps_kmax     = ncm_powspec_get_kmax (ps) * nc_hicosmo_RH_Mpc (cosmo);
  const gdouble nu          = l + 0.5;
  const gdouble xi_max      = nc_distance_comoving (dist, cosmo, xclkg->dn_dz_zmax);
  const gdouble kmin_limber = nu / xi_max;

  *kmin = GSL_MAX (ps_kmin, kmin_limber);
  *kmax = ps_kmax;
}

static void
_nc_xcor_kernel_gal_add_noise (NcXcorKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  NcXcorKernelGal *xclkg = NC_XCOR_KERNEL_GAL (xclk);

  ncm_vector_add_constant (vp1, NOISE_BIAS); /* take noise_bias into account */
  ncm_vector_memcpy (vp2, vp1);
  ncm_vector_add_constant (vp2, xclkg->nbarm1);
}

guint
_nc_xcor_kernel_gal_obs_len (NcXcorKernel *xclk)
{
  return 2;
}

guint
_nc_xcor_kernel_gal_obs_params_len (NcXcorKernel *xclk)
{
  return 1;
}

static void
_nc_xcor_kernel_gal_get_z_range (NcXcorKernel *xclk, gdouble *zmin, gdouble *zmax, gdouble *zmid)
{
  NcXcorKernelGal *xclkg = NC_XCOR_KERNEL_GAL (xclk);

  *zmin = 0.0; /* xclkg->dn_dz_zmin; */
  *zmax = xclkg->dn_dz_zmax;
  *zmid = ncm_vector_get (ncm_spline_get_xv (xclkg->dn_dz), ncm_vector_get_max_index (ncm_spline_get_yv (xclkg->dn_dz)));
}

/**
 * nc_xcor_kernel_gal_set_fast_update:
 * @xclk: a #NcXcorKernelGal
 * @fast_update: a gboolean
 *
 * Set the fast update flag.
 *
 */
void
nc_xcor_kernel_gal_set_fast_update (NcXcorKernelGal *xclk, gboolean fast_update)
{
  xclk->fast_update = fast_update;
}

/**
 * nc_xcor_kernel_gal_get_fast_update:
 * @xclk: a #NcXcorKernelGal
 *
 * Get the fast update flag.
 *
 * Returns: a gboolean
 */
gboolean
nc_xcor_kernel_gal_get_fast_update (NcXcorKernelGal *xclk)
{
  return xclk->fast_update;
}

/**
 * nc_xcor_kernel_gal_set_bias_old:
 * @xclk: a #NcXcorKernelGal
 * @bias_old: a gdouble
 * @noise_bias_old: a gdouble
 *
 * Set the bias parameters.
 *
 */
void
nc_xcor_kernel_gal_set_bias_old (NcXcorKernelGal *xclk, gdouble bias_old, gdouble noise_bias_old)
{
  xclk->bias_old       = bias_old;
  xclk->noise_bias_old = noise_bias_old;
}

/**
 * nc_xcor_kernel_gal_get_bias:
 * @xclk: a #NcXcorKernelGal
 * @bias: (out): a gdouble
 * @bias_old: (out): a gdouble
 * @noise_bias_old: (out): a gdouble
 *
 * Get the bias parameters.
 *
 */
void
nc_xcor_kernel_gal_get_bias (NcXcorKernelGal *xclk, gdouble *bias, gdouble *bias_old, gdouble *noise_bias_old)
{
  *bias           = *(xclk->bias);
  *bias_old       = xclk->bias_old;
  *noise_bias_old = xclk->noise_bias_old;
}

