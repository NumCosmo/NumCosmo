/***************************************************************************
 *            nc_xcor_limber_kernel_weak_lensing.c
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
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
 * NcXcorLimberKernelWeakLensing:
 *
 * Implementation of #NcXcorLimberKernel for galaxy weak lensing.
 *
 * The kernel is given by:
 * \begin{equation}
 *    W^{\kappa_{\mathrm{gal}}} = \frac{3}{2} \frac{\Omega_m H_0^2}{c}
 *    \frac{(1+z)}{H(z)} \chi(z) \int dz^\prime \frac{dn}{dz^\prime} \frac{\chi(^\prime)
 *    - \chi(z)}{\chi(^\prime)}
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
#include "xcor/nc_xcor_limber_kernel_weak_lensing.h"

#include "math/ncm_integrate.h"
#include "math/ncm_spline_gsl.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <cvode/cvode.h>

#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#define SUN_DENSE_ACCESS SM_ELEMENT_D

#include <nvector/nvector_serial.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcorLimberKernelWeakLensing
{
  /*< private >*/
  NcXcorLimberKernel parent_instance;

  /* Kernel integration */
  gpointer cvode;
  SUNContext sunctx;
  N_Vector yv;
  SUNMatrix A;
  SUNLinearSolver LS;
  gdouble reltol;
  gdouble abstol;
  NcmModelCtrl *ctrl_cosmo;

  /* Redshift */
  NcmSpline *dn_dz;
  gdouble dn_dz_zmin;
  gdouble dn_dz_zmax;
  gdouble dn_dz_min;
  gdouble dn_dz_max;

  /* NcmSpline* bias_spline; */
  /* guint nknots; */
  /* gdouble* bias; */

  NcDistance *dist;

  NcmSpline *kernel_W_mz;
  /* gboolean domagbias; */

  /* gboolean fast_update; */
  /* gdouble bias_old; */
  /* gdouble noise_bias_old; */

  gdouble nbar;
  gdouble intr_shear;

  gdouble noise;
};

enum
{
  PROP_0,
  PROP_DN_DZ,
  PROP_NBAR,
  PROP_INTR_SHEAR,
  PROP_DIST,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorLimberKernelWeakLensing, nc_xcor_limber_kernel_weak_lensing, NC_TYPE_XCOR_LIMBER_KERNEL)

#define VECTOR     (NCM_MODEL (xclkg))
#define MAG_BIAS   (ncm_model_orig_param_get (VECTOR, NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_MAG_BIAS))
#define NOISE_BIAS (ncm_model_orig_param_get (VECTOR, NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_NOISE_BIAS))

static void
nc_xcor_limber_kernel_weak_lensing_init (NcXcorLimberKernelWeakLensing *xclkg)
{
  if (SUNContext_Create (SUN_COMM_NULL, &xclkg->sunctx))
    g_error ("ERROR: SUNContext_Create failed\n");

  xclkg->cvode      = NULL;
  xclkg->yv         = N_VNew_Serial (2, xclkg->sunctx);
  xclkg->A          = SUNDenseMatrix (2, 2, xclkg->sunctx);
  xclkg->LS         = SUNLinSol_Dense (xclkg->yv, xclkg->A, xclkg->sunctx);
  xclkg->reltol     = 0.0;
  xclkg->abstol     = 0.0;
  xclkg->ctrl_cosmo = ncm_model_ctrl_new (NULL);

  xclkg->dn_dz       = NULL;
  xclkg->dn_dz_zmin  = 0.0;
  xclkg->dn_dz_zmax  = 0.0;
  xclkg->dn_dz_min   = 0.0;
  xclkg->dn_dz_max   = 0.0;
  xclkg->dist        = NULL;
  xclkg->kernel_W_mz = NULL;
  xclkg->nbar        = 0.0;
  xclkg->intr_shear  = 0.0;
  xclkg->noise       = 0.0;

  NCM_CVODE_CHECK ((gpointer) xclkg->A, "SUNDenseMatrix", 0, );
  NCM_CVODE_CHECK ((gpointer) xclkg->LS, "SUNDenseLinearSolver", 0, );
}

static void
_nc_xcor_limber_kernel_weak_lensing_take_dndz (NcXcorLimberKernelWeakLensing *xclkg, NcmSpline *dn_dz)
{
  NcmVector *z_vec      = ncm_spline_peek_xv (dn_dz);
  const guint z_vec_len = ncm_vector_len (z_vec);

  xclkg->dn_dz      = dn_dz;
  xclkg->dn_dz_zmin = ncm_vector_get (z_vec, 0);
  xclkg->dn_dz_zmax = ncm_vector_get (z_vec, z_vec_len - 1);
}

static void
_nc_xcor_limber_kernel_weak_lensing_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorLimberKernelWeakLensing *xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (object);

  g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL_WEAK_LENSING (object));

  switch (prop_id)
  {
    case PROP_DN_DZ:
      _nc_xcor_limber_kernel_weak_lensing_take_dndz (xclkg, g_value_dup_object (value));
      break;
    case PROP_NBAR:
      xclkg->nbar = g_value_get_double (value);
      break;
    case PROP_INTR_SHEAR:
      xclkg->intr_shear = g_value_get_double (value);
      break;
    case PROP_DIST:
      xclkg->dist = g_value_dup_object (value);
      break;
    case PROP_RELTOL:
      xclkg->reltol = g_value_get_double (value);
      break;
    case PROP_ABSTOL:
      xclkg->abstol = g_value_get_double (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_limber_kernel_weak_lensing_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorLimberKernelWeakLensing *xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (object);

  g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL_WEAK_LENSING (object));

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
    case PROP_DIST:
      g_value_set_object (value, xclkg->dist);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, xclkg->reltol);
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, xclkg->abstol);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_limber_kernel_weak_lensing_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_limber_kernel_weak_lensing_parent_class)->constructed (object);
  {
    NcXcorLimberKernelWeakLensing *xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (object);
    NcXcorLimberKernel *xclk             = NC_XCOR_LIMBER_KERNEL (xclkg);
    gdouble zmin, zmax, zmid;

    nc_xcor_limber_kernel_set_z_range (xclk, 0.0, xclkg->dn_dz_zmax, 0.5 * xclkg->dn_dz_zmax);
    nc_xcor_limber_kernel_get_z_range (xclk, &zmin, &zmax, &zmid);

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

    /* Noise level */
    xclkg->noise = gsl_pow_2 (xclkg->intr_shear) / xclkg->nbar;
  }
}

static void
_nc_xcor_limber_kernel_weak_lensing_dispose (GObject *object)
{
  NcXcorLimberKernelWeakLensing *xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (object);

  ncm_spline_clear (&xclkg->dn_dz);
  ncm_spline_clear (&xclkg->kernel_W_mz);

  nc_distance_clear (&xclkg->dist);
  ncm_model_ctrl_clear (&xclkg->ctrl_cosmo);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_limber_kernel_weak_lensing_parent_class)->dispose (object);
}

static void
_nc_xcor_limber_kernel_weak_lensing_finalize (GObject *object)
{
  NcXcorLimberKernelWeakLensing *xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (object);

  CVodeFree (&xclkg->cvode);
  N_VDestroy (xclkg->yv);

  if (xclkg->A != NULL)
  {
    SUNMatDestroy (xclkg->A);
    xclkg->A = NULL;
  }

  if (xclkg->LS != NULL)
  {
    SUNLinSolFree (xclkg->LS);
    xclkg->LS = NULL;
  }

  SUNContext_Free (&xclkg->sunctx);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_limber_kernel_weak_lensing_parent_class)->finalize (object);
}

static gdouble _nc_xcor_limber_kernel_weak_lensing_eval (NcXcorLimberKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l);
static void _nc_xcor_limber_kernel_weak_lensing_prepare (NcXcorLimberKernel *xclk, NcHICosmo *cosmo);

/*static gdouble _nc_xcor_limber_kernel_weak_lensing_noise_spec (NcXcorLimberKernel* xclk, guint l);*/
static void _nc_xcor_limber_kernel_weak_lensing_add_noise (NcXcorLimberKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin);
guint _nc_xcor_limber_kernel_weak_lensing_obs_len (NcXcorLimberKernel *xclk);
guint _nc_xcor_limber_kernel_weak_lensing_obs_params_len (NcXcorLimberKernel *xclk);

static void
nc_xcor_limber_kernel_weak_lensing_class_init (NcXcorLimberKernelWeakLensingClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcXcorLimberKernelClass *parent_class = NC_XCOR_LIMBER_KERNEL_CLASS (klass);
  NcmModelClass *model_class            = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_xcor_limber_kernel_weak_lensing_set_property;
  model_class->get_property = &_nc_xcor_limber_kernel_weak_lensing_get_property;
  object_class->constructed = &_nc_xcor_limber_kernel_weak_lensing_constructed;
  object_class->dispose     = &_nc_xcor_limber_kernel_weak_lensing_dispose;
  object_class->finalize    = &_nc_xcor_limber_kernel_weak_lensing_finalize;

  ncm_model_class_set_name_nick (model_class, "Xcor limber weak lensing", "Xcor-WL");
  ncm_model_class_add_params (model_class, NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_SPARAM_LEN, NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_VPARAM_LEN, PROP_SIZE);

  g_object_class_install_property (object_class,
                                   PROP_DN_DZ,
                                   g_param_spec_object ("dndz",
                                                        NULL,
                                                        "Source redshift distribution",
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_CONSTRUCT | G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


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

  /**
   * NcXcorLimberKernelWeakLensingC:reltol:
   *
   * Relative tolerance used when integrating the ODE.
   * Default value: $10^{-13}$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance",
                                                        GSL_DBL_EPSILON, 1.0e-1, NC_XCOR_PRECISION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcorLimberKernelWeakLensingC:abstol:
   *
   * Absolute tolerance used when integrating the ODE.
   * Default value: $10^{-50}$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Absolute tolerance tolerance",
                                                        0.0, G_MAXDOUBLE, 1.0e-50,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  parent_class->eval      = &_nc_xcor_limber_kernel_weak_lensing_eval;
  parent_class->prepare   = &_nc_xcor_limber_kernel_weak_lensing_prepare;
  parent_class->add_noise = &_nc_xcor_limber_kernel_weak_lensing_add_noise;

  parent_class->obs_len        = &_nc_xcor_limber_kernel_weak_lensing_obs_len;
  parent_class->obs_params_len = &_nc_xcor_limber_kernel_weak_lensing_obs_params_len;

  ncm_model_class_add_impl_flag (model_class, NC_XCOR_LIMBER_KERNEL_IMPL_ALL);
}

typedef struct _src_int_params
{
  NcXcorLimberKernelWeakLensing *xclkg;
  NcHICosmo *cosmo;
  const gdouble Omega_k0;
} src_int_params;

static gdouble
_kernel_wl_W_U_eval_dndz (NcXcorLimberKernelWeakLensing *xclkg, gdouble z)
{
  const gdouble alpha = 1.0e-2;

  if (z < xclkg->dn_dz_zmin)
    return xclkg->dn_dz_min * exp (-gsl_pow_2 ((z - xclkg->dn_dz_zmin) / alpha));

  if (z > xclkg->dn_dz_zmax)
    return xclkg->dn_dz_max * exp (-gsl_pow_2 ((z - xclkg->dn_dz_zmax) / alpha));

  return ncm_spline_eval (xclkg->dn_dz, z);
}

static gint
_kernel_wl_W_U_f (sunrealtype mz, N_Vector y, N_Vector ydot, gpointer f_data)
{
  src_int_params *ts = (src_int_params *) f_data;
  NcHICosmo *cosmo   = ts->cosmo;
  NcDistance *dist   = ts->xclkg->dist;
  const gdouble z    = -mz;
  const gdouble E    = nc_hicosmo_E (cosmo, z);
  const gdouble dndz = _kernel_wl_W_U_eval_dndz (ts->xclkg, z);
  const gdouble dt   = nc_distance_transverse (dist, cosmo, z);

  NV_Ith_S (ydot, 0) = -NV_Ith_S (y, 1) / E;
  NV_Ith_S (ydot, 1) = -dndz / dt - ts->Omega_k0 * NV_Ith_S (y, 0) / E;

  return 0;
}

static gint
_kernel_wl_W_U_J (sunrealtype mz, N_Vector y, N_Vector fy, SUNMatrix J, void *jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  src_int_params *ts = (src_int_params *) jac_data;
  NcHICosmo *cosmo   = ts->cosmo;
  const gdouble z    = -mz;
  const gdouble E    = nc_hicosmo_E (cosmo, z);

  SUN_DENSE_ACCESS (J, 0, 0) = 0.0;
  SUN_DENSE_ACCESS (J, 0, 1) = -1.0 / E;

  SUN_DENSE_ACCESS (J, 1, 0) = -ts->Omega_k0 / E;
  SUN_DENSE_ACCESS (J, 1, 1) = 0.0;

  return 0;
}

static void
_nc_xcor_limber_kernel_weak_lensing_prepare (NcXcorLimberKernel *xclk, NcHICosmo *cosmo)
{
  NcXcorLimberKernelWeakLensing *xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (xclk);
  gdouble zmin, zmax, zmid;

  nc_xcor_limber_kernel_set_const_factor (xclk, 1.5 * nc_hicosmo_Omega_m0 (cosmo));
  nc_xcor_limber_kernel_get_z_range (xclk, &zmin, &zmax, &zmid);

  /* zmin should always be zero for lensing, so just checking for consistency */
  g_assert_cmpfloat (zmin, ==, 0.0);

  /* Check if the spline includes the boundaries to avoid running into errors when computing */
  g_assert (gsl_finite (ncm_spline_eval (xclkg->dn_dz, zmin)));
  g_assert (gsl_finite (ncm_spline_eval (xclkg->dn_dz, zmax)));

  /* Assign zmid as middle redshift between 0 and max of dn/dz */
  zmid = ncm_vector_get (ncm_spline_get_xv (xclkg->dn_dz), ncm_vector_get_max_index (ncm_spline_get_yv (xclkg->dn_dz))) / 2.0;
  nc_xcor_limber_kernel_set_z_range (xclk, zmin, zmax, zmid);

  nc_distance_prepare_if_needed (xclkg->dist, cosmo);

  {
    src_int_params ts    = {xclkg, cosmo, nc_hicosmo_Omega_k0 (cosmo)};
    gdouble mz_ini       = -zmax;
    const gdouble mz_end = -1.0e-10;
    GArray *x_array, *y_array;
    gint flag;

    if (xclkg->kernel_W_mz != NULL)
    {
      NcmVector *xv = ncm_spline_peek_xv (xclkg->kernel_W_mz);
      NcmVector *yv = ncm_spline_peek_yv (xclkg->kernel_W_mz);

      x_array = ncm_vector_get_array (xv);
      y_array = ncm_vector_get_array (yv);

      g_array_set_size (x_array, 0);
      g_array_set_size (y_array, 0);
    }
    else
    {
      xclkg->kernel_W_mz = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
      x_array            = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
      y_array            = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    }

    {
      /*
       *  Setting up initial conditions based on second order approximation of the
       *  original integral.
       */
      const gdouble dz = zmax * sqrt (GSL_DBL_EPSILON);
      const gdouble f  = 0.5 * gsl_pow_2 (dz) * _kernel_wl_W_U_eval_dndz (xclkg, zmax) / nc_distance_transverse (xclkg->dist, cosmo, zmax) / nc_hicosmo_E (cosmo, zmax);
      const gdouble g  = -dz *_kernel_wl_W_U_eval_dndz (xclkg, zmax) / nc_distance_transverse (xclkg->dist, cosmo, zmax);

      NV_Ith_S (xclkg->yv, 0) = f;
      NV_Ith_S (xclkg->yv, 1) = g;

      mz_ini += dz;
    }

    if (xclkg->cvode == NULL)
    {
      xclkg->cvode = CVodeCreate (CV_BDF, xclkg->sunctx);

      flag = CVodeInit (xclkg->cvode, &_kernel_wl_W_U_f, mz_ini, xclkg->yv);
      NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );

      flag = CVodeSetLinearSolver (xclkg->cvode, xclkg->LS, xclkg->A);
      NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

      flag = CVodeSetJacFn (xclkg->cvode, &_kernel_wl_W_U_J);
      NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );
    }
    else
    {
      flag = CVodeReInit (xclkg->cvode, mz_ini, xclkg->yv);
      NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );

      flag = CVodeSetLinearSolver (xclkg->cvode, xclkg->LS, xclkg->A);
      NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

      flag = CVodeSetJacFn (xclkg->cvode, &_kernel_wl_W_U_J);
      NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );
    }

    flag = CVodeSStolerances (xclkg->cvode, xclkg->reltol, xclkg->abstol);
    NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

    flag = CVodeSetMaxNumSteps (xclkg->cvode, 500000);
    NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

    flag = CVodeSetUserData (xclkg->cvode, &ts);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    flag = CVodeSetStopTime (xclkg->cvode, mz_end);
    NCM_CVODE_CHECK (&flag, "CVodeSetStopTime", 1, );

    flag = CVodeSetMaxStep (xclkg->cvode, 1.0e-1);

    g_array_append_val (x_array, mz_ini);
    g_array_append_val (y_array, NV_Ith_S (xclkg->yv, 0));

    while (TRUE)
    {
      flag = CVode (xclkg->cvode, mz_end, xclkg->yv, &mz_ini, CV_ONE_STEP);

      NCM_CVODE_CHECK (&flag, "CVode", 1, );

      g_array_append_val (x_array, mz_ini);
      g_array_append_val (y_array, NV_Ith_S (xclkg->yv, 0));

      if (mz_ini == mz_end)
        break;
    }

    {
      NcmVector *xv = ncm_vector_new_array (x_array);
      NcmVector *yv = ncm_vector_new_array (y_array);

      ncm_spline_set (xclkg->kernel_W_mz, xv, yv, TRUE);
      ncm_vector_free (xv);
      ncm_vector_free (yv);

      g_array_unref (x_array);
      g_array_unref (y_array);
    }
  }
}

static gdouble
_nc_xcor_limber_kernel_weak_lensing_eval (NcXcorLimberKernel *xclk, NcHICosmo *cosmo, gdouble z, const NcXcorKinetic *xck, gint l)
{
  NcXcorLimberKernelWeakLensing *xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (xclk);
  const gdouble nu                     = l + 0.5;
  const gdouble lfactor                = sqrt ((l + 2.0) * (l + 1.0) * l * (l - 1.0)) / (nu * nu);

  return lfactor * (1.0 + z) * xck->xi_z / xck->E_z * ncm_spline_eval (xclkg->kernel_W_mz, -z);
}

static void
_nc_xcor_limber_kernel_weak_lensing_add_noise (NcXcorLimberKernel *xclk, NcmVector *vp1, NcmVector *vp2, guint lmin)
{
  NcXcorLimberKernelWeakLensing *xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (xclk);

  /* ncm_vector_add_constant (vp1, NOISE_BIAS); / * take noise_bias into account * / */
  ncm_vector_memcpy (vp2, vp1);
  ncm_vector_add_constant (vp2, xclkg->noise);
}

guint
_nc_xcor_limber_kernel_weak_lensing_obs_len (NcXcorLimberKernel *xclk)
{
  NCM_UNUSED (xclk);

  return 2;
}

guint
_nc_xcor_limber_kernel_weak_lensing_obs_params_len (NcXcorLimberKernel *xclk)
{
  NCM_UNUSED (xclk);

  return 1;
}

/**
 * nc_xcor_limber_kernel_weak_lensing_new:
 * @zmin: a gdouble
 * @zmax: a gdouble
 * @dn_dz: a #NcmSpline
 * @nbar: a gdouble, gal density
 * @intr_shear: a gdouble, intrinsic galaxy shear
 * @dist: a #NcDistance
 *
 * Returns: a #NcXcorLimberKernelWeakLensing
 */
NcXcorLimberKernelWeakLensing *
nc_xcor_limber_kernel_weak_lensing_new (gdouble zmin, gdouble zmax, NcmSpline *dn_dz, gdouble nbar, gdouble intr_shear, NcDistance *dist)
{
  NcXcorLimberKernelWeakLensing *xclkg = g_object_new (NC_TYPE_XCOR_LIMBER_KERNEL_WEAK_LENSING,
                                                       "zmin", zmin,
                                                       "zmax", zmax,
                                                       "dndz", dn_dz,
                                                       "nbar", nbar,
                                                       "intr-shear", intr_shear,
                                                       "dist", dist,
                                                       NULL);

  return xclkg;
}

