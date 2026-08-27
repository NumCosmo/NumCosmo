/***************************************************************************
 *            nc_xcor.h
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
 * NcXCor:
 *
 * Angular auto- and cross-spectra.
 *
 * The angular power spectrum between observables $A$ and $B$ with kernels $W^A$ and $W^B$ is
 * \begin{equation}
 * C_{\ell}^{AB} = \int dz W^A(z) \int dz^\prime W^B (z^\prime) \times \int dk \frac{2}{\pi} k^2 P(k, z, z^\prime) j_{\ell}(k\chi(z)) j_{\ell} (k\chi(z^\prime)).
 * \end{equation}
 * In the Limber approximation, it reduces to
 * \begin{equation}
 * C_{\ell}^{AB} = \int_0^{z_*} dz \frac{H(z)}{c \chi^2(z)} W^A(z) W^B (z) P\left(k = \frac{\ell +1/2}{\chi(z)} , z \right),
 * \end{equation}
 * where $P\left(k = \frac{\ell +1/2}{\chi(z)} , z \right)$ is the power spectrum (a #NcmPowspec) at redshift $z$ and $chi(z)$ the comoving distance (a #NcDistance).
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/integration/ncm_integrate.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/integration/ncm_integral_nd.h"
#include "nc/xcor/nc_xcor.h"
#include "ncm/specfunc/ncm_sbessel_integrator_levin.h"
#include "nc/xcor/nc_xcor_priv.h"
#include "nc_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcXcor
{
  /*< private > */
  GObject parent_instance;
  NcDistance *dist;
  NcmPowspec *ps;
  gdouble RH;
  NcXcorMethod meth;
  gdouble reltol;
  guint ell_batch_size;
};

enum
{
  PROP_0,
  PROP_DISTANCE,
  PROP_MATTER_POWER_SPECTRUM,
  PROP_METH,
  PROP_RELTOL,
  PROP_ELL_BATCH_SIZE,
};

G_DEFINE_TYPE (NcXcor, nc_xcor, G_TYPE_OBJECT)

typedef struct _NcXcorArg
{
  NcXcor *xc;
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcmPowspec *ps;

  NcXcorKernel *xclk1;
  NcXcorKernel *xclk2;
  gint *ells;
  guint nells;
  guint comp_offset; /* index of the block component the first output maps to */

  /* Vectorized kernel integrands (for kernel cubature methods) */
  NcXcorKernelIntegrand *xclki1;
  NcXcorKernelIntegrand *xclki2;
  gdouble *W1;
  gdouble *W2;

  gdouble RH;
} NcXcorArg;

static void nc_xcor_auto_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);
static void nc_xcor_auto_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
static void nc_xcor_cross_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);
static void nc_xcor_cross_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
static void nc_xcor_kernel_auto_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);
static void nc_xcor_kernel_auto_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
static void nc_xcor_kernel_cross_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);
static void nc_xcor_kernel_cross_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);

NCM_INTEGRAL_ND_DEFINE_TYPE (NC, XCOR_AUTO, NcXcorAuto, nc_xcor_auto, nc_xcor_auto_dim, nc_xcor_auto_integ, NcXcorArg);
NCM_INTEGRAL_ND_DEFINE_TYPE (NC, XCOR_CROSS, NcXcorCross, nc_xcor_cross, nc_xcor_cross_dim, nc_xcor_cross_integ, NcXcorArg);
NCM_INTEGRAL_ND_DEFINE_TYPE (NC, XCOR_KERNEL_AUTO, NcXcorKernelAuto, nc_xcor_kernel_auto, nc_xcor_kernel_auto_dim, nc_xcor_kernel_auto_integ, NcXcorArg);
NCM_INTEGRAL_ND_DEFINE_TYPE (NC, XCOR_KERNEL_CROSS, NcXcorKernelCross, nc_xcor_kernel_cross, nc_xcor_kernel_cross_dim, nc_xcor_kernel_cross_integ, NcXcorArg);

static void
nc_xcor_init (NcXcor *xc)
{
  xc->ps             = NULL;
  xc->dist           = NULL;
  xc->RH             = 0.0;
  xc->meth           = NC_XCOR_METHOD_LIMBER_Z_GSL;
  xc->ell_batch_size = 8;
}

static void
_nc_xcor_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcor *xc = NC_XCOR (object);

  g_return_if_fail (NC_IS_XCOR (object));

  switch (prop_id)
  {
    case PROP_DISTANCE:
      nc_distance_clear (&xc->dist);
      xc->dist = g_value_dup_object (value);
      break;
    case PROP_MATTER_POWER_SPECTRUM:
      ncm_powspec_clear (&xc->ps);
      xc->ps = g_value_dup_object (value);
      break;
    case PROP_METH:
      xc->meth = g_value_get_enum (value);
      break;
    case PROP_RELTOL:
      nc_xcor_set_reltol (xc, g_value_get_double (value));
      break;
    case PROP_ELL_BATCH_SIZE:
      nc_xcor_set_ell_batch_size (xc, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcor *xc = NC_XCOR (object);

  g_return_if_fail (NC_IS_XCOR (object));

  switch (prop_id)
  {
    case PROP_DISTANCE:
      g_value_set_object (value, xc->dist);
      break;
    case PROP_MATTER_POWER_SPECTRUM:
      g_value_set_object (value, xc->ps);
      break;
    case PROP_METH:
      g_value_set_enum (value, xc->meth);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, nc_xcor_get_reltol (xc));
      break;
    case PROP_ELL_BATCH_SIZE:
      g_value_set_uint (value, nc_xcor_get_ell_batch_size (xc));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_dispose (GObject *object)
{
  NcXcor *xc = NC_XCOR (object);

  nc_distance_clear (&xc->dist);
  ncm_powspec_clear (&xc->ps);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_parent_class)->dispose (object);
}

static void
_nc_xcor_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_parent_class)->finalize (object);
}

static void
nc_xcor_class_init (NcXcorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  /*GObjectClass* parent_class = G_OBJECT_CLASS (klass); */

  object_class->set_property = &_nc_xcor_set_property;
  object_class->get_property = &_nc_xcor_get_property;
  object_class->dispose      = &_nc_xcor_dispose;
  object_class->finalize     = &_nc_xcor_finalize;

  /**
   * NcXcor:distance:
   *
   * This property keeps the distance object.
   */
  g_object_class_install_property (object_class,
                                   PROP_DISTANCE,
                                   g_param_spec_object ("distance",
                                                        NULL,
                                                        "Distance.",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcor:power-spec:
   *
   * This property keeps the matter power spectrum object.
   */
  g_object_class_install_property (object_class,
                                   PROP_MATTER_POWER_SPECTRUM,
                                   g_param_spec_object ("power-spec",
                                                        NULL,
                                                        "Matter power spectrum.",
                                                        NCM_TYPE_POWSPEC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcor:meth:
   *
   * This property keeps the method enumerator to compute the integrals.
   */
  g_object_class_install_property (object_class,
                                   PROP_METH,
                                   g_param_spec_enum ("meth",
                                                      NULL,
                                                      "Method.",
                                                      NC_TYPE_XCOR_METHOD,
                                                      NC_XCOR_METHOD_LIMBER_Z_GSL,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcor:reltol:
   *
   * This property keeps the relative tolerance.
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance.",
                                                        GSL_DBL_EPSILON, 1.0e-1, NC_XCOR_PRECISION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcXcor:ell-batch-size:
   *
   * This property keeps the multipole batch size for the kernel-space block
   * methods. Tuned for 8 or 16; #NC_XCOR_KERNEL_MAX_ELL_BLOCK is only a hard
   * safety cap. See nc_xcor_set_ell_batch_size().
   */
  g_object_class_install_property (object_class,
                                   PROP_ELL_BATCH_SIZE,
                                   g_param_spec_uint ("ell-batch-size",
                                                      NULL,
                                                      "Multipole batch size for the kernel-space block methods.",
                                                      1, NC_XCOR_KERNEL_MAX_ELL_BLOCK, 8,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static void
nc_xcor_auto_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  NcXcorAuto *xcor_auto = NC_XCOR_AUTO (intnd);
  NcXcorArg *xcor_arg   = &xcor_auto->data;

  *dim  = 1;
  *fdim = xcor_arg->nells;
}

static void
nc_xcor_cross_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  NcXcorCross *xcor_cross = NC_XCOR_CROSS (intnd);
  NcXcorArg *xcor_arg     = &xcor_cross->data;

  *dim  = 1;
  *fdim = xcor_arg->nells;
}

static void
nc_xcor_auto_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcXcorAuto *xcor_int    = NC_XCOR_AUTO (intnd);
  NcXcorArg *xcor_int_arg = &xcor_int->data;
  guint i, j;

  for (i = 0; i < npoints; i++)
  {
    const gdouble z         = ncm_vector_fast_get (x, i);
    const gdouble xi_z      = nc_distance_comoving (xcor_int_arg->dist, xcor_int_arg->cosmo, z); /* in units of Hubble radius */
    const gdouble xi_z_phys = xi_z * xcor_int_arg->RH;                                           /* in Mpc */
    const gdouble E_z       = nc_hicosmo_E (xcor_int_arg->cosmo, z);
    const NcXcorKinetic xck = { xi_z, E_z };

    for (j = 0; j < fdim; j++)
    {
      const gint l             = xcor_int_arg->ells[j];
      const gdouble k          = (l + 0.5) / (xi_z_phys); /* in Mpc-1 */
      const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xcor_int_arg->ps), NCM_MODEL (xcor_int_arg->cosmo), z, k);
      const gdouble k1z        = nc_xcor_kernel_eval_limber_z (xcor_int_arg->xclk1, xcor_int_arg->cosmo, z, &xck, l);
      const gdouble res        = gsl_pow_2 (k1z / xi_z) * power_spec / E_z;

      ncm_vector_fast_set (fval, i * fdim + j, res);
    }
  }
}

static void
nc_xcor_cross_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcXcorCross *xcor_cross = NC_XCOR_CROSS (intnd);
  NcXcorArg *xcor_arg     = &xcor_cross->data;
  guint i, j;

  for (i = 0; i < npoints; i++)
  {
    const gdouble z         = ncm_vector_fast_get (x, i);
    const gdouble xi_z      = nc_distance_comoving (xcor_arg->dist, xcor_arg->cosmo, z); /* in units of Hubble radius */
    const gdouble xi_z_phys = xi_z * xcor_arg->RH;                                       /* in Mpc */
    const gdouble E_z       = nc_hicosmo_E (xcor_arg->cosmo, z);
    const NcXcorKinetic xck = { xi_z, E_z };

    for (j = 0; j < fdim; j++)
    {
      const gint l             = xcor_arg->ells[j];
      const gdouble k          = (l + 0.5) / (xi_z_phys); /* in Mpc-1 */
      const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xcor_arg->ps), NCM_MODEL (xcor_arg->cosmo), z, k);
      const gdouble k1z        = nc_xcor_kernel_eval_limber_z (xcor_arg->xclk1, xcor_arg->cosmo, z, &xck, l);
      const gdouble k2z        = nc_xcor_kernel_eval_limber_z (xcor_arg->xclk2, xcor_arg->cosmo, z, &xck, l);
      const gdouble res        = k1z * k2z * power_spec / (xi_z * xi_z * E_z);

      ncm_vector_fast_set (fval, i * fdim + j, res);
    }
  }
}

static void
nc_xcor_kernel_auto_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  NcXcorKernelAuto *xcor_kernel_auto = NC_XCOR_KERNEL_AUTO (intnd);
  NcXcorArg *xcor_arg                = &xcor_kernel_auto->data;

  *dim  = 1;
  *fdim = xcor_arg->nells;
}

static void
nc_xcor_kernel_cross_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  NcXcorKernelCross *xcor_kernel_cross = NC_XCOR_KERNEL_CROSS (intnd);
  NcXcorArg *xcor_arg                  = &xcor_kernel_cross->data;

  *dim  = 1;
  *fdim = xcor_arg->nells;
}

static void
nc_xcor_kernel_auto_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcXcorKernelAuto *xcor_kernel_auto = NC_XCOR_KERNEL_AUTO (intnd);
  NcXcorArg *xcor_arg                = &xcor_kernel_auto->data;
  guint i;

  for (i = 0; i < npoints; i++)
  {
    const gdouble lnk = ncm_vector_fast_get (x, i);
    const gdouble k   = exp (lnk);
    const gdouble k3  = gsl_pow_3 (k);
    guint j;

    /* Evaluate all multipoles at once using pre-computed integrand. The
     * outputs are the block's components starting at comp_offset: a run of
     * multipoles sharing one k range, not necessarily the whole block. */
    nc_xcor_kernel_integrand_eval_comps (xcor_arg->xclki1, k, xcor_arg->comp_offset, fdim, xcor_arg->W1);

    for (j = 0; j < fdim; j++)
    {
      const gdouble W1  = xcor_arg->W1[xcor_arg->comp_offset + j];
      const gdouble res = k3 * W1 * W1;

      ncm_vector_fast_set (fval, i * fdim + j, res);
    }
  }
}

static void
nc_xcor_kernel_cross_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcXcorKernelCross *xcor_kernel_cross = NC_XCOR_KERNEL_CROSS (intnd);
  NcXcorArg *xcor_arg                  = &xcor_kernel_cross->data;
  guint i;

  for (i = 0; i < npoints; i++)
  {
    const gdouble lnk = ncm_vector_fast_get (x, i);
    const gdouble k   = exp (lnk);
    const gdouble k3  = gsl_pow_3 (k);
    guint j;

    /* Evaluate all multipoles at once using pre-computed integrands. The
     * outputs are the block's components starting at comp_offset: a run of
     * multipoles sharing one k range, not necessarily the whole block. */
    nc_xcor_kernel_integrand_eval_comps (xcor_arg->xclki1, k, xcor_arg->comp_offset, fdim, xcor_arg->W1);
    nc_xcor_kernel_integrand_eval_comps (xcor_arg->xclki2, k, xcor_arg->comp_offset, fdim, xcor_arg->W2);

    for (j = 0; j < fdim; j++)
    {
      const gdouble res = k3 * xcor_arg->W1[xcor_arg->comp_offset + j] * xcor_arg->W2[xcor_arg->comp_offset + j];

      ncm_vector_fast_set (fval, i * fdim + j, res);
    }
  }
}

/**
 * nc_xcor_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @meth: a #NcXcorMethod to compute the integrals
 *
 * Two methods are available to compute integrals: independent GSL numerical integration or vector integration using Sundials's CVode algorithm.
 *
 * Returns: (transfer full): a #NcXcor
 *
 */
NcXcor *
nc_xcor_new (NcDistance *dist, NcmPowspec *ps, NcXcorMethod meth)
{
  return g_object_new (NC_TYPE_XCOR,
                       "distance", dist,
                       "power-spec", ps,
                       "meth", meth,
                       NULL);
}

/**
 * nc_xcor_ref:
 * @xc: a #NcXcor
 *
 * Increments the reference count of @xc.
 *
 * Returns: (transfer full): @xc
 */
NcXcor *
nc_xcor_ref (NcXcor *xc)
{
  return g_object_ref (xc);
}

/**
 * nc_xcor_free:
 * @xc: a #NcXcor
 *
 * Decrements the reference count of @xc, and frees it if the count reaches 0.
 *
 */
void
nc_xcor_free (NcXcor *xc)
{
  g_object_unref (xc);
}

/**
 * nc_xcor_clear:
 * @xc: a #NcXcor
 *
 * If *@xc is not %NULL, decrements the reference count of @xc, and frees it if the
 * count reaches 0.
 *
 */
void
nc_xcor_clear (NcXcor **xc)
{
  g_clear_object (xc);
}

/**
 * nc_xcor_get_meth:
 * @xc: a #NcXcor
 *
 * Returns: the #NcXcorMethod used by @xc
 */
NcXcorMethod
nc_xcor_get_meth (NcXcor *xc)
{
  return xc->meth;
}

/**
 * nc_xcor_set_reltol:
 * @xc: a #NcXcor
 * @reltol: a relative tolerance
 *
 * Sets the relative tolerance of @xc.
 *
 */
void
nc_xcor_set_reltol (NcXcor *xc, const gdouble reltol)
{
  xc->reltol = reltol;
}

/**
 * nc_xcor_get_reltol:
 * @xc: a #NcXcor
 *
 * Returns: the relative tolerance of @xc
 */
gdouble
nc_xcor_get_reltol (NcXcor *xc)
{
  return xc->reltol;
}

/**
 * nc_xcor_set_ell_batch_size:
 * @xc: a #NcXcor
 * @ell_batch_size: multipole batch size
 *
 * Sets the multipole batch size used by the kernel-space block methods
 * (%NC_XCOR_METHOD_KERNEL_CUBATURE and %NC_XCOR_METHOD_KERNEL_EXACT): each
 * batch builds one k-space closure per kernel and shares it across the whole
 * batch.
 *
 * The Levin machinery is tuned for 8 (the default) or 16; wider batches are
 * counterproductive, not faster. #NC_XCOR_KERNEL_MAX_ELL_BLOCK is a hard
 * safety cap on a single nc_xcor_kernel_get_eval_vectorized() call, not a
 * setting to aim for.
 *
 */
void
nc_xcor_set_ell_batch_size (NcXcor *xc, const guint ell_batch_size)
{
  g_return_if_fail ((ell_batch_size > 0) && (ell_batch_size <= (guint) NC_XCOR_KERNEL_MAX_ELL_BLOCK));

  xc->ell_batch_size = ell_batch_size;
}

/**
 * nc_xcor_get_ell_batch_size:
 * @xc: a #NcXcor
 *
 * Returns: the multipole batch size
 */
guint
nc_xcor_get_ell_batch_size (NcXcor *xc)
{
  return xc->ell_batch_size;
}

/**
 * nc_xcor_prepare:
 * @xc: a #NcXcor
 * @cosmo: a #NcHICosmo
 *
 * Prepares @xc for computation.
 *
 */
void
nc_xcor_prepare (NcXcor *xc, NcHICosmo *cosmo)
{
  nc_distance_prepare_if_needed (xc->dist, cosmo);
  ncm_powspec_prepare_if_needed (xc->ps, NCM_MODEL (cosmo));

  xc->RH = nc_hicosmo_RH_Mpc (cosmo);
}

typedef struct _xcor_gsl
{
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcmPowspec *ps;

  NcXcorKernel *xclk1;
  NcXcorKernel *xclk2;
  guint l;

  gdouble RH;
} xcor_gsl;

/*
 * QUADPACK's status is a statement about certification, not about the answer:
 * GSL_EROUND means the error estimate stopped improving before reaching the
 * requested tolerance. Where the request carries a safety margin over
 * NcXcor:reltol, the margin is a goal and reltol is the requirement, so the
 * achieved error decides -- aborting on the status alone throws away a result
 * that already meets what the caller asked for.
 */
static void
_nc_xcor_check_qag_status (const gchar *where, gint ret, gdouble reltol, gdouble result, gdouble err)
{
  if (ret == GSL_SUCCESS)
    return;

  if ((ret == GSL_EROUND) && (err <= reltol * fabs (result)))
    return;

  g_error ("%s: %s. Result % 22.15g with estimated error %.8e, against the "
           "requested %.8e (NcXcor:reltol %.8e). The quadrature could not "
           "certify that accuracy on this integrand; loosen NcXcor:reltol to "
           "what the integrand carries.",
           where, gsl_strerror (ret), result, err, reltol * fabs (result), reltol);
}

static gdouble
_xcor_limber_z_gsl_cross_int (gdouble z, gpointer ptr)
{
  xcor_gsl *xclki          = (xcor_gsl *) ptr;
  const gdouble xi_z       = nc_distance_comoving (xclki->dist, xclki->cosmo, z); /* in units of Hubble radius */
  const gdouble xi_z_phys  = xi_z * xclki->RH;                                    /* in Mpc */
  const gdouble E_z        = nc_hicosmo_E (xclki->cosmo, z);
  const NcXcorKinetic xck  = { xi_z, E_z };
  const gdouble k          = (xclki->l + 0.5) / (xi_z_phys); /* in Mpc-1 */
  const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xclki->ps), NCM_MODEL (xclki->cosmo), z, k);

  const gdouble k1z = nc_xcor_kernel_eval_limber_z (xclki->xclk1, xclki->cosmo, z, &xck, xclki->l);
  const gdouble k2z = nc_xcor_kernel_eval_limber_z (xclki->xclk2, xclki->cosmo, z, &xck, xclki->l);

  return k1z * k2z * power_spec / (xi_z * xi_z * E_z);
}

static gdouble
_xcor_limber_z_gsl_auto_int (gdouble z, gpointer ptr)
{
  xcor_gsl *xclki          = (xcor_gsl *) ptr;
  const gdouble xi_z       = nc_distance_comoving (xclki->dist, xclki->cosmo, z); /* in units of Hubble radius */
  const gdouble xi_z_phys  = xi_z * xclki->RH;                                    /* in Mpc */
  const gdouble E_z        = nc_hicosmo_E (xclki->cosmo, z);
  const NcXcorKinetic xck  = { xi_z, E_z };
  const gdouble k          = (xclki->l + 0.5) / (xi_z_phys); /* in Mpc-1 */
  const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xclki->ps), NCM_MODEL (xclki->cosmo), z, k);
  const gdouble k1z        = nc_xcor_kernel_eval_limber_z (xclki->xclk1, xclki->cosmo, z, &xck, xclki->l);

  return gsl_pow_2 (k1z / xi_z) * power_spec / E_z;
}

static void
_nc_xcor_limber_z_gsl (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gdouble zmin, gdouble zmax, gboolean isauto, NcmVector *vp)
{
  xcor_gsl xclki;
  gdouble r, err;
  gsl_function F;
  guint i;
  gint ret;

  xclki.xclk1 = xclk1;
  xclki.xclk2 = xclk2;
  xclki.cosmo = cosmo;
  xclki.dist  = xc->dist;
  xclki.ps    = xc->ps;
  xclki.RH    = xc->RH;

  zmin = (zmin != 0.0) ? zmin : 1.0e-6;

  if (isauto)
    F.function = &_xcor_limber_z_gsl_auto_int;
  else
    F.function = &_xcor_limber_z_gsl_cross_int;

  F.params = &xclki;

  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  for (i = 0; i < lmax - lmin + 1; i++)
  {
    xclki.l = lmin + i;
    /* GSL integration sometimes underestimates the error, so we multiply the relative tolerance by 1e-2 */
    ret = gsl_integration_qag (&F, zmin, zmax, 0.0, xc->reltol * 1.0e-2, NCM_INTEGRAL_PARTITION, 6, *w, &r, &err);

    _nc_xcor_check_qag_status ("_nc_xcor_limber_z_gsl", ret, xc->reltol, r, err);

    ncm_vector_set (vp, i, r);
  }

  ncm_memory_pool_return (w);
}

static void
_nc_xcor_cubature_worker (NcmIntegralND *xcor_int_nd, NcXcorArg *xcor_arg, gdouble zmin, gdouble zmax, guint lmin, guint lmax, NcmVector *vp)
{
  NcmVector *z_min   = ncm_vector_new_data_static (&zmin, 1, 1);
  NcmVector *z_max   = ncm_vector_new_data_static (&zmax, 1, 1);
  const guint size   = lmax - lmin + 1;
  NcmVector *err     = ncm_vector_new (size);
  GArray *ells_array = g_array_new (FALSE, FALSE, sizeof (gint));
  const gint block   = sqrt (size);
  guint i;

  zmin = (zmin != 0.0) ? zmin : 1.0e-6;

  ncm_integral_nd_set_reltol (xcor_int_nd, xcor_arg->xc->reltol);
  ncm_integral_nd_set_abstol (xcor_int_nd, 0.0);
  ncm_integral_nd_set_method (xcor_int_nd, NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V);

  g_array_set_size (ells_array, size);

  for (i = 0; i < size; i++)
  {
    g_array_index (ells_array, gint, i) = lmin + i;
  }

  for (i = 0; i + block < size; i += block)
  {
    NcmVector *vp_i  = ncm_vector_get_subvector (vp, i, block);
    NcmVector *err_i = ncm_vector_get_subvector (err, i, block);

    xcor_arg->ells  = &g_array_index (ells_array, gint, i);
    xcor_arg->nells = block;

    ncm_integral_nd_eval (xcor_int_nd, z_min, z_max, vp_i, err_i);
    ncm_vector_free (vp_i);
    ncm_vector_free (err_i);
  }

  {
    NcmVector *vp_i  = ncm_vector_get_subvector (vp, i, size - i);
    NcmVector *err_i = ncm_vector_get_subvector (err, i, size - i);

    xcor_arg->ells  = &g_array_index (ells_array, gint, i);
    xcor_arg->nells = size - i;

    ncm_integral_nd_eval (xcor_int_nd, z_min, z_max, vp_i, err_i);
    ncm_vector_free (vp_i);
    ncm_vector_free (err_i);
  }

  g_array_unref (ells_array);
  ncm_vector_free (z_min);
  ncm_vector_free (z_max);
  ncm_vector_free (err);
}

static void
_nc_xcor_limber_z_cubature (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gdouble zmin, gdouble zmax, gboolean isauto, NcmVector *vp)
{
  if (isauto)
  {
    NcXcorAuto *xcor_int       = g_object_new (nc_xcor_auto_get_type (), NULL);
    NcXcorArg *xcor_arg        = &xcor_int->data;
    NcmIntegralND *xcor_int_nd = NCM_INTEGRAL_ND (xcor_int);

    xcor_arg->xc    = xc;
    xcor_arg->dist  = xc->dist;
    xcor_arg->ps    = xc->ps;
    xcor_arg->RH    = xc->RH;
    xcor_arg->cosmo = cosmo;
    xcor_arg->xclk1 = xclk1;
    xcor_arg->xclk2 = xclk1;

    _nc_xcor_cubature_worker (xcor_int_nd, xcor_arg, zmin, zmax, lmin, lmax, vp);

    g_object_unref (xcor_int);
  }
  else
  {
    NcXcorCross *xcor_cross    = g_object_new (nc_xcor_cross_get_type (), NULL);
    NcXcorArg *xcor_arg        = &xcor_cross->data;
    NcmIntegralND *xcor_int_nd = NCM_INTEGRAL_ND (xcor_cross);

    xcor_arg->xc    = xc;
    xcor_arg->dist  = xc->dist;
    xcor_arg->ps    = xc->ps;
    xcor_arg->RH    = xc->RH;
    xcor_arg->cosmo = cosmo;
    xcor_arg->xclk1 = xclk1;
    xcor_arg->xclk2 = xclk2;

    _nc_xcor_cubature_worker (xcor_int_nd, xcor_arg, zmin, zmax, lmin, lmax, vp);

    g_object_unref (xcor_cross);
  }
}

static gdouble
_xcor_kernel_gsl_cross_int (gdouble lnk, gpointer ptr)
{
  NcXcorKernelIntegrand **xclki = (NcXcorKernelIntegrand **) ptr;
  const gdouble k               = exp (lnk);
  gdouble W1[1], W2[1];

  nc_xcor_kernel_integrand_eval (xclki[0], k, W1);
  nc_xcor_kernel_integrand_eval (xclki[1], k, W2);

  return gsl_pow_3 (k) * W1[0] * W2[0];
}

static gdouble
_xcor_kernel_gsl_auto_int (gdouble lnk, gpointer ptr)
{
  NcXcorKernelIntegrand **xclki = (NcXcorKernelIntegrand **) ptr;
  const gdouble k               = exp (lnk);
  gdouble W[1];

  nc_xcor_kernel_integrand_eval (xclki[0], k, W);

  return gsl_pow_3 (k) * W[0] * W[0];
}

/*
 * Five-node Gauss-Legendre rule on [-1, 1]. Exact through degree 9, which is
 * one more than the degree-8 polynomial the outer integrand k^2 W_i W_j is on
 * each knot panel of a cubic spline: 2 from k^2 plus 3 from each spline.
 */
#define NC_XCOR_GL5_N 5

static const gdouble _nc_xcor_gl5_x[NC_XCOR_GL5_N] = {
  -0.9061798459386640, -0.5384693101056831, 0.0,
  0.5384693101056831, 0.9061798459386640
};

static const gdouble _nc_xcor_gl5_w[NC_XCOR_GL5_N] = {
  0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
  0.4786286704993665, 0.2369268850561891
};

/*
 * The two GL(5) sweeps over a pair's merged panels. Auto and cross differ only
 * in whether the second integrand is evaluated at all, which is fixed for the
 * whole sweep -- so they are separate functions and the inner loops carry no
 * branch. Each is a flat loop over panel x node x multipole.
 */
/*
 * Everything the error estimate of _nc_xcor_kernel_integrate_block_exact ()
 * needs, accumulated in the same pass as the integral itself. The peaks enter
 * only as multipliers of the accumulated integrals, so a running maximum is
 * enough and no second sweep is required.
 */
typedef struct _NcXcorGL5Err
{
  gdouble *prod;  /* int k^2 |W1 W2| */
  gdouble *abs1;  /* int k^2 |W1|    */
  gdouble *abs2;  /* int k^2 |W2|    */
  gdouble *peak1; /* max |W1|        */
  gdouble *peak2; /* max |W2|        */
} NcXcorGL5Err;

static void
_nc_xcor_gl5_sweep_auto (NcXcorKernelIntegrand *xclki, GArray *edges, guint nell, gdouble *W, gdouble *sum, NcXcorGL5Err *err)
{
  guint ie, ig, il;

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble panel_lo = g_array_index (edges, gdouble, ie);
    const gdouble panel_hi = g_array_index (edges, gdouble, ie + 1);
    const gdouble mid      = 0.5 * (panel_lo + panel_hi);
    const gdouble half     = 0.5 * (panel_hi - panel_lo);

    for (ig = 0; ig < NC_XCOR_GL5_N; ig++)
    {
      const gdouble k = mid + half * _nc_xcor_gl5_x[ig];
      const gdouble w = half * _nc_xcor_gl5_w[ig] * k * k;

      nc_xcor_kernel_integrand_eval (xclki, k, W);

      for (il = 0; il < nell; il++)
      {
        const gdouble term = w * W[il] * W[il];

        sum[il] += term;

        if (err != NULL)
        {
          err->prod[il] += fabs (term);
          err->abs1[il] += fabs (w * W[il]);
          err->peak1[il] = GSL_MAX (err->peak1[il], fabs (W[il]));
        }
      }
    }
  }
}

static void
_nc_xcor_gl5_sweep_cross (NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, GArray *edges, guint nell, gdouble *W1, gdouble *W2, gdouble *sum, NcXcorGL5Err *err)
{
  guint ie, ig, il;

  for (ie = 0; ie + 1 < edges->len; ie++)
  {
    const gdouble panel_lo = g_array_index (edges, gdouble, ie);
    const gdouble panel_hi = g_array_index (edges, gdouble, ie + 1);
    const gdouble mid      = 0.5 * (panel_lo + panel_hi);
    const gdouble half     = 0.5 * (panel_hi - panel_lo);

    for (ig = 0; ig < NC_XCOR_GL5_N; ig++)
    {
      const gdouble k = mid + half * _nc_xcor_gl5_x[ig];
      const gdouble w = half * _nc_xcor_gl5_w[ig] * k * k;

      nc_xcor_kernel_integrand_eval (xclki1, k, W1);
      nc_xcor_kernel_integrand_eval (xclki2, k, W2);

      for (il = 0; il < nell; il++)
      {
        const gdouble term = w * W1[il] * W2[il];

        sum[il] += term;

        if (err != NULL)
        {
          err->prod[il] += fabs (term);
          err->abs1[il] += fabs (w * W1[il]);
          err->abs2[il] += fabs (w * W2[il]);
          err->peak1[il] = GSL_MAX (err->peak1[il], fabs (W1[il]));
          err->peak2[il] = GSL_MAX (err->peak2[il], fabs (W2[il]));
        }
      }
    }
  }
}

/*
 * Common refinement of two knot sets, clipped to [@k_min, @k_max]. Both are
 * sorted, so this is a linear merge; duplicates are dropped so no zero-width
 * panel survives. The result is pre-sized to the exact upper bound of a merge,
 * so the append loop never reallocates: the union of two sorted sets cannot
 * exceed their combined length.
 */
static GArray *
_nc_xcor_merge_knots (NcmVector *knots1, NcmVector *knots2, gdouble k_min, gdouble k_max)
{
  const guint len1 = ncm_vector_len (knots1);
  const guint len2 = ncm_vector_len (knots2);
  GArray *edges    = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), len1 + len2);
  guint i1         = 0;
  guint i2         = 0;

  while ((i1 < len1) || (i2 < len2))
  {
    const gdouble x1 = (i1 < len1) ? ncm_vector_get (knots1, i1) : GSL_POSINF;
    const gdouble x2 = (i2 < len2) ? ncm_vector_get (knots2, i2) : GSL_POSINF;
    const gdouble x  = GSL_MIN (x1, x2);

    if (x1 <= x2)
      i1++;

    if (x2 <= x1)
      i2++;

    if ((x < k_min) || (x > k_max))
      continue;

    if ((edges->len > 0) && (x <= g_array_index (edges, gdouble, edges->len - 1)))
      continue;

    g_array_append_val (edges, x);
  }

  return edges;
}

/*
 * Exact outer quadrature for one pair, on the union of the two kernels' own
 * knot sets.
 *
 * Each kernel is sampled independently -- the same per-kernel closures
 * KERNEL_CUBATURE builds and NcXcorSolver caches -- so the two splines live on
 * different abscissas. On the common refinement of those abscissas each spline
 * is still a single cubic piece per panel, so k^2 W_i W_j is a degree-8
 * polynomial there and GL(5) integrates it exactly. Merging two knot sets is
 * all the coupling the exactness argument needs; a shared abscissa built by
 * sampling the kernels jointly is not required, and costs about twice as much
 * to produce.
 *
 * The range is the intersection of the two integrands' fitted domains, exactly
 * as _nc_xcor_kernel_integrate_block_cubature() uses, because NcmSpline does
 * not range-check and an out-of-domain evaluation returns extrapolation rather
 * than a small number.
 *
 * ## What this means for error control
 *
 * Exact is meant literally: refining every panel fourfold moves the result by
 * 1e-15 to 1e-12, which is rounding. So there is nothing for an embedded
 * quadrature rule to measure here -- a Kronrod extension, or GL(5) against
 * GL(9), would report machine zero on every call and never fire. Do not add
 * one; it would be false confidence rather than error control.
 *
 * The error is entirely in the closure: a spline is a fit, and $C_\ell$ can be
 * far smaller than the integral of the absolute integrand, which amplifies
 * that fit's error. Two disjoint Gaussian bins already reach a cancellation of
 * 1.4e4 by $\ell = 9$, so a closure good to 1e-8 leaves 1e-4 on $C_\ell$.
 * @vp_err reports that product; see nc_xcor_compute_full().
 *
 * The same C^2 kinks are why the adaptive kernel-space methods struggle on
 * this integrand and this one does not. A cubic spline's third derivative
 * jumps at every knot, so an adaptive scheme subdivides at each of them
 * forever, chasing a relative criterion that the fit's own error puts out of
 * reach. Here the knots *are* the panel edges, so the kinks fall between
 * panels and never enter a rule's smoothness assumption.
 */
void
_nc_xcor_kernel_integrate_block_exact (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcmVector *vp, NcmVector *vp_err)
{
  const guint nell           = lmax - lmin + 1;
  const gdouble const_factor = 2.0 / (M_PI * gsl_pow_3 (xc->RH));
  NcmVector *knots1          = nc_xcor_kernel_integrand_peek_knots (xclki1);
  NcmVector *knots2          = isauto ? knots1 : nc_xcor_kernel_integrand_peek_knots (xclki2);
  const gdouble reltol1      = nc_xcor_kernel_integrand_get_reltol (xclki1);
  const gdouble reltol2      = nc_xcor_kernel_integrand_get_reltol (xclki2);
  const gdouble sabs1        = nc_xcor_kernel_integrand_get_scaled_abstol (xclki1);
  const gdouble sabs2        = nc_xcor_kernel_integrand_get_scaled_abstol (xclki2);
  gdouble k_min1, k_max1, k_min2, k_max2, k_min, k_max;
  NcXcorGL5Err err_acc, *err = NULL;
  gdouble *sum, *W1, *W2;
  GArray *edges;
  guint il;

  if (ncm_vector_len (vp) != nell)
    g_error ("_nc_xcor_kernel_integrate_block_exact: vector size does not match multipole limits");

  if ((vp_err != NULL) && (ncm_vector_len (vp_err) != nell))
    g_error ("_nc_xcor_kernel_integrate_block_exact: error vector size does not match multipole limits");

  if ((knots1 == NULL) || (knots2 == NULL))
    g_error ("_nc_xcor_kernel_integrate_block_exact: %s method requires spline-backed "
             "integrands, which report their knots.", "NC_XCOR_METHOD_KERNEL_EXACT");

  nc_xcor_kernel_integrand_get_range (xclki1, &k_min1, &k_max1);
  nc_xcor_kernel_integrand_get_range (xclki2, &k_min2, &k_max2);

  k_min = GSL_MAX (k_min1, k_min2);
  k_max = GSL_MIN (k_max1, k_max2);

  ncm_vector_set_zero (vp);

  if (vp_err != NULL)
    ncm_vector_set_zero (vp_err);

  if (k_min >= k_max)
    return;

  edges = _nc_xcor_merge_knots (knots1, knots2, k_min, k_max);

  if (edges->len < 2)
  {
    g_array_unref (edges);

    return;
  }

  sum = g_new0 (gdouble, nell);
  W1  = g_new0 (gdouble, nc_xcor_kernel_integrand_get_len (xclki1));
  W2  = isauto ? W1 : g_new0 (gdouble, nc_xcor_kernel_integrand_get_len (xclki2));

  if (vp_err != NULL)
  {
    err_acc.prod  = g_new0 (gdouble, nell);
    err_acc.abs1  = g_new0 (gdouble, nell);
    err_acc.abs2  = isauto ? err_acc.abs1 : g_new0 (gdouble, nell);
    err_acc.peak1 = g_new0 (gdouble, nell);
    err_acc.peak2 = isauto ? err_acc.peak1 : g_new0 (gdouble, nell);
    err           = &err_acc;
  }

  /* The auto/cross distinction is fixed for the whole sweep, so it is resolved
   * once here rather than tested at every quadrature node. */
  if (isauto)
    _nc_xcor_gl5_sweep_auto (xclki1, edges, nell, W1, sum, err);
  else
    _nc_xcor_gl5_sweep_cross (xclki1, xclki2, edges, nell, W1, W2, sum, err);

  for (il = 0; il < nell; il++)
    ncm_vector_set (vp, il, const_factor * sum[il]);

  /* The quadrature is exact, so the only error is the closures' own, propagated
   * through d(W1 W2) = |W1| dW2 + |W2| dW1 with dWi the fit criterion of
   * nc_xcor_kernel_integrand_set_tolerances (). Its two halves separate: the
   * relative one rides on the product, the peak-scaled floor on each closure
   * against the other's peak -- which is how a floor set per closure ends up
   * squared where both of them sit on it. */
  if (vp_err != NULL)
  {
    for (il = 0; il < nell; il++)
    {
      const gdouble rel_term   = (reltol1 + reltol2) * err->prod[il];
      const gdouble floor_term = sabs1 * err->peak1[il] * err->abs2[il] +
                                 sabs2 * err->peak2[il] * err->abs1[il];

      ncm_vector_set (vp_err, il, const_factor * (rel_term + floor_term));
    }

    g_free (err->prod);
    g_free (err->abs1);
    g_free (err->peak1);

    if (!isauto)
    {
      g_free (err->abs2);
      g_free (err->peak2);
    }
  }

  g_free (sum);
  g_free (W1);

  if (!isauto)
    g_free (W2);

  g_array_unref (edges);
}

gboolean
_nc_xcor_kernels_limber_disjoint (NcXcorKernel *xclk1, NcXcorKernel *xclk2, gboolean isauto, guint *l_zero)
{
  gdouble zmin, zmax, zmid, zmin_2, zmax_2, zmid_2;
  gint l_limber_1, l_limber_2;

  if (isauto)
    return FALSE;  /* A kernel always overlaps itself. */

  l_limber_1 = nc_xcor_kernel_get_l_limber (xclk1);
  l_limber_2 = nc_xcor_kernel_get_l_limber (xclk2);

  /* A negative l_limber pins the kernel to the non-Limber tier at every
   * multipole, so no threshold exists. */
  if ((l_limber_1 < 0) || (l_limber_2 < 0))
    return FALSE;

  nc_xcor_kernel_get_z_range (xclk1, &zmin, &zmax, &zmid);
  nc_xcor_kernel_get_z_range (xclk2, &zmin_2, &zmax_2, &zmid_2);

  if (GSL_MAX (zmin, zmin_2) < GSL_MIN (zmax, zmax_2))
    return FALSE;

  /* Both kernels are in the Limber tier from the larger of the two thresholds
   * up; below it at least one is still non-Limber and the pair correlates. */
  *l_zero = (guint) GSL_MAX (l_limber_1, l_limber_2);

  return TRUE;
}

static gboolean
_nc_xcor_meth_is_kernel_space (NcXcorMethod meth)
{
  switch (meth)
  {
    case NC_XCOR_METHOD_KERNEL_GSL:
    case NC_XCOR_METHOD_KERNEL_CUBATURE:
    case NC_XCOR_METHOD_KERNEL_EXACT:
      return TRUE;

    default:
      return FALSE;
  }
}

static void
_nc_xcor_kernel_exact (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gboolean isauto, NcmVector *vp, NcmVector *vp_err)
{
  NcmSBesselIntegrator *sbi1 = nc_xcor_kernel_peek_integrator (xclk1);
  NcmSBesselIntegrator *sbi2 = nc_xcor_kernel_peek_integrator (xclk2);
  const guint size           = lmax - lmin + 1;
  const guint block          = xc->ell_batch_size;
  guint i;

  _nc_xcor_check_kernel_tolerance (xc, xclk1);

  if (!isauto)
    _nc_xcor_check_kernel_tolerance (xc, xclk2);

  /* Either kernel's integrator serves a kernel that carries none of its own. */
  if (sbi1 == NULL)
    sbi1 = sbi2;

  if (sbi2 == NULL)
    sbi2 = sbi1;

  /* Batched by NcXcor:ell-batch-size exactly like _nc_xcor_kernel_cubature():
   * one k-space closure per kernel per batch. The batching is not merely an
   * optimization here -- a single closure spanning more than
   * NC_XCOR_KERNEL_MAX_ELL_BLOCK multipoles is a hard error in
   * nc_xcor_kernel_get_eval_vectorized_full(), so an unbatched sweep aborted on
   * any range wider than that. */
  for (i = 0; i < size; i += block)
  {
    const guint block_lmin = lmin + i;
    const guint block_lmax = MIN (block_lmin + block - 1, lmax);
    NcmVector *vp_i        = ncm_vector_get_subvector (vp, i, block_lmax - block_lmin + 1);
    NcmVector *vp_err_i    = (vp_err != NULL) ? ncm_vector_get_subvector (vp_err, i, block_lmax - block_lmin + 1) : NULL;
    NcXcorKernelIntegrand *xclki1;
    NcXcorKernelIntegrand *xclki2;

    xclki1 = nc_xcor_kernel_get_eval_vectorized_full (xclk1, cosmo, block_lmin, block_lmax, sbi1);
    xclki2 = isauto ? NULL : nc_xcor_kernel_get_eval_vectorized_full (xclk2, cosmo, block_lmin, block_lmax, sbi2);

    _nc_xcor_kernel_integrate_block_exact (xc, xclki1, isauto ? xclki1 : xclki2, block_lmin, block_lmax, isauto, vp_i, vp_err_i);

    nc_xcor_kernel_integrand_unref (xclki1);

    if (xclki2 != NULL)
      nc_xcor_kernel_integrand_unref (xclki2);

    ncm_vector_free (vp_i);
    ncm_vector_clear (&vp_err_i);
  }
}

/*
 * Subinterval breakpoints in ln k for one pair, from the integrands' own
 * knots, with the integration limits as the first and last entry. %NULL when
 * either integrand is not spline-backed and so reports no knots.
 */
static GArray *
_nc_xcor_kernel_gsl_breakpoints (NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, gdouble k_min, gdouble k_max)
{
  NcmVector *knots1 = nc_xcor_kernel_integrand_peek_knots (xclki1);
  NcmVector *knots2 = (xclki2 == NULL) ? knots1 : nc_xcor_kernel_integrand_peek_knots (xclki2);
  GArray *edges;
  GArray *pts;
  guint i;

  if ((knots1 == NULL) || (knots2 == NULL))
    return NULL;

  edges = _nc_xcor_merge_knots (knots1, knots2, k_min, k_max);
  pts   = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), edges->len + 2);

  {
    const gdouble ln_k_min = log (k_min);

    g_array_append_val (pts, ln_k_min);
  }

  for (i = 0; i < edges->len; i++)
  {
    const gdouble ln_k = log (g_array_index (edges, gdouble, i));

    if (ln_k > g_array_index (pts, gdouble, pts->len - 1))
      g_array_append_val (pts, ln_k);
  }

  {
    const gdouble ln_k_max = log (k_max);

    if (ln_k_max > g_array_index (pts, gdouble, pts->len - 1))
      g_array_append_val (pts, ln_k_max);
    else
      g_array_index (pts, gdouble, pts->len - 1) = ln_k_max;
  }

  g_array_unref (edges);

  /* Degenerate merge (a domain narrower than one panel): nothing to break on. */
  if ((pts->len < 2) || (g_array_index (pts, gdouble, pts->len - 1) <= g_array_index (pts, gdouble, pts->len - 2)))
  {
    g_array_unref (pts);

    return NULL;
  }

  return pts;
}

void
_nc_xcor_kernel_gsl (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gboolean isauto, NcmVector *vp)
{
  const guint nell              = ncm_vector_len (vp);
  const gdouble const_factor    = 2.0 / (M_PI * gsl_pow_3 (xc->RH));
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  NcXcorKernelIntegrand *xclki_array[2];
  gsl_function F;
  guint i;
  gint ret;

  if (nell != lmax - lmin + 1)
    g_error ("_nc_xcor_kernel_gsl: vector size does not match multipole limits");

  if (lmax < lmin)
    g_error ("_nc_xcor_kernel_gsl: lmax < lmin");

  if (isauto)
    F.function = &_xcor_kernel_gsl_auto_int;
  else
    F.function = &_xcor_kernel_gsl_cross_int;

  F.params = xclki_array;

  for (i = 0; i < nell; i++)
  {
    const guint ell = lmin + i;
    gdouble k_min, k_max, result, err;
    GArray *breakpoints;

    /* Build the integrand(s) first, then read the outer bound off their own
     * fitted domain (get_range) -- NOT the independent Limber-band formula
     * from nc_xcor_kernel_get_k_range(), which has no guarantee of matching
     * it (see plan doc dev-notes/xcor_ultralevin_batching_plan.md §3). */
    xclki_array[0] = nc_xcor_kernel_get_eval (xclk1, cosmo, ell);
    nc_xcor_kernel_integrand_get_range (xclki_array[0], &k_min, &k_max);

    if (isauto)
    {
      xclki_array[1] = NULL;
    }
    else
    {
      gdouble k2_min, k2_max;

      xclki_array[1] = nc_xcor_kernel_get_eval (xclk2, cosmo, ell);
      nc_xcor_kernel_integrand_get_range (xclki_array[1], &k2_min, &k2_max);

      k_min = GSL_MAX (k_min, k2_min);
      k_max = GSL_MIN (k_max, k2_max);
    }

    /* Integrated on the integrands' own knots when they have them, and over
     * the bare interval when they do not. A spline-backed integrand is
     * piecewise cubic: its third derivative jumps at every knot, which a
     * Gauss-Kronrod rule spanning several knots cannot see. Its error estimate
     * then saturates -- measured on the CMB ISW kernel, a panel holding four
     * or five knots reports about 80 times the error the same region reports
     * split on them -- and bisection stops improving it, which is exactly what
     * QUADPACK reports as roundoff (GSL_EROUND). Breaking on the knots leaves
     * every subinterval a single cubic piece: the same integral converges to
     * machine precision instead of stalling near 1e-6, for about twice the
     * evaluations. No safety margin is applied over NcXcor:reltol; the margin
     * that used to be (reltol * 1e-2) only moved the stall.
     */
    breakpoints = _nc_xcor_kernel_gsl_breakpoints (xclki_array[0], isauto ? NULL : xclki_array[1], k_min, k_max);

    if (breakpoints != NULL)
      ret = gsl_integration_qagp (&F, (gdouble *) breakpoints->data, breakpoints->len, 0.0, xc->reltol, NCM_INTEGRAL_PARTITION, *w, &result, &err);
    else
      ret = gsl_integration_qag (&F, log (k_min), log (k_max), 0.0, xc->reltol, NCM_INTEGRAL_PARTITION, 6, *w, &result, &err);

    _nc_xcor_check_qag_status ("_nc_xcor_kernel_gsl", ret, xc->reltol, result, err);

    if (breakpoints != NULL)
      g_array_unref (breakpoints);

    ncm_vector_set (vp, i, const_factor * result);

    nc_xcor_kernel_integrand_unref (xclki_array[0]);

    if (!isauto)
      nc_xcor_kernel_integrand_unref (xclki_array[1]);
  }

  ncm_memory_pool_return (w);
}

/*
 * The outer k-integral cannot resolve the closure more finely than the closure
 * itself is built. pcubature answers an impossible tolerance by exhausting its
 * Clenshaw-Curtis levels and reporting failure, far from the cause, so catch
 * the mismatch here where both numbers are in view.
 */
void
_nc_xcor_check_kernel_tolerance (NcXcor *xc, NcXcorKernel *xclk)
{
  NcmSBesselIntegrator *sbi = nc_xcor_kernel_peek_integrator (xclk);
  gdouble closure_reltol;

  if ((sbi == NULL) || !NCM_IS_SBESSEL_INTEGRATOR_LEVIN (sbi))
    return;

  {
    NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (sbi);

    closure_reltol = GSL_MAX (ncm_sbessel_integrator_levin_get_reltol (sbilv),
                              ncm_sbessel_integrator_levin_get_cheb_reltol (sbilv));
  }

  if (xc->reltol < closure_reltol)
    g_error ("_nc_xcor_check_kernel_tolerance: NcXcor:reltol is %.17g but kernel %s builds "
             "its k-space closure only to %.17g (the looser of the integrator's reltol and "
             "cheb-reltol). The outer integral cannot converge to more precision than the "
             "integrand carries. Loosen NcXcor:reltol to at least %.17g, or construct the "
             "integrator with tighter tolerances.",
             xc->reltol, G_OBJECT_TYPE_NAME (xclk), closure_reltol, closure_reltol);
}

/*
 * The k range one multipole of a pair is supported on: both integrands' own
 * ranges for that component, intersected.
 */
static void
_nc_xcor_kernel_cubature_comp_range (NcXcorArg *xcor_arg, gboolean isauto, guint i, gdouble *k_min, gdouble *k_max)
{
  nc_xcor_kernel_integrand_get_range_comp (xcor_arg->xclki1, i, k_min, k_max);

  if (!isauto)
  {
    gdouble k2_min, k2_max;

    nc_xcor_kernel_integrand_get_range_comp (xcor_arg->xclki2, i, &k2_min, &k2_max);

    *k_min = GSL_MAX (*k_min, k2_min);
    *k_max = GSL_MIN (*k_max, k2_max);
  }
}

/*
 * One ell block, integrated in runs of consecutive multipoles that share a k
 * range.
 *
 * A block's multipoles share one fitted domain but not necessarily one
 * support: under the Limber approximation each is confined to its own band in
 * k and vanishes outside it, so integrating the whole block over the shared
 * domain hands the quadrature a step per multipole. An adaptive rule does not
 * merely work harder there, it can miss the feature outright -- a region whose
 * nodes all land where the component vanishes reports no error and is never
 * subdivided -- and a kernel whose window does not tend to zero at its own
 * band edge (CMB ISW, whose radial integral truncates at recombination) then
 * loses a finite part of the integral under a converged status. Splitting on
 * the range puts each step on an integration limit, where it is not a
 * discontinuity of the integrand at all.
 *
 * Where every multipole shares one range -- the non-Limber case, whose support
 * is the fitted domain itself -- this is one call over the whole block, as
 * before.
 */
static void
_nc_xcor_kernel_cubature_runs (NcmIntegralND *xcor_int_nd, NcXcorArg *xcor_arg, gboolean isauto, guint nells, NcmVector *vp_i, NcmVector *err_i)
{
  NcmVector *lnk_min = ncm_vector_new (1);
  NcmVector *lnk_max = ncm_vector_new (1);
  guint j            = 0;

  while (j < nells)
  {
    gdouble k_min, k_max;
    guint run = 1;

    _nc_xcor_kernel_cubature_comp_range (xcor_arg, isauto, j, &k_min, &k_max);

    while (j + run < nells)
    {
      gdouble k_min_run, k_max_run;

      _nc_xcor_kernel_cubature_comp_range (xcor_arg, isauto, j + run, &k_min_run, &k_max_run);

      if ((k_min_run != k_min) || (k_max_run != k_max))
        break;

      run++;
    }

    if (k_min < k_max)
    {
      NcmVector *vp_j  = ncm_vector_get_subvector (vp_i, j, run);
      NcmVector *err_j = ncm_vector_get_subvector (err_i, j, run);

      xcor_arg->comp_offset = j;
      xcor_arg->nells       = run;

      ncm_vector_set (lnk_min, 0, log (k_min));
      ncm_vector_set (lnk_max, 0, log (k_max));

      ncm_integral_nd_eval (xcor_int_nd, lnk_min, lnk_max, vp_j, err_j);

      ncm_vector_free (vp_j);
      ncm_vector_free (err_j);
    }
    else
    {
      /* Empty support: the pair's bands do not meet inside the domain. */
      guint m;

      for (m = 0; m < run; m++)
      {
        ncm_vector_set (vp_i, j + m, 0.0);
        ncm_vector_set (err_i, j + m, 0.0);
      }
    }

    j += run;
  }

  ncm_vector_free (lnk_min);
  ncm_vector_free (lnk_max);
}

static void
_nc_xcor_kernel_cubature (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gboolean isauto, NcmVector *vp)
{
  const guint size  = lmax - lmin + 1;
  const guint block = xc->ell_batch_size;
  guint i;

  _nc_xcor_check_kernel_tolerance (xc, xclk1);

  if (!isauto)
    _nc_xcor_check_kernel_tolerance (xc, xclk2);

  /* Batched by NcXcor:ell-batch-size exactly like _nc_xcor_kernel_exact():
   * one k-space closure per kernel per batch, built here and handed to the
   * same per-block integrator #NcXcorSolver drives with closures of its own. */
  for (i = 0; i < size; i += block)
  {
    const guint nells      = MIN (block, size - i);
    const guint block_lmin = lmin + i;
    const guint block_lmax = block_lmin + nells - 1;
    NcmVector *vp_i        = ncm_vector_get_subvector (vp, i, nells);
    NcXcorKernelIntegrand *xclki1;
    NcXcorKernelIntegrand *xclki2;

    xclki1 = nc_xcor_kernel_get_eval_vectorized (xclk1, cosmo, block_lmin, block_lmax);
    xclki2 = isauto ? NULL : nc_xcor_kernel_get_eval_vectorized (xclk2, cosmo, block_lmin, block_lmax);

    _nc_xcor_kernel_integrate_block_cubature (xc, xclki1, xclki2, block_lmin, block_lmax, isauto, vp_i);

    nc_xcor_kernel_integrand_unref (xclki1);

    if (xclki2 != NULL)
      nc_xcor_kernel_integrand_unref (xclki2);

    ncm_vector_free (vp_i);
  }
}

/*
 * _nc_xcor_kernel_integrate_block_cubature:
 * @xc: a #NcXcor
 * @xclki1: a pre-built #NcXcorKernelIntegrand for the first kernel,
 * covering multipoles [@lmin, @lmax]
 * @xclki2: (nullable): a pre-built #NcXcorKernelIntegrand for the second
 * kernel, covering the same range, or %NULL for an auto-correlation
 * @lmin: minimum multipole, matching @xclki1's (and @xclki2's) own range
 * @lmax: maximum multipole, matching @xclki1's (and @xclki2's) own range
 * @isauto: %TRUE for an auto-correlation (only @xclki1 is used)
 * @vp: (out): output vector of length (@lmax - @lmin + 1)
 *
 * Computes the outer k-integral for one ell-block from integrand(s) the
 * caller already built, instead of building them internally the way
 * _nc_xcor_kernel_cubature() does. This is #NcXcorSolver's building block
 * for sharing one kernel's k-space closure across every request that needs
 * it in a given block, instead of rebuilding it once per pair (see plan
 * doc dev-notes/xcor_ultralevin_batching_plan.md §5-§6). Not public: see
 * nc_xcor_priv.h.
 *
 * The caller retains ownership of @xclki1/@xclki2 -- they are not unreffed
 * here.
 */
void
_nc_xcor_kernel_integrate_block_cubature (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcmVector *vp)
{
  const guint size           = lmax - lmin + 1;
  const gdouble const_factor = 2.0 / (M_PI * gsl_pow_3 (xc->RH));
  NcmVector *err             = ncm_vector_new (size);
  NcmIntegralND *xcor_int_nd;
  NcXcorArg *xcor_arg;

  if (ncm_vector_len (vp) != size)
    g_error ("_nc_xcor_kernel_integrate_block_cubature: vector size does not match multipole limits");

  if (isauto)
  {
    NcXcorKernelAuto *xcor_kernel_auto = g_object_new (nc_xcor_kernel_auto_get_type (), NULL);

    xcor_int_nd = NCM_INTEGRAL_ND (xcor_kernel_auto);
    xcor_arg    = &xcor_kernel_auto->data;
  }
  else
  {
    NcXcorKernelCross *xcor_kernel_cross = g_object_new (nc_xcor_kernel_cross_get_type (), NULL);

    xcor_int_nd = NCM_INTEGRAL_ND (xcor_kernel_cross);
    xcor_arg    = &xcor_kernel_cross->data;
  }

  xcor_arg->xc     = xc;
  xcor_arg->RH     = xc->RH;
  xcor_arg->xclki1 = xclki1;
  xcor_arg->W1     = g_new (gdouble, size);

  if (!isauto)
  {
    xcor_arg->xclki2 = xclki2;
    xcor_arg->W2     = g_new (gdouble, size);
  }

  ncm_integral_nd_set_reltol (xcor_int_nd, xc->reltol);
  ncm_integral_nd_set_abstol (xcor_int_nd, 0.0);
  ncm_integral_nd_set_method (xcor_int_nd, NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V);

  _nc_xcor_kernel_cubature_runs (xcor_int_nd, xcor_arg, isauto, size, vp, err);

  ncm_vector_scale (vp, const_factor);

  g_free (xcor_arg->W1);

  if (!isauto)
    g_free (xcor_arg->W2);

  g_object_unref (xcor_int_nd);
  ncm_vector_free (err);
}

/*
 * Dispatches one contiguous multipole range to the configured kernel-space
 * method. Split out so nc_xcor_compute() can drive it for a sub-range without
 * duplicating the switch, which is what a range straddling a Limber threshold
 * needs.
 */
static void
_nc_xcor_kernel_space_compute (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gboolean isauto, NcmVector *vp, NcmVector *vp_err)
{
  /* Only the exact method knows its own error budget. The others leave NaN
   * rather than a zero that would read as "no error". */
  if (vp_err != NULL)
    ncm_vector_set_all (vp_err, GSL_NAN);

  switch (xc->meth)
  {
    case NC_XCOR_METHOD_KERNEL_GSL:
      _nc_xcor_kernel_gsl (xc, xclk1, xclk2, cosmo, lmin, lmax, isauto, vp);
      break;

    case NC_XCOR_METHOD_KERNEL_CUBATURE:
      _nc_xcor_kernel_cubature (xc, xclk1, xclk2, cosmo, lmin, lmax, isauto, vp);
      break;

    case NC_XCOR_METHOD_KERNEL_EXACT:
      _nc_xcor_kernel_exact (xc, xclk1, xclk2, cosmo, lmin, lmax, isauto, vp, vp_err);
      break;

    default:
      g_assert_not_reached ();
  }
}

/**
 * nc_xcor_compute_full:
 * @xc: a #NcXcor
 * @xclk1: a #NcXcorKernel
 * @xclk2: (nullable): a #NcXcorKernel, or %NULL for the auto spectrum
 * @cosmo: a #NcHICosmo
 * @lmin: a #guint
 * @lmax: a #guint
 * @vp: a #NcmVector
 * @vp_err: (nullable): a #NcmVector for the error estimate, or %NULL
 *
 * As nc_xcor_compute(), additionally filling @vp_err with an estimate of the
 * absolute error on each $C_\ell$ -- the same length as @vp, and in the same
 * units, so @vp_err[i] / @vp[i] is the relative estimate.
 *
 * Only %NC_XCOR_METHOD_KERNEL_EXACT provides one; every other method leaves
 * NaN, which is deliberate -- a zero there would read as "no error".
 *
 * The estimate is not a quadrature error, and does not belong to the k
 * integral at all. That method integrates its spline closures exactly, so the
 * outer quadrature contributes nothing; what @vp_err reports is a
 * **kernel-building** error -- how well nc_xcor_kernel_get_eval_vectorized_full()
 * fitted $W_\ell(k)$ -- propagated through the conditioning of this particular
 * pair. The quadrature's only contribution is the amplification factor, which
 * is also the one thing that cannot be known before the pair is formed.
 *
 * Propagating $\delta (W^A W^B) = \vert W^A \vert \delta W^B + \vert W^B
 * \vert \delta W^A$ with the fit criterion of
 * nc_xcor_kernel_integrand_set_tolerances(), $\delta W^i \le \epsilon_i \vert
 * W^i \vert + a_i W^i_\mathrm{max}$, splits into two terms:
 *
 * $$ \sigma_\ell \simeq (\epsilon_A + \epsilon_B) \int \mathrm{d}k\, k^2 \vert W^A_\ell W^B_\ell \vert
 *    + a_A W^A_{\ell,\mathrm{max}} \int \mathrm{d}k\, k^2 \vert W^B_\ell \vert
 *    + a_B W^B_{\ell,\mathrm{max}} \int \mathrm{d}k\, k^2 \vert W^A_\ell \vert $$
 *
 * The first rides on the product, so it is the one the pair's cancellation
 * amplifies. The second and third are each closure's peak-scaled floor weighted
 * by the *other* closure's true size, and they are usually the larger of the
 * two: for cluster top-hat bins at the library defaults they dominate the
 * relative term by one to two orders. That is worth stating plainly, because a
 * floor set per closure is often assumed to reach $C_\ell$ only squared. It
 * does so only where both closures sit on their floors at once; wherever just
 * one does, it is linear in $a$ and weighted by the other's real amplitude.
 *
 * ## How conservative, measured
 *
 * Everything above bounds the *criterion*, which the refinement stops at rather
 * than beats, so @vp_err is an upper bound and not an estimate. Against a
 * reference built at reltol $10^{-10}$, cluster top-hat bins over
 * $\ell = 2 \dots 9$ at the library defaults:
 *
 * | pair | true relative error | @vp_err | ratio |
 * |---|---|---|---|
 * | auto | 4e-6 to 1.3e-4 | 1.5e-3 to 3.7e-3 | 12-860 |
 * | cross, adjacent bins | 7e-4 to 0.13 | 0.2 to 4.5 | 35-320 |
 * | cross, separated bins | 0.07 to 17 | 6 to 96 | 4-160 |
 *
 * So read it as a ceiling: a small @vp_err is a strong statement, a large one
 * warrants checking rather than despair. It is loosest exactly where the answer
 * is healthiest. Tightening it needs the refinement's *achieved* residual
 * rather than the tolerance it was asked for, which
 * NcmFunctionSampleSet does not currently report.
 *
 * Those figures are a **worst case**, and deliberately so: a cluster top-hat is
 * discontinuous at its bin edges, so its $W_\ell(k)$ decays only as $1/k$ and is
 * the hardest closure in the library to fit. On the same comoving shells a
 * smooth kernel needs 161 knots against the top-hat's 541, and its cross
 * spectrum is accurate to 7.7e-4 rather than 0.13 -- a factor of 165.
 *
 * That comparison also bounds what @vp_err can do. It reported 4.5 and 2.4 for
 * those two pairs, whose true errors differ by that factor of 165: it is
 * tracking the pair's cancellation, which is similar for both, and not the fit
 * quality, which is not. Read it as a conditioning flag rather than an accuracy
 * figure. Discriminating between them needs the *achieved* residual, as above.
 *
 * One thing pushes the other way, against the conservatism: the criterion is an
 * $L^2$ norm over the whole multipole block at each $k$, not over one
 * multipole, so a multipole that is sub-dominant within its block is held only
 * to the block's norm. For those, $\epsilon \vert W_\ell \vert$ understates the
 * fit error.
 *
 * **What it does not cover**, and it is the same classification again -- the
 * other kernel-building error, which is a range rather than a residual. The
 * outer integral runs over the intersection of
 * the two closures' fitted k-ranges, and @vp_err measures only what is inside
 * that range -- it cannot see what the intersection discarded. Two kernels
 * whose closures are fitted on different k-supports lose the non-overlapping
 * part silently, and that loss grows with separation, the same direction in
 * which the cancellation grows. A small @vp_err therefore means the quadrature
 * and the fit are in hand, not that the $C_\ell$ is right.
 *
 */
void
nc_xcor_compute_full (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, NcmVector *vp, NcmVector *vp_err)
{
  const guint nell           = ncm_vector_len (vp);
  const gboolean isauto      = (xclk2 == NULL) || (xclk2 == xclk1);
  const gdouble const_factor = 1.0 / gsl_pow_3 (xc->RH);
  guint i;
  gdouble zmin, zmax, zmid;

  if (nell != lmax - lmin + 1)
    g_error ("nc_xcor_compute_full: vector size does not match multipole limits");

  if ((vp_err != NULL) && (ncm_vector_len (vp_err) != nell))
    g_error ("nc_xcor_compute_full: error vector size does not match multipole limits");

  if (isauto)
    xclk2 = xclk1;

  /* The kernel-space (non-Limber) methods perform each kernel's radial
   * integral separately and couple the two only through the outer k integral,
   * so two kernels with disjoint redshift support still have a non-zero cross
   * spectrum -- two disjoint radial shells are correlated through the same 3D
   * field. The z-overlap short-circuit below is therefore specific to the
   * Limber-z tier, whose C_l is a single integral over the common support, and
   * must not be applied here.
   *
   * The exception is a kernel-space method running kernels that are themselves
   * in the Limber tier: there W(k) is supported only on its own per-ell band,
   * disjoint bins have disjoint support, and the Cl is zero. Integrating it
   * anyway multiplies the two exponential extrapolation tails, which is a
   * numerical smoothing device rather than physics, and yields a large
   * spurious cross spectrum.
   *
   * The tier is chosen per multipole, so a range straddling the threshold is
   * split: the tail from l_zero up is zeroed, the head below it is integrated
   * normally. */
  if (_nc_xcor_meth_is_kernel_space (xc->meth))
  {
    guint l_zero = 0;

    if (_nc_xcor_kernels_limber_disjoint (xclk1, xclk2, isauto, &l_zero) && (l_zero <= lmax))
    {
      ncm_vector_set_zero (vp);

      /* The zeroed tail is exactly zero, not merely small, so its error is
       * zero too -- the head below overwrites its own entries. */
      if (vp_err != NULL)
        ncm_vector_set_zero (vp_err);

      if (l_zero > lmin)
      {
        NcmVector *vp_head     = ncm_vector_get_subvector (vp, 0, l_zero - lmin);
        NcmVector *vp_err_head = (vp_err != NULL) ? ncm_vector_get_subvector (vp_err, 0, l_zero - lmin) : NULL;

        _nc_xcor_kernel_space_compute (xc, xclk1, xclk2, cosmo, lmin, l_zero - 1, isauto, vp_head, vp_err_head);
        ncm_vector_free (vp_head);
        ncm_vector_clear (&vp_err_head);
      }
    }
    else
    {
      _nc_xcor_kernel_space_compute (xc, xclk1, xclk2, cosmo, lmin, lmax, isauto, vp, vp_err);
    }

    return;
  }

  if (vp_err != NULL)
    ncm_vector_set_all (vp_err, GSL_NAN);

  nc_xcor_kernel_get_z_range (xclk1, &zmin, &zmax, &zmid);

  if (!isauto)
  {
    gdouble zmin_2, zmax_2, zmid_2;

    nc_xcor_kernel_get_z_range (xclk2, &zmin_2, &zmax_2, &zmid_2);
    zmin = GSL_MAX (zmin, zmin_2);
    zmax = GSL_MIN (zmax, zmax_2);
  }

  if (zmin < zmax)
  {
    switch (xc->meth)
    {
      case NC_XCOR_METHOD_LIMBER_Z_GSL:
        _nc_xcor_limber_z_gsl (xc, xclk1, xclk2, cosmo, lmin, lmax, zmin, zmax, isauto, vp);
        break;
      case NC_XCOR_METHOD_LIMBER_Z_CUBATURE:
        _nc_xcor_limber_z_cubature (xc, xclk1, xclk2, cosmo, lmin, lmax, zmin, zmax, isauto, vp);
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
        break;                   /* LCOV_EXCL_LINE */
    }

    if (isauto)
    {
      for (i = 0; i < nell; i++)
      {
        const guint ell                = lmin + i;
        const gdouble const_factor_ell = nc_xcor_kernel_eval_limber_z_prefactor (xclk1, cosmo, ell);

        ncm_vector_mulby (vp, i, const_factor * const_factor_ell * const_factor_ell);
      }
    }
    else
    {
      for (i = 0; i < nell; i++)
      {
        const guint ell                 = lmin + i;
        const gdouble const_factor_ell1 = nc_xcor_kernel_eval_limber_z_prefactor (xclk1, cosmo, ell);
        const gdouble const_factor_ell2 = nc_xcor_kernel_eval_limber_z_prefactor (xclk2, cosmo, ell);

        ncm_vector_mulby (vp, i, const_factor * const_factor_ell1 * const_factor_ell2);
      }
    }
  }
  else
  {
    ncm_vector_set_zero (vp);
  }
}

/**
 * nc_xcor_compute:
 * @xc: a #NcXcor
 * @xclk1: a #NcXcorKernel
 * @xclk2: (nullable): a #NcXcorKernel, or %NULL for the auto spectrum
 * @cosmo: a #NcHICosmo
 * @lmin: a #guint
 * @lmax: a #guint
 * @vp: a #NcmVector
 *
 * Performs the computation of the power spectrum $C_{\ell}^{AB}$. The kernels of
 * observables A and B are @xclk1 and @xclk2. If @xclk2 is NULL, the auto power
 * spectrum is computed. The result for multipoles lmin to lmax (included) is stored in
 * the #NcmVector @vp.
 *
 * Use nc_xcor_compute_full() to also obtain a per-multipole error estimate.
 *
 */
void
nc_xcor_compute (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, NcmVector *vp)
{
  nc_xcor_compute_full (xc, xclk1, xclk2, cosmo, lmin, lmax, vp, NULL);
}

