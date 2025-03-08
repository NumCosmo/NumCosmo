/***************************************************************************
 *            nc_xcor.h
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

#include "math/ncm_integrate.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_cfg.h"
#include "math/ncm_serialize.h"
#include "math/ncm_integral_nd.h"
#include "xcor/nc_xcor.h"
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
  NcXcorLimberMethod meth;
};

enum
{
  PROP_0,
  PROP_DISTANCE,
  PROP_MATTER_POWER_SPECTRUM,
  PROP_METH,
};

G_DEFINE_TYPE (NcXcor, nc_xcor, G_TYPE_OBJECT)

typedef struct _NcXcorLimberArg
{
  NcXcor *xc;
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcmPowspec *ps;

  NcXcorLimberKernel *xclk1;
  NcXcorLimberKernel *xclk2;
  gint *ells;
  guint nells;

  gdouble RH;
} NcXcorLimberArg;

static void nc_xcor_limber_auto_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);
static void nc_xcor_limber_auto_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
static void nc_xcor_limber_cross_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);
static void nc_xcor_limber_cross_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);

NCM_INTEGRAL_ND_DEFINE_TYPE (NC, XCOR_LIMBER_AUTO, NcXcorLimberAuto, nc_xcor_limber_auto, nc_xcor_limber_auto_dim, nc_xcor_limber_auto_integ, NcXcorLimberArg);
NCM_INTEGRAL_ND_DEFINE_TYPE (NC, XCOR_LIMBER_CROSS, NcXcorLimberCross, nc_xcor_limber_cross, nc_xcor_limber_cross_dim, nc_xcor_limber_cross_integ, NcXcorLimberArg);

static void
nc_xcor_init (NcXcor *xc)
{
  xc->ps   = NULL;
  xc->dist = NULL;
  xc->RH   = 0.0;
  xc->meth = NC_XCOR_LIMBER_METHOD_GSL;
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
                                                      NC_TYPE_XCOR_LIMBER_METHOD,
                                                      NC_XCOR_LIMBER_METHOD_GSL,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

static void
nc_xcor_limber_auto_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  NcXcorLimberAuto *xcor_auto = NC_XCOR_LIMBER_AUTO (intnd);
  NcXcorLimberArg *xcor_arg   = &xcor_auto->data;

  *dim  = 1;
  *fdim = xcor_arg->nells;
}

static void
nc_xcor_limber_cross_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  NcXcorLimberCross *xcor_cross = NC_XCOR_LIMBER_CROSS (intnd);
  NcXcorLimberArg *xcor_arg     = &xcor_cross->data;

  *dim  = 1;
  *fdim = xcor_arg->nells;
}

static void
nc_xcor_limber_auto_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcXcorLimberAuto *xcor_int    = NC_XCOR_LIMBER_AUTO (intnd);
  NcXcorLimberArg *xcor_int_arg = &xcor_int->data;
  guint i, j;

  for (i = 0; i < npoints; i++)
  {
    const gdouble z         = ncm_vector_fast_get (x, i);
    const gdouble xi_z      = nc_distance_comoving (xcor_int_arg->dist, xcor_int_arg->cosmo, z); /* in units of Hubble radius */
    const gdouble xi_z_phys = xi_z * xcor_int_arg->RH;                                           /* in Mpc */
    const gdouble E_z       = nc_hicosmo_E (xcor_int_arg->cosmo, z);
    const NcXcorKinetic xck = { xi_z, E_z };

    if (z == 0.0)
    {
      ncm_vector_set_zero (fval);
    }
    else
    {
      for (j = 0; j < fdim; j++)
      {
        const gint l             = xcor_int_arg->ells[j];
        const gdouble k          = (l + 0.5) / (xi_z_phys); /* in Mpc-1 */
        const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xcor_int_arg->ps), NCM_MODEL (xcor_int_arg->cosmo), z, k);
        const gdouble k1z        = nc_xcor_limber_kernel_eval (xcor_int_arg->xclk1, xcor_int_arg->cosmo, z, &xck, l);
        const gdouble res        = E_z * gsl_pow_2 (k1z / xi_z) * power_spec;

        ncm_vector_fast_set (fval, i * fdim + j, res);
      }
    }
  }
}

static void
nc_xcor_limber_cross_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcXcorLimberCross *xcor_cross = NC_XCOR_LIMBER_CROSS (intnd);
  NcXcorLimberArg *xcor_arg     = &xcor_cross->data;
  guint i, j;

  for (i = 0; i < npoints; i++)
  {
    const gdouble z         = ncm_vector_fast_get (x, i);
    const gdouble xi_z      = nc_distance_comoving (xcor_arg->dist, xcor_arg->cosmo, z); /* in units of Hubble radius */
    const gdouble xi_z_phys = xi_z * xcor_arg->RH;                                       /* in Mpc */
    const gdouble E_z       = nc_hicosmo_E (xcor_arg->cosmo, z);
    const NcXcorKinetic xck = { xi_z, E_z };

    if (z == 0.0)
    {
      ncm_vector_set_zero (fval);
    }
    else
    {
      for (j = 0; j < fdim; j++)
      {
        const gint l             = xcor_arg->ells[j];
        const gdouble k          = (l + 0.5) / (xi_z_phys); /* in Mpc-1 */
        const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xcor_arg->ps), NCM_MODEL (xcor_arg->cosmo), z, k);
        const gdouble k1z        = nc_xcor_limber_kernel_eval (xcor_arg->xclk1, xcor_arg->cosmo, z, &xck, l);
        const gdouble k2z        = nc_xcor_limber_kernel_eval (xcor_arg->xclk2, xcor_arg->cosmo, z, &xck, l);
        const gdouble res        = E_z * k1z * k2z * power_spec / (xi_z * xi_z);

        ncm_vector_fast_set (fval, i * fdim + j, res);
      }
    }
  }
}

/**
 * nc_xcor_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @meth: a #NcXcorLimberMethod to compute the Limber integrals
 *
 * Two methods are available to compute Limber-approximated integrals: independent GSL numerical integration or vector integration using Sundials's CVode algorithm.
 *
 * Returns: (transfer full): a #NcXcor
 *
 */
NcXcor *
nc_xcor_new (NcDistance *dist, NcmPowspec *ps, NcXcorLimberMethod meth)
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

typedef struct _xcor_limber_gsl
{
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcmPowspec *ps;

  NcXcorLimberKernel *xclk1;
  NcXcorLimberKernel *xclk2;
  guint l;

  gdouble RH;
} xcor_limber_gsl;

static gdouble
_xcor_limber_gsl_cross_int (gdouble z, gpointer ptr)
{
  xcor_limber_gsl *xclki   = (xcor_limber_gsl *) ptr;
  const gdouble xi_z       = nc_distance_comoving (xclki->dist, xclki->cosmo, z); /* in units of Hubble radius */
  const gdouble xi_z_phys  = xi_z * xclki->RH;                                    /* in Mpc */
  const gdouble E_z        = nc_hicosmo_E (xclki->cosmo, z);
  const NcXcorKinetic xck  = { xi_z, E_z };
  const gdouble k          = (xclki->l + 0.5) / (xi_z_phys); /* in Mpc-1 */
  const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xclki->ps), NCM_MODEL (xclki->cosmo), z, k);

  const gdouble k1z = nc_xcor_limber_kernel_eval (xclki->xclk1, xclki->cosmo, z, &xck, xclki->l);
  const gdouble k2z = nc_xcor_limber_kernel_eval (xclki->xclk2, xclki->cosmo, z, &xck, xclki->l);

  return E_z * k1z * k2z * power_spec / (xi_z * xi_z);
}

static gdouble
_xcor_limber_gsl_auto_int (gdouble z, gpointer ptr)
{
  xcor_limber_gsl *xclki   = (xcor_limber_gsl *) ptr;
  const gdouble xi_z       = nc_distance_comoving (xclki->dist, xclki->cosmo, z); /* in units of Hubble radius */
  const gdouble xi_z_phys  = xi_z * xclki->RH;                                    /* in Mpc */
  const gdouble E_z        = nc_hicosmo_E (xclki->cosmo, z);
  const NcXcorKinetic xck  = { xi_z, E_z };
  const gdouble k          = (xclki->l + 0.5) / (xi_z_phys); /* in Mpc-1 */
  const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xclki->ps), NCM_MODEL (xclki->cosmo), z, k);
  const gdouble k1z        = nc_xcor_limber_kernel_eval (xclki->xclk1, xclki->cosmo, z, &xck, xclki->l);

  return E_z * gsl_pow_2 (k1z / xi_z) * power_spec;
}

static void
_nc_xcor_limber_gsl (NcXcor *xc, NcXcorLimberKernel *xclk1, NcXcorLimberKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gdouble zmin, gdouble zmax, gboolean isauto, NcmVector *vp)
{
  xcor_limber_gsl xclki;
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

  if (isauto)
    F.function = &_xcor_limber_gsl_auto_int;
  else
    F.function = &_xcor_limber_gsl_cross_int;

  F.params = &xclki;

  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  for (i = 0; i < lmax - lmin + 1; i++)
  {
    xclki.l = lmin + i;
    ret     = gsl_integration_qag (&F, zmin, zmax, 0.0, NC_XCOR_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &r, &err);

    if (ret != GSL_SUCCESS)
      g_error ("_nc_xcor_limber_gsl: %s.", gsl_strerror (ret));

    ncm_vector_set (vp, i, r);
  }

  ncm_memory_pool_return (w);
}

static void
_nc_xcor_limber_cubature_worker (NcmIntegralND *xcor_int_nd, NcXcorLimberArg *xcor_arg, gdouble zmin, gdouble zmax, guint lmin, guint lmax, NcmVector *vp)
{
  NcmVector *z_min   = ncm_vector_new_data_static (&zmin, 1, 1);
  NcmVector *z_max   = ncm_vector_new_data_static (&zmax, 1, 1);
  const guint size   = lmax - lmin + 1;
  NcmVector *err     = ncm_vector_new (size);
  GArray *ells_array = g_array_new (FALSE, FALSE, sizeof (gint));
  const gint block   = 30;
  guint i;

  zmin = zmin ? zmin != 0.0 : 1.0e-6;

  ncm_integral_nd_set_reltol (xcor_int_nd, 1.0e-5);
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
_nc_xcor_limber_cubature (NcXcor *xc, NcXcorLimberKernel *xclk1, NcXcorLimberKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gdouble zmin, gdouble zmax, gboolean isauto, NcmVector *vp)
{
  if (isauto)
  {
    NcXcorLimberAuto *xcor_int = g_object_new (nc_xcor_limber_auto_get_type (), NULL);
    NcXcorLimberArg *xcor_arg  = &xcor_int->data;
    NcmIntegralND *xcor_int_nd = NCM_INTEGRAL_ND (xcor_int);

    xcor_arg->dist  = xc->dist;
    xcor_arg->ps    = xc->ps;
    xcor_arg->RH    = xc->RH;
    xcor_arg->cosmo = cosmo;
    xcor_arg->xclk1 = xclk1;
    xcor_arg->xclk2 = xclk1;

    _nc_xcor_limber_cubature_worker (xcor_int_nd, xcor_arg, zmin, zmax, lmin, lmax, vp);

    g_object_unref (xcor_int);
  }
  else
  {
    NcXcorLimberCross *xcor_cross = g_object_new (nc_xcor_limber_cross_get_type (), NULL);
    NcXcorLimberArg *xcor_arg     = &xcor_cross->data;
    NcmIntegralND *xcor_int_nd    = NCM_INTEGRAL_ND (xcor_cross);

    xcor_arg->dist  = xc->dist;
    xcor_arg->ps    = xc->ps;
    xcor_arg->RH    = xc->RH;
    xcor_arg->cosmo = cosmo;
    xcor_arg->xclk1 = xclk1;
    xcor_arg->xclk2 = xclk2;

    _nc_xcor_limber_cubature_worker (xcor_int_nd, xcor_arg, zmin, zmax, lmin, lmax, vp);

    g_object_unref (xcor_cross);
  }
}

/**
 * nc_xcor_limber:
 * @xc: a #NcXcor
 * @xclk1: a #NcXcorLimberKernel
 * @xclk2: a #NcXcorLimberKernel
 * @cosmo: a #NcHICosmo
 * @lmin: a #guint
 * @lmax: a #guint
 * @vp: a #NcmVector
 *
 * Performs the computation of the power spectrum $C_{\ell}^{AB}$ in the Limber
 * approximation. The kernels of observables A and B are @xclk1 and @xclk2. If @xclk2 is
 * NULL, the auto power spectrum is computed. The result for multipoles lmin to lmax
 * (included) is stored in the #NcmVector @vp.
 *
 */
void
nc_xcor_limber (NcXcor *xc, NcXcorLimberKernel *xclk1, NcXcorLimberKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, NcmVector *vp)
{
  const guint nell          = ncm_vector_len (vp);
  const gboolean isauto     = (xclk2 == NULL) || (xclk2 == xclk1);
  const gdouble cons_factor = ((isauto) ?
                               gsl_pow_2 (nc_xcor_limber_kernel_get_const_factor (xclk1)) :
                               nc_xcor_limber_kernel_get_const_factor (xclk1) *
                               nc_xcor_limber_kernel_get_const_factor (xclk2)) / gsl_pow_3 (xc->RH);
  gdouble zmin, zmax, zmid;

  if (nell != lmax - lmin + 1)
    g_error ("nc_xcor_limber: vector size does not match multipole limits");

  nc_xcor_limber_kernel_get_z_range (xclk1, &zmin, &zmax, &zmid);

  if (!isauto)
  {
    gdouble zmin_2, zmax_2, zmid_2;

    nc_xcor_limber_kernel_get_z_range (xclk2, &zmin_2, &zmax_2, &zmid_2);
    zmin = GSL_MAX (zmin, zmin_2);
    zmax = GSL_MIN (zmax, zmax_2);
  }
  else
  {
    xclk2 = xclk1;
  }

  if (zmin < zmax)
  {
    switch (xc->meth)
    {
      case NC_XCOR_LIMBER_METHOD_GSL:
        _nc_xcor_limber_gsl (xc, xclk1, xclk2, cosmo, lmin, lmax, zmin, zmax, isauto, vp);
        break;
      case NC_XCOR_LIMBER_METHOD_CUBATURE:
        _nc_xcor_limber_cubature (xc, xclk1, xclk2, cosmo, lmin, lmax, zmin, zmax, isauto, vp);
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
        break;                   /* LCOV_EXCL_LINE */
    }

    ncm_vector_scale (vp, cons_factor);
  }
  else
  {
    ncm_vector_set_zero (vp);
  }
}

