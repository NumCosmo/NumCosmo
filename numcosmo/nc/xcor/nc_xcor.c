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
   * This property keeps the multipole batch size for cubature methods.
   */
  g_object_class_install_property (object_class,
                                   PROP_ELL_BATCH_SIZE,
                                   g_param_spec_uint ("ell-batch-size",
                                                      NULL,
                                                      "Multipole batch size for cubature methods.",
                                                      1, 64, 8,
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

static gint
_nc_xcor_gdouble_cmp (gconstpointer a, gconstpointer b)
{
  const gdouble da = *(const gdouble *) a;
  const gdouble db = *(const gdouble *) b;

  return (da > db) - (da < db);
}

/*
 * Evaluates the outer k-integral, splitting it at the points where the
 * integrand is discontinuous.
 *
 * Limber collapses the line-of-sight integral to the single point
 * xi = (l + 1/2) / k, so a kernel whose support ends at xi_max without
 * vanishing there turns that edge into a step in k, at a different k for every
 * multipole in the batch. The p-adaptive cubature raises the polynomial order
 * on a fixed interval, which cannot fit a step: it exhausts its Clenshaw-Curtis
 * levels and reports failure. The step locations are known analytically, so
 * integrate between them and sum -- every piece is then smooth for every
 * multipole, including those still identically zero over it.
 *
 * A smooth integrand reports no breakpoints and takes the single-interval path
 * unchanged.
 */
/* Fraction of the largest conceivable value of the outer k-integral below
 * which a piece of it cannot matter. Mirrors NC_XCOR_KERNEL_INTEG_ABSTOL_FRAC
 * on the Levin side: splitting at the discontinuities creates sub-intervals
 * where some multipoles contribute nothing at all, and a purely relative
 * tolerance on a value of zero is unsatisfiable however far the quadrature
 * refines. */
#define NC_XCOR_KERNEL_CUBATURE_ABSTOL_FRAC 1.0e-16

/* Points used to size the integrand before integrating it. The estimate only
 * has to be right to an order of magnitude -- it scales a 1e-16 floor -- and
 * this costs a few dozen spline evaluations against a quadrature that will do
 * far more. */
#define NC_XCOR_KERNEL_CUBATURE_ABSTOL_NPT 64

/*
 * Sizes the outer k-integral of k^3 W1 W2 by sampling it, so a fraction of the
 * result is a floor no piece needs to resolve below.
 *
 * Bounding it analytically by k_max^3 max|W|^2 is useless: |W| peaks at small k
 * while k^3 peaks at large k, so the product overshoots the true integral by
 * orders of magnitude and the floor swallows the answer.
 */
static gdouble
_nc_xcor_kernel_cubature_abstol (NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2,
                                 gdouble k_min, gdouble k_max)
{
  const guint len       = nc_xcor_kernel_integrand_get_len (xclki1);
  const gdouble ln_lo   = log (k_min);
  const gdouble ln_hi   = log (k_max);
  const gdouble width   = ln_hi - ln_lo;
  gdouble *W1           = g_new (gdouble, len);
  gdouble *W2           = (xclki2 != NULL) ? g_new (gdouble, len) : NULL;
  gdouble max_integrand = 0.0;
  guint i, j;

  for (i = 0; i < NC_XCOR_KERNEL_CUBATURE_ABSTOL_NPT; i++)
  {
    const gdouble lnk = ln_lo + width * i / (NC_XCOR_KERNEL_CUBATURE_ABSTOL_NPT - 1.0);
    const gdouble k   = exp (lnk);
    const gdouble k3  = gsl_pow_3 (k);

    nc_xcor_kernel_integrand_eval (xclki1, k, W1);

    if (W2 != NULL)
      nc_xcor_kernel_integrand_eval (xclki2, k, W2);

    for (j = 0; j < len; j++)
    {
      const gdouble f = k3 * W1[j] * ((W2 != NULL) ? W2[j] : W1[j]);

      max_integrand = GSL_MAX (max_integrand, fabs (f));
    }
  }

  g_free (W1);
  g_free (W2);

  {
    const gdouble bound = max_integrand * width;

    if (!gsl_finite (bound) || (bound <= 0.0))
      return 0.0;

    return NC_XCOR_KERNEL_CUBATURE_ABSTOL_FRAC * bound;
  }
}

static void
_nc_xcor_kernel_eval_split (NcmIntegralND *xcor_int_nd,
                            NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2,
                            gdouble k_min, gdouble k_max,
                            NcmVector *vp, NcmVector *err)
{
  GArray *edges     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  NcmVector *lnk_lo = ncm_vector_new (1);
  NcmVector *lnk_hi = ncm_vector_new (1);
  guint i;

  ncm_integral_nd_set_abstol (xcor_int_nd,
                              _nc_xcor_kernel_cubature_abstol (xclki1, xclki2, k_min, k_max));

  {
    NcXcorKernelIntegrand *sources[2] = {xclki1, xclki2};
    guint si;

    for (si = 0; si < 2; si++)
    {
      GArray *bps = (sources[si] != NULL) ? nc_xcor_kernel_integrand_peek_breakpoints (sources[si]) : NULL;

      if (bps == NULL)
        continue;

      for (i = 0; i < bps->len; i++)
      {
        const gdouble bp = g_array_index (bps, gdouble, i);

        if ((bp > k_min) && (bp < k_max))
          g_array_append_val (edges, bp);
      }
    }
  }

  if (edges->len == 0)
  {
    ncm_vector_set (lnk_lo, 0, log (k_min));
    ncm_vector_set (lnk_hi, 0, log (k_max));
    ncm_integral_nd_eval (xcor_int_nd, lnk_lo, lnk_hi, vp, err);
  }
  else
  {
    NcmVector *vp_piece  = ncm_vector_new (ncm_vector_len (vp));
    NcmVector *err_piece = ncm_vector_new (ncm_vector_len (err));
    gdouble lo           = k_min;

    g_array_sort (edges, _nc_xcor_gdouble_cmp);
    ncm_vector_set_zero (vp);
    ncm_vector_set_zero (err);

    for (i = 0; i <= edges->len; i++)
    {
      const gdouble hi = (i < edges->len) ? g_array_index (edges, gdouble, i) : k_max;

      if (hi > lo)
      {
        ncm_vector_set (lnk_lo, 0, log (lo));
        ncm_vector_set (lnk_hi, 0, log (hi));
        ncm_integral_nd_eval (xcor_int_nd, lnk_lo, lnk_hi, vp_piece, err_piece);

        ncm_vector_add (vp, vp_piece);
        ncm_vector_add (err, err_piece);
      }

      lo = hi;
    }

    ncm_vector_free (vp_piece);
    ncm_vector_free (err_piece);
  }

  ncm_vector_free (lnk_lo);
  ncm_vector_free (lnk_hi);
  g_array_unref (edges);
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

    /* Evaluate all multipoles at once using pre-computed integrand */
    nc_xcor_kernel_integrand_eval (xcor_arg->xclki1, k, xcor_arg->W1);

    for (j = 0; j < fdim; j++)
    {
      const gdouble res = k3 * xcor_arg->W1[j] * xcor_arg->W1[j];

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

    /* Evaluate all multipoles at once using pre-computed integrands */
    nc_xcor_kernel_integrand_eval (xcor_arg->xclki1, k, xcor_arg->W1);
    nc_xcor_kernel_integrand_eval (xcor_arg->xclki2, k, xcor_arg->W2);

    for (j = 0; j < fdim; j++)
    {
      const gdouble res = k3 * xcor_arg->W1[j] * xcor_arg->W2[j];

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
 * Sets the multipole batch size for cubature methods.
 *
 */
void
nc_xcor_set_ell_batch_size (NcXcor *xc, const guint ell_batch_size)
{
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

    if (ret != GSL_SUCCESS)
      g_error ("_nc_xcor_limber_z_gsl: %s.", gsl_strerror (ret));

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

    ret = gsl_integration_qag (&F, log (k_min), log (k_max), 0.0, xc->reltol * 1.0e-2, NCM_INTEGRAL_PARTITION, 6, *w, &result, &err);

    if (ret != GSL_SUCCESS)
      g_error ("_nc_xcor_kernel_gsl: %s.", gsl_strerror (ret));

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

static void
_nc_xcor_kernel_cubature (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gboolean isauto, NcmVector *vp)
{
  _nc_xcor_check_kernel_tolerance (xc, xclk1);

  if (!isauto)
    _nc_xcor_check_kernel_tolerance (xc, xclk2);

  const guint size           = lmax - lmin + 1;
  const gdouble const_factor = 2.0 / (M_PI * gsl_pow_3 (xc->RH));
  const gint block           = xc->ell_batch_size;
  GArray *ells_array         = g_array_new (FALSE, FALSE, sizeof (gint));
  guint i;

  /* Outer k-bound is computed per block, after the block's integrand(s) are
   * built, from their own fitted domain (get_range) -- NOT the independent
   * Limber-band formula from nc_xcor_kernel_get_k_range() (see plan doc
   * dev-notes/xcor_ultralevin_batching_plan.md §3). */

  g_array_set_size (ells_array, size);

  for (i = 0; i < size; i++)
  {
    g_array_index (ells_array, gint, i) = lmin + i;
  }

  if (isauto)
  {
    NcXcorKernelAuto *xcor_kernel_auto = g_object_new (nc_xcor_kernel_auto_get_type (), NULL);
    NcXcorArg *xcor_arg                = &xcor_kernel_auto->data;
    NcmIntegralND *xcor_int_nd         = NCM_INTEGRAL_ND (xcor_kernel_auto);
    NcmVector *lnk_min                 = ncm_vector_new (1);
    NcmVector *lnk_max                 = ncm_vector_new (1);
    NcmVector *err                     = ncm_vector_new (size);

    xcor_arg->xc    = xc;
    xcor_arg->dist  = xc->dist;
    xcor_arg->ps    = xc->ps;
    xcor_arg->RH    = xc->RH;
    xcor_arg->cosmo = cosmo;
    xcor_arg->xclk1 = xclk1;
    xcor_arg->xclk2 = xclk1;

    ncm_integral_nd_set_reltol (xcor_int_nd, xc->reltol);
    ncm_integral_nd_set_abstol (xcor_int_nd, 0.0);
    ncm_integral_nd_set_method (xcor_int_nd, NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V);

    for (i = 0; i + block < size; i += block)
    {
      NcmVector *vp_i  = ncm_vector_get_subvector (vp, i, block);
      NcmVector *err_i = ncm_vector_get_subvector (err, i, block);
      gdouble k_min, k_max;

      xcor_arg->ells  = &g_array_index (ells_array, gint, i);
      xcor_arg->nells = block;

      /* Get vectorized kernel integrand for this batch, then read the outer
       * bound off its own fitted domain. */
      xcor_arg->xclki1 = nc_xcor_kernel_get_eval_vectorized (xclk1, cosmo, xcor_arg->ells[0], xcor_arg->ells[block - 1]);
      xcor_arg->W1     = g_new (gdouble, block);

      nc_xcor_kernel_integrand_get_range (xcor_arg->xclki1, &k_min, &k_max);
      ncm_vector_set (lnk_min, 0, log (k_min));
      ncm_vector_set (lnk_max, 0, log (k_max));

      _nc_xcor_kernel_eval_split (xcor_int_nd, xcor_arg->xclki1, xcor_arg->xclki2, k_min, k_max, vp_i, err_i);

      /* Clean up batch resources */
      g_free (xcor_arg->W1);
      nc_xcor_kernel_integrand_unref (xcor_arg->xclki1);
      xcor_arg->xclki1 = NULL;
      xcor_arg->W1     = NULL;

      /* Apply constant factor */
      ncm_vector_scale (vp_i, const_factor);

      ncm_vector_free (vp_i);
      ncm_vector_free (err_i);
    }

    if (i < size)
    {
      NcmVector *vp_i       = ncm_vector_get_subvector (vp, i, size - i);
      NcmVector *err_i      = ncm_vector_get_subvector (err, i, size - i);
      const guint remaining = size - i;
      gdouble k_min, k_max;

      xcor_arg->ells  = &g_array_index (ells_array, gint, i);
      xcor_arg->nells = remaining;

      /* Get vectorized kernel integrand for this batch, then read the outer
       * bound off its own fitted domain. */
      xcor_arg->xclki1 = nc_xcor_kernel_get_eval_vectorized (xclk1, cosmo, xcor_arg->ells[0], xcor_arg->ells[remaining - 1]);
      xcor_arg->W1     = g_new (gdouble, remaining);

      nc_xcor_kernel_integrand_get_range (xcor_arg->xclki1, &k_min, &k_max);
      ncm_vector_set (lnk_min, 0, log (k_min));
      ncm_vector_set (lnk_max, 0, log (k_max));

      _nc_xcor_kernel_eval_split (xcor_int_nd, xcor_arg->xclki1, xcor_arg->xclki2, k_min, k_max, vp_i, err_i);

      /* Clean up batch resources */
      g_free (xcor_arg->W1);
      nc_xcor_kernel_integrand_unref (xcor_arg->xclki1);
      xcor_arg->xclki1 = NULL;
      xcor_arg->W1     = NULL;

      /* Apply constant factor */
      ncm_vector_scale (vp_i, const_factor);

      ncm_vector_free (vp_i);
      ncm_vector_free (err_i);
    }

    ncm_vector_free (lnk_min);
    ncm_vector_free (lnk_max);
    ncm_vector_free (err);
    g_object_unref (xcor_kernel_auto);
  }
  else
  {
    NcXcorKernelCross *xcor_kernel_cross = g_object_new (nc_xcor_kernel_cross_get_type (), NULL);
    NcXcorArg *xcor_arg                  = &xcor_kernel_cross->data;
    NcmIntegralND *xcor_int_nd           = NCM_INTEGRAL_ND (xcor_kernel_cross);
    NcmVector *lnk_min                   = ncm_vector_new (1);
    NcmVector *lnk_max                   = ncm_vector_new (1);
    NcmVector *err                       = ncm_vector_new (size);

    xcor_arg->xc    = xc;
    xcor_arg->dist  = xc->dist;
    xcor_arg->ps    = xc->ps;
    xcor_arg->RH    = xc->RH;
    xcor_arg->cosmo = cosmo;
    xcor_arg->xclk1 = xclk1;
    xcor_arg->xclk2 = xclk2;

    ncm_integral_nd_set_reltol (xcor_int_nd, xc->reltol);
    ncm_integral_nd_set_abstol (xcor_int_nd, 0.0);
    ncm_integral_nd_set_method (xcor_int_nd, NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V);

    for (i = 0; i + block < size; i += block)
    {
      NcmVector *vp_i  = ncm_vector_get_subvector (vp, i, block);
      NcmVector *err_i = ncm_vector_get_subvector (err, i, block);
      gdouble k_min, k_max, k2_min, k2_max;

      xcor_arg->ells  = &g_array_index (ells_array, gint, i);
      xcor_arg->nells = block;

      /* Get vectorized kernel integrands for this batch, then read the
       * outer bound off their own fitted domains (intersection). */
      xcor_arg->xclki1 = nc_xcor_kernel_get_eval_vectorized (xclk1, cosmo, xcor_arg->ells[0], xcor_arg->ells[block - 1]);
      xcor_arg->xclki2 = nc_xcor_kernel_get_eval_vectorized (xclk2, cosmo, xcor_arg->ells[0], xcor_arg->ells[block - 1]);
      xcor_arg->W1     = g_new (gdouble, block);
      xcor_arg->W2     = g_new (gdouble, block);

      nc_xcor_kernel_integrand_get_range (xcor_arg->xclki1, &k_min, &k_max);
      nc_xcor_kernel_integrand_get_range (xcor_arg->xclki2, &k2_min, &k2_max);
      ncm_vector_set (lnk_min, 0, log (GSL_MAX (k_min, k2_min)));
      ncm_vector_set (lnk_max, 0, log (GSL_MIN (k_max, k2_max)));

      _nc_xcor_kernel_eval_split (xcor_int_nd, xcor_arg->xclki1, xcor_arg->xclki2, k_min, k_max, vp_i, err_i);

      /* Clean up batch resources */
      g_free (xcor_arg->W1);
      g_free (xcor_arg->W2);
      nc_xcor_kernel_integrand_unref (xcor_arg->xclki1);
      nc_xcor_kernel_integrand_unref (xcor_arg->xclki2);
      xcor_arg->xclki1 = NULL;
      xcor_arg->xclki2 = NULL;
      xcor_arg->W1     = NULL;
      xcor_arg->W2     = NULL;

      /* Apply constant factor */
      ncm_vector_scale (vp_i, const_factor);

      ncm_vector_free (vp_i);
      ncm_vector_free (err_i);
    }

    if (i < size)
    {
      NcmVector *vp_i       = ncm_vector_get_subvector (vp, i, size - i);
      NcmVector *err_i      = ncm_vector_get_subvector (err, i, size - i);
      const guint remaining = size - i;
      gdouble k_min, k_max, k2_min, k2_max;

      xcor_arg->ells  = &g_array_index (ells_array, gint, i);
      xcor_arg->nells = remaining;

      /* Get vectorized kernel integrands for this batch, then read the
       * outer bound off their own fitted domains (intersection). */
      xcor_arg->xclki1 = nc_xcor_kernel_get_eval_vectorized (xclk1, cosmo, xcor_arg->ells[0], xcor_arg->ells[remaining - 1]);
      xcor_arg->xclki2 = nc_xcor_kernel_get_eval_vectorized (xclk2, cosmo, xcor_arg->ells[0], xcor_arg->ells[remaining - 1]);
      xcor_arg->W1     = g_new (gdouble, remaining);
      xcor_arg->W2     = g_new (gdouble, remaining);

      nc_xcor_kernel_integrand_get_range (xcor_arg->xclki1, &k_min, &k_max);
      nc_xcor_kernel_integrand_get_range (xcor_arg->xclki2, &k2_min, &k2_max);
      ncm_vector_set (lnk_min, 0, log (GSL_MAX (k_min, k2_min)));
      ncm_vector_set (lnk_max, 0, log (GSL_MIN (k_max, k2_max)));

      _nc_xcor_kernel_eval_split (xcor_int_nd, xcor_arg->xclki1, xcor_arg->xclki2, k_min, k_max, vp_i, err_i);

      /* Clean up batch resources */
      g_free (xcor_arg->W1);
      g_free (xcor_arg->W2);
      nc_xcor_kernel_integrand_unref (xcor_arg->xclki1);
      nc_xcor_kernel_integrand_unref (xcor_arg->xclki2);
      xcor_arg->xclki1 = NULL;
      xcor_arg->xclki2 = NULL;
      xcor_arg->W1     = NULL;
      xcor_arg->W2     = NULL;

      /* Apply constant factor */
      ncm_vector_scale (vp_i, const_factor);

      ncm_vector_free (vp_i);
      ncm_vector_free (err_i);
    }

    ncm_vector_free (lnk_min);
    ncm_vector_free (lnk_max);
    ncm_vector_free (err);
    g_object_unref (xcor_kernel_cross);
  }

  g_array_unref (ells_array);
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
  GArray *ells_array         = g_array_new (FALSE, FALSE, sizeof (gint));
  NcmVector *lnk_min         = ncm_vector_new (1);
  NcmVector *lnk_max         = ncm_vector_new (1);
  NcmVector *err             = ncm_vector_new (size);
  gdouble k_min, k_max;
  guint i;

  if (ncm_vector_len (vp) != size)
    g_error ("_nc_xcor_kernel_integrate_block_cubature: vector size does not match multipole limits");

  g_array_set_size (ells_array, size);

  for (i = 0; i < size; i++)
  {
    g_array_index (ells_array, gint, i) = lmin + i;
  }

  nc_xcor_kernel_integrand_get_range (xclki1, &k_min, &k_max);

  if (!isauto)
  {
    gdouble k2_min, k2_max;

    nc_xcor_kernel_integrand_get_range (xclki2, &k2_min, &k2_max);
    k_min = GSL_MAX (k_min, k2_min);
    k_max = GSL_MIN (k_max, k2_max);
  }

  ncm_vector_set (lnk_min, 0, log (k_min));
  ncm_vector_set (lnk_max, 0, log (k_max));

  if (isauto)
  {
    NcXcorKernelAuto *xcor_kernel_auto = g_object_new (nc_xcor_kernel_auto_get_type (), NULL);
    NcXcorArg *xcor_arg                = &xcor_kernel_auto->data;
    NcmIntegralND *xcor_int_nd         = NCM_INTEGRAL_ND (xcor_kernel_auto);

    xcor_arg->xc     = xc;
    xcor_arg->RH     = xc->RH;
    xcor_arg->ells   = &g_array_index (ells_array, gint, 0);
    xcor_arg->nells  = size;
    xcor_arg->xclki1 = xclki1;
    xcor_arg->W1     = g_new (gdouble, size);

    ncm_integral_nd_set_reltol (xcor_int_nd, xc->reltol);
    ncm_integral_nd_set_abstol (xcor_int_nd, 0.0);
    ncm_integral_nd_set_method (xcor_int_nd, NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V);

    _nc_xcor_kernel_eval_split (xcor_int_nd, xclki1, isauto ? NULL : xclki2, k_min, k_max, vp, err);

    g_free (xcor_arg->W1);
    g_object_unref (xcor_kernel_auto);
  }
  else
  {
    NcXcorKernelCross *xcor_kernel_cross = g_object_new (nc_xcor_kernel_cross_get_type (), NULL);
    NcXcorArg *xcor_arg                  = &xcor_kernel_cross->data;
    NcmIntegralND *xcor_int_nd           = NCM_INTEGRAL_ND (xcor_kernel_cross);

    xcor_arg->xc     = xc;
    xcor_arg->RH     = xc->RH;
    xcor_arg->ells   = &g_array_index (ells_array, gint, 0);
    xcor_arg->nells  = size;
    xcor_arg->xclki1 = xclki1;
    xcor_arg->xclki2 = xclki2;
    xcor_arg->W1     = g_new (gdouble, size);
    xcor_arg->W2     = g_new (gdouble, size);

    ncm_integral_nd_set_reltol (xcor_int_nd, xc->reltol);
    ncm_integral_nd_set_abstol (xcor_int_nd, 0.0);
    ncm_integral_nd_set_method (xcor_int_nd, NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V);

    _nc_xcor_kernel_eval_split (xcor_int_nd, xclki1, isauto ? NULL : xclki2, k_min, k_max, vp, err);

    g_free (xcor_arg->W1);
    g_free (xcor_arg->W2);
    g_object_unref (xcor_kernel_cross);
  }

  ncm_vector_scale (vp, const_factor);

  g_array_unref (ells_array);
  ncm_vector_free (lnk_min);
  ncm_vector_free (lnk_max);
  ncm_vector_free (err);
}

/**
 * nc_xcor_compute:
 * @xc: a #NcXcor
 * @xclk1: a #NcXcorKernel
 * @xclk2: a #NcXcorKernel
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
 */
void
nc_xcor_compute (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, NcmVector *vp)
{
  const guint nell           = ncm_vector_len (vp);
  const gboolean isauto      = (xclk2 == NULL) || (xclk2 == xclk1);
  const gdouble const_factor = 1.0 / gsl_pow_3 (xc->RH);
  guint i;
  gdouble zmin, zmax, zmid;

  if (nell != lmax - lmin + 1)
    g_error ("nc_xcor_compute: vector size does not match multipole limits");

  nc_xcor_kernel_get_z_range (xclk1, &zmin, &zmax, &zmid);

  if (!isauto)
  {
    gdouble zmin_2, zmax_2, zmid_2;

    nc_xcor_kernel_get_z_range (xclk2, &zmin_2, &zmax_2, &zmid_2);
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
      case NC_XCOR_METHOD_LIMBER_Z_GSL:
        _nc_xcor_limber_z_gsl (xc, xclk1, xclk2, cosmo, lmin, lmax, zmin, zmax, isauto, vp);
        break;
      case NC_XCOR_METHOD_LIMBER_Z_CUBATURE:
        _nc_xcor_limber_z_cubature (xc, xclk1, xclk2, cosmo, lmin, lmax, zmin, zmax, isauto, vp);
        break;
      case NC_XCOR_METHOD_KERNEL_GSL:
        _nc_xcor_kernel_gsl (xc, xclk1, xclk2, cosmo, lmin, lmax, isauto, vp);

        return;

        break;
      case NC_XCOR_METHOD_KERNEL_CUBATURE:
        _nc_xcor_kernel_cubature (xc, xclk1, xclk2, cosmo, lmin, lmax, isauto, vp);

        return;

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

