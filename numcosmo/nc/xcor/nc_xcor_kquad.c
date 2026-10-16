/***************************************************************************
 *            nc_xcor_kquad.c
 *
 *  Sat August 29 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

/*
 * The outer k quadrature: everything the kernel-space method
 * %NC_XCOR_METHOD_KERNEL_EXACT does with a pair of k-space closures, the exact
 * integral of their panels with NcmSpectralPanelsProduct.
 *
 * Every entry point here takes closures, or the kernels to build them from,
 * and none of them knows about the Limber tier or the multipole policy that
 * chose it -- that lives in nc_xcor.c.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/integration/ncm_integrate.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/algebra/ncm_spectral_panels_product.h"
#include "nc/xcor/nc_xcor.h"
#include "ncm/specfunc/ncm_sbessel_integrator_levin.h"
#include "nc/xcor/nc_xcor_priv.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */

/* The product integrator of k^2 W_1 W_2 dk, one per process: it is reentrant */
static NcmSpectralPanelsProduct *
_nc_xcor_panels_product (void)
{
  static gsize init                     = 0;
  static NcmSpectralPanelsProduct *prod = NULL;

  if (g_once_init_enter (&init))
  {
    prod = ncm_spectral_panels_product_new (2);
    g_once_init_leave (&init, 1);
  }

  return prod;
}

/* The result and error of one multipole's product */
typedef struct _NcXcorTermProduct
{
  NcmVector *res;
  NcmVector *err;
} NcXcorTermProduct;

static gpointer
_nc_xcor_term_product_new (gpointer userdata)
{
  return g_new0 (NcXcorTermProduct, 1);
}

static void
_nc_xcor_term_product_free (gpointer p)
{
  NcXcorTermProduct *tp = p;

  ncm_vector_clear (&tp->res);
  ncm_vector_clear (&tp->err);
  g_free (tp);
}

/* The scratch space of the term products, one per concurrent caller */
static NcmMemoryPool *
_nc_xcor_term_product_pool (void)
{
  static gsize init        = 0;
  static NcmMemoryPool *mp = NULL;

  if (g_once_init_enter (&init))
  {
    mp = ncm_memory_pool_new (_nc_xcor_term_product_new, NULL, _nc_xcor_term_product_free);
    g_once_init_leave (&init, 1);
  }

  return mp;
}

/*
 * Integrates a pair of Chebyshev closures with NcmSpectralPanelsProduct: each
 * closure's panels in their own variables k / s_r. With one term per closure the
 * panel rows are the multipoles and the block is one product. Otherwise a closure of
 * several terms holds each multipole's window summed into one row in k
 * (nc_xcor_kernel_integrand_peek_window()), and the block is one product per
 * multipole: of two such rows, or of one and the other closure's row of that
 * multipole in its own variable.
 */
static void
_nc_xcor_kernel_integrate_block_panels (NcXcor *xc, NcXcorKernelIntegrand *xclki1,
                                        NcXcorKernelIntegrand *xclki2,
                                        gboolean isauto, NcmVector *vp, NcmVector *vp_err)
{
  const gdouble const_factor     = 2.0 / (M_PI * gsl_pow_3 (xc->RH));
  NcmSpectralPanelsProduct *prod = _nc_xcor_panels_product ();
  NcmSpectralPanels *sp1         = nc_xcor_kernel_integrand_peek_panels (xclki1);
  NcmSpectralPanels *sp2         = isauto ? sp1 : nc_xcor_kernel_integrand_peek_panels (xclki2);
  NcmVector *s1                  = nc_xcor_kernel_integrand_peek_scales (xclki1);
  NcmVector *s2                  = isauto ? s1 : nc_xcor_kernel_integrand_peek_scales (xclki2);
  const guint n_terms1           = nc_xcor_kernel_integrand_get_n_terms (xclki1);
  const guint n_terms2           = isauto ? n_terms1 : nc_xcor_kernel_integrand_get_n_terms (xclki2);

  if ((n_terms1 == 1) && (n_terms2 == 1))
  {
    ncm_spectral_panels_product_integrate (prod, sp1, s1, sp2, s2, vp, vp_err);
  }
  else
  {
    const guint len         = ncm_vector_len (vp);
    NcXcorTermProduct **ptr = ncm_memory_pool_get (_nc_xcor_term_product_pool ());
    NcXcorTermProduct *tp   = *ptr;
    NcmVector *err1         = (vp_err != NULL) ? tp->err : NULL;
    guint i;

    if (tp->res == NULL)
    {
      tp->res = ncm_vector_new (1);
      tp->err = ncm_vector_new (1);
    }

    for (i = 0; i < len; i++)
    {
      NcmSpectralPanels *P = nc_xcor_kernel_integrand_peek_window (xclki1, i);
      NcmSpectralPanels *Q = isauto ? P : nc_xcor_kernel_integrand_peek_window (xclki2, i);

      if ((P != NULL) && (Q != NULL))
        ncm_spectral_panels_product_integrate (prod, P, NULL, Q, NULL, tp->res, err1);
      else if (P != NULL)
        ncm_spectral_panels_product_integrate_rows (prod, P, 0, NULL, sp2, i, s2, tp->res, err1);
      else if (Q != NULL)
        ncm_spectral_panels_product_integrate_rows (prod, sp1, i, s1, Q, 0, NULL, tp->res, err1);
      else
        g_error ("_nc_xcor_kernel_integrate_block_panels: a closure of several terms has no summed windows.");

      ncm_vector_set (vp, i, ncm_vector_get (tp->res, 0));

      if (vp_err != NULL)
        ncm_vector_set (vp_err, i, ncm_vector_get (err1, 0));
    }

    ncm_memory_pool_return (ptr);
  }

  ncm_vector_scale (vp, const_factor);

  if (vp_err != NULL)
    ncm_vector_scale (vp_err, const_factor);
}

/*
 * Integrates a pair of closures on one multipole block exactly, as the product of
 * their panels.
 */
void
_nc_xcor_kernel_integrate_block_exact (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcmVector *vp, NcmVector *vp_err)
{
  const guint nell             = lmax - lmin + 1;
  NcXcorKernelIntegrand *side2 = isauto ? xclki1 : xclki2;

  if (ncm_vector_len (vp) != nell)
    g_error ("_nc_xcor_kernel_integrate_block_exact: vector size does not match multipole limits");

  if ((vp_err != NULL) && (ncm_vector_len (vp_err) != nell))
    g_error ("_nc_xcor_kernel_integrate_block_exact: error vector size does not match multipole limits");

  if ((nc_xcor_kernel_integrand_get_len (xclki1) != nell) || (nc_xcor_kernel_integrand_get_len (side2) != nell))
    g_error ("_nc_xcor_kernel_integrate_block_exact: the closures hold %u and %u multipoles, the block %u.",
             nc_xcor_kernel_integrand_get_len (xclki1), nc_xcor_kernel_integrand_get_len (side2), nell);

  _nc_xcor_kernel_integrate_block_panels (xc, xclki1, side2, isauto, vp, vp_err);
}

/*
 * Aborts when NcXcor:reltol is below the reltol of @sbi, the Levin integrator that builds
 * the kernel's k-space closure. The exact quadrature carries no tolerance of its own, so
 * this is a guard on the caller's stated intent: a C_ell cannot carry more precision than
 * its integrand.
 */
void
_nc_xcor_check_kernel_tolerance (NcXcor *xc, NcXcorKernel *xclk, NcmSBesselIntegrator *sbi)
{
  gdouble closure_reltol;

  if ((sbi == NULL) || !NCM_IS_SBESSEL_INTEGRATOR_LEVIN (sbi))
    return;

  {
    NcmSBesselIntegratorLevin *sbilv = NCM_SBESSEL_INTEGRATOR_LEVIN (sbi);

    closure_reltol = ncm_sbessel_integrator_levin_get_reltol (sbilv);
  }

  if (xc->reltol < closure_reltol)
    g_error ("_nc_xcor_check_kernel_tolerance: NcXcor:reltol is %.17g but kernel %s builds "
             "its k-space closure only to %.17g (the integrator's reltol). The outer "
             "integral cannot converge to more precision than the "
             "integrand carries. Loosen NcXcor:reltol to at least %.17g, or construct the "
             "integrator with a tighter reltol.",
             xc->reltol, G_OBJECT_TYPE_NAME (xclk), closure_reltol, closure_reltol);
}

/*
 * The kernel-space methods, one entry each. Adding a quadrature is a line here
 * plus its block function; nothing else selects on the method.
 */
static const NcXcorKQuad _nc_xcor_kquad_table[] = {
  { _nc_xcor_kernel_integrate_block_exact, TRUE, "NC_XCOR_METHOD_KERNEL_EXACT" },
};

const NcXcorKQuad *
_nc_xcor_kquad_for_method (NcXcorMethod meth)
{
  switch (meth)
  {
    case NC_XCOR_METHOD_KERNEL_EXACT:
      return &_nc_xcor_kquad_table[0];

    default:
      return NULL;
  }
}

void
_nc_xcor_kernel_space_run (NcXcor *xc, const NcXcorKQuad *kquad, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gboolean isauto, NcmVector *vp, NcmVector *vp_err)
{
  NcmSBesselIntegrator *sbi1 = nc_xcor_kernel_peek_integrator (xclk1);
  NcmSBesselIntegrator *sbi2 = nc_xcor_kernel_peek_integrator (xclk2);
  const guint size           = lmax - lmin + 1;
  const guint block          = xc->ell_batch_size;
  guint i;

  _nc_xcor_check_kernel_tolerance (xc, xclk1, sbi1);

  if (!isauto)
    _nc_xcor_check_kernel_tolerance (xc, xclk2, sbi2);

  /* Either kernel's integrator serves a kernel that carries none of its own. */
  if (sbi1 == NULL)
    sbi1 = sbi2;

  if (sbi2 == NULL)
    sbi2 = sbi1;

  /* Batched by NcXcor:ell-batch-size: one k-space closure per kernel per
   * batch, built here and handed to the same per-block integrator
   * #NcXcorSolver drives with closures of its own. The batching is not merely
   * an optimization -- a single closure spanning more than
   * NC_XCOR_KERNEL_MAX_ELL_BLOCK multipoles is a hard error in
   * nc_xcor_kernel_get_eval_vectorized_full(), so an unbatched sweep aborted
   * on any range wider than that. */
  for (i = 0; i < size; i += block)
  {
    const guint nells      = MIN (block, size - i);
    const guint block_lmin = lmin + i;
    const guint block_lmax = block_lmin + nells - 1;
    NcmVector *vp_i        = ncm_vector_get_subvector (vp, i, nells);
    NcmVector *vp_err_i    = ((vp_err != NULL) && kquad->has_err) ? ncm_vector_get_subvector (vp_err, i, nells) : NULL;
    NcXcorKernelIntegrand *xclki1;
    NcXcorKernelIntegrand *xclki2;

    xclki1 = nc_xcor_kernel_get_eval_vectorized_full (xclk1, cosmo, block_lmin, block_lmax, sbi1);
    xclki2 = isauto ? NULL : nc_xcor_kernel_get_eval_vectorized_full (xclk2, cosmo, block_lmin, block_lmax, sbi2);

    kquad->block (xc, xclki1, isauto ? xclki1 : xclki2, block_lmin, block_lmax, isauto, vp_i, vp_err_i);

    nc_xcor_kernel_integrand_unref (xclki1);

    if (xclki2 != NULL)
      nc_xcor_kernel_integrand_unref (xclki2);

    ncm_vector_free (vp_i);
    ncm_vector_clear (&vp_err_i);
  }
}

/**
 * nc_xcor_integrate_block:
 * @xc: a #NcXcor
 * @xclki1: a #NcXcorKernelIntegrand covering [@lmin, @lmax]
 * @xclki2: (nullable): the same for the second kernel, %NULL for an auto
 * @lmin: minimum multipole, matching @xclki1's (and @xclki2's) own range
 * @lmax: maximum multipole, matching @xclki1's (and @xclki2's) own range
 * @isauto: %TRUE for an auto spectrum, in which case @xclki2 is ignored
 * @meth: the quadrature to run, which need not be @xc's own #NcXcor:meth
 * @vp: a #NcmVector of length (@lmax - @lmin + 1), filled with the result
 * @vp_err: (nullable): a #NcmVector of the same length for the error
 * estimate, or %NULL
 *
 * Runs one kernel-space quadrature over closures the caller already holds,
 * instead of building them from kernels the way nc_xcor_compute() does.
 *
 * Everything else about a $C_\ell$ -- fitting $W_\ell(k)$, choosing the tier,
 * batching the multipoles -- is an order of magnitude more expensive than the
 * outer integral, so a timing taken through nc_xcor_compute() is a timing of
 * the closure. This is the entry point that separates them: build the pair
 * once with nc_xcor_kernel_get_eval_vectorized_full(), then integrate it here,
 * as #NcXcorSolver does with closures it shares across requests.
 *
 * @meth is given here rather than read from @xc for the same reason. It must
 * be the kernel-space method, %NC_XCOR_METHOD_KERNEL_EXACT. The redshift-space
 * Limber methods have no block quadrature, being a different approximation
 * rather than a different quadrature.
 *
 * @xc still supplies #NcXcor:reltol and the $2/(\pi R_H^3)$ factor, so
 * nc_xcor_prepare() must have been called for the cosmology the closures were
 * built at.
 *
 * @vp_err is filled only by the methods nc_xcor_method_has_error_estimate()
 * reports; the rest fill it with NaN, as nc_xcor_compute_full() does, since a
 * zero there would read as "no error". Ask that question before reading it.
 */
void
nc_xcor_integrate_block (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcXcorMethod meth, NcmVector *vp, NcmVector *vp_err)
{
  const NcXcorKQuad *kquad = _nc_xcor_kquad_for_method (meth);

  g_return_if_fail (NC_IS_XCOR (xc));
  g_return_if_fail (xclki1 != NULL);
  g_return_if_fail (isauto || (xclki2 != NULL));
  g_return_if_fail (lmax >= lmin);
  g_return_if_fail (vp != NULL);

  if (kquad == NULL)
    g_error ("nc_xcor_integrate_block: %s has no block quadrature; it is not a "
             "kernel-space method.",
             nc_xcor_method_get_name (meth));

  if ((vp_err != NULL) && !kquad->has_err)
    ncm_vector_set_all (vp_err, GSL_NAN);

  kquad->block (xc, xclki1, isauto ? xclki1 : xclki2, lmin, lmax, isauto, vp, kquad->has_err ? vp_err : NULL);
}

/**
 * nc_xcor_method_get_name:
 * @meth: a #NcXcorMethod
 *
 * Returns: (transfer none): @meth's enum name, for labelling and error messages.
 */
const gchar *
nc_xcor_method_get_name (NcXcorMethod meth)
{
  const NcXcorKQuad *kquad = _nc_xcor_kquad_for_method (meth);

  if (kquad != NULL)
    return kquad->name;

  switch (meth)
  {
    case NC_XCOR_METHOD_LIMBER_Z_GSL:
      return "NC_XCOR_METHOD_LIMBER_Z_GSL";

    case NC_XCOR_METHOD_LIMBER_Z_CUBATURE:
      return "NC_XCOR_METHOD_LIMBER_Z_CUBATURE";

    default:
      return "NC_XCOR_METHOD_INVALID";
  }
}

/**
 * nc_xcor_method_has_error_estimate:
 * @meth: a #NcXcorMethod
 *
 * Whether @meth fills the @vp_err of nc_xcor_compute_full() and
 * nc_xcor_integrate_block(). Every other method leaves it alone rather than
 * writing a zero that would read as "no error", so a caller that wants an
 * estimate has to ask this first.
 *
 * Returns: %TRUE when @meth reports an error estimate.
 */
gboolean
nc_xcor_method_has_error_estimate (NcXcorMethod meth)
{
  const NcXcorKQuad *kquad = _nc_xcor_kquad_for_method (meth);

  return (kquad != NULL) && kquad->has_err;
}

/**
 * nc_xcor_method_is_kernel_space:
 * @meth: a #NcXcorMethod
 *
 * Whether @meth integrates over $k$ with a pair of fitted closures, as opposed
 * to the redshift-space Limber tier.
 *
 * Returns: %TRUE for the kernel-space methods.
 */
gboolean
nc_xcor_method_is_kernel_space (NcXcorMethod meth)
{
  switch (meth)
  {
    case NC_XCOR_METHOD_KERNEL_EXACT:
      return TRUE;

    default:
      return FALSE;
  }
}

