/***************************************************************************
 *            ncm_sphere_xi.c
 *
 *  Fri Oct 16 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sphere_xi.c
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

/**
 * NcmSphereXi:
 *
 * Angular correlation function of a spin-0 field from its angular power spectrum.
 *
 * For a spin-0 field with angular power spectrum $C_\ell$, the angular correlation
 * function truncated at $\ell_\mathrm{max}$ (#NcmSphereXi:lmax) is
 * $$
 * \xi(\theta) = \sum_{\ell=0}^{\ell_\mathrm{max}} \frac{2\ell + 1}{4\pi} C_\ell P_\ell(\cos\theta).
 * $$
 * The object computes it in two forms.
 *
 * Bin averages. The average over $[\theta_b, \theta_{b+1}]$, weighted by solid angle, is
 * $$
 * \bar{\xi}_b = \sum_{\ell=0}^{\ell_\mathrm{max}} K_{b\ell} C_\ell,
 * \qquad
 * K_{b\ell} = \frac{2\ell + 1}{4\pi} \frac{1}{x_b - x_{b+1}} \int_{x_{b+1}}^{x_b} P_\ell(x)\,\mathrm{d}x,
 * $$
 * where $x_b = \cos\theta_b$ and $\theta_b$ are the entries of #NcmSphereXi:theta-edges.
 * ncm_sphere_xi_prepare() computes the kernel matrix $K$ once, and ncm_sphere_xi_apply()
 * evaluates $\bar{\xi}_b$ for several spectra in one matrix product. The integrals use the
 * antiderivative $-\sin^2\theta\,P_\ell^\prime(\cos\theta) / [\ell(\ell + 1)]$, with the
 * recurrences for $P_\ell$ and $P_\ell^\prime$ written in $1 - \cos\theta$, so the kernel
 * keeps double precision at small angles: against 40-digit values it agrees to
 * $7 \times 10^{-15}$ of $(2\ell + 1) / (4\pi)$ for $\ell \le 3 \times 10^4$, measured on
 * bins from zero to $\pi$.
 *
 * Chebyshev series. ncm_sphere_xi_to_chebyshev() converts the Legendre series to
 * $\xi(\theta) = \sum_k c_k T_k(\cos\theta) = \sum_k c_k \cos(k\theta)$, exactly for the
 * truncated sum, at a cost of $\ell_\mathrm{max}^2 / 4$ multiplications per spectrum.
 * ncm_sphere_xi_eval_chebyshev() evaluates the series at given angles and
 * ncm_sphere_xi_eval_chebyshev_grid() on the uniform grid $\theta_j = j\pi / \ell_\mathrm{max}$
 * by one discrete cosine transform.
 *
 * ncm_sphere_xi_apply() and ncm_sphere_xi_apply_vector() only read the kernel and may be
 * called concurrently on a prepared object; the other methods write to scratch space of the
 * object and may not.
 *
 * See the <a href="../../theory/ncm/sphere/sphere_xi.html">Angular Correlation Function</a>
 * page.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_cblas.h" /* This must be included before any gsl header */
#include "ncm/sphere/ncm_sphere_xi.h"
#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <string.h>
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmSphereXiPrivate
{
  NcmVector *theta_edges;
  guint lmax;
  NcmMatrix *K;
  GArray *G_lo;
  GArray *G_hi;
  gboolean prepared;
  GArray *Lambda;
  GArray *a_parity;
  GArray *w;
  GArray *cos_k;
  fftw_plan dct_plan;
  guint dct_n;
  guint dct_howmany;
  gdouble *dct_in;
  gdouble *dct_out;
} NcmSphereXiPrivate;

struct _NcmSphereXi
{
  GObject parent_instance;
};

enum
{
  PROP_0,
  PROP_THETA_EDGES,
  PROP_LMAX,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmSphereXi, ncm_sphere_xi, G_TYPE_OBJECT)

static void
ncm_sphere_xi_init (NcmSphereXi *sxi)
{
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);

  self->theta_edges = NULL;
  self->lmax        = 0;
  self->K           = NULL;
  self->G_lo        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->G_hi        = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->prepared    = FALSE;
  self->Lambda      = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->a_parity    = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->w           = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->cos_k       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->dct_plan    = NULL;
  self->dct_n       = 0;
  self->dct_howmany = 0;
  self->dct_in      = NULL;
  self->dct_out     = NULL;
}

static void
_ncm_sphere_xi_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSphereXi *sxi = NCM_SPHERE_XI (object);

  g_return_if_fail (NCM_IS_SPHERE_XI (object));

  switch (prop_id)
  {
    case PROP_THETA_EDGES:
      ncm_sphere_xi_set_theta_edges (sxi, g_value_get_object (value));
      break;
    case PROP_LMAX:
      ncm_sphere_xi_set_lmax (sxi, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sphere_xi_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSphereXi *sxi                = NCM_SPHERE_XI (object);
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);

  g_return_if_fail (NCM_IS_SPHERE_XI (object));

  switch (prop_id)
  {
    case PROP_THETA_EDGES:
      g_value_set_object (value, self->theta_edges);
      break;
    case PROP_LMAX:
      g_value_set_uint (value, self->lmax);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sphere_xi_dispose (GObject *object)
{
  NcmSphereXi *sxi                = NCM_SPHERE_XI (object);
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);

  ncm_vector_clear (&self->theta_edges);
  ncm_matrix_clear (&self->K);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sphere_xi_parent_class)->dispose (object);
}

static void
_ncm_sphere_xi_finalize (GObject *object)
{
  NcmSphereXi *sxi                = NCM_SPHERE_XI (object);
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);

  g_clear_pointer (&self->G_lo, g_array_unref);
  g_clear_pointer (&self->G_hi, g_array_unref);
  g_clear_pointer (&self->Lambda, g_array_unref);
  g_clear_pointer (&self->a_parity, g_array_unref);
  g_clear_pointer (&self->w, g_array_unref);
  g_clear_pointer (&self->cos_k, g_array_unref);
  g_clear_pointer (&self->dct_plan, fftw_destroy_plan);
  g_clear_pointer (&self->dct_in, fftw_free);
  g_clear_pointer (&self->dct_out, fftw_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sphere_xi_parent_class)->finalize (object);
}

static void
ncm_sphere_xi_class_init (NcmSphereXiClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_sphere_xi_set_property;
  object_class->get_property = &_ncm_sphere_xi_get_property;
  object_class->dispose      = &_ncm_sphere_xi_dispose;
  object_class->finalize     = &_ncm_sphere_xi_finalize;

  /**
   * NcmSphereXi:theta-edges:
   *
   * The bin edges $\theta_0 < \theta_1 < \dots < \theta_{n}$ in radians, with
   * $0 \le \theta_0$ and $\theta_n \le \pi$; the object keeps a copy.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_THETA_EDGES,
                                   g_param_spec_object ("theta-edges",
                                                        NULL,
                                                        "Bin edges in radians",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSphereXi:lmax:
   *
   * The largest multipole $\ell_\mathrm{max}$ of the sum; zero, the default, is unset.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LMAX,
                                   g_param_spec_uint ("lmax",
                                                      NULL,
                                                      "Largest multipole",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_sphere_xi_new:
 * @theta_edges: (nullable): a #NcmVector with the bin edges in radians
 * @lmax: the largest multipole
 *
 * Creates a new #NcmSphereXi for the bins @theta_edges and multipoles up to @lmax.
 *
 * Returns: (transfer full): a new #NcmSphereXi.
 */
NcmSphereXi *
ncm_sphere_xi_new (NcmVector *theta_edges, guint lmax)
{
  NcmSphereXi *sxi = g_object_new (NCM_TYPE_SPHERE_XI,
                                   "theta-edges", theta_edges,
                                   "lmax", lmax,
                                   NULL);

  return sxi;
}

/**
 * ncm_sphere_xi_ref:
 * @sxi: a #NcmSphereXi
 *
 * Increases the reference count of @sxi by one.
 *
 * Returns: (transfer full): @sxi.
 */
NcmSphereXi *
ncm_sphere_xi_ref (NcmSphereXi *sxi)
{
  return g_object_ref (sxi);
}

/**
 * ncm_sphere_xi_free:
 * @sxi: a #NcmSphereXi
 *
 * Decreases the reference count of @sxi by one.
 *
 */
void
ncm_sphere_xi_free (NcmSphereXi *sxi)
{
  g_object_unref (sxi);
}

/**
 * ncm_sphere_xi_clear:
 * @sxi: a #NcmSphereXi
 *
 * Decreases the reference count of *@sxi by one, and sets the pointer *@sxi to %NULL.
 *
 */
void
ncm_sphere_xi_clear (NcmSphereXi **sxi)
{
  g_clear_object (sxi);
}

/**
 * ncm_sphere_xi_set_theta_edges:
 * @sxi: a #NcmSphereXi
 * @theta_edges: (nullable): a #NcmVector with the bin edges in radians
 *
 * Sets #NcmSphereXi:theta-edges to a copy of @theta_edges. The kernel must then be
 * computed again with ncm_sphere_xi_prepare().
 *
 */
void
ncm_sphere_xi_set_theta_edges (NcmSphereXi *sxi, NcmVector *theta_edges)
{
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);

  ncm_vector_clear (&self->theta_edges);

  if (theta_edges != NULL)
    self->theta_edges = ncm_vector_dup (theta_edges);

  self->prepared = FALSE;
}

/**
 * ncm_sphere_xi_peek_theta_edges:
 * @sxi: a #NcmSphereXi
 *
 * Gets the bin edges of @sxi.
 *
 * Returns: (transfer none) (nullable): the #NcmVector with the bin edges.
 */
NcmVector *
ncm_sphere_xi_peek_theta_edges (NcmSphereXi *sxi)
{
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);

  return self->theta_edges;
}

/**
 * ncm_sphere_xi_set_lmax:
 * @sxi: a #NcmSphereXi
 * @lmax: the largest multipole
 *
 * Sets #NcmSphereXi:lmax. The bin kernel must then be computed again with
 * ncm_sphere_xi_prepare(); the Chebyshev conversion uses the new value directly.
 *
 */
void
ncm_sphere_xi_set_lmax (NcmSphereXi *sxi, guint lmax)
{
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);
  guint m;

  self->lmax     = lmax;
  self->prepared = FALSE;

  /* Lambda(m) = Gamma(m + 1/2) / Gamma(m + 1) for m = 0, ..., lmax, by its ratio recurrence */
  g_array_set_size (self->Lambda, lmax + 1);
  g_array_index (self->Lambda, gdouble, 0) = sqrt (ncm_c_pi ());

  for (m = 0; m < lmax; m++)
    g_array_index (self->Lambda, gdouble, m + 1) = g_array_index (self->Lambda, gdouble, m) * (m + 0.5) / (m + 1.0);
}

/**
 * ncm_sphere_xi_get_lmax:
 * @sxi: a #NcmSphereXi
 *
 * Gets #NcmSphereXi:lmax.
 *
 * Returns: the largest multipole.
 */
guint
ncm_sphere_xi_get_lmax (NcmSphereXi *sxi)
{
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);

  return self->lmax;
}

/**
 * ncm_sphere_xi_get_n_bins:
 * @sxi: a #NcmSphereXi
 *
 * Gets the number of bins, one less than the length of #NcmSphereXi:theta-edges, or zero
 * when the edges are unset.
 *
 * Returns: the number of bins.
 */
guint
ncm_sphere_xi_get_n_bins (NcmSphereXi *sxi)
{
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);

  if ((self->theta_edges == NULL) || (ncm_vector_len (self->theta_edges) < 2))
    return 0;

  return ncm_vector_len (self->theta_edges) - 1;
}

static void _ncm_sphere_xi_check_edges (NcmSphereXiPrivate * const self);
static void _ncm_sphere_xi_edge_integrals (const gdouble theta, const guint lmax, GArray *G_array);

/**
 * ncm_sphere_xi_prepare:
 * @sxi: a #NcmSphereXi
 *
 * Computes the kernel matrix $K_{b\ell}$, with one row per bin and one column per
 * multipole $\ell = 0, \dots, \ell_\mathrm{max}$. Aborts when #NcmSphereXi:lmax is
 * zero, when #NcmSphereXi:theta-edges is unset or has fewer than two entries, and when
 * the edges are not strictly increasing in $[0, \pi]$.
 *
 */
void
ncm_sphere_xi_prepare (NcmSphereXi *sxi)
{
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);
  const gdouble one_4pi           = 1.0 / (4.0 * ncm_c_pi ());
  guint n_bins, b, ell;

  if (self->lmax == 0)
    g_error ("ncm_sphere_xi_prepare: lmax is zero, set it with ncm_sphere_xi_set_lmax() first.");

  _ncm_sphere_xi_check_edges (self);

  n_bins = ncm_vector_len (self->theta_edges) - 1;

  if ((self->K == NULL) || (ncm_matrix_nrows (self->K) != n_bins) || (ncm_matrix_ncols (self->K) != self->lmax + 1))
  {
    ncm_matrix_clear (&self->K);
    self->K = ncm_matrix_new (n_bins, self->lmax + 1);
  }

  g_array_set_size (self->G_lo, self->lmax + 1);
  g_array_set_size (self->G_hi, self->lmax + 1);

  _ncm_sphere_xi_edge_integrals (ncm_vector_get (self->theta_edges, 0), self->lmax, self->G_lo);

  for (b = 0; b < n_bins; b++)
  {
    const gdouble theta_lo = ncm_vector_get (self->theta_edges, b);
    const gdouble theta_hi = ncm_vector_get (self->theta_edges, b + 1);
    const gdouble dx       = 2.0 * sin (0.5 * (theta_hi + theta_lo)) * sin (0.5 * (theta_hi - theta_lo));
    const gdouble *G_lo    = (const gdouble *) self->G_lo->data;
    const gdouble *G_hi;
    gdouble *K_b = ncm_matrix_ptr (self->K, b, 0);
    GArray *swap;

    _ncm_sphere_xi_edge_integrals (theta_hi, self->lmax, self->G_hi);
    G_hi = (const gdouble *) self->G_hi->data;

    K_b[0] = one_4pi;

    for (ell = 1; ell <= self->lmax; ell++)
      K_b[ell] = (2.0 * ell + 1.0) * one_4pi * (G_lo[ell] - G_hi[ell]) / dx;

    swap       = self->G_lo;
    self->G_lo = self->G_hi;
    self->G_hi = swap;
  }

  self->prepared = TRUE;
}

/* Aborts unless the edges hold at least two values strictly increasing in [0, pi] */
static void
_ncm_sphere_xi_check_edges (NcmSphereXiPrivate * const self)
{
  guint n_edges, i;

  if (self->theta_edges == NULL)
    g_error ("ncm_sphere_xi_prepare: theta-edges is unset, set it with ncm_sphere_xi_set_theta_edges() first.");

  n_edges = ncm_vector_len (self->theta_edges);

  if (n_edges < 2)
    g_error ("ncm_sphere_xi_prepare: theta-edges has %u entries, at least two are needed.", n_edges);

  if ((ncm_vector_get (self->theta_edges, 0) < 0.0) || (ncm_vector_get (self->theta_edges, n_edges - 1) > ncm_c_pi ()))
    g_error ("ncm_sphere_xi_prepare: theta-edges spans [% 22.15g, % 22.15g], outside [0, pi].",
             ncm_vector_get (self->theta_edges, 0), ncm_vector_get (self->theta_edges, n_edges - 1));

  for (i = 0; i + 1 < n_edges; i++)
  {
    if (!(ncm_vector_get (self->theta_edges, i) < ncm_vector_get (self->theta_edges, i + 1)))
      g_error ("ncm_sphere_xi_prepare: theta-edges is not strictly increasing at entry %u: % 22.15g, % 22.15g.",
               i, ncm_vector_get (self->theta_edges, i), ncm_vector_get (self->theta_edges, i + 1));
  }
}

/*
 * Writes G[ell] = -sin^2(theta) P'_ell(cos theta) / [ell (ell + 1)] for ell = 1, ..., lmax,
 * the antiderivative of P_ell(x) that vanishes at x = +-1; G[0] is not used and is set to
 * zero. The upward recurrences for P_ell and P'_ell run in u = 1 - x = 2 sin^2(theta / 2) at
 * min(theta, pi - theta), with D_ell = P_ell - P_{ell-1}, so x is never formed: near x = 1 a
 * change of x by one rounding unit changes P_ell by about ell^2 / 2 units. Above pi / 2,
 * G_ell(theta) = (-1)^(ell + 1) G_ell(pi - theta).
 */
static void
_ncm_sphere_xi_edge_integrals (const gdouble theta, const guint lmax, GArray *G_array)
{
  const gboolean reflect  = (theta > 0.5 * ncm_c_pi ());
  const gdouble theta_r   = reflect ? (ncm_c_pi () - theta) : theta;
  const gdouble sign_step = reflect ? -1.0 : 1.0;
  const gdouble sin_h     = sin (0.5 * theta_r);
  const gdouble u         = 2.0 * sin_h * sin_h;
  const gdouble s2        = u * (2.0 - u);
  gdouble * const G       = (gdouble *) G_array->data;
  gdouble sign            = 1.0;
  gdouble P_ell           = 1.0 - u;
  gdouble D_ell           = -u;
  gdouble dP_ell          = 1.0;
  guint ell;

  G[0] = 0.0;

  for (ell = 1; ell <= lmax; ell++)
  {
    const gdouble D_ellp1  = (ell * D_ell - (2.0 * ell + 1.0) * u * P_ell) / (ell + 1.0);
    const gdouble dP_ellp1 = dP_ell - u * dP_ell + (ell + 1.0) * P_ell;

    G[ell] = -sign * s2 * dP_ell / (ell * (ell + 1.0));

    P_ell  = P_ell + D_ellp1;
    D_ell  = D_ellp1;
    dP_ell = dP_ellp1;
    sign  *= sign_step;
  }
}

/**
 * ncm_sphere_xi_peek_kernel:
 * @sxi: a #NcmSphereXi
 *
 * Gets the kernel matrix $K_{b\ell}$ computed by ncm_sphere_xi_prepare(), with one row
 * per bin and one column per multipole $\ell = 0, \dots, \ell_\mathrm{max}$. Aborts when
 * the kernel is not prepared for the current bins and #NcmSphereXi:lmax.
 *
 * Returns: (transfer none): the kernel matrix.
 */
NcmMatrix *
ncm_sphere_xi_peek_kernel (NcmSphereXi *sxi)
{
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);

  if (!self->prepared)
    g_error ("ncm_sphere_xi_peek_kernel: the kernel is not prepared, call ncm_sphere_xi_prepare() first.");

  return self->K;
}

/**
 * ncm_sphere_xi_apply:
 * @sxi: a #NcmSphereXi
 * @C_ells: a #NcmMatrix with one spectrum per row, column $j$ holding $C_{\ell_\mathrm{min} + j}$
 * @lmin: the multipole of the first column of @C_ells
 * @xi: a #NcmMatrix for the result, one row per spectrum and one column per bin
 *
 * Computes $\bar{\xi}_b = \sum_\ell K_{b\ell} C_\ell$ for every row of @C_ells, the
 * multipoles below @lmin taken as zero, with one matrix product. Aborts when the kernel
 * is not prepared, when the last column of @C_ells is above #NcmSphereXi:lmax, and when @xi
 * does not have the rows of @C_ells and one column per bin.
 *
 */
void
ncm_sphere_xi_apply (NcmSphereXi *sxi, NcmMatrix *C_ells, const guint lmin, NcmMatrix *xi)
{
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);
  const guint n_spec              = ncm_matrix_nrows (C_ells);
  const guint n_ell               = ncm_matrix_ncols (C_ells);

  if (!self->prepared)
    g_error ("ncm_sphere_xi_apply: the kernel is not prepared, call ncm_sphere_xi_prepare() first.");

  if (lmin + n_ell > self->lmax + 1)
    g_error ("ncm_sphere_xi_apply: the spectra reach ell = %u, above lmax = %u.", lmin + n_ell - 1, self->lmax);

  if ((ncm_matrix_nrows (xi) != n_spec) || (ncm_matrix_ncols (xi) != ncm_matrix_nrows (self->K)))
    g_error ("ncm_sphere_xi_apply: xi is %u x %u, the result is %u x %u.",
             ncm_matrix_nrows (xi), ncm_matrix_ncols (xi), n_spec, ncm_matrix_nrows (self->K));

  cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasTrans, n_spec, ncm_matrix_nrows (self->K), n_ell,
               1.0, ncm_matrix_data (C_ells), ncm_matrix_tda (C_ells),
               ncm_matrix_ptr (self->K, 0, lmin), ncm_matrix_tda (self->K),
               0.0, ncm_matrix_data (xi), ncm_matrix_tda (xi));
}

/**
 * ncm_sphere_xi_apply_vector:
 * @sxi: a #NcmSphereXi
 * @C_ell: a #NcmVector with entry $j$ holding $C_{\ell_\mathrm{min} + j}$
 * @lmin: the multipole of the first entry of @C_ell
 * @xi: a #NcmVector for the result, one entry per bin
 *
 * Computes $\bar{\xi}_b = \sum_\ell K_{b\ell} C_\ell$ for one spectrum, the multipoles
 * below @lmin taken as zero. Aborts when the kernel is not prepared, when the last entry
 * of @C_ell is above #NcmSphereXi:lmax, and when @xi does not have one entry per bin.
 *
 */
void
ncm_sphere_xi_apply_vector (NcmSphereXi *sxi, NcmVector *C_ell, const guint lmin, NcmVector *xi)
{
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);
  const guint n_ell               = ncm_vector_len (C_ell);

  if (!self->prepared)
    g_error ("ncm_sphere_xi_apply_vector: the kernel is not prepared, call ncm_sphere_xi_prepare() first.");

  if (lmin + n_ell > self->lmax + 1)
    g_error ("ncm_sphere_xi_apply_vector: the spectrum reaches ell = %u, above lmax = %u.", lmin + n_ell - 1, self->lmax);

  if (ncm_vector_len (xi) != ncm_matrix_nrows (self->K))
    g_error ("ncm_sphere_xi_apply_vector: xi has %u entries, the result has %u.",
             ncm_vector_len (xi), ncm_matrix_nrows (self->K));

  cblas_dgemv (CblasRowMajor, CblasNoTrans, ncm_matrix_nrows (self->K), n_ell,
               1.0, ncm_matrix_ptr (self->K, 0, lmin), ncm_matrix_tda (self->K),
               ncm_vector_data (C_ell), ncm_vector_stride (C_ell),
               0.0, ncm_vector_data (xi), ncm_vector_stride (xi));
}

/**
 * ncm_sphere_xi_to_chebyshev:
 * @sxi: a #NcmSphereXi
 * @C_ells: a #NcmMatrix with one spectrum per row, column $j$ holding $C_{\ell_\mathrm{min} + j}$
 * @lmin: the multipole of the first column of @C_ells
 * @c: a #NcmMatrix for the result, one row per spectrum and $\ell_\mathrm{max} + 1$ columns
 *
 * Computes, for every row of @C_ells, the Chebyshev coefficients $c_k$ of the angular
 * correlation function truncated at $\ell_\mathrm{max}$,
 * $$
 * \xi(\theta) = \sum_{\ell=0}^{\ell_\mathrm{max}} \frac{2\ell + 1}{4\pi} C_\ell P_\ell(\cos\theta)
 * = \sum_{k=0}^{\ell_\mathrm{max}} c_k T_k(\cos\theta) = \sum_{k=0}^{\ell_\mathrm{max}} c_k \cos(k\theta),
 * $$
 * the multipoles below @lmin taken as zero. The cost is $\ell_\mathrm{max}^2 / 4$
 * multiplications per spectrum. Evaluate the result with ncm_sphere_xi_eval_chebyshev() or
 * ncm_sphere_xi_eval_chebyshev_grid(), which work in $\theta$; a Clenshaw sum in
 * $x = \cos\theta$, such as ncm_spectral_chebyshev_eval(), forms $x$ and loses accuracy at
 * small $\theta$: measured $2$ to $4 \times 10^{-11}$ of the largest $|\xi|$ on
 * $[2.5', 250']$ for $\ell_\mathrm{max}$ from 2000 to $10^5$. Aborts when #NcmSphereXi:lmax is
 * zero, when the last column of @C_ells is above #NcmSphereXi:lmax, and when @c does not
 * have the rows of @C_ells and $\ell_\mathrm{max} + 1$ columns. The scratch space belongs to
 * @sxi, so concurrent calls on one object are not allowed.
 *
 */
void
ncm_sphere_xi_to_chebyshev (NcmSphereXi *sxi, NcmMatrix *C_ells, const guint lmin, NcmMatrix *c)
{
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);
  const guint n_spec              = ncm_matrix_nrows (C_ells);
  const guint n_ell               = ncm_matrix_ncols (C_ells);
  const guint h                   = self->lmax / 2 + 1;
  const gdouble one_4pi           = 1.0 / (4.0 * ncm_c_pi ());
  const gdouble *Lambda           = (const gdouble *) self->Lambda->data;
  gdouble *a_parity, *w;
  guint r, j, k;

  if (self->lmax == 0)
    g_error ("ncm_sphere_xi_to_chebyshev: lmax is zero, set it with ncm_sphere_xi_set_lmax() first.");

  if (lmin + n_ell > self->lmax + 1)
    g_error ("ncm_sphere_xi_to_chebyshev: the spectra reach ell = %u, above lmax = %u.", lmin + n_ell - 1, self->lmax);

  if ((ncm_matrix_nrows (c) != n_spec) || (ncm_matrix_ncols (c) != self->lmax + 1))
    g_error ("ncm_sphere_xi_to_chebyshev: c is %u x %u, the result is %u x %u.",
             ncm_matrix_nrows (c), ncm_matrix_ncols (c), n_spec, self->lmax + 1);

  /*
   * a_ell = (2 ell + 1) / (4 pi) C_ell, split by the parity q of ell: the block of parity q
   * holds row r at r h and ell at column ell / 2, zero outside [lmin, lmin + n_ell).
   */
  g_array_set_size (self->a_parity, 2 * n_spec * h);
  g_array_set_size (self->w, h);
  a_parity = (gdouble *) self->a_parity->data;
  w        = (gdouble *) self->w->data;
  memset (a_parity, 0, sizeof (gdouble) * 2 * n_spec * h);

  for (r = 0; r < n_spec; r++)
  {
    for (j = 0; j < n_ell; j++)
    {
      const guint ell = lmin + j;

      a_parity[(ell % 2) * n_spec * h + r * h + ell / 2] = (2.0 * ell + 1.0) * one_4pi * ncm_matrix_get (C_ells, r, j);
    }
  }

  /*
   * c_k = (2 - delta_k0) / pi sum_i Lambda(i) Lambda(k + i) a_{k + 2i}: the multipoles
   * k + 2i have the parity of k and sit at columns k / 2 + i of its block.
   */
  for (k = 0; k <= self->lmax; k++)
  {
    const guint m       = (self->lmax - k) / 2 + 1;
    const gdouble scale = ((k == 0) ? 1.0 : 2.0) / ncm_c_pi ();
    guint i;

    for (i = 0; i < m; i++)
      w[i] = scale * Lambda[i] * Lambda[k + i];

    cblas_dgemv (CblasRowMajor, CblasNoTrans, n_spec, m,
                 1.0, &a_parity[(k % 2) * n_spec * h + k / 2], h,
                 w, 1,
                 0.0, ncm_matrix_ptr (c, 0, k), ncm_matrix_tda (c));
  }
}

/**
 * ncm_sphere_xi_eval_chebyshev:
 * @sxi: a #NcmSphereXi
 * @c: a #NcmMatrix of Chebyshev coefficients, one series per row
 * @theta: a #NcmVector of angles in radians
 * @xi: a #NcmMatrix for the result, one row per series and one column per angle
 *
 * Evaluates $\xi(\theta) = \sum_k c_k \cos(k\theta)$ for every row of @c at every entry of
 * @theta, with $\cos(k\theta)$ from the rotation recurrence in $k$. For the coefficients of
 * ncm_sphere_xi_to_chebyshev(), the result agrees with 30-digit sums of the Legendre series
 * to $2 \times 10^{-14}$ of the largest $|\xi|$ on $[2.5', 250']$, measured for
 * $\ell_\mathrm{max} = 2000$ and $2 \times 10^4$. Aborts when @xi does not have the rows of
 * @c and one column per angle. The scratch space belongs to @sxi, so concurrent calls on one
 * object are not allowed.
 *
 */
void
ncm_sphere_xi_eval_chebyshev (NcmSphereXi *sxi, NcmMatrix *c, NcmVector *theta, NcmMatrix *xi)
{
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);
  const guint n_spec              = ncm_matrix_nrows (c);
  const guint n                   = ncm_matrix_ncols (c);
  const guint n_theta             = ncm_vector_len (theta);
  gdouble *cos_k;
  guint t, k;

  if ((ncm_matrix_nrows (xi) != n_spec) || (ncm_matrix_ncols (xi) != n_theta))
    g_error ("ncm_sphere_xi_eval_chebyshev: xi is %u x %u, the result is %u x %u.",
             ncm_matrix_nrows (xi), ncm_matrix_ncols (xi), n_spec, n_theta);

  g_array_set_size (self->cos_k, n);
  cos_k = (gdouble *) self->cos_k->data;

  for (t = 0; t < n_theta; t++)
  {
    const gdouble th     = ncm_vector_get (theta, t);
    const gdouble cos_th = cos (th);
    const gdouble sin_th = sin (th);
    gdouble ck           = 1.0;
    gdouble sk           = 0.0;

    for (k = 0; k < n; k++)
    {
      const gdouble ck1 = ck * cos_th - sk * sin_th;

      cos_k[k] = ck;
      sk       = sk * cos_th + ck * sin_th;
      ck       = ck1;
    }

    cblas_dgemv (CblasRowMajor, CblasNoTrans, n_spec, n,
                 1.0, ncm_matrix_data (c), ncm_matrix_tda (c),
                 cos_k, 1,
                 0.0, ncm_matrix_ptr (xi, 0, t), ncm_matrix_tda (xi));
  }
}

static void _ncm_sphere_xi_dct_plan (NcmSphereXiPrivate * const self, const guint n, const guint howmany);

/**
 * ncm_sphere_xi_eval_chebyshev_grid:
 * @sxi: a #NcmSphereXi
 * @c: a #NcmMatrix of Chebyshev coefficients, one series per row, $n \ge 2$ columns
 * @xi: a #NcmMatrix for the result, the shape of @c
 *
 * Evaluates $\xi(\theta_j) = \sum_k c_k \cos(k\theta_j)$ for every row of @c on the $n$
 * angles $\theta_j = j\pi / (n - 1)$, $j = 0, \dots, n - 1$, with one type-I discrete cosine
 * transform per row. For the coefficients of ncm_sphere_xi_to_chebyshev(), $n - 1$ is
 * $\ell_\mathrm{max}$. The FFTW plan is made on the first call for a given shape and kept.
 * Aborts when @c has fewer than two columns and when @xi does not have the shape of @c. The
 * plan and the buffers belong to @sxi, so concurrent calls on one object are not allowed.
 *
 */
void
ncm_sphere_xi_eval_chebyshev_grid (NcmSphereXi *sxi, NcmMatrix *c, NcmMatrix *xi)
{
  NcmSphereXiPrivate * const self = ncm_sphere_xi_get_instance_private (sxi);
  const guint n_spec              = ncm_matrix_nrows (c);
  const guint n                   = ncm_matrix_ncols (c);
  guint r, k;

  if (n < 2)
    g_error ("ncm_sphere_xi_eval_chebyshev_grid: c has %u columns, at least two are needed.", n);

  if ((ncm_matrix_nrows (xi) != n_spec) || (ncm_matrix_ncols (xi) != n))
    g_error ("ncm_sphere_xi_eval_chebyshev_grid: xi is %u x %u, the result is %u x %u.",
             ncm_matrix_nrows (xi), ncm_matrix_ncols (xi), n_spec, n);

  _ncm_sphere_xi_dct_plan (self, n, n_spec);

  /*
   * FFTW's REDFT00 of X gives Y_j = X_0 + (-1)^j X_{n-1} + 2 sum_{k=1}^{n-2} X_k cos(pi j k / (n - 1)),
   * so the interior coefficients enter halved. The n_spec transforms are interleaved.
   */
  for (k = 0; k < n; k++)
  {
    const gdouble f = ((k == 0) || (k == n - 1)) ? 1.0 : 0.5;

    for (r = 0; r < n_spec; r++)
      self->dct_in[k * n_spec + r] = f * ncm_matrix_get (c, r, k);
  }

  fftw_execute_r2r (self->dct_plan, self->dct_in, self->dct_out);

  for (r = 0; r < n_spec; r++)
    for (k = 0; k < n; k++)
      ncm_matrix_set (xi, r, k, self->dct_out[k * n_spec + r]);
}

/*
 * Keeps one REDFT00 plan of size n for howmany interleaved transforms, with its two buffers,
 * and replaces them when the shape changes. The plan is made on the buffers before they hold
 * any data, and executed through the new-array interface on the same buffers.
 */
static void
_ncm_sphere_xi_dct_plan (NcmSphereXiPrivate * const self, const guint n, const guint howmany)
{
  const fftw_r2r_kind kind[] = { FFTW_REDFT00 };
  const gint nn[]            = { (gint) n };
  gboolean first;

  if ((self->dct_plan != NULL) && (self->dct_n == n) && (self->dct_howmany == howmany))
    return;

  g_clear_pointer (&self->dct_plan, fftw_destroy_plan);
  g_clear_pointer (&self->dct_in, fftw_free);
  g_clear_pointer (&self->dct_out, fftw_free);

  self->dct_in      = fftw_alloc_real (n * howmany);
  self->dct_out     = fftw_alloc_real (n * howmany);
  self->dct_n       = n;
  self->dct_howmany = howmany;

  first          = ncm_cfg_fftw_plan_begin ("ncm_sphere_xi_redft00_%u_%u", n, howmany);
  self->dct_plan = fftw_plan_many_r2r (1, nn, (gint) howmany, self->dct_in, NULL, (gint) howmany, 1,
                                       self->dct_out, NULL, (gint) howmany, 1, kind, ncm_cfg_get_fftw_default_flag ());
  ncm_cfg_fftw_plan_end (first);
}

