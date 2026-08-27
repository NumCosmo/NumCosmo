/***************************************************************************
 *            nc_xcor_ssc_sij.c
 *
 *  Thu August 13 2026
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

/**
 * NcXcorSSCSij:
 *
 * Super-sample covariance $S_{ij}$ matrix for a set of top-hat redshift bins.
 *
 * Each redshift bin becomes a #NcXcorKernelClusterTophat whose radial window is
 * normalized to unit integral in comoving distance, so the angular power
 * spectrum $C^{ij}_\ell$ of the pair is that of the volume-averaged matter
 * density contrast. Given the angular power spectrum $C^{\rm mask}_\ell$ of the
 * survey footprint,
 *
 * $$S_{ij} = \frac{1}{4\pi C^{\rm mask}_0}
 *     \sum_{\ell=0}^{\ell_{\rm max}} (2\ell+1) C^{\rm mask}_\ell C^{ij}_\ell,$$
 *
 * which reduces to the full-sky $S_{ij} = C^{ij}_0 / 4\pi$ for the trivial mask
 * $C^{\rm mask}_\ell = 4\pi \delta_{\ell 0}$ --- the default here, also
 * available as nc_xcor_ssc_sij_mask_cl_fullsky(). See
 * <a href="../../theory/ssc.html">Super-sample covariance</a> for the
 * derivation and the accuracy study.
 *
 * The mask spectrum does not depend on cosmology, so it is supplied once, by
 * nc_xcor_ssc_sij_set_mask_cl(); #NcmSphereMap computes it from a HEALPix
 * footprint. Only the $C^{ij}_\ell$ are recomputed per cosmology, which is what
 * makes this object cheap enough to call once per likelihood step.
 *
 * Every kernel is put in permanent non-Limber mode (`l_limber = -1`): the
 * Limber approximation is meaningless at the low multipoles that dominate
 * $S_{ij}$, and it makes the cross spectrum of two disjoint bins vanish.
 *
 * The default quadrature is %NC_XCOR_METHOD_KERNEL_EXACT, which needs no
 * tolerance and cannot fail to converge. The adaptive alternatives target a
 * tolerance the integrand may not support and abort when they cannot reach it,
 * which would kill a Monte Carlo chain mid-flight.
 *
 * A single #NcXcorSolver is built once and reused across cosmologies, so the
 * per-block spherical Bessel factorizations are paid for only on the first
 * nc_xcor_ssc_sij_prepare().
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/xcor/nc_xcor_ssc_sij.h"
#include "nc/xcor/nc_xcor_kernel_cluster_tophat.h"
#include "nc/xcor/nc_xcor_solver.h"
#include "ncm/core/ncm_c.h"
#include "ncm/model/ncm_model_ctrl.h"
#include "ncm/specfunc/ncm_sbessel_integrator_levin.h"
#include "nc_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* Empirical sweet spot for nc_xcor_solver_plan_blocks(), from
 * dev-notes/xcor_ultralevin_batching_plan.md section 1.3. */
#define NC_XCOR_SSC_SIJ_DEFAULT_BLOCK_SIZE (8)

/* One order tighter than #NcXcorKernel's own 1.0e-4 default. This, not
 * `reltol`, is what limits the accuracy of the off-diagonal S_ij: they are a
 * small residual of a large cancellation, four orders of magnitude below the
 * diagonal for well-separated bins.
 *
 * Deliberately different from NC_XCOR_SSC_SIJ_DEFAULT_RELTOL below. Equal
 * values are the one setting that makes the p-adaptive cubature run out of
 * Clenshaw-Curtis levels, because the refinement then stops exactly at the
 * level the outer rule is trying to resolve. The offset is taken upwards
 * because this object rebuilds S at every likelihood step, where tightening
 * costs ~2x per rebuild for accuracy no forecast can use.
 *
 * Must be kept equal to DEFAULT_SCALED_ABSTOL in numcosmo_py/ssc.py, which is
 * what the frozen path uses: the two are documented to differ only in whether
 * S_ij follows the cosmology, not in how it is computed. */
#define NC_XCOR_SSC_SIJ_DEFAULT_SCALED_ABSTOL (1.0e-5)

#define NC_XCOR_SSC_SIJ_DEFAULT_RELTOL (1.0e-6)

enum
{
  PROP_0,
  PROP_DIST,
  PROP_POWSPEC,
  PROP_Z_EDGES,
  PROP_MASK_CL,
  PROP_AREA,
  PROP_METHOD,
  PROP_BLOCK_SIZE,
  PROP_RELTOL,
  PROP_SCALED_ABSTOL,
  PROP_SIZE,
};

struct _NcXcorSSCSij
{
  /*< private > */
  GObject parent_instance;

  NcDistance *dist;
  NcmPowspec *ps;
  NcmVector *z_edges;
  NcmVector *mask_cl;
  gdouble area;

  NcXcorMethod method;
  guint block_size;
  gdouble reltol;
  gdouble scaled_abstol;

  GPtrArray *kernels; /* element-type NcXcorKernel*, one per bin, owned refs */
  NcXcor *xcor;

  /* Built on the first prepare() and reused: the solver pins one spherical
   * Bessel factorization per ell-block for its whole life, so rebuilding it
   * per cosmology would throw away the only expensive part that is not
   * cosmology dependent. Dropped whenever a setting that changes the block
   * layout or the requested ell range changes. */
  NcXcorSolver *solver;

  NcmMatrix *sij;
  NcmModelCtrl *ctrl;
  gboolean prepared;

  /* `mask-cl` and `area` are mutually exclusive, and GObject sets construct
   * properties in an unspecified order. So the setters only enforce the
   * exclusion once construction is over; until then they just record what
   * they were given and constructed() settles it. */
  gboolean constructed;
};

G_DEFINE_TYPE (NcXcorSSCSij, nc_xcor_ssc_sij, G_TYPE_OBJECT)

static void
nc_xcor_ssc_sij_init (NcXcorSSCSij *ssc_sij)
{
  ssc_sij->dist    = NULL;
  ssc_sij->ps      = NULL;
  ssc_sij->z_edges = NULL;
  ssc_sij->mask_cl = NULL;
  ssc_sij->area    = 0.0;

  ssc_sij->method        = NC_XCOR_METHOD_KERNEL_EXACT;
  ssc_sij->block_size    = NC_XCOR_SSC_SIJ_DEFAULT_BLOCK_SIZE;
  ssc_sij->reltol        = NC_XCOR_SSC_SIJ_DEFAULT_RELTOL;
  ssc_sij->scaled_abstol = NC_XCOR_SSC_SIJ_DEFAULT_SCALED_ABSTOL;

  ssc_sij->kernels     = g_ptr_array_new_with_free_func ((GDestroyNotify) nc_xcor_kernel_free);
  ssc_sij->xcor        = NULL;
  ssc_sij->solver      = NULL;
  ssc_sij->sij         = NULL;
  ssc_sij->ctrl        = ncm_model_ctrl_new (NULL);
  ssc_sij->prepared    = FALSE;
  ssc_sij->constructed = FALSE;
}

/* Any change to the ell range, the block layout or the kernels' tolerances
 * invalidates both the planned solver and the last computed matrix. */
static void
_nc_xcor_ssc_sij_invalidate (NcXcorSSCSij *ssc_sij)
{
  nc_xcor_solver_clear (&ssc_sij->solver);
  ssc_sij->prepared = FALSE;
}

static void
_nc_xcor_ssc_sij_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcorSSCSij *ssc_sij = NC_XCOR_SSC_SIJ (object);

  g_return_if_fail (NC_IS_XCOR_SSC_SIJ (object));

  switch (prop_id)
  {
    case PROP_DIST:
      ssc_sij->dist = g_value_dup_object (value);
      break;
    case PROP_POWSPEC:
      ssc_sij->ps = g_value_dup_object (value);
      break;
    case PROP_Z_EDGES:
      ssc_sij->z_edges = g_value_dup_object (value);
      break;
    case PROP_MASK_CL:
      nc_xcor_ssc_sij_set_mask_cl (ssc_sij, g_value_get_object (value));
      break;
    case PROP_AREA:
      nc_xcor_ssc_sij_set_area (ssc_sij, g_value_get_double (value));
      break;
    case PROP_METHOD:
      nc_xcor_ssc_sij_set_method (ssc_sij, g_value_get_enum (value));
      break;
    case PROP_BLOCK_SIZE:
      nc_xcor_ssc_sij_set_block_size (ssc_sij, g_value_get_uint (value));
      break;
    case PROP_RELTOL:
      nc_xcor_ssc_sij_set_reltol (ssc_sij, g_value_get_double (value));
      break;
    case PROP_SCALED_ABSTOL:
      nc_xcor_ssc_sij_set_scaled_abstol (ssc_sij, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_xcor_ssc_sij_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcorSSCSij *ssc_sij = NC_XCOR_SSC_SIJ (object);

  g_return_if_fail (NC_IS_XCOR_SSC_SIJ (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, ssc_sij->dist);
      break;
    case PROP_POWSPEC:
      g_value_set_object (value, ssc_sij->ps);
      break;
    case PROP_Z_EDGES:
      g_value_set_object (value, ssc_sij->z_edges);
      break;
    case PROP_MASK_CL:
      g_value_set_object (value, ssc_sij->area > 0.0 ? NULL : ssc_sij->mask_cl);
      break;
    case PROP_AREA:
      g_value_set_double (value, ssc_sij->area);
      break;
    case PROP_METHOD:
      g_value_set_enum (value, ssc_sij->method);
      break;
    case PROP_BLOCK_SIZE:
      g_value_set_uint (value, ssc_sij->block_size);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, ssc_sij->reltol);
      break;
    case PROP_SCALED_ABSTOL:
      g_value_set_double (value, ssc_sij->scaled_abstol);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

/* Builds the per-bin kernels and the #NcXcor. Runs once every construct
 * property has been set, so `z-edges`, `dist` and `powspec` can be given in
 * any order and the tolerances are already the final ones.
 *
 * An instance built without them stays empty rather than aborting: a plain
 * g_object_new() with no properties has to survive, since the serialization
 * registry and the introspection tooling both instantiate every registered
 * type that way. Such an instance has no bins, and _check_ready() rejects it
 * on the first prepare. */
static void
_nc_xcor_ssc_sij_constructed (GObject *object)
{
  NcXcorSSCSij *ssc_sij = NC_XCOR_SSC_SIJ (object);
  guint nbins;
  guint i;

  /* Chain up : start */
  G_OBJECT_CLASS (nc_xcor_ssc_sij_parent_class)->constructed (object);

  if ((ssc_sij->mask_cl != NULL) && (ssc_sij->area > 0.0))
    g_error ("_nc_xcor_ssc_sij_constructed: `mask-cl' and `area' are mutually "
             "exclusive; a mask carries its own sky fraction.");

  if (ssc_sij->mask_cl == NULL)
    ssc_sij->mask_cl = nc_xcor_ssc_sij_mask_cl_fullsky ();

  ssc_sij->constructed = TRUE;

  if ((ssc_sij->dist == NULL) || (ssc_sij->ps == NULL) || (ssc_sij->z_edges == NULL))
    return;

  if (ncm_vector_len (ssc_sij->z_edges) < 2)
    g_error ("_nc_xcor_ssc_sij_constructed: z-edges must contain at least two edges, got %u.",
             ncm_vector_len (ssc_sij->z_edges));

  nbins = ncm_vector_len (ssc_sij->z_edges) - 1;

  for (i = 0; i < nbins; i++)
  {
    const gdouble z_lower = ncm_vector_get (ssc_sij->z_edges, i + 0);
    const gdouble z_upper = ncm_vector_get (ssc_sij->z_edges, i + 1);
    NcXcorKernel *kernel;

    if (z_upper <= z_lower)
      g_error ("_nc_xcor_ssc_sij_constructed: z-edges must be strictly increasing, "
               "got %.17g at index %u followed by %.17g.", z_lower, i, z_upper);

    /* #NcXcorKernel refuses non-Limber mode unless it carries an integrator of
     * its own, so one is handed over at construction. The solver ignores it --
     * it copies its own per-block integrators -- but the kernel needs it to
     * accept l_limber = -1. */
    NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (0, ssc_sij->block_size));

    kernel = NC_XCOR_KERNEL (nc_xcor_kernel_cluster_tophat_new_full (ssc_sij->dist, ssc_sij->ps, z_lower, z_upper, sbi));
    ncm_sbessel_integrator_free (sbi);

    /* The Limber approximation is meaningless at the low multipoles dominating
     * S_ij, and makes the cross spectrum of two disjoint bins vanish. */
    nc_xcor_kernel_set_l_limber (kernel, -1);
    nc_xcor_kernel_set_reltol (kernel, ssc_sij->reltol);
    nc_xcor_kernel_set_scaled_abstol (kernel, ssc_sij->scaled_abstol);

    g_ptr_array_add (ssc_sij->kernels, kernel);
  }

  ssc_sij->xcor = nc_xcor_new (ssc_sij->dist, ssc_sij->ps, ssc_sij->method);
  nc_xcor_set_reltol (ssc_sij->xcor, ssc_sij->reltol);

  ssc_sij->sij = ncm_matrix_new (nbins, nbins);
  ncm_matrix_set_zero (ssc_sij->sij);
}

/* Guards the entry points that need the objects `constructed` may have found
 * missing. */
static void
_nc_xcor_ssc_sij_check_ready (NcXcorSSCSij *ssc_sij)
{
  if (ssc_sij->xcor == NULL)
    g_error ("NcXcorSSCSij: cannot be used before `dist', `powspec' and `z-edges' "
             "are set; build it with nc_xcor_ssc_sij_new().");
}

static void
_nc_xcor_ssc_sij_dispose (GObject *object)
{
  NcXcorSSCSij *ssc_sij = NC_XCOR_SSC_SIJ (object);

  nc_distance_clear (&ssc_sij->dist);
  ncm_powspec_clear (&ssc_sij->ps);
  ncm_vector_clear (&ssc_sij->z_edges);
  ncm_vector_clear (&ssc_sij->mask_cl);
  ncm_matrix_clear (&ssc_sij->sij);
  ncm_model_ctrl_clear (&ssc_sij->ctrl);
  nc_xcor_solver_clear (&ssc_sij->solver);
  nc_xcor_clear (&ssc_sij->xcor);

  g_clear_pointer (&ssc_sij->kernels, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_ssc_sij_parent_class)->dispose (object);
}

static void
_nc_xcor_ssc_sij_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_ssc_sij_parent_class)->finalize (object);
}

static void
nc_xcor_ssc_sij_class_init (NcXcorSSCSijClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_xcor_ssc_sij_set_property;
  object_class->get_property = &_nc_xcor_ssc_sij_get_property;
  object_class->constructed  = &_nc_xcor_ssc_sij_constructed;
  object_class->dispose      = &_nc_xcor_ssc_sij_dispose;
  object_class->finalize     = &_nc_xcor_ssc_sij_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_POWSPEC,
                                   g_param_spec_object ("powspec",
                                                        NULL,
                                                        "Linear matter power spectrum",
                                                        NCM_TYPE_POWSPEC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Z_EDGES,
                                   g_param_spec_object ("z-edges",
                                                        NULL,
                                                        "Redshift bin edges",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MASK_CL,
                                   g_param_spec_object ("mask-cl",
                                                        NULL,
                                                        "Angular power spectrum of the survey mask",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_AREA,
                                   g_param_spec_double ("area",
                                                        NULL,
                                                        "Survey area in square degrees for the f_sky rescaling, 0 to disable",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_METHOD,
                                   g_param_spec_enum ("method",
                                                      NULL,
                                                      "Quadrature method used for the angular power spectra",
                                                      NC_TYPE_XCOR_METHOD,
                                                      NC_XCOR_METHOD_KERNEL_EXACT,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BLOCK_SIZE,
                                   g_param_spec_uint ("block-size",
                                                      NULL,
                                                      "Multipole block size for the solver",
                                                      1, G_MAXUINT, NC_XCOR_SSC_SIJ_DEFAULT_BLOCK_SIZE,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance of the kernel spline and the outer k integral",
                                                        GSL_DBL_EPSILON, 1.0, NC_XCOR_SSC_SIJ_DEFAULT_RELTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SCALED_ABSTOL,
                                   g_param_spec_double ("scaled-abstol",
                                                        NULL,
                                                        "Absolute floor of the adaptive refinement of the U_i(k) spline",
                                                        0.0, 1.0, NC_XCOR_SSC_SIJ_DEFAULT_SCALED_ABSTOL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_xcor_ssc_sij_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec, the linear matter power spectrum
 * @z_edges: a #NcmVector of strictly increasing redshift bin edges
 *
 * Creates a new #NcXcorSSCSij over the `len (@z_edges) - 1` top-hat redshift
 * bins delimited by @z_edges. The footprint defaults to the full sky; use
 * nc_xcor_ssc_sij_set_mask_cl() to set a real one.
 *
 * Returns: (transfer full): a new #NcXcorSSCSij
 */
NcXcorSSCSij *
nc_xcor_ssc_sij_new (NcDistance *dist, NcmPowspec *ps, NcmVector *z_edges)
{
  return g_object_new (NC_TYPE_XCOR_SSC_SIJ,
                       "dist", dist,
                       "powspec", ps,
                       "z-edges", z_edges,
                       NULL);
}

/**
 * nc_xcor_ssc_sij_ref:
 * @ssc_sij: a #NcXcorSSCSij
 *
 * Increments the reference count of @ssc_sij.
 *
 * Returns: (transfer full): @ssc_sij
 */
NcXcorSSCSij *
nc_xcor_ssc_sij_ref (NcXcorSSCSij *ssc_sij)
{
  return g_object_ref (ssc_sij);
}

/**
 * nc_xcor_ssc_sij_free:
 * @ssc_sij: a #NcXcorSSCSij
 *
 * Decrements the reference count of @ssc_sij, freeing it if the count reaches 0.
 *
 */
void
nc_xcor_ssc_sij_free (NcXcorSSCSij *ssc_sij)
{
  g_object_unref (ssc_sij);
}

/**
 * nc_xcor_ssc_sij_clear:
 * @ssc_sij: a #NcXcorSSCSij
 *
 * If *@ssc_sij is not %NULL, decrements its reference count, freeing it if the
 * count reaches 0, and sets *@ssc_sij to %NULL.
 *
 */
void
nc_xcor_ssc_sij_clear (NcXcorSSCSij **ssc_sij)
{
  g_clear_object (ssc_sij);
}

/**
 * nc_xcor_ssc_sij_mask_cl_fullsky:
 *
 * Builds the mask spectrum of the full sky, $C^{\rm mask}_\ell = 4\pi
 * \delta_{\ell 0}$, for which $S_{ij}$ reduces to $C^{ij}_0 / 4\pi$. This is
 * the default footprint of a newly created #NcXcorSSCSij.
 *
 * Returns: (transfer full): a length one #NcmVector holding $4\pi$
 */
NcmVector *
nc_xcor_ssc_sij_mask_cl_fullsky (void)
{
  NcmVector *mask_cl = ncm_vector_new (1);

  ncm_vector_set (mask_cl, 0, 4.0 * ncm_c_pi ());

  return mask_cl;
}

/**
 * nc_xcor_ssc_sij_get_nbins:
 * @ssc_sij: a #NcXcorSSCSij
 *
 * Returns: the number of redshift bins
 */
guint
nc_xcor_ssc_sij_get_nbins (NcXcorSSCSij *ssc_sij)
{
  return ssc_sij->kernels->len;
}

/**
 * nc_xcor_ssc_sij_get_lmax:
 * @ssc_sij: a #NcXcorSSCSij
 *
 * Returns: the highest multipole entering the sum, one less than the length of
 * the mask spectrum
 */
guint
nc_xcor_ssc_sij_get_lmax (NcXcorSSCSij *ssc_sij)
{
  return ncm_vector_len (ssc_sij->mask_cl) - 1;
}

/**
 * nc_xcor_ssc_sij_get_fsky:
 * @ssc_sij: a #NcXcorSSCSij
 *
 * Computes the sky fraction of the footprint, $f_{\rm sky} = \sqrt{C^{\rm
 * mask}_0 / 4\pi}$.
 *
 * Returns: the sky fraction
 */
gdouble
nc_xcor_ssc_sij_get_fsky (NcXcorSSCSij *ssc_sij)
{
  if (ssc_sij->area > 0.0)
    return ssc_sij->area * gsl_pow_2 (ncm_c_pi () / 180.0) / (4.0 * ncm_c_pi ());

  return sqrt (ncm_vector_get (ssc_sij->mask_cl, 0) / (4.0 * ncm_c_pi ()));
}

/**
 * nc_xcor_ssc_sij_set_mask_cl:
 * @ssc_sij: a #NcXcorSSCSij
 * @mask_cl: (nullable): the mask angular power spectrum, or %NULL for the full sky
 *
 * Sets the angular power spectrum of the survey footprint, $C^{\rm
 * mask}_\ell$ for $\ell$ in $[0, \ell_{\rm max}]$, so @mask_cl has length
 * $\ell_{\rm max} + 1$ and its truncation sets how many multipoles are summed.
 * #NcmSphereMap computes it from a HEALPix map; it does not depend on
 * cosmology, so it is set once and reused for every cosmology afterwards.
 *
 * Passing %NULL restores the full-sky spectrum, see
 * nc_xcor_ssc_sij_mask_cl_fullsky().
 *
 */
void
nc_xcor_ssc_sij_set_mask_cl (NcXcorSSCSij *ssc_sij, NcmVector *mask_cl)
{
  if ((mask_cl != NULL) && (ncm_vector_len (mask_cl) < 1))
    g_error ("nc_xcor_ssc_sij_set_mask_cl: mask-cl must contain at least the monopole.");

  if ((mask_cl != NULL) && (ncm_vector_get (mask_cl, 0) <= 0.0))
    g_error ("nc_xcor_ssc_sij_set_mask_cl: the mask monopole must be positive, got %.17g.",
             ncm_vector_get (mask_cl, 0));

  ncm_vector_clear (&ssc_sij->mask_cl);

  if (mask_cl != NULL)
  {
    ssc_sij->mask_cl = ncm_vector_ref (mask_cl);

    /* A real footprint already carries its own sky fraction, so the crude
     * area rescaling must not be applied on top of it. */
    ssc_sij->area = 0.0;
  }
  else if (ssc_sij->constructed)
  {
    ssc_sij->mask_cl = nc_xcor_ssc_sij_mask_cl_fullsky ();
  }

  _nc_xcor_ssc_sij_invalidate (ssc_sij);
}

/**
 * nc_xcor_ssc_sij_peek_mask_cl:
 * @ssc_sij: a #NcXcorSSCSij
 *
 * Returns: (transfer none): the mask angular power spectrum in use
 */
NcmVector *
nc_xcor_ssc_sij_peek_mask_cl (NcXcorSSCSij *ssc_sij)
{
  return ssc_sij->mask_cl;
}

/**
 * nc_xcor_ssc_sij_set_area:
 * @ssc_sij: a #NcXcorSSCSij
 * @area: the survey area in square degrees, or 0.0 to disable
 *
 * Selects the crude finite-area estimator, $S_{ij} = S^{\rm fullsky}_{ij} /
 * f_{\rm sky}$ with $f_{\rm sky} = \Omega / 4\pi$: the full-sky matrix simply
 * rescaled by the sky fraction @area subtends. This is **not** a mask
 * deconvolution, and it knows nothing about the shape of the footprint; use
 * nc_xcor_ssc_sij_set_mask_cl() when the shape matters.
 *
 * It is kept because it is what an analysis quoting only a survey area uses,
 * and because it costs a single multipole rather than the whole $\ell$ range a
 * mask needs.
 *
 * Setting a positive @area discards any mask spectrum, the two being mutually
 * exclusive; passing 0.0 restores the plain full-sky matrix.
 *
 */
void
nc_xcor_ssc_sij_set_area (NcXcorSSCSij *ssc_sij, gdouble area)
{
  const gdouble sqd_fullsky = 4.0 * ncm_c_pi () * gsl_pow_2 (180.0 / ncm_c_pi ());

  if (area < 0.0)
    g_error ("nc_xcor_ssc_sij_set_area: area must not be negative, got %.17g.", area);

  if (area >= sqd_fullsky)
    g_error ("nc_xcor_ssc_sij_set_area: area must be smaller than the whole sky "
             "(%.1f square degrees), got %.17g.", sqd_fullsky, area);

  if ((area > 0.0) && ssc_sij->constructed)
  {
    ncm_vector_clear (&ssc_sij->mask_cl);
    ssc_sij->mask_cl = nc_xcor_ssc_sij_mask_cl_fullsky ();
  }

  ssc_sij->area = area;

  _nc_xcor_ssc_sij_invalidate (ssc_sij);
}

/**
 * nc_xcor_ssc_sij_get_area:
 * @ssc_sij: a #NcXcorSSCSij
 *
 * Returns: the survey area in square degrees used for the $f_{\rm sky}$
 * rescaling, or 0.0 when it is disabled
 */
gdouble
nc_xcor_ssc_sij_get_area (NcXcorSSCSij *ssc_sij)
{
  return ssc_sij->area;
}

/**
 * nc_xcor_ssc_sij_set_block_size:
 * @ssc_sij: a #NcXcorSSCSij
 * @block_size: the multipole block size
 *
 * Sets the multipole block size handed to nc_xcor_solver_plan_blocks().
 *
 */
void
nc_xcor_ssc_sij_set_block_size (NcXcorSSCSij *ssc_sij, guint block_size)
{
  g_assert_cmpuint (block_size, >, 0);

  ssc_sij->block_size = block_size;

  _nc_xcor_ssc_sij_invalidate (ssc_sij);
}

/**
 * nc_xcor_ssc_sij_get_block_size:
 * @ssc_sij: a #NcXcorSSCSij
 *
 * Returns: the multipole block size
 */
guint
nc_xcor_ssc_sij_get_block_size (NcXcorSSCSij *ssc_sij)
{
  return ssc_sij->block_size;
}

/**
 * nc_xcor_ssc_sij_set_reltol:
 * @ssc_sij: a #NcXcorSSCSij
 * @reltol: the relative tolerance
 *
 * Sets the relative tolerance of the kernel splines and of the outer $k$
 * integral. This is not the knob limiting the accuracy of the off-diagonal
 * $S_{ij}$, see nc_xcor_ssc_sij_set_scaled_abstol().
 *
 */
void
nc_xcor_ssc_sij_set_reltol (NcXcorSSCSij *ssc_sij, gdouble reltol)
{
  guint i;

  ssc_sij->reltol = reltol;

  for (i = 0; i < ssc_sij->kernels->len; i++)
    nc_xcor_kernel_set_reltol (g_ptr_array_index (ssc_sij->kernels, i), reltol);

  if (ssc_sij->xcor != NULL)
    nc_xcor_set_reltol (ssc_sij->xcor, reltol);

  _nc_xcor_ssc_sij_invalidate (ssc_sij);
}

/**
 * nc_xcor_ssc_sij_get_reltol:
 * @ssc_sij: a #NcXcorSSCSij
 *
 * Returns: the relative tolerance
 */
gdouble
nc_xcor_ssc_sij_get_reltol (NcXcorSSCSij *ssc_sij)
{
  return ssc_sij->reltol;
}

/**
 * nc_xcor_ssc_sij_set_scaled_abstol:
 * @ssc_sij: a #NcXcorSSCSij
 * @scaled_abstol: the absolute floor of the adaptive refinement
 *
 * Sets the absolute floor of the adaptive refinement building the $U_i(k)$
 * spline of every kernel. **This, not the relative tolerance, is what limits
 * the accuracy of the off-diagonal $S_{ij}$**, which are a small residual of a
 * large cancellation: tightening it from #NcXcorKernel's own `1.0e-4` default
 * to the `1.0e-6` used here moves $S_{06}$ by tens of percent for J-PAS-like
 * bins while barely moving $S_{00}$.
 *
 * `1.0e-6` is the end of that road, not a waypoint. The floor is a fraction of
 * the peak of $W_i(k)$ while the integrand is $k^2 W_i W_j$, so it enters
 * squared: `1.0e-6` is already `1.0e-12` there. Below it
 * nc_xcor_kernel_set_scaled_abstol() warns and the accuracy is not recoverable
 * at any cost -- see %NC_XCOR_KERNEL_MIN_USEFUL_SCALED_ABSTOL.
 *
 */
void
nc_xcor_ssc_sij_set_scaled_abstol (NcXcorSSCSij *ssc_sij, gdouble scaled_abstol)
{
  guint i;

  ssc_sij->scaled_abstol = scaled_abstol;

  for (i = 0; i < ssc_sij->kernels->len; i++)
    nc_xcor_kernel_set_scaled_abstol (g_ptr_array_index (ssc_sij->kernels, i), scaled_abstol);

  _nc_xcor_ssc_sij_invalidate (ssc_sij);
}

/**
 * nc_xcor_ssc_sij_get_scaled_abstol:
 * @ssc_sij: a #NcXcorSSCSij
 *
 * Returns: the absolute floor of the adaptive refinement
 */
gdouble
nc_xcor_ssc_sij_get_scaled_abstol (NcXcorSSCSij *ssc_sij)
{
  return ssc_sij->scaled_abstol;
}

/**
 * nc_xcor_ssc_sij_set_method:
 * @ssc_sij: a #NcXcorSSCSij
 * @method: a #NcXcorMethod
 *
 * Sets the quadrature method used for the angular power spectra. The default,
 * %NC_XCOR_METHOD_KERNEL_EXACT, needs no tolerance and cannot fail to
 * converge; the adaptive alternatives abort when they cannot reach their
 * target tolerance, which is fatal inside a Monte Carlo chain.
 *
 */
void
nc_xcor_ssc_sij_set_method (NcXcorSSCSij *ssc_sij, NcXcorMethod method)
{
  ssc_sij->method = method;

  if ((ssc_sij->xcor != NULL) && (nc_xcor_get_meth (ssc_sij->xcor) != method))
  {
    nc_xcor_clear (&ssc_sij->xcor);
    ssc_sij->xcor = nc_xcor_new (ssc_sij->dist, ssc_sij->ps, method);
    nc_xcor_set_reltol (ssc_sij->xcor, ssc_sij->reltol);
  }

  _nc_xcor_ssc_sij_invalidate (ssc_sij);
}

/**
 * nc_xcor_ssc_sij_get_method:
 * @ssc_sij: a #NcXcorSSCSij
 *
 * Returns: the quadrature method in use
 */
NcXcorMethod
nc_xcor_ssc_sij_get_method (NcXcorSSCSij *ssc_sij)
{
  return ssc_sij->method;
}

/* Registers every kernel once and requests the upper triangle of the pair
 * matrix over the whole ell range, then plans the blocks. The result is reused
 * across cosmologies: only nc_xcor_solver_solve() is cosmology dependent. */
static void
_nc_xcor_ssc_sij_ensure_solver (NcXcorSSCSij *ssc_sij)
{
  const guint nbins = ssc_sij->kernels->len;
  const guint lmax  = nc_xcor_ssc_sij_get_lmax (ssc_sij);
  NcmSBesselIntegrator *sbi;
  guint i, j;

  if (ssc_sij->solver != NULL)
    return;

  ssc_sij->solver = nc_xcor_solver_new ();

  for (i = 0; i < nbins; i++)
  {
    const guint kernel_id = nc_xcor_solver_register_kernel (ssc_sij->solver, g_ptr_array_index (ssc_sij->kernels, i));

    /* The kernels are distinct instances, so registration hands them ids in
     * order; the request loop below relies on it, as does the readout. */
    g_assert_cmpuint (kernel_id, ==, i);
  }

  for (i = 0; i < nbins; i++)
  {
    for (j = i; j < nbins; j++)
      nc_xcor_solver_request_cl (ssc_sij->solver, i, j, 0, lmax);
  }

  sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (0, ssc_sij->block_size));
  nc_xcor_solver_set_integrator (ssc_sij->solver, sbi);
  ncm_sbessel_integrator_free (sbi);

  nc_xcor_solver_plan_blocks (ssc_sij->solver, ssc_sij->block_size);
}

/**
 * nc_xcor_ssc_sij_prepare:
 * @ssc_sij: a #NcXcorSSCSij
 * @cosmo: a #NcHICosmo
 *
 * Computes the $S_{ij}$ matrix for @cosmo, unconditionally. Use
 * nc_xcor_ssc_sij_prepare_if_needed() to skip the work when @cosmo has not
 * changed since the last call.
 *
 */
void
nc_xcor_ssc_sij_prepare (NcXcorSSCSij *ssc_sij, NcHICosmo *cosmo)
{
  const guint nbins   = ssc_sij->kernels->len;
  const guint lmax    = nc_xcor_ssc_sij_get_lmax (ssc_sij);
  const gdouble norma = 4.0 * ncm_c_pi () * ncm_vector_get (ssc_sij->mask_cl, 0)
                        * ((ssc_sij->area > 0.0) ? nc_xcor_ssc_sij_get_fsky (ssc_sij) : 1.0);
  guint request = 0;
  guint i, j, l;

  _nc_xcor_ssc_sij_check_ready (ssc_sij);
  _nc_xcor_ssc_sij_ensure_solver (ssc_sij);

  nc_distance_prepare_if_needed (ssc_sij->dist, cosmo);
  ncm_powspec_prepare_if_needed (ssc_sij->ps, NCM_MODEL (cosmo));

  for (i = 0; i < nbins; i++)
    nc_xcor_kernel_prepare (g_ptr_array_index (ssc_sij->kernels, i), cosmo);

  nc_xcor_prepare (ssc_sij->xcor, cosmo);
  nc_xcor_solver_solve (ssc_sij->solver, ssc_sij->xcor, cosmo);

  /* S_ij = sum_l (2l+1) Cl_mask[l] Cl_ij[l] / (4 pi fsky)^2, and since
   * fsky^2 = Cl_mask[0] / (4 pi) the normalization is 4 pi Cl_mask[0]. */
  for (i = 0; i < nbins; i++)
  {
    for (j = i; j < nbins; j++)
    {
      NcmVector *cl = nc_xcor_solver_get_result (ssc_sij->solver, request);
      gdouble Sij   = 0.0;

      for (l = 0; l <= lmax; l++)
        Sij += (2.0 * l + 1.0) * ncm_vector_get (ssc_sij->mask_cl, l) * ncm_vector_get (cl, l);

      Sij /= norma;

      ncm_matrix_set (ssc_sij->sij, i, j, Sij);
      ncm_matrix_set (ssc_sij->sij, j, i, Sij);

      request++;
    }
  }

  ncm_model_ctrl_update (ssc_sij->ctrl, NCM_MODEL (cosmo));
  ssc_sij->prepared = TRUE;
}

/**
 * nc_xcor_ssc_sij_prepare_if_needed:
 * @ssc_sij: a #NcXcorSSCSij
 * @cosmo: a #NcHICosmo
 *
 * Computes the $S_{ij}$ matrix for @cosmo, unless it is already the one held
 * from a previous call. This is the entry point for a likelihood that
 * recomputes $S_{ij}$ per step: it is a no-op whenever @cosmo has not moved.
 *
 */
void
nc_xcor_ssc_sij_prepare_if_needed (NcXcorSSCSij *ssc_sij, NcHICosmo *cosmo)
{
  if (ncm_model_ctrl_update (ssc_sij->ctrl, NCM_MODEL (cosmo)) || !ssc_sij->prepared)
    nc_xcor_ssc_sij_prepare (ssc_sij, cosmo);
}

/**
 * nc_xcor_ssc_sij_peek_matrix:
 * @ssc_sij: a #NcXcorSSCSij
 *
 * Gets the $S_{ij}$ matrix computed by the last nc_xcor_ssc_sij_prepare(). The
 * matrix is owned by @ssc_sij and is overwritten by the next prepare.
 *
 * Returns: (transfer none): the $S_{ij}$ matrix, of size `nbins` by `nbins`
 */
NcmMatrix *
nc_xcor_ssc_sij_peek_matrix (NcXcorSSCSij *ssc_sij)
{
  return ssc_sij->sij;
}

/**
 * nc_xcor_ssc_sij_eval:
 * @ssc_sij: a #NcXcorSSCSij
 * @cosmo: a #NcHICosmo
 *
 * Prepares @ssc_sij for @cosmo if needed and returns a copy of the resulting
 * $S_{ij}$ matrix, which therefore survives the next prepare.
 *
 * Returns: (transfer full): a new #NcmMatrix holding $S_{ij}$
 */
NcmMatrix *
nc_xcor_ssc_sij_eval (NcXcorSSCSij *ssc_sij, NcHICosmo *cosmo)
{
  nc_xcor_ssc_sij_prepare_if_needed (ssc_sij, cosmo);

  return ncm_matrix_dup (ssc_sij->sij);
}

