/***************************************************************************
 *            nc_xcor_solver.c
 *
 *  Thu August 06 2026
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
 * NcXcorSolver:
 *
 * Batched register/request/solve interface for #NcXcorKernel angular power
 * spectra.
 *
 * This is the "companion object" from
 * dev-notes/xcor_ultralevin_batching_plan.md §5: it lets a caller declare
 * every kernel and every desired auto/cross Cl block up front
 * (nc_xcor_solver_register_kernel() / nc_xcor_solver_request_cl()), so a
 * future nc_xcor_solver_solve() can pick an ℓ-block tiling and build each
 * distinct kernel's k-space closure exactly once per block, shared across
 * every request that needs it in that block -- instead of rebuilding it
 * once per pair, as #NcXcor's single-pair nc_xcor_compute() does today.
 *
 * This phase (registration) does no computation at all: it only records
 * which kernels and which (kernel pair, ℓ-range) blocks are wanted.
 * nc_xcor_solver_register_kernel() is idempotent by pointer identity --
 * registering the same #NcXcorKernel instance twice returns the same id.
 *
 * `NcXcor` itself stays untouched and close to stateless per-call; this
 * object exists precisely so bolting registration/request state onto
 * `NcXcor` doesn't change its contract for existing callers of
 * nc_xcor_compute().
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/xcor/nc_xcor_solver.h"
#include "nc/xcor/nc_xcor_priv.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/specfunc/ncm_sbessel_integrator_levin.h"

typedef struct _NcXcorSolverRequest
{
  guint kernel_id_1;
  guint kernel_id_2;
  guint lmin;
  guint lmax;
} NcXcorSolverRequest;

typedef struct _NcXcorSolverBlock
{
  guint lmin;
  guint lmax;
} NcXcorSolverBlock;

struct _NcXcorSolver
{
  /*< private > */
  GObject parent_instance;
  GPtrArray *kernels; /* element-type NcXcorKernel*, owned refs */
  GArray *requests;   /* element-type NcXcorSolverRequest */
  GArray *blocks;     /* element-type NcXcorSolverBlock, from plan_blocks() */
  GPtrArray *results; /* element-type NcmVector*, owned refs, from solve() */

  /* One integrator per ell-block, indexed by block, NULL where the block needs
   * none. Each is pinned to its block's multipole range for its whole life, so
   * the factorization a #NcmSBesselIntegratorLevin builds for that range stays
   * valid and is reused on every later solve(). Cleared by plan_blocks(). */
  NcmSBesselIntegrator *integrator_proto; /* template the per-block ones are copied from, owned */
  GPtrArray *block_integrators;           /* element-type NcmSBesselIntegrator*, owned refs, may hold NULLs */
};

G_DEFINE_TYPE (NcXcorSolver, nc_xcor_solver, G_TYPE_OBJECT)

/* The per-block array is sized to the block count and holds NULL for blocks
 * that need no integrator, so its free function must accept NULL. */
static void
_nc_xcor_solver_integrator_free0 (gpointer sbi)
{
  if (sbi != NULL)
    ncm_sbessel_integrator_free (NCM_SBESSEL_INTEGRATOR (sbi));
}

static void
nc_xcor_solver_init (NcXcorSolver *solver)
{
  solver->kernels  = g_ptr_array_new_with_free_func ((GDestroyNotify) nc_xcor_kernel_free);
  solver->requests = g_array_new (FALSE, FALSE, sizeof (NcXcorSolverRequest));
  solver->blocks   = g_array_new (FALSE, FALSE, sizeof (NcXcorSolverBlock));
  solver->results  = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_vector_free);

  solver->integrator_proto  = NULL;
  solver->block_integrators = g_ptr_array_new_with_free_func (_nc_xcor_solver_integrator_free0);
}

static void
_nc_xcor_solver_dispose (GObject *object)
{
  NcXcorSolver *solver = NC_XCOR_SOLVER (object);

  g_clear_pointer (&solver->kernels, g_ptr_array_unref);
  g_clear_pointer (&solver->requests, g_array_unref);
  g_clear_pointer (&solver->blocks, g_array_unref);
  g_clear_pointer (&solver->results, g_ptr_array_unref);
  g_clear_pointer (&solver->block_integrators, g_ptr_array_unref);
  ncm_sbessel_integrator_clear (&solver->integrator_proto);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_solver_parent_class)->dispose (object);
}

static void
_nc_xcor_solver_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_solver_parent_class)->finalize (object);
}

static void
nc_xcor_solver_class_init (NcXcorSolverClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->dispose  = &_nc_xcor_solver_dispose;
  object_class->finalize = &_nc_xcor_solver_finalize;
}

/**
 * nc_xcor_solver_new:
 *
 * Creates a new, empty #NcXcorSolver.
 *
 * Returns: (transfer full): a new #NcXcorSolver
 */
NcXcorSolver *
nc_xcor_solver_new (void)
{
  return g_object_new (NC_TYPE_XCOR_SOLVER, NULL);
}

/**
 * nc_xcor_solver_ref:
 * @solver: a #NcXcorSolver
 *
 * Increments the reference count of @solver.
 *
 * Returns: (transfer full): @solver
 */
NcXcorSolver *
nc_xcor_solver_ref (NcXcorSolver *solver)
{
  return g_object_ref (solver);
}

/**
 * nc_xcor_solver_free:
 * @solver: a #NcXcorSolver
 *
 * Decrements the reference count of @solver, freeing it if the count
 * reaches 0.
 *
 */
void
nc_xcor_solver_free (NcXcorSolver *solver)
{
  g_object_unref (solver);
}

/**
 * nc_xcor_solver_clear:
 * @solver: a #NcXcorSolver
 *
 * If *@solver is not %NULL, decrements its reference count, freeing it if
 * the count reaches 0, and sets *@solver to %NULL.
 *
 */
void
nc_xcor_solver_clear (NcXcorSolver **solver)
{
  g_clear_object (solver);
}

/**
 * nc_xcor_solver_register_kernel:
 * @solver: a #NcXcorSolver
 * @xclk: a #NcXcorKernel
 *
 * Registers @xclk with @solver, returning an id to use in
 * nc_xcor_solver_request_cl(). Idempotent by pointer identity: registering
 * the same @xclk instance again returns the same id without adding a
 * second entry.
 *
 * Returns: the kernel id, stable for the lifetime of @solver
 */
guint
nc_xcor_solver_register_kernel (NcXcorSolver *solver, NcXcorKernel *xclk)
{
  guint i;

  g_return_val_if_fail (NC_IS_XCOR_SOLVER (solver), 0);
  g_return_val_if_fail (NC_IS_XCOR_KERNEL (xclk), 0);

  for (i = 0; i < solver->kernels->len; i++)
  {
    if (g_ptr_array_index (solver->kernels, i) == xclk)
      return i;
  }

  g_ptr_array_add (solver->kernels, nc_xcor_kernel_ref (xclk));

  return solver->kernels->len - 1;
}

/**
 * nc_xcor_solver_request_cl:
 * @solver: a #NcXcorSolver
 * @kernel_id_1: id of the first kernel, from nc_xcor_solver_register_kernel()
 * @kernel_id_2: id of the second kernel (equal to @kernel_id_1 for an
 * auto-correlation), from nc_xcor_solver_register_kernel()
 * @lmin: minimum multipole
 * @lmax: maximum multipole
 *
 * Declares one desired output block: the auto- or cross-Cl of the kernels
 * identified by @kernel_id_1 and @kernel_id_2, over ℓ in [@lmin, @lmax].
 * Call once per entry in the data vector / covariance block ultimately
 * needed. Registration-only: no computation happens here.
 *
 */
void
nc_xcor_solver_request_cl (NcXcorSolver *solver, guint kernel_id_1, guint kernel_id_2, guint lmin, guint lmax)
{
  NcXcorSolverRequest req;

  g_return_if_fail (NC_IS_XCOR_SOLVER (solver));
  g_return_if_fail (kernel_id_1 < solver->kernels->len);
  g_return_if_fail (kernel_id_2 < solver->kernels->len);

  if (lmax < lmin)
    g_error ("nc_xcor_solver_request_cl: lmax (%u) < lmin (%u)", lmax, lmin);

  req.kernel_id_1 = kernel_id_1;
  req.kernel_id_2 = kernel_id_2;
  req.lmin        = lmin;
  req.lmax        = lmax;

  g_array_append_val (solver->requests, req);
}

/**
 * nc_xcor_solver_get_n_kernels:
 * @solver: a #NcXcorSolver
 *
 * Returns: the number of distinct kernels registered so far
 */
guint
nc_xcor_solver_get_n_kernels (NcXcorSolver *solver)
{
  g_return_val_if_fail (NC_IS_XCOR_SOLVER (solver), 0);

  return solver->kernels->len;
}

/**
 * nc_xcor_solver_get_n_requests:
 * @solver: a #NcXcorSolver
 *
 * Returns: the number of Cl blocks requested so far
 */
guint
nc_xcor_solver_get_n_requests (NcXcorSolver *solver)
{
  g_return_val_if_fail (NC_IS_XCOR_SOLVER (solver), 0);

  return solver->requests->len;
}

/**
 * nc_xcor_solver_peek_kernel:
 * @solver: a #NcXcorSolver
 * @kernel_id: a kernel id, from nc_xcor_solver_register_kernel()
 *
 * Returns: (transfer none): the kernel registered under @kernel_id
 */
NcXcorKernel *
nc_xcor_solver_peek_kernel (NcXcorSolver *solver, guint kernel_id)
{
  g_return_val_if_fail (NC_IS_XCOR_SOLVER (solver), NULL);
  g_return_val_if_fail (kernel_id < solver->kernels->len, NULL);

  return g_ptr_array_index (solver->kernels, kernel_id);
}

/**
 * nc_xcor_solver_get_request:
 * @solver: a #NcXcorSolver
 * @request_index: index of the request, in [0, nc_xcor_solver_get_n_requests())
 * @kernel_id_1: (out): first kernel id of the request
 * @kernel_id_2: (out): second kernel id of the request
 * @lmin: (out): minimum multipole of the request
 * @lmax: (out): maximum multipole of the request
 *
 * Retrieves the details of a previously-made request, in registration
 * order.
 *
 */
void
nc_xcor_solver_get_request (NcXcorSolver *solver, guint request_index, guint *kernel_id_1, guint *kernel_id_2, guint *lmin, guint *lmax)
{
  NcXcorSolverRequest *req;

  g_return_if_fail (NC_IS_XCOR_SOLVER (solver));
  g_return_if_fail (request_index < solver->requests->len);

  req = &g_array_index (solver->requests, NcXcorSolverRequest, request_index);

  *kernel_id_1 = req->kernel_id_1;
  *kernel_id_2 = req->kernel_id_2;
  *lmin        = req->lmin;
  *lmax        = req->lmax;
}

/**
 * nc_xcor_solver_clear_requests:
 * @solver: a #NcXcorSolver
 *
 * Removes all requests, keeping registered kernels intact. Useful to reuse
 * the same kernel registration across several rounds of requests.
 *
 */
void
nc_xcor_solver_clear_requests (NcXcorSolver *solver)
{
  g_return_if_fail (NC_IS_XCOR_SOLVER (solver));

  g_array_set_size (solver->requests, 0);
}

/**
 * nc_xcor_solver_set_integrator:
 * @solver: a #NcXcorSolver
 * @sbi: (nullable): a #NcmSBesselIntegrator to copy per block, or %NULL
 *
 * Sets the integrator @solver copies once per ell-block. Every block gets its
 * own copy of @sbi, pinned to that block's multipole range, so all registered
 * kernels share one integrator per block instead of each carrying its own.
 *
 * Passing %NULL restores the default, which is to copy the integrator of the
 * first registered kernel that has one. Either way the copies, not @sbi
 * itself, are used, and the kernels' own `integrator` properties are ignored
 * during nc_xcor_solver_solve() -- they still apply to nc_xcor_compute().
 *
 * Discards any existing per-block integrators, so the next
 * nc_xcor_solver_solve() rebuilds them.
 *
 */
void
nc_xcor_solver_set_integrator (NcXcorSolver *solver, NcmSBesselIntegrator *sbi)
{
  g_return_if_fail (NC_IS_XCOR_SOLVER (solver));
  g_return_if_fail ((sbi == NULL) || NCM_IS_SBESSEL_INTEGRATOR (sbi));

  ncm_sbessel_integrator_clear (&solver->integrator_proto);

  if (sbi != NULL)
    solver->integrator_proto = ncm_sbessel_integrator_ref (sbi);

  g_ptr_array_set_size (solver->block_integrators, 0);
}

/*
 * A block needs an integrator when any registered kernel is evaluated in
 * non-Limber mode there; this mirrors the branch in
 * nc_xcor_kernel_get_eval_vectorized_full(). Blocks are planned never to
 * straddle an l_limber threshold, so the block's lmin decides for the whole
 * block.
 */
static gboolean
_nc_xcor_solver_block_needs_integrator (NcXcorSolver *solver, NcXcorSolverBlock *block)
{
  guint i;

  for (i = 0; i < solver->kernels->len; i++)
  {
    NcXcorKernel *xclk    = g_ptr_array_index (solver->kernels, i);
    const gint l_limber   = nc_xcor_kernel_get_l_limber (xclk);
    const gboolean limber = (l_limber == 0) || ((l_limber > 0) && ((gint) block->lmin >= l_limber));

    if (!limber)
      return TRUE;
  }

  return FALSE;
}

/*
 * Fills in any per-block integrator not built yet, each a copy of the template
 * pinned to its block's range. Called serially before the parallel block loop:
 * building one is expensive (knots, the j_ell cache and the panel operators)
 * and happens once per block for the lifetime of @solver, not once per solve().
 */
static void
_nc_xcor_solver_prepare_block_integrators (NcXcorSolver *solver)
{
  const guint n_blocks        = solver->blocks->len;
  NcmSBesselIntegrator *proto = solver->integrator_proto;
  NcmSerialize *ser           = NULL;
  guint b;

  if (proto == NULL)
  {
    guint i;

    for (i = 0; (i < solver->kernels->len) && (proto == NULL); i++)
      proto = nc_xcor_kernel_peek_integrator (g_ptr_array_index (solver->kernels, i));
  }

  if (solver->block_integrators->len != n_blocks)
    g_ptr_array_set_size (solver->block_integrators, n_blocks);

  for (b = 0; b < n_blocks; b++)
  {
    NcXcorSolverBlock *block = &g_array_index (solver->blocks, NcXcorSolverBlock, b);
    NcmSBesselIntegrator *sbi;

    if (g_ptr_array_index (solver->block_integrators, b) != NULL)
      continue;

    if (!_nc_xcor_solver_block_needs_integrator (solver, block))
      continue;

    if (proto == NULL)
      g_error ("_nc_xcor_solver_prepare_block_integrators: block [%u, %u] needs an integrator "
               "but no registered kernel has one and none was set with "
               "nc_xcor_solver_set_integrator().", block->lmin, block->lmax);

    if (ser == NULL)
      ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);

    sbi = NCM_SBESSEL_INTEGRATOR (ncm_serialize_dup_obj (ser, G_OBJECT (proto)));
    ncm_serialize_reset (ser, TRUE);

    ncm_sbessel_integrator_set_ell_range (sbi, block->lmin, block->lmax);
    g_ptr_array_index (solver->block_integrators, b) = sbi;
  }

  ncm_serialize_clear (&ser);
}

/**
 * nc_xcor_solver_peek_block_integrator:
 * @solver: a #NcXcorSolver
 * @block_index: index of the block, in [0, nc_xcor_solver_get_n_blocks())
 *
 * Returns: (transfer none) (nullable): the integrator pinned to that block, or
 * %NULL if the block needs none or nc_xcor_solver_solve() has not run yet
 */
NcmSBesselIntegrator *
nc_xcor_solver_peek_block_integrator (NcXcorSolver *solver, guint block_index)
{
  g_return_val_if_fail (NC_IS_XCOR_SOLVER (solver), NULL);

  if (block_index >= solver->block_integrators->len)
    return NULL;

  return g_ptr_array_index (solver->block_integrators, block_index);
}

static gint
_nc_xcor_solver_guint_cmp (gconstpointer a, gconstpointer b)
{
  guint ua = *(const guint *) a;
  guint ub = *(const guint *) b;

  return (ua > ub) - (ua < ub);
}

/**
 * nc_xcor_solver_plan_blocks:
 * @solver: a #NcXcorSolver
 * @default_block_size: preferred number of multipoles per block
 *
 * Tiles the union of all requested ℓ-ranges into contiguous blocks, each
 * covering at most @default_block_size multipoles, further capped by
 * #NC_XCOR_KERNEL_MAX_ELL_BLOCK and by the smallest
 * ncm_sbessel_integrator_levin_get_ell_cache_max() among all registered
 * kernels using a #NcmSBesselIntegratorLevin. A block also never straddles
 * a registered kernel's l_limber threshold, since a kernel evaluated in
 * Limber mode below l_limber and non-Limber mode above it needs the switch
 * to fall on a block boundary.
 *
 * Replaces any blocks from a previous call. Requires at least one request.
 *
 */
void
nc_xcor_solver_plan_blocks (NcXcorSolver *solver, guint default_block_size)
{
  guint max_block_size;
  guint overall_lmin = G_MAXUINT;
  guint overall_lmax = 0;
  GArray *boundaries;
  guint next_boundary_idx = 0;
  guint cursor;
  guint i;

  g_return_if_fail (NC_IS_XCOR_SOLVER (solver));
  g_return_if_fail (default_block_size > 0);

  if (solver->requests->len == 0)
    g_error ("nc_xcor_solver_plan_blocks: no requests, nothing to plan");

  for (i = 0; i < solver->requests->len; i++)
  {
    NcXcorSolverRequest *req = &g_array_index (solver->requests, NcXcorSolverRequest, i);

    overall_lmin = MIN (overall_lmin, req->lmin);
    overall_lmax = MAX (overall_lmax, req->lmax);
  }

  max_block_size = MIN (default_block_size, (guint) NC_XCOR_KERNEL_MAX_ELL_BLOCK);

  for (i = 0; i < solver->kernels->len; i++)
  {
    NcXcorKernel *xclk               = g_ptr_array_index (solver->kernels, i);
    NcmSBesselIntegrator *integrator = nc_xcor_kernel_peek_integrator (xclk);

    if ((integrator != NULL) && NCM_IS_SBESSEL_INTEGRATOR_LEVIN (integrator))
    {
      guint ell_cache_max = ncm_sbessel_integrator_levin_get_ell_cache_max (NCM_SBESSEL_INTEGRATOR_LEVIN (integrator));

      max_block_size = MIN (max_block_size, ell_cache_max);
    }
  }

  /* Every registered kernel's l_limber, if it falls strictly inside the
   * overall range, forces a block boundary there: below it the kernel
   * uses one evaluation mode, at/above it another. */
  boundaries = g_array_new (FALSE, FALSE, sizeof (guint));

  for (i = 0; i < solver->kernels->len; i++)
  {
    NcXcorKernel *xclk = g_ptr_array_index (solver->kernels, i);
    gint l_limber      = nc_xcor_kernel_get_l_limber (xclk);

    if ((l_limber > 0) && ((guint) l_limber > overall_lmin) && ((guint) l_limber <= overall_lmax))
    {
      guint boundary = (guint) l_limber;

      g_array_append_val (boundaries, boundary);
    }
  }

  g_array_sort (boundaries, _nc_xcor_solver_guint_cmp);

  g_array_set_size (solver->blocks, 0);
  g_ptr_array_set_size (solver->block_integrators, 0);

  cursor = overall_lmin;

  while (cursor <= overall_lmax)
  {
    guint block_end = MIN (cursor + max_block_size - 1, overall_lmax);

    while ((next_boundary_idx < boundaries->len) &&
           (g_array_index (boundaries, guint, next_boundary_idx) <= cursor))
      next_boundary_idx++;

    if ((next_boundary_idx < boundaries->len) &&
        (g_array_index (boundaries, guint, next_boundary_idx) - 1 < block_end))
      block_end = g_array_index (boundaries, guint, next_boundary_idx) - 1;

    {
      NcXcorSolverBlock block = { cursor, block_end };

      g_array_append_val (solver->blocks, block);
    }

    cursor = block_end + 1;
  }

  g_array_unref (boundaries);
}

/**
 * nc_xcor_solver_get_n_blocks:
 * @solver: a #NcXcorSolver
 *
 * Returns: the number of blocks from the last nc_xcor_solver_plan_blocks() call
 */
guint
nc_xcor_solver_get_n_blocks (NcXcorSolver *solver)
{
  g_return_val_if_fail (NC_IS_XCOR_SOLVER (solver), 0);

  return solver->blocks->len;
}

/**
 * nc_xcor_solver_get_block:
 * @solver: a #NcXcorSolver
 * @block_index: index of the block, in [0, nc_xcor_solver_get_n_blocks())
 * @lmin: (out): minimum multipole of the block
 * @lmax: (out): maximum multipole of the block
 *
 * Retrieves the ℓ-range of a block planned by nc_xcor_solver_plan_blocks(),
 * in ascending order and non-overlapping.
 *
 */
void
nc_xcor_solver_get_block (NcXcorSolver *solver, guint block_index, guint *lmin, guint *lmax)
{
  NcXcorSolverBlock *block;

  g_return_if_fail (NC_IS_XCOR_SOLVER (solver));
  g_return_if_fail (block_index < solver->blocks->len);

  block = &g_array_index (solver->blocks, NcXcorSolverBlock, block_index);

  *lmin = block->lmin;
  *lmax = block->lmax;
}

/*
 * Computes and stitches in one request's contribution from one block,
 * using an already-built (or freshly cached) pair of integrands. Shared by
 * both the serial and block-parallel KERNEL_CUBATURE loops below; @results
 * is written at request-relative offsets that never overlap between two
 * different blocks (blocks tile disjoint ℓ-ranges), so concurrent calls
 * from different blocks/threads writing into the same @results array are
 * safe as long as each call's own (@kernels, @integrands) are private to
 * its thread.
 */
static void
_nc_xcor_solver_solve_block_request (NcXcor *xc, GPtrArray *kernels, GHashTable *integrands, NcXcorSolverRequest *req, guint result_index, NcXcorSolverBlock *block, NcHICosmo *cosmo, GPtrArray *results, NcmSBesselIntegrator *sbi)
{
  NcXcorKernel *k1, *k2;
  gboolean isauto;
  NcXcorKernelIntegrand *xclki1, *xclki2;
  NcmVector *vp, *block_vp;
  guint lo, hi, i;

  if ((req->lmax < block->lmin) || (req->lmin > block->lmax))
    return;

  k1     = g_ptr_array_index (kernels, req->kernel_id_1);
  k2     = g_ptr_array_index (kernels, req->kernel_id_2);
  isauto = (k1 == k2);

  _nc_xcor_check_kernel_tolerance (xc, k1);

  if (!isauto)
    _nc_xcor_check_kernel_tolerance (xc, k2);

  /* No blanket redshift-overlap short-circuit here: this path is
   * kernel-space only, where the two kernels couple through the outer k
   * integral rather than a shared z integral, so kernels over disjoint
   * redshift ranges still have a non-zero cross spectrum. Kernels that are
   * themselves in the Limber tier are the exception -- see nc_xcor_compute().
   *
   * Testing the block's lmin alone settles the whole block: plan_blocks()
   * forces a block boundary at every registered kernel's l_limber, so no block
   * straddles the threshold and every multipole in it is on the same side. */
  {
    guint l_zero = 0;

    if (_nc_xcor_kernels_limber_disjoint (k1, k2, isauto, &l_zero) && (l_zero <= block->lmin))
      return;  /* zero contribution, result vector already zeroed */
  }

  xclki1 = g_hash_table_lookup (integrands, k1);

  if (xclki1 == NULL)
  {
    xclki1 = nc_xcor_kernel_get_eval_vectorized_full (k1, cosmo, block->lmin, block->lmax, sbi, nc_xcor_get_closure_type (xc));
    g_hash_table_insert (integrands, k1, xclki1);
  }

  if (isauto)
  {
    xclki2 = NULL;
  }
  else
  {
    xclki2 = g_hash_table_lookup (integrands, k2);

    if (xclki2 == NULL)
    {
      xclki2 = nc_xcor_kernel_get_eval_vectorized_full (k2, cosmo, block->lmin, block->lmax, sbi, nc_xcor_get_closure_type (xc));
      g_hash_table_insert (integrands, k2, xclki2);
    }
  }

  block_vp = ncm_vector_new (block->lmax - block->lmin + 1);

  /* The block methods share everything except the outer quadrature: same
   * per-kernel closures, same per-block caching. Which one runs is the
   * table's answer, not an if here. */
  {
    const NcXcorKQuad *kquad = _nc_xcor_kquad_for_method (nc_xcor_get_meth (xc));

    kquad->block (xc, xclki1, xclki2 != NULL ? xclki2 : xclki1, block->lmin, block->lmax, isauto, block_vp, NULL);
  }

  vp = g_ptr_array_index (results, result_index);
  lo = MAX (req->lmin, block->lmin);
  hi = MIN (req->lmax, block->lmax);

  for (i = lo; i <= hi; i++)
    ncm_vector_set (vp, i - req->lmin, ncm_vector_get (block_vp, i - block->lmin));

  ncm_vector_free (block_vp);
}

/**
 * nc_xcor_solver_solve:
 * @solver: a #NcXcorSolver
 * @xc: a #NcXcor
 * @cosmo: a #NcHICosmo
 *
 * Computes every requested Cl block. Requires nc_xcor_solver_plan_blocks()
 * to have been called first. Replaces any results from a previous
 * nc_xcor_solver_solve() call.
 *
 * When @xc's method has a block quadrature -- %NC_XCOR_METHOD_KERNEL_CUBATURE,
 * %NC_XCOR_METHOD_KERNEL_EXACT or %NC_XCOR_METHOD_KERNEL_GSL_BLOCK -- each distinct
 * kernel's k-space closure (nc_xcor_kernel_get_eval_vectorized()) is built
 * once per ℓ-block and shared across every request needing it in that
 * block, instead of rebuilding it once per pair the way nc_xcor_compute()
 * does -- see plan doc dev-notes/xcor_ultralevin_batching_plan.md §5-§6.
 * Different ℓ-blocks share no reusable operator state with each other (a
 * fresh #NcmSBesselIntegratorLevin ℓ-range rebuild is required regardless),
 * so blocks are additionally processed in parallel across an OpenMP team
 * (plan doc §6.3). The registered kernels are *shared* by the team and are
 * prepared against @cosmo before it starts, so they are only ever read inside
 * it; what is private to each thread is the per-block integrator and the
 * integrand cache built from it. No thread touches another's integrand or
 * integrator state, so no locking is needed inside the per-block work itself.
 *
 * They share that whole path and differ only in the outer quadrature --
 * adaptive cubature, exact Gauss-Legendre over the common refinement of each
 * pair's knot sets, or qagp broken on those same knots -- so none of them
 * needs a solve loop of its own.
 *
 * Every other method (%NC_XCOR_METHOD_LIMBER_Z_GSL,
 * %NC_XCOR_METHOD_LIMBER_Z_CUBATURE, %NC_XCOR_METHOD_KERNEL_GSL) has no
 * block-shared closure to reuse -- tier 1 stays untouched by design (plan
 * doc §8), and %NC_XCOR_METHOD_KERNEL_GSL fits its closure per ℓ regardless
 * of pairing -- so those requests are computed directly with
 * nc_xcor_compute(), one call per request, serially, for correctness with no
 * reuse. Keeping %NC_XCOR_METHOD_KERNEL_GSL on that path is also what makes
 * it answer the same here as through nc_xcor_compute();
 * %NC_XCOR_METHOD_KERNEL_GSL_BLOCK is the one to use for its quadrature with
 * the sharing.
 *
 * Results are retrieved with nc_xcor_solver_get_result().
 *
 */
void
nc_xcor_solver_solve (NcXcorSolver *solver, NcXcor *xc, NcHICosmo *cosmo)
{
  guint n_requests;
  guint r;

  g_return_if_fail (NC_IS_XCOR_SOLVER (solver));
  g_return_if_fail (NC_IS_XCOR (xc));
  g_return_if_fail (NC_IS_HICOSMO (cosmo));

  if (solver->blocks->len == 0)
    g_error ("nc_xcor_solver_solve: no blocks planned, call nc_xcor_solver_plan_blocks() first");

  nc_xcor_prepare (xc, cosmo);

  n_requests = solver->requests->len;
  g_ptr_array_set_size (solver->results, 0);

  for (r = 0; r < n_requests; r++)
  {
    NcXcorSolverRequest *req = &g_array_index (solver->requests, NcXcorSolverRequest, r);
    NcmVector *vp            = ncm_vector_new (req->lmax - req->lmin + 1);

    ncm_vector_set_zero (vp);
    g_ptr_array_add (solver->results, vp);
  }

  /* A method with a block quadrature is one this solver can share closures
   * for; everything else goes through nc_xcor_compute(), one call per
   * request. */
  if (_nc_xcor_kquad_for_method (nc_xcor_get_meth (xc)) == NULL)
  {
    for (r = 0; r < n_requests; r++)
    {
      NcXcorSolverRequest *req = &g_array_index (solver->requests, NcXcorSolverRequest, r);
      NcXcorKernel *k1         = g_ptr_array_index (solver->kernels, req->kernel_id_1);
      NcXcorKernel *k2         = g_ptr_array_index (solver->kernels, req->kernel_id_2);
      NcmVector *vp            = g_ptr_array_index (solver->results, r);

      nc_xcor_compute (xc, k1, (k1 == k2) ? NULL : k2, cosmo, req->lmin, req->lmax, vp);
    }

    return;
  }

  {
    const guint n_blocks  = solver->blocks->len;
    GArray *requests      = solver->requests;
    GArray *blocks        = solver->blocks;
    GPtrArray *kernels    = solver->kernels;
    GPtrArray *results    = solver->results;
    GPtrArray *block_sbis = solver->block_integrators;
    guint k;

    /* Kernels are shared by every thread from here on, so prepare them, and
     * build the per-block integrators, before the parallel region. */
    for (k = 0; k < kernels->len; k++)
      nc_xcor_kernel_prepare (g_ptr_array_index (kernels, k), cosmo);

    _nc_xcor_solver_prepare_block_integrators (solver);

    #pragma omp parallel shared (xc, cosmo, n_blocks, requests, blocks, kernels, results, block_sbis)
    {
      GHashTable *integrands = g_hash_table_new_full (NULL, NULL, NULL, (GDestroyNotify) nc_xcor_kernel_integrand_unref);
      guint b;

      #pragma omp for schedule (dynamic)

      for (b = 0; b < n_blocks; b++)
      {
        NcXcorSolverBlock *block  = &g_array_index (blocks, NcXcorSolverBlock, b);
        NcmSBesselIntegrator *sbi = g_ptr_array_index (block_sbis, b);
        guint r_local;

        g_hash_table_remove_all (integrands);

        for (r_local = 0; r_local < requests->len; r_local++)
        {
          NcXcorSolverRequest *req = &g_array_index (requests, NcXcorSolverRequest, r_local);

          _nc_xcor_solver_solve_block_request (xc, kernels, integrands, req, r_local, block, cosmo, results, sbi);
        }
      }

      g_hash_table_unref (integrands);
    }
  }
}

/**
 * nc_xcor_solver_get_result:
 * @solver: a #NcXcorSolver
 * @request_index: index of the request, in [0, nc_xcor_solver_get_n_requests())
 *
 * Returns: (transfer none): the Cl vector computed by the last
 * nc_xcor_solver_solve() call for the request at @request_index, of length
 * (lmax - lmin + 1)
 */
NcmVector *
nc_xcor_solver_get_result (NcXcorSolver *solver, guint request_index)
{
  g_return_val_if_fail (NC_IS_XCOR_SOLVER (solver), NULL);
  g_return_val_if_fail (request_index < solver->results->len, NULL);

  return g_ptr_array_index (solver->results, request_index);
}

