/***************************************************************************
 *            nc_xcor_priv.h
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

#ifndef _NC_XCOR_PRIV_H_
#define _NC_XCOR_PRIV_H_

#include <numcosmo/nc/xcor/nc_xcor.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel.h>

G_BEGIN_DECLS

/*
 * #NcXcor's instance struct and the argument block its integrands carry. Not
 * public API: shared between nc_xcor.c and the two internal translation units
 * it was split into, nc_xcor_limber_z.c and nc_xcor_kquad.c.
 */
struct _NcXcor
{
  /*< private > */
  GObject parent_instance;
  NcDistance *dist;
  NcmPowspec *ps;
  gdouble RH;
  NcXcorMethod meth;
  NcXcorKernelClosure closure_type;
  gdouble reltol;
  guint ell_batch_size;
};

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

/*
 * The redshift-space Limber tier, nc_xcor_limber_z.c. Both take an explicit
 * [@zmin, @zmax] because nc_xcor_compute() narrows it to the kernels' common
 * support before calling.
 */
void _nc_xcor_limber_z_gsl (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gdouble zmin, gdouble zmax, gboolean isauto, NcmVector *vp);
void _nc_xcor_limber_z_cubature (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gdouble zmin, gdouble zmax, gboolean isauto, NcmVector *vp);

/*
 * The kernel-space tier, nc_xcor_kquad.c: one entry point per method, each
 * building its own closures per multipole batch. Only the exact method reports
 * an error estimate, and @vp_err is nullable there.
 */
void _nc_xcor_kernel_gsl (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gboolean isauto, NcmVector *vp);
void _nc_xcor_kernel_cubature (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gboolean isauto, NcmVector *vp);
void _nc_xcor_kernel_exact (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gboolean isauto, NcmVector *vp, NcmVector *vp_err);

/*
 * QUADPACK's status is a statement about certification, not about the answer,
 * so both GSL methods judge it against the error they achieved. Defined in
 * nc_xcor.c because both tiers use it.
 */
void _nc_xcor_check_qag_status (const gchar *where, gint ret, gdouble reltol, gdouble result, gdouble err);

/*
 * Not public API: shared only between nc_xcor.c and nc_xcor_solver.c, so
 * NcXcorSolver can reuse nc_xcor.c's KERNEL_CUBATURE outer-integral
 * machinery with integrand(s) it built and cached itself, instead of
 * letting nc_xcor_compute() rebuild them once per pair.
 */
void _nc_xcor_kernel_integrate_block_cubature (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcmVector *vp);

/*
 * Same for %NC_XCOR_METHOD_KERNEL_EXACT: exact 5-node Gauss-Legendre over the
 * common refinement of the two integrands' own knot sets. Takes the same
 * arguments as the cubature version above and is interchangeable with it, so
 * NcXcorSolver drives both methods through one cached-integrand path.
 */
void _nc_xcor_kernel_integrate_block_exact (NcXcor *xc, NcXcorKernelIntegrand *xclki1, NcXcorKernelIntegrand *xclki2, guint lmin, guint lmax, gboolean isauto, NcmVector *vp, NcmVector *vp_err);

/* Fails loudly when NcXcor:reltol asks for more than the kernel's closure carries. */
void _nc_xcor_check_kernel_tolerance (NcXcor *xc, NcXcorKernel *xclk);

/*
 * TRUE when there is a multipole at or above which both kernels run in the
 * Limber tier while their redshift supports do not overlap, so that their Cl
 * vanishes identically from there up and the kernel-space methods must not
 * integrate it; the lowest such multipole is written to @l_zero. FALSE when no
 * such multipole exists -- the supports overlap, or at least one kernel never
 * enters the Limber tier -- and @l_zero is then left untouched.
 *
 * A Limber kernel is supported only where xi = (l + 1/2) / k lies inside its
 * own radial range, so two disjoint bins have disjoint support in k and their
 * product is zero -- the same statement as the Limber-z tier's overlap test.
 * The non-Limber tier is the opposite case: there the two kernels couple only
 * through the outer k integral and disjoint bins do correlate, which is why
 * this must be asked per tier and not once per method.
 *
 * The tier is chosen per multipole, so this reports a threshold rather than a
 * bare yes/no: a caller whose range straddles @l_zero must zero only the tail
 * and integrate the head normally. Callers whose ranges never straddle a
 * kernel's l_limber (NcXcorSolver, whose plan_blocks() forces a block boundary
 * there) may simply compare @l_zero against the range's lmin.
 */
gboolean _nc_xcor_kernels_limber_disjoint (NcXcorKernel *xclk1, NcXcorKernel *xclk2, gboolean isauto, guint *l_zero);

G_END_DECLS

#endif /* _NC_XCOR_PRIV_H_ */

