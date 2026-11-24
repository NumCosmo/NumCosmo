/***************************************************************************
 *            ncm_vmath_private.h
 *
 *  Tue Nov 17 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_vmath_private.h
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NCM_VMATH_PRIVATE_H_
#define _NCM_VMATH_PRIVATE_H_

#include <glib.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <string.h>

G_BEGIN_DECLS

/*
 * Declares the vector variants of exp(), log(), pow(), sin(), cos() and atan2() that glibc
 * provides (libmvec), so that GCC calls them in the loops it vectorizes. The functions keep
 * their names. A file opts in by including this header, and every loop of that file that
 * GCC vectorizes is affected, with or without omp simd.
 *
 * glibc declares these variants itself only under -ffast-math, which NumCosmo does not use.
 * The declarations below are the same, without the other effects of that option.
 *
 * Requirements, all tested here: GCC (clang ignores the declarations), x86-64, OpenMP, and
 * glibc 2.22 or later for exp(), log(), pow(), sin() and cos(), 2.35 or later for atan2().
 * When one is not met nothing is declared and every call is the one of libm.
 *
 * Effects on the including file, where the variants are used:
 * - a value computed in a vectorized loop can differ from the one of libm in the last bits
 *   (tests/c/ncm/core/test_ncm_vmath.c bounds the difference), so results depend on the
 *   build;
 * - the elements left after the vector body of a loop are computed by libm, so a value can
 *   depend on its position in the loop;
 * - the variants return the zeros, infinities and NaN of libm (/ncm/vmath/special of that
 *   test); they do not set errno nor the floating-point exception flags, and assume the
 *   default rounding mode.
 */
#if defined (__GNUC__) && !defined (__clang__) && defined (__x86_64__) && defined (__GLIBC__) && \
  defined (_OPENMP) && !defined (__FAST_MATH__)

#if __GLIBC_PREREQ (2, 22)
#pragma omp declare simd notinbranch
double exp (double);

#pragma omp declare simd notinbranch
double log (double);

#pragma omp declare simd notinbranch
double pow (double, double);

#pragma omp declare simd notinbranch
double sin (double);

#pragma omp declare simd notinbranch
double cos (double);

#endif /* __GLIBC_PREREQ (2, 22) */

#if __GLIBC_PREREQ (2, 35)
#pragma omp declare simd notinbranch
double atan2 (double, double);

#endif /* __GLIBC_PREREQ (2, 35) */

#endif /* GCC, x86-64, glibc and OpenMP */

/*
 * ncm_vmath_select:
 * @cond: a condition
 * @a: value for a true @cond
 * @b: value for a false @cond
 *
 * Returns @cond ? @a : @b as a bitwise blend of the two values, both of which are
 * computed. GCC keeps a conditional expression on doubles as control flow when one of its
 * branches can raise a floating-point exception (a division, for instance), which leaves
 * the loop scalar; with the blend the loop vectorizes. The value not chosen may be
 * infinity or NaN, and the blend discards it.
 */
static inline gdouble
ncm_vmath_select (const gboolean cond, const gdouble a, const gdouble b)
{
  const guint64 mask = -(guint64) (cond != 0);
  guint64 ua, ub, ur;
  gdouble r;

  memcpy (&ua, &a, sizeof (ua));
  memcpy (&ub, &b, sizeof (ub));
  ur = (ua &mask) | (ub & ~mask);
  memcpy (&r, &ur, sizeof (r));

  return r;
}

G_END_DECLS

#endif /* NUMCOSMO_GIR_SCAN */

#endif /* _NCM_VMATH_PRIVATE_H_ */

