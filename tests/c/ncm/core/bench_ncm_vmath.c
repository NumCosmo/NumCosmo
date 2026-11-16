/***************************************************************************
 *            bench_ncm_vmath.c
 *
 *  Sun Nov 15 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * The time per value of the functions of ncm_vmath_private.h in a loop that vectorizes,
 * against libm in the same loop, for loops of 2 to 1025 values: whether the code of the
 * header is worth using on a target (.github/workflows/arch_check.yml). It needs only
 * GLib and the header, so it is built by the compiler of the target alone:
 *
 *   cc -O3 -fopenmp-simd -DNCM_VMATH_FORCE_VECTOR_EXP_LOG -I numcosmo \
 *     tests/c/ncm/core/bench_ncm_vmath.c $(pkg-config --cflags --libs glib-2.0) -lm
 *
 * Without NCM_VMATH_FORCE_VECTOR_EXP_LOG, exp, log and pow_pos are those the build of
 * the library would use on that target.
 */

#include <glib.h>
#include <math.h>
#include <stdio.h>
#include "ncm/core/ncm_vmath_private.h"

#define BENCH_TOTAL (4000000)
#define BENCH_REPEAT (5)

typedef void (*BenchLoop) (const guint n, const gdouble *x, const gdouble *y, gdouble *out);

#define BENCH_LOOP(name, expr)                                                 \
        static void __attribute__ ((noinline))                                 \
        name (const guint n, const gdouble *x, const gdouble *y, gdouble *out) \
        {                                                                      \
          guint i;                                                             \
          _Pragma("omp simd")                                                  \
          for (i = 0; i < n; i++)                                              \
          out[i] = (expr);                                                     \
          }

BENCH_LOOP (_vm_exp, ncm_vmath_exp (x[i]))
BENCH_LOOP (_lm_exp, exp (x[i]))
BENCH_LOOP (_vm_log, ncm_vmath_log (y[i]))
BENCH_LOOP (_lm_log, log (y[i]))
BENCH_LOOP (_vm_pow, ncm_vmath_pow_pos (y[i], x[i]))
BENCH_LOOP (_lm_pow, pow (y[i], x[i]))
BENCH_LOOP (_vm_atan2, ncm_vmath_atan2 (x[i], y[i]))
BENCH_LOOP (_lm_atan2, atan2 (x[i], y[i]))

/* sin + cos of x: the sum keeps one output per value, as the loops above */
static void __attribute__ ((noinline))

_vm_sincos (const guint n, const gdouble *x, const gdouble *y, gdouble *out)
{
  guint i;

  #pragma omp simd

  for (i = 0; i < n; i++)
  {
    gdouble s, c;

    ncm_vmath_sincos (x[i], &s, &c);
    out[i] = s + c;
  }
}

static void __attribute__ ((noinline))

_lm_sincos (const guint n, const gdouble *x, const gdouble *y, gdouble *out)
{
  guint i;

  #pragma omp simd

  for (i = 0; i < n; i++)
    out[i] = sin (x[i]) + cos (x[i]);
}

/* The least time per value, in ns, of BENCH_REPEAT passes of BENCH_TOTAL values */
static gdouble
_bench (BenchLoop f, const guint n, const gdouble *x, const gdouble *y, gdouble *out, gdouble *sink)
{
  const guint reps = MAX (1, BENCH_TOTAL / n);
  gdouble best     = G_MAXDOUBLE;
  guint k, r;

  for (k = 0; k < BENCH_REPEAT; k++)
  {
    const gint64 t0 = g_get_monotonic_time ();

    for (r = 0; r < reps; r++)
    {
      f (n, x, y, out);
      *sink += out[r % n];
    }

    best = MIN (best, 1.0e3 * (gdouble) (g_get_monotonic_time () - t0) / ((gdouble) reps * n));
  }

  return best;
}

gint
main (void)
{
  const guint ns[]     = {2, 5, 9, 17, 33, 65, 257, 1025};
  const gchar *names[] = {"exp", "log", "pow_pos", "atan2", "sincos"};
  const BenchLoop vm[] = {_vm_exp, _vm_log, _vm_pow, _vm_atan2, _vm_sincos};
  const BenchLoop lm[] = {_lm_exp, _lm_log, _lm_pow, _lm_atan2, _lm_sincos};
  gdouble *x           = g_new (gdouble, 1025);
  gdouble *y           = g_new (gdouble, 1025);
  gdouble *out         = g_new (gdouble, 1025);
  gdouble sink         = 0.0;
  guint f, k, i;

  /* The ranges of the integrals: exponents of a few units, bases in (0, 2] */
  for (i = 0; i < 1025; i++)
  {
    x[i] = -8.0 + 16.0 * ((i * 7919) % 1025) / 1024.0;
    y[i] = pow (10.0, -6.0 + 6.3 * ((i * 104729) % 1025) / 1024.0);
  }

  printf ("compiler: %s\n", __VERSION__);
  printf ("target:%s%s%s%s%s\n",
#if defined (__x86_64__)
          " x86_64",
#else
          "",
#endif
#if defined (__AVX512F__)
          " AVX-512",
#elif defined (__AVX2__)
          " AVX2",
#else
          "",
#endif
#if defined (__aarch64__)
          " aarch64",
#else
          "",
#endif
#if defined (__ARM_NEON)
          " NEON",
#else
          "",
#endif
#if defined (NCM_VMATH_HAVE_VECTOR_EXP_LOG)
          " (exp, log, pow_pos: code of the header)"
#else
          " (exp, log, pow_pos: libm)"
#endif
  );
  printf ("%-8s %6s %10s %10s %8s\n", "function", "n", "vmath ns", "libm ns", "speed-up");

  for (f = 0; f < G_N_ELEMENTS (names); f++)
  {
    for (k = 0; k < G_N_ELEMENTS (ns); k++)
    {
      const gdouble tv = _bench (vm[f], ns[k], x, y, out, &sink);
      const gdouble tl = _bench (lm[f], ns[k], x, y, out, &sink);

      printf ("%-8s %6u %10.3f %10.3f %8.2f\n", names[f], ns[k], tv, tl, tl / tv);
    }
  }

  /* Keeps the loops from being removed */
  if (sink == 0.123456789)
    printf ("%g\n", sink);

  g_free (x);
  g_free (y);
  g_free (out);

  return 0;
}

