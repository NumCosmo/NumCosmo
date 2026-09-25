/***************************************************************************
 *            ncm_complex.c
 *
 *  Thu September 25 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_complex.c
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
 * NcmComplex:
 *
 * Boxed complex double.
 *
 * In C it is `complex double`, so a #NcmComplex pointer can be passed to code using
 * C99 complex numbers or `fftw_complex`. The boxed type makes complex values
 * available through GObject introspection.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_complex.h"

G_DEFINE_BOXED_TYPE (NcmComplex, ncm_complex, ncm_complex_dup, ncm_complex_free)

/**
 * ncm_complex_new:
 *
 * Allocates a complex number set to zero.
 *
 * Returns: (transfer full): a new #NcmComplex.
 */
NcmComplex *
ncm_complex_new ()
{
  return g_new0 (NcmComplex, 1);
}

/**
 * ncm_complex_dup:
 * @c: a #NcmComplex
 *
 * Returns: (transfer full): a newly allocated copy of @c.
 */
NcmComplex *
ncm_complex_dup (NcmComplex *c)
{
  NcmComplex *cc = ncm_complex_new ();

  *cc = *c;

  return cc;
}

/**
 * ncm_complex_free:
 * @c: a #NcmComplex
 *
 * Frees @c, which must come from ncm_complex_new() or ncm_complex_dup().
 */
void
ncm_complex_free (NcmComplex *c)
{
  g_free (c);
}

/**
 * ncm_complex_clear:
 * @c: a #NcmComplex
 *
 * Frees *@c, as ncm_complex_free(), and sets *@c to %NULL.
 */
void
ncm_complex_clear (NcmComplex **c)
{
  g_clear_pointer (c, g_free);
}

/**
 * ncm_complex_set:
 * @c: a #NcmComplex
 * @a: the real part $a$
 * @b: the imaginary part $b$
 *
 * Sets @c to $a + i b$.
 */
/**
 * ncm_complex_set_c: (skip)
 * @c: a #NcmComplex
 * @z: a complex double
 *
 * Sets @c to @z.
 */
/**
 * ncm_complex_set_zero:
 * @c: a #NcmComplex
 *
 * Sets @c to zero.
 */
/**
 * ncm_complex_Re:
 * @c: a #NcmComplex
 *
 * Returns: the real part of @c.
 */
/**
 * ncm_complex_Im:
 * @c: a #NcmComplex
 *
 * Returns: the imaginary part of @c.
 */
/**
 * ncm_complex_Abs:
 * @c: a #NcmComplex
 *
 * Returns: $|c|$.
 */
/**
 * ncm_complex_c: (skip)
 * @c: a #NcmComplex
 *
 * Returns: @c as a complex double.
 */

/**
 * ncm_complex_res_add_mul_real:
 * @c1: a #NcmComplex
 * @c2: a #NcmComplex
 * @v: a double
 *
 * Sets $c_1 \to c_1 + c_2 v$. @c1 and @c2 must not overlap.
 */
/**
 * ncm_complex_res_add_mul:
 * @c1: a #NcmComplex
 * @c2: a #NcmComplex
 * @c3: a #NcmComplex
 *
 * Sets $c_1 \to c_1 + c_2 c_3$. @c1 must not overlap @c2 or @c3.
 */

/**
 * ncm_complex_mul_real:
 * @c: a #NcmComplex
 * @v: a double
 *
 * Sets $c \to c v$.
 */
/**
 * ncm_complex_res_mul:
 * @c1: a #NcmComplex
 * @c2: a #NcmComplex
 *
 * Sets $c_1 \to c_1 c_2$. @c1 and @c2 must not overlap.
 */

