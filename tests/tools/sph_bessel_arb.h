/***************************************************************************
 *            sph_bessel_arb.h
 *
 *  Mon September 28 12:00:00 2026
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
 * Spherical Bessel functions on complex balls, shared by the reference
 * generators in this directory. Header-only and all-static.
 */

#ifndef _SPH_BESSEL_ARB_H
#define _SPH_BESSEL_ARB_H

#if defined (__has_include )
#if __has_include (<flint/acb.h>)
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#else
#include <acb.h>
#include <acb_hypgeom.h>
#endif
#else
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#endif

/*
 * j_ell (z), by whichever form is well conditioned at this z.
 *
 * Near the origin: the entire form
 *
 *   j_ell (z) = sqrt (pi) / 2^(ell+1) * z^ell * 0F1~ (; ell + 3/2; -z^2/4)
 *
 * The textbook sqrt (pi / 2z) J_{ell+1/2} (z) has a removable singularity at
 * z = 0, and a window whose support reaches the observer puts z = 0 inside the
 * domain, where the ball arithmetic sees 1/0 and subdivides without end.
 *
 * Away from it: the textbook form. The 0F1 series needs about z^2/4 terms and
 * loses about that many bits to cancellation, which is affordable at z ~ 100
 * and is not a method at all beyond it -- the C_ell integrand reaches
 * z = k chi ~ 4400 at ell = 200, where it would want 5 million bits. Arb
 * evaluates J_nu by asymptotic expansion out there, at no precision penalty.
 *
 * The switch is on a certified lower bound of |z|, so a ball that straddles
 * the threshold takes the entire form and stays correct either way.
 */
static void
sph_bessel (acb_t out, const acb_t z, long ell, slong prec)
{
  acb_t nu, J, t;
  arb_t az;
  arf_t lb;
  int far;

  arb_init (az);
  arf_init (lb);
  acb_abs (az, z, prec);
  arb_get_lbound_arf (lb, az, prec);
  far = arf_cmp_d (lb, 4.0) > 0;
  arb_clear (az);
  arf_clear (lb);

  acb_init (nu);
  acb_init (J);
  acb_init (t);

  if (far)
  {
    /* sqrt (pi / (2 z)) J_{ell + 1/2} (z) */
    acb_set_si (nu, 2 * ell + 1);
    acb_div_si (nu, nu, 2, prec);
    acb_hypgeom_bessel_j (J, nu, z, prec);
    acb_const_pi (t, prec);
    acb_div (t, t, z, prec);
    acb_div_si (t, t, 2, prec);
    acb_sqrt (t, t, prec);
    acb_mul (out, J, t, prec);
  }
  else
  {
    acb_set_si (nu, 2 * ell + 3);
    acb_div_si (nu, nu, 2, prec); /* ell + 3/2 */
    acb_sqr (J, z, prec);
    acb_div_si (J, J, -4, prec); /* -z^2/4 */
    acb_hypgeom_0f1 (J, nu, J, 1, prec);
    acb_pow_ui (t, z, (ulong) ell, prec);
    acb_mul (J, J, t, prec);
    acb_const_pi (t, prec);
    acb_sqrt (t, t, prec);
    acb_mul (J, J, t, prec);
    acb_mul_2exp_si (J, J, -(ell + 1));
    acb_set (out, J);
  }

  acb_clear (nu);
  acb_clear (J);
  acb_clear (t);
}

#endif /* _SPH_BESSEL_ARB_H */

