/***************************************************************************
 *            nc_xcor_kernel_analytic_arb.c
 *
 *  Tue August 26 12:00:00 2026
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
 * Certified reference values for the radial integral of the analytic xcor
 * windows,
 *
 *   I_ell(k) = int W(chi) j_ell(k chi) dchi ,   W normalized to int W dchi = 1
 *
 * over the window's own truncated support, with chi in Mpc.
 *
 * Re-implements each window of numcosmo/nc/xcor/nc_xcor_kernel_radial_*.c
 * in Arb ball arithmetic, independently of the library. The normalization is
 * *not* copied from the library's closed form -- it is computed here as a
 * certified integral of the unnormalized shape, so comparing against the
 * library checks its normalization too.
 *
 * Working precision is not fixed: every entry is recomputed at doubling
 * precision until its relative radius clears the target. That is what makes
 * the tool usable in the deep sub-band regime, where a fixed-precision
 * reference returns a plausible wrong number instead of reporting failure.
 *
 * Emits TSV on stdout: shape ell k value radius prec.
 *
 * Build (standalone):
 *   gcc -O2 -o xcor_window_arb nc_xcor_kernel_analytic_arb.c \
 *       $(pkg-config --cflags --libs flint) -lm
 */


#include "xcor_window_arb.h"

int
main (int argc, char **argv)
{
  Par p;
  int window_mode;
  double target = 1.0e-25;
  long ell      = 2;
  double k;
  int i;

  par_init (&p);
  window_mode = 0;
  p.shape     = SHAPE_GAUSS;
  p.n_sigma   = 4.0;
  p.n_scale   = 6.0;

  for (i = 1; i < argc; i++)
  {
    const char *a = argv[i];

    if (strncmp (a, "--shape=", 8) == 0)
    {
      p.shape = parse_shape (a + 8);
    }
    else if (strncmp (a, "--ell=", 6) == 0)
    {
      ell = atol (a + 6);
    }
    else if (strncmp (a, "--chi-mean=", 11) == 0)
    {
      p.chi_mean = atof (a + 11);
    }
    else if (strncmp (a, "--chi-sigma=", 12) == 0)
    {
      p.chi_sigma = atof (a + 12);
    }
    else if (strncmp (a, "--n-sigma=", 10) == 0)
    {
      p.n_sigma = atof (a + 10);
    }
    else if (strncmp (a, "--chi-lower=", 12) == 0)
    {
      p.chi_lower = atof (a + 12);
    }
    else if (strncmp (a, "--chi-upper=", 12) == 0)
    {
      p.chi_upper = atof (a + 12);
    }
    else if (strncmp (a, "--chi-scale=", 12) == 0)
    {
      p.chi_scale = atof (a + 12);
    }
    else if (strncmp (a, "--nu=", 5) == 0)
    {
      p.nu = atof (a + 5);
    }
    else if (strncmp (a, "--n-scale=", 10) == 0)
    {
      p.n_scale = atof (a + 10);
    }
    else if (strncmp (a, "--alpha=", 8) == 0)
    {
      p.alpha = atof (a + 8);
    }
    else if (strncmp (a, "--beta=", 7) == 0)
    {
      p.beta = atof (a + 7);
    }
    else if (strncmp (a, "--chi-source-lower=", 19) == 0)
    {
      p.chi_source_lower = atof (a + 19);
    }
    else if (strncmp (a, "--chi-source-upper=", 19) == 0)
    {
      p.chi_source_upper = atof (a + 19);
    }
    else if (strncmp (a, "--mu=", 5) == 0)
    {
      parse_list (a + 5, p.mu, &p.n_bumps);
    }
    else if (strncmp (a, "--sigma=", 8) == 0)
    {
      int m;

      parse_list (a + 8, p.sg, &m);
    }
    else if (strncmp (a, "--weight=", 9) == 0)
    {
      int m;

      parse_list (a + 9, p.wt, &m);
    }
    else if (strncmp (a, "--target-rel=", 13) == 0)
    {
      target = atof (a + 13);
    }
    else if (strcmp (a, "--window") == 0)
    {
      window_mode = 1;
    }
    else
    {
      fprintf (stderr, "unknown argument '%s'\n", a);
      exit (1);
    }
  }

  p.ell = ell;
  shape_support (&p);

  /* The normalization is measured, not copied from the library. */
  {
    acb_t nrm;
    char *s;

    acb_init (nrm);
    p.with_bessel = 0;
    certified (nrm, &p, target);
    s = arb_get_str (acb_realref (nrm), 30, ARB_STR_NO_RADIUS);
    printf ("# shape=%s ell=%ld support=[%.17g,%.17g] norm=%s\n",
            shape_names[p.shape], p.ell, p.chi_min, p.chi_max, s);
    flint_free (s);

    if (window_mode)
    {
      /* W(chi) = W_u(chi) / norm, for chi values on stdin. The window is a
       * closed form, so only the normalization carries a radius. */
      double chi;

      printf ("# shape\tchi\tW\n");

      while (scanf ("%lf", &chi) == 1)
      {
        acb_t x, Wu, Wn;

        acb_init (x);
        acb_init (Wu);
        acb_init (Wn);
        acb_set_d (x, chi);

        if ((chi < p.chi_min) || (chi > p.chi_max))
        {
          acb_zero (Wn);
        }
        else
        {
          window_u (Wu, x, &p, 0, 256);
          acb_div (Wn, Wu, nrm, 256);
        }

        s = arb_get_str (acb_realref (Wn), 30, ARB_STR_NO_RADIUS);
        printf ("%s\t%.17g\t%s\n", shape_names[p.shape], chi, s);
        flint_free (s);
        acb_clear (x);
        acb_clear (Wu);
        acb_clear (Wn);
      }

      acb_clear (nrm);
      flint_cleanup ();

      return 0;
    }

    printf ("# shape\tell\tk\tvalue\tradius\tprec\n");
    p.with_bessel = 1;

    while (scanf ("%lf", &k) == 1)
    {
      acb_t res, val;
      slong prec;

      acb_set_d (p.k, k);
      acb_init (res);
      acb_init (val);
      prec = certified (res, &p, target);
      acb_div (val, res, nrm, prec);
      s = arb_get_str (acb_realref (val), 30, ARB_STR_NO_RADIUS);
      printf ("%s\t%ld\t%.17g\t%s\t%.4e\t%ld\n",
              shape_names[p.shape], p.ell, k, s,
              mag_get_d (arb_radref (acb_realref (val))), (long) prec);
      flint_free (s);
      acb_clear (res);
      acb_clear (val);
    }

    acb_clear (nrm);
  }

  flint_cleanup ();

  return 0;
}

