/***************************************************************************
 *            gobj_itest.c
 *
 *  Mon Jun  6 16:51:56 2011
 *  Copyright  2011 Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * @file
 * @brief FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <glib.h>
#include <math.h>

gdouble
Fsin (gdouble x, gpointer p)
{
  return 1.0 * sin (M_PI * x * 1.0 + 0.0 * M_PI / 2.0) * exp(1.0 * x / 5.0) +
    1.0 * (1.0 * x * x / 2.0 + 5.0 * x + 0.1 * x * x * x + 0.0e-3 * gsl_pow_5 (x));
}

gdouble
Fjl (gdouble x, gpointer p)
{
  return ncm_sf_spherical_bessel (1000, x);
  //return sin (M_PI * x) * exp(x / 5.0);
}

static gdouble
Fxy (gdouble x, gdouble y, gpointer p)
{
  return (  1.0 +  2.0 * x +  3.0 * x * x +  4.0 * x * x * x +
          ( 5.0 +  6.0 * x +  7.0 * x * x +  8.0 * x * x * x ) * y +
          ( 9.0 + 10.0 * x + 11.0 * x * x + 12.0 * x * x * x ) * y * y +
          (13.0 + 14.0 * x + 15.0 * x * x + 16.0 * x * x * x ) * y * y * y
          ) * cos (x * y * 0.01) * exp (x * y);
}

static gdouble
Fx (gdouble x, gpointer p)
{
  return Fxy (x, 7.0, p);
}

static gdouble
Fy (gdouble y, gpointer p)
{
  return Fxy (2.5, y, p);
}

gint
main (gint argc, gchar *argv[])
{
  gint i, j;
  gsl_function F;
  F.function = &Fjl;
  F.params = NULL;
  NcmSpline2d *s2d;
  NcHICosmo *model;
  gdouble max_err = 0.0;

  ncm_cfg_init ();
  if (argc > 1)
  {
    model = nc_hicosmo_new_from_name (NC_TYPE_MODEL_DE, argv[1]);
    ncm_model_free (NCM_MODEL (model));
  }

  if (TRUE)
  {
    NcmSpline *s = ncm_sf_spherical_bessel_spline (2, 0.0, 1000.0, 1e-8);
    gdouble max_e = 0.0;
    guint maxi = 100000;
    GTimer *bench = g_timer_new ();

    for (i = 0; i <= maxi; i++)
    {
      gdouble x = 10.0 / (1.0 * maxi) * i;
      gdouble jl = ncm_sf_spherical_bessel (2, x);
      gdouble sjl = ncm_spline_eval (s, x);
      gdouble e = fabs ((jl - sjl) / jl);
      max_e = GSL_MAX (e, max_e);
    }
    printf ("# max err %e s length %u\n", max_e, ncm_vector_len (s->xv));

    g_timer_start (bench);
    for (i = 0; i <= maxi; i++)
    {
      gdouble x = 1000.0 / (1.0 * maxi) * i;
      max_e += ncm_sf_spherical_bessel (2, x);
    }
    printf ("# elapsed %es\n", g_timer_elapsed (bench, NULL) / maxi);

    g_timer_start (bench);
    for (i = 0; i <= maxi; i++)
    {
      gdouble x = 10.0 / (1.0 * maxi) * i;
      max_e += ncm_spline_eval (s, x);
    }
    printf ("# elapsed %es\n", g_timer_elapsed (bench, NULL) / maxi);


    exit (0);
  }

  if (FALSE)
  {
    //NcHICosmoQSpline *qspline = nc_hicosmo_qspline_new (ncm_spline_type_notaknot (NULL), 6, 2.0);
    NcHICosmoLCDM *lcdm = nc_hicosmo_lcdm_new ();
    NcDistance *dist = nc_distance_new (12.0);
    gdouble sum = 0.0;
    guint j;

    ncm_model_params_print_all (NCM_MODEL (lcdm), stdout);
    for (j = 0; j < 100; j++)
    {
      for (i = 0; i < 1000; i++)
      {
        gdouble zi = 10.0 / (1000 - 1.0) * i;
        gdouble cdi = nc_distance_comoving (dist, NC_HICOSMO (lcdm), zi);
        ncm_spline_prepare (dist->comoving_distance_spline->s);
        //printf ("% 20.15g % 20.15g\n", zi, cdi);
        sum += cdi;
      }
    }
    printf ("% 20.15g % 20.15g\n", 5.34, nc_distance_comoving (dist, NC_HICOSMO (lcdm), 5.34));
  }

  if (TRUE)
  {
    gdouble max_x = 5.0;
    gdouble min_x = 0.0;
    gdouble max_y = 14.0;
    gdouble min_y = 0.0;
    NcmVector *x, *y;
    //		NcMatrix *z;
    gint npx = 100;
    gint npy = 500;
    gsl_rng *rng = ncm_get_rng ();
    GTimer *bench = g_timer_new ();
    gsl_function gFx, gFy;
    gFx.function = &Fx;
    gFy.function = &Fy;

    x = ncm_vector_new (npx);
    y = ncm_vector_new (npy);
    //		z = nc_matrix_new (npy, npx);

    for (i = 0; i < npx; i++)
    {
      gdouble xi = (max_x - min_x) / (npx - 1.0) * i + min_x;
      ncm_vector_set (x, i, xi);
      //printf ("xi = %.5g\n", ncm_vector_get (x, i));
    }
    for (i = 0; i < npy; i++)
    {
      gdouble yi = (max_y - min_y) / (npy - 1.0) * i + min_y;
      ncm_vector_set (y, i, yi);
      //printf ("yi = %.5g\n", ncm_vector_get (y, i));
    }

    s2d = ncm_spline2d_bicubic_notaknot_new ();

    g_timer_start (bench);
    ncm_spline2d_set_function (s2d,
                               NCM_SPLINE_FUNCTION_SPLINE,
                               &gFx, &gFy, min_x, max_x, min_y, max_y, 1.0e-7);
    g_timer_stop (bench);
    x = s2d->xv;
    y = s2d->yv;

    ncm_vector_set (x, ncm_vector_len (x) - 2, ncm_vector_get (x, ncm_vector_len (x) - 3) + 1e-4);
    ncm_vector_set (y, ncm_vector_len (y) - 2, ncm_vector_get (y, ncm_vector_len (y) - 3) + 1e-4);

    npx = ncm_vector_len (x);
    npy = ncm_vector_len (y);
    printf ("# nx %d ny %d\n", npx, npy);

    for (i = 0; i < npy; i++)
    {
      gdouble yp = ncm_vector_get (y, i);
      for (j = 0; j < npx; j++)
      {
        gdouble xp = ncm_vector_get (x, j);
        gdouble g = Fxy (xp, yp, NULL);
        ncm_matrix_set (s2d->zm, i, j, g);
      }
    }

    //s2d = ncm_spline2d_new (ncm_spline2d_type_spsp (ncm_spline_type_gsl (gsl_interp_cspline)), x, y, z, TRUE);
    //s2d = ncm_spline2d_new (ncm_spline2d_type_spsp (ncm_spline_type_gsl (gsl_interp_nc_spline_notaknot)), x, y, z, TRUE);
    //s2d = ncm_spline2d_new (ncm_spline2d_type_spsp (ncm_spline_type_notaknot (NULL)), x, y, z, TRUE);
    //s2d = ncm_spline2d_new (ncm_spline2d_type_gsl (gsl_interp_cspline), x, y, z, TRUE);
    //s2d = ncm_spline2d_new (ncm_spline2d_type_gsl (gsl_interp_nc_spline_notaknot), x, y, z, TRUE);
    //s2d = ncm_spline2d_new (ncm_spline2d_type_bicubic (ncm_spline_type_notaknot (NULL)), x, y, z, FALSE);

    g_timer_continue (bench);
    for (i = 0; i < 10; i++)
      ncm_spline2d_prepare (s2d);
    g_timer_stop (bench);
    printf ("# time per prepare %g\n", g_timer_elapsed (bench, NULL) / 10.0);

    if (TRUE)
    {
      max_err = 0.0;
      for (i = 0; i < NCM_MATRIX_ROW_LEN (s2d->zm) - 1; i++)
      {
        const gdouble xp = (ncm_vector_get (s2d->xv, i + 1) + ncm_vector_get (s2d->xv, i)) * 0.5;
        for (j = 0; j < NCM_MATRIX_COL_LEN (s2d->zm) - 1; j++)
        {
          const gdouble yp = (ncm_vector_get (s2d->yv, j + 1) + ncm_vector_get (s2d->yv, j)) * 0.5;
          gdouble g = Fxy (xp, yp, NULL);
          gdouble gs, err;
          g_timer_continue (bench);
          gs = ncm_spline2d_eval (s2d, xp, yp);
          g_timer_stop (bench);
          err = g != 0.0 ? fabs ((gs-g) / g) : fabs (gs);

          if (gsl_finite (err))
            max_err = GSL_MAX (err, max_err);

          printf ("% 20.15g % 20.15g % 20.15g % 20.15g %.5e %.5e\n", xp, yp, g, gs, err, max_err);
        }
      }
    }

    max_err = 0.0;

    //for (i = 0; i < 10000; i++)
    {
      gdouble xp = (max_x - min_x) * gsl_rng_uniform (rng) + min_x;
      gdouble yp = (max_y - min_y) * gsl_rng_uniform (rng) + min_y;
      //			gdouble g = Fxy (xp, yp, NULL);
      gdouble gs1, gs2;

      gs1 = 0.0; gs2 = 0.0;
      gsl_rng_set (rng, 123);
      g_timer_start (bench);
      for (i = 0; i < 10000; i++)
      {
        gdouble x1 = (max_x - min_x) * gsl_rng_uniform (rng) + min_x;
        gdouble y1 = (max_y - min_y) * gsl_rng_uniform (rng) + min_y;

        gs1 += ncm_spline2d_integ_dx (s2d, min_x, max_x, y1);
        gs2 += ncm_spline2d_integ_dy (s2d, x1, min_y, max_y);
      }
      g_timer_stop (bench);
      printf ("#[ABC] % 20.15g % 20.15g % 20.15g int = % 20.15g % 20.15g\n", xp, yp, g_timer_elapsed (bench, NULL), gs1, gs2);

      gs1 = 0.0; gs2 = 0.0;
      gsl_rng_set (rng, 123);
      g_timer_start (bench);
      for (i = 0; i < 10000; i++)
      {
        gdouble x1 = (max_x - min_x) * gsl_rng_uniform (rng) + min_x;
        gdouble y1 = (max_y - min_y) * gsl_rng_uniform (rng) + min_y;

        gs1 += ncm_spline2d_integ_dx_spline_val (s2d, min_x, max_x, y1);
        gs2 += ncm_spline2d_integ_dy_spline_val (s2d, x1, min_y, max_y);
      }
      g_timer_stop (bench);
      printf ("#[ABC] % 20.15g % 20.15g % 20.15g int = % 20.15g % 20.15g\n", xp, yp, g_timer_elapsed (bench, NULL), gs1, gs2);

    }
    printf ("# max error %.5e %f\n", max_err, g_timer_elapsed (bench, NULL));
    g_timer_destroy (bench);
    ncm_vector_free (x);
    ncm_vector_free (y);
  }

  if (FALSE)
  {
    NcmSpline *s = ncm_spline_cubic_notaknot_new ();
    gdouble max_x = 1200.0;
    gdouble min_x = 1000.0;
    gint np = 1000;

    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, min_x, max_x, 1, 1.0e-4);
    if (FALSE)
    {
      ncm_vector_subfrom (s->xv, 0, 1.0e-1);
      ncm_vector_addto (s->xv, s->len - 1, 1.0e-1);
      ncm_vector_set (s->yv, 0, GSL_FN_EVAL (&F, ncm_vector_get (s->xv, 0)));
      ncm_vector_set (s->yv, s->len - 1, GSL_FN_EVAL (&F, ncm_vector_get (s->xv, s->len - 1)));
      min_x = ncm_vector_get (s->xv, 0);
      max_x = ncm_vector_get (s->xv, s->len - 1);
    }

    ncm_spline_prepare (s);
    fprintf (stderr, "# spline size %zd\n", s->len);
    printf ("# notaknot\n");
    max_err = 0.0;
    for (i = 0; i < np; i++)
    {
      gdouble x = (max_x - min_x) / (np - 1.0) * i + min_x;
      gdouble f = GSL_FN_EVAL (&F, x);
      gdouble fs = ncm_spline_eval (s, x);
      gdouble Ifs = ncm_spline_eval_integ (s, min_x, x);
      gdouble err = f != 0.0 ? fabs ((fs-f) / f) : fabs (fs);

      if (gsl_finite (err))
        max_err = GSL_MAX (err, max_err);

      printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 8.5e %8.5e\n", x, f, fs, Ifs, err, max_err);
    }
    fprintf (stderr, "# max error notaknot %.5e\n", max_err);
    printf ("# cspline\n");
    printf ("\n\n");
    ncm_spline_prepare (s);
    max_err = 0.0;
    for (i = 0; i < np; i++)
    {
      gdouble x = (max_x - min_x) / (np - 1.0) * i + min_x;
      gdouble f = GSL_FN_EVAL (&F, x);
      gdouble fs = ncm_spline_eval (s, x);
      gdouble Ifs = ncm_spline_eval_integ (s, min_x, x);
      gdouble err = f != 0.0 ? fabs ((fs-f) / f) : fabs (fs-f);

      if (gsl_finite (err))
        max_err = GSL_MAX (err, max_err);

      printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 8.5e %8.5e\n", x, f, fs, Ifs, err, max_err);
    }
    fprintf (stderr, "# max error cspline %.5e\n", max_err);
    printf ("# akima\n");
    printf ("\n\n");
    ncm_spline_prepare (s);
    max_err = 0.0;
    for (i = 0; i < np; i++)
    {
      gdouble x = (max_x - min_x) / (np - 1.0) * i + min_x;
      gdouble f = GSL_FN_EVAL (&F, x);
      gdouble fs = ncm_spline_eval (s, x);
      gdouble Ifs = ncm_spline_eval_integ (s, min_x, x);
      gdouble err = f != 0.0 ? fabs ((fs-f) / f) : fabs (fs-f);

      if (gsl_finite (err))
        max_err = GSL_MAX (err, max_err);

      printf ("% 20.15g % 20.15g % 20.15g % 20.15g % 8.5e %8.5e\n", x, f, fs, Ifs, err, max_err);
    }
    fprintf (stderr, "# max error akima %.5e\n", max_err);
  }

  if (FALSE)
  {
    GTimer *bench = g_timer_new ();
    gsl_function F;
    gint i;
    gsize tests = 10000;
    gsl_rng *r = ncm_get_rng ();
    gsl_vector *x = gsl_vector_alloc (tests);
    gsl_vector *val = gsl_vector_alloc (tests);
    F.function = &Fjl;
    for (i = 0; i < tests; i++)
      gsl_vector_set (x, i, gsl_rng_uniform (r) * pow (10.0, 4.0 * gsl_rng_uniform (r)));
    printf ("# Go threads %d\n", NCM_THREAD_POOL_MAX);
    fflush (stdout);
    g_timer_reset (bench);
    for (i = 0; i < x->size; i++)
      val->data[i] = GSL_FN_EVAL (&F, x->data[i]);
    printf ("# Direct took %fs\n", g_timer_elapsed (bench, NULL));
    for (i = 0; i < x->size; i++)
      printf ("% 20.15g % 20.15g\n", x->data[i], val->data[i]);
  }

  return 0;
}
