/***************************************************************************
 *            cmb_maps.c
 *
 *  Thu May  8 11:50:01 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf_coupling.h>
#if defined HAVE_CHEALPIX && defined HAVE_CFITSIO
#include <fitsio.h>
#include <chealpix.h>
#endif /* defined HAVE_CHEALPIX && defined HAVE_CFITSIO */

gint
main (void)
{
#if defined HAVE_CHEALPIX && defined HAVE_CFITSIO
  GTimer *bench          = g_timer_new ();
  NcSphereMapAlm *mapalm = nc_sphere_mapalm_new ();
  NcSphereMap *omap      = nc_sphere_healpix_read_map (argv[1], NULL);
  NcSphereMap *map       = nc_sphere_map_clone (omap);

#ifdef HAVE_FFTW3
  GTimer *global_bench = g_timer_new ();
  NcSphereMapSHT *mapsht;
  gsl_matrix *cls;

#endif /* HAVE_FFTW3 */
  gint i = 1000;
  gdouble noise;

  printf ("# Loading of [%ld] pixels took %fs\n", map->npix, g_timer_elapsed (bench, NULL));

  nc_cfg_init ();

  if (FALSE)
  {
    g_timer_start (bench);
    printf ("# <<%ld>> mean: %g, sd: %g, k: %g, ac: %g, skew: %g\n", map->npix,
            gsl_stats_float_mean (map->dt->data, 1, map->npix),
            gsl_stats_float_sd (map->dt->data, 1, map->npix),
            gsl_stats_float_kurtosis (map->dt->data, 1, map->npix),
            gsl_stats_float_lag1_autocorrelation (map->dt->data, 1, map->npix),
            gsl_stats_float_skew (map->dt->data, 1, map->npix)
           );
    printf ("# Stats of [%ld] pixels took %fs\n", map->npix, g_timer_elapsed (bench, NULL));
  }

  i = 1;

  while (i--)
  {
    g_timer_start (bench);
    nc_sphere_map_set_order (map, NC_SPHERE_MAP_ORDER_RING, FALSE);
    nc_sphere_map_set_order (omap, NC_SPHERE_MAP_ORDER_RING, FALSE);
    printf ("# Convert two maps to %s [%ld] pixels took %fs\n", i % 2 ? "ring" : "nest", map->npix, g_timer_elapsed (bench, NULL));
  }

  g_timer_start (bench);
  noise = nc_sphere_map_homogenize_noise (map, NC_C_WMAP5_COADDED_I_W);
  printf ("# Homogenize map with [%ld] pixels took %fs, noise %g\n", map->npix, g_timer_elapsed (bench, NULL), noise);

  nc_sphere_map_rotate_avg (map, 1000);

  g_timer_start (bench);
  nc_sphere_mapalm_init (mapalm, 1024);
  printf ("# New mapalm with [%ld] pixels and lmax [%d] took %fs\n", map->npix, mapalm->lmax, g_timer_elapsed (bench, NULL));

#ifdef HAVE_FFTW3
  g_timer_start (bench);
  mapsht = nc_sphere_mapsht_new (map, mapalm, FFTW_ESTIMATE); /*FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE */
  printf ("# New SHT with [%ld] pixels and lmax [%d] took %fs\n", map->npix, mapalm->lmax, g_timer_elapsed (bench, NULL));

#define SIMS 1000

  cls = gsl_matrix_alloc (SIMS, mapalm->lmax + 1);

  for (i = 0; i < SIMS; i++)
  {
    glong l;
    gsl_vector_view row = gsl_matrix_row (cls, i);

    g_timer_start (bench);
    nc_sphere_mapsht_map2alm (mapsht, 0.0);
    printf ("# map2alm with [%ld] pixels and lmax [%d] took %fs\n", map->npix, mapalm->lmax, g_timer_elapsed (bench, NULL));

    g_timer_start (bench);
    nc_sphere_map_copy (map, omap);
    nc_sphere_map_homogenize_noise (map, NC_C_WMAP5_COADDED_I_W);
    printf ("# Copy and homogenize map with [%ld] pixels took %fs\n", map->npix, g_timer_elapsed (bench, NULL));

    gsl_vector_memcpy (&row.vector, mapalm->Nc);

    for (l = 0; l <= mapalm->lmax; l++)
    {
      gsl_vector_view col = gsl_matrix_column (cls, l);
      gdouble cl          = gsl_stats_mean (col.vector.data, col.vector.stride, i + 1) * 1.0e6;
      gdouble sigma_cl    = gsl_stats_sd (col.vector.data, col.vector.stride, i + 1) * 1.0e6;

      printf ("%ld %g %g %g\n", l, cl, l * (l + 1.0) / (2.0 * M_PI) * (cl - noise * noise * 0.0), l * (l + 1.0) / (2.0 * M_PI) * sigma_cl);
      fflush (stdout);
    }

    printf ("\n\n");
  }

/*  g_timer_start (bench); */
/*  nc_sphere_mapsht_alm2map (mapsht); */
/*  printf ("# alm2map with [%ld] pixels and lmax [%d] took %fs\n", map->npix, mapalm->lmax, g_timer_elapsed (bench, NULL)); */

  if (FALSE)
  {
    g_timer_start (bench);
    printf ("# <<%ld>> mean: %g, sd: %g, k: %g, ac: %g, skew: %g\n", map->npix,
            gsl_stats_float_mean (map->dt->data, 1, map->npix),
            gsl_stats_float_sd (map->dt->data, 1, map->npix),
            gsl_stats_float_kurtosis (map->dt->data, 1, map->npix),
            gsl_stats_float_lag1_autocorrelation (map->dt->data, 1, map->npix),
            gsl_stats_float_skew (map->dt->data, 1, map->npix)
           );
    printf ("# Stats of [%ld] pixels took %fs\n", map->npix, g_timer_elapsed (bench, NULL));
  }

  nc_sphere_map_set_order (map, NC_SPHERE_MAP_ORDER_NEST, FALSE);
  nc_sphere_healpix_write_map (map, "map_out.fits", TRUE);

  printf ("# whole process with [%ld] pixels and lmax [%d] took %fs\n", map->npix, mapalm->lmax, g_timer_elapsed (global_bench, NULL));

  if (TRUE)
  {
    glong l, m, i, j = 0;

    for (l = 0; l <= mapalm->lmax && TRUE; l++)
    {
      printf ("%.5ld % .12e\n", l, l * (l + 1.0) / (2.0 * M_PI) * gsl_vector_get (mapalm->Nc, l));

      /*printf ("#%.5ld % .12e\n", l, gsl_vector_get (mapalm->Nc, l)); */
      /*printf ("%.4ld % .12e\n", l, gsl_vector_get (mapalm->Nc, l)); */
      for (m = 0; m <= l && FALSE; m++)
      {
        gsl_complex t = gsl_vector_complex_get (mapalm->alm, NC_MAP_ALM_INDEX (mapalm->lmax, l, m));

        printf ("alm[%.7ld][%.4ld, %.4ld] = % .16e + % .16e i\n", j++, l, m,  GSL_REAL (t), GSL_IMAG (t));
      }
    }

    for (i = 0; i < map->npix && FALSE; i++)
      printf ("[%ld] = % .12e\n", i, map->dt->data[i]);

    /*printf ("[%ld] (% .12e, % .12e) = % .12e\n", i, gsl_vector_get (map->theta, i), gsl_vector_get (map->phi, i), map->fmap[i]); */
  }

#endif /* HAVE_FFTW3 */
#else
  g_error ("chealpix or cfitsio not installed.");
#endif /* defined HAVE_CHEALPIX && defined HAVE_CFITSIO */

  return 0;
}

