/***************************************************************************
 *            map.c
 *
 *  Tue Jun 24 16:30:59 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * SECTION:map
 * @title: Spherical Shell Map
 * @short_description: Object representing a spherical shell map
 *
 * Map manipulation algorithms, Ylm decomposition.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "sphere/map.h"
#include "sphere/healpix.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"
#include "math/util.h"

#include <string.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_legendre.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */

/**
 * ncm_sphere_map_new: (skip)
 * @nside: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmSphereMap *
ncm_sphere_map_new (gint nside)
{
  NcmSphereMap *map = g_slice_new (NcmSphereMap);
  memset (map, 0, sizeof(NcmSphereMap));
  map->loaded = FALSE;
  map->is_init_coord = FALSE;
  map->type = 0;
  map->nside = nside;
  map->npix = HEALPIX_NPIX(map->nside);
  map->dt          = NULL;// gsl_vector_float_alloc (map->npix);
  map->qpol        = NULL;// gsl_vector_float_alloc (map->npix);
  map->upol        = NULL;// gsl_vector_float_alloc (map->npix);
  map->spur_signal = NULL;// gsl_vector_float_alloc (map->npix);
  map->nobs        = NULL;// gsl_vector_float_alloc (map->npix);
  return map;
}

/**
 * ncm_sphere_map_copy:
 * @dest: a #NcmSphereMap
 * @orig: a #NcmSphereMap
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_sphere_map_copy (NcmSphereMap *dest, NcmSphereMap *orig)
{
  gsl_vector_float **ovectors[] = {&orig->dt, &orig->qpol, &orig->upol, &orig->spur_signal, &orig->nobs};
  gsl_vector_float **dvectors[] = {&dest->dt, &dest->qpol, &dest->upol, &dest->spur_signal, &dest->nobs};
  gint i;
  g_assert (orig->nside == dest->nside);
  ncm_sphere_map_set_order (dest, orig->order, FALSE);

  for (i = 0; i < 5; i++)
  {
    if (*ovectors[i] != NULL)
    {
      if (*dvectors[i] == NULL)
        *dvectors[i] = gsl_vector_float_alloc (dest->npix);
      gsl_vector_float_memcpy (*dvectors[i], *ovectors[i]);
    }
  }
  dest->loaded = orig->loaded;
  return TRUE;
}

/**
 * ncm_sphere_map_clone: (skip)
 * @map: a #NcmSphereMap
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmSphereMap *
ncm_sphere_map_clone (NcmSphereMap *map)
{
  NcmSphereMap *clone = ncm_sphere_map_new (map->nside);
  ncm_sphere_map_copy (clone, map);
  clone->loaded = TRUE;
  return clone;
}

/**
 * ncm_sphere_map_init_coord:
 * @map: a #NcmSphereMap
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_sphere_map_init_coord (NcmSphereMap *map)
{
  glong i;
  if (map->theta == NULL)
    map->theta = gsl_vector_alloc (map->npix);
  if (map->phi == NULL)
    map->phi = gsl_vector_alloc (map->npix);

  if (map->order == NC_SPHERE_MAP_ORDER_NEST)
    for (i = 0; i < map->npix; i++)
      ncm_sphere_healpix_pix2ang_nest (map->nside, i, gsl_vector_ptr (map->theta, i), gsl_vector_ptr (map->phi, i));
  else
    for (i = 0; i < map->npix; i++)
  {
      ncm_sphere_healpix_pix2ang_ring (map->nside, i, gsl_vector_ptr (map->theta, i), gsl_vector_ptr (map->phi, i));
//    printf ("[%.6ld] %.16f\n", i, map->theta->data[i]);
  }
  map->is_init_coord = TRUE;
  return TRUE;
}

/**
 * ncm_sphere_map_set_order:
 * @map: a #NcmSphereMap
 * @order: a #NcmSphereMapOrder
 * @init_coord: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_sphere_map_set_order (NcmSphereMap *map, NcmSphereMapOrder order, gboolean init_coord)
{
  gsl_vector_float **vectors[] = {&map->dt, &map->qpol, &map->upol, &map->spur_signal, &map->nobs};
  if (map->order != order)
  {
    gsl_vector_float *temp_pix = gsl_vector_float_alloc (map->npix);
    glong (*convert)(gint,glong) = (order == NC_SPHERE_MAP_ORDER_NEST ? ncm_sphere_healpix_ring2nest : ncm_sphere_healpix_nest2ring);
    glong i, j;
    gint v;

    for (v = 0; v < 5; v++)
    {
      if (*vectors[v] != NULL)
      {
        gsl_vector_float_memcpy (temp_pix, *vectors[v]);
        for (i = 0; i < map->npix; i++)
        {
          gfloat val = gsl_vector_float_get (temp_pix, i);
          j = convert (map->nside, i);
          gsl_vector_float_set (*vectors[v], j, val);
        }
      }
    }

    gsl_vector_float_free (temp_pix);
    map->order = order;
  }

  if (init_coord)
    ncm_sphere_map_init_coord (map);
  else
    map->is_init_coord = FALSE;

  return TRUE;
}

/**
 * ncm_sphere_mapalm_new: (skip)
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmSphereMapAlm *
ncm_sphere_mapalm_new (void)
{
  NcmSphereMapAlm *mapalm = g_slice_new (NcmSphereMapAlm);
  mapalm->lmax = 0;
  mapalm->alm = NULL;
  mapalm->Nc = NULL;
  mapalm->loaded = FALSE;
  return mapalm;
}

/**
 * ncm_sphere_mapalm_init:
 * @mapalm: a #NcmSphereMapAlm
 * @lmax: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_sphere_mapalm_init (NcmSphereMapAlm *mapalm, gint lmax)
{
  gint l, m;
  gdouble last_lnpoch = 0.5 * M_LNPI - M_LN2;
  gdouble sgn = -1.0;
  mapalm->lmax = lmax;
  mapalm->alm_size = NCM_MAP_ALM_SIZE(lmax);
  mapalm->alm = gsl_vector_complex_alloc (mapalm->alm_size);
  mapalm->Nc = gsl_vector_alloc (lmax + 1);
  mapalm->sqrt_int = gsl_vector_alloc (2*lmax + 3 + 1);
  mapalm->sphPlm_recur1 = gsl_vector_alloc (mapalm->alm_size);
  mapalm->sphPlm_recur2 = gsl_vector_alloc (mapalm->alm_size);
  mapalm->sphPmm = gsl_vector_alloc (lmax + 1);
  mapalm->lnpoch_m_1_2 = gsl_vector_alloc (lmax + 1);

  for (m = 0; m <= 2*lmax+3; m++)
    gsl_vector_set (mapalm->sqrt_int, m, sqrt(m));

  for (m = 1; m <= lmax; m++)
  {
    gdouble y_mm;
    gsl_vector_set (mapalm->lnpoch_m_1_2, m, 0.5 * last_lnpoch - ncm_c_lnpi_4 ());

    y_mm = ncm_c_sqrt_1_4pi () *
      gsl_vector_get (mapalm->sqrt_int, 2*m + 1) / gsl_vector_get(mapalm->sqrt_int, m) * sgn;
    gsl_vector_set (mapalm->sphPmm, m, y_mm);

    last_lnpoch += gsl_sf_log_1plusx ( 0.5 / ((gdouble)(m)) );
    sgn *= -1.0;
  }

  for (m = 0; m <= lmax; m++)
  {
    gint start_index = NCM_MAP_ALM_M_START(lmax, m);
    for (l = m + 2; l <= lmax; l++)
    {
      gdouble rec1 =
        gsl_vector_get (mapalm->sqrt_int, l-m) *
        gsl_vector_get (mapalm->sqrt_int, 2*l+1) *
        gsl_vector_get (mapalm->sqrt_int, 2*l-1) /
        gsl_vector_get (mapalm->sqrt_int, l + m);
      gdouble rec2 =
        gsl_vector_get (mapalm->sqrt_int, l - m) *
        gsl_vector_get (mapalm->sqrt_int, l - m - 1) *
        gsl_vector_get (mapalm->sqrt_int, 2*l + 1) /
        (gsl_vector_get (mapalm->sqrt_int, l + m) *
         gsl_vector_get (mapalm->sqrt_int, l + m - 1) *
         gsl_vector_get (mapalm->sqrt_int, 2*l - 3));

      gsl_vector_set (mapalm->sphPlm_recur1, start_index + l - m,  rec1 / ((gdouble)(l-m)));
      gsl_vector_set (mapalm->sphPlm_recur2, start_index + l - m,  rec2 * ((gdouble)(l+m-1)) / ((gdouble)(l-m)));
    }
  }

  return TRUE;
}

#ifdef NUMCOSMO_HAVE_FFTW3
/**
 * ncm_sphere_mapsht_new: (skip)
 * @map: a #NcmSphereMap
 * @mapalm: a #NcmSphereMapAlm
 * @fftw_flags: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmSphereMapSHT *
ncm_sphere_mapsht_new (NcmSphereMap *map, NcmSphereMapAlm *mapalm, guint fftw_flags)
{
  NcmSphereMapSHT *mapsht = g_slice_new (NcmSphereMapSHT);
  gint i, j = 0;
  mapsht->n_rings = NCM_MAP_N_RINGS(map->nside);
  mapsht->max_ring_size = NCM_MAP_MAX_RING_SIZE(map->nside);
  mapsht->n_diff_rings = NCM_MAP_N_DIFFERENT_SIZED_RINGS(map->nside);
  mapsht->map = map;
  mapsht->mapalm = mapalm;
  mapsht->ring = gsl_vector_alloc (mapsht->max_ring_size);
  mapsht->fft_ring = gsl_matrix_complex_alloc (mapsht->n_rings, mapsht->max_ring_size);
  mapsht->sphPl0 = gsl_vector_alloc (gsl_sf_legendre_array_size (mapalm->lmax, 0));
  mapsht->sphPlm = gsl_vector_complex_alloc (mapalm->alm_size);
  mapsht->sphPlm_upper_limit = gsl_vector_alloc (mapalm->lmax + 1);

  mapsht->save_wis = TRUE;

  mapsht->forward_plans = g_slice_alloc (mapsht->n_rings * sizeof(fftw_plan));
  mapsht->backward_plans = g_slice_alloc (mapsht->n_rings * sizeof(fftw_plan));

  if (mapsht->save_wis)
  {
    ncm_cfg_init ();
    ncm_cfg_load_fftw_wisdom ("map_nside_%ld.wis", map->nside);
  }

  if (!ncm_cfg_load_vector ("Plm_upper_limit_lmax_%d.dat", mapsht->sphPlm_upper_limit, mapalm->lmax))
  {
    gsl_vector_set_all (mapsht->sphPlm_upper_limit, 1.0);
    for (i = 20; i <= mapalm->lmax; i++)
      gsl_vector_set (mapsht->sphPlm_upper_limit, i, (ncm_sphPlm_x (mapalm->lmax, i, 20)));
    ncm_cfg_save_vector ("Plm_upper_limit_lmax_%d.dat", mapsht->sphPlm_upper_limit, mapalm->lmax);
  }

  for (i = 0; i < mapsht->n_rings; i++)
  {
    gint ring_size = NCM_MAP_RING_SIZE(mapsht->map->nside, i);
    gsl_vector_complex_view fft_ring_view = gsl_matrix_complex_row (mapsht->fft_ring, i);
    gdouble theta, phi;
    ncm_sphere_healpix_pix2ang_ring (mapsht->map->nside, j, &theta, &phi);

    j += ring_size;
    mapsht->forward_plans[i] =
      fftw_plan_dft_r2c_1d (ring_size, mapsht->ring->data, (fftw_complex *)fft_ring_view.vector.data, fftw_flags);

    mapsht->backward_plans[i] =
      fftw_plan_dft_1d (ring_size, (fftw_complex *)fft_ring_view.vector.data, (fftw_complex *)fft_ring_view.vector.data,
                        FFTW_BACKWARD, fftw_flags);
    //mapsht->backward_plans[i] =
    //  fftw_plan_dft_c2r_1d (ring_size, (fftw_complex *)mapsht->fft_ring->data, mapsht->ring->data, fftw_flags);
  }

  if (mapsht->save_wis)
    ncm_cfg_save_fftw_wisdom ("map_nside_%ld.wis", map->nside);
//  if (!lmin_init)
//    ncm_cfg_save_matrix_int ("lmin_%d_%d.dat", mapsht->lmin, mapalm->lmax, map->nside);

  return mapsht;
}

/**
 * ncm_sphere_mapsht_map2alm:
 * @mapsht: a #NcmSphereMapSHT
 * @cut: FIXME
 *
 * Transform the map to alm circle by circle using fft in each one
 *
 * Returns: FIXME
 */
gboolean
ncm_sphere_mapsht_map2alm (NcmSphereMapSHT *mapsht, gdouble cut)
{
  gint i, j = 0;
  gint start_m = 0;
  gint l_size = mapsht->mapalm->lmax;
  gint lmax = mapsht->mapalm->lmax;
  gint block_size = 1024 * 20; // 1024 => 1024 * 8 * 3 = 24kb
  glong last_ring_pix = 0;
  gboolean loop_ctl = TRUE;
  gdouble four_pi_npix = 4*M_PI / mapsht->map->npix;

  g_assert (mapsht->map->order == NC_SPHERE_MAP_ORDER_RING);
  gsl_vector_complex_set_zero (mapsht->mapalm->alm);

  for (i = 0; i < mapsht->n_rings; i++)
  {
    gint ring_size = NCM_MAP_RING_SIZE (mapsht->map->nside, i);
    last_ring_pix += ring_size;
    for (; j < last_ring_pix; j++)
      gsl_vector_set (mapsht->ring, j - last_ring_pix + ring_size,
                      gsl_vector_float_get (mapsht->map->dt, j));
    fftw_execute (mapsht->forward_plans[i]);
  }

  while (loop_ctl)
  {
    gint end_m, mn = 0;
    while (((mn + l_size) < block_size) && (l_size > 0))
      mn += l_size--;
    end_m = (lmax-l_size);
    if (l_size == 0)
      loop_ctl = FALSE;
    j = 0;

    for (i = 0; i < mapsht->n_rings; i++)
    {
      gint ring_size = NCM_MAP_RING_SIZE (mapsht->map->nside, i);
      gdouble theta, phi;
      ncm_sphere_healpix_pix2ang_ring (mapsht->map->nside, j, &theta, &phi);
      j += ring_size;
      if (fabs(theta - M_PI / 2.0) >= ncm_c_degree_to_radian (cut))
        ncm_sphere_mapsht_map2alm_circle (mapsht, i, ring_size, four_pi_npix, theta, phi, start_m, end_m);
    }
    start_m = end_m + 1;
  }

  for (i = 0; i <= mapsht->mapalm->lmax; i++)
  {
    gdouble *Nc = gsl_vector_ptr (mapsht->mapalm->Nc, i);
    *Nc = (gsl_pow_2 (GSL_REAL(gsl_vector_complex_get (mapsht->mapalm->alm, NCM_MAP_ALM_INDEX(mapsht->mapalm->lmax, i, 0))))
        + gsl_pow_2 (GSL_IMAG(gsl_vector_complex_get (mapsht->mapalm->alm, NCM_MAP_ALM_INDEX(mapsht->mapalm->lmax, i, 0)))));
    for (j = 1; j <= i; j++)
    {
      *Nc += 2.0 *
        (gsl_pow_2 (GSL_REAL(gsl_vector_complex_get (mapsht->mapalm->alm, NCM_MAP_ALM_INDEX(mapsht->mapalm->lmax, i, j))))
        + gsl_pow_2 (GSL_IMAG(gsl_vector_complex_get (mapsht->mapalm->alm, NCM_MAP_ALM_INDEX(mapsht->mapalm->lmax, i, j)))));
    }
    *Nc /= (2.0 * i + 1.0);
  }

  return TRUE;
}

/**
 * ncm_sphere_mapsht_map2alm_circle:
 * @mapsht: a #NcmSphereMapSHT
 * @ring: FIXME
 * @ring_size: FIXME
 * @norma: FIXME
 * @theta: FIXME
 * @phi: FIXME
 * @start_m: FIXME
 * @end_m: FIXME
 *
 * Transform the map to alm circle by circle using fft in each one
 * Copied from gsl-1.11 specfunc/legendre_poly.c line 596
 * And then adapted...
 *
 * Returns: FIXME
 */
gboolean
ncm_sphere_mapsht_map2alm_circle (NcmSphereMapSHT *mapsht, gint ring, gint ring_size, gdouble norma, gdouble theta, gdouble phi, gint start_m, gint end_m)
{
  gint lmax = mapsht->mapalm->lmax;
  gint m;
  gdouble y_mm;
  gdouble y_mmp1;
  gdouble x = cos(theta);
  gdouble lncirc = 0.5*gsl_sf_log_1plusx(-x*x);
  gsl_vector_complex *alm = mapsht->mapalm->alm;
  gsl_complex exp_mphi;
  gsl_complex gphase;
  gsl_complex *alm_data = (gsl_complex *)&alm->data[2*NCM_MAP_ALM_M_START(lmax, start_m)];
  gsl_vector_complex_view fft_ring_view = gsl_matrix_complex_row (mapsht->fft_ring, ring);
  gsl_vector_complex *fft_ring = &fft_ring_view.vector;
  GSL_SET_COMPLEX (&exp_mphi, cos(phi), -sin(phi));
  GSL_SET_COMPLEX (&gphase, norma * cos(start_m * phi), - norma * sin (start_m * phi));

  for (m = start_m; m <= lmax && m <= end_m; m++)
  {
    gint start_index = NCM_MAP_ALM_M_START(lmax, m);
    gint fft_ring_index = m % ring_size;
    gsl_complex phase = gphase;
    if (gsl_vector_get (mapsht->sphPlm_upper_limit, m) < fabs(x))
      return TRUE;
    NCM_COMPLEX_MUL(gphase, exp_mphi);

    if (fft_ring_index <= ring_size/2)
      NCM_COMPLEX_MUL (phase, gsl_vector_complex_get (fft_ring, fft_ring_index));
    else
    {
      fft_ring_index = ring_size-fft_ring_index;
      NCM_COMPLEX_MUL_CONJUGATE (phase, gsl_vector_complex_get (fft_ring, fft_ring_index));
    }

    if(m == 0)
    {
      y_mm   = ncm_c_sqrt_1_4pi ();          /* Y00 = 1/sqrt(4pi) */
      y_mmp1 = x * ncm_c_sqrt_3_4pi ();
    }
    else
    {
      const gdouble lnpre = gsl_vector_get(mapsht->mapalm->lnpoch_m_1_2, m) + m*lncirc;
      y_mm = gsl_vector_get (mapsht->mapalm->sphPmm, m) * exp(lnpre);
      y_mmp1 = x * gsl_vector_get(mapsht->mapalm->sqrt_int, 2*m + 3) * y_mm;
    }

    if(lmax == m)
    {
      NCM_COMPLEX_INC_MUL_REAL(alm_data[0], phase, y_mm);
      alm_data = &alm_data[1];
    }
    else if(lmax == m + 1)
    {
      NCM_COMPLEX_INC_MUL_REAL(alm_data[0], phase, y_mm);
      NCM_COMPLEX_INC_MUL_REAL(alm_data[1], phase, y_mmp1);
      alm_data = &alm_data[2];
    }
    else
    {
      gdouble y_ell = 0.0;
      gint ell, rindex, lindex;
      NCM_COMPLEX_INC_MUL_REAL(alm_data[0], phase, y_mm);
      NCM_COMPLEX_INC_MUL_REAL(alm_data[1], phase, y_mmp1);

      /* Compute Y_l^m, l > m+1, upward recursion on l. */
      rindex = start_index + 2;
      lindex = 2;
      for (ell = m + 2; ell <= lmax; ell++)
      {
        const gdouble factor1 = gsl_vector_get (mapsht->mapalm->sphPlm_recur1, rindex);
        const gdouble factor2 = gsl_vector_get (mapsht->mapalm->sphPlm_recur2, rindex);
        y_ell  = (x*y_mmp1*factor1 - y_mm*factor2);
        y_mm   = y_mmp1;
        y_mmp1 = y_ell;
        NCM_COMPLEX_INC_MUL_REAL(alm_data[lindex], phase, y_ell);
        rindex++;
        lindex++;
      }
      alm_data = &alm_data[lindex];
    }
  }
  return TRUE;
}

/**
 * ncm_sphere_mapsht_alm2map:
 * @mapsht: a #NcmSphereMapSHT
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_sphere_mapsht_alm2map (NcmSphereMapSHT *mapsht)
{
  gint i, j = 0;
  glong last_ring_pix = 0;
  glong nerror = 0;

  g_assert (mapsht->map->order == NC_SPHERE_MAP_ORDER_RING);
  //gsl_vector_float_set_zero (&mapsht->map->map_view.vector);

  for (i = 0; i < mapsht->n_rings; i++)
  {
    gint plan_index = NCM_MAP_RING_PLAN_INDEX (mapsht->map->nside, i);
    gint ring_size = 4*(plan_index+1);
    gdouble theta, phi;
    gsl_vector_complex_view fft_ring_view = gsl_matrix_complex_row (mapsht->fft_ring, i);
    gsl_vector_complex *fft_ring = &fft_ring_view.vector;

    ncm_sphere_healpix_pix2ang_ring (mapsht->map->nside, j, &theta, &phi);
    last_ring_pix += ring_size;

    ncm_sphere_mapsht_alm2map_circle (mapsht, i, ring_size, theta, phi);

    fftw_execute (mapsht->backward_plans[i]);

    for (; j < last_ring_pix; j++)
    {
      gint bindex = (j - last_ring_pix + ring_size);// * (1024 / ring_size);
      gdouble error_test = 1.0 - GSL_REAL(gsl_vector_complex_get (fft_ring, bindex))/gsl_vector_float_get (mapsht->map->dt, j);
//      printf ("[%d %d] %g %g | %g %g\n", j, bindex, gsl_vector_float_get (mapsht->map->map, j),
//              GSL_REAL(gsl_vector_complex_get (fft_ring, bindex)),
//              minha, acos(minha));
//      gsl_vector_float_set (mapsht->map->dt, j, GSL_REAL(gsl_vector_complex_get (fft_ring, bindex)));
//      gsl_vector_float_set (mapsht->map->dt, j, -GSL_REAL(gsl_vector_complex_get (fft_ring, bindex))+
//                            gsl_vector_float_get (mapsht->map->dt, j));
      if (fabs(error_test) > 0.001)
        nerror++;
      gsl_vector_float_set (mapsht->map->dt, j, fabs(error_test) > 0.001 ? 4.0 : -4.0);
    }
  }
  printf ("# %ld of %ld == %1.3f%% pixels with error bigger than 100%%\n", nerror, mapsht->map->npix, nerror * 100.0 / (1.0*mapsht->map->npix));
  return TRUE;
}

/**
 * ncm_sphere_mapsht_alm2map_circle:
 * @mapsht: a #NcmSphereMapSHT
 * @ring: FIXME
 * @ring_size: FIXME
 * @theta: FIXME
 * @phi: FIXME
 *
 * Transform the map to alm circle by circle using fft in each one
 * Copied from gsl-1.11 specfunc/legendre_poly.c line 596
 * And then adapted...
 * And then adapted again...
 *
 * Returns: FIXME
 */
gboolean
ncm_sphere_mapsht_alm2map_circle (NcmSphereMapSHT *mapsht, gint ring,  gint ring_size, gdouble theta, gdouble phi)
{
  gint lmax = mapsht->mapalm->lmax;
  gint m;
  gdouble y_mm;
  gdouble y_mmp1;
  gdouble x = cos(theta);
  gdouble lncirc = 0.5*gsl_sf_log_1plusx(-x*x);
  gsl_vector_complex *alm = mapsht->mapalm->alm;
  gsl_complex *alm_data = (gsl_complex *)&alm->data[0];
  gsl_complex exp_mphi;
  gsl_complex gphase;
  gsl_vector_complex_view fft_ring_view = gsl_matrix_complex_row (mapsht->fft_ring, ring);
  gsl_vector_complex *fft_ring = &fft_ring_view.vector;

  GSL_SET_COMPLEX (&exp_mphi, cos(phi), sin(phi));
  GSL_SET_COMPLEX (&gphase, 1.0, 0.0);

  gsl_vector_complex_set_zero (fft_ring);

  for (m = 0; m <= lmax; m++)
  {
    gint start_index = NCM_MAP_ALM_M_START(lmax, m);
    gsl_complex phase = gphase;
    gsl_complex temp_sum;
    gsl_complex *fft_ring_i = gsl_vector_complex_ptr (fft_ring, m % ring_size);
    GSL_SET_COMPLEX (&temp_sum, 0.0, 0.0);
    NCM_COMPLEX_MUL(gphase, exp_mphi);

    if(m == 0)
    {
      y_mm   = ncm_c_sqrt_1_4pi ();          /* Y00 = 1/sqrt(4pi) */
      y_mmp1 = x * ncm_c_sqrt_3_4pi ();
    }
    else
    {
      const gdouble lnpre = gsl_vector_get(mapsht->mapalm->lnpoch_m_1_2, m) + m*lncirc;
      y_mm = gsl_vector_get (mapsht->mapalm->sphPmm, m) * exp(lnpre);
      y_mmp1 = x * gsl_vector_get(mapsht->mapalm->sqrt_int, 2*m + 3) * y_mm;
      GSL_REAL(phase) *= 2;
      GSL_IMAG(phase) *= 2;
    }

    if(lmax == m)
    {
      NCM_COMPLEX_INC_MUL_REAL (temp_sum, alm_data[0], y_mm);
      alm_data = &alm_data[1];
    }
    else if(lmax == m + 1)
    {
      NCM_COMPLEX_INC_MUL_REAL (temp_sum, alm_data[0], y_mm);
      NCM_COMPLEX_INC_MUL_REAL (temp_sum, alm_data[1], y_mmp1);
      alm_data = &alm_data[2];
    }
    else
    {
      gdouble y_ell = 0.0;
      gint ell, rindex, lindex;

      NCM_COMPLEX_INC_MUL_REAL (temp_sum, alm_data[0], y_mm);
      NCM_COMPLEX_INC_MUL_REAL (temp_sum, alm_data[1], y_mmp1);

      /* Compute Y_l^m, l > m+1, upward recursion on l. */
      rindex = start_index + 2;
      lindex = 2;
      for (ell = m + 2; ell <= lmax; ell++)
      {
        const gdouble factor1 = gsl_vector_get (mapsht->mapalm->sphPlm_recur1, rindex);
        const gdouble factor2 = gsl_vector_get (mapsht->mapalm->sphPlm_recur2, rindex);
        y_ell  = (x*y_mmp1*factor1 - y_mm*factor2);
        y_mm   = y_mmp1;
        y_mmp1 = y_ell;
        NCM_COMPLEX_INC_MUL_REAL (temp_sum, alm_data[lindex], y_ell);
        rindex++;
        lindex++;
      }
      alm_data = &alm_data[lindex];
    }
    NCM_COMPLEX_MUL(temp_sum, phase);
    NCM_COMPLEX_ADD(*fft_ring_i, temp_sum);
    //printf ("[%d %d] ==> %g %g\n", m, m % ring_size, GSL_REAL(temp_sum), GSL_IMAG(temp_sum));
  }
  return TRUE;
}
#endif /* NUMCOSMO_HAVE_FFTW3 */

/**
 * ncm_sphere_map_homogenize_noise:
 * @map: a #NcmSphereMap
 * @base_sigma: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_sphere_map_homogenize_noise (NcmSphereMap *map, gdouble base_sigma)
{
  gfloat min;
  gdouble sigma_m;
  glong i;
  gsl_rng *r = ncm_cfg_rng_get ();

  g_assert (map->nobs != NULL);
  min = gsl_vector_float_min (map->nobs);
  sigma_m = base_sigma / sqrt (min);

  for (i = 0; i < map->npix; i++)
  {
    gfloat nobs = gsl_vector_float_get (map->nobs, i);
    gdouble sigma = base_sigma / sqrt (nobs);
    gdouble sigma_c = sqrt (fabs(sigma_m * sigma_m - sigma * sigma));
    gfloat *pix = gsl_vector_float_ptr (map->dt, i);
    /* printf ("%f %f %f\n", sigma_m, sigma, sigma_c); */
    (*pix) += gsl_ran_gaussian_ziggurat (r, sigma_c);
  }
  return sigma_m;
}

/**
 * ncm_sphere_map_rotate_avg:
 * @map: a #NcmSphereMap
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_sphere_map_rotate_avg (NcmSphereMap *map, glong n)
{
  NcmQ *q = ncm_quaternion_new ();
  NcmTriVector v = NC_TRIVEC_NEW;
  NcmSphereMap *ramap = ncm_sphere_map_clone (map);
  glong i, j = n;
  g_assert (n > 0);

  while (j--)
  {
    ncm_quaternion_set_random (q);
    for (i = 0; i < map->npix; i++)
    {
      glong ri = i;
      ncm_sphere_healpix_pix2vec_ring (map->nside, i, v);
      ncm_quaternion_rotate (q, v);
      ncm_sphere_healpix_vec2pix_ring (map->nside, v, &ri);
      gsl_vector_float_set (ramap->dt, i,
                            gsl_vector_float_get (ramap->dt, i) +
                            gsl_vector_float_get (map->dt, ri));
    }
    printf("# nhuc\n");fflush(stdout);
  }

  for (i = 0; i < map->npix; i++)
    gsl_vector_float_set (map->dt, i, gsl_vector_float_get (ramap->dt, i) / (1.0 * n + 1.0));

  //ncm_sphere_map_free (ramap);

  return TRUE;
}
