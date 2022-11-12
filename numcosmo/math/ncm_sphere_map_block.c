/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_sphere_map_block.c
 *
 *  Fri January 19 13:12:06 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sphere_map_block.c
 *
 * Copyright (C) 2018 - Sandro Dias Pinto Vitenti
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

static void 
NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_run_over_l) (NcmSphereMap *pix, NcmSFSphericalHarmonicsY *sphaY, NcmComplex * restrict alm, const NcmComplex *Fim_i)
{
  NcmSphereMapPrivate * const self = pix->priv;
  const gint lmax        = self->lmax;
  const gint lmaxm1      = self->lmax - 1;
	const gint lmaxmstepm2 = self->lmax - NCM_SPHERE_MAP_BLOCK_STEPM2;
  const gint m           = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
  gint l                 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
  gdouble Ylm[NCM_SPHERE_MAP_BLOCK_STEP];

  alm += (l - m);
  
  while (l < lmaxmstepm2) 
  {
    gint i;
    
    ncm_sf_spherical_harmonics_Y_next_l2pn (sphaY, Ylm, NCM_SPHERE_MAP_BLOCK_STEPM2);

    for (i = 0; i < NCM_SPHERE_MAP_BLOCK_STEP; i++)
    {
      ncm_complex_res_add_mul_real (&alm[i], Fim_i, Ylm[i]);
    }

    alm += NCM_SPHERE_MAP_BLOCK_STEP;
    l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
  }
  
  while (l < lmax) 
  {
    gint i;
    ncm_sf_spherical_harmonics_Y_next_l2 (sphaY, Ylm);

    for (i = 0; i < 2; i++)
    {
      ncm_complex_res_add_mul_real (&alm[i], Fim_i, Ylm[i]);
    }

    alm += 2;
    l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
  }
  
	if (l == lmaxm1)
	{
    const gdouble Ylm = ncm_sf_spherical_harmonics_Y_get_lp1m (sphaY);
    ncm_complex_res_add_mul_real (&alm[0], Fim_i, Ylm);
  }  
  else if (l == lmax)
	{
    const gdouble Ylm = ncm_sf_spherical_harmonics_Y_get_lm (sphaY);
    ncm_complex_res_add_mul_real (&alm[0], Fim_i, Ylm);
  }  
}

static void 
NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_run_over_l_array) (NcmSphereMap *pix, NcmSFSphericalHarmonicsYArray *sphaYa, NcmComplex * restrict alm, const NcmComplex * restrict Fim_i)
{
  NcmSphereMapPrivate * const self = pix->priv;
  const gint lmax        = self->lmax;
  const gint lmaxm1      = self->lmax - 1;
	const gint lmaxmstepm2 = self->lmax - NCM_SPHERE_MAP_BLOCK_STEPM2;
  const gint m           = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
  gint l                 = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
	gdouble Ylm[NCM_SPHERE_MAP_BLOCK_STEP * NCM_SPHERE_MAP_BLOCK_NCT];
	/*gdouble *Ylm = fftw_alloc_real (NCM_SPHERE_MAP_BLOCK_STEP * NCM_SPHERE_MAP_BLOCK_NCT);*/ 
  
  alm += (l - m);

  while (l < lmaxmstepm2) 
  {
    gint i, j;
    ncm_sf_spherical_harmonics_Y_array_next_l2pn (sphaYa, NCM_SPHERE_MAP_BLOCK_NCT, Ylm, NCM_SPHERE_MAP_BLOCK_STEPM2);

    for (i = 0; i < NCM_SPHERE_MAP_BLOCK_STEP; i++)
    {
      for (j = 0; j < NCM_SPHERE_MAP_BLOCK_NCT; j++)
      {
        const gdouble Ylm_bj = Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (j, i, NCM_SPHERE_MAP_BLOCK_NCT)];
        ncm_complex_res_add_mul_real (&alm[i], &Fim_i[j], Ylm_bj);
      }
    }
    
    alm += NCM_SPHERE_MAP_BLOCK_STEP;
    l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
  }
  
  while (l < lmax) 
  {
    gint i, j;
    ncm_sf_spherical_harmonics_Y_array_next_l2 (sphaYa, NCM_SPHERE_MAP_BLOCK_NCT, Ylm);

    for (i = 0; i < 2; i++)
    {
      for (j = 0; j < NCM_SPHERE_MAP_BLOCK_NCT; j++)
      {
        const gdouble Ylm_bj = Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (j, i, NCM_SPHERE_MAP_BLOCK_NCT)];
        ncm_complex_res_add_mul_real (&alm[i], &Fim_i[j], Ylm_bj);
      }
    }
    
    alm += 2;
    l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
  }
  
	if (l == lmaxm1)
	{
    gint j;    
    for (j = 0; j < NCM_SPHERE_MAP_BLOCK_NCT; j++)
    {
      const gdouble Ylm = ncm_sf_spherical_harmonics_Y_array_get_lp1m (sphaYa, NCM_SPHERE_MAP_BLOCK_NCT, j);
      ncm_complex_res_add_mul_real (&alm[0], &Fim_i[j], Ylm);
    }
  }  
  else if (l == lmax)
	{
    gint j;    
    for (j = 0; j < NCM_SPHERE_MAP_BLOCK_NCT; j++)
    {
      const gdouble Ylm = ncm_sf_spherical_harmonics_Y_array_get_lm (sphaYa, NCM_SPHERE_MAP_BLOCK_NCT, j);
      ncm_complex_res_add_mul_real (&alm[0], &Fim_i[j], Ylm);      
    }
  }  

	/*fftw_free (Ylm);*/
}

static void
NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_prepare_circle) (NcmSphereMap *pix, const gint64 r_i, const gint64 i, _fft_complex * restrict * restrict Fima, gdouble *theta, gdouble *phi, gint64 *ring_size, gint64 *ring_size_2)
{
  NcmSphereMapPrivate * const self = pix->priv;
  const gint64 ring_fi = ncm_sphere_map_get_ring_first_index (pix, r_i);

  ring_size[i]   = ncm_sphere_map_get_ring_size (pix, r_i);
  ring_size_2[i] = ring_size[i] / 2;

  Fima[i] = &((_fft_complex *)self->fft_pvec)[ring_fi];

  ncm_sphere_map_pix2ang_ring (pix, ring_fi, &theta[i], &phi[i]);
}

static gboolean
NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_get_alm_from_apcircles) (NcmSphereMap *pix, NcmSFSphericalHarmonicsYArray *sphaYa, NcmComplex * restrict alm, const gint64 r_ini, const gint64 nrings, const gint64 mmax)
{
  NcmSphereMapPrivate * const self = pix->priv;
  const gdouble pix_area                   = 4.0 * M_PI / self->npix;
  const NcmSphereMapBlock * restrict block = &g_array_index (self->block_data, NcmSphereMapBlock, r_ini / NCM_SPHERE_MAP_BLOCK_NC);
  
  complex double Fim_i[NCM_SPHERE_MAP_BLOCK_NCT];
  complex double expIphi[NCM_SPHERE_MAP_BLOCK_NCT];
  complex double phase[NCM_SPHERE_MAP_BLOCK_NCT];
  
  gint m, l0;
  gint64 i;

  m  = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
  l0 = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
  
  for (i = 0; i < NCM_SPHERE_MAP_BLOCK_NC; i++)
  {
    const gint64 api   = NCM_SPHERE_MAP_BLOCK_NCT - i - 1;
    
    expIphi[i]   = cexp (-I * block->phi[i]);
    expIphi[api] = cexp (-I * block->phi[api]);
    phase[i]     = pix_area * cexp (-I * block->phi[i] * m);
    phase[api]   = pix_area * cexp (-I * block->phi[api] * m);
  }

  while (TRUE)
  {
    for (i = 0; i < NCM_SPHERE_MAP_BLOCK_NCT; i++)
    {
      const gint fft_ring_index = m % block->ring_size[i];

      if (fft_ring_index <= block->ring_size_2[i])
      {
        Fim_i[i] = block->Fima[i][fft_ring_index] * phase[i];
      }
      else
      {
        Fim_i[i] = conj (block->Fima[i][block->ring_size[i] - fft_ring_index]) * phase[i];
      }
      phase[i] *= expIphi[i];
    }
    
    NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_run_over_l_array) (pix, sphaYa, alm, NCM_COMPLEX (Fim_i));

    if (m < mmax)
    {
      ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, NCM_SPHERE_MAP_BLOCK_NCT);
      l0 = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
      
      if (l0 > self->lmax)
      {
        return FALSE;
        break;
      }
      else
      {
        alm += self->lmax - m + 1;
        m    = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
      }
    }
    else
    {
      if (m < self->lmax)
      {
        ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, NCM_SPHERE_MAP_BLOCK_NCT);
        return TRUE;
      }
      else
        return FALSE;
      break;
    }
  }
}

static void
NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_get_alm_from_circle) (NcmSphereMap *pix, NcmSFSphericalHarmonicsY *sphaY, NcmComplex * restrict alm, gint64 r_i)
{
  NcmSphereMapPrivate * const self = pix->priv;
  const gint64 ring_size      = ncm_sphere_map_get_ring_size (pix, r_i);
  const gint64 ring_size_2    = ring_size / 2;
  const gint64 ring_fi        = ncm_sphere_map_get_ring_first_index (pix, r_i);
  _fft_complex * restrict Fim = &((_fft_complex *)self->fft_pvec)[ring_fi];
  const gdouble pix_area      = 4.0 * M_PI / self->npix;
  gdouble theta_i             = 0.0;
  gdouble phi                 = 0.0; 
  gint m;

  ncm_sphere_map_pix2ang_ring (pix, ring_fi, &theta_i, &phi);

  if (phi != 0.0)
  {
    const complex double expIphi = cexp (-I * phi);
    complex double phase         = pix_area;
    gint l0;

    ncm_sf_spherical_harmonics_start_rec (self->spha, sphaY, theta_i);

    m  = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
    l0 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);

    while (TRUE)
		{
      const gint fft_ring_index = m % ring_size;
      complex double ring_m     = ((fft_ring_index <= ring_size_2) ? Fim[fft_ring_index] : conj (Fim[ring_size - fft_ring_index])) * phase;

      phase *= expIphi;
      NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_run_over_l) (pix, sphaY, alm, NCM_COMPLEX (&ring_m));

      if (m < self->lmax)
      {
        ncm_sf_spherical_harmonics_Y_next_m (sphaY);
        l0 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);

        if (l0 > self->lmax)
          break;
        else
        {
          alm += self->lmax - m + 1;      
          m    = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
        }
      }
      else
        break;
    }
  }
  else
  {
    gint l0;
    ncm_sf_spherical_harmonics_start_rec (self->spha, sphaY, theta_i);

    m  = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
    l0 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
    
    while (TRUE)
		{
      const gint fft_ring_index = m % ring_size;
      complex double ring_m     = ((fft_ring_index <= ring_size_2) ? Fim[fft_ring_index] : conj (Fim[ring_size - fft_ring_index])) * pix_area;

      NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_run_over_l) (pix, sphaY, alm, NCM_COMPLEX (&ring_m));

      if (m < self->lmax)
      {
        ncm_sf_spherical_harmonics_Y_next_m (sphaY);
        l0 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
        
        if (l0 > self->lmax)
          break;
        else
        {
          alm += self->lmax - m + 1;
          m    = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
        }
      }
      else
        break;
    }
  }
}

static void 
NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_map2alm_run) (NcmSphereMap *pix)
{
  NcmSphereMapPrivate * const self = pix->priv;
  NcmSFSphericalHarmonicsY *sphaY  = ncm_sf_spherical_harmonics_Y_new (self->spha, NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL);
  NcmComplex * restrict alm        = NCM_COMPLEX (self->alm);
  const gint64 nrings              = ncm_sphere_map_get_nrings (pix);
  const gint64 lr_i                = self->block_ring_size;
  const gint64 nrleft              = self->last_sing_ring;
  gint64 r_i, m, offset;
  
  memset (alm, 0, sizeof (NcmComplex) * self->alm_len);

  for (r_i = 0; r_i < lr_i; r_i += NCM_SPHERE_MAP_BLOCK_NC)
  {
    const gint64 i                        = r_i / NCM_SPHERE_MAP_BLOCK_NC;
    NcmSFSphericalHarmonicsYArray *sphaYa = g_ptr_array_index (self->sphaYa_array, i);
    ncm_sf_spherical_harmonics_Y_array_reset (sphaYa, NCM_SPHERE_MAP_BLOCK_NCT);
  }

  m      = 0;
  offset = 0;
  
  while (m < self->lmax)
  {
    gboolean c_all = FALSE;
    gint64 chunk = 0;
    do {
      chunk += (self->lmax + 1 - m);
      m++;
    } while ((chunk < (1024 * NCM_SPHERE_MAP_BLOCK_CM)) && (m < self->lmax + 1));
    
    /*printf ("# mmax %ld lmax %u chunk %ld offset %ld\n", m - 1, self->lmax, chunk, offset);*/

    for (r_i = 0; r_i < lr_i; r_i += NCM_SPHERE_MAP_BLOCK_NC)
    {
      const gint64 i                        = r_i / NCM_SPHERE_MAP_BLOCK_NC;
      NcmSFSphericalHarmonicsYArray *sphaYa = g_ptr_array_index (self->sphaYa_array, i);
      const gboolean c                      = NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_get_alm_from_apcircles) (pix, sphaYa, alm + offset, r_i, nrings, m - 1);
      c_all = c_all || c; 
    }
    
    if (!c_all)
      break;

    offset += chunk;
  }

  for (; r_i < nrleft; r_i++)
  {
    NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_get_alm_from_circle) (pix, sphaY, alm, r_i);
  }
}

static void
NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_D_m) (NcmSphereMap *pix, NcmSFSphericalHarmonicsY *sphaY, NcmComplex * restrict * restrict alm_ptr, const gint ring_size, const gint m, NcmComplex *D_m)
{
  NcmSphereMapPrivate * const self = pix->priv;
	const gint lmaxmstepm2 = self->lmax - NCM_SPHERE_MAP_BLOCK_INV_STEPM2;
  gdouble Ylm[NCM_SPHERE_MAP_BLOCK_INV_STEP];
  gint l;

  l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
  alm_ptr[0] += (l - m);

  while (l < lmaxmstepm2)
  {
    gint i;
    ncm_sf_spherical_harmonics_Y_next_l2pn (sphaY, Ylm, NCM_SPHERE_MAP_BLOCK_INV_STEPM2);

    for (i = 0; i < NCM_SPHERE_MAP_BLOCK_INV_STEP; i++)
    {
      ncm_complex_res_add_mul_real (D_m, &alm_ptr[0][i], Ylm[i]);
    }

    alm_ptr[0] += NCM_SPHERE_MAP_BLOCK_INV_STEP;
    l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
  }

  if (l <= self->lmax)
  {
    while (TRUE)
    {
      const gdouble Ylm = ncm_sf_spherical_harmonics_Y_get_lm (sphaY);

      ncm_complex_res_add_mul_real (D_m, &alm_ptr[0][0], Ylm);

      alm_ptr[0]++;

      if (l < self->lmax)
      {
        ncm_sf_spherical_harmonics_Y_next_l (sphaY);
        l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
      }
      else
        break;      
    }
  }
  
  return;
}

static void
NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_D_m_array) (NcmSphereMap *pix, NcmSFSphericalHarmonicsYArray *sphaYa, NcmComplex * restrict * restrict alm_ptr, const gint m, NcmComplex *D_m)
{
  NcmSphereMapPrivate * const self = pix->priv;
	const gint lmaxmstepm2 = self->lmax - NCM_SPHERE_MAP_BLOCK_INV_STEPM2;
  gdouble Ylm[NCM_SPHERE_MAP_BLOCK_INV_STEP * NCM_SPHERE_MAP_BLOCK_INV_NCT];
  gint l;

  l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
  alm_ptr[0] += (l - m);

  while (l < lmaxmstepm2)
  {
    gint i, j;
    ncm_sf_spherical_harmonics_Y_array_next_l2pn (sphaYa, NCM_SPHERE_MAP_BLOCK_INV_NCT, Ylm, NCM_SPHERE_MAP_BLOCK_INV_STEPM2);

		for (i = 0; i < NCM_SPHERE_MAP_BLOCK_INV_STEP; i++)
		{
			for (j = 0; j < NCM_SPHERE_MAP_BLOCK_INV_NCT; j++)
			{
				const gdouble Ylm_bj = Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (j, i, NCM_SPHERE_MAP_BLOCK_INV_NCT)];
        ncm_complex_res_add_mul_real (&D_m[j], &alm_ptr[0][i], Ylm_bj);
			}
		}
		
    alm_ptr[0] += NCM_SPHERE_MAP_BLOCK_INV_STEP;
    l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
  }

  if (l <= self->lmax)
  {
    while (TRUE)
    {
			gint j;
			for (j = 0; j < NCM_SPHERE_MAP_BLOCK_INV_NCT; j++)
			{
				const gdouble Ylm = ncm_sf_spherical_harmonics_Y_array_get_lm (sphaYa, NCM_SPHERE_MAP_BLOCK_INV_NCT, j);
        ncm_complex_res_add_mul_real (&D_m[j], &alm_ptr[0][0], Ylm);
			}
			alm_ptr[0]++;

      if (l < self->lmax)
      {
        ncm_sf_spherical_harmonics_Y_array_next_l (sphaYa, NCM_SPHERE_MAP_BLOCK_INV_NCT);
        l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
      }
      else
        break;      
    }
  }
  
  return;
}

static void
NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_circle_from_alm) (NcmSphereMap *pix, gint64 r_i)
{
  NcmSphereMapPrivate * const self = pix->priv;
  NcmSFSphericalHarmonicsY *sphaY = ncm_sf_spherical_harmonics_Y_new (self->spha, NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL);
  const gint ring_size    = ncm_sphere_map_get_ring_size (pix, r_i);
  const gint ring_fi      = ncm_sphere_map_get_ring_first_index (pix, r_i);
  _fft_complex *Fim       = &((_fft_complex *)self->fft_pvec)[ring_fi];
  gdouble theta_i         = 0.0, phi = 0.0;
  gint ring_size_2        = ring_size / 2;
  complex double *alm_ptr = self->alm;
  complex double expIphi, phase;
  gint m, l0;

  ncm_sphere_map_pix2ang_ring (pix, ring_fi, &theta_i, &phi);

  expIphi = cexp (I * phi);
  phase   = expIphi;
  
  ncm_sf_spherical_harmonics_start_rec (self->spha, sphaY, theta_i);

  memset (Fim, 0, sizeof (_fft_complex) * ring_size);

  m  = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
  l0 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
  {
    complex double D_m = 0.0;
    NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_D_m) (pix, sphaY, NCM_COMPLEX_PTR (&alm_ptr), ring_size, m, NCM_COMPLEX (&D_m));
    Fim[m] += creal (D_m);
  }

  ncm_sf_spherical_harmonics_Y_next_m (sphaY);
  m  = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
  l0 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);

  while (TRUE)
  {
    const gint fft_ring_index = m % ring_size;
    complex double D_m        = 0.0;

		NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_D_m) (pix, sphaY, NCM_COMPLEX_PTR (&alm_ptr), ring_size, m, NCM_COMPLEX (&D_m));

    D_m   *= phase;
    phase *= expIphi;

    if ((fft_ring_index == 0) || (fft_ring_index == ring_size_2))
    {
      Fim[fft_ring_index] += 2.0 * creal (D_m);
    }
    else if (fft_ring_index < ring_size_2)
    {
      Fim[fft_ring_index] += D_m;
    }
    else if (fft_ring_index > ring_size_2)
    {
      Fim[ring_size - fft_ring_index] += conj (D_m);
    }

    if (m < self->lmax)
    {
      ncm_sf_spherical_harmonics_Y_next_m (sphaY);
      m  = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
      l0 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
      if (l0 > self->lmax)
        break;
    }
    else
      break;
  }
  
  ncm_sf_spherical_harmonics_Y_free (sphaY);
}

static void
NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_apcircles_from_alm) (NcmSphereMap *pix, const gint64 r_ini, const gint64 nrings)
{
  NcmSphereMapPrivate * const self = pix->priv;
  NcmSFSphericalHarmonicsYArray *sphaYa = ncm_sf_spherical_harmonics_Y_array_new (self->spha, NCM_SPHERE_MAP_BLOCK_INV_NCT, NCM_SF_SPHERICAL_HARMONICS_ARRAY_DEFAULT_ABSTOL);
  complex double *alm_ptr               = self->alm;

	_fft_complex * restrict Fima[NCM_SPHERE_MAP_BLOCK_INV_NCT];
  gint64 ring_size[NCM_SPHERE_MAP_BLOCK_INV_NCT];
  gint64 ring_size_2[NCM_SPHERE_MAP_BLOCK_INV_NCT];
  gdouble theta[NCM_SPHERE_MAP_BLOCK_INV_NCT];
  gdouble phi[NCM_SPHERE_MAP_BLOCK_INV_NCT];
	complex double D_m[NCM_SPHERE_MAP_BLOCK_INV_NCT];
  complex double expIphi[NCM_SPHERE_MAP_BLOCK_INV_NCT];
  complex double phase[NCM_SPHERE_MAP_BLOCK_INV_NCT];

	gint m, l0, i;

  for (i = 0; i < NCM_SPHERE_MAP_BLOCK_INV_NC; i++)
  {
    const gint64 r_i   = r_ini + i;
    const gint64 apr_i = nrings - r_i - 1;
    const gint64 api   = NCM_SPHERE_MAP_BLOCK_INV_NCT - i - 1;
    
    NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_prepare_circle) (pix, r_i,   i,   Fima, theta, phi, ring_size, ring_size_2);
    NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_prepare_circle) (pix, apr_i, api, Fima, theta, phi, ring_size, ring_size_2);

		{	
		  expIphi[i]   = cexp (I * phi[i]);
		  expIphi[api] = cexp (I * phi[api]);
		  phase[i]     = expIphi[i];
		  phase[api]   = expIphi[api];
	  }

		memset (Fima[i],   0, sizeof (_fft_complex) * ring_size[i]);
		memset (Fima[api], 0, sizeof (_fft_complex) * ring_size[api]);
	}

  ncm_sf_spherical_harmonics_start_rec_array (self->spha, sphaYa, NCM_SPHERE_MAP_BLOCK_INV_NCT, theta);

  m  = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
  l0 = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);

	memset (D_m, 0, sizeof (gsl_complex) * NCM_SPHERE_MAP_BLOCK_INV_NCT);

	NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_D_m_array) (pix, sphaYa, NCM_COMPLEX_PTR (&alm_ptr), m, NCM_COMPLEX (&D_m));

	for (i = 0; i < NCM_SPHERE_MAP_BLOCK_INV_NCT; i++)
  {
		Fima[i][m] += creal (D_m[i]);
  }

  ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, NCM_SPHERE_MAP_BLOCK_INV_NCT);
  m  = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
  l0 = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);

  while (TRUE)
  {
		memset (D_m, 0, sizeof (gsl_complex) * NCM_SPHERE_MAP_BLOCK_INV_NCT);
    NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_D_m_array) (pix, sphaYa, NCM_COMPLEX_PTR (&alm_ptr), m, NCM_COMPLEX (&D_m));

		for (i = 0; i < NCM_SPHERE_MAP_BLOCK_INV_NCT; i++)
		{
			const gint fft_ring_index = m % ring_size[i];

      D_m[i]   *= phase[i];
      phase[i] *= expIphi[i];

			if ((fft_ring_index == 0) || (fft_ring_index == ring_size_2[i]))
			{
				Fima[i][fft_ring_index] += 2.0 * creal (D_m[i]);
			}
			else if (fft_ring_index < ring_size_2[i])
			{
				Fima[i][fft_ring_index] += D_m[i];
			}
			else if (fft_ring_index > ring_size_2[i])
			{
				Fima[i][ring_size[i] - fft_ring_index] += conj (D_m[i]);
			}
		}
		if (m < self->lmax)
    {
      ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, NCM_SPHERE_MAP_BLOCK_INV_NCT);
      m  = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
      l0 = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
      if (l0 > self->lmax)
        break;
    }
    else
      break;
  }
  
  ncm_sf_spherical_harmonics_Y_array_free (sphaYa);
}

static void 
NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_alm2map_run) (NcmSphereMap *pix)
{
  const gint64 nrings   = ncm_sphere_map_get_nrings (pix);
  const gint64 nrings_2 = nrings / 2;
  gint64 r_i, nrleft;

  r_i = 0;
  while (TRUE)
  {
    if (r_i + NCM_SPHERE_MAP_BLOCK_INV_NC < nrings_2)
    {
      NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_apcircles_from_alm) (pix, r_i, nrings);
      r_i += NCM_SPHERE_MAP_BLOCK_INV_NC;
    }
    else
      break;
  }

  nrleft = nrings - r_i;
  for (; r_i < nrleft; r_i++)
  {
    NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_circle_from_alm) (pix, r_i);
  }
}
