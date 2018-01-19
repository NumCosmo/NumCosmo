/***************************************************************************
 *            ncm_sphere_map_block.c
 *
 *  Fri January 19 13:12:06 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
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
NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_run_over_l) (NcmSphereMap *pix, NcmSFSphericalHarmonicsY *sphaY, const gint si, const gint l0, const gint m, const complex double Fim_i)
{
	const gdouble Re_Fim_i                     = creal (Fim_i);
	const gdouble Im_Fim_i                     = cimag (Fim_i);
  const gint lmax                            = pix->lmax;
  const gint lmaxm1                          = pix->lmax - 1;
	const gint lmaxmstepm2                     = pix->lmax - NCM_SPHERE_MAP_BLOCK_STEPM2;
	gsl_complex * restrict alm                 = &g_array_index (pix->alm, gsl_complex, si + l0);
  gdouble Ylm[NCM_SPHERE_MAP_BLOCK_STEP];
  
  gint l = l0;

  while (l < lmaxmstepm2) 
  {
    gint i;
    
    ncm_sf_spherical_harmonics_Y_next_l2pn (sphaY, Ylm, NCM_SPHERE_MAP_BLOCK_STEPM2);

    for (i = 0; i < NCM_SPHERE_MAP_BLOCK_STEP; i++)
    {
      const gdouble Re_alm = Ylm[i] * Re_Fim_i;
      const gdouble Im_alm = Ylm[i] * Im_Fim_i;
      GSL_REAL (alm[i]) += Re_alm;
      GSL_IMAG (alm[i]) += Im_alm;
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
      const gdouble Re_alm = Ylm[i] * Re_Fim_i;
      const gdouble Im_alm = Ylm[i] * Im_Fim_i;
      GSL_REAL (alm[i]) += Re_alm;
      GSL_IMAG (alm[i]) += Im_alm;
    }

    alm += 2;
    l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
  }
  
	if (l == lmaxm1)
	{
    const gdouble Ylm    = ncm_sf_spherical_harmonics_Y_get_lp1m (sphaY);
    const gdouble Re_alm = Ylm * Re_Fim_i;
    const gdouble Im_alm = Ylm * Im_Fim_i;

    GSL_REAL (alm[0]) += Re_alm;
    GSL_IMAG (alm[0]) += Im_alm;
  }  
  else if (l == lmax)
	{
    const gdouble Ylm    = ncm_sf_spherical_harmonics_Y_get_lm (sphaY);
    const gdouble Re_alm = Ylm * Re_Fim_i;
    const gdouble Im_alm = Ylm * Im_Fim_i;

    GSL_REAL (alm[0]) += Re_alm;
    GSL_IMAG (alm[0]) += Im_alm;
  }  
}

static void 
NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_run_over_l_array) (NcmSphereMap *pix, NcmSFSphericalHarmonicsYArray *sphaYa, const gint si, const gint l0, const gint m, const gdouble * restrict Re_Fim_i, const gdouble * restrict Im_Fim_i)
{
  const gint lmax            = pix->lmax;
  const gint lmaxm1          = pix->lmax - 1;
	const gint lmaxmstepm2     = pix->lmax - NCM_SPHERE_MAP_BLOCK_STEPM2;
	gsl_complex * restrict alm = &g_array_index (pix->alm, gsl_complex, si + l0);
  gdouble Ylm[NCM_SPHERE_MAP_BLOCK_STEP * NCM_SPHERE_MAP_BLOCK_NCT];
  gint l = l0;
  
  while (l < lmaxmstepm2) 
  {
    gint i, j;
    ncm_sf_spherical_harmonics_Y_array_next_l2pn (sphaYa, NCM_SPHERE_MAP_BLOCK_NCT, Ylm, NCM_SPHERE_MAP_BLOCK_STEPM2);

    for (i = 0; i < NCM_SPHERE_MAP_BLOCK_STEP; i++)
    {
      for (j = 0; j < NCM_SPHERE_MAP_BLOCK_NCT; j++)
      {
        const gdouble Ylm_bj = Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (j, i, NCM_SPHERE_MAP_BLOCK_NCT, NCM_SPHERE_MAP_BLOCK_STEP)];
        const gdouble Re_alm = Ylm_bj * Re_Fim_i[j];
        const gdouble Im_alm = Ylm_bj * Im_Fim_i[j];

        GSL_REAL (alm[i]) += Re_alm;
        GSL_IMAG (alm[i]) += Im_alm;
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
        const gdouble Ylm_bj = Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (j, i, NCM_SPHERE_MAP_BLOCK_NCT, 2)];
        const gdouble Re_alm = Ylm_bj * Re_Fim_i[j];
        const gdouble Im_alm = Ylm_bj * Im_Fim_i[j];

        GSL_REAL (alm[i]) += Re_alm;
        GSL_IMAG (alm[i]) += Im_alm;
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
      const gdouble Ylm    = ncm_sf_spherical_harmonics_Y_array_get_lp1m (sphaYa, j);
      const gdouble Re_alm = Ylm * Re_Fim_i[j];
      const gdouble Im_alm = Ylm * Im_Fim_i[j];

      GSL_REAL (alm[0]) += Re_alm;
      GSL_IMAG (alm[0]) += Im_alm;
    }
  }  
  else if (l == lmax)
	{
    gint j;    
    for (j = 0; j < NCM_SPHERE_MAP_BLOCK_NCT; j++)
    {
      const gdouble Ylm    = ncm_sf_spherical_harmonics_Y_array_get_lm (sphaYa, j);
      const gdouble Re_alm = Ylm * Re_Fim_i[j];
      const gdouble Im_alm = Ylm * Im_Fim_i[j];

      GSL_REAL (alm[0]) += Re_alm;
      GSL_IMAG (alm[0]) += Im_alm;
    }
  }  
}

static void
NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_prepare_circle) (NcmSphereMap *pix, const gint64 r_i, const gint64 i, _fft_complex * restrict * restrict Fima, gdouble *theta, gdouble *phi, gint64 *ring_size, gint64 *ring_size_2)
{
  const gint64 ring_fi = ncm_sphere_map_get_ring_first_index (pix, r_i);

  ring_size[i]   = ncm_sphere_map_get_ring_size (pix, r_i);
  ring_size_2[i] = ring_size[i] / 2;

  Fima[i] = &((_fft_complex *)pix->fft_pvec)[ring_fi];

  ncm_sphere_map_pix2ang_ring (pix, ring_fi, &theta[i], &phi[i]);
}

static void
NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_get_alm_from_apcircles) (NcmSphereMap *pix, const gint64 r_ini, const gint64 nrings)
{
  const gint64 nc                       = NCM_SPHERE_MAP_BLOCK_NC;
  const gint64 nct                      = NCM_SPHERE_MAP_BLOCK_NCT;
  const gdouble pix_area                = 4.0 * M_PI / pix->npix;
  NcmSFSphericalHarmonicsYArray *sphaYa = ncm_sf_spherical_harmonics_Y_array_new (pix->spha, nct, NCM_SF_SPHERICAL_HARMONICS_ARRAY_DEFAULT_ABSTOL);

  _fft_complex * restrict Fima[NCM_SPHERE_MAP_BLOCK_NCT];
  gint64 ring_size[NCM_SPHERE_MAP_BLOCK_NCT];
  gint64 ring_size_2[NCM_SPHERE_MAP_BLOCK_NCT];
  gdouble theta[NCM_SPHERE_MAP_BLOCK_NCT];
  gdouble phi[NCM_SPHERE_MAP_BLOCK_NCT];
  gdouble Re_Fim_i[NCM_SPHERE_MAP_BLOCK_NCT];
  gdouble Im_Fim_i[NCM_SPHERE_MAP_BLOCK_NCT];
  complex double expIphi_0[NCM_SPHERE_MAP_BLOCK_NCT];
  complex double phase[NCM_SPHERE_MAP_BLOCK_NCT];
  
  gint m, si, l0;
  gint64 i;

  for (i = 0; i < nc; i++)
  {
    const gint64 r_i   = r_ini + i;
    const gint64 apr_i = nrings - r_i - 1;
    const gint64 api   = nct - i - 1;
    
    NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_prepare_circle) (pix, r_i,   i,   Fima, theta, phi, ring_size, ring_size_2);
    NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_prepare_circle) (pix, apr_i, api, Fima, theta, phi, ring_size, ring_size_2);

    expIphi_0[i]   = cexp (-I * phi[i]);
    expIphi_0[api] = cexp (-I * phi[api]);
    phase[i]       = pix_area;
    phase[api]     = pix_area;
  }

  ncm_sf_spherical_harmonics_start_rec_array (pix->spha, sphaYa, nct, theta);

  m  = 0;
  si = 0;
  l0 = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
  
  while (TRUE)
  {
    for (i = 0; i < nct; i++)
    {
      const gint fft_ring_index = m % ring_size[i];

      if (fft_ring_index <= ring_size_2[i])
      {
        const _fft_complex cFim_i = Fima[i][fft_ring_index] * phase[i];
        Re_Fim_i[i] = creal (cFim_i);
        Im_Fim_i[i] = cimag (cFim_i);
      }
      else
      {
        const _fft_complex cFim_i = conj (Fima[i][ring_size[i] - fft_ring_index]) * phase[i];
        Re_Fim_i[i] = creal (cFim_i);
        Im_Fim_i[i] = cimag (cFim_i);
      }
      phase[i] *= expIphi_0[i];
    }
    
    NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_run_over_l_array) (pix, sphaYa, si, l0, m, Re_Fim_i, Im_Fim_i);

    if (m < pix->lmax)
    {
      ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, nct);
      l0 = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
      
      if (l0 > pix->lmax)
        break;
    }
    else
      break;

    m++;
    si += pix->lmax - m + 1;      
  }
  
  ncm_sf_spherical_harmonics_Y_array_free (sphaYa);
}

static void
NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_get_alm_from_circle) (NcmSphereMap *pix, gint64 r_i)
{
  NcmSFSphericalHarmonicsY *sphaY = ncm_sf_spherical_harmonics_Y_new (pix->spha, NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL);
  const gint64 ring_size          = ncm_sphere_map_get_ring_size (pix, r_i);
  const gint64 ring_size_2        = ring_size / 2;
  const gint64 ring_fi            = ncm_sphere_map_get_ring_first_index (pix, r_i);
  _fft_complex * restrict Fim     = &((_fft_complex *)pix->fft_pvec)[ring_fi];
  const gdouble pix_area          = 4.0 * M_PI / pix->npix;
  gdouble theta_i                 = 0.0;
  gdouble phi_0                   = 0.0; 
  gint m, si;

  ncm_sphere_map_pix2ang_ring (pix, ring_fi, &theta_i, &phi_0);

  if (phi_0 != 0.0)
  {
    const complex double expIphi_0 = cexp (-I * phi_0);
    complex double phase           = pix_area;
    gint l0;

    ncm_sf_spherical_harmonics_start_rec (pix->spha, sphaY, theta_i);

    m  = 0;
    si = 0;

    l0 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
    while (TRUE)
		{
      const gint fft_ring_index = m % ring_size;
      complex double ring_m     = ((fft_ring_index <= ring_size_2) ? Fim[fft_ring_index] : conj (Fim[ring_size - fft_ring_index])) * phase;

      phase *= expIphi_0;
      NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_run_over_l) (pix, sphaY, si, l0, m, ring_m);

      if (m < pix->lmax)
      {
        ncm_sf_spherical_harmonics_Y_next_m (sphaY);
        l0 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);

        if (l0 > pix->lmax)
          break;
      }
      else
        break;

      m++;
      si += pix->lmax - m + 1;      
    }
  }
  else
  {
    gint l0;
    ncm_sf_spherical_harmonics_start_rec (pix->spha, sphaY, theta_i);

    m  = 0;
    si = 0;

    l0 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
    while (TRUE)
		{
      const gint fft_ring_index = m % ring_size;
      complex double ring_m     = ((fft_ring_index <= ring_size_2) ? Fim[fft_ring_index] : conj (Fim[ring_size - fft_ring_index])) * pix_area;

      NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_run_over_l) (pix, sphaY, si, l0, m, ring_m);

      if (m < pix->lmax)
      {
        ncm_sf_spherical_harmonics_Y_next_m (sphaY);
        l0 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
        
        if (l0 > pix->lmax)
          break;
      }
      else
        break;

      m++;
      si += pix->lmax - m + 1;      
    }
  }

  ncm_sf_spherical_harmonics_Y_free (sphaY);
}

static void 
NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_map2alm_run) (NcmSphereMap *pix)
{
  const gint64 nrings   = ncm_sphere_map_get_nrings (pix);
  const gint64 nrings_2 = nrings / 2;
  gint64 r_i            = 0;
  gint64 nrleft;

  memset (&g_array_index (pix->alm, gsl_complex, 0), 0, sizeof (gsl_complex) * pix->alm->len);
  
  r_i = 0;
  while (TRUE)
  {
    if (r_i + NCM_SPHERE_MAP_BLOCK_NC < nrings_2)
    {
      NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_get_alm_from_apcircles) (pix, r_i, nrings);
      r_i += NCM_SPHERE_MAP_BLOCK_NC;
    }
    else
      break;
  }

  nrleft = nrings - r_i;
  for (; r_i < nrleft; r_i++)
  {
    NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_get_alm_from_circle) (pix, r_i);
  }
}

static void
NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_D_m) (NcmSphereMap *pix, NcmSFSphericalHarmonicsY *sphaY, gsl_complex * restrict * restrict alm_ptr, const gint ring_size, const gint m, gsl_complex *D_m)
{
	const gint lmaxmstepm2 = pix->lmax - NCM_SPHERE_MAP_BLOCK_INV_STEPM2;
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
      NCM_COMPLEX_INC_MUL_REAL (*D_m, alm_ptr[0][i], Ylm[i]);
    }

    alm_ptr[0] += NCM_SPHERE_MAP_BLOCK_INV_STEP;
    l = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
  }

  if (l <= pix->lmax)
  {
    while (TRUE)
    {
      const gdouble Ylm = ncm_sf_spherical_harmonics_Y_get_lm (sphaY);
      
      NCM_COMPLEX_INC_MUL_REAL (*D_m, alm_ptr[0][0], Ylm);

      alm_ptr[0]++;

      if (l < pix->lmax)
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
NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_D_m_array) (NcmSphereMap *pix, NcmSFSphericalHarmonicsYArray *sphaYa, gsl_complex * restrict * restrict alm_ptr, const gint m, gsl_complex *D_m)
{
	const gint lmaxmstepm2 = pix->lmax - NCM_SPHERE_MAP_BLOCK_INV_STEPM2;
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
				const gdouble Ylm_bj = Ylm[NCM_SF_SPHERICAL_HARMONICS_ARRAY_INDEX (j, i, NCM_SPHERE_MAP_BLOCK_INV_NCT, NCM_SPHERE_MAP_BLOCK_INV_STEP)];
				NCM_COMPLEX_INC_MUL_REAL (D_m[j], alm_ptr[0][i], Ylm_bj);
			}
		}
		
    alm_ptr[0] += NCM_SPHERE_MAP_BLOCK_INV_STEP;
    l = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
  }

  if (l <= pix->lmax)
  {
    while (TRUE)
    {
			gint j;
			for (j = 0; j < NCM_SPHERE_MAP_BLOCK_INV_NCT; j++)
			{
				const gdouble Ylm = ncm_sf_spherical_harmonics_Y_array_get_lm (sphaYa, j);
				NCM_COMPLEX_INC_MUL_REAL (D_m[j], alm_ptr[0][0], Ylm);
			}
			alm_ptr[0]++;

      if (l < pix->lmax)
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
  NcmSFSphericalHarmonicsY *sphaY = ncm_sf_spherical_harmonics_Y_new (pix->spha, NCM_SF_SPHERICAL_HARMONICS_DEFAULT_ABSTOL);
  const gint ring_size    = ncm_sphere_map_get_ring_size (pix, r_i);
  const gint ring_fi      = ncm_sphere_map_get_ring_first_index (pix, r_i);
  _fft_complex *Fim       = &((_fft_complex *)pix->fft_pvec)[ring_fi];
  gdouble theta_i         = 0.0, phi_0 = 0.0;
  gint ring_size_2        = ring_size / 2;
  gsl_complex *alm_ptr    = &g_array_index (pix->alm, gsl_complex, 0);
  gsl_complex Iphi, expIphi, phase;
  gint m, l0;

  ncm_sphere_map_pix2ang_ring (pix, ring_fi, &theta_i, &phi_0);

  Iphi    = gsl_complex_rect (0.0, phi_0);
  expIphi = gsl_complex_exp (Iphi);
  phase   = expIphi;
  
  ncm_sf_spherical_harmonics_start_rec (pix->spha, sphaY, theta_i);

  memset (Fim, 0, sizeof (_fft_complex) * ring_size);

  m  = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
  l0 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
  {
    gsl_complex D_m = NCM_COMPLEX_ZERO;
     NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_D_m) (pix, sphaY, &alm_ptr, ring_size, m, &D_m);
    Fim[m] += GSL_REAL (D_m);
  }

  ncm_sf_spherical_harmonics_Y_next_m (sphaY);
  m  = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
  l0 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);

  while (TRUE)
  {
    const gint fft_ring_index = m % ring_size;
    gsl_complex D_m           = NCM_COMPLEX_ZERO;

    NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_D_m) (pix, sphaY, &alm_ptr, ring_size, m, &D_m);

    NCM_COMPLEX_MUL (D_m, phase);
    NCM_COMPLEX_MUL (phase, expIphi);

    if ((fft_ring_index == 0) || (fft_ring_index == ring_size_2))
    {
      Fim[fft_ring_index] += 2.0 * GSL_REAL (D_m);
    }
    else if (fft_ring_index < ring_size_2)
    {
      Fim[fft_ring_index] += GSL_REAL (D_m) + I * GSL_IMAG (D_m);
    }
    else if (fft_ring_index > ring_size_2)
    {
      Fim[ring_size - fft_ring_index] += GSL_REAL (D_m) - I * GSL_IMAG (D_m);
    }

    if (m < pix->lmax)
    {
      ncm_sf_spherical_harmonics_Y_next_m (sphaY);
      m  = ncm_sf_spherical_harmonics_Y_get_m (sphaY);
      l0 = ncm_sf_spherical_harmonics_Y_get_l (sphaY);
      if (l0 > pix->lmax)
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
  const gint64 nc                       = NCM_SPHERE_MAP_BLOCK_INV_NC;
  const gint64 nct                      = NCM_SPHERE_MAP_BLOCK_INV_NCT;
  NcmSFSphericalHarmonicsYArray *sphaYa = ncm_sf_spherical_harmonics_Y_array_new (pix->spha, nct, NCM_SF_SPHERICAL_HARMONICS_ARRAY_DEFAULT_ABSTOL);
  gsl_complex *alm_ptr                  = &g_array_index (pix->alm, gsl_complex, 0);

	_fft_complex * restrict Fima[NCM_SPHERE_MAP_BLOCK_INV_NCT];
  gint64 ring_size[NCM_SPHERE_MAP_BLOCK_INV_NCT];
  gint64 ring_size_2[NCM_SPHERE_MAP_BLOCK_INV_NCT];
  gdouble theta[NCM_SPHERE_MAP_BLOCK_INV_NCT];
  gdouble phi[NCM_SPHERE_MAP_BLOCK_INV_NCT];
	gsl_complex D_m[NCM_SPHERE_MAP_BLOCK_INV_NCT];
  gsl_complex expIphi[NCM_SPHERE_MAP_BLOCK_INV_NCT];
  gsl_complex phase[NCM_SPHERE_MAP_BLOCK_INV_NCT];

	gint m, l0, i;

  for (i = 0; i < nc; i++)
  {
    const gint64 r_i   = r_ini + i;
    const gint64 apr_i = nrings - r_i - 1;
    const gint64 api   = nct - i - 1;
    
    NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_prepare_circle) (pix, r_i,   i,   Fima, theta, phi, ring_size, ring_size_2);
    NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_prepare_circle) (pix, apr_i, api, Fima, theta, phi, ring_size, ring_size_2);

		{	
		  gsl_complex Iphi_i   = gsl_complex_rect (0.0, phi[i]);
		  gsl_complex Iphi_api = gsl_complex_rect (0.0, phi[api]);
		  expIphi[i]           = gsl_complex_exp (Iphi_i);
		  expIphi[api]         = gsl_complex_exp (Iphi_api);
		  phase[i]             = expIphi[i];
		  phase[api]           = expIphi[api];
	  }

		memset (Fima[i],   0, sizeof (_fft_complex) * ring_size[i]);
		memset (Fima[api], 0, sizeof (_fft_complex) * ring_size[api]);
	}

  ncm_sf_spherical_harmonics_start_rec_array (pix->spha, sphaYa, nct, theta);

  m  = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
  l0 = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);

	memset (D_m, 0, sizeof (gsl_complex) * NCM_SPHERE_MAP_BLOCK_INV_NCT);

	NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_D_m_array) (pix, sphaYa, &alm_ptr, m, D_m);

	for (i = 0; i < nct; i++)
  {
		Fima[i][m] += GSL_REAL (D_m[i]);
  }

  ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, nct);
  m  = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
  l0 = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);

  while (TRUE)
  {
		memset (D_m, 0, sizeof (gsl_complex) * NCM_SPHERE_MAP_BLOCK_INV_NCT);
    NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_get_D_m_array) (pix, sphaYa, &alm_ptr, m, D_m);

		for (i = 0; i < nct; i++)
		{
			const gint fft_ring_index = m % ring_size[i];

			NCM_COMPLEX_MUL (D_m[i], phase[i]);
			NCM_COMPLEX_MUL (phase[i], expIphi[i]);

			if ((fft_ring_index == 0) || (fft_ring_index == ring_size_2[i]))
			{
				Fima[i][fft_ring_index] += 2.0 * GSL_REAL (D_m[i]);
			}
			else if (fft_ring_index < ring_size_2[i])
			{
				Fima[i][fft_ring_index] += GSL_REAL (D_m[i]) + I * GSL_IMAG (D_m[i]);
			}
			else if (fft_ring_index > ring_size_2[i])
			{
				Fima[i][ring_size[i] - fft_ring_index] += GSL_REAL (D_m[i]) - I * GSL_IMAG (D_m[i]);
			}
		}
		if (m < pix->lmax)
    {
      ncm_sf_spherical_harmonics_Y_array_next_m (sphaYa, nct);
      m  = ncm_sf_spherical_harmonics_Y_array_get_m (sphaYa);
      l0 = ncm_sf_spherical_harmonics_Y_array_get_l (sphaYa);
      if (l0 > pix->lmax)
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
