/***************************************************************************
 *            ncm_sphere_map_pix.c
 *
 *  Wed Jul  9 11:09:37 2008 (updated Jul 15 2016)
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_sphere_map_pix.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_sphere_map_pix
 * @title: NcmSphereMapPix
 * @short_description: An re-implementation of Healpix.
 * 
 * Map pixalization/manipulation algorithms, Ylm decomposition.
 * 
 */
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sphere_map_pix.h"
#include "math/ncm_util.h"
#include "math/ncm_cfg.h"
#include "math/ncm_c.h"
#include "ncm_enum_types.h"

#include <gsl/gsl_sf_legendre.h>

#include <math.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */

#ifdef NUMCOSMO_HAVE_FFTW3 
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */

#ifndef HAVE_FFTW3_ALLOC
#define fftwf_alloc_real(n) (double *) fftwf_malloc (sizeof(double) * (n))
#define fftwf_alloc_complex(n) (fftwf_complex *) fftwf_malloc (sizeof(fftw_complex) * (n))
#endif /* HAVE_FFTW3_ALLOC */

#ifdef HAVE_FFTW3F
#  define _fft_vec_alloc fftwf_alloc_real 
#  define _fft_complex complex float 
#  define _fft_vec_alloc_complex fftwf_alloc_complex
#  define _fft_vec_free  fftwf_free
#  define _fft_vec_set_zero(v,s) memset ((v), 0, sizeof (gfloat) * (s))
#  define _fft_vec_memcpy(dest,orig,s) memcpy ((dest), (orig), sizeof (gfloat) * (s))
#  define _fft_vec_ptr(v,i) (&((gfloat *)(v))[i])
#elif defined (HAVE_FFTW3)
#  define _fft_vec_alloc fftw_alloc_real 
#  define _fft_complex complex double 
#  define _fft_vec_alloc_complex fftw_alloc_complex 
#  define _fft_vec_free  fftw_free
#  define _fft_vec_set_zero(v,s) memset ((v), 0, sizeof (gdouble) * (s))
#  define _fft_vec_memcpy(dest,orig,s) memcpy ((dest), (orig), sizeof (gdouble) * (s))
#  define _fft_vec_ptr(v,i) (&((gdouble *)(v))[i])
#else
#  define _fft_vec_alloc gsl_vector_float_alloc 
#  define _fft_vec_free  gsl_vector_float_free
#  define _fft_vec_set_zero(v,s) gsl_vector_float_set_zero (v)
#  define _fft_vec_memcpy(dest,orig,s) gsl_vector_float_memcpy ((dest),(orig))
#  define _fft_vec_ptr(v,i) (gsl_vector_float_ptr ((v),(i)))
#endif
#define _fft_vec_idx(v,i) (*_fft_vec_ptr(v,i))

#include <gsl/gsl_sf_trig.h>

enum
{
  PROP_0,
  PROP_NSIDE,
  PROP_ORDER,
  PROP_COORDSYS,
};

G_DEFINE_TYPE (NcmSphereMapPix, ncm_sphere_map_pix, G_TYPE_OBJECT);

static void
ncm_sphere_map_pix_init (NcmSphereMapPix *pix)
{
  pix->nside             = 0;
  pix->npix              = 0;
  pix->face_size         = 0;
  pix->middle_rings_size = 0;
  pix->cap_size          = 0;
  pix->middle_size       = 0;
  pix->nrings            = 0;
  pix->nrings_cap        = 0;
  pix->nrings_middle     = 0;
  pix->order             = NCM_SPHERE_MAP_PIX_ORDER_RING;
  pix->coordsys          = NCM_SPHERE_MAP_PIX_COORD_SYS_LEN;
  pix->pvec              = NULL;
  pix->fft_pvec          = NULL;
  pix->fft_plan          = g_ptr_array_new ();
#ifdef NUMCOSMO_HAVE_FFTW3
#  ifdef HAVE_FFTW3F
  g_ptr_array_set_free_func (pix->fft_plan, (GDestroyNotify)fftwf_destroy_plan);
#  else
  g_ptr_array_set_free_func (pix->fft_plan, (GDestroyNotify)fftw_destroy_plan);
#  endif
#endif
}

static void
_ncm_sphere_map_pix_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSphereMapPix *pix = NCM_SPHERE_MAP_PIX (object);
  g_return_if_fail (NCM_IS_SPHERE_MAP_PIX (object));

  switch (prop_id)
  {
    case PROP_NSIDE:
      ncm_sphere_map_pix_set_nside (pix, g_value_get_int64 (value));    
      break;
    case PROP_ORDER:
      ncm_sphere_map_pix_set_order (pix, g_value_get_enum (value));
      break;
    case PROP_COORDSYS:
      ncm_sphere_map_pix_set_coordsys (pix, g_value_get_enum (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_sphere_map_pix_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSphereMapPix *pix = NCM_SPHERE_MAP_PIX (object);
  g_return_if_fail (NCM_IS_SPHERE_MAP_PIX (object));

  switch (prop_id)
  {
    case PROP_NSIDE:
      g_value_set_int64 (value, pix->nside);
      break;
    case PROP_ORDER:
      g_value_set_enum (value, ncm_sphere_map_pix_get_order (pix));
      break;
    case PROP_COORDSYS:
      g_value_set_enum (value, ncm_sphere_map_pix_get_coordsys (pix));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_sphere_map_pix_finalize (GObject *object)
{
  NcmSphereMapPix *pix = NCM_SPHERE_MAP_PIX (object);

  ncm_sphere_map_pix_set_nside (pix, 0);
  g_ptr_array_unref (pix->fft_plan);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sphere_map_pix_parent_class)->finalize (object);
}

static void
ncm_sphere_map_pix_class_init (NcmSphereMapPixClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_sphere_map_pix_set_property;
  object_class->get_property = &_ncm_sphere_map_pix_get_property;
  object_class->finalize     = &_ncm_sphere_map_pix_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NSIDE,
                                   g_param_spec_int64 ("nside",
                                                       NULL,
                                                       "nside",
                                                       0, G_MAXINT64, 0,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ORDER,
                                   g_param_spec_enum ("order",
                                                      NULL,
                                                      "Map pixel ordering",
                                                      NCM_TYPE_SPHERE_MAP_PIX_ORDER, NCM_SPHERE_MAP_PIX_ORDER_RING,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_COORDSYS,
                                   g_param_spec_enum ("coordsys",
                                                      NULL,
                                                      "Map coordinate system",
                                                      NCM_TYPE_SPHERE_MAP_PIX_COORD_SYS, NCM_SPHERE_MAP_PIX_COORD_SYS_CELESTIAL,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
}

/**
 * ncm_sphere_map_pix_new:
 * @nside: FIXME
 * 
 * FIXME
 * 
 * 
 * Returns: (transfer full): FIXME
 */
NcmSphereMapPix *
ncm_sphere_map_pix_new (const gint64 nside)
{
  NcmSphereMapPix *pix = g_object_new (NCM_TYPE_SPHERE_MAP_PIX,
                                       "nside", nside,
                                       NULL);

  return pix;
}

/**
 * ncm_sphere_map_pix_ref:
 * @pix: a #NcmSphereMapPix
 * 
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
NcmSphereMapPix *
ncm_sphere_map_pix_ref (NcmSphereMapPix *pix)
{
  return g_object_ref (pix);
}

/**
 * ncm_sphere_map_pix_free:
 * @pix: a #NcmSphereMapPix
 * 
 * FIXME
 * 
 */
void 
ncm_sphere_map_pix_free (NcmSphereMapPix *pix)
{
  g_object_unref (pix);
}

/**
 * ncm_sphere_map_pix_clear:
 * @pix: a #NcmSphereMapPix
 * 
 * FIXME
 * 
 */
void 
ncm_sphere_map_pix_clear (NcmSphereMapPix **pix)
{
  g_clear_object (pix);
}

static gint64 _l_pow_2 (gint64 n) { return n * n; }

/**
 * ncm_sphere_map_pix_set_nside:
 * @pix: a #NcmSphereMapPix
 * @nside: FIXME
 * 
 * FIXME
 * 
 */
void
ncm_sphere_map_pix_set_nside (NcmSphereMapPix *pix, gint64 nside)
{
  if (nside > 0)
    g_assert_cmpint (nside, ==, (gint64) exp2 ((gint64) log2 (nside)));

  if (nside != pix->nside)
  {
    g_clear_pointer (&pix->pvec, (GDestroyNotify) _fft_vec_free);

    pix->nside             = nside;
    pix->npix              = 0;
    pix->face_size         = 0;
    pix->middle_rings_size = 0;
    pix->cap_size          = 0;
    pix->middle_size       = 0;
    pix->nrings            = 0;
    pix->nrings_cap        = 0;
    pix->nrings_middle     = 0;

#ifdef NUMCOSMO_HAVE_FFTW3
#  ifdef HAVE_FFTW3F
    g_clear_pointer (&pix->fft_pvec, (GDestroyNotify) fftwf_free);
#  else
    g_clear_pointer (&pix->fft_pvec, (GDestroyNotify) fftw_free);
#  endif
#endif
    g_ptr_array_set_size (pix->fft_plan, 0);
    
    if (nside > 0)
    {
      pix->npix              = NCM_SPHERE_MAP_PIX_N (pix->nside);
      pix->face_size         = _l_pow_2 (pix->nside);
      pix->middle_rings_size = 4 * nside;
      pix->cap_size          = 2 * nside * (nside - 1);
      pix->nrings            = 4 * nside - 1;
      pix->nrings_cap        = nside - 1;
      pix->nrings_middle     = 2 * nside + 1;
      pix->middle_size       = pix->nrings_middle * pix->middle_rings_size;      
      pix->pvec              = _fft_vec_alloc (pix->npix);
      
      _fft_vec_set_zero (pix->pvec, pix->npix);
      
      g_assert_cmpint (2 * pix->cap_size + pix->middle_size, ==, pix->npix);
      
    }
  }
}

/**
 * ncm_sphere_map_pix_get_nside:
 * @pix: a #NcmSphereMapPix
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
gint64
ncm_sphere_map_pix_get_nside (NcmSphereMapPix *pix)
{
  return pix->nside;
}

/**
 * ncm_sphere_map_pix_get_npix:
 * @pix: a #NcmSphereMapPix
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
gint64
ncm_sphere_map_pix_get_npix (NcmSphereMapPix *pix)
{
  return pix->npix;
}

/**
 * ncm_sphere_map_pix_get_cap_size:
 * @pix: a #NcmSphereMapPix
 *
 * FIXME 
 * 
 * Returns: FIXME
 */
gint64 
ncm_sphere_map_pix_get_cap_size (NcmSphereMapPix *pix)
{
  return pix->cap_size;
}

/**
 * ncm_sphere_map_pix_get_middle_size:
 * @pix: a #NcmSphereMapPix
 *
 * FIXME 
 * 
 * Returns: FIXME
 */
gint64 
ncm_sphere_map_pix_get_middle_size (NcmSphereMapPix *pix)
{
  return pix->middle_size;
}

/**
 * ncm_sphere_map_pix_get_nrings:
 * @pix: a #NcmSphereMapPix
 *
 * FIXME 
 * 
 * Returns: FIXME
 */
gint64 
ncm_sphere_map_pix_get_nrings (NcmSphereMapPix *pix)
{
  return pix->nrings;
}

/**
 * ncm_sphere_map_pix_get_nrings_cap:
 * @pix: a #NcmSphereMapPix
 *
 * FIXME 
 * 
 * Returns: FIXME
 */
gint64 
ncm_sphere_map_pix_get_nrings_cap (NcmSphereMapPix *pix)
{
  return pix->nrings_cap;
}

/**
 * ncm_sphere_map_pix_get_nrings_middle:
 * @pix: a #NcmSphereMapPix
 *
 * FIXME 
 * 
 * Returns: FIXME
 */
gint64 
ncm_sphere_map_pix_get_nrings_middle (NcmSphereMapPix *pix)
{
  return pix->nrings_middle;
}

/**
 * ncm_sphere_map_pix_get_ring_size:
 * @pix: a #NcmSphereMapPix
 * @r_i: ring index
 *
 * FIXME 
 * 
 * Returns: FIXME
 */
gint64 
ncm_sphere_map_pix_get_ring_size (NcmSphereMapPix *pix, gint64 r_i)
{
  if (r_i < pix->nrings_cap) /* North cap (nside - 1 rings) */
  {
    return 4 * (r_i + 1); 
  }
  else if (r_i >= pix->nrings_cap + pix->nrings_middle) /* South cap (it has 3nside rings above) */
  {
    const gint64 rsc_i = r_i - (pix->nrings_cap + pix->nrings_middle); /* South cap index */
    return 4 * (pix->nside - 1 - rsc_i);
  }
  else /* Middle */
  {
    return pix->middle_rings_size;
  }
}

/**
 * ncm_sphere_map_pix_get_ring_first_index:
 * @pix: a #NcmSphereMapPix
 * @r_i: ring index
 *
 * FIXME 
 * 
 * Returns: FIXME
 */
gint64 
ncm_sphere_map_pix_get_ring_first_index (NcmSphereMapPix *pix, gint64 r_i)
{
  if (r_i < pix->nrings_cap) /* North cap (nside - 1 rings) */
  {
    return 2 * r_i * (r_i + 1); 
  }
  else if (r_i >= pix->nrings_cap + pix->nrings_middle) /* South cap (it has 3nside rings above) */
  {
    const gint64 rsc_i = r_i - (pix->nrings_cap + pix->nrings_middle); /* South cap index */
    return pix->cap_size + pix->middle_size + 2 * rsc_i * (2 * pix->nside - 1 - rsc_i);
  }
  else /* Middle */
  {
    const gint64 rm_i = r_i - pix->nrings_cap; /* Middle index */
    return pix->cap_size + pix->middle_rings_size * rm_i;
  }  
}

/**
 * ncm_sphere_map_pix_set_order:
 * @pix: a #NcmSphereMapPix
 * @order: FIXME
 * 
 * FIXME
 * 
 */
void 
ncm_sphere_map_pix_set_order (NcmSphereMapPix *pix, NcmSphereMapPixOrder order)
{
  if (pix->order != order)
  {
    if (pix->nside > 0)
    {
      gpointer temp_pix = _fft_vec_alloc (pix->npix);
      gint64 (*convert) (NcmSphereMapPix *pix, gint64) = (order == NCM_SPHERE_MAP_PIX_ORDER_NEST ? ncm_sphere_map_pix_ring2nest : ncm_sphere_map_pix_nest2ring);
      gint64 i, j;

      _fft_vec_memcpy (temp_pix, pix->pvec, pix->npix);
      for (i = 0; i < pix->npix; i++)
      {
        gfloat val = _fft_vec_idx (temp_pix, i);
        j = convert (pix, i);
        printf ("%ld => %ld val % 20.15g | % 20.15g\n", i, j, val, _fft_vec_idx (pix->pvec, i));
        _fft_vec_idx (pix->pvec, j) = val;
      }
      
      _fft_vec_free (temp_pix);
    }
    pix->order = order;
  }
}

/**
 * ncm_sphere_map_pix_get_order:
 * @pix: a #NcmSphereMapPix
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
NcmSphereMapPixOrder 
ncm_sphere_map_pix_get_order (NcmSphereMapPix *pix)
{
  return pix->order;
}

/**
 * ncm_sphere_map_pix_set_coordsys:
 * @pix: a #NcmSphereMapPix
 * @coordsys: FIXME
 * 
 * FIXME
 * 
 */
void 
ncm_sphere_map_pix_set_coordsys (NcmSphereMapPix *pix, NcmSphereMapPixCoordSys coordsys)
{
  if (pix->coordsys != coordsys)
  {
    switch (coordsys)
    {
      case NCM_SPHERE_MAP_PIX_COORD_SYS_GALACTIC:
      case NCM_SPHERE_MAP_PIX_COORD_SYS_ECLIPTIC:
      case NCM_SPHERE_MAP_PIX_COORD_SYS_CELESTIAL:
        pix->coordsys = coordsys;
        break;
      default:
        g_error ("ncm_sphere_map_pix_set_coordsys: Unknown coordinate system `%d'.", coordsys);
        break;
    }  
  }
}

/**
 * ncm_sphere_map_pix_get_coordsys:
 * @pix: a #NcmSphereMapPix
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
NcmSphereMapPixCoordSys
ncm_sphere_map_pix_get_coordsys (NcmSphereMapPix *pix)
{
  return pix->coordsys;
}

/**
 * ncm_sphere_map_pix_clear_pixels:
 * @pix: a #NcmSphereMapPix
 * 
 * FIXME
 * 
 */
void 
ncm_sphere_map_pix_clear_pixels (NcmSphereMapPix *pix)
{
  if (pix->nside > 0)
    _fft_vec_set_zero (pix->pvec, pix->npix);
}

/**
 * ncm_sphere_map_pix_nest2ring: 
 * @pix: a #NcmSphereMapPix
 * @nest_index: FIXME
 * 
 * FIXME
 *
 * Returns: FIXME
 */
gint64
ncm_sphere_map_pix_nest2ring (NcmSphereMapPix *pix, const gint64 nest_index)
{
  gint64 ring_index;
  gint64 base_pixel;
  gint64 f, h, v;            /* Face number, horizontal coordinate, vertical coordinate */
  gint64 t, p, s, pad, w, l; /* theta, phi, shift, padding, width, local index          */
  gint64 x, y;
  gint64 hf;

  g_assert (nest_index < pix->npix);
  
  f = nest_index / pix->face_size;
  l = nest_index % pix->face_size;
  
  NCM_SPHERE_MAP_PIX_INT_TO_XY (l, x, y);

  h = pix->nside - 1 - y;
  v = 2 * pix->nside - 2 - y - x;
  
  switch (f / 4)
  {
    case 0:
      t = v;
      if (t < (pix->nside - 1))
      {
        w   = (t + 1);
        s   = 0;
        pad = 0;
        hf  = f % 4;
        base_pixel = 2 * (_l_pow_2 (t + 1) - t - 1);
      }
      else
      {
        w   = pix->nside;
        s   = (t - (pix->nside - 1)) / 2;
        pad = 0;
        hf  = f % 4;
        base_pixel = pix->cap_size + (t - pix->nside + 1) * pix->middle_rings_size;
      }
      break;
    case 1:
      t   = v + pix->nside;
      w   = pix->nside;
      s   = (t - (pix->nside - 1)) / 2;
      pad = 0;
      hf  = f % 4;
      base_pixel = pix->cap_size + (v + 1) * pix->middle_rings_size;
      break;
    case 2:
      t   = v + 2 * pix->nside;
      if (v < pix->nside)
      {
        w   = pix->nside;
        s   = (t - (pix->nside - 1)) / 2;
        pad = 0;
        hf  = f % 4 + 1;
        base_pixel = pix->cap_size + (v + 1 + pix->nside) * pix->middle_rings_size;
      }
      else
      {
        w   = 4 * pix->nside - t - 1;
        s   = w;
        pad = (t - 3 * pix->nside + 1);
        hf  = f % 4 + 1;
        base_pixel = (pix->npix - pix->cap_size) + (4 * pix->nside * (v - pix->nside + 1) - 2 * _l_pow_2 (v - pix->nside + 1) + 2 * (v - pix->nside) + 2 - 4 * pix->nside);
      }
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  p = (h + hf * w - s - pad);
  if (p < 0) p += 4 * w;

  ring_index = base_pixel + p;

  return ring_index;
}

/**
 * ncm_sphere_map_pix_ring2nest: 
 * @pix: a #NcmSphereMapPix
 * @ring_index: FIXME
 * 
 * FIXME
 *
 * Returns: FIXME
 */
gint64
ncm_sphere_map_pix_ring2nest (NcmSphereMapPix *pix, const gint64 ring_index)
{
  gint64 nest_index;
  gint64 f, h, v;            /* Face number, horizontal coordinate, vertical coordinate */
  gint64 t, p, s, pad, w, l; /* theta, phi, shift, padding, width, local index          */
  gint64 x, y;
  
  g_assert (ring_index < pix->npix);

  if (ring_index < pix->cap_size)
  {
    t   = (sqrt (1 + 2 * ring_index) - 1) / 2;
    w   = (t + 1);
    p   = ring_index - 2 * (_l_pow_2 (t + 1) - t - 1);
    s   = 0;
    pad = 0;
  }
  else if (ring_index < (pix->npix - pix->cap_size))
  {
    l   = ring_index - pix->cap_size;
    w   = pix->nside;
    t   = (gint64)(l / pix->middle_rings_size) + (pix->nside - 1);
    p   = l % pix->middle_rings_size;
    s   = (t - (pix->nside - 1)) / 2;
    pad = 0;
  }
  else
  {
    l   = ring_index - pix->npix + pix->cap_size;
    t   = 4 * pix->nside - (1 + sqrt (4 * _l_pow_2 (pix->nside) - 4 * pix->nside + 1 - 2 * l)) / 2;
    w   = 4 * pix->nside - t - 1;
    p   = l - (4 * pix->nside * (t - 3 * pix->nside + 1) - 2 * _l_pow_2 (t - 3 * pix->nside + 1) + 2 * (t - 3 * pix->nside) + 2 - 4 * pix->nside);
    s   = w;
    pad = (t - 3 * pix->nside + 1);
  }

  h  = (p + s) % w + pad;
  v  = (t - h) % pix->nside + h;
  f  = ((gint64)((t - h) / pix->nside)) * 4;
  f += (f == 8) ? ((gint64)((p + s) / w) - 1) : ((gint64)((p + s) / w) % 4);
  
  x = pix->nside - 1 - (v - h);
  y = pix->nside - 1 - h;
  
  NCM_SPHERE_MAP_PIX_XY_TO_INT (x, y, nest_index);
  
  nest_index += pix->face_size * f;
  
  return nest_index;
}

static void
_t_p_w_to_theta_phi (const gint64 nside, const gint tm1, const gint pm1, const gint w, gdouble *theta, gdouble *phi)
{
  const gint64 t = tm1 + 1;
  const gint64 p = pm1 + 1;
  
  if (t < nside)
  {
    *theta = acos ((1.0 - t * t * 1.0 / (3.0 * nside * nside)));
    *phi   = (p - 0.5) * M_PI_2 / (1.0 * w);
  }
  else if (t > 3 * nside)
  {
    gint tt = 4 * nside - t;

    *theta = acos (-(1.0 - tt * tt * 1.0 / (3.0 * nside * nside)) );
    *phi   = (p - 0.5) * M_PI_2 / (1.0 * w);
  }
  else
  {
    *theta = acos ((2.0 * nside - t) * 2.0 / (3.0 * nside));
    *phi   = (p - ((t - nside) % 2 + 1.0) * 0.5) * M_PI_2 / (1.0 * w);
  }
}

static void
_t_p_w_to_vector (const gint64 nside, const gint tm1, const gint pm1, const gint w, NcmTriVec *vec)
{
  const gint64 t = tm1 + 1;
  const gint64 p = pm1 + 1;
  gdouble phi, z, sin_theta;

  if (t < nside)
  {
    phi = (p - 0.5) * M_PI_2 / (1.0 * w);
    z   = (1.0 - t * t * 1.0 / (3.0 * nside * nside));
  }
  else if (t > 3 * nside)
  {
    const gint tt = 4 * nside - t;
    phi = (p - 0.5) * M_PI_2 / (1.0 * w);
    z   = -(1.0 - tt * tt * 1.0 / (3.0 * nside * nside));
  }
  else
  {
    phi = (p - ((t - nside) % 2 + 1.0) * 0.5) * M_PI_2 / (1.0 * w);
    z   = (2.0 * nside - t) * 2.0 / (3.0 * nside);
  }
  
  sin_theta = sqrt (1.0 - z * z);
  vec->c[0] = sin_theta * cos (phi);
  vec->c[1] = sin_theta * sin (phi);
  vec->c[2] = z;
}

/**
 * ncm_sphere_map_pix_pix2ang_nest:
 * @pix: a #NcmSphereMapPix
 * @nest_index: FIXME
 * @theta: (out): FIXME
 * @phi: (out): FIXME
 *
 * FIXME 
*/ 
void 
ncm_sphere_map_pix_pix2ang_nest (NcmSphereMapPix *pix, const gint64 nest_index, gdouble *theta, gdouble *phi)
{
  gint64 f, h, v;            /* Face number, horizontal coordinate, vertical coordinate */
  gint64 t, p, s, pad, w, l; /* theta, phi, shift, padding, width, local index          */
  gint64 x, y;
  gint64 hf;

  g_assert (nest_index < pix->npix);
  
  f = nest_index / pix->face_size;
  l = nest_index % pix->face_size;
  
  NCM_SPHERE_MAP_PIX_INT_TO_XY (l, x, y);
 
  h = pix->nside - 1 - y;
  v = 2 * pix->nside - 2 - y - x;
  
  switch (f / 4)
  {
    case 0:
      t = v;
      if (t < (pix->nside - 1))
      {
        w   = (t + 1);
        s   = 0;
        pad = 0;
        hf  = f % 4;
      }
      else
      {
        w   = pix->nside;
        s   = (t - (pix->nside - 1)) / 2;
        pad = 0;
        hf  = f % 4;
      }
      break;
    case 1:
      t   = v + pix->nside;
      w   = pix->nside;
      s   = (t - (pix->nside - 1)) / 2;
      pad = 0;
      hf  = f % 4;
      break;
    case 2:
      t   = v + 2 * pix->nside;
      if (v < pix->nside)
      {
        w   = pix->nside;
        s   = (t - (pix->nside - 1)) / 2;
        pad = 0;
        hf  = f % 4 + 1;
      }
      else
      {
        w   = 4 * pix->nside - t - 1;
        s   = w;
        pad = (t - 3 * pix->nside + 1);
        hf  = f % 4 + 1;
      }
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  p = (h + hf * w - s - pad);
  if (p < 0) p += 4 * w;
  
  _t_p_w_to_theta_phi (pix->nside, t, p, w, theta, phi);
}

/**
 * ncm_sphere_map_pix_pix2ang_ring:
 * @pix: a #NcmSphereMapPix
 * @ring_index: FIXME
 * @theta: (out): FIXME
 * @phi: (out): FIXME
 *
 * FIXME 
*/ 
void
ncm_sphere_map_pix_pix2ang_ring (NcmSphereMapPix *pix, const gint64 ring_index, gdouble *theta, gdouble *phi)
{
  gint64 t, p, w, l; /* theta, phi, shift, padding, width, local index */
  
  g_assert (ring_index < pix->npix);

  if (ring_index < pix->cap_size)
  {
    t   = (sqrt (1 + 2 * ring_index) - 1) / 2;
    w   = (t + 1);
    p   = ring_index - 2 * (_l_pow_2 (t + 1) - t - 1);
  }
  else if (ring_index < (pix->npix - pix->cap_size))
  {
    l   = ring_index - pix->cap_size;
    w   = pix->nside;
    t   = (gint64)(l / pix->middle_rings_size) + (pix->nside - 1);
    p   = l % pix->middle_rings_size;
  }
  else
  {
    l   = ring_index - pix->npix + pix->cap_size;
    t   = 4 * pix->nside - (1 + sqrt (4 * _l_pow_2 (pix->nside) - 4 * pix->nside + 1 - 2 * l)) / 2;
    w   = 4 * pix->nside - t - 1;
    p   = l - (4 * pix->nside * (t - 3 * pix->nside + 1) - 2 * _l_pow_2 (t - 3 * pix->nside + 1) + 2 * (t - 3 * pix->nside) + 2 - 4 * pix->nside);
  }
  
  _t_p_w_to_theta_phi (pix->nside, t, p, w, theta, phi);
}

/**
 * ncm_sphere_map_pix_pix2vec_ring:
 * @pix: a #NcmSphereMapPix
 * @ring_index: FIXME
 * @vec: a #NcmTriVec
 *
 * FIXME 
*/
void 
ncm_sphere_map_pix_pix2vec_ring (NcmSphereMapPix *pix, gint64 ring_index, NcmTriVec *vec)
{
  gint64 t, p, w, l; /* theta, phi, shift, padding, width, local index */
  
  g_assert (ring_index < pix->npix);

  if (ring_index < pix->cap_size)
  {
    t   = (sqrt (1 + 2 * ring_index) - 1) / 2;
    w   = (t + 1);
    p   = ring_index - 2 * (_l_pow_2 (t + 1) - t - 1);
  }
  else if (ring_index < (pix->npix - pix->cap_size))
  {
    l   = ring_index - pix->cap_size;
    w   = pix->nside;
    t   = (gint64)(l / pix->middle_rings_size) + (pix->nside - 1);
    p   = l % pix->middle_rings_size;
  }
  else
  {
    l   = ring_index - pix->npix + pix->cap_size;
    t   = 4 * pix->nside - (1 + sqrt (4 * _l_pow_2 (pix->nside) - 4 * pix->nside + 1 - 2 * l)) / 2;
    w   = 4 * pix->nside - t - 1;
    p   = l - (4 * pix->nside * (t - 3 * pix->nside + 1) - 2 * _l_pow_2 (t - 3 * pix->nside + 1) + 2 * (t - 3 * pix->nside) + 2 - 4 * pix->nside);
  }
  
  _t_p_w_to_vector (pix->nside, t, p, w, vec);
}

/**
 * ncm_sphere_map_pix_pix2vec_nest:
 * @pix: a #NcmSphereMapPix
 * @nest_index: FIXME
 * @vec: a #NcmTriVec
 *
 * FIXME 
*/
void 
ncm_sphere_map_pix_pix2vec_nest (NcmSphereMapPix *pix, gint64 nest_index, NcmTriVec *vec)
{
  gint64 f, h, v;            /* Face number, horizontal coordinate, vertical coordinate */
  gint64 t, p, s, pad, w, l; /* theta, phi, shift, padding, width, local index          */
  gint64 x, y;
  gint64 hf;

  g_assert (nest_index < pix->npix);
  
  f = nest_index / pix->face_size;
  l = nest_index % pix->face_size;

  NCM_SPHERE_MAP_PIX_INT_TO_XY (l, x, y);
 
  h = pix->nside - 1 - y;
  v = 2 * pix->nside - 2 - y - x;
  
  switch (f / 4)
  {
    case 0:
      t = v;
      if (t < (pix->nside - 1))
      {
        w   = (t + 1);
        s   = 0;
        pad = 0;
        hf  = f % 4;
      }
      else
      {
        w   = pix->nside;
        s   = (t - (pix->nside - 1)) / 2;
        pad = 0;
        hf  = f % 4;
      }
      break;
    case 1:
      t   = v + pix->nside;
      w   = pix->nside;
      s   = (t - (pix->nside - 1)) / 2;
      pad = 0;
      hf  = f % 4;
      break;
    case 2:
      t   = v + 2 * pix->nside;
      if (v < pix->nside)
      {
        w   = pix->nside;
        s   = (t - (pix->nside - 1)) / 2;
        pad = 0;
        hf  = f % 4 + 1;
      }
      else
      {
        w   = 4 * pix->nside - t - 1;
        s   = w;
        pad = (t - 3 * pix->nside + 1);
        hf  = f % 4 + 1;
      }
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  p = (h + hf * w - s - pad);
  if (p < 0) p += 4 * w;
  
  _t_p_w_to_vector (pix->nside, t, p, w, vec);
}

static void 
_ncm_sphere_map_pix_zphi2pix_nest (NcmSphereMapPix *pix, const gdouble z, const gdouble onemz2, const gdouble phi, gint64 *nest_index)
{
  const gdouble two_3 = 2.0 / 3.0;
  const gdouble abs_z = fabs (z);
  const gdouble tt    = fmod (phi / M_PI_2, 4.0);
  gint64 x, y, f;

  if (abs_z <= two_3)
  {
    const gdouble temp1 = pix->nside * (0.5 + tt);
    const gdouble temp2 = pix->nside * (z * 0.75);

    const gint64 jp  = (gint64) (temp1 - temp2);
    const gint64 jm  = (gint64) (temp1 + temp2);
    const gint64 ifp = jp / pix->nside;
    const gint64 ifm = jm / pix->nside;

    f = (ifp == ifm) ? (ifp | 4) : ((ifp < ifm) ? ifp : (ifm + 8));
    x = jm & (pix->nside - 1);
    y = pix->nside - (jp & (pix->nside - 1)) - 1;
  }
  else
  {
    const gint64 ntt  = GSL_MIN (3, (gint64) (tt));
    const gdouble tp  = tt - ntt;
    const gdouble tmp = pix->nside * sqrt (3.0 * onemz2 / (1.0 + abs_z));

    const gint64 jp  = (gint64)(tp * tmp);
    const gint64 jm  = (gint64)((1.0 - tp) * tmp);
    const gint64 jpm  = GSL_MIN (jp, pix->nside - 1);
    const gint64 jmm  = GSL_MIN (jm, pix->nside - 1);

    if (z >= 0.0)
    {
      x = pix->nside - jmm - 1;
      y = pix->nside - jpm - 1;
      f = ntt;
    }
    else
    {
      x = jpm;
      y = jmm;
      f = ntt + 8;
    }
  }

  NCM_SPHERE_MAP_PIX_XY_TO_INT (x, y, nest_index[0]);

  nest_index[0] += pix->face_size * f;
}

static void 
_ncm_sphere_map_pix_zphi2pix_ring (NcmSphereMapPix *pix, const gdouble z, const gdouble onemz2, const gdouble phi, gint64 *ring_index)
{
  const gdouble two_3 = 2.0 / 3.0;
  const gdouble abs_z = fabs (z);
  const gdouble tt    = fmod (phi / M_PI_2, 4.0);

  if (abs_z <= two_3)
  {
    const gdouble temp1 = pix->nside * (0.5 + tt);
    const gdouble temp2 = pix->nside * (z * 0.75);

    const gint64 jp  = (gint64) (temp1 - temp2);
    const gint64 jm  = (gint64) (temp1 + temp2);

    const gint64 ir     = pix->nside + 1 + jp - jm;
    const gint64 kshift = 1 - (ir & 1);

    const gint64 t1 = (jp + jm - pix->nside + kshift + 1) / 2;
    const gint64 ip = t1 % pix->middle_rings_size;

    ring_index[0] = pix->cap_size + (ir - 1) * pix->middle_rings_size + ip;
  }
  else  // North & South polar caps
  {
    const gdouble tp  = tt - (gint64)(tt);
    const gdouble tmp = pix->nside * sqrt (3.0 * onemz2 / (1.0 + abs_z));;

    const gint64 jp = (gint64)(tp * tmp);
    const gint64 jm = (gint64)((1.0 - tp) * tmp);

    const gint64 ir = jp + jm + 1;
    const gint64 ip = (gint64) (tt * ir);

    if (z > 0)
    {
      ring_index[0] = 2 * ir * (ir - 1) + ip;
    }
    else
    {
      ring_index[0] = pix->npix - 2 * ir * (ir + 1) + ip;
    }
  }
}

/**
 * ncm_sphere_map_pix_ang2pix_nest:
 * @pix: a #NcmSphereMapPix
 * @theta: FIXME
 * @phi: FIXME
 * @nest_index: (out): FIXME
 *
 * FIXME 
 * 
 */
void 
ncm_sphere_map_pix_ang2pix_nest (NcmSphereMapPix *pix, const gdouble theta, const gdouble phi, gint64 *nest_index)
{
  const gdouble z = cos (theta);
  if (theta > 0.1)
    _ncm_sphere_map_pix_zphi2pix_nest (pix, z, 1.0 - gsl_pow_2 (z), phi, nest_index);
  else
    _ncm_sphere_map_pix_zphi2pix_nest (pix, z, gsl_pow_2 (sin (theta)), phi, nest_index);
}

/**
 * ncm_sphere_map_pix_ang2pix_ring:
 * @pix: a #NcmSphereMapPix
 * @theta: FIXME
 * @phi: FIXME
 * @ring_index: (out): FIXME
 *
 * FIXME 
 * 
 */
void
ncm_sphere_map_pix_ang2pix_ring (NcmSphereMapPix *pix, const gdouble theta, const gdouble phi, gint64 *ring_index)
{
  const gdouble z = cos (theta);
  if (theta > 0.1)
    _ncm_sphere_map_pix_zphi2pix_ring (pix, z, 1.0 - gsl_pow_2 (z), phi, ring_index);
  else
    _ncm_sphere_map_pix_zphi2pix_ring (pix, z, gsl_pow_2 (sin (theta)), phi, ring_index);
}

/**
 * ncm_sphere_map_pix_vec2pix_ring:
 * @pix: a #NcmSphereMapPix
 * @vec: a #NcmTriVec
 * @ring_index: (out): FIXME 
 *
 * FIXME 
 * 
 */
void 
ncm_sphere_map_pix_vec2pix_ring (NcmSphereMapPix *pix, NcmTriVec *vec, gint64 *ring_index)
{
  const gdouble norm = ncm_trivec_norm (vec);
  const gdouble z    = vec->c[2] / norm;
  const gdouble phi  = ncm_trivec_get_phi (vec);
  _ncm_sphere_map_pix_zphi2pix_ring (pix, z, 1.0 - gsl_pow_2 (z), phi, ring_index);
}

/**
 * ncm_sphere_map_pix_vec2pix_nest:
 * @pix: a #NcmSphereMapPix
 * @vec: a #NcmTriVec
 * @nest_index: (out): FIXME 
 *
 * FIXME 
 * 
 */
void 
ncm_sphere_map_pix_vec2pix_nest (NcmSphereMapPix *pix, NcmTriVec *vec, gint64 *nest_index)
{
  const gdouble norm = ncm_trivec_norm (vec);
  const gdouble z    = vec->c[2] / norm;
  const gdouble phi  = ncm_trivec_get_phi (vec);
  _ncm_sphere_map_pix_zphi2pix_nest (pix, z, 1.0 - gsl_pow_2 (z), phi, nest_index);
}

/**
 * ncm_sphere_map_pix_add_to_vec:
 * @pix: a #NcmSphereMapPix
 * @vec: a #NcmTriVec
 * @s: signal
 *
 * FIXME
 * 
 */
void 
ncm_sphere_map_pix_add_to_vec (NcmSphereMapPix *pix, NcmTriVec *vec, const gdouble s)
{
  gint64 index = 0;
  switch (pix->order)
  {
    case NCM_SPHERE_MAP_PIX_ORDER_NEST:
      ncm_sphere_map_pix_vec2pix_nest (pix, vec, &index);
      break;
    case NCM_SPHERE_MAP_PIX_ORDER_RING:
      ncm_sphere_map_pix_vec2pix_ring (pix, vec, &index);
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  _fft_vec_idx (pix->pvec, index) += s;
}

/**
 * ncm_sphere_map_pix_add_to_ang:
 * @pix: a #NcmSphereMapPix
 * @theta: $\theta$
 * @phi: $\phi$
 * @s: signal
 *
 * FIXME
 * 
 */
void 
ncm_sphere_map_pix_add_to_ang (NcmSphereMapPix *pix, const gdouble theta, const gdouble phi, const gdouble s)
{
  gint64 index = 0;
  switch (pix->order)
  {
    case NCM_SPHERE_MAP_PIX_ORDER_NEST:
      ncm_sphere_map_pix_ang2pix_nest (pix, theta, phi, &index);
      break;
    case NCM_SPHERE_MAP_PIX_ORDER_RING:
      ncm_sphere_map_pix_ang2pix_ring (pix, theta, phi, &index);
      break;
    default:
      g_assert_not_reached ();
      break;
  }
  _fft_vec_idx (pix->pvec, index) += s;
}

/**
 * ncm_sphere_map_pix_load_fits:
 * @pix: a #NcmSphereMapPix
 * @fits_file: fits filename
 * @signal_name: (allow-none): signal column name in @fits_file
 *
 * FIXME 
 * 
 */
void
ncm_sphere_map_pix_load_fits (NcmSphereMapPix *pix, const gchar *fits_file, const gchar *signal_name)
{
  gchar comment[FLEN_COMMENT];
  gchar ordering[FLEN_VALUE];
  gchar coordsys[FLEN_VALUE];
  gint  status, hdutype, anynul;   
  glong nside, nfields, naxis2;
  gint signal_i = 0;
  const gchar *sname = signal_name != NULL ?  signal_name : NCM_SPHERE_MAP_PIX_DEFAULT_SIGNAL;
  fitsfile *fptr;

  status = 0;

  fits_open_file (&fptr, fits_file, READONLY, &status); 
  NCM_FITS_ERROR (status);
   
  fits_movabs_hdu (fptr, 2, &hdutype, &status); 
  NCM_FITS_ERROR (status);

  if (hdutype != BINARY_TBL) 
    g_error ("ncm_sphere_map_pix_load_fits: `%s' is not a binary table.", fits_file);
   
  fits_read_key_lng (fptr, "NSIDE", &nside, comment, &status); 
  NCM_FITS_ERROR(status);

  g_assert_cmpint (nside, >, 0);
  ncm_sphere_map_pix_set_nside (pix, nside);
  
  fits_read_key_lng (fptr, "TFIELDS", &nfields, comment, &status); 
  NCM_FITS_ERROR(status);

  g_assert_cmpint (nfields, >, 0);
  
  fits_read_key_lng (fptr, "NAXIS2", &naxis2, comment, &status); 
  NCM_FITS_ERROR(status);

  g_assert_cmpint (naxis2, ==, ncm_sphere_map_pix_get_npix (pix));

  if (fits_get_colnum (fptr, CASESEN, (gchar *)sname, &signal_i, &status))
    g_error ("ncm_sphere_map_pix_load_fits: signal column named `%s' not found in `%s'.",
             sname, fits_file);

#ifdef HAVE_FFTW3F
  fits_read_col_flt (fptr, signal_i, 1, 1, naxis2, NCM_SPHERE_MAP_PIX_HEALPIX_NULLVAL, 
                     _fft_vec_ptr (pix->pvec, 0), &anynul, &status); 
#elif defined (HAVE_FFTW3)
  fits_read_col_dbl (fptr, signal_i, 1, 1, naxis2, NCM_SPHERE_MAP_PIX_HEALPIX_NULLVAL, 
                     _fft_vec_ptr (pix->pvec, 0), &anynul, &status); 
#else
  fits_read_col_flt (fptr, signal_i, 1, 1, naxis2, NCM_SPHERE_MAP_PIX_HEALPIX_NULLVAL, 
                     _fft_vec_ptr (pix->pvec, 0), &anynul, &status); 
#endif  
  NCM_FITS_ERROR (status);

  if (fits_read_key (fptr, TSTRING, "ORDERING", ordering, comment, &status)) 
  {
    g_warning ("ncm_sphere_map_pix_load_fits: Could not find ORDERING in the fits file, assuming RING.");
    status = 0;
  }

  switch (*ordering)
  {
    case 'N':
      pix->order = NCM_SPHERE_MAP_PIX_ORDER_NEST;
      break;
    case 'R':
      pix->order = NCM_SPHERE_MAP_PIX_ORDER_RING;
      break;
    default:
      g_error ("ncm_sphere_map_pix_load_fits: Unknown order type `%s'.", ordering);
      break;
  }  

  if (fits_read_key (fptr, TSTRING, "COORDSYS", coordsys, comment, &status)) 
  {
    g_warning ("ncm_sphere_map_pix_load_fits: Could not find COORDSYS in the fits file, assuming CELESTIAL.");
    coordsys[0] = NCM_SPHERE_MAP_PIX_COORD_SYS_CELESTIAL;
    coordsys[1] = '\0';
    status = 0;
  }

  ncm_sphere_map_pix_set_coordsys (pix, *coordsys);
  
  fits_close_file (fptr, &status);
  NCM_FITS_ERROR (status);
}

/**
 * ncm_sphere_map_pix_save_fits:
 * @pix: a #NcmSphereMapPix
 * @fits_file: fits filename
 * @signal_name: (allow-none): signal column name in @fits_file
 * @overwrite: FIXME
 *
 * FIXME 
 * 
 */
void
ncm_sphere_map_pix_save_fits (NcmSphereMapPix *pix, const gchar *fits_file, const gchar *signal_name, gboolean overwrite)
{    
  const gchar *sname = signal_name != NULL ?  signal_name : NCM_SPHERE_MAP_PIX_DEFAULT_SIGNAL;
  const gchar *ttype[] = { sname };
  const gchar *tform[] = { "1E" };
  const gint64 npix    = ncm_sphere_map_pix_get_npix (pix);
  const gchar extname[] = "BINTABLE";  
  fitsfile *fptr;
  gint status = 0;
 
  if (overwrite && g_file_test (fits_file, G_FILE_TEST_EXISTS))
    g_unlink (fits_file);
  
  fits_create_file (&fptr, fits_file, &status);
  NCM_FITS_ERROR (status);

  fits_create_tbl (fptr, BINARY_TBL, npix, 1, (gchar **)ttype, (gchar **)tform, NULL, extname, &status);
  NCM_FITS_ERROR (status);

  {
    const gchar *order_str = ncm_sphere_map_pix_get_order (pix) == NCM_SPHERE_MAP_PIX_ORDER_NEST ? "NESTED  " : "RING    ";
    fits_write_key (fptr, TSTRING, "ORDERING", (gchar *) order_str,
                    "Pixel ordering scheme, either RING or NESTED", &status);
    NCM_FITS_ERROR (status);
  }

  {
    glong firstpix = 0;
    fits_write_key (fptr, TLONG, "FIRSTPIX", &firstpix,
                    "First pixel # (0 based)", &status);
    NCM_FITS_ERROR (status);
  }

  {
    glong lastpix = npix - 1;
    fits_write_key (fptr, TLONG, "LASTPIX", &lastpix,
                    "Last pixel # (0 based)", &status);
    NCM_FITS_ERROR (status);
  }

  {
    fits_write_key (fptr, TSTRING, "INDXSCHM", (gchar *) "IMPLICIT",
                    "Indexing: IMPLICIT or EXPLICIT", &status);
    NCM_FITS_ERROR (status);
  }
  
  {
    fits_write_key (fptr, TSTRING, "INDXSCHM", (gchar *) "IMPLICIT",
                    "Indexing: IMPLICIT or EXPLICIT", &status);
    NCM_FITS_ERROR (status);
  }

  {
    glong nside = ncm_sphere_map_pix_get_nside (pix);
    fits_write_key (fptr, TLONG, "NSIDE", &nside,
                    "Resolution parameter for HEALPIX", &status);
    NCM_FITS_ERROR (status);
  }

  {
    gchar coordsys_c = (gchar)ncm_sphere_map_pix_get_coordsys (pix);
    gchar *coordsys  = g_strdup_printf ("%c       ", coordsys_c);
    
    fits_write_key(fptr, TSTRING, "COORDSYS", coordsys,
                   "Pixelisation coordinate system", &status);
    NCM_FITS_ERROR (status);
    
    fits_write_comment(fptr,
                       "G = Galactic, E = ecliptic, C = celestial = equatorial", &status);
    NCM_FITS_ERROR (status);

    g_free (coordsys);
  }

#ifdef HAVE_FFTW3F
  fits_write_col (fptr, TFLOAT, 1, 1, 1, npix, _fft_vec_ptr (pix->pvec, 0), &status);
  NCM_FITS_ERROR (status);
#elif defined (HAVE_FFTW3)
  {
    gint64 i;
    for (i = 0; i < npix; i++)
    {
      gfloat v_i = _fft_vec_idx (pix->pvec, i);
      fits_write_col (fptr, TFLOAT, 1, 1 + i, 1, 1, &v_i, &status);
      NCM_FITS_ERROR (status);
    }
  }
#else
  fits_write_col (fptr, TFLOAT, 1, 1, 1, npix, _fft_vec_ptr (pix->pvec, 0), &status);
  NCM_FITS_ERROR (status);
#endif


  fits_close_file (fptr, &status);
  NCM_FITS_ERROR (status);
}

static void 
_ncm_sphere_map_pix_radec_to_ang (const gdouble RA, const gdouble DEC, gdouble *theta, gdouble *phi)
{
  theta[0] = gsl_sf_angle_restrict_pos (ncm_c_degree_to_radian (90.0 - DEC));
  phi[0]   = gsl_sf_angle_restrict_pos (ncm_c_degree_to_radian (RA));
}

/**
 * ncm_sphere_map_pix_load_from_fits_catalog:
 * @pix: a #NcmSphereMapPix
 * @fits_file: fits filename
 * @RA: RA column name in @fits_file
 * @DEC: DEC column name in @fits_file
 * @S: (allow-none): Signal column name in @fits_file
 *
 * FIXME 
 * 
 */
void 
ncm_sphere_map_pix_load_from_fits_catalog (NcmSphereMapPix *pix, const gchar *fits_file, const gchar *RA, const gchar *DEC, const gchar *S)
{
  gchar comment[FLEN_COMMENT];
  gint  status, hdutype, anynul;   
  glong naxis2;
  gint RA_col  = 0;
  gint DEC_col = 0;
  gint S_col   = 0;
  fitsfile *fptr;
  gint64 i;

  status = 0;

  fits_open_file (&fptr, fits_file, READONLY, &status); 
  NCM_FITS_ERROR (status);
   
  fits_movabs_hdu (fptr, 2, &hdutype, &status); 
  NCM_FITS_ERROR (status);

  if (hdutype != BINARY_TBL) 
    g_error ("ncm_sphere_map_pix_load_fits: `%s' is not a binary table.", fits_file);

  if (fits_get_colnum (fptr, CASESEN, (gchar *)RA, &RA_col, &status))
    g_error ("ncm_sphere_map_pix_load_fits: RA column named `%s' not found in `%s'.",
             RA, fits_file);

  if (fits_get_colnum (fptr, CASESEN, (gchar *)DEC, &DEC_col, &status))
    g_error ("ncm_sphere_map_pix_load_fits: DEC column named `%s' not found in `%s'.",
             DEC, fits_file);

  if (S != NULL)
  {
    if (fits_get_colnum (fptr, CASESEN, (gchar *)S, &S_col, &status))
      g_error ("ncm_sphere_map_pix_load_fits: Signal column named `%s' not found in `%s'.",
               S, fits_file);
  }

  fits_read_key_lng (fptr, "NAXIS2", &naxis2, comment, &status); 
  NCM_FITS_ERROR(status);

  if (S != NULL)
  {
    for (i = 0; i < naxis2; i++)
    {
      gdouble RA_i, DEC_i, theta, phi;
      gdouble S_i;

      fits_read_col_dbl (fptr, RA_col, 1 + i, 1, 1, NCM_SPHERE_MAP_PIX_HEALPIX_NULLVAL, 
                         &RA_i, &anynul, &status); 
      NCM_FITS_ERROR (status);

      fits_read_col_dbl (fptr, DEC_col, 1 + i, 1, 1, NCM_SPHERE_MAP_PIX_HEALPIX_NULLVAL, 
                         &DEC_i, &anynul, &status); 
      NCM_FITS_ERROR (status);

      fits_read_col_dbl (fptr, S_col, 1 + i, 1, 1, NCM_SPHERE_MAP_PIX_HEALPIX_NULLVAL, 
                         &S_i, &anynul, &status); 
      NCM_FITS_ERROR (status);

      _ncm_sphere_map_pix_radec_to_ang (RA_i, DEC_i, &theta, &phi);

      ncm_sphere_map_pix_add_to_ang (pix, theta, phi, S_i);
    }
  }
  else
  {
    for (i = 0; i < naxis2; i++)
    {
      gdouble RA_i, DEC_i, theta, phi;

      fits_read_col_dbl (fptr, RA_col, 1 + i, 1, 1, NCM_SPHERE_MAP_PIX_HEALPIX_NULLVAL, 
                         &RA_i, &anynul, &status); 
      NCM_FITS_ERROR (status);

      fits_read_col_dbl (fptr, DEC_col, 1 + i, 1, 1, NCM_SPHERE_MAP_PIX_HEALPIX_NULLVAL, 
                         &DEC_i, &anynul, &status); 
      NCM_FITS_ERROR (status);

      _ncm_sphere_map_pix_radec_to_ang (RA_i, DEC_i, &theta, &phi);

      ncm_sphere_map_pix_add_to_ang (pix, theta, phi, 1.0);
    }
  }
}

static void
_ncm_sphere_map_pix_prepare_fft (NcmSphereMapPix *pix)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  if (pix->fft_plan->len == 0)
  {
    const gint64 npix      = ncm_sphere_map_pix_get_npix (pix);
    const gint64 nring_cap = ncm_sphere_map_pix_get_nrings_cap (pix);
    gint r_i;

    ncm_cfg_load_fftw_wisdom ("ncm_sphere_map_pix_nside_%ld", ncm_sphere_map_pix_get_nside (pix));
#  ifdef HAVE_FFTW3F

    pix->fft_pvec = fftwf_alloc_complex (npix);

    for (r_i = 0; r_i < nring_cap; r_i++)
    {
      const gint ring_size       = ncm_sphere_map_pix_get_ring_size (pix, r_i);
      const gint64 ring_fi_north = ncm_sphere_map_pix_get_ring_first_index (pix, r_i);
      const gint64 ring_fi_south = ncm_sphere_map_pix_get_ring_first_index (pix, ncm_sphere_map_pix_get_nrings (pix) - r_i - 1);
      const gint64 dist          = ring_fi_south - ring_fi_north;
      gfloat *pvec               = pix->pvec;
      complex float *fft_pvec    = pix->fft_pvec;

      
      fftwf_plan plan = fftwf_plan_many_dft_r2c (1, &ring_size, 2,
                                                 &pvec[ring_fi_north], NULL, 
                                                 1, dist,
                                                 &fft_pvec[ring_fi_north], NULL,
                                                 1, dist,
                                                 fftw_default_flags | FFTW_PRESERVE_INPUT);
printf ("Preparing plan for %ld and %ld size %d | npix %ld | %p\n", ring_fi_north, ring_fi_south, ring_size, pix->npix, plan);
      g_ptr_array_add (pix->fft_plan, plan);
    }
    {
      const gint ring_size     = pix->middle_rings_size;
      const gint nrings_middle = ncm_sphere_map_pix_get_nrings_middle (pix);
      const gint cap_size      = ncm_sphere_map_pix_get_cap_size (pix);

      gfloat *pvec               = pix->pvec;
      complex float *fft_pvec    = pix->fft_pvec;


      fftwf_plan plan = fftwf_plan_many_dft_r2c (1, &ring_size, nrings_middle,
                                                 &pvec[cap_size], NULL, 
                                                 1, ring_size,
                                                 &fft_pvec[cap_size], NULL,
                                                 1, ring_size,
                                                 fftw_default_flags | FFTW_PRESERVE_INPUT);
printf ("Preparing plan for %d and %d x %d | npix %ld | %p\n", cap_size, nrings_middle, ring_size, pix->npix, plan);
      g_ptr_array_add (pix->fft_plan, plan);
    }  

#  else

    pix->fft_pvec = fftw_alloc_complex (npix);

    for (r_i = 0; r_i < nring_cap; r_i++)
    {
      const gint ring_size       = ncm_sphere_map_pix_get_ring_size (pix, r_i);
      const gint64 ring_fi_north = ncm_sphere_map_pix_get_ring_first_index (pix, r_i);
      const gint64 ring_fi_south = ncm_sphere_map_pix_get_ring_first_index (pix, ncm_sphere_map_pix_get_nrings (pix) - r_i - 1);
      const gint64 dist          = ring_fi_south - ring_fi_north;
      gdouble *pvec              = pix->pvec;
      complex double *fft_pvec   = pix->fft_pvec;

      fftw_plan plan = fftw_plan_many_dft_r2c (1, &ring_size, 2,
                                               &pvec[ring_fi_north], NULL, 
                                               1, dist,
                                               &fft_pvec[ring_fi_north], NULL,
                                               1, dist,
                                               fftw_default_flags | FFTW_DESTROY_INPUT);
      g_ptr_array_add (pix->fft_plan, plan);
    }
    {
      const gint ring_size     = pix->middle_rings_size;
      const gint nrings_middle = ncm_sphere_map_pix_get_nrings_middle (pix);
      const gint cap_size      = ncm_sphere_map_pix_get_cap_size (pix);

      gdouble *pvec               = pix->pvec;
      complex double *fft_pvec    = pix->fft_pvec;

      fftw_plan plan = fftw_plan_many_dft_r2c (1, &ring_size, nrings_middle,
                                               &pvec[cap_size], NULL, 
                                               1, ring_size,
                                               &fft_pvec[cap_size], NULL,
                                               1, ring_size,
                                               fftw_default_flags | FFTW_DESTROY_INPUT);
      g_ptr_array_add (pix->fft_plan, plan);
    }    
#  endif
    ncm_cfg_save_fftw_wisdom ("ncm_sphere_map_pix_nside_%ld", ncm_sphere_map_pix_get_nside (pix));
  }
#endif
}

#ifdef NUMCOSMO_HAVE_FFTW3
static void
_ncm_sphere_map_pix_get_alm_from_circle (NcmSphereMapPix *pix, gint64 r_i, gint64 lmax, _fft_complex *alm)
{
  const gint64 ring_size  = ncm_sphere_map_pix_get_ring_size (pix, r_i);
  const gint64 ring_fi    = ncm_sphere_map_pix_get_ring_first_index (pix, r_i);
  const _fft_complex *Fim = &((_fft_complex *)pix->fft_pvec)[ring_fi];
  const gdouble pix_area  = 4.0 * M_PI / pix->npix;
  gsize Ylm_size  = gsl_sf_legendre_array_n (lmax);
  static gdouble *Ylm    = NULL;
  gdouble theta_i = 0.0, phi_0 = 0.0, x = 0.0;
  gint64 m;

  ncm_sphere_map_pix_pix2ang_ring (pix, ring_fi, &theta_i, &phi_0);

  x = cos (theta_i);

  if (Ylm == NULL)
    Ylm = g_new (gdouble, Ylm_size);

  gsl_sf_legendre_array_e (GSL_SF_LEGENDRE_SPHARM, lmax, x, -1.0, Ylm);
  
  for (m = 0; m <= lmax; m++)
  {
    gint64 l;
    gint fft_ring_index = m % ring_size;

    _fft_complex phase = cexp (-I * phi_0 * m);
    
    if (fft_ring_index <= ring_size / 2)
      phase *= Fim[fft_ring_index];
    else
      phase *= conj (Fim[ring_size - fft_ring_index]);

    for (l = m; l <= lmax; l++)
    {
      gsize lm_index = gsl_sf_legendre_array_index (l, m);
      alm[lm_index] += Ylm[lm_index] * phase * pix_area;

      /*printf ("Ylm[%ld, %ld](% 20.15g) = % 20.15g\n", l, m, x, Ylm[lm_index]);*/
      
    }
  }
}
#endif

/**
 * ncm_sphere_map_pix_get_alm:
 * @pix: a #NcmSphereMapPix
 *
 * FIXME 
 * 
 */
void 
ncm_sphere_map_pix_get_alm (NcmSphereMapPix *pix)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  gint i;
  _ncm_sphere_map_pix_prepare_fft (pix);

  ncm_sphere_map_pix_set_order (pix, NCM_SPHERE_MAP_PIX_ORDER_RING);
  
  for (i = 0; i < pix->fft_plan->len; i++)
  {
#  ifdef HAVE_FFTW3F
    fftwf_execute (g_ptr_array_index (pix->fft_plan, i));    
#  else
    fftw_execute (g_ptr_array_index (pix->fft_plan, i));
#endif
  }

  if (FALSE)
  {
    gint64 r_i;
    for (r_i = 0; r_i < pix->nrings; r_i++)
    {
      gint64 r_fi = ncm_sphere_map_pix_get_ring_first_index (pix, r_i);
      gint64 r_s  = ncm_sphere_map_pix_get_ring_size (pix, r_i);
      gint64 ir_i;
      for (ir_i = 0; ir_i < r_s; ir_i++)
      {
        if (ir_i <= r_s / 2)
        {
          gint64 j = r_fi + ir_i;
          printf ("r_i %4ld ir_i %4ld p_i %4ld IN % 20.15g OUT % 20.15g % 20.15g\n", 
                  r_i, ir_i, r_fi + ir_i, 
                  ((gfloat *)pix->pvec)[j], 
                  crealf (((complex float *)pix->fft_pvec)[j]), cimagf (((complex float *)pix->fft_pvec)[j]));
        }
        else
        {
          gint64 j = r_fi + (r_s - ir_i);
          printf ("r_i %4ld ir_i %4ld p_i %4ld IN % 20.15g OUT % 20.15g % 20.15g\n", 
                  r_i, ir_i, r_fi + ir_i, 
                  ((gfloat *)pix->pvec)[j], 
                  crealf (((complex float *)pix->fft_pvec)[j]), -cimagf (((complex float *)pix->fft_pvec)[j]));
        }
      }
    }
  }
 
  {
    const gint64 lmax   = 3 * pix->nside - 1;
    const gint64 nrings = ncm_sphere_map_pix_get_nrings (pix);
    gint64 r_i;
    
    _fft_complex *alm = _fft_vec_alloc_complex (NCM_SPHERE_MAP_PIX_ALM_SIZE (lmax));
    memset (alm, 0, sizeof (_fft_complex) * NCM_SPHERE_MAP_PIX_ALM_SIZE (lmax));
    
    for (r_i = 0; r_i < nrings; r_i++)
    {
      _ncm_sphere_map_pix_get_alm_from_circle (pix, r_i, lmax, alm);
    }

    {
      gint64 l, m;

      for (m = 0; m <= lmax; m++)
      {
        for (l = m; l <= lmax; l++)
        {
          gsize lm_index = gsl_sf_legendre_array_index (l, m); 
          printf ("alm[%.4ld, %.4ld] = % 20.15g + i % 20.15g\n", l, m, creal (alm[lm_index]), cimag (alm[lm_index]));
        }
      }
    }
  }
  
  
#else
  g_error ("ncm_sphere_map_pix_get_alm: no fftw3 support, to use this function recompile NumCosmo with fftw.");
#endif
}

