/***************************************************************************
 *            ncm_sphere_map.c
 *
 *  Wed Jul  9 11:09:37 2008 (updated Jul 15 2016)
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sphere_map.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcmSphereMap:
 *
 * Maps on the sphere in the HEALPix pixelization.
 *
 * An independent implementation of the HEALPix pixelization of the sphere, written from
 * its definition in [Górski et al. (2005)](https://arxiv.org/abs/astro-ph/0409513) and
 * not derived from the HEALPix code. For a resolution parameter nside, a power of two,
 * the sphere is divided into $12\,\mathrm{nside}^2$ pixels of equal area on
 * $4\,\mathrm{nside} - 1$ rings of constant latitude. Pixels are numbered in RING or
 * NESTED order (#NcmSphereMapOrder), with the same indices, centres and pixel lookups as
 * HEALPix; the tests compare them with healpy.
 *
 * The object holds a map of values on those pixels, reads and writes it in the HEALPix
 * FITS format, and computes its spherical harmonic coefficients $a_{\ell m}$ and power
 * spectrum $C_\ell$ up to #NcmSphereMap:lmax, and a map from given $a_{\ell m}$.
 *
 */
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/sphere/ncm_sphere_map.h"
#include "ncm/algebra/ncm_vector.h"
#include "ncm/specfunc/ncm_sf_spherical_harmonics.h"
#include "ncm/core/ncm_util.h"
#include "ncm/algebra/ncm_complex.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_timer.h"
#include "ncm/spline/ncm_spline_func.h"
#include "ncm/spline/ncm_spline_cubic_notaknot.h"
#include "ncm_enum_types.h"

/*#define _NCM_SPHERE_MAP_MEASURE 1*/

#ifndef NUMCOSMO_GIR_SCAN
#ifdef HAVE_CFITSIO
#include <fitsio.h>
#endif /* HAVE_CFITSIO */

#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

#  define _fft_vec_alloc fftw_alloc_real
#  define _fft_real gdouble
#  define _fft_complex complex double
#  define _fft_vec_alloc_complex fftw_alloc_complex
#  define _fft_vec_free  fftw_free
#  define _fft_vec_set_zero(v, s) memset ((v), 0, sizeof (gdouble) * (s))
#  define _fft_vec_set_zero_complex(v, s) memset ((v), 0, sizeof (_fft_complex) * (s))
#  define _fft_vec_memcpy(dest, orig, s) memcpy ((dest), (orig), sizeof (gdouble) * (s))
#  define _fft_vec_ptr(v, i) (&((gdouble *) (v))[i])
#define _fft_vec_idx(v, i) (*_fft_vec_ptr (v, i))

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_trig.h>
#endif /* NUMCOSMO_GIR_SCAN */

#ifndef NCM_SPHERE_MAP_BLOCK_NC
#define NCM_SPHERE_MAP_BLOCK_NC 2
#define NCM_SPHERE_MAP_BLOCK_STEP 8
#define NCM_SPHERE_MAP_BLOCK_CM 32
#endif

#define NCM_SPHERE_MAP_BLOCK_NCT (2 * NCM_SPHERE_MAP_BLOCK_NC)
#define NCM_SPHERE_MAP_BLOCK_STEPM2 (NCM_SPHERE_MAP_BLOCK_STEP - 2)

#define NCM_SPHERE_MAP_BLOCK_INV_NC 3
#define NCM_SPHERE_MAP_BLOCK_INV_NCT (2 * NCM_SPHERE_MAP_BLOCK_INV_NC)
#define NCM_SPHERE_MAP_BLOCK_INV_STEP 10
#define NCM_SPHERE_MAP_BLOCK_INV_STEPM2 (NCM_SPHERE_MAP_BLOCK_INV_STEP - 2)

#define NCM_SPHERE_MAP_BLOCK_DEC_XX(name, S, N) name ## _ ## S ## _ ## N
#define NCM_SPHERE_MAP_BLOCK_DEC_X(name, S, N) NCM_SPHERE_MAP_BLOCK_DEC_XX (name, S, N)

#define NCM_SPHERE_MAP_BLOCK_DEC(name) NCM_SPHERE_MAP_BLOCK_DEC_X (name, NCM_SPHERE_MAP_BLOCK_STEP, NCM_SPHERE_MAP_BLOCK_NC)
#define NCM_SPHERE_MAP_BLOCK_INV_DEC(name) NCM_SPHERE_MAP_BLOCK_DEC_X (name, NCM_SPHERE_MAP_BLOCK_INV_STEP, NCM_SPHERE_MAP_BLOCK_INV_NC)

typedef struct _NcmSphereMapPrivate
{
  gint64 nside;
  gint64 npix;
  gint64 face_size;
  gint64 middle_rings_size;
  gint64 cap_size;
  gint64 middle_size;
  gint64 nrings;
  gint64 nrings_cap;
  gint64 nrings_middle;
  gint64 block_ring_size;
  gint64 last_sing_ring;
  NcmSphereMapOrder order;
  NcmSphereMapCoordSys coordsys;
  gpointer pvec;
  gpointer fft_pvec;
  GPtrArray *fft_plan_r2c;
  GPtrArray *fft_plan_c2r;
  guint lmax;
  guint iter;
  _fft_complex *alm;
  gint64 alm_len;
  NcmVector *alm_v;
  NcmVector *Cl;
  gboolean has_Cls;
  NcmTimer *t;
  NcmSFSphericalHarmonics *spha;
  GPtrArray *sphaY_array;
  GPtrArray *sphaYa_array;
  GArray *block_data;
} NcmSphereMapPrivate;

typedef struct _NcmSphereMapBlock
{
  _fft_complex * restrict Fima[NCM_SPHERE_MAP_BLOCK_NCT];
  gint64 ring_size[NCM_SPHERE_MAP_BLOCK_NCT];
  gint64 ring_size_2[NCM_SPHERE_MAP_BLOCK_NCT];
  gdouble theta[NCM_SPHERE_MAP_BLOCK_NCT];
  gdouble phi[NCM_SPHERE_MAP_BLOCK_NCT];
} NcmSphereMapBlock;

enum
{
  PROP_0,
  PROP_NSIDE,
  PROP_ORDER,
  PROP_COORDSYS,
  PROP_LMAX,
  PROP_ITER,
};

struct _NcmSphereMap
{
  GObject parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmSphereMap, ncm_sphere_map, G_TYPE_OBJECT)

static void
ncm_sphere_map_init (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  self->nside             = 0;
  self->npix              = 0;
  self->face_size         = 0;
  self->middle_rings_size = 0;
  self->cap_size          = 0;
  self->middle_size       = 0;
  self->nrings            = 0;
  self->nrings_cap        = 0;
  self->nrings_middle     = 0;
  self->block_ring_size   = 0;
  self->last_sing_ring    = 0;
  self->order             = NCM_SPHERE_MAP_ORDER_RING;
  self->coordsys          = NCM_SPHERE_MAP_COORD_SYS_LEN;
  self->pvec              = NULL;
  self->fft_pvec          = NULL;
  self->fft_plan_r2c      = g_ptr_array_new ();
  self->fft_plan_c2r      = g_ptr_array_new ();
  g_ptr_array_set_free_func (self->fft_plan_r2c, ncm_cfg_fftw_plan_destroy);
  g_ptr_array_set_free_func (self->fft_plan_c2r, ncm_cfg_fftw_plan_destroy);
  self->alm          = NULL;
  self->alm_len      = 0;
  self->lmax         = 0;
  self->iter         = 0;
  self->Cl           = NULL;
  self->has_Cls      = FALSE;
  self->t            = ncm_timer_new ();
  self->spha         = ncm_sf_spherical_harmonics_new (1 << 12);
  self->sphaY_array  = g_ptr_array_new ();
  self->sphaYa_array = g_ptr_array_new ();
  self->block_data   = g_array_new (FALSE, FALSE, sizeof (NcmSphereMapBlock));

  g_ptr_array_set_free_func (self->sphaY_array,  (GDestroyNotify) ncm_sf_spherical_harmonics_Y_free);
  g_ptr_array_set_free_func (self->sphaYa_array, (GDestroyNotify) ncm_sf_spherical_harmonics_Y_array_free);
}

static void
_ncm_sphere_map_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSphereMap *smap = NCM_SPHERE_MAP (object);

  g_return_if_fail (NCM_IS_SPHERE_MAP (object));

  switch (prop_id)
  {
    case PROP_NSIDE:
      ncm_sphere_map_set_nside (smap, g_value_get_int64 (value));
      break;
    case PROP_ORDER:
      ncm_sphere_map_set_order (smap, g_value_get_enum (value));
      break;
    case PROP_COORDSYS:
      ncm_sphere_map_set_coordsys (smap, g_value_get_enum (value));
      break;
    case PROP_LMAX:
      ncm_sphere_map_set_lmax (smap, g_value_get_uint (value));
      break;
    case PROP_ITER:
      ncm_sphere_map_set_iter (smap, g_value_get_uint (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sphere_map_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSphereMap *smap               = NCM_SPHERE_MAP (object);
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  g_return_if_fail (NCM_IS_SPHERE_MAP (object));

  switch (prop_id)
  {
    case PROP_NSIDE:
      g_value_set_int64 (value, self->nside);
      break;
    case PROP_ORDER:
      g_value_set_enum (value, ncm_sphere_map_get_order (smap));
      break;
    case PROP_COORDSYS:
      g_value_set_enum (value, ncm_sphere_map_get_coordsys (smap));
      break;
    case PROP_LMAX:
      g_value_set_uint (value, ncm_sphere_map_get_lmax (smap));
      break;
    case PROP_ITER:
      g_value_set_uint (value, ncm_sphere_map_get_iter (smap));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sphere_map_dispose (GObject *object)
{
  NcmSphereMap *smap               = NCM_SPHERE_MAP (object);
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  g_clear_pointer (&self->alm,  _fft_vec_free);
  ncm_vector_clear (&self->Cl);

  ncm_timer_clear (&self->t);
  ncm_sf_spherical_harmonics_clear (&self->spha);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sphere_map_parent_class)->dispose (object);
}

static void
_ncm_sphere_map_finalize (GObject *object)
{
  NcmSphereMap *smap               = NCM_SPHERE_MAP (object);
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  ncm_sphere_map_set_nside (smap, 0);
  g_ptr_array_unref (self->fft_plan_r2c);
  g_ptr_array_unref (self->fft_plan_c2r);

  g_clear_pointer (&self->alm,  _fft_vec_free);
  ncm_vector_clear (&self->Cl);

  g_clear_pointer (&self->sphaY_array,  g_ptr_array_unref);
  g_clear_pointer (&self->sphaYa_array, g_ptr_array_unref);
  g_clear_pointer (&self->block_data,   g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sphere_map_parent_class)->finalize (object);
}

static void
ncm_sphere_map_class_init (NcmSphereMapClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_sphere_map_set_property;
  object_class->get_property = &_ncm_sphere_map_get_property;
  object_class->dispose      = &_ncm_sphere_map_dispose;
  object_class->finalize     = &_ncm_sphere_map_finalize;

  /**
   * NcmSphereMap:nside:
   *
   * The HEALPix resolution parameter, a power of two, or zero for no pixels.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_NSIDE,
                                   g_param_spec_int64 ("nside",
                                                       NULL,
                                                       "nside",
                                                       0, G_MAXINT64, 0,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSphereMap:order:
   *
   * The pixel ordering; changing it reorders the stored map.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ORDER,
                                   g_param_spec_enum ("order",
                                                      NULL,
                                                      "Map pixel ordering",
                                                      NCM_TYPE_SPHERE_MAP_ORDER, NCM_SPHERE_MAP_ORDER_RING,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSphereMap:coordsys:
   *
   * The coordinate system the map is given in, see ncm_sphere_map_set_coordsys().
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_COORDSYS,
                                   g_param_spec_enum ("coordsys",
                                                      NULL,
                                                      "Map coordinate system",
                                                      NCM_TYPE_SPHERE_MAP_COORD_SYS, NCM_SPHERE_MAP_COORD_SYS_CELESTIAL,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSphereMap:lmax:
   *
   * The largest $\ell$ of the harmonic transforms, see ncm_sphere_map_set_lmax().
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LMAX,
                                   g_param_spec_uint ("lmax",
                                                      NULL,
                                                      "max ell",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmSphereMap:iter:
   *
   * The number of refinement iterations of the map to $a_{\ell m}$ transform, see
   * ncm_sphere_map_set_iter().
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ITER,
                                   g_param_spec_uint ("iter",
                                                      NULL,
                                                      "number of iterations",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_sphere_map_new:
 * @nside: the healpix nside parameter
 *
 * Creates a new #NcmSphereMap for a given @nside.
 *
 * Returns: (transfer full): a new #NcmSphereMap.
 */
NcmSphereMap *
ncm_sphere_map_new (const gint64 nside)
{
  NcmSphereMap *smap = g_object_new (NCM_TYPE_SPHERE_MAP,
                                     "nside", nside,
                                     NULL);

  return smap;
}

/**
 * ncm_sphere_map_ref:
 * @smap: a #NcmSphereMap
 *
 * Increases the reference count of @smap.
 *
 * Returns: (transfer full): @smap.
 */
NcmSphereMap *
ncm_sphere_map_ref (NcmSphereMap *smap)
{
  return g_object_ref (smap);
}

/**
 * ncm_sphere_map_free:
 * @smap: a #NcmSphereMap
 *
 * Decreases the reference count of @smap. When its reference count
 * drops to 0, the object is finalized (i.e. its memory is freed).
 *
 */
void
ncm_sphere_map_free (NcmSphereMap *smap)
{
  g_object_unref (smap);
}

/**
 * ncm_sphere_map_clear:
 * @smap: a #NcmSphereMap
 *
 * If *@smap is not %NULL, decreases the reference count of @smap.
 * When its reference count drops to 0, the object is finalized
 * (i.e. its memory is freed).
 * Set *@smap to %NULL.
 *
 */
void
ncm_sphere_map_clear (NcmSphereMap **smap)
{
  g_clear_object (smap);
}

static gint64
_l_pow_2 (gint64 n)
{
  return n * n;
}

static void
_ncm_sphere_map_check_index (const gint64 index, const gint64 size, const gchar *func, const gchar *what)
{
  if ((index < 0) || (index >= size))
    g_error ("%s: %s index %" G_GINT64_FORMAT " out of range [0, %" G_GINT64_FORMAT ").", func, what, index, size);
}

static void
_ncm_sphere_map_prepare_circle (NcmSphereMap *smap, NcmSphereMapBlock *block, const gint64 r_i, const gint64 i)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  const gint64 ring_fi             = ncm_sphere_map_get_ring_first_index (smap, r_i);

  block->ring_size[i]   = ncm_sphere_map_get_ring_size (smap, r_i);
  block->ring_size_2[i] = block->ring_size[i] / 2;

  block->Fima[i] = &((_fft_complex *) self->fft_pvec)[ring_fi];

  ncm_sphere_map_pix2ang_ring (smap, ring_fi, &block->theta[i], &block->phi[i]);
}

/**
 * ncm_sphere_map_set_nside:
 * @smap: a #NcmSphereMap
 * @nside: the healpix nside parameter
 *
 * Sets the nside parameter of @smap, which must be a power of two, or zero to release
 * the pixels. A change allocates a new map, with every pixel zero.
 *
 */
void
ncm_sphere_map_set_nside (NcmSphereMap *smap, gint64 nside)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  if ((nside < 0) || ((nside & (nside - 1)) != 0))
    g_error ("ncm_sphere_map_set_nside: nside must be zero or a power of two, got %" G_GINT64_FORMAT ".", nside);

  if (nside != self->nside)
  {
    g_clear_pointer (&self->pvec, (GDestroyNotify) _fft_vec_free);

    self->nside             = nside;
    self->npix              = 0;
    self->face_size         = 0;
    self->middle_rings_size = 0;
    self->cap_size          = 0;
    self->middle_size       = 0;
    self->nrings            = 0;
    self->nrings_cap        = 0;
    self->nrings_middle     = 0;
    self->block_ring_size   = 0;
    self->last_sing_ring    = 0;

    g_clear_pointer (&self->fft_pvec, (GDestroyNotify) fftw_free);

    g_ptr_array_set_size (self->fft_plan_r2c, 0);
    g_ptr_array_set_size (self->fft_plan_c2r, 0);

    g_ptr_array_set_size (self->sphaY_array,  0);
    g_ptr_array_set_size (self->sphaYa_array, 0);

    if (nside > 0)
    {
      self->npix              = NCM_SPHERE_MAP_N (self->nside);
      self->face_size         = _l_pow_2 (self->nside);
      self->middle_rings_size = 4 * nside;
      self->cap_size          = 2 * nside * (nside - 1);
      self->nrings            = 4 * nside - 1;
      self->nrings_cap        = nside - 1;
      self->nrings_middle     = 2 * nside + 1;
      self->middle_size       = self->nrings_middle * self->middle_rings_size;
      self->pvec              = _fft_vec_alloc (self->npix);

      _fft_vec_set_zero (self->pvec, self->npix);

      self->fft_pvec = fftw_alloc_complex (self->npix);

      g_assert_cmpint (2 * self->cap_size + self->middle_size, ==, self->npix);

      {
        const gint64 nrings_2 = self->nrings / 2;
        const gint64 nblocks  = nrings_2 / NCM_SPHERE_MAP_BLOCK_NC;
        gint64 i, j;

        self->block_ring_size = nblocks * NCM_SPHERE_MAP_BLOCK_NC;
        self->last_sing_ring  = self->nrings - self->block_ring_size;

        g_array_set_size (self->block_data, nblocks);

        for (j = 0; j < nblocks; j++)
        {
          NcmSFSphericalHarmonicsYArray *sphaYa = ncm_sf_spherical_harmonics_Y_array_new (self->spha, NCM_SPHERE_MAP_BLOCK_NCT, NCM_SF_SPHERICAL_HARMONICS_ARRAY_DEFAULT_ABSTOL);
          NcmSphereMapBlock *block              = &g_array_index (self->block_data, NcmSphereMapBlock, j);
          const gint64 r_ini                    = j * NCM_SPHERE_MAP_BLOCK_NC;

          g_ptr_array_add (self->sphaYa_array, sphaYa);

          for (i = 0; i < NCM_SPHERE_MAP_BLOCK_NC; i++)
          {
            const gint64 r_i   = r_ini + i;
            const gint64 apr_i = self->nrings - r_i - 1;
            const gint64 api   = NCM_SPHERE_MAP_BLOCK_NCT - i - 1;

            _ncm_sphere_map_prepare_circle (smap, block, r_i,   i);
            _ncm_sphere_map_prepare_circle (smap, block, apr_i, api);
          }

          ncm_sf_spherical_harmonics_start_rec_array (self->spha, sphaYa, NCM_SPHERE_MAP_BLOCK_NCT, block->theta);
        }
      }
    }
  }
}

/**
 * ncm_sphere_map_get_nside:
 * @smap: a #NcmSphereMap
 *
 * Get the nside parameter of @smap.
 *
 * Returns: the nside parameter of @smap.
 */
gint64
ncm_sphere_map_get_nside (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  return self->nside;
}

/**
 * ncm_sphere_map_get_npix:
 * @smap: a #NcmSphereMap
 *
 * Gets the number of pixels of @smap.
 *
 * Returns: the number of pixels of @smap.
 */
gint64
ncm_sphere_map_get_npix (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  return self->npix;
}

/**
 * ncm_sphere_map_get_cap_size:
 * @smap: a #NcmSphereMap
 *
 * Get the number of pixels of the north cap of @smap.
 *
 * Returns: the number of pixels of the north cap of @smap.
 */
gint64
ncm_sphere_map_get_cap_size (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  return self->cap_size;
}

/**
 * ncm_sphere_map_get_middle_size:
 * @smap: a #NcmSphereMap
 *
 * Get the number of pixels of the middle of @smap. That is
 * the number of pixels of the middle rings.
 *
 * Returns: the number of pixels of the middle of @smap.
 */
gint64
ncm_sphere_map_get_middle_size (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  return self->middle_size;
}

/**
 * ncm_sphere_map_get_nrings:
 * @smap: a #NcmSphereMap
 *
 * Get the number of rings of @smap.
 *
 * Returns: the number of rings of @smap.
 */
gint64
ncm_sphere_map_get_nrings (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  return self->nrings;
}

/**
 * ncm_sphere_map_get_nrings_cap:
 * @smap: a #NcmSphereMap
 *
 * Get the number of rings of the north cap of @smap.
 *
 * Returns: the number of rings of the north cap of @smap.
 */
gint64
ncm_sphere_map_get_nrings_cap (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  return self->nrings_cap;
}

/**
 * ncm_sphere_map_get_nrings_middle:
 * @smap: a #NcmSphereMap
 *
 * Get the number of rings of the middle of @smap.
 *
 * Returns: the number of rings of the middle of @smap.
 */
gint64
ncm_sphere_map_get_nrings_middle (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  return self->nrings_middle;
}

/**
 * ncm_sphere_map_get_ring_size:
 * @smap: a #NcmSphereMap
 * @r_i: ring index
 *
 * Get the number of pixels of the ring @r_i of @smap, counting rings from the north pole,
 * $0 \le$ @r_i $<$ ncm_sphere_map_get_nrings().
 *
 * Returns: the number of pixels of the ring @r_i of @smap.
 */
gint64
ncm_sphere_map_get_ring_size (NcmSphereMap *smap, gint64 r_i)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  _ncm_sphere_map_check_index (r_i, self->nrings, G_STRFUNC, "ring");

  if (r_i < self->nrings_cap) /* North cap (nside - 1 rings) */
  {
    return 4 * (r_i + 1);
  }
  else if (r_i >= self->nrings_cap + self->nrings_middle) /* South cap (it has 3nside rings above) */
  {
    const gint64 rsc_i = r_i - (self->nrings_cap + self->nrings_middle); /* South cap index */

    return 4 * (self->nside - 1 - rsc_i);
  }
  else /* Middle */
  {
    return self->middle_rings_size;
  }
}

/**
 * ncm_sphere_map_get_ring_first_index:
 * @smap: a #NcmSphereMap
 * @r_i: ring index
 *
 * Get the RING-order index of the first pixel of the ring @r_i of @smap, see
 * ncm_sphere_map_get_ring_size().
 *
 * Returns: the first pixel index of the ring @r_i of @smap.
 */
gint64
ncm_sphere_map_get_ring_first_index (NcmSphereMap *smap, gint64 r_i)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  _ncm_sphere_map_check_index (r_i, self->nrings, G_STRFUNC, "ring");

  if (r_i < self->nrings_cap) /* North cap (nside - 1 rings) */
  {
    return 2 * r_i * (r_i + 1);
  }
  else if (r_i >= self->nrings_cap + self->nrings_middle) /* South cap (it has 3nside rings above) */
  {
    const gint64 rsc_i = r_i - (self->nrings_cap + self->nrings_middle); /* South cap index */

    return self->cap_size + self->middle_size + 2 * rsc_i * (2 * self->nside - 1 - rsc_i);
  }
  else /* Middle */
  {
    const gint64 rm_i = r_i - self->nrings_cap; /* Middle index */

    return self->cap_size + self->middle_rings_size * rm_i;
  }
}

/**
 * ncm_sphere_map_set_order:
 * @smap: a #NcmSphereMap
 * @order: the pixel ordering
 *
 * Sets the pixel ordering of @smap.
 *
 */
void
ncm_sphere_map_set_order (NcmSphereMap *smap, NcmSphereMapOrder order)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  if (self->order != order)
  {
    if (self->nside > 0)
    {
      gpointer temp_pix = _fft_vec_alloc (self->npix);

      gint64 (*convert) (NcmSphereMap *smap, gint64) = (order == NCM_SPHERE_MAP_ORDER_NEST ? ncm_sphere_map_ring2nest : ncm_sphere_map_nest2ring);

      gint64 i, j;

      _fft_vec_memcpy (temp_pix, self->pvec, self->npix);

      /* The copy keeps the map's own precision: a gfloat here rounded every pixel of a
       * double map at 6e-8 on each change of ordering. */
      for (i = 0; i < self->npix; i++)
      {
        const _fft_real val = _fft_vec_idx (temp_pix, i);

        j                            = convert (smap, i);
        _fft_vec_idx (self->pvec, j) = val;
      }

      _fft_vec_free (temp_pix);
    }

    self->order = order;
  }
}

/**
 * ncm_sphere_map_get_order:
 * @smap: a #NcmSphereMap
 *
 * Gets the pixel ordering of @smap.
 *
 * Returns: the pixel ordering of @smap.
 */
NcmSphereMapOrder
ncm_sphere_map_get_order (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  return self->order;
}

/**
 * ncm_sphere_map_set_coordsys:
 * @smap: a #NcmSphereMap
 * @coordsys: the coordinate system
 *
 * Sets the coordinate system of @smap. It labels the map, for instance in the files of
 * ncm_sphere_map_save_fits(); the pixel values are not transformed.
 *
 */
void
ncm_sphere_map_set_coordsys (NcmSphereMap *smap, NcmSphereMapCoordSys coordsys)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  if (self->coordsys != coordsys)
  {
    switch (coordsys)
    {
      case NCM_SPHERE_MAP_COORD_SYS_GALACTIC:
      case NCM_SPHERE_MAP_COORD_SYS_ECLIPTIC:
      case NCM_SPHERE_MAP_COORD_SYS_CELESTIAL:
        self->coordsys = coordsys;
        break;
      default:
        g_error ("ncm_sphere_map_set_coordsys: Unknown coordinate system `%d'.", coordsys);
        break;
    }
  }
}

/**
 * ncm_sphere_map_get_coordsys:
 * @smap: a #NcmSphereMap
 *
 * Gets the coordinate system of @smap.
 *
 * Returns: the coordinate system of @smap.
 */
NcmSphereMapCoordSys
ncm_sphere_map_get_coordsys (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  return self->coordsys;
}

/**
 * ncm_sphere_map_set_lmax:
 * @smap: a #NcmSphereMap
 * @lmax: max value of $\ell$
 *
 * Prepare the object to calculate $a_{\ell{}m}$ and/or $C_\ell$,
 * up to $\ell=$@lmax.
 *
 */
void
ncm_sphere_map_set_lmax (NcmSphereMap *smap, guint lmax)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  if (self->lmax != lmax)
  {
    g_clear_pointer (&self->alm,  _fft_vec_free);
    ncm_vector_clear (&self->Cl);
    self->lmax = lmax;

    if (self->lmax > 0)
    {
      self->Cl      = ncm_vector_new (self->lmax + 1);
      self->alm_len = NCM_SPHERE_MAP_ALM_SIZE (self->lmax);
      self->alm     = _fft_vec_alloc_complex (self->alm_len);
      memset (self->alm, 0, sizeof (_fft_complex) * self->alm_len);

      if (self->lmax + 2 > ncm_sf_spherical_harmonics_get_lmax (self->spha))
        ncm_sf_spherical_harmonics_set_lmax (self->spha, self->lmax + 2);
    }
  }
}

/**
 * ncm_sphere_map_get_lmax:
 * @smap: a #NcmSphereMap
 *
 * Returns: the largest $\ell$ of the harmonic transforms of @smap.
 */
guint
ncm_sphere_map_get_lmax (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  return self->lmax;
}

/**
 * ncm_sphere_map_set_iter:
 * @smap: a #NcmSphereMap
 * @iter: number of iterations for map2alm transform
 *
 * Sets the number of iterations to use in the spherical harmonic
 * transform. Iterations improve the accuracy of the transform for
 * undersampled maps (when lmax < 3 * nside). iter=0 uses only the
 * direct transform; healpy's map2alm defaults to 3 iterations.
 *
 */
void
ncm_sphere_map_set_iter (NcmSphereMap *smap, guint iter)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  self->iter = iter;
}

/**
 * ncm_sphere_map_get_iter:
 * @smap: a #NcmSphereMap
 *
 * Gets the number of iterations used in the spherical harmonic transform.
 *
 * Returns: the number of iterations
 */
guint
ncm_sphere_map_get_iter (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  return self->iter;
}

/**
 * ncm_sphere_map_clear_pixels:
 * @smap: a #NcmSphereMap
 *
 * Clear the pixels of @smap, i.e. set all pixels to zero.
 *
 */
void
ncm_sphere_map_clear_pixels (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  if (self->nside > 0)
    _fft_vec_set_zero (self->pvec, self->npix);
}

/**
 * ncm_sphere_map_nest2ring:
 * @smap: a #NcmSphereMap
 * @nest_index: a pixel index in NESTED ordering
 *
 * Converts a pixel index from NESTED to RING ordering.
 *
 * Returns: the pixel index in RING ordering.
 */
gint64
ncm_sphere_map_nest2ring (NcmSphereMap *smap, const gint64 nest_index)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  gint64 ring_index;
  gint64 base_pixel;
  gint64 f, h, v;            /* Face number, horizontal coordinate, vertical coordinate */
  gint64 t, p, s, pad, w, l; /* theta, phi, shift, padding, width, local index          */
  gint64 x, y;
  gint64 hf;

  _ncm_sphere_map_check_index (nest_index, self->npix, G_STRFUNC, "pixel");

  f = nest_index / self->face_size;
  l = nest_index % self->face_size;

  NCM_SPHERE_MAP_INT_TO_XY (l, x, y);

  h = self->nside - 1 - y;
  v = 2 * self->nside - 2 - y - x;

  switch (f / 4)
  {
    case 0:
      t = v;

      if (t < (self->nside - 1))
      {
        w          = (t + 1);
        s          = 0;
        pad        = 0;
        hf         = f % 4;
        base_pixel = 2 * (_l_pow_2 (t + 1) - t - 1);
      }
      else
      {
        w          = self->nside;
        s          = (t - (self->nside - 1)) / 2;
        pad        = 0;
        hf         = f % 4;
        base_pixel = self->cap_size + (t - self->nside + 1) * self->middle_rings_size;
      }

      break;
    case 1:
      t          = v + self->nside;
      w          = self->nside;
      s          = (t - (self->nside - 1)) / 2;
      pad        = 0;
      hf         = f % 4;
      base_pixel = self->cap_size + (v + 1) * self->middle_rings_size;
      break;
    case 2:
      t = v + 2 * self->nside;

      if (v < self->nside)
      {
        w          = self->nside;
        s          = (t - (self->nside - 1)) / 2;
        pad        = 0;
        hf         = f % 4 + 1;
        base_pixel = self->cap_size + (v + 1 + self->nside) * self->middle_rings_size;
      }
      else
      {
        w          = 4 * self->nside - t - 1;
        s          = w;
        pad        = (t - 3 * self->nside + 1);
        hf         = f % 4 + 1;
        base_pixel = (self->npix - self->cap_size) + (4 * self->nside * (v - self->nside + 1) - 2 * _l_pow_2 (v - self->nside + 1) + 2 * (v - self->nside) + 2 - 4 * self->nside);
      }

      break;
    default:
      g_assert_not_reached ();
      break;
  }

  p = (h + hf * w - s - pad);

  if (p < 0)
    p += 4 * w;

  ring_index = base_pixel + p;

  return ring_index;
}

/**
 * ncm_sphere_map_ring2nest:
 * @smap: a #NcmSphereMap
 * @ring_index: a pixel index in RING ordering
 *
 * Converts a pixel index from RING to NESTED ordering.
 *
 * Returns: the pixel index in NESTED ordering.
 */
gint64
ncm_sphere_map_ring2nest (NcmSphereMap *smap, const gint64 ring_index)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  gint64 nest_index;
  gint64 f, h, v;            /* Face number, horizontal coordinate, vertical coordinate */
  gint64 t, p, s, pad, w, l; /* theta, phi, shift, padding, width, local index          */
  gint64 x, y;

  _ncm_sphere_map_check_index (ring_index, self->npix, G_STRFUNC, "pixel");

  if (ring_index < self->cap_size)
  {
    t   = (sqrt (1 + 2 * ring_index) - 1) / 2;
    w   = (t + 1);
    p   = ring_index - 2 * (_l_pow_2 (t + 1) - t - 1);
    s   = 0;
    pad = 0;
  }
  else if (ring_index < (self->npix - self->cap_size))
  {
    l   = ring_index - self->cap_size;
    w   = self->nside;
    t   = (gint64) (l / self->middle_rings_size) + (self->nside - 1);
    p   = l % self->middle_rings_size;
    s   = (t - (self->nside - 1)) / 2;
    pad = 0;
  }
  else
  {
    l   = ring_index - self->npix + self->cap_size;
    t   = 4 * self->nside - (1 + sqrt (4 * _l_pow_2 (self->nside) - 4 * self->nside + 1 - 2 * l)) / 2;
    w   = 4 * self->nside - t - 1;
    p   = l - (4 * self->nside * (t - 3 * self->nside + 1) - 2 * _l_pow_2 (t - 3 * self->nside + 1) + 2 * (t - 3 * self->nside) + 2 - 4 * self->nside);
    s   = w;
    pad = (t - 3 * self->nside + 1);
  }

  h  = (p + s) % w + pad;
  v  = (t - h) % self->nside + h;
  f  = ((gint64) ((t - h) / self->nside)) * 4;
  f += (f == 8) ? ((gint64) ((p + s) / w) - 1) : ((gint64) ((p + s) / w) % 4);

  x = self->nside - 1 - (v - h);
  y = self->nside - 1 - h;

  NCM_SPHERE_MAP_XY_TO_INT (x, y, nest_index);

  nest_index += self->face_size * f;

  return nest_index;
}

static void _ncm_sphere_map_nest_to_tpw (NcmSphereMapPrivate * const self, const gint64 nest_index, gint64 *t, gint64 *p, gint64 *w);
static void _ncm_sphere_map_ring_to_tpw (NcmSphereMapPrivate * const self, const gint64 ring_index, gint64 *t, gint64 *p, gint64 *w);
static void _ncm_sphere_map_tpw_to_centre (const gint64 nside, const gint64 t0, const gint64 p0, const gint64 w, gdouble *z, gdouble *sin_theta, gdouble *phi);

/**
 * ncm_sphere_map_pix2ang_nest:
 * @smap: a #NcmSphereMap
 * @nest_index: a pixel index in NESTED ordering
 * @theta: (out): the polar angle of the pixel centre, in $[0, \pi]$ (radians)
 * @phi: (out): the azimuth of the pixel centre, in $[0, 2\pi)$ (radians)
 *
 * Gets the centre of the pixel @nest_index. The polar angle is measured from the north
 * pole.
 *
 */
void
ncm_sphere_map_pix2ang_nest (NcmSphereMap *smap, const gint64 nest_index, gdouble *theta, gdouble *phi)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  gint64 t, p, w;
  gdouble z, sin_theta;

  _ncm_sphere_map_check_index (nest_index, self->npix, G_STRFUNC, "pixel");
  _ncm_sphere_map_nest_to_tpw (self, nest_index, &t, &p, &w);
  _ncm_sphere_map_tpw_to_centre (self->nside, t, p, w, &z, &sin_theta, phi);

  *theta = atan2 (sin_theta, z);
}

/**
 * ncm_sphere_map_pix2ang_ring:
 * @smap: a #NcmSphereMap
 * @ring_index: a pixel index in RING ordering
 * @theta: (out): the polar angle of the pixel centre, in $[0, \pi]$ (radians)
 * @phi: (out): the azimuth of the pixel centre, in $[0, 2\pi)$ (radians)
 *
 * Gets the centre of the pixel @ring_index, see ncm_sphere_map_pix2ang_nest().
 *
 */
void
ncm_sphere_map_pix2ang_ring (NcmSphereMap *smap, const gint64 ring_index, gdouble *theta, gdouble *phi)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  gint64 t, p, w;
  gdouble z, sin_theta;

  _ncm_sphere_map_check_index (ring_index, self->npix, G_STRFUNC, "pixel");
  _ncm_sphere_map_ring_to_tpw (self, ring_index, &t, &p, &w);
  _ncm_sphere_map_tpw_to_centre (self->nside, t, p, w, &z, &sin_theta, phi);

  *theta = atan2 (sin_theta, z);
}

/**
 * ncm_sphere_map_pix2vec_ring:
 * @smap: a #NcmSphereMap
 * @ring_index: a pixel index in RING ordering
 * @vec: a #NcmTriVec
 *
 * Sets @vec to the unit vector of the centre of the pixel @ring_index.
 *
 */
void
ncm_sphere_map_pix2vec_ring (NcmSphereMap *smap, gint64 ring_index, NcmTriVec *vec)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  gint64 t, p, w;
  gdouble z, sin_theta, phi;

  _ncm_sphere_map_check_index (ring_index, self->npix, G_STRFUNC, "pixel");
  _ncm_sphere_map_ring_to_tpw (self, ring_index, &t, &p, &w);
  _ncm_sphere_map_tpw_to_centre (self->nside, t, p, w, &z, &sin_theta, &phi);

  vec->c[0] = sin_theta * cos (phi);
  vec->c[1] = sin_theta * sin (phi);
  vec->c[2] = z;
}

/**
 * ncm_sphere_map_pix2vec_nest:
 * @smap: a #NcmSphereMap
 * @nest_index: a pixel index in NESTED ordering
 * @vec: a #NcmTriVec
 *
 * Sets @vec to the unit vector of the centre of the pixel @nest_index.
 *
 */
void
ncm_sphere_map_pix2vec_nest (NcmSphereMap *smap, gint64 nest_index, NcmTriVec *vec)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  gint64 t, p, w;
  gdouble z, sin_theta, phi;

  _ncm_sphere_map_check_index (nest_index, self->npix, G_STRFUNC, "pixel");
  _ncm_sphere_map_nest_to_tpw (self, nest_index, &t, &p, &w);
  _ncm_sphere_map_tpw_to_centre (self->nside, t, p, w, &z, &sin_theta, &phi);

  vec->c[0] = sin_theta * cos (phi);
  vec->c[1] = sin_theta * sin (phi);
  vec->c[2] = z;
}

/* Ring t (from 0 at the north pole), position p in the ring (from 0) and the ring's
 * pixels per face w of a NESTED pixel. */
static void
_ncm_sphere_map_nest_to_tpw (NcmSphereMapPrivate * const self, const gint64 nest_index, gint64 *t_out, gint64 *p_out, gint64 *w_out)
{
  const gint64 f = nest_index / self->face_size; /* Face number */
  const gint64 l = nest_index % self->face_size; /* Index within the face */
  gint64 h, v, t, p, s, pad, w, x, y, hf;        /* Horizontal and vertical coordinates, shift, padding */

  NCM_SPHERE_MAP_INT_TO_XY (l, x, y);

  h = self->nside - 1 - y;
  v = 2 * self->nside - 2 - y - x;

  switch (f / 4)
  {
    case 0:
      t   = v;
      w   = (t < (self->nside - 1)) ? (t + 1) : self->nside;
      s   = (t < (self->nside - 1)) ? 0 : (t - (self->nside - 1)) / 2;
      pad = 0;
      hf  = f % 4;
      break;
    case 1:
      t   = v + self->nside;
      w   = self->nside;
      s   = (t - (self->nside - 1)) / 2;
      pad = 0;
      hf  = f % 4;
      break;
    case 2:
      t = v + 2 * self->nside;

      if (v < self->nside)
      {
        w   = self->nside;
        s   = (t - (self->nside - 1)) / 2;
        pad = 0;
      }
      else
      {
        w   = 4 * self->nside - t - 1;
        s   = w;
        pad = (t - 3 * self->nside + 1);
      }

      hf = f % 4 + 1;
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  p = (h + hf * w - s - pad);

  if (p < 0)
    p += 4 * w;

  *t_out = t;
  *p_out = p;
  *w_out = w;
}

/* Ring t, position p and pixels per face w of a RING pixel, see
 * _ncm_sphere_map_nest_to_tpw(). */
static void
_ncm_sphere_map_ring_to_tpw (NcmSphereMapPrivate * const self, const gint64 ring_index, gint64 *t, gint64 *p, gint64 *w)
{
  if (ring_index < self->cap_size)
  {
    *t = (sqrt (1 + 2 * ring_index) - 1) / 2;
    *w = (*t + 1);
    *p = ring_index - 2 * (_l_pow_2 (*t + 1) - *t - 1);
  }
  else if (ring_index < (self->npix - self->cap_size))
  {
    const gint64 l = ring_index - self->cap_size;

    *w = self->nside;
    *t = (gint64) (l / self->middle_rings_size) + (self->nside - 1);
    *p = l % self->middle_rings_size;
  }
  else
  {
    const gint64 l = ring_index - self->npix + self->cap_size;

    *t = 4 * self->nside - (1 + sqrt (4 * _l_pow_2 (self->nside) - 4 * self->nside + 1 - 2 * l)) / 2;
    *w = 4 * self->nside - *t - 1;
    *p = l - (4 * self->nside * (*t - 3 * self->nside + 1) - 2 * _l_pow_2 (*t - 3 * self->nside + 1) + 2 * (*t - 3 * self->nside) + 2 - 4 * self->nside);
  }
}

/* Centre of the pixel (t, p, w): z = cos(theta), sin(theta) and phi. In the polar caps
 * 1 - |z| = t'^2 / (3 nside^2) is exact from the ring index t' (from 1 at the pole), so
 * sin(theta) = sqrt((1 - |z|)(1 + |z|)) keeps full precision; 1 - z^2 from z loses it as
 * nside^2 (5.5e-11 at nside 1024), and so does acos (z). */
static void
_ncm_sphere_map_tpw_to_centre (const gint64 nside, const gint64 t0, const gint64 p0, const gint64 w, gdouble *z, gdouble *sin_theta, gdouble *phi)
{
  const gint64 t = t0 + 1;
  const gint64 p = p0 + 1;

  if ((t < nside) || (t > 3 * nside))
  {
    const gint64 tc     = (t < nside) ? t : 4 * nside - t;
    const gdouble onemz = tc * tc / (3.0 * nside * nside);
    const gdouble abs_z = 1.0 - onemz;

    *z         = (t < nside) ? abs_z : -abs_z;
    *sin_theta = sqrt (onemz * (1.0 + abs_z));
    *phi       = (p - 0.5) * M_PI_2 / (1.0 * w);
  }
  else
  {
    *z         = (2.0 * nside - t) * 2.0 / (3.0 * nside);
    *sin_theta = sqrt ((1.0 - *z) * (1.0 + *z));
    *phi       = (p - ((t - nside) % 2 + 1.0) * 0.5) * M_PI_2 / (1.0 * w);
  }
}

static void
_ncm_sphere_map_zphi2pix_nest (NcmSphereMap *smap, const gdouble z, const gdouble onemz2, const gdouble phi, gint64 *nest_index)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  const gdouble two_3              = 2.0 / 3.0;
  const gdouble abs_z              = fabs (z);
  const gdouble phi_norm           = (phi < 0.0) ? fmod (phi, 2.0 * M_PI) + 2.0 * M_PI : phi;
  const gdouble tt                 = fmod (phi_norm / M_PI_2, 4.0);
  gint64 x, y, f;

  if (abs_z <= two_3)
  {
    const gdouble temp1 = self->nside * (0.5 + tt);
    const gdouble temp2 = self->nside * (z * 0.75);

    const gint64 jp  = (gint64) (temp1 - temp2);
    const gint64 jm  = (gint64) (temp1 + temp2);
    const gint64 ifp = jp / self->nside;
    const gint64 ifm = jm / self->nside;

    f = (ifp == ifm) ? (ifp | 4) : ((ifp < ifm) ? ifp : (ifm + 8));
    x = jm & (self->nside - 1);
    y = self->nside - (jp & (self->nside - 1)) - 1;
  }
  else
  {
    const gint64 ntt  = GSL_MIN (3, (gint64) (tt));
    const gdouble tp  = tt - ntt;
    const gdouble tmp = self->nside * sqrt (3.0 * onemz2 / (1.0 + abs_z));

    const gint64 jp  = (gint64) (tp * tmp);
    const gint64 jm  = (gint64) ((1.0 - tp) * tmp);
    const gint64 jpm = GSL_MIN (jp, self->nside - 1);
    const gint64 jmm = GSL_MIN (jm, self->nside - 1);

    if (z >= 0.0)
    {
      x = self->nside - jmm - 1;
      y = self->nside - jpm - 1;
      f = ntt;
    }
    else
    {
      x = jpm;
      y = jmm;
      f = ntt + 8;
    }
  }

  NCM_SPHERE_MAP_XY_TO_INT (x, y, nest_index[0]);

  nest_index[0] += self->face_size * f;
}

static void
_ncm_sphere_map_zphi2pix_ring (NcmSphereMap *smap, const gdouble z, const gdouble onemz2, const gdouble phi, gint64 *ring_index)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  const gdouble two_3              = 2.0 / 3.0;
  const gdouble abs_z              = fabs (z);
  const gdouble phi_norm           = (phi < 0.0) ? fmod (phi, 2.0 * M_PI) + 2.0 * M_PI : phi;
  const gdouble tt                 = fmod (phi_norm / M_PI_2, 4.0);

  if (abs_z <= two_3)
  {
    const gdouble temp1 = self->nside * (0.5 + tt);
    const gdouble temp2 = self->nside * (z * 0.75);

    const gint64 jp = (gint64) (temp1 - temp2);
    const gint64 jm = (gint64) (temp1 + temp2);

    const gint64 ir     = self->nside + 1 + jp - jm;
    const gint64 kshift = 1 - (ir & 1);

    const gint64 t1 = (jp + jm - self->nside + kshift + 1) / 2;
    const gint64 ip = t1 % self->middle_rings_size;

    ring_index[0] = self->cap_size + (ir - 1) * self->middle_rings_size + ip;
  }
  else /* North & South polar caps */
  {
    const gdouble tp  = tt - (gint64) (tt);
    const gdouble tmp = self->nside * sqrt (3.0 * onemz2 / (1.0 + abs_z));

    const gint64 jp = (gint64) (tp * tmp);
    const gint64 jm = (gint64) ((1.0 - tp) * tmp);

    const gint64 ir = jp + jm + 1;
    const gint64 ip = (gint64) (tt * ir);

    if (z > 0)
      ring_index[0] = 2 * ir * (ir - 1) + ip;
    else
      ring_index[0] = self->npix - 2 * ir * (ir + 1) + ip;
  }
}

/**
 * ncm_sphere_map_ang2pix_nest:
 * @smap: a #NcmSphereMap
 * @theta: a polar angle, from the north pole (radians)
 * @phi: an azimuth (radians), any value, taken modulo $2\pi$
 * @nest_index: (out): the pixel index in NESTED ordering
 *
 * Gets the NESTED index of the pixel that contains the direction (@theta, @phi).
 *
 */
void
ncm_sphere_map_ang2pix_nest (NcmSphereMap *smap, const gdouble theta, const gdouble phi, gint64 *nest_index)
{
  const gdouble z = cos (theta);

  if (theta > 0.1)
    _ncm_sphere_map_zphi2pix_nest (smap, z, 1.0 - gsl_pow_2 (z), phi, nest_index);
  else
    _ncm_sphere_map_zphi2pix_nest (smap, z, gsl_pow_2 (sin (theta)), phi, nest_index);
}

/**
 * ncm_sphere_map_ang2pix_ring:
 * @smap: a #NcmSphereMap
 * @theta: a polar angle, from the north pole (radians)
 * @phi: an azimuth (radians), any value, taken modulo $2\pi$
 * @ring_index: (out): the pixel index in RING ordering
 *
 * Gets the RING index of the pixel that contains the direction (@theta, @phi).
 *
 */
void
ncm_sphere_map_ang2pix_ring (NcmSphereMap *smap, const gdouble theta, const gdouble phi, gint64 *ring_index)
{
  const gdouble z = cos (theta);

  if (theta > 0.1)
    _ncm_sphere_map_zphi2pix_ring (smap, z, 1.0 - gsl_pow_2 (z), phi, ring_index);
  else
    _ncm_sphere_map_zphi2pix_ring (smap, z, gsl_pow_2 (sin (theta)), phi, ring_index);
}

/**
 * ncm_sphere_map_vec2pix_ring:
 * @smap: a #NcmSphereMap
 * @vec: a #NcmTriVec
 * @ring_index: (out): the pixel index in RING ordering
 *
 * Gets the RING index of the pixel that contains the direction of @vec, any non-zero
 * vector.
 *
 */
void
ncm_sphere_map_vec2pix_ring (NcmSphereMap *smap, NcmTriVec *vec, gint64 *ring_index)
{
  const gdouble norm    = ncm_trivec_norm (vec);
  const gdouble z       = vec->c[2] / norm;
  const gdouble sin_the = hypot (vec->c[0], vec->c[1]) / norm; /* 1 - z^2 would cancel near the poles */
  const gdouble phi     = ncm_trivec_get_phi (vec);

  _ncm_sphere_map_zphi2pix_ring (smap, z, gsl_pow_2 (sin_the), phi, ring_index);
}

/**
 * ncm_sphere_map_vec2pix_nest:
 * @smap: a #NcmSphereMap
 * @vec: a #NcmTriVec
 * @nest_index: (out): the pixel index in NESTED ordering
 *
 * Gets the NESTED index of the pixel that contains the direction of @vec, any non-zero
 * vector.
 *
 */
void
ncm_sphere_map_vec2pix_nest (NcmSphereMap *smap, NcmTriVec *vec, gint64 *nest_index)
{
  const gdouble norm    = ncm_trivec_norm (vec);
  const gdouble z       = vec->c[2] / norm;
  const gdouble sin_the = hypot (vec->c[0], vec->c[1]) / norm; /* 1 - z^2 would cancel near the poles */
  const gdouble phi     = ncm_trivec_get_phi (vec);

  _ncm_sphere_map_zphi2pix_nest (smap, z, gsl_pow_2 (sin_the), phi, nest_index);
}

/**
 * ncm_sphere_map_add_to_vec:
 * @smap: a #NcmSphereMap
 * @vec: a #NcmTriVec
 * @s: signal
 *
 * Adds @s to the pixel of @smap, in its current ordering, that contains the direction
 * of @vec.
 *
 */
void
ncm_sphere_map_add_to_vec (NcmSphereMap *smap, NcmTriVec *vec, const gdouble s)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  gint64 index                     = 0;

  switch (self->order)
  {
    case NCM_SPHERE_MAP_ORDER_NEST:
      ncm_sphere_map_vec2pix_nest (smap, vec, &index);
      break;
    case NCM_SPHERE_MAP_ORDER_RING:
      ncm_sphere_map_vec2pix_ring (smap, vec, &index);
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  _fft_vec_idx (self->pvec, index) += s;
}

/**
 * ncm_sphere_map_add_to_ang:
 * @smap: a #NcmSphereMap
 * @theta: $\theta$
 * @phi: $\phi$
 * @s: signal
 *
 * Adds @s to the pixel of @smap, in its current ordering, that contains the direction
 * (@theta, @phi), see ncm_sphere_map_ang2pix_ring().
 *
 */
void
ncm_sphere_map_add_to_ang (NcmSphereMap *smap, const gdouble theta, const gdouble phi, const gdouble s)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  gint64 index                     = 0;

  switch (self->order)
  {
    case NCM_SPHERE_MAP_ORDER_NEST:
      ncm_sphere_map_ang2pix_nest (smap, theta, phi, &index);
      break;
    case NCM_SPHERE_MAP_ORDER_RING:
      ncm_sphere_map_ang2pix_ring (smap, theta, phi, &index);
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  _fft_vec_idx (self->pvec, index) += s;
}

/**
 * ncm_sphere_map_load_fits:
 * @smap: a #NcmSphereMap
 * @fits_file: fits filename
 * @signal_name: (allow-none): signal column name in @fits_file, %NULL for the first column
 *
 * Loads a full-sky map in the HEALPix FITS format: a binary table in the first extension
 * with the map in one column, any number of pixels per row, and the keys NSIDE, ORDERING
 * (RING if missing, with a warning) and COORDSYS (celestial if missing, with a warning).
 * The nside, ordering and coordinate system of @smap are set from the file. A PIXTYPE other
 * than HEALPIX, a partial-sky map (INDXSCHM = EXPLICIT) or a column whose size does not
 * match NSIDE abort.
 *
 */
void
ncm_sphere_map_load_fits (NcmSphereMap *smap, const gchar *fits_file, const gchar *signal_name)
{
#ifdef HAVE_CFITSIO
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  gchar comment[FLEN_COMMENT];
  gchar ordering[FLEN_VALUE];
  gchar coordsys[FLEN_VALUE];
  gchar value[FLEN_VALUE];
  gint status, hdutype, anynul, typecode;
  glong nside, naxis2, repeat, width;
  gint signal_i = 1;
  fitsfile *fptr;

  G_STATIC_ASSERT (sizeof (_fft_real) == sizeof (gdouble));

  status = 0;

  fits_open_file (&fptr, fits_file, READONLY, &status);
  NCM_FITS_ERROR (status);

  fits_movabs_hdu (fptr, 2, &hdutype, &status);
  NCM_FITS_ERROR (status);

  if (hdutype != BINARY_TBL)
    g_error ("ncm_sphere_map_load_fits: `%s' is not a binary table.", fits_file);

  if (fits_read_key (fptr, TSTRING, "PIXTYPE", value, comment, &status) == 0)
  {
    if (strncmp (value, "HEALPIX", 7) != 0)
      g_error ("ncm_sphere_map_load_fits: `%s' has PIXTYPE `%s', not HEALPIX.", fits_file, value);
  }
  else
  {
    status = 0;
  }

  if (fits_read_key (fptr, TSTRING, "INDXSCHM", value, comment, &status) == 0)
  {
    if (strncmp (value, "EXPLICIT", 8) == 0)
      g_error ("ncm_sphere_map_load_fits: `%s' is a partial-sky map (INDXSCHM = EXPLICIT), which is not supported.", fits_file);
  }
  else
  {
    status = 0;
  }

  fits_read_key_lng (fptr, "NSIDE", &nside, comment, &status);
  NCM_FITS_ERROR (status);

  if (nside <= 0)
    g_error ("ncm_sphere_map_load_fits: `%s' has NSIDE = %ld.", fits_file, nside);

  ncm_sphere_map_set_nside (smap, nside);

  fits_read_key_lng (fptr, "NAXIS2", &naxis2, comment, &status);
  NCM_FITS_ERROR (status);

  if ((signal_name != NULL) && fits_get_colnum (fptr, CASESEN, (gchar *) signal_name, &signal_i, &status))
    g_error ("ncm_sphere_map_load_fits: signal column named `%s' not found in `%s'.",
             signal_name, fits_file);

  fits_get_coltype (fptr, signal_i, &typecode, &repeat, &width, &status);
  NCM_FITS_ERROR (status);

  /* HEALPix writers store several pixels per row (healpy: 1024); cfitsio reads the
   * elements in order across the rows. */
  if (naxis2 * repeat != self->npix)
    g_error ("ncm_sphere_map_load_fits: column %d of `%s' holds %ld values (%ld rows of %ld), "
             "a map of nside %ld has %" G_GINT64_FORMAT " pixels.",
             signal_i, fits_file, naxis2 * repeat, naxis2, repeat, nside, self->npix);

  fits_read_col_dbl (fptr, signal_i, 1, 1, self->npix, NCM_SPHERE_MAP_HEALPIX_NULLVAL,
                     _fft_vec_ptr (self->pvec, 0), &anynul, &status);
  NCM_FITS_ERROR (status);

  if (fits_read_key (fptr, TSTRING, "ORDERING", ordering, comment, &status))
  {
    g_warning ("ncm_sphere_map_load_fits: Could not find ORDERING in the fits file, assuming RING.");
    g_strlcpy (ordering, "RING", FLEN_VALUE);
    status = 0;
  }

  switch (*ordering)
  {
    case 'N':
      self->order = NCM_SPHERE_MAP_ORDER_NEST;
      break;
    case 'R':
      self->order = NCM_SPHERE_MAP_ORDER_RING;
      break;
    default:
      g_error ("ncm_sphere_map_load_fits: Unknown order type `%s'.", ordering);
      break;
  }

  if (fits_read_key (fptr, TSTRING, "COORDSYS", coordsys, comment, &status))
  {
    g_warning ("ncm_sphere_map_load_fits: Could not find COORDSYS in the fits file, assuming CELESTIAL.");
    coordsys[0] = NCM_SPHERE_MAP_COORD_SYS_CELESTIAL;
    coordsys[1] = '\0';
    status      = 0;
  }

  ncm_sphere_map_set_coordsys (smap, *coordsys);

  fits_close_file (fptr, &status);
  NCM_FITS_ERROR (status);
#else
  g_error ("ncm_sphere_map_load_fits: CFITSIO not available.");

#endif
}

/**
 * ncm_sphere_map_save_fits:
 * @smap: a #NcmSphereMap
 * @fits_file: fits filename
 * @signal_name: (allow-none): signal column name in @fits_file, %NULL for
 *   %NCM_SPHERE_MAP_DEFAULT_SIGNAL
 * @overwrite: whether to overwrite @fits_file if it exists
 *
 * Saves @smap in the HEALPix FITS format read by ncm_sphere_map_load_fits() and healpy:
 * a binary table with one double-precision value per row, in the current ordering, and
 * the keys PIXTYPE, ORDERING, NSIDE, FIRSTPIX, LASTPIX, INDXSCHM and COORDSYS. The values
 * are written in double precision, so a save and load returns the map unchanged.
 *
 */
void
ncm_sphere_map_save_fits (NcmSphereMap *smap, const gchar *fits_file, const gchar *signal_name, gboolean overwrite)
{
#ifdef HAVE_CFITSIO
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  const gchar *sname               = signal_name != NULL ?  signal_name : NCM_SPHERE_MAP_DEFAULT_SIGNAL;
  const gchar *ttype[]             = { sname };
  const gchar *tform[]             = { "1D" };
  const gint64 npix                = ncm_sphere_map_get_npix (smap);
  const gchar extname[]            = "BINTABLE";
  fitsfile *fptr;
  gint status = 0;

  if (overwrite && g_file_test (fits_file, G_FILE_TEST_EXISTS))
    g_unlink (fits_file);

  fits_create_file (&fptr, fits_file, &status);
  NCM_FITS_ERROR (status);

  fits_create_tbl (fptr, BINARY_TBL, npix, 1, (gchar **) ttype, (gchar **) tform, NULL, extname, &status);
  NCM_FITS_ERROR (status);

  fits_write_key (fptr, TSTRING, "PIXTYPE", (gchar *) "HEALPIX",
                  "HEALPIX pixelisation", &status);
  NCM_FITS_ERROR (status);

  {
    const gchar *order_str = ncm_sphere_map_get_order (smap) == NCM_SPHERE_MAP_ORDER_NEST ? "NESTED  " : "RING    ";

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
    glong nside = ncm_sphere_map_get_nside (smap);

    fits_write_key (fptr, TLONG, "NSIDE", &nside,
                    "Resolution parameter for HEALPIX", &status);
    NCM_FITS_ERROR (status);
  }

  {
    gchar coordsys_c = (gchar) ncm_sphere_map_get_coordsys (smap);
    gchar *coordsys  = g_strdup_printf ("%c       ", coordsys_c);

    fits_write_key (fptr, TSTRING, "COORDSYS", coordsys,
                    "Pixelization coordinate system", &status);
    NCM_FITS_ERROR (status);

    fits_write_comment (fptr,
                        "G = Galactic, E = ecliptic, C = celestial = equatorial", &status);
    NCM_FITS_ERROR (status);

    g_free (coordsys);
  }

  G_STATIC_ASSERT (sizeof (_fft_real) == sizeof (gdouble));

  fits_write_col (fptr, TDOUBLE, 1, 1, 1, npix, _fft_vec_ptr (self->pvec, 0), &status);
  NCM_FITS_ERROR (status);


  fits_close_file (fptr, &status);
  NCM_FITS_ERROR (status);
#else
  g_error ("ncm_sphere_map_save_fits: CFITSIO not available.");
#endif
}

static void
_ncm_sphere_map_radec_to_ang (const gdouble RA, const gdouble DEC, gdouble *theta, gdouble *phi)
{
  theta[0] = gsl_sf_angle_restrict_pos (ncm_c_degree_to_radian (90.0 - DEC));
  phi[0]   = gsl_sf_angle_restrict_pos (ncm_c_degree_to_radian (RA));
}

/**
 * ncm_sphere_map_load_from_fits_catalog:
 * @smap: a #NcmSphereMap
 * @fits_file: fits filename
 * @RA: right ascension column name in @fits_file (degrees)
 * @DEC: declination column name in @fits_file (degrees)
 * @S: (allow-none): signal column name in @fits_file, %NULL to count objects
 *
 * Adds the objects of a catalog, a binary table in the first extension of @fits_file, to
 * the pixels of @smap that contain them: each object adds its signal, or one when @S is
 * %NULL. The map is not cleared first, see ncm_sphere_map_clear_pixels().
 *
 */
void
ncm_sphere_map_load_from_fits_catalog (NcmSphereMap *smap, const gchar *fits_file, const gchar *RA, const gchar *DEC, const gchar *S)
{
#ifdef HAVE_CFITSIO
  gchar comment[FLEN_COMMENT];
  gint status, hdutype, anynul;
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
    g_error ("ncm_sphere_map_load_from_fits_catalog: `%s' is not a binary table.", fits_file);

  if (fits_get_colnum (fptr, CASESEN, (gchar *) RA, &RA_col, &status))
    g_error ("ncm_sphere_map_load_from_fits_catalog: RA column named `%s' not found in `%s'.",
             RA, fits_file);

  if (fits_get_colnum (fptr, CASESEN, (gchar *) DEC, &DEC_col, &status))
    g_error ("ncm_sphere_map_load_from_fits_catalog: DEC column named `%s' not found in `%s'.",
             DEC, fits_file);

  if (S != NULL)
    if (fits_get_colnum (fptr, CASESEN, (gchar *) S, &S_col, &status))
      g_error ("ncm_sphere_map_load_from_fits_catalog: Signal column named `%s' not found in `%s'.",
               S, fits_file);


  fits_read_key_lng (fptr, "NAXIS2", &naxis2, comment, &status);
  NCM_FITS_ERROR (status);

  if (S != NULL)
  {
    for (i = 0; i < naxis2; i++)
    {
      gdouble RA_i, DEC_i, theta, phi;
      gdouble S_i;

      fits_read_col_dbl (fptr, RA_col, 1 + i, 1, 1, NCM_SPHERE_MAP_HEALPIX_NULLVAL,
                         &RA_i, &anynul, &status);
      NCM_FITS_ERROR (status);

      fits_read_col_dbl (fptr, DEC_col, 1 + i, 1, 1, NCM_SPHERE_MAP_HEALPIX_NULLVAL,
                         &DEC_i, &anynul, &status);
      NCM_FITS_ERROR (status);

      fits_read_col_dbl (fptr, S_col, 1 + i, 1, 1, NCM_SPHERE_MAP_HEALPIX_NULLVAL,
                         &S_i, &anynul, &status);
      NCM_FITS_ERROR (status);

      _ncm_sphere_map_radec_to_ang (RA_i, DEC_i, &theta, &phi);

      ncm_sphere_map_add_to_ang (smap, theta, phi, S_i);
    }
  }
  else
  {
    for (i = 0; i < naxis2; i++)
    {
      gdouble RA_i, DEC_i, theta, phi;

      fits_read_col_dbl (fptr, RA_col, 1 + i, 1, 1, NCM_SPHERE_MAP_HEALPIX_NULLVAL,
                         &RA_i, &anynul, &status);
      NCM_FITS_ERROR (status);

      fits_read_col_dbl (fptr, DEC_col, 1 + i, 1, 1, NCM_SPHERE_MAP_HEALPIX_NULLVAL,
                         &DEC_i, &anynul, &status);
      NCM_FITS_ERROR (status);

      _ncm_sphere_map_radec_to_ang (RA_i, DEC_i, &theta, &phi);

      ncm_sphere_map_add_to_ang (smap, theta, phi, 1.0);
    }
  }

  fits_close_file (fptr, &status);
  NCM_FITS_ERROR (status);
#else
  g_error ("ncm_sphere_map_load_from_fits_catalog: CFITSIO not available.");

#endif
}

static void
_ncm_sphere_map_prepare_fft (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  guint fftw_default_flags         = ncm_cfg_get_fftw_default_flag ();

  if (self->fft_plan_r2c->len == 0)
  {
    const gint64 npix      = ncm_sphere_map_get_npix (smap);
    const gint64 nring_cap = ncm_sphere_map_get_nrings_cap (smap);
    gpointer temp_pix      = _fft_vec_alloc (self->npix);
    gboolean first;
    gint r_i;

    _fft_vec_set_zero_complex (self->fft_pvec, npix);

    _fft_vec_memcpy (temp_pix, self->pvec, self->npix);

    first = ncm_cfg_fftw_plan_begin ("ncm_sphere_map_rings_%ld", ncm_sphere_map_get_nside (smap));

    for (r_i = 0; r_i < nring_cap; r_i++)
    {
      const gint ring_size       = ncm_sphere_map_get_ring_size (smap, r_i);
      const gint64 ring_fi_north = ncm_sphere_map_get_ring_first_index (smap, r_i);
      const gint64 ring_fi_south = ncm_sphere_map_get_ring_first_index (smap, ncm_sphere_map_get_nrings (smap) - r_i - 1);
      const gint64 dist          = ring_fi_south - ring_fi_north;
      gdouble *pvec              = self->pvec;
      complex double *fft_pvec   = self->fft_pvec;

      fftw_plan plan_r2c = fftw_plan_many_dft_r2c (1, &ring_size, 2,
                                                   &pvec[ring_fi_north], NULL,
                                                   1, dist,
                                                   &fft_pvec[ring_fi_north], NULL,
                                                   1, dist,
                                                   fftw_default_flags | FFTW_PRESERVE_INPUT);

      fftw_plan plan_c2r = fftw_plan_many_dft_c2r (1, &ring_size, 2,
                                                   &fft_pvec[ring_fi_north], NULL,
                                                   1, dist,
                                                   &pvec[ring_fi_north], NULL,
                                                   1, dist,
                                                   fftw_default_flags | FFTW_DESTROY_INPUT);

      g_ptr_array_add (self->fft_plan_r2c, plan_r2c);
      g_ptr_array_add (self->fft_plan_c2r, plan_c2r);
    }

    {
      const gint ring_size     = self->middle_rings_size;
      const gint nrings_middle = ncm_sphere_map_get_nrings_middle (smap);
      const gint cap_size      = ncm_sphere_map_get_cap_size (smap);

      gdouble *pvec            = self->pvec;
      complex double *fft_pvec = self->fft_pvec;

      fftw_plan plan_r2c = fftw_plan_many_dft_r2c (1, &ring_size, nrings_middle,
                                                   &pvec[cap_size], NULL,
                                                   1, ring_size,
                                                   &fft_pvec[cap_size], NULL,
                                                   1, ring_size,
                                                   fftw_default_flags | FFTW_PRESERVE_INPUT);
      fftw_plan plan_c2r = fftw_plan_many_dft_c2r (1, &ring_size, nrings_middle,
                                                   &fft_pvec[cap_size], NULL,
                                                   1, ring_size,
                                                   &pvec[cap_size], NULL,
                                                   1, ring_size,
                                                   fftw_default_flags | FFTW_DESTROY_INPUT);

      g_ptr_array_add (self->fft_plan_r2c, plan_r2c);
      g_ptr_array_add (self->fft_plan_c2r, plan_c2r);
    }

    ncm_cfg_fftw_plan_end (first);


    _fft_vec_memcpy (self->pvec, temp_pix, self->npix);
    _fft_vec_free (temp_pix);
  }
}

#include "ncm/sphere/ncm_sphere_map_block.c"

static void
_ncm_sphere_map_map2alm_calc_Cl (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  guint m, l, lm_index = 0;

  ncm_vector_set_zero (self->Cl);

  m = 0;

  for (l = m; l <= self->lmax; l++)
  {
    const gdouble Re_alm = creal (self->alm[lm_index]);
    const gdouble Im_alm = cimag (self->alm[lm_index]);

    lm_index++;
    ncm_vector_fast_addto (self->Cl, l, 1.0 * (Re_alm * Re_alm + Im_alm * Im_alm));
  }

  for (m = 1; m <= self->lmax; m++)
  {
    for (l = m; l <= self->lmax; l++)
    {
      const gdouble Re_alm = creal (self->alm[lm_index]);
      const gdouble Im_alm = cimag (self->alm[lm_index]);

      lm_index++;
      ncm_vector_fast_addto (self->Cl, l, 2.0 * (Re_alm * Re_alm + Im_alm * Im_alm));
    }
  }

  for (l = 0; l <= self->lmax; l++)
  {
    ncm_vector_fast_mulby (self->Cl, l, 1.0 / (2.0 * l + 1.0));
  }

  self->has_Cls = TRUE;
}

/**
 * ncm_sphere_map_prepare_alm:
 * @smap: a #NcmSphereMap
 *
 * Calculates the $a_{\ell{}m}$ from the map @smap, using $\ell_\mathrm{max}$
 * set by ncm_sphere_map_set_lmax(). If $\ell_\mathrm{max} = 0$
 * nothing is done.
 *
 */
void
ncm_sphere_map_prepare_alm (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  guint i;

  if (self->lmax == 0)
  {
    g_warning ("ncm_sphere_map_prepare_alm: lmax equal to zero, returning...");

    return;
  }

#ifdef _NCM_SPHERE_MAP_MEASURE
  printf ("# Optimization control NC=%d STEP=%d CM=%d\n", NCM_SPHERE_MAP_BLOCK_NC, NCM_SPHERE_MAP_BLOCK_STEP, NCM_SPHERE_MAP_BLOCK_CM);
  printf ("# Preparing ffts!\n");
  fflush (stdout);
  ncm_timer_start (self->t);
#endif /* _NCM_SPHERE_MAP_MEASURE */
  _ncm_sphere_map_prepare_fft (smap);

  ncm_sphere_map_set_order (smap, NCM_SPHERE_MAP_ORDER_RING);

#ifdef _NCM_SPHERE_MAP_MEASURE
  printf ("# preparing fft plans, elapsed % 22.15g\n", ncm_timer_elapsed (self->t));
  printf ("# Performing ffts!\n");
  fflush (stdout);
  ncm_timer_start (self->t);
#endif /* _NCM_SPHERE_MAP_MEASURE */

  for (i = 0; i < self->fft_plan_r2c->len; i++)
  {
    fftw_execute (g_ptr_array_index (self->fft_plan_r2c, i));
  }

#ifdef _NCM_SPHERE_MAP_MEASURE
  printf ("# Performing ffts, elapsed % 22.15g\n", ncm_timer_elapsed (self->t));
  printf ("# Transforming rings\n");
  fflush (stdout);
  ncm_timer_start (self->t);
#endif /* _NCM_SPHERE_MAP_MEASURE */

  NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_map2alm_run) (smap);
#ifdef _NCM_SPHERE_MAP_MEASURE
  printf ("# %ld rings transformed, elapsed % 22.15g\n", ncm_sphere_map_get_nrings (smap), ncm_timer_elapsed (self->t));
#endif /* _NCM_SPHERE_MAP_MEASURE */

  /* Iterative refinement if iter > 0 */
  if (self->iter > 0)
  {
    gpointer original_map   = _fft_vec_alloc (self->npix);
    _fft_complex *saved_alm = _fft_vec_alloc_complex (self->alm_len);
    guint iter_i;

    /* Save original map */
    _fft_vec_memcpy (original_map, self->pvec, self->npix);

#ifdef _NCM_SPHERE_MAP_MEASURE
    printf ("# Performing %u iterations\n", self->iter);
    fflush (stdout);
#endif /* _NCM_SPHERE_MAP_MEASURE */

    for (iter_i = 0; iter_i < self->iter; iter_i++)
    {
#ifdef _NCM_SPHERE_MAP_MEASURE
      ncm_timer_start (self->t);
#endif /* _NCM_SPHERE_MAP_MEASURE */

      /* Save current alm */
      memcpy (saved_alm, self->alm, sizeof (_fft_complex) * self->alm_len);

      /* Synthesize: alm -> map */
      NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_alm2map_run) (smap);

      /* Complete synthesis with inverse FFT */
      for (i = 0; i < self->fft_plan_c2r->len; i++)
      {
        fftw_execute (g_ptr_array_index (self->fft_plan_c2r, i));
      }

      /* Compute residual: original_map - synthesized_map */
      for (i = 0; i < self->npix; i++)
      {
        _fft_vec_idx (self->pvec, i) = _fft_vec_idx (original_map, i) - _fft_vec_idx (self->pvec, i);
      }

      /* Transform residual -> alm_correction with forward FFT */
      for (i = 0; i < self->fft_plan_r2c->len; i++)
      {
        fftw_execute (g_ptr_array_index (self->fft_plan_r2c, i));
      }

      /* Compute correction alms from residual */
      NCM_SPHERE_MAP_BLOCK_DEC (_ncm_sphere_map_map2alm_run) (smap);

      /* Accumulate: alm = saved_alm + correction */
      for (i = 0; i < self->alm_len; i++)
      {
        self->alm[i] += saved_alm[i];
      }

#ifdef _NCM_SPHERE_MAP_MEASURE
      printf ("# Iteration %u completed, elapsed % 22.15g\n", iter_i + 1, ncm_timer_elapsed (self->t));
#endif /* _NCM_SPHERE_MAP_MEASURE */
    }

    /* Restore original map */
    _fft_vec_memcpy (self->pvec, original_map, self->npix);

    /* Clean up */
    _fft_vec_free (original_map);
    _fft_vec_free (saved_alm);
  }

  _ncm_sphere_map_map2alm_calc_Cl (smap);
}

/**
 * ncm_sphere_map_update_Cl:
 * @smap: a #NcmSphereMap
 *
 * Updates the values of $C_\ell$ based on the current $a_{lm}$.
 *
 */
void
ncm_sphere_map_update_Cl (NcmSphereMap *smap)
{
  _ncm_sphere_map_map2alm_calc_Cl (smap);
}

/**
 * ncm_sphere_map_compute_cross_Cl:
 * @smap1: first #NcmSphereMap with computed $a_{lm}$
 * @smap2: second #NcmSphereMap with computed $a_{lm}$
 *
 * Computes the cross-spectrum between two maps:
 * $$C_\ell^\mathrm{cross} = \frac{1}{2\ell+1} \sum_{m=-\ell}^\ell \mathrm{Re}(a_{lm}^{(1)} \bar{a}_{lm}^{(2)})$$
 *
 * Both maps must have their spherical harmonic coefficients already computed via
 * ncm_sphere_map_prepare_alm(), and must have the same lmax and nside.
 *
 * Returns: (transfer full): a #NcmVector containing the cross-spectrum $C_\ell$
 * for $\ell = 0, \ldots, \ell_\mathrm{max}$.
 */
NcmVector *
ncm_sphere_map_compute_cross_Cl (NcmSphereMap *smap1, NcmSphereMap *smap2)
{
  NcmSphereMapPrivate * const self1 = ncm_sphere_map_get_instance_private (smap1);
  NcmSphereMapPrivate * const self2 = ncm_sphere_map_get_instance_private (smap2);
  guint m, l, lm_index = 0;
  NcmVector *cross_Cl;

  g_assert (self1->lmax == self2->lmax);
  g_assert (self1->nside == self2->nside);
  g_assert (self1->alm != NULL);
  g_assert (self2->alm != NULL);

  cross_Cl = ncm_vector_new (self1->lmax + 1);
  ncm_vector_set_zero (cross_Cl);

  /* m = 0 term */
  m = 0;

  for (l = m; l <= self1->lmax; l++)
  {
    const gdouble Re_alm1 = creal (self1->alm[lm_index]);
    const gdouble Re_alm2 = creal (self2->alm[lm_index]);

    /* For m=0, a_lm is real, so cross term is just product */
    lm_index++;
    ncm_vector_fast_addto (cross_Cl, l, 1.0 * Re_alm1 * Re_alm2);
  }

  /* m > 0 terms */
  for (m = 1; m <= self1->lmax; m++)
  {
    for (l = m; l <= self1->lmax; l++)
    {
      const gdouble Re_alm1 = creal (self1->alm[lm_index]);
      const gdouble Im_alm1 = cimag (self1->alm[lm_index]);
      const gdouble Re_alm2 = creal (self2->alm[lm_index]);
      const gdouble Im_alm2 = cimag (self2->alm[lm_index]);

      /* Real part of alm1 * conj(alm2) = Re1*Re2 + Im1*Im2 */
      /* Factor of 2 accounts for both +m and -m terms */
      lm_index++;
      ncm_vector_fast_addto (cross_Cl, l, 2.0 * (Re_alm1 * Re_alm2 + Im_alm1 * Im_alm2));
    }
  }

  /* Normalize by (2*ell + 1) */
  for (l = 0; l <= self1->lmax; l++)
  {
    ncm_vector_fast_mulby (cross_Cl, l, 1.0 / (2.0 * l + 1.0));
  }

  return cross_Cl;
}

/**
 * ncm_sphere_map_get_alm:
 * @smap: a #NcmSphereMap
 * @l: value of $l < \ell_\mathrm{max}$
 * @m: value of $m \leq l$.
 * @Re_alm: (out): real part of $a_{lm}$
 * @Im_alm: (out): imaginary part of $a_{lm}$
 *
 * Gets the value of $a_{lm}$ previously calculated by
 * ncm_sphere_map_prepare_alm().
 *
 */
void
ncm_sphere_map_get_alm (NcmSphereMap *smap, guint l, guint m, gdouble *Re_alm, gdouble *Im_alm)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  gint lm_index                    = NCM_SPHERE_MAP_ALM_INDEX (self->lmax, l, m);

  /*Re_alm[0] = ncm_vector_fast_get (self->alm, lm_index + 0);*/
  /*Im_alm[0] = ncm_vector_fast_get (self->alm, lm_index + 1);*/
  Re_alm[0] = creal (self->alm[lm_index]);
  Im_alm[0] = cimag (self->alm[lm_index]);
}

/**
 * ncm_sphere_map_set_alm:
 * @smap: a #NcmSphereMap
 * @l: value of $l < \ell_\mathrm{max}$
 * @m: value of $m \leq l$.
 * @Re_alm: real part of $a_{lm}$
 * @Im_alm: imaginary part of $a_{lm}$
 *
 * Sets the value of $a_{lm}$.
 *
 */
void
ncm_sphere_map_set_alm (NcmSphereMap *smap, guint l, guint m, gdouble Re_alm, gdouble Im_alm)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  gint lm_index                    = NCM_SPHERE_MAP_ALM_INDEX (self->lmax, l, m);

  self->alm[lm_index] = Re_alm + I * Im_alm;
}

/**
 * ncm_sphere_map_get_Cl:
 * @smap: a #NcmSphereMap
 * @l: value of $l < \ell_\mathrm{max}$
 *
 * Gets the value of $C_{\ell}$ previously calculated by
 * ncm_sphere_map_prepare_alm().
 *
 */
gdouble
ncm_sphere_map_get_Cl (NcmSphereMap *smap, guint l)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  return ncm_vector_fast_get (self->Cl, l);
}

/**
 * ncm_sphere_map_get_pix:
 * @smap: a #NcmSphereMap
 * @i: pixel index
 *
 * Gets the value of pixel index by @i.
 *
 */
gdouble
ncm_sphere_map_get_pix (NcmSphereMap *smap, guint i)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  return _fft_vec_idx (self->pvec, i);
}

/**
 * ncm_sphere_map_add_noise:
 * @smap: a #NcmSphereMap
 * @sd: noise standard deviation
 * @rng: a #NcmRNG
 *
 * Adds a Gaussian noise with $\sigma=$ @sd and zero mean to each pixel.
 *
 */
void
ncm_sphere_map_add_noise (NcmSphereMap *smap, const gdouble sd, NcmRNG *rng)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  guint i;

  ncm_rng_lock (rng);

  for (i = 0; i < self->npix; i++)
  {
    const gdouble n_i = ncm_rng_gaussian_gen (rng, 0.0, sd);

    _fft_vec_idx (self->pvec, i) += n_i;
  }

  ncm_rng_unlock (rng);
}

/**
 * ncm_sphere_map_set_map:
 * @smap: a #NcmSphereMap
 * @map: (array) (element-type gdouble): pixels
 *
 * Set map pixels to @map using current ordering.
 *
 */
void
ncm_sphere_map_set_map (NcmSphereMap *smap, GArray *map)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  guint i;

  g_assert_cmpuint (map->len, ==, self->npix);

  for (i = 0; i < self->npix; i++)
  {
    _fft_vec_idx (self->pvec, i) = g_array_index (map, gdouble, i);
  }
}

/**
 * ncm_sphere_map_set_Cls:
 * @smap: a #NcmSphereMap
 * @Cls: a #NcmVector containing the $C_\ell$
 *
 * Set map $C_l$s.
 *
 */
void
ncm_sphere_map_set_Cls (NcmSphereMap *smap, NcmVector *Cls)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  g_assert_cmpint (self->lmax, >, 0);

  ncm_vector_memcpy2 (self->Cl, Cls, 0, 0, self->lmax + 1);

  self->has_Cls = TRUE;
}

/**
 * ncm_sphere_map_alm2map:
 * @smap: a #NcmSphereMap
 *
 * Compute map pixels from current $a_{\ell{}m}$.
 *
 */
void
ncm_sphere_map_alm2map (NcmSphereMap *smap)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);
  guint i;

  /*gfloat *temp_pix = _fft_vec_alloc (self->npix);*/
  /*gfloat *pixels   = self->pvec;*/
  /*_fft_vec_memcpy (temp_pix, self->pvec, self->npix);*/

  g_assert_cmpuint (self->nside, >, 0);

  if (self->lmax == 0)
  {
    g_warning ("ncm_sphere_map_prepare_alm: lmax equal to zero, returning...");

    return;
  }

#ifdef _NCM_SPHERE_MAP_MEASURE
  printf ("# Preparing ffts!\n");
  fflush (stdout);
  ncm_timer_start (self->t);
#endif /* _NCM_SPHERE_MAP_MEASURE */

  _ncm_sphere_map_prepare_fft (smap);

  self->order = NCM_SPHERE_MAP_ORDER_RING;

#ifdef _NCM_SPHERE_MAP_MEASURE
  printf ("# preparing fft plans, elapsed % 22.15g\n", ncm_timer_elapsed (self->t));
  printf ("# Transforming rings\n");
  fflush (stdout);
  ncm_timer_start (self->t);
#endif /* _NCM_SPHERE_MAP_MEASURE */

  NCM_SPHERE_MAP_BLOCK_INV_DEC (_ncm_sphere_map_alm2map_run) (smap);
#ifdef _NCM_SPHERE_MAP_MEASURE
  printf ("# %ld rings transformed, elapsed % 22.15g\n", ncm_sphere_map_get_nrings (smap), ncm_timer_elapsed (self->t));
#endif /* _NCM_SPHERE_MAP_MEASURE */

#ifdef _NCM_SPHERE_MAP_MEASURE
  printf ("# Peforming ffts!\n");
  fflush (stdout);
  ncm_timer_start (self->t);
#endif /* _NCM_SPHERE_MAP_MEASURE */

  for (i = 0; i < self->fft_plan_c2r->len; i++)
  {
    fftw_execute (g_ptr_array_index (self->fft_plan_c2r, i));
  }

#ifdef _NCM_SPHERE_MAP_MEASURE
  printf ("# Peforming ffts, elapsed % 22.15g\n", ncm_timer_elapsed (self->t));
#endif /* _NCM_SPHERE_MAP_MEASURE */
}

static gdouble
_ncm_sphere_map_calc_Ctheta_theta (const gdouble theta, gpointer userdata)
{
  NcmSphereMapPrivate * const self = (NcmSphereMapPrivate * const ) userdata;
  const gdouble x                  = cos (theta);
  gdouble p_ellm2                  = 1.0;
  gdouble p_ellm1                  = x;
  gdouble p_ell;
  guint ell;

  gdouble Ctheta =
    1.0 * (2.0 * 0.0 + 1.0) * ncm_vector_fast_get (self->Cl, 0) +
    x   * (2.0 * 1.0 + 1.0) * ncm_vector_fast_get (self->Cl, 1);

  for (ell = 2; ell <= self->lmax; ell++)
  {
    p_ell   = (x * (2.0 * ell - 1.0) * p_ellm1 - (ell - 1.0) * p_ellm2) / ell;
    p_ellm2 = p_ellm1;
    p_ellm1 = p_ell;

    Ctheta += p_ell * (2.0 * ell + 1.0) * ncm_vector_fast_get (self->Cl, ell);
  }

  return Ctheta / (4.0 * ncm_c_pi ());
}

/**
 * ncm_sphere_map_calc_Ctheta:
 * @smap: a #NcmSphereMap
 * @reltol: required tolerance for $C(\theta)$
 *
 * Computes the two-point correlation function $C(\theta)$ from
 * the precomputed $C_\ell$.
 *
 * Returns: (transfer full): the $C(\theta)$ spline.
 */
NcmSpline *
ncm_sphere_map_calc_Ctheta (NcmSphereMap *smap, const gdouble reltol)
{
  NcmSphereMapPrivate * const self = ncm_sphere_map_get_instance_private (smap);

  if (!self->has_Cls)
  {
    g_error ("ncm_sphere_map_calc_Ctheta: object does not contain Cls.");

    return NULL;
  }
  else
  {
    NcmSpline *s = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
    gsl_function F;

    F.params   = self;
    F.function = &_ncm_sphere_map_calc_Ctheta_theta;

    ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, 0.0, M_PI, -1, reltol);

    return s;
  }
}

