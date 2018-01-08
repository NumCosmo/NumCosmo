/***************************************************************************
 *            ncm_sphere_map_pix.h
 *
 *  Wed Jul  9 11:11:38 2008 (updated Jul 15 2016)
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_sphere_map_pix.h
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

#ifndef _NCM_SPHERE_MAP_PIX_H_
#define _NCM_SPHERE_MAP_PIX_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_quaternion.h>
#include <numcosmo/math/ncm_sf_spherical_harmonics.h>
#include <numcosmo/math/ncm_timer.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_vector_float.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_SPHERE_MAP_PIX             (ncm_sphere_map_pix_get_type ())
#define NCM_SPHERE_MAP_PIX(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SPHERE_MAP_PIX, NcmSphereMapPix))
#define NCM_SPHERE_MAP_PIX_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SPHERE_MAP_PIX, NcmSphereMapPixClass))
#define NCM_IS_SPHERE_MAP_PIX(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SPHERE_MAP_PIX))
#define NCM_IS_SPHERE_MAP_PIX_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SPHERE_MAP_PIX))
#define NCM_SPHERE_MAP_PIX_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SPHERE_MAP_PIX, NcmSphereMapPixClass))

typedef struct _NcmSphereMapPixClass NcmSphereMapPixClass;
typedef struct _NcmSphereMapPix NcmSphereMapPix;

struct _NcmSphereMapPixClass
{
  /*< private >*/
  GObjectClass parent_class;
};

/**
 * NcmSphereMapPixOrder:
 * @NCM_SPHERE_MAP_PIX_ORDER_NEST: FIXME
 * @NCM_SPHERE_MAP_PIX_ORDER_RING: FIXME
 * 
 * FIXME
 * 
 */
typedef enum _NcmSphereMapPixOrder
{
  NCM_SPHERE_MAP_PIX_ORDER_NEST,
  NCM_SPHERE_MAP_PIX_ORDER_RING, 
  /* < private > */
  NCM_SPHERE_MAP_PIX_ORDER_LEN,  /*< skip >*/
} NcmSphereMapPixOrder;

/**
 * NcmSphereMapPixCoordSys:
 * @NCM_SPHERE_MAP_PIX_COORD_SYS_GALACTIC: FIXME
 * @NCM_SPHERE_MAP_PIX_COORD_SYS_ECLIPTIC: FIXME
 * @NCM_SPHERE_MAP_PIX_COORD_SYS_CELESTIAL: FIXME
 * 
 * FIXME
 * 
 */
typedef enum _NcmSphereMapPixCoordSys
{
  NCM_SPHERE_MAP_PIX_COORD_SYS_GALACTIC  = 'G',
  NCM_SPHERE_MAP_PIX_COORD_SYS_ECLIPTIC  = 'E',
  NCM_SPHERE_MAP_PIX_COORD_SYS_CELESTIAL = 'C', 
  /* < private > */
  NCM_SPHERE_MAP_PIX_COORD_SYS_LEN, /*< skip >*/
} NcmSphereMapPixCoordSys;

struct _NcmSphereMapPix
{
  /*< private >*/
  GObject parent_instance;
  gint64 nside;
  gint64 npix;
  gint64 face_size;
  gint64 middle_rings_size;
  gint64 cap_size;
  gint64 middle_size;
  gint64 nrings;
  gint64 nrings_cap;
  gint64 nrings_middle;
  NcmSphereMapPixOrder order;
  NcmSphereMapPixCoordSys coordsys;
  gpointer pvec;
  gpointer fft_pvec;
  GPtrArray *fft_plan_r2c;
  GPtrArray *fft_plan_c2r;
  guint lmax;
  GArray *alm;
  NcmVector *alm_v;
  NcmVector *Cl;
  NcmTimer *t;
	NcmSFSphericalHarmonics *spha;
};

GType ncm_sphere_map_pix_get_type (void) G_GNUC_CONST;

NcmSphereMapPix *ncm_sphere_map_pix_new (const gint64 nside);
NcmSphereMapPix *ncm_sphere_map_pix_ref (NcmSphereMapPix *pix);
void ncm_sphere_map_pix_free (NcmSphereMapPix *pix);
void ncm_sphere_map_pix_clear (NcmSphereMapPix **pix);

void ncm_sphere_map_pix_set_nside (NcmSphereMapPix *pix, gint64 nside);
gint64 ncm_sphere_map_pix_get_nside (NcmSphereMapPix *pix);
gint64 ncm_sphere_map_pix_get_npix (NcmSphereMapPix *pix);
gint64 ncm_sphere_map_pix_get_cap_size (NcmSphereMapPix *pix);
gint64 ncm_sphere_map_pix_get_middle_size (NcmSphereMapPix *pix);

gint64 ncm_sphere_map_pix_get_nrings (NcmSphereMapPix *pix);
gint64 ncm_sphere_map_pix_get_nrings_cap (NcmSphereMapPix *pix);
gint64 ncm_sphere_map_pix_get_nrings_middle (NcmSphereMapPix *pix);
gint64 ncm_sphere_map_pix_get_ring_size (NcmSphereMapPix *pix, gint64 r_i);
gint64 ncm_sphere_map_pix_get_ring_first_index (NcmSphereMapPix *pix, gint64 r_i);

void ncm_sphere_map_pix_set_order (NcmSphereMapPix *pix, NcmSphereMapPixOrder order);
NcmSphereMapPixOrder ncm_sphere_map_pix_get_order (NcmSphereMapPix *pix);

void ncm_sphere_map_pix_set_coordsys (NcmSphereMapPix *pix, NcmSphereMapPixCoordSys coordsys);
NcmSphereMapPixCoordSys ncm_sphere_map_pix_get_coordsys (NcmSphereMapPix *pix);

void ncm_sphere_map_pix_set_lmax (NcmSphereMapPix *pix, guint lmax);
guint ncm_sphere_map_pix_get_lmax (NcmSphereMapPix *pix);

void ncm_sphere_map_pix_clear_pixels (NcmSphereMapPix *pix);

gint64 ncm_sphere_map_pix_nest2ring (NcmSphereMapPix *pix, const gint64 nest_index);
gint64 ncm_sphere_map_pix_ring2nest (NcmSphereMapPix *pix, const gint64 ring_index);
void ncm_sphere_map_pix_pix2ang_nest (NcmSphereMapPix *pix, const gint64 nest_index, gdouble *theta, gdouble *phi);
void ncm_sphere_map_pix_pix2ang_ring (NcmSphereMapPix *pix, const gint64 ring_index, gdouble *theta, gdouble *phi);
void ncm_sphere_map_pix_pix2vec_nest (NcmSphereMapPix *pix, const gint64 nest_index, NcmTriVec *vec);
void ncm_sphere_map_pix_pix2vec_ring (NcmSphereMapPix *pix, const gint64 ring_index, NcmTriVec *vec);
void ncm_sphere_map_pix_ang2pix_nest (NcmSphereMapPix *pix, const gdouble theta, const gdouble phi, gint64 *nest_index);
void ncm_sphere_map_pix_ang2pix_ring (NcmSphereMapPix *pix, const gdouble theta, const gdouble phi, gint64 *ring_index);

void ncm_sphere_map_pix_vec2pix_ring (NcmSphereMapPix *pix, NcmTriVec *vec, gint64 *ring_index);
void ncm_sphere_map_pix_vec2pix_nest (NcmSphereMapPix *pix, NcmTriVec *vec, gint64 *nest_index);

void ncm_sphere_map_pix_add_to_vec (NcmSphereMapPix *pix, NcmTriVec *vec, const gdouble s);
void ncm_sphere_map_pix_add_to_ang (NcmSphereMapPix *pix, const gdouble theta, const gdouble phi, const gdouble s);

void ncm_sphere_map_pix_load_fits (NcmSphereMapPix *pix, const gchar *fits_file, const gchar *signal_name);
void ncm_sphere_map_pix_save_fits (NcmSphereMapPix *pix, const gchar *fits_file, const gchar *signal_name, gboolean overwrite);

void ncm_sphere_map_pix_load_from_fits_catalog (NcmSphereMapPix *pix, const gchar *fits_file, const gchar *RA, const gchar *DEC, const gchar *S);

void ncm_sphere_map_pix_prepare_alm (NcmSphereMapPix *pix);
void ncm_sphere_map_pix_prepare_Cl (NcmSphereMapPix *pix);

void ncm_sphere_map_pix_get_alm (NcmSphereMapPix *pix, guint l, guint m, gdouble *Re_alm, gdouble *Im_alm);
gdouble ncm_sphere_map_pix_get_Cl (NcmSphereMapPix *pix, guint l);

void ncm_sphere_map_pix_add_noise (NcmSphereMapPix *pix, const gdouble sd, NcmRNG *rng);

void ncm_sphere_map_pix_alm2map (NcmSphereMapPix *pix);

#define NCM_SPHERE_MAP_PIX_N(nside) (12 * (nside) * (nside))
#define NCM_SPHERE_MAP_PIX_INT_TO_XY(i,x,y) \
G_STMT_START { \
  gint shift = 0, shifted = i; \
  x = y = 0; \
  do { \
    x |= ((shifted & 1) << shift); \
    y |= (((shifted & 2) >> 1) << shift); \
    shift++; \
  } while (shifted >>= 2); \
} G_STMT_END

#define NCM_SPHERE_MAP_PIX_XY_TO_INT(x,y,i) \
G_STMT_START { \
  gint shift = 0, shifted_x = x, shifted_y = y; \
  g_assert (shifted_x >= 0 && shifted_y >= 0); \
  i = 0; \
  do { \
    i |= ((shifted_x & 1) << (shift + 0)); \
    i |= ((shifted_y & 1) << (shift + 1)); \
    shift += 2; shifted_x >>= 1; shifted_y >>= 1; \
  } while (shifted_x || shifted_y); \
} G_STMT_END

#define NCM_SPHERE_MAP_PIX_HEALPIX_NULLVAL (-1.6375e30)
#define NCM_SPHERE_MAP_PIX_DEFAULT_SIGNAL "SIGNAL"

#define NCM_SPHERE_MAP_PIX_ALM_SIZE(lmax) (((lmax)*(lmax) + 3*(lmax) + 2)/2) 
#define NCM_SPHERE_MAP_PIX_M_START(lmax,m) ((2*(lmax)*(m)-(m)*(m)+3*(m))/2)

/*#define NCM_SPHERE_MAP_PIX_ALM_INDEX(l,m) ((l * (l + 1)) / 2 + m)*/
#define NCM_SPHERE_MAP_PIX_ALM_INDEX(lmax,l,m) ((l) + ((1 + 2 * (lmax) - (m)) * (m)) / 2)

G_END_DECLS

#endif /* _NCM_SPHERE_MAP_PIX_H_ */
