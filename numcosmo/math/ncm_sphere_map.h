/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_sphere_map.h
 *
 *  Wed Jul  9 11:11:38 2008 (updated Jul 15 2016)
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_sphere_map.h
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

#ifndef _NCM_SPHERE_MAP_H_
#define _NCM_SPHERE_MAP_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_quaternion.h>
#include <numcosmo/math/ncm_spline.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_vector_float.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_SPHERE_MAP             (ncm_sphere_map_get_type ())
#define NCM_SPHERE_MAP(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SPHERE_MAP, NcmSphereMap))
#define NCM_SPHERE_MAP_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SPHERE_MAP, NcmSphereMapClass))
#define NCM_IS_SPHERE_MAP(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SPHERE_MAP))
#define NCM_IS_SPHERE_MAP_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SPHERE_MAP))
#define NCM_SPHERE_MAP_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SPHERE_MAP, NcmSphereMapClass))

typedef struct _NcmSphereMapClass NcmSphereMapClass;
typedef struct _NcmSphereMap NcmSphereMap;
typedef struct _NcmSphereMapPrivate NcmSphereMapPrivate;

struct _NcmSphereMapClass
{
  /*< private >*/
  GObjectClass parent_class;
};

/**
 * NcmSphereMapOrder:
 * @NCM_SPHERE_MAP_ORDER_NEST: FIXME
 * @NCM_SPHERE_MAP_ORDER_RING: FIXME
 * 
 * FIXME
 * 
 */
typedef enum _NcmSphereMapOrder
{
  NCM_SPHERE_MAP_ORDER_NEST,
  NCM_SPHERE_MAP_ORDER_RING, 
  /* < private > */
  NCM_SPHERE_MAP_ORDER_LEN,  /*< skip >*/
} NcmSphereMapOrder;

/**
 * NcmSphereMapCoordSys:
 * @NCM_SPHERE_MAP_COORD_SYS_GALACTIC: FIXME
 * @NCM_SPHERE_MAP_COORD_SYS_ECLIPTIC: FIXME
 * @NCM_SPHERE_MAP_COORD_SYS_CELESTIAL: FIXME
 * 
 * FIXME
 * 
 */
typedef enum _NcmSphereMapCoordSys
{
  NCM_SPHERE_MAP_COORD_SYS_GALACTIC  = 'G',
  NCM_SPHERE_MAP_COORD_SYS_ECLIPTIC  = 'E',
  NCM_SPHERE_MAP_COORD_SYS_CELESTIAL = 'C', 
  /* < private > */
  NCM_SPHERE_MAP_COORD_SYS_LEN, /*< skip >*/
} NcmSphereMapCoordSys;

struct _NcmSphereMap
{
  /*< private >*/
  GObject parent_instance;
  NcmSphereMapPrivate *priv;
};

GType ncm_sphere_map_get_type (void) G_GNUC_CONST;

NcmSphereMap *ncm_sphere_map_new (const gint64 nside);
NcmSphereMap *ncm_sphere_map_ref (NcmSphereMap *smap);
void ncm_sphere_map_free (NcmSphereMap *smap);
void ncm_sphere_map_clear (NcmSphereMap **smap);

void ncm_sphere_map_set_nside (NcmSphereMap *smap, gint64 nside);
gint64 ncm_sphere_map_get_nside (NcmSphereMap *smap);
gint64 ncm_sphere_map_get_npix (NcmSphereMap *smap);
gint64 ncm_sphere_map_get_nsmap (NcmSphereMap *smap);
gint64 ncm_sphere_map_get_cap_size (NcmSphereMap *smap);
gint64 ncm_sphere_map_get_middle_size (NcmSphereMap *smap);

gint64 ncm_sphere_map_get_nrings (NcmSphereMap *smap);
gint64 ncm_sphere_map_get_nrings_cap (NcmSphereMap *smap);
gint64 ncm_sphere_map_get_nrings_middle (NcmSphereMap *smap);
gint64 ncm_sphere_map_get_ring_size (NcmSphereMap *smap, gint64 r_i);
gint64 ncm_sphere_map_get_ring_first_index (NcmSphereMap *smap, gint64 r_i);

void ncm_sphere_map_set_order (NcmSphereMap *smap, NcmSphereMapOrder order);
NcmSphereMapOrder ncm_sphere_map_get_order (NcmSphereMap *smap);

void ncm_sphere_map_set_coordsys (NcmSphereMap *smap, NcmSphereMapCoordSys coordsys);
NcmSphereMapCoordSys ncm_sphere_map_get_coordsys (NcmSphereMap *smap);

void ncm_sphere_map_set_lmax (NcmSphereMap *smap, guint lmax);
guint ncm_sphere_map_get_lmax (NcmSphereMap *smap);

void ncm_sphere_map_clear_smapels (NcmSphereMap *smap);

gint64 ncm_sphere_map_nest2ring (NcmSphereMap *smap, const gint64 nest_index);
gint64 ncm_sphere_map_ring2nest (NcmSphereMap *smap, const gint64 ring_index);
void ncm_sphere_map_pix2ang_nest (NcmSphereMap *smap, const gint64 nest_index, gdouble *theta, gdouble *phi);
void ncm_sphere_map_pix2ang_ring (NcmSphereMap *smap, const gint64 ring_index, gdouble *theta, gdouble *phi);
void ncm_sphere_map_pix2vec_nest (NcmSphereMap *smap, const gint64 nest_index, NcmTriVec *vec);
void ncm_sphere_map_pix2vec_ring (NcmSphereMap *smap, const gint64 ring_index, NcmTriVec *vec);
void ncm_sphere_map_ang2pix_nest (NcmSphereMap *smap, const gdouble theta, const gdouble phi, gint64 *nest_index);
void ncm_sphere_map_ang2pix_ring (NcmSphereMap *smap, const gdouble theta, const gdouble phi, gint64 *ring_index);

void ncm_sphere_map_vec2pix_ring (NcmSphereMap *smap, NcmTriVec *vec, gint64 *ring_index);
void ncm_sphere_map_vec2pix_nest (NcmSphereMap *smap, NcmTriVec *vec, gint64 *nest_index);

void ncm_sphere_map_add_to_vec (NcmSphereMap *smap, NcmTriVec *vec, const gdouble s);
void ncm_sphere_map_add_to_ang (NcmSphereMap *smap, const gdouble theta, const gdouble phi, const gdouble s);

void ncm_sphere_map_load_fits (NcmSphereMap *smap, const gchar *fits_file, const gchar *signal_name);
void ncm_sphere_map_save_fits (NcmSphereMap *smap, const gchar *fits_file, const gchar *signal_name, gboolean overwrite);

void ncm_sphere_map_load_from_fits_catalog (NcmSphereMap *smap, const gchar *fits_file, const gchar *RA, const gchar *DEC, const gchar *S);

void ncm_sphere_map_prepare_alm (NcmSphereMap *smap);
void ncm_sphere_map_update_Cl (NcmSphereMap *smap);

void ncm_sphere_map_get_alm (NcmSphereMap *smap, guint l, guint m, gdouble *Re_alm, gdouble *Im_alm);
void ncm_sphere_map_set_alm (NcmSphereMap *smap, guint l, guint m, gdouble Re_alm, gdouble Im_alm);

gdouble ncm_sphere_map_get_Cl (NcmSphereMap *smap, guint l);
gdouble ncm_sphere_map_get_pix (NcmSphereMap *smap, guint i);

void ncm_sphere_map_add_noise (NcmSphereMap *smap, const gdouble sd, NcmRNG *rng);
void ncm_sphere_map_set_map (NcmSphereMap *smap, GArray *map);
void ncm_sphere_map_set_Cls (NcmSphereMap *smap, NcmVector *Cls);

void ncm_sphere_map_alm2map (NcmSphereMap *smap);

NcmSpline *ncm_sphere_map_calc_Ctheta (NcmSphereMap *smap, const gdouble reltol);

#define NCM_SPHERE_MAP_N(nside) (12 * (nside) * (nside))
#define NCM_SPHERE_MAP_INT_TO_XY(i,x,y) \
G_STMT_START { \
  gint shift = 0, shifted = i; \
  x = y = 0; \
  do { \
    x |= ((shifted & 1) << shift); \
    y |= (((shifted & 2) >> 1) << shift); \
    shift++; \
  } while (shifted >>= 2); \
} G_STMT_END

#define NCM_SPHERE_MAP_XY_TO_INT(x,y,i) \
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

#define NCM_SPHERE_MAP_HEALPIX_NULLVAL (-1.6375e30)
#define NCM_SPHERE_MAP_DEFAULT_SIGNAL "SIGNAL"

#define NCM_SPHERE_MAP_ALM_SIZE(lmax) (((lmax)*(lmax) + 3*(lmax) + 2)/2) 
#define NCM_SPHERE_MAP_M_START(lmax,m) ((2*(lmax)*(m)-(m)*(m)+3*(m))/2)

/*#define NCM_SPHERE_MAP_ALM_INDEX(l,m) ((l * (l + 1)) / 2 + m)*/
#define NCM_SPHERE_MAP_ALM_INDEX(lmax,l,m) ((l) + ((1 + 2 * (lmax) - (m)) * (m)) / 2)

G_END_DECLS

#endif /* _NCM_SPHERE_MAP_H_ */
