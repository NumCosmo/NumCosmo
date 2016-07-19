/***************************************************************************
 *            ncm_sphere_map.h
 *
 *  Tue Jun 24 16:35:23 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_sphere_map.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_matrix_complex_double.h>

#ifndef NUMCOSMO_GIR_SCAN
#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */
#endif

G_BEGIN_DECLS

#define NCM_TYPE_SPHERE_MAP             (ncm_sphere_map_get_type ())
#define NCM_SPHERE_MAP(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SPHERE_MAP, NcmSphereMap))
#define NCM_SPHERE_MAP_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SPHERE_MAP, NcmSphereMapClass))
#define NCM_IS_SPHERE_MAP(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SPHERE_MAP))
#define NCM_IS_SPHERE_MAP_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SPHERE_MAP))
#define NCM_SPHERE_MAP_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SPHERE_MAP, NcmSphereMapClass))

typedef struct _NcmSphereMapClass NcmSphereMapClass;
typedef struct _NcmSphereMap NcmSphereMap;

struct _NcmSphereMapClass
{
  /*< private >*/
  GObjectClass parent_class;
};

/**
 * NcmSphereMapOrder:
 * @NC_SPHERE_MAP_ORDER_NEST: FIXME
 * @NC_SPHERE_MAP_ORDER_RING: FIXME
 */
typedef enum _NcmSphereMapOrder
{
  NC_SPHERE_MAP_ORDER_NEST,
  NC_SPHERE_MAP_ORDER_RING
} NcmSphereMapOrder;

/**
 * NcmSphereMapType:
 * @NC_SPHERE_MAP_TYPE_TEMPERATURE: FIXME
 * @NC_SPHERE_MAP_TYPE_Q_POLARIZATION: FIXME
 * @NC_SPHERE_MAP_TYPE_U_POLARISATION: FIXME
 * @NC_SPHERE_MAP_TYPE_SPUR_SIGNAL: FIXME
 * @NC_SPHERE_MAP_TYPE_N_OBS: FIXME
 * 
 * FIXME
 */
typedef enum _NcmSphereMapType
{
  NC_SPHERE_MAP_TYPE_TEMPERATURE    = 1 << 0,
  NC_SPHERE_MAP_TYPE_Q_POLARIZATION = 1 << 1,
  NC_SPHERE_MAP_TYPE_U_POLARISATION = 1 << 2,
  NC_SPHERE_MAP_TYPE_SPUR_SIGNAL    = 1 << 3,
  NC_SPHERE_MAP_TYPE_N_OBS          = 1 << 4
} NcmSphereMapType;

/**
 * NcmSphereMap:
 * 
 * FIXME
 * 
 */
struct _NcmSphereMap
{
  /*< private >*/
  GObject parent_instance;
  gsl_vector_float *dt;
  gsl_vector_float *qpol;
  gsl_vector_float *upol;
  gsl_vector_float *spur_signal;
  gsl_vector_float *nobs;
  NcmSphereMapType type;
  glong npix;
  glong nside;
  glong nrings;
  gboolean loaded;
  NcmSphereMapOrder order;
  gsl_vector *theta;
  gsl_vector *phi;
  gboolean is_init_coord;
};

typedef struct _NcmSphereMapAlm NcmSphereMapAlm;

/**
 * NcmSphereMapAlm:
 * 
 * FIXME
 */
struct _NcmSphereMapAlm
{
  /*< private >*/
  gint lmax;
  gint alm_size;
  gboolean loaded;
  gsl_vector_complex *alm;
  gsl_vector *Nc;
  gsl_vector *sqrt_int;
  gsl_vector *lnpoch_m_1_2;
  gsl_vector *sphPlm_recur1;
  gsl_vector *sphPlm_recur2;
  gsl_vector *sphPmm;
};

#ifdef NUMCOSMO_HAVE_FFTW3

typedef struct _NcmSphereMapSHT NcmSphereMapSHT;

/**
 * NcmSphereMapSHT:
 * 
 * FIXME
 */
struct _NcmSphereMapSHT
{
  /*< private >*/
  NcmSphereMap *map;
  NcmSphereMapAlm *mapalm;
  gint n_rings;
  gint n_diff_rings;
  gint max_ring_size;
  fftw_plan *forward_plans;
  fftw_plan *backward_plans;
  gboolean save_wis;
  gchar *wis_file;
  gchar *lmin_file;
  gsl_matrix_complex *fft_ring;
  gsl_vector *ring;
  gsl_vector_complex *sphPlm;
  gsl_vector *sphPlm_upper_limit;
  gsl_vector *sphPl0;
};
#endif /* NUMCOSMO_HAVE_FFTW3 */

GType ncm_sphere_map_get_type (void) G_GNUC_CONST;

NcmSphereMap *ncm_sphere_map_new (guint nside);
NcmSphereMap *ncm_sphere_map_clone (NcmSphereMap *map);
gboolean ncm_sphere_map_copy (NcmSphereMap *dest, NcmSphereMap *orig);
gboolean ncm_sphere_map_init_coord (NcmSphereMap *map);
gboolean ncm_sphere_map_set_order (NcmSphereMap *map, NcmSphereMapOrder order, gboolean init_coord);

NcmSphereMapAlm *ncm_sphere_mapalm_new (void);
gboolean ncm_sphere_mapalm_init (NcmSphereMapAlm *mapalm, gint lmax);

#ifdef NUMCOSMO_HAVE_FFTW3
NcmSphereMapSHT *ncm_sphere_mapsht_new (NcmSphereMap *map, NcmSphereMapAlm *mapalm, guint fftw_flags);
gboolean ncm_sphere_mapsht_map2alm_circle (NcmSphereMapSHT *mapsht, gint ring, gint ring_size, gdouble norma, gdouble theta, gdouble phi, gint start_m, gint end_m);
gboolean ncm_sphere_mapsht_alm2map_circle (NcmSphereMapSHT *mapsht, gint ring, gint ring_size, gdouble theta, gdouble phi);
gboolean ncm_sphere_mapsht_map2alm (NcmSphereMapSHT *mapsht, gdouble cut);
gboolean ncm_sphere_mapsht_alm2map (NcmSphereMapSHT *mapsht);
#endif /* NUMCOSMO_HAVE_FFTW3 */

gdouble ncm_sphere_map_homogenize_noise (NcmSphereMap *map, gdouble base_sigma);
gdouble ncm_sphere_map_rotate_avg (NcmSphereMap *map, glong n);

#define NCM_SPHERE_MAP_RNG_NAME "sphere_map"

#define NCM_MAP_ALM_SIZE(lmax) (((lmax)*(lmax) + 3*(lmax) + 2)/2)
#define NCM_MAP_N_IND_PLM(lmax) (((lmax)*(lmax) + 3*(lmax) + 2)/2)
#define NCM_MAP_MAX_RING_SIZE(nside) (4*(nside))
#define NCM_MAP_N_DIFFERENT_SIZED_RINGS(nside) (nside)
#define NCM_MAP_RING_PLAN_INDEX(nside,ring_n) (((ring_n) < (nside)) ? (ring_n) : ((ring_n)>=(3*(nside)) ? (4*(nside)-(ring_n)-2) : ((nside)-1)))
#define NCM_MAP_RING_SIZE(nside,ring_n) (4 * (NCM_MAP_RING_PLAN_INDEX (nside,ring_n) + 1))
#define NCM_MAP_N_RINGS(nside) (4*(nside)-1)
#define NCM_MAP_ALM_M_START(lmax,m) ((2*(lmax)*(m)-(m)*(m)+3*(m))/2)
#define NCM_MAP_ALM_INDEX(lmax,l,m) (((l) >= (m)) ? (NCM_MAP_ALM_M_START(lmax,m) + (l) - (m)) : (-1))

G_END_DECLS

#endif /* _NCM_SPHERE_MAP_H_ */
