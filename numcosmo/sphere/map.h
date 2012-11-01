/***************************************************************************
 *            map.h
 *
 *  Tue Jun 24 16:35:23 2008
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

#ifndef _NC_MAP_H
#define _NC_MAP_H

#include <glib.h>
#include <glib-object.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_matrix_complex_double.h>

#ifndef NUMCOSMO_GIR_SCAN
#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */
#endif

G_BEGIN_DECLS

/**
 * NcSphereMapOrder:
 * @NC_SPHERE_MAP_ORDER_NEST: FIXME
 * @NC_SPHERE_MAP_ORDER_RING: FIXME
 */
typedef enum _NcSphereMapOrder
{
  NC_SPHERE_MAP_ORDER_NEST,
  NC_SPHERE_MAP_ORDER_RING
} NcSphereMapOrder;

/**
 * NcSphereMapType:
 * @NC_SPHERE_MAP_TYPE_TEMPERATURE: FIXME
 * @NC_SPHERE_MAP_TYPE_Q_POLARIZATION: FIXME
 * @NC_SPHERE_MAP_TYPE_U_POLARISATION: FIXME
 * @NC_SPHERE_MAP_TYPE_SPUR_SIGNAL: FIXME
 * @NC_SPHERE_MAP_TYPE_N_OBS: FIXME
 * 
 * FIXME
 */
typedef enum _NcSphereMapType
{
  NC_SPHERE_MAP_TYPE_TEMPERATURE    = 1 << 0,
  NC_SPHERE_MAP_TYPE_Q_POLARIZATION = 1 << 1,
  NC_SPHERE_MAP_TYPE_U_POLARISATION = 1 << 2,
  NC_SPHERE_MAP_TYPE_SPUR_SIGNAL    = 1 << 3,
  NC_SPHERE_MAP_TYPE_N_OBS          = 1 << 4
} NcSphereMapType;


typedef struct _NcSphereMap NcSphereMap;

/**
 * NcSphereMap:
 * 
 * FIXME
 * 
 */
struct _NcSphereMap
{
  /*< private >*/
  gsl_vector_float *dt;
  gsl_vector_float *qpol;
  gsl_vector_float *upol;
  gsl_vector_float *spur_signal;
  gsl_vector_float *nobs;
  NcSphereMapType type;
  glong npix;
  glong nside;
  glong nrings;
  gboolean loaded;
  NcSphereMapOrder order;
  gsl_vector *theta;
  gsl_vector *phi;
  gboolean is_init_coord;
};

typedef struct _NcSphereMapAlm NcSphereMapAlm;

/**
 * NcSphereMapAlm:
 * 
 * FIXME
 */
struct _NcSphereMapAlm
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

typedef struct _NcSphereMapSHT NcSphereMapSHT;

/**
 * NcSphereMapSHT:
 * 
 * FIXME
 */
struct _NcSphereMapSHT
{
  /*< private >*/
  NcSphereMap *map;
  NcSphereMapAlm *mapalm;
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

NcSphereMap *nc_sphere_map_new (gint nside);
NcSphereMap *nc_sphere_map_clone (NcSphereMap *map);
gboolean nc_sphere_map_copy (NcSphereMap *dest, NcSphereMap *orig);
gboolean nc_sphere_map_init_coord (NcSphereMap *map);
gboolean nc_sphere_map_set_order (NcSphereMap *map, NcSphereMapOrder order, gboolean init_coord);

NcSphereMapAlm *nc_sphere_mapalm_new (void);
gboolean nc_sphere_mapalm_init (NcSphereMapAlm *mapalm, gint lmax);

#ifdef NUMCOSMO_HAVE_FFTW3
NcSphereMapSHT *nc_sphere_mapsht_new (NcSphereMap *map, NcSphereMapAlm *mapalm, guint fftw_flags);
gboolean nc_sphere_mapsht_map2alm_circle (NcSphereMapSHT *mapsht, gint ring, gint ring_size, gdouble norma, gdouble theta, gdouble phi, gint start_m, gint end_m);
gboolean nc_sphere_mapsht_alm2map_circle (NcSphereMapSHT *mapsht, gint ring, gint ring_size, gdouble theta, gdouble phi);
gboolean nc_sphere_mapsht_map2alm (NcSphereMapSHT *mapsht, gdouble cut);
gboolean nc_sphere_mapsht_alm2map (NcSphereMapSHT *mapsht);
#endif /* NUMCOSMO_HAVE_FFTW3 */
gdouble nc_sphere_map_homogenize_noise (NcSphereMap *map, gdouble base_sigma);
gdouble nc_sphere_map_rotate_avg (NcSphereMap *map, glong n);

G_END_DECLS

#endif /* _MAP_H */
