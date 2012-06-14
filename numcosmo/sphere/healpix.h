/***************************************************************************
 *            healpix.h
 *
 *  Wed Jul  9 11:11:38 2008
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
 
#ifndef _NC_HEALPIX_H
#define _NC_HEALPIX_H

#include <glib.h>

G_BEGIN_DECLS

#ifdef NUMCOSMO_HAVE_CFITSIO
NcSphereMap *nc_sphere_healpix_read_map (gchar *fits_file, NcSphereMap *map);
gboolean nc_sphere_healpix_write_map (NcSphereMap *map, gchar *filename, gboolean overwrite);
#endif /* NUMCOSMO_HAVE_CFITSIO */
glong nc_sphere_healpix_nest2ring (gint nside, glong nest_index);
glong nc_sphere_healpix_ring2nest (gint nside, glong ring_index);
void nc_sphere_healpix_pix2ang_nest (gint nside, glong nest_index, gdouble *theta, gdouble *phi);
void nc_sphere_healpix_pix2ang_ring (gint nside, glong ring_index, gdouble *theta, gdouble *phi);
void nc_sphere_healpix_pix2vec_nest (gint nside, glong nest_index, NcTriVector vec);
void nc_sphere_healpix_pix2vec_ring (gint nside, glong ring_index, NcTriVector v);
void nc_sphere_healpix_vec2pix_ring (gint nside, NcTriVector v, glong *i);

#define HEALPIX_NPIX(nside) (12*(nside)*(nside))
#define HEALPIX_INT_TO_XY(i,x,y) \
do { \
  gint shift = 0, shifted = i; \
  x = y = 0; \
  do { \
    x |= ((shifted & 1) << shift); \
    y |= (((shifted & 2) >> 1) << shift); \
    shift++; \
  } while (shifted >>= 2); \
} while(FALSE)

#define HEALPIX_XY_TO_INT(x,y,i) \
do { \
  gint shift = 0, shifted_x = x, shifted_y = y; \
  g_assert (shifted_x >= 0 && shifted_y >= 0); \
  i = 0; \
  do { \
    i |= ((shifted_x & 1) << (shift + 0)); \
    i |= ((shifted_y & 1) << (shift + 1)); \
    shift += 2; shifted_x >>= 1; shifted_y >>= 1; \
  } while (shifted_x || shifted_y); \
} while(FALSE)

#ifndef NC_HEALPIX_NULLVAL
#define NC_HEALPIX_NULLVAL (-1.6375e30) /* check if its ok to copy it here FIXME*/
#endif /* NC_HEALPIX_NULLVAL */

G_END_DECLS

#endif /* _HEALPIX_H */
