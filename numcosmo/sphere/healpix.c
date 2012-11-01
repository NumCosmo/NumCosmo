/***************************************************************************
 *            healpix.c
 *
 *  Wed Jul  9 11:09:37 2008
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
 * SECTION:healpix
 * @title: Healpix
 * @short_description: Healpix re-implementation
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "sphere/healpix.h"

#include <glib/gstdio.h>
#include <gsl/gsl_math.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */

/***************************************************************************
 *
 *
 ****************************************************************************/

#define NC_FITS_ERROR(status) \
do { \
  gchar errormsg[30]; \
  fits_get_errstatus (status, errormsg); \
  g_error ("FITS: %s", errormsg); \
} while (FALSE);

/**
 * nc_sphere_healpix_str2type: (skip)
 * @str_type: FIXME
 * 
 * FIXME
 *
 * Returns: FIXME
 */
NcSphereMapType
nc_sphere_healpix_str2type (gchar *str_type)
{
  switch (*str_type)
  {
    case 'T':
      return NC_SPHERE_MAP_TYPE_TEMPERATURE; 
      break;
    case 'Q':
      return NC_SPHERE_MAP_TYPE_Q_POLARIZATION; 
      break;
    case 'U':
      return NC_SPHERE_MAP_TYPE_U_POLARISATION; 
      break;
    case 'S':
      return NC_SPHERE_MAP_TYPE_SPUR_SIGNAL; 
      break;
    case 'N':
      return NC_SPHERE_MAP_TYPE_N_OBS; 
      break;
    default:
      g_error ("Unknown type (%s)", str_type);
  }
  return 0;
}

#ifdef NUMCOSMO_HAVE_CFITSIO
/**
 * nc_sphere_healpix_read_map: (skip)
 * @fits_file: FIXME
 * @map: a #NcSphereMap
 * 
 * FIXME
 *
 * Returns: FIXME
 */
NcSphereMap *
nc_sphere_healpix_read_map (gchar *fits_file, NcSphereMap *map) 
{
  glong nside, naxis;
  gint  status, hdutype, anynul, i;   
  gchar comment[FLEN_COMMENT];
  gchar ordering[30];
  gchar col_label[30];    
  fitsfile *fptr;

  status = 0;

  if ( fits_open_file(&fptr, fits_file, READONLY, &status) ) 
    NC_FITS_ERROR(status);
   
  if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) 
    NC_FITS_ERROR(status);

  if (hdutype != BINARY_TBL) 
    g_error ("%s (%d): Extension is not binary!\n", __FILE__, __LINE__);
   
  if ( fits_read_key_lng (fptr, "NSIDE", &nside, comment, &status) ) 
    NC_FITS_ERROR(status);
  
  if (map == NULL)
    map = nc_sphere_map_new (nside);
  else
    g_assert (map->nside == nside);
  
  if ( fits_read_key_lng (fptr, "TFIELDS", &naxis, comment, &status) ) 
    NC_FITS_ERROR(status);
  
  for (i = 0; i < naxis; i++)
  {
    gchar *ttype = g_strdup_printf ("TTYPE%d", i + 1);
    NcSphereMapType type;
    if (fits_read_key_str (fptr, ttype, col_label, comment, &status))
      NC_FITS_ERROR(status);
    type = nc_sphere_healpix_str2type (col_label);
    map->type |= type;
    switch (type)
    {
      case NC_SPHERE_MAP_TYPE_TEMPERATURE:
        map->dt = (map->dt == NULL) ? gsl_vector_float_alloc (map->npix) : map->dt;
        if ( fits_read_col_flt (fptr, i + 1, 1, 1, map->npix, NC_HEALPIX_NULLVAL, map->dt->data, &anynul, &status) ) 
          NC_FITS_ERROR(status);
        break;
      case NC_SPHERE_MAP_TYPE_Q_POLARIZATION:
        map->qpol = (map->qpol == NULL) ? gsl_vector_float_alloc (map->npix) : map->qpol;
        if ( fits_read_col_flt (fptr, i + 1, 1, 1, map->npix, NC_HEALPIX_NULLVAL, map->qpol->data, &anynul, &status) ) 
          NC_FITS_ERROR(status);
        break;
      case NC_SPHERE_MAP_TYPE_U_POLARISATION:
        map->upol = (map->upol == NULL) ? gsl_vector_float_alloc (map->npix) : map->upol;
        if ( fits_read_col_flt (fptr, i + 1, 1, 1, map->npix, NC_HEALPIX_NULLVAL, map->upol->data, &anynul, &status) ) 
          NC_FITS_ERROR(status);
        break;
      case NC_SPHERE_MAP_TYPE_SPUR_SIGNAL:
        map->spur_signal = (map->spur_signal == NULL) ? gsl_vector_float_alloc (map->npix) : map->spur_signal;
        if ( fits_read_col_flt (fptr, i + 1, 1, 1, map->npix, NC_HEALPIX_NULLVAL, map->spur_signal->data, &anynul, &status) ) 
          NC_FITS_ERROR(status);
        break;
      case NC_SPHERE_MAP_TYPE_N_OBS:
        map->nobs = (map->nobs == NULL) ? gsl_vector_float_alloc (map->npix) : map->nobs;
        if ( fits_read_col_flt (fptr, i + 1, 1, 1, map->npix, NC_HEALPIX_NULLVAL, map->nobs->data, &anynul, &status) ) 
          NC_FITS_ERROR(status);
        break;
    }
    g_free (ttype);
  }

  if (fits_read_key(fptr, TSTRING, "ORDERING", ordering, comment, &status)) 
  {
    g_warning ("Could not find ORDERING in the fits file, assuming RING");
    status = 0;
  }
     
  if ( fits_close_file(fptr, &status) ) 
    NC_FITS_ERROR(status);
   
  map->nrings = map->nside * 4 - 1;

  switch (*ordering)
  {
    case 'N':
      map->order = NC_SPHERE_MAP_ORDER_NEST;
      break;
    case 'R':
      map->order = NC_SPHERE_MAP_ORDER_RING;
      break;
    default:
      g_error ("Unknown order type (%s)", ordering);
  }
  map->loaded = TRUE;
  
  /* Later */
  return map;
}

/**
 * nc_sphere_healpix_write_map: 
 * @map: a #NcSphereMap
 * @filename: FIXME 
 * @overwrite: FIXME 
 * 
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_sphere_healpix_write_map (NcSphereMap *map, gchar *filename, gboolean overwrite)
{
  /*******************************************************************/
  /* Create a binary table extension                                 */
  /*******************************************************************/
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  gint status, hdutype;
  glong firstrow, firstelem;
  
  gint bitpix   =  SHORT_IMG;
  glong naxis   =   0;
  glong naxes[] = {0,0};
  
  gint tfields   = 1;
  
  gchar *order;                 /* HEALPix ordering */
  gchar extname[] = "BINTABLE";   /* extension name */  
  gchar *ttype[] = { "SIGNAL" };
  gchar *tform[] = { "1E" };
  gchar *tunit[] = { " " }; 
  
  /* initialize status before calling fitsio routines */
  status = 0;
  
  if (overwrite && g_file_test (filename, G_FILE_TEST_EXISTS))
    g_unlink (filename);
  
  /* create new FITS file */
  if (fits_create_file (&fptr, filename, &status))
    NC_FITS_ERROR(status);

  if (fits_create_img (fptr,  bitpix, naxis, naxes, &status) )
    NC_FITS_ERROR(status);

  if (fits_write_date (fptr, &status) )
    NC_FITS_ERROR(status);
  
  /* move to 1nd HDU  */
  if ( fits_movabs_hdu(fptr, 1, &hdutype, &status) )
    NC_FITS_ERROR(status);
  
  /* append a new empty binary table onto the FITS file */
  if ( fits_create_tbl( fptr, BINARY_TBL, map->npix, tfields, ttype, tform,
                        tunit, extname, &status) )
    NC_FITS_ERROR(status);  
  
  if (fits_write_key(fptr, TSTRING, "PIXTYPE", "HEALPIX", "HEALPIX Pixelisation", &status))   
    NC_FITS_ERROR(status);  
  
  order = map->order == NC_SPHERE_MAP_ORDER_NEST ? "NESTED  " : "RING    ";

  if (fits_write_key( fptr, TSTRING, "ORDERING", order,
                      "Pixel ordering scheme, either RING or NESTED", &status))
    NC_FITS_ERROR(status);
  
  if (fits_write_key (fptr, TLONG, "NSIDE", &map->nside, "Resolution parameter for HEALPIX", &status))
    NC_FITS_ERROR(status);
  
  if (fits_write_key (fptr, TSTRING, "COORDSYS", "G", "Pixelisation coordinate system", &status))
    NC_FITS_ERROR(status);

  if (fits_write_comment(fptr,"           G = Galactic, E = ecliptic, C = celestial = equatorial  ", &status))
    NC_FITS_ERROR(status);

  firstrow  = 1;  /* first row in table to write   */
  firstelem = 1;  /* first element in row  (ignored in ASCII tables)  */

  if (fits_write_col (fptr, TFLOAT, 1, firstrow, firstelem, map->npix, map->dt->data, &status))
    NC_FITS_ERROR(status);

  if ( fits_close_file(fptr, &status) )      
    NC_FITS_ERROR(status);
  
  return TRUE;
}

#endif /* NUMCOSMO_HAVE_CFITSIO */

/**
 * nc_sphere_healpix_nest2ring: 
 * @nside: FIXME
 * @nest_index: FIXME
 * 
 * FIXME
 *
 * Returns: FIXME
 */
glong
nc_sphere_healpix_nest2ring (gint nside, glong nest_index)
{
  glong ring_index;
  glong npix = 12 * gsl_pow_2 (nside);
  glong base_pixel;
  gint middle_rings_size = 4 * nside;  
  glong face_size = nside * nside;          /* Face size in pixels                                     */
  glong cap_size = 2 * nside * (nside - 1); /* Cap size in pixels                                      */
  gint f, h, v;                             /* Face number, horizontal coordinate, vertical coordinate */
  gint t, p, s, pad, w, l;                  /* theta, phi, shift, padding, width, local index          */
  gint x, y;
  gint hf;

  g_assert (nest_index < npix);
  
  f = nest_index / face_size;
  l = nest_index % face_size;
  HEALPIX_INT_TO_XY (l, x, y);

  h = nside - 1 - y;
  v = 2 * nside - 2 - y - x;
  
  switch (f / 4)
  {
    case 0:
      t = v;
      if (t < (nside - 1))
      {
        w   = (t + 1);
        s   = 0;
        pad = 0;
        hf  = f % 4;
        base_pixel = 2 * (gsl_pow_2 (t + 1) - t - 1);
      }
      else
      {
        w   = nside;
        s   = (t - (nside - 1)) / 2;
        pad = 0;
        hf  = f % 4;
        base_pixel = cap_size + (t - nside + 1) * middle_rings_size;
      }
      break;
    case 1:
      t   = v + nside;
      w   = nside;
      s   = (t - (nside - 1)) / 2;
      pad = 0;
      hf  = f % 4;
      base_pixel = cap_size + (v + 1) * middle_rings_size;
      break;
    case 2:
      t   = v + 2 * nside;
      if (v < nside)
      {
        w   = nside;
        s   = (t - (nside - 1)) / 2;
        pad = 0;
        hf  = f % 4 + 1;
        base_pixel = cap_size + (v + 1 + nside) * middle_rings_size;
      }
      else
      {
        w   = 4 * nside - t - 1;
        s   = w;
        pad = (t - 3 * nside + 1);
        hf  = f % 4 + 1;
        base_pixel = (npix - cap_size) + (4 * nside * (v - nside + 1) - 2 * gsl_pow_2 (v - nside + 1) + 2 * (v - nside) + 2 - 4 * nside);
      }
      break;
    default:
      g_assert_not_reached ();
  }

  p = (h + hf * w - s - pad);
  if (p < 0) p += 4 * w;

  ring_index = base_pixel + p;
  return ring_index;
}

/**
 * nc_sphere_healpix_ring2nest: 
 * @nside: FIXME
 * @ring_index: FIXME
 * 
 * FIXME
 *
 * Returns: FIXME
 */
glong
nc_sphere_healpix_ring2nest (gint nside, glong ring_index)
{
  glong nest_index;
  glong npix = 12 * gsl_pow_2 (nside);
  gint middle_rings_size = 4 * nside;
  glong cap_size = 2 * nside * (nside - 1); /* Cap size in pixels                                      */
  gint f, h, v;                             /* Face number, horizontal coordinate, vertical coordinate */
  gint t, p, s, pad, w, l;                  /* theta, phi, shift, padding, width, local index          */
  gint x, y;
  
  g_assert (ring_index < npix);

  if (ring_index < cap_size)
  {
    t   = (sqrt (1 + 2 * ring_index) - 1) / 2;
    w   = (t + 1);
    p   = ring_index - 2 * (gsl_pow_2 (t + 1) - t - 1);
    s   = 0;
    pad = 0;
  }
  else if (ring_index < (npix - cap_size))
  {
    l   = ring_index - cap_size;
    w   = nside;
    t   = (gint)(l / middle_rings_size) + (nside - 1);
    p   = l % middle_rings_size;
    s   = (t - (nside - 1)) / 2;
    pad = 0;
  }
  else
  {
    l   = ring_index - npix + cap_size;
    t   = 4 * nside - (1 + sqrt (4 * gsl_pow_2 (nside) - 4 * nside + 1 - 2 * l)) / 2;
    w   = 4 * nside - t - 1;
    p   = l - (4 * nside * (t - 3 * nside + 1) - 2 * gsl_pow_2 (t - 3 * nside + 1) + 2 * (t - 3 * nside) + 2 - 4 * nside);
    s   = w;
    pad = (t - 3 * nside + 1);
  }

  h  = (p + s) % w + pad;
  v  = (t - h) % nside + h;
  f  = ((gint)((t - h) / nside)) * 4;
  f += (f == 8) ? ((gint)((p + s) / w) - 1) : ((gint)((p + s) / w) % 4);
  
  x = nside - 1 - (v - h);
  y = nside - 1 - h;
  
  HEALPIX_XY_TO_INT (x, y, nest_index);
  
  nest_index += gsl_pow_2 (nside) * f;
  
  return nest_index;
}

static void
_t_p_w_to_theta_phi (const gint nside, const gint tm1, const gint pm1, const gint w, gdouble *theta, gdouble *phi)
{
  const gint t = tm1 + 1;
  const gint p = pm1 + 1;
  if (t < nside)
  {
    *theta = acos((1.0 - t * t * 1.0 / (3.0 * nside * nside)));
    *phi = (p - 0.5) * M_PI_2 / (1.0 * w);
  }
  else if (t > 3 * nside)
  {
    gint tt = 4 * nside - t;
    *theta = acos(-(1.0 - tt * tt * 1.0 / (3.0 * nside * nside)) );
    *phi = (p - 0.5) * M_PI_2 / (1.0 * w);
  }
  else
  {
    *theta = acos((2.0 * nside - t) * 2.0 / (3.0 * nside));
    *phi = (p - ((t - nside) % 2 + 1.0) * 0.5) * M_PI_2 / (1.0 * w);
  }
}

static void
_t_p_w_to_vector (const gint nside, const gint tm1, const gint pm1, const gint w, NcTriVector v)
{
  const gint t = tm1 + 1;
  const gint p = pm1 + 1;
  gdouble phi, z, sin_theta;
  if (t < nside)
  {
    phi = (p - 0.5) * M_PI_2 / (1.0 * w);
    z = (1.0 - t * t * 1.0 / (3.0 * nside * nside));
  }
  else if (t > 3 * nside)
  {
    const gint tt = 4 * nside - t;
    phi = (p - 0.5) * M_PI_2 / (1.0 * w);
    z = -(1.0 - tt * tt * 1.0 / (3.0 * nside * nside));
  }
  else
  {
    phi = (p - ((t - nside) % 2 + 1.0) * 0.5) * M_PI_2 / (1.0 * w);
    z = (2.0 * nside - t) * 2.0 / (3.0 * nside);
  }
  sin_theta = sqrt (1.0 - z * z);
  v.c[0] = sin_theta * cos (phi);
  v.c[1] = sin_theta * sin (phi);
  v.c[2] = z;
}

/**
 * nc_sphere_healpix_pix2ang_nest:
 * @nside: FIXME
 * @nest_index: FIXME
 * @theta: FIXME
 * @phi: FIXME
 *
 * FIXME 
*/ 
void 
nc_sphere_healpix_pix2ang_nest (gint nside, glong nest_index, gdouble *theta, gdouble *phi)
{
  glong npix = 12 * gsl_pow_2 (nside);
  glong face_size = nside * nside;          /* Face size in pixels                                     */
  gint f, h, v;                             /* Face number, horizontal coordinate, vertical coordinate */
  gint t, p, s, pad, w, l;                  /* theta, phi, shift, padding, width, local index          */
  gint x, y;
  gint hf;

  g_assert (nest_index < npix);
  
  f = nest_index / face_size;
  l = nest_index % face_size;
  HEALPIX_INT_TO_XY (l, x, y);
 
  h = nside - 1 - y;
  v = 2 * nside - 2 - y - x;
  
  switch (f / 4)
  {
    case 0:
      t = v;
      if (t < (nside - 1))
      {
        w   = (t + 1);
        s   = 0;
        pad = 0;
        hf  = f % 4;
      }
      else
      {
        w   = nside;
        s   = (t - (nside - 1)) / 2;
        pad = 0;
        hf  = f % 4;
      }
      break;
    case 1:
      t   = v + nside;
      w   = nside;
      s   = (t - (nside - 1)) / 2;
      pad = 0;
      hf  = f % 4;
      break;
    case 2:
      t   = v + 2 * nside;
      if (v < nside)
      {
        w   = nside;
        s   = (t - (nside - 1)) / 2;
        pad = 0;
        hf  = f % 4 + 1;
      }
      else
      {
        w   = 4 * nside - t - 1;
        s   = w;
        pad = (t - 3 * nside + 1);
        hf  = f % 4 + 1;
      }
      break;
    default:
      g_assert_not_reached ();
  }

  p = (h + hf * w - s - pad);
  if (p < 0) p += 4 * w;
  _t_p_w_to_theta_phi (nside, t, p, w, theta, phi);
}

/**
 * nc_sphere_healpix_pix2ang_ring:
 * @nside: FIXME
 * @ring_index: FIXME
 * @theta: FIXME
 * @phi: FIXME
 *
 * FIXME 
*/ 
void
nc_sphere_healpix_pix2ang_ring (gint nside, glong ring_index, gdouble *theta, gdouble *phi)
{
  glong npix = 12 * gsl_pow_2 (nside);
  gint middle_rings_size = 4 * nside;
  glong cap_size = 2 * nside * (nside - 1); /* Cap size in pixels                                      */
  gint t, p, w, l;                  /* theta, phi, shift, padding, width, local index          */
  
  g_assert (ring_index < npix);

  if (ring_index < cap_size)
  {
    t   = (sqrt (1 + 2 * ring_index) - 1) / 2;
    w   = (t + 1);
    p   = ring_index - 2 * (gsl_pow_2 (t + 1) - t - 1);
  }
  else if (ring_index < (npix - cap_size))
  {
    l   = ring_index - cap_size;
    w   = nside;
    t   = (gint)(l / middle_rings_size) + (nside - 1);
    p   = l % middle_rings_size;
  }
  else
  {
    l   = ring_index - npix + cap_size;
    t   = 4 * nside - (1 + sqrt (4 * gsl_pow_2 (nside) - 4 * nside + 1 - 2 * l)) / 2;
    w   = 4 * nside - t - 1;
    p   = l - (4 * nside * (t - 3 * nside + 1) - 2 * gsl_pow_2 (t - 3 * nside + 1) + 2 * (t - 3 * nside) + 2 - 4 * nside);
  }
  _t_p_w_to_theta_phi (nside, t, p, w, theta, phi);
}

/**
 * nc_sphere_healpix_pix2vec_ring:
 * @nside: FIXME
 * @ring_index: FIXME
 * @v: a #NcTriVector
 *
 * FIXME 
*/
void 
nc_sphere_healpix_pix2vec_ring (gint nside, glong ring_index, NcTriVector v)
{
  glong npix = 12 * gsl_pow_2 (nside);
  gint middle_rings_size = 4 * nside;
  glong cap_size = 2 * nside * (nside - 1); /* Cap size in pixels                                      */
  gint t, p, w, l;                  /* theta, phi, shift, padding, width, local index          */
  
  g_assert (ring_index < npix);

  if (ring_index < cap_size)
  {
    t   = (sqrt (1 + 2 * ring_index) - 1) / 2;
    w   = (t + 1);
    p   = ring_index - 2 * (gsl_pow_2 (t + 1) - t - 1);
  }
  else if (ring_index < (npix - cap_size))
  {
    l   = ring_index - cap_size;
    w   = nside;
    t   = (gint)(l / middle_rings_size) + (nside - 1);
    p   = l % middle_rings_size;
  }
  else
  {
    l   = ring_index - npix + cap_size;
    t   = 4 * nside - (1 + sqrt (4 * gsl_pow_2 (nside) - 4 * nside + 1 - 2 * l)) / 2;
    w   = 4 * nside - t - 1;
    p   = l - (4 * nside * (t - 3 * nside + 1) - 2 * gsl_pow_2 (t - 3 * nside + 1) + 2 * (t - 3 * nside) + 2 - 4 * nside);
  }
  _t_p_w_to_vector (nside, t, p, w, v);
}

/**
 * nc_sphere_healpix_pix2vec_nest:
 * @nside: FIXME
 * @nest_index: FIXME
 * @vec: a #NcTriVector
 *
 * FIXME 
*/
void 
nc_sphere_healpix_pix2vec_nest (gint nside, glong nest_index, NcTriVector vec)
{
  glong npix = 12 * gsl_pow_2 (nside);
  glong face_size = nside * nside;          /* Face size in pixels                                     */
  gint f, h, v;                             /* Face number, horizontal coordinate, vertical coordinate */
  gint t, p, s, pad, w, l;                  /* theta, phi, shift, padding, width, local index          */
  gint x, y;
  gint hf;

  g_assert (nest_index < npix);
  
  f = nest_index / face_size;
  l = nest_index % face_size;
  HEALPIX_INT_TO_XY (l, x, y);
 
  h = nside - 1 - y;
  v = 2 * nside - 2 - y - x;
  
  switch (f / 4)
  {
    case 0:
      t = v;
      if (t < (nside - 1))
      {
        w   = (t + 1);
        s   = 0;
        pad = 0;
        hf  = f % 4;
      }
      else
      {
        w   = nside;
        s   = (t - (nside - 1)) / 2;
        pad = 0;
        hf  = f % 4;
      }
      break;
    case 1:
      t   = v + nside;
      w   = nside;
      s   = (t - (nside - 1)) / 2;
      pad = 0;
      hf  = f % 4;
      break;
    case 2:
      t   = v + 2 * nside;
      if (v < nside)
      {
        w   = nside;
        s   = (t - (nside - 1)) / 2;
        pad = 0;
        hf  = f % 4 + 1;
      }
      else
      {
        w   = 4 * nside - t - 1;
        s   = w;
        pad = (t - 3 * nside + 1);
        hf  = f % 4 + 1;
      }
      break;
    default:
      g_assert_not_reached ();
  }

  p = (h + hf * w - s - pad);
  if (p < 0) p += 4 * w;
  _t_p_w_to_vector (nside, t, p, w, vec);
}

/**
 * nc_sphere_healpix_vec2pix_ring:
 * @nside: FIXME
 * @v: a #NcTriVector
 * @i: FIXME 
 *
 * FIXME 
*/
void 
nc_sphere_healpix_vec2pix_ring (gint nside, NcTriVector v, glong *i)
{
  g_assert_not_reached ();
}
