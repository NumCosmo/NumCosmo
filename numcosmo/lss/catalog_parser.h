/***************************************************************************
 *            mass_distribution.h
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_DES_CATALOG_PARSER_H
#define _NC_DES_CATALOG_PARSER_H

#include <gsl/gsl_histogram.h>
#include <fitsio.h>

G_BEGIN_DECLS

/**
 * NcDesCatalogColType:	 
 * @NC_DES_COL_UINT: FIXME
 * @NC_DES_COL_INT: FIXME
 * @NC_DES_COL_FLOAT: FIXME
 * @NC_DES_COL_DOUBLE: FIXME
 * @NC_DES_COL_LONG: FIXME
 *
 * FIXME
 */
typedef enum _NcDesCatalogColType
{
  NC_DES_COL_UINT,
  NC_DES_COL_INT,
  NC_DES_COL_FLOAT,
  NC_DES_COL_DOUBLE,
  NC_DES_COL_LONG,
} NcDesCatalogColType;

#define NC_POINTER_UINT(a)   ((unsigned int *)(a))
#define NC_POINTER_INT(a)    ((int *)(a))
#define NC_POINTER_FLOAT(a)  ((float *)(a))
#define NC_POINTER_DOUBLE(a) ((gdouble *)(a))
#define NC_POINTER_LONG(a)   ((long int *)(a))

typedef struct _NcDesCatalog NcDesCatalog;

struct _NcDesCatalog
{
  /*< private >*/
  gchar *filename;
  fitsfile *fptr;
  FILE *f;
  glong n_cols;
  glong n_rows;
  NcDesCatalogColType *col_type; 
  gchar **col_name;
};

typedef struct _NcDesCatalogCol NcDesCatalogCol;

struct _NcDesCatalogCol
{
  /*< private >*/
  NcDesCatalogColType type;
  size_t index;
  size_t n_rows;
  gpointer col_data;
  int init;
};

NcDesCatalog *nc_des_catalog_open (gchar *filename);
NcDesCatalog *nc_des_catalog_fits_open (gchar *filename);
NcDesCatalogCol *nc_des_catalog_read_cols (NcDesCatalog *cat, gchar **col_name, glong ncols);
NcDesCatalogCol *nc_des_catalog_fits_read_cols (NcDesCatalog *cat, gchar **col_name, glong ncols);
void nc_des_catalog_free (NcDesCatalog *cat);
void nc_des_catalog_col_free (NcDesCatalogCol *cols, glong ncols);

gsl_histogram *nc_des_catalog_col_histogram_create_range (NcDesCatalogCol *col, NcDesCatalogCol *col_w, gdouble *range, size_t n);
gsl_histogram *nc_des_catalog_col_histogram_create_uniform (NcDesCatalogCol *col, NcDesCatalogCol *col_w, gdouble a, gdouble b, size_t n);
void nc_des_catalog_col_histogram_fill (gsl_histogram * h, NcDesCatalogCol *col, NcDesCatalogCol *col_w);
void nc_des_catalog_col_histogram_fill_gt (gsl_histogram * h, NcDesCatalogCol *col, NcDesCatalogCol *col_cmp, gdouble min);

G_END_DECLS

#endif /* _NC_DES_CATALOG_PARSER_H */
