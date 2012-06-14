/***************************************************************************
 *            catalog_parser.c
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

/**
 * SECTION:catalog_parser
 * @title: Catalog Parser
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <glib.h>
#include <glib/gstdio.h>
#include <fitsio.h>

#include <gsl/gsl_histogram.h>


#define MAX_LINE_SIZE 400000
#define MAX_TEMP_STR_SIZE 4000
#define MAX_COLS 100

/**
 * nc_des_catalog_open: (skip)
 * @filename: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcDesCatalog *
nc_des_catalog_open (char *filename)
{
  NcDesCatalog *cat = (NcDesCatalog *)malloc (sizeof(NcDesCatalog));
  char line[MAX_LINE_SIZE];
  char temp_str[MAX_TEMP_STR_SIZE];
  char *col_name[MAX_COLS];
  int readed, total_readed, line_size, i;

  cat->filename = g_strdup (filename); /* comando da glib q duplica uma string */
  cat->f = fopen (cat->filename, "r");
  cat->fptr = NULL;

  if (cat->f == NULL)
  {
    fprintf (stderr, "open_catalog: catalog %s, do not exist\n", cat->filename);
    exit (0);
  }

  fgets (line, MAX_LINE_SIZE, cat->f);
  if (line[0] != '[') /* compara o valor numerico da tabela ASCII entre o lado esquerdo e direito */
  {
    fprintf (stderr, "open_catalog: catalog %s is not a valid catalog\n", cat->filename);
    exit (0);
  }

  /* Parte que le o cabecalho e guarda o numero de colunas e seus nomes */
  readed = 0;
  total_readed = 0;
  line_size = strlen (line) - 2;   /* Essa funcao retorna o tamanho da linha em bytes. O '-2' eh para subtrair os bytes correspondentes aos colchetes. */
  cat->n_cols = 0;

  while ((sscanf (&line[1 + total_readed], " '%[^']', %n", temp_str, &readed) != 0)) /* O q estah entre aspas duplas significa que serah           */
  {										     /* contado um numero qq de espacos, ignorado a primeira aspa, */
    total_readed += readed;                                                          /* atribuido a temp_str o que esta entre aspas, ignorado aspa,*/
                                                                                     /* virgula e um numero qq de espacos e atribuido a readed     */
                                                                                     /* o numero total de caracteres lidos na execucao             */
    col_name[cat->n_cols] = g_strdup(temp_str);
    cat->n_cols++;
    if (total_readed >= line_size)
      break;
  }
  cat->col_name = (char **)malloc (sizeof(char *) * cat->n_cols);
  memcpy (cat->col_name, col_name, sizeof(char *) * cat->n_cols);

  cat->col_type = (NcDesCatalogColType *)malloc(sizeof(NcDesCatalogColType) * cat->n_cols);

  /* Parte que le o numero de linhas e os tipos de cada coluna */

  fgets (line, MAX_LINE_SIZE, cat->f);
  if (line[0] != '(') /* compara o valor numerico da tabela ASCII entre o lado esquerdo e direito */
  {
    fprintf (stderr, "open_catalog: catalog %s is not a valid catalog\n", cat->filename);
    exit (0);
  }

  i = 0;
  readed = 0;
  total_readed = 0;
  line_size = strlen (line) - 2;   /* Essa funcao retorna o tamanho da linha em bytes. O '-2' eh para subtrair os bytes correspondentes aos parenteses. */
  while ((sscanf (&line[1 + total_readed], " %[^,)], %n", temp_str, &readed) != 0))
  {
    size_t temp_str_size = strlen(temp_str);
    int j = 0;
    total_readed += readed;

    cat->col_type[i] = NC_DES_COL_INT;
    while (j < temp_str_size)
      if (temp_str[j++] == '.')
      {
        cat->col_type[i] = NC_DES_COL_DOUBLE;
        break;
      }

    i++;

    if (total_readed >= line_size)
      break;
  }

  cat->n_rows = 1;
  while (fgets (line, MAX_LINE_SIZE, cat->f) != NULL)
    cat->n_rows++;

  return cat;
}

#define NC_FITS_ERROR(status) \
do { \
  if (status) \
  { \
    gchar errormsg[30]; \
    fits_get_errstatus (status, errormsg); \
    g_error ("FITS: %s", errormsg); \
  } \
} while (FALSE);

/**
 * nc_des_catalog_fits_open: (skip)
 * @filename: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcDesCatalog *
nc_des_catalog_fits_open (char *filename)
{
  NcDesCatalog *cat = (NcDesCatalog *)malloc (sizeof(NcDesCatalog));
  char comment[FLEN_COMMENT];
  int i, status = 0, hdutype;

  /*                     */
  cat->filename = g_strdup (filename);
  cat->f = NULL;

  /*                     */
  fits_open_file(&cat->fptr, filename, READONLY, &status);
  NC_FITS_ERROR(status);

  /*                     */
  fits_movabs_hdu(cat->fptr, 2, &hdutype, &status);
  NC_FITS_ERROR(status);
  if (hdutype != BINARY_TBL)
    g_error ("Binary data not found\n");

  /*                     */
  fits_read_key_lng (cat->fptr, "TFIELDS", &cat->n_cols, comment, &status);
  NC_FITS_ERROR(status);
  fits_read_key_lng (cat->fptr, "NAXIS2", &cat->n_rows, comment, &status);
  NC_FITS_ERROR(status);


  /*                     */
  cat->col_name = (char **)malloc (sizeof(char *) * cat->n_cols);
  cat->col_type = (NcDesCatalogColType *)malloc(sizeof(NcDesCatalogColType) * cat->n_cols);

  /*                     */
  for (i = 1; i <= cat->n_cols; i++)
  {
    char *ttype = g_strdup_printf ("TTYPE%d", i);
    char *tform = g_strdup_printf ("TFORM%d", i);
    char col_label[30], col_type[30];
    fits_read_key_str (cat->fptr, ttype, col_label, comment, &status);
    NC_FITS_ERROR(status);
    fits_read_key_str (cat->fptr, tform, col_type,  comment, &status);
    NC_FITS_ERROR(status);

    cat->col_name[i-1] = g_strdup (col_label);

    switch (col_type[0])
    {
      case 'J':
        cat->col_type[i-1] = NC_DES_COL_LONG;
        break;
      case 'E':
        cat->col_type[i-1] = NC_DES_COL_FLOAT;
        break;
      case 'F':
        cat->col_type[i-1] = NC_DES_COL_FLOAT;
        break;
      case 'I':
        cat->col_type[i-1] = NC_DES_COL_INT;
        break;
      case 'D':
        cat->col_type[i-1] = NC_DES_COL_DOUBLE;
        break;
    }

    free (tform);
    free (ttype);
  }

  return cat;
}

/**
 * nc_des_catalog_read_cols: (skip)
 * @cat: a #NcDesCatalog
 * @col_name: FIXME
 * @ncols: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcDesCatalogCol *
nc_des_catalog_read_cols (NcDesCatalog *cat, char **col_name, long int ncols)
{
  int i;
  NcDesCatalogCol *cols = (NcDesCatalogCol *) malloc (sizeof (NcDesCatalogCol) * ncols);
  size_t err = 0;

  memset (cols, 0, sizeof (NcDesCatalogCol) * ncols);

  if (ncols > cat->n_cols)
  {
    fprintf (stderr, "read_cols: ncols[%ld] bigger than cat->n_cols[%ld]\n", ncols, cat->n_cols);
    exit(0);
  }

  for (i = 0; i < cat->n_cols; i++)
  {
    int j;
    for (j = 0; j < ncols; j++)
      if (strcmp(cat->col_name[i], col_name[j]) == 0)
      {
        cols[j].index = i;
        cols[j].type = cat->col_type[i];
        cols[j].init = 1;
        cols[j].n_rows = cat->n_rows;
        if (cols[j].type == NC_DES_COL_INT)
          cols[j].col_data = (int *)malloc(sizeof(int) * cat->n_rows);
        else
          cols[j].col_data = (gdouble *)malloc(sizeof(gdouble) * cat->n_rows);
      }
  }

  for (i = 0; i < ncols; i++)
    if (!cols[i].init)
    {
      fprintf (stderr, "read_cols: cant find col '%s'\n", col_name[i]);
      err++;
    }
  if (err)
    exit(0);

  rewind(cat->f);
  {
    char line[MAX_LINE_SIZE];
    char temp_str[MAX_TEMP_STR_SIZE];
    int readed, total_readed, line_size, line_i;
    fgets (line, MAX_LINE_SIZE, cat->f);

    line_i = 0;
    while (fgets (line, MAX_LINE_SIZE, cat->f) != NULL)
    {
      i = 0;
      readed = 0;
      total_readed = 0;
      line_size = strlen (line) - 2;   /* Essa funcao retorna o tamanho da linha em bytes. O '-2' eh para subtrair os bytes correspondentes aos parenteses. */
      while ((sscanf (&line[1 + total_readed], " %[^,)], %n", temp_str, &readed) != 0))
      {
        int j;
        total_readed += readed;
        for (j = 0; j < ncols; j++)
        {
          if (cols[j].index == i)
          {
            if (cols[j].type == NC_DES_COL_INT)
              sscanf (temp_str, "%d", &(NC_POINTER_INT(cols[j].col_data)[line_i]));
            else if (cols[j].type == NC_DES_COL_DOUBLE)
              sscanf (temp_str, "%lg", &(NC_POINTER_DOUBLE(cols[j].col_data)[line_i]));
          }
        }
        i++;
        if (total_readed >= line_size)
          break;
      }
      line_i++;
    }
  }

  return cols;
}

/**
 * nc_des_catalog_fits_read_cols: (skip)
 * @cat: a #NcDesCatalog
 * @col_name: FIXME
 * @ncols: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcDesCatalogCol *
nc_des_catalog_fits_read_cols (NcDesCatalog *cat, char **col_name, long int ncols)
{
  int i, status = 0;
  NcDesCatalogCol *cols = (NcDesCatalogCol *) malloc (sizeof (NcDesCatalogCol) * ncols);
  size_t err = 0;

  memset (cols, 0, sizeof (NcDesCatalogCol) * ncols);

  if (ncols > cat->n_cols)
  {
    fprintf (stderr, "read_cols: ncols[%ld] bigger than cat->n_cols[%ld]\n", ncols, cat->n_cols);
    exit(0);
  }

  for (i = 0; i < cat->n_cols; i++)
  {
    int j;
    for (j = 0; j < ncols; j++)
      if (strcmp(cat->col_name[i], col_name[j]) == 0)
      {
        cols[j].index = i;
        cols[j].type = cat->col_type[i];
        cols[j].init = 1;
        cols[j].n_rows = cat->n_rows;
        switch (cols[j].type)
        {
          case NC_DES_COL_INT:
            cols[j].col_data = (int *)malloc(sizeof(int) * cat->n_rows);
            fits_read_col_int (cat->fptr, i+1, 1, 1, cat->n_rows, 0, cols[j].col_data, NULL, &status);
            NC_FITS_ERROR(status);
            break;
          case NC_DES_COL_UINT:
            cols[j].col_data = (unsigned int *)malloc(sizeof(unsigned int) * cat->n_rows);
            fits_read_col_uint (cat->fptr, i+1, 1, 1, cat->n_rows, 0, cols[j].col_data, NULL, &status);
            NC_FITS_ERROR(status);
            break;
          case NC_DES_COL_LONG:
            cols[j].col_data = (long int *)malloc(sizeof(long int) * cat->n_rows);
            fits_read_col_lng (cat->fptr, i+1, 1, 1, cat->n_rows, 0, cols[j].col_data, NULL, &status);
            NC_FITS_ERROR(status);
            break;
          case NC_DES_COL_DOUBLE:
            cols[j].col_data = (gdouble *)malloc(sizeof(gdouble) * cat->n_rows);
            fits_read_col_dbl (cat->fptr, i+1, 1, 1, cat->n_rows, 0, cols[j].col_data, NULL, &status);
            NC_FITS_ERROR(status);
            break;
          case NC_DES_COL_FLOAT:
            cols[j].col_data = (float *)malloc(sizeof(float) * cat->n_rows);
            fits_read_col_flt (cat->fptr, i+1, 1, 1, cat->n_rows, 0, cols[j].col_data, NULL, &status);
            NC_FITS_ERROR(status);
            break;
        }
      }
  }

  for (i = 0; i < ncols; i++)
    if (!cols[i].init)
    {
      fprintf (stderr, "read_cols: can't find col '%s'\n", col_name[i]);
      err++;
    }
  if (err)
    exit(0);

  return cols;
}

/**
 * nc_des_catalog_free: (skip)
 * @cat: a #NcDesCatalog
 *
 * FIXME
 *
*/
void
nc_des_catalog_free (NcDesCatalog *cat)
{
  int i, status = 0;
  free (cat->filename);

  if (cat->f != NULL)
    fclose(cat->f);
  else if (cat->fptr != NULL)
  {
    fits_close_file (cat->fptr, &status);
    NC_FITS_ERROR(status);
  }

  for (i = 0; i < cat->n_cols; i++)
    free (cat->col_name[i]);
  free (cat->col_type);

  free (cat);
}

/**
 * nc_des_catalog_col_free:
 * @cols: a #NcDesCatalogCol
 * @ncols: FIXME
 *
 * FIXME
 *
*/
void
nc_des_catalog_col_free (NcDesCatalogCol *cols, long int ncols)
{
  int i;
  for (i = 0; i < ncols; i++)
    free (cols[i].col_data);
  free(cols);
}

/**
 * nc_des_catalog_col_histogram_create_range: (skip)
 * @col: a #NcDesCatalogCol
 * @col_w: a #NcDesCatalogCol
 * @range: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gsl_histogram *
nc_des_catalog_col_histogram_create_range (NcDesCatalogCol *col, NcDesCatalogCol *col_w, gdouble *range, size_t n)
{
  gsl_histogram * h = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges (h, range, n);
  nc_des_catalog_col_histogram_fill (h, col, col_w);

  return h;
}

/**
 * nc_des_catalog_col_histogram_create_uniform: (skip)
 * @col: a #NcDesCatalogCol
 * @col_w: a #NcDesCatalogCol
 * @a: FIXME
 * @b: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gsl_histogram *
nc_des_catalog_col_histogram_create_uniform (NcDesCatalogCol *col, NcDesCatalogCol *col_w, gdouble a, gdouble b, size_t n)
{
  gsl_histogram * h = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h, a, b);
  nc_des_catalog_col_histogram_fill (h, col, col_w);

  return h;
}

/**
 * nc_des_catalog_col_histogram_fill: (skip)
 * @h: FIXME
 * @col: a #NcDesCatalogCol
 * @col_w: a #NcDesCatalogCol
 *
 * FIXME
 *
*/
void
nc_des_catalog_col_histogram_fill (gsl_histogram * h, NcDesCatalogCol *col, NcDesCatalogCol *col_w)
{
  guint i;

  if (col_w == NULL)
  {
    if (col->type == NC_DES_COL_INT)
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_increment (h, NC_POINTER_INT(col->col_data)[i]);
    else if (col->type == NC_DES_COL_DOUBLE)
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_increment (h, NC_POINTER_DOUBLE(col->col_data)[i]);
    else if (col->type == NC_DES_COL_UINT)
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_increment (h, NC_POINTER_UINT(col->col_data)[i]);
    else if (col->type == NC_DES_COL_FLOAT)
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_increment (h, NC_POINTER_FLOAT(col->col_data)[i]);
    else if (col->type == NC_DES_COL_LONG)
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_increment (h, NC_POINTER_LONG(col->col_data)[i]);
  }
  else
  {
    if ((col->type == NC_DES_COL_INT) && (col_w->type == NC_DES_COL_INT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_INT(col->col_data)[i], NC_POINTER_INT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_INT) && (col_w->type == NC_DES_COL_DOUBLE))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_INT(col->col_data)[i], NC_POINTER_DOUBLE(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_INT) && (col_w->type == NC_DES_COL_UINT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_INT(col->col_data)[i], NC_POINTER_UINT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_INT) && (col_w->type == NC_DES_COL_FLOAT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_INT(col->col_data)[i], NC_POINTER_FLOAT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_INT) && (col_w->type == NC_DES_COL_LONG))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_INT(col->col_data)[i], NC_POINTER_LONG(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_DOUBLE) && (col_w->type == NC_DES_COL_INT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_DOUBLE(col->col_data)[i], NC_POINTER_INT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_DOUBLE) && (col_w->type == NC_DES_COL_DOUBLE))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_DOUBLE(col->col_data)[i], NC_POINTER_DOUBLE(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_DOUBLE) && (col_w->type == NC_DES_COL_UINT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_DOUBLE(col->col_data)[i], NC_POINTER_UINT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_DOUBLE) && (col_w->type == NC_DES_COL_FLOAT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_DOUBLE(col->col_data)[i], NC_POINTER_FLOAT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_DOUBLE) && (col_w->type == NC_DES_COL_LONG))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_DOUBLE(col->col_data)[i], NC_POINTER_LONG(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_UINT) && (col_w->type == NC_DES_COL_INT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_UINT(col->col_data)[i], NC_POINTER_INT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_UINT) && (col_w->type == NC_DES_COL_DOUBLE))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_UINT(col->col_data)[i], NC_POINTER_DOUBLE(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_UINT) && (col_w->type == NC_DES_COL_UINT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_UINT(col->col_data)[i], NC_POINTER_UINT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_UINT) && (col_w->type == NC_DES_COL_FLOAT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_UINT(col->col_data)[i], NC_POINTER_FLOAT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_UINT) && (col_w->type == NC_DES_COL_LONG))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_UINT(col->col_data)[i], NC_POINTER_LONG(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_FLOAT) && (col_w->type == NC_DES_COL_INT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_FLOAT(col->col_data)[i], NC_POINTER_INT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_FLOAT) && (col_w->type == NC_DES_COL_DOUBLE))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_FLOAT(col->col_data)[i], NC_POINTER_DOUBLE(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_FLOAT) && (col_w->type == NC_DES_COL_UINT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_FLOAT(col->col_data)[i], NC_POINTER_UINT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_FLOAT) && (col_w->type == NC_DES_COL_FLOAT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_FLOAT(col->col_data)[i], NC_POINTER_FLOAT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_FLOAT) && (col_w->type == NC_DES_COL_LONG))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_FLOAT(col->col_data)[i], NC_POINTER_LONG(col_w->col_data)[i]);
   else if ((col->type == NC_DES_COL_LONG) && (col_w->type == NC_DES_COL_INT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_LONG(col->col_data)[i], NC_POINTER_INT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_LONG) && (col_w->type == NC_DES_COL_DOUBLE))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_LONG(col->col_data)[i], NC_POINTER_DOUBLE(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_LONG) && (col_w->type == NC_DES_COL_UINT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_LONG(col->col_data)[i], NC_POINTER_UINT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_LONG) && (col_w->type == NC_DES_COL_FLOAT))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_LONG(col->col_data)[i], NC_POINTER_FLOAT(col_w->col_data)[i]);
    else if ((col->type == NC_DES_COL_LONG) && (col_w->type == NC_DES_COL_LONG))
      for (i = 0; i < col->n_rows; i++)
        gsl_histogram_accumulate (h, NC_POINTER_LONG(col->col_data)[i], NC_POINTER_LONG(col_w->col_data)[i]);
 }
}

/**
 * nc_des_catalog_col_histogram_fill_gt: (skip)
 * @h: FIXME
 * @col: a #NcDesCatalogCol
 * @col_cmp: a #NcDesCatalogCol
 * @min: FIXME
 *
 * FIXME
 *
*/
void
nc_des_catalog_col_histogram_fill_gt (gsl_histogram * h, NcDesCatalogCol *col, NcDesCatalogCol *col_cmp, gdouble min)
{
  int i;

  if ((col->type == NC_DES_COL_INT) && (col_cmp->type == NC_DES_COL_INT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_INT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_INT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_INT) && (col_cmp->type == NC_DES_COL_DOUBLE))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_INT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_DOUBLE(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_INT) && (col_cmp->type == NC_DES_COL_UINT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_INT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_UINT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_INT) && (col_cmp->type == NC_DES_COL_FLOAT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_INT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_FLOAT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_INT) && (col_cmp->type == NC_DES_COL_LONG))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_INT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_LONG(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_DOUBLE) && (col_cmp->type == NC_DES_COL_INT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_DOUBLE(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_INT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_DOUBLE) && (col_cmp->type == NC_DES_COL_DOUBLE))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_DOUBLE(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_DOUBLE(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_DOUBLE) && (col_cmp->type == NC_DES_COL_UINT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_DOUBLE(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_UINT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_DOUBLE) && (col_cmp->type == NC_DES_COL_FLOAT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_DOUBLE(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_FLOAT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_DOUBLE) && (col_cmp->type == NC_DES_COL_LONG))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_DOUBLE(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_LONG(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_UINT) && (col_cmp->type == NC_DES_COL_INT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_UINT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_INT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_UINT) && (col_cmp->type == NC_DES_COL_DOUBLE))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_UINT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_DOUBLE(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_UINT) && (col_cmp->type == NC_DES_COL_UINT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_UINT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_UINT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_UINT) && (col_cmp->type == NC_DES_COL_FLOAT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_UINT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_FLOAT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_UINT) && (col_cmp->type == NC_DES_COL_LONG))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_UINT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_LONG(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_FLOAT) && (col_cmp->type == NC_DES_COL_INT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_FLOAT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_INT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_FLOAT) && (col_cmp->type == NC_DES_COL_DOUBLE))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_FLOAT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_DOUBLE(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_FLOAT) && (col_cmp->type == NC_DES_COL_UINT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_FLOAT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_UINT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_FLOAT) && (col_cmp->type == NC_DES_COL_FLOAT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_FLOAT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_FLOAT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_FLOAT) && (col_cmp->type == NC_DES_COL_LONG))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_FLOAT(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_LONG(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_LONG) && (col_cmp->type == NC_DES_COL_INT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_LONG(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_INT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_LONG) && (col_cmp->type == NC_DES_COL_DOUBLE))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_LONG(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_DOUBLE(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_LONG) && (col_cmp->type == NC_DES_COL_UINT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_LONG(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_UINT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_LONG) && (col_cmp->type == NC_DES_COL_FLOAT))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_LONG(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_FLOAT(col->col_data)[i]);
    }
  else if ((col->type == NC_DES_COL_LONG) && (col_cmp->type == NC_DES_COL_LONG))
    for (i = 0; i < col->n_rows; i++)
    {
      if (NC_POINTER_LONG(col_cmp->col_data)[i] > min)
        gsl_histogram_increment (h, NC_POINTER_LONG(col->col_data)[i]);
    }
}
