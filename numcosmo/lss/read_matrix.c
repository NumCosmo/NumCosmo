/***************************************************************************
 *            read_matrix.h
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
 * SECTION:read_matrix
 * @title: Read Matrix
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
#ifdef HAVE_LIBCFITSIO
#include <fitsio.h>
#endif

#define MAX_LINE_SIZE 400000
#define MAX_TEMP_STR_SIZE 4000
#define MAX_COLS 100

/**
 * nc_arxive_open: (skip)
 * @filename: FIXME
 * FIXME
 *
 * Returns: FIXME 
*/
gsl_matrix *
nc_arxive_open (char *filename)
{
  FILE *f;
  char line[MAX_LINE_SIZE];
  gdouble temp_dbl;
  int readed, total_readed, line_size;
  gsl_matrix *m;
  size_t nrows,ncols;
  
  f = fopen (filename, "r");
  if (f == NULL)
  {
    fprintf (stderr, "open_arxive: arxive %s, do not exist\n", filename);
    exit (0);
  }
   
  if (fgets (line, MAX_LINE_SIZE, f) == NULL)
    g_error ("nc_arxive_open: io error");
   /* Parte que le o cabecalho e guarda o numero de colunas e seus nomes */
  readed = 0;
  total_readed = 0;
  line_size = strlen (line);   /* Essa funcao retorna o tamanho da linha em bytes. */
  ncols = 0;
  while ((sscanf (&line[total_readed], " %lg %n", &temp_dbl, &readed) != 0))   /* Ignora o numero de espacos em branco e conta o numero de bytes */
  {                                                                                /* do primeiro gdouble (primeiro passo do loop) e guarda em readed. */
    total_readed += readed;                                        
    ncols++;
    if (total_readed >= line_size)
      break;
  }  

  /* Parte que le o numero de linhas e os tipos de cada coluna */
  nrows = 1;
  while (fgets (line, MAX_LINE_SIZE, f) != NULL)
    nrows++;

  rewind (f);

  m = gsl_matrix_alloc (nrows, ncols);
  gsl_matrix_fscanf(f, m);


  return m;
}
