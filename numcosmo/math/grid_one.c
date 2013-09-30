/***************************************************************************
 *            grid_one.c
 *
 *  Mon Feb 22 15:05:25 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
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
 * SECTION:grid_one
 * @title: Unidimensional Grid
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/grid_one.h"
#include "math/ncm_util.h"
#include "math/ncm_cfg.h"

/**
 * ncm_grid_new: (skip)
 * @nnodes: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmGrid *
ncm_grid_new (gulong nnodes)
{
  gulong i;
  NcmGrid *grid = g_slice_new (NcmGrid);

  grid->nodes = g_slice_alloc (nnodes * sizeof(mpq_t));
  grid->nnodes = nnodes;

  grid->sections = g_array_new (TRUE, TRUE, sizeof(NcmGridSection));

  for (i = 0; i < nnodes; i++)
    mpq_init (grid->nodes[i]);
  grid->data = NULL;

  return grid;
}

/**
 * ncm_grid_new_from_sections: (skip)
 * @secs: a #NcmGridSection
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmGrid *
ncm_grid_new_from_sections (NcmGridSection *secs)
{
  NcmGrid *grid;
  NcmGridSection *secs_i = secs;
  g_assert (secs != NULL);

  while (secs_i[1].incl != 0)
    secs_i = &secs_i[1];

  grid = ncm_grid_new (secs_i[0].end);
  ncm_grid_set_sections (grid, secs);

  return grid;
}

/**
 * ncm_grid_free:
 * @grid: a #NcmGrid
 * @free_data: FIXME
 *
 * FIXME
 *
*/
void
ncm_grid_free (NcmGrid *grid, gboolean free_data)
{
  glong i;

  for (i = 0; i < grid->nnodes; i++)
    mpq_clear (grid->nodes[i]);
  g_slice_free1 (grid->nnodes * sizeof(mpq_t), grid->nodes);
  if (free_data)
  {
    g_assert (grid->data != NULL);
    g_slice_free1 (sizeof (gdouble) * grid->nnodes, grid->data);
  }

  g_slice_free (NcmGrid, grid);
}

/**
 * ncm_grid_get_name:
 * @secs: a #NcmGridSection
 *
 * FIXME
 *
 * Returns: FIXME
*/
gchar *
ncm_grid_get_name (NcmGridSection *secs)
{
  gchar *name = g_strdup ("NcGrid ");
  g_assert (secs != NULL);

  do {
    gchar *sec_name = g_strdup_printf ("%u %u %.15e %.15e ", secs[0].start, secs[0].end, secs[0].start_val, secs[0].end_val);
    gchar *temp = g_strconcat (name, sec_name, NULL);
    g_free (name);
    g_free (sec_name);
    name = temp;
  } while ((secs = &secs[1])->incl != 0);

  return name;
}

/**
 * ncm_grid_set_sections:
 * @grid: a #NcmGrid
 * @secs: a #NcmGridSection
 *
 * FIXME
 *
*/
void
ncm_grid_set_sections (NcmGrid *grid, NcmGridSection *secs)
{
  g_assert (secs != NULL);
  do {
    ncm_grid_set_nodes_d (grid, secs[0].incl,
      secs[0].start, secs[0].end,
      secs[0].start_val, secs[0].end_val
      );
  } while ((secs = &secs[1])->incl != 0);
}

/**
 * ncm_grid_set_nodes_d:
 * @grid: a #NcmGrid
 * @incl: a #NcmGridNodesEndPoints
 * @start: FIXME
 * @end: FIXME
 * @start_val: FIXME
 * @end_val: FIXME
 *
 * FIXME
 *
*/
void
ncm_grid_set_nodes_d (NcmGrid *grid, NcmGridNodesEndPoints incl, guint32 start, guint32 end, gdouble start_val, gdouble end_val)
{
  gdouble node;
  gdouble size = end - start;
  glong i;
  glong a = 0, b = 0;
  NcmGridSection sec;

  g_assert ((end > start ) && (end <= grid->nnodes));
  g_assert (end_val > start_val);

  if (grid->sections->len > 0)
  {
    gulong lasti = grid->sections->len - 1;
    g_assert (start == g_array_index(grid->sections,NcmGridSection,lasti).end);
    g_assert (start_val >= g_array_index(grid->sections,NcmGridSection,lasti).end_val);
  }
  else
    g_assert (start == 0);

  sec.incl = incl;
  sec.start = start;
  sec.end = end;
  sec.start_val = start_val;
  sec.end_val = end_val;
  g_array_append_val(grid->sections, sec);

  switch (incl)
  {
    case NCM_GRID_NODES_START:
      a = 0; b =  0; break;
    case NCM_GRID_NODES_END:
      a = 1; b =  0; break;
    case NCM_GRID_NODES_BOTH:
      a = 0; b = -1; break;
    case NCM_GRID_NODES_NONE:
      a = 1; b = +1; break;
  }

  for (i = 0; i < size; i++)
  {
    node = start_val + (end_val - start_val) / (size + b) * (i + a);
    ncm_rational_coarce_double (node, grid->nodes[start + i]);
  }
}

/**
 * ncm_grid_set_nodes_si:
 * @grid: a #NcmGrid
 * @incl: a #NcmGridNodesEndPoints
 * @start: FIXME
 * @end: FIXME
 * @start_val_num: FIXME
 * @start_val_den: FIXME
 * @end_val_num: FIXME
 * @end_val_den: FIXME
 *
 * FIXME
 *
*/
void
ncm_grid_set_nodes_si (NcmGrid *grid, NcmGridNodesEndPoints incl, guint32 start, guint32 end, glong start_val_num, glong start_val_den, glong end_val_num, glong end_val_den)
{
  glong size = end - start;
  glong i;
  glong a = 0, b = 0;
  gdouble start_val = start_val_num * 1.0 / (1.0 * start_val_den);
  gdouble end_val = end_val_num * 1.0 / (1.0 * end_val_den);
  mpz_t numa, numb, den;
  mpz_init (numa);
  mpz_init (numb);
  mpz_init (den);
  NcmGridSection sec;

  g_assert ((end > start ) && (end <= grid->nnodes));
  g_assert (end_val > start_val);

  if (grid->sections->len > 0)
  {
    gulong lasti = grid->sections->len - 1;
    g_assert (start == g_array_index(grid->sections,NcmGridSection,lasti).end);
    g_assert (start_val >= g_array_index(grid->sections,NcmGridSection,lasti).end_val);
  }
  else
    g_assert (start == 0);

  sec.incl = incl;
  sec.start = start;
  sec.end = end;
  sec.start_val = start_val;
  sec.end_val = end_val;
  g_array_append_val(grid->sections, sec);

  switch (incl)
  {
    case NCM_GRID_NODES_START:
      a = 0; b =  0; break;
    case NCM_GRID_NODES_END:
      a = 1; b =  0; break;
    case NCM_GRID_NODES_BOTH:
      a = 0; b = -1; break;
    case NCM_GRID_NODES_NONE:
      a = 1; b = +1; break;
  }

  mpz_set_si (numa, start_val_num);
  mpz_mul_si (numa, numa, end_val_den);
  mpz_mul_si (numa, numa, size + b);

  mpz_set_si (numb, end_val_num);
  mpz_mul_si (numb, numb, start_val_den);
  mpz_set_si (den, start_val_num);
  mpz_mul_si (den, den, end_val_den);
  mpz_sub (numb, numb, den);

  mpz_set_si (den, start_val_den);
  mpz_mul_si (den, den, end_val_den);
  mpz_mul_si (den, den, size + b);

  for (i = 0; i < size; i++)
  {
    const mpz_ptr node_num = mpq_numref (grid->nodes[start + i]);

    mpz_set (node_num, numa);
    mpz_addmul_ui (node_num, numb, i + a);
    mpz_set (mpq_denref (grid->nodes[start + i]), den);
    mpq_canonicalize (grid->nodes[start + i]);
    //mpfr_printf ("# [%d] %Qd %.15f\n", start + i, grid->nodes[start + i], mpq_get_d (grid->nodes[start + i]));
  }

  mpz_clear (numa);
  mpz_clear (numb);
  mpz_clear (den);
}

/**
 * ncm_grid_get_double_array:
 * @grid: a #NcmGrid
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble *
ncm_grid_get_double_array (NcmGrid *grid)
{
  if (grid->data == NULL)
  {
    guint i;
    grid->data = g_slice_alloc (sizeof (gdouble) * grid->nnodes);
    for (i = 0; i < grid->nnodes; i++)
      grid->data[i] = mpq_get_d (grid->nodes[i]);
  }
  return grid->data;
}

/**
 * ncm_grid_get_node_d:
 * @grid: a #NcmGrid
 * @i: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
ncm_grid_get_node_d (NcmGrid *grid, gulong i)
{
  if (grid->data == NULL)
    ncm_grid_get_double_array (grid);
  return grid->data[i];
}

/**
 * ncm_grid_write:
 * @grid: a #NcmGrid
 * @f: FIXME
 *
 * FIXME
 *
*/
void
ncm_grid_write (NcmGrid *grid, FILE *f)
{
  guint32 i;
  NCM_WRITE_UINT32 (f, grid->nnodes);
  for (i = 0; i < grid->nnodes; i++)
    ncm_mpq_out_raw (f, grid->nodes[i]);
}

/**
 * ncm_grid_read: (skip)
 * @f: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmGrid *
ncm_grid_read (FILE *f)
{
  guint32 i;
  NcmGrid *grid;
  NCM_READ_UINT32 (f, i);
  grid = ncm_grid_new (i);
  for (i = 0; i < grid->nnodes; i++)
    ncm_mpq_inp_raw (grid->nodes[i], f);
  return grid;
}
