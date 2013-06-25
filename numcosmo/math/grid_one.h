/***************************************************************************
 *            grid_one.h
 *
 *  Mon Feb 22 15:05:39 2010
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

#ifndef _NC_GRID_ONE_H
#define _NC_GRID_ONE_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <stdio.h>
#ifndef NUMCOSMO_GIR_SCAN
#include <gmp.h>
#endif

G_BEGIN_DECLS

/**
 * NcmGridNodesEndPoints:
 * @NCM_GRID_NODES_START: FIXME
 * @NCM_GRID_NODES_END: FIXME
 * @NCM_GRID_NODES_BOTH: FIXME
 * @NCM_GRID_NODES_NONE: FIXME
 *
 * FIXME
 */
typedef enum _NcmGridNodesEndPoints
{
  NCM_GRID_NODES_START = 1,
  NCM_GRID_NODES_END,
  NCM_GRID_NODES_BOTH,
  NCM_GRID_NODES_NONE,
} NcmGridNodesEndPoints;

typedef struct _NcmGridSection NcmGridSection;

/**
 * NcmGridSection:
 *
 * FIXME
 */
struct _NcmGridSection
{
  /*< private >*/
  NcmGridNodesEndPoints incl;
  guint32 start;
  guint32 end;
  gdouble start_val;
  gdouble end_val;
};

typedef struct _NcmGrid NcmGrid;

/**
 * NcmGrid:
 *
 * FIXME
 */
struct _NcmGrid
{
  /*< private >*/
  GArray *sections;
  guint32 nnodes;
  mpq_t *nodes;
  gdouble *data;
};

#define NCM_GRID_GET_SEC(grid) (&(g_array_index((grid)->sections, NcmGridSection, 0)))

NcmGrid *ncm_grid_new (gulong nnodes);
NcmGrid *ncm_grid_new_from_sections (NcmGridSection *secs);
void ncm_grid_free (NcmGrid *grid, gboolean free_data);
gchar *ncm_grid_get_name (NcmGridSection *secs);
void ncm_grid_set_sections (NcmGrid *grid, NcmGridSection *secs);
void ncm_grid_set_nodes_d (NcmGrid *grid, NcmGridNodesEndPoints incl, guint32 start, guint32 end, gdouble start_val, gdouble end_val);
void ncm_grid_set_nodes_si (NcmGrid *grid, NcmGridNodesEndPoints incl, guint32 start, guint32 end, glong start_val_num, glong start_val_den, glong end_val_num, glong end_val_den);

gdouble ncm_grid_get_node_d (NcmGrid *grid, gulong i);
gdouble *ncm_grid_get_double_array (NcmGrid *grid);

void ncm_grid_write (NcmGrid *grid, FILE *f);
NcmGrid *ncm_grid_read (FILE *f);

G_END_DECLS

#endif /* _NC_GRID_ONE_H */
