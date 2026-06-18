/***************************************************************************
 *            nc_multiplicity_func_st.h
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_MULTIPLICITY_FUNC_ST_H_
#define _NC_MULTIPLICITY_FUNC_ST_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_ST             (nc_multiplicity_func_st_get_type ())

G_DECLARE_FINAL_TYPE (NcMultiplicityFuncST, nc_multiplicity_func_st, NC, MULTIPLICITY_FUNC_ST, NcMultiplicityFunc)


NcMultiplicityFuncST *nc_multiplicity_func_st_new (void);
NcMultiplicityFuncST *nc_multiplicity_func_st_ref (NcMultiplicityFuncST *mst);

void nc_multiplicity_func_st_free (NcMultiplicityFuncST *mst);
void nc_multiplicity_func_st_clear (NcMultiplicityFuncST **mst);

void nc_multiplicity_func_st_set_A (NcMultiplicityFuncST *mst, gdouble A);
gdouble nc_multiplicity_func_st_get_A (const NcMultiplicityFuncST *mst);

void nc_multiplicity_func_st_set_b (NcMultiplicityFuncST *mst, gdouble b);
gdouble nc_multiplicity_func_st_get_b (const NcMultiplicityFuncST *mst);

void nc_multiplicity_func_st_set_p (NcMultiplicityFuncST *mst, gdouble p);
gdouble nc_multiplicity_func_st_get_p (const NcMultiplicityFuncST *mst);

void nc_multiplicity_func_st_set_delta_c (NcMultiplicityFuncST *mst, gdouble delta_c);
gdouble nc_multiplicity_func_st_get_delta_c (const NcMultiplicityFuncST *mst);

G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_ST_H_ */

