/***************************************************************************
 *            nc_multiplicity_func_tinker.h
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

#ifndef _NC_MULTIPLICITY_FUNC_TINKER_H_
#define _NC_MULTIPLICITY_FUNC_TINKER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_TINKER             (nc_multiplicity_func_tinker_get_type ())

G_DECLARE_FINAL_TYPE (NcMultiplicityFuncTinker, nc_multiplicity_func_tinker, NC, MULTIPLICITY_FUNC_TINKER, NcMultiplicityFunc)


NcMultiplicityFuncTinker *nc_multiplicity_func_tinker_new (void);
NcMultiplicityFuncTinker *nc_multiplicity_func_tinker_new_full (NcMultiplicityFuncMassDef mdef, gdouble Delta);
NcMultiplicityFuncTinker *nc_multiplicity_func_tinker_ref (NcMultiplicityFuncTinker *mt);

void nc_multiplicity_func_tinker_free (NcMultiplicityFuncTinker *mt);
void nc_multiplicity_func_tinker_clear (NcMultiplicityFuncTinker **mt);


void nc_multiplicity_func_tinker_set_linear_interp (NcMultiplicityFuncTinker *mulf, gboolean lin_interp);


G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_TINKER_H_ */
