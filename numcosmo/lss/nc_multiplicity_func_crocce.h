/***************************************************************************
 *            nc_multiplicity_func_crocce.h
 *
 *  Wed Feb 15 13:36:09 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2017 <pennalima@gmail.com>
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

#ifndef _NC_MULTIPLICITY_FUNC_CROCCE_H_
#define _NC_MULTIPLICITY_FUNC_CROCCE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_CROCCE             (nc_multiplicity_func_crocce_get_type ())

G_DECLARE_FINAL_TYPE (NcMultiplicityFuncCrocce, nc_multiplicity_func_crocce, NC, MULTIPLICITY_FUNC_CROCCE, NcMultiplicityFunc)


NcMultiplicityFuncCrocce *nc_multiplicity_func_crocce_new (void);
NcMultiplicityFuncCrocce *nc_multiplicity_func_crocce_ref (NcMultiplicityFuncCrocce *mc);

void nc_multiplicity_func_crocce_free (NcMultiplicityFuncCrocce *mc);
void nc_multiplicity_func_crocce_clear (NcMultiplicityFuncCrocce **mc);

G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_CROCCE_H_ */

