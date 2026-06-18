/***************************************************************************
 *            nc_multiplicity_func_bocquet.h
 *
 *  Mon Aug 16 14:39:25 2021
 *  Copyright  2021  Mariana Penna Lima and Cinthia Nunes de Lima
 *  <pennalima@gmail.com>, <cinthia.n.lima@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2021 Mariana Penna Lima and Cinthia Nunes de Lima
 * <pennalima@gmail.com>, <cinthia.n.lima@uel.br>
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

#ifndef _NC_MULTIPLICITY_FUNC_BOCQUET_H_
#define _NC_MULTIPLICITY_FUNC_BOCQUET_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_BOCQUET             (nc_multiplicity_func_bocquet_get_type ())

G_DECLARE_FINAL_TYPE (NcMultiplicityFuncBocquet, nc_multiplicity_func_bocquet, NC, MULTIPLICITY_FUNC_BOCQUET, NcMultiplicityFunc)


/**
 * NcMultiplicityFuncBocquetSim:
 * @NC_MULTIPLICITY_FUNC_BOCQUET_SIM_DM: simulation - Dark Matter only
 * @NC_MULTIPLICITY_FUNC_BOCQUET_SIM_HYDRO: hydrodynamical simulation (baryon effects) 
 *
 */
typedef enum _NcMultiplicityFuncBocquetSim
{
  NC_MULTIPLICITY_FUNC_BOCQUET_SIM_DM = 0,
  NC_MULTIPLICITY_FUNC_BOCQUET_SIM_HYDRO,
  /* < private > */
  NC_MULTIPLICITY_FUNC_BOCQUET_SIM_LEN, /*< skip >*/
} NcMultiplicityFuncBocquetSim;

NcMultiplicityFuncBocquet *nc_multiplicity_func_bocquet_new (void);
NcMultiplicityFuncBocquet *nc_multiplicity_func_bocquet_new_full (NcMultiplicityFuncMassDef mdef, NcMultiplicityFuncBocquetSim sim, gdouble Delta);
NcMultiplicityFuncBocquet *nc_multiplicity_func_bocquet_ref (NcMultiplicityFuncBocquet *mb);

void nc_multiplicity_func_bocquet_free (NcMultiplicityFuncBocquet *mb);
void nc_multiplicity_func_bocquet_clear (NcMultiplicityFuncBocquet **mb);

void nc_multiplicity_func_bocquet_set_sim (NcMultiplicityFuncBocquet *mb, NcMultiplicityFuncBocquetSim sim);
NcMultiplicityFuncBocquetSim nc_multiplicity_func_bocquet_get_sim (const NcMultiplicityFuncBocquet *mb);

G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_BOCQUET_H_ */
