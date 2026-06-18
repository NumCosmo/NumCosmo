/***************************************************************************
 *            nc_multiplicity_func_tinker_mean_normalized.h
 *
 *  Thu July 03 12:31:21 2014
 *  Copyright  2014  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_multiplicity_func_tinker_mean_normalized.h
 * Copyright (C) 2014 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED_H_
#define _NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED             (nc_multiplicity_func_tinker_mean_normalized_get_type ())

G_DECLARE_FINAL_TYPE (NcMultiplicityFuncTinkerMeanNormalized, nc_multiplicity_func_tinker_mean_normalized, NC, MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED, NcMultiplicityFunc)


NcMultiplicityFuncTinkerMeanNormalized *nc_multiplicity_func_tinker_mean_normalized_new (void);
NcMultiplicityFuncTinkerMeanNormalized *nc_multiplicity_func_tinker_mean_normalized_ref (NcMultiplicityFuncTinkerMeanNormalized *mt10);

void nc_multiplicity_func_tinker_mean_normalized_free (NcMultiplicityFuncTinkerMeanNormalized *mt10);
void nc_multiplicity_func_tinker_mean_normalized_clear (NcMultiplicityFuncTinkerMeanNormalized **mt10);


G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED_H_ */

