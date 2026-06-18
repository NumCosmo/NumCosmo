/***************************************************************************
 *            nc_powspec_ml_cbe.h
 *
 *  Tue April 05 10:42:32 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_powspec_ml_cbe.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_POWSPEC_ML_CBE_H_
#define _NC_POWSPEC_ML_CBE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_powspec_ml.h>
#include <numcosmo/nc_cbe.h>

G_BEGIN_DECLS

#define NC_TYPE_POWSPEC_ML_CBE             (nc_powspec_ml_cbe_get_type ())

G_DECLARE_FINAL_TYPE (NcPowspecMLCBE, nc_powspec_ml_cbe, NC, POWSPEC_ML_CBE, NcPowspecML)


NcPowspecMLCBE *nc_powspec_ml_cbe_new (void);
NcPowspecMLCBE *nc_powspec_ml_cbe_new_full (NcCBE *cbe);

void nc_powspec_ml_cbe_set_cbe (NcPowspecMLCBE *ps_cbe, NcCBE *cbe);

NcCBE *nc_powspec_ml_cbe_peek_cbe (NcPowspecMLCBE *ps_cbe);

void nc_powspec_ml_cbe_set_intern_k_min (NcPowspecMLCBE *ps_cbe, const gdouble k_min);
void nc_powspec_ml_cbe_set_intern_k_max (NcPowspecMLCBE *ps_cbe, const gdouble k_max);
gdouble nc_powspec_ml_cbe_get_intern_k_min (NcPowspecMLCBE *ps_cbe);
gdouble nc_powspec_ml_cbe_get_intern_k_max (NcPowspecMLCBE *ps_cbe);

#define NC_POWSPEC_ML_CBE_INTERN_KMIN (1.0e-5)
#define NC_POWSPEC_ML_CBE_INTERN_KMAX (1.0e1)

G_END_DECLS

#endif /* _NC_POWSPEC_ML_CBE_H_ */

