/***************************************************************************
 *            nc_powspec_ml.h
 *
 *  Thu February 18 12:32:25 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_powspec_ml.h
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

#ifndef _NC_POWSPEC_ML_H_
#define _NC_POWSPEC_ML_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_powspec.h>

G_BEGIN_DECLS

#define NC_TYPE_POWSPEC_ML (nc_powspec_ml_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcPowspecML, nc_powspec_ml, NC, POWSPEC_ML, NcmPowspec)

struct _NcPowspecMLClass
{
  /*< private > */
  NcmPowspecClass parent_class;

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[18];
};

NcPowspecML *nc_powspec_ml_ref (NcPowspecML *ps_ml);

void nc_powspec_ml_free (NcPowspecML *ps_ml);
void nc_powspec_ml_clear (NcPowspecML **ps_ml);

G_END_DECLS

#endif /* _NC_POWSPEC_ML_H_ */

