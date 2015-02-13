/***************************************************************************
 *            nc_data_cmb.c
 *
 *  Thu November 22 21:19:23 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

/**
 * SECTION:nc_data_cmb
 * @title: NcDataCMB
 * @short_description: Helper function for instantiating CMB data.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cmb.h"

#include "data/nc_data_cmb_shift_param.h"
#include "data/nc_data_cmb_dist_priors.h"

/**
 * nc_data_cmb_create:
 * @dist: FIXME
 * @id: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmData *
nc_data_cmb_create (NcDistance *dist, NcDataCMBId id)
{
  switch (id)
  {
    case NC_DATA_CMB_SHIFT_PARAM_WMAP3:
    case NC_DATA_CMB_SHIFT_PARAM_WMAP5:
    case NC_DATA_CMB_SHIFT_PARAM_WMAP7:
      return nc_data_cmb_shift_param_new (dist, id);
      break;
    case NC_DATA_CMB_DIST_PRIORS_WMAP5:
    case NC_DATA_CMB_DIST_PRIORS_WMAP7:
    case NC_DATA_CMB_DIST_PRIORS_WMAP9:
      return nc_data_cmb_dist_priors_new (dist, id);
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}
