/***************************************************************************
 *            nc_data_cmb.h
 *
 *  Thu November 22 21:19:23 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_DATA_CMB_H_
#define _NC_DATA_CMB_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data.h>
#include <numcosmo/nc_distance.h>

G_BEGIN_DECLS

/**
 * NcDataCMBDataType:
 * @NC_DATA_CMB_TYPE_TT: Temperature - Temperature correlation data
 * @NC_DATA_CMB_TYPE_EE: Mode E - Mode E correlation data
 * @NC_DATA_CMB_TYPE_BB: Mode B - Mode B correlation data
 * @NC_DATA_CMB_TYPE_TE: Temperature - Mode E correlation data
 * @NC_DATA_CMB_TYPE_TB: Temperature - Mode B correlation data
 * @NC_DATA_CMB_TYPE_EB: Mode E - Mode B correlation data
 * @NC_DATA_CMB_TYPE_PHIPHI: $\phi$ - $\phi$ correlation data
 * @NC_DATA_CMB_TYPE_ALL: All types above
 *
 * Each type corresponds to one CMB two-point correlation function data to be used:
 * TT (temperature only), EE (E-mode only), BB (B-mode only), TE (temperature and E-mode polarization),
 * TB (temperature and B-mode polarization), EB (E- and B-mode polarization) and $\phi\phi$ (CMB lensing).
 * The last option `ALL' selects all data listed above.
 *
 */
typedef enum _NcDataCMBDataType
{
  NC_DATA_CMB_TYPE_TT = 1 << 0,
  NC_DATA_CMB_TYPE_EE = 1 << 1,
  NC_DATA_CMB_TYPE_BB = 1 << 2,
  NC_DATA_CMB_TYPE_TE = 1 << 3,
  NC_DATA_CMB_TYPE_TB = 1 << 4,
  NC_DATA_CMB_TYPE_EB = 1 << 5,
  NC_DATA_CMB_TYPE_PHIPHI = 1 << 6,
  NC_DATA_CMB_TYPE_ALL = (1 << 7) - 1,
  /* < private > */
  NC_DATA_CMB_TYPE_LEN, /*< skip >*/
} NcDataCMBDataType;

/**
 * NcDataCMBId:
 * @NC_DATA_CMB_SHIFT_PARAM_WMAP3: three-year WMAP shift parameter
 * @NC_DATA_CMB_SHIFT_PARAM_WMAP5: five-year WMAP shift parameter
 * @NC_DATA_CMB_SHIFT_PARAM_WMAP7: seven-year WMAP shift parameter
 * @NC_DATA_CMB_DIST_PRIORS_WMAP5: five-year WMAP distance priors
 * @NC_DATA_CMB_DIST_PRIORS_WMAP7: seven-year WMAP distance priors
 * @NC_DATA_CMB_DIST_PRIORS_WMAP9: nine-year WMAP distance priors
 *
 * Wilkinson Microwave Anisotropy Probe (WMAP) distance priors and shift parameter estimates.
 *
 * The so-called distance priors comprises the location of the first acoustic peak $l_A$
 * and the shift parameter $R$, that is,
 * $$l_A \equiv (1 + z_\star) \frac{\pi D_A(z_\star)}{r_s(z_\star)}$$
 * and
 * $$R = \frac{\sqrt{\Omega_m H_0^2}}{c} (1 + z_\star) D_A(z_\star),$$
 * where $z_\star$ is the redshift at decoupling, $D_A(z)$ is the angular diameter distance, and
 * $r_s(z_\star)$ is the sound horizon size.
 *
 */
typedef enum _NcDataCMBId
{
  NC_DATA_CMB_SHIFT_PARAM_WMAP3 = 0,
  NC_DATA_CMB_SHIFT_PARAM_WMAP5,
  NC_DATA_CMB_SHIFT_PARAM_WMAP7,
  NC_DATA_CMB_DIST_PRIORS_WMAP5,
  NC_DATA_CMB_DIST_PRIORS_WMAP7,
  NC_DATA_CMB_DIST_PRIORS_WMAP9,
  /* < private > */
  NC_DATA_CMB_NSAMPLES, /*< skip >*/
} NcDataCMBId;

NcmData *nc_data_cmb_create (NcDistance *dist, NcDataCMBId id);

G_END_DECLS

#endif /* _NC_DATA_CMB_H_ */

