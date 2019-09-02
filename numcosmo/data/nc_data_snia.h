/***************************************************************************
 *            nc_data_snia.h
 *
 *  Mon December 10 00:20:48 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NC_DATA_SNIA_H_
#define _NC_DATA_SNIA_H_

#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data.h>
#include <numcosmo/nc_distance.h>

G_BEGIN_DECLS

/**
 * NcDataSNIAId:
 * @NC_DATA_SNIA_SIMPLE_GOLD_157: FIXME
 * @NC_DATA_SNIA_SIMPLE_GOLD_182: FIXME
 * @NC_DATA_SNIA_SIMPLE_GOLD_182_FULL: FIXME
 * @NC_DATA_SNIA_SIMPLE_ESSENCE: FIXME
 * @NC_DATA_SNIA_SIMPLE_LEGACY: FIXME
 * @NC_DATA_SNIA_SIMPLE_UNION: FIXME
 * @NC_DATA_SNIA_SIMPLE_CfA3: FIXME
 * @NC_DATA_SNIA_SIMPLE_UNION2: FIXME
 * @NC_DATA_SNIA_SIMPLE_UNION2_1: FIXME
 * @NC_DATA_SNIA_SIMPLE_SDSS_EMILLE: FIXME
 * @NC_DATA_SNIA_COV_SNLS3_SYS_STAT: FIXME
 * @NC_DATA_SNIA_COV_SNLS3_STAT_ONLY: FIXME
 * @NC_DATA_SNIA_COV_JLA_SNLS3_SDSS_SYS_STAT: JLA SNIa compilation, reduced covariance already includes intrisic, redshift and peculiar velocity variances. 
 * @NC_DATA_SNIA_COV_JLA_SNLS3_SDSS_SYS_STAT_CMPL: JLA SNIa compilation, complete covariance does *not* include intrisic, redshift and peculiar velocity variances.
 * @NC_DATA_SNIA_COV_PANTHEON: Pantheon SNIa compilation, reduced covariance already includes intrisic, redshift and peculiar velocity variances.
 * 
 * FIXME
 * 
 */
typedef enum _NcDataSNIAId
{
  NC_DATA_SNIA_SIMPLE_GOLD_157 = 0,
  NC_DATA_SNIA_SIMPLE_GOLD_182,
  NC_DATA_SNIA_SIMPLE_GOLD_182_FULL,
  NC_DATA_SNIA_SIMPLE_ESSENCE,
  NC_DATA_SNIA_SIMPLE_LEGACY,
  NC_DATA_SNIA_SIMPLE_UNION,
  NC_DATA_SNIA_SIMPLE_CfA3,
  NC_DATA_SNIA_SIMPLE_UNION2,
  NC_DATA_SNIA_SIMPLE_UNION2_1,
  NC_DATA_SNIA_SIMPLE_SDSS_EMILLE,
  NC_DATA_SNIA_COV_SNLS3_SYS_STAT, 
  NC_DATA_SNIA_COV_SNLS3_STAT_ONLY,
  NC_DATA_SNIA_COV_JLA_SNLS3_SDSS_SYS_STAT,
  NC_DATA_SNIA_COV_JLA_SNLS3_SDSS_SYS_STAT_CMPL, 
  NC_DATA_SNIA_COV_PANTHEON,   
  /* < private > */
  NC_DATA_SNIA_LEN, /*< skip >*/
} NcDataSNIAId;

#define NC_DATA_SNIA_SIMPLE_START NC_DATA_SNIA_SIMPLE_GOLD_157
#define NC_DATA_SNIA_SIMPLE_END NC_DATA_SNIA_SIMPLE_SDSS_EMILLE
#define NC_DATA_SNIA_SIMPLE_LEN ((NC_DATA_SNIA_SIMPLE_END) - (NC_DATA_SNIA_SIMPLE_START) + 1)

#define NC_DATA_SNIA_COV_START NC_DATA_SNIA_COV_SNLS3_SYS_STAT
#define NC_DATA_SNIA_COV_END NC_DATA_SNIA_COV_PANTHEON
#define NC_DATA_SNIA_COV_LEN ((NC_DATA_SNIA_COV_END) - (NC_DATA_SNIA_COV_START) + 1)

G_END_DECLS

#endif /* _NC_DATA_SNIA_H_ */
