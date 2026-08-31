/***************************************************************************
 *            nc_bbn_parthenope.h
 *
 *  Fri August 29 20:45:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_bbn_parthenope.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NC_BBN_PARTHENOPE_H_
#define _NC_BBN_PARTHENOPE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/bbn/nc_bbn.h>

G_BEGIN_DECLS

#define NC_TYPE_BBN_PARTHENOPE (nc_bbn_parthenope_get_type ())

G_DECLARE_FINAL_TYPE (NcBBNParthenope, nc_bbn_parthenope, NC, BBN_PARTHENOPE, NcBBN);

/**
 * NcBBNParthenopeTable:
 * @NC_BBN_PARTHENOPE_TABLE_PLANCK2017: PArthENoPE v1.2, neutron lifetime 880.2 s, the Planck 2017 assumptions
 * @NC_BBN_PARTHENOPE_TABLE_PLANCK2017_MARCUCCI: as above with the Marcucci et al. $d(p,\gamma){}^3\mathrm{He}$ rate
 * @NC_BBN_PARTHENOPE_TABLE_LEGACY: the older table NumCosmo shipped before 2017
 *
 * Which tabulation of $Y_p(\omega_b, \Delta N_\mathrm{eff})$ to interpolate.
 *
 * The choice is not cosmetic: at $\omega_b = 0.02237$, $\Delta N_\mathrm{eff} =
 * 0$ the legacy table gives $Y_p = 0.247840$ against $0.245395$ for
 * %NC_BBN_PARTHENOPE_TABLE_PLANCK2017, a 1% difference. Before this enum existed the
 * table was hardwired per cosmology class, so #NcHICosmoDE and the other
 * background models silently disagreed by that much.
 *
 * %NC_BBN_PARTHENOPE_TABLE_PLANCK2017_MARCUCCI differs from
 * %NC_BBN_PARTHENOPE_TABLE_PLANCK2017 by about $10^{-5}$ in $Y_p$: the
 * $d(p,\gamma){}^3\mathrm{He}$ rate it varies drives Deuterium, not Helium.
 *
 */
typedef enum _NcBBNParthenopeTable /*< prefix=NC_BBN_PARTHENOPE_TABLE >*/
{
  NC_BBN_PARTHENOPE_TABLE_PLANCK2017 = 0,
  NC_BBN_PARTHENOPE_TABLE_PLANCK2017_MARCUCCI,
  NC_BBN_PARTHENOPE_TABLE_LEGACY,
  /* < private > */
} NcBBNParthenopeTable;

NcBBNParthenope *nc_bbn_parthenope_new (void);
NcBBNParthenope *nc_bbn_parthenope_new_from_table (NcBBNParthenopeTable table);
NcBBNParthenope *nc_bbn_parthenope_ref (NcBBNParthenope *bbn_pn);

void nc_bbn_parthenope_free (NcBBNParthenope *bbn_pn);
void nc_bbn_parthenope_clear (NcBBNParthenope **bbn_pn);

void nc_bbn_parthenope_set_table (NcBBNParthenope *bbn_pn, NcBBNParthenopeTable table);
NcBBNParthenopeTable nc_bbn_parthenope_get_table (NcBBNParthenope *bbn_pn);

G_END_DECLS

#endif /* _NC_BBN_PARTHENOPE_H_ */

