/***************************************************************************
 *            nc_data_hubble.h
 *
 *  Thu Apr 22 14:35:37 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
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

#ifndef _NC_DATA_HUBBLE_H_
#define _NC_DATA_HUBBLE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_diag.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_HUBBLE             (nc_data_hubble_get_type ())

G_DECLARE_FINAL_TYPE (NcDataHubble, nc_data_hubble, NC, DATA_HUBBLE, NcmDataGaussDiag);

/**
 * NcDataHubbleId:
 * @NC_DATA_HUBBLE_SIMON2005: FIXME
 * @NC_DATA_HUBBLE_CABRE: FIXME
 * @NC_DATA_HUBBLE_STERN2009: FIXME
 * @NC_DATA_HUBBLE_MORESCO2012_BC03: [Moresco et al. (2012)][XMoresco2012]
 * @NC_DATA_HUBBLE_MORESCO2012_MASTRO: [Moresco et al. (2012)][XMoresco2012]
 * @NC_DATA_HUBBLE_MORESCO2015: [Moresco (2015)][XMoresco2015]
 * @NC_DATA_HUBBLE_MORESCO2016_DR9_BC03: [Moresco et al. (2016)][XMoresco2016]
 * @NC_DATA_HUBBLE_MORESCO2016_DR9_MASTRO: [Moresco et al. (2016)][XMoresco2016]
 * @NC_DATA_HUBBLE_BUSCA2013_BAO_WMAP: FIXME
 * @NC_DATA_HUBBLE_RIESS2008_HST: FIXME
 * @NC_DATA_HUBBLE_ZHANG2012: FIXME
 * @NC_DATA_HUBBLE_RIESS2016_HST_WFC3: FIXME
 * @NC_DATA_HUBBLE_RATSIMBAZAFY2017: Ratsimbazafy et al. 2017 -- arXiv:1702.00418
 * @NC_DATA_HUBBLE_GOMEZ_VALENT_COMP2018: [Gomez-Valent et al. (2018)][XGomez-Valent2018]
 * @NC_DATA_HUBBLE_RIESS2018: FIXME
 * @NC_DATA_HUBBLE_BORGHI2022: Borghi et al. 2022 -- arXiv:2110.04304
 * @NC_DATA_HUBBLE_JIAO2023: Jiao et al. 2023 -- arXiv:2205.05701
 * @NC_DATA_HUBBLE_JIMENEZ2023: Jimenez et al. 2023 -- arXiv:2306.11425
 * @NC_DATA_HUBBLE_TOMASETTI2023: Tomasetti et al. 2023 -- arXiv:2305.16387
 *
 * FIXME
 */
typedef enum _NcDataHubbleId
{
  NC_DATA_HUBBLE_SIMON2005 = 0,
  NC_DATA_HUBBLE_CABRE,
  NC_DATA_HUBBLE_STERN2009,
  NC_DATA_HUBBLE_MORESCO2012_BC03,
  NC_DATA_HUBBLE_MORESCO2012_MASTRO,
  NC_DATA_HUBBLE_MORESCO2015,
  NC_DATA_HUBBLE_MORESCO2016_DR9_BC03,
  NC_DATA_HUBBLE_MORESCO2016_DR9_MASTRO,
  NC_DATA_HUBBLE_BUSCA2013_BAO_WMAP,
  NC_DATA_HUBBLE_RIESS2008_HST,
  NC_DATA_HUBBLE_ZHANG2012,
  NC_DATA_HUBBLE_RIESS2016_HST_WFC3,
  NC_DATA_HUBBLE_RATSIMBAZAFY2017,
  NC_DATA_HUBBLE_GOMEZ_VALENT_COMP2018,
  NC_DATA_HUBBLE_RIESS2018,
  NC_DATA_HUBBLE_BORGHI2022,
  NC_DATA_HUBBLE_JIAO2023,
  NC_DATA_HUBBLE_JIMENEZ2023,
  NC_DATA_HUBBLE_TOMASETTI2023,
  /* < private > */
  NC_DATA_HUBBLE_NSAMPLES, /*< skip >*/
} NcDataHubbleId;

NcDataHubble *nc_data_hubble_new_empty (void);
NcDataHubble *nc_data_hubble_new_from_file (const gchar *filename);
NcDataHubble *nc_data_hubble_new_from_id (NcDataHubbleId id);

void nc_data_hubble_set_sample (NcDataHubble *hubble, NcDataHubbleId id);

G_END_DECLS

#endif /* _NC_DATA_HUBBLE_H_ */

