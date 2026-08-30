/***************************************************************************
 *            nc_bbn_parametrized.h
 *
 *  Sat August 30 10:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_bbn_parametrized.h
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

#ifndef _NC_BBN_PARAMETRIZED_H_
#define _NC_BBN_PARAMETRIZED_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/bbn/nc_bbn.h>

G_BEGIN_DECLS

#define NC_TYPE_BBN_PARAMETRIZED (nc_bbn_parametrized_get_type ())

G_DECLARE_FINAL_TYPE (NcBBNParametrized, nc_bbn_parametrized, NC, BBN_PARAMETRIZED, NcBBN);

/**
 * NcBBNParametrizedSParams:
 * @NC_BBN_PARAMETRIZED_YP_4HE: primordial Helium-4 mass fraction $Y_p$
 *
 * Parameters of #NcBBNParametrized.
 *
 */
typedef enum _NcBBNParametrizedSParams /*< prefix=NC_BBN_PARAMETRIZED >*/
{
  NC_BBN_PARAMETRIZED_YP_4HE = 0,
  /* < private > */
  NC_BBN_PARAMETRIZED_SPARAM_LEN, /*< skip >*/
} NcBBNParametrizedSParams;

#define NC_BBN_PARAMETRIZED_DEFAULT_YP_4HE (0.24)

NcBBNParametrized *nc_bbn_parametrized_new (void);
NcBBNParametrized *nc_bbn_parametrized_ref (NcBBNParametrized *bbn_par);

void nc_bbn_parametrized_free (NcBBNParametrized *bbn_par);
void nc_bbn_parametrized_clear (NcBBNParametrized **bbn_par);

G_END_DECLS

#endif /* _NC_BBN_PARAMETRIZED_H_ */

