/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_data_funnel.h
 *
 *  Wed May 12 21:31:43 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_data_funnel.h
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_DATA_FUNNEL_H_
#define _NCM_DATA_FUNNEL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_cov.h>

G_BEGIN_DECLS

#define NCM_TYPE_DATA_FUNNEL (ncm_data_funnel_get_type ())

G_DECLARE_FINAL_TYPE (NcmDataFunnel, ncm_data_funnel, NCM, DATA_FUNNEL, NcmData)

NcmDataFunnel *ncm_data_funnel_new (void);
NcmDataFunnel *ncm_data_funnel_ref (NcmDataFunnel *dfu);
void ncm_data_funnel_free (NcmDataFunnel *dfu);
void ncm_data_funnel_clear (NcmDataFunnel **dfu);

G_END_DECLS

#endif /* _NCM_DATA_FUNNEL_H_ */

