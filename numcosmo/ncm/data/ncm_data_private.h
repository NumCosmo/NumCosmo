/***************************************************************************
 *            ncm_data_private.h
 *
 *  Fri November 13 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_data_private.h
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

#ifndef _NCM_DATA_PRIVATE_H_
#define _NCM_DATA_PRIVATE_H_

#include <glib.h>
#include "ncm/data/ncm_data.h"

G_BEGIN_DECLS

/*
 * Runs the prepare method of @data for an evaluation (@resampling FALSE), which requires
 * @data to be initialized, or for a resampling (@resampling TRUE), which does not.
 * ncm_data_is_resampling() returns @resampling while the method runs.
 */
void _ncm_data_prepare_full (NcmData *data, NcmMSet *mset, const gboolean resampling);

G_END_DECLS

#endif /* _NCM_DATA_PRIVATE_H_ */

