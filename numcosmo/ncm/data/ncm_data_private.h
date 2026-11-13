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
 * Sets the value returned by ncm_data_is_resampling(). ncm_dataset_resample() uses it to
 * prepare every block as resampling before the blocks are resampled.
 */
void _ncm_data_set_resampling (NcmData *data, gboolean resampling);

G_END_DECLS

#endif /* _NCM_DATA_PRIVATE_H_ */

