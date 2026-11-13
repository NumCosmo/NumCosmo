/***************************************************************************
 *            ncm_spectral_private.h
 *
 *  Thu Nov 12 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

#ifndef _NCM_SPECTRAL_PRIVATE_H_
#define _NCM_SPECTRAL_PRIVATE_H_

#include <glib.h>
#include "ncm/algebra/ncm_spectral.h"

G_BEGIN_DECLS

void _ncm_spectral_transform (NcmSpectral *spectral, guint N, guint howmany, const gdouble *vals, gdouble *coeffs);

G_END_DECLS

#endif /* _NCM_SPECTRAL_PRIVATE_H_ */

