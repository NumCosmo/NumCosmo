/***************************************************************************
 *            ncm_mset_trans_kern_flat.h
 *
 *  Thu October 02 13:36:55 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_trans_kern_flat.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_MSET_TRANS_KERN_FLAT_H_
#define _NCM_MSET_TRANS_KERN_FLAT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset_trans_kern.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET_TRANS_KERN_FLAT (ncm_mset_trans_kern_flat_get_type ())

G_DECLARE_FINAL_TYPE (NcmMSetTransKernFlat, ncm_mset_trans_kern_flat, NCM, MSET_TRANS_KERN_FLAT, NcmMSetTransKern)

NcmMSetTransKernFlat *ncm_mset_trans_kern_flat_new (void);

G_END_DECLS

#endif /* _NCM_MSET_TRANS_KERN_FLAT_H_ */

