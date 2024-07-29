/***************************************************************************
 *            ncm_integral1d_ptr.h
 *
 *  Mon December 05 17:23:53 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_integral1d_ptr.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_INTEGRAL1D_PTR_H_
#define _NCM_INTEGRAL1D_PTR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_integral1d.h>

G_BEGIN_DECLS

#define NCM_TYPE_INTEGRAL1D_PTR (ncm_integral1d_ptr_get_type ())

G_DECLARE_FINAL_TYPE (NcmIntegral1dPtr, ncm_integral1d_ptr, NCM, INTEGRAL1D_PTR, NcmIntegral1d)

typedef gdouble (*NcmIntegral1dPtrF) (gpointer userdata, const gdouble x, const gdouble w);

NcmIntegral1dPtr *ncm_integral1d_ptr_new (NcmIntegral1dPtrF F, GDestroyNotify userfree);
NcmIntegral1dPtr *ncm_integral1d_ptr_new_full (NcmIntegral1dPtrF F, GDestroyNotify userfree, const gdouble reltol, const gdouble abstol, const guint partition, const guint rule);

NcmIntegral1dPtr *ncm_integral1d_ptr_ref (NcmIntegral1dPtr *int1d_ptr);
void ncm_integral1d_ptr_free (NcmIntegral1dPtr *int1d_ptr);
void ncm_integral1d_ptr_clear (NcmIntegral1dPtr **int1d_ptr);

void ncm_integral1d_ptr_set_userdata (NcmIntegral1dPtr *int1d_ptr, gpointer userdata);

G_END_DECLS

#endif /* _NCM_INTEGRAL1D_PTR_H_ */

