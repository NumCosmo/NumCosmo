/***************************************************************************
 *            ncm_powspec_spline2d.h
 *
 *  Tue February 16 17:01:03 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_powspec_spline2d.h
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

#ifndef _NCM_POWSPEC_SPLINE2D_H_
#define _NCM_POWSPEC_SPLINE2D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/math/ncm_integral1d_ptr.h>
#include <numcosmo/math/ncm_powspec.h>
#include <numcosmo/math/ncm_spline2d.h>

G_BEGIN_DECLS

#define NCM_TYPE_POWSPEC_SPLINE2D (ncm_powspec_spline2d_get_type ())

G_DECLARE_FINAL_TYPE (NcmPowspecSpline2d, ncm_powspec_spline2d, NCM, POWSPEC_SPLINE2D, NcmPowspec)

NcmPowspecSpline2d *ncm_powspec_spline2d_new (NcmSpline2d * spline2d);
NcmPowspecSpline2d *ncm_powspec_spline2d_ref (NcmPowspecSpline2d *ps_s2d);

void ncm_powspec_spline2d_free (NcmPowspecSpline2d *ps_s2d);
void ncm_powspec_spline2d_clear (NcmPowspecSpline2d **ps_s2d);

void ncm_powspec_spline2d_set_spline2d (NcmPowspecSpline2d *ps_s2d, NcmSpline2d *spline2d);
NcmSpline2d *ncm_powspec_spline2d_peek_spline2d (NcmPowspecSpline2d *ps_s2d);

G_END_DECLS

#endif /* _NCM_POWSPEC_SPLINE2D_H_ */

