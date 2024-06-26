/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_spline_rbf.h
 *
 *  Fri April 06 20:44:33 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * ncm_spline_rbf.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_SPLINE_RBF_H_
#define _NCM_SPLINE_RBF_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPLINE_RBF (ncm_spline_rbf_get_type ())

G_DECLARE_FINAL_TYPE (NcmSplineRBF, ncm_spline_rbf, NCM, SPLINE_RBF, NcmSpline)

/**
 * NcmSplineRBFType:
 * @NCM_SPLINE_RBF_TYPE_POSDEF_GAUSS: Gaussian [RBF](https://en.wikipedia.org/wiki/Radial_basis_function) interpolation method for positive defined functions.
 * @NCM_SPLINE_RBF_TYPE_GAUSS: Gaussian [RBF](https://en.wikipedia.org/wiki/Radial_basis_function) interpolation method for any kind of function.
 *
 * Enumeration to choose which Gaussian [RBF](https://en.wikipedia.org/wiki/Radial_basis_function) interpolation method to be applied by the object.
 *
 */
typedef enum _NcmSplineRBFType
{
  NCM_SPLINE_RBF_TYPE_POSDEF_GAUSS = 0,
  NCM_SPLINE_RBF_TYPE_GAUSS,
  /* < private > */
  NCM_SPLINE_RBF_TYPE_LEN, /*< skip >*/
} NcmSplineRBFType;

NcmSplineRBF *ncm_spline_rbf_new (NcmSplineRBFType type_id);
NcmSplineRBF *ncm_spline_rbf_ref (NcmSplineRBF *rbf);

void ncm_spline_rbf_free (NcmSplineRBF *rbf);
void ncm_spline_rbf_clear (NcmSplineRBF **rbf);

void ncm_spline_rbf_set_type (NcmSplineRBF *rbf, NcmSplineRBFType type_id);
void ncm_spline_rbf_set_shape_params (NcmSplineRBF *rbf, NcmVector *shape_params);

G_END_DECLS

#endif /* _NCM_SPLINE_RBF_H_ */

