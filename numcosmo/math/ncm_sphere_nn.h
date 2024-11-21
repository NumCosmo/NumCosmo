/***************************************************************************
 *            ncm_sphere_nn.h
 *
 *  Wed Nov 20 19:23:40 2024
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sphere_nn.h
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_SPHERE_NN_H_
#define _NCM_SPHERE_NN_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_quaternion.h>
#include <numcosmo/math/ncm_spline.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_vector_float.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_SPHERE_NN (ncm_sphere_nn_get_type ())

G_DECLARE_FINAL_TYPE (NcmSphereNN, ncm_sphere_nn, NCM, SPHERE_NN, GObject)

NcmSphereNN *ncm_sphere_nn_new ();
NcmSphereNN *ncm_sphere_nn_ref (NcmSphereNN *snn);
void ncm_sphere_nn_free (NcmSphereNN *snn);
void ncm_sphere_nn_clear (NcmSphereNN **snn);

void ncm_sphere_nn_insert (NcmSphereNN *snn, const gdouble r, const gdouble theta, const gdouble phi);
void ncm_sphere_nn_insert_array (NcmSphereNN *snn, GArray *r, GArray *theta, GArray *phi);
void ncm_sphere_nn_get (NcmSphereNN *snn, const gint64 i, gdouble *r, gdouble *theta, gdouble *phi);
gint64 ncm_sphere_nn_get_n (NcmSphereNN *snn);

void ncm_sphere_nn_rebuild (NcmSphereNN *snn);
GArray *ncm_sphere_nn_knn_search (NcmSphereNN *snn, const gdouble r, const gdouble theta, const gdouble phi, const gint64 k);
void ncm_sphere_nn_knn_search_distances (NcmSphereNN *snn, const gdouble r, const gdouble theta, const gdouble phi, const gint64 k, GArray **distances, GArray **indices);

G_END_DECLS

#endif /* _NCM_SPHERE_NN_H_ */

