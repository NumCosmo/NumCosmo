/***************************************************************************
 *            ncm_sphere_xi.h
 *
 *  Fri Oct 16 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sphere_xi.h
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

#ifndef _NCM_SPHERE_XI_H_
#define _NCM_SPHERE_XI_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPHERE_XI (ncm_sphere_xi_get_type ())

G_DECLARE_FINAL_TYPE (NcmSphereXi, ncm_sphere_xi, NCM, SPHERE_XI, GObject)

NcmSphereXi *ncm_sphere_xi_new (NcmVector * theta_edges, guint lmax);
NcmSphereXi *ncm_sphere_xi_ref (NcmSphereXi *sxi);
void ncm_sphere_xi_free (NcmSphereXi *sxi);
void ncm_sphere_xi_clear (NcmSphereXi **sxi);

void ncm_sphere_xi_set_theta_edges (NcmSphereXi *sxi, NcmVector *theta_edges);
NcmVector *ncm_sphere_xi_peek_theta_edges (NcmSphereXi *sxi);
void ncm_sphere_xi_set_lmax (NcmSphereXi *sxi, guint lmax);
guint ncm_sphere_xi_get_lmax (NcmSphereXi *sxi);
guint ncm_sphere_xi_get_n_bins (NcmSphereXi *sxi);

void ncm_sphere_xi_prepare (NcmSphereXi *sxi);
NcmMatrix *ncm_sphere_xi_peek_kernel (NcmSphereXi *sxi);
void ncm_sphere_xi_apply (NcmSphereXi *sxi, NcmMatrix *C_ells, const guint lmin, NcmMatrix *xi);
void ncm_sphere_xi_apply_vector (NcmSphereXi *sxi, NcmVector *C_ell, const guint lmin, NcmVector *xi);

void ncm_sphere_xi_to_chebyshev (NcmSphereXi *sxi, NcmMatrix *C_ells, const guint lmin, NcmMatrix *c);
void ncm_sphere_xi_eval_chebyshev (NcmSphereXi *sxi, NcmMatrix *c, NcmVector *theta, NcmMatrix *xi);
void ncm_sphere_xi_eval_chebyshev_grid (NcmSphereXi *sxi, NcmMatrix *c, NcmMatrix *xi);

G_END_DECLS

#endif /* _NCM_SPHERE_XI_H_ */

