/***************************************************************************
 *            ncm_powspec.h
 *
 *  Tue February 16 17:01:03 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_powspec.h
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

#ifndef _NCM_POWSPEC_H_
#define _NCM_POWSPEC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/math/ncm_integral1d_ptr.h>
#include <numcosmo/math/ncm_spline2d.h>

G_BEGIN_DECLS

#define NCM_TYPE_POWSPEC (ncm_powspec_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmPowspec, ncm_powspec, NCM, POWSPEC, GObject)

struct _NcmPowspecClass
{
  /*< private > */
  GObjectClass parent_class;

  void (*prepare) (NcmPowspec *powspec, NcmModel *model);
  gdouble (*eval) (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k);
  void (*eval_vec) (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk);
  void (*get_nknots) (NcmPowspec *powspec, guint *Nz, guint *Nk);
  NcmSpline2d *(*get_spline_2d) (NcmPowspec *powspec, NcmModel *model);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[13];
};

NcmPowspec *ncm_powspec_ref (NcmPowspec *powspec);

void ncm_powspec_free (NcmPowspec *powspec);
void ncm_powspec_clear (NcmPowspec **powspec);

void ncm_powspec_set_zi (NcmPowspec *powspec, const gdouble zi);
void ncm_powspec_set_zf (NcmPowspec *powspec, const gdouble zf);
void ncm_powspec_set_kmin (NcmPowspec *powspec, const gdouble kmin);
void ncm_powspec_set_kmax (NcmPowspec *powspec, const gdouble kmax);
void ncm_powspec_set_reltol_spline (NcmPowspec *powspec, const gdouble reltol);

void ncm_powspec_require_zi (NcmPowspec *powspec, const gdouble zi);
void ncm_powspec_require_zf (NcmPowspec *powspec, const gdouble zf);
void ncm_powspec_require_kmin (NcmPowspec *powspec, const gdouble kmin);
void ncm_powspec_require_kmax (NcmPowspec *powspec, const gdouble kmax);

gdouble ncm_powspec_get_zi (NcmPowspec *powspec);
gdouble ncm_powspec_get_zf (NcmPowspec *powspec);

gdouble ncm_powspec_get_kmin (NcmPowspec *powspec);
gdouble ncm_powspec_get_kmax (NcmPowspec *powspec);

gdouble ncm_powspec_get_reltol_spline (NcmPowspec *powspec);

void ncm_powspec_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk);

void ncm_powspec_prepare (NcmPowspec *powspec, NcmModel *model);
void ncm_powspec_prepare_if_needed (NcmPowspec *powspec, NcmModel *model);
gdouble ncm_powspec_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k);
void ncm_powspec_eval_vec (NcmPowspec *powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk);
NcmSpline2d *ncm_powspec_get_spline_2d (NcmPowspec *powspec, NcmModel *model);
NcmModelCtrl *ncm_powspec_peek_model_ctrl (NcmPowspec *powspec);

gdouble ncm_powspec_var_tophat_R (NcmPowspec *powspec, NcmModel *model, const gdouble reltol, const gdouble z, const gdouble R);
gdouble ncm_powspec_sigma_tophat_R (NcmPowspec *powspec, NcmModel *model, const gdouble reltol, const gdouble z, const gdouble R);

gdouble ncm_powspec_corr3d (NcmPowspec *powspec, NcmModel *model, const gdouble reltol, const gdouble z, const gdouble r);

gdouble ncm_powspec_sproj (NcmPowspec *powspec, NcmModel *model, const gdouble reltol, const gint ell, const gdouble z1, const gdouble z2, const gdouble xi1, const gdouble xi2);

G_END_DECLS

#endif /* _NCM_POWSPEC_H_ */

