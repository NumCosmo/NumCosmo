/***************************************************************************
 *            nc_scalefactor.h
 *
 *  Wed Nov 12 14:46:40 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NC_SCALEFACTOR_H_
#define _NC_SCALEFACTOR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_distance.h>

G_BEGIN_DECLS

#define NC_TYPE_SCALEFACTOR             (nc_scalefactor_get_type ())
#define NC_SCALEFACTOR(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_SCALEFACTOR, NcScalefactor))
#define NC_SCALEFACTOR_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_SCALEFACTOR, NcScalefactorClass))
#define NC_IS_SCALEFACTOR(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_SCALEFACTOR))
#define NC_IS_SCALEFACTOR_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_SCALEFACTOR))
#define NC_SCALEFACTOR_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_SCALEFACTOR, NcScalefactorClass))

typedef struct _NcScalefactorClass NcScalefactorClass;
typedef struct _NcScalefactor NcScalefactor;
typedef struct _NcScalefactorPrivate NcScalefactorPrivate;

struct _NcScalefactorClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcScalefactor
{
  /*< private >*/
  GObject parent_instance;
  NcScalefactorPrivate *priv;
};

GType nc_scalefactor_get_type (void) G_GNUC_CONST;

NcScalefactor *nc_scalefactor_new (const gdouble zf, NcDistance *dist);
NcScalefactor *nc_scalefactor_ref (NcScalefactor *a);
void nc_scalefactor_free (NcScalefactor *a);
void nc_scalefactor_clear (NcScalefactor **a);

void nc_scalefactor_prepare (NcScalefactor *a, NcHICosmo *cosmo);
void nc_scalefactor_prepare_if_needed (NcScalefactor *a, NcHICosmo *cosmo);

void nc_scalefactor_set_zf (NcScalefactor *a, const gdouble zf);
void nc_scalefactor_require_zf (NcScalefactor *a, const gdouble zf);

void nc_scalefactor_set_a0 (NcScalefactor *a, const gdouble a0);
void nc_scalefactor_set_reltol (NcScalefactor *a, const gdouble reltol);
void nc_scalefactor_set_abstol (NcScalefactor *a, const gdouble abstol);

void nc_scalefactor_set_a0_conformal_normal (NcScalefactor *a, gboolean enable);

gdouble nc_scalefactor_get_zf (NcScalefactor *a);
gdouble nc_scalefactor_get_a0 (NcScalefactor *a);
gdouble nc_scalefactor_get_reltol (NcScalefactor *a);
gdouble nc_scalefactor_get_abstol (NcScalefactor *a);

gdouble nc_scalefactor_eval_z_eta (NcScalefactor *a, const gdouble eta);
gdouble nc_scalefactor_eval_a_eta (NcScalefactor *a, const gdouble eta);

gdouble nc_scalefactor_eval_eta_z (NcScalefactor *a, const gdouble z);
gdouble nc_scalefactor_eval_eta_x (NcScalefactor *a, const gdouble x);

gdouble nc_scalefactor_eval_t_eta (NcScalefactor *a, const gdouble eta);
gdouble nc_scalefactor_eval_eta_t (NcScalefactor *a, const gdouble t);

#define NC_SCALEFACTOR_DEFAULT_ZF (1.0e14)
#define NC_SCALEFACTOR_DEFAULT_A0 (1.0)
#define NC_SCALEFACTOR_DEFAULT_RELTOL (1.0e-13)
#define NC_SCALEFACTOR_DEFAULT_ABSTOL (0.0)
#define NC_SCALEFACTOR_OMEGA_K_ZERO (1.0e-14)
#define NC_SCALEFACTOR_MIN_ETA_STEP (1.0e-11)

G_END_DECLS

#endif /* _NC_SCALEFACTOR_H_ */
