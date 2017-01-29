/***************************************************************************
 *            nc_powspec_mnl_halofit.h
 *
 *  Thu March 17 14:57:27 2016
 *  Copyright  2016  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 ****************************************************************************/
/*
 * nc_powspec_mnl_halofit.h
 * Copyright (C) 2016 Cyrille Doux <cdoux@apc.in2p3.fr>
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

#ifndef _NC_POWSPEC_MNL_HALOFIT_H_
#define _NC_POWSPEC_MNL_HALOFIT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_powspec_filter.h>
#include <numcosmo/nc_powspec_mnl.h>
#include <numcosmo/nc_powspec_ml.h>


G_BEGIN_DECLS

#define NC_TYPE_POWSPEC_MNL_HALOFIT (nc_powspec_mnl_halofit_get_type ())
#define NC_POWSPEC_MNL_HALOFIT(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_POWSPEC_MNL_HALOFIT, NcPowspecMNLHaloFit))
#define NC_POWSPEC_MNL_HALOFIT_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_POWSPEC_MNL_HALOFIT, NcPowspecMNLHaloFitClass))
#define NC_IS_POWSPEC_MNL_HALOFIT(obj) (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_POWSPEC_MNL_HALOFIT))
#define NC_IS_POWSPEC_MNL_HALOFIT_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_POWSPEC_MNL_HALOFIT))
#define NC_POWSPEC_MNL_HALOFIT_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_POWSPEC_MNL_HALOFIT, NcPowspecMNLHaloFitClass))

typedef struct _NcPowspecMNLHaloFitClass NcPowspecMNLHaloFitClass;
typedef struct _NcPowspecMNLHaloFit NcPowspecMNLHaloFit;
typedef struct _NcPowspecMNLHaloFitPrivate NcPowspecMNLHaloFitPrivate;

struct _NcPowspecMNLHaloFitClass
{
  /*< private > */
  NcPowspecMNLClass parent_class;
};

struct _NcPowspecMNLHaloFit
{
  /*< private > */
  NcPowspecMNL parent_instance;
  NcPowspecML* psml;
  gdouble zmaxnl;
  gdouble znl;
  gdouble reltol;
  NcmSpline *Rsigma;
  NcmSpline *neff;
  NcmSpline *Cur;
  NcmPowspecFilter *psml_gauss;
  NcPowspecMNLHaloFitPrivate *priv;
};

GType nc_powspec_mnl_halofit_get_type (void) G_GNUC_CONST;

NcPowspecMNLHaloFit *nc_powspec_mnl_halofit_new (NcPowspecML *psml, gdouble zmaxnl, gdouble reltol);

void nc_powspec_mnl_halofit_set_kbounds_from_ml (NcPowspecMNLHaloFit *pshf);
void nc_powspec_mnl_halofit_pkequal (NcPowspecMNLHaloFit *pshf, gboolean on);

#define NC_POWSPEC_MNL_HALOFIT_F1aPOW   (-0.0732)
#define NC_POWSPEC_MNL_HALOFIT_F2aPOW   (-0.1423)
#define NC_POWSPEC_MNL_HALOFIT_F3aPOW   (0.0725)
#define NC_POWSPEC_MNL_HALOFIT_F1bPOW   (-0.0307)
#define NC_POWSPEC_MNL_HALOFIT_F2bPOW   (-0.0585)
#define NC_POWSPEC_MNL_HALOFIT_F3bPOW   (0.0743)
#define NC_POWSPEC_MNL_HALOFIT_LOGRMIN (-35.0)

G_END_DECLS

#endif /* _NC_POWSPEC_MNL_HALOFIT_H_ */
