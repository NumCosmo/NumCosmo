/***************************************************************************
 *            ncm_fit_esmcmc_walker_aps.h
 *
 *  Sat October 27 13:08:34 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker_aps.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_FIT_ESMCMC_WALKER_APS_H_
#define _NCM_FIT_ESMCMC_WALKER_APS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit_esmcmc_walker.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_mset.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_ESMCMC_WALKER_APS             (ncm_fit_esmcmc_walker_aps_get_type ())
#define NCM_FIT_ESMCMC_WALKER_APS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT_ESMCMC_WALKER_APS, NcmFitESMCMCWalkerAPS))
#define NCM_FIT_ESMCMC_WALKER_APS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT_ESMCMC_WALKER_APS, NcmFitESMCMCWalkerAPSClass))
#define NCM_IS_FIT_ESMCMC_WALKER_APS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT_ESMCMC_WALKER_APS))
#define NCM_IS_FIT_ESMCMC_WALKER_APS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT_ESMCMC_WALKER_APS))
#define NCM_FIT_ESMCMC_WALKER_APS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT_ESMCMC_WALKER_APS, NcmFitESMCMCWalkerAPSClass))

typedef struct _NcmFitESMCMCWalkerAPSClass NcmFitESMCMCWalkerAPSClass;
typedef struct _NcmFitESMCMCWalkerAPS NcmFitESMCMCWalkerAPS;
typedef struct _NcmFitESMCMCWalkerAPSPrivate NcmFitESMCMCWalkerAPSPrivate;

struct _NcmFitESMCMCWalkerAPSClass
{
  /*< private >*/
  NcmFitESMCMCWalkerClass parent_class;
};

struct _NcmFitESMCMCWalkerAPS
{
  /*< private >*/
  NcmFitESMCMCWalker parent_instance;
  NcmFitESMCMCWalkerAPSPrivate *priv;
};

GType ncm_fit_esmcmc_walker_aps_get_type (void) G_GNUC_CONST;

NcmFitESMCMCWalkerAPS *ncm_fit_esmcmc_walker_aps_new (guint nwalkers, guint nparams);

void ncm_fit_esmcmc_walker_aps_set_G (NcmFitESMCMCWalkerAPS *aps, const gdouble G);
gdouble ncm_fit_esmcmc_walker_aps_get_G (NcmFitESMCMCWalkerAPS *aps);

void ncm_fit_esmcmc_walker_aps_set_box (NcmFitESMCMCWalkerAPS *aps, guint n, const gdouble lb, const gdouble ub);
void ncm_fit_esmcmc_walker_aps_set_box_mset (NcmFitESMCMCWalkerAPS *aps, NcmMSet *mset);

G_END_DECLS

#endif /* _NCM_FIT_ESMCMC_WALKER_APS_H_ */
