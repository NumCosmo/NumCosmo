/***************************************************************************
 *            ncm_fit_esmcmc_walker_newton.h
 *
 *  Sat October 27 13:08:34 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker_newton.h
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

#ifndef _NCM_FIT_ESMCMC_WALKER_NEWTON_H_
#define _NCM_FIT_ESMCMC_WALKER_NEWTON_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit_esmcmc_walker.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_mset.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_ESMCMC_WALKER_NEWTON             (ncm_fit_esmcmc_walker_newton_get_type ())
#define NCM_FIT_ESMCMC_WALKER_NEWTON(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT_ESMCMC_WALKER_NEWTON, NcmFitESMCMCWalkerNewton))
#define NCM_FIT_ESMCMC_WALKER_NEWTON_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT_ESMCMC_WALKER_NEWTON, NcmFitESMCMCWalkerNewtonClass))
#define NCM_IS_FIT_ESMCMC_WALKER_NEWTON(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT_ESMCMC_WALKER_NEWTON))
#define NCM_IS_FIT_ESMCMC_WALKER_NEWTON_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT_ESMCMC_WALKER_NEWTON))
#define NCM_FIT_ESMCMC_WALKER_NEWTON_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT_ESMCMC_WALKER_NEWTON, NcmFitESMCMCWalkerNewtonClass))

typedef struct _NcmFitESMCMCWalkerNewtonClass NcmFitESMCMCWalkerNewtonClass;
typedef struct _NcmFitESMCMCWalkerNewton NcmFitESMCMCWalkerNewton;
typedef struct _NcmFitESMCMCWalkerNewtonPrivate NcmFitESMCMCWalkerNewtonPrivate;

struct _NcmFitESMCMCWalkerNewtonClass
{
  /*< private >*/
  NcmFitESMCMCWalkerClass parent_class;
};

struct _NcmFitESMCMCWalkerNewton
{
  /*< private >*/
  NcmFitESMCMCWalker parent_instance;
  NcmFitESMCMCWalkerNewtonPrivate *priv;
};

GType ncm_fit_esmcmc_walker_newton_get_type (void) G_GNUC_CONST;

NcmFitESMCMCWalkerNewton *ncm_fit_esmcmc_walker_newton_new (guint nwalkers, guint nparams);

void ncm_fit_esmcmc_walker_newton_set_G (NcmFitESMCMCWalkerNewton *newton, const gdouble G);
gdouble ncm_fit_esmcmc_walker_newton_get_G (NcmFitESMCMCWalkerNewton *newton);

void ncm_fit_esmcmc_walker_newton_set_box (NcmFitESMCMCWalkerNewton *newton, guint n, const gdouble lb, const gdouble ub);
void ncm_fit_esmcmc_walker_newton_set_box_mset (NcmFitESMCMCWalkerNewton *newton, NcmMSet *mset);

G_END_DECLS

#endif /* _NCM_FIT_ESMCMC_WALKER_NEWTON_H_ */
