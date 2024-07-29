/***************************************************************************
 *            ncm_fit_esmcmc_walker_walk.h
 *
 *  Tue March 29 10:42:03 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker_walk.h
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

#ifndef _NCM_FIT_ESMCMC_WALKER_WALK_H_
#define _NCM_FIT_ESMCMC_WALKER_WALK_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit_esmcmc_walker.h>
#include <numcosmo/math/ncm_matrix.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_ESMCMC_WALKER_WALK (ncm_fit_esmcmc_walker_walk_get_type ())

G_DECLARE_FINAL_TYPE (NcmFitESMCMCWalkerWalk, ncm_fit_esmcmc_walker_walk, NCM, FIT_ESMCMC_WALKER_WALK, NcmFitESMCMCWalker)

NcmFitESMCMCWalkerWalk *ncm_fit_esmcmc_walker_walk_new (guint nwalkers);

void ncm_fit_esmcmc_walker_walk_set_scale (NcmFitESMCMCWalkerWalk *walk, const gdouble a);
gdouble ncm_fit_esmcmc_walker_walk_get_scale (NcmFitESMCMCWalkerWalk *walk);

G_END_DECLS

#endif /* _NCM_FIT_ESMCMC_WALKER_WALK_H_ */

