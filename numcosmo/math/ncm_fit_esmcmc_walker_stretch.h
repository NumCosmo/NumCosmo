/***************************************************************************
 *            ncm_fit_esmcmc_walker_stretch.h
 *
 *  Wed March 16 15:53:15 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker_stretch.h
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

#ifndef _NCM_FIT_ESMCMC_WALKER_STRETCH_H_
#define _NCM_FIT_ESMCMC_WALKER_STRETCH_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit_esmcmc_walker.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_mset.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_ESMCMC_WALKER_STRETCH (ncm_fit_esmcmc_walker_stretch_get_type ())

G_DECLARE_FINAL_TYPE (NcmFitESMCMCWalkerStretch, ncm_fit_esmcmc_walker_stretch, NCM, FIT_ESMCMC_WALKER_STRETCH, NcmFitESMCMCWalker)

NcmFitESMCMCWalkerStretch *ncm_fit_esmcmc_walker_stretch_new (guint nwalkers, guint nparams);

void ncm_fit_esmcmc_walker_stretch_set_scale (NcmFitESMCMCWalkerStretch *stretch, const gdouble a);
gdouble ncm_fit_esmcmc_walker_stretch_get_scale (NcmFitESMCMCWalkerStretch *stretch);

void ncm_fit_esmcmc_walker_stretch_set_box (NcmFitESMCMCWalkerStretch *stretch, guint n, const gdouble lb, const gdouble ub);
void ncm_fit_esmcmc_walker_stretch_set_box_mset (NcmFitESMCMCWalkerStretch *stretch, NcmMSet *mset);

void ncm_fit_esmcmc_walker_stretch_multi (NcmFitESMCMCWalkerStretch *stretch, gboolean multi);

G_END_DECLS

#endif /* _NCM_FIT_ESMCMC_WALKER_STRETCH_H_ */

