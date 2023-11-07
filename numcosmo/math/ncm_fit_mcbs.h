/***************************************************************************
 *            ncm_fit_mcbs.h
 *
 *  Tue February 11 13:54:23 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_fit_mcbs.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_FIT_MCBS_H_
#define _NCM_FIT_MCBS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>
#include <numcosmo/math/ncm_fit_mc.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_MCBS (ncm_fit_mcbs_get_type ())

G_DECLARE_FINAL_TYPE (NcmFitMCBS, ncm_fit_mcbs, NCM, FIT_MCBS, GObject)

NcmFitMCBS *ncm_fit_mcbs_new (NcmFit * fit);
void ncm_fit_mcbs_free (NcmFitMCBS *mcbs);
void ncm_fit_mcbs_clear (NcmFitMCBS **mcbs);

void ncm_fit_mcbs_set_filename (NcmFitMCBS *mcbs, const gchar *filename);
void ncm_fit_mcbs_set_rng (NcmFitMCBS *mcbs, NcmRNG *rng);
void ncm_fit_mcbs_run (NcmFitMCBS *mcbs, NcmMSet *fiduc, guint ni, guint nf, guint nbstraps, NcmFitMCResampleType rtype, NcmFitRunMsgs mtype, guint bsmt);

NcmMSetCatalog *ncm_fit_mcbs_get_catalog (NcmFitMCBS *mcbs);

G_END_DECLS

#endif /* _NCM_FIT_MCBS_H_ */

