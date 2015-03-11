/***************************************************************************
 *            ncm_fit_mcbs.h
 *
 *  Tue February 11 13:54:23 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fit_mcbs.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#define NCM_TYPE_FIT_MCBS             (ncm_fit_mcbs_get_type ())
#define NCM_FIT_MCBS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT_MCBS, NcmFitMCBS))
#define NCM_FIT_MCBS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT_MCBS, NcmFitMCBSClass))
#define NCM_IS_FIT_MCBS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT_MCBS))
#define NCM_IS_FIT_MCBS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT_MCBS))
#define NCM_FIT_MCBS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT_MCBS, NcmFitMCBSClass))

typedef struct _NcmFitMCBSClass NcmFitMCBSClass;
typedef struct _NcmFitMCBS NcmFitMCBS;

struct _NcmFitMCBS
{
  /*< private >*/
  GObject parent_instance;
  NcmFit *fit;
  NcmFitMC *mc_resample;
  NcmFitMC *mc_bstrap;
  NcmMSetCatalog *mcat;
  gchar *base_name;
};

struct _NcmFitMCBSClass
{
  /*< private >*/
  GObjectClass parent_class;
};

GType ncm_fit_mcbs_get_type (void) G_GNUC_CONST;

NcmFitMCBS *ncm_fit_mcbs_new (NcmFit *fit);
void ncm_fit_mcbs_free (NcmFitMCBS *mcbs);
void ncm_fit_mcbs_clear (NcmFitMCBS **mcbs);

void ncm_fit_mcbs_set_filename (NcmFitMCBS *mcbs, const gchar *filename);
void ncm_fit_mcbs_set_rng (NcmFitMCBS *mcbs, NcmRNG *rng);
void ncm_fit_mcbs_run (NcmFitMCBS *mcbs, NcmMSet *fiduc, guint ni, guint nf, guint nbstraps, NcmFitMCResampleType rtype, NcmFitRunMsgs mtype, guint bsmt);

NcmMSetCatalog *ncm_fit_mcbs_get_catalog (NcmFitMCBS *mcbs);

G_END_DECLS

#endif /* _NCM_FIT_MCBS_H_ */
