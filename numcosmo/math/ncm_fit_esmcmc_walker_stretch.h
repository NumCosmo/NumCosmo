/***************************************************************************
 *            ncm_fit_esmcmc_walker_stretch.h
 *
 *  Wed March 16 15:53:15 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker_stretch.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

G_BEGIN_DECLS

#define NCM_TYPE_FIT_ESMCMC_WALKER_STRETCH             (ncm_fit_esmcmc_walker_stretch_get_type ())
#define NCM_FIT_ESMCMC_WALKER_STRETCH(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT_ESMCMC_WALKER_STRETCH, NcmFitESMCMCWalkerStretch))
#define NCM_FIT_ESMCMC_WALKER_STRETCH_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT_ESMCMC_WALKER_STRETCH, NcmFitESMCMCWalkerStretchClass))
#define NCM_IS_FIT_ESMCMC_WALKER_STRETCH(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT_ESMCMC_WALKER_STRETCH))
#define NCM_IS_FIT_ESMCMC_WALKER_STRETCH_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT_ESMCMC_WALKER_STRETCH))
#define NCM_FIT_ESMCMC_WALKER_STRETCH_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT_ESMCMC_WALKER_STRETCH, NcmFitESMCMCWalkerStretchClass))

typedef struct _NcmFitESMCMCWalkerStretchClass NcmFitESMCMCWalkerStretchClass;
typedef struct _NcmFitESMCMCWalkerStretch NcmFitESMCMCWalkerStretch;

struct _NcmFitESMCMCWalkerStretchClass
{
  /*< private >*/
  NcmFitESMCMCWalkerClass parent_class;
};

struct _NcmFitESMCMCWalkerStretch
{
  /*< private >*/
  NcmFitESMCMCWalker parent_instance;
  guint size;
  guint size_2;
  gdouble a;
  NcmVector *z;
  GArray *indices;
};

GType ncm_fit_esmcmc_walker_stretch_get_type (void) G_GNUC_CONST;

NcmFitESMCMCWalkerStretch *ncm_fit_esmcmc_walker_stretch_new (guint nwalkers);

void ncm_fit_esmcmc_walker_stretch_set_scale (NcmFitESMCMCWalkerStretch *stretch, const gdouble a);
gdouble ncm_fit_esmcmc_walker_stretch_get_scale (NcmFitESMCMCWalkerStretch *stretch);

G_END_DECLS

#endif /* _NCM_FIT_ESMCMC_WALKER_STRETCH_H_ */
