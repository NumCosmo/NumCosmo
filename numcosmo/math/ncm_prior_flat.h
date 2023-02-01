/***************************************************************************
 *            ncm_prior_flat.h
 *
 *  Wed August 03 16:58:26 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_prior_flat.h
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

#ifndef _NCM_PRIOR_FLAT_H_
#define _NCM_PRIOR_FLAT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_prior.h>

G_BEGIN_DECLS

#define NCM_TYPE_PRIOR_FLAT             (ncm_prior_flat_get_type ())
#define NCM_PRIOR_FLAT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_PRIOR_FLAT, NcmPriorFlat))
#define NCM_PRIOR_FLAT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_PRIOR_FLAT, NcmPriorFlatClass))
#define NCM_IS_PRIOR_FLAT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_PRIOR_FLAT))
#define NCM_IS_PRIOR_FLAT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_PRIOR_FLAT))
#define NCM_PRIOR_FLAT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_PRIOR_FLAT, NcmPriorFlatClass))

typedef struct _NcmPriorFlatClass NcmPriorFlatClass;
typedef struct _NcmPriorFlat NcmPriorFlat;

typedef gdouble (*NcmPriorFlatMean) (NcmPriorFlat *pf, NcmMSet *mset);

struct _NcmPriorFlatClass
{
  /*< private >*/
  NcmPriorClass parent_class;
  NcmPriorFlatMean mean;  
};

struct _NcmPriorFlat
{
  /*< private >*/
  NcmPrior parent_instance;
  gdouble x_low;
  gdouble x_upp;
  gdouble s;
  gdouble var;
};

GType ncm_prior_flat_get_type (void) G_GNUC_CONST;

NcmPriorFlat *ncm_prior_flat_ref (NcmPriorFlat *pf);
void ncm_prior_flat_free (NcmPriorFlat *pf);
void ncm_prior_flat_clear (NcmPriorFlat **pf);

G_END_DECLS

#endif /* _NCM_PRIOR_FLAT_H_ */
