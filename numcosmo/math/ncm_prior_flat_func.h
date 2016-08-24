/***************************************************************************
 *            ncm_prior_flat_func.h
 *
 *  Wed August 03 16:26:48 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_prior_flat_func.h
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

#ifndef _NCM_PRIOR_FLAT_FUNC_H_
#define _NCM_PRIOR_FLAT_FUNC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_prior_flat.h>

G_BEGIN_DECLS

#define NCM_TYPE_PRIOR_FLAT_FUNC             (ncm_prior_flat_func_get_type ())
#define NCM_PRIOR_FLAT_FUNC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_PRIOR_FLAT_FUNC, NcmPriorFlatFunc))
#define NCM_PRIOR_FLAT_FUNC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_PRIOR_FLAT_FUNC, NcmPriorFlatFuncClass))
#define NCM_IS_PRIOR_FLAT_FUNC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_PRIOR_FLAT_FUNC))
#define NCM_IS_PRIOR_FLAT_FUNC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_PRIOR_FLAT_FUNC))
#define NCM_PRIOR_FLAT_FUNC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_PRIOR_FLAT_FUNC, NcmPriorFlatFuncClass))

typedef struct _NcmPriorFlatFuncClass NcmPriorFlatFuncClass;
typedef struct _NcmPriorFlatFunc NcmPriorFlatFunc;

struct _NcmPriorFlatFuncClass
{
  /*< private >*/
  NcmPriorFlatClass parent_class;
};

struct _NcmPriorFlatFunc
{
  /*< private >*/
  NcmPriorFlat parent_instance;
  NcmMSetFunc *mean_func;
};

GType ncm_prior_flat_func_get_type (void) G_GNUC_CONST;

NcmPriorFlatFunc *ncm_prior_flat_func_new (NcmMSetFunc *mean_func, gdouble x_low, gdouble x_upp, gdouble scale, gdouble variable);
NcmPriorFlatFunc *ncm_prior_flat_func_ref (NcmPriorFlatFunc *pff);

void ncm_prior_flat_func_free (NcmPriorFlatFunc *pff);
void ncm_prior_flat_func_clear (NcmPriorFlatFunc **pff);

G_END_DECLS

#endif /* _NCM_PRIOR_FLAT_FUNC_H_ */
