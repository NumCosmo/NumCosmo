/***************************************************************************
 *            nc_multiplicity_func_tinker_mean.h
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_MULTIPLICITY_FUNC_TINKER_MEAN_H_
#define _NC_MULTIPLICITY_FUNC_TINKER_MEAN_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN             (nc_multiplicity_func_tinker_mean_get_type ())
#define NC_MULTIPLICITY_FUNC_TINKER_MEAN(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN, NcMultiplicityFuncTinkerMean))
#define NC_MULTIPLICITY_FUNC_TINKER_MEAN_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN, NcMultiplicityFuncTinkerMeanClass))
#define NC_IS_MULTIPLICITY_FUNC_TINKER_MEAN(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN))
#define NC_IS_MULTIPLICITY_FUNC_TINKER_MEAN_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN))
#define NC_MULTIPLICITY_FUNC_TINKER_MEAN_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN, NcMultiplicityFuncTinkerMeanClass))

typedef struct _NcMultiplicityFuncTinkerMeanClass NcMultiplicityFuncTinkerMeanClass;
typedef struct _NcMultiplicityFuncTinkerMean NcMultiplicityFuncTinkerMean;

struct _NcMultiplicityFuncTinkerMeanClass
{
  /*< private >*/
  NcMultiplicityFuncClass parent_class;
};

struct _NcMultiplicityFuncTinkerMean
{
  /*< private >*/
  NcMultiplicityFunc parent_instance;
  gdouble Delta;
};

GType nc_multiplicity_func_tinker_mean_get_type (void) G_GNUC_CONST;

NcMultiplicityFunc *nc_multiplicity_func_tinker_mean_new (gdouble Delta);
void nc_multiplicity_func_tinker_mean_set_Delta (NcMultiplicityFuncTinkerMean *mulf_tinker_mean, gdouble Delta);
gdouble nc_multiplicity_func_tinker_mean_get_Delta (const NcMultiplicityFuncTinkerMean *mulf_tinker_mean);

G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_TINKER_MEAN_H_ */
