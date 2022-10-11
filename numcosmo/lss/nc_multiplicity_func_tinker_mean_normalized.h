/***************************************************************************
 *            nc_multiplicity_func_tinker_mean_normalized.h
 *
 *  Thu July 03 12:31:21 2014
 *  Copyright  2014  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_multiplicity_func_tinker_mean_normalized.h
 * Copyright (C) 2014 Mariana Penna Lima <pennalima@gmail.com>
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

#ifndef _NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED_H_
#define _NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED             (nc_multiplicity_func_tinker_mean_normalized_get_type ())
#define NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED, NcMultiplicityFuncTinkerMeanNormalized))
#define NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED, NcMultiplicityFuncTinkerMeanNormalizedClass))
#define NC_IS_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED))
#define NC_IS_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED))
#define NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED, NcMultiplicityFuncTinkerMeanNormalizedClass))

typedef struct _NcMultiplicityFuncTinkerMeanNormalizedClass NcMultiplicityFuncTinkerMeanNormalizedClass;
typedef struct _NcMultiplicityFuncTinkerMeanNormalized NcMultiplicityFuncTinkerMeanNormalized;
typedef struct _NcMultiplicityFuncTinkerMeanNormalizedPrivate NcMultiplicityFuncTinkerMeanNormalizedPrivate;

struct _NcMultiplicityFuncTinkerMeanNormalizedClass
{
  /*< private >*/
  NcMultiplicityFuncClass parent_class;
};

struct _NcMultiplicityFuncTinkerMeanNormalized
{
  /*< private >*/
  NcMultiplicityFunc parent_instance;
  NcMultiplicityFuncTinkerMeanNormalizedPrivate *priv;
};

GType nc_multiplicity_func_tinker_mean_normalized_get_type (void) G_GNUC_CONST;

NcMultiplicityFuncTinkerMeanNormalized *nc_multiplicity_func_tinker_mean_normalized_new (void);
NcMultiplicityFuncTinkerMeanNormalized *nc_multiplicity_func_tinker_mean_normalized_ref (NcMultiplicityFuncTinkerMeanNormalized *mt10);

void nc_multiplicity_func_tinker_mean_normalized_free (NcMultiplicityFuncTinkerMeanNormalized *mt10);
void nc_multiplicity_func_tinker_mean_normalized_clear (NcMultiplicityFuncTinkerMeanNormalized **mt10);


G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED_H_ */

