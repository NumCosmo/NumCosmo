/***************************************************************************
 *            nc_multiplicity_func_tinker_crit.h
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

#ifndef _NC_MULTIPLICITY_FUNC_TINKER_CRIT_H_
#define _NC_MULTIPLICITY_FUNC_TINKER_CRIT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_TINKER_CRIT             (nc_multiplicity_func_tinker_crit_get_type ())
#define NC_MULTIPLICITY_FUNC_TINKER_CRIT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_MULTIPLICITY_FUNC_TINKER_CRIT, NcMultiplicityFuncTinkerCrit))
#define NC_MULTIPLICITY_FUNC_TINKER_CRIT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_MULTIPLICITY_FUNC_TINKER_CRIT, NcMultiplicityFuncTinkerCritClass))
#define NC_IS_MULTIPLICITY_FUNC_TINKER_CRIT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_MULTIPLICITY_FUNC_TINKER_CRIT))
#define NC_IS_MULTIPLICITY_FUNC_TINKER_CRIT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_MULTIPLICITY_FUNC_TINKER_CRIT))
#define NC_MULTIPLICITY_FUNC_TINKER_CRIT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_MULTIPLICITY_FUNC_TINKER_CRIT, NcMultiplicityFuncTinkerCritClass))

typedef struct _NcMultiplicityFuncTinkerCritClass NcMultiplicityFuncTinkerCritClass;
typedef struct _NcMultiplicityFuncTinkerCrit NcMultiplicityFuncTinkerCrit;

struct _NcMultiplicityFuncTinkerCrit
{
  /*< private >*/
  NcMultiplicityFunc parent_instance;
  gdouble Delta;
};

struct _NcMultiplicityFuncTinkerCritClass
{
  /*< private >*/
  NcMultiplicityFuncClass parent_class;
};

GType nc_multiplicity_func_tinker_crit_get_type (void) G_GNUC_CONST;

NcMultiplicityFunc *nc_multiplicity_func_tinker_crit_new (gdouble Delta);
void nc_multiplicity_func_tinker_crit_set_Delta (NcMultiplicityFuncTinkerCrit *mulf_tinker_crit, gdouble Delta);
gdouble nc_multiplicity_func_tinker_crit_get_Delta (const NcMultiplicityFuncTinkerCrit *mulf_tinker_crit);

G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_TINKER_CRIT_H_ */
