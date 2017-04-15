/***************************************************************************
 *            nc_multiplicity_func_crocce.h
 *
 *  Wed Feb 15 13:36:09 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2017 <pennalima@gmail.com>
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

#ifndef _NC_MULTIPLICITY_FUNC_CROCCE_H_
#define _NC_MULTIPLICITY_FUNC_CROCCE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_CROCCE             (nc_multiplicity_func_crocce_get_type ())
#define NC_MULTIPLICITY_FUNC_CROCCE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_MULTIPLICITY_FUNC_CROCCE, NcMultiplicityFuncCrocce))
#define NC_MULTIPLICITY_FUNC_CROCCE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_MULTIPLICITY_FUNC_CROCCE, NcMultiplicityFuncCrocceClass))
#define NC_IS_MULTIPLICITY_FUNC_CROCCE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_MULTIPLICITY_FUNC_CROCCE))
#define NC_IS_MULTIPLICITY_FUNC_CROCCE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_MULTIPLICITY_FUNC_CROCCE))
#define NC_MULTIPLICITY_FUNC_CROCCE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_MULTIPLICITY_FUNC_CROCCE, NcMultiplicityFuncCrocceClass))

typedef struct _NcMultiplicityFuncCrocceClass NcMultiplicityFuncCrocceClass;
typedef struct _NcMultiplicityFuncCrocce NcMultiplicityFuncCrocce;

struct _NcMultiplicityFuncCrocceClass
{
  /*< private >*/
  NcMultiplicityFuncClass parent_class;
};

struct _NcMultiplicityFuncCrocce
{
  /*< private >*/
  NcMultiplicityFunc parent_instance;
  gdouble A0;
  gdouble a0;
  gdouble b0;
  gdouble c0;
};

GType nc_multiplicity_func_crocce_get_type (void) G_GNUC_CONST;

NcMultiplicityFunc *nc_multiplicity_func_crocce_new (gdouble A0, gdouble a0, gdouble b0, gdouble c0);
void nc_multiplicity_func_crocce_set_A0 (NcMultiplicityFuncCrocce *mulf_crocce, gdouble A0);
gdouble nc_multiplicity_func_crocce_get_A0 (const NcMultiplicityFuncCrocce *mulf_crocce);
void nc_multiplicity_func_crocce_set_a0 (NcMultiplicityFuncCrocce *mulf_crocce, gdouble a0);
gdouble nc_multiplicity_func_crocce_get_a0 (const NcMultiplicityFuncCrocce *mulf_crocce);
void nc_multiplicity_func_crocce_set_b0 (NcMultiplicityFuncCrocce *mulf_crocce, gdouble b0);
gdouble nc_multiplicity_func_crocce_get_b0 (const NcMultiplicityFuncCrocce *mulf_crocce);
void nc_multiplicity_func_crocce_set_c0 (NcMultiplicityFuncCrocce *mulf_crocce, gdouble c0);
gdouble nc_multiplicity_func_crocce_get_c0 (const NcMultiplicityFuncCrocce *mulf_crocce);

G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_CROCCE_H_ */
