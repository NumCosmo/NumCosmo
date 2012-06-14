/***************************************************************************
 *            nc_multiplicity_func_st.h
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

#ifndef _NC_MULTIPLICITY_FUNC_ST_H_
#define _NC_MULTIPLICITY_FUNC_ST_H_

#include <glib-object.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_ST             (nc_multiplicity_func_st_get_type ())
#define NC_MULTIPLICITY_FUNC_ST(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_MULTIPLICITY_FUNC_ST, NcMultiplicityFuncST))
#define NC_MULTIPLICITY_FUNC_ST_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_MULTIPLICITY_FUNC_ST, NcMultiplicityFuncSTClass))
#define NC_IS_MULTIPLICITY_FUNC_ST(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_MULTIPLICITY_FUNC_ST))
#define NC_IS_MULTIPLICITY_FUNC_ST_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_MULTIPLICITY_FUNC_ST))
#define NC_MULTIPLICITY_FUNC_ST_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_MULTIPLICITY_FUNC_ST, NcMultiplicityFuncSTClass))

typedef struct _NcMultiplicityFuncSTClass NcMultiplicityFuncSTClass;
typedef struct _NcMultiplicityFuncST NcMultiplicityFuncST;



struct _NcMultiplicityFuncSTClass
{
  /*< private >*/
  NcMultiplicityFuncClass parent_class;
};

struct _NcMultiplicityFuncST
{
  /*< private >*/
  NcMultiplicityFunc parent_instance;
  gdouble A; 
  gdouble b;
  gdouble p;
  gdouble delta_c;
};

GType nc_multiplicity_func_st_get_type (void) G_GNUC_CONST;

NcMultiplicityFunc *nc_multiplicity_func_st_new (gdouble A, gdouble b, gdouble p, gdouble delta_c);
void nc_multiplicity_func_st_set_A (NcMultiplicityFuncST *mulf_st, gdouble A);
gdouble nc_multiplicity_func_st_get_A (const NcMultiplicityFuncST *mulf_st);
void nc_multiplicity_func_st_set_b (NcMultiplicityFuncST *mulf_st, gdouble b);
gdouble nc_multiplicity_func_st_get_b (const NcMultiplicityFuncST *mulf_st);
void nc_multiplicity_func_st_set_p (NcMultiplicityFuncST *mulf_st, gdouble p);
gdouble nc_multiplicity_func_st_get_p (const NcMultiplicityFuncST *mulf_st);
void nc_multiplicity_func_st_set_delta_c (NcMultiplicityFuncST *mulf_st, gdouble delta_c);
gdouble nc_multiplicity_func_st_get_delta_c (const NcMultiplicityFuncST *mulf_st);

G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_ST_H_ */
