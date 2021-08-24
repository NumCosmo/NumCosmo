/***************************************************************************
 *            nc_multiplicity_func_warren.h
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

#ifndef _NC_MULTIPLICITY_FUNC_WARREN_H_
#define _NC_MULTIPLICITY_FUNC_WARREN_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_WARREN             (nc_multiplicity_func_warren_get_type ())
#define NC_MULTIPLICITY_FUNC_WARREN(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_MULTIPLICITY_FUNC_WARREN, NcMultiplicityFuncWarren))
#define NC_MULTIPLICITY_FUNC_WARREN_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_MULTIPLICITY_FUNC_WARREN, NcMultiplicityFuncWarrenClass))
#define NC_IS_MULTIPLICITY_FUNC_WARREN(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_MULTIPLICITY_FUNC_WARREN))
#define NC_IS_MULTIPLICITY_FUNC_WARREN_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_MULTIPLICITY_FUNC_WARREN))
#define NC_MULTIPLICITY_FUNC_WARREN_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_MULTIPLICITY_FUNC_WARREN, NcMultiplicityFuncWarrenClass))

typedef struct _NcMultiplicityFuncWarrenClass NcMultiplicityFuncWarrenClass;
typedef struct _NcMultiplicityFuncWarren NcMultiplicityFuncWarren;
typedef struct _NcMultiplicityFuncWarrenPrivate NcMultiplicityFuncWarrenPrivate;

struct _NcMultiplicityFuncWarrenClass
{
  /*< private >*/
  NcMultiplicityFuncClass parent_class;
};

struct _NcMultiplicityFuncWarren
{
  /*< private >*/
  NcMultiplicityFunc parent_instance;
  NcMultiplicityFuncWarrenPrivate *priv;
};

GType nc_multiplicity_func_warren_get_type (void) G_GNUC_CONST;

NcMultiplicityFuncWarren *nc_multiplicity_func_warren_new (void);
NcMultiplicityFuncWarren *nc_multiplicity_func_warren_ref (NcMultiplicityFuncWarren *mw);

void nc_multiplicity_func_warren_free (NcMultiplicityFuncWarren *mw);
void nc_multiplicity_func_warren_clear (NcMultiplicityFuncWarren **mw);

G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_WARREN_H_ */
