/***************************************************************************
 *            nc_multiplicity_func_watson.h
 *
 *  Sat Sep 11 17:45:13 2021
 *	Copyright  2021  Cinthia Nunes de Lima / Mariana Penna Lima
 *	<cinthia.n.lima@uel.br> / <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 Copyright (C) Cinthia Nunes de Lima <cinthia.n.lima@uel.br>, Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_MULTIPLICITY_FUNC_WATSON_H_
#define _NC_MULTIPLICITY_FUNC_WATSON_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_WATSON             (nc_multiplicity_func_watson_get_type ())
#define NC_MULTIPLICITY_FUNC_WATSON(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_MULTIPLICITY_FUNC_WATSON, NcMultiplicityFuncWatson))
#define NC_MULTIPLICITY_FUNC_WATSON_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_MULTIPLICITY_FUNC_WATSON, NcMultiplicityFuncWatsonClass))
#define NC_IS_MULTIPLICITY_FUNC_WATSON(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_MULTIPLICITY_FUNC_WATSON))
#define NC_IS_MULTIPLICITY_FUNC_WATSON_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_MULTIPLICITY_FUNC_WATSON))
#define NC_MULTIPLICITY_FUNC_WATSON_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_MULTIPLICITY_FUNC_WATSON, NcMultiplicityFuncWatsonClass))

typedef struct _NcMultiplicityFuncWatsonClass NcMultiplicityFuncWatsonClass;
typedef struct _NcMultiplicityFuncWatson NcMultiplicityFuncWatson;
typedef struct _NcMultiplicityFuncWatsonPrivate NcMultiplicityFuncWatsonPrivate;

struct _NcMultiplicityFuncWatsonClass
{
  /*< private >*/
  NcMultiplicityFuncClass parent_class;
};

struct _NcMultiplicityFuncWatson
{
  /*< private >*/
  NcMultiplicityFunc parent_instance;
  NcMultiplicityFuncWatsonPrivate *priv;
};

GType nc_multiplicity_func_watson_get_type (void) G_GNUC_CONST;

NcMultiplicityFuncWatson *nc_multiplicity_func_watson_new (void);
NcMultiplicityFuncWatson *nc_multiplicity_func_watson_ref (NcMultiplicityFuncWatson *mwat);

void nc_multiplicity_func_watson_free (NcMultiplicityFuncWatson *mwat);
void nc_multiplicity_func_watson_clear (NcMultiplicityFuncWatson **mwat);


G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_WATSON_H_ */
