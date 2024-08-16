/***************************************************************************
 *            nc_multiplicity_func_despali.h
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

#ifndef _NC_MULTIPLICITY_FUNC_DESPALI_H_
#define _NC_MULTIPLICITY_FUNC_DESPALI_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_DESPALI             (nc_multiplicity_func_despali_get_type ())
#define NC_MULTIPLICITY_FUNC_DESPALI(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_MULTIPLICITY_FUNC_DESPALI, NcMultiplicityFuncDespali))
#define NC_MULTIPLICITY_FUNC_DESPALI_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_MULTIPLICITY_FUNC_DESPALI, NcMultiplicityFuncDespaliClass))
#define NC_IS_MULTIPLICITY_FUNC_DESPALI(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_MULTIPLICITY_FUNC_DESPALI))
#define NC_IS_MULTIPLICITY_FUNC_DESPALI_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_MULTIPLICITY_FUNC_DESPALI))
#define NC_MULTIPLICITY_FUNC_DESPALI_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_MULTIPLICITY_FUNC_DESPALI, NcMultiplicityFuncDespaliClass))

typedef struct _NcMultiplicityFuncDespaliClass NcMultiplicityFuncDespaliClass;
typedef struct _NcMultiplicityFuncDespali NcMultiplicityFuncDespali;
typedef struct _NcMultiplicityFuncDespaliPrivate NcMultiplicityFuncDespaliPrivate;

struct _NcMultiplicityFuncDespaliClass
{
  /*< private >*/
  NcMultiplicityFuncClass parent_class;
};

struct _NcMultiplicityFuncDespali
{
  /*< private >*/
  NcMultiplicityFunc parent_instance;
  NcMultiplicityFuncDespaliPrivate *priv;
};

GType nc_multiplicity_func_despali_get_type (void) G_GNUC_CONST;

NcMultiplicityFuncDespali *nc_multiplicity_func_despali_new (void);
NcMultiplicityFuncDespali *nc_multiplicity_func_despali_new_full (NcMultiplicityFuncMassDef mdef, gdouble Delta);
NcMultiplicityFuncDespali *nc_multiplicity_func_despali_ref (NcMultiplicityFuncDespali *md);

void nc_multiplicity_func_despali_free (NcMultiplicityFuncDespali *md);
void nc_multiplicity_func_despali_clear (NcMultiplicityFuncDespali **md);

gdouble nc_multiplicity_func_despali_delta_c (NcMultiplicityFuncDespali *md , NcHICosmo *cosmo ,gdouble z);

gdouble nc_multiplicity_func_despali_delta_vir (NcMultiplicityFuncDespali *md , NcHICosmo *cosmo ,gdouble z);
void nc_multiplicity_func_despali_set_eo (NcMultiplicityFuncDespali *md, gboolean on);
gboolean nc_multiplicity_func_despali_get_eo (NcMultiplicityFuncDespali *md);

G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_DESPALI_H_ */
