/***************************************************************************
 *            nc_hireion.h
 *
 *  Tue December 08 14:51:23 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hireion.h
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_HIREION_H_
#define _NC_HIREION_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>

G_BEGIN_DECLS

#define NC_TYPE_HIREION             (nc_hireion_get_type ())
#define NC_HIREION(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIREION, NcHIReion))
#define NC_HIREION_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIREION, NcHIReionClass))
#define NC_IS_HIREION(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIREION))
#define NC_IS_HIREION_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIREION))
#define NC_HIREION_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIREION, NcHIReionClass))

typedef struct _NcHIReionClass NcHIReionClass;
typedef struct _NcHIReion NcHIReion;

struct _NcHIReionClass
{
  /*< private >*/
  NcmModelClass parent_class;
  gdouble (*get_init_x) (NcHIReion *reion, NcHICosmo *cosmo);
  gdouble (*get_Xe) (NcHIReion *reion, NcHICosmo *cosmo, const gdouble lambda, const gdouble Xe_recomb);
  gdouble (*get_tau) (NcHIReion *reion, NcHICosmo *cosmo);
};

struct _NcHIReion
{
  /*< private >*/
  NcmModel parent_instance;
  gdouble prec;
};

GType nc_hireion_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_hireion);

NcHIReion *nc_hireion_new_from_name (GType parent_type, gchar *reion_name);
NcHIReion *nc_hireion_ref (NcHIReion *reion);
void nc_hireion_free (NcHIReion *reion);
void nc_hireion_clear (NcHIReion **reion);

gdouble nc_hireion_get_init_x (NcHIReion *reion, NcHICosmo *cosmo);
gdouble nc_hireion_get_Xe (NcHIReion *reion, NcHICosmo *cosmo, const gdouble lambda, const gdouble Xe_recomb);
gdouble nc_hireion_get_tau (NcHIReion *reion, NcHICosmo *cosmo);

#define NC_HIREION_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_HIREION_H_ */
