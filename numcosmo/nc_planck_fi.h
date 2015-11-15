/***************************************************************************
 *            nc_planck_fi.h
 *
 *  Thu October 22 15:47:48 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_planck_fi.h
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

#ifndef _NC_PLANCK_FI_H_
#define _NC_PLANCK_FI_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_mset.h>

G_BEGIN_DECLS

#define NC_TYPE_PLANCK_FI             (nc_planck_fi_get_type ())
#define NC_PLANCK_FI(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_PLANCK_FI, NcPlanckFI))
#define NC_PLANCK_FI_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_PLANCK_FI, NcPlanckFIClass))
#define NC_IS_PLANCK_FI(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_PLANCK_FI))
#define NC_IS_PLANCK_FI_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_PLANCK_FI))
#define NC_PLANCK_FI_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_PLANCK_FI, NcPlanckFIClass))

typedef struct _NcPlanckFIClass NcPlanckFIClass;
typedef struct _NcPlanckFI NcPlanckFI;

struct _NcPlanckFIClass
{
  /*< private >*/
  NcmModelClass parent_class;
};

struct _NcPlanckFI
{
  /*< private >*/
  NcmModel parent_instance;
  guint version;
};

GType nc_planck_fi_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_planck_fi);

NcPlanckFI *nc_planck_fi_new_from_name (gchar *pfi_name);
NcPlanckFI *nc_planck_fi_ref (NcPlanckFI *pfi);
void nc_planck_fi_free (NcPlanckFI *pfi);
void nc_planck_fi_clear (NcPlanckFI **pfi);

void nc_planck_fi_log_all_models (void);

#define NC_PLANCK_FI_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_PLANCK_FI_H_ */
