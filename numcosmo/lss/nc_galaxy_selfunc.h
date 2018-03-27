/***************************************************************************
 *            nc_galaxy_selfunc.h
 *
 *  Wed March 14 16:30:36 2018
 *  Copyright  2018 Fernando de Simoni 
 *  <fsimoni@id.uff.br>
 ****************************************************************************/
/*
 * nc_galaxy_selfunc.h
 * Copyright (C) 2018 Fernando de Simoni <fsimoni@id.uff.brr>
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

#ifndef _NC_GALAXY_SELFUNC_H_
#define _NC_GALAXY_SELFUNC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_obj_array.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SELFUNC             (nc_galaxy_selfunc_get_type ())
#define NC_GALAXY_SELFUNC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GALAXY_SELFUNC, NcGalaxySelfunc))
#define NC_GALAXY_SELFUNC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GALAXY_SELFUNC, NcGalaxySelfuncClass))
#define NC_IS_GALAXY_SELFUNC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GALAXY_SELFUNC))
#define NC_IS_GALAXY_SELFUNC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GALAXY_SELFUNC))
#define NC_GALAXY_SELFUNC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GALAXY_SELFUNC, NcGalaxySelfuncClass))

typedef struct _NcGalaxySelfuncClass NcGalaxySelfuncClass;
typedef struct _NcGalaxySelfunc NcGalaxySelfunc;
typedef struct _NcGalaxySelfuncPrivate NcGalaxySelfuncPrivate;

struct _NcGalaxySelfuncClass
{
  /*< private > */
  GObjectClass parent_class;
};

struct _NcGalaxySelfunc 
{
  /*< private >*/
  GObject parent_instance;
  NcGalaxySelfuncPrivate *priv;
};

GType nc_galaxy_selfunc_get_type (void) G_GNUC_CONST;

NcGalaxySelfunc *nc_galaxy_selfunc_new (const guint nshells);
NcGalaxySelfunc *nc_galaxy_selfunc_ref (NcGalaxySelfunc *gsf);

void nc_galaxy_selfunc_free (NcGalaxySelfunc *gsf);
void nc_galaxy_selfunc_clear (NcGalaxySelfunc **gsf);

void nc_galaxy_selfunc_set_nshells (NcGalaxySelfunc *gsf, const guint nshells);
guint nc_galaxy_selfunc_get_nshells (NcGalaxySelfunc *gsf);

void nc_galaxy_selfunc_set_shell_splines (NcGalaxySelfunc *gsf, NcmObjArray *dNdz_a);
NcmObjArray *nc_galaxy_selfunc_get_shell_splines (NcGalaxySelfunc *gsf);

void nc_galaxy_selfunc_load_from_txts (NcGalaxySelfunc *gsf, const gchar *prefix, const gchar *suffix);

gdouble nc_galaxy_selfunc_eval (NcGalaxySelfunc *gsf, const guint shell, const gdouble z);
gdouble nc_galaxy_selfunc_get_zmin (NcGalaxySelfunc *gsf, const guint shell);
gdouble nc_galaxy_selfunc_get_zmean (NcGalaxySelfunc *gsf, const guint shell);
gdouble nc_galaxy_selfunc_get_zmax (NcGalaxySelfunc *gsf, const guint shell);

G_END_DECLS

#endif /* _NC_GALAXY_SELFUNC_H_ */
