/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert.h
 *
 *  Tue June 03 15:48:13 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert.h
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_HIPERT_H_
#define _NC_HIPERT_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT             (nc_hipert_get_type ())
#define NC_HIPERT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT, NcHIPert))
#define NC_HIPERT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT, NcHIPertClass))
#define NC_IS_HIPERT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT))
#define NC_IS_HIPERT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT))
#define NC_HIPERT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT, NcHIPertClass))

typedef struct _NcHIPertClass NcHIPertClass;
typedef struct _NcHIPert NcHIPert;
typedef struct _NcHIPertPrivate NcHIPertPrivate;

struct _NcHIPertClass
{
  /*< private >*/
  GObjectClass parent_class;
  void (*set_mode_k) (NcHIPert *pert, gdouble k);
  void (*set_reltol) (NcHIPert *pert, gdouble reltol);
  void (*set_abstol) (NcHIPert *pert, gdouble abstol);
};

struct _NcHIPert
{
  /*< private >*/
  GObject parent_instance;
  NcHIPertPrivate *priv;
};

GType nc_hipert_get_type (void) G_GNUC_CONST;

NCM_INLINE void nc_hipert_set_mode_k (NcHIPert *pert, gdouble k);
NCM_INLINE void nc_hipert_set_reltol (NcHIPert *pert, gdouble reltol);
NCM_INLINE void nc_hipert_set_abstol (NcHIPert *pert, gdouble abstol);

gdouble nc_hipert_get_reltol (NcHIPert *pert);
gdouble nc_hipert_get_abstol (NcHIPert *pert);
gdouble nc_hipert_get_mode_k (NcHIPert *pert);
gboolean nc_hipert_prepared (NcHIPert *pert);
void nc_hipert_set_prepared (NcHIPert *pert, gboolean prepared);

void nc_hipert_set_sys_size (NcHIPert *pert, guint sys_size);
void nc_hipert_set_stiff_solver (NcHIPert *pert, gboolean stiff);
void nc_hipert_reset_solver (NcHIPert *pert);

G_END_DECLS

#endif /* _NC_HIPERT_H_ */

#ifndef _NC_HIPERT_INLINE_H_
#define _NC_HIPERT_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE void 
nc_hipert_set_mode_k (NcHIPert *pert, gdouble k)
{
  NC_HIPERT_GET_CLASS (pert)->set_mode_k (pert, k);
}

NCM_INLINE void 
nc_hipert_set_reltol (NcHIPert *pert, gdouble reltol)
{
  NC_HIPERT_GET_CLASS (pert)->set_reltol (pert, reltol);
}

NCM_INLINE void 
nc_hipert_set_abstol (NcHIPert *pert, gdouble abstol)
{
  NC_HIPERT_GET_CLASS (pert)->set_abstol (pert, abstol);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPERT_INLINE_H_ */
