/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_ode.h
 *
 *  Wed December 12 11:03:16 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_ode.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_ODE_H_
#define _NCM_ODE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>

G_BEGIN_DECLS

#define NCM_TYPE_ODE             (ncm_ode_get_type ())
#define NCM_ODE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_ODE, NcmODE))
#define NCM_ODE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_ODE, NcmODEClass))
#define NCM_IS_ODE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_ODE))
#define NCM_IS_ODE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_ODE))
#define NCM_ODE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_ODE, NcmODEClass))

typedef struct _NcmODEClass NcmODEClass;
typedef struct _NcmODE NcmODE;
typedef struct _NcmODEPrivate NcmODEPrivate;

struct _NcmODEClass
{
  /*< private >*/
  GObjectClass parent_class;
  void (*set_sys_size) (NcmODE *ode, guint sys_size);
};

struct _NcmODE
{
  /*< private >*/
  GObject parent_instance;
  NcmODEPrivate *priv;
};

GType ncm_ode_get_type (void) G_GNUC_CONST;

NcmODE *ncm_ode_ref (NcmODE *ode);

void ncm_ode_free (NcmODE *ode);
void ncm_ode_clear (NcmODE **ode);

void ncm_ode_set_sys_size (NcmODE *ode, guint sys_size);
guint ncm_ode_get_sys_size (NcmODE *ode);

G_END_DECLS

#endif /* _NCM_ODE_H_ */

#ifndef _NCM_ODE_INLINE_H_
#define _NCM_ODE_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS
G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_ODE_INLINE_H_ */
