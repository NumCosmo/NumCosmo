/***************************************************************************
 *            ncm_fftlog_j1pow2.h
 *
 *  Mon July 21 19:59:56 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fftlog_j1pow2.h
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

#ifndef _NCM_FFTLOG_J1POW2_H_
#define _NCM_FFTLOG_J1POW2_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fftlog.h>

#include <glib-object.h>

G_BEGIN_DECLS

#define NCM_TYPE_FFTLOG_J1POW2             (ncm_fftlog_j1pow2_get_type ())
#define NCM_FFTLOG_J1POW2(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FFTLOG_J1POW2, NcmFftlogJ1pow2))
#define NCM_FFTLOG_J1POW2_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FFTLOG_J1POW2, NcmFftlogJ1pow2Class))
#define NCM_IS_FFTLOG_J1POW2(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FFTLOG_J1POW2))
#define NCM_IS_FFTLOG_J1POW2_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FFTLOG_J1POW2))
#define NCM_FFTLOG_J1POW2_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FFTLOG_J1POW2, NcmFftlogJ1pow2Class))

typedef struct _NcmFftlogJ1pow2Class NcmFftlogJ1pow2Class;
typedef struct _NcmFftlogJ1pow2 NcmFftlogJ1pow2;

struct _NcmFftlogJ1pow2Class
{
  NcmFftlogClass parent_class;
};

struct _NcmFftlogJ1pow2
{
  NcmFftlog parent_instance;
};

GType ncm_fftlog_j1pow2_get_type (void) G_GNUC_CONST;

NcmFftlog *ncm_fftlog_j1pow2_new (gdouble r0, gdouble k0, gdouble Lk, guint N);

G_END_DECLS

#endif /* _NCM_FFTLOG_J1POW2_H_ */
