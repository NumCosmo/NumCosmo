/***************************************************************************
 *            ncm_fftlog_sbessel_j.h
 *
 *  Wed July 19 09:58:50 2017
 *  Copyright  2017  Fernando de Simoni
 *  <fernando.saliby@gmail.com>
 ****************************************************************************/
/***************************************************************************
 *            ncm_fftlog_sbessel_j.h
 *
 *  Sat September 02 18:11:24 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fftlog_sbessel_j.h
 *
 * Copyright (C) 2017 - Fernando de Simoni
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NCM_FFTLOG_SBESSEL_J_H_
#define _NCM_FFTLOG_SBESSEL_J_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fftlog.h>

#ifndef NUMCOSMO_GIR_SCAN
#ifdef HAVE_ACB_H
#include <acb.h>
#endif /* HAVE_ACB_H */
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_FFTLOG_SBESSEL_J             (ncm_fftlog_sbessel_j_get_type ())
#define NCM_FFTLOG_SBESSEL_J(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FFTLOG_SBESSEL_J, NcmFftlogSBesselJ))
#define NCM_FFTLOG_SBESSEL_J_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FFTLOG_SBESSEL_J, NcmFftlogSBesselJClass))
#define NCM_IS_FFTLOG_SBESSEL_J(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FFTLOG_SBESSEL_J))
#define NCM_IS_FFTLOG_SBESSEL_J_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FFTLOG_SBESSEL_J))
#define NCM_FFTLOG_SBESSEL_J_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FFTLOG_SBESSEL_J, NcmFftlogSBesselJClass))

typedef struct _NcmFftlogSBesselJClass NcmFftlogSBesselJClass;
typedef struct _NcmFftlogSBesselJ NcmFftlogSBesselJ;
typedef struct _NcmFftlogSBesselJPrivate NcmFftlogSBesselJPrivate;

struct _NcmFftlogSBesselJClass
{
  /*< private >*/
  NcmFftlogClass parent_class;
};

struct _NcmFftlogSBesselJ
{
  /*< private >*/
  NcmFftlog parent_instance;
  NcmFftlogSBesselJPrivate *priv;
};

GType ncm_fftlog_sbessel_j_get_type (void) G_GNUC_CONST;

NcmFftlogSBesselJ *ncm_fftlog_sbessel_j_new (guint ell, gdouble lnr0, gdouble lnk0, gdouble Lk, guint N);

void ncm_fftlog_sbessel_j_set_ell (NcmFftlogSBesselJ *fftlog_jl, const guint ell);
guint ncm_fftlog_sbessel_j_get_ell (NcmFftlogSBesselJ *fftlog_jl);

void ncm_fftlog_sbessel_j_set_q (NcmFftlogSBesselJ *fftlog_jl, const gdouble q);
gdouble ncm_fftlog_sbessel_j_get_q (NcmFftlogSBesselJ *fftlog_jl);

void ncm_fftlog_sbessel_j_set_best_lnr0 (NcmFftlogSBesselJ *fftlog_jl);
void ncm_fftlog_sbessel_j_set_best_lnk0 (NcmFftlogSBesselJ *fftlog_jl);

G_END_DECLS

#endif /* _NCM_FFTLOG_SBESSEL_J_H_ */
