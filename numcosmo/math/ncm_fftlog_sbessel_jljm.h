/***************************************************************************
 *            ncm_fftlog_sbessel_jljm.h
 *
 *  Sun March 24 16:54:15 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fftlog_sbessel_jljm.h
 *
 * Copyright (C) 2017 - Sandro Dias Pinto Vitenti
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

#ifndef _NCM_FFTLOG_SBESSEL_JLJM_H_
#define _NCM_FFTLOG_SBESSEL_JLJM_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fftlog.h>

G_BEGIN_DECLS

#define NCM_TYPE_FFTLOG_SBESSEL_JLJM             (ncm_fftlog_sbessel_jljm_get_type ())
#define NCM_FFTLOG_SBESSEL_JLJM(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FFTLOG_SBESSEL_JLJM, NcmFftlogSBesselJLJM))
#define NCM_FFTLOG_SBESSEL_JLJM_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FFTLOG_SBESSEL_JLJM, NcmFftlogSBesselJLJMClass))
#define NCM_IS_FFTLOG_SBESSEL_JLJM(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FFTLOG_SBESSEL_JLJM))
#define NCM_IS_FFTLOG_SBESSEL_JLJM_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FFTLOG_SBESSEL_JLJM))
#define NCM_FFTLOG_SBESSEL_JLJM_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FFTLOG_SBESSEL_JLJM, NcmFftlogSBesselJLJMClass))

typedef struct _NcmFftlogSBesselJLJMClass NcmFftlogSBesselJLJMClass;
typedef struct _NcmFftlogSBesselJLJM NcmFftlogSBesselJLJM;
typedef struct _NcmFftlogSBesselJLJMPrivate NcmFftlogSBesselJLJMPrivate;

struct _NcmFftlogSBesselJLJMClass
{
  /*< private >*/
  NcmFftlogClass parent_class;
};

struct _NcmFftlogSBesselJLJM
{
  /*< private >*/
  NcmFftlog parent_instance;
  NcmFftlogSBesselJLJMPrivate *priv;
};

GType ncm_fftlog_sbessel_jljm_get_type (void) G_GNUC_CONST;

NcmFftlogSBesselJLJM *ncm_fftlog_sbessel_jljm_new (gint ell, gint dell, gdouble lnw, gdouble lnr0, gdouble lnk0, gdouble Lk, guint N);

void ncm_fftlog_sbessel_jljm_set_ell (NcmFftlogSBesselJLJM *fftlog_jljm, const gint ell);
gint ncm_fftlog_sbessel_jljm_get_ell (NcmFftlogSBesselJLJM *fftlog_jljm);

void ncm_fftlog_sbessel_jljm_set_dell (NcmFftlogSBesselJLJM *fftlog_jljm, const gint dell);
gint ncm_fftlog_sbessel_jljm_get_dell (NcmFftlogSBesselJLJM *fftlog_jljm);

void ncm_fftlog_sbessel_jljm_set_lnw (NcmFftlogSBesselJLJM *fftlog_jljm, const gdouble lnw);
gdouble ncm_fftlog_sbessel_jljm_get_lnw (NcmFftlogSBesselJLJM *fftlog_jljm);

void ncm_fftlog_sbessel_jljm_set_best_lnr0 (NcmFftlogSBesselJLJM *fftlog_jljm);
void ncm_fftlog_sbessel_jljm_set_best_lnk0 (NcmFftlogSBesselJLJM *fftlog_jljm);

G_END_DECLS

#endif /* _NCM_FFTLOG_SBESSEL_JLJM_H_ */
