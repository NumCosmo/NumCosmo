/***************************************************************************
 *            ncm_fftlog.h
 *
 *  Fri May 18 16:44:28 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NCM_FFTLOG_H_
#define _NCM_FFTLOG_H_

#include <glib.h>
#include <glib-object.h>

#ifndef NUMCOSMO_GIR_SCAN
#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */
#endif

G_BEGIN_DECLS

#define NCM_TYPE_FFTLOG             (ncm_fftlog_get_type ())
#define NCM_FFTLOG(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FFTLOG, NcmFftlog))
#define NCM_FFTLOG_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FFTLOG, NcmFftlogClass))
#define NCM_IS_FFTLOG(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FFTLOG))
#define NCM_IS_FFTLOG_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FFTLOG))
#define NCM_FFTLOG_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FFTLOG, NcmFftlogClass))

typedef struct _NcmFftlogClass NcmFftlogClass;
typedef struct _NcmFftlog NcmFftlog;

struct _NcmFftlogClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmFftlog
{
  /*< private >*/
  GObject parent_instance;
  guint N;
  guint N_2;
  gdouble lnk0;
  gdouble lnr0;
  gdouble L;
  gdouble dr;
#ifdef NUMCOSMO_HAVE_FFTW3
  fftw_complex *in;
  fftw_complex *out;
  fftw_plan p_in2out;
  fftw_plan p_out2in;
#endif /* NUMCOSMO_HAVE_FFTW3 */
};

GType ncm_fftlog_get_type (void) G_GNUC_CONST;

G_END_DECLS

#endif /* _NCM_FFTLOG_H_ */
