/***************************************************************************
 *            ncm_fit_levmar.h
 *
 *  Wed Feb 24 21:19:48 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_FIT_LEVMAR_H_
#define _NCM_FIT_LEVMAR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>
#ifdef NUMCOSMO_HAVE_LEVMAR

G_BEGIN_DECLS

#define NCM_TYPE_FIT_LEVMAR             (ncm_fit_levmar_get_type ())
#define NCM_FIT_LEVMAR(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT_LEVMAR, NcmFitLevmar))
#define NCM_FIT_LEVMAR_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT_LEVMAR, NcmFitLevmarClass))
#define NCM_IS_FIT_LEVMAR(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT_LEVMAR))
#define NCM_IS_FIT_LEVMAR_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT_LEVMAR))
#define NCM_FIT_LEVMAR_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT_LEVMAR, NcmFitLevmarClass))

typedef struct _NcmFitLevmarClass NcmFitLevmarClass;
typedef struct _NcmFitLevmar NcmFitLevmar;

/**
 * NcmFitLevmarAlgos:
 * @NCM_FIT_LEVMAR_DER: FIXME
 * @NCM_FIT_LEVMAR_DIF: FIXME
 *
 * FIXME
 */
typedef enum _NcmFitLevmarAlgos
{
  NCM_FIT_LEVMAR_DER = 0,
  NCM_FIT_LEVMAR_DIF,       /*< private >*/
  NCM_FIT_LEVMAR_NUM_ALGOS, /*< skip >*/
} NcmFitLevmarAlgos;

struct _NcmFitLevmarClass
{
  /*< private >*/
  NcmFitClass parent_class;
};

struct _NcmFitLevmar
{
  /*< private >*/
  NcmFit parent_instance;
  gpointer dif;
  gpointer der;
  guint fparam_len;
  guint data_len;
  NcmFitLevmarAlgos algo;
};

GType ncm_fit_levmar_get_type (void) G_GNUC_CONST;

NcmFit *ncm_fit_levmar_new (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, NcmFitLevmarAlgos algo);
NcmFit *ncm_fit_levmar_new_default (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype);
NcmFit *ncm_fit_levmar_new_by_name (NcmLikelihood *lh, NcmMSet *mset, NcmFitGradType gtype, gchar *algo_name);
void ncm_fit_levmar_set_algo (NcmFitLevmar *fit_levmar, NcmFitLevmarAlgos algo);

G_END_DECLS

#endif /* NUMCOSMO_HAVE_LEVMAR */
#endif /* _NCM_FIT_LEVMAR_H_ */

