/***************************************************************************
 *            ncm_hoaa.h
 *
 *  Fri November 04 13:27:40 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_hoaa.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_HOAA_H_
#define _NCM_HOAA_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>

G_BEGIN_DECLS

#define NCM_TYPE_HOAA             (ncm_hoaa_get_type ())
#define NCM_HOAA(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_HOAA, NcmHOAA))
#define NCM_HOAA_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_HOAA, NcmHOAAClass))
#define NCM_IS_HOAA(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_HOAA))
#define NCM_IS_HOAA_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_HOAA))
#define NCM_HOAA_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_HOAA, NcmHOAAClass))

typedef struct _NcmHOAAClass NcmHOAAClass;
typedef struct _NcmHOAA NcmHOAA;
typedef struct _NcmHOAAPrivate NcmHOAAPrivate;

struct _NcmHOAAClass
{
  /*< private >*/
  GObjectClass parent_class;
  gdouble (*eval_m) (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_nu) (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_V) (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
  void (*eval_system) (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu, gdouble *Nt);
  void (*prepare) (NcmHOAA *hoaa, NcmModel *model);
};

struct _NcmHOAA
{
  /*< private >*/
  GObject parent_instance;
  NcmHOAAPrivate *priv;
  gdouble k;
};

GType ncm_hoaa_get_type (void) G_GNUC_CONST;

NcmHOAA *ncm_hoaa_ref (NcmHOAA *hoaa);
void ncm_hoaa_free (NcmHOAA *hoaa);
void ncm_hoaa_clear (NcmHOAA **hoaa);

void ncm_hoaa_set_reltol (NcmHOAA *hoaa, const gdouble reltol);
void ncm_hoaa_set_abstol (NcmHOAA *hoaa, const gdouble abstol);
void ncm_hoaa_set_k (NcmHOAA *hoaa, const gdouble k);
void ncm_hoaa_set_ti (NcmHOAA *hoaa, const gdouble ti);
void ncm_hoaa_set_tf (NcmHOAA *hoaa, const gdouble tf);

void ncm_hoaa_save_evol (NcmHOAA *hoaa, gboolean save_evol);
void ncm_hoaa_prepare (NcmHOAA *hoaa, NcmModel *model);

void ncm_hoaa_get_t0_t1 (NcmHOAA *hoaa, NcmModel *model, gdouble *t0, gdouble *t1);

void ncm_hoaa_eval_adiabatic_approx (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *RQ, gdouble *IQ, gdouble *RLnI, gdouble *ILnI);
void ncm_hoaa_eval_adiabatic_LnI_approx (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble RQ, const gdouble IQ, gdouble *RLnI, gdouble *ILnI);
void ncm_hoaa_eval_AA (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *RQ, gdouble *IQ, gdouble *RLnI, gdouble *ILnI);
void ncm_hoaa_eval_CV (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *Rphi, gdouble *Iphi, gdouble *RPphi, gdouble *IPphi);
void ncm_hoaa_eval_Delta (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *Delta_phi, gdouble *Delta_Pphi);

G_INLINE_FUNC gdouble ncm_hoaa_eval_m (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
G_INLINE_FUNC gdouble ncm_hoaa_eval_nu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
G_INLINE_FUNC gdouble ncm_hoaa_eval_V (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
G_INLINE_FUNC void ncm_hoaa_eval_system (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu, gdouble *Nt);

#define NCM_HOAA_RQi      (0)
#define NCM_HOAA_IQi      (1)
#define NCM_HOAA_RLnIi    (2)
#define NCM_HOAA_ILnIi    (3)
#define NCM_HOAA_SYS_SIZE (4)
#define NCM_HOAA_MAX_ANGLE (100.0)

G_END_DECLS

#endif /* _NCM_HOAA_H_ */

#ifndef _NCM_HOAA_INLINE_H_
#define _NCM_HOAA_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC gdouble 
ncm_hoaa_eval_m (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_HOAA_GET_CLASS (hoaa)->eval_m (hoaa, model, t, k);
}

G_INLINE_FUNC gdouble 
ncm_hoaa_eval_nu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_HOAA_GET_CLASS (hoaa)->eval_nu (hoaa, model, t, k);
}

G_INLINE_FUNC gdouble 
ncm_hoaa_eval_V (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_HOAA_GET_CLASS (hoaa)->eval_V (hoaa, model, t, k);
}

G_INLINE_FUNC void 
ncm_hoaa_eval_system (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu, gdouble *Nt)
{
  NCM_HOAA_GET_CLASS (hoaa)->eval_system (hoaa, model, t, k, nu, dlnmnu, Vnu, Nt);
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_HOAA_INLINE_H_ */
