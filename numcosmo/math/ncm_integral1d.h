/***************************************************************************
 *            ncm_integral1d.h
 *
 *  Sat February 20 14:29:48 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_integral1d.h
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

#ifndef _NCM_INTEGRAL1D_H_
#define _NCM_INTEGRAL1D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NCM_TYPE_INTEGRAL1D             (ncm_integral1d_get_type ())
#define NCM_INTEGRAL1D(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_INTEGRAL1D, NcmIntegral1d))
#define NCM_INTEGRAL1D_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_INTEGRAL1D, NcmIntegral1dClass))
#define NCM_IS_INTEGRAL1D(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_INTEGRAL1D))
#define NCM_IS_INTEGRAL1D_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_INTEGRAL1D))
#define NCM_INTEGRAL1D_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_INTEGRAL1D, NcmIntegral1dClass))

typedef struct _NcmIntegral1dClass NcmIntegral1dClass;
typedef struct _NcmIntegral1d NcmIntegral1d;
typedef struct _NcmIntegral1dPrivate NcmIntegral1dPrivate;

typedef gdouble (*NcmIntegral1dF) (NcmIntegral1d *int1d, const gdouble x, const gdouble w);

struct _NcmIntegral1dClass
{
  /*< private >*/
  GObjectClass parent_class;
  NcmIntegral1dF integrand;
  gpointer padding[9];
};

struct _NcmIntegral1d
{
  /*< private >*/
  GObject parent_instance;
  NcmIntegral1dPrivate *priv;
};

GType ncm_integral1d_get_type (void) G_GNUC_CONST;

NcmIntegral1d *ncm_integral1d_ref (NcmIntegral1d *int1d);
void ncm_integral1d_free (NcmIntegral1d *int1d);
void ncm_integral1d_clear (NcmIntegral1d **int1d);

void ncm_integral1d_set_partition (NcmIntegral1d *int1d, guint partition);
void ncm_integral1d_set_rule (NcmIntegral1d *int1d, guint rule);
void ncm_integral1d_set_reltol (NcmIntegral1d *int1d, gdouble reltol);
void ncm_integral1d_set_abstol (NcmIntegral1d *int1d, gdouble abstol);

guint ncm_integral1d_get_partition (NcmIntegral1d *int1d);
guint ncm_integral1d_get_rule (NcmIntegral1d *int1d);
gdouble ncm_integral1d_get_reltol (NcmIntegral1d *int1d);
gdouble ncm_integral1d_get_abstol (NcmIntegral1d *int1d);

NCM_INLINE gdouble ncm_integral1d_integrand (NcmIntegral1d *int1d, const gdouble x, const gdouble w);

gdouble ncm_integral1d_eval (NcmIntegral1d *int1d, const gdouble xi, const gdouble xf, gdouble *err);
gdouble ncm_integral1d_eval_gauss_hermite_p (NcmIntegral1d *int1d, gdouble *err);
gdouble ncm_integral1d_eval_gauss_hermite (NcmIntegral1d *int1d, gdouble *err);
gdouble ncm_integral1d_eval_gauss_hermite_r_p (NcmIntegral1d *int1d, const gdouble r, gdouble *err);
gdouble ncm_integral1d_eval_gauss_hermite_mur (NcmIntegral1d *int1d, const gdouble r, const gdouble mu, gdouble *err);

gdouble ncm_integral1d_eval_gauss_hermite1_p (NcmIntegral1d *int1d, gdouble *err);
gdouble ncm_integral1d_eval_gauss_hermite1_r_p (NcmIntegral1d *int1d, const gdouble r, gdouble *err);

gdouble ncm_integral1d_eval_gauss_laguerre (NcmIntegral1d *int1d, gdouble *err);
gdouble ncm_integral1d_eval_gauss_laguerre_r (NcmIntegral1d *int1d, const gdouble r, gdouble *err);

#define NCM_INTEGRAL1D_DEFAULT_PARTITION 100000
#define NCM_INTEGRAL1D_DEFAULT_ALG 6
#define NCM_INTEGRAL1D_DEFAULT_RELTOL 1e-13
#define NCM_INTEGRAL1D_DEFAULT_ABSTOL 0.0

G_END_DECLS

#endif /* _NCM_INTEGRAL1D_H_ */

#ifndef _NCM_INTEGRAL1D_INLINE_H_
#define _NCM_INTEGRAL1D_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gdouble 
ncm_integral1d_integrand (NcmIntegral1d *int1d, const gdouble x, const gdouble w)
{
  return NCM_INTEGRAL1D_GET_CLASS (int1d)->integrand (int1d, x, w);  
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_INTEGRAL1D_INLINE_H_ */

