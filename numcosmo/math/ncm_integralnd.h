 /***************************************************************************
 *            ncm_integralnd.h
 *
 *  Thu July 20 08:39:30 2023
 *  Copyright  2023 Eduardo José Barroso
 *  <eduardo.jsbarroso@uel.br>
 ****************************************************************************/
/*
 * ncm_integralnd.h
 * Copyright (C) 2023 Eduardo José Barroso <eduardo.jsbarroso@uel.br>*
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

#ifndef _NCM_INTEGRALND_H_
#define _NCM_INTEGRALND_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NCM_TYPE_INTEGRALND             (ncm_integralnd_get_type ())
#define NCM_INTEGRALND(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_INTEGRALND, NcmIntegralnd))
#define NCM_INTEGRALND_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_INTEGRALND, NcmIntegralndClass))
#define NCM_IS_INTEGRALND(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_INTEGRALND))
#define NCM_IS_INTEGRALND_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_INTEGRALND))
#define NCM_INTEGRALND_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_INTEGRALND, NcmIntegralndClass))

typedef struct _NcmIntegralndClass NcmIntegralndClass;
typedef struct _NcmIntegralnd NcmIntegralnd;
typedef struct _NcmIntegralndPrivate NcmIntegralndPrivate;

typedef gdouble (*NcmIntegralndF) (NcmIntegralnd *intnd, const gdouble x);

struct _NcmIntegralndClass
{
  /*< private >*/
  GObjectClass parent_class;
  NcmIntegralndF integrand;
  gpointer padding[9];
};

struct _NcmIntegralnd
{
  /*< private >*/
  GObject parent_instance;
  NcmIntegralndPrivate *priv;
};

GType ncm_integralnd_get_type (void) G_GNUC_CONST;

NcmIntegralnd *ncm_integralnd_ref (NcmIntegralnd *intnd);
void ncm_integralnd_free (NcmIntegralnd *intnd);
void ncm_integralnd_clear (NcmIntegralnd **intnd);

void ncm_integralnd_set_partition (NcmIntegralnd *intnd, guint partition);
void ncm_integralnd_set_rule (NcmIntegralnd *intnd, guint rule);
void ncm_integralnd_set_reltol (NcmIntegralnd *intnd, gdouble reltol);
void ncm_integralnd_set_abstol (NcmIntegralnd *intnd, gdouble abstol);

guint ncm_integralnd_get_partition (NcmIntegralnd *intnd);
guint ncm_integralnd_get_rule (NcmIntegralnd *intnd);
gdouble ncm_integralnd_get_reltol (NcmIntegralnd *intnd);
gdouble ncm_integralnd_get_abstol (NcmIntegralnd *intnd);

NCM_INLINE gdouble ncm_integralnd_integrand (NcmIntegralnd *intnd, const gdouble x);

gdouble ncm_integralnd_eval_h (NcmIntegralnd *intnd, const gdouble xi, const gdouble xf, gdouble *err);

#define NCM_INTEGRALND_DEFAULT_PARTITION 100000
#define NCM_INTEGRALND_DEFAULT_ALG 6
#define NCM_INTEGRALND_DEFAULT_RELTOL 1e-13
#define NCM_INTEGRALND_DEFAULT_ABSTOL 0.0

G_END_DECLS

#endif /* _NCM_INTEGRALND_H_ */

#ifndef _NCM_INTEGRALND_INLINE_H_
#define _NCM_INTEGRALND_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gdouble 
ncm_integralnd_integrand (NcmIntegralnd *intnd, const gdouble x)
{
  return NCM_INTEGRALND_GET_CLASS (intnd)->integrand (intnd, x);  
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_INTEGRALND_INLINE_H_ */

