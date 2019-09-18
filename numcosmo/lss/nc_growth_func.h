/***************************************************************************
 *            nc_growth_func.h
 *
 *  Tue Apr  6 01:12:58 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_GROWTH_FUNC_H_
#define _NC_GROWTH_FUNC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_GROWTH_FUNC             (nc_growth_func_get_type ())
#define NC_GROWTH_FUNC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_GROWTH_FUNC, NcGrowthFunc))
#define NC_GROWTH_FUNC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_GROWTH_FUNC, NcGrowthFuncClass))
#define NC_IS_GROWTH_FUNC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_GROWTH_FUNC))
#define NC_IS_GROWTH_FUNC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_GROWTH_FUNC))
#define NC_GROWTH_FUNC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_GROWTH_FUNC, NcGrowthFuncClass))

typedef struct _NcGrowthFuncClass NcGrowthFuncClass;
typedef struct _NcGrowthFunc NcGrowthFunc;
typedef struct _NcGrowthFuncPrivate NcGrowthFuncPrivate;

struct _NcGrowthFuncClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcGrowthFunc
{
  /*< private >*/
  GObject parent_instance;
  NcGrowthFuncPrivate *priv;
  NcmSpline *s;
  gdouble Da0;
};

GType nc_growth_func_get_type (void) G_GNUC_CONST;

NcGrowthFunc * nc_growth_func_new (void);
NcGrowthFunc * nc_growth_func_ref (NcGrowthFunc *gf);
void nc_growth_func_free (NcGrowthFunc *gf);
void nc_growth_func_clear (NcGrowthFunc **gf);

void nc_growth_func_prepare (NcGrowthFunc * gf, NcHICosmo *cosmo);
void nc_growth_func_prepare_if_needed (NcGrowthFunc *gf, NcHICosmo *cosmo);

NCM_INLINE gdouble nc_growth_func_eval (NcGrowthFunc *gf, NcHICosmo *cosmo, gdouble z);
NCM_INLINE gdouble nc_growth_func_eval_deriv (NcGrowthFunc *gf, NcHICosmo *cosmo, gdouble z);
NCM_INLINE void nc_growth_func_eval_both (NcGrowthFunc *gf, NcHICosmo *cosmo, gdouble z, gdouble *d, gdouble *f);
NCM_INLINE gdouble nc_growth_func_get_dust_norma_Da0 (NcGrowthFunc *gf);

G_END_DECLS

#endif /* _NC_GROWTH_FUNC_H_ */

#ifndef _NC_GROWTH_FUNC_INLINE_H_
#define _NC_GROWTH_FUNC_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gdouble 
nc_growth_func_eval (NcGrowthFunc *gf, NcHICosmo *cosmo, gdouble z)
{
  const gdouble a = 1.0 / (1.0 + z);
  return ncm_spline_eval (gf->s, a);
}

NCM_INLINE gdouble
nc_growth_func_eval_deriv (NcGrowthFunc *gf, NcHICosmo *cosmo, gdouble z)
{
  const gdouble a = 1.0 / (1.0 + z);
  return - a * a * ncm_spline_eval_deriv (gf->s, a);
}

NCM_INLINE void
nc_growth_func_eval_both (NcGrowthFunc *gf, NcHICosmo *cosmo, gdouble z, gdouble *d, gdouble *f)
{
  const gdouble a = 1.0 / (1.0 + z);
  *d = ncm_spline_eval (gf->s, a);
  *f = - a * a * ncm_spline_eval_deriv (gf->s, a);
}

NCM_INLINE gdouble 
nc_growth_func_get_dust_norma_Da0 (NcGrowthFunc *gf)
{
  return gf->Da0;
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_GROWTH_FUNC_INLINE_H_ */
