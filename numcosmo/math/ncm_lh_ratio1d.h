/***************************************************************************
 *            ncm_lh_ratio1d.h
 *
 *  Fri Aug 15 15:22:57 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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

#ifndef _NCM_LH_RATIO1D_H_
#define _NCM_LH_RATIO1D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/math/ncm_fit.h>

G_BEGIN_DECLS

#define NCM_TYPE_LH_RATIO1D             (ncm_lh_ratio1d_get_type ())
#define NCM_LH_RATIO1D(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_LH_RATIO1D, NcmLHRatio1d))
#define NCM_LH_RATIO1D_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_LH_RATIO1D, NcmLHRatio1dClass))
#define NCM_IS_LH_RATIO1D(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_LH_RATIO1D))
#define NCM_IS_LH_RATIO1D_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_LH_RATIO1D))
#define NCM_LH_RATIO1D_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_LH_RATIO1D, NcmLHRatio1dClass))

typedef struct _NcmLHRatio1dClass NcmLHRatio1dClass;
typedef struct _NcmLHRatio1d NcmLHRatio1d;

/**
 * NcmLHRatio1dRoot:
 * @NCM_LH_RATIO1D_ROOT_BRACKET: FIXME
 * @NCM_LH_RATIO1D_ROOT_NUMDIFF: FIXME
 * 
 */ 
typedef enum _NcmLHRatio1dRoot
{
  NCM_LH_RATIO1D_ROOT_BRACKET = 0,
  NCM_LH_RATIO1D_ROOT_NUMDIFF,
} NcmLHRatio1dRoot;

struct _NcmLHRatio1dClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmLHRatio1d
{
  /*< private >*/
  GObject parent_instance;
  NcmFit *fit;
  NcmFit *constrained;
  NcmFitRunMsgs mtype;
  NcmLHRatio1dRoot rtype;
  NcmMSetPIndex pi;
  gdouble chisquare;
  gdouble lb;
  gdouble ub;
  gdouble bf;
  guint niter;
  guint func_eval;
  guint grad_eval;
};

GType ncm_lh_ratio1d_get_type (void) G_GNUC_CONST;

NcmLHRatio1d *ncm_lh_ratio1d_new (NcmFit *fit, NcmMSetPIndex pi);
void ncm_lh_ratio1d_free (NcmLHRatio1d *lhr1d);
void ncm_lh_ratio1d_clear (NcmLHRatio1d **lhr1d);

void ncm_lh_ratio1d_set_pindex (NcmLHRatio1d *lhr1d, NcmMSetPIndex pi);
void ncm_lh_ratio1d_find_bounds (NcmLHRatio1d *lhr1d, gdouble clevel, NcmFitRunMsgs mtype, gdouble *lb, gdouble *ub);

G_END_DECLS

#endif /* _NCM_LH_RATIO1D_H_ */

