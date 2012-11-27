/***************************************************************************
 *            ncm_lh_ratio2d.h
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

#ifndef _NCM_LH_RATIO2D_H_
#define _NCM_LH_RATIO2D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/math/ncm_fit.h>

G_BEGIN_DECLS

#define NCM_TYPE_LH_RATIO2D             (ncm_lh_ratio2d_get_type ())
#define NCM_LH_RATIO2D(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_LH_RATIO2D, NcmLHRatio2d))
#define NCM_LH_RATIO2D_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_LH_RATIO2D, NcmLHRatio2dClass))
#define NCM_IS_LH_RATIO2D(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_LH_RATIO2D))
#define NCM_IS_LH_RATIO2D_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_LH_RATIO2D))
#define NCM_LH_RATIO2D_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_LH_RATIO2D, NcmLHRatio2dClass))

typedef struct _NcmLHRatio2dClass NcmLHRatio2dClass;
typedef struct _NcmLHRatio2d NcmLHRatio2d;

/**
 * NcmLHRatio2dRoot:
 * @NCM_LH_RATIO2D_ROOT_BRACKET: FIXME
 * @NCM_LH_RATIO2D_ROOT_NUMDIFF: FIXME
 * 
 */ 
typedef enum _NcmLHRatio2dRoot
{
  NCM_LH_RATIO2D_ROOT_BRACKET = 0,
  NCM_LH_RATIO2D_ROOT_NUMDIFF,
} NcmLHRatio2dRoot;

struct _NcmLHRatio2dClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmLHRatio2d
{
  /*< private >*/
  GObject parent_instance;
  NcmFit *fit;
  NcmFit *constrained;
  NcmFitRunMsgs mtype;
  NcmLHRatio2dRoot rtype;
  NcmMSetPIndex pi[2];
  gdouble chisquare;
  gdouble lb[2];
  gdouble ub[2];
  gdouble bf[2];
  NcmMatrix *e_vec;
  NcmVector *e_val;
  gdouble r, theta;
  gdouble shift[2];
  gboolean angular;
  guint niter;
  guint func_eval;
  guint grad_eval;
};

typedef struct _NcmLHRatio2dPoint NcmLHRatio2dPoint;

/**
 * NcmLHRatio2dPoint:
 *
 * FIXME
 */
struct _NcmLHRatio2dPoint
{
  /*< private >*/
  gdouble x;
  gdouble y;
  gdouble theta;
  gdouble p1;
  gdouble p2;
};

GType ncm_lh_ratio2d_get_type (void) G_GNUC_CONST;

NcmLHRatio2d *ncm_lh_ratio2d_new (NcmFit *fit, NcmMSetPIndex pi1, NcmMSetPIndex pi2);
void ncm_lh_ratio2d_free (NcmLHRatio2d *lhr2d);
void ncm_lh_ratio2d_clear (NcmLHRatio2d **lhr2d);

void ncm_lh_ratio2d_set_pindex (NcmLHRatio2d *lhr2d, NcmMSetPIndex pi1, NcmMSetPIndex pi2);

GList *ncm_lh_ratio2d_conf_region (NcmLHRatio2d *lhr2d, gdouble clevel, NcmFitRunMsgs mtype);
GList *ncm_lh_ratio2d_fisher_border (NcmLHRatio2d *lhr2d, gdouble clevel, NcmFitRunMsgs mtype);
void ncm_lh_ratio2d_points_free (GList *points);
void ncm_lh_ratio2d_points_print (GList *points, FILE *out);

G_END_DECLS

#endif /* _NCM_LH_RATIO2D_H_ */

