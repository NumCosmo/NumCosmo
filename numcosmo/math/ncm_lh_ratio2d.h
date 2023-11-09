/***************************************************************************
 *            ncm_lh_ratio2d.h
 *
 *  Fri Aug 15 15:22:57 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_fit.h>
#include <numcosmo/math/ncm_diff.h>

G_BEGIN_DECLS

#define NCM_TYPE_LH_RATIO2D (ncm_lh_ratio2d_get_type ())

G_DECLARE_FINAL_TYPE (NcmLHRatio2d, ncm_lh_ratio2d, NCM, LH_RATIO2D, GObject)

/**
 * NcmLHRatio2dRoot:
 * @NCM_LH_RATIO2D_ROOT_BRACKET: Root finding by bracketing
 * @NCM_LH_RATIO2D_ROOT_NUMDIFF: Root finding by numerical differentiation
 *
 * Root finding methods used by #NcmLHRatio2d.
 *
 */
typedef enum _NcmLHRatio2dRoot
{
  NCM_LH_RATIO2D_ROOT_BRACKET = 0,
  NCM_LH_RATIO2D_ROOT_NUMDIFF,
} NcmLHRatio2dRoot;


/**
 * NcmLHRatio2dPoint:
 *
 * Boxed object containing a point in the 2d parameter space.
 *
 */
typedef struct _NcmLHRatio2dPoint
{
  /*< private >*/
  gdouble x;
  gdouble y;
  gdouble theta;
  gdouble p1;
  gdouble p2;
} NcmLHRatio2dPoint;

typedef struct _NcmLHRatio2dRegion NcmLHRatio2dRegion;

/**
 * NcmLHRatio2dRegion:
 * @np: Number of points.
 * @p1: a #NcmVector containing points of parameter one.
 * @p2: a #NcmVector containing points of parameter two.
 * @clevel: the confidence level represented by the border.
 *
 * Object describing a confidence region.
 *
 */
struct _NcmLHRatio2dRegion
{
  guint np;
  NcmVector *p1;
  NcmVector *p2;
  gdouble clevel;
};

GType ncm_lh_ratio2d_region_get_type (void) G_GNUC_CONST;

NcmLHRatio2d *ncm_lh_ratio2d_new (NcmFit *fit, const NcmMSetPIndex *pi1, const NcmMSetPIndex *pi2, gdouble border_prec);
void ncm_lh_ratio2d_free (NcmLHRatio2d *lhr2d);
void ncm_lh_ratio2d_clear (NcmLHRatio2d **lhr2d);

void ncm_lh_ratio2d_set_pindex (NcmLHRatio2d *lhr2d, NcmMSetPIndex *pi1, NcmMSetPIndex *pi2);

NcmLHRatio2dRegion *ncm_lh_ratio2d_conf_region (NcmLHRatio2d *lhr2d, gdouble clevel, gdouble expected_np, NcmFitRunMsgs mtype);
NcmLHRatio2dRegion *ncm_lh_ratio2d_fisher_border (NcmLHRatio2d *lhr2d, gdouble clevel, gdouble expected_np, NcmFitRunMsgs mtype);
NcmLHRatio2dRegion *ncm_lh_ratio2d_region_dup (NcmLHRatio2dRegion *rg);
void ncm_lh_ratio2d_region_free (NcmLHRatio2dRegion *rg);
void ncm_lh_ratio2d_region_clear (NcmLHRatio2dRegion **rg);
void ncm_lh_ratio2d_region_print (NcmLHRatio2dRegion *rg, FILE *out);

G_END_DECLS

#endif /* _NCM_LH_RATIO2D_H_ */

