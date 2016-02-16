/***************************************************************************
 *            nc_xcor.h
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
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

#ifndef _NC_XCOR_H_
#define _NC_XCOR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/lss/nc_transfer_func.h>
#include <numcosmo/lss/nc_growth_func.h>
#include <numcosmo/xcor/nc_xcor_limber.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR (nc_xcor_get_type ())
#define NC_XCOR(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_XCOR, NcXcor))
#define NC_XCOR_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_XCOR, NcXcorClass))
#define NC_IS_XCOR(obj) (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_XCOR))
#define NC_IS_XCOR_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_XCOR))
#define NC_XCOR_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_XCOR, NcXcorClass))

typedef struct _NcXcorClass NcXcorClass;
typedef struct _NcXcor NcXcor;

struct _NcXcor
{
	/*< private > */
	GObject parent_instance;

	NcmModelCtrl* ctrlcosmo;
	NcmModelCtrl* ctrl1;
	NcmModelCtrl* ctrl2;

	NcDistance* dist;
	NcTransferFunc* tf;
	NcGrowthFunc* gf;

	gdouble normPS;
	gdouble zl, zu;
};

struct _NcXcorClass
{
	/*< private > */
	GObjectClass parent_class;
	gpointer (*alloc) (void);
};

GType nc_xcor_get_type (void) G_GNUC_CONST;


void nc_xcor_prepare (NcXcor *xc, NcHICosmo *model);
NcXcor* nc_xcor_new (NcDistance* dist, NcTransferFunc* tf, NcGrowthFunc* gf, gdouble zl, gdouble zu);
NcXcor* nc_xcor_ref (NcXcor* xcl);
void nc_xcor_free (NcXcor* xcl);
void nc_xcor_clear (NcXcor** xcl);

void nc_xcor_limber_cross_cl (NcXcor* xc, NcXcorLimber* xcl1, NcXcorLimber* xcl2, NcHICosmo* cosmo, NcmVector* ell, NcmVector* vp, guint lmin_idx);
void nc_xcor_limber_auto_cl (NcXcor* xc, NcXcorLimber* xcl, NcHICosmo* cosmo, NcmVector* ell, NcmVector* vp, guint lmin_idx, gboolean withnoise);


G_END_DECLS

#endif /* _NC_XCOR_H_ */
