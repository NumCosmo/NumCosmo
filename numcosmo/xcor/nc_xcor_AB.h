/***************************************************************************
 *            nc_xcor_AB.h
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

#ifndef _NC_XCOR_AB_H_
#define _NC_XCOR_AB_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_matrix.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_AB (nc_xcor_AB_get_type ())
#define NC_XCOR_AB(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_XCOR_AB, NcXcorAB))
#define NC_XCOR_AB_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_XCOR_AB, NcXcorABClass))
#define NC_IS_XCOR_AB(obj) (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_XCOR_AB))
#define NC_IS_XCOR_AB_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_XCOR_AB))
#define NC_XCOR_AB_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_XCOR_AB, NcXcorABClass))

typedef struct _NcXcorABClass NcXcorABClass;
typedef struct _NcXcorAB NcXcorAB;

struct _NcXcorAB
{
  /*< private > */
  GObject parent_instance;
  
  guint a;
  guint b;
  
  guint ell_th_cut_off;
  guint ell_lik_min;
  guint ell_lik_max;
  guint nell_lik;
  
  NcmMatrix *mixing;
  NcmMatrix *cl_th; /*column 0 : C_l^th, 1 : C_l^th+N_l, 2 : mixed C_l */
  NcmVector *cl_obs;
};

struct _NcXcorABClass
{
  /*< private > */
  GObjectClass parent_class;
  
  gpointer (*alloc) (void);
};

GType nc_xcor_AB_get_type (void) G_GNUC_CONST;


/* void nc_xcor_AB_prepare (NcXcorAB* xc, NcHICosmo* model); */
/* NcXcorAB* nc_xcor_AB_new (NcDistance* dist, NcTransferFunc* tf, NcGrowthFunc* gf, gdouble zl, gdouble zu); */
NcXcorAB *nc_xcor_AB_new (guint a, guint b, guint ell_th_cut_off, guint ell_lik_min, guint ell_lik_max, const gchar *clobs_filename, const gchar *mixing_filename, const guint mixing_filelength);
NcXcorAB *nc_xcor_AB_ref (NcXcorAB *xcab);
void nc_xcor_AB_free (NcXcorAB *xcab);
void nc_xcor_AB_clear (NcXcorAB **xcab);

/* gdouble nc_xcor_AB_limber_cross_cl (NcXcorAB* xc, NcXcorABLimberKernel* xclk1, NcXcorABLimberKernel* xclk2, NcHICosmo* cosmo, guint l); */
/* void nc_xcor_AB_limber_cross_cl_range (NcXcorAB* xc, NcXcorABLimberKernel* xclk1, NcXcorABLimberKernel* xclk2, NcHICosmo* cosmo, guint lmin, guint lmax, NcmVector* vp); */
/* gdouble nc_xcor_AB_limber_auto_cl (NcXcorAB* xc, NcXcorABLimberKernel* xclk, NcHICosmo* cosmo, guint l, gboolean withnoise); */
/* void nc_xcor_AB_limber_auto_cl_range (NcXcorAB* xc, NcXcorABLimberKernel* xclk, NcHICosmo* cosmo, guint lmin, guint lmax, NcmVector* vp, gboolean withnoise); */
/* void nc_xcor_AB_limber_auto_cl_clnl (NcXcorAB* xc, NcXcorABLimberKernel* xclk, NcHICosmo* cosmo, guint lmax, NcmVector* vcl, NcmVector* vclnl); */

/* H0/c in h Mpc-1: */
#define H0_c_3 (gsl_pow_3 (1e5 / ncm_c_c ()))

G_END_DECLS

#endif /* _NC_XCOR_AB_H_ */

