/***************************************************************************
 *            nc_xcor_limber.h
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

#ifndef _NC_XCOR_LIMBER_H_
#define _NC_XCOR_LIMBER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/lss/nc_transfer_func.h>
#include <numcosmo/lss/nc_growth_func.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR_LIMBER (nc_xcor_limber_get_type ())
#define NC_XCOR_LIMBER(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_XCOR_LIMBER, NcXcorLimber))
#define NC_XCOR_LIMBER_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_XCOR_LIMBER, NcXcorLimberClass))
#define NC_IS_XCOR_LIMBER(obj) (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_XCOR_LIMBER))
#define NC_IS_XCOR_LIMBER_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_XCOR_LIMBER))
#define NC_XCOR_LIMBER_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_XCOR_LIMBER, NcXcorLimberClass))

typedef struct _NcXcorLimberClass NcXcorLimberClass;
typedef struct _NcXcorLimber NcXcorLimber;

/**
 * NcXcorLimberImpl:
 * @NC_XCOR_LIMBER_EVAL_KERNEL: FIXME
 * @NC_XCOR_LIMBER_PREPARE: FIXME
 * 
 */
typedef enum _NcXcorLimberImpl
{
	NC_XCOR_LIMBER_EVAL_KERNEL = 1 << 0,
	NC_XCOR_LIMBER_PREPARE = 1 << 1,
} NcXcorLimberImpl;

#define NC_XCOR_LIMBER_IMPL_ALL (~0)

struct _NcXcorLimber
{
	/*< private >*/
	NcmModel parent_instance;
	gdouble cons_factor;
	NcmModelCtrl* cosmo_ctrl;
};

struct _NcXcorLimberClass
{
	/*< private >*/
	NcmModelClass parent_class;
	gdouble (*eval_kernel)(NcXcorLimber* xcl, NcHICosmo* cosmo, gdouble z, gint l);
	void (*prepare)(NcXcorLimber* xcl, NcHICosmo* cosmo);
	gdouble (*noise_spec)(NcXcorLimber* xcl, guint l);

	guint (*obs_len)(NcXcorLimber* xcl);
	guint (*obs_params_len)(NcXcorLimber* xcl);
	NcXcorLimberImpl impl;
};

GType nc_xcor_limber_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_xcor_limber);

NcXcorLimber* nc_xcor_limber_new_from_name (gchar* xcor_name);
NcXcorLimber* nc_xcor_limber_ref (NcXcorLimber* xcl);
void nc_xcor_limber_free (NcXcorLimber* xcl);
void nc_xcor_limber_clear (NcXcorLimber** xcl);

NcXcorLimberImpl nc_xcor_limber_impl (NcXcorLimber* xcl);

guint nc_xcor_limber_obs_len (NcXcorLimber* xcl);
guint nc_xcor_limber_obs_params_len (NcXcorLimber* xcl);

gdouble nc_xcor_limber_eval_kernel (NcXcorLimber* xcl, NcHICosmo* cosmo, gdouble z, gint l);
void nc_xcor_limber_prepare (NcXcorLimber* xcl, NcHICosmo* cosmo);
gdouble nc_xcor_limber_noise_spec (NcXcorLimber* xcl, guint l);

void nc_xcor_limber_log_all_models (void);


G_END_DECLS

#endif /* _NC_XCOR_LIMBER_H_ */
