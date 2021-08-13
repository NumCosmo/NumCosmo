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

#include <glib-object.h>
#include <glib.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/xcor/nc_xcor_limber_kernel.h>
#include <numcosmo/math/ncm_powspec.h>

G_BEGIN_DECLS

#define NC_TYPE_XCOR (nc_xcor_get_type ())
#define NC_XCOR(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_XCOR, NcXcor))
#define NC_XCOR_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_XCOR, NcXcorClass))
#define NC_IS_XCOR(obj) (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_XCOR))
#define NC_IS_XCOR_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_XCOR))
#define NC_XCOR_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_XCOR, NcXcorClass))

/**
 * NcXcorLimberMethod:
 * @NC_XCOR_LIMBER_METHOD_GSL: FIXME
 * @NC_XCOR_LIMBER_METHOD_CVODE: FIXME
 * @NC_XCOR_LIMBER_METHOD_SUAVE: FIXME
 *
 * FIXME
 *
 *
 */
typedef enum _NcXcorLimberMethod
{
  NC_XCOR_LIMBER_METHOD_GSL = 0,
  NC_XCOR_LIMBER_METHOD_CVODE,
  NC_XCOR_LIMBER_METHOD_SUAVE,
} NcXcorLimberMethod;

#define NC_XCOR_PRECISION (1e-5)

typedef struct _NcXcorClass NcXcorClass;
typedef struct _NcXcor NcXcor;

struct _NcXcor
{
  /*< private > */
  GObject parent_instance;
  NcDistance *dist;
  NcmPowspec *ps;
  gdouble RH;
  NcXcorLimberMethod meth;
};

struct _NcXcorClass
{
  /*< private > */
  GObjectClass parent_class;
  
  gpointer (*alloc) (void);
};

typedef struct _NcXcorKinetic
{
  gdouble xi_z;
  gdouble E_z;
} NcXcorKinetic;

GType nc_xcor_get_type (void) G_GNUC_CONST;
GType nc_xcor_kinetic_get_type (void) G_GNUC_CONST;

NcXcor *nc_xcor_new (NcDistance *dist, NcmPowspec *ps, NcXcorLimberMethod meth);
NcXcor *nc_xcor_ref (NcXcor *xc);
void nc_xcor_free (NcXcor *xc);
void nc_xcor_clear (NcXcor **xc);

NcXcorKinetic *nc_xcor_kinetic_copy (NcXcorKinetic *xck);
void nc_xcor_kinetic_free (NcXcorKinetic *xck);

void nc_xcor_prepare (NcXcor *xc, NcHICosmo *cosmo);

void nc_xcor_limber (NcXcor *xc, NcXcorLimberKernel *xclk1, NcXcorLimberKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, NcmVector *vp);

G_END_DECLS

#endif /* _NC_XCOR_H_ */

