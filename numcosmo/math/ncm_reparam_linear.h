/***************************************************************************
 *            ncm_reparam_linear.h
 *
 *  Thu March 08 11:05:07 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NCM_REPARAM_LINEAR_H_
#define _NCM_REPARAM_LINEAR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_reparam.h>
#include <numcosmo/math/ncm_model.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_permutation.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_REPARAM_LINEAR             (ncm_reparam_linear_get_type ())
#define NCM_REPARAM_LINEAR(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_REPARAM_LINEAR, NcmReparamLinear))
#define NCM_REPARAM_LINEAR_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_REPARAM_LINEAR, NcmReparamLinearClass))
#define NCM_IS_REPARAM_LINEAR(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_REPARAM_LINEAR))
#define NCM_IS_REPARAM_LINEAR_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_REPARAM_LINEAR))
#define NCM_REPARAM_LINEAR_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_REPARAM_LINEAR, NcmReparamLinearClass))

typedef struct _NcmReparamLinearClass NcmReparamLinearClass;
typedef struct _NcmReparamLinear NcmReparamLinear;

struct _NcmReparamLinearClass
{
  /*< private >*/
  NcmReparamClass parent_class;
};

struct _NcmReparamLinear
{
  /*< private >*/
  NcmReparam parent_instance;
  NcmMatrix *T;
  NcmVector *v;
  NcmVector *vp;
  NcmMatrix *T_LU;
  gsl_permutation *p;
  gint signum;
};

GType ncm_reparam_linear_get_type (void) G_GNUC_CONST;

NcmReparamLinear *ncm_reparam_linear_new (guint size, NcmMatrix *T, NcmVector *v);
void ncm_reparam_linear_set_compat_type (NcmReparamLinear *lin, GType compat_type);

G_END_DECLS

#endif /* _NCM_REPARAM_LINEAR_H_ */
