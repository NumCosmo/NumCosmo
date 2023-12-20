/***************************************************************************
 *            ncm_reparam_linear.h
 *
 *  Thu March 08 11:05:07 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#define NCM_TYPE_REPARAM_LINEAR (ncm_reparam_linear_get_type ())

G_DECLARE_FINAL_TYPE (NcmReparamLinear, ncm_reparam_linear, NCM, REPARAM_LINEAR, NcmReparam)

NcmReparamLinear *ncm_reparam_linear_new (guint size, NcmMatrix * T, NcmVector * v);
void ncm_reparam_linear_set_compat_type (NcmReparamLinear *lin, GType compat_type);

G_END_DECLS

#endif /* _NCM_REPARAM_LINEAR_H_ */

