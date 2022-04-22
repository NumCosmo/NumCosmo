/***************************************************************************
 *            nc_halo_bias_func.h
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_HALO_BIAS_FUNC_H_
#define _NC_HALO_BIAS_FUNC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_halo_bias_type.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/lss/nc_halo_mass_function.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_BIAS_FUNC             (nc_halo_bias_func_get_type ())
#define NC_HALO_BIAS_FUNC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HALO_BIAS_FUNC, NcHaloBiasFunc))
#define NC_HALO_BIAS_FUNC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HALO_BIAS_FUNC, NcHaloBiasFuncClass))
#define NC_IS_HALO_BIAS_FUNC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HALO_BIAS_FUNC))
#define NC_IS_HALO_BIAS_FUNC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HALO_BIAS_FUNC))
#define NC_HALO_BIAS_FUNC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HALO_BIAS_FUNC, NcHaloBiasFuncClass))

typedef struct _NcHaloBiasFuncClass NcHaloBiasFuncClass;
typedef struct _NcHaloBiasFunc NcHaloBiasFunc;

struct _NcHaloBiasFuncClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcHaloBiasFunc
{
  /*< private >*/
  GObject parent_instance;
  NcHaloMassFunction *mfp;
  NcHaloBiasType *biasf;
};

GType nc_halo_bias_func_get_type (void) G_GNUC_CONST;

NcHaloBiasFunc *nc_halo_bias_func_new (NcHaloMassFunction *mfp, NcHaloBiasType *biasf);
NcHaloBiasFunc *nc_halo_bias_func_copy (NcHaloBiasFunc *mbiasf);
void nc_halo_bias_func_free (NcHaloBiasFunc *mbiasf);
void nc_halo_bias_func_clear (NcHaloBiasFunc **mbiasf);

gdouble nc_halo_bias_func_integrand (NcHaloBiasFunc *mbiasf, NcHICosmo *cosmo, gdouble lnM, gdouble z); 
gdouble nc_halo_bias_func_int (NcHaloBiasFunc *mbiasf, NcHICosmo *cosmo);

G_END_DECLS

#endif /* _NC_HALO_BIAS_FUNC_H_ */
