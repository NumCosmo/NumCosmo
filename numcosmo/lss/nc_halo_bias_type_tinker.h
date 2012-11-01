/***************************************************************************
 *            nc_halo_bias_type_tinker.h
 *
 *  Tue June 28 15:41:57 2011
 *  Copyright  2011  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_HALO_BIAS_TYPE_TINKER_H_
#define _NC_HALO_BIAS_TYPE_TINKER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/lss/nc_halo_bias_type.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_BIAS_TYPE_TINKER             (nc_halo_bias_type_tinker_get_type ())
#define NC_HALO_BIAS_TYPE_TINKER(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HALO_BIAS_TYPE_TINKER, NcHaloBiasTypeTinker))
#define NC_HALO_BIAS_TYPE_TINKER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HALO_BIAS_TYPE_TINKER, NcHaloBiasTypeTinkerClass))
#define NC_IS_HALO_BIAS_TYPE_TINKER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HALO_BIAS_TYPE_TINKER))
#define NC_IS_HALO_BIAS_TYPE_TINKER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HALO_BIAS_TYPE_TINKER))
#define NC_HALO_BIAS_TYPE_TINKER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HALO_BIAS_TYPE_TINKER, NcHaloBiasTypeTinkerClass))

typedef struct _NcHaloBiasTypeTinkerClass NcHaloBiasTypeTinkerClass;
typedef struct _NcHaloBiasTypeTinker NcHaloBiasTypeTinker;



struct _NcHaloBiasTypeTinkerClass
{
  /*< private >*/
  NcHaloBiasTypeClass parent_class;
};

struct _NcHaloBiasTypeTinker
{
  /*< private >*/
  NcHaloBiasType parent_instance;
  gdouble delta_c;
  gdouble B;
  gdouble b;
  gdouble c;
  gdouble Delta;   
};

GType nc_halo_bias_type_tinker_get_type (void) G_GNUC_CONST;

NcHaloBiasType *nc_halo_bias_type_tinker_new (gdouble delta_c, gdouble B, gdouble b, gdouble c, gdouble Delta);
void nc_halo_bias_type_tinker_set_delta_c (NcHaloBiasTypeTinker *biasf_tinker, gdouble delta_c);
gdouble nc_halo_bias_type_tinker_get_delta_c (const NcHaloBiasTypeTinker *biasf_tinker);
void nc_halo_bias_type_tinker_set_B (NcHaloBiasTypeTinker *biasf_tinker, gdouble B);
gdouble nc_halo_bias_type_tinker_get_B (const NcHaloBiasTypeTinker *biasf_tinker);
void nc_halo_bias_type_tinker_set_b (NcHaloBiasTypeTinker *biasf_tinker, gdouble b);
gdouble nc_halo_bias_type_tinker_get_b (const NcHaloBiasTypeTinker *biasf_tinker);
void nc_halo_bias_type_tinker_set_c (NcHaloBiasTypeTinker *biasf_tinker, gdouble c);
gdouble nc_halo_bias_type_tinker_get_c (const NcHaloBiasTypeTinker *biasf_tinker);
void nc_halo_bias_type_tinker_set_Delta (NcHaloBiasTypeTinker *biasf_tinker, gdouble Delta);
gdouble nc_halo_bias_type_tinker_get_Delta (const NcHaloBiasTypeTinker *biasf_tinker);

G_END_DECLS

#endif /* _NC_HALO_BIAS_TYPE_TINKER_H_ */
