/***************************************************************************
 *            nc_halo_bias_type_st_ellip.h
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

#ifndef _NC_HALO_BIAS_TYPE_ST_ELLIP_H_
#define _NC_HALO_BIAS_TYPE_ST_ELLIP_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/lss/nc_halo_bias_type.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_BIAS_TYPE_ST_ELLIP             (nc_halo_bias_type_st_ellip_get_type ())
#define NC_HALO_BIAS_TYPE_ST_ELLIP(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HALO_BIAS_TYPE_ST_ELLIP, NcHaloBiasTypeSTEllip))
#define NC_HALO_BIAS_TYPE_ST_ELLIP_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HALO_BIAS_TYPE_ST_ELLIP, NcHaloBiasTypeSTEllipClass))
#define NC_IS_HALO_BIAS_TYPE_ST_ELLIP(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HALO_BIAS_TYPE_ST_ELLIP))
#define NC_IS_HALO_BIAS_TYPE_ST_ELLIP_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HALO_BIAS_TYPE_ST_ELLIP))
#define NC_HALO_BIAS_TYPE_ST_ELLIP_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HALO_BIAS_TYPE_ST_ELLIP, NcHaloBiasTypeSTEllipClass))

typedef struct _NcHaloBiasTypeSTEllipClass NcHaloBiasTypeSTEllipClass;
typedef struct _NcHaloBiasTypeSTEllip NcHaloBiasTypeSTEllip;

struct _NcHaloBiasTypeSTEllipClass
{
  /*< private >*/
  NcHaloBiasTypeClass parent_class;
};

struct _NcHaloBiasTypeSTEllip
{
  /*< private >*/
  NcHaloBiasType parent_instance;
  gdouble delta_c;
  gdouble a;
  gdouble b;
  gdouble c;
};

GType nc_halo_bias_type_st_ellip_get_type (void) G_GNUC_CONST;

NcHaloBiasType *nc_halo_bias_type_st_ellip_new (gdouble delta_c, gdouble a, gdouble b, gdouble c);
void nc_halo_bias_type_st_ellip_set_delta_c (NcHaloBiasTypeSTEllip *biasf_st_ellip, gdouble delta_c);
gdouble nc_halo_bias_type_st_ellip_get_delta_c (const NcHaloBiasTypeSTEllip *biasf_st_ellip);
void nc_halo_bias_type_st_ellip_set_a (NcHaloBiasTypeSTEllip *biasf_st_ellip, gdouble a);
gdouble nc_halo_bias_type_st_ellip_get_a (const NcHaloBiasTypeSTEllip *biasf_st_ellip);
void nc_halo_bias_type_st_ellip_set_b (NcHaloBiasTypeSTEllip *biasf_st_ellip, gdouble b);
gdouble nc_halo_bias_type_st_ellip_get_b (const NcHaloBiasTypeSTEllip *biasf_st_ellip);
void nc_halo_bias_type_st_ellip_set_c (NcHaloBiasTypeSTEllip *biasf_st_ellip, gdouble c);
gdouble nc_halo_bias_type_st_ellip_get_c (const NcHaloBiasTypeSTEllip *biasf_st_ellip);

G_END_DECLS

#endif /* _NC_HALO_BIAS_TYPE_ST_ELLIP_H_ */
