/***************************************************************************
 *            nc_hiprim_atan.h
 *
 *  Thu October 29 15:15:24 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hiprim_atan.h
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_HIPRIM_ATAN_H_
#define _NC_HIPRIM_ATAN_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_c.h>
#include <numcosmo/nc_hiprim.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPRIM_ATAN             (nc_hiprim_atan_get_type ())
#define NC_HIPRIM_ATAN(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPRIM_ATAN, NcHIPrimAtan))
#define NC_HIPRIM_ATAN_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPRIM_ATAN, NcHIPrimAtanClass))
#define NC_IS_HIPRIM_ATAN(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPRIM_ATAN))
#define NC_IS_HIPRIM_ATAN_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPRIM_ATAN))
#define NC_HIPRIM_ATAN_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPRIM_ATAN, NcHIPrimAtanClass))

typedef struct _NcHIPrimAtanClass NcHIPrimAtanClass;
typedef struct _NcHIPrimAtan NcHIPrimAtan;

/**
 * NcHIPrimAtanParams:
 * @NC_HIPRIM_ATAN_LN10E10ASA: Amplitude of the adiabatic scalar mode $\ln(10^10A_{SA})$
 * @NC_HIPRIM_ATAN_N_SA: Adiabatic scalar spectral index
 * @NC_HIPRIM_ATAN_LNKC: Arctan parameter $\ln(k_c)$
 * @NC_HIPRIM_ATAN_C2: Arctan parameter $c_2$
 * @NC_HIPRIM_ATAN_C3: Arctan parameter $c_3$
 *
 * FIXME
 * 
 */
typedef enum _NcHIPrimAtanParams
{
  NC_HIPRIM_ATAN_LN10E10ASA,
  NC_HIPRIM_ATAN_N_SA,
  NC_HIPRIM_ATAN_LNKC,
  NC_HIPRIM_ATAN_C2,
  NC_HIPRIM_ATAN_C3,         /*< private >*/
  NC_HIPRIM_ATAN_SPARAM_LEN, /*< skip >*/
} NcHIPrimAtanParams;

struct _NcHIPrimAtanClass
{
  /*< private >*/
  NcHIPrimClass parent_class;
};

struct _NcHIPrimAtan
{
  /*< private >*/
  NcHIPrim parent_instance;
};

GType nc_hiprim_atan_get_type (void) G_GNUC_CONST;

NcHIPrimAtan *nc_hiprim_atan_new (void);

#define NC_HIPRIM_ATAN_DEFAULT_LN10E10ASA (3.179)
#define NC_HIPRIM_ATAN_DEFAULT_N_SA (0.9742)
#define NC_HIPRIM_ATAN_DEFAULT_LNKC (-5.3)
#define NC_HIPRIM_ATAN_DEFAULT_C2 (1.0)
#define NC_HIPRIM_ATAN_DEFAULT_C3 (1.0)

G_END_DECLS

#endif /* _NC_HIPRIM_ATAN_H_ */
