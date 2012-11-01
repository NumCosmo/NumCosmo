/***************************************************************************
 *            nc_hicosmo_qspline.h
 *
 *  Wed February 15 11:31:28 2012
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

#ifndef _NC_HICOSMO_QSPLINE_H_
#define _NC_HICOSMO_QSPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_ode_spline.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_QSPLINE             (nc_hicosmo_qspline_get_type ())
#define NC_HICOSMO_QSPLINE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_QSPLINE, NcHICosmoQSpline))
#define NC_HICOSMO_QSPLINE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_QSPLINE, NcHICosmoQSplineClass))
#define NC_IS_HICOSMO_QSPLINE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_QSPLINE))
#define NC_IS_HICOSMO_QSPLINE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_QSPLINE))
#define NC_HICOSMO_QSPLINE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_QSPLINE, NcHICosmoQSplineClass))

typedef struct _NcHICosmoQSplineClass NcHICosmoQSplineClass;
typedef struct _NcHICosmoQSpline NcHICosmoQSpline;

/**
 * NcHICosmoQSplineSParams:
 * @NC_HICOSMO_QSPLINE_H0: FIXME
 * @NC_HICOSMO_QSPLINE_OMEGA_T: FIXME
 *
 */
typedef enum _NcHICosmoQSplineSParams
{
  NC_HICOSMO_QSPLINE_H0 = 0,
  NC_HICOSMO_QSPLINE_OMEGA_T,    /*< private >*/
  NC_HICOSMO_QSPLINE_SPARAM_LEN, /*< skip >*/
} NcHICosmoQSplineSParams;

/**
 * NcHICosmoQSplineVParams:
 * @NC_HICOSMO_QSPLINE_Q: FIXME
 *
 */
typedef enum _NcHICosmoQSplineVParams
{
  NC_HICOSMO_QSPLINE_Q,          /*< private >*/
  NC_HICOSMO_QSPLINE_VPARAM_LEN, /*< skip >*/
} NcHICosmoQSplineVParams;

#define NC_HICOSMO_QSPLINE_DEFAULT_H0      NC_C_HUBBLE_CTE_WMAP
#define NC_HICOSMO_QSPLINE_DEFAULT_OMEGA_T ( 1.0)
#define NC_HICOSMO_QSPLINE_DEFAULT_Q       (-0.5)
#define NC_HICOSMO_QSPLINE_DEFAULT_Q_LEN      (6)

struct _NcHICosmoQSplineClass
{
  /*< private >*/
  NcHICosmoClass parent_class;
};

struct _NcHICosmoQSpline
{
  /*< private >*/
  NcHICosmo parent_instance;
  guint nknots;
  guint size;
  gdouble z_f;
  gsize pkey;
  NcmSpline *q_z;
  NcmOdeSpline *E2_z;
};

GType nc_hicosmo_qspline_get_type (void) G_GNUC_CONST;

NcHICosmoQSpline *nc_hicosmo_qspline_new (NcmSpline *s, gsize np, gdouble z_f);

G_END_DECLS

#endif /* _NC_HICOSMO_QSPLINE_H_ */
