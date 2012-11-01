/***************************************************************************
 *            nc_hicosmo_qpw.h
 *
 *  Mon Aug 11 19:52:57 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
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

#ifndef _NC_HICOSMO_QPW_H_
#define _NC_HICOSMO_QPW_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/likelihood/likelihood.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_QPW             (nc_hicosmo_qpw_get_type ())
#define NC_HICOSMO_QPW(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_QPW, NcHICosmoQPW))
#define NC_HICOSMO_QPW_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_QPW, NcHICosmoQPWClass))
#define NC_IS_HICOSMO_QPW(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_QPW))
#define NC_IS_HICOSMO_QPW_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_QPW))
#define NC_HICOSMO_QPW_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_QPW, NcHICosmoQPWClass))

typedef struct _NcHICosmoQPWClass NcHICosmoQPWClass;
typedef struct _NcHICosmoQPW NcHICosmoQPW;

/**
 * NcHICosmoQPWSParams:
 * @NC_HICOSMO_QPW_H0: FIXME
 * @NC_HICOSMO_QPW_OMEGA_T: FIXME
 * @NC_HICOSMO_QPW_Q0: FIXME
 *
 */
typedef enum _NcHICosmoQPWSParams
{
  NC_HICOSMO_QPW_H0 = 0,
  NC_HICOSMO_QPW_OMEGA_T,
  NC_HICOSMO_QPW_Q0,         /*< private >*/
  NC_HICOSMO_QPW_SPARAM_LEN, /*< skip >*/
} NcHICosmoQPWSParams;

/**
 * NcHICosmoQPWVParams:
 * @NC_HICOSMO_QPW_QP: FIXME
 *
 */
typedef enum _NcHICosmoQPWVParams
{
  NC_HICOSMO_QPW_QP,         /*< private >*/
  NC_HICOSMO_QPW_VPARAM_LEN, /*< skip >*/
} NcHICosmoQPWVParams;

#define NC_HICOSMO_QPW_DEFAULT_H0      NC_C_HUBBLE_CTE_WMAP
#define NC_HICOSMO_QPW_DEFAULT_OMEGA_T ( 1.0)
#define NC_HICOSMO_QPW_DEFAULT_Q0      (-0.5)
#define NC_HICOSMO_QPW_DEFAULT_QP      ( 0.0)
#define NC_HICOSMO_QPW_DEFAULT_QP_LEN     (4)
struct _NcHICosmoQPWClass
{
  /*< private >*/
  NcHICosmoClass parent_class;
};

struct _NcHICosmoQPW
{
  /*< private >*/
  NcHICosmo parent_instance;
  guint npieces;
  guint size;
  gdouble z_f;
  gdouble piece;
  gboolean flat;
};

typedef struct _NcHICosmoQPWContPrior NcHICosmoQPWContPrior;
typedef struct _NcHICosmoQPWAsymCDMPrior NcHICosmoQPWAsymCDMPrior;

struct _NcHICosmoQPWContPrior
{
  /*< private >*/
  gint knot;
  gdouble sigma;
};

struct _NcHICosmoQPWAsymCDMPrior
{
  /*< private >*/
  gdouble z;
  gdouble sigma;
  gdouble q;
};

GType nc_hicosmo_qpw_get_type (void) G_GNUC_CONST;

NcHICosmoQPW *nc_hicosmo_qpw_new (guint npieces, gdouble z_f, gboolean flat);
void nc_hicosmo_qpw_add_continuity_prior (NcHICosmoQPW *qpw, NcLikelihood *lh, gint knot, gdouble sigma);
void nc_hicosmo_qpw_add_continuity_priors (NcHICosmoQPW *qpw, NcLikelihood *lh, gdouble sigma);
void nc_hicosmo_qpw_add_asymptotic_cdm_prior (NcHICosmoQPW *qpw, NcLikelihood *lh, gdouble z, gdouble q, gdouble sigma);

void nc_hicosmo_qpw_change_params (NcHICosmoQPW *qpw, gdouble z);
void nc_hicosmo_qpw_change_params_qpp (NcHICosmoQPW *qpw);

guint nc_hicosmo_qpw_index (NcHICosmoQPW *qpw, gdouble z);

G_END_DECLS

#endif /* _NC_HICOSMO_QPW_H_ */
