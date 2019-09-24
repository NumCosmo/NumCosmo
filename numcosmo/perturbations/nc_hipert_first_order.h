/***************************************************************************
 *            nc_hipert_first_order.h
 *
 *  Mon October 09 16:58:26 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_first_order.h
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_HIPERT_FIRST_ORDER_H_
#define _NC_HIPERT_FIRST_ORDER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/perturbations/nc_hipert_boltzmann.h>
#include <numcosmo/perturbations/nc_hipert_grav.h>
#include <numcosmo/perturbations/nc_hipert_comp.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_FIRST_ORDER             (nc_hipert_first_order_get_type ())
#define NC_HIPERT_FIRST_ORDER(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_FIRST_ORDER, NcHIPertFirstOrder))
#define NC_HIPERT_FIRST_ORDER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_FIRST_ORDER, NcHIPertFirstOrderClass))
#define NC_IS_HIPERT_FIRST_ORDER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_FIRST_ORDER))
#define NC_IS_HIPERT_FIRST_ORDER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_FIRST_ORDER))
#define NC_HIPERT_FIRST_ORDER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_FIRST_ORDER, NcHIPertFirstOrderClass))

typedef struct _NcHIPertFirstOrderClass NcHIPertFirstOrderClass;
typedef struct _NcHIPertFirstOrder NcHIPertFirstOrder;
typedef struct _NcHIPertFirstOrderPrivate NcHIPertFirstOrderPrivate;

struct _NcHIPertFirstOrderClass
{
  /*< private >*/
  NcHIPertBoltzmannClass parent_class;
};

struct _NcHIPertFirstOrder
{
  /*< private >*/
  NcHIPertBoltzmann parent_instance;
  NcHIPertFirstOrderPrivate *priv;
};

/**
 * NcHIPertFirstOrderInteg:
 * @NC_HIPERT_FIRST_ORDER_INTEG_CVODE: CVode integrator
 * @NC_HIPERT_FIRST_ORDER_INTEG_ARKODE: ARKode integrator
 * 
 * ODE integrators.
 * 
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_FIRST_ORDER_INTEG >*/
{
  NC_HIPERT_FIRST_ORDER_INTEG_CVODE,
  NC_HIPERT_FIRST_ORDER_INTEG_ARKODE,
  /* < private > */
  NC_HIPERT_FIRST_ORDER_INTEG_LEN, /*< skip >*/
} NcHIPertFirstOrderInteg;

GType nc_hipert_first_order_get_type (void) G_GNUC_CONST;

NcHIPertFirstOrder *nc_hipert_first_order_new (void);
NcHIPertFirstOrder *nc_hipert_first_order_new_full (NcDistance *dist, NcRecomb *recomb, NcScalefactor *a);
NcHIPertFirstOrder *nc_hipert_first_order_ref (NcHIPertFirstOrder *fo);

void nc_hipert_first_order_free (NcHIPertFirstOrder *fo);
void nc_hipert_first_order_clear (NcHIPertFirstOrder **fo);

void nc_hipert_first_order_set_gauge (NcHIPertFirstOrder *fo, NcHIPertGravGauge gauge);
NcHIPertGravGauge nc_hipert_first_order_get_gauge (NcHIPertFirstOrder *fo);

void nc_hipert_first_order_set_reltol (NcHIPertFirstOrder *fo, const gdouble reltol);
void nc_hipert_first_order_set_abstol (NcHIPertFirstOrder *fo, const gdouble abstol);
gdouble nc_hipert_first_order_get_reltol (NcHIPertFirstOrder *fo);
gdouble nc_hipert_first_order_get_abstol (NcHIPertFirstOrder *fo);

void nc_hipert_first_order_set_integ (NcHIPertFirstOrder *fo, NcHIPertFirstOrderInteg integ);
NcHIPertFirstOrderInteg nc_hipert_first_order_get_integ (NcHIPertFirstOrder *fo);

void nc_hipert_first_order_set_grav (NcHIPertFirstOrder *fo, NcHIPertGrav *grav);
NcHIPertGrav *nc_hipert_first_order_get_grav (NcHIPertFirstOrder *fo);
NcHIPertGrav *nc_hipert_first_order_peek_grav (NcHIPertFirstOrder *fo);

void nc_hipert_first_order_add_comp (NcHIPertFirstOrder *fo, NcHIPertComp *comp);

void nc_hipert_first_order_prepare (NcHIPertFirstOrder *fo, NcHICosmo *cosmo);

G_END_DECLS

#endif /* _NC_HIPERT_FIRST_ORDER_H_ */
