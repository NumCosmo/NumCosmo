/***************************************************************************
 *            nc_hicosmo_de_wpoly.h
 *
 *  Mon Sep 28 10:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

#ifndef _NC_HICOSMO_DE_WPOLY_H_
#define _NC_HICOSMO_DE_WPOLY_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc/background/nc_hicosmo_de.h>

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_DE_WPOLY             (nc_hicosmo_de_wpoly_get_type ())
#define NC_HICOSMO_DE_WPOLY(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_DE_WPOLY, NcHICosmoDEWPoly))
#define NC_HICOSMO_DE_WPOLY_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_DE_WPOLY, NcHICosmoDEWPolyClass))
#define NC_IS_HICOSMO_DE_WPOLY(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_DE_WPOLY))
#define NC_IS_HICOSMO_DE_WPOLY_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_DE_WPOLY))
#define NC_HICOSMO_DE_WPOLY_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_DE_WPOLY, NcHICosmoDEWPolyClass))

typedef struct _NcHICosmoDEWPolyClass NcHICosmoDEWPolyClass;
typedef struct _NcHICosmoDEWPoly NcHICosmoDEWPoly;
typedef struct _NcHICosmoDEWPolyPrivate NcHICosmoDEWPolyPrivate;

struct _NcHICosmoDEWPolyClass
{
  /*< private >*/
  NcHICosmoDEClass parent_class;
};

struct _NcHICosmoDEWPoly
{
  /*< private >*/
  NcHICosmoDE parent_instance;
  NcHICosmoDEWPolyPrivate *priv;
};

GType nc_hicosmo_de_wpoly_get_type (void) G_GNUC_CONST;

/**
 * NcHICosmoDEWPolySParams:
 * @NC_HICOSMO_DE_WPOLY_AS_DRAG: Sound horizon at the drag epoch in units of the Hubble distance, $r_s(z_d) H_0 / c$, used when #NcHICosmoDEWPoly:use-as-drag is set
 *
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_DE_WPOLY_SPARAMS,prefix=NC_HICOSMO_DE_WPOLY >*/
{
  NC_HICOSMO_DE_WPOLY_AS_DRAG = NC_HICOSMO_DE_SPARAM_LEN,
  /* < private > */
  NC_HICOSMO_DE_WPOLY_SPARAM_LEN, /*< skip >*/
} NcHICosmoDEWPolySParams;

/**
 * NcHICosmoDEWPolyVParams:
 * @NC_HICOSMO_DE_WPOLY_A: Legendre coefficients $a_\ell$ of $w(z) + 1$ on $[0, z_f]$
 *
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_DE_WPOLY_VPARAMS,prefix=NC_HICOSMO_DE_WPOLY >*/
{
  NC_HICOSMO_DE_WPOLY_A = NC_HICOSMO_DE_VPARAM_LEN,
  /* < private > */
  NC_HICOSMO_DE_WPOLY_VPARAM_LEN, /*< skip >*/
} NcHICosmoDEWPolyVParams;

#define NC_HICOSMO_DE_WPOLY_DEFAULT_A (0.0)
#define NC_HICOSMO_DE_WPOLY_DEFAULT_A_LEN (8)
#define NC_HICOSMO_DE_WPOLY_DEFAULT_AS_DRAG (0.035)

NcHICosmoDEWPoly *nc_hicosmo_de_wpoly_new (gsize ncoeffs, const gdouble z_f);

gdouble nc_hicosmo_de_wpoly_get_z_f (NcHICosmoDEWPoly *wpoly);
guint nc_hicosmo_de_wpoly_get_ncoeffs (NcHICosmoDEWPoly *wpoly);
gdouble nc_hicosmo_de_wpoly_dw_dz (NcHICosmoDEWPoly *wpoly, const gdouble z);
GArray *nc_hicosmo_de_wpoly_turning_points (NcHICosmoDEWPoly *wpoly);

void nc_hicosmo_de_wpoly_set_use_as_drag (NcHICosmoDEWPoly *wpoly, gboolean use_as_drag);
gboolean nc_hicosmo_de_wpoly_get_use_as_drag (NcHICosmoDEWPoly *wpoly);

G_END_DECLS

#endif /* _NC_HICOSMO_DE_WPOLY_H_ */

