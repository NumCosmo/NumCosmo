/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_grav.h
 *
 *  Thu October 12 14:32:08 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_grav.h
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

#ifndef _NC_HIPERT_GRAV_H_
#define _NC_HIPERT_GRAV_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/perturbations/nc_hipert_bg_var.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_GRAV             (nc_hipert_grav_get_type ())
#define NC_HIPERT_GRAV(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_GRAV, NcHIPertGrav))
#define NC_HIPERT_GRAV_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_GRAV, NcHIPertGravClass))
#define NC_IS_HIPERT_GRAV(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_GRAV))
#define NC_IS_HIPERT_GRAV_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_GRAV))
#define NC_HIPERT_GRAV_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_GRAV, NcHIPertGravClass))

typedef struct _NcHIPertGravClass NcHIPertGravClass;
typedef struct _NcHIPertGrav NcHIPertGrav;
typedef struct _NcHIPertGravPrivate NcHIPertGravPrivate;

typedef struct _NcHIPertGravScalar
{
  gdouble phi;
  gdouble dsigma;
  gdouble psi;
} NcHIPertGravScalar;

typedef struct _NcHIPertGravVector
{
  gdouble dsigma[2];
} NcHIPertGravVector;

typedef struct _NcHIPertGravTensor
{
  gdouble h[2];
} NcHIPertGravTensor;

/**
 * NcHIPertCompGauge:
 * @NC_HIPERT_GRAV_GAUGE_SYNCHRONOUS: Synchronous gauge
 * @NC_HIPERT_GRAV_GAUGE_NEWTONIAN: Newtonian gauge
 * @NC_HIPERT_GRAV_GAUGE_CONST_CURV: Constant curvature gauge
 * 
 * FIXME
 * 
 */
typedef enum /*< enum,underscore_name=NC_HIPERT_GRAV_GAUGE >*/
{
  NC_HIPERT_GRAV_GAUGE_SYNCHRONOUS,
  NC_HIPERT_GRAV_GAUGE_NEWTONIAN,
  NC_HIPERT_GRAV_GAUGE_CONST_CURV,
  /* < private > */
  NC_HIPERT_GRAV_GAUGE_LEN, /*< skip >*/
} NcHIPertCompGauge;

typedef void (*NcHIPertGravSetGauge) (NcHIPertGrav *grav, NcHIPertCompGauge gauge);
typedef NcHIPertCompGauge (*NcHIPertGravGetGauge) (NcHIPertGrav *grav);

struct _NcHIPertGravClass
{
  /*< private >*/
  GObjectClass parent_class;
  NcHIPertGravSetGauge set_gauge;
  NcHIPertGravGetGauge get_gauge;
};

struct _NcHIPertGrav
{
  /*< private >*/
  GObject parent_instance;
  NcHIPertGravPrivate *priv;
};

GType nc_hipert_grav_scalar_get_type (void) G_GNUC_CONST;
GType nc_hipert_grav_vector_get_type (void) G_GNUC_CONST;
GType nc_hipert_grav_tensor_get_type (void) G_GNUC_CONST;
GType nc_hipert_grav_get_type (void) G_GNUC_CONST;

G_INLINE_FUNC NcHIPertGravScalar *nc_hipert_grav_scalar_new (void);
G_INLINE_FUNC NcHIPertGravScalar *nc_hipert_grav_scalar_dup (NcHIPertGravScalar *gs);
G_INLINE_FUNC void nc_hipert_grav_scalar_free (NcHIPertGravScalar *gs);
G_INLINE_FUNC void nc_hipert_grav_scalar_set_zero (NcHIPertGravScalar *gs);

G_INLINE_FUNC NcHIPertGravVector *nc_hipert_grav_vector_new (void);
G_INLINE_FUNC NcHIPertGravVector *nc_hipert_grav_vector_dup (NcHIPertGravVector *gv);
G_INLINE_FUNC void nc_hipert_grav_vector_free (NcHIPertGravVector *gv);
G_INLINE_FUNC void nc_hipert_grav_vector_set_zero (NcHIPertGravVector *gv);

G_INLINE_FUNC NcHIPertGravTensor *nc_hipert_grav_tensor_new (void);
G_INLINE_FUNC NcHIPertGravTensor *nc_hipert_grav_tensor_dup (NcHIPertGravTensor *gt);
G_INLINE_FUNC void nc_hipert_grav_tensor_free (NcHIPertGravTensor *gt);
G_INLINE_FUNC void nc_hipert_grav_tensor_set_zero (NcHIPertGravTensor *gt);

NcHIPertGrav *nc_hipert_grav_ref (NcHIPertGrav *grav);
void nc_hipert_grav_free (NcHIPertGrav *grav);
void nc_hipert_grav_clear (NcHIPertGrav **grav);

G_INLINE_FUNC NcHIPertBGVarID nc_hipert_grav_get_id (NcHIPertGrav *grav);

void nc_hipert_grav_set_gauge (NcHIPertGrav *grav, NcHIPertCompGauge gauge);
NcHIPertCompGauge nc_hipert_grav_get_gauge (NcHIPertGrav *grav);

G_END_DECLS

#endif /* _NC_HIPERT_GRAV_H_ */

#ifndef _NC_HIPERT_GRAV_INLINE_H_
#define _NC_HIPERT_GRAV_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC NcHIPertGravScalar *
nc_hipert_grav_scalar_new (void)
{
  return g_new0 (NcHIPertGravScalar, 1);
}

G_INLINE_FUNC NcHIPertGravScalar *
nc_hipert_grav_scalar_dup (NcHIPertGravScalar *gs)
{
  NcHIPertGravScalar *gs_dup = g_new0 (NcHIPertGravScalar, 1);

  gs_dup[0] = gs[0];

  return gs_dup;
}

G_INLINE_FUNC void
nc_hipert_grav_scalar_free (NcHIPertGravScalar *gs)
{
  g_free (gs);
}

G_INLINE_FUNC void
nc_hipert_grav_scalar_set_zero (NcHIPertGravScalar *gs)
{
  gs->phi    = 0.0;
  gs->dsigma = 0.0;
  gs->psi    = 0.0;
}

G_INLINE_FUNC NcHIPertGravVector *
nc_hipert_grav_vector_new (void)
{
  return g_new0 (NcHIPertGravVector, 1);
}

G_INLINE_FUNC NcHIPertGravVector *
nc_hipert_grav_vector_dup (NcHIPertGravVector *gv)
{
  NcHIPertGravVector *gv_dup = g_new0 (NcHIPertGravVector, 1);

  gv_dup[0] = gv[0];

  return gv_dup;
}

G_INLINE_FUNC void
nc_hipert_grav_vector_free (NcHIPertGravVector *gv)
{
  g_free (gv);
}

G_INLINE_FUNC void
nc_hipert_grav_vector_set_zero (NcHIPertGravVector *gv)
{
  gv->dsigma[0] = 0.0;
  gv->dsigma[1] = 0.0;
}

G_INLINE_FUNC NcHIPertGravTensor *
nc_hipert_grav_tensor_new (void)
{
  return g_new0 (NcHIPertGravTensor, 1);
}

G_INLINE_FUNC NcHIPertGravTensor *
nc_hipert_grav_tensor_dup (NcHIPertGravTensor *gt)
{
  NcHIPertGravTensor *gt_dup = g_new0 (NcHIPertGravTensor, 1);

  gt_dup[0] = gt[0];

  return gt_dup;
}

G_INLINE_FUNC void
nc_hipert_grav_tensor_free (NcHIPertGravTensor *gt)
{
  g_free (gt);
}

G_INLINE_FUNC void
nc_hipert_grav_tensor_set_zero (NcHIPertGravTensor *gt)
{
  gt->h[0] = 0.0;
  gt->h[1] = 0.0;
}

G_INLINE_FUNC NcHIPertBGVarID 
nc_hipert_grav_get_id (NcHIPertGrav *grav)
{
  const NcHIPertBGVarID id = nc_hipert_bg_var_class_get_id_by_ns (G_OBJECT_TYPE_NAME (grav));
  g_assert_cmpint (id, >=, 0);
  return id;
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPERT_GRAV_INLINE_H_ */
