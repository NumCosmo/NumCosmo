/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_bg_var.h
 *
 *  Fri October 13 15:56:51 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_bg_var.h
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

#ifndef _NC_HIPERT_BG_VAR_H_
#define _NC_HIPERT_BG_VAR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/nc_recomb.h>
#include <numcosmo/nc_scalefactor.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_BG_VAR             (nc_hipert_bg_var_get_type ())
#define NC_HIPERT_BG_VAR(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_BG_VAR, NcHIPertBGVar))
#define NC_HIPERT_BG_VAR_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_BG_VAR, NcHIPertBGVarClass))
#define NC_IS_HIPERT_BG_VAR(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_BG_VAR))
#define NC_IS_HIPERT_BG_VAR_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_BG_VAR))
#define NC_HIPERT_BG_VAR_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_BG_VAR, NcHIPertBGVarClass))

typedef struct _NcHIPertBGVarClass NcHIPertBGVarClass;
typedef struct _NcHIPertBGVar NcHIPertBGVar;
typedef struct _NcHIPertBGVarPrivate NcHIPertBGVarPrivate;
typedef struct _NcHIPertBGVarDesc NcHIPertBGVarDesc;
typedef struct _NcHIPertBGVarYDY NcHIPertBGVarYDY;

typedef gint32 NcHIPertBGVarID;

struct _NcHIPertBGVarDesc
{
  /*< private >*/
  gboolean init;
  gchar *ns;
  gchar *desc;
  gchar *long_desc;
  guint cstruct_size;
};

struct _NcHIPertBGVarClass
{
  /*< private >*/
  GObjectClass parent_class;
  NcHIPertBGVarID bg_var_id_len;
  GHashTable *ns_table;
  GArray *bg_var_desc_array;
};

struct _NcHIPertBGVar
{
  /*< private >*/
  GObject parent_instance;
  NcHIPertBGVarPrivate *priv;
  GPtrArray *cstructs;
  NcRecomb *recomb;
  NcDistance *dist;
  NcScalefactor *a;
  gdouble t;
  gdouble eta;
  gdouble k;
  gdouble x;
  gdouble E;
};

/**
 * NcHIPertBGVarYDY:
 * 
 * Boxed object containing the current status of the ode system.
 */
struct _NcHIPertBGVarYDY
{
  /*< private >*/
  gdouble *y;
  gdouble *dy;
  guint start_index;
  gint *perm;
  gint *perm_inv;
};

GType nc_hipert_bg_var_ydy_get_type (void) G_GNUC_CONST;
GType nc_hipert_bg_var_get_type (void) G_GNUC_CONST;

void nc_hipert_bg_var_class_register_id (const gchar *ns, const gchar *desc, const gchar *long_desc, guint cstruct_size);
NcHIPertBGVarID nc_hipert_bg_var_class_get_id_by_gtype (GType gt);
NcHIPertBGVarID nc_hipert_bg_var_class_get_id_by_ns (const gchar *ns);

/**
 * NC_HIPERT_BG_VAR_ID_FUNC: (skip)
 * @obj_ns: object namespace
 *
 * FIXME
 *
 */
#define NC_HIPERT_BG_VAR_ID_FUNC(obj_ns) obj_ns##_id

/**
 * NC_HIPERT_BG_VAR_ID_FUNC_DECL: (skip)
 * @obj_ns: object namespace
 *
 * Declare the id function associated with @obj_ns.
 *
 */
#define NC_HIPERT_BG_VAR_ID_FUNC_DECL(obj_ns) NcHIPertBGVarID NC_HIPERT_BG_VAR_ID_FUNC(obj_ns) (void) G_GNUC_CONST

/**
 * NC_HIPERT_BG_VAR_ID_FUNC_IMPL: (skip)
 * @obj_ns: object namespace
 * @ns: namespace
 *
 * The implementation of the id function associated with @obj_ns.
 *
 */
#define NC_HIPERT_BG_VAR_ID_FUNC_IMPL(obj_ns,ns) \
NcHIPertBGVarID NC_HIPERT_BG_VAR_ID_FUNC(obj_ns) (void) \
{ \
  static NcHIPertBGVarID id = -1; \
  if (id == -1) \
  { \
    NcHIPertBGVarClass *bg_var_class = g_type_class_ref (NC_TYPE_HIPERT_BG_VAR); \
    id = nc_hipert_bg_var_class_get_id_by_ns (#ns); \
    g_type_class_unref (bg_var_class); \
  } \
  return id; \
}

NCM_INLINE NcHIPertBGVarYDY *nc_hipert_bg_var_ydy_new (void);
NCM_INLINE NcHIPertBGVarYDY *nc_hipert_bg_var_ydy_dup (NcHIPertBGVarYDY *ydy);
NCM_INLINE void nc_hipert_bg_var_ydy_free (NcHIPertBGVarYDY *ydy);

NCM_INLINE gdouble nc_hipert_bg_var_ydy_get_y_i (NcHIPertBGVarYDY *ydy, guint i);
NCM_INLINE void nc_hipert_bg_var_ydy_set_dy_i (NcHIPertBGVarYDY *ydy, guint i, const gdouble dy_i);
NCM_INLINE gdouble nc_hipert_bg_var_ydy_get_dy_i (NcHIPertBGVarYDY *ydy, guint i);

NcHIPertBGVar *nc_hipert_bg_var_new (void);
NcHIPertBGVar *nc_hipert_bg_var_new_full (NcDistance *dist, NcRecomb *recomb, NcScalefactor *a);
NcHIPertBGVar *nc_hipert_bg_var_ref (NcHIPertBGVar *bg_var);

void nc_hipert_bg_var_free (NcHIPertBGVar *bg_var);
void nc_hipert_bg_var_clear (NcHIPertBGVar **bg_var);

void nc_hipert_bg_var_prepare (NcHIPertBGVar *bg_var, NcHICosmo *cosmo);
void nc_hipert_bg_var_prepare_if_needed (NcHIPertBGVar *bg_var, NcHICosmo *cosmo);

NCM_INLINE void nc_hipert_bg_var_set_dist (NcHIPertBGVar *bg_var, NcDistance *dist);
NCM_INLINE void nc_hipert_bg_var_set_recomb (NcHIPertBGVar *bg_var, NcRecomb *recomb);
NCM_INLINE void nc_hipert_bg_var_set_scalefactor (NcHIPertBGVar *bg_var, NcScalefactor *a);

NCM_INLINE NcDistance *nc_hipert_bg_var_get_dist (NcHIPertBGVar *bg_var);
NCM_INLINE NcRecomb *nc_hipert_bg_var_get_recomb (NcHIPertBGVar *bg_var);
NCM_INLINE NcScalefactor *nc_hipert_bg_var_get_scalefactor (NcHIPertBGVar *bg_var);

NCM_INLINE NcDistance *nc_hipert_bg_var_peek_dist (NcHIPertBGVar *bg_var);
NCM_INLINE NcRecomb *nc_hipert_bg_var_peek_recomb (NcHIPertBGVar *bg_var);
NCM_INLINE NcScalefactor *nc_hipert_bg_var_peek_scalefactor (NcHIPertBGVar *bg_var);

void nc_hipert_bg_var_set_zf (NcHIPertBGVar *bg_var, const gdouble zf);
gdouble nc_hipert_bg_var_get_zf (NcHIPertBGVar *bg_var);

guint nc_hipert_bg_var_len (NcHIPertBGVar *bg_var);

void nc_hipert_bg_var_activate_id (NcHIPertBGVar *bg_var, ...);
void nc_hipert_bg_var_activate_id_array (NcHIPertBGVar *bg_var, GArray *ids);

#define NC_HIPERT_BG_VAR_DEFAULT_ZF (1.0e9)

G_END_DECLS

#endif /* _NC_HIPERT_BG_VAR_H_ */

#ifndef _NC_HIPERT_BG_VAR_INLINE_H_
#define _NC_HIPERT_BG_VAR_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE NcHIPertBGVarYDY *
nc_hipert_bg_var_ydy_new (void)
{
  NcHIPertBGVarYDY *ydy = g_new0 (NcHIPertBGVarYDY, 1);
  return ydy;
}

NCM_INLINE NcHIPertBGVarYDY *
nc_hipert_bg_var_ydy_dup (NcHIPertBGVarYDY *ydy)
{
  NcHIPertBGVarYDY *ydy_dup = nc_hipert_bg_var_ydy_new ();

  ydy_dup[0] = ydy[0];

  return ydy_dup;
}

NCM_INLINE void 
nc_hipert_bg_var_ydy_free (NcHIPertBGVarYDY *ydy)
{
  g_free (ydy);
}

NCM_INLINE gdouble 
nc_hipert_bg_var_ydy_get_y_i (NcHIPertBGVarYDY *ydy, guint i)
{
  return ydy->y[ydy->perm_inv[ydy->start_index + i]];
}

NCM_INLINE void 
nc_hipert_bg_var_ydy_set_dy_i (NcHIPertBGVarYDY *ydy, guint i, const gdouble dy_i)
{
  ydy->dy[ydy->perm_inv[ydy->start_index + i]] = dy_i;
}

NCM_INLINE gdouble 
nc_hipert_bg_var_ydy_get_dy_i (NcHIPertBGVarYDY *ydy, guint i)
{
  return ydy->dy[ydy->perm_inv[ydy->start_index + i]];
}

NCM_INLINE void 
nc_hipert_bg_var_set_dist (NcHIPertBGVar *bg_var, NcDistance *dist)
{
  nc_distance_clear (&bg_var->dist);

  if (dist != NULL)
  {
    bg_var->dist = nc_distance_ref (dist);
    nc_distance_require_zf (dist, nc_hipert_bg_var_get_zf (bg_var));
  }
}

NCM_INLINE void 
nc_hipert_bg_var_set_recomb (NcHIPertBGVar *bg_var, NcRecomb *recomb)
{
  nc_recomb_clear (&bg_var->recomb);

  if (recomb != NULL)
  {
    bg_var->recomb = nc_recomb_ref (recomb);
    nc_recomb_require_zi (recomb, nc_hipert_bg_var_get_zf (bg_var));
  }
}

NCM_INLINE void 
nc_hipert_bg_var_set_scalefactor (NcHIPertBGVar *bg_var, NcScalefactor *a)
{
  nc_scalefactor_clear (&bg_var->a);
  
  if (a != NULL)
  {
    bg_var->a = nc_scalefactor_ref (a);
    nc_scalefactor_require_zf (a, nc_hipert_bg_var_get_zf (bg_var));
  }
}

NCM_INLINE NcDistance *
nc_hipert_bg_var_get_dist (NcHIPertBGVar *bg_var)
{
  return (bg_var->dist != NULL) ? nc_distance_ref (bg_var->dist) : bg_var->dist;
}

NCM_INLINE NcRecomb *
nc_hipert_bg_var_get_recomb (NcHIPertBGVar *bg_var)
{
  return (bg_var->recomb != NULL) ? nc_recomb_ref (bg_var->recomb) : bg_var->recomb;  
}

NCM_INLINE NcScalefactor *
nc_hipert_bg_var_get_scalefactor (NcHIPertBGVar *bg_var)
{
  return (bg_var->a != NULL) ? nc_scalefactor_ref (bg_var->a) : bg_var->a;
}

NCM_INLINE 
NcDistance *nc_hipert_bg_var_peek_dist (NcHIPertBGVar *bg_var)
{
  return bg_var->dist;
}

NCM_INLINE NcRecomb *
nc_hipert_bg_var_peek_recomb (NcHIPertBGVar *bg_var)
{
  return bg_var->recomb;
}

NCM_INLINE NcScalefactor *
nc_hipert_bg_var_peek_scalefactor (NcHIPertBGVar *bg_var)
{
  return bg_var->a;
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPERT_BG_VAR_INLINE_H_ */
