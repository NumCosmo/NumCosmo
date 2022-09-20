/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_iset.h
 *
 *  Sun April 7 16:59:36 2021
 *  Copyright  2021  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_iset.h
 * Copyright (C) 2021 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_ISET_H_
#define _NCM_ISET_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>

G_BEGIN_DECLS

#define NCM_TYPE_ISET             (ncm_iset_get_type ())
#define NCM_ISET(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_ISET, NcmISet))
#define NCM_ISET_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_ISET, NcmISetClass))
#define NCM_IS_ISET(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_ISET))
#define NCM_IS_ISET_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_ISET))
#define NCM_ISET_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_ISET, NcmISetClass))

typedef struct _NcmISetClass NcmISetClass;
typedef struct _NcmISet NcmISet;
typedef struct _NcmISetPrivate NcmISetPrivate;

struct _NcmISetClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmISet
{
  /*< private >*/
  GObject parent_instance;
  NcmISetPrivate *priv;
};

GType ncm_iset_get_type (void) G_GNUC_CONST;

NcmISet *ncm_iset_new (guint n);
NcmISet *ncm_iset_ref (NcmISet *iset);

void ncm_iset_free (NcmISet *iset);
void ncm_iset_clear (NcmISet **iset);

void ncm_iset_add_range (NcmISet *iset, gint ii, gint fi);
void ncm_iset_add (NcmISet *iset, gint i);
void ncm_iset_del (NcmISet *iset, gint i);
void ncm_iset_reset (NcmISet *iset);

void ncm_iset_copy (NcmISet *iset, NcmISet *target);

guint ncm_iset_get_max_size (NcmISet *iset);
guint ncm_iset_get_len (NcmISet *iset);
gdouble ncm_iset_get_vector_max (NcmISet *iset, NcmVector *v, gint *max_i);
NcmVector *ncm_iset_get_subvector (NcmISet *iset, NcmVector *v, NcmVector *v_dup);
GArray *ncm_iset_get_subarray (NcmISet *iset, GArray *a, GArray *a_dup);
NcmMatrix *ncm_iset_get_submatrix (NcmISet *iset, NcmMatrix *M, NcmMatrix *M_dup);
NcmMatrix *ncm_iset_get_submatrix_cols (NcmISet *iset, NcmMatrix *M, NcmMatrix *M_dup);
NcmMatrix *ncm_iset_get_sym_submatrix (NcmISet *iset, gchar UL, NcmMatrix *M, NcmMatrix *M_dup);
void ncm_iset_get_subset_vec_lt (NcmISet *iset, NcmISet *out, NcmVector *v, const gdouble tol);
void ncm_iset_remove_subset (NcmISet *iset, NcmISet *target);
guint ncm_iset_remove_smallest_subset (NcmISet *iset, NcmISet *target, NcmVector *v, guint max_remove);
guint ncm_iset_add_largest_subset (NcmISet *iset, NcmVector *v, const gdouble min, const gdouble add_frac);
void ncm_iset_set_complement (NcmISet *iset, NcmISet *cmplm);

NcmVector *ncm_iset_get_vector_inv_cmp (NcmISet *iset, NcmVector *u, NcmVector *v, NcmVector *v_dup);

void ncm_iset_set_subvector (NcmISet *iset, NcmVector *v, NcmVector *sub);

void ncm_iset_log_vals (NcmISet *iset, const gchar *prefix);

G_END_DECLS

#endif /* _NCM_ISET_H_ */
