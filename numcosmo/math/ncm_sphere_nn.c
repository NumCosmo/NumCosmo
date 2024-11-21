/***************************************************************************
 *            ncm_sphere_nn.c
 *
 *  Wed Nov 20 19:23:40 2024
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_sphere_nn.c
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

/**
 * SECTION:ncm_sphere_nn
 * @title: NcmSphereNN
 * @short_description: An re-implementation of Healpix.
 *
 * NN pixalization/manipulation algorithms, Ylm decomposition.
 *
 */
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sphere_nn.h"
#include "math/ncm_vector.h"
#include "math/ncm_matrix.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_trig.h>
#include "misc/kdtree.h"
#include "misc/rb_knn_list.h"
#endif /* NUMCOSMO_GIR_SCAN */


typedef struct _NcmSphereNNPrivate
{
  struct kdtree *tree;
} NcmSphereNNPrivate;

enum
{
  PROP_0,
  PROP_NOBJS,
};

struct _NcmSphereNN
{
  GObject parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmSphereNN, ncm_sphere_nn, G_TYPE_OBJECT)

static void
ncm_sphere_nn_init (NcmSphereNN *snn)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);

  self->tree = kdtree_init (3);
}

static void
_ncm_sphere_nn_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NCM_IS_SPHERE_NN (object));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sphere_nn_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSphereNN *snn                = NCM_SPHERE_NN (object);
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);

  g_return_if_fail (NCM_IS_SPHERE_NN (object));

  switch (prop_id)
  {
    case PROP_NOBJS:
      g_value_set_int64 (value, self->tree->count);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_sphere_nn_dispose (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sphere_nn_parent_class)->dispose (object);
}

static void
_ncm_sphere_nn_finalize (GObject *object)
{
  NcmSphereNN *snn                = NCM_SPHERE_NN (object);
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);

  g_clear_pointer (&self->tree, kdtree_destroy);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_sphere_nn_parent_class)->finalize (object);
}

static void
ncm_sphere_nn_class_init (NcmSphereNNClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_sphere_nn_set_property;
  object_class->get_property = &_ncm_sphere_nn_get_property;
  object_class->dispose      = &_ncm_sphere_nn_dispose;
  object_class->finalize     = &_ncm_sphere_nn_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NOBJS,
                                   g_param_spec_int64 ("nobjs",
                                                       NULL,
                                                       "number of objects",
                                                       0, G_MAXINT64, 0,
                                                       G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_sphere_nn_new:
 *
 * Creates a new #NcmSphereNN for a given @nobjs.
 *
 * Returns: (transfer full): a new #NcmSphereNN.
 */
NcmSphereNN *
ncm_sphere_nn_new ()
{
  NcmSphereNN *snn = g_object_new (NCM_TYPE_SPHERE_NN,
                                   NULL);

  return snn;
}

/**
 * ncm_sphere_nn_ref:
 * @snn: a #NcmSphereNN
 *
 * Increases the reference count of @snn.
 *
 * Returns: (transfer full): @snn.
 */
NcmSphereNN *
ncm_sphere_nn_ref (NcmSphereNN *snn)
{
  return g_object_ref (snn);
}

/**
 * ncm_sphere_nn_free:
 * @snn: a #NcmSphereNN
 *
 * Decreases the reference count of @snn. When its reference count
 * drops to 0, the object is finalized (i.e. its memory is freed).
 *
 */
void
ncm_sphere_nn_free (NcmSphereNN *snn)
{
  g_object_unref (snn);
}

/**
 * ncm_sphere_nn_clear:
 * @snn: a #NcmSphereNN
 *
 * If *@snn is not %NULL, decreases the reference count of @snn.
 * When its reference count drops to 0, the object is finalized
 * (i.e. its memory is freed).
 * Set *@snn to %NULL.
 *
 */
void
ncm_sphere_nn_clear (NcmSphereNN **snn)
{
  g_clear_object (snn);
}

/**
 * ncm_sphere_nn_insert:
 * @snn: a #NcmSphereNN
 * @theta: the object theta
 * @phi: the object phi
 *
 * Set the object @i of @snn to the given @theta and @phi.
 *
 */
void
ncm_sphere_nn_insert (NcmSphereNN *snn, const gdouble theta, const gdouble phi)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  gdouble coord[3]                = { sin (theta) * cos (phi), sin (theta) * sin (phi), cos (theta) };

  kdtree_insert (self->tree, coord);
}

/**
 * ncm_sphere_nn_get:
 * @snn: a #NcmSphereNN
 * @i: the object index
 * @theta: (out): the object theta
 * @phi: (out): the object phi
 *
 * Get the object @i of @snn.
 *
 */
void
ncm_sphere_nn_get (NcmSphereNN *snn, const gint64 i, gdouble *theta, gdouble *phi)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  gdouble *coord                  = self->tree->coord_table[i];

  *theta = acos (coord[2]);
  *phi   = atan2 (coord[1], coord[0]);
}

/**
 * ncm_sphere_nn_get_n:
 * @snn: a #NcmSphereNN
 *
 * Returns: the number of objects of @snn.
 */
gint64
ncm_sphere_nn_get_n (NcmSphereNN *snn)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);

  return self->tree->count;
}

/**
 * ncm_sphere_nn_rebuild:
 * @snn: a #NcmSphereNN
 *
 * Rebuild the @snn or build it if it is empty.
 *
 */
void
ncm_sphere_nn_rebuild (NcmSphereNN *snn)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);

  kdtree_rebuild (self->tree);
}

/**
 * ncm_sphere_nn_knn_search:
 * @snn: a #NcmSphereNN
 * @theta: the target theta
 * @phi: the target phi
 * @k: the number of nearest neighbors
 *
 * Returns: (transfer full) (element-type glong): a #GArray with the @k nearest neighbors.
 */
GArray *
ncm_sphere_nn_knn_search (NcmSphereNN *snn, const gdouble theta, const gdouble phi, const gint64 k)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  gdouble coord[3]                = { sin (theta) * cos (phi), sin (theta) * sin (phi), cos (theta) };
  rb_knn_list_table_t *table      = kdtree_knn_search (self->tree, coord, k);
  GArray *indices                 = g_array_sized_new (FALSE, FALSE, sizeof (glong), k);
  rb_knn_list_traverser_t trav;
  knn_list_t *p;

  p = rb_knn_list_t_first (&trav, table);

  do {
    g_array_append_val (indices, p->node->coord_index);
  } while ((p = rb_knn_list_t_next (&trav)) != NULL);


  rb_knn_list_destroy (table);

  return indices;
}

/**
 * ncm_sphere_nn_knn_search_distances:
 * @snn: a #NcmSphereNN
 * @theta: the target theta
 * @phi: the target phi
 * @k: the number of nearest neighbors
 * @distances: (out) (transfer full) (element-type gdouble): the distances to the @k nearest neighbors
 * @indices: (out) (transfer full) (element-type glong): the indices of the @k nearest neighbors
 *
 * Computes the @k nearest neighbors of the target point (@theta, @phi) and stores the
 * distances and indices in @distances and @indices, respectively. The distances are
 * Euclidean distances in the 3D space squared.
 *
 */
void
ncm_sphere_nn_knn_search_distances (NcmSphereNN *snn, const gdouble theta, const gdouble phi, const gint64 k, GArray **distances, GArray **indices)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  gdouble coord[3]                = { sin (theta) * cos (phi), sin (theta) * sin (phi), cos (theta) };
  rb_knn_list_table_t *table      = kdtree_knn_search (self->tree, coord, k);
  rb_knn_list_traverser_t trav;
  knn_list_t *p;

  *distances = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), k);
  *indices   = g_array_sized_new (FALSE, FALSE, sizeof (glong), k);

  p = rb_knn_list_t_first (&trav, table);

  do {
    g_array_append_val (*distances, p->distance);
    g_array_append_val (*indices, p->node->coord_index);
  } while ((p = rb_knn_list_t_next (&trav)) != NULL);

  rb_knn_list_destroy (table);
}

