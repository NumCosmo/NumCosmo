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
 * NcmSphereNN:
 *
 * An re-implementation of Healpix.
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

  object_class->dispose  = &_ncm_sphere_nn_dispose;
  object_class->finalize = &_ncm_sphere_nn_finalize;
}

/**
 * ncm_sphere_nn_new:
 *
 * Creates a new #NcmSphereNN for a given @nobjs.
 *
 * Returns: (transfer full): a new #NcmSphereNN.
 */
NcmSphereNN *
ncm_sphere_nn_new (void)
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
 * @r: the object radius
 * @theta: the object theta
 * @phi: the object phi
 *
 * Set the object @i of @snn to the given @theta and @phi.
 *
 */
void
ncm_sphere_nn_insert (NcmSphereNN *snn, const gdouble r, const gdouble theta, const gdouble phi)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  gdouble coord[3]                = { r * sin (theta) * cos (phi), r * sin (theta) * sin (phi), r * cos (theta) };

  kdtree_insert (self->tree, coord);
}

/**
 * ncm_sphere_nn_insert_array:
 * @snn: a #NcmSphereNN
 * @r: (element-type gdouble): the object radius
 * @theta: (element-type gdouble): the object theta
 * @phi: (element-type gdouble): the object phi
 *
 * Inserts an array of objects in @snn.
 *
 */
void
ncm_sphere_nn_insert_array (NcmSphereNN *snn, GArray *r, GArray *theta, GArray *phi)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  gdouble coord[3];
  guint i;

  g_return_if_fail (theta->len == phi->len);
  g_return_if_fail (theta->len == r->len);
  g_assert_cmpuint (g_array_get_element_size (r), ==, sizeof (gdouble));
  g_assert_cmpuint (g_array_get_element_size (theta), ==, sizeof (gdouble));
  g_assert_cmpuint (g_array_get_element_size (phi), ==, sizeof (gdouble));

  for (i = 0; i < theta->len; i++)
  {
    const gdouble r_i = g_array_index (r, gdouble, i);
    gdouble sin_theta, cos_theta, sin_phi, cos_phi;

    sincos (g_array_index (theta, gdouble, i), &sin_theta, &cos_theta);
    sincos (g_array_index (phi, gdouble, i), &sin_phi, &cos_phi);

    coord[0] = r_i * sin_theta * cos_phi;
    coord[1] = r_i * sin_theta * sin_phi;
    coord[2] = r_i * cos_theta;

    kdtree_insert (self->tree, coord);
  }
}

/**
 * ncm_sphere_nn_get:
 * @snn: a #NcmSphereNN
 * @i: the object index
 * @r: (out): the object radius
 * @theta: (out): the object theta
 * @phi: (out): the object phi
 *
 * Get the object @i of @snn.
 *
 */
void
ncm_sphere_nn_get (NcmSphereNN *snn, const gint64 i, gdouble *r, gdouble *theta, gdouble *phi)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  gdouble *coord                  = self->tree->coord_table[i];

  *r     = sqrt (coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]);
  *theta = acos (coord[2] / *r);
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
 * @r: the target radius
 * @theta: the target theta
 * @phi: the target phi
 * @k: the number of nearest neighbors
 *
 * Returns: (transfer full) (element-type glong): a #GArray with the @k nearest neighbors.
 */
GArray *
ncm_sphere_nn_knn_search (NcmSphereNN *snn, const gdouble r, const gdouble theta, const gdouble phi, const gint64 k)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  gdouble coord[3]                = { r * sin (theta) * cos (phi), r * sin (theta) * sin (phi), r * cos (theta) };
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
 * @r: the target radius
 * @theta: the target theta
 * @phi: the target phi
 * @k: the number of nearest neighbors
 * @distances: (out) (transfer full) (element-type gdouble): the distances to the @k nearest neighbors
 * @indices: (out) (transfer full) (element-type glong): the indices of the @k nearest neighbors
 *
 * Computes the @k nearest neighbors of the target point (@theta, @phi) and stores the
 * distances and indices in @distances and @indices, respectively. The distances are
 * Euclidean distances in the 3D space squared. The output distances are sorted in
 * ascending order and are squared. The output indices are sorted in the same order
 * as the distances.
 *
 */
void
ncm_sphere_nn_knn_search_distances (NcmSphereNN *snn, const gdouble r, const gdouble theta, const gdouble phi, const gint64 k, GArray **distances, GArray **indices)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  gdouble coord[3]                = { r * sin (theta) * cos (phi), r * sin (theta) * sin (phi), r * cos (theta) };
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

/**
 * ncm_sphere_nn_knn_search_distances_batch:
 * @snn: a #NcmSphereNN
 * @r: (element-type gdouble): the target radius
 * @theta: (element-type gdouble): the target theta
 * @phi: (element-type gdouble): the target phi
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
ncm_sphere_nn_knn_search_distances_batch (NcmSphereNN *snn, GArray *r, GArray *theta, GArray *phi, const gint64 k, GArray **distances, GArray **indices)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);

  g_assert_cmpuint (theta->len, ==, phi->len);
  g_assert_cmpuint (theta->len, ==, r->len);
  g_assert_cmpuint (g_array_get_element_size (r), ==, sizeof (gdouble));
  g_assert_cmpuint (g_array_get_element_size (theta), ==, sizeof (gdouble));
  g_assert_cmpuint (g_array_get_element_size (phi), ==, sizeof (gdouble));

  *distances = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), k * theta->len);
  *indices   = g_array_sized_new (FALSE, FALSE, sizeof (glong), k * theta->len);

  {
    gdouble coord[3];
    guint i;

    for (i = 0; i < theta->len; i++)
    {
      const gdouble r_i = g_array_index (r, gdouble, i);
      gdouble sin_theta, cos_theta, sin_phi, cos_phi;
      guint n = 0;

      sincos (g_array_index (theta, gdouble, i), &sin_theta, &cos_theta);
      sincos (g_array_index (phi, gdouble, i), &sin_phi, &cos_phi);

      coord[0] = r_i * sin_theta * cos_phi;
      coord[1] = r_i * sin_theta * sin_phi;
      coord[2] = r_i * cos_theta;

      rb_knn_list_table_t *table = kdtree_knn_search (self->tree, coord, k);
      rb_knn_list_traverser_t trav;
      knn_list_t *p;

      p = rb_knn_list_t_first (&trav, table);

      do {
        g_array_append_val (*distances, p->distance);
        g_array_append_val (*indices, p->node->coord_index);
        n++;
      } while ((p = rb_knn_list_t_next (&trav)) != NULL);

      g_assert_cmpuint (n, ==, k);

      rb_knn_list_destroy (table);
    }
  }
}

/**
 * ncm_sphere_nn_dump_tree:
 * @snn: a #NcmSphereNN
 *
 * Print to the standard output the tree structure of @snn.
 *
 */
void
ncm_sphere_nn_dump_tree (NcmSphereNN *snn)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);

  kdtree_dump (self->tree);
  fflush (stdout);
}

