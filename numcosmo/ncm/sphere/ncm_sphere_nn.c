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
 * Nearest-neighbour search for points given in spherical coordinates.
 *
 * Points are inserted as $(r, \theta, \phi)$, with $\theta$ the polar angle from the
 * $z$ axis and $\phi$ the azimuth, and stored as Cartesian points
 * $r(\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)$ in a k-d tree. A point's index
 * is its insertion order, starting at zero.
 *
 * The searches return the $k$ nearest points in 3D Euclidean distance, sorted by
 * increasing distance, and report the distance squared. For points on a sphere of
 * radius $r$ this is the squared chord, $2r^2(1 - \cos\gamma)$ for an angular separation
 * $\gamma$, so the order is that of the angular separation.
 *
 * The tree is built by ncm_sphere_nn_rebuild(): call it after inserting, since a search
 * aborts when points were inserted after the last rebuild. A search needs
 * $1 \le k \le$ ncm_sphere_nn_get_n().
 *
 */
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/sphere/ncm_sphere_nn.h"

#ifndef NUMCOSMO_GIR_SCAN
#include "external/misc/kdtree.h"
#include "external/misc/rb_knn_list.h"
#include <stdio.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmSphereNNPrivate
{
  struct kdtree *tree;
  gint64 n_built;
} NcmSphereNNPrivate;

struct _NcmSphereNN
{
  GObject parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmSphereNN, ncm_sphere_nn, G_TYPE_OBJECT)

static void
ncm_sphere_nn_init (NcmSphereNN *snn)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);

  self->tree    = kdtree_init (3);
  self->n_built = -1;
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

  object_class->finalize = &_ncm_sphere_nn_finalize;
}

static void _ncm_sphere_nn_to_cartesian (const gdouble r, const gdouble theta, const gdouble phi, gdouble coord[3]);
static void _ncm_sphere_nn_check_search (NcmSphereNNPrivate * const self, const gint64 k, const gchar *func);
static void _ncm_sphere_nn_search (NcmSphereNNPrivate * const self, gdouble coord[3], const gint64 k, GArray *distances, GArray *indices);

/**
 * ncm_sphere_nn_new:
 *
 * Creates a new, empty #NcmSphereNN.
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
 * @r: the point radius
 * @theta: the point polar angle $\theta$ (radians)
 * @phi: the point azimuth $\phi$ (radians)
 *
 * Appends the point $(r, \theta, \phi)$ to @snn; its index is the number of points
 * inserted before it. It enters the searches after the next ncm_sphere_nn_rebuild().
 *
 */
void
ncm_sphere_nn_insert (NcmSphereNN *snn, const gdouble r, const gdouble theta, const gdouble phi)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  gdouble coord[3];

  _ncm_sphere_nn_to_cartesian (r, theta, phi, coord);
  kdtree_insert (self->tree, coord);
}

/**
 * ncm_sphere_nn_insert_array:
 * @snn: a #NcmSphereNN
 * @r: (element-type gdouble): the point radii
 * @theta: (element-type gdouble): the point polar angles (radians)
 * @phi: (element-type gdouble): the point azimuths (radians)
 *
 * Appends the points $(r_i, \theta_i, \phi_i)$ to @snn in order, see
 * ncm_sphere_nn_insert(). The three arrays must have the same length.
 *
 */
void
ncm_sphere_nn_insert_array (NcmSphereNN *snn, GArray *r, GArray *theta, GArray *phi)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  guint i;

  if ((theta->len != phi->len) || (theta->len != r->len))
    g_error ("ncm_sphere_nn_insert_array: the arrays have different lengths (%u, %u, %u).", r->len, theta->len, phi->len);

  g_assert_cmpuint (g_array_get_element_size (r), ==, sizeof (gdouble));
  g_assert_cmpuint (g_array_get_element_size (theta), ==, sizeof (gdouble));
  g_assert_cmpuint (g_array_get_element_size (phi), ==, sizeof (gdouble));

  for (i = 0; i < theta->len; i++)
  {
    gdouble coord[3];

    _ncm_sphere_nn_to_cartesian (g_array_index (r, gdouble, i), g_array_index (theta, gdouble, i), g_array_index (phi, gdouble, i), coord);
    kdtree_insert (self->tree, coord);
  }
}

/**
 * ncm_sphere_nn_get:
 * @snn: a #NcmSphereNN
 * @i: the point index
 * @r: (out): the point radius
 * @theta: (out): the point polar angle, in $[0, \pi]$
 * @phi: (out): the point azimuth, in $(-\pi, \pi]$
 *
 * Gets the point @i of @snn, recovered from its Cartesian coordinates: the angles are
 * the inserted ones up to rounding and to the range above.
 *
 */
void
ncm_sphere_nn_get (NcmSphereNN *snn, const gint64 i, gdouble *r, gdouble *theta, gdouble *phi)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  gdouble *coord;

  if ((i < 0) || (i >= (gint64) self->tree->count))
    g_error ("ncm_sphere_nn_get: index %" G_GINT64_FORMAT " out of range, the tree holds %zu points.", i, self->tree->count);

  coord  = self->tree->coord_table[i];
  *r     = sqrt (coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]);
  *theta = acos (coord[2] / *r);
  *phi   = atan2 (coord[1], coord[0]);
}

/**
 * ncm_sphere_nn_get_n:
 * @snn: a #NcmSphereNN
 *
 * Returns: the number of points inserted in @snn.
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
 * Builds the search tree over all the points inserted so far, replacing the previous
 * one.
 *
 */
void
ncm_sphere_nn_rebuild (NcmSphereNN *snn)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);

  kdtree_rebuild (self->tree);
  self->n_built = self->tree->count;
}

/**
 * ncm_sphere_nn_knn_search:
 * @snn: a #NcmSphereNN
 * @r: the target radius
 * @theta: the target polar angle (radians)
 * @phi: the target azimuth (radians)
 * @k: the number of nearest neighbours
 *
 * Finds the @k points of @snn nearest to $(r, \theta, \phi)$.
 *
 * Returns: (transfer full) (element-type glong): the indices of the @k nearest points,
 * nearest first.
 */
GArray *
ncm_sphere_nn_knn_search (NcmSphereNN *snn, const gdouble r, const gdouble theta, const gdouble phi, const gint64 k)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  GArray *distances, *indices;
  gdouble coord[3];

  _ncm_sphere_nn_check_search (self, k, "ncm_sphere_nn_knn_search");
  _ncm_sphere_nn_to_cartesian (r, theta, phi, coord);

  distances = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), k);
  indices   = g_array_sized_new (FALSE, FALSE, sizeof (glong), k);

  _ncm_sphere_nn_search (self, coord, k, distances, indices);
  g_array_unref (distances);

  return indices;
}

/**
 * ncm_sphere_nn_knn_search_distances:
 * @snn: a #NcmSphereNN
 * @r: the target radius
 * @theta: the target polar angle (radians)
 * @phi: the target azimuth (radians)
 * @k: the number of nearest neighbours
 * @distances: (out) (transfer full) (element-type gdouble): the squared distances to the @k nearest points
 * @indices: (out) (transfer full) (element-type glong): the indices of the @k nearest points
 *
 * Finds the @k points of @snn nearest to $(r, \theta, \phi)$, nearest first, with their
 * squared Euclidean distances to it.
 *
 */
void
ncm_sphere_nn_knn_search_distances (NcmSphereNN *snn, const gdouble r, const gdouble theta, const gdouble phi, const gint64 k, GArray **distances, GArray **indices)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  gdouble coord[3];

  _ncm_sphere_nn_check_search (self, k, "ncm_sphere_nn_knn_search_distances");
  _ncm_sphere_nn_to_cartesian (r, theta, phi, coord);

  *distances = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), k);
  *indices   = g_array_sized_new (FALSE, FALSE, sizeof (glong), k);

  _ncm_sphere_nn_search (self, coord, k, *distances, *indices);
}

/**
 * ncm_sphere_nn_knn_search_distances_batch:
 * @snn: a #NcmSphereNN
 * @r: (element-type gdouble): the target radii
 * @theta: (element-type gdouble): the target polar angles (radians)
 * @phi: (element-type gdouble): the target azimuths (radians)
 * @k: the number of nearest neighbours
 * @distances: (out) (transfer full) (element-type gdouble): the squared distances, @k per target
 * @indices: (out) (transfer full) (element-type glong): the indices, @k per target
 *
 * Runs ncm_sphere_nn_knn_search_distances() for each target and concatenates the
 * results: entries $[jk, (j+1)k)$ belong to target $j$.
 *
 */
void
ncm_sphere_nn_knn_search_distances_batch (NcmSphereNN *snn, GArray *r, GArray *theta, GArray *phi, const gint64 k, GArray **distances, GArray **indices)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);
  guint i;

  if ((theta->len != phi->len) || (theta->len != r->len))
    g_error ("ncm_sphere_nn_knn_search_distances_batch: the arrays have different lengths (%u, %u, %u).", r->len, theta->len, phi->len);

  g_assert_cmpuint (g_array_get_element_size (r), ==, sizeof (gdouble));
  g_assert_cmpuint (g_array_get_element_size (theta), ==, sizeof (gdouble));
  g_assert_cmpuint (g_array_get_element_size (phi), ==, sizeof (gdouble));

  _ncm_sphere_nn_check_search (self, k, "ncm_sphere_nn_knn_search_distances_batch");

  *distances = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), k * theta->len);
  *indices   = g_array_sized_new (FALSE, FALSE, sizeof (glong), k * theta->len);

  for (i = 0; i < theta->len; i++)
  {
    gdouble coord[3];

    _ncm_sphere_nn_to_cartesian (g_array_index (r, gdouble, i), g_array_index (theta, gdouble, i), g_array_index (phi, gdouble, i), coord);
    _ncm_sphere_nn_search (self, coord, k, *distances, *indices);
  }
}

/**
 * ncm_sphere_nn_dump_tree:
 * @snn: a #NcmSphereNN
 *
 * Prints the tree structure of @snn to the standard output.
 *
 */
void
ncm_sphere_nn_dump_tree (NcmSphereNN *snn)
{
  NcmSphereNNPrivate * const self = ncm_sphere_nn_get_instance_private (snn);

  kdtree_dump (self->tree);
  fflush (stdout);
}

static void
_ncm_sphere_nn_to_cartesian (const gdouble r, const gdouble theta, const gdouble phi, gdouble coord[3])
{
  gdouble sin_theta, cos_theta, sin_phi, cos_phi;

  sincos (theta, &sin_theta, &cos_theta);
  sincos (phi, &sin_phi, &cos_phi);

  coord[0] = r * sin_theta * cos_phi;
  coord[1] = r * sin_theta * sin_phi;
  coord[2] = r * cos_theta;
}

/* A search on a tree that misses points, or asks for more points than it holds, used to
 * return a wrong answer or dereference an empty result. */
static void
_ncm_sphere_nn_check_search (NcmSphereNNPrivate * const self, const gint64 k, const gchar *func)
{
  if (self->n_built != (gint64) self->tree->count)
    g_error ("%s: %zu points inserted but the tree was built with %" G_GINT64_FORMAT "; "
             "call ncm_sphere_nn_rebuild() after inserting.", func, self->tree->count, self->n_built);

  if ((k < 1) || (k > (gint64) self->tree->count) || (k > G_MAXINT))
    g_error ("%s: k = %" G_GINT64_FORMAT " is out of range, the tree holds %zu points.", func, k, self->tree->count);
}

/* Appends the k nearest points to coord, nearest first. */
static void
_ncm_sphere_nn_search (NcmSphereNNPrivate * const self, gdouble coord[3], const gint64 k, GArray *distances, GArray *indices)
{
  rb_knn_list_table_t *table = kdtree_knn_search (self->tree, coord, (gint) k);
  rb_knn_list_traverser_t trav;
  knn_list_t *p;
  gint64 n = 0;

  for (p = rb_knn_list_t_first (&trav, table); p != NULL; p = rb_knn_list_t_next (&trav))
  {
    g_array_append_val (distances, p->distance);
    g_array_append_val (indices, p->node->coord_index);
    n++;
  }

  rb_knn_list_destroy (table);

  g_assert_cmpint (n, ==, k);
}

