/*
 * Copyright (C) 2017, Leo Ma <begeekmyfriend@gmail.com>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <glib.h>

#include "kdtree.h"
#include "rb_knn_list.h"

static inline int
is_leaf (struct kdnode *node)
{
  return node->left == node->right;
}

static inline void
swap (long *a, long *b)
{
  long tmp = *a;

  *a = *b;
  *b = tmp;
}

static inline double
square (double d)
{
  return d * d;
}

static inline double
distance (double *c1, double *c2, int dim)
{
  double distance = 0;

  while (dim-- > 0)
  {
    distance += square (*c1++ - *c2++);
  }

  return distance;
}

static inline double
D (struct kdtree *tree, long index, int r)
{
  return tree->coord_table[index][r];
}

static inline int
knn_search_on (struct kdtree *tree, rb_knn_list_table_t *table, double *knn_distance, int k, double value, double target)
{
  return table->rb_knn_list_count < k || square (target - value) < *knn_distance;
}

static inline void
coord_index_reset (struct kdtree *tree)
{
  long i;

  for (i = 0; i < tree->capacity; i++)
  {
    tree->coord_indexes[i] = i;
  }
}

static inline void
coord_table_reset (struct kdtree *tree)
{
  long i;

  for (i = 0; i < tree->capacity; i++)
  {
    tree->coord_table[i] = tree->coords + i * tree->dim;
  }
}

static void
insert_sort (struct kdtree *tree, long low, long high, int r)
{
  long i, j;
  long *indexes = tree->coord_indexes;

  for (i = low + 1; i <= high; i++)
  {
    long tmp_idx     = indexes[i];
    double tmp_value = D (tree, indexes[i], r);

    j = i - 1;

    for ( ; j >= low && D (tree, indexes[j], r) > tmp_value; j--)
    {
      indexes[j + 1] = indexes[j];
    }

    indexes[j + 1] = tmp_idx;
  }
}

static void
quicksort (struct kdtree *tree, long low, long high, int r)
{
  if (high - low <= 32)
  {
    insert_sort (tree, low, high, r);

    /*bubble_sort(tree, low, high, r); */
    return;
  }

  long *indexes = tree->coord_indexes;
  /* median of 3 */
  long mid = low + (high - low) / 2;

  if (D (tree, indexes[low], r) > D (tree, indexes[mid], r))
    swap (indexes + low, indexes + mid);

  if (D (tree, indexes[low], r) > D (tree, indexes[high], r))
    swap (indexes + low, indexes + high);

  if (D (tree, indexes[high], r) > D (tree, indexes[mid], r))
    swap (indexes + high, indexes + mid);


  /* D(indexes[low]) <= D(indexes[high]) <= D(indexes[mid]) */
  double pivot = D (tree, indexes[high], r);

  /* 3-way partition
   * +---------+-----------+---------+-------------+---------+
   * |  pivot  |  <=pivot  |   ?     |  >=pivot    |  pivot  |
   * +---------+-----------+---------+-------------+---------+
   * low     lt             i       j               gt    high
   */
  long i  = low - 1;
  long lt = i;
  long j  = high;
  long gt = j;

  for ( ; ;)
  {
    while (D (tree, indexes[++i], r) < pivot)
    {
    }

    while (D (tree, indexes[--j], r) > pivot && j > low)
    {
    }

    if (i >= j)
      break;

    swap (indexes + i, indexes + j);

    if (D (tree, indexes[i], r) == pivot)
      swap (&indexes[++lt], &indexes[i]);

    if (D (tree, indexes[j], r) == pivot)
      swap (&indexes[--gt], &indexes[j]);
  }

  /* i == j or j + 1 == i */
  swap (indexes + i, indexes + high);

  /* Move equal elements to the middle of array */
  long x, y;

  for (x = low, j = i - 1; x <= lt && j > lt; x++, j--)
    swap (indexes + x, indexes + j);

  for (y = high, i = i + 1; y >= gt && i < gt; y--, i++)
    swap (indexes + y, indexes + i);

  quicksort (tree, low, j - lt + x - 1, r);
  quicksort (tree, i + y - gt, high, r);
}

static struct kdnode *

kdnode_alloc (double *coord, long index, int r)
{
  struct kdnode *node = malloc (sizeof (*node));

  if (node != NULL)
  {
    memset (node, 0, sizeof (*node));
    node->coord       = coord;
    node->coord_index = index;
    node->r           = r;
  }

  return node;
}

static void
kdnode_free (struct kdnode *node)
{
  free (node);
}

static int
coord_cmp (double *c1, double *c2, int dim)
{
  int i;
  double ret = 0.0;

  for (i = 0; i < dim; i++)
  {
    ret = *c1++ - *c2++;

    if (fabs (ret) >= DBL_EPSILON)
      return ret > 0 ? 1 : -1;
  }

  if (fabs (ret) < DBL_EPSILON)
    return 0;
  else
    return ret > 0 ? 1 : -1;
}

static void
resize (struct kdtree *tree)
{
  tree->capacity     *= 2;
  tree->coords        = realloc (tree->coords, tree->dim * sizeof (double) * tree->capacity);
  tree->coord_table   = realloc (tree->coord_table, sizeof (double *) * tree->capacity);
  tree->coord_indexes = realloc (tree->coord_indexes, sizeof (long) * tree->capacity);
  coord_table_reset (tree);
  coord_index_reset (tree);
}

static void
kdnode_dump (struct kdnode *node, int dim)
{
  int i;

  if (node->coord != NULL)
  {
    printf ("(");

    for (i = 0; i < dim; i++)
    {
      if (i != dim - 1)
        printf ("%.2f,", node->coord[i]);
      else
        printf ("%.2f)\n", node->coord[i]);
    }
  }
  else
  {
    printf ("(none)\n");
  }
}

void
kdtree_insert (struct kdtree *tree, double *coord)
{
  if (tree->count + 1 > tree->capacity)
    resize (tree);

  memcpy (tree->coord_table[tree->count++], coord, tree->dim * sizeof (double));
}

static void
knn_pickup (struct kdtree *tree, struct kdnode *node, rb_knn_list_table_t *table, double *knn_distance, double *target, int k)
{
  double dist = distance (node->coord, target, tree->dim);

  if (table->rb_knn_list_count < k)
  {
    knn_list_t *log = g_slice_new (knn_list_t);

    log->node     = node;
    log->distance = dist;
    rb_knn_list_insert (table, log);

    *knn_distance = MAX (*knn_distance, dist);
  }
  else if (dist < *knn_distance)
  {
    knn_list_t *log = g_slice_new (knn_list_t);
    rb_knn_list_traverser_t trav;
    knn_list_t *last, *prev;

    log->node     = node;
    log->distance = dist;

    rb_knn_list_insert (table, log);
    last = rb_knn_list_t_last (&trav, table);

    if ((prev = rb_knn_list_t_prev (&trav)) != NULL)
      *knn_distance = prev->distance;

    rb_knn_list_delete (table, last);

    g_slice_free (knn_list_t, last);
  }
}

static void
kdtree_search_recursive (struct kdtree *tree, struct kdnode *node, rb_knn_list_table_t *table, double *knn_distance, double *target, int k, int *pickup, double *min_shift, double min_dist2)
{
  if (node == NULL)
  {
    return;
  }
  else
  {
    int r = node->r;
    double min_shift_new[tree->dim];
    double *min_shift_r;
    double *min_shift_l;
    double min_dist2_new;
    double min_dist2_r;
    double min_dist2_l;

    memcpy (min_shift_new, min_shift, tree->dim * sizeof (double));
    min_shift_new[r] = node->coord[r] - target[r];
    min_dist2_new    = min_dist2 + square (min_shift_new[r]) - square (min_shift[r]);

    if (node->coord[r] >= target[r])
    {
      min_shift_r = min_shift_new;
      min_shift_l = min_shift;
      min_dist2_r = min_dist2_new;
      min_dist2_l = min_dist2;
    }
    else
    {
      min_shift_r = min_shift;
      min_shift_l = min_shift_new;
      min_dist2_r = min_dist2;
      min_dist2_l = min_dist2_new;
    }

    if (!knn_search_on (tree, table, knn_distance, k, node->coord[r], target[r]))
      return;

    /* Testing if the branch can be prunned. */
    if (table->rb_knn_list_count >= k)
      if (min_dist2 > *knn_distance)
        return;


    if (*pickup)
    {
      knn_pickup (tree, node, table, knn_distance, target, k);
      kdtree_search_recursive (tree, node->left, table, knn_distance, target, k, pickup, min_shift_l, min_dist2_l);
      kdtree_search_recursive (tree, node->right, table, knn_distance, target, k, pickup, min_shift_r, min_dist2_r);
    }
    else
    {
      /*printf ("Else!\n"); */

      if (is_leaf (node))
      {
        *pickup = 1;
      }
      else
      {
        if (target[r] <= node->coord[r])
        {
          kdtree_search_recursive (tree, node->left, table, knn_distance, target, k, pickup, min_shift_l, min_dist2_l);
          kdtree_search_recursive (tree, node->right, table, knn_distance, target, k, pickup, min_shift_r, min_dist2_r);
        }
        else
        {
          kdtree_search_recursive (tree, node->right, table, knn_distance, target, k, pickup, min_shift_r, min_dist2_r);
          kdtree_search_recursive (tree, node->left, table, knn_distance, target, k, pickup, min_shift_l, min_dist2_l);
        }
      }

      /* back track and pick up  */
      if (*pickup)
        knn_pickup (tree, node, table, knn_distance, target, k);
    }
  }
}

void *
kdtree_knn_search (struct kdtree *tree, double *target, int k)
{
  if (k > 0)
  {
    rb_knn_list_table_t *table = rb_knn_list_create (NULL);

    int pickup          = 0;
    double knn_distance = 0.0;
    double min_shift[tree->dim];

    memset (min_shift, 0, sizeof (double) * tree->dim);

    kdtree_search_recursive (tree, tree->root, table, &knn_distance, target, k, &pickup, min_shift, 0.0);

    return table;
  }

  return NULL;
}

void
kdtree_delete (struct kdtree *tree, double *coord)
{
  int r                 = 0;
  struct kdnode *node   = tree->root;
  struct kdnode *parent = node;

  while (node != NULL)
  {
    if (node->coord == NULL)
    {
      if (parent->right->coord == NULL)
      {
        break;
      }
      else
      {
        node = parent->right;
        continue;
      }
    }

    if (coord[r] < node->coord[r])
    {
      parent = node;
      node   = node->left;
    }
    else if (coord[r] > node->coord[r])
    {
      parent = node;
      node   = node->right;
    }
    else
    {
      int ret = coord_cmp (coord, node->coord, tree->dim);

      if (ret < 0)
      {
        parent = node;
        node   = node->left;
      }
      else if (ret > 0)
      {
        parent = node;
        node   = node->right;
      }
      else
      {
        node->coord = NULL;
        break;
      }
    }

    r = (r + 1) % tree->dim;
  }
}

static void
kdnode_build (struct kdtree *tree, struct kdnode **nptr, int r, long low, long high)
{
  if (low == high)
  {
    long index = tree->coord_indexes[low];

    *nptr = kdnode_alloc (tree->coord_table[index], index, r);
  }
  else if (low < high)
  {
    /* Sort and fetch the median to build a balanced BST */
    quicksort (tree, low, high, r);

    long median         = low + (high - low) / 2;
    long median_index   = tree->coord_indexes[median];
    struct kdnode *node = *nptr = kdnode_alloc (tree->coord_table[median_index], median_index, r);

    r = (r + 1) % tree->dim;

    kdnode_build (tree, &node->left, r, low, median - 1);
    kdnode_build (tree, &node->right, r, median + 1, high);
  }
}

static void
kdtree_build (struct kdtree *tree)
{
  kdnode_build (tree, &tree->root, 0, 0, tree->count - 1);
}

void
kdtree_rebuild (struct kdtree *tree)
{
  coord_index_reset (tree);
  kdtree_build (tree);
}

struct kdtree *

kdtree_init (int dim)
{
  struct kdtree *tree = malloc (sizeof (*tree));

  if (tree != NULL)
  {
    tree->root          = NULL;
    tree->dim           = dim;
    tree->count         = 0;
    tree->capacity      = 65536;
    tree->coords        = malloc (dim * sizeof (double) * tree->capacity);
    tree->coord_table   = malloc (sizeof (double *) * tree->capacity);
    tree->coord_indexes = malloc (sizeof (long) * tree->capacity);
    coord_index_reset (tree);
    coord_table_reset (tree);
  }

  return tree;
}

static void
kdnode_destroy (struct kdnode *node)
{
  if (node == NULL)
    return;

  kdnode_destroy (node->left);
  kdnode_destroy (node->right);
  kdnode_free (node);
}

void
kdtree_destroy (struct kdtree *tree)
{
  kdnode_destroy (tree->root);
  free (tree->coords);
  free (tree->coord_table);
  free (tree->coord_indexes);
  free (tree);
}

#define _KDTREE_DEBUG

#ifdef _KDTREE_DEBUG
struct kdnode_backlog
{
  struct kdnode *node;
  int next_sub_idx;
};

void
kdtree_dump (struct kdtree *tree)
{
  int level                    = 0;
  struct kdnode *node          = tree->root;
  struct kdnode_backlog *p_nbl = NULL;
  struct kdnode_backlog nbl_stack[KDTREE_MAX_LEVEL];
  struct kdnode_backlog *top = nbl_stack;

  for ( ; ;)
  {
    if (node != NULL)
    {
      /* Fetch the pop-up backlogged node's sub-id.
       * If not backlogged, fetch the first sub-id. */
      int sub_idx = p_nbl != NULL ? p_nbl->next_sub_idx : KDTREE_RIGHT_INDEX;

      /* Backlog should be left in next loop */
      p_nbl = NULL;

      /* Backlog the node */
      if (is_leaf (node) || (sub_idx == KDTREE_LEFT_INDEX))
      {
        top->node         = NULL;
        top->next_sub_idx = KDTREE_RIGHT_INDEX;
      }
      else
      {
        top->node         = node;
        top->next_sub_idx = KDTREE_LEFT_INDEX;
      }

      top++;
      level++;

      /* Draw lines as long as sub_idx is the first one */
      if (sub_idx == KDTREE_RIGHT_INDEX)
      {
        int i;

        for (i = 1; i < level; i++)
        {
          if (i == level - 1)
          {
            printf ("%-8s", "+-------");
          }
          else
          {
            if (nbl_stack[i - 1].node != NULL)
              printf ("%-8s", "|");
            else
              printf ("%-8s", " ");
          }
        }

        kdnode_dump (node, tree->dim);
      }

      /* Move down according to sub_idx */
      node = sub_idx == KDTREE_LEFT_INDEX ? node->left : node->right;
    }
    else
    {
      p_nbl = top == nbl_stack ? NULL : --top;

      if (p_nbl == NULL)
        /* End of traversal */
        break;

      node = p_nbl->node;
      level--;
    }
  }
}

#endif

