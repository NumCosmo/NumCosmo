/* Produced by texiweb from libavl.w. */

/* libavl - library for manipulation of binary trees.
 *  Copyright (C) 1998, 1999, 2000, 2001, 2002, 2004 Free Software
 *  Foundation, Inc.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 *  02110-1301 USA.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include "rb_knn_list.h"
#include "kdtree.h"

int
rb_knn_list_compare (const knn_list_t *rb_a, const knn_list_t *rb_b)
{
  gint ret = (rb_a->distance < rb_b->distance) ? -1 : (rb_a->distance > rb_b->distance);

  if (ret != 0)
    return ret;

  return (rb_a->node->coord_index < rb_b->node->coord_index) ? -1 : (rb_a->node->coord_index > rb_b->node->coord_index);
}

/* Creates and returns a new table
 *  with comparison function |compare| using parameter |param|
 *  and memory allocator |allocator|.
 *  Returns |NULL| if memory allocation failed. */

struct rb_knn_list_table *

rb_knn_list_create (struct libavl_allocator *allocator)
{
  struct rb_knn_list_table *tree;

  if (allocator == NULL)
    allocator = &rb_knn_list_allocator_default;

  tree = allocator->libavl_malloc (allocator, sizeof *tree);

  if (tree == NULL)
    return NULL;

  tree->rb_knn_list_root       = NULL;
  tree->rb_knn_list_alloc      = allocator;
  tree->rb_knn_list_count      = 0;
  tree->rb_knn_list_generation = 0;

  return tree;
}

/* Search |tree| for an item matching |item|, and return it if found.
 *  Otherwise return |NULL|. */
knn_list_t *
rb_knn_list_find (const struct rb_knn_list_table *tree, const knn_list_t *item)
{
  const struct rb_knn_list_node *p;

  assert (tree != NULL && item != NULL);

  for (p = tree->rb_knn_list_root; p != NULL; )
  {
    int cmp = rb_knn_list_compare (item, p->rb_knn_list_data);

    if (cmp < 0)
      p = p->rb_knn_list_link[0];
    else if (cmp > 0)
      p = p->rb_knn_list_link[1];
    else /* |cmp == 0| */
      return p->rb_knn_list_data;
  }

  return NULL;
}

/* Inserts |item| into |tree| and returns a pointer to |item|'s address.
 *  If a duplicate item is found in the tree,
 *  returns a pointer to the duplicate without inserting |item|.
 *  Returns |NULL| in case of memory allocation failure. */
knn_list_t **
rb_knn_list_probe (struct rb_knn_list_table *tree, knn_list_t *item)
{
  struct rb_knn_list_node *pa[RB_MAX_HEIGHT]; /* Nodes on stack. */
  unsigned char da[RB_MAX_HEIGHT];            /* Directions moved from stack nodes. */
  int k;                                      /* Stack height. */

  struct rb_knn_list_node *p; /* Traverses tree looking for insertion point. */
  struct rb_knn_list_node *n; /* Newly inserted node. */

  assert (tree != NULL && item != NULL);

  pa[0] = (struct rb_knn_list_node *) &tree->rb_knn_list_root;
  da[0] = 0;
  k     = 1;

  for (p = tree->rb_knn_list_root; p != NULL; p = p->rb_knn_list_link[da[k - 1]])
  {
    int cmp = rb_knn_list_compare (item, p->rb_knn_list_data);

    if (cmp == 0)
      return &p->rb_knn_list_data;

    pa[k]   = p;
    da[k++] = cmp > 0;
  }

  n = pa[k - 1]->rb_knn_list_link[da[k - 1]] =
    tree->rb_knn_list_alloc->libavl_malloc (tree->rb_knn_list_alloc, sizeof *n);

  if (n == NULL)
    return NULL;

  n->rb_knn_list_data    = item;
  n->rb_knn_list_link[0] = n->rb_knn_list_link[1] = NULL;
  n->rb_knn_list_color   = RB_RED;
  tree->rb_knn_list_count++;
  tree->rb_knn_list_generation++;

  while (k >= 3 && pa[k - 1]->rb_knn_list_color == RB_RED)
  {
    if (da[k - 2] == 0)
    {
      struct rb_knn_list_node *y = pa[k - 2]->rb_knn_list_link[1];

      if ((y != NULL) && (y->rb_knn_list_color == RB_RED))
      {
        pa[k - 1]->rb_knn_list_color = y->rb_knn_list_color = RB_BLACK;
        pa[k - 2]->rb_knn_list_color = RB_RED;
        k                           -= 2;
      }
      else
      {
        struct rb_knn_list_node *x;

        if (da[k - 1] == 0)
        {
          y = pa[k - 1];
        }
        else
        {
          x                              = pa[k - 1];
          y                              = x->rb_knn_list_link[1];
          x->rb_knn_list_link[1]         = y->rb_knn_list_link[0];
          y->rb_knn_list_link[0]         = x;
          pa[k - 2]->rb_knn_list_link[0] = y;
        }

        x                    = pa[k - 2];
        x->rb_knn_list_color = RB_RED;
        y->rb_knn_list_color = RB_BLACK;

        x->rb_knn_list_link[0]                 = y->rb_knn_list_link[1];
        y->rb_knn_list_link[1]                 = x;
        pa[k - 3]->rb_knn_list_link[da[k - 3]] = y;
        break;
      }
    }
    else
    {
      struct rb_knn_list_node *y = pa[k - 2]->rb_knn_list_link[0];

      if ((y != NULL) && (y->rb_knn_list_color == RB_RED))
      {
        pa[k - 1]->rb_knn_list_color = y->rb_knn_list_color = RB_BLACK;
        pa[k - 2]->rb_knn_list_color = RB_RED;
        k                           -= 2;
      }
      else
      {
        struct rb_knn_list_node *x;

        if (da[k - 1] == 1)
        {
          y = pa[k - 1];
        }
        else
        {
          x                              = pa[k - 1];
          y                              = x->rb_knn_list_link[0];
          x->rb_knn_list_link[0]         = y->rb_knn_list_link[1];
          y->rb_knn_list_link[1]         = x;
          pa[k - 2]->rb_knn_list_link[1] = y;
        }

        x                    = pa[k - 2];
        x->rb_knn_list_color = RB_RED;
        y->rb_knn_list_color = RB_BLACK;

        x->rb_knn_list_link[1]                 = y->rb_knn_list_link[0];
        y->rb_knn_list_link[0]                 = x;
        pa[k - 3]->rb_knn_list_link[da[k - 3]] = y;
        break;
      }
    }
  }

  tree->rb_knn_list_root->rb_knn_list_color = RB_BLACK;


  return &n->rb_knn_list_data;
}

/* Inserts |item| into |table|.
 *  Returns |NULL| if |item| was successfully inserted
 *  or if a memory allocation error occurred.
 *  Otherwise, returns the duplicate item. */
knn_list_t *
rb_knn_list_insert (struct rb_knn_list_table *table, knn_list_t *item)
{
  knn_list_t **p = rb_knn_list_probe (table, item);

  return p == NULL || *p == item ? NULL : *p;
}

/* Inserts |item| into |table|, replacing any duplicate item.
 *  Returns |NULL| if |item| was inserted without replacing a duplicate,
 *  or if a memory allocation error occurred.
 *  Otherwise, returns the item that was replaced. */
knn_list_t *
rb_knn_list_replace (struct rb_knn_list_table *table, knn_list_t *item)
{
  knn_list_t **p = rb_knn_list_probe (table, item);

  if ((p == NULL) || (*p == item))
  {
    return NULL;
  }
  else
  {
    knn_list_t *r = *p;

    *p = item;

    return r;
  }
}

/* Deletes from |tree| and returns an item matching |item|.
 *  Returns a null pointer if no matching item found. */
knn_list_t *
rb_knn_list_delete (struct rb_knn_list_table *tree, const knn_list_t *item)
{
  struct rb_knn_list_node *pa[RB_MAX_HEIGHT]; /* Nodes on stack. */
  unsigned char da[RB_MAX_HEIGHT];            /* Directions moved from stack nodes. */
  int k;                                      /* Stack height. */

  struct rb_knn_list_node *p; /* The node to delete, or a node part way to it. */
  int cmp;                    /* Result of comparison between |item| and |p|. */

  assert (tree != NULL && item != NULL);

  k = 0;
  p = (struct rb_knn_list_node *) &tree->rb_knn_list_root;

  for (cmp = -1; cmp != 0;
       cmp = rb_knn_list_compare (item, p->rb_knn_list_data))
  {
    int dir = cmp > 0;

    pa[k]   = p;
    da[k++] = dir;

    p = p->rb_knn_list_link[dir];

    if (p == NULL)
      return NULL;
  }

  item = p->rb_knn_list_data;

  if (p->rb_knn_list_link[1] == NULL)
  {
    pa[k - 1]->rb_knn_list_link[da[k - 1]] = p->rb_knn_list_link[0];
  }
  else
  {
    enum rb_knn_list_color t;
    struct rb_knn_list_node *r = p->rb_knn_list_link[1];

    if (r->rb_knn_list_link[0] == NULL)
    {
      r->rb_knn_list_link[0]                 = p->rb_knn_list_link[0];
      t                                      = r->rb_knn_list_color;
      r->rb_knn_list_color                   = p->rb_knn_list_color;
      p->rb_knn_list_color                   = t;
      pa[k - 1]->rb_knn_list_link[da[k - 1]] = r;
      da[k]                                  = 1;
      pa[k++]                                = r;
    }
    else
    {
      struct rb_knn_list_node *s;
      int j = k++;

      for ( ;;)
      {
        da[k]   = 0;
        pa[k++] = r;
        s       = r->rb_knn_list_link[0];

        if (s->rb_knn_list_link[0] == NULL)
          break;

        r = s;
      }

      da[j]                                  = 1;
      pa[j]                                  = s;
      pa[j - 1]->rb_knn_list_link[da[j - 1]] = s;

      s->rb_knn_list_link[0] = p->rb_knn_list_link[0];
      r->rb_knn_list_link[0] = s->rb_knn_list_link[1];
      s->rb_knn_list_link[1] = p->rb_knn_list_link[1];

      t                    = s->rb_knn_list_color;
      s->rb_knn_list_color = p->rb_knn_list_color;
      p->rb_knn_list_color = t;
    }
  }

  if (p->rb_knn_list_color == RB_BLACK)
  {
    for ( ;;)
    {
      struct rb_knn_list_node *x = pa[k - 1]->rb_knn_list_link[da[k - 1]];

      if ((x != NULL) && (x->rb_knn_list_color == RB_RED))
      {
        x->rb_knn_list_color = RB_BLACK;
        break;
      }

      if (k < 2)
        break;

      if (da[k - 1] == 0)
      {
        struct rb_knn_list_node *w = pa[k - 1]->rb_knn_list_link[1];

        if (w->rb_knn_list_color == RB_RED)
        {
          w->rb_knn_list_color         = RB_BLACK;
          pa[k - 1]->rb_knn_list_color = RB_RED;

          pa[k - 1]->rb_knn_list_link[1]         = w->rb_knn_list_link[0];
          w->rb_knn_list_link[0]                 = pa[k - 1];
          pa[k - 2]->rb_knn_list_link[da[k - 2]] = w;

          pa[k]     = pa[k - 1];
          da[k]     = 0;
          pa[k - 1] = w;
          k++;

          w = pa[k - 1]->rb_knn_list_link[1];
        }

        if (((w->rb_knn_list_link[0] == NULL)
             || (w->rb_knn_list_link[0]->rb_knn_list_color == RB_BLACK))
            && ((w->rb_knn_list_link[1] == NULL)
                || (w->rb_knn_list_link[1]->rb_knn_list_color == RB_BLACK)))
        {
          w->rb_knn_list_color = RB_RED;
        }
        else
        {
          if ((w->rb_knn_list_link[1] == NULL)
              || (w->rb_knn_list_link[1]->rb_knn_list_color == RB_BLACK))
          {
            struct rb_knn_list_node *y = w->rb_knn_list_link[0];

            y->rb_knn_list_color   = RB_BLACK;
            w->rb_knn_list_color   = RB_RED;
            w->rb_knn_list_link[0] = y->rb_knn_list_link[1];
            y->rb_knn_list_link[1] = w;
            w                      = pa[k - 1]->rb_knn_list_link[1] = y;
          }

          w->rb_knn_list_color                      = pa[k - 1]->rb_knn_list_color;
          pa[k - 1]->rb_knn_list_color              = RB_BLACK;
          w->rb_knn_list_link[1]->rb_knn_list_color = RB_BLACK;

          pa[k - 1]->rb_knn_list_link[1]         = w->rb_knn_list_link[0];
          w->rb_knn_list_link[0]                 = pa[k - 1];
          pa[k - 2]->rb_knn_list_link[da[k - 2]] = w;
          break;
        }
      }
      else
      {
        struct rb_knn_list_node *w = pa[k - 1]->rb_knn_list_link[0];

        if (w->rb_knn_list_color == RB_RED)
        {
          w->rb_knn_list_color         = RB_BLACK;
          pa[k - 1]->rb_knn_list_color = RB_RED;

          pa[k - 1]->rb_knn_list_link[0]         = w->rb_knn_list_link[1];
          w->rb_knn_list_link[1]                 = pa[k - 1];
          pa[k - 2]->rb_knn_list_link[da[k - 2]] = w;

          pa[k]     = pa[k - 1];
          da[k]     = 1;
          pa[k - 1] = w;
          k++;

          w = pa[k - 1]->rb_knn_list_link[0];
        }

        if (((w->rb_knn_list_link[0] == NULL)
             || (w->rb_knn_list_link[0]->rb_knn_list_color == RB_BLACK))
            && ((w->rb_knn_list_link[1] == NULL)
                || (w->rb_knn_list_link[1]->rb_knn_list_color == RB_BLACK)))
        {
          w->rb_knn_list_color = RB_RED;
        }
        else
        {
          if ((w->rb_knn_list_link[0] == NULL)
              || (w->rb_knn_list_link[0]->rb_knn_list_color == RB_BLACK))
          {
            struct rb_knn_list_node *y = w->rb_knn_list_link[1];

            y->rb_knn_list_color   = RB_BLACK;
            w->rb_knn_list_color   = RB_RED;
            w->rb_knn_list_link[1] = y->rb_knn_list_link[0];
            y->rb_knn_list_link[0] = w;
            w                      = pa[k - 1]->rb_knn_list_link[0] = y;
          }

          w->rb_knn_list_color                      = pa[k - 1]->rb_knn_list_color;
          pa[k - 1]->rb_knn_list_color              = RB_BLACK;
          w->rb_knn_list_link[0]->rb_knn_list_color = RB_BLACK;

          pa[k - 1]->rb_knn_list_link[0]         = w->rb_knn_list_link[1];
          w->rb_knn_list_link[1]                 = pa[k - 1];
          pa[k - 2]->rb_knn_list_link[da[k - 2]] = w;
          break;
        }
      }

      k--;
    }
  }

  tree->rb_knn_list_alloc->libavl_free (tree->rb_knn_list_alloc, p);
  tree->rb_knn_list_count--;
  tree->rb_knn_list_generation++;

  return (knn_list_t *) item;
}

/* Refreshes the stack of parent pointers in |trav|
 *  and updates its generation number. */
static void
trav_refresh (struct rb_knn_list_traverser *trav)
{
  assert (trav != NULL);

  trav->rb_knn_list_generation = trav->rb_knn_list_table->rb_knn_list_generation;

  if (trav->rb_knn_list_node != NULL)
  {
    struct rb_knn_list_node *node = trav->rb_knn_list_node;
    struct rb_knn_list_node *i;

    trav->rb_knn_list_height = 0;

    for (i = trav->rb_knn_list_table->rb_knn_list_root; i != node; )
    {
      assert (trav->rb_knn_list_height < RB_MAX_HEIGHT);
      assert (i != NULL);

      trav->rb_knn_list_stack[trav->rb_knn_list_height++] = i;
      i                                                   = i->rb_knn_list_link[rb_knn_list_compare (
                                                                                  node->rb_knn_list_data,
                                                                                  i->rb_knn_list_data) > 0];
    }
  }
}

/* Initializes |trav| for use with |tree|
 *  and selects the null node. */
void
rb_knn_list_t_init (struct rb_knn_list_traverser *trav, struct rb_knn_list_table *tree)
{
  trav->rb_knn_list_table      = tree;
  trav->rb_knn_list_node       = NULL;
  trav->rb_knn_list_height     = 0;
  trav->rb_knn_list_generation = tree->rb_knn_list_generation;
}

/* Initializes |trav| for |tree|
 *  and selects and returns a pointer to its least-valued item.
 *  Returns |NULL| if |tree| contains no nodes. */
knn_list_t *
rb_knn_list_t_first (struct rb_knn_list_traverser *trav, struct rb_knn_list_table *tree)
{
  struct rb_knn_list_node *x;

  assert (tree != NULL && trav != NULL);

  trav->rb_knn_list_table      = tree;
  trav->rb_knn_list_height     = 0;
  trav->rb_knn_list_generation = tree->rb_knn_list_generation;

  x = tree->rb_knn_list_root;

  if (x != NULL)
    while (x->rb_knn_list_link[0] != NULL)
    {
      assert (trav->rb_knn_list_height < RB_MAX_HEIGHT);
      trav->rb_knn_list_stack[trav->rb_knn_list_height++] = x;
      x                                                   = x->rb_knn_list_link[0];
    }

  trav->rb_knn_list_node = x;

  return x != NULL ? x->rb_knn_list_data : NULL;
}

/* Initializes |trav| for |tree|
 *  and selects and returns a pointer to its greatest-valued item.
 *  Returns |NULL| if |tree| contains no nodes. */
knn_list_t *
rb_knn_list_t_last (struct rb_knn_list_traverser *trav, struct rb_knn_list_table *tree)
{
  struct rb_knn_list_node *x;

  assert (tree != NULL && trav != NULL);

  trav->rb_knn_list_table      = tree;
  trav->rb_knn_list_height     = 0;
  trav->rb_knn_list_generation = tree->rb_knn_list_generation;

  x = tree->rb_knn_list_root;

  if (x != NULL)
    while (x->rb_knn_list_link[1] != NULL)
    {
      assert (trav->rb_knn_list_height < RB_MAX_HEIGHT);
      trav->rb_knn_list_stack[trav->rb_knn_list_height++] = x;
      x                                                   = x->rb_knn_list_link[1];
    }

  trav->rb_knn_list_node = x;

  return x != NULL ? x->rb_knn_list_data : NULL;
}

/* Searches for |item| in |tree|.
 *  If found, initializes |trav| to the item found and returns the item
 *  as well.
 *  If there is no matching item, initializes |trav| to the null item
 *  and returns |NULL|. */
knn_list_t *
rb_knn_list_t_find (struct rb_knn_list_traverser *trav, struct rb_knn_list_table *tree, knn_list_t *item)
{
  struct rb_knn_list_node *p, *q;

  assert (trav != NULL && tree != NULL && item != NULL);
  trav->rb_knn_list_table      = tree;
  trav->rb_knn_list_height     = 0;
  trav->rb_knn_list_generation = tree->rb_knn_list_generation;

  for (p = tree->rb_knn_list_root; p != NULL; p = q)
  {
    int cmp = rb_knn_list_compare (item, p->rb_knn_list_data);

    if (cmp < 0)
    {
      q = p->rb_knn_list_link[0];
    }
    else if (cmp > 0)
    {
      q = p->rb_knn_list_link[1];
    }
    else /* |cmp == 0| */
    {
      trav->rb_knn_list_node = p;

      return p->rb_knn_list_data;
    }

    assert (trav->rb_knn_list_height < RB_MAX_HEIGHT);
    trav->rb_knn_list_stack[trav->rb_knn_list_height++] = p;
  }

  trav->rb_knn_list_height = 0;
  trav->rb_knn_list_node   = NULL;

  return NULL;
}

/* Attempts to insert |item| into |tree|.
 *  If |item| is inserted successfully, it is returned and |trav| is
 *  initialized to its location.
 *  If a duplicate is found, it is returned and |trav| is initialized to
 *  its location.  No replacement of the item occurs.
 *  If a memory allocation failure occurs, |NULL| is returned and |trav|
 *  is initialized to the null item. */
knn_list_t *
rb_knn_list_t_insert (struct rb_knn_list_traverser *trav, struct rb_knn_list_table *tree, knn_list_t *item)
{
  knn_list_t **p;

  assert (trav != NULL && tree != NULL && item != NULL);

  p = rb_knn_list_probe (tree, item);

  if (p != NULL)
  {
    trav->rb_knn_list_table = tree;
    trav->rb_knn_list_node  =
      ((struct rb_knn_list_node *)
       ((char *) p - offsetof (struct rb_knn_list_node, rb_knn_list_data)));
    trav->rb_knn_list_generation = tree->rb_knn_list_generation - 1;

    return *p;
  }
  else
  {
    rb_knn_list_t_init (trav, tree);

    return NULL;
  }
}

/* Initializes |trav| to have the same current node as |src|. */
knn_list_t *
rb_knn_list_t_copy (struct rb_knn_list_traverser *trav, const struct rb_knn_list_traverser *src)
{
  assert (trav != NULL && src != NULL);

  if (trav != src)
  {
    trav->rb_knn_list_table      = src->rb_knn_list_table;
    trav->rb_knn_list_node       = src->rb_knn_list_node;
    trav->rb_knn_list_generation = src->rb_knn_list_generation;

    if (trav->rb_knn_list_generation == trav->rb_knn_list_table->rb_knn_list_generation)
    {
      trav->rb_knn_list_height = src->rb_knn_list_height;
      memcpy (trav->rb_knn_list_stack, (const void *) src->rb_knn_list_stack,
              sizeof *trav->rb_knn_list_stack * trav->rb_knn_list_height);
    }
  }

  return trav->rb_knn_list_node != NULL ? trav->rb_knn_list_node->rb_knn_list_data : NULL;
}

/* Returns the next data item in inorder
 *  within the tree being traversed with |trav|,
 *  or if there are no more data items returns |NULL|. */
knn_list_t *
rb_knn_list_t_next (struct rb_knn_list_traverser *trav)
{
  struct rb_knn_list_node *x;

  assert (trav != NULL);

  if (trav->rb_knn_list_generation != trav->rb_knn_list_table->rb_knn_list_generation)
    trav_refresh (trav);

  x = trav->rb_knn_list_node;

  if (x == NULL)
  {
    return rb_knn_list_t_first (trav, trav->rb_knn_list_table);
  }
  else if (x->rb_knn_list_link[1] != NULL)
  {
    assert (trav->rb_knn_list_height < RB_MAX_HEIGHT);
    trav->rb_knn_list_stack[trav->rb_knn_list_height++] = x;
    x                                                   = x->rb_knn_list_link[1];

    while (x->rb_knn_list_link[0] != NULL)
    {
      assert (trav->rb_knn_list_height < RB_MAX_HEIGHT);
      trav->rb_knn_list_stack[trav->rb_knn_list_height++] = x;
      x                                                   = x->rb_knn_list_link[0];
    }
  }
  else
  {
    struct rb_knn_list_node *y;

    do {
      if (trav->rb_knn_list_height == 0)
      {
        trav->rb_knn_list_node = NULL;

        return NULL;
      }

      y = x;
      x = trav->rb_knn_list_stack[--trav->rb_knn_list_height];
    } while (y == x->rb_knn_list_link[1]);
  }

  trav->rb_knn_list_node = x;

  return x->rb_knn_list_data;
}

/* Returns the previous data item in inorder
 *  within the tree being traversed with |trav|,
 *  or if there are no more data items returns |NULL|. */
knn_list_t *
rb_knn_list_t_prev (struct rb_knn_list_traverser *trav)
{
  struct rb_knn_list_node *x;

  assert (trav != NULL);

  if (trav->rb_knn_list_generation != trav->rb_knn_list_table->rb_knn_list_generation)
    trav_refresh (trav);

  x = trav->rb_knn_list_node;

  if (x == NULL)
  {
    return rb_knn_list_t_last (trav, trav->rb_knn_list_table);
  }
  else if (x->rb_knn_list_link[0] != NULL)
  {
    assert (trav->rb_knn_list_height < RB_MAX_HEIGHT);
    trav->rb_knn_list_stack[trav->rb_knn_list_height++] = x;
    x                                                   = x->rb_knn_list_link[0];

    while (x->rb_knn_list_link[1] != NULL)
    {
      assert (trav->rb_knn_list_height < RB_MAX_HEIGHT);
      trav->rb_knn_list_stack[trav->rb_knn_list_height++] = x;
      x                                                   = x->rb_knn_list_link[1];
    }
  }
  else
  {
    struct rb_knn_list_node *y;

    do {
      if (trav->rb_knn_list_height == 0)
      {
        trav->rb_knn_list_node = NULL;

        return NULL;
      }

      y = x;
      x = trav->rb_knn_list_stack[--trav->rb_knn_list_height];
    } while (y == x->rb_knn_list_link[0]);
  }

  trav->rb_knn_list_node = x;

  return x->rb_knn_list_data;
}

/* Returns |trav|'s current item. */
knn_list_t *
rb_knn_list_t_cur (struct rb_knn_list_traverser *trav)
{
  assert (trav != NULL);

  return trav->rb_knn_list_node != NULL ? trav->rb_knn_list_node->rb_knn_list_data : NULL;
}

/* Replaces the current item in |trav| by |new| and returns the item replaced.
 |trav| must not have the null item selected.
 *  The new item must not upset the ordering of the tree. */
knn_list_t *
rb_knn_list_t_replace (struct rb_knn_list_traverser *trav, knn_list_t *new)
{
  knn_list_t *old;

  assert (trav != NULL && trav->rb_knn_list_node != NULL && new != NULL);
  old                                      = trav->rb_knn_list_node->rb_knn_list_data;
  trav->rb_knn_list_node->rb_knn_list_data = new;

  return old;
}

/* Destroys |new| with |rb_knn_list_destroy (new, destroy)|,
 *  first setting right links of nodes in |stack| within |new|
 *  to null pointers to avoid touching uninitialized data. */
static void
copy_error_recovery (struct rb_knn_list_node **stack, int height,
                     struct rb_knn_list_table *new)
{
  assert (stack != NULL && height >= 0 && new != NULL);

  for ( ; height > 2; height -= 2)
    stack[height - 1]->rb_knn_list_link[1] = NULL;

  rb_knn_list_destroy (new);
}

/* Copies |org| to a newly created tree, which is returned.
 *  If |copy != NULL|, each data item in |org| is first passed to |copy|,
 *  and the return values are inserted into the tree,
 *  with |NULL| return values taken as indications of failure.
 *  On failure, destroys the partially created new tree,
 *  applying |destroy|, if non-null, to each item in the new tree so far,
 *  and returns |NULL|.
 *  If |allocator != NULL|, it is used for allocation in the new tree.
 *  Otherwise, the same allocator used for |org| is used. */
struct rb_knn_list_table *

rb_knn_list_copy (const struct rb_knn_list_table *org, rb_knn_list_copy_func *copy,
                  struct libavl_allocator *allocator)
{
  struct rb_knn_list_node *stack[2 * (RB_MAX_HEIGHT + 1)];
  int height = 0;

  struct rb_knn_list_table *new;
  const struct rb_knn_list_node *x;
  struct rb_knn_list_node *y;

  assert (org != NULL);
  new = rb_knn_list_create (allocator != NULL ? allocator : org->rb_knn_list_alloc);

  if (new == NULL)
    return NULL;

  new->rb_knn_list_count = org->rb_knn_list_count;

  if (new->rb_knn_list_count == 0)
    return new;

  x = (const struct rb_knn_list_node *) &org->rb_knn_list_root;
  y = (struct rb_knn_list_node *) &new->rb_knn_list_root;

  for ( ;;)
  {
    while (x->rb_knn_list_link[0] != NULL)
    {
      assert (height < 2 * (RB_MAX_HEIGHT + 1));

      y->rb_knn_list_link[0] =
        new->rb_knn_list_alloc->libavl_malloc (new->rb_knn_list_alloc,
                                               sizeof *y->rb_knn_list_link[0]);

      if (y->rb_knn_list_link[0] == NULL)
      {
        if (y != (struct rb_knn_list_node *) &new->rb_knn_list_root)
        {
          y->rb_knn_list_data    = NULL;
          y->rb_knn_list_link[1] = NULL;
        }

        copy_error_recovery (stack, height, new);

        return NULL;
      }

      stack[height++] = (struct rb_knn_list_node *) x;
      stack[height++] = y;
      x               = x->rb_knn_list_link[0];
      y               = y->rb_knn_list_link[0];
    }

    y->rb_knn_list_link[0] = NULL;

    for ( ;;)
    {
      y->rb_knn_list_color = x->rb_knn_list_color;

      if (copy == NULL)
      {
        y->rb_knn_list_data = x->rb_knn_list_data;
      }
      else
      {
        y->rb_knn_list_data = copy (x->rb_knn_list_data);

        if (y->rb_knn_list_data == NULL)
        {
          y->rb_knn_list_link[1] = NULL;
          copy_error_recovery (stack, height, new);

          return NULL;
        }
      }

      if (x->rb_knn_list_link[1] != NULL)
      {
        y->rb_knn_list_link[1] =
          new->rb_knn_list_alloc->libavl_malloc (new->rb_knn_list_alloc,
                                                 sizeof *y->rb_knn_list_link[1]);

        if (y->rb_knn_list_link[1] == NULL)
        {
          copy_error_recovery (stack, height, new);

          return NULL;
        }

        x = x->rb_knn_list_link[1];
        y = y->rb_knn_list_link[1];
        break;
      }
      else
      {
        y->rb_knn_list_link[1] = NULL;
      }

      if (height <= 2)
        return new;

      y = stack[--height];
      x = stack[--height];
    }
  }
}

/* Frees storage allocated for |tree|.
 *  If |destroy != NULL|, applies it to each data item in inorder. */
void
rb_knn_list_destroy (struct rb_knn_list_table *tree)
{
  struct rb_knn_list_node *p, *q;

  assert (tree != NULL);

  for (p = tree->rb_knn_list_root; p != NULL; p = q)
    if (p->rb_knn_list_link[0] == NULL)
    {
      q = p->rb_knn_list_link[1];

      if (p->rb_knn_list_data != NULL)
        g_slice_free (knn_list_t, p->rb_knn_list_data);

      tree->rb_knn_list_alloc->libavl_free (tree->rb_knn_list_alloc, p);
    }
    else
    {
      q                      = p->rb_knn_list_link[0];
      p->rb_knn_list_link[0] = q->rb_knn_list_link[1];
      q->rb_knn_list_link[1] = p;
    }

  tree->rb_knn_list_alloc->libavl_free (tree->rb_knn_list_alloc, tree);
}

/* Allocates |size| bytes of space using |malloc()|.
 *  Returns a null pointer if allocation fails. */
void *
rb_knn_list_malloc (struct libavl_allocator *allocator, size_t size)
{
  assert (allocator != NULL && size > 0);

  return malloc (size);
}

/* Frees |block|. */
void
rb_knn_list_free (struct libavl_allocator *allocator, void *block)
{
  assert (allocator != NULL && block != NULL);
  free (block);
}

/* Default memory allocator that uses |malloc()| and |free()|. */
struct libavl_allocator rb_knn_list_allocator_default =
{
  rb_knn_list_malloc,
  rb_knn_list_free
};

#undef NDEBUG
#include <assert.h>

/* Asserts that |rb_knn_list_insert()| succeeds at inserting |item| into |table|. */
void
rb_knn_list_assert_insert (struct rb_knn_list_table *table, knn_list_t *item)
{
  knn_list_t **p = rb_knn_list_probe (table, item);

  assert (p != NULL && *p == item);
}

/* Asserts that |rb_knn_list_delete()| really removes |item| from |table|,
 *  and returns the removed item. */
knn_list_t *
rb_knn_list_assert_delete (struct rb_knn_list_table *table, knn_list_t *item)
{
  knn_list_t *p = rb_knn_list_delete (table, item);

  assert (p != NULL);

  return p;
}

