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

#ifndef RB_KNN_LIST_H
#define RB_KNN_LIST_H 1

#include <stddef.h>
#include "kdtree.h"

/* Function types. */
typedef knn_list_t *rb_knn_list_copy_func (knn_list_t *rb_knn_list_item);

#ifndef LIBAVL_ALLOCATOR
#define LIBAVL_ALLOCATOR
/* Memory allocator. */
struct libavl_allocator
{
  void *(*libavl_malloc) (struct libavl_allocator *, size_t libavl_size);
  void (*libavl_free) (struct libavl_allocator *, void *libavl_block);
};
#endif

/* Default memory allocator. */
extern struct libavl_allocator rb_knn_list_allocator_default;
void *rb_knn_list_malloc (struct libavl_allocator *, size_t);
void rb_knn_list_free (struct libavl_allocator *, void *);

/* Maximum RB height. */
#ifndef RB_MAX_HEIGHT
#define RB_MAX_HEIGHT 128
#endif

/* Tree data structure. */
typedef struct rb_knn_list_table
{
  struct rb_knn_list_node *rb_knn_list_root;  /* Tree's root. */
  struct libavl_allocator *rb_knn_list_alloc; /* Memory allocator. */
  size_t rb_knn_list_count;                   /* Number of items in tree. */
  unsigned long rb_knn_list_generation;       /* Generation number. */
} rb_knn_list_table_t;

/* Color of a red-black node. */
enum rb_knn_list_color
{
  RB_BLACK, /* Black. */
  RB_RED    /* Red. */
};

/* A red-black tree node. */
struct rb_knn_list_node
{
  struct rb_knn_list_node *rb_knn_list_link[2]; /* Subtrees. */
  knn_list_t *rb_knn_list_data;                 /* Pointer to data. */
  unsigned char rb_knn_list_color;              /* Color. */
};

/* RB traverser structure. */
typedef struct rb_knn_list_traverser
{
  struct rb_knn_list_table *rb_knn_list_table; /* Tree being traversed. */
  struct rb_knn_list_node *rb_knn_list_node;   /* Current node in tree. */
  struct rb_knn_list_node *rb_knn_list_stack[RB_MAX_HEIGHT];
  /* All the nodes above |rb_knn_list_node|. */
  size_t rb_knn_list_height;            /* Number of nodes in |rb_knn_list_parent|. */
  unsigned long rb_knn_list_generation; /* Generation number. */
} rb_knn_list_traverser_t;

/* Table functions. */
struct rb_knn_list_table *rb_knn_list_create (struct libavl_allocator *);
struct rb_knn_list_table *rb_knn_list_copy (const struct rb_knn_list_table *, rb_knn_list_copy_func *,
                                            struct libavl_allocator *);
void rb_knn_list_destroy (struct rb_knn_list_table *);
knn_list_t **rb_knn_list_probe (struct rb_knn_list_table *, knn_list_t *);
knn_list_t *rb_knn_list_insert (struct rb_knn_list_table *, knn_list_t *);
knn_list_t *rb_knn_list_replace (struct rb_knn_list_table *, knn_list_t *);
knn_list_t *rb_knn_list_delete (struct rb_knn_list_table *, const knn_list_t *);
knn_list_t *rb_knn_list_find (const struct rb_knn_list_table *, const knn_list_t *);
void rb_knn_list_assert_insert (struct rb_knn_list_table *, knn_list_t *);
knn_list_t *rb_knn_list_assert_delete (struct rb_knn_list_table *, knn_list_t *);

#define rb_knn_list_count(table) ((size_t) (table)->rb_knn_list_count)

/* Table traverser functions. */
void rb_knn_list_t_init (struct rb_knn_list_traverser *, struct rb_knn_list_table *);
knn_list_t *rb_knn_list_t_first (struct rb_knn_list_traverser *, struct rb_knn_list_table *);
knn_list_t *rb_knn_list_t_last (struct rb_knn_list_traverser *, struct rb_knn_list_table *);
knn_list_t *rb_knn_list_t_find (struct rb_knn_list_traverser *, struct rb_knn_list_table *, knn_list_t *);
knn_list_t *rb_knn_list_t_insert (struct rb_knn_list_traverser *, struct rb_knn_list_table *, knn_list_t *);
knn_list_t *rb_knn_list_t_copy (struct rb_knn_list_traverser *, const struct rb_knn_list_traverser *);
knn_list_t *rb_knn_list_t_next (struct rb_knn_list_traverser *);
knn_list_t *rb_knn_list_t_prev (struct rb_knn_list_traverser *);
knn_list_t *rb_knn_list_t_cur (struct rb_knn_list_traverser *);
knn_list_t *rb_knn_list_t_replace (struct rb_knn_list_traverser *, knn_list_t *);

#endif /* rb.h */

