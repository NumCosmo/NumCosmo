/*
 * Copyright (C) 2017, Leo Ma <begeekmyfriend@gmail.com>
 */

#ifndef _KD_TREE_H
#define _KD_TREE_H

#define KDTREE_MAX_LEVEL 64
#define KDTREE_LEFT_INDEX 0
#define KDTREE_RIGHT_INDEX 1

typedef struct knn_list
{
  struct kdnode *node;
  double distance;
} knn_list_t;

typedef struct kdnode
{
  long coord_index;
  double *coord;
  struct kdnode *left;
  struct kdnode *right;
  int r;
} kdnode_t;

typedef struct kdtree
{
  struct kdnode *root;
  size_t count;
  size_t capacity;
  double *coords;
  double **coord_table;
  long *coord_indexes;
  long *coord_indexes_tmp;
  int dim;
} kdtree_t;

struct kdtree *kdtree_init (int dim);

void kdtree_insert (struct kdtree *tree, double *coord);
void kdtree_rebuild (struct kdtree *tree);
void *kdtree_knn_search (struct kdtree *tree, double *coord, int k);
void kdtree_destroy (struct kdtree *tree);
void kdtree_dump (struct kdtree *tree);

#endif /* _KD_TREE_H */

