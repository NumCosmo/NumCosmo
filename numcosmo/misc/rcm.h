#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

gint     adj_bandwidth   (gint node_num, gint adj_num, gint adj_row[], gint adj[]);
gboolean adj_contains_ij (gint node_num, gint adj_num, gint adj_row[], gint adj[], gint i, gint j);

void adj_insert_ij ( gint node_num, gint adj_max, gint *adj_num, gint adj_row[], gint adj[], gint i, gint j );
gint adj_perm_bandwidth ( gint node_num, gint adj_num, gint adj_row[], gint adj[], gint perm[], gint perm_inv[] );
void adj_perm_bandwidth_low_up ( gint node_num, gint adj_num, gint adj_row[], gint adj[],
  gint perm[], gint perm_inv[], gint *band_hi, gint *band_lo);
void adj_perm_show ( gint node_num, gint adj_num, gint adj_row[], gint adj[], gint perm[], gint perm_inv[] );
void adj_print ( gint node_num, gint adj_num, gint adj_row[], gint adj[], gchar *title );
void adj_print_some ( gint node_num, gint node_lo, gint node_hi, gint adj_num, gint adj_row[], gint adj[], gchar *title );
void adj_set ( gint node_num, gint adj_max, gint *adj_num, gint adj_row[], gint adj[], gint irow, gint jcol );
void adj_show ( gint node_num, gint adj_num, gint adj_row[], gint adj[] );

void degree ( gint root, gint adj_num, gint adj_row[], gint adj[], gint mask[], gint deg[], gint *iccsze, gint ls[], gint node_num );
void genrcm ( gint node_num, gint adj_num, gint adj_row[], gint adj[], gint perm[] );

void graph_01_adj ( gint node_num, gint adj_num, gint adj_row[], gint adj[] );
void graph_01_size ( gint *node_num, gint *adj_num );

gint i4_max ( gint i1, gint i2 );
gint i4_min ( gint i1, gint i2 );
gint i4_sign ( gint i );
void i4_swap ( gint *i, gint *j );
gint i4_uniform ( gint a, gint b, gint *seed );
gint i4col_compare ( gint m, gint n, gint a[], gint i, gint j );
void i4col_sort_a ( gint m, gint n, gint a[] );
void i4col_swap ( gint m, gint n, gint a[], gint irow1, gint irow2 );
void i4mat_print_some ( gint m, gint n, gint a[], gint ilo, gint jlo, gint ihi, gint jhi, gchar *title );
void i4mat_transpose_print ( gint m, gint n, gint a[], gchar *title );
void i4mat_transpose_print_some ( gint m, gint n, gint a[], gint ilo, gint jlo, gint ihi, gint jhi, gchar *title );
void i4vec_heap_d ( gint n, gint a[] );
gint *i4vec_indicator ( gint n );
void i4vec_print ( gint n, gint a[], gchar *title );
void i4vec_reverse ( gint n, gint a[] );
void i4vec_sort_heap_a ( gint n, gint a[] );
void level_set ( gint root, gint adj_num, gint adj_row[], gint adj[], gint mask[], gint *level_num, gint level_row[], gint level[], gint node_num );
void level_set_print ( gint node_num, gint level_num, gint level_row[], gint level[] );
gboolean perm_check ( gint n, gint p[] );
void perm_inverse3 ( gint n, gint perm[], gint perm_inv[] );
gint *perm_uniform ( gint n, gint *seed );
float r4_abs ( float x );
gint r4_nint ( float x );
void r82vec_permute ( gint n, gdouble a[], gint p[] );
void r8mat_print_some ( gint m, gint n, gdouble a[], gint ilo, gint jlo, gint ihi, gint jhi, gchar *title );
void r8mat_transpose_print_some ( gint m, gint n, gdouble a[], gint ilo, gint jlo, gint ihi, gint jhi, gchar *title );
void rcm ( gint root, gint adj_num, gint adj_row[], gint adj[], gint mask[], gint perm[], gint *iccsze, gint node_num );
void root_find ( gint *root, gint adj_num, gint adj_row[], gint adj[], gint mask[], gint *level_num, gint level_row[], gint level[], gint node_num );
void sort_heap_external ( gint n, gint *indx, gint *i, gint *j, gint isgn );
void timestamp ( );
gint *triangulation_neighbor_triangles ( gint triangle_order, gint triangle_num, gint triangle_node[] );
gint triangulation_order3_adj_count ( gint node_num, gint triangle_num, gint triangle_node[], gint triangle_neighbor[], gint adj_col[] );
gint *triangulation_order3_adj_set ( gint node_num, gint triangle_num, gint triangle_node[], gint triangle_neighbor[], gint adj_num, gint adj_col[] );
void triangulation_order3_example2 ( gint node_num, gint triangle_num, gdouble node_xy[], gint triangle_node[], gint triangle_neighbor[] );
void triangulation_order3_example2_size ( gint *node_num, gint *triangle_num, gint *hole_num );
gint triangulation_order6_adj_count ( gint node_num, gint triangle_num, gint triangle_node[], gint triangle_neighbor[], gint adj_col[] );
gint *triangulation_order6_adj_set ( gint node_num, gint triangle_num, gint triangle_node[], gint triangle_neighbor[], gint adj_num, gint adj_col[] );
void triangulation_order6_example2 ( gint node_num, gint triangle_num, gdouble node_xy[], gint triangle_node[], gint triangle_neighbor[] );
void triangulation_order6_example2_size ( gint *node_num, gint *triangle_num, gint *hole_num );
