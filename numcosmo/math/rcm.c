
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include "math/rcm.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//****************************************************************************80

gint adj_bandwidth ( gint node_num, gint adj_num, gint adj_row[], gint adj[] )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_BANDWIDTH computes the bandwidth of an adjacency matrix.
//
//  Discussion:
//
//    Thanks to Man Yuan, of Southeast University, China, for pointing out
//    an inconsistency in the indexing of ADJ_ROW, 06 June 2011.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint ADJ_NUM, the number of adjacency entries.
//
//    Input, gint ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, gint ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Output, gint ADJ_BANDWIDTH, the bandwidth of the adjacency
//    matrix.
//
{
  gint band_hi;
  gint band_lo;
  gint col;
  gint i;
  gint j;
  gint value;

  band_lo = 0;
  band_hi = 0;

  for ( i = 0; i < node_num; i++ )
  {
    for ( j = adj_row[i]; j <= adj_row[i+1]-1; j++ )
    {
      col = adj[j-1] - 1;
      band_lo = i4_max ( band_lo, i - col );
      band_hi = i4_max ( band_hi, col - i );
    }
  }

  value = band_lo + 1 + band_hi;

  return value;
}
//****************************************************************************80

gboolean adj_contains_ij ( gint node_num, gint adj_num, gint adj_row[], gint adj[],
  gint i, gint j )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_CONTAINS_IJ determines if (I,J) is in an adjacency structure.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint ADJ_NUM, the number of adjacency entries.
//
//    Input, gint ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, gint ADJ[ADJ_NUM], the adjacency structure.
//
//    Input, gint I, J, the two nodes, for which we want to know
//    whether I is adjacent to J.
//
//    Output, gboolean ADJ_CONTAINS_IJ, is TRUE if I = J, or the adjacency
//    structure contains the information that I is adjacent to J.
//
{
  gint k;
  gint khi;
  gint klo;
  gboolean value;
//
//  Symmetric entries are not stored.
//
  if ( i == j )
  {
    value = TRUE;
    return value;
  }
//
//  Illegal I, J entries.
//
  if ( node_num < i )
  {
    value = FALSE;
    return value;
  }
  else if ( i < 1 )
  {
    value = FALSE;
    return value;
  }
  else if ( node_num < j )
  {
    value = FALSE;
    return value;
  }
  else if ( j < 1 )
  {
    value = FALSE;
    return value;
  }
//
//  Search the adjacency entries already stored for row I,
//  to see if J has already been stored.
//
  klo = adj_row[i-1];
  khi = adj_row[i]-1;

  for ( k = klo; k <= khi; k++ )
  {
    if ( adj[k-1] == j )
    {
      value = TRUE;
      return value;
    }
  }
  value = FALSE;

  return value;
}
//****************************************************************************80

void adj_insert_ij ( gint node_num, gint adj_max, gint *adj_num, gint adj_row[],
  gint adj[], gint i, gint j )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_INSERT_IJ inserts (I,J) into an adjacency structure.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint ADJ_MAX, the maximum number of adjacency entries.
//
//    Input/output, gint ADJ_NUM, the number of adjacency entries.
//
//    Input/output, gint ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input/output, gint ADJ[ADJ_NUM], the adjacency structure.
//
//    Input, gint I, J, the two nodes which are adjacent.
//
{
  gint j_spot;
  gint k;
//
//  A new adjacency entry must be made.
//  Check that we're not exceeding the storage allocation for ADJ.
//
  if ( adj_max < *adj_num + 1 )
  {
    g_error ("ADJ_INSERT_IJ - Fatal error!\n"
             "  All available storage has been used.\n"
             "  No more information can be stored!\n"
             "  This error occurred for \n"
             "  Row I =    %d\n"
             "  Column J = %d", i, j);
  }
//
//  The action is going to occur between ADJ_ROW(I) and ADJ_ROW(I+1)-1:
//
  j_spot = adj_row[i-1];

  for ( k = adj_row[i-1]; k <= adj_row[i]-1; k++ )
  {
    if ( adj[k-1] == j )
    {
      return;
    }
    else if ( adj[k-1] < j )
    {
      j_spot = k + 1;
    }
    else
    {
      break;
    }
  }

  for ( k = *adj_num; j_spot <= k; k-- )
  {
    adj[k] = adj[k-1];
  }
  adj[j_spot-1] = j;

  for ( k = i; k <= node_num; k++ )
  {
    adj_row[k] = adj_row[k] + 1;
  }

  *adj_num = *adj_num + 1;

  return;
}
//****************************************************************************80

gint adj_perm_bandwidth ( gint node_num, gint adj_num, gint adj_row[], gint adj[],
  gint perm[], gint perm_inv[] )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_PERM_BANDWIDTH computes the bandwidth of a permuted adjacency matrix.
//
//  Discussion:
//
//    The matrix is defined by the adjacency information and a permutation.
//
//    The routine also computes the bandwidth and the size of the envelope.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint ADJ_NUM, the number of adjacency entries.
//
//    Input, gint ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, gint ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input, gint PERM[NODE_NUM], PERM_INV(NODE_NUM), the permutation
//    and inverse permutation.
//
//    Output, gint ADJ_PERM_BANDWIDTH, the bandwidth of the permuted
//    adjacency matrix.
//
{
  gint band_hi;
  gint band_lo;
  gint bandwidth;
  gint col;
  gint i;
  gint j;

  band_lo = 0;
  band_hi = 0;

  for ( i = 0; i < node_num; i++ )
  {
    for ( j = adj_row[perm[i]-1]; j <= adj_row[perm[i]]-1; j++ )
    {
      col = perm_inv[adj[j-1]-1];
      band_lo = i4_max ( band_lo, i - col );
      band_hi = i4_max ( band_hi, col - i );
    }
  }

  bandwidth = band_lo + 1 + band_hi;

  return bandwidth;
}
//****************************************************************************80

void adj_perm_bandwidth_low_up ( gint node_num, gint adj_num, gint adj_row[], gint adj[],
  gint perm[], gint perm_inv[], gint *band_hi, gint *band_lo)
{
  gint col;
  gint i;
  gint j;

  *band_hi = 0;
  *band_lo = 0;

  for ( i = 0; i < node_num; i++ )
  {
    for ( j = adj_row[perm[i]-1]; j <= adj_row[perm[i]]-1; j++ )
    {
      col = perm_inv[adj[j-1]-1];
      *band_lo = i4_max ( *band_lo, i - col );
      *band_hi = i4_max ( *band_hi, col - i );
    }
  }
}

void adj_perm_show ( gint node_num, gint adj_num, gint adj_row[], gint adj[],
  gint perm[], gint perm_inv[] )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_PERM_SHOW displays a symbolic picture of a permuted adjacency matrix.
//
//  Discussion:
//
//    The matrix is defined by the adjacency information and a permutation.
//
//    The routine also computes the bandwidth and the size of the envelope.
//
//    If no permutation has been done, you must set PERM(I) = PERM_INV(I) = I
//    before calling this routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint ADJ_NUM, the number of adjacency entries.
//
//    Input, gint ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, gint ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input, gint PERM[NODE_NUM], PERM_INV[NODE_NUM], the permutation
//    and inverse permutation.
//
{
  char *band;
  gint band_lo;
  gint col;
  gint i;
  gint j;
  gint k;
  gint nonzero_num;

  band = g_new0 (gchar, node_num);

  band_lo = 0;
  nonzero_num = 0;

  g_message ("#\n#  Nonzero structure of matrix:\n");

  for ( i = 0; i < node_num; i++ )
  {
    for ( k = 0; k < node_num; k++ )
    {
      band[k] = '.';
    }

    band[i] = 'D';

    for ( j = adj_row[perm[i]-1]; j <= adj_row[perm[i]]-1; j++ )
    {
      col = perm_inv[adj[j-1]-1] - 1;

      if ( col < i )
      {
        nonzero_num = nonzero_num + 1;
      }

      band_lo = i4_max ( band_lo, i - col );

      if ( col != i )
      {
        band[col] = 'X';
      }
    }
    g_message ("#  %4d ", i + 1);
    for ( j = 0; j < node_num; j++ )
    {
      g_message ("%c", band[j]);
    }
    g_message ("\n");
  }

  g_message ("#\n#  Lower bandwidth = %d\n", band_lo);
  g_message ("#  Lower envelope contains %d nonzeros.\n", nonzero_num);

  return;
}
//****************************************************************************80

void adj_print ( gint node_num, gint adj_num, gint adj_row[], gint adj[],
  gchar *title )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_PRINT prints adjacency information.
//
//  Discussion:
//
//    The list has the form:
//
//    Row   Nonzeros
//
//    1       2   5   9
//    2       7   8   9   15   78   79   81  86  91  99
//          100 103
//    3      48  49  53
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 December 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint ADJ_NUM, the number of adjacency entries.
//
//    Input, gint ADJ_ROW[NODE_NUM+1], organizes the adjacency entries
//    into rows.  The entries for row I are in entries ADJ_ROW(I)
//    through ADJ_ROW(I+1)-1.
//
//    Input, gint ADJ[ADJ_NUM], the adjacency structure, which contains,
//    for each row, the column indices of the nonzero entries.
//
//    Input, gchar *TITLE, a title to be printed.
//
{
  adj_print_some ( node_num, 1, node_num, adj_num, adj_row, adj, title );

  return;
}
//****************************************************************************80

void adj_print_some ( gint node_num, gint node_lo, gint node_hi, gint adj_num,
  gint adj_row[], gint adj[], gchar *title )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_PRINT_SOME prints some adjacency information.
//
//  Discussion:
//
//    The list has the form:
//
//    Row   Nonzeros
//
//    1       2   5   9
//    2       7   8   9   15   78   79   81  86  91  99
//          100 103
//    3      48  49  53
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint NODE_LO, NODE_HI, the first and last nodes for
//    which the adjacency information is to be printed.
//
//    Input, gint ADJ_NUM, the number of adjacency entries.
//
//    Input, gint ADJ_ROW[NODE_NUM+1], organizes the adjacency entries
//    into rows.  The entries for row I are in entries ADJ_ROW(I)
//    through ADJ_ROW(I+1)-1.
//
//    Input, gint ADJ[ADJ_NUM], the adjacency structure, which contains,
//    for each row, the column indices of the nonzero entries.
//
//    Input, gchar *TITLE, a title.
//
{
  gint i;
  gint j;
  gint jhi;
  gint jlo;
  gint jmax;
  gint jmin;

  g_message ("#  %s\n#\n"
             "#  Sparse adjacency structure:\n#\n"
             "#    Number of nodes       = %d\n"
             "#    Number of adjacencies = %d\n"
             "#    Node Min   Max      Nonzeros\n", 
             title, node_num, adj_num);

  for ( i = node_lo; i <= node_hi; i++ )
  {
    jmin = adj_row[i-1];
    jmax = adj_row[i] - 1;

    if ( jmax < jmin )
    {
      g_message ("#    %4d %4d %4d |\n", i, jmin, jmax);
    }
    else
    {
      for ( jlo = jmin; jlo <= jmax; jlo = jlo + 5 )
      {
        jhi = i4_min ( jlo + 4, jmax );

        if ( jlo == jmin )
        {
          g_message ("#    %4d %4d %4d |", i, jmin, jmax);
          for ( j = jlo; j <= jhi; j++ )
          {
            g_message (" %4d", adj[j-1]);
          }
          g_message ("\n");
        }
        else
        {
          g_message ("#                   |");
          for ( j = jlo; j <= jhi; j++ )
          {
            g_message (" %4d", adj[j-1]);
          }
          g_message ("\n");
        }
      }
    }
  }

  return;
}
//****************************************************************************80

void adj_set ( gint node_num, gint adj_max, gint *adj_num, gint adj_row[],
  gint adj[], gint irow, gint jcol )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_SET sets up the adjacency information.
//
//  Discussion:
//
//    The routine records the locations of each nonzero element,
//    one at a time.
//
//    The first call for a given problem should be with IROW or ICOL
//    negative.  This is a signal indicating the data structure should
//    be initialized.
//
//    Then, for each case in which A(IROW,JCOL) is nonzero, or
//    in which IROW is adjacent to JCOL, call this routine once
//    to record that fact.
//
//    Diagonal entries are not to be stored.
//
//    The matrix is assumed to be symmetric, so setting I adjacent to J
//    will also set J adjacent to I.
//
//    Repeated calls with the same values of IROW and JCOL do not
//    actually hurt.  No extra storage will be allocated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint ADJ_MAX, the maximum dimension of the adjacency array.
//
//    Input/output, gint *ADJ_NUM, the number of adjaceny entries.
//
//    Input/output, gint ADJ_ROW[NODE_NUM+1].  Information about
//    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input/output, gint ADJ[ADJ_NUM], the adjacency structure.
//
//    Input, gint IROW, JCOL, the row and column indices of a nonzero
//    entry of the matrix.
//
{
  gint i;
//
//  Negative IROW or JCOL indicates the data structure should be initialized.
//
  if ( irow < 0 || jcol < 0 )
  {
    if (FALSE)
    {
      g_message ("#  ADJ_SET - Note:\n"
                 "#    Initializing adjacency information.\n"
                 "#    Number of nodes NODE_NUM =  %d\n"
                 "#    Maximum adjacency ADJ_MAX = %d\n",
                 node_num, adj_max);
    }
    *adj_num = 0;
    for ( i = 0; i < node_num + 1; i++ )
    {
      adj_row[i] = 1;
    }
    for ( i = 0; i < adj_max; i++ )
    {
      adj[i] = 0;
    }
    return;
  }
//
//  Diagonal entries are not stored.
//
  if ( irow == jcol )
  {
    return;
  }

  if ( node_num < irow )
  {
    g_error ("ADJ_SET - Fatal error!\n"
             "  NODE_NUM < IROW.\n"
             "  IROW =     %d\n"
             "  NODE_NUM = %d",
             irow, node_num
             );
  }
  else if ( irow < 1 )
  {
    g_error ("ADJ_SET - Fatal error!\n"
             "  IROW < 1.\n"
             "  IROW = %d",
             irow);
  }
  else if ( node_num < jcol )
  {
    g_error ("ADJ_SET - Fatal error!\n"
             "  NODE_NUM < JCOL.\n"
             "  JCOL =     %d\n"
             "  NODE_NUM = %d",
             jcol, node_num);
  }
  else if ( jcol < 1 )
  {
    g_error ("ADJ_SET - Fatal error!\n"
             "  JCOL < 1.\n"
             "  JCOL = %d",
             jcol);
  }

  if ( !adj_contains_ij ( node_num, *adj_num, adj_row, adj, irow, jcol ) )
  {
    adj_insert_ij ( node_num, adj_max, adj_num, adj_row, adj, irow, jcol );
  }

  if ( !adj_contains_ij ( node_num, *adj_num, adj_row, adj, jcol, irow ) )
  {
    adj_insert_ij ( node_num, adj_max, adj_num, adj_row, adj, jcol, irow );
  }

  return;
}
//****************************************************************************80

void adj_show ( gint node_num, gint adj_num, gint adj_row[], gint adj[] )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_SHOW displays a symbolic picture of an adjacency matrix.
//
//  Discussion:
//
//    The matrix is defined by the adjacency information and a permutation.
//
//    The routine also computes the bandwidth and the size of the envelope.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint ADJ_NUM, the number of adjacency entries.
//
//    Input, gint ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, gint ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
{
  char *band;
  gint band_lo;
  gint col;
  gint i;
  gint j;
  gint k;
  gint nonzero_num;

  band = g_new0 (gchar, node_num);

  band_lo = 0;
  nonzero_num = 0;

  g_message ("#\n#  Nonzero structure of matrix:\n#\n");

  for ( i = 0; i < node_num; i++ )
  {
    for ( k = 0; k < node_num; k++ )
    {
      band[k] = '.';
    }

    band[i] = 'D';

    for ( j = adj_row[i]; j <= adj_row[i+1]-1; j++ )
    {
      col = adj[j-1] - 1;
      if ( col < i )
      {
        nonzero_num = nonzero_num + 1;
      }
      band_lo = MAX ( band_lo, i - col );
      band[col] = 'X';
    }
    g_message ("#  %8d ", i + 1);
    for ( j = 0; j < node_num; j++ )
    {
      g_message ("%c", band[j]);
    }
    g_message ("\n");
  }

  g_message ("#\n#    Lower bandwidth = %d\n"
             "#    Lower envelope contains %d nonzeros.\n",
             band_lo, nonzero_num);

  g_free (band);

  return;
}
//****************************************************************************80

void degree ( gint root, gint adj_num, gint adj_row[], gint adj[], gint mask[],
  gint deg[], gint *iccsze, gint ls[], gint node_num )

//****************************************************************************80
//
//  Purpose:
//
//    DEGREE computes the degrees of the nodes in the connected component.
//
//  Discussion:
//
//    The connected component is specified by MASK and ROOT.
//    Nodes for which MASK is zero are ignored.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, gint ROOT, the node that defines the connected component.
//
//    Input, gint ADJ_NUM, the number of adjacency entries.
//
//    Input, gint ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, gint ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input, gint MASK[NODE_NUM], is nonzero for those nodes which are
//    to be considered.
//
//    Output, gint DEG[NODE_NUM], contains, for each  node in the connected
//    component, its degree.
//
//    Output, gint *ICCSIZE, the number of nodes in the connected component.
//
//    Output, gint LS[NODE_NUM], stores in entries 1 through ICCSIZE the nodes
//    in the connected component, starting with ROOT, and proceeding
//    by levels.
//
//    Input, gint NODE_NUM, the number of nodes.
//
{
  gint i;
  gint ideg;
  gint j;
  gint jstop;
  gint jstrt;
  gint lbegin;
  gint lvlend;
  gint lvsize;
  gint nbr;
  gint node;
//
//  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
//
  ls[0] = root;
  adj_row[root-1] = -adj_row[root-1];
  lvlend = 0;
  *iccsze = 1;
//
//  LBEGIN is the pointer to the beginning of the current level, and
//  LVLEND points to the end of this level.
//
  for ( ; ; )
  {
    lbegin = lvlend + 1;
    lvlend = *iccsze;
//
//  Find the degrees of nodes in the current level,
//  and at the same time, generate the next level.
//
    for ( i = lbegin; i <= lvlend; i++ )
    {
      node = ls[i-1];
      jstrt = -adj_row[node-1];
      jstop = abs ( adj_row[node] ) - 1;
      ideg = 0;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          ideg = ideg + 1;

          if ( 0 <= adj_row[nbr-1] )
          {
            adj_row[nbr-1] = -adj_row[nbr-1];
            *iccsze = *iccsze + 1;
            ls[*iccsze-1] = nbr;
          }
        }
      }
      deg[node-1] = ideg;
    }
//
//  Compute the current level width.
//
    lvsize = *iccsze - lvlend;
//
//  If the current level width is nonzero, generate another level.
//
    if ( lvsize == 0 )
    {
      break;
    }
  }
//
//  Reset ADJ_ROW to its correct sign and return.
//
  for ( i = 0; i < *iccsze; i++ )
  {
    node = ls[i] - 1;
    adj_row[node] = -adj_row[node];
  }

  return;
}
//****************************************************************************80

void genrcm ( gint node_num, gint adj_num, gint adj_row[], gint adj[], gint perm[] )

//****************************************************************************80
//
//  Purpose:
//
//    GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
//
//  Discussion:
//
//    For each connected component in the graph, the routine obtains
//    an ordering by calling RCM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint ADJ_NUM, the number of adjacency entries.
//
//    Input, gint ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, gint  ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Output, gint  PERM[NODE_NUM], the RCM ordering.
//
//  Local Parameters:
//
//    Local, gint  LEVEL_ROW[NODE_NUM+1], the index vector for a level
//    structure.  The level structure is stored in the currently unused
//    spaces in the permutation vector PERM.
//
//    Local, gint MASK[NODE_NUM], marks variables that have been numbered.
//
{
  gint i;
  gint iccsze;
  gint level_num;
  gint *level_row;
  gint *mask;
  gint num;
  gint root;

  level_row = g_new0 (gint, node_num + 1);
  mask      = g_new0 (gint, node_num);

  for ( i = 0; i < node_num; i++ )
  {
    mask[i] = 1;
  }

  num = 1;

  for ( i = 0; i < node_num; i++ )
  {
//
//  For each masked connected component...
//
    if ( mask[i] != 0 )
    {
      root = i + 1;
//
//  Find a pseudo-peripheral node ROOT.  The level structure found by
//  ROOT_FIND is stored starting at PERM(NUM).
//
      root_find ( &root, adj_num, adj_row, adj, mask, &level_num,
        level_row, perm+num-1, node_num );
//
//  RCM orders the component using ROOT as the starting node.
//
      rcm ( root, adj_num, adj_row, adj, mask, perm+num-1, &iccsze,
        node_num );

      num = num + iccsze;
//
//  We can stop once every node is in one of the connected components.
//
      if ( node_num < num )
      {
        g_free (level_row);
        g_free (mask);
        return;
      }
    }
  }

  g_free (level_row);
  g_free (mask);
  
  return;
}
//****************************************************************************80

void graph_01_adj ( gint node_num, gint adj_num, gint adj_row[], gint adj[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRAPH_01_ADJ returns the adjacency vector for graph 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint ADJ_NUM, the number of adjacencies.
//
//    Output, gint ADJ_ROW[NODE_NUM+1], node pointers into ADJ.
//
//    Output, gint ADJ[ADJ_NUM], the adjacency information.
//
{
# define ADJ_NUM 28
# define NODE_NUM 10

  static gint adj_save[ADJ_NUM] = {
    4, 6,
    3, 5, 7, 10,
    2, 4, 5,
    1, 3, 6, 9,
    2, 3, 7,
    1, 4, 7, 8,
    2, 5, 6, 8,
    6, 7,
    4,
    2 };
  static gint adj_row_save[NODE_NUM+1] = {
    1, 3, 7, 10, 14, 17, 21, 25, 27, 28, 29
  };
  gint i;

  for ( i = 0; i < ADJ_NUM; i++ )
  {
    adj[i] = adj_save[i];
  }

  for ( i = 0; i < NODE_NUM + 1; i++ )
  {
    adj_row[i] = adj_row_save[i];
  }
  return;
# undef ADJ_NUM
# undef NODE_NUM
}
//****************************************************************************80

void graph_01_size ( gint *node_num, gint *adj_num )

//****************************************************************************80
//
//  Purpose:
//
//    GRAPH_01_SIZE returns the number of adjacencies for graph 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Output, gint *NODE_NUM, the number of items that can be adjacent.
//
//    Output, gint *ADJ_NUM, the number of adjacencies.
//
{
  *node_num = 10;
  *adj_num = 28;

  return;
}
//****************************************************************************80

gint i4_max ( gint i1, gint i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint I1, I2, are two integers to be compared.
//
//    Output, gint I4_MAX, the larger of I1 and I2.
//
{
  if ( i2 < i1 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

gint i4_min ( gint i1, gint i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint I1, I2, two integers to be compared.
//
//    Output, gint I4_MIN, the smaller of I1 and I2.
//
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

gint i4_sign ( gint i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SIGN returns the sign of an I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint I, the integer whose sign is desired.
//
//    Output, gint I4_SIGN, the sign of I.
{
  gint value;

  if ( i < 0 )
  {
    value = -1;
  }
  else
  {
    value = 1;
  }
  return value;
}
//****************************************************************************80

void i4_swap ( gint *i, gint *j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SWAP switches two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, gint *I, *J.  On output, the values of I and
//    J have been interchanged.
//
{
  gint k;

  k = *i;
  *i = *j;
  *j = k;

  return;
}
//****************************************************************************80

gint i4_uniform ( gint a, gint b, gint *seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM returns a scaled pseudorandom I4.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, gint A, B, the limits of the interval.
//
//    Input/output, gint *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, gint I4_UNIFORM, a number between A and B.
//
{
  gint k;
  float r;
  gint value;

  if ( *seed == 0 )
  {
    g_error ("I4_UNIFORM - Fatal error!\n"
             "  Input value of SEED = 0.");
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 )
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}
//****************************************************************************80

gint i4col_compare ( gint m, gint n, gint a[], gint i, gint j )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_COMPARE compares columns I and J of an I4COL.
//
//  Example:
//
//    Input:
//
//      M = 3, N = 4, I = 2, J = 4
//
//      A = (
//        1  2  3  4
//        5  6  7  8
//        9 10 11 12 )
//
//    Output:
//
//      I4COL_COMPARE = -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint M, N, the number of rows and columns.
//
//    Input, gint A[M*N], an array of N columns of vectors of length M.
//
//    Input, gint I, J, the columns to be compared.
//    I and J must be between 1 and N.
//
//    Output, gint I4COL_COMPARE, the results of the comparison:
//    -1, column I < column J,
//     0, column I = column J,
//    +1, column J < column I.
//
{
  gint k;
//
//  Check.
//
  if ( i < 1 )
  {
    g_error ("I4COL_COMPARE - Fatal error!\n"
             "  Column index I = %d is less than 1.",
             i);
  }

  if ( n < i )
  {
    g_error ("I4COL_COMPARE - Fatal error!\n"
             "  N = %d is less than column index I = %d.",
             n, i);
  }

  if ( j < 1 )
  {
    g_error ("I4COL_COMPARE - Fatal error!\n"
             "  Column index J = %d is less than 1.",
             j);
  }

  if ( n < j )
  {
    g_error ("I4COL_COMPARE - Fatal error!\n"
             "  N = %d is less than column index J = %d.",
             n, j);
  }

  if ( i == j )
  {
    return 0;
  }

  k = 1;

  while ( k <= m )
  {
    if ( a[k-1+(i-1)*m] < a[k-1+(j-1)*m] )
    {
      return (-1);
    }
    else if ( a[k-1+(j-1)*m] < a[k-1+(i-1)*m] )
    {
      return 1;
    }
    k = k + 1;
  }

  return 0;
}
//****************************************************************************80

void i4col_sort_a ( gint m, gint n, gint a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORT_A ascending sorts the columns of an I4COL.
//
//  Discussion:
//
//    In lexicographic order, the statement "X < Y", applied to two
//    vectors X and Y of length M, means that there is some index I, with
//    1 <= I <= M, with the property that
//
//      X(J) = Y(J) for J < I,
//    and
//      X(I) < Y(I).
//
//    In other words, X is less than Y if, at the first index where they
//    differ, the X value is less than the Y value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint M, the number of rows of A.
//
//    Input, gint N, the number of columns of A.
//
//    Input/output, gint A[M*N].
//    On input, the array of N columns of M vectors;
//    On output, the columns of A have been sorted in ascending
//    lexicographic order.
//
{
  gint i;
  gint indx;
  gint isgn;
  gint j;
//
//  Initialize.
//
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
//
//  Call the external heap sorter.
//
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
    if ( 0 < indx )
    {
      i4col_swap ( m, n, a, i, j );
    }
//
//  Compare the I and J objects.
//
    else if ( indx < 0 )
    {
      isgn = i4col_compare ( m, n, a, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }

  }

  return;
}
//****************************************************************************80

void i4col_swap ( gint m, gint n, gint a[], gint icol1, gint icol2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SWAP swaps two columns of an I4COL.
//
//  Discussion:
//
//    The two dimensional information is stored as a one dimensional
//    array, by columns.
//
//    The row indices are 1 based, NOT 0 based//  However, a preprocessor
//    variable, called OFFSET, can be reset from 1 to 0 if you wish to
//    use 0-based indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint M, N, the number of rows and columns.
//
//    Input/output, gint A[M*N], an array of data.
//
//    Input, gint ICOL1, ICOL2, the two columns to swap.
//    These indices should be between 1 and N.
//
{
# define OFFSET 1

  gint i;
  gint t;
//
//  Check.
//
  if ( icol1 - OFFSET < 0 || n-1 < icol1 - OFFSET )
  {
    g_error ("I4COL_SWAP - Fatal error!\n"
             "  ICOL1 is out of range.");
  }

  if ( icol2 - OFFSET < 0 || n-1 < icol2 - OFFSET )
  {
    g_error ("I4COL_SWAP - Fatal error!\n"
             "  ICOL2 is out of range.");
  }

  if ( icol1 == icol2 )
  {
    return;
  }
  for ( i = 0; i < m; i++ )
  {
    t                     = a[i+(icol1-OFFSET)*m];
    a[i+(icol1-OFFSET)*m] = a[i+(icol2-OFFSET)*m];
    a[i+(icol2-OFFSET)*m] = t;
  }

  return;
# undef OFFSET
}
//****************************************************************************80

void i4mat_print_some ( gint m, gint n, gint a[], gint ilo, gint jlo, gint ihi,
  gint jhi, gchar *title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT_SOME prints some of an I4MAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, gint N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, gint A[M*N], the matrix.
//
//    Input, gint ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, gchar *TITLE, a title.
{
# define INCX 10

  gint i;
  gint i2hi;
  gint i2lo;
  gint j;
  gint j2hi;
  gint j2lo;

  g_message ("#  %s", title);
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    g_message ("\n");
//
//  For each column J in the current range...
//
//  Write the header.
//
    g_message ("#    Col: ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      g_message ("%6d  ", j);
    }
    g_message ("\n#    Row\n#\n");
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to INCX) entries in row I, that lie in the current strip.
//
      g_message ("#  %5d  ", i);
      for ( j = j2lo; j <= j2hi; j++ )
      {
        g_message ("%6d  ", a[i-1+(j-1)*m]);
      }
      g_message ("\n");
    }
  }

  g_message ("#\n");

  return;
# undef INCX
}
//****************************************************************************80

void i4mat_transpose_print ( gint m, gint n, gint a[], gchar *title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint M, the number of rows in A.
//
//    Input, gint N, the number of columns in A.
//
//    Input, gint A[M*N], the M by N matrix.
//
//    Input, gchar *TITLE, a title.
//
{
  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );
  return;
}
//****************************************************************************80

void i4mat_transpose_print_some ( gint m, gint n, gint a[], gint ilo, gint jlo,
  gint ihi, gint jhi, gchar *title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, gint N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, gint A[M*N], the matrix.
//
//    Input, gint ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, gchar *TITLE, a title.
{
# define INCX 10

  gint i;
  gint i2hi;
  gint i2lo;
  gint j;
  gint j2hi;
  gint j2lo;

  g_message ("#  %s", title);
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    g_message ("\n");
//
//  For each row I in the current range...
//
//  Write the header.
//
    g_message ("#    Row: ");
    for ( i = i2lo; i <= i2hi; i++ )
    {
      g_message ("#  %6d  ", i);
    }
    g_message ("\n#    Row\n#\n");
//
//  Determine the range of the rows in this strip.
//
    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
//
//  Print out (up to INCX) entries in column J, that lie in the current strip.
//
      g_message ("#  %5d  ", j);

      for ( i = i2lo; i <= i2hi; i++ )
      {
        g_message ("%6d  ", a[i-1+(j-1)*m]);
      }
      g_message ("\n");
    }

  }

  g_message ("#\n");

  return;
# undef INCX
}
//****************************************************************************80

void i4vec_heap_d ( gint n, gint a[] )

/****************************************************************************80
//
//  Purpose:
//
//    I4VEC_HEAP_D reorders an I4VEC into a descending heap.
//
//  Discussion:
//
//    A heap is an array A with the property that, for every index J,
//    A[J] >= A[2*J+1] and A[J] >= A[2*J+2], (as long as the indices
//    2*J+1 and 2*J+2 are legal).
//
//  Diagram:
//
//                  A(0)
//                /      \
//            A(1)         A(2)
//          /     \        /  \
//      A(3)       A(4)  A(5) A(6)
//      /  \       /   \
//    A(7) A(8)  A(9) A(10)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, gint N, the size of the input array.
//
//    Input/output, gint A[N].
//    On input, an unsorted array.
//    On output, the array has been reordered into a heap.
//*/
{
  gint i;
  gint ifree;
  gint key;
  gint m;
//
//  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
//
  for ( i = (n/2)-1; 0 <= i; i-- )
  {
//
//  Copy the value out of the parent node.
//  Position IFREE is now "open".
//
    key = a[i];
    ifree = i;

    for ( ;; )
    {
//
//  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
//  IFREE.  (One or both may not exist because they equal or exceed N.)
//
      m = 2 * ifree + 1;
//
//  Does the first position exist?
//
      if ( n <= m )
      {
        break;
      }
      else
      {
//
//  Does the second position exist?
//
        if ( m + 1 < n )
        {
//
//  If both positions exist, take the larger of the two values,
//  and update M if necessary.
//
          if ( a[m] < a[m+1] )
          {
            m = m + 1;
          }
        }
//
//  If the large descendant is larger than KEY, move it up,
//  and update IFREE, the location of the free position, and
//  consider the descendants of THIS position.
//
        if ( key < a[m] )
        {
          a[ifree] = a[m];
          ifree = m;
        }
        else
        {
          break;
        }

      }

    }
//
//  When you have stopped shifting items up, return the item you
//  pulled out back to the heap.
//
    a[ifree] = key;

  }

  return;
}
//****************************************************************************80

gint *i4vec_indicator ( gint n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_INDICATOR sets an I4VEC to the indicator vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint N, the number of elements of A.
//
//    Output, gint I4VEC_INDICATOR(N), the initialized array.
//
{
  gint *a;
  gint i;

  a = g_new0 (gint, n);

  for ( i = 0; i < n; i++ )
  {
    a[i] = i + 1;
  }

  return a;
}
//****************************************************************************80

void i4vec_print ( gint n, gint a[], gchar *title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint N, the number of components of the vector.
//
//    Input, gint A[N], the vector to be printed.
//
//    Input, gchar *TITLE, a title.
//
{
  gint i;

  g_message ("#  %s\n#\n#  ", title);
  for ( i = 0; i <= n-1; i++ )
  {
    g_message ("  %8d  %8d", i + 1, a[i]);
  }

  return;
}
//****************************************************************************80

void i4vec_reverse ( gint n, gint a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_REVERSE reverses the elements of an I4VEC.
//
//  Example:
//
//    Input:
//
//      N = 5,
//      A = ( 11, 12, 13, 14, 15 ).
//
//    Output:
//
//      A = ( 15, 14, 13, 12, 11 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint N, the number of entries in the array.
//
//    Input/output, gint A(N), the array to be reversed.
//
{
  gint i;
  gint j;

  for ( i = 0; i < n / 2; i++ )
  {
    j        = a[i];
    a[i]     = a[n-1-i];
    a[n-1-i] = j;
  }

  return;
}
//****************************************************************************80

void i4vec_sort_heap_a ( gint n, gint a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, gint N, the number of entries in the array.
//
//    Input/output, gint A[N].
//    On input, the array to be sorted;
//    On output, the array has been sorted.
//
{
  gint n1;
  gint temp;

  if ( n <= 1 )
  {
    return;
  }
//
//  1: Put A into descending heap form.
//
  i4vec_heap_d ( n, a );
//
//  2: Sort A.
//
//  The largest object in the heap is in A[0].
//  Move it to position A[N-1].
//
  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;
//
//  Consider the diminished heap of size N1.
//
  for ( n1 = n-1; 2 <= n1; n1-- )
  {
//
//  Restore the heap structure of the initial N1 entries of A.
//
    i4vec_heap_d ( n1, a );
//
//  Take the largest object from A[0] and move it to A[N1-1].
//
    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;
  }

  return;
}
//****************************************************************************80

void level_set ( gint root, gint adj_num, gint adj_row[], gint adj[], gint mask[],
  gint *level_num, gint level_row[], gint level[], gint node_num )

//****************************************************************************80
//
//  Purpose:
//
//    LEVEL_SET generates the connected level structure rooted at a given node.
//
//  Discussion:
//
//    Only nodes for which MASK is nonzero will be considered.
//
//    The root node chosen by the user is assigned level 1, and masked.
//    All (unmasked) nodes reachable from a node in level 1 are
//    assigned level 2 and masked.  The process continues until there
//    are no unmasked nodes adjacent to any node in the current level.
//    The number of levels may vary between 2 and NODE_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, gint ROOT, the node at which the level structure
//    is to be rooted.
//
//    Input, gint ADJ_NUM, the number of adjacency entries.
//
//    Input, gint ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, gint ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input/output, gint MASK[NODE_NUM].  On input, only nodes with nonzero
//    MASK are to be processed.  On output, those nodes which were included
//    in the level set have MASK set to 1.
//
//    Output, gint *LEVEL_NUM, the number of levels in the level
//    structure.  ROOT is in level 1.  The neighbors of ROOT
//    are in level 2, and so on.
//
//    Output, gint LEVEL_ROW[NODE_NUM+1], LEVEL[NODE_NUM], the rooted
//    level structure.
//
//    Input, gint NODE_NUM, the number of nodes.
//
{
  gint i;
  gint iccsze;
  gint j;
  gint jstop;
  gint jstrt;
  gint lbegin;
  gint lvlend;
  gint lvsize;
  gint nbr;
  gint node;

  mask[root-1] = 0;
  level[0] = root;
  *level_num = 0;
  lvlend = 0;
  iccsze = 1;
//
//  LBEGIN is the pointer to the beginning of the current level, and
//  LVLEND points to the end of this level.
//
  for ( ; ; )
  {
    lbegin = lvlend + 1;
    lvlend = iccsze;
    *level_num = *level_num + 1;
    level_row[*level_num-1] = lbegin;
//
//  Generate the next level by finding all the masked neighbors of nodes
//  in the current level.
//
    for ( i = lbegin; i <= lvlend; i++ )
    {
      node = level[i-1];
      jstrt = adj_row[node-1];
      jstop = adj_row[node] - 1;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          iccsze = iccsze + 1;
          level[iccsze-1] = nbr;
          mask[nbr-1] = 0;
        }
      }
    }
//
//  Compute the current level width (the number of nodes encountered.)
//  If it is positive, generate the next level.
//
    lvsize = iccsze - lvlend;

    if ( lvsize <= 0 )
    {
      break;
    }
  }

  level_row[*level_num] = lvlend + 1;
//
//  Reset MASK to 1 for the nodes in the level structure.
//
  for ( i = 0; i < iccsze; i++ )
  {
    mask[level[i]-1] = 1;
  }

  return;
}
//****************************************************************************80

void level_set_print ( gint node_num, gint level_num, gint level_row[],
  gint level[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEVEL_SET_PRINT prints level set information.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint LEVEL_NUM, the number of levels.
//
//    Input, gint LEVEL_ROW[LEVEL_NUM+1], organizes the entries of LEVEL.
//    The entries for level I are in entries LEVEL_ROW(I)
//    through LEVEL_ROW(I+1)-1.
//
//    Input, integer LEVEL[NODE_NUM], is simply a list of the nodes in an
//    order induced by the levels.
//
{
  gint i;
  gint j;
  gint jhi;
  gint jlo;
  gint jmax;
  gint jmin;

  g_message ("#  LEVEL_SET_PRINT\n"
             "#    Show the level set structure of a rooted graph.\n"
             "#    The number of nodes is  %d\n"
             "#    The number of levels is %d\n"
             "#\n"
             "#    Level Min Max      Nonzeros\n"
             "#\n",
             node_num, level_num
             );

  for ( i = 0; i < level_num; i++ )
  {
    jmin = level_row[i];
    jmax = level_row[i+1] - 1;

    if ( jmax < jmin )
    {
      g_message ("#    %4d  %4d  %4d", i + 1, jmin, jmax);
    }
    else
    {
      for ( jlo = jmin; jlo <= jmax; jlo = jlo + 5 )
      {
        jhi = i4_min ( jlo + 4, jmax );

        if ( jlo == jmin )
        {
          g_message ("#    %4d  %4d  %4d   ", i + 1, jmin, jmax);
          for ( j = jlo; j <= jhi; j++ )
          {
            g_message ("%8d", level[j-1]);
          }
          g_message ("\n");
        }
        else
        {
          g_message ("#                       ");
          for ( j = jlo; j <= jhi; j++ )
          {
            g_message ("%8d", level[j-1]);
          }
          g_message ("\n");
        }
      }
    }
  }

  return;
}
//****************************************************************************80

gboolean perm_check ( gint n, gint p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CHECK checks that a vector represents a permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from 1
//    to N occurs among the N entries of the permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint N, the number of entries.
//
//    Input, gint P[N], the array to check.
//
//    Output, gboolean PERM_CHECK, is TRUE if the permutation is OK.
//
{
  gboolean found;
  gint i;
  gint seek;

  for ( seek = 1; seek <= n; seek++ )
  {
    found = FALSE;

    for ( i = 0; i < n; i++ )
    {
      if ( p[i] == seek )
      {
        found = TRUE;
        break;
      }
    }

    if ( !found )
    {
      return FALSE;
    }
  }

  return TRUE;
}
//****************************************************************************80

void perm_inverse3 ( gint n, gint perm[], gint perm_inv[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INVERSE3 produces the inverse of a given permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint N, the number of items permuted.
//
//    Input, gint PERM[N], a permutation.
//
//    Output, gint PERM_INV[N], the inverse permutation.
//
{
  gint i;

  for ( i = 0; i < n; i++ )
  {
    perm_inv[perm[i]-1] = i + 1;
  }

  return;
}
//****************************************************************************80

gint *perm_uniform ( gint n, gint *seed )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_UNIFORM selects a random permutation of N objects.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, gint N, the number of objects to be permuted.
//
//    Input/output, gint *SEED, a seed for the random number generator.
//
//    Output, gint PERM_UNIFORM[N], a permutation of (1,, 1, ..., N).
//
{
  gint i;
  gint j;
  gint *p;

  p = g_new0 (gint, n);

  for ( i = 1; i <= n; i++ )
  {
    p[i-1] = i;
  }

  for ( i = 1; i <= n; i++ )
  {
    j = i4_uniform ( i, n, seed );
    i4_swap ( &p[i-1], &p[j-1] );
  }

  return p;
}
//****************************************************************************80

float r4_abs ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ABS returns the absolute value of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the quantity whose absolute value is desired.
//
//    Output, float R4_ABS, the absolute value of X.
//
{
  float value;

  if ( 0.0 <= x )
  {
    value = x;
  }
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

gint r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Examples:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the value.
//
//    Output, gint R4_NINT, the nearest integer to X.
//
{
  gint value;

  if ( x < 0.0 )
  {
    value = - ( gint ) ( r4_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( gint ) ( r4_abs ( x ) + 0.5 );
  }

  return value;
}
//****************************************************************************80

void r82vec_permute ( gint n, gdouble a[], gint p[] )

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_PERMUTE permutes an R82VEC in place.
//
//  Discussion:
//
//    An R82VEC is a vector whose entries are R82's.
//    An R82 is a vector of type gdouble precision with two entries.
//    An R82VEC may be stored as a 2 by N array.
//
//    This routine permutes an array of real "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5
//      P = (   2,    4,    5,    1,    3 )
//      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
//          (11.0, 22.0, 33.0, 44.0, 55.0 )
//
//    Output:
//
//      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
//             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint N, the number of objects.
//
//    Input/output, gdouble A[2*N], the array to be permuted.
//
//    Input, gint P[N], the permutation.  P(I) = J means
//    that the I-th element of the output array should be the J-th
//    element of the input array.  P must be a legal permutation
//    of the integers from 1 to N, otherwise the algorithm will
//    fail catastrophically.
//
{
  gdouble a_temp[2];
  gint i;
  gint iget;
  gint iput;
  gint istart;

  if ( !perm_check ( n, p ) )
  {
    g_error ("R82VEC_PERMUTE - Fatal error!\n"
             "  The input array does not represent\n"
             "  a proper permutation."
             );
    i4vec_print (n, p, "  The faulty permutation:");
  }
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = -p[istart-1];
      continue;
    }
    else
    {
      a_temp[0] = a[0+(istart-1)*2];
      a_temp[1] = a[1+(istart-1)*2];
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = -p[iput-1];

        if ( iget < 1 || n < iget )
        {
          g_error ("R82VEC_PERMUTE - Fatal error!\n"
                   "  Entry IPUT = %d of the permutation has\n"
                   "  an illegal value IGET = %d.\n",
                   iput, iget);
        }

        if ( iget == istart )
        {
          a[0+(iput-1)*2] = a_temp[0];
          a[1+(iput-1)*2] = a_temp[1];
          break;
        }
        a[0+(iput-1)*2] = a[0+(iget-1)*2];
        a[1+(iput-1)*2] = a[1+(iget-1)*2];
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = -p[i];
  }

  return;
}
//****************************************************************************80

void r8mat_print_some ( gint m, gint n, gdouble a[], gint ilo, gint jlo, gint ihi,
  gint jhi, gchar *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, gint N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, gdouble A[M*N], the matrix.
//
//    Input, gint ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, gchar *TITLE, a title.
{
# define INCX 5

  gint i;
  gint i2hi;
  gint i2lo;
  gint j;
  gint j2hi;
  gint j2lo;

  g_message ("#  %s\n", title);
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    g_message ("#  \n");
//
//  For each column J in the current range...
//
//  Write the header.
//
    g_message ("#    Col:    ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      g_message ("%7d       ", j);
    }
    g_message ("\n#    Row\n#\n");
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      g_message ("#  %5d  ", i);
      for ( j = j2lo; j <= j2hi; j++ )
      {
        g_message ("#  %12g  ", a[i-1+(j-1)*m]);
      }
      g_message ("\n");
    }

  }

  g_message ("#\n");

  return;
# undef INCX
}
//****************************************************************************80

void r8mat_transpose_print_some ( gint m, gint n, gdouble a[], gint ilo, gint jlo,
  gint ihi, gint jhi, gchar *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint M, N, the number of rows and columns.
//
//    Input, gdouble A[M*N], an M by N matrix to be printed.
//
//    Input, gint ILO, JLO, the first row and column to print.
//
//    Input, gint IHI, JHI, the last row and column to print.
//
//    Input, gchar *TITLE, a title.
//
{
# define INCX 5

  gint i;
  gint i2;
  gint i2hi;
  gint i2lo;
  gint inc;
  gint j;
  gint j2hi;
  gint j2lo;

  g_message ("#  %s\n", title);

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    g_message ("#\n#    Row: ");
    for ( i = i2lo; i <= i2hi; i++ )
    {
      g_message ("%7d       ", i);
    }
    g_message ("\n#    Col\n#\n");

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      g_message ("#  %5d ", j);
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        g_message ("#  %14g  ", a[(i-1)+(j-1)*m]);
      }
      g_message ("\n");
    }
  }
  g_message ("#\n");

  return;
# undef INCX
}
//****************************************************************************80

void rcm ( gint root, gint adj_num, gint adj_row[], gint adj[], gint mask[],
  gint perm[], gint *iccsze, gint node_num )

//****************************************************************************80
//
//  Purpose:
//
//    RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
//
//  Discussion:
//
//    The connected component is specified by a node ROOT and a mask.
//    The numbering starts at the root node.
//
//    An outline of the algorithm is as follows:
//
//    X(1) = ROOT.
//
//    for ( I = 1 to N-1)
//      Find all unlabeled neighbors of X(I),
//      assign them the next available labels, in order of increasing degree.
//
//    When done, reverse the ordering.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 August 2013
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, gint ROOT, the node that defines the connected component.
//    It is used as the starting point for the RCM ordering.
//
//    Input, gint ADJ_NUM, the number of adjacency entries.
//
//    Input, gint ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, gint ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input/output, gint MASK[NODE_NUM], a mask for the nodes.  Only
//    those nodes with nonzero input mask values are considered by the
//    routine.  The nodes numbered by RCM will have their mask values
//    set to zero.
//
//    Output, gint PERM[NODE_NUM], the RCM ordering.
//
//    Output, gint *ICCSZE, the size of the connected component
//    that has been numbered.
//
//    Input, gint NODE_NUM, the number of nodes.
//
//  Local Parameters:
//
//    Workspace, gint DEG[NODE_NUM], a temporary vector used to hold
//    the degree of the nodes in the section graph specified by mask and root.
//
{
  gint *deg;
  gint fnbr;
  gint i;
  gint j;
  gint jstop;
  gint jstrt;
  gint k;
  gint l;
  gint lbegin;
  gint lnbr;
  gint lperm;
  gint lvlend;
  gint nbr;
  gint node;
//
//  If node_num out of bounds, something is wrong.
//
  if ( node_num < 1 )
  {
    g_error ("RCM - Fatal error!\n"
             "  Unacceptable input value of NODE_NUM = %d\n",
             node_num);
  }
//
//  If the root is out of bounds, something is wrong.
//
  if ( root < 1 || node_num < root )
  {
    g_error ("RCM - Fatal error!\n"
             "  Unacceptable input value of ROOT = %d\n"
             "  Acceptable values are between 1 and %d, inclusive.\n",
             root, node_num);
  }
//
//  Allocate memory for the degree array.
//
  deg = g_new0 (gint, node_num);
//
//  Find the degrees of the nodes in the component specified by MASK and ROOT.
//
  degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, perm, node_num );
//
//  If the connected component size is less than 1, something is wrong.
//
  if ( *iccsze < 1 )
  {
    g_error ("RCM - Fatal error!\n"
             "  Connected component size ICCSZE returned from DEGREE as %d",
             *iccsze);
  }
//
//  Set the mask value for the root.
//
  mask[root-1] = 0;
//
//  If the connected component is a singleton, there is no ordering necessary.
//
  if ( *iccsze == 1 )
  {
    g_free (deg);
    return;
  }
//
//  Carry out the reordering.
//
//  LBEGIN and LVLEND point to the beginning and
//  the end of the current level respectively.
//
  lvlend = 0;
  lnbr = 1;

  while ( lvlend < lnbr )
  {
    lbegin = lvlend + 1;
    lvlend = lnbr;

    for ( i = lbegin; i <= lvlend; i++ )
    {
//
//  For each node in the current level...
//
      node = perm[i-1];
      jstrt = adj_row[node-1];
      jstop = adj_row[node] - 1;
//
//  Find the unnumbered neighbors of NODE.
//
//  FNBR and LNBR point to the first and last neighbors
//  of the current node in PERM.
//
      fnbr = lnbr + 1;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          lnbr = lnbr + 1;
          mask[nbr-1] = 0;
          perm[lnbr-1] = nbr;
        }
      }
//
//  If no neighbors, skip to next node in this level.
//
      if ( lnbr <= fnbr )
      {
        continue;
      }
//
//  Sort the neighbors of NODE in increasing order by degree.
//  Linear insertion is used.
//
      k = fnbr;

      while ( k < lnbr )
      {
        l = k;
        k = k + 1;
        nbr = perm[k-1];

        while ( fnbr < l )
        {
          lperm = perm[l-1];

          if ( deg[lperm-1] <= deg[nbr-1] )
          {
            break;
          }

          perm[l] = lperm;
          l = l - 1;
        }
        perm[l] = nbr;
      }
    }
  }
//
//  We now have the Cuthill-McKee ordering.  
//  Reverse it to get the Reverse Cuthill-McKee ordering.
//
  i4vec_reverse ( *iccsze, perm );
//
//  Free memory.
//
  g_free (deg);

  return;
}
//****************************************************************************80

void root_find ( gint *root, gint adj_num, gint adj_row[], gint adj[], gint mask[],
  gint *level_num, gint level_row[], gint level[], gint node_num )

//****************************************************************************80
//
//  Purpose:
//
//    ROOT_FIND finds a pseudo-peripheral node.
//
//  Discussion:
//
//    The diameter of a graph is the maximum distance (number of edges)
//    between any two nodes of the graph.
//
//    The eccentricity of a node is the maximum distance between that
//    node and any other node of the graph.
//
//    A peripheral node is a node whose eccentricity equals the
//    diameter of the graph.
//
//    A pseudo-peripheral node is an approximation to a peripheral node;
//    it may be a peripheral node, but all we know is that we tried our
//    best.
//
//    The routine is given a graph, and seeks pseudo-peripheral nodes,
//    using a modified version of the scheme of Gibbs, Poole and
//    Stockmeyer.  It determines such a node for the section subgraph
//    specified by MASK and ROOT.
//
//    The routine also determines the level structure associated with
//    the given pseudo-peripheral node; that is, how far each node
//    is from the pseudo-peripheral node.  The level structure is
//    returned as a list of nodes LS, and pointers to the beginning
//    of the list of nodes that are at a distance of 0, 1, 2, ...,
//    NODE_NUM-1 from the pseudo-peripheral node.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//    Norman Gibbs, William Poole, Paul Stockmeyer,
//    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
//    SIAM Journal on Numerical Analysis,
//    Volume 13, pages 236-250, 1976.
//
//    Norman Gibbs,
//    Algorithm 509: A Hybrid Profile Reduction Algorithm,
//    ACM Transactions on Mathematical Software,
//    Volume 2, pages 378-387, 1976.
//
//  Parameters:
//
//    Input/output, gint *ROOT.  On input, ROOT is a node in the
//    the component of the graph for which a pseudo-peripheral node is
//    sought.  On output, ROOT is the pseudo-peripheral node obtained.
//
//    Input, gint ADJ_NUM, the number of adjacency entries.
//
//    Input, gint ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, gint ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input, gint MASK[NODE_NUM], specifies a section subgraph.  Nodes
//    for which MASK is zero are ignored by FNROOT.
//
//    Output, gint *LEVEL_NUM, is the number of levels in the level structure
//    rooted at the node ROOT.
//
//    Output, gint LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the
//    level structure array pair containing the level structure found.
//
//    Input, gint NODE_NUM, the number of nodes.
//
{
  gint iccsze;
  gint j;
  gint jstrt;
  gint k;
  gint kstop;
  gint kstrt;
  gint level_num2;
  gint mindeg;
  gint nabor;
  gint ndeg;
  gint node;
//
//  Determine the level structure rooted at ROOT.
//
  level_set ( *root, adj_num, adj_row, adj, mask, level_num,
    level_row, level, node_num );
//
//  Count the number of nodes in this level structure.
//
  iccsze = level_row[*level_num] - 1;
//
//  Extreme case:
//    A complete graph has a level set of only a single level.
//    Every node is equally good (or bad).
//
  if ( *level_num == 1 )
  {
    return;
  }
//
//  Extreme case:
//    A "line graph" 0--0--0--0--0 has every node in its only level.
//    By chance, we've stumbled on the ideal root.
//
  if ( *level_num == iccsze )
  {
    return;
  }
//
//  Pick any node from the last level that has minimum degree
//  as the starting point to generate a new level set.
//
  for ( ; ; )
  {
    mindeg = iccsze;

    jstrt = level_row[*level_num-1];
    *root = level[jstrt-1];

    if ( jstrt < iccsze )
    {
      for ( j = jstrt; j <= iccsze; j++ )
      {
        node = level[j-1];
        ndeg = 0;
        kstrt = adj_row[node-1];
        kstop = adj_row[node] - 1;

        for ( k = kstrt; k <= kstop; k++ )
        {
          nabor = adj[k-1];
          if ( 0 < mask[nabor-1] )
          {
            ndeg = ndeg + 1;
          }
        }

        if ( ndeg < mindeg )
        {
          *root = node;
          mindeg = ndeg;
        }
      }
    }
//
//  Generate the rooted level structure associated with this node.
//
    level_set ( *root, adj_num, adj_row, adj, mask, &level_num2,
      level_row, level, node_num );
//
//  If the number of levels did not increase, accept the new ROOT.
//
    if ( level_num2 <= *level_num )
    {
      break;
    }

    *level_num = level_num2;
//
//  In the unlikely case that ROOT is one endpoint of a line graph,
//  we can exit now.
//
    if ( iccsze <= *level_num )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

void sort_heap_external ( gint n, gint *indx, gint *i, gint *j, gint isgn )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list is not passed to the routine.  Hence it may
//    consist of integers, reals, numbers, names, etc.  The user,
//    after each return from the routine, will be asked to compare or
//    interchange two items.
//
//    The current version of this code mimics the FORTRAN version,
//    so the values of I and J, in particular, are FORTRAN indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, gint N, the length of the input list.
//
//    Input/output, gint *INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, gint *I, *J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, gint ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
  static gint i_save = 0;
  static gint j_save = 0;
  static gint k = 0;
  static gint k1 = 0;
  static gint n1 = 0;
//
//  INDX = 0: This is the first call.
//
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( *indx < 0 )
  {
    if ( *indx == -2 )
    {
      if ( isgn < 0 )
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn )
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      if ( n1 == 1 )
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }

    k = k - 1;
    k1 = k;

  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( *indx == 1 )
  {
    k1 = k;
  }

  for ( ;; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 )
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 )
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 )
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
  }

  return;
}
//****************************************************************************80

void timestamp ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  /*size_t len;*/
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  /*len =*/ strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  g_message ("#  %s\n", time_buffer);

  return;
# undef TIME_SIZE
}
//****************************************************************************80

gint *triangulation_neighbor_triangles ( gint triangle_order, gint triangle_num,
  gint triangle_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_NEIGHBOR_TRIANGLES determines triangle neighbors.
//
//  Discussion:
//
//    A triangulation of a set of nodes can be completely described by
//    the coordinates of the nodes, and the list of nodes that make up
//    each triangle.  However, in some cases, it is necessary to know
//    triangle adjacency information, that is, which triangle, if any,
//    is adjacent to a given triangle on a particular side.
//
//    This routine creates a data structure recording this information.
//
//    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
//    data items.
//
//    This routine was modified to work with columns rather than rows.
//
//  Example:
//
//    The input information from TRIANGLE_NODE:
//
//    Triangle   Nodes
//    --------   ---------------
//     1         3      4      1
//     2         3      1      2
//     3         3      2      8
//     4         2      1      5
//     5         8      2     13
//     6         8     13      9
//     7         3      8      9
//     8        13      2      5
//     9         9     13      7
//    10         7     13      5
//    11         6      7      5
//    12         9      7      6
//    13        10      9      6
//    14         6      5     12
//    15        11      6     12
//    16        10      6     11
//
//    The output information in TRIANGLE_NEIGHBOR:
//
//    Triangle  Neighboring Triangles
//    --------  ---------------------
//
//     1        -1     -1      2
//     2         1      4      3
//     3         2      5      7
//     4         2     -1      8
//     5         3      8      6
//     6         5      9      7
//     7         3      6     -1
//     8         5      4     10
//     9         6     10     12
//    10         9      8     11
//    11        12     10     14
//    12         9     11     13
//    13        -1     12     16
//    14        11     -1     15
//    15        16     14     -1
//    16        13     15     -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint TRIANGLE_ORDER, the order of the triangles.
//
//    Input, gint TRIANGLE_NUM, the number of triangles.
//
//    Input, gint TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM], the nodes that
//    make up each triangle.
//
//    Output, gint TRIANGLE_ORDER3_NEIGHBOR_TRIANGLES[3*TRIANGLE_NUM],
//    the three triangles
//    that are direct neighbors of a given triangle.  TRIANGLE_NEIGHBOR(1,I)
//    is the index of the triangle which touches side 1, defined by nodes 2
//    and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative if there is no
//    neighbor on that side.  In this case, that side of the triangle lies
//    on the boundary of the triangulation.
//
{
  gint *col;
  gint i;
  gint icol;
  gint j;
  gint k;
  gint side1;
  gint side2;
  gint tri;
  gint tri1;
  gint tri2;
  gint *triangle_neighbor;

  triangle_neighbor = g_new0 (gint, 3 * triangle_num);
  col = g_new0 (gint, 4 * (3 * triangle_num));
//
//  Step 1.
//  From the list of nodes for triangle T, of the form: (I,J,K)
//  construct the three neighbor relations:
//
//    (I,J,3,T) or (J,I,3,T),
//    (J,K,1,T) or (K,J,1,T),
//    (K,I,2,T) or (I,K,2,T)
//
//  where we choose (I,J,3,T) if I < J, or else (J,I,3,T)
//
  for ( tri = 0; tri < triangle_num; tri++ )
  {
    i = triangle_node[0+tri*triangle_order];
    j = triangle_node[1+tri*triangle_order];
    k = triangle_node[2+tri*triangle_order];

    if ( i < j )
    {
      col[0+(3*tri+0)*4] = i;
      col[1+(3*tri+0)*4] = j;
      col[2+(3*tri+0)*4] = 3;
      col[3+(3*tri+0)*4] = tri + 1;
    }
    else
    {
      col[0+(3*tri+0)*4] = j;
      col[1+(3*tri+0)*4] = i;
      col[2+(3*tri+0)*4] = 3;
      col[3+(3*tri+0)*4] = tri + 1;
    }

    if ( j < k )
    {
      col[0+(3*tri+1)*4] = j;
      col[1+(3*tri+1)*4] = k;
      col[2+(3*tri+1)*4] = 1;
      col[3+(3*tri+1)*4] = tri + 1;
    }
    else
    {
      col[0+(3*tri+1)*4] = k;
      col[1+(3*tri+1)*4] = j;
      col[2+(3*tri+1)*4] = 1;
      col[3+(3*tri+1)*4] = tri + 1;
    }

    if ( k < i )
    {
      col[0+(3*tri+2)*4] = k;
      col[1+(3*tri+2)*4] = i;
      col[2+(3*tri+2)*4] = 2;
      col[3+(3*tri+2)*4] = tri + 1;
    }
    else
    {
      col[0+(3*tri+2)*4] = i;
      col[1+(3*tri+2)*4] = k;
      col[2+(3*tri+2)*4] = 2;
      col[3+(3*tri+2)*4] = tri + 1;
    }
  }
//
//  Step 2. Perform an ascending dictionary sort on the neighbor relations.
//  We only intend to sort on rows 1 and 2; the routine we call here
//  sorts on rows 1 through 4 but that won't hurt us.
//
//  What we need is to find cases where two triangles share an edge.
//  Say they share an edge defined by the nodes I and J.  Then there are
//  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
//  we make sure that these two columns occur consecutively.  That will
//  make it easy to notice that the triangles are neighbors.
//
  i4col_sort_a ( 4, 3*triangle_num, col );
//
//  Step 3. Neighboring triangles show up as consecutive columns with
//  identical first two entries.  Whenever you spot this happening,
//  make the appropriate entries in TRIANGLE_NEIGHBOR.
//
  for ( j = 0; j < triangle_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      triangle_neighbor[i+j*3] = -1;
    }
  }

  icol = 1;

  for ( ; ; )
  {
    if ( 3 * triangle_num <= icol )
    {
      break;
    }

    if ( col[0+(icol-1)*4] != col[0+icol*4] ||
         col[1+(icol-1)*4] != col[1+icol*4] )
    {
      icol = icol + 1;
      continue;
    }

    side1 = col[2+(icol-1)*4];
    tri1 =  col[3+(icol-1)*4];
    side2 = col[2+ icol   *4];
    tri2 =  col[3+ icol   *4];

    triangle_neighbor[side1-1+(tri1-1)*3] = tri2;
    triangle_neighbor[side2-1+(tri2-1)*3] = tri1;

    icol = icol + 2;
  }

  g_free (col);

  return triangle_neighbor;
}
//****************************************************************************80

gint triangulation_order3_adj_count ( gint node_num, gint triangle_num,
  gint triangle_node[], gint triangle_neighbor[], gint adj_col[] )

/****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies in a triangulation.
//
//  Discussion:
//
//    This routine is called to count the adjacencies, so that the
//    appropriate amount of memory can be set aside for storage when
//    the adjacency structure is created.
//
//    The triangulation is assumed to involve 3-node triangles.
//
//    Two nodes are "adjacent" if they are both nodes in some triangle.
//    Also, a node is considered to be adjacent to itself.
//
//  Diagram:
//
//       3
//    s  |\
//    i  | \
//    d  |  \
//    e  |   \  side 2
//       |    \
//    3  |     \
//       |      \
//       1-------2
//
//         side 1
//
//    The local node numbering
//
//
//   21-22-23-24-25
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   16-17-18-19-20
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   11-12-13-14-15
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    6--7--8--9-10
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    1--2--3--4--5
//
//    A sample grid.
//
//
//    Below, we have a chart that summarizes the adjacency relationships
//    in the sample grid.  On the left, we list the node, and its neighbors,
//    with an asterisk to indicate the adjacency of the node to itself
//    (in some cases, you want to count this self adjacency and in some
//    you don't).  On the right, we list the number of adjancencies to
//    lower-indexed nodes, to the node itself, to higher-indexed nodes,
//    the total number of adjacencies for this node, and the location
//    of the first and last entries required to list this set of adjacencies
//    in a single list of all the adjacencies.
//
//    N   Adjacencies                Below  Self   Above   Total First  Last
//
//   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
//    1:  *  2  6                        0     1       2       3     1     3
//    2:  1  *  3  6  7                  1     1       3       5     4     8
//    3:  2  *  4  7  8                  1     1       3       5     9    13
//    4:  3  *  5  8  9                  1     1       3       5    14    18
//    5:  4  *  9 10                     1     1       2       4    19    22
//    6:  1  2  *  7 11                  2     1       2       5    23    27
//    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
//    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
//    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
//   10:  5  9  * 14 15                  2     1       2       5    49    53
//   11:  6  7  * 12 16                  2     1       2       5    54    58
//   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
//   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
//   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
//   15: 10 14  * 19 20                  2     1       2       5    80    84
//   16: 11 12  * 17 21                  2     1       2       5    85    89
//   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
//   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
//   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
//   20: 15 19  * 24 25                  2     1       2       5   111   115
//   21: 16 17  * 22                     2     1       1       4   116   119
//   22: 17 18 21  * 23                  3     1       1       5   120   124
//   23: 18 19 22  * 24                  3     1       1       5   125   129
//   24: 19 20 23  * 25                  3     1       1       5   130   134
//   25: 20 24  *                        2     1       0       3   135   137
//   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint TRIANGLE_NUM, the number of triangles.
//
//    Input, gint TRIANGLE_NODE[3*TRIANGLE_NUM], lists the nodes that
//    make up each triangle, in counterclockwise order.
//
//    Input, gint TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
//    a triangle, lists the neighboring triangle, or -1 if there is
//    no neighbor.
//
//    Output, TRIANGULATION_ORDER3_ADJ_COUNT, the number of adjacencies.
//
//    Output, gint ADJ_COL[NODE_NUM+1].  Information about column J is stored
//    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
//*/
{
  gint adj_num;
  gint i;
  gint n1;
  gint n2;
  gint n3;
  gint node;
  gint triangle;
  gint triangle_order = 3;
  gint triangle2;

  adj_num = 0;
//
//  Set every node to be adjacent to itself.
//
  for ( node = 0; node < node_num; node++ )
  {
    adj_col[node] = 1;
  }
//
//  Examine each triangle.
//
  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    n1 = triangle_node[0+triangle*triangle_order];
    n2 = triangle_node[1+triangle*triangle_order];
    n3 = triangle_node[2+triangle*triangle_order];
//
//  Add edge (1,2) if this is the first occurrence,
//  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
//  or if this triangle is the first of the pair in which the edge
//  occurs (TRIANGLE < TRIANGLE2).
//
    triangle2 = triangle_neighbor[0+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[n1-1] = adj_col[n1-1] + 1;
      adj_col[n2-1] = adj_col[n2-1] + 1;
    }
//
//  Add edge (2,3).
//
    triangle2 = triangle_neighbor[1+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[n2-1] = adj_col[n2-1] + 1;
      adj_col[n3-1] = adj_col[n3-1] + 1;
    }
//
//  Add edge (3,1).
//
    triangle2 = triangle_neighbor[2+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[n1-1] = adj_col[n1-1] + 1;
      adj_col[n3-1] = adj_col[n3-1] + 1;
    }
  }
//
//  We used ADJ_COL to count the number of entries in each column.
//  Convert it to pointers into the ADJ array.
//
  for ( node = node_num; 1 <= node; node-- )
  {
    adj_col[node] = adj_col[node-1];
  }
  adj_col[0] = 1;
  for ( i = 1; i <= node_num; i++ )
  {
    adj_col[i]= adj_col[i-1] + adj_col[i];
  }

  adj_num = adj_col[node_num] - 1;

  return adj_num;
}
//****************************************************************************80

gint *triangulation_order3_adj_set ( gint node_num, gint triangle_num,
  gint triangle_node[], gint triangle_neighbor[], gint adj_num, gint adj_col[] )

/****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_ADJ_SET sets adjacencies in a triangulation.
//
//  Discussion:
//
//    This routine is called to count the adjacencies, so that the
//    appropriate amount of memory can be set aside for storage when
//    the adjacency structure is created.
//
//    The triangulation is assumed to involve 3-node triangles.
//
//    Two nodes are "adjacent" if they are both nodes in some triangle.
//    Also, a node is considered to be adjacent to itself.
//
//    This routine can be used to create the compressed column storage
//    for a linear triangle finite element discretization of
//    Poisson's equation in two dimensions.
//
//  Diagram:
//
//       3
//    s  |\
//    i  | \
//    d  |  \
//    e  |   \  side 2
//       |    \
//    3  |     \
//       |      \
//       1-------2
//
//         side 1
//
//    The local node numbering
//
//
//   21-22-23-24-25
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   16-17-18-19-20
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   11-12-13-14-15
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    6--7--8--9-10
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    1--2--3--4--5
//
//    A sample grid
//
//
//    Below, we have a chart that summarizes the adjacency relationships
//    in the sample grid.  On the left, we list the node, and its neighbors,
//    with an asterisk to indicate the adjacency of the node to itself
//    (in some cases, you want to count this self adjacency and in some
//    you don't).  On the right, we list the number of adjancencies to
//    lower-indexed nodes, to the node itself, to higher-indexed nodes,
//    the total number of adjacencies for this node, and the location
//    of the first and last entries required to list this set of adjacencies
//    in a single list of all the adjacencies.
//
//    N   Adjacencies                Below  Self    Above  Total First  Last
//
//   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
//    1:  *  2  6                        0     1       2       3     1     3
//    2:  1  *  3  6  7                  1     1       3       5     4     8
//    3:  2  *  4  7  8                  1     1       3       5     9    13
//    4:  3  *  5  8  9                  1     1       3       5    14    18
//    5:  4  *  9 10                     1     1       2       4    19    22
//    6:  1  2  *  7 11                  2     1       2       5    23    27
//    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
//    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
//    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
//   10:  5  9  * 14 15                  2     1       2       5    49    53
//   11:  6  7  * 12 16                  2     1       2       5    54    58
//   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
//   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
//   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
//   15: 10 14  * 19 20                  2     1       2       5    80    84
//   16: 11 12  * 17 21                  2     1       2       5    85    89
//   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
//   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
//   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
//   20: 15 19  * 24 25                  2     1       2       5   111   115
//   21: 16 17  * 22                     2     1       1       4   116   119
//   22: 17 18 21  * 23                  3     1       1       5   120   124
//   23: 18 19 22  * 24                  3     1       1       5   125   129
//   24: 19 20 23  * 25                  3     1       1       5   130   134
//   25: 20 24  *                        2     1       0       3   135   137
//   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint TRIANGLE_NUM, the number of triangles.
//
//    Input, gint TRIANGLE_NODE[3*TRIANGLE_NUM], lists the nodes that
//    make up each triangle in counterclockwise order.
//
//    Input, gint TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
//    a triangle, lists the neighboring triangle, or -1 if there is
//    no neighbor.
//
//    Input, gint ADJ_NUM, the number of adjacencies.
//
//    Input, gint ADJ_COL[NODE_NUM+1].  Information about column J is stored
//    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
//
//    Output, gint TRIANGULATION_ORDER3_ADJ_SET[ADJ_NUM], the adjacency
//    information.
//*/
{
  gint *adj;
  gint *adj_copy;
  gint k;
  gint k1;
  gint k2;
  gint n1;
  gint n2;
  gint n3;
  gint node;
  gint triangle;
  gint triangle2;
  gint triangle_order = 3;

  adj = g_new0 (gint, adj_num);
  for ( k = 0; k < adj_num; k++ )
  {
    adj[k] = -1;
  }

  adj_copy = g_new0 (gint, node_num);
  for ( node = 0; node < node_num; node++ )
  {
    adj_copy[node] = adj_col[node];
  }
//
//  Set every node to be adjacent to itself.
//
  for ( node = 1; node <= node_num; node++ )
  {
    adj[adj_copy[node-1]-1] = node;
    adj_copy[node-1] = adj_copy[node-1] + 1;
  }
//
//  Examine each triangle.
//
  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    n1 = triangle_node[0+triangle*triangle_order];
    n2 = triangle_node[1+triangle*triangle_order];
    n3 = triangle_node[2+triangle*triangle_order];
//
//  Add edge (1,2) if this is the first occurrence,
//  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
//  or if this triangle is the first of the pair in which the edge
//  occurs (TRIANGLE < TRIANGLE2).
//
    triangle2 = triangle_neighbor[0+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj[adj_copy[n1-1]-1] = n2;
      adj_copy[n1-1] = adj_copy[n1-1] + 1;
      adj[adj_copy[n2-1]-1] = n1;
      adj_copy[n2-1] = adj_copy[n2-1] + 1;
    }
//
//  Add edge (2,3).
//
    triangle2 = triangle_neighbor[1+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj[adj_copy[n2-1]-1] = n3;
      adj_copy[n2-1] = adj_copy[n2-1] + 1;
      adj[adj_copy[n3-1]-1] = n2;
      adj_copy[n3-1] = adj_copy[n3-1] + 1;
    }
//
//  Add edge (3,1).
//
    triangle2 = triangle_neighbor[2+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj[adj_copy[n1-1]-1] = n3;
      adj_copy[n1-1] = adj_copy[n1-1] + 1;
      adj[adj_copy[n3-1]-1] = n1;
      adj_copy[n3-1] = adj_copy[n3-1] + 1;
    }
  }
//
//  Ascending sort the entries for each node.
//
  for ( node = 1; node <= node_num; node++ )
  {
    k1 = adj_col[node-1];
    k2 = adj_col[node]-1;
    i4vec_sort_heap_a ( k2+1-k1, adj+k1-1 );
  }

  g_free (adj_copy);

  return adj;
}
//****************************************************************************80

void triangulation_order3_example2 ( gint node_num, gint triangle_num,
  gdouble node_xy[], gint triangle_node[], gint triangle_neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_EXAMPLE2 sets up a sample triangulation.
//
//  Discussion:
//
//    This triangulation is actually a Delaunay triangulation.
//
//    The appropriate input values of NODE_NUM and TRIANGLE_NUM can be
//    determined by calling TRIANGULATION_ORDER3_EXAMPLE2_SIZE first.
//
//  Diagram:
//
//   21-22-23-24-25
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   16-17-18-19-20
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   11-12-13-14-15
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    6--7--8--9-10
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    1--2--3--4--5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint TRIANGLE_NUM, the number of triangles.
//
//    Output, gdouble NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Output, gint TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the triangles.
//
//    Output, gint TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbors on each side.
//    Negative values indicate edges that lie on the exterior.
//
{
# define DIM_NUM 2
# define NODE_NUM 25
# define TRIANGLE_NUM 32
# define TRIANGLE_ORDER 3

  gint i;
  static gint triangle_neighbor_save[3*TRIANGLE_NUM] = {
    -1,  2, -1,
     9,  1,  3,
    -1,  4,  2,
    11,  3,  5,
    -1,  6,  4,
    13,  5,  7,
    -1,  8,  6,
    15,  7, -1,
     2, 10, -1,
    17,  9, 11,
     4, 12, 10,
    19, 11, 13,
     6, 14, 12,
    21, 13, 15,
     8, 16, 14,
    23, 15, -1,
    10, 18, -1,
    25, 17, 19,
    12, 20, 18,
    27, 19, 21,
    14, 22, 20,
    29, 21, 23,
    16, 24, 22,
    31, 23, -1,
    18, 26, -1,
    -1, 25, 27,
    20, 28, 26,
    -1, 27, 29,
    22, 30, 28,
    -1, 29, 31,
    24, 32, 30,
    -1, 31, -1 };
  static gint triangle_node_save[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     1,  2,  6,
     7,  6,  2,
     2,  3,  7,
     8,  7,  3,
     3,  4,  8,
     9,  8,  4,
     4,  5,  9,
    10,  9,  5,
     6,  7, 11,
    12, 11,  7,
     7,  8, 12,
    13, 12,  8,
     8,  9, 13,
    14, 13,  9,
     9, 10, 14,
    15, 14, 10,
    11, 12, 16,
    17, 16, 12,
    12, 13, 17,
    18, 17, 13,
    13, 14, 18,
    19, 18, 14,
    14, 15, 19,
    20, 19, 15,
    16, 17, 21,
    22, 21, 17,
    17, 18, 22,
    23, 22, 18,
    18, 19, 23,
    24, 23, 19,
    19, 20, 24,
    25, 24, 20 };
  static gdouble node_xy_save[DIM_NUM*NODE_NUM] = {
    0.0, 0.0,
    1.0, 0.0,
    2.0, 0.0,
    3.0, 0.0,
    4.0, 0.0,
    0.0, 1.0,
    1.0, 1.0,
    2.0, 1.0,
    3.0, 1.0,
    4.0, 1.0,
    0.0, 2.0,
    1.0, 2.0,
    2.0, 2.0,
    3.0, 2.0,
    4.0, 2.0,
    0.0, 3.0,
    1.0, 3.0,
    2.0, 3.0,
    3.0, 3.0,
    4.0, 3.0,
    0.0, 4.0,
    1.0, 4.0,
    2.0, 4.0,
    3.0, 4.0,
    4.0, 4.0  };

  for ( i = 0; i < 3 * TRIANGLE_NUM; i++ )
  {
    triangle_neighbor[i] = triangle_neighbor_save[i];
  }

  for ( i = 0; i < TRIANGLE_ORDER * TRIANGLE_NUM; i++ )
  {
    triangle_node[i] = triangle_node_save[i];
  }

  for ( i = 0; i < DIM_NUM * NODE_NUM; i++ )
  {
    node_xy[i] = node_xy_save[i];
  }

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void triangulation_order3_example2_size ( gint *node_num, gint *triangle_num,
  gint *hole_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_EXAMPLE2_SIZE sets sizes for a sample triangulation.
//
//  Diagram:
//
//   21-22-23-24-25
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   16-17-18-19-20
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   11-12-13-14-15
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    6--7--8--9-10
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    1--2--3--4--5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, gint *NODE_NUM, the number of nodes.
//
//    Output, gint *TRIANGLE_NUM, the number of triangles.
//
//    Output, gint *HOLE_NUM, the number of holes.
//
{
  *node_num = 25;
  *triangle_num = 32;
  *hole_num = 0;

  return;
}
//****************************************************************************80

gint triangulation_order6_adj_count ( gint node_num, gint triangle_num,
  gint triangle_node[], gint triangle_neighbor[], gint adj_col[] )

/****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_ADJ_COUNT counts adjacencies in a triangulation.
//
//  Discussion:
//
//    This routine is called to count the adjacencies, so that the
//    appropriate amount of memory can be set aside for storage when
//    the adjacency structure is created.
//
//    The triangulation is assumed to involve 6-node triangles.
//
//    Two nodes are "adjacent" if they are both nodes in some triangle.
//    Also, a node is considered to be adjacent to itself.
//
//  Diagram:
//
//       3
//    s  |\
//    i  | \
//    d  |  \
//    e  6   5  side 2
//       |    \
//    3  |     \
//       |      \
//       1---4---2
//
//         side 1
//
//    The local node numbering
//
//
//   21-22-23-24-25
//    |\    |\    |
//    | \   | \   |
//   16 17 18 19 20
//    |   \ |   \ |
//    |    \|    \|
//   11-12-13-14-15
//    |\    |\    |
//    | \   | \   |
//    6  7  8  9 10
//    |   \ |   \ |
//    |    \|    \|
//    1--2--3--4--5
//
//    A sample grid.
//
//
//    Below, we have a chart that lists the nodes adjacent to each node, with
//    an asterisk to indicate the adjacency of the node to itself
//    (in some cases, you want to count this self adjacency and in some
//    you don't).
//
//    N   Adjacencies
//
//    1:  *  2  3  6  7 11
//    2:  1  *  3  6  7 11
//    3:  1  2  *  4  5  6  7  8  9 11 12 13
//    4:  3  *  5  8  9 13
//    5:  3  4  *  8  9 10 13 14 15
//    6:  1  2  3  *  7 11
//    7:  1  2  3  6  *  8 11 12 13
//    8:  3  4  5  7  *  9 11 12 13
//    9:  3  4  5  8  * 10 13 14 15
//   10:  5  9  * 13 14 15
//   11:  1  2  3  6  7  8  * 12 13 16 17 21
//   12:  3  7  8 11  * 13 16 17 21
//   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
//   14:  5  9 10 13  * 15 18 19 23
//   15:  5  9 10 13 14  * 18 19 20 23 24 25
//   16: 11 12 13  * 17 21
//   17: 11 12 13 16  * 18 21 22 23
//   18: 13 14 15 17  * 19 21 22 23
//   19: 13 14 15 18  * 20 23 24 25
//   20: 15 19  * 23 24 25
//   21: 11 12 13 16 17 18  * 22 23
//   22: 13 17 18 21  * 23
//   23: 13 14 15 17 18 19 20 21 22  * 24 25
//   24: 15 19 20 23  * 25
//   25: 15 19 20 23 24  *
//
//    Below, we list the number of adjancencies to lower-indexed nodes, to
//    the node itself, to higher-indexed nodes, the total number of
//    adjacencies for this node, and the location of the first and last
//    entries required to list this set of adjacencies in a single list
//    of all the adjacencies.
//
//    N   Below  Self   Above   Total First  Last
//
//   --      --    --      --      --   ---     0
//    1:      0     1       5       6     1     6
//    2:      1     1       4       6     7    12
//    3:      2     1       9      12    13    24
//    4:      1     1       4       6    25    30
//    5:      2     1       6       9    31    39
//    6:      3     1       2       6    40    45
//    7:      4     1       4       9    46    54
//    8:      4     1       4       9    55    63
//    9:      4     1       4       9    62    72
//   10:      2     1       3       6    73    78
//   11:      6     1       5      12    79    90
//   12:      4     1       4       9    91    99
//   13:      9     1       9      19   100   118
//   14:      4     1       4       9   119   127
//   15:      5     1       6      12   128   139
//   16:      3     1       2       6   140   145
//   17:      4     1       4       9   146   154
//   18:      4     1       4       9   155   163
//   19:      4     1       4       9   164   172
//   20:      2     1       3       6   173   178
//   21:      6     1       2       9   179   187
//   22:      4     1       1       6   188   193
//   23:      9     1       2      12   194   205
//   24:      4     1       1       6   206   211
//   25:      5     1       0       6   212   217
//   --      --    --      --      --   218   ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint TRIANGLE_NUM, the number of triangles.
//
//    Input, gint TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
//    make up each triangle.  The first three nodes are the vertices,
//    in counterclockwise order.  The fourth value is the midside
//    node between nodes 1 and 2; the fifth and sixth values are
//    the other midside nodes in the logical order.
//
//    Input, gint TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
//    a triangle, lists the neighboring triangle, or -1 if there is
//    no neighbor.
//
//    Output, gint TRIANGULATION_ORDER6_ADJ_COUNT, the number of adjacencies.
//
//    Output, gint ADJ_COL[NODE_NUM+1].  Information about column J is stored
//    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
//*/
{
  gint adj_num;
  gint i;
  gint n1;
  gint n2;
  gint n3;
  gint n4;
  gint n5;
  gint n6;
  gint node;
  gint triangle;
  gint triangle_order = 6;
  gint triangle2;

  adj_num = 0;
//
//  Set every node to be adjacent to itself.
//
  for ( node = 0; node < node_num; node++ )
  {
    adj_col[node] = 1;
  }
//
//  Examine each triangle.
//
  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    n1 = triangle_node[0+triangle*triangle_order];
    n2 = triangle_node[1+triangle*triangle_order];
    n3 = triangle_node[2+triangle*triangle_order];
    n4 = triangle_node[3+triangle*triangle_order];
    n5 = triangle_node[4+triangle*triangle_order];
    n6 = triangle_node[5+triangle*triangle_order];
//
//  For sure, we add the adjacencies:
//    43 / (34)
//    51 / (15)
//    54 / (45)
//    62 / (26)
//    64 / (46)
//    65 / (56)
//
    adj_col[n3-1] = adj_col[n3-1] + 1;
    adj_col[n4-1] = adj_col[n4-1] + 1;
    adj_col[n1-1] = adj_col[n1-1] + 1;
    adj_col[n5-1] = adj_col[n5-1] + 1;
    adj_col[n4-1] = adj_col[n4-1] + 1;
    adj_col[n5-1] = adj_col[n5-1] + 1;
    adj_col[n2-1] = adj_col[n2-1] + 1;
    adj_col[n6-1] = adj_col[n6-1] + 1;
    adj_col[n4-1] = adj_col[n4-1] + 1;
    adj_col[n6-1] = adj_col[n6-1] + 1;
    adj_col[n5-1] = adj_col[n5-1] + 1;
    adj_col[n6-1] = adj_col[n6-1] + 1;
//
//  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
//  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
//  or if this triangle is the first of the pair in which the edge
//  occurs (TRIANGLE < TRIANGLE2).
//
//  Maybe add
//    21 / 12
//    41 / 14
//    42 / 24
//
    triangle2 = triangle_neighbor[0+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[n1-1] = adj_col[n1-1] + 1;
      adj_col[n2-1] = adj_col[n2-1] + 1;
      adj_col[n1-1] = adj_col[n1-1] + 1;
      adj_col[n4-1] = adj_col[n4-1] + 1;
      adj_col[n2-1] = adj_col[n2-1] + 1;
      adj_col[n4-1] = adj_col[n4-1] + 1;
    }
//
//  Maybe add
//    32 / 23
//    52 / 25
//    53 / 35
//
    triangle2 = triangle_neighbor[1+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[n2-1] = adj_col[n2-1] + 1;
      adj_col[n3-1] = adj_col[n3-1] + 1;
      adj_col[n2-1] = adj_col[n2-1] + 1;
      adj_col[n5-1] = adj_col[n5-1] + 1;
      adj_col[n3-1] = adj_col[n3-1] + 1;
      adj_col[n5-1] = adj_col[n5-1] + 1;
    }
//
//  Maybe add
//    31 / 13
//    61 / 16
//    63 / 36
//
    triangle2 = triangle_neighbor[2+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[n1-1] = adj_col[n1-1] + 1;
      adj_col[n3-1] = adj_col[n3-1] + 1;
      adj_col[n1-1] = adj_col[n1-1] + 1;
      adj_col[n6-1] = adj_col[n6-1] + 1;
      adj_col[n3-1] = adj_col[n3-1] + 1;
      adj_col[n6-1] = adj_col[n6-1] + 1;
    }
  }
//
//  We used ADJ_COL to count the number of entries in each column.
//  Convert it to pointers into the ADJ array.
//
  for ( node = node_num; 1 <= node; node-- )
  {
    adj_col[node] = adj_col[node-1];
  }
  adj_col[0] = 1;
  for ( i = 1; i <= node_num; i++ )
  {
    adj_col[i]= adj_col[i-1] + adj_col[i];
  }

  adj_num = adj_col[node_num] - 1;

  return adj_num;
}
//****************************************************************************80

gint *triangulation_order6_adj_set ( gint node_num, gint triangle_num,
  gint triangle_node[], gint triangle_neighbor[], gint adj_num, gint adj_col[] )

/****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_ADJ_SET sets adjacencies in a triangulation.
//
//  Discussion:
//
//    This routine is called to count the adjacencies, so that the
//    appropriate amount of memory can be set aside for storage when
//    the adjacency structure is created.
//
//    The triangulation is assumed to involve 6-node triangles.
//
//    Two nodes are "adjacent" if they are both nodes in some triangle.
//    Also, a node is considered to be adjacent to itself.
//
//    This routine can be used to create the compressed column storage
//    for a quadratic triangle finite element discretization of
//    Poisson's equation in two dimensions.
//
//  Diagram:
//
//       3
//    s  |\
//    i  | \
//    d  |  \
//    e  6   5  side 2
//       |    \
//    3  |     \
//       |      \
//       1---4---2
//
//         side 1
//
//    The local node numbering
//
//
//   21-22-23-24-25
//    |\    |\    |
//    | \   | \   |
//   16 17 18 19 20
//    |   \ |   \ |
//    |    \|    \|
//   11-12-13-14-15
//    |\    |\    |
//    | \   | \   |
//    6  7  8  9 10
//    |   \ |   \ |
//    |    \|    \|
//    1--2--3--4--5
//
//    A sample grid.
//
//
//    Below, we have a chart that lists the nodes adjacent to each node, with
//    an asterisk to indicate the adjacency of the node to itself
//    (in some cases, you want to count this self adjacency and in some
//    you don't).
//
//    N   Adjacencies
//
//    1:  *  2  3  6  7 11
//    2:  1  *  3  6  7 11
//    3:  1  2  *  4  5  6  7  8  9 11 12 13
//    4:  3  *  5  8  9 13
//    5:  3  4  *  8  9 10 13 14 15
//    6:  1  2  3  *  7 11
//    7:  1  2  3  6  *  8 11 12 13
//    8:  3  4  5  7  *  9 11 12 13
//    9:  3  4  5  8  * 10 13 14 15
//   10:  5  9  * 13 14 15
//   11:  1  2  3  6  7  8  * 12 13 16 17 21
//   12:  3  7  8 11  * 13 16 17 21
//   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
//   14:  5  9 10 13  * 15 18 19 23
//   15:  5  9 10 13 14  * 18 19 20 23 24 25
//   16: 11 12 13  * 17 21
//   17: 11 12 13 16  * 18 21 22 23
//   18: 13 14 15 17  * 19 21 22 23
//   19: 13 14 15 18  * 20 23 24 25
//   20: 15 19  * 23 24 25
//   21: 11 12 13 16 17 18  * 22 23
//   22: 13 17 18 21  * 23
//   23: 13 14 15 17 18 19 20 21 22  * 24 25
//   24: 15 19 20 23  * 25
//   25: 15 19 20 23 24  *
//
//    Below, we list the number of adjancencies to lower-indexed nodes, to
//    the node itself, to higher-indexed nodes, the total number of
//    adjacencies for this node, and the location of the first and last
//    entries required to list this set of adjacencies in a single list
//    of all the adjacencies.
//
//    N   Below  Self   Above   Total First  Last
//
//   --      --    --      --      --   ---     0
//    1:      0     1       5       6     1     6
//    2:      1     1       4       6     7    12
//    3:      2     1       9      12    13    24
//    4:      1     1       4       6    25    30
//    5:      2     1       6       9    31    39
//    6:      3     1       2       6    40    45
//    7:      4     1       4       9    46    54
//    8:      4     1       4       9    55    63
//    9:      4     1       4       9    62    72
//   10:      2     1       3       6    73    78
//   11:      6     1       5      12    79    90
//   12:      4     1       4       9    91    99
//   13:      9     1       9      19   100   118
//   14:      4     1       4       9   119   127
//   15:      5     1       6      12   128   139
//   16:      3     1       2       6   140   145
//   17:      4     1       4       9   146   154
//   18:      4     1       4       9   155   163
//   19:      4     1       4       9   164   172
//   20:      2     1       3       6   173   178
//   21:      6     1       2       9   179   187
//   22:      4     1       1       6   188   193
//   23:      9     1       2      12   194   205
//   24:      4     1       1       6   206   211
//   25:      5     1       0       6   212   217
//   --      --    --      --      --   218   ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint TRIANGLE_NUM, the number of triangles.
//
//    Input, gint TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
//    make up each triangle.  The first three nodes are the vertices,
//    in counterclockwise order.  The fourth value is the midside
//    node between nodes 1 and 2; the fifth and sixth values are
//    the other midside nodes in the logical order.
//
//    Input, gint TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
//    a triangle, lists the neighboring triangle, or -1 if there is
//    no neighbor.
//
//    Input, gint ADJ_NUM, the number of adjacencies.
//
//    Input, gint ADJ_COL[NODE_NUM+1].  Information about column J is stored
//    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
//
//    Output, gint TRIANGULATION_ORDER6_ADJ_SET[ADJ_NUM], the adjacency
//    information.
//*/
{
  gint *adj;
  gint *adj_copy;
  gint k;
  gint k1;
  gint k2;
  gint n1;
  gint n2;
  gint n3;
  gint n4;
  gint n5;
  gint n6;
  gint node;
  gint triangle;
  gint triangle2;
  gint triangle_order = 6;

  adj = g_new0 (gint, adj_num);
  for ( k = 0; k < adj_num; k++ )
  {
    adj[k] = -1;
  }

  adj_copy = g_new0 (gint, node_num);
  for ( node = 0; node < node_num; node++ )
  {
    adj_copy[node] = adj_col[node];
  }
//
//  Set every node to be adjacent to itself.
//
  for ( node = 1; node <= node_num; node++ )
  {
    adj[adj_copy[node-1]-1] = node;
    adj_copy[node-1] = adj_copy[node-1] + 1;
  }
//
//  Examine each triangle.
//
  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    n1 = triangle_node[0+triangle*triangle_order];
    n2 = triangle_node[1+triangle*triangle_order];
    n3 = triangle_node[2+triangle*triangle_order];
    n4 = triangle_node[3+triangle*triangle_order];
    n5 = triangle_node[4+triangle*triangle_order];
    n6 = triangle_node[5+triangle*triangle_order];
//
//  For sure, we add the adjacencies:
//    43 / (34)
//    51 / (15)
//    54 / (45)
//    62 / (26)
//    64 / (46)
//    65 / (56)
//
    adj[adj_copy[n3-1]-1] = n4;
    adj_copy[n3-1] = adj_copy[n3-1] + 1;
    adj[adj_copy[n4-1]-1] = n3;
    adj_copy[n4-1] = adj_copy[n4-1] + 1;

    adj[adj_copy[n1-1]-1] = n5;
    adj_copy[n1-1] = adj_copy[n1-1] + 1;
    adj[adj_copy[n5-1]-1] = n1;
    adj_copy[n5-1] = adj_copy[n5-1] + 1;

    adj[adj_copy[n4-1]-1] = n5;
    adj_copy[n4-1] = adj_copy[n4-1] + 1;
    adj[adj_copy[n5-1]-1] = n4;
    adj_copy[n5-1] = adj_copy[n5-1] + 1;

    adj[adj_copy[n2-1]-1] = n6;
    adj_copy[n2-1] = adj_copy[n2-1] + 1;
    adj[adj_copy[n6-1]-1] = n2;
    adj_copy[n6-1] = adj_copy[n6-1] + 1;

    adj[adj_copy[n4-1]-1] = n6;
    adj_copy[n4-1] = adj_copy[n4-1] + 1;
    adj[adj_copy[n6-1]-1] = n4;
    adj_copy[n6-1] = adj_copy[n6-1] + 1;

    adj[adj_copy[n5-1]-1] = n6;
    adj_copy[n5-1] = adj_copy[n5-1] + 1;
    adj[adj_copy[n6-1]-1] = n5;
    adj_copy[n6-1] = adj_copy[n6-1] + 1;
//
//  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
//  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
//  or if this triangle is the first of the pair in which the edge
//  occurs (TRIANGLE < TRIANGLE2).
//
//  Maybe add
//    21 / 12
//    41 / 14
//    42 / 24
//
    triangle2 = triangle_neighbor[0+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj[adj_copy[n1-1]-1] = n2;
      adj_copy[n1-1] = adj_copy[n1-1] + 1;
      adj[adj_copy[n2-1]-1] = n1;
      adj_copy[n2-1] = adj_copy[n2-1] + 1;
      adj[adj_copy[n1-1]-1] = n4;
      adj_copy[n1-1] = adj_copy[n1-1] + 1;
      adj[adj_copy[n4-1]-1] = n1;
      adj_copy[n4-1] = adj_copy[n4-1] + 1;
      adj[adj_copy[n2-1]-1] = n4;
      adj_copy[n2-1] = adj_copy[n2-1] + 1;
      adj[adj_copy[n4-1]-1] = n2;
      adj_copy[n4-1] = adj_copy[n4-1] + 1;
    }
//
//  Maybe add
//    32 / 23
//    52 / 25
//    53 / 35
//
    triangle2 = triangle_neighbor[1+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj[adj_copy[n2-1]-1] = n3;
      adj_copy[n2-1] = adj_copy[n2-1] + 1;
      adj[adj_copy[n3-1]-1] = n2;
      adj_copy[n3-1] = adj_copy[n3-1] + 1;
      adj[adj_copy[n2-1]-1] = n5;
      adj_copy[n2-1] = adj_copy[n2-1] + 1;
      adj[adj_copy[n5-1]-1] = n2;
      adj_copy[n5-1] = adj_copy[n5-1] + 1;
      adj[adj_copy[n3-1]-1] = n5;
      adj_copy[n3-1] = adj_copy[n3-1] + 1;
      adj[adj_copy[n5-1]-1] = n3;
      adj_copy[n5-1] = adj_copy[n5-1] + 1;
    }
//
//  Maybe add
//    31 / 13
//    61 / 16
//    63 / 36
//
    triangle2 = triangle_neighbor[2+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj[adj_copy[n1-1]-1] = n3;
      adj_copy[n1-1] = adj_copy[n1-1] + 1;
      adj[adj_copy[n3-1]-1] = n1;
      adj_copy[n3-1] = adj_copy[n3-1] + 1;
      adj[adj_copy[n1-1]-1] = n6;
      adj_copy[n1-1] = adj_copy[n1-1] + 1;
      adj[adj_copy[n6-1]-1] = n1;
      adj_copy[n6-1] = adj_copy[n6-1] + 1;
      adj[adj_copy[n3-1]-1] = n6;
      adj_copy[n3-1] = adj_copy[n3-1] + 1;
      adj[adj_copy[n6-1]-1] = n3;
      adj_copy[n6-1] = adj_copy[n6-1] + 1;
    }
  }
//
//  Ascending sort the entries for each node.
//
  for ( node = 1; node <= node_num; node++ )
  {
    k1 = adj_col[node-1];
    k2 = adj_col[node]-1;
    i4vec_sort_heap_a ( k2+1-k1, adj+k1-1 );
  }

  g_free (adj_copy);

  return adj;
}
//****************************************************************************80

void triangulation_order6_example2 ( gint node_num, gint triangle_num,
  gdouble node_xy[], gint triangle_node[], gint triangle_neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_EXAMPLE2 sets up a sample triangulation.
//
//  Discussion:
//
//    This triangulation is actually a Delaunay triangulation.
//
//    The appropriate input values of NODE_NUM and TRIANGLE_NUM can be
//    determined by calling TRIANGULATION_ORDER6_EXAMPLE2_SIZE first.
//
//  Diagram:
//
//   21-22-23-24-25
//    |\  6 |\  8 |
//    | \   | \   |
//   16 17 18 19 20
//    |   \ |   \ |
//    | 5  \| 7  \|
//   11-12-13-14-15
//    |\  2 |\  4 |
//    | \   | \   |
//    6  7  8  9 10
//    | 1 \ | 3 \ |
//    |    \|    \|
//    1--2--3--4--5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, gint NODE_NUM, the number of nodes.
//
//    Input, gint TRIANGLE_NUM, the number of triangles.
//
//    Output, gdouble NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Output, gint TRIANGLE_ORDER[6*TRIANGLE_NUM], the nodes that make up
//    the triangles.
//
//    Output, gint TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbors
//    on each side.  Negative values indicate edges that lie on the exterior.
//
{
# define DIM_NUM 2
# define NODE_NUM 25
# define TRIANGLE_NUM 8
# define TRIANGLE_ORDER 6

  gint i;
  gint j;
  static gdouble node_xy_save[DIM_NUM*NODE_NUM] = {
    0.0, 0.0,
    1.0, 0.0,
    2.0, 0.0,
    3.0, 0.0,
    4.0, 0.0,
    0.0, 1.0,
    1.0, 1.0,
    2.0, 1.0,
    3.0, 1.0,
    4.0, 1.0,
    0.0, 2.0,
    1.0, 2.0,
    2.0, 2.0,
    3.0, 2.0,
    4.0, 2.0,
    0.0, 3.0,
    1.0, 3.0,
    2.0, 3.0,
    3.0, 3.0,
    4.0, 3.0,
    0.0, 4.0,
    1.0, 4.0,
    2.0, 4.0,
    3.0, 4.0,
    4.0, 4.0 };
  static gint triangle_node_save[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     1,  3, 11,  2,  7,  6,
    13, 11,  3, 12,  7,  8,
     3,  5, 13,  4,  9,  8,
    15, 13,  5, 14,  9, 10,
    11, 13, 21, 12, 17, 16,
    23, 21, 13, 22, 17, 18,
    13, 15, 23, 14, 19, 18,
    25, 23, 15, 24, 19, 20  };
  static gint triangle_neighbor_save[3*TRIANGLE_NUM] = {
    -1,  2, -1,
     5,  1,  3,
    -1,  4,  2,
     7,  3, -1,
     2,  6, -1,
    -1,  5,  7,
     4,  8,  6,
    -1,  7, -1 };

  for ( j = 0; j < NODE_NUM; j++ )
  {
    for ( i = 0; i < DIM_NUM; i++ )
    {
      node_xy[i+j*DIM_NUM] = node_xy_save[i+j*DIM_NUM];
    }
  }

  for ( j = 0; j < TRIANGLE_NUM; j++ )
  {
    for ( i = 0; i < TRIANGLE_ORDER; i++ )
    {
      triangle_node[i+j*TRIANGLE_ORDER] = triangle_node_save[i+j*TRIANGLE_ORDER];
    }
  }

  for ( j = 0; j < TRIANGLE_NUM; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      triangle_neighbor[i+j*3] = triangle_neighbor_save[i+j*3];
    }
  }

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void triangulation_order6_example2_size ( gint *node_num, gint *triangle_num,
  gint *hole_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_EXAMPLE2_SIZE sets sizes for a sample triangulation.
//
//  Diagram:
//
//   21-22-23-24-25
//    |\  6 |\  8 |
//    | \   | \   |
//   16 17 18 19 20
//    |   \ |   \ |
//    | 5  \| 7  \|
//   11-12-13-14-15
//    |\  2 |\  4 |
//    | \   | \   |
//    6  7  8  9 10
//    | 1 \ | 3 \ |
//    |    \|    \|
//    1--2--3--4--5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, gint *NODE_NUM, the number of nodes.
//
//    Output, gint *TRIANGLE_NUM, the number of triangles.
//
//    Output, gint *HOLE_NUM, the number of holes.
//
{
  *node_num = 25;
  *triangle_num = 8;
  *hole_num = 0;

  return;
}
