/***************************************************************************
 *            mpq_tree.c
 *
 *  Mon Feb 22 13:59:24 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * SECTION:mpq_tree
 * @title: MPQ Tree
 * @short_description: FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <stdio.h>
#include <glib.h>
#include <gmp.h>
#include <mpfr.h>

static void
_slice_mpq_free (gpointer q)
{
  mpq_clear ((mpq_ptr)q);
  g_slice_free (mpq_t, q);
}

static gint
_slice_mpq_cmp (gconstpointer a, gconstpointer b, gpointer user_data)
{
  return mpq_cmp ((mpq_ptr)a, (mpq_ptr)b);
}

/**
 * ncm_mpq_tree_new: (skip)
 *
 * FIXME
 *
 * Returns: FIXME
*/
GTree *
ncm_mpq_tree_new (void)
{
  return g_tree_new_full (&_slice_mpq_cmp, NULL, &_slice_mpq_free, NULL);
 // return g_tree_new ();
}

/**
 * mpq_hash: 
 * @v: FIXME 
 *
 * FIXME
 *
 * Returns: FIXME
*/
guint   
mpq_hash (gconstpointer v)
{
  return mpz_get_ui (mpq_numref((mpq_ptr)v)) + mpz_get_ui (mpq_denref((mpq_ptr)v));
}

/**
 * ncm_mpq_hash_new: 
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
*/
GHashTable *
ncm_mpq_hash_new (void)
{
  return g_hash_table_new_full (&mpq_hash, (GEqualFunc) &mpq_equal, &_slice_mpq_free, NULL);
}
