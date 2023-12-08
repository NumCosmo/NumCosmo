/***************************************************************************
 *            ncm_dtuple.c
 *
 *  Fri December 08 10:28:16 2023
 *  Copyright  2023  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_dtuple.c
 * Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:ncm_dtuple
 * @title: NcmDTuple
 * @short_description: Fixed sized array of double values.
 *
 * In this module we define the #NcmDTuple2, #NcmDTuple3, which are fixed sized arrays
 * of double values.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_dtuple.h"

G_DEFINE_BOXED_TYPE (NcmDTuple2, ncm_dtuple2, ncm_dtuple2_copy, ncm_dtuple2_free)
G_DEFINE_BOXED_TYPE (NcmDTuple3, ncm_dtuple3, ncm_dtuple3_copy, ncm_dtuple3_free)

/**
 * ncm_dtuple2_new:
 * @x: (in): first value
 * @y: (in): second value
 *
 * Creates a new #NcmDTuple2.
 *
 * Returns: (transfer full): a new #NcmDTuple2.
 */
NcmDTuple2 *
ncm_dtuple2_new (const gdouble x, const gdouble y)
{
  NcmDTuple2 *tuple = g_slice_new (NcmDTuple2);

  tuple->elements[0] = x;
  tuple->elements[1] = y;

  return tuple;
}

/**
 * ncm_dtuple2_copy:
 * @tuple: (in): a #NcmDTuple2
 *
 * Creates a new #NcmDTuple2 with the same values of @tuple.
 *
 * Returns: (transfer full): a new #NcmDTuple2.
 */
NcmDTuple2 *
ncm_dtuple2_copy (const NcmDTuple2 *tuple)
{
  NcmDTuple2 *new_tuple = g_slice_new (NcmDTuple2);

  new_tuple->elements[0] = tuple->elements[0];
  new_tuple->elements[1] = tuple->elements[1];

  return new_tuple;
}

/**
 * ncm_dtuple2_free:
 * @tuple: (in): a #NcmDTuple2
 *
 * Frees a #NcmDTuple2.
 */
void
ncm_dtuple2_free (NcmDTuple2 *tuple)
{
  g_slice_free (NcmDTuple2, tuple);
}

/**
 * ncm_dtuple3_new:
 * @x: (in): first value
 * @y: (in): second value
 * @z: (in): third value
 *
 * Creates a new #NcmDTuple3.
 *
 * Returns: (transfer full): a new #NcmDTuple3.
 */

NcmDTuple3 *
ncm_dtuple3_new (const gdouble x, const gdouble y, const gdouble z)
{
  NcmDTuple3 *tuple = g_slice_new (NcmDTuple3);

  tuple->elements[0] = x;
  tuple->elements[1] = y;
  tuple->elements[2] = z;

  return tuple;
}

/**
 * ncm_dtuple3_copy:
 * @tuple: (in): a #NcmDTuple3
 *
 * Creates a new #NcmDTuple3 with the same values of @tuple.
 *
 * Returns: (transfer full): a new #NcmDTuple3.
 */
NcmDTuple3 *
ncm_dtuple3_copy (const NcmDTuple3 *tuple)
{
  NcmDTuple3 *new_tuple = g_slice_new (NcmDTuple3);

  new_tuple->elements[0] = tuple->elements[0];
  new_tuple->elements[1] = tuple->elements[1];
  new_tuple->elements[2] = tuple->elements[2];

  return new_tuple;
}

/**
 * ncm_dtuple3_free:
 * @tuple: (in): a #NcmDTuple3
 *
 * Frees a #NcmDTuple3.
 */
void
ncm_dtuple3_free (NcmDTuple3 *tuple)
{
  g_slice_free (NcmDTuple3, tuple);
}

