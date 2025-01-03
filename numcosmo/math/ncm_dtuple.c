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
 * NcmDTuple2:
 * Fixed sized array of double values.
 *
 * In this module we define the #NcmDTuple2 which are fixed sized arrays of double
 * values.
 *
 */

/**
 * NcmDTuple3:
 * Fixed sized array of double values.
 *
 * In this module we define the #NcmDTuple3 which are fixed sized arrays of double
 * values.
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
  NcmDTuple2 *dt2 = g_slice_new (NcmDTuple2);

  dt2->elements[0] = x;
  dt2->elements[1] = y;

  return dt2;
}

/**
 * ncm_dtuple2_new_from_variant:
 * @var: a #GVariant
 *
 * Creates a new #NcmDTuple2 from a #GVariant.
 * The #GVariant must be of type #NCM_DTUPLE2_TYPE.
 *
 */
NcmDTuple2 *
ncm_dtuple2_new_from_variant (GVariant *var)
{
  NcmDTuple2 *dt2 = g_slice_new (NcmDTuple2);

  g_assert (g_variant_is_of_type (var, G_VARIANT_TYPE (NCM_DTUPLE2_TYPE)));

  g_variant_get (var, NCM_DTUPLE2_TYPE, &dt2->elements[0], &dt2->elements[1]);

  return dt2;
}

/**
 * ncm_dtuple2_copy:
 * @dt2: (in): a #NcmDTuple2
 *
 * Creates a new #NcmDTuple2 with the same values of @dt2.
 *
 * Returns: (transfer full): a new #NcmDTuple2.
 */
NcmDTuple2 *
ncm_dtuple2_copy (const NcmDTuple2 *dt2)
{
  NcmDTuple2 *new_dt2 = g_slice_new (NcmDTuple2);

  new_dt2->elements[0] = dt2->elements[0];
  new_dt2->elements[1] = dt2->elements[1];

  return new_dt2;
}

/**
 * ncm_dtuple2_serialize:
 * @dt2: a #NcmDTuple2
 *
 * Serializes a #NcmDTuple2.
 *
 * Returns: (transfer full): a #GVariant.
 */
GVariant *
ncm_dtuple2_serialize (const NcmDTuple2 *dt2)
{
  return g_variant_ref_sink (g_variant_new (NCM_DTUPLE2_TYPE, dt2->elements[0], dt2->elements[1]));
}

/**
 * ncm_dtuple2_free:
 * @dt2: (in): a #NcmDTuple2
 *
 * Frees a #NcmDTuple2.
 */
void
ncm_dtuple2_free (NcmDTuple2 *dt2)
{
  g_slice_free (NcmDTuple2, dt2);
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
  NcmDTuple3 *dt3 = g_slice_new (NcmDTuple3);

  dt3->elements[0] = x;
  dt3->elements[1] = y;
  dt3->elements[2] = z;

  return dt3;
}

/**
 * ncm_dtuple3_new_from_variant:
 * @var: a #GVariant
 *
 * Creates a new #NcmDTuple3 from a #GVariant.
 *
 * The #GVariant must be of type #NCM_DTUPLE3_TYPE.
 *
 * Returns: (transfer full): a new #NcmDTuple3.
 */
NcmDTuple3 *
ncm_dtuple3_new_from_variant (GVariant *var)
{
  NcmDTuple3 *dt3 = g_slice_new (NcmDTuple3);

  g_assert (g_variant_is_of_type (var, G_VARIANT_TYPE (NCM_DTUPLE3_TYPE)));

  g_variant_get (var, NCM_DTUPLE3_TYPE, &dt3->elements[0], &dt3->elements[1], &dt3->elements[2]);

  return dt3;
}

/**
 * ncm_dtuple3_copy:
 * @dt3: (in): a #NcmDTuple3
 *
 * Creates a new #NcmDTuple3 with the same values of @dt3.
 *
 * Returns: (transfer full): a new #NcmDTuple3.
 */
NcmDTuple3 *
ncm_dtuple3_copy (const NcmDTuple3 *dt3)
{
  NcmDTuple3 *new_dt3 = g_slice_new (NcmDTuple3);

  new_dt3->elements[0] = dt3->elements[0];
  new_dt3->elements[1] = dt3->elements[1];
  new_dt3->elements[2] = dt3->elements[2];

  return new_dt3;
}

/**
 * ncm_dtuple3_serialize:
 * @dt3: a #NcmDTuple3
 *
 * Serializes a #NcmDTuple3.
 *
 * Returns: (transfer full): a #GVariant.
 */
GVariant *
ncm_dtuple3_serialize (const NcmDTuple3 *dt3)
{
  return g_variant_ref_sink (g_variant_new (NCM_DTUPLE3_TYPE, dt3->elements[0], dt3->elements[1], dt3->elements[2]));
}

/**
 * ncm_dtuple3_free:
 * @dt3: (in): a #NcmDTuple3
 *
 * Frees a #NcmDTuple3.
 */
void
ncm_dtuple3_free (NcmDTuple3 *dt3)
{
  g_slice_free (NcmDTuple3, dt3);
}

/**
 * ncm_dtuple2_clear:
 * @dt2: a #NcmDTuple2
 *
 * If *@dt2 is not NULL, frees it and sets it to NULL.
 *
 */
void
ncm_dtuple2_clear (NcmDTuple2 **dt2)
{
  if (*dt2 != NULL)
  {
    ncm_dtuple2_free (*dt2);
    *dt2 = NULL;
  }
}

/**
 * ncm_dtuple3_clear:
 * @dt3: a #NcmDTuple3
 *
 * If *@dt3 is not NULL, frees it and sets it to NULL.
 *
 */
void
ncm_dtuple3_clear (NcmDTuple3 **dt3)
{
  if (*dt3 != NULL)
  {
    ncm_dtuple3_free (*dt3);
    *dt3 = NULL;
  }
}

