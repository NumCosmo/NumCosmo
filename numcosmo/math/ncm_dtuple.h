/***************************************************************************
 *            ncm_dtuple.h
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


#ifndef _NCM_DTUPLE_H_
#define _NCM_DTUPLE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NCM_TYPE_DTUPLE2 (ncm_dtuple2_get_type ())
#define NCM_TYPE_DTUPLE3 (ncm_dtuple3_get_type ())

GType ncm_dtuple2_get_type (void) G_GNUC_CONST;
GType ncm_dtuple3_get_type (void) G_GNUC_CONST;

/**
 * NcmDTuple2:
 * @elements: (array fixed-size=2): The elements of the tuple.
 *
 * A 2-dimensional tuple of double precision floating point numbers.
 */
typedef struct _NcmDTuple2 NcmDTuple2;

/**
 * NcmDTuple3:
 * @elements: (array fixed-size=3): The elements of the tuple.
 *
 * A 3-dimensional tuple of double precision floating point numbers.
 */
typedef struct _NcmDTuple3 NcmDTuple3;

struct _NcmDTuple2
{
  gdouble elements[2];
};

struct _NcmDTuple3
{
  gdouble elements[3];
};

NcmDTuple2 *ncm_dtuple2_new (const gdouble x, const gdouble y);
NcmDTuple2 *ncm_dtuple2_new_from_variant (GVariant *var);

NcmDTuple3 *ncm_dtuple3_new (const gdouble x, const gdouble y, const gdouble z);
NcmDTuple3 *ncm_dtuple3_new_from_variant (GVariant *var);

NcmDTuple2 *ncm_dtuple2_copy (const NcmDTuple2 *dt2);
NcmDTuple3 *ncm_dtuple3_copy (const NcmDTuple3 *dt3);

GVariant *ncm_dtuple2_serialize (const NcmDTuple2 *dt2);
GVariant *ncm_dtuple3_serialize (const NcmDTuple3 *dt3);

void ncm_dtuple2_free (NcmDTuple2 *dt2);
void ncm_dtuple3_free (NcmDTuple3 *dt3);

void ncm_dtuple2_clear (NcmDTuple2 **dt2);
void ncm_dtuple3_clear (NcmDTuple3 **dt3);

/**
 * NCM_DTUPLE2_STATIC_INIT:
 * @x: The first element of the tuple
 * @y: The second element of the tuple
 *
 * Initializes a #NcmDTuple2 with the given elements.
 *
 * For example:
 * |[<!-- language="C" -->
 * NcmDTuple2 tuple = NCM_DTUPLE2_STATIC_INIT (1.0, 2.0);
 * ]|
 *
 */
#define NCM_DTUPLE2_STATIC_INIT(x, y) \
        {                             \
          { x, y }                    \
        }


/**
 * NCM_DTUPLE3_STATIC_INIT:
 * @x: The first element of the tuple
 * @y: The second element of the tuple
 * @z: The third element of the tuple
 *
 * Initializes a #NcmDTuple3 with the given elements.
 *
 * For example:
 * |[<!-- language="C" -->
 * NcmDTuple3 tuple = NCM_DTUPLE3_STATIC_INIT (1.0, 2.0, 3.0);
 * ]|
 *
 */
#define NCM_DTUPLE3_STATIC_INIT(x, y, z) \
        {                                \
          { x, y, z }                    \
        }

/**
 * NCM_DTUPLE2_TYPE:
 *
 * GVariant type string for #NcmDTuple2.
 */
#define NCM_DTUPLE2_TYPE "(dd)"

/**
 * NCM_DTUPLE3_TYPE:
 *
 * GVariant type string for #NcmDTuple3.
 */
#define NCM_DTUPLE3_TYPE "(ddd)"

G_END_DECLS

#endif /* _NCM_DTUPLE_H_ */

