/***************************************************************************
 *            quaternion.h
 *
 *  Fri Aug 22 16:40:54 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NC_QUATERNION_H
#define _NC_QUATERNION_H

#include <glib.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

typedef struct _NcmTriVector NcmTriVector;


/**
 * NcmTriVector:
 *
 * FIXME
 */
struct _NcmTriVector
{
  /*< private >*/
  gdouble c[3];
};

#define NC_TRIVEC_SCALE(a, b) \
do { \
const gdouble __in_value__ = (b); \
(a).c[0] *= __in_value__; \
(a).c[1] *= __in_value__; \
(a).c[2] *= __in_value__; \
} while (FALSE)

#define NC_TRIVEC_NEW {{0.0, 0.0, 0.0}}
#define NC_TRIVEC_SET_0(v) memset(&(v), 0, sizeof(NcmTriVector))
#define NC_TRIVEC_MEMCPY(a, b) memcpy (&(a), &(b), sizeof (NcmTriVector))
#define NC_TRIVEC_NORM(a) sqrt((a).c[0]*(a).c[0] + (a).c[1]*(a).c[1] + (a).c[2]*(a).c[2])
#define NC_TRIVEC_NORMALIZE(a) NC_TRIVEC_SCALE(a, 1.0/NC_TRIVEC_NORM(a))
#define NC_TRIVEC_DOT(a, b) ((a).c[0]*(b).c[0] + (a).c[1]*(b).c[1] + (a).c[2]*(b).c[2])

typedef struct _NcmQ NcmQ;

/**
 * NcmQ:
 *
 * FIXME
 */
struct _NcmQ
{
  /*< private >*/
  gdouble s;
  NcmTriVector x;
};

/**
 * NC_QUATERNION_NEW:
 *
 */
#define NC_QUATERNION_NEW {0.0, {{0.0, 0.0, 0.0}}}

/**
 * NC_QUATERNION_NEW_I:
 *
 */
#define NC_QUATERNION_NEW_I {1.0, {{0.0, 0.0, 0.0}}}

/**
 * NC_QUATERNION_SET_I:
 * @q: FIXME
 *
 * FIXME
 *
 */
#define NC_QUATERNION_SET_I(q) \
do { \
(q)->s = 1.0; \
NC_TRIVEC_SET_0(q->x); \
} while(FALSE)

/**
 * NC_QUATERNION_SET_0:
 * @q: FIXME
 *
 * FIXME
 *
 */
#define NC_QUATERNION_SET_0(q) \
do { \
(q)->s = 0.0; \
NC_TRIVEC_SET_0(q->x); \
} while(FALSE)

/**
 * NC_QUATERNION_NORM:
 * @q: FIXME
 *
 * FIXME
 *
 */
#define NC_QUATERNION_NORM(q) \
(sqrtf(q->s*q->s + q->x.c[0]*q->x.c[0] + q->x.c[1]*q->x.c[1] + q->x.c[2]*q->x.c[2]))

/**
 * NC_QUATERNION_MEMCPY:
 * @a: FIXME
 * @b: FIXME
 *
 * FIXME
 *
 */
#define NC_QUATERNION_MEMCPY(a,b) memcpy (a, b, sizeof(NcmQ))

NcmQ *ncm_quaternion_new ();
NcmQ *ncm_quaternion_new_from_vector (NcmTriVector v);
NcmQ *ncm_quaternion_new_from_data (gdouble x, gdouble y, gdouble z, gdouble theta);

void ncm_quaternion_set_from_data (NcmQ *q, gdouble x, gdouble y, gdouble z, gdouble theta);
void ncm_quaternion_set_random (NcmQ *q);

void ncm_quaternion_free (NcmQ *q);
void ncm_quaternion_normalize (NcmQ *q);
void ncm_quaternion_conjugate (NcmQ *q);

void ncm_quaternion_mul (NcmQ *q, NcmQ *u, NcmQ *res);
void ncm_quaternion_lmul (NcmQ *q, NcmQ *u);
void ncm_quaternion_rmul (NcmQ *q, NcmQ *u);

void ncm_quaternion_conjugate_u_mul (NcmQ *q, NcmQ *u, NcmQ *res);
void ncm_quaternion_conjugate_q_mul (NcmQ *q, NcmQ *u, NcmQ *res);

void ncm_quaternion_rotate (NcmQ *q, NcmTriVector v);
void ncm_quaternion_inv_rotate (NcmQ *q, NcmTriVector v);

#define NCM_QUATERNION_RNG_NAME "quaternion"

G_END_DECLS

#endif /* _NC_QUATERNION_H */
