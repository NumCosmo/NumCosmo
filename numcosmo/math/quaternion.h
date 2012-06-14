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

G_BEGIN_DECLS

typedef struct _NcTriVector NcTriVector;


/**
 * NcTriVector:
 *
 * FIXME
 */
struct _NcTriVector
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
#define NC_TRIVEC_SET_0(v) memset(&(v), 0, sizeof(NcTriVector))
#define NC_TRIVEC_MEMCPY(a, b) memcpy (&(a), &(b), sizeof (NcTriVector))
#define NC_TRIVEC_NORM(a) sqrt((a).c[0]*(a).c[0] + (a).c[1]*(a).c[1] + (a).c[2]*(a).c[2])
#define NC_TRIVEC_NORMALIZE(a) NC_TRIVEC_SCALE(a, 1.0/NC_TRIVEC_NORM(a))
#define NC_TRIVEC_DOT(a, b) ((a).c[0]*(b).c[0] + (a).c[1]*(b).c[1] + (a).c[2]*(b).c[2])

typedef struct _NcQ NcQ;

/**
 * NcQ:
 * 
 * FIXME
 */
struct _NcQ
{
  /*< private >*/
  gdouble s;
  NcTriVector x;
};

/**
 * FIXME
 */
#define NC_QUATERNION_NEW {0.0, {{0.0, 0.0, 0.0}}}

/**
 * FIXME
 */
#define NC_QUATERNION_NEW_I {1.0, {{0.0, 0.0, 0.0}}}

/**
 * FIXME
 */
#define NC_QUATERNION_SET_I(q) \
do { \
(q)->s = 1.0; \
NC_TRIVEC_SET_0(q->x); \
} while(FALSE)

/**
 * FIXME
 */
#define NC_QUATERNION_SET_0(q) \
do { \
(q)->s = 0.0; \
NC_TRIVEC_SET_0(q->x); \
} while(FALSE)

/**
 * FIXME
 */
#define NC_QUATERNION_NORM(q) \
(sqrtf(q->s*q->s + q->x.c[0]*q->x.c[0] + q->x.c[1]*q->x.c[1] + q->x.c[2]*q->x.c[2]))

/**
 * FIXME
 */
#define NC_QUATERNION_MEMCPY(a,b) memcpy (a, b, sizeof(NcQ))

NcQ *nc_quaternion_new ();
NcQ *nc_quaternion_new_from_vector (NcTriVector v);
NcQ *nc_quaternion_new_from_data (gdouble x, gdouble y, gdouble z, gdouble theta);

void nc_quaternion_set_from_data (NcQ *q, gdouble x, gdouble y, gdouble z, gdouble theta);
void nc_quaternion_set_random (NcQ *q);

void nc_quaternion_free (NcQ *q);
void nc_quaternion_normalize (NcQ *q);
void nc_quaternion_conjugate (NcQ *q);

void nc_quaternion_mul (NcQ *q, NcQ *u, NcQ *res);
void nc_quaternion_lmul (NcQ *q, NcQ *u);
void nc_quaternion_rmul (NcQ *q, NcQ *u);

void nc_quaternion_conjugate_u_mul (NcQ *q, NcQ *u, NcQ *res);
void nc_quaternion_conjugate_q_mul (NcQ *q, NcQ *u, NcQ *res);

void nc_quaternion_rotate (NcQ *q, NcTriVector v);
void nc_quaternion_inv_rotate (NcQ *q, NcTriVector v);

G_END_DECLS

#endif /* _NC_QUATERNION_H */
