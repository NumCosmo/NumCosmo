/***************************************************************************
 *            ncm_quaternion.h
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

#ifndef _NCM_QUATERNION_H
#define _NCM_QUATERNION_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_rng.h>

G_BEGIN_DECLS

typedef struct _NcmTriVec NcmTriVec;
typedef struct _NcmQuaternion NcmQuaternion;

struct _NcmTriVec
{
  /*< private >*/
  gdouble c[3];
};

struct _NcmQuaternion
{
  /*< private >*/
  gdouble s;
  NcmTriVec v;
};

GType ncm_trivec_get_type (void) G_GNUC_CONST;
GType ncm_quaternion_get_type (void) G_GNUC_CONST;

NcmTriVec *ncm_trivec_new (void);
NcmTriVec *ncm_trivec_new_full (const gdouble c[3]);
NcmTriVec *ncm_trivec_new_full_c (const gdouble x, const gdouble y, const gdouble z);

NcmTriVec *ncm_trivec_dup (NcmTriVec *v);
void ncm_trivec_free (NcmTriVec *v);

void ncm_trivec_memcpy (NcmTriVec *dest, const NcmTriVec *orig);
void ncm_trivec_set_0 (NcmTriVec *v);
void ncm_trivec_scale (NcmTriVec *v, const gdouble scale);

gdouble ncm_trivec_norm (NcmTriVec *v);
gdouble ncm_trivec_dot (const NcmTriVec *v1, const NcmTriVec *v2);

void ncm_trivec_normalize (NcmTriVec *v);

gdouble ncm_trivec_get_phi (NcmTriVec *v);
void ncm_trivec_set_spherical_coord (NcmTriVec *v, gdouble r, gdouble theta, gdouble phi);
void ncm_trivec_get_spherical_coord (NcmTriVec *v, gdouble *theta, gdouble *phi);

/* Quaternions */

NcmQuaternion *ncm_quaternion_new (void);
NcmQuaternion *ncm_quaternion_new_from_vector (NcmTriVec *v);
NcmQuaternion *ncm_quaternion_new_from_data (gdouble x, gdouble y, gdouble z, gdouble theta);

NcmQuaternion *ncm_quaternion_dup (NcmQuaternion *q);
void ncm_quaternion_free (NcmQuaternion *q);

void ncm_quaternion_memcpy (NcmQuaternion *dest, const NcmQuaternion *orig);

void ncm_quaternion_set_from_data (NcmQuaternion *q, gdouble x, gdouble y, gdouble z, gdouble theta);
void ncm_quaternion_set_random (NcmQuaternion *q, NcmRNG *rng);

void ncm_quaternion_set_I (NcmQuaternion *q);
void ncm_quaternion_set_0 (NcmQuaternion *q);

gdouble ncm_quaternion_norm (NcmQuaternion *q);

void ncm_quaternion_normalize (NcmQuaternion *q);
void ncm_quaternion_conjugate (NcmQuaternion *q);

void ncm_quaternion_mul (NcmQuaternion *q, NcmQuaternion *u, NcmQuaternion *res);
void ncm_quaternion_lmul (NcmQuaternion *q, NcmQuaternion *u);
void ncm_quaternion_rmul (NcmQuaternion *q, NcmQuaternion *u);

void ncm_quaternion_conjugate_u_mul (NcmQuaternion *q, NcmQuaternion *u, NcmQuaternion *res);
void ncm_quaternion_conjugate_q_mul (NcmQuaternion *q, NcmQuaternion *u, NcmQuaternion *res);

void ncm_quaternion_rotate (NcmQuaternion *q, NcmTriVec *v);
void ncm_quaternion_inv_rotate (NcmQuaternion *q, NcmTriVec *v);

/**
 * NCM_QUATERNION_INIT:
 *
 */
#define NCM_QUATERNION_INIT {0.0, {{0.0, 0.0, 0.0}}}

/**
 * NCM_QUATERNION_INIT_I:
 *
 */
#define NCM_QUATERNION_INIT_I {1.0, {{0.0, 0.0, 0.0}}}

G_END_DECLS

#endif /* _NCM_QUATERNION_H */
