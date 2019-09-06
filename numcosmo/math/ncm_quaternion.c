/***************************************************************************
 *            ncm_quaternion.c
 *
 *  Fri Aug 22 16:40:29 2008
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

/**
 * SECTION:ncm_quaternion
 * @title: NcmQuaternion
 * @short_description: Quaternions algebra, three-vectors and mapping to matrix.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_quaternion.h"
#include "math/ncm_rng.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <string.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_BOXED_TYPE (NcmQuaternion, ncm_quaternion, ncm_quaternion_dup, ncm_quaternion_free);
G_DEFINE_BOXED_TYPE (NcmTriVec, ncm_trivec, ncm_trivec_dup, ncm_trivec_free);

/**
 * ncm_trivec_new: (constructor)
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmTriVec *
ncm_trivec_new (void)
{
  NcmTriVec *v = g_new0 (NcmTriVec, 1);
  return v;
}

/**
 * ncm_trivec_new_full: (constructor)
 * @c: components
 * 
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmTriVec *
ncm_trivec_new_full (const gdouble c[3])
{
  NcmTriVec *v = g_new0 (NcmTriVec, 1);
  memcpy (v->c, c, sizeof (gdouble) * 3);
  return v;
}

/**
 * ncm_trivec_new_full_c: (constructor)
 * @x: x-component
 * @y: y-component
 * @z: z-component
 * 
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmTriVec *
ncm_trivec_new_full_c (const gdouble x, const gdouble y, const gdouble z)
{
  NcmTriVec *v = g_new0 (NcmTriVec, 1);

  v->c[0] = x;
  v->c[1] = y;
  v->c[2] = z;
  
  return v;
}

/**
 * ncm_trivec_dup:
 * @v: a #NcmTriVec
 * 
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmTriVec *
ncm_trivec_dup (NcmTriVec *v)
{
  NcmTriVec *cv = ncm_trivec_new ();
  ncm_trivec_memcpy (cv, v);
  return cv;
}

/**
 * ncm_trivec_free:
 * @v: a #NcmTriVec
 * 
 * FIXME
 *
 */
void 
ncm_trivec_free (NcmTriVec *v)
{
  g_free (v);
}

/**
 * ncm_trivec_memcpy:
 * @dest: a #NcmTriVec
 * @orig: a #NcmTriVec
 * 
 * FIXME
 *
 */
void 
ncm_trivec_memcpy (NcmTriVec *dest, const NcmTriVec *orig)
{
  memcpy (dest, orig, sizeof (NcmTriVec));
}

/**
 * ncm_trivec_set_0:
 * @v: a #NcmTriVec
 * 
 * FIXME
 *
 */
void 
ncm_trivec_set_0 (NcmTriVec *v)
{
  memset (v, 0, sizeof (NcmTriVec));
}

/**
 * ncm_trivec_scale:
 * @v: a #NcmTriVec
 * @scale: FIXME
 * 
 * FIXME
 *
 */
void 
ncm_trivec_scale (NcmTriVec *v, const gdouble scale)
{
  v->c[0] *= scale;
  v->c[1] *= scale;
  v->c[2] *= scale;
}

/**
 * ncm_trivec_norm:
 * @v: a #NcmTriVec
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
gdouble
ncm_trivec_norm (NcmTriVec *v)
{
  return gsl_hypot3 (v->c[0], v->c[1], v->c[2]);
}

/**
 * ncm_trivec_dot:
 * @v1: a #NcmTriVec
 * @v2: a #NcmTriVec
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
gdouble
ncm_trivec_dot (const NcmTriVec *v1, const NcmTriVec *v2)
{
  return v1->c[0] * v2->c[0] + v1->c[1] * v2->c[1] + v1->c[2] * v2->c[2];
}

/**
 * ncm_trivec_normalize:
 * @v: a #NcmTriVec
 * 
 * FIXME
 * 
 */
void
ncm_trivec_normalize (NcmTriVec *v)
{
  ncm_trivec_scale (v, 1.0 / ncm_trivec_norm (v));
}

/**
 * ncm_trivec_get_phi:
 * @v: a #NcmTriVec
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
gdouble 
ncm_trivec_get_phi (NcmTriVec *v)
{
  return atan2 (v->c[1], v->c[0]);
}

/**
 * ncm_trivec_set_spherical_coord:
 * @v: a #NcmTriVec
 * @r: FIXME
 * @theta: FIXME
 * @phi: FIXME
 * 
 * FIXME
 * 
 */
void 
ncm_trivec_set_spherical_coord (NcmTriVec *v, gdouble r, gdouble theta, gdouble phi)
{
  v->c[0] = r * sin (theta) * cos (phi);
  v->c[1] = r * sin (theta) * sin (phi);
  v->c[2] = r * cos (theta);
}

/**
 * ncm_trivec_get_spherical_coord:
 * @v: a #NcmTriVec
 * @theta: (out): FIXME
 * @phi: (out): FIXME
 * 
 * FIXME
 * 
 */
void 
ncm_trivec_get_spherical_coord (NcmTriVec *v, gdouble *theta, gdouble *phi)
{
  const gdouble norm = ncm_trivec_norm (v);
  theta[0] = acos (v->c[2] / norm);
  phi[0]   = ncm_trivec_get_phi (v);
}

/**
 * ncm_quaternion_new: (constructor)
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmQuaternion *
ncm_quaternion_new (void)
{
  NcmQuaternion *q = (NcmQuaternion *) g_new0 (NcmQuaternion, 1);
  return q;
}

/**
 * ncm_quaternion_new_from_vector: (constructor)
 * @v: a #NcmTriVec
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
*/
NcmQuaternion *
ncm_quaternion_new_from_vector (NcmTriVec *v)
{
  NcmQuaternion *q = ncm_quaternion_new ();
  ncm_trivec_memcpy (&q->v, v);
  q->s = 0.0;
  return q;
}

/**
 * ncm_quaternion_new_from_data: (constructor)
 * @x: FIXME
 * @y: FIXME
 * @z: FIXME
 * @theta: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmQuaternion *
ncm_quaternion_new_from_data (gdouble x, gdouble y, gdouble z, gdouble theta)
{
  NcmQuaternion *q = ncm_quaternion_new ();
  
  theta /= 2.0;
  q->v.c[0] = x;
  q->v.c[1] = y;
  q->v.c[2] = z;

  ncm_trivec_normalize (&q->v);
  ncm_trivec_scale (&q->v, sin (theta));
  q->s = cos (theta);
  
  return q;
}

/**
 * ncm_quaternion_dup:
 * @q: a #NcmQuaternion
 *
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
NcmQuaternion *
ncm_quaternion_dup (NcmQuaternion *q)
{
  NcmQuaternion *cq = ncm_quaternion_new ();
  ncm_quaternion_memcpy (cq, q);
  
  return cq;
}

/**
 * ncm_quaternion_free:
 * @q: a #NcmQuaternion
 *
 * FIXME
 * 
 */
void
ncm_quaternion_free (NcmQuaternion *q)
{
  g_free (q);
}

/**
 * ncm_quaternion_memcpy:
 * @dest: a #NcmQuaternion
 * @orig: a #NcmQuaternion
 *
 * FIXME
 */
void
ncm_quaternion_memcpy (NcmQuaternion *dest, const NcmQuaternion *orig)
{
  memcpy (dest, orig, sizeof (NcmQuaternion));
}

/**
 * ncm_quaternion_set_from_data:
 * @q: FIXME
 * @x: FIXME
 * @y: FIXME
 * @z: FIXME
 * @theta: FIXME
 *
 * FIXME
 *
 */
void
ncm_quaternion_set_from_data (NcmQuaternion *q, gdouble x, gdouble y, gdouble z, gdouble theta)
{
  theta /= 2.0;
  q->v.c[0] = x;
  q->v.c[1] = y;
  q->v.c[2] = z;

  ncm_trivec_normalize (&q->v);
  ncm_trivec_scale (&q->v, sin (theta));
  
  q->s = cos(theta);
}

/**
 * ncm_quaternion_set_I:
 * @q: FIXME
 *
 * FIXME
 *
 */
void 
ncm_quaternion_set_I (NcmQuaternion *q)
{
  q->s = 1.0;
  ncm_trivec_set_0 (&q->v);
}

/**
 * ncm_quaternion_set_0:
 * @q: FIXME
 *
 * FIXME
 *
 */
void 
ncm_quaternion_set_0 (NcmQuaternion *q)
{
  q->s = 0.0;
  ncm_trivec_set_0 (&q->v);
}

/**
 * ncm_quaternion_norm:
 * @q: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble 
ncm_quaternion_norm (NcmQuaternion *q)
{
  return sqrt (q->s * q->s + q->v.c[0] * q->v.c[0] + q->v.c[1] * q->v.c[1] + q->v.c[2] * q->v.c[2]); 
}

/**
 * ncm_quaternion_set_random:
 * @q: FIXME
 * @rng: a #NcmRNG
 *
 * FIXME
 *
 */
void
ncm_quaternion_set_random (NcmQuaternion *q, NcmRNG *rng)
{
  ncm_rng_lock (rng);
  ncm_quaternion_set_from_data (q,
                                 -1.0 + 2.0 * gsl_rng_uniform_pos (rng->r),
                                 -1.0 + 2.0 * gsl_rng_uniform_pos (rng->r),
                                 -1.0 + 2.0 * gsl_rng_uniform_pos (rng->r),
                                 2.0 * M_PI * gsl_rng_uniform (rng->r));
  ncm_rng_unlock (rng);
  ncm_rng_free (rng);   
}

/**
 * ncm_quaternion_normalize:
 * @q: FIXME
 *
 * FIXME
 *
 */
void
ncm_quaternion_normalize (NcmQuaternion *q)
{
  const gdouble norm = ncm_quaternion_norm (q);

  q->s /= norm;
  q->v.c[0] /= norm;
  q->v.c[1] /= norm;
  q->v.c[2] /= norm;
}

/**
 * ncm_quaternion_conjugate:
 * @q: FIXME
 *
 * FIXME
 *
 */
void
ncm_quaternion_conjugate (NcmQuaternion *q)
{
  ncm_trivec_scale (&q->v, -1.0);
}

/**
 * ncm_quaternion_mul:
 * @q: FIXME
 * @u: FIXME
 * @res: FIXME
 *
 * FIXME
 *
 */
void
ncm_quaternion_mul (NcmQuaternion *q, NcmQuaternion *u, NcmQuaternion *res)
{
  res->s = q->s * u->s - ncm_trivec_dot (&q->v, &u->v);
  res->v.c[0] = q->s * u->v.c[0] + u->s * q->v.c[0] + q->v.c[1] * u->v.c[2] - q->v.c[2] * u->v.c[1];
  res->v.c[1] = q->s * u->v.c[1] + u->s * q->v.c[1] + q->v.c[2] * u->v.c[0] - q->v.c[0] * u->v.c[2];
  res->v.c[2] = q->s * u->v.c[2] + u->s * q->v.c[2] + q->v.c[0] * u->v.c[1] - q->v.c[1] * u->v.c[0];
}

/**
 * ncm_quaternion_lmul:
 * @q: FIXME
 * @u: FIXME
 *
 * FIXME
 *
 */
void
ncm_quaternion_lmul (NcmQuaternion *q, NcmQuaternion *u)
{
  NcmQuaternion t;

  ncm_quaternion_mul (u, q, &t);
  ncm_quaternion_memcpy (q, &t);
}

/**
 * ncm_quaternion_rmul:
 * @q: FIXME
 * @u: FIXME
 *
 * FIXME
 *
 */
void
ncm_quaternion_rmul (NcmQuaternion *q, NcmQuaternion *u)
{
  NcmQuaternion t;

  ncm_quaternion_mul (q, u, &t);
  ncm_quaternion_memcpy (q, &t);
}

/**
 * ncm_quaternion_conjugate_q_mul:
 * @q: FIXME
 * @u: FIXME
 * @res: FIXME
 *
 * FIXME
 *
 */
void
ncm_quaternion_conjugate_q_mul (NcmQuaternion *q, NcmQuaternion *u, NcmQuaternion *res)
{
  res->s = q->s * u->s + ncm_trivec_dot (&q->v, &u->v);
  res->v.c[0] = q->s * u->v.c[0] - u->s * q->v.c[0] - q->v.c[1] * u->v.c[2] + q->v.c[2] * u->v.c[1];
  res->v.c[1] = q->s * u->v.c[1] - u->s * q->v.c[1] - q->v.c[2] * u->v.c[0] + q->v.c[0] * u->v.c[2];
  res->v.c[2] = q->s * u->v.c[2] - u->s * q->v.c[2] - q->v.c[0] * u->v.c[1] + q->v.c[1] * u->v.c[0];
}

/**
 * ncm_quaternion_conjugate_u_mul:
 * @q: FIXME
 * @u: FIXME
 * @res: FIXME
 *
 * FIXME
 *
 */
void
ncm_quaternion_conjugate_u_mul (NcmQuaternion *q, NcmQuaternion *u, NcmQuaternion *res)
{
  res->s = q->s * u->s + ncm_trivec_dot (&q->v, &u->v);
  res->v.c[0] = - q->s * u->v.c[0] + u->s * q->v.c[0] - q->v.c[1] * u->v.c[2] + q->v.c[2] * u->v.c[1];
  res->v.c[1] = - q->s * u->v.c[1] + u->s * q->v.c[1] - q->v.c[2] * u->v.c[0] + q->v.c[0] * u->v.c[2];
  res->v.c[2] = - q->s * u->v.c[2] + u->s * q->v.c[2] - q->v.c[0] * u->v.c[1] + q->v.c[1] * u->v.c[0];
}

/**
 * ncm_quaternion_rotate:
 * @q: FIXME
 * @v: FIXME
 *
 * FIXME
 *
 */
void
ncm_quaternion_rotate (NcmQuaternion *q, NcmTriVec *v)
{
  NcmQuaternion qv = NCM_QUATERNION_INIT;
  NcmQuaternion t;

  ncm_trivec_memcpy (&qv.v, v);
  
  ncm_quaternion_mul (q, &qv, &t);
  ncm_quaternion_conjugate_u_mul (&t, q, &qv);

  ncm_trivec_memcpy (v, &qv.v);
}

/**
 * ncm_quaternion_inv_rotate:
 * @q: FIXME
 * @v: FIXME
 *
 * FIXME
 *
 */
void
ncm_quaternion_inv_rotate (NcmQuaternion *q, NcmTriVec *v)
{
  NcmQuaternion qv = NCM_QUATERNION_INIT;
  NcmQuaternion t;
  
  ncm_trivec_memcpy (&qv.v, v);

  ncm_quaternion_mul (&qv, q, &t);
  ncm_quaternion_conjugate_q_mul (q, &t, &qv);

  ncm_trivec_memcpy (v, &qv.v);
}
