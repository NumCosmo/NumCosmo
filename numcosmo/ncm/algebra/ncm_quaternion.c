/***************************************************************************
 *            ncm_quaternion.c
 *
 *  Fri Aug 22 16:40:29 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
 * NcmQuaternion:
 *
 * Quaternions and three-vectors, for rotations in three dimensions.
 *
 * A quaternion $q = s + \vec v$ has a scalar part $s$ and a vector part $\vec v$; a
 * three-vector, #NcmTriVec, is a quaternion with $s = 0$. The conjugate is
 * $q^\dagger = s - \vec v$, and the norm is $|q| = \sqrt{s^2 + \vec v\cdot\vec v}$. A unit
 * quaternion $q = \cos(\theta/2) + \sin(\theta/2)\,\hat n$ rotates a vector by the angle
 * $\theta$ about $\hat n$, counterclockwise seen from the tip of $\hat n$, as
 * $\vec u \to q\,\vec u\,q^\dagger$.
 *
 * The spherical coordinates of a vector are
 * $(r\sin\theta\cos\phi, r\sin\theta\sin\phi, r\cos\theta)$, with the polar angle $\theta$ from
 * the z-axis and the azimuth $\phi$ from the x-axis. The astronomical ones are
 * $(r\cos\delta\cos\alpha, r\cos\delta\sin\alpha, r\sin\delta)$, with the declination $\delta$
 * from the xy-plane (the equator) and the right ascension $\alpha$ from the x-axis.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_quaternion.h"
#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_rng.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <string.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

static void _ncm_quaternion_set_zy (NcmQuaternion *q, NcmTriVec *v, const gdouble beta);

G_DEFINE_BOXED_TYPE (NcmQuaternion, ncm_quaternion, ncm_quaternion_dup, ncm_quaternion_free)
G_DEFINE_BOXED_TYPE (NcmTriVec, ncm_trivec, ncm_trivec_dup, ncm_trivec_free)

/**
 * ncm_trivec_new: (constructor)
 *
 * Creates a zero #NcmTriVec.
 *
 * Returns: (transfer full): a new #NcmTriVec.
 */
NcmTriVec *
ncm_trivec_new (void)
{
  NcmTriVec *v = g_new0 (NcmTriVec, 1);

  return v;
}

/**
 * ncm_trivec_new_full: (constructor)
 * @c: (array fixed-size=3) (element-type double): the components
 *
 * Returns: (transfer full): a new #NcmTriVec.
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
 * @x: the x component
 * @y: the y component
 * @z: the z component
 *
 * Returns: (transfer full): a new #NcmTriVec.
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
 * ncm_trivec_new_sphere: (constructor)
 * @r: the radius
 * @theta: the polar angle
 * @phi: the azimuth
 *
 * Creates the vector with the given spherical coordinates, see ncm_trivec_set_spherical_coord().
 *
 * Returns: (transfer full): a new #NcmTriVec.
 */
NcmTriVec *
ncm_trivec_new_sphere (gdouble r, gdouble theta, gdouble phi)
{
  NcmTriVec *v = g_new0 (NcmTriVec, 1);

  ncm_trivec_set_spherical_coord (v, r, theta, phi);

  return v;
}

/**
 * ncm_trivec_new_astro_coord: (constructor)
 * @r: the radius
 * @delta: the declination, in radians
 * @alpha: the right ascension, in radians
 *
 * Creates the vector with the given astronomical coordinates, see ncm_trivec_set_astro_coord().
 *
 * Returns: (transfer full): a new #NcmTriVec.
 */
NcmTriVec *
ncm_trivec_new_astro_coord (gdouble r, gdouble delta, gdouble alpha)
{
  NcmTriVec *v = g_new0 (NcmTriVec, 1);

  ncm_trivec_set_astro_coord (v, r, delta, alpha);

  return v;
}

/**
 * ncm_trivec_new_astro_ra_dec: (constructor)
 * @r: the radius
 * @ra: the right ascension, in degrees
 * @dec: the declination, in degrees
 *
 * Creates the vector with the given astronomical coordinates, see ncm_trivec_set_astro_ra_dec().
 *
 * Returns: (transfer full): a new #NcmTriVec.
 */
NcmTriVec *
ncm_trivec_new_astro_ra_dec (gdouble r, gdouble ra, gdouble dec)
{
  NcmTriVec *v = g_new0 (NcmTriVec, 1);

  ncm_trivec_set_astro_ra_dec (v, r, ra, dec);

  return v;
}

/**
 * ncm_trivec_dup:
 * @v: a #NcmTriVec
 *
 * Returns: (transfer full): a copy of @v.
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
 * Frees @v.
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
 * Copies @orig into @dest.
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
 * Sets @v to zero.
 */
void
ncm_trivec_set_0 (NcmTriVec *v)
{
  memset (v, 0, sizeof (NcmTriVec));
}

/**
 * ncm_trivec_scale:
 * @v: a #NcmTriVec
 * @scale: the factor
 *
 * Multiplies @v by @scale.
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
 * Returns: $|\vec v|$.
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
 * Returns: $\vec v_1\cdot\vec v_2$.
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
 * Divides @v by its norm.
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
 * Returns: the azimuth $\phi \in (-\pi, \pi]$ of @v.
 */
gdouble
ncm_trivec_get_phi (NcmTriVec *v)
{
  return atan2 (v->c[1], v->c[0]);
}

/**
 * ncm_trivec_set_spherical_coord:
 * @v: a #NcmTriVec
 * @r: the radius
 * @theta: the polar angle
 * @phi: the azimuth
 *
 * Sets @v to $(r\sin\theta\cos\phi, r\sin\theta\sin\phi, r\cos\theta)$.
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
 * @r: (out): the radius
 * @theta: (out): the polar angle
 * @phi: (out): the azimuth
 *
 * Computes the spherical coordinates of @v, with $\theta \in [0, \pi]$ and $\phi \in (-\pi, \pi]$.
 */
void
ncm_trivec_get_spherical_coord (NcmTriVec *v, gdouble *r, gdouble *theta, gdouble *phi)
{
  const gdouble norm = ncm_trivec_norm (v);

  /* atan2 keeps full precision near the poles, where acos (z / r) does not */
  r[0]     = norm;
  theta[0] = atan2 (hypot (v->c[0], v->c[1]), v->c[2]);
  phi[0]   = ncm_trivec_get_phi (v);
}

/**
 * ncm_trivec_set_astro_coord:
 * @v: a #NcmTriVec
 * @r: the radius
 * @delta: the declination, in radians
 * @alpha: the right ascension, in radians
 *
 * Sets @v to $(r\cos\delta\cos\alpha, r\cos\delta\sin\alpha, r\sin\delta)$.
 */
void
ncm_trivec_set_astro_coord (NcmTriVec *v, gdouble r, gdouble delta, gdouble alpha)
{
  v->c[0] = r * cos (delta) * cos (alpha);
  v->c[1] = r * cos (delta) * sin (alpha);
  v->c[2] = r * sin (delta);
}

/**
 * ncm_trivec_get_astro_coord:
 * @v: a #NcmTriVec
 * @r: (out): the radius
 * @delta: (out): the declination, in radians
 * @alpha: (out): the right ascension, in radians
 *
 * Computes the astronomical coordinates of @v, with $\delta \in [-\pi/2, \pi/2]$ and
 * $\alpha \in (-\pi, \pi]$.
 */
void
ncm_trivec_get_astro_coord (NcmTriVec *v, gdouble *r, gdouble *delta, gdouble *alpha)
{
  const gdouble norm = ncm_trivec_norm (v);

  /* atan2 keeps full precision near the poles, where asin (z / r) does not */
  r[0]     = norm;
  delta[0] = atan2 (v->c[2], hypot (v->c[0], v->c[1]));
  alpha[0] = atan2 (v->c[1], v->c[0]);
}

/**
 * ncm_trivec_set_astro_ra_dec:
 * @v: a #NcmTriVec
 * @r: the radius
 * @ra: the right ascension, in degrees
 * @dec: the declination, in degrees
 *
 * Same as ncm_trivec_set_astro_coord() with the angles in degrees.
 */
void
ncm_trivec_set_astro_ra_dec (NcmTriVec *v, gdouble r, gdouble ra, gdouble dec)
{
  const gdouble delta = ncm_c_degree_to_radian (dec);
  const gdouble alpha = ncm_c_degree_to_radian (ra);

  ncm_trivec_set_astro_coord (v, r, delta, alpha);
}

/**
 * ncm_trivec_get_astro_ra_dec:
 * @v: a #NcmTriVec
 * @r: (out): the radius
 * @ra: (out): the right ascension, in degrees
 * @dec: (out): the declination, in degrees
 *
 * Same as ncm_trivec_get_astro_coord() with the angles in degrees.
 */
void
ncm_trivec_get_astro_ra_dec (NcmTriVec *v, gdouble *r, gdouble *ra, gdouble *dec)
{
  gdouble delta, alpha;

  ncm_trivec_get_astro_coord (v, r, &delta, &alpha);

  dec[0] = ncm_c_radian_to_degree (delta);
  ra[0]  = ncm_c_radian_to_degree (alpha);
}

/**
 * ncm_quaternion_new: (constructor)
 *
 * Creates a zero #NcmQuaternion.
 *
 * Returns: (transfer full): a new #NcmQuaternion.
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
 * Creates the quaternion $0 + \vec v$.
 *
 * Returns: (transfer full): a new #NcmQuaternion.
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
 * @x: the x component of the axis
 * @y: the y component of the axis
 * @z: the z component of the axis
 * @theta: the rotation angle
 *
 * Creates the rotation by @theta about the axis $(x, y, z)$, see
 * ncm_quaternion_set_from_data().
 *
 * Returns: (transfer full): a new #NcmQuaternion.
 */
NcmQuaternion *
ncm_quaternion_new_from_data (gdouble x, gdouble y, gdouble z, gdouble theta)
{
  NcmQuaternion *q = ncm_quaternion_new ();

  ncm_quaternion_set_from_data (q, x, y, z, theta);

  return q;
}

/**
 * ncm_quaternion_dup:
 * @q: a #NcmQuaternion
 *
 * Returns: (transfer full): a copy of @q.
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
 * Frees @q.
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
 * Copies @orig into @dest.
 */
void
ncm_quaternion_memcpy (NcmQuaternion *dest, const NcmQuaternion *orig)
{
  memcpy (dest, orig, sizeof (NcmQuaternion));
}

/**
 * ncm_quaternion_set_from_data:
 * @q: a #NcmQuaternion
 * @x: the x component of the axis
 * @y: the y component of the axis
 * @z: the z component of the axis
 * @theta: the rotation angle
 *
 * Sets @q to the unit quaternion $\cos(\theta/2) + \sin(\theta/2)\,\hat n$, the rotation by
 * @theta about $\hat n = (x, y, z)/|(x, y, z)|$. Aborts if the axis is zero.
 */
void
ncm_quaternion_set_from_data (NcmQuaternion *q, gdouble x, gdouble y, gdouble z, gdouble theta)
{
  theta    /= 2.0;
  q->v.c[0] = x;
  q->v.c[1] = y;
  q->v.c[2] = z;

  g_assert_cmpfloat (ncm_trivec_norm (&q->v), >, 0.0);
  ncm_trivec_normalize (&q->v);
  ncm_trivec_scale (&q->v, sin (theta));

  q->s = cos (theta);
}

/**
 * ncm_quaternion_set_I:
 * @q: a #NcmQuaternion
 *
 * Sets @q to the identity, $1 + \vec 0$.
 */
void
ncm_quaternion_set_I (NcmQuaternion *q)
{
  q->s = 1.0;
  ncm_trivec_set_0 (&q->v);
}

/**
 * ncm_quaternion_set_0:
 * @q: a #NcmQuaternion
 *
 * Sets @q to zero.
 */
void
ncm_quaternion_set_0 (NcmQuaternion *q)
{
  q->s = 0.0;
  ncm_trivec_set_0 (&q->v);
}

/**
 * ncm_quaternion_norm:
 * @q: a #NcmQuaternion
 *
 * Returns: $|q|$.
 */
gdouble
ncm_quaternion_norm (NcmQuaternion *q)
{
  return sqrt (q->s * q->s + q->v.c[0] * q->v.c[0] + q->v.c[1] * q->v.c[1] + q->v.c[2] * q->v.c[2]);
}

/**
 * ncm_quaternion_set_random:
 * @q: a #NcmQuaternion
 * @rng: a #NcmRNG
 *
 * Sets @q to a uniformly distributed rotation: four independent standard Gaussians,
 * normalized, give a point uniform on the unit sphere of quaternions.
 */
void
ncm_quaternion_set_random (NcmQuaternion *q, NcmRNG *rng)
{
  ncm_rng_lock (rng);

  do {
    q->s      = ncm_rng_gaussian_gen (rng, 0.0, 1.0);
    q->v.c[0] = ncm_rng_gaussian_gen (rng, 0.0, 1.0);
    q->v.c[1] = ncm_rng_gaussian_gen (rng, 0.0, 1.0);
    q->v.c[2] = ncm_rng_gaussian_gen (rng, 0.0, 1.0);
  } while (ncm_quaternion_norm (q) == 0.0);

  ncm_rng_unlock (rng);

  ncm_quaternion_normalize (q);
}

/**
 * ncm_quaternion_normalize:
 * @q: a #NcmQuaternion
 *
 * Divides @q by its norm.
 */
void
ncm_quaternion_normalize (NcmQuaternion *q)
{
  const gdouble norm = ncm_quaternion_norm (q);

  q->s      /= norm;
  q->v.c[0] /= norm;
  q->v.c[1] /= norm;
  q->v.c[2] /= norm;
}

/**
 * ncm_quaternion_conjugate:
 * @q: a #NcmQuaternion
 *
 * Sets @q to $q^\dagger$.
 */
void
ncm_quaternion_conjugate (NcmQuaternion *q)
{
  ncm_trivec_scale (&q->v, -1.0);
}

/**
 * ncm_quaternion_mul:
 * @q: a #NcmQuaternion
 * @u: a #NcmQuaternion
 * @res: a #NcmQuaternion, not @q or @u
 *
 * Sets @res to $q\,u$.
 */
void
ncm_quaternion_mul (NcmQuaternion *q, NcmQuaternion *u, NcmQuaternion *res)
{
  res->s      = q->s * u->s - ncm_trivec_dot (&q->v, &u->v);
  res->v.c[0] = q->s * u->v.c[0] + u->s * q->v.c[0] + q->v.c[1] * u->v.c[2] - q->v.c[2] * u->v.c[1];
  res->v.c[1] = q->s * u->v.c[1] + u->s * q->v.c[1] + q->v.c[2] * u->v.c[0] - q->v.c[0] * u->v.c[2];
  res->v.c[2] = q->s * u->v.c[2] + u->s * q->v.c[2] + q->v.c[0] * u->v.c[1] - q->v.c[1] * u->v.c[0];
}

/**
 * ncm_quaternion_lmul:
 * @q: a #NcmQuaternion
 * @u: a #NcmQuaternion
 *
 * Sets @q to $u\,q$.
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
 * @q: a #NcmQuaternion
 * @u: a #NcmQuaternion
 *
 * Sets @q to $q\,u$.
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
 * @q: a #NcmQuaternion
 * @u: a #NcmQuaternion
 * @res: a #NcmQuaternion, not @q or @u
 *
 * Sets @res to $q^\dagger u$.
 */
void
ncm_quaternion_conjugate_q_mul (NcmQuaternion *q, NcmQuaternion *u, NcmQuaternion *res)
{
  res->s      = q->s * u->s + ncm_trivec_dot (&q->v, &u->v);
  res->v.c[0] = q->s * u->v.c[0] - u->s * q->v.c[0] - q->v.c[1] * u->v.c[2] + q->v.c[2] * u->v.c[1];
  res->v.c[1] = q->s * u->v.c[1] - u->s * q->v.c[1] - q->v.c[2] * u->v.c[0] + q->v.c[0] * u->v.c[2];
  res->v.c[2] = q->s * u->v.c[2] - u->s * q->v.c[2] - q->v.c[0] * u->v.c[1] + q->v.c[1] * u->v.c[0];
}

/**
 * ncm_quaternion_conjugate_u_mul:
 * @q: a #NcmQuaternion
 * @u: a #NcmQuaternion
 * @res: a #NcmQuaternion, not @q or @u
 *
 * Sets @res to $q\,u^\dagger$.
 */
void
ncm_quaternion_conjugate_u_mul (NcmQuaternion *q, NcmQuaternion *u, NcmQuaternion *res)
{
  res->s      = q->s * u->s + ncm_trivec_dot (&q->v, &u->v);
  res->v.c[0] = -q->s * u->v.c[0] + u->s * q->v.c[0] - q->v.c[1] * u->v.c[2] + q->v.c[2] * u->v.c[1];
  res->v.c[1] = -q->s * u->v.c[1] + u->s * q->v.c[1] - q->v.c[2] * u->v.c[0] + q->v.c[0] * u->v.c[2];
  res->v.c[2] = -q->s * u->v.c[2] + u->s * q->v.c[2] - q->v.c[0] * u->v.c[1] + q->v.c[1] * u->v.c[0];
}

/**
 * ncm_quaternion_rotate:
 * @q: a #NcmQuaternion
 * @v: a #NcmTriVec
 *
 * Sets @v to $q\,\vec v\,q^\dagger$, its rotation by the unit quaternion @q. For a quaternion
 * that is not unit the result is also scaled by $|q|^2$.
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
 * @q: a #NcmQuaternion
 * @v: a #NcmTriVec
 *
 * Sets @v to $q^\dagger\,\vec v\,q$, the inverse of ncm_quaternion_rotate() for a unit @q.
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

/**
 * ncm_quaternion_set_to_rotate_to_x:
 * @q: a #NcmQuaternion
 * @v: a #NcmTriVec
 *
 * Sets @q to the unit quaternion that takes the direction of @v to the x-axis: a rotation
 * about the z-axis into the xz-plane followed by one about the y-axis. For a zero @v, @q is the
 * identity.
 */
void
ncm_quaternion_set_to_rotate_to_x (NcmQuaternion *q, NcmTriVec *v)
{
  const gdouble rho = hypot (v->c[0], v->c[1]);

  /* Rotation about y by atan2 (v_z, rho), which takes (rho, 0, v_z) to the x-axis */
  _ncm_quaternion_set_zy (q, v, atan2 (v->c[2], rho));
}

/**
 * ncm_quaternion_set_to_rotate_to_z:
 * @q: a #NcmQuaternion
 * @v: a #NcmTriVec
 *
 * Sets @q to the unit quaternion that takes the direction of @v to the z-axis: a rotation
 * about the z-axis into the xz-plane followed by one about the y-axis, so no rotation is added
 * about the final axis. For a zero @v, @q is the identity.
 */
void
ncm_quaternion_set_to_rotate_to_z (NcmQuaternion *q, NcmTriVec *v)
{
  const gdouble rho = hypot (v->c[0], v->c[1]);

  /* Rotation about y by -atan2 (rho, v_z), which takes (rho, 0, v_z) to the z-axis */
  _ncm_quaternion_set_zy (q, v, -atan2 (rho, v->c[2]));
}

/*
 * Sets q to R_y (beta) R_z (-phi), with phi the azimuth of v: the rotation about z takes v
 * into the xz-plane with x >= 0, the one about y by beta then to the target axis. On the
 * z-axis phi is 0, and the zero vector gives the identity.
 */
static void
_ncm_quaternion_set_zy (NcmQuaternion *q, NcmTriVec *v, const gdouble beta)
{
  NcmQuaternion t_z = NCM_QUATERNION_INIT_I;
  NcmQuaternion t_y = NCM_QUATERNION_INIT_I;

  if (ncm_trivec_norm (v) == 0.0)
  {
    ncm_quaternion_set_I (q);

    return;
  }

  {
    const gdouble phi = atan2 (v->c[1], v->c[0]);

    t_z.s      = cos (0.5 * phi);
    t_z.v.c[2] = -sin (0.5 * phi);
    t_y.s      = cos (0.5 * beta);
    t_y.v.c[1] = sin (0.5 * beta);
  }

  ncm_quaternion_mul (&t_y, &t_z, q);
}

