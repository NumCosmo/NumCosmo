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
 * SECTION:ncm_quaternion
 * @title: NcmQuaternion
 * @short_description: Quaternions algebra, three-vectors and mapping to matrix.
 *
 * A quaternion is a four-dimensional vector that can be used to represent rotations in
 * three-dimensional space. The three-dimensional space is represented by the
 * three-dimensional subspace of the quaternions that have zero real part.
 *
 * This object also implements three-dimensional vectors and the mapping of quaternions
 * to rotation matrices.
 *
 * The conjugate of a quaternion is the quaternion with the same real part and the
 * imaginary part negated, that is, if $q = s + \vec{v}$, then the conjugate of $q$ is
 * $q^\dagger = s - \vec{v}$.
 *
 * The norm of a quaternion is the square root of the sum of the squares of its components.
 * The norm of a quaternion is always a positive real number.
 *
 *
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_quaternion.h"
#include "math/ncm_c.h"
#include "math/ncm_rng.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <string.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_DEFINE_BOXED_TYPE (NcmQuaternion, ncm_quaternion, ncm_quaternion_dup, ncm_quaternion_free)
G_DEFINE_BOXED_TYPE (NcmTriVec, ncm_trivec, ncm_trivec_dup, ncm_trivec_free)

/**
 * ncm_trivec_new: (constructor)
 *
 * Creates a new empty #NcmTriVec.
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
 * @c: components
 *
 * Creates a new #NcmTriVec with the given components.
 *
 * Returns: (transfer full): the new #NcmTriVec.
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
 * Creates a new #NcmTriVec with the given components.
 *
 * Returns: (transfer full): the new #NcmTriVec.
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
 * @phi: the azimuthal angle
 *
 * Creates a new #NcmTriVec with the given spherical coordinates. The spherical
 * coordinates are the radius, the polar angle and the azimuthal angle. The polar angle
 * is the angle between the vector and the z-axis, and the azimuthal angle is the angle
 * between the projection of the vector in the xy-plane and the x-axis.
 *
 * Returns: (transfer full): the new #NcmTriVec.
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
 * @delta: the declination
 * @alpha: the right ascension
 *
 * Creates a new #NcmTriVec with the given astronomical coordinates. The astronomical
 * coordinates are the declination and the right ascension. The declination is the angle
 * between the vector and the celestial equator, and the right ascension is the angle
 * between the projection of the vector in the celestial equator and the vernal equinox.
 * The angles are in radians.
 *
 * Returns: (transfer full): the new #NcmTriVec.
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
 * @ra: the right ascension
 * @dec: the declination
 *
 * Creates a new #NcmTriVec with the given astronomical coordinates. The astronomical
 * coordinates are the declination and the right ascension. The declination is the angle
 * between the vector and the celestial equator, and the right ascension is the angle
 * between the projection of the vector in the celestial equator and the vernal equinox.
 * The angles are in degrees.
 *
 * Returns: (transfer full): the new #NcmTriVec.
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
 * Duplicates a #NcmTriVec.
 *
 * Returns: (transfer full): a new #NcmTriVec.
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
 * Frees a #NcmTriVec.
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
 * Copies a #NcmTriVec.
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
 * Sets a #NcmTriVec to zero.
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
 * Scale a #NcmTriVec.
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
 * Calculates the norm of a #NcmTriVec.
 *
 * Returns: the norm of @v.
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
 * Calculates the dot product of two #NcmTriVec.
 *
 * Returns: the dot product of @v1 and @v2.
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
 * Normalize a #NcmTriVec.
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
 * Gets the azimuthal angle of a #NcmTriVec.
 *
 * Returns: the azimuthal angle of @v.
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
 * @phi: the azimuthal angle
 *
 * Sets the spherical coordinates of a #NcmTriVec. The spherical coordinates are the radius,
 * the polar angle and the azimuthal angle. The polar angle is the angle between the vector
 * and the z-axis, and the azimuthal angle is the angle between the projection of the vector
 * in the xy-plane and the x-axis.
 *
 * The vector is defined as:
 * $$\vec{v} = (r \sin(\theta) \cos(\phi), r \sin(\theta) \sin(\phi), r \cos(\theta)).$$
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
 * @r: (out): the radius
 * @theta: (out): the polar angle
 * @phi: (out): the azimuthal angle
 *
 * Computes the spherical coordinates of a #NcmTriVec.
 *
 */
void
ncm_trivec_get_spherical_coord (NcmTriVec *v, gdouble *r, gdouble *theta, gdouble *phi)
{
  const gdouble norm = ncm_trivec_norm (v);

  r[0]     = norm;
  theta[0] = acos (v->c[2] / norm);
  phi[0]   = ncm_trivec_get_phi (v);
}

/**
 * ncm_trivec_set_astro_coord:
 * @v: a #NcmTriVec
 * @r: the radius
 * @delta: the declination
 * @alpha: the right ascension
 *
 * Sets the astronomical coordinates of a #NcmTriVec. The astronomical coordinates are
 * the declination and the right ascension. The vector is defined as:
 * $$\vec{v} = (\cos(\delta) \cos(\alpha), \cos(\delta) \sin(\alpha), \sin(\delta)).$$
 *
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
 * @delta: (out): the declination
 * @alpha: (out): the right ascension
 *
 * Computes the astronomical coordinates of a #NcmTriVec.
 * See ncm_trivec_set_astro_coord() for details.
 *
 */
void
ncm_trivec_get_astro_coord (NcmTriVec *v, gdouble *r, gdouble *delta, gdouble *alpha)
{
  const gdouble norm = ncm_trivec_norm (v);

  r[0]     = norm;
  delta[0] = asin (v->c[2] / norm);
  alpha[0] = atan2 (v->c[1], v->c[0]);
}

/**
 * ncm_trivec_set_astro_ra_dec:
 * @v: a #NcmTriVec
 * @r: the radius
 * @ra: the right ascension (in degrees)
 * @dec: the declination (in degrees)
 *
 * Sets the astronomical coordinates of a #NcmTriVec. The astronomical coordinates are
 * the declination and the right ascension. The declination is the angle between the vector
 * and the z-axis, and the right ascension is the angle between the projection of the vector
 * in the xy-plane and the x-axis. The declination and the right ascension are given in degrees.
 *
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
 * @ra: (out): the right ascension (in degrees)
 * @dec: (out): the declination (in degrees)
 *
 * Computes the astronomical coordinates of a #NcmTriVec.
 * See ncm_trivec_set_astro_ra_dec() for details.
 *
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
 * Creates a new empty #NcmQuaternion.
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
 * Creates a new #NcmQuaternion from a #NcmTriVec.
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
 * @x: the x-component
 * @y: the y-component
 * @z: the z-component
 * @theta: the angle
 *
 * Creates a new #NcmQuaternion from the given components.
 * See ncm_quaternion_set_from_data() for details.
 *
 * Returns: (transfer full): a new #NcmQuaternion.
 */
NcmQuaternion *
ncm_quaternion_new_from_data (gdouble x, gdouble y, gdouble z, gdouble theta)
{
  NcmQuaternion *q = ncm_quaternion_new ();

  theta    /= 2.0;
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
 * Duplicates a #NcmQuaternion.
 *
 * Returns: (transfer full): a new #NcmQuaternion.
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
 * Frees a #NcmQuaternion.
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
 * Copies a #NcmQuaternion.
 *
 */
void
ncm_quaternion_memcpy (NcmQuaternion *dest, const NcmQuaternion *orig)
{
  memcpy (dest, orig, sizeof (NcmQuaternion));
}

/**
 * ncm_quaternion_set_from_data:
 * @q: a #NcmQuaternion
 * @x: the x-component
 * @y: the y-component
 * @z: the z-component
 * @theta: the angle
 *
 * Sets the components of a #NcmQuaternion.
 * The components are the components of a three-dimensional vector and the angle
 * of rotation, the three-dimensional vector is normalized. The final
 * form of the quaternion is:
 * $$q = \cos(\theta/2) + \sin(\theta/2) \hat{v}.$$
 *
 */
void
ncm_quaternion_set_from_data (NcmQuaternion *q, gdouble x, gdouble y, gdouble z, gdouble theta)
{
  theta    /= 2.0;
  q->v.c[0] = x;
  q->v.c[1] = y;
  q->v.c[2] = z;

  ncm_trivec_normalize (&q->v);
  ncm_trivec_scale (&q->v, sin (theta));

  q->s = cos (theta);
}

/**
 * ncm_quaternion_set_I:
 * @q: a #NcmQuaternion
 *
 * Sets a #NcmQuaternion to the identity.
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
 * @q: a #NcmQuaternion
 *
 * Sets a #NcmQuaternion to zero.
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
 * @q: a #NcmQuaternion
 *
 * Calculates the norm of a #NcmQuaternion.
 *
 * Returns: the norm of @q.
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
 * Sets a #NcmQuaternion to a random value, using the given #NcmRNG.
 * The components of the three-dimensional vector are uniformly distributed
 * in the interval [-1, 1] and the angle is uniformly distributed in the
 * interval [0, 2*pi].
 *
 * It represents a random rotation in three-dimensional space.
 *
 */
void
ncm_quaternion_set_random (NcmQuaternion *q, NcmRNG *rng)
{
  ncm_rng_lock (rng);
  ncm_quaternion_set_from_data (q,
                                -1.0 + 2.0 * ncm_rng_uniform01_pos_gen (rng),
                                -1.0 + 2.0 * ncm_rng_uniform01_pos_gen (rng),
                                -1.0 + 2.0 * ncm_rng_uniform01_pos_gen (rng),
                                2.0 * M_PI * ncm_rng_uniform01_gen (rng));
  ncm_rng_unlock (rng);
}

/**
 * ncm_quaternion_normalize:
 * @q: a #NcmQuaternion
 *
 * Normalize a #NcmQuaternion.
 *
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
 * Conjugate a #NcmQuaternion. That is, the vector part is negated.
 *
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
 * @res: a #NcmQuaternion
 *
 * Computes the product of two #NcmQuaternion $r=qu$ where $q$ and $u$ are
 * @q and @u, respectively, and $r$ is @res.
 *
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
 * Computes the product of two #NcmQuaternion and stores the result in @q.
 * That is, @q = @u * @q, where @u and @q are @u and @q, respectively.
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
 * @q: a #NcmQuaternion
 * @u: a #NcmQuaternion
 *
 * Computes the product of two #NcmQuaternion and stores the result in @q.
 * That is, @q = @q * @u, where @q and @u are @q and @u, respectively.
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
 * @q: a #NcmQuaternion
 * @u: a #NcmQuaternion
 * @res: a #NcmQuaternion
 *
 * Computes the product of two #NcmQuaternion and stores the result in @res.
 * The first #NcmQuaternion is conjugated before the multiplication.
 * That is, $r = q^\dagger u$ where $q$ and $u$ are @q and @u, respectively, and $r$ is @res.
 *
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
 * @res: a #NcmQuaternion
 *
 * Computes the product of two #NcmQuaternion and stores the result in @res.
 * The second #NcmQuaternion is conjugated before the multiplication.
 * The result is $r = q u^\dagger$, where $q$ and $u$ are @q and @u, respectively,
 * and $r$ is @res. The conjugation is done by negating the vector part of @u.
 *
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
 * Computes the rotation of a #NcmTriVec by a #NcmQuaternion.
 * The rotation is done by the formula $v' = q v q^\dagger$.
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
 * @q: a #NcmQuaternion
 * @v: a #NcmTriVec
 *
 * Computes the inverse rotation of a #NcmTriVec by a #NcmQuaternion.
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

/**
 * ncm_quaternion_set_to_rotate_to_x:
 * @q: a #NcmQuaternion
 * @v: a #NcmTriVec
 *
 * Sets @q to the rotation that rotates the given #NcmTriVec to the x-axis. It finds
 * first the quaternion that rotates the given vector to the xz-plane and then the
 * quaternion that rotates the vector to the x-axis. Finally, it multiplies the two
 * quaternions and stores the result in @q.
 *
 */
void
ncm_quaternion_set_to_rotate_to_x (NcmQuaternion *q, NcmTriVec *v)
{
  NcmQuaternion t1       = NCM_QUATERNION_INIT_I;
  NcmQuaternion t2       = NCM_QUATERNION_INIT_I;
  const gdouble vx       = v->c[0];
  const gdouble vy       = v->c[1];
  const gdouble vz       = v->c[2];
  const gdouble norma_xy = hypot (vx, vy);
  const gdouble norma    = hypot (vz, norma_xy);

  if (norma > 0.0)
  {
    /*
     * Find the quaternion that rotates the vector to the xz-plane. Skip if the vector
     * is already in the xz-plane or if has no component in the xy-plane.
     */
    const gdouble nx = vx / norma_xy;
    const gdouble ny = vy / norma_xy;
    const gdouble ux = norma_xy / norma;
    const gdouble uz = vz / norma;

    if ((norma_xy > 0.0) && (ny != 0.0))
    {
      /*
       * To avoid division by zero, we use the following formula to find the quaternion
       * that rotates the vector to the xz-plane:
       */
      if ((fabs (ny) < 0.1) && (nx > 0.0))
        t1.v.c[2] = -ny / (1.0 + sqrt (1.0 - ny * ny));
      else
        t1.v.c[2] = -(1.0 - nx) / ny;
    }

    /*
     * Find the quaternion that rotates the vector to the x-axis. Skip if the vector is
     * already in the x-axis.
     */
    if (uz != 0.0)
    {
      /*
       * To avoid division by zero, we use the following formula to find the quaternion
       * that rotates the vector to the x-axis:
       */
      if ((fabs (uz) < 0.1) && (ux > 0.0))
        t2.v.c[1] = uz / (1.0 + sqrt (1.0 - uz * uz));
      else
        t2.v.c[1] = (1.0 - ux) / uz;
    }
  }

  ncm_quaternion_mul (&t2, &t1, q);
  ncm_quaternion_normalize (q);
}

/**
 * ncm_quaternion_set_to_rotate_to_z:
 * @q: a #NcmQuaternion
 * @v: a #NcmTriVec
 *
 * Sets @q to the rotation that rotates the given #NcmTriVec to the z-axis. It finds
 * first the quaternion that rotates the given vector to the xz-plane and then the
 * quaternion that rotates the vector to the z-axis. Finally, it multiplies the two
 * quaternions and stores the result in @q.
 *
 */
void
ncm_quaternion_set_to_rotate_to_z (NcmQuaternion *q, NcmTriVec *v)
{
  NcmQuaternion t1       = NCM_QUATERNION_INIT_I;
  NcmQuaternion t2       = NCM_QUATERNION_INIT_I;
  const gdouble vx       = v->c[0];
  const gdouble vy       = v->c[1];
  const gdouble vz       = v->c[2];
  const gdouble norma_xy = hypot (vx, vy);
  const gdouble norma    = hypot (vz, norma_xy);

  if (norma > 0.0)
  {
    /*
     * Find the quaternion that rotates the vector to the xz-plane. Skip if the vector
     * is already in the xz-plane or if has no component in the xy-plane.
     */
    const gdouble nx = vx / norma_xy;
    const gdouble ny = vy / norma_xy;
    const gdouble ux = norma_xy / norma;
    const gdouble uz = vz / norma;

    if ((norma_xy > 0.0) && (ny != 0.0))
    {
      /*
       * To avoid division by zero, we use the following formula to find the quaternion
       * that rotates the vector to the xz-plane:
       */
      if ((fabs (ny) < 0.1) && (nx > 0.0))
        t1.v.c[2] = -ny / (1.0 + sqrt (1.0 - ny * ny));
      else
        t1.v.c[2] = -(1.0 - nx) / ny;
    }

    /*
     * Find the quaternion that rotates the vector to the z-axis. Skip if the vector is
     * already in the z-axis.
     */
    if (ux != 0.0)
    {
      /*
       * To avoid division by zero, we use the following formula to find the quaternion
       * that rotates the vector to the z-axis:
       */
      if ((fabs (ux) < 0.1) && (uz > 0.0))
        t2.v.c[1] = -ux / (1.0 + sqrt (1.0 - ux * ux));
      else
        t2.v.c[1] = -(1.0 - uz) / ux;
    }
  }

  ncm_quaternion_mul (&t2, &t1, q);
  ncm_quaternion_normalize (q);
}

