/***************************************************************************
 *            quaternion.c
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
 * SECTION:quaternion
 * @title: NcmQ
 * @short_description: Quaternions algebra and mapping to matrix.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/quaternion.h"
#include "math/ncm_rng.h"

#include <string.h>
#include <math.h>

/**
 * ncm_quaternion_new: (skip)
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmQ *
ncm_quaternion_new ()
{
  NcmQ *q = (NcmQ *) g_slice_new (NcmQ);
  NC_QUATERNION_SET_0(q);
  return q;
}

/**
 * ncm_quaternion_free:
 * @q: a #NcmQ
 *
 * FIXME
*/
void
ncm_quaternion_free (NcmQ *q)
{
  g_slice_free (NcmQ, q);
}

/**
 * ncm_quaternion_new_from_vector: (skip)
 * @v: a #NcmTriVector
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmQ *
ncm_quaternion_new_from_vector (NcmTriVector v)
{
  NcmQ *q = (NcmQ *) g_slice_new (NcmQ);
  NC_TRIVEC_MEMCPY(q->x, v);
  q->s = 0.0;
  return q;
}

/**
 * ncm_quaternion_new_from_data: (skip)
 * @x: FIXME
 * @y: FIXME
 * @z: FIXME
 * @theta: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmQ *
ncm_quaternion_new_from_data (gdouble x, gdouble y, gdouble z, gdouble theta)
{
  NcmQ *q = (NcmQ *) g_slice_new (NcmQ);
  theta /= 2.0;
  q->x.c[0] = x;
  q->x.c[1] = y;
  q->x.c[2] = z;
  NC_TRIVEC_NORMALIZE(q->x);
  NC_TRIVEC_SCALE(q->x, sin(theta));
  q->s = cos(theta);
  return q;
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
ncm_quaternion_set_from_data (NcmQ *q, gdouble x, gdouble y, gdouble z, gdouble theta)
{
  theta /= 2.0;
  q->x.c[0] = x;
  q->x.c[1] = y;
  q->x.c[2] = z;
  NC_TRIVEC_NORMALIZE(q->x);
  NC_TRIVEC_SCALE(q->x, sin(theta));
  q->s = cos(theta);
}

/**
 * ncm_quaternion_normalize:
 * @q: FIXME
 *
 * FIXME
 *
 */
void
ncm_quaternion_normalize (NcmQ *q)
{
  gdouble norm = NC_QUATERNION_NORM(q);
  if (fabs(norm - 1.0) > 1e7)
  {
    q->s /= norm;
    q->x.c[0] /= norm;
    q->x.c[1] /= norm;
    q->x.c[2] /= norm;
  }
}

/**
 * ncm_quaternion_conjugate:
 * @q: FIXME
 *
 * FIXME
 *
 */
void
ncm_quaternion_conjugate (NcmQ *q)
{
  NC_TRIVEC_SCALE(q->x, -1.0);
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
ncm_quaternion_rotate (NcmQ *q, NcmTriVector v)
{
  NcmQ qv = NC_QUATERNION_NEW;
  NcmQ t;
  NC_TRIVEC_MEMCPY (qv.x, v);
  ncm_quaternion_mul (q, &qv, &t);
  ncm_quaternion_conjugate_u_mul (&t, q, &qv);
  NC_TRIVEC_MEMCPY(v, qv.x);
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
ncm_quaternion_inv_rotate (NcmQ *q, NcmTriVector v)
{
  NcmQ qv = NC_QUATERNION_NEW;
  NcmQ t;
  NC_TRIVEC_MEMCPY (qv.x, v);

  ncm_quaternion_mul (&qv, q, &t);
  ncm_quaternion_conjugate_q_mul (q, &t, &qv);

  NC_TRIVEC_MEMCPY(v, qv.x);
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
ncm_quaternion_mul (NcmQ *q, NcmQ *u, NcmQ *res)
{
  res->s = q->s * u->s - NC_TRIVEC_DOT(q->x, u->x);
  res->x.c[0] = q->s * u->x.c[0] + u->s * q->x.c[0] + q->x.c[1] * u->x.c[2] - q->x.c[2] * u->x.c[1];
  res->x.c[1] = q->s * u->x.c[1] + u->s * q->x.c[1] + q->x.c[2] * u->x.c[0] - q->x.c[0] * u->x.c[2];
  res->x.c[2] = q->s * u->x.c[2] + u->s * q->x.c[2] + q->x.c[0] * u->x.c[1] - q->x.c[1] * u->x.c[0];
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
ncm_quaternion_lmul (NcmQ *q, NcmQ *u)
{
  NcmQ t;
  ncm_quaternion_mul (u, q, &t);
  NC_QUATERNION_MEMCPY(q, &t);
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
ncm_quaternion_rmul (NcmQ *q, NcmQ *u)
{
  NcmQ t;
  ncm_quaternion_mul (q, u, &t);
  NC_QUATERNION_MEMCPY(q, &t);
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
ncm_quaternion_conjugate_q_mul (NcmQ *q, NcmQ *u, NcmQ *res)
{
  res->s = q->s * u->s + NC_TRIVEC_DOT(q->x, u->x);
  res->x.c[0] = q->s * u->x.c[0] - u->s * q->x.c[0] - q->x.c[1] * u->x.c[2] + q->x.c[2] * u->x.c[1];
  res->x.c[1] = q->s * u->x.c[1] - u->s * q->x.c[1] - q->x.c[2] * u->x.c[0] + q->x.c[0] * u->x.c[2];
  res->x.c[2] = q->s * u->x.c[2] - u->s * q->x.c[2] - q->x.c[0] * u->x.c[1] + q->x.c[1] * u->x.c[0];
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
ncm_quaternion_conjugate_u_mul (NcmQ *q, NcmQ *u, NcmQ *res)
{
  res->s = q->s * u->s + NC_TRIVEC_DOT(q->x, u->x);
  res->x.c[0] = - q->s * u->x.c[0] + u->s * q->x.c[0] - q->x.c[1] * u->x.c[2] + q->x.c[2] * u->x.c[1];
  res->x.c[1] = - q->s * u->x.c[1] + u->s * q->x.c[1] - q->x.c[2] * u->x.c[0] + q->x.c[0] * u->x.c[2];
  res->x.c[2] = - q->s * u->x.c[2] + u->s * q->x.c[2] - q->x.c[0] * u->x.c[1] + q->x.c[1] * u->x.c[0];
}

/**
 * ncm_quaternion_set_random:
 * @q: FIXME
 *
 * FIXME
 *
 */
void
ncm_quaternion_set_random (NcmQ *q)
{
  NcmRNG *rng = ncm_rng_pool_get (NCM_QUATERNION_RNG_NAME);
  ncm_rng_lock (rng);
  ncm_quaternion_set_from_data (q,
                                 -1.0 + 2.0 * gsl_rng_uniform_pos (rng->r),
                                 -1.0 + 2.0 * gsl_rng_uniform_pos (rng->r),
                                 -1.0 + 2.0 * gsl_rng_uniform_pos (rng->r),
                                 2.0 * M_PI * gsl_rng_uniform (rng->r));
  ncm_rng_unlock (rng);
  ncm_rng_free (rng);   
}
