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
 * @title: Quanternions
 * @short_description: Quaternions algebra and mapping to matrix
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <string.h>
#include <glib.h>
#include <glib/gstdio.h>
#include <gsl/gsl_rng.h>


/**
 * nc_quaternion_new: (skip)
 *
 * FIXME
 *
 * Returns: FIXME 
 */
NcQ *
nc_quaternion_new ()
{
  NcQ *q = (NcQ *)g_slice_new (NcQ);
  NC_QUATERNION_SET_0(q);
  return q;
}

/**
 * nc_quaternion_free:
 * @q: a #NcQ
 *
 * FIXME
*/
void
nc_quaternion_free (NcQ *q)
{
  g_slice_free (NcQ, q);
}

/**
 * nc_quaternion_new_from_vector: (skip)
 * @v: a #NcTriVector
 * 
 * FIXME
 *
 * Returns: FIXME
*/
NcQ *
nc_quaternion_new_from_vector (NcTriVector v)
{
  NcQ *q = (NcQ *)g_slice_new (NcQ);
  NC_TRIVEC_MEMCPY(q->x, v);
  q->s = 0.0;
  return q;
}

/**
 * nc_quaternion_new_from_data: (skip)
 * @x: FIXME
 * @y: FIXME
 * @z: FIXME
 * @theta: FIXME
 * 
 * FIXME
 *
 * Returns: FIXME
 */
NcQ *
nc_quaternion_new_from_data (gdouble x, gdouble y, gdouble z, gdouble theta)
{
  NcQ *q = (NcQ *)g_slice_new (NcQ);
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
 * FIXME
 */
void 
nc_quaternion_set_from_data (NcQ *q, gdouble x, gdouble y, gdouble z, gdouble theta)
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
 * FIXME
 */
void
nc_quaternion_normalize (NcQ *q)
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
 * FIXME
 */
void
nc_quaternion_conjugate (NcQ *q)
{
  NC_TRIVEC_SCALE(q->x, -1.0);
}

/**
 * FIXME
 */
void
nc_quaternion_rotate (NcQ *q, NcTriVector v)
{
  NcQ qv = NC_QUATERNION_NEW;
  NcQ t;
  NC_TRIVEC_MEMCPY (qv.x, v);
  nc_quaternion_mul (q, &qv, &t);
  nc_quaternion_conjugate_u_mul (&t, q, &qv);
  NC_TRIVEC_MEMCPY(v, qv.x);
}

/**
 * FIXME
 */
void
nc_quaternion_inv_rotate (NcQ *q, NcTriVector v)
{
  NcQ qv = NC_QUATERNION_NEW;
  NcQ t;
  NC_TRIVEC_MEMCPY (qv.x, v);
  
  nc_quaternion_mul (&qv, q, &t);
  nc_quaternion_conjugate_q_mul (q, &t, &qv);
  
  NC_TRIVEC_MEMCPY(v, qv.x);
}

/**
 * FIXME
 */
void
nc_quaternion_mul (NcQ *q, NcQ *u, NcQ *res)
{
  res->s = q->s * u->s - NC_TRIVEC_DOT(q->x, u->x);
  res->x.c[0] = q->s * u->x.c[0] + u->s * q->x.c[0] + q->x.c[1] * u->x.c[2] - q->x.c[2] * u->x.c[1];
  res->x.c[1] = q->s * u->x.c[1] + u->s * q->x.c[1] + q->x.c[2] * u->x.c[0] - q->x.c[0] * u->x.c[2];
  res->x.c[2] = q->s * u->x.c[2] + u->s * q->x.c[2] + q->x.c[0] * u->x.c[1] - q->x.c[1] * u->x.c[0];
}

/**
 * FIXME
 */
void
nc_quaternion_lmul (NcQ *q, NcQ *u)
{
  NcQ t;
  nc_quaternion_mul (u, q, &t);
  NC_QUATERNION_MEMCPY(q, &t);
}

/**
 * FIXME
 */
void
nc_quaternion_rmul (NcQ *q, NcQ *u)
{
  NcQ t;
  nc_quaternion_mul (q, u, &t);
  NC_QUATERNION_MEMCPY(q, &t);
}

/**
 * FIXME
 */
void
nc_quaternion_conjugate_q_mul (NcQ *q, NcQ *u, NcQ *res)
{
  res->s = q->s * u->s + NC_TRIVEC_DOT(q->x, u->x);
  res->x.c[0] = q->s * u->x.c[0] - u->s * q->x.c[0] - q->x.c[1] * u->x.c[2] + q->x.c[2] * u->x.c[1];
  res->x.c[1] = q->s * u->x.c[1] - u->s * q->x.c[1] - q->x.c[2] * u->x.c[0] + q->x.c[0] * u->x.c[2];
  res->x.c[2] = q->s * u->x.c[2] - u->s * q->x.c[2] - q->x.c[0] * u->x.c[1] + q->x.c[1] * u->x.c[0];
}

/**
 * FIXME
 */
void
nc_quaternion_conjugate_u_mul (NcQ *q, NcQ *u, NcQ *res)
{
  res->s = q->s * u->s + NC_TRIVEC_DOT(q->x, u->x);
  res->x.c[0] = - q->s * u->x.c[0] + u->s * q->x.c[0] - q->x.c[1] * u->x.c[2] + q->x.c[2] * u->x.c[1];
  res->x.c[1] = - q->s * u->x.c[1] + u->s * q->x.c[1] - q->x.c[2] * u->x.c[0] + q->x.c[0] * u->x.c[2];
  res->x.c[2] = - q->s * u->x.c[2] + u->s * q->x.c[2] - q->x.c[0] * u->x.c[1] + q->x.c[1] * u->x.c[0];
}

void 
nc_quaternion_set_random (NcQ *q)
{
  gsl_rng *rand = ncm_get_rng();
  nc_quaternion_set_from_data (q, 
                                 -1.0 + 2.0 * gsl_rng_uniform_pos (rand), 
                                 -1.0 + 2.0 * gsl_rng_uniform_pos (rand), 
                                 -1.0 + 2.0 * gsl_rng_uniform_pos (rand), 
                                 2.0 * M_PI * gsl_rng_uniform (rand));
}
