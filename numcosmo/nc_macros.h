/***************************************************************************
 *            nc_macros.h
 *
 *  Tue Mar 18 13:34:04 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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

#ifndef _NC_MACROS_H
#define _NC_MACROS_H

#include <glib.h>

G_BEGIN_DECLS

#ifdef SUNDIALS_USES_LONG_INT
#define _NCM_SUNDIALS_INT_TYPE glong
#else
#define _NCM_SUNDIALS_INT_TYPE gint
#endif

#if (GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 32)
#define _NCM_MUTEX_LOCK(l) g_static_mutex_lock (l)
#define _NCM_MUTEX_UNLOCK(l) g_static_mutex_unlock (l)
#define _NCM_MUTEX_TRYLOCK(l) g_static_mutex_trylock (l)
#define _NCM_MUTEX_TYPE GStaticMutex
#define _NCM_STATIC_MUTEX_DECL(l) static GStaticMutex l = G_STATIC_MUTEX_INIT
#define _NCM_MUTEX_INIT(l) g_static_mutex_init (l)
#define _NCM_MUTEX_CLEAR(l) g_static_mutex_free (l)
#else
#define _NCM_MUTEX_LOCK(l) g_mutex_lock (l)
#define _NCM_MUTEX_UNLOCK(l) g_mutex_unlock (l)
#define _NCM_MUTEX_TRYLOCK(l) g_mutex_trylock (l)
#define _NCM_MUTEX_TYPE GMutex
#define _NCM_STATIC_MUTEX_DECL(l) static GMutex l
#define _NCM_MUTEX_INIT(l) g_mutex_init (l)
#define _NCM_MUTEX_CLEAR(l) g_mutex_clear (l)
#endif

#define NC_DEGREE_TO_RADIAN(a) ((a) * M_PI/180.0)
#define NC_RADIAN_TO_DEGREE(a) ((a) * 180.0/M_PI)
#define NC_RADIAN_0_2PI(a) ((a)-2.0 * M_PI * floor((a) / (2.0 * M_PI)))
#define NC_SIGN_SIN(a) ((NC_RADIAN_0_2PI(a) < M_PI) ? 1.0 : -1.0)

#define NC_RETURN_IF_INF(a) if (gsl_isinf(a)) return a

#define NC_FLOOR_TRUNC(a,b) (floor ((b) * (a)) / (b))
#define NC_CEIL_TRUNC(a,b) (ceil ((b) * (a)) / (b))
#define NC_ROUND_TRUNC(a,b) (round ((b) * (a)) / (b))

#define NC_TEST_GSL_RESULT(func,ret) if (ret != GSL_SUCCESS) g_error ("%s: %s", func, gsl_strerror (ret))

#define NC_BF_MAX_ITER 100000
#define NC_ZERO_LIMIT 1e-13
#define NC_DEFAULT_PRECISION 1e-7

#ifndef NCM_THREAD_POOL_MAX
#define NCM_THREAD_POOL_MAX 5
#endif

#define NC_MAP_ALM_SIZE(lmax) (((lmax)*(lmax) + 3*(lmax) + 2)/2)
#define NC_MAP_N_IND_PLM(lmax) (((lmax)*(lmax) + 3*(lmax) + 2)/2)
#define NC_MAP_MAX_RING_SIZE(nside) (4*(nside))
#define NC_MAP_N_DIFFERENT_SIZED_RINGS(nside) (nside)
#define NC_MAP_RING_PLAN_INDEX(nside,ring_n) (((ring_n) < (nside)) ? (ring_n) : ((ring_n)>=(3*(nside)) ? (4*(nside)-(ring_n)-2) : ((nside)-1)))
#define NC_MAP_RING_SIZE(nside,ring_n) (4*(NC_MAP_RING_PLAN_INDEX(nside,ring_n)+1))
#define NC_MAP_N_RINGS(nside) (4*(nside)-1)
#define NC_MAP_ALM_M_START(lmax,m) ((2*(lmax)*(m)-(m)*(m)+3*(m))/2)
#define NC_MAP_ALM_INDEX(lmax,l,m) (((l) >= (m)) ? (NC_MAP_ALM_M_START(lmax,m) + (l) - (m)) : (-1))

#ifndef mpz_inits
#define mpz_inits ncm_mpz_inits
#endif /* mpz_inits */

#ifndef mpz_clears
#define mpz_clears ncm_mpz_clears
#endif /* mpz_inits */

#define NC_COMPLEX_INC_MUL_REAL_TEST(a,b,c) \
((fabs(GSL_REAL((b))*(c)/GSL_REAL((a))) < 1e-16) && (fabs(GSL_IMAG((b))*(c)/GSL_IMAG((a))) < 1e-16))

#define NC_COMPLEX_INC_MUL_REAL(a,b,c) \
do { \
  GSL_REAL((a)) += GSL_REAL((b))*(c); \
  GSL_IMAG((a)) += GSL_IMAG((b))*(c); \
} while (FALSE)

#define NC_COMPLEX_INC_MUL(a,b,c) \
do { \
  GSL_REAL((a)) += GSL_REAL((b)) * GSL_REAL((c)) - GSL_IMAG((b)) * GSL_IMAG((c)); \
  GSL_IMAG((a)) += GSL_REAL((b)) * GSL_IMAG((c)) + GSL_IMAG((b)) * GSL_REAL((c)); \
} while (FALSE)

#define NC_COMPLEX_INC_MUL_MUL_REAL(a,b,c,d) \
do { \
  GSL_REAL((a)) += (GSL_REAL((b)) * GSL_REAL((c)) - GSL_IMAG((b)) * GSL_IMAG((c))) * (d); \
  GSL_IMAG((a)) += (GSL_REAL((b)) * GSL_IMAG((c)) + GSL_IMAG((b)) * GSL_REAL((c))) * (d); \
} while (FALSE)

#define NC_COMPLEX_MUL_REAL(a,b,c) \
do { \
  GSL_REAL((a)) = GSL_REAL((b)) * (c); \
  GSL_IMAG((a)) = GSL_IMAG((b)) * (c); \
} while (FALSE)

#define NC_COMPLEX_MUL(a,b) \
do { \
  gdouble temp = GSL_REAL((a)) * GSL_REAL((b)) - GSL_IMAG((a)) * GSL_IMAG((b)); \
  GSL_IMAG(a) = GSL_REAL((a)) * GSL_IMAG((b)) + GSL_IMAG((a)) * GSL_REAL((b)); \
  GSL_REAL(a) = temp; \
} while (FALSE)

#define NC_COMPLEX_ADD(a,b) \
do { \
  GSL_REAL(a) += GSL_REAL((b)); \
  GSL_IMAG(a) += GSL_IMAG((b)); \
} while (FALSE)

#define NC_COMPLEX_MUL_CONJUGATE(a,b) \
do { \
  gdouble temp = GSL_REAL((a)) * GSL_REAL((b)) + GSL_IMAG((a)) * GSL_IMAG((b)); \
  GSL_IMAG(a) = - GSL_REAL((a)) * GSL_IMAG((b)) + GSL_IMAG((a)) * GSL_REAL((b)); \
  GSL_REAL(a) = temp; \
} while (FALSE)

G_END_DECLS

#endif /* _NC_MACROS_H */
