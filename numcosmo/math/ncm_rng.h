/***************************************************************************
 *            ncm_rng.h
 *
 *  Sat August 17 12:39:57 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_rng.h
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_RNG_H_
#define _NCM_RNG_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_RNG             (ncm_rng_get_type ())
#define NCM_RNG(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_RNG, NcmRNG))
#define NCM_RNG_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_RNG, NcmRNGClass))
#define NCM_IS_RNG(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_RNG))
#define NCM_IS_RNG_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_RNG))
#define NCM_RNG_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_RNG, NcmRNGClass))

typedef struct _NcmRNGClass NcmRNGClass;
typedef struct _NcmRNG NcmRNG;

struct _NcmRNGClass
{
  /*< private >*/
  GObjectClass parent_class;
  GRand *seed_gen;
  GHashTable *seed_hash;
};

struct _NcmRNG
{
  /*< private >*/
  GObject parent_instance;
  gsl_rng *r;
  gulong seed_val;
  gboolean seed_set;
  GMutex lock;
};

GType ncm_rng_get_type (void) G_GNUC_CONST;

NcmRNG *ncm_rng_new (const gchar *algo);
NcmRNG *ncm_rng_seeded_new (const gchar *algo, gulong seed);
NcmRNG *ncm_rng_ref (NcmRNG *rng);
void ncm_rng_free (NcmRNG *rng);
void ncm_rng_clear (NcmRNG **rng);

void ncm_rng_lock (NcmRNG *rng);
void ncm_rng_unlock (NcmRNG *rng);

const gchar *ncm_rng_get_algo (NcmRNG *rng);
gchar *ncm_rng_get_state (NcmRNG *rng);
void ncm_rng_set_algo (NcmRNG *rng, const gchar *algo);
void ncm_rng_set_state (NcmRNG *rng, const gchar *state);
gboolean ncm_rng_check_seed (NcmRNG *rng, gulong seed);
void ncm_rng_set_seed (NcmRNG *rng, gulong seed);
gulong ncm_rng_get_seed (NcmRNG *rng);
void ncm_rng_set_random_seed (NcmRNG *rng, gboolean allow_colisions);

NcmRNG *ncm_rng_pool_get (const gchar *name);

NCM_INLINE gdouble ncm_rng_uniform_gen (NcmRNG *rng, const gdouble xl, const gdouble xu); 
NCM_INLINE gdouble ncm_rng_gaussian_gen (NcmRNG *rng, const gdouble mu, const gdouble sigma); 
NCM_INLINE gdouble ncm_rng_ugaussian_gen (NcmRNG *rng); 
NCM_INLINE gdouble ncm_rng_gaussian_tail_gen (NcmRNG *rng, const gdouble a, const gdouble sigma); 
NCM_INLINE gdouble ncm_rng_exponential_gen (NcmRNG *rng, const gdouble mu);
NCM_INLINE gdouble ncm_rng_laplace_gen (NcmRNG *rng, const gdouble a);
NCM_INLINE gdouble ncm_rng_exppow_gen (NcmRNG *rng, const gdouble a, const gdouble b);
NCM_INLINE gdouble ncm_rng_beta_gen (NcmRNG *rng, const gdouble a, const gdouble b);
NCM_INLINE gdouble ncm_rng_gamma_gen (NcmRNG *rng, const gdouble a, const gdouble b);

G_END_DECLS

#endif /* _NCM_RNG_H_ */

#ifndef _NCM_RNG_INLINE_H_
#define _NCM_RNG_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gdouble ncm_rng_uniform_gen (NcmRNG *rng, const gdouble xl, const gdouble xu)         { return gsl_ran_flat (rng->r, xl, xu);            }
NCM_INLINE gdouble ncm_rng_gaussian_gen (NcmRNG *rng, const gdouble mu, const gdouble sigma)     { return gsl_ran_gaussian (rng->r, sigma) + mu;    } 
NCM_INLINE gdouble ncm_rng_ugaussian_gen (NcmRNG *rng)                                           { return gsl_ran_ugaussian (rng->r);               } 
NCM_INLINE gdouble ncm_rng_gaussian_tail_gen (NcmRNG *rng, const gdouble a, const gdouble sigma) { return gsl_ran_gaussian_tail (rng->r, a, sigma); } 
NCM_INLINE gdouble ncm_rng_exponential_gen (NcmRNG *rng, const gdouble mu)                       { return gsl_ran_exponential (rng->r, mu);         }
NCM_INLINE gdouble ncm_rng_laplace_gen (NcmRNG *rng, const gdouble a)                            { return gsl_ran_laplace (rng->r, a);              }
NCM_INLINE gdouble ncm_rng_exppow_gen (NcmRNG *rng, const gdouble a, const gdouble b)            { return gsl_ran_exppow (rng->r, a, b);            }
NCM_INLINE gdouble ncm_rng_beta_gen (NcmRNG *rng, const gdouble a, const gdouble b)              { return gsl_ran_beta (rng->r, a, b);              } 
NCM_INLINE gdouble ncm_rng_gamma_gen (NcmRNG *rng, const gdouble a, const gdouble b)             { return gsl_ran_gamma (rng->r, a, b);             } 

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_RNG_INLINE_H_ */

