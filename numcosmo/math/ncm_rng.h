/***************************************************************************
 *            ncm_rng.h
 *
 *  Sat August 17 12:39:57 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * ncm_rng.h
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

G_BEGIN_DECLS

#define NCM_TYPE_RNG (ncm_rng_get_type ())
#define NCM_TYPE_RNG_DISCRETE (ncm_rng_discrete_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmRNG, ncm_rng, NCM, RNG, GObject)
typedef struct _NcmRNGDiscrete NcmRNGDiscrete;

struct _NcmRNGClass
{
  /*< private >*/
  GObjectClass parent_class;
  GRand *seed_gen;
  GHashTable *seed_hash;

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[16];
};

GType ncm_rng_discrete_get_type (void) G_GNUC_CONST;

NcmRNGDiscrete *ncm_rng_discrete_new (const gdouble *weights, const guint n);
NcmRNGDiscrete *ncm_rng_discrete_copy (NcmRNGDiscrete *rng);
void ncm_rng_discrete_free (NcmRNGDiscrete *rng);

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

gulong ncm_rng_gen_ulong (NcmRNG *rng);
gulong ncm_rng_uniform_int_gen (NcmRNG *rng, gulong n);
gdouble ncm_rng_uniform01_gen (NcmRNG *rng);
gdouble ncm_rng_uniform01_pos_gen (NcmRNG *rng);
gdouble ncm_rng_uniform_gen (NcmRNG *rng, const gdouble xl, const gdouble xu);
gdouble ncm_rng_gaussian_gen (NcmRNG *rng, const gdouble mu, const gdouble sigma);
gdouble ncm_rng_ugaussian_gen (NcmRNG *rng);
gdouble ncm_rng_gaussian_tail_gen (NcmRNG *rng, const gdouble a, const gdouble sigma);
gdouble ncm_rng_exponential_gen (NcmRNG *rng, const gdouble mu);
gdouble ncm_rng_laplace_gen (NcmRNG *rng, const gdouble a);
gdouble ncm_rng_exppow_gen (NcmRNG *rng, const gdouble a, const gdouble b);
gdouble ncm_rng_beta_gen (NcmRNG *rng, const gdouble a, const gdouble b);
gdouble ncm_rng_gamma_gen (NcmRNG *rng, const gdouble a, const gdouble b);
gdouble ncm_rng_chisq_gen (NcmRNG *rng, const gdouble nu);
gdouble ncm_rng_poisson_gen (NcmRNG *rng, const gdouble mu);
gdouble ncm_rng_rayleigh_gen (NcmRNG *rng, const gdouble sigma);
gsize ncm_rng_discrete_gen (NcmRNG *rng, NcmRNGDiscrete *rng_discrete);
void ncm_rng_sample (NcmRNG *rng, void *dest, size_t k, void *src, size_t n, size_t size);
void ncm_rng_choose (NcmRNG *rng, void *dest, size_t k, void *src, size_t n, size_t size);
void ncm_rng_multinomial (NcmRNG *rng, gsize K, guint N, const gdouble *p, guint *n);
void ncm_rng_bivariate_gaussian_gen (NcmRNG *rng, const gdouble sigma_x, const gdouble sigma_y, const gdouble rho, gdouble *x, gdouble *y);

G_END_DECLS

#endif /* _NCM_RNG_H_ */

