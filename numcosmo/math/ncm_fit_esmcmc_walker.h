/***************************************************************************
 *            ncm_fit_esmcmc_walker.h
 *
 *  Wed March 16 13:07:20 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_FIT_ESMCMC_WALKER_H_
#define _NCM_FIT_ESMCMC_WALKER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_rng.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_ESMCMC_WALKER (ncm_fit_esmcmc_walker_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmFitESMCMCWalker, ncm_fit_esmcmc_walker, NCM, FIT_ESMCMC_WALKER, GObject)

struct _NcmFitESMCMCWalkerClass
{
  /*< private >*/
  GObjectClass parent_class;

  void (*set_size) (NcmFitESMCMCWalker *walker, guint size);
  guint (*get_size) (NcmFitESMCMCWalker *walker);
  void (*set_nparams) (NcmFitESMCMCWalker *walker, guint nparams);
  guint (*get_nparams) (NcmFitESMCMCWalker *walker);
  void (*setup) (NcmFitESMCMCWalker *walker, NcmMSet *mset, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng);
  void (*step) (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
  gdouble (*prob) (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star);
  gdouble (*prob_norm) (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
  void (*clean) (NcmFitESMCMCWalker *walker, guint ki, guint kf);
  const gchar *(*desc) (NcmFitESMCMCWalker *walker);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[8];
};

NcmFitESMCMCWalker *ncm_fit_esmcmc_walker_ref (NcmFitESMCMCWalker *walker);
void ncm_fit_esmcmc_walker_free (NcmFitESMCMCWalker *walker);
void ncm_fit_esmcmc_walker_clear (NcmFitESMCMCWalker **walker);

void ncm_fit_esmcmc_walker_set_size (NcmFitESMCMCWalker *walker, guint size);
guint ncm_fit_esmcmc_walker_get_size (NcmFitESMCMCWalker *walker);
void ncm_fit_esmcmc_walker_set_nparams (NcmFitESMCMCWalker *walker, guint nparams);
guint ncm_fit_esmcmc_walker_get_nparams (NcmFitESMCMCWalker *walker);

void ncm_fit_esmcmc_walker_setup (NcmFitESMCMCWalker *walker, NcmMSet *mset, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng);
void ncm_fit_esmcmc_walker_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
gdouble ncm_fit_esmcmc_walker_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star);
gdouble ncm_fit_esmcmc_walker_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
void ncm_fit_esmcmc_walker_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf);
const gchar *ncm_fit_esmcmc_walker_desc (NcmFitESMCMCWalker *walker);

G_END_DECLS

#endif /* _NCM_FIT_ESMCMC_WALKER_H_ */

