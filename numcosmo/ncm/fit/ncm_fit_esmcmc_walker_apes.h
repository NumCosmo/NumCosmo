/***************************************************************************
 *            ncm_fit_esmcmc_walker_apes.h
 *
 *  Sat October 27 13:08:34 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker_apes.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifndef _NCM_FIT_ESMCMC_WALKER_APES_H_
#define _NCM_FIT_ESMCMC_WALKER_APES_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/fit/ncm_fit_esmcmc_walker.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/ncm/stats/ncm_stats_dist.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_ESMCMC_WALKER_APES (ncm_fit_esmcmc_walker_apes_get_type ())

G_DECLARE_FINAL_TYPE (NcmFitESMCMCWalkerAPES, ncm_fit_esmcmc_walker_apes, NCM, FIT_ESMCMC_WALKER_APES, NcmFitESMCMCWalker)

/**
 * NcmFitESMCMCWalkerAPESMethod:
 * @NCM_FIT_ESMCMC_WALKER_APES_METHOD_KDE: Fixed kernel estimation.
 * @NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE: Variable kernel estimation.
 *
 * Posterior estimation method.
 *
 */
typedef enum _NcmFitESMCMCWalkerAPESMethod /*< prefix=NCM_FIT_ESMCMC_WALKER_APES_METHOD >*/
{
  NCM_FIT_ESMCMC_WALKER_APES_METHOD_KDE = 0,
  NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE,
  /* < private > */
  NCM_FIT_ESMCMC_WALKER_APES_METHOD_LEN, /*< skip >*/
} NcmFitESMCMCWalkerAPESMethod;

/**
 * NcmFitESMCMCWalkerAPESKType:
 * @NCM_FIT_ESMCMC_WALKER_APES_KTYPE_CAUCHY: Cauchy kernel.
 * @NCM_FIT_ESMCMC_WALKER_APES_KTYPE_ST3: Student-t kernel with $\nu=3$.
 * @NCM_FIT_ESMCMC_WALKER_APES_KTYPE_GAUSS: Gaussian kernel.
 * @NCM_FIT_ESMCMC_WALKER_APES_KTYPE_AUTO: Student-t kernel whose degrees of freedom are
 * fitted together with the bandwidth, by the same out-of-sample objective; the Gaussian
 * kernel at the upper bound. Needs a #NcmFitESMCMCWalkerAPES:cv-type that fits the
 * bandwidth.
 *
 * Kernel used to build the approximate posterior.
 *
 */
typedef enum _NcmFitESMCMCWalkerAPESKType /*< enum,underscore_name=NCM_FIT_ESMCMC_WALKER_APES_KTYPE,prefix=NCM_FIT_ESMCMC_WALKER_APES_KTYPE >*/
{
  NCM_FIT_ESMCMC_WALKER_APES_KTYPE_CAUCHY = 0,
  NCM_FIT_ESMCMC_WALKER_APES_KTYPE_ST3,
  NCM_FIT_ESMCMC_WALKER_APES_KTYPE_GAUSS,
  NCM_FIT_ESMCMC_WALKER_APES_KTYPE_AUTO,
  /* < private > */
  NCM_FIT_ESMCMC_WALKER_APES_KTYPE_LEN, /*< skip >*/
} NcmFitESMCMCWalkerAPESKType;

NcmFitESMCMCWalkerAPES *ncm_fit_esmcmc_walker_apes_new (guint nwalkers, guint nparams);
NcmFitESMCMCWalkerAPES *ncm_fit_esmcmc_walker_apes_new_full (guint nwalkers, guint nparams, NcmFitESMCMCWalkerAPESMethod method, NcmFitESMCMCWalkerAPESKType k_type, gdouble over_smooth);
NcmFitESMCMCWalkerAPES *ncm_fit_esmcmc_walker_apes_ref (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_free (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_clear (NcmFitESMCMCWalkerAPES **apes);

void ncm_fit_esmcmc_walker_apes_set_method (NcmFitESMCMCWalkerAPES *apes, NcmFitESMCMCWalkerAPESMethod method);
void ncm_fit_esmcmc_walker_apes_set_k_type (NcmFitESMCMCWalkerAPES *apes, NcmFitESMCMCWalkerAPESKType k_type);
void ncm_fit_esmcmc_walker_apes_set_over_smooth (NcmFitESMCMCWalkerAPES *apes, const gdouble os);

NcmFitESMCMCWalkerAPESMethod ncm_fit_esmcmc_walker_apes_get_method (NcmFitESMCMCWalkerAPES *apes);
NcmFitESMCMCWalkerAPESKType ncm_fit_esmcmc_walker_apes_get_k_type (NcmFitESMCMCWalkerAPES *apes);
gdouble ncm_fit_esmcmc_walker_apes_get_over_smooth (NcmFitESMCMCWalkerAPES *apes);

void ncm_fit_esmcmc_walker_apes_set_use_threads (NcmFitESMCMCWalkerAPES *apes, gboolean use_threads);
gboolean ncm_fit_esmcmc_walker_apes_get_use_threads (NcmFitESMCMCWalkerAPES *apes);

void ncm_fit_esmcmc_walker_apes_set_center_shrink (NcmFitESMCMCWalkerAPES *apes, gboolean center_shrink);
gboolean ncm_fit_esmcmc_walker_apes_get_center_shrink (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_set_defensive_frac (NcmFitESMCMCWalkerAPES *apes, const gdouble frac);
gdouble ncm_fit_esmcmc_walker_apes_get_defensive_frac (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_set_defensive_scale (NcmFitESMCMCWalkerAPES *apes, const gdouble scale);
gdouble ncm_fit_esmcmc_walker_apes_get_defensive_scale (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_set_defensive_nu (NcmFitESMCMCWalkerAPES *apes, const gdouble nu);
gdouble ncm_fit_esmcmc_walker_apes_get_defensive_nu (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_set_vkde_points_per_dim (NcmFitESMCMCWalkerAPES *apes, const gdouble points_per_dim);
gdouble ncm_fit_esmcmc_walker_apes_get_vkde_points_per_dim (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_set_uniform_weights (NcmFitESMCMCWalkerAPES *apes, gboolean uniform_weights);
gboolean ncm_fit_esmcmc_walker_apes_get_uniform_weights (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_set_cv_type (NcmFitESMCMCWalkerAPES *apes, NcmStatsDistCV cv_type);
NcmStatsDistCV ncm_fit_esmcmc_walker_apes_get_cv_type (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_set_split_frac (NcmFitESMCMCWalkerAPES *apes, const gdouble split_frac);
gdouble ncm_fit_esmcmc_walker_apes_get_split_frac (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_set_auto_kernel (NcmFitESMCMCWalkerAPES *apes, gboolean auto_kernel);
gboolean ncm_fit_esmcmc_walker_apes_get_auto_kernel (NcmFitESMCMCWalkerAPES *apes);

void ncm_fit_esmcmc_walker_apes_peek_sds (NcmFitESMCMCWalkerAPES *apes, NcmStatsDist **sd0, NcmStatsDist **sd1);

void ncm_fit_esmcmc_walker_apes_set_local_frac (NcmFitESMCMCWalkerAPES *apes, gdouble local_frac);
void ncm_fit_esmcmc_walker_apes_set_cov_fixed_from_mset (NcmFitESMCMCWalkerAPES *apes, NcmMSet *mset);
void ncm_fit_esmcmc_walker_apes_set_cov_robust_diag (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_set_cov_robust (NcmFitESMCMCWalkerAPES *apes);

void ncm_fit_esmcmc_walker_apes_set_exploration (NcmFitESMCMCWalkerAPES *apes, guint exploration);
guint ncm_fit_esmcmc_walker_apes_get_exploration (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_set_exploration_qratio_floor (NcmFitESMCMCWalkerAPES *apes, const gdouble qratio_floor);
gdouble ncm_fit_esmcmc_walker_apes_get_exploration_qratio_floor (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_set_exploration_stop_after (NcmFitESMCMCWalkerAPES *apes, guint stop_after);
guint ncm_fit_esmcmc_walker_apes_get_exploration_stop_after (NcmFitESMCMCWalkerAPES *apes);
gboolean ncm_fit_esmcmc_walker_apes_is_exploring (NcmFitESMCMCWalkerAPES *apes);

G_END_DECLS

#endif /* _NCM_FIT_ESMCMC_WALKER_APES_H_ */

