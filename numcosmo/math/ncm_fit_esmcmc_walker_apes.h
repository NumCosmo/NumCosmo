/***************************************************************************
 *            ncm_fit_esmcmc_walker_apes.h
 *
 *  Sat October 27 13:08:34 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker_apes.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
#include <numcosmo/math/ncm_fit_esmcmc_walker.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_stats_dist.h>

G_BEGIN_DECLS

#define NCM_TYPE_FIT_ESMCMC_WALKER_APES             (ncm_fit_esmcmc_walker_apes_get_type ())
#define NCM_FIT_ESMCMC_WALKER_APES(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FIT_ESMCMC_WALKER_APES, NcmFitESMCMCWalkerAPES))
#define NCM_FIT_ESMCMC_WALKER_APES_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FIT_ESMCMC_WALKER_APES, NcmFitESMCMCWalkerAPESClass))
#define NCM_IS_FIT_ESMCMC_WALKER_APES(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FIT_ESMCMC_WALKER_APES))
#define NCM_IS_FIT_ESMCMC_WALKER_APES_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FIT_ESMCMC_WALKER_APES))
#define NCM_FIT_ESMCMC_WALKER_APES_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FIT_ESMCMC_WALKER_APES, NcmFitESMCMCWalkerAPESClass))

typedef struct _NcmFitESMCMCWalkerAPESClass NcmFitESMCMCWalkerAPESClass;
typedef struct _NcmFitESMCMCWalkerAPES NcmFitESMCMCWalkerAPES;
typedef struct _NcmFitESMCMCWalkerAPESPrivate NcmFitESMCMCWalkerAPESPrivate;

struct _NcmFitESMCMCWalkerAPESClass
{
  /*< private >*/
  NcmFitESMCMCWalkerClass parent_class;
};

struct _NcmFitESMCMCWalkerAPES
{
  /*< private >*/
  NcmFitESMCMCWalker parent_instance;
  NcmFitESMCMCWalkerAPESPrivate *priv;
};

/**
 * NcmFitESMCMCWalkerAPESMethod:
 * @NCM_FIT_ESMCMC_WALKER_APES_METHOD_KDE: Fixed kernel estimation.
 * @NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE: Variable kernel estimation.
 *
 * Posterior estimation method.
 *
 */
typedef enum _NcmFitESMCMCWalkerAPESMethod
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
 *
 * Kernel used to build the approximate posterior.
 *
 */
typedef enum _NcmFitESMCMCWalkerAPESKType /*< enum,underscore_name=NCM_FIT_ESMCMC_WALKER_APES_KTYPE >*/
{
  NCM_FIT_ESMCMC_WALKER_APES_KTYPE_CAUCHY = 0,
  NCM_FIT_ESMCMC_WALKER_APES_KTYPE_ST3,
  NCM_FIT_ESMCMC_WALKER_APES_KTYPE_GAUSS,
  /* < private > */
  NCM_FIT_ESMCMC_WALKER_APES_KTYPE_LEN, /*< skip >*/
} NcmFitESMCMCWalkerAPESKType;

GType ncm_fit_esmcmc_walker_apes_get_type (void) G_GNUC_CONST;

NcmFitESMCMCWalkerAPES *ncm_fit_esmcmc_walker_apes_new (guint nwalkers, guint nparams);
NcmFitESMCMCWalkerAPES *ncm_fit_esmcmc_walker_apes_new_full (guint nwalkers, guint nparams, NcmFitESMCMCWalkerAPESMethod method, NcmFitESMCMCWalkerAPESKType k_type, gdouble over_smooth, gboolean use_interp);
NcmFitESMCMCWalkerAPES *ncm_fit_esmcmc_walker_apes_ref (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_free (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_clear (NcmFitESMCMCWalkerAPES **apes);

void ncm_fit_esmcmc_walker_apes_set_method (NcmFitESMCMCWalkerAPES *apes, NcmFitESMCMCWalkerAPESMethod method);
void ncm_fit_esmcmc_walker_apes_set_k_type (NcmFitESMCMCWalkerAPES *apes, NcmFitESMCMCWalkerAPESKType k_type);
void ncm_fit_esmcmc_walker_apes_set_over_smooth (NcmFitESMCMCWalkerAPES *apes, const gdouble os);

NcmFitESMCMCWalkerAPESMethod ncm_fit_esmcmc_walker_apes_get_method (NcmFitESMCMCWalkerAPES *apes);
NcmFitESMCMCWalkerAPESKType ncm_fit_esmcmc_walker_apes_get_k_type (NcmFitESMCMCWalkerAPES *apes);
gdouble ncm_fit_esmcmc_walker_apes_get_over_smooth (NcmFitESMCMCWalkerAPES *apes);

void ncm_fit_esmcmc_walker_apes_use_interp (NcmFitESMCMCWalkerAPES *apes, gboolean use_interp);
gboolean ncm_fit_esmcmc_walker_apes_interp (NcmFitESMCMCWalkerAPES *apes);

void ncm_fit_esmcmc_walker_apes_peek_sds (NcmFitESMCMCWalkerAPES *apes, NcmStatsDist **sd0, NcmStatsDist **sd1);

void ncm_fit_esmcmc_walker_apes_set_local_frac (NcmFitESMCMCWalkerAPES *apes, gdouble local_frac);
void ncm_fit_esmcmc_walker_apes_set_cov_fixed_from_mset (NcmFitESMCMCWalkerAPES *apes, NcmMSet *mset);
void ncm_fit_esmcmc_walker_apes_set_cov_robust_diag (NcmFitESMCMCWalkerAPES *apes);
void ncm_fit_esmcmc_walker_apes_set_cov_robust (NcmFitESMCMCWalkerAPES *apes);

G_END_DECLS

#endif /* _NCM_FIT_ESMCMC_WALKER_APES_H_ */

