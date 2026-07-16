/***************************************************************************
 *            nc_galaxy_shape_factor_laplace.c
 *
 *  Fri Jul 3 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_laplace.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyShapeFactorLaplace:
 *
 * Laplace-approximate evaluation of the intrinsic-ellipticity marginal.
 *
 * Approximates
 * $$P(\epsilon_\mathrm{obs} \mid g) = \int_{|\chi_I|<1} \mathrm{d}^2\chi_I\,
 *   P_\mathrm{pop}(\chi_I)\, N_2\big(\epsilon_\mathrm{obs} - f_g(\chi_I);
 *   \sigma_\mathrm{noise}^2\big)$$
 * by a single Gaussian expansion around the joint mode of the integrand,
 * found by nc_galaxy_shape_intrinsic_mode_find() (a nested search: since
 * $P_\mathrm{pop}$ has no angular dependence, the whole $\theta$-direction is
 * resolved by the noise term alone at any fixed $\rho$, so only the radial
 * direction needs profiling against $P_\mathrm{pop}$ -- see that file's
 * documentation for why this beats a blind joint 2D minimization). The
 * result is the standard Laplace formula,
 * $$P(\epsilon_\mathrm{obs}\mid g) \approx \frac{2\pi}{\sqrt{\det(-H)}}\,
 *   \exp\big(\ln P_\mathrm{pop}(\chi_{I,\star})+\ln
 *   N_2(\epsilon_\mathrm{obs}-f_g(\chi_{I,\star}))\big),$$
 * with $H$ the Hessian of the log-integrand in Cartesian $\chi_I$ coordinates
 * at the mode $\chi_{I,\star}$.
 *
 * This is a single closed-form evaluation (a handful of Newton steps plus
 * one Hessian, all on cheap arithmetic -- no cubature), orders of magnitude
 * cheaper than #NcGalaxyShapeFactorQuad's ~100ms/evaluation Divonne call, and
 * unlike #NcGalaxyShapeFactorVarAdd it needs no
 * nc_galaxy_shape_pop_get_sigma() support, so it works for any population
 * (e.g. #NcGalaxyShapePopBeta) that implements eval_p().
 *
 * The approximation is exact in the limit of an infinitely sharp,
 * single-peaked integrand and degrades as the population/noise combination
 * gets broader or flatter (verified against #NcGalaxyShapeFactorQuad: from
 * sub-percent for concentrated populations to a genuine few percent for a
 * broad, barely-peaked one, e.g. #NcGalaxyShapePopBeta with a small
 * concentration parameter). nc_galaxy_shape_factor_laplace_eval_marginal()
 * returns NAN (propagating to eval_ln_marginal() as well) when the found
 * point is not a proper local maximum (a non-positive-definite Hessian),
 * which signals that a single local Gaussian is not a meaningful description
 * of the integrand there -- callers needing a guaranteed answer in that
 * regime should fall back to #NcGalaxyShapeFactorQuad.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_laplace.h"
#include "nc/lss/galaxy/nc_galaxy_shape_intrinsic_mode.h"
#include "nc/lss/wl/nc_wl_ellipticity.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyShapeFactorLaplace
{
  NcGalaxyShapeFactor parent_instance;
};

typedef struct _NcGalaxyShapeFactorLaplacePrivate
{
  /* Resolved once at construction, see #NcGalaxyShapeFactorQuad's private
   * struct for why. */
  complex double (*apply_shear) (complex double g, complex double chi);
  complex double (*apply_shear_inv) (complex double g, complex double eps_obs);

  /* TRACE_DET has a closed-form mode-finder (no apply_shear calls in its
   * search loop, see nc_galaxy_shape_intrinsic_mode.c); TRACE still uses
   * the generic finite-difference one. Dispatched once here, not per
   * galaxy -- same "specialize and choose once" convention as apply_shear
   * itself. */
  gboolean use_closed_form_trace_det;

  /* TEMPORARY, for an A/B experiment only -- not a shipped feature, remove
   * after comparing production results with/without the closed-form TRACE
   * mode-finder. Gated on an env var so both arms can be run from the same
   * build. */
  gboolean use_closed_form_trace;
} NcGalaxyShapeFactorLaplacePrivate;

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorLaplace, nc_galaxy_shape_factor_laplace, NC_TYPE_GALAXY_SHAPE_FACTOR)

static void
_nc_galaxy_shape_factor_laplace_ldata_noop (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_nc_galaxy_shape_factor_laplace_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_factor_laplace_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  /* No persistent per-galaxy state: each evaluation is a single, self-
   * contained mode search. */
  data->ldata                  = NULL;
  data->ldata_destroy          = &g_free;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_laplace_ldata_noop;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_laplace_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_laplace_ldata_required_columns;
}

static void
_nc_galaxy_shape_factor_laplace_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset)
{
}

static gdouble
_nc_galaxy_shape_factor_laplace_eval (NcGalaxyShapeFactorLaplace *gsfl, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorLaplacePrivate * const self = nc_galaxy_shape_factor_laplace_get_instance_private (gsfl);
  const complex double g                         = g_1 + I * g_2;
  const complex double eps_obs                   = epsilon_obs_1 + I * epsilon_obs_2;
  NcGalaxyShapeIntrinsicMode mode;

  if (self->use_closed_form_trace_det)
    nc_galaxy_shape_intrinsic_mode_find_trace_det (pop, data->pop_data, g, eps_obs, data->std_noise, &mode);
  else if (self->use_closed_form_trace)
    nc_galaxy_shape_intrinsic_mode_find_trace (pop, data->pop_data, g, eps_obs, data->std_noise, &mode);
  else
    nc_galaxy_shape_intrinsic_mode_find (self->apply_shear, self->apply_shear_inv,
                                         pop, data->pop_data, g, eps_obs, data->std_noise, &mode);

  return nc_galaxy_shape_intrinsic_mode_laplace (&mode);
}

static gdouble
_nc_galaxy_shape_factor_laplace_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return _nc_galaxy_shape_factor_laplace_eval (NC_GALAXY_SHAPE_FACTOR_LAPLACE (gsf), pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2);
}

static gdouble
_nc_galaxy_shape_factor_laplace_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return log (_nc_galaxy_shape_factor_laplace_eval (NC_GALAXY_SHAPE_FACTOR_LAPLACE (gsf), pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2));
}

/* GObject boilerplate --------------------------------------------------- */

static void
nc_galaxy_shape_factor_laplace_init (NcGalaxyShapeFactorLaplace *gsfl)
{
  NcGalaxyShapeFactorLaplacePrivate * const self = nc_galaxy_shape_factor_laplace_get_instance_private (gsfl);

  self->apply_shear               = NULL;
  self->apply_shear_inv           = NULL;
  self->use_closed_form_trace_det = FALSE;
  self->use_closed_form_trace     = FALSE;
}

static void
_nc_galaxy_shape_factor_laplace_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_laplace_parent_class)->constructed (object);
  {
    NcGalaxyShapeFactorLaplace *gsfl               = NC_GALAXY_SHAPE_FACTOR_LAPLACE (object);
    NcGalaxyShapeFactorLaplacePrivate * const self = nc_galaxy_shape_factor_laplace_get_instance_private (gsfl);
    const NcGalaxyWLObsEllipConv ellip_conv        = nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (gsfl));

    switch (ellip_conv)
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
        self->apply_shear           = &nc_wl_ellipticity_apply_shear_trace;
        self->apply_shear_inv       = &nc_wl_ellipticity_apply_shear_inv_trace;
        self->use_closed_form_trace = (g_getenv ("NC_TRACE_CLOSED_FORM") != NULL);
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        self->apply_shear               = &nc_wl_ellipticity_apply_shear_trace_det;
        self->apply_shear_inv           = &nc_wl_ellipticity_apply_shear_inv_trace_det;
        self->use_closed_form_trace_det = TRUE;
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
    }
  }
}

static void
_nc_galaxy_shape_factor_laplace_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_laplace_parent_class)->finalize (object);
}

static void
nc_galaxy_shape_factor_laplace_class_init (NcGalaxyShapeFactorLaplaceClass *klass)
{
  NcGalaxyShapeFactorClass *gsf_class = NC_GALAXY_SHAPE_FACTOR_CLASS (klass);
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);

  object_class->constructed = &_nc_galaxy_shape_factor_laplace_constructed;
  object_class->finalize    = &_nc_galaxy_shape_factor_laplace_finalize;

  gsf_class->data_init        = &_nc_galaxy_shape_factor_laplace_data_init;
  gsf_class->prepare          = &_nc_galaxy_shape_factor_laplace_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_laplace_eval_marginal;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_laplace_eval_ln_marginal;
}

/**
 * nc_galaxy_shape_factor_laplace_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 *
 * Creates a new #NcGalaxyShapeFactorLaplace.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorLaplace.
 */
NcGalaxyShapeFactorLaplace *
nc_galaxy_shape_factor_laplace_new (NcGalaxyWLObsEllipConv ellip_conv)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_LAPLACE,
                       "ellip-conv", ellip_conv,
                       NULL);
}

/**
 * nc_galaxy_shape_factor_laplace_ref:
 * @gsfl: a #NcGalaxyShapeFactorLaplace
 *
 * Increases the reference count of @gsfl by one.
 *
 * Returns: (transfer full): @gsfl.
 */
NcGalaxyShapeFactorLaplace *
nc_galaxy_shape_factor_laplace_ref (NcGalaxyShapeFactorLaplace *gsfl)
{
  return g_object_ref (gsfl);
}

/**
 * nc_galaxy_shape_factor_laplace_free:
 * @gsfl: a #NcGalaxyShapeFactorLaplace
 *
 * Decreases the reference count of @gsfl by one.
 *
 */
void
nc_galaxy_shape_factor_laplace_free (NcGalaxyShapeFactorLaplace *gsfl)
{
  g_object_unref (gsfl);
}

/**
 * nc_galaxy_shape_factor_laplace_clear:
 * @gsfl: a #NcGalaxyShapeFactorLaplace
 *
 * Decreases the reference count of *@gsfl by one, and sets the pointer
 * *@gsfl to NULL.
 *
 */
void
nc_galaxy_shape_factor_laplace_clear (NcGalaxyShapeFactorLaplace **gsfl)
{
  g_clear_object (gsfl);
}

