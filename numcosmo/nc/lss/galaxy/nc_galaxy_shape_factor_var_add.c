/***************************************************************************
 *            nc_galaxy_shape_factor_var_add.c
 *
 *  Thu Jul 2 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_var_add.c
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
 * NcGalaxyShapeFactorVarAdd:
 *
 * Variance-add evaluation of the intrinsic-ellipticity marginal.
 *
 * This is the classic (approximate) Gaussian treatment: instead of
 * marginalizing over $\chi_I$, the observed ellipticity is pulled back through
 * the inverse shear map, $\chi_I^\ast = f_{\tilde g}^{-1}(\epsilon_\mathrm{obs})$,
 * and evaluated on a Gaussian whose variance is the sum of the intrinsic and
 * noise variances,
 * $$P(\epsilon_\mathrm{obs} \mid g) \approx
 *   \frac{e^{-|\chi_I^\ast|^2 / 2(\sigma^2 + \sigma_\mathrm{noise}^2)}}
 *        {2\pi\,(\sigma^2 + \sigma_\mathrm{noise}^2)}\,
 *   \left|\det J_{f^{-1}}(\epsilon_\mathrm{obs})\right|.$$
 * The "variance-add" combination itself (rather than the exact truncated
 * convolution) is exact only in the doubly-linear regime where the shear map
 * is expanded jointly to first order in $g$ AND $\chi_I$ (not just in $g$
 * alone, since the dropped terms mix the two) and the intrinsic integral is
 * extended from the physical unit disc $|\chi_I|<1$ to the whole plane (the
 * only way two independent Gaussians combine in closed form). This
 * implementation removes the map-linearization by using the exact inverse
 * map and its exact Jacobian at a single pulled-back point, but keeps the
 * plane-instead-of-disc approximation (the plain un-truncated Gaussian
 * normalization), matching the legacy implementation's behaviour bit for
 * bit. See the
 * <a href="../../theory/wl_ellipticity.html#the-variance-add-approximation">Variance-Add Approximation</a>
 * section of the theory page for the full derivation.
 *
 * The variance addition is only meaningful for a population parameterized by
 * an (untruncated) Gaussian width sigma, so this method requires the
 * population resolved from the #NcmMSet to support
 * nc_galaxy_shape_pop_get_sigma() (currently #NcGalaxyShapePopGauss or
 * #NcGalaxyShapePopGaussLocal, Global or per-galaxy).
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_var_add.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyShapeFactorVarAdd
{
  NcGalaxyShapeFactor parent_instance;
};

typedef struct _NcGalaxyShapeFactorVarAddPrivate
{
  /* Convention specialization resolved once at construction (ellip-conv is
   * CONSTRUCT_ONLY), keeping the switch out of the per-evaluation path. */
  complex double (*apply_shear_inv) (complex double g, complex double e_obs);

  gdouble (*lndet_jac) (complex double g, complex double e_obs);
} NcGalaxyShapeFactorVarAddPrivate;

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyShapeFactorVarAdd, nc_galaxy_shape_factor_var_add, NC_TYPE_GALAXY_SHAPE_FACTOR)

static void
nc_galaxy_shape_factor_var_add_init (NcGalaxyShapeFactorVarAdd *gsfva)
{
  NcGalaxyShapeFactorVarAddPrivate * const self = nc_galaxy_shape_factor_var_add_get_instance_private (gsfva);

  self->apply_shear_inv = NULL;
  self->lndet_jac       = NULL;
}

static void
_nc_galaxy_shape_factor_var_add_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_var_add_parent_class)->constructed (object);
  {
    NcGalaxyShapeFactorVarAdd *gsfva              = NC_GALAXY_SHAPE_FACTOR_VAR_ADD (object);
    NcGalaxyShapeFactorVarAddPrivate * const self = nc_galaxy_shape_factor_var_add_get_instance_private (gsfva);
    const NcGalaxyWLObsEllipConv ellip_conv       = nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (gsfva));

    switch (ellip_conv)
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
        self->apply_shear_inv = &nc_wl_ellipticity_apply_shear_inv_trace_c;
        self->lndet_jac       = &nc_wl_ellipticity_lndet_jac_trace_c;
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        self->apply_shear_inv = &nc_wl_ellipticity_apply_shear_inv_trace_det_c;
        self->lndet_jac       = &nc_wl_ellipticity_lndet_jac_trace_det_c;
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
    }
  }
}

static void
_nc_galaxy_shape_factor_var_add_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_var_add_parent_class)->finalize (object);
}

typedef struct _NcGalaxyShapeFactorVarAddData
{
  gint placeholder;
} NcGalaxyShapeFactorVarAddData;

static void
_nc_galaxy_shape_factor_var_add_ldata_noop (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_nc_galaxy_shape_factor_var_add_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_factor_var_add_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorVarAddData *ldata = g_new0 (NcGalaxyShapeFactorVarAddData, 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &g_free;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_var_add_ldata_noop;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_var_add_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_var_add_ldata_required_columns;
}

static void
_nc_galaxy_shape_factor_var_add_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  if (data->pop_data->ldata_get_sigma == NULL)
  {
    NcGalaxyShapePop *pop = NC_GALAXY_SHAPE_POP (ncm_mset_peek (mset, nc_galaxy_shape_pop_id ()));

    g_error ("NcGalaxyShapeFactorVarAdd: the variance-add approximation requires a population "
             "supporting nc_galaxy_shape_pop_get_sigma(), got %s.", G_OBJECT_TYPE_NAME (pop));
  }
}

static inline void
_nc_galaxy_shape_factor_var_add_chi2 (NcGalaxyShapeFactorVarAddPrivate * const self, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data,
                                      const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2,
                                      gdouble *total_var, gdouble *chi2, gdouble *lndetjac)
{
  const gdouble sigma      = nc_galaxy_shape_pop_get_sigma (pop, data->pop_data);
  const complex double g   = g_1 + I * g_2;
  const complex double e_o = epsilon_obs_1 + I * epsilon_obs_2;
  const complex double e_s = self->apply_shear_inv (g, e_o);

  /* The two components are reduced separately (chi2_1 + chi2_2, not
   * (e1^2 + e2^2)/var) to stay bit-identical to the legacy implementation. */
  *total_var = gsl_pow_2 (sigma) + gsl_pow_2 (data->std_noise);
  *chi2      = gsl_pow_2 (creal (e_s)) / *total_var + gsl_pow_2 (cimag (e_s)) / *total_var;
  *lndetjac  = self->lndet_jac (g, e_o);
}

static gdouble
_nc_galaxy_shape_factor_var_add_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorVarAddPrivate * const self = nc_galaxy_shape_factor_var_add_get_instance_private (NC_GALAXY_SHAPE_FACTOR_VAR_ADD (gsf));
  gdouble total_var, chi2, lndetjac;

  _nc_galaxy_shape_factor_var_add_chi2 (self, pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2, &total_var, &chi2, &lndetjac);

  return exp (-0.5 * chi2 + lndetjac) / (2.0 * M_PI * total_var);
}

static gdouble
_nc_galaxy_shape_factor_var_add_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorVarAddPrivate * const self = nc_galaxy_shape_factor_var_add_get_instance_private (NC_GALAXY_SHAPE_FACTOR_VAR_ADD (gsf));
  gdouble total_var, chi2, lndetjac;

  _nc_galaxy_shape_factor_var_add_chi2 (self, pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2, &total_var, &chi2, &lndetjac);

  return -0.5 * chi2 + lndetjac - log (2.0 * M_PI * total_var);
}

static void
nc_galaxy_shape_factor_var_add_class_init (NcGalaxyShapeFactorVarAddClass *klass)
{
  NcGalaxyShapeFactorClass *gsf_class = NC_GALAXY_SHAPE_FACTOR_CLASS (klass);
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);

  object_class->constructed = &_nc_galaxy_shape_factor_var_add_constructed;
  object_class->finalize    = &_nc_galaxy_shape_factor_var_add_finalize;

  gsf_class->data_init        = &_nc_galaxy_shape_factor_var_add_data_init;
  gsf_class->prepare          = &_nc_galaxy_shape_factor_var_add_prepare;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_var_add_eval_marginal;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_var_add_eval_ln_marginal;
}

/**
 * nc_galaxy_shape_factor_var_add_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 *
 * Creates a new #NcGalaxyShapeFactorVarAdd.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorVarAdd.
 */
NcGalaxyShapeFactorVarAdd *
nc_galaxy_shape_factor_var_add_new (NcGalaxyWLObsEllipConv ellip_conv)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_VAR_ADD,
                       "ellip-conv", ellip_conv,
                       NULL);
}

/**
 * nc_galaxy_shape_factor_var_add_ref:
 * @gsfva: a #NcGalaxyShapeFactorVarAdd
 *
 * Increases the reference count of @gsfva by one.
 *
 * Returns: (transfer full): @gsfva.
 */
NcGalaxyShapeFactorVarAdd *
nc_galaxy_shape_factor_var_add_ref (NcGalaxyShapeFactorVarAdd *gsfva)
{
  return g_object_ref (gsfva);
}

/**
 * nc_galaxy_shape_factor_var_add_free:
 * @gsfva: a #NcGalaxyShapeFactorVarAdd
 *
 * Decreases the reference count of @gsfva by one.
 *
 */
void
nc_galaxy_shape_factor_var_add_free (NcGalaxyShapeFactorVarAdd *gsfva)
{
  g_object_unref (gsfva);
}

/**
 * nc_galaxy_shape_factor_var_add_clear:
 * @gsfva: a #NcGalaxyShapeFactorVarAdd
 *
 * Decreases the reference count of *@gsfva by one, and sets the pointer
 * *@gsfva to NULL.
 *
 */
void
nc_galaxy_shape_factor_var_add_clear (NcGalaxyShapeFactorVarAdd **gsfva)
{
  g_clear_object (gsfva);
}

