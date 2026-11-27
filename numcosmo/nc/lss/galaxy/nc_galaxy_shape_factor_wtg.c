/***************************************************************************
 *            nc_galaxy_shape_factor_wtg.c
 *
 *  Thu Nov 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_wtg.c
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * NcGalaxyShapeFactorWtG:
 *
 * Tangential shear-estimator likelihood of Applegate et al. (2014).
 *
 * This is the per-galaxy likelihood of the "P(z) method" of Applegate et al.
 * (2014) (Weighing the Giants III, their Eqs. 8-10). Instead of marginalizing
 * over the intrinsic ellipticity, as the other shape factors do, each galaxy's
 * observed ellipticity is turned into a linear estimate of its tangential
 * reduced shear,
 * $$\hat g_t = \frac{\mathrm{Re}\left(\epsilon_\mathrm{obs}\, e^{-2i\phi}\right)}{s\,R},$$
 * with $s = 2$ for the distortion convention
 * (#NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE) and $s = 1$ for the ellipticity
 * convention (#NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET), and $R$ the shear
 * responsivity, read from the catalog column "responsivity". Its density
 * about the calibrated model shear
 * $\mu_t = \mathrm{Re}\left\{\left[(1+m)\,g\,e^{2i\phi} + c\right] e^{-2i\phi}\right\}$
 * (Applegate et al. Eq. 10, the calibration applied to the model) is the
 * one-dimensional kernel of a #NcGalaxyShapePopShearScatter,
 * $$P(\hat g_t \mid \mu_t) = K(\hat g_t - \mu_t),$$
 * a density in $\hat g_t$, as in Applegate et al. The redshift marginalization
 * $\int \mathrm{d}z\, P(z)\, P(\hat g_t \mid \mu_t(z))$ is done, as for every
 * shape factor, by #NcDataClusterWLFactor. The cross component of the
 * ellipticity is not used.
 *
 * The defining approximation is the linear response
 * $\langle \epsilon_t \rangle = s R\, g_t$, with the scatter about it
 * described directly in shear space. The responsivity is a calibration of the
 * shape estimator and of the selected sample, not a model quantity: it is
 * never refitted with the population. For re-Gaussianization shapes (HSC) it
 * is $R = 1 - \langle e_\mathrm{rms}^2 \rangle$, the weighted mean over the
 * analyzed sample (Mandelbaum et al. 2018, Appendix A); for metacalibration
 * or AnaCal it is the ensemble mean response of the sample. A per-galaxy
 * column lets different subsamples carry their own value.
 *
 * The population in the #NcmMSet must be a #NcGalaxyShapePopShearScatter
 * (#NcGalaxyShapePopShearScatterGauss or #NcGalaxyShapePopShearScatterVoigt);
 * the measurement noise std_noise enters it in shear units,
 * $\sigma_\mathrm{noise}/(sR)$.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_wtg.h"
#include "nc/lss/galaxy/nc_galaxy_shape_pop_shear_scatter.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyShapeFactorWtG
{
  NcGalaxyShapeFactor parent_instance;

  /* s of the ellipticity convention, resolved once at construction
   * (ellip-conv is CONSTRUCT_ONLY). */
  gdouble s;
};

G_DEFINE_TYPE (NcGalaxyShapeFactorWtG, nc_galaxy_shape_factor_wtg, NC_TYPE_GALAXY_SHAPE_FACTOR)

typedef struct _NcGalaxyShapeFactorWtGData
{
  gdouble responsivity;
} NcGalaxyShapeFactorWtGData;

static void
nc_galaxy_shape_factor_wtg_init (NcGalaxyShapeFactorWtG *gsfwtg)
{
  gsfwtg->s = 0.0;
}

static void
_nc_galaxy_shape_factor_wtg_constructed (GObject *object)
{
  /* Chain up: start */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_wtg_parent_class)->constructed (object);
  {
    NcGalaxyShapeFactorWtG *gsfwtg = NC_GALAXY_SHAPE_FACTOR_WTG (object);

    switch (nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (gsfwtg)))
    {
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE:
        /* distortion chi: <chi> = 2 R g */
        gsfwtg->s = 2.0;
        break;
      case NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET:
        /* ellipticity epsilon: <epsilon> = R g */
        gsfwtg->s = 1.0;
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
    }
  }
}

static void
_nc_galaxy_shape_factor_wtg_finalize (GObject *object)
{
  /* Chain up: end */
  G_OBJECT_CLASS (nc_galaxy_shape_factor_wtg_parent_class)->finalize (object);
}

static void
_nc_galaxy_shape_factor_wtg_ldata_read_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyShapeFactorWtGData *ldata = (NcGalaxyShapeFactorWtGData *) data->ldata;
  const gdouble responsivity        = ncm_catalog_get (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_WTG_COL_RESPONSIVITY, i, NULL);

  if (!(responsivity > 0.0))
    g_error ("NcGalaxyShapeFactorWtG: row %u has a non-positive responsivity (%g).", i, responsivity);

  ldata->responsivity = responsivity;
}

static void
_nc_galaxy_shape_factor_wtg_ldata_write_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyShapeFactorWtGData *ldata = (NcGalaxyShapeFactorWtGData *) data->ldata;

  ncm_catalog_set (NCM_CATALOG (obs), NC_GALAXY_SHAPE_FACTOR_WTG_COL_RESPONSIVITY, i, ldata->responsivity, NULL);
}

static void
_nc_galaxy_shape_factor_wtg_ldata_required_columns (NcGalaxyShapeFactorData *data, GList **columns)
{
  *columns = g_list_append (*columns, g_strdup (NC_GALAXY_SHAPE_FACTOR_WTG_COL_RESPONSIVITY));
}

static void
_nc_galaxy_shape_factor_wtg_data_init (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorWtGData *ldata = g_new0 (NcGalaxyShapeFactorWtGData, 1);

  /* The responsivity of an unbiased shear estimator, until a row is read
   * or it is set. */
  ldata->responsivity = 1.0;

  data->ldata                  = ldata;
  data->ldata_destroy          = &g_free;
  data->ldata_read_row         = &_nc_galaxy_shape_factor_wtg_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_shape_factor_wtg_ldata_write_row;
  data->ldata_required_columns = &_nc_galaxy_shape_factor_wtg_ldata_required_columns;
}

static gboolean
_nc_galaxy_shape_factor_wtg_accepts_pop (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop)
{
  return NC_IS_GALAXY_SHAPE_POP_SHEAR_SCATTER (pop);
}

/* 1 / (s R): the factor taking an ellipticity to the shear estimate */
static inline gdouble
_nc_galaxy_shape_factor_wtg_eps_to_g (NcGalaxyShapeFactorWtG *gsfwtg, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorWtGData *ldata = (NcGalaxyShapeFactorWtGData *) data->ldata;

  return 1.0 / (gsfwtg->s * ldata->responsivity);
}

static gdouble
_nc_galaxy_shape_factor_wtg_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  NcGalaxyShapeFactorWtG *gsfwtg = NC_GALAXY_SHAPE_FACTOR_WTG (gsf);
  const gdouble eps_to_g         = _nc_galaxy_shape_factor_wtg_eps_to_g (gsfwtg, data);
  const gdouble noise_g          = data->std_noise * eps_to_g;
  gdouble rot_re, rot_im;

  nc_galaxy_shape_factor_data_get_rotation (data, &rot_re, &rot_im);

  {
    /* Tangential components, Re (v e^{-2 i phi}) */
    const gdouble hat_g_t = (epsilon_obs_1 * rot_re + epsilon_obs_2 * rot_im) * eps_to_g;
    const gdouble mu_t    = g_1 * rot_re + g_2 * rot_im;

    return nc_galaxy_shape_pop_shear_scatter_eval_ln_p (NC_GALAXY_SHAPE_POP_SHEAR_SCATTER (pop), data->pop_data, hat_g_t - mu_t, noise_g * noise_g, eps_to_g);
  }
}

static gdouble
_nc_galaxy_shape_factor_wtg_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2)
{
  return exp (_nc_galaxy_shape_factor_wtg_eval_ln_marginal (gsf, pop, data, g_1, g_2, epsilon_obs_1, epsilon_obs_2));
}

/* The shear estimate is mu + delta in the tangential and cross components,
 * each delta drawn independently from the scatter kernel, and the ellipticity
 * is s R times it, rotated to data->coord. There is no intrinsic ellipticity:
 * epsilon_int is set to zero. */
static void
_nc_galaxy_shape_factor_wtg_gen_obs (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, NcmRNG *rng, const gboolean behind,
                                     const gdouble g_1, const gdouble g_2, const gdouble rot_re, const gdouble rot_im)
{
  NcGalaxyShapeFactorWtG *gsfwtg      = NC_GALAXY_SHAPE_FACTOR_WTG (gsf);
  NcGalaxyShapePopShearScatter *gspss = NC_GALAXY_SHAPE_POP_SHEAR_SCATTER (pop);
  const gdouble eps_to_g              = _nc_galaxy_shape_factor_wtg_eps_to_g (gsfwtg, data);
  const gdouble noise_g               = data->std_noise * eps_to_g;
  const gdouble var_noise_g           = noise_g * noise_g;
  const gdouble mu_t                  = g_1 * rot_re + g_2 * rot_im;
  const gdouble mu_x                  = g_2 * rot_re - g_1 * rot_im;
  gdouble hat_g_t, hat_g_x;

  nc_galaxy_shape_pop_prepare (pop, data->pop_data);

  hat_g_t = mu_t + nc_galaxy_shape_pop_shear_scatter_gen_delta (gspss, data->pop_data, rng, var_noise_g, eps_to_g);
  hat_g_x = mu_x + nc_galaxy_shape_pop_shear_scatter_gen_delta (gspss, data->pop_data, rng, var_noise_g, eps_to_g);

  /* (hat_g_t + i hat_g_x) e^{2 i phi} / eps_to_g */
  data->epsilon_obs_1 = (hat_g_t * rot_re - hat_g_x * rot_im) / eps_to_g;
  data->epsilon_obs_2 = (hat_g_t * rot_im + hat_g_x * rot_re) / eps_to_g;
  data->epsilon_int_1 = 0.0;
  data->epsilon_int_2 = 0.0;
}

static void
nc_galaxy_shape_factor_wtg_class_init (NcGalaxyShapeFactorWtGClass *klass)
{
  NcGalaxyShapeFactorClass *gsf_class = NC_GALAXY_SHAPE_FACTOR_CLASS (klass);
  GObjectClass *object_class          = G_OBJECT_CLASS (klass);

  object_class->constructed = &_nc_galaxy_shape_factor_wtg_constructed;
  object_class->finalize    = &_nc_galaxy_shape_factor_wtg_finalize;

  gsf_class->data_init        = &_nc_galaxy_shape_factor_wtg_data_init;
  gsf_class->accepts_pop      = &_nc_galaxy_shape_factor_wtg_accepts_pop;
  gsf_class->eval_marginal    = &_nc_galaxy_shape_factor_wtg_eval_marginal;
  gsf_class->eval_ln_marginal = &_nc_galaxy_shape_factor_wtg_eval_ln_marginal;
  gsf_class->gen_obs          = &_nc_galaxy_shape_factor_wtg_gen_obs;
}

/**
 * nc_galaxy_shape_factor_wtg_new:
 * @ellip_conv: a #NcGalaxyWLObsEllipConv
 *
 * Creates a new #NcGalaxyShapeFactorWtG.
 *
 * Returns: (transfer full): a new #NcGalaxyShapeFactorWtG.
 */
NcGalaxyShapeFactorWtG *
nc_galaxy_shape_factor_wtg_new (NcGalaxyWLObsEllipConv ellip_conv)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_WTG,
                       "ellip-conv", ellip_conv,
                       NULL);
}

/**
 * nc_galaxy_shape_factor_wtg_ref:
 * @gsfwtg: a #NcGalaxyShapeFactorWtG
 *
 * Increases the reference count of @gsfwtg by one.
 *
 * Returns: (transfer full): @gsfwtg.
 */
NcGalaxyShapeFactorWtG *
nc_galaxy_shape_factor_wtg_ref (NcGalaxyShapeFactorWtG *gsfwtg)
{
  return g_object_ref (gsfwtg);
}

/**
 * nc_galaxy_shape_factor_wtg_free:
 * @gsfwtg: a #NcGalaxyShapeFactorWtG
 *
 * Decreases the reference count of @gsfwtg by one.
 *
 */
void
nc_galaxy_shape_factor_wtg_free (NcGalaxyShapeFactorWtG *gsfwtg)
{
  g_object_unref (gsfwtg);
}

/**
 * nc_galaxy_shape_factor_wtg_clear:
 * @gsfwtg: a #NcGalaxyShapeFactorWtG
 *
 * Decreases the reference count of *@gsfwtg by one, and sets the pointer
 * *@gsfwtg to NULL.
 *
 */
void
nc_galaxy_shape_factor_wtg_clear (NcGalaxyShapeFactorWtG **gsfwtg)
{
  g_clear_object (gsfwtg);
}

/**
 * nc_galaxy_shape_factor_wtg_data_set_responsivity:
 * @gsfwtg: a #NcGalaxyShapeFactorWtG
 * @data: a #NcGalaxyShapeFactorData created by @gsfwtg
 * @responsivity: the shear responsivity $R > 0$
 *
 * Sets the shear responsivity of the galaxy in @data, the value written to
 * and read from the "responsivity" column.
 *
 */
void
nc_galaxy_shape_factor_wtg_data_set_responsivity (NcGalaxyShapeFactorWtG *gsfwtg, NcGalaxyShapeFactorData *data, const gdouble responsivity)
{
  NcGalaxyShapeFactorWtGData *ldata = (NcGalaxyShapeFactorWtGData *) data->ldata;

  g_assert_cmpfloat (responsivity, >, 0.0);
  ldata->responsivity = responsivity;
}

/**
 * nc_galaxy_shape_factor_wtg_data_get_responsivity:
 * @gsfwtg: a #NcGalaxyShapeFactorWtG
 * @data: a #NcGalaxyShapeFactorData created by @gsfwtg
 *
 * Returns: the shear responsivity of the galaxy in @data.
 */
gdouble
nc_galaxy_shape_factor_wtg_data_get_responsivity (NcGalaxyShapeFactorWtG *gsfwtg, NcGalaxyShapeFactorData *data)
{
  NcGalaxyShapeFactorWtGData *ldata = (NcGalaxyShapeFactorWtGData *) data->ldata;

  return ldata->responsivity;
}

/**
 * nc_galaxy_shape_factor_wtg_data_get_shear_estimate:
 * @gsfwtg: a #NcGalaxyShapeFactorWtG
 * @data: a #NcGalaxyShapeFactorData created by @gsfwtg
 *
 * Computes the tangential shear estimate $\hat g_t$ of the galaxy in @data,
 * using the rotation cached by the last radius update (see
 * nc_galaxy_shape_factor_update_data_radius()).
 *
 * Returns: $\hat g_t$.
 */
gdouble
nc_galaxy_shape_factor_wtg_data_get_shear_estimate (NcGalaxyShapeFactorWtG *gsfwtg, NcGalaxyShapeFactorData *data)
{
  const gdouble eps_to_g = _nc_galaxy_shape_factor_wtg_eps_to_g (gsfwtg, data);
  gdouble rot_re, rot_im;

  nc_galaxy_shape_factor_data_get_rotation (data, &rot_re, &rot_im);

  return (data->epsilon_obs_1 * rot_re + data->epsilon_obs_2 * rot_im) * eps_to_g;
}

