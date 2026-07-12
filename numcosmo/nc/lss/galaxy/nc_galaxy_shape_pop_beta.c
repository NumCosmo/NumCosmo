/***************************************************************************
 *            nc_galaxy_shape_pop_beta.c
 *
 *  Thu Jun 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop_beta.c
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
 * NcGalaxyShapePopBeta:
 *
 * Beta intrinsic ellipticity distribution (BetaGlobal).
 *
 * The squared intrinsic ellipticity modulus $x = |\chi_I|^2$ follows a Beta
 * distribution $x \sim \mathrm{Beta}(\alpha,\beta)$ with $\alpha = \mu\nu$ and
 * $\beta = (1-\mu)\nu$, where $\mu$ controls the typical ellipticity and $\nu$
 * the concentration:
 * $$P(x) = \frac{x^{\alpha-1}(1-x)^{\beta-1}}{B(\alpha,\beta)}, \qquad x\in[0,1].$$
 * eval_p_rho2() overrides the generic default with the equivalent rational
 * form in $\rho^2$, $x=\rho^2/(1+\rho^2)$:
 * $$P(x(\rho^2)) = \frac{\rho^{2(\alpha-1)}}{B(\alpha,\beta)\,(1+\rho^2)^{\alpha+\beta-2}},$$
 * avoiding forming $1-x$ by subtraction.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_pop_beta.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <complex.h>
#include <gsl/gsl_sf_gamma.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyShapePopBeta
{
  NcGalaxyShapePop parent_instance;
};

typedef struct _NcGalaxyShapePopBetaLData
{
  gdouble lnnorm; /* -ln B(alpha, beta) */
  gdouble alpha;  /* mu nu */
  gdouble beta;   /* (1 - mu) nu */
} NcGalaxyShapePopBetaLData;

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxyShapePopBeta, nc_galaxy_shape_pop_beta, NC_TYPE_GALAXY_SHAPE_POP);

static void
nc_galaxy_shape_pop_beta_init (NcGalaxyShapePopBeta *gspb)
{
}

static void
_nc_galaxy_shape_pop_beta_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_shape_pop_beta_parent_class)->finalize (object);
}

static void _nc_galaxy_shape_pop_beta_data_init (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
static void _nc_galaxy_shape_pop_beta_prepare (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
static gdouble _nc_galaxy_shape_pop_beta_ldata_get_mode_x (NcGalaxyShapePopData *data);
static gdouble _nc_galaxy_shape_pop_beta_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble x);
static gdouble _nc_galaxy_shape_pop_beta_eval_p_rho2 (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble rho2);
static void _nc_galaxy_shape_pop_beta_gen (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2);
static gdouble _nc_galaxy_shape_pop_beta_e_rms (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);

static void
nc_galaxy_shape_pop_beta_class_init (NcGalaxyShapePopBetaClass *klass)
{
  NcGalaxyShapePopClass *gsp_class = NC_GALAXY_SHAPE_POP_CLASS (klass);
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_shape_pop_beta_finalize;

  ncm_model_class_set_name_nick (model_class, "Beta intrinsic ellipticity distribution", "BetaIntrinsic");
  ncm_model_class_add_params (model_class, NC_GALAXY_SHAPE_POP_BETA_SPARAM_LEN, 0, PROP_LEN);

  /**
   * NcGalaxyShapePopBeta:mu:
   *
   * The mean $\mu = \langle x \rangle$ of the Beta distribution of $x = |\chi_I|^2$.
   *
   */
  ncm_model_class_set_sparam (model_class,
                              NC_GALAXY_SHAPE_POP_BETA_MU,
                              "\\mu",
                              "mu", 1.0e-3, 1.0 - 1.0e-3, 1.0e-2,
                              NC_GALAXY_SHAPE_POP_BETA_DEFAULT_PARAMS_ABSTOL,
                              NC_GALAXY_SHAPE_POP_BETA_DEFAULT_MU,
                              NCM_PARAM_TYPE_FIXED);

  /**
   * NcGalaxyShapePopBeta:nu:
   *
   * The concentration $\nu = \alpha + \beta$ of the Beta distribution.
   *
   */
  ncm_model_class_set_sparam (model_class,
                              NC_GALAXY_SHAPE_POP_BETA_NU,
                              "\\nu",
                              "nu", 1.0e-2, 1.0e3, 1.0e-1,
                              NC_GALAXY_SHAPE_POP_BETA_DEFAULT_PARAMS_ABSTOL,
                              NC_GALAXY_SHAPE_POP_BETA_DEFAULT_NU,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  gsp_class->data_init   = &_nc_galaxy_shape_pop_beta_data_init;
  gsp_class->prepare     = &_nc_galaxy_shape_pop_beta_prepare;
  gsp_class->eval_p      = &_nc_galaxy_shape_pop_beta_eval_p;
  gsp_class->eval_p_rho2 = &_nc_galaxy_shape_pop_beta_eval_p_rho2;
  gsp_class->gen         = &_nc_galaxy_shape_pop_beta_gen;
  gsp_class->e_rms       = &_nc_galaxy_shape_pop_beta_e_rms;
}

#define VECTOR (NCM_MODEL (gsp))
#define MU     (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SHAPE_POP_BETA_MU))
#define NU     (ncm_model_orig_param_get (VECTOR, NC_GALAXY_SHAPE_POP_BETA_NU))

static void
_nc_galaxy_shape_pop_beta_ldata_noop (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i)
{
}

static void
_nc_galaxy_shape_pop_beta_ldata_required_columns (NcGalaxyShapePopData *data, GList **columns)
{
}

static void
_nc_galaxy_shape_pop_beta_data_init (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  NcGalaxyShapePopBetaLData *ldata = g_new0 (NcGalaxyShapePopBetaLData, 1);

  data->ldata                  = ldata;
  data->ldata_destroy          = &g_free;
  data->ldata_read_row         = &_nc_galaxy_shape_pop_beta_ldata_noop;
  data->ldata_write_row        = &_nc_galaxy_shape_pop_beta_ldata_noop;
  data->ldata_required_columns = &_nc_galaxy_shape_pop_beta_ldata_required_columns;
  data->ldata_get_mode_x       = &_nc_galaxy_shape_pop_beta_ldata_get_mode_x;
}

static gdouble
_nc_galaxy_shape_pop_beta_ldata_get_mode_x (NcGalaxyShapePopData *data)
{
  NcGalaxyShapePopBetaLData *ldata = (NcGalaxyShapePopBetaLData *) data->ldata;

  /* Mode of a standard Beta(alpha,beta) density: interior stationary point
   * when both shape parameters exceed 1; otherwise the density is monotone
   * (or diverges at both ends), and the mode sits at a disc boundary or is
   * undefined -- 0 is the closest meaningful hint in those cases (it is
   * always in-range and never worse than assuming the Gauss-like default). */
  if ((ldata->alpha > 1.0) && (ldata->beta > 1.0))
    return (ldata->alpha - 1.0) / (ldata->alpha + ldata->beta - 2.0);
  else
    return 0.0;
}

static void
_nc_galaxy_shape_pop_beta_prepare (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  NcGalaxyShapePopBetaLData *ldata = (NcGalaxyShapePopBetaLData *) data->ldata;
  const gdouble mu                 = MU;
  const gdouble nu                 = NU;
  const gdouble alpha              = mu * nu;
  const gdouble beta               = (1.0 - mu) * nu;

  ldata->lnnorm = -gsl_sf_lnbeta (alpha, beta);
  ldata->alpha  = alpha;
  ldata->beta   = beta;
}

static gdouble
_nc_galaxy_shape_pop_beta_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble x)
{
  NcGalaxyShapePopBetaLData *ldata = (NcGalaxyShapePopBetaLData *) data->ldata;

  /* Evaluated in log-space and exponentiated once: for concentrated
   * populations (large alpha, beta) forming pow(x,alpha-1) and the
   * normalization 1/B(alpha,beta) as separate factors overflows/underflows
   * well within the model's own declared nu range (e.g. mu=0.5, nu=1000
   * already pushes -lnbeta past 690, right at exp()'s overflow edge),
   * silently producing NaN (0*inf) with no diagnostic. */
  return exp ((ldata->alpha - 1.0) * log (x) + (ldata->beta - 1.0) * log1p (-x) + ldata->lnnorm);
}

static gdouble
_nc_galaxy_shape_pop_beta_eval_p_rho2 (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble rho2)
{
  NcGalaxyShapePopBetaLData *ldata = (NcGalaxyShapePopBetaLData *) data->ldata;

  /* (1-x)^{beta-1} x^{alpha-1} with 1-x=1/(1+rho2), x=rho2/(1+rho2): avoids
   * ever forming x or 1-x by subtraction (see the class doc), evaluated in
   * log-space for the same overflow reason as eval_p() above. */
  return exp ((ldata->alpha - 1.0) * log (rho2) - (ldata->alpha + ldata->beta - 2.0) * log1p (rho2) + ldata->lnnorm);
}

static void
_nc_galaxy_shape_pop_beta_gen (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2)
{
  NcGalaxyShapePopBetaLData *ldata = (NcGalaxyShapePopBetaLData *) data->ldata;
  const gdouble x                  = ncm_rng_beta_gen (rng, ldata->alpha, ldata->beta);
  const gdouble r                  = sqrt (x);
  const gdouble theta              = ncm_rng_uniform_gen (rng, 0.0, 2.0 * M_PI);

  *e_int_1 = r * cos (theta);
  *e_int_2 = r * sin (theta);
}

static gdouble
_nc_galaxy_shape_pop_beta_e_rms (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  /* <|chi|^2> = <x> = mu, so e_rms = sqrt(mu / 2). */
  return sqrt (0.5 * MU);
}

/**
 * nc_galaxy_shape_pop_beta_new:
 *
 * Creates a new #NcGalaxyShapePopBeta.
 *
 * Returns: (transfer full): a new #NcGalaxyShapePopBeta.
 */
NcGalaxyShapePopBeta *
nc_galaxy_shape_pop_beta_new (void)
{
  NcGalaxyShapePopBeta *gspb = g_object_new (NC_TYPE_GALAXY_SHAPE_POP_BETA,
                                             NULL);

  return gspb;
}

/**
 * nc_galaxy_shape_pop_beta_ref:
 * @gspb: a #NcGalaxyShapePopBeta
 *
 * Increases the reference count of @gspb by one.
 *
 * Returns: (transfer full): @gspb.
 */
NcGalaxyShapePopBeta *
nc_galaxy_shape_pop_beta_ref (NcGalaxyShapePopBeta *gspb)
{
  return g_object_ref (gspb);
}

/**
 * nc_galaxy_shape_pop_beta_free:
 * @gspb: a #NcGalaxyShapePopBeta
 *
 * Decreases the reference count of @gspb by one.
 *
 */
void
nc_galaxy_shape_pop_beta_free (NcGalaxyShapePopBeta *gspb)
{
  g_object_unref (gspb);
}

/**
 * nc_galaxy_shape_pop_beta_clear:
 * @gspb: a #NcGalaxyShapePopBeta
 *
 * Decreases the reference count of *@gspb by one, and sets the pointer *@gspb to
 * NULL.
 *
 */
void
nc_galaxy_shape_pop_beta_clear (NcGalaxyShapePopBeta **gspb)
{
  g_clear_object (gspb);
}

