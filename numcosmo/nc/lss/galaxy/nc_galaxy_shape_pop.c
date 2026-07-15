/***************************************************************************
 *            nc_galaxy_shape_pop.c
 *
 *  Thu Jun 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop.c
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
 * NcGalaxyShapePop:
 *
 * Abstract model for the intrinsic galaxy ellipticity distribution.
 *
 * This small #NcmModel describes the probability density of the intrinsic
 * ellipticity modulus squared $x = |\chi_I|^2 \in [0,1]$. It is a pure
 * parametric law for which the shape calculator (#NcGalaxyWLShapeCalc) performs
 * the exact marginalization over the intrinsic ellipticity.
 *
 * Each concrete model implements its density directly through the
 * nc_galaxy_shape_pop_eval_p() virtual: there is no shared factorization or
 * quadrature scheme imposed at this level, each subclass is simply
 * responsible for its own normalized $P(x)$. nc_galaxy_shape_pop_eval_p_rho2()
 * is the same density parameterized by $\rho^2$ instead of $x$ (related by
 * $x=\rho^2/(1+\rho^2)$, useful for calculators working in a
 * disc-compactified plane, e.g. #NcGalaxyShapeFactorQuad); its default
 * implementation just substitutes into eval_p(), but a subclass may override
 * it with a more direct/better-conditioned form when its density happens to
 * have one (see #NcGalaxyShapePopBeta). Following the #NcGalaxySDShape /
 * #NcGalaxySDPosition idiom, the resolved per-galaxy state lives in a
 * #NcGalaxyShapePopData: typed public fields (@e_rms) plus an opaque @ldata
 * holding the subclass-specific, varying/updatable resolved parameters.
 * nc_galaxy_shape_pop_prepare() resolves the model parameters (and any
 * per-galaxy @e_rms) into that data.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_pop.h"

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_ABSTRACT_TYPE (NcGalaxyShapePop, nc_galaxy_shape_pop, NCM_TYPE_MODEL);
G_DEFINE_BOXED_TYPE (NcGalaxyShapePopData, nc_galaxy_shape_pop_data, nc_galaxy_shape_pop_data_ref, nc_galaxy_shape_pop_data_unref); /* LCOV_EXCL_LINE */

static void
nc_galaxy_shape_pop_init (NcGalaxyShapePop *gsp)
{
}

static void
_nc_galaxy_shape_pop_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_shape_pop_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_shape_pop, NC_TYPE_GALAXY_SHAPE_POP);

/* LCOV_EXCL_START */
static void
_nc_galaxy_shape_pop_data_init (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  g_error ("_nc_galaxy_shape_pop_data_init: method not implemented.");
}

static void
_nc_galaxy_shape_pop_prepare (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  g_error ("_nc_galaxy_shape_pop_prepare: method not implemented.");
}

static gdouble
_nc_galaxy_shape_pop_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble x)
{
  g_error ("_nc_galaxy_shape_pop_eval_p: method not implemented.");

  return 0.0;
}

static void
_nc_galaxy_shape_pop_gen (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2)
{
  g_error ("_nc_galaxy_shape_pop_gen: method not implemented.");
}

static gdouble
_nc_galaxy_shape_pop_e_rms (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  g_error ("_nc_galaxy_shape_pop_e_rms: method not implemented.");

  return 0.0;
}

static void
_nc_galaxy_shape_pop_eval_p_rho2_g_series (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data,
                                           const NcmLaurentSeriesTPS *rho2_series, NcmLaurentSeriesTPS *out)
{
  g_error ("_nc_galaxy_shape_pop_eval_p_rho2_g_series: method not implemented.");
}

/* LCOV_EXCL_STOP */

/* Real (non-stub) default: substitutes x = rho2/(1+rho2) into eval_p(). Any
 * subclass whose density has a more direct/better-conditioned form in rho2
 * (e.g. #NcGalaxyShapePopBeta) may override this instead of inheriting it. */
static gdouble
_nc_galaxy_shape_pop_eval_p_rho2 (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble rho2)
{
  const gdouble x = rho2 / (1.0 + rho2);

  return NC_GALAXY_SHAPE_POP_GET_CLASS (gsp)->eval_p (gsp, data, x);
}

static void
nc_galaxy_shape_pop_class_init (NcGalaxyShapePopClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_shape_pop_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy intrinsic ellipticity distribution", "GalaxyShapePop");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);
  ncm_mset_model_register_id (model_class, "NcGalaxyShapePop", "Galaxy intrinsic ellipticity distribution", NULL, FALSE, NCM_MSET_MODEL_MAIN);
  ncm_model_class_check_params_info (model_class);

  klass->data_init            = &_nc_galaxy_shape_pop_data_init;
  klass->prepare              = &_nc_galaxy_shape_pop_prepare;
  klass->eval_p               = &_nc_galaxy_shape_pop_eval_p;
  klass->eval_p_rho2          = &_nc_galaxy_shape_pop_eval_p_rho2;
  klass->gen                  = &_nc_galaxy_shape_pop_gen;
  klass->e_rms                = &_nc_galaxy_shape_pop_e_rms;
  klass->eval_p_rho2_g_series = &_nc_galaxy_shape_pop_eval_p_rho2_g_series;
}

/**
 * nc_galaxy_shape_pop_ref:
 * @gsp: a #NcGalaxyShapePop
 *
 * Increases the reference count of @gsp by one.
 *
 * Returns: (transfer full): @gsp.
 */
NcGalaxyShapePop *
nc_galaxy_shape_pop_ref (NcGalaxyShapePop *gsp)
{
  return g_object_ref (gsp);
}

/**
 * nc_galaxy_shape_pop_free:
 * @gsp: a #NcGalaxyShapePop
 *
 * Decreases the reference count of @gsp by one.
 *
 */
void
nc_galaxy_shape_pop_free (NcGalaxyShapePop *gsp)
{
  g_object_unref (gsp);
}

/**
 * nc_galaxy_shape_pop_clear:
 * @gsp: a #NcGalaxyShapePop
 *
 * Decreases the reference count of *@gsp by one, and sets the pointer *@gsp to
 * NULL.
 *
 */
void
nc_galaxy_shape_pop_clear (NcGalaxyShapePop **gsp)
{
  g_clear_object (gsp);
}

/**
 * nc_galaxy_shape_pop_data_new:
 * @gsp: a #NcGalaxyShapePop
 *
 * Creates a new per-galaxy #NcGalaxyShapePopData for @gsp.
 *
 * Returns: (transfer full): a new #NcGalaxyShapePopData.
 */
NcGalaxyShapePopData *
nc_galaxy_shape_pop_data_new (NcGalaxyShapePop *gsp)
{
  NcGalaxyShapePopData *data = g_new0 (NcGalaxyShapePopData, 1);

  data->e_rms                  = 0.0;
  data->ldata                  = NULL;
  data->ldata_destroy          = NULL;
  data->ldata_read_row         = NULL;
  data->ldata_write_row        = NULL;
  data->ldata_required_columns = NULL;
  data->ldata_get_sigma        = NULL;
  data->ldata_get_mode_x       = NULL;

  g_atomic_ref_count_init (&data->ref_count);
  NC_GALAXY_SHAPE_POP_GET_CLASS (gsp)->data_init (gsp, data);

  g_assert_nonnull (data->ldata_destroy);
  g_assert_nonnull (data->ldata_read_row);
  g_assert_nonnull (data->ldata_write_row);
  g_assert_nonnull (data->ldata_required_columns);

  return data;
}

/**
 * nc_galaxy_shape_pop_data_ref:
 * @data: a #NcGalaxyShapePopData
 *
 * Increases the reference count of @data by one.
 *
 * Returns: (transfer full): @data.
 */
NcGalaxyShapePopData *
nc_galaxy_shape_pop_data_ref (NcGalaxyShapePopData *data)
{
  g_atomic_ref_count_inc (&data->ref_count);

  return data;
}

/**
 * nc_galaxy_shape_pop_data_unref:
 * @data: a #NcGalaxyShapePopData
 *
 * Decreases the reference count of @data by one. If it reaches 0, @data is freed.
 *
 */
void
nc_galaxy_shape_pop_data_unref (NcGalaxyShapePopData *data)
{
  if (g_atomic_ref_count_dec (&data->ref_count))
  {
    g_assert_nonnull (data->ldata_destroy);
    data->ldata_destroy (data->ldata);
    g_free (data);
  }
}

/**
 * nc_galaxy_shape_pop_data_read_row:
 * @data: a #NcGalaxyShapePopData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Reads the intrinsic-distribution inputs of row @i.
 *
 */
void
nc_galaxy_shape_pop_data_read_row (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i)
{
  data->ldata_read_row (data, obs, i);
}

/**
 * nc_galaxy_shape_pop_data_write_row:
 * @data: a #NcGalaxyShapePopData
 * @obs: a #NcGalaxyWLObs
 * @i: the row index
 *
 * Writes the intrinsic-distribution inputs of row @i.
 *
 */
void
nc_galaxy_shape_pop_data_write_row (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i)
{
  data->ldata_write_row (data, obs, i);
}

/**
 * nc_galaxy_shape_pop_data_required_columns:
 * @data: a #NcGalaxyShapePopData
 *
 * Returns: (element-type utf8) (transfer full): the required columns.
 */
GList *
nc_galaxy_shape_pop_data_required_columns (NcGalaxyShapePopData *data)
{
  GList *columns = NULL;

  data->ldata_required_columns (data, &columns);

  return columns;
}

/**
 * nc_galaxy_shape_pop_prepare: (virtual prepare)
 * @gsp: a #NcGalaxyShapePop
 * @data: a #NcGalaxyShapePopData
 *
 * Resolves the model parameters (and the per-galaxy @data->e_rms, for per-galaxy
 * models) into @data: the subclass-specific parameters stored in @data->ldata
 * that nc_galaxy_shape_pop_eval_p() / nc_galaxy_shape_pop_eval_p_rho2() read.
 *
 */
void
nc_galaxy_shape_pop_prepare (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  NC_GALAXY_SHAPE_POP_GET_CLASS (gsp)->prepare (gsp, data);
}

/**
 * nc_galaxy_shape_pop_eval_p: (virtual eval_p)
 * @gsp: a #NcGalaxyShapePop
 * @data: a resolved #NcGalaxyShapePopData
 * @x: the squared ellipticity modulus $x = |\chi_I|^2 \in [0,1]$
 *
 * Evaluates the normalized density $P(x)$ of the squared intrinsic
 * ellipticity modulus, reading the resolved parameters from @data (no live
 * parameter access).
 *
 * Returns: the density $P(x)$.
 */
gdouble
nc_galaxy_shape_pop_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble x)
{
  return NC_GALAXY_SHAPE_POP_GET_CLASS (gsp)->eval_p (gsp, data, x);
}

/**
 * nc_galaxy_shape_pop_eval_p_rho2: (virtual eval_p_rho2)
 * @gsp: a #NcGalaxyShapePop
 * @data: a resolved #NcGalaxyShapePopData
 * @rho2: $\rho^2 = u^2+v^2 \geq 0$, related to $x=|\chi_I|^2$ by $x = \rho^2/(1+\rho^2)$
 *
 * Evaluates the same density as nc_galaxy_shape_pop_eval_p(), but parameterized
 * directly by $\rho^2$ instead of $x$. The default implementation just
 * substitutes $x=\rho^2/(1+\rho^2)$ into eval_p(); a subclass overrides it
 * only when its density has a genuinely more direct or better-conditioned
 * form in $\rho^2$ (see #NcGalaxyShapePopBeta, whose density is a rational
 * power of $\rho^2$ with no need to ever form $1-x$ by subtraction).
 * Intended for calculators that already work in $(u,v)$ via a
 * disc-compactifying substitution (e.g. #NcGalaxyShapeFactorQuad), where
 * $\rho^2$ is the natural coordinate.
 *
 * Returns: the density $P(x)$.
 */
gdouble
nc_galaxy_shape_pop_eval_p_rho2 (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble rho2)
{
  return NC_GALAXY_SHAPE_POP_GET_CLASS (gsp)->eval_p_rho2 (gsp, data, rho2);
}

/**
 * nc_galaxy_shape_pop_eval_p_rho2_g_series:
 * @gsp: a #NcGalaxyShapePop
 * @data: a resolved #NcGalaxyShapePopData
 * @rho2_series: $\rho^2(g)=|\chi_I(\chi_L,g)|^2$'s own $g$-Taylor
 * coefficients, population-independent
 * @out: this population's $P(\rho^2(g))$ composition, as $g$-Taylor
 * coefficients (same order as @rho2_series)
 *
 * Taylor-in-$g$ analog of nc_galaxy_shape_pop_eval_p_rho2(): composes this
 * population's own normalized density with the (already computed,
 * population-independent) shear-map series $\rho^2(g)$, order by order in
 * $g$. There is no generic default -- every subclass used with
 * #NcGalaxyShapeFactorSeriesLensed must provide its own implementation; the
 * base class errors clearly otherwise.
 */
void
nc_galaxy_shape_pop_eval_p_rho2_g_series (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data,
                                          const NcmLaurentSeriesTPS *rho2_series, NcmLaurentSeriesTPS *out)
{
  NC_GALAXY_SHAPE_POP_GET_CLASS (gsp)->eval_p_rho2_g_series (gsp, data, rho2_series, out);
}

/**
 * nc_galaxy_shape_pop_gen: (virtual gen)
 * @gsp: a #NcGalaxyShapePop
 * @data: a resolved #NcGalaxyShapePopData
 * @rng: a #NcmRNG
 * @e_int_1: (out): first component of the sampled intrinsic ellipticity
 * @e_int_2: (out): second component of the sampled intrinsic ellipticity
 *
 * Samples an intrinsic ellipticity $\chi_I = e_1 + i e_2$ from @data.
 *
 */
void
nc_galaxy_shape_pop_gen (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, NcmRNG *rng, gdouble *e_int_1, gdouble *e_int_2)
{
  NC_GALAXY_SHAPE_POP_GET_CLASS (gsp)->gen (gsp, data, rng, e_int_1, e_int_2);
}

/**
 * nc_galaxy_shape_pop_e_rms: (virtual e_rms)
 * @gsp: a #NcGalaxyShapePop
 * @data: a resolved #NcGalaxyShapePopData
 *
 * Computes the per-component intrinsic RMS ellipticity
 * $e_\mathrm{rms} = \sqrt{\langle |\chi_I|^2 \rangle / 2}$ implied by @data.
 *
 * Returns: the intrinsic RMS ellipticity.
 */
gdouble
nc_galaxy_shape_pop_e_rms (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  return NC_GALAXY_SHAPE_POP_GET_CLASS (gsp)->e_rms (gsp, data);
}

/**
 * nc_galaxy_shape_pop_get_sigma:
 * @gsp: a #NcGalaxyShapePop
 * @data: a resolved #NcGalaxyShapePopData
 *
 * Gets the (untruncated) Gaussian width sigma resolved for this galaxy, for
 * models that parameterize their density this way (the Gauss family, Global
 * or per-galaxy). This is a per-@data capability (@data->ldata_get_sigma),
 * not a virtual method: any concrete model sharing the same internal
 * representation can populate it, regardless of its position in the class
 * hierarchy.
 *
 * Returns: the resolved Gaussian width sigma.
 */
gdouble
nc_galaxy_shape_pop_get_sigma (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  if (data->ldata_get_sigma == NULL)
    g_error ("nc_galaxy_shape_pop_get_sigma: %s does not support a Gaussian width sigma.",
             G_OBJECT_TYPE_NAME (gsp));

  return data->ldata_get_sigma (data);
}

/**
 * nc_galaxy_shape_pop_get_mode_x:
 * @gsp: a #NcGalaxyShapePop
 * @data: a resolved #NcGalaxyShapePopData
 *
 * Gets the mode of $x=|\chi_I|^2$ for this galaxy's resolved population
 * density, i.e. where $P(x)$ peaks. This is a per-@data capability
 * (@data->ldata_get_mode_x), not a virtual method, mirroring
 * nc_galaxy_shape_pop_get_sigma(); unlike that capability, 0 is always a
 * meaningful default (the model is assumed radially symmetric about
 * $\chi_I=0$ unless it says otherwise), so this never errors.
 *
 * Returns: the mode of $x$, or 0 if the model does not override it.
 */
gdouble
nc_galaxy_shape_pop_get_mode_x (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  return (data->ldata_get_mode_x != NULL) ? data->ldata_get_mode_x (data) : 0.0;
}

