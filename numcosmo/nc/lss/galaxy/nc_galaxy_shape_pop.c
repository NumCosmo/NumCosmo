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
 * ellipticity modulus $r = |\chi_I| \in [0,1)$, i.e. its own natural
 * r-marginal density $P_\mathrm{pop}(r)$, normalized so
 * $\int_0^1 P_\mathrm{pop}(r)\,\mathrm{d}r=1$. It is a pure parametric law
 * for which the shape calculator (#NcGalaxyWLShapeCalc) performs the exact
 * marginalization over the intrinsic ellipticity.
 *
 * $P_\mathrm{pop}(r)$ is NOT a 2D area density: a consumer needing the
 * density with respect to $\mathrm{d}^2\chi_I$ (e.g. #NcGalaxyShapeFactorFixedQuad,
 * #NcGalaxyShapeFactorQuad) computes $P_\mathrm{pop}(r)/(2\pi r)$ itself --
 * a genuine, unavoidable pointwise divergence at $r=0$ for any population
 * whose $P_\mathrm{pop}(r)$ does not vanish at least linearly there (e.g.
 * #NcGalaxyShapePopBeta with $\alpha<2$), but a property of those
 * consumers' own (non-polar-in-$\chi_I$) quadrature, not of
 * $P_\mathrm{pop}(r)$ itself (always bounded), nor of the exact 2D
 * integral (stays finite through $r=0$: an integrable singularity,
 * manifestly so in genuine polar coordinates, where the $r$ in
 * $r\,\mathrm{d}r\,\mathrm{d}\theta$ cancels it exactly).
 *
 * Each concrete model implements its density directly through the
 * nc_galaxy_shape_pop_eval_p() virtual: there is no shared factorization or
 * quadrature scheme imposed at this level, each subclass is simply
 * responsible for its own normalized $P_\mathrm{pop}(r)$. Following the
 * #NcGalaxyShapeFactorData /
 * #NcGalaxyPositionFactorData idiom, the resolved per-galaxy state lives in a
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

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_integration.h>
#endif /* NUMCOSMO_GIR_SCAN */

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
_nc_galaxy_shape_pop_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble r)
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
                                           const NcmLaurentSeriesTPS *x_series, NcmLaurentSeriesTPS *out)
{
  g_error ("_nc_galaxy_shape_pop_eval_p_rho2_g_series: method not implemented.");
}

/* LCOV_EXCL_STOP */

/* Real (non-stub) default: plain loop calling eval_p() through the class
 * pointer once (not nc_galaxy_shape_pop_eval_p(), avoiding one extra level
 * of GET_CLASS per element). Correct for every subclass, so this is a safe
 * fallback for populations that do not override eval_p_array(). */
static void
_nc_galaxy_shape_pop_eval_p_array (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const GArray *r, GArray **p)
{
  NcGalaxyShapePopClass *klass = NC_GALAXY_SHAPE_POP_GET_CLASS (gsp);
  gdouble *p_data;
  guint i;

  if (*p == NULL)
    *p = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), r->len);

  g_array_set_size (*p, r->len);

  p_data = (gdouble *) (*p)->data;

  for (i = 0; i < r->len; i++)
    p_data[i] = klass->eval_p (gsp, data, g_array_index (r, gdouble, i));
}

/* Number of fixed Gauss-Legendre nodes for the generic moment_2k() fallback.
 * This path is only exercised by a population with no closed-form override
 * (Gauss, GaussLocal and Beta all have one), so it is sized for headroom
 * rather than hot-loop cost. */
#define NC_GALAXY_SHAPE_POP_MOMENT_2K_DEFAULT_NNODES 64

/* Real (non-stub) default: fixed Gauss-Legendre quadrature of
 * int_0^1 r^2k P_pop(r) dr, batched through eval_p_array() (r^2k is
 * evaluated on the node array, not the population, so this is the second
 * safe generic fallback in the class, after eval_p_array() itself). k=0 is
 * the population's own normalization constant, returned as 1.0 directly
 * rather than integrated -- quadrature there would be a pure error source
 * for a quantity that is exact by construction. */
static gdouble
_nc_galaxy_shape_pop_moment_2k (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const guint k)
{
  const guint n_nodes = NC_GALAXY_SHAPE_POP_MOMENT_2K_DEFAULT_NNODES;

  if (k == 0)
  {
    return 1.0;
  }
  else
  {
    gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc (n_nodes);
    gdouble weight[NC_GALAXY_SHAPE_POP_MOMENT_2K_DEFAULT_NNODES];
    GArray *r_arr  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n_nodes);
    GArray *p_arr  = NULL;
    gdouble moment = 0.0;
    guint i;

    g_array_set_size (r_arr, n_nodes);

    {
      gdouble * const r_data = (gdouble *) r_arr->data;

      for (i = 0; i < n_nodes; i++)
        gsl_integration_glfixed_point (0.0, 1.0, i, &r_data[i], &weight[i], table);
    }

    nc_galaxy_shape_pop_eval_p_array (gsp, data, r_arr, &p_arr);

    {
      const gdouble * const r_data = (const gdouble *) r_arr->data;
      const gdouble * const p_data = (const gdouble *) p_arr->data;

      for (i = 0; i < n_nodes; i++)
        moment += weight[i] * pow (r_data[i], 2.0 * k) * p_data[i];
    }

    g_array_unref (r_arr);
    g_array_unref (p_arr);
    gsl_integration_glfixed_table_free (table);

    return moment;
  }
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
  klass->gen                  = &_nc_galaxy_shape_pop_gen;
  klass->e_rms                = &_nc_galaxy_shape_pop_e_rms;
  klass->eval_p_rho2_g_series = &_nc_galaxy_shape_pop_eval_p_rho2_g_series;
  klass->eval_p_array         = &_nc_galaxy_shape_pop_eval_p_array;
  klass->moment_2k            = &_nc_galaxy_shape_pop_moment_2k;
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
  data->ldata_get_mode_r       = NULL;

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
 * that nc_galaxy_shape_pop_eval_p() reads.
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
 * @r: the ellipticity modulus $r = |\chi_I| \in [0,1)$
 *
 * Evaluates the population's own natural r-marginal density
 * $P_\mathrm{pop}(r)$, normalized so $\int_0^1 P_\mathrm{pop}(r)\,\mathrm{d}r=1$,
 * reading the resolved parameters from @data (no live parameter access).
 * This is NOT a 2D area density -- a consumer needing the density with
 * respect to $\mathrm{d}^2\chi_I$ computes $P_\mathrm{pop}(r)/(2\pi r)$
 * itself (see the class docs for why that division, not this vfunc, is
 * where any $r=0$ divergence for a population like #NcGalaxyShapePopBeta
 * with $\alpha<2$ actually lives).
 *
 * Returns: the density $P_\mathrm{pop}(r)$.
 */
gdouble
nc_galaxy_shape_pop_eval_p (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const gdouble r)
{
  return NC_GALAXY_SHAPE_POP_GET_CLASS (gsp)->eval_p (gsp, data, r);
}

/**
 * nc_galaxy_shape_pop_eval_p_rho2_g_series:
 * @gsp: a #NcGalaxyShapePop
 * @data: a resolved #NcGalaxyShapePopData
 * @x_series: $x(g)=|\chi_I(\chi_L,g)|^2$'s own $g$-Taylor coefficients,
 * population-independent
 * @out: this population's $P(x(g))=\mathrm{eval\_p}(x(g))$ composition, as
 * $g$-Taylor coefficients (same order as @x_series)
 *
 * Taylor-in-$g$ analog of nc_galaxy_shape_pop_eval_p(): composes this
 * population's own fully normalized density with the (already computed,
 * population-independent) shear-map series $x(g)$, order by order in $g$.
 * There is no generic default -- every subclass used with
 * #NcGalaxyShapeFactorSeriesLensed must provide its own implementation; the
 * base class errors clearly otherwise.
 */
void
nc_galaxy_shape_pop_eval_p_rho2_g_series (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data,
                                          const NcmLaurentSeriesTPS *x_series, NcmLaurentSeriesTPS *out)
{
  NC_GALAXY_SHAPE_POP_GET_CLASS (gsp)->eval_p_rho2_g_series (gsp, data, x_series, out);
}

/**
 * nc_galaxy_shape_pop_eval_p_array: (virtual eval_p_array)
 * @gsp: a #NcGalaxyShapePop
 * @data: a resolved #NcGalaxyShapePopData
 * @r: (element-type gdouble): array of $r=|\chi_I|$ values to evaluate
 * @p: (out callee-allocates) (transfer full) (element-type gdouble): output
 * array of $P_\mathrm{pop}(r)$ values, same length as @r
 *
 * Batched form of nc_galaxy_shape_pop_eval_p(): evaluates the population
 * density at every element of @r in one call.
 *
 * The @p parameter supports two usage patterns:
 *
 * - **Python/convenience usage**: pass @p pointing to NULL (`*p == NULL`).
 *   A new #GArray is allocated and returned. The `(out callee-allocates)`
 *   annotation ensures Python bindings automatically use this mode.
 *
 * - **C optimization**: pass @p pointing to a pre-allocated #GArray
 *   (`*p != NULL`). The existing array is resized and refilled, avoiding
 *   repeated allocation/deallocation in performance-critical loops (e.g.
 *   #NcGalaxyShapeFactorFixedQuad, which reuses the same @p across every
 *   $g$ a fit tries).
 *
 * The generic default just loops calling eval_p(), so this is always safe
 * to call; a subclass may override it to batch the work (see
 * #NcGalaxyShapePopBeta).
 */
void
nc_galaxy_shape_pop_eval_p_array (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const GArray *r, GArray **p)
{
  NC_GALAXY_SHAPE_POP_GET_CLASS (gsp)->eval_p_array (gsp, data, r, p);
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
 * nc_galaxy_shape_pop_exponent_at_origin: (virtual exponent_at_origin)
 * @gsp: a #NcGalaxyShapePop
 *
 * Determines how the population density behaves at the origin.
 *
 * Returns: the exponent $\alpha_o$ such that $P_\mathrm{pop}(r) \sim r^{\alpha_o}$ as $r \to 0$.
 */
gdouble
nc_galaxy_shape_pop_exponent_at_origin (NcGalaxyShapePop *gsp)
{
  return NC_GALAXY_SHAPE_POP_GET_CLASS (gsp)->exponent_at_origin (gsp);
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
 * nc_galaxy_shape_pop_get_mode_r:
 * @gsp: a #NcGalaxyShapePop
 * @data: a resolved #NcGalaxyShapePopData
 *
 * Gets the mode of $r=|\chi_I|$ for this galaxy's resolved population
 * density, i.e. where $P_\mathrm{pop}(r)$ peaks. This is a per-@data
 * capability (@data->ldata_get_mode_r), not a virtual method, mirroring
 * nc_galaxy_shape_pop_get_sigma(); unlike that capability, 0 is always a
 * meaningful default (the model is assumed radially symmetric about
 * $\chi_I=0$ unless it says otherwise), so this never errors.
 *
 * Returns: the mode of $r$, or 0 if the model does not override it.
 */
gdouble
nc_galaxy_shape_pop_get_mode_r (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  return (data->ldata_get_mode_r != NULL) ? data->ldata_get_mode_r (data) : 0.0;
}

/**
 * nc_galaxy_shape_pop_moment_2k: (virtual moment_2k)
 * @gsp: a #NcGalaxyShapePop
 * @data: a resolved #NcGalaxyShapePopData
 * @k: the moment order
 *
 * Computes the radial moment
 * $M_{2k} = \langle |\chi_I|^{2k} \rangle = \int_0^1 r^{2k}\,P_\mathrm{pop}(r)\,\mathrm{d}r$,
 * reading the resolved parameters from @data (the
 * nc_galaxy_shape_pop_eval_p() contract: no live parameter access), so
 * @data must have been resolved by nc_galaxy_shape_pop_prepare() since the
 * population's parameters last changed. $k=0$ always returns exactly 1
 * (the population's own normalization), never integrated.
 *
 * The base class provides a real (non-error) fixed-quadrature default valid
 * for any population; concrete models override it with a closed form where
 * one exists.
 *
 * Returns: the radial moment $M_{2k}$.
 */
gdouble
nc_galaxy_shape_pop_moment_2k (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data, const guint k)
{
  return NC_GALAXY_SHAPE_POP_GET_CLASS (gsp)->moment_2k (gsp, data, k);
}

