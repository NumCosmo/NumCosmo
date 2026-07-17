/***************************************************************************
 *            nc_galaxy_redshift_factor_composed.c
 *
 *  Tue Jul 1 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_factor_composed.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcGalaxyRedshiftFactorComposed:
 *
 * Composed redshift calculator scheme: a population convolved with an observable.
 *
 * Builds the per-galaxy joint density as the product of a population
 * distribution $P(z \mid I)$ (a #NcGalaxyRedshiftPop) and a photo-z observable
 * conditional $P(\mathrm{obs} \mid z)$ (a #NcGalaxyRedshiftObs), conditioned on
 * a photometric selection window $[z_{p,\min}, z_{p,\max}]$:
 * $$ p(\mathrm{obs}, z \mid I) = P(z \mid I)\,
 *    \frac{P(\mathrm{obs} \mid z)}{\int_{z_{p,\min}}^{z_{p,\max}}
 *    P(\mathrm{obs}' \mid z)\,\mathrm{d}\,\mathrm{obs}'}. $$
 * The window normalization is the observable's own selection mass
 * (nc_galaxy_redshift_obs_window_mass()), so the scheme never depends on the
 * observable's kernel shape. The scheme hands out this joint as an integrand;
 * the orchestrator performs the single $z$-integral.
 *
 * Following the calculator convention, the scheme does NOT hold the population
 * and observable models: they are resolved from the #NcmMSet passed to each
 * method (both are MAIN models), so a fitter varying them is seen immediately.
 * Only the selection window is scheme configuration.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_redshift_factor_composed.h"
#include "ncm/core/ncm_dtuple.h"
#include "ncm/core/ncm_rng.h"
#include "ncm/integration/ncm_integral1d_ptr.h"
#include "ncm/model/ncm_model.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcGalaxyRedshiftFactorComposedPrivate
{
  gdouble zp_min;
  gdouble zp_max;

  /* Cached by prepare(): a hash of the population/observable models' state,
   * so callers can detect a change without knowing which models are
   * involved -- see #NcGalaxyShapeFactor's analogous radius/optzs/pop
   * hashes for the full rationale. */
  guint64 hash;
} NcGalaxyRedshiftFactorComposedPrivate;

struct _NcGalaxyRedshiftFactorComposed
{
  NcGalaxyRedshiftFactor parent_instance;
};

enum
{
  PROP_0,
  PROP_ZP_LIM,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyRedshiftFactorComposed, nc_galaxy_redshift_factor_composed, NC_TYPE_GALAXY_REDSHIFT_FACTOR);

/* Per-galaxy fragment: embeds the observable model's own {obs} fragment. */
typedef struct _ComposedLData
{
  NcGalaxyRedshiftObsData *obs_data;
} ComposedLData;

/* The population and observable are never held; they are resolved from the mset
 * (both MAIN models) each time a method needs them. */
static void
_composed_peek_models (NcmMSet *mset, NcGalaxyRedshiftPop **population, NcGalaxyRedshiftObs **observable)
{
  *population = NC_GALAXY_REDSHIFT_POP (ncm_mset_peek (mset, nc_galaxy_redshift_pop_id ()));
  *observable = NC_GALAXY_REDSHIFT_OBS (ncm_mset_peek (mset, nc_galaxy_redshift_obs_id ()));

  g_assert_nonnull (*population);
  g_assert_nonnull (*observable);
}

static void
nc_galaxy_redshift_factor_composed_init (NcGalaxyRedshiftFactorComposed *gsdrc)
{
  NcGalaxyRedshiftFactorComposedPrivate * const self = nc_galaxy_redshift_factor_composed_get_instance_private (gsdrc);

  self->zp_min = 0.0;
  self->zp_max = 0.0;
  self->hash   = 0;
}

static void
_nc_galaxy_redshift_factor_composed_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyRedshiftFactorComposed *gsdrc = NC_GALAXY_REDSHIFT_FACTOR_COMPOSED (object);

  g_return_if_fail (NC_IS_GALAXY_REDSHIFT_FACTOR_COMPOSED (gsdrc));

  switch (prop_id)
  {
    case PROP_ZP_LIM:
    {
      NcmDTuple2 *lim = g_value_get_boxed (value);

      if (lim == NULL)
        g_error ("_nc_galaxy_redshift_factor_composed_set_property: zp-lim is NULL");

      nc_galaxy_redshift_factor_composed_set_zp_lim (gsdrc, lim->elements[0], lim->elements[1]);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_redshift_factor_composed_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyRedshiftFactorComposed *gsdrc              = NC_GALAXY_REDSHIFT_FACTOR_COMPOSED (object);
  NcGalaxyRedshiftFactorComposedPrivate * const self = nc_galaxy_redshift_factor_composed_get_instance_private (gsdrc);

  g_return_if_fail (NC_IS_GALAXY_REDSHIFT_FACTOR_COMPOSED (gsdrc));

  switch (prop_id)
  {
    case PROP_ZP_LIM:
      g_value_take_boxed (value, ncm_dtuple2_new (self->zp_min, self->zp_max));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_redshift_factor_composed_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_redshift_factor_composed_parent_class)->finalize (object);
}

static void _nc_galaxy_redshift_factor_composed_data_init (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data);
static void _nc_galaxy_redshift_factor_composed_gen (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng);
static gboolean _nc_galaxy_redshift_factor_composed_gen1 (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng);
static void _nc_galaxy_redshift_factor_composed_prepare (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset);
static guint64 _nc_galaxy_redshift_factor_composed_get_hash (NcGalaxyRedshiftFactor *gsdr);
static NcGalaxyRedshiftFactorIntegrand *_nc_galaxy_redshift_factor_composed_integ (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, gboolean use_lnp);
static void _nc_galaxy_redshift_factor_composed_get_integ_lim (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble *z_min, gdouble *z_max);
static gdouble _nc_galaxy_redshift_factor_composed_norm (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data);
static NcmIntegralFixed *_nc_galaxy_redshift_factor_composed_make_fixed_nodes (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble z_lo, gdouble z_hi, guint n_nodes, guint rule_n);

static void
nc_galaxy_redshift_factor_composed_class_init (NcGalaxyRedshiftFactorComposedClass *klass)
{
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcGalaxyRedshiftFactorClass *gsdr_class = NC_GALAXY_REDSHIFT_FACTOR_CLASS (klass);

  object_class->set_property = &_nc_galaxy_redshift_factor_composed_set_property;
  object_class->get_property = &_nc_galaxy_redshift_factor_composed_get_property;
  object_class->finalize     = &_nc_galaxy_redshift_factor_composed_finalize;

  /**
   * NcGalaxyRedshiftFactorComposed:zp-lim:
   *
   * The photometric selection window $[z_{p,\min}, z_{p,\max}]$.
   */
  g_object_class_install_property (object_class,
                                   PROP_ZP_LIM,
                                   g_param_spec_boxed ("zp-lim",
                                                       NULL,
                                                       "Photometric selection window",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  gsdr_class->data_init        = &_nc_galaxy_redshift_factor_composed_data_init;
  gsdr_class->gen              = &_nc_galaxy_redshift_factor_composed_gen;
  gsdr_class->gen1             = &_nc_galaxy_redshift_factor_composed_gen1;
  gsdr_class->prepare          = &_nc_galaxy_redshift_factor_composed_prepare;
  gsdr_class->integ            = &_nc_galaxy_redshift_factor_composed_integ;
  gsdr_class->get_integ_lim    = &_nc_galaxy_redshift_factor_composed_get_integ_lim;
  gsdr_class->norm             = &_nc_galaxy_redshift_factor_composed_norm;
  gsdr_class->make_fixed_nodes = &_nc_galaxy_redshift_factor_composed_make_fixed_nodes;
  gsdr_class->get_hash         = &_nc_galaxy_redshift_factor_composed_get_hash;
}

/* Fragment vtable ------------------------------------------------------- */

static void
_composed_ldata_free (gpointer ldata)
{
  ComposedLData *cldata = (ComposedLData *) ldata;

  nc_galaxy_redshift_obs_data_unref (cldata->obs_data);
  g_free (cldata);
}

static void
_composed_ldata_read_row (NcGalaxyRedshiftFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  ComposedLData *cldata = (ComposedLData *) data->ldata;

  nc_galaxy_redshift_obs_data_read_row (cldata->obs_data, obs, i);
}

static void
_composed_ldata_write_row (NcGalaxyRedshiftFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  ComposedLData *cldata = (ComposedLData *) data->ldata;

  nc_galaxy_redshift_obs_data_write_row (cldata->obs_data, obs, i);
}

static void
_composed_ldata_required_columns (NcGalaxyRedshiftFactorData *data, GList **columns)
{
  ComposedLData *cldata = (ComposedLData *) data->ldata;
  GList *obs_columns    = nc_galaxy_redshift_obs_data_required_columns (cldata->obs_data);

  *columns = g_list_concat (*columns, obs_columns);
}

static void
_nc_galaxy_redshift_factor_composed_data_init (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data)
{
  ComposedLData *cldata = g_new0 (ComposedLData, 1);
  NcGalaxyRedshiftPop *population;
  NcGalaxyRedshiftObs *observable;

  _composed_peek_models (mset, &population, &observable);
  cldata->obs_data = nc_galaxy_redshift_obs_data_new (observable);

  data->ldata                  = cldata;
  data->ldata_destroy          = &_composed_ldata_free;
  data->ldata_read_row         = &_composed_ldata_read_row;
  data->ldata_write_row        = &_composed_ldata_write_row;
  data->ldata_required_columns = &_composed_ldata_required_columns;
}

/* Sampling -------------------------------------------------------------- */

static void
_nc_galaxy_redshift_factor_composed_gen (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng)
{
  NcGalaxyRedshiftFactorComposed *gsdrc              = NC_GALAXY_REDSHIFT_FACTOR_COMPOSED (gsdr);
  NcGalaxyRedshiftFactorComposedPrivate * const self = nc_galaxy_redshift_factor_composed_get_instance_private (gsdrc);
  ComposedLData *cldata                              = (ComposedLData *) data->ldata;
  guint max_iter                                     = 1000;
  NcGalaxyRedshiftPop *population;
  NcGalaxyRedshiftObs *observable;
  gdouble z, zp;

  _composed_peek_models (mset, &population, &observable);

  do {
    z  = nc_galaxy_redshift_pop_gen (population, rng);
    zp = nc_galaxy_redshift_obs_gen (observable, cldata->obs_data, z, rng);

    if (max_iter-- == 0)
      g_error ("_nc_galaxy_redshift_factor_composed_gen: maximum number of iterations reached.");
  } while ((zp > self->zp_max) || (zp < self->zp_min));

  data->z = z;
}

static gboolean
_nc_galaxy_redshift_factor_composed_gen1 (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng)
{
  NcGalaxyRedshiftFactorComposed *gsdrc              = NC_GALAXY_REDSHIFT_FACTOR_COMPOSED (gsdr);
  NcGalaxyRedshiftFactorComposedPrivate * const self = nc_galaxy_redshift_factor_composed_get_instance_private (gsdrc);
  ComposedLData *cldata                              = (ComposedLData *) data->ldata;
  NcGalaxyRedshiftPop *population;
  NcGalaxyRedshiftObs *observable;
  gdouble z, zp;

  _composed_peek_models (mset, &population, &observable);

  z  = nc_galaxy_redshift_pop_gen (population, rng);
  zp = nc_galaxy_redshift_obs_gen (observable, cldata->obs_data, z, rng);

  data->z = z;

  return (zp <= self->zp_max) && (zp >= self->zp_min);
}

static void
_nc_galaxy_redshift_factor_composed_prepare (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset)
{
  NcGalaxyRedshiftFactorComposed *gsdrc              = NC_GALAXY_REDSHIFT_FACTOR_COMPOSED (gsdr);
  NcGalaxyRedshiftFactorComposedPrivate * const self = nc_galaxy_redshift_factor_composed_get_instance_private (gsdrc);
  NcGalaxyRedshiftPop *population;
  NcGalaxyRedshiftObs *observable;

  /* Validate the required models are present in the mset. */
  _composed_peek_models (mset, &population, &observable);

  /* Neither model's state is otherwise cached (both are always resolved
   * fresh from mset per-call, per the class doc), so there is nothing else
   * to refresh here -- just keep the hash current for get_hash(). */
  self->hash = ncm_model_state_get_pkey (NCM_MODEL (population)) ^
               (ncm_model_state_get_pkey (NCM_MODEL (observable)) * 31U);
}

static guint64
_nc_galaxy_redshift_factor_composed_get_hash (NcGalaxyRedshiftFactor *gsdr)
{
  NcGalaxyRedshiftFactorComposed *gsdrc              = NC_GALAXY_REDSHIFT_FACTOR_COMPOSED (gsdr);
  NcGalaxyRedshiftFactorComposedPrivate * const self = nc_galaxy_redshift_factor_composed_get_instance_private (gsdrc);

  return self->hash;
}

/* Integrand ------------------------------------------------------------- */

struct _IntegData
{
  NcGalaxyRedshiftFactorComposed *gsdrc;
  NcGalaxyRedshiftPop *population;
  NcGalaxyRedshiftObs *observable;
};

/* LCOV_EXCL_START */
static gpointer
_integ_data_copy (gpointer idata)
{
  struct _IntegData *new_idata = g_new0 (struct _IntegData, 1);

  *new_idata = *(struct _IntegData *) idata;

  return new_idata;
}

/* LCOV_EXCL_STOP */

static void
_integ_data_free (gpointer idata)
{
  g_free (idata);
}

static gdouble
_composed_integ_f (gpointer callback_data, const gdouble z, NcGalaxyRedshiftFactorData *data)
{
  const struct _IntegData *int_data                  = (struct _IntegData *) callback_data;
  NcGalaxyRedshiftFactorComposedPrivate * const self = nc_galaxy_redshift_factor_composed_get_instance_private (int_data->gsdrc);
  ComposedLData *cldata                              = (ComposedLData *) data->ldata;
  const gdouble Pz                                   = nc_galaxy_redshift_pop_eval (int_data->population, z);
  const gdouble Pobs                                 = nc_galaxy_redshift_obs_eval (int_data->observable, cldata->obs_data, z);
  const gdouble Wnorm                                = nc_galaxy_redshift_obs_window_mass (int_data->observable, cldata->obs_data, z, self->zp_min, self->zp_max);

  return Pz * Pobs / Wnorm;
}

static gdouble
_composed_ln_integ_f (gpointer callback_data, const gdouble z, NcGalaxyRedshiftFactorData *data)
{
  const struct _IntegData *int_data                  = (struct _IntegData *) callback_data;
  NcGalaxyRedshiftFactorComposedPrivate * const self = nc_galaxy_redshift_factor_composed_get_instance_private (int_data->gsdrc);
  ComposedLData *cldata                              = (ComposedLData *) data->ldata;
  const gdouble ln_Pz                                = nc_galaxy_redshift_pop_ln_eval (int_data->population, z);
  const gdouble Pobs                                 = nc_galaxy_redshift_obs_eval (int_data->observable, cldata->obs_data, z);
  const gdouble Wnorm                                = nc_galaxy_redshift_obs_window_mass (int_data->observable, cldata->obs_data, z, self->zp_min, self->zp_max);

  return ln_Pz + log (Pobs) - log (Wnorm);
}

static NcGalaxyRedshiftFactorIntegrand *
_nc_galaxy_redshift_factor_composed_integ (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, gboolean use_lnp)
{
  NcGalaxyRedshiftFactorComposed *gsdrc  = NC_GALAXY_REDSHIFT_FACTOR_COMPOSED (gsdr);
  struct _IntegData *int_data            = g_new0 (struct _IntegData, 1);
  NcGalaxyRedshiftFactorIntegrand *integ = nc_galaxy_redshift_factor_integrand_new (use_lnp ? _composed_ln_integ_f : _composed_integ_f,
                                                                                    _integ_data_free,
                                                                                    _integ_data_copy,
                                                                                    NULL,
                                                                                    int_data);

  int_data->gsdrc = gsdrc;
  _composed_peek_models (mset, &int_data->population, &int_data->observable);

  return integ;
}

/* Integration limits and normalization ---------------------------------- */

static void
_nc_galaxy_redshift_factor_composed_get_integ_lim (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble *z_min, gdouble *z_max)
{
  ComposedLData *cldata = (ComposedLData *) data->ldata;
  NcGalaxyRedshiftPop *population;
  NcGalaxyRedshiftObs *observable;
  gdouble zi, zf, obs_zi, obs_zf;

  _composed_peek_models (mset, &population, &observable);

  nc_galaxy_redshift_pop_get_lim (population, &zi, &zf);
  nc_galaxy_redshift_obs_get_true_z_lim (observable, cldata->obs_data, &obs_zi, &obs_zf);

  *z_min = MAX (zi, obs_zi);
  *z_max = MIN (zf, obs_zf);
}

typedef struct _ComposedGSLArg
{
  struct _IntegData int_data;
  NcGalaxyRedshiftFactorData *data;
} ComposedGSLArg;

static gdouble
_composed_norm_integ1d_f (gpointer user_data, const gdouble z, const gdouble w)
{
  ComposedGSLArg *arg = (ComposedGSLArg *) user_data;

  return _composed_integ_f (&arg->int_data, z, arg->data);
}

static gdouble
_nc_galaxy_redshift_factor_composed_norm (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data)
{
  NcGalaxyRedshiftFactorComposed *gsdrc = NC_GALAXY_REDSHIFT_FACTOR_COMPOSED (gsdr);
  ComposedGSLArg arg;
  NcmIntegral1dPtr *integrator;
  gdouble z_min, z_max, err, result;

  nc_galaxy_redshift_factor_get_integ_lim (gsdr, mset, data, &z_min, &z_max);

  arg.int_data.gsdrc = gsdrc;
  _composed_peek_models (mset, &arg.int_data.population, &arg.int_data.observable);
  arg.data = data;

  integrator = ncm_integral1d_ptr_new (&_composed_norm_integ1d_f, NULL);
  ncm_integral1d_ptr_set_userdata (integrator, &arg);
  result = ncm_integral1d_eval (NCM_INTEGRAL1D (integrator), z_min, z_max, &err);
  ncm_integral1d_ptr_free (integrator);

  return result;
}

static gdouble
_composed_fixed_nodes_gsl_f (gdouble z, gpointer user_data)
{
  ComposedGSLArg *arg = (ComposedGSLArg *) user_data;

  return _composed_integ_f (&arg->int_data, z, arg->data);
}

static NcmIntegralFixed *
_nc_galaxy_redshift_factor_composed_make_fixed_nodes (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble z_lo, gdouble z_hi, guint n_nodes, guint rule_n)
{
  NcGalaxyRedshiftFactorComposed *gsdrc = NC_GALAXY_REDSHIFT_FACTOR_COMPOSED (gsdr);
  ComposedGSLArg arg;
  NcmIntegralFixed *intf;
  gsl_function F;

  intf = ncm_integral_fixed_new (n_nodes, rule_n, z_lo, z_hi);

  arg.int_data.gsdrc = gsdrc;
  _composed_peek_models (mset, &arg.int_data.population, &arg.int_data.observable);
  arg.data = data;

  F.function = &_composed_fixed_nodes_gsl_f;
  F.params   = &arg;

  ncm_integral_fixed_calc_nodes (intf, &F);

  return intf;
}

/* Public API ------------------------------------------------------------ */

/**
 * nc_galaxy_redshift_factor_composed_new:
 * @zp_min: the minimum photometric redshift of the selection window
 * @zp_max: the maximum photometric redshift of the selection window
 *
 * Creates a new #NcGalaxyRedshiftFactorComposed calculator scheme over the
 * selection window [@zp_min, @zp_max]. The population and observable models are
 * not held; they are resolved from the #NcmMSet passed to each method.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftFactorComposed.
 */
NcGalaxyRedshiftFactorComposed *
nc_galaxy_redshift_factor_composed_new (const gdouble zp_min, const gdouble zp_max)
{
  NcmDTuple2 lim                        = NCM_DTUPLE2_STATIC_INIT (zp_min, zp_max);
  NcGalaxyRedshiftFactorComposed *gsdrc = g_object_new (NC_TYPE_GALAXY_REDSHIFT_FACTOR_COMPOSED,
                                                        "zp-lim", &lim,
                                                        NULL);

  return gsdrc;
}

/**
 * nc_galaxy_redshift_factor_composed_ref:
 * @gsdrc: a #NcGalaxyRedshiftFactorComposed
 *
 * Increases the reference count of @gsdrc by one.
 *
 * Returns: (transfer full): @gsdrc.
 */
NcGalaxyRedshiftFactorComposed *
nc_galaxy_redshift_factor_composed_ref (NcGalaxyRedshiftFactorComposed *gsdrc)
{
  return g_object_ref (gsdrc);
}

/**
 * nc_galaxy_redshift_factor_composed_free:
 * @gsdrc: a #NcGalaxyRedshiftFactorComposed
 *
 * Decreases the reference count of @gsdrc by one.
 *
 */
void
nc_galaxy_redshift_factor_composed_free (NcGalaxyRedshiftFactorComposed *gsdrc)
{
  g_object_unref (gsdrc);
}

/**
 * nc_galaxy_redshift_factor_composed_clear:
 * @gsdrc: a #NcGalaxyRedshiftFactorComposed
 *
 * Decreases the reference count of @gsdrc by one, and sets the pointer *@gsdrc
 * to NULL.
 *
 */
void
nc_galaxy_redshift_factor_composed_clear (NcGalaxyRedshiftFactorComposed **gsdrc)
{
  g_clear_object (gsdrc);
}

/**
 * nc_galaxy_redshift_factor_composed_set_zp_lim:
 * @gsdrc: a #NcGalaxyRedshiftFactorComposed
 * @zp_min: the minimum photometric redshift
 * @zp_max: the maximum photometric redshift
 *
 * Sets the photometric selection window [@zp_min, @zp_max].
 *
 */
void
nc_galaxy_redshift_factor_composed_set_zp_lim (NcGalaxyRedshiftFactorComposed *gsdrc, const gdouble zp_min, const gdouble zp_max)
{
  NcGalaxyRedshiftFactorComposedPrivate * const self = nc_galaxy_redshift_factor_composed_get_instance_private (gsdrc);

  g_assert_cmpfloat (zp_min, <, zp_max);

  self->zp_min = zp_min;
  self->zp_max = zp_max;
}

/**
 * nc_galaxy_redshift_factor_composed_get_zp_lim:
 * @gsdrc: a #NcGalaxyRedshiftFactorComposed
 * @zp_min: (out): the minimum photometric redshift
 * @zp_max: (out): the maximum photometric redshift
 *
 * Gets the photometric selection window.
 *
 */
void
nc_galaxy_redshift_factor_composed_get_zp_lim (NcGalaxyRedshiftFactorComposed *gsdrc, gdouble *zp_min, gdouble *zp_max)
{
  NcGalaxyRedshiftFactorComposedPrivate * const self = nc_galaxy_redshift_factor_composed_get_instance_private (gsdrc);

  *zp_min = self->zp_min;
  *zp_max = self->zp_max;
}

