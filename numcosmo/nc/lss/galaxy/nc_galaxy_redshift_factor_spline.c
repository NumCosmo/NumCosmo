/***************************************************************************
 *            nc_galaxy_redshift_factor_spline.c
 *
 *  Sun Jul 13 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_factor_spline.c
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
 * NcGalaxyRedshiftFactorSpline:
 *
 * Redshift calculator scheme for a per-galaxy $p(z)$ given as data.
 *
 * Each galaxy carries its own redshift density as a #NcGalaxyPz, set with
 * nc_galaxy_redshift_factor_spline_data_set() and stored on the #NcGalaxyWLObs row with
 * nc_galaxy_wl_obs_set_pz()/nc_galaxy_wl_obs_peek_pz(). The scheme has no model: the
 * integrand, the integration limits, the normalization and the draws are those of the
 * #NcGalaxyPz of the galaxy, see #NcGalaxyPzBSplineSqrt for the density of a table
 * $(z_j, p_j)$.
 *
 * Reading a row that has no $p(z)$ aborts.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_redshift_factor_spline.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcGalaxyRedshiftFactorSplinePrivate
{
  gint placeholder;
} NcGalaxyRedshiftFactorSplinePrivate;

struct _NcGalaxyRedshiftFactorSpline
{
  NcGalaxyRedshiftFactor parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyRedshiftFactorSpline, nc_galaxy_redshift_factor_spline, NC_TYPE_GALAXY_REDSHIFT_FACTOR);

/* Per-galaxy fragment: the redshift density of the galaxy, referenced */
typedef struct _SplineLData
{
  NcGalaxyPz *pz;
} SplineLData;

static void
nc_galaxy_redshift_factor_spline_init (NcGalaxyRedshiftFactorSpline *gsdrs)
{
}

static void
_nc_galaxy_redshift_factor_spline_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_redshift_factor_spline_parent_class)->finalize (object);
}

static void _nc_galaxy_redshift_factor_spline_data_init (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data);
static void _nc_galaxy_redshift_factor_spline_gen (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng);
static gboolean _nc_galaxy_redshift_factor_spline_gen1 (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng);
static NcGalaxyRedshiftFactorIntegrand *_nc_galaxy_redshift_factor_spline_integ (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, gboolean use_lnp);
static void _nc_galaxy_redshift_factor_spline_get_integ_lim (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble *z_min, gdouble *z_max);
static gdouble _nc_galaxy_redshift_factor_spline_norm (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data);
static NcmIntegralFixed *_nc_galaxy_redshift_factor_spline_make_fixed_nodes (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble z_lo, gdouble z_hi, guint n_nodes, guint rule_n);

static void
nc_galaxy_redshift_factor_spline_class_init (NcGalaxyRedshiftFactorSplineClass *klass)
{
  GObjectClass *object_class              = G_OBJECT_CLASS (klass);
  NcGalaxyRedshiftFactorClass *gsdr_class = NC_GALAXY_REDSHIFT_FACTOR_CLASS (klass);

  object_class->finalize = &_nc_galaxy_redshift_factor_spline_finalize;

  gsdr_class->data_init        = &_nc_galaxy_redshift_factor_spline_data_init;
  gsdr_class->gen              = &_nc_galaxy_redshift_factor_spline_gen;
  gsdr_class->gen1             = &_nc_galaxy_redshift_factor_spline_gen1;
  gsdr_class->integ            = &_nc_galaxy_redshift_factor_spline_integ;
  gsdr_class->get_integ_lim    = &_nc_galaxy_redshift_factor_spline_get_integ_lim;
  gsdr_class->norm             = &_nc_galaxy_redshift_factor_spline_norm;
  gsdr_class->make_fixed_nodes = &_nc_galaxy_redshift_factor_spline_make_fixed_nodes;

  /*
   * prepare, get_hash and update_data keep the base no-ops: the density depends on no
   * model, and a NcGalaxyPz does not change after construction.
   */
}

/* Fragment vtable --------------------------------------------------------- */

static void
_spline_ldata_free (gpointer ldata)
{
  SplineLData *sldata = (SplineLData *) ldata;

  nc_galaxy_pz_clear (&sldata->pz);
  g_free (sldata);
}

/* References @pz in place of the density the fragment had */
static void
_spline_ldata_set_pz (SplineLData *sldata, NcGalaxyPz *pz)
{
  g_assert_nonnull (pz);

  nc_galaxy_pz_ref (pz);
  nc_galaxy_pz_clear (&sldata->pz);

  sldata->pz = pz;
}

static void
_spline_ldata_read_row (NcGalaxyRedshiftFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  _spline_ldata_set_pz (sldata, nc_galaxy_wl_obs_peek_pz (obs, i));
}

static void
_spline_ldata_write_row (NcGalaxyRedshiftFactorData *data, NcGalaxyWLObs *obs, const guint i)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  nc_galaxy_wl_obs_set_pz (obs, i, sldata->pz);
}

static void
_spline_ldata_required_columns (NcGalaxyRedshiftFactorData *data, GList **columns)
{
  /* The density is carried on the pz slot of the NcGalaxyWLObs row, not on a column */
}

static void
_nc_galaxy_redshift_factor_spline_data_init (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data)
{
  SplineLData *sldata = g_new0 (SplineLData, 1);

  data->ldata                  = sldata;
  data->ldata_destroy          = &_spline_ldata_free;
  data->ldata_read_row         = &_spline_ldata_read_row;
  data->ldata_write_row        = &_spline_ldata_write_row;
  data->ldata_required_columns = &_spline_ldata_required_columns;
}

/* Sampling ----------------------------------------------------------------- */

static void
_nc_galaxy_redshift_factor_spline_gen (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  data->z = nc_galaxy_pz_gen (sldata->pz, rng);
}

static gboolean
_nc_galaxy_redshift_factor_spline_gen1 (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng)
{
  _nc_galaxy_redshift_factor_spline_gen (gsdr, mset, data, rng);

  return TRUE;
}

/* Integrand ---------------------------------------------------------------- */

static gdouble
_spline_integ_f (gpointer callback_data, const gdouble z, NcGalaxyRedshiftFactorData *data)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  return nc_galaxy_pz_eval_p (sldata->pz, z);
}

static gdouble
_spline_ln_integ_f (gpointer callback_data, const gdouble z, NcGalaxyRedshiftFactorData *data)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  return nc_galaxy_pz_eval_ln_p (sldata->pz, z);
}

/* No closure payload needed. */
static gpointer
_spline_integ_data_copy (gpointer callback_data)
{
  return NULL;
}

static void
_spline_integ_data_free (gpointer callback_data)
{
}

static NcGalaxyRedshiftFactorIntegrand *
_nc_galaxy_redshift_factor_spline_integ (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, gboolean use_lnp)
{
  return nc_galaxy_redshift_factor_integrand_new (use_lnp ? _spline_ln_integ_f : _spline_integ_f,
                                                  &_spline_integ_data_free,
                                                  &_spline_integ_data_copy,
                                                  NULL,
                                                  NULL);
}

/* Integration limits and normalization -------------------------------------- */

static void
_nc_galaxy_redshift_factor_spline_get_integ_lim (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble *z_min, gdouble *z_max)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  nc_galaxy_pz_get_range (sldata->pz, z_min, z_max);
}

static gdouble
_nc_galaxy_redshift_factor_spline_norm (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  return nc_galaxy_pz_get_norm (sldata->pz);
}

static gdouble
_spline_fixed_nodes_gsl_f (gdouble z, gpointer user_data)
{
  NcGalaxyRedshiftFactorData *data = (NcGalaxyRedshiftFactorData *) user_data;

  return _spline_integ_f (NULL, z, data);
}

static NcmIntegralFixed *
_nc_galaxy_redshift_factor_spline_make_fixed_nodes (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble z_lo, gdouble z_hi, guint n_nodes, guint rule_n)
{
  NcmIntegralFixed *intf;
  gsl_function F;

  intf = ncm_integral_fixed_new (n_nodes, rule_n, z_lo, z_hi);

  F.function = &_spline_fixed_nodes_gsl_f;
  F.params   = data;

  ncm_integral_fixed_calc_nodes (intf, &F);

  return intf;
}

/* Public API ----------------------------------------------------------------- */

/**
 * nc_galaxy_redshift_factor_spline_new:
 *
 * Creates a new #NcGalaxyRedshiftFactorSpline calculator scheme.
 *
 * Returns: (transfer full): a new #NcGalaxyRedshiftFactorSpline.
 */
NcGalaxyRedshiftFactorSpline *
nc_galaxy_redshift_factor_spline_new (void)
{
  return g_object_new (NC_TYPE_GALAXY_REDSHIFT_FACTOR_SPLINE, NULL);
}

/**
 * nc_galaxy_redshift_factor_spline_ref:
 * @gsdrs: a #NcGalaxyRedshiftFactorSpline
 *
 * Increases the reference count of @gsdrs by one.
 *
 * Returns: (transfer full): @gsdrs.
 */
NcGalaxyRedshiftFactorSpline *
nc_galaxy_redshift_factor_spline_ref (NcGalaxyRedshiftFactorSpline *gsdrs)
{
  return g_object_ref (gsdrs);
}

/**
 * nc_galaxy_redshift_factor_spline_free:
 * @gsdrs: a #NcGalaxyRedshiftFactorSpline
 *
 * Decreases the reference count of @gsdrs by one.
 *
 */
void
nc_galaxy_redshift_factor_spline_free (NcGalaxyRedshiftFactorSpline *gsdrs)
{
  g_object_unref (gsdrs);
}

/**
 * nc_galaxy_redshift_factor_spline_clear:
 * @gsdrs: a #NcGalaxyRedshiftFactorSpline
 *
 * Decreases the reference count of @gsdrs by one, and sets the pointer
 * *@gsdrs to NULL.
 *
 */
void
nc_galaxy_redshift_factor_spline_clear (NcGalaxyRedshiftFactorSpline **gsdrs)
{
  g_clear_object (gsdrs);
}

/**
 * nc_galaxy_redshift_factor_spline_data_set:
 * @gsdrs: a #NcGalaxyRedshiftFactorSpline
 * @data: a #NcGalaxyRedshiftFactorData
 * @pz: a #NcGalaxyPz
 *
 * Sets the redshift density of the galaxy of @data to @pz, which is referenced.
 *
 */
void
nc_galaxy_redshift_factor_spline_data_set (NcGalaxyRedshiftFactorSpline *gsdrs, NcGalaxyRedshiftFactorData *data, NcGalaxyPz *pz)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  _spline_ldata_set_pz (sldata, pz);
}

/**
 * nc_galaxy_redshift_factor_spline_data_peek:
 * @gsdrs: a #NcGalaxyRedshiftFactorSpline
 * @data: a #NcGalaxyRedshiftFactorData
 *
 * Returns: (transfer none): the redshift density of the galaxy of @data, as set.
 */
NcGalaxyPz *
nc_galaxy_redshift_factor_spline_data_peek (NcGalaxyRedshiftFactorSpline *gsdrs, NcGalaxyRedshiftFactorData *data)
{
  SplineLData *sldata = (SplineLData *) data->ldata;

  return sldata->pz;
}

