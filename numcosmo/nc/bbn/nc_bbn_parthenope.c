/***************************************************************************
 *            nc_bbn_parthenope.c
 *
 *  Fri August 29 20:45:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_bbn_parthenope.c
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcBBNParthenope:
 *
 * Primordial Helium from a tabulated PArthENoPE run.
 *
 * Interpolates $Y_p(\omega_b, \Delta N_\mathrm{eff})$ on one of the tables
 * shipped with NumCosmo, selected by #NcBBNParthenope:table. The tables carry
 * $Y_p$ only, so this model implements nc_bbn_Yp_4He() and nc_bbn_get_domain()
 * and nothing else.
 *
 * Outside the tabulated range the model errors rather than extrapolating -- see
 * nc_bbn_check_domain(). All shipped tables cover $\omega_b \in [0.005, 0.04]$
 * and $\Delta N_\mathrm{eff} \in [-3, 7]$.
 *
 * The interpolating splines are deserialized once per table and shared by every
 * instance for the lifetime of the process: a fit or a Monte Carlo run creates
 * many cosmologies, and each one carrying its own copy of the table would mean
 * re-reading it from disk per instance.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/bbn/nc_bbn_parthenope.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/spline/ncm_spline2d.h"
#include "nc_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcBBNParthenope
{
  /*< private >*/
  NcBBN parent_instance;
  NcBBNParthenopeTable table;
  NcmSpline2d *Yp_s2d;
  gdouble cache_wb;
  gdouble cache_DNeff;
  gdouble cache_Yp;
};

G_DEFINE_TYPE (NcBBNParthenope, nc_bbn_parthenope, NC_TYPE_BBN)

enum
{
  PROP_0,
  PROP_TABLE,
  PROP_SIZE,
};

/*
 * One deserialized spline per table, shared by every instance. Held for the
 * lifetime of the process on purpose: the tables are small, immutable and
 * read from disk, and a Monte Carlo run instantiates cosmologies in bulk.
 */
G_LOCK_DEFINE_STATIC (table_cache);

static NcmSpline2d *_table_cache[NC_BBN_PARTHENOPE_TABLE_LEGACY + 1] = { NULL, NULL, NULL };

static NcmSpline2d *
_nc_bbn_parthenope_peek_table (NcBBNParthenopeTable table)
{
  static const gchar *filename[] = {
    "BBN_2017_spline2d.obj",
    "BBN_2017_marcucci_spline2d.obj",
    "BBN_spline2d.obj",
  };
  NcmSpline2d *Yp_s2d;

  g_assert_cmpuint (table, <=, NC_BBN_PARTHENOPE_TABLE_LEGACY);

  G_LOCK (table_cache);

  if (_table_cache[table] == NULL)
  {
    NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
    gchar *full_fname  = ncm_cfg_get_data_filename (filename[table], TRUE);
    GObject *table_obj = ncm_serialize_from_file (ser, full_fname);

    g_assert (NCM_IS_SPLINE2D (table_obj));
    _table_cache[table] = NCM_SPLINE2D (table_obj);

    ncm_serialize_clear (&ser);
    g_free (full_fname);
  }

  Yp_s2d = _table_cache[table];

  G_UNLOCK (table_cache);

  return Yp_s2d;
}

static void
nc_bbn_parthenope_init (NcBBNParthenope *bbn_pn)
{
  bbn_pn->table       = NC_BBN_PARTHENOPE_TABLE_PLANCK2017;
  bbn_pn->Yp_s2d      = NULL;
  bbn_pn->cache_wb    = GSL_NAN;
  bbn_pn->cache_DNeff = GSL_NAN;
  bbn_pn->cache_Yp    = GSL_NAN;
}

static void
_nc_bbn_parthenope_dispose (GObject *object)
{
  NcBBNParthenope *bbn_pn = NC_BBN_PARTHENOPE (object);

  ncm_spline2d_clear (&bbn_pn->Yp_s2d);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_bbn_parthenope_parent_class)->dispose (object);
}

static void
_nc_bbn_parthenope_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_bbn_parthenope_parent_class)->finalize (object);
}

static void
_nc_bbn_parthenope_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcBBNParthenope *bbn_pn = NC_BBN_PARTHENOPE (object);

  g_return_if_fail (NC_IS_BBN_PARTHENOPE (object));

  switch (prop_id)
  {
    case PROP_TABLE:
      nc_bbn_parthenope_set_table (bbn_pn, g_value_get_enum (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_bbn_parthenope_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcBBNParthenope *bbn_pn = NC_BBN_PARTHENOPE (object);

  g_return_if_fail (NC_IS_BBN_PARTHENOPE (object));

  switch (prop_id)
  {
    case PROP_TABLE:
      g_value_set_enum (value, bbn_pn->table);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static gdouble _nc_bbn_parthenope_Yp_4He (NcBBN *bbn, NcHICosmo *cosmo);
static void _nc_bbn_parthenope_get_domain (NcBBN *bbn, gdouble *wb_lb, gdouble *wb_ub, gdouble *DNeff_lb, gdouble *DNeff_ub);

static void
nc_bbn_parthenope_class_init (NcBBNParthenopeClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);
  NcBBNClass *bbn_class      = NC_BBN_CLASS (klass);

  model_class->set_property = &_nc_bbn_parthenope_set_property;
  model_class->get_property = &_nc_bbn_parthenope_get_property;
  object_class->dispose     = &_nc_bbn_parthenope_dispose;
  object_class->finalize    = &_nc_bbn_parthenope_finalize;

  ncm_model_class_set_name_nick (model_class, "Tabulated PArthENoPE primordial Helium.", "NcBBNParthenope");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  /**
   * NcBBNParthenope:table:
   *
   * Which tabulation of $Y_p(\omega_b, \Delta N_\mathrm{eff})$ to interpolate.
   * See #NcBBNParthenopeTable: the tables differ by up to 1% in $Y_p$, so this
   * is a physics choice, not a detail.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_TABLE,
                                   g_param_spec_enum ("table",
                                                      NULL,
                                                      "PArthENoPE tabulation of Yp",
                                                      NC_TYPE_BBN_PARTHENOPE_TABLE,
                                                      NC_BBN_PARTHENOPE_TABLE_PLANCK2017,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_model_class_check_params_info (model_class);

  ncm_model_class_add_impl_opts (model_class, NC_BBN_IMPL_Yp_4He, NC_BBN_IMPL_get_domain, -1);

  bbn_class->Yp_4He     = &_nc_bbn_parthenope_Yp_4He;
  bbn_class->get_domain = &_nc_bbn_parthenope_get_domain;
}

static void
_nc_bbn_parthenope_get_domain (NcBBN *bbn, gdouble *wb_lb, gdouble *wb_ub, gdouble *DNeff_lb, gdouble *DNeff_ub)
{
  NcBBNParthenope *bbn_pn = NC_BBN_PARTHENOPE (bbn);
  NcmVector *wb_v         = ncm_spline2d_peek_xv (bbn_pn->Yp_s2d);
  NcmVector *DNeff_v      = ncm_spline2d_peek_yv (bbn_pn->Yp_s2d);

  wb_lb[0]    = ncm_vector_get (wb_v, 0);
  wb_ub[0]    = ncm_vector_get (wb_v, ncm_vector_len (wb_v) - 1);
  DNeff_lb[0] = ncm_vector_get (DNeff_v, 0);
  DNeff_ub[0] = ncm_vector_get (DNeff_v, ncm_vector_len (DNeff_v) - 1);
}

static gdouble
_nc_bbn_parthenope_Yp_4He (NcBBN *bbn, NcHICosmo *cosmo)
{
  NcBBNParthenope *bbn_pn = NC_BBN_PARTHENOPE (bbn);
  const gdouble wb        = nc_hicosmo_Omega_b0h2 (cosmo);
  const gdouble DNeff     = nc_hicosmo_Neff (cosmo) - 3.046;

  /*
   * Keyed on the two values the table is a function of, not on the cosmology's
   * pkey: a parameter change that moves neither leaves the answer alone.
   */
  if ((wb != bbn_pn->cache_wb) || (DNeff != bbn_pn->cache_DNeff))
  {
    nc_bbn_check_domain (bbn, wb, DNeff);

    bbn_pn->cache_Yp    = ncm_spline2d_eval (bbn_pn->Yp_s2d, wb, DNeff);
    bbn_pn->cache_wb    = wb;
    bbn_pn->cache_DNeff = DNeff;
  }

  return bbn_pn->cache_Yp;
}

/**
 * nc_bbn_parthenope_new:
 *
 * Creates a new #NcBBNParthenope on the default table,
 * %NC_BBN_PARTHENOPE_TABLE_PLANCK2017.
 *
 * Returns: (transfer full): a new #NcBBNParthenope.
 */
NcBBNParthenope *
nc_bbn_parthenope_new (void)
{
  return g_object_new (NC_TYPE_BBN_PARTHENOPE, NULL);
}

/**
 * nc_bbn_parthenope_new_from_table:
 * @table: a #NcBBNParthenopeTable
 *
 * Creates a new #NcBBNParthenope interpolating @table.
 *
 * Returns: (transfer full): a new #NcBBNParthenope.
 */
NcBBNParthenope *
nc_bbn_parthenope_new_from_table (NcBBNParthenopeTable table)
{
  return g_object_new (NC_TYPE_BBN_PARTHENOPE,
                       "table", table,
                       NULL);
}

/**
 * nc_bbn_parthenope_ref:
 * @bbn_pn: a #NcBBNParthenope
 *
 * Increases the reference count of @bbn_pn by one.
 *
 * Returns: (transfer full): @bbn_pn.
 */
NcBBNParthenope *
nc_bbn_parthenope_ref (NcBBNParthenope *bbn_pn)
{
  return g_object_ref (bbn_pn);
}

/**
 * nc_bbn_parthenope_free:
 * @bbn_pn: a #NcBBNParthenope
 *
 * Decreases the reference count of @bbn_pn by one.
 *
 */
void
nc_bbn_parthenope_free (NcBBNParthenope *bbn_pn)
{
  g_object_unref (bbn_pn);
}

/**
 * nc_bbn_parthenope_clear:
 * @bbn_pn: a #NcBBNParthenope
 *
 * If *@bbn_pn is different from NULL, decreases the reference count of
 * *@bbn_pn by one and sets *@bbn_pn to NULL.
 *
 */
void
nc_bbn_parthenope_clear (NcBBNParthenope **bbn_pn)
{
  g_clear_object (bbn_pn);
}

/**
 * nc_bbn_parthenope_set_table:
 * @bbn_pn: a #NcBBNParthenope
 * @table: a #NcBBNParthenopeTable
 *
 * Sets which tabulation of $Y_p$ @bbn_pn interpolates.
 *
 */
void
nc_bbn_parthenope_set_table (NcBBNParthenope *bbn_pn, NcBBNParthenopeTable table)
{
  ncm_spline2d_clear (&bbn_pn->Yp_s2d);

  bbn_pn->table  = table;
  bbn_pn->Yp_s2d = ncm_spline2d_ref (_nc_bbn_parthenope_peek_table (table));

  bbn_pn->cache_wb    = GSL_NAN;
  bbn_pn->cache_DNeff = GSL_NAN;
  bbn_pn->cache_Yp    = GSL_NAN;
}

/**
 * nc_bbn_parthenope_get_table:
 * @bbn_pn: a #NcBBNParthenope
 *
 * Gets which tabulation of $Y_p$ @bbn_pn interpolates.
 *
 * Returns: the #NcBBNParthenopeTable in use.
 */
NcBBNParthenopeTable
nc_bbn_parthenope_get_table (NcBBNParthenope *bbn_pn)
{
  return bbn_pn->table;
}

