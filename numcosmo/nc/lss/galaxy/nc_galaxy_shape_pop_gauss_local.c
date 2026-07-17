/***************************************************************************
 *            nc_galaxy_shape_pop_gauss_local.c
 *
 *  Thu Jul 2 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_pop_gauss_local.c
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
 * NcGalaxyShapePopGaussLocal:
 *
 * Truncated-Gaussian intrinsic ellipticity distribution, per-galaxy width.
 *
 * Same truncated-Gaussian family as #NcGalaxyShapePopGauss (the "Global"
 * variant) — a SIBLING, not a subclass: both derive directly from
 * #NcGalaxyShapePop and share the family's mechanism (data layout,
 * eval_p, gen) as plain functions reused across the two files (see
 * nc_galaxy_shape_pop_gauss_private.h), not through inheritance. The only
 * real difference is where the width $\sigma$ comes from: here it is not a
 * shared model parameter but is resolved per galaxy from an input RMS
 * ellipticity @e_rms carried on #NcGalaxyShapePopData (read from the catalog
 * column "e_rms"). Since $\sigma$ parameterizes an *untruncated* Gaussian
 * while @e_rms is the truncated-disk RMS, resolving $\sigma$ from a target
 * @e_rms requires inverting the (monotonic) forward map; this is done once
 * per galaxy in prepare() via bisection, never in the marginal-likelihood hot
 * loop.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_pop_gauss_local.h"
#include "nc/lss/galaxy/nc_galaxy_shape_pop_gauss_private.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyShapePopGaussLocal
{
  NcGalaxyShapePop parent_instance;
};

enum
{
  PROP_0,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxyShapePopGaussLocal, nc_galaxy_shape_pop_gauss_local, NC_TYPE_GALAXY_SHAPE_POP);

static void
nc_galaxy_shape_pop_gauss_local_init (NcGalaxyShapePopGaussLocal *gspgl)
{
}

static void
_nc_galaxy_shape_pop_gauss_local_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_shape_pop_gauss_local_parent_class)->finalize (object);
}

static void _nc_galaxy_shape_pop_gauss_local_data_init (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
static void _nc_galaxy_shape_pop_gauss_local_prepare (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);
static gdouble _nc_galaxy_shape_pop_gauss_local_e_rms (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data);

static void
nc_galaxy_shape_pop_gauss_local_class_init (NcGalaxyShapePopGaussLocalClass *klass)
{
  NcGalaxyShapePopClass *gsp_class = NC_GALAXY_SHAPE_POP_CLASS (klass);
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

  object_class->finalize = &_nc_galaxy_shape_pop_gauss_local_finalize;

  ncm_model_class_set_name_nick (model_class, "Truncated Gaussian intrinsic ellipticity distribution (per-galaxy)", "GaussLocalIntrinsic");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);
  ncm_model_class_check_params_info (model_class);

  /* eval_p, gen and eval_p_rho2_g_series are the SAME functions
   * NcGalaxyShapePopGauss uses (reused directly, not inherited): all three
   * only ever touch data->ldata, whose {norm, inv_2sigma2, sigma} layout
   * and meaning are identical here. */
  gsp_class->data_init            = &_nc_galaxy_shape_pop_gauss_local_data_init;
  gsp_class->prepare              = &_nc_galaxy_shape_pop_gauss_local_prepare;
  gsp_class->eval_p               = &_nc_galaxy_shape_pop_gauss_eval_p;
  gsp_class->gen                  = &_nc_galaxy_shape_pop_gauss_gen;
  gsp_class->e_rms                = &_nc_galaxy_shape_pop_gauss_local_e_rms;
  gsp_class->eval_p_rho2_g_series = &_nc_galaxy_shape_pop_gauss_eval_p_rho2_g_series;
}

static void
_nc_galaxy_shape_pop_gauss_local_ldata_read_row (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i)
{
  data->e_rms = nc_galaxy_wl_obs_get (obs, NC_GALAXY_SHAPE_POP_GAUSS_LOCAL_COL_E_RMS, i, NULL);
}

static void
_nc_galaxy_shape_pop_gauss_local_ldata_write_row (NcGalaxyShapePopData *data, NcGalaxyWLObs *obs, const guint i)
{
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SHAPE_POP_GAUSS_LOCAL_COL_E_RMS, i, data->e_rms, NULL);
}

static void
_nc_galaxy_shape_pop_gauss_local_ldata_required_columns (NcGalaxyShapePopData *data, GList **columns)
{
  *columns = g_list_append (*columns, g_strdup (NC_GALAXY_SHAPE_POP_GAUSS_LOCAL_COL_E_RMS));
}

static void
_nc_galaxy_shape_pop_gauss_local_data_init (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  /* Plain function call, not vtable chaining: NcGalaxyShapePopGaussLocal is
   * not a subclass of NcGalaxyShapePopGauss, just a sibling reusing its
   * ldata allocation (same layout, see nc_galaxy_shape_pop_gauss_private.h). */
  _nc_galaxy_shape_pop_gauss_data_init (gsp, data);

  /* ldata_get_sigma is already set to _nc_galaxy_shape_pop_gauss_ldata_get_sigma
   * by the call above (same ldata layout), no need to reassign it. */

  data->ldata_read_row         = &_nc_galaxy_shape_pop_gauss_local_ldata_read_row;
  data->ldata_write_row        = &_nc_galaxy_shape_pop_gauss_local_ldata_write_row;
  data->ldata_required_columns = &_nc_galaxy_shape_pop_gauss_local_ldata_required_columns;
}

#define NC_GALAXY_SHAPE_POP_GAUSS_LOCAL_E_RMS_MAX (0.5 - 1.0e-9)

static gdouble
_e_rms_of_sigma (const gdouble sigma)
{
  const gdouble lambda = 0.5 / (sigma * sigma);
  const gdouble exp_ml = exp (-lambda);
  const gdouble mean_x = (1.0 - exp_ml * (1.0 + lambda)) / (lambda * (1.0 - exp_ml));

  return sqrt (0.5 * mean_x);
}

/*
 * e_rms(sigma) is monotonically increasing, from 0 (sigma -> 0) to 0.5
 * (sigma -> infinity, the untruncated-uniform-x limit): bisection is simple
 * and robust here, and this only ever runs once per galaxy in prepare().
 */
static gdouble
_sigma_from_e_rms (const gdouble e_rms)
{
  gdouble sigma_lo = 1.0e-4;
  gdouble sigma_hi = 1.0e2;
  guint iter;

  g_assert_cmpfloat (e_rms, >, 0.0);
  g_assert_cmpfloat (e_rms, <, NC_GALAXY_SHAPE_POP_GAUSS_LOCAL_E_RMS_MAX);

  for (iter = 0; iter < 100; iter++)
  {
    const gdouble sigma_mid = 0.5 * (sigma_lo + sigma_hi);
    const gdouble e_mid     = _e_rms_of_sigma (sigma_mid);

    if (e_mid < e_rms)
      sigma_lo = sigma_mid;
    else
      sigma_hi = sigma_mid;

    if ((sigma_hi - sigma_lo) < 1.0e-14 * sigma_hi)
      break;
  }

  return 0.5 * (sigma_lo + sigma_hi);
}

static void
_nc_galaxy_shape_pop_gauss_local_prepare (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  const gdouble sigma = _sigma_from_e_rms (data->e_rms);

  _nc_galaxy_shape_pop_gauss_ldata_set_sigma (data, sigma);
}

static gdouble
_nc_galaxy_shape_pop_gauss_local_e_rms (NcGalaxyShapePop *gsp, NcGalaxyShapePopData *data)
{
  return data->e_rms;
}

/**
 * nc_galaxy_shape_pop_gauss_local_new:
 *
 * Creates a new #NcGalaxyShapePopGaussLocal.
 *
 * Returns: (transfer full): a new #NcGalaxyShapePopGaussLocal.
 */
NcGalaxyShapePopGaussLocal *
nc_galaxy_shape_pop_gauss_local_new (void)
{
  NcGalaxyShapePopGaussLocal *gspgl = g_object_new (NC_TYPE_GALAXY_SHAPE_POP_GAUSS_LOCAL,
                                                    NULL);

  return gspgl;
}

/**
 * nc_galaxy_shape_pop_gauss_local_ref:
 * @gspgl: a #NcGalaxyShapePopGaussLocal
 *
 * Increases the reference count of @gspgl by one.
 *
 * Returns: (transfer full): @gspgl.
 */
NcGalaxyShapePopGaussLocal *
nc_galaxy_shape_pop_gauss_local_ref (NcGalaxyShapePopGaussLocal *gspgl)
{
  return g_object_ref (gspgl);
}

/**
 * nc_galaxy_shape_pop_gauss_local_free:
 * @gspgl: a #NcGalaxyShapePopGaussLocal
 *
 * Decreases the reference count of @gspgl by one.
 *
 */
void
nc_galaxy_shape_pop_gauss_local_free (NcGalaxyShapePopGaussLocal *gspgl)
{
  g_object_unref (gspgl);
}

/**
 * nc_galaxy_shape_pop_gauss_local_clear:
 * @gspgl: a #NcGalaxyShapePopGaussLocal
 *
 * Decreases the reference count of *@gspgl by one, and sets the pointer
 * *@gspgl to NULL.
 *
 */
void
nc_galaxy_shape_pop_gauss_local_clear (NcGalaxyShapePopGaussLocal **gspgl)
{
  g_clear_object (gspgl);
}

/**
 * nc_galaxy_shape_pop_gauss_local_data_set:
 * @gspgl: a #NcGalaxyShapePopGaussLocal
 * @data: a #NcGalaxyShapePopData
 * @e_rms: the per-galaxy RMS intrinsic ellipticity
 *
 * Sets the per-galaxy RMS intrinsic ellipticity input.
 *
 */
void
nc_galaxy_shape_pop_gauss_local_data_set (NcGalaxyShapePopGaussLocal *gspgl, NcGalaxyShapePopData *data, const gdouble e_rms)
{
  data->e_rms = e_rms;
}

/**
 * nc_galaxy_shape_pop_gauss_local_data_get:
 * @gspgl: a #NcGalaxyShapePopGaussLocal
 * @data: a #NcGalaxyShapePopData
 *
 * Returns: the per-galaxy RMS intrinsic ellipticity input.
 */
gdouble
nc_galaxy_shape_pop_gauss_local_data_get (NcGalaxyShapePopGaussLocal *gspgl, NcGalaxyShapePopData *data)
{
  return data->e_rms;
}

