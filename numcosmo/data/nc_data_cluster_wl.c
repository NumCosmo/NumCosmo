/***************************************************************************
 *            nc_data_cluster_wl.c
 *
 *  Mon Jul 27 16:10:25 2020
 *  Copyright  2020  Mariana Penna Lima
 *  <pennalima@gmail.com>
 *  Tue Jun 15 16:00:13 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_cluster_wl.c
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcDataClusterWL:
 *
 * Cluster weak lensing likelihood.
 *
 * This class implements the weak lensing likelihood for galaxy clusters.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cluster_wl.h"

#include "galaxy/nc_galaxy_sd_shape.h"
#include "galaxy/nc_galaxy_sd_obs_redshift.h"
#include "galaxy/nc_galaxy_sd_obs_redshift_spec.h"
#include "galaxy/nc_galaxy_sd_position.h"
#include "math/ncm_integral_nd.h"
#include <math.h>
#include <gsl/gsl_math.h>

#include "nc_enum_types.h"
#include "nc_hicosmo.h"
#include "lss/nc_halo_position.h"
#include "lss/nc_halo_density_profile.h"
#include "lss/nc_wl_surface_mass_density.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_integrate.h"
#include "math/ncm_memory_pool.h"

#define NC_GALAXY_LOW_PROB 1.0e6
struct _NcDataClusterWLPrivate
{
  NcGalaxyWLObs *obs;
  GPtrArray *shape_data;
  gboolean constructed;
  gdouble r_min;
  gdouble r_max;
  gdouble dr;
  gdouble prec;
  guint len;
  NcmModelCtrl *ctrl_redshift;
  NcmModelCtrl *ctrl_position;
  NcmModelCtrl *ctrl_shape;
  NcDataClusterWLResampleFlag resample_flag;
  gboolean enable_parallel;
  /* Integration temporary variables */
  NcmVector *err;
  NcmVector *zpi;
  NcmVector *zpf;
  NcmVector *res;
  NcmVector *tmp;
  NcHICosmo *cosmo;
  NcWLSurfaceMassDensity *surface_mass_density;
  NcHaloDensityProfile *density_profile;
  NcHaloPosition *halo_position;
  NcGalaxySDShape *galaxy_shape;
  NcGalaxySDObsRedshift *galaxy_redshift;
  NcGalaxySDPosition *galaxy_position;
};

enum
{
  PROP_0,
  PROP_OBS,
  PROP_R_MIN,
  PROP_R_MAX,
  PROP_PREC,
  PROP_LEN,
  PROP_SIZE,
  PROP_RESAMPLE_FLAG,
  PROP_ENABLE_PARALLEL,
};

struct _NcDataClusterWL
{
  /*< private >*/
  NcmData parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcDataClusterWL, nc_data_cluster_wl, NCM_TYPE_DATA);

static void
nc_data_cluster_wl_init (NcDataClusterWL *dcwl)
{
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  self->obs             = NULL;
  self->shape_data      = g_ptr_array_new ();
  self->constructed     = FALSE;
  self->r_max           = 0.0;
  self->r_min           = 0.0;
  self->dr              = 0.0;
  self->prec            = 1.0e-6;
  self->len             = 0;
  self->ctrl_redshift   = ncm_model_ctrl_new (NULL);
  self->ctrl_position   = ncm_model_ctrl_new (NULL);
  self->ctrl_shape      = ncm_model_ctrl_new (NULL);
  self->resample_flag   = NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_ALL;
  self->enable_parallel = FALSE;

  self->err = ncm_vector_new (1);
  self->zpi = ncm_vector_new (1);
  self->zpf = ncm_vector_new (1);
  self->res = ncm_vector_new (1);
  self->tmp = ncm_vector_new (1);

  g_ptr_array_set_free_func (self->shape_data, (GDestroyNotify) nc_galaxy_sd_shape_data_unref);
}

static void
nc_data_cluster_wl_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (object);
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  g_return_if_fail (NC_IS_DATA_CLUSTER_WL (object));

  switch (prop_id)
  {
    case PROP_OBS:
      nc_data_cluster_wl_set_obs (dcwl, g_value_get_object (value));
      break;
    case PROP_R_MIN:
      self->r_min = g_value_get_double (value);

      if (self->constructed)
      {
        g_assert_cmpfloat (self->r_min, <, self->r_max);
        self->dr = self->r_max - self->r_min;
      }

      break;
    case PROP_R_MAX:
      self->r_max = g_value_get_double (value);

      if (self->constructed)
      {
        g_assert_cmpfloat (self->r_min, <, self->r_max);
        self->dr = self->r_max - self->r_min;
      }

      break;
    case PROP_PREC:
      nc_data_cluster_wl_set_prec (dcwl, g_value_get_double (value));
      break;
    case PROP_LEN:
      self->len = g_value_get_uint (value);
      break;
    case PROP_RESAMPLE_FLAG:
      nc_data_cluster_wl_set_resample_flag (dcwl, g_value_get_flags (value));
      break;
    case PROP_ENABLE_PARALLEL:
      self->enable_parallel = g_value_get_boolean (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_cluster_wl_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (object);
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  g_return_if_fail (NC_IS_DATA_CLUSTER_WL (object));

  switch (prop_id)
  {
    case PROP_OBS:
      g_value_set_object (value, nc_data_cluster_wl_peek_obs (dcwl));
      break;
    case PROP_R_MIN:
      g_value_set_double (value, self->r_min);
      break;
    case PROP_R_MAX:
      g_value_set_double (value, self->r_max);
      break;
    case PROP_PREC:
      g_value_set_double (value, self->prec);
      break;
    case PROP_LEN:
      g_value_set_uint (value, self->len);
      break;
    case PROP_RESAMPLE_FLAG:
      g_value_set_flags (value, nc_data_cluster_wl_get_resample_flag (dcwl));
      break;
    case PROP_ENABLE_PARALLEL:
      g_value_set_boolean (value, self->enable_parallel);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_cluster_wl_dispose (GObject *object)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (object);
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  nc_galaxy_wl_obs_clear (&self->obs);

  ncm_model_ctrl_clear (&self->ctrl_redshift);
  ncm_model_ctrl_clear (&self->ctrl_position);
  ncm_model_ctrl_clear (&self->ctrl_shape);

  g_clear_pointer (&self->shape_data, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_wl_parent_class)->dispose (object);
}

static void
nc_data_cluster_wl_finalize (GObject *object)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (object);
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  ncm_vector_clear (&self->err);
  ncm_vector_clear (&self->zpi);
  ncm_vector_clear (&self->zpf);
  ncm_vector_clear (&self->res);
  ncm_vector_clear (&self->tmp);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_wl_parent_class)->finalize (object);
}

static void
_nc_data_cluster_wl_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_cluster_wl_parent_class)->constructed (object);
  {
    NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (object);
    NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

    g_assert_cmpfloat (self->r_min, <, self->r_max);
    self->dr = self->r_max - self->r_min;

    self->constructed = TRUE;
  }
}

static void _nc_data_cluster_wl_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);
static void _nc_data_cluster_wl_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static guint _nc_data_cluster_wl_get_len (NcmData *data);
static void _nc_data_cluster_wl_prepare (NcmData *data, NcmMSet *mset);

static void
nc_data_cluster_wl_class_init (NcDataClusterWLClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->set_property = nc_data_cluster_wl_set_property;
  object_class->get_property = nc_data_cluster_wl_get_property;
  object_class->dispose      = nc_data_cluster_wl_dispose;
  object_class->finalize     = nc_data_cluster_wl_finalize;
  object_class->constructed  = _nc_data_cluster_wl_constructed;

  /**
   * NcDataClusterWL:obs:
   *
   * Galaxy weak lensing observables.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_OBS,
                                   g_param_spec_object ("obs",
                                                        NULL,
                                                        "Galaxy weak lensing observables",
                                                        NC_TYPE_GALAXY_WL_OBS,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:r-min:
   *
   * Minimum radius of the weak lensing observables.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_R_MIN,
                                   g_param_spec_double ("r-min",
                                                        NULL,
                                                        "Minimum radius of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:r-max:
   *
   * Maximum radius of the weak lensing observables.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_R_MAX,
                                   g_param_spec_double ("r-max",
                                                        NULL,
                                                        "Maximum radius of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 5.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:prec:
   *
   * Precision for integral.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_PREC,
                                   g_param_spec_double ("prec",
                                                        NULL,
                                                        "Precision for integral",
                                                        0.0, G_MAXDOUBLE, 1.0e-6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:len:
   *
   * Number of galaxies.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LEN,
                                   g_param_spec_uint ("len",
                                                      NULL,
                                                      "Number of galaxies",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:resample-flag:
   *
   * Resample flag.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RESAMPLE_FLAG,
                                   g_param_spec_flags ("resample-flag",
                                                       NULL,
                                                       "Resample flag",
                                                       NC_TYPE_DATA_CLUSTER_WL_RESAMPLE_FLAG,
                                                       NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_ALL,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:enable-parallel
   *
   * Enable parallelization.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ENABLE_PARALLEL,
                                   g_param_spec_boolean ("enable-parallel",
                                                         NULL,
                                                         "Enable parallelization",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->bootstrap  = TRUE;
  data_class->resample   = &_nc_data_cluster_wl_resample;
  data_class->m2lnL_val  = &_nc_data_cluster_wl_m2lnL_val;
  data_class->get_length = &_nc_data_cluster_wl_get_len;
  data_class->prepare    = &_nc_data_cluster_wl_prepare;
}

struct _NcDataClusterWLIntArg
{
  NcGalaxySDObsRedshiftIntegrand *integrand_redshift;
  NcGalaxySDPositionIntegrand *integrand_position;
  NcGalaxySDShapeIntegrand *integrand_shape;
  NcGalaxySDShapeData *data;
  guint gal_i;
};

static void nc_data_cluster_wl_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
static void nc_data_cluster_wl_int_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);

NCM_INTEGRAL_ND_DEFINE_TYPE (NC, DATA_CLUSTER_WL_INT, NcDataClusterWLInt, nc_data_cluster_wl_integ, nc_data_cluster_wl_int_dim, nc_data_cluster_wl_integ, struct _NcDataClusterWLIntArg);

static void
nc_data_cluster_wl_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcDataClusterWLInt *lh_int                         = NC_DATA_CLUSTER_WL_INT (intnd);
  NcGalaxySDShapeData *data                          = lh_int->data.data;
  NcGalaxySDObsRedshiftIntegrand *integrand_redshift = lh_int->data.integrand_redshift;
  NcGalaxySDPositionIntegrand *integrand_position    = lh_int->data.integrand_position;
  NcGalaxySDShapeIntegrand *integrand_shape          = lh_int->data.integrand_shape;
  guint i;

  for (i = 0; i < npoints; i++)
  {
    const gdouble z         = ncm_vector_fast_get (x, i);
    const gdouble int_z     = nc_galaxy_sd_obs_redshift_integrand_eval (integrand_redshift, z, data->sdpos_data->sdz_data);
    const gdouble int_pos   = nc_galaxy_sd_position_integrand_eval (integrand_position, data->sdpos_data);
    const gdouble int_shape = nc_galaxy_sd_shape_integrand_eval (integrand_shape, z, data);
    const gdouble res       = int_z * int_pos * int_shape;

    ncm_vector_set (fval, i, res);
  }
}

static void
nc_data_cluster_wl_int_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  *dim  = 1;
  *fdim = 1;
}

static gdouble
_nc_data_cluster_wl_eval_m2lnP_weight (NcDataClusterWL *dcwl, const gdouble m2lnP, const gdouble r)
{
  /* NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl); */

/* Removed previous implementation: it was incomplete and introduced bias in the
 * estimates. Leaving this dummy function in place, as it's still unclear whether the
 * radius cuts should be handled at this stage.
 */

  return m2lnP;
}

static gdouble
_nc_data_cluster_wl_eval_m2lnP_integ (NcDataClusterWL *dcwl, NcmMSet *mset, NcmVector *m2lnP_gal)
{
  NcmData *data                       = NCM_DATA (dcwl);
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);
  gdouble result                      = 0;

  #pragma omp parallel reduction(+:result) if (self->enable_parallel)
  {
    NcDataClusterWLInt *likelihood_integral            = g_object_new (nc_data_cluster_wl_integ_get_type (), NULL);
    NcmIntegralND *lh_int                              = NCM_INTEGRAL_ND (likelihood_integral);
    NcGalaxySDObsRedshiftIntegrand *integrand_redshift = nc_galaxy_sd_obs_redshift_integ (self->galaxy_redshift);
    NcGalaxySDPositionIntegrand *integrand_position    = nc_galaxy_sd_position_integ (self->galaxy_position);
    NcGalaxySDShapeIntegrand *integrand_shape          = nc_galaxy_sd_shape_integ (self->galaxy_shape);
    NcmVector *zpi_v                                   = ncm_vector_new (1);
    NcmVector *zpf_v                                   = ncm_vector_new (1);
    NcmVector *res_v                                   = ncm_vector_new (1);
    NcmVector *err_v                                   = ncm_vector_new (1);

    ncm_integral_nd_set_reltol (lh_int, self->prec);
    ncm_integral_nd_set_abstol (lh_int, 0.0);
    ncm_integral_nd_set_method (lh_int, NCM_INTEGRAL_ND_METHOD_CUBATURE_H_V);

    likelihood_integral->data.integrand_redshift = integrand_redshift;
    likelihood_integral->data.integrand_position = integrand_position;
    likelihood_integral->data.integrand_shape    = integrand_shape;

    nc_galaxy_sd_obs_redshift_integrand_prepare (integrand_redshift, mset);
    nc_galaxy_sd_position_integrand_prepare (integrand_position, mset);
    nc_galaxy_sd_shape_integrand_prepare (integrand_shape, mset);

    if (m2lnP_gal != NULL)
      g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

    if (!ncm_data_bootstrap_enabled (data))
    {
      guint gal_i;

      #pragma omp for schedule (dynamic, 5)

      for (gal_i = 0; gal_i < self->len; gal_i++)
      {
        NcGalaxySDShapeData *s_data       = NC_GALAXY_SD_SHAPE_DATA (ncm_obj_array_peek (self->shape_data, gal_i));
        NcGalaxySDPositionData *p_data    = s_data->sdpos_data;
        NcGalaxySDObsRedshiftData *z_data = p_data->sdz_data;
        gdouble m2lnP_gal_i, P_gal_i;
        gdouble zpi;
        gdouble zpf;

        likelihood_integral->data.gal_i = gal_i;

        nc_galaxy_sd_obs_redshift_get_lim (self->galaxy_redshift, z_data, &zpi, &zpf);
        ncm_vector_fast_set (zpi_v, 0, zpi);
        ncm_vector_fast_set (zpf_v, 0, zpf);

        likelihood_integral->data.data = s_data;

        ncm_integral_nd_eval (lh_int, zpi_v, zpf_v, res_v, err_v);

        P_gal_i     = ncm_vector_fast_get (err_v, 0);
        m2lnP_gal_i = P_gal_i > 0.0 ? -2.0 * log (P_gal_i) : NC_GALAXY_LOW_PROB;

        if (!gsl_finite (m2lnP_gal_i))
        {
          g_warning ("_nc_data_cluster_wl_eval_m2lnP_integ: galaxy %d has undefined likelihood [%g]. Skipping it.", gal_i, m2lnP_gal_i);
          continue;
        }

        if (m2lnP_gal != NULL)
          ncm_vector_set (m2lnP_gal, gal_i, m2lnP_gal_i);

        result += _nc_data_cluster_wl_eval_m2lnP_weight (dcwl, m2lnP_gal_i, 0.0);
      }
    }
    else
    {
      NcmBootstrap *bstrap = ncm_data_peek_bootstrap (data);

      const guint bsize = ncm_bootstrap_get_bsize (bstrap);
      guint i;

      #pragma omp for schedule (dynamic, 5)

      for (i = 0; i < bsize; i++)
      {
        guint gal_i = ncm_bootstrap_get (bstrap, i);

        NcGalaxySDShapeData *s_data_i     = NC_GALAXY_SD_SHAPE_DATA (ncm_obj_array_peek (self->shape_data, gal_i));
        NcGalaxySDPositionData *p_data    = s_data_i->sdpos_data;
        NcGalaxySDObsRedshiftData *z_data = p_data->sdz_data;

        gdouble m2lnP_gal_i, P_gal_i;
        gdouble zpi;
        gdouble zpf;

        likelihood_integral->data.gal_i = gal_i;

        nc_galaxy_sd_obs_redshift_get_lim (self->galaxy_redshift, z_data, &zpi, &zpf);
        ncm_vector_fast_set (zpi_v, 0, zpi);
        ncm_vector_fast_set (zpf_v, 0, zpf);

        likelihood_integral->data.data = s_data_i;

        ncm_integral_nd_eval (lh_int, zpi_v, zpf_v, res_v, err_v);

        P_gal_i     = ncm_vector_fast_get (res_v, 0);
        m2lnP_gal_i = P_gal_i > 0.0 ? -2.0 * log (P_gal_i) : NC_GALAXY_LOW_PROB;

        if (!gsl_finite (m2lnP_gal_i))
        {
          g_warning ("_nc_data_cluster_wl_eval_m2lnP_integ: galaxy %d has undefined likelihood [%g]. Skipping it.", gal_i, m2lnP_gal_i);
          continue;
        }

        if (m2lnP_gal != NULL)
          ncm_vector_set (m2lnP_gal, gal_i, m2lnP_gal_i);

        result += _nc_data_cluster_wl_eval_m2lnP_weight (dcwl, m2lnP_gal_i, 0.0);
      }
    }

    ncm_vector_free (zpi_v);
    ncm_vector_free (zpf_v);
    ncm_vector_free (res_v);
    ncm_vector_free (err_v);
    ncm_integral_nd_clear (&lh_int);
    nc_galaxy_sd_shape_integrand_free (integrand_shape);
    nc_galaxy_sd_position_integrand_free (integrand_position);
    nc_galaxy_sd_obs_redshift_integrand_free (integrand_redshift);
  }

  return result;
}

static void
_nc_data_cluster_wl_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcDataClusterWL *dcwl                  = NC_DATA_CLUSTER_WL (data);
  NcDataClusterWLPrivate * const self    = nc_data_cluster_wl_get_instance_private (dcwl);
  NcGalaxySDShape *galaxy_shape          = NC_GALAXY_SD_SHAPE (ncm_mset_peek (mset, nc_galaxy_sd_shape_id ()));
  NcGalaxySDObsRedshift *galaxy_redshift = NC_GALAXY_SD_OBS_REDSHIFT (ncm_mset_peek (mset, nc_galaxy_sd_obs_redshift_id ()));
  NcGalaxySDPosition *galaxy_position    = NC_GALAXY_SD_POSITION (ncm_mset_peek (mset, nc_galaxy_sd_position_id ()));
  NcHaloPosition *halo_position          = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcHICosmo *cosmo                       = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));

  guint gal_i;

  nc_halo_position_prepare_if_needed (halo_position, cosmo);

  for (gal_i = 0; gal_i < self->len; gal_i++)
  {
    NcGalaxySDShapeData *data_i       = NC_GALAXY_SD_SHAPE_DATA (ncm_obj_array_peek (self->shape_data, gal_i));
    NcGalaxySDPositionData *p_data    = data_i->sdpos_data;
    NcGalaxySDObsRedshiftData *z_data = p_data->sdz_data;

    if (self->resample_flag & NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_REDSHIFT)
    {
      nc_galaxy_sd_obs_redshift_prepare (galaxy_redshift, z_data);
      nc_galaxy_sd_obs_redshift_gen (galaxy_redshift, z_data, rng);
    }

    if (self->resample_flag & NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_POSITION)
    {
      gdouble radius;

      do {
        nc_galaxy_sd_position_gen (galaxy_position, p_data, rng);

        radius = nc_halo_position_projected_radius_from_ra_dec (halo_position, cosmo, p_data->ra, p_data->dec);
      } while (radius < self->r_min || radius > self->r_max);
    }

    nc_galaxy_sd_shape_gen (galaxy_shape, mset, data_i, rng);
    nc_galaxy_sd_shape_data_write_row (data_i, self->obs, gal_i);
  }
}

static gdouble
_nc_data_cluster_wl_eval_m2lnP (NcDataClusterWL *dcwl, NcmMSet *mset, NcmVector *m2lnP_gal)
{
  NcmData *data                       = NCM_DATA (dcwl);
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  gdouble result = 0;

  #pragma omp parallel reduction(+:result) if (self->enable_parallel)
  {
    NcGalaxySDObsRedshiftIntegrand *integrand_redshift = nc_galaxy_sd_obs_redshift_integ (self->galaxy_redshift);
    NcGalaxySDPositionIntegrand *integrand_position    = nc_galaxy_sd_position_integ (self->galaxy_position);
    NcGalaxySDShapeIntegrand *integrand_shape          = nc_galaxy_sd_shape_integ (self->galaxy_shape);

    if (m2lnP_gal != NULL)
      g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

    nc_galaxy_sd_obs_redshift_integrand_prepare (integrand_redshift, mset);
    nc_galaxy_sd_position_integrand_prepare (integrand_position, mset);
    nc_galaxy_sd_shape_integrand_prepare (integrand_shape, mset);

    if (!ncm_data_bootstrap_enabled (data))
    {
      guint gal_i;

      #pragma omp for schedule (dynamic, 5)

      for (gal_i = 0; gal_i < self->len; gal_i++)
      {
        NcGalaxySDShapeData *s_data_i = NC_GALAXY_SD_SHAPE_DATA (ncm_obj_array_peek (self->shape_data, gal_i));

        const gdouble z           = s_data_i->sdpos_data->sdz_data->z;
        const gdouble int_z       = nc_galaxy_sd_obs_redshift_integrand_eval (integrand_redshift, z, s_data_i->sdpos_data->sdz_data);
        const gdouble int_pos     = nc_galaxy_sd_position_integrand_eval (integrand_position, s_data_i->sdpos_data);
        const gdouble int_shape   = nc_galaxy_sd_shape_integrand_eval (integrand_shape, z, s_data_i);
        const gdouble P_gal_i     = int_z * int_pos * int_shape;
        const gdouble m2lnP_gal_i = P_gal_i > 0.0 ? -2.0 * log (P_gal_i) : NC_GALAXY_LOW_PROB;

        if (!gsl_finite (m2lnP_gal_i))
        {
          g_warning ("_nc_data_cluster_wl_eval_m2lnP: galaxy %d has undefined likelihood [%g]. Skipping it.", gal_i, m2lnP_gal_i);
          continue;
        }

        if (m2lnP_gal != NULL)
          ncm_vector_set (m2lnP_gal, gal_i, m2lnP_gal_i);

        result += _nc_data_cluster_wl_eval_m2lnP_weight (dcwl, m2lnP_gal_i, 0.0);
      }
    }
    else
    {
      NcmBootstrap *bstrap = ncm_data_peek_bootstrap (data);

      const guint bsize = ncm_bootstrap_get_bsize (bstrap);
      guint i;

      #pragma omp for schedule (dynamic, 5)

      for (i = 0; i < bsize; i++)
      {
        guint gal_i = ncm_bootstrap_get (bstrap, i);

        NcGalaxySDShapeData *s_data_i = NC_GALAXY_SD_SHAPE_DATA (ncm_obj_array_peek (self->shape_data, gal_i));

        const gdouble z           = s_data_i->sdpos_data->sdz_data->z;
        const gdouble int_z       = nc_galaxy_sd_obs_redshift_integrand_eval (integrand_redshift, z, s_data_i->sdpos_data->sdz_data);
        const gdouble int_pos     = nc_galaxy_sd_position_integrand_eval (integrand_position, s_data_i->sdpos_data);
        const gdouble int_shape   = nc_galaxy_sd_shape_integrand_eval (integrand_shape, z, s_data_i);
        const gdouble P_gal_i     = int_z * int_pos * int_shape;
        const gdouble m2lnP_gal_i = P_gal_i > 0.0 ? -2.0 * log (P_gal_i) : NC_GALAXY_LOW_PROB;

        if (!gsl_finite (m2lnP_gal_i))
        {
          g_warning ("_nc_data_cluster_wl_eval_m2lnP: galaxy %d has undefined likelihood [%g]. Skipping it.", gal_i, m2lnP_gal_i);
          continue;
        }

        if (m2lnP_gal != NULL)
          ncm_vector_set (m2lnP_gal, gal_i, m2lnP_gal_i);

        result += _nc_data_cluster_wl_eval_m2lnP_weight (dcwl, m2lnP_gal_i, 0.0);
      }
    }

    nc_galaxy_sd_shape_integrand_free (integrand_shape);
    nc_galaxy_sd_position_integrand_free (integrand_position);
    nc_galaxy_sd_obs_redshift_integrand_free (integrand_redshift);
  }

  return result;
}

static void
_nc_data_cluster_wl_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataClusterWL *dcwl                  = NC_DATA_CLUSTER_WL (data);
  NcGalaxySDObsRedshift *galaxy_redshift = NC_GALAXY_SD_OBS_REDSHIFT (ncm_mset_peek (mset, nc_galaxy_sd_obs_redshift_id ()));

  if (NC_IS_GALAXY_SD_OBS_REDSHIFT_SPEC (galaxy_redshift))
    m2lnL[0] = _nc_data_cluster_wl_eval_m2lnP (dcwl, mset, NULL);
  else
    m2lnL[0] = _nc_data_cluster_wl_eval_m2lnP_integ (dcwl, mset, NULL);

  return;
}

static guint
_nc_data_cluster_wl_get_len (NcmData *data)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (data);
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  if (self->obs == NULL)
    return 0;
  else
    return nc_galaxy_wl_obs_len (self->obs);
}

static void
_nc_data_cluster_wl_load_obs (NcDataClusterWL *dcwl, NcGalaxyWLObs *obs, GPtrArray *shape_data)
{
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  guint i;

  g_ptr_array_set_size (shape_data, 0);

  for (i = 0; i < nc_galaxy_wl_obs_len (obs); i++)
  {
    NcGalaxySDObsRedshiftData *z_data = nc_galaxy_sd_obs_redshift_data_new (self->galaxy_redshift);
    NcGalaxySDPositionData *p_data    = nc_galaxy_sd_position_data_new (self->galaxy_position, z_data);
    NcGalaxySDShapeData *s_data       = nc_galaxy_sd_shape_data_new (self->galaxy_shape, p_data);

    nc_galaxy_sd_shape_data_read_row (s_data, obs, i);
    g_ptr_array_add (shape_data, s_data);

    nc_galaxy_sd_position_data_unref (p_data);
    nc_galaxy_sd_obs_redshift_data_unref (z_data);
  }
}

static void
_nc_data_cluster_wl_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataClusterWL *dcwl                        = NC_DATA_CLUSTER_WL (data);
  NcDataClusterWLPrivate * const self          = nc_data_cluster_wl_get_instance_private (dcwl);
  NcHICosmo *cosmo                             = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcWLSurfaceMassDensity *surface_mass_density = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *density_profile        = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  NcHaloPosition *halo_position                = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcGalaxySDShape *galaxy_shape                = NC_GALAXY_SD_SHAPE (ncm_mset_peek (mset, nc_galaxy_sd_shape_id ()));
  NcGalaxySDObsRedshift *galaxy_redshift       = NC_GALAXY_SD_OBS_REDSHIFT (ncm_mset_peek (mset, nc_galaxy_sd_obs_redshift_id ()));
  NcGalaxySDPosition *galaxy_position          = NC_GALAXY_SD_POSITION (ncm_mset_peek (mset, nc_galaxy_sd_position_id ()));

  g_assert ((cosmo != NULL) && (surface_mass_density != NULL) && (density_profile != NULL) && (halo_position != NULL));
  g_assert ((galaxy_shape != NULL) && (galaxy_redshift != NULL) && (galaxy_position != NULL));

  if (nc_galaxy_wl_obs_get_ellip_conv (self->obs) != nc_galaxy_sd_shape_get_ellip_conv (galaxy_shape))
    g_error ("nc_data_cluster_wl_prepare: ellip_conv mismatch.");

  self->cosmo                = cosmo;
  self->surface_mass_density = surface_mass_density;
  self->density_profile      = density_profile;
  self->halo_position        = halo_position;
  self->galaxy_shape         = galaxy_shape;
  self->galaxy_redshift      = galaxy_redshift;
  self->galaxy_position      = galaxy_position;

  nc_halo_position_prepare_if_needed (halo_position, cosmo);
  nc_wl_surface_mass_density_prepare_if_needed (surface_mass_density, cosmo);
  {
    const gboolean model_update = ncm_model_ctrl_model_update (self->ctrl_position, NCM_MODEL (galaxy_position)) ||
                                  ncm_model_ctrl_model_update (self->ctrl_redshift, NCM_MODEL (galaxy_redshift)) ||
                                  ncm_model_ctrl_model_update (self->ctrl_shape, NCM_MODEL (galaxy_shape));

    if (model_update)
      _nc_data_cluster_wl_load_obs (dcwl, self->obs, self->shape_data);

    nc_galaxy_sd_shape_prepare_data_array (galaxy_shape, mset, self->shape_data);
  }
}

/**
 * nc_data_cluster_wl_new:
 *
 * Creates a new galaxy weak lensing object.
 *
 * Returns: (transfer full): a new NcDataClusterWL.
 */
NcDataClusterWL *
nc_data_cluster_wl_new (void)
{
  NcDataClusterWL *dcwl = g_object_new (NC_TYPE_DATA_CLUSTER_WL,
                                        NULL);

  return dcwl;
}

/**
 * nc_data_cluster_wl_ref:
 * @dcwl: a #NcDataClusterWL
 *
 * Increases the reference count of @dcwl by one.
 *
 * Returns: (transfer full): @dcwl
 */
NcDataClusterWL *
nc_data_cluster_wl_ref (NcDataClusterWL *dcwl)
{
  return g_object_ref (dcwl);
}

/**
 * nc_data_cluster_wl_free:
 * @dcwl: a #NcDataClusterWL
 *
 * Atomically decrements the reference count of @dcwl by one. If the reference count drops to 0,
 * all memory allocated by @dcwl is released.
 *
 */
void
nc_data_cluster_wl_free (NcDataClusterWL *dcwl)
{
  g_object_unref (dcwl);
}

/**
 * nc_data_cluster_wl_clear:
 * @dcwl: a #NcDataClusterWL
 *
 * The reference count of @dcwl is decreased and the pointer is set to NULL.
 *
 */
void
nc_data_cluster_wl_clear (NcDataClusterWL **dcwl)
{
  g_clear_object (dcwl);
}

/**
 * nc_data_cluster_wl_set_prec:
 * @dcwl: a #NcDataClusterWL
 * @prec: precision for integral
 *
 * Sets the number of samples ndata.
 *
 */
void
nc_data_cluster_wl_set_prec (NcDataClusterWL *dcwl, gdouble prec)
{
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  self->prec = prec;
}

/**
 * nc_data_cluster_wl_set_obs:
 * @dcwl: a #NcDataClusterWL
 * @obs: a #NcGalaxyWLObs
 *
 * Sets the observables matrix @obs.
 */
void
nc_data_cluster_wl_set_obs (NcDataClusterWL *dcwl, NcGalaxyWLObs *obs)
{
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  nc_galaxy_wl_obs_clear (&self->obs);
  self->len = nc_galaxy_wl_obs_len (obs);
  self->obs = nc_galaxy_wl_obs_ref (obs);

  ncm_data_set_init (NCM_DATA (dcwl), TRUE);
}

/**
 * nc_data_cluster_wl_set_cut:
 * @dcwl: a #NcDataClusterWL
 * @r_min: minimum angular radius $\theta_\mathrm{min}$
 * @r_max: maximum angular radius $\theta_\mathrm{max}$
 *
 * Sets the radial cut in the observables.
 *
 */
void
nc_data_cluster_wl_set_cut (NcDataClusterWL *dcwl, const gdouble r_min, const gdouble r_max)
{
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  g_assert_cmpfloat (r_min, <, r_max);

  self->r_min = r_min;
  self->r_max = r_max;
  self->dr    = r_max - r_min;
}

/**
 * nc_data_cluster_wl_peek_obs:
 * @dcwl: a #NcDataClusterWL
 *
 * Gets the observables matrix.
 *
 * Returns: (transfer none): the observables matrix.
 */
NcGalaxyWLObs *
nc_data_cluster_wl_peek_obs (NcDataClusterWL *dcwl)
{
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  return self->obs;
}

/**
 * nc_data_cluster_wl_set_resample_flag:
 * @dcwl: a #NcDataClusterWL
 * @resample_flag: resample flag #NcDataClusterWLResampleFlag
 *
 * Sets flag to resample any combination of position, redshift and shape.
 *
 */
void
nc_data_cluster_wl_set_resample_flag (NcDataClusterWL *dcwl, NcDataClusterWLResampleFlag resample_flag)
{
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  if (resample_flag & ~NC_DATA_CLUSTER_WL_RESAMPLE_FLAG_ALL)
    g_error ("nc_data_cluster_wl_set_resample_flag: unknown resample flag %d.", resample_flag);

  self->resample_flag = resample_flag;
}

/**
 * nc_data_cluster_wl_get_resample_flag:
 * @dcwl: a #NcDataClusterWL
 *
 * Gets the resample flag.
 *
 * Returns: the resample flag #NcDataClusterWLResampleFlag.
 */
NcDataClusterWLResampleFlag
nc_data_cluster_wl_get_resample_flag (NcDataClusterWL *dcwl)
{
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  return self->resample_flag;
}

/**
 * nc_data_cluster_wl_peek_data_array:
 * @dcwl: a #NcDataClusterWL
 *
 * Gets the data array.
 *
 * Returns: (transfer none): the data array.
 */
NcmObjArray *
nc_data_cluster_wl_peek_data_array (NcDataClusterWL *dcwl)
{
  NcDataClusterWLPrivate * const self = nc_data_cluster_wl_get_instance_private (dcwl);

  return self->shape_data;
}

