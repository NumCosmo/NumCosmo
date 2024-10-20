/***************************************************************************
 *            nc_data_cluster_wl.c
 *
 *  Mon Jul 27 16:10:25 2020
 *  Copyright  2020  Mariana Penna Lima
 *  <pennalima@gmail.com>
 *  Tue Jun 15 16:00:13 2023
 *  Copyright  2023  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_data_cluster_wl.c
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
 * Copyright (C) 2023 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * SECTION:nc_data_cluster_wl
 * @title: NcDataClusterWL
 * @short_description: Cluster weak lensing likelihood.
 * @stability: Unstable
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
#include "galaxy/nc_galaxy_sd_obs_redshift_gauss.h"
#include "galaxy/nc_galaxy_sd_obs_redshift_spec.h"
#include "galaxy/nc_galaxy_sd_position.h"
#include "math/ncm_integral_nd.h"
#include <math.h>
#include <gsl/gsl_math.h>

#include "nc_hicosmo.h"
#include "lss/nc_halo_position.h"
#include "lss/nc_halo_density_profile.h"
#include "lss/nc_wl_surface_mass_density.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_integrate.h"
#include "math/ncm_memory_pool.h"

struct _NcDataClusterWLPrivate
{
  NcGalaxyWLObs *obs;
  GPtrArray *shape_data;
  gboolean constructed;
  gdouble r_min;
  gdouble r_max;
  gdouble prec;
  guint len;
  NcmModelCtrl *ctrl_redshift;
  NcmModelCtrl *ctrl_position;
  NcmModelCtrl *ctrl_shape;
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

  self->obs           = NULL;
  self->shape_data    = g_ptr_array_new ();
  self->constructed   = FALSE;
  self->r_max         = 0.0;
  self->r_min         = 0.0;
  self->prec          = 1.0e-6;
  self->len           = 0;
  self->ctrl_redshift = ncm_model_ctrl_new (NULL);
  self->ctrl_position = ncm_model_ctrl_new (NULL);
  self->ctrl_shape    = ncm_model_ctrl_new (NULL);

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
      nc_data_cluster_wl_set_obs (dcwl, g_value_dup_object (value));
      break;
    case PROP_R_MIN:
      self->r_min = g_value_get_double (value);

      if (self->constructed)
        g_assert_cmpfloat (self->r_min, <, self->r_max);

      break;
    case PROP_R_MAX:
      self->r_max = g_value_get_double (value);

      if (self->constructed)
        g_assert_cmpfloat (self->r_min, <, self->r_max);

      break;
    case PROP_PREC:
      nc_data_cluster_wl_set_prec (dcwl, g_value_get_double (value));
      break;
    case PROP_LEN:
      self->len = g_value_get_uint (value);
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

    self->constructed = TRUE;
  }
}

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
   * NcDataClusterWL:theta-min:
   *
   * Minimum radius of the weak lensing observables.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_R_MIN,
                                   g_param_spec_double ("theta-min",
                                                        NULL,
                                                        "Minimum radius of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:theta-max:
   *
   * Maximum theta of the weak lensing observables.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_R_MAX,
                                   g_param_spec_double ("r-max",
                                                        NULL,
                                                        "Maximum radius of the weak lensing observables",
                                                        0.0, G_MAXDOUBLE, 10.0,
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
    const gdouble z = ncm_vector_fast_get (x, i);
    gdouble res     = nc_galaxy_sd_obs_redshift_integrand_eval (integrand_redshift, z, data->sdpos_data->sdz_data) *
                      nc_galaxy_sd_position_integrand_eval (integrand_position, data->sdpos_data) *
                      nc_galaxy_sd_shape_integrand_eval (integrand_shape, z, data);

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
_nc_data_cluster_wl_eval_m2lnP_integ (NcDataClusterWL *dcwl, NcmMSet *mset, NcmVector *m2lnP_gal)
{
  NcDataClusterWLPrivate * const self                = nc_data_cluster_wl_get_instance_private (dcwl);
  NcDataClusterWLInt *likelihood_integral            = g_object_new (nc_data_cluster_wl_integ_get_type (), NULL);
  NcmIntegralND *lh_int                              = NCM_INTEGRAL_ND (likelihood_integral);
  NcGalaxySDObsRedshiftIntegrand *integrand_redshift = nc_galaxy_sd_obs_redshift_integ (self->galaxy_redshift);
  NcGalaxySDPositionIntegrand *integrand_position    = nc_galaxy_sd_position_integ (self->galaxy_position);
  NcGalaxySDShapeIntegrand *integrand_shape          = nc_galaxy_sd_shape_integ (self->galaxy_shape);
  gdouble result                                     = 0;
  guint gal_i;

  ncm_vector_fast_set (self->zpi, 0, 0.0);
  ncm_vector_fast_set (self->zpf, 0, 10.0);

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

  for (gal_i = 0; gal_i < self->len; gal_i++)
  {
    NcGalaxySDShapeData *data = NC_GALAXY_SD_SHAPE_DATA (ncm_obj_array_peek (self->shape_data, gal_i));
    gdouble m2lnP_gal_i;

    likelihood_integral->data.data = data;

    ncm_integral_nd_eval (lh_int, self->zpi, self->zpf, self->res, self->err);

    m2lnP_gal_i = -2.0 * log (ncm_vector_fast_get (self->res, 0));

    if (m2lnP_gal != NULL)
      ncm_vector_set (m2lnP_gal, gal_i, m2lnP_gal_i);

    result += m2lnP_gal_i;
  }

  ncm_integral_nd_clear (&lh_int);
  nc_galaxy_sd_shape_integrand_free (integrand_shape);
  nc_galaxy_sd_position_integrand_free (integrand_position);
  nc_galaxy_sd_obs_redshift_integrand_free (integrand_redshift);

  return result;
}

static gdouble
_nc_data_cluster_wl_eval_m2lnP (NcDataClusterWL *dcwl, NcmMSet *mset, NcmVector *m2lnP_gal)
{
  NcDataClusterWLPrivate * const self                = nc_data_cluster_wl_get_instance_private (dcwl);
  NcGalaxySDObsRedshiftIntegrand *integrand_redshift = nc_galaxy_sd_obs_redshift_integ (self->galaxy_redshift);
  NcGalaxySDPositionIntegrand *integrand_position    = nc_galaxy_sd_position_integ (self->galaxy_position);
  NcGalaxySDShapeIntegrand *integrand_shape          = nc_galaxy_sd_shape_integ (self->galaxy_shape);
  gdouble result                                     = 0;
  guint gal_i;

  if (m2lnP_gal != NULL)
    g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

  nc_galaxy_sd_obs_redshift_integrand_prepare (integrand_redshift, mset);
  nc_galaxy_sd_position_integrand_prepare (integrand_position, mset);
  nc_galaxy_sd_shape_integrand_prepare (integrand_shape, mset);

  for (gal_i = 0; gal_i < self->len; gal_i++)
  {
    NcGalaxySDShapeData *data = NC_GALAXY_SD_SHAPE_DATA (ncm_obj_array_peek (self->shape_data, gal_i));
    const gdouble z           = data->sdpos_data->sdz_data->z;
    gdouble m2lnP_gal_i       = nc_galaxy_sd_position_integrand_eval (integrand_position, data->sdpos_data) *
                                nc_galaxy_sd_obs_redshift_integrand_eval (integrand_redshift, z, data->sdpos_data->sdz_data) *
                                nc_galaxy_sd_shape_integrand_eval (integrand_shape, z, data);

    result += m2lnP_gal_i;
  }

  nc_galaxy_sd_shape_integrand_free (integrand_shape);
  nc_galaxy_sd_position_integrand_free (integrand_position);
  nc_galaxy_sd_obs_redshift_integrand_free (integrand_redshift);

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

  self->cosmo                = cosmo;
  self->surface_mass_density = surface_mass_density;
  self->density_profile      = density_profile;
  self->halo_position        = halo_position;
  self->galaxy_shape         = galaxy_shape;
  self->galaxy_redshift      = galaxy_redshift;
  self->galaxy_position      = galaxy_position;

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

  self->len = nc_galaxy_wl_obs_len (obs);
  self->obs = nc_galaxy_wl_obs_ref (obs);
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

