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
  NcmObjArray *s_obs;
  NcmObjArray *s_obs_prep;
  NcmObjArray *z_obs;
  NcmObjArray *p_obs;
  NcGalaxySDShape *s_dist;
  NcGalaxySDObsRedshift *z_dist;
  NcGalaxySDPosition *p_dist;
  gboolean constructed;
  gdouble r_min;
  gdouble r_max;
  gdouble prec;
  guint len;
};

enum
{
  PROP_0,
  PROP_OBS,
  PROP_S_OBS,
  PROP_S_OBS_PREP,
  PROP_Z_OBS,
  PROP_P_OBS,
  PROP_S_DIST,
  PROP_Z_DIST,
  PROP_P_DIST,
  PROP_R_MIN,
  PROP_R_MAX,
  PROP_PREC,
  PROP_LEN,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcDataClusterWL, nc_data_cluster_wl, NCM_TYPE_DATA);

static void
nc_data_cluster_wl_init (NcDataClusterWL *dcwl)
{
  NcDataClusterWLPrivate * const self = dcwl->priv = nc_data_cluster_wl_get_instance_private (dcwl);

  self->obs         = NULL;
  self->s_obs       = NULL;
  self->s_obs_prep  = NULL;
  self->z_obs       = NULL;
  self->p_obs       = NULL;
  self->s_dist      = NULL;
  self->z_dist      = NULL;
  self->p_dist      = NULL;
  self->constructed = FALSE;
  self->r_max       = 0.0;
  self->r_min       = 0.0;
  self->prec        = 1.0e-6;
  self->len         = 0;
}

static void
nc_data_cluster_wl_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (object);
  NcDataClusterWLPrivate * const self = dcwl->priv;

  g_return_if_fail (NC_IS_DATA_CLUSTER_WL (object));

  switch (prop_id)
  {
    case PROP_OBS:
      nc_data_cluster_wl_set_obs (dcwl, g_value_dup_object (value));
      break;
    case PROP_S_OBS:
      self->s_obs = g_value_dup_boxed (value);
      break;
    case PROP_S_OBS_PREP:
      self->s_obs_prep = g_value_dup_boxed (value);
      break;
    case PROP_Z_OBS:
      self->z_obs = g_value_dup_boxed (value);
      break;
    case PROP_P_OBS:
      self->p_obs = g_value_dup_boxed (value);
      break;
    case PROP_S_DIST:
      self->s_dist = g_value_dup_object (value);
      break;
    case PROP_Z_DIST:
      self->z_dist = g_value_dup_object (value);
      break;
    case PROP_P_DIST:
      self->p_dist = g_value_dup_object (value);
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
  NcDataClusterWLPrivate * const self = dcwl->priv;

  g_return_if_fail (NC_IS_DATA_CLUSTER_WL (object));

  switch (prop_id)
  {
    case PROP_OBS:
      g_value_set_object (value, nc_data_cluster_wl_peek_obs (dcwl));
      break;
    case PROP_S_OBS:
      g_value_set_boxed (value, self->s_obs);
      break;
    case PROP_S_OBS_PREP:
      g_value_set_boxed (value, self->s_obs_prep);
      break;
    case PROP_Z_OBS:
      g_value_set_boxed (value, self->z_obs);
      break;
    case PROP_P_OBS:
      g_value_set_boxed (value, self->p_obs);
      break;
    case PROP_S_DIST:
      g_value_set_object (value, self->s_dist);
      break;
    case PROP_Z_DIST:
      g_value_set_object (value, self->z_dist);
      break;
    case PROP_P_DIST:
      g_value_set_object (value, self->p_dist);
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
  NcDataClusterWLPrivate * const self = dcwl->priv;

  if (self->obs != NULL)
    nc_galaxy_wl_obs_clear (&self->obs);

  if (self->s_obs != NULL)
    ncm_obj_array_clear (&self->s_obs);

  if (self->s_obs_prep != NULL)
    ncm_obj_array_clear (&self->s_obs_prep);

  if (self->z_obs != NULL)
    ncm_obj_array_clear (&self->z_obs);

  if (self->p_obs != NULL)
    ncm_obj_array_clear (&self->p_obs);

  nc_galaxy_sd_shape_clear (&self->s_dist);
  nc_galaxy_sd_obs_redshift_clear (&self->z_dist);
  nc_galaxy_sd_position_clear (&self->p_dist);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_wl_parent_class)->dispose (object);
}

static void
nc_data_cluster_wl_finalize (GObject *object)
{
  /*NcDataClusterWL *dcwl = NC_DATA_CLUSTER_WL (object);*/

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
    NcDataClusterWLPrivate * const self = dcwl->priv;

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
   * NcDataClusterWL:s-obs:
   *
   * Galaxy shape observables matrix.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_S_OBS,
                                   g_param_spec_boxed ("s-obs",
                                                       NULL,
                                                       "Galaxy shape observables matrix",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:s-obs-prep:
   *
   * Prepared galaxy shape observables matrix.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_S_OBS_PREP,
                                   g_param_spec_boxed ("s-obs-prep",
                                                       NULL,
                                                       "Prepared galaxy shape observables matrix",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:z-obs:
   *
   * Galaxy redshift observables matrix.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_OBS,
                                   g_param_spec_boxed ("z-obs",
                                                       NULL,
                                                       "Galaxy redshift observables matrix",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:p-obs:
   *
   * Galaxy position observables matrix.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_P_OBS,
                                   g_param_spec_boxed ("p-obs",
                                                       NULL,
                                                       "Galaxy position observables matrix",
                                                       NCM_TYPE_OBJ_ARRAY,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:s-dist:
   *
   * A #NcGalaxySDShape object.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_S_DIST,
                                   g_param_spec_object ("s-dist",
                                                        NULL,
                                                        "Galaxy sample shape distribution",
                                                        NC_TYPE_GALAXY_SD_SHAPE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:z-dist:
   *
   * A #NcGalaxySDObsRedshift object.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_Z_DIST,
                                   g_param_spec_object ("z-dist",
                                                        NULL,
                                                        "Galaxy sample observed redshift distribution",
                                                        NC_TYPE_GALAXY_SD_OBS_REDSHIFT,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcDataClusterWL:p-dist:
   *
   * A #NcGalaxySDPosition object.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_P_DIST,
                                   g_param_spec_object ("p-dist",
                                                        NULL,
                                                        "Galaxy sample position distribution",
                                                        NC_TYPE_GALAXY_SD_POSITION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

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
  NcGalaxySDPosition *p_dist;
  NcGalaxySDObsRedshift *z_dist;
  NcGalaxySDShape *s_dist;
  NcHICosmo *cosmo;
  NcHaloDensityProfile *dp;
  NcWLSurfaceMassDensity *smd;
  NcHaloPosition *hp;
  NcmVector *s_obs;
  NcmVector *z_obs;
  NcmVector *p_obs;
};

static void nc_data_cluster_wl_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
static void nc_data_cluster_wl_int_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);

NCM_INTEGRAL_ND_DEFINE_TYPE (NC, DATA_CLUSTER_WL_INT, NcDataClusterWLInt, nc_data_cluster_wl_integ, nc_data_cluster_wl_int_dim, nc_data_cluster_wl_integ, struct _NcDataClusterWLIntArg);

static void
nc_data_cluster_wl_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcDataClusterWLInt *lh_int = NC_DATA_CLUSTER_WL_INT (intnd);
  guint i;

  for (i = 0; i < npoints; i++)
  {
    const gdouble z = ncm_vector_get (x, i);
    gdouble res     = nc_galaxy_sd_position_integ (lh_int->data.p_dist, lh_int->data.p_obs) *
                      nc_galaxy_sd_obs_redshift_integ (lh_int->data.z_dist, z, lh_int->data.z_obs) *
                      nc_galaxy_sd_shape_integ_optzs (lh_int->data.s_dist, lh_int->data.cosmo, lh_int->data.dp, lh_int->data.smd, lh_int->data.hp, z, lh_int->data.s_obs);

    ncm_vector_set (fval, i, res);
  }
}

static void
nc_data_cluster_wl_int_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  *dim  = 1;
  *fdim = 1;
}

/**
 * nc_data_cluster_wl_eval_m2lnP_integ:
 * @dcwl: a #NcDataClusterWL
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @hp: a #NcHaloPosition
 * @m2lnP_gal: (out) (optional): a #NcmVector
 *
 * Computes the observables probability given the theoretical modeling using
 * integration method.
 *
 *
 * Returns: $-2\ln(P)$.
 */
gdouble
nc_data_cluster_wl_eval_m2lnP_integ (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hp, NcmVector *m2lnP_gal)
{
  NcDataClusterWLPrivate * const self     = dcwl->priv;
  NcDataClusterWLInt *likelihood_integral = g_object_new (nc_data_cluster_wl_integ_get_type (), NULL);
  NcmIntegralND *lh_int                   = NCM_INTEGRAL_ND (likelihood_integral);
  gdouble verr, vzpi, vzpf, vres, vtmp;
  NcmVector *err = ncm_vector_new_data_static (&verr, 1, 1);
  NcmVector *zpi = ncm_vector_new_data_static (&vzpi, 1, 1);
  NcmVector *zpf = ncm_vector_new_data_static (&vzpf, 1, 1);
  NcmVector *res = ncm_vector_new_data_static (&vres, 1, 1);
  NcmVector *tmp = ncm_vector_new_data_static (&vtmp, 1, 1);
  gdouble result = 0;
  guint gal_i;

  vzpi = 0.0;
  vzpf = 10.0;

  ncm_integral_nd_set_reltol (lh_int, self->prec);
  ncm_integral_nd_set_abstol (lh_int, 0.0);
  ncm_integral_nd_set_method (lh_int, NCM_INTEGRAL_ND_METHOD_CUBATURE_H_V);

  nc_galaxy_sd_shape_prepare (self->s_dist, cosmo, hp, nc_galaxy_wl_obs_get_coord (self->obs), FALSE, self->s_obs, self->s_obs_prep);

  if (m2lnP_gal != NULL)
    g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

  for (gal_i = 0; gal_i < self->len; gal_i++)
  {
    gdouble m2lnP_gal_i;

    likelihood_integral->data.p_dist = NC_GALAXY_SD_POSITION (self->p_dist);
    likelihood_integral->data.z_dist = NC_GALAXY_SD_OBS_REDSHIFT (self->z_dist);
    likelihood_integral->data.s_dist = NC_GALAXY_SD_SHAPE (self->s_dist);
    likelihood_integral->data.cosmo  = cosmo;
    likelihood_integral->data.dp     = dp;
    likelihood_integral->data.smd    = smd;
    likelihood_integral->data.hp     = hp;
    likelihood_integral->data.s_obs  = NCM_VECTOR (ncm_obj_array_peek (self->s_obs_prep, gal_i));
    likelihood_integral->data.z_obs  = NCM_VECTOR (ncm_obj_array_peek (self->z_obs, gal_i));
    likelihood_integral->data.p_obs  = NCM_VECTOR (ncm_obj_array_peek (self->p_obs, gal_i));

    nc_galaxy_sd_shape_integ_optzs_prep (self->s_dist, cosmo, dp, smd, hp, likelihood_integral->data.s_obs);

    ncm_integral_nd_eval (lh_int, zpi, zpf, res, err);

    m2lnP_gal_i = -2.0 * log (vres);

    if (m2lnP_gal != NULL)
      ncm_vector_set (m2lnP_gal, gal_i, m2lnP_gal_i);

    result += m2lnP_gal_i;
  }

  ncm_vector_free (err);
  ncm_vector_free (zpi);
  ncm_vector_free (zpf);
  ncm_vector_free (res);
  ncm_vector_free (tmp);
  ncm_integral_nd_clear (&lh_int);

  return result;
}

/**
 * nc_data_cluster_wl_eval_m2lnP:
 * @dcwl: a #NcDataClusterWL
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @hp: a #NcHaloPosition
 * @m2lnP_gal: (out) (optional): a #NcmVector
 *
 * Computes the observables probability given the theoretical modeling using
 * integration method.
 *
 *
 * Returns: $-2\ln(P)$.
 */
gdouble
nc_data_cluster_wl_eval_m2lnP (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hp, NcmVector *m2lnP_gal)
{
  NcDataClusterWLPrivate * const self = dcwl->priv;
  gdouble result                      = 0;
  guint gal_i;

  nc_galaxy_sd_shape_prepare (self->s_dist, cosmo, hp, nc_galaxy_wl_obs_get_coord (self->obs), FALSE, self->s_obs, self->s_obs_prep);

  if (m2lnP_gal != NULL)
    g_assert_cmpuint (ncm_vector_len (m2lnP_gal), ==, self->len);

  for (gal_i = 0; gal_i < self->len; gal_i++)
  {
    gdouble z           = ncm_vector_get (NCM_VECTOR (ncm_obj_array_peek (self->z_obs, gal_i)), 0);
    gdouble m2lnP_gal_i =  nc_galaxy_sd_position_integ (self->p_dist, NCM_VECTOR (ncm_obj_array_peek (self->p_obs, gal_i))) *
                          nc_galaxy_sd_obs_redshift_integ (self->z_dist, z, NCM_VECTOR (ncm_obj_array_peek (self->z_obs, gal_i))) *
                          nc_galaxy_sd_shape_integ_optzs (self->s_dist, cosmo, dp, smd, hp, z, NCM_VECTOR (ncm_obj_array_peek (self->s_obs, gal_i)));

    result += m2lnP_gal_i;
  }

  return result;
}

static void
_nc_data_cluster_wl_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (data);
  NcDataClusterWLPrivate * const self = dcwl->priv;
  NcHICosmo *cosmo                    = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcWLSurfaceMassDensity *smd         = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *dp            = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  NcHaloPosition *hp                  = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  guint i;

  if ((self->s_obs != NULL) && (self->z_obs != NULL) && (self->p_obs != NULL) && (self->s_obs_prep != NULL))
  {
    g_assert_cmpuint (ncm_obj_array_len (self->s_obs), >, 0);
    g_assert_cmpuint (ncm_obj_array_len (self->z_obs), >, 0);
    g_assert_cmpuint (ncm_obj_array_len (self->p_obs), >, 0);
  }
  else
  {
    g_error ("Weak lensing observables matrices are not set.");
  }

  m2lnL[0] = 0.0;

  if (NC_IS_GALAXY_SD_OBS_REDSHIFT_SPEC (self->z_dist))
    m2lnL[0] += nc_data_cluster_wl_eval_m2lnP (dcwl, cosmo, dp, smd, hp, NULL);
  else
    m2lnL[0] += nc_data_cluster_wl_eval_m2lnP_integ (dcwl, cosmo, dp, smd, hp, NULL);

  return;
}

static guint
_nc_data_cluster_wl_get_len (NcmData *data)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (data);
  NcDataClusterWLPrivate * const self = dcwl->priv;

  if (self->obs == NULL)
    return 0;
  else
    return nc_galaxy_wl_obs_len (self->obs);
}

static void
_nc_data_cluster_wl_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataClusterWL *dcwl               = NC_DATA_CLUSTER_WL (data);
  NcDataClusterWLPrivate * const self = dcwl->priv;
  NcHICosmo *cosmo                    = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcWLSurfaceMassDensity *smd         = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *dp            = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  NcHaloPosition *hp                  = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));

  g_assert ((cosmo != NULL) && (smd != NULL) && (dp != NULL) && (hp != NULL));

  nc_wl_surface_mass_density_prepare_if_needed (smd, cosmo);

  nc_galaxy_sd_shape_set_models (self->s_dist, cosmo, hp);

  nc_galaxy_sd_shape_prepare (self->s_dist, cosmo, hp, nc_galaxy_wl_obs_get_coord (self->obs), TRUE, self->s_obs, self->s_obs_prep);
}

/**
 * nc_data_cluster_wl_new:
 * @s_dist: a #NcGalaxySDShape
 * @z_dist: a #NcGalaxySDObsRedshift
 * @p_dist: a #NcGalaxySDPosition
 *
 * Creates a new galaxy weak lensing object.
 * Requires an instance of #NcGalaxySDShape, #NcGalaxySDObsRedshift, and #NcGalaxySDPosition.
 *
 * Returns: (transfer full): a new NcDataClusterWL.
 */
NcDataClusterWL *
nc_data_cluster_wl_new (NcGalaxySDShape *s_dist, NcGalaxySDObsRedshift *z_dist, NcGalaxySDPosition *p_dist)
{
  NcDataClusterWL *dcwl = g_object_new (NC_TYPE_DATA_CLUSTER_WL,
                                        "s-dist", s_dist,
                                        "z-dist", z_dist,
                                        "p-dist", p_dist,
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
  NcDataClusterWLPrivate * const self = dcwl->priv;

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
  NcDataClusterWLPrivate * const self = dcwl->priv;
  NcmVarDict *obs_header              = nc_galaxy_wl_obs_peek_header (obs);
  GStrv s_header                      = nc_galaxy_sd_shape_get_header (self->s_dist);
  GStrv z_header                      = nc_galaxy_sd_obs_redshift_get_header (self->z_dist);
  GStrv p_header                      = nc_galaxy_sd_position_get_header (self->p_dist);
  NcmObjArray *s_obs                  = ncm_obj_array_sized_new (nc_galaxy_wl_obs_len (obs));
  NcmObjArray *s_obs_prep             = ncm_obj_array_sized_new (nc_galaxy_wl_obs_len (obs));
  NcmObjArray *z_obs                  = ncm_obj_array_sized_new (nc_galaxy_wl_obs_len (obs));
  NcmObjArray *p_obs                  = ncm_obj_array_sized_new (nc_galaxy_wl_obs_len (obs));
  gchar **str;
  guint i;

  for (str = s_header; *str != NULL; str++)
  {
    g_assert_true (ncm_var_dict_has_key (obs_header, *str));
  }

  for (str = z_header; *str != NULL; str++)
  {
    g_assert_true (ncm_var_dict_has_key (obs_header, *str));
  }

  for (str = p_header; *str != NULL; str++)
  {
    g_assert_true (ncm_var_dict_has_key (obs_header, *str));
  }

  g_assert_cmpuint (ncm_var_dict_len (obs_header), >=, 4);
  g_assert_cmpuint (nc_galaxy_wl_obs_len (obs), >, 0);

  for (i = 0; i < nc_galaxy_wl_obs_len (obs); i++)
  {
    NcmVector *s_obs_i      = ncm_vector_new (g_strv_length (s_header));
    NcmVector *s_obs_prep_i = ncm_vector_new (nc_galaxy_sd_shape_get_vec_size (self->s_dist));
    NcmVector *z_obs_i      = ncm_vector_new (g_strv_length (z_header));
    NcmVector *p_obs_i      = ncm_vector_new (g_strv_length (p_header));
    gint j                  = 0;

    for (str = s_header; *str != NULL; str++)
    {
      gdouble val = nc_galaxy_wl_obs_get (obs, *str, i);

      ncm_vector_set (s_obs_i, j, val);
      j++;
    }

    j = 0;

    for (str = z_header; *str != NULL; str++)
    {
      gdouble val = nc_galaxy_wl_obs_get (obs, *str, i);

      ncm_vector_set (z_obs_i, j, val);
      j++;
    }

    j = 0;

    for (str = p_header; *str != NULL; str++)
    {
      gdouble val = nc_galaxy_wl_obs_get (obs, *str, i);

      ncm_vector_set (p_obs_i, j, val);
      j++;
    }

    ncm_obj_array_add (s_obs, G_OBJECT (s_obs_i));
    ncm_obj_array_add (s_obs_prep, G_OBJECT (s_obs_prep_i));
    ncm_obj_array_add (z_obs, G_OBJECT (z_obs_i));
    ncm_obj_array_add (p_obs, G_OBJECT (p_obs_i));
  }

  nc_galaxy_wl_obs_clear (&self->obs);
  ncm_obj_array_clear (&self->s_obs);
  ncm_obj_array_clear (&self->s_obs_prep);
  ncm_obj_array_clear (&self->z_obs);
  ncm_obj_array_clear (&self->p_obs);

  self->len        = nc_galaxy_wl_obs_len (obs);
  self->obs        = nc_galaxy_wl_obs_ref (obs);
  self->s_obs      = s_obs;
  self->s_obs_prep = s_obs_prep;
  self->z_obs      = z_obs;
  self->p_obs      = p_obs;
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
  NcDataClusterWLPrivate * const self = dcwl->priv;

  g_assert_cmpfloat (r_min, <, r_max);

  self->r_min = r_min;
  self->r_max = r_max;
}

/**
 * nc_data_cluster_wl_gen_obs:
 * @dcwl: a #NcDataClusterWL
 * @cosmo: a #NcHICosmo
 * @dp: a #NcHaloDensityProfile
 * @smd: a #NcWLSurfaceMassDensity
 * @hp: a #NcHaloPosition
 * @nobs: number of observables to generate
 * @rng: a #NcmRNG
 *
 * Generates @nobs observables.
 */
void
nc_data_cluster_wl_gen_obs (NcDataClusterWL *dcwl, NcHICosmo *cosmo, NcHaloDensityProfile *dp, NcWLSurfaceMassDensity *smd, NcHaloPosition *hp, guint nobs, NcmRNG *rng, NcGalaxyWLObsCoord coord)
{
  NcDataClusterWLPrivate * const self = dcwl->priv;
  GStrv s_header                      = nc_galaxy_sd_shape_get_header (self->s_dist);
  GStrv z_header                      = nc_galaxy_sd_obs_redshift_get_header (self->z_dist);
  GStrv p_header                      = nc_galaxy_sd_position_get_header (self->p_dist);
  gchar *s_string                     = g_strjoinv (" ", s_header);
  gchar *z_string                     = g_strjoinv (" ", z_header);
  gchar *p_string                     = g_strjoinv (" ", p_header);
  gchar *header_string                = g_strconcat (s_string, " ", z_string, " ", p_string, NULL);
  NcGalaxyWLObs *obs                  = nc_galaxy_wl_obs_new (coord, nobs, g_strsplit (header_string, " ", -1));
  guint i;

  for (i = 0; i < nobs; i++)
  {
    NcmVector *p_obs = ncm_vector_new (g_strv_length (p_header));
    NcmVector *z_obs = ncm_vector_new (g_strv_length (z_header));
    NcmVector *s_obs = ncm_vector_new (g_strv_length (s_header));
    gdouble theta    = 0.0;
    gdouble phi      = 0.0;
    gdouble r;
    gchar **str;
    guint j;

    /* do { */
    nc_galaxy_sd_position_gen (self->p_dist, rng, p_obs);
    nc_halo_position_polar_angles (hp, ncm_vector_get (p_obs, 0), ncm_vector_get (p_obs, 1), &theta, &phi);

    r = nc_halo_position_projected_radius (hp, cosmo, theta);
    /* } while (r < self->r_min || r > self->r_max); */

    nc_galaxy_sd_obs_redshift_gen (self->z_dist, rng, z_obs);
    nc_galaxy_sd_shape_gen (self->s_dist, cosmo, dp, smd, hp, rng, coord, p_obs, z_obs, s_obs);

    j = 0;

    for (str = s_header; *str != NULL; str++)
    {
      nc_galaxy_wl_obs_set (obs, *str, i, ncm_vector_get (s_obs, j));

      j++;
    }

    j = 0;

    for (str = z_header; *str != NULL; str++)
    {
      nc_galaxy_wl_obs_set (obs, *str, i, ncm_vector_get (z_obs, j));

      j++;
    }

    j = 0;

    for (str = p_header; *str != NULL; str++)
    {
      nc_galaxy_wl_obs_set (obs, *str, i, ncm_vector_get (p_obs, j));

      j++;
    }
  }

  nc_data_cluster_wl_set_obs (dcwl, obs);
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
  NcDataClusterWLPrivate * const self = dcwl->priv;

  return nc_galaxy_wl_obs_ref (self->obs);
}

