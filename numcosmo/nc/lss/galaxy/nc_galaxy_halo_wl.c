/***************************************************************************
 *            nc_galaxy_halo_wl.c
 *
 *  Fri November 27 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_halo_wl.c
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyHaloWL:
 *
 * The weak-lensing quantities of galaxies with respect to a halo.
 *
 * A galaxy at sky position (ra, dec) lies at a projected radius $R$ and position angle
 * $\phi$ around the halo center (#NcHaloPosition). This object computes, for each
 * galaxy, $R$, the rotation factor $e^{2i\phi}$ in the celestial frame, and the surface
 * density $\Sigma(R)$ and its mean $\overline{\Sigma}(<R)$ through its #NcWLClusterLens,
 * and keeps them in the galaxy's #NcGalaxyHaloWLData. It is not a likelihood factor:
 * #NcDataClusterWLFactor prepares it with its model set, compares its change values and
 * updates the galaxies whose quantities changed.
 *
 * It also provides the quantities of the reachable-range bound: the smallest radius the
 * prior of the halo center allows and the convergence and reduced shear at the corner of
 * the prior box, where the mass is at its upper bound.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_halo_wl.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor_moments_private.h"
#include "nc/lss/halo/nc_halo_position.h"
#include "nc/lss/halo/nc_halo_density_profile.h"
#include "nc/lss/wl/nc_wl_surface_mass_density.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/model/ncm_model_ctrl.h"
#include "ncm/algebra/ncm_quaternion.h"
#include "ncm/core/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyHaloWLData
{
  gdouble ra_v; /* the position of v, NaN before it is set */
  gdouble dec_v;
  NcmTriVec v; /* the unit vector of (ra_v, dec_v) */
  gdouble radius;
  gdouble rot_re;
  gdouble rot_im;
  gdouble sigma;
  gdouble mean_sigma;
  gatomicrefcount ref_count;
};

struct _NcGalaxyHaloWL
{
  /*< private >*/
  GObject parent_instance;
  NcWLClusterLens *lens;
  NcmQuaternion q;
  gdouble pr_prefactor;
  NcmModelCtrl *ctrl_cosmo;
  NcmModelCtrl *ctrl_hp;
  guint64 radius_hash;   /* the prepares that changed the halo position or the cosmology */
  NcmModelCtrl *ctrl_dp; /* the profile object dp_box was copied from */
  NcHaloDensityProfile *dp_box;
  gdouble box_ra_c;
  gdouble box_ra_hw;
  gdouble box_dec_lo;
  gdouble box_dec_hi;
  gboolean c_box;
  guint c_idx;
  gdouble c_lo;
  gdouble c_hi;
  gboolean bound_ready;
  gboolean bound_holds;
  GMutex bound_lock;
};

enum
{
  PROP_0,
  PROP_LENS,
};

G_DEFINE_BOXED_TYPE (NcGalaxyHaloWLData, nc_galaxy_halo_wl_data, nc_galaxy_halo_wl_data_ref, nc_galaxy_halo_wl_data_unref)
G_DEFINE_TYPE (NcGalaxyHaloWL, nc_galaxy_halo_wl, G_TYPE_OBJECT)

static void
nc_galaxy_halo_wl_init (NcGalaxyHaloWL *hwl)
{
  hwl->lens         = NULL;
  hwl->pr_prefactor = GSL_NAN;
  hwl->ctrl_cosmo   = ncm_model_ctrl_new (NULL);
  hwl->ctrl_hp      = ncm_model_ctrl_new (NULL);
  hwl->radius_hash  = 0;
  hwl->ctrl_dp      = ncm_model_ctrl_new (NULL);
  hwl->dp_box       = NULL;
  ncm_quaternion_set_I (&hwl->q);
  hwl->bound_ready = FALSE;
  hwl->bound_holds = FALSE;
  g_mutex_init (&hwl->bound_lock);
}

static void
_nc_galaxy_halo_wl_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyHaloWL *hwl = NC_GALAXY_HALO_WL (object);

  switch (prop_id)
  {
    case PROP_LENS:
      nc_wl_cluster_lens_clear (&hwl->lens);
      hwl->lens = g_value_dup_object (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_halo_wl_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyHaloWL *hwl = NC_GALAXY_HALO_WL (object);

  switch (prop_id)
  {
    case PROP_LENS:
      g_value_set_object (value, hwl->lens);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_halo_wl_dispose (GObject *object)
{
  NcGalaxyHaloWL *hwl = NC_GALAXY_HALO_WL (object);

  nc_wl_cluster_lens_clear (&hwl->lens);
  ncm_model_ctrl_clear (&hwl->ctrl_cosmo);
  ncm_model_ctrl_clear (&hwl->ctrl_hp);
  ncm_model_ctrl_clear (&hwl->ctrl_dp);
  nc_halo_density_profile_clear (&hwl->dp_box);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_halo_wl_parent_class)->dispose (object);
}

static void
_nc_galaxy_halo_wl_finalize (GObject *object)
{
  NcGalaxyHaloWL *hwl = NC_GALAXY_HALO_WL (object);

  g_mutex_clear (&hwl->bound_lock);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_halo_wl_parent_class)->finalize (object);
}

static void
nc_galaxy_halo_wl_class_init (NcGalaxyHaloWLClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_halo_wl_set_property;
  object_class->get_property = &_nc_galaxy_halo_wl_get_property;
  object_class->dispose      = &_nc_galaxy_halo_wl_dispose;
  object_class->finalize     = &_nc_galaxy_halo_wl_finalize;

  /**
   * NcGalaxyHaloWL:lens:
   *
   * The lens of the halo.
   */
  g_object_class_install_property (object_class,
                                   PROP_LENS,
                                   g_param_spec_object ("lens",
                                                        NULL,
                                                        "Cluster lens",
                                                        NC_TYPE_WL_CLUSTER_LENS,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_galaxy_halo_wl_data_new:
 *
 * Creates the weak-lensing quantities of one galaxy, NaN until
 * nc_galaxy_halo_wl_update_data_radius() and nc_galaxy_halo_wl_update_data_profile().
 *
 * Returns: (transfer full): a new #NcGalaxyHaloWLData.
 */
NcGalaxyHaloWLData *
nc_galaxy_halo_wl_data_new (void)
{
  NcGalaxyHaloWLData *data = g_new0 (NcGalaxyHaloWLData, 1);

  data->ra_v       = GSL_NAN;
  data->dec_v      = GSL_NAN;
  data->radius     = GSL_NAN;
  data->rot_re     = GSL_NAN;
  data->rot_im     = GSL_NAN;
  data->sigma      = GSL_NAN;
  data->mean_sigma = GSL_NAN;
  g_atomic_ref_count_init (&data->ref_count);

  return data;
}

/**
 * nc_galaxy_halo_wl_data_ref:
 * @data: a #NcGalaxyHaloWLData
 *
 * Increases the reference count of @data by one.
 *
 * Returns: (transfer full): @data.
 */
NcGalaxyHaloWLData *
nc_galaxy_halo_wl_data_ref (NcGalaxyHaloWLData *data)
{
  g_atomic_ref_count_inc (&data->ref_count);

  return data;
}

/**
 * nc_galaxy_halo_wl_data_unref:
 * @data: a #NcGalaxyHaloWLData
 *
 * Decreases the reference count of @data by one, freeing it at zero.
 */
void
nc_galaxy_halo_wl_data_unref (NcGalaxyHaloWLData *data)
{
  if (g_atomic_ref_count_dec (&data->ref_count))
    g_free (data);
}

/**
 * nc_galaxy_halo_wl_data_get_radius:
 * @data: a #NcGalaxyHaloWLData
 *
 * Returns: the projected radius $R$ of the galaxy, in Mpc.
 */
gdouble
nc_galaxy_halo_wl_data_get_radius (NcGalaxyHaloWLData *data)
{
  return data->radius;
}

/**
 * nc_galaxy_halo_wl_data_get_rotation:
 * @data: a #NcGalaxyHaloWLData
 * @rot_re: (out): $\cos 2\phi$
 * @rot_im: (out): $\sin 2\phi$
 *
 * Gets the rotation factor $e^{2i\phi}$ of the galaxy in the celestial frame.
 */
void
nc_galaxy_halo_wl_data_get_rotation (NcGalaxyHaloWLData *data, gdouble *rot_re, gdouble *rot_im)
{
  *rot_re = data->rot_re;
  *rot_im = data->rot_im;
}

/**
 * nc_galaxy_halo_wl_data_get_sigma:
 * @data: a #NcGalaxyHaloWLData
 * @sigma: (out): the surface density $\Sigma(R)$
 * @mean_sigma: (out): the mean surface density $\overline{\Sigma}(<R)$
 *
 * Gets $\Sigma(R)$ and $\overline{\Sigma}(<R)$ of the galaxy, in $M_\odot/\mathrm{Mpc}^2$.
 */
void
nc_galaxy_halo_wl_data_get_sigma (NcGalaxyHaloWLData *data, gdouble *sigma, gdouble *mean_sigma)
{
  *sigma      = data->sigma;
  *mean_sigma = data->mean_sigma;
}

/**
 * nc_galaxy_halo_wl_new:
 * @lens: a #NcWLClusterLens
 *
 * Creates a new #NcGalaxyHaloWL computing the lensing quantities with @lens.
 *
 * Returns: (transfer full): a new #NcGalaxyHaloWL.
 */
NcGalaxyHaloWL *
nc_galaxy_halo_wl_new (NcWLClusterLens *lens)
{
  return g_object_new (NC_TYPE_GALAXY_HALO_WL,
                       "lens", lens,
                       NULL);
}

/**
 * nc_galaxy_halo_wl_ref:
 * @hwl: a #NcGalaxyHaloWL
 *
 * Increases the reference count of @hwl by one.
 *
 * Returns: (transfer full): @hwl.
 */
NcGalaxyHaloWL *
nc_galaxy_halo_wl_ref (NcGalaxyHaloWL *hwl)
{
  return g_object_ref (hwl);
}

/**
 * nc_galaxy_halo_wl_free:
 * @hwl: a #NcGalaxyHaloWL
 *
 * Decreases the reference count of @hwl by one.
 */
void
nc_galaxy_halo_wl_free (NcGalaxyHaloWL *hwl)
{
  g_object_unref (hwl);
}

/**
 * nc_galaxy_halo_wl_clear:
 * @hwl: a #NcGalaxyHaloWL
 *
 * Decreases the reference count of *@hwl by one and sets the pointer *@hwl to %NULL.
 */
void
nc_galaxy_halo_wl_clear (NcGalaxyHaloWL **hwl)
{
  g_clear_object (hwl);
}

/**
 * nc_galaxy_halo_wl_peek_lens:
 * @hwl: a #NcGalaxyHaloWL
 *
 * Returns: (transfer none): the #NcWLClusterLens of @hwl.
 */
NcWLClusterLens *
nc_galaxy_halo_wl_peek_lens (NcGalaxyHaloWL *hwl)
{
  return hwl->lens;
}

/**
 * nc_galaxy_halo_wl_prepare:
 * @hwl: a #NcGalaxyHaloWL
 * @mset: a #NcmMSet with a #NcHICosmo, a #NcHaloPosition and a #NcHaloDensityProfile
 *
 * Prepares the lens of @hwl with the models of @mset. When the halo position or the
 * cosmology changed, it copies the rotation that takes the halo to the positive z axis
 * (nc_halo_position_get_rotation()) and computes the projected-radius prefactor, the
 * two values the radius updates of the galaxies read. No model is kept.
 */
void
nc_galaxy_halo_wl_prepare (NcGalaxyHaloWL *hwl, NcmMSet *mset)
{
  NcHICosmo *cosmo         = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *hp       = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcHaloDensityProfile *dp = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  gboolean cosmo_changed, hp_changed;

  g_assert_nonnull (cosmo);
  g_assert_nonnull (hp);
  g_assert_nonnull (dp);

  nc_wl_cluster_lens_prepare (hwl->lens, cosmo, hp, dp);

  cosmo_changed = ncm_model_ctrl_update (hwl->ctrl_cosmo, NCM_MODEL (cosmo));
  hp_changed    = ncm_model_ctrl_update (hwl->ctrl_hp, NCM_MODEL (hp));

  if (cosmo_changed || hp_changed)
  {
    nc_halo_position_get_rotation (hp, &hwl->q);
    hwl->pr_prefactor = nc_halo_position_projected_radius_prefactor (hp, cosmo);
    hwl->radius_hash++;
  }
}

/**
 * nc_galaxy_halo_wl_get_radius_hash:
 * @hwl: a #NcGalaxyHaloWL
 *
 * Returns: the number of prepares that changed the halo position or the cosmology, so
 *   $R$ and $e^{2i\phi}$; a different value means a change. A model changed when it is
 *   another object or its state changed (#NcmModelCtrl).
 */
guint64
nc_galaxy_halo_wl_get_radius_hash (NcGalaxyHaloWL *hwl)
{
  return hwl->radius_hash;
}

/**
 * nc_galaxy_halo_wl_get_profile_hash:
 * @hwl: a #NcGalaxyHaloWL
 *
 * Returns: a value that changes when a prepare changed $\Sigma(R)$ and
 *   $\overline{\Sigma}(<R)$ (nc_wl_cluster_lens_get_profile_hash()).
 */
guint64
nc_galaxy_halo_wl_get_profile_hash (NcGalaxyHaloWL *hwl)
{
  return nc_wl_cluster_lens_get_profile_hash (hwl->lens);
}

/**
 * nc_galaxy_halo_wl_update_data_radius:
 * @hwl: a #NcGalaxyHaloWL
 * @data: a #NcGalaxyHaloWLData
 * @pos_data: the #NcGalaxyPositionFactorData of the galaxy
 *
 * Computes the projected radius $R$ and the rotation factor $e^{2i\phi}$ of the galaxy
 * at @pos_data around the halo of the last prepare. The unit vector of the position is
 * kept and recomputed only when the position changed.
 */
void
nc_galaxy_halo_wl_update_data_radius (NcGalaxyHaloWL *hwl, NcGalaxyHaloWLData *data, NcGalaxyPositionFactorData *pos_data)
{
  const gdouble ra  = pos_data->ra;
  const gdouble dec = pos_data->dec;
  NcmTriVec v;

  if ((data->ra_v != ra) || (data->dec_v != dec))
  {
    ncm_trivec_set_astro_ra_dec (&data->v, 1.0, ra, dec);
    data->ra_v  = ra;
    data->dec_v = dec;
  }

  v = data->v;
  ncm_quaternion_rotate (&hwl->q, &v);

  /*
   * With the halo on the z axis, the polar angle theta and the azimuth phi of
   * nc_halo_position_polar_angles() are those of v: sin (theta) = rho / |v| and
   * (cos (2 phi), sin (2 phi)) = (x^2 - y^2, 2 x y) / rho^2, rho^2 = x^2 + y^2, with
   * phi = 0 at rho = 0 as there.
   */
  {
    const gdouble x    = v.c[0];
    const gdouble y    = v.c[1];
    const gdouble z    = v.c[2];
    const gdouble rho2 = x * x + y * y;
    const gdouble rho  = sqrt (rho2);

    data->radius = hwl->pr_prefactor * (rho / sqrt (rho2 + z * z));

    if (rho2 > 0.0)
    {
      data->rot_re = (x * x - y * y) / rho2;
      data->rot_im = 2.0 * x * y / rho2;
    }

    /*
     * LCOV_EXCL_START: rho2 is zero only for a vector exactly on the z axis, which the
     * rounding of the unit vector of (ra, dec) and of the rotation does not produce; it
     * guards the division.
     */
    else
    {
      data->rot_re = 1.0;
      data->rot_im = 0.0;
    }

    /* LCOV_EXCL_STOP */
  }
}

/**
 * nc_galaxy_halo_wl_update_data_profile:
 * @hwl: a #NcGalaxyHaloWL
 * @dp: the #NcHaloDensityProfile of the last prepare
 * @data: a #NcGalaxyHaloWLData, with its radius set
 *
 * Computes $\Sigma(R)$ and $\overline{\Sigma}(<R)$ of the galaxy with the lens and @dp.
 * Galaxies can be updated concurrently once @dp has been evaluated after its last
 * parameter change; a numerically integrated profile builds its splines on that first
 * call.
 */
void
nc_galaxy_halo_wl_update_data_profile (NcGalaxyHaloWL *hwl, NcHaloDensityProfile *dp, NcGalaxyHaloWLData *data)
{
  nc_wl_cluster_lens_sigma (hwl->lens, dp, data->radius, &data->sigma, &data->mean_sigma);
}

/**
 * nc_galaxy_halo_wl_prepare_bound:
 * @hwl: a #NcGalaxyHaloWL, prepared
 * @mset: the #NcmMSet of the last prepare
 *
 * Prepares the reachable-range bound. It sets a copy of the profile of @mset at the corner
 * of the prior box, with mass at its upper bound. It records the centre and concentration
 * boxes. The free parameters may include the centre, mass, concentration and shape
 * population. The
 * bound is not ready when @mset has no #NcWLSurfaceMassDensity.
 */
void
nc_galaxy_halo_wl_prepare_bound (NcGalaxyHaloWL *hwl, NcmMSet *mset)
{
  NcHICosmo *cosmo            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *hp          = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));

  hwl->bound_ready = FALSE;
  hwl->bound_holds = FALSE;

  if (smd == NULL)
    return;

  if (ncm_model_ctrl_model_update (hwl->ctrl_dp, NCM_MODEL (dp)))
  {
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

    nc_halo_density_profile_clear (&hwl->dp_box);
    hwl->dp_box = NC_HALO_DENSITY_PROFILE (ncm_model_dup (NCM_MODEL (dp), ser));
    ncm_serialize_free (ser);
  }

  nc_wl_surface_mass_density_prepare_if_needed (smd, cosmo);

  {
    NcmModel *hp_model      = NCM_MODEL (hp);
    const gboolean ra_free  = ncm_model_param_get_ftype (hp_model, NC_HALO_POSITION_RA) == NCM_PARAM_TYPE_FREE;
    const gboolean dec_free = ncm_model_param_get_ftype (hp_model, NC_HALO_POSITION_DEC) == NCM_PARAM_TYPE_FREE;
    const gdouble ra        = ncm_model_param_get (hp_model, NC_HALO_POSITION_RA);
    const gdouble dec       = ncm_model_param_get (hp_model, NC_HALO_POSITION_DEC);
    const gdouble ra_lo     = ra_free ? ncm_model_param_get_lower_bound (hp_model, NC_HALO_POSITION_RA) : ra;
    const gdouble ra_hi     = ra_free ? ncm_model_param_get_upper_bound (hp_model, NC_HALO_POSITION_RA) : ra;
    const gdouble dec_lo    = dec_free ? ncm_model_param_get_lower_bound (hp_model, NC_HALO_POSITION_DEC) : dec;
    const gdouble dec_hi    = dec_free ? ncm_model_param_get_upper_bound (hp_model, NC_HALO_POSITION_DEC) : dec;

    hwl->box_ra_c   = 0.5 * (ra_lo + ra_hi);
    hwl->box_ra_hw  = 0.5 * fabs (ra_hi - ra_lo);
    hwl->box_dec_lo = MIN (dec_lo, dec_hi);
    hwl->box_dec_hi = MAX (dec_lo, dec_hi);
  }

  hwl->bound_ready = _nc_galaxy_shape_factor_moments_set_box_corner (hwl->dp_box, dp);
  hwl->bound_holds = _nc_galaxy_shape_factor_moments_range_holds (mset, TRUE);
  hwl->c_box       = _nc_galaxy_shape_factor_moments_cdelta_box (dp, &hwl->c_idx, &hwl->c_lo, &hwl->c_hi) &&
                     (hwl->c_hi > hwl->c_lo);
}

/**
 * nc_galaxy_halo_wl_bound_ready:
 * @hwl: a #NcGalaxyHaloWL
 *
 * Returns: whether the last nc_galaxy_halo_wl_prepare_bound() set the corner of the
 *   prior box.
 */
gboolean
nc_galaxy_halo_wl_bound_ready (NcGalaxyHaloWL *hwl)
{
  return hwl->bound_ready;
}

/**
 * nc_galaxy_halo_wl_bound_holds:
 * @hwl: a #NcGalaxyHaloWL
 *
 * Returns: whether the free parameters of the last nc_galaxy_halo_wl_prepare_bound()
 *   are within the supported centre, mass, concentration and shape population parameters.
 */
gboolean
nc_galaxy_halo_wl_bound_holds (NcGalaxyHaloWL *hwl)
{
  return hwl->bound_holds;
}

#define NC_GALAXY_HALO_WL_C_GRID 17
#define NC_GALAXY_HALO_WL_C_GOLDEN 24

/* |g_t| at the corner at concentration c (when c_box) under bound_lock; FALSE when kappa >= 1 */
static gboolean
_nc_galaxy_halo_wl_corner_eval (NcGalaxyHaloWL *hwl, NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo, const gdouble z_cl, const gdouble R, const gdouble z, const gdouble c, gdouble *g_t, gdouble *max_kappa)
{
  gdouble kappa, g;

  if (hwl->c_box)
    ncm_model_param_set (NCM_MODEL (nc_halo_density_profile_peek_mass_summary (hwl->dp_box)), hwl->c_idx, c);

  kappa = nc_wl_surface_mass_density_convergence (smd, hwl->dp_box, cosmo, R, z, z_cl, z_cl);
  g     = fabs (nc_wl_surface_mass_density_reduced_shear (smd, hwl->dp_box, cosmo, R, z, z_cl, z_cl));

  if (!gsl_finite (kappa) || !gsl_finite (g) || !(kappa < 1.0))
    return FALSE;

  *max_kappa = MAX (*max_kappa, kappa);
  g_t[0]     = g;

  return TRUE;
}

/* |g_t| at the corner at concentration c, clearing ok when kappa >= 1 there */
static gdouble
_nc_galaxy_halo_wl_corner_eval_ok (NcGalaxyHaloWL *hwl, NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo, const gdouble z_cl, const gdouble R, const gdouble z, const gdouble c, gboolean *ok, gdouble *max_kappa)
{
  gdouble g = 0.0;

  if (!_nc_galaxy_halo_wl_corner_eval (hwl, smd, cosmo, z_cl, R, z, c, &g, max_kappa))
    ok[0] = FALSE;

  return g;
}

/* The largest |g_t| at the corner over the concentrations of the fit; FALSE when kappa >= 1 at
 * any of the concentrations evaluated. The golden-section search stays between two nodes of
 * the grid where kappa < 1, so it fails only if kappa is not monotonic in c there. */
static gboolean
_nc_galaxy_halo_wl_corner_g_t (NcGalaxyHaloWL *hwl, NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo, const gdouble z_cl, const gdouble R, const gdouble z, gdouble *g_t, gdouble *max_kappa)
{
  const guint n         = NC_GALAXY_HALO_WL_C_GRID;
  const gdouble lc_lo   = log (hwl->c_lo);
  const gdouble dlc     = (log (hwl->c_hi) - lc_lo) / (n - 1.0);
  const gdouble inv_phi = 0.5 * (sqrt (5.0) - 1.0);
  gboolean ok           = TRUE;
  gdouble best          = 0.0, a, b, x1, x2, g1, g2;
  guint k, k_best       = 0, i;

  if (!hwl->c_box)
    return _nc_galaxy_halo_wl_corner_eval (hwl, smd, cosmo, z_cl, R, z, 0.0, g_t, max_kappa);

  for (k = 0; k < n; k++)
  {
    const gdouble g = _nc_galaxy_halo_wl_corner_eval_ok (hwl, smd, cosmo, z_cl, R, z, exp (lc_lo + dlc * k), &ok, max_kappa);

    if (g > best)
    {
      best   = g;
      k_best = k;
    }
  }

  if (!ok)
    return FALSE;

  a  = lc_lo + dlc * ((k_best > 0) ? k_best - 1.0 : 0.0);
  b  = lc_lo + dlc * ((k_best + 1 < n) ? k_best + 1.0 : n - 1.0);
  x1 = b - inv_phi * (b - a);
  x2 = a + inv_phi * (b - a);
  g1 = _nc_galaxy_halo_wl_corner_eval_ok (hwl, smd, cosmo, z_cl, R, z, exp (x1), &ok, max_kappa);
  g2 = _nc_galaxy_halo_wl_corner_eval_ok (hwl, smd, cosmo, z_cl, R, z, exp (x2), &ok, max_kappa);

  for (i = 0; i < NC_GALAXY_HALO_WL_C_GOLDEN; i++)
  {
    best = GSL_MAX (best, GSL_MAX (g1, g2));

    if (g1 > g2)
    {
      b  = x2;
      x2 = x1;
      g2 = g1;
      x1 = b - inv_phi * (b - a);
      g1 = _nc_galaxy_halo_wl_corner_eval_ok (hwl, smd, cosmo, z_cl, R, z, exp (x1), &ok, max_kappa);
    }
    else
    {
      a  = x1;
      x1 = x2;
      g1 = g2;
      x2 = a + inv_phi * (b - a);
      g2 = _nc_galaxy_halo_wl_corner_eval_ok (hwl, smd, cosmo, z_cl, R, z, exp (x2), &ok, max_kappa);
    }
  }

  g_t[0] = GSL_MAX (best, GSL_MAX (g1, g2));

  return ok;
}

/*
 * The smallest angle [rad] between the direction (ra, dec) [deg] and the box of the halo
 * center. With ra inside the right ascensions of the box (modulo 360) the nearest point is
 * on the meridian of the galaxy, at its declination clamped to the box. Otherwise it is on
 * the nearer edge meridian, a great circle, at the declination that maximizes
 * sin(dec) sin(d) + cos(dec) cos(d) cos(Delta), Delta the difference of right ascension:
 * d = atan2 (sin dec, cos dec cos Delta), clamped to the box.
 */
static gdouble
_nc_galaxy_halo_wl_box_min_angle (NcGalaxyHaloWL *hwl, const gdouble ra, const gdouble dec)
{
  const gdouble d_ra = fabs (remainder (ra - hwl->box_ra_c, 360.0));
  const gdouble dr   = dec * M_PI / 180.0;

  if (d_ra <= hwl->box_ra_hw)
  {
    const gdouble d = GSL_MAX (hwl->box_dec_lo, GSL_MIN (hwl->box_dec_hi, dec));

    return fabs (dec - d) * M_PI / 180.0;
  }
  else
  {
    const gdouble Delta = (d_ra - hwl->box_ra_hw) * M_PI / 180.0;
    const gdouble d_opt = atan2 (sin (dr), cos (dr) * cos (Delta)) * 180.0 / M_PI;
    const gdouble d     = GSL_MAX (hwl->box_dec_lo, GSL_MIN (hwl->box_dec_hi, d_opt)) * M_PI / 180.0;
    const gdouble h     = gsl_pow_2 (sin (0.5 * (dr - d))) + cos (dr) * cos (d) * gsl_pow_2 (sin (0.5 * Delta));

    return 2.0 * asin (GSL_MIN (1.0, sqrt (h)));
  }
}

/**
 * nc_galaxy_halo_wl_get_radius_min:
 * @hwl: a #NcGalaxyHaloWL
 * @data: a #NcGalaxyHaloWLData, with its radius set
 *
 * Returns: the smallest projected radius of the galaxy over the prior of the halo
 *   center, using the nearest sky direction in its box.
 */
gdouble
nc_galaxy_halo_wl_get_radius_min (NcGalaxyHaloWL *hwl, NcGalaxyHaloWLData *data)
{
  if ((hwl->box_ra_hw == 0.0) && (hwl->box_dec_lo == hwl->box_dec_hi))
    return data->radius;

  return ncm_util_projected_radius (_nc_galaxy_halo_wl_box_min_angle (hwl, data->ra_v, data->dec_v), hwl->pr_prefactor);
}

/**
 * nc_galaxy_halo_wl_eval_bound:
 * @hwl: a #NcGalaxyHaloWL, with its bound ready
 * @smd: the #NcWLSurfaceMassDensity of the last nc_galaxy_halo_wl_prepare_bound()
 * @cosmo: the #NcHICosmo of the last nc_galaxy_halo_wl_prepare_bound()
 * @R_min: the projected radius, in Mpc
 * @z_s: the source redshift
 * @kappa: (out): the convergence at the corner of the prior box
 * @g_t: (out): the absolute value of the reduced shear there
 *
 * Evaluates the convergence and $|g_t|$ of the profile at the corner of the prior box,
 * at @R_min and @z_s. Thread-safe: a mutex serializes the two evaluations, since the
 * corner profile is shared by the galaxies and updates its state lazily.
 */
void
nc_galaxy_halo_wl_eval_bound (NcGalaxyHaloWL *hwl, NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo, const gdouble R_min, const gdouble z_s, gdouble *kappa, gdouble *g_t)
{
  const gdouble z_cl = nc_wl_cluster_lens_get_z_cl (hwl->lens);

  g_mutex_lock (&hwl->bound_lock);

  if (!hwl->c_box)
  {
    *kappa = nc_wl_surface_mass_density_convergence (smd, hwl->dp_box, cosmo, R_min, z_s, z_cl, z_cl);
    *g_t   = fabs (nc_wl_surface_mass_density_reduced_shear (smd, hwl->dp_box, cosmo, R_min, z_s, z_cl, z_cl));
  }
  else
  {
    *kappa = -GSL_POSINF;

    if (!_nc_galaxy_halo_wl_corner_g_t (hwl, smd, cosmo, z_cl, R_min, z_s, g_t, kappa))
    {
      *kappa = GSL_NAN;
      *g_t   = GSL_NAN;
    }
  }

  g_mutex_unlock (&hwl->bound_lock);
}

/**
 * nc_galaxy_halo_wl_set_shape_bound:
 * @hwl: a #NcGalaxyHaloWL, with its bound prepared
 * @smd: (nullable): the #NcWLSurfaceMassDensity of the last nc_galaxy_halo_wl_prepare_bound()
 * @cosmo: the #NcHICosmo of the last nc_galaxy_halo_wl_prepare_bound()
 * @data: the #NcGalaxyHaloWLData of the galaxy, with its radius set
 * @z_max: the upper end of the redshift support of the galaxy
 * @s_data: the #NcGalaxyShapeFactorData of the galaxy
 *
 * Sets the bound fields of @s_data: whether the corner of the prior box is set and the
 * bound holds, the cluster redshift, the smallest radius of the galaxy over the prior of
 * the center, and the convergence and |g_t| at the corner (R_min, @z_max) when the corner
 * is set, R_min > 0 and @z_max is finite and behind the lens (NaN otherwise).
 */
void
nc_galaxy_halo_wl_set_shape_bound (NcGalaxyHaloWL *hwl, NcWLSurfaceMassDensity *smd, NcHICosmo *cosmo, NcGalaxyHaloWLData *data, const gdouble z_max, NcGalaxyShapeFactorData *s_data)
{
  const gdouble z_cl = nc_wl_cluster_lens_get_z_cl (hwl->lens);

  s_data->bound_ready = hwl->bound_ready;
  s_data->bound_holds = hwl->bound_holds;
  s_data->bound_z_cl  = z_cl;
  s_data->bound_R_min = nc_galaxy_halo_wl_get_radius_min (hwl, data);
  s_data->bound_kappa = GSL_NAN;
  s_data->bound_g_t   = GSL_NAN;

  if (s_data->bound_ready && (s_data->bound_R_min > 0.0) && gsl_finite (z_max) && (z_max > z_cl))
    nc_galaxy_halo_wl_eval_bound (hwl, smd, cosmo, s_data->bound_R_min, z_max, &s_data->bound_kappa, &s_data->bound_g_t);
}

