/***************************************************************************
 *            nc_wl_cluster_lens.c
 *
 *  Fri November 27 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_wl_cluster_lens.c
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
 * NcWLClusterLens:
 *
 * The weak-lensing signal of a cluster at its redshift.
 *
 * A cluster at redshift $z_\mathrm{cl}$ with a halo density profile shears the images
 * of the sources behind it. This calculator computes the reduced shear
 * $$
 * g(R, z_s) = \frac{\left[\overline{\Sigma}(<R) - \Sigma(R)\right]/\Sigma_\mathrm{cr}(z_s)}
 * {1 - \Sigma(R)/\Sigma_\mathrm{cr}(z_s)},
 * $$
 * where $R$ is the projected radius, $z_s$ the source redshift, $\Sigma$ and
 * $\overline{\Sigma}$ the surface density and its mean inside $R$, and
 * $\Sigma_\mathrm{cr}$ the critical surface density; $g = 0$ for $z_s < z_\mathrm{cl}$.
 * See <a href="../../theory/nc/lss/wl/wl_surface_mass_density.html">Weak-Lensing
 * Surface Mass Density</a>.
 *
 * It holds no models: nc_wl_cluster_lens_prepare() receives the cosmology, the halo
 * position and the density profile, and keeps what the evaluations need until the next
 * prepare.
 *
 * A caller that evaluates $\Sigma_\mathrm{cr}^{-1}(z_s)$ at many source redshifts builds a
 * table of it on $[z_\mathrm{cl}, z_\mathrm{max}]$ with
 * nc_wl_cluster_lens_prepare_sigma_crit_inv(), where it gives $z_\mathrm{max}$, and
 * evaluates the #NcmSpectralPanels of nc_wl_cluster_lens_peek_sigma_crit_inv().
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/wl/nc_wl_cluster_lens.h"
#include "ncm/core/ncm_c.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/model/ncm_model_ctrl.h"
#include "ncm/algebra/ncm_spectral_panels_builder.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcWLClusterLens
{
  /*< private >*/
  GObject parent_instance;
  NcDistance *dist;
  gdouble z_cl;
  gint k;
  gdouble sqrt_Omega_k0;
  gdouble dl;
  gdouble sc_over_Dl;
  gdouble r_s;
  gdouble rho_s;
  NcmModelCtrl *ctrl_cosmo;
  NcmModelCtrl *ctrl_hp;
  NcmModelCtrl *ctrl_dp;
  guint64 geometry_hash; /* counts the prepares that changed the geometry */
  guint64 profile_hash;  /* counts the prepares that changed the profile */
  gdouble sc_inv_reltol;
  NcmSpectralPanels *sc_inv; /* 1/Sigma_cr (z_s) on [z_cl, sc_inv_z_max] */
  NcmSpectralPanelsBuilder *builder;
  gdouble sc_inv_z_max;         /* NaN when sc_inv was never built */
  guint64 sc_inv_geometry_hash; /* geometry_hash when sc_inv was built */
};

/*
 * Coefficients of every panel of the 1/Sigma_cr table. N = 8 was the fastest of 5 to 9 in the
 * cluster WL likelihood (NumCosmoDevNotes benches/wl_lens_redesign); the panel count follows
 * from the tolerance.
 */
#define NC_WL_CLUSTER_LENS_SIGMA_CRIT_INV_N (8)

enum
{
  PROP_0,
  PROP_DISTANCE,
  PROP_SIGMA_CRIT_INV_RELTOL,
};

G_DEFINE_TYPE (NcWLClusterLens, nc_wl_cluster_lens, G_TYPE_OBJECT)

static void _nc_wl_cluster_lens_source_distances (NcWLClusterLens *lens, NcHICosmo *cosmo, const gdouble z_s, gdouble *Ds, gdouble *Dls);

static void
nc_wl_cluster_lens_init (NcWLClusterLens *lens)
{
  lens->dist                 = NULL;
  lens->z_cl                 = GSL_NAN;
  lens->k                    = 0;
  lens->sqrt_Omega_k0        = GSL_NAN;
  lens->dl                   = GSL_NAN;
  lens->sc_over_Dl           = GSL_NAN;
  lens->r_s                  = GSL_NAN;
  lens->rho_s                = GSL_NAN;
  lens->ctrl_cosmo           = ncm_model_ctrl_new (NULL);
  lens->ctrl_hp              = ncm_model_ctrl_new (NULL);
  lens->ctrl_dp              = ncm_model_ctrl_new (NULL);
  lens->geometry_hash        = 0;
  lens->profile_hash         = 0;
  lens->sc_inv_reltol        = 0.0;
  lens->sc_inv               = ncm_spectral_panels_new (1);
  lens->builder              = ncm_spectral_panels_builder_new ();
  lens->sc_inv_z_max         = GSL_NAN;
  lens->sc_inv_geometry_hash = 0;
}

static void
_nc_wl_cluster_lens_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcWLClusterLens *lens = NC_WL_CLUSTER_LENS (object);

  switch (prop_id)
  {
    case PROP_DISTANCE:
      nc_distance_clear (&lens->dist);
      lens->dist = g_value_dup_object (value);
      break;
    case PROP_SIGMA_CRIT_INV_RELTOL:
      nc_wl_cluster_lens_set_sigma_crit_inv_reltol (lens, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_wl_cluster_lens_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcWLClusterLens *lens = NC_WL_CLUSTER_LENS (object);

  switch (prop_id)
  {
    case PROP_DISTANCE:
      g_value_set_object (value, lens->dist);
      break;
    case PROP_SIGMA_CRIT_INV_RELTOL:
      g_value_set_double (value, lens->sc_inv_reltol);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_wl_cluster_lens_dispose (GObject *object)
{
  NcWLClusterLens *lens = NC_WL_CLUSTER_LENS (object);

  nc_distance_clear (&lens->dist);
  ncm_model_ctrl_clear (&lens->ctrl_cosmo);
  ncm_model_ctrl_clear (&lens->ctrl_hp);
  ncm_model_ctrl_clear (&lens->ctrl_dp);
  ncm_spectral_panels_clear (&lens->sc_inv);
  ncm_spectral_panels_builder_clear (&lens->builder);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_wl_cluster_lens_parent_class)->dispose (object);
}

static void
nc_wl_cluster_lens_class_init (NcWLClusterLensClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_wl_cluster_lens_set_property;
  object_class->get_property = &_nc_wl_cluster_lens_get_property;
  object_class->dispose      = &_nc_wl_cluster_lens_dispose;

  /**
   * NcWLClusterLens:distance:
   *
   * The distance object of the critical surface density.
   */
  g_object_class_install_property (object_class,
                                   PROP_DISTANCE,
                                   g_param_spec_object ("distance",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcWLClusterLens:sigma-crit-inv-reltol:
   *
   * The relative tolerance of the $\Sigma_\mathrm{cr}^{-1}(z_s)$ table of
   * nc_wl_cluster_lens_prepare_sigma_crit_inv(): on every panel the table differs from
   * $\Sigma_\mathrm{cr}^{-1}$ by at most this value times the larger of its largest value on
   * the panel and $\Sigma_\mathrm{cr}^{-1}(z_\mathrm{max})$, at the points checked by
   * ncm_spectral_panels_builder_build_uniform().
   */
  g_object_class_install_property (object_class,
                                   PROP_SIGMA_CRIT_INV_RELTOL,
                                   g_param_spec_double ("sigma-crit-inv-reltol",
                                                        NULL,
                                                        "Relative tolerance of the inverse critical density table",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-14,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_wl_cluster_lens_new:
 * @dist: a #NcDistance
 *
 * Creates a new #NcWLClusterLens computing its distances with @dist.
 *
 * Returns: (transfer full): a new #NcWLClusterLens.
 */
NcWLClusterLens *
nc_wl_cluster_lens_new (NcDistance *dist)
{
  return g_object_new (NC_TYPE_WL_CLUSTER_LENS,
                       "distance", dist,
                       NULL);
}

/**
 * nc_wl_cluster_lens_ref:
 * @lens: a #NcWLClusterLens
 *
 * Increases the reference count of @lens by one.
 *
 * Returns: (transfer full): @lens.
 */
NcWLClusterLens *
nc_wl_cluster_lens_ref (NcWLClusterLens *lens)
{
  return g_object_ref (lens);
}

/**
 * nc_wl_cluster_lens_free:
 * @lens: a #NcWLClusterLens
 *
 * Decreases the reference count of @lens by one.
 */
void
nc_wl_cluster_lens_free (NcWLClusterLens *lens)
{
  g_object_unref (lens);
}

/**
 * nc_wl_cluster_lens_clear:
 * @lens: a #NcWLClusterLens
 *
 * Decreases the reference count of *@lens by one and sets the pointer *@lens to %NULL.
 */
void
nc_wl_cluster_lens_clear (NcWLClusterLens **lens)
{
  g_clear_object (lens);
}

/**
 * nc_wl_cluster_lens_peek_distance:
 * @lens: a #NcWLClusterLens
 *
 * Returns: (transfer none): the #NcDistance of @lens.
 */
NcDistance *
nc_wl_cluster_lens_peek_distance (NcWLClusterLens *lens)
{
  return lens->dist;
}

/**
 * nc_wl_cluster_lens_prepare:
 * @lens: a #NcWLClusterLens
 * @cosmo: a #NcHICosmo
 * @hp: a #NcHaloPosition, whose redshift is $z_\mathrm{cl}$
 * @dp: a #NcHaloDensityProfile
 *
 * Prepares @lens for the cluster of @hp and @dp in @cosmo: the lens-side distances of
 * $\Sigma_\mathrm{cr}$ when the cosmology or $z_\mathrm{cl}$ changed, and the scale
 * radius and density of @dp at $z_\mathrm{cl}$ when the cosmology, the profile, its mass
 * summary or the halo position changed. A model changed when it is another object or its
 * state changed (#NcmModelCtrl). Prepares the distance and @hp if needed. Call
 * nc_wl_cluster_lens_sigma() with the same @dp.
 */
void
nc_wl_cluster_lens_prepare (NcWLClusterLens *lens, NcHICosmo *cosmo, NcHaloPosition *hp, NcHaloDensityProfile *dp)
{
  const gdouble z_cl = nc_halo_position_get_redshift (hp);
  gboolean cosmo_changed, hp_changed, dp_changed, z_cl_changed;

  nc_distance_prepare_if_needed (lens->dist, cosmo);
  nc_halo_position_prepare_if_needed (hp, cosmo);

  /*
   * The halo position also changes with ra and dec: the geometry compares z_cl by value.
   * The ctrl of the profile also follows its mass summary, a submodel.
   */
  cosmo_changed = ncm_model_ctrl_update (lens->ctrl_cosmo, NCM_MODEL (cosmo));
  hp_changed    = ncm_model_ctrl_update (lens->ctrl_hp, NCM_MODEL (hp));
  dp_changed    = ncm_model_ctrl_update (lens->ctrl_dp, NCM_MODEL (dp));
  z_cl_changed  = !(z_cl == lens->z_cl);

  lens->z_cl = z_cl;

  if (cosmo_changed || z_cl_changed)
  {
    const gdouble a             = ncm_c_c2 () / (4.0 * M_PI * ncm_c_G_mass_solar ()) * ncm_c_Mpc () / nc_hicosmo_RH_Mpc (cosmo);
    const gdouble Omega_k0      = nc_hicosmo_Omega_k0 (cosmo);
    const gdouble sqrt_Omega_k0 = sqrt (fabs (Omega_k0));
    const gint k                = fabs (Omega_k0) < NCM_ZERO_LIMIT ? 0 : (Omega_k0 > 0.0 ? -1 : 1);
    const gdouble dl            = nc_distance_comoving (lens->dist, cosmo, z_cl);
    gdouble Dl;

    switch (k)
    {
      case -1:
        Dl = sinh (sqrt_Omega_k0 * dl) / ((1.0 + z_cl) * sqrt_Omega_k0);
        break;
      case 0:
        Dl = dl / (1.0 + z_cl);
        break;
      case 1:
        Dl = sin (sqrt_Omega_k0 * dl) / ((1.0 + z_cl) * sqrt_Omega_k0);
        break;
      default:                   /* LCOV_EXCL_LINE */
        g_assert_not_reached (); /* LCOV_EXCL_LINE */
        break;                   /* LCOV_EXCL_LINE */
    }

    lens->k             = k;
    lens->sqrt_Omega_k0 = sqrt_Omega_k0;
    lens->dl            = dl;
    lens->sc_over_Dl    = a / Dl;
    lens->geometry_hash++;
  }

  if (cosmo_changed || hp_changed || dp_changed)
  {
    nc_halo_density_profile_r_s_rho_s (dp, cosmo, z_cl, &lens->r_s, &lens->rho_s);
    lens->profile_hash++;
  }
}

/**
 * nc_wl_cluster_lens_get_z_cl:
 * @lens: a #NcWLClusterLens
 *
 * Returns: the cluster redshift $z_\mathrm{cl}$ of the last prepare.
 */
gdouble
nc_wl_cluster_lens_get_z_cl (NcWLClusterLens *lens)
{
  return lens->z_cl;
}

/**
 * nc_wl_cluster_lens_get_geometry_hash:
 * @lens: a #NcWLClusterLens
 *
 * Returns: the number of prepares that changed the cosmology or $z_\mathrm{cl}$, so
 *   $\Sigma_\mathrm{cr}(z_s)$; a different value means a change.
 */
guint64
nc_wl_cluster_lens_get_geometry_hash (NcWLClusterLens *lens)
{
  return lens->geometry_hash;
}

/**
 * nc_wl_cluster_lens_get_profile_hash:
 * @lens: a #NcWLClusterLens
 *
 * Returns: the number of prepares that changed the cosmology, the profile, its mass
 *   summary or the halo position, so $\Sigma(R)$ and $\overline{\Sigma}(<R)$; a different
 *   value means a change.
 */
guint64
nc_wl_cluster_lens_get_profile_hash (NcWLClusterLens *lens)
{
  return lens->profile_hash;
}

/**
 * nc_wl_cluster_lens_set_sigma_crit_inv_reltol:
 * @lens: a #NcWLClusterLens
 * @reltol: the relative tolerance
 *
 * Sets #NcWLClusterLens:sigma-crit-inv-reltol; the next
 * nc_wl_cluster_lens_prepare_sigma_crit_inv() builds the table again.
 */
void
nc_wl_cluster_lens_set_sigma_crit_inv_reltol (NcWLClusterLens *lens, const gdouble reltol)
{
  lens->sc_inv_reltol = reltol;
  lens->sc_inv_z_max  = GSL_NAN;
}

/**
 * nc_wl_cluster_lens_get_sigma_crit_inv_reltol:
 * @lens: a #NcWLClusterLens
 *
 * Returns: #NcWLClusterLens:sigma-crit-inv-reltol.
 */
gdouble
nc_wl_cluster_lens_get_sigma_crit_inv_reltol (NcWLClusterLens *lens)
{
  return lens->sc_inv_reltol;
}

/* The function of the table, 1/Sigma_cr (z_s), and 0 at or in front of the lens */
static void
_nc_wl_cluster_lens_sigma_crit_inv_F (gpointer user_data, guint i, gdouble z_s, NcmVector *fval)
{
  gpointer *args = user_data;

  ncm_vector_set (fval, 0, nc_wl_cluster_lens_sigma_crit_inv (args[0], args[1], z_s));
}

/**
 * nc_wl_cluster_lens_prepare_sigma_crit_inv:
 * @lens: a #NcWLClusterLens
 * @cosmo: the #NcHICosmo of the last nc_wl_cluster_lens_prepare()
 * @z_max: the largest source redshift of the table, above $z_\mathrm{cl}$
 *
 * Builds the table of $\Sigma_\mathrm{cr}^{-1}(z_s)$ (nc_wl_cluster_lens_sigma_crit_inv())
 * on $[z_\mathrm{cl}, z_\mathrm{max}]$ that nc_wl_cluster_lens_peek_sigma_crit_inv()
 * returns: uniform Chebyshev panels (ncm_spectral_panels_builder_build_uniform()) with the
 * relative tolerance #NcWLClusterLens:sigma-crit-inv-reltol and the absolute tolerance
 * that tolerance times $\Sigma_\mathrm{cr}^{-1}(z_\mathrm{max})$. The table is built again
 * only when the geometry (nc_wl_cluster_lens_get_geometry_hash()), @z_max or the tolerance
 * changed since the last build. The caller gives @z_max, the largest source redshift it
 * evaluates, and evaluates the table only on $[z_\mathrm{cl}, z_\mathrm{max}]$. An @z_max
 * not above $z_\mathrm{cl}$ is an error.
 */
void
nc_wl_cluster_lens_prepare_sigma_crit_inv (NcWLClusterLens *lens, NcHICosmo *cosmo, const gdouble z_max)
{
  if ((lens->sc_inv_geometry_hash == lens->geometry_hash) && (z_max == lens->sc_inv_z_max))
    return;

  if (!(z_max > lens->z_cl))
    g_error ("nc_wl_cluster_lens_prepare_sigma_crit_inv: z_max = %.17g is not above z_cl = %.17g.", z_max, lens->z_cl);

  {
    gpointer args[2]    = {lens, cosmo};
    const gdouble scale = nc_wl_cluster_lens_sigma_crit_inv (lens, cosmo, z_max);

    ncm_spectral_panels_builder_build_uniform (lens->builder, _nc_wl_cluster_lens_sigma_crit_inv_F, args,
                                               lens->z_cl, z_max, NC_WL_CLUSTER_LENS_SIGMA_CRIT_INV_N,
                                               lens->sc_inv_reltol, lens->sc_inv_reltol * scale, lens->sc_inv);
  }

  lens->sc_inv_z_max         = z_max;
  lens->sc_inv_geometry_hash = lens->geometry_hash;
}

/**
 * nc_wl_cluster_lens_peek_sigma_crit_inv:
 * @lens: a #NcWLClusterLens
 *
 * Gets the table of $\Sigma_\mathrm{cr}^{-1}(z_s)$ of the last
 * nc_wl_cluster_lens_prepare_sigma_crit_inv(), row 0 of a uniform #NcmSpectralPanels
 * evaluated by ncm_spectral_panels_eval_row(). It is valid on $[z_\mathrm{cl},
 * z_\mathrm{max}]$ until the next prepare.
 *
 * Returns: (transfer none): the table.
 */
NcmSpectralPanels *
nc_wl_cluster_lens_peek_sigma_crit_inv (NcWLClusterLens *lens)
{
  return lens->sc_inv;
}

/**
 * nc_wl_cluster_lens_sigma_crit:
 * @lens: a #NcWLClusterLens
 * @cosmo: the #NcHICosmo of the last prepare
 * @z_s: the source redshift, at least $z_\mathrm{cl}$
 *
 * Computes the critical surface density $\Sigma_\mathrm{cr}(z_s)$ of the last prepare,
 * in $M_\odot/\mathrm{Mpc}^2$.
 *
 * Returns: $\Sigma_\mathrm{cr}(z_s)$.
 */
gdouble
nc_wl_cluster_lens_sigma_crit (NcWLClusterLens *lens, NcHICosmo *cosmo, const gdouble z_s)
{
  gdouble Ds, Dls;

  _nc_wl_cluster_lens_source_distances (lens, cosmo, z_s, &Ds, &Dls);

  return lens->sc_over_Dl * Ds / Dls;
}

/**
 * nc_wl_cluster_lens_sigma_crit_inv:
 * @lens: a #NcWLClusterLens
 * @cosmo: the #NcHICosmo of the last prepare
 * @z_s: the source redshift, at least $z_\mathrm{cl}$
 *
 * Computes the inverse critical surface density $\Sigma_\mathrm{cr}^{-1}(z_s)$ of the last
 * prepare, in $\mathrm{Mpc}^2/M_\odot$; it is 0 at $z_s = z_\mathrm{cl}$.
 *
 * Returns: $\Sigma_\mathrm{cr}^{-1}(z_s)$.
 */
gdouble
nc_wl_cluster_lens_sigma_crit_inv (NcWLClusterLens *lens, NcHICosmo *cosmo, const gdouble z_s)
{
  gdouble Ds, Dls;

  _nc_wl_cluster_lens_source_distances (lens, cosmo, z_s, &Ds, &Dls);

  return Dls / (lens->sc_over_Dl * Ds);
}

/* The angular diameter distances of the source, D_s, and between lens and source, D_ls */
static void
_nc_wl_cluster_lens_source_distances (NcWLClusterLens *lens, NcHICosmo *cosmo, const gdouble z_s, gdouble *Ds, gdouble *Dls)
{
  const gdouble x_i = 1.0 + z_s;
  const gdouble ds  = nc_distance_comoving (lens->dist, cosmo, z_s);

  switch (lens->k)
  {
    case -1:
      *Ds  = sinh (lens->sqrt_Omega_k0 * ds) / (x_i * lens->sqrt_Omega_k0);
      *Dls = sinh (lens->sqrt_Omega_k0 * (ds - lens->dl)) / (x_i * lens->sqrt_Omega_k0);
      break;
    case 0:
      *Ds  = ds / x_i;
      *Dls = (ds - lens->dl) / x_i;
      break;
    case 1:
      *Ds  = sin (lens->sqrt_Omega_k0 * ds) / (x_i * lens->sqrt_Omega_k0);
      *Dls = sin (lens->sqrt_Omega_k0 * (ds - lens->dl)) / (x_i * lens->sqrt_Omega_k0);
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
      break;                   /* LCOV_EXCL_LINE */
  }
}

/**
 * nc_wl_cluster_lens_sigma:
 * @lens: a #NcWLClusterLens
 * @dp: the #NcHaloDensityProfile of the last prepare
 * @R: the projected radius, in Mpc
 * @sigma: (out): the surface density $\Sigma(R)$
 * @mean_sigma: (out): the mean surface density $\overline{\Sigma}(<R)$
 *
 * Computes $\Sigma(R)$ and $\overline{\Sigma}(<R)$ of @dp at $z_\mathrm{cl}$, in
 * $M_\odot/\mathrm{Mpc}^2$, with the scale radius and density of the last prepare.
 * Concurrent calls are safe once @dp has been evaluated after its last parameter change;
 * a numerically integrated profile builds its splines on that first call.
 */
void
nc_wl_cluster_lens_sigma (NcWLClusterLens *lens, NcHaloDensityProfile *dp, const gdouble R, gdouble *sigma, gdouble *mean_sigma)
{
  const gdouble X = R / lens->r_s;

  *mean_sigma = (2.0 * nc_halo_density_profile_eval_dl_cyl_mass (dp, X) / (X * X)) * lens->rho_s * lens->r_s;
  *sigma      = (nc_halo_density_profile_eval_dl_2d_density (dp, X)) * lens->rho_s * lens->r_s;
}

/**
 * nc_wl_cluster_lens_reduced_shear:
 * @sigma: the surface density $\Sigma(R)$
 * @mean_sigma: the mean surface density $\overline{\Sigma}(<R)$
 * @sigma_crit_inv: the inverse critical surface density $\Sigma_\mathrm{cr}^{-1}(z_s)$
 *
 * Computes the tangential reduced shear
 * $g = (\overline{\Sigma} - \Sigma)\Sigma_\mathrm{cr}^{-1}/(1 - \Sigma\Sigma_\mathrm{cr}^{-1})$.
 *
 * Returns: $g$.
 */

/**
 * nc_wl_cluster_lens_reduced_shear_at:
 * @lens: a #NcWLClusterLens
 * @cosmo: the #NcHICosmo of the last prepare
 * @sigma: the surface density $\Sigma(R)$
 * @mean_sigma: the mean surface density $\overline{\Sigma}(<R)$
 * @z_s: the source redshift
 *
 * Computes the tangential reduced shear $g(R, z_s)$ from $\Sigma(R)$ and
 * $\overline{\Sigma}(<R)$ (nc_wl_cluster_lens_sigma()); 0 for $z_s < z_\mathrm{cl}$.
 *
 * Returns: $g(R, z_s)$.
 */
gdouble
nc_wl_cluster_lens_reduced_shear_at (NcWLClusterLens *lens, NcHICosmo *cosmo, const gdouble sigma, const gdouble mean_sigma, const gdouble z_s)
{
  if (z_s < lens->z_cl)
    return 0.0;

  {
    /* The division form of NcWLSurfaceMassDensity, whose values this equals bit for bit */
    const gdouble sigma_crit = nc_wl_cluster_lens_sigma_crit (lens, cosmo, z_s);

    return ((mean_sigma - sigma) / sigma_crit) / (1.0 - sigma / sigma_crit);
  }
}

