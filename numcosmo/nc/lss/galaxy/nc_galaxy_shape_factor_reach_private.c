/***************************************************************************
 *            nc_galaxy_shape_factor_reach_private.c
 *
 *  Mon Nov 16 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor_reach_private.c
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_factor_reach_private.h"
#include "nc/lss/galaxy/nc_galaxy_shape_pop.h"
#include "nc/lss/galaxy/nc_galaxy_position_factor.h"
#include "nc/lss/halo/nc_halo_position.h"
#include "nc/lss/halo/nc_halo_mass_summary.h"
#include "ncm/core/ncm_serialize.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

void
_nc_galaxy_shape_factor_reach_init (NcGalaxyShapeFactorReach *reach, const gboolean sample_c)
{
  memset (reach, 0, sizeof (NcGalaxyShapeFactorReach));
  reach->sample_c = sample_c;
  reach->gen      = 1;
  g_mutex_init (&reach->lock);
}

void
_nc_galaxy_shape_factor_reach_clear (NcGalaxyShapeFactorReach *reach)
{
  nc_hicosmo_clear (&reach->cosmo);
  nc_wl_surface_mass_density_clear (&reach->smd);
  nc_halo_density_profile_clear (&reach->dp_box);
  g_mutex_clear (&reach->lock);
}

static guint64
_nc_galaxy_shape_factor_reach_fnv (guint64 h, gconstpointer p, gsize n)
{
  const guint8 *c = p;
  gsize i;

  for (i = 0; i < n; i++)
  {
    h ^= c[i];
    h *= 1099511628211ULL;
  }

  return h;
}

/* The key of a model of the box: its type, and per parameter its free flag and its bounds or its value */
static guint64
_nc_galaxy_shape_factor_reach_model_key (guint64 h, NcmModel *model)
{
  const guint len = ncm_model_len (model);
  guint pid;

  h = _nc_galaxy_shape_factor_reach_fnv (h, G_OBJECT_TYPE_NAME (model), strlen (G_OBJECT_TYPE_NAME (model)));

  for (pid = 0; pid < len; pid++)
  {
    const gboolean is_free = ncm_model_param_get_ftype (model, pid) == NCM_PARAM_TYPE_FREE;
    gdouble v[2];

    v[0] = is_free ? ncm_model_param_get_lower_bound (model, pid) : ncm_model_param_get (model, pid);
    v[1] = is_free ? ncm_model_param_get_upper_bound (model, pid) : 0.0;

    h = _nc_galaxy_shape_factor_reach_fnv (h, &is_free, sizeof (is_free));
    h = _nc_galaxy_shape_factor_reach_fnv (h, v, sizeof (v));
  }

  return h;
}

/*
 * Whether the bound holds for the free parameters of mset: only log10MDelta, ra, dec, the
 * parameters of the shape population and, with sample_c, cDelta may be free. Read from the
 * models, not from the map of free parameters, which the caller may not have prepared.
 */
static gboolean
_nc_galaxy_shape_factor_reach_holds (NcGalaxyShapeFactorReach *reach, NcmMSet *mset)
{
  const guint n_models = ncm_mset_nmodels (mset);
  guint i, pid;

  for (i = 0; i < n_models; i++)
  {
    NcmModel *model = ncm_mset_peek (mset, ncm_mset_get_mid_array_pos (mset, i));
    const guint len = ncm_model_len (model);

    if (NC_IS_GALAXY_SHAPE_POP (model))
      continue;

    for (pid = 0; pid < len; pid++)
    {
      if (ncm_model_param_get_ftype (model, pid) != NCM_PARAM_TYPE_FREE)
        continue;

      if (NC_IS_HALO_POSITION (model) && ((pid == NC_HALO_POSITION_RA) || (pid == NC_HALO_POSITION_DEC)))
        continue;

      if (NC_IS_HALO_MASS_SUMMARY (model))
      {
        const gchar *name = ncm_model_param_name (model, pid);

        if (g_strcmp0 (name, "log10MDelta") == 0)
          continue;

        if (reach->sample_c && (g_strcmp0 (name, "cDelta") == 0))
          continue;
      }

      return FALSE;
    }
  }

  return TRUE;
}

/* Copies the parameters of src to dest, of the same type */
static void
_nc_galaxy_shape_factor_reach_copy_params (NcmModel *dest, NcmModel *src)
{
  const guint len = ncm_model_len (src);
  guint pid;

  for (pid = 0; pid < len; pid++)
    ncm_model_param_set (dest, pid, ncm_model_param_get (src, pid));
}

/*
 * Takes the box of the fit from mset (see NcGalaxyShapeFactorReach) and changes gen when the
 * box or the fixed models changed since the last call. The profile copy is made again when
 * the type of the profile or of its mass summary changes, and receives all the parameters of
 * the profile of the fit, with the mass at its upper bound.
 */
void
_nc_galaxy_shape_factor_reach_prepare (NcGalaxyShapeFactorReach *reach, NcmMSet *mset)
{
  NcHICosmo *cosmo            = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcHaloPosition *hp          = NC_HALO_POSITION (ncm_mset_peek (mset, nc_halo_position_id ()));
  NcWLSurfaceMassDensity *smd = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcHaloDensityProfile *dp    = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (mset, nc_halo_density_profile_id ()));
  guint64 key                 = 14695981039346656037ULL;

  reach->models_ready = FALSE;
  reach->holds        = FALSE;

  if ((cosmo == NULL) || (hp == NULL) || (smd == NULL) || (dp == NULL))
    return;

  {
    NcmModel *hms_src  = NCM_MODEL (nc_halo_density_profile_peek_mass_summary (dp));
    NcmModel *hp_model = NCM_MODEL (hp);
    GError *error      = NULL;
    NcmModel *hms_box;
    guint m_idx;

    if (!ncm_model_param_index_from_name (hms_src, "log10MDelta", &m_idx, &error))
    {
      /* LCOV_EXCL_START: every NcHaloMassSummary in the library has log10MDelta */
      g_clear_error (&error);

      return;
      /* LCOV_EXCL_STOP */
    }

    if ((reach->dp_box == NULL) || (G_OBJECT_TYPE (reach->dp_box) != G_OBJECT_TYPE (dp)) ||
        (G_OBJECT_TYPE (nc_halo_density_profile_peek_mass_summary (reach->dp_box)) != G_OBJECT_TYPE (hms_src)))
    {
      NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

      nc_halo_density_profile_clear (&reach->dp_box);
      reach->dp_box = NC_HALO_DENSITY_PROFILE (ncm_model_dup (NCM_MODEL (dp), ser));
      ncm_serialize_free (ser);
    }

    hms_box = NCM_MODEL (nc_halo_density_profile_peek_mass_summary (reach->dp_box));
    _nc_galaxy_shape_factor_reach_copy_params (NCM_MODEL (reach->dp_box), NCM_MODEL (dp));
    _nc_galaxy_shape_factor_reach_copy_params (hms_box, hms_src);
    ncm_model_param_set (hms_box, m_idx, ncm_model_param_get_upper_bound (hms_src, m_idx));

    reach->has_c = ncm_model_param_index_from_name (hms_src, "cDelta", &reach->c_idx, &error);
    g_clear_error (&error);

    if (reach->has_c)
    {
      const gboolean c_free = reach->sample_c && (ncm_model_param_get_ftype (hms_src, reach->c_idx) == NCM_PARAM_TYPE_FREE);

      reach->c_lo = c_free ? ncm_model_param_get_lower_bound (hms_src, reach->c_idx) : ncm_model_param_get (hms_src, reach->c_idx);
      reach->c_hi = c_free ? ncm_model_param_get_upper_bound (hms_src, reach->c_idx) : reach->c_lo;
    }

    nc_wl_surface_mass_density_prepare_if_needed (smd, cosmo);
    nc_halo_position_prepare_if_needed (hp, cosmo);

    if (reach->cosmo != cosmo)
    {
      nc_hicosmo_clear (&reach->cosmo);
      reach->cosmo = nc_hicosmo_ref (cosmo);
    }

    if (reach->smd != smd)
    {
      nc_wl_surface_mass_density_clear (&reach->smd);
      reach->smd = nc_wl_surface_mass_density_ref (smd);
    }

    reach->z_cl      = nc_halo_position_get_redshift (hp);
    reach->prefactor = nc_halo_position_projected_radius_prefactor (hp, cosmo);

    /* The box of the halo center: its center and the largest distance from it to a corner */
    {
      const gboolean ra_free  = ncm_model_param_get_ftype (hp_model, NC_HALO_POSITION_RA) == NCM_PARAM_TYPE_FREE;
      const gboolean dec_free = ncm_model_param_get_ftype (hp_model, NC_HALO_POSITION_DEC) == NCM_PARAM_TYPE_FREE;
      const gdouble ra_lo     = ra_free ? ncm_model_param_get_lower_bound (hp_model, NC_HALO_POSITION_RA) : ncm_model_param_get (hp_model, NC_HALO_POSITION_RA);
      const gdouble ra_hi     = ra_free ? ncm_model_param_get_upper_bound (hp_model, NC_HALO_POSITION_RA) : ra_lo;
      const gdouble dec_lo    = dec_free ? ncm_model_param_get_lower_bound (hp_model, NC_HALO_POSITION_DEC) : ncm_model_param_get (hp_model, NC_HALO_POSITION_DEC);
      const gdouble dec_hi    = dec_free ? ncm_model_param_get_upper_bound (hp_model, NC_HALO_POSITION_DEC) : dec_lo;
      const gdouble ra_c[2]   = {ra_lo, ra_hi};
      const gdouble dec_c[2]  = {dec_lo, dec_hi};
      guint i, j;

      ncm_trivec_set_astro_ra_dec (&reach->centre, 1.0, 0.5 * (ra_lo + ra_hi), 0.5 * (dec_lo + dec_hi));
      reach->delta = 0.0;

      for (i = 0; i < 2; i++)
      {
        for (j = 0; j < 2; j++)
        {
          NcmTriVec corner;
          gdouble d2 = 0.0;
          guint k;

          ncm_trivec_set_astro_ra_dec (&corner, 1.0, ra_c[i], dec_c[j]);

          for (k = 0; k < 3; k++)
            d2 += gsl_pow_2 (corner.c[k] - reach->centre.c[k]);

          reach->delta = MAX (reach->delta, reach->prefactor * sin (2.0 * asin (MIN (1.0, 0.5 * sqrt (d2)))));
        }
      }
    }

    /* Everything the bound depends on: the box and the fixed models */
    key = _nc_galaxy_shape_factor_reach_fnv (key, &reach->sample_c, sizeof (reach->sample_c));
    key = _nc_galaxy_shape_factor_reach_model_key (key, NCM_MODEL (cosmo));
    key = _nc_galaxy_shape_factor_reach_model_key (key, NCM_MODEL (smd));
    key = _nc_galaxy_shape_factor_reach_model_key (key, hp_model);
    key = _nc_galaxy_shape_factor_reach_model_key (key, NCM_MODEL (dp));
    key = _nc_galaxy_shape_factor_reach_model_key (key, hms_src);

    if (key != reach->key)
    {
      reach->key = key;
      reach->gen++;
    }

    reach->models_ready = TRUE;
    reach->holds        = _nc_galaxy_shape_factor_reach_holds (reach, mset);
  }
}

/*
 * The largest |g| the galaxy of data can reach over the box, for the redshift support of the
 * galaxy up to z_max: the reduced shear of the profile at the upper bound of the mass, at the
 * smallest projected radius the box of the center allows and at z_max, with the
 * multiplicative and additive bias of the galaxy, |1 + m| |g_t| + |c| (the triangle
 * inequality on the calibrated shear (1 + m) g_t e^{2i phi} + c). At fixed radius, source
 * redshift and concentration both the convergence and the shear grow with the mass and the
 * source redshift, so with a convergence below one the bound is at that corner. The smallest
 * radius is R_c - delta, R_c the projected distance of the galaxy from the center of the box
 * and delta the largest one from that center to a corner (the triangle inequality on the
 * angle, sin being subadditive below pi / 2). The shear is not monotonic in the
 * concentration: a free cDelta (sample_c) is sampled at NC_GALAXY_SHAPE_FACTOR_REACH_NC
 * log-spaced points of its box and the largest value taken times
 * NC_GALAXY_SHAPE_FACTOR_REACH_C_MARGIN.
 *
 * Returns zero for a galaxy with no support behind the halo, and G_MAXDOUBLE when there is
 * nothing to bound: a radius the box of the center takes to the halo, a support with no upper
 * end, or a convergence of one or more at the corner (|g| crossed one inside the box). The
 * value is kept in memo and computed again only when gen or the values of the galaxy change.
 */
gdouble
_nc_galaxy_shape_factor_reach_eval (NcGalaxyShapeFactorReach *reach, NcGalaxyShapeFactorReachMemo *memo, NcGalaxyShapeFactorData *data, const gdouble z_max)
{
  const gdouble ra  = data->pos_data->ra;
  const gdouble dec = data->pos_data->dec;
  gdouble value;

  if ((memo->gen == reach->gen) && (memo->z_max == z_max) && (memo->ra == ra) && (memo->dec == dec) &&
      (memo->m == data->m) && (memo->c1 == data->c1) && (memo->c2 == data->c2))
    return memo->value;

  if (z_max <= reach->z_cl)
  {
    value = 0.0;
  }
  else
  {
    NcmTriVec v;
    gdouble d2 = 0.0, R_min;
    guint k;

    ncm_trivec_set_astro_ra_dec (&v, 1.0, ra, dec);

    for (k = 0; k < 3; k++)
      d2 += gsl_pow_2 (v.c[k] - reach->centre.c[k]);

    R_min = reach->prefactor * sin (2.0 * asin (MIN (1.0, 0.5 * sqrt (d2)))) - reach->delta;

    if ((R_min <= 0.0) || !gsl_finite (z_max))
    {
      value = G_MAXDOUBLE;
    }
    else
    {
      const guint n_c = (reach->has_c && (reach->c_hi > reach->c_lo)) ? NC_GALAXY_SHAPE_FACTOR_REACH_NC : 1;
      NcmModel *hms   = NCM_MODEL (nc_halo_density_profile_peek_mass_summary (reach->dp_box));
      gdouble g_t_max = 0.0;
      guint i;

      g_mutex_lock (&reach->lock);

      for (i = 0; i < n_c; i++)
      {
        gdouble kappa, g_t;

        if (n_c > 1)
          ncm_model_param_set (hms, reach->c_idx, reach->c_lo * pow (reach->c_hi / reach->c_lo, i / (n_c - 1.0)));

        kappa = nc_wl_surface_mass_density_convergence (reach->smd, reach->dp_box, reach->cosmo, R_min, z_max, reach->z_cl, reach->z_cl);
        g_t   = fabs (nc_wl_surface_mass_density_reduced_shear (reach->smd, reach->dp_box, reach->cosmo, R_min, z_max, reach->z_cl, reach->z_cl));

        if (!gsl_finite (kappa) || !gsl_finite (g_t) || (kappa >= 1.0))
        {
          g_t_max = G_MAXDOUBLE;
          break;
        }

        g_t_max = MAX (g_t_max, g_t);
      }

      g_mutex_unlock (&reach->lock);

      if (g_t_max == G_MAXDOUBLE)
        value = G_MAXDOUBLE;
      else
        value = fabs (1.0 + data->m) * g_t_max * ((n_c > 1) ? NC_GALAXY_SHAPE_FACTOR_REACH_C_MARGIN : 1.0) + hypot (data->c1, data->c2);
    }
  }

  memo->gen   = reach->gen;
  memo->z_max = z_max;
  memo->ra    = ra;
  memo->dec   = dec;
  memo->m     = data->m;
  memo->c1    = data->c1;
  memo->c2    = data->c2;
  memo->value = value;

  return value;
}

