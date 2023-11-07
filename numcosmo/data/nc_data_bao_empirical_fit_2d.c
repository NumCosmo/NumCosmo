/***************************************************************************
 *            nc_data_bao_empirical_fit_2d.c
 *
 *  Fri September 1 14:39:23 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_data_bao_empirical_fit_2d.c
 * Copyright (C) 2017 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_data_bao_empirical_fit_2d
 * @title: NcDataBaoEmpiricalFit2d
 * @short_description: Baryon oscillation data -- $D_H / r_d$ and $D_t / r_d$ empirical likelihood.
 *
 * This object implements the BAO data when its likelihood function is provided,
 * e.g., [Bautista et al. (2017)][XBautista2017].
 *
 * The data is stored in a #NcDataBaoEmpiricalFit2d object. The data is stored in a
 * #NcmDataDist2d base class object, which is a subclass of #NcmData.
 * The data represents the likelihood function of the transverse distance $D_t$ and the
 * Hubble distance $D_H$ at the redshift $z$ divided by the sound horizon at the
 * last scattering surface $r_s$. The likelihood function is provided as a
 * #NcmSpline2d object.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_bao_empirical_fit_2d.h"
#include "nc_hicosmo.h"
#include "nc_distance.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_DH_RD_FIDUC,
  PROP_DT_RD_FIDUC,
  PROP_Z,
  PROP_M2LNP,
  PROP_DIST,
  PROP_SIZE
};

G_DEFINE_TYPE (NcDataBaoEmpiricalFit2d, nc_data_bao_empirical_fit_2d, NCM_TYPE_DATA_DIST2D);

static void
nc_data_bao_empirical_fit_2d_init (NcDataBaoEmpiricalFit2d *bao_ef)
{
  bao_ef->Dh_rd_fiduc = 0.0;
  bao_ef->Dt_rd_fiduc = 0.0;
  bao_ef->z           = 0.0;
  bao_ef->m2lnp       = NULL;
  bao_ef->p           = NULL;
  bao_ef->dist        = nc_distance_new (2.0);
}

static void
nc_data_bao_empirical_fit_2d_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_bao_empirical_fit_2d_parent_class)->constructed (object);
  {
    NcDataBaoEmpiricalFit2d *bao_ef = NC_DATA_BAO_EMPIRICAL_FIT_2D (object);

    ncm_stats_dist2d_prepare (bao_ef->p);

    ncm_data_set_init (NCM_DATA (bao_ef), TRUE);
  }
}

static void
nc_data_bao_empirical_fit_2d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataBaoEmpiricalFit2d *bao_ef = NC_DATA_BAO_EMPIRICAL_FIT_2D (object);

  g_return_if_fail (NC_IS_DATA_BAO_EMPIRICAL_FIT_2D (object));

  switch (prop_id)
  {
    case PROP_DH_RD_FIDUC:
      bao_ef->Dh_rd_fiduc = g_value_get_double (value);
      break;
    case PROP_DT_RD_FIDUC:
      bao_ef->Dt_rd_fiduc = g_value_get_double (value);
      break;
    case PROP_Z:
      bao_ef->z = g_value_get_double (value);
      break;
    case PROP_M2LNP:
      ncm_spline2d_clear (&bao_ef->m2lnp);
      ncm_stats_dist2d_clear (&bao_ef->p);
      bao_ef->m2lnp = g_value_dup_object (value);
      bao_ef->p     = NCM_STATS_DIST2D (ncm_stats_dist2d_spline_new (bao_ef->m2lnp));
      break;
    case PROP_DIST:
      nc_data_bao_empirical_fit_2d_set_dist (bao_ef, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_empirical_fit_2d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataBaoEmpiricalFit2d *bao_ef = NC_DATA_BAO_EMPIRICAL_FIT_2D (object);

  g_return_if_fail (NC_IS_DATA_BAO_EMPIRICAL_FIT_2D (object));

  switch (prop_id)
  {
    case PROP_DH_RD_FIDUC:
      g_value_set_double (value, bao_ef->Dh_rd_fiduc);
      break;
    case PROP_DT_RD_FIDUC:
      g_value_set_double (value, bao_ef->Dt_rd_fiduc);
      break;
    case PROP_Z:
      g_value_set_double (value, bao_ef->z);
      break;
    case PROP_M2LNP:
      g_value_set_object (value, bao_ef->m2lnp);
      break;
    case PROP_DIST:
      g_value_set_object (value, bao_ef->dist);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_empirical_fit_2d_dispose (GObject *object)
{
  NcDataBaoEmpiricalFit2d *bao_ef = NC_DATA_BAO_EMPIRICAL_FIT_2D (object);

  ncm_spline2d_clear (&bao_ef->m2lnp);
  ncm_stats_dist2d_clear (&bao_ef->p);

  nc_distance_clear (&bao_ef->dist);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_empirical_fit_2d_parent_class)->dispose (object);
}

static void
nc_data_bao_empirical_fit_2d_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_empirical_fit_2d_parent_class)->finalize (object);
}

static gdouble _nc_data_bao_empirical_fit_2d_m2lnL_val (NcmDataDist2d *dist2d, NcmMSet *mset, gdouble x, gdouble y);
static void _nc_data_bao_empirical_fit_2d_inv_pdf (NcmDataDist2d *dist2d, NcmMSet *mset, gdouble u, gdouble v, gdouble *x, gdouble *y);

static void
nc_data_bao_empirical_fit_2d_class_init (NcDataBaoEmpiricalFit2dClass *klass)
{
  GObjectClass *object_class       = G_OBJECT_CLASS (klass);
  NcmDataDist2dClass *bao_ef_class = NCM_DATA_DIST2D_CLASS (klass);

  object_class->constructed  = nc_data_bao_empirical_fit_2d_constructed;
  object_class->set_property = nc_data_bao_empirical_fit_2d_set_property;
  object_class->get_property = nc_data_bao_empirical_fit_2d_get_property;
  object_class->dispose      = nc_data_bao_empirical_fit_2d_dispose;
  object_class->finalize     = nc_data_bao_empirical_fit_2d_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DH_RD_FIDUC,
                                   g_param_spec_double ("Dh-rd-fiduc",
                                                        NULL,
                                                        "Dh/rd fiducial",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_DT_RD_FIDUC,
                                   g_param_spec_double ("Dt-rd-fiduc",
                                                        NULL,
                                                        "Dt/rd fiducial",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_Z,
                                   g_param_spec_double ("z",
                                                        NULL,
                                                        "Redshift",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_M2LNP,
                                   g_param_spec_object ("m2lnp",
                                                        NULL,
                                                        "Empirical m2lnp",
                                                        NCM_TYPE_SPLINE2D,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  bao_ef_class->m2lnL_val = _nc_data_bao_empirical_fit_2d_m2lnL_val;
  bao_ef_class->inv_pdf   = _nc_data_bao_empirical_fit_2d_inv_pdf;
}

static gdouble
_nc_data_bao_empirical_fit_2d_m2lnL_val (NcmDataDist2d *dist2d, NcmMSet *mset, gdouble x, gdouble y)
{
  NcDataBaoEmpiricalFit2d *bao_ef = NC_DATA_BAO_EMPIRICAL_FIT_2D (dist2d);
  const gdouble alpha_par         = nc_data_bao_empirical_fit_2d_get_alpha_parallel (bao_ef, mset);
  const gdouble alpha_per         = nc_data_bao_empirical_fit_2d_get_alpha_perpendicular (bao_ef, mset);
  const gdouble alphax            = alpha_per - x;
  const gdouble alphay            = alpha_par - y;
  gdouble m2lnL                   = ncm_stats_dist2d_eval_m2lnp (bao_ef->p, alphax, alphay);
  gdouble xl, xu, yl, yu;

  ncm_stats_dist2d_xbounds (bao_ef->p, &xl, &xu);
  ncm_stats_dist2d_ybounds (bao_ef->p, &yl, &yu);

/*
 *  printf ("x, y = % 22.15g, % 22.15g | m2lnL = % 22.15g (% 22.15g % 22.15g) (% 22.15g % 22.15g)\n",
 *     alphax, alphay, m2lnL,
 *     xl, xu, yl, yu);
 */

  if (((alphax < xl) || (alphax > xu)) || ((alphay < yl) || (alphay > yu)))
    return GSL_POSINF;
  else
    return m2lnL;
}

static void
_nc_data_bao_empirical_fit_2d_inv_pdf (NcmDataDist2d *dist2d, NcmMSet *mset, gdouble u, gdouble v, gdouble *x, gdouble *y)
{
  NcDataBaoEmpiricalFit2d *bao_ef = NC_DATA_BAO_EMPIRICAL_FIT_2D (dist2d);

  NCM_UNUSED (bao_ef);
  g_assert_not_reached ();
}

/**
 * nc_data_bao_empirical_fit_2d_new:
 * @m2lnp: a #NcmSpline2d containing $-2\ln (p)$
 * @Dh_rd_fiduc: fiducial $D_H/r_d$
 * @Dt_rd_fiduc: fiducial $D_t/r_d$
 * @z: data redshift
 *
 * Creates a new empirical BAO data object.
 *
 * Returns: the newly created #NcDataBaoEmpiricalFit2d.
 */
NcDataBaoEmpiricalFit2d *
nc_data_bao_empirical_fit_2d_new (NcmSpline2d *m2lnp, gdouble Dh_rd_fiduc, gdouble Dt_rd_fiduc, gdouble z)
{
  NcDataBaoEmpiricalFit2d *bao_ef;
  NcmMatrix *m = ncm_matrix_new (1, 2);

  ncm_matrix_set (m, 0, 0, 0.0);
  ncm_matrix_set (m, 0, 1, 0.0);

  bao_ef = g_object_new (NC_TYPE_DATA_BAO_EMPIRICAL_FIT_2D,
                         "m2lnp", m2lnp,
                         "Dh-rd-fiduc", Dh_rd_fiduc,
                         "Dt-rd-fiduc", Dt_rd_fiduc,
                         "z", z,
                         "n-points", 1,
                         "matrix", m,
                         NULL);
  ncm_matrix_free (m);

  return bao_ef;
}

/**
 * nc_data_bao_empirical_fit_2d_new_from_file:
 * @filename: file containing a serialized #NcDataBaoEmpiricalFit2d.
 *
 * Creates a new #NcDataBaoEmpiricalFit2d from @filename.
 *
 * Returns: (transfer full): the newly created #NcDataBaoEmpiricalFit2d.
 */
NcDataBaoEmpiricalFit2d *
nc_data_bao_empirical_fit_2d_new_from_file (const gchar *filename)
{
  NcDataBaoEmpiricalFit2d *bao_ef = NC_DATA_BAO_EMPIRICAL_FIT_2D (ncm_serialize_global_from_file (filename));

  g_assert (NC_IS_DATA_BAO_EMPIRICAL_FIT_2D (bao_ef));

  return bao_ef;
}

/**
 * nc_data_bao_empirical_fit_2d_new_from_id:
 * @dist: a #NcDistance
 * @id: a #NcDataBaoId
 *
 * Creates a new #NcDataBaoEmpiricalFit2d from @id.
 *
 * Returns: (transfer full): the newly created #NcDataBaoEmpiricalFit2d.
 */
NcDataBaoEmpiricalFit2d *
nc_data_bao_empirical_fit_2d_new_from_id (NcDistance *dist, NcDataBaoId id)
{
  NcDataBaoEmpiricalFit2d *bao_ef;
  gchar *filename;

  switch (id)
  {
    case NC_DATA_BAO_EMPIRICAL_FIT_2D_BAUTISTA2017:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_empirical_fit_2d_bautista2017.obj", TRUE);
      break;
    case NC_DATA_BAO_EMPIRICAL_FIT_2D_SDSS_DR16_LYAUTO_2021:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_empirical_fit_2d_sdss_dr16_lyauto_2021.obj", TRUE);
      break;
    case NC_DATA_BAO_EMPIRICAL_FIT_2D_SDSS_DR16_LYXQSO_2021:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_empirical_fit_2d_sdss_dr16_lyxqso_2021.obj", TRUE);
      break;
    default:
      g_error ("nc_data_bao_empirical_fit_2d_new_from_id: id %d not recognized.", id);
      break;
  }

  bao_ef = nc_data_bao_empirical_fit_2d_new_from_file (filename);
  nc_data_bao_empirical_fit_2d_set_dist (bao_ef, dist);
  g_free (filename);

  g_assert (NC_IS_DATA_BAO_EMPIRICAL_FIT_2D (bao_ef));
  g_assert (NCM_IS_DATA (bao_ef));

  return bao_ef;
}

/**
 * nc_data_bao_empirical_fit_2d_get_alpha_perpendicular:
 * @bao_ef: a #NcDataBaoEmpiricalFit2d
 * @mset: a #NcmMSet
 *
 * Calculates value of $\alpha_{\perp}$ given a #NcmMSet,
 * $$ \alpha_{perp} = \frac{[D_t(z)/r_d]}{[D_t(z)/r_d]_{fid}},$$
 * where $D_t(z)$ is the transverse comoving distance [nc_distance_transverse()], $r_d$ is the sound
 * horizon [nc_distance_sound_horizon()] at the drag epoch, and 'fid' indicates fiducial.
 *
 * Returns: $\alpha_{\perp}$
 */
gdouble
nc_data_bao_empirical_fit_2d_get_alpha_perpendicular (NcDataBaoEmpiricalFit2d *bao_ef, NcmMSet *mset)
{
  NcHICosmo *cosmo         = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  const gdouble Dt_r       = (1.0 + bao_ef->z) * nc_distance_DA_r (bao_ef->dist, cosmo, bao_ef->z);
  const gdouble Dt_r_fiduc = bao_ef->Dt_rd_fiduc;
  const gdouble alpha_per  = Dt_r / Dt_r_fiduc;

  return alpha_per;
}

/**
 * nc_data_bao_empirical_fit_2d_get_alpha_parallel:
 * @bao_ef: a #NcDataBaoEmpiricalFit2d
 * @mset: a #NcmMSet
 *
 * Calculates value of $\alpha_{\parallel}$ given a #NcmMSet.
 *
 * Returns: $\alpha_{\parallel}$
 */
gdouble
nc_data_bao_empirical_fit_2d_get_alpha_parallel (NcDataBaoEmpiricalFit2d *bao_ef, NcmMSet *mset)
{
  NcHICosmo *cosmo         = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  const gdouble Dh_r       = nc_distance_DH_r (bao_ef->dist, cosmo, bao_ef->z);
  const gdouble Dh_r_fiduc = bao_ef->Dh_rd_fiduc;
  const gdouble alpha_par  = Dh_r / Dh_r_fiduc;

  return alpha_par;
}

/**
 * nc_data_bao_empirical_fit_2d_set_dist:
 * @bao_ef: a #NcDataBaoEmpiricalFit2d
 * @dist: a #NcDistance
 *
 * Sets the distance object.
 *
 */
void
nc_data_bao_empirical_fit_2d_set_dist (NcDataBaoEmpiricalFit2d *bao_ef, NcDistance *dist)
{
  nc_distance_clear (&bao_ef->dist);
  bao_ef->dist = nc_distance_ref (dist);
}

