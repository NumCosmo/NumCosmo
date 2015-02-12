/***************************************************************************
 *            nc_data_bao_empirical_fit.c
 *
 *  Wed February 11 13:03:03 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_data_bao_empirical_fit.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_data_bao_empirical_fit
 * @title: Baryon Oscillation data -- Dv/r
 * @short_description: full likelihood
 *
 * This object implements the BAO data when its likelihood function is provided, 
 * e.g., <link linkend="XRoss2014">Ross et al. (2015)</link>. 
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_bao_empirical_fit.h"
#include "nc_hicosmo.h"
#include "nc_distance.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_DV_FIDUC,
  PROP_RS_FIDUC,
  PROP_Z,
  PROP_M2LNP,
  PROP_DIST,
};

G_DEFINE_TYPE (NcDataBaoEmpiricalFit, nc_data_bao_empirical_fit, NCM_TYPE_DATA_DIST1D);

static void
nc_data_bao_empirical_fit_init (NcDataBaoEmpiricalFit *bao_ef)
{
  bao_ef->Dv_fiduc = 0.0;
  bao_ef->rs_fiduc = 0.0;
  bao_ef->z        = 0.0;
  bao_ef->m2lnp    = NULL;
  bao_ef->p        = NULL;
  bao_ef->p_mode   = 0.0;
  bao_ef->dist     = nc_distance_new (2.0);
}

static void
nc_data_bao_empirical_fit_constructed (GObject *object)
{
  NcDataBaoEmpiricalFit *bao_ef = NC_DATA_BAO_EMPIRICAL_FIT (object);

  ncm_stats_dist1d_prepare (bao_ef->p);
  bao_ef->p_mode = ncm_stats_dist1d_mode (bao_ef->p);

  ncm_data_set_init (NCM_DATA (bao_ef), TRUE);
}

static void
nc_data_bao_empirical_fit_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataBaoEmpiricalFit *bao_ef = NC_DATA_BAO_EMPIRICAL_FIT (object);
  g_return_if_fail (NC_IS_DATA_BAO_EMPIRICAL_FIT (object));

  switch (prop_id)
  {
    case PROP_DV_FIDUC:
      bao_ef->Dv_fiduc = g_value_get_double (value);
      break;
    case PROP_RS_FIDUC:
      bao_ef->rs_fiduc = g_value_get_double (value);
      break;
    case PROP_Z:
      bao_ef->z = g_value_get_double (value);
      break;
    case PROP_M2LNP:
      ncm_spline_clear (&bao_ef->m2lnp);
      ncm_stats_dist1d_clear (&bao_ef->p);
      bao_ef->m2lnp = g_value_dup_object (value);
      bao_ef->p = ncm_stats_dist1d_spline_new (bao_ef->m2lnp);
      break;
    case PROP_DIST:
      nc_distance_clear (&bao_ef->dist);
      bao_ef->dist = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_empirical_fit_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataBaoEmpiricalFit *bao_ef = NC_DATA_BAO_EMPIRICAL_FIT (object);
  g_return_if_fail (NC_IS_DATA_BAO_EMPIRICAL_FIT (object));

  switch (prop_id)
  {
    case PROP_DV_FIDUC:
      g_value_set_double (value, bao_ef->Dv_fiduc);
      break;
    case PROP_RS_FIDUC:
      g_value_set_double (value, bao_ef->rs_fiduc);
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
nc_data_bao_empirical_fit_dispose (GObject *object)
{
  NcDataBaoEmpiricalFit *bao_ef = NC_DATA_BAO_EMPIRICAL_FIT (object);

  ncm_spline_clear (&bao_ef->m2lnp);
  ncm_stats_dist1d_clear (&bao_ef->p);
  
  nc_distance_clear (&bao_ef->dist);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_empirical_fit_parent_class)->dispose (object);
}

static void
nc_data_bao_empirical_fit_finalize (GObject *object)
{
  

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_empirical_fit_parent_class)->finalize (object);
}

static gdouble _nc_data_bao_empirical_fit_m2lnL_val (NcmDataDist1d *dist1d, NcmMSet *mset, gdouble x);
static gdouble _nc_data_bao_empirical_fit_inv_pdf (NcmDataDist1d *dist1d, NcmMSet *mset, gdouble u);

static void
nc_data_bao_empirical_fit_class_init (NcDataBaoEmpiricalFitClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataDist1dClass *bao_ef_class = NCM_DATA_DIST1D_CLASS (klass);

  object_class->constructed  = nc_data_bao_empirical_fit_constructed;
  object_class->set_property = nc_data_bao_empirical_fit_set_property;
  object_class->get_property = nc_data_bao_empirical_fit_get_property;
  object_class->dispose      = nc_data_bao_empirical_fit_dispose;
  object_class->finalize     = nc_data_bao_empirical_fit_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DV_FIDUC,
                                   g_param_spec_double ("Dv-fiduc",
                                                        NULL,
                                                        "Dv fiducial",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_RS_FIDUC,
                                   g_param_spec_double ("rs-fiduc",
                                                        NULL,
                                                        "r_s fiducial",
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
                                                        NCM_TYPE_SPLINE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  bao_ef_class->m2lnL_val = _nc_data_bao_empirical_fit_m2lnL_val;
  bao_ef_class->inv_pdf   = _nc_data_bao_empirical_fit_inv_pdf;
}

static gdouble 
_nc_data_bao_empirical_fit_m2lnL_val (NcmDataDist1d *dist1d, NcmMSet *mset, gdouble x)
{
  NcDataBaoEmpiricalFit *bao_ef = NC_DATA_BAO_EMPIRICAL_FIT (dist1d);
  const gdouble alpha  = nc_data_bao_empirical_fit_get_alpha (bao_ef, mset);
  const gdouble alphap = alpha - x;

  return ncm_stats_dist1d_m2lnp (bao_ef->p, alphap);
}

static gdouble 
_nc_data_bao_empirical_fit_inv_pdf (NcmDataDist1d *dist1d, NcmMSet *mset, gdouble u)
{
  NcDataBaoEmpiricalFit *bao_ef = NC_DATA_BAO_EMPIRICAL_FIT (dist1d);
  return ncm_stats_dist1d_inv_pdf (bao_ef->p, u) - bao_ef->p_mode;
}

/**
 * nc_data_bao_empirical_fit_new:
 * @m2lnp: a #NcmSpline containing $-2\ln (p)$
 * @Dv_fiduc: fiducial $D_V$
 * @rs_fiduc: fiducial $r_s$
 * @z: data redshift
 * 
 * Creates a new #NcDataBaoEmpiricalFit.
 * 
 * Returns: the newly created #NcDataBaoEmpiricalFit.
 */
NcDataBaoEmpiricalFit *
nc_data_bao_empirical_fit_new (NcmSpline *m2lnp, gdouble Dv_fiduc, gdouble rs_fiduc, gdouble z)
{
  NcDataBaoEmpiricalFit *bao_ef;
  NcmVector *d = ncm_vector_new (1);

  ncm_vector_set (d, 0, 0.0);

  bao_ef = g_object_new (NC_TYPE_DATA_BAO_EMPIRICAL_FIT,
                         "m2lnp", m2lnp,
                         "Dv-fiduc", Dv_fiduc,
                         "rs-fiduc", rs_fiduc,
                         "z", z,
                         "n-points", 1,
                         "vector", d,
                         NULL);
  ncm_vector_free (d);

  return bao_ef;
}

/**
 * nc_data_bao_empirical_fit_new_from_file:
 * @filename: file containing a serialized #NcDataBaoEmpiricalFit.
 * 
 * Creates a new #NcDataBaoEmpiricalFit from @filename.
 * 
 * Returns: (transfer full): the newly created #NcDataBaoEmpiricalFit.
 */
NcDataBaoEmpiricalFit *
nc_data_bao_empirical_fit_new_from_file (const gchar *filename)
{
  NcDataBaoEmpiricalFit *bao_ef = NC_DATA_BAO_EMPIRICAL_FIT (ncm_serialize_global_from_file (filename));
  g_assert (NC_IS_DATA_BAO_EMPIRICAL_FIT (bao_ef));

  return bao_ef;
}

/**
 * nc_data_bao_empirical_fit_new_from_id:
 * @id: a #NcDataBaoId
 * 
 * Creates a new #NcDataBaoEmpiricalFit from @id.
 * 
 * Returns: (transfer full): the newly created #NcDataBaoEmpiricalFit.
 */
NcDataBaoEmpiricalFit *
nc_data_bao_empirical_fit_new_from_id (NcDataBaoId id)
{
  NcDataBaoEmpiricalFit *bao_ef;
  gchar *filename;
  switch (id)
  {
    case NC_DATA_BAO_EMPIRICAL_FIT_ROSS2015:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_empirical_fit_ross2015.obj", TRUE);
      break;
    default:
      g_error ("nc_data_bao_empirical_fit_new_from_id: id %d not recognized.", id);
      break;
  }

  bao_ef = nc_data_bao_empirical_fit_new_from_file (filename);
  g_free (filename);

  return bao_ef;
}

/**
 * nc_data_bao_empirical_fit_get_mode:
 * @bao_ef: a #NcDataBaoEmpiricalFit
 * 
 * Calculates the mode of the empirical distribution.
 * 
 * Returns: the mode of the distribution.
 */
gdouble 
nc_data_bao_empirical_fit_get_mode (NcDataBaoEmpiricalFit *bao_ef)
{
  return bao_ef->p_mode;
}

/**
 * nc_data_bao_empirical_fit_get_alpha:
 * @bao_ef: a #NcDataBaoEmpiricalFit
 * @mset: a #NcmMSet
 * 
 * Calculates value of $\alpha$ given a #NcmMSet.
 * 
 * Returns: $\alpha$
 */
gdouble 
nc_data_bao_empirical_fit_get_alpha (NcDataBaoEmpiricalFit *bao_ef, NcmMSet *mset)
{
  NcHICosmo *cosmo         = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  const gdouble Dv_r       = 1.0 / nc_distance_bao_r_Dv (bao_ef->dist, cosmo, bao_ef->z);
  const gdouble Dv_r_fiduc = bao_ef->Dv_fiduc / bao_ef->rs_fiduc;
  const gdouble alpha      = Dv_r / Dv_r_fiduc;

  return alpha;
}
