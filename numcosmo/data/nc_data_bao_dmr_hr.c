/***************************************************************************
 *            nc_data_bao_dmr_hr.c
 *
 *  Mon Oct 10 23:29:40 2016
 *  Copyright  2016  Mariana Penna Lima
 *  <pennalima@gmail.com>
 *****************************************************************************/
/*
 * nc_data_bao_dmr_hr.c
 * Copyright (C) 2016 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_data_bao_dmr_hr
 * @title: NcDataBaoDMrHr
 * @short_description: Baryon Oscillation Data -- $(D_M/r,\; H/r)$ data.
 *
 * See: [Ross et al. (2016)][XRoss2016] and
 * [Alam et al. (2016)][XAlam2016].
 *
 * The data is stored in a #NcDataBaoDMrHr object. The data is stored in a
 * #NcDataGaussCov base class object, which is a subclass of #NcmData.
 * The data represents the mean values of the transverse distance $D_M$ and the
 * Hubble parameter $H$ at the redshift $z$ divided by the sound horizon at the
 * last scattering surface $r_s$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_bao_dmr_hr.h"
#include "nc_enum_types.h"
#include "math/ncm_cfg.h"

enum
{
  PROP_0,
  PROP_DIST,
  PROP_Z,
  PROP_RS_FIDUC,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataBaoDMrHr, nc_data_bao_dmr_hr, NCM_TYPE_DATA_GAUSS_COV);

static void
nc_data_bao_dmr_hr_init (NcDataBaoDMrHr *dmh)
{
  dmh->dist     = nc_distance_new (2.0);
  dmh->x        = NULL;
  dmh->rs_fiduc = 0.0;
}

static void
nc_data_bao_dmr_hr_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataBaoDMrHr *dmh = NC_DATA_BAO_DMR_HR (object);

  g_return_if_fail (NC_IS_DATA_BAO_DMR_HR (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_data_bao_dmr_hr_set_dist (dmh, g_value_get_object (value));
      break;
    case PROP_Z:
      ncm_vector_substitute (&dmh->x, g_value_get_object (value), TRUE);
      break;
    case PROP_RS_FIDUC:
      dmh->rs_fiduc = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_dmr_hr_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataBaoDMrHr *dmh = NC_DATA_BAO_DMR_HR (object);

  g_return_if_fail (NC_IS_DATA_BAO_DMR_HR (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, dmh->dist);
      break;
    case PROP_Z:
      g_value_set_object (value, dmh->x);
      break;
    case PROP_RS_FIDUC:
      g_value_set_double (value, dmh->rs_fiduc);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_bao_dmr_hr_dispose (GObject *object)
{
  NcDataBaoDMrHr *dmh = NC_DATA_BAO_DMR_HR (object);

  nc_distance_clear (&dmh->dist);
  ncm_vector_clear (&dmh->x);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_dmr_hr_parent_class)->dispose (object);
}

static void
nc_data_bao_dmr_hr_finalize (GObject *object)
{
  /*NcDataBaoDMrHr *dmh = NC_DATA_BAO_DMR_HR (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_bao_dmr_hr_parent_class)->finalize (object);
}

static void _nc_data_bao_dmr_hr_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_bao_dmr_hr_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);
static void _nc_data_bao_dmr_hr_set_size (NcmDataGaussCov *gauss, guint np);

static void
nc_data_bao_dmr_hr_class_init (NcDataBaoDMrHrClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class          = NCM_DATA_CLASS (klass);
  NcmDataGaussCovClass *gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->set_property = &nc_data_bao_dmr_hr_set_property;
  object_class->get_property = &nc_data_bao_dmr_hr_get_property;
  object_class->dispose      = &nc_data_bao_dmr_hr_dispose;
  object_class->finalize     = &nc_data_bao_dmr_hr_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_Z,
                                   g_param_spec_object ("z",
                                                        NULL,
                                                        "Data redshift",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_RS_FIDUC,
                                   g_param_spec_double ("rs-fiduc",
                                                        NULL,
                                                        "r_s fiducial",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  data_class->prepare    = &_nc_data_bao_dmr_hr_prepare;
  gauss_class->mean_func = &_nc_data_bao_dmr_hr_mean_func;
  gauss_class->set_size  = &_nc_data_bao_dmr_hr_set_size;
}

static void
_nc_data_bao_dmr_hr_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataBaoDMrHr *dmh = NC_DATA_BAO_DMR_HR (data);
  NcHICosmo *cosmo    = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));

  nc_distance_prepare_if_needed (dmh->dist, cosmo);
}

static void
_nc_data_bao_dmr_hr_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataBaoDMrHr *dmh = NC_DATA_BAO_DMR_HR (gauss);
  NcHICosmo *cosmo    = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  const gdouble rd    = nc_distance_r_zd (dmh->dist, cosmo) * nc_hicosmo_RH_Mpc (cosmo);
  const guint np      = ncm_data_gauss_cov_get_size (gauss);
  gint i;

  for (i = 0; i < np; i++)
  {
    if (i % 2 == 0)
    {
      const gdouble z1   = ncm_vector_get (dmh->x, i);
      const gdouble z2   = ncm_vector_get (dmh->x, i + 1);
      const gdouble DM_r = (1.0 + z1) * nc_distance_DA_r (dmh->dist, cosmo, z1) * dmh->rs_fiduc;
      const gdouble Hr   = nc_hicosmo_H (cosmo, z2) * rd / dmh->rs_fiduc;

      ncm_vector_set (vp, i, DM_r);
      ncm_vector_set (vp, i + 1, Hr);
    }
    else
    {
      continue;
    }
  }
}

static void
_nc_data_bao_dmr_hr_set_size (NcmDataGaussCov *gauss, guint np)
{
  NcDataBaoDMrHr *dmh = NC_DATA_BAO_DMR_HR (gauss);
  const guint cnp     = ncm_data_gauss_cov_get_size (gauss);

  if ((np == 0) || (np != cnp))
    ncm_vector_clear (&dmh->x);

  if ((np != 0) && (np != cnp))
    dmh->x = ncm_vector_new (np);

  /* Chain up : end */
  NCM_DATA_GAUSS_COV_CLASS (nc_data_bao_dmr_hr_parent_class)->set_size (gauss, np);
}

/**
 * nc_data_bao_dmr_hr_new_from_file:
 * @filename: file containing a serialized #NcDataBaoDMrHr
 *
 * Creates a new #NcDataBaoDMrHr from @filename.
 *
 * Returns: (transfer full): the newly created #NcDataBaoDMrHr.
 */
NcDataBaoDMrHr *
nc_data_bao_dmr_hr_new_from_file (const gchar *filename)
{
  NcDataBaoDMrHr *dmh = NC_DATA_BAO_DMR_HR (ncm_serialize_global_from_file (filename));

  g_assert (NC_IS_DATA_BAO_DMR_HR (dmh));

  return dmh;
}

/**
 * nc_data_bao_dmr_hr_new_from_id:
 * @dist: a #NcDistance
 * @id: a #NcDataBaoId
 *
 * Creates a new acustic scale data object #NcDataBaoDMrHr from @id.
 * This object requires a #NcDistance object to be set.
 *
 *
 * Returns: (transfer full): a #NcDataBaoDMrHr
 */
NcDataBaoDMrHr *
nc_data_bao_dmr_hr_new_from_id (NcDistance *dist, NcDataBaoId id)
{
  NcDataBaoDMrHr *dmh;
  gchar *filename;

  switch (id)
  {
    case NC_DATA_BAO_DMR_HR_SDSS_DR12_2016:
      filename = ncm_cfg_get_data_filename ("nc_data_bao_dmr_hr_sdss_dr12_2016.obj", TRUE);
      break;
    default:
      g_error ("nc_data_bao_dmr_hr_new_from_id: id %d not recognized.", id);
      break;
  }

  dmh = nc_data_bao_dmr_hr_new_from_file (filename);
  nc_data_bao_dmr_hr_set_dist (dmh, dist);
  g_free (filename);

  return dmh;
}

/**
 * nc_data_bao_dmr_hr_set_dist:
 * @dmh: a #NcDataBaoDMrHr
 * @dist: a #NcDistance
 *
 * Sets the distance object.
 *
 */
void
nc_data_bao_dmr_hr_set_dist (NcDataBaoDMrHr *dmh, NcDistance *dist)
{
  nc_distance_clear (&dmh->dist);
  dmh->dist = nc_distance_ref (dist);
}

