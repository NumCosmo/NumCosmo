/***************************************************************************
 *            nc_data_cluster_counts_box_poisson.c
 *
 *  Mon Feb 20 15:31:47 2017
 *  Copyright  2017  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
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
 * SECTION:nc_data_cluster_counts_box_poisson
 * @title: NcDataClusterCountsBoxPoisson
 * @short_description: Binned (mass) cluster number count data in a box.
 *
 * FIXME Box (simulation), and not redshift space.
 * 
 */

enum
{
  PROP_0,
	PROP_MASS_FUNC, 
  PROP_MASS_KNOTS,
	PROP_REDSHIFT, 
  PROP_VOLUME,
  PROP_SIZE,
};

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_cluster_counts_box_poisson.h"
#include "lss/nc_halo_mass_function.h"
#include "math/ncm_util.h"
#include "math/ncm_spline_cubic_notaknot.h"

G_DEFINE_TYPE (NcDataClusterCountsBoxPoisson, nc_data_cluster_counts_box_poisson, NCM_TYPE_DATA_POISSON);

static void
nc_data_cluster_counts_box_poisson_init (NcDataClusterCountsBoxPoisson *cpoisson)
{
  cpoisson->mfp        = NULL;
	cpoisson->mass_knots = NULL;
	cpoisson->redshift   = 0.0;
	cpoisson->volume     = 0.0;
	cpoisson->dndlog10M  = ncm_spline_cubic_notaknot_new ();
}

static void
nc_data_cluster_counts_box_poisson_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataClusterCountsBoxPoisson *cpoisson = NC_DATA_CLUSTER_COUNTS_BOX_POISSON (object);
  g_return_if_fail (NC_IS_DATA_CLUSTER_COUNTS_BOX_POISSON (object));

  switch (prop_id)
  {
		case PROP_MASS_FUNC:
			cpoisson->mfp = g_value_dup_object (value);
			g_assert (cpoisson->mfp != NULL);
			break;
		case PROP_MASS_KNOTS:
			ncm_vector_memcpy (cpoisson->mass_knots, g_value_get_object (value));
			break;
		case PROP_REDSHIFT:
			cpoisson->redshift = g_value_get_double (value);
			break;	
		case PROP_VOLUME:
			cpoisson->volume = g_value_get_double (value);
			break;	  
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
nc_data_cluster_counts_box_poisson_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
	NcDataClusterCountsBoxPoisson *cpoisson = NC_DATA_CLUSTER_COUNTS_BOX_POISSON (object);
	g_return_if_fail (NC_IS_DATA_CLUSTER_COUNTS_BOX_POISSON (object));

	switch (prop_id)
	{
		case PROP_MASS_FUNC:
			g_value_set_object (value, cpoisson->mfp);
			break;
		case PROP_MASS_KNOTS:
			g_value_set_object (value, cpoisson->mass_knots);
			break;
		case PROP_REDSHIFT:
			g_value_set_double (value, cpoisson->redshift);
			break;	
		case PROP_VOLUME:
			g_value_set_double (value, cpoisson->volume);
			break;	  
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
nc_data_cluster_counts_box_poisson_dispose (GObject *object)
{
  NcDataClusterCountsBoxPoisson *cpoisson = NC_DATA_CLUSTER_COUNTS_BOX_POISSON (object);

	nc_halo_mass_function_clear (&cpoisson->mfp);
	ncm_spline_clear (&cpoisson->dndlog10M); 
	
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_counts_box_poisson_parent_class)->dispose (object);
}

static void
nc_data_cluster_counts_box_poisson_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_cluster_counts_box_poisson_parent_class)->finalize (object);
}

static void _nc_data_cluster_counts_box_poisson_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_cluster_counts_box_poisson_set_size (NcmDataPoisson *poisson, guint nbins);
static gdouble _nc_data_cluster_counts_box_poisson_mean_func (NcmDataPoisson *poisson, NcmMSet *mset, guint n);

static void
nc_data_cluster_counts_box_poisson_class_init (NcDataClusterCountsBoxPoissonClass *klass)
{
  GObjectClass* object_class         = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class           = NCM_DATA_CLASS (klass);
	NcmDataPoissonClass *poisson_class = NCM_DATA_POISSON_CLASS (klass);

  object_class->set_property = &nc_data_cluster_counts_box_poisson_set_property;
  object_class->get_property = &nc_data_cluster_counts_box_poisson_get_property;
  object_class->dispose      = &nc_data_cluster_counts_box_poisson_dispose;
  object_class->finalize     = &nc_data_cluster_counts_box_poisson_finalize;

	g_object_class_install_property (object_class,
	                                 PROP_MASS_FUNC,
	                                 g_param_spec_object ("mass-function",
	                                                      NULL,
	                                                      "Halo mass function",
	                                                      NC_TYPE_HALO_MASS_FUNCTION,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_MASS_KNOTS,
	                                 g_param_spec_object ("mass-knots",
	                                                      NULL,
	                                                      "The mass knots",
	                                                      NCM_TYPE_VECTOR,
	                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
                                   PROP_REDSHIFT,
                                   g_param_spec_double ("redshift",
                                                        NULL,
                                                        "Redshift",
                                                        0.0, 2.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	g_object_class_install_property (object_class,
                                   PROP_VOLUME,
                                   g_param_spec_double ("volume",
                                                        NULL,
                                                        "Volume of the box [Mpc^3 / h^3](simulation)",
                                                        10, 1.0e13, 452984832000.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  data_class->prepare      = &_nc_data_cluster_counts_box_poisson_prepare;
	poisson_class->set_size  = &_nc_data_cluster_counts_box_poisson_set_size;
	poisson_class->mean_func = &_nc_data_cluster_counts_box_poisson_mean_func;
}

typedef struct _NcDataClusterCountsBoxPoissonArg
{
  NcDataClusterCountsBoxPoisson *cpoisson;
	NcHICosmo *cosmo;
} NcDataClusterCountsBoxPoissonArg;

static gdouble 
_nc_data_cluster_counts_box_poisson_dndlog10M (gdouble log10M, gpointer userdata)
{
  NcDataClusterCountsBoxPoissonArg *arg = (NcDataClusterCountsBoxPoissonArg *) userdata;
	const gdouble lnM = log10M * M_LN10;
	return M_LN10 * nc_halo_mass_function_dn_dlnM (arg->cpoisson->mfp, arg->cosmo, lnM, arg->cpoisson->redshift);
}

static void 
_nc_data_cluster_counts_box_poisson_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataClusterCountsBoxPoisson *cpoisson = NC_DATA_CLUSTER_COUNTS_BOX_POISSON (data);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
	NcDataClusterCountsBoxPoissonArg arg = {cpoisson, cosmo};
  const gdouble log10Ml = ncm_vector_get (cpoisson->mass_knots, 0);
	const gdouble log10Mu = ncm_vector_get (cpoisson->mass_knots, ncm_vector_len (cpoisson->mass_knots) - 1);
	gsl_function F;

  g_assert (cpoisson->mfp != NULL);
	
  nc_halo_mass_function_prepare_if_needed (cpoisson->mfp, cosmo);

	F.function = &_nc_data_cluster_counts_box_poisson_dndlog10M;
	F.params   = &arg;
	
  ncm_spline_set_func (cpoisson->dndlog10M, NCM_SPLINE_FUNCTION_SPLINE, &F, log10Ml, log10Mu, 0, 1.0e-9);
}

static void 
_nc_data_cluster_counts_box_poisson_set_size (NcmDataPoisson *poisson, guint nbins)
{
  /* Chain up : start */
  NCM_DATA_POISSON_CLASS (nc_data_cluster_counts_box_poisson_parent_class)->set_size (poisson, nbins);

	{
		NcDataClusterCountsBoxPoisson *cpoisson = NC_DATA_CLUSTER_COUNTS_BOX_POISSON (poisson);
		ncm_vector_clear (&cpoisson->mass_knots);
		if (nbins > 0)
		  cpoisson->mass_knots = ncm_vector_new_data_static (poisson->h->range, nbins + 1, 1);
	}
}

static gdouble 
_nc_data_cluster_counts_box_poisson_mean_func (NcmDataPoisson *poisson, NcmMSet *mset, guint n)
{
  NcDataClusterCountsBoxPoisson *cpoisson = NC_DATA_CLUSTER_COUNTS_BOX_POISSON (poisson);
  /*NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));*/
  const gdouble log10Mn   = ncm_vector_get (cpoisson->mass_knots, n);
	const gdouble log10Mnp1 = ncm_vector_get (cpoisson->mass_knots, n + 1);

  return ncm_spline_eval_integ (cpoisson->dndlog10M, log10Mn, log10Mnp1) * cpoisson->volume;
}

/**
 * nc_data_cluster_counts_box_poisson_new:
 * @mfp: a #NcHaloMassFunction
 * 
 * FIXME
 *
 * Returns: FIXME
 */
NcDataClusterCountsBoxPoisson *
nc_data_cluster_counts_box_poisson_new (NcHaloMassFunction *mfp)
{
	NcDataClusterCountsBoxPoisson *cpoisson = g_object_new (NC_TYPE_DATA_CLUSTER_COUNTS_BOX_POISSON,
	                                                        "mass-function", mfp,
	                                                        NULL);
	return cpoisson;
}

/**
 * nc_data_cluster_counts_box_poisson_init_from_sampling:
 * @cpoisson: a #NcDataClusterCountsBoxPoisson
 * @mset: a #NcmMSet
 * @mass_knots: a #NcmVector containing the histogram knots
 * @volume: box volume
 * @redshift: box redshift
 * @rng: a #NcmRNG
 *
 * FIXME
 *
 */
void
nc_data_cluster_counts_box_poisson_init_from_sampling (NcDataClusterCountsBoxPoisson *cpoisson, NcmMSet *mset, NcmVector *mass_knots, const gdouble volume, const gdouble redshift, NcmRNG *rng)
{
  const gint nbins = ncm_vector_len (mass_knots) - 1;
	
	g_assert_cmpuint (nbins, >, 0);

	ncm_data_poisson_set_size (NCM_DATA_POISSON (cpoisson), nbins);
  ncm_vector_memcpy (cpoisson->mass_knots, mass_knots);

	cpoisson->redshift = redshift;
	cpoisson->volume   = volume;

	ncm_data_resample (NCM_DATA (cpoisson), mset, rng);
}
