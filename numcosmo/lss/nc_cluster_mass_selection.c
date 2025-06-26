/***************************************************************************
 *            nc_cluster_mass_selection.c
 *
 *  Thu Jan 26 18:25:11 2017
 *  Copyright  2017  Mariana Penna Lima and Begoña Selection
 *  <pennalima@gmail.com>, <bego.selection.work@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima and Begoña Selection 2017 <pennalima@gmail.com>
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
 * NcClusterMassSelection:
 *
 * Cluster mass distribution model based on Selection et al.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_cluster_mass_selection.h"
#include "math/ncm_integrate.h"
#include "math/ncm_spline2d_bicubic.h"
#include "math/ncm_spline2d.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_c.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */


#define _NC_CLUSTER_MASS_SELECTION_DEFAULT_INT_KEY 6

struct _NcClusterMassSelectionPrivate
{
	gdouble M0;
	gdouble z0;
	gdouble lnM0;
	gdouble ln1pz0;
	gdouble lnR_max;
	gdouble lnR_min;
	gboolean enable_rejection;
	NcmSpline2d *purity;
  NcmSpline2d *completeness;
};



G_DEFINE_TYPE_WITH_PRIVATE (NcClusterMassSelection, nc_cluster_mass_selection, NC_TYPE_CLUSTER_MASS)

#define VECTOR   (NCM_MODEL (selection))
#define MU_P0    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_SELECTION_MU_P0))
#define MU_P1    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_SELECTION_MU_P1))
#define MU_P2    (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_SELECTION_MU_P2))
#define SIGMA_P0 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_SELECTION_SIGMA_P0))
#define SIGMA_P1 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_SELECTION_SIGMA_P1))
#define SIGMA_P2 (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_SELECTION_SIGMA_P2))
#define LNM_C0   (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_SELECTION_LNM_C0))
#define LNM_CZ   (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_SELECTION_LNM_CZ))
#define A_C0     (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_SELECTION_A_C0))
#define A_CZ     (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_SELECTION_A_CZ))
#define CUT      (ncm_model_orig_param_get (VECTOR, NC_CLUSTER_MASS_SELECTION_CUT))


enum
{
	PROP_0,
	PROP_M0,
	PROP_Z0,
	PROP_LNRICHNESS_MIN,
	PROP_LNRICHNESS_MAX,
	PROP_ENABLE_REJECTION,
	PROP_PURITY,
  PROP_COMPLETENESS,
	PROP_SIZE,
};

static void
nc_cluster_mass_selection_init (NcClusterMassSelection *selection)
{
	NcClusterMassSelectionPrivate * const self = selection->priv = nc_cluster_mass_selection_get_instance_private (selection);

	self->M0               = 0.0;
	self->z0               = 0.0;
	self->lnM0             = 0.0;
	self->ln1pz0           = 0.0;
	self->lnR_min          = GSL_NEGINF;
	self->lnR_max          = GSL_POSINF;
	self->enable_rejection = TRUE;
	self->purity           = NULL;
    self->completeness     = NULL;
}



static void
_nc_cluster_mass_selection_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
	NcClusterMassSelection *selection             = NC_CLUSTER_MASS_SELECTION (object);
	NcClusterMassSelectionPrivate * const self = selection->priv;
    NcmSpline2d *s2d_purity                    = NCM_SPLINE2D (self->purity);
    NcmSpline2d *s2d_completeness              = NCM_SPLINE2D (self->completeness);

	g_return_if_fail (NC_IS_CLUSTER_MASS_SELECTION (object));

	switch (prop_id)
	{
	case PROP_M0:
		self->M0   = g_value_get_double (value);
		self->lnM0 = log (self->M0);
		break;
	case PROP_Z0:
		self->z0     = g_value_get_double (value);
		self->ln1pz0 = log1p (self->z0);
		break;
	case PROP_LNRICHNESS_MIN:
		self->lnR_min = g_value_get_double (value);
		g_assert (self->lnR_min < self->lnR_max);
		break;
	case PROP_LNRICHNESS_MAX:
		self->lnR_max = g_value_get_double (value);
		g_assert (self->lnR_min < self->lnR_max);
		break;
	case PROP_ENABLE_REJECTION:
		nc_cluster_mass_selection_set_enable_rejection (selection, g_value_get_boolean (value));
		break;
	case PROP_PURITY:
		ncm_spline2d_clear (&s2d_purity);
		s2d_purity = g_value_dup_object (value);
		break;
    case PROP_COMPLETENESS:
		ncm_spline2d_clear (&s2d_completeness);
		s2d_completeness  = g_value_dup_object (value);
		break;
	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
_nc_cluster_mass_selection_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
	NcClusterMassSelection *selection             = NC_CLUSTER_MASS_SELECTION (object);
	NcClusterMassSelectionPrivate * const self = selection->priv;

	g_return_if_fail (NC_IS_CLUSTER_MASS_SELECTION (object));

	switch (prop_id)
	{
	case PROP_M0:
		g_value_set_double (value, self->M0);
		break;
	case PROP_Z0:
		g_value_set_double (value, self->z0);
		break;
	case PROP_LNRICHNESS_MIN:
		g_value_set_double (value, self->lnR_min);
		break;
	case PROP_LNRICHNESS_MAX:
		g_value_set_double (value, self->lnR_max);
		break;
	case PROP_ENABLE_REJECTION:
		g_value_set_boolean (value, self->enable_rejection);
		break;
	case PROP_PURITY:
		g_value_set_object (value, self->purity);
		break;
    case PROP_COMPLETENESS:
		g_value_set_object (value, self->completeness);
		break;
	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
_nc_cluster_mass_selection_finalize (GObject *object)
{
	/* Chain up : end */
	G_OBJECT_CLASS (nc_cluster_mass_selection_parent_class)->finalize (object);
}

static void _nc_cluster_mass_selection_completeness(NcClusterMass *clusterm,gdouble lnM, gdouble z, gdouble *completeness);
static void _nc_cluster_mass_selection_purity(NcClusterMass *clusterm,gdouble lnM_obs, gdouble z, gdouble *purity);
static gdouble _nc_cluster_mass_selection_p (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params);
static gdouble _nc_cluster_mass_selection_intp (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z);
static gdouble _nc_cluster_mass_selection_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params);
static gboolean _nc_cluster_mass_selection_resample (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng);
static void _nc_cluster_mass_selection_p_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_selection_p_bin_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper);
static void _nc_cluster_mass_selection_n_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper);
static gdouble _nc_cluster_mass_selection_volume (NcClusterMass *clusterm);

static void _nc_cluster_mass_selection_p_vec_z_lnMobs (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const NcmVector *z, const NcmMatrix *lnM_obs, const NcmMatrix *lnM_obs_params, GArray *res);

static void
nc_cluster_mass_selection_class_init (NcClusterMassSelectionClass *klass)
{
	GObjectClass *object_class       = G_OBJECT_CLASS (klass);
	NcClusterMassClass *parent_class = NC_CLUSTER_MASS_CLASS (klass);
	NcmModelClass *model_class       = NCM_MODEL_CLASS (klass);

	model_class->set_property = &_nc_cluster_mass_selection_set_property;
	model_class->get_property = &_nc_cluster_mass_selection_get_property;
	object_class->finalize    = &_nc_cluster_mass_selection_finalize;

	ncm_model_class_set_name_nick (model_class, "Selection Ln-normal richness distribution with selection function", "Selection");
	ncm_model_class_add_params (model_class, NC_CLUSTER_MASS_SELECTION_SPARAM_LEN, 0, PROP_SIZE);

	/**
	 * NcClusterMassSelection:M0:
	 *
	 * Pivot mass FIXME Set correct values (limits)
	 */
	g_object_class_install_property (object_class,
	                                 PROP_M0,
	                                 g_param_spec_double ("M0",
	                                                      NULL,
	                                                      "Pivot mass",
	                                                      11.0 * M_LN10, G_MAXDOUBLE, 3.0e14 / 0.71,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

/**
 * NcClusterMassSelection:Z0:
 *
 * Pivot redshift FIXME Set correct values (limits)
 */
	g_object_class_install_property (object_class,
	                                 PROP_Z0,
	                                 g_param_spec_double ("z0",
	                                                      NULL,
	                                                      "Pivot redshift",
	                                                      0.0, G_MAXDOUBLE, 0.6,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));



	/**
	 * NcClusterMassSelection:lnRichness_min:
	 *
	 * FIXME Set correct values (limits)
	 */
	g_object_class_install_property (object_class,
	                                 PROP_LNRICHNESS_MIN,
	                                 g_param_spec_double ("lnRichness-min",
	                                                      NULL,
	                                                      "Minimum LnRichness",
	                                                      0.0, G_MAXDOUBLE, M_LN10 * 1.0,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	/**
	 * NcClusterMassSelection:lnRichness_max:
	 *
	 * FIXME Set correct values (limits)intp
	 */
	g_object_class_install_property (object_class,
	                                 PROP_LNRICHNESS_MAX,
	                                 g_param_spec_double ("lnRichness-max",
	                                                      NULL,
	                                                      "Maximum LnRichness",
	                                                      0.0, G_MAXDOUBLE,  M_LN10 * 2.5,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	/**
	 * NcClusterMassSelection:enable_rejection:
	 *
	 * FIXME Set if the objects sampled below CUT are rejected
	 */
	g_object_class_install_property (object_class,
	                                 PROP_ENABLE_REJECTION,
	                                 g_param_spec_boolean ("enable-rejection",
	                                                       NULL,
	                                                       "Whether rejects the sampled objects below the CUT",
	                                                       TRUE,
	                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


	/**
	 * NcClusterMassSelection:MU_P0:
	 *
	 * Distribution's  bias in the mean.
	 * FIXME Set correct values (limits)
	 */
	ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_SELECTION_MU_P0, "mu_p0", "mup0",
	                            0.0,  6.0, 1.0e-1,
	                            NC_CLUSTER_MASS_SELECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_SELECTION_DEFAULT_MU_P0,
	                            NCM_PARAM_TYPE_FIXED);

	/**
	 * NcClusterMassSelection:MU_P1:
	 *
	 * Distribution's slope with respect to the mass in the mean.
	 * FIXME Set correct values (limits)
	 */
	ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_SELECTION_MU_P1, "mu_p1", "mup1",
	                            -10.0,  10.0, 1.0e-2,
	                            NC_CLUSTER_MASS_SELECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_SELECTION_DEFAULT_MU_P1,
	                            NCM_PARAM_TYPE_FIXED);

	/**
	 * NcClusterMassSelection:MU_P2:
	 *
	 * Distribution's slope with respect to the redshift in the mean.
	 * FIXME Set correct values (limits)
	 */
	ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_SELECTION_MU_P2, "mu_p2", "mup2",
	                            -10.0,  10.0, 1.0e-2,
	                            NC_CLUSTER_MASS_SELECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_SELECTION_DEFAULT_MU_P2,
	                            NCM_PARAM_TYPE_FIXED);

	/**
	 * NcClusterMassSelection:sigma_P0:
	 *
	 * Distribution's bias in the standard deviation, $\sigma \in [10^{-4}, 10]$.
	 *
	 * FIXME Set correct values (limits)
	 */
	ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_SELECTION_SIGMA_P0, "\\sigma_p0", "sigmap0",
	                            1.0e-4, 10.0, 1.0e-2,
	                            NC_CLUSTER_MASS_SELECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_SELECTION_DEFAULT_SIGMA_P0,
	                            NCM_PARAM_TYPE_FIXED);

	/**
	 * NcClusterMassSelection:sigma_P1:
	 *
	 * Distribution's slope with respect to the mass in the standard deviation.
	 *
	 * FIXME Set correct values (limits)
	 */
	ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_SELECTION_SIGMA_P1, "\\sigma_p1", "sigmap1",
	                            -10.0, 10.0, 1.0e-2,
	                            NC_CLUSTER_MASS_SELECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_SELECTION_DEFAULT_SIGMA_P1,
	                            NCM_PARAM_TYPE_FIXED);


/**
 * NcClusterMassSelection:sigma_P2:
 *
 * Distribution's slope with respect to the redshift in the standard deviation.
 *
 * FIXME Set correct values (limits)
 */
	ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_SELECTION_SIGMA_P2, "\\sigma_p2", "sigmap2",
	                            -10.0,  10.0, 1.0e-2,
	                            NC_CLUSTER_MASS_SELECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_SELECTION_DEFAULT_SIGMA_P2,
	                            NCM_PARAM_TYPE_FIXED);




/**
 * NcClusterMassSelection:CUT:
 *
 * Cut in richness.
 *
 */
	ncm_model_class_set_sparam (model_class, NC_CLUSTER_MASS_SELECTION_CUT, "CUT", "cut",
	                            0.0,  1.0e16, 1.0e-2,
	                            NC_CLUSTER_MASS_SELECTION_DEFAULT_PARAMS_ABSTOL, NC_CLUSTER_MASS_SELECTION_DEFAULT_CUT,
	                            NCM_PARAM_TYPE_FIXED);


	/* Check for errors in parameters initialization */
	ncm_model_class_check_params_info (model_class);

	parent_class->P               = &_nc_cluster_mass_selection_p;
	parent_class->intP            = &_nc_cluster_mass_selection_intp;
	parent_class->intP_bin        = &_nc_cluster_mass_selection_intp_bin;
	parent_class->resample        = &_nc_cluster_mass_selection_resample;
	parent_class->P_limits        = &_nc_cluster_mass_selection_p_limits;
	parent_class->P_bin_limits    = &_nc_cluster_mass_selection_p_bin_limits;
	parent_class->N_limits        = &_nc_cluster_mass_selection_n_limits;
	parent_class->volume          = &_nc_cluster_mass_selection_volume;
	parent_class->P_vec_z_lnMobs  = &_nc_cluster_mass_selection_p_vec_z_lnMobs;
	parent_class->_obs_len        = 1;
	parent_class->_obs_params_len = 0;

	ncm_model_class_add_impl_flag (model_class, NC_CLUSTER_MASS_IMPL_ALL);
}

static void
_nc_cluster_mass_selection_lnR_sigma (NcClusterMass *clusterm, const gdouble lnM, const gdouble z, gdouble *lnR, gdouble *sigma)
{
	NcClusterMassSelection *selection             = NC_CLUSTER_MASS_SELECTION (clusterm);
	NcClusterMassSelectionPrivate * const self = selection->priv;
	const gdouble DlnM                      = lnM - self->lnM0;
	const gdouble Dln1pz                    = log1p (z) - self->ln1pz0;

	lnR[0]   = MU_P0    + MU_P1    * DlnM + MU_P2    * Dln1pz;
	sigma[0] = SIGMA_P0 + SIGMA_P1 * DlnM + SIGMA_P2 * Dln1pz;
}

void
nc_cluster_mass_selection_set_completeness(NcClusterMassSelection *selection, NcmSpline2dBicubic *completeness)
{

	NcClusterMassSelectionPrivate * const self = selection->priv;
	self->completeness = NCM_SPLINE2D(completeness);
	ncm_spline2d_prepare(self->completeness);
}

NcmSpline2d *
nc_cluster_mass_selection_get_completeness(NcClusterMassSelection *selection)
{
	NcClusterMassSelectionPrivate * const self = selection->priv;

	return self->completeness;
}

static void
_nc_cluster_mass_selection_completeness(NcClusterMass *clusterm, gdouble lnM, gdouble z, gdouble *completeness)
{
	NcClusterMassSelection *selection   = NC_CLUSTER_MASS_SELECTION (clusterm);
    NcClusterMassSelectionPrivate * const self = selection->priv;

    NcmSpline2d *s2d                    = NCM_SPLINE2D (self->completeness);

	completeness[0] =  ncm_spline2d_eval(s2d, lnM, z);


}


void
nc_cluster_mass_selection_set_purity(NcClusterMassSelection *selection, NcmSpline2dBicubic *purity)
{

	NcClusterMassSelectionPrivate * const self = selection->priv;
	self->purity = NCM_SPLINE2D(purity);
	ncm_spline2d_prepare(self->purity);
}

NcmSpline2d *
nc_cluster_mass_selection_get_purity(NcClusterMassSelection *selection)
{
	NcClusterMassSelectionPrivate * const self = selection->priv;

	return self->purity;
}

static void
_nc_cluster_mass_selection_purity(NcClusterMass *clusterm, gdouble lnM_obs, gdouble z,  gdouble *purity)
{
	NcClusterMassSelection *selection   = NC_CLUSTER_MASS_SELECTION (clusterm);
    NcClusterMassSelectionPrivate * const self = selection->priv;

    NcmSpline2d *s2d                    = NCM_SPLINE2D (self->purity);

	purity[0] =  ncm_spline2d_eval(s2d, lnM_obs, z);
}





static gdouble
_nc_cluster_mass_selection_p (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs, const gdouble *lnM_obs_params)
{
	NcClusterMassSelection *selection = NC_CLUSTER_MASS_SELECTION (clusterm);

	gdouble lnR_true, sigma, completeness, ipurity;

	_nc_cluster_mass_selection_lnR_sigma (clusterm, lnM, z, &lnR_true, &sigma);
	_nc_cluster_mass_selection_completeness(clusterm, lnM, z, &completeness);
	_nc_cluster_mass_selection_purity(clusterm, lnM_obs[0], z, &ipurity);

	const gdouble x     = (lnM_obs[0] - lnR_true) / sigma;

    if (lnM_obs[0] < CUT)
        return 0.0;
    else
        return 2.0 / (ncm_c_sqrt_2pi () * sigma) * exp (-0.5 * x * x)/erfc ((CUT - lnR_true) / (M_SQRT2 * sigma)) * completeness *ipurity;;
}


typedef struct _NcClusterMassSelectionInt
{
	NcHICosmo *cosmo;
	NcClusterMass *clusterm;
	gdouble lnM_obs;
	gdouble lnM;
	gdouble z;
	const gdouble *lnM_obs_params;
} NcClusterMassSelectionInt;

static gdouble
_nc_cluster_mass_selection_integrand (gdouble lnM_obs, gpointer userdata)
{
	NcClusterMassSelectionInt *obs_data   = (NcClusterMassSelectionInt *) userdata;
	gdouble lnR_true, sigma, ipurity;
	NcClusterMassSelection *selection             = NC_CLUSTER_MASS_SELECTION (obs_data->clusterm);

	_nc_cluster_mass_selection_lnR_sigma (obs_data->clusterm, obs_data->lnM, obs_data->z, &lnR_true, &sigma);
	_nc_cluster_mass_selection_purity(obs_data->clusterm, lnM_obs, obs_data->z, &ipurity);

	const gdouble x     = (lnM_obs - lnR_true) / sigma;
    if (lnM_obs< CUT)
    {return 0.0; }
    else
    {return 2.0 / (ncm_c_sqrt_2pi () * sigma) * exp (-0.5 * x * x)/erfc ((CUT - lnR_true) / (M_SQRT2 * sigma)) * ipurity;}

}

static gdouble
_nc_cluster_mass_selection_intp (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z)
{
	NcClusterMassSelection *selection = NC_CLUSTER_MASS_SELECTION (clusterm);
	NcClusterMassSelectionPrivate * const self = selection->priv;

	NcClusterMassSelectionInt obs_data;
	gdouble intp, err, completeness;
	gsl_function F;
	gsl_integration_workspace **w = ncm_integral_get_workspace ();

	obs_data.clusterm       = clusterm;
	obs_data.cosmo          = cosmo;
	obs_data.lnM            = lnM;
	obs_data.z              = z;

	F.function = &_nc_cluster_mass_selection_integrand;
	F.params   = &obs_data;

	gsl_integration_qag (&F, CUT, self->lnR_max, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_MASS_SELECTION_DEFAULT_INT_KEY, *w, &intp, &err);
	ncm_memory_pool_return (w);
	_nc_cluster_mass_selection_completeness(clusterm, lnM, z, &completeness);
    return intp * completeness;
}

static gdouble
_nc_cluster_mass_selection_intp_bin (NcClusterMass *clusterm, NcHICosmo *cosmo, gdouble lnM, gdouble z, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params)
{
	NcClusterMassSelection *selection = NC_CLUSTER_MASS_SELECTION (clusterm);

	if ((lnM_obs_lower[0]< CUT) && (lnM_obs_upper[0]< CUT))
	{
		return 0.0;
	}
	else
	{
		NcClusterMassSelectionInt obs_data;
		gdouble intp_bin, err, completeness;
		gsl_function F;
		gsl_integration_workspace **w = ncm_integral_get_workspace ();

		obs_data.clusterm       = clusterm;
		obs_data.cosmo          = cosmo;
		obs_data.lnM            = lnM;
		obs_data.lnM_obs_params = lnM_obs_params;
		obs_data.z              = z;
		_nc_cluster_mass_selection_completeness(clusterm, lnM, z, &completeness);

		F.function = &_nc_cluster_mass_selection_integrand;
		F.params   = &obs_data;

		if ((lnM_obs_lower[0]< CUT) && (lnM_obs_upper[0] >= CUT))
		{
			gsl_integration_qag (&F, CUT, lnM_obs_upper[0], 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_MASS_SELECTION_DEFAULT_INT_KEY, *w, &intp_bin, &err);
			ncm_memory_pool_return (w);
		}
		else
		{
			gsl_integration_qag (&F, lnM_obs_lower[0], lnM_obs_upper[0], 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, _NC_CLUSTER_MASS_SELECTION_DEFAULT_INT_KEY, *w, &intp_bin, &err);
			ncm_memory_pool_return (w);
		}
		if (intp_bin * completeness < 0.0)
        return 0.0;
        else
    	return intp_bin * completeness;
	}
}


static gboolean
_nc_cluster_mass_selection_resample (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble lnM, gdouble z, gdouble *lnM_obs, const gdouble *lnM_obs_params, NcmRNG *rng)
{
	NcClusterMassSelection *selection             = NC_CLUSTER_MASS_SELECTION (clusterm);
	NcClusterMassSelectionPrivate * const self = selection->priv;
	gdouble lnR_true, sigma;

	_nc_cluster_mass_selection_lnR_sigma (clusterm, lnM, z, &lnR_true, &sigma);

	ncm_rng_lock (rng);

	if (self->enable_rejection)
	{
		lnM_obs[0] = ncm_rng_gaussian_gen (rng, lnR_true, sigma);
	}
	else
	{
		lnM_obs[0]  = ncm_rng_gaussian_tail_gen (rng, CUT - lnR_true, sigma);
		lnM_obs[0] += lnR_true;
	}

	ncm_rng_unlock (rng);

	return (lnM_obs[0] <= self->lnR_max) && (lnM_obs[0] >= CUT);
}

static void
_nc_cluster_mass_selection_p_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, const gdouble *lnM_obs, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{

	const gdouble lnMl =  M_LN10 * 13.0;
	const gdouble lnMu =  M_LN10 * 16.0;

	*lnM_lower = lnMl;
	*lnM_upper = lnMu;
	return;
}

static void
_nc_cluster_mass_selection_p_bin_limits (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble *lnM_obs_lower, const gdouble *lnM_obs_upper, const gdouble *lnM_obs_params, gdouble *lnM_lower, gdouble *lnM_upper)
{

	const gdouble lnMl =  M_LN10 * 13.0;
	const gdouble lnMu =  M_LN10 * 16.0;

	*lnM_lower = lnMl;
	*lnM_upper = lnMu;
}

static void
_nc_cluster_mass_selection_n_limits (NcClusterMass *clusterm,  NcHICosmo *cosmo, gdouble *lnM_lower, gdouble *lnM_upper)
{

	const gdouble lnMl =  M_LN10 * 13.0;
	const gdouble lnMu =  M_LN10 * 16.0;

	*lnM_lower = lnMl;
	*lnM_upper = lnMu;

	return;
}

static gdouble
_nc_cluster_mass_selection_volume (NcClusterMass *clusterm)
{
	NcClusterMassSelection *selection             = NC_CLUSTER_MASS_SELECTION (clusterm);
	NcClusterMassSelectionPrivate * const self = selection->priv;

	return (self->lnR_max - CUT);
}

static void
_nc_cluster_mass_selection_p_vec_z_lnMobs (NcClusterMass *clusterm, NcHICosmo *cosmo, const gdouble lnM, const NcmVector *z, const NcmMatrix *lnM_obs, const NcmMatrix *lnM_obs_params, GArray *res)
{
	NcClusterMassSelection *selection             = NC_CLUSTER_MASS_SELECTION (clusterm);
	NcClusterMassSelectionPrivate * const self = selection->priv;

	const gdouble *lnM_obs_ptr = ncm_matrix_const_data (lnM_obs);
	const gdouble *z_ptr       = ncm_vector_const_data (z);
	const guint tda            = ncm_matrix_tda (lnM_obs);
	const guint sz             = ncm_vector_stride (z);
	const guint len            = ncm_vector_len (z);
	const gdouble DlnM         = lnM - self->lnM0;
	const gdouble sqrt_2pi     = ncm_c_sqrt_2pi ();
	const gdouble lnR_pre      = MU_P0    + MU_P1    * DlnM;
	const gdouble sigma_pre    = SIGMA_P0 + SIGMA_P1 * DlnM;
	const gdouble mu_p2        = MU_P2;
	const gdouble sigma_p2     = SIGMA_P2;
	gdouble *res_ptr           = &g_array_index (res, gdouble, 0);
	guint i;

	if ((tda == 1) && (sz == 1))
	{
		for (i = 0; i < len; i++)
		{
			const gdouble Dln1pz = log1p (z_ptr[i]) - self->ln1pz0;
			const gdouble lnR    = lnR_pre + mu_p2 * Dln1pz;
			const gdouble sigma  = sigma_pre + sigma_p2 * Dln1pz;
			const gdouble x      = (lnM_obs_ptr[i] - lnR) / sigma;

			if (lnM_obs_ptr[i] <0.0)
			{
				res_ptr[i]=0.0;
			}
			else
			{
				res_ptr[i] = 1.0 * exp (-0.5 * x * x) / (sqrt_2pi * sigma);
			}
		}
	}

	else
	{
		for (i = 0; i < len; i++)
		{
			const gdouble Dln1pz = log1p (z_ptr[i * sz]) - self->ln1pz0;
			const gdouble lnR    = lnR_pre + mu_p2 * Dln1pz;
			const gdouble sigma  = sigma_pre + sigma_p2 * Dln1pz;
			const gdouble x      = (lnM_obs_ptr[i * tda] - lnR) / sigma;

			if (lnM_obs_ptr[i * tda] <0.0)
			{
				res_ptr[i]=0.0;
			}
			else
			{
				res_ptr[i] = 1.0 * exp (-0.5 * x * x) / (sqrt_2pi * sigma);
			}
		}
	}
}

/**
 * nc_cluster_mass_selection_get_mean_richness:
 * @selection: a #NcClusterMassSelection
 * @lnM: ln of the mass
 * @z: redshift
 *
 * Computes the mean of the richness distribution.
 *
 */
gdouble
nc_cluster_mass_selection_get_mean_richness (NcClusterMassSelection *selection, gdouble lnM, gdouble z)
{
	NcClusterMassSelectionPrivate * const self = selection->priv;
	const gdouble DlnM                      = lnM - self->lnM0;
	const gdouble Dln1pz                    = log1p (z) - self->ln1pz0;

	return MU_P0    + MU_P1    * DlnM + MU_P2    * Dln1pz;
}

/**
 * nc_cluster_mass_selection_get_std_richness:
 * @selection: a #NcClusterMassSelection
 * @lnM: ln of the mass
 * @z: redshift
 *
 * Computes the standard deviation of the richness distribution.
 *
 */
gdouble
nc_cluster_mass_selection_get_std_richness (NcClusterMassSelection *selection, gdouble lnM, gdouble z)
{
	NcClusterMassSelectionPrivate * const self = selection->priv;
	const gdouble DlnM                      = lnM - self->lnM0;
	const gdouble Dln1pz                    = log1p (z) - self->ln1pz0;

	return SIGMA_P0 + SIGMA_P1 * DlnM + SIGMA_P2 * Dln1pz;
}

/**
 * nc_cluster_mass_selection_get_cut:
 * @selection: a #NcClusterMassSelection
 * @lnM: ln of the mass
 * @z: redshift
 *
 * Computes the cut in richness.
 *
 * Returns: the cut in richness.
 */
gdouble
nc_cluster_mass_selection_get_cut (NcClusterMassSelection *selection, gdouble lnM, gdouble z)
{
	return CUT;
}

/**
 * nc_cluster_mass_selection_get_mean:
 * @selection: a #NcClusterMassSelection
 * @lnM: ln of the mass
 * @z: redshift
 *
 * Computes the mean of the richness distribution with the cut correction.
 *
 */
gdouble
nc_cluster_mass_selection_get_mean (NcClusterMassSelection *selection, gdouble lnM, gdouble z)
{
	gdouble lnR_mean, lnR_sigma, A, B, C, mean_correction;

	lnR_mean  = nc_cluster_mass_selection_get_mean_richness (selection, lnM, z);
	lnR_sigma = nc_cluster_mass_selection_get_std_richness  (selection, lnM, z);

	A = (CUT - lnR_mean) / lnR_sigma;

	B = (1.0 / (ncm_c_sqrt_2pi ())) * exp (-0.5 * (A  * A));

	C = 1.0 - 0.5 * (1.0 + erf (A / M_SQRT2));

	mean_correction = (lnR_sigma * B / C);

	return lnR_mean + mean_correction;
}

/**
 * nc_cluster_mass_selection_get_std:
 * @selection: a #NcClusterMassSelection
 * @lnM: ln of the mass
 * @z: redshift
 *
 * Computes the standard deviation of the richness distribution with the cut correction.
 *
 */
gdouble
nc_cluster_mass_selection_get_std (NcClusterMassSelection *selection, gdouble lnM, gdouble z)
{
	gdouble lnR_mean, lnR_sigma, A, B, C, std_correction;


	lnR_mean  = nc_cluster_mass_selection_get_mean_richness (selection, lnM, z);
	lnR_sigma = nc_cluster_mass_selection_get_std_richness  (selection, lnM, z);

	A = (CUT - lnR_mean) / lnR_sigma;

	B = (1.0 / (ncm_c_sqrt_2pi ())) * exp (-0.5 * (A  * A));

	C = 1.0 - 0.5 * (1.0 + erf (A / M_SQRT2));

	std_correction = pow (1.0 + (A * B / C) - (B / C) * (B / C), 0.5);

	return lnR_sigma * std_correction;
}

/**
 * nc_cluster_mass_selection_set_enable_rejection:
 * @selection: a #NcClusterMassSelection
 * @on: a #NcClusterMassSelection
 *
 * Set the enable_rejection property.
 *
 */

void
nc_cluster_mass_selection_set_enable_rejection (NcClusterMassSelection *selection, gboolean on)
{
	NcClusterMassSelectionPrivate * const self = selection->priv;

	self->enable_rejection = on;
}

/**
 * nc_cluster_mass_selection_get_enable_rejection:
 * @selection: a #NcClusterMassSelection
 *
 * Get if the enable_rejection property is on.
 *
 */
gboolean
nc_cluster_mass_selection_get_enable_rejection (NcClusterMassSelection *selection)
{
	NcClusterMassSelectionPrivate * const self = selection->priv;

	return self->enable_rejection;
}
