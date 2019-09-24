/***************************************************************************
 *            nc_cbe.c
 *
 *  Sat October 24 11:56:56 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_cbe.c
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
 * SECTION:nc_cbe
 * @title: NcCBE
 * @short_description: CLASS (Cosmic Linear Anisotropy Solving System) backend
 *
 * This object provides an interface for the CLASS code.
 *
 * If you use this object please cite: [Blas (2011) CLASS II][XBlas2011],
 * see also:
 * - [Lesgourgues (2011) CLASS I][XLesgourgues2011],
 * - [Lesgourgues (2011) CLASS III][XLesgourgues2011a],
 * - [Lesgourgues (2011) CLASS IV][XLesgourgues2011b] and
 * - [CLASS website](http://class-code.net/).
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

/*
 * It must be include before anything else, several symbols clash
 * with the default includes.
 */
#ifndef NUMCOSMO_GIR_SCAN
#include "class/include/class.h"
#endif /* NUMCOSMO_GIR_SCAN */

#include "build_cfg.h"

#include "math/ncm_spline2d_bicubic.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "model/nc_hicosmo_de.h"
#include "model/nc_hicosmo_de_xcdm.h"
#include "model/nc_hicosmo_de_cpl.h"
#include "nc_cbe.h"
#include "nc_enum_types.h"
#include "nc_hiprim.h"
#include "nc_hireion_camb.h"

enum
{
	PROP_0,
	PROP_PREC,
	PROP_TARGET_CLS,
	PROP_CALC_TRANSFER,
	PROP_USE_LENSED_CLS,
	PROP_USE_TENSOR,
	PROP_USE_THERMODYN,
	PROP_SCALAR_LMAX,
	PROP_VECTOR_LMAX,
	PROP_TENSOR_LMAX,
	PROP_MATTER_PK_MAXZ,
	PROP_MATTER_PK_MAXK,
  PROP_USE_PPF,
  PROP_VERBOSE
};

struct _NcCBEPrivate
{
	struct background pba;
	struct thermo pth;
	struct perturbs ppt;
	struct transfers ptr;
	struct primordial ppm;
	struct spectra psp;
	struct nonlinear pnl;
	struct lensing ple;
	struct output pop;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcCBE, nc_cbe, G_TYPE_OBJECT);

static void
nc_cbe_init (NcCBE *cbe)
{
	cbe->priv               = nc_cbe_get_instance_private (cbe);
	cbe->prec               = NULL;
	cbe->ctrl_cosmo         = ncm_model_ctrl_new (NULL);
	cbe->ctrl_prim          = ncm_model_ctrl_new (NULL);
	cbe->a                  = NULL;

	cbe->target_Cls         = 0;
	cbe->calc_transfer      = FALSE;
	cbe->use_lensed_Cls     = FALSE;
	cbe->use_tensor         = FALSE;
	cbe->scalar_lmax        = 0;
	cbe->vector_lmax        = 0;
	cbe->tensor_lmax        = 0;

	cbe->call               = NULL;
	cbe->free               = NULL;
	cbe->allocated          = FALSE;
	cbe->thermodyn_prepared = FALSE;

	/* background structure */
  
	cbe->priv->pba.h                   = 0.0;
	cbe->priv->pba.H0                  = 0.0;
	cbe->priv->pba.T_cmb               = 0.0;
	cbe->priv->pba.Omega0_g            = 0.0;
	cbe->priv->pba.Omega0_ur           = 0.0;
	cbe->priv->pba.Omega0_b            = 0.0;
	cbe->priv->pba.Omega0_cdm          = 0.0;
	cbe->priv->pba.Omega0_dcdmdr       = 0.0;
	cbe->priv->pba.Omega0_dcdm         = 0.0;
	cbe->priv->pba.Gamma_dcdm          = 0.0;
	cbe->priv->pba.N_ncdm              = 0;
	cbe->priv->pba.Omega0_ncdm_tot     = 0.0;
	cbe->priv->pba.ksi_ncdm_default    = 0.0;
	cbe->priv->pba.ksi_ncdm            = NULL;
	cbe->priv->pba.T_ncdm_default      = 0.0;
	cbe->priv->pba.T_ncdm              = NULL;
	cbe->priv->pba.deg_ncdm_default    = 0.0;
	cbe->priv->pba.deg_ncdm            = NULL;
	cbe->priv->pba.ncdm_psd_parameters = NULL;
	cbe->priv->pba.ncdm_psd_files      = NULL;
	cbe->priv->pba.Omega0_scf          = 0.0;
	cbe->priv->pba.attractor_ic_scf    = _FALSE_;
	cbe->priv->pba.scf_parameters      = NULL;
	cbe->priv->pba.scf_parameters_size = 0;
	cbe->priv->pba.scf_tuning_index    = 0;
	cbe->priv->pba.phi_ini_scf         = 0;
	cbe->priv->pba.phi_prime_ini_scf   = 0;
	cbe->priv->pba.Omega0_k            = 0.0;
	cbe->priv->pba.K                   = 0.0;
	cbe->priv->pba.sgnK                = 0;
	cbe->priv->pba.Omega0_lambda       = 0.0;
	cbe->priv->pba.Omega0_fld          = 0.0;
  cbe->priv->pba.fluid_equation_of_state = 0;
	cbe->priv->pba.a_today             = 0.0;
	cbe->priv->pba.w0_fld              = 0.0;
	cbe->priv->pba.wa_fld              = 0.0;
  cbe->priv->pba.Omega_EDE           = 0.0;
	cbe->priv->pba.cs2_fld             = 0.0;
  cbe->priv->pba.use_ppf             = _FALSE_;
  cbe->priv->pba.c_gamma_over_c_fld  = 0.0;
  cbe->priv->pba.shooting_failed     = _FALSE_;

  cbe->priv->pba.got_files                = NULL;
  cbe->priv->pba.m_ncdm_in_eV             = NULL;
  cbe->priv->pba.ncdm_quadrature_strategy = NULL;
  cbe->priv->pba.ncdm_input_q_size        = NULL;
  cbe->priv->pba.ncdm_qmax                = NULL;

	/* thermodynamics structure */

	cbe->priv->pth.YHe                        = 0;
	cbe->priv->pth.recombination              = 0;
	cbe->priv->pth.reio_parametrization       = 0;
	cbe->priv->pth.reio_z_or_tau              = 0;
	cbe->priv->pth.z_reio                     = 0.0;
	cbe->priv->pth.tau_reio                   = 0.0;
	cbe->priv->pth.reionization_exponent      = 0.0;
	cbe->priv->pth.reionization_width         = 0.0;
	cbe->priv->pth.helium_fullreio_redshift   = 0.0;
	cbe->priv->pth.helium_fullreio_width      = 0.0;

	cbe->priv->pth.binned_reio_num            = 0;
	cbe->priv->pth.binned_reio_z              = NULL;
	cbe->priv->pth.binned_reio_xe             = NULL;
	cbe->priv->pth.binned_reio_step_sharpness = 0.0;

	cbe->priv->pth.many_tanh_num              = 0;
	cbe->priv->pth.many_tanh_z                = NULL;
	cbe->priv->pth.many_tanh_xe               = NULL;
	cbe->priv->pth.many_tanh_width            = 0.0;

	cbe->priv->pth.reio_inter_num             = 0;
	cbe->priv->pth.reio_inter_z               = NULL;
	cbe->priv->pth.reio_inter_xe              = NULL;

	cbe->priv->pth.annihilation               = 0.0;
	cbe->priv->pth.decay                      = 0.0;
	cbe->priv->pth.annihilation_variation     = 0.0;
	cbe->priv->pth.annihilation_z             = 0.0;
	cbe->priv->pth.annihilation_zmax          = 0.0;
	cbe->priv->pth.annihilation_zmin          = 0.0;
	cbe->priv->pth.annihilation_f_halo        = 0.0;
	cbe->priv->pth.annihilation_z_halo        = 0.0;
	cbe->priv->pth.has_on_the_spot            = _FALSE_;

	cbe->priv->pth.compute_cb2_derivatives    = _FALSE_;
  cbe->priv->pth.compute_damping_scale      = _FALSE_;

	/* perturbation structure */
  
	cbe->priv->ppt.has_perturbations                   = _FALSE_;
	cbe->priv->ppt.has_cls                             = _FALSE_;

	cbe->priv->ppt.has_cl_cmb_temperature              = _FALSE_;
	cbe->priv->ppt.has_cl_cmb_polarization             = _FALSE_;
	cbe->priv->ppt.has_cl_cmb_lensing_potential        = _FALSE_;
	cbe->priv->ppt.has_cl_number_count                 = _FALSE_;
	cbe->priv->ppt.has_cl_lensing_potential            = _FALSE_;
	cbe->priv->ppt.has_pk_matter                       = _FALSE_;
	cbe->priv->ppt.has_density_transfers               = _FALSE_;
	cbe->priv->ppt.has_velocity_transfers              = _FALSE_;
  cbe->priv->ppt.has_metricpotential_transfers       = _FALSE_;

	cbe->priv->ppt.has_nl_corrections_based_on_delta_m = _FALSE_;

	cbe->priv->ppt.has_nc_density                      = _FALSE_;
	cbe->priv->ppt.has_nc_rsd                          = _FALSE_;
	cbe->priv->ppt.has_nc_lens                         = _FALSE_;
	cbe->priv->ppt.has_nc_gr                           = _FALSE_;

	cbe->priv->ppt.switch_sw                           = 0;
	cbe->priv->ppt.switch_eisw                         = 0;
	cbe->priv->ppt.switch_lisw                         = 0;
	cbe->priv->ppt.switch_dop                          = 0;
	cbe->priv->ppt.switch_pol                          = 0;
	cbe->priv->ppt.eisw_lisw_split_z                   = 0;

	cbe->priv->ppt.has_ad                              = _FALSE_;
	cbe->priv->ppt.has_bi                              = _FALSE_;
	cbe->priv->ppt.has_cdi                             = _FALSE_;
	cbe->priv->ppt.has_nid                             = _FALSE_;
	cbe->priv->ppt.has_niv                             = _FALSE_;
  
	cbe->priv->ppt.has_perturbed_recombination         = _FALSE_;
	cbe->priv->ppt.tensor_method                       = tm_massless_approximation;
	cbe->priv->ppt.evolve_tensor_ur                    = _FALSE_;
	cbe->priv->ppt.evolve_tensor_ncdm                  = _FALSE_;

	cbe->priv->ppt.has_scalars                         = _FALSE_;
	cbe->priv->ppt.has_vectors                         = _FALSE_;
	cbe->priv->ppt.has_tensors                         = _FALSE_;

	cbe->priv->ppt.l_scalar_max                        = 0;
	cbe->priv->ppt.l_vector_max                        = 0;
	cbe->priv->ppt.l_tensor_max                        = 0;
	cbe->priv->ppt.l_lss_max                           = 0;
	cbe->priv->ppt.k_max_for_pk                        = 0.0;

	cbe->priv->ppt.gauge                               = synchronous;

	cbe->priv->ppt.k_output_values_num                 = 0;
	cbe->priv->ppt.store_perturbations                 = _FALSE_;
	cbe->priv->ppt.number_of_scalar_titles             = 0;
	cbe->priv->ppt.number_of_vector_titles             = 0;
	cbe->priv->ppt.number_of_tensor_titles             = 0;

	{
		guint filenum;
		for (filenum = 0; filenum < _MAX_NUMBER_OF_K_FILES_; filenum++)
		{
			cbe->priv->ppt.scalar_perturbations_data[filenum] = NULL;
			cbe->priv->ppt.vector_perturbations_data[filenum] = NULL;
			cbe->priv->ppt.tensor_perturbations_data[filenum] = NULL;
		}
	}
  
	cbe->priv->ppt.index_k_output_values               = NULL;
  
  cbe->priv->ppt.three_ceff2_ur                      = 0.0;
  cbe->priv->ppt.three_cvis2_ur                      = 0.0;

  cbe->priv->ppt.z_max_pk                            = 0.0;

  cbe->priv->ppt.selection_num                       = 0;
  cbe->priv->ppt.selection                           = 0;
  cbe->priv->ppt.selection_mean[0]                   = 0.0;
  cbe->priv->ppt.selection_width[0]                  = 0.0;

	/* primordial structure */
 
	cbe->priv->ppm.primordial_spec_type = analytic_Pk;
	cbe->priv->ppm.k_pivot              = 0.0;
	cbe->priv->ppm.A_s                  = 0.0;
	cbe->priv->ppm.n_s                  = 0.0;
	cbe->priv->ppm.alpha_s              = 0.0;
	cbe->priv->ppm.f_bi                 = 0.0;
	cbe->priv->ppm.n_bi                 = 0.0;
	cbe->priv->ppm.alpha_bi             = 0.0;
	cbe->priv->ppm.f_cdi                = 0.0;
	cbe->priv->ppm.n_cdi                = 0.0;
	cbe->priv->ppm.alpha_cdi            = 0.0;
	cbe->priv->ppm.f_nid                = 0.0;
	cbe->priv->ppm.n_nid                = 0.0;
	cbe->priv->ppm.alpha_nid            = 0.0;
	cbe->priv->ppm.f_niv                = 0.0;
	cbe->priv->ppm.n_niv                = 0.0;
	cbe->priv->ppm.alpha_niv            = 0.0;
	cbe->priv->ppm.c_ad_bi              = 0.0;
	cbe->priv->ppm.n_ad_bi              = 0.0;
	cbe->priv->ppm.alpha_ad_bi          = 0.0;
	cbe->priv->ppm.c_ad_cdi             = 0.0;
	cbe->priv->ppm.n_ad_cdi             = 0.0;
	cbe->priv->ppm.alpha_ad_cdi         = 0.0;
	cbe->priv->ppm.c_ad_nid             = 0.0;
	cbe->priv->ppm.n_ad_nid             = 0.0;
	cbe->priv->ppm.alpha_ad_nid         = 0.0;
	cbe->priv->ppm.c_ad_niv             = 0.0;
	cbe->priv->ppm.n_ad_niv             = 0.0;
	cbe->priv->ppm.alpha_ad_niv         = 0.0;
	cbe->priv->ppm.c_bi_cdi             = 0.0;
	cbe->priv->ppm.n_bi_cdi             = 0.0;
	cbe->priv->ppm.alpha_bi_cdi         = 0.0;
	cbe->priv->ppm.c_bi_nid             = 0.0;
	cbe->priv->ppm.n_bi_nid             = 0.0;
	cbe->priv->ppm.alpha_bi_nid         = 0.0;
	cbe->priv->ppm.c_bi_niv             = 0.0;
	cbe->priv->ppm.n_bi_niv             = 0.0;
	cbe->priv->ppm.alpha_bi_niv         = 0.0;
	cbe->priv->ppm.c_cdi_nid            = 0.0;
	cbe->priv->ppm.n_cdi_nid            = 0.0;
	cbe->priv->ppm.alpha_cdi_nid        = 0.0;
	cbe->priv->ppm.c_cdi_niv            = 0.0;
	cbe->priv->ppm.n_cdi_niv            = 0.0;
	cbe->priv->ppm.alpha_cdi_niv        = 0.0;
	cbe->priv->ppm.c_nid_niv            = 0.0;
	cbe->priv->ppm.n_nid_niv            = 0.0;
	cbe->priv->ppm.alpha_nid_niv        = 0.0;

  cbe->priv->ppm.r                    = 0.0;
	cbe->priv->ppm.n_t                  = 0.0;
	cbe->priv->ppm.alpha_t              = 0.0;
	cbe->priv->ppm.potential            = 0;
	cbe->priv->ppm.phi_end              = 0.0;
  cbe->priv->ppm.phi_pivot_method     = 0;
  cbe->priv->ppm.phi_pivot_target     = 0;
	cbe->priv->ppm.V0                   = 0.0;
	cbe->priv->ppm.V1                   = 0.0;
	cbe->priv->ppm.V2                   = 0.0;
	cbe->priv->ppm.V3                   = 0.0;
	cbe->priv->ppm.V4                   = 0.0;
	cbe->priv->ppm.H0                   = 0.0;
	cbe->priv->ppm.H1                   = 0.0;
	cbe->priv->ppm.H2                   = 0.0;
	cbe->priv->ppm.H3                   = 0.0;
	cbe->priv->ppm.H4                   = 0.0;
  cbe->priv->ppm.behavior             = 0;
	cbe->priv->ppm.command              = NULL;
	cbe->priv->ppm.custom1              = 0.0;
	cbe->priv->ppm.custom2              = 0.0;
	cbe->priv->ppm.custom3              = 0.0;
	cbe->priv->ppm.custom4              = 0.0;
	cbe->priv->ppm.custom5              = 0.0;
	cbe->priv->ppm.custom6              = 0.0;
	cbe->priv->ppm.custom7              = 0.0;
	cbe->priv->ppm.custom8              = 0.0;
	cbe->priv->ppm.custom9              = 0.0;
	cbe->priv->ppm.custom10             = 0.0;

	/* transfer structure */

  cbe->priv->ptr.selection_bias[0]               = 0.0;
  cbe->priv->ptr.selection_magnification_bias[0] = 0.0;
	cbe->priv->ptr.lcmb_rescale                    = 0.0;
	cbe->priv->ptr.lcmb_pivot                      = 0.0;
	cbe->priv->ptr.lcmb_tilt                       = 0.0;
	cbe->priv->ptr.initialise_HIS_cache            = _FALSE_;
	cbe->priv->ptr.has_nz_analytic                 = _FALSE_;
	cbe->priv->ptr.has_nz_file                     = _FALSE_;
	cbe->priv->ptr.has_nz_evo_analytic             = _FALSE_;
	cbe->priv->ptr.has_nz_evo_file                 = _FALSE_;

	/* spectra structure */

	cbe->priv->psp.z_max_pk = 0.0;
	cbe->priv->psp.non_diag = 0;

	/* lensing structure */

	cbe->priv->ple.has_lensed_cls = _FALSE_;

	/* nonlinear structure */

	cbe->priv->pnl.method    = nl_none;
  cbe->priv->pnl.has_pk_eq = _FALSE_;

	/* all verbose parameters */

	cbe->priv->pba.background_verbose     = 0;
	cbe->priv->pth.thermodynamics_verbose = 0;
	cbe->priv->ppt.perturbations_verbose  = 0;
	cbe->priv->ptr.transfer_verbose       = 0;
	cbe->priv->ppm.primordial_verbose     = 0;
	cbe->priv->psp.spectra_verbose        = 0;
	cbe->priv->pnl.nonlinear_verbose      = 0;
	cbe->priv->ple.lensing_verbose        = 0;

	{
		const guint verbosity = 0;
		cbe->bg_verbose       = verbosity;
		cbe->thermo_verbose   = verbosity;
		cbe->pert_verbose     = verbosity;
		cbe->transfer_verbose = verbosity;
		cbe->prim_verbose     = verbosity;
		cbe->spectra_verbose  = verbosity;
		cbe->nonlin_verbose   = verbosity;
		cbe->lensing_verbose  = verbosity;
	}
}

static void
_nc_cbe_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
	NcCBE* cbe = NC_CBE (object);
	g_return_if_fail (NC_IS_CBE (object));

	switch (prop_id)
	{
	case PROP_PREC:
		nc_cbe_set_precision (cbe, g_value_get_object (value));
		break;
	case PROP_TARGET_CLS:
		nc_cbe_set_target_Cls (cbe, g_value_get_flags (value));
		break;
	case PROP_CALC_TRANSFER:
		nc_cbe_set_calc_transfer (cbe, g_value_get_boolean (value));
		break;
	case PROP_USE_LENSED_CLS:
		nc_cbe_set_lensed_Cls (cbe, g_value_get_boolean (value));
		break;
	case PROP_USE_TENSOR:
		nc_cbe_set_tensor (cbe, g_value_get_boolean (value));
		break;
	case PROP_USE_THERMODYN:
		nc_cbe_set_thermodyn (cbe, g_value_get_boolean (value));
		break;
	case PROP_SCALAR_LMAX:
		nc_cbe_set_scalar_lmax (cbe, g_value_get_uint (value));
		break;
	case PROP_VECTOR_LMAX:
		nc_cbe_set_vector_lmax (cbe, g_value_get_uint (value));
		break;
	case PROP_TENSOR_LMAX:
		nc_cbe_set_tensor_lmax (cbe, g_value_get_uint (value));
		break;
	case PROP_MATTER_PK_MAXZ:
		nc_cbe_set_max_matter_pk_z (cbe, g_value_get_double (value));
		break;
	case PROP_MATTER_PK_MAXK:
		nc_cbe_set_max_matter_pk_k (cbe, g_value_get_double (value));
		break;
	case PROP_USE_PPF:
		nc_cbe_use_ppf (cbe, g_value_get_boolean (value));
		break;
	case PROP_VERBOSE:
  {
    const guint verbosity = g_value_get_uint (value);
		cbe->bg_verbose       = verbosity;
		cbe->thermo_verbose   = verbosity;
		cbe->pert_verbose     = verbosity;
		cbe->transfer_verbose = verbosity;
		cbe->prim_verbose     = verbosity;
		cbe->spectra_verbose  = verbosity;
		cbe->nonlin_verbose   = verbosity;
		cbe->lensing_verbose  = verbosity;
		break;
  }
	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
_nc_cbe_get_property (GObject* object, guint prop_id, GValue* value, GParamSpec* pspec)
{
	NcCBE* cbe = NC_CBE (object);
	g_return_if_fail (NC_IS_CBE (object));

	switch (prop_id)
	{
	case PROP_PREC:
		g_value_set_object (value, nc_cbe_peek_precision (cbe));
		break;
	case PROP_TARGET_CLS:
		g_value_set_flags (value, nc_cbe_get_target_Cls (cbe));
		break;
	case PROP_CALC_TRANSFER:
		g_value_set_boolean (value, nc_cbe_calc_transfer (cbe));
		break;
	case PROP_USE_LENSED_CLS:
		g_value_set_boolean (value, nc_cbe_lensed_Cls (cbe));
		break;
	case PROP_USE_TENSOR:
		g_value_set_boolean (value, nc_cbe_tensor (cbe));
		break;
	case PROP_USE_THERMODYN:
		g_value_set_boolean (value, nc_cbe_thermodyn (cbe));
		break;
	case PROP_SCALAR_LMAX:
		g_value_set_uint (value, nc_cbe_get_scalar_lmax (cbe));
		break;
	case PROP_VECTOR_LMAX:
		g_value_set_uint (value, nc_cbe_get_vector_lmax (cbe));
		break;
	case PROP_TENSOR_LMAX:
		g_value_set_uint (value, nc_cbe_get_tensor_lmax (cbe));
		break;
	case PROP_MATTER_PK_MAXZ:
		g_value_set_double (value, nc_cbe_get_max_matter_pk_z (cbe));
		break;
	case PROP_MATTER_PK_MAXK:
		g_value_set_double (value, nc_cbe_get_max_matter_pk_k (cbe));
		break;
	case PROP_USE_PPF:
		g_value_set_boolean (value, cbe->priv->pba.use_ppf);
		break;
	case PROP_VERBOSE:
		g_value_set_uint (value, cbe->bg_verbose);
		break;
	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
_nc_cbe_dispose (GObject* object)
{
	NcCBE* cbe = NC_CBE (object);

	nc_cbe_precision_clear (&cbe->prec);
	nc_scalefactor_clear (&cbe->a);
	ncm_model_ctrl_clear (&cbe->ctrl_cosmo);
	ncm_model_ctrl_clear (&cbe->ctrl_prim);

	/* Chain up : end */
	G_OBJECT_CLASS (nc_cbe_parent_class)->dispose (object);
}

static void _nc_cbe_free_thermo (NcCBE* cbe);

static void
_nc_cbe_finalize (GObject* object)
{
	NcCBE* cbe = NC_CBE (object);

	if (cbe->allocated)
	{
		g_assert (cbe->free != NULL);
		cbe->free (cbe);
		cbe->allocated = FALSE;
	}

	if (cbe->thermodyn_prepared)
	{
		_nc_cbe_free_thermo (cbe);
		cbe->thermodyn_prepared = FALSE;
	}

	/* Chain up : end */
	G_OBJECT_CLASS (nc_cbe_parent_class)->finalize (object);
}

static void
nc_cbe_class_init (NcCBEClass* klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);

	object_class->set_property = &_nc_cbe_set_property;
	object_class->get_property = &_nc_cbe_get_property;
	object_class->dispose      = &_nc_cbe_dispose;
	object_class->finalize     = &_nc_cbe_finalize;

	g_object_class_install_property (object_class,
	                                 PROP_PREC,
	                                 g_param_spec_object ("precision",
	                                                      NULL,
	                                                      "CLASS precision object",
	                                                      NC_TYPE_CBE_PRECISION,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_TARGET_CLS,
	                                 g_param_spec_flags ("target-Cls",
	                                                     NULL,
	                                                     "Target Cls to calculate",
	                                                     NC_TYPE_DATA_CMB_DATA_TYPE, 0,
	                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_CALC_TRANSFER,
	                                 g_param_spec_boolean ("calc-transfer",
	                                                       NULL,
	                                                       "Whether to calculate the transfer function",
	                                                       FALSE,
	                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_USE_LENSED_CLS,
	                                 g_param_spec_boolean ("use-lensed-Cls",
	                                                       NULL,
	                                                       "Whether to use lensed Cls",
	                                                       FALSE,
	                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_USE_TENSOR,
	                                 g_param_spec_boolean ("use-tensor",
	                                                       NULL,
	                                                       "Whether to use tensor contributions",
	                                                       FALSE,
	                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_USE_THERMODYN,
	                                 g_param_spec_boolean ("use-thermodyn",
	                                                       NULL,
	                                                       "Whether to use the thermodynamics module",
	                                                       FALSE,
	                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_SCALAR_LMAX,
	                                 g_param_spec_uint ("scalar-lmax",
	                                                    NULL,
	                                                    "Scalar modes l_max",
	                                                    0, G_MAXUINT, 500,
	                                                    G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_VECTOR_LMAX,
	                                 g_param_spec_uint ("vector-lmax",
	                                                    NULL,
	                                                    "Vector modes l_max",
	                                                    0, G_MAXUINT, 500,
	                                                    G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_TENSOR_LMAX,
	                                 g_param_spec_uint ("tensor-lmax",
	                                                    NULL,
	                                                    "Tensor modes l_max",
	                                                    0, G_MAXUINT, 500,
	                                                    G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_MATTER_PK_MAXZ,
	                                 g_param_spec_double ("matter-pk-maxz",
	                                                      NULL,
	                                                      "Maximum redshift for matter Pk",
	                                                      0.0, G_MAXDOUBLE, 0.0,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_MATTER_PK_MAXK,
	                                 g_param_spec_double ("matter-pk-maxk",
	                                                      NULL,
	                                                      "Maximum mode k for matter Pk",
	                                                      0.0, G_MAXDOUBLE, 1.0,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_USE_PPF,
                                   g_param_spec_boolean ("use-ppf",
                                                         NULL,
                                                         "Whether to use PPF",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
	                                 PROP_VERBOSE,
                                   g_param_spec_uint ("verbosity",
                                                      NULL,
                                                      "Verbosity",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_cbe_new: (constructor)
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcCBEPrecision.
 */
NcCBE*
nc_cbe_new (void)
{
	NcCBEPrecision* prec = nc_cbe_precision_new ();
	NcCBE* cbe = g_object_new (NC_TYPE_CBE,
	                           "precision", prec,
	                           NULL);
	nc_cbe_precision_free (prec);
	return cbe;
}

/**
 * nc_cbe_prec_new: (constructor)
 * @cbe_prec: a #NcCBEPrecision.
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcCBEPrecision.
 */
NcCBE*
nc_cbe_prec_new (NcCBEPrecision* cbe_prec)
{
	NcCBE* cbe = g_object_new (NC_TYPE_CBE,
	                           "precision", cbe_prec,
	                           NULL);
	return cbe;
}

/**
 * nc_cbe_ref:
 * @cbe: a #NcCBE
 *
 * Increases the reference count of @cbe.
 *
 * Returns: (transfer full): @cbe.
 */
NcCBE*
nc_cbe_ref (NcCBE* cbe)
{
	return g_object_ref (cbe);
}

/**
 * nc_cbe_free:
 * @cbe: a #NcCBE
 *
 * Decreases the reference count of @cbe.
 *
 */
void 
nc_cbe_free (NcCBE* cbe)
{
	g_object_unref (cbe);
}

/**
 * nc_cbe_clear:
 * @cbe: a #NcCBE
 *
 * Decreases the reference count of *@cbe and sets *@cbe to NULL.
 *
 */
void 
nc_cbe_clear (NcCBE** cbe)
{
	g_clear_object (cbe);
}

/**
 * nc_cbe_set_precision:
 * @cbe: a #NcCBE
 * @cbe_prec: a #NcCBEPrecision
 *
 * Sets the @cbe_prec as the precision object.
 *
 */
void 
nc_cbe_set_precision (NcCBE* cbe, NcCBEPrecision* cbe_prec)
{
	nc_cbe_precision_clear (&cbe->prec);
	cbe->prec = nc_cbe_precision_ref (cbe_prec);
	ncm_model_ctrl_force_update (cbe->ctrl_cosmo);
}

static void _nc_cbe_update_callbacks (NcCBE* cbe);

/**
 * nc_cbe_set_target_Cls:
 * @cbe: a #NcCBE
 * @target_Cls: a #NcDataCMBDataType.
 *
 * Sets the @target_Cls target.
 *
 */
void 
nc_cbe_set_target_Cls (NcCBE* cbe, NcDataCMBDataType target_Cls)
{
	if (cbe->target_Cls != target_Cls)
	{
		cbe->target_Cls = target_Cls;
		_nc_cbe_update_callbacks (cbe);
	}
}

/**
 * nc_cbe_set_calc_transfer:
 * @cbe: a #NcCBE
 * @calc_transfer: a boolean
 *
 * Sets whether it should calculate the transfer function.
 *
 */
void 
nc_cbe_set_calc_transfer (NcCBE* cbe, gboolean calc_transfer)
{
	if ((calc_transfer && !cbe->calc_transfer) || (!calc_transfer && cbe->calc_transfer))
	{
		cbe->calc_transfer = calc_transfer;
		_nc_cbe_update_callbacks (cbe);
	}
}

/**
 * nc_cbe_set_lensed_Cls:
 * @cbe: a #NcCBE
 * @use_lensed_Cls: a boolean.
 *
 * Sets whether it should use lensed Cl's.
 *
 */
void 
nc_cbe_set_lensed_Cls (NcCBE* cbe, gboolean use_lensed_Cls)
{
	cbe->use_lensed_Cls = use_lensed_Cls;
	_nc_cbe_update_callbacks (cbe);
}

/**
 * nc_cbe_set_tensor:
 * @cbe: a #NcCBE
 * @use_tensor: a boolean
 *
 * Sets whether it should use tensor contribution.
 *
 */
void 
nc_cbe_set_tensor (NcCBE* cbe, gboolean use_tensor)
{
	cbe->use_tensor = use_tensor;
	ncm_model_ctrl_force_update (cbe->ctrl_cosmo);
}

/**
 * nc_cbe_set_thermodyn:
 * @cbe: a #NcCBE
 * @use_thermodyn: a boolean
 *
 * Sets whether it should use the thermodynamics module.
 *
 */
void 
nc_cbe_set_thermodyn (NcCBE* cbe, gboolean use_thermodyn)
{
	cbe->use_thermodyn = use_thermodyn;
	ncm_model_ctrl_force_update (cbe->ctrl_cosmo);
}

/**
 * nc_cbe_set_scalar_lmax:
 * @cbe: a #NcCBE
 * @scalar_lmax: a guint
 *
 * Sets the maximum multipole $\ell_\textrm{max}$ at which the
 * angular power spectrum $C_{\ell}$ of the scalar mode is computed.
 *
 */
void 
nc_cbe_set_scalar_lmax (NcCBE* cbe, guint scalar_lmax)
{
	if (cbe->scalar_lmax != scalar_lmax)
	{
		cbe->scalar_lmax = scalar_lmax;
		ncm_model_ctrl_force_update (cbe->ctrl_cosmo);
	}
}

/**
 * nc_cbe_set_vector_lmax:
 * @cbe: a #NcCBE
 * @vector_lmax: a guint
 *
 * Sets the maximum multipole $\ell_\textrm{max}$ at which the
 * angular power spectrum $C_{\ell}$ of the vector mode is computed.
 *
 */
void 
nc_cbe_set_vector_lmax (NcCBE* cbe, guint vector_lmax)
{
	if (cbe->vector_lmax != vector_lmax)
	{
		cbe->vector_lmax = vector_lmax;
		ncm_model_ctrl_force_update (cbe->ctrl_cosmo);
	}
}

/**
 * nc_cbe_set_tensor_lmax:
 * @cbe: a #NcCBE
 * @tensor_lmax: a guint
 *
 * Sets the maximum multipole $\ell_\textrm{max}$ at which the
 * angular power spectrum $C_{\ell}$ of the tensor mode is computed.
 *
 */
void 
nc_cbe_set_tensor_lmax (NcCBE* cbe, guint tensor_lmax)
{
	if (cbe->tensor_lmax != tensor_lmax)
	{
		cbe->tensor_lmax = tensor_lmax;
		ncm_model_ctrl_force_update (cbe->ctrl_cosmo);
	}
}

/**
 * nc_cbe_set_max_matter_pk_z:
 * @cbe: a #NcCBE
 * @zmax: maximum redshift
 *
 * Sets $z_\mathrm{max}$ for (until?) which the matter power spectrum $P(k, z)$ is evaluated.
 *
 */
void 
nc_cbe_set_max_matter_pk_z (NcCBE* cbe, gdouble zmax)
{
	if (cbe->priv->psp.z_max_pk != zmax)
	{
		cbe->priv->psp.z_max_pk = zmax;
    cbe->priv->ppt.z_max_pk = zmax;
		ncm_model_ctrl_force_update (cbe->ctrl_cosmo);
	}
}

/**
 * nc_cbe_get_max_matter_pk_z:
 * @cbe: a #NcCBE
 *
 * Gets the maximum redshift $z_\mathrm{max}$ for which the matter power spectrum $P(k, z)$ is evaluated.
 *
 * Returns: $z_\mathrm{max}$.
 */
gdouble
nc_cbe_get_max_matter_pk_z (NcCBE* cbe)
{
	return cbe->priv->psp.z_max_pk;
}

/**
 * nc_cbe_set_max_matter_pk_k:
 * @cbe: a #NcCBE
 * @kmax: maximum mode
 *
 * Sets $k_\mathrm{max}$ for which the matter power spectrum $P (k, z)$ is evaluated.
 *
 */
void 
nc_cbe_set_max_matter_pk_k (NcCBE* cbe, gdouble kmax)
{
	if (cbe->priv->ppt.k_max_for_pk != kmax)
	{
		cbe->priv->ppt.k_max_for_pk = kmax;
		ncm_model_ctrl_force_update (cbe->ctrl_cosmo);
	}
}

/**
 * nc_cbe_get_max_matter_pk_k:
 * @cbe: a #NcCBE
 *
 * Gets the maximum mode $k_\mathrm{max}$ for which the matter power spectrum $P (k, z)$ is evaluated.
 *
 * Returns: $k_\mathrm{max}$.
 */
gdouble
nc_cbe_get_max_matter_pk_k (NcCBE* cbe)
{
	return cbe->priv->ppt.k_max_for_pk;
}

/**
 * nc_cbe_peek_precision:
 * @cbe: a #NcCBE
 *
 * Peeks the #NcCBEPrecision object.
 *
 * Returns: (transfer none): the #NcCBEPrecision object.
 */
NcCBEPrecision*
nc_cbe_peek_precision (NcCBE* cbe)
{
	return cbe->prec;
}

/**
 * nc_cbe_get_target_Cls:
 * @cbe: a #NcCBE
 *
 * Gets the target_Cls flags.
 *
 * Returns: the #NcDataCMBDataType flags.
 */
NcDataCMBDataType
nc_cbe_get_target_Cls (NcCBE* cbe)
{
	return cbe->target_Cls;
}

/**
 * nc_cbe_calc_transfer:
 * @cbe: a #NcCBE
 *
 * Gets whether it calculates the transfer function.
 *
 * Returns: a boolean.
 */
gboolean
nc_cbe_calc_transfer (NcCBE* cbe)
{
	return cbe->calc_transfer;
}

/**
 * nc_cbe_lensed_Cls:
 * @cbe: a #NcCBE
 *
 * Gets whether it uses lensed $C_{\ell}$'s.
 *
 * Returns: a boolean.
 */
gboolean
nc_cbe_lensed_Cls (NcCBE* cbe)
{
	return cbe->use_lensed_Cls;
}

/**
 * nc_cbe_tensor:
 * @cbe: a #NcCBE
 *
 * Gets whether it uses tensor contributions.
 *
 * Returns: a boolean.
 */
gboolean
nc_cbe_tensor (NcCBE* cbe)
{
	return cbe->use_tensor;
}

/**
 * nc_cbe_thermodyn:
 * @cbe: a #NcCBE
 *
 * Gets whether it uses the thermodynamics module.
 *
 * Returns: a boolean.
 */
gboolean
nc_cbe_thermodyn (NcCBE* cbe)
{
	return cbe->use_thermodyn;
}

/**
 * nc_cbe_get_scalar_lmax:
 * @cbe: a #NcCBE
 *
 * Gets the maximum multipole $\ell_\textrm{max}$ at which the
 * angular power spectrum $C_{\ell}$ of the scalar mode is computed.
 *
 * Returns: the maximum (scalar) multipole $\ell_\textrm{max}$.
 */
guint nc_cbe_get_scalar_lmax (NcCBE* cbe)
{
	return cbe->scalar_lmax;
}

/**
 * nc_cbe_get_vector_lmax:
 * @cbe: a #NcCBE
 *
 * Gets the maximum multipole $\ell_\textrm{max}$ at which the
 * angular power spectrum $C_{\ell}$ of the vector mode is computed.
 *
 * Returns: the maximum (vector) multipole $\ell_\textrm{max}$.
 */
guint nc_cbe_get_vector_lmax (NcCBE* cbe)
{
	return cbe->vector_lmax;
}

/**
 * nc_cbe_get_tensor_lmax:
 * @cbe: a #NcCBE
 *
 * Gets the maximum multipole $\ell_\textrm{max}$ at which the
 * angular power spectrum $C_{\ell}$ of the tensor mode is computed.
 *
 * Returns: the maximum (tensor) multipole $\ell_\textrm{max}$.
 */
guint nc_cbe_get_tensor_lmax (NcCBE* cbe)
{
	return cbe->tensor_lmax;
}

/**
 * nc_cbe_use_ppf:
 * @cbe: a #NcCBE
 * @use_ppf: whether to use PPF
 * 
 * Sets if PPF should be used.
 * 
 */
void
nc_cbe_use_ppf (NcCBE *cbe, gboolean use_ppf)
{
  cbe->priv->pba.use_ppf = use_ppf ? _TRUE_ : _FALSE_;
}

static void
_nc_cbe_set_bg (NcCBE* cbe, NcHICosmo* cosmo)
{
	if (!g_type_is_a (G_OBJECT_TYPE (cosmo), NC_TYPE_HICOSMO_DE))
		g_error ("_nc_cbe_set_bg: CLASS backend is compatible with darkenergy models only.");

	cbe->priv->pba.h                   = nc_hicosmo_h (cosmo);
	cbe->priv->pba.H0                  = cbe->priv->pba.h * 1.0e5 / ncm_c_c ();
	cbe->priv->pba.T_cmb               = nc_hicosmo_T_gamma0 (cosmo);
	cbe->priv->pba.Omega0_g            = nc_hicosmo_Omega_g0 (cosmo);
	cbe->priv->pba.Omega0_ur           = nc_hicosmo_Omega_nu0 (cosmo);
	cbe->priv->pba.Omega0_b            = nc_hicosmo_Omega_b0 (cosmo);
	cbe->priv->pba.Omega0_cdm          = nc_hicosmo_Omega_c0 (cosmo);
	cbe->priv->pba.Omega0_dcdmdr       = 0.0;
	cbe->priv->pba.Omega0_dcdm         = 0.0;
	cbe->priv->pba.Gamma_dcdm          = 0.0;
	cbe->priv->pba.Omega0_ncdm_tot     = 0.0;
	cbe->priv->pba.ksi_ncdm_default    = 0.0;
	cbe->priv->pba.ksi_ncdm            = NULL;
	cbe->priv->pba.deg_ncdm_default    = 1.0;
	cbe->priv->pba.deg_ncdm            = NULL;
	cbe->priv->pba.ncdm_psd_parameters = NULL;
	cbe->priv->pba.ncdm_psd_files      = NULL;

  {
    const guint N_ncdm = nc_hicosmo_NMassNu (cosmo);

    if (N_ncdm != 0)
    {
      struct background* pba = &cbe->priv->pba;
      struct precision* ppr  = (struct precision*)cbe->prec->priv;
      const gdouble T_gamma0 = nc_hicosmo_T_gamma0 (cosmo);
      guint nu_i;

      cbe->priv->pba.N_ncdm = N_ncdm;

      pba->T_ncdm                   = (gdouble *)malloc (sizeof (gdouble) * N_ncdm);
      pba->ksi_ncdm                 = (gdouble *)malloc (sizeof (gdouble) * N_ncdm);
      pba->deg_ncdm                 = (gdouble *)malloc (sizeof (gdouble) * N_ncdm);
      pba->Omega0_ncdm              = (gdouble *)malloc (sizeof (gdouble) * N_ncdm);
      pba->M_ncdm                   = (gdouble *)malloc (sizeof (gdouble) * N_ncdm);
      pba->m_ncdm_in_eV             = (gdouble *)malloc (sizeof (gdouble) * N_ncdm);
      pba->got_files                = (gboolean *)malloc (sizeof (gboolean) * N_ncdm);

      pba->ncdm_quadrature_strategy = (gint *)malloc (sizeof (gint) * N_ncdm);
      pba->ncdm_input_q_size        = (gint *)malloc (sizeof (gint) * N_ncdm);
      pba->ncdm_qmax                = (gdouble *)malloc (sizeof (gdouble) * N_ncdm);
        
      pba->Omega0_ncdm_tot = 0.0;
      for (nu_i = 0; nu_i < pba->N_ncdm; nu_i++)
      {
        pba->got_files[nu_i] = _FALSE_;

        nc_hicosmo_MassNuInfo (cosmo, nu_i, &pba->m_ncdm_in_eV[nu_i], &pba->T_ncdm[nu_i], &pba->ksi_ncdm[nu_i], &pba->deg_ncdm[nu_i]);
        pba->M_ncdm[nu_i]      = pba->m_ncdm_in_eV[nu_i] * ncm_c_eV () / (ncm_c_kb () * pba->T_ncdm[nu_i] * T_gamma0);
        pba->Omega0_ncdm[nu_i] = nc_hicosmo_Omega_mnu0_n (cosmo, nu_i);

        pba->Omega0_ncdm_tot  += pba->Omega0_ncdm[nu_i];

        pba->ncdm_quadrature_strategy[nu_i] = qm_auto;
        pba->ncdm_input_q_size[nu_i]        = -1;
        pba->ncdm_qmax[nu_i]                = 15.0;
      }

      /* From CLASS input.c */
      /*---------------------------------------------------------*/
      background_ncdm_init (ppr, &cbe->priv->pba);

      /* We must calculate M from omega or vice versa if one of them is missing.
         If both are present, we must update the degeneracy parameter to
         reflect the implicit normalisation of the distribution function.*/
      
      guint n;
      gdouble rho_ncdm;
      for (n = 0; n < N_ncdm; n++)
      {
        if (pba->m_ncdm_in_eV[n] != 0.0)
        {
          /* Case of only mass or mass and Omega/omega: */
          pba->M_ncdm[n] = pba->m_ncdm_in_eV[n] / _k_B_ * _eV_ / pba->T_ncdm[n] / pba->T_cmb;
          background_ncdm_momenta (pba->q_ncdm_bg[n],
                                   pba->w_ncdm_bg[n],
                                   pba->q_size_ncdm_bg[n],
                                   pba->M_ncdm[n],
                                   pba->factor_ncdm[n],
                                   0.,
                                   NULL,
                                   &rho_ncdm,
                                   NULL,
                                   NULL,
                                   NULL);

          pba->Omega0_ncdm[n] = rho_ncdm / pba->H0 / pba->H0;
        }
        else
        {
          /* Case of only Omega/omega: */
          background_ncdm_M_from_Omega (ppr, pba, n);
          //printf("M_ncdm:%g\n",pba->M_ncdm[n]);
          pba->m_ncdm_in_eV[n] = _k_B_ / _eV_ * pba->T_ncdm[n] * pba->M_ncdm[n] * pba->T_cmb;
        }
        pba->Omega0_ncdm_tot += pba->Omega0_ncdm[n];
        //printf("Adding %g to total Omega..\n",pba->Omega0_ncdm[n]);
      }
      /*---------------------------------------------------------*/
    }
    else
    {
      cbe->priv->pba.N_ncdm         = 0;
      cbe->priv->pba.T_ncdm_default = 0.71611;
      cbe->priv->pba.T_ncdm         = NULL;
      cbe->priv->pba.m_ncdm_in_eV   = NULL;
    }
  }

  cbe->priv->pba.Omega0_scf          = 0.0;
	cbe->priv->pba.attractor_ic_scf    = _TRUE_;
	cbe->priv->pba.scf_parameters      = NULL;

  cbe->priv->pba.scf_parameters_size = 0;
	cbe->priv->pba.scf_tuning_index    = 0;
	cbe->priv->pba.phi_ini_scf         = 1;
	cbe->priv->pba.phi_prime_ini_scf   = 1;

	cbe->priv->pba.Omega0_k = nc_hicosmo_Omega_k0 (cosmo);
	if (fabs (cbe->priv->pba.Omega0_k) > NC_HICOSMO_OMEGA_K0_LIMIT)
	{
		cbe->priv->pba.a_today = 1.0;
		cbe->priv->pba.K       = -cbe->priv->pba.Omega0_k * gsl_pow_2 (cbe->priv->pba.a_today * cbe->priv->pba.H0);
		cbe->priv->pba.sgnK    = GSL_SIGN (cbe->priv->pba.K);
  }
	else
	{
		cbe->priv->pba.Omega0_k = 0.0;
		cbe->priv->pba.a_today  = 1.0;
		cbe->priv->pba.K        = 0.0;
		cbe->priv->pba.sgnK     = 0;
	}

  if (NC_IS_HICOSMO_DE_XCDM (cosmo))
  {
    const gdouble w0       = ncm_model_orig_param_get (NCM_MODEL (cosmo), NC_HICOSMO_DE_XCDM_W);
    const gdouble Omega_X0 = ncm_model_orig_param_get (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X);

    if (w0 != -1.0)
    {
      cbe->priv->pba.fluid_equation_of_state = CLP;
      cbe->priv->pba.Omega0_lambda           = 0.0;      /* Disable cosmological constant in CLASS */

      cbe->priv->pba.Omega0_fld              = Omega_X0; /* Enable DE fluid in CLASS */
      cbe->priv->pba.w0_fld                  = w0;
      cbe->priv->pba.wa_fld                  = 0.0;
      cbe->priv->pba.cs2_fld                 = 1.0;      
    }
    else
    {
      cbe->priv->pba.Omega0_lambda           = Omega_X0; /* Enable cosmological constant in CLASS */

      cbe->priv->pba.Omega0_fld              = 0.0;      /* Disable DE fluid in CLASS */
      cbe->priv->pba.w0_fld                  = -1.0;
      cbe->priv->pba.wa_fld                  = 0.0;
      cbe->priv->pba.cs2_fld                 = 1.0;
    }
  }
  else if (NC_IS_HICOSMO_DE_CPL (cosmo))
  {
    const gdouble w0       = ncm_model_orig_param_get (NCM_MODEL (cosmo), NC_HICOSMO_DE_CPL_W0);
    const gdouble w1       = ncm_model_orig_param_get (NCM_MODEL (cosmo), NC_HICOSMO_DE_CPL_W1);
    const gdouble Omega_X0 = ncm_model_orig_param_get (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X);

    cbe->priv->pba.fluid_equation_of_state = CLP;
    cbe->priv->pba.Omega0_lambda           = 0.0;      /* Disable cosmological constant in CLASS */

    cbe->priv->pba.Omega0_fld              = Omega_X0; /* Enable DE fluid in CLASS */
    cbe->priv->pba.w0_fld                  = w0;
    cbe->priv->pba.wa_fld                  = w1;
    cbe->priv->pba.cs2_fld                 = 1.0;
  }
  else
    g_error ("_nc_cbe_set_bg: CLASS in not compatible with the model `%s'.", G_OBJECT_TYPE_NAME (cosmo));

  cbe->priv->pba.Omega_EDE          = 0.0;
  cbe->priv->pba.c_gamma_over_c_fld = 0.4;
  cbe->priv->pba.Omega_ini_dcdm     = 0.0;
  
  cbe->priv->pba.shooting_failed    = _FALSE_;
	cbe->priv->pba.background_verbose = cbe->bg_verbose;
}

static void
_nc_cbe_set_thermo (NcCBE* cbe, NcHICosmo* cosmo)
{
	NcHIReion* reion = nc_hicosmo_peek_reion (cosmo);
	struct precision* ppr = (struct precision*)cbe->prec->priv;

	cbe->priv->pth.YHe                  = nc_hicosmo_Yp_4He (cosmo);
	cbe->priv->pth.recombination        = recfast;
	cbe->priv->pth.reio_parametrization = (reion == NULL) ? reio_none : reio_camb;

  g_assert (NC_IS_HIREION_CAMB (reion));

  cbe->priv->pth.reio_z_or_tau = reio_z;
  cbe->priv->pth.z_reio        = ncm_model_orig_param_get (NCM_MODEL (reion), NC_HIREION_CAMB_HII_HEII_Z);
  cbe->priv->pth.tau_reio      = nc_hireion_get_tau (reion, cosmo);

	cbe->priv->pth.reionization_exponent      = 1.5;
	cbe->priv->pth.reionization_width         = 0.5;
	cbe->priv->pth.helium_fullreio_redshift   = 3.5;
	cbe->priv->pth.helium_fullreio_width      = 0.5;

  cbe->priv->pth.binned_reio_num            = 0;
	cbe->priv->pth.binned_reio_z              = NULL;
	cbe->priv->pth.binned_reio_xe             = NULL;
	cbe->priv->pth.binned_reio_step_sharpness = 0.3;
  
	cbe->priv->pth.annihilation               = 0.0;
	cbe->priv->pth.decay                      = 0.0;
	cbe->priv->pth.annihilation_variation     = 0.0;
	cbe->priv->pth.annihilation_z             = 1000.0;
	cbe->priv->pth.annihilation_zmax          = 2500.0;
	cbe->priv->pth.annihilation_zmin          = 30.0;
	cbe->priv->pth.annihilation_f_halo        = 0.0;
	cbe->priv->pth.annihilation_z_halo        = 30.0;
	cbe->priv->pth.has_on_the_spot            = _TRUE_;

	cbe->priv->pth.compute_cb2_derivatives    = _FALSE_;
  cbe->priv->pth.compute_damping_scale      = _FALSE_;

	cbe->priv->pth.thermodynamics_verbose     = cbe->thermo_verbose;

	if ((ppr->tight_coupling_approximation == (gint)first_order_CLASS) ||
	    (ppr->tight_coupling_approximation == (gint)second_order_CLASS))
		cbe->priv->pth.compute_cb2_derivatives = _TRUE_;
}

static void
_nc_cbe_set_pert (NcCBE* cbe, NcHICosmo* cosmo)
{
	gboolean has_cls = (cbe->target_Cls & NC_DATA_CMB_TYPE_ALL) != 0;
	gboolean has_perturbations = has_cls || cbe->calc_transfer;
	struct precision* ppr = (struct precision*)cbe->prec->priv;

	/*
   * Inside CLASS they compare booleans with _TRUE_ and _FALSE_.
   * This is a bad idea, but to be compatible we must always
   * use their _TRUE_ and _FALSE_.
   */

	if (cbe->target_Cls & (NC_DATA_CMB_TYPE_TB | NC_DATA_CMB_TYPE_EB))
		g_error ("_nc_cbe_set_pert: modes TB and EB are not supported.");

	cbe->priv->ppt.has_perturbations                   = has_perturbations ? _TRUE_ : _FALSE_;
	cbe->priv->ppt.has_cls                             = has_cls ? _TRUE_ : _FALSE_;

	cbe->priv->ppt.has_cl_cmb_temperature              = cbe->target_Cls & NC_DATA_CMB_TYPE_TT ? _TRUE_ : _FALSE_;
	cbe->priv->ppt.has_cl_cmb_polarization             = cbe->target_Cls & (NC_DATA_CMB_TYPE_EE | NC_DATA_CMB_TYPE_BB | NC_DATA_CMB_TYPE_TE) ? _TRUE_ : _FALSE_;
	cbe->priv->ppt.has_cl_cmb_lensing_potential        = ((cbe->target_Cls & NC_DATA_CMB_TYPE_PHIPHI) || (cbe->use_lensed_Cls)) ? _TRUE_ : _FALSE_;
	cbe->priv->ppt.has_cl_number_count                 = _FALSE_;
	cbe->priv->ppt.has_cl_lensing_potential            = _FALSE_;
	cbe->priv->ppt.has_pk_matter                       = cbe->calc_transfer ? _TRUE_ : _FALSE_;
	cbe->priv->ppt.has_density_transfers               = _FALSE_;
	cbe->priv->ppt.has_velocity_transfers              = _FALSE_;
  cbe->priv->ppt.has_metricpotential_transfers       = _FALSE_;

	cbe->priv->ppt.has_nl_corrections_based_on_delta_m = _FALSE_;

	cbe->priv->ppt.has_nc_density                      = _FALSE_;
	cbe->priv->ppt.has_nc_rsd                          = _FALSE_;
	cbe->priv->ppt.has_nc_lens                         = _FALSE_;
	cbe->priv->ppt.has_nc_gr                           = _FALSE_;

	cbe->priv->ppt.switch_sw                           = 1;
	cbe->priv->ppt.switch_eisw                         = 1;
	cbe->priv->ppt.switch_lisw                         = 1;
	cbe->priv->ppt.switch_dop                          = 1;
	cbe->priv->ppt.switch_pol                          = 1;
	cbe->priv->ppt.eisw_lisw_split_z                   = 120.0;

	cbe->priv->ppt.has_ad                              = _TRUE_;
	cbe->priv->ppt.has_bi                              = _FALSE_;
	cbe->priv->ppt.has_cdi                             = _FALSE_;
	cbe->priv->ppt.has_nid                             = _FALSE_;
	cbe->priv->ppt.has_niv                             = _FALSE_;

	cbe->priv->ppt.has_perturbed_recombination         = _FALSE_;
	cbe->priv->ppt.tensor_method                       = tm_massless_approximation;
	cbe->priv->ppt.evolve_tensor_ur                    = _FALSE_;
	cbe->priv->ppt.evolve_tensor_ncdm                  = _FALSE_;

	cbe->priv->ppt.has_scalars                         = _TRUE_;
	cbe->priv->ppt.has_vectors                         = _FALSE_;
	cbe->priv->ppt.has_tensors                         = cbe->use_tensor ? _TRUE_ : _FALSE_;

	cbe->priv->ppt.l_scalar_max                        = cbe->scalar_lmax + (cbe->use_lensed_Cls ? ppr->delta_l_max : 0);

	cbe->priv->ppt.l_vector_max                        = cbe->vector_lmax;
	cbe->priv->ppt.l_tensor_max                        = cbe->tensor_lmax;
	cbe->priv->ppt.l_lss_max                           = 300;

	cbe->priv->ppt.gauge                               = synchronous;

	cbe->priv->ppt.k_output_values_num                 = 0;
  /*cbe->priv->ppt.k_output_values[0]                  = NULL;*/
	cbe->priv->ppt.store_perturbations                 = _FALSE_;
	cbe->priv->ppt.number_of_scalar_titles             = 0;
	cbe->priv->ppt.number_of_vector_titles             = 0;
	cbe->priv->ppt.number_of_tensor_titles             = 0;
	{
		guint filenum;
		for (filenum = 0; filenum < _MAX_NUMBER_OF_K_FILES_; filenum++)
		{
			cbe->priv->ppt.scalar_perturbations_data[filenum] = NULL;
			cbe->priv->ppt.vector_perturbations_data[filenum] = NULL;
			cbe->priv->ppt.tensor_perturbations_data[filenum] = NULL;
		}
	}
	cbe->priv->ppt.index_k_output_values               = NULL;

  cbe->priv->ppt.three_ceff2_ur                      = 1.0;
  cbe->priv->ppt.three_cvis2_ur                      = 1.0;

  /* This is set elsewhere */
  /*cbe->priv->ppt.z_max_pk                            = 0.0;*/
  
	cbe->priv->ppt.selection_num                       = 1;
	cbe->priv->ppt.selection                           = gaussian;
	cbe->priv->ppt.selection_mean[0]                   = 1.0;
	cbe->priv->ppt.selection_width[0]                  = 0.1;

	cbe->priv->ppt.perturbations_verbose = cbe->pert_verbose;
}

static gdouble
_external_Pk_callback_pks (const double lnk, gpointer data)
{
	NcHIPrim* prim = NC_HIPRIM (data);

	return nc_hiprim_lnSA_powspec_lnk (prim, lnk);
}

static gdouble
_external_Pk_callback_pkt (const double lnk, gpointer data)
{
	NcHIPrim* prim = NC_HIPRIM (data);

	return nc_hiprim_lnT_powspec_lnk (prim, lnk);
}

static void
_nc_cbe_set_prim (NcCBE* cbe, NcHICosmo* cosmo)
{
	NcHIPrim* prim = nc_hicosmo_peek_prim (cosmo);

	cbe->priv->ppm.primordial_spec_type = external_Pk_callback;
	/*cbe->priv->ppm.primordial_spec_type = analytic_Pk;*/
	cbe->priv->ppm.external_Pk_callback_pks = &_external_Pk_callback_pks;

  if (cbe->use_tensor)
	{
		g_assert (ncm_model_check_impl_opt (NCM_MODEL (prim), NC_HIPRIM_IMPL_lnT_powspec_lnk));
		cbe->priv->ppm.external_Pk_callback_pkt = &_external_Pk_callback_pkt;
	}
	cbe->priv->ppm.external_Pk_callback_data = prim;

	cbe->priv->ppm.k_pivot          = 0.05;
	cbe->priv->ppm.A_s              = 2.40227188179e-9;
	cbe->priv->ppm.n_s              = 0.9742;
	cbe->priv->ppm.alpha_s          = 0.0;
	cbe->priv->ppm.f_bi             = 1.0;
	cbe->priv->ppm.n_bi             = 1.0;
	cbe->priv->ppm.alpha_bi         = 0.0;
	cbe->priv->ppm.f_cdi            = 1.0;
	cbe->priv->ppm.n_cdi            = 1.0;
	cbe->priv->ppm.alpha_cdi        = 0.0;
	cbe->priv->ppm.f_nid            = 1.0;
	cbe->priv->ppm.n_nid            = 1.0;
	cbe->priv->ppm.alpha_nid        = 0.0;
	cbe->priv->ppm.f_niv            = 1.0;
	cbe->priv->ppm.n_niv            = 1.0;
	cbe->priv->ppm.alpha_niv        = 0.0;
	cbe->priv->ppm.c_ad_bi          = 0.0;
	cbe->priv->ppm.n_ad_bi          = 0.0;
	cbe->priv->ppm.alpha_ad_bi      = 0.0;
	cbe->priv->ppm.c_ad_cdi         = 0.0;
	cbe->priv->ppm.n_ad_cdi         = 0.0;
	cbe->priv->ppm.alpha_ad_cdi     = 0.0;
	cbe->priv->ppm.c_ad_nid         = 0.0;
	cbe->priv->ppm.n_ad_nid         = 0.0;
	cbe->priv->ppm.alpha_ad_nid     = 0.0;
	cbe->priv->ppm.c_ad_niv         = 0.0;
	cbe->priv->ppm.n_ad_niv         = 0.0;
	cbe->priv->ppm.alpha_ad_niv     = 0.0;
	cbe->priv->ppm.c_bi_cdi         = 0.0;
	cbe->priv->ppm.n_bi_cdi         = 0.0;
	cbe->priv->ppm.alpha_bi_cdi     = 0.0;
	cbe->priv->ppm.c_bi_nid         = 0.0;
	cbe->priv->ppm.n_bi_nid         = 0.0;
	cbe->priv->ppm.alpha_bi_nid     = 0.0;
	cbe->priv->ppm.c_bi_niv         = 0.0;
	cbe->priv->ppm.n_bi_niv         = 0.0;
	cbe->priv->ppm.alpha_bi_niv     = 0.0;
	cbe->priv->ppm.c_cdi_nid        = 0.0;
	cbe->priv->ppm.n_cdi_nid        = 0.0;
	cbe->priv->ppm.alpha_cdi_nid    = 0.0;
	cbe->priv->ppm.c_cdi_niv        = 0.0;
	cbe->priv->ppm.n_cdi_niv        = 0.0;
	cbe->priv->ppm.alpha_cdi_niv    = 0.0;
	cbe->priv->ppm.c_nid_niv        = 0.0;
	cbe->priv->ppm.n_nid_niv        = 0.0;
	cbe->priv->ppm.alpha_nid_niv    = 0.0;
	cbe->priv->ppm.r                = 1.0;
	cbe->priv->ppm.n_t              = -cbe->priv->ppm.r / 8.0 * (2.0 - cbe->priv->ppm.r / 8.0 - cbe->priv->ppm.n_s);
	cbe->priv->ppm.alpha_t          = cbe->priv->ppm.r / 8.0 * (cbe->priv->ppm.r / 8.0 + cbe->priv->ppm.n_s - 1.0);
	cbe->priv->ppm.potential        = polynomial;
	cbe->priv->ppm.phi_end          = 0.0;
  cbe->priv->ppm.phi_pivot_method = N_star;
  cbe->priv->ppm.phi_pivot_target = 60;
	cbe->priv->ppm.V0               = 1.25e-13;
	cbe->priv->ppm.V1               = -1.12e-14;
	cbe->priv->ppm.V2               = -6.95e-14;
	cbe->priv->ppm.V3               = 0.0;
	cbe->priv->ppm.V4               = 0.0;
	cbe->priv->ppm.H0               = 3.69e-6;
	cbe->priv->ppm.H1               = -5.84e-7;
	cbe->priv->ppm.H2               = 0.0;
	cbe->priv->ppm.H3               = 0.0;
	cbe->priv->ppm.H4               = 0.0;
  cbe->priv->ppm.behavior         = numerical;
	cbe->priv->ppm.command          = "write here your command for the external Pk";
	cbe->priv->ppm.custom1          = 0.0;
	cbe->priv->ppm.custom2          = 0.0;
	cbe->priv->ppm.custom3          = 0.0;
	cbe->priv->ppm.custom4          = 0.0;
	cbe->priv->ppm.custom5          = 0.0;
	cbe->priv->ppm.custom6          = 0.0;
	cbe->priv->ppm.custom7          = 0.0;
	cbe->priv->ppm.custom8          = 0.0;
	cbe->priv->ppm.custom9          = 0.0;
	cbe->priv->ppm.custom10         = 0.0;

	cbe->priv->ppm.primordial_verbose = cbe->prim_verbose;
}

static void
_nc_cbe_set_transfer (NcCBE* cbe, NcHICosmo* cosmo)
{

  cbe->priv->ptr.selection_bias[0]               = 1.0;
  cbe->priv->ptr.selection_magnification_bias[0] = 0.0;
  
	cbe->priv->ptr.lcmb_rescale                    = 1.0;
	cbe->priv->ptr.lcmb_pivot                      = 0.1;
	cbe->priv->ptr.lcmb_tilt                       = 0.0;
	cbe->priv->ptr.initialise_HIS_cache            = _FALSE_;
	cbe->priv->ptr.has_nz_analytic                 = _FALSE_;
	cbe->priv->ptr.has_nz_file                     = _FALSE_;
	cbe->priv->ptr.has_nz_evo_analytic             = _FALSE_;
	cbe->priv->ptr.has_nz_evo_file                 = _FALSE_;

	cbe->priv->ptr.transfer_verbose                = cbe->transfer_verbose;
}

static void
_nc_cbe_set_spectra (NcCBE* cbe, NcHICosmo* cosmo)
{
	/*cbe->priv->psp.z_max_pk = 0.0;*/
	cbe->priv->psp.non_diag = 0;

	cbe->priv->psp.spectra_verbose = cbe->spectra_verbose;
}

static void
_nc_cbe_set_lensing (NcCBE* cbe, NcHICosmo* cosmo)
{
	cbe->priv->ple.has_lensed_cls  = cbe->use_lensed_Cls ? _TRUE_ : _FALSE_;
	cbe->priv->ple.lensing_verbose = cbe->lensing_verbose;
}

static void
_nc_cbe_set_nonlin (NcCBE* cbe, NcHICosmo* cosmo)
{
	cbe->priv->pnl.method    = nl_none;
  cbe->priv->pnl.has_pk_eq = _FALSE_;
	cbe->priv->pnl.nonlinear_verbose = cbe->nonlin_verbose;
}

static void _nc_cbe_free_bg (NcCBE* cbe);
static void _nc_cbe_free_thermo (NcCBE* cbe);
static void _nc_cbe_free_pert (NcCBE* cbe);
static void _nc_cbe_free_prim (NcCBE* cbe);
static void _nc_cbe_free_nonlin (NcCBE* cbe);
static void _nc_cbe_free_transfer (NcCBE* cbe);
static void _nc_cbe_free_spectra (NcCBE* cbe);
static void _nc_cbe_free_lensing (NcCBE* cbe);

static void
_nc_cbe_call_bg (NcCBE* cbe, NcHICosmo* cosmo)
{
	struct precision* ppr = (struct precision*)cbe->prec->priv;
	_nc_cbe_set_bg (cbe, cosmo);
  if (background_init (ppr, &cbe->priv->pba) == _FAILURE_)
		g_error ("_nc_cbe_call_bg: Error running background_init `%s'\n", cbe->priv->pba.error_message);
}

static void
_nc_cbe_call_thermo (NcCBE* cbe, NcHICosmo* cosmo)
{
	struct precision* ppr = (struct precision*)cbe->prec->priv;

	_nc_cbe_call_bg (cbe, cosmo);

	_nc_cbe_set_thermo (cbe, cosmo);

  if (FALSE)
  {
    struct background *pba = &cbe->priv->pba;
    struct thermo *pth     = &cbe->priv->pth;

    printf ("#NC:  % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g %d % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g %d % 22.15g % 22.15g\n", 
            pba->h, pba->T_cmb, pba->Omega0_g, pba->Omega0_ur, 
            pba->Omega0_b, pba->Omega0_cdm, pba->Omega0_dcdmdr, pba->Omega0_dcdm, 
            pba->Gamma_dcdm, pba->N_ncdm, pba->Omega0_ncdm_tot, pba->Omega0_scf, 
            pba->Omega0_k, pba->Omega0_lambda, pba->Omega0_fld, pba->a_today, 
            pba->w0_fld, pba->wa_fld, pba->cs2_fld, pba->use_ppf, pba->c_gamma_over_c_fld, 
            pth->z_reio);
  }

	if (thermodynamics_init (ppr, &cbe->priv->pba, &cbe->priv->pth) == _FAILURE_)
		g_error ("_nc_cbe_call_thermo: Error running thermodynamics_init `%s'\n", cbe->priv->pth.error_message);
}

static void
_nc_cbe_call_pert (NcCBE* cbe, NcHICosmo* cosmo)
{
	struct precision* ppr = (struct precision*)cbe->prec->priv;

	/*_nc_cbe_call_thermo (cbe, cosmo);*/
	cbe->free = &_nc_cbe_free_pert;

	_nc_cbe_set_pert (cbe, cosmo);
	if (perturb_init (ppr, &cbe->priv->pba, &cbe->priv->pth, &cbe->priv->ppt) == _FAILURE_)
		g_error ("_nc_cbe_call_pert: Error running perturb_init `%s'\n", cbe->priv->ppt.error_message);
}

static void
_nc_cbe_call_prim (NcCBE* cbe, NcHICosmo* cosmo)
{
	struct precision* ppr = (struct precision*)cbe->prec->priv;

	_nc_cbe_call_pert (cbe, cosmo);
	cbe->free = &_nc_cbe_free_prim;

	_nc_cbe_set_prim (cbe, cosmo);
	if (primordial_init (ppr, &cbe->priv->ppt, &cbe->priv->ppm) == _FAILURE_)
		g_error ("_nc_cbe_call_prim: Error running primordial_init `%s'\n", cbe->priv->ppm.error_message);
}

static void
_nc_cbe_call_nonlin (NcCBE* cbe, NcHICosmo* cosmo)
{
	struct precision* ppr = (struct precision*)cbe->prec->priv;

	_nc_cbe_call_prim (cbe, cosmo);
	cbe->free = &_nc_cbe_free_nonlin;

	_nc_cbe_set_nonlin (cbe, cosmo);
	if (nonlinear_init (ppr, &cbe->priv->pba, &cbe->priv->pth, &cbe->priv->ppt, &cbe->priv->ppm, &cbe->priv->pnl) == _FAILURE_)
		g_error ("_nc_cbe_call_nonlin: Error running nonlinear_init `%s'\n", cbe->priv->pnl.error_message);
}

static void
_nc_cbe_call_transfer (NcCBE* cbe, NcHICosmo* cosmo)
{
	struct precision* ppr = (struct precision*)cbe->prec->priv;

	_nc_cbe_call_nonlin (cbe, cosmo);
	cbe->free = &_nc_cbe_free_transfer;

	_nc_cbe_set_transfer (cbe, cosmo);
	if (transfer_init (ppr, &cbe->priv->pba, &cbe->priv->pth, &cbe->priv->ppt, &cbe->priv->pnl, &cbe->priv->ptr) == _FAILURE_)
		g_error ("_nc_cbe_call_transfer: Error running transfer_init `%s'\n", cbe->priv->ptr.error_message);
}

static void
_nc_cbe_call_spectra (NcCBE* cbe, NcHICosmo* cosmo)
{
	struct precision* ppr = (struct precision*)cbe->prec->priv;

	_nc_cbe_call_transfer (cbe, cosmo);
	cbe->free = &_nc_cbe_free_spectra;

	_nc_cbe_set_spectra (cbe, cosmo);
	if (spectra_init (ppr, &cbe->priv->pba, &cbe->priv->ppt, &cbe->priv->ppm, &cbe->priv->pnl, &cbe->priv->ptr, &cbe->priv->psp) == _FAILURE_)
		g_error ("_nc_cbe_call_spectra: Error running spectra_init `%s'\n", cbe->priv->psp.error_message);
}

static void
_nc_cbe_call_lensing (NcCBE* cbe, NcHICosmo* cosmo)
{
	struct precision* ppr = (struct precision*)cbe->prec->priv;

	_nc_cbe_call_spectra (cbe, cosmo);
	cbe->free = &_nc_cbe_free_lensing;

	_nc_cbe_set_lensing (cbe, cosmo);
	if (lensing_init (ppr, &cbe->priv->ppt, &cbe->priv->psp, &cbe->priv->pnl, &cbe->priv->ple) == _FAILURE_)
		g_error ("_nc_cbe_call_lensing: Error running lensing_init `%s'\n", cbe->priv->ple.error_message);
}

static void
_nc_cbe_free_bg (NcCBE* cbe)
{
	if (background_free (&cbe->priv->pba) == _FAILURE_)
		g_error ("_nc_cbe_free_bg: Error running background_free `%s'\n", cbe->priv->pba.error_message);
}

static void
_nc_cbe_free_thermo (NcCBE* cbe)
{
	if (thermodynamics_free (&cbe->priv->pth) == _FAILURE_)
		g_error ("_nc_cbe_free_thermo: Error running thermodynamics_free `%s'\n", cbe->priv->pth.error_message);

	_nc_cbe_free_bg (cbe);
}

static void
_nc_cbe_free_pert (NcCBE* cbe)
{
	if (perturb_free (&cbe->priv->ppt) == _FAILURE_)
		g_error ("_nc_cbe_free_pert: Error running perturb_free `%s'\n", cbe->priv->ppt.error_message);

	/*_nc_cbe_free_thermo (cbe);*/
}

static void
_nc_cbe_free_prim (NcCBE* cbe)
{
	if (primordial_free (&cbe->priv->ppm) == _FAILURE_)
		g_error ("_nc_cbe_free_prim: Error running primordial_free `%s'\n", cbe->priv->ppm.error_message);

	_nc_cbe_free_pert (cbe);
}

static void
_nc_cbe_free_nonlin (NcCBE* cbe)
{
	if (nonlinear_free (&cbe->priv->pnl) == _FAILURE_)
		g_error ("_nc_cbe_free_nonlin: Error running nonlinear_free `%s'\n", cbe->priv->pnl.error_message);

	_nc_cbe_free_prim (cbe);
}

static void
_nc_cbe_free_transfer (NcCBE* cbe)
{
	if (transfer_free (&cbe->priv->ptr) == _FAILURE_)
		g_error ("_nc_cbe_free_transfer: Error running transfer_free `%s'\n", cbe->priv->ptr.error_message);

	_nc_cbe_free_nonlin (cbe);
}

static void
_nc_cbe_free_spectra (NcCBE* cbe)
{
	if (spectra_free (&cbe->priv->psp) == _FAILURE_)
		g_error ("_nc_cbe_free_lensing: Error running spectra_free `%s'\n", cbe->priv->psp.error_message);

	_nc_cbe_free_transfer (cbe);
}

static void
_nc_cbe_free_lensing (NcCBE* cbe)
{
	if (lensing_free (&cbe->priv->ple) == _FAILURE_)
		g_error ("_nc_cbe_free_lensing: Error running lensing_free `%s'\n", cbe->priv->ple.error_message);

	_nc_cbe_free_spectra (cbe);
}

static void
_nc_cbe_update_callbacks (NcCBE* cbe)
{
	gboolean has_Cls = cbe->target_Cls & NC_DATA_CMB_TYPE_ALL;

	ncm_model_ctrl_force_update (cbe->ctrl_cosmo);

	if (cbe->allocated)
	{
		g_assert (cbe->free != NULL);
		cbe->free (cbe);
		cbe->call = NULL;
		cbe->free = NULL;
		cbe->allocated = FALSE;
	}

	if (has_Cls && cbe->use_lensed_Cls)
		cbe->call = _nc_cbe_call_lensing;
	else if (has_Cls || cbe->calc_transfer)
		cbe->call = _nc_cbe_call_spectra;
}

/**
 * nc_cbe_thermodyn_prepare:
 * @cbe: a #NcCBE
 * @cosmo: a #NcHICosmo
 *
 * Prepares the thermodynamic Class structure.
 *
 */
void nc_cbe_thermodyn_prepare (NcCBE* cbe, NcHICosmo* cosmo)
{
	if (cbe->thermodyn_prepared)
	{
		_nc_cbe_free_thermo (cbe);
		cbe->thermodyn_prepared = FALSE;
	}

	_nc_cbe_call_thermo (cbe, cosmo);
	cbe->thermodyn_prepared = TRUE;
}

/**
 * nc_cbe_thermodyn_prepare_if_needed:
 * @cbe: a #NcCBE
 * @cosmo: a #NcHICosmo
 *
 * Prepares the thermodynamic Class structure.
 *
 */
void nc_cbe_thermodyn_prepare_if_needed (NcCBE* cbe, NcHICosmo* cosmo)
{
	if (ncm_model_ctrl_update (cbe->ctrl_cosmo, NCM_MODEL (cosmo)))
	{
		nc_cbe_thermodyn_prepare (cbe, cosmo);
		ncm_model_ctrl_force_update (cbe->ctrl_prim);
	}
}

/**
 * nc_cbe_prepare:
 * @cbe: a #NcCBE
 * @cosmo: a #NcHICosmo
 *
 * Prepares all necessary Class structures.
 *
 */
void nc_cbe_prepare (NcCBE* cbe, NcHICosmo* cosmo)
{
	/*printf ("Preparing CLASS!\n");*/
	if (ncm_model_peek_submodel_by_mid (NCM_MODEL (cosmo), nc_hiprim_id ()) == NULL)
	{
		g_error ("nc_cbe_prepare: cosmo model must contain a NcHIPrim submodel.");
	}

  if (cbe->allocated)
	{
		g_assert (cbe->free != NULL);
		cbe->free (cbe);
		cbe->allocated = FALSE;
	}

	if (cbe->thermodyn_prepared)
	{
		_nc_cbe_free_thermo (cbe);
		cbe->thermodyn_prepared = FALSE;
	}

	_nc_cbe_call_thermo (cbe, cosmo);
	cbe->thermodyn_prepared = TRUE;

	if (cbe->call != NULL)
	{
		cbe->call (cbe, cosmo);
		cbe->allocated = TRUE;
	}
}

/**
 * nc_cbe_prepare_if_needed:
 * @cbe: a #NcCBE
 * @cosmo: a #NcHICosmo
 *
 * Prepares all necessary Class structures.
 *
 */
void nc_cbe_prepare_if_needed (NcCBE* cbe, NcHICosmo* cosmo)
{
	ncm_model_ctrl_update (cbe->ctrl_cosmo, NCM_MODEL (cosmo));

	if (!ncm_model_ctrl_model_has_submodel (cbe->ctrl_cosmo, nc_hiprim_id ()))
	{
		g_error ("nc_cbe_prepare_if_needed: cosmo model must contain a NcHIPrim submodel.");
	}
	else
	{
		gboolean cosmo_up = ncm_model_ctrl_model_last_update (cbe->ctrl_cosmo);
		gboolean prim_up = ncm_model_ctrl_submodel_last_update (cbe->ctrl_cosmo, nc_hiprim_id ());

		/*printf ("cosmo_up %d prim_up %d [%p]\n", cosmo_up, prim_up, cosmo);*/

		if (cosmo_up)
		{
			nc_cbe_prepare (cbe, cosmo);
		}
		else if (prim_up)
		{
			if (cbe->allocated)
			{
				g_assert (cbe->free != NULL);
				cbe->free (cbe);
				cbe->allocated = FALSE;
			}
			if (cbe->call != NULL)
			{
				cbe->call (cbe, cosmo);
				cbe->allocated = TRUE;
			}
		}
	}
}

/**
 * nc_cbe_compare_bg:
 * @cbe: a #NcCBE
 * @cosmo: a #NcHICosmo
 * @log_cmp: whether to print the comparison
 * 
 * Compares CLASS and NumCosmo background calculations and returns the worst discrepancy.
 * 
 * Returns: worst error. 
 */
gdouble
nc_cbe_compare_bg (NcCBE *cbe, NcHICosmo *cosmo, gboolean log_cmp)
{
  nc_cbe_prepare (cbe, cosmo);
  {
    struct precision *ppr  = (struct precision *)cbe->prec->priv;
    struct background *pba = &cbe->priv->pba;
    const gdouble RH       = nc_hicosmo_RH_Mpc (cosmo);
    const gdouble zf       = 1.0 / ppr->a_ini_over_a_today_default;

    gdouble pvecback[pba->bg_size];
    gdouble err = 0.0;
    gboolean isLambda = NC_IS_HICOSMO_DE_XCDM (cosmo) && (ncm_model_orig_param_get (NCM_MODEL (cosmo), NC_HICOSMO_DE_XCDM_W) == -1.0);
    gboolean hasNcdm  = nc_hicosmo_Omega_mnu0 (cosmo) != 0.0;
    guint i;

    if (cbe->a == NULL)
      cbe->a = nc_scalefactor_new (zf, NULL);
    else
      nc_scalefactor_set_zf (cbe->a, zf);

    nc_scalefactor_prepare_if_needed (cbe->a, cosmo);

    for (i = 0; i < pba->bt_size; i++)
    {
      gint last_index = 0;
      background_at_tau (pba,
                         pba->tau_table[i],
                         pba->long_info,
                         pba->inter_normal,
                         &last_index,
                         pvecback);
      {
        const gdouble eta        = pba->tau_table[i] / RH;
        const gdouble a          = pvecback[pba->index_bg_a];
        const gdouble z          = pba->a_today / a - 1.0;    /* nc_scalefactor_z_eta (cbe->a, eta); */
        const gdouble E2         = nc_hicosmo_E2 (cosmo, z);
        const gdouble E          = sqrt (E2);
        const gdouble H          = E / RH;
        const gdouble RH_pow_m2  = 1.0 / (RH * RH);
        const gdouble H_prime    = -0.5 * nc_hicosmo_dE2_dz (cosmo, z) * RH_pow_m2;

        const gdouble rho_g      = nc_hicosmo_E2Omega_g  (cosmo, z) * RH_pow_m2;
        const gdouble rho_ur     = nc_hicosmo_E2Omega_nu (cosmo, z) * RH_pow_m2;
        const gdouble rho_ncdm1  = hasNcdm ? nc_hicosmo_E2Omega_mnu_n (cosmo, 0, z) * RH_pow_m2 : 0.0;
        const gdouble rho_b      = nc_hicosmo_E2Omega_b  (cosmo, z) * RH_pow_m2;
        const gdouble rho_cdm    = nc_hicosmo_E2Omega_c  (cosmo, z) * RH_pow_m2;
        const gdouble rho_Lambda = nc_hicosmo_de_E2Omega_de (NC_HICOSMO_DE (cosmo), z) * RH_pow_m2;
        const gdouble rho_k      = nc_hicosmo_E2Omega_k (cosmo, z) * RH_pow_m2;

        /* CLASS defines Omega_X = rho_X / rho_total */
        const gdouble E2Omega_t  = nc_hicosmo_E2Omega_t (cosmo, z);
        const gdouble Omega_r    = nc_hicosmo_E2Omega_r (cosmo, z) / E2Omega_t;
        const gdouble Omega_m    = nc_hicosmo_E2Omega_m (cosmo, z) / E2Omega_t;

        const gdouble rho_crit   = E2 * RH_pow_m2;

        {
          const gdouble a_diff          = fabs (nc_scalefactor_eval_a_eta (cbe->a, eta) / pvecback[pba->index_bg_a] - 1.0);
          const gdouble H_diff          = fabs (H / pvecback[pba->index_bg_H] - 1.0);
          const gdouble Hprime_diff     = fabs (H_prime / pvecback[pba->index_bg_H_prime] - 1.0);

          const gdouble rho_g_diff      = fabs (rho_g / pvecback[pba->index_bg_rho_g] - 1.0);
          const gdouble rho_ur_diff     = fabs (rho_ur / pvecback[pba->index_bg_rho_ur] - 1.0);
          const gdouble rho_ncdm_diff   = hasNcdm ? fabs (rho_ncdm1 / pvecback[pba->index_bg_rho_ncdm1] - 1.0) : 0.0;
          const gdouble rho_b_diff      = fabs (rho_b / pvecback[pba->index_bg_rho_b] - 1.0);
          const gdouble rho_cdm_diff    = fabs (rho_cdm / pvecback[pba->index_bg_rho_cdm] - 1.0);

          const gdouble rho_Lambda_diff = isLambda ? fabs (rho_Lambda / pvecback[pba->index_bg_rho_lambda] - 1.0) : fabs (rho_Lambda / pvecback[pba->index_bg_rho_fld] - 1.0);

          const gdouble Omega_m0_diff   = fabs (Omega_m / pvecback[pba->index_bg_Omega_m] - 1.0);
          const gdouble Omega_r0_diff   = fabs (Omega_r / pvecback[pba->index_bg_Omega_r] - 1.0);

          const gdouble Omega_k_diff    = (rho_k == -pba->K/a/a) ? 0.0 : fabs ((rho_k - (-pba->K/a/a)) / (-pba->K/a/a + 1.0-2));
          const gdouble rho_crit_diff   = fabs (rho_crit / pvecback[pba->index_bg_rho_crit] - 1.0);

          err = GSL_MAX (err, a_diff);
          err = GSL_MAX (err, H_diff);
          err = GSL_MAX (err, Hprime_diff);
          err = GSL_MAX (err, rho_g_diff);
          err = GSL_MAX (err, rho_ur_diff);
          err = GSL_MAX (err, rho_ncdm_diff);
          err = GSL_MAX (err, rho_b_diff);
          err = GSL_MAX (err, rho_cdm_diff);
          err = GSL_MAX (err, rho_Lambda_diff);
          err = GSL_MAX (err, Omega_r0_diff);
          err = GSL_MAX (err, Omega_m0_diff);
          err = GSL_MAX (err, Omega_k_diff);
          err = GSL_MAX (err, rho_crit_diff);

          if (log_cmp)
            printf ("# eta = % 22.15g a = % 10.5e | % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e % 10.5e [% 10.5e]\n", 
                    eta, pvecback[pba->index_bg_a], 
                    a_diff,          /* 01 */
                    H_diff,          /* 02 */
                    Hprime_diff,     /* 03 */
                    rho_g_diff,      /* 04 */
                    rho_ur_diff,     /* 05 */
                    rho_ncdm_diff,   /* 06 */
                    rho_b_diff,      /* 07 */
                    rho_cdm_diff,    /* 08 */
                    rho_Lambda_diff, /* 09 */
                    Omega_r0_diff,   /* 10 */
                    Omega_m0_diff,   /* 11 */
                    Omega_k_diff,    /* 12 */
                    rho_crit_diff,   /* 13 */
                    err);
        }
      }
    }
    if (log_cmp)
      printf ("# worst error: % 10.5e\n", err);

    return err;
  }
}


/**
 * nc_cbe_thermodyn_get_Xe:
 * @cbe: a #NcCBE
 *
 * Gets the free electrons fraction $X_e$ as a function of the redshift.
 *
 * Returns: (transfer full): a #NcmSpline for Xe.
 */
NcmSpline*
nc_cbe_thermodyn_get_Xe (NcCBE* cbe)
{
	const guint size = cbe->priv->pth.tt_size;
	NcmVector* z_v = ncm_vector_new (size);
	NcmVector* Xe_v = ncm_vector_new (size);
	NcmSpline* Xe_s = ncm_spline_cubic_notaknot_new_full (z_v, Xe_v, FALSE);
	guint i;

	for (i = 0; i < size; i++)
	{
		const gdouble z_i  = cbe->priv->pth.z_table[i];
		const gdouble Xe_i = cbe->priv->pth.thermodynamics_table[cbe->priv->pth.th_size * i + cbe->priv->pth.index_th_xe];

		ncm_vector_fast_set (z_v, size - 1 - i, -log (z_i + 1.0));
		ncm_vector_fast_set (Xe_v, size - 1 - i, Xe_i);
	}

	ncm_vector_clear (&z_v);
	ncm_vector_clear (&Xe_v);

	ncm_spline_prepare (Xe_s);

	return Xe_s;
}

/**
 * nc_cbe_thermodyn_v_tau_max_z:
 * @cbe: a #NcCBE
 *
 * Gets the redshift of the maximum visibility function.
 *
 * Returns: $z_\mathrm{rec}$.
 */
gdouble 
nc_cbe_thermodyn_v_tau_max_z (NcCBE *cbe)
{
  return cbe->priv->pth.z_rec;
}

/**
 * nc_cbe_thermodyn_z_d:
 * @cbe: a #NcCBE
 *
 * Gets drag redshift.
 *
 * Returns: $z_d$.
 */
gdouble 
nc_cbe_thermodyn_z_d (NcCBE *cbe)
{
  return cbe->priv->pth.z_d;
}

/**
 * nc_cbe_get_matter_ps:
 * @cbe: a #NcCBE
 *
 * Gets the logarithm base e of the matter power spectrum as a function of the redshift $z$ and mode $\ln (k)$.
 *
 * Returns: (transfer full): a #NcmSpline2d for the logarithm base e of the matter power spectrum, $\ln P(\ln k, z)$.
 */
NcmSpline2d *
nc_cbe_get_matter_ps (NcCBE* cbe)
{
	const gint z_size  = cbe->priv->psp.ln_tau_size;
	const gint partz   = 4;
	const gint z_tsize = (z_size - 1) * partz + 1;

	NcmVector* lnk_v  = ncm_vector_new (cbe->priv->psp.ln_k_size);
	NcmVector* z_v    = ncm_vector_new (z_tsize);
	NcmMatrix* lnPk   = ncm_matrix_new (z_tsize, cbe->priv->psp.ln_k_size);
	gdouble* pvecback = g_new (gdouble, cbe->priv->pba.bg_size_short);
  gdouble *temp     = g_new (gdouble, cbe->priv->psp.ln_k_size);
	gdouble z_i_a, z_i_b = 0.0;
	guint i, m;

	for (i = 0; i < cbe->priv->psp.ln_k_size; i++)
	{
		ncm_vector_set (lnk_v, i, cbe->priv->psp.ln_k[i]);
	}

	m = 0;
	i = 0;
	{
		gdouble ln_tau_i;
		gint last_index_back;
		guint j;

		z_i_a = 0.0;

		ln_tau_i = cbe->priv->psp.ln_tau[z_size - i - 2];
		background_at_tau (&cbe->priv->pba,
		                   exp (ln_tau_i),
		                   cbe->priv->pba.short_info,
		                   cbe->priv->pba.inter_normal,
		                   &last_index_back,
		                   pvecback);
		z_i_b = cbe->priv->pba.a_today / pvecback[cbe->priv->pba.index_bg_a] - 1.0;

		for (j = 0; j < partz; j++)
		{
			const gdouble z_i = z_i_a + (z_i_b - z_i_a) / (partz * 1.0) * j;

			ncm_vector_set (z_v, m, z_i);
			spectra_pk_at_z (&cbe->priv->pba, &cbe->priv->psp, logarithmic, z_i, ncm_matrix_ptr (lnPk, m, 0), temp, temp, temp);
      
			m++;
		}
	}

	for (i = 1; i < z_size - 1; i++)
	{
		gdouble ln_tau_i;
		gint last_index_back;
		guint j;

		ln_tau_i = cbe->priv->psp.ln_tau[z_size - i - 1];
		background_at_tau (&cbe->priv->pba,
		                   exp (ln_tau_i),
		                   cbe->priv->pba.short_info,
		                   cbe->priv->pba.inter_normal,
		                   &last_index_back,
		                   pvecback);
		z_i_a = cbe->priv->pba.a_today / pvecback[cbe->priv->pba.index_bg_a] - 1.0;

		ln_tau_i = cbe->priv->psp.ln_tau[z_size - i - 2];
		background_at_tau (&cbe->priv->pba,
		                   exp (ln_tau_i),
		                   cbe->priv->pba.short_info,
		                   cbe->priv->pba.inter_normal,
		                   &last_index_back,
		                   pvecback);
		z_i_b = cbe->priv->pba.a_today / pvecback[cbe->priv->pba.index_bg_a] - 1.0;

		for (j = 0; j < partz; j++)
		{
			const gdouble z_i = z_i_a + (z_i_b - z_i_a) / (partz * 1.0) * j;

			ncm_vector_set (z_v, m, z_i);
			spectra_pk_at_z (&cbe->priv->pba, &cbe->priv->psp, logarithmic, z_i, ncm_matrix_ptr (lnPk, m, 0), temp, temp, temp);

			m++;
		}
	}

	{
		const gdouble z_i = (z_i_b * 0.99 + ncm_vector_get (z_v, m - 1) * 0.01); /* Safeguard against interpolation error on CLASS near the end */

		ncm_vector_set (z_v, m, z_i);
		spectra_pk_at_z (&cbe->priv->pba, &cbe->priv->psp, logarithmic, z_i, ncm_matrix_ptr (lnPk, m, 0), temp, temp, temp);

		m++;
	}

  g_free (temp);
	g_free (pvecback);

	{
		NcmSpline2d* lnPk_s = ncm_spline2d_bicubic_notaknot_new ();
		ncm_spline2d_set (lnPk_s, lnk_v, z_v, lnPk, TRUE);

		ncm_vector_free (z_v);
		ncm_vector_free (lnk_v);
		ncm_matrix_free (lnPk);

		return lnPk_s;
	}
}

/**
 * nc_cbe_get_sigma8:
 * @cbe: a #NcCBE
 * 
 * Computes the value of $\sigma_8$ as computed by CLASS, usually with errors $\propto 10^{-4}$.
 * For better precision use: ncm_powspec_sigma_tophat_R() or #NcmPowspecFilter. 
 * 
 * Returns: the value of $\sigma_8$ as computed by CLASS.
 */
gdouble 
nc_cbe_get_sigma8 (NcCBE *cbe)
{
  return cbe->priv->psp.sigma8;
}

/**
 * nc_cbe_get_all_Cls:
 * @cbe: a #NcCBE
 * @PHIPHI_Cls: a #NcmVector
 * @TT_Cls: a #NcmVector
 * @EE_Cls: a #NcmVector
 * @BB_Cls: a #NcmVector
 * @TE_Cls: a #NcmVector
 *
 * Gets and store the angular power spectra $C_l$'s calculated in the vectors @TT_Cls, @EE_Cls,
 * @BB_Cls and @TE_Cls. If any of these vectors are NULL, then it is (they are)
 * ignored.
 *
 */
void 
nc_cbe_get_all_Cls (NcCBE* cbe, NcmVector* PHIPHI_Cls, NcmVector* TT_Cls, NcmVector* EE_Cls, NcmVector* BB_Cls, NcmVector* TE_Cls)
{
	guint all_Cls_size, index_pp, index_tt, index_ee, index_bb, index_te;
	gboolean has_pp, has_tt, has_ee, has_bb, has_te;

	if (cbe->use_lensed_Cls)
	{
		struct lensing *ptr = &cbe->priv->ple;

    all_Cls_size = ptr->lt_size;
		index_pp     = ptr->index_lt_pp;
		index_tt     = ptr->index_lt_tt;
		index_ee     = ptr->index_lt_ee;
		index_bb     = ptr->index_lt_bb;
		index_te     = ptr->index_lt_te;
		has_pp       = ptr->has_pp;
		has_tt       = ptr->has_tt;
		has_ee       = ptr->has_ee;
		has_bb       = ptr->has_bb;
		has_te       = ptr->has_te;
	}
	else
	{
		struct spectra* ptr = &cbe->priv->psp;

    all_Cls_size = ptr->ct_size;
		index_pp     = ptr->index_ct_pp;
		index_tt     = ptr->index_ct_tt;
		index_ee     = ptr->index_ct_ee;
		index_bb     = ptr->index_ct_bb;
		index_te     = ptr->index_ct_te;
		has_pp       = ptr->has_pp;
		has_tt       = ptr->has_tt;
		has_ee       = ptr->has_ee;
		has_bb       = ptr->has_bb;
		has_te       = ptr->has_te;
	}

	{
		const gdouble T_gamma0 = cbe->priv->pba.T_cmb;
		const gdouble Cl_fac   = gsl_pow_2 (1.0e6 * T_gamma0);
		gdouble *all_Cls       = g_new0 (gdouble, all_Cls_size);
		guint PHIPHI_lmax      = has_pp ? (PHIPHI_Cls != NULL ? ncm_vector_len (PHIPHI_Cls) : 0) : 0;
		guint TT_lmax          = has_tt ? (TT_Cls     != NULL ? ncm_vector_len (TT_Cls)     : 0) : 0;
		guint EE_lmax          = has_ee ? (EE_Cls     != NULL ? ncm_vector_len (EE_Cls)     : 0) : 0;
		guint BB_lmax          = has_bb ? (BB_Cls     != NULL ? ncm_vector_len (BB_Cls)     : 0) : 0;
		guint TE_lmax          = has_te ? (TE_Cls     != NULL ? ncm_vector_len (TE_Cls)     : 0) : 0;
		guint l;

		g_assert (!(cbe->target_Cls & NC_DATA_CMB_TYPE_PHIPHI) || has_pp);
		g_assert (!(cbe->target_Cls & NC_DATA_CMB_TYPE_TT)     || has_tt);
		g_assert (!(cbe->target_Cls & NC_DATA_CMB_TYPE_EE)     || has_ee);
		g_assert (!(cbe->target_Cls & NC_DATA_CMB_TYPE_BB)     || has_bb);
		g_assert (!(cbe->target_Cls & NC_DATA_CMB_TYPE_TE)     || has_te);

		for (l = 0; l <= cbe->scalar_lmax; l++)
		{
			if (cbe->use_lensed_Cls)
				lensing_cl_at_l (&cbe->priv->ple, l, all_Cls);
			else
				spectra_cl_at_l (&cbe->priv->psp, l, all_Cls, NULL, NULL);

			if (PHIPHI_lmax > 0)
			{
				const gdouble PHIPHI_Cl = all_Cls[index_pp];
				ncm_vector_set (PHIPHI_Cls, l, PHIPHI_Cl);
				PHIPHI_lmax--;
			}
			if (TT_lmax > 0)
			{
				const gdouble TT_Cl = all_Cls[index_tt];
				ncm_vector_set (TT_Cls, l, Cl_fac * TT_Cl);
				TT_lmax--;
			}
			if (EE_lmax > 0)
			{
				const gdouble EE_Cl = all_Cls[index_ee];
				ncm_vector_set (EE_Cls, l, Cl_fac * EE_Cl);
				EE_lmax--;
			}
			if (BB_lmax > 0)
			{
				const gdouble BB_Cl = all_Cls[index_bb];
				ncm_vector_set (BB_Cls, l, Cl_fac * BB_Cl);
				BB_lmax--;
			}
			if (TE_lmax > 0)
			{
				const gdouble TE_Cl = all_Cls[index_te];
				ncm_vector_set (TE_Cls, l, Cl_fac * TE_Cl);
				TE_lmax--;
			}
		}

		g_free (all_Cls);
	}
}

/**
 * nc_cbe_debug_test:
 * @cbe: a #NcCBE
 *
 * Temporary debug function
 *
 */
void 
nc_cbe_debug_test (NcCBE *cbe)
{
	NCM_UNUSED (cbe);
	g_assert_not_reached ();
}
