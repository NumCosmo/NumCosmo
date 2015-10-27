/***************************************************************************
 *            nc_hipert_boltzmann_cbe.c
 *
 *  Sat October 24 11:56:56 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_boltzmann_cbe.c
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
 * SECTION:nc_hipert_boltzmann_cbe
 * @title: NcHIPertBoltzmannCBE
 * @short_description: CLASS (Cosmic Linear Anisotropy Solving System) backend for perturbations
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
 *
 */

/*
 * It must be include before anything else, several symbols clash
 * with the default includes.
 */
#include "class/include/class.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hiprim.h"
#include "model/nc_hicosmo_de.h"
#include "nc_hipert_boltzmann_cbe.h"

enum
{
  PROP_0,
  PROP_PREC,
};

struct _NcHIPertBoltzmannCBEPrivate
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

G_DEFINE_TYPE (NcHIPertBoltzmannCBE, nc_hipert_boltzmann_cbe, NC_TYPE_HIPERT_BOLTZMANN);

static void
nc_hipert_boltzmann_cbe_init (NcHIPertBoltzmannCBE *cbe)
{
  cbe->priv   = G_TYPE_INSTANCE_GET_PRIVATE (cbe, NC_TYPE_HIPERT_BOLTZMANN_CBE, NcHIPertBoltzmannCBEPrivate);
  cbe->prec   = NULL;
  cbe->lmax   = 0;
  cbe->TT_Cls = NULL;
  cbe->EE_Cls = NULL;
  cbe->BB_Cls = NULL;
  cbe->TE_Cls = NULL;
  cbe->TB_Cls = NULL;
  cbe->EB_Cls = NULL;

  cbe->priv->pba.h                    = 0.0;
  cbe->priv->pba.H0                   = 0.0;
  cbe->priv->pba.T_cmb                = 0.0;
  cbe->priv->pba.Omega0_g             = 0.0;
  cbe->priv->pba.Omega0_ur            = 0.0;
  cbe->priv->pba.Omega0_b             = 0.0;
  cbe->priv->pba.Omega0_cdm           = 0.0;
  cbe->priv->pba.Omega0_dcdmdr        = 0.0;
  cbe->priv->pba.Omega0_dcdm          = 0.0;
  cbe->priv->pba.Gamma_dcdm           = 0.0;
  cbe->priv->pba.N_ncdm               = 0;
  cbe->priv->pba.Omega0_ncdm_tot      = 0.0;
  cbe->priv->pba.ksi_ncdm_default     = 0.0;
  cbe->priv->pba.ksi_ncdm             = NULL;
  cbe->priv->pba.T_ncdm_default       = 0.0;
  cbe->priv->pba.T_ncdm               = NULL;
  cbe->priv->pba.deg_ncdm_default     = 0.0;
  cbe->priv->pba.deg_ncdm             = NULL;
  cbe->priv->pba.ncdm_psd_parameters  = NULL;
  cbe->priv->pba.ncdm_psd_files       = NULL;
  cbe->priv->pba.Omega0_scf           = 0.0;
  cbe->priv->pba.attractor_ic_scf     = _FALSE_;
  cbe->priv->pba.scf_parameters       = NULL;
  cbe->priv->pba.scf_parameters_size  = 0;
  cbe->priv->pba.scf_tuning_index     = 0;
  cbe->priv->pba.phi_ini_scf          = 0;
  cbe->priv->pba.phi_prime_ini_scf    = 0;
  cbe->priv->pba.Omega0_k             = 0.0;
  cbe->priv->pba.K                    = 0.0;
  cbe->priv->pba.sgnK                 = 0;
  cbe->priv->pba.Omega0_lambda        = 0.0;
  cbe->priv->pba.Omega0_fld           = 0.0;
  cbe->priv->pba.a_today              = 0.0;
  cbe->priv->pba.w0_fld               = 0.0;
  cbe->priv->pba.wa_fld               = 0.0;
  cbe->priv->pba.cs2_fld              = 0.0;

  /* thermodynamics structure */

  cbe->priv->pth.YHe                      = 0;
  cbe->priv->pth.recombination            = 0;
  cbe->priv->pth.reio_parametrization     = 0;
  cbe->priv->pth.reio_z_or_tau            = 0;
  cbe->priv->pth.z_reio                   = 0.0;
  cbe->priv->pth.tau_reio                 = 0.0;
  cbe->priv->pth.reionization_exponent    = 0.0;
  cbe->priv->pth.reionization_width       = 0.0;
  cbe->priv->pth.helium_fullreio_redshift = 0.0;
  cbe->priv->pth.helium_fullreio_width    = 0.0;

  cbe->priv->pth.binned_reio_num            = 0;
  cbe->priv->pth.binned_reio_z              = NULL;
  cbe->priv->pth.binned_reio_xe             = NULL;
  cbe->priv->pth.binned_reio_step_sharpness = 0.3;

  cbe->priv->pth.annihilation           = 0.0;
  cbe->priv->pth.decay                  = 0.0;
  cbe->priv->pth.annihilation_variation = 0.0;
  cbe->priv->pth.annihilation_z         = 0.0;
  cbe->priv->pth.annihilation_zmax      = 0.0;
  cbe->priv->pth.annihilation_zmin      = 0.0;
  cbe->priv->pth.annihilation_f_halo    = 0.0;
  cbe->priv->pth.annihilation_z_halo    = 0.0;
  cbe->priv->pth.has_on_the_spot        = _FALSE_;

  cbe->priv->pth.compute_cb2_derivatives = _FALSE_;

  /* perturbation structure */

  cbe->priv->ppt.has_perturbations            = _FALSE_;
  cbe->priv->ppt.has_cls                      = _FALSE_;

  cbe->priv->ppt.has_cl_cmb_temperature       = _FALSE_;
  cbe->priv->ppt.has_cl_cmb_polarization      = _FALSE_;
  cbe->priv->ppt.has_cl_cmb_lensing_potential = _FALSE_;
  cbe->priv->ppt.has_cl_number_count          = _FALSE_;
  cbe->priv->ppt.has_cl_lensing_potential     = _FALSE_;
  cbe->priv->ppt.has_pk_matter                = _FALSE_;
  cbe->priv->ppt.has_density_transfers        = _FALSE_;
  cbe->priv->ppt.has_velocity_transfers       = _FALSE_;

  cbe->priv->ppt.has_nl_corrections_based_on_delta_m = _FALSE_;

  cbe->priv->ppt.has_nc_density = _FALSE_;
  cbe->priv->ppt.has_nc_rsd     = _FALSE_;
  cbe->priv->ppt.has_nc_lens    = _FALSE_;
  cbe->priv->ppt.has_nc_gr      = _FALSE_;

  cbe->priv->ppt.switch_sw         = 0;
  cbe->priv->ppt.switch_eisw       = 0;
  cbe->priv->ppt.switch_lisw       = 0;
  cbe->priv->ppt.switch_dop        = 0;
  cbe->priv->ppt.switch_pol        = 0;
  cbe->priv->ppt.eisw_lisw_split_z = 0;

  cbe->priv->ppt.has_ad  = _FALSE_;
  cbe->priv->ppt.has_bi  = _FALSE_;
  cbe->priv->ppt.has_cdi = _FALSE_;
  cbe->priv->ppt.has_nid = _FALSE_;
  cbe->priv->ppt.has_niv = _FALSE_;

  cbe->priv->ppt.has_perturbed_recombination = _FALSE_;
  cbe->priv->ppt.tensor_method               = tm_massless_approximation;
  cbe->priv->ppt.evolve_tensor_ur            = _FALSE_;
  cbe->priv->ppt.evolve_tensor_ncdm          = _FALSE_;

  cbe->priv->ppt.has_scalars = _FALSE_;
  cbe->priv->ppt.has_vectors = _FALSE_;
  cbe->priv->ppt.has_tensors = _FALSE_;

  cbe->priv->ppt.l_scalar_max = 0;
  cbe->priv->ppt.l_vector_max = 0;
  cbe->priv->ppt.l_tensor_max = 0;
  cbe->priv->ppt.l_lss_max    = 0;
  cbe->priv->ppt.k_max_for_pk = 0.0;

  cbe->priv->ppt.gauge = synchronous;

  cbe->priv->ppt.k_output_values_num     = 0;
  cbe->priv->ppt.store_perturbations     = _FALSE_;
  cbe->priv->ppt.number_of_scalar_titles = 0;
  cbe->priv->ppt.number_of_vector_titles = 0;
  cbe->priv->ppt.number_of_tensor_titles = 0;
  {
    guint filenum;
    for (filenum = 0; filenum<_MAX_NUMBER_OF_K_FILES_; filenum++){
      cbe->priv->ppt.scalar_perturbations_data[filenum] = NULL;
      cbe->priv->ppt.vector_perturbations_data[filenum] = NULL;
      cbe->priv->ppt.tensor_perturbations_data[filenum] = NULL;
    }
  }
  cbe->priv->ppt.index_k_output_values = NULL;

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
  cbe->priv->ppm.ln_aH_ratio          = 0;
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

  cbe->priv->ppt.selection_num        = 0;
  cbe->priv->ppt.selection            = 0;
  cbe->priv->ppt.selection_mean[0]    = 0.0;
  cbe->priv->ppt.selection_width[0]   = 0.0;

  cbe->priv->ptr.lcmb_rescale         = 0.0;
  cbe->priv->ptr.lcmb_pivot           = 0.0;
  cbe->priv->ptr.lcmb_tilt            = 0.0;
  cbe->priv->ptr.initialise_HIS_cache = _FALSE_;
  cbe->priv->ptr.has_nz_analytic      = _FALSE_;
  cbe->priv->ptr.has_nz_file          = _FALSE_;
  cbe->priv->ptr.has_nz_evo_analytic  = _FALSE_;
  cbe->priv->ptr.has_nz_evo_file      = _FALSE_;
  cbe->priv->ptr.bias                 = 0.0;
  cbe->priv->ptr.s_bias               = 0.0;

  /* spectra structure */

  cbe->priv->psp.z_max_pk = 0.0;
  cbe->priv->psp.non_diag = 0;

  /* lensing structure */

  cbe->priv->ple.has_lensed_cls = _FALSE_;

  /* nonlinear structure */

  cbe->priv->pnl.method = nl_none;

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
    guint verbosity = 0;
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
nc_hipert_boltzmann_cbe_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPertBoltzmannCBE *cbe = NC_HIPERT_BOLTZMANN_CBE (object);
  g_return_if_fail (NC_IS_HIPERT_BOLTZMANN_CBE (object));

  switch (prop_id)
  {
    case PROP_PREC:
      nc_cbe_precision_clear (&cbe->prec);
      cbe->prec = g_value_dup_object (value);
      ncm_model_ctrl_force_update (NC_HIPERT_BOLTZMANN (cbe)->ctrl_cosmo);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_boltzmann_cbe_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIPertBoltzmannCBE *cbe = NC_HIPERT_BOLTZMANN_CBE (object);
  g_return_if_fail (NC_IS_HIPERT_BOLTZMANN_CBE (object));

  switch (prop_id)
  {
    case PROP_PREC:
      g_value_set_object (value, cbe->prec);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_boltzmann_cbe_dispose (GObject *object)
{
  NcHIPertBoltzmannCBE *cbe = NC_HIPERT_BOLTZMANN_CBE (object);

  nc_cbe_precision_clear (&cbe->prec);

  ncm_vector_clear (&cbe->TT_Cls);
  ncm_vector_clear (&cbe->EE_Cls);
  ncm_vector_clear (&cbe->BB_Cls);
  ncm_vector_clear (&cbe->TE_Cls);
  ncm_vector_clear (&cbe->TB_Cls);
  ncm_vector_clear (&cbe->EB_Cls);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_boltzmann_cbe_parent_class)->dispose (object);
}

static void
nc_hipert_boltzmann_cbe_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_boltzmann_cbe_parent_class)->finalize (object);
}

static void _nc_hipert_boltzmann_cbe_prepare (NcHIPertBoltzmann *pb, NcHIPrim *prim, NcHICosmo *cosmo);
static void _nc_hipert_boltzmann_cbe_get_TT_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
static void _nc_hipert_boltzmann_cbe_get_EE_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
static void _nc_hipert_boltzmann_cbe_get_BB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
static void _nc_hipert_boltzmann_cbe_get_TE_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
static void _nc_hipert_boltzmann_cbe_get_TB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
static void _nc_hipert_boltzmann_cbe_get_EB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);

static void
nc_hipert_boltzmann_cbe_class_init (NcHIPertBoltzmannCBEClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  g_type_class_add_private (klass, sizeof (NcHIPertBoltzmannCBEPrivate));

  object_class->set_property = nc_hipert_boltzmann_cbe_set_property;
  object_class->get_property = nc_hipert_boltzmann_cbe_get_property;
  object_class->dispose      = nc_hipert_boltzmann_cbe_dispose;
  object_class->finalize     = nc_hipert_boltzmann_cbe_finalize;

  g_object_class_install_property (object_class,
                                   PROP_PREC,
                                   g_param_spec_object ("precision",
                                                        NULL,
                                                        "CLASS precision object",
                                                        NC_TYPE_CBE_PRECISION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  NC_HIPERT_BOLTZMANN_CLASS (klass)->prepare    = &_nc_hipert_boltzmann_cbe_prepare;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_TT_Cls = &_nc_hipert_boltzmann_cbe_get_TT_Cls;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_EE_Cls = &_nc_hipert_boltzmann_cbe_get_EE_Cls;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_BB_Cls = &_nc_hipert_boltzmann_cbe_get_BB_Cls;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_TE_Cls = &_nc_hipert_boltzmann_cbe_get_TE_Cls;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_TB_Cls = &_nc_hipert_boltzmann_cbe_get_TB_Cls;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_EB_Cls = &_nc_hipert_boltzmann_cbe_get_EB_Cls;
}

static void
_nc_hipert_boltzmann_cbe_get_TT_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  NcHIPertBoltzmannCBE *cbe = NC_HIPERT_BOLTZMANN_CBE (pb);
  ncm_vector_memcpy2 (Cls, cbe->TT_Cls, 0, 0, ncm_vector_len (Cls));
}

static void
_nc_hipert_boltzmann_cbe_get_EE_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  NcHIPertBoltzmannCBE *cbe = NC_HIPERT_BOLTZMANN_CBE (pb);
  ncm_vector_memcpy2 (Cls, cbe->EE_Cls, 0, 0, ncm_vector_len (Cls));
}

static void
_nc_hipert_boltzmann_cbe_get_BB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  NcHIPertBoltzmannCBE *cbe = NC_HIPERT_BOLTZMANN_CBE (pb);
  ncm_vector_memcpy2 (Cls, cbe->BB_Cls, 0, 0, ncm_vector_len (Cls));
}

static void
_nc_hipert_boltzmann_cbe_get_TE_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  NcHIPertBoltzmannCBE *cbe = NC_HIPERT_BOLTZMANN_CBE (pb);
  ncm_vector_memcpy2 (Cls, cbe->TE_Cls, 0, 0, ncm_vector_len (Cls));
}

static void
_nc_hipert_boltzmann_cbe_get_TB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  NcHIPertBoltzmannCBE *cbe = NC_HIPERT_BOLTZMANN_CBE (pb);
  ncm_vector_memcpy2 (Cls, cbe->TB_Cls, 0, 0, ncm_vector_len (Cls));
}

static void
_nc_hipert_boltzmann_cbe_get_EB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  NcHIPertBoltzmannCBE *cbe = NC_HIPERT_BOLTZMANN_CBE (pb);
  ncm_vector_memcpy2 (Cls, cbe->EB_Cls, 0, 0, ncm_vector_len (Cls));
}

/**
 * nc_hipert_boltzmann_cbe_new: (constructor)
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcCBEPrecision.
 */
NcHIPertBoltzmannCBE *
nc_hipert_boltzmann_cbe_new (void)
{
  NcCBEPrecision *prec = nc_cbe_precision_new ();
  NcHIPertBoltzmannCBE *cbe = g_object_new (NC_TYPE_HIPERT_BOLTZMANN_CBE,
                                            "precision", prec,
                                            NULL);
  return cbe;
}

/**
 * nc_hipert_boltzmann_cbe_prec_new: (constructor)
 * @cbe_prec: a #NcCBEPrecision.
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcCBEPrecision.
 */
NcHIPertBoltzmannCBE *
nc_hipert_boltzmann_cbe_prec_new (NcCBEPrecision *cbe_prec)
{
  NcHIPertBoltzmannCBE *cbe = g_object_new (NC_TYPE_HIPERT_BOLTZMANN_CBE,
                                            "precision", cbe_prec,
                                            NULL);
  return cbe;
}

/**
 * nc_hipert_boltzmann_cbe_ref:
 * @cbe: a #NcHIPertBoltzmannCBE.
 *
 * Increases the reference count of @cbe.
 *
 * Returns: (transfer full): @cbe.
 */
NcHIPertBoltzmannCBE *
nc_hipert_boltzmann_cbe_ref (NcHIPertBoltzmannCBE *cbe)
{
  return g_object_ref (cbe);
}

/**
 * nc_hipert_boltzmann_cbe_free:
 * @cbe: a #NcHIPertBoltzmannCBE.
 *
 * Decreases the reference count of @cbe.
 *
 */
void
nc_hipert_boltzmann_cbe_free (NcHIPertBoltzmannCBE *cbe)
{
  g_object_unref (cbe);
}

/**
 * nc_hipert_boltzmann_cbe_clear:
 * @cbe: a #NcHIPertBoltzmannCBE.
 *
 * Decreases the reference count of *@cbe and sets *@cbe to NULL.
 *
 */
void
nc_hipert_boltzmann_cbe_clear (NcHIPertBoltzmannCBE **cbe)
{
  g_clear_object (cbe);
}

static void
_nc_hipert_boltzmann_cbe_set_bg (NcHIPertBoltzmannCBE *cbe, NcHICosmo *cosmo)
{
  if (!g_type_is_a (G_OBJECT_TYPE (cosmo), NC_TYPE_HICOSMO_DE))
    g_error ("_nc_hipert_boltzmann_cbe_set_bg: CLASS backend is compatible with darkenergy models only.");

  cbe->priv->pba.h                   = nc_hicosmo_h (cosmo);
  cbe->priv->pba.H0                  = cbe->priv->pba.h * 1.0e5 / ncm_c_c ();
  cbe->priv->pba.T_cmb               = nc_hicosmo_T_gamma0 (cosmo);
  cbe->priv->pba.Omega0_g            = nc_hicosmo_Omega_g (cosmo);
  cbe->priv->pba.Omega0_ur           = nc_hicosmo_Omega_nu (cosmo);
  cbe->priv->pba.Omega0_b            = nc_hicosmo_Omega_b (cosmo);
  cbe->priv->pba.Omega0_cdm          = nc_hicosmo_Omega_c (cosmo);
  cbe->priv->pba.Omega0_dcdmdr       = 0.0;
  cbe->priv->pba.Omega0_dcdm         = 0.0;
  cbe->priv->pba.Gamma_dcdm          = 0.0;
  cbe->priv->pba.N_ncdm              = 0;
  cbe->priv->pba.Omega0_ncdm_tot     = 0.;
  cbe->priv->pba.ksi_ncdm_default    = 0.;
  cbe->priv->pba.ksi_ncdm            = NULL;
  cbe->priv->pba.T_ncdm_default      = 0.71611;
  cbe->priv->pba.T_ncdm              = NULL;
  cbe->priv->pba.deg_ncdm_default    = 1.0;
  cbe->priv->pba.deg_ncdm            = NULL;
  cbe->priv->pba.ncdm_psd_parameters = NULL;
  cbe->priv->pba.ncdm_psd_files      = NULL;

  cbe->priv->pba.Omega0_scf          = 0.0;
  cbe->priv->pba.attractor_ic_scf    = _TRUE_;
  cbe->priv->pba.scf_parameters      = NULL;
  cbe->priv->pba.scf_parameters_size = 0;
  cbe->priv->pba.scf_tuning_index    = 0;
  cbe->priv->pba.phi_ini_scf         = 1;
  cbe->priv->pba.phi_prime_ini_scf   = 1;

  cbe->priv->pba.Omega0_k            = nc_hicosmo_Omega_k (cosmo);
  if (fabs (cbe->priv->pba.Omega0_k) > 1.0e-13)
  {
    cbe->priv->pba.K                   = -GSL_SIGN (cbe->priv->pba.Omega0_k);
    cbe->priv->pba.sgnK                = cbe->priv->pba.K;
    cbe->priv->pba.a_today             = sqrt (1.0 / fabs (cbe->priv->pba.Omega0_k)) / cbe->priv->pba.H0;
  }
  else
  {
    cbe->priv->pba.Omega0_k            = 0.0;
    cbe->priv->pba.K                   = 0.0;
    cbe->priv->pba.sgnK                = 0;
    cbe->priv->pba.a_today             = 1.0;
  }

  cbe->priv->pba.Omega0_lambda       = ncm_model_orig_param_get (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X) ;
  cbe->priv->pba.Omega0_fld          = 0.0;
  cbe->priv->pba.w0_fld              = -1.0;
  cbe->priv->pba.wa_fld              = 0.0;
  cbe->priv->pba.cs2_fld             = 1.0;

  cbe->priv->pba.background_verbose  = cbe->bg_verbose;
}

static void
_nc_hipert_boltzmann_cbe_set_thermo (NcHIPertBoltzmannCBE *cbe, NcHICosmo *cosmo)
{
  cbe->priv->pth.YHe                      = _BBN_;
  cbe->priv->pth.recombination            = recfast;
  cbe->priv->pth.reio_parametrization     = reio_camb;
  cbe->priv->pth.reio_z_or_tau            = reio_z;
  cbe->priv->pth.z_reio                   = 11.357;
  cbe->priv->pth.tau_reio                 = 0.0925;
  cbe->priv->pth.reionization_exponent    = 1.5;
  cbe->priv->pth.reionization_width       = 0.5;
  cbe->priv->pth.helium_fullreio_redshift = 3.5;
  cbe->priv->pth.helium_fullreio_width    = 0.5;

  cbe->priv->pth.binned_reio_num            = 0;
  cbe->priv->pth.binned_reio_z              = NULL;
  cbe->priv->pth.binned_reio_xe             = NULL;
  cbe->priv->pth.binned_reio_step_sharpness = 0.3;

  cbe->priv->pth.annihilation           = 0.0;
  cbe->priv->pth.decay                  = 0.0;
  cbe->priv->pth.annihilation_variation = 0.0;
  cbe->priv->pth.annihilation_z         = 1000.0;
  cbe->priv->pth.annihilation_zmax      = 2500.0;
  cbe->priv->pth.annihilation_zmin      = 30.0;
  cbe->priv->pth.annihilation_f_halo    = 0.0;
  cbe->priv->pth.annihilation_z_halo    = 30.0;
  cbe->priv->pth.has_on_the_spot        = _TRUE_;

  cbe->priv->pth.compute_cb2_derivatives = _FALSE_;

  cbe->priv->pth.thermodynamics_verbose  = cbe->thermo_verbose;
}

static void
_nc_hipert_boltzmann_cbe_set_pert (NcHIPertBoltzmannCBE *cbe, NcHICosmo *cosmo)
{
  NcHIPertBoltzmann *pb = NC_HIPERT_BOLTZMANN (cbe);
  gboolean has_cls = (pb->target_Cls & NC_DATA_CMB_TYPE_ALL) != 0;
  gboolean has_perturbations = has_cls || pb->calc_transfer;

  /*
   * Inside CLASS they compare booleans with _TRUE_ and _FALSE_.
   * This is a bad ideia, but to be compatible we must always
   * use their _TRUE_ and _FALSE_.
   */

  if (pb->target_Cls & (NC_DATA_CMB_TYPE_TB | NC_DATA_CMB_TYPE_EB))
    g_error ("_nc_hipert_boltzmann_cbe_set_pert: modes TB and EB are not supported.");

  cbe->priv->ppt.has_perturbations            = has_perturbations ? _TRUE_ : _FALSE_;
  cbe->priv->ppt.has_cls                      = has_cls ? _TRUE_ : _FALSE_;

  cbe->priv->ppt.has_cl_cmb_temperature       = pb->target_Cls & NC_DATA_CMB_TYPE_TT ? _TRUE_ : _FALSE_;
  cbe->priv->ppt.has_cl_cmb_polarization      =
    pb->target_Cls & (NC_DATA_CMB_TYPE_EE | NC_DATA_CMB_TYPE_BB | NC_DATA_CMB_TYPE_TE) ? _TRUE_ : _FALSE_;
  cbe->priv->ppt.has_cl_cmb_lensing_potential = _FALSE_;
  cbe->priv->ppt.has_cl_number_count          = _FALSE_;
  cbe->priv->ppt.has_cl_lensing_potential     = _FALSE_;
  cbe->priv->ppt.has_pk_matter                = pb->calc_transfer ? _TRUE_ : _FALSE_;
  cbe->priv->ppt.has_density_transfers        = _FALSE_;
  cbe->priv->ppt.has_velocity_transfers       = _FALSE_;

  cbe->priv->ppt.has_nl_corrections_based_on_delta_m = _FALSE_;

  cbe->priv->ppt.has_nc_density = _FALSE_;
  cbe->priv->ppt.has_nc_rsd     = _FALSE_;
  cbe->priv->ppt.has_nc_lens    = _FALSE_;
  cbe->priv->ppt.has_nc_gr      = _FALSE_;

  cbe->priv->ppt.switch_sw         = 1;
  cbe->priv->ppt.switch_eisw       = 1;
  cbe->priv->ppt.switch_lisw       = 1;
  cbe->priv->ppt.switch_dop        = 1;
  cbe->priv->ppt.switch_pol        = 1;
  cbe->priv->ppt.eisw_lisw_split_z = 120;

  cbe->priv->ppt.has_ad  = _TRUE_;
  cbe->priv->ppt.has_bi  = _FALSE_;
  cbe->priv->ppt.has_cdi = _FALSE_;
  cbe->priv->ppt.has_nid = _FALSE_;
  cbe->priv->ppt.has_niv = _FALSE_;

  cbe->priv->ppt.has_perturbed_recombination = _FALSE_;
  cbe->priv->ppt.tensor_method               = tm_massless_approximation;
  cbe->priv->ppt.evolve_tensor_ur            = _FALSE_;
  cbe->priv->ppt.evolve_tensor_ncdm          = _FALSE_;

  cbe->priv->ppt.has_scalars = _TRUE_;
  cbe->priv->ppt.has_vectors = _FALSE_;
  cbe->priv->ppt.has_tensors = _FALSE_;

  {
    guint TT_lmax = nc_hipert_boltzmann_get_TT_lmax (pb);
    guint EE_lmax = nc_hipert_boltzmann_get_EE_lmax (pb);
    guint BB_lmax = nc_hipert_boltzmann_get_BB_lmax (pb);
    guint TE_lmax = nc_hipert_boltzmann_get_TE_lmax (pb);
    guint TB_lmax = nc_hipert_boltzmann_get_TB_lmax (pb);
    guint EB_lmax = nc_hipert_boltzmann_get_EB_lmax (pb);

#define _CHECK_VEC(name) \
G_STMT_START { \
 if (cbe->name##_Cls != NULL) \
 { \
   if (ncm_vector_len (cbe->name##_Cls) != name##_lmax + 1) \
   { \
     ncm_vector_clear (&cbe->name##_Cls); \
     cbe->name##_Cls = ncm_vector_new (name##_lmax + 1); \
   } \
 } \
 else \
   cbe->name##_Cls = ncm_vector_new (name##_lmax + 1); \
} G_STMT_END

    _CHECK_VEC (TT);
    _CHECK_VEC (EE);
    _CHECK_VEC (BB);
    _CHECK_VEC (TE);
    _CHECK_VEC (TB);
    _CHECK_VEC (EB);

    cbe->lmax = 0;
    cbe->lmax = GSL_MAX (cbe->lmax, TT_lmax);
    cbe->lmax = GSL_MAX (cbe->lmax, EE_lmax);
    cbe->lmax = GSL_MAX (cbe->lmax, BB_lmax);
    cbe->lmax = GSL_MAX (cbe->lmax, TE_lmax);
    cbe->lmax = GSL_MAX (cbe->lmax, TB_lmax);
    cbe->lmax = GSL_MAX (cbe->lmax, EB_lmax);

    cbe->priv->ppt.l_scalar_max = cbe->lmax;
    cbe->priv->ppt.l_vector_max = 500;
    cbe->priv->ppt.l_tensor_max = 500;
    cbe->priv->ppt.l_lss_max    = 300;
    cbe->priv->ppt.k_max_for_pk = 0.1;
  }

  cbe->priv->ppt.gauge = synchronous;

  cbe->priv->ppt.k_output_values_num     = 0;
  cbe->priv->ppt.store_perturbations     = _FALSE_;
  cbe->priv->ppt.number_of_scalar_titles = 0;
  cbe->priv->ppt.number_of_vector_titles = 0;
  cbe->priv->ppt.number_of_tensor_titles = 0;
  {
    guint filenum;
    for (filenum = 0; filenum<_MAX_NUMBER_OF_K_FILES_; filenum++)
    {
      cbe->priv->ppt.scalar_perturbations_data[filenum] = NULL;
      cbe->priv->ppt.vector_perturbations_data[filenum] = NULL;
      cbe->priv->ppt.tensor_perturbations_data[filenum] = NULL;
    }
  }
  cbe->priv->ppt.index_k_output_values = NULL;

  cbe->priv->ppt.selection_num      = 1;
  cbe->priv->ppt.selection          = gaussian;
  cbe->priv->ppt.selection_mean[0]  = 1.0;
  cbe->priv->ppt.selection_width[0] = 0.1;

  cbe->priv->ppt.perturbations_verbose = cbe->pert_verbose;
}

gdouble
_external_Pk_callback_pks (const double lnk, gpointer data)
{
  NcHIPrim *prim = NC_HIPRIM (data);

  return nc_hiprim_lnSA_powspec_lnk (prim, lnk);
}

gdouble
_external_Pk_callback_pkt (const double lnk, gpointer data)
{
  NcHIPrim *prim = NC_HIPRIM (data);

  return nc_hiprim_lnT_powspec_lnk (prim, lnk);
}

static void
_nc_hipert_boltzmann_cbe_set_prim (NcHIPertBoltzmannCBE *cbe, NcHIPrim *prim, NcHICosmo *cosmo)
{
  cbe->priv->ppm.primordial_spec_type = external_Pk_callback;
  /*cbe->priv->ppm.primordial_spec_type = analytic_Pk;*/
  cbe->priv->ppm.external_Pk_callback_pks  = &_external_Pk_callback_pks;
  /*cbe->priv->ppm.external_Pk_callback_pkt  = &_external_Pk_callback_pkt;*/
  cbe->priv->ppm.external_Pk_callback_data = prim;

  cbe->priv->ppm.k_pivot       = 0.05;
  cbe->priv->ppm.A_s           = 2.215e-9;
  cbe->priv->ppm.n_s           = 0.9619;
  cbe->priv->ppm.alpha_s       = 0.0;
  cbe->priv->ppm.f_bi          = 1.0;
  cbe->priv->ppm.n_bi          = 1.0;
  cbe->priv->ppm.alpha_bi      = 0.0;
  cbe->priv->ppm.f_cdi         = 1.0;
  cbe->priv->ppm.n_cdi         = 1.0;
  cbe->priv->ppm.alpha_cdi     = 0.0;
  cbe->priv->ppm.f_nid         = 1.0;
  cbe->priv->ppm.n_nid         = 1.0;
  cbe->priv->ppm.alpha_nid     = 0.0;
  cbe->priv->ppm.f_niv         = 1.0;
  cbe->priv->ppm.n_niv         = 1.0;
  cbe->priv->ppm.alpha_niv     = 0.0;
  cbe->priv->ppm.c_ad_bi       = 0.0;
  cbe->priv->ppm.n_ad_bi       = 0.0;
  cbe->priv->ppm.alpha_ad_bi   = 0.0;
  cbe->priv->ppm.c_ad_cdi      = 0.0;
  cbe->priv->ppm.n_ad_cdi      = 0.0;
  cbe->priv->ppm.alpha_ad_cdi  = 0.0;
  cbe->priv->ppm.c_ad_nid      = 0.0;
  cbe->priv->ppm.n_ad_nid      = 0.0;
  cbe->priv->ppm.alpha_ad_nid  = 0.0;
  cbe->priv->ppm.c_ad_niv      = 0.0;
  cbe->priv->ppm.n_ad_niv      = 0.0;
  cbe->priv->ppm.alpha_ad_niv  = 0.0;
  cbe->priv->ppm.c_bi_cdi      = 0.0;
  cbe->priv->ppm.n_bi_cdi      = 0.0;
  cbe->priv->ppm.alpha_bi_cdi  = 0.0;
  cbe->priv->ppm.c_bi_nid      = 0.0;
  cbe->priv->ppm.n_bi_nid      = 0.0;
  cbe->priv->ppm.alpha_bi_nid  = 0.0;
  cbe->priv->ppm.c_bi_niv      = 0.0;
  cbe->priv->ppm.n_bi_niv      = 0.0;
  cbe->priv->ppm.alpha_bi_niv  = 0.0;
  cbe->priv->ppm.c_cdi_nid     = 0.0;
  cbe->priv->ppm.n_cdi_nid     = 0.0;
  cbe->priv->ppm.alpha_cdi_nid = 0.0;
  cbe->priv->ppm.c_cdi_niv     = 0.0;
  cbe->priv->ppm.n_cdi_niv     = 0.0;
  cbe->priv->ppm.alpha_cdi_niv = 0.0;
  cbe->priv->ppm.c_nid_niv     = 0.0;
  cbe->priv->ppm.n_nid_niv     = 0.0;
  cbe->priv->ppm.alpha_nid_niv = 0.0;
  cbe->priv->ppm.r           = 1.0;
  cbe->priv->ppm.n_t         = -cbe->priv->ppm.r / 8.0 * (2.0 - cbe->priv->ppm.r / 8.0 - cbe->priv->ppm.n_s);
  cbe->priv->ppm.alpha_t     = cbe->priv->ppm.r / 8.0 * (cbe->priv->ppm.r / 8.0 + cbe->priv->ppm.n_s - 1.0);
  cbe->priv->ppm.potential   = polynomial;
  cbe->priv->ppm.phi_end     = 0.0;
  cbe->priv->ppm.ln_aH_ratio = 50;
  cbe->priv->ppm.V0          = 1.25e-13;
  cbe->priv->ppm.V1          = -1.12e-14;
  cbe->priv->ppm.V2          = -6.95e-14;
  cbe->priv->ppm.V3          = 0.0;
  cbe->priv->ppm.V4          = 0.0;
  cbe->priv->ppm.H0          = 3.69e-6;
  cbe->priv->ppm.H1          = -5.84e-7;
  cbe->priv->ppm.H2          = 0.0;
  cbe->priv->ppm.H3          = 0.0;
  cbe->priv->ppm.H4          = 0.0;
  cbe->priv->ppm.command     = "write here your command for the external Pk";
  cbe->priv->ppm.custom1     = 0.0;
  cbe->priv->ppm.custom2     = 0.0;
  cbe->priv->ppm.custom3     = 0.0;
  cbe->priv->ppm.custom4     = 0.0;
  cbe->priv->ppm.custom5     = 0.0;
  cbe->priv->ppm.custom6     = 0.0;
  cbe->priv->ppm.custom7     = 0.0;
  cbe->priv->ppm.custom8     = 0.0;
  cbe->priv->ppm.custom9     = 0.0;
  cbe->priv->ppm.custom10    = 0.0;

  cbe->priv->ppm.primordial_verbose = cbe->prim_verbose;
}

static void
_nc_hipert_boltzmann_cbe_set_transfer (NcHIPertBoltzmannCBE *cbe, NcHICosmo *cosmo)
{
  cbe->priv->ptr.lcmb_rescale         = 1.0;
  cbe->priv->ptr.lcmb_pivot           = 0.1;
  cbe->priv->ptr.lcmb_tilt            = 0.0;
  cbe->priv->ptr.initialise_HIS_cache = _FALSE_;
  cbe->priv->ptr.has_nz_analytic      = _FALSE_;
  cbe->priv->ptr.has_nz_file          = _FALSE_;
  cbe->priv->ptr.has_nz_evo_analytic  = _FALSE_;
  cbe->priv->ptr.has_nz_evo_file      = _FALSE_;
  cbe->priv->ptr.bias                 = 1.0;
  cbe->priv->ptr.s_bias               = 0.0;

  cbe->priv->ptr.transfer_verbose = cbe->transfer_verbose;
}

static void
_nc_hipert_boltzmann_cbe_set_spectra (NcHIPertBoltzmannCBE *cbe, NcHICosmo *cosmo)
{
  cbe->priv->psp.z_max_pk = 0.0;
  cbe->priv->psp.non_diag = 0;

  cbe->priv->psp.spectra_verbose = cbe->spectra_verbose;
}

static void
_nc_hipert_boltzmann_cbe_set_lensing (NcHIPertBoltzmannCBE *cbe, NcHICosmo *cosmo)
{
  cbe->priv->ple.has_lensed_cls  = _FALSE_;
  cbe->priv->ple.lensing_verbose = cbe->lensing_verbose;
}

static void
_nc_hipert_boltzmann_cbe_set_nonlin (NcHIPertBoltzmannCBE *cbe, NcHICosmo *cosmo)
{
  cbe->priv->pnl.method = nl_none;
  cbe->priv->pnl.nonlinear_verbose = cbe->nonlin_verbose;

}

static void
_nc_hipert_boltzmann_cbe_prepare (NcHIPertBoltzmann *pb, NcHIPrim *prim, NcHICosmo *cosmo)
{
  NcHIPertBoltzmannCBE *cbe = NC_HIPERT_BOLTZMANN_CBE (pb);
  _nc_hipert_boltzmann_cbe_set_bg (cbe, cosmo);
  _nc_hipert_boltzmann_cbe_set_thermo (cbe, cosmo);
  _nc_hipert_boltzmann_cbe_set_pert (cbe, cosmo);
  _nc_hipert_boltzmann_cbe_set_prim (cbe, prim, cosmo);
  _nc_hipert_boltzmann_cbe_set_transfer (cbe, cosmo);
  _nc_hipert_boltzmann_cbe_set_spectra (cbe, cosmo);
  _nc_hipert_boltzmann_cbe_set_lensing (cbe, cosmo);
  _nc_hipert_boltzmann_cbe_set_nonlin (cbe, cosmo);

  if (background_init ((struct precision *)cbe->prec->priv, &cbe->priv->pba) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running background_init `%s'\n", cbe->priv->pba.error_message);

  if (thermodynamics_init ((struct precision *)cbe->prec->priv, &cbe->priv->pba, &cbe->priv->pth) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running thermodynamics_init `%s'\n", cbe->priv->pth.error_message);

  if (perturb_init ((struct precision *)cbe->prec->priv, &cbe->priv->pba, &cbe->priv->pth, &cbe->priv->ppt) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running perturb_init `%s'\n", cbe->priv->ppt.error_message);

  if (primordial_init ((struct precision *)cbe->prec->priv, &cbe->priv->ppt, &cbe->priv->ppm) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running primordial_init `%s'\n", cbe->priv->ppm.error_message);

  if (nonlinear_init ((struct precision *)cbe->prec->priv, &cbe->priv->pba, &cbe->priv->pth, &cbe->priv->ppt, &cbe->priv->ppm, &cbe->priv->pnl) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running nonlinear_init `%s'\n", cbe->priv->pnl.error_message);

  if (transfer_init ((struct precision *)cbe->prec->priv, &cbe->priv->pba, &cbe->priv->pth, &cbe->priv->ppt, &cbe->priv->pnl, &cbe->priv->ptr) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running transfer_init `%s'\n", cbe->priv->ptr.error_message);

  if (spectra_init ((struct precision *)cbe->prec->priv, &cbe->priv->pba, &cbe->priv->ppt, &cbe->priv->ppm, &cbe->priv->pnl, &cbe->priv->ptr, &cbe->priv->psp) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running spectra_init `%s'\n", cbe->priv->psp.error_message);

  if (lensing_init ((struct precision *)cbe->prec->priv, &cbe->priv->ppt, &cbe->priv->psp, &cbe->priv->pnl, &cbe->priv->ple) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running lensing_init `%s'\n", cbe->priv->ple.error_message);

  {
    guint TT_lmax = cbe->priv->psp.has_tt ? nc_hipert_boltzmann_get_TT_lmax (pb) : 0;
    guint EE_lmax = cbe->priv->psp.has_ee ? nc_hipert_boltzmann_get_EE_lmax (pb) : 0;
    guint BB_lmax = cbe->priv->psp.has_bb ? nc_hipert_boltzmann_get_BB_lmax (pb) : 0;
    guint TE_lmax = cbe->priv->psp.has_te ? nc_hipert_boltzmann_get_TE_lmax (pb) : 0;
    const gdouble T_gamma0 = nc_hicosmo_T_gamma0 (cosmo);
    const gdouble Cl_fac   = gsl_pow_2 (1.0e6 * T_gamma0);
    gdouble *all_Cls       = g_new0 (gdouble, cbe->priv->psp.ct_size);
    guint l;

    g_assert (!(pb->target_Cls & NC_DATA_CMB_TYPE_TT) || cbe->priv->psp.has_tt);
    g_assert (!(pb->target_Cls & NC_DATA_CMB_TYPE_EE) || cbe->priv->psp.has_ee);
    g_assert (!(pb->target_Cls & NC_DATA_CMB_TYPE_BB) || cbe->priv->psp.has_bb);
    g_assert (!(pb->target_Cls & NC_DATA_CMB_TYPE_TE) || cbe->priv->psp.has_te);

    for (l = 2; l <= cbe->lmax; l++)
    {
      spectra_cl_at_l (&cbe->priv->psp, l, all_Cls, NULL, NULL);
      /*printf ("A %d %u %u %u %u\n", l, TT_lmax, EE_lmax, BB_lmax, TE_lmax);fflush (stdout);*/
      if (TT_lmax > 1)
      {
        const gdouble TT_Cl = all_Cls[cbe->priv->psp.index_ct_tt];
        ncm_vector_set (cbe->TT_Cls, l, Cl_fac * TT_Cl);
        TT_lmax--;
      }
      /*printf ("B %d\n", l);fflush (stdout);*/
      if (EE_lmax > 1)
      {
        const gdouble EE_Cl = all_Cls[cbe->priv->psp.index_ct_ee];
        ncm_vector_set (cbe->EE_Cls, l, Cl_fac * EE_Cl);
        EE_lmax--;
      }
      /*printf ("C %d\n", l);fflush (stdout);*/
      if (BB_lmax > 1)
      {
        const gdouble BB_Cl = all_Cls[cbe->priv->psp.index_ct_bb];
        ncm_vector_set (cbe->BB_Cls, l, Cl_fac * BB_Cl);
        BB_lmax--;
      }
      /*printf ("D %d\n", l);fflush (stdout);*/
      if (TE_lmax > 1)
      {
        const gdouble TE_Cl = all_Cls[cbe->priv->psp.index_ct_te];
        ncm_vector_set (cbe->TE_Cls, l, Cl_fac * TE_Cl);
        TE_lmax--;
      }
      /*printf ("E %d\n", l);fflush (stdout);*/
    }
    g_free (all_Cls);
  }

  if (lensing_free (&cbe->priv->ple) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running lensing_free `%s'\n", cbe->priv->ple.error_message);

  if (spectra_free (&cbe->priv->psp) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running spectra_free `%s'\n", cbe->priv->ple.error_message);

  if (transfer_free (&cbe->priv->ptr) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running transfer_free `%s'\n", cbe->priv->ple.error_message);

  if (nonlinear_free (&cbe->priv->pnl) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running nonlinear_free `%s'\n", cbe->priv->ple.error_message);

  if (primordial_free (&cbe->priv->ppm) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running primordial_free `%s'\n", cbe->priv->ple.error_message);

  if (perturb_free (&cbe->priv->ppt) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running perturb_free `%s'\n", cbe->priv->ple.error_message);

  if (thermodynamics_free (&cbe->priv->pth) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running thermodynamics_free `%s'\n", cbe->priv->ple.error_message);

  if (background_free (&cbe->priv->pba) == _FAILURE_)
    g_error ("nc_hipert_boltzmann_cbe_prepare: Error running background_free `%s'\n", cbe->priv->ple.error_message);
}
