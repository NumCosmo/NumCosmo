/***************************************************************************
 *            numcosmo.h
 *
 *  Sun May  6 17:20:29 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#ifndef _NUMCOSMO_H
#define _NUMCOSMO_H

/* Supported libraries and build configuration */
#include <numcosmo/numcosmo-math.h>

/* Macros, Constants and Enums */
#include <numcosmo/nc_enum_types.h>

/* Base types and components */
#include <numcosmo/nc/background/nc_hicosmo.h>
#include <numcosmo/nc/primordial/nc_hiprim.h>
#include <numcosmo/nc/background/nc_distance.h>
#include <numcosmo/nc/background/nc_hicosmo_priors.h>
#include <numcosmo/nc/powspec/nc_powspec_ml.h>
#include <numcosmo/nc/powspec/nc_powspec_ml_spline.h>
#include <numcosmo/nc/powspec/nc_powspec_ml_transfer.h>
#include <numcosmo/nc/powspec/nc_powspec_ml_cbe.h>
#include <numcosmo/nc/powspec/nc_powspec_mnl.h>
#include <numcosmo/nc/powspec/nc_powspec_mnl_halofit.h>
#include <numcosmo/nc/supernova/nc_snia_dist_cov.h>
#include <numcosmo/nc/cmb/nc_planck_fi.h>
#include <numcosmo/nc/cmb/nc_planck_fi_cor_tt.h>
#include <numcosmo/nc/cmb/nc_planck_fi_cor_ttteee.h>
#include <numcosmo/nc/background/nc_scalefactor.h>
#include <numcosmo/nc/cmb/nc_cbe_precision.h>
#include <numcosmo/nc/quantum/nc_hiqg_1d.h>

/* Cosmic thermodynamics */
#include <numcosmo/nc/recomb/nc_recomb.h>
#include <numcosmo/nc/recomb/nc_recomb_seager.h>
#include <numcosmo/nc/reion/nc_hireion.h>
#include <numcosmo/nc/reion/nc_hireion_camb.h>

/* Perturbations */
#include <numcosmo/nc/perturbations/nc_hipert_adiab.h>
#include <numcosmo/nc/perturbations/nc_hipert_bg_var.h>
#include <numcosmo/nc/perturbations/nc_hipert_boltzmann_cbe.h>
#include <numcosmo/nc/perturbations/nc_hipert_boltzmann_std.h>
#include <numcosmo/nc/perturbations/nc_hipert_boltzmann.h>
#include <numcosmo/nc/perturbations/nc_hipert_comp_pb.h>
#include <numcosmo/nc/perturbations/nc_hipert_comp.h>
#include <numcosmo/nc/perturbations/nc_hipert_em.h>
#include <numcosmo/nc/perturbations/nc_hipert_first_order.h>
#include <numcosmo/nc/perturbations/nc_hipert_grav_einstein.h>
#include <numcosmo/nc/perturbations/nc_hipert_grav.h>
#include <numcosmo/nc/perturbations/nc_hipert_gw.h>
#include <numcosmo/nc/perturbations/nc_hipert_itwo_fluids.h>
#include <numcosmo/nc/perturbations/nc_hipert_two_fluids.h>
#include <numcosmo/nc/perturbations/nc_hipert.h>

/* Model implementations */
#include <numcosmo/nc/quantum/nc_de_cont.h>
#include <numcosmo/nc/background/nc_hicosmo_idem2.h>
#include <numcosmo/nc/background/nc_hicosmo_gcg.h>
#include <numcosmo/nc/background/nc_hicosmo_de.h>
#include <numcosmo/nc/background/nc_hicosmo_de_reparam_ok.h>
#include <numcosmo/nc/background/nc_hicosmo_de_reparam_cmb.h>
#include <numcosmo/nc/background/nc_hicosmo_de_cpl.h>
#include <numcosmo/nc/background/nc_hicosmo_de_jbp.h>
#include <numcosmo/nc/background/nc_hicosmo_de_xcdm.h>
#include <numcosmo/nc/background/nc_hicosmo_de_wspline.h>
#include <numcosmo/nc/background/nc_hicosmo_lcdm.h>
#include <numcosmo/nc/background/nc_hicosmo_qconst.h>
#include <numcosmo/nc/background/nc_hicosmo_qlinear.h>
#include <numcosmo/nc/background/nc_hicosmo_qspline.h>
#include <numcosmo/nc/background/nc_hicosmo_qrbf.h>
#include <numcosmo/nc/background/nc_hicosmo_qgrw.h>
#include <numcosmo/nc/background/nc_hicosmo_qgw.h>
#include <numcosmo/nc/background/nc_hicosmo_Vexp.h>
#include <numcosmo/nc/primordial/nc_hiprim_power_law.h>
#include <numcosmo/nc/primordial/nc_hiprim_atan.h>
#include <numcosmo/nc/primordial/nc_hiprim_expc.h>
#include <numcosmo/nc/primordial/nc_hiprim_bpl.h>
#include <numcosmo/nc/primordial/nc_hiprim_sbpl.h>
#include <numcosmo/nc/primordial/nc_hiprim_two_fluids.h>

/* Large Scale Structure / Structure Formation */
#include <numcosmo/nc/powspec/nc_window.h>
#include <numcosmo/nc/powspec/nc_window_tophat.h>
#include <numcosmo/nc/powspec/nc_window_gaussian.h>
#include <numcosmo/nc/powspec/nc_transfer_func.h>
#include <numcosmo/nc/powspec/nc_transfer_func_bbks.h>
#include <numcosmo/nc/powspec/nc_transfer_func_eh.h>
#include <numcosmo/nc/powspec/nc_transfer_func_eh_no_baryon.h>
#include <numcosmo/nc/powspec/nc_transfer_func_camb.h>
#include <numcosmo/nc/powspec/nc_growth_func.h>
#include <numcosmo/nc/lss/halo/nc_halo_catalog.h>
#include <numcosmo/nc/lss/halo/nc_halo_catalog_generator.h>
#include <numcosmo/nc/lss/halo/nc_halo_catalog_member_generator.h>
#include <numcosmo/nc/lss/halo/nc_halo_position.h>
#include <numcosmo/nc/lss/halo/nc_halo_density_profile.h>
#include <numcosmo/nc/lss/halo/nc_halo_density_profile_nfw.h>
#include <numcosmo/nc/lss/halo/nc_halo_density_profile_einasto.h>
#include <numcosmo/nc/lss/halo/nc_halo_density_profile_dk14.h>
#include <numcosmo/nc/lss/halo/nc_halo_density_profile_hernquist.h>
#include <numcosmo/nc/lss/halo/nc_halo_mass_summary.h>
#include <numcosmo/nc/lss/halo/nc_halo_cm_param.h>
#include <numcosmo/nc/lss/halo/nc_halo_cm_duffy08.h>
#include <numcosmo/nc/lss/halo/nc_halo_cm_klypin11.h>
#include <numcosmo/nc/lss/halo/nc_halo_cm_prada12.h>
#include <numcosmo/nc/lss/halo/nc_halo_cm_bhattacharya13.h>
#include <numcosmo/nc/lss/halo/nc_halo_cm_dutton14.h>
#include <numcosmo/nc/lss/halo/nc_halo_cm_diemer15.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_acf.h>
#include <numcosmo/nc/lss/halo/nc_multiplicity_func.h>
#include <numcosmo/nc/lss/halo/nc_multiplicity_func_ps.h>
#include <numcosmo/nc/lss/halo/nc_multiplicity_func_st.h>
#include <numcosmo/nc/lss/halo/nc_multiplicity_func_jenkins.h>
#include <numcosmo/nc/lss/halo/nc_multiplicity_func_warren.h>
#include <numcosmo/nc/lss/halo/nc_multiplicity_func_tinker.h>
#include <numcosmo/nc/lss/halo/nc_multiplicity_func_tinker_mean_normalized.h>
#include <numcosmo/nc/lss/halo/nc_multiplicity_func_crocce.h>
#include <numcosmo/nc/lss/halo/nc_multiplicity_func_bocquet.h>
#include <numcosmo/nc/lss/halo/nc_multiplicity_func_bhattacharya.h>
#include <numcosmo/nc/lss/halo/nc_multiplicity_func_despali.h>
#include <numcosmo/nc/lss/halo/nc_multiplicity_func_watson.h>
#include <numcosmo/nc/lss/halo/nc_halo_mass_function.h>
#include <numcosmo/nc/lss/halo/nc_halo_bias.h>
#include <numcosmo/nc/lss/halo/nc_halo_bias_despali.h>
#include <numcosmo/nc/lss/halo/nc_halo_bias_ps.h>
#include <numcosmo/nc/lss/halo/nc_halo_bias_st_spher.h>
#include <numcosmo/nc/lss/halo/nc_halo_bias_st_ellip.h>
#include <numcosmo/nc/lss/halo/nc_halo_bias_tinker.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_redshift.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_redshift_nodist.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_photoz_gauss.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_photoz_gauss_global.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_mass.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_mass_nodist.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_mass_lnnormal.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_mass_vanderlinde.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_mass_benson.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_mass_benson_xray.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_mass_plcl.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_mass_richness.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_mass_ascaso.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_mass_ext.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_mass_selection.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_abundance.h>
#include <numcosmo/nc/lss/cluster/nc_cluster_pseudo_counts.h>
#include <numcosmo/nc/lss/cluster/nc_cor_cluster_cmb_lens_limber.h>
#include <numcosmo/nc/lss/wl/nc_wl_surface_mass_density.h>
#include <numcosmo/nc/lss/wl/nc_reduced_shear_cluster_mass.h>
#include <numcosmo/nc/lss/wl/nc_reduced_shear_calib.h>
#include <numcosmo/nc/lss/wl/nc_reduced_shear_calib_wtg.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_selfunc.h>

/* Galaxy / Galaxy sample distributions */
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_sd_position.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_sd_position_flat.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_sd_obs_redshift.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_sd_obs_redshift_spec.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_sd_obs_redshift_gauss.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_sd_obs_redshift_pz.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_hod.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_hod_zheng07.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_sd_true_redshift.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_sd_true_redshift_lsst_srd.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_sd_shape.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_sd_shape_hsm_gauss.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_sd_shape_hsm_gauss_global.h>

/* Observable data */
#include <numcosmo/nc/data/nc_data_snia.h>
#include <numcosmo/nc/data/nc_data_dist_mu.h>
#include <numcosmo/nc/data/nc_data_snia_cov.h>
#include <numcosmo/nc/data/nc_data_hubble.h>
#include <numcosmo/nc/data/nc_data_hubble_bao.h>
#include <numcosmo/nc/data/nc_data_bao_a.h>
#include <numcosmo/nc/data/nc_data_bao_dv.h>
#include <numcosmo/nc/data/nc_data_bao_rdv.h>
#include <numcosmo/nc/data/nc_data_bao_dvdv.h>
#include <numcosmo/nc/data/nc_data_bao_empirical_fit.h>
#include <numcosmo/nc/data/nc_data_bao_empirical_fit_2d.h>
#include <numcosmo/nc/data/nc_data_bao_dhr_dar.h>
#include <numcosmo/nc/data/nc_data_bao_dtr_dhr.h>
#include <numcosmo/nc/data/nc_data_bao_dmr_hr.h>
#include <numcosmo/nc/data/nc_data_bao_dvr_dtdh.h>
#include <numcosmo/nc/data/nc_data_bao.h>
#include <numcosmo/nc/data/nc_data_cmb_dist_priors.h>
#include <numcosmo/nc/data/nc_data_cmb_shift_param.h>
#include <numcosmo/nc/data/nc_data_cmb.h>
#include <numcosmo/nc/data/nc_data_cluster_ncount.h>
#include <numcosmo/nc/data/nc_data_cluster_ncounts_gauss.h>
#include <numcosmo/nc/data/nc_data_cluster_pseudo_counts.h>
#include <numcosmo/nc/data/nc_data_cluster_wl.h>
#include <numcosmo/nc/data/nc_data_cluster_mass_rich.h>
#include <numcosmo/nc/data/nc_data_planck_lkl.h>
#include <numcosmo/nc/data/nc_data_xcor.h>

/* Cross-correlations */
#include <numcosmo/nc/xcor/nc_xcor.h>
#include <numcosmo/nc/xcor/nc_xcor_AB.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel_component.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel_gal.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel_cluster.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel_cluster_tophat.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel_CMB_lensing.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel_weak_lensing.h>
#include <numcosmo/nc/xcor/nc_xcor_kernel_tSZ.h>

#endif /* _NUMCOSMO_H */

