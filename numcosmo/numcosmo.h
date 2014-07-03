/***************************************************************************
 *            numcosmo.h
 *
 *  Sun May  6 17:20:29 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/nc_hicosmo_priors.h>
#include <numcosmo/nc_snia_dist_cov.h>
#include <numcosmo/scalefactor.h>
/* Cosmic thermodynamics */
#include <numcosmo/nc_recomb.h>
#include <numcosmo/nc_recomb_seager.h>

/* Perturbations */
#include <numcosmo/perturbations/linear.h>
#include <numcosmo/perturbations/covariance.h>
#include <numcosmo/perturbations/nc_hipert.h>
#include <numcosmo/perturbations/nc_hipert_iadiab.h>
#include <numcosmo/perturbations/nc_hipert_adiab.h>
#include <numcosmo/perturbations/nc_hipert_two_fluids.h>

/* Model implementations */
#include <numcosmo/model/nc_hicosmo_de.h>
#include <numcosmo/model/nc_hicosmo_de_linder.h>
#include <numcosmo/model/nc_hicosmo_de_pad.h>
#include <numcosmo/model/nc_hicosmo_de_qe.h>
#include <numcosmo/model/nc_hicosmo_de_xcdm.h>
#include <numcosmo/model/nc_hicosmo_lcdm.h>
#include <numcosmo/model/nc_hicosmo_qconst.h>
#include <numcosmo/model/nc_hicosmo_qlinear.h>
#include <numcosmo/model/nc_hicosmo_qpw.h>
#include <numcosmo/model/nc_hicosmo_qspline.h>
#include <numcosmo/model/nc_hicosmo_qgrw.h>
#include <numcosmo/model/quantum_gravity.h>

/* Large Scale Structure / Structure Formation */
#include <numcosmo/lss/nc_window.h>
#include <numcosmo/lss/nc_window_tophat.h>
#include <numcosmo/lss/nc_window_gaussian.h>
#include <numcosmo/lss/nc_transfer_func.h>
#include <numcosmo/lss/nc_transfer_func_bbks.h>
#include <numcosmo/lss/nc_transfer_func_eh.h>
#include <numcosmo/lss/nc_transfer_func_camb.h>
#include <numcosmo/lss/nc_transfer_func_pert.h>
#include <numcosmo/lss/nc_growth_func.h>
#include <numcosmo/lss/nc_matter_var.h>
#include <numcosmo/lss/nc_galaxy_acf.h>
#include <numcosmo/lss/nc_multiplicity_func.h>
#include <numcosmo/lss/nc_multiplicity_func_ps.h>
#include <numcosmo/lss/nc_multiplicity_func_st.h>
#include <numcosmo/lss/nc_multiplicity_func_jenkins.h>
#include <numcosmo/lss/nc_multiplicity_func_warren.h>
#include <numcosmo/lss/nc_multiplicity_func_tinker.h>
#include <numcosmo/lss/nc_multiplicity_func_tinker_mean.h>
#include <numcosmo/lss/nc_multiplicity_func_tinker_crit.h>
#include <numcosmo/lss/nc_mass_function.h>
#include <numcosmo/lss/nc_halo_bias_type.h>
#include <numcosmo/lss/nc_halo_bias_type_ps.h>
#include <numcosmo/lss/nc_halo_bias_type_st_spher.h>
#include <numcosmo/lss/nc_halo_bias_type_st_ellip.h>
#include <numcosmo/lss/nc_halo_bias_type_tinker.h>
#include <numcosmo/lss/nc_halo_bias_func.h>
#include <numcosmo/lss/nc_cluster_redshift.h>
#include <numcosmo/lss/nc_cluster_redshift_nodist.h>
#include <numcosmo/lss/nc_cluster_photoz_gauss.h>
#include <numcosmo/lss/nc_cluster_photoz_gauss_global.h>
#include <numcosmo/lss/nc_cluster_mass.h>
#include <numcosmo/lss/nc_cluster_mass_nodist.h>
#include <numcosmo/lss/nc_cluster_mass_lnnormal.h>
#include <numcosmo/lss/nc_cluster_mass_vanderlinde.h>
#include <numcosmo/lss/nc_cluster_mass_benson.h>
#include <numcosmo/lss/nc_cluster_mass_benson_xray.h>
#include <numcosmo/lss/nc_cluster_abundance.h>

/* Observable data */
#include <numcosmo/nc_data_snia.h>
#include <numcosmo/nc_data_dist_mu.h>
#include <numcosmo/nc_data_snia_cov.h>
#include <numcosmo/nc_data_hubble.h>
#include <numcosmo/nc_data_hubble_bao.h>
#include <numcosmo/nc_data_bao_a.h>
#include <numcosmo/nc_data_bao_dv.h>
#include <numcosmo/nc_data_bao_rdv.h>
#include <numcosmo/nc_data_bao_dvdv.h>
#include <numcosmo/nc_data_bao.h>
#include <numcosmo/nc_data_cmb_dist_priors.h>
#include <numcosmo/nc_data_cmb_shift_param.h>
#include <numcosmo/nc_data_cmb.h>
#include <numcosmo/nc_data_cluster_ncount.h>
#include <numcosmo/nc_data_cluster_poisson.h>

#endif /* _NUMCOSMO_H */
