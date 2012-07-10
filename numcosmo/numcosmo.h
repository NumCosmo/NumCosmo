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

#include <string.h>
#include <math.h>
#include <complex.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit_nlin.h>
#include <sundials/sundials_nvector.h>

/* Supported libraries and build configuration */
#include <numcosmo/build_cfg.h>

#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */

/* Macros, Constants and Enums */
#include <numcosmo/ncm_enum_types.h>
#include <numcosmo/nc_macros.h>
#include <numcosmo/nc_constants.h>
#include <numcosmo/nc_enum_types.h>

/* Base types and components */
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_lapack.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_spline_func.h>
#include <numcosmo/math/ncm_spline_gsl.h>
#include <numcosmo/math/ncm_spline_cubic.h>
#include <numcosmo/math/ncm_spline_cubic_notaknot.h>
#include <numcosmo/math/ncm_spline2d.h>
#include <numcosmo/math/ncm_spline2d_spline.h>
#include <numcosmo/math/ncm_spline2d_gsl.h>
#include <numcosmo/math/ncm_spline2d_bicubic.h>
#include <numcosmo/math/ncm_fftlog.h>
#include <numcosmo/math/ncm_sparam.h>
#include <numcosmo/math/ncm_vparam.h>
#include <numcosmo/math/ncm_reparam.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_mset_func.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/math/ncm_reparam_linear.h>
#include <numcosmo/math/function_cache.h>
#include <numcosmo/math/ncm_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/nc_distance.h>
#include <numcosmo/scalefactor.h>

/* Dataset object  */
#include <numcosmo/data/data.h>
#include <numcosmo/data/dataset.h>
#include <numcosmo/data/data_gaussian.h>
#include <numcosmo/data/data_poisson.h>
#include <numcosmo/data/data_onevardist.h>
#include <numcosmo/data/data_distance_modulus.h>
#include <numcosmo/data/data_hubble_function.h>
#include <numcosmo/data/data_baryonic_oscillation.h>
#include <numcosmo/data/data_cosmic_microwave_background.h>

/* Utilities */
#include <numcosmo/math/util.h>

/* Likelihood object */
#include <numcosmo/likelihood/likelihood.h>
#include <numcosmo/likelihood/priors.h>
#include <numcosmo/math/ncm_fit.h>
#include <numcosmo/likelihood/confidence_region.h>
#include <numcosmo/likelihood/multimin.h>
#include <numcosmo/likelihood/multimin_simplex.h>
#include <numcosmo/likelihood/least_squares.h>
#ifdef NUMCOSMO_HAVE_LEVMAR
#include <numcosmo/likelihood/levmar.h>
#endif /* NUMCOSMO_HAVE_LEVMAR */
#ifdef NUMCOSMO_HAVE_NLOPT
#include <numcosmo/likelihood/nc_nlopt.h>
#endif /* NUMCOSMO_HAVE_NLOPT */

/* Utilities */
#include <numcosmo/math/memory_pool.h>
#include <numcosmo/math/cvode_util.h>
#include <numcosmo/math/grid_one.h>
#include <numcosmo/math/mpq_tree.h>
#include <numcosmo/math/quaternion.h>
#include <numcosmo/math/integral.h>
#include <numcosmo/math/poly.h>
#include <numcosmo/math/quadrature.h>
#include <numcosmo/math/matrix_exp.h>
#include <numcosmo/math/magnus_iserles_ode.h>
#include <numcosmo/math/mp_spherical_bessel.h>
#include <numcosmo/math/mp_spherical_bessel_integral.h>
#include <numcosmo/math/spherical_bessel.h>
#include <numcosmo/math/spherical_bessel_integral.h>
#include <numcosmo/math/mp_hypergeometric_0F1.h>
#include <numcosmo/math/trig_integral.h>
#include <numcosmo/math/binsplit.h>
#include <numcosmo/math/dividedifference.h>
#include <numcosmo/math/function_eval.h>
/* #include <numcosmo/math/cvode_util.h> */

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
#include <numcosmo/model/quantum_gravity.h>

/* Cosmic thermodynamics */
#include <numcosmo/thermodyn/recomb.h>

/* Spherical maps, HEALPIX implementation */
#include <numcosmo/sphere/map.h>
#include <numcosmo/sphere/healpix.h>

/* Pert */
#include <numcosmo/perturbations/linear.h>
#include <numcosmo/perturbations/covariance.h>
#include <numcosmo/perturbations/hydrodyn_adiabatic.h>

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
#include <numcosmo/lss/read_matrix.h>
#include <numcosmo/lss/data_cluster_abundance.h>
#include <numcosmo/lss/catalog_parser.h>
#include <numcosmo/lss/print_data.h>

#endif /* _NUMCOSMO_H */
