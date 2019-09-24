/***************************************************************************
 *            numcosmo-math.h
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

#ifndef _NUMCOSMO_MATH_H
#define _NUMCOSMO_MATH_H

/* Supported libraries and build configuration */
#include <numcosmo/build_cfg.h>

/* Macros, Constants and Enums */
#include <numcosmo/ncm_enum_types.h>
#include <numcosmo/math/ncm_c.h>

/* MPI Objects */
#include <numcosmo/math/ncm_mpi_job.h>
#include <numcosmo/math/ncm_mpi_job_test.h>
#include <numcosmo/math/ncm_mpi_job_fit.h>
#include <numcosmo/math/ncm_mpi_job_mcmc.h>

/* Base types and components */
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_serialize.h>
#include <numcosmo/math/ncm_obj_array.h>
#include <numcosmo/math/ncm_integral1d.h>
#include <numcosmo/math/ncm_integral1d_ptr.h>
#include <numcosmo/math/ncm_rng.h>
#include <numcosmo/math/ncm_stats_vec.h>
#include <numcosmo/math/ncm_stats_dist1d.h>
#include <numcosmo/math/ncm_stats_dist1d_spline.h>
#include <numcosmo/math/ncm_stats_dist1d_epdf.h>
#include <numcosmo/math/ncm_stats_dist2d.h>
#include <numcosmo/math/ncm_stats_dist2d_spline.h>
#include <numcosmo/math/ncm_stats_dist_nd.h>
#include <numcosmo/math/ncm_stats_dist_nd_kde_gauss.h>
#include <numcosmo/math/ncm_bootstrap.h>
#include <numcosmo/math/ncm_lapack.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/math/ncm_spline_func.h>
#include <numcosmo/math/ncm_spline_gsl.h>
#include <numcosmo/math/ncm_spline_cubic.h>
#include <numcosmo/math/ncm_spline_cubic_notaknot.h>
#include <numcosmo/math/ncm_spline_rbf.h>
#include <numcosmo/math/ncm_spline2d.h>
#include <numcosmo/math/ncm_spline2d_spline.h>
#include <numcosmo/math/ncm_spline2d_gsl.h>
#include <numcosmo/math/ncm_spline2d_bicubic.h>
#include <numcosmo/math/ncm_powspec.h>
#include <numcosmo/math/ncm_powspec_filter.h>
#include <numcosmo/math/ncm_powspec_corr3d.h>
#include <numcosmo/math/ncm_hoaa.h>
#include <numcosmo/math/ncm_csq1d.h>
#include <numcosmo/math/ncm_func_eval.h>
#include <numcosmo/math/grid_one.h>
#include <numcosmo/math/ncm_mpsf_trig_int.h>
#include <numcosmo/math/ncm_mpsf_sbessel.h>
#include <numcosmo/math/ncm_mpsf_sbessel_int.h>
#include <numcosmo/math/ncm_sf_sbessel.h>
#include <numcosmo/math/ncm_sf_sbessel_int.h>
#include <numcosmo/math/ncm_sf_spherical_harmonics.h>
#include <numcosmo/math/ncm_mpsf_0F1.h>
#include <numcosmo/math/ncm_fftlog.h>
#include <numcosmo/math/ncm_fftlog_sbessel_j.h>
#include <numcosmo/math/ncm_fftlog_tophatwin2.h>
#include <numcosmo/math/ncm_fftlog_gausswin2.h>
#include <numcosmo/math/ncm_sparam.h>
#include <numcosmo/math/ncm_vparam.h>
#include <numcosmo/math/ncm_reparam.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/math/ncm_model_builder.h>
#include <numcosmo/math/ncm_model_mvnd.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_mset_func.h>
#include <numcosmo/math/ncm_mset_func1.h>
#include <numcosmo/math/ncm_mset_func_list.h>
#include <numcosmo/math/ncm_calc.h>
#include <numcosmo/math/ncm_ode_spline.h>
#include <numcosmo/math/ncm_reparam_linear.h>
#include <numcosmo/math/ncm_data.h>
#include <numcosmo/math/ncm_data_dist1d.h>
#include <numcosmo/math/ncm_data_gauss.h>
#include <numcosmo/math/ncm_data_gauss_cov.h>
#include <numcosmo/math/ncm_data_gauss_diag.h>
#include <numcosmo/math/ncm_data_poisson.h>
#include <numcosmo/math/ncm_data_gauss_cov_mvnd.h>
#include <numcosmo/math/ncm_dataset.h>
#include <numcosmo/math/ncm_likelihood.h>
#include <numcosmo/math/ncm_prior.h>
#include <numcosmo/math/ncm_prior_gauss.h>
#include <numcosmo/math/ncm_prior_gauss_param.h>
#include <numcosmo/math/ncm_prior_gauss_func.h>
#include <numcosmo/math/ncm_prior_flat.h>
#include <numcosmo/math/ncm_prior_flat_param.h>
#include <numcosmo/math/ncm_prior_flat_func.h>
#include <numcosmo/math/ncm_function_cache.h>
#include <numcosmo/math/ncm_cfg.h>
#include <numcosmo/math/ncm_util.h>
#include <numcosmo/math/ncm_diff.h>
#include <numcosmo/math/ncm_ode.h>
#include <numcosmo/math/ncm_ode_eval.h>
#include <numcosmo/math/ncm_timer.h>

/* Likelihood object */
#include <numcosmo/math/ncm_fit.h>
#ifdef NUMCOSMO_HAVE_NLOPT
#include <numcosmo/math/ncm_fit_nlopt.h>
#include <numcosmo/ncm_fit_nlopt_enum.h>
#endif /* NUMCOSMO_HAVE_NLOPT */
#include <numcosmo/math/ncm_fit_levmar.h>
#include <numcosmo/math/ncm_fit_gsl_ls.h>
#include <numcosmo/math/ncm_fit_gsl_mm.h>
#include <numcosmo/math/ncm_fit_gsl_mms.h>
#include <numcosmo/math/ncm_mset_catalog.h>
#include <numcosmo/math/ncm_mset_trans_kern.h>
#include <numcosmo/math/ncm_mset_trans_kern_flat.h>
#include <numcosmo/math/ncm_mset_trans_kern_gauss.h>
#include <numcosmo/math/ncm_mset_trans_kern_cat.h>
#include <numcosmo/math/ncm_fit_mc.h>
#include <numcosmo/math/ncm_fit_mcbs.h>
#include <numcosmo/math/ncm_fit_mcmc.h>
#include <numcosmo/math/ncm_fit_esmcmc.h>
#include <numcosmo/math/ncm_fit_esmcmc_walker.h>
#include <numcosmo/math/ncm_fit_esmcmc_walker_stretch.h>
#include <numcosmo/math/ncm_fit_esmcmc_walker_walk.h>
#include <numcosmo/math/ncm_fit_esmcmc_walker_aps.h>
#include <numcosmo/math/ncm_lh_ratio1d.h>
#include <numcosmo/math/ncm_lh_ratio2d.h>
#include <numcosmo/math/ncm_abc.h>
#include <numcosmo/math/ncm_quaternion.h>

/* Utilities */
#include <numcosmo/math/ncm_memory_pool.h>
#include <numcosmo/math/mpq_tree.h>
#include <numcosmo/math/integral.h>
#include <numcosmo/math/poly.h>
#include <numcosmo/math/quadrature.h>
#include <numcosmo/math/matrix_exp.h>
#include <numcosmo/math/magnus_iserles_ode.h>
#include <numcosmo/math/binsplit.h>
#include <numcosmo/math/dividedifference.h>

/* Spherical maps, HEALPIX implementation */
#include <numcosmo/math/ncm_sphere_map.h>

#endif /* _NUMCOSMO_MATH_H */
