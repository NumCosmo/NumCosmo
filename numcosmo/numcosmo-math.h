/***************************************************************************
 *            numcosmo-math.h
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

#ifndef _NUMCOSMO_MATH_H
#define _NUMCOSMO_MATH_H

/* Supported libraries and build configuration */
#include <numcosmo/build_cfg.h>

/* Macros, Constants and Enums */
#include <numcosmo/ncm_enum_types.h>
#include <numcosmo/ncm/core/ncm_c.h>
#include <numcosmo/ncm/data/ncm_catalog.h>
#include <numcosmo/ncm/sphere/ncm_sky_footprint.h>
#include <numcosmo/ncm/sphere/ncm_sky_footprint_rectangular.h>

/* MPI Objects */
#include <numcosmo/ncm/mpi/ncm_mpi_job.h>
#include <numcosmo/ncm/mpi/ncm_mpi_job_test.h>
#include <numcosmo/ncm/mpi/ncm_mpi_job_fit.h>
#include <numcosmo/ncm/mpi/ncm_mpi_job_mcmc.h>
#include <numcosmo/ncm/mpi/ncm_mpi_job_feval.h>

/* Base types and components */
#include <numcosmo/ncm/core/ncm_dtuple.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>
#include <numcosmo/ncm/algebra/ncm_nnls.h>
#include <numcosmo/ncm/algebra/ncm_poly_roots.h>
#include <numcosmo/ncm/core/ncm_serialize.h>
#include <numcosmo/ncm/core/ncm_obj_array.h>
#include <numcosmo/ncm/integration/ncm_integral1d.h>
#include <numcosmo/ncm/integration/ncm_integral1d_ptr.h>
#include <numcosmo/ncm/integration/ncm_integral_nd.h>
#include <numcosmo/ncm/core/ncm_rng.h>
#include <numcosmo/ncm/stats/ncm_stats_vec.h>
#include <numcosmo/ncm/stats/ncm_stats_dist1d.h>
#include <numcosmo/ncm/stats/ncm_stats_dist1d_spline.h>
#include <numcosmo/ncm/stats/ncm_stats_dist1d_epdf.h>
#include <numcosmo/ncm/stats/ncm_stats_dist2d.h>
#include <numcosmo/ncm/stats/ncm_stats_dist2d_spline.h>
#include <numcosmo/ncm/stats/ncm_stats_dist.h>
#include <numcosmo/ncm/stats/ncm_stats_dist_kde.h>
#include <numcosmo/ncm/stats/ncm_stats_dist_vkde.h>
#include <numcosmo/ncm/stats/ncm_stats_dist_kernel.h>
#include <numcosmo/ncm/stats/ncm_stats_dist_kernel_st.h>
#include <numcosmo/ncm/stats/ncm_stats_dist_kernel_gauss.h>
#include <numcosmo/ncm/stats/ncm_bootstrap.h>
#include <numcosmo/ncm/algebra/ncm_lapack.h>
#include <numcosmo/ncm/spline/ncm_spline.h>
#include <numcosmo/ncm/spline/ncm_spline_func.h>
#include <numcosmo/ncm/spline/ncm_spline_func_test.h>
#include <numcosmo/ncm/spline/ncm_spline_gsl.h>
#include <numcosmo/ncm/spline/ncm_spline_cubic.h>
#include <numcosmo/ncm/spline/ncm_spline_cubic_notaknot.h>
#include <numcosmo/ncm/spline/ncm_spline_cubic_d2.h>
#include <numcosmo/ncm/spline/ncm_spline_rbf.h>
#include <numcosmo/ncm/spline/ncm_spline_vec.h>
#include <numcosmo/ncm/stats/ncm_function_sample_set.h>
#include <numcosmo/ncm/spline/ncm_spline2d.h>
#include <numcosmo/ncm/spline/ncm_spline2d_spline.h>
#include <numcosmo/ncm/spline/ncm_spline2d_gsl.h>
#include <numcosmo/ncm/spline/ncm_spline2d_bicubic.h>
#include <numcosmo/ncm/powspec/ncm_powspec_corr3d.h>
#include <numcosmo/ncm/powspec/ncm_powspec_filter.h>
#include <numcosmo/ncm/powspec/ncm_powspec_sphere_proj.h>
#include <numcosmo/ncm/powspec/ncm_powspec_spline2d.h>
#include <numcosmo/ncm/powspec/ncm_powspec.h>
#include <numcosmo/ncm/dynamics/ncm_csq1d.h>
#include <numcosmo/ncm/core/ncm_func_eval.h>
#include <numcosmo/ncm/specfunc/ncm_mpsf_trig_int.h>
#include <numcosmo/ncm/specfunc/ncm_mpsf_sbessel.h>
#include <numcosmo/ncm/specfunc/ncm_sf_sbessel.h>
#include <numcosmo/ncm/specfunc/ncm_sbessel_integrator.h>
#include <numcosmo/ncm/specfunc/ncm_sbessel_integrator_fftl.h>
#include <numcosmo/ncm/specfunc/ncm_sbessel_integrator_gl.h>
#include <numcosmo/ncm/specfunc/ncm_sbessel_integrator_levin.h>
#include <numcosmo/ncm/specfunc/ncm_sbessel_ode_solver.h>
#include <numcosmo/ncm/specfunc/ncm_sf_spherical_harmonics.h>
#include <numcosmo/ncm/specfunc/ncm_mpsf_0F1.h>
#include <numcosmo/ncm/fftlog/ncm_fftlog.h>
#include <numcosmo/ncm/fftlog/ncm_fftlog_sbessel_j.h>
#include <numcosmo/ncm/fftlog/ncm_fftlog_sbessel_jljm.h>
#include <numcosmo/ncm/fftlog/ncm_fftlog_tophatwin2.h>
#include <numcosmo/ncm/fftlog/ncm_fftlog_gausswin2.h>
#include <numcosmo/ncm/sphere/ncm_spectral.h>
#include <numcosmo/ncm/model/ncm_sparam.h>
#include <numcosmo/ncm/model/ncm_vparam.h>
#include <numcosmo/ncm/model/ncm_reparam.h>
#include <numcosmo/ncm/model/ncm_model.h>
#include <numcosmo/ncm/model/ncm_model_ctrl.h>
#include <numcosmo/ncm/model/ncm_model_builder.h>
#include <numcosmo/ncm/model/ncm_model_mvnd.h>
#include <numcosmo/ncm/model/ncm_model_rosenbrock.h>
#include <numcosmo/ncm/model/ncm_model_funnel.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/ncm/model/ncm_mset_func.h>
#include <numcosmo/ncm/model/ncm_mset_func1.h>
#include <numcosmo/ncm/model/ncm_mset_func_list.h>
#include <numcosmo/ncm/spline/ncm_ode_spline.h>
#include <numcosmo/ncm/core/ncm_pln1d.h>
#include <numcosmo/ncm/model/ncm_reparam_linear.h>
#include <numcosmo/ncm/data/ncm_data.h>
#include <numcosmo/ncm/data/ncm_data_dist1d.h>
#include <numcosmo/ncm/data/ncm_data_gauss.h>
#include <numcosmo/ncm/data/ncm_data_gauss_cov.h>
#include <numcosmo/ncm/data/ncm_data_gauss_diag.h>
#include <numcosmo/ncm/data/ncm_data_poisson.h>
#include <numcosmo/ncm/data/ncm_data_gauss_cov_mvnd.h>
#include <numcosmo/ncm/data/ncm_data_rosenbrock.h>
#include <numcosmo/ncm/data/ncm_data_gaussmix2d.h>
#include <numcosmo/ncm/data/ncm_data_funnel.h>
#include <numcosmo/ncm/data/ncm_dataset.h>
#include <numcosmo/ncm/fit/ncm_likelihood.h>
#include <numcosmo/ncm/fit/ncm_prior.h>
#include <numcosmo/ncm/fit/ncm_prior_gauss.h>
#include <numcosmo/ncm/fit/ncm_prior_gauss_param.h>
#include <numcosmo/ncm/fit/ncm_prior_gauss_func.h>
#include <numcosmo/ncm/fit/ncm_prior_flat.h>
#include <numcosmo/ncm/fit/ncm_prior_flat_param.h>
#include <numcosmo/ncm/fit/ncm_prior_flat_func.h>
#include <numcosmo/ncm/core/ncm_function_cache.h>
#include <numcosmo/ncm/core/ncm_cfg.h>
#include <numcosmo/ncm/core/ncm_util.h>
#include <numcosmo/ncm/core/ncm_iset.h>
#include <numcosmo/ncm/integration/ncm_diff.h>
#include <numcosmo/ncm/core/ncm_timer.h>

/* Likelihood object */
#include <numcosmo/ncm/fit/ncm_fit.h>
#ifdef NUMCOSMO_HAVE_NLOPT
#include <numcosmo/ncm/fit/ncm_fit_nlopt.h>
#include <numcosmo/ncm_fit_nlopt_enum.h>
#endif /* NUMCOSMO_HAVE_NLOPT */
#include <numcosmo/ncm/fit/ncm_fit_levmar.h>
#include <numcosmo/ncm/fit/ncm_fit_gsl_ls.h>
#include <numcosmo/ncm/fit/ncm_fit_gsl_mm.h>
#include <numcosmo/ncm/fit/ncm_fit_gsl_mms.h>
#include <numcosmo/ncm/fit/ncm_mset_catalog.h>
#include <numcosmo/ncm/fit/ncm_mset_trans_kern.h>
#include <numcosmo/ncm/fit/ncm_mset_trans_kern_flat.h>
#include <numcosmo/ncm/fit/ncm_mset_trans_kern_gauss.h>
#include <numcosmo/ncm/fit/ncm_mset_trans_kern_cat.h>
#include <numcosmo/ncm/fit/ncm_fit_mc.h>
#include <numcosmo/ncm/fit/ncm_fit_mcbs.h>
#include <numcosmo/ncm/fit/ncm_fit_mcmc.h>
#include <numcosmo/ncm/fit/ncm_fit_esmcmc.h>
#include <numcosmo/ncm/fit/ncm_fit_esmcmc_walker.h>
#include <numcosmo/ncm/fit/ncm_fit_esmcmc_walker_stretch.h>
#include <numcosmo/ncm/fit/ncm_fit_esmcmc_walker_walk.h>
#include <numcosmo/ncm/fit/ncm_fit_esmcmc_walker_apes.h>
#include <numcosmo/ncm/fit/ncm_lh_ratio1d.h>
#include <numcosmo/ncm/fit/ncm_lh_ratio2d.h>
#include <numcosmo/ncm/algebra/ncm_quaternion.h>

/* Utilities */
#include <numcosmo/ncm/specfunc/ncm_binsplit.h>
#include <numcosmo/ncm/core/ncm_memory_pool.h>
#include <numcosmo/ncm/integration/ncm_integrate.h>

/* Spherical maps, HEALPIX implementation */
#include <numcosmo/ncm/sphere/ncm_sphere_map.h>
#include <numcosmo/ncm/sphere/ncm_sphere_nn.h>

#endif /* _NUMCOSMO_MATH_H */

