/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos, Aaron Collier and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNDIALS configuration header file.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_CONFIG_H
#define _SUNDIALS_CONFIG_H

#include "sundials/sundials_export.h"

/* Indicates that the function will not report an error via last_err,
   a return code. In C++, it is just defined as noexcept. */
#if defined(__cplusplus)
#define SUNDIALS_NOEXCEPT noexcept
#else
#define SUNDIALS_NOEXCEPT
#endif

#ifndef SUNDIALS_DEPRECATED_MSG
#  define SUNDIALS_DEPRECATED_MSG(msg) SUNDIALS_DEPRECATED
#endif

#ifndef SUNDIALS_DEPRECATED_EXPORT_MSG
#  define SUNDIALS_DEPRECATED_EXPORT_MSG(msg) SUNDIALS_EXPORT SUNDIALS_DEPRECATED_MSG(msg)
#endif

#ifndef SUNDIALS_DEPRECATED_NO_EXPORT_MSG
#  define SUNDIALS_DEPRECATED_NO_EXPORT_MSG(msg) SUNDIALS_NO_EXPORT SUNDIALS_DEPRECATED_MSG(msg)
#endif

/* ------------------------------------------------------------------
 * Define SUNDIALS version numbers
 * -----------------------------------------------------------------*/


#define SUNDIALS_VERSION "7.3.0"
#define SUNDIALS_VERSION_MAJOR 7
#define SUNDIALS_VERSION_MINOR 3
#define SUNDIALS_VERSION_PATCH 0
#define SUNDIALS_VERSION_LABEL ""
#define SUNDIALS_GIT_VERSION ""


/* ------------------------------------------------------------------
 * SUNDIALS build information
 * -----------------------------------------------------------------*/

#define SUNDIALS_C_COMPILER_HAS_BUILTIN_EXPECT
#define SUNDIALS_C_COMPILER_HAS_ATTRIBUTE_ASSUME
/* #undef SUNDIALS_C_COMPILER_HAS_BUILTIN_ASSUME */
/* #undef SUNDIALS_C_COMPILER_HAS_ASSUME */
#define SUNDIALS_C_COMPILER_HAS_ATTRIBUTE_UNUSED

/* Define precision of SUNDIALS data type 'sunrealtype'
 * Depending on the precision level, one of the following
 * three macros will be defined:
 *     #define SUNDIALS_SINGLE_PRECISION 1
 *     #define SUNDIALS_DOUBLE_PRECISION 1
 *     #define SUNDIALS_EXTENDED_PRECISION 1
 */
#define SUNDIALS_DOUBLE_PRECISION 1

/* Define type of vector indices in SUNDIALS 'sunindextype'.
 * Depending on user choice of index type, one of the following
 * two macros will be defined:
 *     #define SUNDIALS_INT64_T 1
 *     #define SUNDIALS_INT32_T 1
 */
#define SUNDIALS_INT64_T 1

/* Define the type of vector indices in SUNDIALS 'sunindextype'.
 * The macro will be defined with a type of the appropriate size.
 */
#define SUNDIALS_INDEX_TYPE int64_t

/* Define the type used for 'suncountertype'.
 * The macro will be defined with a type of the appropriate size.
 */
#define SUNDIALS_COUNTER_TYPE long int

/* Use POSIX timers if available.
 *     #define SUNDIALS_HAVE_POSIX_TIMERS
 */
#define SUNDIALS_HAVE_POSIX_TIMERS

/* BUILD CVODE with fused kernel functionality */
/* #undef SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS */

/* BUILD SUNDIALS with monitoring functionalities */
/* #undef SUNDIALS_BUILD_WITH_MONITORING */

/* BUILD SUNDIALS with profiling functionalities */
#define SUNDIALS_BUILD_WITH_PROFILING

/* Enable error checking within SUNDIALS */
/* #undef SUNDIALS_ENABLE_ERROR_CHECKS */

/* BUILD SUNDIALS with logging functionalities */
#define SUNDIALS_LOGGING_LEVEL 2

/* Build metadata */
#define SUN_C_COMPILER "GNU"
#define SUN_C_COMPILER_VERSION "15.1.1"
#define SUN_C_COMPILER_FLAGS "-Wmissing-declarations -Wcast-qual -Wno-unknown-warning-option -Wall -Wpedantic -Wextra -Wshadow -Wwrite-strings -Wcast-align -Wdisabled-optimization -Wvla -Walloca -Wduplicated-cond -Wduplicated-branches -Wunused-macros -Wunused-local-typedefs  -Werror"

#define SUN_CXX_COMPILER ""
#define SUN_CXX_COMPILER_VERSION ""
#define SUN_CXX_COMPILER_FLAGS "-Wmissing-declarations -Wcast-qual -Wno-unknown-warning-option -Wall -Wpedantic -Wextra -Wshadow -Wwrite-strings -Wcast-align -Wdisabled-optimization -Wvla -Walloca -Wduplicated-cond -Wduplicated-branches -Wunused-macros -Wunused-local-typedefs  -Werror"

#define SUN_FORTRAN_COMPILER ""
#define SUN_FORTRAN_COMPILER_VERSION ""
#define SUN_FORTRAN_COMPILER_FLAGS ""

#define SUN_BUILD_TYPE "Release"

#define SUN_JOB_ID "20250513172327"
#define SUN_JOB_START_TIME "20250513172327"

#define SUN_TPL_LIST "OPENMP;BLAS_LAPACK;PTHREAD"
#define SUN_TPL_LIST_SIZE ""

#define SUNDIALS_SPACK_VERSION ""

/* ------------------------------------------------------------------
 * SUNDIALS TPL macros
 * -----------------------------------------------------------------*/

/* Caliper */
/* #undef SUNDIALS_CALIPER_ENABLED */

/* Adiak */
/* #undef SUNDIALS_ADIAK_ENABLED */

/* Ginkgo */
/* #undef SUNDIALS_GINKGO_ENABLED */
#define SUN_GINKGO_VERSION ""

/* HYPRE */
/* #undef SUNDIALS_HYPRE_ENABLED */
#define SUN_HYPRE_VERSION ""
#define SUN_HYPRE_VERSION_MAJOR 
#define SUN_HYPRE_VERSION_MINOR 
#define SUN_HYPRE_VERSION_PATCH 

/* KLU */
/* #undef SUNDIALS_KLU_ENABLED */
#define SUN_KLU_VERSION ""

/* KOKKOS */
/* #undef SUNDIALS_KOKKOS_ENABLED */
#define SUN_KOKKOS_VERSION ""

/* KOKKOS_KERNELS */
/* #undef SUNDIALS_KOKKOS_KERNELS_ENABLED */
#define SUN_KOKKOS_KERNELS_VERSION ""

/* LAPACK */
#define SUNDIALS_BLAS_LAPACK_ENABLED
#define SUN_LAPACK_VERSION ""

/* MAGMA */
/* #undef SUNDIALS_MAGMA_ENABLED */
#define SUN_MAGMA_VERSION ""

/* MPI */
#define SUN_MPI_C_COMPILER ""
#define SUN_MPI_C_VERSION ""

#define SUN_MPI_CXX_COMPILER ""
#define SUN_MPI_CXX_VERSION ""

#define SUN_MPI_FORTRAN_COMPILER ""
#define SUN_MPI_FORTRAN_VERSION ""

/* ONEMKL */
/* #undef SUNDIALS_ONEMKL_ENABLED */
#define SUN_ONEMKL_VERSION ""

/* OpenMP */
#define SUNDIALS_OPENMP_ENABLED
#define SUN_OPENMP_VERSION ""

/* PETSC */
/* #undef SUNDIALS_PETSC_ENABLED */
#define SUN_PETSC_VERSION ""

/* PTHREADS */
/* #undef SUNDIALS_PTHREADS_ENABLED */
#define SUN_PTHREADS_VERSION ""

/* RAJA */
/* #undef SUNDIALS_RAJA_ENABLED */
#define SUN_RAJA_VERSION ""

/* SUPERLUDIST */
/* #undef SUNDIALS_SUPERLUDIST_ENABLED */
#define SUN_SUPERLUDIST_VERSION ""

/* SUPERLUMT */
/* #undef SUNDIALS_SUPERLUMT_ENABLED */
#define SUN_SUPERLUMT_VERSION ""

/* TRILLINOS */
/* #undef SUNDIALS_TRILLINOS_ENABLED */
#define SUN_TRILLINOS_VERSION ""

/* XBRAID */
/* #undef SUNDIALS_XBRAID_ENABLED */
#define SUN_XBRAID_VERSION ""

/* RAJA backends */
/* #undef SUNDIALS_RAJA_BACKENDS_CUDA */
/* #undef SUNDIALS_RAJA_BACKENDS_HIP */
/* #undef SUNDIALS_RAJA_BACKENDS_SYCL */

/* Ginkgo backends */
/* #undef SUNDIALS_GINKGO_BACKENDS_CUDA */
/* #undef SUNDIALS_GINKGO_BACKENDS_HIP */
/* #undef SUNDIALS_GINKGO_BACKENDS_OMP */
/* #undef SUNDIALS_GINKGO_BACKENDS_REF */
/* #undef SUNDIALS_GINKGO_BACKENDS_SYCL */

/* MAGMA backends */
/* #undef SUNDIALS_MAGMA_BACKENDS_CUDA */
/* #undef SUNDIALS_MAGMA_BACKENDS_HIP */

/* Set if SUNDIALS is built with MPI support, then
 *     #define SUNDIALS_MPI_ENABLED 1
 * otherwise
 *     #define SUNDIALS_MPI_ENABLED 0
 */
#define SUNDIALS_MPI_ENABLED 0

/* oneMKL interface options */
/* #undef SUNDIALS_ONEMKL_USE_GETRF_LOOP */
/* #undef SUNDIALS_ONEMKL_USE_GETRS_LOOP */

/* SUPERLUMT threading type */
#define SUNDIALS_SUPERLUMT_THREAD_TYPE "PTHREAD"

/* Trilinos with MPI is available, then
 *    #define SUNDIALS_TRILINOS_HAVE_MPI
 */
/* #undef SUNDIALS_TRILINOS_HAVE_MPI */


/* ------------------------------------------------------------------
 * SUNDIALS language macros
 * -----------------------------------------------------------------*/

/* CUDA */
/* #undef SUNDIALS_CUDA_ENABLED */
#define SUN_CUDA_VERSION ""
#define SUN_CUDA_COMPILER ""
#define SUN_CUDA_ARCHITECTURES ""

/* HIP */
/* #undef SUNDIALS_HIP_ENABLED */
#define SUN_HIP_VERSION ""
#define SUN_AMDGPU_TARGETS ""

/* SYCL options */
/* #undef SUNDIALS_SYCL_2020_UNSUPPORTED */


/* ------------------------------------------------------------------
 * SUNDIALS modules enabled
 * -----------------------------------------------------------------*/

#define SUNDIALS_ARKODE 1
#define SUNDIALS_CVODE 1
#define SUNDIALS_NVECTOR_SERIAL 1
#define SUNDIALS_NVECTOR_MANYVECTOR 1
#define SUNDIALS_NVECTOR_OPENMP 1
#define SUNDIALS_NVECTOR_PTHREADS 1
#define SUNDIALS_SUNMATRIX_BAND 1
#define SUNDIALS_SUNMATRIX_DENSE 1
#define SUNDIALS_SUNMATRIX_SPARSE 1
#define SUNDIALS_SUNLINSOL_BAND 1
#define SUNDIALS_SUNLINSOL_DENSE 1
#define SUNDIALS_SUNLINSOL_PCG 1
#define SUNDIALS_SUNLINSOL_SPBCGS 1
#define SUNDIALS_SUNLINSOL_SPFGMR 1
#define SUNDIALS_SUNLINSOL_SPGMR 1
#define SUNDIALS_SUNLINSOL_SPTFQMR 1
#define SUNDIALS_SUNLINSOL_LAPACKBAND 1
#define SUNDIALS_SUNLINSOL_LAPACKDENSE 1
#define SUNDIALS_SUNNONLINSOL_NEWTON 1
#define SUNDIALS_SUNNONLINSOL_FIXEDPOINT 1


#endif /* _SUNDIALS_CONFIG_H */
