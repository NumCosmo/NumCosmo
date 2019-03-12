## Process this file with automake to produce Makefile.in

SUNDIALS_DIR = $(top_srcdir)/numcosmo/sundials

MY_CFLAGS = \
	 -DG_LOG_DOMAIN=\"NUMCOSMO\" \
	 -I$(SUNDIALS_DIR)/include   \
	 $(SUNDIALS_DEFINES)         \
	 $(OPENMP_CFLAGS)            \
	 $(VISIBILITY_HIDDEN_FLAG)

AM_CFLAGS = \
	 -Wall \
	 $(MY_CFLAGS)

AM_LDFLAGS = $(VISIBILITY_HIDDEN_FLAG)

SUBDIRS = .

noinst_LTLIBRARIES = libsundials.la

libsundials_la_LIBADD = \
	$(OPENMP_CFLAGS)

## Sundials modules

sundials_arkode_sources = \
	include/arkode/arkode_arkstep.h      \
	include/arkode/arkode_bandpre.h      \
	include/arkode/arkode_bbdpre.h       \
	include/arkode/arkode_butcher_dirk.h \
	include/arkode/arkode_butcher_erk.h  \
	include/arkode/arkode_butcher.h      \
	include/arkode/arkode_erkstep.h      \
	include/arkode/arkode.h              \
	include/arkode/arkode_ls.h           \
	include/arkode/arkode_mristep.h      \
	src/arkode/arkode_adapt.c            \
	src/arkode/arkode_arkstep_impl.h     \
	src/arkode/arkode_bandpre.c          \
	src/arkode/arkode_bbdpre_impl.h      \
	src/arkode/arkode_butcher_erk.c      \
	src/arkode/arkode_erkstep_impl.h     \
	src/arkode/arkode_interp.c           \
	src/arkode/arkode_ls.c               \
	src/arkode/arkode_mristep_impl.h     \
	src/arkode/arkode_root_impl.h        \
	src/arkode/arkode_adapt_impl.h       \
	src/arkode/arkode_arkstep_io.c       \
	src/arkode/arkode_bandpre_impl.h     \
	src/arkode/arkode_butcher.c          \
	src/arkode/arkode.c                  \
	src/arkode/arkode_erkstep_io.c       \
	src/arkode/arkode_interp_impl.h      \
	src/arkode/arkode_ls_impl.h          \
	src/arkode/arkode_mristep_io.c       \
	src/arkode/arkode_arkstep.c          \
	src/arkode/arkode_arkstep_nls.c      \
	src/arkode/arkode_bbdpre.c           \
	src/arkode/arkode_butcher_dirk.c     \
	src/arkode/arkode_erkstep.c          \
	src/arkode/arkode_impl.h             \
	src/arkode/arkode_io.c               \
	src/arkode/arkode_mristep.c          \
	src/arkode/arkode_root.c

sundials_cvode_sources = \
	include/cvode/cvode_bandpre.h  \
	include/cvode/cvode_bbdpre.h   \
	include/cvode/cvode_diag.h     \
	include/cvode/cvode.h          \
	include/cvode/cvode_ls.h       \
	src/cvode/cvode_bandpre.c      \
	src/cvode/cvode_bandpre_impl.h \
	src/cvode/cvode_bbdpre.c       \
	src/cvode/cvode_bbdpre_impl.h  \
	src/cvode/cvode.c              \
	src/cvode/cvode_diag.c         \
	src/cvode/cvode_diag_impl.h    \
	src/cvode/cvode_impl.h         \
	src/cvode/cvode_io.c           \
	src/cvode/cvode_ls.c           \
	src/cvode/cvode_ls_impl.h      \
	src/cvode/cvode_nls.c

sundials_nvector_sources = \
	include/nvector/nvector_openmp.h        \
	include/nvector/nvector_pthreads.h      \
	include/nvector/nvector_serial.h        \
	src/nvector/openmp/nvector_openmp.c     \
	src/nvector/pthreads/nvector_pthreads.c \
	src/nvector/serial/nvector_serial.c

sundials_sunlinsol_sources = \
	include/sunlinsol/sunlinsol_band.h                \
	include/sunlinsol/sunlinsol_dense.h               \
	include/sunlinsol/sunlinsol_lapackband.h          \
	include/sunlinsol/sunlinsol_lapackdense.h         \
	include/sunlinsol/sunlinsol_pcg.h                 \
	include/sunlinsol/sunlinsol_spbcgs.h              \
	include/sunlinsol/sunlinsol_spfgmr.h              \
	include/sunlinsol/sunlinsol_spgmr.h               \
	include/sunlinsol/sunlinsol_sptfqmr.h             \
	src/sunlinsol/band/sunlinsol_band.c               \
	src/sunlinsol/dense/sunlinsol_dense.c             \
	src/sunlinsol/lapackband/sunlinsol_lapackband.c   \
	src/sunlinsol/lapackdense/sunlinsol_lapackdense.c \
	src/sunlinsol/pcg/sunlinsol_pcg.c                 \
	src/sunlinsol/spbcgs/sunlinsol_spbcgs.c           \
	src/sunlinsol/spfgmr/sunlinsol_spfgmr.c           \
	src/sunlinsol/spgmr/sunlinsol_spgmr.c             \
	src/sunlinsol/sptfqmr/sunlinsol_sptfqmr.c

sundials_sunmatrix_sources = \
	include/sunmatrix/sunmatrix_band.h      \
	include/sunmatrix/sunmatrix_dense.h     \
	include/sunmatrix/sunmatrix_sparse.h    \
	src/sunmatrix/band/sunmatrix_band.c     \
	src/sunmatrix/dense/sunmatrix_dense.c   \
	src/sunmatrix/sparse/sunmatrix_sparse.c

sundials_sunnonlinsol_sources = \
	include/sunnonlinsol/sunnonlinsol_fixedpoint.h        \
	include/sunnonlinsol/sunnonlinsol_newton.h            \
	src/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c \
	src/sunnonlinsol/newton/sunnonlinsol_newton.c

sundials_sources = \
	include/sundials/sundials_band.h                \
	include/sundials/sundials_lapack.h              \
	include/sundials/sundials_nvector.h             \
	include/sundials/sundials_spfgmr.h              \
	include/sundials/sundials_version.h             \
	include/sundials/sundials_config.h              \
	include/sundials/sundials_linearsolver.h        \
	include/sundials/sundials_nvector_senswrapper.h \
	include/sundials/sundials_spgmr.h               \
	include/sundials/sundials_dense.h               \
	include/sundials/sundials_math.h                \
	include/sundials/sundials_pcg.h                 \
	include/sundials/sundials_sptfqmr.h             \
	include/sundials/sundials_direct.h              \
	include/sundials/sundials_matrix.h              \
	include/sundials/sundials_sparse.h              \
	include/sundials/sundials_iterative.h           \
	include/sundials/sundials_nonlinearsolver.h     \
	include/sundials/sundials_spbcgs.h              \
	include/sundials/sundials_types.h               \
	src/sundials/sundials_band.c                    \
	src/sundials/sundials_iterative.c               \
	src/sundials/sundials_matrix.c                  \
	src/sundials/sundials_nvector.c                 \
	src/sundials/sundials_sparse.c                  \
	src/sundials/sundials_spgmr.c                   \
	src/sundials/sundials_dense.c                   \
	src/sundials/sundials_linearsolver.c            \
	src/sundials/sundials_nvector_senswrapper.c     \
	src/sundials/sundials_spbcgs.c                  \
	src/sundials/sundials_sptfqmr.c                 \
	src/sundials/sundials_direct.c                  \
	src/sundials/sundials_math.c                    \
	src/sundials/sundials_nonlinearsolver.c         \
	src/sundials/sundials_pcg.c                     \
	src/sundials/sundials_spfgmr.c                  \
	src/sundials/sundials_version.c

libsundials_la_SOURCES = \
	$(sundials_arkode_sources)       \
	$(sundials_cvode_sources)        \
	$(sundials_nvector_sources)      \
	$(sundials_sunlinsol_sources)    \
	$(sundials_sunmatrix_sources)    \
	$(sundials_sunnonlinsol_sources) \
	$(sundials_sources)

EXTRA_DIST = \
	LICENSE            \
	NOTICE             \
	src/cvode/LICENSE  \
	src/cvode/README   \
	src/cvode/NOTICE   \
	src/arkode/LICENSE \
	src/arkode/README  \
	src/arkode/NOTICE

CLEANFILES =
