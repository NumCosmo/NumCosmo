## Process this file with automake to produce Makefile.in

PLC_DIR = $(top_srcdir)/numcosmo/plc

MY_CFLAGS = \
	 -DG_LOG_DOMAIN=\"NUMCOSMO\" \
	 $(GLIB_CFLAGS)       \
	 $(SUNDIALS_CFLAGS)   \
	 $(GSL_CFLAGS)        \
	 $(FFTW3_CFLAGS)      \
	 $(CFITSIO_CFLAGS)    \
	 $(NLOPT_CFLAGS)      \
	 -I$(top_srcdir)      \
	 -I$(PLC_DIR)         \
	 -I$(PLC_DIR)/cldf    \
	 -I$(PLC_DIR)/minipmc \
	 -I$(PLC_DIR)/plik    \
	 $(PLC_DEFINES)       \
	 $(PLC_NOWARN_FLAGS)  \
	 $(OPENMP_CFLAGS)     \
	 $(VISIBILITY_HIDDEN_FLAG)

AM_CFLAGS = \
	 -Wall \
	 $(MY_CFLAGS)

AM_FCFLAGS = \
	$(PLC_DEFINES) \
	$(OPENMP_FCFLAGS)

AM_LDFLAGS = $(VISIBILITY_HIDDEN_FLAG)

SUBDIRS = .

noinst_LTLIBRARIES = libclik.la

libclik_la_LIBADD = \
	$(OPENMP_CFLAGS)

## Fortran modules
actspt_extra.mod:             actspt/clik_actspt_f90.lo
act_equa_likelihood.mod:      actspt/ACT_equa_likelihood.lo
act_south_likelihood.mod:     actspt/ACT_south_likelihood.lo
bflike.mod:                   bflike/bflike.lo
bflike_extra.mod:             bflike/clik_bflike_f90.lo
bflike_qu.mod:                bflike/bflike_QUonly.lo
bflike_qu_extra.mod:          bflike/clik_bflike_f90.lo
bflike_smw.mod:               bflike/bflike_smw.lo
bflike_smw_extra.mod:         bflike/clik_bflike_f90.lo
bflike_utils_smw.mod:         bflike/bflike_smw.lo
br_mod_dist.mod:              lowlike/br_mod_dist.lo
camspec_extra.mod:            camspec/clik_CAMspec_f90.lo
clik.mod:                     clik.lo
comm_br_mod.mod:              gibbs/comm_br_mod.lo
comm_gauss_br_mod.mod:        gibbs/comm_gauss_br_mod.lo
comm_gauss_br_mod_v3.mod:     gibbs/comm_gauss_br_mod_v3.lo
comm_lowl_extra.mod:          gibbs/clik_gibbs_f90.lo
comm_lowl_mod_dist.mod:       gibbs/comm_lowl_mod_dist.lo
fiducial_cls_mod_smw.mod:     bflike/bflike_smw.lo
fitstools_smw.mod:            bflike/fitstools_smw.lo
foregrounds_loading.mod:      actspt/Foregrounds_loading.lo
gibbs_extra.mod:              gibbs/clik_gibbs_f90.lo
healpix_types.mod:            lowlike/healpix_types.lo
highell_likelihood.mod:       actspt/Highell_likelihood.lo
highell_options.mod:          actspt/Highell_options.lo
highell_subroutines.mod:      actspt/Highell_subroutines.lo
long_intrinsic_smw.mod:       bflike/long_intrinsic_smw.lo
lowlike_extra.mod:            lowlike/clik_lowlike_f90.lo
planck_likelihood.mod:        lowlike/Planck_likelihood.lo
planck_options.mod:           lowlike/Planck_options.lo
planck_teeebb_lowl.mod:       lowlike/Planck_teeebb_pixlike.lo
plik_cmbonly.mod:             cmbonly/plik_cmbonly.lo
plik_cmbonly_extra.mod:       cmbonly/clik_cmbonly_f90.lo
spt_keisler_likelihood.mod:   actspt/SPT_keisler_likelihood.lo
spt_reichardt_likelihood.mod: actspt/SPT_reichardt_likelihood.lo
temp_like.mod:                camspec/CAMspec.lo
temp_like_camspec.mod:        camspec/temp_like.lo
temp_like_camspec3.mod:       camspec/temp_like_v3.lo

## Fortran dependencies

actspt/ACT_equa_likelihood.lo: \
	foregrounds_loading.mod \
	highell_options.mod     \
	highell_subroutines.mod

actspt/ACT_south_likelihood.lo: \
	foregrounds_loading.mod \
	highell_options.mod     \
	highell_subroutines.mod

actspt/clik_actspt_f90.lo: \
	highell_likelihood.mod \
	highell_options.mod

actspt/Foregrounds_loading.lo:  \
	highell_options.mod     \
	highell_subroutines.mod

actspt/Highell_likelihood.lo: \
	act_equa_likelihood.mod      \
	act_south_likelihood.mod     \
	foregrounds_loading.mod      \
	highell_options.mod          \
	spt_keisler_likelihood.mod   \
	spt_reichardt_likelihood.mod

actspt/Inversecov_selection.lo: \
	highell_options.mod

actspt/SPT_keisler_likelihood.lo: \
	foregrounds_loading.mod \
	highell_options.mod     \
	highell_subroutines.mod 

actspt/SPT_reichardt_likelihood.lo: \
	foregrounds_loading.mod \
	highell_options.mod     \
	highell_subroutines.mod 

bflike/bflike.lo: \
	fitstools_smw.mod \
	healpix_types.mod

bflike/bflike_QUonly.lo: \
	fitstools_smw.mod \
	healpix_types.mod

bflike/bflike_smw.lo: \
	fitstools_smw.mod \
	healpix_types.mod

bflike/clik_bflike_f90.lo: \
	bflike.mod           \
	bflike_qu.mod        \
	bflike_smw.mod       \
	bflike_utils_smw.mod

bflike/fitstools_smw.lo: \
	healpix_types.mod \
	long_intrinsic_smw.mod

bflike/long_intrinsic_smw.lo: \
	healpix_types.mod

camspec/clik_CAMspec_f90.lo:   \
	temp_like.mod          \
	temp_like_camspec.mod  \
	temp_like_camspec3.mod

cmbonly/clik_cmbonly_f90.lo: \
	plik_cmbonly.mod

gibbs/clik_gibbs_f90.lo: \
	comm_br_mod.mod          \
	comm_gauss_br_mod.mod    \
	comm_gauss_br_mod_v3.mod \
	comm_lowl_mod_dist.mod

gibbs/validate_comm_lowl_mod.lo: \
	comm_lowl_mod_dist.mod

lowlike/clik_lowlike_f90.lo: \
	healpix_types.mod     \
	planck_likelihood.mod \
	planck_options.mod

lowlike/Planck_likelihood.lo: \
	planck_options.mod     \
	planck_teeebb_lowl.mod

lowlike/Planck_teeebb_pixlike.lo: \
	planck_options.mod

plik/clik_plik_f90.lo: \
	clik.mod

## CLik modules sources

clik_tool_sources =            \
        minipmc/pmc.h          \
        minipmc/maths_base.h   \
	minipmc/errorlist.c    \
	minipmc/errorlist.h    \
	minipmc/io.c           \
	minipmc/io.h           \
	minipmc/distribution.c \
	cldf/cldf.c            \
	cldf/cldf.h            \
	clik_dic.c             \
	clik_dic.h             \
	lapack_clik.h

clik_main_sources = \
	clik.c          \
	clik.h          \
	lklbs.c         \
	lklbs.h         \
	lowly_common.c  \
	lowly_common.h  \
	clik_helper.c   \
	clik_helper.h   \
	clik_lensing.h

clik_act_spt_lkl_sources = \
	actspt/Highell_options.F90          \
	actspt/Highell_subroutines.F90      \
	actspt/Foregrounds_loading.F90      \
	actspt/ACT_equa_likelihood.F90      \
	actspt/SPT_reichardt_likelihood.F90 \
	actspt/ACT_south_likelihood.F90     \
	actspt/SPT_keisler_likelihood.F90   \
	actspt/Highell_likelihood.F90       \
	actspt/clik_actspt_f90.f90          \
	actspt/clik_actspt.c

clik_camspec_lkl_sources =   \
	camspec/clik_CAMspec_f90.F90 \
	camspec/CAMspec.f90          \
	camspec/temp_like.f90        \
	camspec/temp_like_v3.f90     \
	camspec/clik_CAMspec.c

clik_gibbs_lkl_sources = \
	gibbs/comm_br_mod.f90          \
	gibbs/comm_gauss_br_mod.f90    \
	gibbs/comm_gauss_br_mod_v3.f90 \
	gibbs/comm_lowl_mod_dist.f90   \
	gibbs/clik_gibbs_f90.f90       \
	gibbs/clik_gibbs.c

clik_lowlike_lkl_sources = \
	lowlike/healpix_types.f90         \
	lowlike/read_archive_map.f90      \
	lowlike/read_fits.f90             \
	lowlike/br_mod_dist.f90           \
	lowlike/Planck_options.F90        \
	lowlike/Planck_teeebb_pixlike.F90 \
	lowlike/Planck_likelihood.F90     \
	lowlike/clik_lowlike_f90.f90      \
	lowlike/clik_lowlike.c

clik_bflike_lkl_sources =     \
	bflike/long_intrinsic_smw.f90 \
	bflike/fitstools_smw.f90      \
	bflike/bflike_QUonly.f90      \
	bflike/bflike.f90             \
	bflike/bflike_smw.f90         \
	bflike/clik_bflike_f90.f90    \
	bflike/clik_bflike.c

clik_pliklite_lkl_sources =   \
	cmbonly/plik_cmbonly.f90      \
	cmbonly/clik_cmbonly_f90.f90  \
	cmbonly/clik_cmbonly.c

clik_plik_lkl_sources =   \
	plik/smica.c                 \
	plik/smica.h                 \
	plik/clik_hfipack.c          \
	plik/clik_parametric.c       \
	plik/clik_parametric.h       \
	plik/clik_parametric_addon.c \
	plik/clik_parametric_addon.h \
	plik/component_plugin/rel2015/fg2015.c    \
	plik/component_plugin/rel2015/corrnoise.c \
	plik/component_plugin/rel2015/leakage.c

clik_lens_lkl_sources = \
	lenslike/plenslike/plenslike.h          \
	lenslike/plenslike/plenslike_dat_mono.c \
	lenslike/plenslike/plenslike_dat_quad.c \
	lenslike/plenslike/plenslike_dat_qecl.c \
	lenslike/plenslike/plenslike_dat_full.c \
	lenslike/plenslike/qest.c               \
	lenslike/plenslike/wignerd.c            \
	clik_lensing.c

clik_simall_sources = \
	simall/clik_simall.c

libclik_la_SOURCES = \
	$(clik_tool_sources)         \
	$(clik_main_sources)         \
	$(clik_act_spt_lkl_sources)  \
	$(clik_camspec_lkl_sources)  \
	$(clik_gibbs_lkl_sources)    \
	$(clik_lowlike_lkl_sources)  \
	$(clik_bflike_lkl_sources)   \
	$(clik_pliklite_lkl_sources) \
	$(clik_plik_lkl_sources)     \
	$(clik_lens_lkl_sources)     \
	$(clik_simall_sources)

EXTRA_DIST = \
	lowlike/read_archive_map.fh \
	lowlike/read_fits.fh        


CLEANFILES = \
	actspt_extra.mod act_equa_likelihood.mod act_south_likelihood.mod  \
	bflike.mod bflike_extra.mod bflike_qu.mod bflike_qu_extra.mod      \
	bflike_smw.mod bflike_smw_extra.mod bflike_utils_smw.mod           \
	br_mod_dist.mod camspec_extra.mod clik.mod comm_br_mod.mod         \
	comm_gauss_br_mod.mod comm_gauss_br_mod_v3.mod comm_lowl_extra.mod \
	comm_lowl_mod_dist.mod fiducial_cls_mod_smw.mod fitstools_smw.mod  \
	foregrounds_loading.mod gibbs_extra.mod healpix_types.mod          \
	highell_likelihood.mod highell_options.mod highell_subroutines.mod \
	long_intrinsic_smw.mod lowlike_extra.mod planck_likelihood.mod     \
	planck_options.mod planck_teeebb_lowl.mod plik_cmbonly.mod         \
	plik_cmbonly_extra.mod spt_keisler_likelihood.mod                  \
	spt_reichardt_likelihood.mod temp_like.mod temp_like_camspec.mod   \
	temp_like_camspec3.mod
