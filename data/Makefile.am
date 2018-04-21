## Process this file with automake to produce Makefile.in

SUBDIRS =

data_observational_DATA = \
	BBN_spline2d.obj                        \
	nc_data_bao_a_eisenstein2005.obj        \
	nc_data_bao_dhr_dar_sdss_dr11_2015.obj  \
	nc_data_bao_dhr_dar_sdss_dr11_2015_LyaF_auto_cross.obj \
	nc_data_bao_dmr_hr_sdss_dr12_2016.obj   \
	nc_data_bao_dv_eisenstein2005.obj       \
	nc_data_bao_dvdv_percival2007.obj       \
	nc_data_bao_dvdv_percival2010.obj       \
	nc_data_bao_empirical_fit_ross2015.obj  \
	nc_data_bao_empirical_fit_2d_bautista2017.obj \
	nc_data_bao_rdv_percival2007.obj        \
	nc_data_bao_rdv_percival2010.obj        \
	nc_data_bao_rdv_beutler2011.obj         \
	nc_data_bao_rdv_padmanabhan2012.obj     \
	nc_data_bao_rdv_blake2012.obj           \
	nc_data_bao_rdv_anderson2012.obj        \
	nc_data_bao_rdv_kazin2014.obj           \
	nc_data_bao_rdv_boss_qso_ata2017.obj    \
	nc_data_cmb_wmap3_shift_param.obj       \
	nc_data_cmb_wmap5_shift_param.obj       \
	nc_data_cmb_wmap7_shift_param.obj       \
	nc_data_cmb_wmap5_dist_priors.obj       \
	nc_data_cmb_wmap7_dist_priors.obj       \
	nc_data_cmb_wmap9_dist_priors.obj       \
	nc_data_hubble_busca2013_bao_wmap.obj   \
	nc_data_hubble_cabre.obj                \
	nc_data_hubble_moresco2012_bc03.obj     \
	nc_data_hubble_moresco2012_mastro.obj   \
	nc_data_hubble_moresco2015.obj          \
	nc_data_hubble_moresco2016_DR9_bc03.obj \
	nc_data_hubble_moresco2016_DR9_mastro.obj \
	nc_data_hubble_riess2008_hst.obj        \
	nc_data_hubble_stern2009.obj            \
	nc_data_hubble_simon2005.obj            \
	nc_data_hubble_zhang2012.obj            \
	nc_data_hubble_riess2016.obj            \
	nc_data_snia_diag_cfa3.obj              \
	nc_data_snia_diag_essence.obj           \
	nc_data_snia_diag_gold_157.obj          \
	nc_data_snia_diag_gold_182.obj          \
	nc_data_snia_diag_gold_206.obj          \
	nc_data_snia_diag_legacy.obj            \
	nc_data_snia_diag_sdss_emille.obj       \
	nc_data_snia_diag_union.obj             \
	nc_data_snia_diag_union2.obj            \
	nc_data_snia_diag_union2_1.obj

data_observationaldir = \
	$(pkgdatadir)-$(VERSION)/data

class_data_bbn_DATA = \
	class_data/bbn/sBBN.dat

class_data_bbndir = \
	$(pkgdatadir)-$(VERSION)/data/class_data/bbn

class_data_hyrec_DATA = \
	class_data/hyrec/Alpha_inf.dat          \
	class_data/hyrec/input.dat              \
	class_data/hyrec/R_inf.dat              \
	class_data/hyrec/two_photon_tables.dat

class_data_hyrecdir = \
	$(pkgdatadir)-$(VERSION)/data/class_data/hyrec

EXTRA_DIST = \
	$(data_observational_DATA) \
	$(class_data_bbn_DATA)     \
	$(class_data_hyrec_DATA)

