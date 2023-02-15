## Process this file with automake to produce Makefile.in

SUBDIRS = .

nbs_datadir = $(pkgdatadir)-$(VERSION)/notebooks

apes_testsdir = $(nbs_datadir)/apes_tests
apes_tests_DATA = \
	apes_tests/funnel_mcmc.ipynb     \
	apes_tests/gaussmix2d_mcmc.ipynb \
	apes_tests/rosenbrock_mcmc.ipynb


extern_lib_comparisonsdir = $(nbs_datadir)/extern_lib_comparisons
extern_lib_comparisons_DATA = \
	extern_lib_comparisons/ccl_background.ipynb                \
	extern_lib_comparisons/ccl_bocquet.ipynb                   \
	extern_lib_comparisons/ccl_colossus_density_profiles.ipynb \
	extern_lib_comparisons/ccl_mass_functions.ipynb            \
	extern_lib_comparisons/ccl_watson.ipynb                    \
	extern_lib_comparisons/clmm_backend_comp.ipynb            


likelihood_cluster_ncountdir = $(nbs_datadir)/likelihood_cluster_ncount
likelihood_cluster_ncount_DATA = \
	likelihood_cluster_ncount/DC2_Datancount.ipynb                  \
	likelihood_cluster_ncount/NumCosmo_abundance_skysim.ipynb       \
	likelihood_cluster_ncount/cmp_cluster_ccl_numcosmo.ipynb        \
	likelihood_cluster_ncount/cmp_cluster_ccl_numcosmo_proxyZ.ipynb \
	likelihood_cluster_ncount/dc2_fof.ipynb                         \
	likelihood_cluster_ncount/generate_corner.ipynb                 \
	likelihood_cluster_ncount/generate_mcmc_data.ipynb              \
	likelihood_cluster_ncount/mean_bias_ascaso.ipynb                \
	likelihood_cluster_ncount/run_ascaso_mc_binned.ipynb            \
	likelihood_cluster_ncount/run_ascaso_mc_unbinned.ipynb          \
	likelihood_cluster_ncount/run_ascaso_mcmc_binned.ipynb          \
	likelihood_cluster_ncount/run_ascaso_mcmc_unbinned.ipynb       


likelihood_cluster_ncount_cases_comparisondir = $(nbs_datadir)/likelihood_cluster_ncount/cases_comparison
likelihood_cluster_ncount_cases_comparison_DATA = \
	likelihood_cluster_ncount/cases_comparison/comparisons.ipynb        \
	likelihood_cluster_ncount/cases_comparison/corner_comparisons.ipynb


primordial_perturbationsdir = $(nbs_datadir)/primordial_perturbations
primordial_perturbations_DATA = \
	primordial_perturbations/bounce_vector.ipynb          \
	primordial_perturbations/magnetic_dust_bounce.ipynb   \
	primordial_perturbations/vacuum_study.ipynb           \
	primordial_perturbations/vacuum_study_adiabatic.ipynb \
	primordial_perturbations/vexp_bounce.ipynb           


sniadir = $(nbs_datadir)/snia
snia_DATA = \
	snia/pantheon.ipynb   \
	snia/simulation.ipynb


stats_dist_testsdir = $(nbs_datadir)/stats_dist_tests
stats_dist_tests_DATA = \
	stats_dist_tests/stats_dist_gauss.ipynb      \
	stats_dist_tests/stats_dist_rosenbrock.ipynb


stats_toolsdir = $(nbs_datadir)/stats_tools
stats_tools_DATA = \
	stats_tools/fisher_matrix.ipynb         \
	stats_tools/likelihood_ratio_test.ipynb


tutorialsdir = $(nbs_datadir)/tutorials
tutorials_DATA = \
	tutorials/cosmo_distances.ipynb     \
	tutorials/dark_energy_intro.ipynb   \
	tutorials/primordial_universe.ipynb


unstable_testsdir = $(nbs_datadir)/unstable_tests
unstable_tests_DATA = \
	unstable_tests/non_limber.ipynb


EXTRA_DIST = \
	$(apes_tests_DATA)                                 \
	$(extern_lib_comparisons_DATA)                     \
	$(likelihood_cluster_ncount_DATA)                  \
	$(likelihood_cluster_ncount_cases_comparison_DATA) \
	$(primordial_perturbations_DATA)                   \
	$(snia_DATA)                                       \
	$(stats_dist_tests_DATA)                           \
	$(stats_tools_DATA)                                \
	$(tutorials_DATA)                                  \
	$(unstable_tests_DATA)                            
