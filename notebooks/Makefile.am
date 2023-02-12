## Process this file with automake to produce Makefile.in

SUBDIRS = tutorials .

notebooksdir = \
	$(pkgdatadir)-$(VERSION)/notebooks

notebooks_DATA = \
	LRTvsFisher.ipynb                                             \
	primordial_perturbations/MagDustBounce.ipynb                  \
	primordial_perturbations/BounceVecPert.ipynb                  \
	extern_lib_comparisons/NC_CCL_Colossus_density_profiles.ipynb \
	FisherMatrixExample.ipynb                                     \
	stats_dist_tests/stats_dist_gauss.ipynb                       \
	stats_dist_tests/stats_dist_rosenbrock.ipynb

EXTRA_DIST = \
	$(notebooks_DATA)

