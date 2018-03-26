CHANGELOG
----------------------

[Current]

[v0.14.0]
 * Fixed data install path. Tweaked fit tests.

 * NumCosmo version written by autoconf in numcosmo-docs.sgml.

 * Updated print in python scripts.

 * Fixed deploy file.

 * Bumped version and updated ChangeLog.

 * Tweaking tests.

 * Fixed new package (in progress).

 * Implementing ncm_data_voigt object: work in progress

 * Workaround for older sundials bug (2).

 * Tweaking tests and workaround for older sundials bug.

 * Fixed missing prototype and numpy on travis-ci.

 * Tweaking tests.

 * Minor fixes (portability related).

 * Improved error msgs in ncm_spline_func.c. Working on ncm_qm_prop.

 * Working on ncm_qm_prop.

 * Working on ncm_qm_prop.

 * Conditional compilation of NcmQMProp.

 * New test QM object.

 * Created three tests: distance, density profile (NFW) and surface mass density.

     Distance: Implemented functions to compute the comoving and transverse
     distances from z to infinity

     DensityProfileNFW: fixed bugs, tested

     NcWLSurfaceMassDensity : implemented functions like convergence, shear and
     reduced shear... all funtions tested using NFW density profile.

     Included example for NcWLSurfaceMassDensity.

 * Fixed log spacing.

 * Fixed parameter ranges.

 * Improvement in example_hiprim_Tmodes.py. Minor fixes. Added C(theta) calc in
     ncm_sphere_map.h.

 * Updated default alpha in nc_snia_dist_cov.h.

 * Fixed error in reading m2lnp_var from catalog file.

 * Added option to calculate evidence in mcat_analyze.

 * Added volume estimator testing to test_ncm_fit_esmcmc (fixed bug in vol
     estimation when the catalog contains repeated points).

 * New code for Bayesian evidence and posterior volume and its unit tests.

 * Better parameters for test_ncm_sphere_map.c.

 * Renamed test to match the new object name.

 * Encapsulated the NcmMSetCatalog object. New NcmMSetCatalog tests.

 * Removed travis-ci brew science tap.

 * Working on wl related objects (in progress).

 * New test test_ncm_fit_esmcmc.c. Fixed solar/G consistence.

 * Travis-ci backports (removed debug).

 * Travis-ci backports (debug1).

 * Travis-ci backports (debug).

 * Travis-ci backports.

 * Travis-ci backports.

 * Travis-ci backports.

 * Travis-ci backports.

 * Testing backports in travis-ci.

 * New test for NcmFit. GSL requirement changed to 2.0. Old code related to older
     gsl versions removed. New NcmDataset constructor. Updated NcmFitLS.

 * Removed old doc file.

 * Finished sphere_map. Removed old codes. Added support fo sundials 3.x.x.

 * Finished organizing and further folding of the sphere_map map2alm algo.

 * Fixing gcc install issue in macos/travis-ci.

 * Finished organizing code.

 * Fixed old sphere/ scan in docs. Split block algo from main code in
     ncm_sphere_map.c.

 * Renamed new sphere map object: NcmSphereMapPix -> NcmSphereMap.

 * Removed old Spherical Map/Healpix implementation.

 * removed clang support in travis-ci.

 * Travis...

 * Travis...

 * Still fixing travis...

 * Trying to fix openmp+clang problem in travis.

 * Adding openmp cflags to the introspection cflags.

 * Removed march from the options added zhen using --enable-opt-cflags. Removed
     debug messages from ncm_sphere_map_pix.c. Updated ax_gcc_archflag.m4.
     Testing clang options for travis-ci.

 * Finished block opt for sphere map pix.

 * Updates.

 * Improvements in NcmSphereMapPix (optimization of alm2map in progress).

 * Finished alm2map (missing block algo).

 * Finished block algorithms and tests.

 * Finishing new block interface for SF SH. Added tests for SF SH.

 * Updated m4/ax_cc_maxopt.m4 and m4/ax_gcc_archflag.m4. Working on optimization
     of ncm_sphere_map_pix.

 * Optimizing code (unstable).

 * Missing m4 in some platforms.

 * Testing new optimizations, new Spherical Harmonics object, finishing sphere_map
     object.

 * New example NcmDiff.

 * Applied fix (gtkdoc/glib-mkenums scan) to the ncm namespace.

 * Fixed gtkdoc/glib-mkenums scan problem.

 * Added a guard to avoid introspection into wrong headers. Solving the enum
     parsing erros/warnings (glib's bug).

 * Update deprecated glib function.

 * New kinetic w function.

 * Updated scripts/corner.py, new acoustic scale in Mpc function in NcDistance.

 * Documentation.

 * New option to print out functions in mcat_analyze. Fixed minor bug in
     plot2DCorner.

 * Reordered fortran probing in configure.ac (fix problems in some weird gcc
     installations).

 * Removed two prints and p-mode member. Included test of the x and y (splines)
     bounds to compute m2lnL (bao empirical fit 2D).

 * Added Bautista et al. obj in the Makefile.

 * Corrected typo on the documentaion.

 * Removed log printing.

 * Fixing instrospection in macos.

 * cat config.log in travis ci.

 * Fixing macos build.

 * Still trying to fix maxos build.

 * Implemented NcmDataDist2d: data described by two-variable (arbitrary)
     distribution.

     Implemented NcDataBaoEmpiricalFit2d: included Bautista et al. 2017 data
     (SDSS/BOSS DR12).

 * Testing travis ci.

 * Trying fix numpy install using brew.

 * Right order.

 * overwrite option added to numpy at travis ci.

 * Added numpy to brew install in travis ci.

 * Added condition on having fftw to the unit test.

 * Minor tweaks. Improved NcmFftlog and added unit testing. Updated NcmABC and
     NcABCClusterNCount.

 * Update py_sline_gauss.py

     fixing dsymm call to met the new arguments order.
 * Update py_sline_gauss.py

     Fixing dsymm parameters order
 * Finished objects infrastructure.

 * Added tests for inverse distribution computation.

 * Improved ncm_stats_dist1d, including new tests in test_ncm_stats_dist1d_epdf.

 * Comment Evrard's formula for ST multiplicity.

 * Minor tweaks.

 * Improving ncm_stats_dist1 (not ready!)

 * Fix minor bugs.

 * Minor tweaks.

 * New structure for NcDensityProfile and NcDensityProfileNFW objects. They must
     implement more functions, which are used to compute the Weak Lensing (WL)
     surface mass density, WL shear... Work in progress! Not tested.

     New object: NcWLSurfaceMassDensity. Work in progress! Not finalized nor
     tested.

     Property included in NcmStats2dSpline to indicate which marginal to be
     computed.

 * Added arxiv numbers.

 * Fixed doc.

 * Fixed freed null pointer in NcmReparam.

 * Removed debug message in darkenergy. Fixed minor leaks.

 * Implemented framework for Fisher matrix calculation. Finished implementation of
     general Fisher matrix calculation for NcmDataGauss* family. Organized
     methods of NcmFit to calculated covariance through observed or expected
     Fisher matrix. Added new methods for NcmMatrix.

 * Improved error handling in NcmDiff.

 * Removed old ncm_numdiff functions. Updated code to use the new NcmDiff object.

 * Fixed minor typos. Included NcmDiff use in NcmFit.

 * Improved documentation of NcHICosmoVexp. Added paper Bacalhau et al. (2017). 
     Corrected typo on the documention of ncm_stats_dist2d.c.

 * New NcmDiff object that contains all numerical differentiation in a organized
     framework. Improved ncm_assert_cmpdouble test and error message.

 * Included more BAO points in the example. Minor modifications on the plot.

 * Created abstract class and one child to compute reconstruct an arbitrary
     two-dimensional probability distribution. Work in progress!

 * New smooth bpl model.

 * Documentation completed. Included PROP_SIZE in the stats_dist1d enumerator.

 * Removed debug message.

 * New Spherical Bessel FFTLog code.

 * Minor tweaks.

 * Improved the documentation on both hiprim examples.

 * Improved the documentation.

 * NcRecomb documentation is completed.

 * New full C example.

 * Few improvements on the recombinatio figures, new format svg.

 * Documentation and figure improvements.

 * Example better documented, plots added and small bugs fixed.

 * Fixed example name in Makefile.am

 * New example.

 * Modified the file name, included the computation of the halo mass function. All
     figures are in svg format.

 * Fixed mixing http/https in docs. Removed old doc from README.

 * Fixed manual url.

 * Updated automake and some examples, removed old docs (merged into the new
     site).

 * Added the initial and final masses and redsfhits as properties in the
     NcHaloMassFunction object.

 * Renamed the DE model Linder and Pad to CPL and JBP, respectively.

     Documentation improvements.

 * Documentation.

 * Fixed doc typo.

 * Removing old docs (merging all docs in a single place). Fixed bugs in
     recomb_seager.

 * Qlinear and Qconst -- documentation completed.

 * Documentation improvements.

 * Fixed log typo.

 * Added restart run in NcmFit. Added restart option in darkenergy. Reorganized
     NcRecomb and added tau_drag functions.

 * Removed typo from data file.

 * Removed -u option in cp since this flag is missing on macos.

 * Fixed doc typo and adjusted OmegaL range in HICosmoDE.

 * Improved the documentation for the xcor module and data object.

 * Fixing doc typos.

 * Adjusted r range and scale.

 * Implemented the object NcPowspecMLFixSpline: it computes the linear matter
     power spectrum from a file, which contains the knots k and their respective
     P(k) values.

     Included PROP_SIZE in the enumerator of some objects. Small improvements on
     the documentation.

 * Added example with tensor contributions to CMB.

 * Fixed indentation and moved the interface (XHeII and XHII) to NcRecomb (as
     virtual functions).

 * Added functions in NcRecomSeager to evaluate XHII and XHeII.

 * Updated values using data from https://sdss3.org/science/boss_publications.php

 * Fixed data object.

 * Last tweaks after merging with WL branch.

 * Removed before merging with WL.

 * Last tweaks after xcor merge.

 * Removing old files, preparing for the merge with xcor.

 * Fixing indentation before merge.

 * Fixing indentation in ncm_data_gauss_cov.c

 * Tested and fixed the last ensemble check in NcmFitESMCMC.

 * Added last ensemble check to ESMCMC.

 * Removed trailing space in Makefile.am.

 * Fixed the dumb error in .travis.yml and the NcScalefactor GO interface.

 * Trying another way to find the right gcc in travis-ci+macos

 * Fixing travisci build in macos.

 * Checking error in macos build.

 * New interacting dark energy model IDEM2

 * New interacting dark energy model IDEM2

 * Removed spurious - in HOAA. Finished the addition of a new BAO point.

 * Fully working version of HOAA for tensor and scalar modes of the Vexp model.

 * New interaction dakr energy model IDEM2

 * Included new BAO data point: Ata et al. (2017), BOSS DR14 QSO catalog.

 * Code reorganization and examples improvements.

 * Found a good parametrization for HOAA and a method to avoid roundoff during the
     transitions.

 * Still testing parametrizations in HOAA.

 * Testing parametrizations in HOAA>

 * Working on reparametrization of HOAA during singular transitions.

 * NcHICosmoDE documentation (header file).

 * Improve the documentation of some functions such that the bindings can be
     properly created. For instance, @lnM_obs: (array) (element-type gdouble):
     logarithm base e of the observed mass.

 * Finishing HOAA cleaning and testing.

 * Documentation fix.

 * Tweaking examples.

 * New examples and improvements on Vexp and HOAA.

 * Test gcc detection in travis-ci.

 * Cleaning and organizing NcmHOAA.

 * Minor tweaks.

 * Trying travis ci releases.

 * Missing header in toeplitz.

 * Better status control on _ncm_mset_catalog_open_create_file.

 * Fixed bug when reading a catalog with wrong mset fmap.

 * Better output notation for visual HW.

 * Fixed wrong parameter call in visual HW and added an assert to
     ncm_stats_vec_ar_ess to avoid future error like this.

 * Using the ensemble mean when the catalog has more than one chain.

 * New visual HW test added.

 * Fixed ar_fit 0 order case.

 * Missing reset status.

 * Fixed string allocation.

 * Support for reading fits + incompatible mset file.

 * Typo in mcat_analize.

 * Fixed typo in assert.

 * Improved script.

 * Updated .gitignore to include backup files and others.

 * Better asserts.

 * Improved interface with gsl minimizers. Included restarting for mms algorithms.

 * Improved example.

 * Chains diag output fix.

 * Missing refs and typo.

 * Added new diagnostics to NcmMSetCatalog, max ESS and Heidelberger and Welch's
     convergence diagnostic.

     Both can be applied to any NcmStatsVec object or through the NcmMSetCatalog
     interface. In the latter the test can be applied to individual chains, to
     the full catalog and to the ensemble average. Added an option to
     NcmFitESMCMC to automatically trim  the catalog during an ESMCMC run using
     the diagnostics to estimate the best burnin. Added options to run the
     diagnostics to mcat_analize.

 * Not allowing travis to fail on osx.

 * Removed duplicate sundials on travis+osx.

 * Removed klu support on sundials for osx.

 * Added another no-warning flag.

 * Better compilers warnings switches.

 * Added new object to docs.

 * Adding new Toeplitz solvers to docs ignore list.

 * New Toeplitz solvers added. Improved NcmMSetCatalog and NcmFitESMCMC. Added new
     helper functions in several objects.

 * Building docs on linux .travis.yml

 * Trying gcc-6 in .travis.yml

 * Added support for gcov.

 * Trying to rehash in osx.

 * Updated old finite call in levmar, ignoring errors in travis+osx.

 * Fixed all plc warnings and minor bugs.

 * Removed gcc recomp in .travis.yml

 * Reordered commands in .travis.yml

 * Removed CC export in .travis.yml

 * Asserting that gcc will be used in .travis.yml

 * Adding science deps on .travis.yml.

 * Other osx deps.

 * Trying to install gfortran via brew for osx travis.

 * Adding gfortran dep to travis osx build.

 * Adding deps for travis+osx.

 * Removed wrong dist-hook.

 * Log on check and dist.

 * Testing MACOS build.

 * Testing macos support for travis. Fixed minor doc typos.

 * Removed travis log output.

 * Still fixing doc building in travis.

 * Missing texlive package for travis doc compilation.

 * Log try typo.

 * Testing doc building in travis.

 * Better make mensages in travis.

 * Added latex support for travis.

 * Building docs on travis.

 * Fixed type warning in tests/test_ncm_integral1d.c.

 * Fixing last clang related warnings.

 * Fixed another set of minors clang related bugs.

 * Fixed several clang warning related minor bugs.

 * Added return to avoid warnings. Fixed multiple typedefs.

 * Fix plc's Makefile.am.

 * Fixed typo.

 * Conditional use of warning flags depending on the compiler. Fixed abs -> fabs
     bug in Planck likelihood.

 * Changing to make check.

 * Missing deps.

 * Testing trusty.

 * Removed update line.

 * Trying lucid.

 * Testing deps.

 * Testing dependencies .travis.yml

 * Removed debug gtkdocize on .travis.yml

 * Improved autogen.sh to work with old gtkdoc (and without it!).

 * Testing gtk-doc + .travis.yml

 * Testing .travis.yml

 * Including dependencies.

 * Testing travis.yml.

 * Improvements on NcHICosmoGCG. Testing new MCMC diagnostics and NcmStatsVec
     algorithms. Including support for Travis CI.

 * Created functions to obtain the expected means and the observed values.

 * New GCG model. Testing new diagnostics tool for catalogs.

 * Connected the knots vector of Poisson data with mass_knots.

 * Implemented data object for cluster number counts in a box (not redshift
     space). It follows a Poisson distribution.

     Implemented Crocce's et al. 2009 multiplicity function.

     NcmData: Better hooks for begin function.

     NcmDataPoisson: improved.

 * Using a warning instead of a assert in the final optimization test.

 * Testing better optimization finishing clean-up.

 * Created Crocce's 2009 multiplicity function.

     Created Cluster counts data in a box (not redshift space). In progress.

 * Fixed dependency link bug.

 * Removed log from bflike_smw.f90.

 * Moved prepare if needed to nc_hicosmo_sigma8.

 * Increased parameters scale.

 * Fixed typos and increased parameters scales in NcHIPrim*.

 * Fixed typo and increased lambdac range in BPL.

 * Updated c2 variables.

 * Added support for weighted observations in ncm_stats_dist1d_epdf. Added tests
     for ncm_stats_dist1d_epdf. New sampling functions in NcmRNG.

 * Added current time to (ES)MC(MC) logs.

 * Added gtkdocize to autogen.sh.

 * Added autoreset of the acc when splines are reset. Fixed warnings in CLASS
     lensing.c.

 * Added doc.

 * Fixed NcClusterMassAscaso compilation errors.

 * Improved example, added a child of NcmDataGaussCov.

 * Added a new parameter to Atan HIPrim model. Better (de)serialization for
     NcmMatrix. New serialization to binary file. Improved examples.

 * Created new cluster mass (relation provided in Ascaso et al. 2016). Work in
     progress.

 * Included additional parameter at the autocorrelation time calculation.

 * Fixed nlopt search libs.

 * Updated NLOPT library name from PKG_MODULE.

 * Update Dockerfile
 * Update Dockerfile
 * Added support from partial reset (only autosaved objects) for NcmSerialize.

 * Fixed conditional compilation for old GSL.

 * Fixed bug in cubic spline and removed PKEqual debug messages.

 * Fixed typos and commented old code.

 * Fixed typo in references.bib

 * Added PKEqual for HaloFit+Linder parametrization.

 * Imported improvements from xcor branch.

 * New deg2 to steradian convertion factor.

 * Modified NcXcor to select method for Limber integrals at construction.

 * Created function to compute the p-value of a function, giving the upper limits
     of the integral of the probability distribution function.

     Created option in mcat_analyze to compute the p-value of a function at
     different redshift values, giving the upper limits of the integral of the
     probability distribution function.

 * Correction to the Dockerfile for multi-threading.

 * Fixed a bug in Halofit

 * Fixed nc_hicosmo_de_reparam_cmb bug.

 * Switch xcor_limber integrals back to GSL (for now)

 * Updated Halofit (not tested yet)

 * Missing test file.

 * Fixed example neutrino masses. Fixed high-z neutrino calculations at
     NcHICosmoDE (needs improvement).

 * Missing doc tags.

 * Updated implementation flag on xcor.

 * Removed old files.

 * Added tests on CBE background. Pulled improvements on NcmSplineFunc from
     another branch. Improved speed on NcHICosmoDE using splines for massive
     neutrino calculations.

 * Fixed test.

 * Update examples to use massive neutrinos.

 * First tests OK. Working beta.

 * Fixed Omega_m usage.

 * Updated implementation flags code, and improving CLASS/NumCosmo comparison.
     TESTING VERSION!

 * CLASS updated to v2.5.0. Finishing massive neutrino interface and
     implementation on NcHICosmoDE.

 * Documentation NcmFit (in progress).

 * Replaced Omega_m0*(1+z)^3 by nc_hicosmo_E2Omega_m in several objects to take
     neutrinos into account.

 * Imported work in progress on NcHICosmoDE from xcor.

 * Working on singularity crossing.

 * WARNING : unfinished work on neutrinos in nc_hi_cosmo_de

 * Included GObject-introspection in the requirement list.

 * Improved massive neutrino interface. Split NcmIntegral1d.

 * Added missing author.

 * Imported Dockerfile from xcor.

 * Imported neutrino interface improvement from xcor branch.

 * Working on the massive neutrino interface.

 * Update Dockerfile
 * Update Dockerfile
 * Update Dockerfile
 * Update Dockerfile
 * Create Dockerfile
 * Initial commit for weak lensing branch.

 * Some debugging for interfacing neutrinos/ncdm with CLASS...

 * Added a minimal interface for massive neutrinos (will change in the near
     future), testing code. Simple implementation of this interface in
     NcHICosmoDE (not matching the CLASS background yet).

 * Work in progress Vexp, NcHICosmoAdiab and NcmHOAA.

 * Working version Vexp + HOAA + Adiab, it needs structure.

 * Typo corrections.

 * Added sincos detection to configure. New Harmonic Oscillator Action Angle
     variable object. Improvements on NcHICosmoVexp.

 * New Vexp model.

 * New mcat_join tool, it joins different catalogs of the same experiment.

 * Changed safeguard in nc_cbe

 * Missing doc tag.

 * Fixed indentation.

 * Organized and improved (testing phase).

 * Added a safeguard for halofit (Brent solver, in case fdf solver crashes).

 * Corrected a bug in halofit.

 * Corrected leaks in nc_data_xcor.c and modified ncm_data_gauss_cov.c in case of
     singular matrix.

 * Corrected some leaks in NcDataXcor.

 * Increasing maxsteps in xcor.

 * Modified xcor and halofit.

 * Fixed conflicts.

 * Fixed header name.

 * Corrections in xcor and cbe.

 * Fixed the close to the edge bug (emanating from CLASS).

 * Fixed inconsistencies.

 * Small fixes.

 * Increased output sampling of NcPowspecCBE to avoid interpolation errors.

 * Correction in nc_xcor.c and ncm_vector.h

 * Organizing code.

 * Organizing code.

 * Adding Xcor data objects.

 * Organizing and tweaking new Xcor objects.

 * Imported updated XCor codes. First tweaks and documentations fixes.

 * plop

 * Increased output sampling of NcPowspecCBE to avoid interpolation errors.

 * Correction in nc_xcor.c and ncm_vector.h

 * Organizing code.

 * Organizing code.

 * Adding Xcor data objects.

 * Organizing and tweaking new Xcor objects.

 * Imported updated XCor codes. First tweaks and documentations fixes.


[v0.13.3]
 * Updated changelog.

 * New changelog file.

 * Version bumped to 0.13.3.

 * Added smoothing scale to eval by vector function.

 * Added a smooth transition from non-linear to linear power spectrum for high
     redshift in halofit. Added the znl finder to obtain the redshift where we
     should stop applying the halofit.

 * Imported from xcor branch.

 * More stable safeguard.

 * Added get_tau from NcHIReion to MSetFuncList.

 * Added safeguard to minimization in ncm_stats_dist1d.

 * Removed debug printf in mcat_analyze.

 * Missing reference in nocite.

 * Small improvement on ncm_stats_dist1d_epdf and ncm_stats_vec. Included NEC on
     nc_hicosmo. Fixed core detection bug on configure.ac.

 * Included H(z) data: Moresco et al. (2016), arXiv:1601.01701.

 * Testing new algorithm in ncm_stats_dist1d_epdf.

 * Example with zt. Function to get cov from NcmFit.

 * Implemented function to compute the deceleration-acceleration transition
     redshift.

 * Improved reentrancy support when using ifort.

 * Added OPENMP flags log.

 * Included automatic flags for plc compilation. Fixed typo in ncm_fit_esmcmc.c.

 * Improved error handling.

 * Improved linear ps from CBE.

 * Log commented

 * Fixing docs typos and simple bugs.

 * Bug fixing.

 * Fixed setting mset parameters using vectors when some models have no free
     parameters.

 * Added H(z) data: Moresco (2015).

 * Removed verbosity at NcCBE.

 * Included BAO data: SDSS BOSS DR11 -- LyaF auto-correlation and LyaF-QSO
     cross-correlation. Modified nc_data_bao_dhr_dar.c to consider any number of
     data points.

 * New BAO object and new Sundials detection.

     Included detection for the new Sundials version (2.7.0) in configure.ac. 
     New BAO object based on D_M/r_d, H(z)*r_d estimates - SDSS BOSS DR12.

 * Updated the script mass_calibration_planck_clash.py, new funtion to include a
     gaussian prior. Corrected typos in the documentation.

 * Removed broken test for when old gsl is present.

 * Fixed double AC_CONFIG_MACRO_DIR.

 * Removed local link file.

 * Missing header.

 * Fixed g_clear_pointer workaround.

 * Fixing compiling bugs on opensuse.

 * Remove openmp flags from g-ir-scanner.

 * Another conditional compilation bug.

 * Fixed conditional include of ARKode.

 * Updated version.

 * Fixed max redshift in example_ca.py.

 * Updated examples.

 * Improved docs.

 * New ESMCMC example.

 * Fixed conditional threads compiling.

 * Missing ending string null.

 * Align.

 * Fall back to default files.

 * Conditional usage of gsl >= 2.2 functions.

 * Conditional use of gsl_sf_legendre_array_ functions.

 * Fixed docs typos.

 * Fixed doc typos. New catalog sampler NcmMSetTransKernCat. Removing warnings.

 * Added missing GSL support for darkenergy, updated mset_gen to generate mset 
     with models and submodels.

 * Missing gsl link for mcat_analyze.

 * Missing glib link to darkenergy.

 * Explicity link to glib in tools.

 * Test for files in cbe_precision.

 * Add lock to fftw plans.

 * Working on NcHIPertWKB (in progress, unstable).

 * No modification in example_hiprim.py

 * Fix for dlsym on macos

 * Fixed header inclusion (new gsl stuff).

 * Included sigma8 in NcmMSetFuncList. Small adjustements. Working in progress in
     WKB.

 * Updating WKB module (work in progress).

 * Temporary debug prints.

 * Small improvements.

 * Finishing the alm2pix transform.

 * Improving outsource compiling

 * Improving outsource compiling

 * Improving outsource compiling

 * Improving outsource compiling.

 * Removed unecessary comment.

 * Reordered class and instance structs, now all objects declare first the class
     struct and then the instance struct.

 * Support for require maximum redshift in NcDistance.

 * Changed the maxium redshift requirement of Halofit to match the maximum asked
     not the maximum non-linear.

 * mcat_analyse now outputs the full covariance when --info is enabled.

 * Added skip in unbindable functions in nc_cluster_mass_plcl.c. Support for
     including function in ESMCMC analysis through darkenergy. Added support for
     evaluating NcmMSetFunc in fixed points.

 * Many improvements and additions.

     Extended Serialize to work with some instances in different places. 
     Reworked NcmMSetFunc to be an object including description and symbols
     related to the function it calculates. NcmMSetFuncList introduces an
     generic catalog containing NcmMSetFunc from any observable. An initial list
     from NcHICosmo and NcDistance was already created, more to come. Reworked
     NcmPrior on top of NcmMSetFunc, now it can be serialized and saved to disk.
     Now priors can be implemented directly in Python. New generic abstract
     objects NcmPriorFlat and NcmPriorGauss. Improved TwoFluids objects and
     examples to match last paper equations of motion arXiv:1510.06628 (Work in
     progress). Moving all pixalization and spherical harmonics decomposition
     related code to NcmSphereMapPix (Work in progress). New tests on
     NcmSphereMapPix (test_ncm_sphere_map_pix). NcmObjArray now supports saving
     and loading from disk. New NcmObjArray tests.

     Several minor improvements.

 * Create README.md listing and describing the scripts. Included script
     mass_calibration_planck_clash.py (ref. arXiv:1608.05356).

 * Updating TwoFluids perturbation object, working in progress.

 * References included - documentation in progress.

 * Finalizing NcmSphereMapPix (including spherical harmonics decomp).

 * Peakfinder functions were rewritten in terms og GSL functions, therefore the
     objects NcClusterMassPlCL and NcCluster PseudoCounts no longer depend on
     the Levmar library.

     Documentation in progress: README, dependencies.xml, NcmSpline,
     NcmPowspecFilter, NcmFftlog.

 * New reorganized NcmSphereMapPix object.

 * Reorganizing quaternions and spherical map objects.

 * Removed debug messages.

 * Added support for ARB and included ARB calculation of NcmFFTLogTophatwin2.
     Fixed minor bug in FFTLog.

 * Documentation: work in progress.

 * Removed debug print.

 * Finalized the inclusion of NcPowspecMLNHaloFit and adaptating to
     NcmPowspecFilter. Added support for derivatives in NcmSpline2dBicubic.

 * Removed debug print from exaple_ps.py

 * Growth function adjusted in NcPowspecMLTransfer (~1.0e-4 precision comparing
     Class and EH at z = 0).

 * Removed old powerspectrum from NcHICosmo and moved everthing to NcHIPrim. All
     objects were adapted accordingly. (Work in progress!)

 * Documentation: ncm_fftlog, ncm_fftlog_gausswin2, ncm_fftlog_tophatwin2

 * Documentation: nc_cbe, nc_powspec, nc_powspec_ml, nc_powspec_ml_cbe,
     nc_powspec_ml_transfer

 * Updating NcHaloMassFunction to use the new NcmPowspec family.

 * Renamed NcMassFunction to NcHaloMassFunction

 * Reorganizing fftlog object, added calibration method to adjust the number of
     knots. New PowspecFilter object to apply filters (curretly gaussian or
     tophat) to any powerspectrum. Modifying example_ps.py (not ready yet).

 * Functions implemented: nc_cluster)mass_plcl_pdf_only_lognormal and
     nc_cluster_pseudo_counts_mf_lognormal_integral.

 * Minimal README for the python example.

 * Documenting python children objects.

 * Better doc in python example.

 * New Monte Carlo example. Relaxed Serialize to deal with python derived objects.
     New GObject frontend to random number generation functions.

 * Python mcmc example.

 * New external code Faddeeva for error function calc. New python example.
     Improved bandwidth in NcmStatsDist1dEPDF.

 * Missing data file.

 * New Hubble H_0 data Riess2016. Added helper function for CMB reparam.

 * Updating documentation: README, dependencies and compiling

 * New HIPrim models (broken power law and exponential cut).

 * Pseudo counts parametrized in terms of lnMcut, instead of lnTx.

 * Plot scripts update due to numpy modifications.

     Cluster pseudo counts - new variable Tx (substituting the old one - Mcut).

 * First tests with Planck polarization likelihood.

 * Added support for TE EE data from Planck likelihood.

 * Bug fix in ncm_fit_esmcmc_walker_stretch.

 * Added missing object registry.

 * New reparametrization nc_hicosmo_de_reparam_cmb and nc_hiprim_atan. Better
     handling of border cases in ncm_fit_esmcmc_walker_stretch.

 * Improved parallelization of NcMatterVar by removing a mutex. Fix bug in
     NcClusterPseudoCounts. Improved border handling in
     NcFitESMCMCWalkerStretch.

 * Cleaning wrong annotations, added individual shrink factors calculation in
     mcat.

 * Update parameter in ncm_mset_catalog to improve the shrink factor calculation.

 * Implemented selection function considering the relation between X-ray
     temperature and true mass. Defined new parameter:
     NC_CLUSTER_PSEUDO_COUNTS_LNTX_STAR_CUT

     Example example_hiprim.py has a bug.

 * Support for changing kmax and kmin in NcPowspecMLCBE.

 * Missing file.

 * New NcPowspecMLCBE for extracting linear matter power spectrum from CLASS.
     Trying new walkers (and options) for ESMCMC. New options for ESMCMC added
     to darkenergy. New example example_ps.py of how to use new NcPowspecML
     objects. Organized code in NcTransferFuncEH.

 * File simple corner.

 * Removing pyc files.

 * Fixing objects definition order (just cosmetics). Improving NcmMSetCatalog (now
     supports printing ensemble time evolution). Designing new NcmCalc abstract
     object.

 * Added plot scripts.

 * A typo and a leak.

 * Restructured NcmFitESMCMC.

     Now NcmFitESMCMC supports different walkers through NcmFitESMCMCWalker
     interface. Both serial and parallel versions produce the same result.

     mcat_analize calculates now the integrated autocorrelation time with the -i
     option. darkenergy support for Planck data.

 * Fixed bug in mcat_analyze.c.

 * Imposed the same out-of-interval prior in both serial and parallel modes.

 * Moved back the default value of the parameter A of NcmFitESMCMC to 2.

 * Several fixes and improvements.

     New quantile support in ncm_stats_vec (using gsl implementation). Improved
     synchronization robustness for NcmMSetCatalog. Modified
     ncm_stats_dist1d_epdf to use silverman's rule of thumb to calculate the
     bandwidth. Added support to Planck likelihood usage through darkenergy. 
     Made clik_gibbs_f90 thread safe to avoid multithread conflicts.

 * Finished NcmIntegral1d first interfaces for Hermite and Leguerre like
     integrals. Added a new test for NcmIntegral1d.

 * Updated dependency on glib to version 2.32.0 and cleaned old legacy code.

 * Added Gauss-Hermit integration to NcmIntegral1d.

 * New NcPowspecML object for abstract linear matter powerspectrum. New
     NcmIntegral1d object for generic one dimensional integration. Organized
     NcTransferFuncBBKS internally.

 * Fixed documentation.

 * Bug fixed: normalization of the Planck and CLASH masses distributions are now
     implemented considering M_PL >= 0 and M_CL >=0. Normalization is given in
     terms of the error functions. This modification was done for the
     computation of the 3D integral!!!

 * Fixing examples. New example example_epdf1d.py added.

 * Moved NcPowerSpectrum to NcmPowspec (more general base object).

 * Added new abstract class for powerspectrum implementation.

 * Finished gitignore organization.

 * Organizing gitignore to clean the index.

 * Adding .gitignore to the repo.

 * Missing ChangeLog in libcuba.

 * Updated example out filename.

 * Added submodel support for NcmModelCtrl and finished the transition for
     submodels in all derived objects.

     Added tests for submodel in NcmModelCtrl.

 * Added submodel concept in NcmModel. NcHIReion and NcHIPrim are now submodels of
     NcHICosmo.

 * Renamed submodel for stackpos (stack position) in NcmMSet internals.

 * Finished resampling for NcDataPseudoCounts and its tests.

 * Added set_cad function in DataClusterPseudoCounts.


[v0.13.1]
 * Last updates in examples. ChangeLog updated.

 * Fixed docs and updated ChangeLog.

 * Updated and improved PseudoCount related objects. Updated of all examples
     finished. Fixed NcmMSet typo. Added accelerated bsearch option for
     NcmSpline2d.

 * Updating examples and organizing prepare calls in calc objects.

 * Bug fixed: ncm_data_set_init(...) was included in
     nc_data_cluster_pseudo_counts_init_from_sampling().

 * Updated ChangeLog

 * Bumped to v0.13.1

 * Now using ax_cc_maxopt to detect the best optimization flags (removing
     fast-math if included).

 * Removed dependency in Sqlite3.

 * Moved all Hubble data from sqlite3 to .boj files.

 * Moved all SNIa data from SQLite to obj files.

 * Fixed names of nc_data_cmb_wmap?_shift_param.obj files. Moved all distance
     priors data to obj files.

 * Moved shift parameter data to obj files. Fixed bug in NcmLikelihood. Removed
     old shift parameter constants in NcmC.

 * Added stackable and nonstackable models option.

 * Fixed doc not including NcmSplineCubic*. Finished support for NcmSpline2d
     serialization. Fixed typo in NcPlanckFI. Advanced in the class background.c
     replacement. Added 4He Yp from BBN interpolation table in NcHICosmoDE
     models.

 * Reorganized functions names in NcHICosmo, mostly for documentation reasons.

 * Internal reorganization.

 * Fixed typo in nc_recomb_seager.c (missing 1/3 factor).

 * Added more doc in NcRecombSeager and removed old code.

 * Finished all He switches in NcRecombSeager. (Working in progress)

 * Updating recombination code to match the theory used in recfast 1.5.2.

 * Updating constants in NcmC namespace.

     Updated CODATA constants from the last release of 2014. Included concise
     atomic weights from IUPAC. Organized astronomical constants using IAU
     recommendations. Included concise atomic spectra from NIST database. 
     Included references for the constants. Updated constants.txt to reflect the
     2014 CODATA release.

     Updated/added documentation for most functions in NcmC.

     Updated all objects to conform to the new names for the constants in NcmC.

     Updated CLASS thermodynamics.c, arrays.h and arrays.c.

     Working in progress in joining/comparing CLASS thermodynamics and NcRecomb.

 * Implemented function to compute 1-3 sigma error bars for the best fit.

 * Added message to be print when mode_error (mcat_analyze) is called.

 * mcat_analyze: implemented options mode_errors and median_errors. They provide
     the mode (median) and the 1-3 sigma error bars of a parameter.

 * Implemented functions to perform Planck-CLASH analyses considering flat priors 
     for the selection and mass functions.

 * Fixed typo.

 * New NcHIReion* objects. Moved Yp to the cosmological model NcHICosmo.

     Added new NcHIReion* objects to implement reionization models. Included
     CAMB like reionization NcHIReionCamb. The reparametrization object
     NcHIReionCambReparamTau permits the
      usage of tau_reion as the parameter for NcHIReionCamb.

 * New reionization objects (in development). Fixed minor bugs (including bugs in
     libcuba). Unstable boltzmann codes (in development).

 * Fixed sampling function of nc_data_cluster_pseudo_counts (and
     nc_cluster_mass_plcl).

 * Script to perform ESMCMC analysis of the Planck-CLASH clusters.

 * Added build hook to copy modified doc files to the building directory.

 * Fixed types in test_nc_recomb.c.

 * Missing file.

 * Missing files.

 * Improved documentation.

 * Documentation about GObject (basic concepts).

 * Updating recombination code.

 * Documentation

 * Removing support for clapack usage (some headers are broken).

 * Fixed minor bugs,

 * Fixed typos on README.md file. Imrpoved documentation. Function
     nc_data_cluster_pseudo_counts_init_from_sampling created. Sampling of
     cluster pseudo counts is working.

 * Better organization of Bolztmann code options and NcHIPrim implementation
     example.

 * Fixed all virtual functions in abstract classes to be recognized as such by the
     GObject introspection. Created new NcmModelBuilder to create NcmModel from
     binded language. Added a new example for this new feature.

 * Made gtkdoc optional (testing).

 * Improved NcmDataGaussCov tests.

 * Added support to GSL-2.0.

 * Removed spurious print from NcmModelTest. Added control on OPENMP in
     ncm_cfg_init.

 * Fixed parameter name in NcCBEPrecision.

 * Fixed memory leaks and vector parameter allocation in NcmMSet (very obscure
     bugs only active in weird cases). Added name and nick for every Model for
     debug purposes.

 * Fixed leak in numcosmo/nc_cbe_precision.c. And improved tests.

 * Updated macro NCM_TEST_FREE to use a safer method.

 * Add VERBOSE = 1 in make check.

 * New function ncm_mset_trans_kern_gauss_set_cov_from_rescale.

 * Fixed reallocation problem.

 * Added gi.require_version in python examples. Fixed opendir leak in libclik
     (plc-2.0).

 * Fixed typo.

 * Fixed return statement in clik_get_check_param.

 * Fixed fprintf usage in class.

 * Removed data repetition.

 * Fixed reference.

 * new NcHIprimAtan object (primordial spectrum power law x atan). New mset_gen
     tool to generate .mset files. Added flag controling the tensor mode usage
     in NcHIPertBotlzmannCBE. New references (to the atan models). Bumped to
     version 0.13.0.


[v0.13.0]
 * new NcHIprimAtan object (primordial spectrum power law x atan). New mset_gen
     tool to generate .mset files. Added flag controling the tensor mode usage
     in NcHIPertBotlzmannCBE. New references (to the atan models).

 * Better error mensage when trying to de-serialize an invalid string.

 * Fixed parameter name in NcHICosmoDE z_re -> tau_re. Fixed NcHICosmoBoltzmannCBE
     to account correctly the lmax when using lensed Cls. Added free/fixed
     parameter manipulation functions to NcmMSet.

 * Missing HIPrim implementation PowerLaw (nc_hiprim_power_law).

     Working in progress in the Planck+CLASS interface.

 * New objects and support for primordial cosmology NumCosmo <=> CLASS.

     New NcHIPrim object to deal with primordial cosmology. New NcHIPrimPowerLaw
     object to implement simple power
      law primordial spectra.

     Organized main documentation page. Added support for external callback
     function in CLASS so
      it can call NumCosmo to get primordial spectra. Wired CLASS so it uses
     NumCosmo NcHIPrim to get the
      primordial spectra.

     Testing results of NumCosmo + Planck + CLASS.

 * First working version of the Planck+CLASS interface. Minor bug fixes.

     New NcHICosmoDEReparamOk to better handle Omega_x -> Omega_k
     reparametrization. Removed old code and updated NcmReparam. Moved precision
     data for CLASS Backend CBE to NcCBEPrecision. Added tests to check
     consistency of NcHICosmoDE. New parameters ENnu effective number of
     neutrinos to NcHICosmoDE. New methods of NcmHICosmo nc_hicosmo_Omega_g and
     nc_hicosmo_Omega_nu
      to account for electromagnetic density and ultra-relativistic
      neutrinos density.

     Finished the interface NcHIPertBoltzmannCBE which uses CLASS
      as backend for perturbative computation, needs polishing. Finished the
     interface for Planck likelihood NcDataPlanckLKL.

 * Working on the CLASS interface. All precision parameters mapped.

 * Initial phase of the Class backend interface.

     Modified ncm_cfg_get_data_filename so it also search in PACKAGE_SOURCE_DIR
     for data files. That way a non-installed and non-configured numcosmo can
     run using the data from the source directory.

 * Added Class as backend. Documentation fixes. Renamed object NcPlanckFI_TT to
     NcPlanckFICorTT.

 * Missing files in the last commit.

 * New NcPlanckFI objects.

     New NcPlanckFI model type and one implementation NcPlanckFI_TT were
     created. These models implement the necessary parameters to deal with the 
     Planck likelihood.

     Included the conection between NcPlanckFI* model and NcDataPlanckLKL
     likelihoods.

 * Updated test_nc_cluster_pseudo_counts.

 * Bug fix and initial object development.

     Added support for extra flags cheking in Fortran. Added new object
     NcDataPlanckLKL, almost finished laking the
      NcmData implementation.

     Fixed two bugs in bflike_smw.f90, now it works with -O3 + gfortran.

 * Fixed last steps for making releases.

 * Added Planck likelihood 2.0 to the building system.

     Organized the PLC likelihood nested in the NumCosmo building system. Added
     the Fortran prerequisites in order to allow parallel compilation.

 * Updated to internal libcuba 4.2.

 * Resample function of nc_cluster_data_pseudo_counts is a work in progress. 
     NcClusterPseudoCounts object has a new property: ncluster - number of
     clusters.

     NcClusterMass and NcClusterRedshift are no longer properties of
     NcClusterAbundance. Examples still need to be update!

 * Fixing minor bugs.

 * Fixed bug .

 * Support for more sundials 2.6.x versions.

 * Fixed tests and removed debug prints.

 * Testing new parametrizations in hipert_two_fluids.

 * Deleted functions related to the 3-dimensional integral on
     nc_cluster_pseudo_counts.c and renamed all functions of the new
     3-dimensional computation removing the label "_new_variables".

 * Functions to compute the 3-dimensional integral over the true, SZ and lensing
     masses are working. There are two set of functions to compute it
     (independently). The main diference between them is the set of integral
     variables: 1) logarithm base e of the masses and  2) new variables (we
     performed a change of variables). The later provides the best results and
     it is in agreement with the 1+2 integral (integration over true mass and a
     bidimensional integration over SZ and lensing masses) for any values of the
     parameters.

 * Cleaned the code

 * Implemented limber approximation for cross-correlations and likelihood analysis

 * Minor modifications: in progress!

 * Added missing data file.

 * Fixed typo. Working in progress...

 * Several improvements. New sub-fit support.

     Updated configure to detect new versions of SUNDIALS. All dependent code
     updated accordingly.

     New subsidiary fit support. Before each step in a fitting process, it will
     fit a subsidiary likelihood using a subset of the parametric space.

     Automatic convertion between variant vector of doubles and NcmVector and
     variant matrix of doubles and NcmMatrix. It is no longer necessary to
     manually convert NcmVector and NcmMatrix to variant when using it as object
     property.

     Updated all objects (but NcDataClusterPseudoCounts) to use
     NcmVector/NcmMatrix directly as properties instead of their GVariant
     versions.

 * Documentation improvements (in progress): some NcClusterMass' children,
     ncm_abc.c and ncm_lh_ratio1d.c.

     Implemented bidimensional integrations [divonne function (libcuba)] using
     and not using peakfinder. Idem for tridimensional integration with no
     peakfinder.

     Implemented functions to compute the integrand peaks in
     nc_cluster_mass_plcl and nc_cluster_pseudo_counts: testing! Created
     test_nc_cluster_pseudo_counts: in progress!

 * Added support for new version of SUNDIALS.

     Working on NcHIPertBoltzmann.

 * Fixed bootstrap support in NcmDataDist1d

 * Moved BAO data from hardcoded to serialized objects. New BAO data. Minor test
     updates.

     All BAO data now are included as serialized objects. New tests in
     test_ncm_mset. New NcDataBaoDHrDAr object for (D_H/r_zd, D_A/r_zd) data. 
     New field "long-desc" in NcmData to include detailed description of data. 
     Removed old data BAO from NcmC.

 * Fixed bugs in ncm_lh_ratio1d from last update in this object.

     Added shallow copy to NcmMSet and its tests.

 * Added more testing in test_ncm_mset, fixed minor bugs.

 * New support for multples models of the same type in NcmMSet. Minor fixes and
     updates.

     Created test suit for NcmMSet.

 * Documentation: improvements on nc_cluster_redshift, nc_cluster_mass and
     nc_hicosmo.

     Unstable version of nc_cluster_mass_plcl.

 * Minor update.

     Updated libcuba (4.1 => 4.2). Fixed typos and simple bugs. Added test for
     NcmFuncEval. Added tests for finite values of m2lnL in ESMCMC.

 * Added check for missing set/get functions in NcmModel.

 * Pseudo cluster number counts: observable and data objects were created. 
     Integral new function: function to compute tri-dimensional integral
     implemented using cuhre function (libcuba). Documentation: improvement in
     different files.

 * Added support for jerk in DE models.

 * Implementing Planck-CLASH mass function: in progress.

     Documentation: nc_cluster_mass.

 * Bumped to v0.12.2

     Fixed typo in ncm_stats_dist1d.c. Updated Changelog


[v0.12.2]
 * Bumped to v0.12.2

 * Fixed bug in ncm_fit_esmcmc_run_lre.

 * Added new lnsigma_lens parameter to NcSNIADistCov.

     New release of JLA data available on NumCosmo site, full sample without any
     additional variance included (everything now is included by the model).

 * Tools reorganization and several improvements.

     Moved darkenergy to tools directory. Removed from darkenergy options for
     catalog analyze.

     Created mcat_analize to perform analysis on Monte Carlo catalogs
      including MC, MCMC, ESMCMC and Bootstrap MC.

     New NcmStatsDist1dEPDF for 1-dimensional empirical distributions,
      it implements a fast Gaussian basis function interpolation for
      general 1d distributions, recreating the pdf, cdf and inverse cdf.
      Allowing fast resample and quantiles calculation.

     Added hash table interface to access functions from NcHICosmo and
     NcDistance. All these functions are accessible from mcat_analyze, allowing
     it to generate quantiles, distributions, cdf for all these quantities.

 * Added missing docs directives.

 * Improved interface to NcmLHRatio2d in darkenergy. Minor improvements.

 * New cluster mass relation: Planck-CLASH correlated mass-observable relations.
     Planck - SZ signal. CLASH - lensing signal. This object is not finalized.
     Work in progress.

     The function "ncm_model_class_set_name_nick" shall be called before
     "ncm_model_class_add_params". Fixed this order in various models.

     Documentation (improvement): ncm_stats_dist1d.c, ncm_stats_dist1d_spline.c,
     ncm_model.c, ncm_sparam.c, ncm_sparam.h.

     Removed spurious line (sd1->norma = 1.0;) from ncm_stats_dist1d_prepare
     function.

     Improved message error in ncm_model_class_set_sparam and
     ncm_model_class_set_vparam functions.

 * Missing file.

 * Optimization flags.

 * Added and enable flag to include compiler's optimization/warnings flags. Made
     several minor code quality improvements.

 * Added a internal version of Cuba. Fixed minor typos and updated autogen to use
     autoreconf.

     It is no longer necessary to have an installed version of Cuba. If it is
     not available in the system the library compiles its own internal version
     of Cuba. The same was done recently for levmar. Updated README.md to
     explain these points.

     Removed cubature lib support (now it always use Cuba).

     Added option for different inclusing of redshift errors in NcSNIADistCov.

     Removed the call for cholesky_decomp in cholesky_inverse, now the user must 
     call both.

 * Better workaround for the missing fffree/fits_free_memory functions and
     SUNDIALS_USES_LONG_INT macro. Corrected version for g_test_subprocess
     usage.

 * Fixed threads competitions with OpenBLAS or MKL. Finished the NcDataSNIACov
     interface.

 * New minor version. Several improvements.

     Better building tool support for lapack/BLAS. Improved speed using
     different interfaces for Lapack/BLAS. Included levmar internally and
     removed as optional dependency. Levmar intreface now supports box
     constraints. Added virtual functions for lnNorma2 and resample in
     NcmDataGaussCov. Updated interval vector usage in NcmFitNLOpt and
     NcmFitLevmar. Added support for estimating true width and colour in
     NcDataSNIACov. Improved resampling in NcDataSNIACov, now it resamples all
     data m_B, w and c. Several new functions in NcmMatrix (support to writing
     data in col-major order). New ncm_mset_catalog_calc_ci for calculating
     confidence intervals for generic function using a NcmMSetCatalog. Added
     option for confidence intervals from catalogs in darkenergy. Added
     templates as dependency for NcmFitNLOpt enums files.

     NcDataSNIACov now uses the (hopefully) the right normalization (new feature
     still experimental).

     Added some misc statistical functions in NcmUtil (experimental).


[v0.12.1]
 * New minor version. Several improvements.
 * Added option for unordered MC runs. Fixed typo in parameter of NcSNIADistCov.
 * Several minor fixes and improvements.
 * Formated section documentation for all object, some simple doc fixes.
 * Added option to OGbject introspection scanner use the same CC flag.
 * Added checks for file open/close.
 * Documentation: changed titles and short description: BAO
 * Missing files.
 * Added missing references.
 * New test: test_nc_data_bao_dvdv.c Added nc_cor_cluster_cmb_lens_limber.h in
     numcomso.h
 * Fixed typo in name nc_data_bao_empirical_fit. Support for data filenames. New
     BAO data.
 * Improved NcmSpline is now serializable. New ncm_stats_dist1d_* family. New BAO
     data nc_data_bao_empirical_fit. Alpha^3 version of the new Boltzmann
     object.
 * Created test for object nc_data_bao_rdv.
 * Include catalog description in the fits in NcDataSNIACov and fixed minors
     leaks.
 * Added dataset minimum id check in NcDataSNIACov.
 * Fixed documentation.
 * Fixed missing docs tag in nc_snia_dist_cov.h.
 * Added SNIa from SDSS-II/SNLS3 ( arXiv:1401.4064 ).
 * Object ncm_lapck.c was documented.
 * Fixed other compilation problems.
 * Fixed some minor compilation problems.
 * New Ensemble Monte Carlo object added, several minor fixes and MC codes
     organization. Work in progress on nc_hipert_two_fluids family.
 * Testing twofluids_wkb. Added precision property in nc_mass_function, default
     10^{-6}.
 * Better error message in ncm_fit_new and fixed wrong parameters in
     example_simple.py and example_simple.c.
 * Added README to conform with automake. Fixed verbosity in NCM_TEST_FAIL and
     NCM_TEST_PASS.
 * Added links to README.md
 * Corrected file.
 * Added generic INSTALL file.
 * Removed old README.
 * Reorganized data objects. Improved README.md and docs.
 * Fixed bug in g_clear_pointer usage.
 * Added support for old cfitsio >= 3.25.
 * Added traceback for error messages. Testing differents summaries in
     NcABCClusterNCount.
 * Added support for the new version of libcuba 4.0.
 * Trying different ABC summary statistics.
 * Removed old gdarkenergy from building. Scale fisher matrix by two when using in
     the transition kernel in mcmc.
 * Removed debug message.
 * Fixed bug, trying to compile a vala source without vala available.
 * Fixed memory leaks in NcABCClusterNCount and NcClusterAbundance. Some new
     additions and more stable ABC code.
 * Better error message in prepare_base.
 * Updated the threaded evaluation function.
 * Fixed binning methods.
 * Fixed documentation in NcABCClusterNCount.
 * Added more options for binning in nc_data_cluster_ncount and
     nc_abc_cluster_ncount.
 * Fixed bug in ncm_abc.c (wrong number of additional columns). New main-seed and
     nthreads options in darkenergy.
 * Improvement on ABC interface. Added seed and nthreads options to darkenergy.
 * Fixed more compilation problems (when lacking of fftw3l).
 * Fixed some compilation errors for systems without fftw or cfitsio.
 * Missing file in EXTRA.
 * Fix several minor memory leaks, racing conditions. Functional version of NcmABC
     and NcABCClusterNCount.
 * Improvement on NcmSplineCubicNotaknot and NcmABC.
 * Improved blas/lapack search.
 * Fixed non-defined variable for lapackless systems.
 * Fixed lapack header inclusion.
 * Better support for blas and lapack search in configure and code organization
     and new NcmABC.
 * Fixed: property seed is no longer G_PARAM_CONSTRUCT (ncm_rng.c).
 * Fixed wrong casting.
 * Organized mc samplers in NcmMCSampler abstract object. New methods for Vector,
     Matrix, Model and MSet. Fixed bugs in NcClusterPhotozGaussGlobal.
 * Fixed nc_hicosmo_de prototypes for better bindings. New get/set functions by
     parameter names.
 * New example ode_spline. NcClusterRedshift transformed in NcmModel.
 * Minor fixes, release candidate 0.12.0rc1.
 * Fixed tests (using g_test_trap_fork unitl glib < 2.40).
 * Several fixes for release candidate 0.12.0rc0.
 * Fixed references.bib.
 * Updated ChangeLog.
 * Improved nc_data_cluster_ncount and fixed documentation typos.
 * New WKB codes and improvements in the perturbations code. Bumped to 0.12.0.
 * Improving adiabatic perturbations code.
 * New perl bindings example. Fixed bugs in ncm_spline.c. New expermental code in
     adiab.
 * Improved example_ca_sampling.py.
 * Added methods to exam a NcDataClusterNCount contents. Improved
     example_ca_sampling.py.
 * Added new example to Makefile.am
 * Fixed parameters setting order.
 * Removed testing code.
 * Better adaptation of NcCluster* to generating better bindings.
 * Added missing files.
 * Fixed documentation issues.
 * Improved README of examples. Introduced NUMCOSMO_DATA_DIR environment variable
     to allow running darkenergy with data without installing the library.
 * Fixed bug (infinity recursion to compute dE2/dz) in nc_hicosmo_qspline.c.
     Documentation of ncm_fftlog.c is partially done.
 * Finishing conversion of perturbation code to interfaces.
 * Fixed warnings.
 * Fixed header paths.
 * Create objects related to the matter density profile (abstract and NFW) and the
     computation of the cross corelation between clusters and CMB lensing
     potential. These codes are in development.
 * Improved tests. Organizing old code. New perturbations code.
 * Several additions.
 * Fixed many memory leaks.
 * Added ncm_func_eval_threaded_loop_full to run one worker per index.
 * Fixed memory leaks in serialization and minor bugs.
 * Two arXiv references were added as comment. Modify redshift from Beutler et al.
     2011 (before z = 0.1, now z = 0.106).
 * Adapting tests to conform to the new g_test_trap_subprocess. Fixed cvodes/cvode
     usage.
 * Removed gtester support.
 * Fixed bugs in testing.
 * Fixed gtkdoc's and introspection warnings.
 * Lower bound of the Dsz parameter from NcClusterMassBenson model is now 0.01.
 * Several minor improvements.
 * Added support for general gaussian priors in darkenergy.
 * Fixed minor bug (params_max and params_min were not being allocated).
 * Removed old code from nc_cluster_mass_lnnormal.c.
 * New NcmFitCatalog and NcmFitMCBS objects.
 * Comments removed in these files.
 * NcClusterMassLnnormal has now two properties: bias and sigma.
 * Added cpu core counting to set NTHREADS automatically.
 * Added support to libcuba 3.3.
 * Bug fixed.
 * Angular reduction when looking for bounds to avoid infinity repetition.
 * Fixed bug in NcmData which didn't called begin when a sample was generated by a
     resampling.
 * Fixed but in NcDataClusterNCount which discarded old references of
     NcClusterMass and NcClusterRedshift.
 * Fixed bug in darkenergy (always setting params reltol to 1e-5).
 * Added support for setting reltol and params-reltol in the NcmFit object.
 * Added global variables initialization for gsl_rng functions.
 * Fixed bug in ncm_func_eval_threaded_loop.
 * Fixed lock/unlock problem in NcDataCluster resampling.
 * Fixed miscellaneous bugs with valgrind (memcheck and helgrind).
 * Finished support for multithreading montecarlo and bootstraping.
 * Several improvements.
 * Several improvements and new objects.
 * Adding support for bootstrap in ncm_data_gauss_cov. Improved continuity prior
     on nc_hicosmo_qspline (using three knots with five points straight line
     fitting for continuity prior).
 * Inverted Class struct possition to avoid gtk-doc's bug.
 * Modified qspline continuity prior to fit line using n points for each three
     knots.
 * Separeted two types of priors in NcmLikelihood chisq and m2lnL. Modified
     continuity prior to use m2lnL priors. Modified darkenergy to receive the
     snia_cov serialized object.
 * Transformed continuity priors in NcmModel to fit the prior variance. Added
     automatically (set|get)_property to object_class in
     ncm_model_class_add_params and check for right functions. Changed models
     and tests accordingly.
 * Fixed compilation error in 32bits platforms.
 * Minor version increased, adapted glib versioning system.
 * Fixed continuity constraints in hicosmo_qspline. Improved continuity prior by
     using linear fitting to aproximate three points by a straight line.
 * Changed from g_hash_table_contains to g_hash_table_lookup != NULL to work with
     older glib.
 * Testing a new penalty function for overfitting in NcmDataGaussCov. Fixed
     compilation with new libcuba release. Improving documentation.
 * Still testing.
 * Testing new continuity priors in HICosmoQSpline.
 * Fixing documentation issues.
 * Reorganized NcSNIADistCov and NcDataSNIACov. Now all data is allocated in
     NcDataSNIACov and only model parameters stay in NcSNIADistCov. This fix the
     montecarlo with fiducial model issues.
 * Added an assert to NcSNIADistCov to check if the data is loaded.
 * Fixed bug model changing in NcDataSNIACov, now it work with alternating
     NcSNIADistCov models.
 * Testing different continuity priors in HICosmoQSpline. Added missing doc tags.
 * Added support for constraints in NcmFit and multiple algorithms in NcmFitNLOpt.
 * Improved numerical differentiation calculation through better choice of steps.
 * Removed wrong documentation.
 * Better update control in NcHICosmoQSpline object (fixed bug).
 * Fixed parameter name in function documentation.
 * Transformed QSpline continuity prior in an object. Added subdir-objects to
     automake.
 * Included numcosmo/build_cfg.h  in every header to generate bindings for
     conditional compilation functions.
 * Included missing header.
 * Improved catalog functions in nc_data_snia to allow better bindings.
 * Added missing parameters documentations in mset macros. Removed debugging
     g_error.
 * Modified model id framework, now the ids are exported by functions to allow
     better bindings through introspection.
 * Added function to count number of named instances and correct type annotation
     for named instances functions.
 * Improved NcmFitMC, log messages and support for different fitting and fiducial
     models.
 * Fixed free error in darkenergy.
 * Support for named instances, a global object pool. Added serialization of named
     instances. Added mset_load/save method for mset serialization. Added
     save-mset and fiducial options to darkenergy to allow saving   NcmMSet used
     in a analysis and defining a arbitrary fiducial   model for Montecarlo
     studies. Better organization of model registry and id.
 * New test: test_ncm_data_gauss_cov, testing resample and sanity.
 * Working on fftlog. Added kinematic functions for DEC and WEC. Added DEC and WEC
     tests to darkenergy option --kinematics-sigma.
 * Added property maxiter and method to change it in NcmFit. Connected this method
     with the --max-iter option in darkenergy. Included kinematic output in the
     --out option in darkenergy.
 * Imported some code from glib to allow serialization under glib < 2.30. Fixed
     warnings for compilations without sqlite3.
 * Workaround to g_clear_object usage.
 * Added workaround to check for sundials header correctly.
 * Implemented function to compute the inverse of the square normalized Hubble
     function (nc_hicosmo_Em2). Modified q-sigma, q-n and q-z-max DE options to
     kinematics-sigma, kinematics-n and kinematics-z-max. The kinematics-sigma
     option computes the deceleration parameter, the squared normalized Hubble
     function, its inverse and their error bars via Fisher Matrix approach.
 * Fixed default option for qspline continuity priors.
 * Several additions.
 * Implemented function nc_hicosmo_Omega_mh2. Functions
     nc_distance_decoupling_redshift and nc_distance_drag_redshift now use
     nc_hicosmo_Omega_bh2 and nc_hicosmo_Omega_mh2. Function
     nc_distance_dsound_horizon_dz was included in nc_distance.h. Documentation
     of nc_distance.c in progress (approximately 2/3 ready).
 * Updated ChangeLog.
 * Removed functions: nc_distance_curvature, nc_distance_comoving_a0 and
     nc_distance_comoving_a0_lss. nc_distance documentation im progress.
 * New data included. Silent rules by default during make.
 * Move back the functions p_limits and n_limits to be computed with 7 * sigma.
 * Added support for NcmVector and NcmMatrix serialization and their repective
     tests.
 * Updated integration on zeta true. The gap (1 < zeta < 2) was removed and the
     integration is performed in the entire interval. This is different from SPT
     code (their normalization does not take into account this gap), but it is
     consistent with the normalization used.
 * Improving tests and examples.
 * Improving tests and examples.
 * Added support for repeated options in darkenergy. New H(z) and BAO data.
     Support for non darkenergy models in darkenergy application.
 * Updated to glib 2.36, g_type_init no longer required. Added macros for testing
     for older versions.
 * Baryonic density (Omega_b) Gaussian prior from Big Bang Nucleosynthesis (BBN)
     was implemented.
 * Fixed bug: _ncm_data_gauss_resample now correctly uses the lower diagonal of
     the Cholesky decomposition.
 * File nc_data_cmb_dist_priors.c was documented.
 * Included WMAP 9 year distance priors.
 * Splitting namespace in NumCosmo and NumCosmoMath.
 * Testing.
 * Improving documentation, removed old files and fixed msg in autogen.sh.
 * Improvements to allow better bindings.
 * Better code for prereq finding during configure.
 * Copied the values of the magnitude from NcSNIaDistCovt to NcDataSNIaCov.
 * Testing SPT fitting.
 * New example: simple SN Ia model fitting.
 * Fixed g-ir-scanner sources argument.
 * Changed vector parameter in models to GVariant to improve serialization.
 * Added examples to dist and installation.
 * Updated ChangeLog.
 * Better organization of SN Ia data. Added BAO data.
 * Several bugs fixed. New supernovae Ia data with covariance matrix.
 * Fixed bug: Modified the number of bins to build the histogram of -2lnL (Monte
     Carlo). Fixed bug: set has_covar fit member equal to TRUE in function
     ncm_fit_mc_mean_covar.
 * Several improvements in documentation.
 * Removed enum doc.
 * Reorganized the conditional compilation of levmar and nlopt. Better solution
     for nlopt header processing.
 * Removed bugged unnecessary copy in numcosmo/Makefile.am.
 * Translated confidence region in NcmLHRatio1d and NcmLHRatio2d objects. Improved
     two dimensional confidence region algorithm. Added test in NcmOdeSpline to
     detect integration problems. Added test in NcGrowthFunc to detect
     integration problems. Bumped to version 0.9.0. Organized Monte Carlo code
     in NcmFitMC including gof tests. Added new functions on NcmMSet to set/get
     all models parameters.
 * Fixed bug in ncm_cfg_create_from_string which didnt indentify object strings
     with leading whitespaces.
 * Improved comment organization in keyfile generation.
 * Fixed segfault in darkenergy.
 * Fixed tests and clear functions in NcmSpline2d.
 * Fixed bugs in allocation in NcmSpline2d. Fixed wrong property name in
     NcMassFunction. Added support to validity check in NcmModel, NcmMSet and
     added   these checks in the minimization algorithms. Fixed
     NcDataClusterNCount description.
 * Added --fit-list printing all Fit options. Added get_dof to NcmData, returns
     the effective degrees of freedom of that data. Fixed bug in floating
     objects (matrix|vector) now all saved references are sunk. Fixing
     indentation. Fixed bug in Fit numerical differentiation that uses the wrong
     function to obtain   the number of free parameters. New
     ncm_(matrix|vector)_new_gsl_static functions. Fixed bug in levmar (it was
     using wrong measurement vector). Converted NcmFit object and its
     derivatives to GObject framework. Included the nlopt enum to obtain the
     correct list of algorithms. Fixing memory leaks with valgrind memcheck.
     Added ncm_message_ww to suport logging with word-wrap.
 * Removed: old nc_data_cluster_abundance.c.
 * Converted Data object (and its derived objects) to GObject framework. Converted
     Dataset and Likelihood to GObject framework. Fixed dispose method in every
     object. Added clear method for all objects. Reorganized priors in
     ncm_priors.(c|h) and nc_hicosmo_priors.(c|h). Added warning and LU
     decomposition in NcmFit when inverting the covariance matrix.
 * Static function (_ncm_fit_run_empty) was created to compute m2lnL when there is
     no free parameter. This is used to compute profile likelihood confidence
     regions.
 * Added darkenergy.1 to dist. Added test_nc_recomb. New function ncm_cmp to
     compare doubles.
 * Updating recombination object to GObject framework and several improvements.
     Updating documentation.
 * Implemented shift parameter and distance priors for WMAP7. Corrected value of
     WMAP5 shift parameter standard deviation. Message log with models and data
     used are printed for Monte Carlo runs.
 * Improving examples.
 * Testing function to print mass functions data from catalog.
 * New function ncm_model_id_by_type. Added project's URLs in configure.ac.
     Updated glib's threads usage. Fixed several documentation bugs. Fixed
     NCM_TYPE_GMSET -> NCM_TYPE_MSET. Fixed NCM_TYPE_GMSET_FUNC ->
     NCM_TYPE_MSET_FUNC Fixed make check, still missing several tests. Fixed
     lapack functions conditional compilation. Fixed constructors names
     ncm_mset_func_new_hicosmo_func(0|1). Added backward compatibility to
     compile with glib >= 2.26. Fixed backward compatibility with gsl < 1.15.
     Reworked the constants to be compatible with introspection. Working
     examples in C and Python.
 * Testing...
 * Added NUMCOSMO_HAVE_LAPACK test in ncm_lapack.h.
 * Complety rework of headers organization. Removed old lss/Makefile.am. Added a
     compatibility layer for g-ir-scan.
 * Fixed redefinition of NUMCOSMO_HAVE_INLINE.
 * Much simpler method for compiling inlined functions. Removed extra argument in
     darkenergy man page.
 * Added missing file.
 * Fixing inline functions. Added a new *_inline.c to explicity compile the
     inlined functions.
 * Fixed inline functions macros.
 * Fixed (updated) example.
 * Added backward compat for older fftw.
 * Removed INSTALL from installed docs.
 * Added format to fprintf in: ncm_cfg.c, confidence_region.c, util.c, recomb.c,
     darkenergy.c.
 * Fixed entries in darkenergy.xml.
 * Removed old catalog_parser doc
 * Fixed several bugs in conditional building. Removed asciidoc parsing to remove
     this dependency when building from repository.
 * Fixed positivity prior to use the original parameters.
 * Added test to avoid writting comments for empty entries.
 * Added positivity prior for Omega_x.
 * Removed all exit(); calls from the library.
 * Updated NEWS.
 * Including cfitsio via PKG_CHECK_MODULES.
 * Bugs: sizeof format and fgets return fixes.
 * Added tests for fit support in darkenergy.
 * Added return tests for scanf/fread/etc family functions.
 * Removed spurious header fitsio.h from print_data.c.
 * Fixed NCM_(WRITE|READ)_* macros. Added platform independent format when
     printing sizeof.
 * Added correct ifdef for cfitsio presence. Corrected printf types in ncm_cfg.c.
     Added read/write error testing in NCM_(WRITE|READ)_* macros.
 * Removed old catalog_parser.
 * Organized data object in nc_ namespace. Added support to choose data samples by
     name or nick when runnig darkenergy. Added list options in darkenergy to
     list available data options.
 * Corrected requested minimum version of gsl to use
     gsl_integration_glfixed_table.
 * Removed old INCLUDES from tools/Makefile.am. Added atlas libraries when testing
     for atlas_lapack.
 * Improved Tinker multiplicity function (critical density): for Delta_z > 3200
     the multiplicity function coefficients are now computed using the fitting
     formula given in Tinker et al. paper. Previously, when Delta_z > 3200, the
     coefficients were computed assuming Delta_z =  3200.
 * Changed precedence in darkenergy, now command line options takes precedence
     over configuration file.
 * Reworked darkenergy command line options, now each run can be saved as a ini
     file, also the options now can be specified by a .ini file which takes
     precedence over the command line options.
 * Fixed bug: when copying a likelihood the priors were copied without increasing
     their reference count. Extended NcmModel: added new property for each
     parameter describing the fit type. darkenergy not functional, changing from
     --fit-params to directly setting the fit type by setting the parameter
     property.
 * NcClusterMass... (Vanderlinde, BensonXRay, Lnnormal, Nodist) were adapted to
     NcmModel.
 * Fixed bug: HICosmo macros were modified (old: NC_MODEL; new: NC_HICOSMO...).
 * The integrations _nc_cluster_mass_vanderlinde_significance_m_p and
     _nc_cluster_mass_vanderlinde_intp were optimized.
 * The integrals to compute the probability distributions of the Vanderlinde
     mass-observable relation is being optimized.
 * Removed old INCLUES in Makefile.am. Threaded evaluation for real data in
     cluster abundance. Reorganized ncm_func_eval_threaded_loop to simply run
     the loop function when threads are disabled. Added CUBACORES=0 to
     environment in to avoid parallelization in libcuba.
 * Documentation fixes.
 * Memory leak fixed in nc_cluster_mass_benson.c and
     nc_cluster_mass_vanderlinde.c. Debug messages were removed.
 * Moved gobj_itest from bin_PROGRAMS to noinst_PROGRAMS. Added transfer full to
     the return value of ncm_reparam_ref.
 * Bug fixed: when it was set a reparametrization, all other parameters were reset
     to their default values. Now only those parameters modified by the
     reparametrization are set to the new parameter default values.
 * Tinker multiplicity function - critical: for Delta > 3200, it is set Delta =
     3200 and a warning message is provided.
 * Testing resample and montecarlo tools.
 * Fixed bug in ncm_fit_montecarlo_matrix
 * Developing new mass-observable relation.
 * Fixed bugs in reparams. Now --flat works again.
 * Testing cr algorithms.
 * Extended limits in Omega variables.
 * Added test to check if libnlopt exists. Testing nc_galaxy_acf.
 * New implementations of NcClusterMass.
 * Renamed special functions to comply with the library standards.
 * Reorganized NcMassFunction. Adjusted to correct functions names and _prepare
     function usage.
 * Fixed indentation.
 * Included g_assert in Tinker multiplicity functions (mean and critical) to
     assert that Delta <= 3200.
 * Renamed flag plane => flat.
 * Finished NcClusterMassLnnormal.
 * Finished paralelization to compute m2ln of cluster abundance. Still in test.
 * Fixed bug in function_eval lf => lfunc.
 * Updated configure.ac using autoupdate.
 * New organization of NcClusterAbundance and NcDataClusterAbundance.
 * Missing files.
 * Bug fixing and new implementations.
 * Added a simple GObject (de)serialization function set. Added nc_hicosmo_free
     function. Added tests for GObject (de)serialization. This msg is for the
     last commit.
 * Updated manual URL
 * Updated manual URL
 * Fixed documentation build.
 * New repository for savannah upload. Corrected AUTHORS and README, added
     COPYING. Erased old TODO. Version bumped to 0.8.0.
