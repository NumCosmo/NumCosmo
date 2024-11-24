CHANGELOG
----------------------

[Current]

[v0.23.0]
 * V0.23.0 (#187)

     * Updating version number.
     * Fixing documentation glitches.
 * Notebook to generate the plots for the notaknot paper (cosmology sess… (#169)

     * Notebook to generate the plots for the notaknot paper (cosmology
     session).
     * Created functions to get the spline's information (size, number of knots)
     of the halo density profile object.
     * Cleaning notebooks.
     * Updated private access.
     * Updated stubs.
     * Adding tests for new method.
     * Testing workaround for brew pkgconf.
     * Trying uninstalling pkg-config.
     * Another way to remove old pkg-config.
     
     ---------
      Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Add NcmSphereNN for finding nearest neighbors within a spherical shell. (#186)

     * Add NcmSphereNN for finding nearest neighbors within a spherical shell.
     * Uncrustify.
     * Adding support for different radius.
     * Adding missing test for dump.
     * Removed unused properties.
 * Raising error for unknown key in param_set_desc. (#184)


 * Adding support for MC analysis. (#183)

     * Adding support for MC analysis.
     * Reorganizing python modules.
 * Mass and concentration summary  (#180)

     * First draft for Halo Summary.
     * Updates for all dependent objects.
     * Fixing model update.
     * Fixing tests.
     * Updated stubs.
     * Removed unused variable.
     * Fixing leak in constructors.
     * Updating python code.
     * Tests for NcHaloMassSummary.
     * Better calibration for test_ncm_spline.
     * Removed untestable lines.
     * More testing for halo_density_profile.
     
     ---------
      Co-authored-by: Mariana Penna Lima <pennalima@gmail.com>
 * Configuring conda-incubator/setup-miniconda@v3.

 * Updating conda-incubator/setup-miniconda@v3 usage.

 * Updated conda-incubator/setup-miniconda@v3 use.

 * Adding support for version checks in numcosmo. (#179)


 * Testing more parallel tests.

 * Fix leftover merge lines.

 * Fftw config (#178)

     * Improving fftw planner control.
     * Testing ncm_cfg fftw flags manipulation.
     * Adding missing tests.
     * Fixed exception string match.
     * More functions to control fftw planner.
     * Connecting meson option to fftw planner.
 * Improving fftw planner control.

 * Configuring fftw-planner during build.

 * Using FFTW_ESTIMATE by default. Added NC_FFTW_DEFAULT_FLAGS and
     NC_FFTW_TIMELIMIT environment variables.

 * Forcing cache update.

 * Removing use-only-tar-bz2: true.

 * Adding use-only-tar-bz2: true to miniforge action.

 * Updated GHA workflow.

 * Removed old coveralls badge.

 * Twofluids update (#177)

     * Updates to TwoFluids model.
     * More tests for TwoFluids model.
     * Calibrating tests.
     * New S8 MSetFunc.
     * New S8 Gaussian prior.
     * Updated stubs.
 * Updating tests use of Vexp, fixing documentation bugs. (#176)

     * Updating tests use of Vexp.
     * Fixing documentation bugs.
 * Magnetic vexp (#175)

     * Re-parametrized magneto model.
     * Renamed and removed hardcoded paths.
     
     ---------
      Co-authored-by: EFrion <frion.emmanuel@hotmail.fr>
 * New nc_galaxy_wl_obs object  (#167)

     * Redesign of the whole cluster weak lensing analysis.
     * New unit testing for all new objects.
     * New generate command for the numcosmo app.
     
     ---------
      Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Updating to actions/upload-artifact@v4.

 * Restricting setuptools version to avoid gobject-instrospection problems.

 * Improving model interface and error handling (#174)

     * Adapting calls to methods with GError support.
     * Adding support for GError handling to methods.
     * Helper function to set GError or call g_error.
     * Bug: the same reparam passed twice would free.
     * Added GError support for methods.
     * New helper functions for GError handling.
     * Adding tests for exception raising from C.
     * Updating stubs.
     * Removed testing of first derivative when using linear interpolation.
     * Minimal tests for ncm_mset_save/load.
     * Adding missing tests for priors.
     * Improved interface for model.
     * Improved interface for mset.
 * Updated documentation of ncm_m_mass_solar. CODATA 2022.

 * Updated to latest CODATA, NIST and IAU (and others) constants. (#171)

     * Updated to latest CODATA, NIST and IAU (and others) constants.
     
     * Updated compiler versions.
     
     * Updated tests.
     
     * Updated documentation and cross-checked the CODATA, IUPAC and NIST
     values.
     
     * Fixed identation.
     
     ---------
      Co-authored-by: Mariana Penna Lima <pennalima@gmail.com>
 * Xcor cmp (#85)

     * The first version of tSZ kernel is working.
     * Fixed kernel for tSZ.
     * Added tests for tSZ.
     
     ---------
      Co-authored-by: Arthur de Souza Molina <arthur.souza.molina@gmail.com>
 
     Co-authored-by: Mariana Penna Lima <pennalima@gmail.com>
 * Xcor CCL comparisons (#168)

     * Adding the dndz notebook.
     * Adding the file with binned gaussians as dndz.
     * Adding notebooks already running the latest versions of CCL and NumCosmo.
     
     * Cleaning the notebooks outputs.
     * Updated and cleaned CCL/XCor notebooks.
     * Updated CCL precision to avoid roundoff errors.
     * Updated cmb lensing to compute correctly in curved cosmologies.
     * Cleaning notebooks.
     
     ---------
      Co-authored-by: Luigi Lucas de Carvalho Silva <luigi.lcsilva@gmail.com>
 * Implemented Integrated Sachs-Wolfe kernel. (#166)

     * Implemented Integrated Sachs-Wolfe kernel.
     * Implemented the derivative of the growth function with respect to
     redshift.
     * Cleaning notebooks.
     * Uncrustify.
     * Fixing new object nc_xcor_limber_kernel_cmb_isw.
     * Encapsulating xcor objects.
     * New tests for xcor.
     * Made tests a package to allow relative imports.
     * Organized fixtures in a different file.
     * Using Stefan-Boltzmann constant.
     * Adding guard when eval inverse distance.
     * Fixed limits determination.
     * Setting more updated constants to CCL.
     * More high precision parameters.
     * New tests comparing with CCL.
     * Fix l dependent factor.
     * Fixed upper redshift for integration.
     * Removed time limit for some tests.
     * Fixed docstrings.
     * Fixing import.
     * Fixed power-spectrum derivative.
     * More tests for kernels.
     * Ignoring untestable lines.
     * Updated pylint python version to 3.10.
     * Added method to NcDistance to compute distance from z1 to z2 without
     cancellation.
     * New Cosmology python class to hold NumCosmo's cosmology and computation
     tools.
     * Renamed fixture files and reorganizing fixtures.
     * More fixtures.
     * Adding types and using Cosmology to hold NumCosmo outputs.
     * Improving weak-lensing kernel computation.
     * Tests with reorganized fixtures.
     * Computing the weak-lensing kernel in a efficient way.
     * Adding tests for weak-lensing kernel.
     * Adding tests for galaxy counts kernel.
     * Reorganized all fixture and tests. 
     * Improved magnification bias computation.
     * Testing different bias interpolations.
     * Fixed bug in gsl spline serialization.
     * Testing comoving distance small difference.
     * Testing GSL set/get type features.
     * Testing galaxy kernel methods.
     * Testing CLASS powspec derivative.
     * Removed old untested alternative integration methods in xcor.
     
     ---------
      Co-authored-by: Mariana Penna Lima <pennalima@gmail.com>
 * Galaxy WL reformulation (#93)

     * New NcGalaxyWL object and related prototypes (nc_galaxy_sd_position,
     nc_galaxy_sd_z_proxy, nc_galaxy_sd_shape, etc.)
     * Changed naming scheme from GSD to GalaxySD
     * Fixed typos, documentation, and copyright notices
     * Added new observation matrix property and related methods (eval_m2lnP,
     nc_galaxy_wl_likelihood_prepare, etc.)
     * Registered new objects and prototypes (nc_galaxy_sd_position_flat,
     nc_galaxy_sd_z_proxy_gauss, nc_galaxy_sd_shape_gauss, etc.)
     * Added unit tests for new objects and prototypes
     * Improved galaxy weak lensing likelihood implementation
     * Optimized Monte Carlo integration and sampling methods
     * Added support for integration to nc_galaxy_sd objects and
     nc_data_cluster_wl
     * Refactored Galaxy objects to simplify properties
     * Added leave-one-out cross validation method
     * Enabled integral parallelization
     * Added normalization factors and new parameters (z_cluster, true_z_min)
     * Fixed bugs in integration and KDE evaluation methods
     * Added tests for integration and KDE comparison
     * Updated weak lensing cluster mass fitting interface
     * Refactored and cleaned up code
     
     ---------
      Co-authored-by: Caio Lima de Oliveira <caiolimadeoliveira@proton.me>
 
     Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Magnetic Fields in Vexp cosmology (#153)

     * Added spectrum computation to NcHIPertEM and updated notebook.
     * New notebook for vexp_bounce_adiabatic.
     * Updates on vexp_bounce_adiabatic.

[v0.22.0]
 * New version v0.22.0

 * Mix experiments options (#164)

     * Fixing doc-strings.
     * Fixing default lnk0 value.
     * Reorganizing and documenting SNIa objects.
     * Reorganized SNIa serialization.
     * Changed serialization of NULL string to be NULL not an empty string.
     * Added option of using SNIa data to Planck experiments.
     * Adding best fit extraction from catalog.
     * Adding support for lensing plik.
     * Fixed output bug and added support for rasterizing corner plots.
 * CCL background power (#162)

     * Updated notebooks to CCL version 3.0.
     * Cleaning notebooks.
     * Updated numcosmo imports.
     * Updated CCL usage.
     * Updated Colossus usage.
     * Fixed wrong format in print.
     * Updated notebooks: CCL and NumCosmo cross-check  - background and power
     modules.
     * Implemented the comoving volume element function in nc_distance.c. It is
     valid for any curvature.
     * Implemented Sigma critical(infinity) in nc_distance.c. The respective
     funtions in nc_wl_surface_mass_density.c calls it now.
     * nc_halo_mass_function_dv_dzdomega: updated - it calls
     nc_distance_comoving_volume_element.
     * Indenting nc_distance.c
     * Updated test due to fixed volume computation for curved cosmologies.
     * Testing numcosmo generate planck18.
     * Making healpy optional.
     * Minor fixes, added test for nc_halo_mass_function_dv_dzdomega.
     
     ---------
      Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Two Fluids primordial model (#160)

     * New notebook where we compute perturbations for the two fluid model.
     * Updates on two fluids code and notebook.
     * New HIPrimTwoFluid model.
     * New notebook with HIPrimTwoFluid calibration code.
     * Updated CMB experiments to use omega parametrization.
     * Including two_point model in the app experiment generator.
     * Registering two_fluids model.
     * Extending calibration to w = 1.0e-7.
     * Added option to skip catalog validation when continuing a MCMC.
     * Better logic to use skip_check.
     * Fixed bug in calc_param_ensemble_evol.
     * Updated docstring and added more catalog related tools.
     * Save plot option.
     * Calibrating model.
     * mcmc_file is now non-optional.
     * Fixed default value.
     * Adding log_file option to numcosmo.
     * Tests for nc_hiprim_two_fluids.
     * Adding tests for hipert_two_fluids.
     * Testing new ensemble evol methods.
     * Adding tests for output.
     * More tests for ncm_fit_esmcmc.
     * More testing for hipert_two_fluids.
     * Increasing number of point in fftlog jljm.
     * Fixed old edge case.
     * Added extra safe-guard.
     * Removing coveralls upload.
     
     ---------
      Co-authored-by: Luiz Demétrio <demetrio.luizfelipe.fis@gmail.com>
 * Support for multiple corner plots. (#161)


 * Added support for parameter filtering in numcosmo app. (#159)

     * Added support for parameter filtering in numcosmo app.
     
     * More fftlog calibrate.
     
     * Support for different typer behaviors.
 * Adding support for 1d distributions. (#158)


 * Updating uncrustify configuration. (#157)

     * Updating uncrustify configuration.
     * Applying to all files in numcosmo/ numcosmo/math numcosmo/model
     numcosmo/xcor
     * Adding .clang-format for future reference.
     * Testing indentation on CI.
     * Removed uncrustify validation due to outdated uncrustify in conda-forge.
     * Testing another strategy to check indent on CI.
     * Updating package for arb.
     * Added break-system-packages to pip.
 * Calibrating fftlog tests number of knots.

 * Fixing spline constructors to return the right type. (#156)


 * Removing glib version restriction (#155)

     * Removing version restriction.
     * Updated stubs.
     * Using brew to install pygobject instead of pip.
 * NcHIPert Reformulation (#95)

     * Removed old singularity code from CSQ1D. Working version of
     vexp_bounce.ipynb.
     * Removing old singularity interface (deprecated).
     * Updated perturbations to use CSQ1D instead of the old and deprecated
     HOAA.
     * Simplified CSQ1D interface. Now the subclasses are responsible for any
     extra parameter of the system.
     * Redesigned CSQ1D to better organize the output in terms of different
     parametrizations and frames.
     * Updated nc_de_cont, vacuum_study.ipynb and test_py_csq1d.py to use the
     new design.
     * Updated stubs.
     * New method to compute the state at a time and frame.
     * Updated vacuum_study_adiabatic.
     * Updated vacuum_study.
     * Updated unit testing.
     * Added tests for non adiabatic vacuum and its propagation using ODE and
     the propagator.
     * More tests for NcmCSQ1D and NcmCSQ1DState.
     * Cleaning notebooks.
     * Added electromagnetic constants to NcmC. Added unit testing.
     * New NcHIPertEM for free electromagnetic field computation.
     * Improved interface for perturbation objects (work in progress).
     * Updated vexp_bounce.ipynb to use the new interface.
     * Updating primordial_perturbations/magnetic_dust_bounce.ipynb.
     
     ---------
      Co-authored-by: EFrion <frion.emmanuel@hotmail.fr>
 Co-authored-by:
     Eduardo Barroso <eduardojsbarroso@gmail.com>
 * Adding Bayesian evidence support for numcosmo app. (#152)

     * Adding Bayesian evidence support for numcosmo app.
     * Removed black version restriction.
 * Sample variance (#107)

     * SSC comparison with lacasa
     * Added gauss data object
     * Added ncounts data object
     * SSC gauss
     * adding bias crosscheck
     * Fisher matrix notebook
     * fisher matrix for cluster ncounts with super sample covariance with
     numcosmo internal object
     * Removing old ncounts gauss object
     * Fixed example.
     * some files to run PySSC with NumCosmo
     * Cleaning notebooks.
     * Fixing python linter problems.
     * Generating stubs.
     * More tests.
     
     ---------
      Co-authored-by: Henrique Lettieri <henrique.cnl@hotmail.com>
 * Updating codecov to v4. (#151)

     * Updating codecov to v4.
     * Skipping openmpi on macos.
 * Updating requirement versions. (#150)

     * Updating requirement versions.
     * Updating python version in CI.
 * Fixed instrospection error for gobject-instrospection >= 1.80. (#149)

     * Fixed instrospection error for gobject-instrospection >= 1.80.
     * Trying to fix bug in mpi run in MacOS by using brew's openmpi.
     * Adding support to control mpi use in meson.
     * Disabling mpi in macos brew build (bug in the MPI lib).
 * Cmb parametrization (#148)

     * Change the default parametrization for CMB experiments.
     
     * Removing w from parameters to set.
 * Planck data analysis reorganization (#145)

     * Adding stub generation script to gitignore.
     * Fixed planck code linking to allow dlsym to work.
     * Reorganized planck parameters to match baseline defaults.
     * Added option to generate planck 18 experiments.
     * Added automatic download of planck data.
     * Removing protected from symbol visibility in plc and using link_whole to
     allow plc dlsym usage.
     * Set default parameters for experiments.
 * Adding stub generation script to gitignore.


[v0.21.2]
 * Updated stubs.

 * New minor release v0.21.2

 * Adding MPICH support. (#144)

     * Adding MPICH support.
     * Descreasing allowed m2lnL variance for exploration phase (leading to
     overflow during matrix inversion for large dimensions).
     * Testing MPI support.
     * Better names for CI.
 * NumCosmo product file (#143)

     * Introducing the product-file options.
     * Updated Halofit to return linear Pk when all required scales are linear.
     * Support for exploration phase in APES.
     * Adding calibrate option to numcosmo app.
     * Reorganized numcosmo app.
     * Removed calibrate_apes tool.
     * Testing power-spectra with/without halofit.
     * More tests for linear universe.
     * Added tests for APES exploration.
     * Testing APES MPI.
     * Removed unused code and optimizing tests.

[v0.21.1]
 * Moving release v0.21.1.

 * More options to the conversion tool from-cosmosis. (#141)

     * mute-cosmosis makes cosmosis do not print info messages.
     * reltol sets the tolerance for NumCosmo underlying cosmology.
 * New minor release v0.21.1.

 * Updates and tests for NumCosmo app (#140)

     * Fixed minor bugs in levmar.
     * Added serialization for strv.
     * Tests for NumCosmo app.
     * Adding strv to serialization tests.
     * Making all test files in a tmp dir.
     * Conditional tests for numcosmo app (depends on typer and rich).
     * Added optional reqs to pyproject.
     * Fixed conditional compilation of yaml serialization methods.
     * Minimal tests for the complete functionality of NumCosmo app.
     * Finished numcosmo analyze.
     * Fixing flake8 issues.
     * NumCosmo analyze behaves well for number of iterations < 10.
     * Improved message for test failing on all parameters.
 * New bug fix release.

 * Fixed mypy issues.

 * Ran black.

 * Fixed cosmosis required parameters issue due to returning iterator. Fixed
     restart issue on numcosmo run fit.


[v0.21.0]
 * New release v0.21.0

 * numcosmo command line tool (#137)

     Introduced a new command line tool for NumCosmo (experimental):
     
     * `numcosmo from-cosmosis` converts a cosmosis ini file to NumCosmo yaml
     format
     * `numcosmo run fit ` computes the best-fit for an experiment (NumCosmo's
     analysis)
     * `numcosmo run test` test an experiment
     * `numcosmo run fisher` computes a fisher matrix
     * `numcosmo run fisher-bias` computes a fisher matrix
     * `numcosmo run theory-vector` computes the theory vector
     * `numcosmo run mcmc apes` computes the MCMC sampling of the experiment
     using APES
     * Fixed serialization for require_nonlinear_pk.
     * Now NcmFit calls m2lnL just once if no parameters are free.
     * New methods to NcmData and NcmDataset to check if mean_vector is
     available.
     * Updated stubs and requiring black < 24 due to difference in formatting.
     * Unit tests for the newly added code.
 * Variant dictionary support  (#135)

     * New NcmVarDict boxed object describing a dict of str keys and basic types
     values.
     * Added unit testing
     * Support for serialization of NcmVarDict
     * Improved valgrind suppresion file
     * Changing Variant type of object to tuples.
     * Finished update of Object variant type. Added tests for data files.
     * Finished support for VarDict as object properties.
     * Updated conda environment file to use openblas compatible with openmp.
     * Fixing problem with fft wisdow when MKL is being used.
     * Updated priors to use named parameters. Improved Model and MSet objects
     use of full parameter names.
     * Improving reports to codecov.
     * Using conda build for coverage.
     * Fixed wrong signness comparison and coverage .
     * Testing lcov 1.16 options.
     * Removing timeout for coverage tests.
     * Removing external codes from coverage.
     * Disabling documentation build in CI.
     * Removing G_DECLARE_ from coverage.
     * Ignoring G_DEFINE_ in coverage.
     * Using lcov for coveralls.
     * Disabling branch detection.
     * Extra tests for NcmMSet and adding tests back to coverage.
 * Support for object dictionaries, NcmObjDictStr and NcmObjDictInt. (#134)

     * Support for object dictionaries, NcmObjDictStr and NcmObjDictInt.
     * Unit testing
     * Updated stubs
     * Fixed leaks
 * Better python executable finding.

 * Added GSL as a dependency for libmisc (internal library). (#133)

     * Added GSL as a dependency for libmisc (internal library).
     * More missing deps for libmisc.
     * Removed unnecessary includes omp and added missing deps to class.
     * Improving a few includes.

[v0.20.0]
 * Updated stubs.

 * New minor version v0.20.0

 * Support for computing fisher bias vector (#132)

     * Added support for computing fisher bias vector and corresponding unit
     tests.
     * Increased timeout for GaussCov and conditional testing in likelihood
     ratio.
     * Interface for fisher bias computation and tests.
     * Better calibration for ncm_fit tests.
     * Added retry in the hessian computation.
 * Improving tests for NumCosmoMath (#131)

     * Improving tests for NcmFftlog.
     * Testing q=0.5 case for j_l.
     * Added support to make meson use tap protocol for testing.
     * Using g_assert_true instead of g_assert in tests.
     * More tests for NcmFit.
     * Tests for likelihood ratio tests, removed old and unused ncm_fit methods.
     
     * Testing MPI.
     * Adding libopenmpi-dev to ubuntu installations.
     * Adding MPI tests only on supported envs.
     * Created the USE_NCM_MPI flag to use when compiling code that use
     numcosmo's MPI facilities.
     * Disabled TAP when running pytest-tap and mpi.
     * Moved ode_spline from example to unit testing.
 * Adding support for libflint arb usage. (#130)

     * Adding support for libflint arb usage
 * Adding more python based tests using external libs (astropy and scipy). (#129)

     * Adding more python based tests using external libs (astropy and scipy).
     * Adding astropy and scipy to coveralls job.
     * Testing adiabatic solutions on CSQ1D.
     * Tests for NcmDataDist1d.
     * Tests for DataDist2d.
     * Testing DataFunnel.
     * Tests for NcmDataGauss.
     * Added C test for generic garbage collection tests.
     * More tests for DataGauss.
     * Adding LCOV_EXCL_LINE non-testable lines.
     * Tests for NcmDataGaussDiag.
     * Testing bootstrap+wmean.
     * Testing GaussMix2D.
     * Testing DataPoisson and finished Poisson fisher support.
     * Testing NcmDataset.
 * Fixed package name in pyproject.toml.


[v0.19.2]
 * Including typing data into pyproject.toml. Updating changelog.

 * Fixing minor doc glitches. (#128)


 * Updated changelog.

 * Fixed project name in pyproject.toml.

 * New minor version.

 * Using pip to install python modules. (#127)

     * Using pip to install python modules.
     * Adding conda in the CI.
     * Added devel_environment.yml.
     * Update numcosmo_py.
     * Fixed flake8 issues.
 * More objects encapsulation (#126)

     * Removed NcmCalc (unfinished). Updated csq1d.
     * Encapsulated all NumCosmoMath objects.
     * Including CI testing log.
     * Adding documentation to every NumCosmoMath objects.

[v0.19.1]
 * New minor release 0.19.1

 * Updated meson to deal with cross compiling and GI building. Updated ncm.pyi.

 * Removed git ignored files related to autotools and in-source building.

 * Removed unnecessary packages.

 * Yaml implementation (#125)

     * Updated minimum glib version.
     * Complete version of the yaml serialization, including special types.
     * Updated Python stubs.
     * Adding fyaml to CI.
 * Adding fyaml to CI.

 * Updated Python stubs.

 * Complete version of the yaml serialization, including special types.

 * First version of from_yaml and to_yaml serialization. Updated minimum glib
     version.

 * New tuple boxed type (#124)

     * Added new NcmDTuple boxed objects.
     * Added serialization support and tests.
 * Objects encapsulation (#122)

     * Deleting old unnecessary files.
     * Encapsulating and documenting NcmMPI objects.
     * Basic documention for NcmMPIJob and NcmMSetTransKern.
     * Documenting NcmPrior and subclasses. Finished documentation of
     NcmPowspec.
     * More documentation for NcmModelCtrl.
     * Uncrustified NcmModel.
     * Reordered NcmModelCtrl.
     * Encapsulated and documented NcmModel. All subclasses were adapated.
     * Disabling gsl range check by default and enabled inlining in GSL.
     * Refactored models to use ncm_model_orig_param_get.
 * Removed unnecessary header inclusions to avoid propagating depedencies.

 * Fixing warnings in conda build. (#121)


 * Mypy to ignore python scripts inside meson builds.

 * New test for simple vector set/get.

 * Removed old files.


[v0.19.0]
 * Release v0.19.0.

 * Testing before adding warn supp. Testing for isfinite declaration.

 * Adding cfitstio to plc.

 * Addind examples to installation.

 * Added libdl dependency to plc.

 * Added GSL blas definition to avoid double typedefs.

 * Updated changelog.

 * Updated stubs.

 * Moving to meson (#120)

     * Removing trailing ; from G_DECLARE and G_DEFINE.
     * Removing wrong return statement in void functions.
     * Moving NcXcorKinetic boxed object from Xcor to XcorLimberKernel.
     * Removed vector_dot method from inlined methods.
     * Uncrustified.
     * Fixing sign and unsigned mixing.
     * Fixed bug where we used height instead of nnodes.
     * Adding blas header (we are removing an overall cblas include).
     * Adding missing initialization.
     * Fixing ifdef for optional FFTW.
     * Fixing object to function pointer transformation.
     * Better control for optional dependencies.
     * Fixing fallthrough warnings.
     * Removing NUMCOSMO_ prefix from internal macros.
     * Encapsulating NcmStatsVec
     * Moving NUMCOSMO_HAVE_CFITSIO to HAVE_CFITSIO and removing unnecessary
     macros and includes.
     * Organizing ncm_fit headers.
     * Adding support for meson build system.
     * Added support for generating vala binding.
     * Removed outdated factory functions.
     * Removed outdated GSL versions (new requires 2.4).
     * Introduced HAVE_MPI.
     * Updated stubs. Installing python modules.
     * Updated factory functions.
     * Using a fixed version of references.xml.
     * Removing autotools.
     * Removed makefile leftover.
     * Removed makefiles.
     * Removed old autotools files.
     * Comment explaning pkg generation.
     * Including gsl in numcosmo.pc.
     * Forcing interface to avoid broken blas.
     * Try calling gcovr directly.
 * Minor improvements. (#119)

     - Adding new sampler to the notebook with sampler comparions.
     - Added support for weight samples.
     - Fixed bug in ncm_stats_vec.c (first element with zero weight resulted in
     a nan).
 * Added correct prefix for NCM_FIT_GRAD.

 * Kde loocv (#118)

     * Testing new CV types.
     * Removed wrong break.
     * Fixed allocation error.
     * Testing stopping criteria for integration.
     * Testing amise integral.
     * Removing debug printing.
     * Sampling using antithetic variates to improve convergence.
     * Optimizing MC integration for LOO.
 * Creating tests and documentation for n-dimensional integration object (#108)

     * Created initial test
     * Introduced macros to simplify the creation of IntegralND subclasses and
     resolved unit test issues.
     * Improved tests
     * Improved documentation
     * Typo fix
     
     ---------
      Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Improving documentation and encapsulating objects (#116)

     * Now MPI jobs do not require setting nthreads.
     * Using per rank fftw wisdom.
     * Documenting NcmFitESMCMC.
     * More documentation for BAO data objects.
     * More docs for BAO and NcDataClusterNCount.
     * More documentation for Fit* objects and Monte Carlo analysis object.
     * Updated and encapsulated NcmFit and all depending objects.
     * Moving examples to tests.
     * Fixed test name and Makefiles.
     * Moving example_diff.py -> ../tests/test_py_diff.py.
     * Encapsulated FitState and updated all required objects.
     * More tests for NcmFitState.
     * Fixed possible negative precision.
     * Excluing impossible lines from coverage.
     * Included levmar in unit testing.
     * Fixed bug in levmar and gsl_mm.
     * Testing serialization of NcmFit.
     * Removed old analytical derivative support.
     * Removed last link on the analytical derivative support.
     * Testing fit restart.
     * More testing for NcmFit and documentation for NcmMSet.
     * Fixed bug in fit_levmar, more testing for sub fits.
     * Fixed bug in accurate grad (missing matrix transposition).
     * More testing, fixed sub fit testing.
     * Added missing reset states in _gsl_mms.
     * Improved sub-vector manipulation and added necessary tests.
     * Added equality constraint tests.
     * Adding testing for inequality constraints.
     * Testing constraints serialization.
     * Adjusted diff to use a better estimate when errors cannot be estimated.
     * Improving diff computation when never converging.
     * Organized fisher.
     * Testing Fisher matrix and covariances.
     * Removing old unused methods.
     * Removed option to print fisher matrix to file.
 * Now MPI jobs do not require setting nthreads. (#115)

     * Now MPI jobs do not require setting nthreads.
     * Using per rank fftw wisdom.
 * 109 example describing 3d correlation (#113)

     * New tutorial.
     * New version of FFT code with changes in scale and in the xi function.
     Translated from portuguese to english.
     * Added the generalized Fourier Transform function and its inverse in the
     calculations to compare with xi(r,z).
     * Working on the FFT
     
     ---------
      Co-authored-by: Maria Vitoria Lazarin <mvitoria.lazarin@gmail.com>
 * Added pocoMC to rosenbrock_simple.ipynb. (#112)


 * Updated rosenbrock_simple.ipynb.


[v0.18.2]
 * New version v0.18.2.

 * Improving stubs.


[v0.18.1]
 * New minor version v0.18.1.

 * Missing files for python typing

 * Updated changelog.


[v0.18.0]
 * Updated changelog.

 * New minor release 0.18.0.

 * Implementing n-dimensional integration object (#106)

     * Initiating implementation
     * Adding IntegralND to docs.
     * Changed error for testing
     * Adding tests to integralnd.
     * removed unnecessary variables
     * Adding numpy to testing.
     * Renamed object and improved coverage.
     * Added documentation.
     
     ---------
      Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Create SECURITY.md
 * Create CONTRIBUTING.md
 * Create CODE_OF_CONDUCT.md
 * Update issue templates
 * Update bug_report.md
 * Update issue templates (#104)

     * Update issue templates
 * 102 add notebook for gauss constraint tests (#103)

     * New gauss_constraint_mcmc.ipynb notebook. 
     * Minor improvements and tweaks on numcosmo_py.
     * Adding xcdm experiment to the example_apes.py. 
     * Updated default oversmooth to match new VKDE parametrization.
     * Fixed typos.
     * Added burnin option in getdist helper function.
     * Cleaning notebooks.
     * Black on notebooks.
 * Added TMVN sampler. (#101)

     * Added TMVN sampler.
 * Update README.md
 * Update README.md (#99)


 * Updated python interface. (#98)

     * Updated interface and generated stubs for mypy. All examples were
     updated.
     * Adding tests to new methods.
     * Improving tests.
 * Fixed bug that resets the values of use_threads in APES.

 * Several improvements on APES and others. (#92)

     * Several improvements on APES and others.
     
     - Added linters to python code on NumCosmo.
     - Updated severeal examples, more to come.
     - Removed old and/or incomplete examples.
     - Added support for for BLIS (BLAS like library).
     - This commit has timing logs on the APES code (it will be removed soon).
     - Complete parallelization using OpenMP.
     - Combininig different parallelizations OpenBLAS, VKDE, Interpolation, etc.
     
     - Better documentation and error messages. Removed debug and timing prints.
     
     - Updated examples. Organized walker's thread usage.
     - Removed old and unused code.
     - Minimal tests to ncm_cfg.
     - Calibrating lcov exclusions.
 * Removed printf from test.

 * Improving parallelization for APES.

     - Made several methods reentrant (kdtree, ncm_stats_dist_vkde
     ncm_stats_dist_kde).
     - Added parallelization to KDE preparation (VKDE and KDE).
     - Updated multi-thread model to use OpenMP in NcmFitESMCMC.
     - Fixed the associated tests.

 * Fixed setting of max_ess.

 * Added conditional compilation of internal function.

 * Added support in ncm_mset_catalog and mcat_analyze to compute acceptance ratio.

 * Many minor improvements.

     * Added one more numcosmo_py experiment: gauss_constraint.
     * Added source dir tools PATH numcosmo_export.sh to allow non-installed
     version to find python tools.
     * Working in a new primordial model.
     * Added option to mcat_calibrate_apes to whether to plot the calibrated
     results.
     * Minor identation and documentation tweaks.

 * Fixed scripts shebang.

 * New version 0.17.0.

     Minor documentation glitch fixes.


[v0.17.0]
 * New version 0.17.0.

     Minor documentation glitch fixes.

 * New experimental python interface for sampling. New sampler comparisons. (#90)

     * New experimental python interface for sampling. New sampler comparisons.
     * More details on the rosenbrock_simple.ipynb notebook. Ran black.
 * Encapsulating objects (#72)

     * Encapsulated ncm_c, ncm_csq1d, ncm_data and ncm_data_dist1d. Uncrustify
     all ncm_data_* files.
     * Encapsulated ncm_bootstrap and ncm_data_dist2d.
     * Encapsulated and indented ABC objects.
     * Encapsulating ncm_data_funnel, ncm_data_gauss_cov and
     ncm_data_gauss_cov_mvnd.
     * Encapsulated ncm_data_gauss and ncm_data_gauss_diag.
     * Encapsulated ncm_data_poisson. Added documentation to ncm_fit_state.
     * Encapsulated Fftlog objects. Added documentation.
     * Added documentation and fixed flake8 and mypy issues on scripts.
     * Fixed tests to adapt to the encapsulated objects.
     * Uncrustifying files.
 * New features (#81)

     * New support for pytest unit testing.
     * Removed old man pages generations.
     * Added optional support for libfyaml.
     * Added rb_knn_list to documentation ignore list.
     * Cleaned and organized apes tests notebooks.
     * New apes_tests/xcdm_nopert.ipynb.
     * Fixed documentation glitches.
     * Started implementation of Serialization yaml backend.
     * Fixed AC_DEFINE for libfyaml.
     * Finished first version of the YAML serialization.
     * Fixed NcmModel unit tests to handle double sub-indices.
     * Added support for installing numcosmo_py. Added missing docs to
     numcosmo-docs. Renamed mcat related scripts. New interface in
     nc_hicosmo_Vexp.
     * Added command to update apt index before installing prereqs.
     * Added xcdm experiment with cosmological dataset without perturbation
     dependent data.
     * Added support for pytest in CI.
     * Testing example compilation without installing the library.
     * Fixed libfyaml dependent static function in ncm_serialize.c. Added link
     to numcosmo library in test of non-installed library.
     * Better test calibration.
     * Improved doc in ncm_cfg.
     * Improved wisdown handling in ncm_cfg.
     * Improved coverage tweaking.
     * Do not CI all branches.
     * Testing python example running without installing.
 * WL binned likelihood object (#77)

     * New method for likelihood utilizing KDE
     * Added necessary support to compute the galaxy wl likelihood using KDE.
     Fixed leak in nc_galaxy_wl.
     * Cut galaxies with e,g < 0, e,g > 0.05
     * First working version of the KDE likelihood for WL.
     * Fixed behavior for g_i < 0
     * Better control of the border in the galaxy KDE.
     * Renamed reduced shear to ellipticity to focus on the true weak lensing
     observable. Created new object class nc_galaxy_wl_ellipticity_kde and moved
     calculations from nc_galxy_wl to nc_galaxy_wl_ellipticity_kde. Introduced
     new method nc_galaxy_wl_dist_m2lnP_initial_prep. Edited Makefile.am,
     Makefile.in, ncm_cfg.c and numcosmo.h to accomodate changes. Changes have
     made calculations slower but results seem to be in line with previous
     version.
     * Freeing s_kde and g_vec on _nc_galaxy_wl_dist_initial_prep seems to have
     fixed memory leak (?) issue on last commit
     * Fix's fix. Freeing the memory allocated to s_kde and g_vec is what was
     initially causing the segmentation fault.
     * Fixed indentation.
     * Removed comments and added doc.
     * Reorganizing internal of nc_galaxy_wl_ellipticity_kde. Indentation on
     other related objects.
     * Fixed email and unnecessary variable.
     * Fixed object names on docs.
     * Fixed numcosmo/lss/nc_galaxy_wl_ellipticity_kde.c bug. When resetting
     self->kde, epdf_bw_type stayed as FIXED instead of being recast as RoT.
     * First attempt at creating binned object for wl likelihood. Currently
     working on creating a NcmObjArray with the galaxy data belonging to each
     bin. Probably (certainly) very buggy.
     * First version of unit test for nc_galaxy_wl_ellipticity_kde.
     * Added lss/nc_galaxy_wl_ellipticity_kde to numcosmo.h and math/ncm_cfg.c
     * Removed old code necessary for debugging.
     * Fixed copyright notice, fixed reset of self->e_vec, added peek_kde and
     peek_e_vec methods.
     * Fixed file name when adding tests, fixed nc_distance_comoving error by
     preparing dist object, reduced number of tests.
     * Registered new object NcGalaxyWLEllipticityBinned
     * Changed approach to binning, fixed set_bin_obs, object is now compiling
     * Finished binned weak lensing object prototype
     * Started implementation of unit test for wl_ellipticity_binned, setting
     gebin->binobs as NcmObjArray on initialization.
     * Fixed binning behaviour. Bug was caused by casting bin limits as an int
     instead of a float. Tests are all passing now.
     * Update nc_galaxy_wl.c
     * Fixed copyright notices
     
     ---------
      Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 
     Co-authored-by: Caio Lima de Oliveira <caiolimadeoliveira@proton.me>
 * Added check to see if the python interface is available.

 * Improved tests.

 * New method for likelihood utilizing KDE (#65)

     * New method for likelihood utilizing KDE
     * Added necessary support to compute the galaxy wl likelihood using KDE.
     Fixed leak in nc_galaxy_wl.
     * Cut galaxies with e,g < 0, e,g > 0.05
     * First working version of the KDE likelihood for WL.
     * Fixed behavior for g_i < 0
     * Better control of the border in the galaxy KDE.
     * Renamed reduced shear to ellipticity to focus on the true weak lensing
     observable. Created new object class nc_galaxy_wl_ellipticity_kde and moved
     calculations from nc_galxy_wl to nc_galaxy_wl_ellipticity_kde. Introduced
     new method nc_galaxy_wl_dist_m2lnP_initial_prep. Edited Makefile.am,
     Makefile.in, ncm_cfg.c and numcosmo.h to accommodate changes. Changes have
     made calculations slower but results seem to be in line with previous
     version.
     * Freeing s_kde and g_vec on _nc_galaxy_wl_dist_initial_prep seems to have
     fixed memory leak (?) issue on last commit
     * Fix's fix. Freeing the memory allocated to s_kde and g_vec is what was
     initially causing the segmentation fault.
     * Fixed indentation.
     * Removed comments and added doc.
     * Reorganizing internal of nc_galaxy_wl_ellipticity_kde. Indentation on
     other related objects.
     * Fixed email and unnecessary variable.
     * Fixed object names on docs.
     * Fixed numcosmo/lss/nc_galaxy_wl_ellipticity_kde.c bug. When resetting
     self->kde, epdf_bw_type stayed as FIXED instead of being recast as RoT.
     * First version of unit test for nc_galaxy_wl_ellipticity_kde.
     * Added lss/nc_galaxy_wl_ellipticity_kde to numcosmo.h and math/ncm_cfg.c
     * Removed old code necessary for debugging.
     * Fixed copyright notice, fixed reset of self->e_vec, added peek_kde and
     peek_e_vec methods.
     * Fixed file name when adding tests, fixed nc_distance_comoving error by
     preparing dist object, reduced number of tests.
     * Fixed indentation
     
     ---------
      Authored-by: Caio Lima de Oliveira <caiolimadeoliveira@gmail.com>
 * Reordering -I to include first internal sub-packages.

 * Added conditional use of g_tree_remove_all. Removed setting of all threads to
     one. Reintroduced non fatal assertions in test_ncm_stats_dist.

 * Improved ax_code_coverage.m4 to work with newer lcov versions.

 * Organized m4 files and fixed lcov issues.

 * Fixed a few lcov issues.

 * Incresead number of points when testing StatsDist with rubust-diag.

 * Changed divisions to multiplications.

 * Fix bug in AR fitting when only two elements were available.

 * Improved tests, added test to robust covariance computation.

 * Modified vkde to use block triangular system solver.

 * Finished the refactor of kdtree to use a red-black tree and prune impossible
     branches. Significant increase in speed!

 * Removed old tree. Finishing prunning.

 * New red-black BT to improve kdtree. Added support for prunning kdtree during
     search.

 * 75 organizing python modules (#76)

     * First reorganization and cleaning.
     * New directory structure.
     * Refactoring scripts to satisfy linters.
     * Organizing notebooks.
     * Reorganizing notebooks, new MCMC tests with gaussian mixture models.
     * Organizing python scripts and support in NumCosmoPy.
     * Refactored mcat_calibrate_apes.py.
     * New interpolation object.
     * Improving getdist (added asinh filter).
     * Added check for unclean notebooks.
     * Cleaning notebooks.
     * Checking notebooks first in GHA.
     * Normalizing filenames and organizing Makefile.am.
     * Using lower-case names for notebooks.
     * Updated Makefile.am.
 * Updated kdtree and directories in notebooks/Makefile.am.

 * Reorganizing notebooks.

 * Added missing cell.

 * Fixed missing properties (unused). autogen.

 * Notebooks massfunc (#74)

     * Time tests and tests of the Bocquet multiplicity function were added in
     NC_CCL_mass_function notebook.
     * Info about critical Delta was added.
     * Notebook NC_CCL_mass_function was organized.
     * Included data used in Penna-Lima et al. (2017) - Planck-CLASH clusters.
     * Updated plcl script.
     
     ---------
      Co-authored-by: Cinthia Lima <cinthia.n.lima@hotmail.com>
 Co-authored-by:
     Mariana Penna Lima <pennalima@gmail.com>
 * Halo bias tests (#73)

     * Fixed Function Type Definition
     
     * Finished Halo Bias Tests
 * Mean bias (#61)

     * posterior volume
     * Mean halo bias
     * mean bias unbinned plot
     * mean bias binned case
     * Data DC2
     * Redmapper data ncounts with richness
     * Data preparation with richness
     * skysim mcmcm file
     * DC2 tests
     * checking watson on nersc
     * moved set/get Delta methods to parent class. Merged Bias type and Bias
     Func.
     * Updated all bias_func and bias_type to bias.
     * Nodist MCMC analyses on Mock catalog
     * Removed data files from repository.
     * Added nc_halo_mass_function_peek_multiplicity_function to access the
     NcMultiplicityFunc in NcHaloMassFunction.
     * Tinker Bias Delta correction.
     * Including the NcHICosmo in bias_eval function.
     * Refactoring old Bias code to match new design.
     * Tinker and ST_spher correction on new and new_full functions.
     * Updated autotools. Added numcosmo-valgrind.supp for valgrind memcheck.
     Fixed optimization bug in ncm_sphere_map.
     * Uncrustify all halo bias objects.
     * Implementation on the volume element on the bias integrand
     * Improving integrand for bias computation.
     * Bias as function of mass
     * Added minimal documentation to bias objects.
     * Finished Documentations
      Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 
     Co-authored-by: Eduardo Barroso <eduardojsbarroso@gmail.com>
 
     Co-authored-by: root <root@eduardo>
 * Added interface to generate models using an array of NcmSParams

 * Added two missing files to the releases.

 * New minor release v0.16.0.


[v0.16.0]
 * New minor release v0.16.0.

 * 40 numcosmo unit test coverage (#68)

     * Fixed e-mail addresses.
     
     * Updated file/object names. Fixed warnings in libcuba.
     
     * Fixed symbol
     
     * Fixed minor documentation issues.
 * Added new method to set model parameter fit types to their default values.

 * Missing semicolon.

 * Minor changes on NcDistance initialization order.

 * Updated gcc version for macos ci.

 * Better debug messages in GHA.

 * Added more robust testing for power-spectra.

 * Updating e-mails.

 * Updated e-mail in copyright notices.

 * uncrustify code.

 * Minor fixes in documentation. Finished coverage and tests for special
     functions. Removed old code.

 * Added support for namespace search in ncm_mset_func_list. Added plot_corner
     helper script.

 * More cleaning and adding more files to .gitignore.

 * Cleaning autotools files and old unused tools. (#67)


 * Removed old and unused code.

 * Uncrustify sources.

 * Reorganized all ncm_spline2d objects and improved unit testing and coverage.

 * Uncrustify ncm_spline2d_bicubic.

 * Improved coverage of NcmDiff.

 * Uncrustify ncm_diff.c.

 * 60 statsdist1d error (#63)

     * Removed old code causing a bug in ncm_stats_dist1d_epdf_reset. Updated
     autotools. Included more tests in test_ncm_stats_dist1d_epdf.
     
     * Uncrustify sources.
     
     * Adding more tests to test_ncm_stats_dist1d_epdf.
     
     * Fixed retry leak and decreased max_retries.
     
     * More updates in coverage support.
     
     * Removed wrong macro in g_test_trap_subprocess
 * Removed old code causing a bug in ncm_stats_dist1d_epdf_reset.  (#62)

     * Removed old code causing a bug in ncm_stats_dist1d_epdf_reset. Updated
     autotools. Included more tests in test_ncm_stats_dist1d_epdf.
     
     * Uncrustify sources.
     
     * Adding more tests to test_ncm_stats_dist1d_epdf.
     
     * Fixed retry leak and decreased max_retries.
     
     * More updates in coverage support.
     
     * Removed wrong macro in g_test_trap_subprocess
 * Fixed unimportant warnings in class.

 * Removed debug message.

 * Minor fixes in twofluids framework. Updating StatsDist to use only a fraction
     of the sample when computing the bandwidth using a split cross-validation.

 * Multiplicity watson (#59)

     * The files of Watson et al. multiplicity function (.c and .h) has been
     created and included in the files makefile.am, math/ncm_cfg.c and
     numcosmo.h
     
     * watson multiplicity function updated
     
     * complementing the watson et al. multiplicity function
     
     * Fixed bugs and added tests.
     
     * Fixed bugs.
     
     * uptade libtool files
     
     * removed extra files
     
     * missing files
     
     * multiplicity_watson_install
      Co-authored-by: Cinthia Lima <cinthia.n.lima@hotmail.com>
 Co-authored-by:
     Henrique Lettieri <henrique.cnl@hotmail.com>
 * Removed ckern algo.

 * Debug version, do not use it. Version containing the constant kernel option in
     NcmStatsDist.

 * Several minor improvements.

     - Removed configure call from autogen.sh.
     - Updated dataset in examples/example_fit_bao_sdss_dr16.py.
     - Fitting w in examples/example_fit_snia_cov.py.
     - Updating example in examples/pydata_simple.
     - Fixed reading of uninitalized memory in numcosmo/data/nc_data_snia_cov.c.
     - Fixed leak in numcosmo/math/integral.c.
     - Added error testing in numcosmo/math/ncm_csq1d.c.
     - Working in progress in APES and related objects.

 * Removed old CLAPACK and LAPACKE support.

 * W reconstruction (#58)

     * new object WSpline.
     
     * Better support for extrapolation for large redshifts in wspline object.
     
     * New Cosmic Chronometers data objects.
     
     * Added options to use polynomial interpolation for nknots < 6.
     
     * Improving error handling in NcmDiff.
     
     * Adding SDSS DR16 samples.
     
     * Added SDSS DR16 empirical fit objects.
     
     * New examples fitting BAO and Hz data.
     
     * Fixed gsl spline border problem.
     
     * Minor tweaks on documentation.
      Co-authored-by: Sander23 <sander23.ribeiro@uel.br>
 * Trying to find correct path due to broken glib in brew.

 * Debugging missing prereq.

 * New dependency resulting from split package in homebrew.

 * 56 lastest pantheon (#57)

     * Filter implementation for SNIa.
     * New SNIaCov example.
     * Improvements on SNIa constructors and example.
     * Updated autotools.
 * Minor updates in the figures of VacuumStudy and VacuumStudyAdiabatic.

 * Added volume method to nc_cluster_mass_nodist.

 * Adding Minkowski functions to CSQ1D.

 * Cleaning notebooks.

 * Added constructor annotation to ncm_mset_load(). New Vacuum study notebooks.

 * Added method to get the best fit from catalogs.

 * Included more frames for CSQ1D

 * Testing coverage tweaking.

 * Halo bias (#53)

     * halo bias branch
     
     * mean bias in the unbinned and binned case with proxies
     
     * Added the tests to the makefile. Updated the redshift object to one that
     implements p(z).
      Co-authored-by: Henrique Lettieri <henrique.cnl@hotmail.com>
 * Updated autotools.

 * Implementing frames in csq1d.

 * More tests for nc_data_cluster_ncount.c.

 * Removed option to print the mass function (old code).

 * Removed old method nc_data_cluster_ncount_print.

 * More tests for test_nc_data_cluster_ncount.c.

 * Removed inclusion of removed objects documentation.

 * Adding new integration routines to the ignore list in docs.

 * New integration code. Now vector integration used in nc_data_cluster_ncount.
     Fixed bug in NcmFitMC (it was using the bestfit from catalog instead of
     fiducial model to resample). Fixed typos.

 * Added a full corner plot comparing all outputs.

 * Updated generate_corner.ipynb to use ChainConsumer.

 * Fixed bug in catalog_load nc_data_cluster_ncount. New corner plot notebook.

 * Fixed minor leaks in ncm_reparam.c ncm_powspec_filter.c ncm_mset_catalog.c.
     Improved sampling in ncm_fit_esmcmc_walker_apes (now the second half use
     the updated first half when moving the walkers). Support for binning in
     nc_data_cluster_ncount. New notebooks comparing binning vs unbinning.

 * Inclusion of the time function to compare the effiency between CCL and Numcosmo

 * Reorganized binning options in NcDataClusterNCount.

 * Unbinned and binned analisys in the ascaso proxy

 * Reorganizing cluster mass ascaso object.

 * CCL- Numcosmo comparison using a mass proxy, both binned and unbinned analysis

 * Tests with de cluster abundance with a mass proxy

 * Proxy comparation

 * Fixed conflict leftovers.

 * Ascasp changes

 * Removed old data objects all binned versions now reside in NcDataNCount.
     NcABCClusterNCount needs updating. Now lenghts of cluster mass and redshift
     and class properties. Cluster abundance must be instantiated with both mass
     and redshift proxies defined. NcClusterMass/Redshift objects reorganized.

 * New helpers scripts with new tools: a function create pairs of NumCosmo/CCL
     objects, increase CCL precision and notebook plots with comparison between
     NumCosmo and CCL outputs. Updated notebooks to use helper functions.

 * Inclusion of the Cluster Number as a function of mass in the binned case both
     for CosmoSim and Numcosmo

 * Implementation of the inp_bin and p_bin_limits function in the
     gauss_global_photoz redshift proxy

 * Removed checkpoints and output files.

 * Implementation of binning in the lnnormal mass-observable relation

 * Binned and unbinned comparison between Numcosmo and CCL cluster abundace
     objects with no mass or redshift proxies

 * binned and unbinned comparison between CCL and Numcosmo cluster abundance with
     no mass or redshift proxies

 * Working version for binning proxies in NcCluster* family.

 * notebook on cluster mass comparison between CCL and Numcosmo update

 * addition of  binning in cluster_mass.c and cluster_mass.h and unbinning
     comparison between CCL and Numcosmo cluster mass objects(not ready yet)

 * Old modifications on hiqg and updates on NumCosmo vs CCL tests. Starting the
     implementation of binning for cluster mass and redshift.

 * Comparison between numcosmo and ccl cluster abundance objects

 * New spline object for functions with known second derivative. Updated
     nc_multiplicity_func_tinker to use interpolation objects, added option to
     use linear interpolation. Removed old notebook NC_CCL_Bocquet_Test2.ipynb.
     Updated NC_CCL_mass_function.ipynb (fixed bugs).

 * Mass functions comparisons notebook.

 * Updated version to match new interface.

 * Better limits for nc_halo_mass_function. Setting properties through gobject to
     catch out-of-bounds values.

 * Adjusted esmcmc run_lre minimum runs in tests.

 * Calibrated integrals to work on any point of the allowed parametric space.
     Added mores tests.

 * Modified ranges of concentration and alpha (Einasto) parameters.

 * Improved stability in nc_halo_density_profile.c.

 * Smaller lower bounds for ncm_fit_esmcmc_run_lre. Added option for starting
     value of over-smooth in mcat_analize calibration.

 * Option to calibrate over-smooth.

 * New minor version 0.15.4.


[v0.15.4]
 * New minor version 0.15.4.

 * Added missing ncm_cfg_register_obj call.

 * Delete NC_CCL_Bocquet_Test-checkpoint.ipynb
 * Delete .project
 * test of execution time

 * updates

 * New option to use kde instead of interpolation in APES.

 * Notebooks testing.

 * Moved model validating to workers (slaves or threads).

 * Improved fparam set methods.

 * Faster kde sampling.

 * Improved MPI debug messages added timming.

 * Improved control thread avoinding aggressive pooling by MPI.

 * Added conditional compilation of MPI dependent code.

 * Using switch to choose between kernel types.

 * Fixed memory leak.

 * the hydro and dm functions of the CCL were included

 * Finalized tests for kernel class.

 * Better handling of the case where 0 threads are allowed. Fixed limits on
     nc_cluster_photoz_gauss_global. Incresed lower limit in As in
     nc_hiprim_power_law.

 * Fixed leak.

 * New mpi run jobs async (master - slaves).

 * Added test for the kernel sample function.

 * Implentationg of tests for the #NcmStatsDistKernel class.

 * Configuration.

 * More tweaks on omegab range.

 * New notebook NC_CCL_Bocquet_Test has been created

 * Increasing lower limit of Omega_bh^2.

 * Fixed variable types for simulation (sim).

 * Improved bounds on nc_hicosmo_de_reparam_cmb.

 * Clean up Tinker: no need to set some parameters as properties. Delta -
     CONSTRUCT and not CONSTRUCTED_ONLY

 * Fixed bug in Bocquet multiplicity function, e.g., properties are CONSTRUCT not
     CONSTRUCT_ONLY.

     Clean up Crocce, Jenkins and Warren's multiplicity functions: no need to
     set the parameters as properties.

 * Fixing mcat_analize to work with small catalogs.

 * Resolved conflict.

 * Fixed minor warnings.

 * Implemented Bocquet et al. 2016 multiplicity function. Two new functions in
     NcMultiplicityFunc: has_correction_factor and correction_factor. Bocquet
     provides fits for mean and critical mass definitions, but the latter
     depends on the first.

 * Fixed some edge cases in ncm_fit_esmcmc.c and ncm_fit_esmcmc_walker_apes.c.
     Minor reorganization.

 * Add files via upload
 * Testing 10D.

 * Added missing object registry.

 * Removed incomplete tests.

 * Fixed allocation problem in ncm_stats_dist.c. Fixed other minor bugs and
     tweaks.

 * Included tests for the error messages in stats_dist_kernel.c

 * Test if kernel test is implemented right.

 * Finished the documentation of ncm_fit_esmcmc_walker.c and
     ncm_fit_esmcmc_walker_apes.c

 * Implemented documentation of ncm_fit_esmcmc_Walker.c

 * Updated automake file.

 * uncrustify and more tweaks on test_ncm_diff.c removing edge cases.

 * Fixed internal struct access.

 * uncrustify.

 * Updated test, and fixed minor issues.

 * Fixed docs and set nc_multiplicity_func.c to abstract.

 * Refactoring of the multiplicity function object is complete. Main difference:
     included mass definition as a property. Examples were properly updated.

 * Tweaked test_ncm_mset_catalog.c and test_ncm_diff.c. Solved APES offboard
     sampling.

 * Changed the size of figures in docs and improved the documentation of
     StatsDistKernel objects.

 * Improving coverage and fixed casting.

 * Generating graphs with the notebooks.

 * Included over-smooth option in APES. Added the same option to darkenergy's
     command line interface. Improved documentation and coverage.

 * Improved tests and coverage for NcmStatsDist* family.

 * Improving NcmStatsDist* coverage.

 * More tweaks on NcmDiff tests.

 * Tweaking tests to avoid false positives.

 * Improved unit tests for NcCBE, NcCBEPrecision and NcmVector.

 * Improved interface to NcmFitESMCMCWalkerAPES. Included and tweaked unit tests.

 * I am rewriting the multiplicity function objects. Including missing functions
     (e.g., ref, free, clear...), put in the correct order. Add "mass
     definition" as a property.

 * Fix documentation glitches and solve warnings.

 * Documentation for stats dist objects with image problems

 * Unfinished stats dist objects documentation

 * Removed whitespace following trailing backslash.

 * Added missing include directory.

 * Working on stats_dist.c documentation

 * Uncrustify tests. Tweak mcmc tests.

 * Fixed bug in ncm_mset_trans_kern_cat.c (re-preparing for each sampling). Added
     missing files. Added new test to test_ncm_vector.c. Tweaking tests.

 * uncrustify and rename APS to APES.

 * Fixed wrong href when computing IM in VKDE. Fixed over_smooth tweak in
     prepare_interp.

 * Removed unecessary files. Added notebooks.

 * Working version of ncm_stats_dist*. Not yet fully tested.

 * First (incomplete) reorganized version of NcmStatsDist*. Updated mkenums
     templates.

 * Updated notebook. Halo profile uses log10(M) instead of M. Modifying
     Multiplicity function objets: mass definition is a property. Work in
     progress.

 * Removed CNearTree.

 * Working version of vbk.

 * New script to use numcosmo without installing.

 * Testing for fit with no free parameters bug. Fixed the same bug in fit impls.

 * Fixed indentation.

 * Removed unecessary files.

 * vbk_studentt working on notebook. Memory error for rosenbrock. Check slack for
     info.

 * Adding support for non-adiabatic computation.

 * vbk_studentt working for eval and evan_m2lnp. Copy of APS to work with vbk (not
     included in makefile). Copy of gauss to gauss vbk(included in makefile)

 * Missing files from last commit.

 * Functions prepare_args and preapre_interp running. Starting to work on
     eval_m2lnp. Interp.py is the test file.

 * Added more precise delta_c.

 * Working on the examples.

 * Working on VBK.

 * Updated autotools file and removed binnary.

 * example_neartree is the example from documentation, test_neartree is build by
     me and slightly documented.

 * Fixed the includes for CNearTree, inserted a flag in Makefile.am and created a
     test to check.

 * Added gtk-doc to mac os build.

 * removed azure.

 * removed azure.

 * Removed travis-ci.

 * Updated autotools files and removed travis-ci.

 * Added the required files for CNearTree.c library, created copies of stats dist
     to work on, and added the necessary lines in the makefiles.

 * Added gtk-doc to mac os build.

 * Funnel example and notebook.

 * New test likelihood Funnel.

 * Included the RoT for the Student t distributions in
     ncm_stats_dist_nd_kde_studentt (truncated for nu < 3.0 since it is not
     defined for these values).

 * Set default to aps with studentt (Cauchy dist) with no CV and over smooth 1.5.

 * Fixed bug in ncm_stats_dist_nd (it didn't set weights vector to zero before
     fitting).

 * New notebook used to plot Rosenbrock MCMC evolution.

 * New Rosenbrock model/likelihood to check MCMC convergence. New option to thin
     chains. New example to run Rosenbrock MCMC.

 * Typo fix.

 * Working version of the reorganized code (NcmNNLS, NcmISet and NcmStatsDistNd).

 * Working on ncm_stats_dist_nd + ncm_nnls. Working version, finishing code
     reorganization.

 * Working version (not organized yet, full of debug prints...).

 * Added documentation and comentaries in the  notebook TestInterp.ipnb.

 * Improved the description in the documentation of ncm_stats_dist_nd.c,
     ncm_stats_dist_nd_studentt.c and ncm_stats_dist_nd_gauss.c.

 * Moved headers to the right place.

 * Missing Makefile.am.

 * Moved external codes to a new (sub)library to remove these codes from the
     coverage and to make the symbols invisible.

 * Minor release v0.15.3.

 * Added interpolation case where only the most probable point is necessary.

 * Added tests for KDEStudentt.

 * Fixed a few documentation glitches.

 * Reorganized ncm_stats_dist_nd* objects family. Testing different solvers to the
     NNLS problem.

 * Change on the file numcosmo-docs.sgml to include
     ncm_stats_dist_nd_kde_studentt.c. Did not create a studentt HTML as I
     expected.

 * Reupdated m4 and automake stuff.

 * Minor identation/positional tweaks.

 * Implementation of the comentaries from the commit "New implementation of
     studentt function for ncm_stats_nd_kde.".

 * (Re)updated m4 macros and gtk-doc.make.

 * Adding the updated Jupyter notebook

 * New implementation of studentt function for ncm_stats_nd_kde.

 * Updated private instance get function. Fixed doc issues.

 * Better support for arb.

 * New example.

 * Missing file in branch.

 * Testing

 * First version of the ncm_powspec_sphere_proj and ncm_fftlog_sbessel_jljm.

     Computing Cells without RSD is already working.

     Several improvements and extensions in other objects.

 * New FFTLog object to compute the integral with the kernel j_lj_m.


[v0.15.3]
 * Minor release.

 * Included a function to compute numerical integrals of the NFW profile (instead
     of the analytical forms). To be used for testing only!

 * New methods to access Ym values in NcmFftlog.

 * Added a second run to avoid unfinished minimization process.

 * Removed debug msgs from coverage build.

 * Removed coverage flags from introspection build.

 * Debug coverage build.

 * Debug coverage build.

 * Debug coverage build.

 * Removed LDFLAGS for coverage.

 * Debug coveralls build.

 * Moved (all) flags to the right places.

 * Moved flags to the right place.

 * Added explict CODE_COVERAGE_LIBS to introspection build.

 * Debug coveralls build.

 * Test speedups.

 * Allowed reasonable failures.

 * Added 10% allowed test errors when estimating hessian computation error.

 * Testing ncm_stats_dist_nd_kde_gauss.c. Minor modifications to
     ncm_data_gauss_cov_mvnd.c. New notebook to test multidimensional
     interpolation.

 * Debug mac-os GHA

 * Debug mac-os GHA

 * Debug mac-os GHA

 * Debug mac-os GHA

 * Debug mac-os GHA

 * Debug mac-os GHA.

 * Debug mac-os GHA.

 * Debug mac-os GHA build.

 * Trying reinstalling gmp.

 * Testing a solution for GHA on mac-os.

 * Still debugging macos build in GHA.

 * Debug macos build.

 * Conditional use of sincos.

 * Fixed sincos warning.

 * More compiler env.

 * Fixed sincos included warning.

 * Updated example.

 * Setting compilers.

 * Cask install for gfortran in macos build.

 * Testing lib dir in GHA.

 * Added cask install fortran for macos build.

 * Trying lib dirs.

 * Added gfortran req to macos build.

 * Added prefix option to configure in GHA.

 * Fixed example name and moved test.

 * Rolled back autoconf version req.

 * Included missing make install in build check.

 * Updated autotools and deps. New check in GHA. Fixed bug in numcosmo.pc.in.

 * Working on ncm_csq1d.c. New notebook FisherMatrixExample.ipynb.

 * Adding timezone info.

 * New docker image with NumCosmo prereqs.

 * Working on nc_de_cont.

 * Running actions in every branch.

 * Testing GHA

 * Testing GHA

 * Testing GHA

 * Testing GHA

 * Testing coveralls build.

 * Updated CI badge to GHA.

 * Adding missing prereq for the macos build.

 * Better workflow name and removed unnecessary prereq in the macos build.

 * Adding macos build.

 * Removed debug print in c-cpp.yml.

 * Adding references.xml to the repo.

 * Update c-cpp.yml

     Checking xml logs
 * Update c-cpp.yml

     Debug xml build
 * Fixed doc typo.

 * Updated to sundials 5.5.0.

 * Added NumCosmo CCL test notebook.

 * Fixed conditional compilation for system with gsl < 2.4.


[v0.15.2]
 * Minor release 0.15.2.

 * Updated tests and fixed indentation.

 * New framework for Cluster fitting with WL data (in progress).

 * New minor version. Reorganizing WL likelihood (in progress).


[v0.15.1]
 * Default refine set to 1.

 * More options to refine.

 * Add refine as an option.

 * Added vectorized interface for nc_wl_surface_mass_density_reduced_shear. Minor
     other improvements.

 * Improvement in ncm_spline_func to remove outliers.

 * Update c-cpp.yml

     Removed distcheck (split into check and dist)
 * Update c-cpp.yml

     Added --enable-man
 * Update c-cpp.yml

     Added reqs for doc building.
 * Update c-cpp.yml

     Added doc building
 * Update c-cpp.yml

     Added support for gtk-doc
 * Update c-cpp.yml

     added parallel build and upload artifact
 * Update c-cpp.yml

     More deps and removed double configure.
 * Update c-cpp.yml

     More missing deps
 * Update c-cpp.yml

     Added more missing deps
 * Update c-cpp.yml

     Added missing dep: gfortran.
 * Update c-cpp.yml

     Testing install prereqs.
 * Update c-cpp.yml

     Added prereqs
 * Updated and finished support for CCL in Dockerfile-clmm-jupyter.

 * Added python3-yaml support.

 * Added support for CAMB and CCL.

 * Support for CCL and CAMB.

 * Adding support for camb and ccl.

 * Fixed new filename.

 * Notebook comparing Colossus and CCL with NumCosmo: density profiles, surface
     mass density and the excess smd.

 * Added hook between nc_halo_mass_function and ncm_powspec_filter to ensure the
     same redshift range.

 * Fixed minor leak.

 * Minor fixes and improvements in ncm_spline_func_test.*.

 * Update c-cpp.yml

     Testing github actions
 * Create c-cpp.yml

     Trying the github actions.
 * Included new distance functions from z1 to z2.

 * Test suit NcmSplineFuncTest is now stable enough. Next step: add some
     cosmological functions examples.

 * Doc. minor changes.

 * Memory leak - ncm_spline_new_function_4

 * Added new outlier function to last description example.

 * Corrected vector memory lost.

 * Added option to save outliers grid to further analysis.

 * Added test suite to NcmSplineFunc.

 * Tutorial reviewed.

 * Text review (Mari).

 * Added missing ipywidgets from docker build.

 * Fixed makefiles.

 * Reorganized and added copyright notices to notebooks.

 * New tutorial.

 * New code for homogeneous knots.

 * ncm_vector.c documentation improve.

 * Doc. improvement.

        * NcDataCMB
       * NcDataCMBShiftParam
       * NcDataCMBDistPriors

 * Minor doc. modification.

 * Minor doc. changes.

 * Added support for abstol in NcmSplineFunc.

 * Set max order to 3 in NcmODESpline to make the ode integration tolerance agree
     with spline interpolation error.

 * Removed old CCL interface (they no longer have a C api, we are moving to test
     in python since their API is only there).

 * Improve doc. & fixed indentation:

        * NcPowspecML
       * NcPowspecMLFixSpline
       * NcPowspecMLTransfer
       * NcPowspecMLCBE
       * NcPowspecMNL
       * NcPowspecMNLHaloFit

 * Transfer function improve doc.

 * NcTransferFuncEH: fixed indentation.

 * NcTransferFuncEH: improve doc.

 * NcTransferFuncBBKS: corrected minor typo.

 * NcTransferFuncBBKS: fixed indentation.

 * NcTransferFuncBBKS: improve doc. Add BBKS ref.

 * NcTransferFunc: fixed indentation.

 * NcTransferFunc: improve doc.

 * NcWindowGaussian & NcWindowTophat: standardization between both descriptions.

 * NcWindowGaussian: fixed indentation.

 * NcWindowGaussian: improve doc.

 * NcWindow: fixed indentation.

 * NcWindow: improve doc.

 * NcWindowTophat: fixed indentation.

 * NcWindowTophat: improve doc.

     Note: it seems that latex commands "cases" and "array" does not work.

 * More debug messages in MPI.

 * Better debug messages and identation.

 * Documentation.

 * Fixed details in the documentation.

 * NcmPowspecFilter: fixed indentation.

 * NcmPowspecFilter: doc. improvement.

 * NcmPowspec: reference to function NcmPowspecFilter in ncm_powspec_var_tophat_R
     ()

 * NcmPowspec: Fixed indentation.

 * NcmPowspec: doc. improvement.

 * Minor typo.

 * NcmODEEval: Fixed indentation.

 * NcmODEEval: doc. improvement.

 * NcmODE fixed indentation.

 * NcmODE doc. improvement.

 * NcmSpline2dBicubic: fixed indentation and tweak doc.

     *The documentation still needs lots of work.*

 * Fixed minor typos.

 * Fixed indentation:

        * ncm_spline2d_spline.h/c
       * ncm_spline2d_gsl.h/c

 * NcmSpline2dSpline and NcmSpline2dGsl doc tweaks.

 * NcmSpline2d: fixed indentation.

 * NcmSpline2d: Added Include and Stable tags + Minor tweaks.

 * Fixed indentation:

        * ncm_fftlog_tophatwin2.h/c
       * ncm_fftlog_gausswin2.h/c

 * NcmFftlogTophatwin2 and NcmFftlogGausswin2: doc. improvement.

 * Fixed indentation: ncm_powspec_corr3d.c/h.

 * NcmPowspecCorr3d: doc. improvement.

 * NcmFftlogSBesselJ: fixed description and minor tweaks.

 * NcmFftlogSBesselJ: tiny tweaks in the description.

 * Fixed indentantion.

 * NcmFftlogSBesselJ: corrected indentation.

 * NcmFftlogSBesselJ: documentation improved.

 * Fixed wrong lower bound for abstol.

 * Removed old test in autogen.sh and overwritting of gtk-doc.make.

 * Add gtkdoc related files (instead of soft links).

 * Added to repo all necessary m4 files.

 * NcmGrowthFunc: doc tiny tweaks

 * NcmFftlog: corrected indentation.

 * NcmFftlog: documentation's minor improvement.

 * Tweaking NcGrowthFunc documentation and fixed wrong link for NcmSplineFunc.

 * NcGrowthFunc: changed description to a vague explanation on the initial 
     conditions. Added Martinez and Saar book on the references.

 * NcGrowthFunc: corrected indentation.

 * NcGrowthFunc: improved documentation.

 * NcDistance: Standardization of function documentation

        *_free.c
       *_clear.c
       *_ref.c

 * NcDistance: modified two static functions names:

        * comoving_distance_integral_argument -->
     _comoving_distance_integral_argument
       * dcddz --> _dcddz

 * NcDistance: corrected indentation.

 * NcmDistance: improved documentation.

 * Testing support for gcov.

 * ncm_timer.* - correct indentation with uncrustify.

 * NcmTimer: improved documentation.

 * Corrected a broken link in short description.

 * Indentation using uncrustify.

 * Improved documentation from NcmRNG.

 * Changed "abs" --> "abstol" in ncm_ode_spline_class_init. Also some minor
     changes.

 * Improved #NcmOdeSpline documentation.

 * Using different branch in CLMM.

 * Added colossus to Dockerfile-clmm-jupyter.

 * Version 0.15.0

 * Corrected indentation of ncm_ode_spline.*

 * [DOC] Improved main and enum description.

        * ncm_spline_func.h
       * ncm_spline_func.c

 * Removed punctuation from parameter descriptions:

        * ncm_spline.c
       * ncm_spline_cubic_notaknot.c
       * ncm_spline_gsl.c

 * Removed punctuation from parameter descriptions in ncm_matrix.h/c.

 * Removed punctuation from parameter descriptions. Added bindable function to
     NcmSplineFunc.

 * Added description/documentation to ncm_spline_func.h/c

     Corrected missing link in NcmC.

 * Fixed minor bugs.

 * Added part of doc from spline_func module.

 * Included "@stability: Unstable" and "@include: numcosmo/math/ncm_spline_rbf.h".

     Added NcmSplineRBFType enum description.

     Added doc in function ncm_spline_rbf_class_init
     (g_object_class_install_property).

 * Added NcmSplineGslType enum description.

     Added doc in function ncm_spline_gsl_class_init
     (g_object_class_install_property).

     Corrected indentation of ncm_spline_gsl.c & ncm_spline_gsl.h with
     uncrustify.

 * Added "@stability: Stable" and "@include:
     numcosmo/math/ncm_spline_cubic_notaknot.h".

     Corrected indentation of ncm_spline_cubic_notaknot.c & 
     ncm_spline_cubic_notaknot.h with uncrustify.

 * Corrected indentation of ncm_spline_cubic.c and ncm_spline_cubic.h with 
     uncrustify.

 * Added doc to functions:
       * ncm_spline_is_empty
       * ncm_spline_class_init (g_object_class_install_property)

     Correct indentation of ncm_spline.c and ncm_spline.h with uncrustify.

 * Added doc in functions ncm_vector_class_init & ncm_matrix_class_init.

       /**
       * NcmMatrix:values:
       *
       * GVariant representation of the matrix used to serialize the object.
       *
       */

     And the same for NcmVector.

 * Some minor tweaks:

        * ncm_matrix.c
       * ncm_matrix.h
       * ncm_vector.c

 * Passed ncm_matrix.h and ncm_matrix.c through uncrustify to set indentation.

 * Added doc to the following functions of NcmMatrix:

        * ncm_matrix_get_array
       * ncm_matrix_fast_get
       * ncm_matrix_fast_set
       * ncm_matrix_gsl
       * ncm_matrix_const_gsl
       * ncm_matrix_col_len
       * ncm_matrix_row_len
       * ncm_matrix_nrows
       * ncm_matrix_ncols
       * ncm_matrix_add_mul
       * ncm_matrix_cmp
       * ncm_matrix_cmp_diag
       * ncm_matrix_dsymm
       * ncm_matrix_set_colmajor
       * ncm_matrix_get_diag
       * ncm_matrix_set_diag
       * NCM_MATRIX_SLICE
       * NCM_MATRIX_GSL_MATRIX
       * NCM_MATRIX_MALLOC
       * NCM_MATRIX_GARRAY
       * NCM_MATRIX_DERIVED

     Also deleted NcmMatrix struct doc.


[v0.15.0]
 * Version 0.15.0

 * Polishing nc_halo_density_profile. Added NumCosmo x Colossus comparison
     notebook.

 * Fixed missing parameter doc.

 * Updated test test_nc_wl_surface_mass_density.

 * Better integration strategy for NcHaloDensityProfile. Updated test
     test_nc_halo_density_profile.

 * Final tweaks before release.

 * Minor improvements in notebooks/BounceVecPert.ipynb.

 * Added new profile (Hernquist), Einasto implementation is now complete. Added
     documentation.

 * Improved notebooks/BounceVecPert.ipynb.

 * Removed old nlopt header in csq1d.

 * Fixed indentation and conditional load of NLOPT library object.

 * Implemented Einasto profile (just rho, not the integrals). Included the
     funciton to compute the magnification.

 * New refactored NcHaloDensityProfile (working in progress). Added support for
     different internal checkpoints in NcmModel.

 * Added numcosmo's uncrustiify settings.

 * Uniform indentation.

 * Uniform indentation.

 * Homogenization and standardization of the #NcmVector module.

     A bunch of minor changes, .e.g. all #NcmVector are @cv now. Minimal changes
     in some docs.

     Added description at enum NcmVectorInternal.

 * Uniform indentation.

 * Finished first version of #NcmVector documentation. Still needs a careful
     check.

 * Renamed NcDensity* objects to NcHaloDensity*.

 * Added "@stability: Stable" and "@include: numcosmo/math/ncm_c.h" to section in
     ncm_c.c.

 * Fixed use of Planck likelihood without check_param. Fixed typo in
     numcosmo/math/ncm_c.c.

 * 1) Changed function name: ncm_c_hubble_cte_planck_base_2018 to 
     ncm_c_hubble_cte_planck6_base.

     2) Added doc to functions:
       * ncm_c_blackbody_energy_density
       * ncm_c_blackbody_per_crit_density_h2

     3) Corrected units in functions ncm_c_h() and ncm_c_hbar(): Js^{-1} to Js.

     4) Deleted last function ncm_c_radiation_h2Omega_r0_to_temp Not used 
     anywhere.

 * Now the last commit is correct.

 *    * ncm_c_crit_density_h2
       * ncm_c_crit_mass_density_h2

 * 1) Deleted function from #NcmC:
       * ncm_c_hubble_cte_msa - it was not applied anywhere.

     2) Replaced function from #NcmC:
       * ncm_c_hubble_cte_wmap H0 = 72 --> ncm_c_hubble_cte_planck_base_2018 H0
     = 67.36

     3) Replaced ncm_c_hubble_cte_wmap to ncm_c_hubble_cte_planck_base_2018 in
     functions from module #NcHIcosmo*
       * nc_hicosmo_de
       * nc_hicosmo_gcg
       * nc_hicosmo_idem2
       * nc_hicosmo_qconst
       * nc_hicosmo_qgrw
       * nc_hicosmo_qlinear
       * nc_hicosmo_qspline

     4) Added Planck reference in references.bib and references.tex

     5) Added docs in functions from #NcmC:
       * ncm_c_hubble_cte_planck_base_2018
       * ncm_c_hubble_radius_hm1_Mpc
       * ncm_c_crit_number_density_p
       * ncm_c_crit_number_density_n

 * 1) Changed documentation to the following #NcmC module functions:
       * ncm_c_wmap5_coadded_I_K
       * ncm_c_wmap5_coadded_I_Ka
       * ncm_c_hubble_cte_hst

     2) Added documentation to the following #NcmC module functions:
       * ncm_c_wmap5_coadded_I_Q
       * ncm_c_wmap5_coadded_I_V
       * ncm_c_wmap5_coadded_I_W

 * Added documentation to the following #NcmC modules functions:

        * ncm_c_wmap5_coadded_I_K
       * ncm_c_wmap5_coadded_I_Ka
       * ncm_c_hubble_cte_hst

 * Documented function ncm_vector_len.

 * Better expansion for tan(x+d)-tan(x) for small d.

 * Fixed corner case in numcosmo/model/nc_hiprim_atan.c.

 * Fixed MPI in hdf5 incompatibility.

 * Removed update option on homebrew.

 * Fixed error handling, clik returns wrong values when an error occurs (due to a
     wrong usage of forwardError), to fix this we changed the likelihood to
     return m2lnL = 1.0e10 whenever clik returns an error.

 * Workaround to fix travis ci bundle issue.

 * Removed wrong free in CLIK_CHECK_ERROR.

 * Updated old m4 files and building system to keep them updated in the m4/
     folder.

 * New MPI server (in progress). NcDataPlanckLKL no longer kills process when clik
     returns an error (just sends a warning).

 * Better OpenMP (and others) number of threads control.

 * Finished update to 2018 Planck likelihood (in testing).

 * Updates in the notebook.

 * Fixed sprintf related warnings.

 * Removed typo.

 * Fixing new docs build process...

 * Fixing error messages in plik. Working on doc building process.

 * Fixing new docs build process...

 * Reorganized docs building process.

 * Testing travis-ci on osx.

 * Testing travis-ci on osx.

 * Testing travis-ci on osx.

 * Testing travis osx.

 * Fixed doc typos. Testing travis on osx.

 * Fixed doc in NcDensityProfile. Testing travis ci on osx.

 * Updated sundials to version 5.1.0. Fixed somes tests and updated to TAP.

 * Update MagDustBounce.ipynb from Emmanuel Frion.

 * Minor updates and annotation improvements.

 * Copy all examples in Dockerfile-clmm-jupyter.

 * Updated notebooks and Dockerfile-clmm-jupyter.

 * New parametrization for CSQ1D, updates in BounceVecPert.ipynb and
     MagDustBounce.ipynb.

 * Several improvements in MagDustBounce.ipynb.

 * New NcHICosmoQRBF model.

 * Added options to logger function to redirect all library logs.

 * Testing CLMM+NumCosmo notebook.

 * Removed debug print.

 * Fixed documentation error.

 * Added missing scipy for CLMM.

 * Added missing Astropy for CLMM.

 * Cloning the right branch from CLMM.

 * Copying examples from CLMM to work.

 * Added COPY from opt.

 * Fixed build script name...

 * Testing different build order.

 * Testing build from git.

 * New Dockerfile for CLMM comparison.

 * Updated to CODATA 2018. Reorganized density profile objets (in progress).

 * Fixed test test_nc_ccl_dist.c, decreased number of tests in test_ncm_fftlog.c
     and test_ncm_mset_catalog.c. Testing ncm_csq1d.c. Updated
     binder/Dockerfile.

 * Adding missing notebooks to _DATA.

 * Updated binder notebook.

 * Removed debug messages in CSQ1D.

 * Minor fixes and new notebook.

 * Missing notebooks in Makefile.

 * Updated image used by binder.

 * Tweaking notebooks.

 * Included two notebooks.

 * Testing different methods to deal with zero-crossing mass.

 * New tutorial notebooks.

 * Modifying density profile objects. E.g., including more mass definitions.

 * Using the full SHA hash.

 * Testing a mybinder using a Dockerfile.

 * Testing methods to integrate regular singular points.

 * Testing docker jupyter notebooks

 * Testing docker jupyter notebooks

 * Testing docker jupyter notebooks

 * Testing docker jupyter notebooks

 * Testing docker jupyter notebooks

 * Testing docker jupyter notebooks

 * Testing docker jupyter notebooks.

 * Added support for matplotlib and scipy in the docker image.

 * Updated to python3 on Dockerfile.

 * Working on Dockerfile.

 * Working on Dockerfile.

 * Added necessary dist.prepare to examples.

 * Working on Dockerfile.

 * Working on dockerfile.

 * Updating Dockerfile.

 * Fixed lgamma_r declaration presence test. Fixed glong/gint64 mismatch.

 * Added detection for lgamma_r declaration and workaround when it is not declared
     but present (usually implemented by the compiler).

 * Fixed many documentation bugs (in most part by adding __GTK_DOC_IGNORE__ to the
     inline sections).

 * Added the GSL 2.2 guard back to where it was really necessary.

 * Fixing documentation bugs.

 * Removed old GSL guards from tests.

 * Added prepare_ functions on NcWLSurfaceMassDensity object. Added necessary
     prepare calls to tests.

 * Support for gcc 9 in macos.

 * Added missing header.

 * Fixed inlined sincos to use default c keywords.

 * Changed inline macro to the actual keyword in config_extra.h.

 * Several updates and fixes.

     - Updated config.h, now it includes config_extra.h that contains local
     compile only functions necessary to NumCosmo.
     - Updated function_cache to a proper object.
     - New NcmCSQ1D object that implemest the complex structure quantization
     method.
     - New macro to control inlining (Glib's deprecated their).
     - Updated inlining macro at all necessary headers.

 * Added rcm to ignore list in docs.

 * Added missing CPPFLAGS for HDF5.

 * Testing HDF5 in ubuntu.

 * Fixing missing doc. Testing HDF5 in ubuntu.

 * Fixed bug in preparing a fparams_map with no free variables.

 * Removed fitting test in APS.

 * Testing APS.

 * Fix pitch.

 * Improving edge cases and validation.

 * New validation on mset.

 * Removed old debug/test print in class.

 * Updated and reorganized SNIa catalogs support. Added Pantheon.

 * Working on new Boltzmann code.

 * Working on new Boltzmann code.

 * Set up CI with Azure Pipelines

     [skip ci]
 * Fix transport vectors.

 * Fixed instrospection tags.

 * Fix email.

 * Updated to new sundials interface.

 * Updated Planck likelihood code (not tested).

 * Minor changes in HIPrimTutorial.ipynb. Updated ccl interface.

 * Minor modifications related to some tests comparing with "Cluster toolkit".

 * New function to update Cls.

 * Minor improvements.

 * Added missing test file.

 * Fixed indentation.

 * New object NcmPowspecCorr3d (for the moment computes the simples 2point
     function). Moved filter functions from NcPowspecML to NcmPowspec. Improved
     bounds sync in Halofit.

 * Update README.md
 * Converting python examples to jupyter notebooks.

 * New Nonlinear Pk tests. Finished nc_powspec_mnl_halofit encapsulation (private
     members).

 * Fixes xcor to work with sundials 4.0.1.

 * Updating to new sundials API.

 * Updating to the new sundials API.

 * Fixed conditional use of OPENMP (mainly for clang).

 * Fixed merge leftovers.

 * Several improvements in the new Boltzmann code (pre-alpha). Minor fixes.

 * Reorganized the hipert usage of bg_var. Moved gauge enumerator to Grav
     namespace.

 * Added support for system reordering to lower the bandwidth. New component PB
     photon-baryon and gravitation Einstein.

 * New abstract class describing first order problems. Fixed enum type names.

 * New abstract class to describe arbitrary perturbation components.

 * First commit of new perturbation module.


[v0.14.2]
 * Final tweaks for the v0.14.2 release.

 * Working on ccl vs numcosmo unit tests. Minor improvements.

 * Minor version update 0.14.1 => 0.14.2.

     Minor improvements in ncm_ode_spline.

 * Added missing header (in some contexts).

 * Fixed aliasing problem in ncm_matrix_triang_to_sym. Removed log info from
     travis-ci.

 * Fixing doc issues, added missing docs. Better debug message for
     ncm_matrix_sym_posdef_log. Fixing travis-ci.

 * Added debug to travis-ci. Included sundials at the ignore list for
     documentation.

 * Removed no python option in numpy.

 * Lapack now is required, added openblas and lapack to travis-ci. Added
     no-undefined (when available) to libnumcosmo.

 * Added redshift direction tolerance for NcmPowspecFilter. New unit test
     CCLxNumCosmo test_nc_ccl_massfunc.

 * Finished upgrade to CLASS 2.7.1. Added option to hide symbols of dependencies.
     Finished CCL tests for background, distances and Pk (BBKS, EH, CLASS).
     Minor tweaks.

 * Updated CLASS to version 2.7.1.

 * Fixed a typo in the transverse distance. Test distances: Nc and CCL.

 * Missing unit test file.

 * Added warning for initial point in minimization not being finite.

 * Added support for CCL, first unit test for NumCosmo and CCL comparison. Minor
     improvements in NcHIQG1D.

 * Missing file.

 * Moved to xenial in travis ci.

 * Removed backports repository in travis-ci (it no longer exists...), waiting for
     something to break.

 * Removed sundials as dependence in travis-ci.

 * Updated to new sundials version -- 4.1.0.

 * Updated directory structure to match that of the new version (4.1.0).

 * Encapsulating sundials version 4.0.2. Many additions and improvements.

     Reorganized all ode interfaces to match last sundials release. Testing a
     new abstract ODE framework (NcmODE*). New lapack functions added and old
     functions migrated to use NcmLapackWS. New matrix tests (also migrated to
     the newer unit test interface). New matrix log and exp functions added. 
     Testing the fit of the whole covariance matrix in NcmStatsDistNdKDEGauss. 
     Several minor improvements.

 * Added calibration objects for the reduced shear.

 * Working on the new ODE interface.

 * Fixed version mismatch.

 * Added back support for sundials 2.5.0 (as used by Ubuntu trusty).

 * Better support for Sundials versions, now it detects the version automatically.
     Minimum Sundials version is 2.6.0, minor updates to codes using Sundials.

 * New sampling options for NcmMSetTransKernCat object, testing new sampling
     options in test_ncm_fit_esmcmc.

 * Removed non-unsed typedef.

 * Removing unecessary headers.

 * New test to determine the burnin phase (better fitted for low self-correlation
     samplers).

 * Travis OK, removing log.

 * Travis...

 * Still testing travis.

 * Testing travis builds.

 * Fixing travis syntax.

 * Triggering travis.

 * Cleaning and updating NcmHOAA, fixing travis ci bug in mac os.

 * Fixed mac os image.

 * Testing travis ci, mac os config.

 * Debugging travis glitch in mac os.

 * Updated travis.yml to match new environment (again...)

 * Fixed minor documentation bugs.

 * Added new header for fortran lapack functions prototypes.

 * Added conditional macro compilation of the new suave support in Xcor.

 * New walker `Approximate Posterior Sampling' APS based on RBF interpolation
     using (NcmStatsDistNdKDEGauss).

 * Trying a different approach for the new walker.

 * Tests for the multidimensional kernel interpolation/density estimation object.

 * Added two codes for quadratic programming gsl_qp (gsl extras) and LowRankQP
     (borrowed from R). New ESMCMC walker Newton (does not work as expected,
     transforming in another sampler in the next commit). New multidimensional
     kernel interpolation/density estimation for arbitrary distribution
     (abstract interface and gaussian kernel implementation).

 * Removed an exit() in example_diff.py. Added support for new version for
     sundials.

 * Support for the new Sundials version.

 * Updated tests for the new interface for lnnorm computation (including error).

 * New tool for trimming catalogs, fixed typos in parameters names in NcHIPrim*.
     Fixed error estimation in the posterior normalization.

 * Fixed minor issues from codacy.

 * Fixing macos+travis issue.

 * Removed old includes in tests. Added debug in travis.

 * Moved MVND objects to main library code. Improved estimates of the Bayes
     factor, included new unit tests.

 * Change tau range. My tau definition is different from root (CERN) code. Voigt
     profile.

 * Fixed: updated deprecated glib functions. Included new Sundials version in
     configure.

 * Fixed and tested (against CCL) the galaxy weak lensing module inside XCOR.

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * Testing PITCHME.md

 * New PITCHME.md presentation. Tweaks on pydata_simple example. Minor bug fixes.

 * Included missing property "dists" (Photo-z distributions) in class_init.

     Included NcDataReducedShear... and NcGalaxyRedshift... in ncm_cfg.c

 * Included other variable types to be read from hdf5 catalog.

 * Added a fix to avoid setting null sparam-array.

 * Proper serialization of model with modified parameter properties and no
     reparametrization.

 * Fixed stupid bug.

 * Turned backport function static to avoid double definition.

 * Add option to choose parameter by name in mcat_analize.

 * New options on mcat_analize to control dump. Updated parameters range in
     nc_hicosmo_de and nc_hicosmo_lcdm.

 * Fixed missing doc.

 * Finished ESMCMC MPI support.

 * New bindable ncm_mset_func1 abstract class. Included function array in
     example_esmcmc.py. Fixed bug in function array + ncm_fit_esmcmc + MPI.

 * Fixed flags search when using ifort/icc. Added missing objects registry.

 * Added guard to avoid parsing blas header during g-ir-scanner.

 * Moved HDF5 LDFLAGS to LIBS.

 * Re-included gslblas inclusion/removal in configure.ac.

 * Included additional sundials library to LIBS.

 * Inverted header order to avoid clash with CLASS headers. Added backport from
     glib 2.54.

 * Better serialization of constraint tolerance arrays.

 * Safeguards in NcmVector constructors.

 * Fixing minor leaks.

 * Fixed initializer (to remove harmless warning).

 * Added conditional compiling of _nc_hiqg_1d_bohm_f.

 * Fixed MPIJob crash when MPI is not supported.

 * Fixed leaks in ncm_fit_esmcmc.c, ncm_mpi_job_fit.c and ncm_mpi_job_mcmc.c.

 * Updated tests and removed debug info from travis-ci.

 * Fixed blas detection (ax_check_typedef is broken! Now using AC_CHECK_TYPES).

 * Debugging travis-ci build.

 * Debugging travis-ci macos build.

 * Typo in blas enum detection.

 * Detecting lapack xblas functions availability. Fixing cblas headers in
     different scenarios.

 * Fixing BLAS headers compatibility.

 * Included new functions on NcWLSurfaceMassDensity: critical surface mass
     density, shear and convergence new functions are computed when the source
     plane is at infinite redshift.

     Bug fixed in nc_data_reduced_shear_cluster_mass.c: the probability function
     of the reduced convergence is computed using the reduced shear at infinite
     redshift.

 * First working version of ncm_fit_esmcmc + MPI. Example pysimple updated.

 * First working version of nc_hiqg_1d.h (some speed-ups are still necessary).

 * Missing object files.

 * Renaming quantum gravity object.

 * Added a front-end for other lapack functions. Reorganized the blas header
     inclusion. Working in progress in ncm_qm_prop, last commit before removing
     different approaches code.

 * New examples and working in progress for NcmFitESMCMC and NcmMPIJobMCMC.

 * Finished support for complex messages in NcmMPIJob. Two implementations tested
     NcmMPIJobTest and NcmMPIJobFit. Working on NcmFitESMCMC parallelization
     using MPI.

 * Documentation.

 * Finished the first version of MPI support, including the helper objects
     NcmMPIJob*.

 * Trying different interpolation methods in ncm_qm_prop.c.

 * Workaround travis-ci problem.

 * Fixed test test_nc_wl_surface_mass_density.c.

 * Removed leftover headers in ncm_spline_rbf.c.

 * Updated example example_wl_surface_mass_density.py.

 * First working version of nc_data_reduced_shear_cluster_mass. Improvements in
     all related objects.

 * Added support (optional) to HDF5. New objects NcGalaxyRedshift,
     NcGalaxyRedshiftSpec and NcGalaxyRedshiftSpline to describe galaxy redshift
     distributions. Included support for loading hdf5 catalog in
     nc_data_reduced_shear_cluster_mass.

 * Fixed bug: ncm_util_position_angle was returning -(Pi/2 - theta). Corrected to
     return theta.

 * Removed printf in great_circle_distance function.

 * Implemented the position_angle and great_circle_distances functions.

 * Included new lapack encapsulation functions. Trying different methods in
     NcmQMProp. New 1D interpolation object NcmSplineRBF.

 * Removed unnecessary range check in gobject parameter properties. Working on
     ncm_qm_prop.

 * Added function eval_full to nc_xcor_limber_kernel.

 * Added missing conditional compilation of MPI support.

 * First tests of MPI slaves.

 * Testing different autoconf mpi detections.


[v0.14.1]
 * Removed deploy and less verbosity on make.

 * Removed extra header.

 * New version 0.14.1.

 * Included link for parallel linear solvers in sundials. Fixed bug in nc_cbe
     (lensed CMB requirements).

 * Updated examples to python3. New object nc_galaxy_selfunc. Working on
     ncm_qm_prop.

 * Fixed possible (impossible in practice) overflow in background.c. Added support
     for binder.

 * Fixed several documentation glitches.

 * Finished support for Planck lensing likelihood.

 * Created data object NcDataReducedShearClusterMass. Work on progress.

 * Fixed bug in ncm_spline.h. Working on ncm_qm_prop.

 * Included properties. Work in progress.

 * Improved regex that greps SUNDIALS_VERSION in configure.ac.

 * Working on example_qm.c, removed old file from numcosmo-docs.sgml.in, added
     quotes to SUNDIALS_VERSION grep in configure.ac.

 * Implementing object to estimate mass from reduced shear. Work in progress.

 * Fixed data install path. Tweaked fit tests.


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


[]
 * V0.23.0 (#187)
 * Notebook to generate the plots for the notaknot paper (cosmology sess… (#169)
 * Add NcmSphereNN for finding nearest neighbors within a spherical shell. (#186)
 * Raising error for unknown key in param_set_desc. (#184)
 * Adding support for MC analysis. (#183)
 * Mass and concentration summary  (#180)
 * Configuring conda-incubator/setup-miniconda@v3.
 * Updating conda-incubator/setup-miniconda@v3 usage.
 * Updated conda-incubator/setup-miniconda@v3 use.
 * Adding support for version checks in numcosmo. (#179)
 * Testing more parallel tests.
 * Fix leftover merge lines.
 * Fftw config (#178)
 * Improving fftw planner control.
 * Configuring fftw-planner during build.
 * Using FFTW_ESTIMATE by default. Added NC_FFTW_DEFAULT_FLAGS and
     NC_FFTW_TIMELIMIT environment variables.
 * Forcing cache update.
 * Removing use-only-tar-bz2: true.
 * Adding use-only-tar-bz2: true to miniforge action.
 * Updated GHA workflow.
 * Removed old coveralls badge.
 * Twofluids update (#177)
 * Updating tests use of Vexp, fixing documentation bugs. (#176)
 * Magnetic vexp (#175)
 * New nc_galaxy_wl_obs object  (#167)
 * Updating to actions/upload-artifact@v4.
 * Restricting setuptools version to avoid gobject-instrospection problems.
 * Improving model interface and error handling (#174)
 * Updated documentation of ncm_m_mass_solar. CODATA 2022.
 * Updated to latest CODATA, NIST and IAU (and others) constants. (#171)
 * Xcor cmp (#85)
 * Xcor CCL comparisons (#168)
 * Implemented Integrated Sachs-Wolfe kernel. (#166)
 * Galaxy WL reformulation (#93)
 * Magnetic Fields in Vexp cosmology (#153)
 * New version v0.22.0
 * Mix experiments options (#164)
 * CCL background power (#162)
 * Two Fluids primordial model (#160)
 * Support for multiple corner plots. (#161)
 * Added support for parameter filtering in numcosmo app. (#159)
 * Adding support for 1d distributions. (#158)
 * Updating uncrustify configuration. (#157)
 * Calibrating fftlog tests number of knots.
 * Fixing spline constructors to return the right type. (#156)
 * Removing glib version restriction (#155)
 * NcHIPert Reformulation (#95)
 * Adding Bayesian evidence support for numcosmo app. (#152)
 * Sample variance (#107)
 * Updating codecov to v4. (#151)
 * Updating requirement versions. (#150)
 * Fixed instrospection error for gobject-instrospection >= 1.80. (#149)
 * Cmb parametrization (#148)
 * Planck data analysis reorganization (#145)
 * Adding stub generation script to gitignore.
 * Updated stubs.
 * New minor release v0.21.2
 * Adding MPICH support. (#144)
 * NumCosmo product file (#143)
 * Moving release v0.21.1.
 * More options to the conversion tool from-cosmosis. (#141)
 * New minor release v0.21.1.
 * Updates and tests for NumCosmo app (#140)
 * New bug fix release.
 * Fixed mypy issues.
 * Ran black.
 * Fixed cosmosis required parameters issue due to returning iterator. Fixed
     restart issue on numcosmo run fit.
 * New release v0.21.0
 * numcosmo command line tool (#137)
 * Variant dictionary support  (#135)
 * Support for object dictionaries, NcmObjDictStr and NcmObjDictInt. (#134)
 * Better python executable finding.
 * Added GSL as a dependency for libmisc (internal library). (#133)
 * Updated stubs.
 * New minor version v0.20.0
 * Support for computing fisher bias vector (#132)
 * Improving tests for NumCosmoMath (#131)
 * Adding support for libflint arb usage. (#130)
 * Adding more python based tests using external libs (astropy and scipy). (#129)
 * Fixed package name in pyproject.toml.
 * Including typing data into pyproject.toml. Updating changelog.
 * Fixing minor doc glitches. (#128)
 * Updated changelog.
 * Fixed project name in pyproject.toml.
 * New minor version.
 * Using pip to install python modules. (#127)
 * More objects encapsulation (#126)
 * New minor release 0.19.1
 * Updated meson to deal with cross compiling and GI building. Updated ncm.pyi.
 * Removed git ignored files related to autotools and in-source building.
 * Removed unnecessary packages.
 * Yaml implementation (#125)
 * Adding fyaml to CI.
 * Updated Python stubs.
 * Complete version of the yaml serialization, including special types.
 * First version of from_yaml and to_yaml serialization. Updated minimum glib
     version.
 * New tuple boxed type (#124)
 * Objects encapsulation (#122)
 * Removed unnecessary header inclusions to avoid propagating depedencies.
 * Fixing warnings in conda build. (#121)
 * Mypy to ignore python scripts inside meson builds.
 * New test for simple vector set/get.
 * Removed old files.
 * Release v0.19.0.
 * Testing before adding warn supp. Testing for isfinite declaration.
 * Adding cfitstio to plc.
 * Addind examples to installation.
 * Added libdl dependency to plc.
 * Added GSL blas definition to avoid double typedefs.
 * Updated changelog.
 * Updated stubs.
 * Moving to meson (#120)
 * Minor improvements. (#119)
 * Added correct prefix for NCM_FIT_GRAD.
 * Kde loocv (#118)
 * Creating tests and documentation for n-dimensional integration object (#108)
 * Improving documentation and encapsulating objects (#116)
 * Now MPI jobs do not require setting nthreads. (#115)
 * 109 example describing 3d correlation (#113)
 * Added pocoMC to rosenbrock_simple.ipynb. (#112)
 * Updated rosenbrock_simple.ipynb.
 * New version v0.18.2.
 * Improving stubs.
 * New minor version v0.18.1.
 * Missing files for python typing
 * Updated changelog.
 * New minor release 0.18.0.
 * Implementing n-dimensional integration object (#106)
 * Create SECURITY.md
 * Create CONTRIBUTING.md
 * Create CODE_OF_CONDUCT.md
 * Update issue templates
 * Update bug_report.md
 * Update issue templates (#104)
 * 102 add notebook for gauss constraint tests (#103)
 * Added TMVN sampler. (#101)
 * Update README.md
 * Update README.md (#99)
 * Updated python interface. (#98)
 * Fixed bug that resets the values of use_threads in APES.
 * Several improvements on APES and others. (#92)
 * Removed printf from test.
 * Improving parallelization for APES.
 * Fixed setting of max_ess.
 * Added conditional compilation of internal function.
 * Added support in ncm_mset_catalog and mcat_analyze to compute acceptance ratio.
 * Many minor improvements.
 * Fixed scripts shebang.
 * New version 0.17.0.
 * New experimental python interface for sampling. New sampler comparisons. (#90)
 * Encapsulating objects (#72)
 * New features (#81)
 * WL binned likelihood object (#77)
 * Added check to see if the python interface is available.
 * Improved tests.
 * New method for likelihood utilizing KDE (#65)
 * Reordering -I to include first internal sub-packages.
 * Added conditional use of g_tree_remove_all. Removed setting of all threads to
     one. Reintroduced non fatal assertions in test_ncm_stats_dist.
 * Improved ax_code_coverage.m4 to work with newer lcov versions.
 * Organized m4 files and fixed lcov issues.
 * Fixed a few lcov issues.
 * Incresead number of points when testing StatsDist with rubust-diag.
 * Changed divisions to multiplications.
 * Fix bug in AR fitting when only two elements were available.
 * Improved tests, added test to robust covariance computation.
 * Modified vkde to use block triangular system solver.
 * Finished the refactor of kdtree to use a red-black tree and prune impossible
     branches. Significant increase in speed!
 * Removed old tree. Finishing prunning.
 * New red-black BT to improve kdtree. Added support for prunning kdtree during
     search.
 * 75 organizing python modules (#76)
 * Updated kdtree and directories in notebooks/Makefile.am.
 * Reorganizing notebooks.
 * Added missing cell.
 * Fixed missing properties (unused). autogen.
 * Notebooks massfunc (#74)
 * Halo bias tests (#73)
 * Mean bias (#61)
 * Added interface to generate models using an array of NcmSParams
 * Added two missing files to the releases.
 * New minor release v0.16.0.
 * 40 numcosmo unit test coverage (#68)
 * Added new method to set model parameter fit types to their default values.
 * Missing semicolon.
 * Minor changes on NcDistance initialization order.
 * Updated gcc version for macos ci.
 * Better debug messages in GHA.
 * Added more robust testing for power-spectra.
 * Updating e-mails.
 * Updated e-mail in copyright notices.
 * uncrustify code.
 * Minor fixes in documentation. Finished coverage and tests for special
     functions. Removed old code.
 * Added support for namespace search in ncm_mset_func_list. Added plot_corner
     helper script.
 * More cleaning and adding more files to .gitignore.
 * Cleaning autotools files and old unused tools. (#67)
 * Removed old and unused code.
 * Uncrustify sources.
 * Reorganized all ncm_spline2d objects and improved unit testing and coverage.
 * Uncrustify ncm_spline2d_bicubic.
 * Improved coverage of NcmDiff.
 * Uncrustify ncm_diff.c.
 * 60 statsdist1d error (#63)
 * Removed old code causing a bug in ncm_stats_dist1d_epdf_reset.  (#62)
 * Fixed unimportant warnings in class.
 * Removed debug message.
 * Minor fixes in twofluids framework. Updating StatsDist to use only a fraction
     of the sample when computing the bandwidth using a split cross-validation.
 * Multiplicity watson (#59)
 * Removed ckern algo.
 * Debug version, do not use it. Version containing the constant kernel option in
     NcmStatsDist.
 * Several minor improvements.
 * Removed old CLAPACK and LAPACKE support.
 * W reconstruction (#58)
 * Trying to find correct path due to broken glib in brew.
 * Debugging missing prereq.
 * New dependency resulting from split package in homebrew.
 * 56 lastest pantheon (#57)
 * Minor updates in the figures of VacuumStudy and VacuumStudyAdiabatic.
 * Added volume method to nc_cluster_mass_nodist.
 * Adding Minkowski functions to CSQ1D.
 * Cleaning notebooks.
 * Added constructor annotation to ncm_mset_load(). New Vacuum study notebooks.
 * Added method to get the best fit from catalogs.
 * Included more frames for CSQ1D
 * Testing coverage tweaking.
 * Halo bias (#53)
 * Updated autotools.
 * Implementing frames in csq1d.
 * More tests for nc_data_cluster_ncount.c.
 * Removed option to print the mass function (old code).
 * Removed old method nc_data_cluster_ncount_print.
 * More tests for test_nc_data_cluster_ncount.c.
 * Removed inclusion of removed objects documentation.
 * Adding new integration routines to the ignore list in docs.
 * New integration code. Now vector integration used in nc_data_cluster_ncount.
     Fixed bug in NcmFitMC (it was using the bestfit from catalog instead of
     fiducial model to resample). Fixed typos.
 * Added a full corner plot comparing all outputs.
 * Updated generate_corner.ipynb to use ChainConsumer.
 * Fixed bug in catalog_load nc_data_cluster_ncount. New corner plot notebook.
 * Fixed minor leaks in ncm_reparam.c ncm_powspec_filter.c ncm_mset_catalog.c.
     Improved sampling in ncm_fit_esmcmc_walker_apes (now the second half use
     the updated first half when moving the walkers). Support for binning in
     nc_data_cluster_ncount. New notebooks comparing binning vs unbinning.
 * Inclusion of the time function to compare the effiency between CCL and Numcosmo
 * Reorganized binning options in NcDataClusterNCount.
 * Unbinned and binned analisys in the ascaso proxy
 * Reorganizing cluster mass ascaso object.
 * CCL- Numcosmo comparison using a mass proxy, both binned and unbinned analysis
 * Tests with de cluster abundance with a mass proxy
 * Proxy comparation
 * Fixed conflict leftovers.
 * Ascasp changes
 * Removed old data objects all binned versions now reside in NcDataNCount.
     NcABCClusterNCount needs updating. Now lenghts of cluster mass and redshift
     and class properties. Cluster abundance must be instantiated with both mass
     and redshift proxies defined. NcClusterMass/Redshift objects reorganized.
 * New helpers scripts with new tools: a function create pairs of NumCosmo/CCL
     objects, increase CCL precision and notebook plots with comparison between
     NumCosmo and CCL outputs. Updated notebooks to use helper functions.
 * Inclusion of the Cluster Number as a function of mass in the binned case both
     for CosmoSim and Numcosmo
 * Implementation of the inp_bin and p_bin_limits function in the
     gauss_global_photoz redshift proxy
 * Removed checkpoints and output files.
 * Implementation of binning in the lnnormal mass-observable relation
 * Binned and unbinned comparison between Numcosmo and CCL cluster abundace
     objects with no mass or redshift proxies
 * binned and unbinned comparison between CCL and Numcosmo cluster abundance with
     no mass or redshift proxies
 * Working version for binning proxies in NcCluster* family.
 * notebook on cluster mass comparison between CCL and Numcosmo update
 * addition of  binning in cluster_mass.c and cluster_mass.h and unbinning
     comparison between CCL and Numcosmo cluster mass objects(not ready yet)
 * Old modifications on hiqg and updates on NumCosmo vs CCL tests. Starting the
     implementation of binning for cluster mass and redshift.
 * Comparison between numcosmo and ccl cluster abundance objects
 * New spline object for functions with known second derivative. Updated
     nc_multiplicity_func_tinker to use interpolation objects, added option to
     use linear interpolation. Removed old notebook NC_CCL_Bocquet_Test2.ipynb.
     Updated NC_CCL_mass_function.ipynb (fixed bugs).
 * Mass functions comparisons notebook.
 * Updated version to match new interface.
 * Better limits for nc_halo_mass_function. Setting properties through gobject to
     catch out-of-bounds values.
 * Adjusted esmcmc run_lre minimum runs in tests.
 * Calibrated integrals to work on any point of the allowed parametric space.
     Added mores tests.
 * Modified ranges of concentration and alpha (Einasto) parameters.
 * Improved stability in nc_halo_density_profile.c.
 * Smaller lower bounds for ncm_fit_esmcmc_run_lre. Added option for starting
     value of over-smooth in mcat_analize calibration.
 * Option to calibrate over-smooth.
 * New minor version 0.15.4.
 * Added missing ncm_cfg_register_obj call.
 * Delete NC_CCL_Bocquet_Test-checkpoint.ipynb
 * Delete .project
 * test of execution time
 * updates
 * New option to use kde instead of interpolation in APES.
 * Notebooks testing.
 * Moved model validating to workers (slaves or threads).
 * Improved fparam set methods.
 * Faster kde sampling.
 * Improved MPI debug messages added timming.
 * Improved control thread avoinding aggressive pooling by MPI.
 * Added conditional compilation of MPI dependent code.
 * Using switch to choose between kernel types.
 * Fixed memory leak.
 * the hydro and dm functions of the CCL were included
 * Finalized tests for kernel class.
 * Better handling of the case where 0 threads are allowed. Fixed limits on
     nc_cluster_photoz_gauss_global. Incresed lower limit in As in
     nc_hiprim_power_law.
 * Fixed leak.
 * New mpi run jobs async (master - slaves).
 * Added test for the kernel sample function.
 * Implentationg of tests for the #NcmStatsDistKernel class.
 * Configuration.
 * More tweaks on omegab range.
 * New notebook NC_CCL_Bocquet_Test has been created
 * Increasing lower limit of Omega_bh^2.
 * Fixed variable types for simulation (sim).
 * Improved bounds on nc_hicosmo_de_reparam_cmb.
 * Clean up Tinker: no need to set some parameters as properties. Delta -
     CONSTRUCT and not CONSTRUCTED_ONLY
 * Fixed bug in Bocquet multiplicity function, e.g., properties are CONSTRUCT not
     CONSTRUCT_ONLY.
 * Fixing mcat_analize to work with small catalogs.
 * Resolved conflict.
 * Fixed minor warnings.
 * Implemented Bocquet et al. 2016 multiplicity function. Two new functions in
     NcMultiplicityFunc: has_correction_factor and correction_factor. Bocquet
     provides fits for mean and critical mass definitions, but the latter
     depends on the first.
 * Fixed some edge cases in ncm_fit_esmcmc.c and ncm_fit_esmcmc_walker_apes.c.
     Minor reorganization.
 * Add files via upload
 * Testing 10D.
 * Added missing object registry.
 * Removed incomplete tests.
 * Fixed allocation problem in ncm_stats_dist.c. Fixed other minor bugs and
     tweaks.
 * Included tests for the error messages in stats_dist_kernel.c
 * Test if kernel test is implemented right.
 * Finished the documentation of ncm_fit_esmcmc_walker.c and
     ncm_fit_esmcmc_walker_apes.c
 * Implemented documentation of ncm_fit_esmcmc_Walker.c
 * Updated automake file.
 * uncrustify and more tweaks on test_ncm_diff.c removing edge cases.
 * Fixed internal struct access.
 * uncrustify.
 * Updated test, and fixed minor issues.
 * Fixed docs and set nc_multiplicity_func.c to abstract.
 * Refactoring of the multiplicity function object is complete. Main difference:
     included mass definition as a property. Examples were properly updated.
 * Tweaked test_ncm_mset_catalog.c and test_ncm_diff.c. Solved APES offboard
     sampling.
 * Changed the size of figures in docs and improved the documentation of
     StatsDistKernel objects.
 * Improving coverage and fixed casting.
 * Generating graphs with the notebooks.
 * Included over-smooth option in APES. Added the same option to darkenergy's
     command line interface. Improved documentation and coverage.
 * Improved tests and coverage for NcmStatsDist* family.
 * Improving NcmStatsDist* coverage.
 * More tweaks on NcmDiff tests.
 * Tweaking tests to avoid false positives.
 * Improved unit tests for NcCBE, NcCBEPrecision and NcmVector.
 * Improved interface to NcmFitESMCMCWalkerAPES. Included and tweaked unit tests.
 * I am rewriting the multiplicity function objects. Including missing functions
     (e.g., ref, free, clear...), put in the correct order. Add "mass
     definition" as a property.
 * Fix documentation glitches and solve warnings.
 * Documentation for stats dist objects with image problems
 * Unfinished stats dist objects documentation
 * Removed whitespace following trailing backslash.
 * Added missing include directory.
 * Working on stats_dist.c documentation
 * Uncrustify tests. Tweak mcmc tests.
 * Fixed bug in ncm_mset_trans_kern_cat.c (re-preparing for each sampling). Added
     missing files. Added new test to test_ncm_vector.c. Tweaking tests.
 * uncrustify and rename APS to APES.
 * Fixed wrong href when computing IM in VKDE. Fixed over_smooth tweak in
     prepare_interp.
 * Removed unecessary files. Added notebooks.
 * Working version of ncm_stats_dist*. Not yet fully tested.
 * First (incomplete) reorganized version of NcmStatsDist*. Updated mkenums
     templates.
 * Updated notebook. Halo profile uses log10(M) instead of M. Modifying
     Multiplicity function objets: mass definition is a property. Work in
     progress.
 * Removed CNearTree.
 * Working version of vbk.
 * New script to use numcosmo without installing.
 * Testing for fit with no free parameters bug. Fixed the same bug in fit impls.
 * Fixed indentation.
 * Removed unecessary files.
 * vbk_studentt working on notebook. Memory error for rosenbrock. Check slack for
     info.
 * Adding support for non-adiabatic computation.
 * vbk_studentt working for eval and evan_m2lnp. Copy of APS to work with vbk (not
     included in makefile). Copy of gauss to gauss vbk(included in makefile)
 * Missing files from last commit.
 * Functions prepare_args and preapre_interp running. Starting to work on
     eval_m2lnp. Interp.py is the test file.
 * Added more precise delta_c.
 * Working on the examples.
 * Working on VBK.
 * Updated autotools file and removed binnary.
 * example_neartree is the example from documentation, test_neartree is build by
     me and slightly documented.
 * Fixed the includes for CNearTree, inserted a flag in Makefile.am and created a
     test to check.
 * Added gtk-doc to mac os build.
 * removed azure.
 * removed azure.
 * Removed travis-ci.
 * Updated autotools files and removed travis-ci.
 * Added the required files for CNearTree.c library, created copies of stats dist
     to work on, and added the necessary lines in the makefiles.
 * Added gtk-doc to mac os build.
 * Funnel example and notebook.
 * New test likelihood Funnel.
 * Included the RoT for the Student t distributions in
     ncm_stats_dist_nd_kde_studentt (truncated for nu < 3.0 since it is not
     defined for these values).
 * Set default to aps with studentt (Cauchy dist) with no CV and over smooth 1.5.
 * Fixed bug in ncm_stats_dist_nd (it didn't set weights vector to zero before
     fitting).
 * New notebook used to plot Rosenbrock MCMC evolution.
 * New Rosenbrock model/likelihood to check MCMC convergence. New option to thin
     chains. New example to run Rosenbrock MCMC.
 * Typo fix.
 * Working version of the reorganized code (NcmNNLS, NcmISet and NcmStatsDistNd).
 * Working on ncm_stats_dist_nd + ncm_nnls. Working version, finishing code
     reorganization.
 * Working version (not organized yet, full of debug prints...).
 * Added documentation and comentaries in the  notebook TestInterp.ipnb.
 * Improved the description in the documentation of ncm_stats_dist_nd.c,
     ncm_stats_dist_nd_studentt.c and ncm_stats_dist_nd_gauss.c.
 * Moved headers to the right place.
 * Missing Makefile.am.
 * Moved external codes to a new (sub)library to remove these codes from the
     coverage and to make the symbols invisible.
 * Minor release v0.15.3.
 * Added interpolation case where only the most probable point is necessary.
 * Added tests for KDEStudentt.
 * Fixed a few documentation glitches.
 * Reorganized ncm_stats_dist_nd* objects family. Testing different solvers to the
     NNLS problem.
 * Included a function to compute numerical integrals of the NFW profile (instead
     of the analytical forms). To be used for testing only!
 * Change on the file numcosmo-docs.sgml to include
     ncm_stats_dist_nd_kde_studentt.c. Did not create a studentt HTML as I
     expected.
 * Reupdated m4 and automake stuff.
 * Minor identation/positional tweaks.
 * Implementation of the comentaries from the commit "New implementation of
     studentt function for ncm_stats_nd_kde.".
 * (Re)updated m4 macros and gtk-doc.make.
 * New methods to access Ym values in NcmFftlog.
 * Adding the updated Jupyter notebook
 * New implementation of studentt function for ncm_stats_nd_kde.
 * Added a second run to avoid unfinished minimization process.
 * Removed debug msgs from coverage build.
 * Removed coverage flags from introspection build.
 * Debug coverage build.
 * Debug coverage build.
 * Debug coverage build.
 * Removed LDFLAGS for coverage.
 * Debug coveralls build.
 * Moved (all) flags to the right places.
 * Moved flags to the right place.
 * Added explict CODE_COVERAGE_LIBS to introspection build.
 * Debug coveralls build.
 * Test speedups.
 * Allowed reasonable failures.
 * Added 10% allowed test errors when estimating hessian computation error.
 * Testing ncm_stats_dist_nd_kde_gauss.c. Minor modifications to
     ncm_data_gauss_cov_mvnd.c. New notebook to test multidimensional
     interpolation.
 * Debug mac-os GHA
 * Debug mac-os GHA
 * Debug mac-os GHA
 * Debug mac-os GHA
 * Debug mac-os GHA
 * Debug mac-os GHA.
 * Debug mac-os GHA.
 * Debug mac-os GHA build.
 * Trying reinstalling gmp.
 * Testing a solution for GHA on mac-os.
 * Still debugging macos build in GHA.
 * Debug macos build.
 * Conditional use of sincos.
 * Fixed sincos warning.
 * More compiler env.
 * Fixed sincos included warning.
 * Updated example.
 * Setting compilers.
 * Cask install for gfortran in macos build.
 * Testing lib dir in GHA.
 * Added cask install fortran for macos build.
 * Trying lib dirs.
 * Added gfortran req to macos build.
 * Added prefix option to configure in GHA.
 * Fixed example name and moved test.
 * Rolled back autoconf version req.
 * Included missing make install in build check.
 * Updated autotools and deps. New check in GHA. Fixed bug in numcosmo.pc.in.
 * Working on ncm_csq1d.c. New notebook FisherMatrixExample.ipynb.
 * Adding timezone info.
 * New docker image with NumCosmo prereqs.
 * Working on nc_de_cont.
 * Running actions in every branch.
 * Testing GHA
 * Testing GHA
 * Testing GHA
 * Testing GHA
 * Testing coveralls build.
 * Updated CI badge to GHA.
 * Adding missing prereq for the macos build.
 * Better workflow name and removed unnecessary prereq in the macos build.
 * Adding macos build.
 * Removed debug print in c-cpp.yml.
 * Adding references.xml to the repo.
 * Update c-cpp.yml
 * Update c-cpp.yml
 * Fixed doc typo.
 * Updated to sundials 5.5.0.
 * Added NumCosmo CCL test notebook.
 * Fixed conditional compilation for system with gsl < 2.4.
 * Minor release 0.15.2.
 * Updated tests and fixed indentation.
 * New framework for Cluster fitting with WL data (in progress).
 * New minor version. Reorganizing WL likelihood (in progress).
 * Default refine set to 1.
 * More options to refine.
 * Add refine as an option.
 * Added vectorized interface for nc_wl_surface_mass_density_reduced_shear. Minor
     other improvements.
 * Improvement in ncm_spline_func to remove outliers.
 * Update c-cpp.yml
 * Update c-cpp.yml
 * Update c-cpp.yml
 * Update c-cpp.yml
 * Update c-cpp.yml
 * Update c-cpp.yml
 * Update c-cpp.yml
 * Update c-cpp.yml
 * Update c-cpp.yml
 * Update c-cpp.yml
 * Update c-cpp.yml
 * Update c-cpp.yml
 * Updated and finished support for CCL in Dockerfile-clmm-jupyter.
 * Added python3-yaml support.
 * Added support for CAMB and CCL.
 * Support for CCL and CAMB.
 * Adding support for camb and ccl.
 * Fixed new filename.
 * Notebook comparing Colossus and CCL with NumCosmo: density profiles, surface
     mass density and the excess smd.
 * Added hook between nc_halo_mass_function and ncm_powspec_filter to ensure the
     same redshift range.
 * Fixed minor leak.
 * Minor fixes and improvements in ncm_spline_func_test.*.
 * Update c-cpp.yml
 * Create c-cpp.yml
 * Included new distance functions from z1 to z2.
 * Test suit NcmSplineFuncTest is now stable enough. Next step: add some
     cosmological functions examples.
 * Doc. minor changes.
 * Memory leak - ncm_spline_new_function_4
 * Added new outlier function to last description example.
 * Corrected vector memory lost.
 * Added option to save outliers grid to further analysis.
 * Added test suite to NcmSplineFunc.
 * Tutorial reviewed.
 * Text review (Mari).
 * Added missing ipywidgets from docker build.
 * Fixed makefiles.
 * Reorganized and added copyright notices to notebooks.
 * New tutorial.
 * New code for homogeneous knots.
 * ncm_vector.c documentation improve.
 * Doc. improvement.
 * Minor doc. modification.
 * Minor doc. changes.
 * Added support for abstol in NcmSplineFunc.
 * Set max order to 3 in NcmODESpline to make the ode integration tolerance agree
     with spline interpolation error.
 * Removed old CCL interface (they no longer have a C api, we are moving to test
     in python since their API is only there).
 * Improve doc. & fixed indentation:
 * Transfer function improve doc.
 * NcTransferFuncEH: fixed indentation.
 * NcTransferFuncEH: improve doc.
 * NcTransferFuncBBKS: corrected minor typo.
 * NcTransferFuncBBKS: fixed indentation.
 * NcTransferFuncBBKS: improve doc. Add BBKS ref.
 * NcTransferFunc: fixed indentation.
 * NcTransferFunc: improve doc.
 * NcWindowGaussian & NcWindowTophat: standardization between both descriptions.
 * NcWindowGaussian: fixed indentation.
 * NcWindowGaussian: improve doc.
 * NcWindow: fixed indentation.
 * NcWindow: improve doc.
 * NcWindowTophat: fixed indentation.
 * NcWindowTophat: improve doc.
 * More debug messages in MPI.
 * Better debug messages and identation.
 * Documentation.
 * Fixed details in the documentation.
 * NcmPowspecFilter: fixed indentation.
 * NcmPowspecFilter: doc. improvement.
 * NcmPowspec: reference to function NcmPowspecFilter in ncm_powspec_var_tophat_R
     ()
 * NcmPowspec: Fixed indentation.
 * NcmPowspec: doc. improvement.
 * Minor typo.
 * NcmODEEval: Fixed indentation.
 * NcmODEEval: doc. improvement.
 * NcmODE fixed indentation.
 * NcmODE doc. improvement.
 * NcmSpline2dBicubic: fixed indentation and tweak doc.
 * Fixed minor typos.
 * Fixed indentation:
 * NcmSpline2dSpline and NcmSpline2dGsl doc tweaks.
 * NcmSpline2d: fixed indentation.
 * NcmSpline2d: Added Include and Stable tags + Minor tweaks.
 * Fixed indentation:
 * NcmFftlogTophatwin2 and NcmFftlogGausswin2: doc. improvement.
 * Fixed indentation: ncm_powspec_corr3d.c/h.
 * NcmPowspecCorr3d: doc. improvement.
 * NcmFftlogSBesselJ: fixed description and minor tweaks.
 * NcmFftlogSBesselJ: tiny tweaks in the description.
 * Fixed indentantion.
 * NcmFftlogSBesselJ: corrected indentation.
 * NcmFftlogSBesselJ: documentation improved.
 * Fixed wrong lower bound for abstol.
 * Removed old test in autogen.sh and overwritting of gtk-doc.make.
 * Add gtkdoc related files (instead of soft links).
 * Added to repo all necessary m4 files.
 * NcmGrowthFunc: doc tiny tweaks
 * NcmFftlog: corrected indentation.
 * NcmFftlog: documentation's minor improvement.
 * Tweaking NcGrowthFunc documentation and fixed wrong link for NcmSplineFunc.
 * NcGrowthFunc: changed description to a vague explanation on the initial
     conditions. Added Martinez and Saar book on the references.
 * NcGrowthFunc: corrected indentation.
 * NcGrowthFunc: improved documentation.
 * NcDistance: Standardization of function documentation
 * NcDistance: modified two static functions names:
 * NcDistance: corrected indentation.
 * NcmDistance: improved documentation.
 * Testing support for gcov.
 * ncm_timer.* - correct indentation with uncrustify.
 * NcmTimer: improved documentation.
 * Corrected a broken link in short description.
 * Indentation using uncrustify.
 * Improved documentation from NcmRNG.
 * Changed "abs" --> "abstol" in ncm_ode_spline_class_init. Also some minor
     changes.
 * Improved #NcmOdeSpline documentation.
 * Updated private instance get function. Fixed doc issues.
 * Better support for arb.
 * New example.
 * Using different branch in CLMM.
 * Added colossus to Dockerfile-clmm-jupyter.
 * Version 0.15.0
 * Polishing nc_halo_density_profile. Added NumCosmo x Colossus comparison
     notebook.
 * Fixed missing parameter doc.
 * Updated test test_nc_wl_surface_mass_density.
 * Better integration strategy for NcHaloDensityProfile. Updated test
     test_nc_halo_density_profile.
 * Final tweaks before release.
 * Minor improvements in notebooks/BounceVecPert.ipynb.
 * Added new profile (Hernquist), Einasto implementation is now complete. Added
     documentation.
 * Corrected indentation of ncm_ode_spline.*
 * Improved notebooks/BounceVecPert.ipynb.
 * [DOC] Improved main and enum description.
 * Removed punctuation from parameter descriptions:
 * Removed punctuation from parameter descriptions in ncm_matrix.h/c.
 * Removed punctuation from parameter descriptions. Added bindable function to
     NcmSplineFunc.
 * Added description/documentation to ncm_spline_func.h/c
 * Fixed minor bugs.
 * Added part of doc from spline_func module.
 * Removed old nlopt header in csq1d.
 * Fixed indentation and conditional load of NLOPT library object.
 * Implemented Einasto profile (just rho, not the integrals). Included the
     funciton to compute the magnification.
 * New refactored NcHaloDensityProfile (working in progress). Added support for
     different internal checkpoints in NcmModel.
 * Included "@stability: Unstable" and "@include: numcosmo/math/ncm_spline_rbf.h".
 * Added NcmSplineGslType enum description.
 * Added "@stability: Stable" and "@include:
     numcosmo/math/ncm_spline_cubic_notaknot.h".
 * Corrected indentation of ncm_spline_cubic.c and ncm_spline_cubic.h with
     uncrustify.
 * Added doc to functions:    * ncm_spline_is_empty    * ncm_spline_class_init
     (g_object_class_install_property)
 * Added doc in functions ncm_vector_class_init & ncm_matrix_class_init.
 * Some minor tweaks:
 * Passed ncm_matrix.h and ncm_matrix.c through uncrustify to set indentation.
 * Added doc to the following functions of NcmMatrix:
 * Added numcosmo's uncrustiify settings.
 * Uniform indentation.
 * Uniform indentation.
 * Homogenization and standardization of the #NcmVector module.
 * Uniform indentation.
 * Finished first version of #NcmVector documentation. Still needs a careful
     check.
 * Renamed NcDensity* objects to NcHaloDensity*.
 * Added "@stability: Stable" and "@include: numcosmo/math/ncm_c.h" to section in
     ncm_c.c.
 * Fixed use of Planck likelihood without check_param. Fixed typo in
     numcosmo/math/ncm_c.c.
 * 1) Changed function name: ncm_c_hubble_cte_planck_base_2018 to
     ncm_c_hubble_cte_planck6_base.
 * Now the last commit is correct.
 *    * ncm_c_crit_density_h2    * ncm_c_crit_mass_density_h2
 * 1) Deleted function from #NcmC:    * ncm_c_hubble_cte_msa - it was not applied
     anywhere.
 * 1) Changed documentation to the following #NcmC module functions:    *
     ncm_c_wmap5_coadded_I_K    * ncm_c_wmap5_coadded_I_Ka    *
     ncm_c_hubble_cte_hst
 * Added documentation to the following #NcmC modules functions:
 * Documented function ncm_vector_len.
 * Better expansion for tan(x+d)-tan(x) for small d.
 * Fixed corner case in numcosmo/model/nc_hiprim_atan.c.
 * Fixed MPI in hdf5 incompatibility.
 * Removed update option on homebrew.
 * Fixed error handling, clik returns wrong values when an error occurs (due to a
     wrong usage of forwardError), to fix this we changed the likelihood to
     return m2lnL = 1.0e10 whenever clik returns an error.
 * Workaround to fix travis ci bundle issue.
 * Removed wrong free in CLIK_CHECK_ERROR.
 * Updated old m4 files and building system to keep them updated in the m4/
     folder.
 * New MPI server (in progress). NcDataPlanckLKL no longer kills process when clik
     returns an error (just sends a warning).
 * Better OpenMP (and others) number of threads control.
 * Finished update to 2018 Planck likelihood (in testing).
 * Updates in the notebook.
 * Fixed sprintf related warnings.
 * Removed typo.
 * Fixing new docs build process...
 * Fixing error messages in plik. Working on doc building process.
 * Fixing new docs build process...
 * Reorganized docs building process.
 * Testing travis-ci on osx.
 * Testing travis-ci on osx.
 * Testing travis-ci on osx.
 * Testing travis osx.
 * Fixed doc typos. Testing travis on osx.
 * Fixed doc in NcDensityProfile. Testing travis ci on osx.
 * Updated sundials to version 5.1.0. Fixed somes tests and updated to TAP.
 * Update MagDustBounce.ipynb from Emmanuel Frion.
 * Minor updates and annotation improvements.
 * Copy all examples in Dockerfile-clmm-jupyter.
 * Updated notebooks and Dockerfile-clmm-jupyter.
 * New parametrization for CSQ1D, updates in BounceVecPert.ipynb and
     MagDustBounce.ipynb.
 * Several improvements in MagDustBounce.ipynb.
 * New NcHICosmoQRBF model.
 * Added options to logger function to redirect all library logs.
 * Testing CLMM+NumCosmo notebook.
 * Removed debug print.
 * Fixed documentation error.
 * Added missing scipy for CLMM.
 * Added missing Astropy for CLMM.
 * Cloning the right branch from CLMM.
 * Copying examples from CLMM to work.
 * Added COPY from opt.
 * Fixed build script name...
 * Testing different build order.
 * Testing build from git.
 * New Dockerfile for CLMM comparison.
 * Updated to CODATA 2018. Reorganized density profile objets (in progress).
 * Fixed test test_nc_ccl_dist.c, decreased number of tests in test_ncm_fftlog.c
     and test_ncm_mset_catalog.c. Testing ncm_csq1d.c. Updated
     binder/Dockerfile.
 * Adding missing notebooks to _DATA.
 * Updated binder notebook.
 * Removed debug messages in CSQ1D.
 * Minor fixes and new notebook.
 * Missing notebooks in Makefile.
 * Updated image used by binder.
 * Tweaking notebooks.
 * Included two notebooks.
 * Testing different methods to deal with zero-crossing mass.
 * New tutorial notebooks.
 * Modifying density profile objects. E.g., including more mass definitions.
 * Using the full SHA hash.
 * Testing a mybinder using a Dockerfile.
 * Testing methods to integrate regular singular points.
 * Testing docker jupyter notebooks
 * Testing docker jupyter notebooks
 * Testing docker jupyter notebooks
 * Testing docker jupyter notebooks
 * Testing docker jupyter notebooks
 * Testing docker jupyter notebooks
 * Testing docker jupyter notebooks.
 * Added support for matplotlib and scipy in the docker image.
 * Updated to python3 on Dockerfile.
 * Working on Dockerfile.
 * Working on Dockerfile.
 * Added necessary dist.prepare to examples.
 * Working on Dockerfile.
 * Working on dockerfile.
 * Updating Dockerfile.
 * Fixed lgamma_r declaration presence test. Fixed glong/gint64 mismatch.
 * Added detection for lgamma_r declaration and workaround when it is not declared
     but present (usually implemented by the compiler).
 * Fixed many documentation bugs (in most part by adding __GTK_DOC_IGNORE__ to the
     inline sections).
 * Added the GSL 2.2 guard back to where it was really necessary.
 * Fixing documentation bugs.
 * Removed old GSL guards from tests.
 * Added prepare_ functions on NcWLSurfaceMassDensity object. Added necessary
     prepare calls to tests.
 * Support for gcc 9 in macos.
 * Added missing header.
 * Fixed inlined sincos to use default c keywords.
 * Changed inline macro to the actual keyword in config_extra.h.
 * Several updates and fixes.
 * Missing file in branch.
 * Added rcm to ignore list in docs.
 * Added missing CPPFLAGS for HDF5.
 * Testing HDF5 in ubuntu.
 * Fixing missing doc. Testing HDF5 in ubuntu.
 * Fixed bug in preparing a fparams_map with no free variables.
 * Removed fitting test in APS.
 * Testing APS.
 * Fix pitch.
 * Improving edge cases and validation.
 * New validation on mset.
 * Removed old debug/test print in class.
 * Updated and reorganized SNIa catalogs support. Added Pantheon.
 * Working on new Boltzmann code.
 * Working on new Boltzmann code.
 * Set up CI with Azure Pipelines
 * Fix transport vectors.
 * Fixed instrospection tags.
 * Fix email.
 * Updated to new sundials interface.
 * Updated Planck likelihood code (not tested).
 * Minor changes in HIPrimTutorial.ipynb. Updated ccl interface.
 * Minor modifications related to some tests comparing with "Cluster toolkit".
 * New function to update Cls.
 * Minor improvements.
 * Testing
 * First version of the ncm_powspec_sphere_proj and ncm_fftlog_sbessel_jljm.
 * Added missing test file.
 * New FFTLog object to compute the integral with the kernel j_lj_m.
 * Fixed indentation.
 * New object NcmPowspecCorr3d (for the moment computes the simples 2point
     function). Moved filter functions from NcPowspecML to NcmPowspec. Improved
     bounds sync in Halofit.
 * Update README.md
 * Converting python examples to jupyter notebooks.
 * New Nonlinear Pk tests. Finished nc_powspec_mnl_halofit encapsulation (private
     members).
 * Fixes xcor to work with sundials 4.0.1.
 * Updating to new sundials API.
 * Updating to the new sundials API.
 * Fixed conditional use of OPENMP (mainly for clang).
 * Final tweaks for the v0.14.2 release.
 * Working on ccl vs numcosmo unit tests. Minor improvements.
 * Minor version update 0.14.1 => 0.14.2.
 * Added missing header (in some contexts).
 * Fixed aliasing problem in ncm_matrix_triang_to_sym. Removed log info from
     travis-ci.
 * Fixing doc issues, added missing docs. Better debug message for
     ncm_matrix_sym_posdef_log. Fixing travis-ci.
 * Added debug to travis-ci. Included sundials at the ignore list for
     documentation.
 * Removed no python option in numpy.
 * Lapack now is required, added openblas and lapack to travis-ci. Added
     no-undefined (when available) to libnumcosmo.
 * Added redshift direction tolerance for NcmPowspecFilter. New unit test
     CCLxNumCosmo test_nc_ccl_massfunc.
 * Finished upgrade to CLASS 2.7.1. Added option to hide symbols of dependencies.
     Finished CCL tests for background, distances and Pk (BBKS, EH, CLASS).
     Minor tweaks.
 * Updated CLASS to version 2.7.1.
 * Fixed a typo in the transverse distance. Test distances: Nc and CCL.
 * Missing unit test file.
 * Added warning for initial point in minimization not being finite.
 * Added support for CCL, first unit test for NumCosmo and CCL comparison. Minor
     improvements in NcHIQG1D.
 * Missing file.
 * Moved to xenial in travis ci.
 * Removed backports repository in travis-ci (it no longer exists...), waiting for
     something to break.
 * Removed sundials as dependence in travis-ci.
 * Updated to new sundials version -- 4.1.0.
 * Updated directory structure to match that of the new version (4.1.0).
 * Encapsulating sundials version 4.0.2. Many additions and improvements.
 * Added calibration objects for the reduced shear.
 * Working on the new ODE interface.
 * Fixed version mismatch.
 * Added back support for sundials 2.5.0 (as used by Ubuntu trusty).
 * Better support for Sundials versions, now it detects the version automatically.
     Minimum Sundials version is 2.6.0, minor updates to codes using Sundials.
 * New sampling options for NcmMSetTransKernCat object, testing new sampling
     options in test_ncm_fit_esmcmc.
 * Removed non-unsed typedef.
 * Removing unecessary headers.
 * New test to determine the burnin phase (better fitted for low self-correlation
     samplers).
 * Travis OK, removing log.
 * Travis...
 * Still testing travis.
 * Testing travis builds.
 * Fixing travis syntax.
 * Triggering travis.
 * Cleaning and updating NcmHOAA, fixing travis ci bug in mac os.
 * Fixed mac os image.
 * Testing travis ci, mac os config.
 * Debugging travis glitch in mac os.
 * Updated travis.yml to match new environment (again...)
 * Fixed minor documentation bugs.
 * Added new header for fortran lapack functions prototypes.
 * Added conditional macro compilation of the new suave support in Xcor.
 * New walker `Approximate Posterior Sampling' APS based on RBF interpolation
     using (NcmStatsDistNdKDEGauss).
 * Trying a different approach for the new walker.
 * Tests for the multidimensional kernel interpolation/density estimation object.
 * Added two codes for quadratic programming gsl_qp (gsl extras) and LowRankQP
     (borrowed from R). New ESMCMC walker Newton (does not work as expected,
     transforming in another sampler in the next commit). New multidimensional
     kernel interpolation/density estimation for arbitrary distribution
     (abstract interface and gaussian kernel implementation).
 * Removed an exit() in example_diff.py. Added support for new version for
     sundials.
 * Support for the new Sundials version.
 * Updated tests for the new interface for lnnorm computation (including error).
 * New tool for trimming catalogs, fixed typos in parameters names in NcHIPrim*.
     Fixed error estimation in the posterior normalization.
 * Fixed minor issues from codacy.
 * Fixing macos+travis issue.
 * Removed old includes in tests. Added debug in travis.
 * Moved MVND objects to main library code. Improved estimates of the Bayes
     factor, included new unit tests.
 * Change tau range. My tau definition is different from root (CERN) code. Voigt
     profile.
 * Fixed: updated deprecated glib functions. Included new Sundials version in
     configure.
 * Fixed and tested (against CCL) the galaxy weak lensing module inside XCOR.
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * Testing PITCHME.md
 * New PITCHME.md presentation. Tweaks on pydata_simple example. Minor bug fixes.
 * Included missing property "dists" (Photo-z distributions) in class_init.
 * Included other variable types to be read from hdf5 catalog.
 * Added a fix to avoid setting null sparam-array.
 * Proper serialization of model with modified parameter properties and no
     reparametrization.
 * Fixed stupid bug.
 * Turned backport function static to avoid double definition.
 * Add option to choose parameter by name in mcat_analize.
 * New options on mcat_analize to control dump. Updated parameters range in
     nc_hicosmo_de and nc_hicosmo_lcdm.
 * Fixed missing doc.
 * Finished ESMCMC MPI support.
 * New bindable ncm_mset_func1 abstract class. Included function array in
     example_esmcmc.py. Fixed bug in function array + ncm_fit_esmcmc + MPI.
 * Fixed flags search when using ifort/icc. Added missing objects registry.
 * Added guard to avoid parsing blas header during g-ir-scanner.
 * Moved HDF5 LDFLAGS to LIBS.
 * Re-included gslblas inclusion/removal in configure.ac.
 * Included additional sundials library to LIBS.
 * Inverted header order to avoid clash with CLASS headers. Added backport from
     glib 2.54.
 * Better serialization of constraint tolerance arrays.
 * Safeguards in NcmVector constructors.
 * Fixing minor leaks.
 * Fixed initializer (to remove harmless warning).
 * Added conditional compiling of _nc_hiqg_1d_bohm_f.
 * Fixed MPIJob crash when MPI is not supported.
 * Fixed leaks in ncm_fit_esmcmc.c, ncm_mpi_job_fit.c and ncm_mpi_job_mcmc.c.
 * Updated tests and removed debug info from travis-ci.
 * Fixed blas detection (ax_check_typedef is broken! Now using AC_CHECK_TYPES).
 * Debugging travis-ci build.
 * Debugging travis-ci macos build.
 * Typo in blas enum detection.
 * Detecting lapack xblas functions availability. Fixing cblas headers in
     different scenarios.
 * Fixing BLAS headers compatibility.
 * Included new functions on NcWLSurfaceMassDensity: critical surface mass
     density, shear and convergence new functions are computed when the source
     plane is at infinite redshift.
 * First working version of ncm_fit_esmcmc + MPI. Example pysimple updated.
 * First working version of nc_hiqg_1d.h (some speed-ups are still necessary).
 * Missing object files.
 * Renaming quantum gravity object.
 * Added a front-end for other lapack functions. Reorganized the blas header
     inclusion. Working in progress in ncm_qm_prop, last commit before removing
     different approaches code.
 * New examples and working in progress for NcmFitESMCMC and NcmMPIJobMCMC.
 * Finished support for complex messages in NcmMPIJob. Two implementations tested
     NcmMPIJobTest and NcmMPIJobFit. Working on NcmFitESMCMC parallelization
     using MPI.
 * Documentation.
 * Finished the first version of MPI support, including the helper objects
     NcmMPIJob*.
 * Trying different interpolation methods in ncm_qm_prop.c.
 * Workaround travis-ci problem.
 * Fixed test test_nc_wl_surface_mass_density.c.
 * Removed leftover headers in ncm_spline_rbf.c.
 * Updated example example_wl_surface_mass_density.py.
 * First working version of nc_data_reduced_shear_cluster_mass. Improvements in
     all related objects.
 * Added support (optional) to HDF5. New objects NcGalaxyRedshift,
     NcGalaxyRedshiftSpec and NcGalaxyRedshiftSpline to describe galaxy redshift
     distributions. Included support for loading hdf5 catalog in
     nc_data_reduced_shear_cluster_mass.
 * Fixed bug: ncm_util_position_angle was returning -(Pi/2 - theta). Corrected to
     return theta.
 * Removed printf in great_circle_distance function.
 * Implemented the position_angle and great_circle_distances functions.
 * Included new lapack encapsulation functions. Trying different methods in
     NcmQMProp. New 1D interpolation object NcmSplineRBF.
 * Removed unnecessary range check in gobject parameter properties. Working on
     ncm_qm_prop.
 * Removed deploy and less verbosity on make.
 * Removed extra header.
 * New version 0.14.1.
 * Included link for parallel linear solvers in sundials. Fixed bug in nc_cbe
     (lensed CMB requirements).
 * Updated examples to python3. New object nc_galaxy_selfunc. Working on
     ncm_qm_prop.
 * Fixed possible (impossible in practice) overflow in background.c. Added support
     for binder.
 * Fixed several documentation glitches.
 * Finished support for Planck lensing likelihood.
 * Created data object NcDataReducedShearClusterMass. Work on progress.
 * Fixed bug in ncm_spline.h. Working on ncm_qm_prop.
 * Included properties. Work in progress.
 * Improved regex that greps SUNDIALS_VERSION in configure.ac.
 * Working on example_qm.c, removed old file from numcosmo-docs.sgml.in, added
     quotes to SUNDIALS_VERSION grep in configure.ac.
 * Implementing object to estimate mass from reduced shear. Work in progress.
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
 * Fixed merge leftovers.
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
 * Added function eval_full to nc_xcor_limber_kernel.
 * New example NcmDiff.
 * Several improvements in the new Boltzmann code (pre-alpha). Minor fixes.
 * Reorganized the hipert usage of bg_var. Moved gauge enumerator to Grav
     namespace.
 * Added support for system reordering to lower the bandwidth. New component PB
     photon-baryon and gravitation Einstein.
 * New abstract class describing first order problems. Fixed enum type names.
 * New abstract class to describe arbitrary perturbation components.
 * Applied fix (gtkdoc/glib-mkenums scan) to the ncm namespace.
 * Fixed gtkdoc/glib-mkenums scan problem.
 * Added a guard to avoid introspection into wrong headers. Solving the enum
     parsing erros/warnings (glib's bug).
 * Update deprecated glib function.
 * First commit of new perturbation module.
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
 * Testing travis ci.
 * Trying fix numpy install using brew.
 * Right order.
 * overwrite option added to numpy at travis ci.
 * Added numpy to brew install in travis ci.
 * Added condition on having fftw to the unit test.
 * Minor tweaks. Improved NcmFftlog and added unit testing. Updated NcmABC and
     NcABCClusterNCount.
 * Update py_sline_gauss.py
 * Update py_sline_gauss.py
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
 * Added missing conditional compilation of MPI support.
 * First tests of MPI slaves.
 * Better asserts.
 * Improved interface with gsl minimizers. Included restarting for mms algorithms.
 * Testing different autoconf mpi detections.
 * Improved example.
 * Chains diag output fix.
 * Missing refs and typo.
 * Added new diagnostics to NcmMSetCatalog, max ESS and Heidelberger and Welch's
     convergence diagnostic.
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
 * Using a warning instead of a assert in the final optimization test.
 * Testing better optimization finishing clean-up.
 * Created Crocce's 2009 multiplicity function.
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
 * Updated changelog.
 * New changelog file.
 * Version bumped to 0.13.3.
 * Changed safeguard in nc_cbe
 * Missing doc tag.
 * Fixed indentation.
 * Added smoothing scale to eval by vector function.
 * Added a smooth transition from non-linear to linear power spectrum for high
     redshift in halofit. Added the znl finder to obtain the redshift where we
     should stop applying the halofit.
 * Imported from xcor branch.
 * Organized and improved (testing phase).
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
 * Added a safeguard for halofit (Brent solver, in case fdf solver crashes).
 * Updated the script mass_calibration_planck_clash.py, new funtion to include a
     gaussian prior. Corrected typos in the documentation.
 * Removed broken test for when old gsl is present.
 * Fixed double AC_CONFIG_MACRO_DIR.
 * Removed local link file.
 * Corrected a bug in halofit.
 * Corrected leaks in nc_data_xcor.c and modified ncm_data_gauss_cov.c in case of
     singular matrix.
 * Corrected some leaks in NcDataXcor.
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
 * Increasing maxsteps in xcor.
 * Modified xcor and halofit.
 * Fall back to default files.
 * Conditional usage of gsl >= 2.2 functions.
 * Conditional use of gsl_sf_legendre_array_ functions.
 * Fixed docs typos.
 * Fixed conflicts.
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
 * Included sigma8 in NcmMSetFuncList. Small adjustements. Working in progress in
     WKB.
 * Updating WKB module (work in progress).
 * plop
 * Increased output sampling of NcPowspecCBE to avoid interpolation errors.
 * Temporary debug prints.
 * Correction in nc_xcor.c and ncm_vector.h
 * Organizing code.
 * Organizing code.
 * Adding Xcor data objects.
 * Organizing and tweaking new Xcor objects.
 * Imported updated XCor codes. First tweaks and documentations fixes.
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
 * Create README.md listing and describing the scripts. Included script
     mass_calibration_planck_clash.py (ref. arXiv:1608.05356).
 * Updating TwoFluids perturbation object, working in progress.
 * References included - documentation in progress.
 * Finalizing NcmSphereMapPix (including spherical harmonics decomp).
 * Peakfinder functions were rewritten in terms og GSL functions, therefore the
     objects NcClusterMassPlCL and NcCluster PseudoCounts no longer depend on
     the Levmar library.
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
 * Fixed bug in mcat_analyze.c.
 * Imposed the same out-of-interval prior in both serial and parallel modes.
 * Moved back the default value of the parameter A of NcmFitESMCMC to 2.
 * Several fixes and improvements.
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
 * Added submodel concept in NcmModel. NcHIReion and NcHIPrim are now submodels of
     NcHICosmo.
 * Renamed submodel for stackpos (stack position) in NcmMSet internals.
 * Finished resampling for NcDataPseudoCounts and its tests.
 * Added set_cad function in DataClusterPseudoCounts.
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
 * Implemented function to compute 1-3 sigma error bars for the best fit.
 * Added message to be print when mode_error (mcat_analyze) is called.
 * mcat_analyze: implemented options mode_errors and median_errors. They provide
     the mode (median) and the 1-3 sigma error bars of a parameter.
 * Implemented functions to perform Planck-CLASH analyses considering flat priors
     for the selection and mass functions.
 * Fixed typo.
 * New NcHIReion* objects. Moved Yp to the cosmological model NcHICosmo.
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
 * new NcHIprimAtan object (primordial spectrum power law x atan). New mset_gen
     tool to generate .mset files. Added flag controling the tensor mode usage
     in NcHIPertBotlzmannCBE. New references (to the atan models).
 * Better error mensage when trying to de-serialize an invalid string.
 * Fixed parameter name in NcHICosmoDE z_re -> tau_re. Fixed NcHICosmoBoltzmannCBE
     to account correctly the lmax when using lensed Cls. Added free/fixed
     parameter manipulation functions to NcmMSet.
 * Missing HIPrim implementation PowerLaw (nc_hiprim_power_law).
 * New objects and support for primordial cosmology NumCosmo <=> CLASS.
 * First working version of the Planck+CLASS interface. Minor bug fixes.
 * Working on the CLASS interface. All precision parameters mapped.
 * Initial phase of the Class backend interface.
 * Added Class as backend. Documentation fixes. Renamed object NcPlanckFI_TT to
     NcPlanckFICorTT.
 * Missing files in the last commit.
 * New NcPlanckFI objects.
 * Updated test_nc_cluster_pseudo_counts.
 * Bug fix and initial object development.
 * Fixed last steps for making releases.
 * Added Planck likelihood 2.0 to the building system.
 * Updated to internal libcuba 4.2.
 * Resample function of nc_cluster_data_pseudo_counts is a work in progress.
     NcClusterPseudoCounts object has a new property: ncluster - number of
     clusters.
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
 * Documentation improvements (in progress): some NcClusterMass' children,
     ncm_abc.c and ncm_lh_ratio1d.c.
 * Added support for new version of SUNDIALS.
 * Fixed bootstrap support in NcmDataDist1d
 * Moved BAO data from hardcoded to serialized objects. New BAO data. Minor test
     updates.
 * Fixed bugs in ncm_lh_ratio1d from last update in this object.
 * Added more testing in test_ncm_mset, fixed minor bugs.
 * New support for multples models of the same type in NcmMSet. Minor fixes and
     updates.
 * Documentation: improvements on nc_cluster_redshift, nc_cluster_mass and
     nc_hicosmo.
 * Minor update.
 * Added check for missing set/get functions in NcmModel.
 * Pseudo cluster number counts: observable and data objects were created.
     Integral new function: function to compute tri-dimensional integral
     implemented using cuhre function (libcuba). Documentation: improvement in
     different files.
 * Added support for jerk in DE models.
 * Implementing Planck-CLASH mass function: in progress.
 * Bumped to v0.12.2
 * Fixed bug in ncm_fit_esmcmc_run_lre.
 * Added new lnsigma_lens parameter to NcSNIADistCov.
 * Tools reorganization and several improvements.
 * Added missing docs directives.
 * Improved interface to NcmLHRatio2d in darkenergy. Minor improvements.
 * New cluster mass relation: Planck-CLASH correlated mass-observable relations.
     Planck - SZ signal. CLASH - lensing signal. This object is not finalized.
     Work in progress.
 * Missing file.
 * Optimization flags.
 * Added and enable flag to include compiler's optimization/warnings flags. Made
     several minor code quality improvements.
 * Added a internal version of Cuba. Fixed minor typos and updated autogen to use
     autoreconf.
 * Better workaround for the missing fffree/fits_free_memory functions and
     SUNDIALS_USES_LONG_INT macro. Corrected version for g_test_subprocess
     usage.
 * Fixed threads competitions with OpenBLAS or MKL. Finished the NcDataSNIACov
     interface.
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
