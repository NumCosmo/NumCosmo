CHANGELOG
----------------------

[Current]

[v0.28.0]
 * Mark every 0.27 compatibility piece with "dropped in 1.0"

     * Use the same phrase in the Yp sinks, late add_submodel and the Python
     shims
     * Mark the Yp removed-parameter declaration and the legacy ~/.numcosmo code

 * Place stored parameter descriptions by name, so files from 0.27 load correctly

     * Place each NcmModel:sparam-array entry by its parameter name instead of
     its stored index
     * Add ncm_model_class_add_removed_param and skip stored descriptions of
     removed parameters
     * Abort on a stored description of any other unknown parameter
     * Declare Yp removed in NcHICosmo
     * Add 0.27-written fixtures with modified descriptions after the Yp index,
     and their generator
     * Test that each stored description lands on its named parameter and that
     unknown names abort
     * Regenerate ncm.pyi

     Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>

 * Tighten the wording of the data directory and download lock docs

     * State the abort cases in the ncm_cfg_init gtk-doc and fix a missed rename
     * Reduce ncm_cfg_fullpath_base_is_legacy to its Returns line
     * Describe NCM_CFG_HOME_ENV as setting the NumCosmo user data directory
     * Name the owner file and describe the takeover rule plainly in the lock
     docs and comments
     * Rewrite the data-cache action comment and drop an em-dash in
     CONTRIBUTING.md
     * Name the removed variables and state facts in the download test
     docstrings

 * Call the per-user directory the NumCosmo user data directory

     * Rename it in the ncm_cfg gtk-doc, error message and comments
     * Rename it in the legacy-directory notice and the WL catalog download doc
     * Rename it in TESTING.md, CONTRIBUTING.md and the test docstrings
     * Keep "data directory" for the shipped data files under NUMCOSMO_DATA_DIR

 * Add NCM_CFG_HOME_ENV for the NUMCOSMO_HOME variable

     * Define NCM_CFG_HOME_ENV next to NCM_CFG_DATA_DIR_ENV
     * Read it in _ncm_cfg_data_dir and build the error message from it
     * Regenerate ncm.pyi

 * Look up Planck baseline files only in the NumCosmo data directory

     * Drop the hardcoded ~/.numcosmo fallback from find_baseline_file

 * Pin the CI data directory to the XDG default and cache it there

     * Export NUMCOSMO_HOME=$HOME/.local/share/numcosmo from the data-cache
     action
     * Cache the downloaded files under that directory instead of ~/.numcosmo
     * Bump DATA_CACHE_VERSION in build_check.yml and weekly_flaky.yml
     * Describe the pinned directory in CONTRIBUTING.md

 * Point downloads into the legacy ~/.numcosmo at the XDG location

     * Add ncm_cfg_fullpath_base_is_legacy, TRUE when init kept an existing
     ~/.numcosmo
     * Print a one-time deprecation notice when downloading into the legacy
     directory
     * Check the legacy flag in every ncm_cfg path test case
     * Add Python tests for the notice and for the quiet XDG default
     * Regenerate ncm.pyi

 * Keep an existing ~/.numcosmo, otherwise follow the XDG spec

     * Resolve the data directory as NUMCOSMO_HOME, existing ~/.numcosmo,
     $XDG_DATA_HOME/numcosmo, ~/.local/share/numcosmo
     * Treat empty values as unset and ignore a relative XDG_DATA_HOME
     * Abort on a relative NUMCOSMO_HOME and on a data directory that cannot be
     created
     * Document the order in the ncm_cfg_init gtk-doc
     * Make the path tests table-driven and cover the legacy, relative and empty
     cases
     * Look for the wisdom files in the XDG default in test_cfg.py

 * Take over a download lock whose owner is dead on this host

     * Record pid@host in the lock directory's owner file
     * Take the lock over at once when that pid no longer exists on this host
     * Keep the timed takeover for locks with no owner or an owner on another
     host
     * Add test_a_dead_owners_lock_is_taken_over_at_once

 * Release the download lock before a failed fetch aborts

     * Pass the held lockdir to _nc_data_download_file and release it before
     each g_error
     * Assert in test_a_failed_download_leaves_nothing_behind that no lock is
     left

 * Keep the Python test children out of the real data directory

     * Add isolated_env in test_data_download.py, dropping NUMCOSMO_HOME and
     XDG_DATA_HOME
     * Drop NUMCOSMO_HOME and XDG_DATA_HOME from the wisdom child env in
     test_cfg.py

 * Make the ncm_cfg spawn helper robust and its GLib-DEBUG checks real

     * Guard g_spawn_check_wait_status for GLib < 2.70 with
     g_spawn_check_exit_status
     * Return from _test_ncm_cfg_spawn when the spawn fails instead of reading
     NULL output
     * Enable G_MESSAGES_DEBUG=GLib in the child so the setenv/unsetenv checks
     can fire
     * Match GLib's message text, which GLib prefixes with a timestamp
     * Fall back to g_spawn_check_exit_status in test_ncm_mpi_slave.c

 * Make the ncm_cfg child tests able to fail

     * Run the fftw_timelimit_default child through ncm_cfg init again
     * Return g_test_run () from main so a failing child exits non-zero

 * Regenerate the Python stubs for the 0.27 compatibility shims

     * Update HIReionCamb.set_z_from_tau in nc.pyi to the shim signature
     * Update Model.add_submodel in ncm.pyi to the shim signature

 * Keep the 0.27 API that Firecrown uses working until 1.0

     * Allow a deprecated post-construction ncm_model_add_submodel into an
     empty, declared slot
     * Add _ncm_model_class_has_slot_for so the late path still rejects
     unslotted submodel types
     * State in the ncm_model_add_submodel gtk-doc that submodels are
     construction-fixed and the late path is dropped in 1.0
     * Add numcosmo_py/_compat.py with DeprecationWarning shims, imported from
     numcosmo_py/__init__.py
     * Shim late Ncm.Model.add_submodel
     * Shim HICosmo.param_set_by_name("Yp"): set it on an NcBBNParametrized,
     otherwise ignore it as 0.27 did
     * Shim HIReionCamb.set_z_from_tau(cosmo, tau) for unattached, attached and
     tau-reparametrized reion
     * Add /nc/hicosmo_de/late_attach for the late prim and reion attach
     * Expect "declares no submodel slot" in /ncm/mset/submodel/unslotted_attach
     * Add tests/python/numcosmo_py/test_compat.py for the three shims
     * Note the deprecated late attach in the halo mass function example

 * Simplify the CI workflow and its comments

     * Drop the conda activate calls; setup-miniconda activates the environment
     in every login shell
     * Remove the commented-out gcovr and Coveralls steps
     * Rewrite the workflow comments to state the current constraint only; the
     timing measurements move to the developer notes

 * fix(test): isolate ncm_cfg child environments to avoid GLib warnings

 * add(cfg): support configurable data directory

     Use `NUMCOSMO_HOME` as the data directory override, then
     `XDG_DATA_HOME/numcosmo`, with `~/.numcosmo` as the default. Update related
     documentation and test all three path choices.

 * Check the Python stubs in CI and move the mypy config to pyproject.toml

     * Regenerate nc.pyi and ncm.pyi with update_pyi.sh in lint-built and fail
     if they differ from the committed stubs
     * Run the stub check and mypy even when the other fails
     * Move the mypy configuration from .mypy.ini to [tool.mypy] in
     pyproject.toml
     * Describe the stub check in CONTRIBUTING.md

 * Restore the blank lines lost with the pylint comments

     * The removal pattern also matched the newlines before each comment line;
     put back the blank lines it deleted
     * Join the lines black now fits on one line in cluster_wl.py

 * Replace flake8 and pylint with ruff

     * Check numcosmo_py with ruff's default rule set in the lint job instead of
     flake8
     * Configure ruff in pyproject.toml (line length 88, external/ and stubs
     excluded) and delete .flake8 and .pylintrc
     * Swap flake8 and pylint for ruff in environment.yml and the dev extras,
     add black to the dev extras, and regenerate the conda locks
     * Ignore mypy errors in numcosmo_py.external, which mypy 2.4 reports
     through imports despite the exclude
     * Apply ruff's fixes: import order, builtin generics, X | None and X | Y,
     redundant casts, nested ifs, enumerate, itertools.pairwise and the sorted
     __all__ in cluster_richness
     * Keep keys() on NcmVarDict and NcmObjDictStr, which are not iterable
     * Use int(np.rint(...)) in _nside_from_npix so the result is an int under
     numpy 1.x too
     * Move the frozen default generators and integration options in cluster_wl
     to module-level constants
     * Remove the dead try around Ncm.cfg_get_fullpath_base in planck_lite
     * Return Self from the generator enum __new__ methods and mark
     NcmHighlighter.highlights as ClassVar
     * Keep the log file open in logging and the data-error exceptions in
     generate and loading, with noqa and a reason
     * Drop the shebangs from two experiment modules
     * Remove the pylint comments outside vendored code
     * Update CONTRIBUTING.md for black, ruff and mypy

 * Type-check untyped function bodies with mypy and fix what it finds

     * Enable check_untyped_defs in .mypy.ini and move the mypy excludes there
     (meson, external/), so CI runs plain mypy -p numcosmo_py
     * Check generate_stubs.py with mypy again; it passes
     * Replace the removed priors_add_gauss_param in
     mass_calibration_planck_clash with priors_add(PriorGaussParam.new(...))
     * Create an unseeded RNG in MockGenerator when seed is None instead of
     passing None to seeded_new
     * Raise ValueError in MockGenerator when the HMF paths run without hmf or
     cluster_m
     * Type MockGenerator halo_set_size and cluster_set_size as int
     * Type ClusterWL dist as Nc.Distance
     * Drop None from the cluster_redshift_type, cluster_mass_type and
     fitting_sky_cut options in generate
     * Assert theta is set before the APES likelihood evaluates
     * Annotate get_galaxy_coords with array types, the rejection sample list
     and the SkyMatch coordinate dicts in mock.py; build halo_z as an array
     * Pass str(e) to typer.BadParameter and chain the ValueError
     * Ignore the mypy call-arg false positive on the three generator enum
     lookups
     * Replace Optional[X] with X | None outside the stubs and drop the unused
     Optional imports

 * Move the CI quick checks into two lint jobs

     * Add a lint job without a build: conda lock files, doc style, clean
     notebooks, uncrustify, black and flake8, each step running even when an
     earlier one fails
     * Add a lint-built job that builds and installs NumCosmo, then runs mypy
     * Remove black, flake8 and mypy from the build-miniforge legs and
     uncrustify from build-gcc-ubuntu
     * Remove the notebook check from build-gcc-ubuntu and build-gcc-macos
     * Fold the check-conda-locks and check-doc-style jobs into lint
     * Drop the commented-out uncrustify and pylint steps
     * Add uncrustify to environment.yml unpinned, so the weekly lock update
     moves it with black and flake8
     * Remove uncrustify from the apt dependencies
     * Regenerate the conda locks: uncrustify 0.83.0 and flake8 7.4.1
     * Reformat nc_galaxy_shape_intrinsic_mode.h with uncrustify 0.83.0
     * Rename the check-conda-locks job references in CONTRIBUTING.md and
     setup-miniforge

 * Improve NcmDiff robustness, convergence and derivative methods

     * Refactor Richardson convergence around independent ladder controls shared
     by vector, Hessian and dual-series drivers.
     * Make derivative convergence robust to zero derivatives,
     cancellation-dominated steps, non-finite evaluations and poorly scaled
     initial steps.
     * Add domain-aware differentiation with automatic one-sided fallbacks and
     forward second differences near boundaries.
     * Make dual-series Richardson the default method; keep the single-ladder
     scheme available as the faster heuristic alternative.
     * Require cross-ladder agreement and conservative error bounds in
     dual-series mode to reject aliased plateaus.
     * Add adaptive Chebyshev spectral first- and second-derivative methods with
     truncation and cancellation error estimates.
     * Add relative and absolute function-precision controls, reject Richardson
     rows dominated by a single term, and correct second-derivative cancellation
     estimates.
     * Add a derivative benchmark battery and extend C and Python coverage
     across both Richardson schemes, spectral derivatives, domains and difficult
     cancellation cases.
     * Add the NcmDiff theory documentation and align internal names, comments
     and properties with its terminology.
     * Simplify NcmFit numerical Hessian diagnostics by computing the Hessian
     once and reporting the normalized worst error.
     * Regenerate the Python stubs.
 * chore: ignore local conda environment

     Make git ignore the .conda-env/ subdirectory, which a worktrunk workflow
     can use to create a conda environment appropriate for a git worktree.

 * Restore lcov configuration.

 * Make the catalog evidence estimator selectable and fix its short-catalog hang
     and bootstrap

     * Add NcmMSetCatalog:post-lnnorm-method, choosing among the box, bootstrap
     and ellipsoid estimators of ncm_mset_catalog_get_post_lnnorm
     * Give the ellipsoid estimator the error from slices of the rows, through
     the box estimator's sum with a cut
     * Stop the slice loop at one slice; catalogs under 100 rows looped forever
     after dividing by zero
     * Remove the bootstrap's per-resample printf and stop it when the mean is
     known to its tolerance, instead of after about 1e4 resamples
     * Set a NaN error when the ellipsoid's covariance is not positive definite
     * Test the three estimators against the exact evidence of a Gaussian, with
     wide and cutting bounds, the ellipsoid shrunk to the box, and a 50-row
     catalog

 * Cover the catalog, MPI serial path and flat kernel, and warn on a singular
     catalog covariance

     * Exclude the not-implemented default methods from coverage, as
      elsewhere in the library
     * Run ncm_mpi_job at one rank, which takes the serial path
     * Test the catalog distribution, p-value and interval functions with
      messages, the parameter distribution, and the by-chain diagnostics on
      short and long catalogs
     * Factor the catalog covariance before the evidence's ratio estimate, so
      a singular one warns and gives zero instead of aborting; test it with
      identical rows
     * Test the flat transition kernel

 * Give invalid least-squares trial steps a finite residual

     * NcmFitGSLLS fills the residual at an invalid point with sqrt
      (GSL_DBL_MAX / (2 n)) instead of infinity; a BLAS returning a NaN norm
      for an infinite entry made GSL accept the step, and the run looped to
      the iteration limit

 * Report the GSL least-squares solver status when a run fails

     * NcmFitGSLLS logs the GSL status and its info code on failure, as the GSL
     multimin backends do
     * Run the invalid-step test with messages on stderr, so a failure reports
     the solver status outside the TAP stream

 * Return the refill input buffers in the async MPI job loop and finalize MPI
     after its own exit handlers

     * ncm_mpi_job_run_array_async returns the pooled input buffer of every
     refill send; freeing the job hung waiting for it
     * ncm_cfg_init registers its exit handler after MPI_Init, so MPI_Finalize
     runs before the handlers MPI_Init registers; MPICH workers stalled in
     MPI_Finalize and were killed
     * Test a master that runs its inputs slowly, so the workers return while
     inputs are queued

 * Locate a mode at zero in NcmStatsDist1d and draw the EPDF test sample in
     sequence

     * Floor the mode's absolute tolerance at sqrt (reltol) times the grid
     spacing, so a mode at zero stops refining
     * Give up on the mode only after ten iterations with an unchanged bracket;
     Brent may pause for one
     * Report the iterations actually run when the mode's tolerance is not
     achieved
     * Test a mode exactly at zero on three knot sets that failed before
     * Draw the weighted EPDF regression sample in sequence: C leaves the order
     of argument evaluation unspecified, and clang drew another sample

 * Fix clang and GIR warnings, the TAP noise and three-rank MPI tests on two cores

     * Capture the phase-label subprocess output without echoing it into the TAP
     stream
     * Reword the ncm_csq1d_evolve_prop_vector doc so GIR does not read a second
     @frame tag
     * Name the ncm_stats_vec_get_mean_vector and ncm_fit_esmcmc_set_sampler
     parameters as in the headers
     * Zero-initialize NcmFuncEvalCtrl with {0}
     * Spell out every field of the GOptionEntry and spline test sentinels
     * Cast the test_ncm_stats_vec loop bound to guint
     * Allow OpenMPI to oversubscribe in the MPI tests, so three ranks run on
     two-core runners

 * Handle the empty catalog where the best-fit row is used

     * get_bestfit_row is nullable (empty catalog); the catalog app and the
     plotting helper
      used it unchecked and now raise a clear error (found by mypy with the new
     stubs)

 * Fix the NcmCSQ1D propagator shift from a nonzero initial time and test the
     propagator in C

     * evolve_prop_vector changed the initial state to the propagator frame with
     int q m nu^2
      instead of int q^2 m nu^2 (the frame changes use the latter); from ti > 0
     the result
      was wrong (J error 2.9e-2 to 80), from ti = 0 both vanish
     * prepare_prop placed its first knot at ti with the values of tii
     * Document the propagator, its stop, tf_prop, the propagated frames, the
     non-adiabatic
      vacuum and compute_H; the compute_prop_vector doc carried the wrong name
     (stubs
      regenerated); drop two commented-out blocks
     * Tests: /ncm/csq1d/prop/hankel, prop/evolution, nonadiab/vacuum and the
     propagator
      refusals, against the exact non-adiabatic Bessel state
     * Remove the Python propagator tests now covered in C

 * Fix the NcmCSQ1D boost series stop and fail on boosts beyond double precision;
     test the frames in C

     * ncm_util_mln_1mIexpzA_1pIexpmzA stopped its series on the last term
     relative to the
      sum: terms vanish at rho = 0 for some theta (the origin boosted to ADIAB2
     was off by
      1.7e-4) and a zero boost or theta = pi/2 made the stop NaN and the loop
     endless; it
      now stops on the bound of the remaining terms
     * Frame changes abort on a boost of rapidity -ln(epsilon) or more
     (NONADIAB2 at t = -20
      returned garbage); document the frames, their transformations and their
     precision
     * Tests: /ncm/csq1d/frame/changes (round trips and distance invariance at
     1e-15),
      frame/adiab_vacuum, frame/eval_at, both frame aborts, and the failing
     points of the
      boost series; the Bessel model gains the three m nu^2 integrals
     * Remove the Python frame tests now covered in C

 * Fix the NcmCSQ1D adiabatic maximum at an interval end and test the adiabatic
     vacua against Hankel

     * find_adiab_max bracketed a zero-width interval for the borders when the
     minimum of
      |F1| is at an end (GSL abort with the error handler, garbage without); a
     border is
      now the end when |F1 - F1_min| does not reach epsilon on that side
     * Document compute_adiab from the code (no max-order-2 property: the order
     follows the
      initial condition type; the meaning of both error estimates, the
     truncation warning)
      and the two finders
     * Tests: /ncm/csq1d/adiab/hankel (order scaling, estimates above the error,
     the vacuum
      within the tolerance) and /ncm/csq1d/adiab/finders
     * Remove the Python initial-condition tests now covered in C

 * Read the NcmCSQ1D evolution from the time of ad hoc conditions and test the
     evolution in C

     * eval_at and eval_at_frame used the adiabatic expansion before
     vacuum_final_time, which
      ad hoc conditions never set (0): with negative times they returned the
     expansion
      instead of the evolution (distance 49 at t = -0.2 from the same state);
     they now read
      the evolution from the time of the conditions and abort before it
     * Document the regimes of eval_at and the times of get_time_array
     * Tests: /ncm/csq1d/evolution/hankel (the evolved complex structure against
     the exact
      Bessel state, from the second and fourth order vacua),
     /ncm/csq1d/evolution/ad_hoc,
      and eval_at before ad hoc conditions in /ncm/csq1d/prepare/aborts
     * Remove the Python evolution tests now covered in C

 * Give the NcmCSQ1D phase an origin at the start of the evolution, defined on
     both sides of it

     * Drop the pi/4 start of the default int nu; theta = int nu + delta theta
     from the start
      of the numerical evolution t0, where the phase stays small (from ti int nu
     reaches
      1e40 in the QGW background)
     * Before t0 both integrals run backwards on the adiabatic state of
     compute_adiab, so
      phase differences across t0 are defined (they were zero before t0)
     * The phase splines keep machine-precision relative tolerance with an
     absolute
      tolerance of reltol radians; eval_int_nu aborts before
     prepare_phase_splines
     * The saved evolution starts with the initial state; the splines began at
     the first
      step, and eval_at extrapolated between the vacuum time and it
     * Document the phase origin and the evaluated form on the theory page and
     in the API
     * Test /ncm/csq1d/phase/hankel against the exact Hankel phase; remove the
     two Python
      tests of the old convention

 * Document the NcmCSQ1D settings and prepare from the code and drop an unused
     model control

     * The NcmModelCtrl was forced by every setter and never read; remove it
     * Fix the docs of tf (called the initial time), abstol, set_init_cond_adiab
     (does not
      change ti) and NONADIABATIC2 (prepare aborts: not implemented)
     * Document the adiabatic threshold switching, the propagator threshold, the
     vacuum
      options, the initial-condition types, the evolution variables and prepare
     * Tests: /ncm/csq1d/properties and /ncm/csq1d/prepare/aborts; the state
     checks gain the
      explicit half-plane and disc formulas and circles up to r = 1e4
     * Remove the Python tests now covered in C (test_csq1d and the six state
     tests)

 * Make the NcmCSQ1DState distance and circle exact near coincident points and
     give get_phi_Pphi the phase of theta

     * compute_distance returned NaN between equal points with |alpha| > 1 and
     lost small
      distances (1e-12 returned 0); it now sums the non-negative terms of
     cosh(d) - 1 in
      logs
     * get_circle returned NaN when the circle point has alpha = 0
     * get_phi_Pphi returns the phase in which phi is real and positive, the one
     in which
      e^{-i theta} (phi, P_phi) solves the equations of motion (the paper's
     eigenvector
      drifts from it through alpha); J and the Wronskian are unchanged
     * Document the state maps from the code; the theory page gives the
     difference from
      Penna-Lima et al. (2023) and explains why only two-time quantities need
     theta; add
      the paper to references.bib
     * test_de_cont_nonadiab compares the states by their hyperbolic distance,
     which does
      not depend on the phase
     * Tests: /ncm/csq1d/state/maps, distance, distance/exact and circle

 * Fix the propagator solver check in NcmCSQ1D and document its methods from the
     code

     * init checked the wrong linear solver after creating the propagator's
     * NcmCSQ1DClass padding 2 -> 6 (12 virtual functions)
     * The aborting defaults name the method and the class; fix the abstol blurb
     (stubs
      regenerated)
     * Document the virtual functions: definitions, defaults (nu^2, m = e^xi /
     nu, F2 by
      numerical differentiation of F1) and which integrals the non-adiabatic
     parts need
     * New tests/c/ncm/dynamics/test_ncm_csq1d.c with a shared Bessel model;
      /ncm/csq1d/defaults checks the default nu^2, m and F2

 * Document NcmMPIJobTest from the code and drop its unused state

     * Drop an unused return vector, an RNG seeded at random on every
     construction and an
      unused include
     * The buffer checks compared pointers truncated to 32 bits; compare the
     pointers
     * Document the job from the code: the input is an index, the return the
     job's vector
      at it; there is no simulated workload
     * Property blurb for vector (stubs regenerated)

 * Keep unevaluated proposals out of the MPI step statistics and document
     NcmMPIJobMCMC

     * The MPI path of NcmFitESMCMC fed offboard proposals (never evaluated,
     holding an old
      -2 ln L*) and non-finite ones into the step statistics of the diagnostic
     summary;
      it now skips them as the serial path does
     * NcmMPIJobMCMC: the buffer checks compared pointers truncated to 32 bits;
     drop an
      unused timer and property id
     * Document the job from the code: input layout, acceptance rule, negative
     deviates for
      the initial points, bounds checked by NcmFitESMCMC
     * Tests: /ncm/mpi/job/mcmc/run_array and run_array_async check each
     decision against
      the rule; /ncm/fit/esmcmc/parity/stats_serial_vs_mpi compares the
     diagnostic summary
      of serial and MPI runs with offboard proposals

 * Document NcmMPIJobFEval from the code and test its returns across ranks

     * The buffer checks compared pointers truncated to 32 bits; compare the
     pointers
     * Drop an unused property id
     * Document the job from the code: -2 ln L with priors and the functions at
     each input,
      no library user
     * Tests: /ncm/mpi/job/feval/run_array and run_array_async compare each
     return with the
      same evaluation on the master

 * Document NcmMPIJobFit from the code and test its returns across ranks

     * The buffer checks compared pointers truncated to 32 bits; compare the
     pointers
     * Drop an unused property id
     * Document the job from the code: return layout, function requirements, no
     library
      user, convergence not reported
     * Tests: /ncm/mpi/job/fit/run_array and run_array_async compare each return
     with the
      same fit run on the master

 * Move the MPI worker loop to ncm_mpi_slave.c and test it with raw MPI

     * ncm_mpi_slave_serve_job serves one job from INIT to FREE or KILL, one
     function per
      command; ncm_mpi_slave_run serves jobs until KILL
     * ncm_cfg_init runs it on the worker ranks and exits; MPI start-up and the
     public
      ncm_cfg_mpi functions stay in ncm_cfg.c
     * The worker no longer runs a GLib main loop around the blocking receive,
     which also
      drops a 100 ms wait between jobs
     * ncm_mpi_slave.h is internal (not installed, not introspected)
     * test_ncm_mpi_slave.c skips ncm_cfg_init: rank 1 serves jobs, rank 0
     speaks the
      protocol with raw MPI calls (work, reinit, pending returns above the eager
     limit,
      kill)
     * ncm_mpi_slave_abort launches each refusal (INIT twice, WORK before INIT,
     unknown
      command) in its own mpiexec and checks the exit status and the message
     * The shape job moves to a shared test unit; MPI test executables link
     mpi_c_dep

 * Test the MPI job protocol across ranks, shapes and datatypes

     * A test-only job with different input and return lengths through the
     default pooled
      buffers of NcmMPIJob, which no library job uses
     * float_input: inputs as MPI_FLOAT, returns as MPI_DOUBLE; aborts with
     MPI_ERR_TRUNCATE
      when returns are received with the input datatype (the fix in the previous
     commit)
     * Empty input arrays, a single input (fewer inputs than workers), and
     workers serving
      several jobs of different types and shapes in sequence
     * MPI C tests take an optional list of rank counts; ncm_mpi_job runs at two
     and three

 * Receive MPI job returns with the return datatype and run job arrays without MPI
     support

     * run_array and run_array_async received return messages with the input
     datatype
      (harmless while every job uses MPI_DOUBLE for both)
     * Without MPI support run_array and run_array_async run every input on the
     master, as
      with a single rank, instead of aborting; init_all_slaves does nothing
     * NcmMPIJobFit declared a return buffer size of one double without
     functions; unused,
      since the job hands out its return vector's own storage
     * The default vfuncs abort naming the method and the class; drop an empty
     static and
      an unused array in the worker loop
     * Document NcmMPIJob from the code (the master runs the last n/n_ranks
     inputs of
      run_array)
     * Test: tests/c/ncm/mpi/test_ncm_mpi_job.c (run_array, run_array_async), in
     the MPI
      lane at two ranks

 * Keep short catalogs whole in the burn-in criteria and document the catalog
     trimming from the code

     * calc_max_ess_time returned the catalog length below ten iterations, so
     trim_by_type
      (ESMCMC auto trim) dropped every row; it now returns zero with max_ess
     zero
     * calc_heidel_diag returns zero below ten iterations whatever the message
     level (it
      aborted without messages and returned the length with them) and frees its
     p-values
     * heidel_diag_by_chain returned -1 as a guint when no chain passed and
     never set
      wp_pvalue for a single chain; it returns zero, and a single chain goes
     through the
      same loop
     * The by-chain diagnostics return zero below ten iterations instead of
     aborting; their
      out arguments are optional
     * remove_last_ensemble makes a backup only when a file is attached
     * Document trim, trim_p, trim_oob, remove_last_ensemble, the burn-in
     criteria,
      trim_by_type and the by-chain diagnostics (one chain and several were
     swapped)
     * Tests: trim_by_type/short, heidel, heidel/by_chain/fail and
     remove_last_ensemble for
      both fixtures

 * Fix the catalog distribution and interval functions on short catalogs and
     document them from the code

     * calc_ci_interp, calc_pvalue, calc_distrib and the parameter distributions
     divided by
      zero on catalogs with fewer than 100 rows (the progress step ran with
     messages off)
     * The confidence interval and p-value functions accept a single p-value or
     limit
     * calc_pvalue returns one column per limit instead of 2n + 1, n + 1 of them
     uninitialized
     * param_pdf includes the largest value in the last bin; param_pdf_pvalue no
     longer
      aborts at the maximum, no longer counts the bin below the value, and gives
     one below
      the sampled range
     * ensemble_evol asserts at least two grid points
     * Document the shrink factors (the refined MPSRF is not computed; 1e10 on
     failure),
      param_pdf, the confidence interval, p-value, distribution and ensemble
     evolution
      functions
     * Tests: distrib/short, ci/single and param_pdf for both fixtures

 * Return the error of a cached catalog log evidence and document the evidence
     estimators

     * ncm_mset_catalog_get_post_lnnorm left the error unset when it returned
     the cached
      value (every call after the first, including the one in get_post_lnvol);
     the
      hyperbox estimator leaked its generator on a failed Cholesky
     decomposition;
      get_post_lnvol left glnvol unset on its early return
     * Document get_post_lnnorm (ln Z for a flat unit-density prior on the box,
     the
      truncated-Gaussian estimator, the caching), get_post_lnvol, the m2lnL
     quantile, best
      fit, mean and covariance getters; get_bestfit_row is nullable (stubs
     regenerated)
     * Tests: norma checks the error of the cached value

 * Document the NcmMSetCatalog getters and row input from the code

     * The weight of a weighted catalog is the last additional value; the
     zero-mean
      fallback of largest_error applies in [1, 2); get_row_from_time takes a row
     id;
      get_cur_id is the first id minus one when empty; col_by_name also takes
     additional
      value names and column numbers
     * The add functions give the row layout they expect; burn-in, getters and
     the
      read-only blurb condensed; American spelling (stubs regenerated)

 * Fix the catalog reset of the per-ensemble arrays and a generator leak, explain
     K_eff

     * ncm_mset_catalog_reset and reset_stats kept the ensemble variances and
     acceptance
      ratios, so after a reset, trim or remove_last_ensemble peek_e_var_t and
     the
      acceptance array returned stale entries before the new ones
     * set_rng leaked the previous generator of an empty catalog; clearer
     set_file and
      set_rng errors
     * Document the catalog: row and chain layout, first id, Markov chain start,
     burn-in;
      constructors, free and clear (they sync), reset, erase, sync and the
     setters; the
      first-id links pointed to a property that does not exist
     * numcosmo catalog analyze prints what the K_eff column measures
     * Tests: reset with one and many chains, set_rng twice

 * Fix the ESMCMC diagnostic columns and document NcmFitESMCMC from the code

     * The ensemble summary printed cor(p, a) and cor(q, a) under each other's
     headers;
      name the correlations after the statistics they hold
     * set_sampler's doc block carried another function's name; the update and
     ensemble-log
      helpers are static; set_data_file requires a file name; clearer start_run
     error
     * Document NcmFitESMCMC: the ensemble halves, rejections, initial ensemble,
     catalog
      layout, run cycle and resumption, threads and MPI; run counts iterations;
     use_mpi
      does not depend on use-threads; the trimming, min-runs and max-runs-time
     setters and
      validate say what they do

 * Document NcmFitESMCMCWalkerAPES from the code

     * The type doc describes the NcmStatsDist estimators, the
     independence-sampler
      acceptance, the exploration phase and links the theory page; drop the
     outdated
      radial-basis description, nonexistent arguments and the static apes.png
     * new_full, the robust and fixed covariance setters, set_local_frac,
     set_over_smooth
      and the method and kernel enums say what they do; American spelling (stubs
      regenerated)

 * Fix the multi-stretch acceptance and the walk move, set the walkers' dimensions

     * The multi-stretch acceptance omitted the product of the stretches' z^(d -
     1): its
      posterior variances were 0.44 of the true ones in four dimensions
     * NcmFitESMCMCWalkerWalk never got its number of parameters (one from its
     constructor),
      so it moved only the first coordinate; NcmFitESMCMC now sets the walkers'
     size and
      number of parameters to its own
     * The walk step used the walker index as the coordinate index; it is now
     the
      Goodman-Weare walk move, sum_j z_j (X_j - Xbar)
     * Document the walker base (acceptance ratio, prob_norm as ln q; the long
     MCMC review
      becomes a reference), the stretch and walk moves with their proposal
     factors, the
      walk scale range; default methods abort naming the class; padding comment
     * Tests: test_ncm_fit_esmcmc_walkers (stretch, multi-stretch and walk
     against N(y, C),
      every coordinate moving)
     * NcmFitESMCMCWalkerWalk scale default 1 (was 0.2): steps of about one
     standard
      deviation, acceptance 0.49 and tau 20 against 0.88 and 120 at 0.2 in four
     dimensions

 * Fix the catalog kernel with repeated rows and bound its sampling, document the
     kernels

     * NcmMSetTransKernCat RBF sampling aborted when the catalog tail had
     consecutive
      repeated rows (the m2lnL vector kept their slots); RBF and KDE sampling
     redrew
      without limit until inside the bounds; they now draw once and
      ncm_mset_trans_kern_prior_sample draws again
     * NcmMSetTransKernFlat did not lock the RNG; NcmMSetTransKern vfunc
     defaults abort
      instead of being NULL; drop a stray parent_instance; fix two Gauss error
     messages
     * Document the kernels: the proposal of each, symmetry, the catalog kernel
     as a prior
      sampler; free, set_prior, pdf and prior_pdf; padding comment
     * Tests: the catalog kernel samples an RBF interpolation built from a
     catalog with a
      repeated row
     * Remove the unused bernoulli_scheme field of NcmMSetTransKernClass (stubs
     regenerated)

 * Fix NcmFitMCMC acceptance: reject invalid proposals and keep the Gauss kernel
     symmetric

     * A NaN -2 ln L gave acceptance probability one and every later proposal
     was accepted;
      proposals outside the bounds, at invalid parameters or with a non-finite
     -2 ln L are
      rejected, as in NcmFitESMCMC
     * NcmMSetTransKernGauss drew again outside the bounds, which made it
     asymmetric near
      them and biased the chain (half-normal mean 0.1 sigma off); generate draws
     once and
      ncm_mset_trans_kern_prior_sample draws again for initial points (1000
     draws, then
      abort); remove the Gauss max-iter property
     * set_trans_kern leaked the previous kernel and named set_rtype in its
     error; the
      data-file property was never installed; the update helper is static;
     start_run named
      set_data_file in its error
     * Document NcmFitMCMC: the acceptance rule and the symmetric-kernel
     requirement, the
      rejections, the run cycle; stubs regenerated
     * Tests: test_ncm_fit_mcmc (posterior N(y, C), half-normal at a bound,
     invalid region);
      the Gauss kernel tests check one draw per generate and prior-sample
     exhaustion

 * Fix NcmFitMCBS construction and bootstrap type, document and test it

     * ncm_fit_mcbs_new crashed: the fit property read the model set before
     storing the
      fit
     * The bootstrap runs always used the total bootstrap, so BOOTSTRAP_NOMIX
     ran as
      BOOTSTRAP_MIX; the type is now checked before the run starts; drop an
     unused file
      name
     * Document NcmFitMCBS: what the catalog holds, the bootstrap file names,
     the
      realization range; the fiducial model is nullable (stubs regenerated)
     * Tests: test_ncm_fit_mcbs (a run: one bootstrap mean per realization,
     covariance
      set)

 * Document NcmFitMC and test its estimator against the Gaussian resampling

     * Document NcmFitMC: resampling, the catalog of best fits, the run cycle,
     restarts
      from a data file; every public function and NcmFitMCResampleType
     * The catalog update and bootstrap helpers are static; start_run named
      set_data_file in its error; set_rtype checks the type before looking it
     up;
      set_data_file requires a file name; the fiducial model is nullable (stubs
      regenerated)
     * Shorten the keep-order and shared-object comments to the constraints they
     state
     * Tests: test_ncm_fit_mc (4000 realizations from the model: m2lnL of every
     best fit,
      mean and covariance against the fiducial and the data covariance within 5
     sigma,
      mean_covar)

 * Fix NcmLHRatio2d's border, drop its unused parts and test it against the
     Gaussian profile

     * conf_region appended a point one step past the last border point (4.9%
     off the
      border in the test); the region is now the border points closed by the
     first
     * The point list was never freed (points_free freed an empty list); the
     constructor
      checked the first parameter's free status for both and took unsigned fpis
     * One extra constrained fit per root iteration removed (320 instead of 480
     fits in the
      test); "Start found" follows the message level; fisher_border's point
     count no longer
      depends on rounding
     * Remove the unreachable Steffensen root with NcmLHRatio2dRoot and its
     NcmDiff, the
      dead second-try logic and the unused exported point helpers
     (NcmLHRatio2dPoint is
      private); enum nick table and stubs regenerated
     * Document NcmLHRatio2d, the border walk and the regions; NcmLHRatio1d
     names the
      functions that provide the covariance
     * Tests: test_ncm_lh_ratio2d (conf_region and fisher_border against d^T
     M^-1 d =
      chi2_2 to 1e-8, point counts, closure)
     * One extra constrained fit per root iteration removed (320 instead of 480
     fits in the
      test); the border walk aborts after 100 times the expected number of
     steps;
      "Start found" follows the message level; fisher_border's point count no
     longer
      depends on rounding

 * Fix NcmLHRatio1d after a Fisher matrix, drop its unused parts and test it

     * ncm_fit_obs_fisher reset the fit state, clearing is_best_fit and the
     best-fit m2lnL,
      so NcmLHRatio1d aborted after it; it now resets only when the numbers of
     residuals or
      free parameters changed, and the Hessian retry no longer resets
     * NcmLHRatio1d ran one extra constrained fit per root iteration (the root
     solver already
      reports non-finite values; 32 instead of 44 fits in the test) and warned
     twice when
      the interval reached a parameter bound
     * Remove the unused constraint property (it also leaked) and the
     unreachable Steffensen
      root with NcmLHRatio1dRoot; enum nick table and stubs regenerated
     * Document NcmLHRatio1d: the bounds are offsets from the best fit, the fit
     must hold a
      covariance, the search and its precision
     * Tests: test_ncm_lh_ratio1d (bounds against sqrt (chi2_1 C_ii) to 1e-8,
     set_pindex, an
      upper parameter bound)

 * Port NcmFitGSLLS to gsl_multifit_nlinear and fix the fit backends' resets and
     results

     * NcmFitGSLLS uses the GSL trust-region solver (gsl_multifit_nlinear,
     Levenberg-Marquardt)
      instead of the legacy fdfsolver; a trial step to parameters the models
     report invalid
      now gets an infinite residual and is rejected, where the legacy run ended
     with GSL_EDOM
     * GSL least squares, multimin and simplex fits built without free
     parameters crashed or
      hung once parameters were freed (callbacks never set, stale solver,
     infinite first
      step); every reset now recomputes the steps and reallocates only when the
     count changes
     * NcmFitGSLMM returned GSL_EDOM (33) as -2 ln L at invalid points, and fdf
     computed the
      real value there; both now give +inf, so the line minimization backtracks
     * NcmFitGSLMM set_algo freed the description without clearing it (use after
     free and a
      double free) and checked the old algorithm; the fit state stored the last
     line-search
      point instead of the minimum
     * NcmFitGSLMMS restarts from the best vertex instead of the last evaluated
     one
     * NcmFitNLOpt reset installed the local algorithm as the main one, the
     descriptions never
      followed an algorithm change, and the constraint data lived in a GArray
     that was
      reallocated while NLopt held pointers to it
     * GSL multimin, simplex, levmar and NLopt runs report no best fit when
     stopped by the
      iteration limit or an error, as GSL least squares did
     * NcmFitState m2lnL_prec no longer asserts a value below one (and its test
     goes); an
      early-stopped gsl_mm or levmar run aborted there
     * Document the five backends and NcmFitLevmarAlgos; algo_name of the
     new_by_name
      constructors is nullable (stubs regenerated)
     * Tests: /run/empty/free, /run/maxiter, /gsl_ls_invalid_step,
     /gsl_mm_set_algo,
      /gsl_mms_set_algo, /nlopt_set_algo and /constraints/equality/two for every
     algorithm

 * Fix NcmFit covariance and accurate derivatives, the NcmDiff Hessian round-off
     and the levmar stopping tests

     * The covariance of an indefinite Fisher matrix (LU fallback) was halved
     twice: twice the
      inverse, and its log-determinant off by 2 ln 2
     * The accurate least-squares Jacobian aborted unless there were as many
     residuals as free
      parameters; it is now transposed into place
     * The accurate gradient, Jacobian and Hessian set parameters through
      ncm_fit_params_set_vector, so a sub-fit reruns at every point as for
     forward and central
     * NcmDiff Hessian step returned a relative round-off compared with an
     absolute truncation
      error: wrong estimates, and wrong values at a zero coordinate (Rosenbrock
     H12 at (0, 1)
      was 28366 instead of 0)
     * NcmFitLevmar stops at the fit's params_reltol and m2lnL_abstol instead of
     fixed 1e-15;
      the bc algorithms never stopped with a nonzero residual at the minimum
     * Forward val_grad counted one evaluation too many; lr_test documents the
     upper-tail
      p-value and lr_test_range requires two steps; clearer messages and docs
     for the ls and
      gradient functions and get_covar
     * Tests: test_ncm_fit_covar (indefinite Fisher, accurate non-square
     Jacobian, sub-fit with
      accurate derivatives against the closed-form profile); NcmDiff Hessian at
     zero
      coordinates and for Rosenbrock

 * Fix a leak in NcmFit restarts without free parameters and document running and
     logging

     * run_restart returned before freeing its serializer when the fit has no
     free
      parameters
     * Document run, run_restart, reset, the logger and the logging functions;
      priors_m2lnL_val covers the priors alone, not the data
     * _ncm_fit_run_empty is static; the log prints its counters with %u

 * Fix NcmFit covariance lookups with a parameter that was not fit

     * covar_cov and covar_cor checked the first free-parameter index twice, so
     a
      second parameter that was not fit read the covariance at index -1; they
     now
      abort, and the free-parameter accessors assert their indices
     * get_sub_fit returns NULL without a sub-fit instead of referencing NULL
     * set_sub_fit: drop the disabled free-in-both loop, which cannot hold since
     a
      shallow-copied NcmMSet shares the models and their free/fixed types; its
     doc
      states the requirements
     * Document the setters, parameter updates (they rerun the sub-fit) and
      constraints; error messages name their function
     * Tests: /covar/not_fit and /sub_fit/none for every fit algorithm
     * Regenerate ncm.pyi

 * Fix the NcmFit degrees of freedom after a reset and the levmar workspace

     * _ncm_fit_reset, which every run calls, left the priors out of the degrees
     of
      freedom (0 instead of 2 with two priors); construction and reset share one
      computation now
     * NcmFitLevmar kept its workspace sized for the old data length when a
     reset
      changed it (heap overflow after adding priors), and set_algo never
     assigned
      the new algorithm
     * NcmFit aborts in its default copy_new, run and get_desc; class padding 14
     for
      its four virtual methods; drop the parent copy from the private struct
     * Document the class and its properties; drop commented-out code
     * Tests: /dof and /levmar_set_algo for every fit algorithm

 * Fix NcmFitState least-squares m2lnL and its gradient

     * set_ls stored the Euclidean norm |f| as -2 ln L, so GSL least-squares
     fits
      reported sqrt(-2 ln L) (0.8476 against 0.7184); it is now f.f
     * The gradient 2 J^T f used the Jacobian of the previous call (zero on the
      first), since J was copied last; J is copied first
     * The precision is |2 J^T f| / |f.f|, as the other fitters report it,
     instead of
      the square root of the gradient norm
     * Document the state and set_ls
     * Test: /ncm/fit_state/set_ls/values

 * Fix NcmLikelihood annotations and error text and document the posterior

     * priors_take annotated (transfer full) on @lh instead of @prior
     * leastsquares_f aborted with a message naming priors_leastsquares_f
     * priors_peek_f and priors_peek_m2lnL assert their index
     * Document the posterior in the two prior forms, the m2lnL-v layout and the
      properties; drop the trivial finalize and a commented-out log line
     * Tests: new tests/c/ncm/fit/test_ncm_likelihood.c

 * Make the NcmPriorFlat least-squares form match its -2 ln P and document the
     priors

     * The f-form priors enter -2 ln L as f^2 and serve as least-squares
     residuals.
      NcmPriorFlat returned an f whose square was about e^h0 / 4 times the
      documented -2 ln P near a wall, putting the wall at x0 + 0.965 s; it now
      returns f with f^2 = -2 ln P as documented: e^h0 at a limit, one at half a
      width inside, e^-h0 a width inside
     * NcmPriorGauss doc: -2 ln P = (x - mu)^2 / sigma^2
     * NcmPriorGauss and NcmPriorFlat abort in their default mean instead of
     calling
      a NULL vfunc; drop the parent-instance copies from their private structs
     * NcmPriorFlatParam frees its model_ns and param_name
     * Document the properties and the two prior forms; drop trivial finalizers
     and
      a commented-out line; padding comments
     * Tests: new tests/c/ncm/fit/test_ncm_prior.c
     * Regenerate ncm.pyi

 * Fix NcmFunctionSampleSet domain expansion at a hard limit and the
     adaptive-midpoint messages

     * expand_domain added a second copy of an endpoint already at its hard
     limit,
      and the next spline build aborted on the repeated knot; such a side is now
      not expanded. An empty set or a non-positive endpoint aborts, since the
     steps
      are multiplicative
     * adaptive_midpoint reported "Max iterations reached" when the last allowed
      iteration converged; the message now depends on whether every interval
     passed
     * adaptive_midpoint still marks an interval too short to split as passed,
     and
      now reports how many did in a message
     * Share one routine for the range and |y| maxima tracking (four copies) and
     one
      point constructor; drop the trivial finalize
     * Docs: iterator next/prev on the ends, the old-samples precondition of
      adaptive_midpoint, the positive domain of expand_domain, trimmed residual
     docs
     * Tests: new tests/c/ncm/stats/test_ncm_function_sample_set.c

 * Document NcmStatsDistVKDE and drop its dead code

     * Class doc cut to what the class does, with the theory-page link: local
     scale
      matrices from the k nearest neighbors of the whole sample, the bandwidth
      rule, shrinkage by the mean local matrix; the heading gains its colon and
      the outdated scalar-shrinkage description goes
     * Document local-frac and use-rot-href; fix the set_local_frac and
      set_use_rot_href docs and the header's end guard
     * Theory page: the VKDE bandwidth
     * Drop the trivial finalize, a reset override that only chained up, a
      commented-out pragma and OpenMP debug lines; set_property checks its type
     * Regenerate ncm.pyi

 * Fix the NcmStatsDistKDE LSCV closed form under shrinkage and the
     fixed-covariance reprepare

     * The closed-form least-squares cross-validation of a Gaussian kernel took
     the
      integral of the squared mixture over sample points against centers, which
      differ under center shrinkage (-0.1095 against the exact -0.0753). It is
     now
      the mean over the centers of the mixture at sqrt(2) h, exact for the KDE
     and
      the documented approximation for NcmStatsDistVKDE
     * With cov-type FIXED and anisotropic shrinkage every prepare rebuilt the
      unshrunk factor from the shrunk one, compounding the shrinkage; the factor
      made when the matrix is set is kept and copied at each prepare
     * Class doc cut to what the class does, with the theory-page link; theory
     page:
      the closed form and its VKDE approximation; document the properties and
     the
      cov-type enum; drop commented-out code and the trivial finalize
     * Class padding in the stats headers is 18 minus the virtual methods
      (NcmStatsDist 3, NcmStatsDistKDE 18), with the matching comment
     * Remove the obsolete kde.png and vkde.png call-flow images
     * Tests: kde/gauss/lscv_shrink against GSL, kde/gauss/cov_fixed/reprepare;
     free
      the kernel in cov_fixed/nearPD

 * Fix NcmStatsDist evaluation after a defensive-frac change and document the API

     * Evaluation and sampling read defensive-frac directly while the wide
     kernel
      is built at prepare: raising it from 0 after a prepare made eval and
     sample
      dereference a NULL kernel. They now use the fraction applied by the last
      prepare, as documented
     * Document the accessors, evaluators and kernel getters: eval_m2lnp returns
      -2 ln P, get_Ki's cov_i is h^2 Sigma_i, get_rnorm is the squared NNLS
      residual, split-frac is the fraction of kernel centers, the split is in
      insertion order, add_obs takes @y
     * Test: defensive checks the mixture is unchanged, and sampling works,
     between
      set_defensive_frac and the next prepare

 * Make the NcmStatsDist acceptance objective leave-one-out and recover rejected
     fits

     * The acceptance estimate weighted the kernel points by q/pi with q
     containing
      each point's own kernel, which inflated q/pi where pi is small: the
      importance ESS stayed near 3 of 480 at every bandwidth for a Gaussian
     sample
      in d = 10. Kernel points now use the leave-one-out mixture q_{-k}, with
     the
      defensive component mixed in
     * BOBYQA started among rejected (+inf) trials stopped there reporting
     success,
      so the fit kept its start silently. A fit ending on a non-finite value now
      scans ln over-smooth on 21 points, restarts from the best finite one, and
      aborts when none is finite
     * The LOO cross-validation is least-squares cross-validation of the
     integrated
      squared error, not an AMISE estimate: enum doc, theory page and comments
     * Replace the "Gaussian to 1e-4" nu-ceiling comment by the measured
     deviation;
      drop references to private notes; accept and loo_m2lnp definitions static;
      sample2 pairs kernels by n_kernels; drop a commented-out pragma
     * Theory page: the acceptance estimate, its leave-one-out weights and the
     search
     * Tests: accept/rejected_start, accept/rejected_everywhere; free the kernel
     in
      cv_objectives

 * Fix the NcmStatsDist split fallback and size guard, make center shrinkage exact

     * The 0.9/0.1 fallback compared the kernel centres in range with half of
     all
      observations; with split-frac below 0.5 it always triggered and skipped
     the
      NNLS weight fit. Compare with half the kernel centres, as documented
     * Guard on the kernel centres, whose covariance the shapes use: fewer than
      d + 1 aborts instead of a nearPD-repaired degenerate density
     * Center shrinkage solves A (r C + kappa h^2 <Sigma>) A^T = C with
      r = (n - 1) / n: C is the unbiased covariance while the centres of the
      equal-weight mixture have the 1/n one, so the mixture covariance now
     equals C
      to rounding (it was (1 - a^2/n) C); a can reach sqrt(n / (n - 1))
     * Theory page: the eigenbasis transform the code uses and the factor r
     * Tests: center_shrink checks the mixture's own covariance at 1e-12;
      drop_far_points at split-frac 0.5 and 0.3 checking the weight fit ran;
      split/too_few_kernels; defensive asserts the floor eps K and the far-point
      increase only for Gaussian kernels (it failed on 3 of 30 seeds before);
      free the kernels the split tests create

 * Document the NcmStatsDist properties and correct the auto-kernel scope

     * Document the kernel, N, over-smooth, CV-type, use-threads, split-frac and
      print-fit properties; split-frac is the fraction of kernel centers, taken
     in
      insertion order
     * auto-kernel applies to every cross-validation except NONE, not only
      SPLIT_M2LNP; replace "Gaussian to 1e-4 at nu = 1e4" by the measured
      log-density difference [(chi2 - d)^2 - 2d] / (4 nu)
     * State that the defensive component's c C is its scale matrix
     * Theory page: list the four bandwidth objectives; the kernel stays
     Student-t
      at the nu ceiling
     * Restore the g_return_if_fail in set_property, drop the trivial finalize
     * Regenerate ncm.pyi

 * Document the NcmStatsDistKernel family and construct the Student-t nu default

     * Rewrite the NcmStatsDistKernel, NcmStatsDistKernelGauss and
     NcmStatsDistKernelST
      docs: kernel, normalization and chi2 definitions, the scale matrix and its
      upper Cholesky factor, covariance kappa h^2 Sigma, rule-of-thumb
     bandwidths
     * Correct the Student-t density in the docs to (nu pi)^(d/2) and describe
     nu as
      the degrees of freedom
     * Install the nu property with G_PARAM_CONSTRUCT so its default 3 applies;
     a
      kernel built from properties had nu = 0
     * Drop the empty private struct and property handlers, commented-out code
     and
      trivial dispose/finalize; remove the return in void
     ncm_stats_dist_kernel_sample
     * Test normalization, second moment and the AMISE bandwidth by radial
     quadrature,
      samples against the chi-squared and F laws, the aborting defaults through
     a bare
      subclass, and the nu default
     * Regenerate ncm.pyi

 * Document NcmStatsAcorr more precisely and report the cap from every estimator

     * Report the sample-count cap from ncm_stats_acorr_get_tau_method on zero
     variance
     * Name the configurable AR criterion instead of AICc in the docs
     * Document how much of the series the drift and variance-shift halves cover
     * State Geyer's bound as asymptotic
     * Document the two-sample minimum of ncm_stats_acorr_get_ar_fit and
     get_acov
     * Replace the em-dashes in the class doc
     * Correct the ZERO_VARIANCE row of the autocorrelation theory page to the
     cap

 * Fix NcmStatsVec reset and strided input, and document its row semantics

     * Zero the mean in ncm_stats_vec_reset (a row of 0.3 after 1e17 gave a mean
     of 0)
     * Read strided input vectors through a contiguous copy and abort on a wrong
     length
     * Save each row's weight so enabling the quantiles late skips zero-weight
     rows
     * Abort in ncm_stats_vec_get_cov_matrix unless the type is
     NCM_STATS_VEC_COV
     * Document that reset with rm_saved FALSE keeps rows the statistics no
     longer describe
     * Document heidel_diag's values as Cramér-von Mises distribution values (1
     - p-value)
     * Rewrite the NcmStatsVec docs and remove the dead code
     * Add reset, strided, quantile replay and input-check tests

 * Make the NcmStatsDist2d defaults abort and trim NcmStatsDist2dSpline to what it
     implements

     * Give every NcmStatsDist2d method a default that aborts naming the
     subclass
     * Fix the virtual annotation of ncm_stats_dist2d_eval_marginal_cdf()
     * Define ncm_stats_dist2d_eval_marginal_inv_cdf(), declared but missing
     * Abort in the unimplemented NcmStatsDist2dSpline cdf, marginals and
     inv_cond (returned 0)
     * Evaluate the NcmStatsDist2dSpline pdf as exp(-m2lnp/2) and abort in
     prepare without m2lnp
     * Remove the unused marginal-x property and the constant normalization
     fields
     * Rewrite the NcmStatsDist2d and NcmStatsDist2dSpline docs and regenerate
     the stubs
     * Add test_ncm_stats_dist2d_spline.c and ref/clear coverage in
     test_ncm_generic.c

 * Speed up the NcmStatsDist1dEPDF bandwidth selector and document the estimator

     * Precompute the weighted DCT powers and stop each sum where exp underflows
     (same bandwidths)
     * Remove the dead phat-spline paths, their vectors and the unused inverse
     FFTW plan
     * Remove the unused outliers-threshold property
     * Abort in prepare when there are no observations
     * Mark the bandwidth stale in ncm_stats_dist1d_epdf_set_min() and set_max()
     * Rewrite the NcmStatsDist1dEPDF docs with the estimator and its measured
     bandwidth accuracy
     * Add bandwidth regression, Gaussian AMISE, bounds and empty-observation
     tests

 * Continue NcmStatsDist1dSpline tails to zero and add a density spline

     * Continue m2lnp beyond the knots C1 with a curvature that makes it grow to
     +inf
     * Offset the density by the smallest knot value of m2lnp so it cannot
     underflow
     * Add the density property and ncm_stats_dist1d_spline_new_from_density()
     * Abort in prepare unless exactly one of the m2lnp and density splines is
     set
     * Remove the unused tail-sigma property and the dead PROP_SIZE
     * Refine the NcmStatsDist1d mode between two tied grid points
     * Rewrite the NcmStatsDist1dSpline docs and regenerate the stubs
     * Add test_ncm_stats_dist1d_spline.c

 * Rework the NcmStatsDist1d CDF and its inverse, and add Steffen splines

     * Integrate the CDF in units of p(mode)(xf - xi), in at least 1000 steps
     * Invert the CDF knots with a Steffen spline instead of the 1/p ODE
     * Keep the last knot of a zero-density run so no quantile falls in the gap
     * Bracket the eval_mode Brent refinement by the grid neighbours of the
     minimum
     * Abort through default p, m2lnp and get_current_h methods naming the class
     * Return +inf for m2lnp off a point mass and abort when xf < xi
     * Remove the unused max-prob property and the dead PROP_SIZE
     * Rewrite the NcmStatsDist1d docs
     * Add NCM_SPLINE_GSL_STEFFEN and regenerate enum_nicks.json and the stubs
     * Fix ncm_spline_gsl deriv_nmax for splines whose second derivative jumps
     at knots
     * Add test_ncm_stats_dist1d.c with Gaussian, zero-gap and HSC PDR1 P(z)
     cases
     * Add data/hsc_pdr1_pz_multipeak.bin and
     tests/tools/make_hsc_pz_test_data.py
     * Regenerate the NcGalaxyRedshiftFactorSpline gen and WL resample goldens

 * Fix NcmBootstrap realization handling and document it

     * Read the realization property values; a save/load round trip gave all
     zeros
     * Abort on a realization of the wrong length or with an index >= full-size
     * Sort a copy in ncm_bootstrap_get_sortncomp, leaving the realization
     unchanged
     * Return an empty array from ncm_bootstrap_get_sortncomp at bootstrap-size
     0 (segfaulted)
     * Abort in ncm_bootstrap_get_sortncomp when there is no realization
     * Abort in ncm_bootstrap_resample when full-size is 0 (segfaulted)
     * Abort in ncm_bootstrap_remix when bootstrap-size exceeds full-size
     (returned garbage)
     * Discard the realization when ncm_bootstrap_set_fsize or set_bsize changes
     a size
     * Abort in ncm_data_m2lnL_val when the bootstrap has no realization
     (segfaulted)
     * Evaluate NcmDataset blocks through _ncm_dataset_data_m2lnL_val, with the
     same check
     * Resample zero-size blocks in NCM_DATASET_BSTRAP_TOTAL
     * Remove the unused PROP_SIZE
     * Rewrite the NcmBootstrap docs; set_fsize does not set the bootstrap size
     * Add test_ncm_bootstrap.c and the no-realization and empty-block tests for
     NcmData and NcmDataset

 * Allocate the NcDataSNIACov light-curve covariance only when written

     * Allocate cov_full, cov_full_diag and cov_packed lazily in
     _nc_data_snia_cov_alloc_cov_full
     * Free them in set_size and at the start of the FITS and text loaders,
     resetting cov_full_state
     * Ignore cov-full on catalog version 2 so older serializations still load
     * Return NULL from peek_cov_full and peek_cov_packed when no light-curve
     covariance exists
     * Abort in _prep_to_resample and _prep_to_estimate when no light-curve
     covariance exists
     * Remove the unused inv_cov_mm and inv_cov_mm_LU and the Cholesky inverse
     that filled them
     * Document NcDataSNIACov:cov-full and regenerate the stubs
     * Add test_nc_data_snia_cov with lazy-allocation, version 2 and light-curve
     abort tests

 * Warn at initialization when OMP_NUM_THREADS exceeds OMP_THREAD_LIMIT

     * Warn in ncm_cfg_init_full_ptr when omp_get_max_threads() exceeds
     omp_get_thread_limit()
     * Document the warning in the ncm_cfg_init and ncm_cfg_init_full_ptr docs
     * Add /ncm/cfg/omp_thread_limit, which checks the warning in a subprocess
     with 2 threads and limit 1

 * Passing hard priors to getdist.

 * Compare the data tests' doubles exactly only when bit identical

     * test_ncm_data_gauss_cov replace_cov: compare -2lnL and the least-squares
      norm with the sum of the actual squared residuals, to 1e-13 relative;
      the residuals (m + 1) - m are 1 only up to the rounding of the cosine
      mean, and failed on macOS
     * test_ncm_data: check the Fisher matrix symmetry to 1e-14 relative, as
      the BLAS product need not sum both halves in the same order

 * Document the toy data likelihoods and drop their dead code

     * Give NcmDataFunnel, NcmDataRosenbrock and NcmDataGaussMix2D their
      formulas; say that gaussmix2d reads its point from NcmModelRosenbrock
      and that the length and degrees of freedom are a nominal 10
     * Drop the unused private structs, the empty finalize functions and the
      commented printf lines; fix the broken gaussmix2d doc text and the
      free, clear and constructor docs
     * Abort with a message when gaussmix2d has no NcmModelRosenbrock
     * Add test_ncm_data_toys with exact values
     * Regenerate the stubs (RNG seed blurb from #384)

 * Check NcmCatalog rows and data and report bad rows through GError

     * Report a row index not smaller than the length through the existing
      GError, with the new NCM_CATALOG_ERROR_ROW_OUT_OF_RANGE, in set, get and
      their int and bool forms; they read and wrote past the matrix
     * Abort with a message when the data property gets a matrix with another
      number of columns or is unset
     * Return 0 or FALSE from get_int and get_bool on error instead of casting
      NaN
     * Replace the construction asserts with messages; drop the empty
      finalize; update the docs
     * Regenerate enum_nicks.json for the new error code
     * Test the row errors, replacing the data and the construction messages

 * Give NcmDataDist1d and NcmDataDist2d aborting default methods

     * Default dist1d_m2lnL_val, dist2d_m2lnL_val and inv_pdf abort with a
      message naming the class; a subclass without them crashed through a
      NULL method pointer
     * Abort with a message in get_data when there are no points
     * Document both classes, the per-point likelihood and the inv_pdf
      contract, and their properties; drop the empty constructed and finalize
     * Add test_ncm_data_dist with toy subclasses in one and two dimensions

 * Expose and serialize the NcmDataPoisson bin edges and drop dead code

     * Add the bin-edges property, serialized with the data, and
      ncm_data_poisson_get_bin_edges and ncm_data_poisson_get_bin_range, so a
      subclass can compute a bin's mean from its edges; the edges were stored
      but neither exposed nor serialized
     * Give a new size the edges 0, 1, ..., n (gsl_histogram_alloc leaves them
      undefined) and check that given edges are increasing
     * Remove the begin method and log_Nfac, which computed ln N! for nothing
     * Replace the unchecked inputs with messages: empty or mismatched edges
      and counts, the mean property, the mean vector and the bin index
     * Document the deviance form of -2lnL, the properties and the accessors;
      drop the stray parent_instance, the empty constructed and finalize and
      the unused PROP_KNOTS
     * Add test_ncm_data_poisson with a toy subclass that uses the bin edges

 * Pass NcmDataGaussCovMVND covariance updates through set_cov

     * set_cov_mean and gen_cov_mean wrote into the covariance in place, so
      after an evaluation the old Cholesky factor stayed in use (5 instead of
      20); they now go through ncm_data_gauss_cov_set_cov
     * gen reports one realization when there is no bound (it reported 0)
     * Document the class, the constructors and stats_vec's maxiter, which
      counts consecutive rejections
     * Replace bare asserts with messages; drop the empty finalize, the empty
      NUMCOSMO_GIR_SCAN block and a dead commented line
     * Test the replaced covariance, the regenerated one and gen's count

 * Use the profiled likelihood in the NcmDataGaussDiag w-mean Fisher matrix

     * With w-mean, whiten by 1/sigma and project out the weighted constant
      direction in inv_cov_UH and inv_cov_Uf, so the Fisher matrix and bias
      are those of the likelihood profiled over the offset, as least squares
      already was; before, a constant shift of the mean had full information
      (F_aa 261 instead of 0, F_bb 118.4 instead of 53.1), which made
      NcDataDistMu Fisher forecasts overconfident
     * Mark the weights stale whenever sigma_func reports a change, on every
      path, and when sigma is set or the size changes; setting sigma kept the
      old weights and evaluating after a resize segfaulted
     * Document the class, the profiled likelihood, the virtual method
      contracts and the properties; drop the empty constructed and finalize;
      fix the Fisher bootstrap message
     * Add test_ncm_data_gauss_diag

 * Keep the NcmDataGaussCov Cholesky factor in step with its covariance

     * Invalidate the factor in set_cov, the cov property and set_size, and
      whenever cov_func reports a change; after replacing the covariance
      m2lnL and least squares kept the old one (5 instead of 20), and
      resampling after a resize segfaulted
     * Rebuild the factor in one place, with a message naming the data when
      the covariance is not positive definite; remove the unreachable +inf
      branch
     * Check the matrix size in set_cov
     * Count one ln 2pi per bootstrap draw in lnNorma2_bs, not per data point
     * Document the class, the mean_func and cov_func contracts and the
      properties; drop the stray parent_instance and the empty constructed
      and finalize; give bulk_resample and the least-squares bootstrap check
      messages
     * Test the replaced covariance, resize, the bootstrap normalization and
      the new messages

 * Keep the NcmDataGauss Cholesky factor in step with its inverse covariance

     * Mark the factor stale whenever inv_cov_func reports a change, on every
      path; m2lnL consumed the change and later least-squares, bootstrap,
      resampling and Fisher calls used the old factor (2.0 vs 0.5)
     * Mark it stale when inv-cov is set (least squares gave 5 for 20) and on
      set_size (resampling after a resize segfaulted)
     * Rebuild a stale factor in one place, with a message naming the data
      when the inverse covariance is not positive definite
     * Document the class, the mean_func and inv_cov_func contracts and the
      properties; fix an error label and drop the empty finalize
     * Add test_ncm_data_gauss with a parameter-dependent covariance

 * Prepare every NcmDataset member before the mean and Fisher paths

     * Prepare all members before evaluating any in ncm_dataset_mean_vector,
      fisher_matrix and fisher_matrix_bias, as m2lnL and least squares do;
      members sharing a resource were evaluated with only the requirements
      of those prepared before them
     * Make the dataset Fisher IM and delta_theta (inout): reuse a matrix
      passed in (it was ignored and leaked), check sizes and f_true with
      messages, and clear them without free parameters (it segfaulted)
     * Keep the bootstrap type and probabilities in ncm_dataset_copy
     * Rewrite the class doc; fix the property doc titles, parameter docs,
      messages and constructor annotations; drop the empty finalize
     * Add test_ncm_dataset with members sharing a resource

 * Evaluate the NcmData Fisher covariance and bias at the point itself

     * Restore the free parameters right after the numerical derivative in
      fisher_matrix and fisher_matrix_bias; the covariance and the bias mean
      were evaluated at the derivative's last step (a zero shift came out as
      [-0.54, -3.06])
     * Remove the ncm_data_sigma_vector declaration, which had no body
     * Annotate the Fisher IM and delta_theta as (inout), free and clear them
      without free parameters, and check their sizes with messages
     * Run begin again after every resample, since the data changed
     * Make _ncm_data_diff_f static
     * Give prepare on uninitialized data and the bootstrap checks messages
     * Document the evaluation order and fix the property docs
     * Add test_ncm_data

 * Serialize long and ulong properties with their value

     * Transform G_TYPE_LONG and G_TYPE_ULONG values to int64 and uint64 before
     g_dbus_gvalue_to_gvariant(), which reads a 64-bit variant type with the
     int64/uint64 getters and wrote 0 for NcmRNG:seed
     * Test the NcmRNG:seed round trip in test_serialize.py
     * Make a deserialized NcmRNG resume its state and keep its seed
     * Record the seed without reseeding when the seed property is applied after
     the state property, so the loaded state wins in either property order
     * Do not reseed in constructed when only a state was given
     * Factor _ncm_rng_record_seed() out of ncm_rng_set_seed(), whose API
     behavior is unchanged
     * Add /ncm/rng/property_order to test_ncm_rng.c: both property orders,
     state alone, the API reseed, and a serialization round trip
     * Add test_serialize_rng_resumes_stream to test_serialize.py

 * Make two new tests hold on every CI platform

     * numdiff_fparams tests: require the gradient to equal a direct NcmDiff
      call on the same function and NcmDiff's error estimate to bound the
      error against the exact gradient; the 1e-12 bound came from a
      measurement on one machine and failed on three CI jobs (1.01e-12)
     * eval_forms: add an absolute floor scaled by the cancelled terms to the
      second log-derivative check; it vanishes for a power law and a fused
      multiply-subtract on macOS moved it by 2e-16

 * Fix the toy models' docs and names and check the MVND dimension

     * Describe NcmModelFunnel, NcmModelRosenbrock and NcmModelMVND; the first
      two carried the MVND doc
     * Use a descriptive name and a short nick, as the other models do:
      Funnel, Rosenbrock and MVND were registered with swapped strings
     * Abort at construction when NcmModelMVND dim differs from mu-length;
      dim=3 alone gave a mean of length 1
     * Check the vector length in ncm_model_mvnd_mean
     * Drop the placeholder private structs, empty handlers and empty
      NUMCOSMO_GIR_SCAN blocks
     * Set mu-length in the serialize test that built an MVND from dim alone
     * Add test_ncm_model_toys

 * Let NcmModelBuilder extend registered models and re-create types safely

     * Register the new model id only when the parent has none (model_id < 0),
      not whenever main_model_id is -1; building on any registered model
      aborted
     * Place the builder's parameters after the parent's at absolute ids; a
      parent with parameters hit "scalar parameter 0 is already set"
     * create returns an existing type with the same parent and parameters and
      aborts on any other name clash; it returned the invalid type 0
     * Check at construction that the parent exists and is a NcmModel; a
      non-model parent segfaulted
     * Abort with a message when adding parameters after create; drop the
      asserts on construct-only properties and the empty instance_init
     * Add ncm_model_builder_add_vparams, _get_vparams, _free and _clear
     * Document the type and its properties
     * Add test_ncm_model_builder

 * Make NcmMSetFuncList lookups unambiguous and select bindable

     * Register NcmMSetFuncListStruct as a boxed type; select elements no
      longer dangle in bindings after the array is freed
     * has_full_name matches the namespace exactly, as ncm_mset_func_list_new
     * new_ns_name and has_ns_name take the exact namespace first and abort
      when a namespace prefix matches the name in several namespaces
     * Release the previous object when the object property is replaced
     * Clear messages for a missing name, a missing or malformed full name and
      a missing or mistyped object; full-name defaults to NULL
     * Add ncm_mset_func_list_ref, _free and _clear
     * Document the class and its properties; drop the stray parent_instance
      and the empty finalize
     * Test registration, lookups, the object property, select and the errors

 * Make NcmMSetFunc eval-x a checked copy and guard scalar evaluations

     * Copy the eval-x vector (stride-safe) instead of sharing it, check its
      length against nvariables and rebuild the unique names on set
     * ncm_mset_func_set_eval_x takes const x, aborts with a message on a wrong
      length and drops the side-effect peek_desc call
     * Abort in eval_nvar/eval0 for functions with dim != 1, and in eval1 and
      eval_vector unless dim == 1 and nvar <= 1; the scalar evaluations wrote
      dim values into a single double
     * Check that eval_vector's argument and value vectors have the same length
     * Document the type, the argument rule and the three properties; drop the
      stray parent_instance from the private struct
     * Name each MVNDMean in test_fit_mc.py instead of an uninitialised eval_x
     * Test the property copy, stride, clearing and the new aborts

 * Add ncm_mset_func_eval_array and make NcmMSetFunc1 bindable only

     * Add ncm_mset_func_eval_array: returns the dim values as a GArray,
      follows the argument rule and checks the number of arguments
     * Document that ncm_mset_func_eval does not return values to bindings
     * Remove ncm_mset_func1_eval1: it bypassed the argument rule and shadowed
      ncm_mset_func_eval1 in bindings; annotate the NcmMSetFunc1N typedef so
      the eval1 virtual function keeps its element types
     * Make NcmMSetFunc1 abstract with a default eval1 that aborts, and abort
      with a message when eval1 returns the wrong number of values
     * Drop the placeholder private struct and the empty finalize
     * Test eval_array, a C NcmMSetFunc1 subclass and the Python bindings

 * Annotate NcmMSetFunc numdiff out as inout and fix its leak

     * Mark @out of ncm_mset_func_numdiff_fparams as (inout) (allow-none): the
      function reads *out to reuse a vector, which (out) does not allow
     * Free the vector wrapping the saved free parameters, leaked on every call
     * Document that a NULL *out allocates and a non-NULL *out is overwritten
     * Test both paths, the gradient and the restored free parameters

 * Make NcmMSetFunc unique names exact and collision free

     * Write each eval_x component in the shortest locale-independent form that
      reads back to the same double, instead of a locale-dependent %.15g
     * Separate the components by ',' in the unique symbol and by '_' in the
      unique name; (1, 2) and (12) no longer share a name
     * Map '-' to 'm' and '.' to 'p' and drop the exponent '+' in the unique
      name, so 1 and -1, or 1e-05 and 1e+05, no longer collide
     * Rebuild the unique name and symbol in ncm_mset_func_set_meta
     * Add test_ncm_mset_func with the exact names for these cases

 * Accept stacked ids in NcmMSet __setitem__ and validate stack positions

     * Split an integer id into base id and stack position in __setitem__; ids
     with
      a stack position were reported as a mismatch
     * Parse stack positions in peek_by_name and __setitem__ with one strict
     helper:
      an empty or signed position was read as 0 and replaced that model
     * Say which position, name, key type or model is wrong in the
     stack-position,
      key-type and mismatch errors; update the Python patterns
     * Document the accepted keys of __getitem__/__setitem__ and peek_by_name
     * Add a test for getting and setting by id and by name, invalid positions
     and
      a type mismatch

 * Keep the free-parameter map across NcmMSet save/load and fix load cleanup

     * Write the free-parameter map, in order, as fmap in the NcmMSet group and
      restore it on load; files without it get the map prepared as before
     * Free the half-built mset on every ncm_mset_load error
     * Remove the submodel names ncm_mset_load registered in the caller's
      serializer; a second load with the same serializer failed
     * Free the objects held when setting a loaded host fails
     * Fix a use-after-free in ncm_serialize_remove_ser and a leak of the
      instance name in ncm_serialize_from_name_params without autosave
     * Document the mset file format and the load behavior; shorten the
      submodel-unset comment
     * Add tests for the map order, repeated loads, remove_ser and named objects
     * Assert the loaded mset in xcdm_no_perturbations.py; regenerate the stubs

 * Remove ncm_mset_fparam_validate_all and document the free-parameter map

     * Remove ncm_mset_fparam_validate_all and its scratch vector; it changed
     every
      model's pkey to restore the values it tested
     * Check the model and parameter index in ncm_mset_fparam_get_fpi, which
     read
      past the array for an index beyond the model
     * Drop the return of an expression from the void ncm_mset_fparam_set
     * Explain the free-parameter map as a snapshot respected until renewed, and
      the two directions between the map and the fit types
     * Document the NULL and error cases of the name lookups, the vector sizes,
      and fparam_len against fparams_len
     * Add tests for the free-parameter lookups and the fpi range check
     * Regenerate the stubs

 * Make NcmMSet fit types follow a set fmap and fix two printers

     * Set every parameter FIXED or FREE from the map in
      ncm_mset_param_set_ftype_from_fmap and keep the map and its order; it
      left other parameters free and rebuilt the map in model order
     * Remove the space after the header line in ncm_mset_params_pretty_print
     * Count models without parameters in ncm_mset_max_model_nick
     * Document the unchecked model id of the pass-through accessors, where
      pretty_log and params_pretty_print take FREE/FIXED from, ncm_mset_cmp and
      the -1 of ncm_mset_param_get_ftype
     * Add tests for set_fmap with update_models, max_model_nick and the header
     * Regenerate the stubs

 * Keep the NcmMSet free-parameter map consistent and bound stack positions

     * Rebuild the list of models to update in ncm_mset_set_fmap; after moving
     the
      free parameters to another model, setting them left its pkey unchanged
     * Resolve and check every ncm_mset_set_fmap name, repeated ones included,
      before changing the map, so an error keeps the previous one
     * Reject stack positions of NCM_MSET_MAX_STACKSIZE or more in set_pos and
     push;
      position 1000 addressed the next model class
     * Assert a non-negative id in ncm_mset_get_ns_by_id/get_type_by_id
     * Resolve subclasses in ncm_mset_get_id_by_type through
     ncm_model_id_by_type
     * Free the item, with its model reference, when a model is not stackable
     * Document exists/exists_pos, set_fmap, get_fmap and the stacking rules
     * Add tests for the fmap update, a failed fmap, the stack bound, the id
      lookups and a negative id
     * Regenerate the stubs

 * Remove a host's submodels from NcmMSet and reject unknown namespaces

     * Make ncm_mset_remove remove the host's own submodels with it; a
      serialization round trip dropped the ones left behind
     * Return NULL from ncm_mset_peek_by_name for an unregistered namespace; -1
      plus a stack position addressed another model
     * Move the removal out of g_assert in ncm_mset_remove
     * Name ncm_mset_split_full_name in its error and leave its outputs NULL on
      error; update the three Python message patterns
     * Check the namespace and description before registering a model id
     * Remove the stray parent_instance and commented-out code; document the
      properties, the submodel handling and the nullable lookups
     * Check the peeked model in numcosmo_py/app/generate.py
     * Add tests for unknown names, host removal and split_full_name
     * Regenerate the stubs

 * Name the current host in the NcmModel cross-host error

     * Report the type of the host the submodel is attached to, not of the new
      one, when attaching it to a second host
     * Make ncm_model___getitem__/__setitem__ call
     ncm_model_param_get/set_by_name
      and take const names; update the test_model.py message patterns
     * Rewrite the host_wr comment: submodels are set once and never replaced
     * Document the ncm_model_add_submodel contract and the submodel queries
     * Test the cross-host error with hosts of different types
     * Regenerate the stubs

 * Fix NcmModel reparametrization, description and name lookup defects

     * Make ncm_model_set_reparam with NULL a no-op without a reparametrization
      and an error with one; the removal leaked the reparam's vector
     * Make ncm_model_is_equal symmetric in the presence of a reparametrization
     * Remove the undeclared ncm_model_get_reparam
     * Stop ncm_model_param_set_default from marking the description modified
     * Check the current parameters in ncm_model_param_finite and _params_finite
     * Return -1 from ncm_model_id_by_type on error
     * Free the keys and values of the ncm_model_param_get_desc table and copy
      its strings
     * Return an error instead of aborting when an unqualified name was renamed
      by a submodel's reparametrization
     * Stop ncm_model_param_set_desc at the first invalid key and document the
      accepted keys
     * Check the vector length in ncm_model_orig_vparam_set_vector
     * Fix six gtk-doc blocks named after other functions; take const names in
      get_desc/set_desc
     * Include nc_hireion_camb_reparam_tau.h in numcosmo.h
     * Remove the stray parent_instance and commented-out code; document
     NcmModel,
      its properties and the class builders
     * Add tests for all of the above
     * Regenerate the stubs

 * Keep NcmReparam names found after deserialization and reject singular T

     * Rebuild the name table when params-desc is set, so a deserialized
      reparametrization finds its parameters by name
     * Abort when a description's name already describes another parameter
     * Remove the old name outside g_assert in ncm_reparam_set_param_desc
     * Check the LU factorization and solves in NcmReparamLinear and abort on a
      singular matrix
     * Remove the stray parent_instance and an unused include; return
      G_MAXUINT from ncm_reparam_index_from_name when not found
     * Rewrite the NcmReparam and NcmReparamLinear docs
     * Add test_ncm_reparam.c
     * Regenerate the stubs

 * Fix NcmModelCtrl submodel flags and NcmVParam component setters

     * Fold ncm_model_ctrl_update and _model_update into one helper that records
      the main model directly, so a submodel seen first or replaced with its
      host reports a change
     * Take the new reference before releasing the old one in
      ncm_vparam_set_sparam (use-after-free on the held component)
     * Make ncm_vparam_set_fit_type set the fit type; it set the default value
     * Free the old NcmSParam:symbol when it is set again
     * Drop the nc_hicosmo.h include from ncm_sparam.c; default abstol 0
     * Give NcmSParam:scale a positive minimum and a default of 1, matching
      ncm_sparam_set_scale
     * Remove the unused PROP_LEN, take const strings in
      ncm_vparam_set_sparam_full, rename the ctrl property functions
     * Rewrite the NcmModelCtrl, NcmSParam and NcmVParam docs
     * Add the ctrl switch test, the sparam strings/copy test and
      test_ncm_vparam.c; the first update now reports the submodels
     * Regenerate the stubs

 * Accept a NULL model in NcmPowspecFilter and calibrate it at zi

     * Accept a NULL model in prepare and prepare_if_needed; prepare passed it
      to ncm_model_ctrl_update, which dereferences it
     * Calibrate the bias and the k knots at zi instead of z = 0
     * Remove a dead get_nknots call
     * Document the valid (r, z) range of the evaluation functions, which do not
      check their arguments, and rewrite the class, property and eval docs
     * Port test_powspec_filter.py to test_ncm_powspec_filter.c: power-law
      log-derivatives, the Gaussian closed form, the Arb top-hat table, the eval
      forms, nderivs, settings, the r range and a table starting at z = 0.5
     * Regenerate the stubs

 * Continue NcmPowspecCorr3d into the FFTLog padding and fix its calibration

     * Use the smooth padding with the best bias and a fixed starting size, as
      NcmPowspecFilter does; zero padding left r = 100 at -1.1e-4 at any size
     * Calibrate the k knots at zi instead of z = 0
     * Accept a NULL model in prepare and prepare_if_needed
     * Force the model ctrl in set_reltol and set_reltol_z, which
      prepare_if_needed ignored
     * Free the rows taken on the already-calibrated branch of prepare
     * Declare GObject as the parent in the header, as registered
     * Abort in eval_xi outside the grid; get_r_min/get_r_max report the grid
      in use and abort before prepare
     * Remove the unused constructed field, a dead get_nknots call and the
      commented-out printf blocks; rewrite the docs
     * Add test_ncm_powspec_corr3d.c: calibration, the quadrature for
      r >= 100/kmax, D(z)^2 in z, recalibration, the grid ends and the aborts
     * Compare the corr3d test in test_nc_powspec.c from r = 100/kmax, relative
      to the peak of |xi|
     * Regenerate the stubs

 * Continue NcmPowspecSpline2d smoothly and abort outside its z range

     * Continue ln P linearly in ln k with the spline's value and slope at the
      nearest end, replacing the k^3 low end and the -5 dlnk^2 high end
     * Abort in eval for z outside the table and in prepare when zi or zf was
      required beyond it
     * Drop the get_spline_2d override, which returned ln P on (z, ln k)
     * Reference the new table before releasing the old one in set_spline2d
      and force the next prepare
     * Make get_nknots static and remove the unused includes
     * Document the table, the continuation, the property and the constructor
     * Port the Python tests to test_ncm_powspec_spline2d.c and delete them

 * Review NcmPowspec docs and test it against Arb integrals

     * Rewrite the NcmPowspec gtk-doc: correct the <delta delta*> definition,
      units, Returns: lines, the sproj formula, and the concurrency limit
     * Remove the stray GObject parent_instance from NcmPowspecPrivate
     * Set NcmPowspecClass padding to 11 (7 vfuncs), breaking ABI
     * Limit LCOV_EXCL to the prepare/eval g_error stubs; the default
      derivatives are used by halofit, ml_spline, CBE and diemer15
     * Add an --integrals mode to ncm_powspec_analytic_arb: sigma_R^2, xi(r)
      and the two-sphere C_ell of NcmPowspecAnalytic over [k_lo, k_hi]
     * Move sph_bessel from xcor_window_arb.h to the shared sph_bessel_arb.h
     * Add make_powspec_analytic_truth_table.py and the Arb table
      data/truth_tables/powspec/ncm_powspec_analytic_integrals.bin
      (k in [1e-6, 1e2], relative radius below 1e-20)
     * Add test_ncm_powspec.c: the integrals at reltol 1e-9 against the
      table, the default deriv_z/deriv_k against the closed form, the
      default get_spline_2d against its reltol, eval_vec, and the range
      setters and properties

 * Fix build and test warnings found by the PR #382 lanes

     * Format the NcmSphereMap FFTW plan key with G_GINT64_FORMAT (%ld is
      wrong where gint64 is long long)
     * Document NcmLaurentSeries in one block, and match the header
      parameter names of ncm_vector_hypot and ncm_pln1d_set/get_order to
      their docs; regenerate the stubs
     * Capture the NcmFitESMCMC FULL log in the esmcmc test instead of
      printing it into the TAP stream
     * Close the xcor view figure when it is not shown
     * Assert the deprecation warning of --auto-kernel in its test

 * Fix the pip and macOS lanes of PR #382

     * Skip the healpy comparisons at run time instead of at collection: the
      file is the whole sphere_map shard, and a shard with nothing collected
      makes pytest exit with status 5
     * Drop the sphere_map pytest line from the macOS examples step; meson
      runs the shard in Check NumCosmo
     * Keep G_GINT64_FORMAT on the last line of a split format string in
      ncm_sphere_nn.c, which uncrustify 0.78.1 otherwise misaligns
     * Relax the FFTLog bias truth bound to 5e-13 of the peak (measured 2.4e-13
      on macOS arm64, 5.5e-14 on x86-64 Linux)
     * Compare Chebyshev eval_x/deriv_x with their definitions in units of the
      values, with the chain rule in the library's order
     * Build the expected NcmTimer string from the measured time: a 2 ms sleep
      can exceed 10 ms on a busy runner
     * Compare the Levin Dirichlet fallback on the scale of the integral of
      |F j_l|, not relative to results that are rounding noise

 * Remove a stray LCOV_EXCL_START in NcmSpline2dSpline

     * The derivative stubs became real methods in 69cf7ec7 and lost their
      LCOV_EXCL_STOP, leaving an unmatched START that makes lcov abort

 * Test NcmSphereMap in C against frozen healpy truth tables

     * Add tests/tools/make_sphere_healpy_truth_table.py, writing healpy's
      pixel indices and centres, pixels of directions, map2alm/anafast/alm2map
      of seeded maps, a healpy-written NESTED FITS map and four FITS header
      fixtures to data/truth_tables/sphere (about 200 KB)
     * Check pixels, ang2pix, transforms, cross spectra and FITS files against
      those tables in test_ncm_sphere_map.c, with vectors, update_Cl, FITS
      options and header traps
     * Keep ten live healpy comparisons in test_sphere_map_healpy.py
     * Delete test_sphere_map.py and test_sphere_map_basic.py, now covered in C
     * List the new tables and their consumers in data/truth_tables/README.md

 * Check NcmSphereMap C_l, pixel access and noise, and test C(theta)

     * set_lmax clears the C_l flag: after a new lmax calc_Ctheta returned a
     zero C(theta) without complaint
     * set_Cls aborts on a vector shorter than lmax + 1 (it read past its end)
     and on lmax = 0
     * get_pix takes a gint64 index and checks it; add_noise and set_map loop
     over gint64; set_map checks the array
     * Document calc_Ctheta (formula, range, reltol), set_Cls and get_pix
     * C tests: C(theta) against the GSL Legendre sum, the noise moments, and
     traps for the new checks

 * Document NcmSphereMap alm2map and test it against healpy at any lmax

     * alm2map: document the synthesis, the RING order of the result and the
     folding above a ring's Nyquist frequency
     * Abort with a message on a zero nside; drop dead commented code
     * Test random a_lm against healpy's alm2map at lmax 2 nside, 3 nside - 1
     and 4 nside (measured 1.7e-13)

 * Document NcmSphereMap map2alm and check its arguments

     * prepare_alm: document the quadrature (healpy's use_weights=False), the
     iterations, the healpy agreement (1e-14) and the switch to RING order
     * prepare_alm and alm2map abort on lmax = 0; they warned and returned,
     alm2map's warning naming prepare_alm
     * get_alm, set_alm and get_Cl abort on (l, m) out of range; document l <=
     lmax and when C_l are current
     * compute_cross_Cl aborts with messages on mismatched lmax or nside
     * Remove the _NCM_SPHERE_MAP_MEASURE timing code and its NcmTimer, and dead
     code in the block file; describe the blocked scheme there
     * C traps for the new checks

 * Make NcmSphereMap FITS I/O lossless and readable by and from healpy

     * save_fits wrote a single-precision column, changing every pixel by up to
     6e-8 on a round trip; write doubles
     * load_fits read one pixel per row only and aborted on healpy's own files
     (1024 per row); read the elements across rows
     * A missing ORDERING warned "assuming RING" and then aborted on an unset
     buffer; use RING
     * Reject partial-sky maps and other PIXTYPEs with messages; write PIXTYPE;
     drop the duplicate INDXSCHM
     * A NULL signal name reads the first column; document the format, the
     catalog loader and its units
     * Tests: exact round trip in C; healpy interoperability both ways and the
     rejections in Python
     * NcmSphereMap forced HAVE_FFTW3F off and its float branches did not build
     (fftwf_alloc_real returned double *); keep the double code only
     * Drop the dead no-FFTW branches of NcmSphereMap: FFTW3 is mandatory
     * ncm_cfg: remove the float wisdom load and save, the float time limit and
     ncm_cfg_fftwf_plan_destroy
     * meson: remove the optional fftw3f dependency and HAVE_FFTW3F
     * Fixed-seed map2alm, C_l and alm2map outputs are bit-identical before and
     after

 * Compute NcmSphereMap pixel centres without cancellation near the poles

     * Cap centres came from acos (1 - t^2/3nside^2) and sqrt (1 - z^2), losing
     precision as nside^2 (5.8e-11 at nside 1024); compute z and sin(theta)
     exactly and theta = atan2 (sin, z), as in NcmTriVec
     * vec2pix takes sin(theta) from hypot (x, y) / |v|
     * Factor the index to (ring, position, width) code shared by the four
     pix2ang/pix2vec functions
     * Document angle units, ranges and conventions of the conversions
     * NcmSphereNN get: theta = atan2 (hypot (x, y), z); acos (z / r) lost 4e-4
     at theta = 1e-7
     * Tests: cap centres against 2 asin (t / (sqrt(6) nside)) and NcmSphereNN
     get near the poles

 * Fix NcmSphereMap reordering precision and check its indices

     * set_order copied pixels through a gfloat, rounding double maps at 6e-8 on
     every change of ordering
     * Class doc: an independent implementation of the HEALPix pixelization from
     Gorski et al. (2005), matching healpy; document the properties and get_lmax
     * Abort with messages on a non-power-of-two nside and on pixel or ring
     indices out of range (negatives included)
     * Drop the undefined get_nsmap declaration, a duplicate include and
     commented-out code
     * Tests: bit-identical order round trip and the index traps in C; healpy
     comparison at nside 1, 2 and 4

 * Make NcmSphereNN searches safe and document them

     * Class doc described HEALPix; describe the k-d tree search, indices,
     squared chord distances and the rebuild requirement
     * Abort on a search before a rebuild or after new inserts (it dereferenced
     an empty result or missed points)
     * Abort on k outside [1, n] in all searches (single searches returned
     fewer, the batch one aborted), on get out of range and on mismatched array
     lengths
     * Share one result loop between the three searches; drop PROP_NOBJS, the
     empty dispose and unused includes
     * C tests: round trip, brute-force kNN, agreement of the searches, rebuild,
     repeated points, dump and the traps; drop the Python tests

 * Validate and document the rectangular sky footprint

     * Density: document it as the density in dra ddec (degrees), with the
     cos(dec) factor, not per unit solid angle; give the formula
     * Abort on declination limits outside [-90, 90], non-increasing limits or a
     right ascension span outside (0, 360]
     * Document that right ascension is compared as given, without wrapping at 0
     or 360
     * ra-lim and dec-lim are construct properties; a missing limit spans the
     whole sphere in that coordinate
     * C tests: exact density, sampling moments within 4 sigma, the seed-123
     golden draw, the default and the traps; drop the redundant Python tests
     * Fix the class padding comment

 * Place the sbessel_j output grid at the first maximum of j_l

     * set_best_lnr0/lnk0: k0 r0 = x*, the first maximum of j_l by Newton on
     j_l' (1 for l = 0); the docs claimed this rule but the code balanced the
     kernel's size at the ends of the t range
     * Measured for F = k^-1/2 on [1e-4, 1e4]: grid within 1e-6 of the exact
     transform 91.9% vs 77.8% at l = 10, 96.7% vs 74.4% at l = 50
     * Test x* against mpmath for l = 1..200 and the grid coverage at l = 10
     * corr3d: document its k0 r0 = 1 instead of citing the wrong rule
     * Drop unused ARB and gsl_sf_trig includes from sbessel_j

 * Refresh FFTLog docs and tests after the taper

     * Theory page: tapered-padding numbers (unbiased castro stall 2.4e-6, 1e-11
     at 8000 knots), N^-3 down to a roundoff floor with the EH figures
     * Theory page: rounded padding split and the resulting L_T variation;
     midpoint case of the bias choice
     * Theory page: "Reach of the continuation": what the tapered padding leaves
     out, when it matters, and the levers to extend it with their measured
     trade-offs
     * Filter docs: sigma^2 near 1/k_max and 1/k_min is an extrapolation; the r
     range includes those edges
     * no-ringing now shifts the output grid; calibration tolerance is relative
     to the value plus the peak
     * Kernels: document Lk as the length of the fundamental interval
     * Drop the empty ncm_fftlog_array_pos macro and the unused fftw_alloc_real
     * smooth_padding/power_law: measured 6.2e-10 and 4.4e-10; bound tightened
     from 1e-6 to 1e-8

 * Remove NcmFftlogSBesselJLJM and NcmPowspecSphereProj

     * No library code, Python module or test uses them; the only caller was
     examples/example_sphere_proj.py
     * Archived with their tests, the example and the math in NumCosmoDevNotes
     (archive/sphere_proj_jljm)
     * Drop them from the build, numcosmo-math.h, the type registrations in
     ncm_cfg.c and the stubs
     * Drop the jljm fixture, its traps and the generic ref/free test

 * Taper the FFTLog padding ends and enforce the filter's z-knot limit

     * Replace the blend of the two end continuations by a taper of each to zero
     over [0.4 h, 0.8 h], fixed in ln k
     * The blend mixed values of f from both ends, whose k^b differ by e^(b
     L_T): 5.7e-2 of the peak with the chosen bias at padding 0.3
     * A fractional padding rounds the period differently per size; the taper no
     longer moves with it
     * Document the image floor e^(-(b - b_min) L_T) and its dependence on the
     rounded period
     * Pass max-z-knots to the redshift spline and abort the prepare when the
     grid would exceed it
     * Tests: fractional-padding refinement, short-padding calibration trap,
     max-z-knots trap; negative-slope bounds from the new tail measurement

 * Add a transparent bias to FFTLog and make its grid path-independent

     * NcmFftlog:bias: transform F k^-b, scale the output by r^-(1+b),
     derivative factor -(1 + b + a)
     * get_bias_range vfunc (default: bias 0 only); ncm_fftlog_get_end_slopes
     and ncm_fftlog_get_best_bias
     * Tophat coefficients as a Gamma ratio regular on -1 < b < 3; drop the
     unused ACB branch
     * Gaussian and spherical Bessel kernels take the bias; remove the q kernel
     power from sbessel_j and jljm
     * NcmPowspecFilter chooses the bias from the table's end slopes: castro
     row24 converges again
     * Keep the no-ringing shift out of the requested lnr0; restart every filter
     calibration at 100 knots
     * Filter: setters invalidate the calibration; knot-limit getters; r range
     is the calibrated grid
     * Calibration: size helper, max-n checked before allocating, non-finite
     transforms abort
     * Smooth padding: blend fixed in ln k, cut on the growth of F k; noring
     only for even full sizes
     * Tests: bias truth for the three kernels, castro-slope convergence, bias
     rules and traps

 * Set explicit filter tolerances in the castro and golden tests

     * Lower the castro power-law filter reltol to 1e-7, the double-precision
     floor for slope -0.2 is 2e-8
     * Set the golden tests' filter reltol to 1e-6 explicitly
     * Regenerate nc_data_cluster_ncount_golden_seed0.bin with the smooth
     padding

 * Calibrate the power spectrum filter with the smooth padding

     * Turn on the smooth padding in NcmPowspecFilter, whose k^2 P(k) does not
     vanish at the table ends
     * Pass the halofit reltol to its variance filter instead of the filter
     default
     * Set the CCL fixture's tophat filter reltol to a tenth of prec
     * Document that ncm_powspec_var_tophat_R integrates the table only
     * Compare filter and quadrature in test_nc_powspec_filter_tophat two
     decades inside both ends

 * Continue the input smoothly into the FFTLog padding

     * Replace the smooth padding by the power law of each end, anchored at the
     end with value and log-slope from a cubic through the four nearest knots
     * Cut a continuation along which F k grows by a Gaussian in ln k of width
     the inverse of that log-slope
     * Join the two continuations by a C-infinity partition over the middle
     fifth of the padding
     * Fix the phase for an odd full size, which shifted the output by one knot
     * Drop the smooth-padding-scale property and accessors
     * Remove the Nyquist-pair forcing in _ncm_fftlog_eval
     * Make ncm_fftlog_get_Ym clear the prepared flag
     * Enforce a positive period and abort on F <= 0 at the padded ends
     * Compute one transform per size in ncm_fftlog_calibrate_size_gsl and abort
     past max-n
     * Rewrite the class and theory-page docs for the padding as a fraction and
     the continuation
     * Add tests: tophatwin2 truth vs mpmath, power-law law with smooth padding,
     odd full size, get_Ym keeps eval, padding and period traps

 * Review sbessel ODE solver docs

     * Use F(x) for the forcing, as the theory page and Levin do
     * State what operators copy from the solver (tolerance, default constraint,
     tau floor factor) and what reconfiguring resets
     * Describe the first two right-hand side entries as the constraint
     functionals' values, u(a) and u(b) only under Dirichlet
     * Call the third endpoint slot the derivatives' roundoff bound, the value
     get_last_deriv_error returns
     * Fix internal docs that no longer matched their signatures: parameters,
     names, void returns, column counts
     * Describe constraint rows rather than boundary conditions; correct the
     NcmSBesselOdeSolverRow description

 * Review Levin internals; count guard fallbacks once

     * Exclude tau fallbacks from n_locked_eligible_solves, which counted a
     guard rejection again; test it
     * Fix the ell-cache-max doc: a multipole range above it aborts
     * Correct the comments on the forcing fit (x F(x) for order 0), the
     null-panel skip (the Bessel array's threshold, not underflow), the wrapper
     functions, the no-knot path and the dead-junction probe
     * Document _integrate_panel and move the constraint-rule comment to its
     function
     * Merge a duplicated comment and remove a dead variable

 * Review Levin public docs

     * Remove its Simpson path, align defaults, restore K(chi, k)
     * Name the integrand the radial kernel K(chi, k) in the base and GL docs;
     F(x) = K(x/k, k)/k is the x-space forcing
     * Accept result vectors longer than the multipole range, as Levin always
     did; only the leading elements are written
     * Remove Levin's unreachable Simpson path and the scratch array only it
     used
     * Make Levin's property defaults the NCM_SBESSEL_INTEGRATOR_LEVIN_DEFAULT_*
     values, so property construction matches levin_new
     * Document that the range must bracket the kernel: a feature narrower than
     the Chebyshev node spacing is not seen
     * Replace the stale accuracy-floor section with the measurement against Arb
     (1e-12 to 1e-9, better than GL at 25 of 26 k)
     * Fix doc errors in the panel getters, set_max_order and new_full

 * Review sbessel integrator and GL docs

     * Fix GL truncation, remove FFT-Legendre
     * Rewrite NcmSBesselIntegrator and NcmSBesselIntegratorGL docs; scope GL as
     the reference for k b / nu <~ 1.5
     * Default the ell range to [0, 0] instead of aborting, check result lengths
     and negative multipoles in the base class
     * Remove the GL early exit, which replaced the rest of the integral by an
     invalid tail and dropped bumps far from a
     * Place the GL turning-point split in x = k chi and abort on adaptive
     quadrature failure
     * Remove NcmSBesselIntegratorFFTL
     * Add test_ncm_sbessel_integrator.c: base contract, Arb truth table
     entries, exact identity, far-bump regression
     * Reduce test_sbessel_integrator_gl.py to the Python binding tests; delete
     test_sbessel_integrator_fftl.py
     * Exempt GI (closure user_data) annotations in check_doc_style.sh
     * Regenerate stubs

 * Document spherical harmonics accuracy near the poles and test it

     * State that errors are relative to the largest values at the angle, and
     what that means for rows reached after skipped orders near the poles
     * Test the recursion near the poles against a long double reference, scaled
     by the global peak
     * Explain why the GSL comparisons stay away from the poles: GSL takes
     sin(theta) from sqrt(1 - x^2)

 * Review sf_sbessel and spherical harmonics docs

     * Fix overflow and accuracy
     * Rewrite NcmSFSBesselArray and NcmSFSphericalHarmonics docs, define the
     normalization and the skip rule
     * Compute the j_l cutoff from Debye's form: the old estimate overflowed at
     the defaults for x >= 2640
     * Compute L/x directly in the Steed recurrence, raise the threshold minimum
     to 1e-300, use parity for x < 0
     * Key ref_table on the threshold, check ell_max, abort on continued
     fraction failure
     * Convert doubles exactly in ncm_mpsf_sbessel_d, ncm_sf_sin_int and
     ncm_mpsf_0F1_d instead of by a truncated continued fraction
     * Fix the lmax = 0 segfault in NcmSFSphericalHarmonics
     * Remove NcmSFSphericalHarmonicsP, Klm_m and debug comments; skip four
     unbindable inline methods
     * Port test_sf_sbessel.py to C and delete it; add exact-argument tests;
     peak-scaled harmonics tests, fix their stale-l loop and leaks
     * Regenerate stubs

 * Review multiprecision specfunc docs; fix leaks and Si thread safety

     * Rewrite NcmBinSplit, the binsplit evaluator template, NcmMpsf0F1,
     NcmMpsfSBessel and NcmMpsfTrigInt docs
     * Add ncm_binsplit_free and free the pooled NcmBinSplit in the 0F1 and
     spherical Bessel caches
     * Move the sine integral to locked pools, add ncm_mpsf_sin_int_free_cache,
     use Si(-x) = -Si(x)
     * Remove the unused NcmMpsfSBesselRecur and fix the 0F1 macro undefs
     * Add test_ncm_mpsf_sbessel.c and sine integral symmetry and thread tests
     * Regenerate stubs

 * Review integration docs; fix error scaling, leaks and silent failures

     * Rewrite NcmIntegral1d, NcmIntegral1dPtr, NcmIntegralND and ncm_integrate
     docs
     * Use the r^2 convention in ncm_integral1d_eval_gauss_hermite1_r_p like the
     other scaled variants
     * Scale the error estimate of the Gauss-Hermite integrals by sqrt(2 pi);
     remove the dead neval counter
     * Free the user data of NcmIntegral1dPtr on finalization
     * Abort in ncm_integral_nd_eval when cubature stops at maxeval without
     meeting the tolerance
     * Abort on GSL failure in ncm_integral_locked_a_b instead of returning +inf
     * Use the cache tolerances in ncm_integral_cached_0_x
     * Rescale a copy of xgiven in the Divonne wrappers instead of the caller's
     array
     * Return the achieved relative error from ncm_integral_fixed_calibrate
     * Remove ncm_integral_fixed_integ_posdef_mult and the uncalled
     ncm_integrate_3dim, Vegas and Divonne peakfinder wrappers
     * Add C tests for each fix, including test_ncm_integrate.c
     * Regenerate stubs

 * Use KERNEL_EXACT in SijCalculator

     * Build the SijCalculator NcXcor with KERNEL_EXACT, the quadrature of
     NcXcorSSCSij
     * KERNEL_CUBATURE integrated across a Chebyshev closure panel edge where
     the slope of W(k) breaks, and its p-adaptive rule failed on cancelling
     cross pairs
     * Replace the cubature-only notes in the ssc.py docstring with the
     quadrature choice
     * Compare C and Python on KERNEL_EXACT to machine precision, and C
     KERNEL_CUBATURE against it

 * Review ode_spline and 2D spline docs; fix silent and stalling failures

     * Rewrite NcmOdeSpline, NcmSpline2d, NcmSpline2dBicubic, NcmSpline2dGsl and
     NcmSpline2dSpline docs
     * Abort on any failed CVODE call in NcmOdeSpline instead of leaving the
     spline unprepared
     * Detect steps that do not advance x (stop-hnil) and non-finite solutions
     in NcmOdeSpline; a NaN right-hand side looped forever
     * Integrate the yf mode of NcmOdeSpline toward increasing x from any start,
     end the spline at xf, apply min-subdivisions only with xf
     * Abort when auto-abstol gives a zero absolute tolerance from y_i = 0;
     clear a stale root function; remove unused fields
     * Add a default eval_vec_y for every NcmSpline2d and bound the bicubic walk
     at the last cell
     * Prepare NcmSpline2d on evaluation, accept reversed limits in the 2D
     integrals, fill the set_function matrix with NaN
     * Fix the NcmSpline2d property setters leak, the stray parent_instance and
     the NcmSpline2dGsl int_dxdy assert
     * Implement the derivatives of NcmSpline2dSpline
     * Document and skip the bicubic helper functions
     * Add ode_spline and 2D spline C tests; move the NcmPowspecSpline2d Python
     test to the powspec tests
     * Regenerate stubs

 * Review spline docs, remove NcmSplineRBF and 4POINTS, move NcmSplineFuncTest to
     tools

     * Rewrite NcmSplineBSpline, NcmSplineVec and NcmSplineFunc docs; cite
     AutoKnots in NcmSplineFunc
     * Remove NcmSplineRBF (unused, type-id never stored, integral not
     implemented)
     * Remove NCM_SPLINE_FUNCTION_4POINTS and pin the values of the remaining
     NcmSplineFuncType members
     * Move NcmSplineFuncTest out of the library into
     tools/autoknots/autoknots_stress
     * Integrate the extrapolated edge polynomial outside the knots in
     NcmSplineBSpline
     * Copy reltol and abstol in ncm_spline_copy_empty for NcmSplineBSpline
     * Fall back to eval/deriv/integ in the default NcmSpline _idx methods
     * Abort in ncm_spline_vec_set and ncm_spline_vec_set_gpa when there are no
     components
     * Document the shared cache of ncm_function_sample_set_to_spline_vec
     * Port the NcmSplineBSpline Python tests to C and add NcmSplineVec and
     NcmSplineFunc C tests
     * Fix the "cannot achieve requested precision" warning text
     * Regenerate stubs and the enum nick table

 * Improve power spectrum extrapolation and primordial range handling

     * Extrapolate CLASS matter power spectra with the Eisenstein-Hu spectrum
     times the CLASS/EH ratio continued as a power law, with the ratio's mean
     slope over the last computed decade
     * Match at the edges of the k range CLASS actually computed, not the
     requested one
     * Fix the redshift derivative outside the CLASS range, which counted the
     growth twice, and make it consistent with the extrapolated spectrum
     * Add nc_hiprim_get_lnk_range, bounded for the tabulated two-fluid
     spectrum, which aborts when no table is set
     * Restrict transfer-spectrum sampling to the primordial model's k range
     * Fix the NcHIPrimTwoFluids model name
     * Add tests for CLASS extrapolation and primordial spectrum ranges
     * Clarify power spectrum k-range documentation

 * Fix spline integration, validation, and GSL handling

     * integrate curvature norms interval by interval with fixed relative
      tolerance and convergence diagnostics
     * handle reversed limits consistently in spline integration
     * reject non-increasing knot spacing in not-a-knot splines
     * add `ncm_spline_cubic_d2_set_d2` and Python `set_d2`
     * validate GSL interpolation types before use
     * remove unused spline private fields and `PROP_ACC`
     * add C tests covering all spline fixes

 * Reorganizing FFTW interface around the new pattern.

 * Using a more robust angle computation.

 * Fix algebra utilities and numerical edge cases

     * Rewrite quaternion, polynomial-roots, and NNLS documentation.
     * Fix NNLS matrix handling, convergence checks, KKT tolerance, DGELSD
     initialization, and workspace leak.
     * Fix quaternion rotation constructors, uniform random rotations, and
     zero-axis handling.
     * Fix NcmRNG algorithm changes, seed preservation, and full-width
     seed-table keys.
     * Fix threaded function evaluation for small loops and nonpositive pool
     sizes.
     * Make NcmTimer string output robust for unknown timing values.
     * Tighten NcmVarDict getter type handling and integer-to-double conversion.
     * Fix `ncm_util_fact_size` to return the smallest 7-smooth number above the
     requested size.
     * Use the `2J+1`-weighted H-I 2p transition mean and document the
     unweighted He-I triplet convention.
     * Fix `ncm_serialize_global_peek_name` declaration and Python exposure.
     * Refactor Jacobi-Anger evaluation into accumulate and eval operations.
     * Remove `ncm_gsl_blas_types.h` and the deprecated `ncm_c_AR` API.
     * Update Python API annotations and names.
     * Add and extend C and Python regression tests, including NNLS coverage.

 * Documenting and improving Spectral.

     * Removed weighted computation.
     * Moved tests from Python to C.

 * Improving documentation and tests for LaurentSeries.

 * Improving docs and tests for Lapack, Matrix and Vector.

 * Review and improved Serialize (including bug fix).

 * Review improving ncm_cfg documentation.

 * Review constants documentation and tests.

 * Moving documentation to the correct directories.

 * Moving NcmComplex away from NcmUtil.

 * Review NcmUtil docs and added tests.

 * Review ObjArray, RNG and Timer.

     * Adds C tests for all objects.

 * Review and tests for PLN1D, FunctionCache and ISet

     * Moving LamberW work to NcmUtil and adding tests.
     * Fixing bugs in FunctionCache and PLN1D.
     * Moving tests for PLN1D to C and expanding them.
     * Adding a theory page for PLN1D.

 * Review: DTuples, func_eval and MemoryPool.

 * Moving all dev-notes out of the repo.

 * Improving metadata, removed redundant options.

     * Improving sampler metadata on catalogs.
     * Improving sampler metadata presenting in catalog analyze.
     * Removed interpolation option that is now redundant.

 * Rewrite and parallelize VKDE evaluation paths

     * Fuse the triangular solve and squared Mahalanobis-distance calculation
      into `ncm_matrix_chol_chi2_cols()`.
     * Rewrite the batched VKDE evaluator around the fused kernel.
     * Extend the batched path to sample acceptance and leave-one-out
      cross-validation.
     * Rewrite `compute_IM` to use the same fused linear algebra and
      parallelize its full pipeline.
     * Remove repeated work and allocations by precomputing kernel
      weight/normalization terms.
     * Add extensive tests for the new Cholesky-distance primitive across
      dimensions, block sizes, strided inputs, submatrices, and unused
      triangular data.

 * Cut the VKDE evaluation batch into equal tiles of at most 256 points

 * Regenerate the ncm.pyi stub for ncm_matrix_sub_row_vector

 * Optimize NcmStatsDist fitting and expand coverage

     * Warm-start auto-kernel Student-t $\nu$ fitting from the previous optimum,
     substantially reducing repeated objective evaluations and APES runtime.
     * Speed up VKDE tile centering with a new `ncm_matrix_sub_row_vector`
     operation and reuse it in the integration-matrix path.
     * Extend `NcmStatsDist` tests across split filtering, covariance recovery,
     covariance/Cholesky consistency, estimators, kernels, and cross-validation
     objectives.
     * Add complete coverage for the Python `create_stats_dist` builders and
     exercise previously uncovered CLI error and option paths.
     * Raise CI coverage of the reworked `NcmStatsDist`, KDE, and Python
     interpolation code to near-complete levels.

 * Update NcmStatsDist defaults and refactor fitting infrastructure

     * Update APES/NcmStatsDist defaults to use the benchmarked auto-kernel,
     center-shrinkage, and split-validation configuration.
     * Rework NcmStatsDist preparation and fitting to support the new defaults
     consistently across KDE and VKDE implementations.
     * Unify bandwidth and kernel-parameter optimization under the BOBYQA
     driver, including Student-t degrees-of-freedom fitting for auto-kernels.
     * Reorganize center-shrinkage and covariance handling, centralizing
     Cholesky factorization and near-PD recovery.
     * Add vectorized evaluation and reusable scratch storage to reduce repeated
     allocations during distribution evaluation and fitting.
     * Move reusable matrix and algebra operations into the generic algebra
     layer.
     * Rename internal fields, vfuncs, and split options to better reflect their
     roles, and remove obsolete helpers.
     * Record sampler and initial-sampler metadata in `NcmMSetCatalog`, persist
     it through FITS I/O, and expose it in catalog analysis.
     * Expand regression coverage for the new defaults, fitting behavior,
     covariance recovery, sampler metadata, and invalid kernel configurations.

 * Improve catalog analysis and code organization

     * Reorganize imports and improve catalog analysis.
     * Move CosmoSIS-related types to their appropriate module.
     * Expand and improve test coverage.
     * Fix minor spelling issues.

 * Refactor autocorrelation handling and Markovian catalog tracking

     * Move autocorrelation computation out of `NcmStatsVec` into a specialized
     autocorrelation object.
     * Add tests and documentation for the new autocorrelation object.
     * Add enum nicks for the new autocorrelation API.
     * Adapt existing objects to use the new autocorrelation implementation.
     * Update autocorrelation computations and related handling.
     * Add catalog support for tracking when a chain becomes Markovian,
     including the first Markovian row ID.
     * Allow walkers to choose between Markovian and non-Markovian steps based
     on the catalog state.
     * Improve symbol and label usage for consistency.

 * Update APES weighting and covariance controls

     * Remove the shrink and random-walker options and update affected
     interfaces and configuration paths.
     * Add uniform weighting support to APES and `prepare_interp`, and fix
     supplied weights being ignored.
     * Add `points-per-dim` control for local covariance estimation.
     * Improve transition-kernel point matching with a more appropriate
     row-comparison tolerance.
     * Fix MPI bookkeeping by including the missing offboard count.
     * Register the missing synthetic likelihood implementation.
     * Update the Python interface, tests, and documentation for the new and
     removed options.

 * Extend APES shrinkage support and improve ESMCMC robustness

     * Extend `NcmStatsDist` with shrinkage correction, automatic kernel
     selection, and dedicated documentation.
     * Add center shrinkage for VKDE and defensive shrinkage support, including
     the Python interface.
     * Expand `NcmStatsDist` tests, including underflowing RHS cases, and update
     Python stubs.
     * Export `NcmStatsDist` configuration through APES and expose the new APES
     options in Python.
     * Add CLI support for generating synthetic experiments and extend APES
     rebuild/testing coverage.
     * Fix ESMCMC stale-state cleanup and cumulative OpenMP state, and add
     explicit seed configuration.
     * Add MPI versus OpenMP ESMCMC testing.
     * Prepare all dataset components before evaluation and avoid unnecessary
     Boltzmann re-preparation.
     * Improve handling of massive neutrino configurations and match Planck
     neutrino conventions.
     * Fix singular solves in `NcmNNLS`.
     * Update project files and development notes for the new functionality.

 * Add APES centre shrinkage and improve sampler robustness

     * Add covariance-preserving centre shrinkage to `NcmStatsDist` KDE/VKDE
     mixtures and expose it through the APES walker, Python helpers, and catalog
     calibration.
     * Account for kernel covariance in centre shrinkage, including Student-t
     variance factors, and reject kernels without finite covariance.
     * Preserve APES `local-frac`, covariance type, and fixed covariance across
     estimator rebuilds, and fix `METHOD_KDE` to construct the KDE estimator.
     * Defer centre-shrink compatibility checks until the kernel is configured,
     fixing valid Gaussian and Student-t APES configurations.
     * Make the Python APES wrapper apply the supplied initial sample before
     `start_run()`.
     * Harden NNLS interpolation against empty and rank-deficient systems, with
     DGELSD and uniform-weight fallbacks.
     * Fix SNIa absolute-magnitude dataset counting when objects are
     reconfigured or duplicated for threaded runs.
     * Add centre-shrinkage benchmark and validation notes with reproducible
     scripts and update the related documentation.

 * docs: add BAO SDSS DR16 fitting example

 * More coverage tests.

 * Updating UltraLevin project. Improving Xcor CLI.

 * xcor cls command (#374)

     Improve Xcor production C_ell workflows and visibility modeling

     * Add production-oriented `xcor cls` with block processing, multipole
     sampling, pair selection, shared solver state, and reusable CLI/kernel
     setup.
     * Make exact Chebyshev integration the CLI default, rename `fixed` to
     `exact`, and improve generated kernel documentation and metadata handling.
     * Rename `scaled-abstol` to `peak-epsilon` across Xcor/SSC APIs, docs,
     tests, defaults, constants, and stubs, clarifying its role as a
     peak-relative closure-fit floor.
     * Improve RSD evaluation to use P(k) directly, with 2D spline derivatives
     and low-order finite-difference fallback.
     * Add full and minimal visibility-function support for CMB lensing and ISW,
     with selectable thin-screen or visibility treatments.
     * Add animation generation support and expand UltraLevin documentation.
     * Expand tests for CLI solving, quadratures, pairings, kernel construction,
     RSD derivatives, visibility modeling, enum truth tables, and documentation
     completeness.
     * Update generated stubs, formatting, CLI help, and general prose.
 * Documenting ultra levin (#373)

     Improve UltraLevin ODE constraints, notation, and documentation

     * Add free, pinned, and pin-at-peak ODE constraints and propagate them
     through downstream solvers and operators.
     * Certify free and pinned solutions against Arb across analytic windows,
     single and batched solves, including guard and regression tests.
     * Expand UltraLevin documentation and figures to compare Dirichlet, pinned,
     and free constraints, including boundary-term checks, tail behavior, and
     timings.
     * Improve spectral adaptivity with adaptive try, configurable order caps,
     and better wisdom handling.
     * Standardize symbols and notation across the library: use `chi`
     consistently in Xcor while retaining `D_c` in distance APIs, documenting
     the Xcor convention and its connection to the literature.
     * Improve documentation tone and prose across the library, removing
     non-ASCII symbols, indirect narration, and inconsistent terminology.
     * Improve symbol API linking and reorganize UltraLevin files and tests.
     * Expand certification to more multipoles and modes and improve
     tail-by-tail and timing documentation.
     * Fix warnings and linter issues and ensure numcosmo-site depends on
     current typelibs.
 * docs: add H(z) fitting example

 * CI: update conda lock files

 * Skip the last Chebyshev doubling when the level below it predicts failure

     * Add _ncm_spectral_batch_cannot_converge(), which extrapolates a level's
      coefficient envelope decay to the mass the next level would add
     * Skip level k_cap for fatal FALSE callers of
      ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap() when that mass
     is
      predicted above 10x the acceptance tolerance, so the caller splits instead
     * Document why the saving exists: a bisecting caller discards half its
     sampling
      because a child Chebyshev-Lobatto grid is not a subset of its parent's
     * Measured 1.26x at reltol 1e-4 and 1.19x at 1e-6 on a 15-spectrum
      KERNEL_EXACT solve, with all 600 C_ell values bit-identical and the panel
      count unchanged
     * Read the midpoint of chi through a local copy in window_u(), silencing a
     GCC
      16 -Wstringop-overread on every later acb_sub() in the function
     * Clear the Par before flint_cleanup() in nc_xcor_kernel_analytic_arb,
     which
      also puts the previously unused par_clear() to work
     * Drop the unused pass and negligible variables in nc_xcor_kquad_arb

 * Use the attached recombination history for the decoupling redshift (#361)

     Keep CMB distance-prior redshifts consistent and fix CBE recombination
     redshifts

     * Release the `NcRecomb` reference owned by `NcDistance`, fixing a leak.
     * Keep `nc_distance_decoupling_redshift()` in the Hu–Sugiyama convention
     required by CMB distance priors.
     * Derive CBE characteristic redshifts consistently from its optical-depth
     splines, including the tau = 1 and cutoff redshifts.
     * Clarify the visibility-function convention and add cross-backend and
     distance-routing regression tests.
     * Add `nc_recomb_cbe.h` to the umbrella header and clean up recombination
     tests and documentation.
 * Using chebyshev as default (#369)

     Set Chebyshev as the default k-sampling and integration method

     * Make adaptive component boundaries panel edges in the Chebyshev closure,
     sample only components supported on each panel, and validate accepted
     panels against domain-expansion samples to prevent unresolved narrow
     features.
     * Truncate components with touching supports jointly, preserving
     cancellations between component tails and substantially improving closure
     efficiency and accuracy.
     * Correct the certified multi-window Arb reference by zeroing gaps between
     disjoint bump groups and regenerate the affected truth-table entries.
     * Extend the certified lensing reference integration to k = 10/Mpc to
     capture the slow high-k tail from the window's hard edge.
     * Add regression tests covering component boundaries, joint truncation, and
     narrow-feature detection.
     * Add a theory page documenting the full spherical-Bessel projection
     pipeline, including Levin integration, ultraspherical discretisation,
     adaptive fitting, closures, and outer k integration.
     * Review and improve the spherical-Bessel ODE solver, spectral,
     spherical-harmonics, and related theory/API documentation.
 * Drive the CCL bridge through NcXcorKernelTable and NcXcorSolver (#368)

     Support component tables and persistent CCL solves

     * Add multi-component tabulated kernels with density, shear, convergence
     and RSD terms.
     * Route the CCL bridge through `NcXcorKernelTable` and `NcXcorSolver` using
     the exact integration path.
     * Reuse prepared kernels and block integrators across unchanged solves and
     replans.
     * Add in-place table updates and persistent `TracerClSolver` for efficient
     repeated CCL evaluations.
     * Fix spin-2 origin handling and automatically extend the distance range
     for CCL tracers.
     * Extend tests for component tables, CCL mapping, solver reuse, updates and
     refusal paths.
 * Close the RSD follow-ups from the post-merge review (#367)

     Fix RSD follow-ups and CCL bridge Bessel transforms

     * Correct the growth-function derivative docs and tighten the
     spherical-Bessel derivative contracts and Levin order gating.
     * Document the incompatibility of RSD kernels with redshift-space Limber
     methods and expand RSD, serialization, solver, and integrator regression
     coverage.
     * Integrate CCL’s spin-2 `j_l(kχ)/(kχ)^2` weight exactly instead of
     replacing it with its Limber approximation.
     * Handle the spin-2 `1/χ²` factor inside the Levin solve, with linear
     continuation below CCL’s first kernel sample.
     * Replace the fixed outer `k` window with adaptive growth based on bounds
     for the omitted integrated tails.
     * Factor block-transform setup from evaluation so the outer grid can grow
     incrementally without rebuilding the transform.
     * Reduce lensing-bridge errors to `2.4e-5–6.0e-5` at `ell=2–10` while
     cutting the bridge test runtime from 324 s to 37 s.
 * Keep the FIXED_NODES redshift marginal non-negative by construction (#301)

     Fix cluster-WL fixed-node integration stability and warning handling

     * Report achieved calibration error through an optional `relerr_out` and
     suppress per-call warnings when the caller handles reporting.
     * Aggregate failed node calibrations into one warning per `prepare()`,
     including the failure count and worst relative error. 
     * Make FIXED_NODES galaxy probabilities structurally positive using a
     self-normalised panel combination, avoiding spurious `LOW_PROB` walls.
     * Add `nc_data_cluster_wl_factor_get_low_prob_count()` to expose
     probability fallbacks.
     * Improve fixed-grid accuracy with the self-normalised control-variate
     form, including on under-resolved grids. 
     * Keep auto-nodes opt-in with the existing `1e-4` tolerance; document its
     discrete-grid discontinuities and align C and Python defaults. 
     * Update regression, parity, warning, and grid-rebuild tests for the new
     integration behavior and defaults.
 * Close the coverage gaps of the RSD additions

     * Drop the derivative weights from the Levin direct cubature path: the
      path is disabled by the multipole threshold, so the weights were
      untestable; it now fails loudly if it is ever re-enabled with a
      derivative request, matching the GL/FFTL decision.
     * Test the dorsd guard of the redshift-space Limber tier with a
      g_test_trap subprocess.
     * Test the component list of a dorsd kernel (orders 0 and 2), the
      dorsd and bessel-deriv property reads, and component disposal.
     * Test the order-1 Limber peak formula against the exact solve at
      ell = 200 by reconfiguring the density component (measured -8.4e-3,
      the Limber error of this suppressed integral).

 * Removing gtk-doc comment from header.

 * Silence maybe-uninitialized and duplicate-doc warnings in the Levin integrator

     * Zero-initialize the IBP boundary struct in both accumulation paths;
      GCC cannot see that the reads share the deriv > 0 guard with the fill.
     * Drop the header doc block duplicating NcmSBesselIntegratorLevin:,
      which the GIR scanner warned about on every build; the .c block stays.

 * Expose dorsd in numcosmo xcor kernel view

     * Add dorsd (default False) to KernelNumberCountsConfig and pass it to
      NcXcorKernelGal, so the CLI can visualize kernels and C_ell with the
      redshift-space distortion component in both tiers.

 * Add redshift-space distortions to NcXcorKernelGal

     * Add a bessel-deriv property (0-2) to NcXcorKernelComponent: the order
      of the spherical Bessel derivative weighting the component's radial
      integral.
     * Exact tier: the non-Limber closure integrates each component with
      ncm_sbessel_integrator_integrate_deriv() at its own order.
     * Kernel-space Limber tier: replace the three duplicated peak-formula
      sites with one deriv-aware helper; orders 1 and 2 peak-approximate the
      j_l and j_{l+1} pieces of the downward recurrences at nu/k and
      (nu+1)/k, matching CCL's Limber treatment of der_bessel = 2.
     * Add dorsd to NcXcorKernelGal, following CCL's NumberCountsTracer
      convention: a component with kernel -f(z) dn/dz E(z) sqrt(P), growth
      rate from NcGrowthFunc, bessel-deriv = 2. RSD disables the constant-
      bias fast-update shortcut and errors out loudly in the legacy
      redshift-space Limber tier.
     * Extend the CCL bridge (numcosmo_py.ccl.two_point) to der_bessel 1 and
      2 through integrate_deriv.
     * Validation recorded in the tests: NumCosmo growth rate matches CCL to
      2e-9; the kernel-space Limber RSD matches CCL's Limber to ~4e-6 at
      ell = 100-300; at ell = 2-10 the exact tier and the CCL-kernel bridge
      agree to 1e-6 while pyccl's FKEM non-Limber is 38% (ell = 2) and 11%
      (ell = 10) low on a Gaussian bin, so the high-ell Limber comparison is
      the CCL reference used.

 * Add derivative-weighted integrals to NcmSBesselIntegrator

     * Add ncm_sbessel_integrator_integrate_deriv() computing
      int K(x,k) j_l^(d)(kx) dx for d = 0, 1, 2, with a new integrate_deriv
      virtual (default errors out; Levin implements it).
     * Implement the Levin path by double integration by parts: the volume
      solve uses y F'(y) or y F''(y) as forcing, read off the panel fit in
      the ultraspherical basis, plus [F j_l] or [F j_l' - F' j_l] boundary
      terms per panel. No Bessel recurrences in l and no pointwise
      evaluation of F''.
     * Add NcmSpectral maps chebT_deriv_to_gegenbauer_alpha2 (two-band),
      chebT_deriv2_to_gegenbauer_alpha2 (diagonal) and
      gegenbauer_alpha2_xmul (affine multiplication in C^(2)).
     * Add ncm_sf_sbessel_jl_deriv_from_array() and
      ncm_sf_sbessel_xjl_deriv_from_array(); the latter replaces the
      private yj_deriv helper of the Levin integrator.
     * Pad the Levin RHS to the solver's three-coefficient minimum, needed
      when a derivative forcing collapses to fewer coefficients.
     * Test against closed forms (constant K reduces to endpoint values of
      j_l and j_l'), scipy quadrature for Gaussian K across evanescent,
      turning-point and oscillatory regimes, batched-vs-single consistency,
      and the l-recurrence identity as an independent cross-check.

 * Fix the kquad Arb generator applying side a's kdep to side b (#364)

     * Fix the kquad Arb generator applying side a's kdep to side b
     * nc_xcor_kquad_arb.c: stop initializing side b as a copy of side a on the
     first --b: argument; the kdep_on latch has no off switch and leaked across
     * Start side b from its par_init state plus the shared n_sigma/n_scale
     defaults
     * Regenerate X9, which was certified with kdep suppression on both windows
     (0.06%/0.4%/5.3% low at ell 2/10/50); the library now agrees with the
     corrected values to better than 7e-9
     * Regenerate X5 as a null test; it remains bit-identical
     * make_xcor_kquad_truth_table.py: derive the ells header from the entries
     so partial regeneration cannot relabel the table
     * test_k_integral.py: tighten the Chebyshev certified gate 1e-1 -> 3e-2 and
     update the measured constants (closer in 36/43, medians 4.0e-6 vs 1.8e-5)
     * test_k_integral.py: document that a tolerance- and closure-independent
     deviation indicates a wrong reference, not a bad fit
     * Document what gates flipping the closure default to Chebyshev
     * NcXcor:closure-type doc: record the measured Chebyshev case (closer to
     certified truth in 36/43, medians 4.5x smaller at every rung, no
     catastrophic regime, +6% median cost)
     * Record the first flip blocker: KERNEL_EXACT aborts on mixed spline x
     panel pairs produced by differing per-kernel l-limber values, because
     Limber blocks always build splines
     * Record the second flip blocker: the spectral exact path returns no error
     estimate, so flipping the default would silently withdraw error reporting
     * Keep the spline default and state the flip conditions where the property
     is defined
 * Rewrite comments and docs in factual language without dropping content (#363)

     * Restore docs/theory/wl_shape_factor_history.md and the 16 pointers to it
     * Restore Laplace's NaN return on a non-negative-definite Hessian
     * Restore Quad's clamp that prevents a Cuba segfault on a non-finite sample
     * Restore NcXcorKernelRadial's Levin abort on an unclamped chi
     * Restore SeriesLensed's underflow floor and its reason
     * Restore cubacores(0, 0) as what makes Quad safe under OpenMP
     * Restore the NcmMemoryPool checkout contract in SeriesLensed
     * Restore the blocking behaviour of ncm_memory_pool_empty and _free
     * Restore FixedQuad's lower limit on population width, near sigma_pop =
     0.05
     * Restore the alpha >= 2 requirement for SeriesLensed with a Beta
     population
     * Restore the Levin panel-placement accuracy limit and the GL alternative
     * Restore the eval_kernel rule that edges belong in get_limits
     * Restore NC_XCOR_SSC_SIJ_DEFAULT_SCALED_ABSTOL's pairing with
     numcosmo_py/ssc.py
     * Restore why the data cache has no restore-keys and why grep needs -I
     * Restore why FFTW wisdom is excluded from the cached paths
     * Restore the per-test thread policy that avoids the mixed-OpenBLAS
     deadlock
     * Restore the deny-list convention that stops a new suite being skipped
     unnoticed
     * Restore the theory-page links from C and the [[numcosmo|Sym]] links from
     .qmd
     * Restore the vkde.png flowchart reference, which meson still installs
     * Restore the NcmStatsVec online update formulas and usage example
     * Move the vp_err calibration tables to dev-notes/xcor_exact_quadrature.md
     * Point nc_xcor_compute_full() at that section
     * Fix ncm_stats_vec_heidel_diag(): bindex is the smallest qualifying index
     * Fix TESTING.md: accepting a coverage loss differs from losing it
     unrecorded
     * Add .github/scripts/check_doc_style.sh and a CI job that runs it
     * Document the wording rule and placement convention in CONTRIBUTING.md
 * Update TESTING.md for the lane-based CI layout

     * Replace the derived-slice sharding section with the lane-based one
     * Document the tier/lane map and which lanes are instrumented (new 4.1)
     * Add "split by claim, not by check" guidance for mechanics vs statistical
     claims
     * Fix the capability-gate section for the get_closest_marker conftest gate
     * Sync the local coverage recipe: instrumented lane list and tests/
     exclusions

 * Silence the _FORTIFY_SOURCE warning in the coverage build

     * Undefine _FORTIFY_SOURCE alongside the -O0 that conda's -O2 forced us to
     append

 * Name the two test sets in CI instead of subtracting and re-adding

     * Stop excluding py-omp from the fast subset only to run it again on the
     same legs
     * Collect the validation-tier suites in VALIDATION_TEST_ARGS and run them
     in one step

 * Run the C tests as one coverage shard instead of three

     * Replace the three --slice legs with a single --suite c leg
     * Drop the slice branch from the coverage run block

 * Run the C acceptance and statistical tiers optimized, not instrumented

     * Move both suites out of the coverage job into the optimized
     build-miniforge job
     * Leave the coverage job with the unit suite alone

 * Revive dead test code and keep the harness out of the coverage figure

     * Register test_ncm_model_test_finite, which was defined and never wired up
     * Point /nc/halo_bias/integrand at its own function instead of set_get
     * Give that test the mass-function prepare and the prim slot it needs to
     run at all
     * Add a primordial power-law test, covering the tensor spectrum CLASS reads
     * Drop tests/c and tests/python from the published coverage, keeping
     installed helpers

 * Close the reachable xcor gaps in nc_xcor and the kernel component

     * Make the disjoint-pair case actually disjoint, and assert the spectrum is
     zero
     * Read back every NcXcor property
     * Read the component knobs through the property interface, and ref/clear it

 * Exercise both xcor closures in the kernel-space quadrature tests

     * Set the closure type explicitly, the default spline having hidden the
     spectral path
     * Assert the exact method reports no error estimate on a Chebyshev pair, as
     it documents

 * Test the NcmMSetTransKernCat CHOOSE sampler on a synthetic catalog

     * Build the catalog by hand, with rows laid inside the mset bounds and
     distinct m2lnL
     * Check every drawn point is a catalog row and within bounds
     * Check the percentile cut restricts the draw, asking for no more points
     than it allows

 * Test NcmMSetCatalog on a synthetic multi-chain catalog

     * Insert rows directly into a four-chain file-backed catalog, sampling
     nothing
     * Check the shrink factor, which returns 1 outright for a single chain
     * Check column lookup by name against the listing, and the summary
     accessors

 * Test the NcmStatsVec convergence diagnostics on curated chains

     * Add white-noise, AR(1), burn-in and late-break chains from a fixed seed
     * Check heidel_diag, max_ess_time and visual_heidel_diag against them
     * Check that a chain with a burn-in is not reported stationary from the
     start
     * Cover the saved-rows copy, robust diagonal covariance, quantiles and
     correlation

 * Correct the l_limber sense in the xcor integrand and kquad tests

     * Use -1 for never-Limber and 0 for always-Limber, which had been swapped
     * Set the kernel-space methods' kernels to never-Limber, as their names
     claim

 * Gate capability tests on their marker, not on path keywords

     * Check get_closest_marker instead of item.keywords, which also holds
     directory names
     * Collect the marker/option pairs in one table

 * Move the Planck 2018 CLI tests out of the app gate

     * Split test_generate_planck and test_generate_planck_test into their own
     file
     * Place it outside tests/python/numcosmo_py/app, whose name alone gates on
     --run-app

 * Split NcmStatsDist tests into mechanics and distribution recovery

     * Move the shared checks to a common compilation unit with a mode selected
     by main()
     * Build the estimators from a small sample in the instrumented lane
     * Keep the divergence, covariance and sampling-frequency claims in the
     statistical tier

 * Split NcmFitESMCMC tests into mechanics and covariance recovery

     * Move the shared checks to a common compilation unit with a mode selected
     by main()
     * Run the mechanics on a short fixed chain in the instrumented lane
     * Keep the covariance-recovery assertions and their retry loops in the
     statistical tier
     * Copy the data covariance in run() before cov2cor rewrites the peeked
     matrix in place

 * Add C unit tests for the xcor subsystem

     * Cover the seven analytic windows: construction, properties, supports,
     serialization
     * Cover the concrete kernels through the NcXcorKernel surface, with cheap
     inputs
     * Cover the integrand accessors for both closures, with the adaptive
     apparatus capped
     * Cover the Limber methods, the kernel-space block quadratures and NcXcor's
     facade
     * Cover NcXcorSolver registration, the block planner and batched-vs-direct
     agreement
     * Cover NcXcorSSCSij's configuration, including area and mask being
     mutually exclusive

 * Handle k ranges degenerate to within rounding in NcXcorKernelComponent

     * Skip a sampling-grid k range narrower than a few ULP instead of testing
     exact order
     * Warn only when the range is empty by more than rounding, not at the grid
     endpoint
     * Return the endpoint rather than hand gsl_min_fminimizer_set an
     unbracketable interval

 * Give two g_message lines the '# ...\n' shape

     * Prefix and terminate the ncm_fftlog_calibrate_size_gsl size-cap message
     * Prefix and terminate the ncm_function_sample_set_adaptive_midpoint
     max-iter message

 * Split the CLASS-heavy NcCBE checks into their own executable

     * Measured per check, across the eight cosmologies each runs: Cls costs
     12.3 s
      per model and calc_ps 4.4 s, while sanity, serialize, precision, thermodyn
      and compare_bg are 0.4 s or less. Those two are 91% of this file's runtime
      and reach about 2% of the lines nothing else reaches
     * Move them to their own acceptance-tier executable and leave the cheap
     checks
      in the instrumented lane, which drops from one 42 s binary to 3.98 s plus
     a
      separate 38.48 s one -- two units --slice can place in different shards,
      where one executable could not be divided at all
     * Split by check rather than by sharing a prepared state across them.
     Sharing
      would save about 30%, not the several-fold the repeated setup suggests,
     since
      Cls genuinely computes deeper than compare_bg rather than redoing it, and
     it
      would couple subtests: calc_ps mutates the object with
      nc_cbe_set_calc_transfer()
     * Same structure as the fit split: models and checks in
     test_nc_cbe_common.c,
      one small source per group naming the checks it runs
     * Registers exactly the same 57 subtest paths as before, none added or lost

 * Give each NcmFit algorithm its own test executable

     * test_ncm_fit.c was compiled five times with different -DTEST_FIT_* to
     make one
      binary per optimizer group, and inside it TESTS_NCM_ADD token-pasted 36
      g_test_add calls per algorithm. The preprocessor was doing both jobs for
     no
      reason: g_test_add already passes its tdata to setup, test and teardown,
     so a
      descriptor read from fixture data serves every algorithm and the generated
      per-algorithm symbols are unnecessary
     * Move the checks, unchanged, into test_ncm_fit_common.c, which every
     executable
      links, and give each algorithm a source that names it and calls
      test_ncm_fit_main. The registrations become a table rather than a macro,
     so
      the set of checks is readable in one place and greppable by name
     * Fifteen executables instead of five. The longest single fit executable
     falls
      from 61.5 s to 22.2 s, which is what actually bounds the C lane, since
     neither
      --num-processes nor --slice can divide one executable
     * Registers exactly the same 525 subtest paths as before, none added or
     lost
     * Fix a leak this uncovered: test_ncm_fit_equality_constraints and
      test_ncm_fit_inequality_constraints return early when the algorithm is not
      NLOpt:slsqp, which is 13 of the 15, and both returns skipped
      ncm_mset_func_free. They are not dead tests -- each exercises
      ncm_fit_add_*_constraint before skipping the solve

 * Drop the last references to duration-based slicing

     * CONTRIBUTING.md pointed the repository's cache quota at the test-duration
      caches, which no longer exist; the data-file cache is what is left
     * TESTING.md said wall time is balanced by --slice and then that --slice
      balances by count. Say plainly that it balances by count, why weighting by
      duration was removed, and that the remedy for a fat test is to split it

 * Slice the coverage C tests by count and drop the duration machinery

     * The three C legs each computed their own duration-balanced slice, which
     is a
      partition only if all three weigh the tests identically. They restored the
      duration bundle through a restore-keys prefix, which resolves to whichever
      entry is newest when that leg starts, and the legs do not start together,
     so
      a bundle saved by a concurrently finishing run gave later legs different
      weights and three incompatible partitions
     * Seen on run 33457691603: 84 selections covering 67 unique tests, so 17
     ran
      twice and 16 not at all. codecov read that as -3.72% with whole test files
     at
      zero hits, which looked like a regression in the PR under review
     * Use meson's own --slice instead. Measured over the current suite that
     costs
      115s -- a 518s worst slice against 403s -- on legs that take 9-11 min
     while
      this job's critical path is py-xcor at 30-47 min. When one test grows
     enough
      for that to matter the answer is to split it, which is worth doing anyway
     * Removes the duration cache, its version knob, the per-leg extract and
     upload
      steps, and the merge-test-durations job. test_slicer.py loses plan and
      extract, leaving only the timing report, so it becomes test_summary.py
     * Net 98 lines removed from the workflow, and the failure mode is gone
     rather
      than guarded: there is no longer a per-leg input that can disagree

 * Initialize MPI only under a parallel launcher

     * ncm_cfg_init called MPI_Init unconditionally, so every process paid about
     a
      second of OpenMPI/UCX device probing before doing any work: every CLI
      invocation, every one of the 91 test binaries, every pytest worker, and
     every
      CLI subprocess the app tests spawn. Measured 1086 ms -> 53 ms on a test
     binary
      that runs no tests
     * Outside a launcher the call buys nothing. The world is a single rank, so
      _mpi_ctrl keeps the size/rank/nslaves defaults set just above, and the job
      dispatch in ncm_mpi_job.c is reached only when nslaves > 0
      (ncm_fit_esmcmc.c:978), so it never runs
     * Detect the launcher from its own environment rather than initializing and
      backing out: a rank other than the master enters the slave main loop
     during
      init and never returns, so deferring that decision would run the caller's
      whole program once per rank
     * NUMCOSMO_MPI_INIT overrides the detection in both directions, for a
     launcher
      that sets none of the known variables (=1) or to force the serial path
     (=0)
     * Full suite passes 91/91 with MPI uninitialized on every lane but the two
     mpi
      ones, which still run under mpiexec -n 2 and take the unchanged eager path

 * Fail a stuck test instead of hanging its lane

     * Nothing killed a stuck test: pytest-timeout was never configured, and
      conftest.py only re-arms a faulthandler traceback dump, which reports but
     does
      not terminate. A worker that died mid-test therefore left the xdist
     controller
      waiting and hung the whole lane, twice observed at 20+ min before the job
     was
      killed by hand
     * Add --timeout=900 --timeout-method=thread. 900 s is sized from the
     measured
      maxima under -O0 + gcov, the slowest lane by roughly 4x: acceptance 273 s,
      xcor 113 s, omp 75 s, default 73 s, everything else under 10 s
     * thread rather than the default signal: SIGALRM is delivered only between
      Python bytecodes, so it cannot interrupt a test blocked inside a long C
     call,
      which is where this library spends nearly all of its time
     * Declare pytest-timeout in pyproject, environment.yml and the meson module
     list
      so a missing plugin fails at configure time rather than as an unrecognised
      pytest option mid-run, and regenerate the conda locks in the same commit
     so
      the lock check stays green

 * Run the acceptance tier optimized instead of instrumented

     * py-acceptance was excluded from every optimized job by FAST_TEST_ARGS and
     ran
      solely as a coverage leg, which is both the slowest way to run it and the
      least useful one: 48 min on CI against 84 s on an optimized build locally,
      and it was the critical path of the whole coverage job
     * Per-line it is almost entirely redundant there -- it reaches 31021 lines,
     of
      which only 170 are reached by nothing else, measured by running every lane
      alone against a wiped .gcda set
     * Run it as its own step on the optimized Miniforge build and drop the
     coverage
      leg. Gated to ubuntu/openmpi, the one combination it already ran on since
     the
      coverage matrix is ubuntu-only, so this changes where it runs, not what it
      covers, and adds nothing to the macOS job that is the longest of that
     matrix
     * Expected coverage change is -170 lines, 0.18%, inside codecov's 1%
     threshold

 * Revert "Distribute the pytest lanes with worksteal instead of pinning files"

     * This reverts cdb1c638. On CI the xcor coverage lane went from 30m46s to
     47m
     * The pins did more than bound memory: keeping test_k_integral.py on one
     worker
      also built its Frozen cache once, where worksteal spreads that file's
     cases
      and every worker rebuilds what it receives
     * CI runs two workers, not four -- pytest-xdist's `-n auto` uses the
     physical
      core count and the runners have two physical cores -- so there is no spare
      parallelism to hide that recomputation, and the pinned tail is the cheaper
     of
      the two costs
     * Measured on the wrong configuration twice before this: the Optimized
     build at
      twelve workers, then Coverage at four. Both reverse the ordering CI sees
      (worksteal 729 s vs loadgroup 1118 s at four workers), because the fewer
     the
      workers, the less parallelism there is to absorb the rebuild

 * Block pytest-randomly in the test lanes

     * It is not a dependency and is absent from the CI conda locks, but it
     activates
      itself wherever a developer has it installed, shuffling test order and
      reseeding random and numpy.random before every test
     * Local runs were therefore not reproducible while CI runs were, contrary
     to the
      reproducibility TESTING.md describes
     * -p no: is a no-op when the named plugin is absent, so this adds no
     dependency

 * Distribute the pytest lanes with worksteal instead of pinning files

     * Replace --dist loadgroup with --dist worksteal, so a worker that runs out
     of
      tests takes work from one that has not, rather than idling
     * Drop the five xdist_group pins in tests/python/nc/xcor
     * The pins and the frozen-cache LRU bound landed together in 81adac9f;
     measured
      separately, the bound alone holds the memory and the pins cost 5.1x wall
     time
      on the xcor lane (278 s against 54 s), with 2186 tests passing either way
     * At the four workers CI runs, worksteal peaks at 4.5 GB

 * Relax an unreachable tolerance in the xcor view CLI test

     * test_view_kernel_integrator_reltol_reaches_the_computation drove the CLI
     at
      --integrator-cheb-reltol 1e-12, below the integrator's accuracy floor, so
     the
      Levin RHS Chebyshev fit could not converge and reached the fatal max-order
      error; that aborts the process, which kills the xdist worker running it
     and
      leaves the controller waiting, hanging the whole app lane
     * Use 1e-8: the loose and tight runs still differ by 4.6e-5, 46x the
     threshold
      the test asserts, so it checks exactly what it did before
     * 3482d85c removed the absolute floor that previously let this fit
     terminate;
      the accompanying search was by symbol name, which cannot find a caller
     that
      names no abstol symbol

 * Remove the last references to the deleted sbessel abstol

     * test_xcor_window_truth_table.py still called set_abstol, removed in
      3482d85c, so the xcor lane failed on AttributeError before integrating
     * Drop the call and the INTEG_ABSTOL_FRAC constant; the comment claiming
     the
      deep-tail entries abort at max-order without a floor is stale, all 15
     tests
      pass on the relative criterion alone
     * Reword the scaled-abstol paragraph in nc_xcor_kernel.c that described the
      removed constant; the precision limit far below the peak is cancellation

 * Remove the caller abstol from the spherical-Bessel integrator

     * The Levin RHS Chebyshev fit now uses the relative criterion alone; a fit
      that cannot converge relatively is noise and stops at the existing fatal
      max-order error instead of being accepted under a floor
     * Replace the panel abstol computation with an explicit null-panel skip:
     when
      the Bessel weight is exactly zero at both endpoints for every ell in the
      batch, the contribution is identically zero and the panel is skipped
     * Remove ncm_sbessel_integrator_set_abstol/get_abstol; the only caller was
      nc_xcor_kernel.c arming 1e-16 times a running maximum that starts at zero,
      so the floor was smaller than intended until the maximum was measured
     * Remove the integ_max running maximum and its arming from
     nc_xcor_kernel.c;
      the k-space fit tolerances (scaled-abstol, sample-set peak) are unchanged
     * The edge-panel coefficient trim limit keeps only its relative part
     * Remove TestPanelAbstolEvanescent: it fed a synthetic integrand no
     Chebyshev
      order resolves and asserted the floor turned the abort into a finite
     return;
      that behavior is withdrawn on purpose, such input now fails loudly
     * Measured on the full N5K workload (103 ells x 120 spectra, serial and 12
      threads): sn unchanged (2.441384), spectra differ by at most 1.5e-12 of
      block peak, only at ell < 200, in the more-refined direction

 * Collect coverage from unoptimized code

     * conda's compiler activation injects -O2 into CFLAGS and CPPFLAGS, which
      outrank -Dbuildtype=debug, so coverage has been measured from optimized
      builds while meson reported "Optimization level : 0"
     * Inlining detaches a function's line counts from its own record, which
      lcov 2.x rejects outright, and drops inlined functions from the report
     * Applied in the environment rather than per target so it reaches the test
      executables, not just the library

 * Reject undersized two-fluids initial-condition vectors

     * The state has NC_HIPERT_ITWO_FLUIDS_VARS_LEN (8) components, but
      test_evolve_array passed a 6-element vector, so get_init_cond_zetaS wrote
      two doubles past its end and set_init_cond read them back
     * ncm_vector_set does not range check, so this was silent: at -O2 the bytes
      past the end were nonzero and CVODE ran, at -O0 they were zero and its
      error-weight vector became illegal, failing with CV_ILL_INPUT
     * Guard the three public functions taking init_cond and size the test's
      vector like the examples already do

 * Stop instrumenting the bundled external libraries

     * Set b_coverage=false on the 13 static libraries under numcosmo/external
     * Their coverage was captured and then discarded by the --remove step, and
      the Fortran sources in plc broke the lcov 2.x consistency check
     * libfyaml already overrode c_std, so merge both into a single list rather
      than passing the keyword twice, which meson resolves by dropping the first

 * docs: use the house array and DataFrame naming in the SNIa+BAO example

     Rename the plain arrays to the `_a` suffix (centre_a, theta_a, 
     unit_circle_a, points_a) and the best-fit frame to best_fit_pd, matching 
     the naming used across the other examples. The two geom_point calls wrap 
     because the longer name pushes them past 88 columns.

     No numerical change: the three area ratios (1.759, 1.782, 1.814) and the 
     expected/observed sigma_w (0.09969, 0.119) are identical after 
     re-execution.

 * Make NcmSplineBSpline evaluation thread-safe

     * Evaluation called gsl_bspline_calc on the instance's workspace, which
     every
      call uses as scratch; concurrent evaluation of one shared spline corrupted
      both results, and the xcor solver evaluates shared table kernels from its
      OpenMP block loop
     * Evaluate with de Boor's recursion on stack scratch instead, reading only
      state frozen at preparation; values agree with the GSL path to 7.8e-16
      absolute, derivatives and integrals are bit-identical
     * Serialize the derivative, integral and lazy-name paths, which stay on the
      GSL workspace, on a per-instance mutex
     * Add ncm_spline_bspline_threads: concurrent evaluation of one shared
     spline
      must match serial evaluation bit for bit; fails on the unfixed code

 * CI: update conda lock files

 * Add a halo mass function example (#343)

     Add and refine the halo mass function documentation example

     * Add a halo mass function example and wire it into the examples index.
     * Correct `dn/dlnM` labels and use `log10(M/Msun)` for the mass axis.
     * Set primordial parameters explicitly and report `sigma8` for reproducible
     normalization.
     * Explain the FFTLog filter range setup and the physical origin of the
     `dN/dz` peak.
     * Compute and print the counts peak redshift instead of duplicating it in
     prose.
     * Clean up parameter logging, comments, and unused reionization
     configuration.
     * Construct the cosmology with its primordial and reionization submodels,
     following the construction-fixed submodel API.
     * Align the example with the existing documentation style and current
     NumCosmo API.
 * Apply the scale dependence on the radial Limber branch

     * NcXcorKernelRadial's non-Limber component multiplies by its kdep factor
     and
      the Limber branch silently dropped it, so a scale-dependent kernel gave a
      wrong Limber C_ell with the right shape
     * Limber has k = (l + 1/2) / chi in hand, so there is nothing there it
     cannot
      evaluate
     * Assert the Limber and non-Limber branches agree at high ell with a scale
      dependence attached; without the fix the two differ by 0.93%

 * Stop the radial Limber branch from counting growth twice

     * NcXcorKernelRadial carries the growth in W and pairs it with P(k, 0), but
      nc_xcor_limber_z.c multiplies the two kernels by P(k, z); each kernel now
      returns sqrt(P(k,0)/P(k,z)) so the product restores P(k, 0)
     * The shared integrand is untouched, so the physical kernels keep their own
      convention
     * Assert that Limber and non-Limber agree at high ell, the invariant that
      exposed it; without the fix the ratio is D^2(zbar), 0.68 at zbar = 0.37
     * Record the extra factor in the Limber contract test

 * Add a kernel for tabulated radial windows

     * Add NcXcorKernelTable, a NcXcorKernelRadial whose window is supplied as a
      table of (chi, W) samples and reconstructed with a NcmSplineBSpline
     * Default to degree 7 rather than a cubic, which caps a 2000-sample window
      four orders short of what the data supports
     * Add a kind property selecting density or shear, the latter carrying
      1/(k chi)^2 and sqrt((l+2)(l+1)l(l-1)) through the new radial vfuncs
     * Trim leading and trailing zero runs so the support is the interval the
      window occupies
     * Expose the reconstruction's breakpoints for a quadrature to align panels
     to
     * Register the type, add it to the umbrella, record the new enum nicks and
      regenerate nc.pyi

 * Route the last two downloads through the shared fetcher

     * Take the download lock for the curated weak-lensing catalogs
     * Publish the native Planck artifacts by rename instead of writing in place
     * Test both, each failing against the code it replaces

 * Let a radial kernel carry shear factors

     * Add NcXcorKernelRadial::eval_kernel_factor, a (chi, k) factor applied to
      the kernel alongside the scale-dependence hook
     * Add NcXcorKernelRadial::eval_prefactor, an ell-only factor applied to the
      whole kernel
     * Apply both on the non-Limber and the Limber paths, the latter at
      k = (l + 1/2) / chi
     * Both default to 1, so every existing kernel is unchanged

 * Promote the analytic kernel base out of tests

     * Rename NcXcorKernelAnalytic to NcXcorKernelRadial and move it, with
      NcXcorKernelRadialKDep and NcXcorKernelRadialKDepGrowth, from
      nc/xcor/tests/ to nc/xcor/
     * Keep the closed-form subclasses named and located as they were; they now
      derive from NcXcorKernelRadial
     * Rename the base's private component type to NcXcorKernelComponentRadial
     * Update ncm_cfg registration, the numcosmo.h umbrella, meson source and
      header lists, the Arb window tooling and the generic ref/free entry
     * Regenerate nc.pyi

 * Cache the downloaded data files in CI

     * Add .github/actions/data-cache, keyed on DATA_CACHE_VERSION and the
     datafile release tag
     * Add .github/scripts/prefetch_data.py, enumerating every asset from the id
     enums
     * Restore the cache in the jobs whose tests reach a downloaded file
     * Write one entry per OS, from build-gcc-ubuntu and build-gcc-macos only
     * Document the cache in CONTRIBUTING.md

 * Cover the old-catalog read path, and refuse a reparametrization that cannot be
     rebuilt

     * Refuse to rebuild a NcmReparam subclass carrying state of its own;
     NcmReparamLinear's unset matrix was dereferenced
     * Report why a rebuild was refused instead of assuming the descriptors were
     at fault
     * Add make_reparam_catalog_fixtures.py, writing pre-migration catalogs
     without the pre-migration library
     * Add three catalog fixtures: HDU0 vardict, HDU0 bare object, and the
     legacy .mset sidecar
     * Test that each container loads, keeps its reparametrization, maps its
     free parameters and keeps its rows
     * Test the rebuild's branches: grow, missing descriptors, index past the
     end, negative index, extra state
     * Test that deserialization, which has no error channel, aborts with the
     reason
     * Stop .gitignore's blanket *.fits excluding the committed truth-table
     catalogs

 * Fix the Planck baseline download race, and make such failures visible (#344)

     Fix the Planck and SNIa data downloads, and surface fatal test failures

     * Add nc_data_download.c: one locked, atomic download shared by both
     datasets
     * Check wget's and tar's exit status; a failed transfer was reported as
     success
     * Download to a per-process temporary and rename, never straight to the
     final name
     * Extract the Planck tarball to a staging directory renamed into place
     atomically
     * Mark the extracted tree complete after the rename, never by a file the
     tarball provides
     * Replace a partial tree left by an interrupted run instead of trusting it
     * Route fatal messages to the duplicated stderr fd, where output capturing
     cannot lose them
     * Give every pytest lane a real timeout and stop the coverage jobs
     disabling them
     * Test the download on a 310-byte asset: concurrency, contention, 404,
     debris, partial tree
 * Let a model resize a reparametrization it outgrew

     Moving Yp out of NcHICosmoDE took it from eight scalar parameters to seven,
     and every file written before that carries an NcHICosmoDEReparamOk saying
     length is eight. ncm_model_set_reparam() checked the reparametrization's
     compatible type but never its length, so it adopted the eight-long vector
     and died inside ncm_vector_memcpy() -- an abort, not an error, so `numcosmo
     catalog analyze` dumped core rather than reporting a file it could not
     read. Of the catalogs here that embed an mset, this is all of them: 636 of
     636.

     A NcmReparam stores no parameter values. Its length, its descriptors and
     its compatible type are all it carries, and the working vector is allocated
     empty at construction and filled from the model by old2new(). So the length
     is bookkeeping that has to agree with whatever model it lands on, and
     rebuilding one at the model's length loses nothing. NcmReparam:length being
     construct-only is why agreeing means a new instance rather than a resize.

     The descriptors are keyed by index into the model's original parameters, so 
     this holds only while those indices still name the same parameters. Every
     one must land inside the new length or the rebuild is refused and the
     existing NCM_MODEL_ERROR_REPARAM_INCOMPATIBLE is raised. That catches a
     descriptor left pointing past the end -- it cannot catch a parameter
     removed from before one, which nothing here could see. NcHICosmoDEReparamOk
     is safe on both counts: it reparametrizes Omega_x at index 2, and Yp was
     index 4, so nothing it names moved.

     The existing BBN migration kept write-only Yp properties on NcHICosmo and
     is tested across obj, bin and yaml. It missed this because a
     reparametrization is a separate serialized object carrying its own copy of
     the parameter count, so none of the three fixtures had one -- the gap was
     object shape, not file format.

 * Stop the xcor lane from exhausting memory

     An xdist worker is its own pytest session, so a module-scoped fixture is
     built once per worker rather than once. The xcor lane's heaviest files hold
     several GB each that way, and nothing releases it as workers drain --
     memory climbs monotonically to the end of the run. On a 12-core machine
     that killed the lane about a third of the time and on 24 workers every
     time. CI never saw it: its runners have 2-4 cores, so the cost there is a
     sixth of what a developer machine pays, and it will stay invisible to CI
     whatever else is added.

     Measured under a cgroup cap, 2171 tests:

       workers  before          after
      12       OOM (>10 GB)    5.3 GB
      24       OOM (>10 GB)    8.4 GB

     Four changes, in descending order of what they were worth.

     The lane runs --dist loadgroup instead of the default load, and the five 
     heaviest files carry an xdist_group mark. Grouped tests go to one worker; 
     everything else still distributes test by test, so the other ten files keep 
     their parallelism. This was the cheapest change and the largest saving: 9.0
     to 5.3 GB at 12 workers, with no test touched.

     test_k_integral.py's frozen fixture is now an LRU rather than an unbounded 
     cache. It held every key it was ever asked for -- 136 of them at 18 MB for
     an easy one and 105 MB for the worst, 4.0 GB by the end of the file.
     Bounding it is what makes 24 workers survive at all, and it costs 1.9x on
     that file: five test functions each walk the whole key set, so the reuse
     distance is a full pass and no cache size or ordering helps. Measured at
     maxsize 8, 16 and 32 and under deterministic ordering, all within 10% of
     each other. XCOR_FROZEN_CACHE overrides it; 0 restores the old behaviour
     and the old speed.

     A Frozen now shares one NcmSBesselIntegratorLevin across both kernels and
     both closures, which is what nc_xcor_solver_solve() does -- it dups one
     integrator per block and hands the same one to every kernel. Giving each
     kernel its own retained four per key and, more to the point, meant the
     suite never exercised the cross-kernel decomposition reuse that
     dev-notes/xcor_ultralevin_batching_plan.md
     §6.1 measures at 8x. Numerically neutral: bit-identical on every block
     method, and 5.7e-14 on KERNEL_GSL, the one method that fits per multipole
     through the kernel's own integrator. That is 175x inside its tolerance.

     Last, cases.ELLS_SUITE moves from [2, 20, 200] to [2, 6, 10, 50]. The bench 
     sweep over all nine block starts says the methods are worst around l =
     4-10, not at the top: exact/chebyshev is 2.3e-08 at l = 6 against 4.7e-10
     at l = 200, and gsl/chebyshev 1.6e-03 at l = 10 against 1.7e-04. The old
     ladder tested the easy end hard and the hard end not at all. Tolerances and
     REFERENCE_FLOOR are re-derived; the latter moves two orders, because l = 6
     is also where the Level-1 reference is least converged on a Chebyshev
     closure -- 1.3e-08 against the 1.0e-10 the old ladder saw.

     Not fixed here: the lane's peak still scales with worker count, and nothing 
     stops the next ample-scoped fixture doing this again. A memory guard on the 
     test wrapper, so a runaway lane is killed rather than the machine, is a 
     follow-up.

 * Updated stale stubs.

 * Put the kernel-space quadratures behind one table

     The three block methods were selected by an if in nc_xcor.c and another in 
     nc_xcor_solver.c, and each carried its own ell-batching wrapper around an 
     otherwise identical loop. They now share _nc_xcor_kernel_space_run() and
     are chosen from one NcXcorKQuad table, which says what each method's block 
     quadrature is, which closure builder it wants, and whether it reports an
     error estimate. Adding a fourth is a line in that table.

     The signature the table needed carries vp_err for every method, not only
     for the exact one that fills it. A method without an estimate is handed
     NULL and nc_xcor_method_has_error_estimate() says so, rather than it
     answering a zero that would read as "no error".

     NC_XCOR_METHOD_KERNEL_GSL_BLOCK is new: qagp broken on the merged knots,
     the rule KERNEL_GSL already ran, over the block closure KERNEL_EXACT and 
     KERNEL_CUBATURE integrate. KERNEL_GSL keeps its per-multipole closure and
     its own runner, untouched. Folding the two together would have changed its
     numbers and cost the one method whose closure is fitted per multipole --
     and it is that independent fit which showed, in #334, that two closures can
     differ by a factor of three where the quadratures do not.

     nc_xcor_integrate_block() is the public entry onto the table: one block,
     from closures the caller already holds. Without it the outer integral could
     not be timed at all, since every path to it also built the closures. It is
     worth 56 ms of closure against 0.05 to 5.1 ms of quadrature, measured --
     the outer integral is between 0.1% and 9% of a block, and everything else
     is the fit.

     That also let the fourth method say something the other three could not. On
     one shared pair of closures, exact reaches 5.3e-13 of the block's peak and
     gsl_block 8.6e-12, against cubature's 1.7e-5: qagp on the merged knots is
     at GL(5)'s floor, and the gap #334 measured was the closure, not the rule.
     It is a diagnostic and not a production method, at forty times exact's cost
     for the same answer -- QUADPACK is scalar, so it computes a whole block per
     node and keeps one value of it.

     Scratch: every buffer in this file is one double per multipole, and a block
     is capped at NC_XCOR_KERNEL_MAX_ELL_BLOCK, so they are fixed local arrays.
     They are local rather than kept on NcXcor because nc_xcor_solver_solve()
     shares xc across an OpenMP team. The one buffer with no such bound, the
     spectral path's folded coefficients, is a GArray grown per cell instead of
     allocated per cell -- that was the only real churn here, hundreds of
     allocations per block rather than fifteen. No hand-made allocation is left
     in the file.

     Every existing number is unchanged: the sweep over the case matrix
     reproduces all 7344 C_ell bit for bit, and the four pre-existing (method,
     closure) worsts come back as 5.3e-13, 1.1e-12, 4.7e-10, 1.7e-5, 1.7e-4 and
     4.1e-5, matching the table P1 recorded them in.

 * Split nc_xcor.c into its three tiers

     nc_xcor.c held the object, the redshift-space Limber methods and the whole 
     outer k quadrature in 2359 lines. The two tiers share no code -- one
     integrates over z with k pinned to (l + 1/2) / chi(z), the other integrates
     over k with a pair of fitted closures -- so they are now separate
     translation units, nc_xcor_limber_z.c and nc_xcor_kquad.c. What stays
     behind is the object, its properties, nc_xcor_compute(), the tier dispatch
     and the Limber-disjoint policy: 849 lines.

     Pure code motion. The five functions the split promoted from static are 
     declared in nc_xcor_priv.h, which now also carries struct _NcXcor and 
     NcXcorArg, and _nc_xcor_check_qag_status() stays in nc_xcor.c because both 
     tiers use it. Nothing else changed: the sweep over the full k-integral case 
     matrix -- 7344 rows, 17 pairs, 2 closures, 3 methods -- reproduces every
     C_ell bit for bit against master.

 * Construction-fixed submodels (#339)

     Make submodels construction-fixed and rework host/submodel handling

     * Add typed, construction-only submodel slots and migrate all
     submodel-bearing models and constructors to use them.
     * Make every submodel construction-fixed: reject post-construction
     attachment, replacement, cross-host binding, and unconsumed submodels
     during MSet loading.
     * Rework `ncm_mset_load()` into construction-time submodel injection while
     preserving the existing on-disk format.
     * Add `ncm_model_peek_host()` and use the host backpointer for submodel
     reparametrizations, including `NcHIReionCambReparamTau`.
     * Add host-aware parameter-name resolution, including `slot:param`
     qualified names and ambiguity detection.
     * Fix MSet batch-update ordering so host models are committed before their
     submodels.
     * Fix recursive serialization of objects nested inside submodels during
     MSet save/load.
     * Make `NcBBN` a typed cosmology submodel and reduce the deprecated
     `Yp`/`Yp-fit` compatibility properties to inert/erroring sinks.
     * Add full constructors for the HICosmo models accepting `reion`, `prim`,
     and `bbn` submodels and migrate C/Python call sites.
     * Expand C/Python coverage for construction rules, MSet loading, host
     lifecycle, parameter resolution, reparametrization, and error paths.
     * Regenerate Python stubs and enum nick truth tables.
     * Fix the `NC_DE_DATA_SIMPLE_ENTRIES` initializer, formatter-mangled
     `G_DEFINE_QUARK` strings, and a use-after-free in qualified parameter
     lookup. 
 * Rename HWLCatalogID to WLCatalogID

     * Drop the survey from the class name, which its members already carry

     The members became HSC_PDR1_HWL16A_* when the enum's identifiers were made 
     reachable, so naming the class after one survey no longer matched what it 
     holds or what it is for.

 * Drop the early dark energy bound and its BBN name

     * Remove nc_hicosmo_de_new_add_bbn and the darkenergy --BBN option it
     served
     * Register the remaining mset function as early_DE_expansion, not BBN
     * Note that NcCBEPrecision:sBBN-file has no effect

     The helper hardwired a 0.942 +/- 0.03 Gaussian on the expansion speed-up at 
     z = 1e9; anyone wanting that bound can build the prior from the function
     with their own numbers. The function itself constrains expansion, not an
     abundance, so sharing the BBN name with the new NcBBN hierarchy would be
     misleading. The unrelated --BBN-Omega_b prior on omega_b h^2 is untouched.

     CLASS reads its own sBBN table only when pth->YHe is left at the _BBN_ 
     sentinel, and _nc_cbe_set_thermo always assigns a value, so that table is
     never consulted and the two BBN paths cannot disagree.

 * Move Yp from the cosmology to NcBBNParametrized

     * Add NcBBNParametrized, carrying Yp as its own parameter
     * Drop the Yp parameter and the Yp_4He write-back from NcHICosmoDE and
     NcHICosmoLCDM
     * Keep write-only Yp and Yp-fit properties on NcHICosmo to read older files
     * State NcHICosmoLCDM's massive neutrino sector as zero, so it can report
     Neff
     * Test the migration against the frozen fixtures in all three formats

     Yp_4He used to compute the abundance from a BBN spline and write it back
     into the model, bumping the cosmology's pkey in the middle of a Boltzmann
     solve and costing a redundant one per parameter change. It is a submodel's
     answer now, cached on the two values it depends on.

     A file written before the move maps by what its Yp fit type meant: free
     becomes NcBBNParametrized at the stored value, fixed becomes the default
     prediction. NcHICosmoLCDM had no BBN branch, so its fixed Yp was taken at
     face value and now changes -- 0.247800 to 0.245262 for the recorded
     fixture.

 * Let a loaded mset replace a default submodel

     * Check a submodel for duplication against the file's own submodels

     A main model may build submodels of its own on construction, so a
     submodel's slot can already be taken by one of those defaults when
     ncm_mset_load() reaches the file's copy. The attach loop replaces it
     correctly; only the check ahead of it treated the default as a conflict.

 * Let a serialized file set a write-only property

     * Require only G_PARAM_WRITABLE when matching a YAML property to its pspec

     Setting a property needs it writable; the reader was applying the writer's 
     condition, so a write-only property that the GVariant paths accept made the 
     YAML path fail with "object do not have property".

 * Give every cosmology a nucleosynthesis submodel

     * Create a default NcBBNParthenope in nc_hicosmo_constructed()
     * Guard it on there being none, so a deserialized submodel is not discarded
     * Cache the submodel and add nc_hicosmo_peek_bbn()
     * Implement Yp_4He once in NcHICosmo, delegating to the submodel

 * Add the NcBBN submodel and its PArthENoPE implementation

     * Add NcBBN, an abstract NcHICosmo submodel predicting primordial
     abundances
     * Require Yp_4He and get_domain; leave DH, He3H and Li7H optional behind
     impl flags
     * Add nc_bbn_check_domain, which errors instead of extrapolating off-table
     * Add NcBBNParthenope, selecting among the three shipped tables by property
     * Share one deserialized spline per table across instances
     * Key the Yp cache on omega_b and DeltaNeff rather than on the cosmology
     pkey

 * Freeze pre-NcBBN serialization fixtures

     * Add tests/tools/make_bbn_compat_fixtures.py, to be run only before the
     migration
     * Store four cosmologies in GVariant text, GVariant binary and YAML
     * Record each case's Yp_4He, stored Yp, ftype, Omega_b0h2 and Neff in
     golden.json

 * Make the HSC catalog identifiers reachable

     * Pin the enum's prefix at NC_GALAXY_WL_OBS_CATALOG, keeping the survey in
     the nick
     * Bind HWLCatalogID to the members by name instead of by position

     The five members derived to the nicks 002, 007, 060, 064 and 094, which are 
     not valid Python names, so HWLCatalogID had to bind them through raw
     integers and no Python caller could reach them at all. Keeping the survey
     and data release in the identifier also leaves room for catalogs from other
     surveys, where a bare field number would be ambiguous.

     This changes the command-line values to hsc-pdr1-hwl16a-002 and so on. The 
     catalogs are new and have never been released.

 * Pin every enum's nick prefix explicitly

     * Add /*< prefix=... >*/ to the 168 enums that left it implicit
     * Merge it into the trigraph the enum already had, where there was one
     * Record every member's nick in data/truth_tables/enums/enum_nicks.json
     * Guard it with test_enum_nicks.py, keyed on the C identifier

     Left implicit, glib-mkenums infers the stripped prefix from whatever the 
     members of an enum happen to share, so adding one member silently renames 
     every other member. Nicks are not only Python names: numcosmo_py.GEnum is a 
     StrEnum over value_nick, so they are command-line values and serialized 
     experiment-file strings too.

     Verified behaviour-preserving: all 827 nicks are byte-identical before and 
     after.

 * Remove NcHICosmoGCG and NcHICosmoIDEM2

     * Delete nc_hicosmo_gcg and nc_hicosmo_idem2 sources and headers
     * Drop their types and reparams from ncm_cfg_register_objects and
     numcosmo.h
     * Remove the gcg/idem2 branches from the darkenergy tool
     * Delete test_hicosmo.py, which covered only these two models
     * Fix a copy-pasted docstring in test_hireion.py
     * Regenerate nc.pyi

 * Drive the synthetic Planck tests from stored spectra

     * Replace the per-test CBE with FixedClBoltzmann in the data-free Planck
     tests, so nothing is solved where only the likelihood assembly is under
     test.
     * Mirror the C accumulation in the simall table-lookup check, which matched
     numpy's pairwise sum only by luck.
     * Ask for a negative TE explicitly in the simall rejection test instead of
     relying on the sign of the CLASS spectrum.
     * Add test_hipert_boltzmann_cbe.py, covering the CBE spectra directly
     through physical invariants at a modest lmax.
     * Add a solve-free synthetic twin of the lensing Boltzmann
     self-configuration test, which otherwise runs only where plc_3.0 is
     present.
     * Keep one acceptance-marked end-to-end evaluation against a real CBE.

 * Gate the Planck data tests on their own marker

     * Add a planck_data marker and --run-planck-data, declared in pytest.ini
     and wired in conftest.py.
     * Move the eleven real-data Planck tests off app, which means
     application/CLI rather than "needs a plc_3.0 tree".
     * Keep app on the two data-backed CLI tests in test_generate.py and add
     planck_data alongside it.
     * Document the marker in TESTING.md.

 * Planck likelihood reimplementation (#332)

     Add native Planck 2018 likelihoods and clik-free experiment generation

     * Add native Planck 2018 likelihoods for plik_lite, SMICA TT/TTTEEE,
     Commander, SimAll, and CMB lensing, including resampling where supported.
     * Add converters from the public plc_3.0 likelihood data and validate
     native likelihoods against clik to machine precision or the expected
     numerical accuracy.
     * Make native likelihoods self-configure their shared Boltzmann
     requirements in prepare(), including a minimum converged multipole range
     for standalone low-ell blocks.
     * Add clik compatibility details needed for faithful reproduction,
     including Commander single-precision pi and SMICA parameter conversions.
     * Backport the upstream clik lensing CMB-marginalized renormalization
     bugfix and document local PLC divergence from upstream.
     * Add clik-free Planck 2018 experiment generation and self-contained native
     likelihood release artifacts that can be downloaded, cached, serialized,
     reloaded, fitted, and resampled without PLC/clik.
     * Register native Planck data types for fresh-process deserialization and
     add provenance, attribution, citation, release-building, and CLI support.
     * Add live clik comparisons, absolute Planck-data golden references,
     synthetic data-free fixtures, and bit-identical golden snapshots covering
     all native likelihoods and release/generator paths.
 * Split the certified integrals on their oscillation scale

     * Bound the outer k-panel width by the integrand's period, not by an
     octave.
     * Subdivide the inner chi-integration on the same scale.
     * Add --resume so a rerun skips entries already certified.

 * Assert C_ell against certified values

     * Commit the certified C_ell truth table, 43 entries.
     * Add Level-2 assertions, scaled by each pair's own magnitude.
     * Reject a generated entry whose enclosure straddles zero.

 * Pick the conditioned form of j_ell

     * Use sqrt(pi/2z) J_{ell+1/2}(z) away from the origin and the 0F1 form near
     it.
     * Generate one (pair, multipole) per task.
     * Take the k range per multipole rather than as the union over them.

 * State the k truncation instead of bounding it badly

     * Drop the |j_ell| <= 1 tail bound, which gives 6e3 against a C_ell of
     2e-5.
     * Emit k_lo and k_hi with every row.
     * Carry the window parameters as data on KernelSpec.
     * Add the truth-table driver over the case matrix.

 * Certify C_ell against Arb, not just the radial integral

     * Add nc_xcor_kquad_arb.c, computing certified C_ell for a pair of analytic
     windows.
     * Move the shared window code to xcor_window_arb.h.
     * Cap the inner integration and fall back to a finite enclosure when it is
     hit.
     * Walk the k range in panels, escalating the precision per panel.

 * Measure the outer k-integral against its own reference

     * Add cases_k_integral.py with the case matrix, the reference and the
     cancellation ratio.
     * Add test_k_integral.py, asserting each kernel-space method at measured
     tolerances.
     * Add bench_k_integral.py, sweeping accuracy, cost and diagnostics over the
     matrix.
     * Compare each method against the closure it actually integrates.
     * Stop each reference cell against the block's peak rather than its own
     value.
     * Add the narrow hard-edged shells N1 to N3, the case the spectral closure
     exists for.

 * Keep the smoothed top-hat accurate in its own tails

     * Evaluate the window with erfc instead of a difference of erfs.
     * Re-enable tophat_smooth in the Arb closure comparison.
     * Use the cancellation-free form in the unit test's reference too.

 * Adding viewer options for different closure modes.

     * Added --closure-type to the xcor kernel viewer, for plots and C_ell
     alike.
     * Added --compare-closure, drawing the other representation beside the
     chosen one.
     * Generalised the viewer's comparison plots to name either alternative.
     * Added viewer tests for both representations and both comparison modes.

 * * Fix a stale panel coefficient bound in a comment.
     * Moved the closure-type choice from NcXcorKernel to NcXcor.
     * Added closure_type arguments to the three kernel closure entry points.
     * Updated xcor tests, the accuracy tutorial and the xcor viewer for the new
     API.
     * Regenerated Python stubs (nc.pyi).

 * Make the panel order cap configurable

     * Replace the compile-time panel order cap with a kernel property.
     * Keep the fitted order-5 cap as the default while allowing callers to tune
     it for different kernels.
     * Use zero to select the default, bit-identical to setting it explicitly.

 * Cap Chebyshev panels at order 5

     * Reduce the panel cap from order 7 to 5, balancing failed-grid cost
     against panel count.
     * Improve solve time across both kernel families, by up to 1.9x for smooth
     kernels.
     * Preserve production-kernel accuracy and Arb agreement at the sampling
     floor.
     * Recalibrate cap-dependent truth-table tolerances.

 * Pool spectral workspaces instead of sharing them per kernel

     * Borrow `NcmSpectral` workspaces from a pool instead of sharing one per
     kernel.
     * Make concurrent closure construction and restriction safe across ell
     blocks and pairs.
     * Reuse pooled workspaces to avoid repeated FFTW planning costs.
     * Verify concurrent computes are bit-identical to serial for both closure
     types.

 * Choose the block integrator inside the block integrator

     * Move the spline/spectral integration choice into
     `_nc_xcor_kernel_integrate_block_exact()`, shared by both `NcXcorSolver`
     and `nc_xcor_compute()`.
     * Fix `NcXcorSolver` aborts with Chebyshev closures.
     * Preserve the solver's closure cache while using the spectral integration
     path.
     * Test that both entry points select the same integrator.

 * Cover the spectral integration path

     * Add tests exercising spectral closures through the full integration path.
     * Check agreement between `KERNEL_EXACT` and cubature on the same closures.
     * Validate spectral results against spline closures built two orders
     tighter, covering both auto and cross spectra.
     * Account for cancellation in the cross spectrum and cubature's own
     tolerance when setting test bounds.

 * Integrate spectral closures exactly on their common panel refinement

     * Add a spectral `KERNEL_EXACT` route that integrates pairs of Chebyshev
     closures exactly on their common panel refinement.
     * Rebase coefficients onto merged cells instead of refitting or evaluating
     new radial solves.
     * Keep the existing merged-knot GL(5) route unchanged for spline closures.
     * Remove per-pair cubature and make Chebyshev `KERNEL_EXACT` substantially
     faster at matched accuracy.
     * Validate agreement between the spectral and quadrature routes on the same
     closures.

 * Check the Chebyshev closure directly against Arb

     * Validate closure fitting errors directly against certified Arb values
     instead of only checking the resulting integral.
     * Show three-to-five-order accuracy improvements for hard-edge and
     heavy-tail kernels, with comparable accuracy for Gaussian and
     disconnected-support kernels.
     * Expose a pre-existing `tophat_smooth` failure in non-Limber closure
     construction that the direct-integral truth-table test does not exercise.

 * Keep the spline closure for Limber multipoles

     * Apply `NcXcorKernel:closure-type` only to the non-Limber closure.
     * Keep the spline closure under Limber, where per-multipole k-space support
     introduces discontinuities that prevent Chebyshev convergence.
     * Fix the public `set_l_limber()` path with Chebyshev closures, which
     previously reached the minimum-panel guard and aborted.

 * Split the Chebyshev closure into panels

     * Split the Chebyshev closure adaptively into panels with a capped
     per-panel order.
     * Recover knot-placement adaptivity while retaining spectral convergence
     within each panel.
     * Reduce high-multipole expansion sizes and pointwise evaluation cost.
     * Preserve robust convergence for narrow top-hat shells where the spline
     closure fails to converge.
     * Add a capped, non-fatal batch variant to `ncm_spectral` for panel
     refinement.

 * Add a Chebyshev representation for the k-space closure

     * Add a Chebyshev closure selected by `NcXcorKernel:closure-type`, keeping
     the spline as the default.
     * Share sampling seeds and domain expansion between spline and Chebyshev
     closures.
     * Set the Chebyshev order from the total phase `k_max chi_max`.
     * Improve accuracy substantially at comparable sample counts, especially
     for top-hat kernels.
     * Expose the expansion through `get_spectral`, while keeping `KERNEL_GSL`
     and `KERNEL_CUBATURE` unchanged through `eval()`.

 * Reuse caller-supplied coefficient matrix

     * Reuse an existing coefficient matrix when it already has the correct
     shape, matching the scalar path.
     * Keep Python behavior unchanged, with `NULL` always producing a new
     matrix.

 * Expand a set of functions in Chebyshev on a shared grid

     * Add `ncm_spectral_compute_chebyshev_coeffs_batch_adaptive()` to expand
     vector-valued functions on a shared Lobatto grid with batched DCT-I.
     * Reuse function evaluations across components and evaluate only new nodes
     when doubling the nested grid.
     * Check convergence per component while keeping a shared expansion order.
     * Keep the batch implementation separate from the scalar Levin path.
     * Validate against the closed-form coefficients of `exp(alpha x)`, `c_n = 2
     I_n(alpha)`.

 * Report the fit error the closure achieved, not the tolerance it was asked for
     (#326)

     Track closure fit residuals and improve spline error handling

     * Track and report measured closure-fit residuals instead of requested
     tolerances.
     * Use measured residuals in `KERNEL_EXACT` error estimates, with
     tolerance-bound fallback.
     * Warn when refinement tolerances differ by more than two orders of
     magnitude.
     * Restrict the equal-tolerance hazard to `KERNEL_CUBATURE`.
     * Implement `eval_idx` for `NcmSplineBSpline`.
     * Add sample-set fit-error estimation without additional function
     evaluations.
     * Recalibrate `nc_xcor_compute_full()` error-estimate documentation.
 * Separate Tutorials from Examples, and fix the CCL two-point timing comparison
     (#325)

     Clean up documentation structure and examples

     * Set chunk errors and warnings to `false` site-wide.
     * Move the eight Python tutorials under Tutorials, rename Worked Examples
     to Examples, and document the tutorial/example distinction.
     * Add `theory/ssc.qmd` to the sidebar and split the landing-page links to
     the tutorial and example indexes.
     * Improve the CCL two-point benchmark: separate construction/evaluation
     timings, use the supplied `ells` range, and document the Limber
     configurations.
     * Fix the ISW `z_` typo.
     * Rewrite flattened cross-references in the ipynb filter and fix the
     `inspect_data_objects` footnote reference.
 * Document the exact k-quadrature (#324)

     Rename KERNEL_FIXED to KERNEL_EXACT and document its accuracy

     * Rename `KERNEL_FIXED` to `KERNEL_EXACT` to reflect its exact outer
     quadrature.
     * Add `nc_xcor_compute_full()` with propagated kernel-building error
     estimates.
     * Document how relative and scaled-absolute tolerances contribute to the
     estimate.
     * Document conditioning, kernel-edge effects, and interpretation of
     cross-spectrum errors.
     * Add a tutorial on choosing and interpreting xcor accuracy settings.
 * Two profiling fixes: the growth special function at z=0, and angular_cl's
     k-grid (#323)

     Optimize analytic xcor test performance

     * Skip the normalized LCDM growth evaluation at `z = 0`.
     * Expose the `angular_cl` kernel `support_tol` and default it to `1e-6`.
     * Reduce the py-xcor lane from 404 s to 31 s without affecting numerical
     agreement or test coverage.
 * Build test_ncm_fit.c once per optimizer group instead of as one binary (#322)

     Parallelize ncm_fit test groups

     * Split `ncm_fit` into five independently runnable test executables.
     * Preserve all 525 tests and standalone all-groups behavior.
     * Reduce the longest C test binary from 177 s to 57 s.
     * Avoid GLib `-p` splitting, which silently drops test coverage.
 * Add closed-form references for the xcor stack, checked against Arb (#321)

     Add closed-form references for xcor validation

     * Add `NcmPowspecAnalytic` with closed-form transfer and growth functions
     for independent xcor testing.
     * Add Arb-based reference generators for analytic power spectra and xcor
     windows.
     * Validate all analytic shapes and window normalizations against certified
     Arb results.
     * Fix numerical cancellation in the analytic lensing window near its far
     edge.
     * Document known Levin panel-edge, kernel-support, and analytic
     power-spectrum convention limitations.
     * Add comprehensive tests for the new analytic power-spectrum
     infrastructure.
 * Fix the measure of the sbessel convenience integrands, and group truth tables
     (#320)

     Organize truth tables by subsystem

     * Group truth tables under their owning subsystems: `sbessel/`, `halo/`,
     `cluster/`, `sphere/`, and `wl/`.
     * Add a README documenting the directory layout, table formats, provenance,
     and how to add new truth tables.
     * Keep installation unchanged since `install_subdir('data', ...)` already
     handles the new hierarchy.
 * Do not mark delegation-only prepare() methods as current (#319)

     * Revert the model-ctrl update in nc_halo_position_prepare() and
      nc_wl_surface_mass_density_prepare(): both only delegate to
      nc_distance_prepare_if_needed(), which carries its own control, so the
      outer guard can never save work and can suppress a needed re-preparation
      when the NcDistance is shared with a second cosmology.
     * Add /nc/halo_position/shared_distance covering that case, which no
      existing test exercised.
 * Adjusting prepare/prepare_if_needed calls.

 * Improving ultra levin (#317)

     Improve spectral edge continuation and tolerance handling

     * Add Chebyshev rebasing with coefficient-norm bounds for stable edge
     continuation.
     * Improve edge-panel continuation, cleanup, growth limits, and batched
     convergence checks.
     * Fix tolerance handling and operator rebuilds in the Levin integrator,
     with CLI validation and expanded tests.
     * Reorganize spectral and test-support code under `ncm/algebra` and `tests`
     subdirectories; remove obsolete tests.
     * Make tight-tolerance tests robust to the ~1e-7 cancellation floor and
     platform-dependent adaptive refinement.
     * Apply formatting and test-code cleanup.
 * CI: update conda lock files

 * Analytic xcor kernels, and a floor under scaled-abstol (#315)

     Analytic Xcor Kernels and Robustness

     * Add analytic xcor kernels with exact (C_\ell).
     * Add Gaussian, top-hat, multimodal, Student-t, power-exponential, smoothed
     top-hat, and lensing kernels.
     * Add analytic scale dependence via `NcXcorKernelAnalyticKDep`.
     * Handle support boundaries explicitly and clamp round-trip endpoints.
     * Add a `1e-6` floor for `scaled-abstol` and document the rationale.
     * Test closed-form windows, (C_\ell), specifications, serialization, and
     scale dependence.
     * Register new types, update the umbrella header, and regenerate
     introspection stubs.
     * Add remaining test coverage.
 * Fix the two CMB ISW aborts in the kernel-space methods (#297, #298) (#314)

     Fix CMB ISW aborts in the kernel-space methods

     * Handle per-component k ranges and integrate consecutive multipoles
     sharing the same domain.
     * Reuse the per-block NcXcorSolver integrator for the cubature method.
     * Split GSL integration at merged spline knots and apply the requested
     NcXcor:reltol directly.
     * Accept non-success GSL statuses when the achieved error still satisfies
     the requested tolerance.
     * Add ISW kernel-method tests and document their accuracy and performance.
     * Reduce Levin solver memory allocation by sizing the initial storage to
     the resolution floor.
     * Cap the non-Limber test power-spectrum range to reduce memory use and
     runtime without changing the comparisons.
 * * Rename the ell loop variable in compute_kernel's f_ell comprehension, flake8
     E741.

     Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>

 * * Guard the pyccl import with importorskip, so the file skips instead of
     failing collection where pyccl is absent.

     Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>

 * Guard batching invariance to within the measured error budget

 * Clamp the CCL bridge block size to what the solver honours

 * Batch the CCL bridge over ell blocks: one factorisation per block,
     ell-dependent scalars reapplied after

 * Add a CCL-facing non-Limber angular_cl backed by the NumCosmo Levin solver

     compute_kernel gains a tolerance-driven B-spline reconstruction, keeps
     CCL's k-dependent transfer, and rejects unsupported der_bessel explicitly.
     Agrees with pyccl's non-Limber C_ell to 0.1% at ell=2, where it differs
     from Limber by 1.9x.

 * * Cover copy_empty: the order is carried over and the copy is independent.
     * Cover deriv_nmax against (order-1)! for x^(order-1), and against
     eval_deriv/eval_deriv2 at orders 2 and 3.
     * Cover the instance name, its rebuild on an order change, and its use in
     ncm_spline_set's min-size error.
     * Factor the child-process call the abort tests share into _run_child.

     Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>

 * * Restore NCM_FFTW_* with monkeypatch in test_cfg.py, so the deliberately
     invalid values no longer leak into the rest of the session.
     * Drop NCM_FFTW_* from the child environment in
     test_impossible_request_fails_loudly.

     Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>

 * * Extract the apt package list into a .github/actions/apt-deps composite
     action, used by both apt jobs.
     * Group the brew package list and document why the Python packages come
     from brew.
     * Resolve gmp's prefix with brew --prefix instead of a hardcoded Cellar
     path.
     * Drop the cfitsio pkg-config debug step.
     * Align the documented apt/brew lists with CI, and note that Ubuntu 24.04's
     GSL is too old.

     Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>

 * * Require GSL >= 2.8: NcmSplineBSpline uses the rewritten gsl_bspline API.
     * Move the apt build jobs to ubuntu-26.04, which ships GSL 2.8; 24.04 has
     2.7.1.
     * Use libgsl-dev instead of the transitional libgsl0-dev.
     * Bound gsl in environment.yml and regenerate the conda lock files.
     * Update the documented GSL requirement.

     Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>

 * Derive the B-spline order from the requested tolerance, and refuse requests the
     samples cannot support

 * Add NcmSplineBSpline: interpolating B-spline of arbitrary order

     Cubic interpolation caps the accuracy of anything built on a fixed table;
     higher orders reach machine precision. Backed by gsl_bspline with a banded
     O(n) solve. Registered for serialization, tested, stubs regenerated.

 * * Recalibrate the construction-tolerance gap assertion from 1e0 to 1e-3.
     * Add a CLI-level guard that --integrator-reltol reaches the computation.

     Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>

 * xcor kernel view: pass integrator tolerances at construction

     --integrator-reltol and --integrator-cheb-reltol were applied with setters 
     after the integrator was built, which the library documents as a no-op: an
     ODE operator keeps the tolerance in force when it was created. The options
     changed the reported values but never the computation. Add a regression
     test, since asserting the reported values cannot catch this.

     Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>

 * xcor kernel view: compute and plot angular power spectra

     Add --cls, --cls-method and --cls-block-size to the view command; C_ell for 
     every auto- and cross-pair over the --ell/--n-ell range, via NcXcorSolver. 
     With --compare-limber the Limber spectra and their fractional difference
     are plotted alongside. Correct the NcmSBesselIntegratorLevin class
     documentation: the endpoint formula was missing its factors of y, and the
     working variable is y = kx, not x.

     Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>

 * sbessel levin: record per-panel contributions for diagnostics

     The integral is a sum of per-panel boundary terms whose individual sizes
     cannot be recovered from the result, so measuring how much they cancel --
     and whether a single panel's own accuracy is the limit -- requires
     recording them.

     Opt-in and zero cost when off, with typed accessors rather than an
     index-keyed array of generic doubles. Records are cleared at the start of
     each integration.

     First use refuted the standing explanation of the accuracy floor: the
     measured cancellation ratio is 678 where the floor implied 1.2e8, the floor
     anti-correlates with it, and a case with no cancellation at all (ratio
     exactly 1.000) still floors at 3e-10.

     Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>

 * sbessel ode solver: restrict the resolution floor to the oscillatory span

     The floor was set from the whole panel width, but the solution only
     oscillates beyond the turning point y ~ sqrt(ell(ell+1)); below it the
     spherical Bessel functions are evanescent and few coefficients suffice.
     Charging every panel for its full width therefore paid for resolution the
     evanescent panels never needed.

     Use only the span beyond the turning point, taking the smallest ell in the
     batch so the floor stays safe for every member. Accuracy and monotonicity
     are unchanged; a full ell=[0,499] block sweep is 11-28% faster, with the
     gain concentrated at high ell where the evanescent panels dominate.

     Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>

 * sbessel ode solver: floor the truncation order at the oscillation count

     The adaptive QR could declare convergence far below the resolution a panel 
     needs. A panel [a,b] in y carries about (b-a)/pi oscillations, and on such
     a panel the leading Chebyshev coefficients are small and nearly flat, so
     the decay test fired on them: a panel spanning 2162 in y converged at 19
     columns instead of the ~1200 required.

     The symptom was non-monotonic rather than merely inaccurate. Panel
     contributions cancel heavily, so at a loose tolerance every panel was
     under-resolved and the errors partly cancelled, while at an intermediate
     tolerance some panels resolved and others did not, destroying the
     cancellation: requesting 1e-8 was ~50x worse than requesting 1e-6, and
     returned the wrong sign.

     The decay test may now only fire once the column count reaches 2*(b-a)/pi.
     On a strongly oscillatory Gaussian this removes the non-monotonicity and
     improves the loose-tolerance error by up to 2600x; accuracy at tight
     tolerances is unchanged.

     Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>

 * Correct CONTRIBUTING: registration lives in ncm_cfg.c, python tests are
     auto-collected by marker

 * update_pyi.sh: run from its own directory and scope black to the stubs

     Run from anywhere else the script wrote the stubs into the caller's
     directory and then reformatted that whole tree, which rewrote unrelated
     sources. It now cds to its own directory, scopes black to the two generated
     files, and refuses to overwrite a good stub with the empty output of a
     failed generation.

 * Update data object inspection source path (#279)

     Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Ssc cubature fallback (#307)

     SSC: improve cubature fallback and tolerance handling

     * Retry p-adaptive cubature failures with h-adaptive subdivision
     * Bound h-adaptive retries to avoid unbounded memory use
     * Make J-PAS Omega_c prior bounds configurable
     * Tune SSC scaled_abstol using Fisher-parameter convergence
     * Set scaled_abstol to 1e-5 to avoid p-adaptive tolerance coincidence
     * Update NcXcorSSCSij default to match SSC tolerance
     * Test fallback with discontinuous integrands
     * Reset TEST splice version
     * Uncrustify
 * Pin conda toolchains and use committed lock files

     * pin compilers to the pre-15 toolchain required by CAMB
     * add lock files for all CI targets to avoid repeated solver backtracking
     * install CI environments directly from lock files
     * add lock consistency checks and automated lock regeneration
     * remove the CI environment cache
     * document the lock workflow

 * Updating to python 3.13 on CI.

 * Adding more tests. Calibrating old tests.

 * hmf: add Castro halo mass function and bias models with validation and
     documentation

     * Add Castro halo mass function and bias with stable logarithmic
     evaluation.
     * Validate both calibrations against CCToolkit and store reference truth
     tables as NumCosmo binfiles.
     * Compute variance derivatives directly from the power-spectrum transform
     and pass evaluation scales explicitly to multiplicity and bias models.
     * Add documentation and published ADS references for the Castro models.

 * jpas_forecast24: cover the varying-Sij paths

     * Add TestSSCSijCalculator: one case per footprint, both missing-area
     errors,
      the invalid sky cut, and agreement with create_covariance_S at the
     fiducial.
     * Add generator cases for `vary_fitting_Sij` on and off.
     * Takes numcosmo_py/experiments/jpas_forecast24.py to 100% coverage.

 * Update stubs and fix a mypy error

     * Regenerate nc.pyi for NcXcorSSCSij and NcDataClusterNCountsGauss:ssc-sij.
     * Restate the dtype in SijCalculator.partial_sky(): np.tensordot is typed
     as
      returning floating[Any] regardless of its operands.

 * * Add `NcXcorSSCSij` as the native NumCosmo SSC (S_{ij}) calculator, replacing
     the Python implementation.
     * Support full-sky and arbitrary-mask calculations using the mask (C_\ell).
     * Add (f_{\rm sky}) area rescaling through the new `area` property.
     * Add cosmology-dependent (S_{ij}) recalculation for
     `NcmDataClusterNCountsGauss`.
     * Keep the resampling matrix fixed while allowing (S_{ij}) to vary during
     fitting.
     * Add `--vary-fitting-sij` to `jpas_forecast24`.
     * Make `NcXcorSSCSij` safe to default-construct and validate completeness
     in `prepare()`.
     * Fix `set_ssc_sij()` so deserialized data remain initialized.
     * Add the required top-hat kernel constructor with an integrator.
     * Add tests against the Python reference, including full-sky, partial-sky,
     serialization, wiring, and construction cases.
     * Add SSC theory documentation and references.

 * ssc: replace PySSC with NumCosmo implementation (#304)

     * ssc: implement NumCosmo SSC and replace PySSC
     * ssc: add arbitrary-mask and cross-mask support with validation
     * xcor: add and improve KERNEL_FIXED integration
     * xcor: respect Limber tiers and batch shared knot panels
     * tests: update J-PAS forecasts and remove PySSC/healpy dependencies
     * CI: fix Python version and older PyGObject compatibility
     * Update stubs and uncrustify
 * Adding fsky and tests.

 * Adding area dependency to Sij.

 * xcor dev-notes: record the ISW Limber step analysis behind GH #297

     * the discontinuity is a Limber artifact, not physics: the exact W_l(k) is
     smooth, Limber inherits the kernel's support edge as a step at nu/xi_max,
     one per multipole
     * measured discriminant: the lensing efficiency vanishes at its support
     edge (1.5e-25) while ISW is still at its plateau (12.5) when truncated at
     recombination, which is why only ISW aborts
     * cross-spectra are immune -- the intersection with a low-z tracer cuts the
     comb off entirely
     * records what was tried and why each attempt failed, and the per-multipole
     interval approach that avoids the shared-endpoint problem

 * xcor: register NcXcorKernelCMBISW and export the two missing xcor headers

     * ncm_cfg: register NC_TYPE_XCOR_KERNEL_CMB_ISW -- every other xcor kernel
     was registered, so the ISW kernel could not be built from a serialized name
     at all and aborted with "object `NcXcorKernelCMBISW' is not registered"
     * numcosmo.h: export nc_xcor_kernel_cmb_isw.h and
     nc_xcor_lensing_efficiency.h, the public include surface required by
     CONTRIBUTING step 3
     * test_kernel_serialization: check in a fresh subprocess that every
     concrete xcor kernel type resolves by name after cfg_init() alone; the
     existing dup_obj tests cannot catch this because they realize the GType by
     constructing the object first

 * xcor: export NcXcorSolver from the umbrella header, and cover the patch's
     untested lines

     * numcosmo.h: include nc_xcor_solver.h -- the header was installed but
     absent from the umbrella, so NcXcorSolver was unreachable from C
     * test_ncm_generic: NcXcorSolver ref/free/clear coverage, which
     introspection never exercises
     * test_solver: peek_block_integrator before and after solve, per-block
     ell-range pinning, and reuse of the pinned integrators across a second
     solve
     * test_solver: set_integrator clones the prototype rather than sharing it,
     replacing it drops the cache, and clearing it falls back to a registered
     kernel's integrator
     * test_solver: several kernels' l_limber thresholds are sorted before
     tiling; overlapping unordered requests still tile contiguously
     * test_solver: a kernel carrying a non-Levin integrator passes the reltol
     closure check
     * test_xcor_view_app: first tests for the xcor kernel view CLI, asserting
     the precision options reach the integrator
     * test_cfg: FFTW wisdom is written once and read back by a later process,
     which needs its own HOME and planner since the paths are no-ops under
     FFTW_ESTIMATE and the loaded-once cache is per process

 * xcor: expose precision knobs in the CLI, and add the design notes and CCL
     benchmark

     * numcosmo xcor kernel view: l-limber and the Levin integrator's reltol,
     cheb-reltol and max-order are settable from the command line
     * dev-notes/xcor_ultralevin_batching_plan.md: architecture, block-size and
     precision measurements, the CCL comparison, and the defects found along the
     way
     * dev-notes/xcor_ultralevin_solver_vs_ccl_bench.py: the benchmark those
     numbers come from
     * .gitignore: ignore the BenchNative build directory used for compiler-flag
     measurements

 * xcor: add NcXcorSolver, batched angular cross-spectra with per-block Levin
     integrator reuse

     * NcXcorSolver: register N kernels once, request the pairs wanted, and
     solve them together; each kernel's ell-block closure is built once and
     shared by every pair touching that block, so cost is O(N_kernels) instead
     of O(N_pairs)
     * one NcmSBesselIntegratorLevin per ell-block, pinned to that block's
     multipole range and kept across solve() calls so its ODE operators' QR
     factorisation survives; injected into the kernel through
     nc_xcor_kernel_get_eval_vectorized_full()
     * solve() runs its block loop under OpenMP: blocks share no operator state,
     so the axis carries no reuse-versus-parallelism tradeoff
     * nc_xcor_kernel: the k-seed array is call-local, so one kernel can be
     evaluated concurrently for different ell blocks; components declare an
     absolute integral scale to the integrator
     * nc_xcor: fixed the tier-1 Limber lower limit -- "zmin ? zmin != 0.0 :
     1.0e-6" parses as "zmin ? (zmin != 0.0) : 1.0e-6", so any non-zero zmin
     became 1.0, in both the GSL and cubature paths
     * nc_xcor_priv.h: internal header for
     _nc_xcor_kernel_integrate_block_cubature(), which takes pre-built
     integrands

 * ncm_cfg: cache FFTW wisdom file I/O across repeated load/save calls

     * wisdom load/save did full file I/O on every plan creation, ~22% of an
     xcor solve batch; the file contents are now cached process-wide and re-read
     only when they change

 * ncm/nc: give adaptive routines an absolute error scale, and stop evaluating out
     of range

     * ncm_spectral: reaching max-order without converging is a fatal error
     instead of a silent truncation returning unconverged coefficients;
     ncm_spectral_compute_chebyshev_coeffs_adaptive_full() takes a
     caller-supplied abstol so a purely relative criterion is no longer the only
     stopping condition
     * ncm_sbessel_integrator: abstol accessors on the base class, letting a
     caller declare the absolute error it tolerates in an integral because it
     knows the larger quantity that integral feeds
     * ncm_sbessel_integrator_levin: the panel Chebyshev floor scales by the
     actual max (b_p |j_l (b_p)|, a_p |j_l (a_p)|) rather than the |j_l| <= 1
     bound; below the l-th turning point j_l is evanescent and the panel cannot
     move the result, so requiring a relative fit of its integrand there is both
     futile and expensive. The scale never exceeds b_p, so the floor is never
     tighter than the bound it replaces
     * ncm_sbessel_ode_solver: clamp the operator-buffer copy to the smaller of
     the old and new sizes; a shrinking ell-range followed by a higher Chebyshev
     order overflowed the fresh allocation
     * ncm_spline: document that ncm_spline_eval() performs no range check and
     extrapolates the boundary interval's polynomial outside the knot range
     * nc_xcor_lensing_efficiency: evaluate the second-order head on [zmax - dz,
     zmax] instead of extrapolating g(z) below the spline's first knot, which
     returned a floor ten orders of magnitude above the true value; the ODE
     initial condition is unchanged
     * nc_xcor_kernel_gal: bound xi_min by the dn/dz lower edge
     * nc_xcor_kernel_cluster_tophat: nc_xcor_kernel_cluster_tophat_new() passes
     the "dist" property name, so the public constructor no longer aborts
     * ncm_function_sample_set: exclude identically-zero components from
     get_absmaxF_min(), so a vanishing multipole cannot collapse the tolerance
     budget of every other component
     * ncm_integral_nd: report the method, dimensions, bounds and tolerances
     when a cubature fails, instead of a bare assertion

 * Testing new additions.

 * Exposing resample type in CLI for WL.

 * Testing conditional addition of derived for PopBeta.

 * Making Beta pop derived parameters conditional on fitting population.

 * Fixing catalog HDU0 dumping.

 * Adding backwards compatibility with 'std_shape'.

 * Making the overwrite of a RNG state an error on initialized catalogs.

 * Fixing empty ObjectArray crashing experiments.

 * Improving error message on deserializing yaml.

 * * Fixed test_run_mcmc_apes_plot_corner_too_many_plot_names: it asserted on a
     substring straddling "--plot-name", but Rich highlights option-looking
     tokens and injects ANSI codes between their characters when color is forced
     (as CI does for Typer's error panel, unlike a plain local terminal),
     breaking the naive substring match. Assert on a dash-free portion of the
     same message instead.

 * * Closed the two files' patch-coverage gaps flagged by Codecov (loading.py 19
     missing -> 0, catalog.py 8 missing -> 0), verified by intersecting
     coverage.json against the PR diff.
     * Removed dead-code guards in PlotCorner and DerivedQuantityError: Click
     already enforces mcmc_file/--variable/--expr as required, so the manual "at
     least one" checks could never fire.
     * Added tests for the real remaining gaps: --include/--exclude column
     filtering (all three branches), an empty-catalog burnin, --plot-name count
     mismatch, --mark-bestfit, catalog visual-hw and param-evolution (previously
     untested commands), a single-chain (run mc) catalog, a missing-catalog-file
     error, and load_catalog()'s direct negative-tail guard.

 * * Fixed 47 mypy errors: LoadCatalog now declares its load_catalog()-derived
     attributes as typed dataclass fields (mcat, mset, functions, etc.) instead
     of injecting them via self.__dict__.update(), which was invisible to the
     type checker.
     * Widened mcat_to_catalog_data's indices parameter to also accept a plain
     list[int], matching what it already accepted at runtime.
     * mypy --exclude '.*meson.*|numcosmo_py/generate_stubs\.py' -p numcosmo_py
     now passes clean.

 * Adding tests and fixing bugs.

     * Fixed a real bug: ncm_serialize_var_dict_to/from_yaml mangled object and
     object-array values, since their GVariant text form isn't valid YAML; now
     routes them through the existing structured node builder/parser instead of
     a naive print/parse round-trip.
     * Extended the NcmVarDict to_from round-trip test with object and
     object-array entries, covering all 5 formats (variant, yaml,
     variant_binfile, variant_file, yaml_file);
     * Added C tests for ncm_mset_catalog_peek_info_from_file and the
     burnin-exceeds-catalog trap.
     * Added a Python test verifying a run's functions array is embedded in the
     catalog and readable back.
     * Added a Python test for plot-corner with multiple catalogs and
     --plot-name.
     * Added a Python test for derived-error rejecting an unknown --variable
     parameter name.
     * Added Python tests for get-best-fit (success path and missing --output).

 * * Catalogs are now self-sufficient: NcmMSetCatalog embeds the model-set and, if
     used, the functions array in FITS HDU0 as a versioned NcmVarDict, no
     experiment file needed to read one back.
     * NcmVarDict gains typed object/object-array set/get accessors, taking an
     explicit NcmSerialize argument.
     * NcmFitESMCMC/NcmFitMC embed the functions array into the catalog instead
     of writing the never-read .oa sidecar.
     * Added ncm_mset_catalog_peek_info_from_file for cheap nrows/nchains lookup
     without a full catalog load.
     * Added g_assert(NCM_IS_MSET) guards after HDU0 deserialization and
     clarified the burnin-exceeds-catalog error message.
     * numcosmo catalog: LoadCatalog no longer requires an experiment file;
     check-m2lnl is the sole exception since it needs a live likelihood.
     * numcosmo catalog plot-corner: now takes multiple catalogs positionally,
     overlaid in one plot, with --plot-name for legend labels; drops
     --extra-experiment/--extra-mcmc-file/--extra-burnin.
     * numcosmo catalog: --burnin now means iterations (ensemble steps) instead
     of raw rows; added --tail to keep only the last N iterations.
     * numcosmo catalog: user-input errors (bad burnin/tail, missing catalog,
     incompatible experiment, etc.) now raise typer.BadParameter for a clean CLI
     message instead of a traceback.
     * Removed tools/mcat_plot_corner, fully superseded by catalog plot-corner.
     * Added C tests for the new NcmVarDict accessors and NcmMSetCatalog
     HDU0/functions-array round trips; updated Python CLI tests for the new
     signatures.
     * Regenerated ncm.pyi and nc.pyi stubs.

 * Reintroducing probability floor.

 * Adding more tests.

 * Improving tests, accuracy and quadring against negative prob.

 * Increased outer region number of knots.

 * FixedQuad: origin-divergence-aware marginal sum, cheaper tail panel, numerical
     fixes.

     - Add NcGalaxyShapePop::exponent_at_origin vfunc (Beta: alpha-1,
     Gauss/GaussLocal: 1.0) so
      _marginal_two_panel can skip the N0/expm1 correction trick for
     non-divergent pops or narrow
      noise disks, using a plain sum instead.
     - Tighten and symmetrize the expm1 safe bound (50.0 -> 0.1, now on
     |delta_ratio|) to avoid
      precision loss in the correction-sum branch.
     - Two-panel domain: give the smooth tail panel its own small fixed-node GL
     table (5 nodes)
      instead of reusing n_radial, cutting cost with no accuracy loss.
     - nc_wl_ellipticity: rewrite shear_at_origin_trace's 1-|target|^2 as
     (1-|target|)(1+|target|)
      to avoid cancellation near |target|=1.
     - On a non-finite/non-positive marginal, dump node-by-node diagnostics (and
     pop model params)
      to stderr/stdout before g_error, and include the bad result value in the
     error message.

 * Increasing default mass upper bound for WL analysis.

 * Ignoring fyaml compilation warnings.

 * Including _GNU_SOURCE in fyaml compilation.

 * FixedQuad: fix alpha<2 Beta population divergence, plus rotation-covariance and
     caching fixes (#292)

     * Rework weak-lensing shape-factor marginalization around an exact psi
      reparametrization and a hybrid psi/native quadrature strategy, improving
      numerical stability and accuracy across the full parameter space, removing
      the previous box/domain geometry, and fixing invalid behavior when
      |eps_obs|>=1.

     * Redesign the intrinsic shape population interface around a single
      r-native probability density contract, updating all population models,
      quadrature implementations, intrinsic-mode optimization, tests, and
      Python bindings accordingly.

     * Improve FixedQuad support for singular Beta populations by introducing
      safer divergence detection, fixed-knot marginal-spline caching for
      alpha<2 populations, and an optional native polar correction near
      chi_I=0, with accompanying regression tests.

     * Make FixedQuad fully rotation-covariant by rotating the reduced shear
      instead of the observed ellipticity, fixing cache invalidation during
      likelihood evaluation, and making all quadrature branches independent of
      the global coordinate frame.

     * Improve marginal-spline performance by adding CLI controls, memoizing
      repeated gt=0 evaluations, and extending regression coverage.

     * Add diagnostic CLI tools to inspect galaxy-shape integrands and validate
      stored catalog likelihoods against the current implementation.

     * Improve documentation, comments, error handling, and test coverage
      throughout the weak-lensing shape-factor implementation.
 * docs: add SNIa+BAO confidence region example

 * Updates and fixes to allow use of modern C standards gnu or strict C.

 * Fixed mypy.

 * Fixing test.

 * Adding support for multiple expressions.

 * Adding support for derived parameters in catalog analyze. Adding unit tests for
     new features.

 * Fixing tests.

 * Specialized kernel for fixed quad WL computation.

 * Fixing lto related warnings.

 * Reworking beta distribution for ellipticity. Now we model the distribution of
     |chi| or |e| instead of |chi|^2 or |e|^2.

 * Removed guard.

 * * Better calibrating initial beta distribution.
     * Fixing angular border problem with read data analysis.

 * Simplified comments and docs.

 * Improving tests.

 * Testing better testing duration cache/restore.

 * Uncrustify.

 * * Fix ESMCMC OpenMP initialization deadlock by replacing ordered retries with
     serial redraw / parallel evaluation rounds, adding bounded retries and
     deterministic RNG ordering.
     * Fix MPI ESMCMC initialization to redraw walkers in place, preserve walker
     indices, eliminate incorrect compaction, and bound retries.
     * Reset ESMCMC walker acceptance flags before initialization to avoid stale
     state across runs.
     * Add `max-iter` limit to the Gaussian transition kernel to prevent
     infinite retries when proposals remain out of bounds.
     * Add deterministic parity tests comparing serial, OpenMP, and MPI
     initialization results.
     * Detect the number of OpenMP threads at test runtime instead of Meson
     configure time, and update the testing infrastructure and documentation.
     * Add support for loading curated HSC weak-lensing catalogs by catalog ID
     with automatic download and caching.
     * Add catalog-wide metadata support to `NcmCatalog` with serialization.
     * Extend the cluster weak-lensing application to load real catalogs,
     validate catalog coverage, and support metadata defaults.
     * Refactor the cluster weak-lensing CLI to separate mock-data generation
     from real-data loading.
     * Add test coverage for real catalog loading and regenerate Python stub
     files.
     * Improve the CLI interface for loading catalogs.

 * * build_check.yml: add TEST_DURATION_CACHE_VERSION to force-invalidate the
     cached per-test duration files used by the C-suite slicer, bumped to 1 to
     rule out a stale/corrupted cache as the cause of ncm_mset_catalog being
     dropped from all three coverage shards on this PR.
     * build_check.yml: print the computed slice-tests.txt contents in the
     "Compute test slice" step for visibility into what each shard actually
     selects.

 * * test_ncm_mset_catalog.c: add C tests for HDU0 mset round-trip, legacy .mset
     sidecar fallback, and the fatal error when neither is present.

 * * ncm_mset_catalog: embed the mset as GVariant binary in the FITS primary HDU
     (HDU0) instead of a separate .mset GKeyFile sidecar, written once at file
     creation.
     * ncm_mset_catalog: keep read-only support for the legacy .mset sidecar for
     old catalog files.
     * Added "numcosmo catalog dump-mset" CLI command to export a catalog's mset
     as YAML.

 * Update requirements

     * Make FFTW3 (double precision) a hard requirement, removing ~90
     optional-dependency guards across the library, tests, and tools. (#286)
     * Make NLopt a hard requirement, removing its optional-dependency guards.
     * Vendor libfyaml v0.9.6 into numcosmo/external/libfyaml (minimal
     core-parser/emitter subset), replacing the flaky system/conda-forge
     dependency.
     * Fix 4 upstream libfyaml bugs in its non-C11 atomics fallback path,
     exposed by NumCosmo's -std=gnu99 build.
     * Update CI workflows, environment.yml, and docs/install.qmd for the new
     hard requirements and vendored libfyaml.
 * * build_check.yml: don't cache an empty {} duration extract. Observed for real
     on this PR's own first (failed) run: "Save test slice durations" runs with
     if: always(), so a run that fails before any test executes (testlog.json
     never created) still caches an empty durations fragment -- restore-keys
     prefix matching picks the most recent entry regardless of content, so the
     *next* run silently inherited zero real duration data and fell back to a
     count-balanced split. Guard with a has_data check so only non-empty
     extracts get cached.

 * * test_slicer.py: fix --durations argparse config -- it was action="append"
     (expects the flag repeated once per file) but the workflow calls it as one
     --durations flag followed by multiple space-separated filenames, which
     needs nargs="+". Caused the first real CI run's "Compute test slice" step
     to fail outright ("unrecognized arguments"). Re-verified against fresh
     complete local duration data with the exact multi-file invocation the
     workflow uses (252.0s/252.9s/252.8s balance, all 76 tests accounted for)
     before pushing.

 * * test_slicer.py: add a summary subcommand and factor the testlog.json
     JSON-lines parsing shared with extract into _read_testlog(), replacing the
     coverage job's inline Python heredoc ("Test timing summary" step) with a
     real, locally-runnable script call.

 * * .github/actions/setup-miniforge: extract the ~30-line Setup miniforge / Cache
     Conda env / Update environment / Save conda-forge cache sequence
     (previously duplicated between build-miniforge and
     build-miniforge-coverage) into a shared composite action, parameterized by
     python-version and optional mpi.
     * .github/scripts/test_slicer.py: new script replacing meson's
     count-balanced --slice K/N for the coverage job's C-tier shards with a
     duration-aware greedy longest-processing-time-first bin-pack (plan
     subcommand), fed by historical per-test durations extracted from
     testlog.json (extract subcommand). Falls back to a uniform/count-balanced
     split when there's no history yet.
     * build_check.yml: wire the C-tier (c-1/c-2/c-3) legs of
     build-miniforge-coverage to compute their test list via test_slicer.py
     instead of --slice K/N, caching each slice's fresh durations (GitHub
     Actions cache, per-slice keys to avoid the 3 concurrent legs racing) for
     the next run to read back.

 * * build_check.yml: also set UCX_TLS=tcp,self,sm -- OMPI_MCA_pml=ob1 only fixes
     OpenMPI's own PML selection, not MPICH's ch4:ucx netmod (which hit the same
     underlying mana-NIC issue with a different error: MPIDI_UCX_init_worker
     "Address not valid"). UCX_TLS is read by UCX itself under either MPI
     implementation, so this is the fix that actually covers mpich.

 * * build_check.yml: set OMPI_MCA_pml=ob1 globally, bypassing OpenMPI's UCX
     transport -- GitHub-hosted runners' paravirtualized "mana" NIC advertises
     IB-like RDMA verbs it doesn't actually support, causing spurious "Failed to
     create UCP worker" flakes on MPI test jobs (real regardless of code
     changes, e.g. PR #283's ncm_fit_esmcmc_mpi ERROR with exit status 0 and all
     TAP subtests ok).

 * * --bench report: prefix the printed line with '# ' to match NumCosmo's usual
     comment-style log/report output convention.

 * * Add a --bench flag (RunCommonOptions in run_fit.py) reporting wall-clock time
     and peak RSS at the end of a run; covers run fit/test/mc/mcmc since they
     all share this base class.
     * run test: also call end_experiment() at the end (was a pre-existing gap
     -- --output/--log-file were silently ignored, and --bench had nothing to
     hook into).
     * Add test_run_bench covering run test and run fit.

 * * test_fit_mc.py: add a use_threads property/getter round-trip test (was only
     exercised via set_use_threads(), not the GObject property system or
     get_use_threads()).
     * test_ncm_fit_esmcmc.c: extend test_ncm_fit_esmcmc_properties() with the
     same use-threads property/getter round-trip coverage, mirroring the
     existing skip-check/log-time-interval pattern.

 * * tests/python/meson.build: exclude omp-marked tests from the xdist fast lane
     (were running under OMP_NUM_THREADS=1, so when 2-3 landed on one worker it
     ran alone for the tail while the rest of the pool idled); they now run only
     in the dedicated single-process, real-OMP-threads pytest-omp lane.
     * build_check.yml: add a py-omp shard to the coverage job's matrix, since
     it was only getting omp-marked test coverage incidentally (under OMP=1) via
     the plain python suite -- excluding them from that suite would have
     silently dropped coverage.

 * * NcmFitMC/NcmFitESMCMC: replace the dead nthreads guint property with a
     use_threads gboolean (set_use_threads/get_use_threads), matching the
     existing NcmStatsDist/NcmFitESMCMCWalkerAPES convention; real thread count
     remains OMP_NUM_THREADS-driven.
     * NcmFitMCMC: delete nthreads and its entire multi-threaded path outright
     -- it was never implemented (g_assert_not_reached in _ncm_fit_mcmc_mt_eval)
     and would abort if ever triggered.
     * NcmFitESMCMC: move the odd-nwalkers validation out of set_nthreads into
     constructed() (it's a Stretch-move requirement, not a threading concern);
     start-run log now probes the real OpenMP thread count live instead of
     echoing a stored value.
     * NcmFitMCBS: ncm_fit_mcbs_run()'s bsmt param changed guint -> gboolean.
     * darkenergy: --mc-nthreads int flag replaced with --mc-use-threads bool
     flag.
     * Updated all CLI, sampling helper, experiment, and example call sites to
     the new API.
     * Updated C and Python tests for the new API; rewrote
     test_fit_mc_keep_order.py's dead-wiring-bug documentation to describe the
     new design instead.
     * Regenerated numcosmo_py/ncm.pyi.

 * * Add NcDataClusterWLFactor's register_shared override, anchoring its obs
     catalog (including per-galaxy pz splines for the Spline redshift scheme)
     instead of deep-copying it per parallel worker.
     * Move the register_shared call from constructed() into start_run(), so a
     dataset swapped in after construction (e.g. via set_obs()) is still
     anchored correctly.
     * Reset the NcmSerialize instance fully in end_run(), releasing shared
     anchors so a long-lived NcmFitMC/NcmFitESMCMC reused across many runs
     doesn't keep stale data alive between them.
     * Add register_shared regression tests (positive and negative control) to
     test_data_cluster_wl_factor.py.

 * * Add NcmData::register_shared vfunc (default no-op) letting a data subclass
     register its own large read-only internals as NcmSerialize anchors.
     * Add ncm_dataset_register_shared() to call it across every NcmData in a
     dataset.
     * Wire it into NcmFitMC/NcmFitESMCMC's constructed(), using each object's
     own internal NcmSerialize, before any per-worker duplication happens.
     * Regenerate ncm.pyi stubs.

 * * Port docs/tutorials/python/cluster_wl_simul.qmd off the deleted legacy
     NcGalaxySD*/NcDataClusterWL classes to the Factor pipeline (fixes the BDocs
     quarto render failure).
     * Update docs/theory/wl_ellipticity.qmd and galaxy_wl_framework.qmd to drop
     dangling legacy gtk-doc cross-references and stale "still being built"
     status text.
     * Drop remaining dangling NcGalaxySD*/NcDataClusterWL comments across
     nc_data_cluster_wl_factor.c, nc_galaxy_shape_pop.{c,h},
     ncm_laurent_series.{c,h}, and their tests.

 * * Regenerate nc.pyi stubs to drop the removed legacy classes.

 * * Drop dangling "matches/direct translation of legacy" comments in
     nc_data_cluster_wl_factor.c and nc_galaxy_shape_factor_var_add.c now that
     legacy is gone.

 * * Rewrite tests/python/numcosmo_py/experiments/test_wl_app.py against the
     Factor-based cluster-wl CLI schema.
     * Port fixtures_xcor.py's LSST bin fixtures to
     GalaxyRedshiftBinning.lsst_srd_edges/compute_dndz.

 * * Rewrite numcosmo_py/experiments/cluster_wl.py to build the NcGalaxy*Factor
     pipeline instead of the legacy NcGalaxySD* classes.
     * Update numcosmo_py/app/generate.py cluster-wl CLI flags to the new
     z_dist/shape_dist schema.
     * Port examples/example_wl_likelihood.py to the Factor API.
     * Migrate xcor kernels.py/view.py from
     GalaxySDTrueRedshiftLSSTSRDType/new_lsst_srd_bins to
     GalaxyRedshiftPopLSSTSRDType/GalaxyRedshiftBinning.

 * * Delete the legacy NcGalaxySD*/NcDataClusterWL galaxy WL pipeline (24 sources)
     superseded by the NcGalaxy*Factor/NcDataClusterWLFactor pipeline.
     * Remove legacy entries from numcosmo/meson.build and tests/c/meson.build,
     and legacy #include/registration in ncm_cfg.c and numcosmo.h.
     * Delete legacy-only C and Python test files and legacy ref/unref coverage
     in test_ncm_generic.c.

 * * Relocate NcDataClusterWLResampleFlag/IntegMethod enums out of the legacy
     header into nc_data_cluster_wl_factor.h.

 * Make faulthandler dump per-test tracebacks repeatedly instead of once.

 * Pack remaining galaxy WL frozen-fixture dicts into truth_tables/wl binfiles.

     Fix missing atol on test_direct_estimate_parity_legacy's near-zero gt/gx 
     comparison (cross-platform summation-order cancellation, unrelated to the 
     above).

     Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>

 * Packing reference data into bin files. Calibrating tests.

 * More tolerance loosing.

 * Black.

 * Loosen bit-exact frozen-fixture tolerances that route through libm-sensitive
     code

 * Calibrating tests for the new WL framework.

 * Galaxy WL calculators (#277)

     Add galaxy Population/Observable/Factor infrastructure for weak-lensing
     likelihoods Add redshift, position, and shape calculator architecture with
     MSet-based model resolution Implement NcDataClusterWLFactor likelihood
     orchestrating position, redshift, and shape factors Add SeriesLensed
     lensed-frame marginalization based on Laurent-series expansions Generalize
     SeriesLensed through pluggable population g-series composition Implement
     spline-based per-galaxy redshift factorization with
     NcGalaxyRedshiftFactorSpline Promote Laurent-series and weak-lensing series
     machinery into reusable library components Redesign Laurent-series
     containers as self-contained reference-counted boxed types Add exact,
     Laplace, CGF, FixedQuad, and SeriesLensed shape marginalization methods Add
     Beta population support and generalized power-series composition Improve
     GaussLocal compatibility through shared population series evaluation 
     Restore ordered Monte Carlo support and fix unordered ESMCMC initialization 
     Bring FIXED_NODES and CUBATURE integration methods to feature parity Add
     adaptive quadrature calibration and per-galaxy lens-node optimization Fix
     cache invalidation, resampling, serialization, and likelihood state
     consistency Reduce allocations and eliminate unnecessary cache rebuilds in
     hot paths Expand unit tests, bootstrap support, and coverage across new
     components Regenerate Python stubs for newly exposed properties Update
     documentation to describe current behavior and add implementation history 
     Apply formatting fixes for compatibility with CI uncrustify version Add
     missing copyright notices and configuration registrations

     ---------

     Co-authored-by: Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 * Cluster richness poisson lognormal (#276)

     * New MassRichCount for Poisson-LogNormal distributions.
     * Tests for the new data object.
     * Integrating new data object into CLI.
     * Adding fits file generator/support.
 * Refactor wl tools (#275)

     * Add shared WL ellipticity transform module and frame utilities
     * Refactor `NcGalaxySDShape` to use shared WL ellipticity kernels
     * Rename and document WL ellipticity frame API
     * Centralize WL frame parity and position-angle handling
     * Improve `nc_halo_position_polar_angles` documentation
     * Complete WL ellipticity API renaming
     * Add missing documentation links
     * Fix mypy typing issues
     * Apply Uncrustify formatting
 * Add build directory and info files to .gitignore

 * Add optimizations for weak lensing calculations (#269)

     * Add projected-radius prefactor API to `NcHaloPosition` for reuse across
     many angle-to-radius conversions.
     * Add `NcWLSurfaceMassDensityLensCtx` and lens-context prep APIs to cache
     lens-redshift-dependent quantities.
     * Add exact `P(z)` normalization virtual method to `NcGalaxySDObsRedshift`.
     * Expose `NcGalaxySDShape` dispatch slots for direct hot-loop access.
     * Add flag-based cache invalidation to `NcGalaxySDShape` prepare paths.
     * Rework cluster WL fixed-node integration to an exact control-variate
     formulation using analytic foreground normalization.
     * Refactor HSM-Gauss shape preparation to cache lens/radius quantities,
     pre-rotate ellipticities, use direct dispatch calls, and remove the
     `z_cl-\epsilon` workaround.
     * Drive shape-cache refresh decisions from model and lens-redshift changes
     in `_prepare`.
     * Switch `_eval_m2lnP_fixed` OpenMP scheduling from dynamic to static.
     * Refactor `NcDataClusterWL` cache management around a force-rebuild flag.
     * Regenerate Python stubs for the new APIs.
     * Add tests for optimized-source reduced-shear calculations.
     * Add tests for shape-data handling with unprepared caches.
     * Fix uninitialized projected-radius prefactor in HSM-Gauss shape
     preparation and add a regression test.
     * Relax reduced-shear cache test tolerances to account for floating-point
     operation-order differences.
     * Formatting, cleanup, and uncrustify updates.

     ---------

     Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Add data-driven local curvature priors for w(z)/q(z) reconstruction (#271)

     * Add data-driven local curvature priors for w(z)/q(z) reconstruction

     Weighted Lp curvature norm in the spline backend plus wspline/qspline 
     accessors and serializable wlp_* MSetFuncs carrying a weight spline. 
     Fisher-information weight builder (numcosmo_py) and generate.py LOCAL_* 
     prior wiring relax the curvature prior where the data constrains the 
     reconstruction and keep it strong where the data is blind.

     * Restrict local curvature weight to the expected Fisher of Gaussian data

     The data-driven weight needs the expected Fisher J^T C^-1 J, which is only 
     defined for likelihoods with an analytic mean vector. Skip non-Gaussian / 
     empirical data (e.g. DataBaoEmpiricalFit) when building the weight; they
     stay in the fit/MCMC likelihood.

     * Expose knot placement via --knots in de-wspline / qspline generators

     Add a KnotPlacement CLI option (default keeps the model default) wiring the 
     new NcHICosmoSplineKnots property through, so reconstruction experiments
     can select uniform vs Chebyshev knots.

     * Fix mypy: narrow model to HICosmoDE before w_de in eval_w

     * Cover local curvature priors: generate LOCAL_*/knots, wlp_kappa eval

     Adds tests for the data-driven LOCAL_KAPPA/LOCAL_D2 generation paths (both 
     spline models, exercising the q-spline weight builder), explicit UNIFORM/ 
     CHEBYSHEV knot placement, the missing-weight guard, and the wlp_kappa 
     MSetFunc accessor (the geometric-curvature C callback).

     * Cover curvature_weight no-other-block and no-knots guards

     * Mark test_generate.py with pytest.mark.app so the -m app shard selects it

     Without the module-level marker the directory-name keyword made the
     --run-app skip fire (masking the gap) but '-m app' deselected the file, so
     its coverage never counted. Matches the other app test files.
 * Add selectable SNIa resample strategy (#272)

     * Add selectable SNIa resample strategy

     Expose an NcDataSNIACovResample property on NcDataSNIACov letting the user 
     pick how mock realizations are drawn: AUTO (light-curve cov when available, 
     else distance-modulus cov), FROM_COV (always available), or FROM_LIGHTCURVE
     (requires the full SALT2 covariance). Replaces the per-resample warning on 
     mu-only datasets (e.g. Pantheon+) with a silent, correct default; an
     explicit FROM_LIGHTCURVE request on a dataset without the light-curve
     covariance is rejected via NC_DATA_SNIA_COV_ERROR_UNAVAILABLE_RESAMPLE.
     * Use ncm_util_set_or_call_error and GLib.Error in SNIa resample review
     * Drop redundant error precondition guard in set_resample_type
     * ncm_util_set_or_call_error already handles an already-set error.
 * Add selectable knot placement to spline reconstruction models (#270)

     Introduce a shared NcHICosmoSplineKnots enum (UNIFORM / CHEBYSHEV) and a 
     CONSTRUCT_ONLY "knots" property on NcHICosmoDEWSpline and NcHICosmoQSpline, 
     so the w(z) / q(z) knot distribution can be chosen explicitly. Defaults 
     preserve current behaviour: Chebyshev (in alpha) for the w-spline, uniform
     (in z) for the q-spline. Chebyshev clusters knots toward the endpoints; 
     uniform spreads them evenly, giving a less degenerate, lower-variance
     high-z boundary. Regression tests for both models and regenerated stubs.
 * Curvature-prior w(z)/q(z) reconstruction toolkit (#268)

     * Add generic spline-curvature backend with Lp norms.
     * Add curvature-prior functions for wspline and qspline.
     * Make curvature priors configurable in reconstruction generators.
     * Add function-space sampler for reconstruction studies.
     * Add projection-bias and reconstruction-target framework.
     * Add reconstruction bands and q-transition observable.
     * Add fiducial-truth injection for Monte Carlo runs.
     * Add optional H0 fitting in the DE wspline generator.
     * Fix serialization of grid-evaluated derived functions.
     * Improve test reliability and infrastructure.
     * Apply formatting, typing, and mypy fixes.
 * Docs overhaul (#267)

     * Overhaul documentation structure, navigation, and contributor guides
     * Migrate scientific derivations from API documentation to dedicated theory
     pages
     * Add theory documentation for cosmological distances, transfer functions,
     recombination, weak lensing, halo profiles, FFTLog, and spectral methods
     * Streamline API documentation with links to theory references
     * Replace broken gtk-doc citations with inline arXiv and DOI references
     throughout the codebase
     * Standardize bibliography management and add citation-format validation to
     CI
     * Update Read the Docs configuration and documentation build workflow
     * Remove obsolete documentation and completed migration-planning artifacts
 * Directory restructure (#266)

     * Move vendored third-party libraries under numcosmo/external/.
     * Reorganize cosmology sources and tests under numcosmo/nc/.
     * Split LSS components into dedicated halo, cluster, wl, and galaxy
     submodules.
     * Move likelihood data objects and related tests to nc/data/.
     * Refactor the NumCosmoMath library into thematic ncm/ subdirectories.
     * Update CI, coverage, documentation, and source references for the new
     layout.
     * Apply uncrustify formatting consistently across first-party sources and
     expand CI coverage.
     * Move test thread management to Meson and add reproducible/flaky test
     modes.
     * Adjust test tolerances and fix minor build, test, and priority issues.
 * Chore/mechanical improvements (#265)

     * Remove Emacs mode lines from all C/H source files
     * Remove accidentally committed SPECTRAL_INTEGRATION_NOTES.md
     * Remove all FIXMEs from C source files with proper documentation
     * Migrate GObject types with clean-cut parents to G_DECLARE macros
     * Migrated 30 types to G_DECLARE_FINAL_TYPE and 5 to
     G_DECLARE_DERIVABLE_TYPE.
     * Run uncrustify and fix derivable type padding counts
     * Padding now reflects 18 minus the number of virtual function slots used,
     leaving room for future additions without breaking ABI.
     * Uncrustify.
 * Fix .mset save with sub_fit writing sub_fit params instead of main fit params
     (fixes #23)

     Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>

 * Add support for fixed node integration in NcDataClusterWL (#252)

     * Introduce fixed-node integration infrastructure, including reusable
     quadrature nodes, precomputed integrals, node-based likelihood evaluation,
     and fixed-node support for galaxy redshift distributions.
     * Add integration-method support to cluster weak-lensing likelihoods and
     expose per-galaxy likelihood diagnostics for method validation.
     * Unify redshift integration support across all cluster weak-lensing
     integration methods and split integrations at the lens redshift to handle
     the reduced-shear discontinuity accurately.
     * Fix fixed-node cache invalidation and refresh logic when cosmology, halo
     redshift, node count, or dependent models change.
     * Move fixed-node cache ownership to NcDataClusterWL and correct
     model-update bookkeeping to avoid stale caches and unnecessary
     recomputation.
     * Add comprehensive validation of cluster weak-lensing integration methods,
     including per-galaxy comparisons, convergence checks, and deterministic
     truth-table tests.
     * Improve coverage of galaxy shape, redshift, reduced-shear cache,
     strong-lensing, and integration tests, including dedicated tests for
     responsivity factors, node-based likelihoods, and cache updates.
     * Optimize cluster weak-lensing test performance by reducing sample sizes,
     fit counts, and redundant Monte Carlo coverage while preserving statistical
     validation.
     * Fix cluster_wl applications, examples, integration tests, required-column
     checks, strong-lensing tests, coordinate-conversion tests, and
     generated-data workflows.
     * Update Python stubs, typing annotations, mypy compliance, formatting, and
     coverage exclusions. 

     ---------

     Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Organize test directories (#264)

     * Reorganize tests into module directories mirroring current source (Stage
     B). C: tests/c/{math,model,data,lss,galaxy}; remove cluster/ and
     cosmology/. Python:
     tests/python/{math,model,lss,perturbations,galaxy,data,xcor,numcosmo_py/*}. 
     Update tests/c/meson.build paths (cosmology_tests->model_tests,
     cluster_tests->lss_tests). Mechanical move only; tiers/suites/markers
     unchanged.

     * Tier heavy tests as acceptance and move them off the default PR lane
     (Stage C). Mark jpas forecast covariance/generate, cluster-richness MCMC
     analyzer, and SNIa best-fit recovery as @pytest.mark.acceptance
     (heavy-only, per-function). Default pytest lane now runs -m 'not
     statistical and not acceptance' with a real timeout; add pytest-acceptance
     lane (kept in the python suite, like pytest-sphere_map, so CI needs no
     shard change). Add --durations=25 to surface slow tests every run. Rename C
     suites c-slow -> c-acceptance (ncm_fit_esmcmc, nc_data_cluster_wl) /
     c-statistical (ncm_stats_dist); drop c-slow from nc_powspec.

     * Split nc_data_cluster_wl into cheap (gen_obs/m2lnP/serialize) and
     expensive (resample/MC) executables via -D, so meson --slice can balance
     them. Drop hand-curated CI shard suites
     (stats-dist/fit-esmcmc/data-cluster-wl/powspec) from tests/c/meson.build;
     keep tier (c-acceptance/c-statistical) and omp tags. Rewrite the coverage
     CI matrix to shard C tests with `meson test --slice K/3` over the tier
     suites; Python lanes stay their own shards. Add per-test c_args support to
     the C test foreach.

     * Split CI into a fast non-coverage lane and a full coverage lane. 
     Non-coverage jobs run FAST_TEST_ARGS (skip acceptance/statistical tiers and
     opt-in capability lanes); the sharded coverage job runs everything, with
     py-acceptance added as its own shard. Move pytest-acceptance to its own
     py-acceptance suite so the fast lane can exclude it. Split nc_powspec into
     analytic (unit) and CLASS-backed cbe (acceptance) executables via -D,
     mirroring the cluster_wl split. Tier the CCL background cross-validation
     (test_background) as acceptance. Add a per-test wall-time job summary to
     the coverage shards.

     * Stop the cluster mass-selection timing benchmarks from repeating: set
     timeit number=1. These tests have no correctness assertions (pure timing
     smoke); one pass over the nsize grid already covers each path, cutting ~30s
     to ~11s on the PR lane. Left the sbessel integrator accuracy tests intact
     (genuine large-interval/high-ell numerics).

     * Fix nc_powspec split: the analytic transfer half is the slow one, not
     cbe. Measured cbe ~4s (splined once) vs transfer (EH/BBKS + halofit +
     corr3d) >2 min; the previous split had the costs backwards and the transfer
     exec hit the 120s default timeout. Rename to nc_powspec_cbe (fast unit,
     fast lane) and nc_powspec_transfer (c-acceptance, unbounded timeout,
     coverage-only); rename the -D macros to POWSPEC_SPLIT_CBE /
     POWSPEC_SPLIT_TRANSFER.
 * Improving tests (#263)

     * ci: pin BLAS/OpenMP threads to 1 and scale meson test processes to runner
     cores

     Pin OMP/OpenBLAS/BLIS/MKL thread env vars to 1 in build_check.yml: fixes
     the apt/pip OpenBLAS deadlock hanging the (ubuntu, apt, pip) job at
     scipy.linalg.solve in test_sbessel_ode_solver, and prevents cores^2
     oversubscription. Replace hardcoded --num-processes=2 with $(getconf
     _NPROCESSORS_ONLN) in all four meson test invocations (GitHub runners are
     now 4-core Linux / 3-core macOS).

     * fftw: skip wisdom I/O under FFTW_ESTIMATE and add planner-config guard
     test

     Early-return from ncm_cfg_load/save_fftw_wisdom when the planner flag is
     FFTW_ESTIMATE: estimate neither reads nor writes useful wisdom, so the file
     I/O and save/load lock are pure overhead (notably in CI). Add
     /ncm/cfg/fftw_planner test covering flag string round-trip, the
     -Dfftw-planner fallback plumbing, NCM_FFTW_PLANNER env override, and
     invalid-flag errors.

     * docs: add TESTING.md test-organization policy

     Define the policy for test layout and labeling: module->directory
     (mirroring the source tree), tier->marker(py)/suite(c) with tiers
     unit/statistical/acceptance, capability->marker + --run-* opt-in. CI shards
     are derived via meson --slice, not hand-labeled. C and Python share
     identical module dir names; Python uses --import-mode=importlib so a math/
     dir does not shadow stdlib. Tests mirror the current source tree and move
     with it when sources are reorganized.

     * test(python): scaffolding for the test-org policy

     Set --import-mode=importlib (pytest no longer puts test dirs on sys.path,
     so a math/ dir cannot shadow stdlib; basenames are unique). Consolidate all
     marker declarations into pytest.ini as the single source and drop the
     duplicate addinivalue_line block from conftest.py. Add
     statistical/acceptance tier markers; declare sphere_map; document that ccl
     is auto-skipped via importorskip (no --run-ccl). Drop the unused slow
     marker and the dead removed_test_py_xcor.py.

     * test: define the omp capability lane and correct the OMP threading model

     Document in TESTING.md that NumCosmo's OpenMP-parallel paths
     (ESMCMC/MC/cluster_wl/stats_dist use `#pragma omp parallel`, governed by
     OMP_NUM_THREADS — there is no separate thread pool) are serialized by the
     default-lane OMP=1 pin, and are covered instead by a dedicated omp lane:
     process parallelism off (--num-processes=1 / no xdist) with OMP_NUM_THREADS
     = available cores. Declare the omp marker. Fix the build_check.yml comment
     that wrongly claimed the ESMCMC thread pool was unaffected, and note omp
     simd (sbessel) is vectorization, not threads.

     * test(sky_match): fix version-dependent ValueError in test_mask

     Mask is a frozen dataclass, so `mask_a == mask_b` is the auto-generated
     dataclass __eq__ comparing (self.mask,) == (other.mask,) -- a tuple
     comparison over numpy arrays whose truth-value resolution is
     numpy/Python-version dependent (returns an array on numpy 2.4.6, raises
     "truth value ambiguous" on the apt/pip stack). Compare the .array
     attributes explicitly with np.all, matching the safe pattern already used a
     few lines above.

     * test: add a dedicated OMP lane to exercise OpenMP-parallel paths

     Mark the tests that drive the `#pragma omp parallel` paths (Python FitMC
     nthreads tests; C ncm_fit_esmcmc, nc_data_cluster_wl enable-parallel,
     ncm_stats_dist use_threads) with the omp suite/marker. Add a pytest-omp
     meson test (no xdist) and a build-miniforge CI step that re-runs the
     omp/py-omp suites with --num-processes=1 and OMP_NUM_THREADS=cores,
     overriding the default-lane OMP=1 pin so the parallel branches are actually
     executed. The same tests still run serially on the default lane.

     * fix(sky_match): make Mask/BestCandidates dataclasses eq=False

     These frozen dataclasses wrap numpy arrays; the auto-generated dataclass
     __eq__ compares (field,) tuples element-wise, which returns an array or
     raises "truth value ambiguous" depending on the numpy version. Set eq=False
     so equality is well-defined identity comparison; content equality is done
     explicitly on .array where needed.
 * Matching by ID  (#261)

     * Refactor sky_match ID matching and result handling
     * Add Jaccard and shared-fraction ID matching methods
     * Add globally optimal one-to-one distance assignment matching
     * Split SkyMatchIDResult from SkyMatchResult
     * Restore scalable per-component matching for sparse catalogs
     * Expand sky_match test coverage and reorganize tests into a package
     * Move sky_match and mock generation into catalog subpackage
     * Add confusion-matrix metrics module with tests
     * Extract geometry and cosmology helpers from mock generation
     * Fix mock catalog physics and HOD sampling behavior
     * Refactor halo and cluster generation into composable stages
     * Add common NcmCatalog table conversion helpers
     * Add NcHaloCatalog with parent/child linkage support
     * Add NcHaloCatalogGenerator for cluster catalog generation
     * Add NcmSkyFootprint for sky-region sampling and density calculations
     * Add NcGalaxyHOD and Zheng07 HOD implementation
     * Add NcHaloCatalogMemberGenerator for galaxy population synthesis
     * Add MockPipeline for end-to-end mock catalog generation
     * Add optional halo radius output and sky-footprint support
     * Improve error handling, type hints, documentation, and tutorials
     * Remove pandas and other unnecessary dependencies
     * Update notebooks, tests, build configuration, and CI support

     ---------

     Co-authored-by: Cinthia Nunes Lima <cinthia.n.lima@uel.br> Co-authored-by:
     henriquelettieri <henrique.cnl@hotmail.com> Co-authored-by: Sandro Dias
     Pinto Vitenti <vitenti@uel.br>
 * nc_multiplicity_func_bhattacharya: add convention enum selecting the a(z)
     redshift evolution (Bhattacharya 2011 vs Heitmann 2019). Add new_full
     constructor and convention get/set; default keeps Bhattacharya 2011. Add C
     and Python tests covering both conventions and serialization. Regenerated
     Python stubs (nc.pyi).

 * nc_galaxy_sd_shape: HSM shape-measurement models and direct shear estimators.

     Rename NcGalaxySDShapeGauss -> NcGalaxySDShapeHSMGaussGlobal and
     NcGalaxySDShapeGaussHSC -> NcGalaxySDShapeHSMGauss (global vs per-galaxy
     shape noise). Add HSM calibration products to the shape models:
     multiplicative bias m and additive bias c1, c2. Add direct shear estimators
     with convention-correct responsivity (none for trace-det/epsilon, applied
     to both components for trace/distortion). Update cluster_wl experiment,
     generate config and example_wl_likelihood for the new shape naming and
     products. Extend galaxy-shape, cluster-wl and generic tests; add
     deterministic estimator tests. Regenerate Python stubs (nc.pyi).

 * Bt mass function (#260)

     * Bhattacharya mass function was implemented.
     * Bhattacharya model constants updated.
     * Mass definition updated.
     * Bhattacharya multiplicity function: fixed missing a(z) redshift
     evolution.
     * Improved Bhattacharya documentation (formula, reference, property docs)
     and formatting.
     * Added Bhattacharya tests against the paper formula, serialization and HMF
     integration.
     * Regenerated Python stubs (nc.pyi).
     * Adding more tests.

     ---------

     Co-authored-by: cinthia <cinthia.n.lima@hotmail.com> Co-authored-by:
     Cinthia Nunes Lima <cinthia.n.lima@uel.br> Co-authored-by: Sandro Dias
     Pinto Vitenti <vitenti@uel.br>
 * jpas_forecast: make photo-z scatter sigma0 configurable (#259)

     * jpas_forecast: make photo-z scatter sigma0 configurable

     Thread a configurable photo-z scatter normalization (sigma0) through the 
     JPAS forecast generation. Add a --cluster-redshift-sigma0 CLI option
     (default 0.1) that flows into create_cluster_redshift for the GAUSS 
     relation; NODIST ignores it. Also fix a duplicated docstring line.
 * Add data object inspection documentation example (#258)

     * Add data object inspection documentation example
     * Improve data object inspection tables

[v0.27.0]
 * Update files for 0.27.0 release

 * Release version 0.27.0

 * Adding codecov configuration.

 * Avoid create commit status for fork PRs.

 * Fixing minor issues with unused variables. (#257)

     * Fixing minor glitches.
     * Fixing indentation.
 * Full site upload from GHA. (#256)

     * Uploading full site from GHA to RTD.
 * Improving documentaton build time and ipynb generation. (#255)
 * More verbose run in RTD.

 * Settling for requests.

 * Reverted to using requests.

 * Trying token instead of Bearer.

 * Handling RTD stable and latest relases.

 * Build artifacts on GHA and uploading to RTD. (#254)

     * Testing a new RTD scheme and validating the integration between GHA
     artifacts and RTD builds.
     * Installing LaTeX on the GHA build environment and identifying missing
     packages required for documentation compilation.
     * Confirming that artifact downloads from GHA are publicly accessible and
     free.
     * Inspecting hashes and fixing artifact SHA handling to ensure integrity
     verification works correctly.
     * Improving the artifact download script, including using `urllib` instead
     of `requests`, packing related functions together, and testing problematic
     ZIP files and edge cases.
     * Removing temporary debug `echo` statements and applying minor fixes to
     the GHA–RTD workflow.
     * Fixing link generation and validating links pointing to RTD builds and
     pages.
     * Adding, testing, and logging multiple RTD triggering approaches,
     including direct API usage, an older existing GitHub Action, alternative
     triggering strategies, and testing behavior with the RTD application
     disabled.
     * Fixing variable names in the workflow and scripts.
     * Removing problematic RTD trigger approaches when necessary and iterating
     on alternatives.
     * Making the authentication token non-mandatory to support more flexible
     execution paths.
 * Csq1d phase support (#253)

     * Adding support for integrating delta theta.
     * Adding support for phase determination.
     * Improving docs.
     * Improving expression for delta theta.
     * Making phase spline preparation optional.
     * Fixing docs and add call to prepare in tests.
 * Moving wspline example to quarto. (#251)

     * Moving wspline example to quarto.
     * Fixing code-fold everywhere.
 * Update min python version to 3.11 (#250)

     * Updating minimal python version to 3.11.
     * Fixing mypy errors in tests.
 * New inspect command. Removing tap support in pytest. (#249)

     * New inspect command.
     * Adding tests.
     * Removing tap support for python tests.
     * Making CI tests more verbose.
     * Updating docs.
     * Using types-tabulate from conda-forge.
     * Using pygobject-stubs from conda-forge.
 * Calibrating pln1d test. (#248)

     Converted assertions to warnings in the timing tests. These tests were
     failing in CI because they run in parallel, which can skew the timing
     measurements.
 * Downloadable docs (#247)

     * Making all relevant quarto documents downloadable.
     * Downloadable scripts.
 * Fix xcor benchmark document. (#246)

     * Updating stubs.
     * Fixing CCL benchmark.
 * Update firecrown import (#245)

     * Updating firecrown connection.
     * Breaking circular dependency.
 * Improve tests to avoid leaving leftovers. (#244)

     * Improve tests to avoid leaving leftovers.
     * Improving tests type hints.
     * Test file creation in cluster analysis.
 * Improving richness analysis tools (#241)

     * Improving documentation and fixing typos.
     * Removing bootstrap after FitMC; adding Poisson noise; improving resample.
     * Updating stubs and tests.
     * Adding ncm_pln1d (Poisson-LogNormal) with tests and conditional testing.
     * Adding cumulative calculation support and tests.
     * Creating abstract Richness interface; implementing Ascaso and Extended
     models with tests.
     * New cluster richness analysis package; reorganizing tests.
     * Fixing exception handling and logging.
     * Improving diagnostics and integer R diagnostics.
     * Adding support for noise, obs_params (Mobs, zobs), and persistence.
     * Fixing tests and unsupported flags.
     * Code formatting (black) and cleanup.
     * Refining Python requirements and CI:

       - conditional/skipped tests (astropy, getdist, healpy)
      - separating sphere_map suite (no MALLOC_PERTURB_)
      - relaxing timing requirements
      - simplifying/installing reqs
      - macOS adjustments and matplotlib via brew
      - pytest debugging

     * Removing leftovers and debug calls.
     * Removing Amazon LLM files.
 * Extend spectral (#243)

     * Adding weighted decompositon in Spectral.
     * Tests for weighted decomp.
     * Uncrustify.
     * Fixing python warnings.
     * Adjusting docstrings.
 * Update install guide (#239)

     * Update install guide for libflint package rename

     Agent-Logs-Url:
     https://github.com/NumCosmo/NumCosmo/sessions/6068332e-c37f-4953-80cb-cdf009eb73c9

     Co-authored-by: vitenti <7767706+vitenti@users.noreply.github.com>

     * Clarify historical note for libflint package rename

     Agent-Logs-Url:
     https://github.com/NumCosmo/NumCosmo/sessions/6068332e-c37f-4953-80cb-cdf009eb73c9

     Co-authored-by: vitenti <7767706+vitenti@users.noreply.github.com>

     * Restricting libfyaml version.

     * Adding restriction to the environment.yml.

     * More verbose builds.

     * Decreasing verbosity.

     * Makeing Cosmology lazy properties.

     * Improving doc.

     * Adding missing test.

     * Saving cache only when creating.

     ---------

     Co-authored-by: copilot-swe-agent[bot]
     <198982749+Copilot@users.noreply.github.com> Co-authored-by: vitenti
     <7767706+vitenti@users.noreply.github.com> Co-authored-by: Sandro Dias
     Pinto Vitenti <vitenti@uel.br>
 * New xcor (#240)

     * Optimize hypot and memcpy of RHS.
     * Optimize row construction and creation.
     * Optimize Givens rotations.
     * Optimize solve (loop unrolling, alignment).
     * General performance improvements.

     * Improve Chebyshev evaluation and adaptive coefficient computation.
     * Improve spectral object and adaptivity.
     * Fix spectral adaptive bugs and return adaptive order.

     * Refactor solve into diagonalization + backsubstitution.
     * Refactor Spectral to use GArray.
     * Refactor OdeSolve and Levin integrator to use GArray.
     * Simplify operator storage.
     * Simplify memory management in ODE solver.
     * Move allocation of rotations to appropriate location.

     * Implement batched solver and integrator.
     * Reuse previous diagonalization and solution state.
     * Store and reuse rotations.
     * Add support for panel reuse.

     * Improve ODE solver interface and finalize implementation.
     * Remove integration interface from ODE solver.
     * Add direct integration for smooth intervals.

     * Transition fully to Levin integrator.
     * Refactor integrator to receive K(x,k).
     * Add callback-based integrand interface.

     * Introduce kernel component object and new interface.
     * Refactor and reorganize kernel code.
     * Unify Limber and non-Limber kernel construction.
     * Extend Limber to union of bounds with aggressive cutoff.
     * Fix kernel normalization and k-factors.
     * Add cluster kernel and related functionality.

     * Improve FunctionSampleSet and spline vector objects.
     * Improve bucket search and range extensions.

     * Improve documentation and docstrings.
     * Update and fix stubs.
     * Adjust tolerance parameters.

     * Improve and reorganize tests.
     * Split and relocate test files.
     * Add regression and truth-table tests.
     * Improve pytest + meson interaction.
     * Fix parallelization and fixture usage.
     * Skip tests when dependencies (e.g., scipy) are missing.
     * Remove unreliable CI tests.

     * Improve LSST-related functionality and constructors.
     * Add galaxy redshift extensions and LSST bins.

     * Improve logging and remove debug prints.
     * Remove obsolete interfaces and indentation spec.
     * Rename variables (e.g., l → ell).

     * Fix bugs (allocation, precision, FFTL, ODE solver, nlopt, reference
     issues).
     * Fix CI, build, and environment issues (libfyaml, macOS, pip/conda).
     * Improve dependency handling and configuration.

     * Improve download routines (Planck data).
     * Improve coverage exclusion.

     * Code formatting and cleanup (uncrustify, remove redundancy).
 * docs: Replace FIXME placeholders with proper documentation (73/588) (#238)

     * Replace FIXMEs in numcosmo/data directory (27/588 complete)
     * Replace FIXMEs in numcosmo/math directory (10 doc FIXMEs complete)
     * Replace FIXMEs in numcosmo/xcor directory (20 doc FIXMEs complete)
     * Replace FIXMEs in nc_hicosmo_de files (13 doc FIXMEs complete)
     * Replace FIXMEs in model header files (3 more complete - 73 total)

     ---------

     Co-authored-by: copilot-swe-agent[bot]
     <198982749+Copilot@users.noreply.github.com> Co-authored-by: vitenti
     <7767706+vitenti@users.noreply.github.com>
 * docs: Replace FIXME placeholders with proper documentation (170/758) (#237)

     * Replace FIXMEs in nc_hicosmo.c set_impl functions
     * Complete FIXME replacement in nc_hicosmo.c (all 83 FIXMEs resolved)
     * Replace FIXMEs in model enum headers (gcg, idem2)
     * Replace FIXMEs in xcor galaxy kernel enum header
     * Replace FIXMEs in xcor CMB and weak lensing kernel headers
     * Replace FIXMEs in data headers (hubble, snia)
     * New black formating.

     ---------

     Co-authored-by: copilot-swe-agent[bot]
     <198982749+Copilot@users.noreply.github.com> Co-authored-by: vitenti
     <7767706+vitenti@users.noreply.github.com> Co-authored-by: Sandro Dias
     Pinto Vitenti <vitenti@uel.br>
 * Mass concentration bhattacharya13 (#234)

     * Mass concentration duffy08
     * Mass concentration bhattacharya13
     * Corrected issues with mdef
     * Mass concentration dutton14
     * c-M dutton14
     * c-M prada12 added to branch
     * c-M diemer15
     * Mass-concentration relations: fixing structure (properties, struct,
     documentation...) and typo in some equations. 
     * Transfer functions:

       - Created No-baryon EH transfer function, to be used in the Diemer15
     mass-concentration relation.
      - Updated (GObject's sintax) all transfer function objects.

     * Watson multiplicity function: fixed the properties list (included
     PROP_LEN).
     * Unit tests for BBKS, EH and "No baryon" EH transfer functions.
     * Power spectrum - new function: derivative with respect to k (mode).
     * Fixing little-h bug in klypin11.
     * Finishing tests.
     * Adjusting parameters.
     * Testing Bhattacharya13.
     * Updating stubs.
     * More tests for concentration.
     * More testing for Mass Concentration.
     * Adding missing tests for meson.build.
     * Testing virial.

     ---------

     Co-authored-by: thaisornellas <thais.ornellas@uel.br> Co-authored-by:
     Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Jpas forecast bug (#235)

     * pulling spherenn
     * jpas forecast
     * fixing the bug of mcmc for low values of omegac
     * fixing small omegac bug for mcmc, documentation and change from lnM to
     lnM-obs in jpas_forecast24
     * Adding tests for NcmSphereMap.
     * Adding refinement support for SphereMap.
     * Adding support for cross correlations.
     * Tests for JPas.
     * Minor tweaks in formatting and docstrings.
     * Using ncm_util_gaussian.
     * Numpy 2.4 breaks healpy.
     * Adding healpy to envs.
     * Conditional test of sphere map.

     ---------

     Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Create notebook and minor fix in the documentation. (#196)

     Fixing documentation and adding examples.

     Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Twofluid tensor (#236)

     * Adding tensor modes plots.
     * Fixing list and filename.
     * Adding support for saving/loading dicts.
     * Fixing typos.
     * Removed leftover debug print.
     * Updated HIPrimTwoFluids to store and use Pk0.
     * Updating stubs.
     * Adding Pk normalization plot.
     * Final plots and sections in bounce_spectra.
     * Cleaning notebooks.
     * Matching output when more stderr is produced.
     * Updating python version.
     * Separating slower tests to new suites.
 * Jpas forecast (#170)

     * Updating notebooks.
     * Moving psf to psf_tophat.
     * Minor improvements.

     ---------

     Co-authored-by: henriquelettieri <henrique.cnl@hotmail.com>
 * Testing Richness proxies (#117)

     * Added and expanded support for **cluster richness–mass analysis**:
     improvements to `nc_data_cluster_mass_rich`, new resampling/apply_cut
     functions, bootstrap support, stability fixes, and multiple algorithm
     optimizations.
     * Introduced new functions for **Ascaso mass–richness calibration**, fixed
     bugs in the Ascaso model, and added corresponding tests.
     * Improved **cluster mass selection**: numerical integration optimizations,
     new limits handling, bug fixes, and expanded test coverage.
     * Added and reorganized various **tests and benchmarks**, including Despali
     halo bias, photo-z Gaussian model, interpolation tests, and general
     fixture-based refactors.
     * General **code quality updates**: documentation fixes, indentation/style
     cleanup, removal of unused models, interface adjustments for bindings, and
     improved CI caching.
     * Removed all Jupyter notebooks and other unnecessary files during cleanup.

     If you want, I can also draft a final squash commit message.

     ---------

     Co-authored-by: Cinthia Nunes Lima <cinthia.n.lima@uel.br> Co-authored-by:
     cinthia <cinthia.n.lima@hotmail.com> Co-authored-by: Henrique Cardoso Naves
     Lettieri <henrique.cnl@hotmail.com>
 * Removed deprecated option. (#233)

     * Removed deprecated option.
     * Avoiding buggy plotnine 0.15.
 * New version 0.26.0

 * Splitting NumCosmo initialization. (#232)

     * Splitting NumCosmo initialization.
     * Initializing object types and functions automatically at numcosmo_py.
     * Updating to python 3.12 in CI.
 * Bounce tutorial (#212)

     * Implementing adiabatic interface in QGRW.
     * New bouncing model tutorial.
     * Fixing bugs in qgrw and qgw.
     * Including missing factor in the power-spectrum.
     * Updating perturbation code.
     * Updating to sundials 7.3.0.
     * Second order WKB for TwoFluids working.
     * Reorganizing code and WKB approximation for TwoFluids.
     * Support for high level interface for two point observables.
     * Fixing minor bounce related terms in qgw.
     * Implementing adiabatic Psi and drho computation for qgrw.
     * Introducing better error handling.
     * Finished the implementation of compute_spectrum.
     * Updated adiab tests.
     * Renaming variable and adding explicit type conversion.
     * Reorganizing documents and footnotes.
     * Introducing more observables to TwoFluids.
     * Adding Abs interface for Complex.
     * New code to compute spectra at different times.
     * Computing tensor spectrum and improving bounce_spectra.
     * Adding prereqs.
     * Adding more tests for QGRW.
     * Updated Ubuntu build.
     * Removed old nc_hipert_wkb.
     * Testing GW powspec interface.
     * Testing plotting tools.
     * Setting timeout to 0 for Vexp.
     * Updated stubs.

 * DE w(z) spline - experiment (#226)

     * Included experiment in generate.py: Dark Energy - w(z) spline.
     * Fixed documentation (description of the model, nickname)  - wspline.
     * Generate DE wspline: Flat universe.
     * Support for curvature calculation;
     * New test for model_de_wspline.

     ---------

     Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Add option to choose non logarithmic integrand (#230)

     * Add toggle to choose between integrand and log integrand
     * Merge `use_lnp` and `use_lnint` toggles
     * Fix integration calculations in Gaussian redshift model
     * Fix nc_galaxy_sd_true_redshift_ln_integ function
     * Refactor integrand function calls to remove unnecessary de-referencing
     and fix function types
     * Add cubature integration checks to nc_galaxy_sd_true_redshift_integ test
     * Rename and refactor cubature integrand type for clarity and consistency
     * Properly test both integrand interfaces for nc_galaxy_sd_position
     * Properly test both integrand interfaces in nc_galaxy_sd_true_redshift
     * Add TODO comments for refactoring spline usage in
     nc_galaxy_sd_obs_redshift_pz
     * Properly test both integrand interfaces in nc_galaxy_sd_obs_redshift
     objects
     * OBS: spec and gauss integrand interfaces diverge on high z and should be
     better handled later
     * Properly test both integrand interfaces in nc_galaxy_sd_shape objects
     * Add monte_carlo_lnint test for nc_data_cluster_wl
     * Improving test timmings.

     ---------

     Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Restricting NumCosmo version and trying texlive-core.

 * Trying r-tinytex.


[v0.25.0]
 * Fixed version test.

 * New release 0.25.0.

 * Environment for NumCosmo use.

 * Adding support for derived quantities in MC runs. (#231)

     * Adding support for derived quantities in MC runs.
     * Rewrote the thread parallelization using OpenMP.
     * Fixed and testing for NcmFitMC.
     * Removing old parallelization code.
     * Testing the new openmp code.
     * Fixing MC function array.
     * Testing MC function array.
     * Adding missing code for serializing data-file.
     * Testing new property.
     * Updated stubs.
 * Updated stubs.

 * Adding support for random walk in APES. (#229)

     * Adding support for random walk in APES.
     * Adding shrink to KDE in APES.
     * Adding interface options in numcosmo command line.
     * Adding interface for shrink control in APES.
     * Now logging for ensemble stats.
     * Removing old comments and rewrapping text.
     * Support for returning all computed quantiles.
     * Adjusting defaults and re-wrapping text.
     * Add critical section to avoid race conditions in stats update.
     * Splitting tests.

 * Updated stubs.

 * Fixing typos.

 * Update stub formatting to match linter conventions. Fixed typos.

 * Updated stubs.

 * Fixing memory leaks. (#228)

     * Fixing minor leaks.
     * Improving suppression file.
     * Fixing typos.
     * Workaround fyaml bug.
 * Improving galaxy integration (#227)

     * Galaxy distribution objects now compute ln_prob.
     * Fixing typos.
     * Support lnint integrator.
     * Support for ln-int integrator.
     * Minor typo fix and help improvement.
     * Updating stubs.
     * Updated non-integrated likelihood.
     * Updated tests.
     * Fixed leaks.
     * Testing ln_int.
     * Add smooth transition for shear g before lens.
     * Fixing test.
 * Improving stub generation. (#225)

     * Improving stub generation to include numpy arrays.
     * Updating stubs.
     * Removing unnecessary conversions.
     * Fixing typos.
     * Using conda instead of mamba in rtd.
 * Stub update (#224)

     * Adding stubs generating code to numcosmo_py.
     * Minimal script to update stubs.
     * Improving type-hints.
     * Excluding the generate stubs from mypy check in CI.
 * Refactor galaxy redshift limit functions (#222)

     * Refactor galaxy redshift limit functions
     * Renamed and implemented new functions for obtaining redshift limits in
     the galaxy redshift objects:
      - Changed `nc_galaxy_sd_obs_redshift_get_lim` to
     `nc_galaxy_sd_obs_redshift_get_integ_lim` for integration limits.
      - Added `nc_galaxy_sd_obs_redshift_get_zp_lim` for distribution limits.
      - Updated related classes and structures to accommodate the new function
     signatures and ensure proper functionality.
     * Adding helper functions for gaussian integral computation.
     * Moved local implementation to the object where it is related to.
     * Using correct normalization for z-gauss.
     * Fixing method/variable names.
     * Computing the log of the normalization.
     * Updating tests.
     * Improving tests for gaussian integral.
     * Testing missing lines.

     ---------

     Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Adding support for setting seed to a MC. (#221)

     * Adding support for setting seed to a MC.

     * Support for setting fitting precision.

     * Removing timeout of skymatch.
 * Updated stubs.

 * New P(z) redshift object (#190)

     * Add nc_galaxy_sd_obs_redshift_pz object and data access functions
     * Implement nc_galaxy_sd_obs_redshift_pz_gen and add corresponding tests
     * Add nc_galaxy_sd_obs_redshift_pz to documentation and tutorials
     * Add prepare and get_lim methods to nc_galaxy_sd_obs_redshift object
     * Refactor: simplify redshift bounds retrieval using ncm_spline_get_bounds
     * Add sd_obs_redshift_pz tests to test_nc_data_cluster_wl
     * Remove redundant or obsolete pz tests
     * Revamp sd_obs_redshift_pz tests
     * Fix gen, setget, and integration tests for redshift_pz
     * Fix memory leaks and segmentation fault in redshift_pz
     * Fix: update RA field range in HaloPositionData and test fixtures
     * Fix: properly handle coordinate systems in redshift and shape generation
     * Fix: adjust redshift range and improve conversion handling in tests
     * Fix: always update shape on resample
     * Fix: update normalization and integration calculations for redshift
     distributions
     * Fix: remove debug prints and update test normalization
     * Fix: reset mass, ra and dec before resample
     * Fix: use lnint during integration
     * Fix: return GSL_NEGINF for out-of-bounds integrand evaluation
     * Refactor: remove unused functions and headers from sd_true_redshift
     * Add sd_shape_gauss_hsc object and initial tests
     * Add bias, Jacobian, and ellipticity convention support to
     sd_shape_gauss_hsc
     * Add gen1: generate a single galaxy observed redshift
     * Add tests for redshift and shape generation including bad configurations
     * Add strong lensing tests for galaxy shape models
     * Add support for smooth center galaxy cut and truncated redshift
     * Refactor galaxy shape parameters: sigma_int, sigma_true, sigma_obs →
     sigma, std_shape, std_noise
     * Rename shape parameters: e_rms → sigma, e_sigma → std_noise
     * Refactor: make std_shape a NcGalaxySDShapeGauss property
     * Fix: set std_shape from sigma during preparation and generation
     * Fix: adjust handling of Euclidean coordinates in shape calculations
     * Fix: apply bias correction for z < z_cl in shape_gauss_hsc
     * Fix: update tests to cover coordinate conversion and shape noise
     * Fix: update and reorganize test structure for galaxy shape statistics
     * Add cluster WL generation and configuration tests
     * Add validation for cluster model parameters
     * Add signal-to-noise ratio computation to nc_data_cluster_wl
     * Add simplified normal likelihood with upper/lower limit control
     * Add parallelization support and option to disable it
     * Add tools for command-line parsing
     * Add support for coordinate system conversion and related tests
     * Add fixtures for test configuration and bad parameter validation
     * Refactor: simplify resampling logic and eliminate redundant code
     * Refactor: reduce number of Monte Carlo and shape statistics tests
     * Refactor: isolate monte_carlo tests (later reverted)
     * Revert "Reformulate sigma_int"
     * Revert "move monte_carlo tests to its own test suite"
     * Update documentation, tutorial, and environment files
     * Update stubs and fix mypy issues
     * Style: apply black and remove coverage lines
     * Fix: update copyright
     * Fix: CI debugging, typos, and include misuse
     * Fix: improve robustness of Newton iteration
     * Fix: increase sample size and better calibration in MC tests
     * Fix: handling low-probability galaxies and bias estimates
     * Fix: use ncm_model_ctrl_update for model_ctrl
     * Refactor: add LCOV exclusion comments for coverage
     * Ran black

     ---------

     Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * New DES Y5 SNIA support. (#215)

     * New DES Y5 SNIA support.
     * Adding support for qspline experiments.
     * Generate: inclued XCDM experiment.
     * Test: obtain best-fit Flat-wcdm amd compare to DES-Y5 SNeIa result.
     * Dataset: hicosmo.py - Create all_combined_JUN_2025 for BAO data (DESI
     DR2).
     * Aded "@NC_DATA_BAO_DVR_DTDH_DESI_DR2_2025" in the documentation.
     * Added: test for Pantheon+ supernovae data; test_py_generate.py.
     * Updating stubs.
     * Fixing mypy issues.

     ---------

     Co-authored-by: Mariana Penna Lima <pennalima@gmail.com>
 * Adding google analytics to NumCosmo site. (#220)

     * Adding google analytics to NumCosmo site.

     * Adding more control of GA.
 * Improving reltol for numerical int for fftlog test. (#219)

     * Improving reltol for numerical int for fftlog test.

     * Falling back to higher precision integration when failing.
 * Update enums and reqs. (#218)


 * Adding lintegrate (#217)

     * Support for lnint integral using lintegrate.

     * Meson build file.

     * Adding gsl dependency to lintegrate.

     * Testing NcmComplex and removing lintegrate from coverage.

     * Adding more tests for NcmComplex.

     * Adding more tests for NcmComplex.
 * BAO data - DESI  DR2 2025. (#216)

     * BAO data - DESI 2025.
     * Include DESI DR2 BAO data (obj file).
 * * Included new CC data in the enumerator (nc_data_hubble). (#214)

     * Included new CC data in the enumerator (nc_data_hubble).
     * Added test for nc_data_hubble object.
     * Added test_nc_data_hubble

 * Updated stubs.

 * Included new Cosmic Chronometers H(z) data: (#209)

     * Ratsimbazafy et al. (2017) (arXiv:1702.00418):
     nc_data_hubble_ratsimbazafy2017.obj
     * Borghi et al. (2022) (arXiv:2110.04304): nc_data_hubble_borghi2022.obj
     * Jimenez et al. (2023) (arXiv:2306.11425): nc_data_hubble_jimenez2023.obj
     * Jiao et al. (2023) (arXiv:2205.05701): nc_data_hubble_jiao2023.obj
     * Tomasetti (2023) (arXiv:2305.16387): nc_data_hubble_tomasetti2023.obj

 * Bao data desi (#213)

     * Added DESI DR1 BAO data (Adame et al. 2024).
     * Created new data bao object adapted to the LRG and ELG DESI data.
     * Respective tests were created. Updated nc_data_bao_rdv object.
     * Added DESI DR1 BAO data: BGS, QSO, Lyman alpha.
     * Updated object, nc_data_bao_dtr_dhr, and created its respective test.
     * All BAO data objects are updated. 
     * Improved the BAO tests.
 * Updating sundials to 7.2.1. (#211)

     Updated objects to use new sundials interface. Fixed minor bug in
     experiments/planck18.py.
 * Fixing compilation glitches. (#210)

     * Fixing minor compilation warnings.
     * Replacing mamba with conda.
     * Avoiding buggy version of cfitsio.
 * Updated README links.

 * Better calibration for two-point limber. (#208)

     * Better calibration for two-point limber.
     * Testing number counts extrapolation.
     * Fixed interval tests.
 * Benchmark two-point calculations.  (#207)

     * New two-point comparison benchmark document.
     * Adding more unit testing.
     * Added support to cubature vector integration.
     * Integrating ells in blocks.
     * Adjusting tolerance.
     * Adding timing.
     * Support for optimization flags.
     * Creating a intermediate spline for performance.
     * Updating tests and stubs.
     * Improving documentation.
     * Improving benchmark.
     * Calibrating tests.
     * More testing for xcor.
     * Testing setting number of points in CCL tracers.
     * Disabling debug in coverage run.
     * Removing edge case.
     * Modified error treatment.
 * Benchmark PowerSpectra (#206)

     * Reorganizing documents.
     * Initial version of CCL/NumCosmo PS comparison.
     * New ccl_power_spectrum.qmd benchmark.
     * Added CCL compatibility mode to NcHICosmoDE.
     * Updated documents to use shared _setup_models.qmd.
     * New unit tests for CCL power spectrum comparison.
     * Merging unit test files.
     * Excluding unreachable lines from coverage.
 * Removing outdate/unsupported old code. (#205)


 * Benchmarks (#204)

     * Completed ccl_background.qmd implementation  
     * Removed old notebook  
     * Vectorization improvements: renamed NcDistance methods, spline bucket
     search, refactored comparison code  
     * Testing and calibration: spline evaluation, CCL/xcor fixtures, kernel
     adjustments, test intervals  
     * Quarto configuration: fixed getdist.plots conflict (external library
     fix), output logging, Jupyter/MATPLOTLIB adjustments  
     * CCL enhancements: background tests, Omega_x comparisons, high-precision
     calibration  
     * Neutrino model testing: added edge cases, parameter updates (A_s vs
     sigma8)  
     * Output fixes: disabled LaTeX, rc_params updates, matplotlib behavior
     fixes  
     * Maintenance: stub updates, dependency fixes, menu reorganization 
 * Matching algorithm (#203)

     * First version of matching algorithm
     * Adding support for repeated objects.
     * Adding support for batch search.
     * Updated stubs.
     * Testing batch insert and batch search.
     * Adding vectorized distances calculation.
     * Adding new tests for distances.
     * Updating stubs.
     * Refactored SkyMatch.
     * Updated tests.
     * Testing large inserts.
     * Adding missing distance tests.
     * Testing mask and calibrating tolerances.
     * More testing for distances.
     * Removing impossible lines from coverage.
     * Testing dist with model implementing Dc.
     * Testing cross matching.
     * Adding support to assert_never.
     * Testing qconst with curvature.
     * Testing no mask calls.
     * Adding SkyMatch tutorial.
     * Fixing tutorial seed.

     ---------

     Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Mass concentration duffy08 (#202)

     * Mass concentration duffy08
     * Corrected issues with mdef
     * Allowing concentration-mass relations to be redshift dependent.
     * Generalizing Einasto dl_sphere_mass for any radius.
     * Created function: nc_halo_density_profile_eval_spher_mass_delta (compute
     the enclosed mass within radius R_Delta).
     * nc_halo_density_profile_eval_numint_dl_spher_mass: it computes now the
     enclosed mass for any radius (not only R_Delta).
     * Updated the respective tests.
     * Testing concentration mass relations.
     * Adding setting hooks for massdef and Delta for HaloMassSummary.
     * Updating stubs.
     * Adding missing test for dl_spher_mass_s.
     * Adding test to meson.build
     * Adding missing tests.
     * Fixed typo in tests.

     ---------

     Co-authored-by: thaisornellas <thais.ornellas@uel.br> Co-authored-by:
     Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Moving examples to docs (#201)

     * Adding pandas to conda environment.
     * Moving recobination example.
     * Improving docs.
     * Removing old examples and left-overs.
     * Moved hiprim example.
     * Testing new helper functions.
     * Moving example_hiprim_tensor_modes example.
     * Moving example_ps.
     * Moving example_epdf1d.
     * Reorganized plotting, adding auto-thinning.
     * Testing plotting tools.
     * Use tex only if available.
     * Testing numcosmo catalog plot-corner.
 * Bumping to new version in development.


[v0.24.0]
 * Updating version in pyproject.toml

 * Bumping minor version. (#200)

     * Bumping minor version.

     * Updated changelog.
 * Improving docs (#199)

     * Improving doc building.
     * Improving simple example.
 * Imported Despali Mass function from jpas-forecast. (#198)

     * Imported Despali Mass function from jpas-forecast.
     * Adding to numcosmo infra.
 * Moving to a quarto generated documentation (#197)

     * Updated documentation: migrated from `gtkdoc` to `gi-docgen`, reorganized
     structure, and added ReadTheDocs configuration and dependencies.  
     * Improved CI: added environment variables, dependencies, and optimized
     builds for documentation and tests.  
     * Enhanced tests: added mass function tests, fixed Omega_m calculation, and
     optimized with reusable fixtures.  
     * Bug fixes: resolved issues in Bocquet and Watson models, including Delta
     type assertions.  
     * Cleaned up: removed old site, leftovers, and improved notebook and code
     organization.  

 * Coverage python (#185)

     * Test uploading python coverage data.
     * Uploading coverage.xml.
     * Adding pytest-cov to requirements.
     * Fixing multiple files notation for codecov.
     * Adding missing coverage package.
     * Fixing coverage directory.
     * Removing parallel testing when getting python coverage.
     * Adding support for priority.
     * Reorganizing python functions.
     * Working on cosmology refactoring.
     * Tweaking slow tests.
     * Fixed fparam_get_pi_by_name => param_get_by_full_name.
     * Calibrating tests.
     * New tests for Cosmology class.
     * Testing missing lines in Cosmology.
     * Removing external python tools from coverage.
     * Updated stubs.
     * Adding missing tests for sky_match.
     * Renamed class NcSkyMatching => SkyMatch.
 * Matching algorithm (#191)

     * first version of matching class
     * Fixing mypy hints.
     * Adding zlib on the dependency list.
     * Reorganizing matching code.
     * Renaming methods and moving module.
     * Tests for sky_match.
     * More tests.
     * Installing libfabric-devel manually.
     * Black.
     * Updated numcosmo usage.
     
     ---------
      Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Update stubs (#195)

     * Updating stubs.
     
     * Fixing unused variable and comments.
 * Update spline func tests (#194)

     * Updating tests for AutoKnots.
     
     * Fixing pi factor.
 * Fix docs (#193)


 * Mass concentration klypin11 (#192)

     * Created Klypin et al. (2011) concentration-mass relation.
     * Renamed object (nc_halo_cm_param).
     * Update test_ncm_generic.c
     * Fix test.
     * Improved unit test - Klypin11.
     
     ---------
      Co-authored-by: Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Update weak lensing framework (#189)

     * Fixed generated for cluster wl.
     * Added resample support for nc_data_cluster_wl.
     * Adding support for bootstrap in nc_data_cluster_wl.

[v0.23.0]
 * (Re)updating changelog.

 * Better handling of git hash.

 * Updated changelog.

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
 * Mark every 0.27 compatibility piece with "dropped in 1.0"
 * Place stored parameter descriptions by name, so files from 0.27 load correctly
 * Tighten the wording of the data directory and download lock docs
 * Call the per-user directory the NumCosmo user data directory
 * Add NCM_CFG_HOME_ENV for the NUMCOSMO_HOME variable
 * Look up Planck baseline files only in the NumCosmo data directory
 * Pin the CI data directory to the XDG default and cache it there
 * Point downloads into the legacy ~/.numcosmo at the XDG location
 * Keep an existing ~/.numcosmo, otherwise follow the XDG spec
 * Take over a download lock whose owner is dead on this host
 * Release the download lock before a failed fetch aborts
 * Keep the Python test children out of the real data directory
 * Make the ncm_cfg spawn helper robust and its GLib-DEBUG checks real
 * Make the ncm_cfg child tests able to fail
 * Regenerate the Python stubs for the 0.27 compatibility shims
 * Keep the 0.27 API that Firecrown uses working until 1.0
 * Simplify the CI workflow and its comments
 * fix(test): isolate ncm_cfg child environments to avoid GLib warnings
 * add(cfg): support configurable data directory
 * Check the Python stubs in CI and move the mypy config to pyproject.toml
 * Restore the blank lines lost with the pylint comments
 * Replace flake8 and pylint with ruff
 * Type-check untyped function bodies with mypy and fix what it finds
 * Move the CI quick checks into two lint jobs
 * Improve NcmDiff robustness, convergence and derivative methods
 * chore: ignore local conda environment
 * Restore lcov configuration.
 * Make the catalog evidence estimator selectable and fix its short-catalog hang
     and bootstrap
 * Cover the catalog, MPI serial path and flat kernel, and warn on a singular
     catalog covariance
 * Give invalid least-squares trial steps a finite residual
 * Report the GSL least-squares solver status when a run fails
 * Return the refill input buffers in the async MPI job loop and finalize MPI
     after its own exit handlers
 * Locate a mode at zero in NcmStatsDist1d and draw the EPDF test sample in
     sequence
 * Fix clang and GIR warnings, the TAP noise and three-rank MPI tests on two cores
 * Handle the empty catalog where the best-fit row is used
 * Fix the NcmCSQ1D propagator shift from a nonzero initial time and test the
     propagator in C
 * Fix the NcmCSQ1D boost series stop and fail on boosts beyond double precision;
     test the frames in C
 * Fix the NcmCSQ1D adiabatic maximum at an interval end and test the adiabatic
     vacua against Hankel
 * Read the NcmCSQ1D evolution from the time of ad hoc conditions and test the
     evolution in C
 * Give the NcmCSQ1D phase an origin at the start of the evolution, defined on
     both sides of it
 * Document the NcmCSQ1D settings and prepare from the code and drop an unused
     model control
 * Make the NcmCSQ1DState distance and circle exact near coincident points and
     give get_phi_Pphi the phase of theta
 * Fix the propagator solver check in NcmCSQ1D and document its methods from the
     code
 * Document NcmMPIJobTest from the code and drop its unused state
 * Keep unevaluated proposals out of the MPI step statistics and document
     NcmMPIJobMCMC
 * Document NcmMPIJobFEval from the code and test its returns across ranks
 * Document NcmMPIJobFit from the code and test its returns across ranks
 * Move the MPI worker loop to ncm_mpi_slave.c and test it with raw MPI
 * Test the MPI job protocol across ranks, shapes and datatypes
 * Receive MPI job returns with the return datatype and run job arrays without MPI
     support
 * Keep short catalogs whole in the burn-in criteria and document the catalog
     trimming from the code
 * Fix the catalog distribution and interval functions on short catalogs and
     document them from the code
 * Return the error of a cached catalog log evidence and document the evidence
     estimators
 * Document the NcmMSetCatalog getters and row input from the code
 * Fix the catalog reset of the per-ensemble arrays and a generator leak, explain
     K_eff
 * Fix the ESMCMC diagnostic columns and document NcmFitESMCMC from the code
 * Document NcmFitESMCMCWalkerAPES from the code
 * Fix the multi-stretch acceptance and the walk move, set the walkers' dimensions
 * Fix the catalog kernel with repeated rows and bound its sampling, document the
     kernels
 * Fix NcmFitMCMC acceptance: reject invalid proposals and keep the Gauss kernel
     symmetric
 * Fix NcmFitMCBS construction and bootstrap type, document and test it
 * Document NcmFitMC and test its estimator against the Gaussian resampling
 * Fix NcmLHRatio2d's border, drop its unused parts and test it against the
     Gaussian profile
 * Fix NcmLHRatio1d after a Fisher matrix, drop its unused parts and test it
 * Port NcmFitGSLLS to gsl_multifit_nlinear and fix the fit backends' resets and
     results
 * Fix NcmFit covariance and accurate derivatives, the NcmDiff Hessian round-off
     and the levmar stopping tests
 * Fix a leak in NcmFit restarts without free parameters and document running and
     logging
 * Fix NcmFit covariance lookups with a parameter that was not fit
 * Fix the NcmFit degrees of freedom after a reset and the levmar workspace
 * Fix NcmFitState least-squares m2lnL and its gradient
 * Fix NcmLikelihood annotations and error text and document the posterior
 * Make the NcmPriorFlat least-squares form match its -2 ln P and document the
     priors
 * Fix NcmFunctionSampleSet domain expansion at a hard limit and the
     adaptive-midpoint messages
 * Document NcmStatsDistVKDE and drop its dead code
 * Fix the NcmStatsDistKDE LSCV closed form under shrinkage and the
     fixed-covariance reprepare
 * Fix NcmStatsDist evaluation after a defensive-frac change and document the API
 * Make the NcmStatsDist acceptance objective leave-one-out and recover rejected
     fits
 * Fix the NcmStatsDist split fallback and size guard, make center shrinkage exact
 * Document the NcmStatsDist properties and correct the auto-kernel scope
 * Document the NcmStatsDistKernel family and construct the Student-t nu default
 * Document NcmStatsAcorr more precisely and report the cap from every estimator
 * Fix NcmStatsVec reset and strided input, and document its row semantics
 * Make the NcmStatsDist2d defaults abort and trim NcmStatsDist2dSpline to what it
     implements
 * Speed up the NcmStatsDist1dEPDF bandwidth selector and document the estimator
 * Continue NcmStatsDist1dSpline tails to zero and add a density spline
 * Rework the NcmStatsDist1d CDF and its inverse, and add Steffen splines
 * Fix NcmBootstrap realization handling and document it
 * Allocate the NcDataSNIACov light-curve covariance only when written
 * Warn at initialization when OMP_NUM_THREADS exceeds OMP_THREAD_LIMIT
 * Passing hard priors to getdist.
 * Compare the data tests' doubles exactly only when bit identical
 * Document the toy data likelihoods and drop their dead code
 * Check NcmCatalog rows and data and report bad rows through GError
 * Give NcmDataDist1d and NcmDataDist2d aborting default methods
 * Expose and serialize the NcmDataPoisson bin edges and drop dead code
 * Pass NcmDataGaussCovMVND covariance updates through set_cov
 * Use the profiled likelihood in the NcmDataGaussDiag w-mean Fisher matrix
 * Keep the NcmDataGaussCov Cholesky factor in step with its covariance
 * Keep the NcmDataGauss Cholesky factor in step with its inverse covariance
 * Prepare every NcmDataset member before the mean and Fisher paths
 * Evaluate the NcmData Fisher covariance and bias at the point itself
 * Serialize long and ulong properties with their value
 * Make two new tests hold on every CI platform
 * Fix the toy models' docs and names and check the MVND dimension
 * Let NcmModelBuilder extend registered models and re-create types safely
 * Make NcmMSetFuncList lookups unambiguous and select bindable
 * Make NcmMSetFunc eval-x a checked copy and guard scalar evaluations
 * Add ncm_mset_func_eval_array and make NcmMSetFunc1 bindable only
 * Annotate NcmMSetFunc numdiff out as inout and fix its leak
 * Make NcmMSetFunc unique names exact and collision free
 * Accept stacked ids in NcmMSet __setitem__ and validate stack positions
 * Keep the free-parameter map across NcmMSet save/load and fix load cleanup
 * Remove ncm_mset_fparam_validate_all and document the free-parameter map
 * Make NcmMSet fit types follow a set fmap and fix two printers
 * Keep the NcmMSet free-parameter map consistent and bound stack positions
 * Remove a host's submodels from NcmMSet and reject unknown namespaces
 * Name the current host in the NcmModel cross-host error
 * Fix NcmModel reparametrization, description and name lookup defects
 * Keep NcmReparam names found after deserialization and reject singular T
 * Fix NcmModelCtrl submodel flags and NcmVParam component setters
 * Accept a NULL model in NcmPowspecFilter and calibrate it at zi
 * Continue NcmPowspecCorr3d into the FFTLog padding and fix its calibration
 * Continue NcmPowspecSpline2d smoothly and abort outside its z range
 * Review NcmPowspec docs and test it against Arb integrals
 * Fix build and test warnings found by the PR #382 lanes
 * Fix the pip and macOS lanes of PR #382
 * Remove a stray LCOV_EXCL_START in NcmSpline2dSpline
 * Test NcmSphereMap in C against frozen healpy truth tables
 * Check NcmSphereMap C_l, pixel access and noise, and test C(theta)
 * Document NcmSphereMap alm2map and test it against healpy at any lmax
 * Document NcmSphereMap map2alm and check its arguments
 * Make NcmSphereMap FITS I/O lossless and readable by and from healpy
 * Compute NcmSphereMap pixel centres without cancellation near the poles
 * Fix NcmSphereMap reordering precision and check its indices
 * Make NcmSphereNN searches safe and document them
 * Validate and document the rectangular sky footprint
 * Place the sbessel_j output grid at the first maximum of j_l
 * Refresh FFTLog docs and tests after the taper
 * Remove NcmFftlogSBesselJLJM and NcmPowspecSphereProj
 * Taper the FFTLog padding ends and enforce the filter's z-knot limit
 * Add a transparent bias to FFTLog and make its grid path-independent
 * Set explicit filter tolerances in the castro and golden tests
 * Calibrate the power spectrum filter with the smooth padding
 * Continue the input smoothly into the FFTLog padding
 * Review sbessel ODE solver docs
 * Review Levin internals; count guard fallbacks once
 * Review Levin public docs
 * Review sbessel integrator and GL docs
 * Document spherical harmonics accuracy near the poles and test it
 * Review sf_sbessel and spherical harmonics docs
 * Review multiprecision specfunc docs; fix leaks and Si thread safety
 * Review integration docs; fix error scaling, leaks and silent failures
 * Use KERNEL_EXACT in SijCalculator
 * Review ode_spline and 2D spline docs; fix silent and stalling failures
 * Review spline docs, remove NcmSplineRBF and 4POINTS, move NcmSplineFuncTest to
     tools
 * Improve power spectrum extrapolation and primordial range handling
 * Fix spline integration, validation, and GSL handling
 * Reorganizing FFTW interface around the new pattern.
 * Using a more robust angle computation.
 * Fix algebra utilities and numerical edge cases
 * Documenting and improving Spectral.
 * Improving documentation and tests for LaurentSeries.
 * Improving docs and tests for Lapack, Matrix and Vector.
 * Review and improved Serialize (including bug fix).
 * Review improving ncm_cfg documentation.
 * Review constants documentation and tests.
 * Moving documentation to the correct directories.
 * Moving NcmComplex away from NcmUtil.
 * Review NcmUtil docs and added tests.
 * Review ObjArray, RNG and Timer.
 * Review and tests for PLN1D, FunctionCache and ISet
 * Review: DTuples, func_eval and MemoryPool.
 * Moving all dev-notes out of the repo.
 * Improving metadata, removed redundant options.
 * Rewrite and parallelize VKDE evaluation paths
 * Cut the VKDE evaluation batch into equal tiles of at most 256 points
 * Regenerate the ncm.pyi stub for ncm_matrix_sub_row_vector
 * Optimize NcmStatsDist fitting and expand coverage
 * Update NcmStatsDist defaults and refactor fitting infrastructure
 * Improve catalog analysis and code organization
 * Refactor autocorrelation handling and Markovian catalog tracking
 * Update APES weighting and covariance controls
 * Extend APES shrinkage support and improve ESMCMC robustness
 * Add APES centre shrinkage and improve sampler robustness
 * docs: add BAO SDSS DR16 fitting example
 * More coverage tests.
 * Updating UltraLevin project. Improving Xcor CLI.
 * xcor cls command (#374)
 * Documenting ultra levin (#373)
 * docs: add H(z) fitting example
 * CI: update conda lock files
 * Skip the last Chebyshev doubling when the level below it predicts failure
 * Use the attached recombination history for the decoupling redshift (#361)
 * Using chebyshev as default (#369)
 * Drive the CCL bridge through NcXcorKernelTable and NcXcorSolver (#368)
 * Close the RSD follow-ups from the post-merge review (#367)
 * Keep the FIXED_NODES redshift marginal non-negative by construction (#301)
 * Close the coverage gaps of the RSD additions
 * Removing gtk-doc comment from header.
 * Silence maybe-uninitialized and duplicate-doc warnings in the Levin integrator
 * Expose dorsd in numcosmo xcor kernel view
 * Add redshift-space distortions to NcXcorKernelGal
 * Add derivative-weighted integrals to NcmSBesselIntegrator
 * Fix the kquad Arb generator applying side a's kdep to side b (#364)
 * Rewrite comments and docs in factual language without dropping content (#363)
 * Update TESTING.md for the lane-based CI layout
 * Silence the _FORTIFY_SOURCE warning in the coverage build
 * Name the two test sets in CI instead of subtracting and re-adding
 * Run the C tests as one coverage shard instead of three
 * Run the C acceptance and statistical tiers optimized, not instrumented
 * Revive dead test code and keep the harness out of the coverage figure
 * Close the reachable xcor gaps in nc_xcor and the kernel component
 * Exercise both xcor closures in the kernel-space quadrature tests
 * Test the NcmMSetTransKernCat CHOOSE sampler on a synthetic catalog
 * Test NcmMSetCatalog on a synthetic multi-chain catalog
 * Test the NcmStatsVec convergence diagnostics on curated chains
 * Correct the l_limber sense in the xcor integrand and kquad tests
 * Gate capability tests on their marker, not on path keywords
 * Move the Planck 2018 CLI tests out of the app gate
 * Split NcmStatsDist tests into mechanics and distribution recovery
 * Split NcmFitESMCMC tests into mechanics and covariance recovery
 * Add C unit tests for the xcor subsystem
 * Handle k ranges degenerate to within rounding in NcXcorKernelComponent
 * Give two g_message lines the '# ...\n' shape
 * Split the CLASS-heavy NcCBE checks into their own executable
 * Give each NcmFit algorithm its own test executable
 * Drop the last references to duration-based slicing
 * Slice the coverage C tests by count and drop the duration machinery
 * Initialize MPI only under a parallel launcher
 * Fail a stuck test instead of hanging its lane
 * Run the acceptance tier optimized instead of instrumented
 * Revert "Distribute the pytest lanes with worksteal instead of pinning files"
 * Block pytest-randomly in the test lanes
 * Distribute the pytest lanes with worksteal instead of pinning files
 * Relax an unreachable tolerance in the xcor view CLI test
 * Remove the last references to the deleted sbessel abstol
 * Remove the caller abstol from the spherical-Bessel integrator
 * Collect coverage from unoptimized code
 * Reject undersized two-fluids initial-condition vectors
 * Stop instrumenting the bundled external libraries
 * docs: use the house array and DataFrame naming in the SNIa+BAO example
 * Make NcmSplineBSpline evaluation thread-safe
 * CI: update conda lock files
 * Add a halo mass function example (#343)
 * Apply the scale dependence on the radial Limber branch
 * Stop the radial Limber branch from counting growth twice
 * Add a kernel for tabulated radial windows
 * Route the last two downloads through the shared fetcher
 * Let a radial kernel carry shear factors
 * Promote the analytic kernel base out of tests
 * Cache the downloaded data files in CI
 * Cover the old-catalog read path, and refuse a reparametrization that cannot be
     rebuilt
 * Fix the Planck baseline download race, and make such failures visible (#344)
 * Let a model resize a reparametrization it outgrew
 * Stop the xcor lane from exhausting memory
 * Updated stale stubs.
 * Put the kernel-space quadratures behind one table
 * Split nc_xcor.c into its three tiers
 * Construction-fixed submodels (#339)
 * Rename HWLCatalogID to WLCatalogID
 * Drop the early dark energy bound and its BBN name
 * Move Yp from the cosmology to NcBBNParametrized
 * Let a loaded mset replace a default submodel
 * Let a serialized file set a write-only property
 * Give every cosmology a nucleosynthesis submodel
 * Add the NcBBN submodel and its PArthENoPE implementation
 * Freeze pre-NcBBN serialization fixtures
 * Make the HSC catalog identifiers reachable
 * Pin every enum's nick prefix explicitly
 * Remove NcHICosmoGCG and NcHICosmoIDEM2
 * Drive the synthetic Planck tests from stored spectra
 * Gate the Planck data tests on their own marker
 * Planck likelihood reimplementation (#332)
 * Split the certified integrals on their oscillation scale
 * Assert C_ell against certified values
 * Pick the conditioned form of j_ell
 * State the k truncation instead of bounding it badly
 * Certify C_ell against Arb, not just the radial integral
 * Measure the outer k-integral against its own reference
 * Keep the smoothed top-hat accurate in its own tails
 * Adding viewer options for different closure modes.
 * * Fix a stale panel coefficient bound in a comment. * Moved the closure-type
     choice from NcXcorKernel to NcXcor. * Added closure_type arguments to the
     three kernel closure entry points. * Updated xcor tests, the accuracy
     tutorial and the xcor viewer for the new API. * Regenerated Python stubs
     (nc.pyi).
 * Make the panel order cap configurable
 * Cap Chebyshev panels at order 5
 * Pool spectral workspaces instead of sharing them per kernel
 * Choose the block integrator inside the block integrator
 * Cover the spectral integration path
 * Integrate spectral closures exactly on their common panel refinement
 * Check the Chebyshev closure directly against Arb
 * Keep the spline closure for Limber multipoles
 * Split the Chebyshev closure into panels
 * Add a Chebyshev representation for the k-space closure
 * Reuse caller-supplied coefficient matrix
 * Expand a set of functions in Chebyshev on a shared grid
 * Report the fit error the closure achieved, not the tolerance it was asked for
     (#326)
 * Separate Tutorials from Examples, and fix the CCL two-point timing comparison
     (#325)
 * Document the exact k-quadrature (#324)
 * Two profiling fixes: the growth special function at z=0, and angular_cl's
     k-grid (#323)
 * Build test_ncm_fit.c once per optimizer group instead of as one binary (#322)
 * Add closed-form references for the xcor stack, checked against Arb (#321)
 * Fix the measure of the sbessel convenience integrands, and group truth tables
     (#320)
 * Do not mark delegation-only prepare() methods as current (#319)
 * Adjusting prepare/prepare_if_needed calls.
 * Improving ultra levin (#317)
 * CI: update conda lock files
 * Analytic xcor kernels, and a floor under scaled-abstol (#315)
 * Fix the two CMB ISW aborts in the kernel-space methods (#297, #298) (#314)
 * * Rename the ell loop variable in compute_kernel's f_ell comprehension, flake8
     E741.
 * * Guard the pyccl import with importorskip, so the file skips instead of
     failing collection where pyccl is absent.
 * Guard batching invariance to within the measured error budget
 * Clamp the CCL bridge block size to what the solver honours
 * Batch the CCL bridge over ell blocks: one factorisation per block,
     ell-dependent scalars reapplied after
 * Add a CCL-facing non-Limber angular_cl backed by the NumCosmo Levin solver
 * * Cover copy_empty: the order is carried over and the copy is independent. *
     Cover deriv_nmax against (order-1)! for x^(order-1), and against
     eval_deriv/eval_deriv2 at orders 2 and 3. * Cover the instance name, its
     rebuild on an order change, and its use in ncm_spline_set's min-size error.
     * Factor the child-process call the abort tests share into _run_child.
 * * Restore NCM_FFTW_* with monkeypatch in test_cfg.py, so the deliberately
     invalid values no longer leak into the rest of the session. * Drop
     NCM_FFTW_* from the child environment in
     test_impossible_request_fails_loudly.
 * * Extract the apt package list into a .github/actions/apt-deps composite
     action, used by both apt jobs. * Group the brew package list and document
     why the Python packages come from brew. * Resolve gmp's prefix with brew
     --prefix instead of a hardcoded Cellar path. * Drop the cfitsio pkg-config
     debug step. * Align the documented apt/brew lists with CI, and note that
     Ubuntu 24.04's GSL is too old.
 * * Require GSL >= 2.8: NcmSplineBSpline uses the rewritten gsl_bspline API. *
     Move the apt build jobs to ubuntu-26.04, which ships GSL 2.8; 24.04 has
     2.7.1. * Use libgsl-dev instead of the transitional libgsl0-dev. * Bound
     gsl in environment.yml and regenerate the conda lock files. * Update the
     documented GSL requirement.
 * Derive the B-spline order from the requested tolerance, and refuse requests the
     samples cannot support
 * Add NcmSplineBSpline: interpolating B-spline of arbitrary order
 * * Recalibrate the construction-tolerance gap assertion from 1e0 to 1e-3. * Add
     a CLI-level guard that --integrator-reltol reaches the computation.
 * xcor kernel view: pass integrator tolerances at construction
 * xcor kernel view: compute and plot angular power spectra
 * sbessel levin: record per-panel contributions for diagnostics
 * sbessel ode solver: restrict the resolution floor to the oscillatory span
 * sbessel ode solver: floor the truncation order at the oscillation count
 * Correct CONTRIBUTING: registration lives in ncm_cfg.c, python tests are
     auto-collected by marker
 * update_pyi.sh: run from its own directory and scope black to the stubs
 * Update data object inspection source path (#279)
 * Ssc cubature fallback (#307)
 * Pin conda toolchains and use committed lock files
 * Updating to python 3.13 on CI.
 * Adding more tests. Calibrating old tests.
 * hmf: add Castro halo mass function and bias models with validation and
     documentation
 * jpas_forecast24: cover the varying-Sij paths
 * Update stubs and fix a mypy error
 * * Add `NcXcorSSCSij` as the native NumCosmo SSC (S_{ij}) calculator, replacing
     the Python implementation. * Support full-sky and arbitrary-mask
     calculations using the mask (C_\ell). * Add (f_{\rm sky}) area rescaling
     through the new `area` property. * Add cosmology-dependent (S_{ij})
     recalculation for `NcmDataClusterNCountsGauss`. * Keep the resampling
     matrix fixed while allowing (S_{ij}) to vary during fitting. * Add
     `--vary-fitting-sij` to `jpas_forecast24`. * Make `NcXcorSSCSij` safe to
     default-construct and validate completeness in `prepare()`. * Fix
     `set_ssc_sij()` so deserialized data remain initialized. * Add the required
     top-hat kernel constructor with an integrator. * Add tests against the
     Python reference, including full-sky, partial-sky, serialization, wiring,
     and construction cases. * Add SSC theory documentation and references.
 * ssc: replace PySSC with NumCosmo implementation (#304)
 * Adding fsky and tests.
 * Adding area dependency to Sij.
 * xcor dev-notes: record the ISW Limber step analysis behind GH #297
 * xcor: register NcXcorKernelCMBISW and export the two missing xcor headers
 * xcor: export NcXcorSolver from the umbrella header, and cover the patch's
     untested lines
 * xcor: expose precision knobs in the CLI, and add the design notes and CCL
     benchmark
 * xcor: add NcXcorSolver, batched angular cross-spectra with per-block Levin
     integrator reuse
 * ncm_cfg: cache FFTW wisdom file I/O across repeated load/save calls
 * ncm/nc: give adaptive routines an absolute error scale, and stop evaluating out
     of range
 * Testing new additions.
 * Exposing resample type in CLI for WL.
 * Testing conditional addition of derived for PopBeta.
 * Making Beta pop derived parameters conditional on fitting population.
 * Fixing catalog HDU0 dumping.
 * Adding backwards compatibility with 'std_shape'.
 * Making the overwrite of a RNG state an error on initialized catalogs.
 * Fixing empty ObjectArray crashing experiments.
 * Improving error message on deserializing yaml.
 * * Fixed test_run_mcmc_apes_plot_corner_too_many_plot_names: it asserted on a
     substring straddling "--plot-name", but Rich highlights option-looking
     tokens and injects ANSI codes between their characters when color is forced
     (as CI does for Typer's error panel, unlike a plain local terminal),
     breaking the naive substring match. Assert on a dash-free portion of the
     same message instead.
 * * Closed the two files' patch-coverage gaps flagged by Codecov (loading.py 19
     missing -> 0, catalog.py 8 missing -> 0), verified by intersecting
     coverage.json against the PR diff. * Removed dead-code guards in PlotCorner
     and DerivedQuantityError: Click already enforces
     mcmc_file/--variable/--expr as required, so the manual "at least one"
     checks could never fire. * Added tests for the real remaining gaps:
     --include/--exclude column filtering (all three branches), an empty-catalog
     burnin, --plot-name count mismatch, --mark-bestfit, catalog visual-hw and
     param-evolution (previously untested commands), a single-chain (run mc)
     catalog, a missing-catalog-file error, and load_catalog()'s direct
     negative-tail guard.
 * * Fixed 47 mypy errors: LoadCatalog now declares its load_catalog()-derived
     attributes as typed dataclass fields (mcat, mset, functions, etc.) instead
     of injecting them via self.__dict__.update(), which was invisible to the
     type checker. * Widened mcat_to_catalog_data's indices parameter to also
     accept a plain list[int], matching what it already accepted at runtime. *
     mypy --exclude '.*meson.*|numcosmo_py/generate_stubs\.py' -p numcosmo_py
     now passes clean.
 * Adding tests and fixing bugs.
 * * Catalogs are now self-sufficient: NcmMSetCatalog embeds the model-set and, if
     used, the functions array in FITS HDU0 as a versioned NcmVarDict, no
     experiment file needed to read one back. * NcmVarDict gains typed
     object/object-array set/get accessors, taking an explicit NcmSerialize
     argument. * NcmFitESMCMC/NcmFitMC embed the functions array into the
     catalog instead of writing the never-read .oa sidecar. * Added
     ncm_mset_catalog_peek_info_from_file for cheap nrows/nchains lookup without
     a full catalog load. * Added g_assert(NCM_IS_MSET) guards after HDU0
     deserialization and clarified the burnin-exceeds-catalog error message. *
     numcosmo catalog: LoadCatalog no longer requires an experiment file;
     check-m2lnl is the sole exception since it needs a live likelihood. *
     numcosmo catalog plot-corner: now takes multiple catalogs positionally,
     overlaid in one plot, with --plot-name for legend labels; drops
     --extra-experiment/--extra-mcmc-file/--extra-burnin. * numcosmo catalog:
     --burnin now means iterations (ensemble steps) instead of raw rows; added
     --tail to keep only the last N iterations. * numcosmo catalog: user-input
     errors (bad burnin/tail, missing catalog, incompatible experiment, etc.)
     now raise typer.BadParameter for a clean CLI message instead of a
     traceback. * Removed tools/mcat_plot_corner, fully superseded by catalog
     plot-corner. * Added C tests for the new NcmVarDict accessors and
     NcmMSetCatalog HDU0/functions-array round trips; updated Python CLI tests
     for the new signatures. * Regenerated ncm.pyi and nc.pyi stubs.
 * Reintroducing probability floor.
 * Adding more tests.
 * Improving tests, accuracy and quadring against negative prob.
 * Increased outer region number of knots.
 * FixedQuad: origin-divergence-aware marginal sum, cheaper tail panel, numerical
     fixes.
 * Increasing default mass upper bound for WL analysis.
 * Ignoring fyaml compilation warnings.
 * Including _GNU_SOURCE in fyaml compilation.
 * FixedQuad: fix alpha<2 Beta population divergence, plus rotation-covariance and
     caching fixes (#292)
 * docs: add SNIa+BAO confidence region example
 * Updates and fixes to allow use of modern C standards gnu or strict C.
 * Fixed mypy.
 * Fixing test.
 * Adding support for multiple expressions.
 * Adding support for derived parameters in catalog analyze. Adding unit tests for
     new features.
 * Fixing tests.
 * Specialized kernel for fixed quad WL computation.
 * Fixing lto related warnings.
 * Reworking beta distribution for ellipticity. Now we model the distribution of
     |chi| or |e| instead of |chi|^2 or |e|^2.
 * Removed guard.
 * * Better calibrating initial beta distribution. * Fixing angular border problem
     with read data analysis.
 * Simplified comments and docs.
 * Improving tests.
 * Testing better testing duration cache/restore.
 * Uncrustify.
 * * Fix ESMCMC OpenMP initialization deadlock by replacing ordered retries with
     serial redraw / parallel evaluation rounds, adding bounded retries and
     deterministic RNG ordering. * Fix MPI ESMCMC initialization to redraw
     walkers in place, preserve walker indices, eliminate incorrect compaction,
     and bound retries. * Reset ESMCMC walker acceptance flags before
     initialization to avoid stale state across runs. * Add `max-iter` limit to
     the Gaussian transition kernel to prevent infinite retries when proposals
     remain out of bounds. * Add deterministic parity tests comparing serial,
     OpenMP, and MPI initialization results. * Detect the number of OpenMP
     threads at test runtime instead of Meson configure time, and update the
     testing infrastructure and documentation. * Add support for loading curated
     HSC weak-lensing catalogs by catalog ID with automatic download and
     caching. * Add catalog-wide metadata support to `NcmCatalog` with
     serialization. * Extend the cluster weak-lensing application to load real
     catalogs, validate catalog coverage, and support metadata defaults. *
     Refactor the cluster weak-lensing CLI to separate mock-data generation from
     real-data loading. * Add test coverage for real catalog loading and
     regenerate Python stub files. * Improve the CLI interface for loading
     catalogs.
 * * build_check.yml: add TEST_DURATION_CACHE_VERSION to force-invalidate the
     cached per-test duration files used by the C-suite slicer, bumped to 1 to
     rule out a stale/corrupted cache as the cause of ncm_mset_catalog being
     dropped from all three coverage shards on this PR. * build_check.yml: print
     the computed slice-tests.txt contents in the "Compute test slice" step for
     visibility into what each shard actually selects.
 * * test_ncm_mset_catalog.c: add C tests for HDU0 mset round-trip, legacy .mset
     sidecar fallback, and the fatal error when neither is present.
 * * ncm_mset_catalog: embed the mset as GVariant binary in the FITS primary HDU
     (HDU0) instead of a separate .mset GKeyFile sidecar, written once at file
     creation. * ncm_mset_catalog: keep read-only support for the legacy .mset
     sidecar for old catalog files. * Added "numcosmo catalog dump-mset" CLI
     command to export a catalog's mset as YAML.
 * Update requirements
 * * build_check.yml: don't cache an empty {} duration extract. Observed for real
     on this PR's own first (failed) run: "Save test slice durations" runs with
     if: always(), so a run that fails before any test executes (testlog.json
     never created) still caches an empty durations fragment -- restore-keys
     prefix matching picks the most recent entry regardless of content, so the
     *next* run silently inherited zero real duration data and fell back to a
     count-balanced split. Guard with a has_data check so only non-empty
     extracts get cached.
 * * test_slicer.py: fix --durations argparse config -- it was action="append"
     (expects the flag repeated once per file) but the workflow calls it as one
     --durations flag followed by multiple space-separated filenames, which
     needs nargs="+". Caused the first real CI run's "Compute test slice" step
     to fail outright ("unrecognized arguments"). Re-verified against fresh
     complete local duration data with the exact multi-file invocation the
     workflow uses (252.0s/252.9s/252.8s balance, all 76 tests accounted for)
     before pushing.
 * * test_slicer.py: add a summary subcommand and factor the testlog.json
     JSON-lines parsing shared with extract into _read_testlog(), replacing the
     coverage job's inline Python heredoc ("Test timing summary" step) with a
     real, locally-runnable script call.
 * * .github/actions/setup-miniforge: extract the ~30-line Setup miniforge / Cache
     Conda env / Update environment / Save conda-forge cache sequence
     (previously duplicated between build-miniforge and
     build-miniforge-coverage) into a shared composite action, parameterized by
     python-version and optional mpi. * .github/scripts/test_slicer.py: new
     script replacing meson's count-balanced --slice K/N for the coverage job's
     C-tier shards with a duration-aware greedy longest-processing-time-first
     bin-pack (plan subcommand), fed by historical per-test durations extracted
     from testlog.json (extract subcommand). Falls back to a
     uniform/count-balanced split when there's no history yet. *
     build_check.yml: wire the C-tier (c-1/c-2/c-3) legs of
     build-miniforge-coverage to compute their test list via test_slicer.py
     instead of --slice K/N, caching each slice's fresh durations (GitHub
     Actions cache, per-slice keys to avoid the 3 concurrent legs racing) for
     the next run to read back.
 * * build_check.yml: also set UCX_TLS=tcp,self,sm -- OMPI_MCA_pml=ob1 only fixes
     OpenMPI's own PML selection, not MPICH's ch4:ucx netmod (which hit the same
     underlying mana-NIC issue with a different error: MPIDI_UCX_init_worker
     "Address not valid"). UCX_TLS is read by UCX itself under either MPI
     implementation, so this is the fix that actually covers mpich.
 * * build_check.yml: set OMPI_MCA_pml=ob1 globally, bypassing OpenMPI's UCX
     transport -- GitHub-hosted runners' paravirtualized "mana" NIC advertises
     IB-like RDMA verbs it doesn't actually support, causing spurious "Failed to
     create UCP worker" flakes on MPI test jobs (real regardless of code
     changes, e.g. PR #283's ncm_fit_esmcmc_mpi ERROR with exit status 0 and all
     TAP subtests ok).
 * * --bench report: prefix the printed line with '# ' to match NumCosmo's usual
     comment-style log/report output convention.
 * * Add a --bench flag (RunCommonOptions in run_fit.py) reporting wall-clock time
     and peak RSS at the end of a run; covers run fit/test/mc/mcmc since they
     all share this base class. * run test: also call end_experiment() at the
     end (was a pre-existing gap -- --output/--log-file were silently ignored,
     and --bench had nothing to hook into). * Add test_run_bench covering run
     test and run fit.
 * * test_fit_mc.py: add a use_threads property/getter round-trip test (was only
     exercised via set_use_threads(), not the GObject property system or
     get_use_threads()). * test_ncm_fit_esmcmc.c: extend
     test_ncm_fit_esmcmc_properties() with the same use-threads property/getter
     round-trip coverage, mirroring the existing skip-check/log-time-interval
     pattern.
 * * tests/python/meson.build: exclude omp-marked tests from the xdist fast lane
     (were running under OMP_NUM_THREADS=1, so when 2-3 landed on one worker it
     ran alone for the tail while the rest of the pool idled); they now run only
     in the dedicated single-process, real-OMP-threads pytest-omp lane. *
     build_check.yml: add a py-omp shard to the coverage job's matrix, since it
     was only getting omp-marked test coverage incidentally (under OMP=1) via
     the plain python suite -- excluding them from that suite would have
     silently dropped coverage.
 * * NcmFitMC/NcmFitESMCMC: replace the dead nthreads guint property with a
     use_threads gboolean (set_use_threads/get_use_threads), matching the
     existing NcmStatsDist/NcmFitESMCMCWalkerAPES convention; real thread count
     remains OMP_NUM_THREADS-driven. * NcmFitMCMC: delete nthreads and its
     entire multi-threaded path outright -- it was never implemented
     (g_assert_not_reached in _ncm_fit_mcmc_mt_eval) and would abort if ever
     triggered. * NcmFitESMCMC: move the odd-nwalkers validation out of
     set_nthreads into constructed() (it's a Stretch-move requirement, not a
     threading concern); start-run log now probes the real OpenMP thread count
     live instead of echoing a stored value. * NcmFitMCBS: ncm_fit_mcbs_run()'s
     bsmt param changed guint -> gboolean. * darkenergy: --mc-nthreads int flag
     replaced with --mc-use-threads bool flag. * Updated all CLI, sampling
     helper, experiment, and example call sites to the new API. * Updated C and
     Python tests for the new API; rewrote test_fit_mc_keep_order.py's
     dead-wiring-bug documentation to describe the new design instead. *
     Regenerated numcosmo_py/ncm.pyi.
 * * Add NcDataClusterWLFactor's register_shared override, anchoring its obs
     catalog (including per-galaxy pz splines for the Spline redshift scheme)
     instead of deep-copying it per parallel worker. * Move the register_shared
     call from constructed() into start_run(), so a dataset swapped in after
     construction (e.g. via set_obs()) is still anchored correctly. * Reset the
     NcmSerialize instance fully in end_run(), releasing shared anchors so a
     long-lived NcmFitMC/NcmFitESMCMC reused across many runs doesn't keep stale
     data alive between them. * Add register_shared regression tests (positive
     and negative control) to test_data_cluster_wl_factor.py.
 * * Add NcmData::register_shared vfunc (default no-op) letting a data subclass
     register its own large read-only internals as NcmSerialize anchors. * Add
     ncm_dataset_register_shared() to call it across every NcmData in a dataset.
     * Wire it into NcmFitMC/NcmFitESMCMC's constructed(), using each object's
     own internal NcmSerialize, before any per-worker duplication happens. *
     Regenerate ncm.pyi stubs.
 * * Port docs/tutorials/python/cluster_wl_simul.qmd off the deleted legacy
     NcGalaxySD*/NcDataClusterWL classes to the Factor pipeline (fixes the BDocs
     quarto render failure). * Update docs/theory/wl_ellipticity.qmd and
     galaxy_wl_framework.qmd to drop dangling legacy gtk-doc cross-references
     and stale "still being built" status text. * Drop remaining dangling
     NcGalaxySD*/NcDataClusterWL comments across nc_data_cluster_wl_factor.c,
     nc_galaxy_shape_pop.{c,h}, ncm_laurent_series.{c,h}, and their tests.
 * * Regenerate nc.pyi stubs to drop the removed legacy classes.
 * * Drop dangling "matches/direct translation of legacy" comments in
     nc_data_cluster_wl_factor.c and nc_galaxy_shape_factor_var_add.c now that
     legacy is gone.
 * * Rewrite tests/python/numcosmo_py/experiments/test_wl_app.py against the
     Factor-based cluster-wl CLI schema. * Port fixtures_xcor.py's LSST bin
     fixtures to GalaxyRedshiftBinning.lsst_srd_edges/compute_dndz.
 * * Rewrite numcosmo_py/experiments/cluster_wl.py to build the NcGalaxy*Factor
     pipeline instead of the legacy NcGalaxySD* classes. * Update
     numcosmo_py/app/generate.py cluster-wl CLI flags to the new
     z_dist/shape_dist schema. * Port examples/example_wl_likelihood.py to the
     Factor API. * Migrate xcor kernels.py/view.py from
     GalaxySDTrueRedshiftLSSTSRDType/new_lsst_srd_bins to
     GalaxyRedshiftPopLSSTSRDType/GalaxyRedshiftBinning.
 * * Delete the legacy NcGalaxySD*/NcDataClusterWL galaxy WL pipeline (24 sources)
     superseded by the NcGalaxy*Factor/NcDataClusterWLFactor pipeline. * Remove
     legacy entries from numcosmo/meson.build and tests/c/meson.build, and
     legacy #include/registration in ncm_cfg.c and numcosmo.h. * Delete
     legacy-only C and Python test files and legacy ref/unref coverage in
     test_ncm_generic.c.
 * * Relocate NcDataClusterWLResampleFlag/IntegMethod enums out of the legacy
     header into nc_data_cluster_wl_factor.h.
 * Make faulthandler dump per-test tracebacks repeatedly instead of once.
 * Pack remaining galaxy WL frozen-fixture dicts into truth_tables/wl binfiles.
 * Packing reference data into bin files. Calibrating tests.
 * More tolerance loosing.
 * Black.
 * Loosen bit-exact frozen-fixture tolerances that route through libm-sensitive
     code
 * Calibrating tests for the new WL framework.
 * Galaxy WL calculators (#277)
 * Cluster richness poisson lognormal (#276)
 * Refactor wl tools (#275)
 * Add build directory and info files to .gitignore
 * Add optimizations for weak lensing calculations (#269)
 * Add data-driven local curvature priors for w(z)/q(z) reconstruction (#271)
 * Add selectable SNIa resample strategy (#272)
 * Add selectable knot placement to spline reconstruction models (#270)
 * Curvature-prior w(z)/q(z) reconstruction toolkit (#268)
 * Docs overhaul (#267)
 * Directory restructure (#266)
 * Chore/mechanical improvements (#265)
 * Fix .mset save with sub_fit writing sub_fit params instead of main fit params
     (fixes #23)
 * Add support for fixed node integration in NcDataClusterWL (#252)
 * Organize test directories (#264)
 * Improving tests (#263)
 * Matching by ID  (#261)
 * nc_multiplicity_func_bhattacharya: add convention enum selecting the a(z)
     redshift evolution (Bhattacharya 2011 vs Heitmann 2019). Add new_full
     constructor and convention get/set; default keeps Bhattacharya 2011. Add C
     and Python tests covering both conventions and serialization. Regenerated
     Python stubs (nc.pyi).
 * nc_galaxy_sd_shape: HSM shape-measurement models and direct shear estimators.
 * Bt mass function (#260)
 * jpas_forecast: make photo-z scatter sigma0 configurable (#259)
 * Add data object inspection documentation example (#258)
 * Update files for 0.27.0 release
 * Release version 0.27.0
 * Adding codecov configuration.
 * Avoid create commit status for fork PRs.
 * Fixing minor issues with unused variables. (#257)
 * Full site upload from GHA. (#256)
 * Improving documentaton build time and ipynb generation. (#255)
 * More verbose run in RTD.
 * Settling for requests.
 * Reverted to using requests.
 * Trying token instead of Bearer.
 * Handling RTD stable and latest relases.
 * Build artifacts on GHA and uploading to RTD. (#254)
 * Csq1d phase support (#253)
 * Moving wspline example to quarto. (#251)
 * Update min python version to 3.11 (#250)
 * New inspect command. Removing tap support in pytest. (#249)
 * Calibrating pln1d test. (#248)
 * Downloadable docs (#247)
 * Fix xcor benchmark document. (#246)
 * Update firecrown import (#245)
 * Improve tests to avoid leaving leftovers. (#244)
 * Improving richness analysis tools (#241)
 * Extend spectral (#243)
 * Update install guide (#239)
 * New xcor (#240)
 * docs: Replace FIXME placeholders with proper documentation (73/588) (#238)
 * docs: Replace FIXME placeholders with proper documentation (170/758) (#237)
 * Mass concentration bhattacharya13 (#234)
 * Jpas forecast bug (#235)
 * Create notebook and minor fix in the documentation. (#196)
 * Twofluid tensor (#236)
 * Jpas forecast (#170)
 * Testing Richness proxies (#117)
 * Removed deprecated option. (#233)
 * New version 0.26.0
 * Splitting NumCosmo initialization. (#232)
 * Bounce tutorial (#212)
 * DE w(z) spline - experiment (#226)
 * Add option to choose non logarithmic integrand (#230)
 * Restricting NumCosmo version and trying texlive-core.
 * Trying r-tinytex.
 * Fixed version test.
 * New release 0.25.0.
 * Environment for NumCosmo use.
 * Adding support for derived quantities in MC runs. (#231)
 * Updated stubs.
 * Adding support for random walk in APES. (#229)
 * Updated stubs.
 * Fixing typos.
 * Update stub formatting to match linter conventions. Fixed typos.
 * Updated stubs.
 * Fixing memory leaks. (#228)
 * Improving galaxy integration (#227)
 * Improving stub generation. (#225)
 * Stub update (#224)
 * Refactor galaxy redshift limit functions (#222)
 * Adding support for setting seed to a MC. (#221)
 * Updated stubs.
 * New P(z) redshift object (#190)
 * New DES Y5 SNIA support. (#215)
 * Adding google analytics to NumCosmo site. (#220)
 * Improving reltol for numerical int for fftlog test. (#219)
 * Update enums and reqs. (#218)
 * Adding lintegrate (#217)
 * BAO data - DESI  DR2 2025. (#216)
 * * Included new CC data in the enumerator (nc_data_hubble). (#214)
 * Updated stubs.
 * Included new Cosmic Chronometers H(z) data: (#209)
 * Bao data desi (#213)
 * Updating sundials to 7.2.1. (#211)
 * Fixing compilation glitches. (#210)
 * Updated README links.
 * Better calibration for two-point limber. (#208)
 * Benchmark two-point calculations.  (#207)
 * Benchmark PowerSpectra (#206)
 * Removing outdate/unsupported old code. (#205)
 * Benchmarks (#204)
 * Matching algorithm (#203)
 * Mass concentration duffy08 (#202)
 * Moving examples to docs (#201)
 * Bumping to new version in development.
 * Updating version in pyproject.toml
 * Bumping minor version. (#200)
 * Improving docs (#199)
 * Imported Despali Mass function from jpas-forecast. (#198)
 * Moving to a quarto generated documentation (#197)
 * Coverage python (#185)
 * Matching algorithm (#191)
 * Update stubs (#195)
 * Update spline func tests (#194)
 * Fix docs (#193)
 * Mass concentration klypin11 (#192)
 * Update weak lensing framework (#189)
 * (Re)updating changelog.
 * Better handling of git hash.
 * Updated changelog.
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
