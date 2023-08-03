/** @file perturbations.c Documented perturbation module
 *
 * Julien Lesgourgues, 23.09.2010
 *
 * Deals with the perturbation evolution.
 * This module has two purposes:
 *
 * - at the beginning; to initialize the perturbations, i.e. to
 * integrate the perturbation equations, and store temporarily the terms
 * contributing to the source functions as a function of conformal
 * time. Then, to perform a few manipulations of these terms in order to
 * infer the actual source functions \f$ S^{X} (k, \tau) \f$, and to
 * store them as a function of conformal time inside an interpolation
 * table.
 *
 * - at any time in the code; to evaluate the source functions at a
 * given conformal time (by interpolating within the interpolation
 * table).
 *
 * Hence the following functions can be called from other modules:
 *
 * -# perturb_init() at the beginning (but after background_init() and thermodynamics_init())
 * -# perturb_sources_at_tau() at any later time
 * -# perturb_free() at the end, when no more calls to perturb_sources_at_tau() are needed
 */

#include "perturbations.h"
#include <numcosmo/numcosmo.h>

/**
 * Source function \f$ S^{X} (k, \tau) \f$ at a given conformal time tau.
 *
 * Evaluate source functions at given conformal time tau by reading
 * the pre-computed table and interpolating.
 *
 * @param ppt        Input: pointer to perturbation structure containing interpolation tables
 * @param index_md   Input: index of requested mode
 * @param index_ic   Input: index of requested initial condition
 * @param index_type Input: index of requested source function type
 * @param tau        Input: any value of conformal time
 * @param psource    Output: vector (already allocated) of source function as a function of k
 * @return the error status
 */

int
perturb_sources_at_tau (
  struct perturbs *ppt,
  int             index_md,
  int             index_ic,
  int             index_type,
  double          tau,
  double          *psource
                       )
{
  /** Summary: */

  /** - interpolate in pre-computed table contained in ppt */
  class_call (array_interpolate_two_bis (ppt->tau_sampling,
                                         1,
                                         0,
                                         ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index_type],
                                         ppt->k_size[index_md],
                                         ppt->tau_size,
                                         tau,
                                         psource,
                                         ppt->k_size[index_md],
                                         ppt->error_message),
              ppt->error_message,
              ppt->error_message);


  return _SUCCESS_;
}

/**
 * Initialize the perturbs structure, and in particular the table of source functions.
 *
 * Main steps:
 *
 * - given the values of the flags describing which kind of
 *   perturbations should be considered (modes: scalar/vector/tensor,
 *   initial conditions, type of source functions needed...),
 *   initialize indices and wavenumber list
 *
 * - define the time sampling for the output source functions
 *
 * - for each mode (scalar/vector/tensor): initialize the indices of
 *   relevant perturbations, integrate the differential system,
 *   compute and store the source functions.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input: pointer to thermodynamics structure
 * @param ppt Output: Initialized perturbation structure
 * @return the error status
 */

int
perturb_init (
  struct precision        *ppr,
  const struct background *pba,
  struct thermo           *pth,
  struct perturbs         *ppt
             )
{
  /** Summary: */

  /** - define local variables */
  int sgnK;
  double K = nc_hicosmo_curvature_K (pba->cosmo, &sgnK);

  /* running index for modes */
  int index_md;
  /* running index for initial conditions */
  int index_ic;
  /* running index for wavenumbers */
  int index_k;
  /* pointer to one struct perturb_workspace per thread (one if no openmp) */
  struct perturb_workspace **pppw;
  /* background quantities */
  double w_fld_ini, w_fld_0, dw_over_da_fld, integral_fld;
  /* number of threads (always one if no openmp) */
  int number_of_threads = 1;
  /* index of the thread (always 0 if no openmp) */
  int thread = 0;

  /* This code can be optionally compiled with the openmp option for parallel computation.
   *  Inside parallel regions, the use of the command "return" is forbidden.
   *  For error management, instead of "return _FAILURE_", we will set the variable below
   *  to "abort = _TRUE_". This will lead to a "return _FAILURE_" just after leaving the
   *  parallel region. */
  int abort;

  /* unsigned integer that will be set to the size of the workspace */
  size_t sz;

#ifdef _OPENMP
  /* instrumentation times */
  double tstart, tstop, tspent;
#endif

  /** - perform preliminary checks */

  if (ppt->has_perturbations == _FALSE_)
  {
    if (ppt->perturbations_verbose > 0)
      printf ("No sources requested. Perturbation module skipped.\n");

    return _SUCCESS_;
  }
  else
  {
    if (ppt->perturbations_verbose > 0)
      printf ("Computing sources\n");
  }

  class_test ((ppt->gauge == synchronous) && (pba->has_cdm == _FALSE_),
              ppt->error_message,
              "In the synchronous gauge, it is not self-consistent to assume no CDM: the later is used to define the initial timelike hypersurface. You can either add a negligible amount of CDM or switch to newtonian gauge");

  class_test ((ppr->tight_coupling_approximation < first_order_MB) ||
              (ppr->tight_coupling_approximation > compromise_CLASS),
              ppt->error_message,
              "your tight_coupling_approximation is set to %d, out of range defined in perturbations.h", ppr->tight_coupling_approximation);

  class_test ((ppr->radiation_streaming_approximation < rsa_null) ||
              (ppr->radiation_streaming_approximation > rsa_none),
              ppt->error_message,
              "your radiation_streaming_approximation is set to %d, out of range defined in perturbations.h", ppr->radiation_streaming_approximation);

  if (pba->has_ur == _TRUE_)
    class_test ((ppr->ur_fluid_approximation < ufa_mb) ||
                (ppr->ur_fluid_approximation > ufa_none),
                ppt->error_message,
                "your ur_fluid_approximation is set to %d, out of range defined in perturbations.h", ppr->ur_fluid_approximation);

  if (pba->has_ncdm == _TRUE_)
  {
    class_test ((ppr->ncdm_fluid_approximation < ncdmfa_mb) ||
                (ppr->ncdm_fluid_approximation > ncdmfa_none),
                ppt->error_message,
                "your ncdm_fluid_approximation is set to %d, out of range defined in perturbations.h", ppr->ncdm_fluid_approximation);

    if (ppt->has_nc_density == _TRUE_)
      if (ppt->perturbations_verbose > 0)
        fprintf (stdout, " -> [WARNING:] You request the number count Cl's in presence of non-cold dark matter.\n    Like in all previous CLASS and CLASSgal versions, this will be inferred from the total matter density,\n    but it could make much more sense physically to compute it from the CDM+baryon density only.\n    To get the latter behavior you would just need to change one line in transfer.c:\n    search there for a comment starting with 'use here delta_cb'\n");
  }

  if (pba->has_fld == _TRUE_)
  {
    /* check values of w_fld at initial time and today */
    class_call (background_w_fld (pba,     0.,   &w_fld_ini, &dw_over_da_fld, &integral_fld), pba->error_message, ppt->error_message);
    class_call (background_w_fld (pba, pba->a_today, &w_fld_0, &dw_over_da_fld, &integral_fld), pba->error_message, ppt->error_message);

    class_test (w_fld_ini >= 0.,
                ppt->error_message,
                "The fluid is meant to be negligible at early time, and unimportant for defining the initial conditions of other species. You are using parameters for which this assumption may break down, since at early times you have w_fld(a--->0) = %e >= 0", w_fld_ini);

    if (pba->use_ppf == _FALSE_)
    {
      class_test ((w_fld_ini + 1.0) * (w_fld_0 + 1.0) <= 0.0,
                  ppt->error_message,
                  "w crosses -1 between the infinite past and today, and this would lead to divergent perturbation equations for the fluid perturbations. Try to switch to PPF scheme: use_ppf = yes");

      /* the next check is meaningful at least for w(a) = w0 + wa*(1-a/a0); for general formulas and with use_ppf=no, you may prefer to comment it out... */
      class_test ((w_fld_0 == -1.) && (dw_over_da_fld == 0.),
                  ppt->error_message,
                  "Your choice of a fluid with (w0,wa)=(-1,0) is not valid due to instabilities in the unphysical perturbations of such a fluid. Try instead with a plain cosmological constant or with PPF scheme: use_ppf = yes");
    }
  }

  class_test (ppt->has_vectors == _TRUE_,
              ppt->error_message,
              "Vectors not coded yet");

  if ((ppt->has_niv == _TRUE_) && (ppt->perturbations_verbose > 0))
    printf ("Warning: the niv initial conditions in CLASS (and also in CAMB) should still be double-checked: if you want to do it and send feedback, you are welcome!\n");

  if (ppt->has_tensors == _TRUE_)
  {
    ppt->evolve_tensor_ur   = _FALSE_;
    ppt->evolve_tensor_ncdm = _FALSE_;

    switch (ppt->tensor_method)
    {
      case (tm_photons_only):
        break;

      case (tm_massless_approximation):

        if ((pba->has_ur == _TRUE_) || (pba->has_ncdm == _TRUE_))
          ppt->evolve_tensor_ur = _TRUE_;

        break;

      case (tm_exact):

        if (pba->has_ur == _TRUE_)
          ppt->evolve_tensor_ur = _TRUE_;

        if (pba->has_ncdm == _TRUE_)
          ppt->evolve_tensor_ncdm = _TRUE_;

        break;
    }
  }

  /** - initialize all indices and lists in perturbs structure using perturb_indices_of_perturbs() */

  class_call (perturb_indices_of_perturbs (ppr,
                                           pba,
                                           pth,
                                           ppt),
              ppt->error_message,
              ppt->error_message);


  if (ppt->z_max_pk > pth->z_rec)
  {
    class_test (ppt->has_cmb == _TRUE_,
                ppt->error_message,
                "You requested a very high z_pk=%e, higher than z_rec=%e. This works very well when you don't ask for a calculation of the CMB source function(s). Remove any CMB from your output and try e.g. with 'output=mTk' or 'output=mTk,vTk'",
                ppt->z_max_pk,
                pth->z_rec);

    class_test (ppt->has_source_delta_m == _TRUE_,
                ppt->error_message,
                "You requested a very high z_pk=%e, higher than z_rec=%e. This works very well when you ask only transfer functions, e.g. with 'output=mTk' or 'output=mTk,vTk'. But if you need the total matter (e.g. with 'mPk', 'dCl', etc.) there is an issue with the calculation of delta_m at very early times. By default, delta_m is a gauge-invariant variable (the density fluctuation in comoving gauge) and this quantity is hard to get accurately at very early times. The solution is to define delta_m as the density fluctuation in the current gauge, synchronous or newtonian. For the moment this must be done manually by commenting the line 'ppw->delta_m += 3. *ppw->pvecbaack[pba->index_bg_a]*ppw->pvecbaack[pba->index_bg_H] * ppw->theta_m/k2;' in perturb_sources(). In the future there will be an option for doing it in an easier way.",
                ppt->z_max_pk,
                pth->z_rec);
  }



  /** - define the common time sampling for all sources using
   *   perturb_timesampling_for_sources() */

  class_call (perturb_timesampling_for_sources (ppr,
                                                pba,
                                                pth,
                                                ppt),
              ppt->error_message,
              ppt->error_message);

  /** - if we want to store perturbations, write titles and allocate storage */
  class_call (perturb_prepare_output (pba, ppt),
              ppt->error_message,
              ppt->error_message);


  /** - create an array of workspaces in multi-thread case */

#ifdef _OPENMP

#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads ();
  }
#endif

  class_alloc (pppw, number_of_threads * sizeof (struct perturb_workspace *), ppt->error_message);

  /** - loop over modes (scalar, tensors, etc). For each mode: */

  for (index_md = 0; index_md < ppt->md_size; index_md++)
  {
    if (ppt->perturbations_verbose > 1)
      printf ("Evolving mode %d/%d\n", index_md + 1, ppt->md_size);

    abort = _FALSE_;

    sz = sizeof (struct perturb_workspace);

#pragma omp parallel                                             \
    shared(pppw,ppr,pba,pth,ppt,index_md,abort,number_of_threads)  \
    private(thread)                                                \
    num_threads(number_of_threads)

    {
#ifdef _OPENMP
      thread = omp_get_thread_num ();
#endif

      /** - --> (a) create a workspace (one per thread in multi-thread case) */

      class_alloc_parallel (pppw[thread], sz, ppt->error_message);

      /** - --> (b) initialize indices of vectors of perturbations with perturb_indices_of_current_vectors() */

      class_call_parallel (perturb_workspace_init (ppr,
                                                   pba,
                                                   pth,
                                                   ppt,
                                                   index_md,
                                                   pppw[thread]),
                           ppt->error_message,
                           ppt->error_message);
    } /* end of parallel region */

    if (abort == _TRUE_)
      return _FAILURE_;

    /** - --> (c) loop over initial conditions and wavenumbers; for each of them, evolve perturbations and compute source functions with perturb_solve() */

    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++)
    {
      if (ppt->perturbations_verbose > 1)
      {
        printf ("Evolving ic %d/%d\n", index_ic + 1, ppt->ic_size[index_md]);
        printf ("evolving %d wavenumbers\n", ppt->k_size[index_md]);
      }

      abort = _FALSE_;

#pragma omp parallel                                                    \
      shared(pppw,ppr,pba,pth,ppt,index_md,index_ic,abort,number_of_threads) \
      private(index_k,thread,tstart,tstop,tspent)                           \
      num_threads(number_of_threads)

      {
#ifdef _OPENMP
        thread = omp_get_thread_num ();
        tspent = 0.;
#endif

#pragma omp for schedule (dynamic)

        /* integrating backwards is slightly more optimal for parallel runs */
        /*for (index_k = 0; index_k < ppt->k_size; index_k++) { */
        for (index_k = ppt->k_size[index_md] - 1; index_k >= 0; index_k--)
        {
          if ((ppt->perturbations_verbose > 2) && (abort == _FALSE_))
          {
            printf ("evolving mode k=%e /Mpc  (%d/%d)", ppt->k[index_md][index_k], index_k + 1, ppt->k_size[index_md]);

            if (sgnK != 0)
              printf (" (for scalar modes, corresponds to nu=%e)", sqrt (ppt->k[index_md][index_k] * ppt->k[index_md][index_k] + K) / sqrt (sgnK * K));

            printf ("\n");
          }

#ifdef _OPENMP
          tstart = omp_get_wtime ();
#endif

          class_call_parallel (perturb_solve (ppr,
                                              pba,
                                              pth,
                                              ppt,
                                              index_md,
                                              index_ic,
                                              index_k,
                                              pppw[thread]),
                               ppt->error_message,
                               ppt->error_message);

#ifdef _OPENMP
          tstop = omp_get_wtime ();

          tspent += tstop - tstart;
#endif

#pragma omp flush(abort)
        } /* end of loop over wavenumbers */

#ifdef _OPENMP

        if (ppt->perturbations_verbose > 1)
          printf ("In %s: time spent in parallel region (loop over k's) = %e s for thread %d\n",
                  __func__, tspent, omp_get_thread_num ());

#endif
      } /* end of parallel region */

      if (abort == _TRUE_)
        return _FAILURE_;
    } /* end of loop over initial conditions */

    abort = _FALSE_;

#pragma omp parallel                                    \
    shared(pppw,ppt,index_md,abort,number_of_threads)     \
    private(thread)                                       \
    num_threads(number_of_threads)

    {
#ifdef _OPENMP
      thread = omp_get_thread_num ();
#endif

      class_call_parallel (perturb_workspace_free (ppt, index_md, pppw[thread]),
                           ppt->error_message,
                           ppt->error_message);
    } /* end of parallel region */

    if (abort == _TRUE_)
      return _FAILURE_;
  } /* end loop over modes */

  free (pppw);

  return _SUCCESS_;
}

/**
 * Free all memory space allocated by perturb_init().
 *
 * To be called at the end of each run, only when no further calls to
 * perturb_sources_at_tau() are needed.
 *
 * @param ppt Input: perturbation structure to be freed
 * @return the error status
 */

int
perturb_free (
  struct perturbs *ppt
             )
{
  int index_md, index_ic, index_type;
  int filenum;

  if (ppt->has_perturbations == _TRUE_)
  {
    for (index_md = 0; index_md < ppt->md_size; index_md++)
    {
      for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++)
      {
        for (index_type = 0; index_type < ppt->tp_size[index_md]; index_type++)
        {
          free (ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index_type]);
        }
      }

      free (ppt->sources[index_md]);

      free (ppt->k[index_md]);
    }

    free (ppt->tau_sampling);

    free (ppt->tp_size);

    free (ppt->ic_size);

    free (ppt->k);

    free (ppt->k_size_cmb);

    free (ppt->k_size_cl);

    free (ppt->k_size);

    free (ppt->sources);

    /** Stuff related to perturbations output: */

    /** - Free non-NULL pointers */
    if (ppt->index_k_output_values != NULL)
      free (ppt->index_k_output_values);

    for (filenum = 0; filenum < _MAX_NUMBER_OF_K_FILES_; filenum++)
    {
      if (ppt->scalar_perturbations_data[filenum] != NULL)
        free (ppt->scalar_perturbations_data[filenum]);

      if (ppt->vector_perturbations_data[filenum] != NULL)
        free (ppt->vector_perturbations_data[filenum]);

      if (ppt->tensor_perturbations_data[filenum] != NULL)
        free (ppt->tensor_perturbations_data[filenum]);
    }
  }

  return _SUCCESS_;
}

/**
 * Initialize all indices and allocate most arrays in perturbs structure.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input: pointer to thermodynamics structure
 * @param ppt Input/Output: Initialized perturbation structure
 * @return the error status
 */

int
perturb_indices_of_perturbs (
  struct precision        *ppr,
  const struct background *pba,
  struct thermo           *pth,
  struct perturbs         *ppt
                            )
{
  /** Summary: */

  /** - define local variables */

  int index_type;
  int index_md;
  int index_ic;
  int index_type_common;

  /** - count modes (scalar, vector, tensor) and assign corresponding indices */

  index_md = 0;
  class_define_index (ppt->index_md_scalars, ppt->has_scalars, index_md, 1);
  class_define_index (ppt->index_md_vectors, ppt->has_vectors, index_md, 1);
  class_define_index (ppt->index_md_tensors, ppt->has_tensors, index_md, 1);
  ppt->md_size = index_md;

  class_test (index_md == 0,
              ppt->error_message,
              "you should have at least one out of {scalars, vectors, tensors} !!!");

  /** - allocate array of number of types for each mode, ppt->tp_size[index_md] */

  class_alloc (ppt->tp_size, ppt->md_size * sizeof (int), ppt->error_message);

  /** - allocate array of number of initial conditions for each mode, ppt->ic_size[index_md] */

  class_alloc (ppt->ic_size, ppt->md_size * sizeof (int), ppt->error_message);

  /** - allocate array of arrays of source functions for each mode, ppt->source[index_md] */

  class_alloc (ppt->sources, ppt->md_size * sizeof (double *), ppt->error_message);

  /** - initialization of all flags to false (will eventually be set to true later) */

  ppt->has_cmb = _FALSE_;
  ppt->has_lss = _FALSE_;

  ppt->has_source_t            = _FALSE_;
  ppt->has_source_p            = _FALSE_;
  ppt->has_source_delta_m      = _FALSE_;
  ppt->has_source_delta_cb     = _FALSE_;
  ppt->has_source_delta_g      = _FALSE_;
  ppt->has_source_delta_b      = _FALSE_;
  ppt->has_source_delta_cdm    = _FALSE_;
  ppt->has_source_delta_fld    = _FALSE_;
  ppt->has_source_delta_ur     = _FALSE_;
  ppt->has_source_delta_ncdm   = _FALSE_;
  ppt->has_source_theta_m      = _FALSE_;
  ppt->has_source_theta_cb     = _FALSE_;
  ppt->has_source_theta_g      = _FALSE_;
  ppt->has_source_theta_b      = _FALSE_;
  ppt->has_source_theta_cdm    = _FALSE_;
  ppt->has_source_theta_fld    = _FALSE_;
  ppt->has_source_theta_ur     = _FALSE_;
  ppt->has_source_theta_ncdm   = _FALSE_;
  ppt->has_source_phi          = _FALSE_;
  ppt->has_source_phi_prime    = _FALSE_;
  ppt->has_source_phi_plus_psi = _FALSE_;
  ppt->has_source_psi          = _FALSE_;
  ppt->has_source_h            = _FALSE_;
  ppt->has_source_h_prime      = _FALSE_;
  ppt->has_source_eta          = _FALSE_;
  ppt->has_source_eta_prime    = _FALSE_;

  /** - source flags and indices, for sources that all modes have in
   *   common (temperature, polarization, ...). For temperature, the
   *   term t2 is always non-zero, while other terms are non-zero only
   *   for scalars and vectors. For polarization, the term e is always
   *   non-zero, while the term b is only for vectors and tensors. */

  if (ppt->has_cl_cmb_temperature == _TRUE_)
  {
    ppt->has_source_t = _TRUE_;
    ppt->has_cmb      = _TRUE_;
  }

  if (ppt->has_cl_cmb_polarization == _TRUE_)
  {
    ppt->has_source_p = _TRUE_;
    ppt->has_cmb      = _TRUE_;
  }

  index_type = 0;
  class_define_index (ppt->index_tp_t2, ppt->has_source_t, index_type, 1);
  class_define_index (ppt->index_tp_p, ppt->has_source_p, index_type, 1);
  index_type_common = index_type;

  /* indices for perturbed recombination */

  class_define_index (ppt->index_tp_perturbed_recombination_delta_temp, ppt->has_perturbed_recombination, index_type, 1);
  class_define_index (ppt->index_tp_perturbed_recombination_delta_chi, ppt->has_perturbed_recombination, index_type, 1);




  /** - define k values with perturb_get_k_list() */

  class_call (perturb_get_k_list (ppr,
                                  pba,
                                  pth,
                                  ppt),
              ppt->error_message,
              ppt->error_message);

  /** - loop over modes. Initialize flags and indices which are specific to each mode. */

  for (index_md = 0; index_md < ppt->md_size; index_md++)
  {
    /** - (a) scalars */

    if (_scalars_)
    {
      /** - --> source flags and indices, for sources that are specific to scalars */

      if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) || (ppt->has_cl_lensing_potential))
      {
        ppt->has_lss                 = _TRUE_;
        ppt->has_source_phi_plus_psi = _TRUE_;
      }

      if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_nl_corrections_based_on_delta_m))
      {
        ppt->has_lss            = _TRUE_;
        ppt->has_source_delta_m = _TRUE_;

        if (pba->has_ncdm == _TRUE_)
          ppt->has_source_delta_cb = _TRUE_;
      }

      if (ppt->has_density_transfers == _TRUE_)
      {
        ppt->has_lss            = _TRUE_;
        ppt->has_source_delta_g = _TRUE_;
        ppt->has_source_delta_b = _TRUE_;

        if (pba->has_cdm == _TRUE_)
          ppt->has_source_delta_cdm = _TRUE_;

        if (pba->has_fld == _TRUE_)
          ppt->has_source_delta_fld = _TRUE_;

        if (pba->has_ur == _TRUE_)
          ppt->has_source_delta_ur = _TRUE_;

        if (pba->has_ncdm == _TRUE_)
          ppt->has_source_delta_ncdm = _TRUE_;

        /* Thanks to the following lines, (phi,psi) are also stored as sources */
        /* (Obtained directly in newtonian gauge, infereed from (h,eta) in synchronous gauge). */
        /* If density transfer functions are requested in the (default) CLASS format, */
        /* (phi, psi) will be appended to the delta_i's in the final output. */
        ppt->has_source_phi = _TRUE_;
        ppt->has_source_psi = _TRUE_;
      }

      if (ppt->has_velocity_transfers == _TRUE_)
      {
        ppt->has_lss            = _TRUE_;
        ppt->has_source_theta_g = _TRUE_;
        ppt->has_source_theta_b = _TRUE_;

        if ((pba->has_cdm == _TRUE_) && (ppt->gauge != synchronous))
          ppt->has_source_theta_cdm = _TRUE_;

        if (pba->has_fld == _TRUE_)
          ppt->has_source_theta_fld = _TRUE_;

        if (pba->has_ur == _TRUE_)
          ppt->has_source_theta_ur = _TRUE_;

        if (pba->has_ncdm == _TRUE_)
          ppt->has_source_theta_ncdm = _TRUE_;
      }

      if (ppt->has_cl_number_count == _TRUE_)
      {
        ppt->has_lss = _TRUE_;

        if (ppt->has_nc_density == _TRUE_)
          ppt->has_source_delta_m = _TRUE_;

        if (ppt->has_nc_rsd == _TRUE_)
        {
          ppt->has_source_theta_m = _TRUE_;

          if (pba->has_ncdm == _TRUE_)
            /* we may not need theta_cb at all, rsd always defined for
             *  the total matter, but at least this is made
             *  available */
            ppt->has_source_theta_cb = _TRUE_;
        }

        if (ppt->has_nc_lens == _TRUE_)
          ppt->has_source_phi_plus_psi = _TRUE_;

        if (ppt->has_nc_gr == _TRUE_)
        {
          ppt->has_source_phi          = _TRUE_;
          ppt->has_source_psi          = _TRUE_;
          ppt->has_source_phi_prime    = _TRUE_;
          ppt->has_source_phi_plus_psi = _TRUE_;
        }
      }

      if (ppt->has_metricpotential_transfers == _TRUE_)
      {
        if (ppt->gauge == newtonian)
        {
          ppt->has_source_phi       = _TRUE_;
          ppt->has_source_psi       = _TRUE_;
          ppt->has_source_phi_prime = _TRUE_;
        }

        if (ppt->gauge == synchronous)
        {
          ppt->has_source_h         = _TRUE_;
          ppt->has_source_h_prime   = _TRUE_;
          ppt->has_source_eta       = _TRUE_;
          ppt->has_source_eta_prime = _TRUE_;
        }
      }

      index_type = index_type_common;
      class_define_index (ppt->index_tp_t0,         ppt->has_source_t,         index_type, 1);
      class_define_index (ppt->index_tp_t1,         ppt->has_source_t,         index_type, 1);
      class_define_index (ppt->index_tp_delta_m,    ppt->has_source_delta_m,   index_type, 1);
      class_define_index (ppt->index_tp_delta_cb,   ppt->has_source_delta_cb,  index_type, 1);
      class_define_index (ppt->index_tp_delta_g,    ppt->has_source_delta_g,   index_type, 1);
      class_define_index (ppt->index_tp_delta_b,    ppt->has_source_delta_b,   index_type, 1);
      class_define_index (ppt->index_tp_delta_cdm,  ppt->has_source_delta_cdm, index_type, 1);
      class_define_index (ppt->index_tp_delta_fld,  ppt->has_source_delta_fld, index_type, 1);
      class_define_index (ppt->index_tp_delta_ur,   ppt->has_source_delta_ur,  index_type, 1);
      class_define_index (ppt->index_tp_delta_ncdm1, ppt->has_source_delta_ncdm, index_type, pba->N_ncdm);
      class_define_index (ppt->index_tp_theta_m,    ppt->has_source_theta_m,   index_type, 1);
      class_define_index (ppt->index_tp_theta_cb,   ppt->has_source_theta_cb,  index_type, 1);
      class_define_index (ppt->index_tp_theta_g,    ppt->has_source_theta_g,   index_type, 1);
      class_define_index (ppt->index_tp_theta_b,    ppt->has_source_theta_b,   index_type, 1);
      class_define_index (ppt->index_tp_theta_cdm,  ppt->has_source_theta_cdm, index_type, 1);
      class_define_index (ppt->index_tp_theta_fld,  ppt->has_source_theta_fld, index_type, 1);
      class_define_index (ppt->index_tp_theta_ur,   ppt->has_source_theta_ur,  index_type, 1);
      class_define_index (ppt->index_tp_theta_ncdm1, ppt->has_source_theta_ncdm, index_type, pba->N_ncdm);
      class_define_index (ppt->index_tp_phi,        ppt->has_source_phi,       index_type, 1);
      class_define_index (ppt->index_tp_phi_prime,  ppt->has_source_phi_prime, index_type, 1);
      class_define_index (ppt->index_tp_phi_plus_psi, ppt->has_source_phi_plus_psi, index_type, 1);
      class_define_index (ppt->index_tp_psi,        ppt->has_source_psi,       index_type, 1);
      class_define_index (ppt->index_tp_h,          ppt->has_source_h,         index_type, 1);
      class_define_index (ppt->index_tp_h_prime,    ppt->has_source_h_prime,   index_type, 1);
      class_define_index (ppt->index_tp_eta,        ppt->has_source_eta,       index_type, 1);
      class_define_index (ppt->index_tp_eta_prime,  ppt->has_source_eta_prime, index_type, 1);
      ppt->tp_size[index_md] = index_type;

      class_test (index_type == 0,
                  ppt->error_message,
                  "inconsistent input: you asked for scalars, so you should have at least one non-zero scalar source type (temperature, polarization, lensing/gravitational potential, ...). Please adjust your input.");

      /** - --> count scalar initial conditions (for scalars: ad, cdi, nid, niv; for tensors: only one) and assign corresponding indices */

      index_ic = 0;
      class_define_index (ppt->index_ic_ad, ppt->has_ad, index_ic, 1);
      class_define_index (ppt->index_ic_bi, ppt->has_bi, index_ic, 1);
      class_define_index (ppt->index_ic_cdi, ppt->has_cdi, index_ic, 1);
      class_define_index (ppt->index_ic_nid, ppt->has_nid, index_ic, 1);
      class_define_index (ppt->index_ic_niv, ppt->has_niv, index_ic, 1);
      ppt->ic_size[index_md] = index_ic;

      class_test (index_ic == 0,
                  ppt->error_message,
                  "you should have at least one adiabatic or isocurvature initial condition...} !!!");
    }

    /** - (b) vectors */

    if (_vectors_)
    {
      /** - --> source flags and indices, for sources that are specific to vectors */

      index_type = index_type_common;
      class_define_index (ppt->index_tp_t1, ppt->has_source_t, index_type, 1);
      ppt->tp_size[index_md] = index_type;

      /*
       *  class_test(index_type == 0,
       *          ppt->error_message,
       *          "inconsistent input: you asked for vectors, so you should have at least one non-zero vector source type (temperature or polarization). Please adjust your input.");
       */

      /** - --> initial conditions for vectors*/

      index_ic = 0;
      /* not coded yet */
      ppt->ic_size[index_md] = index_ic;
    }

    /** - (c) tensors */
    if (_tensors_)
    {
      /** - --> source flags and indices, for sources that are specific to tensors */

      index_type = index_type_common;
      /* nothing specific, unlike for vectors and scalars! */
      ppt->tp_size[index_md] = index_type;

      /*
       *  class_test(index_type == 0,
       *          ppt->error_message,
       *          "inconsistent input: you asked for tensors, so you should have at least one non-zero tensor source type (temperature or polarization). Please adjust your input.");
       */

      /** - --> only one initial condition for tensors*/

      index_ic = 0;
      class_define_index (ppt->index_ic_ten, _TRUE_, index_ic, 1);
      ppt->ic_size[index_md] = index_ic;
    }

    /** - (d) for each mode, allocate array of arrays of source functions for each initial conditions and wavenumber, (ppt->source[index_md])[index_ic][index_type] */

    class_alloc (ppt->sources[index_md],
                 ppt->ic_size[index_md] * ppt->tp_size[index_md] * sizeof (double *),
                 ppt->error_message);
  }

  return _SUCCESS_;
}

/**
 * Define time sampling for source functions.
 *
 * For each type, compute the list of values of tau at which sources
 * will be sampled.  Knowing the number of tau values, allocate all
 * arrays of source functions.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param pth Input: pointer to thermodynamics structure
 * @param ppt Input/Output: Initialized perturbation structure
 * @return the error status
 */

int
perturb_timesampling_for_sources (
  struct precision        *ppr,
  const struct background *pba,
  struct thermo           *pth,
  struct perturbs         *ppt
                                 )
{
  /** Summary: */

  /** - define local variables */

  const double conformal_age = nc_scalefactor_eval_eta_Mpc_z (pba->scalefactor, 0.0);

  int counter;
  int index_md;
  int index_type;
  int index_ic;
  int last_index_back;
  int last_index_thermo;
  int first_index_back;
  int first_index_thermo;

  double tau;
  double tau_ini;
  double tau_lower;
  double tau_upper;
  double tau_mid;

  double timescale_source;
  double rate_thermo;
  double rate_isw_squared;
  double a_prime_over_a;
  double a_primeprime_over_a;
  double *pvecback;
  double *pvecthermo;

  /** - allocate background/thermodynamics vectors */

  class_alloc (pvecback, pba->bg_size_short * sizeof (double), ppt->error_message);
  class_alloc (pvecthermo, pth->th_size * sizeof (double), ppt->error_message);

  /** - first, just count the number of sampling points in order to allocate the array containing all values */

  /** - (a) if CMB requested, first sampling point = when the universe
   *   stops being opaque; otherwise, start sampling gravitational
   *   potential at recombination [however, if perturbed recombination
   *   is requested, we also need to start the system before
   *   recombination. Otherwise, the initial conditions for gas
   *   temperature and ionization fraction perturbations (delta_T = 1/3
   *   delta_b, delta_x_e) are not valid]. */

  if ((ppt->has_cmb == _TRUE_) || (ppt->has_perturbed_recombination == _TRUE_))
  {
    /* using bisection, search time tau such that the ratio of thermo
     *  to Hubble time scales tau_c/tau_h=aH/kappa' is equal to
     *  start_sources_at_tau_c_over_tau_h */

    tau_lower = pth->tau_ini;

    class_call (background_at_tau (pba,
                                   tau_lower,
                                   pba->short_info,
                                   pba->inter_normal,
                                   &first_index_back,
                                   pvecback),
                pba->error_message,
                ppt->error_message);

    class_call (thermodynamics_at_z (pba,
                                     pth,
                                     1. / pvecback[pba->index_bg_a] - 1., /* redshift z=1/a-1 */
                                     pth->inter_normal,
                                     &first_index_thermo,
                                     pvecback,
                                     pvecthermo),
                pth->error_message,
                ppt->error_message);

    class_test (pvecback[pba->index_bg_a] *
                pvecback[pba->index_bg_H] /
                pvecthermo[pth->index_th_dkappa] >
                ppr->start_sources_at_tau_c_over_tau_h,
                ppt->error_message,
                "your choice of initial time for computing sources is inappropriate: it corresponds to an earlier time than the one at which the integration of thermodynamical variables started (tau=%g). You should increase either 'start_sources_at_tau_c_over_tau_h' or 'recfast_z_initial'\n",
                tau_lower);


    tau_upper = pth->tau_rec;

    class_call (background_at_tau (pba,
                                   tau_upper,
                                   pba->short_info,
                                   pba->inter_normal,
                                   &first_index_back,
                                   pvecback),
                pba->error_message,
                ppt->error_message);

    class_call (thermodynamics_at_z (pba,
                                     pth,
                                     1. / pvecback[pba->index_bg_a] - 1., /* redshift z=1/a-1 */
                                     pth->inter_normal,
                                     &first_index_thermo,
                                     pvecback,
                                     pvecthermo),
                pth->error_message,
                ppt->error_message);

    class_test (pvecback[pba->index_bg_a] *
                pvecback[pba->index_bg_H] /
                pvecthermo[pth->index_th_dkappa] <
                ppr->start_sources_at_tau_c_over_tau_h,
                ppt->error_message,
                "your choice of initial time for computing sources is inappropriate: it corresponds to a time after recombination. You should decrease 'start_sources_at_tau_c_over_tau_h'\n");

    tau_mid = 0.5 * (tau_lower + tau_upper);

    while (tau_upper - tau_lower > ppr->tol_tau_approx)
    {
      class_call (background_at_tau (pba,
                                     tau_mid,
                                     pba->short_info,
                                     pba->inter_normal,
                                     &first_index_back,
                                     pvecback),
                  pba->error_message,
                  ppt->error_message);

      class_call (thermodynamics_at_z (pba,
                                       pth,
                                       1. / pvecback[pba->index_bg_a] - 1., /* redshift z=1/a-1 */
                                       pth->inter_normal,
                                       &first_index_thermo,
                                       pvecback,
                                       pvecthermo),
                  pth->error_message,
                  ppt->error_message);


      if (pvecback[pba->index_bg_a] *
          pvecback[pba->index_bg_H] /
          pvecthermo[pth->index_th_dkappa] >
          ppr->start_sources_at_tau_c_over_tau_h)

        tau_upper = tau_mid;
      else
        tau_lower = tau_mid;

      tau_mid = 0.5 * (tau_lower + tau_upper);
    }

    tau_ini = tau_mid;
  }
  else
  {
    /* check the time corresponding to the highest redshift requested in output plus one */
    class_call (background_tau_of_z (pba,
                                     ppt->z_max_pk + 1,
                                     &tau_ini),
                pba->error_message,
                ppt->error_message);

    /* obsolete: previous choice was to start always at recombination time */
    /* tau_ini = pth->tau_rec; */

    /* set values of first_index_back/thermo */
    class_call (background_at_tau (pba,
                                   tau_ini,
                                   pba->short_info,
                                   pba->inter_normal,
                                   &first_index_back,
                                   pvecback),
                pba->error_message,
                ppt->error_message);

    class_call (thermodynamics_at_z (pba,
                                     pth,
                                     1. / pvecback[pba->index_bg_a] - 1., /* redshift z=1/a-1 */
                                     pth->inter_normal,
                                     &first_index_thermo,
                                     pvecback,
                                     pvecthermo),
                pth->error_message,
                ppt->error_message);
  }

  /** - (b) next sampling point = previous + ppr->perturb_sampling_stepsize * timescale_source, where:
   *   - --> if CMB requested:
   *   timescale_source1 = \f$ |g/\dot{g}| = |\dot{\kappa}-\ddot{\kappa}/\dot{\kappa}|^{-1} \f$;
   *   timescale_source2 = \f$ |2\ddot{a}/a-(\dot{a}/a)^2|^{-1/2} \f$ (to sample correctly the late ISW effect; and
   *   timescale_source=1/(1/timescale_source1+1/timescale_source2); repeat till today.
   *   - --> if CMB not requested:
   *   timescale_source = 1/aH; repeat till today.
   */

  counter           = 1;
  last_index_back   = first_index_back;
  last_index_thermo = first_index_thermo;
  tau               = tau_ini;

  while (tau < conformal_age)
  {
    class_call (background_at_tau (pba,
                                   tau,
                                   pba->short_info,
                                   pba->inter_closeby,
                                   &last_index_back,
                                   pvecback),
                pba->error_message,
                ppt->error_message);

    class_call (thermodynamics_at_z (pba,
                                     pth,
                                     1. / pvecback[pba->index_bg_a] - 1., /* redshift z=1/a-1 */
                                     pth->inter_closeby,
                                     &last_index_thermo,
                                     pvecback,
                                     pvecthermo),
                pth->error_message,
                ppt->error_message);

    if (ppt->has_cmb == _TRUE_)
    {
      /* variation rate of thermodynamics variables */
      rate_thermo = pvecthermo[pth->index_th_rate];

      /* variation rate of metric due to late ISW effect (important at late times) */
      a_prime_over_a      = pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a];
      a_primeprime_over_a = pvecback[pba->index_bg_H_prime] * pvecback[pba->index_bg_a]
                            + 2. * a_prime_over_a * a_prime_over_a;
      rate_isw_squared = fabs (2. * a_primeprime_over_a - a_prime_over_a * a_prime_over_a);

      /* compute rate */
      timescale_source = sqrt (rate_thermo * rate_thermo + rate_isw_squared);
    }
    else
    {
      /* variation rate given by Hubble time */
      a_prime_over_a = pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a];

      timescale_source = a_prime_over_a;
    }

    /* check it is non-zero */
    class_test (timescale_source == 0.,
                ppt->error_message,
                "null evolution rate, integration is diverging");

    /* compute inverse rate */
    timescale_source = 1. / timescale_source;

    class_test (fabs (ppr->perturb_sampling_stepsize * timescale_source / tau) < ppr->smallest_allowed_variation,
                ppt->error_message,
                "integration step =%e < machine precision : leads either to numerical error or infinite loop", ppr->perturb_sampling_stepsize * timescale_source);

    tau = tau + ppr->perturb_sampling_stepsize * timescale_source;
    counter++;
  }

  /** - --> infer total number of time steps, ppt->tau_size */
  ppt->tau_size = counter;

  /** - --> allocate array of time steps, ppt->tau_sampling[index_tau] */
  class_alloc (ppt->tau_sampling, ppt->tau_size * sizeof (double), ppt->error_message);

  /** - --> repeat the same steps, now filling the array with each tau value: */

  /** - --> (b.1.) first sampling point = when the universe stops being opaque */

  counter                    = 0;
  ppt->tau_sampling[counter] = tau_ini;

  /** - --> (b.2.) next sampling point = previous + ppr->perturb_sampling_stepsize * timescale_source, where
   *   timescale_source1 = \f$ |g/\dot{g}| = |\dot{\kappa}-\ddot{\kappa}/\dot{\kappa}|^{-1} \f$;
   *   timescale_source2 = \f$ |2\ddot{a}/a-(\dot{a}/a)^2|^{-1/2} \f$ (to sample correctly the late ISW effect; and
   *   timescale_source=1/(1/timescale_source1+1/timescale_source2); repeat till today.
   *   If CMB not requested:
   *   timescale_source = 1/aH; repeat till today.  */

  last_index_back   = first_index_back;
  last_index_thermo = first_index_thermo;
  tau               = tau_ini;

  while (tau < conformal_age)
  {
    class_call (background_at_tau (pba,
                                   tau,
                                   pba->short_info,
                                   pba->inter_closeby,
                                   &last_index_back,
                                   pvecback),
                pba->error_message,
                ppt->error_message);

    class_call (thermodynamics_at_z (pba,
                                     pth,
                                     1. / pvecback[pba->index_bg_a] - 1., /* redshift z=1/a-1 */
                                     pth->inter_closeby,
                                     &last_index_thermo,
                                     pvecback,
                                     pvecthermo),
                pth->error_message,
                ppt->error_message);

    if (ppt->has_cmb == _TRUE_)
    {
      /* variation rate of thermodynamics variables */
      rate_thermo = pvecthermo[pth->index_th_rate];

      /* variation rate of metric due to late ISW effect (important at late times) */
      a_prime_over_a      = pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a];
      a_primeprime_over_a = pvecback[pba->index_bg_H_prime] * pvecback[pba->index_bg_a]
                            + 2. * a_prime_over_a * a_prime_over_a;
      rate_isw_squared = fabs (2. * a_primeprime_over_a - a_prime_over_a * a_prime_over_a);

      /* compute rate */
      timescale_source = sqrt (rate_thermo * rate_thermo + rate_isw_squared);
    }
    else
    {
      a_prime_over_a   = pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a];
      timescale_source = a_prime_over_a;
    }

    /* check it is non-zero */
    class_test (timescale_source == 0.,
                ppt->error_message,
                "null evolution rate, integration is diverging");

    /* compute inverse rate */
    timescale_source = 1. / timescale_source;

    class_test (fabs (ppr->perturb_sampling_stepsize * timescale_source / tau) < ppr->smallest_allowed_variation,
                ppt->error_message,
                "integration step =%e < machine precision : leads either to numerical error or infinite loop", ppr->perturb_sampling_stepsize * timescale_source);

    tau = tau + ppr->perturb_sampling_stepsize * timescale_source;
    counter++;
    ppt->tau_sampling[counter] = tau;
  }

  /** - last sampling point = exactly today */
  ppt->tau_sampling[counter] = conformal_age;

  free (pvecback);
  free (pvecthermo);

  /** - loop over modes, initial conditions and types. For each of
   *   them, allocate array of source functions. */

  for (index_md = 0; index_md < ppt->md_size; index_md++)
  {
    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++)
    {
      for (index_type = 0; index_type < ppt->tp_size[index_md]; index_type++)
      {
        class_alloc (ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index_type],
                     ppt->k_size[index_md] * ppt->tau_size * sizeof (double),
                     ppt->error_message);
      }
    }
  }

  return _SUCCESS_;
}

/**
 * Define the number of comoving wavenumbers using the information
 * passed in the precision structure.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to thermodynamics structure
 * @param ppt        Input: pointer to perturbation structure
 * @return the error status
 */

int
perturb_get_k_list (
  struct precision        *ppr,
  const struct background *pba,
  struct thermo           *pth,
  struct perturbs         *ppt
                   )
{
  const double H0            = 1.0 / nc_hicosmo_RH_Mpc (pba->cosmo);
  const double conformal_age = nc_scalefactor_eval_eta_Mpc_z (pba->scalefactor, 0.0);
  int sgnK;
  double K = nc_hicosmo_curvature_K (pba->cosmo, &sgnK);

  int index_k, index_k_output, index_mode;
  double k, k_min = 0., k_rec, step, tau1;
  double *k_max_cmb;
  double *k_max_cl;
  double k_max = 0.;
  double scale2;
  double *tmp_k_list;
  int newk_size, index_newk, add_k_output_value;

  /** Summary: */

  class_test (ppr->k_step_transition == 0.,
              ppt->error_message,
              "stop to avoid division by zero");

  class_test (pth->rs_rec == 0.,
              ppt->error_message,
              "stop to avoid division by zero");

  /** - allocate arrays related to k list for each mode */

  class_alloc (ppt->k_size_cmb,
               ppt->md_size * sizeof (int),
               ppt->error_message);
  class_alloc (ppt->k_size_cl,
               ppt->md_size * sizeof (int),
               ppt->error_message);
  class_alloc (ppt->k_size,
               ppt->md_size * sizeof (int),
               ppt->error_message);
  class_alloc (ppt->k,
               ppt->md_size * sizeof (double *),
               ppt->error_message);

  class_calloc (k_max_cmb,
                ppt->md_size,
                sizeof (double),
                ppt->error_message);
  class_calloc (k_max_cl,
                ppt->md_size,
                sizeof (double),
                ppt->error_message);

  /** - scalar modes */

  if (ppt->has_scalars == _TRUE_)
  {
    /* first value */
    if (sgnK == 0)
      /* K<0 (flat)  : start close to zero */
      k_min = ppr->k_min_tau0 / conformal_age;
    else if (sgnK == -1)
      /* K<0 (open)  : start close to sqrt(-K)
       *  (in transfer modules, for scalars, this will correspond to q close to zero;
       *  for vectors and tensors, this value is even smaller than the minimum necessary value) */
      k_min = sqrt (-K + pow (ppr->k_min_tau0 / conformal_age / pth->angular_rescaling, 2));

    else if (sgnK == 1)
      /* K>0 (closed): start from q=sqrt(k2+(1+m)K) equal to 3sqrt(K), i.e. k=sqrt((8-m)K) */
      k_min = sqrt ((8. - 1.e-4) * K);

    /** - --> find k_max (as well as k_max_cmb[ppt->index_md_scalars], k_max_cl[ppt->index_md_scalars]) */

    k_rec = 2. * _PI_ / pth->rs_rec; /* comoving scale corresponding to sound horizon at recombination */

    k_max_cmb[ppt->index_md_scalars] = k_min;
    k_max_cl[ppt->index_md_scalars]  = k_min;
    k_max                            = k_min;

    if (ppt->has_cls == _TRUE_)
    {
      /* find k_max_cmb[ppt->index_md_scalars] : */

      /* choose a k_max_cmb[ppt->index_md_scalars] corresponding to a wavelength on the last
       *  scattering surface seen today under an angle smaller than
       *  pi/lmax: this is equivalent to
       *  k_max_cl[ppt->index_md_scalars]*[comvoving.ang.diameter.distance] > l_max */

      k_max_cmb[ppt->index_md_scalars] = ppr->k_max_tau0_over_l_max * ppt->l_scalar_max
                                         / conformal_age / pth->angular_rescaling;
      k_max_cl[ppt->index_md_scalars] = k_max_cmb[ppt->index_md_scalars];
      k_max                           = k_max_cmb[ppt->index_md_scalars];

      /* find k_max_cl[ppt->index_md_scalars] : */

      /* if we need density/lensing Cl's, we must impose a stronger condition,
       *  such that the minimum wavelength on the shell corresponding
       *  to the center of smallest redshift bin is seen under an
       *  angle smaller than pi/lmax. So we must multiply our previous
       *  k_max_cl[ppt->index_md_scalars] by the ratio tau0/(tau0-tau[center of smallest
       *  redshift bin]). Note that we could do the same with the
       *  lensing potential if we needed a very precise C_l^phi-phi at
       *  large l. We don't do it by default, because the lensed ClT,
       *  ClE would be marginally affected. */

      if ((ppt->has_cl_number_count == _TRUE_) || (ppt->has_cl_lensing_potential == _TRUE_))
      {
        class_call (background_tau_of_z (pba,
                                         ppt->selection_mean[0],
                                         &tau1),
                    pba->error_message,
                    ppt->error_message);

        k_max_cl[ppt->index_md_scalars] = MAX (k_max_cl[ppt->index_md_scalars], ppr->k_max_tau0_over_l_max * ppt->l_lss_max / (conformal_age - tau1)); /* to be very accurate we should use angular diameter distance to given redshift instead of comoving radius: would implement corrections depending on curvature */
        k_max                           = k_max_cl[ppt->index_md_scalars];
      }
    }

    /* find k_max: */

    if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_))
      k_max = MAX (k_max, ppt->k_max_for_pk);

    if (ppt->has_nl_corrections_based_on_delta_m == _TRUE_)
      k_max = MAX (k_max, ppr->halofit_min_k_max);

    /** - --> test that result for k_min, k_max make sense */

    class_test (k_min < 0.,
                ppt->error_message,
                "buggy definition of k_min");

    class_test (k_max < 0.,
                ppt->error_message,
                "buggy definition of k_max");

    class_test (k_max < k_min,
                ppt->error_message,
                "buggy definition of k_min and/or k_max");

    /* if K>0, the transfer function will be calculated for discrete
     *  integer values of nu=3,4,5,... where nu=sqrt(k2+(1+m)K) and
     *  m=0,1,2 for scalars/vectors/tensors. However we are free to
     *  define in the perturbation module some arbitrary values of k:
     *  later on, the transfer module will interpolate at values of k
     *  corresponding exactly to integer values of nu. Hence, apart
     *  from the value of k_min and the step size in the vicinity of
     *  k_min, we define exactly the same sampling in the three cases
     *  K=0, K<0, K>0 */

    /* allocate array with, for the moment, the largest possible size */
    class_alloc (ppt->k[ppt->index_md_scalars],
                 ((int) ((k_max_cmb[ppt->index_md_scalars] - k_min) / k_rec / MIN (ppr->k_step_super, ppr->k_step_sub)) +
                  (int) (MAX (ppr->k_per_decade_for_pk, ppr->k_per_decade_for_bao) * log (k_max / k_min) / log (10.)) + 3)
                 * sizeof (double), ppt->error_message);

    /* first value */

    index_k                                = 0;
    k                                      = k_min;
    ppt->k[ppt->index_md_scalars][index_k] = k;
    index_k++;

    /* values until k_max_cmb[ppt->index_md_scalars] */

    while (k < k_max_cmb[ppt->index_md_scalars])
    {
      /* the linear step is not constant, it has a step-like shape,
       *  centered around the characteristic scale set by the sound
       *  horizon at recombination (associated to the comoving wavenumber
       *  k_rec) */

      step = (ppr->k_step_super
              + 0.5 * (tanh ((k - k_rec) / k_rec / ppr->k_step_transition) + 1.)
              * (ppr->k_step_sub - ppr->k_step_super)) * k_rec;

      /* there is one other thing to take into account in the step
       *  size. There are two other characteristic scales that matter for
       *  the sampling: the Hubble scale today, k0=a0H0, and eventually
       *  curvature scale sqrt(|K|). We define "scale2" as the sum of the
       *  squared Hubble radius and squared curvature radius. We need to
       *  increase the sampling for k<sqrt(scale2), in order to get the
       *  first mutipoles accurate enough. The formula below reduces it
       *  gradually in the k-->0 limit, by up to a factor 10. The actual
       *  stepsize is still fixed by k_step_super, this is just a
       *  reduction factor. */

      scale2 = pow (pba->a_today * H0, 2) + fabs (K);

      step *= (k * k / scale2 + 1.) / (k * k / scale2 + 1. / ppr->k_step_super_reduction);

      class_test (step / k < ppr->smallest_allowed_variation,
                  ppt->error_message,
                  "k step =%e < machine precision : leads either to numerical error or infinite loop",
                  step * k_rec);

      k += step;

      class_test (k <= ppt->k[ppt->index_md_scalars][index_k - 1],
                  ppt->error_message,
                  "consecutive values of k should differ and should be in growing order");

      ppt->k[ppt->index_md_scalars][index_k] = k;

      index_k++;
    }

    ppt->k_size_cmb[ppt->index_md_scalars] = index_k;

    /* values until k_max_cl[ppt->index_md_scalars] */

    while (k < k_max_cl[ppt->index_md_scalars])
    {
      k *= pow (10., 1. / (ppr->k_per_decade_for_pk
                           + (ppr->k_per_decade_for_bao - ppr->k_per_decade_for_pk)
                           * (1. - tanh (pow ((log (k) - log (ppr->k_bao_center * k_rec)) / log (ppr->k_bao_width), 4)))));

      ppt->k[ppt->index_md_scalars][index_k] = k;
      index_k++;
    }

    ppt->k_size_cl[ppt->index_md_scalars] = index_k;

    /* values until k_max */

    while (k < k_max)
    {
      k *= pow (10., 1. / (ppr->k_per_decade_for_pk
                           + (ppr->k_per_decade_for_bao - ppr->k_per_decade_for_pk)
                           * (1. - tanh (pow ((log (k) - log (ppr->k_bao_center * k_rec)) / log (ppr->k_bao_width), 4)))));

      ppt->k[ppt->index_md_scalars][index_k] = k;
      index_k++;
    }

    ppt->k_size[ppt->index_md_scalars] = index_k;

    class_realloc (ppt->k[ppt->index_md_scalars],
                   ppt->k[ppt->index_md_scalars],
                   ppt->k_size[ppt->index_md_scalars] * sizeof (double),
                   ppt->error_message);
  }

  /** - vector modes */

  if (ppt->has_vectors == _TRUE_)
  {
    /* first value */
    if (sgnK == 0)
      /* K<0 (flat)  : start close to zero */
      k_min = ppr->k_min_tau0 / conformal_age;
    else if (sgnK == -1)
      /* K<0 (open)  : start close to sqrt(-K)
       *  (in transfer modules, for scalars, this will correspond to q close to zero;
       *  for vectors and tensors, this value is even smaller than the minimum necessary value) */
      k_min = sqrt (-K + pow (ppr->k_min_tau0 / conformal_age / pth->angular_rescaling, 2));

    else if (sgnK == 1)
      /* K>0 (closed): start from q=sqrt(k2+(1+m)K) equal to 3sqrt(K), i.e. k=sqrt((8-m)K) */
      k_min = sqrt ((7. - 1.e-4) * K);

    /** - --> find k_max (as well as k_max_cmb[ppt->index_md_vectors], k_max_cl[ppt->index_md_vectors]) */

    k_rec = 2. * _PI_ / pth->rs_rec; /* comoving scale corresponding to sound horizon at recombination */

    k_max_cmb[ppt->index_md_vectors] = k_min;
    k_max_cl[ppt->index_md_vectors]  = k_min;
    k_max                            = k_min;

    if (ppt->has_cls == _TRUE_)
    {
      /* find k_max_cmb: */

      /* choose a k_max_cmb corresponding to a wavelength on the last
       *  scattering surface seen today under an angle smaller than
       *  pi/lmax: this is equivalent to
       *  k_max_cl*[comvoving.ang.diameter.distance] > l_max */

      k_max_cmb[ppt->index_md_vectors] = ppr->k_max_tau0_over_l_max * ppt->l_vector_max
                                         / conformal_age / pth->angular_rescaling;
      k_max_cl[ppt->index_md_vectors] = k_max_cmb[ppt->index_md_vectors];
      k_max                           = k_max_cmb[ppt->index_md_vectors];
    }

    /** - --> test that result for k_min, k_max make sense */

    class_test (k_min < 0.,
                ppt->error_message,
                "buggy definition of k_min");

    class_test (k_max < 0.,
                ppt->error_message,
                "buggy definition of k_max");

    class_test (k_max < k_min,
                ppt->error_message,
                "buggy definition of k_min and/or k_max");

    /* if K>0, the transfer function will be calculated for discrete
     *  integer values of nu=3,4,5,... where nu=sqrt(k2+(1+m)K) and
     *  m=0,1,2 for scalars/vectors/tensors. However we are free to
     *  define in the perturbation module some arbitrary values of k:
     *  later on, the transfer module will interpolate at values of k
     *  corresponding exactly to integer values of nu. Hence, apart
     *  from the value of k_min and the step size in the vicinity of
     *  k_min, we define exactly the same sampling in the three cases
     *  K=0, K<0, K>0 */

    /* allocate array with, for the moment, the largest possible size */
    class_alloc (ppt->k[ppt->index_md_vectors],
                 ((int) ((k_max_cmb[ppt->index_md_vectors] - k_min) / k_rec / MIN (ppr->k_step_super, ppr->k_step_sub)) + 1)
                 * sizeof (double), ppt->error_message);

    /* first value */

    index_k                                = 0;
    k                                      = k_min;
    ppt->k[ppt->index_md_vectors][index_k] = k;
    index_k++;

    /* values until k_max_cmb[ppt->index_md_vectors] */

    while (k < k_max_cmb[ppt->index_md_vectors])
    {
      /* the linear step is not constant, it has a step-like shape,
       *  centered around the characteristic scale set by the sound
       *  horizon at recombination (associated to the comoving wavenumber
       *  k_rec) */

      step = (ppr->k_step_super
              + 0.5 * (tanh ((k - k_rec) / k_rec / ppr->k_step_transition) + 1.)
              * (ppr->k_step_sub - ppr->k_step_super)) * k_rec;

      /* there is one other thing to take into account in the step
       *  size. There are two other characteristic scales that matter for
       *  the sampling: the Hubble scale today, k0=a0H0, and eventually
       *  curvature scale sqrt(|K|). We define "scale2" as the sum of the
       *  squared Hubble radius and squared curvature radius. We need to
       *  increase the sampling for k<sqrt(scale2), in order to get the
       *  first mutipoles accurate enough. The formula below reduces it
       *  gradually in the k-->0 limit, by up to a factor 10. The actual
       *  stepsize is still fixed by k_step_super, this is just a
       *  reduction factor. */

      scale2 = pow (pba->a_today * H0, 2) + fabs (K);

      step *= (k * k / scale2 + 1.) / (k * k / scale2 + 1. / ppr->k_step_super_reduction);

      class_test (step / k < ppr->smallest_allowed_variation,
                  ppt->error_message,
                  "k step =%e < machine precision : leads either to numerical error or infinite loop",
                  step * k_rec);

      k += step;

      class_test (k <= ppt->k[ppt->index_md_scalars][index_k - 1],
                  ppt->error_message,
                  "consecutive values of k should differ and should be in growing order");

      ppt->k[ppt->index_md_vectors][index_k] = k;

      index_k++;
    }

    ppt->k_size_cmb[ppt->index_md_vectors] = index_k;
    ppt->k_size_cl[ppt->index_md_vectors]  = index_k;
    ppt->k_size[ppt->index_md_vectors]     = index_k;

    class_realloc (ppt->k[ppt->index_md_vectors],
                   ppt->k[ppt->index_md_vectors],
                   ppt->k_size[ppt->index_md_vectors] * sizeof (double),
                   ppt->error_message);
  }

  /** - tensor modes */

  if (ppt->has_tensors == _TRUE_)
  {
    /* first value */
    if (sgnK == 0)
      /* K<0 (flat)  : start close to zero */
      k_min = ppr->k_min_tau0 / conformal_age;
    else if (sgnK == -1)
      /* K<0 (open)  : start close to sqrt(-K)
       *  (in transfer modules, for scalars, this will correspond to q close to zero;
       *  for vectors and tensors, this value is even smaller than the minimum necessary value) */
      k_min = sqrt (-K + pow (ppr->k_min_tau0 / conformal_age / pth->angular_rescaling, 2));

    else if (sgnK == 1)
      /* K>0 (closed): start from q=sqrt(k2+(1+m)K) equal to 3sqrt(K), i.e. k=sqrt((8-m)K) */
      k_min = sqrt ((6. - 1.e-4) * K);

    /** - --> find k_max (as well as k_max_cmb[ppt->index_md_tensors], k_max_cl[ppt->index_md_tensors]) */

    k_rec = 2. * _PI_ / pth->rs_rec; /* comoving scale corresponding to sound horizon at recombination */

    k_max_cmb[ppt->index_md_tensors] = k_min;
    k_max_cl[ppt->index_md_tensors]  = k_min;
    k_max                            = k_min;

    if (ppt->has_cls == _TRUE_)
    {
      /* find k_max_cmb[ppt->index_md_tensors]: */

      /* choose a k_max_cmb[ppt->index_md_tensors] corresponding to a wavelength on the last
       *  scattering surface seen today under an angle smaller than
       *  pi/lmax: this is equivalent to
       *  k_max_cl[ppt->index_md_tensors]*[comvoving.ang.diameter.distance] > l_max */

      k_max_cmb[ppt->index_md_tensors] = ppr->k_max_tau0_over_l_max * ppt->l_tensor_max
                                         / conformal_age / pth->angular_rescaling;
      k_max_cl[ppt->index_md_tensors] = k_max_cmb[ppt->index_md_tensors];
      k_max                           = k_max_cmb[ppt->index_md_tensors];
    }

    /** - --> test that result for k_min, k_max make sense */

    class_test (k_min < 0.,
                ppt->error_message,
                "buggy definition of k_min");

    class_test (k_max < 0.,
                ppt->error_message,
                "buggy definition of k_max");

    class_test (k_max < k_min,
                ppt->error_message,
                "buggy definition of k_min and/or k_max");

    /* if K>0, the transfer function will be calculated for discrete
     *  integer values of nu=3,4,5,... where nu=sqrt(k2+(1+m)K) and
     *  m=0,1,2 for scalars/vectors/tensors. However we are free to
     *  define in the perturbation module some arbitrary values of k:
     *  later on, the transfer module will interpolate at values of k
     *  corresponding exactly to integer values of nu. Hence, apart
     *  from the value of k_min and the step size in the vicinity of
     *  k_min, we define exactly the same sampling in the three cases
     *  K=0, K<0, K>0 */

    /* allocate array with, for the moment, the largest possible size */
    class_alloc (ppt->k[ppt->index_md_tensors],
                 ((int) ((k_max_cmb[ppt->index_md_tensors] - k_min) / k_rec / MIN (ppr->k_step_super, ppr->k_step_sub)) + 1)
                 * sizeof (double), ppt->error_message);

    /* first value */

    index_k                                = 0;
    k                                      = k_min;
    ppt->k[ppt->index_md_tensors][index_k] = k;
    index_k++;

    /* values until k_max_cmb[ppt->index_md_tensors] */

    while (k < k_max_cmb[ppt->index_md_tensors])
    {
      /* the linear step is not constant, it has a step-like shape,
       *  centered around the characteristic scale set by the sound
       *  horizon at recombination (associated to the comoving wavenumber
       *  k_rec) */

      step = (ppr->k_step_super
              + 0.5 * (tanh ((k - k_rec) / k_rec / ppr->k_step_transition) + 1.)
              * (ppr->k_step_sub - ppr->k_step_super)) * k_rec;

      /* there is one other thing to take into account in the step
       *  size. There are two other characteristic scales that matter for
       *  the sampling: the Hubble scale today, k0=a0H0, and eventually
       *  curvature scale sqrt(|K|). We define "scale2" as the sum of the
       *  squared Hubble radius and squared curvature radius. We need to
       *  increase the sampling for k<sqrt(scale2), in order to get the
       *  first mutipoles accurate enough. The formula below reduces it
       *  gradually in the k-->0 limit, by up to a factor 10. The actual
       *  stepsize is still fixed by k_step_super, this is just a
       *  reduction factor. */

      scale2 = pow (pba->a_today * H0, 2) + fabs (K);

      step *= (k * k / scale2 + 1.) / (k * k / scale2 + 1. / ppr->k_step_super_reduction);

      class_test (step / k < ppr->smallest_allowed_variation,
                  ppt->error_message,
                  "k step =%e < machine precision : leads either to numerical error or infinite loop",
                  step * k_rec);

      k += step;

      class_test (k <= ppt->k[ppt->index_md_tensors][index_k - 1],
                  ppt->error_message,
                  "consecutive values of k should differ and should be in growing order");

      ppt->k[ppt->index_md_tensors][index_k] = k;

      index_k++;
    }

    ppt->k_size_cmb[ppt->index_md_tensors] = index_k;
    ppt->k_size_cl[ppt->index_md_tensors]  = index_k;
    ppt->k_size[ppt->index_md_tensors]     = index_k;

    class_realloc (ppt->k[ppt->index_md_tensors],
                   ppt->k[ppt->index_md_tensors],
                   ppt->k_size[ppt->index_md_tensors] * sizeof (double),
                   ppt->error_message);
  }

  /** - If user asked for k_output_values, add those to all k lists: */
  if (ppt->k_output_values_num > 0)
  {
    /* Allocate storage */
    class_alloc (ppt->index_k_output_values, sizeof (double) * ppt->md_size * ppt->k_output_values_num, ppt->error_message);

    /** - --> Find indices in ppt->k[index_md] corresponding to 'k_output_values'.
     *   We are assuming that ppt->k is sorted and growing, and we have made sure
     *   that ppt->k_output_values is also sorted and growing.*/
    for (index_mode = 0; index_mode < ppt->md_size; index_mode++)
    {
      newk_size = ppt->k_size[index_mode] + ppt->k_output_values_num;

      class_alloc (tmp_k_list, sizeof (double) * newk_size, ppt->error_message);

      index_k        = 0;
      index_k_output = 0;

      for (index_newk = 0; index_newk < newk_size; index_newk++)
      {
        /** - --> Decide if we should add k_output_value now. This has to be this complicated, since we
         *   can only compare the k-values when both indices are in range.*/
        if (index_k >= ppt->k_size[index_mode])
          add_k_output_value = _TRUE_;
        else if (index_k_output >= ppt->k_output_values_num)
          add_k_output_value = _FALSE_;
        else if (ppt->k_output_values[index_k_output] < ppt->k[index_mode][index_k])
          add_k_output_value = _TRUE_;
        else
          add_k_output_value = _FALSE_;

        if (add_k_output_value == _TRUE_)
        {
          tmp_k_list[index_newk]                                                             = ppt->k_output_values[index_k_output];
          ppt->index_k_output_values[index_mode * ppt->k_output_values_num + index_k_output] = index_newk;
          index_k_output++;
        }
        else
        {
          tmp_k_list[index_newk] = ppt->k[index_mode][index_k];
          index_k++;
        }
      }

      free (ppt->k[index_mode]);
      ppt->k[index_mode]      = tmp_k_list;
      ppt->k_size[index_mode] = newk_size;

      index_k = newk_size - 1;

      while (ppt->k[index_mode][index_k] > k_max_cl[index_mode])
        index_k--;

      ppt->k_size_cl[index_mode] = MIN (index_k + 2, ppt->k_size[index_mode]);

      index_k = newk_size - 1;

      while (ppt->k[index_mode][index_k] > k_max_cmb[index_mode])
        index_k--;

      ppt->k_size_cmb[index_mode] = MIN (index_k + 2, ppt->k_size[index_mode]);

      /** - --> The two MIN statements are here because in a normal run, the cl and cmb
       *   arrays contain a single k value larger than their respective k_max.
       *   We are mimicking this behavior. */
    }
  }

  /* For testing, can be useful to print the k list in a file:
   *
   *  FILE * out=fopen("output/k","w");
   *
   *  for (index_k=0; index_k < ppt->k_size[0]; index_k++) {
   *
   *  fprintf(out,"%e\n",ppt->k[0][index_k],K);
   *
   *  }
   *  fclose(out);
   */

  /** - finally, find the global k_min and k_max for the ensemble of all modes 9scalars, vectors, tensors) */

  ppt->k_min = _HUGE_;
  ppt->k_max = 0.;

  if (ppt->has_scalars == _TRUE_)
  {
    ppt->k_min = MIN (ppt->k_min, ppt->k[ppt->index_md_scalars][0]);                                      /* first value, inferred from perturbations structure */
    ppt->k_max = MAX (ppt->k_max, ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars] - 1]); /* last value, inferred from perturbations structure */
  }

  if (ppt->has_vectors == _TRUE_)
  {
    ppt->k_min = MIN (ppt->k_min, ppt->k[ppt->index_md_vectors][0]);                                      /* first value, inferred from perturbations structure */
    ppt->k_max = MAX (ppt->k_max, ppt->k[ppt->index_md_vectors][ppt->k_size[ppt->index_md_vectors] - 1]); /* last value, inferred from perturbations structure */
  }

  if (ppt->has_tensors == _TRUE_)
  {
    ppt->k_min = MIN (ppt->k_min, ppt->k[ppt->index_md_tensors][0]);                                      /* first value, inferred from perturbations structure */
    ppt->k_max = MAX (ppt->k_max, ppt->k[ppt->index_md_tensors][ppt->k_size[ppt->index_md_tensors] - 1]); /* last value, inferred from perturbations structure */
  }

  free (k_max_cmb);
  free (k_max_cl);

  return _SUCCESS_;
}

/**
 * Initialize a perturb_workspace structure. All fields are allocated
 * here, with the exception of the perturb_vector '-->pv' field, which
 * is allocated separately in perturb_vector_init. We allocate one
 * such perturb_workspace structure per thread and per mode
 * (scalar/../tensor). Then, for each thread, all initial conditions
 * and wavenumbers will use the same workspace.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_md Input: index of mode under consideration (scalar/.../tensor)
 * @param ppw        Input/Output: pointer to perturb_workspace structure which fields are allocated or filled here
 * @return the error status
 */

int
perturb_workspace_init (
  struct precision         *ppr,
  const struct background  *pba,
  struct thermo            *pth,
  struct perturbs          *ppt,
  int                      index_md,
  struct perturb_workspace *ppw
                       )
{
  /** Summary: */

  /** - define local variables */

  int index_mt = 0;
  int index_ap;
  int l;

  /** - Compute maximum l_max for any multipole */;

  if (_scalars_)
  {
    ppw->max_l_max = MAX (ppr->l_max_g, ppr->l_max_pol_g);

    if (pba->has_ur == _TRUE_)
      ppw->max_l_max = MAX (ppw->max_l_max, ppr->l_max_ur);

    if (pba->has_ncdm == _TRUE_)
      ppw->max_l_max = MAX (ppw->max_l_max, ppr->l_max_ncdm);
  }

  if (_tensors_)
  {
    ppw->max_l_max = MAX (ppr->l_max_g_ten, ppr->l_max_pol_g_ten);

    if (pba->has_ur == _TRUE_)
      ppw->max_l_max = MAX (ppw->max_l_max, ppr->l_max_ur);

    if (pba->has_ncdm == _TRUE_)
      ppw->max_l_max = MAX (ppw->max_l_max, ppr->l_max_ncdm);
  }

  /** - Allocate \f$ s_l\f$[ ] array for freestreaming of multipoles (see arXiv:1305.3261) and initialize
   *   to 1.0, which is the K=0 value. */
  class_alloc (ppw->s_l, sizeof (double) * (ppw->max_l_max + 1), ppt->error_message);

  for (l = 0; l <= ppw->max_l_max; l++)
  {
    ppw->s_l[l] = 1.0;
  }

  /** - define indices of metric perturbations obeying constraint
   *   equations (this can be done once and for all, because the
   *   vector of metric perturbations is the same whatever the
   *   approximation scheme, unlike the vector of quantities to
   *   be integrated, which is allocated separately in
   *   perturb_vector_init) */

  if (_scalars_)
  {
    /* newtonian gauge */

    if (ppt->gauge == newtonian)
    {
      class_define_index (ppw->index_mt_psi, _TRUE_, index_mt, 1);       /* psi */
      class_define_index (ppw->index_mt_phi_prime, _TRUE_, index_mt, 1); /* phi' */
    }

    /* synchronous gauge (note that eta is counted in the vector of
     *  quantities to be integrated, while here we only consider
     *  quantities obeying to constraint equations) */

    if (ppt->gauge == synchronous)
    {
      class_define_index (ppw->index_mt_h_prime, _TRUE_, index_mt, 1);       /* h' */
      class_define_index (ppw->index_mt_h_prime_prime, _TRUE_, index_mt, 1); /* h'' */
      class_define_index (ppw->index_mt_eta_prime, _TRUE_, index_mt, 1);     /* eta' */
      class_define_index (ppw->index_mt_alpha, _TRUE_, index_mt, 1);         /* alpha = (h' + 6 tau') / (2 k**2) */
      class_define_index (ppw->index_mt_alpha_prime, _TRUE_, index_mt, 1);   /* alpha' */
    }
  }

  if (_vectors_)
  {
    /* newtonian gauge */

    if (ppt->gauge == newtonian)
      class_define_index (ppw->index_mt_V_prime, _TRUE_, index_mt, 1);


    if (ppt->gauge == synchronous)
      class_define_index (ppw->index_mt_hv_prime_prime, _TRUE_, index_mt, 1);
  }

  if (_tensors_)
    class_define_index (ppw->index_mt_gw_prime_prime, _TRUE_, index_mt, 1);

  ppw->mt_size = index_mt;

  /** - allocate some workspace in which we will store temporarily the
   *   values of background, thermodynamics, metric and source
   *   quantities at a given time */

  class_alloc (ppw->pvecback, pba->bg_size_normal * sizeof (double), ppt->error_message);
  class_alloc (ppw->pvecthermo, pth->th_size * sizeof (double), ppt->error_message);
  class_alloc (ppw->pvecmetric, ppw->mt_size * sizeof (double), ppt->error_message);

  /** - count number of approximations, initialize their indices, and allocate their flags */
  index_ap = 0;

  class_define_index (ppw->index_ap_tca, _TRUE_, index_ap, 1);
  class_define_index (ppw->index_ap_rsa, _TRUE_, index_ap, 1);

  if (_scalars_)
  {
    class_define_index (ppw->index_ap_ufa, pba->has_ur, index_ap, 1);
    class_define_index (ppw->index_ap_ncdmfa, pba->has_ncdm, index_ap, 1);
  }

  ppw->ap_size = index_ap;

  if (ppw->ap_size > 0)
    class_alloc (ppw->approx, ppw->ap_size * sizeof (int), ppt->error_message);

  /** - For definiteness, initialize approximation flags to arbitrary
   *   values (correct values are overwritten in
   *   pertub_find_approximation_switches) */

  if (_scalars_)
  {
    ppw->approx[ppw->index_ap_tca] = (int) tca_on;
    ppw->approx[ppw->index_ap_rsa] = (int) rsa_off;

    if (pba->has_ur == _TRUE_)
      ppw->approx[ppw->index_ap_ufa] = (int) ufa_off;

    if (pba->has_ncdm == _TRUE_)
      ppw->approx[ppw->index_ap_ncdmfa] = (int) ncdmfa_off;
  }

  if (_tensors_)
  {
    ppw->approx[ppw->index_ap_tca] = (int) tca_on;
    ppw->approx[ppw->index_ap_rsa] = (int) rsa_off;
  }

  /** - allocate fields where some of the perturbations are stored */

  if (_scalars_)
  {
    if ((ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_) || (ppt->has_source_delta_m == _TRUE_))
    {
      class_alloc (ppw->delta_ncdm, pba->N_ncdm * sizeof (double), ppt->error_message);
      class_alloc (ppw->theta_ncdm, pba->N_ncdm * sizeof (double), ppt->error_message);
      class_alloc (ppw->shear_ncdm, pba->N_ncdm * sizeof (double), ppt->error_message);
    }
  }

  return _SUCCESS_;
}

/**
 * Free the perturb_workspace structure (with the exception of the
 * perturb_vector '-->pv' field, which is freed separately in
 * perturb_vector_free).
 *
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_md Input: index of mode under consideration (scalar/.../tensor)
 * @param ppw        Input: pointer to perturb_workspace structure to be freed
 * @return the error status
 */

int
perturb_workspace_free (
  struct perturbs          *ppt,
  int                      index_md,
  struct perturb_workspace *ppw
                       )
{
  free (ppw->s_l);
  free (ppw->pvecback);
  free (ppw->pvecthermo);
  free (ppw->pvecmetric);

  if (ppw->ap_size > 0)
    free (ppw->approx);

  if (_scalars_)
  {
    if ((ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_) || (ppt->has_source_delta_m == _TRUE_))
    {
      free (ppw->delta_ncdm);
      free (ppw->theta_ncdm);
      free (ppw->shear_ncdm);
    }
  }

  free (ppw);

  return _SUCCESS_;
}

/**
 * Solve the perturbation evolution for a given mode, initial
 * condition and wavenumber, and compute the corresponding source
 * functions.
 *
 * For a given mode, initial condition and wavenumber, this function
 * finds the time ranges over which the perturbations can be described
 * within a given approximation. For each such range, it initializes
 * (or redistributes) perturbations using perturb_vector_init(), and
 * integrates over time. Whenever a "source sampling time" is passed,
 * the source terms are computed and stored in the source table using
 * perturb_sources().
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppt        Input/Output: pointer to the perturbation structure (output source functions S(k,tau) written here)
 * @param index_md Input: index of mode under consideration (scalar/.../tensor)
 * @param index_ic   Input: index of initial condition under consideration (ad, iso...)
 * @param index_k    Input: index of wavenumber
 * @param ppw        Input: pointer to perturb_workspace structure containing index values and workspaces
 * @return the error status
 */

int
perturb_solve (
  struct precision         *ppr,
  const struct background  *pba,
  struct thermo            *pth,
  struct perturbs          *ppt,
  int                      index_md,
  int                      index_ic,
  int                      index_k,
  struct perturb_workspace *ppw
              )
{
  /** Summary: */

  /** - define local variables */

  int sgnK;
  double K = nc_hicosmo_curvature_K (pba->cosmo, &sgnK);

  /* contains all fixed parameters, indices and workspaces used by the perturb_derivs function */
  struct perturb_parameters_and_workspace ppaw;

  /* conformal time */
  double tau, tau_lower, tau_upper, tau_mid;

  /* multipole */
  int l;

  /* index running over time */
  int index_tau;

  /* number of values in the tau_sampling array that should be considered for a given mode */
  int tau_actual_size;

  /* running index over types (temperature, etc) */
  int index_type;

  /* Fourier mode */
  double k;

  /* number of time intervals where the approximation scheme is uniform */
  int interval_number;

  /* index running over such time intervals */
  int index_interval;

  /* number of time intervals where each particular approximation is uniform */
  int *interval_number_of;

  /* edge of intervals where approximation scheme is uniform: tau_ini, tau_switch_1, ..., tau_end */
  double *interval_limit;

  /* array of approximation scheme within each interval: interval_approx[index_interval][index_ap] */
  int **interval_approx;

  /* index running over approximations */
  int index_ap;

  /* approximation scheme within previous interval: previous_approx[index_ap] */
  int *previous_approx;

  int n_ncdm, is_early_enough;

  /* function pointer to ODE evolver and names of possible evolvers */

  extern int evolver_rk ();
  extern int evolver_ndf15 ();

  int (*generic_evolver) ();


  /* Related to the perturbation output */
  int (*perhaps_print_variables) ();

  int index_ikout;

  /** - initialize indices relevant for back/thermo tables search */
  ppw->last_index_back   = 0;
  ppw->last_index_thermo = 0;
  ppw->inter_mode        = pba->inter_normal;

  /** - get wavenumber value */
  k = ppt->k[index_md][index_k];

  class_test (k == 0.,
              ppt->error_message,
              "stop to avoid division by zero");

  /** - If non-zero curvature, update array of free-streaming coefficients ppw->s_l */
  if (pba->has_curvature == _TRUE_)
  {
    for (l = 0; l <= ppw->max_l_max; l++)
    {
      ppw->s_l[l] = sqrt (MAX (1.0 - K * (l * l - 1.0) / k / k, 0.));
    }
  }

  /** - maximum value of tau for which sources are calculated for this wavenumber */

  /* by default, today */
  tau_actual_size = ppt->tau_size;

  /** - using bisection, compute minimum value of tau for which this
   *   wavenumber is integrated */

  /* will be at least the first time in the background table */
  tau_lower = pba->tau_table[0];

  class_call (background_at_tau (pba,
                                 tau_lower,
                                 pba->normal_info,
                                 pba->inter_normal,
                                 &(ppw->last_index_back),
                                 ppw->pvecback),
              pba->error_message,
              ppt->error_message);

  class_call (thermodynamics_at_z (pba,
                                   pth,
                                   1. / ppw->pvecback[pba->index_bg_a] - 1.,
                                   pth->inter_normal,
                                   &(ppw->last_index_thermo),
                                   ppw->pvecback,
                                   ppw->pvecthermo),
              pth->error_message,
              ppt->error_message);

  /* check that this initial time is indeed OK given imposed
   *  conditions on kappa' and on k/aH */

  class_test (ppw->pvecback[pba->index_bg_a] *
              ppw->pvecback[pba->index_bg_H] /
              ppw->pvecthermo[pth->index_th_dkappa] >
              ppr->start_small_k_at_tau_c_over_tau_h, ppt->error_message, "your choice of initial time for integrating wavenumbers is inappropriate: it corresponds to a time before that at which the background has been integrated. You should increase 'start_small_k_at_tau_c_over_tau_h' up to at least %g, or decrease 'a_ini_over_a_today_default'\n",
              ppw->pvecback[pba->index_bg_a] *
              ppw->pvecback[pba->index_bg_H] /
              ppw->pvecthermo[pth->index_th_dkappa]);

  class_test (k / ppw->pvecback[pba->index_bg_a] / ppw->pvecback[pba->index_bg_H] >
              ppr->start_large_k_at_tau_h_over_tau_k,
              ppt->error_message,
              "your choice of initial time for integrating wavenumbers is inappropriate: it corresponds to a time before that at which the background has been integrated. You should increase 'start_large_k_at_tau_h_over_tau_k' up to at least %g, or decrease 'a_ini_over_a_today_default'\n",
              ppt->k[index_md][ppt->k_size[index_md] - 1] / ppw->pvecback[pba->index_bg_a] / ppw->pvecback[pba->index_bg_H]);

  if (pba->has_ncdm == _TRUE_)
  {
    for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
    {
      class_test (fabs (ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm] / ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm] - 1. / 3.) > ppr->tol_ncdm_initial_w,
                  ppt->error_message,
                  "your choice of initial time for integrating wavenumbers is inappropriate: it corresponds to a time at which the ncdm species number %d is not ultra-relativistic anymore, with w=%g, p=%g and rho=%g\n",
                  n_ncdm,
                  ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm] / ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm],
                  ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm],
                  ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm]);
    }
  }

  /* is at most the time at which sources must be sampled */
  tau_upper = ppt->tau_sampling[0];

  /* start bisection */
  tau_mid = 0.5 * (tau_lower + tau_upper);

  while ((tau_upper - tau_lower) / tau_lower > ppr->tol_tau_approx)
  {
    is_early_enough = _TRUE_;

    class_call (background_at_tau (pba,
                                   tau_mid,
                                   pba->normal_info,
                                   pba->inter_normal,
                                   &(ppw->last_index_back),
                                   ppw->pvecback),
                pba->error_message,
                ppt->error_message);

    /* if there are non-cold relics, check that they are relativistic enough */
    if (pba->has_ncdm == _TRUE_)
    {
      for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
      {
        if (fabs (ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm] / ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm] - 1. / 3.) > ppr->tol_ncdm_initial_w)
          is_early_enough = _FALSE_;
      }
    }

    /* also check that the two conditions on (aH/kappa') and (aH/k) are fulfilled */
    if (is_early_enough == _TRUE_)
    {
      class_call (thermodynamics_at_z (pba,
                                       pth,
                                       1. / ppw->pvecback[pba->index_bg_a] - 1., /* redshift z=1/a-1 */
                                       pth->inter_normal,
                                       &(ppw->last_index_thermo),
                                       ppw->pvecback,
                                       ppw->pvecthermo),
                  pth->error_message,
                  ppt->error_message);

      if ((ppw->pvecback[pba->index_bg_a] *
           ppw->pvecback[pba->index_bg_H] /
           ppw->pvecthermo[pth->index_th_dkappa] >
           ppr->start_small_k_at_tau_c_over_tau_h) ||
          (k / ppw->pvecback[pba->index_bg_a] / ppw->pvecback[pba->index_bg_H] >
           ppr->start_large_k_at_tau_h_over_tau_k))

        is_early_enough = _FALSE_;
    }

    if (is_early_enough == _TRUE_)
      tau_lower = tau_mid;
    else
      tau_upper = tau_mid;

    tau_mid = 0.5 * (tau_lower + tau_upper);
  }

  tau = tau_mid;

  /** - find the number of intervals over which approximation scheme is constant */

  class_alloc (interval_number_of, ppw->ap_size * sizeof (int), ppt->error_message);

  ppw->inter_mode = pba->inter_normal;

  class_call (perturb_find_approximation_number (ppr,
                                                 pba,
                                                 pth,
                                                 ppt,
                                                 index_md,
                                                 k,
                                                 ppw,
                                                 tau,
                                                 ppt->tau_sampling[tau_actual_size - 1],
                                                 &interval_number,
                                                 interval_number_of),
              ppt->error_message,
              ppt->error_message);

  class_alloc (interval_limit, (interval_number + 1) * sizeof (double), ppt->error_message);

  class_alloc (interval_approx, interval_number * sizeof (int *), ppt->error_message);

  for (index_interval = 0; index_interval < interval_number; index_interval++)
    class_alloc (interval_approx[index_interval], ppw->ap_size * sizeof (int), ppt->error_message);

  class_call (perturb_find_approximation_switches (ppr,
                                                   pba,
                                                   pth,
                                                   ppt,
                                                   index_md,
                                                   k,
                                                   ppw,
                                                   tau,
                                                   ppt->tau_sampling[tau_actual_size - 1],
                                                   ppr->tol_tau_approx,
                                                   interval_number,
                                                   interval_number_of,
                                                   interval_limit,
                                                   interval_approx),
              ppt->error_message,
              ppt->error_message);

  free (interval_number_of);

  /** - fill the structure containing all fixed parameters, indices
   *   and workspaces needed by perturb_derivs */

  ppaw.ppr                    = ppr;
  ppaw.pba                    = pba;
  ppaw.pth                    = pth;
  ppaw.ppt                    = ppt;
  ppaw.index_md               = index_md;
  ppaw.index_ic               = index_ic;
  ppaw.index_k                = index_k;
  ppaw.k                      = k;
  ppaw.ppw                    = ppw;
  ppaw.ppw->inter_mode        = pba->inter_closeby;
  ppaw.ppw->last_index_back   = 0;
  ppaw.ppw->last_index_thermo = 0;

  /** - check whether we need to print perturbations to a file for this wavenumber */

  perhaps_print_variables = NULL;
  ppw->index_ikout        = -1;

  for (index_ikout = 0; index_ikout < ppt->k_output_values_num; index_ikout++)
  {
    if (ppt->index_k_output_values[index_md * ppt->k_output_values_num + index_ikout] == index_k)
    {
      ppw->index_ikout        = index_ikout;
      perhaps_print_variables = perturb_print_variables;

      /* class_call(perturb_prepare_output_file(
       *  pba,ppt,ppw,index_ikout,index_md),
       *  ppt->error_message,
       *  ppt->error_message);
       */
    }
  }

  /** - loop over intervals over which approximation scheme is uniform. For each interval: */

  for (index_interval = 0; index_interval < interval_number; index_interval++)
  {
    /** - --> (a) fix the approximation scheme */

    for (index_ap = 0; index_ap < ppw->ap_size; index_ap++)
      ppw->approx[index_ap] = interval_approx[index_interval][index_ap];

    /** - --> (b) get the previous approximation scheme. If the current
     *   interval starts from the initial time tau_ini, the previous
     *   approximation is set to be a NULL pointer, so that the
     *   function perturb_vector_init() knows that perturbations must
     *   be initialized */

    if (index_interval == 0)
      previous_approx = NULL;
    else
      previous_approx = interval_approx[index_interval - 1];


    /** - --> (c) define the vector of perturbations to be integrated
     *   over. If the current interval starts from the initial time
     *   tau_ini, fill the vector with initial conditions for each
     *   mode. If it starts from an approximation switching point,
     *   redistribute correctly the perturbations from the previous to
     *   the new vector of perturbations. */

    class_call (perturb_vector_init (ppr,
                                     pba,
                                     pth,
                                     ppt,
                                     index_md,
                                     index_ic,
                                     k,
                                     interval_limit[index_interval],
                                     ppw,
                                     previous_approx),
                ppt->error_message,
                ppt->error_message);

    /** - --> (d) integrate the perturbations over the current interval. */

    if (ppr->evolver == rk)
      generic_evolver = evolver_rk;
    else
      generic_evolver = evolver_ndf15;

    class_call (generic_evolver (perturb_derivs,
                                 interval_limit[index_interval],
                                 interval_limit[index_interval + 1],
                                 ppw->pv->y,
                                 ppw->pv->used_in_sources,
                                 ppw->pv->pt_size,
                                 &ppaw,
                                 ppr->tol_perturb_integration,
                                 ppr->smallest_allowed_variation,
                                 perturb_timescale,
                                 ppr->perturb_integration_stepsize,
                                 ppt->tau_sampling,
                                 tau_actual_size,
                                 perturb_sources,
                                 perhaps_print_variables,
                                 ppt->error_message),
                ppt->error_message,
                ppt->error_message);
  }

  /** - if perturbations were printed in a file, close the file */

  /*if (perhaps_print_variables != NULL) */
  /*  fclose(ppw->perturb_output_file); */

  /** - fill the source terms array with zeros for all times between
   *   the last integrated time tau_max and tau_today. */

  for (index_tau = tau_actual_size; index_tau < ppt->tau_size; index_tau++)
  {
    for (index_type = 0; index_type < ppt->tp_size[index_md]; index_type++)
    {
      ppt->sources[index_md]
      [index_ic * ppt->tp_size[index_md] + index_type]
      [index_tau * ppt->k_size[index_md] + index_k] = 0.;
    }
  }

  /** - free quantities allocated at the beginning of the routine */

  class_call (perturb_vector_free (ppw->pv),
              ppt->error_message,
              ppt->error_message);

  for (index_interval = 0; index_interval < interval_number; index_interval++)
    free (interval_approx[index_interval]);

  free (interval_approx);

  free (interval_limit);

  return _SUCCESS_;
}

int
perturb_prepare_output (const struct background *pba,
                        struct perturbs         *ppt)
{
  int n_ncdm;
  char tmp[40];

  ppt->scalar_titles[0] = '\0';
  ppt->vector_titles[0] = '\0';
  ppt->tensor_titles[0] = '\0';


  if (ppt->k_output_values_num > 0)
  {
    /** Write titles for all perturbations that we would like to print/store. */
    if (ppt->has_scalars == _TRUE_)
    {
      class_store_columntitle (ppt->scalar_titles, "tau [Mpc]", _TRUE_);
      class_store_columntitle (ppt->scalar_titles, "a", _TRUE_);
      class_store_columntitle (ppt->scalar_titles, "delta_g", _TRUE_);
      class_store_columntitle (ppt->scalar_titles, "theta_g", _TRUE_);
      class_store_columntitle (ppt->scalar_titles, "shear_g", _TRUE_);
      class_store_columntitle (ppt->scalar_titles, "pol0_g", _TRUE_);
      class_store_columntitle (ppt->scalar_titles, "pol1_g", _TRUE_);
      class_store_columntitle (ppt->scalar_titles, "pol2_g", _TRUE_);
      class_store_columntitle (ppt->scalar_titles, "delta_b", _TRUE_);
      class_store_columntitle (ppt->scalar_titles, "theta_b", _TRUE_);
      class_store_columntitle (ppt->scalar_titles, "psi", _TRUE_);
      class_store_columntitle (ppt->scalar_titles, "phi", _TRUE_);
      /* Perturbed recombination */
      class_store_columntitle (ppt->scalar_titles, "delta_Tb", ppt->has_perturbed_recombination);
      class_store_columntitle (ppt->scalar_titles, "delta_chi", ppt->has_perturbed_recombination);
      /* Ultrarelativistic species */
      class_store_columntitle (ppt->scalar_titles, "delta_ur", pba->has_ur);
      class_store_columntitle (ppt->scalar_titles, "theta_ur", pba->has_ur);
      class_store_columntitle (ppt->scalar_titles, "shear_ur", pba->has_ur);
      /* Cold dark matter */
      class_store_columntitle (ppt->scalar_titles, "delta_cdm", pba->has_cdm);
      class_store_columntitle (ppt->scalar_titles, "theta_cdm", pba->has_cdm);

      /* Non-cold dark matter */
      if ((pba->has_ncdm == _TRUE_) && ((ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_) || (ppt->has_source_delta_m == _TRUE_)))
      {
        for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
        {
          sprintf (tmp, "delta_ncdm[%d]", n_ncdm);
          class_store_columntitle (ppt->scalar_titles, tmp, _TRUE_);
          sprintf (tmp, "theta_ncdm[%d]", n_ncdm);
          class_store_columntitle (ppt->scalar_titles, tmp, _TRUE_);
          sprintf (tmp, "shear_ncdm[%d]", n_ncdm);
          class_store_columntitle (ppt->scalar_titles, tmp, _TRUE_);
          sprintf (tmp, "cs2_ncdm[%d]", n_ncdm);
          class_store_columntitle (ppt->scalar_titles, tmp, _TRUE_);
        }
      }

      ppt->number_of_scalar_titles =
        get_number_of_titles (ppt->scalar_titles);
    }

    if (ppt->has_tensors == _TRUE_)
    {
      class_store_columntitle (ppt->tensor_titles, "tau [Mpc]", _TRUE_);
      class_store_columntitle (ppt->tensor_titles, "a", _TRUE_);
      class_store_columntitle (ppt->tensor_titles, "delta_g", _TRUE_);
      class_store_columntitle (ppt->tensor_titles, "shear_g", _TRUE_);
      class_store_columntitle (ppt->tensor_titles, "l4_g", _TRUE_);
      class_store_columntitle (ppt->tensor_titles, "pol0_g", _TRUE_);
      class_store_columntitle (ppt->tensor_titles, "pol2_g", _TRUE_);
      class_store_columntitle (ppt->tensor_titles, "pol4_g", _TRUE_);
      class_store_columntitle (ppt->tensor_titles, "H (gw)", _TRUE_);
      class_store_columntitle (ppt->tensor_titles, "Hdot (gwdot)", _TRUE_);

      class_store_columntitle (ppt->tensor_titles, "delta_ur", ppt->evolve_tensor_ur);
      class_store_columntitle (ppt->tensor_titles, "shear_ur", ppt->evolve_tensor_ur);
      class_store_columntitle (ppt->tensor_titles, "l4_ur", ppt->evolve_tensor_ur);

      if (ppt->evolve_tensor_ncdm == _TRUE_)
      {
        for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
        {
          sprintf (tmp, "delta_ncdm[%d]", n_ncdm);
          class_store_columntitle (ppt->tensor_titles, tmp, _TRUE_);
          sprintf (tmp, "theta_ncdm[%d]", n_ncdm);
          class_store_columntitle (ppt->tensor_titles, tmp, _TRUE_);
          sprintf (tmp, "shear_ncdm[%d]", n_ncdm);
          class_store_columntitle (ppt->tensor_titles, tmp, _TRUE_);
        }
      }

      ppt->number_of_tensor_titles =
        get_number_of_titles (ppt->tensor_titles);
    }
  }

  return _SUCCESS_;
}

/**
 * For a given mode and wavenumber, find the number of intervals of
 * time between tau_ini and tau_end such that the approximation
 * scheme (and the number of perturbation equations) is uniform.
 *
 * @param ppr                Input: pointer to precision structure
 * @param pba                Input: pointer to background structure
 * @param pth                Input: pointer to the thermodynamics structure
 * @param ppt                Input: pointer to the perturbation structure
 * @param index_md           Input: index of mode under consideration (scalar/.../tensor)
 * @param k                  Input: index of wavenumber
 * @param ppw                Input: pointer to perturb_workspace structure containing index values and workspaces
 * @param tau_ini            Input: initial time of the perturbation integration
 * @param tau_end            Input: final time of the perturbation integration
 * @param interval_number    Output: total number of intervals
 * @param interval_number_of Output: number of intervals with respect to each particular approximation
 * @return the error status
 */

int
perturb_find_approximation_number (
  struct precision         *ppr,
  const struct background  *pba,
  struct thermo            *pth,
  struct perturbs          *ppt,
  int                      index_md,
  double                   k,
  struct perturb_workspace *ppw,
  double                   tau_ini,
  double                   tau_end,
  int                      *interval_number,
  int                      *interval_number_of /* interval_number_of[index_ap] (already allocated) */
                                  )
{
  /** Summary: */
  /* index running over approximations */
  int index_ap;

  /* value of a given approximation at tau_ini and tau_end */
  int flag_ini, flag_end;

  /** - fix default number of intervals to one (if no approximation switch) */

  *interval_number = 1;

  /** - loop over each approximation and add the number of approximation switching times */

  for (index_ap = 0; index_ap < ppw->ap_size; index_ap++)
  {
    class_call (perturb_approximations (ppr,
                                        pba,
                                        pth,
                                        ppt,
                                        index_md,
                                        k,
                                        tau_ini,
                                        ppw),
                ppt->error_message,
                ppt->error_message);

    flag_ini = ppw->approx[index_ap];

    class_call (perturb_approximations (ppr,
                                        pba,
                                        pth,
                                        ppt,
                                        index_md,
                                        k,
                                        tau_end,
                                        ppw),
                ppt->error_message,
                ppt->error_message);

    flag_end = ppw->approx[index_ap];

    class_test (flag_end < flag_ini,
                ppt->error_message,
                "For each approximation scheme, the declaration of approximation labels in the enumeration must follow chronological order, e.g: enum approx_flags {flag1, flag2, flag3} with flag1 being the initial one and flag3 the final one");

    *interval_number += flag_end - flag_ini;

    interval_number_of[index_ap] = flag_end - flag_ini + 1;
  }

  return _SUCCESS_;
}

/**
 * For a given mode and wavenumber, find the values of time at which
 * the approximation changes.
 *
 * @param ppr                Input: pointer to precision structure
 * @param pba                Input: pointer to background structure
 * @param pth                Input: pointer to the thermodynamics structure
 * @param ppt                Input: pointer to the perturbation structure
 * @param index_md           Input: index of mode under consideration (scalar/.../tensor)
 * @param k                  Input: index of wavenumber
 * @param ppw                Input: pointer to perturb_workspace structure containing index values and workspaces
 * @param tau_ini            Input: initial time of the perturbation integration
 * @param tau_end            Input: final time of the perturbation integration
 * @param precision          Input: tolerance on output values
 * @param interval_number    Input: total number of intervals
 * @param interval_number_of Input: number of intervals with respect to each particular approximation
 * @param interval_limit     Output: value of time at the boundary of the intervals: tau_ini, tau_switch1, ..., tau_end
 * @param interval_approx    Output: value of approximations in each interval
 * @return the error status
 */

int
perturb_find_approximation_switches (
  struct precision         *ppr,
  const struct background  *pba,
  struct thermo            *pth,
  struct perturbs          *ppt,
  int                      index_md,
  double                   k,
  struct perturb_workspace *ppw,
  double                   tau_ini,
  double                   tau_end,
  double                   precision,
  int                      interval_number,
  int                      *interval_number_of,
  double                   *interval_limit,  /* interval_limit[index_interval] (already allocated) */
  int                      **interval_approx /* interval_approx[index_interval][index_ap] (already allocated) */
                                    )
{
  /** Summary: */

  int index_ap;
  int index_switch;
  int index_switch_tot;
  int num_switch;
  double tau_min, lower_bound, upper_bound;
  double mid = 0;
  double *unsorted_tau_switch;
  double next_tau_switch;
  int flag_ini;
  int num_switching_at_given_time;

  /** - write in output arrays the initial time and approximation */

  interval_limit[0] = tau_ini;

  class_call (perturb_approximations (ppr,
                                      pba,
                                      pth,
                                      ppt,
                                      index_md,
                                      k,
                                      tau_ini,
                                      ppw),
              ppt->error_message,
              ppt->error_message);

  for (index_ap = 0; index_ap < ppw->ap_size; index_ap++)
    interval_approx[0][index_ap] = ppw->approx[index_ap];

  /** - if there are no approximation switches, just write final time and return */

  if (interval_number == 1)
  {
    interval_limit[1] = tau_end;
  }

  /** - if there are switches, consider approximations one after each
   *   other.  Find switching time by bisection. Store all switches in
   *   arbitrary order in array unsorted_tau_switch[ ] */

  else
  {
    class_alloc (unsorted_tau_switch, (interval_number - 1) * sizeof (double), ppt->error_message);

    index_switch_tot = 0;

    for (index_ap = 0; index_ap < ppw->ap_size; index_ap++)
    {
      if (interval_number_of[index_ap] > 1)
      {
        num_switch = interval_number_of[index_ap] - 1;

        tau_min = tau_ini;

        flag_ini = interval_approx[0][index_ap];

        for (index_switch = 0; index_switch < num_switch; index_switch++)
        {
          lower_bound = tau_min;
          upper_bound = tau_end;
          mid         = 0.5 * (lower_bound + upper_bound);

          while (upper_bound - lower_bound > precision)
          {
            class_call (perturb_approximations (ppr,
                                                pba,
                                                pth,
                                                ppt,
                                                index_md,
                                                k,
                                                mid,
                                                ppw),
                        ppt->error_message,
                        ppt->error_message);

            if (ppw->approx[index_ap] > flag_ini + index_switch)
              upper_bound = mid;
            else
              lower_bound = mid;

            mid = 0.5 * (lower_bound + upper_bound);
          }

          unsorted_tau_switch[index_switch_tot] = mid;
          index_switch_tot++;

          tau_min = mid;
        }
      }
    }

    class_test (index_switch_tot != (interval_number - 1),
                ppt->error_message,
                "bug in approximation switch search routine: should have %d = %d",
                index_switch_tot, interval_number - 1);

    /** - now sort interval limits in correct order */

    index_switch_tot = 1;

    while (index_switch_tot < interval_number)
    {
      next_tau_switch = tau_end;

      for (index_switch = 0; index_switch < interval_number - 1; index_switch++)
      {
        if ((unsorted_tau_switch[index_switch] > interval_limit[index_switch_tot - 1]) &&
            (unsorted_tau_switch[index_switch] < next_tau_switch))
          next_tau_switch = unsorted_tau_switch[index_switch];
      }

      interval_limit[index_switch_tot] = next_tau_switch;
      index_switch_tot++;
    }

    interval_limit[index_switch_tot] = tau_end;

    class_test (index_switch_tot != interval_number,
                ppt->error_message,
                "most probably two approximation switching time were found to be equal, which cannot be handled\n");

    /** - store each approximation in chronological order */

    for (index_switch = 1; index_switch < interval_number; index_switch++)
    {
      class_call (perturb_approximations (ppr,
                                          pba,
                                          pth,
                                          ppt,
                                          index_md,
                                          k,
                                          0.5 * (interval_limit[index_switch] + interval_limit[index_switch + 1]),
                                          ppw),

                  ppt->error_message,
                  ppt->error_message);

      for (index_ap = 0; index_ap < ppw->ap_size; index_ap++)
      {
        interval_approx[index_switch][index_ap] = ppw->approx[index_ap];

        /* check here that approximation does not go backward (remember
         *  that by definition the value of an approximation can only
         *  increase) */
        class_test (interval_approx[index_switch][index_ap] < interval_approx[index_switch - 1][index_ap],
                    ppt->error_message,
                    "The approximation with label %d is not defined correctly: it goes backward (from %d to %d) for k=%e and between tau=%e and %e; this cannot be handled\n",
                    index_ap,
                    interval_approx[index_switch - 1][index_ap],
                    interval_approx[index_switch][index_ap],
                    k,
                    0.5 * (interval_limit[index_switch - 1] + interval_limit[index_switch]),
                    0.5 * (interval_limit[index_switch] + interval_limit[index_switch + 1])
                   );
      }

      /* check here that more than one approximation is not switched on at a given time */
      num_switching_at_given_time = 0;

      for (index_ap = 0; index_ap < ppw->ap_size; index_ap++)
      {
        if (interval_approx[index_switch][index_ap] != interval_approx[index_switch - 1][index_ap])
          num_switching_at_given_time++;
      }

      class_test (num_switching_at_given_time != 1,
                  ppt->error_message,
                  "for k=%e, at tau=%g, you switch %d approximations at the same time, this cannot be handled. Usually happens in two cases: triggers for different approximations coincide, or one approx is reversible\n",
                  k,
                  interval_limit[index_switch],
                  num_switching_at_given_time);

      if (ppt->perturbations_verbose > 2)
      {
        if (_scalars_)
        {
          if ((interval_approx[index_switch - 1][ppw->index_ap_tca] == (int) tca_on) &&
              (interval_approx[index_switch][ppw->index_ap_tca] == (int) tca_off))
            fprintf (stdout, "Mode k=%e: will switch off tight-coupling approximation at tau=%e\n", k, interval_limit[index_switch]);

          /*fprintf(stderr,"Mode k=%e: will switch off tight-coupling approximation at tau=%e\n",k,interval_limit[index_switch]);  //TBC */

          if ((interval_approx[index_switch - 1][ppw->index_ap_rsa] == (int) rsa_off) &&
              (interval_approx[index_switch][ppw->index_ap_rsa] == (int) rsa_on))
            fprintf (stdout, "Mode k=%e: will switch on radiation streaming approximation at tau=%e\n", k, interval_limit[index_switch]);

          if (pba->has_ur == _TRUE_)
            if ((interval_approx[index_switch - 1][ppw->index_ap_ufa] == (int) ufa_off) &&
                (interval_approx[index_switch][ppw->index_ap_ufa] == (int) ufa_on))
              fprintf (stdout, "Mode k=%e: will switch on ur fluid approximation at tau=%e\n", k, interval_limit[index_switch]);


          if (pba->has_ncdm == _TRUE_)
            if ((interval_approx[index_switch - 1][ppw->index_ap_ncdmfa] == (int) ncdmfa_off) &&
                (interval_approx[index_switch][ppw->index_ap_ncdmfa] == (int) ncdmfa_on))
              fprintf (stdout, "Mode k=%e: will switch on ncdm fluid approximation at tau=%e\n", k, interval_limit[index_switch]);
        }

        if (_tensors_)
        {
          if ((interval_approx[index_switch - 1][ppw->index_ap_tca] == (int) tca_on) &&
              (interval_approx[index_switch][ppw->index_ap_tca] == (int) tca_off))
            fprintf (stdout, "Mode k=%e: will switch off tight-coupling approximation for tensors at tau=%e\n", k, interval_limit[index_switch]);

          if ((interval_approx[index_switch - 1][ppw->index_ap_rsa] == (int) rsa_off) &&
              (interval_approx[index_switch][ppw->index_ap_rsa] == (int) rsa_on))
            fprintf (stdout, "Mode k=%e: will switch on radiation streaming approximation for tensors at tau=%e\n", k, interval_limit[index_switch]);
        }
      }
    }

    free (unsorted_tau_switch);

    class_call (perturb_approximations (ppr,
                                        pba,
                                        pth,
                                        ppt,
                                        index_md,
                                        k,
                                        tau_end,
                                        ppw),

                ppt->error_message,
                ppt->error_message);
  }

  return _SUCCESS_;
}

/**
 * Initialize the field '-->pv' of a perturb_workspace structure, which
 * is a perturb_vector structure. This structure contains indices and
 * values of all quantities which need to be integrated with respect
 * to time (and only them: quantities fixed analytically or obeying
 * constraint equations are NOT included in this vector). This routine
 * distinguishes between two cases:
 *
 * --> the input pa_old is set to the NULL pointer:
 *
 * This happens when we start integrating over a new wavenumber and we
 * want to set initial conditions for the perturbations. Then, it is
 * assumed that ppw-->pv is not yet allocated. This routine allocates
 * it, defines all indices, and then fills the vector ppw-->pv-->y with
 * the initial conditions defined in perturb_initial_conditions.
 *
 * --> the input pa_old is not set to the NULL pointer and describes
 * some set of approximations:
 *
 * This happens when we need to change approximation scheme while
 * integrating over a given wavenumber. The new approximation
 * described by ppw-->pa is then different from pa_old. Then, this
 * routine allocates a new vector with a new size and new index
 * values; it fills this vector with initial conditions taken from the
 * previous vector passed as an input in ppw-->pv, and eventually with
 * some analytic approximations for the new variables appearing at
 * this time; then the new vector comes in replacement of the old one,
 * which is freed.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to the thermodynamics structure
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_md Input: index of mode under consideration (scalar/.../tensor)
 * @param index_ic   Input: index of initial condition under consideration (ad, iso...)
 * @param k          Input: wavenumber
 * @param tau        Input: conformal time
 * @param ppw        Input/Output: workspace containing in input the approximation scheme, the background/thermodynamics/metric quantities, and eventually the previous vector y; and in output the new vector y.
 * @param pa_old     Input: NULL is we need to set y to initial conditions for a new wavenumber; points towards a perturb_approximations if we want to switch of approximation.
 * @return the error status
 */

int
perturb_vector_init (
  struct precision         *ppr,
  const struct background  *pba,
  struct thermo            *pth,
  struct perturbs          *ppt,
  int                      index_md,
  int                      index_ic,
  double                   k,
  double                   tau,
  struct perturb_workspace *ppw, /* ppw->pv unallocated if pa_old = NULL, allocated and filled otherwise */
  int                      *pa_old
                    )
{
  /** Summary: */

  /** - define local variables */

  struct perturb_vector *ppv;

  int index_pt;
  int l;
  int n_ncdm, index_q, ncdm_l_size;
  double rho_plus_p_ncdm, q, q2, epsilon, a, factor;

  /** - allocate a new perturb_vector structure to which ppw-->pv will point at the end of the routine */

  class_alloc (ppv, sizeof (struct perturb_vector), ppt->error_message);

  /** - initialize pointers to NULL (they will be allocated later if
   *   needed), relevant for perturb_vector_free() */
  ppv->l_max_ncdm  = NULL;
  ppv->q_size_ncdm = NULL;

  /** - define all indices in this new vector (depends on approximation scheme, described by the input structure ppw-->pa) */

  index_pt = 0;

  if (_scalars_)
  {
    /* reject inconsistent values of the number of mutipoles in photon temperature hierarchy */
    class_test (ppr->l_max_g < 4,
                ppt->error_message,
                "ppr->l_max_g should be at least 4, i.e. we must integrate at least over photon density, velocity, shear, third and fourth momentum");

    /* reject inconsistent values of the number of mutipoles in photon polarization hierarchy */
    class_test (ppr->l_max_pol_g < 4,
                ppt->error_message,
                "ppr->l_max_pol_g should be at least 4");

    /* reject inconsistent values of the number of mutipoles in ultra relativistic neutrino hierarchy */
    if (pba->has_ur == _TRUE_)
      class_test (ppr->l_max_ur < 4,
                  ppt->error_message,
                  "ppr->l_max_ur should be at least 4, i.e. we must integrate at least over neutrino/relic density, velocity, shear, third and fourth momentum");

    /* photons */

    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off) /* if radiation streaming approximation is off */
    /* temperature */
    {
      ppv->l_max_g = ppr->l_max_g;

      class_define_index (ppv->index_pt_delta_g, _TRUE_, index_pt, 1); /* photon density */
      class_define_index (ppv->index_pt_theta_g, _TRUE_, index_pt, 1); /* photon velocity */

      if (ppw->approx[ppw->index_ap_tca] == (int) tca_off)
      {
        class_define_index (ppv->index_pt_shear_g, _TRUE_, index_pt, 1);             /* photon shear */
        class_define_index (ppv->index_pt_l3_g, _TRUE_, index_pt, ppv->l_max_g - 2); /* higher momenta */

        /* polarization */

        ppv->l_max_pol_g = ppr->l_max_pol_g;

        class_define_index (ppv->index_pt_pol0_g, _TRUE_, index_pt, 1);
        class_define_index (ppv->index_pt_pol1_g, _TRUE_, index_pt, 1);
        class_define_index (ppv->index_pt_pol2_g, _TRUE_, index_pt, 1);
        class_define_index (ppv->index_pt_pol3_g, _TRUE_, index_pt, ppv->l_max_pol_g - 2);
      }
    }

    /* baryons */

    class_define_index (ppv->index_pt_delta_b, _TRUE_, index_pt, 1); /* baryon density */
    class_define_index (ppv->index_pt_theta_b, _TRUE_, index_pt, 1); /* baryon velocity */

    /* cdm */

    class_define_index (ppv->index_pt_delta_cdm, pba->has_cdm, index_pt, 1);                              /* cdm density */
    class_define_index (ppv->index_pt_theta_cdm, pba->has_cdm && (ppt->gauge == newtonian), index_pt, 1); /* cdm velocity */

    /* fluid */

    if (pba->use_ppf == _FALSE_)
    {
      class_define_index (ppv->index_pt_delta_fld, pba->has_fld, index_pt, 1); /* fluid density */
      class_define_index (ppv->index_pt_theta_fld, pba->has_fld, index_pt, 1); /* fluid velocity */
    }
    else
    {
      class_define_index (ppv->index_pt_Gamma_fld, pba->has_fld, index_pt, 1); /* Gamma variable of PPF scheme */
    }

    /* perturbed recombination: the indices are defined once tca is off. */
    if ((ppt->has_perturbed_recombination == _TRUE_) && (ppw->approx[ppw->index_ap_tca] == (int) tca_off))
    {
      class_define_index (ppv->index_pt_perturbed_recombination_delta_temp, _TRUE_, index_pt, 1);
      class_define_index (ppv->index_pt_perturbed_recombination_delta_chi, _TRUE_, index_pt, 1);
    }

    /* ultra relativistic neutrinos */

    if (pba->has_ur && (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off))
    {
      class_define_index (ppv->index_pt_delta_ur, _TRUE_, index_pt, 1); /* density of ultra-relativistic neutrinos/relics */
      class_define_index (ppv->index_pt_theta_ur, _TRUE_, index_pt, 1); /* velocity of ultra-relativistic neutrinos/relics */
      class_define_index (ppv->index_pt_shear_ur, _TRUE_, index_pt, 1); /* shear of ultra-relativistic neutrinos/relics */

      if (ppw->approx[ppw->index_ap_ufa] == (int) ufa_off)
      {
        ppv->l_max_ur = ppr->l_max_ur;
        class_define_index (ppv->index_pt_l3_ur, _TRUE_, index_pt, ppv->l_max_ur - 2); /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */
      }
    }

    /* non-cold dark matter */

    if (pba->has_ncdm == _TRUE_)
    {
      ppv->index_pt_psi0_ncdm1 = index_pt; /* density of ultra-relativistic neutrinos/relics */
      ppv->N_ncdm              = pba->N_ncdm;
      class_alloc (ppv->l_max_ncdm, ppv->N_ncdm * sizeof (double), ppt->error_message);
      class_alloc (ppv->q_size_ncdm, ppv->N_ncdm * sizeof (double), ppt->error_message);

      for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
      {
        /* Set value of ppv->l_max_ncdm: */
        if (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_off)
        {
          /* reject inconsistent values of the number of mutipoles in ultra relativistic neutrino hierarchy */
          class_test (ppr->l_max_ncdm < 4,
                      ppt->error_message,
                      "ppr->l_max_ncdm=%d should be at least 4, i.e. we must integrate at least over first four momenta of non-cold dark matter perturbed phase-space distribution", n_ncdm);
          /*Copy value from precision parameter: */
          ppv->l_max_ncdm[n_ncdm]  = ppr->l_max_ncdm;
          ppv->q_size_ncdm[n_ncdm] = pba->q_size_ncdm[n_ncdm];
        }
        else
        {
          /* In the fluid approximation, hierarchy is cut at lmax = 2 and q dependence is integrated out: */
          ppv->l_max_ncdm[n_ncdm]  = 2;
          ppv->q_size_ncdm[n_ncdm] = 1;
        }

        index_pt += (ppv->l_max_ncdm[n_ncdm] + 1) * ppv->q_size_ncdm[n_ncdm];
      }
    }

    /* metric (only quantities to be integrated, not those obeying constraint equations) */

    /* metric perturbation eta of synchronous gauge */
    class_define_index (ppv->index_pt_eta, ppt->gauge == synchronous, index_pt, 1);

    /* metric perturbation phi of newtonian gauge ( we could fix it
     *  using Einstein equations as a constraint equation for phi, but
     *  integration is numerically more stable if we actually evolve
     *  phi) */
    class_define_index (ppv->index_pt_phi, ppt->gauge == newtonian, index_pt, 1);
  }

  if (_vectors_)
  {
    /* Vector baryon velocity: v_b^{(1)}. */
    class_define_index (ppv->index_pt_theta_b, _TRUE_, index_pt, 1);

    /* eventually reject inconsistent values of the number of mutipoles in photon temperature hierarchy and polarization*/

    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off) /* if radiation streaming approximation is off */
    {
      if (ppw->approx[ppw->index_ap_tca] == (int) tca_off) /* if tight-coupling approximation is off */
      {
        ppv->l_max_g = ppr->l_max_g_ten;

        class_define_index (ppv->index_pt_delta_g, _TRUE_, index_pt, 1);             /* photon density */
        class_define_index (ppv->index_pt_theta_g, _TRUE_, index_pt, 1);             /* photon velocity */
        class_define_index (ppv->index_pt_shear_g, _TRUE_, index_pt, 1);             /* photon shear */
        class_define_index (ppv->index_pt_l3_g, _TRUE_, index_pt, ppv->l_max_g - 2); /* photon l=3 */

        ppv->l_max_pol_g = ppr->l_max_pol_g_ten;

        class_define_index (ppv->index_pt_pol0_g, _TRUE_, index_pt, 1);                    /* photon polarization, l=0 */
        class_define_index (ppv->index_pt_pol1_g, _TRUE_, index_pt, 1);                    /* photon polarization, l=1 */
        class_define_index (ppv->index_pt_pol2_g, _TRUE_, index_pt, 1);                    /* photon polarization, l=2 */
        class_define_index (ppv->index_pt_pol3_g, _TRUE_, index_pt, ppv->l_max_pol_g - 2); /* photon polarization, l=3 */
      }
    }

    /** - (a) metric perturbations V or \f$ h_v \f$ depending on gauge */
    if (ppt->gauge == synchronous)
      class_define_index (ppv->index_pt_hv_prime, _TRUE_, index_pt, 1);

    if (ppt->gauge == newtonian)
      class_define_index (ppv->index_pt_V, _TRUE_, index_pt, 1);
  }

  if (_tensors_)
  {
    /* reject inconsistent values of the number of mutipoles in photon temperature hierarchy */
    class_test (ppr->l_max_g_ten < 4,
                ppt->error_message,
                "ppr->l_max_g_ten should be at least 4, i.e. we must integrate at least over photon density, velocity, shear, third momentum");

    /* reject inconsistent values of the number of mutipoles in photon polarization hierarchy */
    class_test (ppr->l_max_pol_g_ten < 4,
                ppt->error_message,
                "ppr->l_max_pol_g_ten should be at least 4");

    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off) /* if radiation streaming approximation is off */
    {
      if (ppw->approx[ppw->index_ap_tca] == (int) tca_off) /* if tight-coupling approximation is off */
      {
        ppv->l_max_g = ppr->l_max_g_ten;

        class_define_index (ppv->index_pt_delta_g, _TRUE_, index_pt, 1);             /* photon density */
        class_define_index (ppv->index_pt_theta_g, _TRUE_, index_pt, 1);             /* photon velocity */
        class_define_index (ppv->index_pt_shear_g, _TRUE_, index_pt, 1);             /* photon shear */
        class_define_index (ppv->index_pt_l3_g, _TRUE_, index_pt, ppv->l_max_g - 2); /* photon l=3 */

        ppv->l_max_pol_g = ppr->l_max_pol_g_ten;

        class_define_index (ppv->index_pt_pol0_g, _TRUE_, index_pt, 1);                    /* photon polarization, l=0 */
        class_define_index (ppv->index_pt_pol1_g, _TRUE_, index_pt, 1);                    /* photon polarization, l=1 */
        class_define_index (ppv->index_pt_pol2_g, _TRUE_, index_pt, 1);                    /* photon polarization, l=2 */
        class_define_index (ppv->index_pt_pol3_g, _TRUE_, index_pt, ppv->l_max_pol_g - 2); /* photon polarization, l=3 */
      }
    }

    /* ultra relativistic neutrinos */

    class_define_index (ppv->index_pt_delta_ur, ppt->evolve_tensor_ur, index_pt, 1); /* ur density  */
    class_define_index (ppv->index_pt_theta_ur, ppt->evolve_tensor_ur, index_pt, 1); /* ur velocity */
    class_define_index (ppv->index_pt_shear_ur, ppt->evolve_tensor_ur, index_pt, 1); /* ur shear */
    ppv->l_max_ur = ppr->l_max_ur;
    class_define_index (ppv->index_pt_l3_ur, ppt->evolve_tensor_ur, index_pt, ppv->l_max_ur - 2); /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3) */

    if (ppt->evolve_tensor_ncdm == _TRUE_)
    {
      ppv->index_pt_psi0_ncdm1 = index_pt;
      ppv->N_ncdm              = pba->N_ncdm;
      class_alloc (ppv->l_max_ncdm, ppv->N_ncdm * sizeof (double), ppt->error_message);
      class_alloc (ppv->q_size_ncdm, ppv->N_ncdm * sizeof (double), ppt->error_message);

      for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
      {
        /* Set value of ppv->l_max_ncdm: */
        class_test (ppr->l_max_ncdm < 4,
                    ppt->error_message,
                    "ppr->l_max_ncdm=%d should be at least 4, i.e. we must integrate at least over first four momenta of non-cold dark matter perturbed phase-space distribution", n_ncdm);
        /*Copy value from precision parameter: */
        ppv->l_max_ncdm[n_ncdm]  = ppr->l_max_ncdm;
        ppv->q_size_ncdm[n_ncdm] = pba->q_size_ncdm[n_ncdm];

        index_pt += (ppv->l_max_ncdm[n_ncdm] + 1) * ppv->q_size_ncdm[n_ncdm];
      }
    }


    /** - (b) metric perturbation h is a propagating degree of freedom, so h and hdot are included
     *   in the vector of ordinary perturbations, no in that of metric perturbations */

    class_define_index (ppv->index_pt_gw, _TRUE_, index_pt, 1);    /* tensor metric perturbation h (gravitational waves) */
    class_define_index (ppv->index_pt_gwdot, _TRUE_, index_pt, 1); /* its time-derivative */
  }

  ppv->pt_size = index_pt;

  /** - allocate vectors for storing the values of all these
   *   quantities and their time-derivatives at a given time */

  class_calloc (ppv->y, ppv->pt_size, sizeof (double), ppt->error_message);
  class_alloc (ppv->dy, ppv->pt_size * sizeof (double), ppt->error_message);
  class_alloc (ppv->used_in_sources, ppv->pt_size * sizeof (int), ppt->error_message);

  /** - specify which perturbations are needed in the evaluation of source terms */

  /* take all of them by default */
  for (index_pt = 0; index_pt < ppv->pt_size; index_pt++)
    ppv->used_in_sources[index_pt] = _TRUE_;

  /* indicate which ones are not needed (this is just for saving time,
   *  omitting perturbations in this list will not change the
   *  results!) */

  if (_scalars_)
  {
    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
    {
      if (ppw->approx[ppw->index_ap_tca] == (int) tca_off)
      {
        /* we don't need temperature multipoles above l=2 (but they are
         *  defined only when rsa and tca are off) */

        for (index_pt = ppv->index_pt_l3_g; index_pt <= ppv->index_pt_delta_g + ppv->l_max_g; index_pt++)
          ppv->used_in_sources[index_pt] = _FALSE_;

        /* for polarization, we only need l=0,2 (but l =1,3, ... are
         *  defined only when rsa and tca are off) */

        ppv->used_in_sources[ppv->index_pt_pol1_g] = _FALSE_;

        for (index_pt = ppv->index_pt_pol3_g; index_pt <= ppv->index_pt_pol0_g + ppv->l_max_pol_g; index_pt++)
          ppv->used_in_sources[index_pt] = _FALSE_;
      }
    }

    if (pba->has_ur == _TRUE_)
      if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
        if (ppw->approx[ppw->index_ap_ufa] == (int) ufa_off)
          /* we don't need ur multipoles above l=2 (but they are
           *  defined only when rsa and ufa are off) */

          for (index_pt = ppv->index_pt_l3_ur; index_pt <= ppv->index_pt_delta_ur + ppv->l_max_ur; index_pt++)
            ppv->used_in_sources[index_pt] = _FALSE_;


    if (pba->has_ncdm == _TRUE_)
    {
      /* we don't need ncdm multipoles above l=2 (but they are
       *  defined only when ncdmfa is off) */

      index_pt = ppv->index_pt_psi0_ncdm1;

      for (n_ncdm = 0; n_ncdm < ppv->N_ncdm; n_ncdm++)
      {
        for (index_q = 0; index_q < ppv->q_size_ncdm[n_ncdm]; index_q++)
        {
          for (l = 0; l <= ppv->l_max_ncdm[n_ncdm]; l++)
          {
            if (l > 2)
              ppv->used_in_sources[index_pt] = _FALSE_;

            index_pt++;
          }
        }
      }
    }
  }

  if (_tensors_)
  {
    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off) /* if radiation streaming approximation is off */
    {
      if (ppw->approx[ppw->index_ap_tca] == (int) tca_off)
      {
        /* we don't need temperature multipoles above except l=0,2,4 */

        ppv->used_in_sources[ppv->index_pt_theta_g] = _FALSE_;
        ppv->used_in_sources[ppv->index_pt_l3_g]    = _FALSE_;

        for (index_pt = ppv->index_pt_delta_g + 5; index_pt <= ppv->index_pt_delta_g + ppv->l_max_g; index_pt++)
          ppv->used_in_sources[index_pt] = _FALSE_;

        /* same for polarization, we only need l=0,2,4 */

        ppv->used_in_sources[ppv->index_pt_pol1_g] = _FALSE_;
        ppv->used_in_sources[ppv->index_pt_pol3_g] = _FALSE_;

        for (index_pt = ppv->index_pt_pol0_g + 5; index_pt <= ppv->index_pt_pol0_g + ppv->l_max_pol_g; index_pt++)
          ppv->used_in_sources[index_pt] = _FALSE_;
      }
    }

    /* we need h' but not h */
    ppv->used_in_sources[ppv->index_pt_gw] = _FALSE_;
  }

  /** - case of setting initial conditions for a new wavenumber */

  if (pa_old == NULL)
  {
    if (ppt->perturbations_verbose > 2)
      fprintf (stdout, "Mode k=%e: initializing vector at tau=%e\n", k, tau);

    if (_scalars_)
    {
      /** - --> (a) check that current approximation scheme is consistent
       *   with initial conditions */

      class_test (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on,
                  ppt->error_message,
                  "scalar initial conditions assume radiation streaming approximation turned off");


      if (pba->has_ur == _TRUE_)
        class_test (ppw->approx[ppw->index_ap_ufa] == (int) ufa_on,
                    ppt->error_message,
                    "scalar initial conditions assume ur fluid approximation turned off");


      if (pba->has_ncdm == _TRUE_)
        class_test (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on,
                    ppt->error_message,
                    "scalar initial conditions assume ncdm fluid approximation turned off");


      class_test (ppw->approx[ppw->index_ap_tca] == (int) tca_off,
                  ppt->error_message,
                  "scalar initial conditions assume tight-coupling approximation turned on");
    }

    if (_tensors_)
    {
      class_test (ppw->approx[ppw->index_ap_tca] == (int) tca_off,
                  ppt->error_message,
                  "tensor initial conditions assume tight-coupling approximation turned on");

      class_test (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on,
                  ppt->error_message,
                  "tensor initial conditions assume radiation streaming approximation turned off");
    }

    /** - --> (b) let ppw-->pv points towards the perturb_vector structure
     *   that we just created */

    ppw->pv = ppv;

    /** - --> (c) fill the vector ppw-->pv-->y with appropriate initial conditions */

    class_call (perturb_initial_conditions (ppr,
                                            pba,
                                            ppt,
                                            index_md,
                                            index_ic,
                                            k,
                                            tau,
                                            ppw),
                ppt->error_message,
                ppt->error_message);
  }
  /** - case of switching approximation while a wavenumber is being integrated */

  else
  {
    /** - --> (a) for the scalar mode: */

    if (_scalars_)
    {
      /** - ---> (a.1.) check that the change of approximation scheme makes
       *   sense (note: before calling this routine there is already a
       *   check that we wish to change only one approximation flag at
       *   a time) */

      class_test ((pa_old[ppw->index_ap_tca] == (int) tca_off) && (ppw->approx[ppw->index_ap_tca] == (int) tca_on),
                  ppt->error_message,
                  "at tau=%g: the tight-coupling approximation can be switched off, not on", tau);

      /** - ---> (a.2.) some variables (b, cdm, fld, ...) are not affected by
       *   any approximation. They need to be reconducted whatever
       *   the approximation switching is. We treat them here. Below
       *   we will treat other variables case by case. */

      ppv->y[ppv->index_pt_delta_b] =
        ppw->pv->y[ppw->pv->index_pt_delta_b];

      ppv->y[ppv->index_pt_theta_b] =
        ppw->pv->y[ppw->pv->index_pt_theta_b];

      if (pba->has_cdm == _TRUE_)
      {
        ppv->y[ppv->index_pt_delta_cdm] =
          ppw->pv->y[ppw->pv->index_pt_delta_cdm];

        if (ppt->gauge == newtonian)
          ppv->y[ppv->index_pt_theta_cdm] =
            ppw->pv->y[ppw->pv->index_pt_theta_cdm];
      }

      if (pba->has_fld == _TRUE_)
      {
        if (pba->use_ppf == _FALSE_)
        {
          ppv->y[ppv->index_pt_delta_fld] =
            ppw->pv->y[ppw->pv->index_pt_delta_fld];

          ppv->y[ppv->index_pt_theta_fld] =
            ppw->pv->y[ppw->pv->index_pt_theta_fld];
        }
        else
        {
          ppv->y[ppv->index_pt_Gamma_fld] =
            ppw->pv->y[ppw->pv->index_pt_Gamma_fld];
        }
      }

      if (ppt->gauge == synchronous)
        ppv->y[ppv->index_pt_eta] =
          ppw->pv->y[ppw->pv->index_pt_eta];

      if (ppt->gauge == newtonian)
        ppv->y[ppv->index_pt_phi] =
          ppw->pv->y[ppw->pv->index_pt_phi];

      /* -- case of switching off tight coupling
       *  approximation. Provide correct initial conditions to new set
       *  of variables */

      if ((pa_old[ppw->index_ap_tca] == (int) tca_on) && (ppw->approx[ppw->index_ap_tca] == (int) tca_off))
      {
        if (ppt->perturbations_verbose > 2)
          fprintf (stdout, "Mode k=%e: switch off tight-coupling approximation at tau=%e\n", k, tau);

        ppv->y[ppv->index_pt_delta_g] =
          ppw->pv->y[ppw->pv->index_pt_delta_g];

        ppv->y[ppv->index_pt_theta_g] =
          ppw->pv->y[ppw->pv->index_pt_theta_g];

        /* tight-coupling approximation for shear_g (previously
         *  computed in perturb_derivs: perturb_derivs is always
         *  called at the end of generic_evolver, in order to update
         *  all quantities in ppw to the time at which the
         *  approximation is switched off) */
        ppv->y[ppv->index_pt_shear_g] = ppw->tca_shear_g;

        ppv->y[ppv->index_pt_l3_g] = 6. / 7. * k / ppw->pvecthermo[pth->index_th_dkappa] * ppw->s_l[3] * ppv->y[ppv->index_pt_shear_g]; /* second-order tight-coupling approximation for l=3 */

        ppv->y[ppv->index_pt_pol0_g] = 2.5 * ppv->y[ppv->index_pt_shear_g];                                                                      /* first-order tight-coupling approximation for polarization, l=0 */
        ppv->y[ppv->index_pt_pol1_g] = k / ppw->pvecthermo[pth->index_th_dkappa] * (5. - 2. * ppw->s_l[2]) / 6. * ppv->y[ppv->index_pt_shear_g]; /* second-order tight-coupling approximation for polarization, l=1 */
        ppv->y[ppv->index_pt_pol2_g] = 0.5 * ppv->y[ppv->index_pt_shear_g];                                                                      /* first-order tight-coupling approximation for polarization, l=2 */
        ppv->y[ppv->index_pt_pol3_g] = k / ppw->pvecthermo[pth->index_th_dkappa] * 3. * ppw->s_l[3] / 14. * ppv->y[ppv->index_pt_shear_g];       /* second-order tight-coupling approximation for polarization, l=3 */

        if (pba->has_ur == _TRUE_)
        {
          ppv->y[ppv->index_pt_delta_ur] =
            ppw->pv->y[ppw->pv->index_pt_delta_ur];

          ppv->y[ppv->index_pt_theta_ur] =
            ppw->pv->y[ppw->pv->index_pt_theta_ur];

          ppv->y[ppv->index_pt_shear_ur] =
            ppw->pv->y[ppw->pv->index_pt_shear_ur];

          if (ppw->approx[ppw->index_ap_ufa] == (int) ufa_off)
          {
            ppv->y[ppv->index_pt_l3_ur] =
              ppw->pv->y[ppw->pv->index_pt_l3_ur];

            for (l = 4; l <= ppv->l_max_ur; l++)
              ppv->y[ppv->index_pt_delta_ur + l] =
                ppw->pv->y[ppw->pv->index_pt_delta_ur + l];
          }
        }

        if (pba->has_ncdm == _TRUE_)
        {
          index_pt = 0;

          for (n_ncdm = 0; n_ncdm < ppv->N_ncdm; n_ncdm++)
          {
            for (index_q = 0; index_q < ppv->q_size_ncdm[n_ncdm]; index_q++)
            {
              for (l = 0; l <= ppv->l_max_ncdm[n_ncdm]; l++)
              {
                /* This is correct with or without ncdmfa, since ppv->lmax_ncdm is set accordingly. */
                ppv->y[ppv->index_pt_psi0_ncdm1 + index_pt] =
                  ppw->pv->y[ppw->pv->index_pt_psi0_ncdm1 + index_pt];
                index_pt++;
              }
            }
          }
        }

        /* perturbed recombination */
        /* the initial conditions are set when tca is switched off (current block) */
        if (ppt->has_perturbed_recombination == _TRUE_)
        {
          ppv->y[ppv->index_pt_perturbed_recombination_delta_temp] = 1. / 3. * ppv->y[ppw->pv->index_pt_delta_b];
          ppv->y[ppv->index_pt_perturbed_recombination_delta_chi]  = 0.;
        }
      } /* end of block tca ON -> tca OFF */

      /* perturbed recombination */
      /* For any other transition in the approximation scheme, we should just copy the value of the perturbations, provided tca is already off (otherwise the indices are not yet allocated). For instance, we do not want to copy the values in the (k,tau) region where both UFA and TCA are engaged.*/

      if ((ppt->has_perturbed_recombination == _TRUE_) && (pa_old[ppw->index_ap_tca] == (int) tca_off))
      {
        ppv->y[ppv->index_pt_perturbed_recombination_delta_temp] =
          ppw->pv->y[ppw->pv->index_pt_perturbed_recombination_delta_temp];
        ppv->y[ppv->index_pt_perturbed_recombination_delta_chi] =
          ppw->pv->y[ppw->pv->index_pt_perturbed_recombination_delta_chi];
      }


      /* -- case of switching on radiation streaming
       *  approximation. Provide correct initial conditions to new set
       *  of variables */

      if ((pa_old[ppw->index_ap_rsa] == (int) rsa_off) && (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on))
      {
        if (ppt->perturbations_verbose > 2)
          fprintf (stdout, "Mode k=%e: switch on radiation streaming approximation at tau=%e with Omega_r=%g\n", k, tau, ppw->pvecback[pba->index_bg_Omega_r]);

        if (pba->has_ncdm == _TRUE_)
        {
          index_pt = 0;

          for (n_ncdm = 0; n_ncdm < ppv->N_ncdm; n_ncdm++)
          {
            for (index_q = 0; index_q < ppv->q_size_ncdm[n_ncdm]; index_q++)
            {
              for (l = 0; l <= ppv->l_max_ncdm[n_ncdm]; l++)
              {
                ppv->y[ppv->index_pt_psi0_ncdm1 + index_pt] =
                  ppw->pv->y[ppw->pv->index_pt_psi0_ncdm1 + index_pt];
                index_pt++;
              }
            }
          }
        }
      }

      /* -- case of switching on ur fluid
       *  approximation. Provide correct initial conditions to new set
       *  of variables */

      if (pba->has_ur == _TRUE_)
      {
        if ((pa_old[ppw->index_ap_ufa] == (int) ufa_off) && (ppw->approx[ppw->index_ap_ufa] == (int) ufa_on))
        {
          if (ppt->perturbations_verbose > 2)
            fprintf (stdout, "Mode k=%e: switch on ur fluid approximation at tau=%e\n", k, tau);

          if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
          {
            ppv->y[ppv->index_pt_delta_g] =
              ppw->pv->y[ppw->pv->index_pt_delta_g];

            ppv->y[ppv->index_pt_theta_g] =
              ppw->pv->y[ppw->pv->index_pt_theta_g];
          }

          if ((ppw->approx[ppw->index_ap_tca] == (int) tca_off) && (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off))
          {
            ppv->y[ppv->index_pt_shear_g] =
              ppw->pv->y[ppw->pv->index_pt_shear_g];

            ppv->y[ppv->index_pt_l3_g] =
              ppw->pv->y[ppw->pv->index_pt_l3_g];

            for (l = 4; l <= ppw->pv->l_max_g; l++)
            {
              ppv->y[ppv->index_pt_delta_g + l] =
                ppw->pv->y[ppw->pv->index_pt_delta_g + l];
            }

            ppv->y[ppv->index_pt_pol0_g] =
              ppw->pv->y[ppw->pv->index_pt_pol0_g];

            ppv->y[ppv->index_pt_pol1_g] =
              ppw->pv->y[ppw->pv->index_pt_pol1_g];

            ppv->y[ppv->index_pt_pol2_g] =
              ppw->pv->y[ppw->pv->index_pt_pol2_g];

            ppv->y[ppv->index_pt_pol3_g] =
              ppw->pv->y[ppw->pv->index_pt_pol3_g];

            for (l = 4; l <= ppw->pv->l_max_pol_g; l++)
            {
              ppv->y[ppv->index_pt_pol0_g + l] =
                ppw->pv->y[ppw->pv->index_pt_pol0_g + l];
            }
          }

          if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
          {
            ppv->y[ppv->index_pt_delta_ur] =
              ppw->pv->y[ppw->pv->index_pt_delta_ur];

            ppv->y[ppv->index_pt_theta_ur] =
              ppw->pv->y[ppw->pv->index_pt_theta_ur];

            ppv->y[ppv->index_pt_shear_ur] =
              ppw->pv->y[ppw->pv->index_pt_shear_ur];
          }

          if (pba->has_ncdm == _TRUE_)
          {
            index_pt = 0;

            for (n_ncdm = 0; n_ncdm < ppv->N_ncdm; n_ncdm++)
            {
              for (index_q = 0; index_q < ppv->q_size_ncdm[n_ncdm]; index_q++)
              {
                for (l = 0; l <= ppv->l_max_ncdm[n_ncdm]; l++)
                {
                  /* This is correct even when ncdmfa == off, since ppv->l_max_ncdm and
                   *   ppv->q_size_ncdm is updated.*/
                  ppv->y[ppv->index_pt_psi0_ncdm1 + index_pt] =
                    ppw->pv->y[ppw->pv->index_pt_psi0_ncdm1 + index_pt];
                  index_pt++;
                }
              }
            }
          }
        }
      }

      /* -- case of switching on ncdm fluid
       *  approximation. Provide correct initial conditions to new set
       *  of variables */

      if (pba->has_ncdm == _TRUE_)
      {
        if ((pa_old[ppw->index_ap_ncdmfa] == (int) ncdmfa_off) && (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on))
        {
          if (ppt->perturbations_verbose > 2)
            fprintf (stdout, "Mode k=%e: switch on ncdm fluid approximation at tau=%e\n", k, tau);

          if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
          {
            ppv->y[ppv->index_pt_delta_g] =
              ppw->pv->y[ppw->pv->index_pt_delta_g];

            ppv->y[ppv->index_pt_theta_g] =
              ppw->pv->y[ppw->pv->index_pt_theta_g];
          }

          if ((ppw->approx[ppw->index_ap_tca] == (int) tca_off) && (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off))
          {
            ppv->y[ppv->index_pt_shear_g] =
              ppw->pv->y[ppw->pv->index_pt_shear_g];

            ppv->y[ppv->index_pt_l3_g] =
              ppw->pv->y[ppw->pv->index_pt_l3_g];

            for (l = 4; l <= ppw->pv->l_max_g; l++)
            {
              ppv->y[ppv->index_pt_delta_g + l] =
                ppw->pv->y[ppw->pv->index_pt_delta_g + l];
            }

            ppv->y[ppv->index_pt_pol0_g] =
              ppw->pv->y[ppw->pv->index_pt_pol0_g];

            ppv->y[ppv->index_pt_pol1_g] =
              ppw->pv->y[ppw->pv->index_pt_pol1_g];

            ppv->y[ppv->index_pt_pol2_g] =
              ppw->pv->y[ppw->pv->index_pt_pol2_g];

            ppv->y[ppv->index_pt_pol3_g] =
              ppw->pv->y[ppw->pv->index_pt_pol3_g];

            for (l = 4; l <= ppw->pv->l_max_pol_g; l++)
            {
              ppv->y[ppv->index_pt_pol0_g + l] =
                ppw->pv->y[ppw->pv->index_pt_pol0_g + l];
            }
          }

          if (pba->has_ur == _TRUE_)
          {
            if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
            {
              ppv->y[ppv->index_pt_delta_ur] =
                ppw->pv->y[ppw->pv->index_pt_delta_ur];

              ppv->y[ppv->index_pt_theta_ur] =
                ppw->pv->y[ppw->pv->index_pt_theta_ur];

              ppv->y[ppv->index_pt_shear_ur] =
                ppw->pv->y[ppw->pv->index_pt_shear_ur];

              if (ppw->approx[ppw->index_ap_ufa] == (int) ufa_off)
              {
                ppv->y[ppv->index_pt_l3_ur] =
                  ppw->pv->y[ppw->pv->index_pt_l3_ur];

                for (l = 4; l <= ppv->l_max_ur; l++)
                  ppv->y[ppv->index_pt_delta_ur + l] =
                    ppw->pv->y[ppw->pv->index_pt_delta_ur + l];
              }
            }
          }

          a        = ppw->pvecback[pba->index_bg_a];
          index_pt = ppw->pv->index_pt_psi0_ncdm1;

          for (n_ncdm = 0; n_ncdm < ppv->N_ncdm; n_ncdm++)
          {
            /* We are in the fluid approximation, so ncdm_l_size is always 3. */
            ncdm_l_size     = ppv->l_max_ncdm[n_ncdm] + 1;
            rho_plus_p_ncdm = ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm] +
                              ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm];

            for (l = 0; l <= 2; l++)
            {
              ppv->y[ppv->index_pt_psi0_ncdm1 + ncdm_l_size * n_ncdm + l] = 0.0;
            }

            factor = pba->factor_ncdm[n_ncdm] * pow (pba->a_today / a, 4);

            for (index_q = 0; index_q < ppw->pv->q_size_ncdm[n_ncdm]; index_q++)
            {
              /* Integrate over distributions: */
              q                                                        = pba->q_ncdm[n_ncdm][index_q];
              q2                                                       = q * q;
              epsilon                                                  = sqrt (q2 + a * a * pba->M_ncdm[n_ncdm] * pba->M_ncdm[n_ncdm]);
              ppv->y[ppv->index_pt_psi0_ncdm1 + ncdm_l_size * n_ncdm] +=
                pba->w_ncdm[n_ncdm][index_q] * q2 * epsilon *
                ppw->pv->y[index_pt];

              ppv->y[ppv->index_pt_psi0_ncdm1 + ncdm_l_size * n_ncdm + 1] +=
                pba->w_ncdm[n_ncdm][index_q] * q2 * q *
                ppw->pv->y[index_pt + 1];

              ppv->y[ppv->index_pt_psi0_ncdm1 + ncdm_l_size * n_ncdm + 2] +=
                pba->w_ncdm[n_ncdm][index_q] * q2 * q2 / epsilon *
                ppw->pv->y[index_pt + 2];

              /*Jump to next momentum bin in ppw->pv->y: */
              index_pt += (ppw->pv->l_max_ncdm[n_ncdm] + 1);
            }

            ppv->y[ppv->index_pt_psi0_ncdm1 + ncdm_l_size * n_ncdm]     *= factor / ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm];
            ppv->y[ppv->index_pt_psi0_ncdm1 + ncdm_l_size * n_ncdm + 1] *= k * factor / rho_plus_p_ncdm;
            ppv->y[ppv->index_pt_psi0_ncdm1 + ncdm_l_size * n_ncdm + 2] *= 2.0 / 3.0 * factor / rho_plus_p_ncdm;
          }
        }
      }
    }

    /** - --> (b) for the vector mode */

    if (_vectors_)
    {
      /** - ---> (b.1.) check that the change of approximation scheme makes
       *   sense (note: before calling this routine there is already a
       *   check that we wish to change only one approximation flag at
       *   a time) */

      class_test ((pa_old[ppw->index_ap_tca] == (int) tca_off) && (ppw->approx[ppw->index_ap_tca] == (int) tca_on),
                  ppt->error_message,
                  "at tau=%g: the tight-coupling approximation can be switched off, not on", tau);

      /** - ---> (b.2.) some variables (gw, gwdot, ...) are not affected by
       *   any approximation. They need to be reconducted whatever
       *   the approximation switching is. We treat them here. Below
       *   we will treat other variables case by case. */

      if (ppt->gauge == synchronous)
        ppv->y[ppv->index_pt_hv_prime] =
          ppw->pv->y[ppw->pv->index_pt_hv_prime];



      if (ppt->gauge == newtonian)
        ppv->y[ppv->index_pt_V] =
          ppw->pv->y[ppw->pv->index_pt_V];


      ppv->y[ppv->index_pt_theta_b] =
        ppw->pv->y[ppw->pv->index_pt_theta_b];


      /* -- case of switching off tight coupling
       *  approximation. Provide correct initial conditions to new set
       *  of variables */

      if ((pa_old[ppw->index_ap_tca] == (int) tca_on) && (ppw->approx[ppw->index_ap_tca] == (int) tca_off))
      {
        if (ppt->perturbations_verbose > 2)
          fprintf (stdout, "Mode k=%e: switch off tight-coupling approximation at tau=%e\n", k, tau);

        ppv->y[ppv->index_pt_delta_g] = 0.0; /*TBC */
        /*-4./3.*ppw->pv->y[ppw->pv->index_pt_gwdot]/ppw->pvecthermo[pth->index_th_dkappa]; */

        ppv->y[ppv->index_pt_pol0_g] = 0.0; /*TBC */
        /*1./3.*ppw->pv->y[ppw->pv->index_pt_gwdot]/ppw->pvecthermo[pth->index_th_dkappa]; */
      }

      /* -- case of switching on radiation streaming
       *  approximation. Provide correct initial conditions to new set
       *  of variables */

      if ((pa_old[ppw->index_ap_rsa] == (int) rsa_off) && (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on))
        if (ppt->perturbations_verbose > 2)
          fprintf (stdout, "Mode k=%e: switch on radiation streaming approximation at tau=%e with Omega_r=%g\n", k, tau, ppw->pvecback[pba->index_bg_Omega_r]);
    }

    /** - --> (c) for the tensor mode */

    if (_tensors_)
    {
      /** - ---> (c.1.) check that the change of approximation scheme makes
       *   sense (note: before calling this routine there is already a
       *   check that we wish to change only one approximation flag at
       *   a time) */

      class_test ((pa_old[ppw->index_ap_tca] == (int) tca_off) && (ppw->approx[ppw->index_ap_tca] == (int) tca_on),
                  ppt->error_message,
                  "at tau=%g: the tight-coupling approximation can be switched off, not on", tau);

      /** - ---> (c.2.) some variables (gw, gwdot, ...) are not affected by
       *   any approximation. They need to be reconducted whatever
       *   the approximation switching is. We treat them here. Below
       *   we will treat other variables case by case. */


      ppv->y[ppv->index_pt_gw] =
        ppw->pv->y[ppw->pv->index_pt_gw];

      ppv->y[ppv->index_pt_gwdot] =
        ppw->pv->y[ppw->pv->index_pt_gwdot];

      if (ppt->evolve_tensor_ur == _TRUE_)
      {
        /* For now, neutrinos go here. */
        ppv->y[ppv->index_pt_delta_ur] =
          ppw->pv->y[ppw->pv->index_pt_delta_ur];

        ppv->y[ppv->index_pt_theta_ur] =
          ppw->pv->y[ppw->pv->index_pt_theta_ur];

        ppv->y[ppv->index_pt_shear_ur] =
          ppw->pv->y[ppw->pv->index_pt_shear_ur];

        ppv->y[ppv->index_pt_l3_ur] =
          ppw->pv->y[ppw->pv->index_pt_l3_ur];

        for (l = 4; l <= ppv->l_max_ur; l++)
          ppv->y[ppv->index_pt_delta_ur + l] =
            ppw->pv->y[ppw->pv->index_pt_delta_ur + l];
      }

      if (ppt->evolve_tensor_ncdm == _TRUE_)
      {
        index_pt = 0;

        for (n_ncdm = 0; n_ncdm < ppv->N_ncdm; n_ncdm++)
        {
          for (index_q = 0; index_q < ppv->q_size_ncdm[n_ncdm]; index_q++)
          {
            for (l = 0; l <= ppv->l_max_ncdm[n_ncdm]; l++)
            {
              /* This is correct with or without ncdmfa, since ppv->lmax_ncdm is set accordingly. */
              ppv->y[ppv->index_pt_psi0_ncdm1 + index_pt] =
                ppw->pv->y[ppw->pv->index_pt_psi0_ncdm1 + index_pt];
              index_pt++;
            }
          }
        }
      }

      /* -- case of switching off tight coupling
       *  approximation. Provide correct initial conditions to new set
       *  of variables */

      if ((pa_old[ppw->index_ap_tca] == (int) tca_on) && (ppw->approx[ppw->index_ap_tca] == (int) tca_off))
      {
        if (ppt->perturbations_verbose > 2)
          fprintf (stdout, "Mode k=%e: switch off tight-coupling approximation at tau=%e\n", k, tau);

        ppv->y[ppv->index_pt_delta_g] = -4. / 3. * ppw->pv->y[ppw->pv->index_pt_gwdot] / ppw->pvecthermo[pth->index_th_dkappa];

        ppv->y[ppv->index_pt_pol0_g] = 1. / 3. * ppw->pv->y[ppw->pv->index_pt_gwdot] / ppw->pvecthermo[pth->index_th_dkappa];
      }

      /* -- case of switching on radiation streaming
       *  approximation. Provide correct initial conditions to new set
       *  of variables */

      if ((pa_old[ppw->index_ap_rsa] == (int) rsa_off) && (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on))
        if (ppt->perturbations_verbose > 2)
          fprintf (stdout, "Mode k=%e: switch on radiation streaming approximation at tau=%e with Omega_r=%g\n", k, tau, ppw->pvecback[pba->index_bg_Omega_r]);
    }

    /** - --> (d) free the previous vector of perturbations */

    class_call (perturb_vector_free (ppw->pv),
                ppt->error_message,
                ppt->error_message);

    /** - --> (e) let ppw-->pv points towards the perturb_vector structure
     *   that we just created */

    ppw->pv = ppv;
  }

  return _SUCCESS_;
}

/**
 * Free the perturb_vector structure.
 *
 * @param pv        Input: pointer to perturb_vector structure to be freed
 * @return the error status
 */

int
perturb_vector_free (
  struct perturb_vector *pv
                    )
{
  if (pv->l_max_ncdm != NULL)
    free (pv->l_max_ncdm);

  if (pv->q_size_ncdm != NULL)
    free (pv->q_size_ncdm);

  free (pv->y);
  free (pv->dy);
  free (pv->used_in_sources);
  free (pv);

  return _SUCCESS_;
}

/**
 * For each mode, wavenumber and initial condition, this function
 * initializes in the vector all values of perturbed variables (in a
 * given gauge). It is assumed here that all values have previously been
 * set to zero, only non-zero values are set here.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_md   Input: index of mode under consideration (scalar/.../tensor)
 * @param index_ic   Input: index of initial condition under consideration (ad, iso...)
 * @param k          Input: wavenumber
 * @param tau        Input: conformal time
 * @param ppw        Input/Output: workspace containing in input the approximation scheme, the background/thermodynamics/metric quantities, and eventually the previous vector y; and in output the new vector y.
 * @return the error status
 */

int
perturb_initial_conditions (struct precision         *ppr,
                            const struct background  *pba,
                            struct perturbs          *ppt,
                            int                      index_md,
                            int                      index_ic,
                            double                   k,
                            double                   tau,
                            struct perturb_workspace *ppw
                           )
{
  /** Summary: */

  /** --> Declare local variables */
  /* const double H0 = 1.0 / nc_hicosmo_RH_Mpc (pba->cosmo); */
  int sgnK;
  double K = nc_hicosmo_curvature_K (pba->cosmo, &sgnK);

  double a, a_prime_over_a;
  double w_fld, dw_over_da_fld, integral_fld;
  double delta_ur = 0., theta_ur = 0., shear_ur = 0., l3_ur = 0., eta = 0., delta_cdm = 0., alpha;
  double q, epsilon, k2;
  int index_q, n_ncdm, idx;
  double rho_r, rho_m, rho_nu, rho_m_over_rho_r;
  double fracnu, fracg, fracb, fraccdm, om;
  double ktau_two, ktau_three;

  double delta_tot;
  double velocity_tot;
  double s2_squared;

  /** --> For scalars */

  if (_scalars_)
  {
    /** - (a) compute relevant background quantities: compute rho_r,
     *   rho_m, rho_nu (= all relativistic except photons), and their
     *   ratio. */

    class_call (background_at_tau (pba,
                                   tau,
                                   pba->normal_info,
                                   pba->inter_normal,
                                   &(ppw->last_index_back),
                                   ppw->pvecback),
                pba->error_message,
                ppt->error_message);

    a = ppw->pvecback[pba->index_bg_a];

    a_prime_over_a = ppw->pvecback[pba->index_bg_H] * a;

    /* 8piG/3 rho_r(t_i) */
    rho_r = ppw->pvecback[pba->index_bg_rho_g];

    /* 8piG/3 rho_m(t_i) */
    rho_m = ppw->pvecback[pba->index_bg_rho_b];

    /* 8piG/3 rho_nu(t_i) (all neutrinos and collisionless relics being relativistic at that time) */
    rho_nu = 0.;

    if (pba->has_cdm == _TRUE_)
      rho_m += ppw->pvecback[pba->index_bg_rho_cdm];

    if (pba->has_ur == _TRUE_)
    {
      rho_r  += ppw->pvecback[pba->index_bg_rho_ur];
      rho_nu += ppw->pvecback[pba->index_bg_rho_ur];
    }

    if (pba->has_ncdm == _TRUE_)
    {
      for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
      {
        rho_r  += ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm];
        rho_nu += ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm];
      }
    }

    class_test (rho_r == 0.,
                ppt->error_message,
                "stop to avoid division by zero");

    /* f_nu = Omega_nu(t_i) / Omega_r(t_i) */
    fracnu = rho_nu / rho_r;

    /* f_g = Omega_g(t_i) / Omega_r(t_i) */
    fracg = ppw->pvecback[pba->index_bg_rho_g] / rho_r;

    /* f_b = Omega_b(t_i) / Omega_m(t_i) */
    fracb = ppw->pvecback[pba->index_bg_rho_b] / rho_m;

    /* f_cdm = Omega_cdm(t_i) / Omega_m(t_i) */
    fraccdm = 1. - fracb;

    /* Omega_m(t_i) / Omega_r(t_i) */
    rho_m_over_rho_r = rho_m / rho_r;

    /* omega = Omega_m(t_i) a(t_i) H(t_i) / sqrt(Omega_r(t_i))
     *  = Omega_m(t_0) a(t_0) H(t_0) / sqrt(Omega_r(t_0)) assuming rho_m in a-3 and rho_r in a^-4
     *  = (8piG/3 rho_m(t_i)) a(t_i) / sqrt(8piG/3 rho_r(t_i))  in Mpc-1
     *  This (a priori strange) parameter is the relevant one for expressing a
     *  as a function of tau during radiation and matter domination (but not DE domination).
     *  Indeed the exact solution of Friedmann when there is only radiation and matter in
     *  the universe is
     *  a = [H(t_0)^2 Omega_m(t_0) a(t_0)^3 / 4] x [tau^2 + 4 tau / omega]
     */
    om = a * rho_m / sqrt (rho_r);

    /* (k tau)^2, (k tau)^3 */
    ktau_two   = k * k * tau * tau;
    ktau_three = k * tau * ktau_two;


    /* curvature-dependent factors */

    s2_squared = 1. - 3. * K / k / k;

    /** - (b) starts by setting everything in synchronous gauge. If
     *   another gauge is needed, we will perform a gauge
     *   transformation below. */

    /** - --> (b.1.) adiabatic */

    if ((ppt->has_ad == _TRUE_) && (index_ic == ppt->index_ic_ad))
    {
      /* The following formulas are valid at leading order in
       *  (k*tau) and (om*tau), and order zero in
       *  tight-coupling. Identical to first order terms in CRS,
       *  except for normalization (when ppr->curvature_ini=1, tau=1:
       *  leads to factor 1/2 difference between CRS formulas with
       *  beta1=0). Identical to CAMB when om set to zero in theta_g,
       *  theta_ur, shear_ur, tau
       *
       *  In the non-flat case the relation R=eta is still valid
       *  outside the horizon for adiabatic IC. Hence eta is still
       *  set to ppr->curvature_ini at leading order.  Factors s2
       *  appear through the solution of Einstein equations and
       *  equations of motion. */

      /* photon density */
      ppw->pv->y[ppw->pv->index_pt_delta_g] = -ktau_two / 3. * (1. - om * tau / 5.)
                                              * ppr->curvature_ini * s2_squared;

      /* photon velocity */
      ppw->pv->y[ppw->pv->index_pt_theta_g] = -k * ktau_three / 36. * (1. - 3. * (1. + 5. * fracb - fracnu) / 20. / (1. - fracnu) * om * tau)
                                              * ppr->curvature_ini * s2_squared;

      /* tighly-coupled baryons */
      ppw->pv->y[ppw->pv->index_pt_delta_b] = 3. / 4. * ppw->pv->y[ppw->pv->index_pt_delta_g]; /* baryon density */
      ppw->pv->y[ppw->pv->index_pt_theta_b] = ppw->pv->y[ppw->pv->index_pt_theta_g];           /* baryon velocity */

      if (pba->has_cdm == _TRUE_)
        ppw->pv->y[ppw->pv->index_pt_delta_cdm] = 3. / 4. * ppw->pv->y[ppw->pv->index_pt_delta_g];  /* cdm density */

      /* cdm velocity vanishes in the synchronous gauge */

      /* fluid (assumes wa=0, if this is not the case the
       *  fluid will catch anyway the attractor solution) */
      if (pba->has_fld == _TRUE_)
      {
        class_call (background_w_fld (pba, a, &w_fld, &dw_over_da_fld, &integral_fld), pba->error_message, ppt->error_message);

        if (pba->use_ppf == _FALSE_)
        {
          ppw->pv->y[ppw->pv->index_pt_delta_fld] = -ktau_two / 4. * (1. + w_fld) * (4. - 3. * pba->cs2_fld) / (4. - 6. * w_fld + 3. * pba->cs2_fld) * ppr->curvature_ini * s2_squared; /* from 1004.5509 */ /*TBC: curvature */

          ppw->pv->y[ppw->pv->index_pt_theta_fld] = -k * ktau_three / 4. * pba->cs2_fld / (4. - 6. * w_fld + 3. * pba->cs2_fld) * ppr->curvature_ini * s2_squared; /* from 1004.5509 */ /*TBC:curvature */
        }

        /* if use_ppf == _TRUE_, y[ppw->pv->index_pt_Gamma_fld] will be automatically set to zero, and this is what we want (although one could probably work out some small nonzero initial conditions: TODO) */
      }

      /* all relativistic relics: ur, early ncdm, dr */

      if ((pba->has_ur == _TRUE_) || (pba->has_ncdm == _TRUE_))
      {
        delta_ur = ppw->pv->y[ppw->pv->index_pt_delta_g]; /* density of ultra-relativistic neutrinos/relics */

        theta_ur = -k * ktau_three / 36. / (4. * fracnu + 15.) * (4. * fracnu + 11. + 12. * s2_squared - 3. * (8. * fracnu * fracnu + 50. * fracnu + 275.) / 20. / (2. * fracnu + 15.) * tau * om) * ppr->curvature_ini * s2_squared; /* velocity of ultra-relativistic neutrinos/relics */ /*TBC */

        shear_ur = ktau_two / (45. + 12. * fracnu) * (3. * s2_squared - 1.) * (1. + (4. * fracnu - 5.) / 4. / (2. * fracnu + 15.) * tau * om) * ppr->curvature_ini; /*TBC /s2_squared; / * shear of ultra-relativistic neutrinos/relics * /  //TBC:0 */

        l3_ur = ktau_three * 2. / 7. / (12. * fracnu + 45.) * ppr->curvature_ini; /*TBC */
      }

      /* synchronous metric perturbation eta */
      /*eta = ppr->curvature_ini * (1.-ktau_two/12./(15.+4.*fracnu)*(5.+4.*fracnu - (16.*fracnu*fracnu+280.*fracnu+325)/10./(2.*fracnu+15.)*tau*om)) /  s2_squared; */
      /*eta = ppr->curvature_ini * s2_squared * (1.-ktau_two/12./(15.+4.*fracnu)*(15.*s2_squared-10.+4.*s2_squared*fracnu - (16.*fracnu*fracnu+280.*fracnu+325)/10./(2.*fracnu+15.)*tau*om)); */
      eta = ppr->curvature_ini * (1. - ktau_two / 12. / (15. + 4. * fracnu) * (5. + 4. * s2_squared * fracnu - (16. * fracnu * fracnu + 280. * fracnu + 325) / 10. / (2. * fracnu + 15.) * tau * om));
    }

    /* isocurvature initial conditions taken from Bucher, Moodely,
     *  Turok 99, with just a different normalization convention for
     *  tau and the scale factor. [k tau] from BMT99 is left invariant
     *  because it is the ratio [k/aH]. But [Omega_i,0 tau] from BMT99
     *  must be replaced by [frac_i*om*tau/4]. Some doubts remain about
     *  the niv formulas, that should be recheked at some point. We
     *  also checked that for bi,cdi,nid, everything coincides exactly
     *  with the CAMB formulas. */

    /** - --> (b.2.) Cold dark matter Isocurvature */

    if ((ppt->has_cdi == _TRUE_) && (index_ic == ppt->index_ic_cdi))
    {
      class_test (pba->has_cdm == _FALSE_,
                  ppt->error_message,
                  "not consistent to ask for CDI in absence of CDM!");

      ppw->pv->y[ppw->pv->index_pt_delta_g] = ppr->entropy_ini * fraccdm * om * tau * (-2. / 3. + om * tau / 4.);
      ppw->pv->y[ppw->pv->index_pt_theta_g] = -ppr->entropy_ini * fraccdm * om * ktau_two / 12.;

      ppw->pv->y[ppw->pv->index_pt_delta_b] = 3. / 4. * ppw->pv->y[ppw->pv->index_pt_delta_g];
      ppw->pv->y[ppw->pv->index_pt_theta_b] = ppw->pv->y[ppw->pv->index_pt_theta_g];

      ppw->pv->y[ppw->pv->index_pt_delta_cdm] = ppr->entropy_ini + 3. / 4. * ppw->pv->y[ppw->pv->index_pt_delta_g];

      if ((pba->has_ur == _TRUE_) || (pba->has_ncdm == _TRUE_))
      {
        delta_ur = ppw->pv->y[ppw->pv->index_pt_delta_g];
        theta_ur = ppw->pv->y[ppw->pv->index_pt_theta_g];
        shear_ur = -ppr->entropy_ini * fraccdm * ktau_two * tau * om / 6. / (2. * fracnu + 15.);
      }

      eta = -ppr->entropy_ini * fraccdm * om * tau * (1. / 6. - om * tau / 16.);
    }

    /** - --> (b.3.) Baryon Isocurvature */

    if ((ppt->has_bi == _TRUE_) && (index_ic == ppt->index_ic_bi))
    {
      ppw->pv->y[ppw->pv->index_pt_delta_g] = ppr->entropy_ini * fracb * om * tau * (-2. / 3. + om * tau / 4.);
      ppw->pv->y[ppw->pv->index_pt_theta_g] = -ppr->entropy_ini * fracb * om * ktau_two / 12.;

      ppw->pv->y[ppw->pv->index_pt_delta_b] = ppr->entropy_ini + 3. / 4. * ppw->pv->y[ppw->pv->index_pt_delta_g];
      ppw->pv->y[ppw->pv->index_pt_theta_b] = ppw->pv->y[ppw->pv->index_pt_theta_g];

      if (pba->has_cdm == _TRUE_)
        ppw->pv->y[ppw->pv->index_pt_delta_cdm] = 3. / 4. * ppw->pv->y[ppw->pv->index_pt_delta_g];


      if ((pba->has_ur == _TRUE_) || (pba->has_ncdm == _TRUE_))
      {
        delta_ur = ppw->pv->y[ppw->pv->index_pt_delta_g];
        theta_ur = ppw->pv->y[ppw->pv->index_pt_theta_g];
        shear_ur = -ppr->entropy_ini * fracb * ktau_two * tau * om / 6. / (2. * fracnu + 15.);
      }

      eta = -ppr->entropy_ini * fracb * om * tau * (1. / 6. - om * tau / 16.);
    }

    /** - --> (b.4.) Neutrino density Isocurvature */

    if ((ppt->has_nid == _TRUE_) && (index_ic == ppt->index_ic_nid))
    {
      class_test ((pba->has_ur == _FALSE_) && (pba->has_ncdm == _FALSE_),
                  ppt->error_message,
                  "not consistent to ask for NID in absence of ur or ncdm species!");

      ppw->pv->y[ppw->pv->index_pt_delta_g] = ppr->entropy_ini * fracnu / fracg * (-1. + ktau_two / 6.);
      ppw->pv->y[ppw->pv->index_pt_theta_g] = -ppr->entropy_ini * fracnu / fracg * k * k * tau * (1. / 4. - fracb / fracg * 3. / 16. * om * tau);

      ppw->pv->y[ppw->pv->index_pt_delta_b] = ppr->entropy_ini * fracnu / fracg / 8. * ktau_two;
      ppw->pv->y[ppw->pv->index_pt_theta_b] = ppw->pv->y[ppw->pv->index_pt_theta_g];

      if (pba->has_cdm == _TRUE_)
        ppw->pv->y[ppw->pv->index_pt_delta_cdm] = -ppr->entropy_ini * fracnu * fracb / fracg / 80. * ktau_two * om * tau;


      delta_ur = ppr->entropy_ini * (1. - ktau_two / 6.);
      theta_ur = ppr->entropy_ini * k * k * tau / 4.;
      shear_ur = ppr->entropy_ini * ktau_two / (4. * fracnu + 15.) / 2.;

      eta = -ppr->entropy_ini * fracnu / (4. * fracnu + 15.) / 6. * ktau_two;
    }

    /** - --> (b.5.) Neutrino velocity Isocurvature */

    if ((ppt->has_niv == _TRUE_) && (index_ic == ppt->index_ic_niv))
    {
      class_test ((pba->has_ur == _FALSE_) && (pba->has_ncdm == _FALSE_),
                  ppt->error_message,
                  "not consistent to ask for NIV in absence of ur or ncdm species!");

      ppw->pv->y[ppw->pv->index_pt_delta_g] = ppr->entropy_ini * k * tau * fracnu / fracg *
                                              (1. - 3. / 16. * fracb * (2. + fracg) / fracg * om * tau); /* small diff wrt camb */

      ppw->pv->y[ppw->pv->index_pt_theta_g] = ppr->entropy_ini * fracnu / fracg * 3. / 4. * k *
                                              (-1. + 3. / 4. * fracb / fracg * om * tau + 3. / 16. * om * om * tau * tau * fracb / fracg / fracg * (fracg - 3. * fracb) + ktau_two / 6.);

      ppw->pv->y[ppw->pv->index_pt_delta_b] = 3. / 4. * ppw->pv->y[ppw->pv->index_pt_delta_g]; /* small diff wrt camb */
      ppw->pv->y[ppw->pv->index_pt_theta_b] = ppw->pv->y[ppw->pv->index_pt_theta_g];

      if (pba->has_cdm == _TRUE_)
        ppw->pv->y[ppw->pv->index_pt_delta_cdm] = -ppr->entropy_ini * 9. / 64. * fracnu * fracb / fracg * k * tau * om * tau;


      delta_ur = -ppr->entropy_ini * k * tau * (1. + 3. / 16. * fracb * fracnu / fracg * om * tau); /* small diff wrt camb */
      theta_ur = ppr->entropy_ini * 3. / 4. * k * (1. - 1. / 6. * ktau_two * (4. * fracnu + 9.) / (4. * fracnu + 5.));
      shear_ur = ppr->entropy_ini / (4. * fracnu + 15.) * k * tau * (1. + 3. * om * tau * fracnu / (4. * fracnu + 15.)); /* small diff wrt camb */

      eta = ppr->entropy_ini * fracnu * k * tau * (-1. / (4. * fracnu + 5.) + (-3. / 64. * fracb / fracg + 15. / 4. / (4. * fracnu + 15.) / (4. * fracnu + 5.) * om * tau)); /* small diff wrt camb */
    }

    /** - (c) If the needed gauge is really the synchronous gauge, we need to affect the previously computed value of eta to the actual variable eta */

    if (ppt->gauge == synchronous)
      ppw->pv->y[ppw->pv->index_pt_eta] = eta;


    /** - (d) If the needed gauge is the newtonian gauge, we must compute alpha and then perform a gauge transformation for each variable */

    if (ppt->gauge == newtonian)
    {
      /* alpha is like in Ma & Bertschinger: (h'+6 eta')/(2k^2). We obtain it from the first two Einstein equations:
       *
       *  alpha = [eta + 3/2 (a'/a)^2 (delta_rho/rho_c) / k^2 /s_2^2 + 3/2 (a'/a)^3 3 ((rho+p)theta/rho_c) / k^4 / s_2^2] / (a'/a)
       *  = [eta + 3/2 (a'/a)^2 / k^2 /s_2^2 {delta_tot + 3 (a'/a) /k^2 velocity_tot}] / (a'/a)
       *
       *  with
       *
       *  delta_tot = (delta_rho/rho_c)
       *  = [rho_r delta_r + rho_m delta_m] / (rho_r + rho_m)
       *  = [delta_r + (rho_m/rho_r) delta_m] / (1 + rho_m/rho_r)
       *  = [(f_g delta_g + f_nu delta_nu) + (rho_m/rho_r) (f_b delta_b + f_cdm delta_cdm)] / (1 + rho_m/rho_r)
       *
       *  velocity_tot = ((rho+p)theta/rho_c)
       *  = [(4/3) rho_r theta_r + rho_m theta_m] / (rho_r + rho_m)
       *  = [(4/3) theta_r + (rho_m/rho_r) theta_m] / (1 + rho_m/rho_r)
       *  = [(4/3) (f_g theta_g + f_nu theta_nu) + (rho_m/rho_r) (f_b delta_b + f_cdm 0)] / (1 + rho_m/rho_r)
       */

      if (pba->has_cdm == _TRUE_)
        delta_cdm = ppw->pv->y[ppw->pv->index_pt_delta_cdm];
      else
        delta_cdm = 0.;

      /* note: if there are no neutrinos, fracnu, delta_ur and theta_ur below will consistently be zero. */

      delta_tot = (fracg * ppw->pv->y[ppw->pv->index_pt_delta_g] + fracnu * delta_ur + rho_m_over_rho_r * (fracb * ppw->pv->y[ppw->pv->index_pt_delta_b] + fraccdm * delta_cdm)) / (1. + rho_m_over_rho_r);

      velocity_tot = ((4. / 3.) * (fracg * ppw->pv->y[ppw->pv->index_pt_theta_g] + fracnu * theta_ur) + rho_m_over_rho_r * fracb * ppw->pv->y[ppw->pv->index_pt_theta_b]) / (1. + rho_m_over_rho_r);

      alpha = (eta + 3. / 2. * a_prime_over_a * a_prime_over_a / k / k / s2_squared * (delta_tot + 3. * a_prime_over_a / k / k * velocity_tot)) / a_prime_over_a;

      ppw->pv->y[ppw->pv->index_pt_phi] = eta - a_prime_over_a * alpha;

      ppw->pv->y[ppw->pv->index_pt_delta_g] -= 4. * a_prime_over_a * alpha;
      ppw->pv->y[ppw->pv->index_pt_theta_g] += k * k * alpha;

      ppw->pv->y[ppw->pv->index_pt_delta_b] -= 3. * a_prime_over_a * alpha;
      ppw->pv->y[ppw->pv->index_pt_theta_b] += k * k * alpha;

      if (pba->has_cdm == _TRUE_)
      {
        ppw->pv->y[ppw->pv->index_pt_delta_cdm] -= 3. * a_prime_over_a * alpha;
        ppw->pv->y[ppw->pv->index_pt_theta_cdm]  = k * k * alpha;
      }

      /* fluid */
      if ((pba->has_fld == _TRUE_) && (pba->use_ppf == _FALSE_))
      {
        class_call (background_w_fld (pba, a, &w_fld, &dw_over_da_fld, &integral_fld), pba->error_message, ppt->error_message);

        ppw->pv->y[ppw->pv->index_pt_delta_fld] += 3 * (1. + w_fld) * a_prime_over_a * alpha;
        ppw->pv->y[ppw->pv->index_pt_theta_fld] += k * k * alpha;
      }

      if ((pba->has_ur == _TRUE_) || (pba->has_ncdm == _TRUE_))
      {
        delta_ur -= 4. * a_prime_over_a * alpha;
        theta_ur += k * k * alpha;
        /* shear and l3 are gauge invariant */
      }
    } /* end of gauge transformation to newtonian gauge */

    /** - (e) In any gauge, we should now implement the relativistic initial conditions in ur and ncdm variables */

    if (pba->has_ur == _TRUE_)
    {
      ppw->pv->y[ppw->pv->index_pt_delta_ur] = delta_ur;

      ppw->pv->y[ppw->pv->index_pt_theta_ur] = theta_ur;

      ppw->pv->y[ppw->pv->index_pt_shear_ur] = shear_ur;

      ppw->pv->y[ppw->pv->index_pt_l3_ur] = l3_ur;
    }

    if (pba->has_ncdm == _TRUE_)
    {
      idx = ppw->pv->index_pt_psi0_ncdm1;

      for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
      {
        for (index_q = 0; index_q < ppw->pv->q_size_ncdm[n_ncdm]; index_q++)
        {
          q = pba->q_ncdm[n_ncdm][index_q];

          epsilon = sqrt (q * q + a * a * pba->M_ncdm[n_ncdm] * pba->M_ncdm[n_ncdm]);

          ppw->pv->y[idx] = -0.25 * delta_ur * pba->dlnf0_dlnq_ncdm[n_ncdm][index_q];

          ppw->pv->y[idx + 1] =  -epsilon / 3. / q / k * theta_ur * pba->dlnf0_dlnq_ncdm[n_ncdm][index_q];

          ppw->pv->y[idx + 2] = -0.5 * shear_ur * pba->dlnf0_dlnq_ncdm[n_ncdm][index_q];

          ppw->pv->y[idx + 3] = -0.25 * l3_ur * pba->dlnf0_dlnq_ncdm[n_ncdm][index_q];

          /*Jump to next momentum bin: */
          idx += (ppw->pv->l_max_ncdm[n_ncdm] + 1);
        }
      }
    }
  }

  /** --> For tensors */

  if (_tensors_)
  {
    /** tensor initial conditions take into account the fact that
     *  scalar (resp. tensor) \f$ C_l\f$'s are related to the real space
     *  power spectrum of curvature (resp. of the tensor part of
     *  metric perturbations)
     *
     *  \f[ <R(x) R(x)>  \ \  \sum_{ij} <h_{ij}(x) h^{ij}(x)> \f]
     *
     *  In momentum space it is conventional to use the modes R(k)
     *  and h(k) where the quantity h obeying to the equation of
     *  propagation:
     *
     *  \f[ h'' + \frac{2a'}{a} h + [k2+2K] h = 12\pi Ga2 (\rho+p) \sigma = 8\pi Ga2 p \pi \f]
     *
     *  and the power spectra in real space and momentum space are related through:
     *
     *  \f[ <R(x) R(x)> = \int \frac{dk}{k} \left[ \frac{k^3}{2\pi^2} <R(k)R(k)^*>\right] = \int \frac{dk}{k} \mathcal{P}_R(k) \f]
     *  \f[\sum_{ij} <h_{ij}(x) h^{ij}(x)> = \frac{dk}{k} \left[ \frac{k^3}{2\pi^2} F\left(\frac{k^2}{K}\right) <h(k)h(k)^*>\right] = \int \frac{dk}{k} F\left(\frac{k^2}{K}\right) \mathcal{P}_h(k) \f]
     *
     *  where \f$ \mathcal{P}_R\f$ and \f$ \mathcal{P}_h\f$ are the dimensionless spectrum of
     *  curvature R, and F is a function of k2/K, where K is the curvature
     *  parameter. F is equal to one in flat space (K=0), and coming
     *  from the contraction of the laplacian eigentensor \f$ Q_{ij}\f$ with
     *  itself. We will give F explicitly below.
     *
     *  Similarly the scalar (S) and tensor (T) \f$ C_l\f$'s are given by
     *
     *  \f[ C_l^S = 4\pi \int \frac{dk}{k} [\Delta_l^S(q)]^2 \mathcal{P}_R(k) \f]
     *  \f[ C_l^T = 4\pi \int \frac{dk}{k} [\Delta_l^T(q)]^2 F\left(\frac{k^2}{K}\right) \mathcal{P}_h(k) \f]
     *
     *  The usual convention for the tensor-to-scalar ratio
     *  \f$ r = A_t / A_s \f$ at pivot scale
     *  = 16 epsilon in single-field inflation
     *  is such that for constant \f$ \mathcal{P}_R(k)\f$ and \f$ \mathcal{P}_h(k)\f$,
     *
     *  \f[ r = 6 \frac{\mathcal{P}_h(k)}{\mathcal{P}_R(k)} \f]
     *
     *  so
     *
     *  \f[ \mathcal{P}_h(k) = \frac{\mathcal{P}_R(k) r}{6} = \frac{A_s r}{6} = \frac{A_t}{6} \f]
     *
     *  A priori it would make sense to say that for a power-law
     *  primordial spectrum there is an extra factor \f$ (k/k_{pivot})^{n_t} \f$
     *  (and eventually running and so on and so forth...)
     *
     *  However it has been shown that the minimal models of
     *  inflation in a negatively curved bubble lead to
     *  \f$ \mathcal{P}_h(k)=\tanh(\pi*\nu/2)\f$. In open models it is customary to
     *  define the tensor tilt in a non-flat universe as a deviation
     *  from this behavior rather than from true scale-invariance in
     *  the above sense.
     *
     *  Hence we should have
     *
     *  \f[ \mathcal{P}_h(k) = \frac{A_t}{6} [ \tanh(\pi*\frac{\nu}{2})]  (k/k_{pivot})^{(n_t+...)}\f]
     *
     *  where the brackets \f[ [...] \f] mean "if K<0"
     *
     *  Then
     *
     *  \f[ C_l^T = 4\pi \int \frac{dk}{k} [\Delta_l^T(q)]^2 F\left(\frac{k^2}{K}\right) \frac{A_t}{6} [\tanh(\pi*\frac{\nu}{2})] (k/k_{pivot})^{(n_t+...)} \f]
     *
     *  In the code, it is then a matter of choice to write:
     *
     *  - In the primordial module: \f$ \mathcal{P}_h(k) = \frac{A_t}{6} \tanh{(\pi*\frac{\nu}{2})} (k/k^*)^{n_T}\f$
     *  - In the perturbation initial conditions: \f$ h = 1\f$
     *  - In the spectra module: \f$ C_l^T = \frac{4}{\pi} \int \frac{dk}{k} [\Delta_l^T(q)]^2 F\left(\frac{k^2}{K}\right) \mathcal{P}_h(k) \f$
     *
     *  or:
     *
     *  - In the primordial module: \f$ \mathcal{P}_h(k) = A_t (k/k^*)^{n_T} \f$
     *  - In the perturbation initial conditions: \f$ h = \sqrt{[F\left(\frac{k^2}{K}\right) / 6] \tanh{(\pi*\frac{\nu}{2})}} \f$
     *  - In the spectra module: \f$ C_l^T = \frac{4}{\pi} \int \frac{dk}{k} [\Delta_l^T(q)]^2 \mathcal{P}_h(k) \f$
     *
     *  We choose this last option, such that the primordial and
     *  spectra module differ minimally in flat and non-flat space. Then we must impose
     *
     *  \f[ h = \sqrt{\left(\frac{F}{6}\right) \tanh{(\pi*\frac{\nu}{2})}} \f]
     *
     *  The factor F is found to be given by:
     *
     *  \f[ \sum_{ij}<h_{ij}(x) h^{ij}(x)> = \int \frac{dk}{k}  \frac{k2(k2-K)}{(k2+3K)(k2+2K)} \mathcal{P}_h(k) \f]
     *
     *  Introducing as usual \f$ q2 = k2 - 3K \f$  and using qdq = kdk this gives
     *
     *  \f[ \sum_{ij}<h_{ij}(x) h^{ij}(x)> = \int \frac{dk}{k} \frac{(q2-3K)(q2-4K)}{q2(q2-K)} \mathcal{P}_h(k) \f]
     *
     *  Using qdq = kdk this is equivalent to
     *
     *  \f[ \sum_{ij}<h_{ij}(x) h^{ij}(x)> = \int \frac{dq}{q} \frac{q2-4K}{q2-K} \mathcal{P}_h(k(q)) \f]
     *
     *  Finally, introducing \f$ \nu=q/\sqrt{|K|}\f$ and sgnK=SIGN(k)\f$=\pm 1\f$, this could also be written
     *
     *  \f[ \sum_{ij}<h_{ij}(x) h^{ij}(x)> = \int \frac{d\nu}{\nu} \frac{(\nu2-4sgnK)}{(\nu2-sgnK)} \mathcal{P}_h(k(\nu)) \f]
     *
     *  Equation (43,44) of Hu, Seljak, White, Zaldarriaga is
     *  equivalent to absorbing the above factor
     *  \f$ (\nu2-4sgnK)/(\nu2-sgnK)\f$ in the definition of the primordial
     *  spectrum. Since the initial condition should be written in terms of k rather than nu, they should read
     *
     *  \f[ h = \sqrt{ [k2(k2-K)]/[(k2+3K)(k2+2K)] / 6 * \tanh{(\pi*\frac{\nu}{2})} } \f]
     *
     *  We leave the freedom to multiply by an arbitrary number
     *  ppr->gw_ini. The standard convention corresponding to
     *  standard definitions of r, \f$ A_T\f$, \f$ n_T\f$ is however ppr->gw_ini=1.
     *
     */

    if (index_ic == ppt->index_ic_ten)
      ppw->pv->y[ppw->pv->index_pt_gw] = ppr->gw_ini / _SQRT6_;

    k2 = k * k;

    if (sgnK != 0)
      ppw->pv->y[ppw->pv->index_pt_gw] *= sqrt (k2 * (k2 - K) / (k2 + 3. * K) / (k2 + 2. * K));

    if (sgnK == -1)
    {
      if (k * k + 3 * K >= 0.)
        ppw->pv->y[ppw->pv->index_pt_gw] *= sqrt (tanh (_PI_ / 2. * sqrt (k2 + 3 * K) / sqrt (-K)));
      else
        ppw->pv->y[ppw->pv->index_pt_gw] = 0.;
    }
  }

  return _SUCCESS_;
}

/**
 * Evaluate background/thermodynamics at \f$ \tau \f$, infer useful flags / time scales for integrating perturbations.
 *
 * Evaluate background quantities at \f$ \tau \f$, as well as thermodynamics for scalar mode; infer useful flags and time scales for integrating the perturbations:
 * - check whether tight-coupling approximation is needed.
 * - check whether radiation (photons, massless neutrinos...) perturbations are needed.
 * - choose step of integration: step = ppr->perturb_integration_stepsize * min_time_scale, where min_time_scale = smallest time scale involved in the equations. There are three time scales to compare:
 *     -# that of recombination, \f$ \tau_c = 1/\kappa' \f$
 *     -# Hubble time scale, \f$ \tau_h = a/a' \f$
 *     -# Fourier mode, \f$ \tau_k = 1/k \f$
 *
 * So, in general, min_time_scale = \f$ \min(\tau_c, \tau_b, \tau_h, \tau_k) \f$.
 *
 * However, if \f$ \tau_c \ll \tau_h \f$ and \f$ \tau_c
 * \ll \tau_k \f$, we can use the tight-coupling regime for photons
 * and write equations in such way that the time scale \f$
 * \tau_c \f$ becomes irrelevant (no effective mass term in \f$
 * 1/\tau_c \f$).  Then, the smallest
 * scale in the equations is only \f$ \min(\tau_h, \tau_k) \f$.
 * In practise, it is sufficient to use only the condition \f$ \tau_c \ll \tau_h \f$.
 *
 * Also, if \f$ \rho_{matter} \gg \rho_{radiation} \f$ and \f$ k \gg
 * aH \f$, we can switch off radiation perturbations (i.e. switch on
 * the free-streaming approximation) and then the smallest scale is
 * simply \f$ \tau_h \f$.
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to thermodynamics structure
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_md   Input: index of mode under consideration (scalar/.../tensor)
 * @param k          Input: wavenumber
 * @param tau        Input: conformal time
 * @param ppw        Input/Output: in output contains the approximation to be used at this time
 * @return the error status
 */

int
perturb_approximations (
  struct precision         *ppr,
  const struct background  *pba,
  struct thermo            *pth,
  struct perturbs          *ppt,
  int                      index_md,
  double                   k,
  double                   tau,
  struct perturb_workspace *ppw
                       )
{
  /** Summary: */

  /** - define local variables */

  /* (a) time scale of Fourier mode, \f$ \tau_k = 1/k \f$ */
  double tau_k;
  /* (b) time scale of expansion, \f$ \tau_h = a/a' \f$ */
  double tau_h;
  /* (c) time scale of recombination, \f$ \tau_{\gamma} = 1/\kappa' \f$ */
  double tau_c;

  /** - compute Fourier mode time scale = \f$ \tau_k = 1/k \f$ */

  class_test (k == 0.,
              ppt->error_message,
              "stop to avoid division by zero");

  tau_k = 1. / k;

  /** - evaluate background quantities with background_at_tau() and
   *   Hubble time scale \f$ \tau_h = a/a' \f$ */

  class_call (background_at_tau (pba, tau, pba->normal_info, ppw->inter_mode, &(ppw->last_index_back), ppw->pvecback),
              pba->error_message,
              ppt->error_message);

  class_test (ppw->pvecback[pba->index_bg_H] * ppw->pvecback[pba->index_bg_a] == 0.,
              ppt->error_message,
              "aH=0, stop to avoid division by zero");

  tau_h = 1. / (ppw->pvecback[pba->index_bg_H] * ppw->pvecback[pba->index_bg_a]);

  /** - for scalar modes: */

  if (_scalars_)
  {
    /** - --> (a) evaluate thermodynamical quantities with thermodynamics_at_z() */

    class_call (thermodynamics_at_z (pba,
                                     pth,
                                     1. / ppw->pvecback[pba->index_bg_a] - 1., /* redshift z=1/a-1 */
                                     ppw->inter_mode,
                                     &(ppw->last_index_thermo),
                                     ppw->pvecback,
                                     ppw->pvecthermo),
                pth->error_message,
                ppt->error_message);

    /** - ---> (b.1.) if \f$ \kappa'=0 \f$, recombination is finished; tight-coupling approximation must be off */

    if (ppw->pvecthermo[pth->index_th_dkappa] == 0.)
    {
      ppw->approx[ppw->index_ap_tca] = (int) tca_off;
    }
    /** - ---> (b.2.) if \f$ \kappa' \neq 0 \f$, recombination is not finished: check tight-coupling approximation */

    else
    {
      /** - ----> (b.2.a) compute recombination time scale for photons, \f$ \tau_{\gamma} = 1/ \kappa' \f$ */
      tau_c = 1. / ppw->pvecthermo[pth->index_th_dkappa];

      class_test (tau_c < 0.,
                  ppt->error_message,
                  "tau_c = 1/kappa' should always be positive unless there is something wrong in the thermodynamics module. However you have here tau_c=%e at z=%e, conformal time=%e x_e=%e. (This could come from the interpolation of a too poorly sampled reionisation history?).\n",
                  tau_c,
                  1. / ppw->pvecback[pba->index_bg_a] - 1.,
                  tau,
                  ppw->pvecthermo[pth->index_th_xe]);

      /** - ----> (b.2.b) check whether tight-coupling approximation should be on */

      if ((tau_c / tau_h < ppr->tight_coupling_trigger_tau_c_over_tau_h) &&
          (tau_c / tau_k < ppr->tight_coupling_trigger_tau_c_over_tau_k))
        ppw->approx[ppw->index_ap_tca] = (int) tca_on;
      else
        ppw->approx[ppw->index_ap_tca] = (int) tca_off;
    }

    /** - --> (c) free-streaming approximations */

    if ((tau / tau_k > ppr->radiation_streaming_trigger_tau_over_tau_k) &&
        (tau > pth->tau_free_streaming) &&
        (ppr->radiation_streaming_approximation != rsa_none))

      ppw->approx[ppw->index_ap_rsa] = (int) rsa_on;
    else
      ppw->approx[ppw->index_ap_rsa] = (int) rsa_off;

    if (pba->has_ur == _TRUE_)
    {
      if ((tau / tau_k > ppr->ur_fluid_trigger_tau_over_tau_k) &&
          (ppr->ur_fluid_approximation != ufa_none))

        ppw->approx[ppw->index_ap_ufa] = (int) ufa_on;
      else
        ppw->approx[ppw->index_ap_ufa] = (int) ufa_off;
    }

    if (pba->has_ncdm == _TRUE_)
    {
      if ((tau / tau_k > ppr->ncdm_fluid_trigger_tau_over_tau_k) &&
          (ppr->ncdm_fluid_approximation != ncdmfa_none))

        ppw->approx[ppw->index_ap_ncdmfa] = (int) ncdmfa_on;
      else
        ppw->approx[ppw->index_ap_ncdmfa] = (int) ncdmfa_off;
    }
  }

  /** - for tensor modes: */

  if (_tensors_)
  {
    /** - --> (a) evaluate thermodynamical quantities with thermodynamics_at_z() */

    class_call (thermodynamics_at_z (pba,
                                     pth,
                                     1. / ppw->pvecback[pba->index_bg_a] - 1., /* redshift z=1/a-1 */
                                     ppw->inter_mode,
                                     &(ppw->last_index_thermo),
                                     ppw->pvecback,
                                     ppw->pvecthermo),
                pth->error_message,
                ppt->error_message);

    /** - ---> (b.1.) if \f$ \kappa'=0 \f$, recombination is finished; tight-coupling approximation must be off */

    if (ppw->pvecthermo[pth->index_th_dkappa] == 0.)
    {
      ppw->approx[ppw->index_ap_tca] = (int) tca_off;
    }
    /** - ---> (b.2.) if \f$ \kappa' \neq 0 \f$, recombination is not finished: check tight-coupling approximation */

    else
    {
      /** - ----> (b.2.a) compute recombination time scale for photons, \f$ \tau_{\gamma} = 1/ \kappa' \f$ */
      tau_c = 1. / ppw->pvecthermo[pth->index_th_dkappa];

      /** - ----> (b.2.b) check whether tight-coupling approximation should be on */
      if ((tau_c / tau_h < ppr->tight_coupling_trigger_tau_c_over_tau_h) &&
          (tau_c / tau_k < ppr->tight_coupling_trigger_tau_c_over_tau_k))
        ppw->approx[ppw->index_ap_tca] = (int) tca_on;
      else
        ppw->approx[ppw->index_ap_tca] = (int) tca_off;
    }

    if ((tau / tau_k > ppr->radiation_streaming_trigger_tau_over_tau_k) &&
        (tau > pth->tau_free_streaming) &&
        (ppr->radiation_streaming_approximation != rsa_none))

      ppw->approx[ppw->index_ap_rsa] = (int) rsa_on;
    else
      ppw->approx[ppw->index_ap_rsa] = (int) rsa_off;
  }

  return _SUCCESS_;
}

/**
 * Compute typical timescale over which the perturbation equations
 * vary. Some integrators (e.g. Runge-Kunta) benefit from calling this
 * routine at each step in order to adapt the next step.
 *
 * This is one of the few functions in the code which is passed to the generic_integrator() routine.
 * Since generic_integrator() should work with functions passed from various modules, the format of the arguments
 * is a bit special:
 * - fixed parameters and workspaces are passed through a generic pointer.
 *   generic_integrator() doesn't know the content of this pointer.
 * - the error management is a bit special: errors are not written as usual to pth->error_message, but to a generic
 *   error_message passed in the list of arguments.
 *
 * @param tau                      Input: conformal time
 * @param parameters_and_workspace Input: fixed parameters (e.g. indices), workspace, approximation used, etc.
 * @param timescale                Output: perturbation variation timescale (given the approximation used)
 * @param error_message            Output: error message
 */

int
perturb_timescale (
  double   tau,
  void     *parameters_and_workspace,
  double   *timescale,
  ErrorMsg error_message
                  )
{
  /** Summary: */

  /** - define local variables */

  /* (a) time scale of Fourier mode, \f$ \tau_k = 1/k \f$ */
  double tau_k;
  /* (b) time scale of expansion, \f$ \tau_h = a/a' \f$ */
  double tau_h;
  /* (c) time scale of recombination, \f$ \tau_{\gamma} = 1/\kappa' \f$ */
  double tau_c;

  /* various pointers allowing to extract the fields of the
   *  parameter_and_workspace input structure */
  struct perturb_parameters_and_workspace *pppaw;
  const struct background *pba;
  struct thermo *pth;
  struct perturbs *ppt;
  struct perturb_workspace *ppw;
  double *pvecback;
  double *pvecthermo;

  /** - extract the fields of the parameter_and_workspace input structure */
  pppaw      = parameters_and_workspace;
  pba        = pppaw->pba;
  pth        = pppaw->pth;
  ppt        = pppaw->ppt;
  ppw        = pppaw->ppw;
  pvecback   = ppw->pvecback;
  pvecthermo = ppw->pvecthermo;

  /** - compute Fourier mode time scale = \f$ \tau_k = 1/k \f$ */

  class_test (pppaw->k == 0.,
              ppt->error_message,
              "stop to avoid division by zero");

  tau_k = 1. / pppaw->k;

  /** - evaluate background quantities with background_at_tau() and
   *   Hubble time scale \f$ \tau_h = a/a' \f$ */

  class_call (background_at_tau (pba, tau, pba->normal_info, ppw->inter_mode, &(ppw->last_index_back), pvecback),
              pba->error_message,
              error_message);

  class_test (pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a] == 0.,
              error_message,
              "aH=0, stop to avoid division by zero");

  tau_h = 1. / (pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a]);

  /** - for scalars modes: */

  if ((ppt->has_scalars == _TRUE_) && (pppaw->index_md == ppt->index_md_scalars))
  {
    *timescale = tau_h;

    if ((ppw->approx[ppw->index_ap_rsa] == (int) rsa_off) || (pba->has_ncdm == _TRUE_))
      *timescale = MIN (tau_k, *timescale);

    if (ppw->approx[ppw->index_ap_tca] == (int) tca_off)
    {
      class_call (thermodynamics_at_z (pba,
                                       pth,
                                       1. / pvecback[pba->index_bg_a] - 1., /* redshift z=1/a-1 */
                                       ppw->inter_mode,
                                       &(ppw->last_index_thermo),
                                       pvecback,
                                       pvecthermo),
                  pth->error_message,
                  error_message);

      if (pvecthermo[pth->index_th_dkappa] != 0.)
      {
        /** - -->  compute recombination time scale for photons, \f$ \tau_{\gamma} = 1/ \kappa' \f$ */

        tau_c = 1. / pvecthermo[pth->index_th_dkappa];

        *timescale = MIN (tau_c, *timescale);
      }
    }
  }

  /** - for vector modes: */

  if ((ppt->has_vectors == _TRUE_) && (pppaw->index_md == ppt->index_md_vectors))
  {
    *timescale = MIN (tau_h, tau_k);

    if (ppw->approx[ppw->index_ap_tca] == (int) tca_off)
    {
      class_call (thermodynamics_at_z (pba,
                                       pth,
                                       1. / pvecback[pba->index_bg_a] - 1., /* redshift z=1/a-1 */
                                       ppw->inter_mode,
                                       &(ppw->last_index_thermo),
                                       pvecback,
                                       pvecthermo),
                  pth->error_message,
                  error_message);

      if (pvecthermo[pth->index_th_dkappa] != 0.)
      {
        /** - -->  compute recombination time scale for photons, \f$ \tau_{\gamma} = 1/ \kappa' \f$ */

        tau_c = 1. / pvecthermo[pth->index_th_dkappa];

        *timescale = MIN (tau_c, *timescale);
      }
    }
  }

  /** - for tensor modes: */

  if ((ppt->has_tensors == _TRUE_) && (pppaw->index_md == ppt->index_md_tensors))
  {
    *timescale = MIN (tau_h, tau_k);

    if (ppw->approx[ppw->index_ap_tca] == (int) tca_off)
    {
      class_call (thermodynamics_at_z (pba,
                                       pth,
                                       1. / pvecback[pba->index_bg_a] - 1., /* redshift z=1/a-1 */
                                       ppw->inter_mode,
                                       &(ppw->last_index_thermo),
                                       pvecback,
                                       pvecthermo),
                  pth->error_message,
                  error_message);

      if (pvecthermo[pth->index_th_dkappa] != 0.)
      {
        /** - --> compute recombination time scale for photons, \f$ \tau_{\gamma} = 1/ \kappa' \f$ */

        tau_c = 1. / pvecthermo[pth->index_th_dkappa];

        *timescale = MIN (tau_c, *timescale);
      }
    }
  }

  return _SUCCESS_;
}

/**
 * Compute metric perturbations (those not integrated over time) using Einstein equations
 *
 * @param ppr        Input: pointer to precision structure
 * @param pba        Input: pointer to background structure
 * @param pth        Input: pointer to thermodynamics structure
 * @param ppt        Input: pointer to the perturbation structure
 * @param index_md   Input: index of mode under consideration (scalar/.../tensor)
 * @param k          Input: wavenumber
 * @param tau        Input: conformal time
 * @param y          Input: vector of perturbations (those integrated over time) (already allocated)
 * @param ppw        Input/Output: in output contains the updated metric perturbations
 * @return the error status
 */

int
perturb_einstein (
  struct precision         *ppr,
  const struct background  *pba,
  struct thermo            *pth,
  struct perturbs          *ppt,
  int                      index_md,
  double                   k,
  double                   tau,
  double                   *y,
  struct perturb_workspace *ppw
                 )
{
  /** Summary: */

  /** - define local variables */

  int sgnK;
  double K = nc_hicosmo_curvature_K (pba->cosmo, &sgnK);

  double k2, a, a2, a_prime_over_a;
  double s2_squared;
  double shear_g = 0.;

  /** - define wavenumber and scale factor related quantities */

  k2             = k * k;
  a              = ppw->pvecback[pba->index_bg_a];
  a2             = a * a;
  a_prime_over_a = ppw->pvecback[pba->index_bg_H] * a;
  s2_squared     = 1. - 3. * K / k2;

  /** - sum up perturbations from all species */
  class_call (perturb_total_stress_energy (ppr, pba, pth, ppt, index_md, k, y, ppw),
              ppt->error_message,
              ppt->error_message);

  /** - for scalar modes: */

  if (_scalars_)
  {
    /** - --> infer metric perturbations from Einstein equations */

    /* newtonian gauge */
    if (ppt->gauge == newtonian)
    {
      /* in principle we could get phi from the constrain equation:
       *
       *  ppw->pvecmetric[ppw->index_mt_phi] = -1.5 * (a2/k2/k2/s2/s2) * (k2 * delta_rho + 3.*a_prime_over_a * rho_plus_p_theta);
       *
       *  with s2_squared = sqrt(1-3K/k2) = ppw->s_l[2]*ppw->s_l[2]
       *
       *  This was the case in class v1.3. However the integration is
       *  more stable is we treat phi as a dynamical variable
       *  y[ppw->pv->index_pt_phi], which derivative is given by the
       *  second equation below (credits to Guido Walter Pettinari). */

      /* equation for psi */
      ppw->pvecmetric[ppw->index_mt_psi] = y[ppw->pv->index_pt_phi] - 4.5 * (a2 / k2) * ppw->rho_plus_p_shear;

      /* equation for phi' */
      ppw->pvecmetric[ppw->index_mt_phi_prime] = -a_prime_over_a * ppw->pvecmetric[ppw->index_mt_psi] + 1.5 * (a2 / k2) * ppw->rho_plus_p_theta;

      /* eventually, infer radiation streaming approximation for
       *  gamma and ur (this is exactly the right place to do it
       *  because the result depends on h_prime) */

      if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
        class_call (perturb_rsa_delta_and_theta (ppr, pba, pth, ppt, k, y, a_prime_over_a, ppw->pvecthermo, ppw),
                    ppt->error_message,
                    ppt->error_message);
    }

    /* synchronous gauge */
    if (ppt->gauge == synchronous)
    {
      /* first equation involving total density fluctuation */
      ppw->pvecmetric[ppw->index_mt_h_prime] =
        (k2 * s2_squared * y[ppw->pv->index_pt_eta] + 1.5 * a2 * ppw->delta_rho) / (0.5 * a_prime_over_a); /* h' */

      /* eventually, infer radiation streaming approximation for
       *  gamma and ur (this is exactly the right place to do it
       *  because the result depends on h_prime) */

      if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
      {
        class_call (perturb_rsa_delta_and_theta (ppr, pba, pth, ppt, k, y, a_prime_over_a, ppw->pvecthermo, ppw),
                    ppt->error_message,
                    ppt->error_message);

        /* update total theta given rsa approximation results */

        ppw->rho_plus_p_theta += 4. / 3. * ppw->pvecback[pba->index_bg_rho_g] * ppw->rsa_theta_g;

        if (pba->has_ur == _TRUE_)
          ppw->rho_plus_p_theta += 4. / 3. * ppw->pvecback[pba->index_bg_rho_ur] * ppw->rsa_theta_ur;
      }

      /* second equation involving total velocity */
      ppw->pvecmetric[ppw->index_mt_eta_prime] = (1.5 * a2 * ppw->rho_plus_p_theta + 0.5 * K * ppw->pvecmetric[ppw->index_mt_h_prime]) / k2 / s2_squared; /* eta' */

      /* third equation involving total pressure */
      ppw->pvecmetric[ppw->index_mt_h_prime_prime] =
        -2. * a_prime_over_a * ppw->pvecmetric[ppw->index_mt_h_prime]
        + 2. * k2 * s2_squared * y[ppw->pv->index_pt_eta]
        - 9. * a2 * ppw->delta_p;

      /* alpha = (h'+6eta')/2k^2 */
      ppw->pvecmetric[ppw->index_mt_alpha] = (ppw->pvecmetric[ppw->index_mt_h_prime] + 6. * ppw->pvecmetric[ppw->index_mt_eta_prime]) / 2. / k2;

      /* eventually, infer first-order tight-coupling approximation for photon
       *  shear, then correct the total shear */
      if (ppw->approx[ppw->index_ap_tca] == (int) tca_on)
      {
        shear_g = 16. / 45. / ppw->pvecthermo[pth->index_th_dkappa] * (y[ppw->pv->index_pt_theta_g] + k2 * ppw->pvecmetric[ppw->index_mt_alpha]);

        ppw->rho_plus_p_shear += 4. / 3. * ppw->pvecback[pba->index_bg_rho_g] * shear_g;
      }

      /* fourth equation involving total shear */
      ppw->pvecmetric[ppw->index_mt_alpha_prime] = /*TBC */
                                                   -2. * a_prime_over_a * ppw->pvecmetric[ppw->index_mt_alpha]
                                                   + y[ppw->pv->index_pt_eta]
                                                   - 4.5 * (a2 / k2) * ppw->rho_plus_p_shear;
    }

    /* transform (delta_m, theta_m) of the current gauge into
     *  gauge-independent variables (you could comment this out if you
     *  really want gauge-dependent results) */

    if (ppt->has_source_delta_m == _TRUE_)
      ppw->delta_m += 3. * ppw->pvecback[pba->index_bg_a] * ppw->pvecback[pba->index_bg_H] * ppw->theta_m / k2;

    /* note: until 2.4.3 there was a typo, the factor was (-2 H'/H) instead */
    /* of (3 aH). There is the same typo in the CLASSgal paper */
    /* 1307.1459v1,v2,v3. It came from a confusion between (1+w_total) */
    /* and (1+w_matter)=1 [the latter is the relevant one here]. */
    /* */
    /* note2: at this point this gauge-invariant variable is only */
    /* valid if all matter components are pressureless and */
    /* stable. This relation will be generalized soon to the case */
    /* of decaying dark matter. */

    if (ppt->has_source_delta_cb == _TRUE_)
      ppw->delta_cb += 3. * ppw->pvecback[pba->index_bg_a] * ppw->pvecback[pba->index_bg_H] * ppw->theta_cb / k2;  /*check gauge transformation */

    if (ppt->has_source_theta_m == _TRUE_)
      if  (ppt->gauge == synchronous)
        ppw->theta_m += ppw->pvecmetric[ppw->index_mt_alpha] * k2;


    if (ppt->has_source_theta_cb == _TRUE_)
      if  (ppt->gauge == synchronous)
        ppw->theta_cb += ppw->pvecmetric[ppw->index_mt_alpha] * k2;


    /*check gauge transformation */
  }

  /** - for vector modes */

  if (_vectors_)
  {
    if (ppt->gauge == newtonian)
      ppw->pvecmetric[ppw->index_mt_V_prime] = -2. * a_prime_over_a * y[ppw->pv->index_pt_V] - 3. * ppw->vector_source_pi / k;


    if (ppt->gauge == synchronous)
      /* assuming    vector_source_pi = p_class a^2 pi_T^{(1)} and  vector_source_v = (rho_class+p_class)a^2 v^{(1)} */

      /* from Hu and White: */
      ppw->pvecmetric[ppw->index_mt_hv_prime_prime] = -2. * a_prime_over_a * y[ppw->pv->index_pt_hv_prime] - 3. * ppw->vector_source_pi / k2;

    /* what we suspect: */
    /*ppw->pvecmetric[ppw->index_mt_hv_prime_prime] = -2.*a_prime_over_a*y[ppw->pv->index_pt_hv_prime] - 3.*ppw->vector_source_pi; */

    /* if we use the other equation: */
    /*ppw->pvecmetric[ppw->index_mt_hv_prime] = -2./k/ (1.-2.*K/k2) * 3. * ppw->vector_source_v; */
  }

  /** - for tensor modes */

  if (_tensors_)
    /* single einstein equation for tensor perturbations */
    ppw->pvecmetric[ppw->index_mt_gw_prime_prime] = -2. * a_prime_over_a * y[ppw->pv->index_pt_gwdot] - (k2 + 2. * K) * y[ppw->pv->index_pt_gw] + ppw->gw_source;



  return _SUCCESS_;
}

int
perturb_total_stress_energy (
  struct precision         *ppr,
  const struct background  *pba,
  struct thermo            *pth,
  struct perturbs          *ppt,
  int                      index_md,
  double                   k,
  double                   *y,
  struct perturb_workspace *ppw
                            )
{
  /** Summary: */

  /** - define local variables */
  /* const double H0 = 1.0 / nc_hicosmo_RH_Mpc (pba->cosmo); */
  int sgnK;
  double K = nc_hicosmo_curvature_K (pba->cosmo, &sgnK);

  double a, a2, a_prime_over_a, k2;
  double rho_plus_p_tot        = 0.;
  double delta_g               = 0.;
  double theta_g               = 0.;
  double shear_g               = 0.;
  double delta_ur              = 0.;
  double theta_ur              = 0.;
  double shear_ur              = 0.;
  double rho_delta_ncdm        = 0.;
  double rho_plus_p_theta_ncdm = 0.;
  double rho_plus_p_shear_ncdm = 0.;
  double delta_p_ncdm          = 0.;
  double factor;
  double rho_plus_p_ncdm;
  int index_q, n_ncdm, idx;
  double epsilon, q, q2, cg2_ncdm, w_ncdm, rho_ncdm_bg, p_ncdm_bg, pseudo_p_ncdm;
  double rho_m, delta_rho_m, rho_plus_p_m, rho_plus_p_theta_m;
  double w_fld, dw_over_da_fld, integral_fld;
  double gwncdm;
  double rho_relativistic;
  double c_gamma_k_H_square;
  double Gamma_prime_plus_a_prime_over_a_Gamma, alpha = 0., s2sq = 1.;

  /** - wavenumber and scale factor related quantities */

  a              = ppw->pvecback[pba->index_bg_a];
  a2             = a * a;
  a_prime_over_a = ppw->pvecback[pba->index_bg_H] * a;
  k2             = k * k;

  /** - for scalar modes */

  if (_scalars_)
  {
    /** - --> (a) deal with approximation schemes */

    /** - ---> (a.1.) photons */

    if (ppw->approx[ppw->index_ap_tca] == (int) tca_off)
    {
      if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
      {
        /** - ----> (a.1.1.) no approximation */

        delta_g = y[ppw->pv->index_pt_delta_g];
        theta_g = y[ppw->pv->index_pt_theta_g];
        shear_g = y[ppw->pv->index_pt_shear_g];
      }
      else
      {
        /** - ----> (a.1.2.) radiation streaming approximation */

        delta_g = 0.; /* actual free streaming approximation imposed after evaluation of einstein equations */
        theta_g = 0.; /* actual free streaming approximation imposed after evaluation of einstein equations */
        shear_g = 0.; /* shear always neglected in radiation streaming approximation */
      }
    }
    else
    {
      /** - ----> (a.1.3.) tight coupling approximation */

      delta_g = y[ppw->pv->index_pt_delta_g];
      theta_g = y[ppw->pv->index_pt_theta_g];

      /* first-order tight-coupling approximation for photon shear */
      if (ppt->gauge == newtonian)
        shear_g = 16. / 45. / ppw->pvecthermo[pth->index_th_dkappa] * y[ppw->pv->index_pt_theta_g];
      else
        shear_g = 0.;  /* in the synchronous gauge, the expression of
                        *  shear_g (at first-order in a tight-coupling
                        *  expansion) is a function of h' and eta'; but h'
                        *  and eta' are calculated in perturb_einstein()
                        *  as a function of delta_g and theta_g.  Hence,
                        *  we set shear_g temporarily to zero, and set it
                        *  to the right first-order value in
                        *  perturb_einstein(), just before using the
                        *  Einstein equation for the shear. */
    }

    /** - ---> (a.2.) ur */

    if (pba->has_ur == _TRUE_)
    {
      if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
      {
        delta_ur = y[ppw->pv->index_pt_delta_ur];
        theta_ur = y[ppw->pv->index_pt_theta_ur];
        shear_ur = y[ppw->pv->index_pt_shear_ur];
      }
      else
      {
        delta_ur = 0.; /* actual free streaming approximation imposed after evaluation of 1st einstein equation */
        theta_ur = 0.; /* actual free streaming approximation imposed after evaluation of 1st einstein equation */
        shear_ur = 0.; /* shear always neglected in free streaming approximation */
      }
    }

    /** - --> (b) compute the total density, velocity and shear perturbations */

    /* photon and baryon contribution */
    ppw->delta_rho = ppw->pvecback[pba->index_bg_rho_g] * delta_g
                     + ppw->pvecback[pba->index_bg_rho_b] * y[ppw->pv->index_pt_delta_b];
    ppw->rho_plus_p_theta = 4. / 3. * ppw->pvecback[pba->index_bg_rho_g] * theta_g
                            + ppw->pvecback[pba->index_bg_rho_b] * y[ppw->pv->index_pt_theta_b];
    ppw->rho_plus_p_shear = 4. / 3. * ppw->pvecback[pba->index_bg_rho_g] * shear_g;
    ppw->delta_p          = 1. / 3. * ppw->pvecback[pba->index_bg_rho_g] * delta_g
                            + ppw->pvecthermo[pth->index_th_cb2] * ppw->pvecback[pba->index_bg_rho_b] * y[ppw->pv->index_pt_delta_b];
    rho_plus_p_tot = 4. / 3. * ppw->pvecback[pba->index_bg_rho_g] + ppw->pvecback[pba->index_bg_rho_b];

    /* cdm contribution */
    if (pba->has_cdm == _TRUE_)
    {
      ppw->delta_rho = ppw->delta_rho + ppw->pvecback[pba->index_bg_rho_cdm] * y[ppw->pv->index_pt_delta_cdm];

      if (ppt->gauge == newtonian)
        ppw->rho_plus_p_theta = ppw->rho_plus_p_theta + ppw->pvecback[pba->index_bg_rho_cdm] * y[ppw->pv->index_pt_theta_cdm];

      rho_plus_p_tot += ppw->pvecback[pba->index_bg_rho_cdm];
    }

    /* ultra-relativistic neutrino/relics contribution */

    if (pba->has_ur == _TRUE_)
    {
      ppw->delta_rho        = ppw->delta_rho + ppw->pvecback[pba->index_bg_rho_ur] * delta_ur;
      ppw->rho_plus_p_theta = ppw->rho_plus_p_theta + 4. / 3. * ppw->pvecback[pba->index_bg_rho_ur] * theta_ur;
      ppw->rho_plus_p_shear = ppw->rho_plus_p_shear + 4. / 3. * ppw->pvecback[pba->index_bg_rho_ur] * shear_ur;
      ppw->delta_p         += 1. / 3. * ppw->pvecback[pba->index_bg_rho_ur] * delta_ur;
      rho_plus_p_tot       += 4. / 3. * ppw->pvecback[pba->index_bg_rho_ur];
    }

    /* non-cold dark matter contribution */
    if (pba->has_ncdm == _TRUE_)
    {
      idx = ppw->pv->index_pt_psi0_ncdm1;

      if (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on)
      {
        /* The perturbations are evolved integrated: */
        for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
        {
          rho_ncdm_bg   = ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm];
          p_ncdm_bg     = ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm];
          pseudo_p_ncdm = ppw->pvecback[pba->index_bg_pseudo_p_ncdm1 + n_ncdm];

          rho_plus_p_ncdm = rho_ncdm_bg + p_ncdm_bg;
          w_ncdm          = p_ncdm_bg / rho_ncdm_bg;
          cg2_ncdm        = w_ncdm * (1.0 - 1.0 / (3.0 + 3.0 * w_ncdm) * (3.0 * w_ncdm - 2.0 + pseudo_p_ncdm / p_ncdm_bg));

          if ((ppt->has_source_delta_ncdm == _TRUE_) || (ppt->has_source_theta_ncdm == _TRUE_) || (ppt->has_source_delta_m == _TRUE_))
          {
            ppw->delta_ncdm[n_ncdm] = y[idx];
            ppw->theta_ncdm[n_ncdm] = y[idx + 1];
            ppw->shear_ncdm[n_ncdm] = y[idx + 2];
          }

          ppw->delta_rho        += rho_ncdm_bg * y[idx];
          ppw->rho_plus_p_theta += rho_plus_p_ncdm * y[idx + 1];
          ppw->rho_plus_p_shear += rho_plus_p_ncdm * y[idx + 2];
          ppw->delta_p          += cg2_ncdm * rho_ncdm_bg * y[idx];
          rho_plus_p_tot        += rho_plus_p_ncdm;

          idx += ppw->pv->l_max_ncdm[n_ncdm] + 1;
        }
      }
      else
      {
        /* We must integrate to find perturbations: */
        for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
        {
          rho_delta_ncdm        = 0.0;
          rho_plus_p_theta_ncdm = 0.0;
          rho_plus_p_shear_ncdm = 0.0;
          delta_p_ncdm          = 0.0;
          factor                = pba->factor_ncdm[n_ncdm] * pow (pba->a_today / a, 4);

          for (index_q = 0; index_q < ppw->pv->q_size_ncdm[n_ncdm]; index_q++)
          {
            q       = pba->q_ncdm[n_ncdm][index_q];
            q2      = q * q;
            epsilon = sqrt (q2 + pba->M_ncdm[n_ncdm] * pba->M_ncdm[n_ncdm] * a2);

            rho_delta_ncdm        += q2 * epsilon * pba->w_ncdm[n_ncdm][index_q] * y[idx];
            rho_plus_p_theta_ncdm += q2 * q * pba->w_ncdm[n_ncdm][index_q] * y[idx + 1];
            rho_plus_p_shear_ncdm += q2 * q2 / epsilon * pba->w_ncdm[n_ncdm][index_q] * y[idx + 2];
            delta_p_ncdm          += q2 * q2 / epsilon * pba->w_ncdm[n_ncdm][index_q] * y[idx];

            /*Jump to next momentum bin: */
            idx += (ppw->pv->l_max_ncdm[n_ncdm] + 1);
          }

          rho_delta_ncdm        *= factor;
          rho_plus_p_theta_ncdm *= k * factor;
          rho_plus_p_shear_ncdm *= 2.0 / 3.0 * factor;
          delta_p_ncdm          *= factor / 3.;

          if ((ppt->has_source_delta_ncdm == _TRUE_) || (ppt->has_source_theta_ncdm == _TRUE_) || (ppt->has_source_delta_m == _TRUE_))
          {
            ppw->delta_ncdm[n_ncdm] = rho_delta_ncdm / ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm];
            ppw->theta_ncdm[n_ncdm] = rho_plus_p_theta_ncdm /
                                      (ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm] + ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm]);
            ppw->shear_ncdm[n_ncdm] = rho_plus_p_shear_ncdm /
                                      (ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm] + ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm]);
          }

          ppw->delta_rho        += rho_delta_ncdm;
          ppw->rho_plus_p_theta += rho_plus_p_theta_ncdm;
          ppw->rho_plus_p_shear += rho_plus_p_shear_ncdm;
          ppw->delta_p          += delta_p_ncdm;
          rho_plus_p_tot        += ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm] + ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm];
        }
      }
    }

    /* add your extra species here */

    /* fluid contribution */
    if (pba->has_fld == _TRUE_)
    {
      class_call (background_w_fld (pba, a, &w_fld, &dw_over_da_fld, &integral_fld), pba->error_message, ppt->error_message);

      if (pba->use_ppf == _FALSE_)
      {
        ppw->delta_rho_fld        = ppw->pvecback[pba->index_bg_rho_fld] * y[ppw->pv->index_pt_delta_fld];
        ppw->rho_plus_p_theta_fld = (1. + w_fld) * ppw->pvecback[pba->index_bg_rho_fld] * y[ppw->pv->index_pt_theta_fld];
      }
      else
      {
        s2sq = ppw->s_l[2] * ppw->s_l[2];

        if (ppt->gauge == synchronous)
          alpha = (y[ppw->pv->index_pt_eta] + 1.5 * a2 / k2 / s2sq * (ppw->delta_rho + a_prime_over_a / k2 * ppw->rho_plus_p_theta) - y[ppw->pv->index_pt_Gamma_fld]) / a_prime_over_a;
        else
          alpha = 0.;

        ppw->S_fld = ppw->pvecback[pba->index_bg_rho_fld] * (1. + w_fld) * 1.5 * a2 / k2 / a_prime_over_a *
                     (ppw->rho_plus_p_theta / rho_plus_p_tot + k2 * alpha);
        /* note that the last terms in the ratio do not include fld, that's correct, it's the whole point of the PPF scheme */
        c_gamma_k_H_square                    = pow (pba->c_gamma_over_c_fld * k / a_prime_over_a, 2) * pba->cs2_fld;
        ppw->Gamma_prime_fld                  = a_prime_over_a * (ppw->S_fld / (1. + c_gamma_k_H_square) - (1. + c_gamma_k_H_square) * y[ppw->pv->index_pt_Gamma_fld]);
        Gamma_prime_plus_a_prime_over_a_Gamma = ppw->Gamma_prime_fld + a_prime_over_a * y[ppw->pv->index_pt_Gamma_fld];
        /* delta and theta in both gauges gauge: */
        ppw->rho_plus_p_theta_fld = ppw->pvecback[pba->index_bg_rho_fld] * (1. + w_fld) * ppw->rho_plus_p_theta / rho_plus_p_tot -
                                    k2 * 2. / 3. * a_prime_over_a / a2 / (1 + 4.5 * a2 / k2 / s2sq * rho_plus_p_tot) *
                                    (ppw->S_fld - Gamma_prime_plus_a_prime_over_a_Gamma / a_prime_over_a);
        ppw->delta_rho_fld = -2. / 3. * k2 * s2sq / a2 * y[ppw->pv->index_pt_Gamma_fld] - 3 * a_prime_over_a / k2 * ppw->rho_plus_p_theta_fld;
      }

      ppw->delta_rho        += ppw->delta_rho_fld;
      ppw->rho_plus_p_theta += ppw->rho_plus_p_theta_fld;
      ppw->delta_p          += pba->cs2_fld * ppw->delta_rho_fld;
    }

    /* don't add species here, add them before the fluid contribution: because of the PPF scheme that one must be the last one! */

    /* store delta_m in the current gauge. In perturb_einstein, this
     *  will be transformed later on into the gauge-independent variable D
     *  = delta_m - 2H'/H \theta_m/k^2 .  */

    if (ppt->has_source_delta_m == _TRUE_)
    {
      /* include baryons and cold dark matter */

      delta_rho_m = ppw->pvecback[pba->index_bg_rho_b] * y[ppw->pv->index_pt_delta_b];
      rho_m       = ppw->pvecback[pba->index_bg_rho_b];

      if (pba->has_cdm == _TRUE_)
      {
        delta_rho_m += ppw->pvecback[pba->index_bg_rho_cdm] * y[ppw->pv->index_pt_delta_cdm];
        rho_m       += ppw->pvecback[pba->index_bg_rho_cdm];
      }

      /* infer delta_cb */
      if (ppt->has_source_delta_cb)
        ppw->delta_cb = delta_rho_m / rho_m;

      /* include any other species non-relativistic today (like ncdm species) */

      if (pba->has_ncdm == _TRUE_)
      {
        for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
        {
          delta_rho_m += ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm] * ppw->delta_ncdm[n_ncdm];
          rho_m       += ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm];
        }
      }

      /* infer delta_m */

      ppw->delta_m = delta_rho_m / rho_m;
    }

    /* store theta_m in the current gauge. In perturb_einstein, this
     *  will be transformed later on into the gauge-independent variable
     *  Theta . Note that computing theta_m is necessary also if we want
     *  the delta_m source only, because the gauge-invariant delta_m
     *  involves theta_m in the current gauge. */

    if ((ppt->has_source_delta_m == _TRUE_) || (ppt->has_source_theta_m == _TRUE_))
    {
      /* include baryons and cold dark matter */

      rho_plus_p_theta_m = ppw->pvecback[pba->index_bg_rho_b] * y[ppw->pv->index_pt_theta_b];
      rho_plus_p_m       = ppw->pvecback[pba->index_bg_rho_b];

      if (pba->has_cdm == _TRUE_)
      {
        if (ppt->gauge == newtonian)
          rho_plus_p_theta_m += ppw->pvecback[pba->index_bg_rho_cdm] * y[ppw->pv->index_pt_theta_cdm];

        rho_plus_p_m += ppw->pvecback[pba->index_bg_rho_cdm];
      }

      if ((ppt->has_source_delta_cb == _TRUE_) || (ppt->has_source_theta_cb == _TRUE_))
        ppw->theta_cb = rho_plus_p_theta_m / rho_plus_p_m;

      /* include any other species non-relativistic today (like ncdm species) */

      if (pba->has_ncdm == _TRUE_)
      {
        for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
        {
          rho_plus_p_theta_m += (ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm] + ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm]) * ppw->theta_ncdm[n_ncdm];
          rho_plus_p_m       += (ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm] + ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm]);
        }
      }

      /* infer theta_m */

      ppw->theta_m = rho_plus_p_theta_m / rho_plus_p_m;
    }
  }

  /** - for vector modes */

  if (_vectors_)
  {
    ppw->vector_source_pi = 0.;
    ppw->vector_source_v  = 0.;

    /** - --> photon contribution to vector sources: */
    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off) /* if radiation streaming approximation is off */
    {
      if (ppw->approx[ppw->index_ap_tca] == (int) tca_off) /* if tight-coupling approximation is off */
      {
        ppw->vector_source_v += 4. / 3. * a2 * ppw->pvecback[pba->index_bg_rho_g]
                                * (-1. / 4. * _SQRT2_)
                                * (y[ppw->pv->index_pt_delta_g] + 2. * y[ppw->pv->index_pt_delta_g] + y[ppw->pv->index_pt_shear_g]);

        ppw->vector_source_pi += 1. / 3. * a2 * ppw->pvecback[pba->index_bg_rho_g]
                                 * (6. * _SQRT2_ / 5. / sqrt (1. - 2. * K / k / k))
                                 * (4. / 3. / k * y[ppw->pv->index_pt_theta_g] + y[ppw->pv->index_pt_l3_g]);
      }
    }

    /** - --> baryons */
  }

  /** - for tensor modes */

  if (_tensors_)
  {
    ppw->gw_source = 0.0;

    /** - --> photon contribution to gravitational wave source: */
    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off) /* if radiation streaming approximation is off */
      if (ppw->approx[ppw->index_ap_tca] == (int) tca_off) /* if tight-coupling approximation is off */
        ppw->gw_source += (-_SQRT6_ * 4 * a2 * ppw->pvecback[pba->index_bg_rho_g] *
                           (1. / 15. * y[ppw->pv->index_pt_delta_g] +
                            4. / 21. * y[ppw->pv->index_pt_shear_g] +
                            1. / 35. * y[ppw->pv->index_pt_l3_g + 1]));


    /** - --> ur contribution to gravitational wave source: */
    if (ppt->evolve_tensor_ur == _TRUE_)
    {
      rho_relativistic = 0.;

      if (ppt->tensor_method == tm_exact)
        rho_relativistic += ppw->pvecback[pba->index_bg_rho_ur];

      if (ppt->tensor_method == tm_massless_approximation)
      {
        if (pba->has_ur == _TRUE_)
          rho_relativistic += ppw->pvecback[pba->index_bg_rho_ur];

        if (pba->has_ncdm == _TRUE_)
        {
          for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
          {
            /* (3 p_ncdm1) is the "relativistic" contribution to rho_ncdm1 */
            rho_relativistic += 3. * ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm];
          }
        }
      }

      ppw->gw_source += (-_SQRT6_ * 4 * a2 * rho_relativistic *
                         (1. / 15. * y[ppw->pv->index_pt_delta_ur] +
                          4. / 21. * y[ppw->pv->index_pt_shear_ur] +
                          1. / 35. * y[ppw->pv->index_pt_l3_ur + 1]));
    }

    /** - --> ncdm contribution to gravitational wave source: */
    if (ppt->evolve_tensor_ncdm == _TRUE_)
    {
      idx = ppw->pv->index_pt_psi0_ncdm1;

      /* We must integrate to find perturbations: */
      for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
      {
        gwncdm = 0.;

        factor = pba->factor_ncdm[n_ncdm] * pow (pba->a_today / a, 4);

        for (index_q = 0; index_q < ppw->pv->q_size_ncdm[n_ncdm]; index_q++)
        {
          q       = pba->q_ncdm[n_ncdm][index_q];
          q2      = q * q;
          epsilon = sqrt (q2 + pba->M_ncdm[n_ncdm] * pba->M_ncdm[n_ncdm] * a2);

          gwncdm += q2 * q2 / epsilon * pba->w_ncdm[n_ncdm][index_q] * (1. / 15. * y[idx] + 2. / 21. * y[idx + 2] + 1. / 35. * y[idx + 4]);

          /*Jump to next momentum bin: */
          idx += (ppw->pv->l_max_ncdm[n_ncdm] + 1);
        }

        gwncdm *= -_SQRT6_ * 4 * a2 * factor;

        ppw->gw_source += gwncdm;
      }
    }
  }

  return _SUCCESS_;
}

/**
 * Compute the source functions (three terms for temperature, one for
 * E or B modes, etc.)
 *
 * This is one of the few functions in the code which is passed to
 * the generic_integrator() routine. Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed parameters and workspaces are passed through a generic
 * pointer.  generic_integrator() doesn't know the content of this
 * pointer.
 *
 * - the error management is a bit special: errors are not written as
 * usual to pth->error_message, but to a generic error_message passed
 * in the list of arguments.
 *
 * @param tau                      Input: conformal time
 * @param y                        Input: vector of perturbations
 * @param dy                       Input: vector of time derivative of perturbations
 * @param index_tau                Input: index in the array tau_sampling
 * @param parameters_and_workspace Input/Output: in input, all parameters needed by perturb_derivs, in output, source terms
 * @param error_message            Output: error message
 * @return the error status
 */

int
perturb_sources (
  double   tau,
  double   *y,
  double   *dy,
  int      index_tau,
  void     *parameters_and_workspace,
  ErrorMsg error_message
                )
{
  /** Summary: */

  /** - define local variables */

  /* const double H0 = 1.0 / nc_hicosmo_RH_Mpc (((struct perturb_parameters_and_workspace *) parameters_and_workspace)->pba->cosmo); */

  double P;
  int index_type;

  struct perturb_parameters_and_workspace *pppaw;
  struct precision *ppr;
  const struct background *pba;
  struct thermo *pth;
  struct perturbs *ppt;
  int index_md;
  int index_ic;
  int index_k;
  double k;
  double z;
  struct perturb_workspace *ppw;
  double *pvecback;
  double *pvecthermo;
  double *pvecmetric;

  double delta_g;
  double a_prime_over_a = 0.;       /* (a'/a) */
  double a_prime_over_a_prime = 0.; /* (a'/a)' */
  double w_fld, dw_over_da_fld, integral_fld;
  int switch_isw = 1;

  double a_rel;

  /** - rename structure fields (just to avoid heavy notations) */

  pppaw    = parameters_and_workspace;
  ppr      = pppaw->ppr;
  pba      = pppaw->pba;
  pth      = pppaw->pth;
  ppt      = pppaw->ppt;
  index_md = pppaw->index_md;
  index_ic = pppaw->index_ic;
  index_k  = pppaw->index_k;
  k        = pppaw->k;
  ppw      = pppaw->ppw;

  pvecback   = ppw->pvecback;
  pvecthermo = ppw->pvecthermo;
  pvecmetric = ppw->pvecmetric;

  /** - get background/thermo quantities in this point */

  class_call (background_at_tau (pba,
                                 tau,
                                 pba->normal_info,
                                 pba->inter_closeby,
                                 &(ppw->last_index_back),
                                 pvecback),
              pba->error_message,
              error_message);

  z = pba->a_today / pvecback[pba->index_bg_a] - 1.;

  class_call (thermodynamics_at_z (pba,
                                   pth,
                                   z, /* redshift z=1/a-1 */
                                   pth->inter_closeby,
                                   &(ppw->last_index_thermo),
                                   pvecback,
                                   pvecthermo),
              pth->error_message,
              error_message);

  a_rel = ppw->pvecback[pba->index_bg_a] / pba->a_today;

  /* derived background quantities, useful only in synchronous gauge */
  if (ppt->gauge == synchronous)
  {
    a_prime_over_a       = pvecback[pba->index_bg_a] * pvecback[pba->index_bg_H];                                                                        /* (a'/a)=aH */
    a_prime_over_a_prime = pvecback[pba->index_bg_H_prime] * pvecback[pba->index_bg_a] + pow (pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a], 2); /* (a'/a)' = aH'+(aH)^2 */
  }

  /** - for scalars */
  if (_scalars_)
  {
    /** - --> compute metric perturbations */

    class_call (perturb_einstein (ppr,
                                  pba,
                                  pth,
                                  ppt,
                                  index_md,
                                  k,
                                  tau,
                                  y,
                                  ppw),
                ppt->error_message,
                error_message);

    /** - --> compute quantities depending on approximation schemes */

    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    {
      delta_g = ppw->rsa_delta_g;
      P       = 0.;
    }
    else
    {
      delta_g = y[ppw->pv->index_pt_delta_g];

      if (ppw->approx[ppw->index_ap_tca] == (int) tca_on)
        P = 5. * ppw->s_l[2] * ppw->tca_shear_g / 8.;  /* (2.5+0.5+2)shear_g/8 */
      else
        P = (y[ppw->pv->index_pt_pol0_g] + y[ppw->pv->index_pt_pol2_g] + 2. * ppw->s_l[2] * y[ppw->pv->index_pt_shear_g]) / 8.;
    }

    /** - --> for each type, compute source terms */

    /* scalar temperature */
    if (ppt->has_source_t == _TRUE_)
    {
      /* check whether integrated Sachs-Wolf term should be included */
      if ((ppt->switch_eisw == 0) && (z >= ppt->eisw_lisw_split_z))
        switch_isw = 0;

      if ((ppt->switch_lisw == 0) && (z < ppt->eisw_lisw_split_z))
        switch_isw = 0;

      /* newtonian gauge: simplest form, not efficient numerically */

      /*
       *  if (ppt->gauge == newtonian) {
       *  _set_source_(ppt->index_tp_t0) = pvecthermo[pth->index_th_exp_m_kappa] * pvecmetric[ppw->index_mt_phi_prime] + pvecthermo[pth->index_th_g] * delta_g / 4.;
       *  _set_source_(ppt->index_tp_t1) = pvecthermo[pth->index_th_exp_m_kappa] * k* pvecmetric[ppw->index_mt_psi] + pvecthermo[pth->index_th_g] * y[ppw->pv->index_pt_theta_b]/k;
       *  _set_source_(ppt->index_tp_t2) = pvecthermo[pth->index_th_g] * P;
       *  }
       */

      /* newtonian gauge: slightly more complicated form, but more efficient numerically */

      if (ppt->gauge == newtonian)
      {
        _set_source_ (ppt->index_tp_t0) =
          ppt->switch_sw * pvecthermo[pth->index_th_g] * (delta_g / 4. + pvecmetric[ppw->index_mt_psi])
          + switch_isw * (pvecthermo[pth->index_th_g] * (y[ppw->pv->index_pt_phi] - pvecmetric[ppw->index_mt_psi])
                          + pvecthermo[pth->index_th_exp_m_kappa] * 2. * pvecmetric[ppw->index_mt_phi_prime])
          + ppt->switch_dop / k / k * (pvecthermo[pth->index_th_g] * dy[ppw->pv->index_pt_theta_b]
                                       + pvecthermo[pth->index_th_dg] * y[ppw->pv->index_pt_theta_b]);

        _set_source_ (ppt->index_tp_t1) = switch_isw * pvecthermo[pth->index_th_exp_m_kappa] * k * (pvecmetric[ppw->index_mt_psi] - y[ppw->pv->index_pt_phi]);

        _set_source_ (ppt->index_tp_t2) = ppt->switch_pol * pvecthermo[pth->index_th_g] * P;
      }


      /* synchronous gauge: simplest form, not efficient numerically */

      /*
       *  if (ppt->gauge == synchronous) {
       *  _set_source_(ppt->index_tp_t0) = - pvecthermo[pth->index_th_exp_m_kappa] * pvecmetric[ppw->index_mt_h_prime] / 6. + pvecthermo[pth->index_th_g] / 4. * delta_g;
       *  _set_source_(ppt->index_tp_t1) = pvecthermo[pth->index_th_g] * y[ppw->pv->index_pt_theta_b] / k;
       *  _set_source_(ppt->index_tp_t2) = pvecthermo[pth->index_th_exp_m_kappa] * k*k* 2./3. * ppw->s_l[2] * pvecmetric[ppw->index_mt_alpha] + pvecthermo[pth->index_th_g] * P;
       *  }
       */

      /* synchronous gauge: slightly more complicated form, but more efficient numerically */

      if (ppt->gauge == synchronous)
      {
        _set_source_ (ppt->index_tp_t0) =
          ppt->switch_sw * pvecthermo[pth->index_th_g] * (delta_g / 4. + pvecmetric[ppw->index_mt_alpha_prime])
          + switch_isw * (pvecthermo[pth->index_th_g] * (y[ppw->pv->index_pt_eta]
                                                         - pvecmetric[ppw->index_mt_alpha_prime]
                                                         - 2 * a_prime_over_a * pvecmetric[ppw->index_mt_alpha])
                          + pvecthermo[pth->index_th_exp_m_kappa] * 2. * (pvecmetric[ppw->index_mt_eta_prime]
                                                                          - a_prime_over_a_prime * pvecmetric[ppw->index_mt_alpha]
                                                                          - a_prime_over_a * pvecmetric[ppw->index_mt_alpha_prime]))
          + ppt->switch_dop * (pvecthermo[pth->index_th_g] * (dy[ppw->pv->index_pt_theta_b] / k / k + pvecmetric[ppw->index_mt_alpha_prime])
                               + pvecthermo[pth->index_th_dg] * (y[ppw->pv->index_pt_theta_b] / k / k + pvecmetric[ppw->index_mt_alpha]));

        _set_source_ (ppt->index_tp_t1) =
          switch_isw * pvecthermo[pth->index_th_exp_m_kappa] * k * (pvecmetric[ppw->index_mt_alpha_prime]
                                                                    + 2. * a_prime_over_a * pvecmetric[ppw->index_mt_alpha]
                                                                    - y[ppw->pv->index_pt_eta]);

        _set_source_ (ppt->index_tp_t2) =
          ppt->switch_pol * pvecthermo[pth->index_th_g] * P;
      }
    }

    /* scalar polarization */
    if (ppt->has_source_p == _TRUE_)
      /* all gauges. Note that the correct formula for the E source
       *  should have a minus sign, as shown in Hu & White. We put a
       *  plus sign to comply with the 'historical convention'
       *  established in CMBFAST and CAMB. */

      _set_source_ (ppt->index_tp_p) = sqrt (6.) * pvecthermo[pth->index_th_g] * P;


    /* now, non-CMB sources */

    /* Bardeen potential -PHI_H = phi in Newtonian gauge */
    if (ppt->has_source_phi == _TRUE_)
    {
      if (ppt->gauge == newtonian)
        _set_source_ (ppt->index_tp_phi) = y[ppw->pv->index_pt_phi];

      if (ppt->gauge == synchronous)
        _set_source_ (ppt->index_tp_phi) = y[ppw->pv->index_pt_eta] - a_prime_over_a * pvecmetric[ppw->index_mt_alpha];
    }

    /* its derivative phi' */
    if (ppt->has_source_phi_prime == _TRUE_)
    {
      if (ppt->gauge == newtonian)
        _set_source_ (ppt->index_tp_phi_prime) = dy[ppw->pv->index_pt_phi];

      if (ppt->gauge == synchronous)
        _set_source_ (ppt->index_tp_phi_prime) = dy[ppw->pv->index_pt_eta]
                                                 - a_prime_over_a_prime * pvecmetric[ppw->index_mt_alpha]
                                                 - a_prime_over_a * pvecmetric[ppw->index_mt_alpha_prime];
    }

    /* diff of Bardeen potentials PHI_A-PHI_H = psi + phi in newtonian gauge */
    if (ppt->has_source_phi_plus_psi == _TRUE_)
    {
      if (ppt->gauge == newtonian)
        _set_source_ (ppt->index_tp_phi_plus_psi) =
          y[ppw->pv->index_pt_phi] + pvecmetric[ppw->index_mt_psi];

      if (ppt->gauge == synchronous)
        _set_source_ (ppt->index_tp_phi_plus_psi) =
          y[ppw->pv->index_pt_eta] + pvecmetric[ppw->index_mt_alpha_prime];
    }

    /* Bardeen potential PHI_A = psi in newtonian gauge */
    if (ppt->has_source_psi == _TRUE_)
    {
      if (ppt->gauge == newtonian)
        _set_source_ (ppt->index_tp_psi) =
          pvecmetric[ppw->index_mt_psi];

      if (ppt->gauge == synchronous)
        _set_source_ (ppt->index_tp_psi) =
          a_prime_over_a * pvecmetric[ppw->index_mt_alpha] + pvecmetric[ppw->index_mt_alpha_prime];
    }

    /* the metric potentials h and eta in synchronous gauge */
    if (ppt->gauge == synchronous)
    {
      /* cdm is always on in synchronous gauge, see error message above that checks gauge and has_cdm */
      if (ppt->has_source_h == _TRUE_)
        _set_source_ (ppt->index_tp_h) = -2 * y[ppw->pv->index_pt_delta_cdm];

      if (ppt->has_source_h_prime == _TRUE_)
        _set_source_ (ppt->index_tp_h_prime) = pvecmetric[ppw->index_mt_h_prime];

      if (ppt->has_source_eta == _TRUE_)
        _set_source_ (ppt->index_tp_eta) = y[ppw->pv->index_pt_eta];

      if (ppt->has_source_eta_prime == _TRUE_)
        _set_source_ (ppt->index_tp_eta_prime) = dy[ppw->pv->index_pt_eta];
    }

    /* total matter over density (gauge-invariant, defined as in arXiv:1307.1459) */
    if (ppt->has_source_delta_m == _TRUE_)
      _set_source_ (ppt->index_tp_delta_m) = ppw->delta_m;

    /* cdm and baryon over density */
    if (ppt->has_source_delta_cb == _TRUE_)
      _set_source_ (ppt->index_tp_delta_cb) = ppw->delta_cb;

    /* delta_g */
    if (ppt->has_source_delta_g == _TRUE_)
      _set_source_ (ppt->index_tp_delta_g) = delta_g;

    /* delta_baryon */
    if (ppt->has_source_delta_b == _TRUE_)
      _set_source_ (ppt->index_tp_delta_b) = y[ppw->pv->index_pt_delta_b];

    /* delta_cdm */
    if (ppt->has_source_delta_cdm == _TRUE_)
      _set_source_ (ppt->index_tp_delta_cdm) = y[ppw->pv->index_pt_delta_cdm];

    /* delta_fld */
    if (ppt->has_source_delta_fld == _TRUE_)
      _set_source_ (ppt->index_tp_delta_fld) = ppw->delta_rho_fld / pvecback[pba->index_bg_rho_fld];

    /* delta_ur */
    if (ppt->has_source_delta_ur == _TRUE_)
    {
      if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
        _set_source_ (ppt->index_tp_delta_ur) = y[ppw->pv->index_pt_delta_ur];
      else
        _set_source_ (ppt->index_tp_delta_ur) = ppw->rsa_delta_ur;
    }

    /* delta_ncdm1 */
    if (ppt->has_source_delta_ncdm == _TRUE_)
    {
      for (index_type = ppt->index_tp_delta_ncdm1; index_type < ppt->index_tp_delta_ncdm1 + pba->N_ncdm; index_type++)
      {
        _set_source_ (index_type) = ppw->delta_ncdm[index_type - ppt->index_tp_delta_ncdm1];
      }
    }

    /* total velocity (gauge-invariant, defined as in arXiv:1307.1459) */
    if (ppt->has_source_theta_m == _TRUE_)
      _set_source_ (ppt->index_tp_theta_m) = ppw->theta_m;

    /* cdm and baryon velocity */
    if (ppt->has_source_theta_cb == _TRUE_)
      _set_source_ (ppt->index_tp_theta_cb) = ppw->theta_cb;

    /* theta_g */
    if (ppt->has_source_theta_g == _TRUE_)
    {
      if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
        _set_source_ (ppt->index_tp_theta_g) = y[ppw->pv->index_pt_theta_g];
      else
        _set_source_ (ppt->index_tp_theta_g) = ppw->rsa_theta_g;
    }

    /* theta_baryon */
    if (ppt->has_source_theta_b == _TRUE_)
      _set_source_ (ppt->index_tp_theta_b) = y[ppw->pv->index_pt_theta_b];

    /* theta_cdm */
    if (ppt->has_source_theta_cdm == _TRUE_)
      _set_source_ (ppt->index_tp_theta_cdm) = y[ppw->pv->index_pt_theta_cdm];

    /* theta_fld */
    if (ppt->has_source_theta_fld == _TRUE_)
    {
      class_call (background_w_fld (pba, a_rel * pba->a_today, &w_fld, &dw_over_da_fld, &integral_fld), pba->error_message, ppt->error_message);

      _set_source_ (ppt->index_tp_theta_fld) = ppw->rho_plus_p_theta_fld / (1. + w_fld) / pvecback[pba->index_bg_rho_fld];
    }

    /* theta_ur */
    if (ppt->has_source_theta_ur == _TRUE_)
    {
      if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
        _set_source_ (ppt->index_tp_theta_ur) = y[ppw->pv->index_pt_theta_ur];
      else
        _set_source_ (ppt->index_tp_theta_ur) = ppw->rsa_theta_ur;
    }

    /* theta_ncdm1 */
    if (ppt->has_source_theta_ncdm == _TRUE_)
    {
      for (index_type = ppt->index_tp_theta_ncdm1; index_type < ppt->index_tp_theta_ncdm1 + pba->N_ncdm; index_type++)
      {
        _set_source_ (index_type) = ppw->theta_ncdm[index_type - ppt->index_tp_theta_ncdm1];
      }
    }
  }

  /** - for tensors */
  if (_tensors_)
  {
    /** - --> compute quantities depending on approximation schemes */
    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
    {
      if (ppw->approx[ppw->index_ap_tca] == (int) tca_off)
        P = -(1. / 10. * y[ppw->pv->index_pt_delta_g]
              + 2. / 7. * y[ppw->pv->index_pt_shear_g]
              + 3. / 70. * y[ppw->pv->index_pt_delta_g + 4]
              - 3. / 5. * y[ppw->pv->index_pt_pol0_g]
              + 6. / 7. * y[ppw->pv->index_pt_pol2_g]
              - 3. / 70. * y[ppw->pv->index_pt_pol0_g + 4])
            / sqrt (6.);

      else
        P = 2. / 5. * _SQRT6_ * y[ppw->pv->index_pt_gwdot] / ppw->pvecthermo[pth->index_th_dkappa];  /*TBC */
    }
    else
    {
      P = 0.;
    }

    /* tensor temperature */
    if (ppt->has_source_t == _TRUE_)
      _set_source_ (ppt->index_tp_t2) = -y[ppw->pv->index_pt_gwdot] * pvecthermo[pth->index_th_exp_m_kappa] + pvecthermo[pth->index_th_g] * P;

    /* tensor polarization */
    if (ppt->has_source_p == _TRUE_)
      /* Note that the correct formula for the polarization source
       *  should have a minus sign, as shown in Hu & White. We put a
       *  plus sign to comply with the 'historical convention'
       *  established in CMBFAST and CAMB. */

      _set_source_ (ppt->index_tp_p) = sqrt (6.) * pvecthermo[pth->index_th_g] * P;
  }

  return _SUCCESS_;
}

/**
 * When testing the code or a cosmological model, it can be useful to
 * output perturbations at each step of integration (and not just the
 * delta's at each source sampling point, which is achieved simply by
 * asking for matter transfer functions). Then this function can be
 * passed to the generic_evolver routine.
 *
 * By default, instead of passing this function to generic_evolver,
 * one passes a null pointer. Then this function is just not used.
 *
 * @param tau                      Input: conformal time
 * @param y                        Input: vector of perturbations
 * @param dy                       Input: vector of its derivatives (already allocated)
 * @param parameters_and_workspace Input: fixed parameters (e.g. indices)
 * @param error_message            Output: error message
 *
 */

int
perturb_print_variables (double   tau,
                         double   *y,
                         double   *dy,
                         void     *parameters_and_workspace,
                         ErrorMsg error_message
                        )
{
  struct perturb_parameters_and_workspace *pppaw;
  /** Summary: */

  /** - define local variables */
  /* const double H0 = 1.0 / nc_hicosmo_RH_Mpc (((struct perturb_parameters_and_workspace *) parameters_and_workspace)->pba->cosmo); */

  double k;
  int index_md;
  struct precision *ppr;
  const struct background *pba;
  struct thermo *pth;
  struct perturbs *ppt;
  struct perturb_workspace *ppw;
  double *pvecback;
  double *pvecthermo;
  double *pvecmetric;

  double delta_g, theta_g, shear_g, l4_g, pol0_g, pol1_g, pol2_g, pol4_g;
  double delta_b, theta_b;
  double delta_cdm = 0., theta_cdm = 0.;
  double delta_ur = 0., theta_ur = 0., shear_ur = 0., l4_ur = 0.;
  /** - ncdm sector begins */
  int n_ncdm;
  double *delta_ncdm = NULL, *theta_ncdm = NULL, *shear_ncdm = NULL, *delta_p_over_delta_rho_ncdm = NULL;
  double rho_ncdm_bg, p_ncdm_bg, pseudo_p_ncdm, w_ncdm;
  double rho_delta_ncdm = 0.0;
  double rho_plus_p_theta_ncdm = 0.0;
  double rho_plus_p_shear_ncdm = 0.0;
  double delta_p_ncdm = 0.0;
  double factor = 0.0;
  double q, q2, epsilon;
  /** - ncdm sector ends */
  double phi = 0., psi = 0., alpha = 0.;
  double delta_temp = 0., delta_chi = 0.;

  double a, a2;
  int idx, index_q, storeidx;
  double *dataptr;


  /** - rename structure fields (just to avoid heavy notations) */

  pppaw      = parameters_and_workspace;
  k          = pppaw->k;
  index_md   = pppaw->index_md;
  ppr        = pppaw->ppr;
  pba        = pppaw->pba;
  pth        = pppaw->pth;
  ppt        = pppaw->ppt;
  ppw        = pppaw->ppw;
  pvecback   = ppw->pvecback;
  pvecthermo = ppw->pvecthermo;
  pvecmetric = ppw->pvecmetric;

  /** - update background/thermo quantities in this point */

  class_call (background_at_tau (pba,
                                 tau,
                                 pba->normal_info,
                                 pba->inter_closeby,
                                 &(ppw->last_index_back),
                                 pvecback),
              pba->error_message,
              error_message);

  class_call (thermodynamics_at_z (pba,
                                   pth,
                                   1. / pvecback[pba->index_bg_a] - 1.,
                                   pth->inter_closeby,
                                   &(ppw->last_index_thermo),
                                   pvecback,
                                   pvecthermo),
              pth->error_message,
              error_message);

  /** - update metric perturbations in this point */

  class_call (perturb_einstein (ppr,
                                pba,
                                pth,
                                ppt,
                                index_md,
                                k,
                                tau,
                                y,
                                ppw),
              ppt->error_message,
              error_message);

  a  = pvecback[pba->index_bg_a];
  a2 = a * a;

  if (pba->has_ncdm == _TRUE_)
  {
    class_alloc (delta_ncdm, sizeof (double) * pba->N_ncdm, error_message);
    class_alloc (theta_ncdm, sizeof (double) * pba->N_ncdm, error_message);
    class_alloc (shear_ncdm, sizeof (double) * pba->N_ncdm, error_message);
    class_alloc (delta_p_over_delta_rho_ncdm, sizeof (double) * pba->N_ncdm, error_message);
  }

  /** - calculate perturbed recombination */

  if ((ppt->has_perturbed_recombination == _TRUE_) && (ppw->approx[ppw->index_ap_tca] == (int) tca_off))
  {
    delta_temp = y[ppw->pv->index_pt_perturbed_recombination_delta_temp];
    delta_chi  = y[ppw->pv->index_pt_perturbed_recombination_delta_chi];
  }

  /** - for scalar modes */
  if (_scalars_)
  {
    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
    {
      delta_g = y[ppw->pv->index_pt_delta_g];
      theta_g = y[ppw->pv->index_pt_theta_g];
    }
    else
    {
      delta_g = ppw->rsa_delta_g;
      theta_g = ppw->rsa_theta_g;
    }

    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
    {
      if (ppw->approx[ppw->index_ap_tca] == (int) tca_on)
      {
        shear_g = ppw->tca_shear_g;
        /*l3_g = 6./7.*k/pvecthermo[pth->index_th_dkappa]*ppw->tca_shear_g; */
        pol0_g = 2.5 * ppw->tca_shear_g;
        pol1_g = 7. / 12. * 6. / 7. * k / pvecthermo[pth->index_th_dkappa] * ppw->tca_shear_g;
        pol2_g = 0.5 * ppw->tca_shear_g;
        /*pol3_g = 0.25*6./7.*k/pvecthermo[pth->index_th_dkappa]*ppw->tca_shear_g; */
      }
      else
      {
        shear_g = y[ppw->pv->index_pt_shear_g];
        /*l3_g = y[ppw->pv->index_pt_l3_g]; */
        pol0_g = y[ppw->pv->index_pt_pol0_g];
        pol1_g = y[ppw->pv->index_pt_pol1_g];
        pol2_g = y[ppw->pv->index_pt_pol2_g];
        /*pol3_g = y[ppw->pv->index_pt_pol3_g]; */
      }
    }
    else
    {
      shear_g = 0;
      /*l3_g = 0; */
      pol0_g = 0;
      pol1_g = 0;
      pol2_g = 0;
      /*pol3_g = 0.; */
    }

    if (pba->has_ur == _TRUE_)
    {
      if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
      {
        delta_ur = y[ppw->pv->index_pt_delta_ur];
        theta_ur = y[ppw->pv->index_pt_theta_ur];
        shear_ur = y[ppw->pv->index_pt_shear_ur];
      }
      else
      {
        delta_ur = ppw->rsa_delta_ur;
        theta_ur = ppw->rsa_theta_ur;
        shear_ur = 0.;
      }
    }

    delta_b = y[ppw->pv->index_pt_delta_b];
    theta_b = y[ppw->pv->index_pt_theta_b];

    if (pba->has_cdm == _TRUE_)
    {
      delta_cdm = y[ppw->pv->index_pt_delta_cdm];

      if (ppt->gauge == synchronous)
        theta_cdm = 0.;
      else
        theta_cdm = y[ppw->pv->index_pt_theta_cdm];
    }

    /* gravitational potentials */
    if (ppt->gauge == synchronous)
    {
      alpha = pvecmetric[ppw->index_mt_alpha];

      psi = pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a] * alpha + pvecmetric[ppw->index_mt_alpha_prime];
      phi = y[ppw->pv->index_pt_eta] - pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a] * alpha;
    }
    else if (ppt->gauge == newtonian)
    {
      psi = pvecmetric[ppw->index_mt_psi];
      phi = y[ppw->pv->index_pt_phi];
    }
    else
    {
      psi = 0.0;
      phi = 0.0;
    }

    if (pba->has_ncdm == _TRUE_)
    {
      /** - --> Get delta, deltaP/rho, theta, shear and store in array */
      idx = ppw->pv->index_pt_psi0_ncdm1;

      if (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on)
      {
        /* The perturbations are evolved integrated: */
        for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
        {
          rho_ncdm_bg   = pvecback[pba->index_bg_rho_ncdm1 + n_ncdm];
          p_ncdm_bg     = pvecback[pba->index_bg_p_ncdm1 + n_ncdm];
          pseudo_p_ncdm = pvecback[pba->index_bg_pseudo_p_ncdm1 + n_ncdm];
          w_ncdm        = p_ncdm_bg / rho_ncdm_bg;

          delta_ncdm[n_ncdm] = y[idx];
          theta_ncdm[n_ncdm] = y[idx + 1];
          shear_ncdm[n_ncdm] = y[idx + 2];
          /*This is the adiabatic sound speed: */
          delta_p_over_delta_rho_ncdm[n_ncdm] = w_ncdm * (1.0 - 1.0 / (3.0 + 3.0 * w_ncdm) * (3.0 * w_ncdm - 2.0 + pseudo_p_ncdm / p_ncdm_bg));
          idx                                += ppw->pv->l_max_ncdm[n_ncdm] + 1;
        }
      }
      else
      {
        /* We must integrate to find perturbations: */
        for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
        {
          rho_delta_ncdm        = 0.0;
          rho_plus_p_theta_ncdm = 0.0;
          rho_plus_p_shear_ncdm = 0.0;
          delta_p_ncdm          = 0.0;
          factor                = pba->factor_ncdm[n_ncdm] * pow (pba->a_today / a, 4);

          for (index_q = 0; index_q < ppw->pv->q_size_ncdm[n_ncdm]; index_q++)
          {
            q       = pba->q_ncdm[n_ncdm][index_q];
            q2      = q * q;
            epsilon = sqrt (q2 + pba->M_ncdm[n_ncdm] * pba->M_ncdm[n_ncdm] * a2);

            rho_delta_ncdm        += q2 * epsilon * pba->w_ncdm[n_ncdm][index_q] * y[idx];
            rho_plus_p_theta_ncdm += q2 * q * pba->w_ncdm[n_ncdm][index_q] * y[idx + 1];
            rho_plus_p_shear_ncdm += q2 * q2 / epsilon * pba->w_ncdm[n_ncdm][index_q] * y[idx + 2];
            delta_p_ncdm          += q2 * q2 / epsilon * pba->w_ncdm[n_ncdm][index_q] * y[idx];

            /*Jump to next momentum bin: */
            idx += (ppw->pv->l_max_ncdm[n_ncdm] + 1);
          }

          rho_delta_ncdm        *= factor;
          rho_plus_p_theta_ncdm *= k * factor;
          rho_plus_p_shear_ncdm *= 2.0 / 3.0 * factor;
          delta_p_ncdm          *= factor / 3.;

          delta_ncdm[n_ncdm] = rho_delta_ncdm / ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm];
          theta_ncdm[n_ncdm] = rho_plus_p_theta_ncdm /
                               (ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm] + ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm]);
          shear_ncdm[n_ncdm] = rho_plus_p_shear_ncdm /
                               (ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm] + ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm]);
          delta_p_over_delta_rho_ncdm[n_ncdm] = delta_p_ncdm / rho_delta_ncdm;
        }
      }
    }

    /* converting synchronous variables to newtonian ones */
    if (ppt->gauge == synchronous)
    {
      /* density and velocity perturbations (comment out if you wish to keep synchronous variables) */

      delta_g -= 4. * pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a] * alpha;
      theta_g += k * k * alpha;

      delta_b -= 3. * pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a] * alpha;
      theta_b += k * k * alpha;

      if (pba->has_ur == _TRUE_)
      {
        delta_ur -= 4. * pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a] * alpha;
        theta_ur += k * k * alpha;
      }

      if (pba->has_cdm == _TRUE_)
      {
        delta_cdm -= 3. * pvecback[pba->index_bg_H] * pvecback[pba->index_bg_a] * alpha;
        theta_cdm += k * k * alpha;
      }

      if (pba->has_ncdm == _TRUE_)
      {
        for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
        {
          /** - --> Do gauge transformation of delta, deltaP/rho (?) and theta using -= 3aH(1+w_ncdm) alpha for delta. */
        }
      }
    }

    /*    fprintf(ppw->perturb_output_file," "); */
    /** - --> Handle (re-)allocation */
    if (ppt->scalar_perturbations_data[ppw->index_ikout] == NULL)
    {
      class_alloc (ppt->scalar_perturbations_data[ppw->index_ikout],
                   sizeof (double) * ppt->number_of_scalar_titles,
                   error_message);
      ppt->size_scalar_perturbation_data[ppw->index_ikout] = 0;
    }
    else
    {
      ppt->scalar_perturbations_data[ppw->index_ikout] =
        realloc (ppt->scalar_perturbations_data[ppw->index_ikout],
                 sizeof (double) * (ppt->size_scalar_perturbation_data[ppw->index_ikout] + ppt->number_of_scalar_titles));
    }

    storeidx = 0;
    dataptr  = ppt->scalar_perturbations_data[ppw->index_ikout] +
               ppt->size_scalar_perturbation_data[ppw->index_ikout];
    ppt->size_scalar_perturbation_data[ppw->index_ikout] += ppt->number_of_scalar_titles;

    class_store_double (dataptr, tau, _TRUE_, storeidx);
    class_store_double (dataptr, pvecback[pba->index_bg_a], _TRUE_, storeidx);
    class_store_double (dataptr, delta_g, _TRUE_, storeidx);
    class_store_double (dataptr, theta_g, _TRUE_, storeidx);
    class_store_double (dataptr, shear_g, _TRUE_, storeidx);
    class_store_double (dataptr, pol0_g, _TRUE_, storeidx);
    class_store_double (dataptr, pol1_g, _TRUE_, storeidx);
    class_store_double (dataptr, pol2_g, _TRUE_, storeidx);
    class_store_double (dataptr, delta_b, _TRUE_, storeidx);
    class_store_double (dataptr, theta_b, _TRUE_, storeidx);
    class_store_double (dataptr, psi, _TRUE_, storeidx);
    class_store_double (dataptr, phi, _TRUE_, storeidx);
    /* perturbed recombination */
    class_store_double (dataptr, delta_temp, ppt->has_perturbed_recombination, storeidx);
    class_store_double (dataptr, delta_chi, ppt->has_perturbed_recombination, storeidx);
    /* Ultra relativistic species */
    class_store_double (dataptr, delta_ur, pba->has_ur, storeidx);
    class_store_double (dataptr, theta_ur, pba->has_ur, storeidx);
    class_store_double (dataptr, shear_ur, pba->has_ur, storeidx);
    /* Cold dark matter */
    class_store_double (dataptr, delta_cdm, pba->has_cdm, storeidx);
    class_store_double (dataptr, theta_cdm, pba->has_cdm, storeidx);

    /* Non-cold Dark Matter */
    if ((pba->has_ncdm == _TRUE_) && ((ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_) || (ppt->has_source_delta_m == _TRUE_)))
    {
      for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
      {
        class_store_double (dataptr, delta_ncdm[n_ncdm], _TRUE_, storeidx);
        class_store_double (dataptr, theta_ncdm[n_ncdm], _TRUE_, storeidx);
        class_store_double (dataptr, shear_ncdm[n_ncdm], _TRUE_, storeidx);
        class_store_double (dataptr, delta_p_over_delta_rho_ncdm[n_ncdm],  _TRUE_, storeidx);
      }
    }

    /*fprintf(ppw->perturb_output_file,"\n"); */
  }

  /** - for tensor modes: */

  if (_tensors_)
  {
    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
    {
      if (ppw->approx[ppw->index_ap_tca] == (int) tca_off)
      {
        delta_g = y[ppw->pv->index_pt_delta_g];
        shear_g = y[ppw->pv->index_pt_shear_g];
        l4_g    = y[ppw->pv->index_pt_delta_g + 4];
        pol0_g  = y[ppw->pv->index_pt_pol0_g];
        pol2_g  = y[ppw->pv->index_pt_pol2_g];
        pol4_g  = y[ppw->pv->index_pt_pol0_g + 4];
      }
      else
      {
        delta_g = -4. / 3. * ppw->pv->y[ppw->pv->index_pt_gwdot] / pvecthermo[pth->index_th_dkappa]; /*TBC */
        shear_g = 0.;
        l4_g    = 0.;
        pol0_g  = 1. / 3. * ppw->pv->y[ppw->pv->index_pt_gwdot] / pvecthermo[pth->index_th_dkappa]; /*TBC */
        pol2_g  = 0.;
        pol4_g  = 0.;
      }
    }
    else
    {
      delta_g = 0.;
      shear_g = 0.;
      l4_g    = 0.;
      pol0_g  = 0.;
      pol2_g  = 0.;
      pol4_g  = 0.;
    }

    if (ppt->evolve_tensor_ur == _TRUE_)
    {
      delta_ur = y[ppw->pv->index_pt_delta_ur];
      shear_ur = y[ppw->pv->index_pt_shear_ur];
      l4_ur    = y[ppw->pv->index_pt_delta_ur + 4];
    }

    /** - --> Handle (re-)allocation */
    if (ppt->tensor_perturbations_data[ppw->index_ikout] == NULL)
    {
      class_alloc (ppt->tensor_perturbations_data[ppw->index_ikout],
                   sizeof (double) * ppt->number_of_tensor_titles,
                   error_message);
      ppt->size_tensor_perturbation_data[ppw->index_ikout] = 0;
    }
    else
    {
      ppt->tensor_perturbations_data[ppw->index_ikout] =
        realloc (ppt->tensor_perturbations_data[ppw->index_ikout],
                 sizeof (double) * (ppt->size_tensor_perturbation_data[ppw->index_ikout] + ppt->number_of_tensor_titles));
    }

    storeidx = 0;
    dataptr  = ppt->tensor_perturbations_data[ppw->index_ikout] +
               ppt->size_tensor_perturbation_data[ppw->index_ikout];
    ppt->size_tensor_perturbation_data[ppw->index_ikout] += ppt->number_of_tensor_titles;

    /*fprintf(ppw->perturb_output_file," "); */
    class_store_double (dataptr, tau, _TRUE_, storeidx);
    class_store_double (dataptr, pvecback[pba->index_bg_a], _TRUE_, storeidx);
    class_store_double (dataptr, delta_g, _TRUE_, storeidx);
    class_store_double (dataptr, shear_g, _TRUE_, storeidx);
    class_store_double (dataptr, l4_g, _TRUE_, storeidx);
    class_store_double (dataptr, pol0_g, _TRUE_, storeidx);
    class_store_double (dataptr, pol2_g, _TRUE_, storeidx);
    class_store_double (dataptr, pol4_g, _TRUE_, storeidx);
    class_store_double (dataptr, y[ppw->pv->index_pt_gw], _TRUE_, storeidx);
    class_store_double (dataptr, y[ppw->pv->index_pt_gwdot], _TRUE_, storeidx);

    class_store_double (dataptr, delta_ur, ppt->evolve_tensor_ur, storeidx);
    class_store_double (dataptr, shear_ur, ppt->evolve_tensor_ur, storeidx);
    class_store_double (dataptr, l4_ur, ppt->evolve_tensor_ur, storeidx);
    /*printf("index_pt_delta+ur = %d\n",ppw->pv->index_pt_delta_ur); */

    /* Non-cold Dark Matter */
    if (ppt->evolve_tensor_ncdm == _TRUE_)
    {
      idx = ppw->pv->index_pt_psi0_ncdm1;

      for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
      {
        rho_delta_ncdm        = 0.0;
        rho_plus_p_theta_ncdm = 0.0;
        rho_plus_p_shear_ncdm = 0.0;
        delta_p_ncdm          = 0.0;
        factor                = pba->factor_ncdm[n_ncdm] * pow (pba->a_today / a, 4);

        for (index_q = 0; index_q < ppw->pv->q_size_ncdm[n_ncdm]; index_q++)
        {
          q       = pba->q_ncdm[n_ncdm][index_q];
          q2      = q * q;
          epsilon = sqrt (q2 + pba->M_ncdm[n_ncdm] * pba->M_ncdm[n_ncdm] * a2);

          rho_delta_ncdm        += q2 * epsilon * pba->w_ncdm[n_ncdm][index_q] * y[idx];
          rho_plus_p_theta_ncdm += q2 * q * pba->w_ncdm[n_ncdm][index_q] * y[idx + 1];
          rho_plus_p_shear_ncdm += q2 * q2 / epsilon * pba->w_ncdm[n_ncdm][index_q] * y[idx + 2];
          delta_p_ncdm          += q2 * q2 / epsilon * pba->w_ncdm[n_ncdm][index_q] * y[idx];

          /*Jump to next momentum bin: */
          idx += (ppw->pv->l_max_ncdm[n_ncdm] + 1);
        }

        rho_delta_ncdm        *= factor;
        rho_plus_p_theta_ncdm *= k * factor;
        rho_plus_p_shear_ncdm *= 2.0 / 3.0 * factor;
        delta_p_ncdm          *= factor / 3.;

        delta_ncdm[n_ncdm] = rho_delta_ncdm / ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm];
        theta_ncdm[n_ncdm] = rho_plus_p_theta_ncdm /
                             (ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm] + ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm]);
        shear_ncdm[n_ncdm] = rho_plus_p_shear_ncdm /
                             (ppw->pvecback[pba->index_bg_rho_ncdm1 + n_ncdm] + ppw->pvecback[pba->index_bg_p_ncdm1 + n_ncdm]);

        class_store_double (dataptr, delta_ncdm[n_ncdm], _TRUE_, storeidx);
        class_store_double (dataptr, theta_ncdm[n_ncdm], _TRUE_, storeidx);
        class_store_double (dataptr, shear_ncdm[n_ncdm], _TRUE_, storeidx);
      }
    }

    /*    fprintf(ppw->perturb_output_file,"\n"); */
  }

  if (pba->has_ncdm == _TRUE_)
  {
    free (delta_ncdm);
    free (theta_ncdm);
    free (shear_ncdm);
    free (delta_p_over_delta_rho_ncdm);
  }

  return _SUCCESS_;
}

/**
 * Compute derivative of all perturbations to be integrated
 *
 * For each mode (scalar/vector/tensor) and each wavenumber k, this
 * function computes the derivative of all values in the vector of
 * perturbed variables to be integrated.
 *
 * This is one of the few functions in the code which is passed to the generic_integrator() routine.
 * Since generic_integrator() should work with functions passed from various modules, the format of the arguments
 * is a bit special:
 * - fixed parameters and workspaces are passed through a generic pointer.
 *   generic_integrator() doesn't know what the content of this pointer is.
 * - errors are not written as usual in pth->error_message, but in a generic
 *   error_message passed in the list of arguments.
 *
 * @param tau                      Input: conformal time
 * @param y                        Input: vector of perturbations
 * @param dy                       Output: vector of its derivatives (already allocated)
 * @param parameters_and_workspace Input/Output: in input, fixed parameters (e.g. indices); in output, background and thermo quantities evaluated at tau.
 * @param error_message            Output: error message
 */

int
perturb_derivs (double   tau,
                double   *y,
                double   *dy,
                void     *parameters_and_workspace,
                ErrorMsg error_message
               )
{
  /** Summary: */

  /** - define local variables */

  const double T_cmb    = nc_hicosmo_T_gamma0 (((struct perturb_parameters_and_workspace *) parameters_and_workspace)->pba->cosmo);
  const double Omega0_b = nc_hicosmo_Omega_b0 (((struct perturb_parameters_and_workspace *) parameters_and_workspace)->pba->cosmo);
  int sgnK;
  double K = nc_hicosmo_curvature_K (((struct perturb_parameters_and_workspace *) parameters_and_workspace)->pba->cosmo, &sgnK);

  /* multipole */
  int l;

  /* scale factor and other background quantities */
  double a, a2, a_prime_over_a, R;

  /* short-cut names for the fields of the input structure */
  struct perturb_parameters_and_workspace *pppaw;
  double k, k2;
  int index_md;
  struct precision *ppr;
  const struct background *pba;
  struct thermo *pth;
  struct perturbs *ppt;
  struct perturb_workspace *ppw;
  double *pvecback;
  double *pvecthermo;
  double *pvecmetric;
  double *s_l;
  struct perturb_vector *pv;

  /* short-cut notations for the perturbations */
  double delta_g = 0., theta_g = 0., shear_g = 0.;
  double delta_b, theta_b;
  double cb2, cs2, ca2;
  double metric_continuity = 0., metric_euler = 0., metric_shear = 0., metric_ufa_class = 0.;

  /* perturbed recombination (just to simplify the notation) */

  double H0 = 0., Nnow = 0., n_H = 0., fHe = 0.;
  double delta_temp = 0., delta_chi = 0., chi = 0.;
  double alpha_rec = 0., delta_alpha_rec = 0.;
  double a_rad = 0., Compton_CR = 0.;
  double Tb_in_K = 0.;


  /* Non-metric source terms for photons, i.e. \mathcal{P}^{(m)} from arXiv:1305.3261  */
  double P0, P1, P2;

  /* for use with fluid (fld): */
  double w_fld, dw_over_da_fld, w_prime_fld, integral_fld;

  /* for use with non-cold dark matter (ncdm): */
  int index_q, n_ncdm, idx;
  double q, epsilon, dlnf0_dlnq, qk_div_epsilon;
  double rho_ncdm_bg, p_ncdm_bg, pseudo_p_ncdm, w_ncdm, ca2_ncdm, ceff2_ncdm = 0., cvis2_ncdm = 0.;

  /* for use with curvature */
  double cotKgen, sqrt_absK;
  double s2_squared, ssqrt3;

  /** - rename the fields of the input structure (just to avoid heavy notations) */

  pppaw = parameters_and_workspace;

  k        = pppaw->k;
  k2       = k * k;
  index_md = pppaw->index_md;
  ppr      = pppaw->ppr;
  pba      = pppaw->pba;
  pth      = pppaw->pth;
  ppt      = pppaw->ppt;
  ppw      = pppaw->ppw;

  s_l        = ppw->s_l;
  pvecback   = ppw->pvecback;
  pvecthermo = ppw->pvecthermo;
  pvecmetric = ppw->pvecmetric;
  pv         = ppw->pv;

  /** - get background/thermo quantities in this point */

  class_call (background_at_tau (pba,
                                 tau,
                                 pba->normal_info,
                                 pba->inter_closeby,
                                 &(ppw->last_index_back),
                                 pvecback),
              pba->error_message,
              error_message);

  class_call (thermodynamics_at_z (pba,
                                   pth,
                                   1. / pvecback[pba->index_bg_a] - 1., /* redshift z=1/a-1 */
                                   pth->inter_closeby,
                                   &(ppw->last_index_thermo),
                                   pvecback,
                                   pvecthermo),
              pth->error_message,
              error_message);

  /** - get metric perturbations with perturb_einstein() */
  class_call (perturb_einstein (ppr,
                                pba,
                                pth,
                                ppt,
                                index_md,
                                k,
                                tau,
                                y,
                                ppw),
              ppt->error_message,
              error_message);

  /** - compute related background quantities */

  a              = pvecback[pba->index_bg_a];
  a2             = a * a;
  a_prime_over_a = pvecback[pba->index_bg_H] * a;
  R              = 4. / 3. * pvecback[pba->index_bg_rho_g] / pvecback[pba->index_bg_rho_b];

  /** - Compute 'generalised cotK function of argument \f$ \sqrt{|K|}*\tau \f$, for closing hierarchy.
   *   (see equation 2.34 in arXiv:1305.3261): */
  if (pba->has_curvature == _FALSE_)
  {
    cotKgen = 1.0 / (k * tau);
  }
  else
  {
    sqrt_absK = sqrt (fabs (K));

    if (K < 0)
      cotKgen = sqrt_absK / k / tanh (sqrt_absK * tau);
    else
      cotKgen = sqrt_absK / k / tan (sqrt_absK * tau);
  }

  s2_squared = 1. - 3. * K / k2;

  /** - for scalar modes: */
  if (_scalars_)
  {
    /** - --> (a) define short-cut notations for the scalar perturbations */
    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
    {
      delta_g = y[pv->index_pt_delta_g];
      theta_g = y[pv->index_pt_theta_g];
    }

    delta_b = y[pv->index_pt_delta_b];
    theta_b = y[pv->index_pt_theta_b];
    cb2     = pvecthermo[pth->index_th_cb2];

    /** - --> (b) perturbed recombination **/

    if ((ppt->has_perturbed_recombination == _TRUE_) && (ppw->approx[ppw->index_ap_tca] == (int) tca_off))
    {
      delta_temp = y[ppw->pv->index_pt_perturbed_recombination_delta_temp];
      delta_chi  = y[ppw->pv->index_pt_perturbed_recombination_delta_chi];
      chi        = pvecthermo[pth->index_th_xe];

      /* Conversion of H0 in inverse seconds (H0 is [H0/c] in inverse Mpcs) */
      H0 = H0 * _c_ / _Mpc_over_m_;

      /*Computation of Nnow in SI units */
      Nnow = 3. * H0 * H0 * Omega0_b * (1. - pth->YHe) / (8. * _PI_ * _G_ * _m_H_);

      /* total amount of hydrogen today */
      n_H = (pba->a_today / a) * (pba->a_today / a) * (pba->a_today / a) * Nnow;

      /* Helium-to-hydrogen ratio */
      fHe = pth->YHe / (_not4_ * (1 - pth->YHe));

      /* The constant such that rho_gamma = a_rad * T^4 */
      a_rad = 8. / 15. * pow (_PI_, 5) * pow (_k_B_, 4) / pow (_c_ * _h_P_, 3);

      /* Compton cooling rate in Mpc^(-1) */
      Compton_CR = 8. / 3. * _sigma_ * a_rad / (_m_e_ * _c_ * _c_) * _Mpc_over_m_;

      /* Temperature is already in Kelvin */
      Tb_in_K = pvecthermo[pth->index_th_Tb];

      /* Alpha in m^3/s, cf. Recfast paper */
      alpha_rec = 1.14 * 4.309e-19 * pow ((Tb_in_K * 1e-4), -0.6166) / (1 + 0.6703 * pow ((Tb_in_K * 1e-4), 0.53));

      /* delta alpha, dimensionless */
      delta_alpha_rec = (-0.6166 + 0.6703 * pow ((Tb_in_K * 1e-4), 0.53) * (-0.6166 - 0.53)) / (1 + 0.6703 * pow ((Tb_in_K * 1e-4), 0.53)) * delta_temp;
    } /* end of perturbed recombination related quantities */

    /** - --> (c) compute metric-related quantities (depending on gauge; additional gauges can be coded below)
     *
     *       - Each continuity equation contains a term in (theta+metric_continuity) with
     *        metric_continuity = (h_prime/2) in synchronous gauge, (-3 phi_prime) in newtonian gauge
     *
     *       - Each Euler equation contains a source term metric_euler with
     *        metric_euler = 0 in synchronous gauge, (k2 psi) in newtonian gauge
     *
     *       - Each shear derivative equation contains a source term metric_shear equal to
     *        metric_shear = (h_prime+6eta_prime)/2 in synchronous gauge, 0 in newtonian gauge
     *
     *       - metric_shear_prime is the derivative of metric_shear
     *
     *       - In the ufa_class approximation, the leading-order source term is (h_prime/2) in synchronous gauge,
     *        (-3 (phi_prime+psi_prime)) in newtonian gauge: we approximate the later by (-6 phi_prime) */

    if (ppt->gauge == synchronous)
    {
      metric_continuity = pvecmetric[ppw->index_mt_h_prime] / 2.;
      metric_euler      = 0.;
      metric_shear      = k2 * pvecmetric[ppw->index_mt_alpha];
      /*metric_shear_prime = k2 * pvecmetric[ppw->index_mt_alpha_prime]; */
      metric_ufa_class = pvecmetric[ppw->index_mt_h_prime] / 2.;
    }

    if (ppt->gauge == newtonian)
    {
      metric_continuity = -3. * pvecmetric[ppw->index_mt_phi_prime];
      metric_euler      = k2 * pvecmetric[ppw->index_mt_psi];
      metric_shear      = 0.;
      /*metric_shear_prime = 0.; */
      metric_ufa_class = -6. * pvecmetric[ppw->index_mt_phi_prime];
    }

    /** - --> (d) if some approximation schemes are turned on, enforce a few y[] values computed in perturb_einstein */

    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    {
      delta_g = ppw->rsa_delta_g;
      theta_g = ppw->rsa_theta_g;
    }

    /** - --> (e) BEGINNING OF ACTUAL SYSTEM OF EQUATIONS OF EVOLUTION */

    /* Note concerning perturbed recombination: $cb2*delta_b$ must be replaced everywhere by $cb2*(delta_b+delta_temp)$. If perturbed recombination is not required, delta_temp is equal to zero. */

    /** - ---> photon temperature density */

    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
      dy[pv->index_pt_delta_g] = -4. / 3. * (theta_g + metric_continuity);


    /** - ---> baryon density */

    dy[pv->index_pt_delta_b] = -(theta_b + metric_continuity);

    /** - ---> baryon velocity (depends on tight-coupling approximation=tca) */

    if (ppw->approx[ppw->index_ap_tca] == (int) tca_off)
    {
      /* without tca */

      /** - ----> perturbed recombination has an impact **/
      dy[pv->index_pt_theta_b] =
        -a_prime_over_a * theta_b
        + metric_euler
        + k2 * cb2 * (delta_b + delta_temp)
        + R * pvecthermo[pth->index_th_dkappa] * (theta_g - theta_b);
    }
    else
    {
      /* with tca */
      class_call (perturb_tca_slip_and_shear (y, pppaw, error_message),
                  error_message,
                  error_message);

      /* perturbed recombination has an impact **/
      dy[pv->index_pt_theta_b] =
        (-a_prime_over_a * theta_b
         + k2 * (cb2 * (delta_b + delta_temp) + R * (delta_g / 4. - s2_squared * ppw->tca_shear_g))
         + R * ppw->tca_slip) / (1. + R)
        + metric_euler;
    }

    /** - ---> photon temperature higher momenta and photon polarization (depend on tight-coupling approximation) */

    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
    {
      /** - ----> if photon tight-coupling is off */
      if (ppw->approx[ppw->index_ap_tca] == (int) tca_off)
      {
        /** - -----> define \f$ \Pi = G_{\gamma 0} + G_{\gamma 2} + F_{\gamma 2} \f$ */
        P0 = (y[pv->index_pt_pol0_g] + y[pv->index_pt_pol2_g] + 2. * s_l[2] * y[pv->index_pt_shear_g]) / 8.;

        /** - -----> photon temperature velocity */

        dy[pv->index_pt_theta_g] =
          k2 * (delta_g / 4. - s2_squared * y[pv->index_pt_shear_g])
          + metric_euler
          + pvecthermo[pth->index_th_dkappa] * (theta_b - theta_g);

        /** - -----> photon temperature shear */
        dy[pv->index_pt_shear_g] =
          0.5 * (8. / 15. * (theta_g + metric_shear)
                 - 3. / 5. * k * s_l[3] / s_l[2] * y[pv->index_pt_l3_g]
                 - pvecthermo[pth->index_th_dkappa] * (2. * y[pv->index_pt_shear_g] - 4. / 5. / s_l[2] * P0));

        /** - -----> photon temperature l=3 */

        l                     = 3;
        dy[pv->index_pt_l3_g] = k / (2.0 * l + 1.0) *
                                (l * s_l[l] * 2. * s_l[2] * y[pv->index_pt_shear_g] - (l + 1.) * s_l[l + 1] * y[pv->index_pt_l3_g + 1])
                                - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_l3_g];

        /** - -----> photon temperature l>3 */
        for (l = 4; l < pv->l_max_g; l++)
        {
          dy[pv->index_pt_delta_g + l] = k / (2.0 * l + 1.0) *
                                         (l * s_l[l] * y[pv->index_pt_delta_g + l - 1] - (l + 1) * s_l[l + 1] * y[pv->index_pt_delta_g + l + 1])
                                         - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_delta_g + l];
        }

        /** - -----> photon temperature lmax */
        l                            = pv->l_max_g; /* l=lmax */
        dy[pv->index_pt_delta_g + l] =
          k * (s_l[l] * y[pv->index_pt_delta_g + l - 1] - (1. + l) * cotKgen * y[pv->index_pt_delta_g + l])
          - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_delta_g + l];

        /** - -----> photon polarization l=0 */

        dy[pv->index_pt_pol0_g] =
          -k * y[pv->index_pt_pol0_g + 1]
          - pvecthermo[pth->index_th_dkappa] * (y[pv->index_pt_pol0_g] - 4. * P0);

        /** - -----> photon polarization l=1 */

        dy[pv->index_pt_pol1_g] =
          k / 3. * (y[pv->index_pt_pol1_g - 1] - 2. * s_l[2] * y[pv->index_pt_pol1_g + 1])
          - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_pol1_g];

        /** - -----> photon polarization l=2 */

        dy[pv->index_pt_pol2_g] =
          k / 5. * (2. * s_l[2] * y[pv->index_pt_pol2_g - 1] - 3. * s_l[3] * y[pv->index_pt_pol2_g + 1])
          - pvecthermo[pth->index_th_dkappa] * (y[pv->index_pt_pol2_g] - 4. / 5. * P0);

        /** - -----> photon polarization l>2 */

        for (l = 3; l < pv->l_max_pol_g; l++)
          dy[pv->index_pt_pol0_g + l] = k / (2. * l + 1) *
                                        (l * s_l[l] * y[pv->index_pt_pol0_g + l - 1] - (l + 1.) * s_l[l + 1] * y[pv->index_pt_pol0_g + l + 1])
                                        - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_pol0_g + l];

        /** - -----> photon polarization lmax_pol */

        l                           = pv->l_max_pol_g;
        dy[pv->index_pt_pol0_g + l] =
          k * (s_l[l] * y[pv->index_pt_pol0_g + l - 1] - (l + 1) * cotKgen * y[pv->index_pt_pol0_g + l])
          - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_pol0_g + l];
      }
      /** - ----> if photon tight-coupling is on: */

      else
      {
        /** - -----> in that case, only need photon velocity */


        /* perturbed recombination has an impact **/
        dy[pv->index_pt_theta_g] =
          -(dy[pv->index_pt_theta_b] + a_prime_over_a * theta_b - cb2 * k2 * (delta_b + delta_temp)) / R
          + k2 * (0.25 * delta_g - s2_squared * ppw->tca_shear_g) + (1. + R) / R * metric_euler;
      }
    }

    /** - ---> cdm */

    if (pba->has_cdm == _TRUE_)
    {
      /** - ----> newtonian gauge: cdm density and velocity */

      if (ppt->gauge == newtonian)
      {
        dy[pv->index_pt_delta_cdm] = -(y[pv->index_pt_theta_cdm] + metric_continuity); /* cdm density */

        dy[pv->index_pt_theta_cdm] = -a_prime_over_a * y[pv->index_pt_theta_cdm] + metric_euler; /* cdm velocity */
      }

      /** - ----> synchronous gauge: cdm density only (velocity set to zero by definition of the gauge) */

      if (ppt->gauge == synchronous)
        dy[pv->index_pt_delta_cdm] = -metric_continuity;  /* cdm density */
    }

    /* perturbed recombination */
    /* computes the derivatives of delta x_e and delta T_b */

    if ((ppt->has_perturbed_recombination == _TRUE_) && (ppw->approx[ppw->index_ap_tca] == (int) tca_off))
    {
      /* alpha * n_H is in inverse seconds, so we have to multiply it by Mpc_in_sec */
      dy[ppw->pv->index_pt_perturbed_recombination_delta_chi] = -alpha_rec * a * chi * n_H  * (delta_alpha_rec + delta_chi + delta_b) * _Mpc_over_m_ / _c_;

      /* see the documentation for this formula */
      dy[ppw->pv->index_pt_perturbed_recombination_delta_temp] =  2. / 3. * dy[ppw->pv->index_pt_delta_b] - a * Compton_CR * pow (T_cmb / a, 4) * chi / (1. + chi + fHe) * ((1. - T_cmb * pba->a_today / a / pvecthermo[pth->index_th_Tb]) * (delta_g + delta_chi * (1. + fHe) / (1. + chi + fHe)) + T_cmb * pba->a_today / a / pvecthermo[pth->index_th_Tb] * (delta_temp - 1. / 4. * delta_g));
    }

    /** - ---> fluid (fld) */

    if (pba->has_fld == _TRUE_)
    {
      if (pba->use_ppf == _FALSE_)
      {
        /** - ----> factors w, w_prime, adiabatic sound speed ca2 (all three background-related),
         *   plus actual sound speed in the fluid rest frame cs2 */

        class_call (background_w_fld (pba, a, &w_fld, &dw_over_da_fld, &integral_fld), pba->error_message, ppt->error_message);
        w_prime_fld = dw_over_da_fld * a_prime_over_a * a;

        ca2 = w_fld - w_prime_fld / 3. / (1. + w_fld) / a_prime_over_a;
        cs2 = pba->cs2_fld;

        /** - ----> fluid density */

        dy[pv->index_pt_delta_fld] =
          -(1 + w_fld) * (y[pv->index_pt_theta_fld] + metric_continuity)
          - 3. * (cs2 - w_fld) * a_prime_over_a * y[pv->index_pt_delta_fld]
          - 9. * (1 + w_fld) * (cs2 - ca2) * a_prime_over_a * a_prime_over_a * y[pv->index_pt_theta_fld] / k2;

        /** - ----> fluid velocity */

        dy[pv->index_pt_theta_fld] = /* fluid velocity */
                                     -(1. - 3. * cs2) * a_prime_over_a * y[pv->index_pt_theta_fld]
                                     + cs2 * k2 / (1. + w_fld) * y[pv->index_pt_delta_fld]
                                     + metric_euler;
      }
      else
      {
        dy[pv->index_pt_Gamma_fld] = ppw->Gamma_prime_fld; /* Gamma variable of PPF formalism */
      }
    }

    /** - ---> ultra-relativistic neutrino/relics (ur) */

    if (pba->has_ur == _TRUE_)
    {
      /** - ----> if radiation streaming approximation is off */

      if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
      {
        /** - -----> ur density */
        dy[pv->index_pt_delta_ur] =
          /* standard term */
          -4. / 3. * (y[pv->index_pt_theta_ur] + metric_continuity)
          /* non-standard term, non-zero if if ceff2_ur not 1/3 */
          + (1. - ppt->three_ceff2_ur) * a_prime_over_a * (y[pv->index_pt_delta_ur] + 4. * a_prime_over_a * y[pv->index_pt_theta_ur] / k / k);

        /** - -----> ur velocity */
        dy[pv->index_pt_theta_ur] =
          /* standard term with extra coefficient (3 ceff2_ur), normally equal to one */
          k2 * (ppt->three_ceff2_ur * y[pv->index_pt_delta_ur] / 4. - s2_squared * y[pv->index_pt_shear_ur]) + metric_euler
          /* non-standard term, non-zero if ceff2_ur not 1/3 */
          - (1. - ppt->three_ceff2_ur) * a_prime_over_a * y[pv->index_pt_theta_ur];

        if (ppw->approx[ppw->index_ap_ufa] == (int) ufa_off)
        {
          /** - -----> exact ur shear */
          dy[pv->index_pt_shear_ur] =
            0.5 * (
              /* standard term */
              8. / 15. * (y[pv->index_pt_theta_ur] + metric_shear) - 3. / 5. * k * s_l[3] / s_l[2] * y[pv->index_pt_shear_ur + 1]
              /* non-standard term, non-zero if cvis2_ur not 1/3 */
              - (1. - ppt->three_cvis2_ur) * (8. / 15. * (y[pv->index_pt_theta_ur] + metric_shear)));

          /** - -----> exact ur l=3 */
          l                      = 3;
          dy[pv->index_pt_l3_ur] = k / (2. * l + 1.) *
                                   (l * 2. * s_l[l] * s_l[2] * y[pv->index_pt_shear_ur] - (l + 1.) * s_l[l + 1] * y[pv->index_pt_l3_ur + 1]);

          /** - -----> exact ur l>3 */
          for (l = 4; l < pv->l_max_ur; l++)
          {
            dy[pv->index_pt_delta_ur + l] = k / (2. * l + 1) *
                                            (l * s_l[l] * y[pv->index_pt_delta_ur + l - 1] - (l + 1.) * s_l[l + 1] * y[pv->index_pt_delta_ur + l + 1]);
          }

          /** - -----> exact ur lmax_ur */
          l                             = pv->l_max_ur;
          dy[pv->index_pt_delta_ur + l] =
            k * (s_l[l] * y[pv->index_pt_delta_ur + l - 1] - (1. + l) * cotKgen * y[pv->index_pt_delta_ur + l]);
        }
        else
        {
          /** - -----> in fluid approximation (ufa): only ur shear needed */
          /*TBC: curvature? */
          /* a la Ma & Bertschinger */
          if (ppr->ur_fluid_approximation == ufa_mb)
            dy[pv->index_pt_shear_ur] =
              -3. / tau * y[pv->index_pt_shear_ur]
              + 2. / 3. * (y[pv->index_pt_theta_ur] + metric_shear);


          /* a la Hu */
          if (ppr->ur_fluid_approximation == ufa_hu)
            dy[pv->index_pt_shear_ur] =
              -3. * a_prime_over_a * y[pv->index_pt_shear_ur]
              + 2. / 3. * (y[pv->index_pt_theta_ur] + metric_shear);


          /* a la CLASS */
          if (ppr->ur_fluid_approximation == ufa_CLASS)
            dy[pv->index_pt_shear_ur] =
              -3. / tau * y[pv->index_pt_shear_ur]
              + 2. / 3. * (y[pv->index_pt_theta_ur] + metric_ufa_class);
        }
      }
    }

    /** - ---> non-cold dark matter (ncdm): massive neutrinos, WDM, etc. */
    /*TBC: curvature in all ncdm */
    if (pba->has_ncdm == _TRUE_)
    {
      idx = pv->index_pt_psi0_ncdm1;

      /** - ----> first case: use a fluid approximation (ncdmfa) */
      /*TBC: curvature */
      if (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on)
      {
        /** - -----> loop over species */

        for (n_ncdm = 0; n_ncdm < pv->N_ncdm; n_ncdm++)
        {
          /** - -----> define intermediate quantitites */

          rho_ncdm_bg   = pvecback[pba->index_bg_rho_ncdm1 + n_ncdm];                        /* background density */
          p_ncdm_bg     = pvecback[pba->index_bg_p_ncdm1 + n_ncdm];                          /* background pressure */
          pseudo_p_ncdm = pvecback[pba->index_bg_pseudo_p_ncdm1 + n_ncdm];                   /* pseudo-pressure (see CLASS IV paper) */
          w_ncdm        = p_ncdm_bg / rho_ncdm_bg;                                           /* equation of state parameter */
          ca2_ncdm      = w_ncdm / 3.0 / (1.0 + w_ncdm) * (5.0 - pseudo_p_ncdm / p_ncdm_bg); /* adiabatic sound speed */

          /* c_eff is (delta p / delta rho) in the gauge under
           *  consideration (not in the gauge comoving with the
           *  fluid) */

          /* c_vis is introduced in order to close the system */

          /* different ansatz for sound speed c_eff and viscosity speed c_vis */
          if (ppr->ncdm_fluid_approximation == ncdmfa_mb)
          {
            ceff2_ncdm = ca2_ncdm;
            cvis2_ncdm = 3. * w_ncdm * ca2_ncdm;
          }

          if (ppr->ncdm_fluid_approximation == ncdmfa_hu)
          {
            ceff2_ncdm = ca2_ncdm;
            cvis2_ncdm = w_ncdm;
          }

          if (ppr->ncdm_fluid_approximation == ncdmfa_CLASS)
          {
            ceff2_ncdm = ca2_ncdm;
            cvis2_ncdm = 3. * w_ncdm * ca2_ncdm;
          }

          /** - -----> exact continuity equation */

          dy[idx] = -(1.0 + w_ncdm) * (y[idx + 1] + metric_continuity) -
                    3.0 * a_prime_over_a * (ceff2_ncdm - w_ncdm) * y[idx];

          /** - -----> exact euler equation */

          dy[idx + 1] = -a_prime_over_a * (1.0 - 3.0 * ca2_ncdm) * y[idx + 1] +
                        ceff2_ncdm / (1.0 + w_ncdm) * k2 * y[idx] - k2 * y[idx + 2]
                        + metric_euler;

          /** - -----> different ansatz for approximate shear derivative */

          if (ppr->ncdm_fluid_approximation == ncdmfa_mb)
            dy[idx + 2] = -3.0 * (a_prime_over_a * (2. / 3. - ca2_ncdm - pseudo_p_ncdm / p_ncdm_bg / 3.) + 1. / tau) * y[idx + 2]
                          + 8.0 / 3.0 * cvis2_ncdm / (1.0 + w_ncdm) * s_l[2] * (y[idx + 1] + metric_shear);


          if (ppr->ncdm_fluid_approximation == ncdmfa_hu)
            dy[idx + 2] = -3.0 * a_prime_over_a * ca2_ncdm / w_ncdm * y[idx + 2]
                          + 8.0 / 3.0 * cvis2_ncdm / (1.0 + w_ncdm) * s_l[2] * (y[idx + 1] + metric_shear);


          if (ppr->ncdm_fluid_approximation == ncdmfa_CLASS)
            dy[idx + 2] = -3.0 * (a_prime_over_a * (2. / 3. - ca2_ncdm - pseudo_p_ncdm / p_ncdm_bg / 3.) + 1. / tau) * y[idx + 2]
                          + 8.0 / 3.0 * cvis2_ncdm / (1.0 + w_ncdm) * s_l[2] * (y[idx + 1] + metric_ufa_class);


          /** - -----> jump to next species */

          idx += pv->l_max_ncdm[n_ncdm] + 1;
        }
      }
      /** - ----> second case: use exact equation (Boltzmann hierarchy on momentum grid) */

      else
      {
        /** - -----> loop over species */

        for (n_ncdm = 0; n_ncdm < pv->N_ncdm; n_ncdm++)
        {
          /** - -----> loop over momentum */

          for (index_q = 0; index_q < pv->q_size_ncdm[n_ncdm]; index_q++)
          {
            /** - -----> define intermediate quantities */

            dlnf0_dlnq     = pba->dlnf0_dlnq_ncdm[n_ncdm][index_q];
            q              = pba->q_ncdm[n_ncdm][index_q];
            epsilon        = sqrt (q * q + a2 * pba->M_ncdm[n_ncdm] * pba->M_ncdm[n_ncdm]);
            qk_div_epsilon = k * q / epsilon;

            /** - -----> ncdm density for given momentum bin */

            dy[idx] = -qk_div_epsilon * y[idx + 1] + metric_continuity * dlnf0_dlnq / 3.;

            /** - -----> ncdm velocity for given momentum bin */

            dy[idx + 1] = qk_div_epsilon / 3.0 * (y[idx] - 2 * s_l[2] * y[idx + 2])
                          - epsilon * metric_euler / (3 * q * k) * dlnf0_dlnq;

            /** - -----> ncdm shear for given momentum bin */

            dy[idx + 2] = qk_div_epsilon / 5.0 * (2 * s_l[2] * y[idx + 1] - 3. * s_l[3] * y[idx + 3])
                          - s_l[2] * metric_shear * 2. / 15. * dlnf0_dlnq;

            /** - -----> ncdm l>3 for given momentum bin */

            for (l = 3; l < pv->l_max_ncdm[n_ncdm]; l++)
            {
              dy[idx + l] = qk_div_epsilon / (2. * l + 1.0) * (l * s_l[l] * y[idx + (l - 1)] - (l + 1.) * s_l[l + 1] * y[idx + (l + 1)]);
            }

            /** - -----> ncdm lmax for given momentum bin (truncation as in Ma and Bertschinger)
             *   but with curvature taken into account a la arXiv:1305.3261 */

            dy[idx + l] = qk_div_epsilon * y[idx + l - 1] - (1. + l) * k * cotKgen * y[idx + l];

            /** - -----> jump to next momentum bin or species */

            idx += (pv->l_max_ncdm[n_ncdm] + 1);
          }
        }
      }
    }

    /** - ---> metric */

    /** - ---> eta of synchronous gauge */

    if (ppt->gauge == synchronous)
      dy[pv->index_pt_eta] = pvecmetric[ppw->index_mt_eta_prime];


    if (ppt->gauge == newtonian)
      dy[pv->index_pt_phi] = pvecmetric[ppw->index_mt_phi_prime];
  }

  /** - vector mode */
  if (_vectors_)
  {
    fprintf (stderr, "we are in vectors\n");

    ssqrt3 = sqrt (1. - 2. * K / k2);
    cb2    = pvecthermo[pth->index_th_cb2];

    /** - --> baryon velocity */

    if (ppt->gauge == synchronous)
      dy[pv->index_pt_theta_b] = -(1 - 3. * cb2) * a_prime_over_a * y[pv->index_pt_theta_b]
                                 - pvecthermo[pth->index_th_dkappa] * (_SQRT2_ / 4. * delta_g + y[pv->index_pt_theta_b]);


    else if (ppt->gauge == newtonian)
      dy[pv->index_pt_theta_b] = -(1 - 3. * cb2) * a_prime_over_a * y[pv->index_pt_theta_b]
                                 - _SQRT2_ / 4. * pvecthermo[pth->index_th_dkappa] * (delta_g + 2. * _SQRT2_ * y[pv->index_pt_theta_b])
                                 + pvecmetric[ppw->index_mt_V_prime] + (1. - 3. * cb2) * a_prime_over_a * y[pv->index_pt_V];



    /*
     *  if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off) {
     *  if (ppw->approx[ppw->index_ap_tca]==(int)tca_off) {
     */

    /* short-cut notations for the tensor perturbations */
    delta_g = y[pv->index_pt_delta_g];
    theta_g = y[pv->index_pt_theta_g];
    shear_g = y[pv->index_pt_shear_g];


    /* (P^{(1)}) (see Eq. B.23 in 1305.3261)*/
    P1 = -_SQRT6_ / 40. * (
      4. / (3. * k) * theta_g /*F1 */
      + y[pv->index_pt_delta_g + 3]
      + 2. * y[pv->index_pt_pol0_g]
      + 10. / 7. * y[pv->index_pt_pol2_g]
      - 4. / 7. * y[pv->index_pt_pol0_g + 4]);

    if (ppt->gauge == synchronous)
    {
      /* photon density (delta_g = F_0) */
      dy[pv->index_pt_delta_g] =
        -4. / 3. * theta_g
        - pvecthermo[pth->index_th_dkappa] * (delta_g + 2. * _SQRT2_ * y[pv->index_pt_theta_b]);

      /* photon velocity (theta_g = (3k/4)*F_1) */
      dy[pv->index_pt_theta_g] =
        k2 * (delta_g / 4. - s_l[2] * shear_g)
        - pvecthermo[pth->index_th_dkappa] * (theta_g + 4.0 / _SQRT6_ * P1)
        + 4.0 / (3.0 * _SQRT2_) * ssqrt3 * y[pv->index_pt_hv_prime];
    }
    else if (ppt->gauge == newtonian)
    {
      /* photon density (delta_g = F_0) */
      dy[pv->index_pt_delta_g] =
        -4. / 3. * theta_g
        - pvecthermo[pth->index_th_dkappa] * (delta_g + 2. * _SQRT2_ * y[pv->index_pt_theta_b])
        - 2. * _SQRT2_ * pvecmetric[ppw->index_mt_V_prime];

      /* photon velocity (theta_g = (3k/4)*F_1) */
      dy[pv->index_pt_theta_g] =
        k2 * (delta_g / 4. - s_l[2] * shear_g)
        - pvecthermo[pth->index_th_dkappa] * (theta_g + 4.0 / _SQRT6_ * P1);
    }

    /* photon shear (shear_g = F_2/2) */
    dy[pv->index_pt_shear_g] =
      4. / 15. * s_l[2] * theta_g - 3. / 10. * k * s_l[3] * y[pv->index_pt_shear_g + 1]
      - pvecthermo[pth->index_th_dkappa] * shear_g;

    /* photon l=3 */
    dy[pv->index_pt_l3_g] =
      k / 7. * (6. * s_l[3] * shear_g - 4. * s_l[4] * y[pv->index_pt_l3_g + 1])
      - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_l3_g];

    /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3,4) */
    for (l = 4; l < pv->l_max_g; l++)
      dy[pv->index_pt_delta_g + l] =
        k / (2. * l + 1.) * (l * s_l[l] * y[pv->index_pt_delta_g + l - 1]
                             - (l + 1.) * s_l[l + 1] * y[pv->index_pt_delta_g + l + 1])
        - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_delta_g + l];

    /* l=lmax */
    l                            = pv->l_max_g;
    dy[pv->index_pt_delta_g + l] =
      k * (s_l[l] * y[pv->index_pt_delta_g + l - 1]
           - (1. + l) * cotKgen * y[pv->index_pt_delta_g + l])
      - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_delta_g + l];

    /* photon polarization, l=0 (pol0_g = G_0)*/
    dy[pv->index_pt_pol0_g] =
      -k * y[pv->index_pt_pol0_g + 1]
      - pvecthermo[pth->index_th_dkappa] * (y[pv->index_pt_pol0_g] - _SQRT6_ * P1);

    /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3,4) */
    for (l = 1; l < pv->l_max_pol_g; l++)
      dy[pv->index_pt_pol0_g + l] =
        k / (2. * l + 1.) * (l * s_l[l] * y[pv->index_pt_pol0_g + l - 1]
                             - (l + 1.) * s_l[l + 1] * y[pv->index_pt_pol0_g + l + 1])
        - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_pol0_g + l];

    /* l=lmax */
    l                           = pv->l_max_pol_g;
    dy[pv->index_pt_pol0_g + l] =
      k * (s_l[l] * y[pv->index_pt_pol0_g + l - 1]
           - (l + 1.) * cotKgen * y[pv->index_pt_pol0_g + l])
      - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_pol0_g + l];

    /*
     *  }
     *  }
     */

    if (ppt->gauge == synchronous)
      /* Vector metric perturbation in synchronous gauge: */
      dy[pv->index_pt_hv_prime] = pvecmetric[ppw->index_mt_hv_prime_prime];

    else if (ppt->gauge == newtonian)
      /* Vector metric perturbation in Newtonian gauge: */
      dy[pv->index_pt_V] = pvecmetric[ppw->index_mt_V_prime];
  }


  /** - tensor modes: */
  if (_tensors_)
  {
    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
    {
      if (ppw->approx[ppw->index_ap_tca] == (int) tca_off)
      {
        /* short-cut notations for the tensor perturbations */
        delta_g = y[pv->index_pt_delta_g];
        theta_g = y[pv->index_pt_theta_g];
        shear_g = y[pv->index_pt_shear_g];


        /* (P^{(2)}) */
        P2 = -1.0 / _SQRT6_ * (
          1. / 10. * delta_g
          + 2. / 7. * shear_g
          + 3. / 70. * y[pv->index_pt_delta_g + 4]
          - 3. / 5. * y[pv->index_pt_pol0_g]
          + 6. / 7. * y[pv->index_pt_pol2_g]
          - 3. / 70. * y[pv->index_pt_pol0_g + 4]);

        /* above expression from paper, expression below matches old class but is not correct
         *  P2 = -1.0/_SQRT6_*(
         *  1./10.*delta_g
         +2./35.*shear_g
         +1./210.*y[pv->index_pt_delta_g+4]
         *  -3./5.*y[pv->index_pt_pol0_g]
         +6./35.*y[pv->index_pt_pol2_g]
         *  -1./210.*y[pv->index_pt_pol0_g+4]
         *  );
         */

        /* photon density (delta_g = F_0) */
        dy[pv->index_pt_delta_g] =
          -4. / 3. * theta_g
          - pvecthermo[pth->index_th_dkappa] * (delta_g + _SQRT6_ * P2)
          /*+y[pv->index_pt_gwdot]; */
          + _SQRT6_ * y[pv->index_pt_gwdot]; /*TBC */

        /* photon velocity (theta_g = (3k/4)*F_1) */
        dy[pv->index_pt_theta_g] =
          k2 * (delta_g / 4. - s_l[2] * shear_g)
          - pvecthermo[pth->index_th_dkappa] * theta_g;

        /* photon shear (shear_g = F_2/2) */
        dy[pv->index_pt_shear_g] =
          4. / 15. * s_l[2] * theta_g - 3. / 10. * k * s_l[3] * y[pv->index_pt_shear_g + 1]
          - pvecthermo[pth->index_th_dkappa] * shear_g;

        /* photon l=3 */
        dy[pv->index_pt_l3_g] =
          k / 7. * (6. * s_l[3] * shear_g - 4. * s_l[4] * y[pv->index_pt_l3_g + 1])
          - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_l3_g];

        /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3,4) */
        for (l = 4; l < pv->l_max_g; l++)
          dy[pv->index_pt_delta_g + l] =
            k / (2. * l + 1.) * (l * s_l[l] * y[pv->index_pt_delta_g + l - 1]
                                 - (l + 1.) * s_l[l + 1] * y[pv->index_pt_delta_g + l + 1])
            - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_delta_g + l];

        /* l=lmax */
        l                            = pv->l_max_g;
        dy[pv->index_pt_delta_g + l] =
          k * (s_l[l] * y[pv->index_pt_delta_g + l - 1]
               - (1. + l) * cotKgen * y[pv->index_pt_delta_g + l])
          - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_delta_g + l];

        /* photon polarization, l=0 (pol0_g = G_0)*/
        dy[pv->index_pt_pol0_g] =
          -k * y[pv->index_pt_pol0_g + 1]
          - pvecthermo[pth->index_th_dkappa] * (y[pv->index_pt_pol0_g] - _SQRT6_ * P2);

        /* additional momenta in Boltzmann hierarchy (beyond l=0,1,2,3,4) */
        for (l = 1; l < pv->l_max_pol_g; l++)
          dy[pv->index_pt_pol0_g + l] =
            k / (2. * l + 1.) * (l * s_l[l] * y[pv->index_pt_pol0_g + l - 1]
                                 - (l + 1.) * s_l[l + 1] * y[pv->index_pt_pol0_g + l + 1])
            - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_pol0_g + l];

        /* l=lmax */
        l                           = pv->l_max_pol_g;
        dy[pv->index_pt_pol0_g + l] =
          k * (s_l[l] * y[pv->index_pt_pol0_g + l - 1]
               - (l + 1.) * cotKgen * y[pv->index_pt_pol0_g + l])
          - pvecthermo[pth->index_th_dkappa] * y[pv->index_pt_pol0_g + l];
      }
    }

    if (ppt->evolve_tensor_ur == _TRUE_)
    {
      dy[pv->index_pt_delta_ur] = -4. / 3. * y[pv->index_pt_theta_ur] + _SQRT6_ * y[pv->index_pt_gwdot];

      dy[pv->index_pt_theta_ur] = k2 * (y[pv->index_pt_delta_ur] / 4. - s2_squared * y[pv->index_pt_shear_ur]);

      dy[pv->index_pt_shear_ur] = (4. / 15. * y[pv->index_pt_theta_ur]
                                   - 3. / 10. * k * s_l[3] / s_l[2] * y[pv->index_pt_shear_ur + 1]);

      l                      = 3;
      dy[pv->index_pt_l3_ur] = k / (2. * l + 1.) *
                               (l * 2. * s_l[l] * s_l[2] * y[pv->index_pt_shear_ur] - (l + 1.) * s_l[l + 1] * y[pv->index_pt_l3_ur + 1]);

      for (l = 4; l < pv->l_max_ur; l++)
      {
        dy[pv->index_pt_delta_ur + l] = k / (2. * l + 1) *
                                        (l * s_l[l] * y[pv->index_pt_delta_ur + l - 1] - (l + 1.) * s_l[l + 1] * y[pv->index_pt_delta_ur + l + 1]);
      }

      l                             = pv->l_max_ur;
      dy[pv->index_pt_delta_ur + l] =
        k * (s_l[l] * y[pv->index_pt_delta_ur + l - 1] - (1. + l) * cotKgen * y[pv->index_pt_delta_ur + l]);
    }

    /** - --> non-cold dark matter (ncdm): massive neutrinos, WDM, etc. */
    /*TBC: curvature in all ncdm */
    if (ppt->evolve_tensor_ncdm == _TRUE_)
    {
      idx = pv->index_pt_psi0_ncdm1;

      /** - ---> loop over species */

      for (n_ncdm = 0; n_ncdm < pv->N_ncdm; n_ncdm++)
      {
        /** - ----> loop over momentum */

        for (index_q = 0; index_q < pv->q_size_ncdm[n_ncdm]; index_q++)
        {
          /** - ----> define intermediate quantities */

          dlnf0_dlnq     = pba->dlnf0_dlnq_ncdm[n_ncdm][index_q];
          q              = pba->q_ncdm[n_ncdm][index_q];
          epsilon        = sqrt (q * q + a2 * pba->M_ncdm[n_ncdm] * pba->M_ncdm[n_ncdm]);
          qk_div_epsilon = k * q / epsilon;

          /** - ----> ncdm density for given momentum bin */

          dy[idx] = -qk_div_epsilon * y[idx + 1] - 0.25 * _SQRT6_ * y[pv->index_pt_gwdot] * dlnf0_dlnq;

          /** - ----> ncdm l>0 for given momentum bin */

          for (l = 1; l < pv->l_max_ncdm[n_ncdm]; l++)
          {
            dy[idx + l] = qk_div_epsilon / (2. * l + 1.0) * (l * s_l[l] * y[idx + (l - 1)] - (l + 1.) * s_l[l + 1] * y[idx + (l + 1)]);
          }

          /** - ----> ncdm lmax for given momentum bin (truncation as in Ma and Bertschinger)
           *   but with curvature taken into account a la arXiv:1305.3261 */

          dy[idx + l] = qk_div_epsilon * y[idx + l - 1] - (1. + l) * k * cotKgen * y[idx + l];

          /** - ----> jump to next momentum bin or species */

          idx += (pv->l_max_ncdm[n_ncdm] + 1);
        }
      }
    }

    /** - --> tensor metric perturbation h (gravitational waves) */
    dy[pv->index_pt_gw] = y[pv->index_pt_gwdot];

    /** - --> its time-derivative */
    dy[pv->index_pt_gwdot] = pvecmetric[ppw->index_mt_gw_prime_prime];
  }

  return _SUCCESS_;
}

int
perturb_tca_slip_and_shear (double   *y,
                            void     *parameters_and_workspace,
                            ErrorMsg error_message
                           )
{
  /** Summary: */

  /** - define local variables */
  int sgnK;
  double K = nc_hicosmo_curvature_K (((struct perturb_parameters_and_workspace *) parameters_and_workspace)->pba->cosmo, &sgnK);

  /* scale factor and other background quantities */
  double a, a_prime_over_a, a_primeprime_over_a, R;

  /* useful terms for tight-coupling approximation */
  double slip = 0.;
  double tau_c = 0., dtau_c = 0.;
  double theta_prime, shear_g_prime = 0., theta_prime_prime;
  double g0, g0_prime, g0_prime_prime;
  double F = 0., F_prime = 0., F_prime_prime = 0.;

  /* short-cut names for the fields of the input structure */
  struct perturb_parameters_and_workspace *pppaw;
  double k, k2;
  struct precision *ppr;
  const struct background *pba;
  struct thermo *pth;
  struct perturbs *ppt;
  struct perturb_workspace *ppw;
  double *pvecback;
  double *pvecthermo;
  double *pvecmetric;
  struct perturb_vector *pv;

  /* short-cut notations for the perturbations */
  double delta_g = 0., theta_g = 0., shear_g = 0.;
  double delta_b, theta_b;
  double Delta;
  double cb2;
  double metric_continuity = 0., metric_euler = 0., metric_shear = 0., metric_shear_prime = 0.;

  /* perturbed recombination */
  double delta_temp = 0.;

  /* for use with curvature */
  double s2_squared;

  /** - rename the fields of the input structure (just to avoid heavy notations) */

  pppaw = parameters_and_workspace;

  k  = pppaw->k;
  k2 = k * k;

  ppr = pppaw->ppr;
  pba = pppaw->pba;
  pth = pppaw->pth;
  ppt = pppaw->ppt;
  ppw = pppaw->ppw;

  pvecback   = ppw->pvecback;
  pvecthermo = ppw->pvecthermo;
  pvecmetric = ppw->pvecmetric;
  pv         = ppw->pv;

  /** - compute related background quantities */

  a                   = pvecback[pba->index_bg_a];
  a_prime_over_a      = pvecback[pba->index_bg_H] * a;
  a_primeprime_over_a = pvecback[pba->index_bg_H_prime] * a + 2. * a_prime_over_a * a_prime_over_a;
  /*z = pba->a_today-1.; */
  R          = 4. / 3. * pvecback[pba->index_bg_rho_g] / pvecback[pba->index_bg_rho_b];
  s2_squared = 1. - 3. * K / k2;

  /** - --> (a) define short-cut notations for the scalar perturbations */
  if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
  {
    delta_g = y[pv->index_pt_delta_g];
    theta_g = y[pv->index_pt_theta_g];
  }

  delta_b = y[pv->index_pt_delta_b];
  theta_b = y[pv->index_pt_theta_b];
  cb2     = pvecthermo[pth->index_th_cb2];

  /* perturbed recombination */
  if ((ppt->has_perturbed_recombination == _TRUE_) && (ppw->approx[ppw->index_ap_tca] == (int) tca_off))
    delta_temp = y[pv->index_pt_perturbed_recombination_delta_temp];

  /** - --> (b) define short-cut notations used only in tight-coupling approximation */
  tau_c  = 1. / pvecthermo[pth->index_th_dkappa];              /* inverse of opacity */
  dtau_c = -pvecthermo[pth->index_th_ddkappa] * tau_c * tau_c; /* its first derivative wrt conformal time */
  F      = tau_c / (1 + R);                                    /* F = tau_c/(1+R) */

  if (ppr->tight_coupling_approximation >= (int) second_order_CLASS)
  {
    F_prime = dtau_c / (1 + R) + tau_c * a_prime_over_a * R / (1 + R) / (1 + R); /*F' needed by second_order_CLASS and compromise_CLASS */

    if (ppr->tight_coupling_approximation == (int) second_order_CLASS)
      F_prime_prime = (-pvecthermo[pth->index_th_dddkappa] * tau_c * tau_c /* F'' needed by second_order_CLASS only */
                       + 2. * pvecthermo[pth->index_th_ddkappa] * pvecthermo[pth->index_th_ddkappa] * tau_c * tau_c * tau_c) / (1 + R)
                      + 2. * dtau_c * a_prime_over_a * R / (1 + R) / (1 + R)
                      + tau_c * ((a_primeprime_over_a - 2. * a_prime_over_a * a_prime_over_a) + 2. * a_prime_over_a * a_prime_over_a * R / (1 + R)) * R / (1 + R) / (1 + R);
  }

  /** - --> (c) compute metric-related quantities (depending on gauge; additional gauges can be coded below)
   *
   *       - Each continuity equation contains a term in (theta+metric_continuity) with
   *        metric_continuity = (h_prime/2) in synchronous gauge, (-3 phi_prime) in newtonian gauge
   *
   *       - Each Euler equation contains a source term metric_euler with
   *       metric_euler = 0 in synchronous gauge, (k2 psi) in newtonian gauge
   *
   *       - Each shear derivative equation contains a source term metric_shear equal to
   *        metric_shear = (h_prime+6eta_prime)/2 in synchronous gauge, 0 in newtonian gauge
   *
   *       - metric_shear_prime is the derivative of metric_shear
   *
   *       - In the ufa_class approximation, the leading-order source term is (h_prime/2) in synchronous gauge,
   *        (-3 (phi_prime+psi_prime)) in newtonian gauge: we approximate the later by (-6 phi_prime) */

  if (ppt->gauge == synchronous)
  {
    metric_continuity  = pvecmetric[ppw->index_mt_h_prime] / 2.;
    metric_euler       = 0.;
    metric_shear       = k2 * pvecmetric[ppw->index_mt_alpha];
    metric_shear_prime = k2 * pvecmetric[ppw->index_mt_alpha_prime];
  }

  if (ppt->gauge == newtonian)
  {
    metric_continuity  = -3. * pvecmetric[ppw->index_mt_phi_prime];
    metric_euler       = k2 * pvecmetric[ppw->index_mt_psi];
    metric_shear       = 0.;
    metric_shear_prime = 0.;
  }

  /** - --> (d) if some approximation schemes are turned on, enforce a few y[ ] values computed in perturb_einstein */

  /* free-streaming photon velocity */
  if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    theta_g = ppw->rsa_theta_g;


  /** - ---> like Ma & Bertschinger */
  if (ppr->tight_coupling_approximation == (int) first_order_MB)
    slip = 2. * R / (1. + R) * a_prime_over_a * (theta_b - theta_g)
           + F * (-a_primeprime_over_a * theta_b
                  + k2 * (-a_prime_over_a * delta_g / 2.
                          + cb2 * (-theta_b - metric_continuity)
                          - 4. / 3. * (-theta_g - metric_continuity) / 4.)
                  - a_prime_over_a * metric_euler);


  /** - ---> relax assumption dkappa~a\f$^{-2}\f$ (like in CAMB) */
  if ((ppr->tight_coupling_approximation == (int) first_order_CAMB) || (ppr->tight_coupling_approximation == (int) compromise_CLASS))
    slip = (dtau_c / tau_c - 2. * a_prime_over_a / (1. + R)) * (theta_b - theta_g)
           + F * (-a_primeprime_over_a * theta_b
                  + k2 * (-a_prime_over_a * delta_g / 2.
                          + cb2 * (-theta_b - metric_continuity)
                          - 4. / 3. * (-theta_g - metric_continuity) / 4.)
                  - a_prime_over_a * metric_euler);

  /** - ---> also relax assumption cb2~a\f$^{-1}\f$ */
  if ((ppr->tight_coupling_approximation == (int) first_order_CLASS) || (ppr->tight_coupling_approximation == (int) second_order_CLASS))
    slip = (dtau_c / tau_c - 2. * a_prime_over_a / (1. + R)) * (theta_b - theta_g)
           + F * (-a_primeprime_over_a * theta_b
                  + k2 * (-a_prime_over_a * delta_g / 2.
                          + pvecthermo[pth->index_th_dcb2] * delta_b
                          + cb2 * (-theta_b - metric_continuity)
                          - 4. / 3. * (-theta_g - metric_continuity) / 4.)
                  - a_prime_over_a * metric_euler);

  /** - ---> intermediate quantities for 2nd order tca: shear_g at first order in tight-coupling */
  shear_g = 16. / 45. * tau_c * (theta_g + metric_shear);

  /* (Ma & Bertschinger give (1/9)*(4/3) instead of (2/15)*(4/3)
   *  because they didn't include the contribution of G_gamma0
   *  and G_gamma2, which are of the same order as sigma_g. This
   *  was already consistently included in CAMB) */

  /** - ---> intermediate quantities for 2nd order tca: zero order for theta_b' = theta_g' */
  /** - ----> perturbed recombination has an impact **/
  theta_prime = (-a_prime_over_a * theta_b + k2 * (cb2 * (delta_b + delta_temp) + R / 4. * delta_g)) / (1. + R) + metric_euler;

  /** - ---> intermediate quantities for 2nd order tca: shear_g_prime at first order in tight-coupling */
  shear_g_prime = 16. / 45. * (tau_c * (theta_prime + metric_shear_prime) + dtau_c * (theta_g + metric_shear));

  /** - ---> 2nd order as in CRS*/
  if (ppr->tight_coupling_approximation == (int) second_order_CRS)
  {
    if (ppt->gauge == newtonian)
      class_stop (error_message,
                  "the second_order_CRS approach to tight-coupling is coded in synchronous gauge, not newtonian: change gauge or try another tight-coupling scheme");


    if (ppt->gauge == synchronous)
    {
      class_test (sgnK != 0,
                  ppt->error_message,
                  "the second_order_CRS approach to tight-coupling is coded in the flat case only: for non-flat try another tight-coupling scheme");

      /* infer Delta from h'' using Einstein equation */

      Delta = 2 * k2 * y[pv->index_pt_eta]
              - 2 * a_prime_over_a * pvecmetric[ppw->index_mt_h_prime]
              - pvecmetric[ppw->index_mt_h_prime_prime];

      /* monster expression for slip at second-order in tight-coupling */
      slip = (-2. / (1. + R) * a_prime_over_a - pvecthermo[pth->index_th_ddkappa] / pvecthermo[pth->index_th_dkappa]) * (theta_b - theta_g)
             + (-a_primeprime_over_a * theta_b
                - k2 * a_prime_over_a * (delta_g / 2. - 2. * shear_g)
                + k2 * (cb2 * (-theta_b - metric_continuity)
                        - 4. / 3. * (-theta_g - metric_continuity) / 4.
                        + shear_g_prime)
               ) / pvecthermo[pth->index_th_dkappa] / (1. + R)
             - 2. * R * (3. * a_prime_over_a * a_prime_over_a * cb2 + (1. + R) * (a_primeprime_over_a - a_prime_over_a * a_prime_over_a) - 3. * a_prime_over_a * a_prime_over_a)
             / (1. + R) / (1. + R) / (1. + R) * (theta_b - theta_g) / pvecthermo[pth->index_th_dkappa]
             + (
        a_primeprime_over_a * a_prime_over_a * ((2. - 3. * cb2) * R - 2.) * theta_b / (1. + R)
        + a_prime_over_a * k2 * (1. - 3. * cb2) * theta_b / 3. / (1. + R)
        /* perturbed recombination has an impact (next two lines) */
        + a_primeprime_over_a * k2 * cb2 * (delta_b + delta_temp) / (1. + R)
        + k2 * k2 * (3. * cb2 - 1.) * cb2 * (delta_b + delta_temp) / 3. / (1. + R)
        + k2 * k2 * R * (3. * cb2 - 1.) * delta_g / 12. / (1. + R)
        + a_primeprime_over_a * k2 * (2. + 3. * R) * delta_g / 4. / (1. + R)
        + a_prime_over_a * a_prime_over_a * k2 * ((2. - 3. * cb2) * R - 1.) * delta_g / 2. / (1. + R)
        + a_prime_over_a * k2 * cb2 * (1. + (3. * cb2 - 2.) * R) * (-theta_b - metric_continuity) / (1. + R)
        + a_prime_over_a * k2 * (2. + (5. - 3. * cb2) * R) * 4. / 3. * (-theta_g - metric_continuity) / 4. / (1. + R)
        + a_prime_over_a * (1. - 3. * cb2) * k2 * 2. * metric_shear / 3.
        + k2 * k2 * (3. * cb2 - 1.) * y[pv->index_pt_eta] / 3.
        + 2. * a_prime_over_a * k2 * (3. * cb2 - 1.) * pvecmetric[ppw->index_mt_eta_prime]
        + k2 * (1. - 3. * cb2) * Delta / 6.
               ) / pvecthermo[pth->index_th_dkappa] / pvecthermo[pth->index_th_dkappa] / (1. + R) / (1. + R)
             - (4. * a_primeprime_over_a * theta_b - 4. * k2 * cb2 * (-theta_b - metric_continuity) + 2. * a_prime_over_a * k2 * delta_g + k2 * 4. / 3. * (-theta_g - metric_continuity)) / 2. / (1. + R) / (1. + R) * pvecthermo[pth->index_th_ddkappa] / pvecthermo[pth->index_th_dkappa] / pvecthermo[pth->index_th_dkappa] / pvecthermo[pth->index_th_dkappa]
             + 4. * a_prime_over_a * R / (1. + R) / (1. + R) * pvecthermo[pth->index_th_ddkappa] / pvecthermo[pth->index_th_dkappa] / pvecthermo[pth->index_th_dkappa] * (theta_b - theta_g);

      /* second-order correction to shear */
      shear_g = (1. - 11. / 6. * dtau_c) * shear_g - 11. / 6. * tau_c * 16. / 45. * tau_c * (theta_prime + k2 * pvecmetric[ppw->index_mt_alpha_prime]);
    }
  }

  /** - ---> 2nd order like in CLASS paper */
  if (ppr->tight_coupling_approximation == (int) second_order_CLASS)
  {
    if (ppt->gauge == newtonian)
      class_stop (error_message,
                  "the second_order_CLASS approach to tight-coupling is coded in synchronous gauge, not newtonian: change gauge or try another tight-coupling scheme");


    if (ppt->gauge == synchronous)
    {
      /* zero order for theta_b'' = theta_g'' */
      theta_prime_prime = ((R - 1.) * a_prime_over_a * theta_prime - (a_primeprime_over_a - a_prime_over_a * a_prime_over_a) * theta_b
                           + k2 * (pvecthermo[pth->index_th_dcb2] * delta_b + cb2 * (-theta_b - metric_continuity) - a_prime_over_a * R / 4. * delta_g + R / 4. * 4. / 3. * (-theta_g - metric_continuity))) / (1. + R);

      /* zero-order quantities g0, g0', go'' */
      g0             = -a_prime_over_a * theta_b + k2 * (cb2 * delta_b - delta_g / 4.);
      g0_prime       = -a_prime_over_a * theta_prime - (a_primeprime_over_a - a_prime_over_a * a_prime_over_a) * theta_b + k2 * (pvecthermo[pth->index_th_dcb2] * delta_b + (1. / 3. - cb2) * (theta_b + 0.5 * pvecmetric[ppw->index_mt_h_prime]));
      g0_prime_prime = -a_prime_over_a * theta_prime_prime - 2. * (a_primeprime_over_a - a_prime_over_a * a_prime_over_a) * theta_prime
                       - (2. * a_prime_over_a * a_prime_over_a * a_prime_over_a - 3. * a_primeprime_over_a * a_prime_over_a) * theta_b
                       + k2 * (pvecthermo[pth->index_th_ddcb2] * delta_b - 2. * pvecthermo[pth->index_th_dcb2] * (theta_b + 0.5 * pvecmetric[ppw->index_mt_h_prime]) + (1. / 3. - cb2) * (theta_prime + 0.5 * pvecmetric[ppw->index_mt_h_prime_prime]));

      /* slip at second order */
      slip = (1. - 2 * a_prime_over_a * F) * slip + F * k2 * s2_squared * (2. * a_prime_over_a * shear_g + shear_g_prime)
             - F * (F_prime_prime * g0 + 2. * F_prime * g0_prime + F * g0_prime_prime);

      /* second-order correction to shear */
      shear_g = (1. - 11. / 6. * dtau_c) * shear_g - 11. / 6. * tau_c * 16. / 45. * tau_c * (theta_prime + metric_shear_prime);
    }
  }

  /** - ---> add only the most important 2nd order terms */
  if (ppr->tight_coupling_approximation == (int) compromise_CLASS)
  {
    /* slip at second order (only leading second-order terms) */
    slip = (1. - 2. * a_prime_over_a * F) * slip + F * k2 * (2. * a_prime_over_a * s2_squared * shear_g + s2_squared * shear_g_prime - (1. / 3. - cb2) * (F * theta_prime + 2. * F_prime * theta_b));

    /* second-order correction to shear */
    shear_g = (1. - 11. / 6. * dtau_c) * shear_g - 11. / 6. * tau_c * 16. / 45. * tau_c * (theta_prime + metric_shear_prime);
  }

  /** - ---> store tight-coupling values of photon shear and its derivative */

  ppw->tca_shear_g = shear_g;
  ppw->tca_slip    = slip;


  return _SUCCESS_;
}

int
perturb_rsa_delta_and_theta (
  struct precision         *ppr,
  const struct background  *pba,
  struct thermo            *pth,
  struct perturbs          *ppt,
  double                   k,
  double                   *y,
  double                   a_prime_over_a,
  double                   *pvecthermo,
  struct perturb_workspace *ppw
                            )
{
  /* - define local variables */

  double k2;

  k2 = k * k;

  /* formulas below TBC for curvaturema */

  /* newtonian gauge */
  if (ppt->gauge == newtonian)
  {
    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    {
      if (ppr->radiation_streaming_approximation == rsa_null)
      {
        ppw->rsa_delta_g = 0.;
        ppw->rsa_theta_g = 0.;
      }
      else
      {
        ppw->rsa_delta_g = -4. * y[ppw->pv->index_pt_phi];

        ppw->rsa_theta_g = 6. * ppw->pvecmetric[ppw->index_mt_phi_prime];
      }

      if (ppr->radiation_streaming_approximation == rsa_MD_with_reio)
      {
        ppw->rsa_delta_g +=
          -4. / k2 * ppw->pvecthermo[pth->index_th_dkappa] * y[ppw->pv->index_pt_theta_b];

        ppw->rsa_theta_g +=
          3. / k2 * (ppw->pvecthermo[pth->index_th_ddkappa] * y[ppw->pv->index_pt_theta_b]
                     + ppw->pvecthermo[pth->index_th_dkappa] *
                     (-a_prime_over_a * y[ppw->pv->index_pt_theta_b]
                      + ppw->pvecthermo[pth->index_th_cb2] * k2 * y[ppw->pv->index_pt_delta_b]
                      + k2 * y[ppw->pv->index_pt_phi]));
      }

      if (pba->has_ur == _TRUE_)
      {
        if (ppr->radiation_streaming_approximation == rsa_null)
        {
          ppw->rsa_delta_ur = 0.;
          ppw->rsa_theta_ur = 0.;
        }
        else
        {
          ppw->rsa_delta_ur = -4. * y[ppw->pv->index_pt_phi];
          ppw->rsa_theta_ur = 6. * ppw->pvecmetric[ppw->index_mt_phi_prime];
        }
      }
    }
  }

  /* synchronous gauge */
  if (ppt->gauge == synchronous)
  {
    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    {
      if (ppr->radiation_streaming_approximation == rsa_null)
      {
        ppw->rsa_delta_g = 0.;
        ppw->rsa_theta_g = 0.;
      }
      else
      {
        ppw->rsa_delta_g = 4. / k2 * (a_prime_over_a * ppw->pvecmetric[ppw->index_mt_h_prime]
                                      - k2 * y[ppw->pv->index_pt_eta]);
        ppw->rsa_theta_g = -0.5 * ppw->pvecmetric[ppw->index_mt_h_prime];
      }

      if (ppr->radiation_streaming_approximation == rsa_MD_with_reio)
      {
        ppw->rsa_delta_g +=
          -4. / k2 * ppw->pvecthermo[pth->index_th_dkappa] * (y[ppw->pv->index_pt_theta_b] + 0.5 * ppw->pvecmetric[ppw->index_mt_h_prime]);

        ppw->rsa_theta_g +=
          3. / k2 * (ppw->pvecthermo[pth->index_th_ddkappa] *
                     (y[ppw->pv->index_pt_theta_b]
                      + 0.5 * ppw->pvecmetric[ppw->index_mt_h_prime])
                     + ppw->pvecthermo[pth->index_th_dkappa] *
                     (-a_prime_over_a * y[ppw->pv->index_pt_theta_b]
                      + ppw->pvecthermo[pth->index_th_cb2] * k2 * y[ppw->pv->index_pt_delta_b]
                      - a_prime_over_a * ppw->pvecmetric[ppw->index_mt_h_prime]
                      + k2 * y[ppw->pv->index_pt_eta]));
      }

      if (pba->has_ur == _TRUE_)
      {
        if (ppr->radiation_streaming_approximation == rsa_null)
        {
          ppw->rsa_delta_ur = 0.;
          ppw->rsa_theta_ur = 0.;
        }
        else
        {
          ppw->rsa_delta_ur = 4. / k2 * (a_prime_over_a * ppw->pvecmetric[ppw->index_mt_h_prime]
                                         - k2 * y[ppw->pv->index_pt_eta]);
          ppw->rsa_theta_ur = -0.5 * ppw->pvecmetric[ppw->index_mt_h_prime];
        }
      }
    }
  }

  return _SUCCESS_;
}

