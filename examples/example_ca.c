#include <glib.h>
#include <numcosmo/numcosmo.h>

gint
main (gint argc, gchar *argv[])
{
  NcHIPrim  *prim;
  NcHIReion *reion;
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcTransferFunc *tf;
  NcMultiplicityFunc *mulf;
  NcHaloMassFunction *mf;
  NcPowspecML *psml;
  NcmPowspecFilter *psf;
  guint np = 50;
  gint i;

  /**************************************************************************** 
   * Initializing the library objects, this must be called before 
   * any other library function.
   ****************************************************************************/  
  ncm_cfg_init ();
  
  /**************************************************************************** 
   * New homogeneous and isotropic cosmological model NcHICosmoDEXcdm.
   ****************************************************************************/  
  cosmo = nc_hicosmo_new_from_name (NC_TYPE_HICOSMO, "NcHICosmoDEXcdm");

  /**************************************************************************** 
   * New homogeneous and isotropic reionization object.
   ****************************************************************************/  
  reion = NC_HIREION (nc_hireion_camb_new ());

  /**************************************************************************** 
   * New homogeneous and isotropic primordial object.
   ****************************************************************************/  
  prim  = NC_HIPRIM (nc_hiprim_power_law_new ());

  /**************************************************************************** 
   * Adding the submodels to the main cosmology model.
   ****************************************************************************/  
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (reion));
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (prim));

  /**************************************************************************** 
   * New cosmological distance objects optimizied to perform calculations
   * up to redshift 2.0.
   ****************************************************************************/  
  dist = nc_distance_new (2.0);
 
  /**************************************************************************** 
   * New transfer function 'NcTransferFuncEH' using the Einsenstein and Hu
   * fitting formula.
   ****************************************************************************/  
  tf = nc_transfer_func_new_from_name ("NcTransferFuncEH");

  /**************************************************************************** 
   * New linear matter power spectrum object based of the EH transfer function.
   ****************************************************************************/  
  psml = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  ncm_powspec_require_kmin (NCM_POWSPEC (psml), 1.0e-3);
  ncm_powspec_require_kmax (NCM_POWSPEC (psml), 1.0e3);

  /**************************************************************************** 
   * Apply a tophat filter to the psml object, set best output interval.
   ****************************************************************************/     
  psf = ncm_powspec_filter_new (NCM_POWSPEC (psml), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  ncm_powspec_filter_set_best_lnr0 (psf);

  /**************************************************************************** 
   * New multiplicity function 'NcMultiplicityFuncTinkerMean'
   ****************************************************************************/  
  mulf = nc_multiplicity_func_tinker_new_full (NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN, 200.0);
  
  /**************************************************************************** 
   * New mass function object using the objects defined above.
   ****************************************************************************/  
  mf = nc_halo_mass_function_new (dist, psf, mulf);

  /**************************************************************************** 
   * Setting values for the cosmological model, those not set stay in the
   * default values. Remeber to use the _orig_ version to set the original
   * parameters in case when a reparametrization is used.
   ****************************************************************************/ 
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_H0,        70.0);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_C,    0.25);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X,    0.7);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_T_GAMMA0,   1.0);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_B,    0.05);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_XCDM_W,    -1.1);

  /**************************************************************************** 
   * Printing the parameters used.
   ****************************************************************************/
  printf ("# Model parameters:\n#"); 
  ncm_model_params_log_all (NCM_MODEL (cosmo));

  /**************************************************************************** 
   * Printing the growth function and its derivative with respect to z 
   * up to redshift 1.
   ****************************************************************************/ 
  ncm_powspec_prepare (NCM_POWSPEC (psml), NCM_MODEL (cosmo));

  {
    NcGrowthFunc *gf = nc_powspec_ml_transfer_peek_gf (NC_POWSPEC_ML_TRANSFER (psml));
    for (i = 0; i < np; i++)
    {
      gdouble z = 1.0 / (np - 1.0) * i;
      gdouble gfz  = nc_growth_func_eval (gf, cosmo, z);
      gdouble dgfz = nc_growth_func_eval_deriv (gf, cosmo, z);
      printf ("% 10.8f % 21.15g % 21.15g\n", z, gfz, dgfz);
    }
    printf ("\n\n");
  }

  /**************************************************************************** 
   * Printing the matter power spectrum in the k (in unities of Mpc^-1) 
   * interval [1e-3, 1e3] 
   ****************************************************************************/

  for (i = 0; i < np; i++)
  {
    gdouble lnk = log (1e-3) + log (1e6) / (np - 1.0) * i;
    gdouble Pmk = ncm_powspec_eval (NCM_POWSPEC (psml), NCM_MODEL (cosmo), 0.0, exp (lnk));
    printf ("% 10.8f % 21.15g\n", exp (lnk), Pmk);
  }
  printf ("\n\n");

  /**************************************************************************** 
   * Printing the variance filtered with the tophat windown function using
   * scales R in the interval [5, 50] at redshift 0.3.
   * First calculates the growth function at z = 0.3 and then the spectrum
   * amplitude from the sigma8 parameter.
   ****************************************************************************/
  ncm_powspec_filter_prepare (psf, NCM_MODEL (cosmo));
  {
    for (i = 0; i < np; i++)
    {
      gdouble lnR = log (5.0) + log (10.0) / (np - 1.0) * i;
      gdouble sigma2       = ncm_powspec_filter_eval_var_lnr (psf, 0.0, lnR);
      gdouble dsigma2_dlnR = ncm_powspec_filter_eval_dvar_dlnr (psf, 0.0, lnR);
      printf ("% 10.8f % 21.15g % 21.15g\n", exp (lnR), sigma2, dsigma2_dlnR);
    }
    printf ("\n\n");
  }

  /**************************************************************************** 
   * Printing the mass function integrated in the mass interval [1e14, 1e16]
   * for the redhshifts in the interval [0, 2.0] and area 200 squared degree.
   ****************************************************************************/
  nc_halo_mass_function_set_area_sd (mf, 200.0);
  nc_halo_mass_function_set_eval_limits (mf, cosmo, log (1e14), log (1e16), 0.0, 2.0);
  nc_halo_mass_function_prepare (mf, cosmo);

  for (i = 0; i < np; i++)
  {
    gdouble z = 2.0 / (np - 1.0) * i;
    gdouble dndz = nc_halo_mass_function_dn_dz (mf, cosmo, log (1.0e14), log (1.0e16), z, FALSE);
    printf ("% 10.8f % 21.15g\n", z, dndz);
  }
  printf ("\n\n");


  /**************************************************************************** 
   * Freeing objects.
   ****************************************************************************/ 
  nc_distance_free (dist);
  ncm_model_free (NCM_MODEL (cosmo));
  nc_hireion_free (reion);
  nc_transfer_func_free (tf);
  nc_powspec_ml_free (psml);
  ncm_powspec_filter_free (psf);
  nc_multiplicity_func_free (mulf);
  nc_halo_mass_function_free (mf);

  return 0;
}
