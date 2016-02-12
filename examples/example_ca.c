#include <glib.h>
#include <numcosmo/numcosmo.h>

gint
main (gint argc, gchar *argv[])
{
  NcHIReion *reion;
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcWindow *wp;
  NcTransferFunc *tf;
  NcMatterVar *vp;
  NcGrowthFunc *gf;
  NcMultiplicityFunc *mulf;
  NcMassFunction *mf;
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
   * New cosmological distance objects optimizied to perform calculations
   * up to redshift 2.0.
   ****************************************************************************/  
  dist = nc_distance_new (2.0);
 
  /**************************************************************************** 
   * New windown function 'NcWindowTophat'
   ****************************************************************************/  
  wp = nc_window_new_from_name ("NcWindowTophat");

  /**************************************************************************** 
   * New transfer function 'NcTransferFuncEH' using the Einsenstein, Hu
   * fitting formula.
   ****************************************************************************/  
  tf = nc_transfer_func_new_from_name ("NcTransferFuncEH");

  /**************************************************************************** 
   * New matter variance object using FFT method for internal calculations and
   * the window and transfer functions defined above.
   ****************************************************************************/  
  vp = nc_matter_var_new (NC_MATTER_VAR_FFT, wp, tf);

  /**************************************************************************** 
   * New growth function
   ****************************************************************************/  
  gf = nc_growth_func_new ();

  /**************************************************************************** 
   * New multiplicity function 'NcMultiplicityFuncTinkerMean'
   ****************************************************************************/  
  mulf = nc_multiplicity_func_new_from_name ("NcMultiplicityFuncTinkerMean");

  /**************************************************************************** 
   * New mass function object using the objects defined above.
   ****************************************************************************/  
  mf = nc_mass_function_new (dist, vp, gf, mulf);

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
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_SPECINDEX,  1.0);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_SIGMA8,     0.9);
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
  nc_growth_func_prepare (gf, cosmo);

  for (i = 0; i < np; i++)
  {
    gdouble z = 1.0 / (np - 1.0) * i;
    gdouble gfz = nc_growth_func_eval (gf, cosmo, z);
    gdouble dgfz = nc_growth_func_eval_deriv (gf, cosmo, z);
    printf ("% 10.8f % 20.15g % 20.15g\n", z, gfz, dgfz);
  }
  printf ("\n\n");

  /**************************************************************************** 
   * Printing the transfer function and the matter power spectrum in the 
   * kh (in unities of h/Mpc) interval [1e-3, 1e3] 
   ****************************************************************************/
  nc_transfer_func_prepare (tf, reion, cosmo);

  for (i = 0; i < np; i++)
  {
    gdouble lnkh = log (1e-3) + log (1e6) / (np - 1.0) * i;
    gdouble tfkh = nc_transfer_func_eval (tf, cosmo, exp (lnkh));
    gdouble Pmkh = nc_transfer_func_matter_powerspectrum (tf, cosmo, exp (lnkh));
    printf ("% 10.8f % 20.15g % 20.15g\n", exp (lnkh), tfkh, Pmkh);
  }
  printf ("\n\n");

  /**************************************************************************** 
   * Printing the variance filtered with the tophat windown function using
   * scales R in the interval [5, 50] at redshift 0.3.
   * First calculates the growth function at z = 0.3 and then the spectrum
   * amplitude from the sigma8 parameter.
   ****************************************************************************/
  nc_matter_var_prepare (vp, reion, cosmo);
  {
    gdouble Dz = nc_growth_func_eval (gf, cosmo, 0.3);
    gdouble A = nc_matter_var_sigma8_sqrtvar0 (vp, cosmo);
    gdouble prefac = Dz * Dz * A * A;
    
    for (i = 0; i < np; i++)
    {
      gdouble lnR = log (5.0) + log (10.0) / (np - 1.0) * i;
      gdouble sigma2 = prefac * nc_matter_var_var0 (vp, cosmo, lnR);
      gdouble dsigma2_dlnR = nc_matter_var_dlnvar0_dlnR (vp, cosmo, lnR);
      printf ("% 10.8f % 20.15g % 20.15g\n", exp (lnR), sigma2, dsigma2_dlnR);
    }
    printf ("\n\n");
  }

  /**************************************************************************** 
   * Printing the mass function integrated in the mass interval [1e14, 1e16]
   * for the redhshifts in the interval [0, 2.0] and area 200 squared degree.
   ****************************************************************************/
  nc_mass_function_set_area_sd (mf, 200.0);
  nc_mass_function_set_eval_limits (mf, cosmo, log (1e14), log (1e16), 0.0, 2.0);
  nc_mass_function_prepare (mf, reion, cosmo);

  for (i = 0; i < np; i++)
  {
    gdouble z = 2.0 / (np - 1.0) * i;
    gdouble dndz = nc_mass_function_dn_dz (mf, cosmo, log(1e14), log(1e16), z, FALSE);
    printf ("% 10.8f % 20.15g\n", z, dndz);
  }
  printf ("\n\n");


  /**************************************************************************** 
   * Freeing objects.
   ****************************************************************************/ 
  nc_distance_free (dist);
  ncm_model_free (NCM_MODEL (cosmo));
  nc_hireion_free (reion);
  nc_window_free (wp);
  nc_transfer_func_free (tf);
  nc_matter_var_free (vp);
  nc_growth_func_free (gf);
  nc_multiplicity_func_free (mulf);
  nc_mass_function_free (mf);

  return 0;
}
