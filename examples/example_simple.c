#include <glib.h>
#include <numcosmo/numcosmo.h>

gint
main (gint argc, gchar *argv[])
{
  NcHICosmo *cosmo;
  NcDistance *dist;

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
   * New cosmological distance objects optimizied to perform calculations
   * up to redshift 2.0.
   ****************************************************************************/  
  dist = nc_distance_new (2.0);
 
  /**************************************************************************** 
   * Setting values for the cosmological model, those not set stay in the
   * default values. Remeber to use the _orig_ version to set the original
   * parameters in case when a reparametrization is used.
   ****************************************************************************/ 
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_H0,       70.00);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_C,   0.25);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X,   0.70);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.72);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_B,   0.05);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_XCDM_W,   -1.10);

  /**************************************************************************** 
   * Printing the parameters used.
   ****************************************************************************/
  printf ("# Model parameters:\n"); 
  ncm_model_params_log_all (NCM_MODEL (cosmo));

  /**************************************************************************** 
   * Printing some distances up to redshift 1.0.
   ****************************************************************************/ 
  {
    const guint N        = 20;
    const gdouble RH_Mpc = nc_hicosmo_RH_Mpc (cosmo); /* Hubble radius in Mpc */
    gint i;
    
    for (i = 0; i < N; i++)
    {
      const gdouble z  = 1.0 / (N - 1.0) * i;
      const gdouble Dc = nc_distance_comoving (dist, cosmo, z); /* Comoving distance in Hubble radius */
      const gdouble dc = RH_Mpc * Dc;                           /* Comoving distance in Mpc           */
      printf ("% 10.8f % 22.15g [c/H0] % 22.15g [Mpc]\n", z, Dc, dc);
    }
  }

  /**************************************************************************** 
   * Freeing objects.
   ****************************************************************************/ 
  nc_distance_free (dist);
  nc_hicosmo_free (cosmo);

  return 0;
}
