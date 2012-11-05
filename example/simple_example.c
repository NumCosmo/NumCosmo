#include <glib.h>
#include <numcosmo/numcosmo.h>

gint
main (gint argc, gchar *argv[])
{
  NcHICosmo *cosmo;
  NcDistance *dist;
  gint i;

  /**************************************************************************** 
   * Initialize the library objects, this must be called before 
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
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_H0, 70.0);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_C, 0.25);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X, 0.7);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_T_GAMMA0, 1.0);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_B, 0.05);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_SPECINDEX, 1.0);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_SIGMA8, 0.9);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_XCDM_W, -1.1);

  /**************************************************************************** 
   * Printing the parameters used.
   ****************************************************************************/
  printf ("# Model parameters:\n"); 
  ncm_model_params_log_all (NCM_MODEL (cosmo));

  /**************************************************************************** 
   * Printing some distances up to redshift 1.0.
   ****************************************************************************/ 
  for (i = 0; i < 10; i++)
  {
    gdouble z = 1.0 / 9.0 * i;
    gdouble cd = ncm_c_hubble_radius () * nc_distance_comoving (dist, cosmo, z);
    printf ("% 10.8f % 20.15g\n", z, cd);
  }

  /**************************************************************************** 
   * Freeing objects.
   ****************************************************************************/ 
  nc_distance_free (dist);
  ncm_model_free (NCM_MODEL (cosmo));

  return 0;
}
