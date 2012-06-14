#include <glib.h>
#include <numcosmo/numcosmo.h>

int
main()
{
  ClModel *model;
  ClParams *cp;
  gint i;
  
  cl_cfg_init ();
  
  /* model = cl_model_xcdm_new (); */
  model = cl_model_linder_new ();
  cp = cl_model_get_params (model);
 
  cl_params_set (cp, CL_MODEL_DE_H0, 70.0);
  cl_params_set (cp, CL_MODEL_DE_OMEGA_C, 0.25);
  cl_params_set (cp, CL_MODEL_DE_OMEGA_X, 0.7);
  cl_params_set (cp, CL_MODEL_DE_T_GAMMA0, 1.0);
  cl_params_set (cp, CL_MODEL_DE_OMEGA_B, 0.05);
  cl_params_set (cp, CL_MODEL_DE_SPECINDEX, 1.0);
  cl_params_set (cp, CL_MODEL_DE_SIGMA8, 1.0);
  //cl_params_set (cp, CL_MODEL_DE_XCDM_W, -1.1);
  cl_params_set (cp, CL_MODEL_DE_LINDER_W0, -1.0);
  cl_params_set (cp, CL_MODEL_DE_LINDER_W1, 1.0);

  for (i = 0; i < 10; i++)
  {
    gdouble z = 1.0 / 9.0 * i;
    gdouble cd = CL_C_HUBBLE_RADIUS * cl_distance_comoving (cp, z, NULL);
    printf ("% 10.8f % 20.15g\n", z, cd);
  }

  cl_params_free (cp);
  cl_model_free (model);

  return 0;
}