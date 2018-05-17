#include <glib.h>
#include <numcosmo/numcosmo.h>

gint
main (gint argc, gchar *argv[])
{
  ncm_cfg_init ();
  
  {
    NcmQMProp *p         = ncm_qm_prop_new_full (200, 5.0);
    //NcmQMPropGauss *psi0 = ncm_qm_prop_gauss_new (0.0, 1.0, 1.0, -0.5);
    NcmQMPropExp *psi0   = ncm_qm_prop_exp_new (ncm_qm_prop_spec_get_acs_a (p), 1.0, -1.5);
    gint i;

    //ncm_qm_prop_set_init_cond_gauss (p, psi0, 1.0e-4, 8.0);
    ncm_qm_prop_set_init_cond_exp (p, psi0, 0.0, 50.0);
    ncm_qm_prop_spec_prepare (p);
      
    for (i = 0; i < 1000; i++)
    {
      const gdouble t = 1.0e-4 * (i + 1);
      gint n = ncm_qm_prop_spec_nBohm (p);
      gint j;
      
      //ncm_qm_prop_evolve_spec (p, 0.01 * (i + 1));
      //ncm_qm_prop_evolve (p, 0.01 * (i + 1));
      ncm_qm_prop_spec_evol (p, t);
      printf ("# STEP %3d, t = %f:", i, t);
      for (j = 0; j < n; j++)
      {
        printf (" % 20.13g", ncm_qm_prop_spec_Bohm (p, j));
      }
      printf ("\n");
    }
    
    ncm_qm_prop_exp_free (psi0);
    ncm_qm_prop_free (p);
  }
}


