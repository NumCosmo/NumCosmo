#include <glib.h>
#include <numcosmo/numcosmo.h>

gint
main (gint argc, gchar *argv[])
{
  ncm_cfg_init ();
  
  {
    NcmQMProp *p         = ncm_qm_prop_new ();
    NcmQMPropGauss *psi0 = ncm_qm_prop_gauss_new (0.0, 1.0, 1.0, 1.0);
    gint i;
    
    g_object_set (p, 
      "nknots", 1500, 
      "lambda", 0.0, 
      NULL);
      
    ncm_qm_prop_set_init_cond_gauss (p, psi0, 0.0, 30.0);
    
    for (i = 0; i < 10; i++)
    {
      ncm_qm_prop_evolve (p, 0.01 * (i + 1));
      printf ("# STEP %d\n", i);
    }
  }
}


