/***************************************************************************
 *            test_ncm_sbessel_ode_solver_ivp.c
 *
 *  Mon Aug 24 2026
 *  Copyright 2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

static gdouble
_polynomial_forcing (gpointer user_data, gdouble x)
{
  const gdouble dx = x - 1.0;

  return x * x * (2.0 + dx * dx);
}

static void
test_ncm_sbessel_ode_solver_ivp_solve (void)
{
  NcmSBesselOdeSolverIVP *solver;
  gdouble u, du;

  solver = ncm_sbessel_ode_solver_ivp_new ();
  ncm_sbessel_ode_solver_ivp_set_reltol (solver, 1.0e-12);
  ncm_sbessel_ode_solver_ivp_set_abstol (solver, 1.0e-14);

  /* For ell = 0, f(x) = x^2 [2 + (x - 1)^2] gives u(x) = (x - 1)^2. */
  ncm_sbessel_ode_solver_ivp_solve (solver, &_polynomial_forcing, 1.0, 2.0, 0, NULL, &u, &du);
  ncm_assert_cmpdouble_e (u, ==, 1.0, 1.0e-10, 1.0e-12);
  ncm_assert_cmpdouble_e (du, ==, 2.0, 1.0e-10, 1.0e-12);

  ncm_sbessel_ode_solver_ivp_set_method (solver, NCM_SBESSEL_ODE_SOLVER_IVP_METHOD_FROBENIUS);
  g_assert_cmpint (ncm_sbessel_ode_solver_ivp_get_method (solver), ==, NCM_SBESSEL_ODE_SOLVER_IVP_METHOD_FROBENIUS);
  ncm_sbessel_ode_solver_ivp_solve (solver, &_polynomial_forcing, 1.0, 2.0, 0, NULL, &u, &du);
  ncm_assert_cmpdouble_e (u, ==, 1.0, 1.0e-10, 1.0e-12);
  ncm_assert_cmpdouble_e (du, ==, 2.0, 1.0e-10, 1.0e-12);

  ncm_sbessel_ode_solver_ivp_solve (solver, NULL, 1.0, 10.0, 5, NULL, &u, &du);
  ncm_assert_cmpdouble_e (u, ==, 0.0, 0.0, 1.0e-14);
  ncm_assert_cmpdouble_e (du, ==, 0.0, 0.0, 1.0e-14);

  ncm_sbessel_ode_solver_ivp_free (solver);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);

  g_test_add_func ("/ncm/specfunc/sbessel_ode_solver_ivp/solve", test_ncm_sbessel_ode_solver_ivp_solve);

  return g_test_run ();
}

