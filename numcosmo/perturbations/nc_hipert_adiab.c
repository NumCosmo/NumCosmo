/***************************************************************************
 *            nc_hipert_adiab.c
 *
 *  Tue June 03 17:20:42 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_adiab.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
 *
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION:nc_hipert_adiab
 * @title: Adiabatic Perturbation Object
 * @short_description: Perturbation object for adiabatic mode only 
 *
 * This object provides the computation of the adiabatic mode for the cosmological 
 * perturbations. It solves the equation of motion for the gauge invariant variable
 * (see <link linkend="XVitenti2013">Vitenti (2013)</link> for notation and details)
 * $$\zeta \equiv \Psi - \frac{2\bar{K}}{\kappa(\bar{\rho} + \bar{p})} + H\mathcal{V}.$$
 * Its conjugated momentum is give by
 * \begin{split}
 * P_\zeta &= \frac{2\bar{D}^2_\bar{K}\Psi}{x^3H},
 * \end{split}
 * 
 * The equations of motion in their first order form are
 * \begin{align}
 * \zeta^\prime &= \frac{P_\zeta}{m_\zeta}, \\\\
 * P_\zeta^\prime &= -m_\zeta\mu_\zeta^2\zeta.
 * \end{align}
 * The mass $m_\zeta$ and the frequency $\mu_\zeta$ are defined by
 * \begin{align}
 * m_\zeta     &= \frac{3\Delta_\bar{K}(\bar{\rho} + \bar{p})}{\rho_\text{crit0} N x^3 c_s^2 E^2}, \\\\
 * \mu_\zeta^2 &= x^2N^2c_s^2k^2,
 * \end{align}
 * where $\bar{\rho} + \bar{p}$ is the background total energy density plus pressure defined by nc_hicosmo_rhopp(),
 * $E^2 = H^2/H_0^2$ is the adimensional Hubble function squared (nc_hicosmo_E2()), $c_s^2$ the speed of sound (nc_hicosmo_cs2()),
 * $N$ is the lapse function that in this case (using $\alpha$ as time variable) is $N \equiv \vert{}E\vert^{-1}$, $\rho_\text{crit0}$
 * is the critical density today defined by $\rho_\text{crit0} \equiv 3H_0^2/\kappa$ and $$\Delta_\bar{K} \equiv \frac{k^2}{k^2 + \Omega_k}.$$
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic_notaknot.h"
#include "perturbations/nc_hipert_iadiab.h"
#include "perturbations/nc_hipert_adiab.h"

#include "math/cvode_util.h"

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <nvector/nvector_serial.h> 
#include <gsl/gsl_roots.h>

G_DEFINE_TYPE (NcHIPertAdiab, nc_hipert_adiab, NC_TYPE_HIPERT);

enum {
  PROP_0,
  PROP_SIZE,
};

static gdouble _nc_hipert_adiab_wkb_phase (gdouble y, gdouble x, gpointer userdata);

typedef struct _NcHIPertAdiabArg
{
  NcHICosmo *cosmo;
  NcHIPertAdiab *pa;
} NcHIPertAdiabArg;

static void
nc_hipert_adiab_init (NcHIPertAdiab *pa)
{
  NcmSpline *s     = ncm_spline_cubic_notaknot_new ();
  pa->wkb_phase    = ncm_ode_spline_new (s, &_nc_hipert_adiab_wkb_phase);

  ncm_ode_spline_set_abstol (pa->wkb_phase, 0.0);
  ncm_spline_free (s);

  pa->mzetanuA         = N_VNew_Serial (3);
  pa->cvode_phase      = CVodeCreate (CV_BDF, CV_NEWTON);
  pa->cvode_phase_init = FALSE;

  pa->wkb_prep = FALSE;
  pa->ode_prep = FALSE;

  pa->lnphase  = ncm_spline_cubic_notaknot_new ();
  pa->lnFzeta  = ncm_spline_cubic_notaknot_new ();
  pa->dlnFzeta = ncm_spline_cubic_notaknot_new ();
}

static void
nc_hipert_adiab_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /* NcHIPertAdiab *pa = NC_HIPERT_ADIAB (object); */
  g_return_if_fail (NC_IS_HIPERT_ADIAB (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_adiab_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /* NcHIPertAdiab *pa = NC_HIPERT_ADIAB (object); */
  g_return_if_fail (NC_IS_HIPERT_ADIAB (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_adiab_dispose (GObject *object)
{
  NcHIPertAdiab *pa = NC_HIPERT_ADIAB (object);

  ncm_ode_spline_clear (&pa->wkb_phase);

  ncm_spline_clear (&pa->lnphase);
  ncm_spline_clear (&pa->lnFzeta);
  ncm_spline_clear (&pa->dlnFzeta);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_adiab_parent_class)->dispose (object);
}


static void
nc_hipert_adiab_finalize (GObject *object)
{
  NcHIPertAdiab *pa = NC_HIPERT_ADIAB (object);

  if (pa->cvode_phase != NULL)
  {
    CVodeFree (&pa->cvode_phase);
    pa->cvode_phase = NULL;
  }
  
  if (pa->mzetanuA != NULL)
  {
    N_VDestroy (pa->mzetanuA);
    pa->mzetanuA = NULL;
  }

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_adiab_parent_class)->finalize (object);
}

static void
nc_hipert_adiab_class_init (NcHIPertAdiabClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->dispose      = &nc_hipert_adiab_dispose;
  object_class->finalize     = &nc_hipert_adiab_finalize;
  object_class->set_property = &nc_hipert_adiab_set_property;
  object_class->get_property = &nc_hipert_adiab_get_property;
  
}

/**
 * nc_hipert_adiab_new:
 * 
 * Creates a new #NcHIPertAdiab object.
 * 
 * Returns: (transfer full): a new #NcHIPertAdiab.
 */
NcHIPertAdiab *
nc_hipert_adiab_new (void)
{
  NcHIPertAdiab *pa = g_object_new (NC_TYPE_HIPERT_ADIAB, 
                                        NULL);

  return pa;
}

/**
 * nc_hipert_adiab_ref:
 * @pa: a #NcHIPertAdiab.
 * 
 * Increases the reference count of @pa.
 * 
 * Returns: (transfer full): @pa. 
 */
NcHIPertAdiab *
nc_hipert_adiab_ref (NcHIPertAdiab *pa)
{
  return g_object_ref (pa);
}

/**
 * nc_hipert_adiab_free:
 * @pa: a #NcHIPertAdiab.
 * 
 * Decreases the reference count of @pa.
 * 
 */
void 
nc_hipert_adiab_free (NcHIPertAdiab *pa)
{
  g_object_unref (pa);
}

/**
 * nc_hipert_adiab_clear:
 * @pa: a #NcHIPertAdiab.
 * 
 * Decreases the reference count of *@pa and sets *@pa to NULL.
 * 
 */
void 
nc_hipert_adiab_clear (NcHIPertAdiab **pa)
{
  g_clear_object (pa);
}

static gdouble 
_nc_hipert_adiab_wkb_phase (gdouble y, gdouble x, gpointer userdata)
{
  NcHIPertAdiabArg *arg = (NcHIPertAdiabArg *) userdata;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (arg->cosmo), x, NC_HIPERT (arg->pa)->k);
  return sqrt (nuA2);
}

/**
 * nc_hipert_adiab_prepare_wkb:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha_i: initial log-redshift time.
 * @alpha_f: final log-redshift time.
 * 
 * Prepare the object for WKB calculations using the cosmology @cosmo.
 * 
 */
void 
nc_hipert_adiab_prepare_wkb (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha_i, gdouble alpha_f)
{
  NcHIPert *pert = NC_HIPERT (pa);
  if (!pert->prepared)
  {
    pa->wkb_prep   = FALSE;
    pa->ode_prep   = FALSE;
    pert->prepared = TRUE;
  }  
  
  if (!pa->wkb_prep)
  {
    NcHIPertAdiabArg arg;

    arg.cosmo = cosmo;
    arg.pa    = pa;

    ncm_ode_spline_set_interval (pa->wkb_phase, M_PI * 0.25, alpha_i, alpha_f);

    ncm_ode_spline_prepare (pa->wkb_phase, &arg);

    pa->wkb_prep = TRUE;
  }
}

static gint
_nc_hipert_adiab_phase_f (realtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertAdiabArg *arg   = (NcHIPertAdiabArg *) f_data;
  NcHIPertIAdiabEOM *eom  = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (arg->cosmo), alpha, NC_HIPERT (arg->pa)->k);
  const gdouble dlnmzeta  = nc_hipert_iadiab_dlnmzeta (NC_HIPERT_IADIAB (arg->cosmo), alpha, NC_HIPERT (arg->pa)->k);
  const gdouble lnFzeta   = NV_Ith_S (y, 0);
  const gdouble dlnFzeta  = NV_Ith_S (y, 1);
  const gdouble lnphase   = NV_Ith_S (y, 2);
  const gdouble dlnFzeta2 = dlnFzeta * dlnFzeta;
  const gdouble mzeta2    = eom->m * eom->m; 
  
  NV_Ith_S (ydot, 0) = dlnFzeta;
  NV_Ith_S (ydot, 1) = 0.5 * dlnFzeta2 - dlnFzeta * dlnmzeta + 2.0 * (eom->nu2 - exp (2.0 * lnFzeta) / mzeta2);
  NV_Ith_S (ydot, 2) = exp (lnFzeta - lnphase) / eom->m;


/*printf ("% 20.15e % 20.15e % 20.15e %e\n", nc_hicosmo_x_alpha (arg->cosmo, alpha), eom->nu2, exp (2.0 * lnFzeta) / mzeta2, fabs ((exp (2.0 * lnFzeta) / mzeta2 - eom->nu2)/eom->nu2) );*/
  
  return 0;
}

static gint
_nc_hipert_adiab_phase_J (_NCM_SUNDIALS_INT_TYPE N, realtype alpha, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertAdiabArg *arg  = (NcHIPertAdiabArg *) jac_data;
  NcHIPertIAdiabEOM *eom = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (arg->cosmo), alpha, NC_HIPERT (arg->pa)->k);
  const gdouble dlnmzeta  = nc_hipert_iadiab_dlnmzeta (NC_HIPERT_IADIAB (arg->cosmo), alpha, NC_HIPERT (arg->pa)->k);
  const gdouble lnFzeta  = NV_Ith_S (y, 0);
  const gdouble dlnFzeta = NV_Ith_S (y, 1);
  const gdouble lnphase  = NV_Ith_S (y, 2);
  const gdouble mzeta2   = eom->m * eom->m; 

  DENSE_ELEM (J, 0, 1) = 1.0;

  DENSE_ELEM (J, 1, 0) = - 4.0 * exp (2.0 * lnFzeta) / mzeta2;
  DENSE_ELEM (J, 1, 1) = dlnFzeta - dlnmzeta;

  DENSE_ELEM (J, 2, 0) = exp (lnFzeta - lnphase) / eom->m;
  DENSE_ELEM (J, 2, 2) = -exp (lnFzeta - lnphase) / eom->m;
  
  return 0;
}

/**
 * nc_hipert_adiab_prepare_ode:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha_i: initial log-redshift time.
 * @alpha_f: final log-redshift time.
 * 
 * Prepare the object for WKB calculations using the cosmology @cosmo. Instead
 * of using the wkb approximation as in @nc_hipert_adiab_prepare_wkb it solves
 * the non-linear equation of motion for $\nu_A$. 
 * 
 */
void 
nc_hipert_adiab_prepare_ode (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha_i, gdouble alpha_f)
{
  NcHIPert *pert = NC_HIPERT (pa);
  if (!pert->prepared)
  {
    pa->wkb_prep   = FALSE;
    pa->ode_prep   = FALSE;
    pert->prepared = TRUE;
  }  

  if (!pa->ode_prep)
  {
    gint flag;
    const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (cosmo), alpha_i, pert->k);
    const gdouble nuA = sqrt (nuA2);
    const gdouble dmzetanuA_nuA = nc_hipert_iadiab_dmzetanuA_nuA (NC_HIPERT_IADIAB (cosmo), alpha_i, pert->k);
    NcHIPertIAdiabEOM *eom = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (cosmo), alpha_i, pert->k);
    const gdouble dmzetanuA = dmzetanuA_nuA * nuA;
    gdouble alpha = alpha_i;
    GArray *alpha_a;
    GArray *lnphase;
    GArray *lnFzeta;
    GArray *dlnFzeta;
    
    NcmVector *v = NULL;

    if ((v = ncm_spline_get_xv (pa->lnFzeta)) != NULL)
    {
      alpha_a = ncm_vector_get_array (v);
      ncm_vector_free (v);

      v = ncm_spline_get_yv (pa->lnFzeta);
      lnFzeta = ncm_vector_get_array (v);
      ncm_vector_free (v);

      v = ncm_spline_get_yv (pa->dlnFzeta);
      dlnFzeta = ncm_vector_get_array (v);
      ncm_vector_free (v);

      v = ncm_spline_get_yv (pa->lnphase);
      lnphase = ncm_vector_get_array (v);
      ncm_vector_clear (&v);

      g_array_set_size (alpha_a, 0);
      g_array_set_size (lnphase, 0);
      g_array_set_size (lnFzeta, 0);
      g_array_set_size (dlnFzeta, 0);
    }
    else
    {
      alpha_a  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
      lnphase  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
      lnFzeta  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
      dlnFzeta = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    }
 
    NV_Ith_S (pa->mzetanuA, 0) = log (nuA * eom->m);
    NV_Ith_S (pa->mzetanuA, 1) = dmzetanuA / (nuA * eom->m);
    NV_Ith_S (pa->mzetanuA, 2) = log (M_PI * 0.25);    
    
    if (!pa->cvode_phase_init)
    {
      flag = CVodeInit (pa->cvode_phase, &_nc_hipert_adiab_phase_f, alpha_i, pa->mzetanuA);
      CVODE_CHECK (&flag, "CVodeInit", 1, );
      pa->cvode_phase_init = TRUE;
    }
    else
    {
      flag = CVodeReInit (pa->cvode_phase, alpha_i, pa->mzetanuA);
      CVODE_CHECK (&flag, "CVodeReInit", 1, );
    }

    flag = CVodeSStolerances (pa->cvode_phase, pert->reltol, pert->abstol);
    CVODE_CHECK(&flag, "CVodeSStolerances", 1,);

    flag = CVodeSetMaxNumSteps (pa->cvode_phase, 1000000);
    CVODE_CHECK(&flag, "CVodeSetMaxNumSteps", 1, );

    flag = CVDense (pa->cvode_phase, 3);
    CVODE_CHECK(&flag, "CVDense", 1, );

    flag = CVDlsSetDenseJacFn (pa->cvode_phase, &_nc_hipert_adiab_phase_J);
    CVODE_CHECK(&flag, "CVDlsSetDenseJacFn", 1, );  

    {
      gdouble last = alpha_i;
      NcHIPertAdiabArg arg;
      arg.cosmo = cosmo;
      arg.pa = pa;
      
      flag = CVodeSetUserData (pa->cvode_phase, &arg);
      CVODE_CHECK (&flag, "CVodeSetFdata", 1, );

      while (alpha < alpha_f)
      {
        flag = CVode (pa->cvode_phase, alpha_f, pa->mzetanuA, &alpha, CV_ONE_STEP);
        CVODE_CHECK (&flag, "CVode[nc_hipert_adiab_prepare_ode]", 1, );

        if (fabs (2.0 * (alpha - last) / (alpha + last)) > 1.0e-6)
        {
          g_array_append_val (alpha_a, alpha);
          g_array_append_val (lnFzeta, NV_Ith_S (pa->mzetanuA, 0));
          g_array_append_val (dlnFzeta, NV_Ith_S (pa->mzetanuA, 1));
          g_array_append_val (lnphase, NV_Ith_S (pa->mzetanuA, 2));
          last = alpha;
        }
      }
    }

    ncm_spline_set_array (pa->lnphase, alpha_a, lnphase, TRUE);
    ncm_spline_set_array (pa->lnFzeta, alpha_a, lnFzeta, TRUE);
    ncm_spline_set_array (pa->dlnFzeta, alpha_a, dlnFzeta, TRUE);

    g_array_unref (alpha_a);
    g_array_unref (lnphase);
    g_array_unref (lnFzeta);
    g_array_unref (dlnFzeta);
    pa->ode_prep = TRUE;
  }
}

/**
 * nc_hipert_adiab_wkb_zeta:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @Re_zeta: (out caller-allocates): Real part of the wkb solution.
 * @Im_zeta: (out caller-allocates): Imaginary part of the wkb solution.
 * 
 * Computes the WKB solution $\zeta_\text{WKB}$ for the mode $k$ at the time $\alpha$. 
 * 
 */
void
nc_hipert_adiab_wkb_zeta (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta)
{
  NcHIPert *pert = NC_HIPERT (pa);
  NcHIPertIAdiabEOM *eom = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  complex double zeta;
  gdouble int_nuA;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  const gdouble nuA = sqrt (nuA2); 

  g_assert (pa->wkb_prep);

  int_nuA = ncm_spline_eval (pa->wkb_phase->s, alpha); 
  
  zeta = cexp (-I * int_nuA) / sqrt (2.0 * eom->m * nuA);

  *Re_zeta = creal (zeta);
  *Im_zeta = cimag (zeta);
}

/**
 * nc_hipert_adiab_wkb_zeta_Pzeta:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @Re_zeta: (out caller-allocates): Real part of the wkb solution.
 * @Im_zeta: (out caller-allocates): Imaginary part of the wkb solution.
 * @Re_Pzeta: (out caller-allocates): Real part of the wkb solution momentum.
 * @Im_Pzeta: (out caller-allocates): Imaginary part of the wkb solution momentum.
 * 
 * Computes the WKB solution $\zeta_\text{WKB}$ and its momentum for the mode $k$ at the time $\alpha$. 
 * 
 */
void
nc_hipert_adiab_wkb_zeta_Pzeta (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta, gdouble *Re_Pzeta, gdouble *Im_Pzeta)
{
  NcHIPert *pert = NC_HIPERT (pa);
  NcHIPertIAdiabEOM *eom = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  complex double zeta, Pzeta;
  gdouble int_nuA;
  gdouble dmzetanuA_nuA;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  const gdouble nuA = sqrt (nuA2); 

  g_assert (pa->wkb_prep);

  int_nuA = ncm_spline_eval (pa->wkb_phase->s, alpha);

  dmzetanuA_nuA = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);

  zeta = cexp (-I * int_nuA) / sqrt (2.0 * eom->m * nuA);

  Pzeta = -I * cexp (-I * int_nuA) * sqrt (0.5 * eom->m * nuA) - 0.5 * dmzetanuA_nuA * zeta;

  *Re_zeta = creal (zeta);
  *Im_zeta = cimag (zeta);

  *Re_Pzeta = creal (Pzeta);
  *Im_Pzeta = cimag (Pzeta);  
}

/**
 * nc_hipert_adiab_ode_zeta:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @Re_zeta: (out caller-allocates): Real part of the ode solution.
 * @Im_zeta: (out caller-allocates): Imaginary part of the ode solution.
 * 
 * Computes the solution $\zeta_\text{ode}$ for the mode $k$ at the time $\alpha$. 
 * 
 */
void
nc_hipert_adiab_ode_zeta (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta)
{
  complex double zeta;
  gdouble phase;
  gdouble lnFzeta;
  const gdouble one_sqrt2 = 1.0 / sqrt (2.0);

  g_assert (pa->ode_prep);

  phase = exp (ncm_spline_eval (pa->lnphase, alpha));
  lnFzeta = ncm_spline_eval (pa->lnFzeta, alpha);
  
  zeta = cexp (-I * phase - lnFzeta * 0.5) * one_sqrt2;

  *Re_zeta = creal (zeta);
  *Im_zeta = cimag (zeta);
}

/**
 * nc_hipert_adiab_ode_zeta_Pzeta:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @Re_zeta: (out caller-allocates): Real part of the ode solution.
 * @Im_zeta: (out caller-allocates): Imaginary part of the ode solution.
 * @Re_Pzeta: (out caller-allocates): Real part of the ode solution momentum.
 * @Im_Pzeta: (out caller-allocates): Imaginary part of the ode solution momentum.
 * 
 * Computes the ode solution $\zeta_\text{ode}$ and its momentum for the mode $k$ at the time $\alpha$. 
 * 
 */
void
nc_hipert_adiab_ode_zeta_Pzeta (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_zeta, gdouble *Im_zeta, gdouble *Re_Pzeta, gdouble *Im_Pzeta)
{
  NcHIPert *pert = NC_HIPERT (pa);
  NcHIPertIAdiabEOM *eom = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  complex double zeta, Pzeta;
  gdouble phase;
  gdouble lnFzeta;
  gdouble dlnFzeta;
  const gdouble one_sqrt2 = 1.0 / sqrt (2.0);

  g_assert (pa->ode_prep);

  phase    = exp (ncm_spline_eval (pa->lnphase, alpha));
  lnFzeta  = ncm_spline_eval (pa->lnFzeta, alpha);
  dlnFzeta = ncm_spline_eval (pa->dlnFzeta, alpha);

  zeta = cexp (-I * phase - lnFzeta * 0.5) * one_sqrt2;

  Pzeta = -I * cexp (-I * phase + lnFzeta * 0.5) * one_sqrt2 - 0.5 * eom->m * dlnFzeta * zeta;

  *Re_zeta = creal (zeta);
  *Im_zeta = cimag (zeta);

  *Re_Pzeta = creal (Pzeta);
  *Im_Pzeta = cimag (Pzeta);  
}

static gdouble 
_nc_hipert_adiab_wkb_prec (gdouble alpha, gpointer userdata)
{
  NcHIPertAdiabArg *arg = (NcHIPertAdiabArg *) userdata;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (arg->cosmo), alpha, NC_HIPERT (arg->pa)->k);
  const gdouble test = log (fabs (nuA2 * NC_HIPERT (arg->pa)->reltol / (alpha * alpha)));
  return test;
}

static gdouble 
_nc_hipert_adiab_wkb_nuA2 (gdouble alpha, gpointer userdata)
{
  NcHIPertAdiabArg *arg = (NcHIPertAdiabArg *) userdata;
  const gdouble nuA2 = nc_hipert_iadiab_nuA2 (NC_HIPERT_IADIAB (arg->cosmo), alpha, NC_HIPERT (arg->pa)->k);
  NcHIPertIAdiabEOM *eom  = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (arg->cosmo), alpha, NC_HIPERT (arg->pa)->k);
  return nuA2 / eom->nu2;
}

static gint
_nc_hipert_adiab_wkb_nu2_root (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha0, gdouble *alpha, gdouble (*f) (gdouble, gpointer))
{
  NcHIPert *pert = NC_HIPERT (pa);
  gint status;
  gint iter = 0, max_iter = 1000000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  gdouble prec = pert->reltol, alpha1 = *alpha;
  NcHIPertAdiabArg arg;

  arg.cosmo = cosmo;
  arg.pa    = pa;

  F.function = f;
  F.params   = &arg;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, alpha0, alpha1);

  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    if (status)
    {
      gsl_root_fsolver_free (s);
      return status;
    }

    *alpha = gsl_root_fsolver_root (s);
    alpha0 = gsl_root_fsolver_x_lower (s);
    alpha1 = gsl_root_fsolver_x_upper (s);

    status = gsl_root_test_residual (GSL_FN_EVAL ((&F), *alpha), prec);
    
    if (status == GSL_CONTINUE && (gsl_root_test_interval (alpha0, alpha1, 0, prec) == GSL_SUCCESS))
    {
      gsl_root_fsolver_free (s);
      return status;
    }
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  if (iter >= max_iter)
    return -1;

  gsl_root_fsolver_free (s);

  return 0;
}

/**
 * nc_hipert_adiab_wkb_maxtime:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha0: the initial log-redshift time.
 * @alpha1: the final log-redshift time.
 * 
 * Search for the root of $\nu_A^2$ between $\alpha_0$ and $\alpha_1$. 
 * 
 * Returns: the root of $\nu_A^2$ between $\alpha_0$ and $\alpha_1$ or NaN if not found.
 */
gdouble 
nc_hipert_adiab_wkb_maxtime (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha0, gdouble alpha1)
{
  gdouble alpha = alpha1;
  guint ret = _nc_hipert_adiab_wkb_nu2_root (pa, cosmo, alpha0, &alpha, &_nc_hipert_adiab_wkb_nuA2);

  if (ret == 0)
    return alpha;
  else
    return GSL_NAN;
}

/**
 * nc_hipert_adiab_wkb_maxtime_prec:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha0: the initial log-redshift time.
 * @alpha1: the final log-redshift time.
 * 
 * Search for the instant at which the WKB approximation starts to fails within the asked precision. 
 * 
 * Returns: the instant $\alpha$ between $\alpha_0$ and $\alpha_1$ or NaN if not found.
 */
gdouble 
nc_hipert_adiab_wkb_maxtime_prec (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha0, gdouble alpha1)
{
  gdouble alpha = alpha1;
  guint ret = _nc_hipert_adiab_wkb_nu2_root (pa, cosmo, alpha0, &alpha, &_nc_hipert_adiab_wkb_nuA2);

  if (ret == 0)
  {
    ret = _nc_hipert_adiab_wkb_nu2_root (pa, cosmo, alpha0, &alpha, &_nc_hipert_adiab_wkb_prec);
    if (ret == 0)
      return alpha;
    else
      return GSL_NAN;
  }
  else
    return GSL_NAN;
}

/**
 * nc_hipert_adiab_wkb_v:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha: the log-redshift time.
 * @Re_v: (out caller-allocates): Real part of the wkb solution.
 * @Im_v: (out caller-allocates): Imaginary part of the wkb solution.
 * 
 * Computes the WKB solution $v_\text{WKB}$ for the mode $k$ at the time $\alpha$. 
 * 
 * Returns: the value of $v_\text{WKB}(k, \alpha)$.
 */
void 
nc_hipert_adiab_wkb_v (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha, gdouble *Re_v, gdouble *Im_v)
{
  NcHIPert *pert = NC_HIPERT (pa);
  NcHIPertIAdiabEOM *eom = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (cosmo), alpha, pert->k);
  complex double v;
  gdouble int_theta;

  g_assert (pa->wkb_prep);

  int_theta = ncm_spline_eval (pa->wkb_phase->s, alpha); 
  
  v = cexp (-I * int_theta) / sqrt (2.0 * sqrt (eom->nu2));

  *Re_v = creal (v);
  *Im_v = cimag (v);
}

static gint
_nc_hipert_adiab_f (realtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertAdiabArg *arg = (NcHIPertAdiabArg *) f_data;
  NcHIPertIAdiabEOM *eom = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (arg->cosmo), alpha, NC_HIPERT (arg->pa)->k);
  
  NV_Ith_S (ydot, NC_HIPERT_ADIAB_RE_ZETA) = NV_Ith_S (y, NC_HIPERT_ADIAB_RE_PZETA) / eom->m;
  NV_Ith_S (ydot, NC_HIPERT_ADIAB_IM_ZETA) = NV_Ith_S (y, NC_HIPERT_ADIAB_IM_PZETA) / eom->m;

  NV_Ith_S (ydot, NC_HIPERT_ADIAB_RE_PZETA) = - eom->m * eom->nu2 * NV_Ith_S (y, NC_HIPERT_ADIAB_RE_ZETA);
  NV_Ith_S (ydot, NC_HIPERT_ADIAB_IM_PZETA) = - eom->m * eom->nu2 * NV_Ith_S (y, NC_HIPERT_ADIAB_IM_ZETA);

  return 0;
}

static gint
_nc_hipert_adiab_J (_NCM_SUNDIALS_INT_TYPE N, realtype alpha, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertAdiabArg *arg = (NcHIPertAdiabArg *) jac_data;
  NcHIPertIAdiabEOM *eom = nc_hipert_iadiab_eom (NC_HIPERT_IADIAB (arg->cosmo), alpha, NC_HIPERT (arg->pa)->k);

  DENSE_ELEM (J, NC_HIPERT_ADIAB_RE_ZETA, NC_HIPERT_ADIAB_RE_PZETA) = 1.0 / eom->m;
  DENSE_ELEM (J, NC_HIPERT_ADIAB_IM_ZETA, NC_HIPERT_ADIAB_IM_PZETA) = 1.0 / eom->m;
  DENSE_ELEM (J, NC_HIPERT_ADIAB_RE_PZETA, NC_HIPERT_ADIAB_RE_ZETA) = - eom->m * eom->nu2;
  DENSE_ELEM (J, NC_HIPERT_ADIAB_IM_PZETA, NC_HIPERT_ADIAB_IM_ZETA) = - eom->m * eom->nu2;

  return 0;
}

/**
 * nc_hipert_adiab_set_init_cond:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha_i: the log-redshift time.
 * @Re_zeta: Real part of the wkb solution.
 * @Im_zeta: Imaginary part of the wkb solution.
 * @Re_Pzeta: Real part of the wkb solution momentum.
 * @Im_Pzeta: Imaginary part of the wkb solution momentum.
 * 
 * Sets the initial conditions for the zeta evolution. 
 * 
 */
void 
nc_hipert_adiab_set_init_cond (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha_i, gdouble Re_zeta, gdouble Im_zeta, gdouble Re_Pzeta, gdouble Im_Pzeta)
{
  NcHIPert *pert = NC_HIPERT (pa);
  gint flag;

  pert->alpha0 = alpha_i;

  if (pert->y == NULL)
    pert->y = N_VNew_Serial (4);

  NV_Ith_S (pert->y, NC_HIPERT_ADIAB_RE_ZETA) = Re_zeta;
  NV_Ith_S (pert->y, NC_HIPERT_ADIAB_IM_ZETA) = Im_zeta;
  NV_Ith_S (pert->y, NC_HIPERT_ADIAB_RE_PZETA) = Re_Pzeta;
  NV_Ith_S (pert->y, NC_HIPERT_ADIAB_IM_PZETA) = Im_Pzeta;

  if (!pert->cvode_init)
  {
    flag = CVodeInit (pert->cvode, &_nc_hipert_adiab_f, alpha_i, pert->y);
    CVODE_CHECK (&flag, "CVodeInit", 1, );
    pert->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (pert->cvode, alpha_i, pert->y);
    CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }

  flag = CVodeSStolerances (pert->cvode, pert->reltol, pert->abstol);
  CVODE_CHECK(&flag, "CVodeSStolerances", 1,);

  flag = CVodeSetMaxNumSteps (pert->cvode, 1000000);
  CVODE_CHECK(&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVDense (pert->cvode, 4);
  CVODE_CHECK(&flag, "CVDense", 1, );

  flag = CVDlsSetDenseJacFn (pert->cvode, &_nc_hipert_adiab_J);
  CVODE_CHECK(&flag, "CVDlsSetDenseJacFn", 1, );  
}

/**
 * nc_hipert_adiab_set_init_cond_wkb:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alpha_i: the log-redshift time.
 * 
 * Sets the initial conditions for the zeta evolution using the value of the WKB solution at @alpha_i. 
 * 
 */
void 
nc_hipert_adiab_set_init_cond_wkb (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alpha_i)
{
  gdouble Re_zeta, Im_zeta, Re_Pzeta, Im_Pzeta;

  nc_hipert_adiab_wkb_zeta_Pzeta (pa, cosmo, alpha_i, &Re_zeta, &Im_zeta, &Re_Pzeta, &Im_Pzeta);
  nc_hipert_adiab_set_init_cond (pa, cosmo, alpha_i, Re_zeta, Im_zeta, Re_Pzeta, Im_Pzeta);
}

/**
 * nc_hipert_adiab_evolve:
 * @pa: a #NcHIPertAdiab.
 * @cosmo: a #NcHICosmo.
 * @alphaf: the final log-redshift time.
 * 
 * Evolve the system until @alphaf.
 * 
 */
void 
nc_hipert_adiab_evolve (NcHIPertAdiab *pa, NcHICosmo *cosmo, gdouble alphaf)
{
  NcHIPert *pert = NC_HIPERT (pa);
  NcHIPertAdiabArg arg;
  gdouble alpha_i = 0.0;
  gint flag;
  
  arg.cosmo = cosmo;
  arg.pa = pa;

  flag = CVodeSetUserData (pert->cvode, &arg);
  CVODE_CHECK (&flag, "CVodeSetFdata", 1, );

  flag = CVode (pert->cvode, alphaf, pert->y, &alpha_i, CV_NORMAL);
  CVODE_CHECK (&flag, "CVode[nc_hipert_adiab_evolve]", 1, );

  pert->alpha0 = alpha_i;
}

/**
 * nc_hipert_adiab_get_values:
 * @pa: a #NcHIPertAdiab.
 * @alpha_i: (out caller-allocates): Current time.
 * @Re_zeta: (out caller-allocates): Real part of the solution.
 * @Im_zeta: (out caller-allocates): Imaginary part of the solution.
 * @Re_Pzeta: (out caller-allocates): Real part of the solution momentum.
 * @Im_Pzeta: (out caller-allocates): Imaginary part of the solution momentum.
 * 
 * Get the current time and values of the numerical solution.
 * 
 */
void
nc_hipert_adiab_get_values (NcHIPertAdiab *pa, gdouble *alpha_i, gdouble *Re_zeta, gdouble *Im_zeta, gdouble *Re_Pzeta, gdouble *Im_Pzeta)
{
  NcHIPert *pert = NC_HIPERT (pa);
  
  *alpha_i = pert->alpha0;
  *Re_zeta  = NV_Ith_S (pert->y, NC_HIPERT_ADIAB_RE_ZETA);
  *Im_zeta  = NV_Ith_S (pert->y, NC_HIPERT_ADIAB_IM_ZETA);
  *Re_Pzeta = NV_Ith_S (pert->y, NC_HIPERT_ADIAB_RE_PZETA);
  *Im_Pzeta = NV_Ith_S (pert->y, NC_HIPERT_ADIAB_IM_PZETA);
}
