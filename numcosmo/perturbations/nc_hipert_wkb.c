/***************************************************************************
 *            nc_hipert_wkb.c
 *
 *  Sun August 03 20:39:05 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_wkb.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcHIPertWKB:
 *
 * WKB perturbation object.
 *
 * Generic implementation of WKB analysis.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline_func.h"
#include "nc_hipert_wkb.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <cvode/cvode.h>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

#include <gsl/gsl_roots.h>
#define SUN_DENSE_ACCESS SM_ELEMENT_D
#define SUN_BAND_ACCESS SM_ELEMENT_D
#endif /* NUMCOSMO_GIR_SCAN */

#include "perturbations/nc_hipert_private.h"

enum
{
  PROP_0,
  PROP_IMPL_TYPE,
  PROP_NUA2,
  PROP_V,
  PROP_DMNUA_NUA,
  PROP_EOM,
};

G_DEFINE_ABSTRACT_TYPE (NcHIPertWKB, nc_hipert_wkb, NC_TYPE_HIPERT)

static gdouble _nc_hipert_wkb_phase (gdouble x, gpointer userdata);

typedef struct _NcHIPertWKBArg
{
  NcmModel *model;
  NcHIPertWKB *wkb;
} NcHIPertWKBArg;

static void
nc_hipert_wkb_init (NcHIPertWKB *wkb)
{
  wkb->nuA = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());

  wkb->alpha_i = 0.0;
  wkb->alpha_f = 0.0;
  wkb->alpha_p = 0.0;

  wkb->alpha_phase = 0.0;
  wkb->cur_phase   = 0.0;

  wkb->lnF  = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  wkb->dlnF = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
}

static void
nc_hipert_wkb_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hipert_wkb_parent_class)->constructed (object);
  {
    nc_hipert_set_stiff_solver (NC_HIPERT (object), TRUE);
  }
}

static void
nc_hipert_wkb_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  /*NcHIPertWKB *wkb = NC_HIPERT_WKB (object);*/
  g_return_if_fail (NC_IS_HIPERT_WKB (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_wkb_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  /*NcHIPertWKB *wkb = NC_HIPERT_WKB (object);*/
  g_return_if_fail (NC_IS_HIPERT_WKB (object));

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_wkb_dispose (GObject *object)
{
  NcHIPertWKB *wkb = NC_HIPERT_WKB (object);

  ncm_spline_clear (&wkb->nuA);

  ncm_spline_clear (&wkb->lnF);
  ncm_spline_clear (&wkb->dlnF);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_wkb_parent_class)->dispose (object);
}

static void
nc_hipert_wkb_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_wkb_parent_class)->finalize (object);
}

static void _nc_hipert_wkb_set_mode_k (NcHIPert *pert, gdouble k);
static void _nc_hipert_wkb_set_abstol (NcHIPert *pert, gdouble abstol);
static void _nc_hipert_wkb_set_reltol (NcHIPert *pert, gdouble reltol);

static void
_nc_hipert_wkb_get_nu_V (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k, gdouble *nu, gdouble *V)
{
  g_error ("_nc_hipert_wkb_get_nu_V: not implemented by `%s'.", G_OBJECT_TYPE_NAME (wkb));
}

static void
_nc_hipert_wkb_get_mnu_dmnu (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k, gdouble *mnu, gdouble *dmnu)
{
  g_error ("_nc_hipert_wkb_get_mnu_dmnu: not implemented by `%s'.", G_OBJECT_TYPE_NAME (wkb));
}

static gdouble
_nc_hipert_wkb_get_m (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k)
{
  g_error ("_nc_hipert_wkb_get_m: not implemented by `%s'.", G_OBJECT_TYPE_NAME (wkb));

  return 0.0;
}

static gdouble
_nc_hipert_wkb_get_nu2 (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k)
{
  g_error ("_nc_hipert_wkb_get_nu2: not implemented by `%s'.", G_OBJECT_TYPE_NAME (wkb));

  return 0.0;
}

static gdouble
_nc_hipert_wkb_get_dVnu2 (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble k)
{
  g_error ("_nc_hipert_wkb_get_dVnu2: not implemented by `%s'.", G_OBJECT_TYPE_NAME (wkb));

  return 0.0;
}

static void
nc_hipert_wkb_class_init (NcHIPertWKBClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = nc_hipert_wkb_constructed;
  object_class->set_property = nc_hipert_wkb_set_property;
  object_class->get_property = nc_hipert_wkb_get_property;
  object_class->dispose      = nc_hipert_wkb_dispose;
  object_class->finalize     = nc_hipert_wkb_finalize;

  NC_HIPERT_CLASS (klass)->set_mode_k = &_nc_hipert_wkb_set_mode_k;
  NC_HIPERT_CLASS (klass)->set_abstol = &_nc_hipert_wkb_set_abstol;
  NC_HIPERT_CLASS (klass)->set_reltol = &_nc_hipert_wkb_set_reltol;

  klass->get_nu_V     = &_nc_hipert_wkb_get_nu_V;
  klass->get_mnu_dmnu = &_nc_hipert_wkb_get_mnu_dmnu;
  klass->get_m        = &_nc_hipert_wkb_get_m;
  klass->get_nu2      = &_nc_hipert_wkb_get_nu2;
  klass->get_dVnu2    = &_nc_hipert_wkb_get_dVnu2;
}

static void
_nc_hipert_wkb_set_mode_k (NcHIPert *pert, gdouble k)
{
  NC_HIPERT_CLASS (nc_hipert_wkb_parent_class)->set_mode_k (pert, k);
  /* Chain up : start */

/*
 *  if (!nc_hipert_prepared (pert))
 *  {
 *   NcHIPertWKB *wkb = NC_HIPERT_WKB (pert);
 *   nc_hipert_set_prepared (pert, TRUE);
 *  }
 */
}

static void
_nc_hipert_wkb_set_abstol (NcHIPert *pert, gdouble abstol)
{
  NC_HIPERT_CLASS (nc_hipert_wkb_parent_class)->set_abstol (pert, abstol);
  /* Chain up : start */

/*
 *  if (!nc_hipert_prepared (pert))
 *  {
 *   NcHIPertWKB *wkb = NC_HIPERT_WKB (pert);
 *   nc_hipert_set_prepared (pert, TRUE);
 *  }
 */
}

static void
_nc_hipert_wkb_set_reltol (NcHIPert *pert, gdouble reltol)
{
  NC_HIPERT_CLASS (nc_hipert_wkb_parent_class)->set_reltol (pert, reltol);
  /* Chain up : start */

/*
 *  if (!nc_hipert_prepared (pert))
 *  {
 *   NcHIPertWKB *wkb = NC_HIPERT_WKB (pert);
 *   nc_hipert_set_prepared (pert, TRUE);
 *  }
 */
}

/**
 * nc_hipert_wkb_ref:
 * @wkb: a #NcHIPertWKB.
 *
 * Increases the reference count of @wkb.
 *
 * Returns: (transfer full): @wkb.
 */
NcHIPertWKB *
nc_hipert_wkb_ref (NcHIPertWKB *wkb)
{
  return g_object_ref (wkb);
}

/**
 * nc_hipert_wkb_free:
 * @wkb: a #NcHIPertWKB.
 *
 * Decreases the reference count of @wkb.
 *
 */
void
nc_hipert_wkb_free (NcHIPertWKB *wkb)
{
  g_object_unref (wkb);
}

/**
 * nc_hipert_wkb_clear:
 * @wkb: a #NcHIPertWKB.
 *
 * Decreases the reference count of *@wkb and sets *@wkb to NULL.
 *
 */
void
nc_hipert_wkb_clear (NcHIPertWKB **wkb)
{
  g_clear_object (wkb);
}

/**
 * nc_hipert_wkb_set_interval:
 * @wkb: a #NcHIPertWKB
 * @alpha_i: initial log-redshift time
 * @alpha_f: final log-redshift time
 *
 * Sets the interval to calculate the WKB modes to $(\alpha_i,\,\alpha_f)$.
 *
 */
void
nc_hipert_wkb_set_interval (NcHIPertWKB *wkb, gdouble alpha_i, gdouble alpha_f)
{
  if ((wkb->alpha_i != alpha_i) || (wkb->alpha_f != alpha_f))
  {
    g_assert_cmpfloat (alpha_f, >, alpha_i);

    wkb->alpha_i = alpha_i;
    wkb->alpha_f = alpha_f;

    nc_hipert_set_prepared (NC_HIPERT (wkb), FALSE);
  }
}

/**
 * nc_hipert_wkb_get_nu_V: (virtual get_nu_V)
 * @wkb: a #NcHIPertWKB
 * @model: a #NcmModel
 * @alpha: log-redshift time
 * @k: mode $k$
 * @nu: (out): frequency $\nu$
 * @V: (out): WKB potential
 *
 * FIXME
 *
 */
/**
 * nc_hipert_wkb_get_mnu_dmnu: (virtual get_mnu_dmnu)
 * @wkb: a #NcHIPertWKB
 * @model: a #NcmModel
 * @alpha: log-redshift time
 * @k: mode $k$
 * @mnu: (out): mass-frequency $m\nu$
 * @dmnu: (out):  mass-frequency derivative $\mathrm{d}m\nu/\mathrm{d}\alpha$
 *
 * FIXME
 *
 */
/**
 * nc_hipert_wkb_get_m: (virtual get_m)
 * @wkb: a #NcHIPertWKB
 * @model: a #NcmModel
 * @alpha: log-redshift time
 * @k: mode $k$
 *
 * FIXME
 *
 * Returns: mass $m$.
 */
/**
 * nc_hipert_wkb_get_nu2: (virtual get_nu2)
 * @wkb: a #NcHIPertWKB
 * @model: a #NcmModel
 * @alpha: log-redshift time
 * @k: mode $k$
 *
 * FIXME
 *
 * Returns: frequency squared $m^2$.
 */
/**
 * nc_hipert_wkb_get_dVnu2: (virtual get_dVnu2)
 * @wkb: a #NcHIPertWKB
 * @model: a #NcmModel
 * @alpha: log-redshift time
 * @k: mode $k$
 *
 * FIXME
 *
 * Returns: derivative of the WKB potential $\mathrm{d}(V/\nu^2)/\mathrm{d}\alpha$.
 */

static gdouble
_nc_hipert_wkb_phase (gdouble alpha, gpointer userdata)
{
  NcHIPertWKBArg *arg = (NcHIPertWKBArg *) userdata;
  gdouble nu, V;

  nc_hipert_wkb_get_nu_V (arg->wkb, arg->model, alpha, nc_hipert_get_mode_k (NC_HIPERT (arg->wkb)), &nu, &V);

  return nu;
}

void
_nc_hipert_wkb_prepare_approx (NcHIPertWKB *wkb, NcmModel *model)
{
  NcHIPertWKBArg arg;
  gsl_function F;

  arg.model = model;
  arg.wkb   = wkb;

  F.params   = &arg;
  F.function = &_nc_hipert_wkb_phase;

  ncm_spline_set_func (wkb->nuA, NCM_SPLINE_FUNCTION_SPLINE, &F, wkb->alpha_i, wkb->alpha_f, 1000000000, nc_hipert_get_reltol (NC_HIPERT (wkb)));
}

static gint
_nc_hipert_wkb_phase_f (sunrealtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
/*  NcHIPertWKBArg *arg   = (NcHIPertWKBArg *) f_data; */
/*  gdouble nu = 0.0, V = 0.0; */

  /*nc_hipert_wkb_get_nu_V (arg->wkb, arg->model, alpha, nc_hipert_get_mode_k (NC_HIPERT (arg->wkb)), &nu, &dnu_nu, &V); */
  {
    /*const gdouble rnu   = NV_Ith_S (y, 0); */
    const gdouble U = NV_Ith_S (y, 1);

/*    const gdouble Rnu   = exp (rnu); */
/*    const gdouble Rnunu = Rnu / nu; */
/*    const gdouble nu2   = nu * nu; */

    NV_Ith_S (ydot, 0) = U;
    NV_Ith_S (ydot, 1) = U * U;

    return 0;
  }
}

static gint
_nc_hipert_wkb_phase_J (sunrealtype alpha, N_Vector y, N_Vector fy, SUNMatrix J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertWKBArg *arg = (NcHIPertWKBArg *) jac_data;
  gdouble nu = 0.0, V = 0.0;

  nc_hipert_wkb_get_nu_V (arg->wkb, arg->model, alpha, nc_hipert_get_mode_k (NC_HIPERT (arg->wkb)), &nu, &V);
  {
    const gdouble rnu   = NV_Ith_S (y, 0);
    const gdouble U     = NV_Ith_S (y, 1);
    const gdouble Rnu   = exp (rnu);
    const gdouble Rnunu = Rnu / nu;
    const gdouble nu2   = nu * nu;

    SUN_DENSE_ACCESS (J, 0, 0) = -Rnunu * U;
    SUN_DENSE_ACCESS (J, 0, 1) = -Rnunu;

    SUN_DENSE_ACCESS (J, 1, 0) = -(1.0 / Rnunu) * (nu2 * expm1 (4.0 * rnu) + V) + (4.0 * nu2 / Rnunu) * gsl_pow_4 (Rnu);
    SUN_DENSE_ACCESS (J, 1, 1) = 0.0;

    return 0;
  }
}

void
_nc_hipert_wkb_prepare_exact (NcHIPertWKB *wkb, NcmModel *model)
{
  NcHIPert *pert                = NC_HIPERT (wkb);
  NcHIPertPrivate * const pself = nc_hipert_get_private (pert);
  gint flag;
  gdouble alpha = wkb->alpha_i;
  GArray *alpha_a;
  GArray *lnF;
  GArray *dlnF;
  GArray *alpha_nuA_a;
  GArray *nuA_a;
  NcmVector *v = NULL;

  v = ncm_spline_get_xv (wkb->nuA);
  g_assert (v != NULL);

  alpha_nuA_a = ncm_vector_get_array (v);
  ncm_vector_free (v);

  v = ncm_spline_get_yv (wkb->nuA);
  g_assert (v != NULL);

  nuA_a = ncm_vector_get_array (v);
  ncm_vector_free (v);

  g_assert_cmpfloat (g_array_index (alpha_nuA_a, gdouble, alpha_nuA_a->len - 1), ==, wkb->alpha_i);

  if ((v = ncm_spline_get_xv (wkb->lnF)) != NULL)
  {
    alpha_a = ncm_vector_get_array (v);
    ncm_vector_free (v);

    v   = ncm_spline_get_yv (wkb->lnF);
    lnF = ncm_vector_get_array (v);
    ncm_vector_free (v);

    v    = ncm_spline_get_yv (wkb->dlnF);
    dlnF = ncm_vector_get_array (v);
    ncm_vector_free (v);

    g_array_set_size (alpha_a, 0);
    g_array_set_size (lnF, 0);
    g_array_set_size (dlnF, 0);
  }
  else
  {
    alpha_a = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    lnF     = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    dlnF    = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  }


  nc_hipert_set_sys_size (pert, 2);

  {
    gdouble nu, nu2, V, dVnu2;

    nc_hipert_wkb_get_nu_V (wkb, model, alpha, nc_hipert_get_mode_k (pert), &nu, &V);

    nu2   = nu * nu;
    dVnu2 = nc_hipert_wkb_get_dVnu2 (wkb, model, alpha, nc_hipert_get_mode_k (pert));

    NV_Ith_S (pself->y, 0) = 0.25 * log1p (V / nu2);
    NV_Ith_S (pself->y, 1) = 0.25 * (dVnu2 / (1.0 + V / nu2));
  }

  g_array_append_val (alpha_a, wkb->alpha_i);
  g_array_append_val (lnF, NV_Ith_S (pself->y, 0));
  g_array_append_val (dlnF, NV_Ith_S (pself->y, 1));

  if (!pself->cvode_init)
  {
    flag = CVodeInit (pself->cvode, &_nc_hipert_wkb_phase_f, wkb->alpha_i, pself->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
    pself->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (pself->cvode, wkb->alpha_i, pself->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }

  flag = CVodeSStolerances (pself->cvode, pself->reltol, pself->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );

  flag = CVodeSetMaxNumSteps (pself->cvode, 1000000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVodeSetLinearSolver (pself->cvode, pself->LS, pself->A);
  NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );

  flag = CVodeSetJacFn (pself->cvode, &_nc_hipert_wkb_phase_J);
  NCM_CVODE_CHECK (&flag, "CVodeSetJacFn", 1, );

  {
    gdouble last = wkb->alpha_i;
    NcHIPertWKBArg arg;

    arg.model = model;
    arg.wkb   = wkb;

    flag = CVodeSetUserData (pself->cvode, &arg);
    NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );

    while (alpha < wkb->alpha_f)
    {
      gdouble mnu, dmnu;

      flag = CVode (pself->cvode, wkb->alpha_f, pself->y, &alpha, CV_ONE_STEP);
      NCM_CVODE_CHECK (&flag, "CVode[_nc_hipert_wkb_prepare_exact]", 1, );

      if (fabs (2.0 * (alpha - last) / (alpha + last)) > 1.0e-6)
      {
        const gdouble rnu   = NV_Ith_S (pself->y, 0);
        const gdouble Unu   = NV_Ith_S (pself->y, 1);
        const gdouble Rnu   = exp (rnu);
        const gdouble m_i   = nc_hipert_wkb_get_m (wkb, model, alpha, nc_hipert_get_mode_k (pert));
        const gdouble nu_i  = sqrt (nc_hipert_wkb_get_nu2 (wkb, model, alpha, nc_hipert_get_mode_k (pert)));
        const gdouble nuA_i = Rnu * Rnu * nu_i;

        const gdouble lnF_i = log (m_i * nuA_i);
        gdouble dlnF_i;

        nc_hipert_wkb_get_mnu_dmnu (wkb, model, alpha, nc_hipert_get_mode_k (pert), &mnu, &dmnu);

        dlnF_i = dmnu / mnu - 2.0 * Rnu * Unu / nu_i;

        g_array_append_val (alpha_a, alpha);
        g_array_append_val (lnF, lnF_i);
        g_array_append_val (dlnF, dlnF_i);

        g_array_append_val (alpha_nuA_a, alpha);
        g_array_append_val (nuA_a, nuA_i);

        last = alpha;
      }
    }
  }

  ncm_spline_set_array (wkb->lnF, alpha_a, lnF, TRUE);
  ncm_spline_set_array (wkb->dlnF, alpha_a, dlnF, TRUE);
  ncm_spline_set_array (wkb->nuA, alpha_nuA_a, nuA_a, TRUE);

  g_array_unref (alpha_a);
  g_array_unref (lnF);
  g_array_unref (dlnF);
  g_array_unref (alpha_nuA_a);
  g_array_unref (nuA_a);
}

/**
 * nc_hipert_wkb_prepare:
 * @wkb: a #NcHIPertWKB
 * @model: a #NcmModel
 *
 * Prepare the object for WKB calculations using the model @model. It uses the wkb
 * approximation until @wkb->reltol is reached and then it solves the non-linear equation of motion
 * for $\nu_A$ for the rest of the interval.
 *
 */
void
nc_hipert_wkb_prepare (NcHIPertWKB *wkb, NcmModel *model)
{
  NcHIPert *pert     = NC_HIPERT (wkb);
  const gdouble prec = nc_hipert_get_reltol (pert);
  gdouble nu_i, V_i, nu_f, V_f;

  if (!nc_hipert_prepared (pert))
  {
    gdouble nu2_i, nu2_f;

    nc_hipert_wkb_get_nu_V (wkb, model, wkb->alpha_i, nc_hipert_get_mode_k (pert), &nu_i, &V_i);
    nc_hipert_wkb_get_nu_V (wkb, model, wkb->alpha_f, nc_hipert_get_mode_k (pert), &nu_f, &V_f);

    nu2_i = nu_i * nu_i;
    nu2_f = nu_f * nu_f;

    if (fabs (nu2_i / V_i) < 1.0)
    {
      g_error ("nc_hipert_wkb_prepare: cannot prepare wkb solution in the interval [% 7.5g % 7.5g], the WKB approximation is not valid at % 7.5g (nuA2_i = % 7.5g, nuA2_i = % 7.5g).",
               wkb->alpha_i, wkb->alpha_f, wkb->alpha_i, nu2_i, nu2_f);
    }
    else
    {
      if (prec < fabs (V_i / nu2_i))
      {
        g_error ("nc_hipert_wkb_prepare: cannot prepare wkb solution in the interval [% 7.5g % 7.5g], this interval is beyond the desired precision [%7.5e].",
                 wkb->alpha_i, wkb->alpha_f, prec);
      }
      else
      {
        if (fabs (V_f / nu2_f) < prec)
        {
          _nc_hipert_wkb_prepare_approx (wkb, model);
          wkb->alpha_p = wkb->alpha_f;
        }
        else
        {
          const gdouble alpha_p = nc_hipert_wkb_maxtime_prec (wkb, model, NC_HIPERT_WKB_CMP_POTENTIAL, wkb->alpha_i, wkb->alpha_f);

          if (!gsl_finite (alpha_p))
          {
            g_error ("nc_hipert_wkb_prepare: cannot find the precision [%7.5e] point in the interval [% 7.5g % 7.5g].",
                     prec, wkb->alpha_i, wkb->alpha_f);
          }
          else
          {
            wkb->alpha_p = alpha_p;
            /*printf ("# Preparing approx [% 21.16g % 21.16g]\n", alpha_i, alpha_p);*/
            _nc_hipert_wkb_prepare_approx (wkb, model);
            /*printf ("# Preparing exact [% 21.16g % 21.16g]\n", alpha_p, alpha_f);*/
            _nc_hipert_wkb_prepare_exact (wkb, model);
          }
        }
      }
    }

    wkb->alpha_phase = wkb->alpha_i;
    wkb->cur_phase   = M_PI * 0.25;
    nc_hipert_set_prepared (pert, TRUE);
  }
}

/**
 * nc_hipert_wkb_q:
 * @wkb: a #NcHIPertWKB
 * @model: a #NcmModel
 * @alpha: the log-redshift time
 * @Re_q: (out caller-allocates): Real part of the wkb solution
 * @Im_q: (out caller-allocates): Imaginary part of the wkb solution
 *
 * Computes the WKB solution $q_\text{WKB}$ for the mode $k$ at the time $\alpha$.
 *
 */
void
nc_hipert_wkb_q (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble *Re_q, gdouble *Im_q)
{
  NcHIPert *pert = NC_HIPERT (wkb);
  complex double q;
  const gdouble int_nuA = nc_hipert_wkb_phase (wkb, model, alpha);

  if (alpha < wkb->alpha_p)
  {
    const gdouble nuA = ncm_spline_eval (wkb->nuA, alpha);
    gdouble m         = nc_hipert_wkb_get_m (wkb, model, alpha, nc_hipert_get_mode_k (pert));

    q = cexp (-I * int_nuA) / sqrt (2.0 * m * nuA);
  }
  else
  {
    gdouble lnF;
    const gdouble one_sqrt2 = 1.0 / sqrt (2.0);

    lnF = ncm_spline_eval (wkb->lnF, alpha);

    q = cexp (-I * int_nuA - lnF * 0.5) * one_sqrt2;
  }

  *Re_q = creal (q);
  *Im_q = cimag (q);
}

/**
 * nc_hipert_wkb_q_p:
 * @wkb: a #NcHIPertWKB.
 * @model: a #NcmModel
 * @alpha: the log-redshift time.
 * @Re_q: (out caller-allocates): Real part of the wkb solution.
 * @Im_q: (out caller-allocates): Imaginary part of the wkb solution.
 * @Re_p: (out caller-allocates): Real part of the wkb solution momentum.
 * @Im_p: (out caller-allocates): Imaginary part of the wkb solution momentum.
 *
 * Computes the WKB solution $q_\text{WKB}$ and its momentum for the mode $k$ at the time $\alpha$.
 *
 */
void
nc_hipert_wkb_q_p (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha, gdouble *Re_q, gdouble *Im_q, gdouble *Re_p, gdouble *Im_p)
{
  NcHIPert *pert = NC_HIPERT (wkb);
  complex double q, p;
  gdouble int_nuA = nc_hipert_wkb_phase (wkb, model, alpha);

  if (alpha < wkb->alpha_p)
  {
    gdouble mnu, dmnu;
    const gdouble nuA = ncm_spline_eval (wkb->nuA, alpha);
    const gdouble m   = nc_hipert_wkb_get_m (wkb, model, alpha, nc_hipert_get_mode_k (pert));

    nc_hipert_wkb_get_mnu_dmnu (wkb, model, alpha, nc_hipert_get_mode_k (pert), &mnu, &dmnu);

    q = cexp (-I * int_nuA) / sqrt (2.0 * m * nuA);
    p = -I *cexp (-I *int_nuA) * sqrt (0.5 * m * nuA) - 0.5 * m * (dmnu / mnu) * q;
  }
  else
  {
    const gdouble m   = nc_hipert_wkb_get_m (wkb, model, alpha, nc_hipert_get_mode_k (pert));
    const gdouble lnF = ncm_spline_eval (wkb->lnF, alpha);

    ;

    const gdouble dlnF      = ncm_spline_eval (wkb->dlnF, alpha);
    const gdouble one_sqrt2 = 1.0 / sqrt (2.0);

    q = cexp (-I * int_nuA - lnF * 0.5) * one_sqrt2;
    p = -I *cexp (-I *int_nuA + lnF * 0.5) * one_sqrt2 - 0.5 * m * dlnF * q;
  }

  *Re_q = creal (q);
  *Im_q = cimag (q);

  *Re_p = creal (p);
  *Im_p = cimag (p);
}

static gdouble
_nc_hipert_wkb_prec (gdouble alpha, gpointer userdata)
{
  NcHIPertWKBArg *arg = (NcHIPertWKBArg *) userdata;
  gdouble nu, V;

  nc_hipert_wkb_get_nu_V (arg->wkb, arg->model, alpha, nc_hipert_get_mode_k (NC_HIPERT (arg->wkb)), &nu, &V);
  {
    const gdouble test = log (fabs (nu * nu * nc_hipert_get_reltol (NC_HIPERT (arg->wkb)) / V));

    return test;
  }
}

static gdouble
_nc_hipert_wkb_prec_alpha2 (gdouble alpha, gpointer userdata)
{
  NcHIPertWKBArg *arg = (NcHIPertWKBArg *) userdata;
  gdouble nu, V;

  nc_hipert_wkb_get_nu_V (arg->wkb, arg->model, alpha, nc_hipert_get_mode_k (NC_HIPERT (arg->wkb)), &nu, &V);

  {
    const gdouble nuA2 = nu * nu - V;
    const gdouble test = log (fabs (nuA2 * nc_hipert_get_reltol (NC_HIPERT (arg->wkb)) / (alpha * alpha)));

    return test;
  }
}

static gdouble
_nc_hipert_wkb_nuA2 (gdouble alpha, gpointer userdata)
{
  NcHIPertWKBArg *arg = (NcHIPertWKBArg *) userdata;
  gdouble nu, V;

  nc_hipert_wkb_get_nu_V (arg->wkb, arg->model, alpha, nc_hipert_get_mode_k (NC_HIPERT (arg->wkb)), &nu, &V);

  {
    const gdouble nu2  = nu * nu;
    const gdouble nuA2 = nu2 - V;

    return nuA2 / nu2 - nc_hipert_get_reltol (NC_HIPERT (arg->wkb));
  }
}

static gint
_nc_hipert_wkb_nu2_root (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha0, gdouble *alpha, gsl_function *F)
{
  NcHIPert *pert                = NC_HIPERT (wkb);
  NcHIPertPrivate * const pself = nc_hipert_get_private (pert);
  gint iter                     = 0;
  gint max_iter                 = 1000000;
  gdouble prec                  = pself->reltol;
  gdouble alpha1                = *alpha;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gint status;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, F, alpha0, alpha1);

  do {
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

    status = gsl_root_test_residual (GSL_FN_EVAL (F, *alpha), prec);

    if ((status == GSL_CONTINUE) && (gsl_root_test_interval (alpha0, alpha1, 0, prec) == GSL_SUCCESS))
    {
      gsl_root_fsolver_free (s);

      return status;
    }
  } while (status == GSL_CONTINUE && iter < max_iter);

  if (iter >= max_iter)
    return -1;

  gsl_root_fsolver_free (s);

  return 0;
}

/**
 * nc_hipert_wkb_maxtime:
 * @wkb: a #NcHIPertWKB
 * @model: a #NcmModel
 * @alpha0: the initial log-redshift time
 * @alpha1: the final log-redshift time
 *
 * Search for the root of $\nu_A^2$ between $\alpha_0$ and $\alpha_1$.
 *
 * Returns: the root of $\nu_A^2$ between $\alpha_0$ and $\alpha_1$ or NaN if not found.
 */
gdouble
nc_hipert_wkb_maxtime (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha0, gdouble alpha1)
{
  gdouble alpha = alpha1;
  gsl_function F;
  NcHIPertWKBArg arg;
  guint ret;

  arg.model = model;
  arg.wkb   = wkb;

  F.function = &_nc_hipert_wkb_nuA2;
  F.params   = &arg;

  ret = _nc_hipert_wkb_nu2_root (wkb, model, alpha0, &alpha, &F);

  if (ret == 0)
    return alpha;
  else
    return GSL_NAN;
}

/**
 * nc_hipert_wkb_maxtime_prec:
 * @wkb: a #NcHIPertWKB
 * @model: a #NcmModel
 * @cmp: Comparison type
 * @alpha0: FIXME
 * @alpha1: FIXME
 *
 * Search for the instant at which the WKB approximation starts to fails within the asked precision.
 *
 * Returns: the instant $\alpha$ between $\alpha_0$ and $\alpha_1$ or NaN if not found.
 */
gdouble
nc_hipert_wkb_maxtime_prec (NcHIPertWKB *wkb, NcmModel *model, NcHIPertWKBCmp cmp, gdouble alpha0, gdouble alpha1)
{
  gdouble alpha = alpha1;
  gsl_function F;
  NcHIPertWKBArg arg;
  guint ret;

  arg.model = model;
  arg.wkb   = wkb;

  F.function = &_nc_hipert_wkb_nuA2;
  F.params   = &arg;

  ret = _nc_hipert_wkb_nu2_root (wkb, model, alpha0, &alpha, &F);

  if (ret == 0)
  {
    switch (cmp)
    {
      case NC_HIPERT_WKB_CMP_POTENTIAL:
        F.function = &_nc_hipert_wkb_prec;
        break;
      case NC_HIPERT_WKB_CMP_ALPHA2:
        F.function = &_nc_hipert_wkb_prec_alpha2;
        break;
      default:
        g_assert_not_reached ();
        break;
    }

    ret = _nc_hipert_wkb_nu2_root (wkb, model, alpha0, &alpha, &F);

    if (ret == 0)
      return alpha;
    else
      return GSL_NAN;
  }
  else
  {
    return GSL_NAN;
  }
}

/**
 * nc_hipert_wkb_nuA:
 * @wkb: a #NcHIPertWKB
 * @model: a #NcmModel
 * @alpha: the log-redshift time
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hipert_wkb_nuA (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha)
{
  g_assert (nc_hipert_prepared (NC_HIPERT (wkb)));

  return ncm_spline_eval (wkb->nuA, alpha);
}

/**
 * nc_hipert_wkb_phase:
 * @wkb: a #NcHIPertWKB
 * @model: a #NcmModel
 * @alpha: the log-redshift time
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
nc_hipert_wkb_phase (NcHIPertWKB *wkb, NcmModel *model, gdouble alpha)
{
  g_assert (nc_hipert_prepared (NC_HIPERT (wkb)));

  if (alpha == wkb->alpha_phase)
  {
    return wkb->cur_phase;
  }
  else
  {
    NcmVector *nuA_xv = ncm_spline_peek_xv (wkb->nuA);
    gdouble alpha_i   = wkb->alpha_phase;
    gdouble alpha_f   = alpha;
    gdouble sign      = 1.0;
    gdouble Dphase    = 0.0;
    gdouble dphase, alpha_k;
    guint cur, tar, k;

    if (alpha_i > alpha_f)
    {
      gdouble tmp = alpha_i;

      alpha_i = alpha_f;
      alpha_f = tmp;
      sign    = -1.0;
    }

    cur = ncm_spline_get_index (wkb->nuA, alpha_i);
    tar = ncm_spline_get_index (wkb->nuA, alpha_f);

    g_assert_cmpuint (cur, <=, tar);

    alpha_k = alpha_i;

    for (k = cur; k < tar; k++)
    {
      gdouble alpha_kp1 = ncm_vector_get (nuA_xv, k + 1);

      dphase  = ncm_spline_eval_integ (wkb->nuA, alpha_k, alpha_kp1);
      Dphase  = fmod (Dphase + dphase, 2.0 * M_PI);
      alpha_k = alpha_kp1;
    }

    dphase = ncm_spline_eval_integ (wkb->nuA, alpha_k, alpha_f);
    Dphase = fmod (Dphase + dphase, 2.0 * M_PI);

    wkb->cur_phase = fmod (wkb->cur_phase + sign * Dphase, 2.0 * M_PI);

    if (wkb->cur_phase < 0.0)
      wkb->cur_phase += 2.0 * M_PI;
  }

  wkb->alpha_phase = alpha;

  return wkb->cur_phase;
}

