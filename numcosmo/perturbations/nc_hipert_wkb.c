/***************************************************************************
 *            nc_hipert_wkb.c
 *
 *  Sun August 03 20:39:05 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_wkb.c
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
 * SECTION:nc_hipert_wkb
 * @title: NcHIPertWKB
 * @short_description: WKB perturbation object. 
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

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <nvector/nvector_serial.h> 
#include <gsl/gsl_roots.h>

enum
{
  PROP_0,
  PROP_IMPL_TYPE,
  PROP_NUA2,
  PROP_V,
  PROP_DMNUA_NUA,
  PROP_EOM,
};

G_DEFINE_TYPE (NcHIPertWKB, nc_hipert_wkb, NC_TYPE_HIPERT);

static gdouble _nc_hipert_wkb_phase (gdouble x, gpointer userdata);

typedef struct _NcHIPertWKBArg
{
  GObject *obj;
  NcHIPertWKB *wkb;
  gdouble prec;
} NcHIPertWKBArg;


static void
nc_hipert_wkb_init (NcHIPertWKB *wkb)
{
  wkb->impl_type   = G_TYPE_INVALID;
  wkb->obj         = NULL;
  wkb->nuA         = ncm_spline_cubic_notaknot_new ();

  wkb->alpha_i     = 0.0;
  wkb->alpha_f     = 0.0;
  wkb->alpha_p     = 0.0;

  wkb->alpha_phase = 0.0;
  wkb->cur_phase   = 0.0;

  wkb->lnF         = ncm_spline_cubic_notaknot_new ();
  wkb->dlnF        = ncm_spline_cubic_notaknot_new ();

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
  NcHIPertWKB *wkb = NC_HIPERT_WKB (object); 
  g_return_if_fail (NC_IS_HIPERT_WKB (object));

  switch (prop_id)
  {
    case PROP_IMPL_TYPE:
      wkb->impl_type = g_value_get_gtype (value);
      break;
    case PROP_NUA2:
      wkb->nuA2 = g_value_get_pointer (value);
      break;
    case PROP_V:
      wkb->V = g_value_get_pointer (value);
      break;
    case PROP_DMNUA_NUA:
      wkb->dmnuA_nuA = g_value_get_pointer (value);
      break;
    case PROP_EOM:
      wkb->eom = g_value_get_pointer (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_wkb_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIPertWKB *wkb = NC_HIPERT_WKB (object);
  g_return_if_fail (NC_IS_HIPERT_WKB (object));

  switch (prop_id)
  {
    case PROP_IMPL_TYPE:
      g_value_set_gtype (value, wkb->impl_type);
      break;
    case PROP_NUA2:
      g_value_set_pointer (value, wkb->nuA2);
      break;
    case PROP_V:
      g_value_set_pointer (value, wkb->V);
      break;
    case PROP_DMNUA_NUA:
      g_value_set_pointer (value, wkb->dmnuA_nuA);
      break;
    case PROP_EOM:
      g_value_set_pointer (value, wkb->eom);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_wkb_dispose (GObject *object)
{
  NcHIPertWKB *wkb = NC_HIPERT_WKB (object);

  g_clear_object (&wkb->obj);
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
nc_hipert_wkb_class_init (NcHIPertWKBClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = nc_hipert_wkb_constructed;
  object_class->set_property = nc_hipert_wkb_set_property;
  object_class->get_property = nc_hipert_wkb_get_property;
  object_class->dispose      = nc_hipert_wkb_dispose;
  object_class->finalize     = nc_hipert_wkb_finalize;

  g_object_class_install_property (object_class,
                                   PROP_IMPL_TYPE,
                                   g_param_spec_gtype ("impl-type",
                                                       NULL,
                                                       "Implementation GType",
                                                       G_TYPE_NONE,
                                                       G_PARAM_WRITABLE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  g_object_class_install_property (object_class,
                                   PROP_NUA2,
                                   g_param_spec_pointer ("nuA2",
                                                         NULL,
                                                         "nu_A^2",
                                                         G_PARAM_WRITABLE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_V,
                                   g_param_spec_pointer ("V",
                                                         NULL,
                                                         "Potential",
                                                         G_PARAM_WRITABLE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_DMNUA_NUA,
                                   g_param_spec_pointer ("dmnuA-nuA",
                                                         NULL,
                                                         "dm\nu_A/\nu_A",
                                                         G_PARAM_WRITABLE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_EOM,
                                   g_param_spec_pointer ("eom",
                                                         NULL,
                                                         "Equation of motion",
                                                         G_PARAM_WRITABLE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  NC_HIPERT_CLASS (klass)->set_mode_k = &_nc_hipert_wkb_set_mode_k;
  NC_HIPERT_CLASS (klass)->set_abstol = &_nc_hipert_wkb_set_abstol;
  NC_HIPERT_CLASS (klass)->set_reltol = &_nc_hipert_wkb_set_reltol;
}

static void 
_nc_hipert_wkb_set_mode_k (NcHIPert *pert, gdouble k) 
{
  NC_HIPERT_CLASS (nc_hipert_wkb_parent_class)->set_mode_k (pert, k);
  /* Chain up : start */
/*
  if (!pert->prepared)
  {
    NcHIPertWKB *wkb = NC_HIPERT_WKB (pert);
    pert->prepared    = TRUE;
  }
*/  
}

static void 
_nc_hipert_wkb_set_abstol (NcHIPert *pert, gdouble abstol) 
{
  NC_HIPERT_CLASS (nc_hipert_wkb_parent_class)->set_abstol (pert, abstol);
  /* Chain up : start */
/*
  if (!pert->prepared)
  {
    NcHIPertWKB *wkb = NC_HIPERT_WKB (pert);
    pert->prepared    = TRUE;
  }
*/  
}

static void 
_nc_hipert_wkb_set_reltol (NcHIPert *pert, gdouble reltol) 
{
  NC_HIPERT_CLASS (nc_hipert_wkb_parent_class)->set_reltol (pert, reltol);
  /* Chain up : start */
/*
  if (!pert->prepared)
  {
    NcHIPertWKB *wkb = NC_HIPERT_WKB (pert);
    pert->prepared    = TRUE;
  }
*/
}

/**
 * nc_hipert_wkb_new:
 * @impl_type: FIXME
 * @nuA2: (scope notified): FIXME
 * @V: (scope notified): FIXME
 * @dmnuA_nuA: (scope notified): FIXME 
 * @eom: (scope notified): FIXME
 * 
 * Creates a new #NcHIPertWKB object.
 * 
 * Returns: (transfer full): a new #NcHIPertWKB.
 */
NcHIPertWKB *
nc_hipert_wkb_new (GType impl_type, NcHIPertWKBFunc nuA2, NcHIPertWKBFunc V, NcHIPertWKBFunc dmnuA_nuA, NcHIPertWKBEom eom)
{
  NcHIPertWKB *wkb = g_object_new (NC_TYPE_HIPERT_WKB,
                                   "impl-type", impl_type,
                                   "nuA2", nuA2,
                                   "V", V,
                                   "dmnuA_nuA", dmnuA_nuA,
                                   "eom", eom,
                                   NULL);
  return wkb;
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

static gdouble 
_nc_hipert_wkb_phase (gdouble alpha, gpointer userdata)
{
  NcHIPertWKBArg *arg = (NcHIPertWKBArg *) userdata;
  const gdouble nuA2 = arg->wkb->nuA2 (arg->obj, alpha, NC_HIPERT (arg->wkb)->k);
  return sqrt (nuA2);
}

void 
_nc_hipert_wkb_prepare_approx (NcHIPertWKB *wkb, GObject *obj, gdouble alpha_i, gdouble alpha_f)
{
  NcHIPertWKBArg arg;
  gsl_function F;

  arg.obj = obj;
  arg.wkb = wkb;

  F.params = &arg;
  F.function = &_nc_hipert_wkb_phase;

  ncm_spline_set_func (wkb->nuA, NCM_SPLINE_FUNCTION_SPLINE, &F, alpha_i, alpha_f, 1000000000, NC_HIPERT (wkb)->reltol);

  wkb->alpha_i = alpha_i;
  wkb->alpha_f = alpha_f;
}

static gint
_nc_hipert_wkb_phase_f (realtype alpha, N_Vector y, N_Vector ydot, gpointer f_data)
{
  NcHIPertWKBArg *arg   = (NcHIPertWKBArg *) f_data;
  gdouble m = 0.0, nu2 = 0.0, dlnm = 0.0;
  arg->wkb->eom (arg->obj, alpha, NC_HIPERT (arg->wkb)->k, &nu2, &m, &dlnm);
  {
    const gdouble lnF     = NV_Ith_S (y, 0);
    const gdouble dlnF    = NV_Ith_S (y, 1);
    const gdouble dlnF2   = dlnF * dlnF;
    const gdouble m2      = m * m;
    
    NV_Ith_S (ydot, 0) = dlnF;
    NV_Ith_S (ydot, 1) = (0.5 * dlnF2 + 2.0 * (nu2 - exp (2.0 * lnF) / m2)) - dlnF * dlnm;

    return 0;
  }
}

static gint
_nc_hipert_wkb_phase_J (_NCM_SUNDIALS_INT_TYPE N, realtype alpha, N_Vector y, N_Vector fy, DlsMat J, gpointer jac_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  NcHIPertWKBArg *arg  = (NcHIPertWKBArg *) jac_data;
  gdouble m = 0.0, nu2 = 0.0, dlnm = 0.0;
  arg->wkb->eom (arg->obj, alpha, NC_HIPERT (arg->wkb)->k, &nu2, &m, &dlnm);
  {
    const gdouble lnF     = NV_Ith_S (y, 0);
    const gdouble dlnF    = NV_Ith_S (y, 1);
    const gdouble m2      = m * m; 

    DENSE_ELEM (J, 0, 1) = 1.0;

    DENSE_ELEM (J, 1, 0) = - 4.0 * exp (2.0 * lnF) / m2;
    DENSE_ELEM (J, 1, 1) = dlnF - dlnm;
    
    return 0;
  }
}

void 
_nc_hipert_wkb_prepare_exact (NcHIPertWKB *wkb, GObject *obj, gdouble alpha_i, gdouble alpha_f)
{
  NcHIPert *pert = NC_HIPERT (wkb);
  gint flag;
  gdouble alpha = alpha_i;
  gdouble m = 0.0, nu2 = 0.0, dlnm = 0.0;
  const gdouble nuA2 = wkb->nuA2 (obj, alpha, pert->k);
  const gdouble nuA = sqrt (nuA2);
  const gdouble dmnuA_nuA = wkb->dmnuA_nuA (obj, alpha, pert->k);
  const gdouble dmnuA = dmnuA_nuA * nuA;
  GArray *alpha_a;
  GArray *lnF;
  GArray *dlnF;
  GArray *alpha_nuA_a;
  GArray *nuA_a;
  
  NcmVector *v = NULL;

  wkb->eom (obj, alpha, pert->k, &nu2, &m, &dlnm);

  v = ncm_spline_get_xv (wkb->nuA);
  g_assert (v != NULL);
  alpha_nuA_a = ncm_vector_get_array (v);
  ncm_vector_free (v);

  v = ncm_spline_get_yv (wkb->nuA);
  g_assert (v != NULL);
  nuA_a = ncm_vector_get_array (v);
  ncm_vector_free (v);

  g_assert_cmpfloat (g_array_index (alpha_nuA_a, gdouble, alpha_nuA_a->len - 1), ==, alpha_i);
  
  if ((v = ncm_spline_get_xv (wkb->lnF)) != NULL)
  {
    alpha_a = ncm_vector_get_array (v);
    ncm_vector_free (v);

    v = ncm_spline_get_yv (wkb->lnF);
    lnF = ncm_vector_get_array (v);
    ncm_vector_free (v);

    v = ncm_spline_get_yv (wkb->dlnF);
    dlnF = ncm_vector_get_array (v);
    ncm_vector_free (v);

    g_array_set_size (alpha_a, 0);
    g_array_set_size (lnF, 0);
    g_array_set_size (dlnF, 0);
  }
  else
  {
    alpha_a  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    lnF  = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
    dlnF = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), 1000);
  }

  if (pert->y == NULL)
    pert->y = N_VNew_Serial (2);

  NV_Ith_S (pert->y, 0) = log (nuA * m);
  NV_Ith_S (pert->y, 1) = dmnuA / (nuA * m);

  g_array_append_val (alpha_a, alpha_i);
  g_array_append_val (lnF, NV_Ith_S (pert->y, 0));
  g_array_append_val (dlnF, NV_Ith_S (pert->y, 1));

  if (!pert->cvode_init)
  {
    flag = CVodeInit (pert->cvode, &_nc_hipert_wkb_phase_f, alpha_i, pert->y);
    NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
    pert->cvode_init = TRUE;
  }
  else
  {
    flag = CVodeReInit (pert->cvode, alpha_i, pert->y);
    NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
  }

  flag = CVodeSStolerances (pert->cvode, pert->reltol, pert->abstol);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1,);

  flag = CVodeSetMaxNumSteps (pert->cvode, 1000000);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );

  flag = CVDense (pert->cvode, 2);
  NCM_CVODE_CHECK (&flag, "CVDense", 1, );

  flag = CVDlsSetDenseJacFn (pert->cvode, &_nc_hipert_wkb_phase_J);
  NCM_CVODE_CHECK (&flag, "CVDlsSetDenseJacFn", 1, );  

  {
    gdouble last = alpha_i;
    NcHIPertWKBArg arg;
    arg.obj = obj;
    arg.wkb = wkb;

    flag = CVodeSetUserData (pert->cvode, &arg);
    NCM_CVODE_CHECK (&flag, "CVodeSetFdata", 1, );

    while (alpha < alpha_f)
    {
      flag = CVode (pert->cvode, alpha_f, pert->y, &alpha, CV_ONE_STEP);
      NCM_CVODE_CHECK (&flag, "CVode[_nc_hipert_wkb_prepare_exact]", 1, );

      if (fabs (2.0 * (alpha - last) / (alpha + last)) > 1.0e-6)
      {
        gdouble nuA_k;
        wkb->eom (obj, alpha, pert->k, &nu2, &m, &dlnm);
        nuA_k = exp (NV_Ith_S (pert->y, 0)) / m;
        
        g_array_append_val (alpha_a, alpha);
        g_array_append_val (lnF, NV_Ith_S (pert->y, 0));
        g_array_append_val (dlnF, NV_Ith_S (pert->y, 1));

        g_array_append_val (alpha_nuA_a, alpha);
        g_array_append_val (nuA_a, nuA_k);

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
 * @wkb: a #NcHIPertWKB.
 * @obj: the WKB implementation.
 * @prec: Required precision.
 * @alpha_i: initial log-redshift time.
 * @alpha_f: final log-redshift time.
 * 
 * Prepare the object for WKB calculations using the implementation @obj. It uses the wkb 
 * approximation until @prec is reached and then it solves the non-linear equation of motion 
 * for $\nu_A$ for the rest of the interval. 
 * 
 */
void 
nc_hipert_wkb_prepare (NcHIPertWKB *wkb, GObject *obj, gdouble prec, gdouble alpha_i, gdouble alpha_f)
{
  NcHIPert *pert = NC_HIPERT (wkb);
  if (wkb->obj != obj)
  {
    g_clear_object (&wkb->obj);
    g_assert (g_type_is_a (G_OBJECT_TYPE (obj), wkb->impl_type));
    wkb->obj = g_object_ref (obj);
  }

  if (!pert->prepared || (alpha_i != wkb->alpha_i || alpha_f != wkb->alpha_f))
  {
    const gdouble nuA2_i = wkb->nuA2 (obj, alpha_i, pert->k);
    const gdouble nuA2_f = wkb->nuA2 (obj, alpha_f, pert->k);
    if (nuA2_i < 0.0)
    {
      g_error ("nc_hipert_wkb_prepare_patched: cannot prepare wkb solution in the interval [% 7.5g % 7.5g], the WKB approximation is not valid at % 7.5g (nuA2_i = % 7.5g, nuA2_i = % 7.5g).",
               alpha_i, alpha_f, alpha_i, nuA2_i, nuA2_f);
    }
    else
    {
      gdouble m = 0.0, nu2 = 0.0, dlnm = 0.0;
      const gdouble V_i = wkb->V (obj, alpha_i, pert->k);
      wkb->eom (obj, alpha_i, pert->k, &nu2, &m, &dlnm);
      
      if (prec < fabs (V_i / nu2))
      {
        g_error ("nc_hipert_wkb_prepare_patched: cannot prepare wkb solution in the interval [% 7.5g % 7.5g], this interval is beyond the desired precision [%7.5e].",
                 alpha_i, alpha_f, prec);
      }
      else
      {
        const gdouble V_f = wkb->V (obj, alpha_f, pert->k);
        wkb->eom (obj, alpha_f, pert->k, &nu2, &m, &dlnm);
        if (fabs (V_f / nu2) < prec)
        {
          _nc_hipert_wkb_prepare_approx (wkb, obj, alpha_i, alpha_f);
          wkb->alpha_p = alpha_f;
        }
        else
        {
          const gdouble alpha_p = nc_hipert_wkb_maxtime_prec (wkb, obj, NC_HIPERT_WKB_CMP_POTENTIAL, prec, alpha_i, alpha_f);
          if (!gsl_finite (alpha_p))
          {
            g_error ("nc_hipert_wkb_prepare_patched: cannot find the precision [%7.5e] point in the interval [% 7.5g % 7.5g].", 
                     prec, alpha_i, alpha_f);
          }
          else
          {
            wkb->alpha_p = alpha_p;
            /*printf ("# Preparing approx [% 21.16g % 21.16g]\n", alpha_i, alpha_p);*/
            _nc_hipert_wkb_prepare_approx (wkb, obj, alpha_i, alpha_p);
            /*printf ("# Preparing exact [% 21.16g % 21.16g]\n", alpha_p, alpha_f);*/
            _nc_hipert_wkb_prepare_exact (wkb, obj, alpha_p, alpha_f);
            printf ("# Prepared [% 21.16g % 21.16g] %zu\n", alpha_i, alpha_f, wkb->nuA->len);
          }
        }
      }
      
    } 

    wkb->alpha_phase = alpha_i;
    wkb->cur_phase   = M_PI * 0.25;
    wkb->alpha_i     = alpha_i;
    wkb->alpha_f     = alpha_f;
    pert->prepared   = TRUE;
  }
}

/**
 * nc_hipert_wkb_q:
 * @wkb: a #NcHIPertWKB.
 * @obj: the WKB implementation.
 * @alpha: the log-redshift time.
 * @Re_q: (out caller-allocates): Real part of the wkb solution.
 * @Im_q: (out caller-allocates): Imaginary part of the wkb solution.
 * 
 * Computes the WKB solution $q_\text{WKB}$ for the mode $k$ at the time $\alpha$. 
 * 
 */
void
nc_hipert_wkb_q (NcHIPertWKB *wkb, GObject *obj, gdouble alpha, gdouble *Re_q, gdouble *Im_q)
{
  if (wkb->obj != obj)
  {
    g_clear_object (&wkb->obj);
    g_assert (g_type_is_a (G_OBJECT_TYPE (obj), wkb->impl_type));
    wkb->obj = g_object_ref (obj);
  }
  
  {
    NcHIPert *pert = NC_HIPERT (wkb);
    complex double q;
    const gdouble int_nuA = nc_hipert_wkb_phase (wkb, obj, alpha);
    
    if (alpha < wkb->alpha_p)
    {
      const gdouble nuA = ncm_spline_eval (wkb->nuA, alpha); 
      gdouble m = 0.0, nu2 = 0.0, dlnm = 0.0;
      wkb->eom (obj, alpha, pert->k, &nu2, &m, &dlnm);

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
}

/**
 * nc_hipert_wkb_q_p:
 * @wkb: a #NcHIPertWKB.
 * @obj: the WKB implementation.
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
nc_hipert_wkb_q_p (NcHIPertWKB *wkb, GObject *obj, gdouble alpha, gdouble *Re_q, gdouble *Im_q, gdouble *Re_p, gdouble *Im_p)
{
  if (wkb->obj != obj)
  {
    g_clear_object (&wkb->obj);
    g_assert (g_type_is_a (G_OBJECT_TYPE (obj), wkb->impl_type));
    wkb->obj = g_object_ref (obj);
  }
  {
    NcHIPert *pert = NC_HIPERT (wkb);
    complex double q, p;
    gdouble int_nuA = nc_hipert_wkb_phase (wkb, obj, alpha);
    gdouble m = 0.0, nu2 = 0.0, dlnm = 0.0;
    wkb->eom (obj, alpha, pert->k, &nu2, &m, &dlnm);

    if (alpha < wkb->alpha_p)
    {
      const gdouble nuA2 = wkb->nuA2 (obj, alpha, pert->k);
      const gdouble nuA = sqrt (nuA2); 
      const gdouble dmnuA_nuA = wkb->dmnuA_nuA (obj, alpha, pert->k);

      q = cexp (-I * int_nuA) / sqrt (2.0 * m * nuA);
      p = -I * cexp (-I * int_nuA) * sqrt (0.5 * m * nuA) - 0.5 * dmnuA_nuA * q;
    }
    else
    {
      const gdouble lnF       = ncm_spline_eval (wkb->lnF, alpha);;
      const gdouble dlnF      = ncm_spline_eval (wkb->dlnF, alpha);
      const gdouble one_sqrt2 = 1.0 / sqrt (2.0);

      q = cexp (-I * int_nuA - lnF * 0.5) * one_sqrt2;
      p = -I * cexp (-I * int_nuA + lnF * 0.5) * one_sqrt2 - 0.5 * m * dlnF * q;
    }
    
    *Re_q = creal (q);
    *Im_q = cimag (q);

    *Re_p = creal (p);
    *Im_p = cimag (p);    
  }
}

static gdouble 
_nc_hipert_wkb_prec (gdouble alpha, gpointer userdata)
{
  NcHIPertWKBArg *arg = (NcHIPertWKBArg *) userdata;
  const gdouble nuA2 = arg->wkb->nuA2 (arg->obj, alpha, NC_HIPERT (arg->wkb)->k);
  const gdouble V = arg->wkb->V (arg->obj, alpha, NC_HIPERT (arg->wkb)->k);
  const gdouble test = log (fabs (nuA2 * arg->prec / V));

  return test;
}

static gdouble 
_nc_hipert_wkb_prec_alpha2 (gdouble alpha, gpointer userdata)
{
  NcHIPertWKBArg *arg = (NcHIPertWKBArg *) userdata;
  const gdouble nuA2 = arg->wkb->nuA2 (arg->obj, alpha, NC_HIPERT (arg->wkb)->k);
  const gdouble test = log (fabs (nuA2 * arg->prec / (alpha * alpha)));

  return test;
}

static gdouble 
_nc_hipert_wkb_nuA2 (gdouble alpha, gpointer userdata)
{
  NcHIPertWKBArg *arg = (NcHIPertWKBArg *) userdata;
  const gdouble nuA2 = arg->wkb->nuA2 (arg->obj, alpha, NC_HIPERT (arg->wkb)->k);
  gdouble m = 0.0, nu2 = 0.0, dlnm = 0.0;
  arg->wkb->eom (arg->obj, alpha, NC_HIPERT (arg->wkb)->k, &nu2, &m, &dlnm);

  return nuA2 / nu2 - NC_HIPERT (arg->wkb)->reltol;
}

static gint
_nc_hipert_wkb_nu2_root (NcHIPertWKB *wkb, GObject *obj, gdouble alpha0, gdouble *alpha, gsl_function *F)
{
  NcHIPert *pert = NC_HIPERT (wkb);
  gint status;
  gint iter = 0, max_iter = 1000000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gdouble prec = pert->reltol, alpha1 = *alpha;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, F, alpha0, alpha1);

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

    status = gsl_root_test_residual (GSL_FN_EVAL (F, *alpha), prec);
    
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
 * nc_hipert_wkb_maxtime:
 * @wkb: a #NcHIPertWKB.
 * @obj: the WKB implementation.
 * @alpha0: the initial log-redshift time.
 * @alpha1: the final log-redshift time.
 * 
 * Search for the root of $\nu_A^2$ between $\alpha_0$ and $\alpha_1$. 
 * 
 * Returns: the root of $\nu_A^2$ between $\alpha_0$ and $\alpha_1$ or NaN if not found.
 */
gdouble 
nc_hipert_wkb_maxtime (NcHIPertWKB *wkb, GObject *obj, gdouble alpha0, gdouble alpha1)
{
  if (wkb->obj != obj)
  {
    g_clear_object (&wkb->obj);
    g_assert (g_type_is_a (G_OBJECT_TYPE (obj), wkb->impl_type));
    wkb->obj = g_object_ref (obj);
  }
  {
    gdouble alpha = alpha1;
    gsl_function F;
    NcHIPertWKBArg arg;
    guint ret;

    arg.obj = obj;
    arg.wkb = wkb;

    F.function = &_nc_hipert_wkb_nuA2;
    F.params   = &arg;
    
    ret = _nc_hipert_wkb_nu2_root (wkb, obj, alpha0, &alpha, &F);

    if (ret == 0)
      return alpha;
    else
      return GSL_NAN;
  }
}

/**
 * nc_hipert_wkb_maxtime_prec:
 * @wkb: a #NcHIPertWKB.
 * @obj: the WKB implementation.
 * @cmp: Comparison type.
 * @prec: Required precision.
 * @alpha0: the initial log-redshift time.
 * @alpha1: the final log-redshift time.
 * 
 * Search for the instant at which the WKB approximation starts to fails within the asked precision. 
 * 
 * Returns: the instant $\alpha$ between $\alpha_0$ and $\alpha_1$ or NaN if not found.
 */
gdouble 
nc_hipert_wkb_maxtime_prec (NcHIPertWKB *wkb, GObject *obj, NcHIPertWKBCmp cmp, gdouble prec, gdouble alpha0, gdouble alpha1)
{
  if (wkb->obj != obj)
  {
    g_clear_object (&wkb->obj);
    g_assert (g_type_is_a (G_OBJECT_TYPE (obj), wkb->impl_type));
    wkb->obj = g_object_ref (obj);
  }
  {  
    gdouble alpha = alpha1;
    gsl_function F;
    NcHIPertWKBArg arg;
    guint ret;

    arg.obj  = obj;
    arg.wkb  = wkb;
    arg.prec = prec;

    F.function = &_nc_hipert_wkb_nuA2;
    F.params   = &arg;

    ret = _nc_hipert_wkb_nu2_root (wkb, obj, alpha0, &alpha, &F);

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
      ret = _nc_hipert_wkb_nu2_root (wkb, obj, alpha0, &alpha, &F);
      if (ret == 0)
        return alpha;
      else
        return GSL_NAN;
    }
    else
      return GSL_NAN;
  }
}

/**
 * nc_hipert_wkb_nuA:
 * @wkb: a #NcHIPertWKB.
 * @obj: the WKB implementation.
 * @alpha: the log-redshift time.
 * 
 * FIXME 
 * 
 * Returns: FIXME
 */
gdouble 
nc_hipert_wkb_nuA (NcHIPertWKB *wkb, GObject *obj, gdouble alpha)
{
  g_assert (NC_HIPERT (wkb)->prepared);
  return ncm_spline_eval (wkb->nuA, alpha);
}

/**
 * nc_hipert_wkb_phase:
 * @wkb: a #NcHIPertWKB.
 * @obj: the WKB implementation.
 * @alpha: the log-redshift time.
 * 
 * FIXME 
 * 
 * Returns: FIXME
 */
gdouble 
nc_hipert_wkb_phase (NcHIPertWKB *wkb, GObject *obj, gdouble alpha)
{
  g_assert (NC_HIPERT (wkb)->prepared);
  if (alpha == wkb->alpha_phase)
    return wkb->cur_phase;
  else
  {
    gdouble alpha_i = wkb->alpha_phase;
    gdouble alpha_f = alpha;
    gdouble sign = 1.0;
    gdouble Dphase = 0.0;
    gdouble dphase, alpha_k;
    guint cur, tar, k;

    if (alpha_i > alpha_f)
    {
      gdouble tmp = alpha_i;
      alpha_i = alpha_f;
      alpha_f = tmp;
      sign = -1.0;
    }

    cur = ncm_spline_get_index (wkb->nuA, alpha_i);
    tar = ncm_spline_get_index (wkb->nuA, alpha_f);

    g_assert_cmpuint (cur, <=, tar);

    alpha_k = alpha_i;
    for (k = cur; k < tar; k++)
    {
      gdouble alpha_kp1 = ncm_vector_get (wkb->nuA->xv, k + 1);
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
