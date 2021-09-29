/***************************************************************************
 *            nc_xcor.h
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
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
 * SECTION:nc_xcor
 * @title: Cross-correlations
 * @short_description: Angular auto- and cross-spectra.
 *
 * The angular power spectrum between observables $A$ and $B$ with kernels $W^A$ and $W^B$ is
 * \begin{equation}
 * C_{\ell}^{AB} &= \int dz W^A(z) \int dz^\prime W^B (z^\prime) \times \int dk \frac{2}{\pi} k^2 P(k, z, z^\prime) j_{\ell}(k\chi(z)) j_{\ell} (k\chi(z^\prime)).
 * \end{equation}
 * In the Limber approximation, it reduces to
 * \begin{equation}
 * C_{\ell}^{AB} = \int_0^{z_*} dz \frac{H(z)}{c \chi^2(z)} W^A(z) W^B (z) P\left(k = \frac{\ell +1/2}{\chi(z)} , z \right),
 * \end{equation}
 * where $P\left(k = \frac{\ell +1/2}{\chi(z)} , z \right)$ is the power spectrum (a #NcmPowspec) at redshift $z$ and $chi(z)$ the comoving distance (a #NcDistance).
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/integral.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_cfg.h"
#include "math/ncm_serialize.h"
#include "xcor/nc_xcor.h"
#include "nc_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <cuba.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_DISTANCE,
  PROP_MATTER_POWER_SPECTRUM,
  PROP_METH,
};

G_DEFINE_TYPE (NcXcor, nc_xcor, G_TYPE_OBJECT);
G_DEFINE_BOXED_TYPE (NcXcorKinetic, nc_xcor_kinetic, nc_xcor_kinetic_copy, nc_xcor_kinetic_free);

static void
nc_xcor_init (NcXcor *xc)
{
  xc->ps   = NULL;
  xc->dist = NULL;
  xc->RH   = 0.0;
  xc->meth = NC_XCOR_LIMBER_METHOD_GSL;
}

static void
_nc_xcor_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcXcor *xc = NC_XCOR (object);
  
  g_return_if_fail (NC_IS_XCOR (object));
  
  switch (prop_id)
  {
    case PROP_DISTANCE:
      nc_distance_clear (&xc->dist);
      xc->dist = g_value_dup_object (value);
      break;
    case PROP_MATTER_POWER_SPECTRUM:
      ncm_powspec_clear (&xc->ps);
      xc->ps = g_value_dup_object (value);
      break;
    case PROP_METH:
      xc->meth = g_value_get_enum (value);
      break;
      
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_xcor_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcXcor *xc = NC_XCOR (object);
  
  g_return_if_fail (NC_IS_XCOR (object));
  
  switch (prop_id)
  {
    case PROP_DISTANCE:
      g_value_set_object (value, xc->dist);
      break;
    case PROP_MATTER_POWER_SPECTRUM:
      g_value_set_object (value, xc->ps);
      break;
    case PROP_METH:
      g_value_set_enum (value, xc->meth);
      break;
      
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_xcor_dispose (GObject *object)
{
  NcXcor *xc = NC_XCOR (object);
  
  nc_distance_clear (&xc->dist);
  ncm_powspec_clear (&xc->ps);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_parent_class)->dispose (object);
}

static void
_nc_xcor_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_xcor_parent_class)->finalize (object);
}

static void
nc_xcor_class_init (NcXcorClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  
  /*GObjectClass* parent_class = G_OBJECT_CLASS (klass); */
  
  object_class->set_property = &_nc_xcor_set_property;
  object_class->get_property = &_nc_xcor_get_property;
  object_class->dispose      = &_nc_xcor_dispose;
  object_class->finalize     = &_nc_xcor_finalize;
  
  /**
   * NcXcor:distance:
   *
   * This property keeps the distance object.
   */
  g_object_class_install_property (object_class,
                                   PROP_DISTANCE,
                                   g_param_spec_object ("distance",
                                                        NULL,
                                                        "Distance.",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcXcor:power-spec:
   *
   * This property keeps the matter power spectrum object.
   */
  g_object_class_install_property (object_class,
                                   PROP_MATTER_POWER_SPECTRUM,
                                   g_param_spec_object ("power-spec",
                                                        NULL,
                                                        "Matter power spectrum.",
                                                        NCM_TYPE_POWSPEC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcXcor:meth:
   *
   * This property keeps the method enumerator to compute the integrals.
   */
  g_object_class_install_property (object_class,
                                   PROP_METH,
                                   g_param_spec_enum ("meth",
                                                      NULL,
                                                      "Method.",
                                                      NC_TYPE_XCOR_LIMBER_METHOD,
                                                      NC_XCOR_LIMBER_METHOD_GSL,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_xcor_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 * @meth: a #NcXcorLimberMethod to compute the Limber integrals
 *
 * Two methods are available to compute Limber-approximated integrals: independent GSL numerical integration or vector integration using Sundials's CVode algorithm.
 *
 * Returns: FIXME
 *
 */
NcXcor *
nc_xcor_new (NcDistance *dist, NcmPowspec *ps, NcXcorLimberMethod meth)
{
  return g_object_new (NC_TYPE_XCOR,
                       "distance", dist,
                       "power-spec", ps,
                       "meth", meth,
                       NULL);
}

/**
 * nc_xcor_ref:
 * @xc: a #NcXcor
 *
 * Returns: (transfer full): @xc
 */
NcXcor *
nc_xcor_ref (NcXcor *xc)
{
  return g_object_ref (xc);
}

/**
 * nc_xcor_free:
 * @xc: a #NcXcor
 *
 * FIXME
 *
 */
void
nc_xcor_free (NcXcor *xc)
{
  g_object_unref (xc);
}

/**
 * nc_xcor_clear:
 * @xc: a #NcXcor
 *
 * FIXME
 *
 */
void
nc_xcor_clear (NcXcor **xc)
{
  g_clear_object (xc);
}

/**
 * nc_xcor_kinetic_copy:
 * @xck: a #NcXcorKinetic
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcXcorKinetic *
nc_xcor_kinetic_copy (NcXcorKinetic *xck)
{
  NcXcorKinetic *xck_copy = g_new (NcXcorKinetic, 1);
  
  xck_copy[0] = xck[0];
  
  return xck_copy;
}

/**
 * nc_xcor_kinetic_free:
 * @xck: a #NcXcorKinetic
 *
 * FIXME
 *
 */
void
nc_xcor_kinetic_free (NcXcorKinetic *xck)
{
  g_free (xck);
}

/**
 * nc_xcor_prepare:
 * @xc: a #NcXcor
 * @cosmo: a #NcHICosmo
 *
 * FIXME
 *
 */
void
nc_xcor_prepare (NcXcor *xc, NcHICosmo *cosmo)
{
  nc_distance_prepare_if_needed (xc->dist, cosmo);
  ncm_powspec_prepare_if_needed (xc->ps, NCM_MODEL (cosmo));
  
  xc->RH = nc_hicosmo_RH_Mpc (cosmo);
}

typedef struct _xcor_limber_cvode
{
  gboolean isauto;
  guint lmin, lmax;
  guint nell;
  NcmVector *k;
  NcmVector *Pk;
  NcXcor *xc;
  NcXcorLimberKernel *xclk1;
  NcXcorLimberKernel *xclk2;
  NcHICosmo *cosmo;
} xcor_limber_cvode;

static gint
_xcor_limber_cvode_int (realtype z, N_Vector y, N_Vector ydot, gpointer params)
{
  xcor_limber_cvode *xclc = (xcor_limber_cvode *) params;
  const gdouble xi_z      = nc_distance_comoving (xclc->xc->dist, xclc->cosmo, z); /* in units of Hubble radius */
  const gdouble xi_z_phys = xi_z * xclc->xc->RH;                                   /* in Mpc */
  const gdouble E_z       = nc_hicosmo_E (xclc->cosmo, z);
  const NcXcorKinetic xck = { xi_z, E_z };
  const gdouble k1z       = nc_xcor_limber_kernel_eval (xclc->xclk1, xclc->cosmo, z, &xck, 0);
  
  gdouble geoW1W2;
  guint i, l;
  
  if (G_UNLIKELY (z == 0.0))
  {
    N_VConst (0.0, ydot);
    
    return 0;
  }
  
  if (xclc->isauto)
  {
    geoW1W2 = E_z * gsl_pow_2 (k1z / xi_z);
  }
  else
  {
    const gdouble k2z = nc_xcor_limber_kernel_eval (xclc->xclk2, xclc->cosmo, z, &xck, 0);
    
    geoW1W2 = E_z * k1z * k2z / (xi_z * xi_z);
  }
  
  for (i = 0; i < xclc->nell; i++)
  {
    l = i + xclc->lmin;
    ncm_vector_fast_set (xclc->k, i, ((gdouble) l + 0.5) / xi_z_phys); /* in Mpc-1 */
  }
  
  ncm_powspec_eval_vec (xclc->xc->ps, NCM_MODEL (xclc->cosmo), z, xclc->k, xclc->Pk);
  
  for (i = 0; i < xclc->nell; i++)
  {
    NV_Ith_S (ydot, i) = ncm_vector_fast_get (xclc->Pk, i) * geoW1W2;
  }
  
  return 0;
}

static void
_nc_xcor_limber_cvode (NcXcor *xc, NcXcorLimberKernel *xclk1, NcXcorLimberKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gdouble zmin, gdouble zmax, gboolean isauto, NcmVector *vp)
{
  const guint nell = lmax - lmin + 1;
  const gdouble dz = (zmax - zmin) * NC_XCOR_PRECISION;
  
  N_Vector yv        = N_VNew_Serial (nell);
  N_Vector yv0       = N_VNew_Serial (nell);
  SUNMatrix A        = SUNDenseMatrix (nell, nell);
  SUNLinearSolver LS = SUNDenseLinearSolver (yv, A);
  gpointer cvode     = CVodeCreate (CV_ADAMS);
  gpointer cvodefunc = &_xcor_limber_cvode_int;
  NcmVector *Pk      = ncm_vector_new (nell);
  NcmVector *k       = ncm_vector_new (nell);
  gdouble z;
  gint flag;
  guint i;
  
  ncm_vector_set_zero (k);
  ncm_vector_set_zero (Pk);
  
  for (i = 0; i < nell; i++)
  {
    NV_Ith_S (yv, i)  = 0.0;
    NV_Ith_S (yv0, i) = 0.0;
  }
  
  xcor_limber_cvode xclc = { isauto, lmin, lmax, nell, k, Pk, xc, xclk1, xclk2, cosmo };
  
  /* Find initial value = y'(zmid) * dz */
  const gdouble zmid = expm1 ((log1p (xclk1->zmid) + log1p (xclk2->zmid)) / 2.0);
  
  _xcor_limber_cvode_int (zmid, yv, yv0, &xclc);
  N_VScale (dz, yv0, yv0);
  
  /* Init and set flags */
  flag = CVodeInit (cvode, cvodefunc, zmid, yv0);
  NCM_CVODE_CHECK (&flag, "CVodeInit", 1, );
  
  flag = CVodeSStolerances (cvode, NC_XCOR_PRECISION, 0.0);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );
  flag = CVodeSetMaxNumSteps (cvode, G_MAXUINT32);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );
  flag = CVodeSetUserData (cvode, &xclc);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );
  flag = CVodeSetLinearSolver (cvode, LS, A);
  NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );
  
  /* First integrate from zmid to zmin */
  flag = CVode (cvode, zmin, yv, &z, CV_NORMAL);
  NCM_CVODE_CHECK (&flag, "CVode", 1, );
  
  /* Put int[zmin->zmid] in yv */
  for (i = 0; i < nell; i++)
  {
    NV_Ith_S (yv, i) = -(NV_Ith_S (yv, i) - NV_Ith_S (yv0, i));
  }
  
  /* Then integrate from zmid to zmax */
  flag = CVodeReInit (cvode, zmid, yv);
  NCM_CVODE_CHECK (&flag, "CVodeReInit", 1, );
  flag = CVodeSStolerances (cvode, NC_XCOR_PRECISION, 0.0);
  NCM_CVODE_CHECK (&flag, "CVodeSStolerances", 1, );
  flag = CVodeSetMaxNumSteps (cvode, G_MAXUINT32);
  NCM_CVODE_CHECK (&flag, "CVodeSetMaxNumSteps", 1, );
  flag = CVodeSetUserData (cvode, &xclc);
  NCM_CVODE_CHECK (&flag, "CVodeSetUserData", 1, );
  flag = CVodeSetLinearSolver (cvode, LS, A);
  NCM_CVODE_CHECK (&flag, "CVodeSetLinearSolver", 1, );
  
  flag = CVode (cvode, zmax, yv, &z, CV_NORMAL);
  NCM_CVODE_CHECK (&flag, "CVode", 1, );
  
  /* Finally put the result in vp */
  for (i = 0; i < nell; i++)
  {
    ncm_vector_set (vp, i, NV_Ith_S (yv, i));
  }
  
  CVodeFree (&cvode);
  N_VDestroy (yv);
  N_VDestroy (yv0);
  SUNMatDestroy (A);
  SUNLinSolFree (LS);
  
  ncm_vector_free (k);
  ncm_vector_free (Pk);
}

#ifdef HAVE_LIBCUBA

/* typedef int (*integrand_t)(const int *ndim, const double x[], const int *ncomp, double f[], void *userdata); */

typedef struct _xcor_limber_suave
{
  gboolean isauto;
  guint lmin, lmax;
  guint nell;
  NcmVector *k;
  NcmVector *Pk;
  NcXcor *xc;
  NcXcorLimberKernel *xclk1;
  NcXcorLimberKernel *xclk2;
  NcHICosmo *cosmo;
  gdouble zmin;
  gdouble zmax_zmin;
  guint counter;
} xcor_limber_suave;

gint
_nc_xcor_limber_suave_integrand (const gint *ndim, const gdouble x[], const gint *ncomp, gdouble f[], gpointer userdata)
{
  xcor_limber_suave *xcls = (xcor_limber_suave *) userdata;
  
  const gdouble z = xcls->zmin + x[0] * xcls->zmax_zmin; /* x goes from 0 to 1 */
  
  /* printf("counter = %15i, z = %g\n", xcls->counter, z); */
  
  const gdouble xi_z      = nc_distance_comoving (xcls->xc->dist, xcls->cosmo, z); /* in units of Hubble radius */
  const gdouble xi_z_phys = xi_z * xcls->xc->RH;                                   /* in Mpc */
  const gdouble E_z       = nc_hicosmo_E (xcls->cosmo, z);
  const NcXcorKinetic xck = { xi_z, E_z };
  const gdouble k1z       = nc_xcor_limber_kernel_eval (xcls->xclk1, xcls->cosmo, z, &xck, 0);
  
  gdouble geoW1W2;
  guint i, l;
  
  if (G_UNLIKELY (z == 0.0))
  {
    for (i = 0; i < xcls->nell; i++)
    {
      f[i] = 0.0;
    }
    
    return 0; /*supposes f is initialized to zero */
  }
  
  if (xcls->isauto)
  {
    geoW1W2 = E_z * gsl_pow_2 (k1z / xi_z);
  }
  else
  {
    const gdouble k2z = nc_xcor_limber_kernel_eval (xcls->xclk2, xcls->cosmo, z, &xck, 0);
    
    geoW1W2 = E_z * k1z * k2z / (xi_z * xi_z);
  }
  
  for (i = 0; i < xcls->nell; i++)
  {
    l = i + xcls->lmin;
    ncm_vector_fast_set (xcls->k, i, ((gdouble) l + 0.5) / xi_z_phys); /* in Mpc-1 */
  }
  
  ncm_powspec_eval_vec (xcls->xc->ps, NCM_MODEL (xcls->cosmo), z, xcls->k, xcls->Pk);
  
  for (i = 0; i < xcls->nell; i++)
  {
    f[i] = ncm_vector_fast_get (xcls->Pk, i) * geoW1W2;
  }
  
  xcls->counter += 1;
  
  return 0;
}

static void
_nc_xcor_limber_suave (NcXcor *xc, NcXcorLimberKernel *xclk1, NcXcorLimberKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gdouble zmin, gdouble zmax, gboolean isauto, NcmVector *vp)
{
  const guint nell = lmax - lmin + 1;
  NcmVector *Pk    = ncm_vector_new (nell);
  NcmVector *k     = ncm_vector_new (nell);
  
  ncm_vector_set_zero (k);
  ncm_vector_set_zero (Pk);
  
  guint counter = 0;
  
  xcor_limber_suave xcls = { isauto, lmin, lmax, nell, k, Pk, xc, xclk1, xclk2, cosmo, zmin, zmax - zmin, counter };
  
  /*  in */
  const gdouble flatness = 100.0;
  const guint mineval    = 1;
  const guint maxeval    = 1000000;
  const guint nnew       = 1000;
  const guint nmin       = 500;
  const gint flags       = 0;
  
  /* out */
  gint nregions, neval, fail;
  
  gdouble integral[nell];
  gdouble error[nell];
  gdouble prob[nell];
  
  Suave (1, nell,
         &_nc_xcor_limber_suave_integrand, &xcls, 1,
         NC_XCOR_PRECISION, 0.0,
         flags, 0,
         mineval, maxeval,
         nnew, nmin,
         flatness, "", NULL,
         &nregions, &neval, &fail,
         integral, error, prob);
  
  /* printf("nregions = %i, neval = %i, fail = %i \n", nregions, neval, fail ); */
  
  guint i;
  
  for (i = 0; i < nell; i++)
  {
    ncm_vector_set (vp, i, integral[i]);
    /* printf("%i, integral = %14.10g, error = %14.10g, prob = %14.10g \n", i, integral[i], error[i], prob[i]); */
  }
  
  /* suave integrates over the hypercube, so scaling is necessary */
  ncm_vector_scale (vp, zmax - zmin);
}

#endif /* HAVE_LIBCUBA */


typedef struct _xcor_limber_gsl
{
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcmPowspec *ps;
  
  NcXcorLimberKernel *xclk1;
  NcXcorLimberKernel *xclk2;
  guint l;
  
  gdouble RH;
} xcor_limber_gsl;

static gdouble
_xcor_limber_gsl_cross_int (gdouble z, gpointer ptr)
{
  xcor_limber_gsl *xclki   = (xcor_limber_gsl *) ptr;
  const gdouble xi_z       = nc_distance_comoving (xclki->dist, xclki->cosmo, z); /* in units of Hubble radius */
  const gdouble xi_z_phys  = xi_z * xclki->RH;                                    /* in Mpc */
  const gdouble E_z        = nc_hicosmo_E (xclki->cosmo, z);
  const NcXcorKinetic xck  = { xi_z, E_z };
  const gdouble k          = (xclki->l + 0.5) / (xi_z_phys); /* in Mpc-1 */
  const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xclki->ps), NCM_MODEL (xclki->cosmo), z, k);
  
  const gdouble k1z = nc_xcor_limber_kernel_eval (xclki->xclk1, xclki->cosmo, z, &xck, xclki->l);
  const gdouble k2z = nc_xcor_limber_kernel_eval (xclki->xclk2, xclki->cosmo, z, &xck, xclki->l);
  
  return E_z * k1z * k2z * power_spec / (xi_z * xi_z);
}

static gdouble
_xcor_limber_gsl_auto_int (gdouble z, gpointer ptr)
{
  xcor_limber_gsl *xclki   = (xcor_limber_gsl *) ptr;
  const gdouble xi_z       = nc_distance_comoving (xclki->dist, xclki->cosmo, z); /* in units of Hubble radius */
  const gdouble xi_z_phys  = xi_z * xclki->RH;                                    /* in Mpc */
  const gdouble E_z        = nc_hicosmo_E (xclki->cosmo, z);
  const NcXcorKinetic xck  = { xi_z, E_z };
  const gdouble k          = (xclki->l + 0.5) / (xi_z_phys); /* in Mpc-1 */
  const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xclki->ps), NCM_MODEL (xclki->cosmo), z, k);
  const gdouble k1z        = nc_xcor_limber_kernel_eval (xclki->xclk1, xclki->cosmo, z, &xck, xclki->l);
  
  return E_z * gsl_pow_2 (k1z / xi_z) * power_spec;
}

static void
_nc_xcor_limber_gsl (NcXcor *xc, NcXcorLimberKernel *xclk1, NcXcorLimberKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gdouble zmin, gdouble zmax, gboolean isauto, NcmVector *vp)
{
  xcor_limber_gsl xclki;
  gdouble r, err;
  gsl_function F;
  guint i;
  gint ret;
  
  xclki.xclk1 = xclk1;
  xclki.xclk2 = xclk2;
  xclki.cosmo = cosmo;
  xclki.dist  = xc->dist;
  xclki.ps    = xc->ps;
  xclki.RH    = xc->RH;
  
  if (isauto)
    F.function = &_xcor_limber_gsl_auto_int;
  else
    F.function = &_xcor_limber_gsl_cross_int;

  F.params = &xclki;
  
  gsl_integration_workspace **w = ncm_integral_get_workspace ();
  
  for (i = 0; i < lmax - lmin + 1; i++)
  {
    xclki.l = lmin + i;
    ret     = gsl_integration_qag (&F, zmin, zmax, 0.0, NC_XCOR_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &r, &err);
    
    if (ret != GSL_SUCCESS)
      g_error ("_nc_xcor_limber_gsl: %s.", gsl_strerror (ret));
    
    ncm_vector_set (vp, i, r);
  }
  
  ncm_memory_pool_return (w);
}

/**
 * nc_xcor_limber:
 * @xc: a #NcXcor
 * @xclk1: a #NcXcorLimberKernel
 * @xclk2: a #NcXcorLimberKernel
 * @cosmo: a #NcHICosmo
 * @lmin: a #guint
 * @lmax: a #guint
 * @vp: a #NcmVector
 *
 * Performs the computation of the power spectrum $C_{\ell}^{AB}$ in the Limber approximation. The kernels of observables A and B are @xclk1 and @xclk2. If @xclk2 is NULL, the auto power spectrum is computed. The result for multipoles lmin to lmax (included) is stored in the #NcmVector @vp.
 *
 */
void
nc_xcor_limber (NcXcor *xc, NcXcorLimberKernel *xclk1, NcXcorLimberKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, NcmVector *vp)
{
  const guint nell = ncm_vector_len (vp);
  const gboolean isauto = (xclk2 == xclk1);
  const gdouble cons_factor = ((isauto) ? gsl_pow_2 (xclk1->cons_factor) : xclk1->cons_factor * xclk2->cons_factor) / gsl_pow_3 (xc->RH);
  gdouble zmin, zmax;
  
  if (nell != lmax - lmin + 1)
    g_error ("nc_xcor_limber: vector size does not match multipole limits");
  
  if (isauto)
  {
    zmin = xclk1->zmin;
    zmax = xclk1->zmax;
  }
  else
  {
    zmin = GSL_MAX (xclk1->zmin, xclk2->zmin);
    zmax = GSL_MIN (xclk1->zmax, xclk2->zmax);
  }
  
  if (zmin < zmax)
  {
    switch (xc->meth)
    {
      case NC_XCOR_LIMBER_METHOD_CVODE:
        _nc_xcor_limber_cvode (xc, xclk1, xclk2, cosmo, lmin, lmax, zmin, zmax, isauto, vp);
        break;
      case NC_XCOR_LIMBER_METHOD_GSL:
        _nc_xcor_limber_gsl (xc, xclk1, xclk2, cosmo, lmin, lmax, zmin, zmax, isauto, vp);
        break;
#ifdef HAVE_LIBCUBA
      case NC_XCOR_LIMBER_METHOD_SUAVE:
        _nc_xcor_limber_suave (xc, xclk1, xclk2, cosmo, lmin, lmax, zmin, zmax, isauto, vp);
        break;
#endif /* HAVE_LIBCUBA */
      default:
        g_assert_not_reached ();
        break;
    }
    
    ncm_vector_scale (vp, cons_factor);
  }
  else
  {
    ncm_vector_set_zero (vp);
  }
}

