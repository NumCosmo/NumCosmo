/***************************************************************************
 *            nc_xcor_limber_z.c
 *
 *  Sat August 29 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

/*
 * The redshift-space Limber tier: both methods behind
 * %NC_XCOR_METHOD_LIMBER_Z_GSL and %NC_XCOR_METHOD_LIMBER_Z_CUBATURE, and the
 * two #NcmIntegralND subtypes they drive.
 *
 * This is a different approximation from the rest of #NcXcor, not a different
 * quadrature of the same integral: it integrates over z with k fixed at
 * (l + 1/2) / chi(z), so it shares no closure, no knot set and no code with
 * the kernel-space methods in nc_xcor_kquad.c.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/integration/ncm_integrate.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/integration/ncm_integral_nd.h"
#include "nc/xcor/nc_xcor.h"
#include "nc/xcor/nc_xcor_priv.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */

static void nc_xcor_auto_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);
static void nc_xcor_auto_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);
static void nc_xcor_cross_dim (NcmIntegralND *intnd, guint *dim, guint *fdim);
static void nc_xcor_cross_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval);

NCM_INTEGRAL_ND_DEFINE_TYPE (NC, XCOR_AUTO, NcXcorAuto, nc_xcor_auto, nc_xcor_auto_dim, nc_xcor_auto_integ, NcXcorArg);
NCM_INTEGRAL_ND_DEFINE_TYPE (NC, XCOR_CROSS, NcXcorCross, nc_xcor_cross, nc_xcor_cross_dim, nc_xcor_cross_integ, NcXcorArg);

static void
nc_xcor_auto_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  NcXcorAuto *xcor_auto = NC_XCOR_AUTO (intnd);
  NcXcorArg *xcor_arg   = &xcor_auto->data;

  *dim  = 1;
  *fdim = xcor_arg->nells;
}

static void
nc_xcor_cross_dim (NcmIntegralND *intnd, guint *dim, guint *fdim)
{
  NcXcorCross *xcor_cross = NC_XCOR_CROSS (intnd);
  NcXcorArg *xcor_arg     = &xcor_cross->data;

  *dim  = 1;
  *fdim = xcor_arg->nells;
}

static void
nc_xcor_auto_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcXcorAuto *xcor_int    = NC_XCOR_AUTO (intnd);
  NcXcorArg *xcor_int_arg = &xcor_int->data;
  guint i, j;

  for (i = 0; i < npoints; i++)
  {
    const gdouble z         = ncm_vector_fast_get (x, i);
    const gdouble xi_z      = nc_distance_comoving (xcor_int_arg->dist, xcor_int_arg->cosmo, z); /* in units of Hubble radius */
    const gdouble xi_z_phys = xi_z * xcor_int_arg->RH;                                           /* in Mpc */
    const gdouble E_z       = nc_hicosmo_E (xcor_int_arg->cosmo, z);
    const NcXcorKinetic xck = { xi_z, E_z };

    for (j = 0; j < fdim; j++)
    {
      const gint l             = xcor_int_arg->ells[j];
      const gdouble k          = (l + 0.5) / (xi_z_phys); /* in Mpc-1 */
      const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xcor_int_arg->ps), NCM_MODEL (xcor_int_arg->cosmo), z, k);
      const gdouble k1z        = nc_xcor_kernel_eval_limber_z (xcor_int_arg->xclk1, xcor_int_arg->cosmo, z, &xck, l);
      const gdouble res        = gsl_pow_2 (k1z / xi_z) * power_spec / E_z;

      ncm_vector_fast_set (fval, i * fdim + j, res);
    }
  }
}

static void
nc_xcor_cross_integ (NcmIntegralND *intnd, NcmVector *x, guint dim, guint npoints, guint fdim, NcmVector *fval)
{
  NcXcorCross *xcor_cross = NC_XCOR_CROSS (intnd);
  NcXcorArg *xcor_arg     = &xcor_cross->data;
  guint i, j;

  for (i = 0; i < npoints; i++)
  {
    const gdouble z         = ncm_vector_fast_get (x, i);
    const gdouble xi_z      = nc_distance_comoving (xcor_arg->dist, xcor_arg->cosmo, z); /* in units of Hubble radius */
    const gdouble xi_z_phys = xi_z * xcor_arg->RH;                                       /* in Mpc */
    const gdouble E_z       = nc_hicosmo_E (xcor_arg->cosmo, z);
    const NcXcorKinetic xck = { xi_z, E_z };

    for (j = 0; j < fdim; j++)
    {
      const gint l             = xcor_arg->ells[j];
      const gdouble k          = (l + 0.5) / (xi_z_phys); /* in Mpc-1 */
      const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xcor_arg->ps), NCM_MODEL (xcor_arg->cosmo), z, k);
      const gdouble k1z        = nc_xcor_kernel_eval_limber_z (xcor_arg->xclk1, xcor_arg->cosmo, z, &xck, l);
      const gdouble k2z        = nc_xcor_kernel_eval_limber_z (xcor_arg->xclk2, xcor_arg->cosmo, z, &xck, l);
      const gdouble res        = k1z * k2z * power_spec / (xi_z * xi_z * E_z);

      ncm_vector_fast_set (fval, i * fdim + j, res);
    }
  }
}

typedef struct _xcor_gsl
{
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcmPowspec *ps;

  NcXcorKernel *xclk1;
  NcXcorKernel *xclk2;
  guint l;

  gdouble RH;
} xcor_gsl;

static gdouble
_xcor_limber_z_gsl_cross_int (gdouble z, gpointer ptr)
{
  xcor_gsl *xclki          = (xcor_gsl *) ptr;
  const gdouble xi_z       = nc_distance_comoving (xclki->dist, xclki->cosmo, z); /* in units of Hubble radius */
  const gdouble xi_z_phys  = xi_z * xclki->RH;                                    /* in Mpc */
  const gdouble E_z        = nc_hicosmo_E (xclki->cosmo, z);
  const NcXcorKinetic xck  = { xi_z, E_z };
  const gdouble k          = (xclki->l + 0.5) / (xi_z_phys); /* in Mpc-1 */
  const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xclki->ps), NCM_MODEL (xclki->cosmo), z, k);

  const gdouble k1z = nc_xcor_kernel_eval_limber_z (xclki->xclk1, xclki->cosmo, z, &xck, xclki->l);
  const gdouble k2z = nc_xcor_kernel_eval_limber_z (xclki->xclk2, xclki->cosmo, z, &xck, xclki->l);

  return k1z * k2z * power_spec / (xi_z * xi_z * E_z);
}

static gdouble
_xcor_limber_z_gsl_auto_int (gdouble z, gpointer ptr)
{
  xcor_gsl *xclki          = (xcor_gsl *) ptr;
  const gdouble xi_z       = nc_distance_comoving (xclki->dist, xclki->cosmo, z); /* in units of Hubble radius */
  const gdouble xi_z_phys  = xi_z * xclki->RH;                                    /* in Mpc */
  const gdouble E_z        = nc_hicosmo_E (xclki->cosmo, z);
  const NcXcorKinetic xck  = { xi_z, E_z };
  const gdouble k          = (xclki->l + 0.5) / (xi_z_phys); /* in Mpc-1 */
  const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xclki->ps), NCM_MODEL (xclki->cosmo), z, k);
  const gdouble k1z        = nc_xcor_kernel_eval_limber_z (xclki->xclk1, xclki->cosmo, z, &xck, xclki->l);

  return gsl_pow_2 (k1z / xi_z) * power_spec / E_z;
}

void
_nc_xcor_limber_z_gsl (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gdouble zmin, gdouble zmax, gboolean isauto, NcmVector *vp)
{
  xcor_gsl xclki;
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

  zmin = (zmin != 0.0) ? zmin : 1.0e-6;

  if (isauto)
    F.function = &_xcor_limber_z_gsl_auto_int;
  else
    F.function = &_xcor_limber_z_gsl_cross_int;

  F.params = &xclki;

  gsl_integration_workspace **w = ncm_integral_get_workspace ();

  for (i = 0; i < lmax - lmin + 1; i++)
  {
    xclki.l = lmin + i;
    /* GSL integration sometimes underestimates the error, so we multiply the relative tolerance by 1e-2 */
    ret = gsl_integration_qag (&F, zmin, zmax, 0.0, xc->reltol * 1.0e-2, NCM_INTEGRAL_PARTITION, 6, *w, &r, &err);

    _nc_xcor_check_qag_status ("_nc_xcor_limber_z_gsl", ret, xc->reltol, r, err);

    ncm_vector_set (vp, i, r);
  }

  ncm_memory_pool_return (w);
}

static void
_nc_xcor_cubature_worker (NcmIntegralND *xcor_int_nd, NcXcorArg *xcor_arg, gdouble zmin, gdouble zmax, guint lmin, guint lmax, NcmVector *vp)
{
  NcmVector *z_min   = ncm_vector_new_data_static (&zmin, 1, 1);
  NcmVector *z_max   = ncm_vector_new_data_static (&zmax, 1, 1);
  const guint size   = lmax - lmin + 1;
  NcmVector *err     = ncm_vector_new (size);
  GArray *ells_array = g_array_new (FALSE, FALSE, sizeof (gint));
  const gint block   = sqrt (size);
  guint i;

  zmin = (zmin != 0.0) ? zmin : 1.0e-6;

  ncm_integral_nd_set_reltol (xcor_int_nd, xcor_arg->xc->reltol);
  ncm_integral_nd_set_abstol (xcor_int_nd, 0.0);
  ncm_integral_nd_set_method (xcor_int_nd, NCM_INTEGRAL_ND_METHOD_CUBATURE_P_V);

  g_array_set_size (ells_array, size);

  for (i = 0; i < size; i++)
  {
    g_array_index (ells_array, gint, i) = lmin + i;
  }

  for (i = 0; i + block < size; i += block)
  {
    NcmVector *vp_i  = ncm_vector_get_subvector (vp, i, block);
    NcmVector *err_i = ncm_vector_get_subvector (err, i, block);

    xcor_arg->ells  = &g_array_index (ells_array, gint, i);
    xcor_arg->nells = block;

    ncm_integral_nd_eval (xcor_int_nd, z_min, z_max, vp_i, err_i);
    ncm_vector_free (vp_i);
    ncm_vector_free (err_i);
  }

  {
    NcmVector *vp_i  = ncm_vector_get_subvector (vp, i, size - i);
    NcmVector *err_i = ncm_vector_get_subvector (err, i, size - i);

    xcor_arg->ells  = &g_array_index (ells_array, gint, i);
    xcor_arg->nells = size - i;

    ncm_integral_nd_eval (xcor_int_nd, z_min, z_max, vp_i, err_i);
    ncm_vector_free (vp_i);
    ncm_vector_free (err_i);
  }

  g_array_unref (ells_array);
  ncm_vector_free (z_min);
  ncm_vector_free (z_max);
  ncm_vector_free (err);
}

void
_nc_xcor_limber_z_cubature (NcXcor *xc, NcXcorKernel *xclk1, NcXcorKernel *xclk2, NcHICosmo *cosmo, guint lmin, guint lmax, gdouble zmin, gdouble zmax, gboolean isauto, NcmVector *vp)
{
  if (isauto)
  {
    NcXcorAuto *xcor_int       = g_object_new (nc_xcor_auto_get_type (), NULL);
    NcXcorArg *xcor_arg        = &xcor_int->data;
    NcmIntegralND *xcor_int_nd = NCM_INTEGRAL_ND (xcor_int);

    xcor_arg->xc    = xc;
    xcor_arg->dist  = xc->dist;
    xcor_arg->ps    = xc->ps;
    xcor_arg->RH    = xc->RH;
    xcor_arg->cosmo = cosmo;
    xcor_arg->xclk1 = xclk1;
    xcor_arg->xclk2 = xclk1;

    _nc_xcor_cubature_worker (xcor_int_nd, xcor_arg, zmin, zmax, lmin, lmax, vp);

    g_object_unref (xcor_int);
  }
  else
  {
    NcXcorCross *xcor_cross    = g_object_new (nc_xcor_cross_get_type (), NULL);
    NcXcorArg *xcor_arg        = &xcor_cross->data;
    NcmIntegralND *xcor_int_nd = NCM_INTEGRAL_ND (xcor_cross);

    xcor_arg->xc    = xc;
    xcor_arg->dist  = xc->dist;
    xcor_arg->ps    = xc->ps;
    xcor_arg->RH    = xc->RH;
    xcor_arg->cosmo = cosmo;
    xcor_arg->xclk1 = xclk1;
    xcor_arg->xclk2 = xclk2;

    _nc_xcor_cubature_worker (xcor_int_nd, xcor_arg, zmin, zmax, lmin, lmax, vp);

    g_object_unref (xcor_cross);
  }
}

