/***************************************************************************
 *            ncm_sf_sbessel.c
 *
 *  Wed Mar 10 17:15:25 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:ncm_sf_sbessel
 * @title: Spherical Bessel -- Double Precision
 * @short_description: Double precision spherical bessel implementation
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sf_sbessel.h"
#include "math/ncm_mpsf_sbessel.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include <gsl/gsl_math.h>
#include <mpfr.h>

/**
 * ncm_sf_sbessel_recur_new: (skip)
 * @x_grid: a #NcmGrid
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmSFSBesselRecur *
ncm_sf_sbessel_recur_new (NcmGrid *x_grid)
{
  NcmSFSBesselRecur *jlrec = g_slice_new (NcmSFSBesselRecur);
  jlrec->prepared = FALSE;
  jlrec->x_grid = x_grid;
  jlrec->jl = g_slice_alloc (x_grid->nnodes * sizeof(gdouble));
  jlrec->jlp1 = g_slice_alloc (x_grid->nnodes * sizeof(gdouble));

  return jlrec;
}

/**
 * ncm_sf_sbessel_recur_set:
 * @jlrec: a #NcmSFSBesselRecur
 * @l: FIXME
 *
 * FIXME
*/
void
ncm_sf_sbessel_recur_set (NcmSFSBesselRecur *jlrec, glong l)
{
  guint i;
  NcmGrid *x_grid = jlrec->x_grid;
  jlrec->l = l;

  for (i = 0; i < x_grid->nnodes; i++)
  {
    gdouble x = ncm_grid_get_node_d (x_grid, i);
    jlrec->jl[i] = ncm_sf_sbessel (l, x);
    jlrec->jlp1[i] = ncm_sf_sbessel (l + 1, x);
  }
  jlrec->prepared = TRUE;
}

/**
 * ncm_sf_sbessel_recur_free:
 * @jlrec: a #NcmSFSBesselRecur
 * @free_grid: FIXME
 *
 * FIXME
*/
void
ncm_sf_sbessel_recur_free (NcmSFSBesselRecur *jlrec, gboolean free_grid)
{
  NcmGrid *x_grid = jlrec->x_grid;
  g_slice_free1 (x_grid->nnodes * sizeof(gdouble), jlrec->jl);
  g_slice_free1 (x_grid->nnodes * sizeof(gdouble), jlrec->jlp1);
  g_slice_free (NcmSFSBesselRecur, jlrec);
  if (free_grid)
    ncm_grid_free (x_grid, TRUE);
}

/**
 * ncm_sf_sbessel_recur_next:
 * @jlrec: a #NcmSFSBesselRecur
 *
 * FIXME
*/
void
ncm_sf_sbessel_recur_next (NcmSFSBesselRecur *jlrec)
{
  NcmGrid *x_grid = jlrec->x_grid;
  guint i;
  
  for (i = 0; i < x_grid->nnodes; i++)
  {
    gdouble x = ncm_grid_get_node_d (x_grid, i);
    if (x != 0)
    {
      const gdouble temp = jlrec->jlp1[i] * (2.0 * jlrec->l + 3.0) / x - jlrec->jl[i];
      jlrec->jl[i] = jlrec->jlp1[i];
      jlrec->jl[i] = temp;
    }
  }
  jlrec->l++;
}

/**
 * ncm_sf_sbessel_recur_previous:
 * @jlrec: a #NcmSFSBesselRecur
 *
 * FIXME
*/
void
ncm_sf_sbessel_recur_previous (NcmSFSBesselRecur *jlrec)
{
  NcmGrid *x_grid = jlrec->x_grid;
  guint i;

  for (i = 0; i < x_grid->nnodes; i++)
  {
    gdouble x = ncm_grid_get_node_d (x_grid, i);
    if (x != 0)
    {
      const gdouble temp = jlrec->jl[i] * (2.0 * jlrec->l + 1.0) / x - jlrec->jlp1[i];
      jlrec->jlp1[i] = jlrec->jl[i];
      jlrec->jl[i] = temp;
    }
  }
  jlrec->l--;
}

/**
 * ncm_sf_sbessel_recur_goto:
 * @jlrec: a #NcmSFSBesselRecur
 * @l: FIXME
 *
 * FIXME
*/
void
ncm_sf_sbessel_recur_goto (NcmSFSBesselRecur *jlrec, glong l)
{
  glong sign = GSL_SIGN (l - jlrec->l);
  glong sub = labs(l - jlrec->l);
  glong i;
  if (sub == 0)
    return;
  if (sign == 1)
    for (i = 0; i < sub; i++)
      ncm_sf_sbessel_recur_next (jlrec);
  else
    for (i = 0; i < sub; i++)
      ncm_sf_sbessel_recur_previous (jlrec);
}


static void
_taylor_jl (const glong l, const gdouble x, const gdouble x2, const gdouble x3, const gdouble jl, const gdouble jlp1, gdouble *deriv)
{
  const gdouble llm1 = l * (l - 1.0);
  const gdouble llm1lm2 = llm1 * (l - 2);

  deriv[0] = jl;
  deriv[1] = (l * jl - x * jlp1) / x;
  deriv[2] = (((llm1 - x2) * jl + 2.0 * x * jlp1) / x2) / (1.0 * 2.0);
  deriv[3] = (((llm1lm2 - (l - 2.0) * x2) * jl - x * (l * (l + 1.0) + 6.0 - x2) * jlp1) / x3) / (1.0 * 2.0 * 3.0);
}

/**
 * ncm_sf_sbessel_taylor_coeff_jl_jlp1:
 * @jlrec: a #NcmSFSBesselRecur
 * @n: FIXME
 * @djl: FIXME
 * @djlp1: FIXME
 *
 * FIXME
*/
void
ncm_sf_sbessel_taylor_coeff_jl_jlp1 (NcmSFSBesselRecur *jlrec, guint n, gdouble *djl, gdouble *djlp1)
{
  if (jlrec->x_grid->data[n] != 0.0)
  {
    const gdouble x  = jlrec->x_grid->data[n];
    const gdouble x2 = x * x;
    const gdouble x3 = x2 * x;
    const glong l = jlrec->l;
    const gdouble jl   = jlrec->jl[n];
    const gdouble jlp1 = jlrec->jlp1[n];
    const gdouble jlp2 = (2.0 * l + 3.0) * jlp1 / x - jl;

    _taylor_jl (l,     x, x2, x3, jl,   jlp1, djl);
    _taylor_jl (l + 1, x, x2, x3, jlp1, jlp2, djlp1);
  }
  else if (jlrec->l > 3)
    djl[0] = djlp1[0] = djl[1] = djlp1[1] = djl[2] = djlp1[2] = djl[3] = djlp1[3] = 0.0;
  else
  {
    switch (jlrec->l)
    {
      case 0:
        djl[0] = 1.0;
        djl[1] = 0.0;
        djl[2] = -1.0 / 3.0 / 2.0;
        djl[3] = 0.0;
        djlp1[0] = 0.0;
        djlp1[1] = 1.0 / 3.0;
        djlp1[2] = 0.0;
        djlp1[3] = -1.0 / 5.0 / 6.0;
        break;
      case 1:
        djl[0] = 0.0;
        djl[1] = 1.0 / 3.0;
        djl[2] = 0.0;
        djl[3] = -1.0 / 5.0 / 6.0;
        djlp1[0] = 0.0;
        djlp1[1] = 0.0;
        djlp1[2] = 2.0 / 15.0 / 2.0;
        djlp1[3] = 0.0;
        break;
      case 2:
        djl[0] = 0.0;
        djl[1] = 0.0;
        djl[2] = 2.0 / 15.0 / 2.0;
        djl[3] = 0.0;
        djlp1[0] = 0.0;
        djlp1[1] = 0.0;
        djlp1[2] = 0.0;
        djlp1[3] = 2.0 / 35.0 / 6.0;
        break;
      case 3:
        djl[0] = 0.0;
        djl[1] = 0.0;
        djl[2] = 0.0;
        djl[3] = 2.0 / 35.0 / 6.0;
        djlp1[0] = 0.0;
        djlp1[1] = 0.0;
        djlp1[2] = 0.0;
        djlp1[3] = 0.0;
        break;
    }
  }
}

/**
 * ncm_sf_sbessel_recur_write: (skip)
 * @jlrec: a #NcmSFSBesselRecur
 * @f: FIXME
 *
 * FIXME
*/
void
ncm_sf_sbessel_recur_write (NcmSFSBesselRecur *jlrec, FILE *f)
{
  guint i;
  NCM_WRITE_INT32(f, jlrec->l);
  ncm_grid_write (jlrec->x_grid, f);
  for (i = 0; i < jlrec->x_grid->nnodes; i++)
    NCM_WRITE_DOUBLE(f, jlrec->jl[i]);
  for (i = 0; i < jlrec->x_grid->nnodes; i++)
    NCM_WRITE_DOUBLE(f, jlrec->jlp1[i]);
}

/**
 * ncm_sf_sbessel_recur_read: (skip)
 * @f: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcmSFSBesselRecur *
ncm_sf_sbessel_recur_read (FILE *f)
{
  guint i;
  gint32 l;
  NcmSFSBesselRecur *jlrec;
  NcmGrid *x_grid;

  NCM_READ_INT32(f, l);
  x_grid = ncm_grid_read (f);
  jlrec = ncm_sf_sbessel_recur_new (x_grid);
  jlrec->l = l;

  for (i = 0; i < jlrec->x_grid->nnodes; i++)
    NCM_READ_DOUBLE(f, jlrec->jl[i]);
  for (i = 0; i < jlrec->x_grid->nnodes; i++)
    NCM_READ_DOUBLE(f, jlrec->jlp1[i]);

  jlrec->prepared = TRUE;

  return jlrec;
}

/**
 * ncm_sf_sbessel:
 * @l: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
ncm_sf_sbessel (gulong l, gdouble x)
{
  MPFR_DECL_INIT (res, 53); /* Should it be 53? FIXME */
  gdouble res_d;
  ncm_mpsf_sbessel_d (l, x, res, GMP_RNDN);
  res_d = mpfr_get_d (res, GMP_RNDN);
  return res_d;
}

/**
 * ncm_sf_sbessel_taylor:
 * @l: FIXME
 * @x: FIXME
 * @djl: FIXME
 *
 * FIXME
*/
void
ncm_sf_sbessel_taylor (gulong l, gdouble x, gdouble *djl)
{
  const gdouble jl = ncm_sf_sbessel (l, x);
  const gdouble jlp1 = ncm_sf_sbessel (l + 1, x);
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;

  _taylor_jl (l, x, x2, x3, jl, jlp1, djl);
  return;
}

/**
 * ncm_sf_sbessel_deriv:
 * @l: FIXME
 * @x: FIXME
 * @jl: FIXME
 * @jlp1: FIXME
 * @djl: FIXME
 *
 * FIXME
*/
void
ncm_sf_sbessel_deriv (gulong l, gdouble x, gdouble jl, gdouble jlp1, gdouble *djl)
{
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;
  const gdouble x4 = x2 * x2;
  const gdouble x5 = x3 * x2;
  const gdouble llm1 = l * (l - 1.0);
  const gdouble llm1lm2 = llm1 * (l - 2);
  const gdouble llm1lm2lm3 = llm1lm2 * (l - 3);
  const gdouble llm1lm2lm3lm4 = llm1lm2lm3 * (l - 4);

  if (x != 0.0)
  {
    djl[0] = jl;
    djl[1] = (l * jl - x * jlp1) / x;
    djl[2] = ((llm1 - x2) * jl + 2.0 * x * jlp1) / x2;
    djl[3] = ((llm1lm2 - (l - 2.0) * x2) * jl - x * (l * (l + 1.0) + 6.0 - x2) * jlp1) / x3;
    djl[4] = ((llm1lm2lm3 - 2.0 * (4.0 + llm1) * x2 + x4) * jl + (4.0 * x * (6.0 + 2.0 * l * (l + 1.0) - x2)) * jlp1) / x4;
    djl[5] = ((llm1lm2lm3lm4 - 2.0 * (l * (2.0 + (l - 7.0) * l) - 20.0) * x2 + (l - 4.0) * x4) * jl - (x * (l * (l + 1) * (58.0 + l * (1.0 + l)) - 2.0 * l * (l + 1.0) * x2 + x4 - 20.0 * (x2 - 6.0))) * jlp1) / x5;
  }
  else if (l > 5)
    djl[0] = djl[1] = djl[2] = djl[3] = djl[4] = djl[5] = 0.0;
  else
  {
    switch (l)
    {
      case 0:
        djl[0] = 1.0;
        djl[1] = 0.0;
        djl[2] = -1.0 / 3.0;
        djl[3] = 0.0;
        djl[4] = 1.0 / 5.0;
        djl[5] = 0.0;
        break;
      case 1:
        djl[0] = 0.0;
        djl[1] = 1.0 / 3.0;
        djl[2] = 0.0;
        djl[3] = -1.0 / 5.0;
        djl[4] = 0.0;
        djl[5] = 1.0 / 7.0;
        break;
      case 2:
        djl[0] = 0.0;
        djl[1] = 0.0;
        djl[2] = 2.0 / 15.0;
        djl[3] = 0.0;
        djl[4] = -4.0 / 35.0;
        djl[5] = 0.0;
        break;
      case 3:
        djl[0] = 0.0;
        djl[1] = 0.0;
        djl[2] = 0.0;
        djl[3] = 2.0 / 35.0;
        djl[4] = 0.0;
        djl[5] = -4.0 / 63.0;
        break;
      case 4:
        djl[0] = 0.0;
        djl[1] = 0.0;
        djl[2] = 0.0;
        djl[3] = 0.0;
        djl[4] = 8.0 / 315.0;
        djl[5] = 0.0;
        break;
      case 5:
        djl[0] = 0.0;
        djl[1] = 0.0;
        djl[2] = 0.0;
        djl[3] = 0.0;
        djl[4] = 0.0;
        djl[5] = 8.0 / 693.0;
        break;
    }
  }
  return;
}

static gdouble
_ncm_sf_sbessel_spline_calc (gdouble x, gpointer data)
{
	gulong *l = (gulong *)data;
	return ncm_sf_sbessel (*l, x);
}

/**
 * ncm_sf_sbessel_spline:
 * @l: FIXME
 * @xi: FIXME
 * @xf: FIXME
 * @reltol: FIXME
 *
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmSpline *
ncm_sf_sbessel_spline (gulong l, gdouble xi, gdouble xf, gdouble reltol)
{
	NcmSpline *s = ncm_spline_cubic_notaknot_new ();
	gsl_function F;

	F.function = &_ncm_sf_sbessel_spline_calc;
	F.params = &l;

	ncm_spline_set_func (s, NCM_SPLINE_FUNCTION_SPLINE, &F, xi, xf, 0, reltol);
	return s;
}
