/***************************************************************************
 *            ncm_sf_sbessel_int.c
 *
 *  Wed Mar 10 17:16:19 2010
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
 * SECTION:ncm_sf_sbessel_int
 * @title: Spherical Bessel Integral -- Double Precision
 * @short_description: Double precision spherical bessel integrals implementation
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_sf_sbessel_int.h"
#include "math/ncm_mpsf_sbessel_int.h"
#include "math/ncm_spline_gsl.h"
#include "math/dividedifference.h"
#include "math/ncm_cfg.h"

#include <gsl/gsl_math.h>
#include <mpfr.h>

/**
 * ncm_sf_sbessel_jl_xj_integral_recur_new: (skip)
 * @jlrec: a #NcSFSBesselRecur
 * @x_grid: a #NcmGrid
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcSFSphericalBesselIntegRecur *
ncm_sf_sbessel_jl_xj_integral_recur_new (NcSFSBesselRecur *jlrec, NcmGrid *x_grid)
{
  NcSFSphericalBesselIntegRecur *xnjlrec = g_slice_new (NcSFSphericalBesselIntegRecur);
  gint i;

  g_assert (((jlrec != NULL) && (x_grid == NULL)) || ((jlrec == NULL) && (x_grid != NULL)));

  if (jlrec == NULL)
    xnjlrec->jlrec = ncm_sf_sbessel_recur_new (x_grid);
  else
    xnjlrec->jlrec = jlrec;

  for (i = 0; i < 4; i++)
  {
    xnjlrec->int_jl_xn[i] = g_slice_alloc (xnjlrec->jlrec->x_grid->nnodes * sizeof(gdouble));
    xnjlrec->int_jlp1_xn[i] = g_slice_alloc (xnjlrec->jlrec->x_grid->nnodes * sizeof(gdouble));
  }

  xnjlrec->prepared = FALSE;

  return xnjlrec;
}

/**
 * ncm_sf_sbessel_jl_xj_integral_recur_new_from_section: (skip)
 * @x_sec: a #NcmGridSection
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcSFSphericalBesselIntegRecur *
ncm_sf_sbessel_jl_xj_integral_recur_new_from_section (NcmGridSection *x_sec)
{
  NcmGrid *x_grid = ncm_grid_new_from_sections (x_sec);
  return ncm_sf_sbessel_jl_xj_integral_recur_new (NULL, x_grid);
}

/**
 * ncm_sf_sbessel_jl_xj_integral_recur_cached_new: (skip)
 * @l: FIXME
 * @x_sec: a #NcmGridSection
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcSFSphericalBesselIntegRecur *
ncm_sf_sbessel_jl_xj_integral_recur_cached_new (glong l, NcmGridSection *x_sec)
{
  NcSFSphericalBesselIntegRecur *xnjlrec;
  gchar *name_x = ncm_grid_get_name (x_sec);
  gchar *name_x_hash = g_compute_checksum_for_string (G_CHECKSUM_MD5, name_x, strlen (name_x));
  gchar *filename = g_strdup_printf ("xnjl_rule_double_%ld_%s.dat", l, name_x_hash);
  //printf ("# looking for cache (%s) name (%s)\n", filename, name_x);

  if (ncm_cfg_exists (filename))
    xnjlrec = ncm_sf_sbessel_jl_xj_integral_recur_load (filename);
  else
  {
    xnjlrec = ncm_sf_sbessel_jl_xj_integral_recur_new_from_section (x_sec);
    ncm_sf_sbessel_jl_xj_integral_recur_set (xnjlrec, l);
    ncm_sf_sbessel_jl_xj_integral_recur_save (xnjlrec, filename);
  }
  g_free (name_x);
  g_free (name_x_hash);
  g_free (filename);

  return xnjlrec;
}

/**
 * ncm_sf_sbessel_jl_xj_integral_recur_set:
 * @xnjlrec: a #NcSFSphericalBesselIntegRecur
 * @l: FIXME
 *
 * FIXME
 *
*/
void
ncm_sf_sbessel_jl_xj_integral_recur_set (NcSFSphericalBesselIntegRecur *xnjlrec, glong l)
{
  gint i;
  NcmGrid *x_grid = xnjlrec->jlrec->x_grid;

  ncm_sf_sbessel_recur_set (xnjlrec->jlrec, l);

  for (i = 0; i < 4; i++)
  {
    gint j;
    for (j = 0; j < x_grid->nnodes; j++)
    {
      gdouble x = ncm_grid_get_node_d (x_grid, j);
      xnjlrec->int_jl_xn[i][j]   = ncm_sf_sbessel_jl_xj_integral (l + 0, i, x);
      xnjlrec->int_jlp1_xn[i][j] = ncm_sf_sbessel_jl_xj_integral (l + 1, i, x);
      printf ("%.15g %.15g %.15g %ld %d\n", x, xnjlrec->int_jl_xn[i][j], xnjlrec->int_jlp1_xn[i][j], l, i);
    }
    printf ("\n\n");
  }
  xnjlrec->prepared = TRUE;
}

/**
 * ncm_sf_sbessel_jl_xj_integral_recur_free:
 * @xnjlrec: a #NcSFSphericalBesselIntegRecur
 * @free_grid: FIXME
 *
 * FIXME
 *
*/
void
ncm_sf_sbessel_jl_xj_integral_recur_free (NcSFSphericalBesselIntegRecur *xnjlrec, gboolean free_grid)
{
  ncm_sf_sbessel_recur_free (xnjlrec->jlrec, free_grid);
  g_slice_free (NcSFSphericalBesselIntegRecur, xnjlrec);
}

glong
ncm_sf_sbessel_jl_xj_integral_recur_next (NcSFSphericalBesselIntegRecur *xnjlrec)
{
  gint i, j;
  NcSFSBesselRecur *jlrec = xnjlrec->jlrec;
  NcmGrid *x_grid = jlrec->x_grid;
  const glong l = jlrec->l;

  for (j = 0; j < x_grid->nnodes; j++)
  {
    gdouble x = ncm_grid_get_node_d (x_grid, j);
    gdouble xn = 1.0;

    if (x != 0)
    {
      const gdouble temp = jlrec->jlp1[j] * (2.0 * jlrec->l + 3.0) / x - jlrec->jl[j];
      jlrec->jl[j] = jlrec->jlp1[j];
      jlrec->jl[j] = temp;
    }

    for (i = 0; i < 4; i++)
    {
      const gdouble temp1 = xn * xnjlrec->jlrec->jl[j] * (2.0 * l + 3.0);
      const gdouble temp2 = xnjlrec->int_jl_xn[i][j] * (i + l + 1.0);

      xnjlrec->int_jl_xn[i][j] = xnjlrec->int_jlp1_xn[i][j];
      xnjlrec->int_jlp1_xn[i][j] = (temp1 - temp2) / (i - l - 2.0);

      xn *= x;
    }
  }

  jlrec->l++;
  return jlrec->l;
}

/**
 * ncm_sf_sbessel_jl_xj_integral_recur_previous:
 * @xnjlrec: a #NcSFSphericalBesselIntegRecur
 *
 * FIXME
 *
 * Returns: FIXME
*/
glong
ncm_sf_sbessel_jl_xj_integral_recur_previous (NcSFSphericalBesselIntegRecur *xnjlrec)
{
  gint i, j;
  NcSFSBesselRecur *jlrec = xnjlrec->jlrec;
  NcmGrid *x_grid = jlrec->x_grid;
  const glong l = jlrec->l;

  for (j = 0; j < x_grid->nnodes; j++)
  {
    gdouble x = ncm_grid_get_node_d (x_grid, j);
    gdouble xn = 1.0;

    if (x != 0)
    {
      const gdouble temp = jlrec->jl[j] * (2.0 * jlrec->l + 1.0) / x - jlrec->jlp1[j];
      jlrec->jlp1[j] = jlrec->jl[j];
      jlrec->jl[j] = temp;
    }

    for (i = 0; i < 4; i++)
    {
      const gdouble temp1 = xn * jlrec->jlp1[j] * (2.0 * l + 1.0);
      const gdouble temp2 = xnjlrec->int_jlp1_xn[i][j] * (i - l - 1.0);
      xnjlrec->int_jlp1_xn[i][j] = xnjlrec->int_jl_xn[i][j];
      xnjlrec->int_jl_xn[i][j] = (temp1 - temp2) / (1.0 * (i + l));

      xn *= x;
    }
  }

  jlrec->l--;
  return jlrec->l;
}

/**
 * ncm_sf_sbessel_jl_xj_integral_recur_goto:
 * @xnjlrec: a #NcSFSphericalBesselIntegRecur
 * @l: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
glong
ncm_sf_sbessel_jl_xj_integral_recur_goto (NcSFSphericalBesselIntegRecur *xnjlrec, glong l)
{
  glong sign = GSL_SIGN (l - xnjlrec->jlrec->l);
  glong sub = labs(l - xnjlrec->jlrec->l);
  glong i;

  if (sub != 0)
    printf ("# Moving %ld positions %s\n", sub, sign == -1 ? "backward" : "forward");

  if (sub == 0)
    return 0;
  if (sign == 1)
    for (i = 0; i < sub; i++)
      ncm_sf_sbessel_jl_xj_integral_recur_next (xnjlrec);
  else
    for (i = 0; i < sub; i++)
      ncm_sf_sbessel_jl_xj_integral_recur_previous (xnjlrec);
  return sub;
}

/**
 * ncm_sf_sbessel_jl_xj_integral_recur_write: (skip)
 * @xnjlrec: a #NcSFSphericalBesselIntegRecur
 * @f: FIXME
 *
 * FIXME
 *
*/
void
ncm_sf_sbessel_jl_xj_integral_recur_write (NcSFSphericalBesselIntegRecur *xnjlrec, FILE *f)
{
  gint i, j;
  NcmGrid *x_grid = xnjlrec->jlrec->x_grid;

  ncm_sf_sbessel_recur_write (xnjlrec->jlrec, f);

  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < x_grid->nnodes; j++)
      NCM_WRITE_DOUBLE (f, xnjlrec->int_jl_xn[i][j]);
    for (j = 0; j < x_grid->nnodes; j++)
      NCM_WRITE_DOUBLE (f, xnjlrec->int_jlp1_xn[i][j]);
  }
}

/**
 * ncm_sf_sbessel_jl_xj_integral_recur_save:
 * @xnjlrec: a #NcSFSphericalBesselIntegRecur
 * @filename: FIXME
 *
 * FIXME
 *
*/
void
ncm_sf_sbessel_jl_xj_integral_recur_save (NcSFSphericalBesselIntegRecur *xnjlrec, gchar *filename)
{
  FILE *f = ncm_cfg_fopen (filename, "w");
  ncm_sf_sbessel_jl_xj_integral_recur_write (xnjlrec, f);
  fclose (f);
}

/**
 * ncm_sf_sbessel_jl_xj_integral_recur_read: (skip)
 * @f: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcSFSphericalBesselIntegRecur *
ncm_sf_sbessel_jl_xj_integral_recur_read (FILE *f)
{
  gint i, j;
  NcSFSBesselRecur *jlrec = ncm_sf_sbessel_recur_read (f);
  NcSFSphericalBesselIntegRecur *xnjlrec = ncm_sf_sbessel_jl_xj_integral_recur_new (jlrec, NULL);
  NcmGrid *x_grid = xnjlrec->jlrec->x_grid;

  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < x_grid->nnodes; j++)
      NCM_READ_DOUBLE (f, xnjlrec->int_jl_xn[i][j]);
    for (j = 0; j < x_grid->nnodes; j++)
      NCM_READ_DOUBLE (f, xnjlrec->int_jlp1_xn[i][j]);
  }
  xnjlrec->prepared = TRUE;

  return xnjlrec;
}

/**
 * ncm_sf_sbessel_jl_xj_integral_recur_load: (skip)
 * @filename: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcSFSphericalBesselIntegRecur *
ncm_sf_sbessel_jl_xj_integral_recur_load (gchar *filename)
{
  FILE *f = ncm_cfg_fopen (filename, "r");
  NcSFSphericalBesselIntegRecur *xnjlrec = ncm_sf_sbessel_jl_xj_integral_recur_read (f);
  fclose (f);
  return xnjlrec;
}

/**
 * ncm_sf_sbessel_jl_xj_integral:
 * @l: FIXME
 * @j: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
ncm_sf_sbessel_jl_xj_integral (gint l, gint j, gdouble x)
{
  MPFR_DECL_INIT (res, 53); /* Should it be 53? FIXME */
  gdouble res_d;
  ncm_mpsf_sbessel_jl_xj_integral (l, j, x, res, GMP_RNDN);
  res_d = mpfr_get_d (res, GMP_RNDN);
  return res_d;
}

/*********************************************************************************************************
 *
 *
 *********************************************************************************************************/

/**
 * ncm_sf_sbessel_jl_xj_integrate_spline_new: (skip)
 * @xnjlrec: a #NcSFSphericalBesselIntegRecur
 * @init: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcSFSphericalBesselIntSpline *
ncm_sf_sbessel_jl_xj_integrate_spline_new (NcSFSphericalBesselIntegRecur *xnjlrec, gboolean init)
{
  NcSFSphericalBesselIntSpline *int_jlspline = g_slice_new (NcSFSphericalBesselIntSpline);
  gint i;

  int_jlspline->xnjlrec = xnjlrec;
  int_jlspline->x_grid = xnjlrec->jlrec->x_grid;

  int_jlspline->jl = NULL;
  int_jlspline->jlp1 = NULL;

  for (i = 0; i < 4; i++)
    int_jlspline->int_jl_xn[i] = NULL;

  int_jlspline->prepared = FALSE;

  if (init)
    ncm_sf_sbessel_jl_xj_integrate_spline_reset (int_jlspline);

  return int_jlspline;
}

/**
 * ncm_sf_sbessel_jl_xj_integrate_spline_cached_new: (skip)
 * @l: FIXME
 * @x_sec: a #NcmGridSection
 * @init: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
NcSFSphericalBesselIntSpline *
ncm_sf_sbessel_jl_xj_integrate_spline_cached_new (glong l, NcmGridSection *x_sec, gboolean init)
{
  static NcSFSphericalBesselIntSpline *int_jlspline[3000];
  static gboolean tonull = TRUE;
  if (tonull)
  {
    memset (int_jlspline, 0, sizeof (NcSFSphericalBesselIntSpline *) * 3000);
    tonull = FALSE;
  }
  g_assert (l < 3000);

  if (int_jlspline[l] == NULL)
  {
    NcSFSphericalBesselIntegRecur *xnjlrec = ncm_sf_sbessel_jl_xj_integral_recur_cached_new (l, x_sec);
    int_jlspline[l] = ncm_sf_sbessel_jl_xj_integrate_spline_new (xnjlrec, init);
  }

  return int_jlspline[l];
}

/**
 * ncm_sf_sbessel_jl_xj_integral_recur_taylor_coeff:
 * @xnjlrec: a #NcSFSphericalBesselIntegRecur
 * @n: FIXME
 * @djl: FIXME
 * @djlp1: FIXME
 * @dint_jl_x0: FIXME
 * @dint_jl_x1: FIXME
 * @dint_jl_x2: FIXME
 * @dint_jl_x3: FIXME
 *
 * FIXME
 *
*/
void
ncm_sf_sbessel_jl_xj_integral_recur_taylor_coeff (NcSFSphericalBesselIntegRecur *xnjlrec,
                                                          guint n, gdouble *djl, gdouble *djlp1,
                                                          gdouble *dint_jl_x0, gdouble *dint_jl_x1,
                                                          gdouble *dint_jl_x2, gdouble *dint_jl_x3)
{
  const gdouble x = xnjlrec->jlrec->x_grid->data[n];
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;
  ncm_sf_sbessel_taylor_coeff_jl_jlp1 (xnjlrec->jlrec, n, djl, djlp1);
  {
    dint_jl_x0[0] = xnjlrec->int_jl_xn[0][n];
    dint_jl_x0[1] = djl[0];
    dint_jl_x0[2] = djl[1] / 2.0;
    dint_jl_x0[3] = djl[2] / 3.0;

    dint_jl_x1[0] = xnjlrec->int_jl_xn[1][n];
    dint_jl_x1[1] = x * djl[0];
    dint_jl_x1[2] = (djl[0] + x * djl[1]) / 2.0;
    dint_jl_x1[3] = (djl[1] + x * djl[2]) / 3.0;

    dint_jl_x2[0] = xnjlrec->int_jl_xn[2][n];
    dint_jl_x2[1] = x2 * djl[0];
    dint_jl_x2[2] = (2.0 * x * djl[0] + x2 * djl[1]) / 2.0;
    dint_jl_x2[3] = (djl[0] + 2.0 * x * djl[1] + x2 * djl[2]) / 3.0;

    dint_jl_x3[0] = xnjlrec->int_jl_xn[3][n];
    dint_jl_x3[1] = x3 * djl[0];
    dint_jl_x3[2] = (3.0 * x2 * djl[0] + x3 * djl[1]) / 2.0;
    dint_jl_x3[3] = (3.0 * x * djl[0] + 3.0 * x2 * djl[1] + x3 * djl[2]) / 3.0;
  }
}

/**
 * ncm_sf_sbessel_jl_xj_integrate_spline_reset:
 * @int_jlspline: a #NcSFSphericalBesselIntSpline
 *
 * FIXME
*/
void
ncm_sf_sbessel_jl_xj_integrate_spline_reset (NcSFSphericalBesselIntSpline *int_jlspline)
{
  gint i;
  gulong n_dd = int_jlspline->x_grid->nnodes - 1;
	NcmSpline *s = ncm_spline_gsl_new (gsl_interp_cspline);

  g_assert (int_jlspline->xnjlrec->prepared);

  if (!int_jlspline->prepared)
  {
    int_jlspline->dd_jl = g_new (gdouble, n_dd * 8);
    int_jlspline->dd_jlp1 = g_new (gdouble, n_dd * 8);
    int_jlspline->dd_int_jl_x0 = g_new (gdouble, n_dd * 8);
    int_jlspline->dd_int_jl_x1 = g_new (gdouble, n_dd * 8);
    int_jlspline->dd_int_jl_x2 = g_new (gdouble, n_dd * 8);
    int_jlspline->dd_int_jl_x3 = g_new (gdouble, n_dd * 8);
    int_jlspline->x_data = ncm_grid_get_double_array (int_jlspline->x_grid);

    int_jlspline->jl =
      ncm_spline_new_data (s,
          int_jlspline->x_data, int_jlspline->xnjlrec->jlrec->jl,
          int_jlspline->x_grid->nnodes, FALSE);

    int_jlspline->jlp1 =
      ncm_spline_new_data (s,
          int_jlspline->x_data, int_jlspline->xnjlrec->jlrec->jlp1,
          int_jlspline->x_grid->nnodes, FALSE);

    for (i = 0; i < 4; i++)
    {
      int_jlspline->int_jl_xn[i] =
        ncm_spline_new_data (s,
            int_jlspline->x_data, int_jlspline->xnjlrec->int_jl_xn[i],
            int_jlspline->x_grid->nnodes, FALSE);
    }
    int_jlspline->prepared = TRUE;
  }

  {
    ncm_grid_get_double_array (int_jlspline->xnjlrec->jlrec->x_grid);
    ncm_sf_sbessel_jl_xj_integral_recur_taylor_coeff (int_jlspline->xnjlrec, 0,
                                                              &int_jlspline->dd_jl[0], &int_jlspline->dd_jlp1[0],
                                                              &int_jlspline->dd_int_jl_x0[0], &int_jlspline->dd_int_jl_x1[0],
                                                              &int_jlspline->dd_int_jl_x2[0], &int_jlspline->dd_int_jl_x3[0]);

    for (i = 0; i < n_dd; i++)
    {
      const gsize rs = 4 * sizeof(gdouble);
      const guint ni = (2 * i + 1) * 4;
      const guint ci = (2 * i + 2) * 4;

      ncm_sf_sbessel_jl_xj_integral_recur_taylor_coeff (int_jlspline->xnjlrec, i + 1,
                                                                &int_jlspline->dd_jl[ni], &int_jlspline->dd_jlp1[ni],
                                                                &int_jlspline->dd_int_jl_x0[ni], &int_jlspline->dd_int_jl_x1[ni],
                                                                &int_jlspline->dd_int_jl_x2[ni], &int_jlspline->dd_int_jl_x3[ni]);
      memcpy (&int_jlspline->dd_jl[ci],        &int_jlspline->dd_jl[ni], rs);
      memcpy (&int_jlspline->dd_jlp1[ci],      &int_jlspline->dd_jlp1[ni], rs);
      memcpy (&int_jlspline->dd_int_jl_x0[ci], &int_jlspline->dd_int_jl_x0[ni], rs);
      memcpy (&int_jlspline->dd_int_jl_x1[ci], &int_jlspline->dd_int_jl_x1[ni], rs);
      memcpy (&int_jlspline->dd_int_jl_x2[ci], &int_jlspline->dd_int_jl_x2[ni], rs);
      memcpy (&int_jlspline->dd_int_jl_x3[ci], &int_jlspline->dd_int_jl_x3[ni], rs);

      nc_interp_dd_init_2_4 (&int_jlspline->x_data[i], &int_jlspline->dd_jl[i * 8]);
      nc_interp_dd_init_2_4 (&int_jlspline->x_data[i], &int_jlspline->dd_jlp1[i * 8]);
      nc_interp_dd_init_2_4 (&int_jlspline->x_data[i], &int_jlspline->dd_int_jl_x0[i * 8]);
      nc_interp_dd_init_2_4 (&int_jlspline->x_data[i], &int_jlspline->dd_int_jl_x1[i * 8]);
      nc_interp_dd_init_2_4 (&int_jlspline->x_data[i], &int_jlspline->dd_int_jl_x2[i * 8]);
      nc_interp_dd_init_2_4 (&int_jlspline->x_data[i], &int_jlspline->dd_int_jl_x3[i * 8]);
    }
  }
	ncm_spline_free (s);
}

/**
 * ncm_sf_sbessel_jl_xj_integrate_spline_set:
 * @int_jlspline: a #NcSFSphericalBesselIntSpline
 * @l: FIXME
 *
 * FIXME
*/
void
ncm_sf_sbessel_jl_xj_integrate_spline_set (NcSFSphericalBesselIntSpline *int_jlspline, glong l)
{
  ncm_sf_sbessel_jl_xj_integral_recur_set (int_jlspline->xnjlrec, l);
  ncm_sf_sbessel_jl_xj_integrate_spline_reset (int_jlspline);
}

/**
 * ncm_sf_sbessel_jl_xj_integrate_spline_next:
 * @int_jlspline: a #NcSFSphericalBesselIntSpline
 *
 * FIXME
*/
void
ncm_sf_sbessel_jl_xj_integrate_spline_next (NcSFSphericalBesselIntSpline *int_jlspline)
{
  ncm_sf_sbessel_jl_xj_integral_recur_next (int_jlspline->xnjlrec);
  ncm_sf_sbessel_jl_xj_integrate_spline_reset (int_jlspline);
}

/**
 * ncm_sf_sbessel_jl_xj_integrate_spline_previous:
 * @int_jlspline: a #NcSFSphericalBesselIntSpline
 *
 * FIXME
*/
void
ncm_sf_sbessel_jl_xj_integrate_spline_previous (NcSFSphericalBesselIntSpline *int_jlspline)
{
  ncm_sf_sbessel_jl_xj_integral_recur_previous (int_jlspline->xnjlrec);
  ncm_sf_sbessel_jl_xj_integrate_spline_reset (int_jlspline);
}

/**
 * ncm_sf_sbessel_jl_xj_integrate_spline_goto:
 * @int_jlspline: a #NcSFSphericalBesselIntSpline
 * @l: FIXME
 *
 * FIXME
*/
void
ncm_sf_sbessel_jl_xj_integrate_spline_goto (NcSFSphericalBesselIntSpline *int_jlspline, glong l)
{
  if (ncm_sf_sbessel_jl_xj_integral_recur_goto (int_jlspline->xnjlrec, l))
    ncm_sf_sbessel_jl_xj_integrate_spline_reset (int_jlspline);
}

inline static gdouble
_calc_xndjl (glong n, gdouble xn, gdouble jl, gdouble int_jl_xnm1)
{
  return (xn * jl - n * int_jl_xnm1);
}

inline static gdouble
_calc_xnd2jl (glong n, glong l, gdouble xnm1, gdouble xn, gdouble jl, gdouble jlp1, gdouble int_jl_xnm2)
{
  gdouble rule = 0.0;

  if (l != n)
    rule = xnm1 * jl * (l - n);

  rule -= xn * jlp1;

  if (n > 1)
    rule += int_jl_xnm2 * n * (n - 1);

  return rule;
}

inline static void
_integral_a_b_center (gdouble d, gdouble *r)
{
  const gdouble d2 = d * d;
  const gdouble d3 = d2 * d;

  r[1] -= d * r[0];
  r[2] -= 2.0 * d * r[1] + d2 * r[0];
  r[3] -= 3.0 * d * r[2] + 3.0 * d2 * r[1] + d3 * r[0];
}

/**
 * ncm_sf_sbessel_jl_xj_integral_a_b:
 * @int_jlspline: a #NcSFSphericalBesselIntSpline
 * @x0: FIXME
 * @x1: FIXME
 * @w: FIXME
 * @xnjl_rules: FIXME
 * @xndjl_rules: FIXME
 * @xnd2jl_rules: FIXME
 *
 * FIXME
*/
void
ncm_sf_sbessel_jl_xj_integral_a_b (NcSFSphericalBesselIntSpline *int_jlspline, gdouble x0, gdouble x1, gdouble w, gdouble *xnjl_rules, gdouble *xndjl_rules, gdouble *xnd2jl_rules)
{
  const glong l = int_jlspline->xnjlrec->jlrec->l;
  const gdouble x1w = x1 * w;
  const gdouble x0w = x0 * w;
  const guint nn1 = ncm_spline_get_index (int_jlspline->jl, x1w); /* gsl_interp_accel_find (int_jlspline->jl->acc, int_jlspline->jl->x, int_jlspline->jl->len, x1w); */
  const guint nn0 = ncm_spline_get_index (int_jlspline->jl, x0w); /* gsl_interp_accel_find (int_jlspline->jl->acc, int_jlspline->jl->x, int_jlspline->jl->len, x0w); */
  const gdouble jlx0   = nc_interp_dd_eval_2_4 (&int_jlspline->x_data[nn0], &int_jlspline->dd_jl[nn0 * 8], x0w);
  const gdouble jlp1x0 = nc_interp_dd_eval_2_4 (&int_jlspline->x_data[nn0], &int_jlspline->dd_jlp1[nn0 * 8], x0w);
  if (fabs(x1w - x0w) < 1e-1)
  {
    const gdouble dx = x1 - x0;
    const gdouble dx2 = dx * dx;
    const gdouble dx3 = dx2 * dx;
    const gdouble dx4 = dx2 * dx2;
    const gdouble dxw = dx * w;
    const gdouble dxw2 = dxw * dxw;
    const gdouble dxw3 = dxw2 * dxw;
    gdouble djl[6];
    ncm_sf_sbessel_deriv (l, x0w, jlx0, jlp1x0, djl);

    xnjl_rules[0]   = dx  * (djl[0] / 1.0 + djl[1] * dxw / 2.0 + djl[2] * dxw2 /  6.0 + djl[3] * dxw3 / 24.0);
    xnjl_rules[1]   = dx2 * (djl[0] / 2.0 + djl[1] * dxw / 3.0 + djl[2] * dxw2 /  8.0 + djl[3] * dxw3 / 30.0);
    xnjl_rules[2]   = dx3 * (djl[0] / 3.0 + djl[1] * dxw / 4.0 + djl[2] * dxw2 / 10.0 + djl[3] * dxw3 / 36.0);
    xnjl_rules[3]   = dx4 * (djl[0] / 4.0 + djl[1] * dxw / 5.0 + djl[2] * dxw2 / 12.0 + djl[3] * dxw3 / 42.0);

    xndjl_rules[0]  = dx  * (djl[1] / 1.0 + djl[2] * dxw / 2.0 + djl[3] * dxw2 /  6.0 + djl[4] * dxw3 / 24.0);
    xndjl_rules[1]  = dx2 * (djl[1] / 2.0 + djl[2] * dxw / 3.0 + djl[3] * dxw2 /  8.0 + djl[4] * dxw3 / 30.0);
    xndjl_rules[2]  = dx3 * (djl[1] / 3.0 + djl[2] * dxw / 4.0 + djl[3] * dxw2 / 10.0 + djl[4] * dxw3 / 36.0);
    xndjl_rules[3]  = dx4 * (djl[1] / 4.0 + djl[2] * dxw / 5.0 + djl[3] * dxw2 / 12.0 + djl[4] * dxw3 / 42.0);

    xnd2jl_rules[0] = dx  * (djl[2] / 1.0 + djl[3] * dxw / 2.0 + djl[4] * dxw2 /  6.0 + djl[5] * dxw3 / 24.0);
    xnd2jl_rules[1] = dx2 * (djl[2] / 2.0 + djl[3] * dxw / 3.0 + djl[4] * dxw2 /  8.0 + djl[5] * dxw3 / 30.0);
    xnd2jl_rules[2] = dx3 * (djl[2] / 3.0 + djl[3] * dxw / 4.0 + djl[4] * dxw2 / 10.0 + djl[5] * dxw3 / 36.0);
    xnd2jl_rules[3] = dx4 * (djl[2] / 4.0 + djl[3] * dxw / 5.0 + djl[4] * dxw2 / 12.0 + djl[5] * dxw3 / 42.0);
  }
  else
  {
    const gdouble jlx1   = nc_interp_dd_eval_2_4 (&int_jlspline->x_data[nn1], &int_jlspline->dd_jl[nn1 * 8], x1w);
    const gdouble jlp1x1 = nc_interp_dd_eval_2_4 (&int_jlspline->x_data[nn1], &int_jlspline->dd_jlp1[nn1 * 8], x1w);
    gdouble w_pow_ip1 = w;
    gdouble int_jl_x1nm1 = 0.0, int_jl_x0nm1 = 0.0;
    gdouble int_jl_x1nm2 = 0.0, int_jl_x0nm2 = 0.0;
    gdouble x0wn = 1;
    gdouble x1wn = 1;
    gdouble x0wnm1 = (x0w == 0) ? 0.0 : 1.0 / x0w;
    gdouble x1wnm1 = (x1w == 0) ? 0.0 : 1.0 / x1w;
    gdouble int_jl_x1wn[4];
    gdouble int_jl_x0wn[4];
    gint n;

    int_jl_x1wn[0] = nc_interp_dd_eval_2_4 (&int_jlspline->x_data[nn1], &int_jlspline->dd_int_jl_x0[nn1 * 8], x1w);
    int_jl_x0wn[0] = nc_interp_dd_eval_2_4 (&int_jlspline->x_data[nn0], &int_jlspline->dd_int_jl_x0[nn0 * 8], x0w);

    int_jl_x1wn[1] = nc_interp_dd_eval_2_4 (&int_jlspline->x_data[nn1], &int_jlspline->dd_int_jl_x1[nn1 * 8], x1w);
    int_jl_x0wn[1] = nc_interp_dd_eval_2_4 (&int_jlspline->x_data[nn0], &int_jlspline->dd_int_jl_x1[nn0 * 8], x0w);

    int_jl_x1wn[2] = nc_interp_dd_eval_2_4 (&int_jlspline->x_data[nn1], &int_jlspline->dd_int_jl_x2[nn1 * 8], x1w);
    int_jl_x0wn[2] = nc_interp_dd_eval_2_4 (&int_jlspline->x_data[nn0], &int_jlspline->dd_int_jl_x2[nn0 * 8], x0w);

    int_jl_x1wn[3] = nc_interp_dd_eval_2_4 (&int_jlspline->x_data[nn1], &int_jlspline->dd_int_jl_x3[nn1 * 8], x1w);
    int_jl_x0wn[3] = nc_interp_dd_eval_2_4 (&int_jlspline->x_data[nn0], &int_jlspline->dd_int_jl_x3[nn0 * 8], x0w);

    for (n = 0; n < 4; n++)
    {
      gdouble int_jl_x1n = int_jl_x1wn[n];
      gdouble int_jl_x0n = int_jl_x0wn[n];

      xnjl_rules[n]  =
        (
         int_jl_x1n -
         int_jl_x0n
         ) / w_pow_ip1;

      xndjl_rules[n] =
        (
         _calc_xndjl (n, x1wn, jlx1, int_jl_x1nm1) -
         _calc_xndjl (n, x0wn, jlx0, int_jl_x0nm1)
         ) / w_pow_ip1;

      xnd2jl_rules[n] =
        (
         _calc_xnd2jl (n, l, x1wnm1, x1wn, jlx1, jlp1x1, int_jl_x1nm2) -
         _calc_xnd2jl (n, l, x0wnm1, x0wn, jlx0, jlp1x0, int_jl_x0nm2)
         ) / w_pow_ip1;

      int_jl_x1nm2 = int_jl_x1nm1;
      int_jl_x1nm1 = int_jl_x1n;

      int_jl_x0nm2 = int_jl_x0nm1;
      int_jl_x0nm1 = int_jl_x0n;

      w_pow_ip1 *= w;
      x1wnm1 = x1wn;
      x1wn *= x1w;
      x0wnm1 = x0wn;
      x0wn *= x0w;
    }

    _integral_a_b_center (x0, xnjl_rules);
    _integral_a_b_center (x0, xndjl_rules);
    _integral_a_b_center (x0, xnd2jl_rules);
  }
}

typedef struct
{
  double * c;
  double * g;
  double * diag;
  double * offdiag;
} cspline_state_t;

static inline void
coeff_calc (const double c_array[], double dy, double dx, size_t index, double *b, double *c, double *d)
{
  const double c_i = c_array[index];
  const double c_ip1 = c_array[index + 1];
  *b = (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0;
  *c = c_i;
  *d = (c_ip1 - c_i) / (3.0 * dx);
}

inline static void
_get_spline_coeff (NcmSpline *s, gdouble *c, gint i)
{
	NcmSplineGsl *sg = NCM_SPLINE_GSL (s);
  const cspline_state_t *state = (const cspline_state_t *) (sg->interp)->state;
  const gdouble a = ncm_vector_get (s->xv, i);
  const gdouble b = ncm_vector_get (s->xv, i + 1);
  const gdouble y_lo = ncm_vector_get (s->yv, i);
  const gdouble y_hi = ncm_vector_get (s->yv, i + 1);
  const gdouble dx = b - a;
  const gdouble dy = y_hi - y_lo;

  c[0] = y_lo;
  coeff_calc(state->c, dy, dx, i,  &c[1], &c[2], &c[3]);
}

/**
 * ncm_sf_sbessel_jl_xj_integral_spline:
 * @int_jlspline: a #NcSFSphericalBesselIntSpline
 * @s0: a #NcmSpline
 * @s1: a #NcmSpline
 * @s2: a #NcmSpline
 * @w: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
ncm_sf_sbessel_jl_xj_integral_spline (NcSFSphericalBesselIntSpline *int_jlspline, NcmSpline *s0, NcmSpline *s1, NcmSpline *s2, gdouble w)
{
  gdouble res[3] = {0.0, 0.0, 0.0};
  gdouble c[4];
  gdouble r0[4];
  gdouble r1[4];
  gdouble r2[4];
  gdouble *r;
  gint i;

  for (i = s0->len - 2; i >= 0; i--)
  {
    const gdouble a = ncm_vector_get (s0->xv, i);
    const gdouble b = ncm_vector_get (s0->xv, i + 1);
    gdouble p[3];

    ncm_sf_sbessel_jl_xj_integral_a_b (int_jlspline, a, b, w, r0, r1, r2);

    r = r0;
    _get_spline_coeff (s0, c, i);
    p[0] = (c[0] * r[0]) + (c[1] * r[1]) + (c[2] * r[2]) + (c[3] * r[3]);
    res[0] += p[0];

    r = r1;
    _get_spline_coeff (s1, c, i);
    p[1] = (c[0] * r[0]) + (c[1] * r[1]) + (c[2] * r[2]) + (c[3] * r[3]);
    res[1] += p[1];

    r = r2;
    _get_spline_coeff (s2, c, i);
    p[2] = (c[0] * r[0]) + (c[1] * r[1]) + (c[2] * r[2]) + (c[3] * r[3]);
    res[2] += p[2];

    if ((fabs((p[0] + p[1] + p[2])/(res[0] + res[1] + res[2])) < GSL_DBL_EPSILON) || ((p[0] + p[1] + p[2]) == 0))
      break;
  }
//printf ("%.15g %.15g %.15g %.15g\n", w, res[0], res[1], res[2]);
  return res[0] + res[1] + res[2];
}

/**
 * ncm_sf_sbessel_jl_xj_integrate_spline_eval:
 * @int_jlspline: a #NcSFSphericalBesselIntSpline
 * @d: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
ncm_sf_sbessel_jl_xj_integrate_spline_eval (NcSFSphericalBesselIntSpline *int_jlspline, gint d, gdouble x)
{
  return ncm_spline_eval (int_jlspline->int_jl_xn[d], x);
}
