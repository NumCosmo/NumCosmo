/***************************************************************************
 *            dividedifference.c
 *
 *  Wed Mar 17 16:20:30 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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
 * SECTION:dividedifference
 * @title: Divided Difference
 * @short_description: Divided difference methods for function interpolation with derivatives
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <string.h>
#include <stdio.h>
#include <glib.h>
#include <gmp.h>
#include <mpfr.h>

static void
_interp_dd_direct (const gdouble *vx, const gint nt, const gint nf, const gint n, register gint xbi, register gint xbf, gdouble *dd)
{
  register gint i;
  register gdouble last = dd[n];
  register gint xai = 0;
  register gint xaf = 0;
  const gint next_xbi = xbi;
  const gint next_xbf = xbf;

  for (i = n + 1; i < nt;)
  {
    const gint k = i++;
    const gdouble new_dd_k = (dd[k] - last) / (vx[xbi] - vx[xai]);
    last = dd[k];
    dd[k] = new_dd_k;
    xaf++; xbf++;
    if (xaf >= nf) { xaf = 0; xai++; }
    if (xbf >= nf) { xbf = 0; xbi++; }
  }

  if (n < nt)
  {
    if (next_xbf + 1 >= nf)
      _interp_dd_direct (vx, nt, nf, n + 1, next_xbi + 1, 0, dd);
    else
      _interp_dd_direct (vx, nt, nf, n + 1, next_xbi, next_xbf + 1, dd);
  }
  
}

static void
_interp_dd (const gdouble *vx, const gint np, const gint nf, const gint n, gdouble *dd)
{
  register gint i;
  register gdouble last = dd[n];
  
  for (i = 0; i < np - 1; i++)
  {
    register gint j;
    const gdouble dxi = vx[i + 1] - vx[i];
    
    for (j = 0; j <= n; j++)
    {
      const gint k = (i + 1) * nf + j;
      const gdouble new_dd_k = (dd[k] - last) / dxi;
      last = dd[k];
      dd[k] = new_dd_k;
    }
  }

  if (n + 2 >= nf)
    _interp_dd_direct (vx, np * nf, nf, n + 1, 1, 0, dd);
  else
    _interp_dd (vx, np, nf, n + 1, dd);

}

/***********************************************************************************************************/

static void
_interp_dd_2_4 (const gdouble *vx, gdouble *dd)
{
  const gdouble one_dxi1  = 1.0 / (vx[1] - vx[0]);
  const gdouble one_dxi2  = one_dxi1 * one_dxi1;
  const gdouble one_dxi3  = one_dxi2 * one_dxi1;
  const gdouble one_dxi4  = one_dxi2 * one_dxi2;
  const gdouble one_dxi5  = one_dxi3 * one_dxi2;
  const gdouble one_dxi6  = one_dxi3 * one_dxi3;
  const gdouble one_dxi7  = one_dxi4 * one_dxi3;

  dd[7] = 
    20.0 * one_dxi7 * dd[0] +
    10.0 * one_dxi6 * dd[1] +
     4.0 * one_dxi5 * dd[2] +
     1.0 * one_dxi4 * dd[3] -
    20.0 * one_dxi7 * dd[4] + 
    10.0 * one_dxi6 * dd[5] -
     4.0 * one_dxi5 * dd[6] +
     1.0 * one_dxi4 * dd[7];

  dd[6] = 
    -(
      10.0 * one_dxi6 * dd[0] + 
       6.0 * one_dxi5 * dd[1] + 
       3.0 * one_dxi4 * dd[2] +
       1.0 * one_dxi3 * dd[3]
      ) +
    10.0 * one_dxi6 * dd[4] -
     4.0 * one_dxi5 * dd[5] +
     1.0 * one_dxi4 * dd[6];

  dd[5] = 
    4.0 * one_dxi5 * dd[0] + 
    3.0 * one_dxi4 * dd[1] + 
    2.0 * one_dxi3 * dd[2] + 
    1.0 * one_dxi2 * dd[3] -
    4.0 * one_dxi5 * dd[4] + 
    1.0 * one_dxi4 * dd[5];

  dd[4] = 
    -(
      one_dxi4 * dd[0] + 
      one_dxi3 * dd[1] + 
      one_dxi2 * dd[2] + 
      one_dxi1 * dd[3]
      )+
    one_dxi4 * dd[4];
}

/************************************************************************************************************/

void
nc_interp_dd_init (const gdouble *vx, gdouble *dd, const gint np, const gint nf)
{
  g_assert (nf >= 1);
  
  if (nf == 1)
    _interp_dd_direct (vx, np, nf, 0, 1, 0, dd);
  else if (nf == 4 && np == 2 && TRUE)
    _interp_dd_2_4 (vx, dd);
  else
    _interp_dd (vx, np, nf, 0, dd);
  
  return;
}

void
nc_interp_dd_init_2_4 (const gdouble *vx, gdouble *dd)
{
  _interp_dd_2_4 (vx, dd);
}

/************************************************************************************************************/

static gdouble
_interp_dd_eval (const gdouble *vx, const gdouble *dd, const gdouble x, const gint np, const gint nf, 
    register gint i, register gint j, register gint k)
{
  if (j < nf - 1)
    return dd[i] + (x - vx[k]) * _interp_dd_eval (vx, dd, x, np, nf, i + 1, j + 1, k);
  else if (k < np - 1)
    return dd[i] + (x - vx[k]) * _interp_dd_eval (vx, dd, x, np, nf, i + 1, 0, k + 1);
  else
    return dd[i];
}

gdouble 
nc_interp_dd_eval (const gdouble *vx, const gdouble *dd, const gdouble x, const gint np, const gint nf)
{
  return _interp_dd_eval (vx, dd, x, np, nf, 0, 0, 0);
}

gdouble 
nc_interp_dd_eval_2_4 (const gdouble *vx, const gdouble *dd, const gdouble x)
{
  gdouble dx0 = (x - vx[0]);
  gdouble dx1 = (x - vx[1]);
  return dd[0] + dx0 * (dd[1] + dx0 * (dd[2] + dx0 * (dd[3] + dx0 * (dd[4] + dx1 * (dd[5] + dx1 * (dd[6] + dx1 * dd[7]))))));
}
