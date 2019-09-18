/***************************************************************************
 *            config_extra.h
 *
 *  Tue September 17 17:57:01 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * config_extra.h
 *
 * Copyright (C) 2019 - Sandro Dias Pinto Vitenti
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _CONFIG_EXTRA_H_
#define _CONFIG_EXTRA_H_

#ifndef NUMCOSMO_GIR_SCAN

#ifndef HAVE_EXP10
#define exp10(x) (exp ((x) * M_LN10))
#endif /* HAVE_EXP10 */

#ifndef HAVE_SINCOS
#include <math.h>

static inline void 
sincos (double x, double *s, double *c)
{
  s[0] = sin (x);
  c[0] = cos (x);
}
#endif /* HAVE_SINCOS */

#if HAVE_DECL_LGAMMA_R == 0
double lgamma_r(double x, int *signp);
#endif /* HAVE_DECL_LGAMMA_R == 0 */

#endif /* NUMCOSMO_GIR_SCAN */

#endif /* _CONFIG_EXTRA_H_ */
