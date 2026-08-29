/***************************************************************************
 *            config_extra.h
 *
 *  Tue September 17 17:57:01 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
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

#include <math.h>

/* glibc only declares exp10() in <math.h> under __USE_GNU (_GNU_SOURCE) or
 * in C23 mode; HAVE_EXP10 (a link-level check) does not guarantee the
 * declaration is visible. NumCosmo does not require _GNU_SOURCE, so provide
 * the fallback whenever the real declaration isn't in scope. */
#if !defined (__USE_GNU) && (!defined (__STDC_VERSION__) || __STDC_VERSION__ < 202311L)
#define exp10(x) (exp ((x) * M_LN10))
#endif

/* sincos() is a glibc/BSD extension only declared under __USE_GNU (i.e. with
 * _GNU_SOURCE), independently of the C standard version and of HAVE_SINCOS
 * (a link-level check that doesn't guarantee the declaration is visible).
 * NumCosmo does not require _GNU_SOURCE, so provide the portable fallback
 * whenever the real declaration isn't in scope. */
#ifndef __USE_GNU

static inline void
sincos (double x, double *s, double *c)
{
  s[0] = sin (x);
  c[0] = cos (x);
}

#endif /* __USE_GNU */

#if HAVE_DECL_LGAMMA_R == 0
double lgamma_r (double x, int *signp);

#endif /* HAVE_DECL_LGAMMA_R == 0 */

#endif /* NUMCOSMO_GIR_SCAN */

#endif /* _CONFIG_EXTRA_H_ */

