/***************************************************************************
 *            ncm_complex.h
 *
 *  Thu September 25 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_complex.h
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

#ifndef _NCM_COMPLEX_H_
#define _NCM_COMPLEX_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#ifndef NUMCOSMO_GIR_SCAN
typedef complex double NcmComplex;
#else /* NUMCOSMO_GIR_SCAN */
typedef struct _NcmComplexShouldNeverAppear NcmComplex;
#endif /* NUMCOSMO_GIR_SCAN */

GType ncm_complex_get_type (void) G_GNUC_CONST;

NcmComplex *ncm_complex_new (void);
NcmComplex *ncm_complex_dup (NcmComplex *c);
void ncm_complex_free (NcmComplex *c);
void ncm_complex_clear (NcmComplex **c);

NCM_INLINE void ncm_complex_set (NcmComplex *c, const gdouble a, const gdouble b);
NCM_INLINE void ncm_complex_set_zero (NcmComplex *c);

NCM_INLINE gdouble ncm_complex_Re (const NcmComplex *c);
NCM_INLINE gdouble ncm_complex_Im (const NcmComplex *c);
NCM_INLINE gdouble ncm_complex_Abs (const NcmComplex *c);

#ifndef NUMCOSMO_GIR_SCAN
NCM_INLINE void ncm_complex_set_c (NcmComplex *c, const complex double z);
NCM_INLINE complex double ncm_complex_c (const NcmComplex *c);

#endif /* NUMCOSMO_GIR_SCAN */

NCM_INLINE void ncm_complex_res_add_mul_real (NcmComplex * restrict c1, const NcmComplex * restrict c2, const gdouble v);
NCM_INLINE void ncm_complex_res_add_mul (NcmComplex * restrict c1, const NcmComplex * restrict c2, const NcmComplex * restrict c3);

NCM_INLINE void ncm_complex_mul_real (NcmComplex *c, const gdouble v);
NCM_INLINE void ncm_complex_res_mul (NcmComplex * restrict c1, const NcmComplex * restrict c2);

#define NCM_COMPLEX_ZERO (0.0)
#define NCM_COMPLEX(p) ((NcmComplex *) (p))
#define NCM_COMPLEX_PTR(p) ((NcmComplex **) (p))
#define NCM_COMPLEX_INIT(z) (z)
#define NCM_COMPLEX_INIT_REAL(z) (z)

G_END_DECLS

#endif /* _NCM_COMPLEX_H_ */

#ifndef _NCM_COMPLEX_INLINE_H_
#define _NCM_COMPLEX_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE void
ncm_complex_set (NcmComplex *c, const gdouble a, const gdouble b)
{
  *c = a + I * b;
}

NCM_INLINE void
ncm_complex_set_zero (NcmComplex *c)
{
  *c = 0.0;
}

NCM_INLINE gdouble
ncm_complex_Re (const NcmComplex *c)
{
  return creal (*c);
}

NCM_INLINE gdouble
ncm_complex_Im (const NcmComplex *c)
{
  return cimag (*c);
}

NCM_INLINE gdouble
ncm_complex_Abs (const NcmComplex *c)
{
  return cabs (*c);
}

#ifndef NUMCOSMO_GIR_SCAN

NCM_INLINE void
ncm_complex_set_c (NcmComplex *c, const complex double z)
{
  *c = z;
}

NCM_INLINE complex double
ncm_complex_c (const NcmComplex *c)
{
  return *c;
}

#endif /* NUMCOSMO_GIR_SCAN */

NCM_INLINE void
ncm_complex_res_add_mul_real (NcmComplex * restrict c1, const NcmComplex * restrict c2, const gdouble v)
{
  *c1 += (*c2) * v;
}

NCM_INLINE void
ncm_complex_res_add_mul (NcmComplex * restrict c1, const NcmComplex * restrict c2, const NcmComplex * restrict c3)
{
  *c1 += (*c2) * (*c3);
}

NCM_INLINE void
ncm_complex_mul_real (NcmComplex *c, const gdouble v)
{
  *c *= v;
}

NCM_INLINE void
ncm_complex_res_mul (NcmComplex * restrict c1, const NcmComplex * restrict c2)
{
  *c1 *= *c2;
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_COMPLEX_INLINE_H_ */

