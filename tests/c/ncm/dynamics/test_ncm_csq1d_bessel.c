/***************************************************************************
 *            test_ncm_csq1d_bessel.c
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_csq1d_bessel.c
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>
#include <math.h>

#include "test_ncm_csq1d_bessel.h"

/* m = (s t)^(1 + 2 a), xi = ln (m nu), F1 = xi' / (2 nu), F2 = F1' / (2 nu). */

gdouble
test_csq1d_bessel_m (gdouble a, gdouble k, gdouble s, gdouble t)
{
  return pow (s * t, 1.0 + 2.0 * a);
}

gdouble
test_csq1d_bessel_xi (gdouble a, gdouble k, gdouble s, gdouble t)
{
  return log (k) + (1.0 + 2.0 * a) * log (s * t);
}

gdouble
test_csq1d_bessel_F1 (gdouble a, gdouble k, gdouble s, gdouble t)
{
  return 0.5 * (1.0 + 2.0 * a) / (k * t);
}

gdouble
test_csq1d_bessel_F2 (gdouble a, gdouble k, gdouble s, gdouble t)
{
  return -0.25 * (1.0 + 2.0 * a) / gsl_pow_2 (k * t);
}

struct _TestCSQ1DBessel
{
  NcmCSQ1D parent_instance;
  gdouble a;
  gdouble k;
  gdouble s;
};

G_DEFINE_TYPE (TestCSQ1DBessel, test_csq1d_bessel, NCM_TYPE_CSQ1D)

static void
test_csq1d_bessel_init (TestCSQ1DBessel *b)
{
  b->a = 0.0;
  b->k = 0.0;
  b->s = 0.0;
}

#define _BESSEL(c) TEST_CSQ1D_BESSEL (c)

static gdouble
_test_csq1d_bessel_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return test_csq1d_bessel_xi (_BESSEL (csq1d)->a, _BESSEL (csq1d)->k, _BESSEL (csq1d)->s, t);
}

static gdouble
_test_csq1d_bessel_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return _BESSEL (csq1d)->k;
}

static gdouble
_test_csq1d_bessel_eval_nu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return gsl_pow_2 (_BESSEL (csq1d)->k);
}

static gdouble
_test_csq1d_bessel_eval_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return test_csq1d_bessel_m (_BESSEL (csq1d)->a, _BESSEL (csq1d)->k, _BESSEL (csq1d)->s, t);
}

static gdouble
_test_csq1d_bessel_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return test_csq1d_bessel_F1 (_BESSEL (csq1d)->a, _BESSEL (csq1d)->k, _BESSEL (csq1d)->s, t);
}

static gdouble
_test_csq1d_bessel_eval_F2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return test_csq1d_bessel_F2 (_BESSEL (csq1d)->a, _BESSEL (csq1d)->k, _BESSEL (csq1d)->s, t);
}

/* int dt / m = -s (s t)^(-2 a) / (2 a) */
static gdouble
_test_csq1d_bessel_eval_int_1_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  const gdouble a = _BESSEL (csq1d)->a;
  const gdouble s = _BESSEL (csq1d)->s;

  return -s *pow (s *t, -2.0 *a) / (2.0 * a);
}

/* With q = int dt / m: int m nu^2 = -k^2 (s t)^(2 + 2a) / (2a + 2),
 * int q m nu^2 = -k^2 t^2 / (4 a), int q^2 m nu^2 = -k^2 (s t)^(2 - 2a) / (8 a^2 (1 - a)). */
static gdouble
_test_csq1d_bessel_eval_int_mnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  const gdouble a = _BESSEL (csq1d)->a;
  const gdouble k = _BESSEL (csq1d)->k;
  const gdouble s = _BESSEL (csq1d)->s;

  return -k *k *pow (s *t, 2.0 + 2.0 *a) / (2.0 * a + 2.0);
}

static gdouble
_test_csq1d_bessel_eval_int_qmnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  const gdouble a = _BESSEL (csq1d)->a;
  const gdouble k = _BESSEL (csq1d)->k;

  return -gsl_pow_2 (k * t) / (4.0 * a);
}

static gdouble
_test_csq1d_bessel_eval_int_q2mnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  const gdouble a = _BESSEL (csq1d)->a;
  const gdouble k = _BESSEL (csq1d)->k;
  const gdouble s = _BESSEL (csq1d)->s;

  return -k *k *pow (s *t, 2.0 - 2.0 *a) / (8.0 * a * a * (1.0 - a));
}

static void
test_csq1d_bessel_class_init (TestCSQ1DBesselClass *klass)
{
  NcmCSQ1DClass *csq1d_class = NCM_CSQ1D_CLASS (klass);

  csq1d_class->eval_xi         = &_test_csq1d_bessel_eval_xi;
  csq1d_class->eval_nu         = &_test_csq1d_bessel_eval_nu;
  csq1d_class->eval_nu2        = &_test_csq1d_bessel_eval_nu2;
  csq1d_class->eval_m          = &_test_csq1d_bessel_eval_m;
  csq1d_class->eval_F1         = &_test_csq1d_bessel_eval_F1;
  csq1d_class->eval_F2         = &_test_csq1d_bessel_eval_F2;
  csq1d_class->eval_int_1_m    = &_test_csq1d_bessel_eval_int_1_m;
  csq1d_class->eval_int_mnu2   = &_test_csq1d_bessel_eval_int_mnu2;
  csq1d_class->eval_int_qmnu2  = &_test_csq1d_bessel_eval_int_qmnu2;
  csq1d_class->eval_int_q2mnu2 = &_test_csq1d_bessel_eval_int_q2mnu2;
}

TestCSQ1DBessel *
test_csq1d_bessel_new (gdouble a, gdouble k, gboolean adiab)
{
  TestCSQ1DBessel *b = g_object_new (TEST_TYPE_CSQ1D_BESSEL, NULL);

  b->a = a;
  b->k = k;
  b->s = adiab ? -1.0 : 1.0;

  return b;
}

struct _TestCSQ1DBesselMin
{
  NcmCSQ1D parent_instance;
  gdouble a;
  gdouble k;
  gdouble s;
};

G_DEFINE_TYPE (TestCSQ1DBesselMin, test_csq1d_bessel_min, NCM_TYPE_CSQ1D)

static void
test_csq1d_bessel_min_init (TestCSQ1DBesselMin *b)
{
  b->a = 0.0;
  b->k = 0.0;
  b->s = 0.0;
}

#define _BESSEL_MIN(c) TEST_CSQ1D_BESSEL_MIN (c)

static gdouble
_test_csq1d_bessel_min_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return test_csq1d_bessel_xi (_BESSEL_MIN (csq1d)->a, _BESSEL_MIN (csq1d)->k, _BESSEL_MIN (csq1d)->s, t);
}

static gdouble
_test_csq1d_bessel_min_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return _BESSEL_MIN (csq1d)->k;
}

static gdouble
_test_csq1d_bessel_min_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t)
{
  return test_csq1d_bessel_F1 (_BESSEL_MIN (csq1d)->a, _BESSEL_MIN (csq1d)->k, _BESSEL_MIN (csq1d)->s, t);
}

static void
test_csq1d_bessel_min_class_init (TestCSQ1DBesselMinClass *klass)
{
  NcmCSQ1DClass *csq1d_class = NCM_CSQ1D_CLASS (klass);

  csq1d_class->eval_xi = &_test_csq1d_bessel_min_eval_xi;
  csq1d_class->eval_nu = &_test_csq1d_bessel_min_eval_nu;
  csq1d_class->eval_F1 = &_test_csq1d_bessel_min_eval_F1;
}

TestCSQ1DBesselMin *
test_csq1d_bessel_min_new (gdouble a, gdouble k, gboolean adiab)
{
  TestCSQ1DBesselMin *b = g_object_new (TEST_TYPE_CSQ1D_BESSEL_MIN, NULL);

  b->a = a;
  b->k = k;
  b->s = adiab ? -1.0 : 1.0;

  return b;
}

