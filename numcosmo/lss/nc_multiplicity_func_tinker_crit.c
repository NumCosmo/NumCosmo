/***************************************************************************
 *            nc_multiplicity_func_tinker_crit.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_multiplicity_func_tinker_crit
 * @title: NcMultiplicityFuncTinkerCrit
 * @short_description: Dark matter halo -- Tinker multiplicity function critical density.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_tinker_crit.h"

G_DEFINE_TYPE (NcMultiplicityFuncTinkerCrit, nc_multiplicity_func_tinker_crit, NC_TYPE_MULTIPLICITY_FUNC);

enum
{
  PROP_0,
  PROP_DELTA
};

/**
 * nc_multiplicity_func_tinker_crit_new:
 * @Delta: FIXME
 *
 * FIXME
 *
 * Returns: A new #NcMultiplicityFunc.
 */
NcMultiplicityFunc *
nc_multiplicity_func_tinker_crit_new (gdouble Delta)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_TINKER_CRIT,
                       "Delta", Delta,
                       NULL);
}

static gdouble
calc_polynomial (const gdouble *d, gdouble x)
{
  const gdouble y = x - d[0];
  const gdouble dx = d[1] - d[0];
  const gdouble dx2 = dx * dx;
  const gdouble p0 = d[2];
  const gdouble p2 = d[3];
  const gdouble p3 = (d[5] - p2) / dx;
  const gdouble p1 = (d[4] - p0 - (p2 + p3 * dx / 3.0) * dx2 / 2.0) / dx;
  const gdouble P = p0 + y * (p1 + y * (p2 + p3 * y / 3.0) / 2.0);

  //printf ("% 20.15g % 20.15g % 20.15g\n% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n% 20.15g % 20.15g % 20.15g % 20.15g\n", x, y, dx, d[0], d[1], d[2], d[3], d[4], d[5],
  //        p0, y * p1, y*y/2.0*p2, y*y*y/6.0*p3);

  return P;
}

static gdouble
_nc_multiplicity_func_tinker_crit_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)   /* $f(\sigma)$ Tinker: astro-ph/0803.2706 */
{
  NcMultiplicityFuncTinkerCrit *mulf_tinker_crit = NC_MULTIPLICITY_FUNC_TINKER_CRIT (mulf);
  const gdouble E2      = nc_hicosmo_E2 (cosmo, z);
  const gdouble Omega_m = nc_hicosmo_E2Omega_m (cosmo, z) / E2;
  const gdouble Delta_z = mulf_tinker_crit->Delta / Omega_m;
  const gdouble log10_Delta_z = log10 (Delta_z);
/*
  const gdouble log10_200  = 2.301029995663981195213738894724493026768; //log10 (200.0);
  const gdouble log10_300  = 2.477121254719662437295027903255115309202; //log10 (300.0);
  const gdouble log10_400  = 2.602059991327962390427477789448986053536; //log10 (400.0);
  const gdouble log10_600  = 2.778151250383643632508766797979608335970; //log10 (600.0);
  const gdouble log10_800  = 2.903089986991943585641216684173479080304; //log10 (800.0);
  const gdouble log10_1200 = 3.079181246047624827722505692704101362737; //log10 (1200.0);
  const gdouble log10_1600 = 3.204119982655924780854955578897972107072; //log10 (1600.0);
  const gdouble log10_2400 = 3.380211241711606022936244587428594389505; //log10 (2400.0);
  const gdouble log10_3200 = 3.505149978319905976068694473622465133840; //log10 (3200.0);
*/

#define log10_200  2.301029995663981195213738894724493026768 //log10 (200.0);
#define log10_300  2.477121254719662437295027903255115309202 //log10 (300.0);
#define log10_400  2.602059991327962390427477789448986053536 //log10 (400.0);
#define log10_600  2.778151250383643632508766797979608335970 //log10 (600.0);
#define log10_800  2.903089986991943585641216684173479080304 //log10 (800.0);
#define log10_1200 3.079181246047624827722505692704101362737 //log10 (1200.0);
#define log10_1600 3.204119982655924780854955578897972107072 //log10 (1600.0);
#define log10_2400 3.380211241711606022936244587428594389505 //log10 (2400.0);
#define log10_3200 3.505149978319905976068694473622465133840 //log10 (3200.0);

  static const gdouble coef_A[8][6] = {
    {log10_200, log10_300, 0.186, 0.0, 0.2, 0.5},
    {log10_300, log10_400, 0.2, 0.5, 0.212, -1.56},
    {log10_400, log10_600, 0.212, -1.56, 0.218, 3.05},
    {log10_600, log10_800, 0.218, 3.05, 0.248, -2.95},
    {log10_800, log10_1200, 0.248, -2.95, 0.255, 1.07},
    {log10_1200, log10_1600, 0.255, 1.07, 0.26, -0.71},
    {log10_1600, log10_2400, 0.26, -0.71, 0.26, 0.21},
    {log10_2400, log10_3200, 0.26, 0.21, 0.26, 0.0}};
  static const gdouble coef_a[8][6] = {
    {log10_200, log10_300, 1.47, 0.0, 1.52, 1.19},
    {log10_300, log10_400, 1.52, 1.19, 1.46, -6.34},
    {log10_400, log10_600, 1.46, -6.34, 1.61, 21.36},
    {log10_600, log10_800, 1.61, 21.36, 1.87, -10.95},
    {log10_800, log10_1200, 1.87, -10.95, 2.13, 2.59},
    {log10_1200, log10_1600, 2.13, 2.59, 2.3, -0.85},
    {log10_1600, log10_2400, 2.3, -0.85, 2.53, -2.07},
    {log10_2400, log10_3200, 2.53, -2.07, 2.66, 0.0}};
  static const gdouble coef_b[8][6] = {
    {log10_200, log10_300, 2.57, 0.0, 2.25, -1.08},
    {log10_300, log10_400, 2.25, -1.08, 2.05, 12.61},
    {log10_400, log10_600, 2.05, 12.61, 1.87, -20.96},
    {log10_600, log10_800, 1.87, -20.96, 1.59, 24.08},
    {log10_800, log10_1200, 1.59, 24.08, 1.51, -6.64},
    {log10_1200, log10_1600, 1.51, -6.64, 1.46, 3.84},
    {log10_1600, log10_2400, 1.46, 3.84, 1.44, -2.09},
    {log10_2400, log10_3200, 1.44, -2.09, 1.41, 0.0}};
  static const gdouble coef_c[8][6] = {
    {log10_200, log10_300, 1.19, 0.0, 1.27, 0.94},
    {log10_300, log10_400, 1.27, 0.94, 1.34, -0.43},
    {log10_400, log10_600, 1.34, -0.43, 1.45, 4.61},
    {log10_600, log10_800, 1.45, 4.61, 1.58, 0.01},
    {log10_800, log10_1200, 1.58, 0.01, 1.8, 1.21},
    {log10_1200, log10_1600, 1.8, 1.21, 1.97, 1.43},
    {log10_1600, log10_2400, 1.97, 1.43, 2.24, 0.33},
    {log10_2400, log10_3200, 2.24, 0.33, 2.44, 0.0}};
  gint i;
  gdouble A0 = 0.0, a0 = 0.0, b0 = 0.0, c = 0.0;

  if (log10_Delta_z > log10_3200)
  {
    const gdouble a0_fit_3200 = 1.43 + pow (log10_3200 - 2.3, 1.5);
    const gdouble b0_fit_3200 = 1.0 + pow (log10_3200 - 1.6, -1.5);
    const gdouble c_fit_3200 = 1.2 + pow (log10_3200 - 2.35, 1.6);
    const gdouble A0_s_3200 = calc_polynomial (coef_A[7], log10_3200);
    const gdouble a0_s_3200 = calc_polynomial (coef_a[7], log10_3200);
    const gdouble b0_s_3200 = calc_polynomial (coef_b[7], log10_3200);
    const gdouble c_s_3200  = calc_polynomial (coef_c[7], log10_3200);

    A0 = A0_s_3200;
    a0 = a0_s_3200 / a0_fit_3200 * (1.43 + pow (log10_Delta_z - 2.3, 1.5));
    b0 = b0_s_3200 / b0_fit_3200 * (1.0 + pow (log10_3200 - 1.6, -1.5));
    c = c_s_3200 / c_fit_3200 * (1.2 + pow (log10_3200 - 2.35, 1.6));
  }
  else
  {
    for (i = 0; i < 8; i++)
    {
      if (log10_Delta_z >= coef_A[i][1])
        continue;
      else
      {
        A0 = calc_polynomial (coef_A[i], log10_Delta_z);
        a0 = calc_polynomial (coef_a[i], log10_Delta_z);
        b0 = calc_polynomial (coef_b[i], log10_Delta_z);
        c = calc_polynomial (coef_c[i], log10_Delta_z);
        break;
      }
    }
  }

  {
    const gdouble A          = A0 * pow(1.0 + z, -0.14);
    const gdouble a          = a0 * pow(1.0 + z, -0.06);
    const gdouble log10alpha = - pow (0.75 / log10 (Delta_z / 75.0), 1.2);
    const gdouble alpha      = pow (10.0, log10alpha);
    const gdouble b          = b0 * pow(1.0 + z, -alpha);
    const gdouble f_Tinker_Delta_spline = A * (pow (sigma / b, -a) + 1.0) * exp(-c / (sigma * sigma));

    return f_Tinker_Delta_spline;
  }
}

/**
 * nc_multiplicity_func_tinker_crit_set_Delta:
 * @mulf_tinker_crit: a #NcMultiplicityFuncTinkerCrit.
 * @Delta: value of #NcMultiplicityFuncTinkerCrit:Delta.
 *
 * Sets the value @Delta to the #NcMultiplicityFuncTinkerCrit:Delta property.
 *
 */
void
nc_multiplicity_func_tinker_crit_set_Delta (NcMultiplicityFuncTinkerCrit *mulf_tinker_crit, gdouble Delta)
{
  g_assert (Delta >= 0);
  mulf_tinker_crit->Delta = Delta;
}

/**
 * nc_multiplicity_func_tinker_crit_get_Delta:
 * @mulf_tinker_crit: a #NcMultiplicityFuncTinkerCrit.
 *
 * Returns: the value of #NcMultiplicityFuncTinkerCrit:Delta property.
 */
gdouble
nc_multiplicity_func_tinker_crit_get_Delta (const NcMultiplicityFuncTinkerCrit *mulf_tinker_crit)
{
  return mulf_tinker_crit->Delta;
}

static void
nc_multiplicity_func_tinker_crit_init (NcMultiplicityFuncTinkerCrit *mulf_tinker_crit)
{
  mulf_tinker_crit->Delta = 200.0;
}

static void
_nc_multiplicity_func_tinker_crit_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_multiplicity_func_tinker_crit_parent_class)->finalize (object);
}

static void
_nc_multiplicity_func_tinker_crit_set_property (GObject * object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncTinkerCrit *mulf_tinker_crit = NC_MULTIPLICITY_FUNC_TINKER_CRIT (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_TINKER_CRIT (object));

  switch (prop_id)
  {
    case PROP_DELTA:
      nc_multiplicity_func_tinker_crit_set_Delta (mulf_tinker_crit, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_tinker_crit_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncTinkerCrit *mulf_tinker_crit = NC_MULTIPLICITY_FUNC_TINKER_CRIT (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_TINKER_CRIT (object));

  switch (prop_id)
  {
    case PROP_DELTA:
      g_value_set_double (value, mulf_tinker_crit->Delta);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}


static void
nc_multiplicity_func_tinker_crit_class_init (NcMultiplicityFuncTinkerCritClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  parent_class->eval = &_nc_multiplicity_func_tinker_crit_eval;

  object_class->finalize = _nc_multiplicity_func_tinker_crit_finalize;
  object_class->set_property = _nc_multiplicity_func_tinker_crit_set_property;
  object_class->get_property = _nc_multiplicity_func_tinker_crit_get_property;

  /**
   * NcMultiplicityFuncTinkerCrit:Delta:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA,
                                   g_param_spec_double ("Delta",
                                                        NULL,
                                                        "Delta",
                                                        200.0, 3200.0, 200.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}
