/***************************************************************************
 *            nc_multiplicity_func_tinker_mean.c
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
 * SECTION:nc_multiplicity_func_tinker_mean
 * @title: NcMultiplicityFuncTinkerMean
 * @short_description: Dark matter halo -- Tinker multiplicity function mean matter density. 
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_tinker_mean.h"

G_DEFINE_TYPE (NcMultiplicityFuncTinkerMean, nc_multiplicity_func_tinker_mean, NC_TYPE_MULTIPLICITY_FUNC);

enum
{
  PROP_0,
  PROP_DELTA
};

/**
 * nc_multiplicity_func_tinker_mean_new:
 * @Delta: FIXME 
 *   
 * FIXME
 *
 * Returns: A new #NcMultiplicityFunc.
 */
NcMultiplicityFunc *
nc_multiplicity_func_tinker_mean_new (gdouble Delta)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN,
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
_nc_multiplicity_func_tinker_mean_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)   /* $f(\sigma)$ Tinker: astro-ph/0803.2706 */
{
  /* This function is a copy of the nc_multiplicity_function_tinker_critical, adapted to mean matter density.*/
  NcMultiplicityFuncTinkerMean *mulf_tinker_mean = NC_MULTIPLICITY_FUNC_TINKER_MEAN (mulf);
  //const gdouble Omega_m0 = nc_hicosmo_Omega_m0 (NC_HICOSMO (model));
  //const gdouble E2 = nc_hicosmo_E2 (NC_HICOSMO (model), z);
  //const gdouble Delta_z = mulf_tinker_mean->Delta * E2 / (Omega_m0 * gsl_pow_3 (1.0 + z));
  const gdouble Delta = mulf_tinker_mean->Delta;
  const gdouble log10_Delta = log10 (Delta);
  const gdouble log10_200 = log10 (200.0);
  const gdouble log10_300 = log10 (300.0);
  const gdouble log10_400 = log10 (400.0);
  const gdouble log10_600 = log10 (600.0);
  const gdouble log10_800 = log10 (800.0);
  const gdouble log10_1200 = log10 (1200.0);
  const gdouble log10_1600 = log10 (1600.0);
  const gdouble log10_2400 = log10 (2400.0);
  const gdouble log10_3200 = log10 (3200.0);

  const gdouble coef_A[8][6] = {
    {log10_200, log10_300, 0.186, 0.0, 0.2, 0.5},
    {log10_300, log10_400, 0.2, 0.5, 0.212, -1.56},
    {log10_400, log10_600, 0.212, -1.56, 0.218, 3.05},
    {log10_600, log10_800, 0.218, 3.05, 0.248, -2.95},
    {log10_800, log10_1200, 0.248, -2.95, 0.255, 1.07},
    {log10_1200, log10_1600, 0.255, 1.07, 0.26, -0.71},
    {log10_1600, log10_2400, 0.26, -0.71, 0.26, 0.21},
    {log10_2400, log10_3200, 0.26, 0.21, 0.26, 0.0}};
  const gdouble coef_a[8][6] = {
    {log10_200, log10_300, 1.47, 0.0, 1.52, 1.19},
    {log10_300, log10_400, 1.52, 1.19, 1.46, -6.34},
    {log10_400, log10_600, 1.46, -6.34, 1.61, 21.36},
    {log10_600, log10_800, 1.61, 21.36, 1.87, -10.95},
    {log10_800, log10_1200, 1.87, -10.95, 2.13, 2.59},
    {log10_1200, log10_1600, 2.13, 2.59, 2.3, -0.85},
    {log10_1600, log10_2400, 2.3, -0.85, 2.53, -2.07},
    {log10_2400, log10_3200, 2.53, -2.07, 2.66, 0.0}};
  const gdouble coef_b[8][6] = {
    {log10_200, log10_300, 2.57, 0.0, 2.25, -1.08},
    {log10_300, log10_400, 2.25, -1.08, 2.05, 12.61},
    {log10_400, log10_600, 2.05, 12.61, 1.87, -20.96},
    {log10_600, log10_800, 1.87, -20.96, 1.59, 24.08},
    {log10_800, log10_1200, 1.59, 24.08, 1.51, -6.64},
    {log10_1200, log10_1600, 1.51, -6.64, 1.46, 3.84},
    {log10_1600, log10_2400, 1.46, 3.84, 1.44, -2.09},
    {log10_2400, log10_3200, 1.44, -2.09, 1.41, 0.0}};
  const gdouble coef_c[8][6] = {
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

  g_assert (Delta <= 3200.0);

  NCM_UNUSED (cosmo);

  for (i = 0; i < 8; i++)
  {
    if (log10_Delta >= coef_A[i][1])
      continue;
    else 
    {
      A0 = calc_polynomial (coef_A[i], log10_Delta);
      a0 = calc_polynomial (coef_a[i], log10_Delta);
      b0 = calc_polynomial (coef_b[i], log10_Delta);
      c = calc_polynomial (coef_c[i], log10_Delta);
      break;
    }
  }

  {
    const gdouble A = A0 * pow(1.0 + z, -0.14);
    const gdouble a = a0 * pow(1.0 + z, -0.06);
    const gdouble log10alpha = -pow(0.75 / log10 (Delta / 75.0), 1.2);
    const gdouble alpha = pow(10.0, log10alpha);
    const gdouble b = b0 * pow(1.0 + z, -alpha);
    const gdouble f_Tinker_mean = A * (pow(sigma / b, -a) + 1.0) * exp(-c / (sigma * sigma));

    //printf ("NEW % 5.5g % 5.5g % 20.15g % 20.15g % 20.15g % 20.15g\n", z, mulf_tinker_mean->Delta, Delta_z, E2, Omega_m0, gsl_pow_3 (1.0 + z));
    //printf ("NC: A = %22.15g a = %22.15g b = %22.15g c = %22.15g sigma = %22.15g\n", A, a, b, c, sigma);
		
    return f_Tinker_mean;
  }
}

/**
 * nc_multiplicity_func_tinker_mean_set_Delta:
 * @mulf_tinker_mean: a #NcMultiplicityFuncTinkerMean.
 * @Delta: value of #NcMultiplicityFuncTinkerMean:Delta.
 *
 * Sets the value @Delta to the #NcMultiplicityFuncTinkerMean:Delta property.
 *
 */
void
nc_multiplicity_func_tinker_mean_set_Delta (NcMultiplicityFuncTinkerMean *mulf_tinker_mean, gdouble Delta)
{
  g_assert (Delta >= 0);
  mulf_tinker_mean->Delta = Delta;
}

/**
 * nc_multiplicity_func_tinker_mean_get_Delta:
 * @mulf_tinker_mean: a #NcMultiplicityFuncTinkerMean.
 *
 * Returns: the value of #NcMultiplicityFuncTinkerMean:Delta property.
 */
gdouble
nc_multiplicity_func_tinker_mean_get_Delta (const NcMultiplicityFuncTinkerMean *mulf_tinker_mean)
{
  return mulf_tinker_mean->Delta;
}

static void
nc_multiplicity_func_tinker_mean_init (NcMultiplicityFuncTinkerMean *mulf_tinker_mean)
{
  mulf_tinker_mean->Delta = 200.0;
}

static void
_nc_multiplicity_func_tinker_mean_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_multiplicity_func_tinker_mean_parent_class)->finalize (object);
}

static void
_nc_multiplicity_func_tinker_mean_set_property (GObject * object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncTinkerMean *mulf_tinker_mean = NC_MULTIPLICITY_FUNC_TINKER_MEAN (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_TINKER_MEAN (object));

  switch (prop_id)
  {
    case PROP_DELTA:
      nc_multiplicity_func_tinker_mean_set_Delta (mulf_tinker_mean, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_tinker_mean_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncTinkerMean *mulf_tinker_mean = NC_MULTIPLICITY_FUNC_TINKER_MEAN (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_TINKER_MEAN (object));

  switch (prop_id)
  {
    case PROP_DELTA:
      g_value_set_double (value, mulf_tinker_mean->Delta);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_multiplicity_func_tinker_mean_class_init (NcMultiplicityFuncTinkerMeanClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  parent_class->eval = &_nc_multiplicity_func_tinker_mean_eval;

  object_class->finalize = _nc_multiplicity_func_tinker_mean_finalize;
  object_class->set_property = _nc_multiplicity_func_tinker_mean_set_property;
  object_class->get_property = _nc_multiplicity_func_tinker_mean_get_property;

  /**
   * NcMultiplicityFuncTinkerMean:Delta:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA,
                                   g_param_spec_double ("Delta",
                                                        NULL,
                                                        "Delta",
                                                        200.0, 3200.0, 200.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

