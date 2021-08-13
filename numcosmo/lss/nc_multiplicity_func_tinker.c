/***************************************************************************
 *            nc_multiplicity_func_tinker.c
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
 * SECTION:nc_multiplicity_func_tinker
 * @title: NcMultiplicityFuncTinker
 * @short_description: Dark matter halo -- Tinker multiplicity function.
 *
 * FIXME
 * Reference: arxiv:0803.2706
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_tinker.h"

struct _NcMultiplicityFuncTinkerPrivate
{
  NcMultiplicityFuncMassDef mdef;
  gdouble (*eval) (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);
  gdouble A0;
  gdouble a0;
  gdouble b0;
  gdouble c;
  gdouble Delta;
};

enum
{
  PROP_0,
  PROP_A0,
  PROP_A1,
  PROP_B0,
  PROP_C,
  PROP_DELTA,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncTinker, nc_multiplicity_func_tinker, NC_TYPE_MULTIPLICITY_FUNC);

static gdouble _nc_multiplicity_func_tinker_eval_error (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z) { g_error ("method eval not correctly initialized by %s.", G_OBJECT_TYPE_NAME (mulf)); return 0.0;}

static void
nc_multiplicity_func_tinker_init (NcMultiplicityFuncTinker *mt)
{
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv = nc_multiplicity_func_tinker_get_instance_private (mt);

  self->mdef  = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;
  self->eval  = &_nc_multiplicity_func_tinker_eval_error;
  self->A0    = 0.0;
  self->a0    = 0.0;
  self->b0    = 0.0;
  self->c     = 0.0;
  self->Delta = 0.0;
}

static void
_nc_multiplicity_func_tinker_set_property (GObject * object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncTinker *mt = NC_MULTIPLICITY_FUNC_TINKER (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_TINKER (object));
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;

  switch (prop_id)
  {
    case PROP_DELTA:
      nc_multiplicity_func_tinker_set_Delta (mt, g_value_get_double (value));
      break;
    case PROP_A0:
      self->A0 = g_value_get_double (value);
      break;
    case PROP_A1:
      self->a0 = g_value_get_double (value);
      break;
    case PROP_B0:
      self->b0 = g_value_get_double (value);
      break;
    case PROP_C:
      self->c = g_value_get_double (value);
      break;        
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_tinker_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncTinker *mt = NC_MULTIPLICITY_FUNC_TINKER (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_TINKER (object));
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;

  switch (prop_id)
  {
    case PROP_DELTA:
      g_value_set_double (value, nc_multiplicity_func_tinker_get_Delta (mt));
      break;
    case PROP_A0:
      g_value_set_double (value, self->A0);
      break;
    case PROP_A1:
      g_value_set_double (value, self->a0);
      break;
    case PROP_B0:
      g_value_set_double (value, self->b0);
      break; 
    case PROP_C:
      g_value_set_double (value, self->c);
      break;     
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_tinker_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_tinker_parent_class)->finalize (object);
}

static void _nc_multiplicity_func_tinker_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef); 
static NcMultiplicityFuncMassDef _nc_multiplicity_func_tinker_get_mdef (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_tinker_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

static void
nc_multiplicity_func_tinker_class_init (NcMultiplicityFuncTinkerClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = _nc_multiplicity_func_tinker_set_property;
  object_class->get_property = _nc_multiplicity_func_tinker_get_property;
  object_class->finalize = _nc_multiplicity_func_tinker_finalize;

  /**
   * NcMultiplicityFuncTinker:A0:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A0,
                                   g_param_spec_double ("A0",
                                                        NULL,
                                                        "A0",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.186,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncTinker:a0:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A1,
                                   g_param_spec_double ("a0",
                                                        NULL,
                                                        "a0",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.47,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcMultiplicityFuncTinker:b0:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_B0,
                                   g_param_spec_double ("b0",
                                                        NULL,
                                                        "b0",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 2.57,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncTinker:c:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_C,
                                   g_param_spec_double ("c",
                                                        NULL,
                                                        "c",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.19,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncTinker:Delta:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA,
                                   g_param_spec_double ("Delta",
                                                        NULL,
                                                        "Delta",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 200.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  parent_class->set_mdef = &_nc_multiplicity_func_tinker_set_mdef;
  parent_class->get_mdef = &_nc_multiplicity_func_tinker_get_mdef;
  parent_class->eval     = &_nc_multiplicity_func_tinker_eval;
}

/* Eval mean and critical */

#define log10_200  2.301029995663981195213738894724493026768 //log10 (200.0);
#define log10_300  2.477121254719662437295027903255115309202 //log10 (300.0);
#define log10_400  2.602059991327962390427477789448986053536 //log10 (400.0);
#define log10_600  2.778151250383643632508766797979608335970 //log10 (600.0);
#define log10_800  2.903089986991943585641216684173479080304 //log10 (800.0);
#define log10_1200 3.079181246047624827722505692704101362737 //log10 (1200.0);
#define log10_1600 3.204119982655924780854955578897972107072 //log10 (1600.0);
#define log10_2400 3.380211241711606022936244587428594389505 //log10 (2400.0);
#define log10_3200 3.505149978319905976068694473622465133840 //log10 (3200.0);

static const gdouble
_nc_multiplicity_func_tinker_coef_A[8][6] = {
    {log10_200, log10_300, 0.186, 0.0, 0.2, 0.5},
    {log10_300, log10_400, 0.2, 0.5, 0.212, -1.56},
    {log10_400, log10_600, 0.212, -1.56, 0.218, 3.05},
    {log10_600, log10_800, 0.218, 3.05, 0.248, -2.95},
    {log10_800, log10_1200, 0.248, -2.95, 0.255, 1.07},
    {log10_1200, log10_1600, 0.255, 1.07, 0.26, -0.71},
    {log10_1600, log10_2400, 0.26, -0.71, 0.26, 0.21},
    {log10_2400, log10_3200, 0.26, 0.21, 0.26, 0.0}};

static const gdouble 
_nc_multiplicity_func_tinker_coef_a[8][6] = {
    {log10_200, log10_300, 1.47, 0.0, 1.52, 1.19},
    {log10_300, log10_400, 1.52, 1.19, 1.46, -6.34},
    {log10_400, log10_600, 1.46, -6.34, 1.61, 21.36},
    {log10_600, log10_800, 1.61, 21.36, 1.87, -10.95},
    {log10_800, log10_1200, 1.87, -10.95, 2.13, 2.59},
    {log10_1200, log10_1600, 2.13, 2.59, 2.3, -0.85},
    {log10_1600, log10_2400, 2.3, -0.85, 2.53, -2.07},
    {log10_2400, log10_3200, 2.53, -2.07, 2.66, 0.0}};

static const gdouble 
_nc_multiplicity_func_tinker_coef_b[8][6] = {
    {log10_200, log10_300, 2.57, 0.0, 2.25, -1.08},
    {log10_300, log10_400, 2.25, -1.08, 2.05, 12.61},
    {log10_400, log10_600, 2.05, 12.61, 1.87, -20.96},
    {log10_600, log10_800, 1.87, -20.96, 1.59, 24.08},
    {log10_800, log10_1200, 1.59, 24.08, 1.51, -6.64},
    {log10_1200, log10_1600, 1.51, -6.64, 1.46, 3.84},
    {log10_1600, log10_2400, 1.46, 3.84, 1.44, -2.09},
    {log10_2400, log10_3200, 1.44, -2.09, 1.41, 0.0}};

static const gdouble 
_nc_multiplicity_func_tinker_coef_c[8][6] = {
    {log10_200, log10_300, 1.19, 0.0, 1.27, 0.94},
    {log10_300, log10_400, 1.27, 0.94, 1.34, -0.43},
    {log10_400, log10_600, 1.34, -0.43, 1.45, 4.61},
    {log10_600, log10_800, 1.45, 4.61, 1.58, 0.01},
    {log10_800, log10_1200, 1.58, 0.01, 1.8, 1.21},
    {log10_1200, log10_1600, 1.8, 1.21, 1.97, 1.43},
    {log10_1600, log10_2400, 1.97, 1.43, 2.24, 0.33},
    {log10_2400, log10_3200, 2.24, 0.33, 2.44, 0.0}};

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

  //printf ("% 20.15g % 20.15g % 20.15g\n% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n% 20.15g % 20.15g % 20.15g % 20.15g\n", x, y, dx, d[0], d[1], 
  //        p0, y * p1, y*y/2.0*p2, y*y*y/6.0*p3);

  return P;
}

static gdouble
_nc_multiplicity_func_tinker_mean_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)   /* $f(\sigma)$ Tinker: astro-ph/0803.2706 */
{
  NcMultiplicityFuncTinker *mt = NC_MULTIPLICITY_FUNC_TINKER (mulf);
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;

  const gdouble Delta = self->Delta;
  const gdouble log10_Delta = log10 (Delta);
  gint i;
  gdouble A0 = 0.0, a0 = 0.0, b0 = 0.0, c = 0.0;

  g_assert (Delta <= 3200.0);

  NCM_UNUSED (cosmo);

  for (i = 0; i < 8; i++)
  {
    if (log10_Delta >= _nc_multiplicity_func_tinker_coef_A[i][1])
      continue;
    else 
    {
      A0 = calc_polynomial (_nc_multiplicity_func_tinker_coef_A[i], log10_Delta);
      a0 = calc_polynomial (_nc_multiplicity_func_tinker_coef_a[i], log10_Delta);
      b0 = calc_polynomial (_nc_multiplicity_func_tinker_coef_b[i], log10_Delta);
      c  = calc_polynomial (_nc_multiplicity_func_tinker_coef_c[i], log10_Delta);
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

static gdouble
_nc_multiplicity_func_tinker_crit_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)   /* $f(\sigma)$ Tinker: astro-ph/0803.2706 */
{
  NcMultiplicityFuncTinker *mt = NC_MULTIPLICITY_FUNC_TINKER (mulf);
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;
  const gdouble E2      = nc_hicosmo_E2 (cosmo, z);
  const gdouble Omega_m = nc_hicosmo_E2Omega_m (cosmo, z) / E2;
  const gdouble Delta_z = self->Delta / Omega_m;
  const gdouble log10_Delta_z = log10 (Delta_z);
  gint i;
  gdouble A0 = 0.0, a0 = 0.0, b0 = 0.0, c = 0.0;

  if (log10_Delta_z > log10_3200)
  {
    const gdouble a0_fit_3200 = 1.43 + pow (log10_3200 - 2.3, 1.5);
    const gdouble b0_fit_3200 = 1.0 + pow (log10_3200 - 1.6, -1.5);
    const gdouble c_fit_3200 = 1.2 + pow (log10_3200 - 2.35, 1.6);
    const gdouble A0_s_3200 = calc_polynomial (_nc_multiplicity_func_tinker_coef_A[7], log10_3200);
    const gdouble a0_s_3200 = calc_polynomial (_nc_multiplicity_func_tinker_coef_a[7], log10_3200);
    const gdouble b0_s_3200 = calc_polynomial (_nc_multiplicity_func_tinker_coef_b[7], log10_3200);
    const gdouble c_s_3200  = calc_polynomial (_nc_multiplicity_func_tinker_coef_c[7], log10_3200);

    A0 = A0_s_3200;
    a0 = a0_s_3200 / a0_fit_3200 * (1.43 + pow (log10_Delta_z - 2.3, 1.5));
    b0 = b0_s_3200 / b0_fit_3200 * (1.0 + pow (log10_3200 - 1.6, -1.5));
    c = c_s_3200 / c_fit_3200 * (1.2 + pow (log10_3200 - 2.35, 1.6));
  }
  else
  {
    for (i = 0; i < 8; i++)
    {
      if (log10_Delta_z >= _nc_multiplicity_func_tinker_coef_A[i][1])
        continue;
      else
      {
        A0 = calc_polynomial (_nc_multiplicity_func_tinker_coef_A[i], log10_Delta_z);
        a0 = calc_polynomial (_nc_multiplicity_func_tinker_coef_a[i], log10_Delta_z);
        b0 = calc_polynomial (_nc_multiplicity_func_tinker_coef_b[i], log10_Delta_z);
        c  = calc_polynomial (_nc_multiplicity_func_tinker_coef_c[i], log10_Delta_z);
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

static void 
_nc_multiplicity_func_tinker_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncTinker *mt = NC_MULTIPLICITY_FUNC_TINKER (mulf);
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;

  switch (mdef)
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
      self->eval = &_nc_multiplicity_func_tinker_mean_eval;
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
      self->eval = &_nc_multiplicity_func_tinker_crit_eval;
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
      g_error ("NcMultiplicityFuncTinker does not support virial mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_FOF:
      g_error ("NcMultiplicityFuncTinker does not support fof mass def");
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  self->mdef = mdef;
}

static NcMultiplicityFuncMassDef 
_nc_multiplicity_func_tinker_get_mdef (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncTinker *mt = NC_MULTIPLICITY_FUNC_TINKER (mulf);
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;

  return self->mdef;
}

static gdouble
_nc_multiplicity_func_tinker_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)   /* $f(\sigma)$ Tinker: astro-ph/0803.2706 */
{
  NcMultiplicityFuncTinker *mt = NC_MULTIPLICITY_FUNC_TINKER (mulf);
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;
  
  return self->eval(mulf, cosmo, sigma, z);
}


/**
 * nc_multiplicity_func_tinker_new:
 *   
 * FIXME
 *
 * Returns: A new #NcMultiplicityFuncTinker.
 */
NcMultiplicityFuncTinker *
nc_multiplicity_func_tinker_new (void)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_TINKER,
                       "mass-def", NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN,
                       NULL);
}

/**
 * nc_multiplicity_func_tinker_new_full:
 * @mdef: a #NcMultiplicityFuncMassDef
 * @Delta: parameter that multiplies the background mass density (mean ou critical)   
 * 
 * FIXME
 *
 * Returns: A new #NcMultiplicityFuncTinker.
 */
NcMultiplicityFuncTinker *
nc_multiplicity_func_tinker_new_full (NcMultiplicityFuncMassDef mdef, gdouble Delta)
{
  
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_TINKER,
                       "mass-def", mdef,
                       "Delta",    Delta,
                       NULL);
}

/**
 * nc_multiplicity_func_tinker_ref:
 * @mt: a #NcMultiplicityFuncTinker
 *
 * Increases the reference count of @mt by one.
 *
 * Returns: (transfer full): @mt
 */
NcMultiplicityFuncTinker *
nc_multiplicity_func_tinker_ref (NcMultiplicityFuncTinker *mt)
{
  return g_object_ref (mt);
}

/**
 * nc_multiplicity_func_tinker_free:
 * @mt: a #NcMultiplicityFuncTinker
 *
 * Atomically decrements the reference count of @mt by one. If the reference count drops to 0,
 * all memory allocated by @mt is released.
 *
 */
void
nc_multiplicity_func_tinker_free (NcMultiplicityFuncTinker *mt)
{
  g_object_unref (mt);
}

/**
 * nc_multiplicity_func_tinker_clear:
 * @mt: a #NcMultiplicityFuncTinker
 *
 * Atomically decrements the reference count of @mt by one. If the reference count drops to 0,
 * all memory allocated by @mt is released. Set the pointer to NULL;
 *
 */
void
nc_multiplicity_func_tinker_clear (NcMultiplicityFuncTinker **mt)
{
  g_clear_object (mt);
}

/**
 * nc_multiplicity_func_tinker_set_Delta:
 * @mt: a #NcMultiplicityFuncTinker.
 * @Delta: value of #NcMultiplicityFuncTinker:Delta.
 *
 * Sets the value @Delta to the #NcMultiplicityFuncTinker:Delta property.
 *
 */
void
nc_multiplicity_func_tinker_set_Delta (NcMultiplicityFuncTinker *mt, gdouble Delta)
{
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;

  g_assert (Delta >= 0);

  self->Delta = Delta;
}

/**
 * nc_multiplicity_func_tinker_get_Delta:
 * @mt: a #NcMultiplicityFuncTinker.
 *
 * Returns: the value of #NcMultiplicityFuncTinker:Delta property.
 */
gdouble
nc_multiplicity_func_tinker_get_Delta (const NcMultiplicityFuncTinker *mt)
{
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;
  
  return self->Delta;
}



