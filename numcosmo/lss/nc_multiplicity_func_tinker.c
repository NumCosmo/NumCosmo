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
 * NcMultiplicityFuncTinker:
 *
 * Dark matter halo -- Tinker multiplicity function.
 *
 * The Tinker multiplicity function is a parametric form for the mass function of
 * dark matter halos.
 * Reference: arxiv:0803.2706
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_tinker.h"
#include "math/ncm_spline_cubic_d2.h"
#include "math/ncm_spline_gsl.h"

struct _NcMultiplicityFuncTinkerPrivate
{
  NcMultiplicityFuncMassDef mdef;
  gdouble (*eval) (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);
  gdouble A0;
  gdouble a0;
  gdouble b0;
  gdouble c;
  gdouble Delta;
  NcmSpline *A_s;
  NcmSpline *a_s;
  NcmSpline *b_s;
  NcmSpline *c_s;
};

enum
{
  PROP_0,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncTinker, nc_multiplicity_func_tinker, NC_TYPE_MULTIPLICITY_FUNC)

static gdouble
_nc_multiplicity_func_tinker_eval_error (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  g_error ("method eval not correctly initialized by %s.", G_OBJECT_TYPE_NAME (mulf));

  return 0.0;
}

#define log10_200  2.301029995663981195213738894724493026768 /*log10 (200.0); */
#define log10_300  2.477121254719662437295027903255115309202 /*log10 (300.0); */
#define log10_400  2.602059991327962390427477789448986053536 /*log10 (400.0); */
#define log10_600  2.778151250383643632508766797979608335970 /*log10 (600.0); */
#define log10_800  2.903089986991943585641216684173479080304 /*log10 (800.0); */
#define log10_1200 3.079181246047624827722505692704101362737 /*log10 (1200.0); */
#define log10_1600 3.204119982655924780854955578897972107072 /*log10 (1600.0); */
#define log10_2400 3.380211241711606022936244587428594389505 /*log10 (2400.0); */
#define log10_3200 3.505149978319905976068694473622465133840 /*log10 (3200.0); */

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
  self->A_s   = NULL;
  self->a_s   = NULL;
  self->b_s   = NULL;
  self->c_s   = NULL;

  nc_multiplicity_func_tinker_set_linear_interp (mt, FALSE);
}

static void
_nc_multiplicity_func_tinker_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_TINKER (object));
  /*NcMultiplicityFuncTinkerPrivate * const self = mt->priv;*/

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_tinker_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_TINKER (object));
  /* NcMultiplicityFuncTinkerPrivate * const self = mt->priv;*/

  switch (prop_id)
  {
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_tinker_dispose (GObject *object)
{
  NcMultiplicityFuncTinker *mt                 = NC_MULTIPLICITY_FUNC_TINKER (object);
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;

  ncm_spline_clear (&self->A_s);
  ncm_spline_clear (&self->a_s);
  ncm_spline_clear (&self->b_s);
  ncm_spline_clear (&self->c_s);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_tinker_parent_class)->dispose (object);
}

static void
_nc_multiplicity_func_tinker_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_tinker_parent_class)->finalize (object);
}

static void _nc_multiplicity_func_tinker_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef);
static void _nc_multiplicity_func_tinker_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta);
static NcMultiplicityFuncMassDef _nc_multiplicity_func_tinker_get_mdef (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_tinker_get_Delta (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_tinker_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

static void
nc_multiplicity_func_tinker_class_init (NcMultiplicityFuncTinkerClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass *parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = &_nc_multiplicity_func_tinker_set_property;
  object_class->get_property = &_nc_multiplicity_func_tinker_get_property;
  object_class->dispose      = &_nc_multiplicity_func_tinker_dispose;
  object_class->finalize     = &_nc_multiplicity_func_tinker_finalize;


  parent_class->set_mdef  = &_nc_multiplicity_func_tinker_set_mdef;
  parent_class->set_Delta = &_nc_multiplicity_func_tinker_set_Delta;
  parent_class->get_mdef  = &_nc_multiplicity_func_tinker_get_mdef;
  parent_class->get_Delta = &_nc_multiplicity_func_tinker_get_Delta;
  parent_class->eval      = &_nc_multiplicity_func_tinker_eval;
}

static gdouble
_nc_multiplicity_func_tinker_mean_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z) /* $f(\sigma)$ Tinker: astro-ph/0803.2706 */
{
  NcMultiplicityFuncTinker *mt                 = NC_MULTIPLICITY_FUNC_TINKER (mulf);
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;

  const gdouble Delta = self->Delta;

  const gdouble A             = self->A0 * pow (1.0 + z, -0.14);
  const gdouble a             = self->a0 * pow (1.0 + z, -0.06);
  const gdouble log10alpha    = -pow (0.75 / log10 (Delta / 75.0), 1.2);
  const gdouble alpha         = pow (10.0, log10alpha);
  const gdouble b             = self->b0 * pow (1.0 + z, -alpha);
  const gdouble f_Tinker_mean = A * (pow (sigma / b, -a) + 1.0) * exp (-self->c / (sigma * sigma));

  return f_Tinker_mean;
}

static gdouble
_nc_multiplicity_func_tinker_crit_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z) /* $f(\sigma)$ Tinker: astro-ph/0803.2706 */
{
  NcMultiplicityFuncTinker *mt                 = NC_MULTIPLICITY_FUNC_TINKER (mulf);
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;
  const gdouble E2                             = nc_hicosmo_E2 (cosmo, z);
  const gdouble Omega_m                        = nc_hicosmo_E2Omega_m (cosmo, z) / E2;
  const gdouble Delta_z                        = self->Delta / Omega_m;
  const gdouble log10_Delta_z                  = log10 (Delta_z);
  gdouble A0 = 0.0, a0 = 0.0, b0 = 0.0, c = 0.0;

  if (log10_Delta_z > log10_3200)
  {
    const gdouble a0_fit_3200 = 1.43 + pow (log10_3200 - 2.3, 1.5);
    const gdouble b0_fit_3200 = 1.0 + pow (log10_3200 - 1.6, -1.5);
    const gdouble c_fit_3200  = 1.2 + pow (log10_3200 - 2.35, 1.6);
    const gdouble A0_s_3200   = ncm_spline_eval (self->A_s, log10_3200);
    const gdouble a0_s_3200   = ncm_spline_eval (self->a_s, log10_3200);
    const gdouble b0_s_3200   = ncm_spline_eval (self->b_s, log10_3200);
    const gdouble c_s_3200    = ncm_spline_eval (self->c_s, log10_3200);

    A0 = A0_s_3200;
    a0 = a0_s_3200 / a0_fit_3200 * (1.43 + pow (log10_Delta_z - 2.3, 1.5));
    b0 = b0_s_3200 / b0_fit_3200 * (1.0 + pow (log10_3200 - 1.6, -1.5));
    c  = c_s_3200 / c_fit_3200 * (1.2 + pow (log10_3200 - 2.35, 1.6));
  }
  else
  {
    A0 = ncm_spline_eval (self->A_s, log10_Delta_z);
    a0 = ncm_spline_eval (self->a_s, log10_Delta_z);
    b0 = ncm_spline_eval (self->b_s, log10_Delta_z);
    c  = ncm_spline_eval (self->c_s, log10_Delta_z);
  }

  {
    const gdouble A                     = A0 * pow (1.0 + z, -0.14);
    const gdouble a                     = a0 * pow (1.0 + z, -0.06);
    const gdouble log10alpha            = -pow (0.75 / log10 (Delta_z / 75.0), 1.2);
    const gdouble alpha                 = pow (10.0, log10alpha);
    const gdouble b                     = b0 * pow (1.0 + z, -alpha);
    const gdouble f_Tinker_Delta_spline = A * (pow (sigma / b, -a) + 1.0) * exp (-c / (sigma * sigma));

    return f_Tinker_Delta_spline;
  }
}

static void
_nc_multiplicity_func_tinker_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncTinker *mt                 = NC_MULTIPLICITY_FUNC_TINKER (mulf);
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
  NcMultiplicityFuncTinker *mt                 = NC_MULTIPLICITY_FUNC_TINKER (mulf);
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;

  return self->mdef;
}

static gdouble
_nc_multiplicity_func_tinker_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z) /* $f(\sigma)$ Tinker: astro-ph/0803.2706 */
{
  NcMultiplicityFuncTinker *mt                 = NC_MULTIPLICITY_FUNC_TINKER (mulf);
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;

  return self->eval (mulf, cosmo, sigma, z);
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
static void
_nc_multiplicity_func_tinker_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta)
{
  NcMultiplicityFuncTinker *mt                 = NC_MULTIPLICITY_FUNC_TINKER (mulf);
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;
  const gdouble log10_Delta                    = log10 (Delta);

  g_assert (Delta >= 0);
  g_assert (Delta <= 3200.0);

  self->Delta = Delta;

  self->A0 = ncm_spline_eval (self->A_s, log10_Delta);
  self->a0 = ncm_spline_eval (self->a_s, log10_Delta);
  self->b0 = ncm_spline_eval (self->b_s, log10_Delta);
  self->c  = ncm_spline_eval (self->c_s, log10_Delta);
}

/**
 * nc_multiplicity_func_tinker_get_Delta:
 * @mt: a #NcMultiplicityFuncTinker.
 *
 * Returns: the value of #NcMultiplicityFuncTinker:Delta property.
 */
gdouble
_nc_multiplicity_func_tinker_get_Delta (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncTinker *mt                 = NC_MULTIPLICITY_FUNC_TINKER (mulf);
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;

  return self->Delta;
}

/**
 * nc_multiplicity_func_tinker_set_linear_interp:
 * @mulf: a #NcMultiplicityFuncTinker.
 * @lin_interp: a @gboolean
 *
 * If @lin_interp is true uses linear interpolation to compute the
 * coefficients A, a0, b0 and c. Otherwise it uses cubic interpolation
 * as described in arxiv:0803.2706.
 *
 */
void
nc_multiplicity_func_tinker_set_linear_interp (NcMultiplicityFuncTinker *mt, gboolean lin_interp)
{
  NcMultiplicityFuncTinkerPrivate * const self = mt->priv;
  NcMultiplicityFunc *mulf                     = NC_MULTIPLICITY_FUNC (mt);
  const gdouble Delta                          = nc_multiplicity_func_get_Delta (mulf);

  gdouble log10_delta[9] = {log10_200, log10_300,  log10_400,  log10_600, log10_800, log10_1200, log10_1600, log10_2400, log10_3200};
  gdouble coeff_A[9]     = {0.186, 0.200, 0.212, 0.218, 0.248, 0.255, 0.260, 0.260, 0.260};
  gdouble coeff_a[9]     = {1.470, 1.520, 1.560, 1.610, 1.870, 2.130, 2.300, 2.530, 2.660};
  gdouble coeff_b[9]     = {2.570, 2.250, 2.050, 1.870, 1.590, 1.510, 1.460, 1.440, 1.410};
  gdouble coeff_c[9]     = {1.190, 1.270, 1.340, 1.450, 1.580, 1.800, 1.970, 2.240, 2.440};

  gdouble d2f_A[9] = {0.000,  0.500, -1.560,   3.050,  -2.950,  1.070, -0.710,  0.210, 0.000};
  gdouble d2f_a[9] = {0.000,  1.190, -6.340,  21.360, -10.950,  2.590, -0.850, -2.070, 0.000};
  gdouble d2f_b[9] = {0.000, -1.080, 12.610, -20.960,  24.080, -6.640,  3.840, -2.090, 0.000};
  gdouble d2f_c[9] = {0.000,  0.940, -0.430,   4.610,   0.010,  1.210,  1.430,  0.330, 0.000};

  NcmVector *log10_delta_v = ncm_vector_new_data_dup (log10_delta, 9, 1);
  NcmVector *coeff_A_v     = ncm_vector_new_data_dup (coeff_A, 9, 1);
  NcmVector *coeff_a_v     = ncm_vector_new_data_dup (coeff_a, 9, 1);
  NcmVector *coeff_b_v     = ncm_vector_new_data_dup (coeff_b, 9, 1);
  NcmVector *coeff_c_v     = ncm_vector_new_data_dup (coeff_c, 9, 1);
  NcmVector *d2f_A_v       = ncm_vector_new_data_dup (d2f_A, 9, 1);
  NcmVector *d2f_a_v       = ncm_vector_new_data_dup (d2f_a, 9, 1);
  NcmVector *d2f_b_v       = ncm_vector_new_data_dup (d2f_b, 9, 1);
  NcmVector *d2f_c_v       = ncm_vector_new_data_dup (d2f_c, 9, 1);

  ncm_spline_clear (&self->A_s);
  ncm_spline_clear (&self->a_s);
  ncm_spline_clear (&self->b_s);
  ncm_spline_clear (&self->c_s);

  if (!lin_interp)
  {
    self->A_s = NCM_SPLINE (ncm_spline_cubic_d2_new (log10_delta_v, coeff_A_v, d2f_A_v, TRUE));
    self->a_s = NCM_SPLINE (ncm_spline_cubic_d2_new (log10_delta_v, coeff_a_v, d2f_a_v, TRUE));
    self->b_s = NCM_SPLINE (ncm_spline_cubic_d2_new (log10_delta_v, coeff_b_v, d2f_b_v, TRUE));
    self->c_s = NCM_SPLINE (ncm_spline_cubic_d2_new (log10_delta_v, coeff_c_v, d2f_c_v, TRUE));
  }
  else
  {
    self->A_s = NCM_SPLINE (ncm_spline_gsl_new_full_by_id (NCM_SPLINE_GSL_LINEAR, log10_delta_v, coeff_A_v, TRUE));
    self->a_s = NCM_SPLINE (ncm_spline_gsl_new_full_by_id (NCM_SPLINE_GSL_LINEAR, log10_delta_v, coeff_a_v, TRUE));
    self->b_s = NCM_SPLINE (ncm_spline_gsl_new_full_by_id (NCM_SPLINE_GSL_LINEAR, log10_delta_v, coeff_b_v, TRUE));
    self->c_s = NCM_SPLINE (ncm_spline_gsl_new_full_by_id (NCM_SPLINE_GSL_LINEAR, log10_delta_v, coeff_c_v, TRUE));
  }

  ncm_vector_free (log10_delta_v);
  ncm_vector_free (coeff_A_v);
  ncm_vector_free (coeff_a_v);
  ncm_vector_free (coeff_b_v);
  ncm_vector_free (coeff_c_v);
  ncm_vector_free (d2f_A_v);
  ncm_vector_free (d2f_a_v);
  ncm_vector_free (d2f_b_v);
  ncm_vector_free (d2f_c_v);

  if (self->Delta > 0.0)
    _nc_multiplicity_func_tinker_set_Delta (mulf, Delta);
}

