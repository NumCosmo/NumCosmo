/***************************************************************************
 *            nc_multiplicity_func_bocquet.c
 *
 *  Mon Aug 16 14:39:25 2021
 *  Copyright  2021  Mariana Penna Lima and Cinthia Nunes de Lima
 *  <pennalima@gmail.com>, <cinthia.n.lima@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2021 Mariana Penna Lima and Cinthia Nunes de Lima 
 * <pennalima@gmail.com>, <cinthia.n.lima@uel.br>
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
 * SECTION:nc_multiplicity_func_bocquet
 * @title: NcMultiplicityFuncBocquet
 * @short_description: Dark matter halo -- Bocquet multiplicity function.
 *
 * FIXME
 *
 * The multiplicity function $ f(\sigma) $ is
 * \begin{equation}\label{3}
 * f(\sigma) = A \left[ \left( \frac{\sigma}{b}\right)^{-a} + 1 \right] exp \left( - \frac{c}{\sigma^{2}} \right), 
 * \end{equation}
 * where the parameters A, a, b and c are calibrated through the simulations.
 *
 * For $\Delta_{200m}$:
 * \begin{equation}\label{4}
 * \begin{split}
 * &A(z) = A_{0} (1 + z)^{A_{z}} \\
 * &a(z) = a_{0} (1 + z)^{a_{z}} \\
 * &b(z) = b_{0} (1 + z)^{b_{z}} \\
 * &c(z) = c_{0} (1 + z)^{c_{z}}.
 * \end{split}
 * \end{equation}
 *
 * Obs.: the subscript 0 denotes $ z = 0$.
 *
 * Reference: arxiv:1502.07357
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_bocquet.h"
#include "numcosmo/nc_enum_types.h"

struct _NcMultiplicityFuncBocquetPrivate
{
  NcMultiplicityFuncMassDef mdef;
  NcMultiplicityFuncBocquetSim sim;
  gdouble A0;
  gdouble a0;
  gdouble b0;
  gdouble c0;
  gdouble Az;
  gdouble az;
  gdouble bz;
  gdouble cz;
  gdouble Delta;
  gboolean constructed;
};

enum
{
  PROP_0,
  PROP_SIM,
  PROP_DELTA,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncBocquet, nc_multiplicity_func_bocquet, NC_TYPE_MULTIPLICITY_FUNC);

static void
nc_multiplicity_func_bocquet_init (NcMultiplicityFuncBocquet *mb)
{
  NcMultiplicityFuncBocquetPrivate * const self = mb->priv = nc_multiplicity_func_bocquet_get_instance_private (mb);

  self->mdef  = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;
  self->sim   = NC_MULTIPLICITY_FUNC_BOCQUET_SIM_LEN;
  self->Delta = 0.0;
  self->A0    = 0.0;
  self->a0    = 0.0;
  self->b0    = 0.0;
  self->c0    = 0.0;
  self->Az    = 0.0;
  self->az    = 0.0;
  self->bz    = 0.0;
  self->cz    = 0.0; 
  self->constructed = FALSE;
}

static void
_nc_multiplicity_func_bocquet_set_property (GObject * object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncBocquet *mb = NC_MULTIPLICITY_FUNC_BOCQUET (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_BOCQUET (object));

  switch (prop_id)
  {
    case PROP_DELTA:
      nc_multiplicity_func_bocquet_set_Delta (mb, g_value_get_double (value));
      break;
    case PROP_SIM:
      nc_multiplicity_func_bocquet_set_sim (mb, g_value_get_enum (value));
      break;  
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_bocquet_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncBocquet *mb = NC_MULTIPLICITY_FUNC_BOCQUET (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_BOCQUET (object));

  switch (prop_id)
  {
    case PROP_DELTA:
      g_value_set_double (value, nc_multiplicity_func_bocquet_get_Delta (mb));
      break;
    case PROP_SIM:
      g_value_set_enum (value, nc_multiplicity_func_bocquet_get_sim (mb));
      break;  
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void _nc_multiplicity_func_bocquet_set_all (NcMultiplicityFuncBocquet *mb);

static void
_nc_multiplicity_func_bocquet_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_multiplicity_func_bocquet_parent_class)->constructed (object);

  {
    NcMultiplicityFuncBocquet *mb = NC_MULTIPLICITY_FUNC_BOCQUET (object);
    NcMultiplicityFuncBocquetPrivate * const self = mb->priv;

    _nc_multiplicity_func_bocquet_set_all (mb);
    
    self->constructed = TRUE;
  }
}

static void
_nc_multiplicity_func_bocquet_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_bocquet_parent_class)->finalize (object);
}

static void _nc_multiplicity_func_bocquet_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef); 
static NcMultiplicityFuncMassDef _nc_multiplicity_func_bocquet_get_mdef (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_bocquet_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);
static gboolean _nc_multiplicity_func_bocquet_has_correction_factor (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_bocquet_correction_factor (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z, gdouble lnM);

static void
nc_multiplicity_func_bocquet_class_init (NcMultiplicityFuncBocquetClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = _nc_multiplicity_func_bocquet_set_property;
  object_class->get_property = _nc_multiplicity_func_bocquet_get_property;
  object_class->constructed  = _nc_multiplicity_func_bocquet_constructed;
  object_class->finalize     = _nc_multiplicity_func_bocquet_finalize;

  /**
   * NcMultiplicityFuncBocquet:Delta:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA,
                                   g_param_spec_double ("Delta",
                                                        NULL,
                                                        "Delta",
                                                        200.0, G_MAXDOUBLE, 200.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncBocquet:sim:
   *
   * Simulation used: dark matter only (DM) or hydrodynamical (HYDRO)
   */
  g_object_class_install_property (object_class,
                                   PROP_SIM,
                                   g_param_spec_enum ("sim",
                                                      NULL,
                                                      "Simulation type",
                                                      NC_TYPE_MULTIPLICITY_FUNC_BOCQUET_SIM, NC_MULTIPLICITY_FUNC_BOCQUET_SIM_DM,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  

  parent_class->set_mdef = &_nc_multiplicity_func_bocquet_set_mdef;
  parent_class->get_mdef = &_nc_multiplicity_func_bocquet_get_mdef;
  parent_class->eval     = &_nc_multiplicity_func_bocquet_eval;
  parent_class->has_correction_factor = &_nc_multiplicity_func_bocquet_has_correction_factor;
  parent_class->correction_factor     = &_nc_multiplicity_func_bocquet_correction_factor;
}

static void
_nc_multiplicity_func_bocquet_set_all (NcMultiplicityFuncBocquet *mb)
{
  NcMultiplicityFuncBocquetPrivate * const self = mb->priv;

  switch (self->sim)
  {
    case NC_MULTIPLICITY_FUNC_BOCQUET_SIM_DM:
      switch (self->mdef)
      {
        case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
          if (self->Delta != 200.0)
            g_error ("NcMultiplicityFuncBocquet does not support Delta != 200 (mass def - mean density).");
          self->A0 = 0.175; 
          self->a0 = 1.53;
          self->b0 = 2.55;
          self->c0 = 1.19;
          self->Az = -0.012; 
          self->az = -0.040;
          self->bz = -0.194;
          self->cz = -0.021;
          break;
        case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
          if (self->Delta != 200.0 && self->Delta != 500.0)
            g_error ("NcMultiplicityFuncBocquet does not support Delta != 200 or 500 (mass def - critical density).");
          else if (self->Delta == 200.0)
          {
            self->A0 = 0.222; 
            self->a0 = 1.71;
            self->b0 = 2.24;
            self->c0 = 1.46;
            self->Az = 0.269; 
            self->az = 0.321;
            self->bz = -0.621;
            self->cz = -0.153;  
          }
          else
          {
            self->A0 = 0.241; 
            self->a0 = 2.18;
            self->b0 = 2.35;
            self->c0 = 2.02;
            self->Az = 0.370; 
            self->az = 0.251;
            self->bz = -0.698;
            self->cz = -0.310;  
          }
          break;
        case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
        case NC_MULTIPLICITY_FUNC_MASS_DEF_FOF:
          g_error ("NcMultiplicityFuncBocquet does not support virial or FOF mass def.");
          break;
        default:
          g_error ("NcMultiplicityFuncMassDef does not have this option (%d).", self->mdef);
          break;  
      }
      break;
    case NC_MULTIPLICITY_FUNC_BOCQUET_SIM_HYDRO:
      switch (self->mdef)
      {
        case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
          if (self->Delta != 200.0)
            g_error ("NcMultiplicityFuncBocquet does not support Delta != 200 (mass def - mean density).");
          self->A0 = 0.228; 
          self->a0 = 2.15;
          self->b0 = 1.69;
          self->c0 = 1.30;
          self->Az = 0.285; 
          self->az = -0.058;
          self->bz = -0.366;
          self->cz = -0.045;
          break;
        case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
          if (self->Delta != 200.0 && self->Delta != 500.0)
            g_error ("NcMultiplicityFuncBocquet does not support Delta != 200 or 500 (mass def - critical density).");
          else if (self->Delta == 200.0)
          {
            self->A0 = 0.202; 
            self->a0 = 2.21;
            self->b0 = 2.00;
            self->c0 = 1.57;
            self->Az = 1.147; 
            self->az = 0.375;
            self->bz = -1.074;
            self->cz = -0.196;  
          }
          else
          {
            self->A0 = 0.180; 
            self->a0 = 2.29;
            self->b0 = 2.44;
            self->c0 = 1.97;
            self->Az = 1.088; 
            self->az = 0.150;
            self->bz = -1.008;
            self->cz = -0.322;  
          }
          break;
        case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
        case NC_MULTIPLICITY_FUNC_MASS_DEF_FOF:
          g_error ("NcMultiplicityFuncBocquet does not support virial or FOF mass def.");
          break;
        default:
          g_error ("NcMultiplicityFuncMassDef does not have this option.");
          break;  
      }
      break;
      default:
        g_error ("NcMultiplicityFuncBocquetSim does not have this option (%d).", self->sim);
        break; 
  }  
}


static void 
_nc_multiplicity_func_bocquet_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncBocquet *mb = NC_MULTIPLICITY_FUNC_BOCQUET (mulf);
  NcMultiplicityFuncBocquetPrivate * const self = mb->priv;

  self->mdef = mdef;
  if (self->constructed)
    _nc_multiplicity_func_bocquet_set_all (mb);  
}

static NcMultiplicityFuncMassDef 
_nc_multiplicity_func_bocquet_get_mdef (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncBocquet *mb = NC_MULTIPLICITY_FUNC_BOCQUET (mulf);
  NcMultiplicityFuncBocquetPrivate * const self = mb->priv;

  return self->mdef;
}

static gboolean 
_nc_multiplicity_func_bocquet_has_correction_factor (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncBocquet *mb = NC_MULTIPLICITY_FUNC_BOCQUET (mulf);
  NcMultiplicityFuncBocquetPrivate * const self = mb->priv;
  
  if (self->mdef == NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL)
    return TRUE;
  else
    return FALSE;   
}

static gdouble
_nc_multiplicity_func_bocquet_correction_factor (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z, gdouble lnM)
{
  NcMultiplicityFuncBocquet *mt = NC_MULTIPLICITY_FUNC_BOCQUET (mulf);
  NcMultiplicityFuncBocquetPrivate * const self = mt->priv;
  const gdouble Omega_m = nc_hicosmo_E2Omega_m (cosmo, z);
  
  if (self->Delta == 200.0)
  {
    const gdouble gamma0 = 3.54e-2 + pow(Omega_m, 0.09);
    const gdouble gamma1 = 4.56e-2 + 2.68e-2 / Omega_m;
    const gdouble gamma2 = 0.721 + 3.50e-2 / Omega_m;
    const gdouble gamma3 = 0.628 + 0.164 / Omega_m;
    const gdouble delta0 = -1.67e-2 + 2.18e-2 * Omega_m;
    const gdouble delta1 = 6.52e-3 - 6.86e-3 * Omega_m;
    const gdouble gamma = gamma0 + gamma1 * exp(- gsl_pow_2((gamma2 - z) / gamma3));
    const gdouble delta = delta0 + delta1 * z;
    const gdouble M200c_M200m = gamma + delta * lnM;
    
    return M200c_M200m;
  }
  else
  {
    g_assert (self->Delta == 500.0);
    const gdouble alpha0 = 0.880 + 0.329 * Omega_m;
    const gdouble alpha1 = 1.00 + 4.31e-2 / Omega_m;
    const gdouble alpha2 = -0.365 + 0.254 / Omega_m;
    const gdouble alpha = alpha0 * (alpha1 * z + alpha2) / (z + alpha2);
    const gdouble beta = -1.7e-2 + 3.74e-3 * Omega_m;
    const gdouble M500c_M200m = alpha + beta * lnM;

    return M500c_M200m;
  }
}

static gdouble
_nc_multiplicity_func_bocquet_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)   /* $f(\sigma)$ Bocquet: astro-ph/0803.2706 */
{
  NcMultiplicityFuncBocquet *mt = NC_MULTIPLICITY_FUNC_BOCQUET (mulf);
  NcMultiplicityFuncBocquetPrivate * const self = mt->priv;

  const gdouble A = self->A0 * pow(1.0 + z, self->Az);
  const gdouble a = self->a0 * pow(1.0 + z, self->az);
  const gdouble b = self->b0 * pow(1.0 + z, self->bz);
  const gdouble c = self->c0 * pow(1.0 + z, self->cz);

  const gdouble f_Bocquet = A * (pow(sigma / b, -a) + 1.0) * exp(-c / (sigma * sigma));
  
  return f_Bocquet;
}


/**
 * nc_multiplicity_func_bocquet_new:
 *   
 * FIXME
 *
 * Returns: A new #NcMultiplicityFuncBocquet.
 */
NcMultiplicityFuncBocquet *
nc_multiplicity_func_bocquet_new (void)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_BOCQUET,
                       NULL);
}

/**
 * nc_multiplicity_func_bocquet_new_full:
 * @mdef: a #NcMultiplicityFuncMassDef
 * @sim: a #NcMultiplicityFuncBocquetSim
 * @Delta: parameter that multiplies the background mass density (mean ou critical)   
 * 
 * FIXME
 *
 * Returns: A new #NcMultiplicityFuncBocquet.
 */
NcMultiplicityFuncBocquet *
nc_multiplicity_func_bocquet_new_full (NcMultiplicityFuncMassDef mdef, NcMultiplicityFuncBocquetSim sim, gdouble Delta)
{
  
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_BOCQUET,
                       "mass-def", mdef,
                       "sim",      sim,
                       "Delta",    Delta,
                       NULL);
}

/**
 * nc_multiplicity_func_bocquet_ref:
 * @mb: a #NcMultiplicityFuncBocquet
 *
 * Increases the reference count of @mb by one.
 *
 * Returns: (transfer full): @mb
 */
NcMultiplicityFuncBocquet *
nc_multiplicity_func_bocquet_ref (NcMultiplicityFuncBocquet *mb)
{
  return g_object_ref (mb);
}

/**
 * nc_multiplicity_func_bocquet_free:
 * @mb: a #NcMultiplicityFuncBocquet
 *
 * Atomically decrements the reference count of @mb by one. If the reference count drops to 0,
 * all memory allocated by @mb is released.
 *
 */
void
nc_multiplicity_func_bocquet_free (NcMultiplicityFuncBocquet *mb)
{
  g_object_unref (mb);
}

/**
 * nc_multiplicity_func_bocquet_clear:
 * @mb: a #NcMultiplicityFuncBocquet
 *
 * Atomically decrements the reference count of @mt by one. If the reference count drops to 0,
 * all memory allocated by @mb is released. Set the pointer to NULL;
 *
 */
void
nc_multiplicity_func_bocquet_clear (NcMultiplicityFuncBocquet **mb)
{
  g_clear_object (mb);
}

/**
 * nc_multiplicity_func_bocquet_set_Delta:
 * @mb: a #NcMultiplicityFuncBocquet.
 * @Delta: value of #NcMultiplicityFuncBocquet:Delta.
 *
 * Sets the value @Delta to the #NcMultiplicityFuncBocquet:Delta property.
 *
 */
void
nc_multiplicity_func_bocquet_set_Delta (NcMultiplicityFuncBocquet *mb, gdouble Delta)
{
  NcMultiplicityFuncBocquetPrivate * const self = mb->priv;

  self->Delta = Delta;
  if (self->constructed)
    _nc_multiplicity_func_bocquet_set_all (mb);
}

/**
 * nc_multiplicity_func_bocquet_get_Delta:
 * @mb: a #NcMultiplicityFuncBocquet.
 *
 * Returns: the value of #NcMultiplicityFuncBocquet:Delta property.
 */
gdouble
nc_multiplicity_func_bocquet_get_Delta (const NcMultiplicityFuncBocquet *mb)
{
  NcMultiplicityFuncBocquetPrivate * const self = mb->priv;
  
  return self->Delta;
}


/**
 * nc_multiplicity_func_bocquet_set_sim:
 * @mb: a #NcMultiplicityFuncBocquet
 * @sim: a #NcMultiplicityFuncBocquetSim
 *
 * Sets the value @sim to the #NcMultiplicityFuncBocquet:sim property.
 *
 */
void
nc_multiplicity_func_bocquet_set_sim (NcMultiplicityFuncBocquet *mb, NcMultiplicityFuncBocquetSim sim)
{
  NcMultiplicityFuncBocquetPrivate * const self = mb->priv;

  self->sim = sim;
  if (self->constructed)
    _nc_multiplicity_func_bocquet_set_all (mb);
}

/**
 * nc_multiplicity_func_bocquet_get_sim:
 * @mb: a #NcMultiplicityFuncBocquet
 *
 * Returns: a #NcMultiplicityFuncBocquetSim.
 */
NcMultiplicityFuncBocquetSim
nc_multiplicity_func_bocquet_get_sim (const NcMultiplicityFuncBocquet *mb)
{
  NcMultiplicityFuncBocquetPrivate * const self = mb->priv;

  return self->sim;
}
