/***************************************************************************
 *            nc_observer_velocity.c
 *
 *  Sun Oct 18 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_observer_velocity.c
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

/**
 * NcObserverVelocity:
 *
 * Velocity of the observer relative to the rest frame of the cosmological model.
 *
 * The parameters are the Cartesian components $v_x$, $v_y$, $v_z$ [km/s] of the observer's
 * velocity in the frame of the sky directions of the data (equatorial for the catalogs of
 * NumCosmo), $\vec\beta = \vec v / c$. A source seen in the direction $\hat n$ at the
 * observed redshift $z_\mathrm{obs}$ has, in the rest frame,
 * $$
 * 1 + z_\mathrm{rest} = (1 + z_\mathrm{obs})\, \mathcal{D}, \qquad
 * \mathcal{D} = \frac{1}{\gamma\,(1 - \vec\beta\cdot\hat n)}, \qquad
 * \gamma = \frac{1}{\sqrt{1 - \beta^2}},
 * $$
 * and its luminosity distance transforms as
 * $D_{L,\mathrm{obs}} = D_{L,\mathrm{rest}} / \mathcal{D}$, so the observed distance modulus is
 * that of the rest frame at $z_\mathrm{rest}$ minus $5\log_{10}\mathcal{D}$. To first order
 * in $\beta$ this is the dipole $\delta\mu \simeq -\frac{5}{\ln 10}\,\vec\beta\cdot\hat n\,
 * \frac{1 + z}{z\,(1 + z/2)}$ of a bulk flow of the sources with velocity $-\vec v$, which
 * the data cannot tell from the observer's motion. The aberration of $\hat n$ itself is of
 * order $\beta$ in the direction and changes the modulus only at second order, so the
 * observed direction is used in the rest frame.
 *
 * The redshifts of a catalog in the CMB frame already carry the Sun's velocity relative to
 * the CMB, and the model velocity is then that of the CMB frame relative to the rest frame
 * of the cosmological model. The defaults, $\vec v = 0$, leave a likelihood unchanged.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_observer_velocity.h"
#include "ncm/core/ncm_c.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcObserverVelocity
{
  /*< private >*/
  NcmModel parent_instance;
};

enum
{
  PROP_0,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcObserverVelocity, nc_observer_velocity, NCM_TYPE_MODEL)

#define VECTOR (NCM_MODEL (ov))
#define VX     (ncm_model_orig_param_get (VECTOR, NC_OBSERVER_VELOCITY_VX))
#define VY     (ncm_model_orig_param_get (VECTOR, NC_OBSERVER_VELOCITY_VY))
#define VZ     (ncm_model_orig_param_get (VECTOR, NC_OBSERVER_VELOCITY_VZ))

static void
nc_observer_velocity_init (NcObserverVelocity *ov)
{
}

static void
_nc_observer_velocity_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_OBSERVER_VELOCITY (object));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_observer_velocity_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  g_return_if_fail (NC_IS_OBSERVER_VELOCITY (object));

  switch (prop_id)
  {
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_observer_velocity_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_observer_velocity_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_observer_velocity, NC_TYPE_OBSERVER_VELOCITY);

static void
nc_observer_velocity_class_init (NcObserverVelocityClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_observer_velocity_set_property;
  model_class->get_property = &_nc_observer_velocity_get_property;
  object_class->finalize    = &_nc_observer_velocity_finalize;

  ncm_model_class_set_name_nick (model_class, "Observer velocity", "ObserverVelocity");
  ncm_model_class_add_params (model_class, NC_OBSERVER_VELOCITY_SPARAM_LEN, 0, PROP_SIZE);

  ncm_mset_model_register_id (model_class,
                              "NcObserverVelocity",
                              "Velocity of the observer relative to the rest frame of the model.",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_set_sparam (model_class, NC_OBSERVER_VELOCITY_VX, "v_x", "vx",
                              -NC_OBSERVER_VELOCITY_MAX_V, NC_OBSERVER_VELOCITY_MAX_V, 10.0,
                              NC_OBSERVER_VELOCITY_DEFAULT_PARAMS_ABSTOL, NC_OBSERVER_VELOCITY_DEFAULT_V,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_OBSERVER_VELOCITY_VY, "v_y", "vy",
                              -NC_OBSERVER_VELOCITY_MAX_V, NC_OBSERVER_VELOCITY_MAX_V, 10.0,
                              NC_OBSERVER_VELOCITY_DEFAULT_PARAMS_ABSTOL, NC_OBSERVER_VELOCITY_DEFAULT_V,
                              NCM_PARAM_TYPE_FIXED);
  ncm_model_class_set_sparam (model_class, NC_OBSERVER_VELOCITY_VZ, "v_z", "vz",
                              -NC_OBSERVER_VELOCITY_MAX_V, NC_OBSERVER_VELOCITY_MAX_V, 10.0,
                              NC_OBSERVER_VELOCITY_DEFAULT_PARAMS_ABSTOL, NC_OBSERVER_VELOCITY_DEFAULT_V,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);
}

/**
 * nc_observer_velocity_new:
 *
 * Creates a new #NcObserverVelocity at rest.
 *
 * Returns: (transfer full): a new #NcObserverVelocity
 */
NcObserverVelocity *
nc_observer_velocity_new (void)
{
  NcObserverVelocity *ov = g_object_new (NC_TYPE_OBSERVER_VELOCITY, NULL);

  return ov;
}

/**
 * nc_observer_velocity_ref:
 * @ov: a #NcObserverVelocity
 *
 * Increases the reference count of @ov by one.
 *
 * Returns: (transfer full): @ov
 */
NcObserverVelocity *
nc_observer_velocity_ref (NcObserverVelocity *ov)
{
  return g_object_ref (ov);
}

/**
 * nc_observer_velocity_free:
 * @ov: a #NcObserverVelocity
 *
 * Decreases the reference count of @ov by one.
 *
 */
void
nc_observer_velocity_free (NcObserverVelocity *ov)
{
  g_object_unref (ov);
}

/**
 * nc_observer_velocity_clear:
 * @ov: a #NcObserverVelocity
 *
 * If *@ov is not %NULL, decreases its reference count and sets *@ov to %NULL.
 *
 */
void
nc_observer_velocity_clear (NcObserverVelocity **ov)
{
  g_clear_object (ov);
}

/**
 * nc_observer_velocity_speed:
 * @ov: a #NcObserverVelocity
 *
 * Returns: $|\vec v|$ [km/s]
 */
gdouble
nc_observer_velocity_speed (NcObserverVelocity *ov)
{
  return sqrt (VX * VX + VY * VY + VZ * VZ);
}

/**
 * nc_observer_velocity_beta_n:
 * @ov: a #NcObserverVelocity
 * @nx: x component of the unit vector $\hat n$ towards the source
 * @ny: y component of $\hat n$
 * @nz: z component of $\hat n$
 *
 * Returns: $\vec\beta\cdot\hat n$
 */
gdouble
nc_observer_velocity_beta_n (NcObserverVelocity *ov, const gdouble nx, const gdouble ny, const gdouble nz)
{
  return (VX * nx + VY * ny + VZ * nz) / (1.0e-3 * ncm_c_c ());
}

/**
 * nc_observer_velocity_doppler:
 * @ov: a #NcObserverVelocity
 * @nx: x component of the unit vector $\hat n$ towards the source
 * @ny: y component of $\hat n$
 * @nz: z component of $\hat n$
 *
 * Returns: the Doppler factor $\mathcal{D} = (1 + z_\mathrm{rest}) / (1 + z_\mathrm{obs})$, see #NcObserverVelocity
 */
gdouble
nc_observer_velocity_doppler (NcObserverVelocity *ov, const gdouble nx, const gdouble ny, const gdouble nz)
{
  const gdouble beta2  = (VX * VX + VY * VY + VZ * VZ) / gsl_pow_2 (1.0e-3 * ncm_c_c ());
  const gdouble beta_n = nc_observer_velocity_beta_n (ov, nx, ny, nz);

  return sqrt (1.0 - beta2) / (1.0 - beta_n);
}

/**
 * nc_observer_velocity_z_rest:
 * @ov: a #NcObserverVelocity
 * @z_obs: observed redshift
 * @nx: x component of the unit vector $\hat n$ towards the source
 * @ny: y component of $\hat n$
 * @nz: z component of $\hat n$
 *
 * Returns: the rest-frame redshift $z_\mathrm{rest} = (1 + z_\mathrm{obs})\,\mathcal{D} - 1$
 */
gdouble
nc_observer_velocity_z_rest (NcObserverVelocity *ov, const gdouble z_obs, const gdouble nx, const gdouble ny, const gdouble nz)
{
  return (1.0 + z_obs) * nc_observer_velocity_doppler (ov, nx, ny, nz) - 1.0;
}

