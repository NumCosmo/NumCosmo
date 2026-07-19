/***************************************************************************
 *            nc_galaxy_redshift_pop.c
 *
 *  Wed Jul 31 20:52:43 2024
 *  Copyright  2024  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_pop.c
 * Copyright (C) 2024 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyRedshiftPop:
 *
 * Class describing galaxy sample redshift distributions.
 *
 * This class describes a galaxy sample redshift distributions.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_redshift_pop.h"
#include "ncm/core/ncm_dtuple.h"

typedef struct _NcGalaxyRedshiftPopPrivate
{
  gint placeholder;
} NcGalaxyRedshiftPopPrivate;

enum
{
  PROP_0,
  PROP_LIM,
  PROP_LEN,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyRedshiftPop, nc_galaxy_redshift_pop, NCM_TYPE_MODEL)

static void
nc_galaxy_redshift_pop_init (NcGalaxyRedshiftPop *gsdrp)
{
}

static void
_nc_galaxy_redshift_pop_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyRedshiftPop *gsdrp = NC_GALAXY_REDSHIFT_POP (object);

  g_return_if_fail (NC_IS_GALAXY_REDSHIFT_POP (gsdrp));

  switch (prop_id)
  {
    case PROP_LIM:
    {
      NcmDTuple2 *lim = g_value_get_boxed (value);

      if (lim == NULL)
        g_error ("_nc_galaxy_redshift_pop_set_property: lim is NULL");

      nc_galaxy_redshift_pop_set_lim (gsdrp, lim->elements[0], lim->elements[1]);
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_redshift_pop_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyRedshiftPop *gsdrp = NC_GALAXY_REDSHIFT_POP (object);

  g_return_if_fail (NC_IS_GALAXY_REDSHIFT_POP (gsdrp));

  switch (prop_id)
  {
    case PROP_LIM:
    {
      gdouble z_min, z_max;

      nc_galaxy_redshift_pop_get_lim (gsdrp, &z_min, &z_max);

      g_value_take_boxed (value, ncm_dtuple2_new (z_min, z_max));
      break;
    }
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_redshift_pop_finalize (GObject *object)
{
  G_OBJECT_CLASS (nc_galaxy_redshift_pop_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_galaxy_redshift_pop, NC_TYPE_GALAXY_REDSHIFT_POP);

/* LCOV_EXCL_START */
static gdouble
_nc_galaxy_redshift_pop_gen (NcGalaxyRedshiftPop *gsdrp, NcmRNG *rng)
{
  g_error ("_nc_galaxy_redshift_pop_gen_z: not implemented");

  return 0.0;
}

static gdouble
_nc_galaxy_redshift_pop_eval (NcGalaxyRedshiftPop *gsdrp, gdouble z)
{
  g_error ("_nc_galaxy_redshift_pop_eval: not implemented");

  return 0.0;
}

static gdouble
_nc_galaxy_redshift_pop_ln_eval (NcGalaxyRedshiftPop *gsdrp, gdouble z)
{
  g_error ("_nc_galaxy_redshift_pop_ln_eval: not implemented");

  return 0.0;
}

static void
_nc_galaxy_redshift_pop_set_lim (NcGalaxyRedshiftPop *gsdrp, const gdouble z_min, const gdouble z_max)
{
  g_error ("_nc_galaxy_redshift_pop_set_lim: not implemented");
}

static void
_nc_galaxy_redshift_pop_get_lim (NcGalaxyRedshiftPop *gsdrp, gdouble *z_min, gdouble *z_max)
{
  g_error ("_nc_galaxy_redshift_pop_get_lim: not implemented");
}

/* LCOV_EXCL_STOP */

static void
nc_galaxy_redshift_pop_class_init (NcGalaxyRedshiftPopClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_nc_galaxy_redshift_pop_set_property;
  model_class->get_property = &_nc_galaxy_redshift_pop_get_property;
  object_class->finalize    = &_nc_galaxy_redshift_pop_finalize;

  ncm_model_class_set_name_nick (model_class, "Galaxy sample redshift distribution", "GalaxyRedshiftPop");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  /**
   * NcGalaxyRedshiftPop:lim:
   *
   * Galaxy sample redshift distribution limits.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LIM,
                                   g_param_spec_boxed ("lim",
                                                       NULL,
                                                       "Galaxy sample redshift distribution limits",
                                                       NCM_TYPE_DTUPLE2,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_mset_model_register_id (model_class, "NcGalaxyRedshiftPop", "Galaxy sample redshift distribution", NULL, FALSE, NCM_MSET_MODEL_MAIN);
  ncm_model_class_check_params_info (model_class);

  klass->gen     = &_nc_galaxy_redshift_pop_gen;
  klass->eval    = &_nc_galaxy_redshift_pop_eval;
  klass->ln_eval = &_nc_galaxy_redshift_pop_ln_eval;
  klass->set_lim = &_nc_galaxy_redshift_pop_set_lim;
  klass->get_lim = &_nc_galaxy_redshift_pop_get_lim;
}

/**
 * nc_galaxy_redshift_pop_ref:
 * @gsdrp: a #NcGalaxyRedshiftPop
 *
 * Increases the reference count of @gsdrp by one.
 *
 * Returns: (transfer full): @gsdrp.
 */
NcGalaxyRedshiftPop *
nc_galaxy_redshift_pop_ref (NcGalaxyRedshiftPop *gsdrp)
{
  return g_object_ref (gsdrp);
}

/**
 * nc_galaxy_redshift_pop_free:
 * @gsdrp: a #NcGalaxyRedshiftPop
 *
 * Decreases the reference count of @gsdrp by one.
 *
 */
void
nc_galaxy_redshift_pop_free (NcGalaxyRedshiftPop *gsdrp)
{
  g_object_unref (gsdrp);
}

/**
 * nc_galaxy_redshift_pop_clear:
 * @gsdrp: a #NcGalaxyRedshiftPop
 *
 * Decreases the reference count of @gsdrp by one, and sets the
 * pointer @gsdrp to NULL.
 *
 */
void
nc_galaxy_redshift_pop_clear (NcGalaxyRedshiftPop **gsdrp)
{
  g_clear_object (gsdrp);
}

/**
 * nc_galaxy_redshift_pop_set_lim:
 * @gsdrp: a #NcGalaxyRedshiftPop
 * @z_min: a #gdouble representing minimum redshift
 * @z_max: a #gdouble representing maximum redshift
 *
 * Sets the redshift limits of the galaxy sample redshift distribution.
 *
 */
void
nc_galaxy_redshift_pop_set_lim (NcGalaxyRedshiftPop *gsdrp, const gdouble z_min, const gdouble z_max)
{
  NC_GALAXY_REDSHIFT_POP_GET_CLASS (gsdrp)->set_lim (gsdrp, z_min, z_max);
}

/**
 * nc_galaxy_redshift_pop_get_lim:
 * @gsdrp: a #NcGalaxyRedshiftPop
 * @z_min: (out): a #gdouble pointer representing minimum redshift
 * @z_max: (out): a #gdouble pointer representing maximum redshift
 *
 * Gets the redshift limits of the galaxy sample redshift distribution.
 *
 */
void
nc_galaxy_redshift_pop_get_lim (NcGalaxyRedshiftPop *gsdrp, gdouble *z_min, gdouble *z_max)
{
  NC_GALAXY_REDSHIFT_POP_GET_CLASS (gsdrp)->get_lim (gsdrp, z_min, z_max);
}

/**
 * nc_galaxy_redshift_pop_gen:
 * @gsdrp: a #NcGalaxyRedshiftPop
 * @rng: a #NcmRNG
 *
 * Generates a redshift value from the galaxy sample redshift distribution.
 *
 * Returns: the generated redshift.
 */
gdouble
nc_galaxy_redshift_pop_gen (NcGalaxyRedshiftPop *gsdrp, NcmRNG *rng)
{
  return NC_GALAXY_REDSHIFT_POP_GET_CLASS (gsdrp)->gen (gsdrp, rng);
}

/**
 * nc_galaxy_redshift_pop_eval:
 * @gsdrp: a #NcGalaxyRedshiftPop
 * @z: the redshift
 *
 * Evaluates the galaxy sample redshift distribution at redshift @z.
 *
 * Returns: the probability density at $z$, $P(z)$.
 */
gdouble
nc_galaxy_redshift_pop_eval (NcGalaxyRedshiftPop *gsdrp, gdouble z)
{
  return NC_GALAXY_REDSHIFT_POP_GET_CLASS (gsdrp)->eval (gsdrp, z);
}

/**
 * nc_galaxy_redshift_pop_ln_eval:
 * @gsdrp: a #NcGalaxyRedshiftPop
 * @z: the redshift
 *
 * Evaluates the galaxy sample redshift distribution at redshift @z.
 *
 * Returns: the probability density at $z$, $P(z)$.
 */
gdouble
nc_galaxy_redshift_pop_ln_eval (NcGalaxyRedshiftPop *gsdrp, gdouble z)
{
  return NC_GALAXY_REDSHIFT_POP_GET_CLASS (gsdrp)->ln_eval (gsdrp, z);
}

