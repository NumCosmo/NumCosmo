/***************************************************************************
 *            test_nc_galaxy_position_factor.c
 *
 *  Wed Jul 22 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
#include <glib.h>
#include <glib-object.h>

/*
 * A minimal test-only NcGalaxyPositionFactor subclass overriding only the
 * three vfuncs with no default (data_init/gen/integ -- unimplemented aborts
 * unless overridden, see nc_galaxy_position_factor.c's LCOV_EXCL'd stubs),
 * left otherwise bare so get_desc() falls through to the abstract base's own
 * default implementation. NcGalaxyPositionFactorFlat, the only real
 * subclass, overrides get_desc() with its own footprint-reporting version,
 * so the base default is otherwise unreachable through any production code
 * path -- exactly the kind of gap this test-double idiom exists for (see
 * e.g. ncm_model_test.c).
 */

#define TEST_TYPE_GALAXY_POSITION_FACTOR_BARE (test_galaxy_position_factor_bare_get_type ())

G_DECLARE_FINAL_TYPE (TestGalaxyPositionFactorBare, test_galaxy_position_factor_bare, TEST, GALAXY_POSITION_FACTOR_BARE, NcGalaxyPositionFactor)

struct _TestGalaxyPositionFactorBare
{
  NcGalaxyPositionFactor parent_instance;
};

G_DEFINE_TYPE (TestGalaxyPositionFactorBare, test_galaxy_position_factor_bare, NC_TYPE_GALAXY_POSITION_FACTOR)

static void
test_galaxy_position_factor_bare_init (TestGalaxyPositionFactorBare *self)
{
}

static void
_test_galaxy_position_factor_bare_data_init (NcGalaxyPositionFactor *gspf, NcmMSet *mset, NcGalaxyPositionFactorData *data)
{
}

static void
_test_galaxy_position_factor_bare_gen (NcGalaxyPositionFactor *gspf, NcmMSet *mset, NcGalaxyPositionFactorData *data, NcmRNG *rng)
{
}

static NcGalaxyPositionFactorIntegrand *
_test_galaxy_position_factor_bare_integ (NcGalaxyPositionFactor *gspf, NcmMSet *mset, gboolean use_lnp)
{
  return NULL;
}

static void
test_galaxy_position_factor_bare_class_init (TestGalaxyPositionFactorBareClass *klass)
{
  NcGalaxyPositionFactorClass *position_factor_class = NC_GALAXY_POSITION_FACTOR_CLASS (klass);

  position_factor_class->data_init = &_test_galaxy_position_factor_bare_data_init;
  position_factor_class->gen       = &_test_galaxy_position_factor_bare_gen;
  position_factor_class->integ     = &_test_galaxy_position_factor_bare_integ;
}

void test_nc_galaxy_position_factor_default_get_desc (void);

int
main (int argc, char *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/nc/galaxy_position_factor/default_get_desc", test_nc_galaxy_position_factor_default_get_desc);

  g_test_run ();
}

void
test_nc_galaxy_position_factor_default_get_desc (void)
{
  NcGalaxyPositionFactor *gspf = g_object_new (TEST_TYPE_GALAXY_POSITION_FACTOR_BARE, NULL);
  gchar *desc                  = nc_galaxy_position_factor_get_desc (gspf);

  g_assert_cmpstr (desc, ==, G_OBJECT_TYPE_NAME (gspf));

  g_free (desc);
  g_object_unref (gspf);
}

