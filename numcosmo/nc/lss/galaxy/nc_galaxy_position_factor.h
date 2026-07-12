/***************************************************************************
 *            nc_galaxy_position_factor.h
 *
 *  Wed Jul 2 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_position_factor.h
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _NC_GALAXY_POSITION_FACTOR_H_
#define _NC_GALAXY_POSITION_FACTOR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/core/ncm_rng.h>
#include <numcosmo/ncm/core/ncm_util.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_POSITION_FACTOR           (nc_galaxy_position_factor_get_type ())
#define NC_TYPE_GALAXY_POSITION_FACTOR_DATA      (nc_galaxy_position_factor_data_get_type ())
#define NC_TYPE_GALAXY_POSITION_FACTOR_INTEGRAND (nc_galaxy_position_factor_integrand_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxyPositionFactor, nc_galaxy_position_factor, NC, GALAXY_POSITION_FACTOR, GObject)

typedef struct _NcGalaxyPositionFactorData NcGalaxyPositionFactorData;

/*
 * Integrand callback: the per-galaxy position density p(ra, dec | I) evaluated
 * at the galaxy's (ra, dec) stored in @data. The position is observed directly
 * (no scatter kernel, no marginalization), so there is no integration variable
 * here — the callback is evaluated at the fixed measured position.
 */
NCM_UTIL_DECLARE_CALLBACK (NcGalaxyPositionFactorIntegrand,
                           NC_GALAXY_POSITION_FACTOR_INTEGRAND,
                           nc_galaxy_position_factor_integrand,
                           gdouble,
                           NCM_UTIL_CALLBACK_ARGS (NcGalaxyPositionFactorData * data))

#define NC_GALAXY_POSITION_FACTOR_DATA(obj) ((NcGalaxyPositionFactorData *) (obj))

struct _NcGalaxyPositionFactorClass
{
  /*< private >*/
  GObjectClass parent_class;

  void (*data_init) (NcGalaxyPositionFactor *gspf, NcmMSet *mset, NcGalaxyPositionFactorData *data);
  void (*gen) (NcGalaxyPositionFactor *gspf, NcmMSet *mset, NcGalaxyPositionFactorData *data, NcmRNG *rng);
  void (*prepare) (NcGalaxyPositionFactor *gspf, NcmMSet *mset);
  NcGalaxyPositionFactorIntegrand *(*integ) (NcGalaxyPositionFactor *gspf, NcmMSet *mset, gboolean use_lnp);

  /* Factory-level change-detection: get_hash() returns an opaque value that
   * changes whenever prepare() refreshed something relevant (default:
   * a constant, "never changes"); update_data() unconditionally refreshes
   * one galaxy's cached state from what the last prepare() call resolved
   * (default: no-op). Callers (e.g. the cluster-WL orchestrator) call
   * prepare() once, compare get_hash() against their own last-seen value,
   * and call update_data() per galaxy only when it changed -- see
   * #NcGalaxyShapeFactor's analogous (but concrete, since that machinery is
   * identical across all its subclasses) radius/optzs/pop hashes for the
   * full rationale. Unlike Shape, these stay virtual: different concrete
   * Position schemes may cache completely different things. */
  guint64 (*get_hash) (NcGalaxyPositionFactor *gspf);
  void (*update_data) (NcGalaxyPositionFactor *gspf, NcGalaxyPositionFactorData *data);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[12];
};

/*
 * Per-galaxy data fragment. `ra`/`dec` are the inline base fields (the measured
 * sky position); `ldata` is the scheme's opaque fragment, packed to / unpacked
 * from an NcGalaxyWLObs via the fragment vtable below.
 */
struct _NcGalaxyPositionFactorData
{
  gdouble ra;
  gdouble dec;
  gpointer ldata;
  GDestroyNotify ldata_destroy;
  void (*ldata_read_row) (NcGalaxyPositionFactorData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_write_row) (NcGalaxyPositionFactorData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_required_columns) (NcGalaxyPositionFactorData *data, GList **columns);
  gatomicrefcount ref_count;
};

GType nc_galaxy_position_factor_data_get_type (void) G_GNUC_CONST;

NcGalaxyPositionFactorData *nc_galaxy_position_factor_data_ref (NcGalaxyPositionFactorData *data);
void nc_galaxy_position_factor_data_unref (NcGalaxyPositionFactorData *data);
void nc_galaxy_position_factor_data_read_row (NcGalaxyPositionFactorData *data, NcGalaxyWLObs *obs, const guint i);
void nc_galaxy_position_factor_data_write_row (NcGalaxyPositionFactorData *data, NcGalaxyWLObs *obs, const guint i);
GList *nc_galaxy_position_factor_data_required_columns (NcGalaxyPositionFactorData *data);

NcGalaxyPositionFactor *nc_galaxy_position_factor_ref (NcGalaxyPositionFactor *gspf);
void nc_galaxy_position_factor_free (NcGalaxyPositionFactor *gspf);
void nc_galaxy_position_factor_clear (NcGalaxyPositionFactor **gspf);

NcGalaxyPositionFactorData *nc_galaxy_position_factor_data_new (NcGalaxyPositionFactor *gspf, NcmMSet *mset);
void nc_galaxy_position_factor_gen (NcGalaxyPositionFactor *gspf, NcmMSet *mset, NcGalaxyPositionFactorData *data, NcmRNG *rng);
void nc_galaxy_position_factor_prepare (NcGalaxyPositionFactor *gspf, NcmMSet *mset);
NcGalaxyPositionFactorIntegrand *nc_galaxy_position_factor_integ (NcGalaxyPositionFactor *gspf, NcmMSet *mset, gboolean use_lnp);
guint64 nc_galaxy_position_factor_get_hash (NcGalaxyPositionFactor *gspf);
void nc_galaxy_position_factor_update_data (NcGalaxyPositionFactor *gspf, NcGalaxyPositionFactorData *data);

#define NC_GALAXY_POSITION_FACTOR_COL_RA "ra"
#define NC_GALAXY_POSITION_FACTOR_COL_DEC "dec"

G_END_DECLS

#endif /* _NC_GALAXY_POSITION_FACTOR_H_ */

