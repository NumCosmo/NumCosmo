/***************************************************************************
 *            nc_galaxy_redshift_factor.h
 *
 *  Tue Jul 1 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_factor.h
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
#ifndef _NC_GALAXY_REDSHIFT_FACTOR_H_
#define _NC_GALAXY_REDSHIFT_FACTOR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/core/ncm_rng.h>
#include <numcosmo/ncm/core/ncm_util.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/ncm/integration/ncm_integrate.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_REDSHIFT_FACTOR           (nc_galaxy_redshift_factor_get_type ())
#define NC_TYPE_GALAXY_REDSHIFT_FACTOR_DATA      (nc_galaxy_redshift_factor_data_get_type ())
#define NC_TYPE_GALAXY_REDSHIFT_FACTOR_INTEGRAND (nc_galaxy_redshift_factor_integrand_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxyRedshiftFactor, nc_galaxy_redshift_factor, NC, GALAXY_REDSHIFT_FACTOR, GObject)

typedef struct _NcGalaxyRedshiftFactorData NcGalaxyRedshiftFactorData;

/*
 * Integrand callback: the per-galaxy JOINT density p(z_phot, z | I) as a
 * function of the true redshift z (the calculator NEVER integrates z itself;
 * the orchestrator integrates this against every other z-dependent factor).
 */
NCM_UTIL_DECLARE_CALLBACK (NcGalaxyRedshiftFactorIntegrand,
                           NC_GALAXY_REDSHIFT_FACTOR_INTEGRAND,
                           nc_galaxy_redshift_factor_integrand,
                           gdouble,
                           NCM_UTIL_CALLBACK_ARGS (const gdouble z, NcGalaxyRedshiftFactorData * data))

#define NC_GALAXY_REDSHIFT_FACTOR_DATA(obj) ((NcGalaxyRedshiftFactorData *) (obj))

struct _NcGalaxyRedshiftFactorClass
{
  /*< private >*/
  GObjectClass parent_class;

  void (*data_init) (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data);
  void (*gen) (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng);
  gboolean (*gen1) (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng);
  void (*prepare) (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset);
  NcGalaxyRedshiftFactorIntegrand *(*integ) (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, gboolean use_lnp);
  void (*get_integ_lim) (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble *z_min, gdouble *z_max);
  gdouble (*norm) (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data);
  NcmIntegralFixed *(*make_fixed_nodes) (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble z_lo, gdouble z_hi, guint n_nodes, guint rule_n);

  /* Factory-level change-detection, same rationale/contract as
   * #NcGalaxyPositionFactor's analogous pair (see there): get_hash()
   * reflects what the last prepare() call resolved (default: a constant);
   * update_data() unconditionally refreshes one galaxy's cached state
   * (default: no-op). */
  guint64 (*get_hash) (NcGalaxyRedshiftFactor *gsdr);
  void (*update_data) (NcGalaxyRedshiftFactor *gsdr, NcGalaxyRedshiftFactorData *data);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[8];
};

/*
 * Per-galaxy data fragment. `z` is the inline base field (the true redshift,
 * sampled or the integration variable); `ldata` is the scheme's opaque fragment
 * (embedding the Observable model's own {z_phot, sigma0} fragment for Composed),
 * packed to / unpacked from an NcGalaxyWLObs via the fragment vtable below.
 */
struct _NcGalaxyRedshiftFactorData
{
  gdouble z;
  gpointer ldata;
  GDestroyNotify ldata_destroy;
  void (*ldata_read_row) (NcGalaxyRedshiftFactorData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_write_row) (NcGalaxyRedshiftFactorData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_required_columns) (NcGalaxyRedshiftFactorData *data, GList **columns);
  gatomicrefcount ref_count;
};

GType nc_galaxy_redshift_factor_data_get_type (void) G_GNUC_CONST;

NcGalaxyRedshiftFactorData *nc_galaxy_redshift_factor_data_ref (NcGalaxyRedshiftFactorData *data);
void nc_galaxy_redshift_factor_data_unref (NcGalaxyRedshiftFactorData *data);
void nc_galaxy_redshift_factor_data_read_row (NcGalaxyRedshiftFactorData *data, NcGalaxyWLObs *obs, const guint i);
void nc_galaxy_redshift_factor_data_write_row (NcGalaxyRedshiftFactorData *data, NcGalaxyWLObs *obs, const guint i);
GList *nc_galaxy_redshift_factor_data_required_columns (NcGalaxyRedshiftFactorData *data);

NcGalaxyRedshiftFactor *nc_galaxy_redshift_factor_ref (NcGalaxyRedshiftFactor *gsdr);
void nc_galaxy_redshift_factor_free (NcGalaxyRedshiftFactor *gsdr);
void nc_galaxy_redshift_factor_clear (NcGalaxyRedshiftFactor **gsdr);

NcGalaxyRedshiftFactorData *nc_galaxy_redshift_factor_data_new (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset);
void nc_galaxy_redshift_factor_gen (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng);
gboolean nc_galaxy_redshift_factor_gen1 (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, NcmRNG *rng);
void nc_galaxy_redshift_factor_prepare (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset);
NcGalaxyRedshiftFactorIntegrand *nc_galaxy_redshift_factor_integ (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, gboolean use_lnp);
void nc_galaxy_redshift_factor_get_integ_lim (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble *z_min, gdouble *z_max);
gdouble nc_galaxy_redshift_factor_norm (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data);
NcmIntegralFixed *nc_galaxy_redshift_factor_make_fixed_nodes (NcGalaxyRedshiftFactor *gsdr, NcmMSet *mset, NcGalaxyRedshiftFactorData *data, gdouble z_lo, gdouble z_hi, guint n_nodes, guint rule_n);
guint64 nc_galaxy_redshift_factor_get_hash (NcGalaxyRedshiftFactor *gsdr);
void nc_galaxy_redshift_factor_update_data (NcGalaxyRedshiftFactor *gsdr, NcGalaxyRedshiftFactorData *data);

#define NC_GALAXY_REDSHIFT_FACTOR_COL_Z "z"

G_END_DECLS

#endif /* _NC_GALAXY_REDSHIFT_FACTOR_H_ */

