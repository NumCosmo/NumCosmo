/***************************************************************************
 *            nc_galaxy_shape_factor.h
 *
 *  Thu Jul 2 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_shape_factor.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
#ifndef _NC_GALAXY_SHAPE_FACTOR_H_
#define _NC_GALAXY_SHAPE_FACTOR_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/core/ncm_rng.h>
#include <numcosmo/ncm/core/ncm_util.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_pop.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_position_factor.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_redshift_factor.h>
#include <numcosmo/nc/lss/wl/nc_wl_ellipticity.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_FACTOR           (nc_galaxy_shape_factor_get_type ())
#define NC_TYPE_GALAXY_SHAPE_FACTOR_DATA      (nc_galaxy_shape_factor_data_get_type ())
#define NC_TYPE_GALAXY_SHAPE_FACTOR_INTEGRAND (nc_galaxy_shape_factor_integrand_get_type ())
#define NC_GALAXY_SHAPE_FACTOR_ERROR          (nc_galaxy_shape_factor_error_quark ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxyShapeFactor, nc_galaxy_shape_factor, NC, GALAXY_SHAPE_FACTOR, GObject)

/**
 * NcGalaxyShapeFactorError:
 * @NC_GALAXY_SHAPE_FACTOR_ERROR_ELLIP_CONV_MISMATCH: the observation's ellipticity convention does not match the factor's.
 *
 * Error codes returned by #NcGalaxyShapeFactor methods.
 */
typedef enum _NcGalaxyShapeFactorError
{
  NC_GALAXY_SHAPE_FACTOR_ERROR_ELLIP_CONV_MISMATCH,
} NcGalaxyShapeFactorError;

GQuark nc_galaxy_shape_factor_error_quark (void);

typedef struct _NcGalaxyShapeFactorData NcGalaxyShapeFactorData;

/*
 * Integrand callback: the per-galaxy shape likelihood factor
 * P(epsilon_obs | z, data) as a function of the (source) redshift z. The
 * ellipticity itself is never an integration variable at this level: the
 * intrinsic-ellipticity marginalization is performed *inside* the callback by
 * the integration-method subclass (see eval_marginal below), and the single
 * outer z-integral against the other pipeline factors belongs to the
 * orchestrator.
 */
NCM_UTIL_DECLARE_CALLBACK (NcGalaxyShapeFactorIntegrand,
                           NC_GALAXY_SHAPE_FACTOR_INTEGRAND,
                           nc_galaxy_shape_factor_integrand,
                           gdouble,
                           NCM_UTIL_CALLBACK_ARGS (const gdouble z, NcGalaxyShapeFactorData * data))

#define NC_GALAXY_SHAPE_FACTOR_DATA(obj) ((NcGalaxyShapeFactorData *) (obj))

/*
 * NcGalaxyShapeFactor is the whole HSM measurement engine: it owns the
 * generative model (intrinsic shape from the NcGalaxyShapePop in the mset ->
 * deterministic shear map with (1+m) g + c calibration bias -> additive
 * per-galaxy Gaussian pixel noise), the per-galaxy geometry caches and the
 * frame bookkeeping. Subclasses supply ONLY the intrinsic-ellipticity
 * marginalization
 *
 *   P(e_obs | g) = \int_{|chi_I|<1} d^2chi_I P_pop(chi_I) N_2(e_obs - f_g(chi_I); std_noise^2)
 *
 * through eval_marginal / eval_ln_marginal: the evaluation strategy
 * (variance-add approximation, 2D quadrature, Laplace, truncated series) is
 * the only axis that varies. Both hooks receive the reduced shear g (bias
 * already applied) and
 * the observed ellipticity in the SAME frame (the tangential/cross frame,
 * where the frame rotation has been handled by the engine); the population
 * model and its per-galaxy fragment come via @pop and @data->pop_data.
 */
struct _NcGalaxyShapeFactorClass
{
  /*< private >*/
  GObjectClass parent_class;

  void (*data_init) (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data);
  void (*prepare) (NcGalaxyShapeFactor *gsf, NcmMSet *mset);
  gdouble (*eval_marginal) (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
  gdouble (*eval_ln_marginal) (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);

  /* Human-readable one-line description of this scheme's own configuration
   * (default: the concrete type name plus :ellip-conv, the one property
   * shared by every subclass), same rationale as #NcGalaxyPositionFactor's
   * analogous vfunc (see there) -- override to append scheme-specific
   * configuration (e.g. SeriesLensed's truncation order), chaining to the
   * parent class implementation rather than duplicating the ellip-conv
   * part. */
  gchar *(*get_desc) (NcGalaxyShapeFactor *gsf);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[13];
};

/*
 * Per-galaxy data. Upstream fragments are held as flat references: the shape
 * factor needs ra/dec for the lens geometry and z at generation time. The
 * intrinsic-ellipticity fragment @pop_data belongs to the NcGalaxyShapePop
 * resolved from the mset at data_new time. @cdata holds the engine-owned
 * geometry caches (opaque, defined in the .c); @ldata is the
 * integration-method subclass scratch.
 */
struct _NcGalaxyShapeFactorData
{
  NcGalaxyPositionFactorData *pos_data;
  NcGalaxyRedshiftFactorData *z_data;
  NcGalaxyShapePopData *pop_data;

  /* Handedness frame in which this galaxy's ellipticity components are
   * expressed (see #NcWLEllipticityFrame). Sky positions are always RA/Dec, so
   * only the ellipticity basis is selected by this field. */
  NcWLEllipticityFrame coord;
  gdouble epsilon_int_1;
  gdouble epsilon_int_2;
  gdouble epsilon_obs_1;
  gdouble epsilon_obs_2;
  gdouble std_noise;
  gdouble c1;
  gdouble c2;
  gdouble m;
  gpointer cdata;
  gpointer ldata;
  GDestroyNotify ldata_destroy;
  void (*ldata_read_row) (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_write_row) (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i);
  void (*ldata_required_columns) (NcGalaxyShapeFactorData *data, GList **columns);
  gatomicrefcount ref_count;
};

GType nc_galaxy_shape_factor_data_get_type (void) G_GNUC_CONST;

NcGalaxyShapeFactorData *nc_galaxy_shape_factor_data_ref (NcGalaxyShapeFactorData *data);
void nc_galaxy_shape_factor_data_unref (NcGalaxyShapeFactorData *data);
void nc_galaxy_shape_factor_data_read_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i);
void nc_galaxy_shape_factor_data_write_row (NcGalaxyShapeFactorData *data, NcGalaxyWLObs *obs, const guint i);
GList *nc_galaxy_shape_factor_data_required_columns (NcGalaxyShapeFactorData *data);
gdouble nc_galaxy_shape_factor_data_get_radius (NcGalaxyShapeFactorData *data);

NcGalaxyShapeFactor *nc_galaxy_shape_factor_ref (NcGalaxyShapeFactor *gsf);
void nc_galaxy_shape_factor_free (NcGalaxyShapeFactor *gsf);
void nc_galaxy_shape_factor_clear (NcGalaxyShapeFactor **gsf);

NcGalaxyWLObsEllipConv nc_galaxy_shape_factor_get_ellip_conv (NcGalaxyShapeFactor *gsf);

gboolean nc_galaxy_shape_factor_check_obs (NcGalaxyShapeFactor *gsf, NcGalaxyWLObs *obs, GError **error);

gchar *nc_galaxy_shape_factor_get_desc (NcGalaxyShapeFactor *gsf);

void nc_galaxy_shape_factor_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset);
guint64 nc_galaxy_shape_factor_get_radius_hash (NcGalaxyShapeFactor *gsf);
guint64 nc_galaxy_shape_factor_get_optzs_hash (NcGalaxyShapeFactor *gsf);
guint64 nc_galaxy_shape_factor_get_crit_hash (NcGalaxyShapeFactor *gsf);
guint64 nc_galaxy_shape_factor_get_pop_hash (NcGalaxyShapeFactor *gsf);
void nc_galaxy_shape_factor_update_data_radius (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data);
void nc_galaxy_shape_factor_update_data_optzs (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data);
void nc_galaxy_shape_factor_update_data_pop (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data);
void nc_galaxy_shape_factor_update_data_at_nodes_crit (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, const NcmVector *z_nodes);
void nc_galaxy_shape_factor_update_data_at_nodes_sigma (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data);

void nc_galaxy_shape_factor_apply_shear (NcGalaxyShapeFactor *gsf, const NcmComplex *g, const NcmComplex *E, NcmComplex *E_obs);
void nc_galaxy_shape_factor_apply_shear_inv (NcGalaxyShapeFactor *gsf, const NcmComplex *g, const NcmComplex *E_obs, NcmComplex *E);
gdouble nc_galaxy_shape_factor_lndet_jac (NcGalaxyShapeFactor *gsf, const NcmComplex *g, const NcmComplex *E_obs);

NcGalaxyShapeFactorData *nc_galaxy_shape_factor_data_new (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyPositionFactorData *pos_data, NcGalaxyRedshiftFactorData *z_data);
void nc_galaxy_shape_factor_data_set (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble std_noise, const gdouble c1, const gdouble c2, const gdouble m, NcWLEllipticityFrame coord);
void nc_galaxy_shape_factor_data_get (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, gdouble *epsilon_obs_1, gdouble *epsilon_obs_2, gdouble *std_noise, gdouble *c1, gdouble *c2, gdouble *m);

void nc_galaxy_shape_factor_gen (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data, NcmRNG *rng);
NcGalaxyShapeFactorIntegrand *nc_galaxy_shape_factor_integ (NcGalaxyShapeFactor *gsf, NcmMSet *mset, gboolean use_lnp);
gboolean nc_galaxy_shape_factor_prepare_data_array (NcGalaxyShapeFactor *gsf, NcmMSet *mset, GPtrArray *data_array, gboolean update_radius, gboolean update_optzs);
gboolean nc_galaxy_shape_factor_prepare_data_array_at_nodes (NcGalaxyShapeFactor *gsf, NcmMSet *mset, GPtrArray *data_array, const GPtrArray *z_nodes_per_galaxy, gboolean update_radius, gboolean update_crit, gboolean update_sigma);
void nc_galaxy_shape_factor_eval_at_nodes (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data, const NcmVector *z_nodes, NcmVector *out);
void nc_galaxy_shape_factor_direct_estimate (NcGalaxyShapeFactor *gsf, NcmMSet *mset, GPtrArray *data_array, gdouble *gt, gdouble *gx, gdouble *sigma_t, gdouble *sigma_x, gdouble *rho);

gdouble nc_galaxy_shape_factor_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
gdouble nc_galaxy_shape_factor_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);

#define NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_INT_1 "epsilon_int_1"
#define NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_INT_2 "epsilon_int_2"
#define NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_1 "epsilon_obs_1"
#define NC_GALAXY_SHAPE_FACTOR_COL_EPSILON_OBS_2 "epsilon_obs_2"
#define NC_GALAXY_SHAPE_FACTOR_COL_STD_NOISE "std_noise"
#define NC_GALAXY_SHAPE_FACTOR_COL_C1 "c1"
#define NC_GALAXY_SHAPE_FACTOR_COL_C2 "c2"
#define NC_GALAXY_SHAPE_FACTOR_COL_M "m"

G_END_DECLS

#endif /* _NC_GALAXY_SHAPE_FACTOR_H_ */

