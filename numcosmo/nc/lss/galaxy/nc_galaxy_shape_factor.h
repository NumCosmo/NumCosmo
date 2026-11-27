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
#include <numcosmo/ncm/algebra/ncm_complex.h>
#include <numcosmo/ncm/model/ncm_mset.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_wl_obs.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_pop.h>
#include <numcosmo/nc/lss/wl/nc_wl_ellipticity.h>

G_BEGIN_DECLS

#define NC_TYPE_GALAXY_SHAPE_FACTOR       (nc_galaxy_shape_factor_get_type ())
#define NC_TYPE_GALAXY_SHAPE_FACTOR_DATA  (nc_galaxy_shape_factor_data_get_type ())
#define NC_GALAXY_SHAPE_FACTOR_ERROR      (nc_galaxy_shape_factor_error_quark ())

G_DECLARE_DERIVABLE_TYPE (NcGalaxyShapeFactor, nc_galaxy_shape_factor, NC, GALAXY_SHAPE_FACTOR, GObject)

/**
 * NcGalaxyShapeFactorError:
 * @NC_GALAXY_SHAPE_FACTOR_ERROR_ELLIP_CONV_MISMATCH: the observation's ellipticity convention does not match the factor's.
 *
 * Error codes returned by #NcGalaxyShapeFactor methods.
 */
typedef enum _NcGalaxyShapeFactorError /*< prefix=NC_GALAXY_SHAPE_FACTOR_ERROR_ELLIP_CONV >*/
{
  NC_GALAXY_SHAPE_FACTOR_ERROR_ELLIP_CONV_MISMATCH,
} NcGalaxyShapeFactorError;

GQuark nc_galaxy_shape_factor_error_quark (void);

typedef struct _NcGalaxyShapeFactorData NcGalaxyShapeFactorData;

#define NC_GALAXY_SHAPE_FACTOR_DATA(obj) ((NcGalaxyShapeFactorData *) (obj))

/*
 * NcGalaxyShapeFactor is the whole HSM measurement engine: it owns the
 * generative model (intrinsic shape from the NcGalaxyShapePop in the mset ->
 * deterministic shear map -> additive per-galaxy Gaussian pixel noise) and the
 * evaluation of the marginal at a given reduced shear. Subclasses supply ONLY the
 * intrinsic-ellipticity marginalization
 *
 *   P(e_obs | g) = \int_{|chi_I|<1} d^2chi_I P_pop(chi_I) N_2(e_obs - f_g(chi_I); std_noise^2)
 *
 * through eval_marginal / eval_ln_marginal: the evaluation strategy
 * (variance-add approximation, 2D quadrature, Laplace, truncated series) is
 * the only axis that varies. Both hooks receive the reduced shear g, with the
 * calibration bias (1+m) g + c applied by the caller, and the observed
 * ellipticity in the same frame, the frame of @data; the population model and
 * its per-galaxy fragment come via @pop and @data->pop_data.
 */
struct _NcGalaxyShapeFactorClass
{
  /*< private >*/
  GObjectClass parent_class;

  void (*data_init) (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data);
  void (*prepare) (NcGalaxyShapeFactor *gsf, NcmMSet *mset);
  gdouble (*eval_marginal) (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
  gdouble (*eval_ln_marginal) (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);

  /* One-line config description (default: type name + ellip-conv); override
   * to append and chain to the parent. */
  gchar *(*get_desc) (NcGalaxyShapeFactor *gsf);

  /*
   * Per-galaxy precomputation, called once per galaxy by #NcDataClusterWLFactor
   * after the population update and the bound fields, and before any redshift
   * quadrature evaluates the marginal of the galaxy. @z_max bounds the
   * galaxy's redshift support. The default does nothing: a subclass needs
   * this only when it caches something that depends on the galaxy and on
   * the models in @mset, which neither eval_marginal (no mset) nor prepare
   * (no data) can compute. What it builds must also be built lazily from
   * eval_marginal, since nothing guarantees this ran. */
  void (*data_prepare) (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data, const gdouble z_max);

  /* Prefetch hint for the subclass-owned per-galaxy state (@data->ldata and
   * what it points to), called by nc_galaxy_shape_factor_data_prefetch()
   * one stage at a time; see there for what each stage may assume is
   * already cached. It issues prefetches only, reads only pointers the
   * previous stage fetched, and changes no state. NULL (the default)
   * prefetches nothing. */
  void (*data_prefetch) (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, const guint stage);

  /*
   * The values of eval_marginal (want_log FALSE) or eval_ln_marginal (want_log TRUE) of one
   * galaxy at the n shears (g_1[j], g_2[j]), written to out[j]. The default calls them at
   * each shear; a subclass overrides it to evaluate the shears together. Used by
   * nc_galaxy_shape_factor_eval_shear_nodes(). @workspace comes from workspace_get, or is
   * NULL, and then the subclass takes its own for the call.
   */
  void (*eval_nodes) (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, gpointer workspace, const guint n, const gdouble *g_1, const gdouble *g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gboolean want_log, gdouble *out);

  /*
   * The values of eval_nodes at the n shears b + t[j] u of one line, written to out[j] (out
   * may be t), when the subclass computes them along the line; FALSE, with out unchanged
   * when out is not t, when it does not. The default returns FALSE. @workspace as in
   * eval_nodes.
   */
  gboolean (*eval_line) (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, gpointer workspace, const guint n, const gdouble *t, const gdouble u_1, const gdouble u_2, const gdouble b_1, const gdouble b_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gboolean want_log, gdouble *out);

  /*
   * The scratch space of eval_nodes and eval_line, taken once by a caller that evaluates
   * many galaxies and used by one thread at a time; workspace_return gives it back. The
   * default returns NULL and does nothing.
   */
  gpointer (*workspace_get) (NcGalaxyShapeFactor *gsf);
  void (*workspace_return) (NcGalaxyShapeFactor *gsf, gpointer workspace);

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[7];
};

/*
 * Per-galaxy shape data. The intrinsic-ellipticity fragment @pop_data belongs to the
 * NcGalaxyShapePop resolved from the mset at data_new time. @ldata is the
 * integration-method subclass scratch. The position and redshift of the galaxy are in
 * the fragments of the position and redshift factors, held beside this one by
 * #NcDataClusterWLFactor.
 */
struct _NcGalaxyShapeFactorData
{
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

  /*
   * The bound of the reduced shear of the galaxy over the prior box, set by
   * nc_galaxy_halo_wl_set_shape_bound() before nc_galaxy_shape_factor_data_prepare():
   * bound_ready, the corner of the box is set; bound_holds, only the halo center,
   * log10MDelta and shape-population parameters are free; bound_z_cl, the cluster
   * redshift; bound_R_min, the smallest projected radius of the galaxy over the center's
   * prior; bound_kappa and bound_g_t, the convergence and |g_t| at the corner
   * (bound_R_min, z_max), NaN when not evaluated. A galaxy whose bound was not set has
   * bound_ready FALSE.
   */
  gboolean bound_ready;
  gboolean bound_holds;
  gdouble bound_z_cl;
  gdouble bound_R_min;
  gdouble bound_kappa;
  gdouble bound_g_t;
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

NcGalaxyShapeFactor *nc_galaxy_shape_factor_ref (NcGalaxyShapeFactor *gsf);
void nc_galaxy_shape_factor_free (NcGalaxyShapeFactor *gsf);
void nc_galaxy_shape_factor_clear (NcGalaxyShapeFactor **gsf);

NcGalaxyWLObsEllipConv nc_galaxy_shape_factor_get_ellip_conv (NcGalaxyShapeFactor *gsf);

gboolean nc_galaxy_shape_factor_check_obs (NcGalaxyShapeFactor *gsf, NcGalaxyWLObs *obs, GError **error);

gchar *nc_galaxy_shape_factor_get_desc (NcGalaxyShapeFactor *gsf);

void nc_galaxy_shape_factor_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset);
void nc_galaxy_shape_factor_data_prepare (NcGalaxyShapeFactor *gsf, NcmMSet *mset, NcGalaxyShapeFactorData *data, const gdouble z_max);
void nc_galaxy_shape_factor_data_prefetch (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, const guint stage);
guint64 nc_galaxy_shape_factor_get_pop_hash (NcGalaxyShapeFactor *gsf);
void nc_galaxy_shape_factor_update_data_pop (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data);

void nc_galaxy_shape_factor_apply_shear (NcGalaxyShapeFactor *gsf, const NcmComplex *g, const NcmComplex *E, NcmComplex *E_obs);
void nc_galaxy_shape_factor_apply_shear_inv (NcGalaxyShapeFactor *gsf, const NcmComplex *g, const NcmComplex *E_obs, NcmComplex *E);
gdouble nc_galaxy_shape_factor_lndet_jac (NcGalaxyShapeFactor *gsf, const NcmComplex *g, const NcmComplex *E_obs);

NcGalaxyShapeFactorData *nc_galaxy_shape_factor_data_new (NcGalaxyShapeFactor *gsf, NcmMSet *mset);
void nc_galaxy_shape_factor_data_set (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2, const gdouble std_noise, const gdouble c1, const gdouble c2, const gdouble m, NcWLEllipticityFrame coord);
void nc_galaxy_shape_factor_data_get (NcGalaxyShapeFactor *gsf, NcGalaxyShapeFactorData *data, gdouble *epsilon_obs_1, gdouble *epsilon_obs_2, gdouble *std_noise, gdouble *c1, gdouble *c2, gdouble *m);

gboolean nc_galaxy_shape_factor_prepare_data_array (NcGalaxyShapeFactor *gsf, NcmMSet *mset, GPtrArray *data_array);
gdouble nc_galaxy_shape_factor_eval_shear (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gboolean want_log);
gpointer nc_galaxy_shape_factor_workspace_get (NcGalaxyShapeFactor *gsf);
void nc_galaxy_shape_factor_workspace_return (NcGalaxyShapeFactor *gsf, gpointer workspace);
void nc_galaxy_shape_factor_eval_shear_nodes (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, gpointer workspace, const guint n, const gdouble *g_1, const gdouble *g_2, const gboolean want_log, gdouble *out);
gboolean nc_galaxy_shape_factor_eval_shear_line (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, gpointer workspace, const guint n, const gdouble *t, const gdouble u_1, const gdouble u_2, const gdouble b_1, const gdouble b_2, const gboolean want_log, gdouble *out);
void nc_galaxy_shape_factor_gen_shear (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gboolean lensed, const gdouble g_1, const gdouble g_2, NcmRNG *rng);
void nc_galaxy_shape_factor_direct_estimate (NcGalaxyShapeFactor *gsf, NcmMSet *mset, GPtrArray *data_array, GArray *rotations, gdouble *gt, gdouble *gx, gdouble *sigma_t, gdouble *sigma_x, gdouble *rho);

gdouble nc_galaxy_shape_factor_eval_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
gdouble nc_galaxy_shape_factor_eval_ln_marginal (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const gdouble g_1, const gdouble g_2, const gdouble epsilon_obs_1, const gdouble epsilon_obs_2);
void nc_galaxy_shape_factor_eval_ln_marginal_array (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const NcmVector *g_1, const NcmVector *g_2, const NcmVector *epsilon_obs_1, const NcmVector *epsilon_obs_2, NcmVector *out);

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

