/***************************************************************************
 *            nc_galaxy_wl_obs.c
 *
 *  Tue Jul 16 06:43:45 2024
 *  Copyright  2024 Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_wl_obs.c
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * NcGalaxyWLObs:
 *
 * Galaxy weak lensing observation data.
 *
 * A class to store galaxy weak lensing observation data and information about its
 * coordinate system.
 *
 * A catalog may also carry the interpolants of $\ln P(\epsilon_\mathrm{obs} \mid g)$ of its
 * galaxies with the record of what built them, a #NcGalaxyWLObsShapeInterp
 * (#NcGalaxyWLObs:shape-interp), either in memory or backed by a file. They describe the rows
 * whose data hash (ncm_catalog_get_data_hash()) they carry. A catalog that receives them
 * through the property, as by deserialization, compares that hash with the one of its rows
 * at the first call of nc_galaxy_wl_obs_peek_shape_interp() or
 * nc_galaxy_wl_obs_peek_shape_interp_source(), and aborts with an error when they differ.
 * From then on they describe the rows as they are, by the key of #NcmCatalog
 * (ncm_catalog_key_is_current()): a later change of the rows makes the two functions return
 * %NULL and the property return %NULL, so a serialized catalog never carries interpolants of
 * other rows.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"
#include "stdarg.h"
#include "nc/lss/galaxy/nc_galaxy_wl_obs.h"
#include "ncm/algebra/ncm_matrix.h"
#include "ncm/core/ncm_obj_array.h"
#include "ncm/spline/ncm_spline.h"
#include "ncm/algebra/ncm_vector.h"
#include "ncm/algebra/ncm_spectral_tensor.h"
#include "nc_enum_types.h"
#include "nc/lss/galaxy/nc_galaxy_wl_obs.h"
#include "nc/background/nc_hicosmo.h"
#include "nc/lss/halo/nc_halo_density_profile.h"
#include "nc/lss/halo/nc_halo_position.h"
#include "nc/lss/wl/nc_wl_surface_mass_density.h"
#include "nc/data/nc_data_download_priv.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/core/ncm_util.h"

/* *INDENT-OFF* */
G_DEFINE_QUARK (nc-galaxy-wl-obs-error, nc_galaxy_wl_obs_error) 
/* *INDENT-ON* */

struct _NcGalaxyWLObs
{
  /*< private >*/
  NcmCatalog parent_instance;
};

struct _NcGalaxyWLObsPrivate
{
  NcmObjDictInt *pz;
  NcWLEllipticityFrame coord;
  NcGalaxyWLObsEllipConv ellip_conv;
  NcGalaxyWLObsShapeInterp *shape_interp;
  guint64 shape_interp_key; /* catalog key of the rows the interpolants describe */
  NcGalaxyWLObsChebGrid *cheb_grid;
  guint64 cheb_grid_key; /* catalog key of the rows the grid describes */
};

enum
{
  PROP_PZ = 1000,
  PROP_COORD,
  PROP_ELLIP_CONV,
  PROP_SHAPE_INTERP,
  PROP_CHEB_GRID,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyWLObs, nc_galaxy_wl_obs, NCM_TYPE_CATALOG)

static void
nc_galaxy_wl_obs_init (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  self->pz               = ncm_obj_dict_int_new ();
  self->coord            = NC_GALAXY_WL_OBS_COORD_EUCLIDEAN;
  self->ellip_conv       = 0;
  self->shape_interp     = NULL;
  self->shape_interp_key = NCM_CATALOG_KEY_PENDING;
  self->cheb_grid        = NULL;
  self->cheb_grid_key    = NCM_CATALOG_KEY_PENDING;
}

static gboolean _nc_galaxy_wl_obs_shape_interp_current (NcGalaxyWLObs *obs);
static gboolean _nc_galaxy_wl_obs_cheb_grid_current (NcGalaxyWLObs *obs);

static void
_nc_galaxy_wl_obs_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLObs *obs                = NC_GALAXY_WL_OBS (object);
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  switch (prop_id)
  {
    case PROP_PZ:
      g_value_set_boxed (value, self->pz);
      break;
    case PROP_COORD:
      g_value_set_enum (value, nc_galaxy_wl_obs_get_coord (obs));
      break;
    case PROP_ELLIP_CONV:
      g_value_set_enum (value, nc_galaxy_wl_obs_get_ellip_conv (obs));
      break;
    case PROP_SHAPE_INTERP:

      if (_nc_galaxy_wl_obs_shape_interp_current (obs))
        g_value_set_object (value, self->shape_interp);

      break;
    case PROP_CHEB_GRID:

      if (_nc_galaxy_wl_obs_cheb_grid_current (obs))
        g_value_set_object (value, self->cheb_grid);

      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void _nc_galaxy_wl_obs_take_pz (NcGalaxyWLObsPrivate *self, NcmObjDictInt *pz);
static void _nc_galaxy_wl_obs_set_ellip_conv (NcGalaxyWLObs *obs, NcGalaxyWLObsEllipConv type);

static void
_nc_galaxy_wl_obs_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLObs *obs         = NC_GALAXY_WL_OBS (object);
  NcGalaxyWLObsPrivate *self = nc_galaxy_wl_obs_get_instance_private (obs);

  switch (prop_id)
  {
    case PROP_PZ:
      _nc_galaxy_wl_obs_take_pz (self, g_value_dup_boxed (value));
      break;
    case PROP_COORD:
      nc_galaxy_wl_obs_set_coord (obs, g_value_get_enum (value));
      break;
    case PROP_ELLIP_CONV:
      _nc_galaxy_wl_obs_set_ellip_conv (obs, g_value_get_enum (value));
      break;
    case PROP_SHAPE_INTERP:
      nc_galaxy_wl_obs_shape_interp_clear (&self->shape_interp);
      self->shape_interp     = g_value_dup_object (value);
      self->shape_interp_key = NCM_CATALOG_KEY_PENDING;
      break;
    case PROP_CHEB_GRID:
      nc_galaxy_wl_obs_cheb_grid_clear (&self->cheb_grid);
      self->cheb_grid     = g_value_dup_object (value);
      self->cheb_grid_key = NCM_CATALOG_KEY_PENDING;
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_wl_obs_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_galaxy_wl_obs_parent_class)->constructed (object);
}

static void
_nc_galaxy_wl_obs_dispose (GObject *object)
{
  NcGalaxyWLObs *obs                = NC_GALAXY_WL_OBS (object);
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  ncm_obj_dict_int_clear (&self->pz);
  nc_galaxy_wl_obs_shape_interp_clear (&self->shape_interp);
  nc_galaxy_wl_obs_cheb_grid_clear (&self->cheb_grid);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_obs_parent_class)->dispose (object);
}

/**
 * _nc_galaxy_wl_obs_finalize:
 * @object: a #GObject.
 *
 * Finalizes the #NcGalaxyWLObs object.
 *
 */
static void
_nc_galaxy_wl_obs_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_obs_parent_class)->finalize (object);
}

static void
nc_galaxy_wl_obs_class_init (NcGalaxyWLObsClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->get_property = &_nc_galaxy_wl_obs_get_property;
  object_class->set_property = &_nc_galaxy_wl_obs_set_property;
  object_class->constructed  = &_nc_galaxy_wl_obs_constructed;
  object_class->dispose      = &_nc_galaxy_wl_obs_dispose;
  object_class->finalize     = &_nc_galaxy_wl_obs_finalize;

  /**
   * NcGalaxyWLObs:pz:
   *
   * The redshift densities of the galaxies: a #NcGalaxyPz for each row that has one,
   * keyed by the row.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_PZ,
                                   g_param_spec_boxed ("pz",
                                                       "P(z)",
                                                       "Redshift density of each galaxy",
                                                       NCM_TYPE_OBJ_DICT_INT,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObs:coord:
   *
   * Coordinate system used to store the data.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_COORD,
                                   g_param_spec_enum ("coord",
                                                      "Coordinate system",
                                                      "Coordinate system used to store the data",
                                                      NC_TYPE_GALAXY_WL_OBS_COORD,
                                                      NC_GALAXY_WL_OBS_COORD_EUCLIDEAN,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

/**
 * NcGalaxyWLObs:ellip-conv:
 *
 * Ellipticity convention.
 *
 */
  g_object_class_install_property (object_class,
                                   PROP_ELLIP_CONV,
                                   g_param_spec_enum ("ellip-conv",
                                                      "Ellipticity convention",
                                                      "Weak lensing observables ellipticity convention",
                                                      NC_TYPE_GALAXY_WL_OBS_ELLIP_CONV,
                                                      NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObs:shape-interp:
   *
   * The interpolants of $\ln P(\epsilon_\mathrm{obs} \mid g)$ of the galaxies and the
   * record of what built them; %NULL when there are none or when the rows changed after
   * they were set.
   */
  g_object_class_install_property (object_class,
                                   PROP_SHAPE_INTERP,
                                   g_param_spec_object ("shape-interp",
                                                        "Shape interpolants",
                                                        "Interpolants of ln P in g of the galaxies",
                                                        NC_TYPE_GALAXY_WL_OBS_SHAPE_INTERP,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  /**
   * NcGalaxyWLObs:cheb-grid:
   *
   * The CHEB_WEIGHTS redshift grid of the galaxies and the record of what built it; %NULL
   * when there is none, after the values of the rows changed or after a p(z) was set.
   */
  g_object_class_install_property (object_class,
                                   PROP_CHEB_GRID,
                                   g_param_spec_object ("cheb-grid",
                                                        "Chebyshev grid",
                                                        "CHEB_WEIGHTS redshift grid of the galaxies",
                                                        NC_TYPE_GALAXY_WL_OBS_CHEB_GRID,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));
}

/*
 * Whether the interpolants and their record describe the current rows. Interpolants set
 * through the property are checked, at the first call, against the data hash of the rows;
 * a mismatch aborts.
 */
static gboolean
_nc_galaxy_wl_obs_shape_interp_current (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  if ((self->shape_interp == NULL) || (nc_galaxy_wl_obs_shape_interp_peek_source (self->shape_interp) == NULL))
    return FALSE;

  if (self->shape_interp_key == NCM_CATALOG_KEY_PENDING)
  {
    const gchar *data_hash = nc_galaxy_wl_obs_shape_interp_peek_data_hash (self->shape_interp);
    gchar *rows_hash       = ncm_catalog_get_data_hash (NCM_CATALOG (obs));

    if (g_strcmp0 (data_hash, rows_hash) != 0)
      g_error ("nc_galaxy_wl_obs: the interpolants describe the rows of data hash %s, not these rows of data hash %s.",
               data_hash, rows_hash);

    g_free (rows_hash);
  }

  return ncm_catalog_key_is_current (NCM_CATALOG (obs), &self->shape_interp_key);
}

/*
 * Whether the grid describes the current rows. A grid set through the property is checked,
 * at the first call, against the data hash of the rows; a mismatch aborts. After a change
 * of the rows the data hash is compared again: rows written back with the same values, as
 * by nc_data_cluster_wl_factor_write_shape_interp(), keep the grid; other values drop it.
 */
static gboolean
_nc_galaxy_wl_obs_cheb_grid_current (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  if (self->cheb_grid == NULL)
    return FALSE;

  if (self->cheb_grid_key == NCM_CATALOG_KEY_PENDING)
  {
    const gchar *data_hash = nc_galaxy_wl_obs_cheb_grid_peek_data_hash (self->cheb_grid);
    gchar *rows_hash       = ncm_catalog_get_data_hash (NCM_CATALOG (obs));

    if (g_strcmp0 (data_hash, rows_hash) != 0)
      g_error ("nc_galaxy_wl_obs: the Chebyshev grid describes the rows of data hash %s, not these rows of data hash %s.",
               data_hash, rows_hash);

    g_free (rows_hash);
  }

  if (!ncm_catalog_key_is_current (NCM_CATALOG (obs), &self->cheb_grid_key))
  {
    const gchar *data_hash = nc_galaxy_wl_obs_cheb_grid_peek_data_hash (self->cheb_grid);
    gchar *rows_hash       = ncm_catalog_get_data_hash (NCM_CATALOG (obs));
    const gboolean same    = (g_strcmp0 (data_hash, rows_hash) == 0);

    g_free (rows_hash);

    if (!same)
    {
      nc_galaxy_wl_obs_cheb_grid_clear (&self->cheb_grid);

      return FALSE;
    }

    self->cheb_grid_key = ncm_catalog_get_key (NCM_CATALOG (obs));
  }

  return TRUE;
}

void
_nc_galaxy_wl_obs_take_pz (NcGalaxyWLObsPrivate *self, NcmObjDictInt *pz)
{
  guint i;

  ncm_obj_dict_int_clear (&self->pz);
  self->pz = pz;

  if (self->pz != NULL)
  {
    GArray *keys  = ncm_obj_dict_int_keys (self->pz);
    const guint n = keys->len;

    for (i = 0; i < n; i++)
    {
      const gint key = g_array_index (keys, gint, i);
      GObject *obj   = ncm_obj_dict_int_peek (self->pz, key);

      /* Catalogs written by NumCosmo 0.x before NcGalaxyPz: this check is dropped in 1.0 */
      if (NCM_IS_SPLINE (obj))
        g_error ("nc_galaxy_wl_obs: the catalog holds the p(z) of row %d as a %s, the format written before NcGalaxyPz. "
                 "Update the file with tools/nc_galaxy_pz_update_catalog.py. A catalog downloaded by NumCosmo "
                 "(wl_obs_*.gvar in the NumCosmo data directory) can instead be removed, so that it is downloaded again.",
                 key, G_OBJECT_TYPE_NAME (obj));

      g_assert (NC_IS_GALAXY_PZ (obj));
    }

    g_array_unref (keys);
  }
  else
  {
    self->pz = ncm_obj_dict_int_new ();
  }
}

/**
 * nc_galaxy_wl_obs_new:
 * @ellip_conv: the ellipticity convention
 * @coord: the coordinate system used to store the data
 * @nrows: the number of data rows
 * @col_names: the column names
 * @col_types: (array length=n_types): the type of each column of @col_names
 * @n_types: the number of entries of @col_types, the number of @col_names
 *
 * Creates a new #NcGalaxyWLObs object with @nrows rows and the typed columns
 * @col_names, as ncm_catalog_new_full(). The columns are read and written with the
 * #NcmCatalog functions.
 *
 * Returns: (transfer full): a new #NcGalaxyWLObs object.
 */
NcGalaxyWLObs *
nc_galaxy_wl_obs_new (NcGalaxyWLObsEllipConv ellip_conv, NcWLEllipticityFrame coord, guint nrows, GStrv col_names,
                      const NcmCatalogColType *col_types, guint n_types)
{
  GArray *types = g_array_sized_new (FALSE, TRUE, sizeof (guint32), n_types);
  NcGalaxyWLObs *obs;
  GVariant *var;
  guint i;

  for (i = 0; i < n_types; i++)
  {
    const guint32 t = col_types[i];

    g_array_append_val (types, t);
  }

  var = ncm_cfg_array_to_variant (types, G_VARIANT_TYPE ("u"));
  obs = g_object_new (NC_TYPE_GALAXY_WL_OBS,
                      "ellip-conv", ellip_conv,
                      "columns", col_names,
                      "col-types", var,
                      "coord", coord,
                      "len", nrows,
                      NULL);

  g_variant_unref (var);
  g_array_unref (types);

  return obs;
}

/* Filenames of the curated Subaru HSC-SSP PDR1 catalogs in the NumCosmo
 * datafile-release-v1.0.0 GitHub release, indexed by #NcGalaxyWLObsCatalogId. */
static const gchar *_nc_galaxy_wl_obs_catalog_files[] = {
  "wl_obs_HWL16a-002.gvar",
  "wl_obs_HWL16a-007.gvar",
  "wl_obs_HWL16a-060.gvar",
  "wl_obs_HWL16a-064.gvar",
  "wl_obs_HWL16a-094.gvar",
};

/**
 * nc_galaxy_wl_obs_catalog_id_get_filename:
 * @id: a #NcGalaxyWLObsCatalogId
 *
 * Downloads (if not already cached) and returns the local filename of the
 * catalog identified by @id. The file is fetched from the NumCosmo
 * datafile-release-v1.0.0 GitHub release into the NumCosmo user data directory
 * (see ncm_cfg_get_fullpath_base()) the first time it is requested; later
 * calls reuse the cached copy.
 *
 * Returns: (transfer full): Full path for the catalog file.
 */
gchar *
nc_galaxy_wl_obs_catalog_id_get_filename (NcGalaxyWLObsCatalogId id)
{
  g_return_val_if_fail (id < NC_GALAXY_WL_OBS_CATALOG_LEN, NULL);

  {
    const gchar *filename = _nc_galaxy_wl_obs_catalog_files[id];
    gchar *full_filename  = ncm_cfg_get_fullpath (filename);

    if (!g_file_test (full_filename, G_FILE_TEST_EXISTS))
    {
      gchar *url_str = g_strdup_printf ("https://github.com/NumCosmo/NumCosmo/releases/download/datafile-release-v1.0.0/%s", filename);
      gchar *lockdir = NULL;

      /* Serialize with any other process fetching the same catalog, and let it
       * publish by rename: writing straight to full_filename let a concurrent
       * reader open a half-written catalog, and an interrupted fetch left a
       * truncated one that every later run took for complete. */
      if (_nc_data_download_lock (full_filename, full_filename, 900, &lockdir))
      {
        _nc_data_download_file (url_str, full_filename, filename, lockdir);
        _nc_data_download_unlock (lockdir);
      }

      g_free (url_str);
    }

    return full_filename;
  }
}

/**
 * nc_galaxy_wl_obs_new_from_catalog_id: (constructor)
 * @id: a #NcGalaxyWLObsCatalogId
 *
 * Downloads (if necessary) and loads the catalog identified by @id.
 *
 * Returns: (transfer full): a new #NcGalaxyWLObs object.
 */
NcGalaxyWLObs *
nc_galaxy_wl_obs_new_from_catalog_id (NcGalaxyWLObsCatalogId id)
{
  gchar *full_filename = nc_galaxy_wl_obs_catalog_id_get_filename (id);
  NcmSerialize *ser    = ncm_serialize_new (NCM_SERIALIZE_OPT_KEEP_SHARED);
  GObject *obj         = ncm_serialize_from_binfile (ser, full_filename);

  g_free (full_filename);
  ncm_serialize_free (ser);

  if (!NC_IS_GALAXY_WL_OBS (obj))
    g_error ("nc_galaxy_wl_obs_new_from_catalog_id: file for catalog id %d does not contain a NcGalaxyWLObs.", id);

  return NC_GALAXY_WL_OBS (obj);
}

/**
 * nc_galaxy_wl_obs_new_from_file: (constructor)
 * @filename: the name of a file with a serialized #NcGalaxyWLObs
 * @error: a #GError
 *
 * Reads the catalog of the binary file @filename (ncm_serialize_from_binfile()). A file
 * that holds another type of object sets %NC_GALAXY_WL_OBS_ERROR_WRONG_TYPE; a file that
 * cannot be read aborts with an error.
 *
 * Returns: (transfer full): a new #NcGalaxyWLObs object.
 */
NcGalaxyWLObs *
nc_galaxy_wl_obs_new_from_file (const gchar *filename, GError **error)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_KEEP_SHARED);
  GObject *obj      = ncm_serialize_from_binfile (ser, filename);

  ncm_serialize_free (ser);

  if (!NC_IS_GALAXY_WL_OBS (obj))
  {
    ncm_util_set_or_call_error (error, NC_GALAXY_WL_OBS_ERROR, NC_GALAXY_WL_OBS_ERROR_WRONG_TYPE,
                                "nc_galaxy_wl_obs_new_from_file: the file %s holds a %s, not a NcGalaxyWLObs.",
                                filename, G_OBJECT_TYPE_NAME (obj));
    g_clear_object (&obj);

    return NULL;
  }

  return NC_GALAXY_WL_OBS (obj);
}

/**
 * nc_galaxy_wl_obs_ref:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Increases the reference count of the #NcGalaxyWLObs object.
 *
 * Returns: (transfer full): the #NcGalaxyWLObs object.
 */
NcGalaxyWLObs *
nc_galaxy_wl_obs_ref (NcGalaxyWLObs *obs)
{
  return g_object_ref (obs);
}

/**
 * nc_galaxy_wl_obs_free:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Frees the #NcGalaxyWLObs object.
 *
 */
void
nc_galaxy_wl_obs_free (NcGalaxyWLObs *obs)
{
  g_object_unref (obs);
}

/**
 * nc_galaxy_wl_obs_clear:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Clears the #NcGalaxyWLObs object.
 *
 */
void
nc_galaxy_wl_obs_clear (NcGalaxyWLObs **obs)
{
  g_clear_object (obs);
}

/**
 * nc_galaxy_wl_obs_set_pz:
 * @obs: a #NcGalaxyWLObs object.
 * @i: row index.
 * @pz: a #NcGalaxyPz.
 *
 * Sets the redshift density of row @i to @pz, which is referenced. A #NcGalaxyPz does
 * not change after construction, so one object may be set in several rows and catalogs.
 *
 */
void
nc_galaxy_wl_obs_set_pz (NcGalaxyWLObs *obs, const guint i, NcGalaxyPz *pz)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  ncm_obj_dict_int_add (self->pz, i, G_OBJECT (pz));
  nc_galaxy_wl_obs_cheb_grid_clear (&self->cheb_grid);
}

/**
 * nc_galaxy_wl_obs_set_cheb_grid:
 * @obs: a #NcGalaxyWLObs object
 * @grid: (nullable): a #NcGalaxyWLObsChebGrid of the current rows
 *
 * Sets the CHEB_WEIGHTS redshift grid of the galaxies (#NcGalaxyWLObs:cheb-grid), or
 * removes it when @grid is %NULL. A grid of another number of galaxies or of rows of
 * another data hash aborts.
 */
void
nc_galaxy_wl_obs_set_cheb_grid (NcGalaxyWLObs *obs, NcGalaxyWLObsChebGrid *grid)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  nc_galaxy_wl_obs_cheb_grid_clear (&self->cheb_grid);
  self->cheb_grid_key = NCM_CATALOG_KEY_PENDING;

  if (grid == NULL)
    return;

  if (nc_galaxy_wl_obs_cheb_grid_len (grid) != ncm_catalog_len (NCM_CATALOG (obs)))
    g_error ("nc_galaxy_wl_obs_set_cheb_grid: the grid has %u galaxies, the catalog %u.",
             nc_galaxy_wl_obs_cheb_grid_len (grid), ncm_catalog_len (NCM_CATALOG (obs)));

  self->cheb_grid = nc_galaxy_wl_obs_cheb_grid_ref (grid);
  _nc_galaxy_wl_obs_cheb_grid_current (obs);
}

/**
 * nc_galaxy_wl_obs_peek_cheb_grid:
 * @obs: a #NcGalaxyWLObs object
 *
 * Gets the CHEB_WEIGHTS redshift grid of the galaxies; %NULL when there is none, after
 * the values of the rows changed or after a p(z) was set.
 *
 * Returns: (transfer none) (nullable): the #NcGalaxyWLObsChebGrid.
 */
NcGalaxyWLObsChebGrid *
nc_galaxy_wl_obs_peek_cheb_grid (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  return _nc_galaxy_wl_obs_cheb_grid_current (obs) ? self->cheb_grid : NULL;
}

/**
 * nc_galaxy_wl_obs_peek_pz:
 * @obs: a #NcGalaxyWLObs object.
 * @i: row index.
 *
 * Gets the redshift density of row @i.
 *
 * Returns: (transfer none) (nullable): the #NcGalaxyPz of row @i, or %NULL when the row
 * has none.
 *
 */
NcGalaxyPz *
nc_galaxy_wl_obs_peek_pz (NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  return (NcGalaxyPz *) ncm_obj_dict_int_peek (self->pz, i);
}

/* The interpolants of obs, created empty in memory when there are none */
static NcGalaxyWLObsShapeInterp *
_nc_galaxy_wl_obs_peek_or_new_shape_interp (NcGalaxyWLObsPrivate *self)
{
  if (self->shape_interp == NULL)
  {
    self->shape_interp     = nc_galaxy_wl_obs_shape_interp_new ();
    self->shape_interp_key = NCM_CATALOG_KEY_PENDING;
  }

  return self->shape_interp;
}

/**
 * nc_galaxy_wl_obs_set_shape_interp:
 * @obs: a #NcGalaxyWLObs
 * @i: row index
 * @st: the interpolant of row @i
 *
 * Sets the interpolant of $\ln P(\epsilon_\mathrm{obs} \mid g)$ of row @i, referencing
 * @st (nc_galaxy_wl_obs_shape_interp_set_interp()). It is returned by
 * nc_galaxy_wl_obs_peek_shape_interp() only after nc_galaxy_wl_obs_set_shape_interp_source()
 * and while the rows are unchanged.
 */
void
nc_galaxy_wl_obs_set_shape_interp (NcGalaxyWLObs *obs, const guint i, NcmSpectralTensor *st)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  nc_galaxy_wl_obs_shape_interp_set_interp (_nc_galaxy_wl_obs_peek_or_new_shape_interp (self), i, st, NULL);
}

/**
 * nc_galaxy_wl_obs_peek_shape_interp:
 * @obs: a #NcGalaxyWLObs
 * @i: row index
 *
 * Gets the interpolant of row @i (nc_galaxy_wl_obs_shape_interp_peek_interp(), which reads
 * the file of #NcGalaxyWLObs:shape-interp at the first call). A failure to read the file
 * aborts with the error of nc_galaxy_wl_obs_shape_interp_peek_interp().
 *
 * Returns: (transfer none) (nullable): the interpolant of row @i, or %NULL when the row has
 * none, there is no record or the rows changed after the record was set.
 */
NcmSpectralTensor *
nc_galaxy_wl_obs_peek_shape_interp (NcGalaxyWLObs *obs, const guint i)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  if (!_nc_galaxy_wl_obs_shape_interp_current (obs))
    return NULL;

  return nc_galaxy_wl_obs_shape_interp_peek_interp (self->shape_interp, i, NULL);
}

/**
 * nc_galaxy_wl_obs_set_shape_interp_source:
 * @obs: a #NcGalaxyWLObs
 * @source: the record of what built the interpolants
 *
 * Sets the record of the interpolants set with nc_galaxy_wl_obs_set_shape_interp(),
 * referencing @source, with the data hash of the rows (ncm_catalog_get_data_hash()), and
 * records ncm_catalog_get_key(): the interpolants describe the rows as they are now. Set it
 * after the last change of the rows.
 */
void
nc_galaxy_wl_obs_set_shape_interp_source (NcGalaxyWLObs *obs, NcmObjDictStr *source)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);
  NcGalaxyWLObsShapeInterp *si      = _nc_galaxy_wl_obs_peek_or_new_shape_interp (self);
  gchar *data_hash                  = ncm_catalog_get_data_hash (NCM_CATALOG (obs));

  nc_galaxy_wl_obs_shape_interp_set_source (si, source, NULL);
  nc_galaxy_wl_obs_shape_interp_set_data_hash (si, data_hash, NULL);
  self->shape_interp_key = ncm_catalog_get_key (NCM_CATALOG (obs));

  g_free (data_hash);
}

/**
 * nc_galaxy_wl_obs_peek_shape_interp_source:
 * @obs: a #NcGalaxyWLObs
 *
 * Gets the record of the interpolants, without reading the file of
 * #NcGalaxyWLObs:shape-interp.
 *
 * Returns: (transfer none) (nullable): the record of the interpolants, or %NULL when there is
 * none or the rows changed after it was set.
 */
NcmObjDictStr *
nc_galaxy_wl_obs_peek_shape_interp_source (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  if (!_nc_galaxy_wl_obs_shape_interp_current (obs))
    return NULL;

  return nc_galaxy_wl_obs_shape_interp_peek_source (self->shape_interp);
}

/**
 * nc_galaxy_wl_obs_set_coord:
 * @obs: a #NcGalaxyWLObs object.
 * @coord: the ellipticity handedness frame #NcWLEllipticityFrame.
 *
 * Sets the frame in which the stored ellipticity components are expressed. Sky
 * positions are always RA/Dec; this only selects the ellipticity basis (see
 * #NcWLEllipticityFrame).
 *
 */
void
nc_galaxy_wl_obs_set_coord (NcGalaxyWLObs *obs, NcWLEllipticityFrame coord)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  self->coord = coord;
}

static void
_nc_galaxy_wl_obs_set_ellip_conv (NcGalaxyWLObs *obs, NcGalaxyWLObsEllipConv type)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  self->ellip_conv = type;
}

/**
 * nc_galaxy_wl_obs_get_ellip_conv:
 * @obs: a #NcGalaxyWLObs object
 *
 * Gets the ellipse type used to store the data.
 *
 * Returns: the ellipse type #NcGalaxyWLObsEllipConv
 */
NcGalaxyWLObsEllipConv
nc_galaxy_wl_obs_get_ellip_conv (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  return self->ellip_conv;
}

/**
 * nc_galaxy_wl_obs_get_coord:
 * @obs: a #NcGalaxyWLObs object.
 *
 * Gets the frame in which the stored ellipticity components are expressed (see
 * #NcWLEllipticityFrame).
 *
 * Returns: the ellipticity handedness frame #NcWLEllipticityFrame.
 *
 */
NcWLEllipticityFrame
nc_galaxy_wl_obs_get_coord (NcGalaxyWLObs *obs)
{
  NcGalaxyWLObsPrivate * const self = nc_galaxy_wl_obs_get_instance_private (obs);

  return self->coord;
}

