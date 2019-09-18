/***************************************************************************
 *            nc_data_reduced_shear_cluster_mass.c
 *
 *  Thu Mar 22 16:10:25 2018
 *  Copyright  2018  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_data_reduced_shear_cluster_mass.c
 * Copyright (C) 2018 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_data_reduced_shear_cluster_mass
 * @title: NcDataReducedShearClusterMass
 * @short_description: Galaxy clusters data -- pseudo number counts likelihood.
 * 
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_reduced_shear_cluster_mass.h"
#include "lss/nc_reduced_shear_calib.h"
#include "nc_hicosmo.h"
#include "lss/nc_cluster_redshift.h"
#include "lss/nc_cluster_mass.h"
#include "lss/nc_cluster_abundance.h"
#include "lss/nc_galaxy_redshift_spec.h"
#include "lss/nc_galaxy_redshift_spline.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/integral.h"
#include "math/ncm_memory_pool.h"

#ifndef NUMCOSMO_GIR_SCAN
#ifdef HAVE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif /* HAVE_HDF5  */
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcDataReducedShearClusterMassPrivate
{
	NcDistance *dist;
	NcmObjArray *photoz_array;
	NcmMatrix *gal_obs;
	NcmVector *rh_vec;
	gboolean has_rh;
	gdouble psf_size;
	gdouble z_cluster;
	gdouble ra_cluster;
	gdouble dec_cluster;
};

enum
{
  PROP_0,
	PROP_DIST,
  PROP_PHOTOZ_ARRAY,
  PROP_GAL_OBS,
	PROP_HAS_RH,
	PROP_PSF_SIZE,
  PROP_Z_CLUSTER, 
  PROP_RA_CLUSTER, 
  PROP_DEC_CLUSTER, 
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcDataReducedShearClusterMass, nc_data_reduced_shear_cluster_mass, NCM_TYPE_DATA);

static void
nc_data_reduced_shear_cluster_mass_init (NcDataReducedShearClusterMass *drs)
{
	NcDataReducedShearClusterMassPrivate * const self = drs->priv = nc_data_reduced_shear_cluster_mass_get_instance_private (drs);
	self->dist         = NULL;
  self->photoz_array = ncm_obj_array_new ();
  self->gal_obs      = NULL;
  self->z_cluster    = 0.0;
  self->ra_cluster   = 0.0;
  self->dec_cluster  = 0.0;
	self->has_rh       = FALSE;
	self->psf_size     = 0.0;
}

static void
nc_data_reduced_shear_cluster_mass_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (object);
	NcDataReducedShearClusterMassPrivate * const self = drs->priv;
  g_return_if_fail (NC_IS_DATA_REDUCED_SHEAR_CLUSTER_MASS (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_data_reduced_shear_cluster_mass_set_dist (drs, g_value_get_object (value));
      break;
    case PROP_PHOTOZ_ARRAY:
		{
			NcmObjArray *photoz_array = g_value_get_boxed (value);
			
			g_clear_pointer (&self->photoz_array, ncm_obj_array_unref);
			self->photoz_array = ncm_obj_array_ref (photoz_array);

			if (self->gal_obs != NULL)
			{
				g_assert_cmpuint (self->photoz_array->len, ==, ncm_matrix_nrows (self->gal_obs));
			}

      break;
		}
    case PROP_GAL_OBS:
		{
			NcmMatrix *gal_obs = g_value_get_object (value);

			ncm_matrix_clear (&self->gal_obs);
			self->gal_obs = ncm_matrix_ref (gal_obs);

			if (self->photoz_array != NULL)
			{
				g_assert_cmpuint (self->photoz_array->len, ==, ncm_matrix_nrows (self->gal_obs));
			}

			break;
		}
		case PROP_HAS_RH:
			self->has_rh = g_value_get_boolean (value);
      break;
		case PROP_PSF_SIZE:
			self->psf_size = g_value_get_double (value);
      break;
    case PROP_Z_CLUSTER:
      self->z_cluster = g_value_get_double (value);
      break;
    case PROP_RA_CLUSTER:
      self->ra_cluster = g_value_get_double (value);
      break;
    case PROP_DEC_CLUSTER:
      self->dec_cluster = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_reduced_shear_cluster_mass_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (object);
	NcDataReducedShearClusterMassPrivate * const self = drs->priv;
  g_return_if_fail (NC_IS_DATA_REDUCED_SHEAR_CLUSTER_MASS (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, self->dist);
      break;
		case PROP_PHOTOZ_ARRAY:
      g_value_set_boxed (value, self->photoz_array);
      break;
		case PROP_GAL_OBS:
      g_value_set_object (value, self->gal_obs);
      break;
		case PROP_HAS_RH:
			g_value_set_boolean (value, self->has_rh);
      break;
		case PROP_PSF_SIZE:
			g_value_set_double (value, self->psf_size);
      break;
    case PROP_Z_CLUSTER:
      g_value_set_double (value, self->z_cluster);
      break;
    case PROP_RA_CLUSTER:
      g_value_set_double (value, self->ra_cluster);
      break;
    case PROP_DEC_CLUSTER:
      g_value_set_double (value, self->dec_cluster);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_reduced_shear_cluster_mass_dispose (GObject *object)
{
  NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (object);
	NcDataReducedShearClusterMassPrivate * const self = drs->priv;

  ncm_matrix_clear (&self->gal_obs);
	ncm_obj_array_clear (&self->photoz_array);
	nc_distance_clear (&self->dist);
	
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_reduced_shear_cluster_mass_parent_class)->dispose (object);
}

static void
nc_data_reduced_shear_cluster_mass_finalize (GObject *object)
{
  /*NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (object);*/

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_reduced_shear_cluster_mass_parent_class)->finalize (object);
}

static void _nc_data_reduced_shear_cluster_mass_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);
static guint _nc_data_reduced_shear_cluster_mass_get_len (NcmData *data);
static void _nc_data_reduced_shear_cluster_mass_prepare (NcmData *data, NcmMSet *mset);

static void
nc_data_reduced_shear_cluster_mass_class_init (NcDataReducedShearClusterMassClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class = NCM_DATA_CLASS (klass);
	
  object_class->set_property = nc_data_reduced_shear_cluster_mass_set_property;
  object_class->get_property = nc_data_reduced_shear_cluster_mass_get_property;
  object_class->dispose      = nc_data_reduced_shear_cluster_mass_dispose;
  object_class->finalize     = nc_data_reduced_shear_cluster_mass_finalize;

	g_object_class_install_property (object_class,
	                                 PROP_DIST,
	                                 g_param_spec_object ("dist",
	                                                      NULL,
	                                                      "Distance object",
	                                                      NC_TYPE_DISTANCE,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_PHOTOZ_ARRAY,
	                                 g_param_spec_boxed ("photoz-array",
	                                                     NULL,
	                                                     "Array of photometric redshift objects",
	                                                     NCM_TYPE_OBJ_ARRAY,
	                                                     G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_GAL_OBS,
	                                 g_param_spec_object ("gal-obs",
	                                                      NULL,
	                                                      "Matrix containing galaxy observables",
	                                                      NCM_TYPE_MATRIX,
	                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_HAS_RH,
                                   g_param_spec_boolean ("has-rh",
                                                         NULL,
                                                         "Has the galaxy size (rh) information",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
                                   PROP_PSF_SIZE,
                                   g_param_spec_double ("psf-size",
                                                         NULL,
                                                         "PSF size",
                                                         0.0, G_MAXDOUBLE, 0.0,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	
	g_object_class_install_property (object_class,
                                   PROP_Z_CLUSTER,
                                   g_param_spec_double ("z-cluster",
                                                         NULL,
                                                         "Cluster (halo) redshift",
                                                         0.0, G_MAXDOUBLE, 0.0,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
                                   PROP_RA_CLUSTER,
                                   g_param_spec_double ("ra-cluster",
                                                         NULL,
                                                         "Cluster (halo) RA",
                                                         -360.0, 360.0, 0.0,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
                                   PROP_DEC_CLUSTER,
                                   g_param_spec_double ("dec-cluster",
                                                         NULL,
                                                         "Cluster (halo) DEC",
                                                         -360.0, 360.0, 0.0,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
 
  data_class->m2lnL_val  = &_nc_data_reduced_shear_cluster_mass_m2lnL_val;
  data_class->get_length = &_nc_data_reduced_shear_cluster_mass_get_len;
  data_class->prepare    = &_nc_data_reduced_shear_cluster_mass_prepare;
}

typedef struct _NcDataReducedShearClusterMassInteg
{
	NcGalaxyRedshift *gz;
	NcReducedShearClusterMass *rs;
	NcReducedShearCalib *rs_calib;
	NcDistance *dist;
	NcHICosmo *cosmo;
	NcWLSurfaceMassDensity *smd;
	NcDensityProfile *dp;
	gdouble z_cluster;
	gdouble R_Mpc;
	gdouble dt_cluster;
	gdouble g_obs;
	guint interval_index;
	gdouble psf_size;
	gdouble rh;
	} NcDataReducedShearClusterMassInteg;

static gdouble
_nc_data_reduced_shear_cluster_mass_Pgal (gdouble z_gal, NcDataReducedShearClusterMassInteg *integ)
{
	const gdouble g_th = nc_wl_surface_mass_density_reduced_shear_infinity (integ->smd, integ->dp, integ->cosmo, integ->R_Mpc, z_gal, integ->z_cluster, integ->z_cluster);
	const gdouble Pgal = nc_reduced_shear_cluster_mass_P_z_gth_gobs (integ->rs, integ->cosmo, z_gal, g_th, integ->g_obs);
	
	return Pgal;
}

static gdouble
_nc_data_reduced_shear_cluster_mass_Pgal_calib (gdouble z_gal, NcDataReducedShearClusterMassInteg *integ)
{
	const gdouble g_th       = nc_wl_surface_mass_density_reduced_shear_infinity (integ->smd, integ->dp, integ->cosmo, integ->R_Mpc, z_gal, integ->z_cluster, integ->z_cluster);
	const gdouble g_th_calib = nc_reduced_shear_calib_eval (integ->rs_calib, g_th, integ->psf_size, integ->rh);
	const gdouble Pgal       = nc_reduced_shear_cluster_mass_P_z_gth_gobs (integ->rs, integ->cosmo, z_gal, g_th_calib, integ->g_obs);
	return Pgal;
}

static gdouble
_nc_data_reduced_shear_cluster_mass_PgalPz (gdouble z_gal, NcDataReducedShearClusterMassInteg *integ)
{
	const gdouble Pgal = _nc_data_reduced_shear_cluster_mass_Pgal (z_gal, integ);
	const gdouble Pz   = nc_galaxy_redshift_pdf (integ->gz, integ->interval_index, z_gal);

	return Pgal * Pz;
}

static gdouble
_nc_data_reduced_shear_cluster_mass_PgalPz_calib (gdouble z_gal, NcDataReducedShearClusterMassInteg *integ)
{
	const gdouble Pgal = _nc_data_reduced_shear_cluster_mass_Pgal_calib (z_gal, integ);
	const gdouble Pz   = nc_galaxy_redshift_pdf (integ->gz, integ->interval_index, z_gal);

	return Pgal * Pz;
}

static void
_nc_data_reduced_shear_cluster_mass_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (data);
	NcDataReducedShearClusterMassPrivate * const self = drs->priv;
  NcHICosmo *cosmo              = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcWLSurfaceMassDensity *smd   = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcDensityProfile *dp          = NC_DENSITY_PROFILE (ncm_mset_peek (mset, nc_density_profile_id ()));
	NcReducedShearClusterMass *rs = NC_REDUCED_SHEAR_CLUSTER_MASS (ncm_mset_peek (mset, nc_reduced_shear_cluster_mass_id ()));
	NcReducedShearCalib *rs_calib = NC_REDUCED_SHEAR_CALIB (ncm_mset_peek (mset, nc_reduced_shear_calib_id ()));
	const guint ngal              = self->photoz_array->len;
	const gdouble RH              = nc_hicosmo_RH_Mpc (cosmo);
	const gdouble dA              = nc_distance_angular_diameter (self->dist, cosmo, self->z_cluster) * RH;
	const gdouble dt              = nc_distance_transverse (self->dist, cosmo, self->z_cluster);
	gsl_integration_workspace **w = ncm_integral_get_workspace ();
  	
	NcDataReducedShearClusterMassInteg integ_data;
	gint i;

	gdouble (*Pgal) (gdouble z_gal, NcDataReducedShearClusterMassInteg *integ);
  gdouble (*PgalPz) (gdouble z_gal, NcDataReducedShearClusterMassInteg *integ);
	
  g_assert (cosmo != NULL);
  g_assert (smd != NULL);
  g_assert (dp != NULL);
  g_assert (rs != NULL);

	integ_data.rs         = rs;
	integ_data.rs_calib   = rs_calib;
	integ_data.dist       = self->dist;
	integ_data.cosmo      = cosmo;
	integ_data.smd        = smd;
	integ_data.dp         = dp;
	integ_data.z_cluster  = self->z_cluster;
	integ_data.dt_cluster = dt;
  integ_data.psf_size   = self->psf_size;
	
	m2lnL[0] = 0.0;

	if (self->has_rh == TRUE)
	{
		Pgal   = &_nc_data_reduced_shear_cluster_mass_Pgal_calib;
		PgalPz = &_nc_data_reduced_shear_cluster_mass_PgalPz_calib;		
		g_assert (rs_calib != NULL);
	}
	else
	{
		Pgal   = &_nc_data_reduced_shear_cluster_mass_Pgal;
		PgalPz = &_nc_data_reduced_shear_cluster_mass_PgalPz;
	}
		
	for (i = 0; i < ngal; i++)
	{
		NcGalaxyRedshift *gz      = NC_GALAXY_REDSHIFT (ncm_obj_array_peek (self->photoz_array, i));
		const gdouble r_arcmin    = ncm_matrix_get (self->gal_obs, i, 0);
		const gdouble g_obs       = ncm_matrix_get (self->gal_obs, i, 1);
		const gdouble rh          = self->has_rh ? ncm_matrix_get (self->gal_obs, i, 2) : 0.0;
		const gdouble R_Mpc       = (r_arcmin / 60.0) * (M_PI / 180.0) * dA;

		const gdouble z_gal       = nc_galaxy_redshift_mode (gz);

		if (R_Mpc < 0.75 || R_Mpc > 3.0 || z_gal < self->z_cluster + 0.1 || z_gal > 1.25)
		{
			//printf ("r_arcmin = %.5g R_Mpc = %.5g\n", r_arcmin, R_Mpc);
			//printf ("z_gal = %.5g\n", z_gal);
			continue;
		}

    integ_data.gz    = gz;
		integ_data.R_Mpc = R_Mpc;
		integ_data.g_obs = g_obs;
		integ_data.rh    = rh;
		
		if (!nc_galaxy_redshift_has_dist (gz))
		{
			//const gdouble z_gal = nc_galaxy_redshift_mode (gz);
			const gdouble P_i   = Pgal (z_gal, &integ_data);
			//printf ("zgal = %.5g\n", z_gal);
			m2lnL[0] += log (P_i);
		}
		else
		{
			const guint ndists = nc_galaxy_redshift_nintervals (gz);
			gint j;

			for (j = 0; j < ndists; j++)
			{
				gdouble z_gal_lower, z_gal_upper, P_ij, abserr;
				gsl_function F;
				
				nc_galaxy_redshift_pdf_limits (gz, j, &z_gal_lower, &z_gal_upper);

				//printf ("zgal_l = %.5g zgal_u = %.5g\n", z_gal_lower, z_gal_upper);
				integ_data.interval_index = j;

				F.params   = &integ_data;
				F.function = (gdouble (*)(gdouble,gpointer)) PgalPz;
				
				gsl_integration_qag (&F, z_gal_lower, z_gal_upper, 0.0, 1.0e-5, NCM_INTEGRAL_PARTITION, 1, *w, &P_ij, &abserr);

				m2lnL[0] += log (P_ij);
			}
		}
	}

	m2lnL[0] = -2.0 * m2lnL[0];
	
	ncm_memory_pool_return (w);
  return;
}

static guint 
_nc_data_reduced_shear_cluster_mass_get_len (NcmData *data)
{ 
	NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (data);
	NcDataReducedShearClusterMassPrivate * const self = drs->priv;

	return (self->photoz_array != NULL) ? self->photoz_array->len : 0; 
}

static void 
_nc_data_reduced_shear_cluster_mass_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (data);
	NcDataReducedShearClusterMassPrivate * const self = drs->priv;
	NcHICosmo *cosmo              = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcWLSurfaceMassDensity *smd   = NC_WL_SURFACE_MASS_DENSITY (ncm_mset_peek (mset, nc_wl_surface_mass_density_id ()));
  NcDensityProfile *dp          = NC_DENSITY_PROFILE (ncm_mset_peek (mset, nc_density_profile_id ()));
  NcReducedShearClusterMass *rs = NC_REDUCED_SHEAR_CLUSTER_MASS (ncm_mset_peek (mset, nc_reduced_shear_cluster_mass_id ()));
	
  g_assert ((cosmo != NULL) && (smd != NULL) && (dp != NULL) && (rs != NULL));

  nc_distance_prepare_if_needed (self->dist, cosmo);
}


/**
 * nc_data_reduced_shear_cluster_mass_new:
 * @dist: a #NcDistance 
 * 
 * Creates a new #NcDataReducedShearClusterMass.
 * 
 * Returns: the newly created #NcDataReducedShearClusterMass.
 */
NcDataReducedShearClusterMass *
nc_data_reduced_shear_cluster_mass_new (NcDistance *dist)
{
	NcDataReducedShearClusterMass *drs = g_object_new (NC_TYPE_DATA_REDUCED_SHEAR_CLUSTER_MASS,
	                                                   "dist", dist,
	                                                   NULL);
	
  return drs;
}

/**
 * nc_data_reduced_shear_cluster_mass_new_from_file:
 * @filename: file containing a serialized #NcDataReducedShearClusterMass
 * 
 * Creates a new #NcDataReducedShearClusterMass from @filename.
 * 
 * Returns: (transfer full): the newly created #NcDataReducedShearClusterMass.
 */
NcDataReducedShearClusterMass *
nc_data_reduced_shear_cluster_mass_new_from_file (const gchar *filename)
{
  NcDataReducedShearClusterMass *drs = NC_DATA_REDUCED_SHEAR_CLUSTER_MASS (ncm_serialize_global_from_file (filename));
  g_assert (NC_IS_DATA_REDUCED_SHEAR_CLUSTER_MASS (drs));

  return drs;
}

/**
 * nc_data_reduced_shear_cluster_mass_ref:
 * @drs: a #NcDataReducedShearClusterMass
 *
 * Increases the reference count of @drs by one.
 * 
 * Returns: (transfer full): @drs
 */
NcDataReducedShearClusterMass *
nc_data_reduced_shear_cluster_mass_ref (NcDataReducedShearClusterMass *drs)
{
  return g_object_ref (drs);
}

/**
 * nc_data_reduced_shear_cluster_mass_free:
 * @drs: a #NcDataReducedShearClusterMass
 *
 * Atomically decrements the reference count of @drs by one. If the reference count drops to 0,
 * all memory allocated by @drs is released.
 * 
 */
void
nc_data_reduced_shear_cluster_mass_free (NcDataReducedShearClusterMass *drs)
{
  g_object_unref (drs);
}

/**
 * nc_data_reduced_shear_cluster_mass_clear:
 * @drs: a #NcDataReducedShearClusterMass
 *
 * The reference count of @drs is decreased and the pointer is set to NULL.
 * 
 */
void
nc_data_reduced_shear_cluster_mass_clear (NcDataReducedShearClusterMass **drs)
{
  g_clear_object (drs);
}

/**
 * nc_data_reduced_shear_cluster_mass_set_dist:
 * @drs: a NcDataReducedShearClusterMass
 * @dist: a #NcDistance
 * 
 * Sets the distance object.
 * 
 */
void 
nc_data_reduced_shear_cluster_mass_set_dist (NcDataReducedShearClusterMass *drs, NcDistance *dist)
{
	NcDataReducedShearClusterMassPrivate * const self = drs->priv;
	
  nc_distance_clear (&self->dist);
  self->dist = nc_distance_ref (dist);
}

#ifdef HAVE_HDF5
typedef struct _NcmHDF5Table
{
	gchar *name;
	hsize_t nfields;
	hsize_t nrecords;
	GPtrArray *field_names;
	GArray *field_sizes;
	GArray *field_offsets;
	size_t type_size;
	hid_t h5f;
	hid_t table_h5d;
	hid_t table_h5t;
} NcmHDF5Table;

NcmHDF5Table *
ncm_hdf5_table_new (hid_t h5f, const gchar *tname)
{
	NcmHDF5Table *h5tb = g_new (NcmHDF5Table, 1);
	gint ret, i;

	/* Table info: init */
	h5tb->h5f           = h5f;
	h5tb->name          = g_strdup (tname);
	h5tb->field_names   = g_ptr_array_new_with_free_func (g_free);
	h5tb->field_sizes   = g_array_new (TRUE, TRUE, sizeof (size_t));
	h5tb->field_offsets = g_array_new (TRUE, TRUE, sizeof (size_t));
		
	/* Table info: read table info */
	ret = H5TBget_table_info (h5f, tname, &h5tb->nfields, &h5tb->nrecords);
	g_assert_cmpint (ret, !=, -1);

	g_ptr_array_set_size (h5tb->field_names,   h5tb->nfields);
	g_array_set_size     (h5tb->field_sizes,   h5tb->nfields);
	g_array_set_size     (h5tb->field_offsets, h5tb->nfields);
	
	for (i = 0; i < h5tb->nfields; i++) 
		g_ptr_array_index (h5tb->field_names, i) = g_new (gchar, 255);

	/* Table info: read fields info */
	ret = H5TBget_field_info (h5f, tname, (gchar **)h5tb->field_names->pdata, (size_t *)h5tb->field_sizes->data, (size_t *)h5tb->field_offsets->data, &h5tb->type_size);
	g_assert_cmpint (ret, !=, -1);

	/* Table datatype */
	h5tb->table_h5d = H5Dopen2 (h5f, tname, H5P_DEFAULT);
	h5tb->table_h5t = H5Dget_type (h5tb->table_h5d);

	/* Consistence check */
	{
		H5T_class_t table_h5t_class = H5Tget_class (h5tb->table_h5t);
		if (table_h5t_class != H5T_COMPOUND)
			g_error ("ncm_hdf5_table_new: only support tables with compound elements.");
	}

	return h5tb;
}

void
ncm_hdf5_table_free (NcmHDF5Table *h5tb)
{
	g_free (h5tb->name);
	g_ptr_array_unref (h5tb->field_names);
	g_array_unref (h5tb->field_sizes);
	g_array_unref (h5tb->field_offsets);

	H5Dclose (h5tb->table_h5d);
	H5Tclose (h5tb->table_h5t);

	g_free (h5tb);
}

#if !GLIB_CHECK_VERSION(2,54,0)
static gboolean
g_ptr_array_find_with_equal_func (GPtrArray     *haystack,
                                  gconstpointer  needle,
                                  GEqualFunc     equal_func,
                                  guint         *index_)
{
  guint i;

  g_return_val_if_fail (haystack != NULL, FALSE);

  if (equal_func == NULL)
    equal_func = g_direct_equal;

  for (i = 0; i < haystack->len; i++)
    {
      if (equal_func (g_ptr_array_index (haystack, i), needle))
        {
          if (index_ != NULL)
            *index_ = i;
          return TRUE;
        }
    }

  return FALSE;
}
#endif

gboolean
ncm_hdf5_table_has_col (NcmHDF5Table *h5tb, const gchar *col)
{
	return g_ptr_array_find_with_equal_func (h5tb->field_names, col, g_str_equal, NULL);
}

void
ncm_hdf5_table_read_col_as_vec (NcmHDF5Table *h5tb, const gchar *col, NcmVector **vec)
{
	size_t field_sizes_sd[1] = { sizeof (gdouble) };
	size_t field_offset[1]   = { 0 };
	hid_t mtype, mntype;
	guint _cindex = 0;
	gint cindex;
	gint mi, ret;
	
	if (*vec != NULL)
		g_assert_cmpuint (ncm_vector_len (*vec), ==, h5tb->nrecords);
	else
		*vec = ncm_vector_new (h5tb->nrecords);

	mi      = H5Tget_member_index (h5tb->table_h5t, col);
	mtype   = H5Tget_member_type (h5tb->table_h5t, mi);
	mntype  = H5Tget_native_type (mtype, H5T_DIR_ASCEND);

	if (!g_ptr_array_find_with_equal_func (h5tb->field_names, col, g_str_equal, &_cindex))
		g_error ("ncm_hdf5_table_read_col_as_vec: column `%s' not found in table `%s'.", col, h5tb->name);

	cindex = _cindex;

	g_assert_cmpint (cindex, ==, mi);
	
	if (H5Tequal (mntype, H5T_NATIVE_FLOAT) > 0)
	{
		GArray *tmp = g_array_new (FALSE, FALSE, sizeof (gfloat)); 
		gint i;
		
		g_array_set_size (tmp, h5tb->nrecords);

		field_sizes_sd[0] = sizeof (gfloat);
		ret = H5TBread_fields_index (h5tb->h5f, h5tb->name, 1, &cindex, 0, h5tb->nrecords, sizeof (gfloat), field_offset, field_sizes_sd, tmp->data);
		g_assert_cmpint (ret, !=, -1);

		for (i = 0; i < h5tb->nrecords; i++)
		{
			ncm_vector_set (*vec, i, g_array_index (tmp, gfloat, i));
		}
		
		g_array_unref (tmp);
	}
	else
	{
		g_assert_cmpint (H5Tequal (mntype, H5T_NATIVE_DOUBLE), >, 0);

		ret = H5TBread_fields_index (h5tb->h5f, h5tb->name, 1, &cindex, 0, h5tb->nrecords, sizeof (gdouble), field_offset, field_sizes_sd, ncm_vector_data (*vec));
		g_assert_cmpint (ret, !=, -1);
	}


	ret = H5Tclose (mntype);
	g_assert_cmpint (ret, !=, -1);
	ret = H5Tclose (mtype);
	g_assert_cmpint (ret, !=, -1);
}

void
ncm_hdf5_table_read_col_as_longarray (NcmHDF5Table *h5tb, const gchar *col, GArray **a)
{
	size_t field_sizes_sd[1] = { sizeof (glong) };
	size_t field_offset[1]   = { 0 };
	hid_t mtype, mntype;
	guint _cindex = 0;
	gint cindex;
	gint mi, ret;
	
	if (*a != NULL)
	{
		g_array_set_size (*a, h5tb->nrecords);
	}
	else
	{
		*a = g_array_new (TRUE, TRUE, sizeof (glong));
		g_array_set_size (*a, h5tb->nrecords);
	}

	mi      = H5Tget_member_index (h5tb->table_h5t, col);
	mtype   = H5Tget_member_type (h5tb->table_h5t, mi);
	mntype  = H5Tget_native_type (mtype, H5T_DIR_ASCEND);

  if (!g_ptr_array_find_with_equal_func (h5tb->field_names, col, g_str_equal, &_cindex))
		g_error ("ncm_hdf5_table_read_col_as_longarray: column `%s' not found in table `%s'.", col, h5tb->name);
	cindex = _cindex;

	g_assert_cmpint (cindex, ==, mi);


	if (H5Tequal (mntype, H5T_NATIVE_INT) > 0)
	{
		GArray *tmp = g_array_new (TRUE, TRUE, sizeof (gint)); 
		gint i;

		g_array_set_size (tmp, h5tb->nrecords);

		//printf ("name '%s', nrecords %llu\n", h5tb->name, h5tb->nrecords);
		field_sizes_sd[0] = sizeof (gint);
		ret = H5TBread_fields_index (h5tb->h5f, h5tb->name, 1, &cindex, 0, h5tb->nrecords, sizeof (gint), field_offset, field_sizes_sd, tmp->data);
		g_assert_cmpint (ret, !=, -1);

		for (i = 0; i < h5tb->nrecords; i++)
		{
			g_array_index (*a, glong, i) = g_array_index (tmp, gint, i);
		}

		g_array_unref (tmp);
	}
  else
	{
		g_assert_cmpint (H5Tequal (mntype, H5T_NATIVE_LONG), >, 0);
		
		ret = H5TBread_fields_index (h5tb->h5f, h5tb->name, 1, &cindex, 0, h5tb->nrecords, sizeof (glong), field_offset, field_sizes_sd, (*a)->data);
	  g_assert_cmpint (ret, !=, -1);
	}

	ret = H5Tclose (mntype);
	g_assert_cmpint (ret, !=, -1);
	ret = H5Tclose (mtype);
	g_assert_cmpint (ret, !=, -1);
}

void
ncm_hdf5_table_read_col_as_mat (NcmHDF5Table *h5tb, const gchar *col, NcmMatrix **mat)
{
	size_t field_offset[1] = { 0 };
	size_t field_sizes_sd[1];
	hid_t mtype, btype, mntype;
	guint _cindex = 0;
	gint cindex;
	gint mi, ret;
	H5T_class_t mtype_class;
	hsize_t ncols;
	
	mi          = H5Tget_member_index (h5tb->table_h5t, col);
	mtype       = H5Tget_member_type (h5tb->table_h5t, mi);
	mtype_class = H5Tget_class (mtype);

	g_assert_cmpint (mtype_class, ==, H5T_ARRAY);

	btype  = H5Tget_super (mtype);
	mntype = H5Tget_native_type (btype, H5T_DIR_ASCEND);

	g_assert_cmpint (H5Tget_array_ndims (mtype), ==, 1);

	ret = H5Tget_array_dims2 (mtype, &ncols); 
	g_assert_cmpint (ret, !=, -1);
	
	if (*mat != NULL)
	{
		g_assert_cmpuint (ncm_matrix_nrows (*mat), ==, h5tb->nrecords);
		g_assert_cmpuint (ncm_matrix_ncols (*mat), ==, ncols);
	}
	else
		*mat = ncm_matrix_new (h5tb->nrecords, ncols);

	if (!g_ptr_array_find_with_equal_func (h5tb->field_names, col, g_str_equal, &_cindex))
		g_error ("ncm_hdf5_table_read_col_as_mat: column `%s' not found in table `%s'.", col, h5tb->name);
	cindex = _cindex;

	g_assert_cmpint (cindex, ==, mi);

	if (H5Tequal (mntype, H5T_NATIVE_FLOAT) > 0)
	{
		GArray *tmp = g_array_new (FALSE, FALSE, sizeof (gfloat)); 
		gint i;
		
		g_array_set_size (tmp, h5tb->nrecords * ncols);

		field_sizes_sd[0] = sizeof (gfloat) * ncols;
		ret = H5TBread_fields_index (h5tb->h5f, h5tb->name, 1, &cindex, 0, h5tb->nrecords, sizeof (gfloat) * ncols, field_offset, field_sizes_sd, tmp->data);
		g_assert_cmpint (ret, !=, -1);

		for (i = 0; i < h5tb->nrecords; i++)
		{
			gint j;
			gfloat *row = &g_array_index (tmp, gfloat, i * ncols);
			for (j = 0; j < ncols; j++)
			{
				ncm_matrix_set (*mat, i, j, row[j]);
			}
		}
		
		g_array_unref (tmp);
	}
	else
	{
		g_assert_cmpint (H5Tequal (mntype, H5T_NATIVE_DOUBLE), >, 0);

	  field_sizes_sd[0] = sizeof (gdouble) * ncols;
	  ret = H5TBread_fields_index (h5tb->h5f, h5tb->name, 1, &cindex, 0, h5tb->nrecords, sizeof (gdouble) * ncols, field_offset, field_sizes_sd, ncm_matrix_data (*mat));
	  g_assert_cmpint (ret, !=, -1);
	}

	ret = H5Tclose (mntype);
	g_assert_cmpint (ret, !=, -1);

	ret = H5Tclose (btype);
	g_assert_cmpint (ret, !=, -1);

	ret = H5Tclose (mtype);
	g_assert_cmpint (ret, !=, -1);
}

void
ncm_hdf5_table_read_col_as_fixstr (NcmHDF5Table *h5tb, const gchar *col, GArray **a)
{
	size_t field_offset[1] = { 0 };
	size_t field_sizes_sd[1];
	hid_t mtype;
	guint _cindex = 0;
	gint cindex;
	gint mi, ret;
	H5T_class_t mtype_class;
	size_t len;
	
	mi          = H5Tget_member_index (h5tb->table_h5t, col);
	mtype       = H5Tget_member_type (h5tb->table_h5t, mi);
	mtype_class = H5Tget_class (mtype);

	g_assert_cmpint (mtype_class, ==, H5T_STRING);
	g_assert_cmpint (H5Tis_variable_str (mtype), ==, 0);
	g_assert_cmpint (H5Tget_cset (mtype), ==, H5T_CSET_ASCII);

	len = H5Tget_size (mtype);
	
	if (*a != NULL)
		g_assert_cmpuint (g_array_get_element_size (*a), ==, len);
	else
		*a = g_array_new (TRUE, TRUE, len);
	g_array_set_size (*a, h5tb->nrecords);
	
	if (!g_ptr_array_find_with_equal_func (h5tb->field_names, col, g_str_equal, &_cindex))
		g_error ("ncm_hdf5_table_read_col_as_fixstr: column `%s' not found in table `%s'.", col, h5tb->name);
	cindex = _cindex;

	g_assert_cmpint (cindex, ==, mi);
	
	field_sizes_sd[0] = sizeof (gchar) * len;
	ret = H5TBread_fields_index (h5tb->h5f, h5tb->name, 1, &cindex, 0, h5tb->nrecords, field_sizes_sd[0], field_offset, field_sizes_sd, (*a)->data);
	g_assert_cmpint (ret, !=, -1);

	ret = H5Tclose (mtype);
	g_assert_cmpint (ret, !=, -1);
}

#endif /* HAVE_HDF5 */

static gboolean
_g_long_equal (gconstpointer v1, gconstpointer v2)
{
  return *((const glong*) v1) == *((const glong*) v2);
}

static guint
_g_long_hash (gconstpointer v)
{
  return (guint) *(const glong*) v;
}

/**
 * nc_data_reduced_shear_cluster_mass_load_hdf5:
 * @drs: a #NcDataReducedShearClusterMass
 * @hdf5_file: a file containing #NcDataReducedShearClusterMass data
 * @ftype: filter type ('u', 'g', 'r', 'i', 'z')
 * @z_cluster: cluster redshift
 * @ra_cluster: cluster RA
 * @dec_cluster: cluster DEC
 * 
 * Loads from a HDF5 file the data, filters by @ftype when available, using
 * the cluster information in @z_cluster, @ra_cluster and @dec_cluster.
 * 
 */
void 
nc_data_reduced_shear_cluster_mass_load_hdf5 (NcDataReducedShearClusterMass *drs, const gchar *hdf5_file, const gchar ftype, const gdouble z_cluster, const gdouble ra_cluster, const gdouble dec_cluster)
{
#ifdef HAVE_HDF5
	NcDataReducedShearClusterMassPrivate * const self = drs->priv;
	GHashTable *fdata     = g_hash_table_new_full (_g_long_hash, _g_long_equal, NULL, NULL); 
	NcmVector *vecs[5]    = {NULL, NULL, NULL, NULL, NULL};
	NcmMatrix *mats[2]    = {NULL, NULL};
	NcmVector *rh_vec     = NULL;
	GArray *gal_id        = NULL;
	GArray *photz_id      = NULL;
	GArray *filter        = NULL;
	GArray *gal_obs_array = NULL;
	NcmHDF5Table *photz_table;
	NcmHDF5Table *gal_table;
	hid_t h5f;
	herr_t ret;
	gint i;
	
	h5f = H5Fopen (hdf5_file, H5F_ACC_RDONLY, H5P_DEFAULT); 
	g_assert_cmpint (h5f, !=, -1);

	gal_table   = ncm_hdf5_table_new (h5f, "deepCoadd_meas");
	photz_table = ncm_hdf5_table_new (h5f, "zphot_ref");
		
	{
		const gchar *cols[] = {"coord_ra_deg", "coord_dec_deg", "ext_shapeHSM_HsmShapeRegauss_e1", "ext_shapeHSM_HsmShapeRegauss_e2"};
		for (i = 0; i < 4; i++)
		{
			vecs[i] = NULL;
			ncm_hdf5_table_read_col_as_vec (gal_table, cols[i], &vecs[i]);
		}

		if (ncm_hdf5_table_has_col (gal_table, "id"))
			ncm_hdf5_table_read_col_as_longarray (gal_table, "id", &gal_id);
		else
			g_error ("nc_data_reduced_shear_cluster_mass_load_hdf5: cannot find id in galaxy data table.");
		
		if (ncm_hdf5_table_has_col (gal_table, "filter"))
		{
			ncm_hdf5_table_read_col_as_fixstr (gal_table, "filter", &filter);
			
			for (i = 0; i < gal_table->nrecords; i++)
			{
				if (g_array_index (filter, gchar, i) == ftype)
				{
					glong *key_ptr = &g_array_index (gal_id, glong, i);
					glong key      = g_array_index (gal_id, glong, i);

					if (g_hash_table_contains (fdata, &key))
					{
						g_error ("nc_data_reduced_shear_cluster_mass_load_hdf5: duplicated key %ld\n", key);
					}
					else
					{
						g_hash_table_insert (fdata, key_ptr, GINT_TO_POINTER (i));
					}
				}
			}
		}
		else
		{
			for (i = 0; i < gal_table->nrecords; i++)
			{
				glong *key_ptr = &g_array_index (gal_id, glong, i);
				glong key      = g_array_index (gal_id, glong, i);

				if (g_hash_table_contains (fdata, &key))
				{
					g_error ("nc_data_reduced_shear_cluster_mass_load_hdf5: duplicated key %ld\n", key);
				}
				else
				{
					g_hash_table_insert (fdata, key_ptr, GINT_TO_POINTER (i));
				}
			}
		}

		if (ncm_hdf5_table_has_col (gal_table, "rh"))
			ncm_hdf5_table_read_col_as_vec (gal_table, "rh", &rh_vec);
		else
			rh_vec = NULL; /* Table does not contain galaxy size information - required to calibrate. */
	}

	if (g_hash_table_size (fdata) != photz_table->nrecords)
	{
		g_warning ("nc_data_reduced_shear_cluster_mass_load_hdf5: filtered galaxy data and photz tables do not match in size (%u %llu). Using only galaxy with data in both tables.", 
		           g_hash_table_size (fdata), photz_table->nrecords);
	}

	{
		const gchar *cols[]  = {"Z_BEST"};
		const gchar *acols[] = {"zbins", "pdz"};
		for (i = 0; i < 1; i++)
		{
			vecs[i+4] = NULL;
			ncm_hdf5_table_read_col_as_vec (photz_table, cols[i], &vecs[i+4]);
		}
		for (i = 0; i < 2; i++)
		{
			mats[i] = NULL;
			ncm_hdf5_table_read_col_as_mat (photz_table, acols[i], &mats[i]);
		}

		if (ncm_hdf5_table_has_col (photz_table, "id"))
			ncm_hdf5_table_read_col_as_longarray (photz_table, "id", &photz_id);
		else if (ncm_hdf5_table_has_col (photz_table, "objectId"))
			ncm_hdf5_table_read_col_as_longarray (photz_table, "objectId", &photz_id);
		else
			g_error ("nc_data_reduced_shear_cluster_mass_load_hdf5: cannot find id in photoz table.");
	}

	g_ptr_array_set_size ((GPtrArray *) self->photoz_array, 0);
	self->z_cluster   = z_cluster;
	self->ra_cluster  = ra_cluster;
	self->dec_cluster = dec_cluster;
	gal_obs_array     = g_array_new (FALSE, FALSE, sizeof (gdouble));

	for (i = 0; i < photz_table->nrecords; i++)
	{
		glong photz_id_i = g_array_index (photz_id, glong, i);
		if (!g_hash_table_contains (fdata, &photz_id_i))
			continue;
		else
		{
			const gint gindex = GPOINTER_TO_INT (g_hash_table_lookup (fdata, &photz_id_i));
			const gdouble ra  = ncm_vector_get (vecs[0], gindex);
			const gdouble dec = ncm_vector_get (vecs[1], gindex);
			const gdouble e1  = ncm_vector_get (vecs[2], gindex);
			const gdouble e2  = ncm_vector_get (vecs[3], gindex);
			const gdouble rh  = (rh_vec != NULL) ? ncm_vector_get (rh_vec, gindex) : 0.0;
			const gdouble zb  = ncm_vector_get (vecs[4], i);
			const guint ncols = ncm_matrix_ncols (mats[0]);

			gint first_nz     = -1;
			gint last_nz      = -1;
			gint j;

			g_assert_cmpuint (ncm_matrix_ncols (mats[0]), ==, ncm_matrix_ncols (mats[1]));
			g_assert_cmpuint (ncm_matrix_nrows (mats[0]), ==, ncm_matrix_nrows (mats[1]));

			for (j = 0; j < ncols; j++)
			{
				const gdouble Pz_j = ncm_matrix_get (mats[1], i, j);

				if (Pz_j > 0.0)
				{
					if (first_nz == -1)
					{
						first_nz = j;
					}
					last_nz = j;
				}
			}

			if (first_nz == -1)
			{
				g_error ("nc_data_reduced_shear_cluster_mass_load_hdf5: galaxy %d, empty P_z.", i);
			}

			{
				NcGalaxyRedshift *gz = NULL;
				gboolean is_delta = (first_nz == last_nz);
				if (is_delta)
				{
					gz = NC_GALAXY_REDSHIFT (nc_galaxy_redshift_spec_new (zb));
				}
				else
				{
					NcmVector *zv               = ncm_matrix_get_row (mats[0], i);
					NcmVector *Pzv              = ncm_matrix_get_row (mats[1], i);
					NcGalaxyRedshiftSpline *gzs = nc_galaxy_redshift_spline_new ();

					gz = NC_GALAXY_REDSHIFT (gzs);
					
					nc_galaxy_redshift_spline_init_from_vectors (gzs, zv, Pzv);
					
					ncm_vector_clear (&zv);
					ncm_vector_clear (&Pzv);
				}

				if (gsl_finite (ra) && gsl_finite (dec) && gsl_finite (e1) && gsl_finite (e2))
				{
					const gdouble r_arcmin = ncm_util_great_circle_distance (ra, dec, ra_cluster, dec_cluster) * 60.0;
					const gdouble posangle = -(0.5 * M_PI - ncm_util_position_angle (ra, dec, ra_cluster, dec_cluster));
					const gdouble cos2phi  = cos (2.0 * posangle);
					const gdouble sin2phi  = sin (2.0 * posangle);
					const gdouble g_obs    = - (e1 * cos2phi + e2 * sin2phi); 
					
					ncm_obj_array_add (self->photoz_array, G_OBJECT (gz));

					g_array_append_val (gal_obs_array, r_arcmin);
					g_array_append_val (gal_obs_array, g_obs);
					if (rh_vec != NULL)
						g_array_append_val (gal_obs_array, rh);
				}

				if (FALSE)
				{
					printf ("[%ld %ld] % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g %5d %5d", 
					        g_array_index (gal_id, glong, gindex), 
					        g_array_index (photz_id, glong, i), 
					        ra, dec, e1, e2, zb, first_nz, last_nz);

					if (nc_galaxy_redshift_has_dist (gz))
					{
						printf (" DISTR Z_BEST % 22.15g NINTERVALS %d", nc_galaxy_redshift_mode (gz), nc_galaxy_redshift_nintervals (gz));
						if (nc_galaxy_redshift_nintervals (gz) > 1)
						{
							gint k;
							for (k = 0; k < nc_galaxy_redshift_nintervals (gz); k++)
							{
								printf (" % 22.15g", nc_galaxy_redshift_interval_weight (gz, k));
							}
						}
						printf ("\n");
					}
					else
						printf (" DELTA Z_BEST % 22.15g\n", nc_galaxy_redshift_mode (gz));
				}

				nc_galaxy_redshift_free (gz);
			}
		}
	}

	ncm_matrix_clear (&self->gal_obs);
	if (rh_vec != NULL)
	{
		self->gal_obs = ncm_matrix_new_array (gal_obs_array, 3);
		self->has_rh  = TRUE;
	}
	else
	{
		self->gal_obs = ncm_matrix_new_array (gal_obs_array, 2);
		self->has_rh  = FALSE;
	}
	
	g_assert_cmpuint (self->photoz_array->len, ==, ncm_matrix_nrows (self->gal_obs));

	ncm_data_set_init (NCM_DATA (drs), TRUE);
	
	ncm_hdf5_table_free (gal_table);
	ncm_hdf5_table_free (photz_table);
	g_hash_table_unref (fdata);

	for (i = 0; i < 5; i++)
		ncm_vector_clear (&vecs[i]);

	ncm_vector_clear (&rh_vec);
	
	for (i = 0; i < 2; i++)
		ncm_matrix_clear (&mats[i]);

	g_array_unref (gal_id);
	g_array_unref (photz_id);
	g_array_unref (gal_obs_array);
	if (filter != NULL)
		g_array_unref (filter);

	ret = H5Fclose (h5f);
	g_assert_cmpint (ret, !=, -1);
#else  /* HAVE_HDF5 */
	g_error ("nc_data_reduced_shear_cluster_mass_load_hdf5: numcosmo built without support for HDF5 files.");
#endif /* HAVE_HDF5 */
}
