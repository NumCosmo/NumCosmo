/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            nc_galaxy_redshift_gauss.c
 *
 *  Thu July 27 14:38:59 2020
 *  Copyright  2020  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_gauss.c
 * Copyright (C) 2020 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
 * Copyright (C) 2020 Mariana Penna Lima <pennalima@gmail.com>
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
 * SECTION:nc_galaxy_redshift_gauss
 * @title: NcGalaxyRedshiftGauss
 * @short_description: Class describing Gaussian photometric galaxy redshifts.
 * @stability: Unstable
 *
 * Class used to define a generic galaxy redshift probability distribution
 * $P^g_i(z)$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_galaxy_redshift_gauss.h"
#include "math/ncm_vector.h"
#include "math/ncm_obj_array.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

struct _NcGalaxyRedshiftGaussPrivate
{
  guint len;
  NcmMatrix *obs;
  GPtrArray *fpws;
  GPtrArray *nodes;
  GPtrArray *weights;
  GArray *nnodes;
  GArray *norma;
  gboolean constructed;
};

enum
{
  PROP_0,
  PROP_OBS,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyRedshiftGauss, nc_galaxy_redshift_gauss, NC_TYPE_GALAXY_REDSHIFT);

static void
nc_galaxy_redshift_gauss_init (NcGalaxyRedshiftGauss *gzg)
{
  NcGalaxyRedshiftGaussPrivate * const self = gzg->priv = nc_galaxy_redshift_gauss_get_instance_private (gzg);
  
  self->len         = 0;
  self->obs         = NULL;
  self->fpws        = g_ptr_array_new ();
  self->nodes       = g_ptr_array_new ();
  self->weights     = g_ptr_array_new ();
  self->nnodes      = g_array_new (FALSE, FALSE, sizeof (guint));
  self->norma       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->constructed = FALSE;
  
  g_ptr_array_set_free_func (self->fpws, (GDestroyNotify) gsl_integration_fixed_free);
}

static void
_nc_galaxy_redshift_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyRedshiftGauss *gzg = NC_GALAXY_REDSHIFT_GAUSS (object);
  
  /*NcGalaxyRedshiftGaussPrivate * const self = gzg->priv;*/
  
  g_return_if_fail (NC_IS_GALAXY_REDSHIFT_GAUSS (object));
  
  switch (prop_id)
  {
    case PROP_OBS:
      nc_galaxy_redshift_gauss_set_obs (gzg, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_redshift_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyRedshiftGauss *gzg = NC_GALAXY_REDSHIFT_GAUSS (object);
  
  /* NcGalaxyRedshiftGaussPrivate * const self = gzg->priv; */
  
  g_return_if_fail (NC_IS_GALAXY_REDSHIFT_GAUSS (object));
  
  switch (prop_id)
  {
    case PROP_OBS:
      g_value_set_object (value, nc_galaxy_redshift_gauss_peek_obs (gzg));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_redshift_gauss_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_galaxy_redshift_gauss_parent_class)->constructed (object);
  {
    NcGalaxyRedshiftGauss *gzg                = NC_GALAXY_REDSHIFT_GAUSS (object);
    NcGalaxyRedshiftGaussPrivate * const self = gzg->priv;
    
    self->constructed = TRUE;
  }
}

static void
_nc_galaxy_redshift_gauss_dispose (GObject *object)
{
  NcGalaxyRedshiftGauss *gzg                = NC_GALAXY_REDSHIFT_GAUSS (object);
  NcGalaxyRedshiftGaussPrivate * const self = gzg->priv;
  
  g_clear_pointer (&self->fpws,    g_ptr_array_unref);
  g_clear_pointer (&self->nodes,   g_ptr_array_unref);
  g_clear_pointer (&self->weights, g_ptr_array_unref);
  
  g_clear_pointer (&self->nnodes, g_array_unref);
  g_clear_pointer (&self->norma,  g_array_unref);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_redshift_gauss_parent_class)->dispose (object);
}

static void
_nc_galaxy_redshift_gauss_finalize (GObject *object)
{
  /*NcGalaxyRedshiftGauss *gzg = NC_GALAXY_REDSHIFT_GAUSS (object);*/
  /*NcGalaxyRedshiftGaussPrivate * const self = gzg->priv;*/
  
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_redshift_gauss_parent_class)->finalize (object);
}

static gboolean _nc_galaxy_redshift_gauss_has_dist (NcGalaxyRedshift *gz);
static gdouble _nc_galaxy_redshift_gauss_mode (NcGalaxyRedshift *gz);
static guint _nc_galaxy_redshift_gauss_nintervals (NcGalaxyRedshift *gz);
static gdouble _nc_galaxy_redshift_gauss_interval_weight (NcGalaxyRedshift *gz, const guint di);
static void _nc_galaxy_redshift_gauss_pdf_limits (NcGalaxyRedshift *gz, const guint di, gdouble *zmin, gdouble *zmax);
static gdouble _nc_galaxy_redshift_gauss_pdf (NcGalaxyRedshift *gz, const guint di, const gdouble z);
static gdouble _nc_galaxy_redshift_gauss_gen (NcGalaxyRedshift *gz, NcmRNG *rng);
static gdouble _nc_galaxy_redshift_gauss_compute_mean_m2lnf (NcGalaxyRedshift *gz, guint gal_i, NcGalaxyRedshiftF m2lnf, gpointer userdata);
static guint _nc_galaxy_redshift_gauss_len (NcGalaxyRedshift *gz);

static void
nc_galaxy_redshift_gauss_class_init (NcGalaxyRedshiftGaussClass *klass)
{
  GObjectClass *object_class      = G_OBJECT_CLASS (klass);
  NcGalaxyRedshiftClass *gz_class = NC_GALAXY_REDSHIFT_CLASS (klass);
  
  object_class->set_property = &_nc_galaxy_redshift_gauss_set_property;
  object_class->get_property = &_nc_galaxy_redshift_gauss_get_property;
  object_class->constructed  = &_nc_galaxy_redshift_gauss_constructed;
  object_class->dispose      = &_nc_galaxy_redshift_gauss_dispose;
  object_class->finalize     = &_nc_galaxy_redshift_gauss_finalize;
  
  g_object_class_install_property (object_class,
                                   PROP_OBS,
                                   g_param_spec_object ("obs",
                                                        NULL,
                                                        "Redshift observables",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  gz_class->has_dist           = &_nc_galaxy_redshift_gauss_has_dist;
  gz_class->mode               = &_nc_galaxy_redshift_gauss_mode;
  gz_class->nintervals         = &_nc_galaxy_redshift_gauss_nintervals;
  gz_class->interval_weight    = &_nc_galaxy_redshift_gauss_interval_weight;
  gz_class->pdf_limits         = &_nc_galaxy_redshift_gauss_pdf_limits;
  gz_class->pdf                = &_nc_galaxy_redshift_gauss_pdf;
  gz_class->gen                = &_nc_galaxy_redshift_gauss_gen;
  gz_class->compute_mean_m2lnf = &_nc_galaxy_redshift_gauss_compute_mean_m2lnf;
  gz_class->len                = &_nc_galaxy_redshift_gauss_len;
}

static gboolean
_nc_galaxy_redshift_gauss_has_dist (NcGalaxyRedshift *gz)
{
  return TRUE;
}

static gdouble
_nc_galaxy_redshift_gauss_mode (NcGalaxyRedshift *gz)
{
  /*NcGalaxyRedshiftGauss *gzg                = NC_GALAXY_REDSHIFT_GAUSS (gz);*/
  /*NcGalaxyRedshiftGaussPrivate * const self = gzg->priv;*/
  
  return 0.0;
}

static guint
_nc_galaxy_redshift_gauss_nintervals (NcGalaxyRedshift *gz)
{
  return 1;
}

static gdouble
_nc_galaxy_redshift_gauss_interval_weight (NcGalaxyRedshift *gz, const guint di)
{
  return 1.0;
}

static void
_nc_galaxy_redshift_gauss_pdf_limits (NcGalaxyRedshift *gz, const guint di, gdouble *zmin, gdouble *zmax)
{
  /*NcGalaxyRedshiftGauss *gzg                = NC_GALAXY_REDSHIFT_GAUSS (gz);*/
  /*NcGalaxyRedshiftGaussPrivate * const self = gzg->priv;*/
  
  zmin[0] = 0.0;
  zmax[0] = 1.0;
}

static gdouble
_nc_galaxy_redshift_gauss_pdf (NcGalaxyRedshift *gz, const guint di, const gdouble z)
{
  /*NcGalaxyRedshiftGauss *gzg                = NC_GALAXY_REDSHIFT_GAUSS (gz);*/
  /*NcGalaxyRedshiftGaussPrivate * const self = gzg->priv;*/
  
  return 0.0;
}

static gdouble
_nc_galaxy_redshift_gauss_gen (NcGalaxyRedshift *gz, NcmRNG *rng)
{
/*
 *  NcGalaxyRedshiftGauss *gzg                = NC_GALAXY_REDSHIFT_GAUSS (gz);
 *  NcGalaxyRedshiftGaussPrivate * const self = gzg->priv;
 *  gdouble z;
 *
 *  do {
 *   z = ncm_rng_gaussian_gen (rng, self->z_obs, self->sigma_z);
 *  } while (z < 0.0);
 */
  
  return 0.0;
}

static gdouble
_nc_galaxy_redshift_gauss_compute_mean_m2lnf (NcGalaxyRedshift *gz, guint gal_i, NcGalaxyRedshiftF m2lnf, gpointer userdata)
{
  NcGalaxyRedshiftGauss *gzg                = NC_GALAXY_REDSHIFT_GAUSS (gz);
  NcGalaxyRedshiftGaussPrivate * const self = gzg->priv;
  const guint nnodes                        = g_array_index (self->nnodes, guint, gal_i);
  const gdouble norma                       = g_array_index (self->norma, gdouble, gal_i);
  const gdouble *nodes                      = g_ptr_array_index (self->nodes, gal_i);
  const gdouble *weights                    = g_ptr_array_index (self->weights, gal_i);
  gdouble res                               = 0.0;
  gint i;
  
  for (i = 0; i < nnodes; i++)
  {
    res += weights[i] * exp (-0.5 * m2lnf (nodes[i], userdata));
  }
  
  return -2.0 * log (res / norma);
}

static guint
_nc_galaxy_redshift_gauss_len (NcGalaxyRedshift *gz)
{
  NcGalaxyRedshiftGauss *gzg                = NC_GALAXY_REDSHIFT_GAUSS (gz);
  NcGalaxyRedshiftGaussPrivate * const self = gzg->priv;
  
  return self->len;
}

/**
 * nc_galaxy_redshift_gauss_new:
 *
 * Creates a new #NcGalaxyRedshiftGauss.
 *
 * Returns: (transfer full): The newly created #NcGalaxyRedshiftGauss.
 */
NcGalaxyRedshiftGauss *
nc_galaxy_redshift_gauss_new (void)
{
  NcGalaxyRedshiftGauss *gzg = g_object_new (NC_TYPE_GALAXY_REDSHIFT_GAUSS,
                                             NULL);
  
  return gzg;
}

/**
 * nc_galaxy_redshift_gauss_ref:
 * @gzg: a #NcGalaxyRedshiftGauss
 *
 * Increase the reference of @gzg by one.
 *
 * Returns: (transfer full): @gzg.
 */
NcGalaxyRedshiftGauss *
nc_galaxy_redshift_gauss_ref (NcGalaxyRedshiftGauss *gzg)
{
  return g_object_ref (gzg);
}

/**
 * nc_galaxy_redshift_gauss_free:
 * @gzg: a #NcGalaxyRedshiftGauss
 *
 * Decrease the reference count of @gzg by one.
 *
 */
void
nc_galaxy_redshift_gauss_free (NcGalaxyRedshiftGauss *gzg)
{
  g_object_unref (gzg);
}

/**
 * nc_galaxy_redshift_gauss_clear:
 * @gzg: a #NcGalaxyRedshiftGauss
 *
 * Decrease the reference count of @gzg by one, and sets the pointer *@gzg to
 * NULL.
 *
 */
void
nc_galaxy_redshift_gauss_clear (NcGalaxyRedshiftGauss **gzg)
{
  g_clear_object (gzg);
}

/**
 * nc_galaxy_redshift_gauss_set_obs:
 * @gzg: a #NcGalaxyRedshiftGauss
 * @obs: a #NcmMatrix
 *
 * Sets the observables of the redshift distribution through
 * the matrix @obs. @obs should contain two columns where the
 * first contain the mean redshift $\bar{z}$ and the second
 * the standard deviation $\sigma_z$.
 *
 */
void
nc_galaxy_redshift_gauss_set_obs (NcGalaxyRedshiftGauss *gzg, NcmMatrix *obs)
{
  NcGalaxyRedshiftGaussPrivate * const self = gzg->priv;
  guint i;
  
  g_assert_cmpuint (ncm_matrix_ncols (obs), ==, 2);
  g_assert_cmpuint (ncm_matrix_nrows (obs), >, 0);
  
  ncm_matrix_clear (&self->obs);
  
  self->len = ncm_matrix_nrows (obs);
  self->obs = ncm_matrix_ref (obs);
  
  g_ptr_array_set_size (self->fpws, 0);
  g_ptr_array_set_size (self->fpws,    self->len);
  g_ptr_array_set_size (self->nodes,   self->len);
  g_ptr_array_set_size (self->weights, self->len);
  
  g_array_set_size (self->nnodes, self->len);
  g_array_set_size (self->norma,  self->len);
  
  for (i = 0; i < self->len; i++)
  {
    const gdouble z_obs                   = ncm_matrix_get (obs, i, 0);
    const gdouble sigma_z                 = ncm_matrix_get (obs, i, 1);
    const gint N                          = 120;
    gsl_integration_fixed_workspace *fpws = gsl_integration_fixed_alloc (gsl_integration_fixed_hermite, N, z_obs, 0.5 / (sigma_z * sigma_z), 0.0, 0.0);
    gdouble *nodes                        = gsl_integration_fixed_nodes (fpws);
    gdouble *weights                      = gsl_integration_fixed_weights (fpws);
    const gdouble norma                   = sqrt (2.0 * M_PI) * sigma_z * 0.5 * (1.0 + erf (z_obs / (sqrt (2.0) * sigma_z)));
    gint j;
    
    for (j = 0; j < N; j++)
    {
      if (nodes[j] > 0)
        break;
    }
    
    g_ptr_array_index (self->fpws, i)    = fpws;
    g_ptr_array_index (self->nodes, i)   = &nodes[j];
    g_ptr_array_index (self->weights, i) = &weights[j];
    
    g_array_index (self->nnodes, guint, i)   = N - j;
    g_array_index (self->norma,  gdouble, i) = norma;
  }
}

/**
 * nc_galaxy_redshift_gauss_peek_obs:
 * @gzg: a #NcGalaxyRedshiftGauss
 *
 * Gets observations matrix.
 *
 * Returns: (transfer none): observations matrix.
 */
NcmMatrix *
nc_galaxy_redshift_gauss_peek_obs (NcGalaxyRedshiftGauss *gzg)
{
  NcGalaxyRedshiftGaussPrivate * const self = gzg->priv;
  
  return self->obs;
}

