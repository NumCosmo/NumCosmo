/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_galaxy_redshift_spline.c
 *
 *  Thu April 19 14:38:59 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br>, <pennalima@gmail.com>
 ****************************************************************************/
/*
 * nc_galaxy_redshift_spline.c
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_galaxy_redshift_spline
 * @title: NcGalaxyRedshiftSpline
 * @short_description: Class describing spectroscopic galaxy redshifts.
 *
 * Class used to define a generic galaxy redshift probability distribution 
 * $P^g_i(z)$ using splines for each disjoint interval $i$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_vector.h"
#include "math/ncm_obj_array.h"
#include "math/ncm_stats_dist1d_spline.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline_gsl.h"
#include "lss/nc_galaxy_redshift_spline.h"

struct _NcGalaxyRedshiftSplinePrivate
{
	gdouble z_best;
	GArray *normas;
	NcmObjArray *dists;
};

enum
{
	PROP_0,
	PROP_Z_BEST,
	PROP_DISTS,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxyRedshiftSpline, nc_galaxy_redshift_spline, NC_TYPE_GALAXY_REDSHIFT);

static void
nc_galaxy_redshift_spline_init (NcGalaxyRedshiftSpline *gzs)
{
	NcGalaxyRedshiftSplinePrivate * const self = gzs->priv = nc_galaxy_redshift_spline_get_instance_private (gzs);

	self->z_best = 0.0;
	self->normas = g_array_new (TRUE, TRUE, sizeof (gdouble));
	
	self->dists  = ncm_obj_array_new ();
}

static void
_nc_galaxy_redshift_spline_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
	NcGalaxyRedshiftSpline *gzs = NC_GALAXY_REDSHIFT_SPLINE (object);
	NcGalaxyRedshiftSplinePrivate * const self = gzs->priv;
	g_return_if_fail (NC_IS_GALAXY_REDSHIFT_SPLINE (object));

	switch (prop_id)
	{
		case PROP_Z_BEST:
			nc_galaxy_redshift_spline_set_z_best (gzs, g_value_get_double (value));
			break;
		case PROP_DISTS:
		{
			NcmObjArray *dists = g_value_get_boxed (value);
			gdouble norma_t    = 0.0;
			gint i;
			
			g_clear_pointer (&self->dists, ncm_obj_array_unref);
			self->dists = ncm_obj_array_ref (dists);

			g_array_set_size (self->normas, 0);
			
			for (i = 0; i < dists->len; i++)
			{
				NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (ncm_obj_array_peek (self->dists, i));
				ncm_stats_dist1d_prepare (sd1);
				
				norma_t	+= ncm_stats_dist1d_eval_norma (sd1);
			}

			for (i = 0; i < dists->len; i++)
			{
				NcmStatsDist1d *sd1   = NCM_STATS_DIST1D (ncm_obj_array_peek (self->dists, i));				
				const gdouble norma_i = ncm_stats_dist1d_eval_norma (sd1) / norma_t;
				g_array_append_val (self->normas, norma_i);
			}

			break;
		}
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
_nc_galaxy_redshift_spline_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
	NcGalaxyRedshiftSpline *gzs = NC_GALAXY_REDSHIFT_SPLINE (object);
	NcGalaxyRedshiftSplinePrivate * const self = gzs->priv;
	g_return_if_fail (NC_IS_GALAXY_REDSHIFT_SPLINE (object));

	switch (prop_id)
	{
		case PROP_Z_BEST:
			g_value_set_double (value, nc_galaxy_redshift_spline_get_z_best (gzs));
			break;
		case PROP_DISTS:
			g_value_set_boxed (value, self->dists);
			break;
		default:
			G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
			break;
	}
}

static void
_nc_galaxy_redshift_spline_dispose (GObject *object)
{
	NcGalaxyRedshiftSpline *gzs = NC_GALAXY_REDSHIFT_SPLINE (object);
	NcGalaxyRedshiftSplinePrivate * const self = gzs->priv;

	g_clear_pointer (&self->normas, g_array_unref);
	ncm_obj_array_clear (&self->dists);

	/* Chain up : end */
	G_OBJECT_CLASS (nc_galaxy_redshift_spline_parent_class)->dispose (object);
}

static void
_nc_galaxy_redshift_spline_finalize (GObject *object)
{
	/*NcGalaxyRedshiftSpline *gzs = NC_GALAXY_REDSHIFT_SPLINE (object);*/
	/*NcGalaxyRedshiftSplinePrivate * const self = gzs->priv;*/


	/* Chain up : end */
	G_OBJECT_CLASS (nc_galaxy_redshift_spline_parent_class)->finalize (object);
}

static gboolean _nc_galaxy_redshift_spline_has_dist (NcGalaxyRedshift *gz);
static gdouble _nc_galaxy_redshift_spline_mode (NcGalaxyRedshift *gz);
static guint _nc_galaxy_redshift_spline_nintervals (NcGalaxyRedshift *gz);
static gdouble _nc_galaxy_redshift_spline_interval_weight (NcGalaxyRedshift *gz, const guint di);
static void _nc_galaxy_redshift_spline_pdf_limits (NcGalaxyRedshift *gz, const guint di, gdouble *zmin, gdouble *zmax);
static gdouble _nc_galaxy_redshift_spline_pdf (NcGalaxyRedshift *gz, const guint di, const gdouble z);
static gdouble _nc_galaxy_redshift_spline_gen (NcGalaxyRedshift *gz, NcmRNG *rng);

static void
nc_galaxy_redshift_spline_class_init (NcGalaxyRedshiftSplineClass *klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);
	NcGalaxyRedshiftClass *gz_class = NC_GALAXY_REDSHIFT_CLASS (klass);

	object_class->set_property = &_nc_galaxy_redshift_spline_set_property;
	object_class->get_property = &_nc_galaxy_redshift_spline_get_property;
	object_class->dispose      = &_nc_galaxy_redshift_spline_dispose;
	object_class->finalize     = &_nc_galaxy_redshift_spline_finalize;

	g_object_class_install_property (object_class,
	                                 PROP_Z_BEST,
	                                 g_param_spec_double ("z-best",
	                                                      NULL,
	                                                      "Distributions mode",
	                                                      0.0, G_MAXDOUBLE, 0.0,
	                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_DISTS,
	                                 g_param_spec_boxed ("dists",
	                                                     NULL,
	                                                     "Distribution objects",
	                                                     NCM_TYPE_OBJ_ARRAY,
	                                                     G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));	
	
	gz_class->has_dist        = &_nc_galaxy_redshift_spline_has_dist;
	gz_class->mode            = &_nc_galaxy_redshift_spline_mode;
	gz_class->nintervals      = &_nc_galaxy_redshift_spline_nintervals;
	gz_class->interval_weight = &_nc_galaxy_redshift_spline_interval_weight;
	gz_class->pdf_limits      = &_nc_galaxy_redshift_spline_pdf_limits;
  gz_class->pdf             = &_nc_galaxy_redshift_spline_pdf;
  gz_class->gen             = &_nc_galaxy_redshift_spline_gen;

}

static gboolean 
_nc_galaxy_redshift_spline_has_dist (NcGalaxyRedshift *gz)
{
	return TRUE;
}

static gdouble 
_nc_galaxy_redshift_spline_mode (NcGalaxyRedshift *gz)
{
	NcGalaxyRedshiftSpline *gzs = NC_GALAXY_REDSHIFT_SPLINE (gz);
	NcGalaxyRedshiftSplinePrivate * const self = gzs->priv;

	return self->z_best;
}

static guint 
_nc_galaxy_redshift_spline_nintervals (NcGalaxyRedshift *gz)
{
	NcGalaxyRedshiftSpline *gzs = NC_GALAXY_REDSHIFT_SPLINE (gz);
	NcGalaxyRedshiftSplinePrivate * const self = gzs->priv;

	return self->dists->len;
}

static gdouble 
_nc_galaxy_redshift_spline_interval_weight (NcGalaxyRedshift *gz, const guint di)
{
	NcGalaxyRedshiftSpline *gzs = NC_GALAXY_REDSHIFT_SPLINE (gz);
	NcGalaxyRedshiftSplinePrivate * const self = gzs->priv;

	g_assert_cmpuint (di, <, self->normas->len);
	
	return g_array_index (self->normas, gdouble, di);
}

static void 
_nc_galaxy_redshift_spline_pdf_limits (NcGalaxyRedshift *gz, const guint di, gdouble *zmin, gdouble *zmax)
{
	NcGalaxyRedshiftSpline *gzs = NC_GALAXY_REDSHIFT_SPLINE (gz);
	NcGalaxyRedshiftSplinePrivate * const self = gzs->priv;

	{
		NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (ncm_obj_array_peek (self->dists, di));
		zmin[0] = ncm_stats_dist1d_get_xi (sd1);
		zmax[0] = ncm_stats_dist1d_get_xf (sd1);
	}
}

static gdouble 
_nc_galaxy_redshift_spline_pdf (NcGalaxyRedshift *gz, const guint di, const gdouble z)
{
	NcGalaxyRedshiftSpline *gzs = NC_GALAXY_REDSHIFT_SPLINE (gz);
	NcGalaxyRedshiftSplinePrivate * const self = gzs->priv;

	{
		NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (ncm_obj_array_peek (self->dists, di));
		return g_array_index (self->normas, gdouble, di) * ncm_stats_dist1d_eval_p (sd1, z);
	}
}

static gdouble 
_nc_galaxy_redshift_spline_gen (NcGalaxyRedshift *gz, NcmRNG *rng)
{
	NcGalaxyRedshiftSpline *gzs = NC_GALAXY_REDSHIFT_SPLINE (gz);
	NcGalaxyRedshiftSplinePrivate * const self = gzs->priv;
	guint *n = g_new0 (guint, self->normas->len);
	gboolean done = FALSE;
	gdouble z = 0.0;
	guint i;

	gsl_ran_multinomial (rng->r, self->normas->len, 1, (gdouble *)self->normas->data, n);
	for (i = 0; i < self->normas->len; i++)
	{
		if (n[i] == 1)
		{
			NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (ncm_obj_array_peek (self->dists, i));
			g_assert (!done);
			z = ncm_stats_dist1d_gen (sd1, rng); 
			done = TRUE;
		}
	}

	g_assert (done);
	g_free (n);

	return z;
}	

/**
 * nc_galaxy_redshift_spline_new:
 * 
 * Creates a new #NcGalaxyRedshiftSpline.
 *
 * Returns: (transfer full): The newly created #NcGalaxyRedshiftSpline.
 */
NcGalaxyRedshiftSpline *
nc_galaxy_redshift_spline_new (void)
{
	NcGalaxyRedshiftSpline *gzs = g_object_new (NC_TYPE_GALAXY_REDSHIFT_SPLINE, 
	                                            NULL);

	return gzs;
}

/**
 * nc_galaxy_redshift_spline_ref:
 * @gzs: a #NcGalaxyRedshiftSpline
 *
 * Increase the reference of @gzs by one.
 *
 * Returns: (transfer full): @gzs.
 */
NcGalaxyRedshiftSpline *
nc_galaxy_redshift_spline_ref (NcGalaxyRedshiftSpline *gzs)
{
  return g_object_ref (gzs);
}

/**
 * nc_galaxy_redshift_spline_free:
 * @gzs: a #NcGalaxyRedshiftSpline
 *
 * Decrease the reference count of @gzs by one.
 *
 */
void
nc_galaxy_redshift_spline_free (NcGalaxyRedshiftSpline *gzs)
{
  g_object_unref (gzs);
}

/**
 * nc_galaxy_redshift_spline_clear:
 * @gzs: a #NcGalaxyRedshiftSpline
 *
 * Decrease the reference count of @gzs by one, and sets the pointer *@gzs to
 * NULL.
 *
 */
void
nc_galaxy_redshift_spline_clear (NcGalaxyRedshiftSpline **gzs)
{
  g_clear_object (gzs);
}


/**
 * nc_galaxy_redshift_spline_set_z_best:
 * @gzs: a #NcGalaxyRedshiftSpline
 * @z_best: the mode of the redshift distribution
 * 
 * Sets the mode of the redshift distribution @z_best.
 *
 */
void
nc_galaxy_redshift_spline_set_z_best (NcGalaxyRedshiftSpline *gzs, const gdouble z_best)
{
	NcGalaxyRedshiftSplinePrivate * const self = gzs->priv;
	
	self->z_best = z_best;
}

/**
 * nc_galaxy_redshift_spline_get_z_best:
 * @gzs: a #NcGalaxyRedshiftSpline
 * 
 * Gets $z_\mathrm{best}$, the mode of the redshift distribution.
 *
 * Returns: $z_\mathrm{best}$.
 */
gdouble
nc_galaxy_redshift_spline_get_z_best (NcGalaxyRedshiftSpline *gzs)
{
	NcGalaxyRedshiftSplinePrivate * const self = gzs->priv;
	
	return self->z_best;
}

/**
 * nc_galaxy_redshift_spline_init_from_vectors:
 * @gzs: a #NcGalaxyRedshiftSpline
 * @zv: a #NcmVector containing the redshift knots
 * @Pzv: a #NcmVector containing the $P(z)$ at the redshift knots
 * 
 * Initialize @gzs (cleaning any previous data) using the data from
 * @z and @Pzv.
 *
 */
void
nc_galaxy_redshift_spline_init_from_vectors (NcGalaxyRedshiftSpline *gzs, NcmVector *zv, NcmVector *Pzv)
{
	NcGalaxyRedshiftSplinePrivate * const self = gzs->priv;
	gint first_nz    = -1;
	guint len        = ncm_vector_len (Pzv);
	GArray *z_a      = NULL;
	GArray *m2lnPz_a = NULL;
	gdouble z_best   = 0.0;
	gdouble Pz_best  = 0.0;
	gint i;

	g_assert_cmpuint (ncm_vector_len (zv), ==, ncm_vector_len (Pzv));

	g_ptr_array_set_size ((GPtrArray *)self->dists, 0);

	for (i = 0; i < len; i++)
	{
		const gdouble z_i  = ncm_vector_get (zv, i);
		const gdouble Pz_i = ncm_vector_get (Pzv, i);
		g_assert_cmpfloat (Pz_i, >=, 0.0);

		if (Pz_i > Pz_best)
		{
			z_best  = z_i;
			Pz_best = Pz_i;
		}

		if (Pz_i > 0.0)
		{
			const gdouble m2lnPz_i = -2.0 * log (Pz_i);
			
			if (first_nz == -1)
			{
				first_nz = i;

				z_a      = g_array_new (TRUE, TRUE, sizeof (gdouble));
				m2lnPz_a = g_array_new (TRUE, TRUE, sizeof (gdouble));

				if (i > 0)
				{
					const gdouble z_im1      = ncm_vector_get (zv, i - 1);
					const gdouble m2lnPz_im1 = -2.0 * NC_GALAXY_REDSHIFT_SPLINE_LKNOT_DROP + m2lnPz_i;

					g_array_append_val (z_a, z_im1);
					g_array_append_val (m2lnPz_a, m2lnPz_im1);
				}
			}

			g_array_append_val (z_a, z_i);
			g_array_append_val (m2lnPz_a, m2lnPz_i);
		}
		else
		{
			if (first_nz != -1)
			{
				NcmStatsDist1dSpline *dist = NULL;
				const gdouble m2lnPz_i     = -2.0 * (NC_GALAXY_REDSHIFT_SPLINE_LKNOT_DROP + log (ncm_vector_get (Pzv, i - 1)));

				g_array_append_val (z_a, z_i);
				g_array_append_val (m2lnPz_a, m2lnPz_i);

				{
					NcmSpline *s = (z_a->len >= 6) ? NCM_SPLINE (ncm_spline_cubic_notaknot_new ()) : NCM_SPLINE (ncm_spline_gsl_new (gsl_interp_polynomial));
					ncm_spline_set_array (s, z_a, m2lnPz_a, TRUE);

					dist = ncm_stats_dist1d_spline_new (s);

					/*ncm_vector_log_vals (s->xv, "XV", "% 22.15g", TRUE);*/
					/*ncm_vector_log_vals (s->yv, "YV", "% 22.15g", TRUE);*/
					ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (dist));
					
					ncm_spline_free (s);
				}

				ncm_obj_array_add (self->dists, G_OBJECT (dist));

				ncm_stats_dist1d_free (NCM_STATS_DIST1D (dist));
				g_array_unref (z_a);
				g_array_unref (m2lnPz_a);
				
				first_nz = -1;
			}
		}
	}

	if (first_nz != -1)
	{
		NcmStatsDist1dSpline *dist = NULL;
		NcmSpline *s               = (z_a->len >= 6) ? NCM_SPLINE (ncm_spline_cubic_notaknot_new ()) : NCM_SPLINE (ncm_spline_gsl_new (gsl_interp_polynomial));

		ncm_spline_set_array (s, z_a, m2lnPz_a, TRUE);
		dist = ncm_stats_dist1d_spline_new (s);
		ncm_stats_dist1d_prepare (NCM_STATS_DIST1D (dist));

		ncm_obj_array_add (self->dists, G_OBJECT (dist));

		ncm_stats_dist1d_free (NCM_STATS_DIST1D (dist));
		ncm_spline_free (s);
		g_array_unref (z_a);
		g_array_unref (m2lnPz_a);

		first_nz = -1;
	}
	
	if (self->dists->len == 0)
	{
		g_error ("nc_galaxy_redshift_spline_init_from_vectors: empty P_z.");
	}

	{
		gdouble norma_t = 0.0;
		g_array_set_size (self->normas, 0);

		for (i = 0; i < self->dists->len; i++)
		{
			NcmStatsDist1d *sd1 = NCM_STATS_DIST1D (ncm_obj_array_peek (self->dists, i));				
			norma_t	+= ncm_stats_dist1d_eval_norma (sd1);
		}
		
		for (i = 0; i < self->dists->len; i++)
		{
			NcmStatsDist1d *sd1   = NCM_STATS_DIST1D (ncm_obj_array_peek (self->dists, i));				
			const gdouble norma_i = ncm_stats_dist1d_eval_norma (sd1) / norma_t;
			g_array_append_val (self->normas, norma_i);
		}

		self->z_best = z_best;
	}
}
