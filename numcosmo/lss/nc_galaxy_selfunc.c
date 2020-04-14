/***************************************************************************
 *            nc_galaxy_selfunc.c
 *
 *  Wed March 14 16:30:36 2018
 *  Copyright  2018 Fernando de Simoni 
 *  <fsimoni@id.uff.br>
 ****************************************************************************/
/*
 * nc_galaxy_selfunc.h
 * Copyright (C) 2018 Fernando de Simoni <fsimoni@id.uff.brr>
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
 * SECTION:nc_galaxy_selfunc
 * @title: NcGalaxySelfunc
 * @short_description: Galaxy phenomelogical selection function.
 * @include: nc_galaxy_selfunc.h
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_galaxy_selfunc.h"
#include "math/ncm_vector.h"
#include "math/ncm_spline.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_obj_array.h"

struct _NcGalaxySelfuncPrivate
{
  guint nshells;

  gdouble zmin_all;
  gdouble zmax_all;

  NcmVector *norm;
  NcmVector *zmin;
  NcmVector *zmean;
  NcmVector *zmax;

  NcmObjArray *dNdz_a;
};

enum
{
  PROP_0,
  PROP_NSHELLS,
	PROP_SHELL_SPLINES,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcGalaxySelfunc, nc_galaxy_selfunc, G_TYPE_OBJECT);

static void
nc_galaxy_selfunc_init (NcGalaxySelfunc *gsf)
{
	NcGalaxySelfuncPrivate * const self = gsf->priv = nc_galaxy_selfunc_get_instance_private (gsf);
	
  self->nshells  = 0;
  self->zmin_all = 0;
  self->zmax_all = 0;
  self->norm     = NULL;
  self->zmin     = NULL;
  self->zmean    = NULL;
  self->zmax     = NULL;
  self->dNdz_a   = ncm_obj_array_new ();
}

static void
_nc_galaxy_selfunc_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxySelfunc *gsf = NC_GALAXY_SELFUNC (object);
  g_return_if_fail (NC_IS_GALAXY_SELFUNC (object));

  switch (prop_id)
  {
    case PROP_NSHELLS:
      nc_galaxy_selfunc_set_nshells (gsf, g_value_get_uint (value));
      break;
    case PROP_SHELL_SPLINES:
      nc_galaxy_selfunc_set_shell_splines (gsf, g_value_get_boxed (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_selfunc_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxySelfunc *gsf = NC_GALAXY_SELFUNC (object);
  g_return_if_fail (NC_IS_GALAXY_SELFUNC (object));

  switch (prop_id)
  {
    case PROP_NSHELLS:
      g_value_set_uint (value, nc_galaxy_selfunc_get_nshells (gsf));
      break;
    case PROP_SHELL_SPLINES:
      g_value_take_boxed (value, nc_galaxy_selfunc_get_shell_splines (gsf));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_galaxy_selfunc_dispose (GObject *object)
{
  NcGalaxySelfunc *gsf = NC_GALAXY_SELFUNC (object);
	NcGalaxySelfuncPrivate * const self = gsf->priv;

  ncm_vector_clear (&self->norm);
  ncm_vector_clear (&self->zmin);
  ncm_vector_clear (&self->zmean);
  ncm_vector_clear (&self->zmax);
	ncm_obj_array_clear (&self->dNdz_a);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_selfunc_parent_class)->dispose (object);
}

static void
_nc_galaxy_selfunc_finalize (GObject *object)
{
	
  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_selfunc_parent_class)->finalize (object);
}

static void
nc_galaxy_selfunc_class_init (NcGalaxySelfuncClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_galaxy_selfunc_set_property;
  object_class->get_property = &_nc_galaxy_selfunc_get_property;
  object_class->dispose      = &_nc_galaxy_selfunc_dispose;
  object_class->finalize     = &_nc_galaxy_selfunc_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NSHELLS,
                                   g_param_spec_uint ("nshells",
                                                      NULL,
                                                      "Galaxy survey number of redshift shells",
                                                      1, G_MAXUINT, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_SHELL_SPLINES,
	                                 g_param_spec_boxed ("shell-splines",
	                                                     NULL,
	                                                     "Galaxy survey shell splines",
	                                                     NCM_TYPE_OBJ_ARRAY,
	                                                     G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_galaxy_selfunc_new:
 * @nshells: number of redshift shells
 * 
 * Creates a new instance of #NcGalaxySelfunc.
 * 
 * Returns: (transfer full): a new #NcGalaxySelfunc.
 */ 
NcGalaxySelfunc *
nc_galaxy_selfunc_new (const guint nshells)
{
	NcGalaxySelfunc *gsf = g_object_new (NC_TYPE_GALAXY_SELFUNC,
	                                     "nshells", nshells,
	                                     NULL);
	return gsf;
}

/**
 * nc_galaxy_selfunc_ref:
 * @gsf: a #NcGalaxySelfunc
 * 
 * Increase the reference of @diff by one.
 * 
 * Returns: (transfer full): @gsf.
 */
NcGalaxySelfunc *
nc_galaxy_selfunc_ref (NcGalaxySelfunc *gsf)
{
  return g_object_ref (gsf);
}

/**
 * nc_galaxy_selfunc_free:
 * @gsf: a #NcGalaxySelfunc
 * 
 * Decrease the reference count of @diff by one.
 * 
 */
void
nc_galaxy_selfunc_free (NcGalaxySelfunc *gsf)
{
  g_object_unref (gsf);
}

/**
 * nc_galaxy_selfunc_clear:
 * @gsf: a #NcGalaxySelfunc
 * 
 * Decrease the reference count of @diff by one, and sets the pointer *@diff to
 * NULL.
 * 
 */
void
nc_galaxy_selfunc_clear (NcGalaxySelfunc **gsf)
{
  g_clear_object (gsf);
}

/**
 * nc_galaxy_selfunc_set_nshells:
 * @gsf: a #NcGalaxySelfunc
 * @nshells: number of shells
 * 
 * Sets the number of shells in @gsf. This also delete all 
 * data inside @gsf if @nshells is different from the current
 * value.
 * 
 */
void
nc_galaxy_selfunc_set_nshells (NcGalaxySelfunc *gsf, const guint nshells)
{
	NcGalaxySelfuncPrivate * const self = gsf->priv;
		
  if (self->nshells != nshells)
  {
    ncm_vector_clear (&self->norm);
    ncm_vector_clear (&self->zmin);
    ncm_vector_clear (&self->zmean);
    ncm_vector_clear (&self->zmax);

    g_ptr_array_set_size ((GPtrArray *) self->dNdz_a, 0);

    self->nshells = nshells;

    if (self->nshells > 0)
    {
      guint iz;

      self->norm  = ncm_vector_new (self->nshells);
      self->zmin  = ncm_vector_new (self->nshells);
      self->zmean = ncm_vector_new (self->nshells);
      self->zmax  = ncm_vector_new (self->nshells);

      for (iz = 0; iz < self->nshells; iz++)
			{
				NcmSpline *s = ncm_spline_cubic_notaknot_new();
				ncm_obj_array_add (self->dNdz_a, G_OBJECT (s));
				ncm_spline_free (s);
			}
    }
    else
    {
      g_error ("nc_galaxy_selfunc_set_nshells: number of shells must be greater than 0.");
    }
  }
}

static void _nc_galaxy_selfunc_set_zmean (NcGalaxySelfunc *gsf, const guint shell);
static void _nc_galaxy_selfunc_set_minmax_norma (NcGalaxySelfunc *gsf, const guint shell);
static void _nc_galaxy_selfunc_set_survey_redshift_limits (NcGalaxySelfunc *gsf);

/**
 * nc_galaxy_selfunc_set_shell_splines:
 * @gsf: a #NcGalaxySelfunc
 * @dNdz_a: a #NcmObjArray
 * 
 * Sets splines for each shell using the splines in @dNdz_a.
 * 
 */
void 
nc_galaxy_selfunc_set_shell_splines (NcGalaxySelfunc *gsf, NcmObjArray *dNdz_a)
{
	NcGalaxySelfuncPrivate * const self = gsf->priv;
  guint iz;

  for (iz = 0; iz < self->nshells; iz++)
  {
		NcmSpline *dNdz = NCM_SPLINE (ncm_obj_array_peek (dNdz_a, iz));
    ncm_obj_array_set (self->dNdz_a, iz, G_OBJECT (dNdz));
		ncm_spline_prepare (dNdz);
		
		_nc_galaxy_selfunc_set_minmax_norma  (gsf, iz);
    _nc_galaxy_selfunc_set_zmean (gsf, iz);
  }

  _nc_galaxy_selfunc_set_survey_redshift_limits (gsf);
}

/**
 * nc_galaxy_selfunc_get_shell_splines:
 * @gsf: a #NcGalaxySelfunc
 * 
 * Gets the shell splines.
 * 
 * Returns: the current shell splines array.
 */
NcmObjArray *
nc_galaxy_selfunc_get_shell_splines (NcGalaxySelfunc *gsf)
{
	NcGalaxySelfuncPrivate * const self = gsf->priv;
	return ncm_obj_array_ref (self->dNdz_a);
}

/**
 * nc_galaxy_selfunc_get_nshells:
 * @gsf: a #NcGalaxySelfunc
 * 
 * Gets the current number of shells in @gsf.
 * 
 * Returns: the current number of shells.
 */
guint
nc_galaxy_selfunc_get_nshells (NcGalaxySelfunc *gsf)
{
	NcGalaxySelfuncPrivate * const self = gsf->priv;
  return self->nshells;
}

static void
_nc_galaxy_selfunc_set_zmean (NcGalaxySelfunc *gsf, const guint shell)
{
	NcGalaxySelfuncPrivate * const self = gsf->priv;
	
  guint iz;
  gdouble zmin, zmean, zmax;

  NcmSpline *dNdz = NCM_SPLINE (ncm_obj_array_peek (self->dNdz_a, shell));
  guint len       = ncm_spline_get_len (dNdz);
  NcmVector *z    = ncm_spline_get_xv (dNdz);
  NcmVector *phi  = ncm_spline_get_yv (dNdz);
  NcmVector *zxs  = ncm_vector_new (len);
  NcmSpline *spl  = ncm_spline_cubic_notaknot_new ();

  for (iz = 0; iz < len; iz++)
  {
    const gdouble m = ncm_vector_get (z, iz) * ncm_vector_get (phi, iz);
    ncm_vector_set (zxs, iz, m);
  }

  ncm_spline_set (spl, z, zxs, TRUE);
  ncm_spline_prepare (spl);

  ncm_spline_get_bounds (spl, &zmin, &zmax);
  zmean = ncm_spline_eval_integ (spl, zmin, zmax);
  ncm_vector_set (self->zmean, shell, zmean);

  ncm_vector_free (z);
  ncm_vector_free (phi);
  ncm_vector_free (zxs);
  ncm_spline_free (spl);
}

static void
_nc_galaxy_selfunc_set_minmax_norma (NcGalaxySelfunc *gsf, const guint shell)
{
	NcGalaxySelfuncPrivate * const self = gsf->priv;
	NcmSpline *dNdz = NCM_SPLINE (ncm_obj_array_peek (self->dNdz_a, shell));
	NcmVector *z    = ncm_spline_get_xv (dNdz);
	NcmVector *phi  = ncm_spline_get_yv (dNdz);	
  gdouble norm, zmin, zmax;

	zmin = ncm_vector_get (z, 0);
	zmax = ncm_vector_get (z, ncm_vector_len (z) - 1);
		
  ncm_vector_set (self->zmin, shell, zmin);
  ncm_vector_set (self->zmax, shell, zmax);

  norm = ncm_spline_eval_integ (dNdz, zmin, zmax);
  ncm_vector_set (self->norm, shell, norm);
  ncm_vector_scale (phi, 1.0 / norm);
  ncm_spline_prepare (dNdz);

  ncm_vector_free (z);
  ncm_vector_free (phi);
}

static void
_nc_galaxy_selfunc_read_file_to_spline (const gchar *filename, NcmSpline *dNdz_s)
{
  GArray *z_a    = g_array_new (FALSE, FALSE, sizeof(gdouble));
	GArray *dNdz_a = g_array_new (FALSE, FALSE, sizeof(gdouble));
  guint nl = 0;
  gdouble x, y;
  FILE *selfile;

	if (!g_file_test (filename, G_FILE_TEST_EXISTS))
    g_error ("_nc_galaxy_selfunc_read_file_to_garray: selection function file \"%s\" does not exist in folder.\n", filename);

  selfile = fopen (filename, "r");

  nl = 0;
  while (fscanf (selfile, "%le %le", &x, &y ) != EOF)
  {
    g_array_append_val (z_a, x);
    g_array_append_val (dNdz_a, y);
    nl++;
  }
  fclose (selfile);

	ncm_spline_set_array (dNdz_s, z_a, dNdz_a, TRUE);

	g_array_unref (z_a);
	g_array_unref (dNdz_a);
}

static void
_nc_galaxy_selfunc_set_dNdz_from_file (NcGalaxySelfunc *gsf, const gchar *prefix, const gchar *suffix, const guint shell)
{
	NcGalaxySelfuncPrivate * const self = gsf->priv;
  NcmSpline *dNdz = NCM_SPLINE (ncm_obj_array_peek (self->dNdz_a, shell));
	gchar *filename;

	g_assert_cmpuint (shell, <, self->nshells);
	filename = (suffix == NULL) ? g_strdup_printf ("%s.%u", prefix, shell) : g_strdup_printf ("%s.%u.%s", prefix, shell, suffix);

  _nc_galaxy_selfunc_read_file_to_spline (filename, dNdz);
	g_free (filename);
}

static void
_nc_galaxy_selfunc_set_survey_redshift_limits (NcGalaxySelfunc *gsf)
{
	NcGalaxySelfuncPrivate * const self = gsf->priv;
	
  self->zmin_all = ncm_vector_get_min (self->zmin);
  self->zmax_all = ncm_vector_get_max (self->zmax);
}

/**
 * nc_galaxy_selfunc_load_from_txts:
 * @gsf: a #NcGalaxySelfunc
 * @prefix: files prefix
 * @suffix: (nullable): files prefix
 * 
 * Reads each shell selection function from files.
 * The files are named as @prefix.n.@suffix or 
 * @prefix.n if @suffix is NULL. The value of n
 * goes from 0 to nshells - 1. Each file must contain
 * pairs of doubles describing the redshift its
 * corresponding selection function value.
 * 
 */
void
nc_galaxy_selfunc_load_from_txts (NcGalaxySelfunc *gsf, const gchar *prefix, const gchar *suffix)
{
	NcGalaxySelfuncPrivate * const self = gsf->priv;
  guint iz;

  for (iz=0; iz < self->nshells; iz++)
  {
    _nc_galaxy_selfunc_set_dNdz_from_file (gsf, prefix, suffix, iz);
		_nc_galaxy_selfunc_set_minmax_norma  (gsf, iz);
    _nc_galaxy_selfunc_set_zmean (gsf, iz);
  }

  _nc_galaxy_selfunc_set_survey_redshift_limits (gsf);
}

/**
 * nc_galaxy_selfunc_eval:
 * @gsf: a #NcGalaxySelfunc
 * @shell: shell number
 * @z: redshift $z$
 * 
 * Evaluates the @shell-th shell at redshift $z=$@z.
 * 
 * Returns: the selection function at @shell and @z. 
 */
gdouble
nc_galaxy_selfunc_eval (NcGalaxySelfunc *gsf, const guint shell, const gdouble z)
{
	NcGalaxySelfuncPrivate * const self = gsf->priv;
  return ncm_spline_eval (NCM_SPLINE (ncm_obj_array_peek (self->dNdz_a, shell)), z);
}

/**
 * nc_galaxy_selfunc_get_zmin:
 * @gsf: a #NcGalaxySelfunc
 * @shell: shell number
 * 
 * Gets the @shell-th shell lowest redshift.
 * 
 * Returns: the lowest redshift at @shell. 
 */
gdouble
nc_galaxy_selfunc_get_zmin (NcGalaxySelfunc *gsf, guint shell)
{
	NcGalaxySelfuncPrivate * const self = gsf->priv;
	return ncm_vector_get (self->zmin, shell);
}

/**
 * nc_galaxy_selfunc_get_zmean:
 * @gsf: a #NcGalaxySelfunc
 * @shell: shell number
 * 
 * Gets the @shell-th shell mean redshift.
 * 
 * Returns: the mean redshift at @shell. 
 */
gdouble
nc_galaxy_selfunc_get_zmean (NcGalaxySelfunc *gsf, guint shell)
{
	NcGalaxySelfuncPrivate * const self = gsf->priv;
	return ncm_vector_get (self->zmean, shell);
}

/**
 * nc_galaxy_selfunc_get_zmax:
 * @gsf: a #NcGalaxySelfunc
 * @shell: shell number
 * 
 * Gets the @shell-th shell maximum redshift.
 * 
 * Returns: the maximum redshift at @shell. 
 */
gdouble
nc_galaxy_selfunc_get_zmax (NcGalaxySelfunc *gsf, guint shell)
{
	NcGalaxySelfuncPrivate * const self = gsf->priv;
	return ncm_vector_get (self->zmax, shell);
}
