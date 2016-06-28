/***************************************************************************
 *            ncm_powspec_filter.c
 *
 *  Fri June 17 10:12:06 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/* excerpt from: */
/***************************************************************************
 *            nc_window_gaussian.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/***************************************************************************
 *            nc_window_tophat.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * ncm_powspec_filter.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_powspec_filter
 * @title: NcmPowspecFilter
 * @short_description: Abstract class for implementing powerspectrum filters.
 * 
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_powspec_filter.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline2d_bicubic.h"
#include "ncm_enum_types.h"

enum
{
  PROP_0,
  PROP_TYPE,
  PROP_LNR0,
  PROP_RELTOL,
  PROP_POWERSPECTRUM
};

G_DEFINE_TYPE (NcmPowspecFilter, ncm_powspec_filter, G_TYPE_OBJECT);

static void
ncm_powspec_filter_init (NcmPowspecFilter *psf)
{
  psf->ps         = NULL;
  psf->lnr0       = 0.0;
  psf->reltol     = 0.0;
  psf->type       = NCM_POWSPEC_FILTER_TYPE_LEN;
  psf->fftlog     = NULL;
  psf->calibrated = FALSE;
  psf->var        = ncm_spline2d_bicubic_notaknot_new ();
  psf->dvar       = ncm_spline2d_bicubic_notaknot_new ();
}

static void
ncm_powspec_filter_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPowspecFilter *psf = NCM_POWSPEC_FILTER (object);
  g_return_if_fail (NCM_IS_POWSPEC_FILTER (object));

  switch (prop_id)
  {
    case PROP_TYPE:
      ncm_powspec_filter_set_type (psf, g_value_get_enum (value));
      break;
    case PROP_LNR0:
      ncm_powspec_filter_set_lnr0 (psf, g_value_get_double (value));
      break;
    case PROP_RELTOL:
      psf->reltol = g_value_get_double (value);
      break;
    case PROP_POWERSPECTRUM:
      psf->ps = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_powspec_filter_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPowspecFilter *psf = NCM_POWSPEC_FILTER (object);
  g_return_if_fail (NCM_IS_POWSPEC_FILTER (object));

  switch (prop_id)
  {
    case PROP_TYPE:
      g_value_set_enum (value, psf->type);
      break;
    case PROP_LNR0:
      g_value_set_double (value, psf->lnr0);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, psf->reltol);
      break;
    case PROP_POWERSPECTRUM:
      g_value_set_object (value, psf->ps);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_powspec_filter_dispose (GObject *object)
{
  NcmPowspecFilter *psf = NCM_POWSPEC_FILTER (object);

  ncm_powspec_clear (&psf->ps);
  ncm_fftlog_clear (&psf->fftlog);

  ncm_spline2d_clear (&psf->var);
  ncm_spline2d_clear (&psf->dvar);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_filter_parent_class)->dispose (object);
}

static void
ncm_powspec_filter_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_filter_parent_class)->finalize (object);
}

static void
ncm_powspec_filter_class_init (NcmPowspecFilterClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = ncm_powspec_filter_set_property;
  object_class->get_property = ncm_powspec_filter_get_property;
  object_class->dispose      = ncm_powspec_filter_dispose;
  object_class->finalize     = ncm_powspec_filter_finalize;

  g_object_class_install_property (object_class,
                                   PROP_LNR0,
                                   g_param_spec_double ("lnr0",
                                                        NULL,
                                                        "Output center value",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance for calibration",
                                                        0, G_MAXDOUBLE, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TYPE,
                                   g_param_spec_enum ("type",
                                                      NULL,
                                                      "Filter type",
                                                      NCM_TYPE_POWSPEC_FILTER_TYPE, NCM_POWSPEC_FILTER_TYPE_TOPHAT,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_POWERSPECTRUM,
                                   g_param_spec_object ("powerspectrum",
                                                        NULL,
                                                        "NcmPowspec object",
                                                        NCM_TYPE_POWSPEC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_powspec_filter_new:
 * @ps: a #NcmPowspec
 * @type: a type from #NcmPowspecFilterType.
 * 
 * Creates a new #NcPowspecFilter from the powerspectrum @ps. 
 * 
 * Returns: (transfer full): the newly created #NcPowspecFilter.
 */
NcmPowspecFilter *
ncm_powspec_filter_new (NcmPowspec *ps, NcmPowspecFilterType type)
{
  NcmPowspecFilter *psf = g_object_new (NCM_TYPE_POWSPEC_FILTER,
                                        "type", type,
                                        "powerspectrum", ps,
                                        NULL);
  
  return psf;
}

void
ncm_powspec_filter_set_type (NcmPowspecFilter *psf, NcmPowspecFilterType type)
{
  if (type != psf->type)
  {
    const gdouble lnk_min = log (psf->ps->kmin);
    const gdouble lnk_max = log (psf->ps->kmax);
    const gdouble lnk0    = 0.5 * (lnk_max + lnk_min);
    const gdouble Lk0     = (lnk_max - lnk_min);

    ncm_fftlog_clear (&psf->fftlog);
    psf->type = type;
    psf->lnr0 = -lnk0;

    switch (psf->type)
    {
      case NCM_POWSPEC_FILTER_TYPE_TOPHAT:
        psf->fftlog = NCM_FFTLOG (ncm_fftlog_tophatwin2_new (psf->lnr0, lnk0, Lk0, 200));
        break;
      case NCM_POWSPEC_FILTER_TYPE_GAUSS:
        psf->fftlog = NCM_FFTLOG (ncm_fftlog_gausswin2_new (psf->lnr0, lnk0, Lk0, 200));
        break;
      default:
        g_assert_not_reached ();
        break;
    }

    ncm_fftlog_set_padding (psf->fftlog, 2.0);
  }
}

typedef struct _NcmPowspecFilterArg
{
  NcmPowspecFilter *psf;
  NcmModel *model;
  gdouble z;
} NcmPowspecFilterArg;

static gdouble 
_ncm_powspec_filter_k2Pk (gdouble k, gpointer userdata)
{
  NcmPowspecFilterArg *arg = (NcmPowspecFilterArg *) userdata;
  const gdouble k2 = k * k;
  const gdouble Pk = ncm_powspec_eval (arg->psf->ps, arg->model, arg->z, k);
  const gdouble f  = Pk * k2 / (2.0 * M_PI * M_PI);
  
  return f;
}

/**
 * ncm_powspec_filter_prepare:
 * @psf: a #NcmPowspecFilter
 * @model: a #NcmModel
 * 
 * Prepares the object applying the filter to the powerspectrum.
 * 
 */
void
ncm_powspec_filter_prepare (NcmPowspecFilter *psf, NcmModel *model)
{
  NcmPowspecFilterArg arg;
  gsl_function F;
  const gdouble zmin = psf->ps->zi;
  const gdouble zmax = psf->ps->zf;

  F.function = &_ncm_powspec_filter_k2Pk;
  F.params   = &arg;

  arg.psf   = psf;
  arg.model = model;
  arg.z     = 0.0;

  ncm_powspec_prepare_if_needed (psf->ps, model);

/*printf ("# zmin % 20.15g zmax % 20.15g\n", zmin, zmax);*/
  
  if (!psf->calibrated)
  {
    NcmMatrix *lnvar, *dlnvar;
    NcmVector *z_vec, *lnr_vec;
    guint N_k = 0, N_z = 0;
    guint i;

    ncm_powspec_get_nknots (psf->ps, &N_z, &N_k);
    ncm_fftlog_calibrate_size (psf->fftlog, &F, psf->reltol);
    N_k = ncm_fftlog_get_size (psf->fftlog);

    g_assert_cmpuint (N_z, >, 0);
    g_assert_cmpuint (N_k, >, 0);
    
    lnvar   = ncm_matrix_new (N_z, N_k);
    dlnvar  = ncm_matrix_new (N_z, N_k);
    z_vec   = ncm_vector_new (N_z);
    lnr_vec = ncm_fftlog_get_vector_lnr (psf->fftlog);

    for (i = 0; i < N_z; i++)
    {
      NcmVector *var_z  = ncm_matrix_get_row (lnvar, i);
      NcmVector *dvar_z = ncm_matrix_get_row (dlnvar, i);

      arg.z = zmin + i * (zmax - zmin) / (N_z - 1.0);
      ncm_fftlog_eval_by_function (psf->fftlog, &F);

      ncm_vector_set (z_vec, i, arg.z);
      ncm_vector_memcpy (var_z, ncm_fftlog_peek_output_vector (psf->fftlog, 0));
      ncm_vector_memcpy (dvar_z, ncm_fftlog_peek_output_vector (psf->fftlog, 1));
    }

    ncm_spline2d_set (psf->var, lnr_vec, z_vec, lnvar, TRUE);
    ncm_spline2d_set (psf->dvar, lnr_vec, z_vec, dlnvar, TRUE);

    ncm_vector_free (z_vec);
    ncm_vector_free (lnr_vec);
    ncm_matrix_free (lnvar);
    ncm_matrix_free (dlnvar);
    
    psf->calibrated = TRUE;
  }
}

/**
 * ncm_powspec_filter_set_lnr0:
 * @psf: a #NcmPowspecFilter
 * @lnr0: the output center value $\ln(r_0)$
 * 
 * FIXME
 * 
 */
void 
ncm_powspec_filter_set_lnr0 (NcmPowspecFilter *psf, gdouble lnr0)
{
  psf->lnr0 = lnr0;
}

/**
 * ncm_powspec_filter_get_r_min:
 * @psf: a #NcmPowspecFilter
 * 
 * FIXME
 * 
 * Returns: FIXME 
 */
gdouble
ncm_powspec_filter_get_r_min (NcmPowspecFilter *psf)
{
  return exp (psf->lnr0 - ncm_fftlog_get_length (psf->fftlog) * 0.5);
}

/**
 * ncm_powspec_filter_get_r_max:
 * @psf: a #NcmPowspecFilter
 * 
 * FIXME
 * 
 * Returns: FIXME 
 */
gdouble
ncm_powspec_filter_get_r_max (NcmPowspecFilter *psf)
{
  return exp (psf->lnr0 + ncm_fftlog_get_length (psf->fftlog) * 0.5);
}

/**
 * ncm_powspec_filter_eval_lnvar_lnr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @lnr: FIXME
 * 
 * Evaluate the filtered variance at @lnr.
 * 
 * Returns: FIXME 
 */
gdouble
ncm_powspec_filter_eval_lnvar_lnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return ncm_spline2d_eval (psf->var, lnr, z);
}

/**
 * ncm_powspec_filter_eval_var:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @r: FIXME
 * 
 * Evaluate the filtered variance at @r.
 * 
 * Returns: FIXME 
 */
gdouble
ncm_powspec_filter_eval_var (NcmPowspecFilter *psf, const gdouble z, const gdouble r)
{
  return exp (ncm_powspec_filter_eval_lnvar_lnr (psf, z, log (r)));
}

/**
 * ncm_powspec_filter_eval_sigma:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @r: FIXME
 * 
 * Evaluate the filtered variance at @r.
 * 
 * Returns: FIXME 
 */
gdouble
ncm_powspec_filter_eval_sigma (NcmPowspecFilter *psf, const gdouble z, const gdouble r)
{
  return exp (0.5 * ncm_powspec_filter_eval_lnvar_lnr (psf, z, log (r)));
}

/**
 * ncm_powspec_filter_eval_dlnvar_dlnr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @lnr: FIXME
 * 
 * Evaluate the filtered variance at @lnr.
 * 
 * Returns: FIXME 
 */
gdouble
ncm_powspec_filter_eval_dlnvar_dlnr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return ncm_spline2d_eval (psf->dvar, lnr, z);
}

/**
 * ncm_powspec_filter_eval_dlnvar_dr:
 * @psf: a #NcmPowspecFilter
 * @z: redshift $z$
 * @lnr: FIXME
 * 
 * Evaluate the filtered variance at @lnr.
 * 
 * Returns: FIXME 
 */
gdouble
ncm_powspec_filter_eval_dlnvar_dr (NcmPowspecFilter *psf, const gdouble z, const gdouble lnr)
{
  return ncm_powspec_filter_eval_dlnvar_dlnr (psf, z, lnr) * exp (-lnr);
}

/**
 * ncm_powspec_filter_volume_rm3:
 * @psf: a #NcmPowspecFilter
 * 
 * Calculates the volume of the filter over $r^3$.
 * 
 * Returns: Filter's volume over the radius squared $V r^{-3}$. 
 */
gdouble 
ncm_powspec_filter_volume_rm3 (NcmPowspecFilter *psf)
{
  const gdouble tophat_volumeRm3 = 4.0 * M_PI / 3.0;
  const gdouble gauss_volumeRm3  = sqrt (2.0 * M_PI) * sqrt(2.0 * M_PI) * sqrt(2.0 * M_PI);
  switch (psf->type)
  {
    case NCM_POWSPEC_FILTER_TYPE_TOPHAT:
      return tophat_volumeRm3;
      break;
    case NCM_POWSPEC_FILTER_TYPE_GAUSS:
      return gauss_volumeRm3;
      break;
    default:
      g_assert_not_reached ();
      break;
  }
  return 0.0;
}

