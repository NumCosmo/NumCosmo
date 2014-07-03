/***************************************************************************
 *            ncm_fftlog.c
 *
 *  Fri May 18 16:44:23 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:ncm_fftlog
 * @title: Logarithm Fast Fourier Algorithm
 * @short_description: Object implementing logarithm fast fourier transform
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fftlog.h"
#include "math/ncm_cfg.h"

#include <math.h>
#ifdef NUMCOSMO_HAVE_FFTW3 
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */
#ifndef HAVE_FFTW3_ALLOC
#define fftw_alloc_real(n) (double *) fftw_malloc(sizeof(double) * (n))
#define fftw_alloc_complex(n) (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (n))
#endif /* HAVE_FFTW3_ALLOC */

enum
{
  PROP_0,
  PROP_R0,
  PROP_K0,
  PROP_LR,
  PROP_N,
  PROP_NCOMP,
  PROP_NAME,
};

G_DEFINE_ABSTRACT_TYPE (NcmFftlog, ncm_fftlog, G_TYPE_OBJECT);

static void
ncm_fftlog_init (NcmFftlog *fftlog)
{
  fftlog->lnr0  = 0.0;
  fftlog->lnk0  = 0.0;
  fftlog->Lk    = 0.0;
  fftlog->Lk_N  = 0.0;
  fftlog->N     = 0;
  fftlog->prepared  = FALSE;
  fftlog->evaluated = FALSE;

  fftlog->lnr_vec = NULL;
  fftlog->Gr_vec  = NULL;

#ifdef NUMCOSMO_HAVE_FFTW3
  fftlog->Fk        = NULL;
  fftlog->Cm        = NULL;
  fftlog->Gr        = NULL;
  fftlog->Ym        = NULL;
  fftlog->CmYm      = NULL;
  fftlog->p_Fk2Cm   = NULL;
  fftlog->p_CmYm2Gr = NULL;
#endif /* NUMCOSMO_HAVE_FFTW3 */
}

static void
_ncm_fftlog_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFftlog *fftlog = NCM_FFTLOG (object);
  g_return_if_fail (NCM_IS_FFTLOG (object));

  switch (prop_id)
  {
    case PROP_R0:
      fftlog->lnr0 = log (g_value_get_double (value));
      break;
    case PROP_K0:
      fftlog->lnk0 = log (g_value_get_double (value));
      break;
    case PROP_LR:
      fftlog->Lk = g_value_get_double (value);
      break;
    case PROP_N:
      ncm_fftlog_set_size (fftlog, g_value_get_uint (value));
      break;
    case PROP_NCOMP:
      g_assert_not_reached ();
      break;
    case PROP_NAME:
      g_assert_not_reached ();
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fftlog_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFftlog *fftlog = NCM_FFTLOG (object);

  g_return_if_fail (NCM_IS_FFTLOG (object));

  switch (prop_id)
  {
    case PROP_R0:
      g_value_set_double (value, exp (fftlog->lnr0));
      break;
    case PROP_K0:
      g_value_set_double (value, exp (fftlog->lnk0));
      break;
    case PROP_LR:
      g_value_set_double (value, fftlog->Lk);
      break;
    case PROP_N:
      g_value_set_uint (value, fftlog->N);
      break;
    case PROP_NAME:
      g_value_set_string (value, NCM_FFTLOG_GET_CLASS (fftlog)->name);
      break;
    case PROP_NCOMP:
      g_value_set_uint (value, NCM_FFTLOG_GET_CLASS (fftlog)->ncomp);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fftlog_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fftlog_parent_class)->constructed (object);
  {
    NcmFftlog *fftlog = NCM_FFTLOG (object);
    guint ncomp = NCM_FFTLOG_GET_CLASS (fftlog)->ncomp;

    g_assert_cmpuint (ncomp, >, 0);

    fftlog->Gr_vec    = g_new0 (NcmVector *, ncomp);

    fftlog->Gr        = g_new0 (fftw_complex *, ncomp);
    fftlog->Ym        = g_new0 (fftw_complex *, ncomp);
    fftlog->CmYm      = g_new0 (fftw_complex *, ncomp);
    fftlog->p_CmYm2Gr = g_new0 (fftw_plan, ncomp);
#ifndef NUMCOSMO_HAVE_FFTW3
    g_error ("Cannot construct FFTLog object, fftw3 library is not present, recompile numcosmo with fftw3 support.");
#endif /* NUMCOSMO_HAVE_FFTW3 */
  }
}

static void
_ncm_fftlog_free_all (NcmFftlog *fftlog)
{
  guint ncomp = NCM_FFTLOG_GET_CLASS (fftlog)->ncomp;
  guint i;

  g_clear_pointer (&fftlog->Fk, fftw_free);
  g_clear_pointer (&fftlog->Cm, fftw_free);
  g_clear_pointer (&fftlog->p_Fk2Cm, fftw_destroy_plan);

  ncm_vector_clear (&fftlog->lnr_vec);
  
  for (i = 0; i < ncomp; i++)
  {
    g_clear_pointer (&fftlog->Gr[i], fftw_free);
    g_clear_pointer (&fftlog->Ym[i], fftw_free);
    g_clear_pointer (&fftlog->CmYm[i], fftw_free);
    
    g_clear_pointer (&fftlog->p_CmYm2Gr[i], fftw_destroy_plan);
    ncm_vector_clear (&fftlog->Gr_vec[i]);
  }
}

static void
ncm_fftlog_finalize (GObject *object)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  NcmFftlog *fftlog = NCM_FFTLOG (object);
  _ncm_fftlog_free_all (fftlog);

  g_free (fftlog->Gr_vec);
  
  g_free (fftlog->Gr);
  g_free (fftlog->Ym);
  g_free (fftlog->CmYm);
  g_free (fftlog->p_CmYm2Gr);

#endif /* NUMCOSMO_HAVE_FFTW3 */
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fftlog_parent_class)->finalize (object);
}

static void
ncm_fftlog_class_init (NcmFftlogClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_fftlog_constructed;
  object_class->set_property = &_ncm_fftlog_set_property;
  object_class->get_property = &_ncm_fftlog_get_property;
  object_class->finalize     = &ncm_fftlog_finalize;

  g_object_class_install_property (object_class,
                                   PROP_R0,
                                   g_param_spec_double ("r0",
                                                        NULL,
                                                        "Center value for r",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_K0,
                                   g_param_spec_double ("k0",
                                                        NULL,
                                                        "Center value for k",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_LR,
                                   g_param_spec_double ("Lk",
                                                        NULL,
                                                        "Function log-period",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_N,
                                   g_param_spec_uint ("N",
                                                      NULL,
                                                      "Number of knots",
                                                      0, G_MAXUINT, 10,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_NCOMP,
                                   g_param_spec_string ("ncomp",
                                                      NULL,
                                                      "Number of components",
                                                      0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_NAME,
                                   g_param_spec_string ("name",
                                                      NULL,
                                                      "FFTW Plan wisdown name",
                                                      "fftlog_default_wisdown",
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

}


/**
 * ncm_fftlog_ref:
 * @fftlog: a #NcmFftlog
 * 
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
NcmFftlog *ncm_fftlog_ref (NcmFftlog *fftlog)
{
  return g_object_ref (fftlog);
}

/**
 * ncm_fftlog_free:
 * @fftlog: a #NcmFfftlog
 * 
 * FIXME
 * 
 */
void 
ncm_fftlog_free (NcmFftlog *fftlog)
{
  g_object_unref (fftlog);
}

/**
 * ncm_fftlog_clear:
 * @fftlog: a #NcmFfftlog
 * 
 * FIXME
 * 
 */
void 
ncm_fftlog_clear (NcmFftlog **fftlog)
{
  g_clear_object (fftlog);
}

/**
 * ncm_fftlog_peek_name:
 * @fftlog: FIXME
 * 
 * FIXME
 * 
 * Returns: (transfer none): FIXME
 */
gchar *
ncm_fftlog_peek_name (NcmFftlog *fftlog)
{
  return NCM_FFTLOG_GET_CLASS (fftlog)->name;
}

/**
 * ncm_fftlog_set_size:
 * @fftlog: FIXME
 * @n: FIXME
 * 
 * FIXME
 * 
 */
void
ncm_fftlog_set_size (NcmFftlog *fftlog, guint n)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  if (n != fftlog->N)
  {
    guint i;
    
    fftlog->N = n;
    if ((fftlog->N % 2) == 0)
      fftlog->N++;
    
    fftlog->Lk_N = fftlog->Lk / (1.0 * fftlog->N);
    fftlog->N_2  = (fftlog->N - 1) / 2;

    _ncm_fftlog_free_all (fftlog);

    fftlog->Fk = fftw_alloc_complex (fftlog->N);
    fftlog->Cm = fftw_alloc_complex (fftlog->N);

    fftlog->lnr_vec = ncm_vector_new (fftlog->N);

    ncm_cfg_load_fftw_wisdom (NCM_FFTLOG_GET_CLASS (fftlog)->name);
    
    fftlog->p_Fk2Cm = fftw_plan_dft_1d (fftlog->N, fftlog->Fk, fftlog->Cm, FFTW_FORWARD, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);

    for (i = 0; i < NCM_FFTLOG_GET_CLASS (fftlog)->ncomp; i++)
    {
      fftlog->Gr_vec[i] = ncm_vector_new (fftlog->N);
      
      fftlog->CmYm[i]   = fftw_alloc_complex (fftlog->N);
      fftlog->Gr[i]     = fftw_alloc_complex (fftlog->N);
      fftlog->Ym[i]     = fftw_alloc_complex (fftlog->N);
      
      fftlog->p_CmYm2Gr[i] = fftw_plan_dft_1d (fftlog->N, fftlog->CmYm[i], fftlog->Gr[i], FFTW_FORWARD, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
    }

    ncm_cfg_save_fftw_wisdom (NCM_FFTLOG_GET_CLASS (fftlog)->name);

    fftlog->prepared  = FALSE;
    fftlog->evaluated = FALSE;
  }
#endif /* NUMCOSMO_HAVE_FFTW3 */
}

/**
 * ncm_fftlog_get_size:
 * @fftlog: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
guint
ncm_fftlog_get_size (NcmFftlog *fftlog)
{
  return fftlog->N;
}

/**
 * ncm_fftlog_set_length:
 * @fftlog: FIXME
 * @Lk: FIXME
 * 
 * FIXME
 * 
 */
void
ncm_fftlog_set_length (NcmFftlog *fftlog, gdouble Lk)
{
  if (fftlog->Lk != Lk)
  {
    fftlog->Lk = Lk;
    fftlog->Lk_N = fftlog->Lk / (1.0 * fftlog->N);
    fftlog->prepared  = FALSE;
    fftlog->evaluated = FALSE;
  }
}

/**
 * ncm_fftlog_get_length:
 * @fftlog: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
gdouble
ncm_fftlog_get_length (NcmFftlog *fftlog)
{
  return fftlog->Lk;
}

void
_ncm_fftlog_eval (NcmFftlog *fftlog)
{
  guint comp;
  gint ii;
  guint ncomp = NCM_FFTLOG_GET_CLASS (fftlog)->ncomp;

  fftw_execute (fftlog->p_Fk2Cm);

  if (!fftlog->prepared)
  {
    NCM_FFTLOG_GET_CLASS (fftlog)->get_Ym (fftlog);
    fftlog->prepared = TRUE;
  }
  
  for (comp = 0; comp < ncomp; comp++)
  {
    for (ii = 0; ii < fftlog->N; ii++)
    {
      fftlog->CmYm[comp][ii] = fftlog->Cm[ii] * fftlog->Ym[comp][ii]; 
    }

    fftw_execute (fftlog->p_CmYm2Gr[comp]);
  }

  NCM_FFTLOG_GET_CLASS (fftlog)->generate_Gr (fftlog);

  fftlog->evaluated = TRUE;
}

/**
 * ncm_fftlog_eval_by_vector:
 * @fftlog: a #NcmFftlog
 * @Fk: Fk function vector.
 * 
 * FIXME
 * 
 */
void 
ncm_fftlog_eval_by_vector (NcmFftlog *fftlog, NcmVector *Fk)
{
  gint i;

  for (i = -fftlog->N_2; i <= fftlog->N_2; i++)
  {
    const gint ii = (i < 0) ? i + fftlog->N : i;
    fftlog->Fk[ii] = ncm_vector_get (Fk, ii);
  }
  _ncm_fftlog_eval (fftlog);
}

/**
 * ncm_fftlog_eval_by_function: (skip)
 * @fftlog: a #NcmFftlog
 * @Fk: Fk function pointer
 * 
 * FIXME
 * 
 */
void 
ncm_fftlog_eval_by_function (NcmFftlog *fftlog, gsl_function *Fk)
{
  gint i;
  for (i = -fftlog->N_2; i <= fftlog->N_2; i++)
  {
    const gint ii = (i < 0) ? i + fftlog->N : i;
    const gdouble k = exp (fftlog->lnk0 + fftlog->Lk_N * i);
    const gdouble Fk_i = GSL_FN_EVAL (Fk, k);
    fftlog->Fk[ii] = Fk_i;
  }
  _ncm_fftlog_eval (fftlog);
}

/**
 * ncm_fftlog_get_vector_lnr:
 * @fftlog: a #NcmFftlog
 * 
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
NcmVector *
ncm_fftlog_get_vector_lnr (NcmFftlog *fftlog)
{
  return ncm_vector_ref (fftlog->lnr_vec);
}

/**
 * ncm_fftlog_get_vector_Gr:
 * @fftlog: a #NcmFftlog
 * @comp: component number
 * 
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
NcmVector *
ncm_fftlog_get_vector_Gr (NcmFftlog *fftlog, guint comp)
{
  g_assert (fftlog->evaluated);
  return ncm_vector_ref (fftlog->Gr_vec[comp]);
}
