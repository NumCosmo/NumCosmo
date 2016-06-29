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
 * @title: NcmFftlog
 * @short_description: Abstract class for implementing logarithm fast fourier transform.
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fftlog.h"
#include "math/ncm_cfg.h"
#include "math/ncm_spline_cubic_notaknot.h"

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
  PROP_LNR0,
  PROP_LNK0,
  PROP_LR,
  PROP_N,
  PROP_PAD,
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
  fftlog->pad_p = 0.0;
  fftlog->Nr    = 0;
  fftlog->N     = 0;
  fftlog->N_2   = 0;
  fftlog->Nf    = 0;
  fftlog->Nf_2  = 0;
  fftlog->pad   = 0;
  fftlog->prepared  = FALSE;
  fftlog->evaluated = FALSE;

  fftlog->lnr_vec = NULL;
  fftlog->Gr_vec  = NULL;
  fftlog->Gr_s    = NULL;

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

#define ncm_fftlog_array_pos(fftlog,array)  

static void
_ncm_fftlog_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFftlog *fftlog = NCM_FFTLOG (object);
  g_return_if_fail (NCM_IS_FFTLOG (object));

  switch (prop_id)
  {
    case PROP_LNR0:
      ncm_fftlog_set_lnr0 (fftlog, g_value_get_double (value));
      break;
    case PROP_LNK0:
      ncm_fftlog_set_lnk0 (fftlog, g_value_get_double (value));
      break;
    case PROP_LR:
      ncm_fftlog_set_length (fftlog, g_value_get_double (value));
      break;
    case PROP_N:
      ncm_fftlog_set_size (fftlog, g_value_get_uint (value));
      break;
    case PROP_PAD:
      ncm_fftlog_set_padding (fftlog, g_value_get_double (value));
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
    case PROP_LNR0:
      g_value_set_double (value, ncm_fftlog_get_lnr0 (fftlog));
      break;
    case PROP_LNK0:
      g_value_set_double (value, ncm_fftlog_get_lnk0 (fftlog));
      break;
    case PROP_LR:
      g_value_set_double (value, ncm_fftlog_get_length (fftlog));
      break;
    case PROP_N:
      g_value_set_uint (value, ncm_fftlog_get_size (fftlog));
      break;
    case PROP_PAD:
      g_value_set_double (value, ncm_fftlog_get_padding (fftlog));
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
    fftlog->Gr_s      = g_new0 (NcmSpline *, ncomp);

#ifdef NUMCOSMO_HAVE_FFTW3
    fftlog->Gr        = g_new0 (fftw_complex *, ncomp);
    fftlog->Ym        = g_new0 (fftw_complex *, ncomp);
    fftlog->CmYm      = g_new0 (fftw_complex *, ncomp);
    fftlog->p_CmYm2Gr = g_new0 (fftw_plan, ncomp);
#else
    g_error ("Cannot construct FFTLog object, fftw3 library is not present, recompile numcosmo with fftw3 support.");
#endif /* NUMCOSMO_HAVE_FFTW3 */
  }
}

#ifdef NUMCOSMO_HAVE_FFTW3
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
    ncm_spline_clear (&fftlog->Gr_s[i]);
  }
}
#endif /* NUMCOSMO_HAVE_FFTW3 */

static void
_ncm_fftlog_finalize (GObject *object)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  NcmFftlog *fftlog = NCM_FFTLOG (object);
  _ncm_fftlog_free_all (fftlog);

  g_free (fftlog->Gr_vec);
  g_free (fftlog->Gr_s);
  
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
  object_class->finalize     = &_ncm_fftlog_finalize;

  g_object_class_install_property (object_class,
                                   PROP_LNR0,
                                   g_param_spec_double ("lnr0",
                                                        NULL,
                                                        "Center value for ln(r)",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LNK0,
                                   g_param_spec_double ("lnk0",
                                                        NULL,
                                                        "Center value for ln(k)",
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
                                   PROP_PAD,
                                   g_param_spec_double ("padding",
                                                      NULL,
                                                      "Padding percentage",
                                                      0, G_MAXDOUBLE, 1.0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
 * Increases the reference count of @fftlog by one.
 * 
 * Returns: (transfer full): @fftlog
 */
NcmFftlog *ncm_fftlog_ref (NcmFftlog *fftlog)
{
  return g_object_ref (fftlog);
}

/**
 * ncm_fftlog_free:
 * @fftlog: a #NcmFftlog
 * 
 * Decreases the reference count of @fftlog by one.
 * 
 */
void 
ncm_fftlog_free (NcmFftlog *fftlog)
{
  g_object_unref (fftlog);
}

/**
 * ncm_fftlog_clear:
 * @fftlog: a #NcmFftlog
 * 
 * If @fftlog is different from NULL, decreases the reference count of 
 * @fftlog by one and sets @fftlog to NULL.
 * 
 */
void 
ncm_fftlog_clear (NcmFftlog **fftlog)
{
  g_clear_object (fftlog);
}

/**
 * ncm_fftlog_peek_name:
 * @fftlog: a #NcmFftlog
 * 
 * 
 * Returns: (transfer none): The internal string describing #NcmFftlog.
 */
const gchar *
ncm_fftlog_peek_name (NcmFftlog *fftlog)
{
  return NCM_FFTLOG_GET_CLASS (fftlog)->name;
}

/**
 * ncm_fftlog_set_lnr0:
 * @fftlog: a #NcmFftlog
 * @lnr0: output center $\ln(r_0)$
 * 
 * Sets the center of the transform output $\ln(r_0)$.
 * 
 */
void 
ncm_fftlog_set_lnr0 (NcmFftlog *fftlog, const gdouble lnr0)
{
  if (lnr0 != fftlog->lnr0)
  {
    fftlog->lnr0      = lnr0;
    fftlog->prepared  = FALSE;
    fftlog->evaluated = FALSE;
  }
}

/**
 * ncm_fftlog_get_lnr0:
 * @fftlog: a #NcmFftlog
 * 
 * Gets the center of the transform output.
 * 
 * Returns: the output center $\ln(r_0)$.
 */
gdouble 
ncm_fftlog_get_lnr0 (NcmFftlog *fftlog)
{
  return fftlog->lnr0;
}

/**
 * ncm_fftlog_set_lnk0:
 * @fftlog: a #NcmFftlog
 * @lnk0: input center $\ln(k_0)$
 * 
 * Sets the center of the transform input $\ln(k_0)$.
 * 
 */
void 
ncm_fftlog_set_lnk0 (NcmFftlog *fftlog, const gdouble lnk0)
{
  if (lnk0 != fftlog->lnk0)
  {
    fftlog->lnk0      = lnk0;
    fftlog->prepared  = FALSE;
    fftlog->evaluated = FALSE;
  }
}

/**
 * ncm_fftlog_get_k0:
 * @fftlog: a #NcmFftlog
 * 
 * Gets the center of the transform input $\ln(k_0)$.
 * 
 * Returns: the input center $\ln(k_0)$.
 */
gdouble 
ncm_fftlog_get_lnk0 (NcmFftlog *fftlog)
{
  return fftlog->lnk0;
}

static gulong
_ncm_fftlog_fact_size (gulong n)
{
  if (n == 1)
    return 0;
  else
  {
    gulong r3 = n % 3;
    gulong r5 = n % 5;
    gulong r7 = n % 7;
    gulong m = 1;

    if (r3 == 0)
      m *= 3;
    if (r5 == 0)
      m *= 5;
    if (r7 == 0)
      m *= 7;

    if (m != 1)
    {
      if (n / m == 1)
        return m;
      else
        return m * _ncm_fftlog_fact_size (n / m);
    }
    else
    {
      const gulong d3 = 3 - r3;
      const gulong d5 = 5 - r5;
      const gulong d7 = 7 - r7;

      return _ncm_fftlog_fact_size (n + GSL_MIN (d3, GSL_MIN (d5, d7)));
    }
  }
}


/**
 * ncm_fftlog_set_size:
 * @fftlog: a #NcmFftlog
 * @n: number of knots
 * 
 * Sets the number of knots where the integrated function is evaluated.
 * 
 */
void
ncm_fftlog_set_size (NcmFftlog *fftlog, guint n)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  guint nt = n * (1.0 + fftlog->pad_p);

  fftlog->Nr = n;

  nt = _ncm_fftlog_fact_size (nt);
  fftlog->pad = nt * fftlog->pad_p * 0.5 / (1.0 + fftlog->pad_p);
  n = nt - 2 * fftlog->pad;
  
  if ((n != fftlog->N) || (n + 2 * fftlog->pad != fftlog->Nf))
  {
    guint i;
    
    fftlog->N    = n;
    fftlog->Lk_N = fftlog->Lk / (1.0 * fftlog->N - 1.0);
    fftlog->N_2  = (fftlog->N - 1) / 2;

    fftlog->Nf   = fftlog->N + 2 * fftlog->pad;
    fftlog->Nf_2 = fftlog->N_2 + fftlog->pad;

    _ncm_fftlog_free_all (fftlog);

    fftlog->Fk = fftw_alloc_complex (fftlog->Nf);
    fftlog->Cm = fftw_alloc_complex (fftlog->Nf);

    fftlog->lnr_vec = ncm_vector_new (fftlog->N);

    ncm_cfg_load_fftw_wisdom (NCM_FFTLOG_GET_CLASS (fftlog)->name);
    
    fftlog->p_Fk2Cm = fftw_plan_dft_1d (fftlog->Nf, fftlog->Fk, fftlog->Cm, FFTW_FORWARD, fftw_default_flags | FFTW_DESTROY_INPUT);

    for (i = 0; i < NCM_FFTLOG_GET_CLASS (fftlog)->ncomp; i++)
    {
      fftlog->Gr_vec[i] = ncm_vector_new (fftlog->N);
      fftlog->Gr_s[i]   = ncm_spline_cubic_notaknot_new ();
      
      fftlog->CmYm[i]   = fftw_alloc_complex (fftlog->Nf);
      fftlog->Gr[i]     = fftw_alloc_complex (fftlog->Nf);
      fftlog->Ym[i]     = fftw_alloc_complex (fftlog->Nf);
      
      fftlog->p_CmYm2Gr[i] = fftw_plan_dft_1d (fftlog->Nf, fftlog->CmYm[i], fftlog->Gr[i], FFTW_FORWARD, fftw_default_flags | FFTW_DESTROY_INPUT);
    }

    ncm_cfg_save_fftw_wisdom (NCM_FFTLOG_GET_CLASS (fftlog)->name);

    fftlog->prepared  = FALSE;
    fftlog->evaluated = FALSE;
  }
#endif /* NUMCOSMO_HAVE_FFTW3 */
}

/**
 * ncm_fftlog_set_padding:
 * @fftlog: a #NcmFftlog
 * @pad_p: padding percentage
 * 
 * Sets the size of the padding in percetange of the interval.
 * 
 */
void
ncm_fftlog_set_padding (NcmFftlog *fftlog, gdouble pad_p)
{
  if (fftlog->pad_p != pad_p)
  {
    fftlog->pad_p = pad_p;

    if (fftlog->Nr != 0)
      ncm_fftlog_set_size (fftlog, fftlog->Nr);

    fftlog->prepared  = FALSE;
    fftlog->evaluated = FALSE;
  }
}

/**
 * ncm_fftlog_get_padding:
 * @fftlog: a #NcmFftlog
 * 
 * Gets the padding percentage.
 * 
 * Returns: the padding percentage.
 */
gdouble
ncm_fftlog_get_padding (NcmFftlog *fftlog)
{
  return fftlog->pad_p;
}

/**
 * ncm_fftlog_set_length:
 * @fftlog: a #NcmFftlog
 * @Lk: period in the logarithmic space
 * 
 * Sets the length of the period @Lk, where the function is periodic in logarithmic space log10 (r).  
 * 
 */
void
ncm_fftlog_set_length (NcmFftlog *fftlog, gdouble Lk)
{
  if (fftlog->Lk != Lk)
  {
    fftlog->Lk        = Lk;
    fftlog->Lk_N      = fftlog->Lk / (1.0 * fftlog->N);
    fftlog->prepared  = FALSE;
    fftlog->evaluated = FALSE;
  }
}

#ifdef NUMCOSMO_HAVE_FFTW3
static void
_ncm_fftlog_eval (NcmFftlog *fftlog)
{
  guint ncomp = NCM_FFTLOG_GET_CLASS (fftlog)->ncomp;
  guint comp;
  gint i;

  fftw_execute (fftlog->p_Fk2Cm);
  
  if (!fftlog->prepared)
  {
    NCM_FFTLOG_GET_CLASS (fftlog)->get_Ym (fftlog);
    fftlog->prepared = TRUE;
  }
  
  for (comp = 0; comp < ncomp; comp++)
  {
    for (i = 0; i < fftlog->Nf; i++)
    {
      fftlog->CmYm[comp][i] = fftlog->Cm[i] * fftlog->Ym[comp][i]; 
    }

    fftw_execute (fftlog->p_CmYm2Gr[comp]);
  }

  NCM_FFTLOG_GET_CLASS (fftlog)->generate_Gr (fftlog);

  for (comp = 0; comp < ncomp; comp++)
  {
    ncm_spline_set (fftlog->Gr_s[comp], fftlog->lnr_vec, fftlog->Gr_vec[comp], FALSE);
  }

  fftlog->evaluated = TRUE;
}
#endif /* NUMCOSMO_HAVE_FFTW3 */

/**
 * ncm_fftlog_eval_by_vector:
 * @fftlog: a #NcmFftlog
 * @Fk: Fk function vector
 * 
 * FIXME
 * 
 */
void 
ncm_fftlog_eval_by_vector (NcmFftlog *fftlog, NcmVector *Fk)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  gint i;

  memset (fftlog->Fk, 0, sizeof (complex double) * fftlog->Nf);
  
  for (i = 0; i < fftlog->N; i++)
  {
    fftlog->Fk[fftlog->pad + i] = ncm_vector_get (Fk, i);
  }

  _ncm_fftlog_eval (fftlog);
#endif /* NUMCOSMO_HAVE_FFTW3 */
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
#ifdef NUMCOSMO_HAVE_FFTW3
  gint i;
  
  memset (fftlog->Fk, 0, sizeof (complex double) * (fftlog->Nf));

  for (i = 0; i < fftlog->N; i++)
  {
    const gint phys_i   = i - fftlog->N_2;
    const gdouble lnk_i = fftlog->lnk0 + fftlog->Lk_N * phys_i;
    const gdouble k_i   = exp (lnk_i);
    const gdouble Fk_i  = GSL_FN_EVAL (Fk, k_i);

    fftlog->Fk[fftlog->pad + i] = Fk_i;
  }
  
  _ncm_fftlog_eval (fftlog);
#endif /* NUMCOSMO_HAVE_FFTW3 */
}

/**
 * ncm_fftlog_prepare_splines:
 * @fftlog: a #NcmFftlog
 * 
 * FIXME
 * 
 */
void 
ncm_fftlog_prepare_splines (NcmFftlog *fftlog)
{
  guint ncomp = NCM_FFTLOG_GET_CLASS (fftlog)->ncomp;
  guint comp;

  g_assert (fftlog->evaluated);
  
  for (comp = 0; comp < ncomp; comp++)
    ncm_spline_prepare (fftlog->Gr_s[comp]);
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

/**
 * ncm_fftlog_peek_spline_Gr:
 * @fftlog: a #NcmFftlog
 * @comp: component number
 * 
 * FIXME
 * 
 * Returns: (transfer none): FIXME
 */
NcmSpline *
ncm_fftlog_peek_spline_Gr (NcmFftlog *fftlog, guint comp)
{
  g_assert (fftlog->evaluated);
  return fftlog->Gr_s[comp];
}

/**
 * ncm_fftlog_eval_output:
 * @fftlog: a #NcmFftlog
 * @comp: component number
 * @lnr: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
gdouble 
ncm_fftlog_eval_output (NcmFftlog *fftlog, guint comp, const gdouble lnr)
{
  return ncm_spline_eval (ncm_fftlog_peek_spline_Gr (fftlog, comp), lnr);
}

/**
 * ncm_fftlog_calibrate_size: (skip)
 * @fftlog: a #NcmFftlog
 * @Fk: FIXME
 * @reltol: FIXME
 * 
 * FIXME
 * 
 */
void 
ncm_fftlog_calibrate_size (NcmFftlog *fftlog, gsl_function *Fk, gdouble reltol)
{
  const guint ncomp = NCM_FFTLOG_GET_CLASS (fftlog)->ncomp;
  NcmSpline **s = g_new0 (NcmSpline *, ncomp);
  guint comp, size;
  gdouble lreltol = 0.0;

  ncm_fftlog_eval_by_function (fftlog, Fk);
  for (comp = 0; comp < ncomp; comp++)
  {
    s[comp] = fftlog->Gr_s[comp];
    fftlog->Gr_s[comp] = ncm_spline_copy_empty (fftlog->Gr_s[comp]);
  }

  /*printf ("# Initial size %u.\n", fftlog->N);*/
  ncm_fftlog_set_size (fftlog, fftlog->N * 1.2);
  /*printf ("# Trying size %u.\n", fftlog->N);*/
  ncm_fftlog_eval_by_function (fftlog, Fk);

  size = ncm_spline_get_len (fftlog->Gr_s[0]);
  for (comp = 0; comp < ncomp; comp++)
  {
    guint i;
    gdouble absmin, absmax;
    ncm_vector_get_absminmax (fftlog->Gr_vec[comp], &absmin, &absmax);
    
    /*printf ("# Testing component %u [% 20.15g, % 20.15g].\n", comp, absmin, absmax);*/
    for (i = 0; i < size; i++)
    {
      const gdouble lnr_i = ncm_vector_get (fftlog->lnr_vec, i);
      const gdouble lnG_i = ncm_vector_get (fftlog->Gr_vec[comp], i);
      const gdouble lnS_i = ncm_spline_eval (s[comp], lnr_i);
      const gdouble lreltol_i = fabs ((lnG_i - lnS_i) / (fabs (lnG_i) + absmax));
      
      /*printf ("% 20.15g % 20.15g % 20.15g | % 20.15e\n", lnr_i, lnG_i, lnS_i, lreltol_i);*/
      lreltol = GSL_MAX (lreltol_i, lreltol);
    }

    ncm_spline_clear (&s[comp]);
    /*printf ("# Largest error up to component %u is %e.\n", comp, lreltol);*/
  }

  g_clear_pointer (&s, g_free);

  if (lreltol > reltol)
    ncm_fftlog_calibrate_size (fftlog, Fk, reltol);
}

/**
 * ncm_fftlog_get_size:
 * @fftlog: a #NcmFftlog
 * 
 * Gets the number of knots N where the integrated function is evaluated.
 * 
 * Returns: the number of knots N.
 */
/**
 * ncm_fftlog_get_full_size:
 * @fftlog: a #NcmFftlog
 * 
 * Gets the number of knots N where the integrated function is evaluated
 * plus padding.
 * 
 * Returns: the total number of knots Nf.
 */
/**
 * ncm_fftlog_get_norma:
 * @fftlog: a #NcmFftlog
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
/**
 * ncm_fftlog_get_length:
 * @fftlog: a #NcmFftlog
 * 
 * Gets the value of the period, where the function is periodic in logarithmic space $\ln(r)$.
 * 
 * Returns: the period
 */
/**
 * ncm_fftlog_get_full_length:
 * @fftlog: a #NcmFftlog
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
/**
 * ncm_fftlog_get_mode_index:
 * @fftlog: a #NcmFftlog
 * @i: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
/**
 * ncm_fftlog_get_variable_index:
 * @fftlog: a #NcmFftlog
 * @i: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
/**
 * ncm_fftlog_get_output_index:
 * @fftlog: a #NcmFftlog
 * @i: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
/**
 * ncm_fftlog_peek_output_vector:
 * @fftlog: a #NcmFftlog
 * @comp: FIXME
 * 
 * FIXME
 * 
 * Returns: (transfer none): FIXME
 */
