/***************************************************************************
 *            ncm_fftlog.c
 *
 *  Fri May 18 16:44:23 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
 * NcmFftlog:
 *
 * Base class for FFTLog transforms.
 *
 * This class computes the Fast Fourier Transform of a function assumed to be a
 * periodic sequence of logarithmically spaced points, following the FFTLog
 * approach of Hamilton (2000) with extensions. A function $G(r)$ is decomposed as
 * \begin{equation*}
 * G(r) = \int_0^\infty F(k)\,K(kr)\,\mathrm{d}k,
 * \end{equation*}
 * with $F(k)$ expanded in its lowest Fourier modes over a logarithmic
 * fundamental interval and $K(kr)$ a kernel. The kernel-dependent coefficients
 * are provided by the concrete implementations (e.g. #NcmFftlogTophatwin2,
 * #NcmFftlogGausswin2).
 *
 * For the full derivation, discretization, and padding scheme, see the
 * theoretical background page:
 * <a href="../../theory/ncm/fftlog/fftlog.html">FFTLog</a>.
 * Reference: [Hamilton (2000)](https://arxiv.org/abs/astro-ph/9905191).
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fftlog/ncm_fftlog.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_util.h"
#include "ncm/spline/ncm_spline_cubic_notaknot.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <complex.h>
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

#define fftw_alloc_real(n) (double *) fftw_malloc (sizeof (double) * (n))
#define fftw_alloc_complex(n) (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * (n))

typedef struct _NcmFftlogPrivate
{
  gint Nr;
  gint N;
  gint N_2;
  gint Nf;
  gint Nf_2;
  guint max_n;
  guint nderivs;
  guint pad;
  gdouble lnk0;
  gdouble lnr0;
  gdouble eval_r_min;
  gdouble eval_r_max;
  gdouble Lk;
  gdouble Lk_N;
  gdouble pad_p;
  gboolean smooth_padding;
  gboolean use_eval_int;
  gboolean noring;
  gboolean prepared;
  gboolean evaluated;
  NcmVector *lnr_vec;
  GPtrArray *Gr_vec;
  GPtrArray *Gr_s;
  GPtrArray *Ym;

  fftw_complex *Fk;
  fftw_complex *Cm;
  fftw_complex *Gr;
  fftw_complex *CmYm;
  fftw_plan p_Fk2Cm;
  fftw_plan p_CmYm2Gr;
} NcmFftlogPrivate;

enum
{
  PROP_0,
  PROP_NDERIV,
  PROP_LNR0,
  PROP_LNK0,
  PROP_LR,
  PROP_N,
  PROP_MAX_N,
  PROP_PAD,
  PROP_NORING,
  PROP_NAME,
  PROP_USE_EVAL_INT,
  PROP_SMOOTH_PADDING,
  PROP_EVAL_R_MIN,
  PROP_EVAL_R_MAX,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmFftlog, ncm_fftlog, G_TYPE_OBJECT)

static void
ncm_fftlog_init (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  self->lnr0           = 0.0;
  self->use_eval_int   = FALSE;
  self->smooth_padding = FALSE;
  self->eval_r_min     = 0.0;
  self->eval_r_max     = 0.0;
  self->lnk0           = 0.0;
  self->Lk             = 0.0;
  self->Lk_N           = 0.0;
  self->pad_p          = 0.0;
  self->Nr             = 0;
  self->N              = 0;
  self->N_2            = 0;
  self->Nf             = 0;
  self->Nf_2           = 0;
  self->max_n          = 0;
  self->pad            = 0;
  self->noring         = FALSE;
  self->prepared       = FALSE;
  self->evaluated      = FALSE;

  self->lnr_vec = NULL;
  self->Gr_vec  = g_ptr_array_new ();
  self->Gr_s    = g_ptr_array_new ();

  g_ptr_array_set_free_func (self->Gr_vec, (GDestroyNotify) ncm_vector_free);
  g_ptr_array_set_free_func (self->Gr_s, (GDestroyNotify) ncm_spline_free);

  self->Fk        = NULL;
  self->Cm        = NULL;
  self->Gr        = NULL;
  self->Ym        = g_ptr_array_new ();
  self->CmYm      = NULL;
  self->p_Fk2Cm   = NULL;
  self->p_CmYm2Gr = NULL;
  g_ptr_array_set_free_func (self->Ym, (GDestroyNotify) fftw_free);
}

#define ncm_fftlog_array_pos(fftlog, array)

static void
_ncm_fftlog_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFftlog *fftlog = NCM_FFTLOG (object);

  g_return_if_fail (NCM_IS_FFTLOG (object));

  switch (prop_id)
  {
    case PROP_NDERIV:
      ncm_fftlog_set_nderivs (fftlog, g_value_get_uint (value));
      break;
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
    case PROP_MAX_N:
      ncm_fftlog_set_max_size (fftlog, g_value_get_uint (value));
      break;
    case PROP_PAD:
      ncm_fftlog_set_padding (fftlog, g_value_get_double (value));
      break;
    case PROP_NORING:
      ncm_fftlog_set_noring (fftlog, g_value_get_boolean (value));
      break;
    case PROP_NAME:
      g_assert_not_reached ();
      break;
    case PROP_USE_EVAL_INT:
      ncm_fftlog_use_eval_interval (fftlog, g_value_get_boolean (value));
      break;
    case PROP_SMOOTH_PADDING:
      ncm_fftlog_use_smooth_padding (fftlog, g_value_get_boolean (value));
      break;
    case PROP_EVAL_R_MIN:
      ncm_fftlog_set_eval_r_min (fftlog, g_value_get_double (value));
      break;
    case PROP_EVAL_R_MAX:
      ncm_fftlog_set_eval_r_max (fftlog, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_fftlog_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFftlog *fftlog             = NCM_FFTLOG (object);
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  g_return_if_fail (NCM_IS_FFTLOG (object));

  switch (prop_id)
  {
    case PROP_NDERIV:
      g_value_set_uint (value, ncm_fftlog_get_nderivs (fftlog));
      break;
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
    case PROP_MAX_N:
      g_value_set_uint (value, ncm_fftlog_get_max_size (fftlog));
      break;
    case PROP_PAD:
      g_value_set_double (value, ncm_fftlog_get_padding (fftlog));
      break;
    case PROP_NORING:
      g_value_set_boolean (value, ncm_fftlog_get_noring (fftlog));
      break;
    case PROP_NAME:
      g_value_set_string (value, NCM_FFTLOG_GET_CLASS (fftlog)->name);
      break;
    case PROP_USE_EVAL_INT:
      g_value_set_boolean (value, self->use_eval_int);
      break;
    case PROP_SMOOTH_PADDING:
      g_value_set_boolean (value, self->smooth_padding);
      break;
    case PROP_EVAL_R_MIN:
      g_value_set_double (value, ncm_fftlog_get_eval_r_min (fftlog));
      break;
    case PROP_EVAL_R_MAX:
      g_value_set_double (value, ncm_fftlog_get_eval_r_max (fftlog));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_fftlog_free_all (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  g_clear_pointer (&self->Fk, fftw_free);
  g_clear_pointer (&self->Cm, fftw_free);
  g_clear_pointer (&self->CmYm, fftw_free);
  g_clear_pointer (&self->Gr, fftw_free);

  g_clear_pointer (&self->p_Fk2Cm, ncm_cfg_fftw_plan_destroy);
  g_clear_pointer (&self->p_CmYm2Gr, ncm_cfg_fftw_plan_destroy);

  ncm_vector_clear (&self->lnr_vec);

  g_ptr_array_set_size (self->Gr_vec, 0);
  g_ptr_array_set_size (self->Gr_s, 0);
  g_ptr_array_set_size (self->Ym, 0);
}

static void
_ncm_fftlog_finalize (GObject *object)
{
  NcmFftlog *fftlog             = NCM_FFTLOG (object);
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  _ncm_fftlog_free_all (fftlog);

  g_clear_pointer (&self->Gr_vec, g_ptr_array_unref);
  g_clear_pointer (&self->Gr_s,   g_ptr_array_unref);
  g_clear_pointer (&self->Ym,     g_ptr_array_unref);


  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fftlog_parent_class)->finalize (object);
}

static void
ncm_fftlog_class_init (NcmFftlogClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_fftlog_set_property;
  object_class->get_property = &_ncm_fftlog_get_property;
  object_class->finalize     = &_ncm_fftlog_finalize;

  /**
   * NcmFftlog:nderivs:
   *
   * The number of derivatives to be estimated.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_NDERIV,
                                   g_param_spec_uint ("nderivs",
                                                      NULL,
                                                      "Number of derivatives",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFftlog:lnr0:
   *
   * Center $\ln r_0$ of the output grid.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LNR0,
                                   g_param_spec_double ("lnr0",
                                                        NULL,
                                                        "Center value for ln(r)",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFftlog:lnk0:
   *
   * Center $\ln k_0$ of the input grid.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LNK0,
                                   g_param_spec_double ("lnk0",
                                                        NULL,
                                                        "Center value for ln(k)",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFftlog:Lk:
   *
   * Length $L > 0$ of the fundamental interval in $\ln k$, the period of $F$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LR,
                                   g_param_spec_double ("Lk",
                                                        NULL,
                                                        "Function log-period",
                                                        G_MINDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFftlog:N:
   *
   * The number of knots in the fundamental interval.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_N,
                                   g_param_spec_uint ("N",
                                                      NULL,
                                                      "Number of knots",
                                                      0, G_MAXUINT, 10,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFftlog:max-n:
   *
   * The maximum number of knots in the fundamental interval.
   * This limit is used when calibrating the number of knots.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_N,
                                   g_param_spec_uint ("max-n",
                                                      NULL,
                                                      "Maximum number of knots",
                                                      0, G_MAXUINT, 100000,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFftlog:padding:
   *
   * Padding as a fraction of the number of knots $N$: the transform uses
   * $N_f = N(1 + \mathrm{padding})$ points, the extra ones split between the two ends.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_PAD,
                                   g_param_spec_double ("padding",
                                                        NULL,
                                                        "Padding fraction",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFftlog:no-ringing:
   *
   * True to use the no-ringing adjustment of $\ln(r_0)$ and False otherwise.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_NORING,
                                   g_param_spec_boolean ("no-ringing",
                                                         NULL,
                                                         "No ringing",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFftlog:name:
   *
   * FFTW Plan wisdom's name to perform the transformation.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_NAME,
                                   g_param_spec_string ("name",
                                                        NULL,
                                                        "FFTW Plan wisdom name",
                                                        "fftlog_default_wisdom",
                                                        G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFftlog:use-eval-int:
   *
   * Whether to use evaluation interval
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_USE_EVAL_INT,
                                   g_param_spec_boolean ("use-eval-int",
                                                         NULL,
                                                         "Whether to use evaluation interval",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFftlog:use-smooth-padding:
   *
   * Whether the padding continues the input smoothly, see
   * ncm_fftlog_use_smooth_padding().
   */
  g_object_class_install_property (object_class,
                                   PROP_SMOOTH_PADDING,
                                   g_param_spec_boolean ("use-smooth-padding",
                                                         NULL,
                                                         "Whether the padding continues the input smoothly",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFftlog:eval-r-min:
   *
   * The minimum value of the evaluation interval.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_EVAL_R_MIN,
                                   g_param_spec_double ("eval-r-min",
                                                        NULL,
                                                        "Evaluation r_min",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFftlog:eval-r-max:
   *
   * The maximum value of the evaluation interval.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_EVAL_R_MAX,
                                   g_param_spec_double ("eval-r-max",
                                                        NULL,
                                                        "Evaluation r_max",
                                                        0.0, G_MAXDOUBLE, G_MAXDOUBLE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_fftlog_ref:
 * @fftlog: a #NcmFftlog
 *
 * Increases the reference count of @fftlog by one.
 *
 * Returns: (transfer full): @fftlog
 */
NcmFftlog *
ncm_fftlog_ref (NcmFftlog *fftlog)
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
 * If *@fftlog is not NULL, decreases its reference count by one and sets *@fftlog to NULL.
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
 * Returns: (transfer none): the #NcmFftlog:name, used for the FFTW wisdom
 */
const gchar *
ncm_fftlog_peek_name (NcmFftlog *fftlog)
{
  return NCM_FFTLOG_GET_CLASS (fftlog)->name;
}

/**
 * ncm_fftlog_reset:
 * @fftlog: a #NcmFftlog
 *
 * Reset the evaluation and internal coefficients forcing
 * their recomputation.
 *
 */
void
ncm_fftlog_reset (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  self->prepared  = FALSE;
  self->evaluated = FALSE;
}

/**
 * ncm_fftlog_set_nderivs:
 * @fftlog: a #NcmFftlog
 * @nderivs: the number of derivatives
 *
 * Sets @nderivs as the number of derivatives to calculate.
 *
 */
void
ncm_fftlog_set_nderivs (NcmFftlog *fftlog, guint nderivs)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  if (self->nderivs != nderivs)
  {
    if (nderivs < self->nderivs)
    {
      g_ptr_array_set_size (self->Gr_vec, nderivs + 1);
      g_ptr_array_set_size (self->Gr_s, nderivs + 1);
      g_ptr_array_set_size (self->Ym, nderivs + 1);
    }
    else
    {
      if (self->N != 0)
      {
        guint i;

        for (i = self->nderivs + 1; i <= nderivs; i++)
        {
          NcmVector *Gr_vec_i = ncm_vector_new (self->N);
          NcmSpline *Gr_s_i   = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (self->lnr_vec, Gr_vec_i, FALSE));
          fftw_complex *Ym_i  = fftw_alloc_complex (self->Nf);

          g_ptr_array_add (self->Gr_vec, Gr_vec_i);
          g_ptr_array_add (self->Gr_s, Gr_s_i);
          g_ptr_array_add (self->Ym, Ym_i);
        }
      }

      ncm_fftlog_reset (fftlog);
    }

    self->nderivs = nderivs;
  }
}

/**
 * ncm_fftlog_get_nderivs:
 * @fftlog: a #NcmFftlog
 *
 * Gets the number of derivatives the object is currently
 * calculating.
 *
 * Returns: the number of derivatives calculated.
 */
guint
ncm_fftlog_get_nderivs (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return self->nderivs;
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
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  if (lnr0 != self->lnr0)
  {
    self->lnr0 = lnr0;
    ncm_fftlog_reset (fftlog);
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
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return self->lnr0;
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
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  if (lnk0 != self->lnk0)
  {
    self->lnk0 = lnk0;
    ncm_fftlog_reset (fftlog);
  }
}

/**
 * ncm_fftlog_get_lnk0:
 * @fftlog: a #NcmFftlog
 *
 * Gets the center of the transform input $\ln(k_0)$.
 *
 * Returns: the input center $\ln(k_0)$.
 */
gdouble
ncm_fftlog_get_lnk0 (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return self->lnk0;
}

/**
 * ncm_fftlog_set_size:
 * @fftlog: a #NcmFftlog
 * @n: number of knots $N$
 *
 * Sets the size of the transform from @n knots: the full size $N_f^\prime$ is the
 * smallest $2^a 3^b 5^c 7^d \geq N(1 + \mathrm{padding})$, and the fundamental interval
 * keeps $N^\prime = N_f^\prime - 2 N_\mathrm{pad}$ of them, the rest being padding.
 *
 */
void
ncm_fftlog_set_size (NcmFftlog *fftlog, guint n)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);
  gint nt                       = n * (1.0 + self->pad_p);
  gint n_new;

  self->Nr = n;

  nt        = ncm_util_fact_size (nt);
  self->pad = nt * self->pad_p * 0.5 / (1.0 + self->pad_p);
  n_new     = nt - 2 * self->pad;

  if ((n_new != self->N) || (n_new + 2 * (gint) self->pad != self->Nf))
  {
    guint fftw_default_flags = ncm_cfg_get_fftw_default_flag ();
    gboolean first_plan;
    gint i;

    self->N    = n_new;
    self->N_2  = self->N / 2;
    self->Lk_N = self->Lk / (1.0 * self->N);

    self->Nf   = self->N + 2 * self->pad;
    self->Nf_2 = self->N_2 + self->pad;

    _ncm_fftlog_free_all (fftlog);

    self->Fk   = fftw_alloc_complex (self->Nf);
    self->Cm   = fftw_alloc_complex (self->Nf);
    self->CmYm = fftw_alloc_complex (self->Nf);
    self->Gr   = fftw_alloc_complex (self->Nf);

    self->lnr_vec = ncm_vector_new (self->N);

    first_plan = ncm_cfg_fftw_plan_begin ("ncm_fftlog_dft_1d_%d", self->Nf);

    self->p_Fk2Cm   = fftw_plan_dft_1d (self->Nf, self->Fk,   self->Cm, FFTW_FORWARD, fftw_default_flags | FFTW_DESTROY_INPUT);
    self->p_CmYm2Gr = fftw_plan_dft_1d (self->Nf, self->CmYm, self->Gr, FFTW_FORWARD, fftw_default_flags | FFTW_DESTROY_INPUT);

    for (i = 0; i <= (gint) self->nderivs; i++)
    {
      NcmVector *Gr_vec_i = ncm_vector_new (self->N);
      NcmSpline *Gr_s_i   = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (self->lnr_vec, Gr_vec_i, FALSE));
      fftw_complex *Ym_i  = fftw_alloc_complex (self->Nf);

      g_ptr_array_add (self->Gr_vec, Gr_vec_i);
      g_ptr_array_add (self->Gr_s, Gr_s_i);
      g_ptr_array_add (self->Ym, Ym_i);
    }

    ncm_cfg_fftw_plan_end (first_plan);

    ncm_fftlog_reset (fftlog);
  }
}

/**
 * ncm_fftlog_set_max_size:
 * @fftlog: a #NcmFftlog
 * @max_n: maximum number of knots
 *
 * Sets the maximum number of knots in the fundamental interval.
 *
 */
void
ncm_fftlog_set_max_size (NcmFftlog *fftlog, guint max_n)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  self->max_n = max_n;
}

/**
 * ncm_fftlog_get_max_size:
 * @fftlog: a #NcmFftlog
 *
 * Gets the maximum number of knots in the fundamental interval.
 *
 * Returns: the #NcmFftlog:max-n
 */
guint
ncm_fftlog_get_max_size (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return self->max_n;
}

/**
 * ncm_fftlog_set_padding:
 * @fftlog: a #NcmFftlog
 * @pad_p: padding fraction
 *
 * Sets #NcmFftlog:padding.
 *
 */
void
ncm_fftlog_set_padding (NcmFftlog *fftlog, gdouble pad_p)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  if (self->pad_p != pad_p)
  {
    self->pad_p = pad_p;

    if (self->Nr != 0)
      ncm_fftlog_set_size (fftlog, self->Nr);

    ncm_fftlog_reset (fftlog);
  }
}

/**
 * ncm_fftlog_get_padding:
 * @fftlog: a #NcmFftlog
 *
 * Returns: the #NcmFftlog:padding fraction
 */
gdouble
ncm_fftlog_get_padding (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return self->pad_p;
}

/**
 * ncm_fftlog_set_noring:
 * @fftlog: a #NcmFftlog
 * @active: whether to use the no-ringing adjustment of $\ln(r_0)$
 *
 * Sets whether to use the no-ringing adjustment of $\ln(r_0)$.
 *
 */
void
ncm_fftlog_set_noring (NcmFftlog *fftlog, gboolean active)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  if ((!self->noring && active) || (self->noring && !active))
  {
    self->noring = active;
    ncm_fftlog_reset (fftlog);
  }
}

/**
 * ncm_fftlog_get_noring:
 * @fftlog: a #NcmFftlog
 *
 * Returns: whether the no-ringing adjustment of $\ln r_0$ is active
 */
gboolean
ncm_fftlog_get_noring (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return self->noring;
}

/**
 * ncm_fftlog_set_length:
 * @fftlog: a #NcmFftlog
 * @Lk: period in the logarithmic space
 *
 * Sets the length of the period @Lk, where the function is periodic in logarithmic space $\ln k$.
 *
 */
void
ncm_fftlog_set_length (NcmFftlog *fftlog, gdouble Lk)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  if (!(Lk > 0.0))
    g_error ("ncm_fftlog_set_length: the period must be positive, got %g.", Lk);

  if (self->Lk != Lk)
  {
    self->Lk   = Lk;
    self->Lk_N = self->Lk / (1.0 * self->N);
    ncm_fftlog_reset (fftlog);
  }
}

/**
 * ncm_fftlog_use_eval_interval:
 * @fftlog: a #NcmFftlog
 * @use_eval_interval: whether to restrict the output
 *
 * Sets whether to use a restricted evaluation interval $[r_\mathrm{min}, r_\mathrm{max}]$.
 * See ncm_fftlog_set_eval_r_min() and ncm_fftlog_set_eval_r_max().
 *
 */
void
ncm_fftlog_use_eval_interval (NcmFftlog *fftlog, gboolean use_eval_interval)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  if (use_eval_interval)
  {
    self->use_eval_int = TRUE;
  }
  else
  {
    if (self->use_eval_int)
    {
      guint nd;

      for (nd = 0; nd <= self->nderivs; nd++)
      {
        NcmVector *Gr_vec_nd = g_ptr_array_index (self->Gr_vec, nd);

        ncm_spline_set (g_ptr_array_index (self->Gr_s, nd), self->lnr_vec, Gr_vec_nd, FALSE);
      }
    }

    self->use_eval_int = FALSE;
  }
}

/**
 * ncm_fftlog_use_smooth_padding:
 * @fftlog: a #NcmFftlog
 * @use_smooth_padding: whether to pad smoothly
 *
 * Sets whether the padding continues the input instead of holding zeros. Zeros put a
 * step at each end of the interval, and the transform of a step rings at $r$ near the
 * inverse of that end at a level that falls only as $1/N$. The continuation is the power
 * law of each end, with the value and log-slope of $F$ there taken from a cubic through
 * the four nearest knots, anchored at the end itself so that it does not move with the
 * knots. Where $F k$, the integrand in $\ln k$, would grow along the continuation, the
 * power law is cut by a Gaussian in $\ln k$ of width the inverse of that log-slope, so
 * that the padding never holds much more integral than the interval does. The two
 * continuations are joined by a $C^\infty$ partition of unity that
 * switches over the middle fifth of the padding, so each side keeps its own continuation
 * over two fifths of the padding, the periodic input is smooth everywhere, and the
 * transform converges as $N^{-3}$ over the whole output grid. The continuation is fitted
 * in log space, so $F$ must be positive at the four knots nearest each end.
 *
 * The result within a few e-foldings of $1/k_\mathrm{max}$ or $1/k_\mathrm{min}$ depends
 * on the continuation, which is an extrapolation of the input beyond its interval. Far
 * below the peak of the output, the periodic images of the padded input set a floor of
 * order $e^{-L_T} \int F \, \mathrm{d}k$; a padding fraction of one puts it at
 * $e^{-2L}$ times that integral.
 *
 */
void
ncm_fftlog_use_smooth_padding (NcmFftlog *fftlog, gboolean use_smooth_padding)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  self->smooth_padding = use_smooth_padding;
  self->evaluated      = FALSE;
}

/**
 * ncm_fftlog_set_eval_r_min:
 * @fftlog: a #NcmFftlog
 * @eval_r_min: the value of $r_\mathrm{min}$
 *
 * Sets $r_\mathrm{min}$ to @eval_r_min.
 *
 */
void
ncm_fftlog_set_eval_r_min (NcmFftlog *fftlog, const gdouble eval_r_min)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  if (self->eval_r_min != eval_r_min)
    self->eval_r_min = eval_r_min;
}

/**
 * ncm_fftlog_set_eval_r_max:
 * @fftlog: a #NcmFftlog
 * @eval_r_max: the value of $r_\mathrm{max}$
 *
 * Sets $r_\mathrm{max}$ to @eval_r_max.
 *
 */
void
ncm_fftlog_set_eval_r_max (NcmFftlog *fftlog, const gdouble eval_r_max)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  if (self->eval_r_max != eval_r_max)
    self->eval_r_max = eval_r_max;
}

/**
 * ncm_fftlog_get_eval_r_min:
 * @fftlog: a #NcmFftlog
 *
 * Returns: the value of $r_\mathrm{min}$
 */
gdouble
ncm_fftlog_get_eval_r_min (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return self->eval_r_min;
}

/**
 * ncm_fftlog_get_eval_r_max:
 * @fftlog: a #NcmFftlog
 *
 * Returns: the value of $r_\mathrm{max}$
 */
gdouble
ncm_fftlog_get_eval_r_max (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return self->eval_r_max;
}

static void
_ncm_fftlog_eval (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);
  guint nd;
  gint i;

  fftw_execute (self->p_Fk2Cm);

  if (!self->prepared)
  {
    const gdouble Lt       = ncm_fftlog_get_full_length (fftlog);
    const gdouble twopi_Lt = 2.0 * M_PI / Lt;
    fftw_complex *Ym_0     = g_ptr_array_index (self->Ym, 0);

    /* Knots sit at (i - Nf_2) Lk_N on both sides, so the two offsets cancel only when
     * 2 Nf_2 = Nf. For an odd Nf one knot is left over and enters the phase here. */
    gdouble lnr0k0 = self->lnk0 + self->lnr0 + (self->Nf - 2 * self->Nf_2) * self->Lk_N;
    fftw_complex *Ym_ndm1;

    NCM_FFTLOG_GET_CLASS (fftlog)->compute_Ym (fftlog, Ym_0);

    if (self->noring)
    {
      gint i;

      for (i = 0; i < 5; i++)
      {
        fftw_complex YNf_2_0 = Ym_0[self->Nf / 2];
        const gdouble theta  = carg (YNf_2_0);
        const gdouble M      = (self->Nf / Lt) * lnr0k0 - theta / M_PI;
        const glong M_round  = M;
        const gdouble dM     = M - M_round;

        lnr0k0     -= (Lt / self->Nf) * dM;
        self->lnr0 -= (Lt / self->Nf) * dM;
      }
    }

    for (i = 0; i < self->Nf; i++)
    {
      const gint phys_i      = ncm_fftlog_get_mode_index (fftlog, i);
      const complex double a = twopi_Lt * phys_i * I;

      Ym_0[i] *= cexp (-a * lnr0k0);

      Ym_ndm1 = Ym_0;

      for (nd = 1; nd <= self->nderivs; nd++)
      {
        fftw_complex *Ym_nd = g_ptr_array_index (self->Ym, nd);

        Ym_nd[i] = -(1.0 + a) * Ym_ndm1[i];
        Ym_ndm1  = Ym_nd;
      }
    }

    if ((self->Nf % 2) == 0)
    {
      const gint Nf_2_index = ncm_fftlog_get_array_index (fftlog, +self->Nf / 2);

      for (nd = 0; nd <= self->nderivs; nd++)
      {
        fftw_complex *Ym_nd = g_ptr_array_index (self->Ym, nd);

        Ym_nd[Nf_2_index] = creal (Ym_nd[Nf_2_index]);
      }
    }

    self->prepared = TRUE;
  }

  for (i = 0; i < self->N; i++)
  {
    const gint phys_i = i - self->N_2;
    const gdouble lnr = self->lnr0 + phys_i * self->Lk_N;

    ncm_vector_set (self->lnr_vec, i, lnr);
  }

  for (nd = 0; nd <= self->nderivs; nd++)
  {
    const gdouble norma = ncm_fftlog_get_norma (fftlog);
    NcmVector *Gr_nd    = g_ptr_array_index (self->Gr_vec, nd);
    fftw_complex *Ym_nd = g_ptr_array_index (self->Ym, nd);

    /* For real F and a real kernel C_m Y_m is conjugate symmetric, so the transform is
     * real; the Nyquist entry of an even Nf is already real, since both factors are. No
     * entry is altered: making one member of a +-m pair real drops half of that mode's
     * sine part. */
    for (i = 0; i < self->Nf; i++)
    {
      self->CmYm[i] = self->Cm[i] * Ym_nd[i];
    }

    fftw_execute (self->p_CmYm2Gr);

    for (i = 0; i < self->N; i++)
    {
      const gdouble lnr     = ncm_vector_get (self->lnr_vec, i);
      const gdouble rm1     = exp (-lnr);
      const gdouble Gr_nd_i = creal (self->Gr[i + self->pad]) * rm1 / norma;

      ncm_vector_set (Gr_nd, i, Gr_nd_i);
    }
  }

  self->evaluated = TRUE;
}

/**
 * ncm_fftlog_get_Ym:
 * @fftlog: a #NcmFftlog
 * @size: (out): return size
 *
 * Computes the kernel coefficients $Y_m$, before the phase factor of the grid centres,
 * as interleaved real and imaginary parts of the ncm_fftlog_get_full_size() modes.
 * They are computed into the buffer the transform uses, so the next evaluation prepares
 * its own copy again.
 *
 * Returns: (transfer none) (array length=size): $Y_m$
 */
gdouble *
ncm_fftlog_get_Ym (NcmFftlog *fftlog, guint *size)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  fftw_complex *Ym_0 = g_ptr_array_index (self->Ym, 0);

  NCM_FFTLOG_GET_CLASS (fftlog)->compute_Ym (fftlog, Ym_0);
  self->prepared = FALSE;

  size[0] = ncm_fftlog_get_full_size (fftlog) * 2;

  return (gdouble *) Ym_0;
}

/**
 * ncm_fftlog_get_lnk_vector:
 * @fftlog: a #NcmFftlog
 * @lnk: output, of length ncm_fftlog_get_size()
 *
 * Fills @lnk with the knots $\ln k_m$ of the fundamental interval.
 *
 */
void
ncm_fftlog_get_lnk_vector (NcmFftlog *fftlog, NcmVector *lnk)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);
  gint i;

  g_assert_cmpuint (self->N, ==, ncm_vector_len (lnk));

  for (i = 0; i < self->N; i++)
  {
    const gint phys_i   = i - self->N_2;
    const gdouble lnk_i = self->lnk0 + self->Lk_N * phys_i;

    ncm_vector_set (lnk, i, lnk_i);
  }
}

static void _ncm_fftlog_end_power_law (const gdouble u[4], const gdouble lnF[4], gdouble *A, gdouble *s);
static gdouble _ncm_fftlog_continuation (const gdouble A, const gdouble s, const gdouble sigma, const gdouble u);
static gdouble _ncm_fftlog_smooth_step (const gdouble t);

/* Fills the padding with the continuation described in ncm_fftlog_use_smooth_padding().
 * The padding is one stretch of 2 pad slots in the periodic array, running from just
 * above ln k_max (slot pad + N) around the wrap to just below ln k_min (slot pad - 1);
 * t runs from 0 to 1 along it and the partition switches over its middle fifth. */
static void
_ncm_fftlog_add_smooth_padding (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);
  const gdouble lnk_max         = self->lnk0 + 0.5 * self->Lk;
  const gdouble lnk_min         = self->lnk0 - 0.5 * self->Lk;
  const gdouble lnk_first       = self->lnk0 - self->N_2 * self->Lk_N;
  const gdouble lnk_last        = self->lnk0 + (self->N - 1 - self->N_2) * self->Lk_N;
  const gint stretch            = 2 * self->pad;
  gdouble u_hi[4], lnF_hi[4], u_lo[4], lnF_lo[4];
  gdouble A_hi, s_hi, A_lo, s_lo;
  gint i;

  g_assert_cmpint (self->N, >=, 4);

  /* u is the distance from the end along ln k, negative at the knots, positive in the
   * padding. */
  for (i = 0; i < 4; i++)
  {
    const gdouble F_hi = creal (self->Fk[self->pad + self->N - 4 + i]);
    const gdouble F_lo = creal (self->Fk[self->pad + i]);

    if (!(F_hi > 0.0) || !(F_lo > 0.0))
      g_error ("ncm_fftlog: smooth padding needs F > 0 at the four knots nearest each end of the interval, got %g and %g.", F_lo, F_hi);

    u_hi[i]   = lnk_last - (3 - i) * self->Lk_N - lnk_max;
    lnF_hi[i] = log (F_hi);

    u_lo[i]   = lnk_min - (lnk_first + i * self->Lk_N);
    lnF_lo[i] = log (F_lo);
  }

  _ncm_fftlog_end_power_law (u_hi, lnF_hi, &A_hi, &s_hi);
  _ncm_fftlog_end_power_law (u_lo, lnF_lo, &A_lo, &s_lo);

  for (i = 0; i < stretch; i++)
  {
    const gdouble x_hi = lnk_last + (i + 1) * self->Lk_N - lnk_max;
    const gdouble x_lo = lnk_min - (lnk_first - (stretch - i) * self->Lk_N);
    const gdouble t    = (i + 1.0) / (stretch + 1.0);
    const gdouble w    = _ncm_fftlog_smooth_step ((t - 0.4) / 0.2);
    const gdouble F_i  = (1.0 - w) * _ncm_fftlog_continuation (A_hi, s_hi, s_hi + 1.0, x_hi) + w * _ncm_fftlog_continuation (A_lo, s_lo, s_lo, x_lo);

    if (i < (gint) self->pad)
      self->Fk[self->pad + self->N + i] = F_i;
    else
      self->Fk[i - self->pad] = F_i;
  }
}

/* Value A and log-slope s of ln F at u = 0 from the cubic through four (u, ln F) knots. */
static void
_ncm_fftlog_end_power_law (const gdouble u[4], const gdouble lnF[4], gdouble *A, gdouble *s)
{
  gdouble dd[4], c[4], w[4];

  gsl_poly_dd_init (dd, u, lnF, 4);
  gsl_poly_dd_taylor (c, 0.0, dd, u, 4, w);

  *A = c[0];
  *s = c[1];
}

/* The power law exp (A + s u) at distance u > 0 beyond an end, cut by a Gaussian in ln k
 * of width 1 / sigma when sigma > 0. sigma is the log-slope of F k, what the transform
 * integrates in ln k: s + 1 above the interval, s below it (there k shrinks as u grows).
 * A cut on the growth of F alone lets F k grow for several e-foldings and its periodic
 * image floods the output far below its peak. */
static gdouble
_ncm_fftlog_continuation (const gdouble A, const gdouble s, const gdouble sigma, const gdouble u)
{
  const gdouble su = sigma * u;

  return exp (A + s * u - ((sigma > 0.0) ? 0.5 * su * su : 0.0));
}

/* C-infinity step from 0 at t <= 0 to 1 at t >= 1, all derivatives zero at both ends. */
static gdouble
_ncm_fftlog_smooth_step (const gdouble t)
{
  if (t <= 0.0)
    return 0.0;

  if (t >= 1.0)
    return 1.0;

  {
    const gdouble a = exp (-1.0 / t);
    const gdouble b = exp (-1.0 / (1.0 - t));

    return a / (a + b);
  }
}

/**
 * ncm_fftlog_eval_by_vector:
 * @fftlog: a #NcmFftlog
 * @Fk: values of $F$ at the knots $\ln k_m$, see ncm_fftlog_get_lnk_vector()
 *
 * Computes the transform from the values of $F$ at the knots.
 *
 */
void
ncm_fftlog_eval_by_vector (NcmFftlog *fftlog, NcmVector *Fk)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);
  gint i;

  g_assert_cmpuint (self->N, ==, ncm_vector_len (Fk));

  memset (self->Fk, 0, sizeof (complex double) * self->Nf);

  for (i = 0; i < self->N; i++)
  {
    self->Fk[self->pad + i] = ncm_vector_get (Fk, i);
  }

  if (self->smooth_padding)
    _ncm_fftlog_add_smooth_padding (fftlog);

  _ncm_fftlog_eval (fftlog);
}

/**
 * ncm_fftlog_eval_by_gsl_function: (skip)
 * @fftlog: a #NcmFftlog
 * @Fk: Fk function pointer
 *
 * Evaluates the function @Fk at each knot $\ln k_m$.
 *
 */
void
ncm_fftlog_eval_by_gsl_function (NcmFftlog *fftlog, gsl_function *Fk)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);
  gint i;

  memset (self->Fk, 0, sizeof (complex double) * (self->Nf));

  for (i = 0; i < self->N; i++)
  {
    const gint phys_i   = i - self->N_2;
    const gdouble lnk_i = self->lnk0 + self->Lk_N * phys_i;
    const gdouble k_i   = exp (lnk_i);
    const gdouble Fk_i  = GSL_FN_EVAL (Fk, k_i);

    self->Fk[self->pad + i] = Fk_i;
  }

  if (self->smooth_padding)
    _ncm_fftlog_add_smooth_padding (fftlog);

  _ncm_fftlog_eval (fftlog);
}

/**
 * ncm_fftlog_eval_by_function:
 * @fftlog: a #NcmFftlog
 * @Fk: (scope call): a #NcmFftlogFunc
 * @user_data: @Fk user data
 *
 * Evaluates the function @Fk at each knot $\ln k_m$.
 *
 */
void
ncm_fftlog_eval_by_function (NcmFftlog *fftlog, NcmFftlogFunc Fk, gpointer user_data)
{
  gsl_function F;

  F.function = Fk;
  F.params   = user_data;

  ncm_fftlog_eval_by_gsl_function (fftlog, &F);
}

/**
 * ncm_fftlog_prepare_splines:
 * @fftlog: a #NcmFftlog
 *
 * Prepares the set of splines respective to the function $G(r)$
 * and, if required, its n-order derivatives.
 *
 */
void
ncm_fftlog_prepare_splines (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);
  guint nd;

  g_assert (self->evaluated);

  if (self->use_eval_int)
  {
    g_assert_cmpfloat (self->eval_r_min, <, self->eval_r_max);
    {
      const gint i0           = ncm_vector_find_closest_index (self->lnr_vec, log (self->eval_r_min));
      const gint i1           = ncm_vector_find_closest_index (self->lnr_vec, log (self->eval_r_max)) + 1;
      const gint size         = i1 - i0 + 1;
      NcmVector *eval_lnr_vec = ncm_vector_get_subvector (self->lnr_vec, i0, size);

      for (nd = 0; nd <= self->nderivs; nd++)
      {
        NcmVector *Gr_vec_nd      = g_ptr_array_index (self->Gr_vec, nd);
        NcmVector *eval_Gr_vec_nd = ncm_vector_get_subvector (Gr_vec_nd, i0, size);

        ncm_spline_set (g_ptr_array_index (self->Gr_s, nd), eval_lnr_vec, eval_Gr_vec_nd, TRUE);

        ncm_vector_free (eval_Gr_vec_nd);
      }

      ncm_vector_free (eval_lnr_vec);
    }
  }
  else
  {
    for (nd = 0; nd <= self->nderivs; nd++)
      ncm_spline_prepare (g_ptr_array_index (self->Gr_s, nd));
  }
}

/**
 * ncm_fftlog_get_vector_lnr:
 * @fftlog: a #NcmFftlog
 *
 * Gets the vector of the $\ln r$ knots.
 *
 * Returns: (transfer full): the $\ln r$ knots
 */
NcmVector *
ncm_fftlog_get_vector_lnr (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  if (!self->prepared)
    g_warning ("ncm_fftlog_get_vector_lnr: returning the lnr vector without preparing evaluating, the vector may change after evaluation.");

  return (self->lnr_vec != NULL) ? ncm_vector_ref (self->lnr_vec) : NULL;
}

/**
 * ncm_fftlog_get_vector_Gr:
 * @fftlog: a #NcmFftlog
 * @nderiv: derivative number
 *
 * Gets the vector of the transformed function $G(r)$, @nderiv = 0, or
 * its @nderiv-th derivative with respect to $\ln r$.
 *
 * Returns: (transfer full): a vector of $G(r)$ values or its @nderiv-th derivative.
 */
NcmVector *
ncm_fftlog_get_vector_Gr (NcmFftlog *fftlog, guint nderiv)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  g_assert (self->evaluated);

  return ncm_vector_ref (g_ptr_array_index (self->Gr_vec, nderiv));
}

/**
 * ncm_fftlog_peek_spline_Gr:
 * @fftlog: a #NcmFftlog
 * @nderiv: derivative number
 *
 * Peeks the spline of $G(r)$, @nderiv = 0,
 * or the spline of the @nderiv-th derivative of $G(r)$ with
 * respect to $\ln r$.
 *
 * Returns: (transfer none): the @nderiv component of the spline.
 */
NcmSpline *
ncm_fftlog_peek_spline_Gr (NcmFftlog *fftlog, guint nderiv)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  g_assert (self->evaluated);

  return g_ptr_array_index (self->Gr_s, nderiv);
}

/**
 * ncm_fftlog_eval_output:
 * @fftlog: a #NcmFftlog
 * @nderiv: derivative number
 * @lnr: logarithm base e of $r$
 *
 * Evaluates the function $G(r)$, or the @nderiv-th derivative,
 * at the point @lnr.
 *
 * Returns: $\mathrm{d}^n G / (\mathrm{d}\ln r)^n$ at @lnr, with $n$ = @nderiv
 */
gdouble
ncm_fftlog_eval_output (NcmFftlog *fftlog, guint nderiv, const gdouble lnr)
{
  return ncm_spline_eval (ncm_fftlog_peek_spline_Gr (fftlog, nderiv), lnr);
}

/**
 * ncm_fftlog_calibrate_size_gsl: (skip)
 * @fftlog: a #NcmFftlog
 * @Fk: Fk function pointer
 * @reltol: relative tolerance
 *
 * Increases the number of knots by 20% at a time until $G(r)$ and its derivatives change
 * by less than @reltol, relative to each component's peak, from one size to the next.
 * Each step computes one transform. Aborts if #NcmFftlog:max-n is passed first. Leaves
 * @fftlog evaluated at the final size.
 *
 */
void
ncm_fftlog_calibrate_size_gsl (NcmFftlog *fftlog, gsl_function *Fk, const gdouble reltol)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);
  NcmSpline **prev              = g_new0 (NcmSpline *, self->nderivs + 1);

  ncm_fftlog_eval_by_gsl_function (fftlog, Fk);
  ncm_fftlog_prepare_splines (fftlog);

  while (TRUE)
  {
    const gint N_prev = self->N;
    gdouble err       = 0.0;
    NcmVector *lnr;
    guint nd, i, size, n_try;

    /* Keep the current splines: set_size() replaces them when the size changes. */
    for (nd = 0; nd <= self->nderivs; nd++)
      prev[nd] = ncm_spline_ref (g_ptr_array_index (self->Gr_s, nd));

    /* Grow by 20%; the factorable full size can round back to the same N, so grow
     * further until it moves. */
    for (n_try = 1; self->N == N_prev; n_try++)
      ncm_fftlog_set_size (fftlog, (guint) ceil (N_prev * pow (1.2, n_try)));

    ncm_fftlog_eval_by_gsl_function (fftlog, Fk);
    ncm_fftlog_prepare_splines (fftlog);

    /* Largest difference between the two sizes on the new knots, relative to the value
     * plus the component's peak. */
    lnr  = ncm_spline_get_xv (g_ptr_array_index (self->Gr_s, 0));
    size = ncm_spline_get_len (g_ptr_array_index (self->Gr_s, 0));

    for (nd = 0; nd <= self->nderivs; nd++)
    {
      NcmVector *Gr = ncm_spline_get_yv (g_ptr_array_index (self->Gr_s, nd));
      gdouble absmin, absmax;

      ncm_vector_get_absminmax (Gr, &absmin, &absmax);

      for (i = 0; i < size; i++)
      {
        const gdouble G_new  = ncm_vector_get (Gr, i);
        const gdouble G_prev = ncm_spline_eval (prev[nd], ncm_vector_get (lnr, i));

        err = GSL_MAX (err, fabs (G_new - G_prev) / (fabs (G_new) + absmax));
      }

      ncm_spline_clear (&prev[nd]);
      ncm_vector_free (Gr);
    }

    ncm_vector_free (lnr);

    if (err <= reltol)
      break;

    if (self->N > (gint) self->max_n)
      g_error ("ncm_fftlog_calibrate_size_gsl: the maximum number of knots (%u) was reached "
               "at relative accuracy %e, the requested one is %e.", self->max_n, err, reltol);
  }

  g_free (prev);
}

/**
 * ncm_fftlog_calibrate_size:
 * @fftlog: a #NcmFftlog
 * @Fk: (scope call): a #NcmFftlogFunc
 * @user_data: @Fk user data
 * @reltol: relative tolerance
 *
 * Increases the number of knots by 20% at a time until $G(r)$ and its derivatives change
 * by less than @reltol, relative to each component's peak, from one size to the next.
 * Each step computes one transform. Aborts if #NcmFftlog:max-n is passed first. Leaves
 * @fftlog evaluated at the final size.
 *
 */
void
ncm_fftlog_calibrate_size (NcmFftlog *fftlog, NcmFftlogFunc Fk, gpointer user_data, const gdouble reltol)
{
  gsl_function F;

  F.function = Fk;
  F.params   = user_data;

  ncm_fftlog_calibrate_size_gsl (fftlog, &F, reltol);
}

/**
 * ncm_fftlog_get_size:
 * @fftlog: a #NcmFftlog
 *
 * Gets the number of knots $N^\prime$ where the integrated function is evaluated.
 *
 * Returns: the number of knots $N^\prime$.
 */
guint
ncm_fftlog_get_size (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return self->N;
}

/**
 * ncm_fftlog_get_full_size:
 * @fftlog: a #NcmFftlog
 *
 * Gets the number of knots $N_f^\prime$ where the integrated function is evaluated
 * plus padding.
 *
 * Returns: the total number of knots $N_f^\prime$.
 */
gint
ncm_fftlog_get_full_size (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return self->Nf;
}

/**
 * ncm_fftlog_get_norma:
 * @fftlog: a #NcmFftlog
 *
 * Gets the normalization of the discrete transform, the full size $N_f^\prime$.
 *
 * Returns: $N_f^\prime$ as a double
 */
gdouble
ncm_fftlog_get_norma (NcmFftlog *fftlog)
{
  return ncm_fftlog_get_full_size (fftlog);
}

/**
 * ncm_fftlog_get_length:
 * @fftlog: a #NcmFftlog
 *
 * Gets the value of the ``physical'' period, i.e., period of the fundamental interval.
 *
 * Returns: the period $L$.
 */
gdouble
ncm_fftlog_get_length (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return self->Lk;
}

/**
 * ncm_fftlog_get_full_length:
 * @fftlog: a #NcmFftlog
 *
 * Gets the value of the total period, i.e., period defined by the fundamental interval plus the padding size.
 *
 * Returns: the total period $L_T$.
 */
gdouble
ncm_fftlog_get_full_length (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return self->Lk + 2.0 * self->Lk_N * self->pad;
}

/**
 * ncm_fftlog_get_mode_index:
 * @fftlog: a #NcmFftlog
 * @i: index
 *
 * Gets the mode $n$ held in slot @i of the FFT array: $n = i$ up to $N_f^\prime/2$ and
 * $n = i - N_f^\prime$ above, the $n$ of the decomposition on the theory page.
 *
 * Returns: the mode $n$
 */
gint
ncm_fftlog_get_mode_index (NcmFftlog *fftlog, gint i)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return (i > self->Nf_2) ? i - self->Nf : i;
}

/**
 * ncm_fftlog_get_array_index:
 * @fftlog: a #NcmFftlog
 * @phys_i: index
 *
 * Gets the slot of the FFT array holding mode @phys_i, the inverse of
 * ncm_fftlog_get_mode_index().
 *
 * Returns: the array index of mode @phys_i
 */
gint
ncm_fftlog_get_array_index (NcmFftlog *fftlog, gint phys_i)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return (phys_i < 0) ? phys_i + self->Nf : phys_i;
}

/**
 * ncm_fftlog_peek_output_vector:
 * @fftlog: a #NcmFftlog
 * @nderiv: derivative number
 *
 * Peeks the output vector of $G(r)$, @nderiv = 0, or of its @nderiv-th derivative with
 * respect to $\ln r$.
 *
 * Returns: (transfer none): the output vector
 */
NcmVector *
ncm_fftlog_peek_output_vector (NcmFftlog *fftlog, guint nderiv)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return g_ptr_array_index (self->Gr_vec, nderiv);
}

