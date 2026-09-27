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
 * With a bias $b$ (#NcmFftlog:bias) the expansion is of $F(k) k^{-b}$ instead, and
 * the kernel becomes $(kr)^b K(kr)$. The result is the same $G(r)$; the bias only
 * chooses which function the discrete transform represents, see
 * ncm_fftlog_set_bias().
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
  gdouble lnr0_grid;
  gdouble eval_r_min;
  gdouble eval_r_max;
  gdouble Lk;
  gdouble Lk_N;
  gdouble pad_p;
  gdouble bias;
  gdouble end_slope_min;
  gdouble end_slope_max;
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
  PROP_BIAS,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmFftlog, ncm_fftlog, G_TYPE_OBJECT)

static void
ncm_fftlog_init (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  self->lnr0           = 0.0;
  self->lnr0_grid      = 0.0;
  self->use_eval_int   = FALSE;
  self->smooth_padding = FALSE;
  self->eval_r_min     = 0.0;
  self->eval_r_max     = 0.0;
  self->lnk0           = 0.0;
  self->Lk             = 0.0;
  self->Lk_N           = 0.0;
  self->pad_p          = 0.0;
  self->bias           = 0.0;
  self->end_slope_min  = GSL_NAN;
  self->end_slope_max  = GSL_NAN;
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
    case PROP_BIAS:
      ncm_fftlog_set_bias (fftlog, g_value_get_double (value));
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
    case PROP_BIAS:
      g_value_set_double (value, ncm_fftlog_get_bias (fftlog));
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

static void _ncm_fftlog_get_bias_range_none (NcmFftlog *fftlog, gdouble *bias_min, gdouble *bias_max);

static void
ncm_fftlog_class_init (NcmFftlogClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_fftlog_set_property;
  object_class->get_property = &_ncm_fftlog_get_property;
  object_class->finalize     = &_ncm_fftlog_finalize;

  klass->get_bias_range = &_ncm_fftlog_get_bias_range_none;

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

  /**
   * NcmFftlog:bias:
   *
   * The bias $b$: the transform expands $F(k) k^{-b}$, see ncm_fftlog_set_bias().
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_BIAS,
                                   g_param_spec_double ("bias",
                                                        NULL,
                                                        "Bias",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/* A kernel that does not declare its range supports no bias. */
static void
_ncm_fftlog_get_bias_range_none (NcmFftlog *fftlog, gdouble *bias_min, gdouble *bias_max)
{
  *bias_min = 0.0;
  *bias_max = 0.0;
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
 * Gets the requested center of the transform output. With the no-ringing adjustment
 * (ncm_fftlog_set_noring()) the output grid is centred up to one knot away from it, see
 * ncm_fftlog_get_vector_lnr().
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

/* The fundamental and padding knot counts set_size() gives @n knots, without allocating. */
static void
_ncm_fftlog_size_for (NcmFftlogPrivate * const self, guint n, gint *n_new, guint *pad)
{
  const gint nt = ncm_util_fact_size ((gint) (n * (1.0 + self->pad_p)));

  *pad   = nt * self->pad_p * 0.5 / (1.0 + self->pad_p);
  *n_new = nt - 2 * (gint) * pad;
}

void
ncm_fftlog_set_size (NcmFftlog *fftlog, guint n)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);
  gint n_new;

  self->Nr = n;

  _ncm_fftlog_size_for (self, n, &n_new, &self->pad);

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
 * Sets whether to use Hamilton's low-ringing adjustment, which moves $\ln r_0$ by less
 * than one knot so that the kernel coefficient of the Nyquist mode is real. It applies
 * only when ncm_fftlog_get_full_size() is even, since an odd size has no Nyquist mode.
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
 * ncm_fftlog_set_bias:
 * @fftlog: a #NcmFftlog
 * @bias: the bias $b$
 *
 * Sets the bias $b$. The transform then expands $f(k) = F(k) k^{-b}$ and uses the kernel
 * $(kr)^b K(kr)$, so that
 * \begin{equation*}
 * G(r) = r^{-b} \int_0^\infty f(k)\,(kr)^b K(kr)\,\mathrm{d}k
 * \end{equation*}
 * is the same function for any $b$. The coefficients become
 * $Y_n = \int_0^\infty t^{b + 2\pi i n / L_T} K(t)\,\mathrm{d}t$, which exist only for
 * $b$ inside the range of ncm_fftlog_get_bias_range(); the evaluation aborts otherwise.
 *
 * What the bias changes is the discrete representation: the padding holds $f$, not
 * $F$, and the roundoff of the transform scales with the dynamic range of $f$ over the
 * interval and padding. A bias close to the log-slope of $F$ at the ends keeps $f$ from
 * growing into the padding; for a power law $F \propto k^b$ it makes $f$ constant.
 *
 */
void
ncm_fftlog_set_bias (NcmFftlog *fftlog, const gdouble bias)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  if (self->bias != bias)
  {
    self->bias = bias;
    ncm_fftlog_reset (fftlog);
  }
}

/**
 * ncm_fftlog_get_bias:
 * @fftlog: a #NcmFftlog
 *
 * Returns: the bias $b$, see ncm_fftlog_set_bias()
 */
gdouble
ncm_fftlog_get_bias (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  return self->bias;
}

/**
 * ncm_fftlog_get_bias_range:
 * @fftlog: a #NcmFftlog
 * @bias_min: (out): lower end of the range
 * @bias_max: (out): upper end of the range
 *
 * Gets the open range $(b_\mathrm{min}, b_\mathrm{max})$ of biases for which the kernel
 * coefficients exist, set by the behaviour of $t^b K(t)$ at $t \to 0$ and $t \to \infty$.
 * A kernel that supports no bias returns $b_\mathrm{min} = b_\mathrm{max} = 0$, and then
 * only $b = 0$ is accepted.
 *
 */
void
ncm_fftlog_get_bias_range (NcmFftlog *fftlog, gdouble *bias_min, gdouble *bias_max)
{
  NCM_FFTLOG_GET_CLASS (fftlog)->get_bias_range (fftlog, bias_min, bias_max);
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
 * inverse of that end at a level that falls only as $1/N$. The padding continues the
 * function the transform expands, $f = F k^{-b}$ with the bias of ncm_fftlog_set_bias().
 * The continuation is the power law of each end, with the value and log-slope of $f$
 * there taken from a cubic through the four nearest knots, anchored at the end itself so
 * that it does not move with the knots. Where $F k$, the integrand in $\ln k$, would grow along the continuation, the
 * power law is cut by a Gaussian in $\ln k$ of width the inverse of that log-slope, so
 * that the padding never holds much more integral than the interval does. The two
 * continuations are joined by a $C^\infty$ partition of unity that
 * switches over the middle fifth of the padding, so each side keeps its own continuation
 * over two fifths of the padding. The periodic input is then continuous with its first
 * derivative at each end of the interval, where the continuation matches the value and
 * log-slope of $F$ but not its higher derivatives, and smooth elsewhere; the transform
 * converges as $N^{-3}$ over the whole output grid. The continuation is fitted
 * in log space, so $F$ must be positive at the four knots nearest each end. The fitted
 * slopes are kept, see ncm_fftlog_get_end_slopes().
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

static void _ncm_fftlog_check_bias (NcmFftlog *fftlog);

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
    gdouble shift  = 0.0;
    fftw_complex *Ym_ndm1;

    _ncm_fftlog_check_bias (fftlog);
    NCM_FFTLOG_GET_CLASS (fftlog)->compute_Ym (fftlog, Ym_0);

    /* Hamilton's low-ringing condition sets the phase of the Nyquist mode, which only an
     * even Nf has. theta does not depend on ln r0, so one step makes M an integer. */
    if (self->noring && ((self->Nf % 2) == 0))
    {
      const gdouble theta = carg (Ym_0[self->Nf / 2]);
      const gdouble M     = (self->Nf / Lt) * lnr0k0 - theta / M_PI;
      const gdouble dM    = M - (glong) M;

      shift   = -(Lt / self->Nf) * dM;
      lnr0k0 += shift;
    }

    /* The shift applies to the grid, never to the requested ln r0: otherwise each new size
     * shifts from the last one and the grid depends on the sizes visited before. */
    self->lnr0_grid = self->lnr0 + shift;

    for (i = 0; i < self->Nf; i++)
    {
      const gint phys_i      = ncm_fftlog_get_mode_index (fftlog, i);
      const complex double a = twopi_Lt * phys_i * I;

      Ym_0[i] *= cexp (-a * lnr0k0);

      Ym_ndm1 = Ym_0;

      /* G(r) is r^(-1 - b) times the sum over the modes r^(-a), so each derivative in
       * ln r brings down -(1 + b + a). */
      for (nd = 1; nd <= self->nderivs; nd++)
      {
        fftw_complex *Ym_nd = g_ptr_array_index (self->Ym, nd);

        Ym_nd[i] = -(1.0 + self->bias + a) * Ym_ndm1[i];
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
    const gdouble lnr = self->lnr0_grid + phys_i * self->Lk_N;

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
      const gdouble rm1mb   = exp (-(1.0 + self->bias) * lnr);
      const gdouble Gr_nd_i = creal (self->Gr[i + self->pad]) * rm1mb / norma;

      ncm_vector_set (Gr_nd, i, Gr_nd_i);
    }
  }

  self->evaluated = TRUE;
}

static void
_ncm_fftlog_check_bias (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);
  gdouble bias_min, bias_max;

  if (self->bias == 0.0)
    return;

  ncm_fftlog_get_bias_range (fftlog, &bias_min, &bias_max);

  if ((bias_min == bias_max) || !((self->bias > bias_min) && (self->bias < bias_max)))
    g_error ("ncm_fftlog: the kernel of `%s' has coefficients only for a bias in (%g, %g), got %g.",
             G_OBJECT_TYPE_NAME (fftlog), bias_min, bias_max, self->bias);
}

/**
 * ncm_fftlog_get_Ym:
 * @fftlog: a #NcmFftlog
 * @size: (out): return size
 *
 * Computes the kernel coefficients $Y_m$ at the current bias, before the phase factor of
 * the grid centres, as interleaved real and imaginary parts of the
 * ncm_fftlog_get_full_size() modes.
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

  _ncm_fftlog_check_bias (fftlog);
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

/**
 * ncm_fftlog_get_end_slopes:
 * @fftlog: a #NcmFftlog
 * @slope_min: (out): $\mathrm{d}\ln F / \mathrm{d}\ln k$ at $\ln k_0 - L/2$
 * @slope_max: (out): $\mathrm{d}\ln F / \mathrm{d}\ln k$ at $\ln k_0 + L/2$
 *
 * Gets the log-slopes of $F$ at the two ends of the interval, as the smooth padding
 * fitted them in the last evaluation (see ncm_fftlog_use_smooth_padding()). They are
 * slopes of $F$, not of $F k^{-b}$, so they do not depend on the bias. The evaluation
 * must have used the smooth padding.
 *
 */
void
ncm_fftlog_get_end_slopes (NcmFftlog *fftlog, gdouble *slope_min, gdouble *slope_max)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);

  if (!self->evaluated || !self->smooth_padding)
    g_error ("ncm_fftlog_get_end_slopes: the slopes come from the last evaluation with the smooth padding.");

  *slope_min = self->end_slope_min;
  *slope_max = self->end_slope_max;
}

/**
 * ncm_fftlog_get_best_bias:
 * @fftlog: a #NcmFftlog
 *
 * Chooses a bias from the end slopes of the last evaluation, see
 * ncm_fftlog_get_end_slopes(). With $s_\mathrm{min}$ and $s_\mathrm{max}$ the slopes at
 * $k_\mathrm{min}$ and $k_\mathrm{max}$, $f = F k^{-b}$ grows into neither padding when
 * $s_\mathrm{max} \le b \le s_\mathrm{min}$; the choice is the $b$ closest to zero in
 * that interval, so an input that already decays at both ends keeps $b = 0$, and a power
 * law $F \propto k^s$ gets $b = s$. When the interval is empty, $F$ grows at both ends
 * and the choice is its midpoint, which splits the growth between them.
 *
 * The periodic image of the input one period $L_T$ away enters with a weight
 * $e^{-(b - b_\mathrm{min}) L_T}$ at one end of ncm_fftlog_get_bias_range() and
 * $e^{-(b_\mathrm{max} - b) L_T}$ at the other, so the choice is then kept
 * $\ln(1/\epsilon) / L_T$ inside the range, with $\epsilon$ the double precision
 * epsilon, where both weights are below roundoff. A range narrower than that gives its
 * midpoint, and a kernel that supports no bias gives zero.
 *
 * Returns: the bias
 */
gdouble
ncm_fftlog_get_best_bias (NcmFftlog *fftlog)
{
  const gdouble margin = -log (GSL_DBL_EPSILON) / ncm_fftlog_get_full_length (fftlog);
  gdouble slope_min, slope_max, bias_min, bias_max, bias;

  ncm_fftlog_get_end_slopes (fftlog, &slope_min, &slope_max);
  ncm_fftlog_get_bias_range (fftlog, &bias_min, &bias_max);

  if (bias_min == bias_max)
    return 0.0;

  if (slope_max <= slope_min)
    bias = GSL_MIN (GSL_MAX (0.0, slope_max), slope_min);
  else
    bias = 0.5 * (slope_min + slope_max);

  if (bias_max - bias_min <= 2.0 * margin)
    return 0.5 * (bias_min + bias_max);

  return GSL_MIN (GSL_MAX (bias, bias_min + margin), bias_max - margin);
}

static void _ncm_fftlog_end_power_law (const gdouble u[4], const gdouble lnF[4], gdouble *A, gdouble *s);
static gdouble _ncm_fftlog_continuation (const gdouble A, const gdouble s, const gdouble sigma, const gdouble u);
static gdouble _ncm_fftlog_smooth_step (const gdouble t);

/* Fills the padding with the continuation described in ncm_fftlog_use_smooth_padding().
 * The padding is one stretch of 2 pad slots in the periodic array, running from just
 * above ln k_max (slot pad + N) around the wrap to just below ln k_min (slot pad - 1);
 * t is the position along it in ln k, from 0 to 1, and the partition switches over its
 * middle fifth. Measuring t in ln k, not in slots, keeps the partition where it is when
 * the grid is refined: an index ratio moves it by O(1/N) and so does the result. */
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

  /* The knots hold f = F k^(-b): its slope in u is that of F less b above the interval
   * and plus b below it, where u runs down in ln k. */
  self->end_slope_max = s_hi + self->bias;
  self->end_slope_min = self->bias - s_lo;

  for (i = 0; i < stretch; i++)
  {
    const gdouble x_hi = lnk_last + (i + 1) * self->Lk_N - lnk_max;
    const gdouble x_lo = lnk_min - (lnk_first - (stretch - i) * self->Lk_N);
    const gdouble t    = x_hi / (stretch * self->Lk_N); /* x_hi + x_lo = stretch Lk_N */
    const gdouble w    = _ncm_fftlog_smooth_step ((t - 0.4) / 0.2);
    const gdouble F_i  = (1.0 - w) * _ncm_fftlog_continuation (A_hi, s_hi, s_hi + 1.0 + self->bias, x_hi)
                         + w * _ncm_fftlog_continuation (A_lo, s_lo, s_lo - 1.0 - self->bias, x_lo);

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
 * of width 1 / sigma when sigma > 0. The knots hold f = F k^(-b) and s is its log-slope
 * in u; sigma is that of F k = f k^(1 + b), what the transform integrates in ln k:
 * s + 1 + b above the interval, where ln k = ln k_max + u, and s - 1 - b below it, where
 * ln k = ln k_min - u. It does not depend on b. A cut on the growth of f alone lets F k
 * grow for several e-foldings and its periodic image floods the output far below its
 * peak; one on s below the interval also cuts tails where f grows but F k decays. */
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
    const gint phys_i   = i - self->N_2;
    const gdouble lnk_i = self->lnk0 + self->Lk_N * phys_i;

    self->Fk[self->pad + i] = ncm_vector_get (Fk, i) * exp (-self->bias * lnk_i);
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

    self->Fk[self->pad + i] = Fk_i * exp (-self->bias * lnk_i);
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
 * Each step computes one transform; the next size is checked against #NcmFftlog:max-n
 * before anything is allocated, and the calibration aborts if it would pass it first. A
 * component that is zero at both sizes counts as converged, and a non-finite transform
 * aborts. Leaves @fftlog evaluated at the final size.
 *
 */
void
ncm_fftlog_calibrate_size_gsl (NcmFftlog *fftlog, gsl_function *Fk, const gdouble reltol)
{
  NcmFftlogPrivate * const self = ncm_fftlog_get_instance_private (fftlog);
  NcmSpline **prev              = g_new0 (NcmSpline *, self->nderivs + 1);
  gdouble err                   = GSL_NAN; /* Change at the last step, NaN before the first */

  ncm_fftlog_eval_by_gsl_function (fftlog, Fk);
  ncm_fftlog_prepare_splines (fftlog);

  while (TRUE)
  {
    const gint N_prev = self->N;
    gint N_next       = N_prev;
    guint n_req       = N_prev;
    guint pad_next, nd, i, size, n_try;
    NcmVector *lnr;

    /* Grow by 20%; the factorable full size can round back to the same N, so grow further
     * until it moves. The limit is checked on the size before anything is allocated. */
    for (n_try = 1; N_next <= N_prev; n_try++)
    {
      n_req = (guint) ceil (N_prev * pow (1.2, n_try));
      _ncm_fftlog_size_for (self, n_req, &N_next, &pad_next);
    }

    if ((N_next > (gint) self->max_n) && isnan (err))
      g_error ("ncm_fftlog_calibrate_size_gsl: the next size, %d knots, exceeds the maximum (%u) "
               "before any comparison; the requested relative accuracy is %e.",
               N_next, self->max_n, reltol);
    else if (N_next > (gint) self->max_n)
      g_error ("ncm_fftlog_calibrate_size_gsl: the next size, %d knots, exceeds the maximum (%u) "
               "before the requested relative accuracy %e was reached; the last change was %e.",
               N_next, self->max_n, reltol, err);

    /* Keep the current splines: set_size() replaces them when the size changes. */
    for (nd = 0; nd <= self->nderivs; nd++)
      prev[nd] = ncm_spline_ref (g_ptr_array_index (self->Gr_s, nd));

    ncm_fftlog_set_size (fftlog, n_req);
    ncm_fftlog_eval_by_gsl_function (fftlog, Fk);
    ncm_fftlog_prepare_splines (fftlog);

    /* Largest difference between the two sizes on the new knots, relative to the value
     * plus the component's peak. A component zero at both sizes has not changed. */
    lnr  = ncm_spline_get_xv (g_ptr_array_index (self->Gr_s, 0));
    size = ncm_spline_get_len (g_ptr_array_index (self->Gr_s, 0));
    err  = 0.0;

    for (nd = 0; nd <= self->nderivs; nd++)
    {
      NcmVector *Gr = ncm_spline_get_yv (g_ptr_array_index (self->Gr_s, nd));
      gdouble absmin, absmax;

      ncm_vector_get_absminmax (Gr, &absmin, &absmax);

      for (i = 0; i < size; i++)
      {
        const gdouble G_new  = ncm_vector_get (Gr, i);
        const gdouble G_prev = ncm_spline_eval (prev[nd], ncm_vector_get (lnr, i));
        const gdouble scale  = fabs (G_new) + absmax;

        if (!isfinite (G_new) || !isfinite (G_prev))
          g_error ("ncm_fftlog_calibrate_size_gsl: the transform is not finite (component %u, "
                   "ln r = %g, %g at %d knots and %g before).",
                   nd, ncm_vector_get (lnr, i), G_new, self->N, G_prev);

        if (scale > 0.0)
          err = GSL_MAX (err, fabs (G_new - G_prev) / scale);
        else if (G_prev != 0.0)
          err = GSL_MAX (err, 1.0);
      }

      ncm_spline_clear (&prev[nd]);
      ncm_vector_free (Gr);
    }

    ncm_vector_free (lnr);

    if (err <= reltol)
      break;
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
 * Each step computes one transform; the next size is checked against #NcmFftlog:max-n
 * before anything is allocated, and the calibration aborts if it would pass it first. A
 * component that is zero at both sizes counts as converged, and a non-finite transform
 * aborts. Leaves @fftlog evaluated at the final size.
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

