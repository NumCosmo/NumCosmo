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
 * @stability: Stable
 * @include: numcosmo/math/ncm_fftlog.h
 *
 * This class provides the tools to compute the Fast Fourier Transform
 * of any function, which is assumed to be a periodic
 * sequence of logarithmically spaced points.
 * It is inspired on the approach FFTLog developed by
 * [Hamilton (2000)][XHamilton2000] [[arXiv](https://arxiv.org/abs/astro-ph/9905191)],
 * which was extended as described below.
 *
 * A function $G(r)$ is written as
 * \begin{equation}\label{eq:Gr} G(r) = \int_0^\infty F(k) \ K(kr) dk, \end{equation}
 * where $F(k)$ is defined in the fundamental interval $[\ln k_0 - L/2, \ln k_0 + L/2]$, $L$ is the period,
 * $\ln k_0$ is the center value and $K(kr)$ is a kernel function. Assuming that $F(k)$ can be written in terms of the
 * $N$ lowest Fourier modes, we have
 *
 * $$F(k) = \sum_{n} c_n e^{\frac{2\pi i n}{L} \ln\left(\frac{k}{k_0}\right)}.$$
 * Substituting $F(k)$ in Eq. \eqref{eq:Gr} and changing the variable $k \rightarrow t = kr$, thus
 * \begin{align}\label{eq:Gr_decomp}
 * r G(r) &= \sum_n c_n \int_0^\infty \frac{k}{k_0}^{\frac{2\pi i n}{L}} K(kr)^2 d(kr) \\
 * &= \sum_n c_n  \int_0^\infty \frac{t}{k_0 r}^{\frac{2\pi i n}{L}} K(t) dt \\
 * &= \sum_n c_n e^{-\frac{2\pi i n}{L} \ln\left(\frac{r}{r_0}\right)} e^{-\frac{2\pi i n}{L} \ln(k_0 r_0)} Y_n,
 * \end{align}
 * where
 * $$Y_n = \int_0^\infty t^{\frac{2\pi i n}{L}} K(t) dt,$$
 * and the Fourier coefficients are
 * $$c_n = \frac{1}{N} \sum_m F(k_m) e^{- \frac{2\pi i nm}{N}}.$$
 * The total number of points $N$ corresponds to the number of knots in the fundamental interval, which is equally spaced.
 *
 * The variables discretization is different depending whether $N$ is even or odd, in general
 * $$ k_n = k_0 \mathrm{e}^{n L / N}, \qquad r_m = r_0 \mathrm{e}^{m L / N}. $$
 * If $N$ is odd, $n$ and $m$ runs from $[-\lfloor N/2\rfloor, \lfloor N/2\rfloor]$,
 * where $\lfloor N/2\rfloor$  is the round-down (largest integer smaller than $N/2$)
 * of $N/2$. In this case
 * \begin{align}
 * \ln\left(k_{-\lfloor N/2\rfloor}\right) = \ln(k_0) - \frac{(N-1)}{N} \frac{L}{2}, \\
 * \ln\left(k_{+\lfloor N/2\rfloor}\right) = \ln(k_0) + \frac{(N-1)}{N} \frac{L}{2}.
 * \end{align}
 * This means that for odd $N$ the values of $k_n$ (and $r_n$) never touches the borders $\ln(k_0) \pm L/2$
 * and $\ln(r_0) \pm L/2$. On the other hand if $N$ is even ($\lfloor N/2\rfloor = N/2$)
 * \begin{align}
 * \ln\left(k_{-\lfloor N/2\rfloor}\right) = \ln(k_0) - \frac{L}{2}, \\
 * \ln\left(k_{+\lfloor N/2\rfloor}\right) = \ln(k_0) + \frac{L}{2}.
 * \end{align}
 * However, since we are assuming that these functions are periodic with period $L$ these two
 * points refer to the same value of the functions. Thus, we do not need to include both points and in
 * the case of even $N$ we include the point $\ln(k_0) - \frac{L}{2}$ only. In the original [FFTLog][XHamilton2000]
 * paper they include both points but give them a $1/2$ weight, here we avoid this complication by using the lower
 * end point only.
 *
 * The user must provide the following input values: $\ln k_0$ - ncm_fftlog_set_lnk0(), $\ln r_0$ - ncm_fftlog_set_lnr0(),
 * $L$ - ncm_fftlog_set_length(), padding percentage - ncm_fftlog_set_padding(), $N$ - ncm_fftlog_set_size(),
 * $F(k)$ (or $F(k_m)$ -- see description below).
 *
 * - Since the algorithm assumes that the function to be decomposed is periodic, it is worth extending the interval in $\ln k$ such that
 * $F(k) \equiv 0$ in the intervals $\left[\ln k_0 -\frac{L_T}{2}, \ln k_0 - \frac{L}{2} \right)$ and
 * $ \left(\ln k_0 + \frac{L}{2}, \ln k_0 + \frac{L_T}{2}\right]$, where the total period $L_T$ is defined by the final
 * number of knots, i.e., $N_f = N (1 + \mathrm{padding})$.
 * - $N$ knots are equally distributed in the fundamental interval and $N \times \mathrm{padding}$ knots are distributed
 * in the two simetric intervals as mentioned above.
 * - For the sake of optimization, the final number of points $N_f$ is substituted by the smallest number $N_f^\prime$ (bigger than $N_f$)
 * which can be decomposed as $N_f \leq N_f^\prime = N^\prime (1 + \mathrm{padding}) = 2^a 3^b 5^c 7^d$, where $a$,
 * $b$, $c$ and $d$ are positive integers.
 * - The function $F(k)$ can be provided as:
 * 1. a gsl_function - ncm_fftlog_eval_by_gsl_function() - whose values are computed at the knots
 * within the fundamental interval, and set to zero within the padding intervals.
 * 2. as a vector - ncm_fftlog_eval_by_vector() - first one must get the vector of $\ln k$ knots, ncm_fftlog_get_lnk_vector(),
 * and then pass a vector containing the values of the function computed at each knot.
 *
 * - Regarding $Y_n$, see the different implementations of #NcmFftlog, e.g., #NcmFftlogTophatwin2 and #NcmFftlogGausswin2.
 *
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fftlog.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"
#include "math/ncm_spline_cubic_notaknot.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <complex.h>
#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */
#endif /* NUMCOSMO_GIR_SCAN */

#ifndef HAVE_FFTW3_ALLOC
#define fftw_alloc_real(n) (double *) fftw_malloc (sizeof (double) * (n))
#define fftw_alloc_complex(n) (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * (n))
#endif /* HAVE_FFTW3_ALLOC */

struct _NcmFftlogPrivate
{
  gint Nr;
  gint N;
  gint N_2;
  gint Nf;
  gint Nf_2;
  guint nderivs;
  guint pad;
  gdouble lnk0;
  gdouble lnr0;
  gdouble eval_r_min;
  gdouble eval_r_max;
  gdouble Lk;
  gdouble Lk_N;
  gdouble pad_p;
  gdouble smooth_padding_scale;
  gboolean smooth_padding;
  gboolean use_eval_int;
  gboolean noring;
  gboolean prepared;
  gboolean evaluated;
  NcmVector *lnr_vec;
  GPtrArray *Gr_vec;
  GPtrArray *Gr_s;
  
#ifdef NUMCOSMO_HAVE_FFTW3
  fftw_complex *Fk;
  fftw_complex *Cm;
  fftw_complex *Gr;
  fftw_complex *CmYm;
  GPtrArray *Ym;
  fftw_plan p_Fk2Cm;
  fftw_plan p_CmYm2Gr;
#endif /* NUMCOSMO_HAVE_FFTW3 */
};

enum
{
  PROP_0,
  PROP_NDERIV,
  PROP_LNR0,
  PROP_LNK0,
  PROP_LR,
  PROP_N,
  PROP_PAD,
  PROP_NORING,
  PROP_NAME,
  PROP_USE_EVAL_INT,
  PROP_SMOOTH_PADDING,
  PROP_SMOOTH_PADDING_SCALE,
  PROP_EVAL_R_MIN,
  PROP_EVAL_R_MAX,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmFftlog, ncm_fftlog, G_TYPE_OBJECT);

static void
ncm_fftlog_init (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = fftlog->priv = ncm_fftlog_get_instance_private (fftlog);
  
  self->lnr0                 = 0.0;
  self->use_eval_int         = FALSE;
  self->smooth_padding       = FALSE;
  self->smooth_padding_scale = 0.0;
  self->eval_r_min           = 0.0;
  self->eval_r_max           = 0.0;
  self->lnk0                 = 0.0;
  self->Lk                   = 0.0;
  self->Lk_N                 = 0.0;
  self->pad_p                = 0.0;
  self->Nr                   = 0;
  self->N                    = 0;
  self->N_2                  = 0;
  self->Nf                   = 0;
  self->Nf_2                 = 0;
  self->pad                  = 0;
  self->noring               = FALSE;
  self->prepared             = FALSE;
  self->evaluated            = FALSE;
  
  self->lnr_vec = NULL;
  self->Gr_vec  = g_ptr_array_new ();
  self->Gr_s    = g_ptr_array_new ();
  
  g_ptr_array_set_free_func (self->Gr_vec, (GDestroyNotify) ncm_vector_free);
  g_ptr_array_set_free_func (self->Gr_s, (GDestroyNotify) ncm_spline_free);
  
#ifdef NUMCOSMO_HAVE_FFTW3
  self->Fk        = NULL;
  self->Cm        = NULL;
  self->Gr        = NULL;
  self->Ym        = g_ptr_array_new ();
  self->CmYm      = NULL;
  self->p_Fk2Cm   = NULL;
  self->p_CmYm2Gr = NULL;
  g_ptr_array_set_free_func (self->Ym, (GDestroyNotify) fftw_free);
#endif /* NUMCOSMO_HAVE_FFTW3 */
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
    case PROP_SMOOTH_PADDING_SCALE:
      ncm_fftlog_set_smooth_padding_scale (fftlog, g_value_get_double (value));
      break;
    case PROP_EVAL_R_MIN:
      ncm_fftlog_set_eval_r_min (fftlog, g_value_get_double (value));
      break;
    case PROP_EVAL_R_MAX:
      ncm_fftlog_set_eval_r_max (fftlog, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fftlog_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFftlog *fftlog             = NCM_FFTLOG (object);
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
    case PROP_SMOOTH_PADDING_SCALE:
      g_value_set_double (value, ncm_fftlog_get_smooth_padding_scale (fftlog));
      break;
    case PROP_EVAL_R_MIN:
      g_value_set_double (value, ncm_fftlog_get_eval_r_min (fftlog));
      break;
    case PROP_EVAL_R_MAX:
      g_value_set_double (value, ncm_fftlog_get_eval_r_max (fftlog));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

#ifdef NUMCOSMO_HAVE_FFTW3

static void
_ncm_fftlog_free_all (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  
  g_clear_pointer (&self->Fk, fftw_free);
  g_clear_pointer (&self->Cm, fftw_free);
  g_clear_pointer (&self->CmYm, fftw_free);
  g_clear_pointer (&self->Gr, fftw_free);
  
  g_clear_pointer (&self->p_Fk2Cm, fftw_destroy_plan);
  g_clear_pointer (&self->p_CmYm2Gr, fftw_destroy_plan);
  
  ncm_vector_clear (&self->lnr_vec);
  
  g_ptr_array_set_size (self->Gr_vec, 0);
  g_ptr_array_set_size (self->Gr_s, 0);
  g_ptr_array_set_size (self->Ym, 0);
}

#endif /* NUMCOSMO_HAVE_FFTW3 */

static void
_ncm_fftlog_finalize (GObject *object)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  NcmFftlog *fftlog             = NCM_FFTLOG (object);
  NcmFftlogPrivate * const self = fftlog->priv;
  
  _ncm_fftlog_free_all (fftlog);
  
  g_clear_pointer (&self->Gr_vec, g_ptr_array_unref);
  g_clear_pointer (&self->Gr_s,   g_ptr_array_unref);
  g_clear_pointer (&self->Ym,     g_ptr_array_unref);
  
#endif /* NUMCOSMO_HAVE_FFTW3 */
  
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
   * The Center value for $\ln(r)$.
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
   * The Center value for $\ln(k)$.
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
   * The function $F(k)$'s period in natural logarithm base.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LR,
                                   g_param_spec_double ("Lk",
                                                        NULL,
                                                        "Function log-period",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.0,
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
   * NcmFftlog:padding:
   *
   * The padding percentage of the number of knots $N$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_PAD,
                                   g_param_spec_double ("padding",
                                                        NULL,
                                                        "Padding percentage",
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
   * FFTW Plan wisdown's name to perform the transformation.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_NAME,
                                   g_param_spec_string ("name",
                                                        NULL,
                                                        "FFTW Plan wisdown name",
                                                        "fftlog_default_wisdown",
                                                        G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_USE_EVAL_INT,
                                   g_param_spec_boolean ("use-eval-int",
                                                         NULL,
                                                         "Whether to use evaluation interval",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SMOOTH_PADDING,
                                   g_param_spec_boolean ("use-smooth-padding",
                                                         NULL,
                                                         "Whether to use a smooth padding",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SMOOTH_PADDING_SCALE,
                                   g_param_spec_double ("smooth-padding-scale",
                                                        NULL,
                                                        "Log10 of the smoothing scale",
                                                        -G_MAXDOUBLE, +G_MAXDOUBLE, -200.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_EVAL_R_MIN,
                                   g_param_spec_double ("eval-r-min",
                                                        NULL,
                                                        "Evaluation r_min",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
 * This function peeks the @fftlog's associated name.
 *
 * Returns: (transfer none): The internal string describing #NcmFftlog.
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
          NcmSpline *Gr_s_i   = ncm_spline_cubic_notaknot_new_full (self->lnr_vec, Gr_vec_i, FALSE);
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
  return self->lnk0;
}

/**
 * ncm_fftlog_set_size:
 * @fftlog: a #NcmFftlog
 * @n: number of knots
 *
 * Sets the number of knots $N_f^\prime$ where the integrated function is evaluated,
 * given the input number of knots @n, plus padding.
 *
 */
void
ncm_fftlog_set_size (NcmFftlog *fftlog, guint n)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  NcmFftlogPrivate * const self = fftlog->priv;
  guint nt                      = n * (1.0 + self->pad_p);
  
  self->Nr = n;
  
  nt        = ncm_util_fact_size (nt);
  self->pad = nt * self->pad_p * 0.5 / (1.0 + self->pad_p);
  n         = nt - 2 * self->pad;
  
  if ((n != self->N) || (n + 2 * self->pad != self->Nf))
  {
    guint i;
    
    self->N    = n;
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
    
    ncm_cfg_load_fftw_wisdom ("ncm_fftlog_%s", NCM_FFTLOG_GET_CLASS (fftlog)->name);
    
    ncm_cfg_lock_plan_fftw ();
    
    self->p_Fk2Cm   = fftw_plan_dft_1d (self->Nf, self->Fk,   self->Cm, FFTW_FORWARD, fftw_default_flags | FFTW_DESTROY_INPUT);
    self->p_CmYm2Gr = fftw_plan_dft_1d (self->Nf, self->CmYm, self->Gr, FFTW_FORWARD, fftw_default_flags | FFTW_DESTROY_INPUT);
    
    for (i = 0; i <= self->nderivs; i++)
    {
      NcmVector *Gr_vec_i = ncm_vector_new (self->N);
      NcmSpline *Gr_s_i   = ncm_spline_cubic_notaknot_new_full (self->lnr_vec, Gr_vec_i, FALSE);
      fftw_complex *Ym_i  = fftw_alloc_complex (self->Nf);
      
      g_ptr_array_add (self->Gr_vec, Gr_vec_i);
      g_ptr_array_add (self->Gr_s, Gr_s_i);
      g_ptr_array_add (self->Ym, Ym_i);
    }
    
    ncm_cfg_unlock_plan_fftw ();
    
    ncm_cfg_save_fftw_wisdom ("ncm_fftlog_%s", NCM_FFTLOG_GET_CLASS (fftlog)->name);
    
    ncm_fftlog_reset (fftlog);
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
 * Gets the padding percentage.
 *
 * Returns: the padding percentage.
 */
gdouble
ncm_fftlog_get_padding (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
 *
 * Returns: whether no-ringing condition is activated.
 */
gboolean
ncm_fftlog_get_noring (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
 * @use_eval_interal: a gboolean
 *
 * Sets whether to use a restricted evaluation interval $[r_\mathrm{min}, r_\mathrm{max}]$.
 * See ncm_fftlog_set_eval_r_min() and ncm_fftlog_set_eval_r_max().
 *
 */
void
ncm_fftlog_use_eval_interval (NcmFftlog *fftlog, gboolean use_eval_interal)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  
  if (use_eval_interal)
  {
    self->use_eval_int = TRUE;
  }
  else
  {
    if (self->use_eval_int)
    {
      gint nd;
      
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
 * @use_smooth_padding: a gboolean
 *
 * Sets whether to use pad the fft using a power-law continuation of the
 * input function which is continuous at the border and drops to a scale
 * determined by ncm_fftlog_set_smooth_padding_scale().
 *
 */
void
ncm_fftlog_use_smooth_padding (NcmFftlog *fftlog, gboolean use_smooth_padding)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  
  self->smooth_padding = use_smooth_padding;
  self->evaluated      = FALSE;
}

/**
 * ncm_fftlog_set_smooth_padding_scale:
 * @fftlog: a #NcmFftlog
 * @log10sc: a gdouble containing the $\log_{10}(s)$ of the smoothing scale $s$
 *
 * Sets the value of the smoothing scale $s$ which is used if ncm_fftlog_use_smooth_padding()
 * is turned on.
 *
 */
void
ncm_fftlog_set_smooth_padding_scale (NcmFftlog *fftlog, gdouble log10sc)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  
  self->smooth_padding_scale = log10sc;
  self->evaluated            = FALSE;
}

/**
 * ncm_fftlog_get_smooth_padding_scale:
 * @fftlog: a #NcmFftlog
 *
 * Gets the log10 of the current value of the smoothing scale $s$.
 *
 * Returns: $\log_{10}(s)$.
 */
gdouble
ncm_fftlog_get_smooth_padding_scale (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  
  return self->smooth_padding_scale;
}

/**
 * ncm_fftlog_set_eval_r_min:
 * @fftlog: a #NcmFftlog
 * @eval_r_min: the value of $r_\mathrm{min}$
 *
 * Sets $r_\mathrm{min}$ to @r_min.
 *
 */
void
ncm_fftlog_set_eval_r_min (NcmFftlog *fftlog, const gdouble eval_r_min)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  
  if (self->eval_r_min != eval_r_min)
    self->eval_r_min = eval_r_min;
}

/**
 * ncm_fftlog_set_eval_r_max:
 * @fftlog: a #NcmFftlog
 * @eval_r_max: the value of $r_\mathrm{max}$
 *
 * Sets $r_\mathrm{max}$ to @r_max.
 *
 */
void
ncm_fftlog_set_eval_r_max (NcmFftlog *fftlog, const gdouble eval_r_max)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
  return self->eval_r_max;
}

#ifdef NUMCOSMO_HAVE_FFTW3

static void
_ncm_fftlog_eval (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  guint nd;
  gint i;
  
  fftw_execute (self->p_Fk2Cm);
  
  if (!self->prepared)
  {
    const gdouble Lt       = ncm_fftlog_get_full_length (fftlog);
    const gdouble twopi_Lt = 2.0 * M_PI / Lt;
    fftw_complex *Ym_0     = g_ptr_array_index (self->Ym, 0);
    gdouble lnr0k0         = self->lnk0 + self->lnr0;
    fftw_complex *Ym_ndm1;
    
    NCM_FFTLOG_GET_CLASS (fftlog)->get_Ym (fftlog, Ym_0);
    
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
    
    for (i = 0; i < self->Nf; i++)
    {
      self->CmYm[i] = self->Cm[i] * Ym_nd[i];
    }
    
    self->CmYm[self->Nf_2]     = creal (self->CmYm[self->Nf_2]);
    self->CmYm[self->Nf_2 + 1] = creal (self->CmYm[self->Nf_2 + 1]);
    
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

#endif /* NUMCOSMO_HAVE_FFTW3 */

/**
 * ncm_fftlog_get_Ym:
 * @fftlog: a #NcmFftlog
 * @size: (out): return size
 *
 * Computes the $Y_m$ vector.
 *
 * Returns: (transfer none) (array length=size): $Y_m$.
 */
gdouble *
ncm_fftlog_get_Ym (NcmFftlog *fftlog, guint *size)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  
  fftw_complex *Ym_0 = g_ptr_array_index (self->Ym, 0);
  
  NCM_FFTLOG_GET_CLASS (fftlog)->get_Ym (fftlog, Ym_0);
  
  size[0] = ncm_fftlog_get_full_size (fftlog) * 2;
  
  return (gdouble *) Ym_0;
}

/**
 * ncm_fftlog_get_lnk_vector:
 * @fftlog: a #NcmFftlog
 * @lnk: a #NcmVector
 *
 * Computes the $\ln k$ vector.
 *
 */
void
ncm_fftlog_get_lnk_vector (NcmFftlog *fftlog, NcmVector *lnk)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  gint i;
  
  g_assert_cmpuint (self->N, ==, ncm_vector_len (lnk));
  
  for (i = 0; i < self->N; i++)
  {
    const gint phys_i   = i - self->N_2;
    const gdouble lnk_i = self->lnk0 + self->Lk_N * phys_i;
    
    ncm_vector_set (lnk, i, lnk_i);
  }
}

static void
_ncm_fftlog_add_smooth_padding (NcmFftlog *fftlog)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  NcmFftlogPrivate * const self = fftlog->priv;
  gint size = 3;
  gdouble dd1[3], xa1[3], ya1[3];
  gdouble dd2[3], xa2[3], ya2[3];
  gint i;
  
  for (i = 0; i < size; i++)
  {
    xa1[i] = log1p (self->pad + 1.0 * i);
    ya1[i] = log (creal (self->Fk[self->pad + i]));
    
    xa2[i] = log1p (self->pad + self->N + 1.0 * i - size);
    ya2[i] = log (creal (self->Fk[self->pad + self->N - size + i]));
  }
  
  xa1[2] = 0.0;
  ya1[2] = self->smooth_padding_scale * M_LN10;
  
  xa2[0] = log1p (2.0 * self->pad + self->N - 1.0);
  ya2[0] = self->smooth_padding_scale * M_LN10;
  
  gsl_poly_dd_init (dd1, xa1, ya1, size);
  gsl_poly_dd_init (dd2, xa2, ya2, size);
  
  for (i = 0; i < self->pad; i++)
  {
    self->Fk[i]                       = exp (gsl_poly_dd_eval (dd1, xa1, size, log1p (1.0 * i)));
    self->Fk[self->pad + self->N + i] = exp (gsl_poly_dd_eval (dd2, xa2, size, log1p (self->pad + self->N + i)));
    
    /*printf ("%8d % 22.15g %8d % 22.15g\n", i, creal (self->Fk[i]), self->pad + self->N + i, creal (self->Fk[self->pad + self->N + i]));*/
  }
  
#endif /* NUMCOSMO_HAVE_FFTW3 */
}

/**
 * ncm_fftlog_eval_by_vector:
 * @fftlog: a #NcmFftlog
 * @Fk: a #NcmVector
 *
 * @Fk is a vector which contains the values of the function at each knot $\ln k_m$.
 *
 */
void
ncm_fftlog_eval_by_vector (NcmFftlog *fftlog, NcmVector *Fk)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  NcmFftlogPrivate * const self = fftlog->priv;
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
#endif /* NUMCOSMO_HAVE_FFTW3 */
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
#ifdef NUMCOSMO_HAVE_FFTW3
  NcmFftlogPrivate * const self = fftlog->priv;
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
#endif /* NUMCOSMO_HAVE_FFTW3 */
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
  NcmFftlogPrivate * const self = fftlog->priv;
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
 * Returns: (transfer full):
 */
NcmVector *
ncm_fftlog_get_vector_lnr (NcmFftlog *fftlog)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
 * Returns: $\frac{\mathrm{d}^nG(r)}{\mathrm{d}\ln r}$ value computed at @lnr.
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
 * Increases the original (input) number of knots until the $G(r)$ splines reach
 * the required precision @reltol.
 *
 */
void
ncm_fftlog_calibrate_size_gsl (NcmFftlog *fftlog, gsl_function *Fk, const gdouble reltol)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  NcmSpline **s                 = g_new0 (NcmSpline *, self->nderivs + 1);
  gdouble lreltol               = 0.0;
  NcmVector *eval_lnr_vec;
  guint nd, size;
  
  ncm_fftlog_eval_by_gsl_function (fftlog, Fk);
  ncm_fftlog_prepare_splines (fftlog);
  
  for (nd = 0; nd <= self->nderivs; nd++)
  {
    s[nd] = ncm_spline_ref (g_ptr_array_index (self->Gr_s, nd));
  }
  
  /*printf ("# Initial size %u [%u].\n", self->N, self->pad);*/
  ncm_fftlog_set_size (fftlog, self->N * 1.2);
  /*printf ("# Trying size %u [%u].\n", self->N, self->pad);*/
  ncm_fftlog_eval_by_gsl_function (fftlog, Fk);
  ncm_fftlog_prepare_splines (fftlog);
  
  eval_lnr_vec = ncm_spline_get_xv (g_ptr_array_index (self->Gr_s, 0));
  size         = ncm_spline_get_len (g_ptr_array_index (self->Gr_s, 0));
  
  for (nd = 0; nd <= self->nderivs; nd++)
  {
    NcmVector *eval_Gr_vec_nd = ncm_spline_get_yv (g_ptr_array_index (self->Gr_s, nd));
    gdouble absmin, absmax;
    guint i /*, i_max = 0*/;
    
    ncm_vector_get_absminmax (eval_Gr_vec_nd, &absmin, &absmax);
    
    /*printf ("# Testing component %u [% 20.15g, % 20.15g].\n", nd, absmin, absmax);*/
    for (i = 0; i < size; i++)
    {
      const gdouble lnr_i     = ncm_vector_get (eval_lnr_vec, i);
      const gdouble lnG_i     = ncm_vector_get (eval_Gr_vec_nd, i);
      const gdouble lnS_i     = ncm_spline_eval (s[nd], lnr_i);
      const gdouble lreltol_i = fabs ((lnG_i - lnS_i) / (fabs (lnG_i) + absmax));
      
      /*printf ("% 20.15g % 20.15e % 20.15e % 20.15e | % 20.15e % 20.15e\n", lnr_i, exp (lnr_i), lnG_i, lnS_i, lreltol_i, fabs ((lnG_i - lnS_i) / fabs (lnG_i)));*/
      if (lreltol_i > lreltol)
        lreltol = lreltol_i;
      
      /*i_max   = i;*/
    }
    
    ncm_spline_clear (&s[nd]);
    ncm_vector_free (eval_Gr_vec_nd);
    /*printf ("# Largest error up to component %u is [%u] %e.\n", nd, i_max, lreltol); fflush (stdout);*/
  }
  
  ncm_vector_free (eval_lnr_vec);
  g_clear_pointer (&s, g_free);
  
  if (lreltol > reltol)
    ncm_fftlog_calibrate_size_gsl (fftlog, Fk, reltol);
}

/**
 * ncm_fftlog_calibrate_size:
 * @fftlog: a #NcmFftlog
 * @Fk: (scope call): a #NcmFftlogFunc
 * @user_data: @Fk user data
 * @reltol: relative tolerance
 *
 * Increases the original (input) number of knots until the $G(r)$ splines reach
 * the required precision @reltol.
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
  return self->Nf;
}

/**
 * ncm_fftlog_get_norma:
 * @fftlog: a #NcmFftlog
 *
 * Gets the number of knots $N_f^\prime$ where the integrated function is evaluated
 * plus padding.
 *
 * Returns: the total number of knots $N_f^\prime$ (double).
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
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
  NcmFftlogPrivate * const self = fftlog->priv;
  
  return self->Lk + 2.0 * self->Lk_N * self->pad;
}

/**
 * ncm_fftlog_get_mode_index:
 * @fftlog: a #NcmFftlog
 * @i: index
 *
 * Gets the index of the mode @i of the Fourier decomposition. This index corresponds
 * to the label $n$ in Eq. \eqref{eq:Gr_decomp}.
 *
 * Returns: the index of the mode
 */
gint
ncm_fftlog_get_mode_index (NcmFftlog *fftlog, gint i)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  
  return (i > self->Nf_2) ? i - self->Nf : i;
}

/**
 * ncm_fftlog_get_array_index:
 * @fftlog: a #NcmFftlog
 * @phys_i: index
 *
 * Gets the array index @i of the Fourier decomposition. This index corresponds the position
 * in the fft array of the element $n$ in Eq. \eqref{eq:Gr_decomp}.
 *
 * Returns: the array index corresponding to @phys_i
 */
gint
ncm_fftlog_get_array_index (NcmFftlog *fftlog, gint phys_i)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  
  return (phys_i < 0) ? phys_i + self->Nf : phys_i;
}

/**
 * ncm_fftlog_peek_output_vector:
 * @fftlog: a #NcmFftlog
 * @nderiv: derivative number
 *
 * Peeks the output vector respective to $G(r)$, @nderiv = 0, or
 * its @comp-th derivative with respect to $\ln r$.
 *
 * Returns: (transfer none): the output vector $G(r)$ or its @comp-th derivative.
 */
NcmVector *
ncm_fftlog_peek_output_vector (NcmFftlog *fftlog, guint nderiv)
{
  NcmFftlogPrivate * const self = fftlog->priv;
  
  return g_ptr_array_index (self->Gr_vec, nderiv);
}

