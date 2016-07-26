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
 * This class provides the tools to compute the Fast Fourier Transform of any function, which is assumed to be a periodic 
 * sequence of logarithmically spaced points. It is inspired on the approach [FFTLog][XHamilton2000], which we extended as described below.
 * 
 * A function $G(r)$ is written as 
 * \begin{equation}\label{eq:Gr} G(r) = \int_0^\infty F(k) \ K(kr) dk, \end{equation}
 * where $F(k)$ is defined in the fundamental interval $[\ln k_0 - L/2, \ln k_0 + L/2]$, $L$ is the period, 
 * $\ln k_0$ is the center value and $K(kr)$ is a kernel function. Assuming that $F(k)$ can be written in terms of the 
 * $N$ lowest Fourier modes, we have 
 *
 * $$F(k) = \sum_{n} c_n e^{\frac{2\pi i n}{L} \ln\left(\frac{k}{k_0}\right)}.$$
 * Substituting $F(k)$ in Eq. \eqref{eq:Gr} and changing the variable $k \rightarrow t = kr$, thus 
 * \begin{eqnarray}\label{eq:Gr_decomp}
 * r G(r) &=& \sum_n c_n \int_0^\infty \frac{k}{k_0}^{\frac{2\pi i n}{L}} K(kr)^2 d(kr) \\
 * &=& \sum_n c_n  \int_0^\infty \frac{t}{k_0 r}^{\frac{2\pi i n}{L}} K(t) dt \\
 * &=& \sum_n c_n e^{-\frac{2\pi i n}{L} \ln\left(\frac{r}{r_0}\right)} e^{-\frac{2\pi i n}{L} \ln(k_0 r_0)} Y_n, 
 * \end{eqnarray}
 * where 
 * $$Y_n = \int_0^\infty t^{\frac{2\pi i n}{L}} K(t) dt,$$
 * and the Fourier coefficients are
 * $$c_n = \frac{1}{N} \sum_m F(k_m) e^{- \frac{2\pi i nm}{N}}.$$
 * The total number of points $N$ corresponds to the number of knots in the fundamental interval, which is equally spaced.    
 * 
 * The user must provide the following input values: $\ln k_0$ - ncm_fftlog_set_lnk0(), $\ln r_0$ - ncm_fftlog_set_lnr0(), 
 * $L$ - ncm_fftlog_set_length(), padding percentage - ncm_fftlog_set_padding(), $N$ - ncm_fftlog_set_size(), 
 * $F(k)$ (or $F(k_m)$ -- see description below). 
 * 
 * - Since the algorithm assumes that the function to be decomposed is periodic, it is worth extending the interval in $\ln k$ such that 
 * $F(k) \equiv 0$ in the intervals $\left[\ln k_0 -\frac{L_T}{2}, \ln k_0 - \frac{L}{2} \right)$ and 
 * $ \left(\ln k_0 + \frac{L}{2}, \ln k_0 + \frac{L_T}{2}\right]$, where the total period $L_T$ is defined by the final 
 * number of knots, i.e., $N_f = N (1 + \mathrm{padding})$. 
 * - $N$ knots are equally distributed in the fundamental interval and $N \times \mathrm{padding}$ knots are distributed in 
 * in the two simetric intervals as mentioned above. 
 * - For the sake of optimization, the final number of points $N_f$ is substituted by the smallest number $N_f^\prime$ (bigger than $N_f$) 
 * which can be decomposed as $N_f \leq N_f^\prime = N^\prime (1 + \mathrm{padding}) = 3^a 5^b 7^c$, where $a$, 
 * $b$ and $c$ are positive integers. 
 * - The function $F(k)$ can be provided as:
 * 1. a gsl_function - ncm_fftlog_eval_by_function() - whose values are computed at the knots 
 * within the fundamental interval, and set to zero within the padding intervals. 
 * 2. as a vector - ncm_fftlog_eval_by_vector() - first one must get the vector of $\ln k$ knots, ncm_fftlog_get_vector_lnr(), 
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
  PROP_NDERIV,
  PROP_LNR0,
  PROP_LNK0,
  PROP_LR,
  PROP_N,
  PROP_PAD,
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
  fftlog->Gr_vec  = g_ptr_array_new ();
  fftlog->Gr_s    = g_ptr_array_new ();

  g_ptr_array_set_free_func (fftlog->Gr_vec, (GDestroyNotify)ncm_vector_free);
  g_ptr_array_set_free_func (fftlog->Gr_s, (GDestroyNotify)ncm_spline_free);
    
#ifdef NUMCOSMO_HAVE_FFTW3
  fftlog->Fk        = NULL;
  fftlog->Cm        = NULL;
  fftlog->Gr        = NULL;
  fftlog->Ym        = g_ptr_array_new ();
  fftlog->CmYm      = NULL;
  fftlog->p_Fk2Cm   = NULL;
  fftlog->p_CmYm2Gr = NULL;
  g_ptr_array_set_free_func (fftlog->Ym, (GDestroyNotify)fftw_free);
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
    case PROP_NAME:
      g_value_set_string (value, NCM_FFTLOG_GET_CLASS (fftlog)->name);
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
  g_clear_pointer (&fftlog->Fk, fftw_free);
  g_clear_pointer (&fftlog->Cm, fftw_free);
  g_clear_pointer (&fftlog->CmYm, fftw_free);
  g_clear_pointer (&fftlog->Gr, fftw_free);
  
  g_clear_pointer (&fftlog->p_Fk2Cm, fftw_destroy_plan);
  g_clear_pointer (&fftlog->p_CmYm2Gr, fftw_destroy_plan);

  ncm_vector_clear (&fftlog->lnr_vec);

  g_ptr_array_set_size (fftlog->Gr_vec, 0);
  g_ptr_array_set_size (fftlog->Gr_s, 0);
  g_ptr_array_set_size (fftlog->Ym, 0);  
}
#endif /* NUMCOSMO_HAVE_FFTW3 */

static void
_ncm_fftlog_finalize (GObject *object)
{
#ifdef NUMCOSMO_HAVE_FFTW3
  NcmFftlog *fftlog = NCM_FFTLOG (object);
  _ncm_fftlog_free_all (fftlog);

  g_clear_pointer (&fftlog->Gr_vec, (GDestroyNotify)g_ptr_array_unref);
  g_clear_pointer (&fftlog->Gr_s, (GDestroyNotify)g_ptr_array_unref);
  g_clear_pointer (&fftlog->Ym, (GDestroyNotify)g_ptr_array_unref);
  
#endif /* NUMCOSMO_HAVE_FFTW3 */
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fftlog_parent_class)->finalize (object);
}

static void
ncm_fftlog_class_init (NcmFftlogClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_fftlog_set_property;
  object_class->get_property = &_ncm_fftlog_get_property;
  object_class->finalize     = &_ncm_fftlog_finalize;

  g_object_class_install_property (object_class,
                                   PROP_NDERIV,
                                   g_param_spec_uint ("nderivs",
                                                      NULL,
                                                      "Number of derivatives",
                                                      0, G_MAXUINT32, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
 * 
 * Returns: (transfer none): The internal string describing #NcmFftlog.
 */
const gchar *
ncm_fftlog_peek_name (NcmFftlog *fftlog)
{
  return NCM_FFTLOG_GET_CLASS (fftlog)->name;
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
  if (fftlog->nderivs != nderivs)
  {
    if (nderivs < fftlog->nderivs)
    {
      g_ptr_array_set_size (fftlog->Gr_vec, nderivs + 1);
      g_ptr_array_set_size (fftlog->Gr_s, nderivs + 1);
      g_ptr_array_set_size (fftlog->Ym, nderivs + 1);
    }
    else
    {
      if (fftlog->N != 0)
      {
        guint i;
        for (i = fftlog->nderivs + 1; i <= nderivs; i++)
        {
          NcmVector *Gr_vec_i = ncm_vector_new (fftlog->N);
          NcmSpline *Gr_s_i   = ncm_spline_cubic_notaknot_new_full (fftlog->lnr_vec, Gr_vec_i, FALSE);
          fftw_complex *Ym_i  = fftw_alloc_complex (fftlog->Nf);

          g_ptr_array_add (fftlog->Gr_vec, Gr_vec_i);
          g_ptr_array_add (fftlog->Gr_s, Gr_s_i);
          g_ptr_array_add (fftlog->Ym, Ym_i);
        }
      }

      fftlog->prepared  = FALSE;
      fftlog->evaluated = FALSE;
    }

    fftlog->nderivs = nderivs;
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
  return fftlog->nderivs;
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
 * Sets the number of knots $N_f^\prime$ where the integrated function is evaluated, 
 * given the input number of knots @n, plus padding.  
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

    fftlog->Fk      = fftw_alloc_complex (fftlog->Nf);
    fftlog->Cm      = fftw_alloc_complex (fftlog->Nf);
    fftlog->CmYm    = fftw_alloc_complex (fftlog->Nf);
    fftlog->Gr      = fftw_alloc_complex (fftlog->Nf);

    fftlog->lnr_vec = ncm_vector_new (fftlog->N);

    ncm_cfg_load_fftw_wisdom ("ncm_fftlog_%s.fftw", NCM_FFTLOG_GET_CLASS (fftlog)->name);
    
    fftlog->p_Fk2Cm   = fftw_plan_dft_1d (fftlog->Nf, fftlog->Fk, fftlog->Cm, FFTW_FORWARD, fftw_default_flags | FFTW_DESTROY_INPUT);
    fftlog->p_CmYm2Gr = fftw_plan_dft_1d (fftlog->Nf, fftlog->CmYm, fftlog->Gr, FFTW_FORWARD, fftw_default_flags | FFTW_DESTROY_INPUT);

    for (i = 0; i <= fftlog->nderivs; i++)
    {
      NcmVector *Gr_vec_i = ncm_vector_new (fftlog->N);
      NcmSpline *Gr_s_i   = ncm_spline_cubic_notaknot_new_full (fftlog->lnr_vec, Gr_vec_i, FALSE);
      fftw_complex *Ym_i  = fftw_alloc_complex (fftlog->Nf);

      g_ptr_array_add (fftlog->Gr_vec, Gr_vec_i);
      g_ptr_array_add (fftlog->Gr_s, Gr_s_i);
      g_ptr_array_add (fftlog->Ym, Ym_i);
    }

    ncm_cfg_save_fftw_wisdom ("ncm_fftlog_%s.fftw", NCM_FFTLOG_GET_CLASS (fftlog)->name);

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
 * Sets the length of the period @Lk, where the function is periodic in logarithmic space $\ln k$.  
 * 
 */
void
ncm_fftlog_set_length (NcmFftlog *fftlog, gdouble Lk)
{
  if (fftlog->Lk != Lk)
  {
    fftlog->Lk        = Lk;
    fftlog->Lk_N      = fftlog->Lk / (1.0 * fftlog->N - 1.0);
    fftlog->prepared  = FALSE;
    fftlog->evaluated = FALSE;
  }
}

#ifdef NUMCOSMO_HAVE_FFTW3
static void
_ncm_fftlog_eval (NcmFftlog *fftlog)
{
  guint nd;
  gint i;

  fftw_execute (fftlog->p_Fk2Cm);
  
  if (!fftlog->prepared)
  {
    const gdouble twopi_Lt = 2.0 * M_PI / ncm_fftlog_get_full_length (fftlog);
    fftw_complex *Ym_0     = g_ptr_array_index (fftlog->Ym, 0);
    fftw_complex *Ym_ndm1;
    
    NCM_FFTLOG_GET_CLASS (fftlog)->get_Ym (fftlog, Ym_0);

    for (i = 0; i < fftlog->Nf; i++)
    {
      const gint phys_i      = ncm_fftlog_get_mode_index (fftlog, i);
      const complex double a = twopi_Lt * phys_i * I;
      
      Ym_0[i] *= cexp (- a * (fftlog->lnk0 + fftlog->lnr0));

      Ym_ndm1 = Ym_0;
      for (nd = 1; nd <= fftlog->nderivs; nd++)
      {
        fftw_complex *Ym_nd = g_ptr_array_index (fftlog->Ym, nd);
        Ym_nd[i] = -(1.0 + a) * Ym_ndm1[i];
        Ym_ndm1  = Ym_nd;
      }
    }
    fftlog->prepared = TRUE;
  }

  for (i = 0; i < fftlog->N; i++)
  {
    const gint phys_i      = i - fftlog->N_2;
    const gdouble lnr      = fftlog->lnr0 + phys_i * fftlog->Lk_N;
    
    ncm_vector_set (fftlog->lnr_vec, i, lnr);
  }
  
  for (nd = 0; nd <= fftlog->nderivs; nd++)
  {
    const gdouble norma = ncm_fftlog_get_norma (fftlog);
    NcmVector *Gr_nd    = g_ptr_array_index (fftlog->Gr_vec, nd);
    fftw_complex *Ym_nd = g_ptr_array_index (fftlog->Ym, nd);
    
    for (i = 0; i < fftlog->Nf; i++)
    {
      fftlog->CmYm[i] = fftlog->Cm[i] * Ym_nd[i]; 
    }
/*    
    printf ("% 20.15g % 20.15g | % 20.15g % 20.15g\n", 
            creal (fftlog->CmYm[fftlog->Nf_2]),
            cimag (fftlog->CmYm[fftlog->Nf_2]),
            creal (fftlog->CmYm[fftlog->Nf_2 + 1]),
            cimag (fftlog->CmYm[fftlog->Nf_2 + 1])
            );
*/
    fftlog->CmYm[fftlog->Nf_2]     = creal (fftlog->CmYm[fftlog->Nf_2]);
    fftlog->CmYm[fftlog->Nf_2 + 1] = creal (fftlog->CmYm[fftlog->Nf_2 + 1]);

    fftw_execute (fftlog->p_CmYm2Gr);

    for (i = 0; i < fftlog->N; i++)
    {
      const gdouble lnr      = ncm_vector_get (fftlog->lnr_vec, i);
      const gdouble rm1      = exp (-lnr);
      const gdouble Gr_nd_i  = creal (fftlog->Gr[i + fftlog->pad]) * rm1 / norma;

      ncm_vector_set (Gr_nd, i, Gr_nd_i);
    }
  }

  fftlog->evaluated = TRUE;
}
#endif /* NUMCOSMO_HAVE_FFTW3 */

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
 * Evaluates the function @Fk at each knot $\ln k_m$.
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
 * Prepares the set of splines respective to the function $G(r)$ 
 * and, if required, its n-order derivatives. 
 * 
 */
void 
ncm_fftlog_prepare_splines (NcmFftlog *fftlog)
{
  guint nd;

  g_assert (fftlog->evaluated);
  
  for (nd = 0; nd <= fftlog->nderivs; nd++)
    ncm_spline_prepare (g_ptr_array_index (fftlog->Gr_s, nd));
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
  return ncm_vector_ref (fftlog->lnr_vec);
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
  g_assert (fftlog->evaluated);
  return ncm_vector_ref (g_ptr_array_index (fftlog->Gr_vec, nderiv));
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
  g_assert (fftlog->evaluated);
  return g_ptr_array_index (fftlog->Gr_s, nderiv);
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
 * ncm_fftlog_calibrate_size: (skip)
 * @fftlog: a #NcmFftlog
 * @Fk: Fk function pointer
 * @reltol: relative tolerance
 * 
 * Increases the original (input) number of knots until the $G(r)$ splines reach 
 * the required precision @reltol.  
 * 
 */
void 
ncm_fftlog_calibrate_size (NcmFftlog *fftlog, gsl_function *Fk, gdouble reltol)
{
  NcmSpline **s = g_new0 (NcmSpline *, fftlog->nderivs + 1);
  guint nd, size;
  gdouble lreltol = 0.0;

  ncm_fftlog_eval_by_function (fftlog, Fk);
  ncm_fftlog_prepare_splines (fftlog);

  for (nd = 0; nd <= fftlog->nderivs; nd++)
  {
    s[nd] = ncm_spline_ref (g_ptr_array_index (fftlog->Gr_s, nd));
  }

  /*printf ("# Initial size %u [%u].\n", fftlog->N, fftlog->pad);*/
  ncm_fftlog_set_size (fftlog, fftlog->N * 1.2);
  /*printf ("# Trying size %u [%u].\n", fftlog->N, fftlog->pad);*/
  ncm_fftlog_eval_by_function (fftlog, Fk);
  ncm_fftlog_prepare_splines (fftlog);

  size = ncm_spline_get_len (g_ptr_array_index (fftlog->Gr_s, 0));
  for (nd = 0; nd <= fftlog->nderivs; nd++)
  {
    guint i;
    gdouble absmin, absmax;
    NcmVector *Gr_vec_nd = g_ptr_array_index (fftlog->Gr_vec, nd);
    ncm_vector_get_absminmax (Gr_vec_nd, &absmin, &absmax);
    
    /*printf ("# Testing component %u [% 20.15g, % 20.15g].\n", nd, absmin, absmax);*/
    for (i = 0; i < size; i++)
    {
      const gdouble lnr_i = ncm_vector_get (fftlog->lnr_vec, i);
      const gdouble lnG_i = ncm_vector_get (Gr_vec_nd, i);
      const gdouble lnS_i = ncm_spline_eval (s[nd], lnr_i);
      const gdouble lreltol_i = fabs ((lnG_i - lnS_i) / (fabs (lnG_i) + absmax));
      
      /*printf ("% 20.15g% 20.15e % 20.15g % 20.15g | % 20.15e\n", lnr_i, exp (lnr_i), lnG_i, lnS_i, lreltol_i);*/
      lreltol = GSL_MAX (lreltol_i, lreltol);
    }

    ncm_spline_clear (&s[nd]);
    /*printf ("# Largest error up to component %u is %e.\n", nd, lreltol);*/
    /*fflush (stdout);*/
  }

  g_clear_pointer (&s, g_free);

  if (lreltol > reltol)
    ncm_fftlog_calibrate_size (fftlog, Fk, reltol);
}

/**
 * ncm_fftlog_get_size:
 * @fftlog: a #NcmFftlog
 * 
 * Gets the number of knots $N^\prime$ where the integrated function is evaluated.
 * 
 * Returns: the number of knots $N^\prime$.
 */
/**
 * ncm_fftlog_get_full_size:
 * @fftlog: a #NcmFftlog
 * 
 * Gets the number of knots $N_f^\prime$ where the integrated function is evaluated
 * plus padding.
 * 
 * Returns: the total number of knots $N_f^\prime$.
 */
/**
 * ncm_fftlog_get_norma:
 * @fftlog: a #NcmFftlog
 * 
 * Gets the number of knots $N_f^\prime$ where the integrated function is evaluated
 * plus padding.
 * 
 * Returns: the total number of knots $N_f^\prime$ (double).
 */
/**
 * ncm_fftlog_get_length:
 * @fftlog: a #NcmFftlog
 * 
 * Gets the value of the ``physical'' period, i.e., period of the fundamental interval.
 * 
 * Returns: the period $L$.
 */
/**
 * ncm_fftlog_get_full_length:
 * @fftlog: a #NcmFftlog
 * 
 * Gets the value of the total period, i.e., period defined by the fundamental interval plus the padding size.
 * 
 * Returns: the total period $L_T$.
 */
/**
 * ncm_fftlog_get_mode_index:
 * @fftlog: a #NcmFftlog
 * @i: index
 * 
 * Gets the index of the mode @i of the Fourier decomposition. This index corresponds 
 * to the lable $n$ in Eq. \eqref{eq:Gr_decomp}.
 * 
 * Returns: the index of the mode
 */
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
