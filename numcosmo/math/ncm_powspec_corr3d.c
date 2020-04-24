/***************************************************************************
 *            ncm_powspec_corr3d.c
 *
 *  Wed March 20 15:22:51 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_powspec_corr3d.c
 * Copyright (C) 2019 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_powspec_corr3d
 * @title: NcmPowspecCorr3d
 * @short_description: Class to compute filtered power spectrum
 * @stability: Stable
 * @include: numcosmo/math/ncm_powspec_corr3d.h
 * 
 * This class computes the 3d spatial correlation function $\xi(r, z)$ from the power spectrum 
 * using the FFTLog approach (see #NcmFftlog and #NcmFftlogSBesselJ),
 * \begin{equation}\label{eq:variance}
 * \xi(r, z) = \sigma^2(r, z) = \frac{1}{2\pi^2} \int_0^\infty k^2 \ P(k, z) j_0(kr) \ \mathrm{d}k, 
 * \end{equation}
 * where $P(k, z)$ is the power spectrum at mode $k$ and redshift $z$ and $j_0(kr)$ is the order zero spherical bessel function.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_powspec_corr3d.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline2d_bicubic.h"
#include "math/ncm_fftlog_sbessel_j.h"
#include "math/ncm_c.h"
#include "ncm_enum_types.h"

enum
{
  PROP_0,
  PROP_LNR0,
  PROP_ZI,
  PROP_ZF,
  PROP_RELTOL,
  PROP_RELTOL_Z,
  PROP_POWERSPECTRUM,
	PROP_SIZE,
};

G_DEFINE_TYPE (NcmPowspecCorr3d, ncm_powspec_corr3d, G_TYPE_OBJECT);

static void
ncm_powspec_corr3d_init (NcmPowspecCorr3d *psc)
{
  psc->ps          = NULL;
  psc->lnr0        = 0.0;
  psc->lnk0        = 0.0;
  psc->Lk          = 0.0;
  psc->zi          = 0.0;
  psc->zf          = 0.0;
  psc->reltol      = 0.0;
  psc->fftlog      = NULL;
  psc->calibrated  = FALSE;
  psc->xi          = ncm_spline2d_bicubic_notaknot_new ();
  psc->ctrl        = ncm_model_ctrl_new (NULL);
  psc->constructed = FALSE;
}

static void
_ncm_powspec_corr3d_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmPowspecCorr3d *psc = NCM_POWSPEC_CORR3D (object);
  g_return_if_fail (NCM_IS_POWSPEC_CORR3D (object));

  switch (prop_id)
  {
    case PROP_LNR0:
      ncm_powspec_corr3d_set_lnr0 (psc, g_value_get_double (value));
      break;
    case PROP_ZI:
      ncm_powspec_corr3d_set_zi (psc, g_value_get_double (value));
      break;
    case PROP_ZF:
      ncm_powspec_corr3d_set_zf (psc, g_value_get_double (value));
      break;
    case PROP_RELTOL:
      psc->reltol = g_value_get_double (value);
      break;
    case PROP_RELTOL_Z:
      psc->reltol_z = g_value_get_double (value);
      break;
    case PROP_POWERSPECTRUM:
      psc->ps = g_value_dup_object (value);
      psc->zi = ncm_powspec_get_zi (psc->ps);
      psc->zf = ncm_powspec_get_zf (psc->ps);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_powspec_corr3d_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmPowspecCorr3d *psc = NCM_POWSPEC_CORR3D (object);
  g_return_if_fail (NCM_IS_POWSPEC_CORR3D (object));

  switch (prop_id)
  {
    case PROP_LNR0:
      g_value_set_double (value, psc->lnr0);
      break;
    case PROP_ZI:
      g_value_set_double (value, psc->zi);
      break;
    case PROP_ZF:
      g_value_set_double (value, psc->zf);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, psc->reltol);
      break;
    case PROP_RELTOL_Z:
      g_value_set_double (value, psc->reltol_z);
      break;
    case PROP_POWERSPECTRUM:
      g_value_set_object (value, psc->ps);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_powspec_corr3d_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_powspec_corr3d_parent_class)->constructed (object);
  {
    NcmPowspecCorr3d *psc = NCM_POWSPEC_CORR3D (object);
    const gdouble lnk_min = log (ncm_powspec_get_kmin (psc->ps));
    const gdouble lnk_max = log (ncm_powspec_get_kmax (psc->ps));

    psc->lnk0 = 0.5 * (lnk_max + lnk_min);
    psc->Lk   = (lnk_max - lnk_min);

    ncm_fftlog_clear (&psc->fftlog);
    psc->fftlog = NCM_FFTLOG (ncm_fftlog_sbessel_j_new (0, psc->lnr0, psc->lnk0, psc->Lk, 100));
    /*ncm_fftlog_sbessel_j_set_q (NCM_FFTLOG_SBESSEL_J (psc->fftlog), 0.5);*/
    
    ncm_fftlog_set_padding (psc->fftlog, 1.0);
    ncm_fftlog_set_nderivs (psc->fftlog, 0);

    ncm_powspec_corr3d_set_best_lnr0 (psc);
    
    ncm_model_ctrl_force_update (psc->ctrl);
    psc->calibrated = FALSE;
  }
}

static void
_ncm_powspec_corr3d_dispose (GObject *object)
{
  NcmPowspecCorr3d *psc = NCM_POWSPEC_CORR3D (object);

  ncm_powspec_clear (&psc->ps);
  ncm_fftlog_clear (&psc->fftlog);

  ncm_spline2d_clear (&psc->xi);

  ncm_model_ctrl_clear (&psc->ctrl);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_corr3d_parent_class)->dispose (object);
}

static void
_ncm_powspec_corr3d_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_powspec_corr3d_parent_class)->finalize (object);
}

static void
ncm_powspec_corr3d_class_init (NcmPowspecCorr3dClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_powspec_corr3d_set_property;
  object_class->get_property = &_ncm_powspec_corr3d_get_property;
  object_class->constructed  = &_ncm_powspec_corr3d_constructed;
  object_class->dispose      = &_ncm_powspec_corr3d_dispose;
  object_class->finalize     = &_ncm_powspec_corr3d_finalize;

  /**
   * NcmPowspecCorr3d:lnr0:
   *
   * The Center value for $\ln(r)$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_LNR0,
                                   g_param_spec_double ("lnr0",
                                                        NULL,
                                                        "Output center value",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcmPowspecCorr3d:zi:
   *
   * The output inital time $z_i$ of the correlation function, $\xi(r,z)$. 
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ZI,
                                   g_param_spec_double ("zi",
                                                        NULL,
                                                        "Output initial time",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcmPowspecCorr3d:zf:
   *
   * The output final time $z_f$ of the correlation function, $\xi(r,z)$. 
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ZF,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "Output final time",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcmPowspecCorr3d:reltol:
   *
   * The relative tolerance for calibration in the distance direction. 
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance for calibration",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcmPowspecCorr3d:reltol-z:
   *
   * The relative tolerance for calibration in the redshift direction. 
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RELTOL_Z,
                                   g_param_spec_double ("reltol-z",
                                                        NULL,
                                                        "Relative tolerance for calibration in the redshift direction",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-6,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcmPowspecCorr3d:powerspectrum:
   *
   * The #NcmPowspec object to be transformed into $\xi(r, z)$. 
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_POWERSPECTRUM,
                                   g_param_spec_object ("powerspectrum",
                                                        NULL,
                                                        "NcmPowspec object",
                                                        NCM_TYPE_POWSPEC,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_powspec_corr3d_new:
 * @ps: a #NcmPowspec
 * 
 * Creates a new #NcmPowspecCorr3d from the power spectrum @ps. 
 * 
 * Returns: (transfer full): the newly created #NcmPowspecCorr3d.
 */
NcmPowspecCorr3d *
ncm_powspec_corr3d_new (NcmPowspec *ps)
{
  NcmPowspecCorr3d *psc = g_object_new (NCM_TYPE_POWSPEC_CORR3D,
                                        "powerspectrum", ps,
                                        NULL);
  
  return psc;
}

/**
 * ncm_powspec_corr3d_ref:
 * @psc: a #NcmPowspecCorr3d
 * 
 * Increases the reference count of @psc by one.
 * 
 * Returns: (transfer full): @psc
 */
NcmPowspecCorr3d *
ncm_powspec_corr3d_ref (NcmPowspecCorr3d *psc)
{
  return g_object_ref (psc);
}

/**
 * ncm_powspec_corr3d_free:
 * @psc: a #NcmPowspecCorr3d
 * 
 * Decreases the reference count of @psc by one.
 * 
 */
void
ncm_powspec_corr3d_free (NcmPowspecCorr3d *psc)
{
  g_object_unref (psc);
}

/**
 * ncm_powspec_corr3d_clear:
 * @psc: a #NcmPowspecCorr3d
 * 
 * If @psc is different from NULL, decreases the reference count of 
 * @psc by one and sets @fftlog to NULL.
 * 
 */
void
ncm_powspec_corr3d_clear (NcmPowspecCorr3d **psc)
{
  g_clear_object (psc);
}

typedef struct _NcmPowspecCorr3dArg
{
  NcmPowspecCorr3d *psc;
  NcmModel *model;
  gdouble z;
} NcmPowspecCorr3dArg;

static gdouble 
_ncm_powspec_corr3d_k2Pk (gdouble k, gpointer userdata)
{
  NcmPowspecCorr3dArg *arg = (NcmPowspecCorr3dArg *) userdata;
  const gdouble k2 = k * k;
  const gdouble Pk = ncm_powspec_eval (arg->psc->ps, arg->model, arg->z, k);
  const gdouble f  = Pk * k2 / ncm_c_2_pi_2 ();
  
  return f;
}

static gdouble 
_ncm_powspec_corr3d_dummy_z (gdouble z, gpointer userdata)
{
  NcmPowspecCorr3dArg *arg = (NcmPowspecCorr3dArg *) userdata;
  gsl_function F;
  
  F.function = &_ncm_powspec_corr3d_k2Pk;
  F.params   = arg;

  arg->z = z;
  ncm_fftlog_eval_by_gsl_function (arg->psc->fftlog, &F);
  /*printf ("# z-knots % 20.15g % 20.15g\n", z, ncm_vector_get (ncm_fftlog_peek_output_vector (arg->psc->fftlog, 0), 0));*/
  return ncm_vector_get (ncm_fftlog_peek_output_vector (arg->psc->fftlog, 0), 0);
}

/**
 * ncm_powspec_corr3d_prepare:
 * @psc: a #NcmPowspecCorr3d
 * @model: a #NcmModel
 * 
 * Prepares the object applying the filter to the power spectrum.
 * 
 */
void
ncm_powspec_corr3d_prepare (NcmPowspecCorr3d *psc, NcmModel *model)
{
  NcmPowspecCorr3dArg arg;
  gsl_function F;

  F.function = &_ncm_powspec_corr3d_k2Pk;
  F.params   = &arg;

  arg.psc   = psc;
  arg.model = model;
  arg.z     = 0.0;

  ncm_powspec_prepare_if_needed (psc->ps, model);
  
  {
    const gdouble lnk_min = log (ncm_powspec_get_kmin (psc->ps));
    const gdouble lnk_max = log (ncm_powspec_get_kmax (psc->ps));

    psc->lnk0 = 0.5 * (lnk_max + lnk_min);
    psc->Lk   = (lnk_max - lnk_min);

    if ((psc->lnk0 != ncm_fftlog_get_lnk0 (psc->fftlog)) || (psc->Lk != ncm_fftlog_get_length (psc->fftlog)))
    {
      ncm_fftlog_set_lnk0 (psc->fftlog, psc->lnk0);
      ncm_fftlog_set_length (psc->fftlog, psc->Lk);
      psc->calibrated = FALSE;
    }
  }

  if (!psc->calibrated)
  {
    NcmMatrix *xi;
    NcmVector *z_vec, *lnr_vec;
    guint N_k = 0, N_z = 0;
    guint i;

    ncm_powspec_get_nknots (psc->ps, &N_z, &N_k);
    
    ncm_fftlog_calibrate_size_gsl (psc->fftlog, &F, psc->reltol);
    N_k = ncm_fftlog_get_size (psc->fftlog);

    {
      NcmSpline *dummy_z = ncm_spline_cubic_notaknot_new ();
      gsl_function Fdummy_z;

      Fdummy_z.function = &_ncm_powspec_corr3d_dummy_z;
      Fdummy_z.params   = &arg;

      ncm_spline_set_func (dummy_z, NCM_SPLINE_FUNCTION_SPLINE, &Fdummy_z, psc->zi, psc->zf, 0, psc->reltol_z);

      z_vec = ncm_spline_get_xv (dummy_z);
      N_z = ncm_vector_len (z_vec);

      ncm_spline_clear (&dummy_z);
    }

    g_assert_cmpuint (N_z, >, 0);
    g_assert_cmpuint (N_k, >, 0);
/*    
    printf ("# Calibrating in zmin % 20.15g zmax % 20.15g, rmin % 20.15g rmax % 20.15g, N_z = %u, N_k = %u\n",
            psc->zi, psc->zf, 
            ncm_powspec_corr3d_get_r_min (psc), 
            ncm_powspec_corr3d_get_r_max (psc), 
            N_z, N_k);
*/    
    xi      = ncm_matrix_new (N_z, N_k);
    lnr_vec = ncm_fftlog_get_vector_lnr (psc->fftlog);
        
    for (i = 0; i < N_z; i++)
    {
      NcmVector *xi_z  = ncm_matrix_get_row (xi, i);

      arg.z = ncm_vector_get (z_vec, i);
      ncm_fftlog_eval_by_gsl_function (psc->fftlog, &F);

      ncm_vector_memcpy (xi_z, ncm_fftlog_peek_output_vector (psc->fftlog, 0));

      ncm_vector_free (xi_z);
    }

    ncm_spline2d_set (psc->xi, lnr_vec, z_vec, xi, TRUE);

    ncm_vector_free (z_vec);
    ncm_vector_free (lnr_vec);
    ncm_matrix_free (xi);

    psc->calibrated = TRUE;
  }
  else
  {
    NcmMatrix *xi = psc->xi->zm;
    
    guint N_z = ncm_matrix_nrows (xi);
    guint i;

    for (i = 0; i < N_z; i++)
    {
      NcmVector *xi_z  = ncm_matrix_get_row (xi, i);

      arg.z = ncm_vector_get (psc->xi->yv, i);
      ncm_fftlog_eval_by_gsl_function (psc->fftlog, &F);

      ncm_vector_memcpy (xi_z, ncm_fftlog_peek_output_vector (psc->fftlog, 0));
    }

    ncm_spline2d_prepare (psc->xi);
  }
}

/**
 * ncm_powspec_corr3d_prepare_if_needed:
 * @psc: a #NcmPowspecCorr3d
 * @model: a #NcmModel
 * 
 * Prepares (if necessary) the object applying the filter to the power spectrum.
 * 
 */
void
ncm_powspec_corr3d_prepare_if_needed (NcmPowspecCorr3d *psc, NcmModel *model)
{
  gboolean model_up = ncm_model_ctrl_update (psc->ctrl, model);

  if (model_up)
    ncm_powspec_corr3d_prepare (psc, model);
}

/**
 * ncm_powspec_corr3d_set_lnr0:
 * @psc: a #NcmPowspecCorr3d
 * @lnr0: the output center value $\ln(r_0)$
 * 
 * Sets the center of the transform output $\ln(r_0)$ (see ncm_fftlog_set_lnr0()). 
 * 
 */
void 
ncm_powspec_corr3d_set_lnr0 (NcmPowspecCorr3d *psc, gdouble lnr0)
{
  if (psc->lnr0 != lnr0)
  {
    const gdouble lnk_min = log (ncm_powspec_get_kmin (psc->ps));
    const gdouble lnk_max = log (ncm_powspec_get_kmax (psc->ps));

    psc->lnk0 = 0.5 * (lnk_max + lnk_min);
    psc->Lk   = (lnk_max - lnk_min);

    ncm_fftlog_set_lnk0 (psc->fftlog, psc->lnk0);
    ncm_fftlog_set_length (psc->fftlog, psc->Lk);
    
    psc->lnr0 = lnr0;
    ncm_model_ctrl_force_update (psc->ctrl);
    psc->calibrated = FALSE;
    ncm_fftlog_set_lnr0 (psc->fftlog, psc->lnr0);

    if (psc->lnr0 < -psc->lnk0)
      g_warning ("ncm_powspec_corr3d_set_lnr0: the requested center of the output does not satisfy r0k0 > 1.");
  }
}

/**
 * ncm_powspec_corr3d_set_best_lnr0:
 * @psc: a #NcmPowspecCorr3d
 * 
 * Sets the value of $\ln(r_0)$ which gives the best results for
 * the transformation based on the current value of $\ln(k_0)$,
 * this is based in the rule of thumb $\mathrm{max}_{x^*}(j_l)$
 * where $ x^* \approx l + 1$ (see ncm_fftlog_sbessel_j_set_best_lnr0()).
 * 
 */
void 
ncm_powspec_corr3d_set_best_lnr0 (NcmPowspecCorr3d *psc)
{
  const gdouble lnk_min = log (ncm_powspec_get_kmin (psc->ps));
  const gdouble lnk_max = log (ncm_powspec_get_kmax (psc->ps));

  psc->lnk0 = 0.5 * (lnk_max + lnk_min);
  psc->Lk   = (lnk_max - lnk_min);

  ncm_fftlog_set_lnk0 (psc->fftlog, psc->lnk0);
  ncm_fftlog_set_length (psc->fftlog, psc->Lk);
  
  ncm_powspec_corr3d_set_lnr0 (psc, -psc->lnk0);
}

/**
 * ncm_powspec_corr3d_set_zi:
 * @psc: a #NcmPowspecCorr3d
 * @zi: the output initial time $z_i$
 * 
 * Sets the inital time $z_i$. 
 * 
 */
void 
ncm_powspec_corr3d_set_zi (NcmPowspecCorr3d *psc, gdouble zi)
{
  if (psc->zi != zi)
  {
    psc->zi = zi;
    ncm_model_ctrl_force_update (psc->ctrl);
    psc->calibrated = FALSE;
    ncm_powspec_require_zi (psc->ps, zi);
  }
}

/**
 * ncm_powspec_corr3d_set_zf:
 * @psc: a #NcmPowspecCorr3d
 * @zf: the output final time $z_f$
 * 
 * Sets the final time $z_f$. 
 * 
 */
void 
ncm_powspec_corr3d_set_zf (NcmPowspecCorr3d *psc, gdouble zf)
{
  if (psc->zf != zf)
  {
    psc->zf = zf;
    ncm_model_ctrl_force_update (psc->ctrl);
    psc->calibrated = FALSE;
    ncm_powspec_require_zf (psc->ps, zf);
  }
}

/**
 * ncm_powspec_corr3d_get_r_min:
 * @psc: a #NcmPowspecCorr3d
 * 
 * This function returns $\xi(r, z)$'s minimum evaluated distance. 
 * 
 * Returns: the minimum distance $r_{\mathrm{min}}$.  
 */
gdouble
ncm_powspec_corr3d_get_r_min (NcmPowspecCorr3d *psc)
{
  return exp (psc->lnr0 - psc->Lk * 0.5);
}

/**
 * ncm_powspec_corr3d_get_r_max:
 * @psc: a #NcmPowspecCorr3d
 * 
 * This function returns $\xi(r, z)$'s maximum evaluated distance. 
 * 
 * Returns: the maximum distance $r_{\mathrm{max}}$.  
 */
gdouble
ncm_powspec_corr3d_get_r_max (NcmPowspecCorr3d *psc)
{
  return exp (psc->lnr0 + psc->Lk * 0.5);
}

/**
 * ncm_powspec_corr3d_eval_xi_lnr:
 * @psc: a #NcmPowspecCorr3d
 * @z: redshift
 * @lnr: logarithm base e of $r$
 * 
 * Evaluates the function $\xi(z, r)$ at @lnr and @z.
 * 
 * Returns: $\xi(\ln r, z)$. 
 */
gdouble
ncm_powspec_corr3d_eval_xi_lnr (NcmPowspecCorr3d *psc, const gdouble z, const gdouble lnr)
{
  return ncm_spline2d_eval (psc->xi, lnr, z);
}

/**
 * ncm_powspec_corr3d_eval_xi:
 * @psc: a #NcmPowspecCorr3d
 * @z: redshift $z$
 * @r: distance $r$ 
 * 
 * Evaluate the function $\xi(z, r)$ at @r and @z.
 * 
 * Returns: $\xi(r, z)$. 
 */
gdouble
ncm_powspec_corr3d_eval_xi (NcmPowspecCorr3d *psc, const gdouble z, const gdouble r)
{
  return ncm_powspec_corr3d_eval_xi_lnr (psc, z, log (r));
}
