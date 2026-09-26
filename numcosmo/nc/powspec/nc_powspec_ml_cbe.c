/***************************************************************************
 *            nc_powspec_ml_cbe.c
 *
 *  Tue April 05 10:42:14 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_powspec_ml_cbe.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcPowspecMLCBE:
 *
 * Linear matter power spectrum from CLASS.
 *
 * Computes the linear matter power spectrum with the
 * [CLASS](https://lesgourg.github.io/class_public/class.html) backend #NcCBE, up to
 * #NcPowspecMLCBE:intern-k-max. Outside the range CLASS computed, $P$ is the
 * Eisenstein-Hu spectrum, see #NcTransferFuncEH, times the ratio $r = P_\mathrm{CLASS} /
 * P_\mathrm{EH}$ continued as a power law in $k$ from the nearest edge $k_e$,
 * $$P(k, z) = r(k_e, z) \left(\frac{k}{k_e}\right)^{\beta(z)} P_\mathrm{EH}(k, z),$$
 * where $\beta$ is the mean slope of $\ln r$ over the decade of computed modes next to
 * $k_e$, or over the whole computed range when it spans less than a decade.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/powspec/nc_powspec_ml_cbe.h"
#include "nc/powspec/nc_powspec_ml_transfer.h"
#include "nc/powspec/nc_transfer_func_eh.h"
#include "nc/primordial/nc_hiprim.h"

enum
{
  PROP_0,
  PROP_CBE,
  PROP_CBE_K_MIN,
  PROP_CBE_K_MAX,
  PROP_SIZE
};

typedef struct _NcPowspecMLCBEPrivate
{
  NcCBE *cbe;
  NcmSpline2d *lnPk;
  NcPowspecML *eh;
  gdouble intern_k_min;
  gdouble calc_k_min;
  gdouble calc_k_max;
  gdouble intern_k_max;
} NcPowspecMLCBEPrivate;


struct _NcPowspecMLCBE
{
  NcPowspecML parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcPowspecMLCBE, nc_powspec_ml_cbe, NC_TYPE_POWSPEC_ML)

static void
nc_powspec_ml_cbe_init (NcPowspecMLCBE *ps_cbe)
{
  NcPowspecMLCBEPrivate * const self = nc_powspec_ml_cbe_get_instance_private (ps_cbe);
  NcTransferFunc *tf                 = nc_transfer_func_eh_new ();

  self->cbe          = NULL;
  self->lnPk         = NULL;
  self->eh           = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  self->intern_k_min = 0.0;
  self->calc_k_min   = 0.0;
  self->calc_k_max   = 0.0;
  self->intern_k_max = 0.0;

  nc_transfer_func_free (tf);
}

static void
_nc_powspec_ml_cbe_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcPowspecMLCBE *ps_cbe = NC_POWSPEC_ML_CBE (object);

  g_return_if_fail (NC_IS_POWSPEC_ML_CBE (object));

  switch (prop_id)
  {
    case PROP_CBE:
      nc_powspec_ml_cbe_set_cbe (ps_cbe, g_value_get_object (value));
      break;
    case PROP_CBE_K_MIN:
      nc_powspec_ml_cbe_set_intern_k_min (ps_cbe, g_value_get_double (value));
      break;
    case PROP_CBE_K_MAX:
      nc_powspec_ml_cbe_set_intern_k_max (ps_cbe, g_value_get_double (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_powspec_ml_cbe_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcPowspecMLCBE *ps_cbe = NC_POWSPEC_ML_CBE (object);

  g_return_if_fail (NC_IS_POWSPEC_ML_CBE (object));

  switch (prop_id)
  {
    case PROP_CBE:
      g_value_set_object (value, nc_powspec_ml_cbe_peek_cbe (ps_cbe));
      break;
    case PROP_CBE_K_MIN:
      g_value_set_double (value, nc_powspec_ml_cbe_get_intern_k_min (ps_cbe));
      break;
    case PROP_CBE_K_MAX:
      g_value_set_double (value, nc_powspec_ml_cbe_get_intern_k_max (ps_cbe));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_powspec_ml_cbe_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_powspec_ml_cbe_parent_class)->constructed (object);
  {
    NcPowspecMLCBE *ps_cbe             = NC_POWSPEC_ML_CBE (object);
    NcPowspecMLCBEPrivate * const self = nc_powspec_ml_cbe_get_instance_private (ps_cbe);

    if (self->cbe == NULL)
      self->cbe = nc_cbe_new ();

    g_assert_cmpfloat (self->intern_k_min, <, self->intern_k_max);
  }
}

static void
_nc_powspec_ml_cbe_dispose (GObject *object)
{
  NcPowspecMLCBE *ps_cbe             = NC_POWSPEC_ML_CBE (object);
  NcPowspecMLCBEPrivate * const self = nc_powspec_ml_cbe_get_instance_private (ps_cbe);

  nc_cbe_clear (&self->cbe);
  ncm_spline2d_clear (&self->lnPk);
  nc_powspec_ml_clear (&self->eh);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_powspec_ml_cbe_parent_class)->dispose (object);
}

static void
_nc_powspec_ml_cbe_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_powspec_ml_cbe_parent_class)->finalize (object);
}

static void _nc_powspec_ml_cbe_prepare (NcmPowspec *powspec, NcmModel *model);
static gdouble _nc_powspec_ml_cbe_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k);
static gdouble _nc_powspec_ml_cbe_deriv_z (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k);
static void _nc_powspec_ml_cbe_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk);

static void
nc_powspec_ml_cbe_class_init (NcPowspecMLCBEClass *klass)
{
  GObjectClass *object_class     = G_OBJECT_CLASS (klass);
  NcmPowspecClass *powspec_class = NCM_POWSPEC_CLASS (klass);

  object_class->set_property = &_nc_powspec_ml_cbe_set_property;
  object_class->get_property = &_nc_powspec_ml_cbe_get_property;

  object_class->constructed = &_nc_powspec_ml_cbe_constructed;
  object_class->dispose     = &_nc_powspec_ml_cbe_dispose;
  object_class->finalize    = &_nc_powspec_ml_cbe_finalize;

  /**
   * NcPowspecMLCBE:cbe:
   *
   * Class backend object.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_CBE,
                                   g_param_spec_object ("cbe",
                                                        NULL,
                                                        "Class backend object",
                                                        NC_TYPE_CBE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcPowspecMLCBE:intern-k-min:
   *
   * The smallest mode $k$ requested from CLASS, in $\mathrm{Mpc}^{-1}$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_CBE_K_MIN,
                                   g_param_spec_double ("intern-k-min",
                                                        NULL,
                                                        "Class minimum mode k",
                                                        G_MINDOUBLE, G_MAXDOUBLE, NC_POWSPEC_ML_CBE_INTERN_KMIN,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcPowspecMLCBE:intern-k-max:
   *
   * The largest mode $k$ computed by CLASS, in $\mathrm{Mpc}^{-1}$.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_CBE_K_MAX,
                                   g_param_spec_double ("intern-k-max",
                                                        NULL,
                                                        "Class maximum mode k",
                                                        G_MINDOUBLE, G_MAXDOUBLE, NC_POWSPEC_ML_CBE_INTERN_KMAX,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  powspec_class->prepare    = &_nc_powspec_ml_cbe_prepare;
  powspec_class->eval       = &_nc_powspec_ml_cbe_eval;
  powspec_class->deriv_z    = &_nc_powspec_ml_cbe_deriv_z;
  powspec_class->get_nknots = &_nc_powspec_ml_cbe_get_nknots;
}

static void
_nc_powspec_ml_cbe_prepare (NcmPowspec *powspec, NcmModel *model)
{
  NcHICosmo *cosmo                   = NC_HICOSMO (model);
  NcPowspecMLCBE *ps_cbe             = NC_POWSPEC_ML_CBE (powspec);
  NcPowspecMLCBEPrivate * const self = nc_powspec_ml_cbe_get_instance_private (ps_cbe);

  g_assert (NC_IS_HICOSMO (model));
  g_assert (ncm_model_peek_submodel_by_mid (model, nc_hiprim_id ()) != NULL);

  nc_cbe_set_calc_transfer (self->cbe, TRUE);
  nc_cbe_set_max_matter_pk_z (self->cbe, ncm_powspec_get_zf (powspec));

  nc_cbe_set_max_matter_pk_k (self->cbe, self->intern_k_max);

  nc_cbe_prepare_if_needed (self->cbe, cosmo);

  ncm_spline2d_clear (&self->lnPk);

  self->lnPk = nc_cbe_get_matter_ps (self->cbe);

  /* The range CLASS computed, from its grid; intern_k_min/max are the range requested */
  {
    NcmVector *lnk_v = ncm_spline2d_peek_xv (self->lnPk);

    self->calc_k_min = exp (ncm_vector_get (lnk_v, 0));
    self->calc_k_max = exp (ncm_vector_get (lnk_v, ncm_vector_len (lnk_v) - 1));
  }

  ncm_powspec_prepare_if_needed (NCM_POWSPEC (self->eh), model);
}

/*
 * The edge k_e nearest to k and the point k_1 one decade inside the computed range from
 * it, clamped to the other edge, over which the slope of ln r = ln P_CLASS - ln P_EH is
 * measured.
 */
static void
_nc_powspec_ml_cbe_match_points (NcPowspecMLCBEPrivate * const self, const gdouble k, gdouble *lnk_e, gdouble *lnk_1)
{
  const gdouble lnk_min = log (self->calc_k_min);
  const gdouble lnk_max = log (self->calc_k_max);

  if (k < self->calc_k_min)
  {
    *lnk_e = lnk_min;
    *lnk_1 = GSL_MIN_DBL (lnk_min + M_LN10, lnk_max);
  }
  else
  {
    *lnk_e = lnk_max;
    *lnk_1 = GSL_MAX_DBL (lnk_max - M_LN10, lnk_min);
  }
}

static gdouble
_nc_powspec_ml_cbe_eval (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  NcPowspecMLCBE *ps_cbe             = NC_POWSPEC_ML_CBE (powspec);
  NcPowspecMLCBEPrivate * const self = nc_powspec_ml_cbe_get_instance_private (ps_cbe);

  if ((k < self->calc_k_min) || (k > self->calc_k_max))
  {
    NcmPowspec *eh = NCM_POWSPEC (self->eh);
    gdouble lnk_e, lnk_1;

    _nc_powspec_ml_cbe_match_points (self, k, &lnk_e, &lnk_1);
    {
      const gdouble lnr_e = ncm_spline2d_eval (self->lnPk, lnk_e, z) - log (ncm_powspec_eval (eh, model, z, exp (lnk_e)));
      const gdouble lnr_1 = ncm_spline2d_eval (self->lnPk, lnk_1, z) - log (ncm_powspec_eval (eh, model, z, exp (lnk_1)));
      const gdouble beta  = (lnr_e - lnr_1) / (lnk_e - lnk_1);

      return exp (lnr_e + beta * (log (k) - lnk_e)) * ncm_powspec_eval (eh, model, z, k);
    }
  }
  else
  {
    return exp (ncm_spline2d_eval (self->lnPk, log (k), z));
  }
}

static gdouble
_nc_powspec_ml_cbe_deriv_z (NcmPowspec *powspec, NcmModel *model, const gdouble z, const gdouble k)
{
  NcPowspecMLCBE *ps_cbe             = NC_POWSPEC_ML_CBE (powspec);
  NcPowspecMLCBEPrivate * const self = nc_powspec_ml_cbe_get_instance_private (ps_cbe);

  if ((k < self->calc_k_min) || (k > self->calc_k_max))
  {
    NcmPowspec *eh = NCM_POWSPEC (self->eh);
    gdouble lnk_e, lnk_1;

    _nc_powspec_ml_cbe_match_points (self, k, &lnk_e, &lnk_1);
    {
      const gdouble k_e       = exp (lnk_e);
      const gdouble k_1       = exp (lnk_1);
      const gdouble Peh_e     = ncm_powspec_eval (eh, model, z, k_e);
      const gdouble Peh_1     = ncm_powspec_eval (eh, model, z, k_1);
      const gdouble Peh       = ncm_powspec_eval (eh, model, z, k);
      const gdouble lnr_e     = ncm_spline2d_eval (self->lnPk, lnk_e, z) - log (Peh_e);
      const gdouble lnr_1     = ncm_spline2d_eval (self->lnPk, lnk_1, z) - log (Peh_1);
      const gdouble dlnr_e_dz = ncm_spline2d_deriv_dzdy (self->lnPk, lnk_e, z) - ncm_powspec_deriv_z (eh, model, z, k_e) / Peh_e;
      const gdouble dlnr_1_dz = ncm_spline2d_deriv_dzdy (self->lnPk, lnk_1, z) - ncm_powspec_deriv_z (eh, model, z, k_1) / Peh_1;
      const gdouble beta      = (lnr_e - lnr_1) / (lnk_e - lnk_1);
      const gdouble dbeta_dz  = (dlnr_e_dz - dlnr_1_dz) / (lnk_e - lnk_1);
      const gdouble dlnk      = log (k) - lnk_e;
      const gdouble Pk        = exp (lnr_e + beta * dlnk) * Peh;

      return Pk * (dlnr_e_dz + dbeta_dz * dlnk + ncm_powspec_deriv_z (eh, model, z, k) / Peh);
    }
  }
  else
  {
    return exp (ncm_spline2d_eval (self->lnPk, log (k), z)) * ncm_spline2d_deriv_dzdy (self->lnPk, log (k), z);
  }
}

static void
_nc_powspec_ml_cbe_get_nknots (NcmPowspec *powspec, guint *Nz, guint *Nk)
{
  NcPowspecMLCBE *ps_cbe             = NC_POWSPEC_ML_CBE (powspec);
  NcPowspecMLCBEPrivate * const self = nc_powspec_ml_cbe_get_instance_private (ps_cbe);

  Nz[0] = ncm_vector_len (ncm_spline2d_peek_yv (self->lnPk));
  Nk[0] = ncm_vector_len (ncm_spline2d_peek_xv (self->lnPk));
}

/**
 * nc_powspec_ml_cbe_new:
 *
 * Creates a new #NcPowspecMLCBE from a new #NcCBE.
 *
 * Returns: (transfer full): the newly created #NcPowspecMLCBE.
 */
NcPowspecMLCBE *
nc_powspec_ml_cbe_new (void)
{
  NcPowspecMLCBE *ps_cbe = g_object_new (NC_TYPE_POWSPEC_ML_CBE,
                                         NULL);

  return ps_cbe;
}

/**
 * nc_powspec_ml_cbe_new_full:
 * @cbe: a #NcCBE
 *
 * Creates a new #NcPowspecMLCBE from @cbe.
 *
 * Returns: (transfer full): the newly created #NcPowspecMLCBE.
 */
NcPowspecMLCBE *
nc_powspec_ml_cbe_new_full (NcCBE *cbe)
{
  NcPowspecMLCBE *ps_cbe = g_object_new (NC_TYPE_POWSPEC_ML_CBE,
                                         "cbe", cbe,
                                         NULL);

  return ps_cbe;
}

/**
 * nc_powspec_ml_cbe_set_cbe:
 * @ps_cbe: a #NcPowspecMLCBE
 * @cbe: a #NcCBE
 *
 * Sets the #NcCBE to @cbe.
 *
 */
void
nc_powspec_ml_cbe_set_cbe (NcPowspecMLCBE *ps_cbe, NcCBE *cbe)
{
  NcPowspecMLCBEPrivate * const self = nc_powspec_ml_cbe_get_instance_private (ps_cbe);

  g_clear_object (&self->cbe);
  self->cbe = nc_cbe_ref (cbe);
}

/**
 * nc_powspec_ml_cbe_peek_cbe:
 * @ps_cbe: a #NcPowspecMLCBE
 *
 * Peeks the #NcCBE inside @ps_cbe.
 *
 * Returns: (transfer none): the #NcCBE inside @ps_cbe.
 */
NcCBE *
nc_powspec_ml_cbe_peek_cbe (NcPowspecMLCBE *ps_cbe)
{
  NcPowspecMLCBEPrivate * const self = nc_powspec_ml_cbe_get_instance_private (ps_cbe);

  return self->cbe;
}

/**
 * nc_powspec_ml_cbe_set_intern_k_min:
 * @ps_cbe: a #NcPowspecMLCBE
 * @k_min: the smallest $k$ requested from CLASS
 *
 * Sets #NcPowspecMLCBE:intern-k-min.
 *
 */
void
nc_powspec_ml_cbe_set_intern_k_min (NcPowspecMLCBE *ps_cbe, const gdouble k_min)
{
  NcPowspecMLCBEPrivate * const self = nc_powspec_ml_cbe_get_instance_private (ps_cbe);

  self->intern_k_min = k_min;

  if (self->intern_k_max > 0.0)
    g_assert_cmpfloat (self->intern_k_min, <, self->intern_k_max);
}

/**
 * nc_powspec_ml_cbe_set_intern_k_max:
 * @ps_cbe: a #NcPowspecMLCBE
 * @k_max: the largest $k$ computed by CLASS
 *
 * Sets #NcPowspecMLCBE:intern-k-max; the spectrum is extrapolated beyond it, see
 * #NcPowspecMLCBE.
 *
 */
void
nc_powspec_ml_cbe_set_intern_k_max (NcPowspecMLCBE *ps_cbe, const gdouble k_max)
{
  NcPowspecMLCBEPrivate * const self = nc_powspec_ml_cbe_get_instance_private (ps_cbe);

  self->intern_k_max = k_max;

  if (self->intern_k_min > 0.0)
    g_assert_cmpfloat (self->intern_k_min, <, self->intern_k_max);
}

/**
 * nc_powspec_ml_cbe_get_intern_k_min:
 * @ps_cbe: a #NcPowspecMLCBE
 *
 * Returns: the current value of the minimum mode $k$ computed by CLASS.
 */
gdouble
nc_powspec_ml_cbe_get_intern_k_min (NcPowspecMLCBE *ps_cbe)
{
  NcPowspecMLCBEPrivate * const self = nc_powspec_ml_cbe_get_instance_private (ps_cbe);

  return self->intern_k_min;
}

/**
 * nc_powspec_ml_cbe_get_intern_k_max :
 * @ps_cbe: a #NcPowspecMLCBE
 *
 * Returns: the current value of the maximum mode $k$ computed by CLASS.
 */
gdouble
nc_powspec_ml_cbe_get_intern_k_max (NcPowspecMLCBE *ps_cbe)
{
  NcPowspecMLCBEPrivate * const self = nc_powspec_ml_cbe_get_instance_private (ps_cbe);

  return self->intern_k_max;
}

