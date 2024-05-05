/***************************************************************************
 *            nc_hiprim_two_fluids.c
 *
 *  Tue April 27 11:14:53 2024
 *  Copyright  2024  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hiprim_two_fluids.c
 * Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:nc_hiprim_two_fluids
 * @title: NcHIPrimTwoFluids
 * @short_description: Two Fluids implementation for primordial spectra.
 *
 * Primordial adiabatic scalar power spectrum obtained from a two fluids model.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hiprim_two_fluids.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"

typedef struct _NcHIPrimTwoFluidsPrivate
{
  NcmSpline2d *lnSA_powspec_lnk_lnw;
  gboolean use_default_calib;
} NcHIPrimTwoFluidsPrivate;


G_DEFINE_TYPE_WITH_PRIVATE (NcHIPrimTwoFluids, nc_hiprim_two_fluids, NC_TYPE_HIPRIM)

enum
{
  PROP_0,
  PROP_LNK_LNW_SPLINE,
  PROP_USE_DEFAULT_CALIB,
  PROP_SIZE,
};

static void
nc_hiprim_two_fluids_init (NcHIPrimTwoFluids *two_fluids)
{
  NcHIPrimTwoFluidsPrivate * const self = nc_hiprim_two_fluids_get_instance_private (two_fluids);

  self->lnSA_powspec_lnk_lnw = NULL;
  self->use_default_calib    = FALSE;
}

static void
_nc_hiprim_two_fluids_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPrimTwoFluids *two_fluids         = NC_HIPRIM_TWO_FLUIDS (object);
  NcHIPrimTwoFluidsPrivate * const self = nc_hiprim_two_fluids_get_instance_private (two_fluids);

  switch (prop_id)
  {
    case PROP_LNK_LNW_SPLINE:
      nc_hiprim_two_fluids_set_lnk_lnw_spline (two_fluids, g_value_get_object (value));
      break;
    case PROP_USE_DEFAULT_CALIB:
      nc_hiprim_two_fluids_set_use_default_calib (two_fluids, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hiprim_two_fluids_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIPrimTwoFluids *two_fluids         = NC_HIPRIM_TWO_FLUIDS (object);
  NcHIPrimTwoFluidsPrivate * const self = nc_hiprim_two_fluids_get_instance_private (two_fluids);

  switch (prop_id)
  {
    case PROP_LNK_LNW_SPLINE:
      g_value_set_object (value, nc_hiprim_two_fluids_peek_lnk_lnw_spline (two_fluids));
      break;
    case PROP_USE_DEFAULT_CALIB:
      g_value_set_boolean (value, nc_hiprim_two_fluids_get_use_default_calib (two_fluids));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_hiprim_two_fluids_dispose (GObject *object)
{
  NcHIPrimTwoFluids *two_fluids         = NC_HIPRIM_TWO_FLUIDS (object);
  NcHIPrimTwoFluidsPrivate * const self = nc_hiprim_two_fluids_get_instance_private (two_fluids);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hiprim_two_fluids_parent_class)->dispose (object);
}

static void
_nc_hiprim_two_fluids_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hiprim_two_fluids_parent_class)->finalize (object);
}

static gdouble _nc_hiprim_two_fluids_lnSA_powespec_lnk (NcHIPrim *prim, const gdouble lnk);
static gdouble _nc_hiprim_two_fluids_lnT_powespec_lnk (NcHIPrim *prim, const gdouble lnk);

static void
nc_hiprim_two_fluids_class_init (NcHIPrimTwoFluidsClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcHIPrimClass *prim_class  = NC_HIPRIM_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize    = &_nc_hiprim_two_fluids_finalize;
  object_class->dispose     = &_nc_hiprim_two_fluids_dispose;
  model_class->set_property = &_nc_hiprim_two_fluids_set_property;
  model_class->get_property = &_nc_hiprim_two_fluids_get_property;

  ncm_model_class_set_name_nick (model_class, "Power Law model for primordial spectra", "TwoFluids");
  ncm_model_class_add_params (model_class, NC_HIPRIM_TWO_FLUIDS_SPARAM_LEN, 0, PROP_SIZE);

  /* Set ln10e10ASA param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_TWO_FLUIDS_LN10E10ASA, "\\log(10^{10}A_{\\mathrm{SA}})", "ln10e10ASA",
                              2.0, 5.0, 1.0e0,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_TWO_FLUIDS_DEFAULT_LN10E10ASA,
                              NCM_PARAM_TYPE_FIXED);

  /* Set T_SA_ratio param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_TWO_FLUIDS_T_SA_RATIO, "A_T/A_{\\mathrm{SA}}", "T_SA_ratio",
                              0.0, 10.0, 1.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_TWO_FLUIDS_DEFAULT_T_SA_RATIO,
                              NCM_PARAM_TYPE_FIXED);

  /* Set ln(k_0) param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_TWO_FLUIDS_LNK0, "\\ln(k_0)", "lnk0",
                              -5.0 * M_LN10, -2.0 * M_LN10, 1.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_TWO_FLUIDS_DEFAULT_LNK0,
                              NCM_PARAM_TYPE_FIXED);

  /* Set ln(w) param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_TWO_FLUIDS_LNW, "\\ln(w)", "lnw",
                              -7.0 * M_LN10, -1.0 * M_LN10, 1.0e-1,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_TWO_FLUIDS_DEFAULT_LNW,
                              NCM_PARAM_TYPE_FIXED);

  /* Set N_T param info */
  ncm_model_class_set_sparam (model_class, NC_HIPRIM_TWO_FLUIDS_N_T, "n_{\\mathrm{T}}", "n_T",
                              -0.5, 0.5, 1.0e-2,
                              NC_HIPRIM_DEFAULT_PARAMS_ABSTOL, NC_HIPRIM_TWO_FLUIDS_DEFAULT_N_T,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  g_object_class_install_property (object_class,
                                   PROP_LNK_LNW_SPLINE,
                                   g_param_spec_object ("lnk-lnw-spline",
                                                        NULL,
                                                        "Spline for the primordial adiabatic scalar power spectrum as a function of ln(k) and ln(w)",
                                                        NCM_TYPE_SPLINE2D,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_USE_DEFAULT_CALIB,
                                   g_param_spec_boolean ("use-default-calib",
                                                         NULL,
                                                         "Use default calibration",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  nc_hiprim_set_lnSA_powspec_lnk_impl (prim_class, &_nc_hiprim_two_fluids_lnSA_powespec_lnk);
  nc_hiprim_set_lnT_powspec_lnk_impl  (prim_class, &_nc_hiprim_two_fluids_lnT_powespec_lnk);
}

#define VECTOR     (NCM_MODEL (prim))
#define LN10E10ASA (ncm_model_orig_param_get (VECTOR, NC_HIPRIM_TWO_FLUIDS_LN10E10ASA))
#define T_SA_RATIO (ncm_model_orig_param_get (VECTOR, NC_HIPRIM_TWO_FLUIDS_T_SA_RATIO))
#define LNK0       (ncm_model_orig_param_get (VECTOR, NC_HIPRIM_TWO_FLUIDS_LNK0))
#define LNW        (ncm_model_orig_param_get (VECTOR, NC_HIPRIM_TWO_FLUIDS_LNW))
#define N_T        (ncm_model_orig_param_get (VECTOR, NC_HIPRIM_TWO_FLUIDS_N_T))

/****************************************************************************
 * Power spectra
 ****************************************************************************/

static gdouble
_nc_hiprim_two_fluids_lnSA_powespec_lnk (NcHIPrim *prim, const gdouble lnk)
{
  NcHIPrimTwoFluids *two_fluids         = NC_HIPRIM_TWO_FLUIDS (prim);
  NcHIPrimTwoFluidsPrivate * const self = nc_hiprim_two_fluids_get_instance_private (two_fluids);
  const gdouble ln_ka                   = lnk - LNK0;
  const gdouble lnw                     = LNW;

  return LN10E10ASA - 10.0 * M_LN10 + ncm_spline2d_eval (self->lnSA_powspec_lnk_lnw, ln_ka, lnw);
}

static gdouble
_nc_hiprim_two_fluids_lnT_powespec_lnk (NcHIPrim *prim, const gdouble lnk)
{
  const gdouble ln_ka = lnk - prim->lnk_pivot;

  return N_T * ln_ka + LN10E10ASA - 10.0 * M_LN10 + log (T_SA_RATIO);
}

/**
 * nc_hiprim_two_fluids_new: (constructor)
 *
 * This function instantiates a new object of type #NcHIPrimTwoFluids.
 *
 * Returns: (transfer full): A new #NcHIPrimTwoFluids
 */
NcHIPrimTwoFluids *
nc_hiprim_two_fluids_new (void)
{
  NcHIPrimTwoFluids *prim_pl = g_object_new (NC_TYPE_HIPRIM_TWO_FLUIDS,
                                             NULL);

  return prim_pl;
}

/**
 * nc_hiprim_two_fluids_set_use_default_calib:
 * @two_fluids: a #NcHIPrimTwoFluids
 * @use_default_calib: a #gboolean
 *
 * Set the use_default_calib flag.
 *
 */
void
nc_hiprim_two_fluids_set_use_default_calib (NcHIPrimTwoFluids *two_fluids, gboolean use_default_calib)
{
  NcHIPrimTwoFluidsPrivate * const self = nc_hiprim_two_fluids_get_instance_private (two_fluids);

  if ((use_default_calib && self->use_default_calib) || (!use_default_calib && !self->use_default_calib))
    return;

  if (use_default_calib)
  {
    NcmSerialize *ser         = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    gchar *default_calib_file = ncm_cfg_get_data_filename ("hiprim_2f_spline.bin", TRUE);

    ncm_spline2d_clear (&self->lnSA_powspec_lnk_lnw);

    self->lnSA_powspec_lnk_lnw = NCM_SPLINE2D (ncm_serialize_from_binfile (ser, default_calib_file));
    g_assert_nonnull (self->lnSA_powspec_lnk_lnw);
    g_assert (NCM_IS_SPLINE2D (self->lnSA_powspec_lnk_lnw));
    g_assert (ncm_spline2d_is_init (self->lnSA_powspec_lnk_lnw));

    g_free (default_calib_file);
    ncm_serialize_free (ser);
  }

  self->use_default_calib = use_default_calib;
}

/**
 * nc_hiprim_two_fluids_get_use_default_calib:
 * @two_fluids: a #NcHIPrimTwoFluids
 *
 * Get the use_default_calib flag.
 *
 * Returns: a #gboolean
 */
gboolean
nc_hiprim_two_fluids_get_use_default_calib (NcHIPrimTwoFluids *two_fluids)
{
  NcHIPrimTwoFluidsPrivate * const self = nc_hiprim_two_fluids_get_instance_private (two_fluids);

  return self->use_default_calib;
}

/**
 * nc_hiprim_two_fluids_set_lnk_lnw_spline:
 * @two_fluids: a #NcHIPrimTwoFluids
 * @lnSA_powspec_lnk_lnw: a #NcmSpline2d
 *
 * Set the spline for the primordial adiabatic scalar power spectrum as a function of ln(k) and ln(w).
 *
 */
void
nc_hiprim_two_fluids_set_lnk_lnw_spline (NcHIPrimTwoFluids *two_fluids, NcmSpline2d *lnSA_powspec_lnk_lnw)
{
  NcHIPrimTwoFluidsPrivate * const self = nc_hiprim_two_fluids_get_instance_private (two_fluids);

  if (self->use_default_calib)
  {
    g_error ("nc_hiprim_two_fluids_set_lnk_lnw_spline: Cannot set spline when use_default_calib is TRUE.");

    return;
  }

  ncm_spline2d_clear (&self->lnSA_powspec_lnk_lnw);

  self->lnSA_powspec_lnk_lnw = lnSA_powspec_lnk_lnw;
}

/**
 * nc_hiprim_two_fluids_peek_lnk_lnw_spline:
 * @two_fluids: a #NcHIPrimTwoFluids
 *
 * Get the spline for the primordial adiabatic scalar power spectrum as a function of
 * $\ln(k)$ and $\ln(w)$.
 *
 * Returns: (transfer none): the spline for the primordial power spectrum.
 */
NcmSpline2d *
nc_hiprim_two_fluids_peek_lnk_lnw_spline (NcHIPrimTwoFluids *two_fluids)
{
  NcHIPrimTwoFluidsPrivate * const self = nc_hiprim_two_fluids_get_instance_private (two_fluids);

  if (self->use_default_calib)
    return NULL;

  return self->lnSA_powspec_lnk_lnw;
}

