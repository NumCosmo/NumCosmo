/***************************************************************************
 *            nc_hicosmo_de.c
 *
 *  Tue Mar 18 15:33:13 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
 * SECTION:nc_hicosmo_de
 * @title: NcHICosmoDE
 * @short_description: Abstract class for implementing dark energy models.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_de.h"
#include "model/nc_hicosmo_de_reparam_ok.h"

G_DEFINE_ABSTRACT_TYPE (NcHICosmoDE, nc_hicosmo_de, NC_TYPE_HICOSMO);

#define VECTOR    (NCM_MODEL (cosmo)->params)
#define MACRO_H0  (ncm_vector_get (VECTOR, NC_HICOSMO_DE_H0))
#define OMEGA_C   (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_C))
#define OMEGA_X   (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_X))
#define T_GAMMA0  (ncm_vector_get (VECTOR, NC_HICOSMO_DE_T_GAMMA0))
#define HE_YP     (ncm_vector_get (VECTOR, NC_HICOSMO_DE_HE_YP))
#define ENNU      (ncm_vector_get (VECTOR, NC_HICOSMO_DE_ENNU))
#define OMEGA_R   nc_hicosmo_Omega_r0 (NC_HICOSMO (cosmo))
#define OMEGA_B   (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_B))
#define SPECINDEX (ncm_vector_get (VECTOR, NC_HICOSMO_DE_SPECINDEX))
#define SIGMA8    (ncm_vector_get (VECTOR, NC_HICOSMO_DE_SIGMA8))

#define OMEGA_M (OMEGA_B + OMEGA_C)
#define OMEGA_K (1.0 - (OMEGA_B + OMEGA_C + OMEGA_R + OMEGA_X))

/****************************************************************************
 * Normalized Hubble function
 ****************************************************************************/

static gdouble
_nc_hicosmo_de_E2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble omega_k = OMEGA_K;
  const gdouble x = 1.0 + z;
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;
  const gdouble x4 = x3 * x;
  const gdouble E2 = OMEGA_R * x4 + OMEGA_M * x3 + omega_k * x2 + nc_hicosmo_de_weff (NC_HICOSMO_DE (cosmo), z);
  return E2;
}

/****************************************************************************
 * dE2_dz
 ****************************************************************************/

static gdouble
_nc_hicosmo_de_dE2_dz (NcHICosmo *cosmo, gdouble z)
{
  const gdouble omega_k = OMEGA_K;
  const gdouble x = 1.0 + z;
  const gdouble x2 = x*x;
  const gdouble x3 = x2*x;

  return (4.0 * OMEGA_R * x3 + 3.0 * OMEGA_M * x2 + 2.0 * omega_k * x + nc_hicosmo_de_dweff_dz (NC_HICOSMO_DE (cosmo), z));
}

/****************************************************************************
 * d2E2_dz2
 ****************************************************************************/

static gdouble
_nc_hicosmo_de_d2E2_dz2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble omega_k = OMEGA_K;
  const gdouble x = 1.0 + z;
  const gdouble x2 = x*x;

  return (12.0 * OMEGA_R * x2 + 6.0 * OMEGA_M * x + 2.0 * omega_k + nc_hicosmo_de_d2weff_dz2 (NC_HICOSMO_DE (cosmo), z));
}

/****************************************************************************
 * Simple functions
 ****************************************************************************/
static gdouble _nc_hicosmo_de_H0 (NcHICosmo *cosmo) { return MACRO_H0; }
static gdouble _nc_hicosmo_de_Omega_t0 (NcHICosmo *cosmo) { return OMEGA_M + OMEGA_X + OMEGA_R; }
static gdouble _nc_hicosmo_de_Omega_c0 (NcHICosmo *cosmo) { return OMEGA_C; }
static gdouble _nc_hicosmo_de_T_gamma0 (NcHICosmo *cosmo) { return T_GAMMA0; }

static gdouble 
_nc_hicosmo_de_Yp_4He (NcHICosmo *cosmo) 
{
  NcmModel *model       = NCM_MODEL (cosmo);
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (cosmo);
  if (ncm_model_param_get_ftype (model, NC_HICOSMO_DE_DEFAULT_HE_YP) == NCM_PARAM_TYPE_FIXED)
  {
    const gdouble wb     = nc_hicosmo_Omega_b0h2 (cosmo);
    const gdouble DENNU  = ENNU - 3.046;

    if (FALSE)
    {
      const gdouble wb2    = wb * wb;
      const gdouble DENNU2 = DENNU * DENNU;

      const gdouble Yp = 0.2311 + 0.9502 * wb - 11.27 * wb2 + 
        DENNU * (0.01356 + 0.008581 * wb - 0.1810 * wb2) +
        DENNU2 * (-0.0009795 - 0.001370 * wb + 0.01746 * wb2);

      return Yp;
    }
    else
    {
      
      if (model->pkey != cosmo_de->HE4_Yp_key)
      {
        const gdouble Yp = ncm_spline2d_eval (NC_HICOSMO_DE (cosmo)->BBN_spline2d, wb, DENNU);
        ncm_vector_set (VECTOR, NC_HICOSMO_DE_HE_YP, Yp);
        cosmo_de->HE4_Yp_key = model->pkey;
        /*printf ("# omega_b % 20.15g DeltaNnu % 20.15g Yp % 20.15g\n",  wb, DENNU, Yp);*/
      }
    }
  }
  
  return HE_YP;
}

static gdouble _nc_hicosmo_de_Omega_g0 (NcHICosmo *cosmo)
{
  const gdouble h = MACRO_H0 / 100.0;
  const gdouble h2 = h * h;
  return ncm_c_radiation_temp_to_h2omega_r (T_GAMMA0) / h2;
}
static gdouble
_nc_hicosmo_de_Omega_nu0 (NcHICosmo *cosmo)
{
  const gdouble conv = 7.0 / 8.0 * pow (4.0 / 11.0, 4.0 / 3.0);
  return ENNU * conv * _nc_hicosmo_de_Omega_g0 (cosmo);
}
static gdouble
_nc_hicosmo_de_Omega_r0 (NcHICosmo *cosmo)
{
  const gdouble conv = 7.0 / 8.0 * pow (4.0 / 11.0, 4.0 / 3.0);
  return (1.0 + ENNU * conv) * _nc_hicosmo_de_Omega_g0 (cosmo);
}
static gdouble _nc_hicosmo_de_Omega_b0 (NcHICosmo *cosmo) { return OMEGA_B; }
static gdouble
_nc_hicosmo_de_bgp_cs2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble x        = 1.0 + z;
  const gdouble Omega_g0 = _nc_hicosmo_de_Omega_g0 (cosmo);
  const gdouble Omega_b0 = _nc_hicosmo_de_Omega_b0 (cosmo);
  const gdouble nine_4   = 9.0 / 4.0;
  
  return 1.0 / (3.0 + nine_4 * Omega_b0 / (Omega_g0 * x));
}

static gdouble _nc_hicosmo_de_sigma_8 (NcHICosmo *cosmo) { return SIGMA8; }
static gdouble _nc_hicosmo_de_powspec (NcHICosmo *cosmo, gdouble k) { return pow (k, SPECINDEX); }

void
nc_hicosmo_de_set_wmap5_params (NcHICosmoDE *cosmo_de)
{
  g_assert (NC_IS_HICOSMO_DE (cosmo_de));
  ncm_model_orig_param_set (NCM_MODEL (cosmo_de), NC_HICOSMO_DE_H0,       72.4000);
  ncm_model_orig_param_set (NCM_MODEL (cosmo_de), NC_HICOSMO_DE_OMEGA_C,   0.2060);
  ncm_model_orig_param_set (NCM_MODEL (cosmo_de), NC_HICOSMO_DE_OMEGA_X,   0.7510);
  ncm_model_orig_param_set (NCM_MODEL (cosmo_de), NC_HICOSMO_DE_OMEGA_B,   0.0432);
  ncm_model_orig_param_set (NCM_MODEL (cosmo_de), NC_HICOSMO_DE_T_GAMMA0,  2.7250);
  ncm_model_orig_param_set (NCM_MODEL (cosmo_de), NC_HICOSMO_DE_ENNU,      3.046);
  ncm_model_orig_param_set (NCM_MODEL (cosmo_de), NC_HICOSMO_DE_SPECINDEX, 0.9610);
  ncm_model_orig_param_set (NCM_MODEL (cosmo_de), NC_HICOSMO_DE_SIGMA8,    0.7870);
}

/**
 * nc_hicosmo_de_omega_x2omega_k:
 * @cosmo_de: FIXME
 *
 * FIXME
 *
 */
void
nc_hicosmo_de_omega_x2omega_k (NcHICosmoDE *cosmo_de)
{
  NcHICosmoDEReparamOk *de_reparam_ok = nc_hicosmo_de_reparam_ok_new (ncm_model_len (NCM_MODEL (cosmo_de)));
  ncm_model_set_reparam (NCM_MODEL (cosmo_de), NCM_REPARAM (de_reparam_ok));
  return;
}

static void
bbn_prior (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (ncm_mset_peek (mset, nc_hicosmo_id ()));
  gdouble z_bbn = 1.0e9;
  gdouble bbn, a;

  NCM_UNUSED (obj);
  NCM_UNUSED (x);

  a = nc_hicosmo_de_weff (cosmo_de, z_bbn) / nc_hicosmo_E2 (NC_HICOSMO (cosmo_de), z_bbn);
  bbn = 1.0 / sqrt(1.0 - a);
  f[0] = (bbn - 0.942) / 0.03;
}

/**
 * nc_hicosmo_de_new_add_bbn:
 * @lh: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_hicosmo_de_new_add_bbn (NcmLikelihood *lh)
{
  NcmMSetFunc *func = ncm_mset_func_new (bbn_prior, 0, 1, NULL, NULL);
  ncm_likelihood_priors_add (lh, func, FALSE);
  return TRUE;
}

enum {
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_de_init (NcHICosmoDE *cosmo_de)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *filename = ncm_cfg_get_data_filename ("BBN_spline2d.obj", TRUE);

  cosmo_de->BBN_spline2d = NCM_SPLINE2D (ncm_serialize_from_file (ser, filename));
  g_assert (NCM_IS_SPLINE2D (cosmo_de->BBN_spline2d));

  ncm_serialize_clear (&ser);

  cosmo_de->HE4_Yp_key = NCM_MODEL (cosmo_de)->pkey - 1;
}

static void
nc_hicosmo_de_dispose (GObject *object)
{
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (object);

  ncm_spline2d_clear (&cosmo_de->BBN_spline2d);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_de_parent_class)->finalize (object);
}

static void
nc_hicosmo_de_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_de_parent_class)->finalize (object);
}

static void
nc_hicosmo_de_class_init (NcHICosmoDEClass *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass* parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (klass);

  object_class->dispose      = &nc_hicosmo_de_dispose;
  object_class->finalize     = &nc_hicosmo_de_finalize;

  ncm_model_class_set_name_nick (model_class, "Darkenergy models abstract class", "NcHICosmoDE");
  ncm_model_class_add_params (model_class, NC_HICOSMO_DE_SPARAM_LEN, 0, PROP_SIZE);
  /* Set H_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_H0, "H_0", "H0",
                              40.0, 120.0, 1.0,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_H0,
                              NCM_PARAM_TYPE_FIXED);
  /* Set Omega_c0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_C, "\\Omega_{c0}", "Omegac",
                              1e-8,  1.2, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_OMEGA_C,
                              NCM_PARAM_TYPE_FREE);
  /* Set Omega_x0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_X, "\\Omega_{x0}", "Omegax",
                              1e-8,  2.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_OMEGA_X,
                              NCM_PARAM_TYPE_FREE);
  /* Set T_gamma0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_T_GAMMA0, "T_{\\gamma0}", "Tgamma0",
                              2.0,  3.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_T_GAMMA0,
                              NCM_PARAM_TYPE_FIXED);
  /* Set He Yp param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_HE_YP, "Y_p", "Yp",
                              0.0,  1.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_HE_YP,
                              NCM_PARAM_TYPE_FIXED);
  /* Set ENnu param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_ENNU, "N_\\nu", "ENnu",
                              0.0,  4.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_ENNU,
                              NCM_PARAM_TYPE_FIXED);
  /* Set Omega_b0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_B, "\\Omega_{b0}", "Omegab",
                              0.03,  0.05, 5.0e-4,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_OMEGA_B,
                              NCM_PARAM_TYPE_FIXED);
  /* Set n_s param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_SPECINDEX, "n_s", "ns",
                              0.5,   1.5, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_SPECINDEX,
                              NCM_PARAM_TYPE_FIXED);
  /* Set sigma_8 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_SIGMA8, "\\sigma_8", "sigma8",
                              0.2,   1.8, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_SIGMA8,
                              NCM_PARAM_TYPE_FIXED);
  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hicosmo_set_H0_impl        (parent_class, &_nc_hicosmo_de_H0);
  nc_hicosmo_set_E2_impl        (parent_class, &_nc_hicosmo_de_E2);
  nc_hicosmo_set_Omega_c0_impl  (parent_class, &_nc_hicosmo_de_Omega_c0);
  nc_hicosmo_set_Omega_r0_impl  (parent_class, &_nc_hicosmo_de_Omega_r0);
  nc_hicosmo_set_Omega_b0_impl  (parent_class, &_nc_hicosmo_de_Omega_b0);
  nc_hicosmo_set_Omega_g0_impl  (parent_class, &_nc_hicosmo_de_Omega_g0);
  nc_hicosmo_set_Omega_nu0_impl (parent_class, &_nc_hicosmo_de_Omega_nu0);
  nc_hicosmo_set_Omega_t0_impl  (parent_class, &_nc_hicosmo_de_Omega_t0);
  nc_hicosmo_set_sigma_8_impl   (parent_class, &_nc_hicosmo_de_sigma_8);
  nc_hicosmo_set_T_gamma0_impl  (parent_class, &_nc_hicosmo_de_T_gamma0);
  nc_hicosmo_set_Yp_4He_impl    (parent_class, &_nc_hicosmo_de_Yp_4He);
  nc_hicosmo_set_dE2_dz_impl    (parent_class, &_nc_hicosmo_de_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl  (parent_class, &_nc_hicosmo_de_d2E2_dz2);
  nc_hicosmo_set_bgp_cs2_impl   (parent_class, &_nc_hicosmo_de_bgp_cs2);
  nc_hicosmo_set_powspec_impl   (parent_class, &_nc_hicosmo_de_powspec);
}

#define NC_HICOSMO_DE_SET_IMPL_FUNC(name) \
void \
nc_hicosmo_de_set_##name##_impl (NcHICosmoDEClass *cosmo_de_class, NcmFuncF f, NcmFuncPF pf, NcmFuncDF df) \
{ \
NCM_MODEL_CLASS (cosmo_de_class)->impl |= NC_HICOSMO_DE_IMPL_##name; \
g_assert (f != NULL); \
cosmo_de_class->name = *ncm_func_stub; \
cosmo_de_class->name.f = f; \
cosmo_de_class->name.pf = pf; \
cosmo_de_class->name.df = df; \
}

/**
 * nc_hicosmo_de_set_weff_impl: (skip)
 * @cosmo_de_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO_DE,NcHICosmoDE,nc_hicosmo_de,NcHICosmoDEFunc1,weff)
/**
 * nc_hicosmo_de_set_dweff_dz_impl: (skip)
 * @cosmo_de_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO_DE,NcHICosmoDE,nc_hicosmo_de,NcHICosmoDEFunc1,dweff_dz)
/**
 * nc_hicosmo_de_set_d2weff_dz2_impl: (skip)
 * @cosmo_de_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO_DE,NcHICosmoDE,nc_hicosmo_de,NcHICosmoDEFunc1,d2weff_dz2)
/**
 * nc_hicosmo_weff:
 * @cosmo_de: a #NcHICosmoDE
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_dweff_dz:
 * @cosmo_de: a #NcHICosmoDE
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_d2weff_dz2:
 * @cosmo_de: a #NcHICosmoDE
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
