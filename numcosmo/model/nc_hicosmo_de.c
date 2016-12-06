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
#include "model/nc_hicosmo_de_reparam_cmb.h"
#include "math/ncm_mset_func_list.h"
#include "math/ncm_prior_gauss_func.h"

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

#define OMEGA_M (OMEGA_B + OMEGA_C)
#define OMEGA_K (1.0 - (OMEGA_B + OMEGA_C + OMEGA_R + OMEGA_X))

/****************************************************************************
 * Normalized Hubble function
 ****************************************************************************/

static gdouble
_nc_hicosmo_de_E2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble Omega_k = OMEGA_K;
  const gdouble x = 1.0 + z;
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;
  const gdouble x4 = x3 * x;
  const gdouble E2 = OMEGA_R * x4 + OMEGA_M * x3 + Omega_k * x2 + nc_hicosmo_de_E2Omega_de (NC_HICOSMO_DE (cosmo), z);
  return E2;
}

/****************************************************************************
 * dE2_dz
 ****************************************************************************/

static gdouble
_nc_hicosmo_de_dE2_dz (NcHICosmo *cosmo, gdouble z)
{
  const gdouble Omega_k = OMEGA_K;
  const gdouble x = 1.0 + z;
  const gdouble x2 = x*x;
  const gdouble x3 = x2*x;

  return (4.0 * OMEGA_R * x3 + 3.0 * OMEGA_M * x2 + 2.0 * Omega_k * x + nc_hicosmo_de_dE2Omega_de_dz (NC_HICOSMO_DE (cosmo), z));
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

  return (12.0 * OMEGA_R * x2 + 6.0 * OMEGA_M * x + 2.0 * omega_k + nc_hicosmo_de_d2E2Omega_de_dz2 (NC_HICOSMO_DE (cosmo), z));
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
        ncm_model_orig_param_set (model, NC_HICOSMO_DE_HE_YP, Yp);
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

static guint
_nc_hicosmo_de_NMassNu (NcHICosmo *cosmo)
{
  return ncm_model_vparam_len (NCM_MODEL (cosmo), NC_HICOSMO_DE_MASSNU_M);
}

static void
_nc_hicosmo_de_MassNuInfo (NcHICosmo *cosmo, guint nu_i, gdouble *mass_eV, gdouble *T_0, gdouble *xi, gdouble *g)
{
  NcmModel *model = NCM_MODEL (cosmo);

  g_assert_cmpuint (nu_i, <, ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_T));

  mass_eV[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_M, nu_i);
  
  if (ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_T) == 1)
  {
    T_0[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_T, 0);
  }
  else
  {
    T_0[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_T, nu_i);
  }

  if (ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_XI) == 1)
  {
    xi[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_XI, 0);
  }
  else
  {
    xi[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_XI, nu_i);
  }
  
  if (ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_G) == 1)
  {
    g[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_G, 0);
  }
  else
  {
    g[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_G, nu_i);
  }  
}

static gdouble 
_nc_hicosmo_Omega_mnu0 (NcHICosmo *cosmo, const guint nu_i, const gdouble z) 
{ 
  
  return 0.0; 
}

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
}

/**
 * nc_hicosmo_de_omega_x2omega_k:
 * @cosmo_de: a #NcHICosmoDE
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

/**
 * nc_hicosmo_de_cmb_params:
 * @cosmo_de: a #NcHICosmoDE
 *
 * FIXME
 *
 */
void
nc_hicosmo_de_cmb_params (NcHICosmoDE *cosmo_de)
{
  NcHICosmoDEReparamCMB *de_reparam_cmb = nc_hicosmo_de_reparam_cmb_new (ncm_model_len (NCM_MODEL (cosmo_de)));
  ncm_model_set_reparam (NCM_MODEL (cosmo_de), NCM_REPARAM (de_reparam_cmb));
  return;
}

/**
 * nc_hicosmo_de_new_add_bbn:
 * @lh: a #NcmLikelihood
 *
 * FIXME
 *
 */
void
nc_hicosmo_de_new_add_bbn (NcmLikelihood *lh)
{
  NcmMSetFunc *bbn       = NCM_MSET_FUNC (ncm_mset_func_list_new ("NcHICosmoDE:BBN", NULL));
  NcmPriorGaussFunc *pgf = ncm_prior_gauss_func_new (bbn, 0.942, 0.03, 0.0);
  
  ncm_likelihood_priors_add (lh, NCM_PRIOR (pgf));

  ncm_mset_func_clear (&bbn);
  ncm_prior_gauss_func_clear (&pgf);
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

  g_free (filename);
}

static void
_nc_hicosmo_de_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hicosmo_de_parent_class)->constructed (object);
  {
    NcmModel *model = NCM_MODEL (object);
    const guint m_len  = ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_M);
    const guint T_len  = ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_T);
    const guint xi_len = ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_XI);
    const guint g_len  = ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_G);
    
    if (!((m_len == T_len) || (m_len > 0 && T_len == 1)))
    {
      g_error ("NcHICosmoDE: number of neutrinos masses must match the number of massive neutrino temperatures,\n"
               " or the neutrino temperature vector must be of size one to use the same value for all massive neutrinos.");
    }
    if (!((m_len == xi_len) || (m_len > 0 && xi_len == 1)))
    {
      g_error ("NcHICosmoDE: number of neutrinos masses must match the number of massive neutrino relative chemical potential,\n"
               " or the neutrino relative chemical potential vector must be of size one to use the same value for all massive neutrinos.");
    }
    if (!((m_len == g_len) || (m_len > 0 && g_len == 1)))
    {
      g_error ("NcHICosmoDE: number of neutrinos masses must match the number of massive neutrino degeneracy,\n"
               " or the neutrino degeneracy vector must be of size one to use the same value for all massive neutrinos.");
    }
  }
}

static void
_nc_hicosmo_de_dispose (GObject *object)
{
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (object);

  ncm_spline2d_clear (&cosmo_de->BBN_spline2d);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_de_parent_class)->dispose (object);
}

static void
_nc_hicosmo_de_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_de_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_de_E2Omega_de (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_dE2Omega_de_dz (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_d2E2Omega_de_dz2 (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_w_de (NcHICosmoDE *cosmo_de, gdouble z);

static void
nc_hicosmo_de_class_init (NcHICosmoDEClass *klass)
{
  GObjectClass* object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass* parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (klass);

  object_class->constructed = &_nc_hicosmo_de_constructed;
  object_class->dispose     = &_nc_hicosmo_de_dispose;
  object_class->finalize    = &_nc_hicosmo_de_finalize;

  ncm_model_class_set_name_nick (model_class, "Darkenergy models abstract class", "NcHICosmoDE");
  ncm_model_class_add_params (model_class, NC_HICOSMO_DE_SPARAM_LEN, NC_HICOSMO_DE_VPARAM_LEN, PROP_SIZE);

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

  /* Set massive neutrinos mass vector param */
  ncm_model_class_set_vparam (model_class, NC_HICOSMO_DE_MASSNU_M, 0, "m_\\nu", "massnu", 
                              0.0, 10.0, 0.01, 
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_NU_MASS, 
                              NCM_PARAM_TYPE_FIXED);

  /* Set massive neutrinos temperature vector param */
  ncm_model_class_set_vparam (model_class, NC_HICOSMO_DE_MASSNU_T, 0, "T_{\\nu0}", "Tnu", 
                              0.0, 10.0, 0.01, 
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_NU_T, 
                              NCM_PARAM_TYPE_FIXED);

  /* Set massive neutrinos relative chemical potential vector param */
  ncm_model_class_set_vparam (model_class, NC_HICOSMO_DE_MASSNU_XI, 0, "\\xi_{\\nu}", "xinu", 
                              -10.0, 10.0, 0.01, 
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_NU_XI, 
                              NCM_PARAM_TYPE_FIXED);

  /* Set massive neutrinos degeneracy vector param */
  ncm_model_class_set_vparam (model_class, NC_HICOSMO_DE_MASSNU_G, 0, "g_{\\nu}", "gnu", 
                              0.0, 10.0, 0.01, 
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_NU_G, 
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hicosmo_set_H0_impl         (parent_class, &_nc_hicosmo_de_H0);
  nc_hicosmo_set_E2_impl         (parent_class, &_nc_hicosmo_de_E2);
  nc_hicosmo_set_Omega_c0_impl   (parent_class, &_nc_hicosmo_de_Omega_c0);
  nc_hicosmo_set_Omega_r0_impl   (parent_class, &_nc_hicosmo_de_Omega_r0);
  nc_hicosmo_set_Omega_b0_impl   (parent_class, &_nc_hicosmo_de_Omega_b0);
  nc_hicosmo_set_Omega_g0_impl   (parent_class, &_nc_hicosmo_de_Omega_g0);
  nc_hicosmo_set_Omega_nu0_impl  (parent_class, &_nc_hicosmo_de_Omega_nu0);
  nc_hicosmo_set_Omega_t0_impl   (parent_class, &_nc_hicosmo_de_Omega_t0);
  nc_hicosmo_set_T_gamma0_impl   (parent_class, &_nc_hicosmo_de_T_gamma0);
  nc_hicosmo_set_Yp_4He_impl     (parent_class, &_nc_hicosmo_de_Yp_4He);
  nc_hicosmo_set_dE2_dz_impl     (parent_class, &_nc_hicosmo_de_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl   (parent_class, &_nc_hicosmo_de_d2E2_dz2);
  nc_hicosmo_set_bgp_cs2_impl    (parent_class, &_nc_hicosmo_de_bgp_cs2);
  nc_hicosmo_set_NMassNu_impl    (parent_class, &_nc_hicosmo_de_NMassNu);
  nc_hicosmo_set_MassNuInfo_impl (parent_class, &_nc_hicosmo_de_MassNuInfo);
  nc_hicosmo_set_Omega_mnu0_impl (parent_class, &_nc_hicosmo_Omega_mnu0);

  klass->E2Omega_de       = &_nc_hicosmo_de_E2Omega_de;
  klass->dE2Omega_de_dz   = &_nc_hicosmo_de_dE2Omega_de_dz;
  klass->d2E2Omega_de_dz2 = &_nc_hicosmo_de_d2E2Omega_de_dz2;
  klass->w_de             = &_nc_hicosmo_de_w_de;
}

static gdouble _nc_hicosmo_de_E2Omega_de (NcHICosmoDE *cosmo_de, gdouble z)       { g_error ("nc_hicosmo_de_E2Omega_de: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo_de)); return 0.0;  }
static gdouble _nc_hicosmo_de_dE2Omega_de_dz (NcHICosmoDE *cosmo_de, gdouble z)   { g_error ("nc_hicosmo_de_dE2Omega_de_dz: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo_de)); return 0.0;  }
static gdouble _nc_hicosmo_de_d2E2Omega_de_dz2 (NcHICosmoDE *cosmo_de, gdouble z) { g_error ("nc_hicosmo_de_d2E2Omega_de_dz2: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo_de)); return 0.0;  }
static gdouble _nc_hicosmo_de_w_de (NcHICosmoDE *cosmo_de, gdouble z)             { g_error ("nc_hicosmo_de_w_de: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo_de)); return 0.0;  }

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
 * nc_hicosmo_de_set_E2Omega_de_impl: (skip)
 * @cosmo_de_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO_DE,NcHICosmoDE,nc_hicosmo_de,NcHICosmoDEFunc1,E2Omega_de)
/**
 * nc_hicosmo_de_set_dE2Omega_de_dz_impl: (skip)
 * @cosmo_de_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO_DE,NcHICosmoDE,nc_hicosmo_de,NcHICosmoDEFunc1,dE2Omega_de_dz)
/**
 * nc_hicosmo_de_set_d2E2Omega_de_dz2_impl: (skip)
 * @cosmo_de_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO_DE,NcHICosmoDE,nc_hicosmo_de,NcHICosmoDEFunc1,d2E2Omega_de_dz2)
/**
 * nc_hicosmo_de_set_w_de_impl: (skip)
 * @cosmo_de_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO_DE,NcHICosmoDE,nc_hicosmo_de,NcHICosmoDEFunc1,w_de)

/**
 * nc_hicosmo_E2Omega_de:
 * @cosmo_de: a #NcHICosmoDE
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_dE2Omega_de_dz:
 * @cosmo_de: a #NcHICosmoDE
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_d2E2Omega_de_dz2:
 * @cosmo_de: a #NcHICosmoDE
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_w_de:
 * @cosmo_de: a #NcHICosmoDE
 * @z: redshift $z$
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_de_E2Omega_de_onepw:
 * @cosmo_de: a #NcHICosmoDE
 * @z: redshift $z$
 *
 * $E^2\Omega_\mathrm{de}(1+w)$.
 *
 * Returns: FIXME
 */

static void
_nc_hicosmo_de_flist_w0 (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *f)
{
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (ncm_mset_peek (mset, nc_hicosmo_id ()));

  g_assert (NC_IS_HICOSMO_DE (cosmo_de));
  
  f[0] = nc_hicosmo_de_w_de (cosmo_de, 0.0);
}

static void
_nc_hicosmo_de_flist_BBN (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *f)
{
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (ncm_mset_peek (mset, nc_hicosmo_id ()));
  gdouble z_bbn = 1.0e9;
  gdouble Omega_de;

  g_assert (NC_IS_HICOSMO_DE (cosmo_de));

  NCM_UNUSED (x);

  Omega_de = nc_hicosmo_de_E2Omega_de (cosmo_de, z_bbn) / nc_hicosmo_E2 (NC_HICOSMO (cosmo_de), z_bbn);

  f[0] = 1.0 / sqrt (1.0 - Omega_de);
}

void
_nc_hicosmo_de_register_functions (void)
{
  ncm_mset_func_list_register ("wDE", "\\omega_\\mathrm{de}", "NcHICosmoDE", "Darkenergy equation of state today", G_TYPE_NONE, _nc_hicosmo_de_flist_w0,  0, 1);
  ncm_mset_func_list_register ("BBN", "BBN",                  "NcHICosmoDE", "BBN",                                G_TYPE_NONE, _nc_hicosmo_de_flist_BBN, 0, 1);
}
