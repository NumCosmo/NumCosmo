/***************************************************************************
 *            nc_hicosmo_gcg.c
 *
 *  Wed March 08 13:48:46 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2017 <sandro@isoftware.com.br>
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
 * SECTION:nc_hicosmo_gcg
 * @title: NcHICosmoGCG
 * @short_description: Generalized Chaplygin Gas
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_mset_func_list.h"
#include "math/ncm_prior_gauss_func.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "model/nc_hicosmo_gcg.h"

struct _NcHICosmoGCGPrivate
{
  NcmSpline2d *BBN_spline2d;
  guint64 HE4_Yp_key;
  NcmIntegral1dPtr *nu_rho;
  NcmIntegral1dPtr *nu_p;
  NcmSpline *nu_rho_s[10];
  NcmSpline *nu_p_s[10];
  gdouble zmax;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcHICosmoGCG, nc_hicosmo_gcg, NC_TYPE_HICOSMO);

enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_gcg_init (NcHICosmoGCG *cosmo_de)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *filename   = ncm_cfg_get_data_filename ("BBN_spline2d.obj", TRUE);
  gint i;

  cosmo_de->priv               = nc_hicosmo_gcg_get_instance_private (cosmo_de);
  cosmo_de->priv->BBN_spline2d = NCM_SPLINE2D (ncm_serialize_from_file (ser, filename));

  g_assert (NCM_IS_SPLINE2D (cosmo_de->priv->BBN_spline2d));
  ncm_serialize_clear (&ser);
  cosmo_de->priv->HE4_Yp_key = NCM_MODEL (cosmo_de)->pkey - 1;

  cosmo_de->priv->nu_rho   = NULL;
  cosmo_de->priv->nu_p     = NULL;

  for (i = 0; i < 10; i++)
  {
    cosmo_de->priv->nu_rho_s[i] = NULL;
    cosmo_de->priv->nu_p_s[i]   = NULL;
  }

  cosmo_de->priv->zmax = 1.0e12;

  g_free (filename);
}

static gdouble _nc_hicosmo_gcg_neutrino_rho_integrand (gpointer userdata, const gdouble v, const gdouble w);
static gdouble _nc_hicosmo_gcg_neutrino_p_integrand (gpointer userdata, const gdouble v, const gdouble w);

#define _NC_HICOSMO_GCG_MNU_PREC (1.0e-6)

static void
_nc_hicosmo_gcg_constructed (GObject *object)
{
  /* Before model's construct!!!! */
  /* 
   * If other massive neutrinos parameters are not set, 
   * set it here before constructing the model. 
   */
  {
    NcmModel *model = NCM_MODEL (object);
    const guint nmassnu = ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_M);
  
    if (nmassnu > 0)
    {
      if (ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_T) == 0)
        g_array_index (model->vparam_len, guint, NC_HICOSMO_GCG_MASSNU_T) = nmassnu;

      if (ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_MU) == 0)
        g_array_index (model->vparam_len, guint, NC_HICOSMO_GCG_MASSNU_MU) = nmassnu;
      
      if (ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_G) == 0)
        g_array_index (model->vparam_len, guint, NC_HICOSMO_GCG_MASSNU_G) = nmassnu;      
    }
  }
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hicosmo_gcg_parent_class)->constructed (object);
  {
    NcmModel *model    = NCM_MODEL (object);
    const guint m_len  = ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_M);
    const guint T_len  = ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_T);
    const guint mu_len = ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_MU);
    const guint g_len  = ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_G);

    if (! ( (m_len == T_len) || (m_len > 0 && T_len == 1)))
    {
      g_error ("NcHICosmoGCG: number of neutrinos masses must match the number of massive neutrino temperatures,\n"
               " or the neutrino temperature vector must be of size one to use the same value for all massive neutrinos.");
    }

    if (! ( (m_len == mu_len) || (m_len > 0 && mu_len == 1)))
    {
      g_error ("NcHICosmoGCG: number of neutrinos masses must match the number of massive neutrino relative chemical potential,\n"
               " or the neutrino relative chemical potential vector must be of size one to use the same value for all massive neutrinos.");
    }

    if (! ( (m_len == g_len) || (m_len > 0 && g_len == 1)))
    {
      g_error ("NcHICosmoGCG: number of neutrinos masses must match the number of massive neutrino degeneracy,\n"
               " or the neutrino degeneracy vector must be of size one to use the same value for all massive neutrinos.");
    }

    if (m_len != 0)
    {
      NcHICosmoGCG *cosmo_de  = NC_HICOSMO_GCG (model);
      gint i;
      
      cosmo_de->priv->nu_rho = ncm_integral1d_ptr_new (&_nc_hicosmo_gcg_neutrino_rho_integrand, NULL);
      cosmo_de->priv->nu_p   = ncm_integral1d_ptr_new (&_nc_hicosmo_gcg_neutrino_p_integrand, NULL);

      ncm_integral1d_set_reltol (NCM_INTEGRAL1D (cosmo_de->priv->nu_rho), _NC_HICOSMO_GCG_MNU_PREC);
      ncm_integral1d_set_reltol (NCM_INTEGRAL1D (cosmo_de->priv->nu_p),   _NC_HICOSMO_GCG_MNU_PREC);

      ncm_integral1d_set_rule (NCM_INTEGRAL1D (cosmo_de->priv->nu_rho), 1);
      ncm_integral1d_set_rule (NCM_INTEGRAL1D (cosmo_de->priv->nu_p), 1);

      for (i = 0; i < m_len; i++)
      {
        cosmo_de->priv->nu_rho_s[i] = ncm_spline_cubic_notaknot_new ();
        cosmo_de->priv->nu_p_s[i]   = ncm_spline_cubic_notaknot_new ();
      }
    }
  }
}

static void
_nc_hicosmo_gcg_dispose (GObject *object)
{
  NcHICosmoGCG *cosmo_de = NC_HICOSMO_GCG (object);
  gint i;

  ncm_spline2d_clear (&cosmo_de->priv->BBN_spline2d);

  ncm_integral1d_ptr_clear (&cosmo_de->priv->nu_rho);
  ncm_integral1d_ptr_clear (&cosmo_de->priv->nu_p);

  for (i = 0; i < 10; i++)
  {
    ncm_spline_clear (&cosmo_de->priv->nu_rho_s[i]);
    ncm_spline_clear (&cosmo_de->priv->nu_p_s[i]);
  }
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_gcg_parent_class)->dispose (object);
}

static void
_nc_hicosmo_gcg_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_gcg_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_gcg_E2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_gcg_dE2_dz (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_gcg_d2E2_dz2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_gcg_H0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_gcg_Omega_t0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_gcg_Omega_c0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_gcg_T_gamma0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_gcg_Yp_4He (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_gcg_Omega_g0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_gcg_Omega_nu0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_gcg_Omega_m0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_gcg_Omega_r0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_gcg_Omega_b0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_gcg_bgp_cs2 (NcHICosmo *cosmo, gdouble z);
static guint _nc_hicosmo_gcg_NMassNu (NcHICosmo *cosmo);
static void _nc_hicosmo_gcg_MassNuInfo (NcHICosmo *cosmo, guint nu_i, gdouble *mass_eV, gdouble *T_0, gdouble *mu, gdouble *g);

static gdouble _nc_hicosmo_gcg_Omega_mnu0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_gcg_Press_mnu0 (NcHICosmo *cosmo);

static gdouble _nc_hicosmo_gcg_Omega_mnu0_n (NcHICosmo *cosmo, const guint n);
static gdouble _nc_hicosmo_gcg_Press_mnu0_n (NcHICosmo *cosmo, const guint n);

static gdouble _nc_hicosmo_gcg_E2Omega_mnu (NcHICosmo *cosmo, const gdouble z);
static gdouble _nc_hicosmo_gcg_E2Press_mnu (NcHICosmo *cosmo, const gdouble z);

static gdouble _nc_hicosmo_gcg_E2Omega_mnu_n (NcHICosmo *cosmo, const guint n, const gdouble z);
static gdouble _nc_hicosmo_gcg_E2Press_mnu_n (NcHICosmo *cosmo, const guint n, const gdouble z);

static gdouble _nc_hicosmo_gcg_E2Omega_m (NcHICosmo *cosmo, const gdouble z);
static gdouble _nc_hicosmo_gcg_E2Omega_r (NcHICosmo *cosmo, const gdouble z);

static void
nc_hicosmo_gcg_class_init (NcHICosmoGCGClass *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass *parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (klass);

  object_class->constructed = &_nc_hicosmo_gcg_constructed;
  object_class->dispose     = &_nc_hicosmo_gcg_dispose;
  object_class->finalize    = &_nc_hicosmo_gcg_finalize;

  ncm_model_class_set_name_nick (model_class, "GCG class", "NcHICosmoGCG");

  ncm_model_class_add_params (model_class,
                              NC_HICOSMO_GCG_SPARAM_LEN, NC_HICOSMO_GCG_VPARAM_LEN, PROP_SIZE);

  /* Set H_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_GCG_H0, "H_0", "H0",
                              40.0, 120.0, 1.0,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_GCG_DEFAULT_H0, NCM_PARAM_TYPE_FIXED);

  /* Set Omega_c0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_GCG_OMEGA_C, "\\Omega_{c0}", "Omegac",
                              1e-8, 0.90, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_GCG_DEFAULT_OMEGA_C,
                              NCM_PARAM_TYPE_FREE);

  /* Set Omega_x0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_GCG_OMEGA_X, "\\Omega_{x0}", "Omegax",
                              1e-8, 2.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_GCG_DEFAULT_OMEGA_X,
                              NCM_PARAM_TYPE_FREE);

  /* Set T_gamma0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_GCG_T_GAMMA0, "T_{\\gamma0}", "Tgamma0",
                              2.0, 3.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_GCG_DEFAULT_T_GAMMA0,
                              NCM_PARAM_TYPE_FIXED);

  /* Set He Yp param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_GCG_HE_YP, "Y_p", "Yp",
                              0.0, 1.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_GCG_DEFAULT_HE_YP,
                              NCM_PARAM_TYPE_FIXED);

  /* Set ENnu param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_GCG_ENNU, "N_\\nu", "ENnu",
                              0.0, 4.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_GCG_DEFAULT_ENNU,
                              NCM_PARAM_TYPE_FIXED);

  /* Set Omega_b0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_GCG_OMEGA_B, "\\Omega_{b0}", "Omegab",
                              0.03, 0.05, 5.0e-4,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_GCG_DEFAULT_OMEGA_B,
                              NCM_PARAM_TYPE_FIXED);

  /* Set Gamma param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_GCG_GAMMA, "\\gamma", "gamma",
                              -1.5, +1.5, 0.05,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_GCG_DEFAULT_GAMMA,
                              NCM_PARAM_TYPE_FIXED);

  /* Set massive neutrinos mass vector param */
  ncm_model_class_set_vparam (model_class, NC_HICOSMO_GCG_MASSNU_M, 0, "m_\\nu", "massnu",
                              0.0, 10.0, 0.01,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_GCG_DEFAULT_NU_MASS,
                              NCM_PARAM_TYPE_FIXED);

  /* Set massive neutrinos temperature vector param */
  ncm_model_class_set_vparam (model_class, NC_HICOSMO_GCG_MASSNU_T, 0, "T_{\\nu0}", "Tnu",
                              0.0, 10.0, 0.01,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_GCG_DEFAULT_NU_T,
                              NCM_PARAM_TYPE_FIXED);

  /* Set massive neutrinos relative chemical potential vector param */
  ncm_model_class_set_vparam (model_class, NC_HICOSMO_GCG_MASSNU_MU, 0, "\\mu_{\\nu}", "munu",
                              -10.0, 10.0, 0.01,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_GCG_DEFAULT_NU_MU,
                              NCM_PARAM_TYPE_FIXED);

  /* Set massive neutrinos degeneracy vector param */
  ncm_model_class_set_vparam (model_class, NC_HICOSMO_GCG_MASSNU_G, 0, "g_{\\nu}", "gnu",
                              0.0, 10.0, 0.01,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_GCG_DEFAULT_NU_G,
                              NCM_PARAM_TYPE_FIXED);

  /* Check for errors in parameters initialization */
  ncm_model_class_check_params_info (model_class);

  nc_hicosmo_set_H0_impl         (parent_class, &_nc_hicosmo_gcg_H0);
  nc_hicosmo_set_E2_impl         (parent_class, &_nc_hicosmo_gcg_E2);
  nc_hicosmo_set_Omega_c0_impl   (parent_class, &_nc_hicosmo_gcg_Omega_c0);
  nc_hicosmo_set_Omega_b0_impl   (parent_class, &_nc_hicosmo_gcg_Omega_b0);
  nc_hicosmo_set_Omega_g0_impl   (parent_class, &_nc_hicosmo_gcg_Omega_g0);
  nc_hicosmo_set_Omega_nu0_impl  (parent_class, &_nc_hicosmo_gcg_Omega_nu0);
  nc_hicosmo_set_Omega_m0_impl   (parent_class, &_nc_hicosmo_gcg_Omega_m0);
  nc_hicosmo_set_Omega_r0_impl   (parent_class, &_nc_hicosmo_gcg_Omega_r0);
  nc_hicosmo_set_Omega_t0_impl   (parent_class, &_nc_hicosmo_gcg_Omega_t0);
  nc_hicosmo_set_T_gamma0_impl   (parent_class, &_nc_hicosmo_gcg_T_gamma0);
  nc_hicosmo_set_Yp_4He_impl     (parent_class, &_nc_hicosmo_gcg_Yp_4He);
  nc_hicosmo_set_dE2_dz_impl     (parent_class, &_nc_hicosmo_gcg_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl   (parent_class, &_nc_hicosmo_gcg_d2E2_dz2);
  nc_hicosmo_set_bgp_cs2_impl    (parent_class, &_nc_hicosmo_gcg_bgp_cs2);

  /* Massive neutrino related implementations */
  nc_hicosmo_set_NMassNu_impl    (parent_class, &_nc_hicosmo_gcg_NMassNu);
  nc_hicosmo_set_MassNuInfo_impl (parent_class, &_nc_hicosmo_gcg_MassNuInfo);

  nc_hicosmo_set_Omega_mnu0_impl (parent_class, &_nc_hicosmo_gcg_Omega_mnu0);
  nc_hicosmo_set_Press_mnu0_impl (parent_class, &_nc_hicosmo_gcg_Press_mnu0);

  nc_hicosmo_set_Omega_mnu0_n_impl (parent_class, &_nc_hicosmo_gcg_Omega_mnu0_n);
  nc_hicosmo_set_Press_mnu0_n_impl (parent_class, &_nc_hicosmo_gcg_Press_mnu0_n);
  
  nc_hicosmo_set_E2Omega_mnu_impl (parent_class, &_nc_hicosmo_gcg_E2Omega_mnu);
  nc_hicosmo_set_E2Press_mnu_impl (parent_class, &_nc_hicosmo_gcg_E2Press_mnu);

  nc_hicosmo_set_E2Omega_mnu_n_impl (parent_class, &_nc_hicosmo_gcg_E2Omega_mnu_n);
  nc_hicosmo_set_E2Press_mnu_n_impl (parent_class, &_nc_hicosmo_gcg_E2Press_mnu_n);

  nc_hicosmo_set_E2Omega_m_impl   (parent_class, &_nc_hicosmo_gcg_E2Omega_m);
  nc_hicosmo_set_E2Omega_r_impl   (parent_class, &_nc_hicosmo_gcg_E2Omega_r);
}

static gdouble _nc_hicosmo_gcg_Omega_mnu0_n (NcHICosmo *cosmo, const guint n);
static gdouble _nc_hicosmo_gcg_Omega_gnu0 (NcHICosmo *cosmo);

#define VECTOR (NCM_MODEL (cosmo)->params)
#define MACRO_H0 (ncm_vector_get (VECTOR, NC_HICOSMO_GCG_H0))
#define OMEGA_C (ncm_vector_get (VECTOR, NC_HICOSMO_GCG_OMEGA_C))
#define OMEGA_X (ncm_vector_get (VECTOR, NC_HICOSMO_GCG_OMEGA_X))
#define T_GAMMA0 (ncm_vector_get (VECTOR, NC_HICOSMO_GCG_T_GAMMA0))
#define HE_YP (ncm_vector_get (VECTOR, NC_HICOSMO_GCG_HE_YP))
#define ENNU (ncm_vector_get (VECTOR, NC_HICOSMO_GCG_ENNU))
#define OMEGA_R (_nc_hicosmo_gcg_Omega_gnu0 (NC_HICOSMO (cosmo)))
#define OMEGA_B (ncm_vector_get (VECTOR, NC_HICOSMO_GCG_OMEGA_B))
#define GAMMA (ncm_vector_get (VECTOR, NC_HICOSMO_GCG_GAMMA))

#define OMEGA_M (OMEGA_B + OMEGA_C)
#define OMEGA_K (1.0 - (OMEGA_B + OMEGA_C + OMEGA_R + OMEGA_X + _nc_hicosmo_gcg_Omega_mnu0 (cosmo)))

typedef struct _NcHICosmoGCGNuInt
{
  NcHICosmoGCG *cosmo_de;
  gdouble xi2_0;
  gdouble yi;
} NcHICosmoGCGNuInt;

typedef struct _neutrino_int
{
  const gdouble xi2;
  const gdouble yi;
} neutrino_int;

static gdouble
_nc_hicosmo_gcg_nu_rho_f (gdouble z, gpointer userdata)
{
  NcHICosmoGCGNuInt *nu_int = (NcHICosmoGCGNuInt *) userdata;
  neutrino_int nudata = { nu_int->xi2_0 / gsl_pow_2 (1.0 + z), nu_int->yi };

  ncm_integral1d_ptr_set_userdata (nu_int->cosmo_de->priv->nu_rho, &nudata);

  {
    gdouble err = 0.0;
    const gdouble int_rho = ncm_integral1d_eval_gauss_laguerre (NCM_INTEGRAL1D (nu_int->cosmo_de->priv->nu_rho), &err);

    return int_rho;
  }
}

static gdouble
_nc_hicosmo_gcg_nu_p_f (gdouble z, gpointer userdata)
{
  NcHICosmoGCGNuInt *nu_int = (NcHICosmoGCGNuInt *) userdata;
  neutrino_int nudata = { nu_int->xi2_0 / gsl_pow_2 (1.0 + z), nu_int->yi };

  ncm_integral1d_ptr_set_userdata (nu_int->cosmo_de->priv->nu_p, &nudata);

  {
    gdouble err = 0.0;
    const gdouble int_p = ncm_integral1d_eval_gauss_laguerre (NCM_INTEGRAL1D (nu_int->cosmo_de->priv->nu_p), &err);

    return int_p;
  }
}

static void
_nc_hicosmo_gcg_prepare (NcHICosmoGCG *cosmo_de)
{
  NcmModel *model = NCM_MODEL (cosmo_de);
  if (!ncm_model_state_is_update (model))
  {
    const guint m_len  = ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_M);

    if (m_len > 0)
    {
      gint n;

      for (n = 0; n < m_len; n++)
      {
        const gdouble m       = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_M, n);
        const gdouble Tgamma0 = nc_hicosmo_T_gamma0 (NC_HICOSMO (cosmo_de));
        const gdouble T0      = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_T, n) * Tgamma0;
        const gdouble yi      = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_MU, n);
        const gdouble xi_0    = m * ncm_c_eV () / (ncm_c_kb () * T0);
        const gdouble xi2_0   = xi_0 * xi_0;

        NcHICosmoGCGNuInt nu_int = {cosmo_de, xi2_0, yi};
        gsl_function F;

        F.params   = &nu_int;

        F.function = &_nc_hicosmo_gcg_nu_rho_f;
        ncm_spline_set_func (cosmo_de->priv->nu_rho_s[n], NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT, &F, 0.0, cosmo_de->priv->zmax, 0, _NC_HICOSMO_GCG_MNU_PREC);

        F.function = &_nc_hicosmo_gcg_nu_p_f;
        ncm_spline_set_func (cosmo_de->priv->nu_p_s[n], NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT, &F, 0.0, cosmo_de->priv->zmax, 0, _NC_HICOSMO_GCG_MNU_PREC);
      }
    }

    ncm_model_state_set_update (model);
  }
  else
    return;
}


/****************************************************************************
 *Normalized Hubble function
 ****************************************************************************/

static gdouble
_nc_hicosmo_gcg_E2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble Omega_k         = OMEGA_K;
  const gdouble Omega_c         = OMEGA_C;
  const gdouble Omega_x         = OMEGA_X;
  const gdouble x               = 1.0 + z;
  const gdouble x2              = x * x;
  const gdouble x3              = x2 * x;
  const gdouble x4              = x3 * x;

  const gdouble gamma           = GAMMA;
  const gdouble x3_1_p_gamma    = pow (x3, 1.0 + gamma);
  const gdouble arg             = (Omega_x + Omega_c * x3_1_p_gamma) / (Omega_x + Omega_c);
  
  const gdouble E2 = 
    OMEGA_R * x4 + OMEGA_B * x3 + Omega_k * x2 +
    (Omega_x + Omega_c) * pow (arg, 1.0 / (1.0 + gamma))
    + _nc_hicosmo_gcg_E2Omega_mnu (cosmo, z);
  
  return E2;
}

/****************************************************************************
 * dE2_dz
 ****************************************************************************/

static gdouble
_nc_hicosmo_gcg_dE2_dz (NcHICosmo *cosmo, gdouble z)
{
  const gdouble Omega_k         = OMEGA_K;
  const gdouble Omega_c         = OMEGA_C;
  const gdouble Omega_x         = OMEGA_X;
  const gdouble x               = 1.0 + z;
  const gdouble x2              = x * x;
  const gdouble x3              = x2 * x;
  const gdouble dE2Omega_mnu_dz = 3.0 * (_nc_hicosmo_gcg_E2Omega_mnu (cosmo, z) + _nc_hicosmo_gcg_E2Press_mnu (cosmo, z)) / x;

  const gdouble gamma           = GAMMA;
  const gdouble x3_1_p_gamma    = pow (x3, 1.0 + gamma);
  const gdouble arg             = (Omega_x + Omega_c * x3_1_p_gamma) / (Omega_x + Omega_c);

  return (4.0 * OMEGA_R * x3 + 3.0 * OMEGA_B * x2 + 2.0 * Omega_k * x)
    + 3.0 * Omega_c * (x3_1_p_gamma / x) * pow (arg, -gamma / (1.0 + gamma))
    + dE2Omega_mnu_dz;
}

/****************************************************************************
 * d2E2_dz2
 ****************************************************************************/

static gdouble
_nc_hicosmo_gcg_d2E2_dz2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble omega_k = OMEGA_K;
  const gdouble Omega_c = OMEGA_C;
  const gdouble Omega_x = OMEGA_X;
  const gdouble x       = 1.0 + z;
  const gdouble x2      = x * x;
  const gdouble x3      = x2 * x;

  const gdouble gamma           = GAMMA;
  const gdouble x3_1_p_gamma    = pow (x3, 1.0 + gamma);
  const gdouble arg             = (Omega_x + Omega_c * x3_1_p_gamma) / (Omega_x + Omega_c);

g_assert_not_reached ();
  
  return (12.0 * OMEGA_R * x2 + 6.0 * OMEGA_B * x + 2.0 * omega_k)
    - 9.0 * (pow(x3_1_p_gamma, 2.0) / x2) * gamma * pow(Omega_c, 2.0) * pow(arg, -2.0 + 1.0 / (1.0 + gamma)) / (Omega_c + Omega_x)
    + 3.0 * (x3_1_p_gamma / x2) * (2.0 + 3.0 * gamma) * Omega_c * pow(arg, -gamma / (1.0 + gamma));
}

/****************************************************************************
 * Simple functions
 ****************************************************************************/
static gdouble
_nc_hicosmo_gcg_H0 (NcHICosmo *cosmo)
{
  return MACRO_H0;
}

static gdouble
_nc_hicosmo_gcg_Omega_t0 (NcHICosmo *cosmo)
{
  return OMEGA_M + OMEGA_X + OMEGA_R + _nc_hicosmo_gcg_Omega_mnu0 (cosmo);
}

static gdouble
_nc_hicosmo_gcg_Omega_c0 (NcHICosmo *cosmo)
{
  return OMEGA_C;
}

static gdouble
_nc_hicosmo_gcg_T_gamma0 (NcHICosmo *cosmo)
{
  return T_GAMMA0;
}

static gdouble
_nc_hicosmo_gcg_Yp_4He (NcHICosmo *cosmo)
{
  NcmModel *model = NCM_MODEL (cosmo);
  NcHICosmoGCG *cosmo_de = NC_HICOSMO_GCG (cosmo);

  if (ncm_model_param_get_ftype (model, NC_HICOSMO_GCG_HE_YP) == NCM_PARAM_TYPE_FIXED)
  {
    const gdouble wb    = nc_hicosmo_Omega_b0h2 (cosmo);
    const gdouble DNeff = nc_hicosmo_Neff (cosmo) - 3.046;

    if (FALSE)
    {
      const gdouble wb2 = wb * wb;
      const gdouble DNeff2 = DNeff * DNeff;
      const gdouble Yp = 0.2311 + 0.9502 * wb - 11.27 * wb2 + DNeff * (0.01356 + 0.008581 * wb - 0.1810 * wb2) + DNeff2 * (-0.0009795 - 0.001370 * wb + 0.01746 * wb2);
      return Yp;
    }

    else
    {
      if (model->pkey != cosmo_de->priv->HE4_Yp_key)
      {
        const gdouble Yp = ncm_spline2d_eval (NC_HICOSMO_GCG (cosmo)->priv->BBN_spline2d,
                                              wb, DNeff);
        ncm_model_orig_param_set (model, NC_HICOSMO_GCG_HE_YP, Yp);
        cosmo_de->priv->HE4_Yp_key = model->pkey;
        /*printf ("# omega_b % 20.15g DeltaNnu % 20.15g Yp % 20.15g\n",  wb, DNeff, Yp); */
      }
    }
  }

  return HE_YP;
}

static gdouble
_nc_hicosmo_gcg_Omega_g0 (NcHICosmo *cosmo)
{
  const gdouble h = MACRO_H0 / 100.0;
  const gdouble h2 = h * h;
  return ncm_c_radiation_temp_to_h2omega_r (T_GAMMA0) / h2;
}

static gdouble
_nc_hicosmo_gcg_Omega_nu0 (NcHICosmo *cosmo)
{
  const gdouble conv = 7.0 / 8.0 * pow (4.0 / 11.0, 4.0 / 3.0);
  return ENNU * conv * _nc_hicosmo_gcg_Omega_g0 (cosmo);
}

static gdouble
_nc_hicosmo_gcg_Omega_gnu0 (NcHICosmo *cosmo)
{
  const gdouble conv        = 7.0 / 8.0 * pow (4.0 / 11.0, 4.0 / 3.0);
  
  return (1.0 + ENNU * conv) * _nc_hicosmo_gcg_Omega_g0 (cosmo);
}

static gdouble
_nc_hicosmo_gcg_Omega_r0 (NcHICosmo *cosmo)
{
  const gdouble Omega_mnu_r = 3.0 * _nc_hicosmo_gcg_Press_mnu0 (cosmo);  
  return _nc_hicosmo_gcg_Omega_gnu0 (cosmo) + Omega_mnu_r;
}

static gdouble
_nc_hicosmo_gcg_E2Omega_r (NcHICosmo *cosmo, const gdouble z)
{
  const gdouble x4          = gsl_pow_4 (1.0 + z);
  const gdouble conv        = 7.0 / 8.0 * pow (4.0 / 11.0, 4.0 / 3.0);
  const gdouble Omega_mnu_r = 3.0 * _nc_hicosmo_gcg_E2Press_mnu (cosmo, z);
  
  return (1.0 + ENNU * conv) * _nc_hicosmo_gcg_Omega_g0 (cosmo) * x4 + Omega_mnu_r;
}

static gdouble
_nc_hicosmo_gcg_Omega_m0 (NcHICosmo *cosmo)
{
  const gdouble Omega_mnu_r = 3.0 * _nc_hicosmo_gcg_Press_mnu0 (cosmo);
  const gdouble Omega_mnu_d = _nc_hicosmo_gcg_Omega_mnu0 (cosmo) - Omega_mnu_r;
  
  return OMEGA_M + Omega_mnu_d;
}

static gdouble
_nc_hicosmo_gcg_E2Omega_m (NcHICosmo *cosmo, const gdouble z)
{
  const gdouble x3          = gsl_pow_3 (1.0 + z);
  const gdouble Omega_mnu_r = 3.0 * _nc_hicosmo_gcg_E2Press_mnu (cosmo, z);
  const gdouble Omega_mnu_d = _nc_hicosmo_gcg_E2Omega_mnu (cosmo, z) - Omega_mnu_r;
  const gdouble Omega_c     = OMEGA_C;
  const gdouble Omega_x     = OMEGA_X;
  const gdouble gamma       = GAMMA;

  return OMEGA_B * x3 + OMEGA_C * x3 * pow ((Omega_c + Omega_x * pow (x3, (1.0 + gamma)))/ (Omega_c + Omega_x), -gamma / (1.0 + gamma)) + Omega_mnu_d;
}

static gdouble
_nc_hicosmo_gcg_Omega_b0 (NcHICosmo *cosmo)
{
  return OMEGA_B;
}

static gdouble
_nc_hicosmo_gcg_bgp_cs2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble x        = 1.0 + z;
  const gdouble Omega_g0 = _nc_hicosmo_gcg_Omega_g0 (cosmo);
  const gdouble Omega_b0 = _nc_hicosmo_gcg_Omega_b0 (cosmo);
  const gdouble nine_4   = 9.0 / 4.0;
  
  return 1.0 / (3.0 + nine_4 * Omega_b0 / (Omega_g0 * x));
}

static guint
_nc_hicosmo_gcg_NMassNu (NcHICosmo *cosmo)
{
  return ncm_model_vparam_len (NCM_MODEL (cosmo), NC_HICOSMO_GCG_MASSNU_M);
}

static void
_nc_hicosmo_gcg_MassNuInfo (NcHICosmo *cosmo, guint nu_i, gdouble *mass_eV, gdouble *T_0, gdouble *mu, gdouble *g)
{
  NcmModel *model = NCM_MODEL (cosmo);
  g_assert_cmpuint (nu_i, <, ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_T));
  mass_eV[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_M, nu_i);

  if (ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_T) == 1)
  {
    T_0[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_T, 0);
  }
  else
  {
    T_0[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_T, nu_i);
  }

  if (ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_MU) == 1)
  {
    mu[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_MU, 0);
  }
  else
  {
    mu[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_MU, nu_i);
  }

  if (ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_G) == 1)
  {
    g[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_G, 0);
  }
  else
  {
    g[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_G, nu_i);
  }
}

static gdouble
_nc_hicosmo_gcg_neutrino_rho_integrand (gpointer userdata, const gdouble u, const gdouble w)
{
  /* see equations 4.11 to 4.16 of Peter & Uzan, Primordial Cosmology with change of variable v = u - x_i and Gauss-Laguerre integration */
  neutrino_int *nudata  = (neutrino_int *) userdata;
  const gdouble u2      = u * u;
  const gdouble exp_yi  = exp (nudata->yi);
  const gdouble exp_myi = 1.0 / exp_yi;

  return u2 * sqrt (u2 + nudata->xi2) * (1.0 / (w + exp_myi) + 1.0 / (w + exp_yi));
}

static gdouble
_nc_hicosmo_gcg_neutrino_p_integrand (gpointer userdata, const gdouble u, const gdouble w)
{
  /* see equations 4.11 to 4.16 of Peter & Uzan, Primordial Cosmology with change of variable v = u - x_i and Gauss-Laguerre integration */
  neutrino_int *nudata  = (neutrino_int *) userdata;
  const gdouble u2      = u * u;
  const gdouble u4      = u2 * u2;
  const gdouble exp_yi  = exp (nudata->yi);
  const gdouble exp_myi = 1.0 / exp_yi;

  return u4 / (3.0 * sqrt (u2 + nudata->xi2)) * (1.0 / (w + exp_myi) + 1.0 / (w + exp_yi));
}

typedef struct _NcHICosmoGCGNeutrinoparams
{
  NcHICosmoGCG *cosmo_de;
  NcmIntegral1dPtr *int1d_ptr;
} NcHICosmoGCGNeutrinoparams;

static gdouble
_nc_hicosmo_gcg_Omega_mnu0_n (NcHICosmo *cosmo, const guint n)
{
  return _nc_hicosmo_gcg_E2Omega_mnu_n (cosmo, n, 0.0);
}

static gdouble
_nc_hicosmo_gcg_Press_mnu0_n (NcHICosmo *cosmo, const guint n)
{
  return _nc_hicosmo_gcg_E2Press_mnu_n (cosmo, n, 0.0);
}

static gdouble
_nc_hicosmo_gcg_Omega_mnu0 (NcHICosmo *cosmo)
{
  return _nc_hicosmo_gcg_E2Omega_mnu (cosmo, 0.0);
}

static gdouble
_nc_hicosmo_gcg_Press_mnu0 (NcHICosmo *cosmo)
{
  return _nc_hicosmo_gcg_E2Press_mnu (cosmo, 0.0);
}

static gdouble
_nc_hicosmo_gcg_E2Omega_mnu_n (NcHICosmo *cosmo, const guint n, const gdouble z)
{
  NcHICosmoGCG *cosmo_de = NC_HICOSMO_GCG (cosmo);
  NcmModel *model       = NCM_MODEL (cosmo);
  const gdouble h2      = nc_hicosmo_h2 (cosmo);
  const gdouble ffac    = 15.0 / (2.0 * gsl_pow_4 (M_PI)); 
  const gdouble Tgamma  = (1.0 + z) * nc_hicosmo_T_gamma0 (cosmo);
  const gdouble T       = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_T, n) * Tgamma;
  const gdouble g       = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_G, n);

  _nc_hicosmo_gcg_prepare (NC_HICOSMO_GCG (cosmo_de));

  {
    const gdouble int_Ef       = ncm_spline_eval (cosmo_de->priv->nu_rho_s[n], z > cosmo_de->priv->zmax ? cosmo_de->priv->zmax : z);
    const gdouble Omega_mnu0_n = ffac * g * gsl_pow_4 (T) * ncm_c_blackbody_per_crit_density_h2 () / h2 * int_Ef;
    
    return Omega_mnu0_n;
  }
}

static gdouble
_nc_hicosmo_gcg_E2Press_mnu_n (NcHICosmo *cosmo, const guint n, const gdouble z)
{
  NcHICosmoGCG *cosmo_de = NC_HICOSMO_GCG (cosmo);
  NcmModel *model       = NCM_MODEL (cosmo);
  const gdouble h2      = nc_hicosmo_h2 (cosmo);
  const gdouble ffac    = 15.0 / (2.0 * gsl_pow_4 (M_PI)); 
  const gdouble Tgamma  = (1.0 + z) * nc_hicosmo_T_gamma0 (cosmo);
  const gdouble T       = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_T, n) * Tgamma;
  const gdouble g       = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_G, n);

  _nc_hicosmo_gcg_prepare (NC_HICOSMO_GCG (cosmo_de));

  {
    const gdouble int_pf       = ncm_spline_eval (cosmo_de->priv->nu_p_s[n], z > cosmo_de->priv->zmax ? cosmo_de->priv->zmax : z);
    const gdouble Press_mnu0_n = ffac * g * gsl_pow_4 (T) * ncm_c_blackbody_per_crit_density_h2 () / h2 * int_pf;

    return Press_mnu0_n;
  }
}

static gdouble
_nc_hicosmo_gcg_E2Omega_mnu (NcHICosmo *cosmo, const gdouble z)
{
  NcHICosmoGCG *cosmo_de = NC_HICOSMO_GCG (cosmo);
  NcmModel *model       = NCM_MODEL (cosmo);
  const guint nmassnu   = ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_M);
  const gdouble h2      = nc_hicosmo_h2 (cosmo);
  const gdouble ffac    = 15.0 / (2.0 * gsl_pow_4 (M_PI)); 
  const gdouble Tgamma  = (1.0 + z) * nc_hicosmo_T_gamma0 (cosmo);
  gdouble Omega_mnu0    = 0.0;

  guint n;

  _nc_hicosmo_gcg_prepare (NC_HICOSMO_GCG (cosmo_de));
  
  for (n = 0; n < nmassnu; n++)
  {
    const gdouble T            = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_T, n) * Tgamma;
    const gdouble g            = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_G, n);
    const gdouble int_Ef       = ncm_spline_eval (cosmo_de->priv->nu_rho_s[n], z > cosmo_de->priv->zmax ? cosmo_de->priv->zmax : z);
    const gdouble Omega_mnu0_n = ffac * g * gsl_pow_4 (T) * ncm_c_blackbody_per_crit_density_h2 () / h2 * int_Ef;

    Omega_mnu0 += Omega_mnu0_n;
  }

  return Omega_mnu0;
}

static gdouble
_nc_hicosmo_gcg_E2Press_mnu (NcHICosmo *cosmo, const gdouble z)
{
  NcHICosmoGCG *cosmo_de = NC_HICOSMO_GCG (cosmo);
  NcmModel *model       = NCM_MODEL (cosmo);
  const guint nmassnu   = ncm_model_vparam_len (model, NC_HICOSMO_GCG_MASSNU_M);
  const gdouble Tgamma  = (1.0 + z) * nc_hicosmo_T_gamma0 (cosmo);
  const gdouble h2      = nc_hicosmo_h2 (cosmo);
  const gdouble ffac    = 15.0 / (2.0 * gsl_pow_4 (M_PI)); 
  gdouble Press_mnu0    = 0.0;
  guint n;

  _nc_hicosmo_gcg_prepare (NC_HICOSMO_GCG (cosmo_de));
  
  for (n = 0; n < nmassnu; n++)
  {
    const gdouble T            = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_T, n) * Tgamma;
    const gdouble g            = ncm_model_orig_vparam_get (model, NC_HICOSMO_GCG_MASSNU_G, n);
    const gdouble int_pf       = ncm_spline_eval (cosmo_de->priv->nu_p_s[n], z > cosmo_de->priv->zmax ? cosmo_de->priv->zmax : z);
    const gdouble Press_mnu0_n = ffac * g * gsl_pow_4 (T) * ncm_c_blackbody_per_crit_density_h2 () / h2 * int_pf;

    Press_mnu0 += Press_mnu0_n;
  }
  return Press_mnu0;
}

/**
 * nc_hicosmo_gcg_omega_x2omega_k:
 * @cosmo_gcg: a #NcHICosmoGCG
 *
 * FIXME
 *
 */
void
nc_hicosmo_gcg_omega_x2omega_k (NcHICosmoGCG *cosmo_de)
{
  NcHICosmoGCGReparamOk *gcg_reparam_ok = nc_hicosmo_gcg_reparam_ok_new (ncm_model_len (NCM_MODEL (cosmo_de)));
  ncm_model_set_reparam (NCM_MODEL (cosmo_de), NCM_REPARAM (gcg_reparam_ok));
  return;
}

/**
 * nc_hicosmo_gcg_cmb_params:
 * @cosmo_gcg: a #NcHICosmoGCG
 *
 * FIXME
 *
 */
void
nc_hicosmo_gcg_cmb_params (NcHICosmoGCG *cosmo_de)
{
  NcHICosmoGCGReparamCMB *gcg_reparam_cmb = nc_hicosmo_gcg_reparam_cmb_new (ncm_model_len (NCM_MODEL (cosmo_de)));
  ncm_model_set_reparam (NCM_MODEL (cosmo_de), NCM_REPARAM (gcg_reparam_cmb));
  return;
}

/***********************************************************************/
/* Reparam CMB                                                         */
/***********************************************************************/

G_DEFINE_TYPE (NcHICosmoGCGReparamCMB, nc_hicosmo_gcg_reparam_cmb, NCM_TYPE_REPARAM);

static void
nc_hicosmo_gcg_reparam_cmb_init (NcHICosmoGCGReparamCMB *gcg_reparam_ok)
{
}

static void
nc_hicosmo_gcg_reparam_cmb_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hicosmo_gcg_reparam_cmb_parent_class)->constructed (object);
  {
    NcHICosmoGCGReparamCMB *reparam_cmb = NC_HICOSMO_GCG_REPARAM_CMB (object);

    ncm_reparam_set_param_desc_full (NCM_REPARAM (reparam_cmb), NC_HICOSMO_GCG_OMEGA_C,
                                     "omegac","\\omega_{c0}", 1.0e-8, 5.0e-1, 1.0e-2,
                                     NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, 0.0, NCM_PARAM_TYPE_FIXED);

    ncm_reparam_set_param_desc_full (NCM_REPARAM (reparam_cmb), NC_HICOSMO_GCG_OMEGA_B,
                                     "omegab","\\omega_{b0}", 1.0e-8, 5.0e-1, 1.0e-2,
                                     NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, 0.0, NCM_PARAM_TYPE_FIXED);

    ncm_reparam_set_param_desc_full (NCM_REPARAM (reparam_cmb), NC_HICOSMO_GCG_OMEGA_X,
                                     "Omegak","\\omega_{k0}", -5.0e-1, 5.0e-1, 1.0e-2,
                                     NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, 0.0, NCM_PARAM_TYPE_FIXED);

    NCM_REPARAM (reparam_cmb)->compat_type = NC_TYPE_HICOSMO_GCG;
  }
}

static void
nc_hicosmo_gcg_reparam_cmb_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_gcg_reparam_cmb_parent_class)->finalize (object);
}

static gboolean _nc_hicosmo_gcg_reparam_cmb_old2new (NcmReparam *reparam, NcmModel *model);
static gboolean _nc_hicosmo_gcg_reparam_cmb_new2old (NcmReparam *reparam, NcmModel *model);
static gboolean _nc_hicosmo_gcg_reparam_cmb_jac (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac);

static void
nc_hicosmo_gcg_reparam_cmb_class_init (NcHICosmoGCGReparamCMBClass *klass)
{
  GObjectClass* object_class     = G_OBJECT_CLASS (klass);
  NcmReparamClass *reparam_class = NCM_REPARAM_CLASS (klass);

  object_class->constructed = nc_hicosmo_gcg_reparam_cmb_constructed;
  object_class->finalize    = nc_hicosmo_gcg_reparam_cmb_finalize;

  reparam_class->old2new = &_nc_hicosmo_gcg_reparam_cmb_old2new;
  reparam_class->new2old = &_nc_hicosmo_gcg_reparam_cmb_new2old;
  reparam_class->jac     = &_nc_hicosmo_gcg_reparam_cmb_jac;
}

static gboolean
_nc_hicosmo_gcg_reparam_cmb_old2new (NcmReparam *reparam, NcmModel *model)
{
  NcHICosmo *cosmo  = NC_HICOSMO (model);
  NcmVector *params = ncm_model_orig_params_peek_vector (model);

  const gdouble Omega_k0 = nc_hicosmo_Omega_k0 (cosmo);
  const gdouble omega_c0 = nc_hicosmo_Omega_c0h2 (cosmo);
  const gdouble omega_b0 = nc_hicosmo_Omega_b0h2 (cosmo);

  ncm_vector_memcpy (reparam->new_params, params);

  ncm_vector_set (reparam->new_params, NC_HICOSMO_GCG_OMEGA_B, omega_b0);
  ncm_vector_set (reparam->new_params, NC_HICOSMO_GCG_OMEGA_C, omega_c0);
  ncm_vector_set (reparam->new_params, NC_HICOSMO_GCG_OMEGA_X, Omega_k0);
  
  return TRUE;
}

static gboolean
_nc_hicosmo_gcg_reparam_cmb_new2old (NcmReparam *reparam, NcmModel *model)
{
  NcmVector *params = ncm_model_orig_params_peek_vector (model);
  ncm_vector_memcpy (params, reparam->new_params);

  {
    NcHICosmo *cosmo       = NC_HICOSMO (model);
    const gdouble h2       = nc_hicosmo_h2 (cosmo);
    const gdouble Omega_c0 = ncm_vector_get (reparam->new_params, NC_HICOSMO_GCG_OMEGA_C) / h2;
    const gdouble Omega_b0 = ncm_vector_get (reparam->new_params, NC_HICOSMO_GCG_OMEGA_B) / h2;
    const gdouble Omega_k0 = ncm_vector_get (reparam->new_params, NC_HICOSMO_GCG_OMEGA_X);

    ncm_vector_set (params, NC_HICOSMO_GCG_OMEGA_C, Omega_c0);
    ncm_vector_set (params, NC_HICOSMO_GCG_OMEGA_B, Omega_b0);

    {
      const gdouble Omega_r0 = nc_hicosmo_Omega_r0 (cosmo);
      const gdouble Omega_m0 = nc_hicosmo_Omega_m0 (cosmo);
      const gdouble Omega_x0 = 1.0 - (Omega_m0 + Omega_r0 + Omega_k0);

      ncm_vector_set (params, NC_HICOSMO_GCG_OMEGA_X, Omega_x0);
    }
  }

  return TRUE;
}

static gboolean
_nc_hicosmo_gcg_reparam_cmb_jac (NcmReparam *reparam, NcmModel *model, NcmMatrix *jac)
{
  g_assert_not_reached ();
}

/**
 * nc_hicosmo_gcg_reparam_cmb_new: (constructor)
 * @length: number of parameters
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcHICosmoGCGReparamCMB
 */
NcHICosmoGCGReparamCMB *
nc_hicosmo_gcg_reparam_cmb_new (guint length)
{
  NcHICosmoGCGReparamCMB *gcg_reparam_cmb = g_object_new (NC_TYPE_HICOSMO_GCG_REPARAM_CMB,
                                                        "length", length,
                                                        NULL);
  return gcg_reparam_cmb;
}

/***********************************************************************/
/* Reparam Omega_x -> Omega_k                                          */
/***********************************************************************/

G_DEFINE_TYPE (NcHICosmoGCGReparamOk, nc_hicosmo_gcg_reparam_ok, NCM_TYPE_REPARAM);

static void
nc_hicosmo_gcg_reparam_ok_init (NcHICosmoGCGReparamOk *gcg_reparam_ok)
{
}

static void
nc_hicosmo_gcg_reparam_ok_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hicosmo_gcg_reparam_ok_parent_class)->constructed (object);
  {
    NcHICosmoGCGReparamOk *reparam_Ok = NC_HICOSMO_GCG_REPARAM_OK (object);
    ncm_reparam_set_param_desc_full (NCM_REPARAM (reparam_Ok), NC_HICOSMO_GCG_OMEGA_X,
                                     "Omegak","\\Omega_{k0}", -5.0e-1, 5.0e-1, 1.0e-2,
                                     NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, 0.0, NCM_PARAM_TYPE_FIXED);
    NCM_REPARAM (reparam_Ok)->compat_type = NC_TYPE_HICOSMO_GCG;
  }
}

static void
nc_hicosmo_gcg_reparam_ok_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_gcg_reparam_ok_parent_class)->finalize (object);
}

static gboolean _nc_hicosmo_gcg_reparam_ok_old2new (NcmReparam *reparam, NcmModel *model);
static gboolean _nc_hicosmo_gcg_reparam_ok_new2old (NcmReparam *reparam, NcmModel *model);
static gboolean _nc_hicosmo_gcg_reparam_ok_jac (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac);

static void
nc_hicosmo_gcg_reparam_ok_class_init (NcHICosmoGCGReparamOkClass *klass)
{
  GObjectClass* object_class     = G_OBJECT_CLASS (klass);
  NcmReparamClass *reparam_class = NCM_REPARAM_CLASS (klass);

  object_class->constructed = nc_hicosmo_gcg_reparam_ok_constructed;
  object_class->finalize    = nc_hicosmo_gcg_reparam_ok_finalize;

  reparam_class->old2new = &_nc_hicosmo_gcg_reparam_ok_old2new;
  reparam_class->new2old = &_nc_hicosmo_gcg_reparam_ok_new2old;
  reparam_class->jac     = &_nc_hicosmo_gcg_reparam_ok_jac;
}

static gboolean
_nc_hicosmo_gcg_reparam_ok_old2new (NcmReparam *reparam, NcmModel *model)
{
  NcmVector *params = ncm_model_orig_params_peek_vector (model);
  const gdouble Omega_k0 = nc_hicosmo_Omega_k0 (NC_HICOSMO (model));
  ncm_vector_memcpy (reparam->new_params, params);
  ncm_vector_set (reparam->new_params, NC_HICOSMO_GCG_OMEGA_X, Omega_k0);
  return TRUE;
}

static gboolean
_nc_hicosmo_gcg_reparam_ok_new2old (NcmReparam *reparam, NcmModel *model)
{
  NcmVector *params = ncm_model_orig_params_peek_vector (model);
  ncm_vector_memcpy (params, reparam->new_params);
  {
    NcHICosmo *cosmo = NC_HICOSMO (model);
    const gdouble Omega_x0 = 1.0 - (nc_hicosmo_Omega_m0 (cosmo) + nc_hicosmo_Omega_r0 (cosmo) + ncm_vector_get (reparam->new_params, NC_HICOSMO_GCG_OMEGA_X));
    ncm_vector_set (params, NC_HICOSMO_GCG_OMEGA_X, Omega_x0);
  }
  return TRUE;
}

static gboolean
_nc_hicosmo_gcg_reparam_ok_jac (NcmReparam *reparam, NcmModel *model, NcmMatrix *jac)
{
  g_assert_not_reached ();
}

/**
 * nc_hicosmo_gcg_reparam_ok_new: (constructor)
 * @length: number of parameters
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcHICosmoGCGReparamOk
 */
NcHICosmoGCGReparamOk *
nc_hicosmo_gcg_reparam_ok_new (guint length)
{
  NcHICosmoGCGReparamOk *gcg_reparam_ok = g_object_new (NC_TYPE_HICOSMO_GCG_REPARAM_OK,
                                                      "length", length,
                                                      NULL);
  return gcg_reparam_ok;
}
