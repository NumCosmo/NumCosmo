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
 * @short_description: Abstract class for implementing dark energy models
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
#include "model/nc_hicosmo_de.h"
#include "model/nc_hicosmo_de_reparam_cmb.h"
#include "model/nc_hicosmo_de_reparam_ok.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_min.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcHICosmoDEPrivate
{
  NcmSpline2d *BBN_spline2d;
  guint64 HE4_Yp_key;
  NcmIntegral1dPtr *nu_rho;
  NcmIntegral1dPtr *nu_p;
  NcmSpline *nu_rho_s[10];
  NcmSpline *nu_p_s[10];
  gdouble zmax;
  gsl_min_fminimizer *min;
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcHICosmoDE, nc_hicosmo_de, NC_TYPE_HICOSMO);

enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_de_init (NcHICosmoDE *cosmo_de)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *filename   = ncm_cfg_get_data_filename ("BBN_2017_spline2d.obj", TRUE);
  gint i;

  cosmo_de->priv               = nc_hicosmo_de_get_instance_private (cosmo_de);
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

  cosmo_de->priv->min = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);
}

static gdouble _nc_hicosmo_de_neutrino_rho_integrand (gpointer userdata, const gdouble v, const gdouble w);
static gdouble _nc_hicosmo_de_neutrino_p_integrand (gpointer userdata, const gdouble v, const gdouble w);

#define _NC_HICOSMO_DE_MNU_PREC (1.0e-6)

static void
_nc_hicosmo_de_constructed (GObject *object)
{
  /* Before model's construct!!!! */
  /* 
   * If other massive neutrinos parameters are not set, 
   * set it here before constructing the model. 
   */
  {
    NcmModel *model = NCM_MODEL (object);
    const guint nmassnu = ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_M);
  
    if (nmassnu > 0)
    {
      if (ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_T) == 0)
        g_array_index (model->vparam_len, guint, NC_HICOSMO_DE_MASSNU_T) = nmassnu;

      if (ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_MU) == 0)
        g_array_index (model->vparam_len, guint, NC_HICOSMO_DE_MASSNU_MU) = nmassnu;
      
      if (ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_G) == 0)
        g_array_index (model->vparam_len, guint, NC_HICOSMO_DE_MASSNU_G) = nmassnu;      
    }
  }
  /* Chain up : start */
  G_OBJECT_CLASS (nc_hicosmo_de_parent_class)->constructed (object);
  {
    NcmModel *model    = NCM_MODEL (object);
    const guint m_len  = ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_M);
    const guint T_len  = ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_T);
    const guint mu_len = ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_MU);
    const guint g_len  = ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_G);

    if (! ( (m_len == T_len) || (m_len > 0 && T_len == 1)))
    {
      g_error ("NcHICosmoDE: number of neutrinos masses must match the number of massive neutrino temperatures,\n"
               " or the neutrino temperature vector must be of size one to use the same value for all massive neutrinos.");
    }

    if (! ( (m_len == mu_len) || (m_len > 0 && mu_len == 1)))
    {
      g_error ("NcHICosmoDE: number of neutrinos masses must match the number of massive neutrino relative chemical potential,\n"
               " or the neutrino relative chemical potential vector must be of size one to use the same value for all massive neutrinos.");
    }

    if (! ( (m_len == g_len) || (m_len > 0 && g_len == 1)))
    {
      g_error ("NcHICosmoDE: number of neutrinos masses must match the number of massive neutrino degeneracy,\n"
               " or the neutrino degeneracy vector must be of size one to use the same value for all massive neutrinos.");
    }

    if (m_len != 0)
    {
      NcHICosmoDE *cosmo_de  = NC_HICOSMO_DE (model);
      gint i;
      
      cosmo_de->priv->nu_rho = ncm_integral1d_ptr_new (&_nc_hicosmo_de_neutrino_rho_integrand, NULL);
      cosmo_de->priv->nu_p   = ncm_integral1d_ptr_new (&_nc_hicosmo_de_neutrino_p_integrand, NULL);

      ncm_integral1d_set_reltol (NCM_INTEGRAL1D (cosmo_de->priv->nu_rho), _NC_HICOSMO_DE_MNU_PREC);
      ncm_integral1d_set_reltol (NCM_INTEGRAL1D (cosmo_de->priv->nu_p),   _NC_HICOSMO_DE_MNU_PREC);

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
_nc_hicosmo_de_dispose (GObject *object)
{
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (object);
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
  G_OBJECT_CLASS (nc_hicosmo_de_parent_class)->dispose (object);
}

static void
_nc_hicosmo_de_finalize (GObject *object)
{
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (object);

  g_clear_pointer (&cosmo_de->priv->min, gsl_min_fminimizer_free);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hicosmo_de_parent_class)->finalize (object);
}

static gdouble _nc_hicosmo_de_E2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_de_dE2_dz (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_de_d2E2_dz2 (NcHICosmo *cosmo, gdouble z);
static gdouble _nc_hicosmo_de_H0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_de_Omega_t0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_de_Omega_c0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_de_T_gamma0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_de_Yp_4He (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_de_Omega_g0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_de_Omega_nu0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_de_Omega_m0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_de_Omega_r0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_de_Omega_b0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_de_bgp_cs2 (NcHICosmo *cosmo, gdouble z);
static guint _nc_hicosmo_de_NMassNu (NcHICosmo *cosmo);
static void _nc_hicosmo_de_MassNuInfo (NcHICosmo *cosmo, guint nu_i, gdouble *mass_eV, gdouble *T_0, gdouble *mu, gdouble *g);

static gdouble _nc_hicosmo_de_Omega_mnu0 (NcHICosmo *cosmo);
static gdouble _nc_hicosmo_de_Press_mnu0 (NcHICosmo *cosmo);

static gdouble _nc_hicosmo_de_Omega_mnu0_n (NcHICosmo *cosmo, const guint n);
static gdouble _nc_hicosmo_de_Press_mnu0_n (NcHICosmo *cosmo, const guint n);

static gdouble _nc_hicosmo_de_E2Omega_mnu (NcHICosmo *cosmo, const gdouble z);
static gdouble _nc_hicosmo_de_E2Press_mnu (NcHICosmo *cosmo, const gdouble z);

static gdouble _nc_hicosmo_de_E2Omega_mnu_n (NcHICosmo *cosmo, const guint n, const gdouble z);
static gdouble _nc_hicosmo_de_E2Press_mnu_n (NcHICosmo *cosmo, const guint n, const gdouble z);

static gdouble _nc_hicosmo_de_E2Omega_m (NcHICosmo *cosmo, const gdouble z);
static gdouble _nc_hicosmo_de_E2Omega_r (NcHICosmo *cosmo, const gdouble z);

static gdouble _nc_hicosmo_de_E2Omega_de (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_dE2Omega_de_dz (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_d2E2Omega_de_dz2 (NcHICosmoDE *cosmo_de, gdouble z);
static gdouble _nc_hicosmo_de_w_de (NcHICosmoDE *cosmo_de, gdouble z);

static void _nc_hicosmo_de_get_bg_var (NcHICosmo *cosmo, const gdouble t, NcHIPertBGVar *bg_var);

static gboolean _nc_hicosmo_de_valid (NcmModel *model);

static void
nc_hicosmo_de_class_init (NcHICosmoDEClass *klass)
{
  GObjectClass *object_class   = G_OBJECT_CLASS (klass);
  NcHICosmoClass *parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class   = NCM_MODEL_CLASS (klass);

  object_class->constructed = &_nc_hicosmo_de_constructed;
  object_class->dispose     = &_nc_hicosmo_de_dispose;
  object_class->finalize    = &_nc_hicosmo_de_finalize;

  ncm_model_class_set_name_nick (model_class, "Darkenergy models abstract class", "NcHICosmoDE");

  ncm_model_class_add_params (model_class,
                              NC_HICOSMO_DE_SPARAM_LEN, NC_HICOSMO_DE_VPARAM_LEN, PROP_SIZE);

  /* Set H_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_H0, "H_0", "H0",
                              31.0, 99.0, 1.0,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_H0, NCM_PARAM_TYPE_FIXED);

  /* Set Omega_c0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_C, "\\Omega_{c0}", "Omegac",
                              0.01, 0.9, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_OMEGA_C,
                              NCM_PARAM_TYPE_FREE);

  /* Set Omega_x0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_X, "\\Omega_{x0}", "Omegax",
                              0.01, 1.4, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_OMEGA_X,
                              NCM_PARAM_TYPE_FREE);

  /* Set T_gamma0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_T_GAMMA0, "T_{\\gamma0}", "Tgamma0",
                              2.0, 3.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_T_GAMMA0,
                              NCM_PARAM_TYPE_FIXED);

  /* Set He Yp param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_HE_YP, "Y_p", "Yp",
                              0.0, 1.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_HE_YP,
                              NCM_PARAM_TYPE_FIXED);

  /* Set ENnu param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_ENNU, "N_\\nu", "ENnu",
                              0.0, 4.0, 1.0e-2,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_ENNU,
                              NCM_PARAM_TYPE_FIXED);

  /* Set Omega_b0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_B, "\\Omega_{b0}", "Omegab",
                              0.03, 0.05, 5.0e-4,
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
  ncm_model_class_set_vparam (model_class, NC_HICOSMO_DE_MASSNU_MU, 0, "\\mu_{\\nu}", "munu",
                              -10.0, 10.0, 0.01,
                              NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_NU_MU,
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
  nc_hicosmo_set_Omega_b0_impl   (parent_class, &_nc_hicosmo_de_Omega_b0);
  nc_hicosmo_set_Omega_g0_impl   (parent_class, &_nc_hicosmo_de_Omega_g0);
  nc_hicosmo_set_Omega_nu0_impl  (parent_class, &_nc_hicosmo_de_Omega_nu0);
  nc_hicosmo_set_Omega_m0_impl   (parent_class, &_nc_hicosmo_de_Omega_m0);
  nc_hicosmo_set_Omega_r0_impl   (parent_class, &_nc_hicosmo_de_Omega_r0);
  nc_hicosmo_set_Omega_t0_impl   (parent_class, &_nc_hicosmo_de_Omega_t0);
  nc_hicosmo_set_T_gamma0_impl   (parent_class, &_nc_hicosmo_de_T_gamma0);
  nc_hicosmo_set_Yp_4He_impl     (parent_class, &_nc_hicosmo_de_Yp_4He);
  nc_hicosmo_set_dE2_dz_impl     (parent_class, &_nc_hicosmo_de_dE2_dz);
  nc_hicosmo_set_d2E2_dz2_impl   (parent_class, &_nc_hicosmo_de_d2E2_dz2);
  nc_hicosmo_set_bgp_cs2_impl    (parent_class, &_nc_hicosmo_de_bgp_cs2);

  /* Massive neutrino related implementations */
  nc_hicosmo_set_NMassNu_impl    (parent_class, &_nc_hicosmo_de_NMassNu);
  nc_hicosmo_set_MassNuInfo_impl (parent_class, &_nc_hicosmo_de_MassNuInfo);

  nc_hicosmo_set_Omega_mnu0_impl (parent_class, &_nc_hicosmo_de_Omega_mnu0);
  nc_hicosmo_set_Press_mnu0_impl (parent_class, &_nc_hicosmo_de_Press_mnu0);

  nc_hicosmo_set_Omega_mnu0_n_impl (parent_class, &_nc_hicosmo_de_Omega_mnu0_n);
  nc_hicosmo_set_Press_mnu0_n_impl (parent_class, &_nc_hicosmo_de_Press_mnu0_n);
  
  nc_hicosmo_set_E2Omega_mnu_impl (parent_class, &_nc_hicosmo_de_E2Omega_mnu);
  nc_hicosmo_set_E2Press_mnu_impl (parent_class, &_nc_hicosmo_de_E2Press_mnu);

  nc_hicosmo_set_E2Omega_mnu_n_impl (parent_class, &_nc_hicosmo_de_E2Omega_mnu_n);
  nc_hicosmo_set_E2Press_mnu_n_impl (parent_class, &_nc_hicosmo_de_E2Press_mnu_n);

  nc_hicosmo_set_E2Omega_m_impl   (parent_class, &_nc_hicosmo_de_E2Omega_m);
  nc_hicosmo_set_E2Omega_r_impl   (parent_class, &_nc_hicosmo_de_E2Omega_r);

  nc_hicosmo_set_get_bg_var_impl (parent_class, &_nc_hicosmo_de_get_bg_var);
  
  klass->E2Omega_de       = &_nc_hicosmo_de_E2Omega_de;
  klass->dE2Omega_de_dz   = &_nc_hicosmo_de_dE2Omega_de_dz;
  klass->d2E2Omega_de_dz2 = &_nc_hicosmo_de_d2E2Omega_de_dz2;
  klass->w_de             = &_nc_hicosmo_de_w_de;

  model_class->valid = _nc_hicosmo_de_valid;
}

static gdouble _nc_hicosmo_de_Omega_mnu0_n (NcHICosmo *cosmo, const guint n);
static gdouble _nc_hicosmo_de_Omega_gnu0 (NcHICosmo *cosmo);

#define VECTOR (NCM_MODEL (cosmo)->params)
#define MACRO_H0 (ncm_vector_get (VECTOR, NC_HICOSMO_DE_H0))
#define OMEGA_C (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_C))
#define OMEGA_X (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_X))
#define T_GAMMA0 (ncm_vector_get (VECTOR, NC_HICOSMO_DE_T_GAMMA0))
#define HE_YP (ncm_vector_get (VECTOR, NC_HICOSMO_DE_HE_YP))
#define ENNU (ncm_vector_get (VECTOR, NC_HICOSMO_DE_ENNU))
#define OMEGA_R (_nc_hicosmo_de_Omega_gnu0 (NC_HICOSMO (cosmo)))
#define OMEGA_B (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_B))

#define OMEGA_M (OMEGA_B + OMEGA_C)
#define OMEGA_K (1.0 - (OMEGA_B + OMEGA_C + OMEGA_R + OMEGA_X + _nc_hicosmo_de_Omega_mnu0 (cosmo)))

typedef struct _NcHICosmoDENuInt
{
  NcHICosmoDE *cosmo_de;
  gdouble xi2_0;
  gdouble yi;
} NcHICosmoDENuInt;

typedef struct _neutrino_int
{
  const gdouble xi2;
  const gdouble yi;
} neutrino_int;

static gdouble
_nc_hicosmo_de_nu_rho_f (gdouble z, gpointer userdata)
{
  NcHICosmoDENuInt *nu_int = (NcHICosmoDENuInt *) userdata;
  neutrino_int nudata = { nu_int->xi2_0 / gsl_pow_2 (1.0 + z), nu_int->yi };

  ncm_integral1d_ptr_set_userdata (nu_int->cosmo_de->priv->nu_rho, &nudata);

  {
    gdouble err = 0.0;
    const gdouble int_rho = ncm_integral1d_eval_gauss_laguerre (NCM_INTEGRAL1D (nu_int->cosmo_de->priv->nu_rho), &err);

    return int_rho;
  }
}

static gdouble
_nc_hicosmo_de_nu_p_f (gdouble z, gpointer userdata)
{
  NcHICosmoDENuInt *nu_int = (NcHICosmoDENuInt *) userdata;
  neutrino_int nudata = { nu_int->xi2_0 / gsl_pow_2 (1.0 + z), nu_int->yi };

  ncm_integral1d_ptr_set_userdata (nu_int->cosmo_de->priv->nu_p, &nudata);

  {
    gdouble err = 0.0;
    const gdouble int_p = ncm_integral1d_eval_gauss_laguerre (NCM_INTEGRAL1D (nu_int->cosmo_de->priv->nu_p), &err);

    return int_p;
  }
}

static void
_nc_hicosmo_de_prepare (NcHICosmoDE *cosmo_de)
{
  NcmModel *model = NCM_MODEL (cosmo_de);
  if (!ncm_model_state_is_update (model))
  {
    const guint m_len  = ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_M);
    if (m_len > 0)
    {
      gint n;

      for (n = 0; n < m_len; n++)
      {
        const gdouble m       = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_M, n);
        const gdouble Tgamma0 = nc_hicosmo_T_gamma0 (NC_HICOSMO (cosmo_de));
        const gdouble T0      = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_T, n) * Tgamma0;
        const gdouble yi      = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_MU, n);
        const gdouble xi_0    = m * ncm_c_eV () / (ncm_c_kb () * T0);
        const gdouble xi2_0   = xi_0 * xi_0;

        NcHICosmoDENuInt nu_int = {cosmo_de, xi2_0, yi};
        gsl_function F;

        F.params   = &nu_int;

        F.function = &_nc_hicosmo_de_nu_rho_f;
        ncm_spline_set_func (cosmo_de->priv->nu_rho_s[n], NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT, &F, 0.0, cosmo_de->priv->zmax, 0, _NC_HICOSMO_DE_MNU_PREC);

        F.function = &_nc_hicosmo_de_nu_p_f;
        ncm_spline_set_func (cosmo_de->priv->nu_p_s[n], NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT, &F, 0.0, cosmo_de->priv->zmax, 0, _NC_HICOSMO_DE_MNU_PREC);
      }
    }

    ncm_model_state_set_update (model);
  }
  else
    return;
}


/****************************************************************************
 * Normalized Hubble function
 ****************************************************************************/

static gdouble
_nc_hicosmo_de_E2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble Omega_k     = OMEGA_K;
  const gdouble x           = 1.0 + z;
  const gdouble x2          = x * x;
  const gdouble x3          = x2 * x;
  const gdouble x4          = x3 * x;
  const gdouble E2Omega_de  = nc_hicosmo_de_E2Omega_de (NC_HICOSMO_DE (cosmo), z);
  const gdouble E2Omega_mnu = _nc_hicosmo_de_E2Omega_mnu (cosmo, z);
  
  const gdouble E2          = OMEGA_R * x4 + OMEGA_M * x3 + Omega_k * x2 + E2Omega_de + E2Omega_mnu;
  
  return E2;
}

/****************************************************************************
 * dE2_dz
 ****************************************************************************/

static gdouble
_nc_hicosmo_de_dE2_dz (NcHICosmo *cosmo, gdouble z)
{
  const gdouble Omega_k         = OMEGA_K;
  const gdouble x               = 1.0 + z;
  const gdouble x2              = x * x;
  const gdouble x3              = x2 * x;
  const gdouble dE2Omega_mnu_dz = 3.0 * (_nc_hicosmo_de_E2Omega_mnu (cosmo, z) + _nc_hicosmo_de_E2Press_mnu (cosmo, z)) / x;

  return (4.0 * OMEGA_R * x3 + 3.0 * OMEGA_M * x2 + 2.0 * Omega_k * x 
    + nc_hicosmo_de_dE2Omega_de_dz (NC_HICOSMO_DE (cosmo), z))
    + dE2Omega_mnu_dz;
}

/****************************************************************************
 * d2E2_dz2
 ****************************************************************************/

static gdouble
_nc_hicosmo_de_d2E2_dz2 (NcHICosmo *cosmo, gdouble z)
{
  const gdouble omega_k = OMEGA_K;
  const gdouble x       = 1.0 + z;
  const gdouble x2      = x * x;
  
  return (12.0 * OMEGA_R * x2 + 6.0 * OMEGA_M * x + 2.0 * omega_k + nc_hicosmo_de_d2E2Omega_de_dz2 (NC_HICOSMO_DE (cosmo), z));
}

/****************************************************************************
 * Simple functions
 ****************************************************************************/
static gdouble
_nc_hicosmo_de_H0 (NcHICosmo *cosmo)
{
  return MACRO_H0;
}

static gdouble
_nc_hicosmo_de_Omega_t0 (NcHICosmo *cosmo)
{
  return OMEGA_M + OMEGA_X + OMEGA_R + _nc_hicosmo_de_Omega_mnu0 (cosmo);
}

static gdouble
_nc_hicosmo_de_Omega_c0 (NcHICosmo *cosmo)
{
  return OMEGA_C;
}

static gdouble
_nc_hicosmo_de_T_gamma0 (NcHICosmo *cosmo)
{
  return T_GAMMA0;
}

static gdouble
_nc_hicosmo_de_Yp_4He (NcHICosmo *cosmo)
{
  NcmModel *model = NCM_MODEL (cosmo);
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (cosmo);

  if (ncm_model_param_get_ftype (model, NC_HICOSMO_DE_HE_YP) == NCM_PARAM_TYPE_FIXED)
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
        const gdouble Yp = ncm_spline2d_eval (NC_HICOSMO_DE (cosmo)->priv->BBN_spline2d,
                                              wb, DNeff);
        ncm_model_orig_param_set (model, NC_HICOSMO_DE_HE_YP, Yp);
        cosmo_de->priv->HE4_Yp_key = model->pkey;
        /*printf ("# omega_b % 20.15g DeltaNnu % 20.15g Yp % 20.15g\n",  wb, DNeff, Yp);*/
      }
    }
  }

  return HE_YP;
}

static gdouble
_nc_hicosmo_de_Omega_g0 (NcHICosmo *cosmo)
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
_nc_hicosmo_de_Omega_gnu0 (NcHICosmo *cosmo)
{
  const gdouble conv = 7.0 / 8.0 * pow (4.0 / 11.0, 4.0 / 3.0);
  
  return (1.0 + ENNU * conv) * _nc_hicosmo_de_Omega_g0 (cosmo);
}

static gdouble
_nc_hicosmo_de_Omega_r0 (NcHICosmo *cosmo)
{
  const gdouble Omega_mnu_r = 3.0 * _nc_hicosmo_de_Press_mnu0 (cosmo);  
  return _nc_hicosmo_de_Omega_gnu0 (cosmo) + Omega_mnu_r;
}

static gdouble
_nc_hicosmo_de_E2Omega_r (NcHICosmo *cosmo, const gdouble z)
{
  const gdouble x4          = gsl_pow_4 (1.0 + z);
  const gdouble conv        = 7.0 / 8.0 * pow (4.0 / 11.0, 4.0 / 3.0);
  const gdouble Omega_mnu_r = 3.0 * _nc_hicosmo_de_E2Press_mnu (cosmo, z);
  
  return (1.0 + ENNU * conv) * _nc_hicosmo_de_Omega_g0 (cosmo) * x4 + Omega_mnu_r;
}

static gdouble
_nc_hicosmo_de_Omega_m0 (NcHICosmo *cosmo)
{
  const gdouble Omega_mnu_r = 3.0 * _nc_hicosmo_de_Press_mnu0 (cosmo);
  const gdouble Omega_mnu_d = _nc_hicosmo_de_Omega_mnu0 (cosmo) - Omega_mnu_r;
  
  return OMEGA_M + Omega_mnu_d;
}

static gdouble
_nc_hicosmo_de_E2Omega_m (NcHICosmo *cosmo, const gdouble z)
{
  const gdouble x3          = gsl_pow_3 (1.0 + z);
  const gdouble Omega_mnu_r = 3.0 * _nc_hicosmo_de_E2Press_mnu (cosmo, z);
  const gdouble Omega_mnu_d = _nc_hicosmo_de_E2Omega_mnu (cosmo, z) - Omega_mnu_r;
  
  return OMEGA_M * x3 + Omega_mnu_d;
}

static gdouble
_nc_hicosmo_de_Omega_b0 (NcHICosmo *cosmo)
{
  return OMEGA_B;
}

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
_nc_hicosmo_de_MassNuInfo (NcHICosmo *cosmo, guint nu_i, gdouble *mass_eV, gdouble *T_0, gdouble *mu, gdouble *g)
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

  if (ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_MU) == 1)
  {
    mu[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_MU, 0);
  }
  else
  {
    mu[0] = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_MU, nu_i);
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
_nc_hicosmo_de_neutrino_rho_integrand (gpointer userdata, const gdouble u, const gdouble w)
{
  /* see equations 4.11 to 4.16 of Peter & Uzan, Primordial Cosmology with change of variable v = u - x_i and Gauss-Laguerre integration */
  neutrino_int *nudata  = (neutrino_int *) userdata;
  const gdouble u2      = u * u;
  const gdouble exp_yi  = exp (nudata->yi);
  const gdouble exp_myi = 1.0 / exp_yi;

  return u2 * sqrt (u2 + nudata->xi2) * (1.0 / (w + exp_myi) + 1.0 / (w + exp_yi));
}

static gdouble
_nc_hicosmo_de_neutrino_p_integrand (gpointer userdata, const gdouble u, const gdouble w)
{
  /* see equations 4.11 to 4.16 of Peter & Uzan, Primordial Cosmology with change of variable v = u - x_i and Gauss-Laguerre integration */
  neutrino_int *nudata  = (neutrino_int *) userdata;
  const gdouble u2      = u * u;
  const gdouble u4      = u2 * u2;
  const gdouble exp_yi  = exp (nudata->yi);
  const gdouble exp_myi = 1.0 / exp_yi;

  return u4 / (3.0 * sqrt (u2 + nudata->xi2)) * (1.0 / (w + exp_myi) + 1.0 / (w + exp_yi));
}

typedef struct _NcHICosmoDENeutrinoparams
{
  NcHICosmoDE *cosmo_de;
  NcmIntegral1dPtr *int1d_ptr;
} NcHICosmoDENeutrinoparams;

static gdouble
_nc_hicosmo_de_Omega_mnu0_n (NcHICosmo *cosmo, const guint n)
{
  return _nc_hicosmo_de_E2Omega_mnu_n (cosmo, n, 0.0);
}

static gdouble
_nc_hicosmo_de_Press_mnu0_n (NcHICosmo *cosmo, const guint n)
{
  return _nc_hicosmo_de_E2Press_mnu_n (cosmo, n, 0.0);
}

static gdouble
_nc_hicosmo_de_Omega_mnu0 (NcHICosmo *cosmo)
{
  return _nc_hicosmo_de_E2Omega_mnu (cosmo, 0.0);
}

static gdouble
_nc_hicosmo_de_Press_mnu0 (NcHICosmo *cosmo)
{
  return _nc_hicosmo_de_E2Press_mnu (cosmo, 0.0);
}

static gdouble
_nc_hicosmo_de_E2Omega_mnu_n (NcHICosmo *cosmo, const guint n, const gdouble z)
{
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (cosmo);
  NcmModel *model       = NCM_MODEL (cosmo);
  const gdouble h2      = nc_hicosmo_h2 (cosmo);
  const gdouble ffac    = 15.0 / (2.0 * gsl_pow_4 (M_PI)); 
  const gdouble Tgamma  = (1.0 + z) * nc_hicosmo_T_gamma0 (cosmo);
  const gdouble T       = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_T, n) * Tgamma;
  const gdouble g       = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_G, n);

  _nc_hicosmo_de_prepare (NC_HICOSMO_DE (cosmo_de));

  {
    const gdouble int_Ef       = ncm_spline_eval (cosmo_de->priv->nu_rho_s[n], z > cosmo_de->priv->zmax ? cosmo_de->priv->zmax : z);
    const gdouble Omega_mnu0_n = ffac * g * gsl_pow_4 (T) * ncm_c_blackbody_per_crit_density_h2 () / h2 * int_Ef;
    
    return Omega_mnu0_n;
  }
}

static gdouble
_nc_hicosmo_de_E2Press_mnu_n (NcHICosmo *cosmo, const guint n, const gdouble z)
{
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (cosmo);
  NcmModel *model       = NCM_MODEL (cosmo);
  const gdouble h2      = nc_hicosmo_h2 (cosmo);
  const gdouble ffac    = 15.0 / (2.0 * gsl_pow_4 (M_PI)); 
  const gdouble Tgamma  = (1.0 + z) * nc_hicosmo_T_gamma0 (cosmo);
  const gdouble T       = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_T, n) * Tgamma;
  const gdouble g       = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_G, n);

  _nc_hicosmo_de_prepare (NC_HICOSMO_DE (cosmo_de));

  {
    const gdouble int_pf       = ncm_spline_eval (cosmo_de->priv->nu_p_s[n], z > cosmo_de->priv->zmax ? cosmo_de->priv->zmax : z);
    const gdouble Press_mnu0_n = ffac * g * gsl_pow_4 (T) * ncm_c_blackbody_per_crit_density_h2 () / h2 * int_pf;

    return Press_mnu0_n;
  }
}

static gdouble
_nc_hicosmo_de_E2Omega_mnu (NcHICosmo *cosmo, const gdouble z)
{
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (cosmo);
  NcmModel *model       = NCM_MODEL (cosmo);
  const guint nmassnu   = ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_M);
  const gdouble h2      = nc_hicosmo_h2 (cosmo);
  const gdouble ffac    = 15.0 / (2.0 * gsl_pow_4 (M_PI)); 
  const gdouble Tgamma  = (1.0 + z) * nc_hicosmo_T_gamma0 (cosmo);
  gdouble Omega_mnu0    = 0.0;

  guint n;

  _nc_hicosmo_de_prepare (NC_HICOSMO_DE (cosmo_de));
  
  for (n = 0; n < nmassnu; n++)
  {
    const gdouble T            = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_T, n) * Tgamma;
    const gdouble g            = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_G, n);
    const gdouble int_Ef       = ncm_spline_eval (cosmo_de->priv->nu_rho_s[n], z > cosmo_de->priv->zmax ? cosmo_de->priv->zmax : z);
    const gdouble Omega_mnu0_n = ffac * g * gsl_pow_4 (T) * ncm_c_blackbody_per_crit_density_h2 () / h2 * int_Ef;

    Omega_mnu0 += Omega_mnu0_n;
  }

  return Omega_mnu0;
}

static gdouble
_nc_hicosmo_de_E2Press_mnu (NcHICosmo *cosmo, const gdouble z)
{
  NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (cosmo);
  NcmModel *model       = NCM_MODEL (cosmo);
  const guint nmassnu   = ncm_model_vparam_len (model, NC_HICOSMO_DE_MASSNU_M);
  const gdouble Tgamma  = (1.0 + z) * nc_hicosmo_T_gamma0 (cosmo);
  const gdouble h2      = nc_hicosmo_h2 (cosmo);
  const gdouble ffac    = 15.0 / (2.0 * gsl_pow_4 (M_PI)); 
  gdouble Press_mnu0    = 0.0;
  guint n;

  _nc_hicosmo_de_prepare (NC_HICOSMO_DE (cosmo_de));
  
  for (n = 0; n < nmassnu; n++)
  {
    const gdouble T            = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_T, n) * Tgamma;
    const gdouble g            = ncm_model_orig_vparam_get (model, NC_HICOSMO_DE_MASSNU_G, n);
    const gdouble int_pf       = ncm_spline_eval (cosmo_de->priv->nu_p_s[n], z > cosmo_de->priv->zmax ? cosmo_de->priv->zmax : z);
    const gdouble Press_mnu0_n = ffac * g * gsl_pow_4 (T) * ncm_c_blackbody_per_crit_density_h2 () / h2 * int_pf;

    Press_mnu0 += Press_mnu0_n;
  }
  return Press_mnu0;
}

static void 
_nc_hicosmo_de_get_bg_var (NcHICosmo *cosmo, const gdouble t, NcHIPertBGVar *bg_var)
{
}

static gdouble
min_E2 (gdouble z, gpointer params)
{
  NcHICosmo *cosmo = NC_HICOSMO (params);
  return _nc_hicosmo_de_E2 (cosmo, z);
}

static gboolean 
_nc_hicosmo_de_valid (NcmModel *model)
{
  if (nc_hicosmo_Omega_k0 (NC_HICOSMO (model)) < 0.0)
  {
    NcHICosmoDE *cosmo_de = NC_HICOSMO_DE (model);
    gdouble E2_min = 1.0e300;
    gdouble z_min  = 0.0;
    gdouble z_max  = 4.0;
    gint i;

    for (i = 0 ; i <= 40; i++)
    {
      const gdouble z  = 0.05 * i;
      const gdouble E2 = _nc_hicosmo_de_E2 (NC_HICOSMO (model), z);

      if (E2 < E2_min)
      {
        E2_min = E2;
        z_min  = z;
      }
    }

    if (E2_min <= 0.0)
      return FALSE;

    z_max = 1.0e4;
    
    {
      gsl_function F;
      gint max_iter = 10000;
      gint iter = 0;
      gint status;
      gdouble a, b;

      F.function = &min_E2;
      F.params   = model;

      gsl_min_fminimizer_set (cosmo_de->priv->min, &F, z_min, 0.0, z_max);

      do
      {
        iter++;
        status = gsl_min_fminimizer_iterate (cosmo_de->priv->min);
        z_min  = gsl_min_fminimizer_x_minimum (cosmo_de->priv->min);
        a      = gsl_min_fminimizer_x_lower (cosmo_de->priv->min);
        b      = gsl_min_fminimizer_x_upper (cosmo_de->priv->min);

        status = gsl_min_test_interval (a, b, 0.01, 0.0);        
      }
      while (status == GSL_CONTINUE && iter < max_iter);
/*
      if (_nc_hicosmo_de_E2 (NC_HICOSMO (model), z_min) <= 0.0)
        printf ("FOUND % 22.15g % 22.15g % 22.15g\n", 
                z_min, 
                nc_hicosmo_Omega_k0 (NC_HICOSMO (model)), 
                _nc_hicosmo_de_E2 (NC_HICOSMO (model), z_min));
*/      
      return _nc_hicosmo_de_E2 (NC_HICOSMO (model), z_min) > 0.0;
    }
  }
  else
    return TRUE;
}

void
nc_hicosmo_de_set_wmap5_params (NcHICosmoDE *cosmo_de)
{
  g_assert (NC_IS_HICOSMO_DE (cosmo_de));
  ncm_model_orig_param_set (NCM_MODEL (cosmo_de), NC_HICOSMO_DE_H0, 72.4000);
  ncm_model_orig_param_set (NCM_MODEL (cosmo_de), NC_HICOSMO_DE_OMEGA_C, 0.2060);
  ncm_model_orig_param_set (NCM_MODEL (cosmo_de), NC_HICOSMO_DE_OMEGA_X, 0.7510);
  ncm_model_orig_param_set (NCM_MODEL (cosmo_de), NC_HICOSMO_DE_OMEGA_B, 0.0432);
  ncm_model_orig_param_set (NCM_MODEL (cosmo_de), NC_HICOSMO_DE_T_GAMMA0, 2.7250);
  ncm_model_orig_param_set (NCM_MODEL (cosmo_de), NC_HICOSMO_DE_ENNU, 3.046);
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

  ncm_reparam_free (NCM_REPARAM (de_reparam_ok));
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

  ncm_reparam_free (NCM_REPARAM (de_reparam_cmb));
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

static gdouble
_nc_hicosmo_de_E2Omega_de (NcHICosmoDE *cosmo_de, gdouble z)
{
  g_error ("nc_hicosmo_de_E2Omega_de: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo_de));
  return 0.0;
}

static gdouble
_nc_hicosmo_de_dE2Omega_de_dz (NcHICosmoDE *cosmo_de, gdouble z)
{
  g_error ("nc_hicosmo_de_dE2Omega_de_dz: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo_de));
  return 0.0;
}

static gdouble
_nc_hicosmo_de_d2E2Omega_de_dz2 (NcHICosmoDE *cosmo_de, gdouble z)
{
  g_error ("nc_hicosmo_de_d2E2Omega_de_dz2: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo_de));
  return 0.0;
}

static gdouble
_nc_hicosmo_de_w_de (NcHICosmoDE *cosmo_de, gdouble z)
{
  g_error ("nc_hicosmo_de_w_de: model `%s' does not implement this function.", G_OBJECT_TYPE_NAME (cosmo_de));
  return 0.0;
}

#define NC_HICOSMO_DE_SET_IMPL_FUNC(name)                                                                    \
	void                                                                                                       \
	nc_hicosmo_de_set_##name##_impl (NcHICosmoDEClass *cosmo_de_class, NcmFuncF f, NcmFuncPF pf, NcmFuncDF df) \
	{                                                                                                          \
		ncm_model_class_add_impl_opts (NCM_MODEL_CLASS (cosmo_de_class), NC_HICOSMO_DE_IMPL_##name, -1);         \
		g_assert (f != NULL);                                                                                    \
		cosmo_de_class->name = *ncm_func_stub;                                                                   \
		cosmo_de_class->name.f = f;                                                                              \
		cosmo_de_class->name.pf = pf;                                                                            \
		cosmo_de_class->name.df = df;                                                                            \
	}

/**
 * nc_hicosmo_de_set_E2Omega_de_impl: (skip)
 * @cosmo_de_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC (NC_HICOSMO_DE, NcHICosmoDE, nc_hicosmo_de, NcHICosmoDEFunc1, E2Omega_de)
/**
 * nc_hicosmo_de_set_dE2Omega_de_dz_impl: (skip)
 * @cosmo_de_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC (NC_HICOSMO_DE, NcHICosmoDE, nc_hicosmo_de, NcHICosmoDEFunc1, dE2Omega_de_dz)
/**
 * nc_hicosmo_de_set_d2E2Omega_de_dz2_impl: (skip)
 * @cosmo_de_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC (NC_HICOSMO_DE, NcHICosmoDE, nc_hicosmo_de, NcHICosmoDEFunc1, d2E2Omega_de_dz2)
/**
 * nc_hicosmo_de_set_w_de_impl: (skip)
 * @cosmo_de_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC (NC_HICOSMO_DE, NcHICosmoDE, nc_hicosmo_de, NcHICosmoDEFunc1, w_de)
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
  ncm_mset_func_list_register ("wDE", "\\omega_\\mathrm{de}", "NcHICosmoDE", "Darkenergy equation of state today", G_TYPE_NONE, _nc_hicosmo_de_flist_w0, 0, 1);
  ncm_mset_func_list_register ("BBN", "BBN", "NcHICosmoDE", "BBN", G_TYPE_NONE, _nc_hicosmo_de_flist_BBN, 0, 1);
}
