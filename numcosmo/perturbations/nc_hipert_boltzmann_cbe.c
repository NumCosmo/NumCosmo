/***************************************************************************
 *            nc_hipert_boltzmann_cbe.c
 *
 *  Sat October 24 11:56:56 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_boltzmann_cbe.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_hipert_boltzmann_cbe
 * @title: NcHIPertBoltzmannCBE
 * @short_description: CLASS (Cosmic Linear Anisotropy Solving System) backend for perturbations
 *
 * FIXME
 *
 * If you use this object please cite: [Blas (2011) CLASS II][XBlas2011],
 * see also:
 * - [Lesgourgues (2011) CLASS I][XLesgourgues2011],
 * - [Lesgourgues (2011) CLASS III][XLesgourgues2011a],
 * - [Lesgourgues (2011) CLASS IV][XLesgourgues2011b] and
 * - [CLASS website](http://class-code.net/).
 *
 *
 */

/*
 * It must be include before anything else, several symbols clash
 * with the default includes.
 */
#include "class/include/class.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hipert_boltzmann_cbe.h"

enum
{
  PROP_0,
  PROP_A_INI_A_0,
  PROP_BACK_INT_STEP,
  PROP_BACK_TOL,
  PROP_INITIAL_OMEGA_R_TOL,
  PROP_M_NCDM_TOL,
  PROP_NCDM_TOL,
  PROP_NCDM_SYNCHRONOUS_TOL,
  PROP_NCDM_NEWTONIAN_TOL,
  PROP_NCDM_BG_TOL,
  PROP_NCDM_INITIAL_W_TOL,
  PROP_SBBN_FILE,
  PROP_RECFAST_Z_INI,
  PROP_RECFAST_NZ0,
  PROP_THERMO_INTEGRATION_TOL,
  PROP_RECFAST_HE_SWITCH,
  PROP_RECFAST_FUDGE_HE,
  PROP_RECFAST_H_SWITCH,
  PROP_RECFAST_FUDGE_H,
  PROP_RECFAST_DELTA_FUDGE_H,
  PROP_RECFAST_A_GAUSS_1,
  PROP_RECFAST_A_GAUSS_2,
  PROP_RECFAST_Z_GAUSS_1,
  PROP_RECFAST_Z_GAUSS_2,
  PROP_RECFAST_W_GAUSS_1,
  PROP_RECFAST_W_GAUSS_2,
 };

struct _NcHIPertBoltzmannCBEPrivate
{
  struct precision ppr;
  struct background pba;
  struct thermo pth;
  struct perturbs ppt;
  struct transfers ptr;
  struct primordial ppm;
  struct spectra psp;
  struct nonlinear pnl;
  struct lensing ple;
  struct output pop;
};

G_DEFINE_TYPE (NcHIPertBoltzmannCBE, nc_hi_pert_boltzmann_cbe, NC_TYPE_HIPERT_BOLTZMANN);

static void
nc_hi_pert_boltzmann_cbe_init (NcHIPertBoltzmannCBE *cbe)
{
  cbe->priv = G_TYPE_INSTANCE_GET_PRIVATE (cbe, NC_TYPE_HIPERT_BOLTZMANN_CBE, NcHIPertBoltzmannCBEPrivate);
}

static void
nc_hi_pert_boltzmann_cbe_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPertBoltzmannCBE *pclass = NC_HIPERT_BOLTZMANN_CBE (object);
  g_return_if_fail (NC_IS_HIPERT_BOLTZMANN_CBE (object));

  switch (prop_id)
  {
    case PROP_A_INI_A_0:
      pclass->priv->ppr.a_ini_over_a_today_default = g_value_get_double (value);
      break;
    case PROP_BACK_INT_STEP:
      pclass->priv->ppr.back_integration_stepsize  = g_value_get_double (value);
      break;
    case PROP_BACK_TOL:
      pclass->priv->ppr.tol_background_integration = g_value_get_double (value);
      break;
    case PROP_INITIAL_OMEGA_R_TOL:
      pclass->priv->ppr.tol_initial_Omega_r        = g_value_get_double (value);
      break;
    case PROP_M_NCDM_TOL:
      pclass->priv->ppr.tol_M_ncdm                 = g_value_get_double (value);
      break;
    case PROP_NCDM_TOL:
      pclass->priv->ppr.tol_ncdm                   = g_value_get_double (value);
      break;
    case PROP_NCDM_SYNCHRONOUS_TOL:
      pclass->priv->ppr.tol_ncdm_synchronous       = g_value_get_double (value);
      break;
    case PROP_NCDM_NEWTONIAN_TOL:
      pclass->priv->ppr.tol_ncdm_newtonian         = g_value_get_double (value);
      break;
    case PROP_NCDM_BG_TOL:
      pclass->priv->ppr.tol_ncdm_bg                = g_value_get_double (value);
      break;
    case PROP_NCDM_INITIAL_W_TOL:
      pclass->priv->ppr.tol_ncdm_initial_w         = g_value_get_double (value);
      break;
    case PROP_SBBN_FILE:
    {
      const gchar *sBBN_file = g_value_get_string (value);
      const guint sBBN_file_len = strlen (sBBN_file);
      g_assert_cmpuint (sBBN_file_len, <= ,_FILENAMESIZE_);
      memcpy (pclass->priv->ppr.sBBN_file, sBBN_file, sBBN_file_len);
      break;
    }
    case PROP_RECFAST_Z_INI:
      pclass->priv->ppr.recfast_z_initial          = g_value_get_double (value);
      break;
    case PROP_RECFAST_NZ0:
      pclass->priv->ppr.recfast_Nz0                = g_value_get_double (value);
      break;
    case PROP_THERMO_INTEGRATION_TOL:
      pclass->priv->ppr.tol_thermo_integration     = g_value_get_double (value);
      break;
    case PROP_RECFAST_HE_SWITCH:
      pclass->priv->ppr.recfast_Heswitch           = g_value_get_int (value);
      break;
    case PROP_RECFAST_FUDGE_HE:
      pclass->priv->ppr.recfast_fudge_He           = g_value_get_double (value);
      break;
    case PROP_RECFAST_H_SWITCH:
      pclass->priv->ppr.recfast_Hswitch            = g_value_get_int (value);
      break;
    case PROP_RECFAST_FUDGE_H:
      pclass->priv->ppr.recfast_fudge_H            = g_value_get_double (value);
      break;
    case PROP_RECFAST_DELTA_FUDGE_H:
      pclass->priv->ppr.recfast_delta_fudge_H      = g_value_get_double (value);
      break;
    case PROP_RECFAST_A_GAUSS_1:
      pclass->priv->ppr.recfast_AGauss1            = g_value_get_double (value);
      break;
    case PROP_RECFAST_A_GAUSS_2:
      pclass->priv->ppr.recfast_AGauss2            = g_value_get_double (value);
      break;
    case PROP_RECFAST_Z_GAUSS_1:
      pclass->priv->ppr.recfast_zGauss1            = g_value_get_double (value);
      break;
    case PROP_RECFAST_Z_GAUSS_2:
      pclass->priv->ppr.recfast_zGauss2            = g_value_get_double (value);
      break;
    case PROP_RECFAST_W_GAUSS_1:
      pclass->priv->ppr.recfast_wGauss1            = g_value_get_double (value);
      break;
    case PROP_RECFAST_W_GAUSS_2:
      pclass->priv->ppr.recfast_wGauss2            = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hi_pert_boltzmann_cbe_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIPertBoltzmannCBE *pclass = NC_HIPERT_BOLTZMANN_CBE (object);
  g_return_if_fail (NC_IS_HIPERT_BOLTZMANN_CBE (object));

  switch (prop_id)
  {
    case PROP_A_INI_A_0:
      g_value_set_double (value, pclass->priv->ppr.a_ini_over_a_today_default);
      break;
    case PROP_BACK_INT_STEP:
      g_value_set_double (value, pclass->priv->ppr.back_integration_stepsize);
      break;
    case PROP_BACK_TOL:
      g_value_set_double (value, pclass->priv->ppr.tol_background_integration);
      break;
    case PROP_INITIAL_OMEGA_R_TOL:
      g_value_set_double (value, pclass->priv->ppr.tol_initial_Omega_r);
      break;
    case PROP_M_NCDM_TOL:
      g_value_set_double (value, pclass->priv->ppr.tol_M_ncdm);
      break;
    case PROP_NCDM_TOL:
      g_value_set_double (value, pclass->priv->ppr.tol_ncdm);
      break;
    case PROP_NCDM_SYNCHRONOUS_TOL:
      g_value_set_double (value, pclass->priv->ppr.tol_ncdm_synchronous);
      break;
    case PROP_NCDM_NEWTONIAN_TOL:
      g_value_set_double (value, pclass->priv->ppr.tol_ncdm_newtonian);
      break;
    case PROP_NCDM_BG_TOL:
      g_value_set_double (value, pclass->priv->ppr.tol_ncdm_bg);
      break;
    case PROP_NCDM_INITIAL_W_TOL:
      g_value_set_double (value, pclass->priv->ppr.tol_ncdm_initial_w);
      break;
    case PROP_SBBN_FILE:
      g_value_set_string (value, pclass->priv->ppr.sBBN_file);
      break;
    case PROP_RECFAST_Z_INI:
      g_value_set_double (value, pclass->priv->ppr.recfast_z_initial);
      break;
    case PROP_RECFAST_NZ0:
      g_value_set_int (value, pclass->priv->ppr.recfast_Nz0);
      break;
    case PROP_THERMO_INTEGRATION_TOL:
      g_value_set_double (value, pclass->priv->ppr.tol_thermo_integration);
      break;
    case PROP_RECFAST_HE_SWITCH:
      g_value_set_int (value, pclass->priv->ppr.recfast_Heswitch);
      break;
    case PROP_RECFAST_FUDGE_HE:
      g_value_set_double (value, pclass->priv->ppr.recfast_fudge_He);
      break;
    case PROP_RECFAST_H_SWITCH:
      g_value_set_int (value, pclass->priv->ppr.recfast_Hswitch);
      break;
    case PROP_RECFAST_FUDGE_H:
      g_value_set_double (value, pclass->priv->ppr.recfast_fudge_H);
      break;
    case PROP_RECFAST_DELTA_FUDGE_H:
      g_value_set_double (value, pclass->priv->ppr.recfast_delta_fudge_H);
      break;
    case PROP_RECFAST_A_GAUSS_1:
      g_value_set_double (value, pclass->priv->ppr.recfast_AGauss1);
      break;
    case PROP_RECFAST_A_GAUSS_2:
      g_value_set_double (value, pclass->priv->ppr.recfast_AGauss2);
      break;
    case PROP_RECFAST_Z_GAUSS_1:
      g_value_set_double (value, pclass->priv->ppr.recfast_zGauss1);
      break;
    case PROP_RECFAST_Z_GAUSS_2:
      g_value_set_double (value, pclass->priv->ppr.recfast_zGauss2);
      break;
    case PROP_RECFAST_W_GAUSS_1:
      g_value_set_double (value, pclass->priv->ppr.recfast_wGauss1);
      break;
    case PROP_RECFAST_W_GAUSS_2:
      g_value_set_double (value, pclass->priv->ppr.recfast_wGauss2);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hi_pert_boltzmann_cbe_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hi_pert_boltzmann_cbe_parent_class)->finalize (object);
}

static void
nc_hi_pert_boltzmann_cbe_class_init (NcHIPertBoltzmannCBEClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  g_type_class_add_private (klass, sizeof (NcHIPertBoltzmannCBEPrivate));

  object_class->set_property = nc_hi_pert_boltzmann_cbe_set_property;
  object_class->get_property = nc_hi_pert_boltzmann_cbe_get_property;
  object_class->finalize     = nc_hi_pert_boltzmann_cbe_finalize;

  /*
   * Background related parameters.
   */
  g_object_class_install_property (object_class,
                                   PROP_A_INI_A_0,
                                   g_param_spec_double ("a-ini-over-a-today-default",
                                                        NULL,
                                                        "a_ini / a_0",
                                                        0.0, G_MAXDOUBLE, 1.0e-14,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BACK_INT_STEP,
                                   g_param_spec_double ("back-integration-stepsize",
                                                        NULL,
                                                        "Background integration step size",
                                                        0.0, G_MAXDOUBLE, 7.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BACK_TOL,
                                   g_param_spec_double ("tol-background-integration",
                                                        NULL,
                                                        "Tolerance in background integration",
                                                        0.0, G_MAXDOUBLE, 1.0e-2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_INITIAL_OMEGA_R_TOL,
                                   g_param_spec_double ("tol-initial-Omega-r",
                                                        NULL,
                                                        "Tolerance in initial Omega_r",
                                                        0.0, G_MAXDOUBLE, 1.0e-4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_M_NCDM_TOL,
                                   g_param_spec_double ("tol-M-ncdm",
                                                        NULL,
                                                        "Tolerance M_ncdm",
                                                        0.0, G_MAXDOUBLE, 1.0e-7,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NCDM_TOL,
                                   g_param_spec_double ("tol-ncdm",
                                                        NULL,
                                                        "Tolerance ncdm",
                                                        0.0, G_MAXDOUBLE, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NCDM_SYNCHRONOUS_TOL,
                                   g_param_spec_double ("tol-ncdm-synchronous",
                                                        NULL,
                                                        "Tolerance ncdm synchronous",
                                                        0.0, G_MAXDOUBLE, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NCDM_NEWTONIAN_TOL,
                                   g_param_spec_double ("tol-ncdm-newtonian",
                                                        NULL,
                                                        "Tolerance ncdm newtonian",
                                                        0.0, G_MAXDOUBLE, 1.0e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NCDM_BG_TOL,
                                   g_param_spec_double ("tol-ncdm-bg",
                                                        NULL,
                                                        "Tolerance ncdm bg",
                                                        0.0, G_MAXDOUBLE, 1.0e-5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NCDM_INITIAL_W_TOL,
                                   g_param_spec_double ("tol-ncdm-initial-w",
                                                        NULL,
                                                        "Tolerance ncdm initial w",
                                                        0.0, G_MAXDOUBLE, 1.0e-3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /*
   * Thermodynamics related parameters
   */
  {
    gchar *sBBN_file = ncm_cfg_get_data_filename ("class_data/bbn/sBBN.dat", TRUE);
    g_object_class_install_property (object_class,
                                     PROP_SBBN_FILE,
                                     g_param_spec_string ("sBBN-file",
                                                          NULL,
                                                          "SBBN filename",
                                                          sBBN_file,
                                                          G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  }
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_Z_INI,
                                   g_param_spec_double ("recfast-z-initial",
                                                        NULL,
                                                        "Recfast initial redshift",
                                                        0.0, G_MAXDOUBLE, 1.0e4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_NZ0,
                                   g_param_spec_int ("recfast-Nz0",
                                                     NULL,
                                                     "Recfast Nz0",
                                                     0, G_MAXINT, 20000,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_THERMO_INTEGRATION_TOL,
                                   g_param_spec_double ("tol-thermo-integration",
                                                        NULL,
                                                        "Thermodynamic integration tolerance",
                                                        0.0, G_MAXDOUBLE, 1.0e-2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_HE_SWITCH,
                                   g_param_spec_int ("recfast-Heswitch",
                                                     NULL,
                                                     "Recfast He switch",
                                                     0, G_MAXINT, 6,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_FUDGE_HE,
                                   g_param_spec_double ("recfast-fudge-He",
                                                        NULL,
                                                        "Recfast fudge He",
                                                        0.0, G_MAXDOUBLE, 0.86,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_H_SWITCH,
                                   g_param_spec_int ("recfast-Hswitch",
                                                     NULL,
                                                     "Recfast He switch",
                                                     0, G_MAXINT, _TRUE_,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_FUDGE_H,
                                   g_param_spec_double ("recfast-fudge-H",
                                                        NULL,
                                                        "Recfast fudge H",
                                                        0.0, G_MAXDOUBLE, 1.14,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_DELTA_FUDGE_H,
                                   g_param_spec_double ("recfast-delta-fudge-H",
                                                        NULL,
                                                        "Recfast delta fudge H",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, -0.015,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_A_GAUSS_1,
                                   g_param_spec_double ("recfast-AGauss1",
                                                        NULL,
                                                        "Recfast A Gauss 1",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, -0.14,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_A_GAUSS_2,
                                   g_param_spec_double ("recfast-AGauss2",
                                                        NULL,
                                                        "Recfast A Gauss 2",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.079,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_Z_GAUSS_1,
                                   g_param_spec_double ("recfast-zGauss1",
                                                        NULL,
                                                        "Recfast z Gauss 1",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 7.28,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_Z_GAUSS_2,
                                   g_param_spec_double ("recfast-zGauss2",
                                                        NULL,
                                                        "Recfast z Gauss 2",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 6.73,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_W_GAUSS_1,
                                   g_param_spec_double ("recfast-wGauss1",
                                                        NULL,
                                                        "Recfast w Gauss 1",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.18,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECFAST_W_GAUSS_2,
                                   g_param_spec_double ("recfast-wGauss2",
                                                        NULL,
                                                        "Recfast w Gauss 2",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.33,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

}
