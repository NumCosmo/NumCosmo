/***************************************************************************
 *            nc_hipert_boltzmann.c
 *
 *  Sat Oct 25 21:02:36 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_boltzmann.c
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
 * SECTION:nc_hipert_boltzmann
 * @title: NcHIPertBoltzmann
 * @short_description: Abstract class for perturbative Boltzmann hierarchy.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_enum_types.h"
#include "perturbations/nc_hipert_boltzmann.h"

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <nvector/nvector_serial.h>
#include <gsl/gsl_roots.h>

G_DEFINE_ABSTRACT_TYPE (NcHIPertBoltzmann, nc_hipert_boltzmann, NC_TYPE_HIPERT);

enum
{
  PROP_0,
  PROP_RECOMB,
  PROP_TARGET_CLS,
  PROP_CALC_TRANSFER,
  PROP_USE_LENSED_CLS,
  PROP_USE_TENSOR,
  PROP_TT_LMAX,
  PROP_EE_LMAX,
  PROP_BB_LMAX,
  PROP_TE_LMAX,
  PROP_TB_LMAX,
  PROP_EB_LMAX,
  PROP_SIZE,
};

static void
nc_hipert_boltzmann_init (NcHIPertBoltzmann *pb)
{
  pb->recomb                 = NULL;
  pb->cosmo                  = NULL;
  pb->a                      = NULL;
  pb->eta0                   = 0.0;
  pb->lambdai                = 0.0;
  pb->lambdaf                = 0.0;
  pb->lambda_opt_cutoff      = 0.0;
  pb->lambda_rec             = 0.0;
  pb->lambda_rec_10m2_max[0] = 0.0;
  pb->lambda                 = 0.0;

  pb->target_Cls             = 0;
  pb->use_lensed_Cls         = FALSE;
  pb->use_tensor             = FALSE;
  pb->calc_transfer          = FALSE;
  pb->TT_lmax                = 0;
  pb->EE_lmax                = 0;
  pb->BB_lmax                = 0;
  pb->TE_lmax                = 0;
  pb->TB_lmax                = 0;
  pb->EB_lmax                = 0;
  pb->tight_coupling         = FALSE;

  pb->ctrl_cosmo             = ncm_model_ctrl_new (NULL);
  pb->ctrl_prim              = ncm_model_ctrl_new (NULL);
}

static void
_nc_hipert_boltzmann_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPertBoltzmann *pb = NC_HIPERT_BOLTZMANN (object);
  g_return_if_fail (NC_IS_HIPERT_BOLTZMANN (object));

  switch (prop_id)
  {
    case PROP_RECOMB:
      nc_hipert_boltzmann_set_recomb (pb, g_value_get_object (value));
      break;
    case PROP_TARGET_CLS:
      nc_hipert_boltzmann_set_target_Cls (pb, g_value_get_flags (value));
      break;
    case PROP_CALC_TRANSFER:
      nc_hipert_boltzmann_set_calc_transfer (pb, g_value_get_boolean (value));
      break;
    case PROP_USE_LENSED_CLS:
      nc_hipert_boltzmann_set_lensed_Cls (pb, g_value_get_boolean (value));
      break;
    case PROP_USE_TENSOR:
      nc_hipert_boltzmann_set_tensor (pb, g_value_get_boolean (value));
      break;
    case PROP_TT_LMAX:
      nc_hipert_boltzmann_set_TT_lmax (pb, g_value_get_uint (value));
      break;
    case PROP_EE_LMAX:
      nc_hipert_boltzmann_set_EE_lmax (pb, g_value_get_uint (value));
      break;
    case PROP_BB_LMAX:
      nc_hipert_boltzmann_set_BB_lmax (pb, g_value_get_uint (value));
      break;
    case PROP_TE_LMAX:
      nc_hipert_boltzmann_set_TE_lmax (pb, g_value_get_uint (value));
      break;
    case PROP_TB_LMAX:
      nc_hipert_boltzmann_set_TB_lmax (pb, g_value_get_uint (value));
      break;
    case PROP_EB_LMAX:
      nc_hipert_boltzmann_set_EB_lmax (pb, g_value_get_uint (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_boltzmann_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIPertBoltzmann *pb = NC_HIPERT_BOLTZMANN (object);
  g_return_if_fail (NC_IS_HIPERT_BOLTZMANN (object));

  switch (prop_id)
  {
    case PROP_RECOMB:
      g_value_set_object (value, pb->recomb);
      break;
    case PROP_TARGET_CLS:
      g_value_set_flags (value, nc_hipert_boltzmann_get_target_Cls (pb));
      break;
    case PROP_CALC_TRANSFER:
      g_value_set_boolean (value, nc_hipert_boltzmann_get_calc_transfer (pb));
      break;
    case PROP_USE_LENSED_CLS:
      g_value_set_boolean (value, nc_hipert_boltzmann_lensed_Cls (pb));
      break;
    case PROP_USE_TENSOR:
      g_value_set_boolean (value, nc_hipert_boltzmann_tensor (pb));
      break;
    case PROP_TT_LMAX:
      g_value_set_uint (value, nc_hipert_boltzmann_get_TT_lmax (pb));
      break;
    case PROP_EE_LMAX:
      g_value_set_uint (value, nc_hipert_boltzmann_get_EE_lmax (pb));
      break;
    case PROP_BB_LMAX:
      g_value_set_uint (value, nc_hipert_boltzmann_get_BB_lmax (pb));
      break;
    case PROP_TE_LMAX:
      g_value_set_uint (value, nc_hipert_boltzmann_get_TE_lmax (pb));
      break;
    case PROP_TB_LMAX:
      g_value_set_uint (value, nc_hipert_boltzmann_get_TB_lmax (pb));
      break;
    case PROP_EB_LMAX:
      g_value_set_uint (value, nc_hipert_boltzmann_get_EB_lmax (pb));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_boltzmann_dispose (GObject *object)
{
  NcHIPertBoltzmann *pb = NC_HIPERT_BOLTZMANN (object);

  nc_recomb_clear (&pb->recomb);
  nc_hicosmo_clear (&pb->cosmo);
  nc_scalefactor_clear (&pb->a);

  ncm_model_ctrl_clear (&pb->ctrl_cosmo);
  ncm_model_ctrl_clear (&pb->ctrl_prim);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_boltzmann_parent_class)->dispose (object);
}

static void
_nc_hipert_boltzmann_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_boltzmann_parent_class)->finalize (object);
}

static void _nc_hipert_boltzmann_set_mode_k (NcHIPert *pert, gdouble k);
static void _nc_hipert_boltzmann_set_abstol (NcHIPert *pert, gdouble abstol);
static void _nc_hipert_boltzmann_set_reltol (NcHIPert *pert, gdouble reltol);
static void _nc_hipert_boltzmann_prepare (NcHIPertBoltzmann *pb, NcHIPrim *prim, NcHIReion *reion, NcHICosmo *cosmo);
static void _nc_hipert_boltzmann_prepare_if_needed (NcHIPertBoltzmann *pb, NcHIPrim *prim, NcHIReion *reion, NcHICosmo *cosmo);

static void
nc_hipert_boltzmann_class_init (NcHIPertBoltzmannClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_hipert_boltzmann_set_property;
  object_class->get_property = &_nc_hipert_boltzmann_get_property;
  object_class->dispose      = &_nc_hipert_boltzmann_dispose;
  object_class->finalize     = &_nc_hipert_boltzmann_finalize;

  g_object_class_install_property (object_class,
                                   PROP_RECOMB,
                                   g_param_spec_object ("recomb",
                                                        NULL,
                                                        "Recombination object",
                                                        NC_TYPE_RECOMB,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TARGET_CLS,
                                   g_param_spec_flags ("target-Cls",
                                                       NULL,
                                                       "Which Cls must be calculated",
                                                       NC_TYPE_DATA_CMB_DATA_TYPE, NC_DATA_CMB_TYPE_TT,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_CALC_TRANSFER,
                                   g_param_spec_boolean ("calc-transfer",
                                                         NULL,
                                                         "Whether to calculate the matter transfer function",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_USE_LENSED_CLS,
                                   g_param_spec_boolean ("use-lensed-Cls",
                                                         NULL,
                                                         "Whether use the lensed corrected Cls",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_USE_TENSOR,
                                   g_param_spec_boolean ("use-tensor",
                                                         NULL,
                                                         "Whether use tensor contribution",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TT_LMAX,
                                   g_param_spec_uint ("TT-l-max",
                                                      NULL,
                                                      "Last multipole in the TT correlation",
                                                      4, G_MAXUINT32, 30,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_EE_LMAX,
                                   g_param_spec_uint ("EE-l-max",
                                                      NULL,
                                                      "Last multipole in the EE correlation",
                                                      4, G_MAXUINT32, 4,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_BB_LMAX,
                                   g_param_spec_uint ("BB-l-max",
                                                      NULL,
                                                      "Last multipole in the BB correlation",
                                                      4, G_MAXUINT32, 4,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TE_LMAX,
                                   g_param_spec_uint ("TE-l-max",
                                                      NULL,
                                                      "Last multipole in the TE correlation",
                                                      4, G_MAXUINT32, 4,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_TB_LMAX,
                                   g_param_spec_uint ("TB-l-max",
                                                      NULL,
                                                      "Last multipole in the TB correlation",
                                                      4, G_MAXUINT32, 4,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_EB_LMAX,
                                   g_param_spec_uint ("EB-l-max",
                                                      NULL,
                                                      "Last multipole in the EB correlation",
                                                      4, G_MAXUINT32, 4,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  NC_HIPERT_CLASS (klass)->set_mode_k        = &_nc_hipert_boltzmann_set_mode_k;
  NC_HIPERT_CLASS (klass)->set_abstol        = &_nc_hipert_boltzmann_set_abstol;
  NC_HIPERT_CLASS (klass)->set_reltol        = &_nc_hipert_boltzmann_set_reltol;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->prepare           = &_nc_hipert_boltzmann_prepare;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->prepare_if_needed = &_nc_hipert_boltzmann_prepare_if_needed;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_TT_Cls        = NULL;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_EE_Cls        = NULL;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_BB_Cls        = NULL;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_TE_Cls        = NULL;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_TB_Cls        = NULL;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_EB_Cls        = NULL;
}

static void
_nc_hipert_boltzmann_set_mode_k (NcHIPert *pert, gdouble k)
{
  NC_HIPERT_CLASS (nc_hipert_boltzmann_parent_class)->set_mode_k (pert, k);
  /* Chain up : start */
  {
    NcHIPertBoltzmann *pb = NC_HIPERT_BOLTZMANN (pert);
    ncm_model_ctrl_force_update (pb->ctrl_cosmo);
  }
}

static void
_nc_hipert_boltzmann_set_abstol (NcHIPert *pert, gdouble abstol)
{
  NC_HIPERT_CLASS (nc_hipert_boltzmann_parent_class)->set_abstol (pert, abstol);
  /* Chain up : start */
  {
    NcHIPertBoltzmann *pb = NC_HIPERT_BOLTZMANN (pert);
    ncm_model_ctrl_force_update (pb->ctrl_cosmo);
  }
}

static void
_nc_hipert_boltzmann_set_reltol (NcHIPert *pert, gdouble reltol)
{
  NC_HIPERT_CLASS (nc_hipert_boltzmann_parent_class)->set_reltol (pert, reltol);
  /* Chain up : start */
  {
    NcHIPertBoltzmann *pb = NC_HIPERT_BOLTZMANN (pert);
    ncm_model_ctrl_force_update (pb->ctrl_cosmo);
  }
}

static void
_nc_hipert_boltzmann_prepare (NcHIPertBoltzmann *pb, NcHIPrim *prim, NcHIReion *reion, NcHICosmo *cosmo)
{
  nc_recomb_prepare_if_needed (pb->recomb, reion, cosmo);

  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->init (pb, cosmo);
  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->reset (pb);
  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->set_opts (pb);
  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->evol (pb, pb->lambdaf);
}

static void
_nc_hipert_boltzmann_prepare_if_needed (NcHIPertBoltzmann *pb, NcHIPrim *prim, NcHIReion *reion, NcHICosmo *cosmo)
{
  gboolean cosmo_up = ncm_model_ctrl_update (pb->ctrl_cosmo, NCM_MODEL (cosmo));
  gboolean prim_up = ncm_model_ctrl_update (pb->ctrl_prim, NCM_MODEL (prim));

  if (cosmo_up || prim_up)
    nc_hipert_boltzmann_prepare (pb, prim, reion, cosmo);
}

/**
 * nc_hipert_boltzmann_ref:
 * @pb: a #NcHIPertBoltzmann.
 *
 * Increases the reference count of @pb.
 *
 * Returns: (transfer full): @pb.
 */
NcHIPertBoltzmann *
nc_hipert_boltzmann_ref (NcHIPertBoltzmann *pb)
{
  return g_object_ref (pb);
}

/**
 * nc_hipert_boltzmann_free:
 * @pb: a #NcHIPertBoltzmann.
 *
 * Decreases the reference count of @pb.
 *
 */
void
nc_hipert_boltzmann_free (NcHIPertBoltzmann *pb)
{
  g_object_unref (pb);
}

/**
 * nc_hipert_boltzmann_clear:
 * @pb: a #NcHIPertBoltzmann.
 *
 * Decreases the reference count of *@pb and sets *@pb to NULL.
 *
 */
void
nc_hipert_boltzmann_clear (NcHIPertBoltzmann **pb)
{
  g_clear_object (pb);
}

/**
 * nc_hipert_boltzmann_set_target_Cls:
 * @pb: a #NcHIPertBoltzmann
 * @tCls: Cls targets
 *
 * FIXME
 *
 */
void
nc_hipert_boltzmann_set_target_Cls (NcHIPertBoltzmann *pb, NcDataCMBDataType tCls)
{
  pb->target_Cls = tCls;
  ncm_model_ctrl_force_update (pb->ctrl_cosmo);
}

/**
 * nc_hipert_boltzmann_get_target_Cls:
 * @pb: a #NcHIPertBoltzmann.
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcDataCMBDataType
nc_hipert_boltzmann_get_target_Cls (NcHIPertBoltzmann *pb)
{
  return pb->target_Cls;
}

/**
 * nc_hipert_boltzmann_set_calc_transfer:
 * @pb: a #NcHIPertBoltzmann
 * @calc_transfer: a boolean
 *
 * FIXME
 *
 */
void
nc_hipert_boltzmann_set_calc_transfer (NcHIPertBoltzmann *pb, gboolean calc_transfer)
{
  pb->calc_transfer = calc_transfer;
  ncm_model_ctrl_force_update (pb->ctrl_cosmo);
}

/**
 * nc_hipert_boltzmann_get_calc_transfer:
 * @pb: a #NcHIPertBoltzmann.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_hipert_boltzmann_get_calc_transfer (NcHIPertBoltzmann *pb)
{
  return pb->calc_transfer;
}

/**
 * nc_hipert_boltzmann_set_lensed_Cls:
 * @pb: a #NcHIPertBoltzmann
 * @use_lensed_Cls: a boolean
 *
 * FIXME
 *
 */
void
nc_hipert_boltzmann_set_lensed_Cls (NcHIPertBoltzmann *pb, gboolean use_lensed_Cls)
{
  pb->use_lensed_Cls = use_lensed_Cls;
  ncm_model_ctrl_force_update (pb->ctrl_cosmo);
}

/**
 * nc_hipert_boltzmann_lensed_Cls:
 * @pb: a #NcHIPertBoltzmann.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_hipert_boltzmann_lensed_Cls (NcHIPertBoltzmann *pb)
{
  return pb->use_lensed_Cls;
}

/**
 * nc_hipert_boltzmann_set_tensor:
 * @pb: a #NcHIPertBoltzmann
 * @use_tensor: a boolean
 *
 * FIXME
 *
 */
void
nc_hipert_boltzmann_set_tensor (NcHIPertBoltzmann *pb, gboolean use_tensor)
{
  pb->use_tensor = use_tensor;
  ncm_model_ctrl_force_update (pb->ctrl_cosmo);
}

/**
 * nc_hipert_boltzmann_tensor:
 * @pb: a #NcHIPertBoltzmann.
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
nc_hipert_boltzmann_tensor (NcHIPertBoltzmann *pb)
{
  return pb->use_tensor;
}

/**
 * nc_hipert_boltzmann_set_TT_lmax:
 * @pb: a #NcHIPertBoltzmann.
 * @lmax: last mutipole.
 *
 * FIXME
 *
 */
void
nc_hipert_boltzmann_set_TT_lmax (NcHIPertBoltzmann *pb, guint lmax)
{
  g_assert_cmpuint (lmax, >=, 4);

  pb->TT_lmax = lmax;
  nc_hipert_set_sys_size (NC_HIPERT (pb), NC_HIPERT_BOLTZMANN_LEN + 2 * (lmax + 1 - 3));
  ncm_model_ctrl_force_update (pb->ctrl_cosmo);
}

/**
 * nc_hipert_boltzmann_set_EE_lmax:
 * @pb: a #NcHIPertBoltzmann.
 * @lmax: last mutipole.
 *
 * FIXME
 *
 */
void
nc_hipert_boltzmann_set_EE_lmax (NcHIPertBoltzmann *pb, guint lmax)
{
  g_assert_cmpuint (lmax, >=, 4);
  pb->EE_lmax = lmax;
  ncm_model_ctrl_force_update (pb->ctrl_cosmo);
}

/**
 * nc_hipert_boltzmann_set_BB_lmax:
 * @pb: a #NcHIPertBoltzmann.
 * @lmax: last mutipole.
 *
 * FIXME
 *
 */
void
nc_hipert_boltzmann_set_BB_lmax (NcHIPertBoltzmann *pb, guint lmax)
{
  g_assert_cmpuint (lmax, >=, 4);
  pb->BB_lmax = lmax;
  ncm_model_ctrl_force_update (pb->ctrl_cosmo);
}

/**
 * nc_hipert_boltzmann_set_TE_lmax:
 * @pb: a #NcHIPertBoltzmann.
 * @lmax: last mutipole.
 *
 * FIXME
 *
 */
void
nc_hipert_boltzmann_set_TE_lmax (NcHIPertBoltzmann *pb, guint lmax)
{
  g_assert_cmpuint (lmax, >=, 4);
  pb->TE_lmax = lmax;
  ncm_model_ctrl_force_update (pb->ctrl_cosmo);
}

/**
 * nc_hipert_boltzmann_set_TB_lmax:
 * @pb: a #NcHIPertBoltzmann.
 * @lmax: last mutipole.
 *
 * FIXME
 *
 */
void
nc_hipert_boltzmann_set_TB_lmax (NcHIPertBoltzmann *pb, guint lmax)
{
  g_assert_cmpuint (lmax, >=, 4);
  pb->TB_lmax = lmax;
  ncm_model_ctrl_force_update (pb->ctrl_cosmo);
}

/**
 * nc_hipert_boltzmann_set_EB_lmax:
 * @pb: a #NcHIPertBoltzmann.
 * @lmax: last mutipole.
 *
 * FIXME
 *
 */
void
nc_hipert_boltzmann_set_EB_lmax (NcHIPertBoltzmann *pb, guint lmax)
{
  g_assert_cmpuint (lmax, >=, 4);
  pb->EB_lmax = lmax;
  ncm_model_ctrl_force_update (pb->ctrl_cosmo);
}

/**
 * nc_hipert_boltzmann_get_TT_lmax:
 * @pb: a #NcHIPertBoltzmann.
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint nc_hipert_boltzmann_get_TT_lmax (NcHIPertBoltzmann *pb) { return pb->TT_lmax; }
/**
 * nc_hipert_boltzmann_get_EE_lmax:
 * @pb: a #NcHIPertBoltzmann.
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint nc_hipert_boltzmann_get_EE_lmax (NcHIPertBoltzmann *pb) { return pb->EE_lmax; }
/**
 * nc_hipert_boltzmann_get_BB_lmax:
 * @pb: a #NcHIPertBoltzmann.
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint nc_hipert_boltzmann_get_BB_lmax (NcHIPertBoltzmann *pb) { return pb->BB_lmax; }
/**
 * nc_hipert_boltzmann_get_TE_lmax:
 * @pb: a #NcHIPertBoltzmann.
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint nc_hipert_boltzmann_get_TE_lmax (NcHIPertBoltzmann *pb) { return pb->TE_lmax; }
/**
 * nc_hipert_boltzmann_get_TB_lmax:
 * @pb: a #NcHIPertBoltzmann.
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint nc_hipert_boltzmann_get_TB_lmax (NcHIPertBoltzmann *pb) { return pb->TB_lmax; }
/**
 * nc_hipert_boltzmann_get_EB_lmax:
 * @pb: a #NcHIPertBoltzmann.
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint nc_hipert_boltzmann_get_EB_lmax (NcHIPertBoltzmann *pb) { return pb->EB_lmax; }

/**
 * nc_hipert_boltzmann_set_recomb:
 * @pb: a #NcHIPertBoltzmann.
 * @recomb: a #NcRecomb.
 *
 * Sets the #NcRecomb object to be used in the Boltzmann evolution.
 *
 */
void
nc_hipert_boltzmann_set_recomb (NcHIPertBoltzmann *pb, NcRecomb *recomb)
{
  if (pb->recomb != recomb)
  {
    nc_recomb_clear (&pb->recomb);
    pb->recomb = nc_recomb_ref (recomb);
    NC_HIPERT (pb)->prepared = FALSE;
    ncm_model_ctrl_force_update (pb->ctrl_cosmo);
  }
}

/**
 * nc_hipert_boltzmann_prepare: (virtual prepare)
 * @pb: a #NcHIPertBoltzmann
 * @prim: a #NcHIPrim
 * @reion: a #NcHIReion
 * @cosmo: a NcHICosmo
 *
 * Prepares the #NcHIPertBoltzmann object.
 *
 */
void
nc_hipert_boltzmann_prepare (NcHIPertBoltzmann *pb, NcHIPrim *prim, NcHIReion *reion, NcHICosmo *cosmo)
{
  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->prepare (pb, prim, reion, cosmo);
}

/**
 * nc_hipert_boltzmann_prepare_if_needed: (virtual prepare_if_needed)
 * @pb: a #NcHIPertBoltzmann
 * @prim: a #NcHIPrim
 * @reion: a #NcHIReion
 * @cosmo: a #NcHICosmo
 *
 * Prepares the #NcHIPertBoltzmann object if not prepared for @cosmo.
 *
 */
void
nc_hipert_boltzmann_prepare_if_needed (NcHIPertBoltzmann *pb, NcHIPrim *prim, NcHIReion *reion, NcHICosmo *cosmo)
{
  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->prepare_if_needed (pb, prim, reion, cosmo);
}

/**
 * nc_hipert_boltzmann_get_TT_Cls: (virtual get_TT_Cls)
 * @pb: a #NcHIPertBoltzmann.
 * @Cls: a #NcmVector
 *
 * Prepares the #NcHIPertBoltzmann object.
 *
 */
void
nc_hipert_boltzmann_get_TT_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  guint lmax = ncm_vector_len (Cls);
  g_assert_cmpuint (lmax, >, 1);
  lmax--;

  if (!(pb->target_Cls & NC_DATA_CMB_TYPE_TT))
    g_error ("nc_hipert_boltzmann_get_TT_Cls: TT was not calculated, include it on the targets.");
  if (lmax > pb->TT_lmax)
    g_error ("nc_hipert_boltzmann_get_TT_Cls: TT was calculated up to ell = %u, but ell = %u was requested.",
             pb->TT_lmax, lmax);

  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->get_TT_Cls (pb, Cls);
}

/**
 * nc_hipert_boltzmann_get_EE_Cls: (virtual get_EE_Cls)
 * @pb: a #NcHIPertBoltzmann.
 * @Cls: a #NcmVector
 *
 * Prepares the #NcHIPertBoltzmann object.
 *
 */
void
nc_hipert_boltzmann_get_EE_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  guint lmax = ncm_vector_len (Cls);
  g_assert_cmpuint (lmax, >, 1);
  lmax--;

  if (!(pb->target_Cls & NC_DATA_CMB_TYPE_EE))
    g_error ("nc_hipert_boltzmann_get_EE_Cls: EE was not calculated, include it on the targets.");
  if (lmax > pb->EE_lmax)
    g_error ("nc_hipert_boltzmann_get_EE_Cls: EE was calculated up to ell = %u, but ell = %u was requested.",
             pb->EE_lmax, lmax);

  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->get_EE_Cls (pb, Cls);
}

/**
 * nc_hipert_boltzmann_get_BB_Cls: (virtual get_BB_Cls)
 * @pb: a #NcHIPertBoltzmann.
 * @Cls: a #NcmVector
 *
 * Prepares the #NcHIPertBoltzmann object.
 *
 */
void
nc_hipert_boltzmann_get_BB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  guint lmax = ncm_vector_len (Cls);
  g_assert_cmpuint (lmax, >, 1);
  lmax--;

  if (!(pb->target_Cls & NC_DATA_CMB_TYPE_BB))
    g_error ("nc_hipert_boltzmann_get_BB_Cls: BB was not calculated, include it on the targets.");
  if (lmax > pb->BB_lmax)
    g_error ("nc_hipert_boltzmann_get_BB_Cls: BB was calculated up to ell = %u, but ell = %u was requested.",
             pb->BB_lmax, lmax);

  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->get_BB_Cls (pb, Cls);
}

/**
 * nc_hipert_boltzmann_get_TE_Cls: (virtual get_TE_Cls)
 * @pb: a #NcHIPertBoltzmann
 * @Cls: a #NcmVector
 *
 * Prepares the #NcHIPertBoltzmann object.
 *
 */
void
nc_hipert_boltzmann_get_TE_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  guint lmax = ncm_vector_len (Cls);
  g_assert_cmpuint (lmax, >, 1);
  lmax--;

  if (!(pb->target_Cls & NC_DATA_CMB_TYPE_TE))
    g_error ("nc_hipert_boltzmann_get_TE_Cls: TE was not calculated, include it on the targets.");
  if (lmax > pb->TE_lmax)
    g_error ("nc_hipert_boltzmann_get_TE_Cls: TE was calculated up to ell = %u, but ell = %u was requested.",
             pb->TE_lmax, lmax);

  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->get_TE_Cls (pb, Cls);
}

/**
 * nc_hipert_boltzmann_get_TB_Cls: (virtual get_TB_Cls)
 * @pb: a #NcHIPertBoltzmann.
 * @Cls: a #NcmVector
 *
 * Prepares the #NcHIPertBoltzmann object.
 *
 */
void
nc_hipert_boltzmann_get_TB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  guint lmax = ncm_vector_len (Cls);
  g_assert_cmpuint (lmax, >, 1);
  lmax--;

  if (!(pb->target_Cls & NC_DATA_CMB_TYPE_TB))
    g_error ("nc_hipert_boltzmann_get_TB_Cls: TB was not calculated, include it on the targets.");
  if (lmax > pb->TB_lmax)
    g_error ("nc_hipert_boltzmann_get_TB_Cls: TB was calculated up to ell = %u, but ell = %u was requested.",
             pb->TB_lmax, lmax);

  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->get_TB_Cls (pb, Cls);
}

/**
 * nc_hipert_boltzmann_get_EB_Cls: (virtual get_EB_Cls)
 * @pb: a #NcHIPertBoltzmann
 * @Cls: a #NcmVector
 *
 * Prepares the #NcHIPertBoltzmann object.
 *
 */
void
nc_hipert_boltzmann_get_EB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  guint lmax = ncm_vector_len (Cls);
  g_assert_cmpuint (lmax, >, 1);
  lmax--;

  if (!(pb->target_Cls & NC_DATA_CMB_TYPE_EB))
    g_error ("nc_hipert_boltzmann_get_EB_Cls: EB was not calculated, include it on the targets.");
  if (lmax > pb->EB_lmax)
    g_error ("nc_hipert_boltzmann_get_EB_Cls: EB was calculated up to ell = %u, but ell = %u was requested.",
             pb->EB_lmax, lmax);

  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->get_EB_Cls (pb, Cls);
}
