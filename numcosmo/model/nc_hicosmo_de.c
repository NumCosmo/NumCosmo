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
 * @title: Darkenergy Virtual Model
 * @short_description: Generic model for darkenergy parametric models
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "model/nc_hicosmo_de.h"
#include "math/ncm_reparam_linear.h"

G_DEFINE_ABSTRACT_TYPE (NcHICosmoDE, nc_hicosmo_de, NC_TYPE_HICOSMO);

#define VECTOR    (model->params)
#define MACRO_H0  (ncm_vector_get (VECTOR, NC_HICOSMO_DE_H0))
#define OMEGA_C   (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_C))
#define OMEGA_X   (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_X))
#define T_GAMMA0  (ncm_vector_get (VECTOR, NC_HICOSMO_DE_T_GAMMA0))
#define OMEGA_R   nc_hicosmo_Omega_r (NC_HICOSMO (model))
#define OMEGA_B   (ncm_vector_get (VECTOR, NC_HICOSMO_DE_OMEGA_B))
#define SPECINDEX (ncm_vector_get (VECTOR, NC_HICOSMO_DE_SPECINDEX))
#define SIGMA8    (ncm_vector_get (VECTOR, NC_HICOSMO_DE_SIGMA8))

#define OMEGA_M (OMEGA_B + OMEGA_C)
#define OMEGA_K (1.0 - (OMEGA_B + OMEGA_C + OMEGA_R + OMEGA_X))

/****************************************************************************
 * Normalized Hubble function
 ****************************************************************************/

static gdouble
_nc_hicosmo_de_E2 (NcmModel *model, gdouble z)
{
  const gdouble omega_k = OMEGA_K;
  const gdouble x = 1.0 + z;
  const gdouble x2 = x * x;
  const gdouble x3 = x2 * x;
  const gdouble x4 = x3 * x;
  return (OMEGA_R * x4 + OMEGA_M * x3 + omega_k * x2 + nc_hicosmo_de_weff (NC_HICOSMO_DE (model), z));
}

/****************************************************************************
 * dE_dz
 ****************************************************************************/

static gdouble
_nc_hicosmo_de_dE2_dz (NcmModel *model, gdouble z)
{
  const gdouble omega_k = OMEGA_K;
  const gdouble x = 1.0 + z;
  const gdouble x2 = x*x;
  const gdouble x3 = x2*x;

  return (4.0 * OMEGA_R * x3 + 3.0 * OMEGA_M * x2 + 2.0 * omega_k * x + nc_hicosmo_de_dweff_dz (NC_HICOSMO_DE (model), z));
}

/****************************************************************************
 * Simple functions
 ****************************************************************************/
static gdouble _nc_hicosmo_de_H0 (NcmModel *model) { return MACRO_H0; }
static gdouble _nc_hicosmo_de_Omega_t (NcmModel *model) { return OMEGA_M + OMEGA_X + OMEGA_R; }
static gdouble _nc_hicosmo_de_Omega_c (NcmModel *model) { return OMEGA_C; }
static gdouble _nc_hicosmo_de_T_gamma0 (NcmModel *model) { return T_GAMMA0; }
static gdouble
_nc_hicosmo_de_Omega_r (NcmModel *model)
{
  const gdouble h = MACRO_H0 / 100.0;
  const gdouble h2 = h * h;
  return (1.0 + 0.2271 * NC_C_NEUTRINO_N_EFF) * NC_C_RADIATION_TEMP_TO_h2OMEGA_R (T_GAMMA0) / h2;
}
static gdouble _nc_hicosmo_de_Omega_b (NcmModel *model) { return OMEGA_B; }
static gdouble _nc_hicosmo_de_sigma_8 (NcmModel *model) { return SIGMA8; }
static gdouble _nc_hicosmo_de_powspec (NcmModel *model, gdouble k) { return pow (k, SPECINDEX); }

/**
 * FIXME
 */
void
nc_hicosmo_de_omega_x2omega_k (NcHICosmo *model)
{
  guint size = ncm_model_len (NCM_MODEL (model));
  NcmMatrix *T = ncm_matrix_new (size, size);
  NcmVector *v = ncm_vector_new (size);
  NcmReparamLinear *relin;

  ncm_matrix_set_identity (T);
  ncm_vector_set_zero (v);

  ncm_matrix_set (T, NC_HICOSMO_DE_OMEGA_X, NC_HICOSMO_DE_OMEGA_C, -1.0);
  ncm_matrix_set (T, NC_HICOSMO_DE_OMEGA_X, NC_HICOSMO_DE_OMEGA_X, -1.0);
  ncm_matrix_set (T, NC_HICOSMO_DE_OMEGA_X, NC_HICOSMO_DE_OMEGA_B, -1.0);

  ncm_vector_set (v, NC_HICOSMO_DE_OMEGA_X, 1.0);

  relin = ncm_reparam_linear_new (ncm_model_len (NCM_MODEL (model)), T, v);
  {
	NcmReparam *reparam = NCM_REPARAM (relin);
	NcmSParam *sp_Omegak = ncm_sparam_new ("Omegak","\\Omega_k", -5.0, 5.0, 1.0e-2, NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, 0.0, NCM_PARAM_TYPE_FIXED);
	g_ptr_array_index (reparam->sparams, NC_HICOSMO_DE_OMEGA_X) = sp_Omegak;
  }

  ncm_model_set_reparam (NCM_MODEL (model), NCM_REPARAM (relin));

  return;
}

static void
bbn_prior (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f)
{
  NcmModel *model = ncm_mset_peek (mset, NC_HICOSMO_ID);
  gdouble z_bbn = 1.0e9;
  gdouble bbn, a;

  a = nc_hicosmo_de_weff (NC_HICOSMO_DE (model), z_bbn) / nc_hicosmo_E2 (NC_HICOSMO (model), z_bbn);
  bbn = 1.0 / sqrt(1.0 - a);
  f[0] = (bbn - 0.942) / 0.03;
}

/**
 * FIXME
 */
gboolean
nc_hicosmo_de_new_add_bbn (NcLikelihood *lh)
{
  NcmMSetFunc *func = ncm_mset_func_new (bbn_prior, 0, 1, NULL, NULL);
  nc_likelihood_priors_add (lh, func);
  return TRUE;
}

enum {
  PROP_0,
  PROP_SIZE,
};

static void
nc_hicosmo_de_init (NcHICosmoDE *object)
{
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
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcHICosmoClass* parent_class = NC_HICOSMO_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);
  object_class->set_property = &ncm_model_class_set_property;
  object_class->get_property = &ncm_model_class_get_property;
  object_class->finalize     = &nc_hicosmo_de_finalize;

  ncm_model_class_add_params (model_class, 7, 0, PROP_SIZE);
  /* Set H_0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_H0, "H_0", "H0",
                               10.0, 500.0, 1.0,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_H0,
                               NCM_PARAM_TYPE_FIXED);
  /* Set Omega_c param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_C, "\\Omega_c", "Omegac",
                               1e-8,  10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_OMEGA_C,
                               NCM_PARAM_TYPE_FREE);
  /* Set Omega_x param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_X, "\\Omega_x", "Omegax",
                               1e-8,  10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_OMEGA_X,
                               NCM_PARAM_TYPE_FREE);
  /* Set T_gamma0 param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_T_GAMMA0, "T_{\\gamma0}", "Tgamma0",
                               1e-8,  10.0, 1.0e-2,
                               NC_HICOSMO_DEFAULT_PARAMS_ABSTOL, NC_HICOSMO_DE_DEFAULT_T_GAMMA0,
                               NCM_PARAM_TYPE_FIXED);
  /* Set Omega_b param info */
  ncm_model_class_set_sparam (model_class, NC_HICOSMO_DE_OMEGA_B, "\\Omega_b", "Omegab",
                               1e-8,  10.0, 1.0e-2,
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
  nc_hicosmo_set_Omega_c_impl   (parent_class, &_nc_hicosmo_de_Omega_c);
  nc_hicosmo_set_Omega_r_impl   (parent_class, &_nc_hicosmo_de_Omega_r);
  nc_hicosmo_set_Omega_b_impl   (parent_class, &_nc_hicosmo_de_Omega_b);
  nc_hicosmo_set_Omega_t_impl   (parent_class, &_nc_hicosmo_de_Omega_t);
  nc_hicosmo_set_sigma_8_impl   (parent_class, &_nc_hicosmo_de_sigma_8);
  nc_hicosmo_set_T_gamma0_impl  (parent_class, &_nc_hicosmo_de_T_gamma0);
  nc_hicosmo_set_dE2_dz_impl    (parent_class, &_nc_hicosmo_de_dE2_dz);
  nc_hicosmo_set_powspec_impl   (parent_class, &_nc_hicosmo_de_powspec);
}

#define NC_HICOSMO_DE_SET_IMPL_FUNC(name) \
void \
nc_hicosmo_de_set_##name##_impl (NcHICosmoDEClass *model_de_class, NcmFuncF f, NcmFuncPF pf, NcmFuncDF df) \
{ \
  NCM_MODEL_CLASS (model_de_class)->impl |= NC_HICOSMO_DE_IMPL_##name; \
  g_assert (f != NULL); \
  model_de_class->name = *ncm_func_stub; \
  model_de_class->name.f = f; \
  model_de_class->name.pf = pf; \
  model_de_class->name.df = df; \
}

/**
 * nc_hicosmo_de_set_weff_impl: (skip)
 * @model_de_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO_DE,NcHICosmoDE,nc_hicosmo_de,NcmModelFunc1,weff)
/**
 * nc_hicosmo_de_set_dweff_dz_impl: (skip)
 * @model_de_class: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 */
NCM_MODEL_SET_IMPL_FUNC(NC_HICOSMO_DE,NcHICosmoDE,nc_hicosmo_de,NcmModelFunc1,dweff_dz)
/**
 * nc_hicosmo_weff:
 * @model: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_weff_pf:
 * @model: FIXME
 * @x: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_weff_df:
 * @model: FIXME
 * @pt: FIXME
 * @x: FIXME
 * @v: FIXME
 *
 * FIXME
 *
 */
/**
 * nc_hicosmo_dweff_dz:
 * @model: FIXME
 * @x: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_dweff_dz_pf:
 * @model: FIXME
 * @x: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
/**
 * nc_hicosmo_dweff_dz_df:
 * @model: FIXME
 * @pt: FIXME
 * @x: FIXME
 * @v: FIXME
 *
 * FIXME
 *
 */
