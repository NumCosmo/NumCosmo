/***************************************************************************
 *            nc_bbn_parametrized.c
 *
 *  Sat August 30 10:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_bbn_parametrized.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * NcBBNParametrized:
 *
 * Primordial Helium as a free parameter.
 *
 * Reports $Y_p$ from its own #NcBBNParametrized:Yp parameter, predicting
 * nothing. This is the model for an analysis that constrains the Helium
 * abundance from the data rather than assuming Big Bang Nucleosynthesis --
 * the Planck "$Y_p$ free" runs, for instance -- and it is where the $Y_p$
 * parameter lives now that it is no longer a parameter of the background
 * models themselves.
 *
 * Being a parameter and not a prediction, it is defined for any cosmology:
 * nc_bbn_get_domain() reports the whole plane.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/bbn/nc_bbn_parametrized.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcBBNParametrized
{
  /*< private >*/
  NcBBN parent_instance;
};

G_DEFINE_TYPE (NcBBNParametrized, nc_bbn_parametrized, NC_TYPE_BBN)

#define VECTOR (NCM_MODEL (bbn_par))
#define YP_4HE (ncm_model_orig_param_get (VECTOR, NC_BBN_PARAMETRIZED_YP_4HE))

enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_bbn_parametrized_init (NcBBNParametrized *bbn_par)
{
}

static void
_nc_bbn_parametrized_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_bbn_parametrized_parent_class)->finalize (object);
}

static gdouble _nc_bbn_parametrized_Yp_4He (NcBBN *bbn, NcHICosmo *cosmo);
static void _nc_bbn_parametrized_get_domain (NcBBN *bbn, gdouble *wb_lb, gdouble *wb_ub, gdouble *DNeff_lb, gdouble *DNeff_ub);

static void
nc_bbn_parametrized_class_init (NcBBNParametrizedClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);
  NcBBNClass *bbn_class      = NC_BBN_CLASS (klass);

  object_class->finalize = &_nc_bbn_parametrized_finalize;

  ncm_model_class_set_name_nick (model_class, "Free primordial Helium.", "NcBBNParametrized");
  ncm_model_class_add_params (model_class, NC_BBN_PARAMETRIZED_SPARAM_LEN, 0, PROP_SIZE);

  ncm_model_class_set_sparam (model_class, NC_BBN_PARAMETRIZED_YP_4HE, "Y_p", "Yp",
                              0.0, 1.0, 1.0e-2,
                              NC_BBN_DEFAULT_PARAMS_ABSTOL, NC_BBN_PARAMETRIZED_DEFAULT_YP_4HE,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  ncm_model_class_add_impl_opts (model_class, NC_BBN_IMPL_Yp_4He, NC_BBN_IMPL_get_domain, -1);

  bbn_class->Yp_4He     = &_nc_bbn_parametrized_Yp_4He;
  bbn_class->get_domain = &_nc_bbn_parametrized_get_domain;
}

static gdouble
_nc_bbn_parametrized_Yp_4He (NcBBN *bbn, NcHICosmo *cosmo)
{
  NcBBNParametrized *bbn_par = NC_BBN_PARAMETRIZED (bbn);

  return YP_4HE;
}

static void
_nc_bbn_parametrized_get_domain (NcBBN *bbn, gdouble *wb_lb, gdouble *wb_ub, gdouble *DNeff_lb, gdouble *DNeff_ub)
{
  wb_lb[0]    = 0.0;
  wb_ub[0]    = GSL_POSINF;
  DNeff_lb[0] = GSL_NEGINF;
  DNeff_ub[0] = GSL_POSINF;
}

/**
 * nc_bbn_parametrized_new:
 *
 * Creates a new #NcBBNParametrized.
 *
 * Returns: (transfer full): a new #NcBBNParametrized.
 */
NcBBNParametrized *
nc_bbn_parametrized_new (void)
{
  return g_object_new (NC_TYPE_BBN_PARAMETRIZED, NULL);
}

/**
 * nc_bbn_parametrized_ref:
 * @bbn_par: a #NcBBNParametrized
 *
 * Increases the reference count of @bbn_par by one.
 *
 * Returns: (transfer full): @bbn_par.
 */
NcBBNParametrized *
nc_bbn_parametrized_ref (NcBBNParametrized *bbn_par)
{
  return g_object_ref (bbn_par);
}

/**
 * nc_bbn_parametrized_free:
 * @bbn_par: a #NcBBNParametrized
 *
 * Decreases the reference count of @bbn_par by one.
 *
 */
void
nc_bbn_parametrized_free (NcBBNParametrized *bbn_par)
{
  g_object_unref (bbn_par);
}

/**
 * nc_bbn_parametrized_clear:
 * @bbn_par: a #NcBBNParametrized
 *
 * If *@bbn_par is different from NULL, decreases the reference count of
 * *@bbn_par by one and sets *@bbn_par to NULL.
 *
 */
void
nc_bbn_parametrized_clear (NcBBNParametrized **bbn_par)
{
  g_clear_object (bbn_par);
}

