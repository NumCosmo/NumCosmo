/***************************************************************************
 *            nc_bbn.c
 *
 *  Fri August 29 20:30:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_bbn.c
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
 * NcBBN:
 *
 * Abstract class for primordial nucleosynthesis models.
 *
 * A #NcBBN answers what Big Bang Nucleosynthesis predicts for a cosmology,
 * given the two quantities BBN actually depends on: the baryon density
 * $\omega_b = \Omega_{b0}h^2$ and the extra relativistic degrees of freedom
 * $\Delta N_\mathrm{eff} = N_\mathrm{eff} - 3.046$. Both are read from the
 * #NcHICosmo passed to each method, never stored: a #NcBBN is a submodel of
 * #NcHICosmo and holding a reference back to its parent would be a cycle.
 *
 * The only prediction every implementation must provide is the Helium-4 mass
 * fraction $Y_p$, nc_bbn_Yp_4He(). The tables shipped with NumCosmo tabulate
 * $Y_p$ alone, so nc_bbn_DH(), nc_bbn_He3H() and nc_bbn_Li7H() are optional --
 * check ncm_model_check_impl_opt() against #NcBBNImpl before calling them.
 *
 * Implementations that interpolate a table are valid only inside it.
 * nc_bbn_get_domain() reports that range and nc_bbn_check_domain() turns a
 * cosmology outside it into an error: a table asked about a baryon density it
 * never covered has no answer, and returning an extrapolated one silently is
 * how a chain ends up converged on invented helium.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/bbn/nc_bbn.h"
#include "ncm/core/ncm_cfg.h"

typedef struct _NcBBNPrivate
{
  gint placeholder;
} NcBBNPrivate;

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcBBN, nc_bbn, NCM_TYPE_MODEL)

enum
{
  PROP_0,
  PROP_SIZE,
};

static void
nc_bbn_init (NcBBN *bbn)
{
  NcBBNPrivate * const self = nc_bbn_get_instance_private (bbn);

  self->placeholder = 0;
}

static void
nc_bbn_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_bbn_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_bbn, NC_TYPE_BBN);

static gdouble _nc_bbn_Yp_4He (NcBBN *bbn, NcHICosmo *cosmo);
static void _nc_bbn_get_domain (NcBBN *bbn, gdouble *wb_lb, gdouble *wb_ub, gdouble *DNeff_lb, gdouble *DNeff_ub);
static gdouble _nc_bbn_DH (NcBBN *bbn, NcHICosmo *cosmo);
static gdouble _nc_bbn_He3H (NcBBN *bbn, NcHICosmo *cosmo);
static gdouble _nc_bbn_Li7H (NcBBN *bbn, NcHICosmo *cosmo);

static void
nc_bbn_class_init (NcBBNClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  object_class->finalize = nc_bbn_finalize;

  ncm_model_class_set_name_nick (model_class, "Abstract class for primordial nucleosynthesis models.", "NcBBN");
  ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);

  ncm_mset_model_register_id (model_class,
                              "NcBBN",
                              "Primordial nucleosynthesis models.",
                              NULL,
                              FALSE,
                              nc_hicosmo_id ());

  ncm_model_class_check_params_info (model_class);

  klass->Yp_4He     = &_nc_bbn_Yp_4He;
  klass->get_domain = &_nc_bbn_get_domain;
  klass->DH         = &_nc_bbn_DH;
  klass->He3H       = &_nc_bbn_He3H;
  klass->Li7H       = &_nc_bbn_Li7H;
}

/* LCOV_EXCL_START */

static gdouble
_nc_bbn_Yp_4He (NcBBN *bbn, NcHICosmo *cosmo)
{
  g_error ("_nc_bbn_Yp_4He: object `%s' does not implement this virtual function.",
           g_type_name (G_OBJECT_TYPE (bbn)));

  return 0.0;
}

static void
_nc_bbn_get_domain (NcBBN *bbn, gdouble *wb_lb, gdouble *wb_ub, gdouble *DNeff_lb, gdouble *DNeff_ub)
{
  g_error ("_nc_bbn_get_domain: object `%s' does not implement this virtual function.",
           g_type_name (G_OBJECT_TYPE (bbn)));
}

static gdouble
_nc_bbn_DH (NcBBN *bbn, NcHICosmo *cosmo)
{
  g_error ("_nc_bbn_DH: object `%s' does not predict the Deuterium abundance.",
           g_type_name (G_OBJECT_TYPE (bbn)));

  return 0.0;
}

static gdouble
_nc_bbn_He3H (NcBBN *bbn, NcHICosmo *cosmo)
{
  g_error ("_nc_bbn_He3H: object `%s' does not predict the Helium-3 abundance.",
           g_type_name (G_OBJECT_TYPE (bbn)));

  return 0.0;
}

static gdouble
_nc_bbn_Li7H (NcBBN *bbn, NcHICosmo *cosmo)
{
  g_error ("_nc_bbn_Li7H: object `%s' does not predict the Lithium-7 abundance.",
           g_type_name (G_OBJECT_TYPE (bbn)));

  return 0.0;
}

/* LCOV_EXCL_STOP */

/**
 * nc_bbn_ref:
 * @bbn: a #NcBBN
 *
 * Increases the reference count of @bbn by one.
 *
 * Returns: (transfer full): @bbn.
 */
NcBBN *
nc_bbn_ref (NcBBN *bbn)
{
  return g_object_ref (bbn);
}

/**
 * nc_bbn_free:
 * @bbn: a #NcBBN
 *
 * Decreases the reference count of @bbn by one.
 *
 */
void
nc_bbn_free (NcBBN *bbn)
{
  g_object_unref (bbn);
}

/**
 * nc_bbn_clear:
 * @bbn: a #NcBBN
 *
 * If *@bbn is different from NULL, decreases the reference count of
 * *@bbn by one and sets *@bbn to NULL.
 *
 */
void
nc_bbn_clear (NcBBN **bbn)
{
  g_clear_object (bbn);
}

/**
 * nc_bbn_Yp_4He: (virtual Yp_4He)
 * @bbn: a #NcBBN
 * @cosmo: a #NcHICosmo
 *
 * Computes the primordial Helium-4 mass fraction $Y_p$ predicted for @cosmo.
 *
 * Returns: $Y_p$.
 */
gdouble
nc_bbn_Yp_4He (NcBBN *bbn, NcHICosmo *cosmo)
{
  return NC_BBN_GET_CLASS (bbn)->Yp_4He (bbn, cosmo);
}

/**
 * nc_bbn_get_domain: (virtual get_domain)
 * @bbn: a #NcBBN
 * @wb_lb: (out): lower bound on $\omega_b$
 * @wb_ub: (out): upper bound on $\omega_b$
 * @DNeff_lb: (out): lower bound on $\Delta N_\mathrm{eff}$
 * @DNeff_ub: (out): upper bound on $\Delta N_\mathrm{eff}$
 *
 * Gets the range of $(\omega_b, \Delta N_\mathrm{eff})$ over which @bbn is
 * defined. A model with no intrinsic range reports the widest one it accepts.
 *
 */
void
nc_bbn_get_domain (NcBBN *bbn, gdouble *wb_lb, gdouble *wb_ub, gdouble *DNeff_lb, gdouble *DNeff_ub)
{
  NC_BBN_GET_CLASS (bbn)->get_domain (bbn, wb_lb, wb_ub, DNeff_lb, DNeff_ub);
}

/**
 * nc_bbn_DH: (virtual DH)
 * @bbn: a #NcBBN
 * @cosmo: a #NcHICosmo
 *
 * Computes the primordial Deuterium abundance $\mathrm{D}/\mathrm{H}$ predicted
 * for @cosmo. Optional: check ncm_model_check_impl_opt() for %NC_BBN_IMPL_DH first.
 *
 * Returns: $\mathrm{D}/\mathrm{H}$.
 */
gdouble
nc_bbn_DH (NcBBN *bbn, NcHICosmo *cosmo)
{
  return NC_BBN_GET_CLASS (bbn)->DH (bbn, cosmo);
}

/**
 * nc_bbn_He3H: (virtual He3H)
 * @bbn: a #NcBBN
 * @cosmo: a #NcHICosmo
 *
 * Computes the primordial Helium-3 abundance ${}^3\mathrm{He}/\mathrm{H}$
 * predicted for @cosmo. Optional: check ncm_model_check_impl_opt() for
 * %NC_BBN_IMPL_He3H first.
 *
 * Returns: ${}^3\mathrm{He}/\mathrm{H}$.
 */
gdouble
nc_bbn_He3H (NcBBN *bbn, NcHICosmo *cosmo)
{
  return NC_BBN_GET_CLASS (bbn)->He3H (bbn, cosmo);
}

/**
 * nc_bbn_Li7H: (virtual Li7H)
 * @bbn: a #NcBBN
 * @cosmo: a #NcHICosmo
 *
 * Computes the primordial Lithium-7 abundance ${}^7\mathrm{Li}/\mathrm{H}$
 * predicted for @cosmo. Optional: check ncm_model_check_impl_opt() for
 * %NC_BBN_IMPL_Li7H first.
 *
 * Returns: ${}^7\mathrm{Li}/\mathrm{H}$.
 */
gdouble
nc_bbn_Li7H (NcBBN *bbn, NcHICosmo *cosmo)
{
  return NC_BBN_GET_CLASS (bbn)->Li7H (bbn, cosmo);
}

/**
 * nc_bbn_check_domain:
 * @bbn: a #NcBBN
 * @wb: the baryon density $\omega_b = \Omega_{b0}h^2$
 * @DNeff: the extra relativistic degrees of freedom $\Delta N_\mathrm{eff}$
 *
 * Errors out unless $(@wb, @DNeff)$ lies inside the domain @bbn reports through
 * nc_bbn_get_domain(). Implementations that interpolate a table call this
 * before evaluating: outside the table there is no prediction to give, and an
 * extrapolated value would be indistinguishable from a real one.
 *
 */
void
nc_bbn_check_domain (NcBBN *bbn, const gdouble wb, const gdouble DNeff)
{
  gdouble wb_lb, wb_ub, DNeff_lb, DNeff_ub;

  nc_bbn_get_domain (bbn, &wb_lb, &wb_ub, &DNeff_lb, &DNeff_ub);

  if ((wb < wb_lb) || (wb > wb_ub))
    g_error ("nc_bbn_check_domain: `%s' is defined for omega_b in [% 20.15g, % 20.15g], got % 20.15g.",
             g_type_name (G_OBJECT_TYPE (bbn)), wb_lb, wb_ub, wb);

  if ((DNeff < DNeff_lb) || (DNeff > DNeff_ub))
    g_error ("nc_bbn_check_domain: `%s' is defined for DeltaNeff in [% 20.15g, % 20.15g], got % 20.15g.",
             g_type_name (G_OBJECT_TYPE (bbn)), DNeff_lb, DNeff_ub, DNeff);
}

