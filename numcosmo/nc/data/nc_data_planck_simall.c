/***************************************************************************
 *            nc_data_planck_simall.c
 *
 *  Wed July 23 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_data_planck_simall.c
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
 * NcDataPlanckSimall:
 *
 * Planck low-ℓ simall (SimAll) polarization likelihood, native port.
 *
 * A per-multipole tabulated likelihood: each theory $D_\ell = C_\ell
 * \ell(\ell+1)/2\pi$ (with $C_\ell$ divided by $A_\mathrm{planck}^2$) indexes a
 * precomputed log-probability table, $\ln L = \sum_\ell
 * \mathrm{prob}_\ell[\lfloor D_\ell/\mathrm{step}\rfloor]$, summed over the
 * present spectra (EE, BB, TE). This #NcmData reports $-2\ln L$.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_data_planck_simall.h"
#include "nc/cmb/nc_planck_fi.h"
#include "ncm/core/ncm_serialize.h"

#include <math.h>

typedef struct _NcSimallSpectrum
{
  NcmVector *prob; /* nell*nsteps flattened, or NULL if absent */
  gdouble step;
  guint nsteps;
  const gdouble *pdata;
  NcmVector *cl; /* Cl cache */
} NcSimallSpectrum;

struct _NcDataPlanckSimall
{
  /*< private >*/
  NcmData parent_instance;
  NcHIPertBoltzmann *pb;
  guint lmin;
  guint lmax;
  guint nell;
  gchar *calib_name;
  gint calib_pi;
  NcSimallSpectrum ee;
  NcSimallSpectrum bb;
  NcSimallSpectrum te;
};

enum
{
  PROP_0,
  PROP_PB,
  PROP_LMIN,
  PROP_LMAX,
  PROP_CALIB_NAME,
  PROP_STEP_EE,
  PROP_PROB_EE,
  PROP_STEP_BB,
  PROP_PROB_BB,
  PROP_STEP_TE,
  PROP_PROB_TE,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataPlanckSimall, nc_data_planck_simall, NCM_TYPE_DATA)

static void
nc_data_planck_simall_init (NcDataPlanckSimall *simall)
{
  simall->pb         = NULL;
  simall->lmin       = 0;
  simall->lmax       = 0;
  simall->nell       = 0;
  simall->calib_name = NULL;
  simall->calib_pi   = -1;
  simall->ee         = (NcSimallSpectrum) {
    NULL, 0.0, 0, NULL, NULL
  };
  simall->bb = (NcSimallSpectrum) {
    NULL, 0.0, 0, NULL, NULL
  };
  simall->te = (NcSimallSpectrum) {
    NULL, 0.0, 0, NULL, NULL
  };
}

static void
_nc_data_planck_simall_spectrum_setup (NcDataPlanckSimall *simall, NcSimallSpectrum *s)
{
  if (s->prob == NULL)
    return;

  g_assert_cmpuint (ncm_vector_len (s->prob) % simall->nell, ==, 0);
  s->nsteps = ncm_vector_len (s->prob) / simall->nell;
  g_assert_cmpuint (ncm_vector_stride (s->prob), ==, 1);
  s->pdata = ncm_vector_data (s->prob);
  s->cl    = ncm_vector_new (simall->lmax + 1);
}

static void
nc_data_planck_simall_constructed (GObject *object)
{
  NcDataPlanckSimall *simall = NC_DATA_PLANCK_SIMALL (object);

  g_assert_cmpuint (simall->lmax, >=, simall->lmin);
  simall->nell = simall->lmax - simall->lmin + 1;

  _nc_data_planck_simall_spectrum_setup (simall, &simall->ee);
  _nc_data_planck_simall_spectrum_setup (simall, &simall->bb);
  _nc_data_planck_simall_spectrum_setup (simall, &simall->te);

  g_assert ((simall->ee.prob != NULL) || (simall->bb.prob != NULL) || (simall->te.prob != NULL));

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_simall_parent_class)->constructed (object);
}

static void
nc_data_planck_simall_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataPlanckSimall *simall = NC_DATA_PLANCK_SIMALL (object);

  g_return_if_fail (NC_IS_DATA_PLANCK_SIMALL (object));

  switch (prop_id)
  {
    case PROP_PB:
      nc_data_planck_simall_set_hipert_boltzmann (simall, g_value_get_object (value));
      break;
    case PROP_LMIN:
      simall->lmin = g_value_get_uint (value);
      break;
    case PROP_LMAX:
      simall->lmax = g_value_get_uint (value);
      break;
    case PROP_CALIB_NAME:
      g_free (simall->calib_name);
      simall->calib_name = g_value_dup_string (value);
      break;
    case PROP_STEP_EE:
      simall->ee.step = g_value_get_double (value);
      break;
    case PROP_PROB_EE:
      ncm_vector_substitute (&simall->ee.prob, g_value_get_object (value), TRUE);
      break;
    case PROP_STEP_BB:
      simall->bb.step = g_value_get_double (value);
      break;
    case PROP_PROB_BB:
      ncm_vector_substitute (&simall->bb.prob, g_value_get_object (value), TRUE);
      break;
    case PROP_STEP_TE:
      simall->te.step = g_value_get_double (value);
      break;
    case PROP_PROB_TE:
      ncm_vector_substitute (&simall->te.prob, g_value_get_object (value), TRUE);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_planck_simall_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataPlanckSimall *simall = NC_DATA_PLANCK_SIMALL (object);

  g_return_if_fail (NC_IS_DATA_PLANCK_SIMALL (object));

  switch (prop_id)
  {
    case PROP_PB:
      g_value_set_object (value, simall->pb);
      break;
    case PROP_LMIN:
      g_value_set_uint (value, simall->lmin);
      break;
    case PROP_LMAX:
      g_value_set_uint (value, simall->lmax);
      break;
    case PROP_CALIB_NAME:
      g_value_set_string (value, simall->calib_name);
      break;
    case PROP_STEP_EE:
      g_value_set_double (value, simall->ee.step);
      break;
    case PROP_PROB_EE:
      g_value_set_object (value, simall->ee.prob);
      break;
    case PROP_STEP_BB:
      g_value_set_double (value, simall->bb.step);
      break;
    case PROP_PROB_BB:
      g_value_set_object (value, simall->bb.prob);
      break;
    case PROP_STEP_TE:
      g_value_set_double (value, simall->te.step);
      break;
    case PROP_PROB_TE:
      g_value_set_object (value, simall->te.prob);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
nc_data_planck_simall_dispose (GObject *object)
{
  NcDataPlanckSimall *simall = NC_DATA_PLANCK_SIMALL (object);

  nc_hipert_boltzmann_clear (&simall->pb);
  ncm_vector_clear (&simall->ee.prob);
  ncm_vector_clear (&simall->ee.cl);
  ncm_vector_clear (&simall->bb.prob);
  ncm_vector_clear (&simall->bb.cl);
  ncm_vector_clear (&simall->te.prob);
  ncm_vector_clear (&simall->te.cl);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_simall_parent_class)->dispose (object);
}

static void
nc_data_planck_simall_finalize (GObject *object)
{
  NcDataPlanckSimall *simall = NC_DATA_PLANCK_SIMALL (object);

  g_clear_pointer (&simall->calib_name, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_planck_simall_parent_class)->finalize (object);
}

static guint _nc_data_planck_simall_get_length (NcmData *data);
static void _nc_data_planck_simall_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_planck_simall_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL);

static void
nc_data_planck_simall_class_init (NcDataPlanckSimallClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);

  object_class->set_property = &nc_data_planck_simall_set_property;
  object_class->get_property = &nc_data_planck_simall_get_property;
  object_class->constructed  = &nc_data_planck_simall_constructed;
  object_class->dispose      = &nc_data_planck_simall_dispose;
  object_class->finalize     = &nc_data_planck_simall_finalize;

  g_object_class_install_property (object_class,
                                   PROP_PB,
                                   g_param_spec_object ("hipert-boltzmann",
                                                        NULL,
                                                        "Perturbations (Cls source)",
                                                        NC_TYPE_HIPERT_BOLTZMANN,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LMIN,
                                   g_param_spec_uint ("lmin", NULL, "Minimum multipole",
                                                      0, G_MAXUINT, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_LMAX,
                                   g_param_spec_uint ("lmax", NULL, "Maximum multipole",
                                                      0, G_MAXUINT, 29,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_CALIB_NAME,
                                   g_param_spec_string ("calib-name", NULL,
                                                        "Free-calibration parameter name (NULL => none)",
                                                        "A_planck",
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_STEP_EE,
                                   g_param_spec_double ("step-ee", NULL, "EE table Dl step",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PROB_EE,
                                   g_param_spec_object ("prob-ee", NULL,
                                                        "EE log-probability table (flattened nell*nsteps)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_STEP_BB,
                                   g_param_spec_double ("step-bb", NULL, "BB table Dl step",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PROB_BB,
                                   g_param_spec_object ("prob-bb", NULL,
                                                        "BB log-probability table (flattened nell*nsteps)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_STEP_TE,
                                   g_param_spec_double ("step-te", NULL, "TE table Dl step",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_PROB_TE,
                                   g_param_spec_object ("prob-te", NULL,
                                                        "TE log-probability table (flattened nell*nsteps)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->get_length = &_nc_data_planck_simall_get_length;
  data_class->prepare    = &_nc_data_planck_simall_prepare;
  data_class->m2lnL_val  = &_nc_data_planck_simall_m2lnL_val;
}

static guint
_nc_data_planck_simall_get_length (NcmData *data)
{
  NcDataPlanckSimall *simall = NC_DATA_PLANCK_SIMALL (data);
  guint n                    = 0;

  if (simall->ee.prob != NULL)
    n += simall->nell;

  if (simall->bb.prob != NULL)
    n += simall->nell;

  if (simall->te.prob != NULL)
    n += simall->nell;

  return n;
}

static void
_nc_data_planck_simall_prepare (NcmData *data, NcmMSet *mset)
{
  NcDataPlanckSimall *simall = NC_DATA_PLANCK_SIMALL (data);
  NcHICosmo *cosmo           = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));

  if (simall->pb == NULL)
    g_error ("_nc_data_planck_simall_prepare: no NcHIPertBoltzmann set.");

  if (cosmo == NULL)
    g_error ("_nc_data_planck_simall_prepare: no NcHICosmo in mset.");

  nc_hipert_boltzmann_prepare_if_needed (simall->pb, cosmo);
}

/* Sum the tabulated log-prob over @s for the given Cls (already in @s->cl),
 * with Dl = Cl/calib2 * l(l+1)/2pi and position = floor(Dl/step). */
static gdouble
_nc_data_planck_simall_spectrum_lnL (NcDataPlanckSimall *simall, NcSimallSpectrum *s, gdouble calib2, gboolean *bad)
{
  const gdouble twopi = 2.0 * M_PI;
  gdouble res         = 0.0;
  guint il;

  for (il = 0; il < simall->nell; il++)
  {
    const guint ell     = il + simall->lmin;
    const gdouble dl    = ncm_vector_get (s->cl, ell) / calib2 * (ell * (ell + 1.0)) / twopi;
    const gint position = (gint) (dl / s->step);

    if ((position < 0) || ((guint) position >= s->nsteps))
    {
      *bad = TRUE;

      return 0.0;
    }

    res += s->pdata[il * s->nsteps + position];
  }

  return res;
}

static void
_nc_data_planck_simall_m2lnL_val (NcmData *data, NcmMSet *mset, gdouble *m2lnL)
{
  NcDataPlanckSimall *simall = NC_DATA_PLANCK_SIMALL (data);
  gdouble calib2             = 1.0;
  gdouble lnL                = 0.0;
  gboolean bad               = FALSE;

  if (simall->calib_name != NULL)
  {
    NcPlanckFI *pfi = NC_PLANCK_FI (ncm_mset_peek (mset, nc_planck_fi_id ()));
    gdouble cal;

    if (pfi == NULL)
      g_error ("_nc_data_planck_simall_m2lnL_val: calibration '%s' requested but no NcPlanckFI in mset.", simall->calib_name);

    if (simall->calib_pi < 0)
    {
      guint idx;

      if (!ncm_model_param_index_from_name (NCM_MODEL (pfi), simall->calib_name, &idx, NULL))
        g_error ("_nc_data_planck_simall_m2lnL_val: calibration parameter '%s' not in model.", simall->calib_name);

      simall->calib_pi = idx;
    }

    cal    = ncm_model_param_get (NCM_MODEL (pfi), simall->calib_pi);
    calib2 = cal * cal;
  }

  if (simall->ee.prob != NULL)
  {
    nc_hipert_boltzmann_get_EE_Cls (simall->pb, simall->ee.cl);
    lnL += _nc_data_planck_simall_spectrum_lnL (simall, &simall->ee, calib2, &bad);
  }

  if (!bad && (simall->bb.prob != NULL))
  {
    nc_hipert_boltzmann_get_BB_Cls (simall->pb, simall->bb.cl);
    lnL += _nc_data_planck_simall_spectrum_lnL (simall, &simall->bb, calib2, &bad);
  }

  if (!bad && (simall->te.prob != NULL))
  {
    nc_hipert_boltzmann_get_TE_Cls (simall->pb, simall->te.cl);
    lnL += _nc_data_planck_simall_spectrum_lnL (simall, &simall->te, calib2, &bad);
  }

  if (bad)
    *m2lnL = 1.0e30;
  else
    *m2lnL = -2.0 * lnL;
}

/**
 * nc_data_planck_simall_new:
 * @lmin: lowest multipole
 * @lmax: highest multipole
 * @step_ee: EE table Dl step (0 if absent)
 * @prob_ee: (transfer none) (nullable): EE log-prob table (flattened nell*nsteps)
 * @step_bb: BB table Dl step (0 if absent)
 * @prob_bb: (transfer none) (nullable): BB log-prob table
 * @step_te: TE table Dl step (0 if absent)
 * @prob_te: (transfer none) (nullable): TE log-prob table
 *
 * Builds a #NcDataPlanckSimall (calibration defaults to A_planck). The theory
 * $C_\ell$ source is set separately.
 *
 * Returns: (transfer full): the new #NcDataPlanckSimall.
 */
NcDataPlanckSimall *
nc_data_planck_simall_new (guint     lmin,
                           guint     lmax,
                           gdouble   step_ee,
                           NcmVector *prob_ee,
                           gdouble   step_bb,
                           NcmVector *prob_bb,
                           gdouble   step_te,
                           NcmVector *prob_te)
{
  return g_object_new (NC_TYPE_DATA_PLANCK_SIMALL,
                       "lmin", lmin,
                       "lmax", lmax,
                       "step-ee", step_ee,
                       "prob-ee", prob_ee,
                       "step-bb", step_bb,
                       "prob-bb", prob_bb,
                       "step-te", step_te,
                       "prob-te", prob_te,
                       NULL);
}

/**
 * nc_data_planck_simall_new_from_file:
 * @filename: a serialized #NcDataPlanckSimall file.
 *
 * Returns: (transfer full): the new #NcDataPlanckSimall.
 */
NcDataPlanckSimall *
nc_data_planck_simall_new_from_file (const gchar *filename)
{
  NcDataPlanckSimall *simall = NC_DATA_PLANCK_SIMALL (ncm_serialize_global_from_file (filename));

  g_assert (NC_IS_DATA_PLANCK_SIMALL (simall));

  return simall;
}

/**
 * nc_data_planck_simall_set_hipert_boltzmann:
 * @simall: a #NcDataPlanckSimall
 * @pb: a #NcHIPertBoltzmann
 *
 * Sets the theory $C_\ell$ source.
 */
void
nc_data_planck_simall_set_hipert_boltzmann (NcDataPlanckSimall *simall, NcHIPertBoltzmann *pb)
{
  nc_hipert_boltzmann_clear (&simall->pb);

  if (pb != NULL)
    simall->pb = nc_hipert_boltzmann_ref (pb);
}

/**
 * nc_data_planck_simall_peek_hipert_boltzmann:
 * @simall: a #NcDataPlanckSimall
 *
 * Returns: (transfer none): the perturbations object, or %NULL.
 */
NcHIPertBoltzmann *
nc_data_planck_simall_peek_hipert_boltzmann (NcDataPlanckSimall *simall)
{
  return simall->pb;
}

