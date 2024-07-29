/***************************************************************************
 *            nc_multiplicity_func_st.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_multiplicity_func_st
 * @title: NcMultiplicityFuncST
 * @short_description: Dark matter halo -- Sheth-Tormen multiplicity function.
 *
 * Computes the multiplicity function of dark matter halos using the
 * Sheth-Tormen model.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_st.h"

struct _NcMultiplicityFuncSTPrivate
{
  NcMultiplicityFuncMassDef mdef;
  gdouble A;
  gdouble b;
  gdouble p;
  gdouble delta_c;
  gdouble Delta;
};

enum
{
  PROP_0,
  PROP_A,
  PROP_B,
  PROP_P,
  PROP_DELTA_C,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcMultiplicityFuncST, nc_multiplicity_func_st, NC_TYPE_MULTIPLICITY_FUNC)

static void
nc_multiplicity_func_st_init (NcMultiplicityFuncST *mst)
{
  NcMultiplicityFuncSTPrivate * const self = mst->priv = nc_multiplicity_func_st_get_instance_private (mst);

  self->mdef    = NC_MULTIPLICITY_FUNC_MASS_DEF_LEN;
  self->A       = 0.0;
  self->b       = 0.0;
  self->p       = 0.0;
  self->delta_c = 0.0;
  self->Delta   = 0.0;
}

static void
_nc_multiplicity_func_st_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncST *mst = NC_MULTIPLICITY_FUNC_ST (object);

  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_ST (object));

  switch (prop_id)
  {
    case PROP_A:
      nc_multiplicity_func_st_set_A (mst, g_value_get_double (value));
      break;
    case PROP_B:
      nc_multiplicity_func_st_set_b (mst, g_value_get_double (value));
      break;
    case PROP_P:
      nc_multiplicity_func_st_set_p (mst, g_value_get_double (value));
      break;
    case PROP_DELTA_C:
      nc_multiplicity_func_st_set_delta_c (mst, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_st_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcMultiplicityFuncST *mst = NC_MULTIPLICITY_FUNC_ST (object);

  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_ST (object));

  switch (prop_id)
  {
    case PROP_A:
      g_value_set_double (value, nc_multiplicity_func_st_get_A (mst));
      break;
    case PROP_B:
      g_value_set_double (value, nc_multiplicity_func_st_get_b (mst));
      break;
    case PROP_P:
      g_value_set_double (value, nc_multiplicity_func_st_get_p (mst));
      break;
    case PROP_DELTA_C:
      g_value_set_double (value, nc_multiplicity_func_st_get_delta_c (mst));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_st_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_st_parent_class)->finalize (object);
}

static void _nc_multiplicity_func_st_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef);
static void _nc_multiplicity_func_st_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta);
static NcMultiplicityFuncMassDef _nc_multiplicity_func_st_get_mdef (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_st_get_Delta (NcMultiplicityFunc *mulf);
static gdouble _nc_multiplicity_func_st_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z);

static void
nc_multiplicity_func_st_class_init (NcMultiplicityFuncSTClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass *parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  object_class->set_property = _nc_multiplicity_func_st_set_property;
  object_class->get_property = _nc_multiplicity_func_st_get_property;
  object_class->finalize     = _nc_multiplicity_func_st_finalize;

  /**
   * NcMultiplicityFuncST:A:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A,
                                   g_param_spec_double ("A",
                                                        NULL,
                                                        "A",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.3222,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcMultiplicityFuncST:b:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_B,
                                   g_param_spec_double ("b",
                                                        NULL,
                                                        "b",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.707,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcMultiplicityFuncST:p:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_P,
                                   g_param_spec_double ("p",
                                                        NULL,
                                                        "p",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.3,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcMultiplicityFuncST:critical_delta:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_DELTA_C,
                                   g_param_spec_double ("critical-delta",
                                                        NULL,
                                                        "Critical delta",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, NC_MULTIPLICITY_FUNC_DELTA_C0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  parent_class->set_mdef  = &_nc_multiplicity_func_st_set_mdef;
  parent_class->get_mdef  = &_nc_multiplicity_func_st_get_mdef;
  parent_class->set_Delta = &_nc_multiplicity_func_st_set_Delta;
  parent_class->get_Delta = &_nc_multiplicity_func_st_get_Delta;
  parent_class->eval      = &_nc_multiplicity_func_st_eval;
}

static void
_nc_multiplicity_func_st_set_mdef (NcMultiplicityFunc *mulf, NcMultiplicityFuncMassDef mdef)
{
  NcMultiplicityFuncST *mst                = NC_MULTIPLICITY_FUNC_ST (mulf);
  NcMultiplicityFuncSTPrivate * const self = mst->priv;

  switch (mdef)
  {
    case NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN:
      /* nothing to do */
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL:
      g_error ("NcMultiplicityFuncST does not support critical mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL:
      g_error ("NcMultiplicityFuncST does not support virial mass def");
      break;
    case NC_MULTIPLICITY_FUNC_MASS_DEF_FOF:
      g_error ("NcMultiplicityFuncST does not support fof mass def");
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  self->mdef = mdef;
}

static NcMultiplicityFuncMassDef
_nc_multiplicity_func_st_get_mdef (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncST *mst                = NC_MULTIPLICITY_FUNC_ST (mulf);
  NcMultiplicityFuncSTPrivate * const self = mst->priv;

  return self->mdef;
}

static void
_nc_multiplicity_func_st_set_Delta (NcMultiplicityFunc *mulf, gdouble Delta)
{
  NcMultiplicityFuncST *mst                = NC_MULTIPLICITY_FUNC_ST (mulf);
  NcMultiplicityFuncSTPrivate * const self = mst->priv;

  self->Delta = Delta;
}

static gdouble
_nc_multiplicity_func_st_get_Delta (NcMultiplicityFunc *mulf)
{
  NcMultiplicityFuncST *mst                = NC_MULTIPLICITY_FUNC_ST (mulf);
  NcMultiplicityFuncSTPrivate * const self = mst->priv;

  return self->Delta;
}

static gdouble
_nc_multiplicity_func_st_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z) /* f(\sigma) - Sheth \& Tormen (ST) */
{
  NcMultiplicityFuncST *mst                = NC_MULTIPLICITY_FUNC_ST (mulf);
  NcMultiplicityFuncSTPrivate * const self = mst->priv;

/*  const gdouble bc1 = sqrt(b/(2.0 * M_PI)); */
  const gdouble A   = self->A;
  const gdouble b   = self->b;
  const gdouble bc1 = sqrt (2.0 * b / M_PI);
  const gdouble p   = self->p;
  gdouble x         = self->delta_c / sigma;
  gdouble x2        = x * x;
  /*gdouble b2 = b * b; */
  gdouble f_ST = A * bc1 * (1.0 + pow (x2 * b, -p)) * exp (-(b * x2) / 2.0) * x; /* Jenkin's and Crocce's paper */

  /*gdouble f_ST = A * bc1 * (1.0 + pow(x * b, -p)) * exp(-(b2 * x2) / 2.0) * x; // Evrard' s paper */

  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);
  /*  printf ("A = %.5g, b=%.5g, p=%.5g, delta_c= %.5g\n", A, b, p, mulf_st->delta_c); */

  return f_ST;
}

/**
 * nc_multiplicity_func_st_new:
 *
 * FIXME
 *
 * Returns: A new #NcMultiplicityFuncST.
 */
NcMultiplicityFuncST *
nc_multiplicity_func_st_new (void)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_ST,
                       "mass-def", NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN,
                       NULL);
}

/**
 * nc_multiplicity_func_st_ref:
 * @mst: a #NcMultiplicityFuncST
 *
 * Increases the reference count of @mst by one.
 *
 * Returns: (transfer full): @mst
 */
NcMultiplicityFuncST *
nc_multiplicity_func_st_ref (NcMultiplicityFuncST *mst)
{
  return g_object_ref (mst);
}

/**
 * nc_multiplicity_func_st_free:
 * @mst: a #NcMultiplicityFuncST
 *
 * Atomically decrements the reference count of @mst by one. If the reference count drops to 0,
 * all memory allocated by @mst is released.
 *
 */
void
nc_multiplicity_func_st_free (NcMultiplicityFuncST *mst)
{
  g_object_unref (mst);
}

/**
 * nc_multiplicity_func_st_clear:
 * @mst: a #NcMultiplicityFuncST
 *
 * Atomically decrements the reference count of @mst by one. If the reference count drops to 0,
 * all memory allocated by @mst is released. Set the pointer to NULL;
 *
 */
void
nc_multiplicity_func_st_clear (NcMultiplicityFuncST **mst)
{
  g_clear_object (mst);
}

/**
 * nc_multiplicity_func_st_set_A:
 * @mst: a #NcMultiplicityFuncST.
 * @A: value of #NcMultiplicityFuncST:A.
 *
 * Sets the value @A to the #NcMultiplicityFuncST:A property.
 *
 */
void
nc_multiplicity_func_st_set_A (NcMultiplicityFuncST *mst, gdouble A)
{
  NcMultiplicityFuncSTPrivate * const self = mst->priv;

  g_assert (A >= 0);

  self->A = A;
}

/**
 * nc_multiplicity_func_st_get_A:
 * @mst: a #NcMultiplicityFuncST.
 *
 * Returns: the value of #NcMultiplicityFuncST:A property.
 */
gdouble
nc_multiplicity_func_st_get_A (const NcMultiplicityFuncST *mst)
{
  NcMultiplicityFuncSTPrivate * const self = mst->priv;

  return self->A;
}

/**
 * nc_multiplicity_func_st_set_b:
 * @mst: a #NcMultiplicityFuncST.
 * @b: value of #NcMultiplicityFuncST:b.
 *
 * Sets the value @b to the #NcMultiplicityFuncST:b property.
 *
 */
void
nc_multiplicity_func_st_set_b (NcMultiplicityFuncST *mst, gdouble b)
{
  NcMultiplicityFuncSTPrivate * const self = mst->priv;

  g_assert (b >= 0);

  self->b = b;
}

/**
 * nc_multiplicity_func_st_get_b:
 * @mst: a #NcMultiplicityFuncST.
 *
 * Returns: the value of #NcMultiplicityFuncST:b property.
 */
gdouble
nc_multiplicity_func_st_get_b (const NcMultiplicityFuncST *mst)
{
  NcMultiplicityFuncSTPrivate * const self = mst->priv;

  return self->b;
}

/**
 * nc_multiplicity_func_st_set_p:
 * @mst: a #NcMultiplicityFuncST.
 * @p: value of #NcMultiplicityFuncST:p.
 *
 * Sets the value @p to the #NcMultiplicityFuncST:p property.
 *
 */
void
nc_multiplicity_func_st_set_p (NcMultiplicityFuncST *mst, gdouble p)
{
  NcMultiplicityFuncSTPrivate * const self = mst->priv;

  g_assert (p >= 0);

  self->p = p;
}

/**
 * nc_multiplicity_func_st_get_p:
 * @mst: a #NcMultiplicityFuncST.
 *
 * Returns: the value of #NcMultiplicityFuncST:p property.
 */
gdouble
nc_multiplicity_func_st_get_p (const NcMultiplicityFuncST *mst)
{
  NcMultiplicityFuncSTPrivate * const self = mst->priv;

  return self->p;
}

/**
 * nc_multiplicity_func_st_set_delta_c:
 * @mst: a #NcMultiplicityFuncST.
 * @delta_c: value of #NcMultiplicityFuncST:critical-delta.
 *
 * Sets the value @delta_c to the #NcMultiplicityFuncST:critical-delta property.
 *
 */
void
nc_multiplicity_func_st_set_delta_c (NcMultiplicityFuncST *mst, gdouble delta_c)
{
  NcMultiplicityFuncSTPrivate * const self = mst->priv;

  g_assert (delta_c >= 0);

  self->delta_c = delta_c;
}

/**
 * nc_multiplicity_func_st_get_delta_c:
 * @mst: a #NcMultiplicityFuncST.
 *
 * Returns: the value of #NcMultiplicityFuncST:critical_delta property.
 */
gdouble
nc_multiplicity_func_st_get_delta_c (const NcMultiplicityFuncST *mst)
{
  NcMultiplicityFuncSTPrivate * const self = mst->priv;

  return self->delta_c;
}

