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
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_st.h"

G_DEFINE_TYPE (NcMultiplicityFuncST, nc_multiplicity_func_st, NC_TYPE_MULTIPLICITY_FUNC);

enum
{
  PROP_0,
  PROP_A,
  PROP_B,
  PROP_P,
  PROP_DELTA_C
};

/**
 * nc_multiplicity_func_st_new:
 * @A: FIXME
 * @b: FIXME
 * @p: FIXME
 * @delta_c: FIXME
 *
 * FIXME
 *
 * Returns: A new #NcMultiplicityFunc.
 */
NcMultiplicityFunc *
nc_multiplicity_func_st_new (gdouble A, gdouble b, gdouble p, gdouble delta_c)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_ST,
                       "A", A,
                       "b", b,
                       "p", p,
                       "critical-delta", delta_c,
                       NULL);
}

static gdouble
_nc_multiplicity_func_st_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)          /* f(\sigma) - Sheth \& Tormen (ST) */
{
  NcMultiplicityFuncST *mulf_st = NC_MULTIPLICITY_FUNC_ST (mulf);
//  const gdouble bc1 = sqrt(b/(2.0 * M_PI));
  const gdouble A = mulf_st->A;
  const gdouble b = mulf_st->b;
  const gdouble bc1 = sqrt (2.0 * b / M_PI);
  const gdouble p = mulf_st->p;
  gdouble x = mulf_st->delta_c / sigma;
  gdouble x2 = x * x;
  gdouble b2 = b * b;
  gdouble f_ST = A * bc1 * (1.0 + pow(x2 * b, -p)) * exp(-(b * x2) / 2.0) * x; // Jenkin's and Crocce's paper
  //gdouble f_ST = A * bc1 * (1.0 + pow(x * b, -p)) * exp(-(b2 * x2) / 2.0) * x; // Evrard' s paper

  NCM_UNUSED (cosmo);
  NCM_UNUSED (z);
 //  printf ("A = %.5g, b=%.5g, p=%.5g, delta_c= %.5g\n", A, b, p, mulf_st->delta_c);

  return f_ST;
}

/**
 * nc_multiplicity_func_st_set_A:
 * @mulf_st: a #NcMultiplicityFuncST.
 * @A: value of #NcMultiplicityFuncST:A.
 *
 * Sets the value @A to the #NcMultiplicityFuncST:A property.
 *
 */
void
nc_multiplicity_func_st_set_A (NcMultiplicityFuncST *mulf_st, gdouble A)
{
  g_assert (A >= 0);
  mulf_st->A = A;
}

/**
 * nc_multiplicity_func_st_get_A:
 * @mulf_st: a #NcMultiplicityFuncST.
 *
 * Returns: the value of #NcMultiplicityFuncST:A property.
 */
gdouble
nc_multiplicity_func_st_get_A (const NcMultiplicityFuncST *mulf_st)
{
  return mulf_st->A;
}

/**
 * nc_multiplicity_func_st_set_b:
 * @mulf_st: a #NcMultiplicityFuncST.
 * @b: value of #NcMultiplicityFuncST:b.
 *
 * Sets the value @b to the #NcMultiplicityFuncST:b property.
 *
 */
void
nc_multiplicity_func_st_set_b (NcMultiplicityFuncST *mulf_st, gdouble b)
{
  g_assert (b >= 0);
  mulf_st->b = b;
}

/**
 * nc_multiplicity_func_st_get_b:
 * @mulf_st: a #NcMultiplicityFuncST.
 *
 * Returns: the value of #NcMultiplicityFuncST:b property.
 */
gdouble
nc_multiplicity_func_st_get_b (const NcMultiplicityFuncST *mulf_st)
{
  return mulf_st->b;
}

/**
 * nc_multiplicity_func_st_set_p:
 * @mulf_st: a #NcMultiplicityFuncST.
 * @p: value of #NcMultiplicityFuncST:p.
 *
 * Sets the value @p to the #NcMultiplicityFuncST:p property.
 *
 */
void
nc_multiplicity_func_st_set_p (NcMultiplicityFuncST *mulf_st, gdouble p)
{
  g_assert (p >= 0);
  mulf_st->p = p;
}

/**
 * nc_multiplicity_func_st_get_p:
 * @mulf_st: a #NcMultiplicityFuncST.
 *
 * Returns: the value of #NcMultiplicityFuncST:p property.
 */
gdouble
nc_multiplicity_func_st_get_p (const NcMultiplicityFuncST *mulf_st)
{
  return mulf_st->p;
}

/**
 * nc_multiplicity_func_st_set_delta_c:
 * @mulf_st: a #NcMultiplicityFuncST.
 * @delta_c: value of #NcMultiplicityFuncST:critical-delta.
 *
 * Sets the value @delta_c to the #NcMultiplicityFuncST:critical-delta property.
 *
 */
void
nc_multiplicity_func_st_set_delta_c (NcMultiplicityFuncST *mulf_st, gdouble delta_c)
{
  g_assert (delta_c >= 0);
  mulf_st->delta_c = delta_c;
}

/**
 * nc_multiplicity_func_st_get_delta_c:
 * @mulf_st: a #NcMultiplicityFuncST.
 *
 * Returns: the value of #NcMultiplicityFuncST:critical_delta property.
 */
gdouble
nc_multiplicity_func_st_get_delta_c (const NcMultiplicityFuncST *mulf_st)
{
  return mulf_st->delta_c;
}

static void
nc_multiplicity_func_st_init (NcMultiplicityFuncST *mulf_st)
{
  /* TODO: Add initialization code here */
  mulf_st->A = 0.3222;
  mulf_st->b = 0.707;
  mulf_st->p = 0.3;
  mulf_st->delta_c = 1.686;
}

static void
_nc_multiplicity_func_st_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_multiplicity_func_st_parent_class)->finalize (object);
}

static void
_nc_multiplicity_func_st_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncST *mulf_st = NC_MULTIPLICITY_FUNC_ST (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_ST (object));

  switch (prop_id)
  {
	case PROP_A:
	  mulf_st->A = g_value_get_double (value);
	  break;
	case PROP_B:
	  mulf_st->b = g_value_get_double (value);
	  break;
	case PROP_P:
	  mulf_st->p = g_value_get_double (value);
	  break;
	case PROP_DELTA_C:
	  mulf_st->delta_c = g_value_get_double (value);
	  break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_st_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncST *mulf_st = NC_MULTIPLICITY_FUNC_ST (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_ST (object));

  switch (prop_id)
  {
	case PROP_A:
	  g_value_set_double (value, mulf_st->A);
	  break;
	case PROP_B:
	  g_value_set_double (value, mulf_st->b);
	  break;
	case PROP_P:
	  g_value_set_double (value, mulf_st->p);
	  break;
    case PROP_DELTA_C:
      g_value_set_double (value, mulf_st->delta_c);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_multiplicity_func_st_class_init (NcMultiplicityFuncSTClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);

  parent_class->eval = &_nc_multiplicity_func_st_eval;

  object_class->finalize = _nc_multiplicity_func_st_finalize;
  object_class->set_property = _nc_multiplicity_func_st_set_property;
  object_class->get_property = _nc_multiplicity_func_st_get_property;

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
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
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
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.686,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

