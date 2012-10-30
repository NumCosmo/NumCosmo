/***************************************************************************
 *            nc_multiplicity_func_jenkins.c
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
 * SECTION:nc_multiplicity_func_jenkins
 * @title: Jenkins Multiplicity Function
 * @short_description: Dark Matter Halo FIXME
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func_jenkins.h"

G_DEFINE_TYPE (NcMultiplicityFuncJenkins, nc_multiplicity_func_jenkins, NC_TYPE_MULTIPLICITY_FUNC);

enum
{
  PROP_0,
  PROP_A,
  PROP_A_TCDM,
  PROP_B,
  PROP_B_TCDM,
  PROP_EPSILON,
  PROP_EPSILON_TCDM
};

/**
 * nc_multiplicity_func_jenkins_new:
 * @A: FIXME
 * @A_tCDM: FIXME 
 * @B: FIXME
 * @B_tCDM: FIXME
 * @epsilon: FIXME 
 * @epsilon_tCDM: FIXME
 *   
 * FIXME
 *
 * Returns: A new #NcMultiplicityFunc.
 */
NcMultiplicityFunc *
nc_multiplicity_func_jenkins_new (gdouble A, gdouble A_tCDM, gdouble B, gdouble B_tCDM, gdouble epsilon, gdouble epsilon_tCDM)
{
  return g_object_new (NC_TYPE_MULTIPLICITY_FUNC_JENKINS,
                       "A", A,
                       "A-tCDM", A_tCDM,
                       "B", B, 
                       "B-tCDM", B_tCDM,
                       "epsilon", epsilon,
                       "epsilon-tCDM", epsilon_tCDM,
                       NULL);
}

/* Simulacao FoF: A = 0.315, B = 0.61, epsilon = 3.8. */
/* Simulacao SO(200): A = 0.22, B = 0.73, epsilon = 3.86. Ref. The Astrophys. Journal, 573:7-36, 2002 July 1. */

static gdouble
_nc_multiplicity_func_jenkins_eval (NcMultiplicityFunc *mulf, NcHICosmo *model, gdouble sigma, gdouble z)
{
  NcMultiplicityFuncJenkins *mulf_jenkins = NC_MULTIPLICITY_FUNC_JENKINS (mulf);
  
//  gdouble Omega_m0 = mulf->model->wm / gsl_pow_2(mulf->model->h);
/*
  gdouble x = (1.0 - Omega_m0 * gsl_pow_3(1.0 + z))/(1.0 - Omega_m0);
  gdouble A = (1.0 - x) * mulf_jenkins->A_tCDM + x * mulf_jenkins->A;
  gdouble B = (1.0 - x) * mulf_jenkins->B_tCDM + x * mulf_jenkins->B;
  gdouble epsilon = (1.0 - x) * mulf_jenkins->epsilon_tCDM + x * mulf_jenkins->epsilon;  
  gdouble f_Jenkins = A * exp(-pow(fabs(-log(sigma) + B), epsilon));
*/

/*  Evrard (the incomplete article I have).
  gdouble A = mulf_jenkins->A_tCDM - 0.07 * Omega_m0 * gsl_pow_3(1.0 + z);
  gdouble B = mulf_jenkins->B_tCDM + 0.11 * Omega_m0 * gsl_pow_3(1.0 + z);
  gdouble epsilon = 3.8;
*/

  /* Simulacao FoF */
  mulf_jenkins->A = 0.315;
  mulf_jenkins->B = 0.61;
  mulf_jenkins->epsilon = 3.8;

  gdouble f_Jenkins = mulf_jenkins->A * exp(-pow(fabs(-log(sigma) + mulf_jenkins->B), mulf_jenkins->epsilon));
    
  return f_Jenkins;
}

/**
 * nc_multiplicity_func_jenkins_set_A:
 * @mulf_jenkins: a #NcMultiplicityFuncJenkins.
 * @A: value of #NcMultiplicityFuncJenkins:A.
 *
 * Sets the value @A to the #NcMultiplicityFuncJenkins:A property.
 *
 */
void
nc_multiplicity_func_jenkins_set_A (NcMultiplicityFuncJenkins *mulf_jenkins, gdouble A)
{
  g_assert (A >= 0);
  mulf_jenkins->A = A;
}

/**
 * nc_multiplicity_func_jenkins_get_A:
 * @mulf_jenkins: a #NcMultiplicityFuncJenkins.
 *
 * Returns: the value of #NcMultiplicityFuncJenkins:A property.
 */
gdouble
nc_multiplicity_func_jenkins_get_A (const NcMultiplicityFuncJenkins *mulf_jenkins)
{
  return mulf_jenkins->A;
}

/**
 * nc_multiplicity_func_jenkins_set_A_tCDM:
 * @mulf_jenkins: a #NcMultiplicityFuncJenkins.
 * @A_tCDM: value of #NcMultiplicityFuncJenkins:A-tCDM.
 *
 * Sets the value @A_tCDM to the #NcMultiplicityFuncJenkins:A-tCDM property.
 *
 */
void
nc_multiplicity_func_jenkins_set_A_tCDM (NcMultiplicityFuncJenkins *mulf_jenkins, gdouble A_tCDM)
{
  g_assert (A_tCDM >= 0);
  mulf_jenkins->A_tCDM = A_tCDM;
}

/**
 * nc_multiplicity_func_jenkins_get_A_tCDM:
 * @mulf_jenkins: a #NcMultiplicityFuncJenkins.
 *
 * Returns: the value of #NcMultiplicityFuncJenkins:A-tCDM property.
 */
gdouble
nc_multiplicity_func_jenkins_get_A_tCDM (const NcMultiplicityFuncJenkins *mulf_jenkins)
{
  return mulf_jenkins->A_tCDM;
}

/**
 * nc_multiplicity_func_jenkins_set_B:
 * @mulf_jenkins: a #NcMultiplicityFuncJenkins.
 * @B: value of #NcMultiplicityFuncJenkins:B.
 *
 * Sets the value @B to the #NcMultiplicityFuncJenkins:B property.
 *
 */
void
nc_multiplicity_func_jenkins_set_B (NcMultiplicityFuncJenkins *mulf_jenkins, gdouble B)
{
  g_assert (B >= 0);
  mulf_jenkins->B = B;
}

/**
 * nc_multiplicity_func_jenkins_get_B:
 * @mulf_jenkins: a #NcMultiplicityFuncJenkins.
 *
 * Returns: the value of #NcMultiplicityFuncJenkins:B property.
 */
gdouble
nc_multiplicity_func_jenkins_get_B (const NcMultiplicityFuncJenkins *mulf_jenkins)
{
  return mulf_jenkins->B;
}

/**
 * nc_multiplicity_func_jenkins_set_B_tCDM:
 * @mulf_jenkins: a #NcMultiplicityFuncJenkins.
 * @B_tCDM: value of #NcMultiplicityFuncJenkins:B-tCDM.
 *
 * Sets the value @B_tCDM to the #NcMultiplicityFuncJenkins:B-tCDM property.
 *
 */
void
nc_multiplicity_func_jenkins_set_B_tCDM (NcMultiplicityFuncJenkins *mulf_jenkins, gdouble B_tCDM)
{
  g_assert (B_tCDM >= 0);
  mulf_jenkins->B_tCDM = B_tCDM;
}

/**
 * nc_multiplicity_func_jenkins_get_B_tCDM:
 * @mulf_jenkins: a #NcMultiplicityFuncJenkins.
 *
 * Returns: the value of #NcMultiplicityFuncJenkins:B-tCDM property.
 */
gdouble
nc_multiplicity_func_jenkins_get_B_tCDM (const NcMultiplicityFuncJenkins *mulf_jenkins)
{
  return mulf_jenkins->B_tCDM;
}

/**
 * nc_multiplicity_func_jenkins_set_epsilon:
 * @mulf_jenkins: a #NcMultiplicityFuncJenkins.
 * @epsilon: value of #NcMultiplicityFuncJenkins:epsilon.
 *
 * Sets the value @epsilon to the #NcMultiplicityFuncJenkins:epsilon property.
 *
 */
void
nc_multiplicity_func_jenkins_set_epsilon (NcMultiplicityFuncJenkins *mulf_jenkins, gdouble epsilon)
{
  g_assert (epsilon >= 0);
  mulf_jenkins->epsilon_tCDM = epsilon;
}

/**
 * nc_multiplicity_func_jenkins_get_epsilon:
 * @mulf_jenkins: a #NcMultiplicityFuncJenkins.
 *
 * Returns: the value of #NcMultiplicityFuncJenkins:epsilon property.
 */
gdouble
nc_multiplicity_func_jenkins_get_epsilon (const NcMultiplicityFuncJenkins *mulf_jenkins)
{
  return mulf_jenkins->epsilon;
}

/**
 * nc_multiplicity_func_jenkins_set_epsilon_tCDM:
 * @mulf_jenkins: a #NcMultiplicityFuncJenkins.
 * @epsilon_tCDM: value of #NcMultiplicityFuncJenkins:epsilon-tCDM.
 *
 * Sets the value @epsilon_tCDM to the #NcMultiplicityFuncJenkins:epsilon-tCDM property.
 *
 */
void
nc_multiplicity_func_jenkins_set_epsilon_tCDM (NcMultiplicityFuncJenkins *mulf_jenkins, gdouble epsilon_tCDM)
{
  g_assert (epsilon_tCDM >= 0);
  mulf_jenkins->epsilon_tCDM = epsilon_tCDM;
}

/**
 * nc_multiplicity_func_jenkins_get_epsilon_tCDM:
 * @mulf_jenkins: a #NcMultiplicityFuncJenkins.
 *
 * Returns: the value of #NcMultiplicityFuncJenkins:epsilon-tCDM property.
 */
gdouble
nc_multiplicity_func_jenkins_get_epsilon_tCDM (const NcMultiplicityFuncJenkins *mulf_jenkins)
{
  return mulf_jenkins->epsilon_tCDM;
}

/* The numbers correspond to the astro-ph/number. */
/* I have to verify if the coeficients change with z. FIX */
// _NC_MULTIPLICITY_FUNCTION_JENKINS_DATASET_FOF_0005260 = {0.315, 0.0, 0.61, 0.0, 3.8, 0.0};
// _NC_MULTIPLICITY_FUNCTION_JENKINS_DATASET_SO_0110246 = {0.22, 0.27, 0.73, 0.65, 3.86, 3.77};

static void
nc_multiplicity_func_jenkins_init (NcMultiplicityFuncJenkins *mulf_jenkins)
{
  /* TODO: Add initialization code here */
  mulf_jenkins->A = 0.315;
  mulf_jenkins->A_tCDM = 0.0;
  mulf_jenkins->B = 0.61;
  mulf_jenkins->B_tCDM = 0.0;
  mulf_jenkins->epsilon = 3.8;
  mulf_jenkins->epsilon_tCDM = 0.0;
}

static void
_nc_multiplicity_func_jenkins_finalize (GObject *object)
{
  /* TODO: Add deinitalization code here */

  G_OBJECT_CLASS (nc_multiplicity_func_jenkins_parent_class)->finalize (object);
}

static void
_nc_multiplicity_func_jenkins_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncJenkins *mulf_jenkins = NC_MULTIPLICITY_FUNC_JENKINS (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_JENKINS (object));

  switch (prop_id)
  {
	case PROP_A:
	  mulf_jenkins->A = g_value_get_double (value);
	  break;
	case PROP_A_TCDM:
	  mulf_jenkins->A_tCDM = g_value_get_double (value);
	  break;
	case PROP_B:
	  mulf_jenkins->B = g_value_get_double (value);
	  break;
	case PROP_B_TCDM:
	  mulf_jenkins->B_tCDM = g_value_get_double (value);
	  break;
	case PROP_EPSILON:
	  mulf_jenkins->epsilon = g_value_get_double (value);
	  break;
	case PROP_EPSILON_TCDM:
	  mulf_jenkins->epsilon_tCDM = g_value_get_double (value);
	  break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_multiplicity_func_jenkins_get_property (GObject * object, guint prop_id, GValue * value, GParamSpec * pspec)
{
  NcMultiplicityFuncJenkins *mulf_jenkins = NC_MULTIPLICITY_FUNC_JENKINS (object);
  g_return_if_fail (NC_IS_MULTIPLICITY_FUNC_JENKINS (object));

  switch (prop_id)
  {
	case PROP_A:
	  g_value_set_double (value, mulf_jenkins->A);
	  break;
	case PROP_A_TCDM:
	  g_value_set_double (value, mulf_jenkins->A_tCDM);
	  break;
	case PROP_B:
	  g_value_set_double (value, mulf_jenkins->B);
	  break;
	case PROP_B_TCDM:
	  g_value_set_double (value, mulf_jenkins->B_tCDM);
	  break;
	case PROP_EPSILON:
	  g_value_set_double (value, mulf_jenkins->epsilon);
	  break;
	case PROP_EPSILON_TCDM:
	  g_value_set_double (value, mulf_jenkins->epsilon_tCDM);
	  break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_multiplicity_func_jenkins_class_init (NcMultiplicityFuncJenkinsClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcMultiplicityFuncClass* parent_class = NC_MULTIPLICITY_FUNC_CLASS (klass);
  
  parent_class->eval = &_nc_multiplicity_func_jenkins_eval;
  
  object_class->finalize = _nc_multiplicity_func_jenkins_finalize;
  object_class->set_property = _nc_multiplicity_func_jenkins_set_property;
  object_class->get_property = _nc_multiplicity_func_jenkins_get_property;

  /**
   * NcMultiplicityFuncJenkins:A:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A,
                                   g_param_spec_double ("A",
                                                        NULL,
                                                        "A",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.315,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncJenkins:A-tCDM:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_A_TCDM,
                                   g_param_spec_double ("A-tCDM",
                                                        NULL,
                                                        "A-tCDM",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncJenkins:B:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_B,
                                   g_param_spec_double ("B",
                                                        NULL,
                                                        "B",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.61,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncJenkins:B-tCDM:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_B_TCDM,
                                   g_param_spec_double ("B-tCDM",
                                                        NULL,
                                                        "B-tCDM",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncJenkins:epsilon:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_EPSILON,
                                   g_param_spec_double ("epsilon",
                                                        NULL,
                                                        "epsilon",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 3.8,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /**
   * NcMultiplicityFuncJenkins:epsilon-tCDM:
   *
   * FIXME Set correct values (limits)
   */
  g_object_class_install_property (object_class,
                                   PROP_EPSILON_TCDM,
                                   g_param_spec_double ("epsilon-tCDM",
                                                        NULL,
                                                        "epsilon-tCDM",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

