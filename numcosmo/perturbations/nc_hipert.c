/***************************************************************************
 *            nc_hipert.c
 *
 *  Tue June 03 15:47:45 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert.c
 * Copyright (C) 2014 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_hipert
 * @title: General Perturbation Object
 * @short_description: Perturbation object for homogeneous and isotropic cosmologies 
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/nc_hipert.h"

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <nvector/nvector_serial.h> 

enum {
  PROP_0,
  PROP_K,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_ALPHAI,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE (NcHIPert, nc_hipert, G_TYPE_OBJECT);

static void
nc_hipert_init (NcHIPert *pert)
{
  pert->alpha0      = 0.0;
  pert->reltol      = 0.0;
  pert->k           = 0.0;
  pert->y           = NULL;
  pert->cvode       = CVodeCreate (CV_ADAMS, CV_FUNCTIONAL);
  pert->cvode_init  = FALSE;
  pert->cvode_stiff = FALSE;
}

static void
nc_hipert_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPert *pert = NC_HIPERT (object);
  g_return_if_fail (NC_IS_HIPERT (object));

  switch (prop_id)
  {
    case PROP_K:
      nc_hipert_set_mode_k (pert, g_value_get_double (value));
      break;
    case PROP_RELTOL:
      nc_hipert_set_reltol (pert, g_value_get_double (value));
      break;    
    case PROP_ABSTOL:
      nc_hipert_set_abstol (pert, g_value_get_double (value));
      break;    
    case PROP_ALPHAI:
      pert->alpha0 = g_value_get_double (value);
      break;    
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIPert *pert = NC_HIPERT (object);
  g_return_if_fail (NC_IS_HIPERT (object));

  switch (prop_id)
  {
    case PROP_K:
      g_value_set_double (value, pert->k);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, nc_hipert_get_reltol (pert));
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, nc_hipert_get_abstol (pert));
      break;
    case PROP_ALPHAI:
      g_value_set_double (value, pert->alpha0);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_finalize (GObject *object)
{
  NcHIPert *pert = NC_HIPERT (object);

  if (pert->cvode != NULL)
  {
    CVodeFree (&pert->cvode);
    pert->cvode = NULL;
  }
  
  if (pert->y != NULL)
  {
    N_VDestroy (pert->y);
    pert->y = NULL;
  }

  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_parent_class)->finalize (object);
}

static void _nc_hipert_set_mode_k (NcHIPert *pert, gdouble k);
static void _nc_hipert_set_abstol (NcHIPert *pert, gdouble abstol);
static void _nc_hipert_set_reltol (NcHIPert *pert, gdouble reltol);

static void
nc_hipert_class_init (NcHIPertClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &nc_hipert_set_property;
  object_class->get_property = &nc_hipert_get_property;
  object_class->finalize     = &nc_hipert_finalize;

  g_object_class_install_property (object_class,
                                   PROP_K,
                                   g_param_spec_double ("k",
                                                        NULL,
                                                        "Mode k",
                                                        0.0, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance",
                                                        0.0, 1.0, 1e-7,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_double ("abstol",
                                                        NULL,
                                                        "Absolute tolerance",
                                                        0.0, 1.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_ALPHAI,
                                   g_param_spec_double ("alphai",
                                                        NULL,
                                                        "Initial time",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  klass->set_mode_k = &_nc_hipert_set_mode_k;
  klass->set_abstol = &_nc_hipert_set_abstol;
  klass->set_reltol = &_nc_hipert_set_reltol;
}

/**
 * nc_hipert_set_mode_k:
 * @pert: a #NcHIPert.
 * @k: the mode value.
 * 
 * Sets the value of the mode to be computed.
 * 
 */
static void 
_nc_hipert_set_mode_k (NcHIPert *pert, gdouble k)
{
  if (pert->k != k)
  {
    pert->k        = k;
    pert->prepared = FALSE;
  }
}

/**
 * nc_hipert_set_stiff_solver:
 * @pert: a #NcHIPert.
 * @stiff: whenever to enable or disable a stiff solver.
 * 
 * Sets the ode algorithm to use.
 * 
 */
void 
nc_hipert_set_stiff_solver (NcHIPert *pert, gboolean stiff)
{
  guint a = stiff ? 1 : 0;
  guint b = pert->cvode_stiff ? 1 : 0;
  if (a != b)
  {
    if (pert->cvode != NULL)
    {
      CVodeFree (&pert->cvode);
      pert->cvode = NULL;
    }
    pert->cvode_stiff = stiff;

    if (stiff)
      pert->cvode = CVodeCreate (CV_BDF, CV_NEWTON);
    else
      pert->cvode = CVodeCreate (CV_ADAMS, CV_FUNCTIONAL);
    
    pert->cvode_init = FALSE;
  }
}


/**
 * nc_hipert_set_reltol:
 * @pert: a #NcHIPert.
 * @reltol: the relative tolarance.
 * 
 * Sets the value of the relative tolerance.
 * 
 */
static void 
_nc_hipert_set_reltol (NcHIPert *pert, gdouble reltol)
{
  if (pert->reltol != reltol)
  {
    pert->reltol   = reltol;
    pert->prepared = FALSE;
  }
}

/**
 * nc_hipert_set_abstol:
 * @pert: a #NcHIPert.
 * @abstol: the absolute tolarance.
 * 
 * Sets the value of the absolute tolerance.
 * 
 */
static void 
_nc_hipert_set_abstol (NcHIPert *pert, gdouble abstol)
{
  if (pert->abstol != abstol)
  {
    pert->abstol   = abstol;
    pert->prepared = FALSE;
  }
}

/**
 * nc_hipert_get_reltol:
 * @pert: a #NcHIPert.
 * 
 * Gets the value of the relative tolerance.
 * 
 * Returns: reltol
 */
gdouble 
nc_hipert_get_reltol (NcHIPert *pert)
{
  return pert->reltol;
}

/**
 * nc_hipert_get_abstol:
 * @pert: a #NcHIPert.
 * 
 * Gets the value of the relative tolerance.
 * 
 */
gdouble
nc_hipert_get_abstol (NcHIPert *pert)
{
  return pert->abstol;
}
