/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
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
 * @title: NcHIPert
 * @short_description: Abstract class for perturbation in homogeneous and isotropic cosmologies.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include <stdio.h>
#include "perturbations/nc_hipert.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

#define SUN_DENSE_ACCESS SM_ELEMENT_D
#define SUN_BAND_ACCESS SM_ELEMENT_D
#endif /* NUMCOSMO_GIR_SCAN */

#include "perturbations/nc_hipert_private.h"

enum {
  PROP_0,
  PROP_K,
  PROP_SYS_SIZE,
  PROP_RELTOL,
  PROP_ABSTOL,
  PROP_ALPHAI,
  PROP_SIZE,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcHIPert, nc_hipert, G_TYPE_OBJECT);

static void
nc_hipert_init (NcHIPert *pert)
{
  NcHIPertPrivate * const self = pert->priv = nc_hipert_get_instance_private (pert);
  self->alpha0      = 0.0;
  self->reltol      = 0.0;
  self->k           = 0.0;
  self->sys_size    = 0;
  self->y           = NULL;
  self->A           = NULL;
  self->LS          = NULL;

  self->vec_abstol  = NULL;
  self->cvode       = CVodeCreate (CV_ADAMS);
  self->cvode_init  = FALSE;
  self->cvode_stiff = FALSE;
}

static void
nc_hipert_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPert *pert = NC_HIPERT (object);
  NcHIPertPrivate * const self = pert->priv;
  g_return_if_fail (NC_IS_HIPERT (object));

  switch (prop_id)
  {
    case PROP_K:
      nc_hipert_set_mode_k (pert, g_value_get_double (value));
      break;
    case PROP_SYS_SIZE:
      nc_hipert_set_sys_size (pert, g_value_get_uint (value));
      break;
    case PROP_RELTOL:
      nc_hipert_set_reltol (pert, g_value_get_double (value));
      break;
    case PROP_ABSTOL:
      nc_hipert_set_abstol (pert, g_value_get_double (value));
      break;
    case PROP_ALPHAI:
      self->alpha0 = g_value_get_double (value);
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
  NcHIPertPrivate * const self = pert->priv;
  g_return_if_fail (NC_IS_HIPERT (object));

  switch (prop_id)
  {
    case PROP_K:
      g_value_set_double (value, self->k);
      break;
    case PROP_SYS_SIZE:
      g_value_set_uint (value, self->sys_size);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, nc_hipert_get_reltol (pert));
      break;
    case PROP_ABSTOL:
      g_value_set_double (value, nc_hipert_get_abstol (pert));
      break;
    case PROP_ALPHAI:
      g_value_set_double (value, self->alpha0);
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
  NcHIPertPrivate * const self = pert->priv;

  if (self->cvode != NULL)
  {
    CVodeFree (&self->cvode);
    self->cvode = NULL;
  }

  if (self->y != NULL)
  {
    N_VDestroy (self->y);
    self->y = NULL;
  }

  if (self->vec_abstol != NULL)
  {
    N_VDestroy (self->vec_abstol);
    self->vec_abstol = NULL;
  }
  
  if (self->A != NULL)
  {
    SUNMatDestroy (self->A);
    self->A = NULL;
  }
  
  if (self->LS != NULL)
  {
    SUNLinSolFree (self->LS);
    self->LS = NULL;
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
                                   PROP_SYS_SIZE,
                                   g_param_spec_uint ("sys-size",
                                                      NULL,
                                                      "System size",
                                                      1, G_MAXUINT32, 1,
                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_RELTOL,
                                   g_param_spec_double ("reltol",
                                                        NULL,
                                                        "Relative tolerance",
                                                        0.0, 1.0, NCM_DEFAULT_PRECISION,
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
  NcHIPertPrivate * const self = pert->priv;
  if (self->k != k)
  {
    self->k        = k;
    self->prepared = FALSE;
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
  NcHIPertPrivate * const self = pert->priv;
  if (self->reltol != reltol)
  {
    self->reltol   = reltol;
    self->prepared = FALSE;
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
  NcHIPertPrivate * const self = pert->priv;
  if (self->abstol != abstol)
  {
    self->abstol   = abstol;
    self->prepared = FALSE;
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
  NcHIPertPrivate * const self = pert->priv;
  return self->reltol;
}

/**
 * nc_hipert_get_abstol:
 * @pert: a #NcHIPert.
 *
 * Gets the value of the relative tolerance.
 * 
 * Returns: abstol
 */
gdouble 
nc_hipert_get_abstol (NcHIPert *pert)
{
  NcHIPertPrivate * const self = pert->priv;
  return self->abstol;
}

/**
 * nc_hipert_get_mode_k:
 * @pert: a #NcHIPert.
 *
 * Gets the value of the mode $k$.
 * 
 * Returns: $k$
 */
gdouble 
nc_hipert_get_mode_k (NcHIPert *pert)
{
  NcHIPertPrivate * const self = pert->priv;
  return self->k;
}

/**
 * nc_hipert_prepared:
 * @pert: a #NcHIPert.
 *
 * Returns: Whether the object is prepared.
 */
gboolean 
nc_hipert_prepared (NcHIPert *pert)
{
  NcHIPertPrivate * const self = pert->priv;
  return self->prepared;
}

/**
 * nc_hipert_set_prepared:
 * @pert: a #NcHIPert
 * @prepared: a boolean
 * 
 * Sets the object to @prepared.
 * 
 */
void 
nc_hipert_set_prepared (NcHIPert *pert, gboolean prepared)
{
  NcHIPertPrivate * const self = pert->priv;
  self->prepared = prepared;
}

/**
 * nc_hipert_set_sys_size:
 * @pert: a #NcHIPert.
 * @sys_size: the system size.
 *
 * Sets the system size.
 *
 */
void
nc_hipert_set_sys_size (NcHIPert *pert, guint sys_size)
{
  NcHIPertPrivate * const self = pert->priv;
  if (self->sys_size != sys_size)
  {
    if (self->y != NULL)
    {
      N_VDestroy (self->y);
      self->y = NULL;
    }
    if (self->A != NULL)
    {
      SUNMatDestroy (self->A);
      self->A = NULL;
    }

    if (self->LS != NULL)
    {
      SUNLinSolFree (self->LS);
      self->LS = NULL;
    }
    
    if (self->vec_abstol != NULL)
    {
      N_VDestroy (self->vec_abstol);
      self->vec_abstol = NULL;
    }

    self->sys_size = sys_size;
    if (self->sys_size > 0)
    {
      self->y          = N_VNew_Serial (sys_size);
      self->vec_abstol = N_VNew_Serial (sys_size);

      self->A          = SUNDenseMatrix (sys_size, sys_size);
      self->LS         = SUNDenseLinearSolver (self->y, self->A);

      NCM_CVODE_CHECK ((gpointer)self->A, "SUNDenseMatrix", 0, );
      NCM_CVODE_CHECK ((gpointer)self->LS, "SUNDenseLinearSolver", 0, );
    }
    self->prepared = FALSE;

    nc_hipert_reset_solver (pert);
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
  NcHIPertPrivate * const self = pert->priv;
  guint a = stiff ? 1 : 0;
  guint b = self->cvode_stiff ? 1 : 0;
  if (a != b)
  {
    if (self->cvode != NULL)
    {
      CVodeFree (&self->cvode);
      self->cvode = NULL;
    }
    self->cvode_stiff = stiff;

    if (stiff)
      self->cvode = CVodeCreate (CV_BDF);
    else
      self->cvode = CVodeCreate (CV_ADAMS);

    self->cvode_init = FALSE;
  }
}

/**
 * nc_hipert_reset_solver:
 * @pert: a #NcHIPert
 *
 * Destroy and recreates the solver.
 *
 */
void
nc_hipert_reset_solver (NcHIPert *pert)
{
  NcHIPertPrivate * const self = pert->priv;
  if (self->cvode != NULL)
  {
    CVodeFree (&self->cvode);
    self->cvode = NULL;
  }

  if (self->cvode_stiff)
    self->cvode = CVodeCreate (CV_BDF);
  else
    self->cvode = CVodeCreate (CV_ADAMS);

  self->cvode_init = FALSE;
}

