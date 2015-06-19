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
  PROP_LMAX,
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
  pb->lmax                   = 0;
  pb->tight_coupling         = FALSE;
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
    case PROP_LMAX:
      nc_hipert_boltzmann_set_lmax (pb, g_value_get_uint (value));
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
    case PROP_LMAX:
      g_value_set_uint (value, pb->lmax);
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
  nc_scale_factor_clear (&pb->a);
  
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
                                   PROP_LMAX,
                                   g_param_spec_uint ("l-max",
                                                      NULL,
                                                      "Last multipole",
                                                      2, G_MAXUINT32, 2,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  
  NC_HIPERT_CLASS (klass)->set_mode_k = &_nc_hipert_boltzmann_set_mode_k;
  NC_HIPERT_CLASS (klass)->set_abstol = &_nc_hipert_boltzmann_set_abstol;
  NC_HIPERT_CLASS (klass)->set_reltol = &_nc_hipert_boltzmann_set_reltol;  
}

static void 
_nc_hipert_boltzmann_set_mode_k (NcHIPert *pert, gdouble k) 
{
  NC_HIPERT_CLASS (nc_hipert_boltzmann_parent_class)->set_mode_k (pert, k);
  /* Chain up : start */
  if (!pert->prepared)
  {
    NcHIPertBoltzmann *pb = NC_HIPERT_BOLTZMANN (pert);
    NCM_UNUSED (pb);
  }
}

static void 
_nc_hipert_boltzmann_set_abstol (NcHIPert *pert, gdouble abstol) 
{
  NC_HIPERT_CLASS (nc_hipert_boltzmann_parent_class)->set_abstol (pert, abstol);
  /* Chain up : start */
  if (!pert->prepared)
  {
    NcHIPertBoltzmann *pb = NC_HIPERT_BOLTZMANN (pert);
    NCM_UNUSED (pb);
  }
}

static void 
_nc_hipert_boltzmann_set_reltol (NcHIPert *pert, gdouble reltol) 
{
  NC_HIPERT_CLASS (nc_hipert_boltzmann_parent_class)->set_reltol (pert, reltol);
  /* Chain up : start */
  if (!pert->prepared)
  {
    NcHIPertBoltzmann *pb = NC_HIPERT_BOLTZMANN (pert);
    NCM_UNUSED (pb);
  }
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
 * nc_hipert_boltzmann_set_lmax:
 * @pb: a #NcHIPertBoltzmann.
 * @lmax: last mutipole.
 * 
 * FIXME
 * 
 */
void
nc_hipert_boltzmann_set_lmax (NcHIPertBoltzmann *pb, guint lmax)
{
  g_assert_cmpuint (lmax, >=, 2);

  pb->lmax = lmax;
  nc_hipert_set_sys_size (NC_HIPERT (pb), NC_HIPERT_BOLTZMANN_LEN + 2 * (lmax + 1 - 3));
}

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
  }
}

void 
nc_hipert_boltzmann_prepare (NcHIPertBoltzmann *pb, NcHICosmo *cosmo)
{
  /*NcHIPert *pert = NC_HIPERT (pb);*/

  nc_recomb_prepare_if_needed (pb->recomb, cosmo);
  
  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->init (pb, cosmo);
  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->reset (pb);
  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->set_opts (pb);
  NC_HIPERT_BOLTZMANN_GET_CLASS (pb)->evol (pb, pb->lambdaf);
}


