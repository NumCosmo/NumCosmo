/***************************************************************************
 *            nc_hipert_boltzmann_cbe.c
 *
 *  Sat October 24 11:56:56 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_boltzmann_cbe.c
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:nc_hipert_boltzmann_cbe
 * @title: NcHIPertBoltzmannCBE
 * @short_description: CLASS (Cosmic Linear Anisotropy Solving System) backend for perturbations
 *
 * If you use this object please cite: [Blas (2011) CLASS II][XBlas2011],
 * see also:
 * - [Lesgourgues (2011) CLASS I][XLesgourgues2011],
 * - [Lesgourgues (2011) CLASS III][XLesgourgues2011a],
 * - [Lesgourgues (2011) CLASS IV][XLesgourgues2011b] and
 * - [CLASS website](http://class-code.net/).
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */

/*
 * It must be include before anything else, several symbols clash
 * with the default includes.
 */
#ifndef NUMCOSMO_GIR_SCAN
#include "class/include/class.h"
#endif /* NUMCOSMO_GIR_SCAN */

#include "build_cfg.h"

#include "nc_hiprim.h"
#include "model/nc_hicosmo_de.h"
#include "nc_hipert_boltzmann_cbe.h"

enum
{
  PROP_0,
  PROP_CBE,
};

G_DEFINE_TYPE (NcHIPertBoltzmannCBE, nc_hipert_boltzmann_cbe, NC_TYPE_HIPERT_BOLTZMANN)

static void
nc_hipert_boltzmann_cbe_init (NcHIPertBoltzmannCBE *boltzmann_cbe)
{
  boltzmann_cbe->cbe        = NULL;
  boltzmann_cbe->PHIPHI_Cls = NULL;
  boltzmann_cbe->TT_Cls     = NULL;
  boltzmann_cbe->EE_Cls     = NULL;
  boltzmann_cbe->BB_Cls     = NULL;
  boltzmann_cbe->TE_Cls     = NULL;
}

static void
nc_hipert_boltzmann_cbe_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPertBoltzmannCBE *boltzmann_cbe = NC_HIPERT_BOLTZMANN_CBE (object);
  g_return_if_fail (NC_IS_HIPERT_BOLTZMANN_CBE (object));

  switch (prop_id)
  {
    case PROP_CBE:
      boltzmann_cbe->cbe = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_boltzmann_cbe_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIPertBoltzmannCBE *boltzmann_cbe = NC_HIPERT_BOLTZMANN_CBE (object);
  g_return_if_fail (NC_IS_HIPERT_BOLTZMANN_CBE (object));

  switch (prop_id)
  {
    case PROP_CBE:
      g_value_set_object (value, boltzmann_cbe->cbe);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_hipert_boltzmann_cbe_dispose (GObject *object)
{
  NcHIPertBoltzmannCBE *boltzmann_cbe = NC_HIPERT_BOLTZMANN_CBE (object);

  nc_cbe_clear (&boltzmann_cbe->cbe);
  ncm_vector_clear (&boltzmann_cbe->PHIPHI_Cls);
  ncm_vector_clear (&boltzmann_cbe->TT_Cls);
  ncm_vector_clear (&boltzmann_cbe->EE_Cls);
  ncm_vector_clear (&boltzmann_cbe->BB_Cls);
  ncm_vector_clear (&boltzmann_cbe->TE_Cls);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_boltzmann_cbe_parent_class)->dispose (object);
}

static void
nc_hipert_boltzmann_cbe_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_boltzmann_cbe_parent_class)->finalize (object);
}

static void _nc_hipert_boltzmann_cbe_prepare (NcHIPertBoltzmann *pb, NcHICosmo *cosmo);
static void _nc_hipert_boltzmann_cbe_get_PHIPHI_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
static void _nc_hipert_boltzmann_cbe_get_TT_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
static void _nc_hipert_boltzmann_cbe_get_EE_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
static void _nc_hipert_boltzmann_cbe_get_BB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
static void _nc_hipert_boltzmann_cbe_get_TE_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
static void _nc_hipert_boltzmann_cbe_get_TB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);
static void _nc_hipert_boltzmann_cbe_get_EB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls);

static void
nc_hipert_boltzmann_cbe_class_init (NcHIPertBoltzmannCBEClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);


  object_class->set_property = nc_hipert_boltzmann_cbe_set_property;
  object_class->get_property = nc_hipert_boltzmann_cbe_get_property;
  object_class->dispose      = nc_hipert_boltzmann_cbe_dispose;
  object_class->finalize     = nc_hipert_boltzmann_cbe_finalize;

  g_object_class_install_property (object_class,
                                   PROP_CBE,
                                   g_param_spec_object ("cbe",
                                                        NULL,
                                                        "CLASS backend object",
                                                        NC_TYPE_CBE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  NC_HIPERT_BOLTZMANN_CLASS (klass)->prepare        = &_nc_hipert_boltzmann_cbe_prepare;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_PHIPHI_Cls = &_nc_hipert_boltzmann_cbe_get_PHIPHI_Cls;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_TT_Cls     = &_nc_hipert_boltzmann_cbe_get_TT_Cls;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_EE_Cls     = &_nc_hipert_boltzmann_cbe_get_EE_Cls;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_BB_Cls     = &_nc_hipert_boltzmann_cbe_get_BB_Cls;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_TE_Cls     = &_nc_hipert_boltzmann_cbe_get_TE_Cls;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_TB_Cls     = &_nc_hipert_boltzmann_cbe_get_TB_Cls;
  NC_HIPERT_BOLTZMANN_CLASS (klass)->get_EB_Cls     = &_nc_hipert_boltzmann_cbe_get_EB_Cls;
}

static void
_nc_hipert_boltzmann_cbe_get_PHIPHI_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  NcHIPertBoltzmannCBE *boltzmann_cbe = NC_HIPERT_BOLTZMANN_CBE (pb);
  ncm_vector_memcpy2 (Cls, boltzmann_cbe->PHIPHI_Cls, 0, 0, ncm_vector_len (Cls));
}

static void
_nc_hipert_boltzmann_cbe_get_TT_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  NcHIPertBoltzmannCBE *boltzmann_cbe = NC_HIPERT_BOLTZMANN_CBE (pb);
  ncm_vector_memcpy2 (Cls, boltzmann_cbe->TT_Cls, 0, 0, ncm_vector_len (Cls));
}

static void
_nc_hipert_boltzmann_cbe_get_EE_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  NcHIPertBoltzmannCBE *boltzmann_cbe = NC_HIPERT_BOLTZMANN_CBE (pb);
  ncm_vector_memcpy2 (Cls, boltzmann_cbe->EE_Cls, 0, 0, ncm_vector_len (Cls));
}

static void
_nc_hipert_boltzmann_cbe_get_BB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  NcHIPertBoltzmannCBE *boltzmann_cbe = NC_HIPERT_BOLTZMANN_CBE (pb);
  ncm_vector_memcpy2 (Cls, boltzmann_cbe->BB_Cls, 0, 0, ncm_vector_len (Cls));
}

static void
_nc_hipert_boltzmann_cbe_get_TE_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  NcHIPertBoltzmannCBE *boltzmann_cbe = NC_HIPERT_BOLTZMANN_CBE (pb);
  ncm_vector_memcpy2 (Cls, boltzmann_cbe->TE_Cls, 0, 0, ncm_vector_len (Cls));
}

static void
_nc_hipert_boltzmann_cbe_get_TB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  NcHIPertBoltzmannCBE *boltzmann_cbe = NC_HIPERT_BOLTZMANN_CBE (pb);
  ncm_vector_memcpy2 (Cls, boltzmann_cbe->TB_Cls, 0, 0, ncm_vector_len (Cls));
}

static void
_nc_hipert_boltzmann_cbe_get_EB_Cls (NcHIPertBoltzmann *pb, NcmVector *Cls)
{
  NcHIPertBoltzmannCBE *boltzmann_cbe = NC_HIPERT_BOLTZMANN_CBE (pb);
  ncm_vector_memcpy2 (Cls, boltzmann_cbe->EB_Cls, 0, 0, ncm_vector_len (Cls));
}

/**
 * nc_hipert_boltzmann_cbe_new: (constructor)
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcCBEPrecision.
 */
NcHIPertBoltzmannCBE *
nc_hipert_boltzmann_cbe_new (void)
{
  NcCBE *cbe = nc_cbe_new ();
  NcHIPertBoltzmannCBE *boltzmann_cbe = g_object_new (NC_TYPE_HIPERT_BOLTZMANN_CBE,
                                                      "cbe", cbe,
                                                      NULL);
  return boltzmann_cbe;
}

/**
 * nc_hipert_boltzmann_cbe_full_new: (constructor)
 * @cbe: a #NcCBE.
 *
 * FIXME
 *
 * Returns: (transfer full): a new #NcHIPertBoltzmannCBE.
 */
NcHIPertBoltzmannCBE *
nc_hipert_boltzmann_cbe_full_new (NcCBE *cbe)
{
  NcHIPertBoltzmannCBE *boltzmann_cbe = g_object_new (NC_TYPE_HIPERT_BOLTZMANN_CBE,
                                                      "cbe", cbe,
                                                      NULL);
  return boltzmann_cbe;
}

/**
 * nc_hipert_boltzmann_cbe_ref:
 * @boltzmann_cbe: a #NcHIPertBoltzmannCBE.
 *
 * Increases the reference count of @boltzmann_cbe.
 *
 * Returns: (transfer full): @boltzmann_cbe.
 */
NcHIPertBoltzmannCBE *
nc_hipert_boltzmann_cbe_ref (NcHIPertBoltzmannCBE *boltzmann_cbe)
{
  return g_object_ref (boltzmann_cbe);
}

/**
 * nc_hipert_boltzmann_cbe_free:
 * @boltzmann_cbe: a #NcHIPertBoltzmannCBE.
 *
 * Decreases the reference count of @boltzmann_cbe.
 *
 */
void
nc_hipert_boltzmann_cbe_free (NcHIPertBoltzmannCBE *boltzmann_cbe)
{
  g_object_unref (boltzmann_cbe);
}

/**
 * nc_hipert_boltzmann_cbe_clear:
 * @boltzmann_cbe: a #NcHIPertBoltzmannCBE.
 *
 * Decreases the reference count of *@boltzmann_cbe and sets *@boltzmann_cbe to NULL.
 *
 */
void
nc_hipert_boltzmann_cbe_clear (NcHIPertBoltzmannCBE **boltzmann_cbe)
{
  g_clear_object (boltzmann_cbe);
}

static void
_nc_hipert_boltzmann_cbe_prepare (NcHIPertBoltzmann *pb, NcHICosmo *cosmo)
{
  NcHIPertBoltzmannCBE *boltzmann_cbe = NC_HIPERT_BOLTZMANN_CBE (pb);

  guint PHIPHI_lmax = nc_hipert_boltzmann_get_PHIPHI_lmax (pb);
  guint TT_lmax     = nc_hipert_boltzmann_get_TT_lmax (pb);
  guint EE_lmax     = nc_hipert_boltzmann_get_EE_lmax (pb);
  guint BB_lmax     = nc_hipert_boltzmann_get_BB_lmax (pb);
  guint TE_lmax     = nc_hipert_boltzmann_get_TE_lmax (pb);
  guint TB_lmax     = nc_hipert_boltzmann_get_TB_lmax (pb);
  guint EB_lmax     = nc_hipert_boltzmann_get_EB_lmax (pb);
  guint scalar_lmax = 0;
  
#define _CHECK_VEC(name) \
G_STMT_START { \
 if (boltzmann_cbe->name##_Cls != NULL) \
 { \
   if (ncm_vector_len (boltzmann_cbe->name##_Cls) != name##_lmax + 1) \
   { \
     ncm_vector_clear (&boltzmann_cbe->name##_Cls); \
     boltzmann_cbe->name##_Cls = ncm_vector_new (name##_lmax + 1); \
   } \
 } \
 else \
   boltzmann_cbe->name##_Cls = ncm_vector_new (name##_lmax + 1); \
} G_STMT_END

	_CHECK_VEC (PHIPHI);
  _CHECK_VEC (TT);
  _CHECK_VEC (EE);
  _CHECK_VEC (BB);
  _CHECK_VEC (TE);
  _CHECK_VEC (TB);
  _CHECK_VEC (EB);

  scalar_lmax = 0;
  scalar_lmax = GSL_MAX (scalar_lmax, PHIPHI_lmax);
  scalar_lmax = GSL_MAX (scalar_lmax, TT_lmax);
  scalar_lmax = GSL_MAX (scalar_lmax, EE_lmax);
  scalar_lmax = GSL_MAX (scalar_lmax, BB_lmax);
  scalar_lmax = GSL_MAX (scalar_lmax, TE_lmax);
  scalar_lmax = GSL_MAX (scalar_lmax, TB_lmax);
  scalar_lmax = GSL_MAX (scalar_lmax, EB_lmax);

  nc_cbe_set_target_Cls (boltzmann_cbe->cbe, nc_hipert_boltzmann_get_target_Cls (pb));
  nc_cbe_set_lensed_Cls (boltzmann_cbe->cbe, nc_hipert_boltzmann_lensed_Cls (pb));
  nc_cbe_set_tensor (boltzmann_cbe->cbe, nc_hipert_boltzmann_tensor (pb));

  nc_cbe_set_scalar_lmax (boltzmann_cbe->cbe, scalar_lmax);

  nc_cbe_prepare_if_needed (boltzmann_cbe->cbe, cosmo);
  nc_cbe_get_all_Cls (boltzmann_cbe->cbe,
                      boltzmann_cbe->PHIPHI_Cls,
                      boltzmann_cbe->TT_Cls, 
                      boltzmann_cbe->EE_Cls, 
                      boltzmann_cbe->BB_Cls, 
                      boltzmann_cbe->TE_Cls);  
}

