/***************************************************************************
 *            nc_hipert_itwo_fluids.h
 *
 *  Tue July 22 17:37:11 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_itwo_fluids.h
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

#ifndef _NC_HIPERT_ITWO_FLUIDS_H_
#define _NC_HIPERT_ITWO_FLUIDS_H_

#include <glib-object.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/perturbations/nc_hipert_wkb.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_ITWO_FLUIDS               (nc_hipert_itwo_fluids_get_type ())
#define NC_HIPERT_ITWO_FLUIDS(obj)               (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_ITWO_FLUIDS, NcHIPertITwoFluids))
#define NC_IS_HIPERT_ITWO_FLUIDS(obj)            (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_ITWO_FLUIDS))
#define NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE(obj) (G_TYPE_INSTANCE_GET_INTERFACE ((obj), NC_TYPE_HIPERT_ITWO_FLUIDS, NcHIPertITwoFluidsInterface))

typedef struct _NcHIPertITwoFluids NcHIPertITwoFluids;
typedef struct _NcHIPertITwoFluidsInterface NcHIPertITwoFluidsInterface;
typedef struct _NcHIPertITwoFluidsEOM NcHIPertITwoFluidsEOM;

typedef gdouble (*NcHIPertITwoFluidsFuncNuA2) (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
typedef gdouble (*NcHIPertITwoFluidsFuncNuB2) (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
typedef gdouble (*NcHIPertITwoFluidsFuncVA) (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
typedef gdouble (*NcHIPertITwoFluidsFuncVB) (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
typedef gdouble (*NcHIPertITwoFluidsFuncDlnmzeta) (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
typedef gdouble (*NcHIPertITwoFluidsFuncDlnmS) (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
typedef gdouble (*NcHIPertITwoFluidsFuncDmzetanuAnuA) (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
typedef gdouble (*NcHIPertITwoFluidsFuncDmSnuBnuB) (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
typedef NcHIPertITwoFluidsEOM *(*NcHIPertITwoFluidsFuncEOM) (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);

struct _NcHIPertITwoFluidsInterface
{
  /*< private >*/
  GTypeInterface parent;

  NcHIPertITwoFluidsFuncNuA2 nuA2;
  NcHIPertITwoFluidsFuncNuB2 nuB2;
  NcHIPertITwoFluidsFuncVA VA;
  NcHIPertITwoFluidsFuncVB VB;
  NcHIPertITwoFluidsFuncDlnmzeta dlnmzeta;
  NcHIPertITwoFluidsFuncDlnmS dlnmS;
  NcHIPertITwoFluidsFuncDmzetanuAnuA dmzetanuA_nuA;
  NcHIPertITwoFluidsFuncDmSnuBnuB dmSnuB_nuB;
  NcHIPertITwoFluidsFuncEOM eom;
  NcHIPertITwoFluidsFuncEOM eom_full;
  NcHIPertWKBEom wkb_zeta_eom;
  NcHIPertWKBEom wkb_S_eom;
};

/**
 * NcHICosmoEOMTwoFluids:
 * 
 * FIXME
 * 
 */
struct _NcHIPertITwoFluidsEOM
{
  /*< private >*/
  guint64 skey;
  gdouble alpha;
  gdouble k;
  gdouble mzeta;
  gdouble mS;
  gdouble nuzeta2;
  gdouble nuS2;
  gdouble Y;
  gdouble nu_plus;
  gdouble nu_minus;
  gdouble lambda_zeta;
  gdouble lambda_s;
  gdouble Yt;
  gdouble Uplus;
  gdouble Uminus;
  gdouble Wplus;
  gdouble Wminus;
  gdouble m_plus;
  gdouble m_minus;
  gdouble US[3];
  gdouble UA;
};

GType nc_hipert_itwo_fluids_eom_get_type (void) G_GNUC_CONST;
GType nc_hipert_itwo_fluids_get_type (void) G_GNUC_CONST;

NcHIPertITwoFluidsEOM *nc_hipert_itwo_fluids_eom_dup (NcHIPertITwoFluidsEOM *tf_eom);
void nc_hipert_itwo_fluids_eom_free (NcHIPertITwoFluidsEOM *tf_eom);

G_INLINE_FUNC gdouble nc_hipert_itwo_fluids_nuA2 (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
G_INLINE_FUNC gdouble nc_hipert_itwo_fluids_nuB2 (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);

G_INLINE_FUNC gdouble nc_hipert_itwo_fluids_VA (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
G_INLINE_FUNC gdouble nc_hipert_itwo_fluids_VB (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);

G_INLINE_FUNC gdouble nc_hipert_itwo_fluids_dlnmzeta (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
G_INLINE_FUNC gdouble nc_hipert_itwo_fluids_dlnmS (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);

G_INLINE_FUNC gdouble nc_hipert_itwo_fluids_dmzetanuA_nuA (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
G_INLINE_FUNC gdouble nc_hipert_itwo_fluids_dmSnuB_nuB (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);

G_INLINE_FUNC NcHIPertITwoFluidsEOM *nc_hipert_itwo_fluids_eom (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);
G_INLINE_FUNC NcHIPertITwoFluidsEOM *nc_hipert_itwo_fluids_eom_full (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k);

G_INLINE_FUNC void nc_hipert_itwo_fluids_wkb_zeta_eom (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k, gdouble *nu2, gdouble *m, gdouble *dlnm);
G_INLINE_FUNC void nc_hipert_itwo_fluids_wkb_S_eom (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k, gdouble *nu2, gdouble *m, gdouble *dlnm);

G_END_DECLS

#endif /* _NC_HIPERT_ITWO_FLUIDS_H_ */

#ifndef _NC_HIPERT_ITWO_FLUIDS_INLINE_H_
#define _NC_HIPERT_ITWO_FLUIDS_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC gdouble 
nc_hipert_itwo_fluids_nuA2 (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE (itf)->nuA2 (itf, alpha, k);
}

G_INLINE_FUNC gdouble 
nc_hipert_itwo_fluids_nuB2 (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE (itf)->nuB2 (itf, alpha, k);
}

G_INLINE_FUNC gdouble 
nc_hipert_itwo_fluids_VA (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE (itf)->VA (itf, alpha, k);
}

G_INLINE_FUNC gdouble 
nc_hipert_itwo_fluids_VB (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE (itf)->VB (itf, alpha, k);
}

G_INLINE_FUNC gdouble 
nc_hipert_itwo_fluids_dlnmzeta (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE (itf)->dlnmzeta (itf, alpha, k);
}

G_INLINE_FUNC gdouble 
nc_hipert_itwo_fluids_dlnmS (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE (itf)->dlnmS (itf, alpha, k);
}

G_INLINE_FUNC gdouble 
nc_hipert_itwo_fluids_dmzetanuA_nuA (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE (itf)->dmzetanuA_nuA (itf, alpha, k);
}

G_INLINE_FUNC gdouble 
nc_hipert_itwo_fluids_dmSnuB_nuB (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE (itf)->dmSnuB_nuB (itf, alpha, k);
}

G_INLINE_FUNC NcHIPertITwoFluidsEOM *
nc_hipert_itwo_fluids_eom (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE (itf)->eom (itf, alpha, k);
}

G_INLINE_FUNC NcHIPertITwoFluidsEOM *
nc_hipert_itwo_fluids_eom_full (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE (itf)->eom_full (itf, alpha, k);
}

G_INLINE_FUNC void
nc_hipert_itwo_fluids_wkb_zeta_eom (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k, gdouble *nu2, gdouble *m, gdouble *dlnm)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE (itf)->wkb_zeta_eom (G_OBJECT (itf), alpha, k, nu2, m, dlnm);
}

G_INLINE_FUNC void
nc_hipert_itwo_fluids_wkb_S_eom (NcHIPertITwoFluids *itf, gdouble alpha, gdouble k, gdouble *nu2, gdouble *m, gdouble *dlnm)
{
  return NC_HIPERT_ITWO_FLUIDS_GET_INTERFACE (itf)->wkb_S_eom (G_OBJECT (itf), alpha, k, nu2, m, dlnm);
}


G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPERT_ITWO_FLUIDS_INLINE_H_ */
