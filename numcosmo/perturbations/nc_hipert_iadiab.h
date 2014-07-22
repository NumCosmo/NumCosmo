/***************************************************************************
 *            nc_hipert_iadiab.h
 *
 *  Fri July 18 15:10:14 2014
 *  Copyright  2014  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_iadiab.h
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

#ifndef _NC_HIPERT_IADIAB_H_
#define _NC_HIPERT_IADIAB_H_

#include <glib-object.h>
#include <numcosmo/nc_hicosmo.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_IADIAB               (nc_hipert_iadiab_get_type ())
#define NC_HIPERT_IADIAB(obj)               (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_IADIAB, NcHIPertIAdiab))
#define NC_IS_HIPERT_IADIAB(obj)            (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_IADIAB))
#define NC_HIPERT_IADIAB_GET_INTERFACE(obj) (G_TYPE_INSTANCE_GET_INTERFACE ((obj), NC_TYPE_HIPERT_IADIAB, NcHIPertIAdiabInterface))

typedef struct _NcHIPertIAdiab NcHIPertIAdiab;
typedef struct _NcHIPertIAdiabInterface NcHIPertIAdiabInterface;
typedef struct _NcHIPertIAdiabEOM NcHIPertIAdiabEOM;

typedef gdouble (*NcHIPertIAdiabFuncNuA2) (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k);
typedef gdouble (*NcHIPertIAdiabFuncDmzetanuAnuA) (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k);
typedef NcHIPertIAdiabEOM *(*NcHIPertIAdiabFuncEOM) (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k);

struct _NcHIPertIAdiabInterface
{
  /*< private >*/
  GTypeInterface parent;

  NcHIPertIAdiabFuncNuA2 nuA2;
  NcHIPertIAdiabFuncDmzetanuAnuA dmzetanuA_nuA;
  NcHIPertIAdiabFuncEOM eom;
};

/**
 * NcHIPertIAdiabEOM:
 * 
 * FIXME
 * 
 */
struct _NcHIPertIAdiabEOM
{
  /*< private >*/
  guint64 skey;
  gdouble alpha;
  gdouble k;
  gdouble m;
  gdouble nu2;
};

GType nc_hipert_iadiab_eom_get_type (void) G_GNUC_CONST;
GType nc_hipert_iadiab_get_type (void) G_GNUC_CONST;

NcHIPertIAdiabEOM *nc_hipert_iadiab_eom_dup (NcHIPertIAdiabEOM *adiab_eom);
void nc_hipert_iadiab_eom_free (NcHIPertIAdiabEOM *adiab_eom);

G_INLINE_FUNC gdouble nc_hipert_iadiab_nuA2 (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k);
G_INLINE_FUNC gdouble nc_hipert_iadiab_dmzetanuA_nuA (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k);
G_INLINE_FUNC NcHIPertIAdiabEOM *nc_hipert_iadiab_eom (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k);

G_END_DECLS

#endif /* _NC_HIPERT_IADIAB_H_ */

#ifndef _NC_HIPERT_IADIAB_INLINE_H_
#define _NC_HIPERT_IADIAB_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC gdouble 
nc_hipert_iadiab_nuA2 (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k)
{
  return NC_HIPERT_IADIAB_GET_INTERFACE (iadiab)->nuA2 (iadiab, alpha, k);
}

G_INLINE_FUNC gdouble 
nc_hipert_iadiab_dmzetanuA_nuA (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k)
{
  return NC_HIPERT_IADIAB_GET_INTERFACE (iadiab)->dmzetanuA_nuA (iadiab, alpha, k);
}

G_INLINE_FUNC NcHIPertIAdiabEOM *
nc_hipert_iadiab_eom (NcHIPertIAdiab *iadiab, gdouble alpha, gdouble k)
{
  return NC_HIPERT_IADIAB_GET_INTERFACE (iadiab)->eom (iadiab, alpha, k);
}



G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_HIPERT_IADIAB_INLINE_H_ */
