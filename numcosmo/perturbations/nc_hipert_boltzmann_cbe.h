/***************************************************************************
 *            nc_hipert_boltzmann_cbe.h
 *
 *  Sat October 24 11:57:37 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_hipert_boltzmann_cbe.h
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

#ifndef _NC_HIPERT_BOLTZMANN_CBE_H_
#define _NC_HIPERT_BOLTZMANN_CBE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/nc_cbe.h>
#include <numcosmo/perturbations/nc_hipert_boltzmann.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_BOLTZMANN_CBE             (nc_hipert_boltzmann_cbe_get_type ())
#define NC_HIPERT_BOLTZMANN_CBE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_BOLTZMANN_CBE, NcHIPertBoltzmannCBE))
#define NC_HIPERT_BOLTZMANN_CBE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_BOLTZMANN_CBE, NcHIPertBoltzmannCBEClass))
#define NC_IS_HIPERT_BOLTZMANN_CBE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_BOLTZMANN_CBE))
#define NC_IS_HIPERT_BOLTZMANN_CBE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_BOLTZMANN_CBE))
#define NC_HIPERT_BOLTZMANN_CBE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_BOLTZMANN_CBE, NcHIPertBoltzmannCBEClass))

typedef struct _NcHIPertBoltzmannCBEClass NcHIPertBoltzmannCBEClass;
typedef struct _NcHIPertBoltzmannCBE NcHIPertBoltzmannCBE;

struct _NcHIPertBoltzmannCBEClass
{
  /*< private >*/
  NcHIPertBoltzmannClass parent_class;
};

struct _NcHIPertBoltzmannCBE
{
  /*< private >*/
  NcHIPertBoltzmann parent_instance;
  NcCBE *cbe;
  NcmVector *PHIPHI_Cls;
  NcmVector *TT_Cls;
  NcmVector *EE_Cls;
  NcmVector *BB_Cls;
  NcmVector *TE_Cls;
  NcmVector *TB_Cls;
  NcmVector *EB_Cls;
};

GType nc_hipert_boltzmann_cbe_get_type (void) G_GNUC_CONST;

NcHIPertBoltzmannCBE *nc_hipert_boltzmann_cbe_new (void);
NcHIPertBoltzmannCBE *nc_hipert_boltzmann_cbe_full_new (NcCBE *cbe);
NcHIPertBoltzmannCBE *nc_hipert_boltzmann_cbe_ref (NcHIPertBoltzmannCBE *boltzmann_cbe);
void nc_hipert_boltzmann_cbe_free (NcHIPertBoltzmannCBE *boltzmann_cbe);
void nc_hipert_boltzmann_cbe_clear (NcHIPertBoltzmannCBE **boltzmann_cbe);

NcCBE *nc_hipert_boltzmann_cbe_peek_cbe (NcHIPertBoltzmannCBE *boltzmann_cbe);

G_END_DECLS

#endif /* _NC_HIPERT_BOLTZMANN_CBE_H_ */

