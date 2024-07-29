/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/*
 * nc_hipert_boltzmann_std.h
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

#ifndef _NC_HIPERT_BOLTZMANN_STD_H_
#define _NC_HIPERT_BOLTZMANN_STD_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/perturbations/nc_hipert_boltzmann.h>

G_BEGIN_DECLS

#define NC_TYPE_HIPERT_BOLTZMANN_STD             (nc_hipert_boltzmann_std_get_type ())
#define NC_HIPERT_BOLTZMANN_STD(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIPERT_BOLTZMANN_STD, NcHIPertBoltzmannStd))
#define NC_HIPERT_BOLTZMANN_STD_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIPERT_BOLTZMANN_STD, NcHIPertBoltzmannStdClass))
#define NC_IS_HIPERT_BOLTZMANN_STD(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIPERT_BOLTZMANN_STD))
#define NC_IS_HIPERT_BOLTZMANN_STD_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIPERT_BOLTZMANN_STD))
#define NC_HIPERT_BOLTZMANN_STD_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIPERT_BOLTZMANN_STD, NcHIPertBoltzmannStdClass))

typedef struct _NcHIPertBoltzmannStdClass NcHIPertBoltzmannStdClass;
typedef struct _NcHIPertBoltzmannStd NcHIPertBoltzmannStd;

struct _NcHIPertBoltzmannStdClass
{
  /*< private >*/
  NcHIPertBoltzmannClass parent_class;
};

struct _NcHIPertBoltzmannStd
{
  /*< private >*/
  NcHIPertBoltzmann parent_instance;
};

GType nc_hipert_boltzmann_std_get_type (void) G_GNUC_CONST;

NcHIPertBoltzmannStd *nc_hipert_boltzmann_std_new (NcRecomb *recomb, guint lmax);

G_END_DECLS

#endif /* _NC_HIPERT_BOLTZMANN_STD_H_ */

