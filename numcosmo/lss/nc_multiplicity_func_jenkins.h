/***************************************************************************
 *            nc_multiplicity_func_st.h
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

#ifndef _NC_MULTIPLICITY_FUNC_JENKINS_H_
#define _NC_MULTIPLICITY_FUNC_JENKINS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_multiplicity_func.h>

G_BEGIN_DECLS

#define NC_TYPE_MULTIPLICITY_FUNC_JENKINS             (nc_multiplicity_func_jenkins_get_type ())
#define NC_MULTIPLICITY_FUNC_JENKINS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_MULTIPLICITY_FUNC_JENKINS, NcMultiplicityFuncJenkins))
#define NC_MULTIPLICITY_FUNC_JENKINS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_MULTIPLICITY_FUNC_JENKINS, NcMultiplicityFuncJenkinsClass))
#define NC_IS_MULTIPLICITY_FUNC_JENKINS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_MULTIPLICITY_FUNC_JENKINS))
#define NC_IS_MULTIPLICITY_FUNC_JENKINS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_MULTIPLICITY_FUNC_JENKINS))
#define NC_MULTIPLICITY_FUNC_JENKINS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_MULTIPLICITY_FUNC_JENKINS, NcMultiplicityFuncJenkinsClass))

typedef struct _NcMultiplicityFuncJenkinsClass NcMultiplicityFuncJenkinsClass;
typedef struct _NcMultiplicityFuncJenkins NcMultiplicityFuncJenkins;

struct _NcMultiplicityFuncJenkinsClass
{
  /*< private >*/
  NcMultiplicityFuncClass parent_class;
};

struct _NcMultiplicityFuncJenkins
{
  /*< private >*/
  NcMultiplicityFunc parent_instance;
  gdouble A;
  gdouble A_tCDM;
  gdouble B;
  gdouble B_tCDM;
  gdouble epsilon;
  gdouble epsilon_tCDM;
};

GType nc_multiplicity_func_jenkins_get_type (void) G_GNUC_CONST;

NcMultiplicityFunc *nc_multiplicity_func_jenkins_new (gdouble A, gdouble A_tCDM, gdouble B, gdouble B_tCDM, gdouble epsilon, gdouble epsilon_tCDM);
void nc_multiplicity_func_jenkins_set_A (NcMultiplicityFuncJenkins *mulf_jenkins, gdouble A);
gdouble nc_multiplicity_func_jenkins_get_A (const NcMultiplicityFuncJenkins *mulf_jenkins);
void nc_multiplicity_func_jenkins_set_A_tCDM (NcMultiplicityFuncJenkins *mulf_jenkins, gdouble A_tCDM);
gdouble nc_multiplicity_func_jenkins_get_A_tCDM (const NcMultiplicityFuncJenkins *mulf_jenkins);
void nc_multiplicity_func_jenkins_set_B (NcMultiplicityFuncJenkins *mulf_jenkins, gdouble B);
gdouble nc_multiplicity_func_jenkins_get_B (const NcMultiplicityFuncJenkins *mulf_jenkins);
void nc_multiplicity_func_jenkins_set_B_tCDM (NcMultiplicityFuncJenkins *mulf_jenkins, gdouble B_tCDM);
gdouble nc_multiplicity_func_jenkins_get_B_tCDM (const NcMultiplicityFuncJenkins *mulf_jenkins);
void nc_multiplicity_func_jenkins_set_epsilon (NcMultiplicityFuncJenkins *mulf_jenkins, gdouble epsilon);
gdouble nc_multiplicity_func_jenkins_get_epsilon (const NcMultiplicityFuncJenkins *mulf_jenkins);
void nc_multiplicity_func_jenkins_set_epsilon_tCDM (NcMultiplicityFuncJenkins *mulf_jenkins, gdouble epsilon_tCDM);
gdouble nc_multiplicity_func_jenkins_get_epsilon_tCDM (const NcMultiplicityFuncJenkins *mulf_jenkins);

G_END_DECLS

#endif /* _NC_MULTIPLICITY_FUNC_JENKINS_H_ */
