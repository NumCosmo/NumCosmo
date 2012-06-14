/***************************************************************************
 *            linear_internal.h
 *
 *  Fri Nov  6 15:27:13 2009
 *  Copyright  2009  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@lapsandro>
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

#define _NC_USE_CUTOFF 1
//#define _NC_USE_CUTOFF_TEST (g > pert->g_opt_cutoff)
//#define _NC_USE_CUTOFF_TEST (!pert->pws->tight_coupling)
#define _NC_USE_CUTOFF_TEST TRUE
#define _NC_TIGHT_COUPLING_END 1.0e5

#define _NC_PHI (LINEAR_VEC_COMP(y, NC_PERT_PHI))
#define _NC_C0 (LINEAR_VEC_COMP(y, NC_PERT_C0))

#define _NC_B0 (LINEAR_VEC_COMP(y, NC_PERT_B0))
#define _NC_dB0 (LINEAR_VEC_COMP(y, NC_PERT_dB0))

#define _NC_THETA0 (LINEAR_VEC_COMP(y, NC_PERT_THETA0))
#define _NC_dTHETA0 (LINEAR_VEC_COMP(y, NC_PERT_dTHETA0))

#define _NC_C1 (LINEAR_VEC_COMP(y, NC_PERT_C1))
#define _NC_V (LINEAR_VEC_COMP(y, NC_PERT_V))

#define _NC_U (LINEAR_VEC_COMP(y, NC_PERT_U))
#define _NC_THETA1 (LINEAR_VEC_COMP(y, NC_PERT_THETA1))

#define _NC_T (LINEAR_VEC_COMP(y, NC_PERT_T))
#define _NC_B1 (LINEAR_VEC_COMP(y, NC_PERT_B1))

#define _NC_THETA2 (LINEAR_VEC_COMP(y, NC_PERT_THETA2))
#define _NC_THETA(n) (LINEAR_VEC_COMP(y, NC_PERT_THETA(n)))

#define _NC_THETA_P0 (LINEAR_VEC_COMP(y, NC_PERT_THETA_P0))
#define _NC_THETA_P1 (LINEAR_VEC_COMP(y, NC_PERT_THETA_P1))
#define _NC_THETA_P2 (LINEAR_VEC_COMP(y, NC_PERT_THETA_P2))
#define _NC_THETA_P(n) (LINEAR_VEC_COMP(y, NC_PERT_THETA_P(n)))

#define _NC_DPHI (LINEAR_VEC_COMP(ydot, NC_PERT_PHI))
#define _NC_DC0 (LINEAR_VEC_COMP(ydot, NC_PERT_C0))

#define _NC_DB0 (LINEAR_VEC_COMP(ydot, NC_PERT_B0))
#define _NC_DdB0 (LINEAR_VEC_COMP(ydot, NC_PERT_dB0))

#define _NC_DTHETA0 (LINEAR_VEC_COMP(ydot, NC_PERT_THETA0))
#define _NC_DdTHETA0 (LINEAR_VEC_COMP(ydot, NC_PERT_dTHETA0))

#define _NC_DV (LINEAR_VEC_COMP(ydot, NC_PERT_V))
#define _NC_DC1 (LINEAR_VEC_COMP(ydot, NC_PERT_C1))

#define _NC_DU (LINEAR_VEC_COMP(ydot, NC_PERT_U))
#define _NC_DTHETA1 (LINEAR_VEC_COMP(ydot, NC_PERT_THETA1))

#define _NC_DT (LINEAR_VEC_COMP(ydot, NC_PERT_T))
#define _NC_DB1 (LINEAR_VEC_COMP(ydot, NC_PERT_B1))

#define _NC_DTHETA2 (LINEAR_VEC_COMP(ydot, NC_PERT_THETA2))
#define _NC_DTHETA(n) (LINEAR_VEC_COMP(ydot, NC_PERT_THETA(n)))

#define _NC_DTHETA_P0 (LINEAR_VEC_COMP(ydot, NC_PERT_THETA_P0))
#define _NC_DTHETA_P1 (LINEAR_VEC_COMP(ydot, NC_PERT_THETA_P1))
#define _NC_DTHETA_P2 (LINEAR_VEC_COMP(ydot, NC_PERT_THETA_P2))
#define _NC_DTHETA_P(n) (LINEAR_VEC_COMP(ydot, NC_PERT_THETA_P(n)))
