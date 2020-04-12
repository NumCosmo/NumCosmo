/***************************************************************************
 *            ncm_spline_func.h
 *
 *  Wed Aug 13 21:13:59 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NCM_SPLINE_FUNC_H
#define _NCM_SPLINE_FUNC_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

/**
 * NcmSplineFuncType:
 * @NCM_SPLINE_FUNCTION_4POINTS: This function initially uses 4 knots and applies a cubic polynomial interpolation. It compares the polynomial interpolation and @F at the intermediate points between the knots. If the @rel_error is achieved in one of the points than it is saved. If it is not then it's division is again splitted in 4 knots and another polynomial interpolation is applied and the same comparison is done as before. This procedure is repeated until @rel_error is achieved for each division. Note that this procedure does not create a homogeneous grid of knots.   
 * @NCM_SPLINE_FUNCTION_SPLINE: FIXME
 * @NCM_SPLINE_FUNCTION_SPLINE_LNKNOT: FIXME
 * @NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT: FIXME
 *
 * Enumeration to choose which of the functions to be applied when interpolating the input #gsl_function *@F, $f$, with the desired @rel_error in the range [@xi, @xf]. The interpolation knots, $\mathbf{x}$, are automatically defined internally by the functions. 
 * 
 * All available algorithms initially start with $n_0$ knots, $\mathbf{x}_0$, equally spaced in the range [@xi, @xf], including both limiting points. The value of $n_0$ depends on the chosen interpolation method given by @s, except for @NCM_SPLINE_FUNCTION_4POINTS that starts with $n_0 = 4$. Below we show the common procedure applied by the functions @NCM_SPLINE_FUNCTION_SPLINE, @NCM_SPLINE_FUNCTION_SPLINE_LNKNOT and @NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT. After, this procedure will be adapted to the function @NCM_SPLINE_FUNCTION_4POINTS which distributes knots and creates the spline function in a slightly different manner.
 *
 * The function $f$ is first interpolated at the $\mathbf{x}_0$ knots, producing the $\hat{f}_0$ interpolated function. Next, the existing $n_0 - 1$ bins, $\Delta \mathbf{x}_0 = \mathbf{x}_0^{\mathrm{r}} - \mathbf{x}_0^{\mathrm{l}}$, between the $\mathbf{x}_0$ knots are divided in half and a new set of knots are created, $\overline{\mathbf{x}}_0$. Then the following tests are done for each bin $\Delta \mathbf{x}_0$,
 * \begin{equation*}
 *   \left| \frac{ \hat{f}_0(\overline{\mathbf{x}}_0) - f(\overline{\mathbf{x}}_0)}{f(\overline{\mathbf{x}}_0)} \right| < \mathrm{rel \\_ error}
 * \end{equation*}
 * and
 * \begin{equation*}
 *   \left| \frac{ \int_{\Delta \mathbf{x}_0} \hat{f}_0  - \int_{\Delta \mathbf{x}_0} f }{ \int_{\Delta \mathbf{x}_0} f } \right| < \mathrm{rel \\_ error}.
 * \end{equation*}
 * Where $\int_{\Delta \mathbf{x}_0} f$ is the integral of the input function $f$ evaluated using [Simpson's rule](https://en.wikipedia.org/wiki/Simpson%27s_rule)
 *\begin{equation*}
 * \int_{\Delta \mathbf{x}_0} f = \frac{\Delta \mathbf{x}_0}{6} \left[ f(\mathbf{x}_0^{\mathrm{l}}) + 4 f(\overline{\mathbf{x}}_0) + f(\mathbf{x}_0^{\mathrm{r}}) \right]
 *\end{equation*}
 *
 * and the interpolated function $\hat{f}_0$ is integrated by applying ncm_spline_eval_integ() function. Both conditions are verified for each one of the $\Delta \mathbf{x}_0$ bins, separately. If any bin passes those relations, then its associated $\mathbf{x}_0 \cup  \overline{\mathbf{x}}_0$ knots are defined as a good representation of the function $f$ and this specific bin does not need to be refined anymore. If not, then $\mathbf{x}_0 \cup \overline{\mathbf{x}}_0$ is splitted once again into two more symmetric knots around $\overline{\mathbf{x}}_0$. All the knots from the previous fase defines a new set of knots $\mathbf{x}_1 = \mathbf{x}_0 \cup \overline{\mathbf{x}}_0$ and the spllited ones, which did not pass the tests, define another set of knots $\overline{\mathbf{x}}_1$ that lies between the $\mathbf{x}_0\cup \overline{\mathbf{x}}_0$. The new set of knots $\mathbf{x}_1$ are used to create a new interpolated function $\hat{f}_1$, always for the full range [@xi, @xf]. The same tests are performed as before, but now with $\hat{f}_1(\overline{\mathbf{x}}_1)$, $f(\overline{\mathbf{x}}_1)$ and the integral has limits $\Delta \mathbf{x}_1 = \Delta \mathbf{x}_0/2$, but only for the bins that did not pass the previous test. This procedure is repeated until the desired accuracy is met across the range [@xi, @xf]. Note that it will most probably create a inhomogeneous set of knots. 
 *
 * <inlinegraphic fileref="/home/fsimoni/cosmology/NumCosmo/docs/images/spline_func_knots_evolution.png" format="PNG" scale="95" align="center"/>
 *
 *
 * The figure on the left shows a schematic evolution. First it starts with 6 knots, $\mathbf{x}_0$ (black filled circles), used to create the interpolated function $\hat{f}_0$.
 * The $\overline{\mathbf{x}}_0$ are created (blue squares) and the first tests are performed. In this case only the first and the fourth original bins did not pass the test. 
 * A new set of knots is create $\mathbf{x}_1$ (second line) and also $\overline{\mathbf{x}}_1$ (red diamonds). Again the tests are done and only two passes, the first and the fourth 
 * bins with $\overline{\mathbf{x}}_1$. A new set of knots are created $\mathbf{x}_2 = \mathbf{x}_1 \cup \overline{\mathbf{x}}_1$ and the procedure is repeated once more for 2 bins. 
 * In this schematic example, the final set of knots is given by the last line $\mathbf{x}_3$ with 19 knots in total, also showing that the final distribution is not homogeneous. 
 *
 *
 *
 * The specifics of each function are described below.
 *
 */
typedef enum _NcmSplineFuncType
{
  NCM_SPLINE_FUNCTION_4POINTS,
  NCM_SPLINE_FUNCTION_SPLINE,
  NCM_SPLINE_FUNCTION_SPLINE_LNKNOT,
  NCM_SPLINE_FUNCTION_SPLINE_SINHKNOT,
} NcmSplineFuncType;

void ncm_spline_set_func (NcmSpline *s, NcmSplineFuncType ftype, gsl_function *F, gdouble xi, gdouble xf, gsize max_nodes, gdouble rel_error);

#define NCM_SPLINE_FUNC_DEFAULT_MAX_NODES 10000000
#define NCM_SPLINE_KNOT_DIFF_TOL (GSL_DBL_EPSILON * 1.0e2)

G_END_DECLS

#endif /* _NCM_SPLINE_FUNC_H */

