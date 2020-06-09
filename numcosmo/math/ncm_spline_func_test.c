/***************************************************************************
 *            ncm_spline_func_test.c
 *
 *  Wed March 14 16:30:36 2020
 *  Copyright  2020 Fernando de Simoni
 *  <fsimoni@id.uff.br>
 ****************************************************************************/
/*
 * ncm_spline_func_test.c
 * Copyright (C) 2020 Fernando de Simoni <fsimoni@id.uff.brr>
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
 * SECTION:ncm_spline_func_test
 * @title: NcmSplineFuncTest
 * @short_description: test suite to analyze the NcmSplineFunc's knots distribution.
 * @stability: Private and Unstable
 * @include: numcosmo/math/ncm_spline_func_test.h
 *
 * This module is intended to be a test suite for the [NcmSplineFunc](numcosmo-NcmSplineFunc.html) object. 
 * It performs a brute force approach to check if the required tolerance (#NcmSplineFuncTest:rel-error and #NcmSplineFuncTest:scale) 
 * was achieved throughout the entire desired range [#NcmSplineFuncTest:xi, #NcmSplineFuncTest:xf] 
 * to evaluate a base function $f(x)$. To perform the test, it is created an evenly distributed grid 
 * with #NcmSplineFuncTest:ngrid knots, including both limiting points.
 *
 * The criteria applied by [NcmSplineFunc](numcosmo-NcmSplineFunc.html) to include a knot is given by the condition,
 * \begin{equation}\label{eq:condition}
 * |f(x) - \hat{f}(x)| < \mathrm{rel \\_ error} \times \left[ |f(x)| + \mathrm{scale} \right] \,\, .
 * \end{equation}
 * Where, $f(x)$ is the base (true) function to be analyzed and $\hat{f}(x)$ is the estimation 
 * given by [NcmSplineFunc](numcosmo-NcmSplineFunc.html) using the interpolation method #NcmSplineCubicNotaknot.
 *
 * The test can be done on any function provided by the user (#NCM_SPLINE_FUNC_TEST_TYPE_USER), 
 * but it is also a stress test tool with some built in base functions (see #NcmSplineFuncTestType):
 *
 * - #NCM_SPLINE_FUNC_TEST_TYPE_POLYNOMIAL: it applies a polynomial interpolation with the desired degree upon a given set of points. 
 *   The polynomial degree is given by the number of points provided by the user subtracted by one.
 *   It uses [GSL](https://www.gnu.org/software/gsl/)'s polynomial interpolation method encapsulated with #NcmSplineGsl.
 *
 * - #NCM_SPLINE_FUNC_TEST_TYPE_COSINE: it applies a summation of cosine functions given by
 *   \begin{equation*}
 *     f(x) = \sum_{i=0}^{N-1} A_{i} \cos \left( 2 \pi \, \nu_i \, x  \right) \,\, .
 *   \end{equation*}
 *   The user provides the amplitudes $A_i$ and the frequencies $\nu_i$.
 *
 * - #NCM_SPLINE_FUNC_TEST_TYPE_RBF: it is similar to the polynomial, but now it applies a RBF interpolation upon a given set of points. 
 *
 * For both interpolated base functions, polynomial and RBF, all ordinate 
 * points $y$ can be drawn from a flat or a normal (gaussian) distribution.
 * The abscissa points $x$ are always drawn from a flat distribution, 
 * and its limiting points are fixed at #NcmSplineFuncTest:xi and #NcmSplineFuncTest:xf.
 *
 * The randomness goal is to produce a base function with several shapes. 
 * The #NcmSplineFuncTestTypePDF type provides the switch between both available PDFs. 
 *
 * # Setting the parameters # {#set-par}
 * In order to provide the base function parameters, the user has two options: 
 * 
 * 1. #ncm_spline_func_test_set_params_info(): the user must supply a #NcmMatrix with the number of rows as the number of parameters and two columns.
 *    If it is a flat PDF, each column represents the minimum and maximum parameter value possible to be drawn. 
 *    If it is a normal (gaussian) PDF, each column represents its mean and standard deviation. The supplied matrix is copied to #NcmSplineFuncTest:par-info. 
 *
 * 2. #ncm_spline_func_test_set_params_info_all(): the user must supply the number of parameters and the same value for each column. 
 *    Therefore, all parameters will have the same statistical properties. It creates the matrix #NcmSplineFuncTest:par-info internally.
 *
 * <emphasis>Fixing a parameter:</emphasis> the user also have the option to fix some or all the parameters. 
 * 
 * * #NCM_SPLINE_FUNC_TEST_TYPE_PDF_FLAT: set the parameter's minimum and maximum to the same value (both columns with the same value).
 *
 * * #NCM_SPLINE_FUNC_TEST_TYPE_PDF_NORMAL: set the parameter's standard deviation to zero (second column equal to zero).
 *
 * <emphasis>Built in functions parameters:</emphasis>
 *
 * * For the two interpolation methods, polynomial and RBF, the provided parameters are the ordinate points $y$.
 *
 * * For the cosine summation method, the provided parameters are the amplitude and the frequency, alternating in that order. 
 *   If someone wants to perform a sum of two cosine functions, must provide a matrix with four rows: $A_0$, $\nu_0$, $A_1$ and $\nu_1$, for example. 
 * <note>
 *   <para>
 *     The #NCM_SPLINE_FUNC_TEST_TYPE_COSINE input parameters matrix must have an even number of rows. 
 *   </para>
 * </note>
 *
 *
 * # One grid examples # {#grid-ex}
 *
 * To perform one grid statistics, it is needed to place ncm_spline_func_test_set_one_grid_stats() after preparing it (ncm_spline_func_test_prepare()).
 * Two examples are shown below together with their results.
 *
 * ## Example: 7th degree polynomial interpolation. # {#7th-ex}
 * |[<!-- language="C" -->
 * #include <numcosmo/numcosmo.h>
 *
 * int
 * main (void)
 * {
 *   guint npar = 8, seed = 1, ngrid = 100000;
 *
 *   gdouble rel_error = 1.e-10, scale = 1.0, mean = 0.0, sigma = 1.0;
 *
 *   NcmSplineFuncTest *sft = ncm_spline_func_test_new ();
 *
 *   ncm_spline_func_test_set_seed (sft, seed);
 *
 *   ncm_spline_func_test_set_ngrid (sft, ngrid);
 *
 *   ncm_spline_func_test_set_scale (sft, scale);
 *
 *   ncm_spline_func_test_set_rel_error (sft, rel_error);
 *
 *   ncm_spline_func_test_set_params_info_all (sft, npar, mean, sigma);
 *
 *   ncm_spline_func_test_prepare (sft, NCM_SPLINE_FUNCTION_SPLINE, NCM_SPLINE_FUNC_TEST_TYPE_PDF_NORMAL);
 *
 *   ncm_spline_func_test_set_one_grid_stats (sft);
 *
 *   ncm_spline_func_test_log_vals_one_grid_stats (sft);
 *
 *   ncm_spline_func_test_save_grid_functions_to_txt (sft, "functions.txt");
 *
 *   ncm_spline_func_test_save_knots_to_txt (sft, "knots.txt");
 *
 *   ncm_spline_func_test_unref (sft);
 *
 *   return 0;
 * }
 * ]|
 *
 * The function ncm_spline_func_test_set_params_info_all() creates the 8 rows matrix, 
 * “#NcmSplineFuncTest:par-info”, the number of ordinate points $y$ (consequently 8 abscissa $x$), and two columns.
 * The first column is $\mu=0$ and the second $\sigma = 1$ with 
 * the same value for all parameters (rows) to be used by the gaussian PDF.
 *
 * The function ncm_spline_func_test_prepare() sets the [NcmSplineFunc](numcosmo-NcmSplineFunc.html) method
 * used, [NCM_SPLINE_FUNCTION_SPLINE](numcosmo-NcmSplineFunc.html), and also sets the PDF, #NCM_SPLINE_FUNC_TEST_TYPE_PDF_NORMAL.  
 *
 * The functions ncm_spline_func_test_save_grid_functions_to_txt() and ncm_spline_func_test_save_knots_to_txt() 
 * create the ascii files "functions.txt", with all the information about the functions at all grid, 
 * and "knots.txt", which saves [NcmSplineFunc](numcosmo-NcmSplineFunc.html)'s knots.
 *
 * The function ncm_spline_func_test_log_vals_one_grid_stats() prints 
 * the statistics evaluated by ncm_spline_func_test_set_one_grid_stats().
 * In the above example it displays:
 * <informalexample>
 *   <programlisting>
 *
 * ##### Grid statistics  #####
 *
 * NcmSplineFunc number of knots = 3847 
 *
 * (ncm) diff. = -7.377e-13 +/- 2.582e-11 (abs. max. diff. = 2.917e-10) | outliers =  0.00 %
 * (lin) diff. = -3.826e-13 +/- 3.091e-11 (abs. max. diff. = 1.211e-09) | outliers =  0.09 %
 *
 *   </programlisting>
 * </informalexample>
 * * "NcmSplineFunc number of knots = 3847": is self explanatory.
 * * "diff." first value (-7.377e-13 and -3.826e-13): it is the [mean signed difference](https://en.wikipedia.org/wiki/Mean_signed_deviation) 
 *   between the base and approximated functions, $f(x)$ and $\hat{f}(x)$. 
 *   The signed difference is defined as $\Delta f(x)= f(x) - \hat{f}(x)$. Therefore its mean value through all grid knots should
 *   be as close to zero as possible, regardless of the required tolerance.
 * * "diff." second value (2.582e-11 and 3.091e-11): it is the [standard deviation](https://en.wikipedia.org/wiki/Standard_deviation) 
 *   of $\Delta f(x)$. This value is related to the required tolerance.
 *   If $\Delta f(x)$ distribution is assumed to be normal, which is not, we have $5\sigma \approx 10^{-10}$, the required #NcmSplineFuncTest:rel-error in this example.   
 * * "abs. max. diff.": it is the maximum absolute difference between the entire grid knots, $|\Delta f(x)|_{\mathrm{max}}$. 
 *   In this example, it appears that "ncm" has a higher value than the required tolerance, 
 *   but in fact it is not true, because @scale is not zero.
 * * "outliers": it is defined as the knots in the linear grid that did not passed the criteria given by Eq. \eqref{eq:condition}. 
 *   It is given in percentage of the number of knots (#NcmSplineFuncTest:ngrid).
 *
 *
 * ## Example: cosine summation with 50 terms.
 * |[<!-- language="C" -->
 * #include <numcosmo/numcosmo.h>
 *
 * int
 * main (void)
 * {
 *   NcmVector *c           = ncm_vector_new (50);
 *   NcmMatrix *params      = ncm_matrix_new (50, 2);
 *   NcmSplineFuncTest *sft = ncm_spline_func_test_new ();
 *
 *   ncm_vector_set_all (c, -5.0);
 *   ncm_matrix_set_col (params, 0, c);
 *   ncm_vector_set_all (c, +5.0);
 *   ncm_matrix_set_col (params, 1, c);
 *
 *   ncm_matrix_set (params, 0, 0, 100.); // Sets the first
 *   ncm_matrix_set (params, 0, 1, 100.); // amplitude fixed: A_0 = 100.
 *   ncm_matrix_set (params, 1, 0, 0.);   // Sets the first         
 *   ncm_matrix_set (params, 1, 1, 0.);   // frequency fixed: \nu_0 = 0. 
 *
 *   ncm_spline_func_test_set_type (sft, NCM_SPLINE_FUNC_TEST_TYPE_COSINE);
 *
 *   ncm_spline_func_test_set_params_info (sft, params);
 *
 *   ncm_spline_func_test_set_ngrid (sft, 1000000);
 *
 *   ncm_spline_func_test_set_seed (sft, 73649);
 *
 *   ncm_spline_func_test_prepare (sft, NCM_SPLINE_FUNCTION_4POINTS, NCM_SPLINE_FUNC_TEST_TYPE_PDF_FLAT);
 *
 *   ncm_spline_func_test_set_one_grid_stats (sft);
 *
 *   ncm_spline_func_test_log_vals_one_grid_stats (sft);
 *
 *   ncm_vector_free (c);
 *   ncm_matrix_free (params);
 *   ncm_spline_func_test_unref (sft);
 *
 *   return 0;
 * }
 * ]|
 *
 * In this example the user created a matrix with 50 parameters, 25 amplitudes and 25 frequencies. 
 * But note that the first amplitude is fixed to $A_0 = 100$ and the first frequency to $\nu_0 = 0$.
 * Therefore, creating a base function with some kind of oscilatory behaviour added to a constant factor,
 * \begin{equation*}
 *  f(x) = A_0 + \sum_{i=1}^{49} A_i \cos \left( 2 \pi \, \nu_i \, x  \right) \,\,\, , 
 * \end{equation*}
 * with $A_i$ and $\nu_i$ drawn from a flat PDF between the values [-5, 5]. 
 * Below is shown the output from ncm_spline_func_test_log_vals_one_grid_stats()
 * <informalexample>
 *   <programlisting>
 *
 * ##### Grid statistics  #####
 *
 * NcmSplineFunc number of knots = 7318 
 *
 * (ncm) diff. =  2.050e-13 +/- 3.565e-12 (abs. max. diff. = 5.193e-11) | outliers =  0.00 %
 * (lin) diff. = -1.088e-14 +/- 2.517e-12 (abs. max. diff. = 7.323e-11) | outliers =  0.00 %
 *
 *   </programlisting>
 * </informalexample>
 *
 * Those statistics are remarkably better than the required tolerance given by Eq. \eqref{eq:condition}. 
 * In this case we have the default values, $\mathrm{rel \\_ error} = 10^{-8}$ and $\mathrm{scale} = 0$.
 * The result condition is around $|f(x)| \times \mathrm{rel \\_ error} \approx 10^{-6}$. 
 * Note that both maximum absolute error are given by $|\Delta f(x)|_{\mathrm{max}} \approx 10^{-10}$, four orders of magnitude better than expected.
 * This fact is due to the [NcmSplineFunc](numcosmo-NcmSplineFunc.html) method applied in this example, #NCM_SPLINE_FUNCTION_4POINTS.
 * It is a much more conservative approach compared to #NCM_SPLINE_FUNCTION_SPLINE, applied in the previous example.   
 *
 * # Monte Carlo examples # {#mc-ex}
 *
 * To perform a Monte Carlo statistics it is needed to place ncm_spline_func_test_monte_carlo_and_save_to_txt() 
 * or ncm_spline_func_test_monte_carlo() after preparing it (ncm_spline_func_test_prepare()).
 *
 * Two examples are shown below, together with their results.
 *
 * ## Example: 5th degree polynomial interpolation with 1 million realizations.
 *
 * |[<!-- language="C" -->
 * #include <numcosmo/numcosmo.h>
 *
 * int
 * main (void)
 * {
 *   guint npar = 6, nsim = 1000000;
 *
 *   NcmSplineFuncTest *sft = ncm_spline_func_test_new ();
 *
 *   ncm_spline_func_test_set_params_info_all (sft, npar, -10.0, 10.0);
 *
 *   ncm_spline_func_test_set_scale (sft, 1.0);
 *
 *   ncm_spline_func_test_set_ngrid (sft, 10000);
 *
 *   ncm_spline_func_test_prepare (sft, NCM_SPLINE_FUNCTION_SPLINE, NCM_SPLINE_FUNC_TEST_TYPE_PDF_FLAT);
 *
 *   ncm_spline_func_test_monte_carlo_and_save_to_txt (sft, nsim, "mc.txt");
 *
 *   ncm_spline_func_test_log_vals_mc_stats (sft);
 *
 *   ncm_spline_func_test_unref (sft);
 *
 *   return 0;
 * }
 * ]|
 *
 * The function ncm_spline_func_test_monte_carlo_and_save_to_txt() creats the ascii file "mc.txt".
 * It saves the same statistics as printed by the function ncm_spline_func_test_log_vals_one_grid_stats () for each realization.
 * Therefore it is going to be quite a big file. In this example it has 165 megabytes. 
 * That is the reason why it is not allowed to save all the grids informations and also because the file is filled on the fly.
 *
 * Below is shown the output of ncm_spline_func_test_log_vals_mc_stats():
 * <informalexample>
 *   <programlisting>
 *
 * ####  Monte Carlo statistics for 1000000 simulations  ####
 *
 *  * NcmSplineFunc number of knots:   750.58 +/-   114.83
 *
 *  * Ncm clean grid (%): 97.59
 *  * Lin clean grid (%): 27.78
 *
 *  * Ncm outliers (%):  0.03 +/-  0.24
 *  * Lin outliers (%):  0.25 +/-  0.56
 *
 *  * Ncm diff. : 1.624e-07 +/- 6.353e-06
 *  * Lin diff. : 7.880e-09 +/- 2.385e-07
 *
 *  * Ncm abs. max. diff. : 3.891e-05 +/- 9.496e-03
 *  * Lin abs. max. diff. : 5.886e-06 +/- 1.099e-03
 *
 *   </programlisting>
 * </informalexample>
 *
 * First, we need to define two different statistics, "one grid" and "Monte Carlo". 
 * The mean value for each will be represented by $\overline{X}$ and $\langle X \rangle$, respectively.
 * The standard deviation $\sigma$ is applied to the Monte Carlo realizations.
 *
 * * "NcmSplineFunc number of knots:" $\langle \mathrm{spline \\_ length}  \rangle \pm  \sigma \left( \mathrm{spline \\_ length} \right)$.
 * * "clean grid (%):" it is the proportion of the grids realizations with no outliers at all.
 * * "outliers (%):" $\Big \langle \overline{\mathrm{outlier}} \Big  \rangle \pm  \sigma \left( \overline{\mathrm{outlier}} \right)$.
 * * "diff.:" $\Big \langle \overline{\Delta f(x)} \Big \rangle \pm  \sigma \left( \overline{\Delta f(x)} \right)$.
 * * "abs. max. diff.:" $\Big \langle |\Delta f(x)|_{\mathrm{max}} \Big \rangle \pm  \sigma \left( |\Delta f(x)|_{\mathrm{max}} \right)$.
 * 
 * Note that that the "clean grid" shows a much better result for [NcmSplineFunc](numcosmo-NcmSplineFunc.html) over the linear grid.
 * Also "outliers" shows the same trend but not at the same level.
 * Nevertheless, the "diff." and "abs. max. diff." seems to contradict both of them. 
 * To check this apparent contradiction, it is needed to look into the "mc.txt" file, which should show that the ~2% "non-clean"
 * [NcmSplineFunc](numcosmo-NcmSplineFunc.html) $\hat{f}(x)$ should have outliers with high values, probably indicating some issue with those functions.
 *
 * ## Example: user supplied function with 100 thousand realizations.
 *
 * |[<!-- language="C" -->
 * #include <numcosmo/numcosmo.h>
 *
 * gdouble
 * f_user (gdouble x, gpointer pin)
 * {
 *   NcmSplineFuncTest *sft = NCM_SPLINE_FUNC_TEST (pin);
 *
 *   NcmVector *params = ncm_spline_func_test_peek_current_params (sft);
 *
 *   gdouble p1 = ncm_vector_fast_get (params, 0);
 *   gdouble p2 = ncm_vector_fast_get (params, 1);
 *
 *   return exp (sin (p1 * x)) * sin (x) * sin (x) + p2;
 * }
 *
 * int 
 * main (void)
 * {
 *   gsl_function F;
 *
 *   NcmSplineFuncTest *sft = ncm_spline_func_test_new ();
 *
 *   ncm_spline_func_test_set_type (sft, NCM_SPLINE_FUNC_TEST_TYPE_USER);
 *  
 *   ncm_spline_func_test_set_rel_error (sft, 1.e-10);
 *   ncm_spline_func_test_set_out_threshold (sft, 10.0);
 * 
 *   ncm_spline_func_test_set_xi (sft, -0.5);
 *   ncm_spline_func_test_set_xf (sft, +0.5);
 *
 *   ncm_spline_func_test_set_params_info_all (sft, 2, 0., 10.);
 *
 *   F.function = &f_user;
 *   F.params   = sft;
 *
 *   ncm_spline_func_test_set_user_gsl_function (sft, &F);
 *
 *   ncm_spline_func_test_prepare (sft, NCM_SPLINE_FUNCTION_SPLINE, NCM_SPLINE_FUNC_TEST_TYPE_PDF_NORMAL);
 *
 *   ncm_spline_func_test_monte_carlo (sft, 100000);
 *
 *   ncm_spline_func_test_log_vals_mc_stats (sft);
 *
 *   ncm_spline_func_test_unref (sft);
 *
 *   return 0;
 * }
 * ]|
 *
 * The user function "f_user" has two parameters, $p_1$ and $p_2$, and is given by
 * \begin{equation*}
 *  f(x) = \exp \left[ \sin(p_1 \, x)  \right] \sin^{2}(x) + p_2 \,\, .
 * \end{equation*}
 * Unlike the previous example, the Monte Carlo statistic is not saved in a file. 
 * Instead, it is used the ncm_spline_func_test_monte_carlo() function.
 * This procedure has a faster execution time, but the drawback is information waste.
 * The function ncm_spline_func_test_set_out_threshold() saves information only
 * for the functions with outilers above 10\% of the threshold given in Eq. \eqref{eq:condition}.
 *
 * Below is shown the output of ncm_spline_func_test_log_vals_mc_stats():
 * <informalexample>
 *   <programlisting>
 *
 * ####  Monte Carlo statistics for 100000 simulations  ####
 *
 *  * NcmSplineFunc number of knots:   738.92 +/-   451.61
 *
 *  * Ncm clean grid (%): 83.33
 *  * Lin clean grid (%): 96.18
 *
 *  * Ncm outliers (%):  0.16 +/-  0.54
 *  * Lin outliers (%):  0.00 +/-  0.02
 *
 *  * Ncm diff. : 1.926e-11 +/- 3.332e-10
 *  * Lin diff. : -9.121e-13 +/- 2.690e-11
 *
 *  * Ncm abs. max. diff. : 2.619e-09 +/- 3.357e-07
 *  * Lin abs. max. diff. : 4.122e-10 +/- 4.484e-10
 *
 *   </programlisting>
 * </informalexample>
 *
 * The results shows a better approximation by the linear grid 
 * compared to the [NcmSplineFunc](numcosmo-NcmSplineFunc.html) method.
 * In order to try to understand this behaviour, the user should look into the files
 * created by ncm_spline_func_test_set_out_threshold(), but note that "Ncm abs. max. diff." has mean 
 * near the desired tolerance. 
 *
 *
 * # Future improvements
 *
 * 1. Add OpenMP on grid and/or Monte Carlo loops.
 *
 * 2. Add @nsim parameter as an object property.
 *
 * 3. Add option for the user to set the $x$ values as input for the spline methods, polynomial and RBF.
 *
 * 4. Add a function that returns to the user the parameters that created the base function.
 *
 * # Known bugs
 *
 * 1. In the last example with the user function, if the method from [NcmSplineFunc](numcosmo-NcmSplineFunc.html)
 *    is changed to #NCM_SPLINE_FUNCTION_4POINTS, after some realizations it gives "Segmentation fault (core dumped)". 
 *    But does not have memory issue anymore.
 *
 * 2. When #NCM_SPLINE_FUNC_TEST_TYPE_POLYNOMIAL is used with different values from 
 *    the default for #NcmSplineFuncTest:xi and #NcmSplineFuncTest:xf it become very unstable.
 *    Rarely it even can not perform the ncm_spline_func_test_set_one_grid_stats(), but if it
 *    is able to enter the Monte Carlo simulation, it crashes within very few realizations.
 *    It happens with both [NcmSplineFunc](numcosmo-NcmSplineFunc.html) methods and PDFs.
 *    Curiously, with the default values, 0 and 1, it is stable.
 *    The displayed error is "gsl: interp.c:150: ERROR: interpolation error".
 *       
 * 3. When #NCM_SPLINE_FUNC_TEST_TYPE_RBF is used it becomes very unstable. Rarely it can pass the function
 *    ncm_spline_func_test_set_one_grid_stats(). But when it is able to, it can produce a series of errors and warnings listed below:
 *
 *   * NUMCOSMO-WARNING: ncm_spline_new_function_spline: cannot archive requested precision with at most 5000000 nodes.
 *     But it seems to not happen with #NCM_SPLINE_FUNCTION_4POINTS method.
 *   * Segmentation fault (core dumped): sometimes without even passing by ncm_spline_func_test_set_one_grid_stats(),
 *     but when enter the Monte Carlo, it happens after very few realizations (the highest value checked was ~150 realizations).  
 *   * NUMCOSMO-ERROR: _ncm_spline_rbf_prepare[ncm_matrix_cholesky_solve]: 7. (seems to happen very rarely).
 *   * NUMCOSMO-ERROR: Tolerance of the difference between knots was reached. Interpolated function is probably discontinuous (seems to happen rarely).
 *   * It is interesting to note that this method seems to have no memory issue.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spline_func_test.h"
#include "math/ncm_spline_func.h"
#include "math/ncm_spline.h"
#include "math/ncm_spline_gsl.h"
#include "math/ncm_spline_rbf.h"
#include "math/ncm_util.h"
#include "math/ncm_rng.h"
#include "math/ncm_vector.h"
#include "math/ncm_matrix.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_stats_vec.h"
#include <gsl/gsl_sort_vector.h>
#include "ncm_enum_types.h"
#include "math/ncm_cfg.h"

struct _NcmSplineFuncTestPrivate
{
  NcmSplineFuncTestType type;
  NcmSplineFuncTestTypePDF pdftype;
  NcmSplineFuncType ftype;
  
  guint ngrid;
  gdouble xi;
  gdouble xf;
  NcmVector *xgrid;
  NcmVector *ygrid_true;
  NcmVector *ygrid_ncm;
  NcmVector *ygrid_lin;
  NcmVector *diff_ncm;
  NcmVector *diff_lin;
  
  gdouble rel_error;
  gdouble scale;
  
  guint npar;
  NcmVector *x_spl;
  NcmVector *params;
  NcmMatrix *par_info;
  GArray *par_is_fixed;
  gboolean all_fixed;
  
  guint ncm_good;
  guint lin_good;
  guint ncm_out;
  guint lin_out;
  gdouble out_threshold;
  NcmStatsVec *stats_grid;
  NcmStatsVec *stats_mc;
  
  guint len;
  NcmSpline *ncm;
  NcmSpline *lin;
  
  guint nsim;
  gulong seed;
  NcmRNG *rng;
  
  gsl_function F;
};

enum
{
  PROP_0,
  PROP_TYPE,
  PROP_NGRID,
  PROP_SEED,
  PROP_PAR_INFO,
  PROP_XI,
  PROP_XF,
  PROP_REL_ERROR,
  PROP_SCALE,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmSplineFuncTest, ncm_spline_func_test, G_TYPE_OBJECT);

static void
ncm_spline_func_test_init (NcmSplineFuncTest *sft)
{
  NcmSplineFuncTestPrivate * const self = sft->priv = ncm_spline_func_test_get_instance_private (sft);
  
  self->type    = 0;
  self->pdftype = 0;
  self->ftype   = 0;
  
  self->ngrid      = 0;
  self->xi         = 0.0;
  self->xf         = 0.0;
  self->xgrid      = NULL;
  self->ygrid_true = NULL;
  self->ygrid_ncm  = NULL;
  self->ygrid_lin  = NULL;
  self->diff_ncm   = NULL;
  self->diff_lin   = NULL;
  
  self->rel_error = 0.0;
  self->scale     = 0.0;
  
  self->npar         = 0;
  self->x_spl        = NULL;
  self->params       = NULL;
  self->par_info     = NULL;
  self->par_is_fixed = g_array_new (FALSE, FALSE, sizeof (gboolean));
  self->all_fixed    = TRUE;
  
  self->ncm_good      = 0;
  self->lin_good      = 0;
  self->ncm_out       = 0;
  self->lin_out       = 0;
  self->out_threshold = 0.0;
  self->stats_grid    = NULL;
  self->stats_mc      = NULL;
  
  self->len = 0;
  self->ncm = ncm_spline_cubic_notaknot_new ();
  self->lin = ncm_spline_cubic_notaknot_new ();
  
  self->nsim = 0;
  self->seed = 0;
  self->rng  = ncm_rng_new (NULL);
  
  self->F.function = NULL;
  self->F.params   = NULL;
}

static void
_ncm_spline_func_test_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSplineFuncTest *sft = NCM_SPLINE_FUNC_TEST (object);
  
  g_return_if_fail (NCM_IS_SPLINE_FUNC_TEST (object));
  
  switch (prop_id)
  {
    case PROP_TYPE:
      ncm_spline_func_test_set_type (sft, g_value_get_enum (value));
      break;
    case PROP_NGRID:
      ncm_spline_func_test_set_ngrid (sft, g_value_get_uint (value));
      break;
    case PROP_SEED:
      ncm_spline_func_test_set_seed (sft, g_value_get_ulong (value));
      break;
    case PROP_PAR_INFO:
      ncm_spline_func_test_set_params_info (sft, g_value_get_object (value));
      break;
    case PROP_XI:
      ncm_spline_func_test_set_xi (sft, g_value_get_double (value));
      break;
    case PROP_XF:
      ncm_spline_func_test_set_xf (sft, g_value_get_double (value));
      break;
    case PROP_REL_ERROR:
      ncm_spline_func_test_set_rel_error (sft, g_value_get_double (value));
      break;
    case PROP_SCALE:
      ncm_spline_func_test_set_scale (sft, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_spline_func_test_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSplineFuncTest *sft                = NCM_SPLINE_FUNC_TEST (object);
  NcmSplineFuncTestPrivate * const self = ncm_spline_func_test_get_instance_private (sft);
  
  g_return_if_fail (NCM_IS_SPLINE_FUNC_TEST (object));
  
  switch (prop_id)
  {
    case PROP_TYPE:
      g_value_set_enum (value, self->type);
      break;
    case PROP_NGRID:
      g_value_set_uint (value, ncm_spline_func_test_get_ngrid (sft));
      break;
    case PROP_SEED:
      g_value_set_ulong (value, ncm_spline_func_test_get_seed (sft));
      break;
    case PROP_PAR_INFO:
      g_value_take_object (value, ncm_spline_func_test_get_params_info (sft));
      break;
    case PROP_XI:
      g_value_set_double (value, ncm_spline_func_test_get_xi (sft));
      break;
    case PROP_XF:
      g_value_set_double (value, ncm_spline_func_test_get_xf (sft));
      break;
    case PROP_REL_ERROR:
      g_value_set_double (value, ncm_spline_func_test_get_rel_error (sft));
      break;
    case PROP_SCALE:
      g_value_set_double (value, ncm_spline_func_test_get_scale (sft));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_spline_func_test_dispose (GObject *object)
{
  NcmSplineFuncTest *sft                = NCM_SPLINE_FUNC_TEST (object);
  NcmSplineFuncTestPrivate * const self = sft->priv;
  
  ncm_vector_clear (&self->xgrid);
  ncm_vector_clear (&self->ygrid_true);
  ncm_vector_clear (&self->ygrid_ncm);
  ncm_vector_clear (&self->ygrid_lin);
  ncm_vector_clear (&self->diff_ncm);
  ncm_vector_clear (&self->diff_lin);
  
  ncm_vector_clear (&self->x_spl);
  ncm_vector_clear (&self->params);
  ncm_matrix_clear (&self->par_info);
  
  ncm_stats_vec_clear (&self->stats_mc);
  ncm_stats_vec_clear (&self->stats_grid);
  
  if (self->type == NCM_SPLINE_FUNC_TEST_TYPE_POLYNOMIAL)
  {
    NcmSpline *pol = NCM_SPLINE (self->F.params);
    
    ncm_spline_clear (&pol);
  }
  else if (self->type == NCM_SPLINE_FUNC_TEST_TYPE_RBF)
  {
    NcmSplineRBF *rbf = NCM_SPLINE_RBF (self->F.params);
    
    ncm_spline_rbf_clear (&rbf);
  }
  
  ncm_spline_clear (&self->ncm);
  ncm_spline_clear (&self->lin);
  
  ncm_rng_clear (&self->rng);
  
  g_array_unref (self->par_is_fixed);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline_func_test_parent_class)->dispose (object);
}

static void
_ncm_spline_func_test_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_spline_func_test_parent_class)->finalize (object);
}

static void
ncm_spline_func_test_class_init (NcmSplineFuncTestClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);
  
  object_class->set_property = &_ncm_spline_func_test_set_property;
  object_class->get_property = &_ncm_spline_func_test_get_property;
  object_class->dispose      = &_ncm_spline_func_test_dispose;
  object_class->finalize     = &_ncm_spline_func_test_finalize;
  
  /**
   * NcmSplineFuncTest:type:
   *
   * The method applied to test the [NcmSplineFunc](numcosmo-NcmSplineFunc.html) object.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_TYPE,
                                   g_param_spec_enum ("type",
                                                      NULL,
                                                      "Type",
                                                      NCM_TYPE_SPLINE_FUNC_TEST_TYPE, NCM_SPLINE_FUNC_TEST_TYPE_POLYNOMIAL,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcmSplineFuncTest:ngrid:
   *
   * The number of knots to evaluate the comparison between the base function with the
   * the spline from [NcmSplineFunc](numcosmo-NcmSplineFunc.html) and the spline created with a homogeneous 
   * distribution of knots but with the same number as [NcmSplineFunc](numcosmo-NcmSplineFunc.html).
   *
   * The grid knots are evenly distributed in the range [@xi, @xf].
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_NGRID,
                                   g_param_spec_uint ("ngrid",
                                                      NULL,
                                                      "Number of grid nodes",
                                                      10, G_MAXUINT32, 1000,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcmSplineFuncTest:seed:
   *
   * The seed for the random number generator @rng.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SEED,
                                   g_param_spec_ulong ("seed",
                                                       NULL,
                                                       "RNG seed",
                                                       0, G_MAXUINT32, 0,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcmSplineFuncTest:par-info:
   *
   * #NcmMatrix with the test function parameters information.
   *
   * - Rows = number of parameters.
   * - Columns = 2, with the information to perform Monte Carlo: (min, max) for a flat PDF or (mean, sigma) for a normal PDF.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_PAR_INFO,
                                   g_param_spec_object ("par-info",
                                                        NULL,
                                                        "Test function parameters information",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcmSplineFuncTest:xi:
   *
   * The initial abscissa value @xi to performe the comparison.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_XI,
                                   g_param_spec_double ("xi",
                                                        NULL,
                                                        "Initial abscissa value",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  /**
   * NcmSplineFuncTest:xf:
   *
   * The final abscissa value @xf to performe the comparison.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_XF,
                                   g_param_spec_double ("xf",
                                                        NULL,
                                                        "Final abscissa value",
                                                        -G_MAXDOUBLE, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
/**
 * NcmSplineFuncTest:rel-error:
 *
 * Relative error (tolerance) used by [NcmSplineFunc](numcosmo-NcmSplineFunc.html) object when defining its knots.
 *
 */
  g_object_class_install_property (object_class,
                                   PROP_REL_ERROR,
                                   g_param_spec_double ("rel-error",
                                                        NULL,
                                                        "Relative error",
                                                        GSL_DBL_EPSILON, 1.0, 1.0e-8,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
/**
 * NcmSplineFuncTest:scale:
 *
 * Scale of function. It is used to compute the absolute tolerance.
 *
 */
  g_object_class_install_property (object_class,
                                   PROP_SCALE,
                                   g_param_spec_double ("scale",
                                                        NULL,
                                                        "Scale",
                                                        0.0, G_MAXDOUBLE, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_spline_func_test_new:
 *
 * Allocates memory for a new #NcmSplineFuncTest suite
 * with the default parameters values.
 *
 * Returns: a new #NcmSplineFuncTest.
 */
NcmSplineFuncTest *
ncm_spline_func_test_new (void)
{
  NcmSplineFuncTest *sft = g_object_new (NCM_TYPE_SPLINE_FUNC_TEST,
                                         NULL);
  
  return sft;
}

/**
 * ncm_spline_func_test_ref:
 * @sft: a #NcmSplineFuncTest
 *
 * Increases the reference count of @sft by one atomically.
 *
 * Returns: (transfer full): @sft.
 */
NcmSplineFuncTest *
ncm_spline_func_test_ref (NcmSplineFuncTest *sft)
{
  return g_object_ref (sft);
}

/**
 * ncm_spline_func_test_unref:
 * @sft: a #NcmSplineFuncTest
 *
 * Atomically decrements the reference count of @sft by one.
 * If the reference count drops to 0, all memory allocated by @sft is released.
 *
 */
void
ncm_spline_func_test_unref (NcmSplineFuncTest *sft)
{
  g_object_unref (sft);
}

/**
 * ncm_spline_func_test_clear:
 * @sft: a #NcmSplineFuncTest
 *
 * If @sft is different from NULL,
 * atomically decrements the reference count of @sft by one.
 * If the reference count drops to 0,
 * all memory allocated by @sft is released and @sft is set to NULL.
 *
 */
void
ncm_spline_func_test_clear (NcmSplineFuncTest **sft)
{
  g_clear_object (sft);
}

/**
 * ncm_spline_func_test_set_type:
 * @sft: a #NcmSplineFuncTest
 * @type: a #NcmSplineFuncTestType
 *
 * Sets the type method to generate the base function
 * to be used as test to the [NcmSplineFunc](numcosmo-NcmSplineFunc.html) object.
 *
 */
void
ncm_spline_func_test_set_type (NcmSplineFuncTest *sft, NcmSplineFuncTestType type)
{
  NcmSplineFuncTestPrivate * const self = sft->priv;
  
  self->type = type;
}

/**
 * ncm_spline_func_test_set_ngrid:
 * @sft: a #NcmSplineFuncTest
 * @ngrid: the number of knots in the grid used to compare the base function with [NcmSplineFunc](numcosmo-NcmSplineFunc.html) 
 *
 * Sets the number of knots in the linear grid used to compare the base function with [NcmSplineFunc](numcosmo-NcmSplineFunc.html).
 *
 * The knots are evenly distributed in the range [@xi, @xf], including both limiting points.
 *
 */
void
ncm_spline_func_test_set_ngrid (NcmSplineFuncTest *sft, const guint ngrid)
{
  NcmSplineFuncTestPrivate * const self = sft->priv;
  
  self->ngrid = ngrid;
}

/**
 * ncm_spline_func_test_set_seed:
 * @sft: a #NcmSplineFuncTest
 * @seed: the seed for the random number generator
 *
 * Sets the seed for the random number generator. Must be a positive integer.
 * If it is set to zero, it generates a seed internally.
 *
 * See also #ncm_rng_set_seed and #ncm_rng_set_random_seed from #NcmRNG.
 *
 */
void
ncm_spline_func_test_set_seed (NcmSplineFuncTest *sft, const gulong seed)
{
  NcmSplineFuncTestPrivate * const self = sft->priv;
  
  self->seed = seed;
}

/**
 * ncm_spline_func_test_set_params_info:
 * @sft: a #NcmSplineFuncTest
 * @par_info: a #NcmMatrix
 *
 * Sets the matrix @par_info in order to create the base function and to perform statistics.
 *
 * The matrix must have dimensions (rows = number of parameters, cols = 2).
 * - columns: (min, max) for a flat PDF or (mean, sigma) for a normal PDF.
 *
 */
void
ncm_spline_func_test_set_params_info (NcmSplineFuncTest *sft, NcmMatrix *par_info)
{
  NcmSplineFuncTestPrivate * const self = sft->priv;
  
  g_assert_nonnull (par_info);
  g_assert_cmpuint (ncm_matrix_nrows (par_info),  >, 0);
  g_assert_cmpuint (ncm_matrix_ncols (par_info), ==, 2);
  
  self->par_info = ncm_matrix_ref (par_info);
}

/**
 * ncm_spline_func_test_set_params_info_all:
 * @sft: a #NcmSplineFuncTest
 * @npar: number of parameters
 * @p1: parameter to be used by the PDF
 * @p2: parameter to be used by the PDF
 *
 * Sets all values of the matrix @par_info to @p1 and @p2
 * in order to create the base function and to perform statistics.
 *
 * The matrix must have dimensions (rows = @npar, cols = 2).
 * - columns: (@p1 = min. value , @p2 = max. value) for a flat PDF or (@p1 = mean, @p2 = sigma) for a normal PDF.
 *
 */
void
ncm_spline_func_test_set_params_info_all (NcmSplineFuncTest *sft, const guint npar, const gdouble p1, const gdouble p2)
{
  NcmSplineFuncTestPrivate * const self = sft->priv;
  
  guint i;
  
  self->par_info = ncm_matrix_new (npar, 2);
  
  for (i = 0; i < npar; i++)
  {
    ncm_matrix_set (self->par_info, i, 0, p1);
    ncm_matrix_set (self->par_info, i, 1, p2);
  }
}

/**
 * ncm_spline_func_test_set_xi:
 * @sft: a #NcmSplineFuncTest
 * @xi: initial abscissa value
 *
 * Sets the initial abscissa value to @xi.
 *
 */
void
ncm_spline_func_test_set_xi (NcmSplineFuncTest *sft, const gdouble xi)
{
  NcmSplineFuncTestPrivate * const self = sft->priv;
  
  self->xi = xi;
}

/**
 * ncm_spline_func_test_set_xf:
 * @sft: a #NcmSplineFuncTest
 * @xf: final abscissa value
 *
 * Sets the final abscissa value to @xf.
 *
 */
void
ncm_spline_func_test_set_xf (NcmSplineFuncTest *sft, const gdouble xf)
{
  NcmSplineFuncTestPrivate * const self = sft->priv;
  
  self->xf = xf;
}

/**
 * ncm_spline_func_test_set_rel_error:
 * @sft: a #NcmSplineFuncTest
 * @rel_error: relative error
 *
 * Sets the relative error (tolerance) used by [NcmSplineFunc](numcosmo-NcmSplineFunc.html) object when defining its knots.
 *
 * See also #ncm_spline_set_func and #ncm_spline_set_func_scale from [NcmSplineFunc](numcosmo-NcmSplineFunc.html).
 *
 */
void
ncm_spline_func_test_set_rel_error (NcmSplineFuncTest *sft, const gdouble rel_error)
{
  NcmSplineFuncTestPrivate * const self = sft->priv;
  
  self->rel_error = rel_error;
}

/**
 * ncm_spline_func_test_set_scale:
 * @sft: a #NcmSplineFuncTest
 * @scale: scale of function
 *
 * Sets the @scale. It is used to compute the absolute tolerance.
 * See Eq. \eqref{eq:condition}. 
 *
 * See also #ncm_spline_set_func_scale from [NcmSplineFunc](numcosmo-NcmSplineFunc.html).
 *
 */
void
ncm_spline_func_test_set_scale (NcmSplineFuncTest *sft, const gdouble scale)
{
  NcmSplineFuncTestPrivate * const self = sft->priv;
  
  self->scale = scale;
}

/**
 * ncm_spline_func_test_set_out_threshold:
 * @sft: a #NcmSplineFuncTest
 * @out_threshold: threshold of outliers in a grid (\%)
 *
 * Sets the outliers threshold in a grid in order to save the functions informations 
 * for further analysis.ncm_spline_func_test_save_grid_functions_to_txt() 
 * and ncm_spline_func_test_save_knots_to_txt().
 * The files names are "functions_with_outlier_above_threshold.*" and "knots_with_outlier_above_threshold.*". 
 * Where * stands for the realization number.
 *
 * If the #NcmSplineFuncTestType is an interpolation, then saves one more file with
 * the interpolated points with the function , "interpolated_points_with_outlier_above_threshold.*".
 *
 * The threshold is given in percentage of @ngrid.
 *
 */
void
ncm_spline_func_test_set_out_threshold (NcmSplineFuncTest *sft, const gdouble out_threshold)
{
  NcmSplineFuncTestPrivate * const self = sft->priv;

  self->out_threshold = out_threshold;
}

/**
 * ncm_spline_func_test_set_user_gsl_function:
 * @sft: a #NcmSplineFuncTest
 * @F: a pointer to a GSL function
 *
 * Sets the user supplied [GSL](https://www.gnu.org/software/gsl/) function as base for the test.
 *
 * See also [gsl_function](https://www.gnu.org/software/gsl/doc/html/roots.html?highlight=gsl_function#c.gsl_function) for more information.
 *
 */
void
ncm_spline_func_test_set_user_gsl_function (NcmSplineFuncTest *sft, gsl_function *F)
{
  NcmSplineFuncTestPrivate * const self = sft->priv;
  
  g_assert_cmpuint (self->type, ==, NCM_SPLINE_FUNC_TEST_TYPE_USER);
  
  self->F.function = F->function;
  self->F.params   = F->params;
}

/**
 * ncm_spline_func_test_get_ngrid:
 * @sft: a #NcmSplineFuncTest
 *
 * Returns: the number of knots in the range [@xi, @xf], including both limiting points.
 */
guint
ncm_spline_func_test_get_ngrid (NcmSplineFuncTest *sft)
{
  return sft->priv->ngrid;
}

/**
 * ncm_spline_func_test_get_seed:
 * @sft: a #NcmSplineFuncTest
 *
 * Returns: the random number generator @seed.
 */
gulong
ncm_spline_func_test_get_seed (NcmSplineFuncTest *sft)
{
  return sft->priv->seed;
}

/**
 * ncm_spline_func_test_get_params_info:
 * @sft: a #NcmSplineFuncTest
 *
 * Returns: (transfer full): a #NcmMatrix with parameters information.
 */
NcmMatrix *
ncm_spline_func_test_get_params_info (NcmSplineFuncTest *sft)
{
  return ncm_matrix_ref (sft->priv->par_info);
}

/**
 * ncm_spline_func_test_peek_current_params:
 * @sft: a #NcmSplineFuncTest
 *
 * A #NcmVector with the currrent parameters value in the same order
 * as in the user provided in #NcmMatrix @par_info.
 *
 * This function is necessary when #NCM_SPLINE_FUNC_TEST_TYPE_USER is applied.
 *
 * Returns: (transfer none): a pointer to the currrent parameters vector.
 */
NcmVector *
ncm_spline_func_test_peek_current_params (NcmSplineFuncTest *sft)
{
  return sft->priv->params;
}

/**
 * ncm_spline_func_test_get_xi:
 * @sft: a #NcmSplineFuncTest
 *
 * Returns: the initial abscissa value @xi.
 */
gdouble
ncm_spline_func_test_get_xi (NcmSplineFuncTest *sft)
{
  return sft->priv->xi;
}

/**
 * ncm_spline_func_test_get_xf:
 * @sft: a #NcmSplineFuncTest
 *
 * Returns: the final abscissa value @xf.
 */
gdouble
ncm_spline_func_test_get_xf (NcmSplineFuncTest *sft)
{
  return sft->priv->xf;
}

/**
 * ncm_spline_func_test_get_rel_error:
 * @sft: a #NcmSplineFuncTest
 *
 * Returns: the relative error @rel_error.
 */
gdouble
ncm_spline_func_test_get_rel_error (NcmSplineFuncTest *sft)
{
  return sft->priv->rel_error;
}

/**
 * ncm_spline_func_test_get_scale:
 * @sft: a #NcmSplineFuncTest
 *
 * Returns: the @scale.
 */
gdouble
ncm_spline_func_test_get_scale (NcmSplineFuncTest *sft)
{
  return sft->priv->scale;
}

/**
 * ncm_spline_func_test_get_out_threshold:
 * @sft: a #NcmSplineFuncTest
 *
 * Returns: @out_threshold.
 */
gdouble
ncm_spline_func_test_get_out_threshold (NcmSplineFuncTest *sft)
{
  return sft->priv->out_threshold;
}

static void _ncm_spline_func_test_drawn_params (NcmSplineFuncTestPrivate *self);

static void _ncm_spline_func_test_set_params_is_fixed_normal (NcmSplineFuncTestPrivate *self, guint i, const gdouble p1, const gdouble p2);
static void _ncm_spline_func_test_set_params_is_fixed_flat (NcmSplineFuncTestPrivate *self, guint i, const gdouble p1, const gdouble p2);
static void _ncm_spline_func_test_set_params_is_fixed (NcmSplineFuncTestPrivate *self);
static void _ncm_spline_func_test_set_grid (NcmSplineFuncTestPrivate *self);
static void _ncm_spline_func_test_set_rng (NcmSplineFuncTestPrivate *self);

static void _ncm_spline_func_test_prepare_spl (NcmSplineFuncTestPrivate *self);
static void _ncm_spline_func_test_prepare_to_loop (NcmSplineFuncTestPrivate *self);

static void _ncm_spline_func_test_save_interpolated_points (NcmSplineFuncTestPrivate *self, gchar *fname);

static gdouble _ncm_spline_func_test_gsl_eval_cos (gdouble x, gpointer p);
static gdouble _ncm_spline_func_test_gsl_eval_spl (gdouble x, gpointer p);

/**
 * ncm_spline_func_test_prepare:
 * @sft: a #NcmSplineFuncTest
 * @ftype: a #NcmSplineFuncType
 * @pdftype: a #NcmSplineFuncTestTypePDF
 *
 * Prepares #NcmSplineFuncTest suite in order to evaluate
 * one grid statistics, ncm_spline_func_test_set_one_grid_stats(),
 *  and also Monte Carlo statistics, ncm_spline_func_test_monte_carlo() and ncm_spline_func_test_monte_carlo_and_save_to_txt().
 */
void
ncm_spline_func_test_prepare (NcmSplineFuncTest *sft, NcmSplineFuncType ftype, NcmSplineFuncTestTypePDF pdftype)
{
  NcmSplineFuncTestPrivate *self = sft->priv;
  
  self->pdftype = pdftype;
  self->ftype   = ftype;
  self->npar    = ncm_matrix_nrows (self->par_info);
  self->params  = ncm_vector_new (self->npar);
  self->x_spl   = ncm_vector_new (self->npar);
  
  
  if (self->type == NCM_SPLINE_FUNC_TEST_TYPE_COSINE)
    g_assert (GSL_IS_EVEN (self->npar));
  
  _ncm_spline_func_test_set_grid (self);
  
  _ncm_spline_func_test_set_params_is_fixed (self);
  
  _ncm_spline_func_test_set_rng (self);
  
  switch (self->type)
  {
    case NCM_SPLINE_FUNC_TEST_TYPE_POLYNOMIAL:
      self->F.function = &_ncm_spline_func_test_gsl_eval_spl;
      self->F.params   = ncm_spline_gsl_new (gsl_interp_polynomial);
      break;
    case NCM_SPLINE_FUNC_TEST_TYPE_COSINE:
      self->F.function = &_ncm_spline_func_test_gsl_eval_cos;
      self->F.params   = self->params;
      break;
    case NCM_SPLINE_FUNC_TEST_TYPE_RBF:
      self->F.function = &_ncm_spline_func_test_gsl_eval_spl;
      self->F.params   = ncm_spline_rbf_new (NCM_SPLINE_RBF_TYPE_GAUSS);
      break;
    case NCM_SPLINE_FUNC_TEST_TYPE_USER:
      g_assert_nonnull (self->F.function);
      break;
    default:
      printf ("This option is not available yet.\n");
      g_assert_not_reached ();
      break;
  }
}

/**
 * ncm_spline_func_test_set_one_grid_stats:
 * @sft: a #NcmSplineFuncTest
 *
 * Sets the #NcmSplineFuncTest suite in order to perform one grid statistics.
 *
 * <note>
 *   <para>
 *    It must be called after ncm_spline_func_test_prepare().
 *   </para>
 * </note>
 */
void
ncm_spline_func_test_set_one_grid_stats (NcmSplineFuncTest *sft)
{
  NcmSplineFuncTestPrivate *self = sft->priv;
  
  guint i;
  
  self->ncm_out = 0;
  self->lin_out = 0;
  
  ncm_stats_vec_reset (self->stats_grid, TRUE);
  
  _ncm_spline_func_test_prepare_to_loop (self);
  
  for (i = 0; i < self->ngrid; i++)
  {
    const gdouble x     = ncm_vector_fast_get (self->xgrid, i);
    const gdouble ytrue = GSL_FN_EVAL (&self->F, x);
    const gdouble yncm  = ncm_spline_eval (self->ncm, x);
    const gdouble ylin  = ncm_spline_eval (self->lin, x);
    
    const gdouble diff_ncm  = ytrue - yncm;
    const gdouble diff_lin  = ytrue - ylin;
    const gdouble threshold = (fabs (ytrue) + self->scale) * self->rel_error;
    
    if (fabs (diff_ncm) > threshold)
      self->ncm_out++;
    
    if (fabs (diff_lin) > threshold)
      self->lin_out++;
    
    ncm_vector_fast_set (self->ygrid_true, i, ytrue);
    ncm_vector_fast_set (self->ygrid_ncm,  i, yncm);
    ncm_vector_fast_set (self->ygrid_lin,  i, ylin);
    ncm_vector_fast_set (self->diff_ncm,   i, diff_ncm);
    ncm_vector_fast_set (self->diff_lin,   i, diff_lin);
    
    ncm_stats_vec_set (self->stats_grid, 0, diff_ncm);
    ncm_stats_vec_set (self->stats_grid, 1, diff_lin);
    
    ncm_stats_vec_update (self->stats_grid);
  }
}

/**
 * ncm_spline_func_test_monte_carlo:
 * @sft: a #NcmSplineFuncTest
 * @nsim: the number of Monte Carlo simulations to be performed
 *
 * Performs a Monte Carlo simulation on the base function.
 *
 * See #NcmSplineFuncTestType and #NcmSplineFuncTestTypePDF for the available options.
 *
 * <note>
 *   <para>
 *     It must be called after ncm_spline_func_test_prepare().
 *   </para>
 * </note>
 */
void
ncm_spline_func_test_monte_carlo (NcmSplineFuncTest *sft, guint nsim)
{
  NcmSplineFuncTestPrivate *self = sft->priv;

  guint ncount = 0;

  if (self->all_fixed)
    nsim = 1;

  self->nsim = nsim;

  self->stats_mc = ncm_stats_vec_new (9, NCM_STATS_VEC_VAR, FALSE);

  do {
    gdouble ncm_min, ncm_max, lin_min, lin_max, ncm_out, lin_out;

    ncm_spline_func_test_set_one_grid_stats (sft);

    if (self->ncm_out == 0)
      self->ncm_good++;

    if (self->lin_out == 0)
      self->lin_good++;

    ncm_out = (100. * self->ncm_out) / (1.0 * self->ngrid);
    lin_out = (100. * self->lin_out) / (1.0 * self->ngrid);

    if ((ncm_out > self->out_threshold) && (self->out_threshold > G_MINDOUBLE))
    {
      gchar *f1 = g_strdup_printf ("%s.%u", "functions_with_outlier_above_threshold", ncount + 1);
      gchar *f2 = g_strdup_printf ("%s.%u", "knots_with_outlier_above_threshold", ncount + 1);

      ncm_spline_func_test_save_grid_functions_to_txt (sft, f1);
      ncm_spline_func_test_save_knots_to_txt (sft, f2);

      g_free (f1);
      g_free (f2);

      if ((self->type == NCM_SPLINE_FUNC_TEST_TYPE_POLYNOMIAL) || self->type == NCM_SPLINE_FUNC_TEST_TYPE_RBF)
      {
        gchar *f3 = g_strdup_printf ("%s.%u", "interpolated_points_with_outlier_above_threshold", ncount + 1);
        _ncm_spline_func_test_save_interpolated_points (self, f3);
        g_free (f3);
      }
    }

    ncm_vector_get_absminmax (self->diff_ncm, &ncm_min, &ncm_max);
    ncm_vector_get_absminmax (self->diff_lin, &lin_min, &lin_max);

    ncm_stats_vec_set (self->stats_mc, 0, (1.0 * self->len));
    ncm_stats_vec_set (self->stats_mc, 1, ncm_out);
    ncm_stats_vec_set (self->stats_mc, 2, lin_out);
    ncm_stats_vec_set (self->stats_mc, 3, ncm_stats_vec_get_mean (self->stats_grid, 0));
    ncm_stats_vec_set (self->stats_mc, 4, ncm_stats_vec_get_mean (self->stats_grid, 1));
    ncm_stats_vec_set (self->stats_mc, 5, ncm_stats_vec_get_sd (self->stats_grid, 0));
    ncm_stats_vec_set (self->stats_mc, 6, ncm_stats_vec_get_sd (self->stats_grid, 1));
    ncm_stats_vec_set (self->stats_mc, 7, ncm_max);
    ncm_stats_vec_set (self->stats_mc, 8, lin_max);

    ncm_stats_vec_update (self->stats_mc);

    ncount++;
  } while (ncount < nsim);
}

/**
 * ncm_spline_func_test_monte_carlo_and_save_to_txt:
 * @sft: a #NcmSplineFuncTest
 * @nsim: the number of Monte Carlo simulations to be performed
 * @fname: the name of the text file to save mc statistics
 *
 * Performs a Monte Carlo simulation on the base function
 * and saves some statistics in the @fname file.
 *
 * Same as ncm_spline_func_test_monte_carlo() but creating a ascii file with the statistics on the fly. 
 *
 * See #NcmSplineFuncTestType and #NcmSplineFuncTestTypePDF for the available options.
 *
 * <note>
 *   <para>
 *     It must be called after ncm_spline_func_test_prepare().
 *   </para>
 * </note>
 */
void
ncm_spline_func_test_monte_carlo_and_save_to_txt (NcmSplineFuncTest *sft, guint nsim, gchar *fname)
{
  NcmSplineFuncTestPrivate *self = sft->priv;
  
  FILE *fout = fopen (fname, "w");
  
  fprintf (fout, "#### len, ncm_out (%%), lin_out (%%), ncm mean diff., lin mean diff., ncm stdev. diff., lin stdev. diff., ncm max. diff., lin max. diff. ####\n");
  
  guint ncount = 0;
  
  if (self->all_fixed)
    nsim = 1;
  
  self->nsim = nsim;
  
  self->stats_mc = ncm_stats_vec_new (9, NCM_STATS_VEC_VAR, FALSE);
  
  do {
    gdouble ncm_min, ncm_max, lin_min, lin_max, ncm_out, lin_out;
    
    ncm_spline_func_test_set_one_grid_stats (sft);

    if (self->ncm_out == 0)
      self->ncm_good++;

    if (self->lin_out == 0)
      self->lin_good++;

    ncm_out = (100. * self->ncm_out) / (1.0 * self->ngrid);
    lin_out = (100. * self->lin_out) / (1.0 * self->ngrid);

    if ((ncm_out > self->out_threshold) && (self->out_threshold > G_MINDOUBLE))
    {
      gchar *f1 = g_strdup_printf ("%s.%u", "functions_with_outlier_above_threshold", ncount + 1);
      gchar *f2 = g_strdup_printf ("%s.%u", "knots_with_outlier_above_threshold", ncount + 1);

      ncm_spline_func_test_save_grid_functions_to_txt (sft, f1);
      ncm_spline_func_test_save_knots_to_txt (sft, f2);

      g_free (f1);
      g_free (f2);

      if ((self->type == NCM_SPLINE_FUNC_TEST_TYPE_POLYNOMIAL) || self->type == NCM_SPLINE_FUNC_TEST_TYPE_RBF)
      {
        gchar *f3 = g_strdup_printf ("%s.%u", "interpolated_points_with_outlier_above_threshold", ncount + 1);
        _ncm_spline_func_test_save_interpolated_points (self, f3);
        g_free (f3);
      }
    }

    ncm_vector_get_absminmax (self->diff_ncm, &ncm_min, &ncm_max);
    ncm_vector_get_absminmax (self->diff_lin, &lin_min, &lin_max);
    
    fprintf (fout, "%6u  %7.3e  %7.3e  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e  %22.15e\n", 
             self->len,
             (100. * self->ncm_out) / (1. * self->ngrid),
             (100. * self->lin_out) / (1. * self->ngrid),
             ncm_stats_vec_get_mean (self->stats_grid, 0),
             ncm_stats_vec_get_mean (self->stats_grid, 1),
             ncm_stats_vec_get_sd (self->stats_grid, 0),
             ncm_stats_vec_get_sd (self->stats_grid, 1),
             ncm_max, lin_max);
    
    ncm_stats_vec_set (self->stats_mc, 0, (1.0 * self->len));
    ncm_stats_vec_set (self->stats_mc, 1, ncm_out);
    ncm_stats_vec_set (self->stats_mc, 2, lin_out);
    ncm_stats_vec_set (self->stats_mc, 3, ncm_stats_vec_get_mean (self->stats_grid, 0));
    ncm_stats_vec_set (self->stats_mc, 4, ncm_stats_vec_get_mean (self->stats_grid, 1));
    ncm_stats_vec_set (self->stats_mc, 5, ncm_stats_vec_get_sd (self->stats_grid, 0));
    ncm_stats_vec_set (self->stats_mc, 6, ncm_stats_vec_get_sd (self->stats_grid, 1));
    ncm_stats_vec_set (self->stats_mc, 7, ncm_max);
    ncm_stats_vec_set (self->stats_mc, 8, lin_max);
    
    ncm_stats_vec_update (self->stats_mc);
    
    ncount++;
  } while (ncount < nsim);
  
  fclose (fout);
}

/**
 * ncm_spline_func_test_log_vals_one_grid_stats:
 * @sft: a #NcmSplineFuncTest
 *
 * Prints the current grid statistics. 
 *
 * See [description][NcmSplineFuncTest.description] for information about the parameters.
 *
 * <note>
 *   <para>
 *    It must be called after ncm_spline_func_test_set_one_grid_stats(). 
 *   </para>
 * </note>
 */
void
ncm_spline_func_test_log_vals_one_grid_stats (NcmSplineFuncTest *sft)
{
  NcmSplineFuncTestPrivate *self = sft->priv;
  
  gdouble ncm_min, ncm_max, lin_min, lin_max;
  
  ncm_vector_get_absminmax (self->diff_ncm, &ncm_min, &ncm_max);
  ncm_vector_get_absminmax (self->diff_lin, &lin_min, &lin_max);
  
  printf ("\n##### Grid statistics  #####\n\n");
  
  printf (" * NcmSplineFunc number of knots = %u\n\n", self->len);
  
  printf (" * (ncm) diff. = %8.3e +/- %8.3e (abs. max. diff. = %8.3e) | outliers = %5.2f %%\n", 
          ncm_stats_vec_get_mean (self->stats_grid, 0),
          ncm_stats_vec_get_sd (self->stats_grid, 0),
          ncm_max,
          (100. * self->ncm_out) / (1. * self->ngrid));
  
  printf (" * (lin) diff. = %8.3e +/- %8.3e (abs. max. diff. = %8.3e) | outliers = %5.2f %%\n\n", 
          ncm_stats_vec_get_mean (self->stats_grid, 1),
          ncm_stats_vec_get_sd (self->stats_grid, 1),
          lin_max,
          (100. * self->lin_out) / (1. * self->ngrid));
}

/**
 * ncm_spline_func_test_save_grid_functions_to_txt:
 * @sft: a #NcmSplineFuncTest
 * @fname: text file name to save one grid informations
 *
 * Saves one grid functions in the text @fname file. The colums are:
 * - $x$.
 * - the base function $f(x)$.
 * - the [NcmSplineFunc](numcosmo-NcmSplineFunc.html) estimation of $f(x)$.
 * - the linear grid estimation of $f(x)$ with the same number of knots as the [NcmSplineFunc](numcosmo-NcmSplineFunc.html).
 *
 * <note>
 *   <para>
 *    It must be called after ncm_spline_func_test_set_one_grid_stats(). 
 *   </para>
 * </note>
 */
void
ncm_spline_func_test_save_grid_functions_to_txt (NcmSplineFuncTest *sft, gchar *fname)
{
  NcmSplineFuncTestPrivate *self = sft->priv;
  
  guint i;
  
  FILE *fout = fopen (fname, "w");
  
  fprintf (fout, "##### x , true,  ncm,  linear  #####\n");
  
  for (i = 0; i < self->ngrid; i++)
  {
    const gdouble x    = ncm_vector_fast_get (self->xgrid, i);
    //const gdouble dfdx = ncm_spline_eval_deriv2 (NCM_SPLINE (self->F.params), x); 

    //fprintf (fout, "%22.15e  %22.15e  %22.15e  %22.15e  %22.15e\n", 
    fprintf (fout, "%22.15e  %22.15e  %22.15e  %22.15e\n", 
             x,
             ncm_vector_fast_get (self->ygrid_true, i),
             //dfdx,
             ncm_vector_fast_get (self->ygrid_ncm, i),
             ncm_vector_fast_get (self->ygrid_lin, i));
  }
  
  fclose (fout);
}

/**
 * ncm_spline_func_test_save_knots_to_txt:
 * @sft: a #NcmSplineFuncTest
 * @fname: text file name to save the [NcmSplineFunc](numcosmo-NcmSplineFunc.html) knots
 *
 * Saves the knots defined by [NcmSplineFunc](numcosmo-NcmSplineFunc.html) in the text file @fname (only the current grid).
 *
 * <note>
 *   <para>
 *    It must be called after ncm_spline_func_test_set_one_grid_stats(). 
 *   </para>
 * </note>
 */
void
ncm_spline_func_test_save_knots_to_txt (NcmSplineFuncTest *sft, gchar *fname)
{
  NcmSplineFuncTestPrivate * const self = sft->priv;
  
  guint i;
  
  NcmVector *knots = ncm_spline_get_xv (self->ncm);
  
  FILE *fout = fopen (fname, "w");
  
  for (i = 0; i < ncm_spline_get_len (self->ncm); i++)
    fprintf (fout, "%22.15e\n", ncm_vector_fast_get (knots, i));
  
  fclose (fout);
  
  ncm_vector_free (knots);
}

/**
 * ncm_spline_func_test_log_vals_mc_stats:
 * @sft: a #NcmSplineFuncTest
 *
 * Prints the Monte Carlo statistics. 
 *
 * See [description][NcmSplineFuncTest.description] for information about the parameters.
 *
 * <note>
 *   <para>
 *    It must be called after ncm_spline_func_test_monte_carlo() or ncm_spline_func_test_monte_carlo_and_save_to_txt(). 
 *   </para>
 * </note>
 */
void
ncm_spline_func_test_log_vals_mc_stats (NcmSplineFuncTest *sft)
{
  NcmSplineFuncTestPrivate * const self = sft->priv;
  
  printf ("\n####  Monte Carlo statistics for %u simulations  ####\n\n", self->nsim);
  printf ("   * NcmSplineFunc number of knots: %8.2f +/- %8.2f\n", ncm_stats_vec_get_mean (self->stats_mc, 0), ncm_stats_vec_get_sd (self->stats_mc, 0));
  printf ("\n");
  printf ("   * Ncm clean grid (%%): %5.2f\n", 100.*self->ncm_good / (1.*self->nsim));
  printf ("   * Lin clean grid (%%): %5.2f\n", 100.*self->lin_good / (1.*self->nsim));
  printf ("\n");
  printf ("   * Ncm outliers (%%): %5.2f +/- %5.2f\n", ncm_stats_vec_get_mean (self->stats_mc, 1), ncm_stats_vec_get_sd (self->stats_mc, 1));
  printf ("   * Lin outliers (%%): %5.2f +/- %5.2f\n", ncm_stats_vec_get_mean (self->stats_mc, 2), ncm_stats_vec_get_sd (self->stats_mc, 2));
  printf ("\n");
  printf ("   * Ncm diff. : %7.3e +/- %7.3e\n", ncm_stats_vec_get_mean (self->stats_mc, 3), ncm_stats_vec_get_mean (self->stats_mc, 5));
  printf ("   * Lin diff. : %7.3e +/- %7.3e\n", ncm_stats_vec_get_mean (self->stats_mc, 4), ncm_stats_vec_get_mean (self->stats_mc, 6));
  printf ("\n");
  printf ("   * Ncm abs. max. diff. : %7.3e +/- %7.3e\n", ncm_stats_vec_get_mean (self->stats_mc, 7), ncm_stats_vec_get_sd (self->stats_mc, 7));
  printf ("   * Lin abs. max. diff. : %7.3e +/- %7.3e\n", ncm_stats_vec_get_mean (self->stats_mc, 8), ncm_stats_vec_get_sd (self->stats_mc, 8));
  printf ("\n");
}

static void
_ncm_spline_func_test_drawn_params (NcmSplineFuncTestPrivate *self)
{
  guint i;
  
  for (i = 0; i < self->npar; i++)
  {
    const gboolean fixed = g_array_index (self->par_is_fixed, gboolean, i);
    
    if (fixed)
    {
      ncm_vector_fast_set (self->params, i, ncm_matrix_get (self->par_info, i, 0));
    }
    else
    {
      const gdouble p1 = ncm_matrix_get (self->par_info, i, 0);
      const gdouble p2 = ncm_matrix_get (self->par_info, i, 1);
      
      switch (self->pdftype)
      {
        case NCM_SPLINE_FUNC_TEST_TYPE_PDF_FLAT:
          ncm_vector_fast_set (self->params, i, ncm_rng_uniform_gen (self->rng, p1, p2));
          break;
        case NCM_SPLINE_FUNC_TEST_TYPE_PDF_NORMAL:
          ncm_vector_fast_set (self->params, i, ncm_rng_gaussian_gen (self->rng, p1, p2));
          break;
        default:
          printf ("This PDF option is not available yet.\n");
          g_assert_not_reached ();
          break;
      }
    }
  }
}

static void
_ncm_spline_func_test_set_grid (NcmSplineFuncTestPrivate *self)
{
  guint i;
  
  const gdouble dx = (self->xf - self->xi) / (self->ngrid - 1.0);
  
  gint cmp = ncm_cmp (self->xi, self->xf, 1.e-12, 0.0);
  
  ncm_assert_cmpdouble_e (self->xf, !=, self->xi, 1.e-12, 0.0);
  
  if (cmp == 1)
  {
    gdouble xf = self->xi;
    
    self->xi = self->xf;
    self->xf = xf;
  }
  
  self->xgrid      = ncm_vector_new (self->ngrid);
  self->ygrid_true = ncm_vector_new (self->ngrid);
  self->ygrid_ncm  = ncm_vector_new (self->ngrid);
  self->ygrid_lin  = ncm_vector_new (self->ngrid);
  self->diff_ncm   = ncm_vector_new (self->ngrid);
  self->diff_lin   = ncm_vector_new (self->ngrid);
  
  self->stats_grid = ncm_stats_vec_new (2, NCM_STATS_VEC_VAR, FALSE);
  
  for (i = 0; i < self->ngrid; i++)
    ncm_vector_fast_set (self->xgrid, i, (dx * i + self->xi));
}

static void
_ncm_spline_func_test_set_params_is_fixed_normal (NcmSplineFuncTestPrivate *self, guint i, const gdouble p1, const gdouble p2)
{
  const gboolean f = FALSE;
  const gboolean t = TRUE;
  
  if (fabs (p2) > G_MINDOUBLE)
  {
    g_array_append_val (self->par_is_fixed, f);
    self->all_fixed = FALSE;
    
    ncm_matrix_set (self->par_info, i, 0, p1);
    ncm_matrix_set (self->par_info, i, 1, fabs (p2));
  }
  else
  {
    g_array_append_val (self->par_is_fixed, t);
  }
}

static void
_ncm_spline_func_test_set_params_is_fixed_flat (NcmSplineFuncTestPrivate *self, guint i, const gdouble p1, const gdouble p2)
{
  const gboolean f = FALSE;
  const gboolean t = TRUE;
  
  if (ncm_cmp (p1, p2, 1.e-12, 0.0))
  {
    g_array_append_val (self->par_is_fixed, f);
    self->all_fixed = FALSE;
    
    ncm_matrix_set (self->par_info, i, 0, GSL_MIN (p1, p2));
    ncm_matrix_set (self->par_info, i, 1, GSL_MAX (p1, p2));
  }
  else
  {
    g_array_append_val (self->par_is_fixed, t);
  }
}

static void
_ncm_spline_func_test_set_params_is_fixed (NcmSplineFuncTestPrivate *self)
{
  guint i;
  
  for (i = 0; i < self->npar; i++)
  {
    const gdouble p1 = ncm_matrix_get (self->par_info, i, 0);
    const gdouble p2 = ncm_matrix_get (self->par_info, i, 1);
    
    switch (self->pdftype)
    {
      case NCM_SPLINE_FUNC_TEST_TYPE_PDF_FLAT:
        _ncm_spline_func_test_set_params_is_fixed_flat (self, i, p1, p2);
        break;
      case NCM_SPLINE_FUNC_TEST_TYPE_PDF_NORMAL:
        _ncm_spline_func_test_set_params_is_fixed_normal (self, i, p1, p2);
        break;
      default:
        printf ("This PDF option is not available yet.\n");
        g_assert_not_reached ();
        break;
    }
  }
}

static void
_ncm_spline_func_test_set_rng (NcmSplineFuncTestPrivate *self)
{
  if (self->seed == 0)
    ncm_rng_set_random_seed (self->rng, FALSE);
  else
    ncm_rng_set_seed (self->rng, self->seed);
}

static gdouble
_ncm_spline_func_test_gsl_eval_cos (gdouble x, gpointer p)
{
  guint i;
  
  const gdouble twopi = 2.0 * M_PI;
  gdouble val         = 0.0;
  NcmVector *v        = NCM_VECTOR (p);
  
  for (i = 0; i < ncm_vector_len (v); i += 2)
  {
    const gdouble amp = ncm_vector_fast_get (v,     i);
    const gdouble f   = ncm_vector_fast_get (v, i + 1);
    
    val += amp * cos (twopi * f * x);
  }
  
  return val;
}

static gdouble
_ncm_spline_func_test_gsl_eval_spl (gdouble x, gpointer p)
{
  return ncm_spline_eval (NCM_SPLINE (p), x);
}

static void
_ncm_spline_func_test_prepare_spl (NcmSplineFuncTestPrivate *self)
{
  guint i;
  
  for (i = 0; i < self->npar; i++)
    ncm_vector_fast_set (self->x_spl, i, ncm_rng_uniform_gen (self->rng, self->xi, self->xf));
  
  gsl_sort_vector (ncm_vector_gsl (self->x_spl));

  ncm_vector_fast_set (self->x_spl,              0, self->xi);
  ncm_vector_fast_set (self->x_spl, self->npar - 1, self->xf);
  
  if (self->all_fixed)
  {
    NcmVector *v = ncm_matrix_get_col (self->par_info, 0);

    self->params = ncm_vector_dup (v);
    ncm_vector_free (v);
  }
  else
  {
    _ncm_spline_func_test_drawn_params (self);
  }
  
  ncm_spline_set (NCM_SPLINE (self->F.params), self->x_spl, self->params, FALSE);
  
  ncm_spline_prepare (NCM_SPLINE (self->F.params));
}

static void
_ncm_spline_func_test_prepare_to_loop (NcmSplineFuncTestPrivate *self)
{
  switch (self->type)
  {
    case NCM_SPLINE_FUNC_TEST_TYPE_POLYNOMIAL:
      _ncm_spline_func_test_prepare_spl (self);
      break;
    case NCM_SPLINE_FUNC_TEST_TYPE_COSINE:
      _ncm_spline_func_test_drawn_params (self);
      break;
    case NCM_SPLINE_FUNC_TEST_TYPE_RBF:
      _ncm_spline_func_test_prepare_spl (self);
      break;
    case NCM_SPLINE_FUNC_TEST_TYPE_USER:
      _ncm_spline_func_test_drawn_params (self);
      break;
    default:
      printf ("This option is not available yet.\n");
      g_assert_not_reached ();
      break;
  }
  
  ncm_spline_set_func_scale (self->ncm, self->ftype, &self->F, self->xi, self->xf, 5000000, self->rel_error, self->scale);
  
  self->len = ncm_spline_get_len (self->ncm);
  
  ncm_spline_set_func_grid (self->lin, NCM_SPLINE_FUNC_GRID_LINEAR, &self->F, self->xi, self->xf, self->len);
}

static void
_ncm_spline_func_test_save_interpolated_points (NcmSplineFuncTestPrivate *self, gchar *fname)
{
  guint i;

  FILE *fout = fopen (fname, "w");

  for (i = 0; i < self->npar; i++)
    fprintf (fout, "%22.15e  %22.15e\n", ncm_vector_fast_get (self->x_spl, i), 
                                         ncm_vector_fast_get (self->params, i));

  fclose (fout);
}
