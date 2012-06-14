#!/usr/bin/env python

import sys
import configparser
import os ###The os module provides dozens of functions for interacting with the operating system
import subprocess
import matplotlib
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
import datetime
import shutil
import time

params_ids = ['clusterid', 'Fit_H0', 'Prior_H0', 'H0', 'plane', 'Omega_k', 'Fit_Omega_c', 'Omega_c', 'Fit_Omega_x', 'Omega_x',
              'Fit_Omega_b', 'Omega_b', 'Fit_ns', 'n_s_val', 'Fit_sigma_8', 'Sigma_8', 'DE_model', 'Fit_w0', 'eqos_w0', 
              'Fit_w1', 'eqos_w1', 'multiplicity', 'mf_ds_index', 'area', 'z_initial', 'z_final', 'Mi', 'Mf', 'photoz', 
              'photoz_sigma', 'Mobs', 'Mobs_sigma', 'estimator', 'zbin', 'statist_method', 'cr_x_axis', 'cr_y_axis', 
              'darkenergy', 'massf_print', 'fiducial_available']

def build_command(run_params_m):
    # building command line
    comando = []
    run_params = {}
    # put the full path for darkenergy while it is not in tawala. ask to be tawalaed.
    for param in run_params_m:
        run_params[param] = str(run_params_m[param])

    comando += [run_params['darkenergy']]
    comando += ['--cluster-id']
    comando += [run_params['clusterid']]
    comando += ['--catalog']
    comando += [run_params['halo_cat']]
    if ([run_params['plane']] == '0'):
        comando += ['--Omega_k']
        comando += [run_params['Omega_k']]
    else:
        comando += [str('--plane')]    
    comando += ['--H0']
    comando += [run_params['H0']]
    comando += ['--Omega_c']
    comando += [run_params['Omega_c']]
    comando += ['--Omega_x']
    comando += [run_params['Omega_x']]
    comando += ['--Omega_b']
    comando += [run_params['Omega_b']]
    comando += ['--n_s']
    comando += [run_params['n_s_val']]
    comando += ['--sigma_8']
    comando += [run_params['Sigma_8']]
    comando += ['--model']
    comando += [run_params['DE_model']]
    fit_params_str = 'H_0='+run_params['Fit_H0']+',\\Omega_c='+run_params['Fit_Omega_c']+',\\Omega_x='+run_params['Fit_Omega_x']+',\\Omega_b='+run_params['Fit_Omega_b']+',n_s='+run_params['Fit_ns']+',\\sigma_8='+run_params['Fit_sigma_8']
    if (run_params['DE_model'] == '0'):
        comando += ['--omega_0'] + ['-1.0']
        comando += ['--omega_1'] + ['0.0']
    elif (run_params['DE_model'] == '1'):
        comando += ['--omega_0']
        comando += [run_params['eqos_w0']]
        comando += ['--omega_1'] + ['0.0']
        fit_params_str += ',w_0='+run_params['Fit_w0']
    else:
        comando += ['--omega_0']
        comando += [run_params['eqos_w0']]
        comando += ['--omega_1']
        comando += [run_params['eqos_w1']]
        fit_params_str += ',w_0='+run_params['Fit_w0']+',w_1='+run_params['Fit_w1']
    comando += ['--fit-params']
    comando += [fit_params_str]         
    if (run_params['Prior_H0'] == '1'):
        comando += ['--H0_Hst']
    comando += ['--multiplicity']
    comando += [run_params['multiplicity']]
    comando += ['--mf_ds_index']
    comando += [run_params['mf_ds_index']]
    comando += ['--area']
    comando += [run_params['area']]
    comando += ['--z_initial']
    comando += [run_params['z_initial']]
    comando += ['--z_final']
    comando += [run_params['z_final']]
    comando += ['--Mi']
    comando += [run_params['Mi']]
    comando += ['--Mf']
    comando += [run_params['Mf']]
    if (run_params['photoz'] == '1'):
        comando += ['--photoz']
        comando += ['--photoz_sigma0']
        comando += [run_params['photoz_sigma']]
    if (run_params['Mobs'] == '1'):
        comando += ['--Mobs']
        comando += ['--lnM_sigma0']
        comando += [run_params['Mobs_sigma']]
    if (run_params['estimator'] == '0'):
        comando += ['']
    elif (run_params['estimator'] == '1'):
        comando += ['--binmass']
    else:
        comando += ['--binned']
        comando += ['--n_bins']
        comando += [run_params['zbin']]  
    if (run_params['statist_method'] == '0'):
        comando += ['--fisher']
        comando += ['--save-fisher']
    else:
        comando += ['--n-sigma'] + ['0']
    comando += ['--cr-x']
    comando += [run_params['cr_x_axis']]
    comando += ['--cr-y']
    comando += [run_params['cr_y_axis']]
    if (run_params['massf_print'] == '0'):
        comando += ['--print_mf']
    comando += ['--out']
    comando += [run_params['outFile']]
    
    return comando

#---------------------------------------------------------------------------------------------------
def plot_mass_function (run_params, outFile3):

    list_plots = {}
    list_plots['massfunc'] = []
    file = open(outFile3).read().splitlines()
    mf_vs_mass = []
    count_enter = 0
    z_bin = 0
    saved = 0
    x1 = []  
    y1 = []
    z1 = []
    u1 = []  
    for row in file:
        col = row.split()
        if not col:
            count_enter = count_enter + 1
            if count_enter == 2:
                t1 = x1
                t2 = y1
                t3 = z1
                t4 = u1
                mf_vs_mass.append([t1, t2, t3, t4])
                x1 = []
                y1 = []
                z1 = []
                u1 = []
                count_enter = 0
                saved = 1
        elif col[0].startswith('#'):
            continue
        else:
            x1.append(float(col[0]))
            y1.append(float(col[1]))
            z1.append(float(col[2]))
            u1.append(float(col[3]))
            saved = 0
    if saved == 0:
        mf_vs_mass.append([x1, y1, z1, u1])

    # plotting mass function vs. mass
    to_plot = [5, 11, 17]
    #ax = plt.subplot(111)
    for to_plot_b in to_plot:
        plot_bin = mf_vs_mass[to_plot_b]
        # plot theoretical mass function vs. mass   
        p1 = plt.plot(plot_bin[1], plot_bin[3], 'r')
        #plt.legend (p1, 'model', 'upper right')
        # plot catalog  mass function vs. mass
        p2 = plt.bar (plot_bin[1], plot_bin[2], fc='k', log=True)
        ax = plt.subplot(111)
        ax.set_yscale('log')
        ax.set_xscale('log')
        legend( (p1, p2), ('model', 'catalog'), 'best')
        xlabel = '\log_{10}[M / (h^{-1} M_{\odot})]'
        ylabel = 'dn / d\log_{10} M'
        plt.xlabel(r'$' + xlabel + '$')
        plt.ylabel('$' + ylabel + '$')
        name = 'mass_function_' + str(to_plot_b + 1) + 'zbin.png'
        #caption_filename = 'caption.txt'
        zmean = "%.2f" % float(plot_bin[0][0])
        title = 'Mean redshift: $\\bar{z} = $' + str(zmean)
        plt.title(title)
        plt.savefig(name)
        #build_caption (caption_filename, run_params, fitted_params, i, j)
        list_plots['massfunc'].append (name)
        plt.clf() 

    # plotting mass function vs. redshift for minimum mass threshold
    z = []
    mf_catalog = []
    mf_theory = []
    for i in range (len(mf_vs_mass) - 1):
        z.append(mf_vs_mass[i][0][0])
        mf_catalog.append(mf_vs_mass[i][2][0])
        mf_theory.append(mf_vs_mass[i][3][0])
        
    # plot theoretical mass function vs. redshift    
    q1 = plt.plot (z, mf_theory)
    # plot catalog mass function vs. redshift
    delta_z = ((float(run_params['z_final']) - float(run_params['z_initial'])) / 20.0) # redshift bins = 20
    # z1 list will give the initial value of each redshift bin, z list gives the mean value
    z1 = []
    for value in range(len(z)):
        z1_item = (z[value] - (delta_z/2.0))
        z1.append(z1_item)
    q2 = plt.bar (z1, mf_catalog, fc='r', alpha=0.4, width=delta_z, log=False)
    bx = plt.subplot(111)
    legend( (q1, q2), ('model', 'catalog'), 'best')
    ylabel = 'dn / d\log_{10} M'
    plt.xlabel('$ z $')
    plt.ylabel('$' + ylabel + '$')
    name_mf_z = 'mass_function_vs_z.png'
    M_min = "%.2e" % float(mf_vs_mass[0][1][0])
    title_mf_z = '$M_{min} = $' + str(M_min)
    plt.title(title_mf_z)
    plt.savefig(name_mf_z)
    list_plots['massfunc'].append (name_mf_z)
    plt.clf()

#---------------------------------------------------------------------------------------------------
##################################################################
# Function to generate the figure s caption (confidence regions) #
##################################################################
def build_caption (file_caption, run_params, fitted_params, param_x, param_y):
    #print 'chegou no caption'
    f = open (file_caption, 'w')
    
    x = 0
    y = 0
    if (isinstance(fitted_params[param_x], (list, tuple))):
        x = fitted_params[param_x][0]
        y = fitted_params[param_y][0]
    else: 
        x = fitted_params[param_x]
        y = fitted_params[param_y]

    caption = '1-3$\sigma$ confidence regions on the ($'+str(x)+'$, $'+str(y)+'$) plane obtained with '

    if (int(run_params['statist_method']) == 0):
        caption += 'Fisher Matrix.\n'
        n = len (fitted_params)
        if (n > 2):
            caption += 'Marginalized over '
            for param in range (0, n):
                if (param == param_x or param == param_y):
                    continue
                else:
                    z = fitted_params[param]
                    caption += str(z)+' '
            caption += '.\n'
            if (int(run_params['Prior_H0']) == 1):
                caption += 'It was used a Gaussian prior on $\H_0$ with mean 73.8 and $\sigma$ = 2.4.\n'
        else:
            caption += 'All other cosmological parameters were kept fixed.\n'
    else:
        caption += 'Profile Likelihood method.\n' 
        n = len (fitted_params)
        param_free = []
        i = 0
        for param in range (0, n):
            if ("FREE" == fitted_params[i][1]):
                param_free.append(fitted_params[i][0])
            i = i + 1
        m = len (param_free)
        if (m > 2):
            caption += 'Marginalized over '
            for param2 in range (0, m):
                if (param2 == param_x or param2 == param_y):
                    continue
                else:
                    z = param_free [param2]  
                    caption += '$'+str(z)+'$'
            caption += '.\n'
            if (int(run_params['Prior_H0']) == 1):
                caption += 'It was used a Gaussian prior on $\H_0$ with mean 73.8 and $\sigma$ = 2.4.\n'         
        else:
            caption += 'All other cosmological parameters were kept fixed.\n'
 
    caption += 'These constraints were obtained with halo abundance using '
    if (int(run_params['estimator']) == 0):
        caption += 'an unbinned likelihood in both mass and redshift.\n'
    elif (int(run_params['estimator']) == 1):
        caption += 'an unbinned likelihood in redshift.\n'
    else:
        caption += 'Poisson likelihood with ' + str(run_params['zbin']) + 'redshift bins.\n'

    if (int(run_params['photoz']) == 1):
        caption += 'Photo-z uncertainty was taken into account considering a Gaussian distribution with $\sigma$ = ' + str(run_params['photoz_sigma']) +'.\n'
    if (int(run_params['Mobs']) == 1):
        caption += 'Mass-observable uncertainty was taken into account considering a log-normal distribution with $\sigma_{\ln M}$ = ' + str(run_params['Mobs_sigma'])+'.\n'    

    caption += 'Area = '+ str(run_params['area']) +' $\\text{deg}^2$, $\M_{\\text{min}}$ = '+ str(run_params['Mi']) +' $\h^{-1} \\text{M}_{\odot}$ and redshift interval (' + str(run_params['z_initial']) + ', ' + str(run_params['z_final']) +').'

    f.write (caption)

#--------------------------------------------------------------------------------------------------------------------
##########################################
# Function to plot Fisher Matrix figures #
########################################## 
def plot_Fisher_Matrix (run_params, outFile1, outFile2, outFile3):
    # The function below is useful when you have Latex installed in your machine.
    #plt.rc('text', usetex=True)

    cosmo_params = {}
    cosmo_params[r'H_0'] = float(run_params['H0'])
    cosmo_params[r'\Omega_x'] = float(run_params['Omega_x'])
    cosmo_params[r'\Omega_b'] = float(run_params['Omega_b'])
    cosmo_params[r'\Omega_c'] = float(run_params['Omega_c'])
    cosmo_params[r'\sigma_8'] = float(run_params['Sigma_8'])
    cosmo_params[r'n_s'] = float(run_params['n_s_val'])
    cosmo_params[r'w_0'] = float(run_params['eqos_w0'])
    cosmo_params[r'w_1'] = float(run_params['eqos_w1'])
    
    list_plots = {}
    list_plots['fig'] = []
    file_bf = open(outFile1).read().splitlines()
    n = 0
    fitted_params = []
    bf_values = []
    for line in file_bf:
        pieces = line.split()
        if pieces[0]!= '#':
            if ("FREE" == pieces[1]):
                n = n + 1
                fitted_params.append(pieces[0])
                bf_values.append(pieces[2])

    if (run_params['fiducial_available'] == '0'):
       fiducial_params = []
       for item in range(len(fitted_params)):
           fiducial_params.append(cosmo_params[fitted_params[item]])
        
    file_MF = open(outFile2).read().splitlines()
    mfisher = []
    mfisher_row = []
    for line in file_MF:
        element = line.split()
        if element[0]!= '#':
            for k in xrange(0, n):
                mfisher_row.append(element[k])
            mfisher.append(mfisher_row)
            mfisher_row=[]
    
    k = 1
    fm_1sigma = 1.52
    fm_2sigma = 2.48
    fm_3sigma = 3.44
    for i in range(0, n):
        for j in range (i, n):
            plt.clf()
            plt.cla()
            if j == i:
                C11 = float(mfisher[i][i])
            else:
                C12 = float(mfisher[i][j])
                C22 = float(mfisher[j][j])
                p1_bf = [bf_values[i]]
                p2_bf = [bf_values[j]]
                
                a = np.sqrt((C11 + C22) / 2.0 + np.sqrt((C11 - C22) **2 / 4.0 + C12 * C12))
                b = np.sqrt((C11 + C22) / 2.0 - np.sqrt((C11 - C22) **2 / 4.0 + C12 * C12))
                theta = 0.5 * math.atan2(2.0 * C12, (C11 - C22))
                theta1 = theta * 180 / np.pi
                
                ellipse_1sigma = mpatches.Ellipse((p1_bf[0], p2_bf[0]), 2.0 * fm_1sigma * a, 2.0 * fm_1sigma * b, theta1)
                ellipse_2sigma = mpatches.Ellipse([p1_bf[0], p2_bf[0]], 2.0 * fm_2sigma * a, 2.0 * fm_2sigma * b, theta1)
                ellipse_3sigma = mpatches.Ellipse([p1_bf[0], p2_bf[0]], 2.0 * fm_3sigma * a, 2.0 * fm_3sigma * b, theta1)
                
                p1 = plt.plot(p1_bf, p2_bf, 'ko', ms=6)
                legend([p1], ['best-fit'], 'best')
                if (run_params['fiducial_available'] == '0'):
                    p2 = plt.plot(fiducial_params[i], fiducial_params[j], 'ro', ms=6)
                    legend([p1, p2], ['best-fit', 'fiducial'], 'best')
                ax = plt.subplot(111)
                
                ellipse_1sigma.set_alpha(0.1)
                ax.add_artist(ellipse_1sigma)
                ellipse_2sigma.set_alpha(0.1)
                ax.add_artist(ellipse_2sigma)
                ellipse_3sigma.set_alpha(0.1)
                ax.add_artist(ellipse_3sigma)
                
                margin_scale = 1.1 * fm_3sigma
                
                if theta >= 0.0 and theta < 1.570796327:
                    xl_lim = float(p1_bf[0]) - margin_scale * (a * np.cos(theta) + b * np.sin(theta))
                    xu_lim = float(p1_bf[0]) + margin_scale * (a * np.cos(theta) + b * np.sin(theta))
                    yl_lim = float(p2_bf[0]) - margin_scale * (a * np.sin(theta) + b * np.cos(theta))
                    yu_lim = float(p2_bf[0]) + margin_scale * (a * np.sin(theta) + b * np.cos(theta))
                    
                    ax.set_xlim(xl_lim, xu_lim)
                    ax.set_ylim(yl_lim, yu_lim)
                    
                elif theta >= 1.570796327 and theta < 3.141592654:
                    xl_lim = float(p1_bf[0]) + margin_scale * (a * np.cos(theta) - b * np.sin(theta))
                    xu_lim = float(p1_bf[0]) - margin_scale * (a * np.cos(theta) - b * np.sin(theta))
                    yl_lim = float(p2_bf[0]) - margin_scale * (a * np.sin(theta) - b * np.cos(theta))
                    yu_lim = float(p2_bf[0]) + margin_scale * (a * np.sin(theta) - b * np.cos(theta))
                    
                    ax.set_xlim(xl_lim, xu_lim)
                    ax.set_ylim(yl_lim, yu_lim)
                    
                elif theta >= -1.570796327 and theta < 0.0:
                    xl_lim = float(p1_bf[0]) - margin_scale * (a * np.cos(theta) - b * np.sin(theta))
                    xu_lim = float(p1_bf[0]) + margin_scale * (a * np.cos(theta) - b * np.sin(theta))
                    yl_lim = float(p2_bf[0]) + margin_scale * (a * np.sin(theta) - b * np.cos(theta))
                    yu_lim = float(p2_bf[0]) - margin_scale * (a * np.sin(theta) - b * np.cos(theta))
                    
                    ax.set_xlim(xl_lim, xu_lim)
                    ax.set_ylim(yl_lim, yu_lim)
                    
                else:
                    xl_lim = float(p1_bf[0]) + margin_scale * (a * np.cos(theta) + b * np.sin(theta))
                    xu_lim = float(p1_bf[0]) - margin_scale * (a * np.cos(theta) + b * np.sin(theta))
                    yl_lim = float(p2_bf[0]) + margin_scale * (a * np.sin(theta) + b * np.cos(theta))
                    yu_lim = float(p2_bf[0]) - margin_scale * (a * np.sin(theta) + b * np.cos(theta))
                    
                    ax.set_xlim(xl_lim, xu_lim)
                    ax.set_ylim(yl_lim, yu_lim)
                    
                xlabel = fitted_params[i]            
                ylabel = fitted_params[j]
                plt.xlabel('$' + xlabel + '$')
                plt.ylabel('$' + ylabel + '$')
                t = plt.title('Fisher Matrix')
                name = 'fisher_' + str(k) + '.png'
                caption_filename = 'caption_' +str(k) + '.txt'
                k = k + 1
                plt.savefig(name)
                build_caption (caption_filename, run_params, fitted_params, i, j)
                list_plots['fig'].append ([name, caption_filename])
                
                ellipse_1sigma = 0
                ellipse_2sigma = 0
                ellipse_3sigma = 0
    
    return list_plots          
              
#--------------------------------------------------------------------------------------------

###############################################
# Function to plot Profile Likelihood figures #
###############################################
def plot_Profile_Likelihood (run_params, param1, param2, outFile1, outFile2, outFile3):
    #plt.rc('text', usetex=True) 

    cosmo_params = {}
    cosmo_params[r'H_0'] = float(run_params['H0'])
    cosmo_params[r'T_{\gamma}'] = 0.0
    cosmo_params[r'\Omega_x'] = float(run_params['Omega_x'])
    cosmo_params[r'\Omega_b'] = float(run_params['Omega_b'])
    cosmo_params[r'\Omega_c'] = float(run_params['Omega_c'])
    cosmo_params[r'\sigma_8'] = float(run_params['Sigma_8'])
    cosmo_params[r'n_s'] = float(run_params['n_s_val'])
    cosmo_params[r'w_0'] = float(run_params['eqos_w0'])
    cosmo_params[r'w_1'] = float(run_params['eqos_w1'])

    list_plots = {}
    list_plots['fig'] = []
    file_bf = open(outFile1).read().splitlines()
    fitted_params = []
    bf_values = []
    for line in file_bf:
        pieces = line.split()
        if pieces[0] != '#':
            fitted_params.append([pieces[0], pieces[1]])
            bf_values.append(pieces[2])
    
    if (run_params['fiducial_available'] == '0'):
       fiducial_params = []
       for item in range(len(fitted_params)):
           fiducial_params.append(cosmo_params[fitted_params[item][0]])

    file = open(outFile2).read().splitlines()
    regions = []
    count_reg = 0
    saved = 0
    x1 = []
    y1 = []
    for row in file:
        col = row.split()
        if not col:
            count_reg = count_reg + 1
            if count_reg == 2:
                t1 = x1
                t2 = y1
                regions.append([t1, t2])
                x1 = []
                y1 = []
                count_reg = 0
                saved = 1
        elif col[0].startswith('#'):
            continue
        else:
            x1.append(col[0])
            y1.append(col[1])
            saved = 0
    if saved == 0:
        regions.append([x1, y1])
 
    p1_bf = [bf_values[int(param1)]]
    p2_bf = [bf_values[int(param2)]]
    
    # plot best-fit point #
    p1 = plt.plot(p1_bf, p2_bf, 'ko', ms=6)
    legend([p1], ['best-fit'], 'best')
    if (run_params['fiducial_available'] == '0'):
        p2 = plt.plot(fiducial_params[int(param1)], fiducial_params[int(param2)], 'ro', ms=6)
        legend([p1, p2], ['best-fit', 'fiducial'], 'best')
    #ax = plt.subplot(111)
    # plot 1 sigma confidence region #
    plt.plot(regions[0][0], regions[0][1], 'b')
    # plot 2 sigma confidence region #
    plt.plot(regions[1][0], regions[1][1], 'b')
    # plot 3 sigma confidence region #
    plt.plot(regions[2][0], regions[2][1], 'b')
    ax = plt.subplot(111)
    i = int (param1)
    j = int (param2)
    xlabel = fitted_params[i][0]
    ylabel = fitted_params[j][0]
    plt.xlabel('$' + xlabel + '$')
    plt.ylabel('$' + ylabel + '$')
    name = 'profile.png'
    caption_filename = 'caption.txt'
    plt.title('Profile Likelihood')
    plt.savefig(name)
    build_caption (caption_filename, run_params, fitted_params, i, j)
    list_plots['fig'].append ([name, caption_filename])
    plt.clf()

    return list_plots    
#--------------------------------------------------------------------------------------------

################################
# Run in the DES-Brazil portal #
################################
def run():
    import lib.xmlutil as lx
    import tools.xml2ascii as x2a ### import package that translates XML to ASCII
    import Ft.Xml as xml ### import package for the creation of the XML product log
    import orchestration.io as cpio ### import package that handles input, output and configuration XML files
    import lib.log
    import lib.config

    print (" ---------------------------------------------------------")
    print ("  CosmoLibGui: Halo Abundance")
    print ("  Fit cosmologial parameters using halo number counts.")
    print ("  Author: Mariana Penna-Lima and Sandro D. P. Vitenti")
    print ("  Wrappist monk: Mariana Penna-Lima")
    print ("  Thanks to Angelo Fausti, Ricardo Ogando, Alexandre Pretto, Fernando Simoni and Sandro Vitenti")
    print (" ---------------------------------------------------------")
    
    logger = lib.log.get_logger ()
    logger.info ("Creating ComponentConfig and ComponentIO...")
    # initializing access to IO
    io = cpio.ComponentIO()
    # initializing access to Configuration
    conf = cpio.ComponentConfig()
    
    # getting the input catalog name
    """
    cfcmp should change to something like simulations/halo
    and you should ask for a Serch Engine of Halos only
    """
    run_params = {}

    # getting configuration parameters
    for id in params_ids:
        run_params[id] = conf.getScalarById (id)

    zi = run_params['z_initial']

    zf = run_params['z_final']

    mi = run_params['Mi']

    mf = run_params['Mf']

    area = run_params['area']
    print ('iniciando output')
    # output.xml: parameters
    io.putParam (_value=str(area), _id='area', _name='Area [square degrees]', _type='float',  _section='Sample Selection', _publish='True')
    io.putParam (_value=str(zi), _id='z_initial', _name='Minimum redshift', _type='float', _section='Sample Selection', _publish='True')
    io.putParam (_value=str(zf), _id='z_final', _name='Maximum redshift', _type='float', _section='Sample Selection', _publish='True')
    io.putParam (_value=str(mi), _id='Mi', _name='Lower mass threshold', _type='float', _section='Sample Selection', _publish='True')
    io.putParam (_value=str(mf), _id='Mf', _name='Upper mass threshold', _type='float', _section='Sample Selection', _publish='True')
    print ('terminando output')

    if (int(run_params['Mobs']) == 0):
        mobs = 'F'
        mobs_value = 'Not included'
        io.putParam (_value=str(mobs_value), _id='Mobs', _name='Mass-Observable probability distribution', _type='string',  _section='Observational Effects', _publish='True')
    else:
        mobs = 'T'
        mobs_value = 'Included'
        io.putParam (_value=str(mobs_value), _id='Mobs', _name='Mass-Observable probability distribution', _type='string',  _section='Observational Effects', _publish='True')
        io.putParam (_value=str(run_params['Mobs_sigma']), _id='Mobs_sigma', _name='Standard deviation', _type='float',  _section='Observational Effects', _publish='True')
    mobs_sd = run_params['Mobs_sigma']
    if (int(run_params['photoz']) == 0):  
        zphot = 'F'
        zphot_value = 'Not included'
        io.putParam (_value=str(zphot_value), _id='photoz', _name='Photo-z probability distribution', _type='string',  _section='Observational Effects', _publish='True')
    else:
        zphot = 'T'
        zphot_value = 'Included'
        io.putParam (_value=str(zphot_value), _id='photoz', _name='Photo-z probability distribution', _type='string',  _section='Observational Effects', _publish='True')
        io.putParam (_value=str(run_params['photoz_sigma']), _id='photoz_sigma', _name='Standard deviation', _type='float',  _section='Observational Effects', _publish='True')
    zphot_sd = run_params['photoz_sigma']
   
    if (int(run_params['estimator']) == 0):
        io.putParam (_value=str('Unbinned in z and M'), _id='estimator_nc', _name='Estimator', _type='float',  _section='Number Counts Estimator', _publish='True')
    elif (int(run_params['estimator']) == 0):
        io.putParam (_value=str('Unbinned in z'), _id='estimator_nc', _name='Estimator', _type='float',  _section='Number Counts Estimator', _publish='True')
    else:
        io.putParam (_value=str('Poisson (binned in z)'), _id='estimator_nc', _name='Estimator', _type='float',  _section='Number Counts Estimator', _publish='True')
        io.putParam (_value=str(zbins), _id='zbins_nc', _name='Number of redshift bins', _type='float',  _section='Number Counts Estimator', _publish='True')
    
    #catalog = os.getcwd() + '/'
    catalog = io.getFileByMimeClass('application/fits-table','simulation/halos')
    catalog += '[col #Z_I=' + str(zi)
    catalog += '; #Z_F=' + str(zf) 
    catalog += '; #M_I=' + str(mi) 
    catalog += '; #M_F=' + str(mf) 
    catalog += '; #AREA=' + str(area) 
    catalog += '; #M_OBS=' + str(mobs) 
    catalog += '; #M_OBS_S0=' + str(mobs_sd) 
    catalog += '; #PHOTOZ=' + str(zphot) 
    catalog += '; #PHOTOZS0=' + str(zphot_sd)  + ']'    

    print (catalog)

    run_params['halo_cat'] = catalog
    run_params['outFile'] = 'halo_abundance.log'
    
    comando = build_command(run_params)
                 
    list_files = run_darkenergy (comando, False, run_params)

    print (os.getcwd())
    
    #i = 0
    # output.xml: files
    #os.system("tar -cvf halo_abundance.tar "+" ".join(list_files['txt']))
    all_files = list_files['txt']
    if (int(len(list_files['fig']) > 0)):
        for fig in list_files['fig']:
            all_files += fig
    print (all_files)
    cmd_to_call = "tar -cvf halo_abundance.tar " + " ".join(all_files)
    print (cmd_to_call)
    os.system(cmd_to_call)
    #os.system("tar -cvf halo_abundance.tar "+" ".join(list_files['fig']))
    #making the output file for product-log for portal
    io.addOutput(_value='halo_abundance.tar',_class='science',_mimetype='application/x-tar',_name='',_output_type='file',_publish="True")
    #for file in list_files['txt']:
        #filename = str(file)
        #io.putFile (_value=filename,_class='science',_mimetype='text/plain',_publish="True")
        #io.putFile (_value=filename0,_class='science',_mimetype='text/plain',_publish="True")
    #for figure in list_files['fig']:
        #print 'chegou em fig'
        #filename0 = str(figure[0])
        #filename1 = str(figure[1])
        #io.putFile (_value=filename0,_class='cosmohalo',_mimetype='image/png',_publish="True")
        #io.putFile (_value=filename1,_class='cosmohalo',_mimetype='text/plain',_publish="True")
   
    ##############################
    #  XML Product Log Creation  #
    ##############################

    print ('chegou productlog')

    io.beginLog('Cosmological Fits')
   
    io.beginSection('Data vs. Model')

    io.beginSubSection('Mass Function')
    io.beginStats()
    io.addStat('Multiplicity Function', 'qq coisa')
    io.endStats()
    io.endSubSection()
    io.endSection()

    io.beginSection('Cosmological Constraints')
    
    latex2utf8 = {}
    latex2utf8[r'H_0'] = 'H0 [km/(s Mpc)]'
    latex2utf8[r'T_{\gamma}'] = 'T&#947; [K]'
    latex2utf8[r'\Omega_x'] = '&#x3A9;x'
    latex2utf8[r'\Omega_b'] = '&#x3A9;b'
    latex2utf8[r'\Omega_c'] = '&#x3A9;c'
    latex2utf8[r'\sigma_8'] = '&#x3C3;8'
    latex2utf8[r'n_s'] = 'n_s'
    latex2utf8[r'w_0'] = 'w0'
    latex2utf8[r'w_1'] = 'w1'
    latex2utf8[r'$'] = ''
    latex2utf8[r'\sigma'] = '&#x3C3;'
    latex2utf8[r'\text{deg}^2'] = 'deg&#xB2;'
    latex2utf8[r'\M_{\text{min}}'] = 'M_min'
    latex2utf8[r'\h^{-1} \text{M}_{\odot}'] = 'h{-1} M&#8857;'

    for file_txt in list_files['txt']:
        name_pieces = file_txt.split('_')
        if str(name_pieces[0]) == 'best':
            best_file = open(file_txt)
            io.beginSubSection('Fixed Cosmological Parameters')
            io.beginStats()
            for line_fixed in best_file:
                if not line_fixed.startswith("#"):
                    pieces = line_fixed.split()
                    if (str(pieces[1]) == 'FIXED'):
                        parameter = latex2utf8[pieces[0]]
                        val = "%.3f" % float(pieces[2])
                        io.addStat(str(parameter), val)
            io.endStats()
            io.endSubSection()
            best_file.seek (0)
            
            if (int(run_params['statist_method']) == 0):
                i = 0
                err_list = []
                for file_txt in list_files['txt']:
                    name_pieces = file_txt.split('_')
                    if str(name_pieces[0]) == 'MF':
                        MF_file = open(file_txt)
                        for line_mf in MF_file:
                            if not line_mf.startswith("#"):
                                pieces_err = line_mf.split()
                                err_list.append("%.3f" % math.sqrt(float(pieces_err[i])))
                                print (err_list[i])
                                i = i + 1
                print (err_list)
                io.beginSubSection('Fitted Cosmological Parameters')
                io.beginStats()
                i = 0
                for line_free in best_file:
                    if not line_free.startswith("#"):
                        pieces = line_free.split()
                        if (str(pieces[1]) == 'FREE'):
                            parameter = latex2utf8[pieces[0]]
                            val = "%.3f" % float(pieces[2])
                            err = err_list[i]
                            io.addStat(str(parameter), str(val) +  ' &#xB1; ' + str(err))
                            i = i + 1
                            print (i)
                io.endStats()
                io.endSubSection()
            else:
                io.beginSubSection('Fitted cosmological parameters')
                io.beginStats()
                for line_free in best_file:
                    if not line_free.startswith("#"):
                        pieces = line_free.split()
                        if (str(pieces[1]) == 'FREE'):
                            parameter = latex2utf8[pieces[0]]
                            val = "%.3f" % float(pieces[2])
                            io.addStat(str(parameter), val)
                io.endStats()
                io.endSubSection()
    
    io.endSection()

    if (int(len(list_files['fig'])) > 0):
    
        io.beginSection('Confidence Regions')

        io.beginSubSection()
    
        for figure in list_files['fig']:
            cap = open(figure[1]) 
            caption = cap.read()
            for k in latex2utf8:
                caption = caption.replace(k, latex2utf8[k])
            print (caption)
            io.addPlot(figure[0],'', caption)
            print (figure[0])

        io.endSubSection()
        io.endSection()
    
    io.endLog()

    print ('The End')
#--------------------------------------------------------------------------------------------------------------------------

#######################################################################################################################
# This function calls darkenergy with especifications given by "comando" and also call the function to generate plots.#
#######################################################################################################################
def run_darkenergy (comando, create_subdir, run_params):
    # Open a directory and change the current dir to it.
    mycurdir = os.getcwd ()
    d2 = '.'
    if create_subdir:
        d1 = datetime.datetime.now()
        d2 = d1.strftime("%Y%m%dT%H%M%S")
        os.mkdir(d2)
        os.chdir(d2)
    
    # list with names' files to be published at the portal
    list_files = {}
    list_files['fig'] = []
    list_files['txt'] = []
    list_files['massfunc'] = []
    res = subprocess.Popen(comando).wait()
    list_files['txt'].append(run_params['outFile'])
    list_files['outdir'] = d2
    if (res):
        cmdline = ""
        for co in comando:
            cmdline += "'" + co + "' "
        print (cmdline)
        raise [Exception, 'Could not run command: ' + comando[0]]
    if (not (run_params['cr_y_axis'] == '' or run_params['outFile'] == '')):
        lines = open(run_params['outFile']).read().splitlines()
        outFile1 = ''
        outFile2 = ''
        outFile3 = ''
        for line in lines:
            if line.startswith('# Params file:'):
                pieces = line.split()
                outFile1 = pieces[3]
                list_files['txt'].append(outFile1)
            if (int(run_params['statist_method']) == 0):
                if line.startswith('# FM file:'):
                    pieces = line.split()
                    outFile2 = pieces[3]
                    list_files['txt'].append(outFile2)
                    list_files2 = plot_Fisher_Matrix (run_params, outFile1, outFile2, run_params['outFile'])
                    for k in list_files2.keys():
                        list_files[k].extend(list_files2[k])
            elif (int(run_params['statist_method']) == 1):
                if line.startswith('# PL file:'):
                    pieces = line.split()
                    outFile2 = pieces[3]
                    list_files['txt'].append(outFile2)
                    list_files2 = plot_Profile_Likelihood (run_params, run_params['cr_x_axis'], run_params['cr_y_axis'], outFile1, outFile2, run_params['outFile'])
                    for k in list_files2.keys():
                        list_files[k].extend(list_files2[k])
            if (int(run_params['massf_print']) == 0):
                if line.startswith('# MassFunction file:'):
                    pieces = line.split()
                    outFile3 = pieces[3]
                    list_files['txt'].append(outFile3)
                    list_files3 = plot_mass_function (run_params, outFile3)
    
    if create_subdir:
        os.chdir (mycurdir)
    return list_files

#--------------------------------------------------------------------------------------------------------------------------

try:
    import orchestration.io
except ImportError:
    pass
    try:
        filename = sys.argv[1]
    except IndexError:
        print ('Usage ' + sys.argv[0] + 'config_halo_abundance.ini')
        sys.exit(1)

    config = ConfigParser.SafeConfigParser()
    file = open(filename)
    config.readfp(file)

    run_params = {}
    for id in params_ids:
        run_params[id] = config.get('RunParams', id)
    run_params['halo_cat'] = config.get('RunParams', 'halo_cat')
    run_params['outFile'] = config.get('RunParams', 'outFile')

    comando = build_command(run_params)
    print (" ".join(comando))

    run_darkenergy (comando, True, run_params)


