#!/usr/bin/env python

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import sys
import time
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.interpolate import CubicSpline

from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm
from scipy.interpolate import interp1d

matplotlib.rcParams.update({'font.size': 11})

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

#
#  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm 
#
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm{'massnu-length':<1>}")
cosmo.omega_x2omega_k ()
cosmo.param_set_by_name ("Omegak",    0.0)
cosmo.param_set_by_name ("w",        -1.0)
cosmo.param_set_by_name ("Omegab",    0.022)
cosmo.param_set_by_name ("Omegac",    0.12)
cosmo.param_set_by_name ("massnu_0",  0.06)
cosmo.param_set_by_name ("ENnu",      2.0328)

reion = Nc.HIReionCamb.new ()
prim  = Nc.HIPrimPowerLaw.new ()

cosmo.param_set_by_name ("H0", 67)

prim.param_set_by_name ("n_SA", 0.96)
prim.param_set_by_name ("ln10e10ASA", 3.0904)

reion.param_set_by_name ("z_re", 9.9999)


cosmo.add_submodel (reion)
cosmo.add_submodel (prim)


#
#  Printing the parameters used.
#
print ("# Model parameters: ", end=' ') 
cosmo.params_log_all ()
print ("# Omega_X0: % 22.15g" % (cosmo.E2Omega_de (0.0)))


ps_lin = Nc.PowspecMLCBE.new ()

# Redshift bounds
z_min = 0.0
z_max = 2.0
zdiv  = 0.49999999999

# Mode bounds
k_min = 1.0e-8
k_max = 1.0e6

ps_lin.set_kmin (k_min)
ps_lin.set_kmax (k_max)
ps_lin.require_zi (z_min)
ps_lin.require_zf (z_max)

ps_lin.set_intern_k_min (k_min)
ps_lin.set_intern_k_max (50.0)
psf = Ncm.PowspecFilter.new (ps_lin, Ncm.PowspecFilterType.TOPHAT)
psf.set_best_lnr0 ()

ps = Nc.PowspecMNLHaloFit.new (ps_lin, 2.0, 1.0e-8)

ps.set_kmin (k_min)
ps.set_kmax (k_max)
ps.require_zi (z_min)
ps.require_zf (z_max)

old_amplitude = math.exp (prim.props.ln10e10ASA)
prim.props.ln10e10ASA = math.log ((0.81/ cosmo.sigma8(psf))**2 * old_amplitude)

ell_min = 0
ell_max = 1

sproj = Ncm.PowspecSphereProj.new (ps, ell_min, ell_max)

sproj.props.reltol = 1.0e-8
sproj.props.reltol_z = 1.0e-8
xi_i = 1.0e-1
xi_f = 1.0e7

sproj.set_xi_i (xi_i)
sproj.set_xi_f (xi_f)

sproj.prepare (cosmo)

dist = Nc.Distance.new (1.0e11)
dist.compute_inv_comoving (True)

scal = Nc.Scalefactor.new (1.0e11, dist)
gf   = Nc.GrowthFunc.new ()

dist.prepare (cosmo)
scal.prepare (cosmo)
gf.prepare (cosmo)

RH_Mpc = cosmo.RH_Mpc ()

#for ell in range (ell_min, ell_max + 1):
  #(lnxi, Cell) = sproj.get_ell (ell)
  #xi_a = np.exp (np.array (lnxi))

xi_a   = np.geomspace (xi_i, xi_f, num = 5000)
z_a    = np.array ([dist.inv_comoving (cosmo, xi / RH_Mpc) for xi in xi_a])
index  = np.logical_and ((z_a > 0.0), (z_a < 10.0))

z_a    = np.linspace (0.2, 2, num = 2001)[1:]




#print (xi_a)
#print (z_a)

for ell in range (ell_min, ell_min + 1):
  
  limber = np.array ([ps.eval (cosmo, z, (0.5 + ell) / xi) / (xi * xi * RH_Mpc / cosmo.E (z)) for (xi, z) in zip (xi_a, z_a)])
  
  #plt.plot (xi_a[index], limber[index], label = r'$\ell$ = %d, limber' % ell)
  
  #for w_i in range (20):
    #w      = sproj.get_w (w_i)
    #Cell_s = sproj.peek_ell_spline (w_i, ell)
    #Cell   = [Cell_s.eval (math.log (xi)) for xi in xi_a]

    #for xi, Cell_i in zip (xi_a, Cell):
      #sCell_i = sproj.eval_Cell_xi1_xi2 (ell, xi / w, xi * w)
      #print ("xi % 22.15g w % 22.15g Cell % 22.15g | % 22.15g %e" % (xi, w, Cell_i, sCell_i, sCell_i / Cell_i - 1.0))
      #Cell_int_i = ps.sproj (cosmo, 1.0e-5, ell, 0.0, 0.0, xi / w, xi * w)
      #print ("% 22.15g <% 22.15g> [% 22.15g <% 22.15g> % 22.15g <% 22.15g> : % 22.15g] % 22.15g % 22.15g %e % 22.15g" % (xi, dist.inv_comoving (cosmo, xi / RH_Mpc), xi * w, dist.inv_comoving (cosmo, xi * w / RH_Mpc), xi / w, dist.inv_comoving (cosmo, xi / w / RH_Mpc), xi * (1.0 / w - w), Cell_int_i, Cell_i, Cell_i / Cell_int_i - 1.0, RH_Mpc * dist.comoving (cosmo, dist.inv_comoving (cosmo, xi / RH_Mpc))))
  
    #z_a    = np.array (z_a)
    #Cell   = np.array ([(gf.eval (cosmo, z)**2 * Cell_i) for (Cell_i, z) in zip (Cell, z_a)])
  

    #plt.plot (xi_a, Cell, label = r'$\ell$ = %d, $w$ = %f, full' % (ell, w))
  m    = [[sproj.eval_Cell_xi1_xi2 (cosmo, ell, z1, z2, RH_Mpc * dist.comoving (cosmo, z1), RH_Mpc * dist.comoving (cosmo, z2)) for z1 in z_a] for z2 in z_a]
  cov  = np.asanyarray (m)
  #print (np.diag (cov))
  std_ = np.sqrt (np.diag (cov))
  corr = cov / np.outer(std_, std_)
  #print (m)
  #print (corr)  
  #plt.matshow (np.abs (corr), norm = LogNorm (vmin=1.0e-8, vmax = 1.0))
  #plt.matshow (corr, extent = np.array ([z_a[0], z_a[-1], z_a[0], z_a[-1]]), origin = "lower")
  plt.matshow (cov, extent = np.array ([z_a[0], z_a[-1], z_a[0], z_a[-1]]), origin = "lower", norm = SymLogNorm (1.0e-6))
  plt.xlabel(r'$z_2$')
  plt.ylabel(r'$z_1$')
  #plt.matshow (corr)
  #plt.matshow (corr, norm = SymLogNorm (1.0e-4))
  plt.title (r'$C_{%d}(z_1, z_2)$' % (ell))  
  plt.colorbar ()
  plt.show ()
#plt.xticks(xi_a)

#plt.xlabel (r'$\xi \; [\mathrm{Mpc}]$')
#plt.ylabel (r'$C_\ell(\xi, \xi)$')
#plt.legend (loc = "best")
#plt.xscale ('log')
#plt.yscale ('symlog', linthreshy=1.0e-8)
#plt.xlim ([200.0, 4000.0])

#plt.show ()
print(cov/np.log(10))


'''''''''
import sys
sys.path.insert(0, '../../Programas_Cosmologia/')
import PySSC
from PySSC import PySSC
zstakes = np.array([0.2,0.3,0.4,0.5,0.6,0.7,0.8,  ,1,1.1,1.2,1.3,1.4,1.5])

nz       = 2000
z_arr    = np.linspace(0.2,2,num=nz+1)[1:] 



nbins_T   = len(zstakes)-1
kernels_T = np.zeros((nbins_T,nz))
for i in range(nbins_T):
    zminbin = zstakes[i] ; zmaxbin = zstakes[i+1] ; Dz = zmaxbin-zminbin
    for iz in range(nz):
        z = z_arr[iz]
        if ((z>zminbin) and (z<=zmaxbin)):
            kernels_T[i,iz] = 1/Dz



Sijw_T_alt = PySSC.Sij_alt_fullsky(z_arr,kernels_T)
'''''''''

































"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Compute necessary FFTs and make interpolation functions
def lacasa_cll(z):
    nz = len(z)
    keq         = 0.02/cosmo.h()                            #Equality matter radiation in 1/Mpc (more or less)
    klogwidth   = 10                                        #Factor of width of the integration range.
                                                            #10 seems ok ; going higher needs to increase nk_fft to reach convergence (fine cancellation issue noted in Lacasa & Grain)
    kmin        = 1.0e-6
    kmax        = 1.0e2
    nk_fft      = 2**11                                     #seems to be enough. Increase to test precision, reduce to speed up.
    k_4fft      = np.linspace(kmin,kmax,nk_fft)             #linear grid on k, as we need to use an FFT
    Deltak      = kmax - kmin
    Dk          = Deltak/nk_fft
    Pk_4fft     = np.zeros(nk_fft)
    for ik in range(nk_fft):
        Pk_4fft[ik] = ps.eval(cosmo,0.,k_4fft[ik])               #In Mpc^3
    dr_fft      = np.linspace(0,nk_fft//2,nk_fft//2+1)*2*np.pi/Deltak

    
    fft0        = np.fft.rfft(Pk_4fft)*Dk
    dct0        = fft0.real ; dst0 = -fft0.imag
    Pk_dct      = interp1d(dr_fft,dct0,kind='cubic')
    
    # Compute sigma^2(z1,z2)
    sigma2_nog = np.zeros((nz,nz))
    #First with P(k,z=0) and z1<=z2
    for iz in range(nz):
        r1 =   dist.comoving (cosmo, z[iz])
        for jz in range(iz,nz):
            r2                =  dist.comoving (cosmo, z[jz])
            rsum              = r1+r2
            rdiff             = abs(r1-r2)
            Icp0              = Pk_dct(rsum) ; Icm0 = Pk_dct(rdiff)
            sigma2_nog[iz,jz] = (Icm0-Icp0)/(4*np.pi**2 * r1 * r2)
    #Now fill by symmetry and put back growth functions
    sigma2      = np.zeros((nz,nz))
    for iz in range(nz):
        growth1 = ps.eval(cosmo,z[iz],0.1)/ps.eval(cosmo,0.,0.1)
        for jz in range(nz):
            growth2       = ps.eval(cosmo,z[jz],0.1)/ps.eval(cosmo,0.,0.1)
            sigma2[iz,jz] = sigma2_nog[min(iz,jz),max(iz,jz)]*growth1*growth2

    return sigma2

print(lacasa_cll(z_a))
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
