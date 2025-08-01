# -*- coding: utf-8 -*-
"""Test2.py


"""

# Basic Libraries
import sys
import math
import numpy as np                                                # imports the Numpy Library
from tqdm import tqdm
from scipy.integrate import quad

from numcosmo_py import Nc, Ncm                                   # imports the NumCosmo library
from numcosmo_py.plotting.tools import set_rc_params_article      # imports Numcosmo plotting tools

# Background Quantities
class model_evol:
  def __init__(self,model,alpha_list,k=1.0):
    self.model = model
    self.time = alpha_list
    self.k = k
    self.w = model.props.w
    w = model.props.w
    self.evol = np.array( [ model.eom_eval(alpha,k) for alpha in alpha_list] )
    self.y = np.array( [ eqmotion.y for eqmotion in self.evol ] )
    #self.misc = [ print(1) for alpha in alpha_list ]

    # Masses
    self.ms = np.array( [ eqmotion.m_s for eqmotion in self.evol ] )
    self.mz = np.array( [ eqmotion.m_zeta for eqmotion in self.evol ] )

    self.mnu2_z = np.array( [ eqmotion.mnu2_zeta for eqmotion in self.evol ] )
    self.mnu2_s = np.array( [ eqmotion.mnu2_s for eqmotion in self.evol] )

    # Frequencies
    self.nu1 = np.array( [ eqmotion.nu1 for eqmotion in self.evol ] )
    self.nu2 = np.array( [ eqmotion.nu2 for eqmotion in self.evol ] )
    self.nuz = np.sqrt( self.mnu2_z / self.mz )
    self.nus = np.sqrt( self.mnu2_s / self.ms )

    self.F = self.nu1 / np.sqrt(1/3)

    # Coupling Matrices
    self.gamma11 = np.array( [ eqmotion.gammabar11 for eqmotion in self.evol ] )
    self.gamma12 = np.array( [ eqmotion.gammabar12 for eqmotion in self.evol ] )
    self.gamma22 = np.array( [ eqmotion.gammabar22 for eqmotion in self.evol ] )
    self.tau12 = np.array( [ eqmotion.taubar for eqmotion in self.evol ] )

    # Angles
    #self.cos2phi = ( (self.nu1 ** 2) * ( self.nuz ** 2 ) - (self.nu2 ** 2) * ( self.nus ** 2 ) ) / ( self.nu1**4 - self.nu2**4 )
    #self.sin2phi = ( (self.nu1 ** 2) * ( self.nus ** 2 ) - (self.nu2 ** 2) * ( self.nuz ** 2 ) ) / ( self.nu1**4 - self.nu2**4 )
    self.cos2phi = np.array( [ eqmotion.cos2phi for eqmotion in self.evol ] )
    self.sin2phi = np.array( [ eqmotion.sin2phi for eqmotion in self.evol ] )

    # Sound Velocities
    c1 = np.sqrt(1/3)
    c2 = np.sqrt( model.props.w )
    #self.cs2 = (c1 ** 2) * self.cos2phi + (c2 ** 2) * self.sin2phi
    #self.cm2 = (c2 ** 2) * self.cos2phi + (c1 ** 2) * self.sin2phi
    self.cs2 = np.array( [ eqmotion.cs2 for eqmotion in self.evol ] )
    self.cm2 = np.array( [ eqmotion.cm2 for eqmotion in self.evol ] )

    #
    self.nubar1 = self.nu1 * ( self.nu1 - self.gamma11 )
    self.nubar2 = self.nu2 * ( self.nu2 - self.gamma22 )

    #
    #self.sigma1 = prim( self.nu1 - self.gamma11/2, alpha_list )
    #self.sigma2 = prim( self.nu2 -self.gamma22/2, alpha_list  )

    # Effective Couplings
    varrho = 1 / ( 1 - (self.y ** 2) * self.ms * self.mz )

    self.Rh_z = varrho * ( np.gradient( self.mz, alpha_list )/ self.mz + self.y * self.mz * self.ms * np.gradient(self.y,alpha_list) )
    self.Rh_s = varrho * ( np.gradient( self.ms, alpha_list )/ self.ms + self.y * self.mz * self.ms * np.gradient(self.y,alpha_list))
    self.aleph_zs = self.y * self.mnu2_s
    self.aleph_sz = self.y * self.mnu2_z
    self.beth_zs = varrho * (self.ms / self.mz) * np.gradient( self.mz * self.y , alpha_list )
    self.beth_sz = varrho * (self.mz / self.ms) * np.gradient( self.ms * self.y , alpha_list)
    self.wz = (self.nuz**2) - (1/2) * ( np.gradient( self.Rh_z, alpha_list ) + (self.Rh_z)**2 )
    self.ws = (self.nus**2) - (1/2) * ( np.gradient( self.Rh_s, alpha_list ) + (self.Rh_s)**2 )

    # critical times
    tanphi = np.sqrt( self.sin2phi / self.cos2phi )
    tanphiz = np.sqrt(1/(3*w))
    tanphiS = np.sqrt(3*w)

    self.phi = np.arctan( tanphi )
    self.phiz = np.arctan( np.sqrt(1/(3*w)) )
    self.phiS = np.arctan( np.sqrt(3*w) )

    for i in range(len(alpha_list)):
      if np.abs(tanphi[i] - tanphiz) == np.min(np.abs(tanphi - tanphiz)):
        tz = alpha_list[i]
        phiz = self.phi[i]
      if np.abs(tanphi[i] - 1) == np.min(np.abs(tanphi - 1)):
        tc = alpha_list[i]
      if np.abs(tanphi[i] - tanphiS) == np.min(np.abs(tanphi - tanphiS)):
        ts = alpha_list[i]
        phiS = self.phi[i]

    self.phiz = phiz
    self.phiS = phiS

    self.tz = tz
    self.ts = ts
    self.tc = tc

  # Relevant time instants
  def trans(self):
    alpha_list = self.time
    t_list = []
    for i in range(len(self.cos2phi )):
      d = np.abs(  self.cos2phi[i]  - self.sin2phi[i]   )
      if d < 0.99:
        t_list.append(alpha_list[i])
    t_i = t_list[0]
    t_f = t_list[-1]
    return [t_i,t_f]

  def transH(self):
    alpha_list = self.time
    t_list = []
    Δ = 0.1

    for i in range(len(self.Rh_z )):
      d = np.abs(  ( self.Rh_z[i]  - self.Rh_z[0] ) / self.Rh_z[i]    )
      if d > Δ:
        t1 = alpha_list[i]
        break

    alpha_list_flip = np. flip(self.time)
    for i in range(len(self.Rh_z )):
      d = np.abs(  ( self.Rh_z[i]  - self.Rh_z[-1] ) / self.Rh_z[i]    )
      #if d >  6 * Δ * np.abs( self.Rh_z[-1]/self.Rh_z[0] ):
      if d >  0.399:
        t2 = alpha_list_flip[i]
        break

    return [t1,t2]

  def eq(self):
    alpha_list = self.time
    d = np.abs( ( self.mz  - self.ms ) / (self.mz)   )
    for i in range(len( alpha_list )):
      if d[i] == np.min(d):
        t_eq = alpha_list[i]
        break
    return t_eq

  def cross_z(self,k_value):
    alpha_list = self.time
    d1 = np.abs( np.log( np.abs( self.Rh_z/ ( ( self.nuz**2 ) ) * ( (self.k/k_value)**2 ) ) )   )
    for i in range(len( self.Rh_z )):
      if d1[i] == np.min(d1):
        t_cross1 = alpha_list[i]
        break

    t_cross = t_cross1

    #d2 = np.abs( np.log( np.abs( ( self.Rh_z / self.aleph_zs )  ) ) )
    #for i in range(len( self.Rh_z )):
    #  if d2[i] == np.min(d2):
    #    t_cross2 = alpha_list[i]
    #    break
    #  else:
    #    t_cross2 = alpha_list[0]

    #t_cross = np.max([t_cross1, t_cross2])

    return t_cross

  def cross_s(self, k_value):
    alpha_list = self.time
    d1 = np.abs( ( self.Rh_s - ( self.nus**2 ) * ( (k_value/self.k)**2 ) ) / (self.Rh_s)   )
    for i in range(len( self.Rh_s )):
      if d1[i] == np.min(d1):
        t_cross = alpha_list[i]
        break

    return t_cross

  def cross_za(self):
    alpha_list = self.time
    d = np.abs( np.log( np.abs( ( self.Rh_s / self.aleph_zs )  ) ) )
    for i in range(len( self.Rh_s )):
      if d[i] == np.min(d):
        t_cross = alpha_list[i]
        break

    return t_cross

  def cross_zb(self):
    alpha_list = self.time
    d = np.abs( np.log( np.abs( ( self.Rh_s / self.beth_zs )  ) ) )
    for i in range(len( self.Rh_s )):
      if d[i] == np.min(d):
        t_cross = alpha_list[i]
        break

    return t_cross

  def cross_sa(self):
    alpha_list = self.time
    d = np.abs( np.log( np.abs( ( self.Rh_s / self.aleph_sz )  ) ) )
    for i in range(len( self.Rh_s )):
      if d[i] == np.min(d):
        t_cross = alpha_list[i]
        break

    return t_cross

  def cross_sb(self):
    alpha_list = self.time
    d = np.abs( np.log( np.abs( ( self.Rh_s / self.beth_sz )  ) ) )
    for i in range(len( self.Rh_s )):
      if d[i] == np.min(d):
        t_cross = alpha_list[i]
        break

    return t_cross

  # horizon crossing time for various modes
  def cross_Z(self,k_interval):
    cross_z_list = []
    for k_value in k_interval:
      cross_z_list.append(self.cross_z(k_value))
    return cross_z_list

  def cross_S(self,k_interval):
    cross_s_list = []
    for k in k_interval:
      cross_s_list.append(self.cross_s(k))
    return cross_s_list

  # critical k values

  def kz1(self,k_interval=np.geomspace(1e1,1e3,1000)):
    difz1 = np.abs( self.cross_Z(k_interval) - self.tz )
    m = np.min(difz1)
    for i in range(len(k_interval)):
      if difz1[i] == m:
        kzeta1 = k_interval[i]
        break
    return kzeta1

  def kz2(self,k_interval=np.geomspace(1e6,1e10,1000)):
    difz2 = np.abs( self.cross_Z(k_interval) - self.ts )
    m = np.min(difz2)
    for i in range(len(k_interval)):
      if difz2[i] == m:
        kzeta2 = k_interval[i]
        break
    return kzeta2

  def ks1(self,k_interval):
    difs1 = np.abs( self.cross_S(k_interval)- self.tz )
    m = np.min(difs1)
    for i in range(len(k_interval)):
      if difs1[i] == m:
        kS1 = k_interval[i]
        break
    return kS1

  def ks2(self,k_interval):
    difs2 = np.abs( self.cross_S(k_interval) - self.ts )
    m = np.min(difs2)
    for i in range(len(k_interval)):
      if difs2[i] == m:
        kS2 = k_interval[i]
        break
    return kS2

  def kcrit(self, k_interval= np.geomspace(1e3,1e6,1000) ):
    dif_c = np.abs( np.log( np.array( self.cross_Z(k_interval) ) / np.array( self.cross_S(k_interval) ) ) )
    for i in range(len(k_interval)):
      if dif_c[i] == np.min(dif_c):
        kc = k_interval[i]
        break
    return kc

####################################################################################################################################################

# 2. Perturbative Quantities

def get_zeta(v):
    return v.get(Nc.HIPertITwoFluidsVars.ZETA_R) + 1.0j * v.get(Nc.HIPertITwoFluidsVars.ZETA_I)

def get_S(v):
    return v.get(Nc.HIPertITwoFluidsVars.S_R) + 1.0j * v.get(Nc.HIPertITwoFluidsVars.S_I)

def get_Pzeta(v):
    return v.get(Nc.HIPertITwoFluidsVars.PZETA_R) + 1.0j * v.get(Nc.HIPertITwoFluidsVars.PZETA_I)

def get_PS(v):
    return v.get(Nc.HIPertITwoFluidsVars.PS_R) + 1.0j * v.get(Nc.HIPertITwoFluidsVars.PS_I)

# OLD integrate system

def integrate_system(cosmo,k):
    # Defining relative tolerance for integration
    prec = 1.0e-10
    # Ratio potential frequency to define cross time
    cross_size = 1.0e-9
    # cross_size = 1.0e-9 # original choice
    # cross_size = 1.0e-8 # leads to very long runtimes

    # New perturbations object
    pert1 = Nc.HIPertTwoFluids.new()
    pert2 = Nc.HIPertTwoFluids.new()
    # Setting reltol
    pert1.props.reltol = prec
    pert2.props.reltol = prec
    # Setting k
    pert1.set_mode_k(k)
    pert2.set_mode_k(k)

    # Choose an initial condition
    alpha_try = -cosmo.abs_alpha(1.0e-14 * k**2)

    # New vector to store initial conditions
    # 8 dimensional (Q_1, Q_2, P_1, P_2), real and imaginary parts
    ci1 = Ncm.Vector.new(8)
    ci2 = Ncm.Vector.new(8)

    alphai1 = pert1.get_cross_time(cosmo, Nc.HIPertTwoFluidsCross.MODE1MAIN, alpha_try, cross_size)
    alphai2 = pert2.get_cross_time(cosmo, Nc.HIPertTwoFluidsCross.MODE2MAIN, alpha_try, cross_size)
    # Sandro tried cross_size x 10-2 for mode 2

    # Compute initial conditions at alpha_try store at ci and use normalization factor
    # pi/4
    pert1.get_init_cond_zetaS(cosmo, alphai1, 1, 0.25 * math.pi, ci1)
    pert2.get_init_cond_zetaS(cosmo, alphai2, 2, 0.25 * math.pi, ci2)

    # Use the previously computed initial conditions to start the system at alpha_try
    pert1.set_init_cond(cosmo, alphai1, 30, False, ci1)
    pert2.set_init_cond(cosmo, alphai2, 30, False, ci2)
    # print(f"Setting initial conditions for zeta1 and S1 at {alphai1}")
    # print(f"Setting initial conditions for zeta2 and S2 at {alphai2}")

    if alphai2 > alphai1:
        pert1.evolve(cosmo, alphai2)
        ci1, _ = pert1.peek_state(cosmo)
    else:
        pert2.evolve(cosmo, alphai1)
        ci2, _ = pert2.peek_state(cosmo)

    alphai = max(alphai1, alphai2)

    # Create a array of times to integrate the system over
    alpha_evol = np.linspace(alphai, -1.0e-1, 1000)

    # Integrate the system by stepping through alpha_evol using .evolve
    zeta1_a = [get_zeta(ci1)]
    S1_a = [get_S(ci1)]
    Pzeta1_a = [get_Pzeta(ci1)]
    PS1_a = [get_PS(ci1)]

    zeta2_a = [get_zeta(ci2)]
    S2_a = [get_S(ci2)]
    Pzeta2_a = [get_Pzeta(ci2)]
    PS2_a = [get_PS(ci2)]

    for alpha in tqdm(alpha_evol[1:], desc="Time evolution", position=1, leave=False):
        pert1.evolve(cosmo, alpha)
        pert2.evolve(cosmo, alpha)
        v1, _alphac1 = pert1.peek_state(cosmo)
        v2, _alphac2 = pert2.peek_state(cosmo)

        zeta1_a.append(get_zeta(v1))
        S1_a.append(get_S(v1))
        Pzeta1_a.append(get_Pzeta(v1))
        PS1_a.append(get_PS(v1))

        zeta2_a.append(get_zeta(v2))
        S2_a.append(get_S(v2))
        Pzeta2_a.append(get_Pzeta(v2))
        PS2_a.append(get_PS(v2))

    zeta1 = np.array(zeta1_a)
    S1 = np.array(S1_a)
    Pzeta1 = np.array(Pzeta1_a)
    PS1 = np.array(PS1_a)

    zeta2 = np.array(zeta2_a)
    S2 = np.array(S2_a)
    Pzeta2 = np.array(Pzeta2_a)
    PS2 = np.array(PS2_a)

    #cosmo_e = model_evol(cosmo, alpha_evol ,k)

    #Q11 = ( 1/np.sqrt( nu1(cosmo_e) ) ) * ( np.sqrt( mz(cosmo_e) ) * nu_z(cosmo_e) * np.sqrt( cos2phi(cosmo_e) ) * zeta1 + np.sqrt( ms(cosmo_e) ) * nu_s(cosmo_e) * np.sqrt( sin2phi(cosmo_e) ) * S1 )
    #Q21 = ( 1/np.sqrt( nu2(cosmo_e) ) ) * ( -np.sqrt( mz(cosmo_e) ) * nu_z(cosmo_e) * np.sqrt( sin2phi(cosmo_e) ) * zeta1 + np.sqrt( ms(cosmo_e) ) * nu_s(cosmo_e) * np.sqrt( cos2phi(cosmo_e) ) * S1 )
    #Q12 = ( 1/np.sqrt( nu1(cosmo_e) ) ) * ( np.sqrt( mz(cosmo_e) ) * nu_z(cosmo_e) * np.sqrt( cos2phi(cosmo_e) ) * zeta2 + np.sqrt( ms(cosmo_e) ) * nu_s(cosmo_e) * np.sqrt( sin2phi(cosmo_e) ) * S2 )
    #Q22 = ( 1/np.sqrt( nu2(cosmo_e) ) ) * ( -np.sqrt( mz(cosmo_e) ) * nu_z(cosmo_e) * np.sqrt( sin2phi(cosmo_e) ) * zeta2 + np.sqrt( ms(cosmo_e) ) * nu_s(cosmo_e) * np.sqrt( cos2phi(cosmo_e) ) * S2 )

    Q11 = []
    Q21 = []
    Q12 = []
    Q22 = []
    for i in range(len(alpha_evol)):
      dummy = cosmo.tv_eval(alpha_evol[i],k)
      A = dummy.zeta[0]
      B = dummy.zeta[1]
      C = dummy.s[0]
      D = dummy.s[1]

      det = A * D - B * C
      Anew = D/det
      Bnew = -B/det
      Cnew = -C/det
      Dnew = A/det

      Q11.append( Anew * zeta1[i] * Bnew * S1[i] )
      Q21.append( Cnew * zeta1[i] * Dnew * S1[i] )
      Q12.append( Anew * zeta2[i] * Bnew * S2[i] )
      Q22.append( Cnew * zeta2[i] * Dnew * S2[i] )

    Q11 = np.array(Q11)
    Q21 = np.array(Q21)
    Q12 = np.array(Q12)
    Q22 = np.array(Q22)

    return (alpha_evol, zeta1, S1, Pzeta1, PS1, zeta2, S2, Pzeta2, PS2, Q11, Q21, Q12, Q22)
    #return (alpha_evol, zeta1, S1, Pzeta1, PS1, zeta2, S2, Pzeta2, PS2, Q11, Q21, Q12, Q22, ci1.dup_array(), ci2.dup_array())

def spectrum(model,k_range):
  PI_zeta1 = []
  PI_zeta2 = []
  PI_S1 = []
  PI_S2 = []

  for k in tqdm(k_range, desc= "Mode evolution", position=0):
    alpha_evol, zeta1, S1, Pzeta1, PS1, zeta2, S2, Pzeta2, PS2, Q11, Q21, Q12, Q22 = integrate_system(model,k)

    PI_zeta1.append(np.abs(zeta1[-1])**2)
    PI_zeta2.append(np.abs(zeta2[-1])**2)
    PI_S1.append(np.abs(S1[-1])**2)
    PI_S2.append(np.abs(S2[-1])**2)

  PI_zeta1 = np.array(PI_zeta1)
  PI_zeta2 = np.array(PI_zeta2)
  PI_S1 = np.array(PI_S1)
  PI_S2 = np.array(PI_S2)

  As1 = np.polyfit(np.log(k_range),np.log(k_range**3 * PI_zeta1),1)[1]
  ns1 = np.polyfit(np.log(k_range),np.log(k_range**3 * PI_zeta1),1)[0] + 1
  As2 = np.polyfit(np.log(k_range),np.log(k_range**3 * PI_zeta2),1)[1]
  ns2 = np.polyfit(np.log(k_range),np.log(k_range**3 * PI_zeta2),1)[0] + 1
  AsT = np.polyfit(np.log(k_range),np.log(k_range**3 * np.sqrt( PI_zeta1**2 + PI_zeta2**2 ) ),1)[1]
  nsT = np.polyfit(np.log(k_range),np.log(k_range**3 * np.sqrt( PI_zeta1**2 + PI_zeta2**2 ) ),1)[0] + 1

  return [ns1, As1, ns2, As2, PI_zeta1, PI_S1, PI_zeta2, PI_S2, alpha_evol, nsT, AsT]

# new integrate system
def integrate_system(cosmo,k):
    # Defining relative tolerance for integration
    prec = 1.0e-10
    # Ratio potential frequency to define cross time
    cross_size = 1.0e-9

    # New perturbations object
    pert1 = Nc.HIPertTwoFluids.new()
    pert2 = Nc.HIPertTwoFluids.new()
    # Setting reltol
    pert1.props.reltol = prec
    pert2.props.reltol = prec
    # Setting k
    pert1.set_mode_k(k)
    pert2.set_mode_k(k)

    # Choose an initial condition
    alpha_try = -cosmo.abs_alpha(1.0e-14 * k**2)

    # New vector to store initial conditions
    # 8 dimensional (Q_1, Q_2, P_1, P_2), real and imaginary parts
    ci1 = Ncm.Vector.new(8)
    ci2 = Ncm.Vector.new(8)

    alphai1 = pert1.get_cross_time(cosmo, Nc.HIPertTwoFluidsCross.MODE1MAIN, alpha_try, cross_size)
    alphai2 = pert2.get_cross_time(cosmo, Nc.HIPertTwoFluidsCross.MODE2MAIN, alpha_try, cross_size)

    # Compute initial conditions at alpha_try store at ci and use normalization factor
    # pi/4
    pert1.get_init_cond_zetaS(cosmo, alphai1, 1, 0.25 * math.pi, ci1)
    pert2.get_init_cond_zetaS(cosmo, alphai2, 2, 0.25 * math.pi, ci2)

    # Use the previously computed initial conditions to start the system at alpha_try
    pert1.set_init_cond(cosmo, alphai1, 30, False, ci1)
    pert2.set_init_cond(cosmo, alphai2, 30, False, ci2)
    # print(f"Setting initial conditions for zeta1 and S1 at {alphai1}")
    # print(f"Setting initial conditions for zeta2 and S2 at {alphai2}")

    if alphai2 > alphai1:
        pert1.evolve(cosmo, alphai2)
        ci1, _ = pert1.peek_state(cosmo)
    else:
        pert2.evolve(cosmo, alphai1)
        ci2, _ = pert2.peek_state(cosmo)

    alphai = max(alphai1, alphai2)

    # Create a array of times to integrate the system over
    alpha_evol = np.linspace(alphai, -1.0e-1, 1000)

    # Integrate the system by stepping through alpha_evol using .evolve
    zeta1_a = [get_zeta(ci1)]
    S1_a = [get_S(ci1)]
    Pzeta1_a = [get_Pzeta(ci1)]
    PS1_a = [get_PS(ci1)]

    zeta2_a = [get_zeta(ci2)]
    S2_a = [get_S(ci2)]
    Pzeta2_a = [get_Pzeta(ci2)]
    PS2_a = [get_PS(ci2)]

    for alpha in tqdm(alpha_evol[1:], desc="Time evolution", position=1, leave=False):
        pert1.evolve(cosmo, alpha)
        pert2.evolve(cosmo, alpha)
        v1, _alphac1 = pert1.peek_state(cosmo)
        v2, _alphac2 = pert2.peek_state(cosmo)

        zeta1_a.append(get_zeta(v1))
        S1_a.append(get_S(v1))
        Pzeta1_a.append(get_Pzeta(v1))
        PS1_a.append(get_PS(v1))

        zeta2_a.append(get_zeta(v2))
        S2_a.append(get_S(v2))
        Pzeta2_a.append(get_Pzeta(v2))
        PS2_a.append(get_PS(v2))

    zeta1 = np.array(zeta1_a)
    S1 = np.array(S1_a)
    Pzeta1 = np.array(Pzeta1_a)
    PS1 = np.array(PS1_a)

    zeta2 = np.array(zeta2_a)
    S2 = np.array(S2_a)
    Pzeta2 = np.array(Pzeta2_a)
    PS2 = np.array(PS2_a)

    return (alpha_evol, zeta1, S1, Pzeta1, PS1, zeta2, S2, Pzeta2, PS2)

# 1. Integrates the first initial condition
def integrate_mode1(k, cosmo, max_time=-1.0):
    # Defining relative tolerance for integration
    prec = 1.0e-7
    # Ratio potential frequency to define cross time
    cross_size = 1.0e-13

    pert = Nc.HIPertTwoFluids.new()
    pert.set_stiff_solver(True)
    pert.props.reltol = prec
    pert.set_mode_k(k)
    alpha_try = -cosmo.abs_alpha(1.0e-14 * k**2)
    ci = Ncm.Vector.new(8)
    alphai = pert.get_cross_time(cosmo, Nc.HIPertTwoFluidsCross.MODE1MAIN, alpha_try, cross_size)
    pert.get_init_cond_zetaS(cosmo, alphai, 1, 0.25 * math.pi, ci)
    pert.set_init_cond(cosmo, alphai, 1, False, ci)

    res = pert.evolve_array(cosmo, max_time)
    res_a = np.array(res.dup_array()).reshape(-1,9)

    return res_a

# 2. Integrates the second initial condition
def integrate_mode2(k, cosmo, max_time=-1.0):
    # Defining relative tolerance for integration
    prec = 1.0e-7
    # Ratio potential frequency to define cross time
    cross_size = 1.0e-5

    pert = Nc.HIPertTwoFluids.new()
    pert.set_stiff_solver(True)
    pert.props.reltol = prec
    pert.set_mode_k(k)
    alpha_try = -cosmo.abs_alpha(1.0e-14 * k**2)
    ci = Ncm.Vector.new(8)
    alphai = pert.get_cross_time(cosmo, Nc.HIPertTwoFluidsCross.MODE2MAIN, alpha_try, cross_size)
    pert.get_init_cond_zetaS(cosmo, alphai, 2, 0.25 * math.pi, ci)
    pert.set_init_cond(cosmo, alphai, 2, False, ci)

    res = pert.evolve_array(cosmo, max_time)
    res_a = np.array(res.dup_array()).reshape(-1,9)

    #del pert
    #gc.collect()
    return res_a

# Slower integrating functions, considered for testing purposes

def integrate_mode1_qp(k, cosmo, max_time=-1.0):
    # Defining relative tolerance for integration
    prec = 1.0e-11
    # Ratio potential frequency to define cross time
    cross_size = 1.0e-13

    pert = Nc.HIPertTwoFluids.new()
    pert.set_stiff_solver(True)
    pert.props.reltol = prec
    pert.set_mode_k(k)
    alpha_try = -cosmo.abs_alpha(1.0e-14 * k**2)
    ci = Ncm.Vector.new(8)
    alphai = pert.get_cross_time(cosmo, Nc.HIPertTwoFluidsCross.MODE1MAIN, alpha_try, cross_size)
    alphaf = pert.get_cross_time(cosmo, Nc.HIPertTwoFluidsCross.MODE1MAIN, alpha_try, 1.0)
    pert.get_init_cond_QP(cosmo, alphai, 1, 0.25 * math.pi, ci)
    pert.set_init_cond(cosmo, alphai, 1, True, ci)

    res = pert.evolve_array(cosmo, min(alphaf, max_time))
    res_a = np.array(res.dup_array()).reshape(-1,9)
    for row in res_a:
        alpha = row[0]
        ci.set_array(row[1:])
        pert.to_zeta_s(cosmo, alpha, ci)
        row[1:] = np.array(ci.dup_array())
    #print(f"A {max_time} {alphaf} {res_a[-1,0]}")

    if max_time > alphaf:
        pert.set_init_cond(cosmo, alphaf, 1, False, ci)
        res_c = pert.evolve_array(cosmo, max_time)
        res_c_a = np.array(res_c.dup_array()).reshape(-1,9)
        res_a = np.concatenate((res_a, res_c_a))

    #del pert
    #gc.collect()

    return res_a, alphaf


def integrate_mode2_qp(k, cosmo, max_time=-1.0):
    # Defining relative tolerance for integration
    prec = 1.0e-11
    # Ratio potential frequency to define cross time
    cross_size = 1.0e-5

    pert = Nc.HIPertTwoFluids.new()
    pert.set_stiff_solver(True)
    pert.props.reltol = prec
    pert.set_mode_k(k)
    alpha_try = -cosmo.abs_alpha(1.0e-14 * k**2)
    ci = Ncm.Vector.new(8)
    alphai = pert.get_cross_time(cosmo, Nc.HIPertTwoFluidsCross.MODE2MAIN, alpha_try, cross_size)
    alphaf = pert.get_cross_time(cosmo, Nc.HIPertTwoFluidsCross.MODE2MAIN, alpha_try, 1.0)
    pert.get_init_cond_QP(cosmo, alphai, 2, 0.25 * math.pi, ci)
    pert.set_init_cond(cosmo, alphai, 2, True, ci)

    res = pert.evolve_array(cosmo, min(alphaf, max_time))
    res_a = np.array(res.dup_array()).reshape(-1,9)
    for row in res_a:
        alpha = row[0]
        ci.set_array(row[1:])
        pert.to_zeta_s(cosmo, alpha, ci)
        row[1:] = np.array(ci.dup_array())
    #print(f"B {max_time} {alphaf} {res_a[-1,0]}")

    if max_time > alphaf:
        pert.set_init_cond(cosmo, alphaf, 2, False, ci)
        res_c = pert.evolve_array(cosmo, max_time)
        res_c_a = np.array(res_c.dup_array()).reshape(-1,9)
        res_a = np.concatenate((res_a, res_c_a))

    #del pert
    #gc.collect()
    return res_a, alphaf



# Definitions of the OLD modes class

class mode1:
  def __init__(self, k, cosmo, max_time=-1):
    mode1 = integrate_mode1(k, cosmo, max_time)
    self.t = mode1[:,0]
    self.zeta_r = mode1[:,zr_index]
    self.s_r = mode1[:,sr_index]
    self.zeta_im = mode1[:,zi_index]
    self.s_im = mode1[:,si_index]

class mode1_qp:
  def __init__(self, k, cosmo, max_time=-1):
    mode1, alpha = integrate_mode1_qp(k, cosmo, max_time)
    self.t = mode1[:,0]
    self.zeta_r = mode1[:,zr_index]
    self.zeta_im = mode1[:,zi_index]
    self.s_r = mode1[:,sr_index]
    self.s_im = mode1[:,si_index]

class mode2:
  def __init__(self, k, cosmo, max_time=-1):
    mode2 = integrate_mode2(k, cosmo, max_time)
    self.t = mode2[:,0]
    self.zeta_r = mode2[:,zr_index]
    self.s_r = mode2[:,sr_index]
    self.zeta_im = mode2[:,zi_index]
    self.s_im = mode2[:,si_index]

class mode2_qp:
  def __init__(self, k, cosmo, max_time=-1):
    mode2, alpha = integrate_mode2_qp(k, cosmo, max_time)
    self.t = mode2[:,0]
    self.zeta_r = mode2[:,zr_index]
    self.zeta_im = mode2[:,zi_index]
    self.s_r = mode2[:,sr_index]
    self.s_im = mode2[:,si_index]

# SANDRO's functions to solve for the modes

def solve_mode1(k, cosmo, max_time=-1.0):
    # Defining relative tolerance for integration
    prec = 1.0e-7
    # Ratio potential frequency to define cross time
    cross_size = 1.0e-13

    pert = Nc.HIPertTwoFluids.new()
    pert.set_stiff_solver(True)
    pert.props.reltol = prec
    pert.set_mode_k(k)
    alpha_try = -cosmo.abs_alpha(1.0e-14 * k**2)
    ci = Ncm.Vector.new(8)
    alphai = pert.get_cross_time(cosmo, Nc.HIPertTwoFluidsCross.MODE1MAIN, alpha_try, cross_size)
    pert.get_init_cond_zetaS(cosmo, alphai, 1, 0.25 * math.pi, ci)
    pert.set_init_cond(cosmo, alphai, 1, False, ci)

    res = pert.evolve_array(cosmo, max_time)

    #del pert
    #gc.collect()
    return res

def solve_mode2(k, cosmo, max_time=-1.0):
    # Defining relative tolerance for integration
    prec = 1.0e-7
    # Ratio potential frequency to define cross time
    cross_size = 1.0e-5

    pert = Nc.HIPertTwoFluids.new()
    pert.set_stiff_solver(True)
    pert.props.reltol = prec
    pert.set_mode_k(k)
    alpha_try = -cosmo.abs_alpha(1.0e-14 * k**2)
    ci = Ncm.Vector.new(8)
    alphai = pert.get_cross_time(cosmo, Nc.HIPertTwoFluidsCross.MODE2MAIN, alpha_try, cross_size)
    pert.get_init_cond_zetaS(cosmo, alphai, 2, 0.25 * math.pi, ci)
    pert.set_init_cond(cosmo, alphai, 2, False, ci)

    res = pert.evolve_array(cosmo, max_time)

    #del pert
    #gc.collect()
    return res

# Demétrio's NEW classes for the modes

# Initial condition 1
class mode_1:
  def __init__(self, k, cosmo, max_time=-1):
    self.NcmMatrix = solve_mode1(k, cosmo, max_time)
    self.t = np.array(self.NcmMatrix.get_col(0).dup_array())
    self.model = cosmo
    self.k = k

  def cross_times(self):
    time_interval = np.geomspace( self.t[0], -1, int( ( self.t[0]/(-1) ) * 10 ) )
    dynamic = model_evol(self.model, time_interval, self.k)
    return [dynamic.cross_z(self.k), dynamic.cross_s(self.k)]

  def zeta_r(self):
    return np.array(self.NcmMatrix.get_col(1).dup_array())

  def s_r(self):
    return np.array(self.NcmMatrix.get_col(2).dup_array())

  def pz_r(self):
    return np.array(self.NcmMatrix.get_col(3).dup_array())

  def ps_r(self):
    return np.array(self.NcmMatrix.get_col(4).dup_array())

  def zeta_im(self):
    return np.array(self.NcmMatrix.get_col(5).dup_array())

  def s_im(self):
    return np.array(self.NcmMatrix.get_col(6).dup_array())

  def pz_im(self):
    return np.array(self.NcmMatrix.get_col(7).dup_array())

  def ps_im(self):
    return np.array(self.NcmMatrix.get_col(8).dup_array())

# Initial condition 2
class mode_2:
  def __init__(self, k, cosmo, max_time=-1):
    self.NcmMatrix = solve_mode2(k, cosmo, max_time)
    self.t = np.array(self.NcmMatrix.get_col(0).dup_array())
    self.model = cosmo
    self.k = k

  def cross_times(self):
    time_interval = np.geomspace( self.t[0], -1, int( ( self.t[0]/(-1) ) * 10 ) )
    dynamic = model_evol(self.model, time_interval, self.k)
    return [dynamic.cross_z(self.k), dynamic.cross_s(self.k)]

  def zeta_r(self):
    return np.array(self.NcmMatrix.get_col(1).dup_array())

  def s_r(self):
    return np.array(self.NcmMatrix.get_col(2).dup_array())

  def pz_r(self):
    return np.array(self.NcmMatrix.get_col(3).dup_array())

  def ps_r(self):
    return np.array(self.NcmMatrix.get_col(4).dup_array())

  def zeta_im(self):
    return np.array(self.NcmMatrix.get_col(5).dup_array())

  def s_im(self):
    return np.array(self.NcmMatrix.get_col(6).dup_array())

  def pz_im(self):
    return np.array(self.NcmMatrix.get_col(7).dup_array())

  def ps_im(self):
    return np.array(self.NcmMatrix.get_col(8).dup_array())

###################################################################################################################################################

# C. Spectra commands

def spec_params(Omegars, w, xb,  E0, k_interval):
    cosmo = Nc.HICosmoQGRW() 
    cosmo.props.w = w
    cosmo.props.Omegar = E0 * Omegars
    cosmo.props.Omegaw = E0 * (1.0 - Omegars)
    cosmo.props.xb = xb

    pert = Nc.HIPertTwoFluids.new()
    pert.props.reltol = 1.0e-9

    k0 = k_interval[0]
    kf = k_interval[-1]
    L = len(k_interval)

    spec1 = pert.compute_zeta_spectrum(cosmo, 1, -cosmo.abs_alpha(1.0e-14), -1.0, k0, kf, L)
    spec2 = pert.compute_zeta_spectrum(cosmo, 2, -cosmo.abs_alpha(1.0e-14), -1.0, k0, kf, L)

    return spec1, spec2

# careful with this function
def entropy_spec_params(Omegars, w, E0, k_interval):
    cosmo.props.w = w
    cosmo.props.Omegar = E0 * Omegars
    cosmo.props.Omegaw = E0 * (1.0 - Omegars)

    pert = Nc.HIPertTwoFluids.new()
    pert.props.reltol = 1.0e-9

    k0 = k_interval[0]
    kf = k_interval[-1]
    L = len(k_interval)

    spec1 = pert.computes(cosmo, 1, -cosmo.abs_alpha(1.0e-14), -1.0, k0, kf, L)
    spec2 = pert.computes(cosmo, 2, -cosmo.abs_alpha(1.0e-14), -1.0, k0, kf, L)

    return spec1, spec2

def spec(cosmo, k_interval, E0 = 1):
  w = cosmo.props.w
  Omegars = cosmo.props.Omegar
  Omegarm  = E0 - cosmo.props.Omegaw
  xb = cosmo.props.xb
  spec1, spec2 = spec_params(Omegars, w, xb,  E0, k_interval)

  Pk_a1 = np.exp( np.array(spec1.peek_yv().dup_array()) )
  Pk_a2 = np.exp( np.array(spec2.peek_yv().dup_array()) )

  return ( Pk_a1 + Pk_a2 )

def specs(cosmo, k_interval, E0 = 1):
  w = cosmo.props.w
  Omegars = cosmo.props.Omegar
  Omegars  = 1.0 - cosmo.props.Omegaw
  xb = cosmo.props.xb
  spec1, spec2 = spec_params(Omegars, w, xb,  E0, k_interval)

  Pk_a1 = np.exp( np.array(spec1.peek_yv().dup_array()) )
  Pk_a2 = np.exp( np.array(spec2.peek_yv().dup_array()) )

  return [ Pk_a1 , Pk_a2 ]

def spec_index(cosmo):
  CMB_range = np.geomspace(1e-3,1e3,1000)
  power_spectrum = spec(cosmo,CMB_range)

  ms, As = np.polyfit(np.log(CMB_range) , np.log(power_spectrum), 1 )
  return ms+1, As

###################################################################################################################################################

# D. Correlation Functions

class spectra:
  def __init__(self, k_interval, cosmo, max_time=-1):
    self.k = k_interval
    self.model = cosmo
    self.tf = max_time

  # Final values of the modes for the Initial Condition 1
  def values1k(self,k):
    mode = mode_1(k,self.model)
    n = mode.NcmMatrix.col_len() - 1
    values = []
    for j in range(1,9):
      values.append(mode.NcmMatrix.get(n,j))
    #del mode
    #gc.collect()
    return values

  # Final values of the modes for the Initial Condition 2
  def values2k(self,k):
    mode = mode_2(k,self.model)
    n = mode.NcmMatrix.col_len() - 1
    values = []
    for j in range(1,9):
      values.append(mode.NcmMatrix.get(n,j))
    #del mode
    #gc.collect()
    return values

  def zz1(self):
    result = []
    for k in self.k:
      values = self.values1k(k)
      result.append( values[0]**2 + values[4]**2 )
    return np.array(result)

  def zz2(self):
    result = []
    for k in self.k:
      values = self.values2k(k)
      result.append( values[0]**2 + values[4]**2 )
    return np.array(result)

  def corr1(self):
    #corrs = []

    P_zz = []
    P_sz = []
    P_pzz = []
    P_psz = []

    P_ss = []
    P_pzs = []
    P_pss = []

    P_pzpz = []
    P_pspz = []

    P_psps = []

    for k in self.k:
      values1 = self.values1k(k)

      zr1 = values1[0]
      sr1 = values1[1]
      pzr1 = values1[2]
      psr1 = values1[3]
      zim1 = values1[4]
      sim1 = values1[5]
      pzim1 = values1[6]
      psim1 = values1[7]

      z1 = complex(zr1,zim1)
      s1 = complex(sr1,sim1)
      pz1 = complex(pzr1,pzim1)
      ps1 = complex(psr1,psim1)

      # z column
      p_zz = z1 * np.conjugate(z1)
      p_sz = s1 * np.conjugate(z1)
      p_pzz = pz1 * np.conjugate(z1)
      p_psz = ps1 * np.conjugate(z1)

      # s column
      p_ss = s1 * np.conjugate(s1)
      p_pzs = pz1 * np.conjugate(s1)
      p_pss = ps1 * np.conjugate(s1)

      # pz column
      p_pzpz = pz1 * np.conjugate(pz1)
      p_pspz = ps1 * np.conjugate(pz1)

      # ps column
      p_psps = ps1 * np.conjugate(ps1)

      P_zz.append(p_zz)
      P_sz.append(p_sz)
      P_pzz.append(p_pzz)
      P_psz.append(p_psz)

      P_ss.append(p_ss)
      P_pzs.append(p_pzs)
      P_pss.append(p_pss)

      P_pzpz.append(p_pzpz)
      P_pspz.append(p_pspz)

      P_psps.append(p_psps)

    corr =  (self.k**3) * np.array( [ np.array(P_zz), np.array(P_sz),  np.array(P_pzz), np.array(P_psz),
                                                      np.array(P_ss), np.array(P_pzs),  np.array(P_pss),
                                                                      np.array(P_pzpz), np.array(P_pspz),
                                                                                        np.array(P_psps)] )
    return corr

  def corr2(self):
    #corrs = []

    P_zz = []
    P_sz = []
    P_pzz = []
    P_psz = []

    P_ss = []
    P_pzs = []
    P_pss = []

    P_pzpz = []
    P_pspz = []

    P_psps = []

    for k in self.k:
      values2 = self.values2k(k)

      zr2 = values2[0]
      sr2 = values2[1]
      pzr2 = values2[2]
      psr2 = values2[3]
      zim2 = values2[4]
      sim2 = values2[5]
      pzim2 = values2[6]
      psim2 = values2[7]

      z2 = complex(zr2,zim2)
      s2 = complex(sr2,sim2)
      pz2 = complex(pzr2,pzim2)
      ps2 = complex(psr2,psim2)

      #result.append( values[0]**2 + values[4]**2 )

      # z column
      p_zz = z2 * np.conjugate(z2)
      p_sz = s2 * np.conjugate(z2)
      p_pzz = pz2 * np.conjugate(z2)
      p_psz = ps2 * np.conjugate(z2)

      # s column
      p_ss = s2 * np.conjugate(s2)
      p_pzs = pz2 * np.conjugate(s2)
      p_pss = ps2 * np.conjugate(s2)

      # pz column
      p_pzpz = pz2 * np.conjugate(pz2)
      p_pspz = ps2 * np.conjugate(pz2)

      # ps column
      p_psps = ps2 * np.conjugate(ps2)

      P_zz.append(p_zz)
      P_sz.append(p_sz)
      P_pzz.append(p_pzz)
      P_psz.append(p_psz)

      P_ss.append(p_ss)
      P_pzs.append(p_pzs)
      P_pss.append(p_pss)

      P_pzpz.append(p_pzpz)
      P_pspz.append(p_pspz)

      P_psps.append(p_psps)

    corr =  (self.k**3) * np.array( [ np.array(P_zz), np.array(P_sz),  np.array(P_pzz), np.array(P_psz),
                                                      np.array(P_ss), np.array(P_pzs),  np.array(P_pss),
                                                                      np.array(P_pzpz), np.array(P_pspz),
                                                                                        np.array(P_psps)] )
    return corr

# Correlation Matrix

class corrmatrix:
  def __init__(self,corr_list):
    self.corrs = corr_list

    # z column
    self.zz = corr_list[0]
    self.sz = corr_list[1]
    self.pzz = corr_list[2]
    self.psz = corr_list[3]

    # s column
    self.zs = np.conjugate( self.sz )
    self.ss = corr_list[4]
    self.pzs = corr_list[5]
    self.pss = corr_list[6]

    # pz column
    self.zpz = np.conjugate( self.pzz )
    self.spz = np.conjugate( self.pzs )
    self.pzpz = corr_list[7]
    self.pspz = corr_list[8]

    # ps column
    self.zps = np.conjugate( self.psz )
    self.sps = np.conjugate( self.pss )
    self.pzps = np.conjugate( self.pspz )
    self.psps = corr_list[9]

# Covariance Matrix

class covmatrix:
  def __init__(self,corr_list):
    self.corrs = corr_list

    # z column
    self.zz = corr_list[0]
    self.sz = np.absolute(corr_list[1]) / (corr_list[0] * corr_list[4])
    self.pzz = np.absolute(corr_list[2]) / ( corr_list[0] * corr_list[7] )
    self.psz = np.absolute(corr_list[3]) / (  corr_list[0] * corr_list[9] )

    # s column
    self.zs = self.sz
    self.ss = corr_list[4]
    self.pzs = np.absolute(corr_list[5]) / (corr_list[4] * corr_list[7])
    self.pss = np.absolute(corr_list[6]) / (corr_list[4] * corr_list[9])

    # pz column
    self.zpz = self.pzz
    self.spz = self.pzs
    self.pzpz = corr_list[7]
    self.pspz = np.absolute(corr_list[8]) / (corr_list[7] * corr_list[9])

    # ps column
    self.zps = self.psz
    self.sps = self.pss
    self.pzps = self.pspz
    self.psps = corr_list[9]
