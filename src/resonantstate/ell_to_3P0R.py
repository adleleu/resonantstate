#Plots a 3-planet 0th degree MMR p*n1 - (p+q)*n2 + q*n3 into the pendulum defined by the Hamiltonian H = 1/(2*delta)*(Sig - delta)^2 + cos(sig)

# Author : Jérémy COUTURIER

# Functions : 
# - ell2pendulum(p, q, m1, m2, m3, T1, T2, T3, lbd1, lbd2, lbd3) -> Returns the coordinates (Sig, sig, delta) of the pendulum Hamiltonian

import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib as mpl
import numpy as np
import celeries as cl
import celeries.perham3pla as h3
import celeries.perham as ph
import celeries.series as se
import celeries.mpfrac as mf
import celeries.laplace as lp
import celeries.ellipseries as es
import celeries.prime as pm
se.setdisplay('expi')
se.setseparator('\n')
ctx = se.mf.mp
ctx.dps = 128 #Guaranteeing correct computation of Xi for large p,q or atypical n1/n2 and n2/n3
se.mf.setctx(ctx)


def ell_to_3P0R(p, q, m1, m2, m3, T1, T2, T3, lbd1, lbd2, lbd3, Xi = 0., verbose=True):
      # Converts from elliptic elements to the coordinates (Sig, sig, delta) of the pendulum of a three-planet 0th degree resonance
      # p and q are such that the triplet is close from MMR p*n1 - (p+q)*n2 + q*n3 = 0
      # The planetary masses are relative to the star's mass
      # The lbd_j are in radians
      # The T_j can be given in any units
      # If Xi is given in argument, it is not recomputed. If given, can be an array
      # All elliptic elements given in argument can be 1D arrays
      # Returns [Sig, sig, delta, (m1 + m2 + m3)**(4./3.)*delta, Xi] such that the Hamiltonian is 1/(2*delta)*(Sig-delta)**2 + cos(sig). See Sect. 3 of Couturier et al. (2026) for details
      
      #Trivial validity check
      if isinstance(T2, np.ndarray):
            if (np.median(T1) > np.median(T2) or np.median(T2) > np.median(T3)):
                  raise Exception("      Error: The periods must verify T1 < T2 < T3")
      else:
            if (T1 > T2 or T2 > T3):
                  raise Exception("      Error: The periods must verify T1 < T2 < T3")
            
      if (p <= 0 or q <= 0):
            raise Exception("      Error: p and q must both be strictly positive integers")

      # Period of inner planet is normalized to 1
      T3 = T3/T1
      T2 = T2/T1
      T1 = 1.

      # Getting semi-major axes and Lambda
      G     = 4.*np.pi**2
      beta1 = m1/(1. + m1)
      beta2 = m2/(1. + m2)
      beta3 = m3/(1. + m3)
      mu1   = G*(1. + m1)
      mu2   = G*(1. + m2)
      mu3   = G*(1. + m3)
      n1    = 2.*np.pi
      n2    = 2.*np.pi/T2
      n3    = 2.*np.pi/T3      
      lambda_over_two = (p*n1 - (p+q)*n2 + q*n3)/(p**2 + q**2 + (p+q)**2)
      n10   = n1 - p*lambda_over_two #n10, n20 and n30 are chosen by minimizing (n1 - n10)**2 + (n2 - n20)**2 + (n3 - n30)**2, and by forcing p*n10 - (p+q)*n20 + q*n30 = 0
      n20   = n2 + (p+q)*lambda_over_two
      n30   = n3 - q*lambda_over_two
      a1    = (mu1/n1**2)**(1./3.)
      a2    = (mu2/n2**2)**(1./3.)
      a3    = (mu3/n3**2)**(1./3.)
      Lbd1  = beta1*np.sqrt(mu1*a1)
      Lbd2  = beta2*np.sqrt(mu2*a2)
      Lbd3  = beta3*np.sqrt(mu3*a3)

      #Defining the exact resonance
      a10   = (mu1/n10**2)**(1./3.)
      a20   = (mu2/n20**2)**(1./3.)
      a30   = (mu3/n30**2)**(1./3.)
      Lbd10 = beta1*np.sqrt(mu1*a10)
      Lbd20 = beta2*np.sqrt(mu2*a20)
      Lbd30 = beta3*np.sqrt(mu3*a30)
      
      #Normalizing
      Gamma0 = q*(p + 1)/p*Lbd10 + q*Lbd20 + (q - 1)*Lbd30
      Gamma = (q*(p + 1)/p*Lbd1 + q*Lbd2 + (q - 1)*Lbd3)/Gamma0
      Upsilon = (Lbd1 + Lbd2 + Lbd3)/Gamma0
      Phi = (Lbd1/p)/Gamma0
      
      #Getting resonance coefficients
      Lbd10 = Lbd10/Gamma0
      Lbd20 = Lbd20/Gamma0
      Lbd30 = Lbd30/Gamma0
      K = n10*p**2/Lbd10 + n20*(p+q)**2/Lbd20 + n30*q**2/Lbd30
      Phi_eq = ((p+q)*n20/Lbd20*(Gamma - (q - 1)*Upsilon) + q*n30/Lbd30*(Gamma - q*Upsilon))/K
      if (Xi == 0.):
          if (isinstance(lambda_over_two, np.ndarray) and np.std(n1/n2) > .0002):
                if (verbose):
                      print("      Standard deviation of the mean motions is large. Xi will be computed for each point.")
                H3 = h3.PerHam3pla(degree = 0, keplerian=True)
                Xis = -2.*H3.angle((p, -p-q, q))
                Xi = np.zeros(len(lambda_over_two))
                for i in range(len(lambda_over_two)):
                      if (i % (len(lambda_over_two)//100) == 0 and verbose):
                            print('\r', f'      progress = {(i*100)//len(lambda_over_two)}%', end='')
                      Xi[i] = H3.eval((p, -p-q, q), Xis, (n10[i], n20[i], n30[i])).toConst()
                print()
          else:
                if (isinstance(lambda_over_two, np.ndarray)):
                      H3 = h3.PerHam3pla(degree = 0, n0 = (np.median(n10), np.median(n20), np.median(n30)), ev = True, keplerian=True)
                else:
                      H3 = h3.PerHam3pla(degree = 0, n0 = (n10, n20, n30), ev = True, keplerian=True)
                Xi = -2.*H3.angle((p, -p-q, q))
                Xi = Xi.toConst()
      Rpq = n30*Xi*Lbd30

      #Getting delta, Sig and sig
      delta = 3.*K*Phi_eq**2/(m1*m2*Rpq)
      Sig = delta*Phi/Phi_eq
      sig = np.mod(p*lbd1 - (p+q)*lbd2 + q*lbd3, 2.*np.pi)
      if isinstance(sig, np.ndarray):
            for i in range(len(sig)):
                  if (sig[i] < 0):
                        sig[i] += 2.*np.pi
      else:
            if (sig < 0):
                  sig += 2.*np.pi
      return [Sig, sig, delta, (m1 + m2 + m3)**(4./3.)*delta, Xi]


def plot_topology(ax1, linewidth=4, alpha=1, grid=False):

      ### Plots the topology of the phase space (separatrices and fixed points) of the pendulum on the axis ax1 ###

      
      ax1.autoscale(axis='both')
      X_lim = ax1.get_xlim()
      Y_lim = ax1.get_ylim()

      ax1.tick_params(axis='both', which='major', labelsize=16)
      #ax1.set_xlabel(xlabel=r"$\varepsilon^{4/3}\delta$", fontsize=16)
      #ax1.set_ylabel(ylabel=r"$\left(\Sigma-\delta\right)/\sqrt{\delta}$", fontsize=16)

      #Plotting the separatrix
      ax1.hlines(y = 2., xmin = 5000., xmax = 1.e8, colors='red', linestyles='solid', data=None, linewidth = 3, alpha = 1, label = 'Separatrix')
      ax1.hlines(y = -2.,xmin = 5000., xmax = 1.e8, colors='red', linestyles='solid', data=None, linewidth = 3, alpha = 1)
      ax1.hlines(y = 2., xmin = -1.e8, xmax = -5000, colors='red', linestyles='solid', data=None, linewidth = 3, alpha = 1)
      ax1.hlines(y = -2.,xmin = -1.e8, xmax = -5000, colors='red', linestyles='solid', data=None, linewidth = 3, alpha = 1)
      #ax1.hlines(y = 2., xmin = 1000., xmax = X_lim[1], colors='red', linestyles='solid', data=None, linewidth = 3, alpha = 1, label = 'Separatrix')
      #ax1.hlines(y =-2., xmin = 1000., xmax = X_lim[1], colors='red', linestyles='solid', data=None, linewidth = 3, alpha = 1)
      for x in np.linspace(100., 5000., 4901):
            alpha = (x - 100.)/4900.
            ax1.hlines(y = 2., xmin = x, xmax = x+1., colors='red', linestyles='solid', data=None, linewidth = 3, alpha = alpha)
            ax1.hlines(y =-2., xmin = x, xmax = x+1., colors='red', linestyles='solid', data=None, linewidth = 3, alpha = alpha)
      for x in np.linspace(-5000., -100., 4901):
            alpha = (-x - 100.)/4900.
            ax1.hlines(y = 2., xmin = x, xmax = x+1., colors='red', linestyles='solid', data=None, linewidth = 3, alpha = alpha)
            ax1.hlines(y =-2., xmin = x, xmax = x+1., colors='red', linestyles='solid', data=None, linewidth = 3, alpha = alpha)


      if grid:
            ax1.grid(linewidth=0.3, alpha = 0.5)


