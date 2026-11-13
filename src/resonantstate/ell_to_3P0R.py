#Plots a 3-planet 0th degree MMR p*n1 - (p+q)*n2 + q*n3 into the pendulum defined by the Hamiltonian H = Sig**2/2 + sgn(Xi) * cos(sig)

# Author : Jérémy COUTURIER

# Functions : 
# - ell_to_3P0R(p, q, m1, m2, m3, T1, T2, T3, lbd1, lbd2, lbd3, Xi = 0., verbose=False) -> Returns the coordinates (Sig, sig, D3P0R, D2P1R, Xi) of the pendulum Hamiltonian.

# If you use this file, cite Couturier et al. (2027)

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


def Lbdj0(Lbdj, wj, cj, prec = 1.e-13):
      #Returns Lbdj0 = [\Lambda_{1,0}, \Lambda_{2,0}, \Lambda_{3,0}], from Lbdj = [\Lambda_1, \Lambda_2, \Lambda_3], wj = [p, -p-q, q] and cj = \beta_j^3\mu_j^2
      #The returned \Lambda_j,0 comply with pn_1,0 - (p+q)n_2,0 + qn_3,0=0, as well as \Delta\Gamma = 0 and \Delta\Upsilon = 0 (Couturier et al. 2026, Sect. 3)
      #Denoting \Delta\Lambda_j = \Lambda_j - \Lambda_j,0, the solution reads G(s) = 0 with \Delta\Lambda_j = w_j*s and G(s) = \sum_j w_j*c_j/(\Lambda_j - w_j*s)**3
      #Unfortunately, G(s) is a sextic in s and no analytical expression of its roots can be given. Therefore, I find s with a Newton-Raphson method
      #Since s belongs to the range [Lbd2/w_2, min(Lbd1/w_1, Lbd3/w_3)], I check after each step of the Newton-Raphson whether the solution is in the range.
      #G(s) is an increasing function of s in that range, so at each new evaluation of G(s), I keep track of a range [low, high] inside which s must be.
      #If an iteration gets out of range, I restart from the middle of [low, high].
      
      if not isinstance(Lbdj, np.ndarray):
            Lbdj = np.array(Lbdj)
      if not isinstance(wj, np.ndarray):
            wj = np.array(wj)
      if not isinstance(cj, np.ndarray):
            cj = np.array(cj)
      
      low = Lbdj[1]/wj[1]
      high = min(Lbdj[0]/wj[0], Lbdj[2]/wj[2])
      prec = prec*(high - low)
      s0 = 0.
      step = 1.e300
      it = 0
      while(abs(step) > prec):
            Gs = np.sum(wj*cj/(Lbdj-wj*s0)**3)
            if (Gs > 0):
                  high = s0
            else:
                  low = s0
            dGs = np.sum(3.*wj**2*cj/(Lbdj-wj*s0)**4)
            step = -Gs/dGs
            s = s0 + step
            s0 = s
            if (s > high or s < low):
                  s0 = .5*(low+high)
                  step = 1.e300
            it = it + 1
            if (it > 10000):
                  raise Exception("Newton-Raphson method does not converge in function Lbdj0")
      return Lbdj - wj*s


def ell_to_3P0R(p, q, m1, m2, m3, T1, T2, T3, lbd1, lbd2, lbd3, Xi = 0., verbose=False):
      # Converts from elliptic elements to the coordinates (Sig, sig, D3P0R, D2P1R) of the pendulum of a three-planet 0th degree resonance. See Couturier et al. 2026, Sects 2 and 3
      # p and q are such that the triplet is close from MMR p*n1 - (p+q)*n2 + q*n3 = 0
      # The planetary masses are relative to the star's mass
      # The lbd_j are in radians
      # The T_j can be given in any units
      # If Xi is given in argument, it is not recomputed. If given, can be an array
      # All elliptic elements given in argument can be 1D arrays
      # Returns [Sig, sig, D3P0R, D2P1R, Xi] such that the Hamiltonian is Sig**2/2 + sgn(Xi) * cos(sig).
      
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
      a1    = (mu1/n1**2)**(1./3.)
      a2    = (mu2/n2**2)**(1./3.)
      a3    = (mu3/n3**2)**(1./3.)
      Lbd1  = beta1*np.sqrt(mu1*a1)
      Lbd2  = beta2*np.sqrt(mu2*a2)
      Lbd3  = beta3*np.sqrt(mu3*a3)
      
      #Defining the exact resonance such that pn10 - (p+q)n20 + qn30 = 0, \Delta\Gamma = 0, \Delta\Upsilon = 0
      wj = np.array([p, -p-q, q])
      if isinstance(Lbd2, np.ndarray):
            N = len(Lbd2)
            Lbd10 = np.zeros(N)
            Lbd20 = np.zeros(N)
            Lbd30 = np.zeros(N)
            for i in range(N):
                  Lbd_j0 = Lbdj0(np.array([Lbd1[i], Lbd2[i], Lbd3[i]]), wj, np.array([beta1[i]**3*mu1[i]**2, beta2[i]**3*mu2[i]**2, beta3[i]**3*mu3[i]**2]))
                  Lbd10[i] = Lbd_j0[0]
                  Lbd20[i] = Lbd_j0[1]
                  Lbd30[i] = Lbd_j0[2]
      else:
            Lbd_j0 = Lbdj0(np.array([Lbd1, Lbd2, Lbd3]), wj, np.array([beta1**3*mu1**2, beta2**3*mu2**2, beta3**3*mu3**2]))
            Lbd10 = Lbd_j0[0]
            Lbd20 = Lbd_j0[1]
            Lbd30 = Lbd_j0[2]      
      a10 = Lbd10**2/(mu1*beta1**2)
      a20 = Lbd20**2/(mu2*beta2**2)
      a30 = Lbd30**2/(mu3*beta3**2)
      n10 = np.sqrt(mu1/a10**3)
      n20 = np.sqrt(mu2/a20**3)
      n30 = np.sqrt(mu3/a30**3)     
      
      #Computing Xi if not given in argument
      if (Xi == 0.):
          if (isinstance(n20, np.ndarray) and np.std(n1/n2) > .0002):
                if (verbose):
                      print("      Standard deviation of the mean motions is large. Xi will be computed for each point.")
                H3 = h3.PerHam3pla(degree = 0, keplerian=True)
                Xis = -2*H3.angle((p, -p-q, q))
                Xi = np.zeros(len(n20))
                for i in range(len(n20)):
                      if (i % (len(n20)//100) == 0 and verbose):
                            print('\r', f'      progress = {(i*100)//len(n20)}%', end='')
                      Xi[i] = H3.eval((p, -p-q, q), Xis, (n10[i], n20[i], n30[i]), dps=128).toConst()
                print()
          else:
                if (isinstance(n20, np.ndarray)):
                      H3 = h3.PerHam3pla(degree = 0, n0 = (np.median(n10), np.median(n20), np.median(n30)), ev = True, keplerian=True, dps=128)
                else:
                      H3 = h3.PerHam3pla(degree = 0, n0 = (n10, n20, n30), ev = True, keplerian=True, dps=128)
                Xi = -2.*H3.angle((p, -p-q, q))
                Xi = Xi.toConst()
      
      
      #Computing Sig and sig
      sgnXi = np.abs(Xi)/Xi
      DLbd1 = Lbd1 - Lbd10
      DLbd2 = Lbd2 - Lbd20
      DLbd3 = Lbd3 - Lbd30
      DPhi = DLbd1/p
      DGamma = (p+q)/p*DLbd1 + DLbd2 #Should be zero
      DUpsilon = DLbd1 + DLbd2 + DLbd3 #Should be zero
      if (np.max(np.abs(DGamma)) > 1.e-9 or np.max(np.abs(DUpsilon)) > 1.e-9 or np.max(np.abs(p*n10/n20-(p+q)+q*n30/n20)) > 1.e-9):
            print('(p*n10-(p+q)*n20+q*n30)/n20 =', np.max(np.abs(p*n10/n20-(p+q)+q*n30/n20)))
            print('DGamma =', np.max(np.abs(DGamma)))
            print('DUpsilon =', np.max(np.abs(DUpsilon)))
            raise Exception("A problem occured while computing the n_j,0 in function ell_to_3P0R")
      K = n10*p**2/Lbd10 + n20*(p+q)**2/Lbd20 + n30*q**2/Lbd30
      omega = sgnXi*m1*m2*n30*Lbd30*Xi
      B = np.sqrt(omega/(3.*K))
      Sig = DPhi/B
      sig = np.mod(p*lbd1 - (p+q)*lbd2 + q*lbd3, 2.*np.pi)
      if isinstance(sig, np.ndarray):
            for i in range(len(sig)):
                  if (sig[i] < 0):
                        sig[i] += 2.*np.pi
      else:
            if (sig < 0):
                  sig += 2.*np.pi
      
      #Computing D3P0R
      H = .5*Sig**2 + sgnXi*np.cos(sig)
      D3P0R = np.sqrt(2. + 2.*H)
      if isinstance(H, np.ndarray):
            for i in range(len(H)):
                  if (H[i] > 1. and Sig[i] < 0.):
                        D3P0R[i] = -D3P0R[i]
      else:
            if (H > 1. and Sig < 0.):
                  D3P0R = -D3P0R
      
      #Computing D2P1R
      D2P1R = (m1+m2+m3)**(-2./3.)*((p/(p+1)*n10/n20)**(1./3.) - 1.)

      return [Sig, sig, D3P0R, D2P1R, Xi]


def plot_topology(ax1, linewidth=4, alpha=1, grid=False):

      ### Plots the topology of the phase space (separatrices and fixed points) of the pendulum on the axis ax1 ###

      
      ax1.autoscale(axis='both')
      X_lim = ax1.get_xlim()
      Y_lim = ax1.get_ylim()

      ax1.tick_params(axis='both', which='major', labelsize=16)
      #ax1.set_xlabel(xlabel=r"$\varepsilon^{4/3}\delta$", fontsize=16)
      #ax1.set_ylabel(ylabel=r"$\left(\Sigma-\delta\right)/\sqrt{\delta}$", fontsize=16)

      #Plotting the separatrix
      ax1.hlines(y = 2., xmin = 2.,    xmax = 1.e5, colors='red', linestyles='solid', data=None, linewidth = 3, alpha = 1, label = 'Separatrix')
      ax1.hlines(y = -2.,xmin = 2.,    xmax = 1.e5, colors='red', linestyles='solid', data=None, linewidth = 3, alpha = 1)
      ax1.hlines(y = 2., xmin = -1.e5, xmax = -2.,  colors='red', linestyles='solid', data=None, linewidth = 3, alpha = 1)
      ax1.hlines(y = -2.,xmin = -1.e5, xmax = -2.,  colors='red', linestyles='solid', data=None, linewidth = 3, alpha = 1)
      #ax1.hlines(y = 2., xmin = 1000., xmax = X_lim[1], colors='red', linestyles='solid', data=None, linewidth = 3, alpha = 1, label = 'Separatrix')
      #ax1.hlines(y =-2., xmin = 1000., xmax = X_lim[1], colors='red', linestyles='solid', data=None, linewidth = 3, alpha = 1)
      for x in np.linspace(1., 2., 1001):
            alpha = 1000.*(x - 1.)/1000.
            ax1.hlines(y = 2., xmin = x, xmax = x+.001, colors='red', linestyles='solid', data=None, linewidth = 3, alpha = alpha)
            ax1.hlines(y =-2., xmin = x, xmax = x+.001, colors='red', linestyles='solid', data=None, linewidth = 3, alpha = alpha)
      for x in np.linspace(-2., -1., 1001):
            alpha = -1000.*(x + 1.)/1000.
            ax1.hlines(y = 2., xmin = x, xmax = x+.001, colors='red', linestyles='solid', data=None, linewidth = 3, alpha = alpha)
            ax1.hlines(y =-2., xmin = x, xmax = x+.001, colors='red', linestyles='solid', data=None, linewidth = 3, alpha = alpha)


      if grid:
            ax1.grid(linewidth=0.3, alpha = 0.5)


