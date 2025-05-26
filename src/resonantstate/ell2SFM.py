# This notebook manipulates 1-dimensional numpy arrays of elliptic elements to convert and plot into the Second Fundamental Model of resonance (SFM).
# The definition of the SFM used here is detailed at https://jeremycouturier.com/img/SFM.pdf
# The numpy arrays can come from posterior data of MCMC analysis or from numerical simulations.
# The user can (optionally) import a sample from the DACE table of the Geneva Resonant State Workshop (GRSW) for analysis

# Author : Jérémy COUTURIER

# Functions : 
# - ell2SFM(p, e1, e2, vp1, vp2, m1, m2, T1, T2, lbd1, lbd2) -> Returns the coordinates (X, Y, X2, Y2, delta) of the SFM. (X, Y) corresponds to the unique degree of freedom
#                                                               of the SFM and (X2, Y2) to the first integral. delta is the unique parameter of the SFM
# - SFM2useful(X, Y, X2, Y2, delta) -> Returns [sig, Sig, sig2, Sig2, x1, x2, IsResonant] where X+iY = sqrt(2*Sig)*e^(i*sig) and X2+iY2 = sqrt(2*Sig2)*e^(i*sig2)
#                                      x1 and x2 are such that (x1, 0) and (x2, 0) are on the same level line as (X, Y)
#                                      IsResonant is 1 if the system is in the resonance, and 0 else.

# import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib as mpl
import numpy as np
from itertools import cycle

import pandas as pd

from resonantstate.simulations_resonance_analysis import *


# import os 

#Hamiltonian coefficients
f1s = [ 1.190493697849547, 2.025222689938346, 2.840431856715441, 3.649618244089652, 4.456142785027623, 5.261253831849899, 6.065523627097718, 6.869251916417852, 7.672611034475267, 8.475707148201764, 9.278609253466129, 10.08136413618922, 10.88400465063751, 11.68655455857515, 12.48903144896030, 13.29144865274429, 14.09381638467312, 14.89614276587963, 15.69843412935734, 16.50069558620453]
f2s = [-0.428389834143869,-2.484005183303907,-3.283256721821090,-4.083705371769611,-4.884706297511002,-5.686007411626633,-6.487489727907814,-7.289089771453291,-8.090770598035306,-8.892509248107672,-9.694290707819164,-10.49610474146903,-11.29794413782656,-12.09980367496610,-12.90167944505811,-13.70356853306293,-14.50546862185001,-15.30737796425819,-16.10929512977600,-16.91121894121170]


def ell2SFM(p, e1, e2, vp1, vp2, m1, m2, T1, T2, lbd1, lbd2):
      # Converts from elliptic elements to the coordinates (X, Y, X2, Y2, delta) of the Second Fondamental Model of resonance (SFM).
      # (X, Y) corresponds to the unique degree of freedom of the SFM and (X2, Y2) to the first integral. delta is the unique parameter of the SFM

      # Assumptions : e1, e2, vp1, vp2, m1, m2, T1, T2, lbd1, lbd2 are 1-dimensional numpy arrays. (e.g. outputs of a MCMC posterior analysis)
      # Subscript 1 (resp. 2) is for the inner (resp. outer) planet
      # The masses m1 and m2 are relative to the star. e1 and e2 are the eccentricities.
      # vp1 and vp2 are the longitudes of the periapsis, lbd1 and lbd2 are the mean longitudes, in radians
      # T1 and T2 are the periods. They can be given in any unit because they will be renormalized by the inner period.
      # p is an integer such that the pair is close to the resonance p : p + 1

      if (p > 20 or p < 1):
            raise Exception('The index p of the resonance must be between 1 and 20 included')

      # Period of inner planet is normalized to 1
      T2    = T2/T1
      T1    = T1/T1

      # Getting semi-major axes and Lambda
      G     = 4.*np.pi**2
      beta1 = m1/(1. + m1)
      beta2 = m2/(1. + m2)
      mu1   = G*(1. + m1)
      mu2   = G*(1. + m2) 
      n1    = 2.*np.pi/T1
      n2    = 2.*np.pi/T2
      n10   = 2.*np.pi
      n20   = p*n10/(p + 1)
      a1    = (mu1/n1**2)**(1./3.)
      a2    = (mu2/n2**2)**(1./3.)
      Lbd1  = beta1*np.sqrt(mu1*a1)
      Lbd2  = beta2*np.sqrt(mu2*a2)
      D1    = Lbd1*(1. - np.sqrt(1. - e1**2))
      D2    = Lbd2*(1. - np.sqrt(1. - e2**2))

      #Defining the exact resonance
      a10   = (mu1/n10**2)**(1./3.)
      a20   = (mu2/n20**2)**(1./3.)
      Lbd10 = beta1*np.sqrt(mu1*a10)
      Lbd20 = beta2*np.sqrt(mu2*a20)

      # Getting G and Gamma and normalizing
      G     = Lbd1 + Lbd2 - D1 - D2
      Gamma = (p + 1)*Lbd1 + p*Lbd2
      g     = G/Gamma
      d1    = D1/Gamma
      d2    = D2/Gamma
      C1    = Gamma/Lbd10
      C2    = Gamma/Lbd20

      #Getting alpha, beta, gamma, delta, R and S
      f1    = f1s[p - 1]
      f2    = f2s[p - 1]
      R     = (f1**2*C1*d1 + f2**2*C2*d2 + 2.*f1*f2*np.sqrt(C1*d1*C2*d2)*np.cos(vp1 - vp2))/(f1**2*C1 + f2**2*C2)
      S     = (f1**2*C1*d2 + f2**2*C2*d1 - 2.*f1*f2*np.sqrt(C1*d1*C2*d2)*np.cos(vp1 - vp2))/(f1**2*C1 + f2**2*C2)
      alpha = -3.*n10*p*((g + S)*(p*C1 + (p + 1)*C2) - C1 - C2)
      beta  = 1.5*n10*p*(p*C1 + (p + 1)*C2)
      gamma = m1*n20/C2*np.sqrt(f1**2*C1 + f2**2*C2)
      delta = alpha*(4./(27.*beta*gamma**2))**(1./3.)

      #Getting X and Y
      K     = (2.*beta/gamma)**(-2./3.)
      omega = beta*(2.*beta/gamma)**(-4./3.)
      Sigma = R/K
      Sigma2= S/K
      xi    = -p*lbd1 + (p + 1)*lbd2
      sig1  = xi - vp1
      sig2  = xi - vp2
      u1    = np.sqrt(2.*d1)*np.cos(sig1)
      u2    = np.sqrt(2.*d2)*np.cos(sig2)
      v1    = np.sqrt(2.*d1)*np.sin(sig1)
      v2    = np.sqrt(2.*d2)*np.sin(sig2)
      z     = f2*np.sqrt(C2)/(f1*np.sqrt(C1))
      cophi = 1./np.sqrt(1. + z**2)
      siphi = z /np.sqrt(1. + z**2)
      x1    = cophi*u1 + siphi*u2
      y1    = cophi*v1 + siphi*v2
      x2    = cophi*u2 - siphi*u1
      y2    = cophi*v2 - siphi*v1
      cossig= x1/np.sqrt(2.*R)
      sinsig= y1/np.sqrt(2.*R)
      cosig2= x2/np.sqrt(2.*S)
      sisig2= y2/np.sqrt(2.*S)
      X     = np.sqrt(2.*Sigma)*cossig
      Y     = np.sqrt(2.*Sigma)*sinsig
      X2    = np.sqrt(2.*Sigma2)*cosig2
      Y2    = np.sqrt(2.*Sigma2)*sisig2
      return [X, Y, X2, Y2, delta]

def cubic(a_3, a_2, a_1, a_0):
      # Returns the real roots of a_3*X^3 + a_2*X^2 + a_1*X + a_0 = 0 using analytical expressions of Cardan's method
      if (a_3 == 0.):
            if (a_2 == 0.):
                  if (a_1 == 0.):
                        return []
                  return [-a_0/a_1]
            Delta = a_1**2 - 4.*a_2*a_0
            if (Delta >= 0.):
                  return [(-a_1 + np.sqrt(Delta))/(2.*a_2), (-a_1 - np.sqrt(Delta))/(2.*a_2)]
            return []
      if (a_3 != 1.):
            return cubic(1., a_2/a_3, a_1/a_3, a_0/a_3)
      ### Equation is Y^3 + p*Y + q = 0 with Y = X + s ###
      p = a_1 - a_2**2/3.
      q = a_0 - a_1*a_2/3. + 2.*a_2**3/27.
      s = a_2/3.
      D = q**2/4. + p**3/27.
      if (D < 0.): # 3 real solutions
            u3 = -q/2. + cm.sqrt(D)
            v3 = -q/2. - cm.sqrt(D)
            [mod_u3, arg_u3] = cm.polar(u3)
            [mod_v3, arg_v3] = cm.polar(v3)
            u  = mod_u3**(1./3.)*cm.exp(1j*arg_u3/3.)
            v  = mod_v3**(1./3.)*cm.exp(1j*arg_v3/3.)
            j  = cm.exp( 2.*1j*np.pi/3.)
            jb = cm.exp(-2.*1j*np.pi/3.)
            S1 = (u + v).real
            S2 = (j*u + jb*v).real
            S3 = (jb*u + j*v).real
            return [S1 - s, S2 - s, S3 - s]
      else: # 1 real solution
            u3 = -q/2. + np.sqrt(D)
            v3 = -q/2. - np.sqrt(D)
            if (u3 < 0.):
                  u = -(-u3)**(1./3.)
            else:
                  u  = u3**(1./3.)
            if (v3 < 0.):
                  v = -(-v3)**(1./3.)
            else:
                  v  = v3**(1./3.)
            return [u + v - s]

def quartic(a_4, a_3, a_2, a_1, a_0):
      # Returns the real roots of a_4*X^4 + a_3*X^3 + a_2*X^2 + a_1*X + a_0 = 0 using analytical expressions of Ferrari's method
      if (a_4 == 0.):
            return cubic(a_3, a_2, a_1, a_0)
      if (a_4 != 1.):
            return quartic(1., a_3/a_4, a_2/a_4, a_1/a_4, a_0/a_4)
      ### Equation is Y^4 + p*Y^2 + q*Y + r = 0 with Y = X + s ###
      p = a_2 - 3./8.*a_3**2
      q = a_1 + a_3**3/8. - 0.5*a_2*a_3
      r = a_0 + a_2*a_3**2/16. - a_1*a_3/4. - 3.*a_3**4/256.
      s = a_3/4.
      ### Taking care of bi-quartic case ###
      if (abs(q) < 1.e-14):
            Sol = cubic(0., 1., p, r)
            if (len(Sol) == 0):
                  return []
            [S1, S2] = Sol
            Sol = []
            # First pair
            if (S1 >= -1.e-14):
                  Sol.append( np.sqrt(abs(S1)) - s)
                  Sol.append(-np.sqrt(abs(S1)) - s)
            if (S2 >= -1.e-14):
                  Sol.append( np.sqrt(abs(S2)) - s)
                  Sol.append(-np.sqrt(abs(S2)) - s)
            return Sol
      ### Getting a solution of the resolving cubic ###
      if (abs(q) < 1.e-6): #Too close from bi-quartic. Must be done differently
            if (abs(2.*p**2 - 8.*r) > 1.e-10):
                  if (2.*p**2 - 8.*r < 0.):
                        return []
                  sqM = abs(q)/np.sqrt(2.*p**2 - 8.*r)
            else:
                  if (abs(8.*p) > 1.e-10):
                        if (p < 0.):
                              sqM = np.sqrt(-p)
                        else:
                              sqM = abs(q)/np.sqrt(8.*p)
                  else:
                        sqM = (q**2/8.)**(1./6.)
      else:
            Sol = cubic(1., p, p**2/4. - r, -q**2/8.)
            Sol.sort()
            if (Sol[-1] < 0.): #There are no real solutions
                  return []
            sqM = np.sqrt(Sol[-1])
      ### Getting first pair of solutions ###
      Sol = []
      D = q/(2.*np.sqrt(2.)*sqM)-(p + sqM**2)/2.
      if (D >= -1.e-14):
            S1 = -sqM/np.sqrt(2.) + np.sqrt(abs(D))
            S2 = -sqM/np.sqrt(2.) - np.sqrt(abs(D))
            Sol.append(S1 - s)
            Sol.append(S2 - s)
      ### Getting second pair of solutions ###
      D = -q/(2.*np.sqrt(2.)*sqM)-(p + sqM**2)/2.
      if (D >= -1.e-14):
            S3 = sqM/np.sqrt(2.) + np.sqrt(abs(D))
            S4 = sqM/np.sqrt(2.) - np.sqrt(abs(D))
            Sol.append(S3 - s)
            Sol.append(S4 - s)
      return Sol

def X1X2(X, Y, delta):
      # Returns X1 and X2 such that (X1, 0) and (X2, 0) are on the same level line as (X, Y)
      H   = 1.5*delta*(X**2 + Y**2) - 0.25*(X**2 + Y**2)**2 + 2.*X
      Sol = quartic(-0.25, 0., 1.5*delta, 2., -H)
      Sol.sort()
      if (len(Sol) == 0):
            print("Warning: Problem with quartic.")
            return [0., 0.]
      if (len(Sol) == 2):
            return Sol
      #Four solutions. Either [Sol[0], Sol[3]] or [Sol[1], Sol[2]] should be returned
      if (delta < 1.): #Four solutions should be impossible when delta < 1
            print("Warning: Four solutions were found even though delta < 1 in function X1X2")
      topo = topology_light(delta)
      if (len(topo) == 1):
            print("Warning: Could not find xhyp and xint in function X1X2")
            return [Sol[0], Sol[3]]
      [xres, xint, xhyp] = topo
      # A very simple criterion proposed by Max Goldberg to determine which pair of solution should be returned 
      if ((X - xint)**2 + Y**2 > (xhyp - xint)**2):
            return [Sol[0], Sol[3]]
      else:
            return [Sol[1], Sol[2]]

def SFM2useful(X, Y, X2, Y2, delta):
      # Returns [sig, Sig, sig2, Sig2, x1, x2] where X+iY = sqrt(2*Sig)*e^(i*sig) and X2+iY2 = sqrt(2*Sig2)*e^(i*sig2)
      # x1 and x2 are such that (x1, 0) and (x2, 0) are on the same level line as (X, Y)
      
      Sig  = (X**2  + Y**2) /2.
      Sig2 = (X2**2 + Y2**2)/2.
      sig  = np.arctan2(Y,  X)
      sig2 = np.arctan2(Y2, X2)
      
      x1 = []
      x2 = []
      IR = []
      
      n = len(delta)
      for i in range(n):
            [xx1, xx2] = X1X2(X[i], Y[i], delta[i])
            x1.append(xx1)
            x2.append(xx2)
            [xres, xint, xhyp] = topology_light(delta[i])
            if (delta[i] < 1.):
                  IR.append(0)
            else:
                  Hseparatrix = 1.5*delta[i]*xhyp**2 - 0.25*xhyp**4 + 2.*xhyp
                  H           = 1.5*delta[i]*(X[i]**2 + Y[i]**2) - 0.25*(X[i]**2 + Y[i]**2)**2 + 2.*X[i]
                  if (H > Hseparatrix):
                        IR.append(1)
                  else:
                        IR.append(0)
      x1 = np.array(x1)
      x2 = np.array(x2)
      IR = np.array(IR)
      return [sig, Sig, sig2, Sig2, x1, x2, IR]

def topology(delta):
      #Returns [Xmin, Xmax, Xres, Xint, Xhyp] from analytical expressions instead of reading from file
      #Xmin and Xmax are the lower and upper separatrices respectively, Xres is the resonance center, Xint is the center of the internal circulation and Xhyp is the hyperbolic fixed point.
      #If delta < 1, returns [0, 0, Xres, 0, 0]
      #When delta >= 1; then Xhyp <= Xint <= Xmin <= Xres <= Xmax
      if (delta == 1.):
            return [-1., 3., 2., -1., -1.]
      Sol = cubic(1., 0., -3.*delta, -2.) #Getting Xres, Xint and Xhyp
      if (len(Sol) == 1):
            return [0., 0., Sol[0], 0., 0.]
      else:
            [S1, S2, S3] = Sol
            if (S2**2 < 3.*delta and S2**2 > delta):
                  xhyp = S2
                  xint = S3
            else:
                  xhyp = S3
                  xint = S2
            H  = 1.5*delta*xhyp**2 - 0.25*xhyp**4 + 2.*xhyp
            Sl = quartic(-0.25, 0., 1.5*delta, 2., -H) #Getting Xmin and Xmax
            Sl.sort()
            if (len(Sl) < 2):
                  print("Warning in function topology : The separatrix could not be obtained at delta =", delta)
                  return [0., 0., S1, xint, xhyp]
            if (len(Sl) == 2):
                  return [Sl[0], Sl[1], S1, xint, xhyp]
            return [Sl[2], Sl[3], S1, xint, xhyp]
            
def topology_light(delta):
      #Same as topology but only returns [Xres, Xint, Xhyp]
      if (delta == 1.):
            return [2., -1., -1.]
      Sol = cubic(1., 0., -3.*delta, -2.)
      if (len(Sol) == 1):
            return [Sol[0], 0., 0.]
      else:
            [S1, S2, S3] = Sol
            if (S2**2 < 3.*delta and S2**2 > delta):
                  xhyp = S2
                  xint = S3
            else:
                  xhyp = S3
                  xint = S2
            return [S1, xint, xhyp]

#Plotting
#if (plot_DACE_data):


def samples2ell_twoplanets(sample, pair):
      I    = pair[0]
      J    = pair[1]
      if isinstance(sample, pd.DataFrame) and isinstance(I, str) and isinstance(J, str):
            lbd1 = sample[f'mean_longitude_deg_{I}'].values
            lbd2 = sample[f'mean_longitude_deg_{J}'].values
            P1   = sample[f'period_days_{I}'].values
            P2   = sample[f'period_days_{J}'].values
            k1 = sample[f'k_{I}'].values
            k2 = sample[f'k_{J}'].values
            h1 = sample[f'h_{I}'].values
            h2 = sample[f'h_{J}'].values
            m1 = sample[f'planet_star_mass_ratio_{I}'].values
            m2 = sample[f'planet_star_mass_ratio_{J}'].values
      elif isinstance(sample, pd.DataFrame) and isinstance(I, (int, np.integer)) and isinstance(J, (int, np.integer)):
            lbd1 = sample.iloc[:, 1 + 8*I].values
            lbd2 = sample.iloc[:, 1 + 8*J].values
            P1   = sample.iloc[:, 2 + 8*I].values
            P2   = sample.iloc[:, 2 + 8*J].values
            k1   = sample.iloc[:, 3 + 8*I].values
            k2   = sample.iloc[:, 3 + 8*J].values
            h1   = sample.iloc[:, 4 + 8*I].values
            h2   = sample.iloc[:, 4 + 8*J].values
            m1   = sample.iloc[:, 7 + 8*I].values
            m2   = sample.iloc[:, 7 + 8*J].values
      else:
            lbd1 = sample[1 + 8*I,:]
            lbd2 = sample[1 + 8*J,:]
            P1   = sample[2 + 8*I,:]
            P2   = sample[2 + 8*J,:]
            k1   = sample[3 + 8*I,:]
            k2   = sample[3 + 8*J,:]
            h1   = sample[4 + 8*I,:]
            h2   = sample[4 + 8*J,:]
            m1   = sample[7 + 8*I,:]
            m2   = sample[7 + 8*J,:]
      e1   = np.sqrt(k1**2 + h1**2)
      e2   = np.sqrt(k2**2 + h2**2)
      vp1  = np.arctan2(h1, k1)
      vp2  = np.arctan2(h2, k2)
      lbd2 = lbd2*np.pi/180.
      lbd1 = lbd1*np.pi/180.

      return [e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2]




def plot_SFM(fig, ax1, Ds, x1s, x2s, pair, p, colors, label='', color_min=None, color_max=None, alpha=0.7, markersize=80):
      I    = pair[0]
      J    = pair[1]

      if (isinstance(colors, np.ndarray)):
            #Plotting
            ax1.scatter(Ds, x1s, c = colors, cmap='hsv', vmin=color_min, vmax=color_max, marker = 'o',  s = markersize, alpha = alpha, label = label + f' pair {I} {J}')
            ax1.scatter(Ds, x2s, c = colors, cmap='hsv', vmin=color_min, vmax=color_max, marker = 'o',  s = markersize, alpha = alpha)
      else:
            ax1.scatter(Ds, x1s, color = colors, marker = 'o',  s = markersize, alpha = alpha, label = label + f' pair {I} {J}')
            ax1.scatter(Ds, x2s, color = colors, marker = 'o',  s = markersize, alpha = alpha)

      if (isinstance(colors, np.ndarray)):
            cbar=fig.colorbar(mpl.cm.ScalarMappable(cmap=mpl.cm.hsv, norm=mpl.colors.Normalize(color_min, color_max)), ax=ax1, aspect=40, pad=0.01)
            cbar.ax.tick_params()


def samples2SFM(sample, pair, p):
      [e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2] = samples2ell_twoplanets(sample, pair)
      [X, Y, X2, Y2, delta] = ell2SFM(p, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2)
      return [X, Y, X2, Y2, delta]


def samples2usefull(sample, pair, p):
      [X, Y, X2, Y2, Ds] = samples2SFM(sample, pair, p)
      [sig, Sig, sig2, Sig2, x1, x2, IR] = SFM2useful(X, Y, X2, Y2, Ds)
      return [sig, Sig, sig2, Sig2, x1, x2, IR]

def plot_ell(fig, ax1, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2, pair, p, colors, label):

      [X, Y, X2, Y2, Ds] = ell2SFM(p, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2)
      [sig, Sig, sig2, Sig2, x1s, x2s, IsResonant] = SFM2useful(X, Y, X2, Y2, Ds)

      print('pair',pair, ':', 100*np.mean(IsResonant), '% within resonance.')
      plot_SFM(fig, ax1, Ds, x1s, x2s, pair, p, colors, label)


def plot_samples_SFM(fig, ax1, sample, pairs, p_indexes, colors, label):

      ### Plots the sample in the phase space of the Second Fundamental Model ###

      if isinstance(p_indexes, (list, np.ndarray)):
            if isinstance(colors, (list, np.ndarray)):
                  colors = cycle(colors)
            for pair, p, color in zip(pairs, p_indexes, colors):
                  [e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2] = samples2ell_twoplanets(sample, pair)
                  plot_ell(fig, ax1, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2, pair, p, color, label)
      else:
            [e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2] = samples2ell_twoplanets(sample, pairs)
            plot_ell(fig, ax1, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2, pairs, p_indexes, colors, label)


#def plot_samples_SFM(fig, ax1, sample, pairs, p_indexes, colors, color_min=None, color_max=None):
#
#      ### Plots the sample in the phase space of the Second Fundamental Model ###
#
#      if isinstance(p_indexes, (list, np.ndarray)):
#            #if len(pairs) != len(p_indexes):
#            #      raise ValueError('The number of planet pairs must match the number of resonances.')
#            #if isinstance(colors, (list, np.ndarray)):
#            #      colors_ = cycle(colors)
#            #else:
#            #      colors_ = [colors] * len(pairs)
#            for pair, p, color in zip(pairs, p_indexes, colors):
#                  [e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2] = samples2ell_twoplanets(sample, pair)
#                  plot_ell(fig, ax1, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2, pair, p, color, color_min, color_max)
#      else:
#            [e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2] = samples2ell_twoplanets(sample, pairs)
#            plot_ell(fig, ax1, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2, pairs, p_indexes, colors, color_min, color_max)



def plot_topology(ax1, linewidth=4, alpha=1, grid=True):

      ### Plots the topology of the phase space (separatrices and fixed points) of the Second Fundamental Model on the axis ax1 ###

      #delta_min, delta_max = delta_lim
      #X_min, X_max = X_lim

      #ax1.set_xlim(xmin = delta_min, xmax = delta_max)
      #ax1.set_ylim(ymin = X_min,     ymax = X_max)

      ax1.autoscale(axis='x')
      delta_lim = ax1.get_xlim()
      if delta_lim[0] > -3:
            delta_lim = (-3, delta_lim[1])
      if delta_lim[1] < 5:
            delta_lim = (delta_lim[0], 5)
      ax1.set_xlim(delta_lim)

      ax1.autoscale(axis='y')
      X_lim = ax1.get_ylim()
      if X_lim[0] > -5:
            X_lim = (-5, X_lim[1])
      if X_lim[1] < 5:
            X_lim = (X_lim[0], 5)
      ax1.set_ylim(X_lim)

      ax1.tick_params(axis='both', which='major')
      ax1.set_xlabel(xlabel=r"$\delta$", labelpad = 3)
      ax1.set_ylabel(ylabel=r"$X$",      labelpad = 4, rotation = 0)
      delt = np.linspace(delta_lim[0], delta_lim[1], 512)
      Xmin = np.zeros(512)
      Xmax = np.zeros(512)
      Xres = np.zeros(512)
      Xint = np.zeros(512)
      Xhyp = np.zeros(512)
      count = 0
      for delta in delt:
            [xmin, xmax, xres, xint, xhyp] = topology(delta)
            Xmin[count] = xmin
            Xmax[count] = xmax
            Xres[count] = xres
            Xint[count] = xint
            Xhyp[count] = xhyp
            count = count + 1            
      ax1.plot(delt[delt >= 1.], Xint[delt >= 1.], color = 'black', linewidth = linewidth, linestyle = '-', alpha = alpha)
      ax1.plot(delt[delt >= 1.], Xhyp[delt >= 1.], color = 'red',   linewidth = linewidth, linestyle = ':', alpha = alpha, label = 'Hyperbolic')
      ax1.plot(delt,             Xres,             color = 'black', linewidth = linewidth, linestyle = '-', alpha = alpha, label = 'Elliptic')
      ax1.plot(delt[delt >= 1.], Xmin[delt >= 1.], color = 'red',   linewidth = linewidth, linestyle = '-', alpha = alpha, label = 'Separatrix')
      ax1.plot(delt[delt >= 1.], Xmax[delt >= 1.], color = 'red',   linewidth = linewidth, linestyle = '-', alpha = alpha)
      ax1.fill_between(delt[delt >= 1.], Xmin[delt >= 1.], Xmax[delt >= 1.], color = 'red', alpha = 0.1)
      if grid:
            ax1.grid(linewidth=0.3, alpha = 0.5)



def plot_ell2SFM(data, colors=None):
      fig, ax = py.subplots(1, 1, figsize=(9,9))


      if isinstance(data, dict):
            samples = data['samples']
            analysis_id = data['samples_name']
            fig.suptitle(f'Analysis {analysis_id}', fontsize=16)
      elif isinstance(data, pd.DataFrame) or isinstance(data, np.ndarray):
            samples = data
      else:
            raise TypeError('Unsupported data type. Input has to be a pandas DataFrame, a numpy array, or a dictionary containing the "samples" key.')

      # Get the pairs of planets in resonance
      pairs_resonances = get_near_resonant_pairs(samples, 0)
      first_order_pairs = [(pair.tolist(), res, order) for pair, res, order in pairs_resonances if order == 1]
      print('Found', len(first_order_pairs), 'first order pairs.')

      # Print the pairs and their resonances
      planet_pairs, p_indexes, orders = map(list, zip(*first_order_pairs))
      print('Pairs:', planet_pairs)
      print('Resonances:', p_indexes)

      # Use a colormap if colors is not provided
      if colors is None:
            cmap = py.get_cmap('tab10')
            colors = [cmap(i % cmap.N) for i in range(len(planet_pairs))]

      # Plot the samples
      plot_samples_SFM(fig, ax, samples, planet_pairs, p_indexes, colors, label='')
      plot_topology(ax)
      py.legend()
      py.tight_layout()
      py.show()
      return fig, ax

def plot_ell2SFM_comparison(data_list, planet_pair, resonance, colors=None):
    
      fig, ax = py.subplots(figsize=(9, 9))

      cmap = py.get_cmap('tab10')
      colors = [cmap(i % cmap.N) for i in range(len(data_list))]

      for i, data in enumerate(data_list):
            samples = data['samples']
            label = f'analysis {data['analysis_id']}'
            print('Analysis ID:', data['analysis_id'])
            color = colors[i]

            if isinstance(planet_pair, np.ndarray):
                  pair = planet_pair.tolist()
            else:
                  pair = planet_pair

            # Plot for this analysis
            plot_samples_SFM(fig, ax, samples, pair, resonance, colors=color, label=label)
      plot_topology(ax)

      py.legend()
      py.tight_layout()
      py.show()
      return fig, ax






#def plot_ell2SFM(data, planet_pairs=(0,1), resonances=2, colors='green', 
#                 delta_lim=(-3,5), X_lim=(-5,5), color_lim=(None, None), check_resonance=False, 
#                 grid=True, markersize=80, alpha=0.7, linewidth=4):
#      """
#      Converts and plots the elliptic elements into the Second Fundamental Model of resonance (SFM).
#
#      Parameters
#      ----------
#      data : pandas.DataFrame or np.ndarray or dict or list of dict
#            Input data containing the posterior samples.
#            - If DataFrame or np.ndarray: used directly as sample input.
#            - If dict: must contain keys 'sample' (DataFrame) and 'samples_name' (str).
#            - If list: a list of the above dictionaries.
#      planet_pairs : list
#            Pairs of planets to be considered in the sample
#      resonance : list 
#            Resonance of the corresponding pair (p such that resonance is p:p+1).
#      colors : list 
#            List of color values to use for plotting each pair/analysis.
#            - Each entry in main list corresponds to one analysis.
#            - Each entry in nested list corresponds to one pair.
#            - Nested list entries can be strings or numpy arrays to be colormapped. 
#      delta_lim : tuple 
#            Lower and upper limits of the x-axis (delta). Set to 'auto' for automatic scaling.
#      X_lim : tuple
#            Lower and upper limits of the y-axis (X). Set to 'auto' for automatic scaling.
#      check_resonance : bool
#            If True, prints the percentage of samples within the resonance.
#      grid : bool
#            If True, adds a grid to the plot.
#      markersize : float
#            Size of the markers in the scatterplot of samples.
#      alpha : float
#            Transparency of the markers in the scatterplot of samples.
#      linewidth : float
#            Width of the lines representing the equilibrium points and separatrix. 
#      
#      Returns
#      -------
#      fig, ax : matplotlib.figure.Figure, matplotlib.axes.Axes
#            The figure and axes objects of the plot.
#      """
#      color_min, color_max = color_lim
#
#      fig, ax = py.subplots(1, 1, figsize=(9,9))
#
#      if isinstance(data, list):
#            for df_dict, color in zip(data, cycle(colors)):
#                  if check_resonance:
#                        print('Analysis', df_dict['samples_name'], ':')
#                  plot_samples_SFM(fig, ax, df_dict['samples'], planet_pairs, resonances, colors=color, color_min=color_min, color_max=color_max,
#                               label_name=df_dict['samples_name'], check_resonance=check_resonance, markersize=markersize, alpha=alpha)
#
#      elif isinstance(data, dict):
#            plot_samples_SFM(fig, ax, data['samples'], planet_pairs, resonances, colors=colors, color_min=color_min, color_max=color_max, label_name=data['samples_name'],
#                          check_resonance=check_resonance, markersize=markersize, alpha=alpha)
#
#      elif isinstance(data, pd.DataFrame) or isinstance(data, np.ndarray):
#            plot_samples_SFM(fig, ax, data, planet_pairs, resonances, colors=colors, color_min=color_min, color_max=color_max,
#                         check_resonance=check_resonance, markersize=markersize, alpha=alpha)
#
#      else:
#            raise TypeError('Unsupported data type. Input has to be a pandas DataFrame, a dictionary, or a list of dictionaries.')
#
#      if delta_lim == 'auto':
#            ax.autoscale(axis='x')
#            delta_lim = ax.get_xlim()
#            if delta_lim[0] > -3:
#                  delta_lim = (-3, delta_lim[1])
#            if delta_lim[1] < 5:
#                  delta_lim = (delta_lim[0], 5)
#            ax.set_xlim(delta_lim)
#      if X_lim == 'auto':
#            ax.autoscale(axis='y')
#            X_lim = ax.get_ylim()
#            if X_lim[0] > -5:
#                  X_lim = (-5, X_lim[1])
#            if X_lim[1] < 5:
#                  X_lim = (X_lim[0], 5)
#            ax.set_ylim(X_lim)
#      plot_topology(ax, delta_lim, X_lim)
#
#      py.legend()
#      py.tight_layout()
#      py.show()
#      return fig, ax
