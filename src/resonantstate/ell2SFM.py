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

import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib as mpl
import numpy as np

import pandas as pd


import os 


dir_path = os.path.dirname(os.path.realpath(__file__))



#plot_DACE_data = 1 # Determines if the data from a table of the GRSW are plotted in the SFM

# Defining sample to be used (only if plot_DACE_data is 1)
#path2sample = './Kepler-54_2_samples.csv'
#pairs       = [[0,1]] # Pairs of planets to be considered in the sample
#ps          = [2]     # Resonance of the corresponding pair (p such that resonance is p:p+1)
#sample      = np.loadtxt(path2sample, dtype = np.float64, delimiter=',', unpack=True)
#colors      = [sample[7,:] + sample[15,:]] # Either a list of numpy arrays or a list of strings like ['green','red']. Here the color is the total mass of the pair over the stellar mass
#if (isinstance(colors[0], np.ndarray)):
#      if (len(colors) >= 2):
#            full_color_array = np.concatenate((colors[0], colors[1]))
#            for cll in range (2, len(colors)):
#                  full_color_array = np.concatenate((full_color_array, colors[cll]))
#            color_min = min(full_color_array)
#            color_max = max(full_color_array)
#      else:
#            color_min = min(colors[0])
#            color_max = max(colors[0])

#delta_min  = -3. #To be chosen by trial and error. Irrelevant if plot_DACE_data is 0
#delta_max  = 5.
#X_min      = -5.
#X_max      = 5.

#Hamiltonian coefficients
f1s = [ 1.190493697849547, 2.025222689938346, 2.840431856715441, 3.649618244089652, 4.456142785027623, 5.261253831849899, 6.065523627097718, 6.869251916417852, 7.672611034475267, 8.475707148201764, 9.278609253466129, 10.08136413618922, 10.88400465063751, 11.68655455857515, 12.48903144896030, 13.29144865274429, 14.09381638467312, 14.89614276587963, 15.69843412935734, 16.50069558620453]
f2s = [-0.428389834143869,-2.484005183303907,-3.283256721821090,-4.083705371769611,-4.884706297511002,-5.686007411626633,-6.487489727907814,-7.289089771453291,-8.090770598035306,-8.892509248107672,-9.694290707819164,-10.49610474146903,-11.29794413782656,-12.09980367496610,-12.90167944505811,-13.70356853306293,-14.50546862185001,-15.30737796425819,-16.10929512977600,-16.91121894121170]

txt_file = os.path.join(dir_path,'continuedSeparatrix.txt')
delt, Xmin, Xmax, Xint, Xext, Xhyp = np.loadtxt(txt_file, dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 1, 2, 3, 4, 5]))

#delt, Xmin, Xmax, Xint, Xext, Xhyp = np.loadtxt('./continuedSeparatrix.txt', dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 1, 2, 3, 4, 5]))

delt = np.flip(delt)
Xmin = np.flip(Xmin)
Xmax = np.flip(Xmax)
Xint = np.flip(Xint)
Xext = np.flip(Xext)
Xhyp = np.flip(Xhyp)


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
      G     = 4.*m.pi**2
      beta1 = m1/(1. + m1)
      beta2 = m2/(1. + m2)
      mu1   = G*(1. + m1)
      mu2   = G*(1. + m2) 
      n1    = 2.*m.pi/T1
      n2    = 2.*m.pi/T2
      n10   = 2.*m.pi
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

def quartic(a_4, a_3, a_2, a_1, a_0):
      # Returns the real solutions of a_4*X^4 + a_3*X^3 + a_2*X^2 + a_1*X + a_0 = 0 from analytical expressions of Ferrari's solution
      Sol = []
      if (a_4 == 0.):
            print("Warning: The term of 4th degree must be non-zero in function quartic.")
            return []
      if (a_4 != 1.):
            return quartic(1., a_3/a_4, a_2/a_4, a_1/a_4, a_0/a_4)
      ### Equation is Y^4 + p*Y^2 + q*Y + r = 0 with Y = X + s ###
      p = a_2 - 3./8.*a_3**2
      q = a_1 + a_3**3/8. - 0.5*a_2*a_3
      r = a_0 + a_2*a_3**2/16. - a_1*a_3/4. - 3.*a_3**4/256.
      s = a_3/4.
      ### Getting M solution of the resolving cubic ###
      P  = -(r + p**2/12.)
      Q  = p*r/3. - q**2/8. - p**3/108.
      D  = Q**2+4./27.*P**3
      if (D < 0.):
            u3 = 0.5*(-Q + cm.sqrt(D))
            v3 = 0.5*(-Q - cm.sqrt(D))
            [mod_u3, arg_u3] = cm.polar(u3)
            [mod_v3, arg_v3] = cm.polar(v3)
            u  = mod_u3**(1./3.)*cm.exp(1j*arg_u3/3.)
            v  = mod_v3**(1./3.)*cm.exp(1j*arg_v3/3.)
            M  = (u + v).real - p/3.
      else:
            u3 = 0.5*(-Q + m.sqrt(D))
            v3 = 0.5*(-Q - m.sqrt(D))
            if (u3 < 0.):
                  u = -(-u3)**(1./3.)
            else:
                  u  = u3**(1./3.)
            if (v3 < 0.):
                  v = -(-v3)**(1./3.)
            else:
                  v  = v3**(1./3.)
            M  = u + v - p/3.
      if (abs(M) < 1.e-15):
            print("Warning: Problem with quartic. It may be bi-quartic.")
            return []
      if (M < 0.): #There are no real solutions
            return []
      ### Getting first pair of solutions ###
      D = q/(2.*m.sqrt(2.*M))-(p + M)/2.
      if (D >= 0.):
            S1 = -m.sqrt(2.*M)/2. + m.sqrt(D)
            S2 = -m.sqrt(2.*M)/2. - m.sqrt(D)
            Sol.append(S1 - s)
            Sol.append(S2 - s)
      ### Getting second pair of solutions ###
      D = -q/(2.*m.sqrt(2.*M))-(p + M)/2.
      if (D >= 0.):
            S3 = m.sqrt(2.*M)/2. + m.sqrt(D)
            S4 = m.sqrt(2.*M)/2. - m.sqrt(D)
            Sol.append(S3 - s)
            Sol.append(S4 - s)
      return Sol

def SFM2useful(X, Y, X2, Y2, delta):
      # Returns [sig, Sig, sig2, Sig2, nu, x1, x2] where X+iY = sqrt(2*Sig)*e^(i*sig) and X2+iY2 = sqrt(2*Sig2)*e^(i*sig2)
      # nu is the frequency of the orbit starting at (X, Y). x1 and x2 are such that (x1, 0) and (x2, 0) are on the same level line as (X, Y)
      
      Sig  = (X**2  + Y**2) /2.
      Sig2 = (X2**2 + Y2**2)/2.
      sig  = np.arctan2(Y,  X)
      sig2 = np.arctan2(Y2, X2)
      
      x1 = []
      x2 = []
      IR = []
      
      n = len(delta)
      for i in range(n):
            H   = 1.5*delta[i]*(X[i]**2 + Y[i]**2) - 0.25*(X[i]**2 + Y[i]**2)**2 + 2.*X[i]
            Sol = quartic(-0.25, 0., 1.5*delta[i], 2., -H)
            if   (len(Sol) == 0):
                  print("Warning: Problem with quartic.")
                  [xx1, xx2] = [1.e300, 1.e300]
            elif (len(Sol) == 2):
                  [xx1, xx2] = Sol
            else: #Four solutions. Retaining the two whose distance to the origin is closest to that of (X,Y)
                  if (delta[i] < 1.): #Four solutions should be impossible when delta < 1
                        print("Warning: Four solutions were found even though delta < 1")
                  [S1, S2, S3, S4] = Sol
                  D                = m.sqrt(X[i]**2 + Y[i]**2)
                  [D1, D2, D3, D4] = [abs(S1), abs(S2), abs(S3), abs(S4)]
                  DD               = np.sort(np.array([abs(D - D1), abs(D - D2), abs(D - D3), abs(D - D4)]))
                  if   (DD[0] == abs(D - D1)):
                        xx1 = S1
                  elif (DD[0] == abs(D - D2)):
                        xx1 = S2
                  elif (DD[0] == abs(D - D3)):
                        xx1 = S3
                  else:
                        xx1 = S4
                  if   (DD[1] == abs(D - D1)):
                        xx2 = S1
                  elif (DD[1] == abs(D - D2)):
                        xx2 = S2
                  elif (DD[1] == abs(D - D3)):
                        xx2 = S3
                  else:
                        xx2 = S4
            x1.append(xx1)
            x2.append(xx2)
            [xmin, xmax, xint, xext, xhyp] = topologie(delta[i])
            if (delta[i] > 1. and max(xx1,xx2) <= xmax and min(xx1,xx2) >= xmin):
                  IR.append(1)
            else:
                  IR.append(0)
      x1 = np.array(x1)
      x2 = np.array(x2)
      IR = np.array(IR)
      return [sig, Sig, sig2, Sig2, x1, x2, IR]

def topologie(delta):
      #Returns [Xmin, Xmax, Xint, Xext, Xhyp] as a function of delta
      #Instead of a direct calculation, extrapolates from file './continuedSeparatrix.txt'
      
      if (delta < -1000.):
            print("Warning : delta < -1000. Maybe p is ill-chosen.")
      
      if (delta < -20.):
            xint = -2./(3.*delta)
            xmin = xint - np.sqrt(6.)
            xmax = xint + np.sqrt(6.)
            return [xmin, xmax, xint, 0., 0.]
      	
      if (delta > 1000.):
            print("Warning : delta > 1000. Boolean IsResonant might be incorrect. Maybe p is ill-chosen.")
            return [0., 0., 0., 0., 0.]
            
      N       = len(delt[delt < delta])
      xminmin = Xmin[N - 1]
      xminmax = Xmin[N]
      xmaxmin = Xmax[N - 1]
      xmaxmax = Xmax[N]
      xintmin = Xint[N - 1]
      xintmax = Xint[N]
      xextmin = Xext[N - 1]
      xextmax = Xext[N]
      xhypmin = Xhyp[N - 1]
      xhypmax = Xhyp[N]
      Dmin    = delt[N - 1]
      Dmax    = delt[N]
      if (delta > Dmax or delta < Dmin):
            print("N = ", N)
            print("delta = ", delta)
            print("Dmin = ", Dmin)
            print("Dmax = ", Dmax)
            raise Exception("Problem with Dmin or Dmax\n")
      t    = (Dmax - delta)/(Dmax - Dmin)
      xmin = t*xminmin + (1. - t)*xminmax
      xmax = t*xmaxmin + (1. - t)*xmaxmax
      xint = t*xintmin + (1. - t)*xintmax
      xext = t*xextmin + (1. - t)*xextmax
      xhyp = t*xhypmin + (1. - t)*xhypmax
      return [xmin, xmax, xint, xext, xhyp]


#Plotting
#if (plot_DACE_data):

def plot_samples(fig, ax1, sample, pairs, ps, colors, label_name='', check_resonance=False):

      if (isinstance(colors[0], np.ndarray)):
            if (len(colors) >= 2):
                  full_color_array = np.concatenate((colors[0], colors[1]))
                  for cll in range (2, len(colors)):
                        full_color_array = np.concatenate((full_color_array, colors[cll]))
                  color_min = min(full_color_array)
                  color_max = max(full_color_array)
            else:
                  color_min = min(colors[0])
                  color_max = max(colors[0])

      #fig, ((ax1)) = py.subplots(1, 1, sharex=True, sharey=True, gridspec_kw={'width_ratios': [1]}, constrained_layout=False)
      #py.subplots_adjust(left=0.26, right=0.71, bottom=0.1, top=0.95)

      N_pairs = len(ps)
      for pair in range(N_pairs):
            I    = pairs[pair][0]
            J    = pairs[pair][1]
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
            lbd1 = lbd1*np.pi/180.
            lbd2 = lbd2*np.pi/180.
            e1   = np.sqrt(k1**2 + h1**2)
            e2   = np.sqrt(k2**2 + h2**2)
            vp1  = np.arctan2(h1, k1)
            vp2  = np.arctan2(h2, k2)
            p    = ps[pair]
            n10  = 2.*m.pi
            n20  = p*n10/(p + 1)
            n    = len(m1)
            [X, Y, X2, Y2, Ds] = ell2SFM(p, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2)
            [sig, Sig, sig2, Sig2, x1s, x2s, IsResonant] = SFM2useful(X, Y, X2, Y2, Ds)

            if check_resonance:
                  resonant = np.where(IsResonant == 1)
                  percent = (resonant[0].size / IsResonant.size)*100
                  print('pair', pairs[pair], ':', percent, '% within resonance.')

            if (isinstance(colors[pair], np.ndarray)):
                  #Making sure that all plots use the same colorbar
                  Ds    = np.concatenate((Ds,  np.array([1.e300, 1.e300])))
                  x1s   = np.concatenate((x1s, np.array([1.e300, 1.e300])))
                  x2s   = np.concatenate((x2s, np.array([1.e300, 1.e300])))
                  color = np.concatenate((colors[pair], np.array([color_min, color_max])))
                  #Plotting
                  ax1.scatter(Ds, x1s, c = color, cmap='hsv', marker = 'o',  s = 80, alpha = 0.7, label = label_name + ' ' + 'pair ' + str(I) + str(J))
                  ax1.scatter(Ds, x2s, c = color, cmap='hsv', marker = 'o',  s = 80, alpha = 0.7)
            else:
                  ax1.scatter(Ds, x1s, c = colors[pair], marker = 'o',  s = 80, alpha = 0.7, label = label_name + ' ' + 'pair ' + str(I) + str(J))
                  ax1.scatter(Ds, x2s, c = colors[pair], marker = 'o',  s = 80, alpha = 0.7)

      if (isinstance(colors[0], np.ndarray)):
            cbar=fig.colorbar(mpl.cm.ScalarMappable(cmap=mpl.cm.hsv, norm=mpl.colors.Normalize(color_min, color_max)), ax=ax1, aspect=40, pad=0.01)
            cbar.ax.tick_params()


def plot_auxiliary(ax1, delta_lim, X_lim):
      delta_min, delta_max = delta_lim
      X_min, X_max = X_lim

      ax1.set_xlim(xmin = delta_min, xmax = delta_max)
      ax1.set_ylim(ymin = X_min,     ymax = X_max)
      ax1.tick_params(axis='both', which='major')
      ax1.set_xlabel(xlabel=r"$\delta$", labelpad = 3)
      ax1.set_ylabel(ylabel=r"$X$",      labelpad = 4, rotation = 0)
      ax1.plot(delt[delt >= 1.], Xext[delt >= 1.], color = 'black', linewidth = 4, linestyle = '-', alpha = 1)
      ax1.plot(delt[delt >= 1.], Xhyp[delt >= 1.], color = 'red',   linewidth = 4, linestyle = ':', alpha = 1, label = 'Hyperbolic')
      ax1.plot(delt,             Xint,             color = 'black', linewidth = 4, linestyle = '-', alpha = 1, label = 'Elliptic')
      ax1.plot(delt[delt >= 1.], Xmin[delt >= 1.], color = 'red',   linewidth = 4, linestyle = '-', alpha = 1, label = 'Separatrix')
      ax1.plot(delt[delt >= 1.], Xmax[delt >= 1.], color = 'red',   linewidth = 4, linestyle = '-', alpha = 1)
      ax1.fill_between(delt[delt >= 1.], Xmin[delt >= 1.], Xmax[delt >= 1.], color = 'red', alpha = 0.1)
      ax1.grid(linewidth=0.3, alpha = 0.5)


def plot_ell2SFM(data, planet_pairs=[[0,1]], resonances=[2], colors=[['green']], delta_lim=(-3,5), X_lim=(-5,5), check_resonance=False):
      """
      Converts and plots the elliptic elements into the Second Fundamental Model of resonance (SFM).

      Parameters
      ----------
      data : pandas.DataFrame or dict or list of dict
            Input data containing the posterior samples.
            - If DataFrame: used directly as sample input.
            - If dict: must contain keys 'sample' (DataFrame) and 'sample_name' (str).
            - If list: a list of the above dictionaries.
      planet_pairs : list
            Pairs of planets to be considered in the sample
      resonance : list 
            Resonance of the corresponding pair (p such that resonance is p:p+1).
      colors : list 
            List of color values to use for plotting each pair/analysis.
            - Each entry in main list corresponds to one analysis.
            - Each entry in nested list corresponds to one pair.
      delta_lim : tuple 
            Lower and upper limits of the x-axis (delta).
      X_lim : tuple
            Lower and upper limits of the y-axis (X).
      """

      plot_params = planet_pairs, resonances
      ax_limits = delta_lim, X_lim

      fig, ax = py.subplots(1, 1, figsize=(7,7))
      plot_auxiliary(ax, *ax_limits)

      if isinstance(data, list):
            for df_dict, color in zip(data, colors):
                  if check_resonance:
                        print('Analysis', df_dict['sample_name'], ':')
                  sample = np.vstack([df_dict['sample'][col] for col in df_dict['sample'].columns])
                  plot_samples(fig, ax, sample, *plot_params, colors=color, 
                               label_name=df_dict['sample_name'], check_resonance=check_resonance)

      elif isinstance(data, dict):
            sample = np.vstack([data['sample'][col] for col in data['sample'].columns])
            plot_samples(fig, ax, sample, *plot_params, colors=colors[0], label_name=data['sample_name'],
                          check_resonance=check_resonance)

      elif isinstance(data, pd.DataFrame):
            sample = np.vstack([data[col] for col in data.columns])
            plot_samples(fig, ax, sample, *plot_params, colors=colors[0], check_resonance=check_resonance)

      else:
            raise TypeError('Unsupported data type. Input has to be a pandas DataFrame, a dictionary, or a list of dictionaries.')

      py.legend()
      py.tight_layout()
      py.show()                
