# This notebook manipulates 1-dimensional numpy arrays of elliptic elements to convert and plot into the Second Fundamental Model of resonance (SFM).
# The definition of the SFM used here is detailed at https://jeremycouturier.com/img/SFM.pdf
# The numpy arrays can come from posterior data of MCMC analysis or from numerical simulations.
# The user can (optionally) import a sample from the DACE table of the Geneva Resonant State Workshop (GRSW) for analysis

# Author : Jérémy COUTURIER

# Functions : 
# - ell2SFM(p, e1, e2, vp1, vp2, m1, m2, T1, T2, lbd1, lbd2) -> Returns the coordinates (X, Y, X2, Y2, delta) of the SFM. (X, Y) corresponds to the unique degree of freedom
#                                                               of the SFM and (X2, Y2) to the first integral. delta is the unique parameter of the SFM
# - SFM2useful(X, Y, X2, Y2, delta) -> Returns [sig, Sig, sig2, Sig2, nu, x1, x2] where X+iY = sqrt(2*Sig)*e^(i*sig) and X2+iY2 = sqrt(2*Sig2)*e^(i*sig2)
#                                      nu is the frequency of the orbit starting at (X, Y). x1 and x2 are such that (x1, 0) and (x2, 0) are on the same level line as (X, Y)

import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib as mpl
import numpy as np

plot_DACE_data = 1 # Determines if the data from a table of the GRSW are plotted in the SFM

# Defining sample to be used (only if plot_DACE_data is 1)
path2sample = './Kepler-54_2_samples.csv'
pairs       = [[0,1]] # Pairs of planets to be considered in the sample
ps          = [2]     # Resonance of the corresponding pair (p such that resonance is p:p+1)
colors      = ['green']
sample      = np.loadtxt(path2sample, dtype = np.float64, delimiter=',', unpack=True)

delta_min  = -3. #To be chosen by trial and error. Irrelevant if plot_DACE_data is 0
delta_max  = 5.
X_min      = -5.
X_max      = 5.

#Hamiltonian coefficients
f1s = [ 1.1904936978,  2.0252226899,  2.8404318567,  3.6496182441,  4.4561427851]
f2s = [-0.4283898341, -2.4840051833, -3.2832567218, -4.0837053718, -4.8847062975]

delt, Xmin, Xmax, Xint, Xext, Xhyp = np.loadtxt('./continuedSeparatrix.txt', dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 1, 2, 3, 4, 5]))
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
            
def rk4(X, Y, delta, dt):
      #Starting from (X,Y), integrates with a Runge-Kutta 4 until Y has changed sign twice
      #Returns the two values of X corresponding to Y=0 and the frequency of the orbit

      X0    = X
      Y0    = 0.
      count = 0
      go_on = 1
      first_time = 1
      
      while(go_on):
            Zx    = X0
            Zy    = Y0
            k1_x  = -3.*delta*Zy + Zy*(Zx**2 + Zy**2)
            k1_y  =  3.*delta*Zx - Zx*(Zx**2 + Zy**2) + 2.
            Zx    = X0 + 0.5*dt*k1_x
            Zy    = Y0 + 0.5*dt*k1_y
            k2_x  = -3.*delta*Zy + Zy*(Zx**2 + Zy**2)
            k2_y  =  3.*delta*Zx - Zx*(Zx**2 + Zy**2) + 2.
            Zx    = X0 + 0.5*dt*k2_x
            Zy    = Y0 + 0.5*dt*k2_y
            k3_x  = -3.*delta*Zy + Zy*(Zx**2 + Zy**2)
            k3_y  =  3.*delta*Zx - Zx*(Zx**2 + Zy**2) + 2.
            Zx    = X0 + dt*k3_x
            Zy    = Y0 + dt*k3_y
            k4_x  = -3.*delta*Zy + Zy*(Zx**2 + Zy**2)
            k4_y  =  3.*delta*Zx - Zx*(Zx**2 + Zy**2) + 2.
            oldY  = Y0
            oldX  = X0
            X0    = X0 + dt/6.*(k1_x + 2.*k2_x + 2.*k3_x + k4_x)
            Y0    = Y0 + dt/6.*(k1_y + 2.*k2_y + 2.*k3_y + k4_y)
            count = count + 1
            if (count >= 10000):
                  return [0., 0., 1000.]
            if (Y0 == 0. or oldY*Y0 < 0.):
                  t    = abs(Y0/(Y0 - oldY))
                  tOld = (count - 1)*dt
                  tNew = count*dt
                  if (first_time):
                        first_time = 0
                        X1 = t*oldX + (1. - t)*X0
                        t1 = t*tOld + (1. - t)*tNew
                  else:
                        go_on = 0
                        X2 = t*oldX + (1. - t)*X0
                        t2 = t*tOld + (1. - t)*tNew
      Period    = 2.*(t2 - t1)
      frequency = 2.*np.pi/Period
      return [X1, X2, frequency]

def SFM2useful(X, Y, X2, Y2, delta):
      # Returns [sig, Sig, sig2, Sig2, nu, x1, x2] where X+iY = sqrt(2*Sig)*e^(i*sig) and X2+iY2 = sqrt(2*Sig2)*e^(i*sig2)
      # nu is the frequency of the orbit starting at (X, Y). x1 and x2 are such that (x1, 0) and (x2, 0) are on the same level line as (X, Y)
      
      Sig  = (X**2  + Y**2) /2.
      Sig2 = (X2**2 + Y2**2)/2.
      sig  = np.arctan2(Y,  X)
      sig2 = np.arctan2(Y2, X2)
      
      nu_est = 3.*delta - 2.*Sig + 2.*np.cos(sig)/np.sqrt(2.*Sig)
      T_est  = 2.*np.pi/nu_est
      dt     = T_est/192.
      
      nu = []
      x1 = []
      x2 = []
      
      n = len(delta)
      for i in range(n):
            [xx1, xx2, frequency] = rk4(X[i], Y[i], delta[i], dt[i])
            nu.append(frequency)
            x1.append(xx1)
            x2.append(xx2)
      nu = np.array(nu)
      x1 = np.array(x1)
      x2 = np.array(x2)
      return [sig, Sig, sig2, Sig2, nu, x1, x2]

def topologie(delta):
      #Returns [Xmin, Xmax, Xint, Xext, Xhyp] as a function of delta
      #Instead of a direct calculation, extrapolates from file './continuedSeparatrix.txt'
            
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
if (plot_DACE_data):

      fig, ((ax1)) = py.subplots(1, 1, sharex=True, sharey=True, gridspec_kw={'width_ratios': [1]}, constrained_layout=False)
      py.subplots_adjust(left=0.26, right=0.71, bottom=0.1, top=0.95)

      row     = sample[:,0]
      n       = len(row)
      k       = n%8
      N       = (n - k)//8 #Number of planets
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
            [sig, Sig, sig2, Sig2, nus, x1s, x2s] = SFM2useful(X, Y, X2, Y2, Ds)
            ax1.scatter(Ds, x1s, c = colors[pair], marker = 'o',  s = 80, alpha = 0.4, label = 'pair ' + str(I) + str(J))
            ax1.scatter(Ds, x2s, c = colors[pair], marker = 'o',  s = 80, alpha = 0.4)


      ax1.set_xlim(xmin = delta_min, xmax = delta_max)
      ax1.set_ylim(ymin = X_min,     ymax = X_max)
      ax1.tick_params(axis='both', which='major', labelsize=25)
      ax1.set_xlabel(xlabel="$\delta$", fontsize=25, labelpad=3)
      ax1.set_ylabel(ylabel="$X$",      fontsize=25, labelpad=4, rotation=0)
      ax1.plot(delt[delt >= 1.], Xext[delt >= 1.], color = 'black',     linewidth = 4, linestyle = '-',  alpha = 1)
      ax1.plot(delt[delt >= 1.], Xhyp[delt >= 1.], color = 'lightpink', linewidth = 4, linestyle = '--', alpha = 1, label = 'Hyperbolic')
      ax1.plot(delt,             Xint,             color = 'black',     linewidth = 4, linestyle = '-',  alpha = 1, label = 'Elliptic')
      ax1.plot(delt[delt >= 1.], Xmin[delt >= 1.], color = 'lightpink', linewidth = 4, linestyle = '-',  alpha = 1, label = 'Separatrix')
      ax1.plot(delt[delt >= 1.], Xmax[delt >= 1.], color = 'lightpink', linewidth = 4, linestyle = '-',  alpha = 1)
      ax1.grid(linewidth=0.3, alpha = 0.5)
      py.legend(fontsize = 20)
      py.show()
