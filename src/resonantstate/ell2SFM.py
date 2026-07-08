# import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib as mpl
import numpy as np
from itertools import cycle
import pandas as pd
from resonantstate.constants import *
from resonantstate.simulations_resonance_analysis import *


Tp=np.arange(1,20)
def _dist_Dpp1p(Pr):


    if Pr>2.5 or Pr<1.1:
        return np.inf,-1,False
    
    else :
        dist_H= -2*(np.abs(f1s[Tp-1])/3/Tp**(1/2))**(-2/3)*(1-(Pr*Tp/(Tp+1))**(1/3))
        Id_closest=np.abs(dist_H).argmin()

        dist_H_best,p_best=dist_H[Id_closest],Tp[Id_closest]

        if dist_H_best>0:
            Pr_test=(2-(p_best/(1+p_best)*Pr)**(1/3))**3*(p_best+1)/p_best
            dist_H_test= -2*(np.abs(f1s[Tp-1])/3/Tp**(1/2))**(-2/3)*(1-(Pr_test*Tp/(Tp+1))**(1/3))
            Id_closest_test=np.abs(dist_H_test).argmin()
            
            if Tp[Id_closest_test]==p_best:
                sym=True
            else:
                sym=False
        else:
            sym=True

        
        return dist_H[Id_closest],Tp[Id_closest],sym

def dist_Dpp1p(Pr):
    """deal with different entry shape for Pr

    Parameters
    ----------
    Pr : float or array
        period ratios

    Returns
    -------
    dist : float

    p : integer

    sym: boolean
        
    """

    if isinstance(Pr,np.ndarray):
        if len(Pr.shape)==1:
            dist_H,Tp,sym=np.zeros(Pr.size),np.zeros(Pr.size),np.zeros(Pr.size,dtype=bool)
            for k,pr in enumerate(Pr):
                dist_H[k],Tp[k],sym[k]=_dist_Dpp1p(pr)
            return dist_H,Tp,sym
    else:
        return  _dist_Dpp1p(Pr)
    
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

def ell2SFM(p, e1, e2, vp1, vp2, m1, m2, T1, T2, lbd1, lbd2):    
      r"""
      Converts from elliptic elements to the coordinates :math:`\left(X, Y, X_2, Y_2, \delta\right)` of the Second Fondamental Model of resonance (SFM).
      :math:`\left(X, Y\right)` corresponds to the unique degree of freedom of the SFM and :math:`\left(X_2, Y_2\right)` to the first integral. :math:`\delta` is the unique parameter of the SFM
      The SFM is described byt the Hamiltonian :math:`\mathcal{H}(\Sigma,\sigma)=3\delta\Sigma-\Sigma^2+2\sqrt{2\Sigma}\cos\sigma` where :math:`X+iY=\sqrt{2\Sigma}e^{i\sigma}`.
      
      Author : Jeremy Couturier. https://jeremycouturier.com
      
      Parameters
      ----------
      p: int
            Integer such that the two planets are in resonance :math:`p` : :math:`p+1`.
      e1: float or 1-dimensional numpy array of floats
            The eccentricity :math:`e_1` of the inner planet.
      e2: float or 1-dimensional numpy array of floats
            The eccentricity :math:`e_2` of the outer planet.
      vp1: float or 1-dimensional numpy array of floats
            The longitude of periapsis :math:`\varpi_1` of the inner planet, in radians.
      vp2: float or 1-dimensional numpy array of floats
            The longitude of periapsis :math:`\varpi_2` of the outer planet, in radians.
      m1: float or 1-dimensional numpy array of floats
            The mass :math:`m_1` of the inner planet in units of the stellar mass.
      m2: float or 1-dimensional numpy array of floats
            The mass :math:`m_2` of the outer planet in units of the stellar mass.
      T1: float or 1-dimensional numpy array of floats
            The period :math:`2\pi/n_1` of the inner planet in any units. Such that :math:`\mathcal{G}\left(m_0+m_1\right)=n_1^2a_1^3`.
      T2: float or 1-dimensional numpy array of floats
            The period :math:`2\pi/n_2` of the outer planet in the same units as T1. Such that :math:`\mathcal{G}\left(m_0+m_2\right)=n_2^2a_2^3`.
      lbd1: float or 1-dimensional numpy array of floats
            The mean longitude :math:`\lambda_1` of the inner planet, in radians.
      lbd2: float or 1-dimensional numpy array of floats
            The mean longitude :math:`\lambda_2` of the outer planet, in radians.
            
      Returns
      -------
      l : List
            The list :math:`\left[X, Y, X_2, Y_2, \delta\right]` of elements of the Second Fondamental Model of resonance.
      """

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
      DG    = Lbd1 - Lbd10 + Lbd2 - Lbd20 - D1 - D2
      Gamma = (p + 1)*Lbd1 + p*Lbd2
      DGamma= Gamma - ((p + 1)*Lbd10 + p*Lbd20)
      g     = G/Gamma
      Dg    = DG/Gamma
      Dgamma= DGamma/Gamma
      d1    = D1/Gamma
      d2    = D2/Gamma
      C1    = Gamma/Lbd10
      C2    = Gamma/Lbd20

      #Getting alpha, beta, gamma, delta, R and S
      f1    = f1s[p - 1]
      f2    = f2s[p - 1]
      R     = (f1**2*C1*d1 + f2**2*C2*d2 + 2.*f1*f2*np.sqrt(C1*d1*C2*d2)*np.cos(vp1 - vp2))/(f1**2*C1 + f2**2*C2)
      S     = (f1**2*C1*d2 + f2**2*C2*d1 - 2.*f1*f2*np.sqrt(C1*d1*C2*d2)*np.cos(vp1 - vp2))/(f1**2*C1 + f2**2*C2)
      #alpha = -3.*n10*p*((g + S)*(p*C1 + (p + 1)*C2) - (C1 + C2)) #Old expression. Equal to the new one although written differently
      alpha = -3.*n10*p*((Dg + S)*(p*C1 + (p + 1)*C2) - (C1 + C2)*Dgamma)
      beta  = 1.5*n10*p*(p*C1 + (p + 1)*C2)
      gamma = m1*n20/C2*np.sqrt(f1**2*C1 + f2**2*C2)
      delta = alpha*(4./(27.*beta*gamma**2))**(1./3.) - 1.

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
      H   = 1.5*(delta + 1.)*(X**2 + Y**2) - 0.25*(X**2 + Y**2)**2 + 2.*X
      Sol = quartic(-0.25, 0., 1.5*(delta + 1.), 2., -H)
      Sol.sort()
      if (len(Sol) == 0):
            print("Warning: Problem with quartic.")
            return [0., 0.]
      if (len(Sol) == 2):
            return Sol
      #Four solutions. Either [Sol[0], Sol[3]] or [Sol[1], Sol[2]] should be returned
      if (delta < 0.): #Four solutions should be impossible when delta < 0
            print("Warning: Four solutions were found even though delta < 0 in function X1X2")
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
            [lib1,lib2]= X1X2(0, 0, delta[i])
            if (delta[i] < 1.):
                  
                  if  max(xx1,xx2)<max(lib1,lib2):    
                        IR.append(2) #external libration
                  else:
                        IR.append(0) #external circulation
            else:
                  Hseparatrix = 1.5*(delta[i] + 1.)*xhyp**2 - 0.25*xhyp**4 + 2.*xhyp
                  H           = 1.5*(delta[i] + 1.)*(X[i]**2 + Y[i]**2) - 0.25*(X[i]**2 + Y[i]**2)**2 + 2.*X[i]
                  if (H > Hseparatrix):
                        IR.append(1) #formally resonant
                  else:
                        if max(xx1,xx2)<xres:
                              if min(xx1,xx2)>min(lib1,lib2):
                                    IR.append(-2) #internal libration
                              else:
                                    IR.append(-1) #internal circulation
                        else:
                              IR.append(0) #external circulation

      x1 = np.array(x1)
      x2 = np.array(x2)
      IR = np.array(IR)
      return [sig, Sig, sig2, Sig2, x1, x2, IR]

def topology(delta):
      #Returns [Xmin, Xmax, Xres, Xint, Xhyp] from analytical expressions instead of reading from file
      #Xmin and Xmax are the lower and upper separatrices respectively, Xres is the resonance center, Xint is the center of the internal circulation and Xhyp is the hyperbolic fixed point.
      #If delta < 0, returns [0, 0, Xres, 0, 0]
      #When delta >= 0; then Xhyp <= Xint <= Xmin <= Xres <= Xmax
      if (delta == 0.):
            return [-1., 3., 2., -1., -1.]
      Sol = cubic(1., 0., -3.*(delta + 1.), -2.) #Getting Xres, Xint and Xhyp
      if (len(Sol) == 1):
            return [0., 0., Sol[0], 0., 0.]
      else:
            [S1, S2, S3] = Sol
            if (S2**2 < 3.*(delta + 1.) and S2**2 > (delta + 1.)):
                  xhyp = S2
                  xint = S3
            else:
                  xhyp = S3
                  xint = S2
            H  = 1.5*(delta + 1.)*xhyp**2 - 0.25*xhyp**4 + 2.*xhyp
            Sl = quartic(-0.25, 0., 1.5*(delta + 1.), 2., -H) #Getting Xmin and Xmax
            Sl.sort()
            if (len(Sl) < 2):
                  print("Warning in function topology : The separatrix could not be obtained at delta =", delta)
                  return [0., 0., S1, xint, xhyp]
            if (len(Sl) == 2):
                  return [Sl[0], Sl[1], S1, xint, xhyp]
            return [Sl[2], Sl[3], S1, xint, xhyp]
            
def topology_light(delta):
      #Same as topology but only returns [Xres, Xint, Xhyp]
      if (delta == 0.):
            return [2., -1., -1.]
      Sol = cubic(1., 0., -3.*(delta + 1.), -2.)
      if (len(Sol) == 1):
            return [Sol[0], 0., 0.]
      else:
            [S1, S2, S3] = Sol
            if (S2**2 < 3.*(delta + 1.) and S2**2 > (delta + 1.)):
                  xhyp = S2
                  xint = S3
            else:
                  xhyp = S3
                  xint = S2
            return [S1, xint, xhyp]

def samples2SFM(sample, pair, p):
      [e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2] = samples2ell_twoplanets(sample, pair)
      [X, Y, X2, Y2, delta] = ell2SFM(p, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2)
      return [X, Y, X2, Y2, delta]

def samples2usefull(sample, pair, p):
      [X, Y, X2, Y2, Ds] = samples2SFM(sample, pair, p)
      [sig, Sig, sig2, Sig2, x1, x2, IR] = SFM2useful(X, Y, X2, Y2, Ds)
      return [sig, Sig, sig2, Sig2, x1, x2, IR]

def get_SFM_quantities(samples,pair,p=None,sampling=1):
    """call functions of the package to extract the quantities used in the first paper : delta, DMMR, er, Der, es, IR, p

    Parameters
    ----------
    samples : pandas dataframe
        posterior of a dynamical analysis
    pair : tuple of ints
        identifiant of the pair to consider
    p : int, optional
        p+1:p first order MMR to consider
    sampling : int, optional
        thinning of the samples

    Returns
    -------
    delta : np.array of floats of the lengh of the thinned samples
        delta parameter of the second fundamental model for resonance (Henrard and Lemaitre 1983)
    DMMR : np.array of floats of the lengh of the thinned samples
        er^2-3*delta
    er : np.array of floats of the lengh of the thinned samples
        rer parameter of the second fundamental model for resonance (see GRSW paper 1 - Leleu et al 2026) also called resonant eccentricities
    Der : np.array of floats of the lengh of the thinned samples
        Der parameter of the second fundamental model for resonance (see GRSW paper 1 - Leleu et al 2026) 
    es : np.array of floats of the lengh of the thinned samples
        secular eccentricities (see GRSW paper 1 - Leleu et al 2026) 
    IR : np.array of floats of the lengh of the thinned samples
        position of the sample in the second fundamental model:
            2 : upper libration
            1 : formally resonant
            0 : upper circulation
            -1 : lower circulation
            -2 : lower libration
    p : int
        pair attributed to the p+1:p MMR
    
    """

    P1=np.mean(samples['period_days_'+str(pair[0])].values)
    P2=np.mean(samples['period_days_'+str(pair[1])].values)

    if isinstance(samples, pd.DataFrame):
            samples=np.vstack([samples[col] for col in samples.columns])
    

    if p==None:
        dist,p,sym=dist_Dpp1p(P2/P1)
    
    [X, Y, X2, Y2, delta]=samples2SFM( samples[:,::sampling], pair, p)
    [sig, Sig, sig2, Sig2, x1, x2, IR] = samples2usefull( samples[:,::sampling], pair,p)
    
    es=np.sqrt(2*Sig2)

    Der=np.zeros(delta.size)
    er=np.zeros(delta.size)


    for kl in range(delta.size):
        [fixp,o1,o2]=topology_light(delta[kl])
        Vx=np.array([x1[kl], x2[kl]])
        Idmaxx=abs(Vx).argmax()
        Idxlib=abs(Vx[Idmaxx]-np.array([fixp,o1,o2])).argmin()
        Der[kl]=Vx[Idmaxx]-np.array([fixp,o1,o2])[Idxlib]
        er[kl]=Vx[Idmaxx]



    DMMR=er**2-3*(delta+1)


    return delta, DMMR, er, Der, es, IR, p



#Plotting

def plot_topology(ax1, delta_lim=None, X_lim=None, linewidth=4, alpha=1, grid=False,dark=False,alpha_fill=0.1,legend=True):

      ### Plots the topology of the phase space (separatrices and fixed points) of the Second Fundamental Model on the axis ax1 ###

      #delta_min, delta_max = delta_lim
      #X_min, X_max = X_lim

      #ax1.set_xlim(xmin = delta_min, xmax = delta_max)
      #ax1.set_ylim(ymin = X_min,     ymax = X_max)

      if delta_lim==None:
            ax1.autoscale(axis='x')
            delta_lim = ax1.get_xlim()
            if delta_lim[0] > -3:
                  delta_lim = (-3, delta_lim[1])
            if delta_lim[1] < 5:
                  delta_lim = (delta_lim[0], 5)


      ax1.set_xlim(delta_lim)

      if X_lim==None:
            ax1.autoscale(axis='y')
            X_lim = ax1.get_ylim()
            if X_lim[0] > -5:
                  X_lim = (-5, X_lim[1])
            if X_lim[1] < 5:
                  X_lim = (X_lim[0], 5)
      ax1.set_ylim(X_lim)

      ax1.tick_params(axis='both', which='major')
      ax1.set_xlabel(xlabel=r"$\delta$", labelpad = 3)
      ax1.set_ylabel(ylabel=r"$x_r$",    labelpad = 4, rotation = 0)
      delt = np.linspace(delta_lim[0], delta_lim[1], 1024)
      Xmin = np.zeros(1024)
      Xmax = np.zeros(1024)
      Xres = np.zeros(1024)
      Xint = np.zeros(1024)
      Xhyp = np.zeros(1024)
      count = 0
      for delta in delt:
            [xmin, xmax, xres, xint, xhyp] = topology(delta)
            Xmin[count] = xmin
            Xmax[count] = xmax
            Xres[count] = xres
            Xint[count] = xint
            Xhyp[count] = xhyp
            count = count + 1            
      if dark:        
            ax1.plot(delt[delt >= 0.], Xint[delt >= 0.], color = 'white', linewidth = linewidth, linestyle = '-', alpha = alpha)
            if legend:
                  ax1.plot(delt[delt >= 0.], Xhyp[delt >= 0.], color = 'pink',  linewidth = linewidth, linestyle = ':', alpha = alpha, label = 'Hyperbolic')
                  ax1.plot(delt,             Xres,             color = 'white', linewidth = linewidth, linestyle = '-', alpha = alpha, label = 'Elliptic')
                  ax1.plot(delt[delt >= 0.], Xmin[delt >= 0.], color = 'pink',  linewidth = linewidth, linestyle = '-', alpha = alpha, label = 'Separatrix')
            else:                        
                  ax1.plot(delt[delt >= 0.], Xhyp[delt >= 0.], color = 'pink',  linewidth = linewidth, linestyle = ':', alpha = alpha)
                  ax1.plot(delt,             Xres,             color = 'white', linewidth = linewidth, linestyle = '-', alpha = alpha)
                  ax1.plot(delt[delt >= 0.], Xmin[delt >= 0.], color = 'pink',  linewidth = linewidth, linestyle = '-', alpha = alpha)
            ax1.plot(      delt[delt >= 0.], Xmax[delt >= 0.], color = 'pink',  linewidth = linewidth, linestyle = '-', alpha = alpha)
            ax1.fill_between(delt[delt >= 0.], Xmin[delt >= 0.], Xmax[delt >= 0.], color = 'pink', alpha = alpha_fill)
      else:
            ax1.plot(      delt[delt >= 0.], Xint[delt >= 0.], color = 'black', linewidth = linewidth, linestyle = '-', alpha = alpha)
            if legend:
                  ax1.plot(delt[delt >= 0.], Xhyp[delt >= 0.], color = 'red',   linewidth = linewidth, linestyle = ':', alpha = alpha, label = 'Hyperbolic')
                  ax1.plot(delt,             Xres,             color = 'black', linewidth = linewidth, linestyle = '-', alpha = alpha, label = 'Elliptic')
                  ax1.plot(delt[delt >= 0.], Xmin[delt >= 0.], color = 'red',   linewidth = linewidth, linestyle = '-', alpha = alpha, label = 'Separatrix')
            else:
                  ax1.plot(delt[delt >= 0.], Xhyp[delt >= 0.], color = 'red',   linewidth = linewidth, linestyle = ':', alpha = alpha)
                  ax1.plot(delt,             Xres,             color = 'black', linewidth = linewidth, linestyle = '-', alpha = alpha)
                  ax1.plot(delt[delt >= 0.], Xmin[delt >= 0.], color = 'red',   linewidth = linewidth, linestyle = '-', alpha = alpha)    
            ax1.plot(delt[delt >= 0.], Xmax[delt >= 0.], color = 'red',   linewidth = linewidth, linestyle = '-', alpha = alpha)
            ax1.fill_between(delt[delt >= 0.], Xmin[delt >= 0.], Xmax[delt >= 0.], color = 'red', alpha = alpha_fill)

      if grid:
            ax1.grid(linewidth=0.3, alpha = 0.5)


def plot_topology_DMMR(ax1, DMMR_lim=None, X_lim=None, linewidth=4, alpha=1, grid=False,dark=False,legend=True):

      ### Plots the topology of the phase space (separatrices and fixed points) of the Second Fundamental Model on the axis ax1 ###

      if DMMR_lim==None:
            ax1.autoscale(axis='x')
            DMMR_lim = ax1.get_xlim()
            if DMMR_lim[0] > -3:
                  DMMR_lim = (-3, DMMR_lim[1])
            if DMMR_lim[1] < 5:
                  DMMR_lim = (DMMR_lim[0], 5)


      ax1.set_xlim(DMMR_lim)

      if X_lim==None:
            ax1.autoscale(axis='y')
            X_lim = ax1.get_ylim()
            if X_lim[0] > -5:
                  X_lim = (-5, X_lim[1])
            if X_lim[1] < 5:
                  X_lim = (X_lim[0], 5)

      
      bound=(max(np.abs(X_lim))**2+max(np.abs(DMMR_lim)))/3
      ax1.set_ylim(X_lim)

      ax1.tick_params(axis='both', which='major')
      ax1.set_xlabel(xlabel=r"$\Delta_{MMR}$", labelpad = 3)
      ax1.set_ylabel(ylabel=r"$e_r$",          labelpad = 4, rotation = 0)
      #delt = np.linspace(delta_lim[0]-X2max, delta_lim[1]+X2max, 512)
      delt = np.arange(-bound, bound, .01)
      Xmin = np.zeros(delt.size)
      Xmax = np.zeros(delt.size)
      Xres = np.zeros(delt.size)
      Xint = np.zeros(delt.size)
      Xhyp = np.zeros(delt.size)
      count = 0
      for count,delta in enumerate(delt):
            [xmin, xmax, xres, xint, xhyp] = topology(delta)
            Xmin[count] = xmin
            Xmax[count] = xmax
            Xres[count] = xres
            Xint[count] = xint
            Xhyp[count] = xhyp

      if dark:        
            
            if legend:
                  ax1.plot(Xint[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xint[delt >= 0.], color = 'grey', linewidth = linewidth, linestyle = '-', alpha = alpha, label = 'secondary')
                  ax1.plot(Xhyp[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xhyp[delt >= 0.], color = 'pink',   linewidth = linewidth, linestyle = ':', alpha = alpha, label = 'hyperbolic')
                  ax1.plot(Xres**2-3*(delt + 1.),             Xres,             color = 'white', linewidth = linewidth, linestyle = '-', alpha = alpha, label = 'primary')
                  ax1.plot(Xmin[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xmin[delt >= 0.], color = 'pink',   linewidth = linewidth, linestyle = '-', alpha = alpha, label = 'separatrix')
            else:
                  ax1.plot(Xint[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xint[delt >= 0.], color = 'grey', linewidth = linewidth, linestyle = '-', alpha = alpha)
                  ax1.plot(Xhyp[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xhyp[delt >= 0.], color = 'pink',   linewidth = linewidth, linestyle = ':', alpha = alpha)
                  ax1.plot(Xres**2-3*(delt + 1.),             Xres,             color = 'white', linewidth = linewidth, linestyle = '-', alpha = alpha)
                  ax1.plot(Xmin[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xmin[delt >= 0.], color = 'pink',   linewidth = linewidth, linestyle = '-', alpha = alpha)
            ax1.plot(Xmax[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xmax[delt >= 0.], color = 'pink',   linewidth = linewidth, linestyle = '-', alpha = alpha)
            #ax1.fill_between(delt[delt >= 0.], Xmin[delt >= 0.], Xmax[delt >= 0.], color = 'pink', alpha = alpha_fill)
      else:
            
            if legend:
                  ax1.plot(Xint[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xint[delt >= 0.], color = 'grey', linewidth = linewidth, linestyle = '-', alpha = alpha, label = 'secondary')
                  ax1.plot(Xhyp[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xhyp[delt >= 0.], color = 'red',   linewidth = linewidth, linestyle = ':', alpha = alpha, label = 'hyperbolic')
                  ax1.plot(Xres**2-3*(delt + 1.),             Xres,             color = 'black', linewidth = linewidth, linestyle = '-', alpha = alpha, label = 'primary')
                  ax1.plot(Xmin[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xmin[delt >= 0.], color = 'red',   linewidth = linewidth, linestyle = '-', alpha = alpha, label = 'separatrix')
            else:
                  ax1.plot(Xint[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xint[delt >= 0.], color = 'grey', linewidth = linewidth, linestyle = '-', alpha = alpha)
                  ax1.plot(Xhyp[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xhyp[delt >= 0.], color = 'red',   linewidth = linewidth, linestyle = ':', alpha = alpha)
                  ax1.plot(Xres**2-3*(delt + 1.),             Xres,             color = 'black', linewidth = linewidth, linestyle = '-', alpha = alpha)
                  ax1.plot(Xmin[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xmin[delt >= 0.], color = 'red',   linewidth = linewidth, linestyle = '-', alpha = alpha)    
            ax1.plot(Xmax[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xmax[delt >= 0.], color = 'red',   linewidth = linewidth, linestyle = '-', alpha = alpha)
            #ax1.fill_between(delt[delt >= 0.], Xmin[delt >= 0.], Xmax[delt >= 0.], color = 'red', alpha = alpha_fill)

      if grid:
            ax1.grid(linewidth=0.3, alpha = 0.5)


def plot_SFM(fig, ax1, Ds, x1s, x2s, pair, p, colors, color_lim=None, label='', alpha=0.7, markersize=80, marker= 'o'):
      I    = pair[0]
      J    = pair[1]

      if (isinstance(colors, np.ndarray)):
            #Plotting
            ax1.scatter(Ds, x1s, c = colors, cmap='hsv', vmin=color_lim[0], vmax=color_lim[1], marker = marker,  s = markersize, alpha = alpha, label = label + f' pair {I} {J}')
            ax1.scatter(Ds, x2s, c = colors, cmap='hsv', vmin=color_lim[0], vmax=color_lim[1], marker = marker,  s = markersize, alpha = alpha)
      else:
            ax1.scatter(Ds, x1s, color = colors, marker = marker,  s = markersize, alpha = alpha, label = label + f' pair {I} {J}')
            ax1.scatter(Ds, x2s, color = colors, marker = marker,  s = markersize, alpha = alpha)


def plot_ell(fig, ax1, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2, pair, p, colors, color_lim, label,alpha = 1,markersize=20, marker = 'o'):

      [X, Y, X2, Y2, Ds] = ell2SFM(p, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2)
      [sig, Sig, sig2, Sig2, x1s, x2s, IsResonant] = SFM2useful(X, Y, X2, Y2, Ds)

      print('pair',pair, ':', 100*np.mean(IsResonant), '% within resonance.')
      plot_SFM(fig, ax1, Ds, x1s, x2s, pair, p, colors, color_lim, label,alpha =alpha,markersize=markersize, marker = marker)


def plot_samples_SFM(fig, ax1, sample, pairs, p_indexes, colors, color_lim=None, label='',alpha = 1,markersize=20,marker='o'):

      ### Plots the sample in the phase space of the Second Fundamental Model ###

      if isinstance(p_indexes, (list, np.ndarray)):
            if isinstance(colors, (list, np.ndarray)):
                  colors = cycle(colors)
            for pair, p, color in zip(pairs, p_indexes, colors):
                  [e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2] = samples2ell_twoplanets(sample, pair)
                  plot_ell(fig, ax1, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2, pair, p, color, color_lim, label,alpha =alpha, markersize=markersize, marker=marker)
      else:
            [e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2] = samples2ell_twoplanets(sample, pairs)
            plot_ell(fig, ax1, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2, pairs, p_indexes, colors, color_lim, label,alpha = alpha, markersize=markersize, marker=marker)



def get_p_by_pair(samples):

      #Returns (pairs, ps) where pairs = [(i1, j1), (i2, j2), ...] is a list of pairs close to a resonance p:p+1 and ps = [p1, p2, ...] contains the associated values of p
      #samples is a dataframe
      
      data = samples[samples.columns[samples.columns.str.contains('period')]]
      
      stdd = data['period_days_0'].std() #If the standard deviation of the periods is very low, the sample is likely an observation, otherwise, it is likely a simulation

      if (stdd < 0.05): #The sample is probably an observation, we take the median
            data = data.median().sort_values()
      else:            #The sample is probably a simulation, we take the last row
            last_row = samples.index[-1]
            data = data.iloc[last_row].sort_values()
            
      periods = np.array(data.values)
      N       = len(periods)
      
      pairs = []
      ps    = []
      for j in range(N):
            for i in range(j):
                  #Pair is (i,j)                   
                  PeriodRatio = periods[j]/periods[i]
                  if (PeriodRatio < 2.1 and PeriodRatio > 1.0488):
                        distance = 1.e300
                        for p in range(1, 21): #Resonance 21:22 and onward are unsupported
                              if (abs(PeriodRatio - (p + 1)/p)/((p + 1)/p) < distance):
                                    distance = abs(PeriodRatio - (p + 1)/p)/((p + 1)/p)
                                    the_p    = p
                        if (distance < 0.05): 
                              pairs.append((i,j))
                              ps.append(the_p)
      return (pairs, ps)





