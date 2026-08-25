# import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib as mpl
import numpy as np
from itertools import cycle
import pandas as pd
from resonantstate.constants import *
from resonantstate.simulations_resonance_analysis import *
from resonantstate.quartic import quartic
from resonantstate.sample_conversion import Sample2cart, Cart2sample
import celeries as cl
import celeries.perham3pla as h3
import celeries.perhamavg as hv
import celeries.perhamavg_2pla as hv2
import celeries.perham as ph
import celeries.series as se
import celeries.mpfrac as mf
import celeries.laplace as lp
import celeries.ellipseries as es
import celeries.prime as pm

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
      The SFM is described by the Hamiltonian :math:`\mathcal{H}(\Sigma,\sigma)=3\delta\Sigma-\Sigma^2+2\sqrt{2\Sigma}\cos\sigma` where :math:`X+iY=\sqrt{2\Sigma}e^{i\sigma}`.
      
      Author: Jeremy Couturier. https://jeremycouturier.com
      
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

      #Period of inner planet is normalized to 1
      T2    = T2/T1
      T1    = T1/T1

      #Getting semi-major axes and Lambda
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
      
      #Converting from osculating to mean variables
      #if Nharm:
      #      X1 = np.sqrt(2.*D1/Lbd1)*np.exp(1j*vp1)
      #      X2 = np.sqrt(2.*D2/Lbd2)*np.exp(1j*vp2)
      #      if (isinstance(e1, np.ndarray)):
      #            for i in range(len(e1)):
      #                  Hv = hv2.PerHamavg(degree=1,ang2pla=((1,2,p,-p-1),),n0=(n1[i],n2[i]),masses=(1.,m1[i],m2[i]),Lbd=(Lbd1[i],Lbd2[i]),lbd=(lbd1[i],lbd2[i]),X=(X1[i],X2[i]),N=Nharm)
      #                  dLbd1 = np.real(Hv.mean('Lbd_i').toConst()) # Lbd1' - Lbd1
      #                  dLbd2 = np.real(Hv.mean('Lbd_j').toConst()) # Lbd2' - Lbd2
      #                  dlbd1 = np.real(Hv.mean('lbd_i').toConst()) # lbd1' - lbd1
      #                  dlbd2 = np.real(Hv.mean('lbd_j').toConst()) # lbd2' - lbd2
      #                  dX1 = Hv.mean('X_i').toConst() # X1' - X1
      #                  dX2 = Hv.mean('X_j').toConst() # X2' - X2
      #                  Lbd1[i] = Lbd1[i] + dLbd1
      #                  Lbd2[i] = Lbd2[i] + dLbd2
      #                  lbd1[i] = lbd1[i] + dlbd1
      #                  lbd2[i] = lbd2[i] + dlbd2
      #                  X1[i] = X1[i] + dX1
      #                  X2[i] = X2[i] + dX2
      #      else:
      #            Hv = hv2.PerHamavg(degree=1,ang2pla=((1,2,p,-p-1),),n0=(n1,n2),masses=(1.,m1,m2),Lbd=(Lbd1,Lbd2),lbd=(lbd1,lbd2),X=(X1,X2),N=Nharm)
      #            dLbd1 = np.real(Hv.mean('Lbd_i').toConst()) # Lbd1' - Lbd1
      #            dLbd2 = np.real(Hv.mean('Lbd_j').toConst()) # Lbd2' - Lbd2
      #            dlbd1 = np.real(Hv.mean('lbd_i').toConst()) # lbd1' - lbd1
      #            dlbd2 = np.real(Hv.mean('lbd_j').toConst()) # lbd2' - lbd2
      #            dX1 = Hv.mean('X_i').toConst() # X1' - X1
      #            dX2 = Hv.mean('X_j').toConst() # X2' - X2
      #            Lbd1 = Lbd1 + dLbd1
      #            Lbd2 = Lbd2 + dLbd2
      #            lbd1 = lbd1 + dlbd1
      #            lbd2 = lbd2 + dlbd2
      #            X1 = X1 + dX1
      #            X2 = X2 + dX2
      #      vp1 = np.atan2(np.imag(X1), np.real(X1))
      #      vp2 = np.atan2(np.imag(X2), np.real(X2))
      #      sq2DL1 = np.abs(X1)
      #      sq2DL2 = np.abs(X2)
      #      D1 = Lbd1*sq2DL1**2/2.
      #      D2 = Lbd2*sq2DL2**2/2.
            
      #Defining the exact resonance
      a10   = (mu1/n10**2)**(1./3.)
      a20   = (mu2/n20**2)**(1./3.)
      Lbd10 = beta1*np.sqrt(mu1*a10)
      Lbd20 = beta2*np.sqrt(mu2*a20)

      #Getting G and Gamma and normalizing
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
      
def toCanonical(v, Npla):
      r"""
      Converts the numpy array v from non-canonical (astrocentric speeds, astrocentric positions) to canonical (barycentric speeds, astrocentric positions) coordinates
      
      Author: Jeremy Couturier. https://jeremycouturier.com
      
      Parameters
      ----------
      v: Numpy array in the format of the Geneva Resonant State Workshop (GRSW).
            The columns are Id or timestamp, lbd_j(°), P_j (days), k_j, h_j, I_j(°), Omega_j(°), m_j/m_0, R_j/R_0, ..., m_0(Solar masses), R_0(Solar radii)
      Npla: Integer
            Number of planets in the planetary system. Must be consistent with v
            
      Returns
      -------
      v' : Numpy array in the format of the Geneva Resonant State Workshop (GRSW)
            The numpy array given in input, converted to canonical coordinates.
      """
      
      #### Converting to cartesian heliocentric non-canonical ####
      v_cart = Sample2cart(v, 'Heliocentric', v.shape[1]%8 - 3)
      
      #### Converting to cartesian heliocentric canonical ####
      m = 0.
      v0x = np.zeros(v.shape[0])
      v0y = np.zeros(v.shape[0])
      v0z = np.zeros(v.shape[0])
      vjx = [v_cart[:,4+8*j] for j in range(Npla)]
      vjy = [v_cart[:,5+8*j] for j in range(Npla)]
      vjz = [v_cart[:,6+8*j] for j in range(Npla)]
      mj = [v_cart[:,7+8*j] for j in range(Npla)]
      for i in range(Npla):
            v0x = v0x + mj[i]*vjx[i]
            v0y = v0y + mj[i]*vjy[i]
            v0z = v0z + mj[i]*vjz[i]
            m = m + mj[i]
      m = m + 1.
      v0x = -v0x/m
      v0y = -v0y/m
      v0z = -v0z/m
      for i in range(Npla):
            vjx[i] = vjx[i] + v0x
            vjy[i] = vjy[i] + v0y
            vjz[i] = vjz[i] + v0z
            v_cart[:,4+8*i] = vjx[i]
            v_cart[:,5+8*i] = vjy[i]
            v_cart[:,6+8*i] = vjz[i]
            
      #### Converting to elliptic heliocentric canonical ####
      out = Cart2sample(v_cart, 'Heliocentric', v.shape[1]%8 - 3)
      
      return out
      
      
      
def circular_median(a):  ##Written by Claude Opus 5##
      r"""
      Median of a sample of angles, valid across the 0 / 2*pi discontinuity.
      
      The ordinary median of angles reduced modulo 2*pi is meaningless when the sample straddles 0, since sorting then puts the discontinuity in the middle of the sorted array.
      The sample is therefore rotated so that its circular mean lies at 0, the ordinary median is taken in (-pi, pi], and the result is rotated back.
      The returned angle is not necessarily in [0, 2*pi), which does not matter since the mean longitudes only ever appear through their imaginary exponential.
      
      Parameters
      ----------
      a: 1-dimensional numpy array of floats
            The angles, in radians.
            
      Returns
      -------
      m : float
            The circular median of the angles, in radians.
      """
      
      mu = np.angle(np.mean(np.exp(1j*a)))
      return mu + np.median(np.mod(a - mu + np.pi, 2.*np.pi) - np.pi)


def ell2SFM_Npla(pair, p, v, Npla=2, Nharm=0, ang2pla=(), canonical=False, ev_flag='median'):
      r"""
      Converts from elliptic elements to the coordinates :math:`\left(X, Y, X_2, Y_2, \delta\right)` of the Second Fondamental Model of resonance (SFM).
      :math:`\left(X, Y\right)` corresponds to the unique degree of freedom of the SFM and :math:`\left(X_2, Y_2\right)` to the first integral. :math:`\delta` is the unique parameter of the SFM
      The SFM is described by the Hamiltonian :math:`\mathcal{H}(\Sigma,\sigma)=3\delta\Sigma-\Sigma^2+2\sqrt{2\Sigma}\cos\sigma` where :math:`X+iY=\sqrt{2\Sigma}e^{i\sigma}`.
      
      Author: Jeremy Couturier. https://jeremycouturier.com
      
      Parameters
      ----------
      pair: tuple of integers
            The pair (i,j) to be plotted in the SFM. The innermost planet of the system has index 1, the outermost has index Npla
      p: Integer
            Integer such that the two planets are in resonance :math:`p` : :math:`p+1`.
      v: Numpy array in the format of the Geneva Resonant State Workshop (GRSW).
            The columns are Id or timestamp, lbd_j(°), P_j (days), k_j, h_j, I_j(°), Omega_j(°), m_j/m_0, R_j/R_0, ..., m_0(Solar masses), R_0(Solar radii)
      Npla: Integer
            Number of planets in the planetary system. Must be consistent with v
      Nharm, ang2pla: Refer to class PerHamavg of celeries.
            These arguments allow the osculating coordinates to be converted into the averaged coordinates associated with the Lie serie expansion.
      canonical: Boolean
            If True, the non-canonical coordinates (astrocentric speeds, astrocentric positions) are converted to canonical coordinates (barycentric speeds, astrocentric positions).
      ev_flag: String among 'median', 'each' and 'scratch'. Matters only when the osculating coordinates are converted into mean coordinates, that is, when Nharm is not 0.
            If 'median', the difference mean - osculating is computed numerically for the median of the sample and applied as is to the full sample. Valid for compact samples only. By far fastest.
            If 'each', the difference mean - osculating is computed analytically once, and then evaluated for each line of v. Always valid.
            If 'scratch', the difference mean - osculating is computed numerically from scratch for each line of v. Always valid.
            For large samples, 'each' is faster than 'scratch'. For small samples, 'scratch' is faster than 'each' because evaluating during computation makes manipulated series smaller.
            'median' must not be used for samples of numerical simulations or spread out posterior analyses.            
            
      Returns
      -------
      l : List
            The list :math:`\left[X, Y, X_2, Y_2, \delta, e_r\right]` of numpy arrays of elements of the Second Fondamental Model of resonance.
      """

      if (p > 20 or p < 1):
            raise Exception('The index p of the resonance must be between 1 and 20 included')

      #Converting to canonical coordinates
      if canonical:
            v = toCanonical(v, Npla)

      #Extracting the elements
      G    = 4.*np.pi**2
      I,J  = pair
      lbd  = [v[:,1+8*j]*np.pi/180. for j in range(Npla)]
      lbd  = [np.mod(lbd[j],2.*np.pi) for j in range(Npla)]
      T    = [v[:,2+8*j] for j in range(Npla)]
      TI   = T[I-1]
      T    = [T[j]/TI for j in range(Npla)] #Normalizing periods
      k    = [v[:,3+8*j] for j in range(Npla)]
      h    = [v[:,4+8*j] for j in range(Npla)]
      m    = [v[:,7+8*j] for j in range(Npla)]
      e    = [np.sqrt(k[j]**2 + h[j]**2) for j in range(Npla)]
      vp   = [np.atan2(h[j], k[j]) for j in range(Npla)]
      beta = [m[j]/(1. + m[j]) for j in range(Npla)]
      mu   = [G*(1. + m[j]) for j in range(Npla)]
      n    = [2.*np.pi/T[j] for j in range(Npla)]
      a    = [(mu[j]/n[j]**2)**(1./3.) for j in range(Npla)]
      Lbd  = [beta[j]*np.sqrt(mu[j]*a[j]) for j in range(Npla)]
      D    = [Lbd[j]*(1. - np.sqrt(1. - e[j]**2)) for j in range(Npla)]
      X    = [np.sqrt(2.*D[j]/Lbd[j])*np.exp(1j*vp[j]) for j in range(Npla)]
      Nsamp= v.shape[0]
      
      #Converting from osculating to mean variables
      if Nharm:
            if ev_flag == 'median': ############################# Numerical computation at the sample's median once and for all #############################
                  Xmed = tuple([np.median(np.real(X[j])) + 1j*np.median(np.imag(X[j])) for j in range(Npla)])  ##Modified by Claude Opus 5##
                  mmed = tuple([np.median(m[j]) for j in range(Npla)])
                  Lbdmed = tuple([np.median(Lbd[j]) for j in range(Npla)])
                  lbdmed = tuple([circular_median(lbd[j]) for j in range(Npla)])  ##Modified by Claude Opus 5##
                  nmed = tuple([np.median(n[j]) for j in range(Npla)])
                  Hv = hv.PerHamavg(degree=1,Npla=Npla,ang2pla=ang2pla,n0=nmed,masses=(1.,)+mmed,Lbd=Lbdmed,lbd=lbdmed,X=Xmed,Nharm=Nharm)
                  dLbdI = np.real(Hv.mean(f'Lbd_{I}').toConst()) # Numerical value of Lbd_I' - Lbd_I at the sample's median
                  dLbdJ = np.real(Hv.mean(f'Lbd_{J}').toConst()) # Numerical value of Lbd_J' - Lbd_J at the sample's median
                  dlbdI = np.real(Hv.mean(f'lbd_{I}').toConst()) # Numerical value of lbd_I' - lbd_I at the sample's median
                  dlbdJ = np.real(Hv.mean(f'lbd_{J}').toConst()) # Numerical value of lbd_J' - lbd_J at the sample's median
                  dXI   = Hv.mean(f'X_{I}').toConst()   # Numerical value of X_I' - X_I at the sample's median
                  dXJ   = Hv.mean(f'X_{J}').toConst()   # Numerical value of X_J' - X_J at the sample's median
                  Lbd[I-1] = Lbd[I-1] + dLbdI
                  Lbd[J-1] = Lbd[J-1] + dLbdJ
                  lbd[I-1] = lbd[I-1] + dlbdI
                  lbd[J-1] = lbd[J-1] + dlbdJ
                  X[I-1] = X[I-1] + dXI
                  X[J-1] = X[J-1] + dXJ
            elif ev_flag == 'each': ############################# Analytical computation once and for all #############################
                  Hv = hv.PerHamavg(degree=1,Npla=Npla,ang2pla=ang2pla,Nharm=Nharm)
                  dLbdI = Hv.mean(f'Lbd_{I}') # Analytical expression of Lbd_I' - Lbd_I
                  dLbdJ = Hv.mean(f'Lbd_{J}') # Analytical expression of Lbd_J' - Lbd_J
                  dlbdI = Hv.mean(f'lbd_{I}') # Analytical expression of lbd_I' - lbd_I
                  dlbdJ = Hv.mean(f'lbd_{J}') # Analytical expression of lbd_J' - lbd_J
                  dXI   = Hv.mean(f'X_{I}')   # Analytical expression of X_I' - X_I
                  dXJ   = Hv.mean(f'X_{J}')   # Analytical expression of X_J' - X_J
                  for j in range(Nsamp):
                        Xj = tuple([X[i][j] for i in range(Npla)])
                        nj = tuple([n[i][j] for i in range(Npla)])
                        mj = tuple([m[i][j] for i in range(Npla)])
                        Lbdj = tuple([Lbd[i][j] for i in range(Npla)])
                        lbdj = tuple([lbd[i][j] for i in range(Npla)])
                        dLbdIj = np.real(Hv.eval(f'Lbd_{I}',dLbdI,n0=nj,masses=(1.,)+mj,Lbd=Lbdj,lbd=lbdj,X=Xj).toConst()) # Numerical evaluation of Lbd_I' - Lbd_I for each sample point
                        dLbdJj = np.real(Hv.eval(f'Lbd_{J}',dLbdJ,n0=nj,masses=(1.,)+mj,Lbd=Lbdj,lbd=lbdj,X=Xj).toConst()) # Numerical evaluation of Lbd_J' - Lbd_J for each sample point
                        dlbdIj = np.real(Hv.eval(f'lbd_{I}',dlbdI,n0=nj,masses=(1.,)+mj,Lbd=Lbdj,lbd=lbdj,X=Xj).toConst()) # Numerical evaluation of lbd_I' - lbd_I for each sample point
                        dlbdJj = np.real(Hv.eval(f'lbd_{J}',dlbdJ,n0=nj,masses=(1.,)+mj,Lbd=Lbdj,lbd=lbdj,X=Xj).toConst()) # Numerical evaluation of lbd_J' - lbd_J for each sample point
                        dXIj = Hv.eval(f'X_{I}',dXI,n0=nj,masses=(1.,)+mj,Lbd=Lbdj,lbd=lbdj,X=Xj).toConst() # Numerical evaluation of X_I' - X_I for each sample point
                        dXJj = Hv.eval(f'X_{J}',dXJ,n0=nj,masses=(1.,)+mj,Lbd=Lbdj,lbd=lbdj,X=Xj).toConst() # Numerical evaluation of X_J' - X_J for each sample point
                        Lbd[I-1][j] = Lbd[I-1][j] + dLbdIj
                        Lbd[J-1][j] = Lbd[J-1][j] + dLbdJj
                        lbd[I-1][j] = lbd[I-1][j] + dlbdIj
                        lbd[J-1][j] = lbd[J-1][j] + dlbdJj
                        X[I-1][j] = X[I-1][j] + dXIj
                        X[J-1][j] = X[J-1][j] + dXJj
            elif ev_flag == 'scratch': ############################# Numerical computation from scratch at each sample's point #############################
                  for j in range(Nsamp):
                        Xj = tuple([X[i][j] for i in range(Npla)])
                        nj = tuple([n[i][j] for i in range(Npla)])
                        mj = tuple([m[i][j] for i in range(Npla)])
                        Lbdj = tuple([Lbd[i][j] for i in range(Npla)])
                        lbdj = tuple([lbd[i][j] for i in range(Npla)])
                        Hv = hv.PerHamavg(degree=1,Npla=Npla,ang2pla=ang2pla,n0=nj,masses=(1.,)+mj,Lbd=Lbdj,lbd=lbdj,X=Xj,Nharm=Nharm)
                        dLbdIj = np.real(Hv.mean(f'Lbd_{I}').toConst()) # Numerical evaluation of Lbd_I' - Lbd_I for each sample point
                        dLbdJj = np.real(Hv.mean(f'Lbd_{J}').toConst()) # Numerical evaluation of Lbd_J' - Lbd_J for each sample point
                        dlbdIj = np.real(Hv.mean(f'lbd_{I}').toConst()) # Numerical evaluation of lbd_I' - lbd_I for each sample point
                        dlbdJj = np.real(Hv.mean(f'lbd_{J}').toConst()) # Numerical evaluation of lbd_J' - lbd_J for each sample point
                        dXIj = Hv.mean(f'X_{I}').toConst() # Numerical evaluation of X_I' - X_I for each sample point
                        dXJj = Hv.mean(f'X_{J}').toConst() # Numerical evaluation of X_J' - X_J for each sample point
                        Lbd[I-1][j] = Lbd[I-1][j] + dLbdIj
                        Lbd[J-1][j] = Lbd[J-1][j] + dLbdJj
                        lbd[I-1][j] = lbd[I-1][j] + dlbdIj
                        lbd[J-1][j] = lbd[J-1][j] + dlbdJj
                        X[I-1][j] = X[I-1][j] + dXIj
                        X[J-1][j] = X[J-1][j] + dXJj
            else:
                  raise Exception("'ev_flag' must be a string among 'median', 'each' and 'scratch' in function ell2SFM_Npla.")
                  
      #Extracting coordinates of the pair needed for ell2SFM
      X1 = X[I-1]
      X2 = X[J-1]
      Lbd1 = Lbd[I-1]
      Lbd2 = Lbd[J-1]
      lbd1 = lbd[I-1]
      lbd2 = lbd[J-1]
      m1 = m[I-1]
      m2 = m[J-1]
      vp1 = np.atan2(np.imag(X1), np.real(X1))
      vp2 = np.atan2(np.imag(X2), np.real(X2))
      sq2DL1 = np.abs(X1)
      sq2DL2 = np.abs(X2)
      D1 = Lbd1*sq2DL1**2/2.
      D2 = Lbd2*sq2DL2**2/2.
      beta1 = m1/(1. + m1)
      beta2 = m2/(1. + m2)
      mu1   = G*(1. + m1)
      mu2   = G*(1. + m2)
      a1 = Lbd1**2/(mu1*beta1**2)
      a2 = Lbd2**2/(mu2*beta2**2)
      n1 = np.sqrt(mu1/a1**3)
      n2 = np.sqrt(mu2/a2**3)
      if (abs(np.median(n1) - 2.*np.pi)/(2.*np.pi) > .2):
            print("median(n1) = ", np.median(n1))
            raise Exception("Problem with n1 in function ell2SFM_Npla.")
      n2 = n2/n1*2.*np.pi
      n1 = 2.*np.pi*np.ones(len(n2))
      n10 = 2.*np.pi
      n20 = p*n10/(p + 1)
            
      #Defining the exact resonance
      a10   = (mu1/n10**2)**(1./3.)
      a20   = (mu2/n20**2)**(1./3.)
      Lbd10 = beta1*np.sqrt(mu1*a10)
      Lbd20 = beta2*np.sqrt(mu2*a20)

      #Getting G and Gamma and normalizing
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
      e_r = np.zeros(len(delta))
      for i in range(len(delta)):
            [x1, x2] = X1X2(X[i], Y[i], delta[i])
            x_r = max(x1, x2)
            if delta[i] >= 0:
                  [xmin, xmax, xres, xint, xhyp]=topology(delta[i])
                  if x_r < xres: #Internal circulation
                        x_r = min(x1, x2) #Leleu's convention
            e_r[i] = x_r
      return [X, Y, X2, Y2, delta, e_r]

def X1X2(X, Y, delta):
      # Returns X1 and X2 such that (X1, 0) and (X2, 0) are on the same level line as (X, Y)
      H   = 1.5*(delta + 1.)*(X**2 + Y**2) - 0.25*(X**2 + Y**2)**2 + 2.*X
      Sol = quartic(-0.25, 0., 1.5*(delta + 1.), 2., -H)
      Sol.sort()
      if (len(Sol) == 0):
            print(f"Warning: No roots found by quartic in function X1X2 with delta={delta}, X={X}, Y={Y} and H={H}")
            return [0., 0.]
      elif(np.isnan(Sol[0])):
            return [np.nan, np.nan]
      if (len(Sol) == 2):
            return Sol
      #Four solutions. Either [Sol[0], Sol[3]] or [Sol[1], Sol[2]] should be returned
      if (delta < 0.): #Four solutions should be impossible when delta < 0
            print("Warning: Four solutions were found even though delta < 0 in function X1X2")
      [xmin, xmax, xres, xint, xhyp] = topology(delta)
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
            [xmin, xmax, xres, xint, xhyp] = topology(delta[i])
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
      Sol = quartic(0., 1., 0., -3.*(delta + 1.), -2.) #Getting Xres, Xint and Xhyp
      if (len(Sol) == 1):
            return [0., 0., Sol[0], 0., 0.]
      else:
            Sol.sort()
            [xhyp, xint, xres] = Sol
            # Xmin and Xmax can be obtained directly as a function of Hhyp, as shown in Petit et al. (2017), Appendix C.1.
            xmax = - xhyp + 2./np.sqrt(-xhyp)
            xmin = - xhyp - 2./np.sqrt(-xhyp)
            # Ugly but temporary bug fix
            #H3 = 1.5*(delta+1.)*xhyp**2-.25*xhyp**4+2.*xhyp
            #Sol = quartic(-.25, 0., 1.5*(delta+1.), 2., -H3)
            #if (len(Sol) == 0):
            #      raise SystemExit("Problem in topology. len(Sol) == 0")
            #Sol.sort()
            #xmax = Sol[-1]
            #xmin = Sol[-2]
            # Sanity check
            #H3 = 1.5*(delta+1.)*xhyp**2-.25*xhyp**4+2.*xhyp
            #if (abs(1.5*(delta+1.)*xmin**2-.25*xmin**4+2.*xmin - H3) > 1.e-9):
            #      print("abs(1.5*(delta+1.)*xmin**2-.25*xmin**4+2.*xmin - H3) =", abs(1.5*(delta+1.)*xmin**2-.25*xmin**4+2.*xmin - H3))
            #      raise SystemExit("problem in topology")
            #if (abs(1.5*(delta+1.)*xmax**2-.25*xmax**4+2.*xmax - H3) > 1.e-9):
            #      print("abs(1.5*(delta+1.)*xmax**2-.25*xmax**4+2.*xmax - H3) =", abs(1.5*(delta+1.)*xmax**2-.25*xmax**4+2.*xmax - H3))
            #      raise SystemExit("problem in topology")
            return [xmin, xmax, xres, xint, xhyp]

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
    [sig, Sig, sig2, Sig2, x1, x2, IR] = SFM2useful(X, Y, X2, Y2, delta)
    #[sig, Sig, sig2, Sig2, x1, x2, IR] = samples2usefull( samples[:,::sampling], pair,p)
    
    es=np.sqrt(2*Sig2)

    Der=np.zeros(delta.size)
    er=np.zeros(delta.size)


    for kl in range(delta.size):
        [_1,_2,fixp,o1,o2]=topology(delta[kl])
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


def plot_samples_SFM_pair( ax1, sample, pair, p, color, color_lim=None, label='',alpha = 1,markersize=20,marker='o',flag_DMMR=True, sampling=1):
      delta, DMMR, er, Der, es, IR, p= get_SFM_quantities(sample,pair,p=p,sampling=sampling)
      I    = pair[0]
      J    = pair[1]

      if flag_DMMR==True:
            px=DMMR
      else:
            px=delta

      if (isinstance(color, np.ndarray)):
            #Plotting
            ax1.scatter(px, er, c = color, cmap='hsv', vmin=color_lim[0], vmax=color_lim[1], marker = marker,  s = markersize, alpha = alpha, label = label + f' pair {I} {J}')
            
      else:
            ax1.scatter(px, er, color = color, marker = marker,  s = markersize, alpha = alpha, label = label + f' pair {I} {J}')


def plot_samples_SFM(ax1, sample, pairs, p_indexes, colors, color_lim=None, label='',alpha = 1,markersize=20,marker='o',flag_DMMR=True, sampling=1, Nharm=0, ang2pla=(), canonical=False, ev_flag='median'):
      ### Plots the sample in the phase space of the Second Fundamental Model ###
      if isinstance(sample, pd.DataFrame):
            v = sample.to_numpy()
      else:
            v = sample
      v = v[::sampling,:]
      if isinstance(p_indexes, (list, np.ndarray)):
            if isinstance(colors, (list, np.ndarray)):
                  colors = cycle(colors)
            for pair, p, color in zip(pairs, p_indexes, colors):
                  if (Nharm or canonical):
                        I,J=pair
                        pair = (pair[0]+1, pair[1]+1)
                        [X, Y, X2, Y2, delta, er] = ell2SFM_Npla(pair, p, v, Npla=(v.shape[1]-v.shape[1]%8)//8, Nharm=Nharm, ang2pla=ang2pla, canonical=canonical, ev_flag=ev_flag)
                        if flag_DMMR:
                              ax1.scatter(er**2-3*(delta+1), er, marker=marker, s=markersize, c=color, label=label + f' pair {I} {J}', alpha=alpha)
                        else:
                              ax1.scatter(delta, er, marker=marker, s=markersize, c=color, label=label + f' pair {I} {J}', alpha=alpha)
                  else:
                        plot_samples_SFM_pair( ax1, sample, pair, p, color, color_lim=color_lim, label=label,alpha =alpha,markersize=markersize,marker=marker,flag_DMMR=flag_DMMR, sampling=sampling)
                           
      else:
            if (Nharm or canonical):
                  I,J=pairs
                  pairs = (pairs[0]+1, pairs[1]+1)
                  [X, Y, X2, Y2, delta, er] = ell2SFM_Npla(pairs, p_indexes, v, Npla=(v.shape[1]-v.shape[1]%8)//8, Nharm=Nharm, ang2pla=ang2pla, canonical=canonical, ev_flag=ev_flag)
                  if flag_DMMR:
                        ax1.scatter(er**2-3*(delta+1), er, marker=marker, s=markersize, c=colors, label=label + f' pair {I} {J}', alpha=alpha)
                  else:
                        ax1.scatter(delta, er, marker=marker, s=markersize, c=colors, label=label + f' pair {I} {J}', alpha=alpha)
            else:
                  plot_samples_SFM_pair( ax1, sample, pairs, p_indexes, colors, color_lim=color_lim, label=label,alpha =alpha,markersize=markersize,marker=marker,flag_DMMR=flag_DMMR, sampling=sampling)
                  

