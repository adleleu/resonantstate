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


Tp=np.arange(1,20)



def _dist_Dpp1p(Pr):
    r"""
      Compute the distance to the closest first-order resonance, and verify if it lands in the symetrical domain around that resonance

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
    
    Tp=np.arange(1,21)

    if Pr>2.5 or Pr<1.1:
        return np.inf,-1,False
    
    else :
        pspp1=Tp/(Tp+1)
        Kp=-2*(np.abs(f1s[Tp-1])/3/Tp**(1/2))**(-2/3)
        dist_H= (1-(Pr*Tp/(Tp+1))**(1/3))*Kp*(np.sqrt(pspp1**(2/3)+pspp1**(5/3)*(f2s/f1s)**2+pspp1**(4/3)+pspp1*(f2s/f1s)**2)/2)**(-2/3)
        Id_closest=np.abs(dist_H).argmin()

        dist_H_best,p_best=dist_H[Id_closest],Tp[Id_closest]

        if dist_H_best>0:
            Pr_test=(2-(p_best/(1+p_best)*Pr)**(1/3))**3*(p_best+1)/p_best
            dist_H_test= (1-(Pr_test*Tp/(Tp+1))**(1/3))*Kp*(np.sqrt(pspp1**(2/3)+pspp1**(5/3)*(f2s/f1s)**2+pspp1**(4/3)+pspp1*(f2s/f1s)**2)/2)**(-2/3)#-2*(np.abs(f1s[Tp-1])/3/Tp**(1/2))**(-2/3)*(1-(Pr_test*Tp/(Tp+1))**(1/3))
            Id_closest_test=np.abs(dist_H_test).argmin()
            
            if Tp[Id_closest_test]==p_best:
                sym=True
            else:
                sym=False
        else:
            sym=True

        
        return dist_H[Id_closest],Tp[Id_closest],sym


def dist_Dpp1p(Pr):
    r"""
    Compute the distance to the closest first-order resonance, and verify if it lands in the symetrical domain around that resonance

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

    
def _samples2ell_twoplanets(sample, pair):
      r"""convert the sample into numpy array of orbital elements and masses

      Parameters
      ----------
      sample : np.array or pandaframe
          samples of dynamical analysis of observed system or evolution of simulated system
      pair : tuple of ints
          (i,j) indices of the pair to consider

      Returns
      -------
      np.arrays of floats
          e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2
      """
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
            lbd1 = sample[:,1 + 8*I]
            lbd2 = sample[:,1 + 8*J]
            P1   = sample[:,2 + 8*I]
            P2   = sample[:,2 + 8*J]
            k1   = sample[:,3 + 8*I]
            k2   = sample[:,3 + 8*J]
            h1   = sample[:,4 + 8*I]
            h2   = sample[:,4 + 8*J]
            m1   = sample[:,7 + 8*I]
            m2   = sample[:,7 + 8*J]
      e1   = np.sqrt(k1**2 + h1**2)
      e2   = np.sqrt(k2**2 + h2**2)
      vp1  = np.arctan2(h1, k1)
      vp2  = np.arctan2(h2, k2)
      lbd2 = lbd2*np.pi/180.
      lbd1 = lbd1*np.pi/180.

      return [e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2]

def ell2SFMR(p, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2, return_all=False):    
      r"""
      Converts from elliptic elements to the coordinates :math:`\left(X, Y, X_2, Y_2, \delta\right)` of the Second Fondamental Model of resonance (SFM).
      :math:`\left(X, Y\right)` corresponds to the unique degree of freedom of the SFM and :math:`\left(X_2, Y_2\right)` to the first integral. :math:`\delta` is the unique parameter of the SFM
      The SFM is described byt the Hamiltonian :math:`\mathcal{H}(\Sigma,\sigma)=3(\delta+1)\Sigma-\Sigma^2+2\sqrt{2\Sigma}\cos\sigma` where :math:`X+iY=\sqrt{2\Sigma}e^{i\sigma}`.
      
      Author : Jeremy Couturier. https://jeremycouturier.com
      
      Update 02/10/2026 By A. Leleu to refect the changes in the paper (removed the normalisation by Gamma, and chose \Delta\gamma=0, rename the variables)

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
      P1: float or 1-dimensional numpy array of floats
            The period :math:`2\pi/n_1` of the inner planet in any units. Such that :math:`\mathcal{G}\left(m_0+m_1\right)=n_1^2a_1^3`.
      P2: float or 1-dimensional numpy array of floats
            The period :math:`2\pi/n_2` of the outer planet in the same units as P1. Such that :math:`\mathcal{G}\left(m_0+m_2\right)=n_2^2a_2^3`.
      lbd1: float or 1-dimensional numpy array of floats
            The mean longitude :math:`\lambda_1` of the inner planet, in radians.
      lbd2: float or 1-dimensional numpy array of floats
            The mean longitude :math:`\lambda_2` of the outer planet, in radians.
      return_all: bool
            if return_all==False, the function returns :math:`\left[X, Y, X_2, Y_2, \delta\right]`, if return_all==True it return a dictionary with additional variables instead.
            
      Returns
      -------
      l : List
            The list :math:`\left[X, Y, X_2, Y_2, \delta\right]` of elements of the Second Fondamental Model of Resonance.
      """

      if (p > 20 or p < 1):
            raise Exception('The index p of the resonance must be between 1 and 20 included')


      # Getting semi-major axes and Lambda
      cG     = 4.*np.pi**2
      beta1 = m1/(1. + m1)
      beta2 = m2/(1. + m2)
      mu1   = cG*(1. + m1)
      mu2   = cG*(1. + m2) 
      n1    = 2.*np.pi/P1
      n2    = 2.*np.pi/P2
      a1    = (mu1/n1**2)**(1./3.)
      a2    = (mu2/n2**2)**(1./3.)
      Lbd1  = beta1*np.sqrt(mu1*a1)
      Lbd2  = beta2*np.sqrt(mu2*a2)
      D1    = Lbd1*(1. - np.sqrt(1. - e1**2))
      D2    = Lbd2*(1. - np.sqrt(1. - e2**2))

      #we construct the reference mean motions such that they are at exact commensurability, and Gamma0=Gamma
      Lbd10=((p+1)*Lbd1+p*Lbd2)/((p+1)+p*beta2/beta1*((p+1)/p)**(1/3)*(mu2/mu1)**(2/3) )
      Lbd20=((p+1)*Lbd1+p*Lbd2-(p+1)*Lbd10)/p
      n10=beta1**3*mu1**2/Lbd10**3
      n20=beta2**3*mu2**2/Lbd20**3


      # Getting DG 
      DG    = Lbd1 - Lbd10 + Lbd2 - Lbd20 - D1 - D2
      C1    = 1/Lbd10
      C2    = 1/Lbd20

      #Getting alpha, beta, gamma, delta, R and S
      f1    = f1s[p - 1]
      f2    = f2s[p - 1]
      Rp     = (f1**2*C1*D1 + f2**2*C2*D2 + 2.*f1*f2*np.sqrt(C1*D1*C2*D2)*np.cos(vp1 - vp2))/(f1**2*C1 + f2**2*C2)
      Sp     = (f1**2*C1*D2 + f2**2*C2*D1 - 2.*f1*f2*np.sqrt(C1*D1*C2*D2)*np.cos(vp1 - vp2))/(f1**2*C1 + f2**2*C2)

      
      beta  = 1.5*n10*p*(p*C1 + (p + 1)*C2)
      alpha = -2.*beta*(DG + Sp)
      gamma = m1*n20/C2*np.sqrt(f1**2*C1 + f2**2*C2)
      delta = alpha*(4./(27.*beta*gamma**2))**(1./3.) - 1.



      #Getting X and Y
      K     = (2.*beta/gamma)**(-2./3.)
      
      R = Rp/K
      S = Sp/K
      phi2    = -p*lbd1 + (p + 1)*lbd2
      sig1  = phi2 - vp1
      sig2  = phi2 - vp2
      u1    = np.sqrt(2.*D1)*np.cos(sig1)
      u2    = np.sqrt(2.*D2)*np.cos(sig2)
      v1    = np.sqrt(2.*D1)*np.sin(sig1)
      v2    = np.sqrt(2.*D2)*np.sin(sig2)
      
      
      alpha1 = -f1/np.sqrt(Lbd10)
      alpha2 = -f2/np.sqrt(Lbd20)
      cophi = alpha1/np.sqrt(alpha1**2+alpha2**2)
      siphi = alpha2/np.sqrt(alpha1**2+alpha2**2)


      x1    = cophi*u1 + siphi*u2
      y1    = cophi*v1 + siphi*v2
      x2    = cophi*u2 - siphi*u1
      y2    = cophi*v2 - siphi*v1
      cosr= x1/np.sqrt(2.*Rp)
      sinr= y1/np.sqrt(2.*Rp)
      coss= x2/np.sqrt(2.*Sp)
      sins= y2/np.sqrt(2.*Sp)
      X     = np.sqrt(2.*R)*cosr
      Y     = np.sqrt(2.*R)*sinr
      X2    = np.sqrt(2.*S)*coss
      Y2    = np.sqrt(2.*S)*sins

      if return_all:
            return {'X':X, 
                    'Y':Y,
                    'X2':X2,
                    'Y2':Y2, 
                    'delta':delta, 
                    'Lbd1':Lbd1,
                    'Lbd2':Lbd2, 
                    'Rp':Rp,
                    'Sp':Sp,
                    'cosr':cosr, 
                    'sinr':sinr,         
                    'coss':coss, 
                    'sins':sins,    
                    'R':R,
                    'S':S,   
                    'cophi':cophi,
                    'siphi':siphi,        
                    }
      else:
            return [X, Y, X2, Y2, delta]



def samples2erjesj(samples, pair, p):
    r"""compute the :math:`k=e\cos(\varpi)` and :math:`h=e\sin(\varpi)` components of the resonant and secular eccentricities of a pair of planets

      Parameters
      ----------
      sample : np.array or pandaframe
                samples of dynamical analysis of observed system or evolution of simulated system
      pair : tuple of ints
          (i,j) indices of the pair to consider
      p : int
          p+1:p resonance

      Returns
      -------
      list of np arrays of floats
          kr1,hr1,kr2,hr2,ks1,hs1,ks2,hs2
    """
      
    
    dic=samples2SFMR(samples, pair, p,return_all=True)


    phi2=-p*samples[:,1+8*pair[0]]*pi/180+(p+1)*samples[:,1+8*pair[1]]*pi/180

    r=np.arctan2(dic['sinr'],dic['cosr'])
    s=np.arctan2(dic['sins'],dic['coss'])

    kr1=np.sqrt(1/dic['Lbd1'])*dic['cophi']*np.sqrt(2*dic['Rp'])*np.cos(phi2-r)
    hr1=np.sqrt(1/dic['Lbd1'])*dic['cophi']*np.sqrt(2*dic['Rp'])*np.sin(phi2-r)

    ks1=-np.sqrt(1/dic['Lbd1'])*dic['siphi']*np.sqrt(2*dic['Sp'])*np.cos(phi2-s)
    hs1=-np.sqrt(1/dic['Lbd1'])*dic['siphi']*np.sqrt(2*dic['Sp'])*np.sin(phi2-s)

    kr2=np.sqrt(1/dic['Lbd2'])*dic['siphi']*np.sqrt(2*dic['Rp'])*np.cos(phi2-r)
    hr2=np.sqrt(1/dic['Lbd2'])*dic['siphi']*np.sqrt(2*dic['Rp'])*np.sin(phi2-r)

    ks2=np.sqrt(1/dic['Lbd2'])*dic['cophi']*np.sqrt(2*dic['Sp'])*np.cos(phi2-s)
    hs2=np.sqrt(1/dic['Lbd2'])*dic['cophi']*np.sqrt(2*dic['Sp'])*np.sin(phi2-s)
    
    return kr1,hr1,kr2,hr2,ks1,hs1,ks2,hs2


def X1X2(X, Y, delta):
      r"""Returns :math:`X_1` and :math:`X_2` such that :math:`(X_1, 0)` and :math:`(X2, 0)` are on the same level line as :math:`(X, Y)` for a given value of :math:`\delta` in the SMFR

      Parameters
      ----------
      X : float
          :math:`\sqrt{2R}\cos r`
      Y : float
          :math:`\sqrt{2R}\sin r`
      delta : float
          :math:`\delta` parameter of the SFMR

      Returns
      -------
      _type_
          :math:`X_1` (float), :math:`X_2` (float)
      """
      # 
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

def SFMR2useful(X, Y, X2, Y2, delta):
      r""" Returns [r, R, s, S, x1, x2] where :math:`X+iY = \sqrt{2R}e^{ir}` and :math:`X_2+iY_2 = \sqrt{2S}e^{is}`
      :math:`x_1` and :math:`x_2` are such that :math:`(x1, 0)` and :math:`(x2, 0)` are on the same level line as :math:`(X, Y)`

      Parameters
      ----------
      X : np.array
          Cartesian coordinate :math:`\sqrt{2R}\cos r`
      Y : np.array
          Cartesian coordinate :math:`\sqrt{2R}\sin r`
      X2 : np.array
          Cartesian coordinate :math:`\sqrt{2S}\cos s`
      Y2 : np.array
          Cartesian coordinate :math:`\sqrt{2S}\sin s`
      delta : np.array
          paramter of the SFMR

      Returns
      -------
      list of np.arrays
          [r, R, s, S, x1, x2, RS, RT]
          RS = 1 if formally resonant, 2 for upper circulation, 0 for upper libration, -1 for lower circulation and -2 for lower libration
          RT = True if the pair is in the resonance-trapped domain, False otherwise 
      """

      
      R  = (X**2  + Y**2) /2.
      S = (X2**2 + Y2**2)/2.
      r  = np.arctan2(Y,  X)
      s = np.arctan2(Y2, X2)

      n = len(delta)
      x1 = np.zeros(n)
      x2 = np.zeros(n)
      RS = np.zeros(n)
      RT = np.zeros(n)
      
      
      for i in range(n):
            [xx1, xx2] = X1X2(X[i], Y[i], delta[i])
            x1[i]=xx1
            x2[i]=xx2
            [xmin, xmax, xres, xint, xhyp] = topology(delta[i])
            [lib1,lib2]= X1X2(0, 0, delta[i])
            if (delta[i] < 0.):
                  
                  if  max(xx1,xx2)<max(lib1,lib2):    
                        RS[i]=2 #upper libration
                  else:
                        RS[i]=0 #upper circulation
            else:
                  Hseparatrix = 1.5*(delta[i] + 1.)*xhyp**2 - 0.25*xhyp**4 + 2.*xhyp
                  H           = 1.5*(delta[i] + 1.)*(X[i]**2 + Y[i]**2) - 0.25*(X[i]**2 + Y[i]**2)**2 + 2.*X[i]
                  if (H > Hseparatrix):
                        RS[i]=1 #formally resonant
                  else:
                        if max(xx1,xx2)<xres:
                              if min(xx1,xx2)>min(lib1,lib2):
                                    RS[i]=-2 #lower libration
                              else:
                                    RS[i]=-1 #lower circulation
                        else:
                              RS[i]=0 #upper circulation

            if (RS[i]==1) or ((delta[i]<0) and (delta[i]>-2) and ((max(xx1,xx2)-xres)<1)):
                  RT[i]=True 
            else:
                  RT[i]=False 
      return [r, R, s, S, x1, x2, RS, RT]

def topology(delta):
      r"""Returns :math:`[Xmin, Xmax, Xres, Xint, Xhyp]` where :math:`Xmin, Xmax` are the position of the separatrix, the :math:`Xres` the position of the primary fixed point, 
      :math:`Xint` the position of the secondary fixed point, and :math:`Xhyp` the position of the hyperbolic fixed point, all in the :math:`sin r=0` plane. 
      If :math:`\delta < 0`, returns :math:`[0, 0, Xres, 0, 0]`
      When :math:`\delta >= 0`; then :math:`Xhyp <= Xint <= Xmin <= Xres <= Xmax`

      Parameters
      ----------
      delta : float
          paramter of the SFMR

      Returns
      -------
      list of floats
          [xmin, xmax, xres, xint, xhyp]
      """
      

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
            return [xmin, xmax, xres, xint, xhyp]

def samples2SFMR(sample, pair, p,return_all=False):
      """_summary_

      Parameters
      ----------
      sample : np.array or pandaframe
          samples of dynamical analysis of observed system or evolution of simulated system
      pair : tuple of ints
          (i,j) indices of the pair to consider
      p : int
          indice of the p+1:p resonance
      return_all : bool, optional
          if return_all==False, the function returns :math:`\left[X, Y, X_2, Y_2, \delta\right]`, if return_all==True it return a dictionary with additional variables instead.

      Returns
      -------
      l : List
            The list :math:`\left[X, Y, X_2, Y_2, \delta\right]` of elements of the Second Fondamental Model of Resonance.
      """
      [e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2] = _samples2ell_twoplanets(sample, pair)
      return ell2SFMR(p, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2,return_all=return_all)

def samples2usefull(sample, pair, p):
      r""" Returns [r, R, s, S, x1, x2] where :math:`X+iY = \sqrt{2R}e^{ir}` and :math:`X_2+iY_2 = \sqrt{2S}e^{is}`
      :math:`x_1` and :math:`x_2` are such that :math:`(x1, 0)` and :math:`(x2, 0)` are on the same level line as :math:`(X, Y)`

      Parameters
      ----------
      sample : np.array or pandaframe
      samples of dynamical analysis of observed system or evolution of simulated system
      pair : tuple of ints
      (i,j) indices of the pair to consider
      p : int
      indice of the p+1:p resonance

      Returns
      -------
      list of np.arrays
            [r, R, s, S, x1, x2, RS, RT]
            RS = 1 if formally resonant, 2 for upper circulation, 0 for upper libration, -1 for lower circulation and -2 for lower libration
            RT = True if the pair is in the resonance-trapped domain, False otherwise 
      """
      [X, Y, X2, Y2, delta] = samples2SFMR(sample, pair, p)
      [sig, Sig, sig2, Sig2, x1, x2, RS, RT] = SFMR2useful(X, Y, X2, Y2, delta)
      return [sig, Sig, sig2, Sig2, x1, x2, RS, RT]

class _AttrDict(dict):
      __getattr__ = dict.__getitem__
      __setattr__ = dict.__setitem__

def get_SFMR_quantities(samples,pair,p=None,sampling=1):
    r"""call functions of the package to extract the quantities used in the first paper : delta, DMMR, er, Der, es, RS, RT, p

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
    Dictionary with the following keys:

    delta : np.array of floats of the lengh of the thinned samples
        :math:`\delta` parameter of the second fundamental model for resonance (Henrard and Lemaitre 1983)
    DMMR : np.array of floats of the lengh of the thinned samples
        :math:`\\xi_r^2-3(\delta+1)`
    er : np.array of floats of the lengh of the thinned samples
        :math:`\\xi_r` parameter of the second fundamental model for resonance (see GRSW paper 1 - Leleu et al 2026) also called resonant eccentricities
    Der : np.array of floats of the lengh of the thinned samples
        :math:`\Delta \\xi_r` parameter of the second fundamental model for resonance (see GRSW paper 1 - Leleu et al 2026) 
    es : np.array of floats of the lengh of the thinned samples
        :math:`\\xi_s` secular eccentricities (see GRSW paper 1 - Leleu et al 2026) 
    RS : np.array of floats of the lengh of the thinned samples
        position of the sample in the second fundamental model:
            2 : upper libration
            1 : formally resonant
            0 : upper circulation
            -1 : lower circulation
            -2 : lower libration
    p : int
        pair attributed to the p+1:p MMR
    
    """


    if isinstance(samples, pd.DataFrame):
            P1=np.mean(samples['period_days_'+str(pair[0])].values)
            P2=np.mean(samples['period_days_'+str(pair[1])].values)
            samples=samples.to_numpy()#np.vstack([samples[col] for col in samples.columns])
    else:
            P1=np.mean(samples[:,2+pair[0]*8])
            P2=np.mean(samples[:,2+pair[1]*8])


    if p==None:
        dist,p,sym=dist_Dpp1p(P2/P1)
    
    [X, Y, X2, Y2, delta]=samples2SFMR( samples[::sampling,:], pair, p)
    [sig, Sig, sig2, Sig2, x1, x2, RS, RT] = samples2usefull( samples[::sampling,:], pair,p)
    
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

    dic={
          'delta':delta,
          'DMMR':DMMR,
          'er':er,
          'es':es,
          'Der':Der,
          'RS':RS,
          'RT':RT,
          'p':p
    }
    

    return _AttrDict(dic)




#Plotting

def plot_topology(ax1, delta_lim=None, X_lim=None, linewidth=4, alpha=1, grid=False,dark=False,alpha_fill=0.1,legend=True):
      r"""Plots the topology of the phase space (separatrices and fixed points) of the Second Fundamental Model of Resonance  

      Parameters
      ----------
      ax1 : ax
          ax of the figure
      delta_lim : tuple of floats, optional
          (:math: `\delta_{min},\delta_{max}`) boundaries of the figure :math: `\delta`, by default None
      X_lim :  tuple of floats, optional
          (:math: `\\xi_r_{min},\xi_r_{max}`) boundaries of the figure :math: `\\xi_r`, by default None
      linewidth : int, optional
          _description_, by default 4
      alpha : int, optional
          _description_, by default 1
      grid : bool, optional
          _description_, by default False
      dark : bool, optional
          _description_, by default False
      alpha_fill : float, optional
          _description_, by default 0.1
      legend : bool, optional
          _description_, by default True
      """
      

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
      ax1.set_ylabel(ylabel=r"$\xi_r$",    labelpad = 4, rotation = 0)
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


def plot_topology_DMMR(ax1, DMMR_lim=None, X_lim=None, linewidth=4, alpha=1, grid=False,dark=False,legend=True,zorder=1,resonance_trapped=False):
      r"""Plots the topology of the phase space (separatrices and fixed points) of the Second Fundamental Model of Resonance  

      Parameters
      ----------
      ax1 : ax
          ax of the figure
      DMMR_lim : tuple of floats, optional
          (:math:`\Delta_{MMR,min},\Delta_{MMR,max}`) boundaries of the figure for the x axis :math:`\Delta_{MMR}`, by default None
      X_lim :  tuple of floats, optional
          (:math:`\xi_{r,min},\xi_{r,max}`) boundaries of the figure for the y axis :math:`\xi_r`, by default None
      linewidth : int, optional
          _description_, by default 4
      alpha : int, optional
          _description_, by default 1
      grid : bool, optional
          _description_, by default False
      dark : bool, optional
          _description_, by default False
      alpha_fill : float, optional
          _description_, by default 0.1
      legend : bool, optional
          _description_, by default True
      """

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
      ax1.set_ylabel(ylabel=r"$\xi_r$",          labelpad = 4, rotation = 0)
      #delt = np.linspace(delta_lim[0]-X2max, delta_lim[1]+X2max, 512)
      delt = np.arange(-bound, bound, .01)
      Xmin = np.zeros(delt.size)
      Xmax = np.zeros(delt.size)
      Xres = np.zeros(delt.size)
      Xint = np.zeros(delt.size)
      Xhyp = np.zeros(delt.size)
      if resonance_trapped:
            Xrestrapped = np.ones(delt.size)*np.nan
      count = 0
      for count,delta in enumerate(delt):
            [xmin, xmax, xres, xint, xhyp] = topology(delta)
            Xmin[count] = xmin
            Xmax[count] = xmax
            Xres[count] = xres
            Xint[count] = xint
            Xhyp[count] = xhyp
            if resonance_trapped:
                  if delta >=-2 and delta<0:
                        Xrestrapped[count] = xres+1
                        if count< Xrestrapped.size:
                              Xrestrapped[count+1] = xres

      if resonance_trapped:
            Id_b=np.searchsorted(delt,-2)
            Id_e=np.searchsorted(delt,0)

            delt_restrap=delt[Id_b:Id_e]
            X_restrap=Xres[Id_b:Id_e]+1
            X_restrap=np.append(np.linspace(X_restrap[0]-1,X_restrap[0],20),X_restrap)
            delt_restrap=np.append(np.ones(20)*-2,delt_restrap)
  
                              
                  
      if dark:        
            
            if legend:
                  ax1.plot(Xint[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xint[delt >= 0.], color = 'grey', linewidth = linewidth, zorder=zorder, linestyle = '-', alpha = alpha, label = 'secondary')
                  ax1.plot(Xhyp[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xhyp[delt >= 0.], color = 'pink',   linewidth = linewidth, zorder=zorder, linestyle = ':', alpha = alpha, label = 'hyperbolic')
                  ax1.plot(Xres**2-3*(delt + 1.),             Xres,             color = 'white', linewidth = linewidth, zorder=zorder, linestyle = '-', alpha = alpha, label = 'primary')
                  ax1.plot(Xmin[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xmin[delt >= 0.], color = 'pink',   linewidth = linewidth, zorder=zorder, linestyle = '-', alpha = alpha, label = 'separatrix')
            else:
                  ax1.plot(Xint[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xint[delt >= 0.], color = 'grey', linewidth = linewidth, zorder=zorder, linestyle = '-', alpha = alpha)
                  ax1.plot(Xhyp[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xhyp[delt >= 0.], color = 'pink',   linewidth = linewidth, zorder=zorder, linestyle = ':', alpha = alpha)
                  ax1.plot(Xres**2-3*(delt + 1.),             Xres,             color = 'white', linewidth = linewidth, zorder=zorder, linestyle = '-', alpha = alpha)
                  ax1.plot(Xmin[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xmin[delt >= 0.], color = 'pink',   linewidth = linewidth, zorder=zorder, linestyle = '-', alpha = alpha)
            ax1.plot(Xmax[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xmax[delt >= 0.], color = 'pink',   linewidth = linewidth, zorder=zorder, linestyle = '-', alpha = alpha)
            #ax1.fill_between(delt[delt >= 0.], Xmin[delt >= 0.], Xmax[delt >= 0.], color = 'pink', alpha = alpha_fill)
      else:
            
            if legend:
                  ax1.plot(Xint[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xint[delt >= 0.], color = 'grey', linewidth = linewidth, zorder=zorder, linestyle = '-', alpha = alpha, label = 'secondary')
                  ax1.plot(Xhyp[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xhyp[delt >= 0.], color = 'red',   linewidth = linewidth, zorder=zorder, linestyle = ':', alpha = alpha, label = 'hyperbolic')
                  ax1.plot(Xres**2-3*(delt + 1.),             Xres,             color = 'black', linewidth = linewidth, zorder=zorder, linestyle = '-', alpha = alpha, label = 'primary')
                  ax1.plot(Xmin[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xmin[delt >= 0.], color = 'red',   linewidth = linewidth, zorder=zorder, linestyle = '-', alpha = alpha, label = 'separatrix')
                  if resonance_trapped:
                        ax1.plot(X_restrap**2-3*(delt_restrap + 1.),X_restrap, color = 'red',   linewidth = linewidth, zorder=zorder, linestyle = '--', alpha = alpha, label = 'resonance-\ntrapped')
            else:
                  ax1.plot(Xint[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xint[delt >= 0.], color = 'grey', linewidth = linewidth, zorder=zorder, linestyle = '-', alpha = alpha)
                  ax1.plot(Xhyp[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xhyp[delt >= 0.], color = 'red',   linewidth = linewidth, zorder=zorder, linestyle = ':', alpha = alpha)
                  ax1.plot(Xres**2-3*(delt + 1.),             Xres,             color = 'black', linewidth = linewidth, zorder=zorder, linestyle = '-', alpha = alpha)
                  ax1.plot(Xmin[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xmin[delt >= 0.], color = 'red',   linewidth = linewidth, zorder=zorder, linestyle = '-', alpha = alpha)   
                  if resonance_trapped:
                        ax1.plot(X_restrap**2-3*(delt_restrap + 1.),X_restrap, color = 'red',   linewidth = linewidth, zorder=zorder, linestyle = '--', alpha = alpha) 
            ax1.plot(Xmax[delt >= 0.]**2-3*(delt[delt >= 0.] + 1.), Xmax[delt >= 0.], color = 'red',   linewidth = linewidth, zorder=zorder, linestyle = '-', alpha = alpha)
            #ax1.fill_between(delt[delt >= 0.], Xmin[delt >= 0.], Xmax[delt >= 0.], color = 'red', alpha = alpha_fill)

      if grid:
            ax1.grid(linewidth=0.3, alpha = 0.5)


def plot_samples_SFMR_pair( ax1, sample, pair, p, color, color_lim=None, label='',alpha = 1,markersize=20,marker='o',flag_DMMR=True):
      delta, DMMR, er, Der, es, RS, RT, p= get_SFMR_quantities(sample,pair,p=p)
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


def plot_samples_SFMR( ax1, sample, pairs, p_indexes, colors, color_lim=None, label='',alpha = 1,markersize=20,marker='o',flag_DMMR=True):
      ### Plots the sample in the phase space of the Second Fundamental Model ###

      if isinstance(p_indexes, (list, np.ndarray)):
            if isinstance(colors, (list, np.ndarray)):
                  colors = cycle(colors)
            for pair, p, color in zip(pairs, p_indexes, colors):

                  plot_samples_SFMR_pair( ax1, sample, pair, p, color, color_lim=color_lim, label=label,alpha =alpha,markersize=markersize,marker=marker,flag_DMMR=flag_DMMR)
                           
      else:
                  plot_samples_SFMR_pair( ax1, sample, pairs, p_indexes, colors, color_lim=color_lim, label=label,alpha =alpha,markersize=markersize,marker=marker,flag_DMMR=flag_DMMR)
                  

