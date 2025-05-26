                                                ###### Conversion of samples ######

# We provide here 6 conversion functions for the samples of the workshop.
# A row of sample is : Id or timestamp, lbd(°), Period(days), k, h, Inclination(°), Omega(°), m/M*, R/R*, ... , M* (Msun), R* (Rsun)
# Columns Id or timestamp, mj/Mj*, Rj/Rj*, M* (Msun), R* (Rsun) are NOT modified by these functions.
# Only columns lbd(°), Period(days), k, h, Inclination(°), Omega(°) are modified for each planet.


# - Sample2cart(sample, typeOfCoordinates, adcol)   -->  Converts the sample to cartesian coordinates.
#                                                        The columns are now Id or timestamp, X, Y, Z, vX, vY, vZ, m/M*, R/R*, ... , M* (Msun), R* (Rsun)
#                                                        typeOfCoordinates is either 'Jacobi', 'JacobiWisdomHolman' or 'Heliocentric'
#                                                        Jacobi differs from JacobiWisdomHolman in the sense that the mu is defined differently. See Wisdom & Holman 1991
#                                                        adcol is the number of additional columns in the sample (often 0 unless otherwise specified by the sample's author)

# - Sample2aeiMoO(sample, typeOfCoordinates, adcol) -->  Converts the sample to elliptic coordinates (a, e, i, M, omega, Omega) = (semi-major axis, eccentricity, inclination, 
#                                                        mean anomaly, argument of periapsis, longitude of ascending node).
#                                                        The columns are now Id or timestamp, a, e, i(rad), M(rad), omega(rad), Omega(rad), m/M*, R/R*, ... , M* (Msun), R* (Rsun)
#                                                        typeOfCoordinates is either 'Jacobi', 'JacobiWisdomHolman' or 'Heliocentric'
#                                                        adcol is the number of additional columns in the sample (often 0 unless otherwise specified by the sample's author)

# - Sample2alkhqp(sample, typeOfCoordinates, adcol) -->  Converts the sample to elliptic coordinates (a, lbd, k, h, q, p) = (semi-major axis, mean longitude,
#                                                        e*cos(varpi), e*sin(varpi), sin(i/2)*cos(Omega), sin(i/2)*sin(Omega)).
#                                                        The columns are now Id or timestamp, a, lbd(rad), k, h, q, p, m/M*, R/R*, ... , M* (Msun), R* (Rsun)
#                                                        typeOfCoordinates is either 'Jacobi', 'JacobiWisdomHolman' or 'Heliocentric'
#                                                        adcol is the number of additional columns in the sample (often 0 unless otherwise specified by the sample's author)

# - Cart2sample(sample, typeOfCoordinates, adcol)   -->  Inverse of Sample2cart.

# CAUTION : typeOfCoordinates indicates the type of coordinates the sample is given in. It cannot be used to change from Jacobi to Heliocentric or vice-versa.
#           For example, if your sample is in Jacobi coordinates, Sample2cart(sample, 'Jacobi') will convert it into Jacobi cartesian coordinates, whereas
#           Sample2cart(sample, 'Heliocentric') will have undefined behavior. Use functions Jac2Hel and Hel2Jac for conversion Jacobi <--> Heliocentric

# - Jac2Hel(sample, adcol)   --> Converts the sample from Jacobi coordinates to Heliocentric coordinates. The regular Jacobi convention mu = G(M* + m1 + ... + mj) is used
#                                adcol is the number of additional columns in the sample (often 0 unless otherwise specified by the sample's author)

# - Hel2Jac(sample, adcol)   --> Converts the sample from Heliocentric coordinates to Jacobi coordinates. The regular Jacobi convention mu = G(M* + m1 + ... + mj) is used
#                                adcol is the number of additional columns in the sample (often 0 unless otherwise specified by the sample's author)

# - JacWH2Hel(sample, adcol) --> Converts the sample from Jacobi coordinates to Heliocentric coordinates. The WisdomHolman convention mu = GM*(M* + m1 + ... + mj)/(M* + m1 + ... + m_{j-1}) is used
#                                adcol is the number of additional columns in the sample (often 0 unless otherwise specified by the sample's author)

# - Hel2JacWH(sample, adcol) --> Converts the sample from Heliocentric coordinates to Jacobi coordinates. The WisdomHolman convention mu = GM*(M* + m1 + ... + mj)/(M* + m1 + ... + m_{j-1}) is used
#                                adcol is the number of additional columns in the sample (often 0 unless otherwise specified by the sample's author)

# These functions convert into a system of units where : (see https://iau-a3.gitlab.io/NSFA/IAU2009_consts.html)
# - The unit of length is the astronomical unit (1.49597870700e11 m, as recommended by the IAU)
# - The unit of mass   is the mass of the Sun (1.98841583e30 kg, taking G and G*Msun as recommended by the IAU)
# - The unit of time   is the day (86400 s, as recommended by the IAU)


# Author : Jérémy COUTURIER

                                         ###### See the very bottom of this page for examples of use ######

import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib
import numpy as np
import random


##########################################################################################################################
##########################################################################################################################
##########################################################################################################################


#Defining the gravitational constant
G = 0.0002959122 # AU^3 Msun^-1 day^-2 (Using the value G=6.67428e-11 m^3 kg^-1 s^-2 recommended by the IAU)


def ell2cart_true(aeinuoO, mass):

      #Returns the cartesian coordinates. This is an auxiliary function

      #aeinuoO = [a, e, i, nu, omega, Omega] = [semi-major axis, eccentricity, inclination, true anomaly, argument of the periapsis, longitude of the ascending node]
      #i, nu, omega, Omega are in radians
      #mass is the central mass. e.g. mass = M* + m_j for heliocentric, M* + m_1 + ... + m_j for regular Jacobi and M*(M* + m1 + ... + mj)/(M* + m1 + ... + m_{j-1}) for WisdomHolman Jacobi
      
      #Extracted from NcorpiON software and translated from C to python
      
      a        = aeinuoO[0]
      e        = aeinuoO[1]
      i        = aeinuoO[2]
      nu       = aeinuoO[3]
      omega    = aeinuoO[4]
      Omega    = aeinuoO[5]
      
      mu       = G*mass
      cosnu    = np.cos(nu)
      sinnu    = np.sin(nu)
      cosvarpi = np.cos(omega + Omega)
      sinvarpi = np.sin(omega + Omega)
      q        = np.sin(i/2.)*np.cos(Omega)
      p        = np.sin(i/2.)*np.sin(Omega)
      chi      = np.cos(i/2.)
      pp       = 1. - 2.*p*p
      qq       = 1. - 2.*q*q
      dpq      = 2.*p*q
      
      ## In the orbital plane (see e.g. Laskar & Robutel 1995) ##
      r       = a*(1. - e*e)/(1. + e*cosnu)
      g       = np.sqrt(mu*a*(1. - e*e))
      dnudt   = g/(r*r)
      drdt    = a*e*dnudt*sinnu*(1. - e*e)/((1. + e*cosnu)*(1. + e*cosnu))
      X_buff  = r*cosnu
      Y_buff  = r*sinnu
      vX_buff = drdt*cosnu - r*dnudt*sinnu
      vY_buff = drdt*sinnu + r*dnudt*cosnu

      ## Rotations to convert to reference plane (see e.g. Laskar & Robutel 1995) ##      
      X  =  X_buff*(pp*cosvarpi + dpq*sinvarpi) +  Y_buff*(dpq*cosvarpi - pp*sinvarpi)
      vX = vX_buff*(pp*cosvarpi + dpq*sinvarpi) + vY_buff*(dpq*cosvarpi - pp*sinvarpi)
      Y  =  X_buff*(qq*sinvarpi + dpq*cosvarpi) +  Y_buff*(qq*cosvarpi  - dpq*sinvarpi)
      vY = vX_buff*(qq*sinvarpi + dpq*cosvarpi) + vY_buff*(qq*cosvarpi  - dpq*sinvarpi)
      Z  =  X_buff*(2.*q*chi*sinvarpi - 2.*p*chi*cosvarpi) +  Y_buff*(2.*p*chi*sinvarpi + 2.*q*chi*cosvarpi)
      vZ = vX_buff*(2.*q*chi*sinvarpi - 2.*p*chi*cosvarpi) + vY_buff*(2.*p*chi*sinvarpi + 2.*q*chi*cosvarpi)
      
      return [X, Y, Z, vX, vY, vZ]
      
      
def ell2cart_eccentric(aeiEoO, mass):

      #Returns the cartesian coordinates. This is an auxiliary function

      #aeiEoO = [a, e, i, E, omega, Omega] = [semi-major axis, eccentricity, inclination, eccentric anomaly, argument of the periapsis, longitude of the ascending node]
      #i, E, omega, Omega are in radians
      #mass is the central mass. e.g. mass = M* + m_j for heliocentric, M* + m_1 + ... + m_j for regular Jacobi and M*(M* + m1 + ... + mj)/(M* + m1 + ... + m_{j-1}) for WisdomHolman Jacobi
      
      #Extracted from NcorpiON software and translated from C to python
      
      a        = aeiEoO[0]
      e        = aeiEoO[1]
      i        = aeiEoO[2]
      E        = aeiEoO[3]
      omega    = aeiEoO[4]
      Omega    = aeiEoO[5]
      
      mu       = G*mass
      cosE     = np.cos(E)
      sinE     = np.sin(E)
      cosvarpi = np.cos(omega + Omega)
      sinvarpi = np.sin(omega + Omega)
      q        = np.sin(i/2.)*np.cos(Omega)
      p        = np.sin(i/2.)*np.sin(Omega)
      chi      = np.cos(i/2.)
      pp       = 1. - 2.*p*p
      qq       = 1. - 2.*q*q
      dpq      = 2.*p*q
      
      ## In the orbital plane (see e.g. Laskar & Robutel 1995) ##
      na      = np.sqrt(mu/a)
      X_buff  = a*(cosE - e)
      Y_buff  = a*sinE*np.sqrt(1. - e*e)
      vX_buff = -na*sinE/(1. - e*cosE)
      vY_buff = na*np.sqrt(1. - e*e)*cosE/(1. - e*cosE)

      ## Rotations to convert to reference plane (see e.g. Laskar & Robutel 1995) ##      
      X  =  X_buff*(pp*cosvarpi + dpq*sinvarpi) +  Y_buff*(dpq*cosvarpi - pp*sinvarpi)
      vX = vX_buff*(pp*cosvarpi + dpq*sinvarpi) + vY_buff*(dpq*cosvarpi - pp*sinvarpi)
      Y  =  X_buff*(qq*sinvarpi + dpq*cosvarpi) +  Y_buff*(qq*cosvarpi  - dpq*sinvarpi)
      vY = vX_buff*(qq*sinvarpi + dpq*cosvarpi) + vY_buff*(qq*cosvarpi  - dpq*sinvarpi)
      Z  =  X_buff*(2.*q*chi*sinvarpi - 2.*p*chi*cosvarpi) +  Y_buff*(2.*p*chi*sinvarpi + 2.*q*chi*cosvarpi)
      vZ = vX_buff*(2.*q*chi*sinvarpi - 2.*p*chi*cosvarpi) + vY_buff*(2.*p*chi*sinvarpi + 2.*q*chi*cosvarpi)
      
      return [X, Y, Z, vX, vY, vZ]
      

def cart2ell(XYZvXvYvZ, mass):

      #Converts the cartesian coordinates into the elliptic elements [a, lambda, k, h, q, p] where k + ih = e*exp(i*varpi) and q + ip = sin(I/2)*exp(i*Omega)
      #mass is the central mass. e.g. mass = M* + m_j for heliocentric, M* + m_1 + ... + m_j for regular Jacobi and M*(M* + m1 + ... + mj)/(M* + m1 + ... + m_{j-1}) for WisdomHolman Jacobi
      
      #This function was originally provided by Mickaël Gastineau of the IMMCE lab for NcorpiON software and converted from C to python for the workshop

      mu = G*mass

      ##### Getting the cartesian coordinates #####
      X   = XYZvXvYvZ[0]
      Y   = XYZvXvYvZ[1]
      Z   = XYZvXvYvZ[2]
      vX  = XYZvXvYvZ[3]
      vY  = XYZvXvYvZ[4]
      vZ  = XYZvXvYvZ[5]

      ##### Computing the semi-major axis #####
      R  = np.sqrt(X*X + Y*Y + Z*Z)
      V2 = vX*vX + vY*vY + vZ*vZ
      AA = R*mu/(2.*mu - R*V2) #Division by zero if the trajectory is perfectly parabolic.

      ##### Normalizing the velocities (Adopting the convention of J. Laskar's 2004 lectures notes) #####
      SMU = np.sqrt(mu)
      vX  = vX/SMU
      vY  = vY/SMU
      vZ  = vZ/SMU

      ##### Computing the angular momentum #####
      V2 = vX*vX + vY*vY + vZ*vZ
      RV = X*vX  + Y*vY  + Z*vZ
      C1 = Y*vZ  - Z*vY
      C2 = Z*vX  - X*vZ
      C3 = X*vY  - Y*vX
      CC = C1*C1 + C2*C2 + C3*C3
      DC = np.sqrt(CC)

      ##### Computing (q, p) #####
      aux0 = np.sqrt(2.*(CC + DC*C3))
      if (aux0 == 0.):
            q = 0.
            p = 0.
      else:
            q = -C2/aux0
            p =  C1/aux0

      ##### Computing (k, h) #####
      if (R == 0. or V2 == 0.):
            K = 0.
            H = 0.

      else:
            ##### Computing the matrix coefficients needed for (k, h) #####
            a11 = V2 - 1./R
            a12 = RV/(R*DC)
            a21 = -RV
            a22 = DC - R/DC

            ##### Computing (k, h) #####
            c11  =  X*a11  + vX*a21
            c12  =  X*a12  + vX*a22
            c21  =  Y*a11  + vY*a21
            c22  =  Y*a12  + vY*a22
            FAC1 =  C1/(DC + C3)
            FAC2 =  C2/(DC + C3)
            K1   =  c11 - FAC1*(Z*a11 + vZ*a21)
            H1   = -c12 + FAC1*(Z*a12 + vZ*a22)
            H2   =  c21 - FAC2*(Z*a11 + vZ*a21) #Should be equal to H1
            K2   =  c22 - FAC2*(Z*a12 + vZ*a22) #Should be equal to K1
            if (abs(H1 - H2) + abs(K1 - K2) > 1.e-6):
                  print("Warning : Bad computation of (k,h) in function cart2ell. (K1 - K2, H1 - H2) = ", K1 - K2, H1 - H2)
            K    =  0.5*(K1 + K2)
            H    =  0.5*(H1 + H2)

      ##### Computing lambda #####
      if (R == 0. or V2 == 0.):
            lbd = 0.

      else:
            ##### Computing the mean longitude l = M + varpi #####
            USQA = np.sqrt(2./R - V2)
            if ((USQA) >= 0.): #Elliptic case
                  b12  = vX - FAC1*vZ
                  b22  = vY - FAC2*vZ
                  aux1 = (R*V2 - 1.)/(1. + DC*USQA)
                  sinF = -b12*R*USQA + H*aux1
                  cosF =  b22*R*USQA + K*aux1
                  F    =  np.arctan2(sinF, cosF)
                  lbd  = F - RV*USQA

            else: #Hyperbolic case
                  USQA = np.sqrt(-(2./R - V2))
                  E    = np.arctanh(RV*USQA/(R*V2 - 1.))
                  M    = np.sqrt(K*K + H*H) * np.sinh(E) - E
                  lbd  = M + np.arctan2(H, K)

      return [AA, lbd, K, H, q, p]


def get_true_anomaly(aeiMoO, mass):

      ##### Computes the true anomaly from the differential equation                         #####
      ##### dnu/dt = sqrt(mu)*(1 + e cos nu)^2/(a(1-e^2))^(3/2) using a Runge-Kutta 4 method #####
      ##### This is equivalent to solving Kepler's equation.                                 #####
      
      #mass is the central mass. e.g. mass = M* + m_j for heliocentric, M* + m_1 + ... + m_j for regular Jacobi and M*(M* + m1 + ... + mj)/(M* + m1 + ... + m_{j-1}) for WisdomHolman Jacobi

      #Extracted from NcorpiON software and translated from C to python

      a     = aeiMoO[0]
      e     = aeiMoO[1]
      i     = aeiMoO[2]
      M     = aeiMoO[3]
      o     = aeiMoO[4]
      O     = aeiMoO[5]

      mu    = G*mass
      sq_mu = np.sqrt(mu)
      denom = (a*(1. - e*e))**1.5
      n     = np.sqrt(mu/abs(a)**3)
      M     = np.mod(M, 2.*np.pi)
      if (M < 0.):
            M = M + 2.*np.pi
      time  = M/n
      
      previous_tra = 0.
      period       = 2.*np.pi/n
      t            = time
      if (e < 1.):
            N_step = np.floor(512.*t/period) + 1
      else:
            N_step = 512
      if (e > 0.8 and e < 1.):
            N_step = N_step*4
      dt = t/N_step
      N_step = int(N_step)
      
      ##### Integrating #####
      for j in range(N_step):
            partial_tra   = previous_tra
            num           = 1. + e*np.cos(partial_tra)
            K1            = sq_mu*num*num/denom
            partial_tra   = previous_tra + 0.5*K1*dt
            num           = 1. + e*np.cos(partial_tra)
            K2            = sq_mu*num*num/denom
            partial_tra   = previous_tra + 0.5*K2*dt
            num           = 1. + e*np.cos(partial_tra)
            K3            = sq_mu*num*num/denom
            partial_tra   = previous_tra + K3*dt
            num           = 1. + e*np.cos(partial_tra)
            K4            = sq_mu*num*num/denom
            previous_tra  = previous_tra + dt*(K1 + 2.*K2 + 2.*K3 + K4)/6.
      
      nu = previous_tra
      return [a, e, i, nu, o, O]


def l2F(l,k,h):
      #Solves the Kepler equation without a numerical integration. Much faster than the above function
      #Returns F = E + varpi
      #(l, k, h) = (Mean longitude, e*cos(varpi), e*sin(varpi))
      #Written by ASD Team LTE lab (former IMCCE)

      eps = 2*2.26e-16
      imax = 20
      # depart methode d'ordre 3
      a=l
      ca=np.cos(a)
      sa=np.sin(a)
      se=k*sa-h*ca
      ce=k*ca+h*sa
      fa=a-se-l
      f1a=1.0-ce
      f2a=se/2.0
      f3a=ce/6.0
      d1=-fa/f1a
      d2=-fa/(f1a-d1*f2a)
      d3 =-fa/(f1a+d2*(f2a+d2*f3a))
      a=a+d3
      # methode d'ordre 6
      ca=np.cos(a)
      sa=np.sin(a)
      se=k*sa-h*ca
      ce=k*ca+h*sa
      fa=a-se-l
      f1a=1.0-ce
      f2a=se/2.0
      f3a=ce/6.0
      f4a=-se/24.0
      f5a=-ce/120.0
      d1=-fa/f1a
      d2=-fa/(f1a-d1*f2a)
      d3=-fa/(f1a+d2*(f2a+d2*f3a))
      d4=-fa/(f1a+d3*(f2a+d3*(f3a+d3*f4a)))
      d5=-fa/( f1a+d4*(f2a+d4*(f3a+d4*(f4a+d4*f5a))))
      a=a+d5
      #     return
      #     on teste la precision obtenue
      i=0
      while True:
            i=i+1
            ca=np.cos(a)
            sa=np.sin(a)
            se=k*sa-h*ca
            fa=a-se-l
            ce=k*ca+h*sa
            f1a=1.0-ce
            d1=-fa/f1a
            #     si la precison n'est pas bonne, on continue les calculs
            #     en iterant la methode d'ordre 1
            if (abs(d1)/max(1.0,abs(a)) > eps):
                  if (i > imax):
                        #     write(*,*) 'erreur fatale dans elliptid:keplkh2'
                        #     write(*,*) 'erreur :',abs(d1)/dmax1(1.0,abs(a))
                        return(a)
                  a=a+d1
            else:
                  return(a)


def lPkhiO_to_aeiMoO(lPkhiO, mass, daysInUOT):

      #Converts from the workshop format to elliptic elements [semi-major axis, eccentricity, inclination, mean anomaly, argument of periapsis, longitude of ascending node]
      #lPkhiO = [mean longitude (deg), Period (days), k, h, inclination(deg), longitude of the ascending node(deg)]
      #mass is the central mass. e.g. mass = M* + m_j for heliocentric, M* + m_1 + ... + m_j for regular Jacobi and M*(M* + m1 + ... + mj)/(M* + m1 + ... + m_{j-1}) for WisdomHolman Jacobi

      
      mu    = G*mass
      lbd   = lPkhiO[0]
      P     = lPkhiO[1]
      k     = lPkhiO[2]
      h     = lPkhiO[3]
      i     = lPkhiO[4]
      Omega = lPkhiO[5]
      
      lbd   = lbd  *m.pi/180.
      i     = i    *m.pi/180.
      Omega = Omega*m.pi/180.
      P     = P/daysInUOT

      e     = np.sqrt(h**2 + k**2)
      varpi = np.arctan2(h, k)
      M     = lbd - varpi
      omega = varpi - Omega
      a     = (mu*P**2/(4.*np.pi**2))**(1./3.)

      return [a, e, i, M, omega, Omega]
      
      
def aeiMoO_to_alkhqp(aeiMoO, mass, daysInUOT):

      #Converts from elliptic elements [semi-major axis, eccentricity, inclination, mean anomaly, argument of periapsis, longitude of ascending node]
      #to elliptic elements [semi-major axis, mean longitude (deg), k, h, q, p]
      #mass is the central mass. e.g. mass = M* + m_j for heliocentric, M* + m_1 + ... + m_j for regular Jacobi and M*(M* + m1 + ... + mj)/(M* + m1 + ... + m_{j-1}) for WisdomHolman Jacobi

      
      mu    = G*mass
      a     = aeiMoO[0]
      e     = aeiMoO[1]
      i     = aeiMoO[2]
      M     = aeiMoO[3]
      omega = aeiMoO[4]
      Omega = aeiMoO[5]
      
      varpi = omega + Omega
      k     = e*np.cos(varpi)
      h     = e*np.sin(varpi)
      q     = np.sin(i/2.)*np.cos(Omega)
      p     = np.sin(i/2.)*np.sin(Omega)
      lbd   = M + varpi

      return [a, lbd, k, h, q, p]


def lPkhiO_to_cart(lPkhiO, mass, daysInUOT):

      #Converts from the workshop format to cartesian coordinates
      #lPkhiO = [mean longitude (deg), Period (days), k, h, inclination(deg), longitude of the ascending node(deg)]
      #mass is the central mass. e.g. mass = M* + m_j for heliocentric, M* + m_1 + ... + m_j for regular Jacobi and M*(M* + m1 + ... + mj)/(M* + m1 + ... + m_{j-1}) for WisdomHolman Jacobi

      aeiMoO = lPkhiO_to_aeiMoO(lPkhiO, mass, daysInUOT)
      a = aeiMoO[0]
      e = aeiMoO[1]
      i = aeiMoO[2]
      M = aeiMoO[3]
      o = aeiMoO[4]
      O = aeiMoO[5]
      #aeinuoO = get_true_anomaly(aeiMoO, mass)
      #cart    = ell2cart_true(aeinuoO, mass)
      vrp = o + O
      l   = M + vrp
      k   = e*np.cos(vrp)
      h   = e*np.sin(vrp)
      aeiEoO  = [a, e, i, l2F(l,k,h) - vrp, o, O]
      cart    = ell2cart_eccentric(aeiEoO, mass)

      return cart;


def sample2cart_row(row, typeOfCoordinates, adcol):

      #Converts one row of the workshop into the same row in cartesian coordinates
      #typeOfCoordinates is either 'Jacobi', 'JacobiWisdomHolman' or 'Heliocentric'
      #typeOfCoordinates must be the type of coordinates into which the sample is given !
      #This function cannot be used to change the type of coordinates ! Its purpose is only to convert to cartesian
      #Use functions Jac2Hel_row or Hel2Jac_row to change the type of coordinates of the sample row
      #adcol is the number of additional columns in the sample (often 0). (e.g. adcol = 1 for Hadden's posteriors)

      n = len(row)
      if (n%8 != 3 + adcol):
            raise Exception("The sample line does not have 8*k + 3 + adcol columns.")

      N = (len(row) - 3 - adcol)//8 #Number of planets
      
      output = np.array([row[0]])
      
      daysInUOT = 1.

      for i in range(N):
            lPkhiO = row[1 + 8*i : 7 + 8*i]
            if (typeOfCoordinates == 'Jacobi'):
                  mass = 1.
                  for j in range(i + 1):
                        mass = mass + row[7 + 8*j]
            elif (typeOfCoordinates == 'Heliocentric'):
                  mass = 1. + row[7 + 8*i]
            elif (typeOfCoordinates == 'JacobiWisdomHolman'):
                  massum = 1. # M* + m1 + ... + mi
                  for j in range(i + 1):
                        massum = massum + row[7 + 8*j]
                  mass = massum/(massum - row[7 + 8*i]) # M*(M* + m1 + ... + mi)/(M* + m1 + ... + m_{i-1})
            else:
                  raise Exception("typeOfCoordinates can be either 'Jacobi' or 'Heliocentric' or 'JacobiWisdomHolman'.")
            mass = mass*row[-2-adcol]
            cart   = np.array(lPkhiO_to_cart(lPkhiO, mass, daysInUOT))
            output = np.concatenate((output, cart, np.array([row[7 + 8*i], row[8 + 8*i]])))

      output = np.concatenate((output, np.array([row[-2-adcol], row[-1-adcol]])))
      for i in range(adcol):
            output = np.concatenate((output, np.array([row[-adcol+i]])))
      return output


def sample2aeiMoO_row(row, typeOfCoordinates, adcol):

      #Converts one row of the workshop into the same row in coordinates a, e, i, M, o, O
      #[a, e, i, M, o, O] = [semi-major axis, eccentricity, inclination, mean anomaly, argument of periapsis, longitude of ascending node]
      #typeOfCoordinates is either 'Jacobi', 'JacobiWisdomHolman' or 'Heliocentric'
      #typeOfCoordinates must be the type of coordinates into which the sample is given !
      #This function cannot be used to change the type of coordinates ! Its purpose is only to convert to a, e, i, M, o, O
      #Use functions Jac2Hel_row or Hel2Jac_row to change the type of coordinates of the sample row
      #adcol is the number of additional columns in the sample (often 0). (e.g. adcol = 1 for Hadden's posteriors)

      n = len(row)
      if (n%8 != 3 + adcol):
            raise Exception("The sample line does not have 8*k + 3 + adcol columns.")

      N = (len(row) - 3 - adcol)//8 #Number of planets
      
      output = np.array([row[0]])
      
      daysInUOT = 1.

      for i in range(N):
            lPkhiO = row[1 + 8*i : 7 + 8*i]
            if (typeOfCoordinates == 'Jacobi'):
                  mass = 1.
                  for j in range(i + 1):
                        mass = mass + row[7 + 8*j]
            elif (typeOfCoordinates == 'Heliocentric'):
                  mass = 1. + row[7 + 8*i]
            elif (typeOfCoordinates == 'JacobiWisdomHolman'):
                  massum = 1. # M* + m1 + ... + mi
                  for j in range(i + 1):
                        massum = massum + row[7 + 8*j]
                  mass = massum/(massum - row[7 + 8*i]) # M*(M* + m1 + ... + mi)/(M* + m1 + ... + m_{i-1})
            else:
                  raise Exception("typeOfCoordinates can be either 'Jacobi' or 'Heliocentric' or 'JacobiWisdomHolman'.")
            mass = mass*row[-2-adcol]
            aeiMoO = np.array(lPkhiO_to_aeiMoO(lPkhiO, mass, daysInUOT))
            output = np.concatenate((output, aeiMoO, np.array([row[7 + 8*i], row[8 + 8*i]])))

      output = np.concatenate((output, np.array([row[-2-adcol], row[-1-adcol]])))
      for i in range(adcol):
            output = np.concatenate((output, np.array([row[-adcol+i]])))
      return output
      
      
      
def sample2alkhqp_row(row, typeOfCoordinates, adcol):

      #Converts one row of the workshop into the same row in coordinates a, lbd, k, h, q, p
      #[a, lbd, k, h, q, p] = [semi-major axis, mean longitude, e*cos(varpi), e*sin(varpi), sin(i/2)*cos(Omega), sin(i/2)*sin(Omega)]
      #typeOfCoordinates is either 'Jacobi', 'JacobiWisdomHolman' or 'Heliocentric'
      #typeOfCoordinates must be the type of coordinates into which the sample is given !
      #This function cannot be used to change the type of coordinates ! Its purpose is only to convert to a, lbd, k, h, q, p
      #Use functions Jac2Hel_row or Hel2Jac_row to change the type of coordinates of the sample row
      #adcol is the number of additional columns in the sample (often 0). (e.g. adcol = 1 for Hadden's posteriors)

      n = len(row)
      if (n%8 != 3 + adcol):
            raise Exception("The sample line does not have 8*k + 3 + adcol columns.")

      N = (len(row) - 3 - adcol)//8 #Number of planets
      
      output = np.array([row[0]])
      
      daysInUOT = 1.

      for i in range(N):
            lPkhiO = row[1 + 8*i : 7 + 8*i]
            if (typeOfCoordinates == 'Jacobi'):
                  mass = 1.
                  for j in range(i + 1):
                        mass = mass + row[7 + 8*j]
            elif (typeOfCoordinates == 'Heliocentric'):
                  mass = 1. + row[7 + 8*i]
            elif (typeOfCoordinates == 'JacobiWisdomHolman'):
                  massum = 1. # M* + m1 + ... + mi
                  for j in range(i + 1):
                        massum = massum + row[7 + 8*j]
                  mass = massum/(massum - row[7 + 8*i]) # M*(M* + m1 + ... + mi)/(M* + m1 + ... + m_{i-1})
            else:
                  raise Exception("typeOfCoordinates can be either 'Jacobi' or 'Heliocentric' or 'JacobiWisdomHolman'.")
            mass = mass*row[-2-adcol]
            aeiMoO = lPkhiO_to_aeiMoO(lPkhiO, mass, daysInUOT)
            alkhqp = aeiMoO_to_alkhqp(aeiMoO, mass, daysInUOT)
            output = np.concatenate((output, alkhqp, np.array([row[7 + 8*i], row[8 + 8*i]])))

      output = np.concatenate((output, np.array([row[-2-adcol], row[-1-adcol]])))
      for i in range(adcol):
            output = np.concatenate((output, np.array([row[-adcol+i]])))
      return output
      

def cart2sample_row(row, typeOfCoordinates, adcol):

      #Inverse of function sample2cart_row

      n = len(row)
      if (n%8 != 3 + adcol):
            raise Exception("The sample line does not have 8*k + 3 + adcol columns.")

      N = (len(row) - 3 - adcol)//8 #Number of planets
      
      output = np.array([row[0]])
      
      daysInUOT = 1.

      for i in range(N):
            XYZvXvYvZ = row[1 + 8*i : 7 + 8*i]
            if (typeOfCoordinates == 'Jacobi'):
                  mass = 1.
                  for j in range(i + 1):
                        mass = mass + row[7 + 8*j]
            elif (typeOfCoordinates == 'Heliocentric'):
                  mass = 1. + row[7 + 8*i]
            elif (typeOfCoordinates == 'JacobiWisdomHolman'):
                  massum = 1. # M* + m1 + ... + mi
                  for j in range(i + 1):
                        massum = massum + row[7 + 8*j]
                  mass = massum/(massum - row[7 + 8*i]) # M*(M* + m1 + ... + mi)/(M* + m1 + ... + m_{i-1})
            else:
                  raise Exception("typeOfCoordinates can be either 'Jacobi' or 'Heliocentric' or 'JacobiWisdomHolman'.")
            mass = mass*row[-2-adcol]
            alkhqp = cart2ell(XYZvXvYvZ, mass)
            a      = alkhqp[0]
            lbd    = alkhqp[1]
            k      = alkhqp[2]
            h      = alkhqp[3]
            q      = alkhqp[4]
            p      = alkhqp[5]
            n      = np.sqrt(G*mass/a**3)
            P      = 2.*np.pi/n
            P      = P*daysInUOT
            lbd    = 180./np.pi*lbd
            Omega  = np.arctan2(p, q)
            sini2  = np.sqrt(p**2 + q**2)
            I      = 2.*np.arcsin(sini2)
            I      = 180./np.pi*I
            Omega  = 180./np.pi*Omega
            lPkhiO = np.array([lbd, P, k, h, I, Omega])
            output = np.concatenate((output, lPkhiO, np.array([row[7 + 8*i], row[8 + 8*i]])))

      output = np.concatenate((output, np.array([row[-2-adcol], row[-1-adcol]])))
      for i in range(adcol):
            output = np.concatenate((output, np.array([row[-adcol+i]])))
      return output


def Jac2Hel_rowCart(row, adcol):

      #To be called on the return value of sample2cart_row
      #Converts from Jacobi to Heliocentric
      
      n = len(row)
      if (n%8 != 3 + adcol):
            raise Exception("The sample line does not have 8*k + 3 + adcol columns.")

      N = (len(row) - 3 - adcol)//8 #Number of planets
      
      output = np.zeros((len(row)))
      for i in range(len(row)):
            output[i] = row[i]
      
      masses     = [1.]
      total_mass = [1.]
      
      for i in range(N):
            masses.append(row[7 + 8*i])
      for i in range(N):
            m = 1.
            for j in range(i + 1):
                  m = m + row[7 + 8*j]
            total_mass.append(m)
      masses     = np.array(masses)
      total_mass = np.array(total_mass)
      
      M_J = np.zeros((N + 1, N + 1)) #Conversion matrix from Barycentric to Jacobi
      M_H = np.zeros((N + 1, N + 1)) #Conversion matrix from Barycentric to Heliocentric

      for i in range(N + 1):
            M_H[i, i] = 1.
            M_J[i, i] = 1.

      for i in range(N):
            M_H[i + 1, 0] = -1.

      M_J[1,0] = -1.
      for i in range(2, N + 1):
            denom = total_mass[i - 1]
            for j in range(i):
                  num = masses[j]
                  M_J[i,j] = -num/denom
      
      M = M_H @ np.linalg.inv(M_J) #Conversion matrix from Jacobi to Heliocentric
      
      Xs  = np.zeros((N + 1))
      Ys  = np.zeros((N + 1))
      Zs  = np.zeros((N + 1))
      vXs = np.zeros((N + 1))
      vYs = np.zeros((N + 1))
      vZs = np.zeros((N + 1))
      
      for i in range(N):
            Xs [i + 1] = row[1 + 8*i]
            Ys [i + 1] = row[2 + 8*i]
            Zs [i + 1] = row[3 + 8*i]
            vXs[i + 1] = row[4 + 8*i]
            vYs[i + 1] = row[5 + 8*i]
            vZs[i + 1] = row[6 + 8*i]
      
      newX  = np.transpose(M @ np.transpose(Xs))
      newY  = np.transpose(M @ np.transpose(Ys))
      newZ  = np.transpose(M @ np.transpose(Zs))
      newvX = np.transpose(M @ np.transpose(vXs))
      newvY = np.transpose(M @ np.transpose(vYs))
      newvZ = np.transpose(M @ np.transpose(vZs))
      
      for i in range(N):
            output[1 + 8*i] = newX[i + 1]
            output[2 + 8*i] = newY[i + 1]
            output[3 + 8*i] = newZ[i + 1]
            output[4 + 8*i] = newvX[i + 1]
            output[5 + 8*i] = newvY[i + 1]
            output[6 + 8*i] = newvZ[i + 1]
            
      return output
            

def Hel2Jac_rowCart(row, adcol):

      #To be called on the return value of sample2cart_row
      #Converts from Heliocentric to Jacobi
      
      n = len(row)
      if (n%8 != 3 + adcol):
            raise Exception("The sample line does not have 8*k + 3 + adcol columns.")

      N = (len(row) - 3 - adcol)//8 #Number of planets
      
      output = np.zeros((len(row)))
      for i in range(len(row)):
            output[i] = row[i]
      
      masses     = [1.]
      total_mass = [1.]
      
      for i in range(N):
            masses.append(row[7 + 8*i])
      for i in range(N):
            m = 1.
            for j in range(i + 1):
                  m = m + row[7 + 8*j]
            total_mass.append(m)
      masses     = np.array(masses)
      total_mass = np.array(total_mass)
      
      M_J = np.zeros((N + 1, N + 1)) #Conversion matrix from Barycentric to Jacobi
      M_H = np.zeros((N + 1, N + 1)) #Conversion matrix from Barycentric to Heliocentric

      for i in range(N + 1):
            M_H[i, i] = 1.
            M_J[i, i] = 1.

      for i in range(N):
            M_H[i + 1, 0] = -1.

      M_J[1,0] = -1.
      for i in range(2, N + 1):
            denom = total_mass[i - 1]
            for j in range(i):
                  num = masses[j]
                  M_J[i,j] = -num/denom
      
      M = M_J @ np.linalg.inv(M_H) #Conversion matrix from Heliocentric to Jacobi
      
      Xs  = np.zeros((N + 1))
      Ys  = np.zeros((N + 1))
      Zs  = np.zeros((N + 1))
      vXs = np.zeros((N + 1))
      vYs = np.zeros((N + 1))
      vZs = np.zeros((N + 1))
      
      for i in range(N):
            Xs [i + 1] = row[1 + 8*i]
            Ys [i + 1] = row[2 + 8*i]
            Zs [i + 1] = row[3 + 8*i]
            vXs[i + 1] = row[4 + 8*i]
            vYs[i + 1] = row[5 + 8*i]
            vZs[i + 1] = row[6 + 8*i]
      
      newX  = np.transpose(M @ np.transpose(Xs))
      newY  = np.transpose(M @ np.transpose(Ys))
      newZ  = np.transpose(M @ np.transpose(Zs))
      newvX = np.transpose(M @ np.transpose(vXs))
      newvY = np.transpose(M @ np.transpose(vYs))
      newvZ = np.transpose(M @ np.transpose(vZs))
      
      for i in range(N):
            output[1 + 8*i] = newX[i + 1]
            output[2 + 8*i] = newY[i + 1]
            output[3 + 8*i] = newZ[i + 1]
            output[4 + 8*i] = newvX[i + 1]
            output[5 + 8*i] = newvY[i + 1]
            output[6 + 8*i] = newvZ[i + 1]
            
      return output
            

def Jac2Hel_row(row, adcol):

      #Converts one row of sample from Jacobi to Heliocentric
      
      rowCartJac       = sample2cart_row(row, 'Jacobi', adcol)
      rowCartHel       = Jac2Hel_rowCart(rowCartJac, adcol)
      rowlPkhiO        = cart2sample_row(rowCartHel, 'Heliocentric', adcol)
      return rowlPkhiO
      

def Hel2Jac_row(row, adcol):

      #Converts one row of sample from Heliocentric to Jacobi
      
      rowCartHel       = sample2cart_row(row, 'Heliocentric', adcol)
      rowCartJac       = Hel2Jac_rowCart(rowCartHel, adcol)
      rowlPkhiO        = cart2sample_row(rowCartJac, 'Jacobi', adcol)
      return rowlPkhiO
      
      
def JacWH2Hel_row(row, adcol):

      #Converts one row of sample from JacobiWisdomHolman to Heliocentric
      
      rowCartJac       = sample2cart_row(row, 'JacobiWisdomHolman', adcol)
      rowCartHel       = Jac2Hel_rowCart(rowCartJac, adcol)
      rowlPkhiO        = cart2sample_row(rowCartHel, 'Heliocentric', adcol)
      return rowlPkhiO
      

def Hel2JacWH_row(row, adcol):

      #Converts one row of sample from Heliocentric to JacobiWisdomHolman
      
      rowCartHel       = sample2cart_row(row, 'Heliocentric', adcol)
      rowCartJac       = Hel2Jac_rowCart(rowCartHel, adcol)
      rowlPkhiO        = cart2sample_row(rowCartJac, 'JacobiWisdomHolman', adcol)
      return rowlPkhiO


##########################################################################################################################
############################ End of auxiliary functions. Defining main functions #########################################
##########################################################################################################################


def Sample2cart(sample, typeOfCoordinates, adcol):

      n = sample.shape[1]
      output = np.copy(sample)
      print("Converting sample into Cartesian coordinates")
      print("Progress = ", 0., "%")
      K = n // 100
      for i in range(n):
            row    = sample[:,i]
            Newrow = sample2cart_row(row, typeOfCoordinates, adcol)
            output[:,i] = Newrow
            if ((i+1)%K == 0):
                  progress = 100.*(i + 1)/n
                  print("Progress = ", progress, "%")
      return output


def Cart2sample(sample, typeOfCoordinates, adcol):

      n = sample.shape[1]
      output = np.copy(sample)
      print("Converting sample from Cartesian to workshop format")
      print("Progress = ", 0., "%")
      K = n // 100
      for i in range(n):
            row    = sample[:,i]
            Newrow = cart2sample_row(row, typeOfCoordinates, adcol)
            output[:,i] = Newrow
            if ((i+1)%K == 0):
                  progress = 100.*(i + 1)/n
                  print("Progress = ", progress, "%")
      return output


def Sample2aeiMoO(sample, typeOfCoordinates, adcol):

      n = sample.shape[1]
      output = np.copy(sample)
      print("Converting sample into (a, e, i, M, o, O) elliptic elements")
      print("Progress = ", 0., "%")
      K = n // 100
      for i in range(n):
            row    = sample[:,i]
            Newrow = sample2aeiMoO_row(row, typeOfCoordinates, adcol)
            output[:,i] = Newrow
            if ((i+1)%K == 0):
                  progress = 100.*(i + 1)/n
                  print("Progress = ", progress, "%")
      return output


def Sample2alkhqp(sample, typeOfCoordinates, adcol):

      n = sample.shape[1]
      output = np.copy(sample)
      print("Converting sample into (a, lbd, k, h, q, p) elliptic elements")
      print("Progress = ", 0., "%")
      K = n // 100
      for i in range(n):
            row    = sample[:,i]
            Newrow = sample2alkhqp_row(row, typeOfCoordinates, adcol)
            output[:,i] = Newrow
            if ((i+1)%K == 0):
                  progress = 100.*(i + 1)/n
                  print("Progress = ", progress, "%")
      return output
      
      
def Jac2Hel(sample, adcol):

      n = sample.shape[1]
      output = np.copy(sample)
      print("Converting sample from Jacobi to Heliocentric coordinates")
      print("Progress = ", 0., "%")
      K = n // 100
      for i in range(n):
            row    = sample[:,i]
            Newrow = Jac2Hel_row(row, adcol)
            output[:,i] = Newrow
            if ((i+1)%K == 0):
                  progress = 100.*(i + 1)/n
                  print("Progress = ", progress, "%")
      return output
      
      
def Hel2Jac(sample, adcol):

      n = sample.shape[1]
      output = np.copy(sample)
      print("Converting sample from Heliocentric to Jacobi coordinates")
      print("Progress = ", 0., "%")
      K = n // 100
      for i in range(n):
            row    = sample[:,i]
            Newrow = Hel2Jac_row(row, adcol)
            output[:,i] = Newrow
            if ((i+1)%K == 0):
                  progress = 100.*(i + 1)/n
                  print("Progress = ", progress, "%")
      return output
      
def JacWH2Hel(sample, adcol):

      n = sample.shape[1]
      output = np.copy(sample)
      print("Converting sample from JacobiWisdomHolman to Heliocentric coordinates")
      print("Progress = ", 0., "%")
      K = n // 100
      for i in range(n):
            row    = sample[:,i]
            Newrow = JacWH2Hel_row(row, adcol)
            output[:,i] = Newrow
            if ((i+1)%K == 0):
                  progress = 100.*(i + 1)/n
                  print("Progress = ", progress, "%")
      return output
      
      
def Hel2JacWH(sample, adcol):

      n = sample.shape[1]
      output = np.copy(sample)
      print("Converting sample from Heliocentric to JacobiWisdomHolman coordinates")
      print("Progress = ", 0., "%")
      K = n // 100
      for i in range(n):
            row    = sample[:,i]
            Newrow = Hel2JacWH_row(row, adcol)
            output[:,i] = Newrow
            if ((i+1)%K == 0):
                  progress = 100.*(i + 1)/n
                  print("Progress = ", progress, "%")
      return output
      
      
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################

      
## Examples of use ##

'''
path = 'path_towards_sample.csv'
sample = np.loadtxt(path, dtype = np.float64, delimiter=',', unpack=True) #A sample in Heliocentric coordinates


S1 = Sample2cart(sample, 'Heliocentric', 0)         # Converts the sample into Heliocentric cartesian coordinates

S2 = Sample2cart(Hel2Jac(sample, 0), 'Jacobi', 0)   # Converts the sample into Jacobi cartesian coordinates

S3 = Sample2aeiMoO(sample, 'Heliocentric', 0)       # Converts the sample into Heliocentric elliptic elements (a, e, i, M, omega, Omega)

S4 = Sample2alkhqp(sample, 'Heliocentric', 0)       # Converts the sample into Heliocentric elliptic elements (a, l, k, h, q, p)

S5 = Sample2aeiMoO(Hel2Jac(sample, 0), 'Jacobi', 0) # Converts the sample into Jacobi elliptic elements (a, e, i, M, omega, Omega)

S6 = Sample2alkhqp(Hel2Jac(sample, 0), 'Jacobi', 0) # Converts the sample into Jacobi elliptic elements (a, l, k, h, q, p)

S7 = Sample2cart(Hel2JacWH(sample, 0), 'Jacobi', 0) # Converts the sample into Jacobi cartesian coordinates with the Wisdom-Holman convention for mu
'''


## Sanity checks :

# - Jac2Hel(Hel2Jac(sample, adcol), adcol) or Hel2Jac(Jac2Hel(sample, adcol), adcol) should leave the sample mostly untouched (within errors due to Kepler equation solving)
# - Sample2cart(Cart2sample(sample, 'Heliocentric', adcol), 'Heliocentric', adcol) and variations should leave the sample mostly untouched as well.



