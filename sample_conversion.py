                                                ###### Conversion of samples ######

# We provide here 6 conversion functions for the samples of the workshop.
# A row of sample is : Id or timestamp, lbd(°), Period(days), k, h, Inclination(°), Omega(°), m/M*, R/R*, ... , M* (Msun), R* (Rsun)
# Columns Id or timestamp, m/M*, R/R*, M* (Msun), R* (Rsun) are NOT modified by these functions.
# Only columns lbd(°), Period(days), k, h, Inclination(°), Omega(°) are modified for each planet.


# - Sample2cart(sample, typeOfCoordinates)   -->  Converts the sample to cartesian coordinates.
#                                                 The columns are now Id or timestamp, X, Y, Z, vX, vY, vZ, m/M*, R/R*, ... , M* (Msun), R* (Rsun)
#                                                 typeOfCoordinates is either 'Jacobi' or 'Heliocentric'

# - Sample2aeiMoO(sample, typeOfCoordinates) -->  Converts the sample to elliptic coordinates (a, e, i, M, omega, Omega) = (semi-major axis, eccentricity, inclination, 
#                                                 mean anomaly, argument of periapsis, longitude of ascending node).
#                                                 The columns are now Id or timestamp, a, e, i(rad), M(rad), omega(rad), Omega(rad), m/M*, R/R*, ... , M* (Msun), R* (Rsun)
#                                                 typeOfCoordinates is either 'Jacobi' or 'Heliocentric'

# - Sample2alkhqp(sample, typeOfCoordinates) -->  Converts the sample to elliptic coordinates (a, lbd, k, h, q, p) = (semi-major axis, mean longitude,
#                                                 e*cos(varpi), e*sin(varpi), sin(i/2)*cos(Omega), sin(i/2)*sin(Omega)).
#                                                 The columns are now Id or timestamp, a, lbd(rad), k, h, q, p, m/M*, R/R*, ... , M* (Msun), R* (Rsun)
#                                                 typeOfCoordinates is either 'Jacobi' or 'Heliocentric'

# - Cart2sample(sample, typeOfCoordinates)   -->  Inverse of Sample2cart.

# CAUTION : typeOfCoordinates indicates the type of coordinates the sample is given in. It cannot be used to change from Jacobi to Heliocentric or vice-versa.
#           For example, if your sample is in Jacobi coordinates, Sample2cart(sample, 'Jacobi') will convert it into Jacobi cartesian coordinates, whereas
#           Sample2cart(sample, 'Heliocentric') will have undefined behavior. Use functions Jac2Hel and Hel2Jac for conversion Jacobi <--> Heliocentric

# - Jac2Hel(sample) --> Converts the sample from Jacobi coordinates to Heliocentric coordinates

# - Hel2Jac(sample) --> Converts the sample from Heliocentric coordinates to Jacobi coordinates


# By default, these functions convert into an absolute system of units where :
# - The unit of length is the Sun radius
# - The unit of mass   is the Sun mass
# - The unit of time   is the day

# If it is more convenient for you, you can convert into a relative system of units where
# - The unit of length is the semi-major axis of the first planet of the current row
# - The unit of mass   is the stellar mass of the current row
# - The unit of time   is the orbital period of a massless particle with same semi-major axis that the first planet of the current row

# To do so, set this parameter to 0
convert2absolute = 1

# Conversion functions written by Jérémy COUTURIER

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
G = 4.*np.pi**2
if (convert2absolute):
      G = 0.034053


def ell2cart_true(aeinuoO, mass):

      #Returns the cartesian coordinates. This is an auxiliary function

      #aeinuoO = [a, e, i, nu, omega, Omega] = [semi-major axis, eccentricity, inclination, true longitude, argument of the periapsis, longitude of the ascending node]
      #i, nu, omega, Omega are in radians
      #mass is the central mass. e.g. mass = M_star + m_j for heliocentric or M_star + m_1 + ... + m_j for Jacobi.
      
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
      r = a*(1. - e*e)/(1. + e*cosnu)
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
      

def cart2ell(XYZvXvYvZ, mass):

      #Converts the cartesian coordinates into the elliptic elements [a, lambda, k, h, q, p] where k + ih = e*exp(i*varpi) and q + ip = sin(I/2)*exp(i*Omega)
      #mass is the central mass. e.g. mass = M_star + m_j for heliocentric or M_star + m_1 + ... + m_j for Jacobi.
      
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
      
      #mass is the central mass. e.g. mass = M_star + m_j for heliocentric or M_star + m_1 + ... + m_j for Jacobi.

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


def lPkhiO_to_aeiMoO(lPkhiO, mass, daysInUOT):

      #Converts from the workshop format to elliptic elements [semi-major axis, eccentricity, inclination, mean longitude, argument of periapsis, longitude of ascending node]
      #lPkhiO = [mean longitude (deg), Period (days), k, h, inclination(deg), longitude of the ascending node(deg)]
      #mass is the central mass. e.g. mass = M_star + m_j for heliocentric or M_star + m_1 + ... + m_j for Jacobi.

      
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

      #Converts from the workshop format to elliptic elements [semi-major axis, eccentricity, inclination, mean longitude, argument of periapsis, longitude of ascending node]
      #lPkhiO = [mean longitude (deg), Period (days), k, h, inclination(deg), longitude of the ascending node(deg)]
      #mass is the central mass. e.g. mass = M_star + m_j for heliocentric or M_star + m_1 + ... + m_j for Jacobi.

      
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
      #mass is the central mass. e.g. mass = M_star + m_j for heliocentric or M_star + m_1 + ... + m_j for Jacobi.

      aeiMoO  = lPkhiO_to_aeiMoO(lPkhiO, mass, daysInUOT)
      aeinuoO = get_true_anomaly(aeiMoO, mass)
      cart    = ell2cart_true(aeinuoO, mass)

      return cart;


def sample2cart_row(row, typeOfCoordinates):

      #Converts one row of the workshop into the same row in cartesian coordinates
      #typeOfCoordinates can be either 'Jacobi' or 'Heliocentric'
      #typeOfCoordinates must be the type of coordinates into which the sample is given !
      #This function cannot be used to change the type of coordinates ! Its purpose is only to convert to cartesian
      #Use functions Jac2Hel_row or Hel2Jac_row to change the type of coordinates of the sample row

      n = len(row)
      if (n%8 != 3):
            raise Exception("The sample line does not have 8*k + 3 columns.")

      N = (len(row) - 3)//8 #Number of planets
      
      output = np.array([row[0]])
      
      if (convert2absolute):
            daysInUOT = 1.
      else:
            daysInUOT = row[2]*np.sqrt(1. + row[7]) #Orbital period of inner planet is (1 + m1/m0)**(-1/2) since it is not massless

      for i in range(N):
            lPkhiO = row[1 + 8*i : 7 + 8*i]
            if (typeOfCoordinates == 'Jacobi'):
                  mass = 1.
                  for j in range(i + 1):
                        mass = mass + row[7 + 8*j]
            elif (typeOfCoordinates == 'Heliocentric'):
                  mass = 1. + row[7 + 8*i]
            else:
                  raise Exception("typeOfCoordinates can be either 'Jacobi' or 'Heliocentric'.")
            if (convert2absolute):
                  mass = mass*row[-2]
            cart   = np.array(lPkhiO_to_cart(lPkhiO, mass, daysInUOT))
            output = np.concatenate((output, cart, np.array([row[7 + 8*i], row[8 + 8*i]])))

      output = np.concatenate((output, np.array([row[-2], row[-1]])))
      return output


def sample2aeiMoO_row(row, typeOfCoordinates):

      #Converts one row of the workshop into the same row in coordinates a, e, i, M, o, O
      #[a, e, i, M, o, O] = [semi-major axis, eccentricity, inclination, mean anomaly, argument of periapsis, longitude of ascending node]
      #typeOfCoordinates can be either 'Jacobi' or 'Heliocentric'
      #typeOfCoordinates must be the type of coordinates into which the sample is given !
      #This function cannot be used to change the type of coordinates ! Its purpose is only to convert to a, e, i, M, o, O
      #Use functions Jac2Hel_row or Hel2Jac_row to change the type of coordinates of the sample row

      n = len(row)
      if (n%8 != 3):
            raise Exception("The sample line does not have 8*k + 3 columns.")

      N = (len(row) - 3)//8 #Number of planets
      
      output = np.array([row[0]])
      
      if (convert2absolute):
            daysInUOT = 1.
      else:
            daysInUOT = row[2]*np.sqrt(1. + row[7]) #Orbital period of inner planet is (1 + m1/m0)**(-1/2) since it is not massless

      for i in range(N):
            lPkhiO = row[1 + 8*i : 7 + 8*i]
            if (typeOfCoordinates == 'Jacobi'):
                  mass = 1.
                  for j in range(i + 1):
                        mass = mass + row[7 + 8*j]
            elif (typeOfCoordinates == 'Heliocentric'):
                  mass = 1. + row[7 + 8*i]
            else:
                  raise Exception("typeOfCoordinates can be either 'Jacobi' or 'Heliocentric'.")
            if (convert2absolute):
                  mass = mass*row[-2]
            aeiMoO = np.array(lPkhiO_to_aeiMoO(lPkhiO, mass, daysInUOT))
            output = np.concatenate((output, aeiMoO, np.array([row[7 + 8*i], row[8 + 8*i]])))

      output = np.concatenate((output, np.array([row[-2], row[-1]])))
      return output
      
      
      
def sample2alkhqp_row(row, typeOfCoordinates):

      #Converts one row of the workshop into the same row in coordinates a, lbd, k, h, q, p
      #[a, lbd, k, h, q, p] = [semi-major axis, mean longitude, e*cos(varpi), e*sin(varpi), sin(i/2)*cos(Omega), sin(i/2)*sin(Omega)]
      #typeOfCoordinates can be either 'Jacobi' or 'Heliocentric'
      #typeOfCoordinates must be the type of coordinates into which the sample is given !
      #This function cannot be used to change the type of coordinates ! Its purpose is only to convert to a, lbd, k, h, q, p
      #Use functions Jac2Hel_row or Hel2Jac_row to change the type of coordinates of the sample row

      n = len(row)
      if (n%8 != 3):
            raise Exception("The sample line does not have 8*k + 3 columns.")

      N = (len(row) - 3)//8 #Number of planets
      
      output = np.array([row[0]])
      
      
      if (convert2absolute):
            daysInUOT = 1.
      else:
            daysInUOT = row[2]*np.sqrt(1. + row[7]) #Orbital period of inner planet is (1 + m1/m0)**(-1/2) since it is not massless

      for i in range(N):
            lPkhiO = row[1 + 8*i : 7 + 8*i]
            if (typeOfCoordinates == 'Jacobi'):
                  mass = 1.
                  for j in range(i + 1):
                        mass = mass + row[7 + 8*j]
            elif (typeOfCoordinates == 'Heliocentric'):
                  mass = 1. + row[7 + 8*i]
            else:
                  raise Exception("typeOfCoordinates can be either 'Jacobi' or 'Heliocentric'.")
            if (convert2absolute):
                  mass = mass*row[-2]
            aeiMoO = lPkhiO_to_aeiMoO(lPkhiO, mass, daysInUOT)
            alkhqp = aeiMoO_to_alkhqp(aeiMoO, mass, daysInUOT)
            output = np.concatenate((output, alkhqp, np.array([row[7 + 8*i], row[8 + 8*i]])))

      output = np.concatenate((output, np.array([row[-2], row[-1]])))
      return output
      

def cart2sample_row(row, typeOfCoordinates):

      #Inverse of function sample2cart_row

      n = len(row)
      if (n%8 != 3):
            raise Exception("The sample line does not have 8*k + 3 columns.")

      N = (len(row) - 3)//8 #Number of planets
      
      output = np.array([row[0]])
      
      #Getting the period of the inner planet
      XYZvXvYvZ = row[1 : 7]
      mass = 1. + row[7]
      if (convert2absolute):
            mass = mass*row[-2]
      alkhqp = cart2ell(XYZvXvYvZ, mass)
      n      = np.sqrt(G*mass/alkhqp[0]**3)
      Pin    = 2.*np.pi/n
      if (convert2absolute):
            daysInUOT = 1.
      else:
            daysInUOT = Pin*np.sqrt(1. + row[7]) #Orbital period of inner planet is (1 + m1/m0)**(-1/2) since it is not massless

      for i in range(N):
            XYZvXvYvZ = row[1 + 8*i : 7 + 8*i]
            if (typeOfCoordinates == 'Jacobi'):
                  mass = 1.
                  for j in range(i + 1):
                        mass = mass + row[7 + 8*j]
            elif (typeOfCoordinates == 'Heliocentric'):
                  mass = 1. + row[7 + 8*i]
            else:
                  raise Exception("typeOfCoordinates can be either 'Jacobi' or 'Heliocentric'.")
            if (convert2absolute):
                  mass = mass*row[-2]
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

      output = np.concatenate((output, np.array([row[-2], row[-1]])))
      return output


def Jac2Hel_rowCart(row):

      #To be called on the return value of sample2cart_row
      #Converts from Jacobi to Heliocentric
      
      n = len(row)
      if (n%8 != 3):
            raise Exception("The sample line does not have 8*k + 3 columns.")

      N = (len(row) - 3)//8 #Number of planets
      
      output = np.zeros((len(row)))
      for i in range(len(row)):
            output[i] = row[i]
      
      masses = [1.]
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
            

def Hel2Jac_rowCart(row):

      #To be called on the return value of sample2cart_row
      #Converts from Heliocentric to Jacobi
      
      n = len(row)
      if (n%8 != 3):
            raise Exception("The sample line does not have 8*k + 3 columns.")

      N = (len(row) - 3)//8 #Number of planets
      
      output = np.zeros((len(row)))
      for i in range(len(row)):
            output[i] = row[i]
      
      masses = [1.]
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
            

def Jac2Hel_row(row):

      #Converts one row of sample from Jacobi to Heliocentric
      
      C                = convert2absolute
      convert2absolute = 1
      rowCartJac       = sample2cart_row(row, 'Jacobi')
      rowCartHel       = Jac2Hel_rowCart(rowCartJac)
      rowlPkhiO        = cart2sample_row(rowCartHel, 'Heliocentric')
      convert2absolute = C
      return rowlPkhiO
      

def Hel2Jac_row(row):

      #Converts one row of sample from Jacobi to Heliocentric
      
      C                = convert2absolute
      convert2absolute = 1
      rowCartHel       = sample2cart_row(row, 'Heliocentric')
      rowCartJac       = Hel2Jac_rowCart(rowCartHel)
      rowlPkhiO        = cart2sample_row(rowCartJac, 'Jacobi')
      convert2absolute = C
      return rowlPkhiO


##########################################################################################################################
############################ End of auxiliary functions. Defining main functions #########################################
##########################################################################################################################


def Sample2cart(sample, typeOfCoordinates):

      n = sample.shape[1]
      print("Converting sample into Cartesian coordinates")
      print("Progress = ", 0., "%")
      K = n // 100
      for i in range(n):
            row    = sample[:,i]
            Newrow = sample2cart_row(row, typeOfCoordinates)
            sample[:,i] = Newrow
            if ((i+1)%K == 0):
                  progress = 100.*(i + 1)/n
                  print("Progress = ", progress, "%")
      return sample


def Cart2sample(sample, typeOfCoordinates):

      n = sample.shape[1]
      print("Converting sample from Cartesian to workshop format")
      print("Progress = ", 0., "%")
      K = n // 100
      for i in range(n):
            row    = sample[:,i]
            Newrow = cart2sample_row(row, typeOfCoordinates)
            sample[:,i] = Newrow
            if ((i+1)%K == 0):
                  progress = 100.*(i + 1)/n
                  print("Progress = ", progress, "%")
      return sample


def Sample2aeiMoO(sample, typeOfCoordinates):

      n = sample.shape[1]
      print("Converting sample into (a, e, i, M, o, O) elliptic elements")
      print("Progress = ", 0., "%")
      K = n // 100
      for i in range(n):
            row    = sample[:,i]
            Newrow = sample2aeiMoO_row(row, typeOfCoordinates)
            sample[:,i] = Newrow
            if ((i+1)%K == 0):
                  progress = 100.*(i + 1)/n
                  print("Progress = ", progress, "%")
      return sample


def Sample2alkhqp(sample, typeOfCoordinates):

      n = sample.shape[1]
      print("Converting sample into (a, lbd, k, h, q, p) elliptic elements")
      print("Progress = ", 0., "%")
      K = n // 100
      for i in range(n):
            row    = sample[:,i]
            Newrow = sample2alkhqp_row(row, typeOfCoordinates)
            sample[:,i] = Newrow
            if ((i+1)%K == 0):
                  progress = 100.*(i + 1)/n
                  print("Progress = ", progress, "%")
      return sample
      
      
def Jac2Hel(sample):

      n = sample.shape[1]
      print("Converting sample from Jacobi to Heliocentric coordinates")
      print("Progress = ", 0., "%")
      K = n // 100
      for i in range(n):
            row    = sample[:,i]
            Newrow = Jac2Hel_row(row)
            sample[:,i] = Newrow
            if ((i+1)%K == 0):
                  progress = 100.*(i + 1)/n
                  print("Progress = ", progress, "%")
      return sample
      
      
def Hel2Jac(sample):

      n = sample.shape[1]
      print("Converting sample from Heliocentric to Jacobi coordinates")
      print("Progress = ", 0., "%")
      K = n // 100
      for i in range(n):
            row    = sample[:,i]
            Newrow = Hel2Jac_row(row)
            sample[:,i] = Newrow
            if ((i+1)%K == 0):
                  progress = 100.*(i + 1)/n
                  print("Progress = ", progress, "%")
      return sample
      
      
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################

      
## Examples of use ##

path = './Agol/Agol_TRAPPIST1_0_samples.csv' # This sample is given in Jacobi coordinates
sample = np.loadtxt(path, dtype = np.float64, delimiter=',', unpack=True)

S1 = Sample2cart(sample, 'Jacobi')                  # Converts the sample into Jacobi cartesian coordinates

S2 = Sample2cart(Jac2Hel(sample), 'Heliocentric')   # Converts the sample into Heliocentric cartesian coordinates

S3 = Sample2aeiMoO(sample, 'Jacobi')                # Converts the sample into Jacobi elliptic elements (a, e, i, M, omega, Omega)

S4 = Sample2alkhqp(sample, 'Jacobi')                # Converts the sample into Jacobi elliptic elements (a, l, k, h, q, p)

S5 = Sample2aeiMoO(Jac2Hel(sample), 'Heliocentric') # Converts the sample into Heliocentric elliptic elements (a, e, i, M, omega, Omega)

S6 = Sample2alkhqp(Jac2Hel(sample), 'Heliocentric') # Converts the sample into Heliocentric elliptic elements (a, l, k, h, q, p)


## Sanity checks :

# If convert2absolute is 1, then : 
# - Jac2Hel(Hel2Jac(sample)) and Hel2Jac(Jac2Hel(sample)) should leave the sample mostly untouched (within numerical errors)
# - Sample2cart(Cart2sample(sample, 'Jacobi'), 'Jacobi') and Sample2cart(Cart2sample(sample, 'Heliocentric'), 'Heliocentric') should leave the sample mostly untouched as well.


