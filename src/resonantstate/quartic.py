import numpy as np
import cmath as cm
import mpmath as mpm
from mpmath import mp, fp, mpc, mpf, pi


def _cubic(a_3, a_2, a_1, a_0, is_mp):
      # Returns the real roots of a_3*X^3 + a_2*X^2 + a_1*X + a_0 = 0 using analytical expressions of Cardan's method
      if (a_3 == 0.):
            if (a_2 == 0.):
                  if (a_1 == 0.):
                        return []
                  return [-a_0/a_1]
            Delta = a_1**2. - 4.*a_2*a_0
            if (Delta >= 0.):
                  if not is_mp:
                        return [(-a_1 + np.sqrt(Delta))/(2.*a_2), (-a_1 - np.sqrt(Delta))/(2.*a_2)]
                  else:
                        return [(-a_1 + mpm.sqrt(Delta))/(2.*a_2), (-a_1 - mpm.sqrt(Delta))/(2.*a_2)]
            return []
      if (a_3 != 1.):
            return _cubic(a_3/a_3, a_2/a_3, a_1/a_3, a_0/a_3, is_mp)
            
      ### Equation is Y^3 + p*Y + q = 0 with Y = X + s ###
      p = a_1 - a_2**2./3.
      q = a_0 - a_1*a_2/3. + 2.*a_2**3./27.
      s = a_2/3.
      D = q**2./4. + p**3./27.
      if (D < 0.): # 3 real solutions
            if not is_mp:
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
            else:
                  u3 = -q/2. + mpm.sqrt(D)
                  v3 = -q/2. - mpm.sqrt(D)
                  [mod_u3, arg_u3] = mpm.polar(u3)
                  [mod_v3, arg_v3] = mpm.polar(v3)
                  u  = mod_u3**(mpf(1.)/mpf(3.))*mpm.exp(mpc(0.,1.)*arg_u3/3.)
                  v  = mod_v3**(mpf(1.)/mpf(3.))*mpm.exp(mpc(0.,1.)*arg_v3/3.)
                  j  = mpm.exp( 2.*mpc(0.,1.)*pi/3.)
                  jb = mpm.exp(-2.*mpc(0.,1.)*pi/3.)
                  S1 = mpm.re(u + v)
                  S2 = mpm.re(j*u + jb*v)
                  S3 = mpm.re(jb*u + j*v)
            sol = [S1 - s, S2 - s, S3 - s]
      else: # 1 real solution
            if not is_mp:
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
            else:
                  u3 = -q/2. + mpm.sqrt(D)
                  v3 = -q/2. - mpm.sqrt(D)
                  if (u3 < 0.):
                        u = -(-u3)**(mpf(1.)/mpf(3.))
                  else:
                        u  = u3**(mpf(1.)/mpf(3.))
                  if (v3 < 0.):
                        v = -(-v3)**(mpf(1.)/mpf(3.))
                  else:
                        v  = v3**(mpf(1.)/mpf(3.))
            sol = [u + v - s]
      return sol
      
def quartic(a_4, a_3, a_2, a_1, a_0, acceptable_error = 1.e-11, dps = 16):
      r"""
      Computes the real roots of the quartic :math:`a_4X^4+a_3X^3+a_2X^2+a_1X+a_0` using analytical expressions of
      Ferrari's method when :math:`a_4\neq0` and Cardan's method when :math:`a_4=0`
      
      Author : Jeremy Couturier. https://jeremycouturier.com
      
      Parameters
      ----------
      a_j: Floats or mp floats (mpf).
            The coefficients of the quartic polynomial
      acceptable_error: Positive float much smaller than 1
            Any returned root :math:`X` verifies :math:`|a_4X^4+a_3X^3+a_2X^2+a_1X+a_0|/(|X|+|a_0|+|a_1|+|a_2|+|a_3|+|a_4|) <` acceptable_error. Default is 1.e-11. No minimum.
      dps: Positive integer
            Number of decimal places of floating point operations. Default and minimum is 16. No maximum. 
            When 16, floating point operations do not use mpmath, unless the :math:`a_j` are mp floats, or unless it is required. 
            The user generally does not need to specify this parameter, since decreasing acceptable_error will increase dps accordingly.
            
      Returns
      -------
      l : List of floats or mp floats
            The list of real roots of the quartic
            
      Edge cases
      ----------
      Closeness to bi-quartic : If :math:`a_1 + a_3^3/8 - a_2a_3/2 \ll a_j`, calculations must be done with more decimal places. 
            This is handled automatically and does not require any action from the user. However in that case, mpmath will always be used. 
      Closeness of two roots : When two roots are close, or when a root is multiple, the analytical expressions of Ferrari's method are numerically hard to evaluate. 
            This is because it involves substracting very large and nearly identical quantities. The problem can be solved by decreasing acceptable_error or by increasing dps (but see below). 
      The precision on the roots is bounded by the precision on the coefficients : Sometimes, it is required to compute the :math:`a_j` with more decimal places before calling this function.
            In this case, decreasing acceptable_error or increasing dps is useless. The :math:`a_j` must be computed using mpmath and be given to this function as mp floats (mpf)
      Some real roots can be missed : This is directly related to the above problem.
            With approximate :math:`a_j`, the quartic can barely miss the x-axis, while it would have touched it or crossed it with precise enough :math:`a_j`
      """
      if (dps < 16):
            dps = 16
      if (acceptable_error < 1.e-15):
            current_prec = mp.prec
            mp.prec = max(int(-np.log2(acceptable_error)) + 4, current_prec)
      if (dps > 16):
            current_dps = mp.dps
            mp.dps = max(dps, current_dps)
      if (mp.dps > 16):
            if not isinstance(a_4, mpf):
                  a_4 = mpf(a_4)
            if not isinstance(a_3, mpf):
                  a_3 = mpf(a_3)
            if not isinstance(a_2, mpf):
                  a_2 = mpf(a_2)
            if not isinstance(a_1, mpf):
                  a_1 = mpf(a_1)
            if not isinstance(a_0, mpf):
                  a_0 = mpf(a_0)
      if (a_4 != 1. and a_4 != 0.):
            return quartic(a_4/a_4, a_3/a_4, a_2/a_4, a_1/a_4, a_0/a_4, acceptable_error = acceptable_error, dps = dps)
      if (a_4 == 0.):
            sol = _cubic(a_3, a_2, a_1, a_0, isinstance(a_0, mpf) or isinstance(a_1, mpf) or isinstance(a_2, mpf) or isinstance(a_3, mpf))
      else:
            ### Equation is Y^4 + p*Y^2 + q*Y + r = 0 with Y = X + s ###
            p = a_2 - 3./8.*a_3**2.
            q = a_1 + a_3**3./8. - .5*a_2*a_3
            r = a_0 + a_2*a_3**2./16. - a_1*a_3/4. - 3.*a_3**4./256.
            s = a_3/4.
            ### Taking care of bi-quartic case ###
            if (abs(q)/(abs(p)+abs(r)) < 10.**(-dps+1)):
                  Sol = _cubic(0., p/p, p, r, isinstance(p, mpf) or isinstance(r, mpf))
                  if (len(Sol) == 0):
                        return []
                  [S1, S2] = Sol
                  sol = []
                  # First pair
                  if (S1 >= -10.**(-dps+1)):
                        if (dps <= 16):
                              sol.append( np.sqrt(abs(S1)) - s)
                              sol.append(-np.sqrt(abs(S1)) - s)
                        else:
                              sol.append( mpm.sqrt(abs(S1)) - s)
                              sol.append(-mpm.sqrt(abs(S1)) - s)
                  if (S2 >= -10.**(-dps+1)):
                        if (dps <= 16):
                              sol.append( np.sqrt(abs(S2)) - s)
                              sol.append(-np.sqrt(abs(S2)) - s)
                        else:
                              sol.append( mpm.sqrt(abs(S2)) - s)
                              sol.append(-mpm.sqrt(abs(S2)) - s)
            elif (abs(q)/(abs(p)+abs(r)) < 10.**(-dps//2+4)): ##Too close to bi-quartic
                  return quartic(a_4, a_3, a_2, a_1, a_0, acceptable_error = acceptable_error, dps = 2*dps)
            else:
                  ### Getting a solution of the resolving cubic ###
                  Sol = _cubic(p/p, p, p**2./4. - r, -q**2./8., isinstance(p, mpf) or isinstance(r, mpf) or isinstance(q, mpf))
                  Sol.sort()
                  if (Sol[-1] < 0.): #There are no real solutions
                        return []
                  if (dps <= 16):
                        sqM = np.sqrt(Sol[-1])
                        ### Getting first pair of solutions ###
                        sol = []
                        D = q/(2.*np.sqrt(2.)*sqM)-(p + sqM**2.)/2.
                        if (D >= -10.**(-dps+1)):
                              S1 = -sqM/np.sqrt(2.) + np.sqrt(abs(D))
                              S2 = -sqM/np.sqrt(2.) - np.sqrt(abs(D))
                              sol.append(S1 - s)
                              sol.append(S2 - s)
                        ### Getting second pair of solutions ###
                        D = -q/(2.*np.sqrt(2.)*sqM)-(p + sqM**2.)/2.
                        if (D >= -10.**(-dps+1)):
                              S3 = sqM/np.sqrt(2.) + np.sqrt(abs(D))
                              S4 = sqM/np.sqrt(2.) - np.sqrt(abs(D))
                              sol.append(S3 - s)
                              sol.append(S4 - s)
                  else:
                        sqM = mpm.sqrt(Sol[-1])
                        ### Getting first pair of solutions ###
                        sol = []
                        D = q/(2.*mpm.sqrt(2.)*sqM)-(p + sqM**2.)/2.
                        if (D >= -10.**(-dps+1)):
                              S1 = -sqM/mpm.sqrt(2.) + mpm.sqrt(abs(D))
                              S2 = -sqM/mpm.sqrt(2.) - mpm.sqrt(abs(D))
                              sol.append(S1 - s)
                              sol.append(S2 - s)
                        ### Getting second pair of solutions ###
                        D = -q/(2.*mpm.sqrt(2.)*sqM)-(p + sqM**2.)/2.
                        if (D >= -10.**(-dps+1)):
                              S3 = sqM/mpm.sqrt(2.) + mpm.sqrt(abs(D))
                              S4 = sqM/mpm.sqrt(2.) - mpm.sqrt(abs(D))
                              sol.append(S3 - s)
                              sol.append(S4 - s)
      #Verifying the solution
      for X in sol:
            error = abs(a_4*X**4 + a_3*X**3 + a_2*X**2 + a_1*X + a_0)/(abs(a_0) + abs(a_1) + abs(a_2) + abs(a_3) + abs(a_4) + abs(X))
            if (error > acceptable_error):
                  return quartic(a_4, a_3, a_2, a_1, a_0, acceptable_error = acceptable_error, dps = 2*dps) #Increasing number of decimal places if error > acceptable_error for one of the roots
      if (dps > 16):
            mp.prec = max(int(-np.log2(acceptable_error)) + 3, 53) 
            if (acceptable_error >= 1.e-15):
                  for i in range(len(sol)):
                        sol[i] = float(sol[i])
      return sol



#Example. Finding the roots of P(X) = (X - 1)(X - 3.4)(X - 2)(X - 2.000000000001) with at least 40 decimal places

#mp.dps = 60
#A1 = 1.
#A2 = 3.+mpf(4.)/10.
#A3 = 2.
#A4 = 2.+mpf(10.)**(-12.)
#a4 = 1.
#a3 = -A4-A3-A2-A1
#a2 = A3*A4+A2*A4+A1*A4+A2*A3+A1*A3+A1*A2
#a1 = -(A2*A3*A4)-A1*A3*A4-A1*A2*A4-A1*A2*A3
#a0 = A1*A2*A3*A4
#sol = quartic(a4, a3, a2, a1, a0, acceptable_error = 1.e-60, dps = 16)
#print(sol)

