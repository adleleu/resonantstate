The columns of continuedSeparatrix.txt are : delta, Xmin, Xmax, Xint, Xext, Xhyp
where Xint is the X value of the resonance center, Xext is the external circulation, Xmin and Xmax are the X value of the separatrix around the resonance center and
Xhyp is the hyperbolic fixed point.
When delta < 1, Xmin and Xmax are defined as the two X values of the unique level line whose area is 6*pi (area of the cardioid when delta = 1).
Xext and Xhyp are undefined when delta < 1 and the two corresponding columns are 0.
The file only spans -20 <= delta <= 50.
For delta < -20, Xint = -2/(3*delta) + O(delta^-4), Xmax = Xint + sqrt(6) and Xmin = Xint - sqrt(6)
For delta > 50,  Xint = sqrt(3*delta) + 1/(3*delta) + O(delta^(-5/2))
The Hamiltonian is H(X,Y) = 3/2*delta*(X^2 + Y^2) - 1/4*(X^2 + Y^2)^2 + 2X

The frequency on the internal branch goes to infinity when delta goes to either infinity of minus infinity.
It reaches a minimum of sqrt(3)*2^(2/3) at delta = 0 and is 1 at delta = 1 when the separatrices appear.
For delta << -1, nu = -3*delta + 24/(27*delta^2) + O(delta^-5)
For delta >>  1, nu = 2*(3*delta)^(1/4) + O(delta^(-5/4))

For delta << -1, the frequency on the continued separatrix is nu = 6 - 3*delta + 24/(27*delta^2) + O(delta^-5)
