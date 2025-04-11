                                                ###### Conversion of samples ######
                                                ###### Running sanity checks ######


# Author : Jérémy COUTURIER


import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib
import numpy as np
import random
import sample_conversion as sc

## Sanity checks :

# Jac2Hel(Hel2Jac(sample, adcol), adcol) or Hel2Jac(Jac2Hel(sample, adcol), adcol) should leave the sample mostly untouched (within errors due to Kepler equation solving)
# Sample2cart(Cart2sample(sample, 'Heliocentric', adcol), 'Heliocentric', adcol) and variations should leave the sample mostly untouched as well.


## Running sanity checks
path = './Kepler-60_2_samples.csv' #A sample in Heliocentric coordinates with no additional column
sample = np.loadtxt(path, dtype = np.float64, delimiter=',', unpack=True)

S1 = sc.Jac2Hel(sc.Hel2Jac(sample, 0), 0)
S2 = sc.Cart2sample(sc.Sample2cart(sample, 'Heliocentric', 0), 'Heliocentric', 0)

print(np.max(np.absolute(sample - S1))) # Should print a small value
print(np.max(np.absolute(sample - S2))) # Should print a small value as well
