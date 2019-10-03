# Author : Shahram Talei @ 2017 
# this code creates random Rs for exponential disk
# see : https://github.com/jbailinua/probability/blob/master/Transformation%20of%20Probability%20-%20Generating%20Distributions.ipynb
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import random

from scipy.interpolate import interp1d
# make sure these numbers are the same as Fz.c
R0=2.8
N = 200000
############
R_sample = np.arange(0,5*R0,0.2)
xi_sample = (1-(1+R_sample/R0)*np.exp(-(R_sample/R0)))
r_from_xi = interp1d(xi_sample, R_sample)
xi = random(N)*(max(xi_sample))
print "max xi:",max(xi_sample)
r = r_from_xi(xi)
print len(r)
file = open('R.bin', 'wb')

r.tofile(file)
file.close()

fig = plt.figure(1)
figsize = fig.get_size_inches()
figsize[0] *= 2
fig.set_size_inches(figsize)
ax1 = fig.add_subplot(121)
ax1.hist(xi)
ax1.set_xlabel('$\\xi$')
ax1.set_ylabel('$P(\\xi)$')
ax1.set_xlim(-0.2,1.2)
ax2 = fig.add_subplot(122)
ax2.hist(r)
ax2.set_xlabel('$a$')
ax2.set_ylabel('$P(a)$')
ax2.set_xlim(0,5*R0)

plt.figure(2)
plt.plot(xi,r,'b.')
plt.title('Fitting')
plt.xlabel('$\\xi$')
plt.ylabel('$R(\\xi)$')
plt.show()


