"""

  Plotting routines for StormFlash convergence tests.
  Nicole Beisiegel, August 2012

"""
import matplotlib.pyplot as plt
import numpy as np
#from numpy import *

# Adjust
SFDIR = "/home/nicole/Development/StormFlash2d/compile/linux_g64"
# For convergence tests
min = 9
max = 9
#a=range(min,max)
a=9
plt.ioff()
fig = plt.figure()

ax1 = fig.add_subplot(221) # Aufteilung Gebiet 2x2, 1 Quadrant
#plt.xlim(0,1)
#plt.ylim(-2,2)
plt.xlabel("time [s]")
plt.ylabel("$L_1$-Norm")
plt.title("Mass")

ax2 = fig.add_subplot(222)
#plt.xlim(0,1)
#plt.ylim(0,1)
plt.xlabel("time [s]")
plt.ylabel("$L_1$-Norm")
plt.title("Momentum")

ax3 = fig.add_subplot(223)
#plt.xlim(0,1)
#plt.ylim(0,1)
plt.xlabel("time [s]")
plt.ylabel("$L_2$-Norm")
plt.title("Kinetic Energy")

ax4 = fig.add_subplot(224)
#plt.xlim(0,1)
#plt.ylim(0,1)
plt.xlabel("time [s]")
plt.ylabel("$L_2$-Norm")
plt.title("Potential Energy")

for x in a:
    z  = np.genfromtxt(SFDIR+"/output-conservation "+ str(x)+".txt")
    zz = np.genfromtxt(SFDIR+"/output-energy "+ str(x)+".txt")
    m = z.size
    n = m/7
    z.reshape(n,7)
    ax1.plot(z[:,0],z[:,1])
    ax2.plot(z[:,0],z[:,2])
    ax3.plot(zz[:,0],zz[:,1])
    ax4.plot(zz[:,0],zz[:,2])

plt.tight_layout()

plt.savefig("Test.png")
plt.show()

# plt.plot(x,y,linestyle =' ',marker ='+',label ="$x ^2$ ( Markers )")

