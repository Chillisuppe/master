"""

  Plotting routines for StormFlash convergence tests if analytical solutions
  are known
  Nicole Beisiegel, August 2012

"""
import matplotlib.pyplot as plt
import numpy as np
import math
import os

# Adjust
SFDIR = os.getcwd() #"/home/nicole/Development/StormFlash2d/compile/linux_g32"
# For convergence tests please adjust
eingabe = input("Geben Sie den Verfeinerungsgrad ein \n")
x= eingabe
if x < 10:
    t = np.genfromtxt(SFDIR+"/output-conservation "+str(x)+".txt")
else:
    t = np.genfromtxt(SFDIR+"/output-conservation"+str(x)+".txt") 
m = t.size
n = m/10
t.reshape(n,10)
tmax = t[:,0].max()

plt.ioff()
fig = plt.figure(1)

ax1 = fig.add_subplot(231) # Aufteilung Gebiet 2x2, 1 Quadrant
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_1$-Norm")
plt.title("Mass")

ax2 = fig.add_subplot(232)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_1$-Norm")
plt.title("Momentum 1")

ax22 = fig.add_subplot(233)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_1$-Norm")
plt.title("Momentum 2")

ax3 = fig.add_subplot(234)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_2$-Norm")
plt.title("Kinetic Energy")

ax4 = fig.add_subplot(235)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_2$-Norm")
plt.title("Potential Energy")

plt.tight_layout()
plt.savefig("Test_conservation.png")

fig2 = plt.figure(2)

axx2= fig2.add_subplot(231)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_2$-Error")
plt.title("Relative $L_2$-Error for ssh")

axx3= fig2.add_subplot(232)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_2$-Error")
plt.title("Relative $L_2$-Error for u")

axx4= fig2.add_subplot(233)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_2$-Error")
plt.title("Relative $L_2$-Error for v")


axx5= fig2.add_subplot(234)
plt.xlabel("time [s]")
plt.ylabel("$L_{\inf}$-Error")
plt.title("Max Error for ssh")

axx6= fig2.add_subplot(235)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_{\inf}$-Error")
plt.title("Max Error for u")

axx7= fig2.add_subplot(236)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_{\inf}$-Error")
plt.title("Max Error for v")

plt.tight_layout()
plt.savefig("Test.png")


if x < 10:
    z  = np.genfromtxt(SFDIR+"/output-conservation "+ str(x)+".txt")
    zz = np.genfromtxt(SFDIR+"/output-energy "+ str(x)+".txt")
    y  = np.genfromtxt(SFDIR+"/output-absolute "+ str(x)+".txt")
else:
    z  = np.genfromtxt(SFDIR+"/output-conservation"+ str(x)+".txt")
    zz = np.genfromtxt(SFDIR+"/output-energy"+ str(x)+".txt")
    y  = np.genfromtxt(SFDIR+"/output-absolute"+ str(x)+".txt")

r = tmax/z[0,0]
m = z.size
n = m/10
z.reshape(n,10)
zz.reshape(n,3)
y.reshape(n,10)
ax1.plot(z[:,0],z[:,7], label="N="+str(x))
ax2.plot(z[:,0],z[:,8])
ax22.plot(z[:,0],z[:,9])
ax3.plot(zz[:,0],zz[:,1])
ax4.plot(zz[:,0],zz[:,2])
axx2.plot(y[:,0], np.true_divide(y[:,1],z[:,1]))
axx3.plot(y[:,0], np.true_divide(y[:,2],z[:,2]))
axx4.plot(y[:,0], np.true_divide(y[:,3],z[:,3]))
axx5.plot(y[:,0],y[:,4])
axx6.plot(y[:,0],y[:,5])
axx7.plot(y[:,0],y[:,6])
   
ax1.legend(loc='upper center', bbox_to_anchor=(1.15,0.55))

plt.show()





