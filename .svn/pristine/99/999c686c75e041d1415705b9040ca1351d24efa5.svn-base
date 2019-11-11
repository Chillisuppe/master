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
min = input("Choose minimal level of refinement\n")
max = input("Choose maximal level of refinement\n")
a = range(min,max)
if float(max) < 10:
   t = np.genfromtxt(SFDIR+"/output-conservation "+str(max)+".txt")
else:
    t = np.genfromtxt(SFDIR+"/output-conservation"+str(max)+".txt")
m = t.size
n = m/10
t.reshape(n,10)
tmax = t[:,0].max()

mm = []
xx = []
mvec= []
hh = []

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
plt.title("Momentum $\phi u$")

ax22 = fig.add_subplot(233)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_1$-Norm")
plt.title("Momentum $\phi v$")

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
axx1= fig2.add_subplot(311)
plt.xlabel("Shortest edge")
plt.ylabel("$L_2$-Error")
plt.title("Manufactured Solution Convergence Order")

axx2= fig2.add_subplot(334)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_2$-Error")
plt.title("Relative $L_2$-Error for ssh")

axx3= fig2.add_subplot(335)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_2$-Error")
plt.title("Relative $L_2$-Error for u")

axx4= fig2.add_subplot(336)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_2$-Error")
plt.title("Relative $L_2$-Error for v")


axx5= fig2.add_subplot(337)
plt.xlabel("time [s]")
plt.ylabel("$L_{\inf}$-Error")
plt.title("Max Error for ssh")

axx6= fig2.add_subplot(338)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_{\inf}$-Error")
plt.title("Max Error for u")

axx7= fig2.add_subplot(339)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_{\inf}$-Error")
plt.title("Max Error for v")

plt.tight_layout()
plt.savefig("Test1.png")

fig3=plt.figure(3)
axxx1=fig3.add_subplot(211)
plt.xlim(0,tmax)
plt.xlabel("computational time [s]")
plt.ylabel("$L_2$-Error")
plt.title("Max Error for v")

axxx2=fig3.add_subplot(212)
plt.xlim(0,tmax)
plt.xlabel("time [s]")
plt.ylabel("$L_{\inf}$-Error")
plt.title("Max Error for v")

plt.tight_layout()
plt.savefig("Test.png")

fig4=plt.figure(4)
az1=fig4.add_subplot(111)
plt.xlabel("Shortest edge")
plt.ylabel("$L_2$-Error")
plt.title("Logplot for Convergence")

plt.tight_layout()
plt.savefig("Logplot.png")

for x in a:
    if x < 10:
        z  = np.genfromtxt(SFDIR+"/output-conservation "+ str(x)+".txt")
       # zz = np.genfromtxt(SFDIR+"/output-energy "+ str(x)+".txt")
        y  = np.genfromtxt(SFDIR+"/output-absolute "+ str(x)+".txt")
        yy = np.genfromtxt(SFDIR+"/output-time "+ str(x)+".txt")
    else:
        z  = np.genfromtxt(SFDIR+"/output-conservation"+ str(x)+".txt")
       # zz = np.genfromtxt(SFDIR+"/output-energy"+ str(x)+".txt")
        y  = np.genfromtxt(SFDIR+"/output-absolute"+ str(x)+".txt")
        yy = np.genfromtxt(SFDIR+"/output-time"+ str(x)+".txt")
    r = tmax/z[0,0]
    m = z.size
    n = m/10
    z.reshape(n,10)
    #zz.reshape(n,3)
    y.reshape(n,10)
    #yy.reshape(n,7)
    ax1.plot(z[:,0],z[:,7], label="N="+str(x))
    ax2.plot(z[:,0],z[:,8])
    ax22.plot(z[:,0],z[:,9])
    #ax3.plot(zz[:,0],zz[:,1])
    #ax4.plot(zz[:,0],zz[:,2])
    axx2.plot(y[:,0], np.true_divide(y[:,1],z[:,1]))
    axx3.plot(y[:,0], np.true_divide(y[:,2],z[:,2]))
    axx4.plot(y[:,0], np.true_divide(y[:,3],z[:,3]))
    axx5.plot(y[:,0],y[:,4])
    axx6.plot(y[:,0],y[:,5])
    axx7.plot(y[:,0],y[:,6])
    #axxx1.plot(yy[:,3]+yy[:,4]+yy[:,5]+yy[:,6], np.true_divide(y[:,3],z[:,3]))
    maxim = np.true_divide(y[1:r,1],z[1:r,1]).max()
    minim = 1/(math.pow(2,(0.5*(x-1))))
    mm.append(np.true_divide(y[1:r,1],z[1:r,1]).max())
    mvec.append(math.log(maxim,2))
    xx.append(x)
    hh.append/(math.log(minim,2))
    print mvec,hh

#mvec = math.log(mm)


ax1.legend(loc='upper center', bbox_to_anchor=(1.15,0.55))

axx1.plot(xx,mm)
#plt.show()

az1.plot(hh,mvec)
plt.show()




# plt.plot(x,y,linestyle =' ',marker ='+',label ="$x ^2$ ( Markers )")

