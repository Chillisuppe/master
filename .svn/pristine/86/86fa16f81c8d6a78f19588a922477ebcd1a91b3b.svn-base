# This module computes the solution to the 1-dimensional shocktube riemann problem for given initial data
# for a wet bed acc. to Toro adapted from Cristobal Castro's code.
import os
import numpy as np
import matplotlib.pyplot as plt
import math 

#Set initial conditions - to be altered by user:
hl = 1.0
hr = 0.0
ul = 0.0
ur = 0.0
x0 = 25.0
xlen = 50

g = 9.81
tol = 10e-6
t = 2.0
n = 1000
x = []
y = []
yy= []
for i in range(n):
    x.append(-x0 +(xlen)*(float(i)-1.0)/float(n))

def sf(h):
    erg = math.sqrt(g*h) 
    return erg

def swguess(hl,ul,al,hr,ur,ar):
    # Minimum h
    minh = min(hl,hr)
    # Two rarefaction guess
    erg = (0.5*(al+ar)+0.25*(ul-ur))**2.0/g
    if (erg > minh):
        #Use two shock guess
        h0 = erg
        gl = math.sqrt(0.5*g * (h0+hl)/(h0*hl))
        gr = math.sqrt(0.5*g * (h0+hr)/(h0*hr))
        erg = (gl*hl + gr*hr + ul - ur)/(gl + gr)
    return erg

def swfunc(h, hk):
    a = sf(h)
    ak = sf(hk)
    if (h <= hk):
        erg = 2.0*(a-ak)
    else:
        erg = (h-hk)*math.sqrt(0.5*g*((h+hk)/(h*hk)))
    return erg

def swfunc2(h, hk):
    a=sf(h)
    if (h <= hk):
        erg = g/a
    else:
        gk    = math.sqrt(0.5*g*(h+hk)/(h*hk))
        erg = gk - g*(h-hk)/(4.0*h**2.0*gk)
    return erg


#Compute sound speed
al = sf(hl)
ar = sf(hr)
h0 = 0.0
u0 = 0.0
it = 0

if (hl > 0.0) & (hr > 0.0):
    Delta_u_crit = 2.0*(al+ar)	       
    Delta_u = ur-ul
    Delta_h = tol + 1.0
    
    if (Delta_u_crit <= Delta_u):		#Data produce vacuum
        hstar = 0.0
        ustar = 0.0
    else:
        # Two rarefaction guess
        h_l = swguess(hl,ul,al,hr,ur,ar)
        while ((Delta_h>tol) & (it<100)):
            f  = swfunc(h_l,hl) + swfunc(h_l,hr) + Delta_u
            df = swfunc2(h_l,hl) + swfunc2(h_l,hr)
            hstar = h_l - f/df
            Delta_h = 2.0*math.fabs(hstar-h_l)/(hstar+h_l)
            it = it + 1
            h_l = hstar
	    ustar = 0.5*(ul+ur)+0.5*(swfunc(hstar,hr)-swfunc(hstar,hl))
            astar = sf(hstar) #Sound speed in star region
else:
    hstar = 0.0
    ustar = 0.0

print hstar, ustar , ustar*hstar
 #erg is  hstar and  hstar*ustar 

#Sample the function
for k in range(n):
    xt = x[k]/t
    if ((hl>0.0) & (hr>0.0)):				       # There is not a dry/wet case but data may produce dry bed in the middle
        al = sf(hl)						       # Sound speed in left region
        ar = sf(hr)						       # Sound speed in right region
           
        Delta_u = ur-ul
        Delta_u_crit = 2.0*(al+ar)
    
        if (Delta_u >= Delta_u_crit):			       # Data produce vacuum. Dry bed in the middle
            S_headl = ul - al
            S_taill = ul + 2.0*al
            S_headr = ur + ar
            S_tailr = ur - 2.0*ar
            if (xt<S_headl):
                h0 = hl
                u0 = ul
            elif (xt<S_taill):
                h0 = (ul+2.0*al-xt)**2/(9.0*g)
                u0 = (ul+2.0*al+2.0*xt)/3.0
            elif (xt<S_tailr):
                h0 = 0.0
                u0 = 0.0
            elif (xt<S_headr):
                h0 = (-ur+2.0*ar+xt)**2/(9.0*g)
                u0 = (ur-2.0*ar+2.0*xt)/3.0
            else:
                h0 = hr
                u0 = ur
        else:
            astar = sf(hstar)                                  # Sound speed in star region
	# Compute wave speeds
            if (hstar<=hl):                            # Left rarefaction wave
                S_headl = ul - al
                S_taill = ustar - astar
            else:
                ql      = math.sqrt(0.5*((hstar+hl)*hstar/hl**2))
                S_l     = ul - al*ql
                
            if (hstar<=hr):                            # Right rarefaction wave
                S_headr = ur + ar
                S_tailr = ustar + astar
            else:
                qr      = math.sqrt(0.5*((hstar+hr)*hstar/hr**2))
                S_r     = ur + ar*qr

            if (xt<ustar):
                  if (hstar<=hl):			   # Left rarefaction wave
                       if (xt<S_headl):
                           h0 = hl
                           u0 = ul
                       elif (xt<=S_taill):
                           a0 = (ul+2.0*al-xt)/3.
                           h0 = a0**2.0/g
                           u0 = (ul+2.0*al+2.0*xt)/3.0
                       else:
                           h0 = hstar
                           u0 = ustar
                  else:				   # Left shock wave
                       if (xt<S_l):
                           h0 = hl
                           u0 = ul
                       else:
                           h0 = hstar
                           u0 = ustar
            else:
                if (hstar<=hr):			   # Right rarefaction wave
                    if (xt>S_headr):
                        h0 = hr
                        u0 = ur
                    elif (xt>=S_tailr):
                        a0 = (-ur+2.0*ar+xt)/3.0
                        h0 = a0**2/g
                        u0 = (ur-2.0*ar+2.0*xt)/3.0
                    else:
                        h0 = hstar
                        u0 = ustar
                else:				   # Right shock wave
                    if (xt>S_r):
                        h0 = hr
                        u0 = ur
                    else:
                        h0 = hstar
                        u0 = ustar
    elif ((hl <=0) &(hr>0)):			   # Dry bed region on the left
        ar = sf(hr)					   # Sound speed in right region
        S_headr = ur + ar
        S_tailr = ur - 2.0*ar			
        if (xt>S_headr):
            h0 = hr
            u0 = ur
        elif (xt>S_tailr):
            a0 = (-ur+2.0*ar+xt)/3.0
            h0 = a0**2/g
            u0 = (ur-2.0*ar+2.0*xt)/3.0
        else:
            h0 = 0.0
            u0 = 0.0
    elif ((hr<=0) & (hl>0)):			    # Dry bed region on the right
        al = sf(hl)					    # Sound speed in left region
        S_headl = ul - al
        S_taill = ul + 2.0*al
        if (xt<S_headl):
            h0 = hl
            u0 = ul
        elif (xt<S_taill):
            a0 = (ul+2.0*al-xt)/3.0
            h0 = a0**2.0/g
            u0 = (ul+2.0*al+2.0*xt)/3.0
        else:
            h0 = 0.0
            u0 = 0.0
    else:
        h0 = 0.0
        u0 = 0.0
    y.append(h0)
    yy.append(u0)

#Plots for time t
plt.ioff()
fig = plt.figure(1)

for i in range(n):
    x[i]=x[i]+x0

ax1 = fig.add_subplot(211) # Aufteilung Gebiet 2x2, 1 Quadrant
plt.xlim(0,xlen)
plt.xlabel("length [m]")
plt.ylabel("$\phi$")
plt.title("Sea Surface Height")
plt.grid()

ax2 = fig.add_subplot(212) # Aufteilung Gebiet 2x2, 1 Quadrant
plt.xlim(0,xlen)
plt.xlabel("length [m]")
plt.ylabel("u")
plt.title("Velocity")
plt.grid()

plt.tight_layout()
plt.savefig("Shocktest_at_"+str(t)+"_sec.png")

ax1.plot(x,y)
ax2.plot(x,yy)

plt.show()
  
