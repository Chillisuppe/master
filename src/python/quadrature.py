"""
  functions for quadrature in 1D and on the triangle

  Nicole Beisiegel, April 2013
  Stefan Vater, October 2013
"""

import os
import numpy as np

def JacobiGQ(alpha, beta, N):
    """
    Compute the N'th order Gauss quadrature points, x, and weights, w,
    associated with the Jacobi polynomial of type (alpha,beta) with
    alpha, beta > -1 (!= -0.5).

    Note: adapted from Hesthaven and Warburton, 2008
    """

    from math import sqrt, gamma

    f = np.finfo(float)
    x = np.zeros(N+1)
    w = np.zeros(N+1)

    if (N==0):
        x[0] = (alpha-beta) / (alpha+beta+2.0)
        w[0] = 2.0
        return x,w

    # form symmetric matrix from recurrence.
    h1 = 2.0*np.arange(N+1) + alpha + beta
    J  = np.diag(-0.5 * (alpha**2 - beta**2)*h1 / (h1+2.0))

    for i in range(N):
        J[i,i+1] = 2.0/(h1[i]+2.0) * sqrt((i+1.0) * (i+1.0+alpha+beta) * (i+1.0+alpha) * \
                                          (i+1.0+beta) / ((h1[i]+1.0) * (h1[i]+3.0)))

    if ((alpha+beta) < f.resolution):
        J[0,0] = 0.0

    J = J + J.T

    # compute quadrature by eigenvalue solve
    [x,V] = np.linalg.eig(J)
    for j in range(N+1):
        w[j] = V[0,j]**2 * 2**(alpha+beta+1) / (alpha+beta+1.0) * \
               gamma(alpha+1) * gamma(beta+1) / gamma(alpha+beta+1)

    # obtain ascending order of Gauss points
    ind = np.argsort(x)

    return x[ind], w[ind]

#print JacobiGQ(0,0,2)


def JacobiGL(alpha, beta, N):
    """
    Compute the N'th order Gauss Lobatto quadrature points, x,
    associated with the Jacobi polynomial of type (alpha,beta) with
    alpha, beta > -1 (!= -0.5).

    Note: adapted from Hesthaven and Warburton, 2008
    """

    x = np.zeros(N+1)
    x[0] = -1.0
    x[N] =  1.0

    if (N==1):
        return x

    [xint,w] = JacobiGQ(alpha+1, beta+1, N-2)
    x[1:N] = xint
    return x

#print JacobiGL(-1,-1,5)


def readGaussTriangle(N):
    """
    obtain Gauss points and weights from file

    @param  N   exactness of quadrature devided by 2 (in terms of polynomial degree)
    """

    with open(os.path.join(os.path.dirname(__file__),'../../data/gausspoints','gausspoints{0:02d}.dat'.format(N)), 'r') as f:
      while True:
        line = f.readline().rstrip()
        if not line: break
        if (line == 'NUM_GAUSSPOINTS'):
          no = int(f.readline().rstrip())
          omega = np.zeros(no)
          gl    = np.zeros((no,2))
        if (line == 'GAUSS_POINTS'):
          for i in range(no):
            line = f.readline().rstrip()
            gl[i,0], gl[i,1], omega[i] = [float(x) for x in line.split(',')]

    return omega, gl


def legendre_gauss(ngl):

    from math import cos, pi
    from interpolation import legendre_poly

    kmax = 20
    n    = ngl-1
    nh   = (n+1)/2
    xgl  = np.zeros(ngl)
    wgl  = np.zeros(ngl)

    for i in range(1,nh+1):
        x = cos((2.*i - 1.) / (2.*n + 1.) * pi)

        for k in range(1,kmax+1):
            (p0,p0_1,p0_2,p1,p1_1,p1_2,p2,p2_1,p2_2,p00,p00_1,p00_2) = legendre_poly(n,x)
            dx = -p00/p00_1
            x  = x + dx
            if (abs(dx) < 1.0e-20): break

        xgl[n+2-i-1] = x
        wgl[n+2-i-1] = 2.0 / ((1.0-x**2) * p00_1**2)

    # check for Zero
    if (n+1 != 2*nh):
        x = 0
        (p0,p0_1,p0_2,p1,p1_1,p1_2,p2,p2_1,p2_2,p00,p00_1,p00_2) = legendre_poly(n,x)
        xgl[nh] = x
        wgl[nh] = 2.0 / ((1.0-x**2) * p00_1**2)

    # rest is symmetry
    for i in range(1,nh+1):
        xgl[i-1] = -xgl[n+2-i-1]
        wgl[i-1] = +wgl[n+2-i-1]

    return xgl, wgl


def legendre_gauss_lobatto(ngl):

    from math import cos, pi
    from interpolation import legendre_poly

    kmax = 20
    xgl  = np.zeros(ngl)
    wgl  = np.zeros(ngl)
    n    = ngl-1
    nh   = (n+1)/2

    # First find half of the Roots
    for i in range(1,nh+1):
        x = cos((2.*i - 1.) / (2.*n + 1.) * pi)

        for k in range(1,kmax+1):
            # Construct Legendre Polynomial and Derivatives
            (p0,p0_1,p0_2,p1,p1_1,p1_2,p2,p2_1,p2_2,p00,p00_1,p00_2) = legendre_poly(n,x)
            # Get next Newton Iterative
            dx = -(1.0 - x**2) * p0_1 / (-2.*x*p0_1 + (1.0 - x**2)*p0_2)
            x  = x + dx

            if (abs(dx) < 1.0e-20): break

        xgl[n+2-i-1] = x
        wgl[n+2-i-1] = 2.0 / (float(n*(n+1)) * p0**2)

        # Check for Zero
        if (n+1 != 2*nh):
            x = 0
            (p0,p0_1,p0_2,p1,p1_1,p1_2,p2,p2_1,p2_2,p00,p00_1,p00_2) = legendre_poly(n,x,)
            xgl[nh] = x
            wgl[nh] = 2.0 / (float(n*(n+1)) * p0**2)

    # Rest of the roots
    for i in range(1,nh+1):
        xgl[i-1] = -xgl[n+2-i-1]
        wgl[i-1] = +wgl[n+2-i-1]

    return xgl, wgl

#print legendre_gauss(3)
#print legendre_gauss_lobatto(3)
