"""
  functions for mapping between different coordinate systems in the triangle

  Nicole Beisiegel, April 2013
  Stefan Vater, October 2013
"""

import numpy as np

def xytors(x,y):
    """
    Map from (x,y) coordinates on equilateral triangle to (r,s) coordinates on
    right triangle
    """

    from math import sqrt

    L1 = (sqrt(3.0)*y + 1.0)/3.0
    L2 = (-3.0*x - sqrt(3.0)*y + 2.0)/6.0
    L3 = ( 3.0*x - sqrt(3.0)*y + 2.0)/6.0

    r = -L2 + L3 - L1
    s = -L2 - L3 + L1

    return r,s


def xytobaryequilat(x, y):
    """
    Map from (x,y) coordinates on equilateral triangle to barycentric
    coordinates
    """

    from math import sqrt

    f = np.finfo(float)
    L = np.zeros((3, x.size))

    xv = (-1.0, 1.0, 0.0)
    yv = (-1.0/sqrt(3.0), -1.0/sqrt(3.0), 2.0/sqrt(3.0))

    detT =  (yv[1]-yv[2])*(xv[0]-xv[2])+(xv[2]-xv[1])*(yv[0]-yv[2])
    L[0] = ((yv[1]-yv[2])*(x    -xv[2])+(xv[2]-xv[1])*(y    -yv[2])) / detT
    L[1] = ((yv[2]-yv[0])*(x    -xv[2])+(xv[0]-xv[2])*(y    -yv[2])) / detT
    L[2]   = 1.0 - L[0] - L[1]

    # set small values really to zero
    L[abs(L) < f.resolution] = 0.0

    return L.T


def rstobaryright(x, y):
    """
    Map from (r,s) coordinates on right triangle to barycentric
    coordinates
    """

    from math import sqrt

    L = np.zeros((3, x.size))

    xv = (-1.0,  1.0, -1.0)
    yv = (-1.0, -1.0,  1.0)

    detT =  (yv[1]-yv[2])*(xv[0]-xv[2])+(xv[2]-xv[1])*(yv[0]-yv[2])
    L[0] = ((yv[1]-yv[2])*(x    -xv[2])+(xv[2]-xv[1])*(y    -yv[2])) / detT
    L[1] = ((yv[2]-yv[0])*(x    -xv[2])+(xv[0]-xv[2])*(y    -yv[2])) / detT
    L[2]   = 1.0 - L[0] - L[1]
    return L.T


