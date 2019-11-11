"""
  functions for interpolation in 1D and on the triangle

  Nicole Beisiegel, April 2013
  Stefan Vater, October 2013
"""

import numpy as np


def _rstoab(r, s):
    """
    Map from (r,s) coordinates on right triangle to (a,b) coordinates on
    square [-1,1]^2 for evaluation of Jacobi polynomials

    Note: adapted from Hesthaven and Warburton, 2008
    """

    s.dtype = np.float
    a       = -np.ones(s.size)
    mask    = (s != 1.0)
    a[mask] = 2.0*(1.0+r[mask])/(1.0-s[mask]) - 1.0

    return (a, s)


def JacobiP(x, alpha, beta, N):
    """
    Evaluate the Jacobi polynomial of type (alpha, beta), alpha, beta > -1,
    (alpha+beta) != -1, at points x for order N and returns P(x). The
    polynomials are normalized to be orthonormal.

    Note: adapted from Hesthaven and Warburton, 2008
    """

    from math import sqrt, gamma

    PL = np.zeros((N+1, np.size(x)))

    # initial values P_0(x) and P_1(x)
    gamma0  = 2**(alpha+beta+1) / (alpha+beta+1.0) * \
              gamma(alpha+1) * gamma(beta+1) / gamma(alpha+beta+1)
    PL[0,:] = 1.0 / sqrt(gamma0)

    if(N == 0):
        return PL[0]

    gamma1  = (alpha+1.0) * (beta+1.0) / (alpha+beta+3.0) * gamma0
    PL[1,:] = ((alpha+beta+2.0)*x/2.0 + (alpha-beta)/2.0) / sqrt(gamma1)

    if(N == 1):
        return PL[1]

    # repeat value in recurrence
    aold = 2.0 / (2.0+alpha+beta) * sqrt((alpha+1.0)*(beta+1.0) / (alpha+beta+3.0))

    # forward recurrence using the symmetry of the recurrence
    for i in range(1,N):
        h1   = 2.0*i + alpha + beta
        anew = 2.0/(h1+2.0) * sqrt((i+1.0) * (i+1.0+alpha+beta) * (i+1.0+alpha) * \
                                   (i+1.0+beta) / ((h1+1.0)*(h1 + 3.0)))
        bnew = - (alpha**2 - beta**2) / (h1*(h1+2.0))
        PL[i+1] = 1.0/anew * (-aold*PL[i-1] + (x-bnew)*PL[i])
        aold = anew

    return PL[N]


def Vandermonde1D(N,r):
    """
    Initialize the 1D Vandermonde Matrix, V_{ij} = phi_j(r_i)

    Note: adapted from Hesthaven and Warburton, 2008
    """

    V1D = np.zeros((np.size(r), N+1))

    for j in range(N+1):
        V1D[:,j]= JacobiP(r,0,0,j)

    return V1D


def GradJacobiP(r, alpha, beta, N):
    """
    Evaluate the derivative of the Jacobi polynomial of type (alpha,beta) with
    alpha, beta > -1, at points r for order N and returns dP[1:length(r))]

    Note: adapted from Hesthaven and Warburton, 2008
    """

    from math import sqrt

    if(N == 0):
        return np.zeros(np.size(r))
    else:
        return sqrt(N * (N+alpha+beta+1.0)) * JacobiP(r, alpha+1, beta+1, N-1)


def Simplex2DP(r, s, i, j):
    """
    Evaluate 2D orthonormal polynomial of order (i,j) on simplex at
    points (r,s).

    Note: adapted from Hesthaven and Warburton, 2008
    """

    from math import sqrt

    (a,b) = _rstoab(r,s)
    h1 = JacobiP(a,0,0,i)
    h2 = JacobiP(b,2*i+1,0,j)
    P  = sqrt(2.0) * h1 * h2 * (1-b)**i

    return P


def Warpfactor(N, rout):
    """
    Compute scaled warp function at order N based on rout interpolation nodes

    Note: adapted from Hesthaven and Warburton, 2008
    """

    from quadrature import JacobiGL

    # compute LGL and equidistant node distribution
    LGLr = JacobiGL(0,0,N)
    req  = np.linspace(-1,1,N+1)

    # compute V based on req
    Veq = Vandermonde1D(N, req)

    # evaluate Lagrange polynomial at rout
    Nr = np.size(rout)
    Pmat = np.zeros((N+1,Nr))
    for i in range(N+1):
        Pmat[i] = JacobiP(rout,0,0,i)

    Lmat = np.linalg.solve(Veq.T, Pmat)

    # Compute warp factor
    warp = np.dot(Lmat.T, LGLr-req)

    # Scale factor
    zerof = (np.absolute(rout)<1.0-1.0E-10)
    sf    = 1.0 - (zerof*rout)**2
    warp  = warp/sf + warp*(zerof-1.0)

    return warp

#print Warpfactor(3,(-1,0,1))


def Nodes2D(N):
    """
    Compute (x,y) Legendre-Gauss-Lobatto type nodes on equilateral
    triangle for polynomial interpolation of order N

    Note: adapted from Hesthaven and Warburton, 2008
    """

    from math import sqrt, cos, sin, pi

    alpopt = ( 0.0000, 0.0000, 0.0000, 1.4152, 0.1001, 0.2751, 0.9800, 1.0999, \
               1.2832, 1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258)

    # set optimized parameter, alpha, depending on order
    if (N < 16):
        alpha = alpopt[N]
    else:
        alpha = 5.0 / 3.0

    # total number of nodes
    Np=(N+1)*(N+2)/2

    # create equidistributed nodes on equilateral triangle
    L1 = np.zeros(Np)
    L3 = np.zeros(Np)
    sk = 2
    for k in range(N):
        for m in range(N+1-k):
            if(m!=N or k!=0):
              L1[sk] = float(k)/N
              L3[sk] = float(m)/N
              sk = sk+1

    L1[2] = 1.0
    L3[1] = 1.0
    L2    = 1.0 - L1 - L3

    x = -L2 + L3
    y = (-L2 - L3 + 2*L1) / sqrt(3.0)

    # compute blending function at each node for each edge
    blend1 = 4*L2*L3
    blend2 = 4*L1*L3
    blend3 = 4*L1*L2

    # amount of warp for each node, for each edge
    warpf1 = Warpfactor(N, L3-L2)
    warpf2 = Warpfactor(N, L1-L3)
    warpf3 = Warpfactor(N, L2-L1)

    # combine blend & warp
    warp1 = blend1 * warpf1 * (1 + (alpha*L1)**2)
    warp2 = blend2 * warpf2 * (1 + (alpha*L2)**2)
    warp3 = blend3 * warpf3 * (1 + (alpha*L3)**2)

    # accumulate deformations associated with each edge
    x = x + 1.0*warp1 + cos(2.0*pi/3.0)*warp2 + cos(4.0*pi/3.0)*warp3
    y = y + 0.0*warp1 + sin(2.0*pi/3.0)*warp2 + sin(4.0*pi/3.0)*warp3

    return x,y


def Vandermonde2D(N, r, s):
    """
    Initialize the 2D Vandermonde matrix, V_{ij} = phi_j(r_i, s_i)

    Note: adapted from Hesthaven and Warburton, 2008
    """

    m  = (N+1)*(N+2)/2
    V  = np.zeros((r.size,m))

    # build the Vandermonde matrix
    sk = 0
    for i in range(N+1):
        for j in range(N+1-i):
            V[:,sk] = Simplex2DP(r, s, i, j)
            sk = sk+1
    return V


def GradSimplex2DP(a, b, id, jd):
    """
    Return the derivatives of the modal basis (id,jd) on the 2D simplex at (a,b).

    Note: adapted from Hesthaven and Warburton, 2008
    """

    fa  = JacobiP(a, 0, 0, id)
    dfa = GradJacobiP(a, 0, 0, id)
    gb  = JacobiP(b, 2*id+1.0, 0, jd)
    dgb = GradJacobiP(b, 2*id+1.0, 0, jd)

    # r-derivative
    # d/dr = da/dr d/da + db/dr d/db = (2/(1-s)) d/da = (2/(1-b)) d/da
    dmodedr = dfa * gb
    if(id>0):
        dmodedr = dmodedr * ((0.5*(1.0-b))**(id-1))

    # s-derivative
    # d/ds = ((1+a)/2)/((1-b)/2) d/da + d/db
    dmodeds = dfa * (gb*(0.5*(1.0+a)))
    if(id>0):
        dmodeds = dmodeds * ((0.5*(1.0-b))**(id-1))

    tmp = dgb*((0.5*(1.0-b))**id)
    if(id>0):
        tmp = tmp - 0.5*id*gb * ((0.5*(1.0-b))**(id-1))

    dmodeds = dmodeds + fa*tmp

    # Normalize
    dmodedr = 2**(id+0.5) * dmodedr
    dmodeds = 2**(id+0.5) * dmodeds

    return dmodedr, dmodeds


def GradVandermonde2D(N, r, s):
    """
    Initialize the gradient of the modal basis (i,j) at (r,s) at order N

    Note: adapted from Hesthaven and Warburton, 2008
    """

    m    = (N+1)*(N+2)/2
    V2Dr = np.zeros((r.size, m))
    V2Ds = np.zeros((r.size, m))

    # find tensor-product coordinates
    [a,b] = _rstoab(r,s)

    # initialize matrices
    sk = 0
    for i in range(N+1):
        for j in range(N+1-i):
            V2Dr[:,sk], V2Ds[:,sk] = GradSimplex2DP(a,b,i,j)
            sk = sk+1

    return V2Dr, V2Ds


def Dmatrices2D(N, r, s, V):
    """
    Initialize the (r,s) differentiation matrices on the simplex, evaluated
    at (r,s) at order N (using the Vandermonde matrix V)

    Note: adapted from Hesthaven and Warburton, 2008
    """

    Vr, Vs = GradVandermonde2D(N, r, s)
    Dr     = np.linalg.solve(V.T, Vr.T)
    Ds     = np.linalg.solve(V.T, Vs.T)

    return Dr.T, Ds.T


def legendre_poly(n,x):

    p1   = 0
    p1_1 = 0
    p1_2 = 0
    p0   = 1
    p0_1 = 0
    p0_2 = 0

    # The nth order polynomials
    for j in range(1,n+1):
        p2   = p1
        p2_1 = p1_1
        p2_2 = p1_2
        p1   = p0
        p1_1 = p0_1
        p1_2 = p0_2
        a = (2.0*j-1.0) / float(j)
        b = (j-1.0) / float(j)
        p0   = a*x*p1 - b*p2
        p0_1 = a*( p1 + x*p1_1 ) - b*p2_1
        p0_2 = a*( 2.0*p1_1 + x*p1_2 ) - b*p2_2
        a = (2.0*j+1.0) / (j+1.0)
        b = float(j) / (j+1.0)
        p00   = a*x*p0 - b*p1
        p00_1 = a*( p0 + x*p0_1 ) - b*p1_1
        p00_2 = a*( 2.0*p0_1 + x*p0_2 ) - b*p1_2

    return p0, p0_1, p0_2, p1, p1_1, p1_2, p2, p2_1, p2_2, p00, p00_1,p00_2


def printtriangle(n):
    """
    print interpolation points with outline of equilateral and
    computational right triangle
    """

    import matplotlib.pyplot as plt
    from math import sqrt
    from trianglemappings import xytors

    fig = plt.figure(1)

    # nodes of equilateral triangle
    xt = np.array([-1.0, 1.0, 0.0, -1.0])
    yt = np.array([-1.0, -1.0, 2.0, -1.0])/sqrt(3.0)
    (x,y) = Nodes2D(n)

    # plot nodes on equilateral triangle
    ax1 = fig.add_subplot(211)
    ax1.set_xlim(-1.1, 1.1)
    ax1.set_ylim(-1.1/sqrt(3.0), 2.1/sqrt(3.0))
    ax1.set_xlabel("")
    ax1.set_aspect('equal')
    #ax1.tight_layout()
    ax1.plot(xt, yt)
    ax1.plot(x , y, '*')

    # nodes of computational right triangle
    xt2 = np.array([-1.0,  1.0, -1.0, -1.0])
    yt2 = np.array([-1.0, -1.0,  1.0, -1.0])
    (x2,y2) = xytors(x,y)

    ax2 = fig.add_subplot(212)
    ax2.set_xlim(-1.1,1.1)
    ax2.set_ylim(-1.1,1.1)
    ax2.set_aspect('equal')
    ax2.plot(xt2, yt2)
    ax2.plot(x2 , y2, '*')

    plt.show()
