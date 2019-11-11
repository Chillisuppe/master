"""
  program for generating signature files

  Version 0.1

  Nicole Beisiegel, February 2013
  Stefan Vater, October 2013

  Note: To start the main program in interactive mode, just write "main()"
"""

import numpy as np
import interpolation as interp
import quadrature as quad
import trianglemappings as trmap


def write_signature(fname, N, m, edofcoo, gquad, equad, Dpsi, rot, ref, prolo, restr):
    """
    writes the signature to file
    """

    #import tkMessageBox

    numfrmt  = '% .14e'
    rowfrmtg = ['%2u']+[numfrmt]*2
    rowfrmte = ['%2u']+[numfrmt]*3
    rowfrmtw = ['%2u']+[numfrmt]

    stormflash_input = open(fname, "w")
    stormflash_input.write("################################################################################\n")
    stormflash_input.write("#Signature for DG Method\n")
    stormflash_input.write("#Created with signature.py version 0.1\n")
    stormflash_input.write("#Authors: Nicole Beisiegel, Stefan Vater\n")
    stormflash_input.write("################################################################################\n")
    stormflash_input.write("SIGNATURE_DESCRIPTION\n"
                         + "DG Method" + "\n"
                         + "POLY_DEGREE\n"
                         + str(N) +"\n"
                         + "NODE_DOFS\n"
                         + str(0) +"\n"
                         + "EDGE_DOFS\n"
                         + str(0) +"\n"
                         + "ELMT_DOFS\n"
                         + str(m) +"\n"
                         + "EDGE_QUADPOINTS\n"
                         + str(gquad.npts) +"\n"
                         + "ELMT_QUADPOINTS\n"
                         + str(equad.npts) +"\n"
                         + "EDGE_PSINONZERO\n"
                         + str(gquad.psinonzero) +"\n")

    stormflash_input.write("ELMT_DOFCOORDINATES\n")
    tmpdata = np.column_stack((np.arange(1,m+1),edofcoo))
    np.savetxt(stormflash_input, tmpdata, fmt=rowfrmte)

    stormflash_input.write("EDGE_QUADCOORDINATES\n")
    tmpdata = np.column_stack((np.arange(1,gquad.npts+1),gquad.coo))
    np.savetxt(stormflash_input, tmpdata, fmt=rowfrmtg)

    stormflash_input.write("EDGE_QUADWEIGHTS\n")
    tmpdata = np.column_stack((np.arange(1,gquad.npts+1),gquad.wei))
    np.savetxt(stormflash_input, tmpdata, fmt=rowfrmtw)

    stormflash_input.write("EDGE_PSIINDEX\n")
    np.savetxt(stormflash_input, [gquad.psiidx], fmt=['%2u']*gquad.psinonzero)

    stormflash_input.write("EDGE_PSIQUADPOINTS\n")
    tmpdata = np.column_stack((np.arange(1,gquad.npts+1),gquad.psi))
    np.savetxt(stormflash_input, tmpdata, fmt=['%2u']+[numfrmt]*gquad.psinonzero)

    stormflash_input.write("ELMT_QUADCOORDINATES\n")
    tmpdata = np.column_stack((np.arange(1,equad.npts+1),equad.coo))
    np.savetxt(stormflash_input, tmpdata, fmt=rowfrmte)

    stormflash_input.write("ELMT_QUADWEIGHTS\n")
    tmpdata = np.column_stack((np.arange(1,equad.npts+1),equad.wei))
    np.savetxt(stormflash_input, tmpdata, fmt=rowfrmtw)

    stormflash_input.write("ELMT_PSIQUADPOINTS\n")
    tmpdata = np.column_stack((np.arange(1,equad.npts+1),equad.psi))
    np.savetxt(stormflash_input, tmpdata, fmt=['%2u']+[numfrmt]*m)

    stormflash_input.write("DPSIDXI_MATRIX\n")
    np.savetxt(stormflash_input, Dpsi.Dxi, fmt=numfrmt)

    stormflash_input.write("DPSIDETA_MATRIX\n")
    np.savetxt(stormflash_input, Dpsi.Deta, fmt=numfrmt)

    stormflash_input.write("!--- Note: This is experimental and has to be checked!\n")
    stormflash_input.write("DOF_ROTATIONMATRIX\n")
    np.savetxt(stormflash_input, rot, fmt='%2u')

    stormflash_input.write("!--- Note: This is experimental and has to be checked!\n")
    stormflash_input.write("DOF_REFLECTIONMATRIX\n")
    np.savetxt(stormflash_input, ref, fmt='%2u')

    stormflash_input.write("!--- Note: This is experimental and has to be checked!\n")
    stormflash_input.write("ADAPT_PROLONGATION\n")
    np.savetxt(stormflash_input, prolo, fmt=numfrmt)

    stormflash_input.write("!--- Note: you have to adjust some zeros!\n")
    stormflash_input.write("ADAPT_RESTRICTION\n")
    np.savetxt(stormflash_input, restr, fmt=numfrmt)

    stormflash_input.close()

    print("The output file has been written.")
    # tkMessageBox.showinfo("Success", "The output file has been written.")


def main():
    """
    main program
    """

    f = np.finfo(float)

    # obtain desired interpolation order from commandline
    while True:
        try:
            N = int(raw_input("Please choose the interpolation order: "))
            if N > 0:
                break
            else:
                raise ValueError
        except ValueError:
            print("Oops!  That was no valid number.  Try again...")

    # total number of nodes in element
    m = (N+1)*(N+2)/2

    # obtain number of Gauss points for edge quadrature
    while True:
        try:
            ngpts = raw_input("Please choose the number of Gauss points for edge quadrature (default: "+str(m)+"): ")
            if (ngpts == ""):
              ngpts = m
            else:
              ngpts = int(ngpts)
            if ngpts > 0:
                break
            else:
                raise ValueError
        except ValueError:
            print("Oops!  That was no valid number.  Try again...")

    # obtain exactness of elemental Gauss quadrature
    while True:
        try:
            equadN = raw_input("Please choose the exactness of element quadrature devided by 2 (default: "+str(N)+"): ")
            if (equadN == ""):
              equadN = N
            else:
              equadN = int(equadN)
            if equadN > 0:
                break
            else:
                raise ValueError
        except ValueError:
            print("Oops!  That was no valid number.  Try again...")

    fname = raw_input("Please choose the signature file name (default: Signature.ftf): ")
    if (fname == ""):
        fname = "Signature.ftf"

    # coordinates of the right triangle
    xv  = np.array([-1.0,  1.0, -1.0])
    yv  = np.array([-1.0, -1.0,  1.0])

    # obtain Lagrange points and map them from (x,y) on equilateral triangle to
    # barycentric coordinates / (r,s) on right triangle
    xt, yt  = interp.Nodes2D(N)
    edofcoo = trmap.xytobaryequilat(xt, yt)
    r, s    = trmap.xytors(xt, yt)

    # compute Vandermonde matrix and its inverse, cuttoff values smaller than
    # machine resolution
    vander    = interp.Vandermonde2D(N, r, s)
    vander[abs(vander) < f.resolution] = 0.0

    vanderinv = np.linalg.inv(vander)
    vanderinv[abs(vanderinv) < f.resolution] = 0.0

    # element Gauss points, weights and nodal basis function evaluations at quad points
    class equad():
        eN      = equadN
        wei, gl = quad.readGaussTriangle(eN)    # weights and coordinates on right triangle
        npts    = wei.size                      # number of Gauss points
        # map to barycentric coordinates
        coo     = trmap.rstobaryright(gl[:,0], gl[:,1])
        # compute node values of nodal (Lagrange) basis function
        psi     = np.dot(interp.Vandermonde2D(N, gl[:,0], gl[:,1]), vanderinv)
        psi[abs(psi) < f.resolution] = 0.0

    # edge Gauss points, weights and nodal basis function evaluations at quad points
    class gquad():
        npts = ngpts                            # number of Gauss points
        xlgl, wei = quad.JacobiGQ(0, 0, npts-1) # coordinates and weights on interval [-1,1]

        # map to barycentric coordinates on the triangle
        blgl      = np.zeros((npts,3))
        blgl[:,1] = -(xlgl - 1.0)/2.0
        blgl[:,2] = 1.0 - blgl[:,1]
        # extract barycentric coordinates on the reference edge
        coo  = blgl[:,1:]

        # map to (r,s) coordinates on the right trangle
        lgl = np.zeros((npts,2))
        lgl[:,0] = np.dot(blgl, xv)
        lgl[:,1] = np.dot(blgl, yv)
        # compute node values of nodal (Lagrange) basis function
        psi = np.dot(interp.Vandermonde2D(N, lgl[:,0], lgl[:,1]), vanderinv)

        # determine basis functions nonzero on reference edge
        psiidx     = np.where(abs(psi).max(axis=0) > f.resolution)
        psi        = psi[:,psiidx[0]]
        psiidx     = psiidx[0]+1
        psinonzero = psiidx.size

    # differentiation matrices
    class Dpsi():
        Dxi, Deta = interp.Dmatrices2D(N, r, s, vander)
        Dxi[ abs(Dxi)  < f.resolution] = 0.0
        Deta[abs(Deta) < f.resolution] = 0.0

    # rotation matrix
    rot = np.zeros((m,4), dtype=int)
    rot[:,0] = np.arange(1,m+1) # first column has node indices
    rot[:,1] = np.arange(1,m+1) # second column is always the identity

    # other columns are obtained by circular permutation of barycentric coordinates
    for iperm in range(2):
        rotcoo = np.roll(edofcoo, shift=iperm+1, axis=1)
        for irow in range(m):
            rotwhere = np.all(np.abs(edofcoo - np.array([rotcoo[irow,:],]*m)) < f.resolution*10, axis=1)
            rot[irow,2+iperm] = np.where(rotwhere)[0]+1

    # reflection matrix
    ref = np.zeros((m,3), dtype=int)
    ref[:,0] = np.arange(1,m+1) # first column has node indices
    ref[:,1] = np.arange(1,m+1) # second column is always the identity

    # third column is obtained from switching second and third barycentric coordinates
    refcoo = edofcoo[:,(0,2,1)]
    for irow in range(m):
        refwhere = np.all(np.abs(edofcoo - np.array([refcoo[irow,:],]*m)) < f.resolution*10, axis=1)
        ref[irow,2] = np.where(refwhere)[0]+1

    # prolongation matrix
    ptsprolo = np.copy(edofcoo)
    ptsprolo[:,2] = ptsprolo[:,2]*0.5
    ptsprolo[:,1] = ptsprolo[:,1]+ptsprolo[:,2]

    # map to (r,s) coordinates on the right trangle
    ptsrs = np.zeros((m,2))
    ptsrs[:,0] = np.dot(ptsprolo, xv)
    ptsrs[:,1] = np.dot(ptsprolo, yv)
    # compute node values of nodal (Lagrange) basis function
    prolo = np.dot(interp.Vandermonde2D(N, ptsrs[:,0], ptsrs[:,1]), vanderinv)
    prolo[abs(prolo) < f.resolution] = 0.0

    # restriction matrix (no algorithm yet)
    restr = np.zeros((m,m))

    # write signature to file
    write_signature(fname, N, m, edofcoo, gquad, equad, Dpsi, rot, ref, prolo, restr)

if __name__ == "__main__":
    main()
