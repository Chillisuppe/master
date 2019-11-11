# author: Stefan Vater, 02/2017
#
# Note: two non-standard python packages are needed implicitly for these functions:
#       pyugrid and cell_tree2d

import copy
import numpy as np
import types

import interpolation as interp
import quadrature as quad
import trianglemappings as trmap

import pyugrid

def xytobary(x, y, xv, yv):
    """
    Map from (x,y) coordinates barycentric coordinates with vertex
    coordinates (xv, yv)
    """

    from math import sqrt

    f = np.finfo(float)
    #L = np.zeros((3, x.size))
    L = np.zeros((3,) + x.shape)

    detT =  (yv[1]-yv[2])*(xv[0]-xv[2])+(xv[2]-xv[1])*(yv[0]-yv[2])
    L[0] = ((yv[1]-yv[2])*(x    -xv[2])+(xv[2]-xv[1])*(y    -yv[2])) / detT
    L[1] = ((yv[2]-yv[0])*(x    -xv[2])+(xv[0]-xv[2])*(y    -yv[2])) / detT
    L[2]   = 1.0 - L[0] - L[1]

    # set small values really to zero
    L[abs(L) < f.resolution] = 0.0

    return L.T


def get_crosssection(x1, x2, ugrid, varkws, time=0):
  x = np.linspace(x1[0], x2[0], 1201)
  y = np.linspace(x1[1], x2[1], 1201)

  faces = ugrid.faces
  coo   = ugrid.nodes

  vpts = {}
  for kw in varkws:
    vpts[kw] = np.zeros(x.shape)

  ptstriang = ugrid.locate_faces(np.array([x, y]).T)
  triangs   = np.unique(ptstriang)
  discont = np.empty(0, dtype=int)

  for tr_i in range(len(triangs)):
    pts_i       = np.where(ptstriang == triangs[tr_i])[0]
    pts_bary    = xytobary(x[pts_i], y[pts_i], coo[faces[triangs[tr_i]],0],
                           coo[faces[triangs[tr_i]],1])
    for kw in varkws:
      v = ugrid.data[kw]
      vpts[kw][pts_i] = np.dot(pts_bary, v.data[time,:][faces[triangs[tr_i]]])
    discont     = np.append(discont, pts_i[0])

  discont = np.sort(discont)[::-1]

  for di in range(len(discont)-1):
    x    = np.insert(x   , discont[di], np.nan)
    y    = np.insert(y   , discont[di], np.nan)
    for kw in varkws:
      vpts[kw] = np.insert(vpts[kw], discont[di], np.nan)

  return x, y, vpts


def getcommonmesh_remesh(ugrid0, ugrid1, varkws):
  """
  interpolate two UGRID meshes (create with IO_netcdfplot_dg_remesh) and the
  given variables to common fine "base" mesh. Locally one mesh must
  be a submesh in the other.
  """

  # compute triangle midpoints on both grids
  midpointsf0 = np.mean(ugrid0.nodes[ugrid0.faces],axis=1)
  midpointsf1 = np.mean(ugrid1.nodes[ugrid1.faces],axis=1)

  # locate which triangles from on grid are in which triangles from the other
  f1mf0 = ugrid1.locate_faces(midpointsf0)
  f0mf1 = ugrid0.locate_faces(midpointsf1)

  # compute areas of triangles to decide, which triangles are common, and which are
  # sub-triangles in the other grid
  nodesf0 = ugrid0.nodes[ugrid0.faces]
  areaf0  = 0.5*np.abs((nodesf0[:,0,0]-nodesf0[:,2,0])*(nodesf0[:,1,1]-nodesf0[:,0,1]) -
                       (nodesf0[:,0,0]-nodesf0[:,1,0])*(nodesf0[:,2,1]-nodesf0[:,0,1]))
  nodesf1 = ugrid1.nodes[ugrid1.faces]
  areaf1  = 0.5*np.abs((nodesf1[:,0,0]-nodesf1[:,2,0])*(nodesf1[:,1,1]-nodesf1[:,0,1]) -
                       (nodesf1[:,0,0]-nodesf1[:,1,0])*(nodesf1[:,2,1]-nodesf1[:,0,1]))

  # compute common nodes
  idxf0same = np.isclose(areaf0, areaf1[f1mf0])
  idxn0same = ugrid0.faces[idxf0same].ravel()
  idxn1same = ugrid1.faces[idxf0same[f0mf1]].ravel()
  newnodes = ugrid0.nodes[idxn0same]

  # setup variables
  vars0    = {}
  vars1    = {}
  vars0new = {}
  vars1new = {}

  for kw in varkws:
    vars0[kw]    = ugrid0.data[kw].data[0,:]
    vars1[kw]    = ugrid1.data[kw].data[0,:]
    vars0new[kw] = vars0[kw][idxn0same]
    vars1new[kw] = vars1[kw][idxn1same]

  # walk through all triangle in first mesh, which are not other mesh
  for f0i in np.where(~idxf0same)[0]:
    face0i   = ugrid0.faces[f0i]
    nodesf0i = ugrid0.nodes[face0i]
    areaf0i  = areaf0[f0i]

    # find corresponding triangle other mesh and compare areas
    f1i      = f1mf0[f0i]
    face1i   = ugrid1.faces[f1i]
    areaf1i  = areaf1[f1i]
    if (areaf0i < areaf1i):
      # if area is smaller, add nodes to new mesh and linearly
      # interpolate variables to finer mesh
      nodesf1i = ugrid1.nodes[face1i]
      newnodes = np.vstack([newnodes, ugrid0.nodes[face0i]])
      nodesf0b = xytobary(nodesf0i[:,0], nodesf0i[:,1], nodesf1i[:,0], nodesf1i[:,1])
      for kw in varkws:
        vars0new[kw] = np.append(vars0new[kw], vars0[kw][face0i])
        vars1new[kw] = np.append(vars1new[kw], np.dot(nodesf0b, vars1[kw][face1i]))
    else:
      # if are is larger, add nodes from other mesh to new mesh and interpolate
      # variables from first mesh
      f1add = np.where(f0mf1 == f0i)[0]
      nodesf1 = ugrid1.nodes[ugrid1.faces[f1add]]
      n1add   = ugrid1.faces[f1add].ravel()
      newnodes = np.vstack([newnodes, ugrid1.nodes[n1add]])
      nodesf1b = np.swapaxes(xytobary(nodesf1[:,:,0], nodesf1[:,:,1], nodesf0i[:,0], nodesf0i[:,1]), 0, 1)
      for kw in varkws:
        vars0new[kw] = np.append(vars0new[kw], np.dot(nodesf1b, vars0[kw][face0i]))
        vars1new[kw] = np.append(vars1new[kw], vars1[kw][ugrid1.faces[f1add]])

  # create new fine common mesh
  newfaces = np.reshape(np.arange(newnodes.shape[0]),(newnodes.shape[0]/3,3))
  ugridnew = pyugrid.UGrid(nodes     = newnodes,
                           faces     = newfaces,
                           mesh_name = ugrid0.mesh_name)

  # copy variables to new fine mesh
  for kw in varkws:
    uvar0new = copy.deepcopy(ugrid0.data[kw])
    uvar0new.data = np.array([vars0new[kw]])
    uvar0new.name = kw+'_0'
    ugridnew.add_data(uvar0new)

    uvar1new = copy.deepcopy(ugrid1.data[kw])
    uvar1new.data = np.array([vars1new[kw]])
    uvar1new.name = kw+'_1'
    ugridnew.add_data(uvar1new)

  return ugridnew


class QuadUGRID:
  """
  2D quadrature on UGRID grids (using pyugrid)
  """

  def __init__(self, DGdegree, QuadExactness):

    self.N    = DGdegree
    self.NQ   = QuadExactness

    f = np.finfo(float)

    # obtain Lagrange points and map them from (x,y) on equilateral triangle to
    # barycentric coordinates / (r,s) on right triangle
    xt, yt = interp.Nodes2D(self.N)
    r, s   = trmap.xytors(xt, yt)
    self.dofcoo = trmap.xytobaryequilat(xt, yt)

    # compute Vandermonde matrix and its inverse, cuttoff values smaller than
    # machine resolution
    vander    = interp.Vandermonde2D(self.N, r, s)
    vander[abs(vander) < f.resolution] = 0.0

    vanderinv = np.linalg.inv(vander)
    vanderinv[abs(vanderinv) < f.resolution] = 0.0

    # compute quadrature weights and barycentric coordinates
    self.wei, gl = quad.readGaussTriangle(self.NQ) # weights and coordinates on right triangle
    self.npts    = self.wei.size                   # number of Gauss points
    # map to barycentric coordinates
    self.quadcoo = trmap.rstobaryright(gl[:,0], gl[:,1])

    # compute node values of nodal (Lagrange) basis function
    self.psi     = np.dot(interp.Vandermonde2D(self.N, gl[:,0], gl[:,1]), vanderinv)
    self.psi[abs(self.psi) < f.resolution] = 0.0

  def __call__(self, ugrid, arg, fcn=None, timeslice=0):
    """
    compute intregal by quadrature rule
    """
    coonod = ugrid.nodes
    faces  = ugrid.faces

    if type(arg) is types.StringType:
      argdata = ugrid.data[arg].data[timeslice]
      # check that this is a nodal variable!
    elif type(arg) is types.ListType:
      argdata = {}
      for i in range(len(arg)):
        argdata[arg[i]] = ugrid.data[arg[i]].data[timeslice]
      # check that this is a nodal variable!

    if type(fcn) is types.NoneType:
      fcn = lambda x: x

    quadf = 0.0
    for iface in range(faces.shape[0]):
      ctri = coonod[faces[iface]]
      area = 0.5*np.abs((ctri[0,0]-ctri[2,0])*(ctri[1,1]-ctri[0,1]) -
                        (ctri[0,0]-ctri[1,0])*(ctri[2,1]-ctri[0,1]))

      if type(arg) is types.FunctionType:
        cooquad = np.dot(self.quadcoo, ctri).T
        dataqcoo = arg(cooquad)
      elif type(arg) is types.StringType:
        dataqcoo = np.dot(self.psi, argdata[faces[iface]])
      elif type(arg) is types.ListType:
        dataqcoo = {}
        for i in range(len(arg)):
          dataqcoo[arg[i]] = np.dot(self.psi, argdata[arg[i]][faces[iface]])

      quadf = quadf + np.dot(self.wei, fcn(dataqcoo))/self.wei.sum()*area

    return quadf
