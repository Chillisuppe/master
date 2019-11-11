function grid = ugrid_loadgrid(fname)

% Open netCDF file
ncid = netcdf.open(fname, 'NC_NOWRITE');

% Get global attributes
grid.title       = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'Title');
grid.conventions = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'Conventions');
grid.institution = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'Institution');
grid.references  = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'References');

% Get time step information
grid.time = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'time'));

% Get grid dimensions
[~, grid.nnodes] = netcdf.inqDim(ncid, netcdf.inqDimID(ncid, 'nMesh2_node'));
[~, grid.nelmts] = netcdf.inqDim(ncid, netcdf.inqDimID(ncid, 'nMesh2_face'));

% Get node coordinates and element-node connectivity
grid.nodecoor = [netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'Mesh2_node_x'), 'double'), ...
                 netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'Mesh2_node_y'), 'double')];
grid.face_nodes = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'Mesh2_face_nodes'))' + 1;

% Compute cell (triangle) centers
for ielmt = 1:grid.nelmts
  grid.elmtcenter(:,ielmt) = mean(grid.nodecoor(grid.face_nodes(ielmt,:),:), 1)';
end

% Close netCDF file
netcdf.close(ncid);
