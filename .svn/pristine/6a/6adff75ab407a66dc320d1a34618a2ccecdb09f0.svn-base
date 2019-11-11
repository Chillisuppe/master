function var = ugrid_loadvar(fname, varname)

% Open netCDF file
ncid = netcdf.open(fname, 'NC_NOWRITE');

% Get variable ID
varid = netcdf.inqVarID(ncid, varname);

% Get global attributes
var.standard_name     = netcdf.getAtt(ncid, varid, 'standard_name');
var.units             = netcdf.getAtt(ncid, varid, 'units');
var.mesh              = netcdf.getAtt(ncid, varid, 'mesh');
var.location          = netcdf.getAtt(ncid, varid, 'location');

% Get variable data
var.data = netcdf.getVar(ncid, varid, 'double');

% Close netCDF file
netcdf.close(ncid);
