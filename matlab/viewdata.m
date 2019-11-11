clear all

%% define data to load
basepath = '~/rechnen/amatos/StormFlash2d/compile/linux_g64';
datapath = 'vortexStrong08';
% datapath = '';
fileprefix = 'StormFlash2d_';

%% set min/max for plotting
x_lim = [0 4];
y_lim = [0 2];
minhgt = 0.9;
maxhgt = 1.0;
minu   = -0.5;
maxu   =  0.5;

ncfiles = dir(fullfile(basepath, datapath, [fileprefix '*.nc']));

figure(1)
for istep = 1:length(ncfiles)

  %% load data
  grfname = fullfile(basepath, datapath, ncfiles(istep).name);
  gr      = ugrid_loadgrid(grfname);

  hgt = ugrid_loadvar(grfname, 'Mesh2_height');
  bat = ugrid_loadvar(grfname, 'Mesh2_bathy');
  v_x = ugrid_loadvar(grfname, 'Mesh2_v_x');
  v_y = ugrid_loadvar(grfname, 'Mesh2_v_y');
  lev = ugrid_loadvar(grfname, 'Mesh2_level');

  %% plot data
  clf
  
  dummy  = -9999*ones([gr.nnodes, 1]);
  trdata = TriRep(double(gr.face_nodes), [gr.nodecoor, dummy]);

%% plot height
  subplot(321)
  trisurf(trdata, 'CData', hgt.data);
  view(0,90);
  % shading interp;
  colormap(jet);
  axis equal
  set(gca, 'XLim', x_lim, 'YLim', y_lim)
  caxis([minhgt maxhgt])
  colorbar;
  title('height')
  text(0, 1.1, ['time = ' num2str(gr.time)], 'Units', 'normalized')

  %% plot bathymetry
  subplot(322)
  trisurf(trdata, 'CData', bat.data);
  view(0,90);
  % shading interp;
  colormap(jet);
  colorbar;
  axis equal
  set(gca, 'XLim', x_lim, 'YLim', y_lim)
  %set(gca, 'DataAspectRatio', [1 1 1])
  title('bathymetry')

  %% plot velocity field (vector field)
  subplot(323)
  h = trimesh(trdata);
  set(h, 'EdgeColor', 0.8*[1 1 1])
  hold on
  quiver3(gr.elmtcenter(1,:), gr.elmtcenter(2,:), ones([1, gr.nelmts]), ...
          v_x.data', v_y.data', zeros([1, gr.nelmts]))
  hold off
  view(0,90);
  axis equal
  set(gca, 'XLim', x_lim, 'YLim', y_lim)
  title('velocity')

  %% plot velocity field (x-component)
  subplot(324)
  trisurf(trdata, 'CData', v_x.data);
  view(0,90);
  axis equal
  set(gca, 'XLim', x_lim, 'YLim', y_lim)
  caxis([minu maxu]+1)
  colorbar;
  title('x-velocity')

  %% plot velocity field (y-component)
  subplot(325)
  trisurf(trdata, 'CData', v_y.data);
  view(0,90);
  axis equal
  set(gca, 'XLim', x_lim, 'YLim', y_lim)
  caxis([minu maxu])
  colorbar;
  title('y-velocity')

  %% plot refinement levels
  subplot(326)
  trisurf(trdata, 'CData', lev.data);
  view(0,90);
  axis equal
  colorbar;
  set(gca, 'XLim', x_lim, 'YLim', y_lim)
  title('refinement level')

  if istep == 1
    pause()
  else
    pause(0.1)
  end
end
