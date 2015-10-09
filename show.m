function show(dirichlet,neumann,coordinates)
  % Visual representaion of the cylinder deformation
  % Input:
    % Dirichlet surfaces: dirichlet
    % Neumann surfaces: neumann
    % Coordinates of all nodes: coordinates

  E=[dirichlet;neumann];

  map = repmat([0.1:0.05:1]',[1,3]);
  map = [map;flipud(map)];
  colormap(map)
  C = coordinates(:,1);
  trisurf(E,coordinates(:,1),coordinates(:,2),coordinates(:,3), C,...
          'EdgeColor','none');
  view(-50,35)
  axis equal;
  xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
