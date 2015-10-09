function laminate = stackSim(stacks,E,nu,G,meshDim,coordinates,elements,dirichlet,neumann,force)
  % Run FEM for set of stacks
  % Input:
    % Stack configuration (ply angles in radians): stacks
    % Elastic modulus: E = [E_r, E_phi, E_z] in GPa
    % Poisson's ratio: nu = [nu_phi_z, nu_r_z, nu_r_phi]
    % Shear modulus: G = [G_phi_z, G_r_z, G_r_phi] in GPa
    % Mesh dimensions: meshdim = [r_o,r_i,dr,h,n]
      % where meshdim(1:4) is in meters and n is the number of equally spaced nodes
    % Coordinates of element nodes: coordinates
    % Element numbering: elements
    % Elements with Dirichlet conditions: dirichlet
    % Elements with Neumann conditions: neumann
    % Surface force applied to an individual element in Newtons: force
  % Output:
    % Structure containing fields with stack configuration and resulting maximum deflection: laminate

  nStacks = size(stacks,3);
  nAngles = size(stacks,1);
  laminate = struct;
  legendEntries = cell(1,size(stacks,3));
  r_o = meshDim(1);
  r_i = meshDim(2);
  dr = meshDim(3);

  for k = 1:nStacks
    display('-------------- NEW TYPE OF STACKING SEQUENCE --------------')
    laminate(k).stack = stacks(:,:,k);
    maxDisplacement = zeros(nAngles,1);

    for i = 1:nAngles
      phi = stacks(i,:,k);
      display(['New ply angles: ', num2str(radtodeg(phi))])
      display('Computing effective material stiffness matrix...')
      S = effCompMatrix(phi,meshDim,E,nu,G);
      u = FEMcylinder(S,coordinates,elements,neumann,dirichlet,meshDim,force);
      maxDisplacement(i) = max(u(2:3:end))*10^3;
      display(['Maximum displacement: ', num2str(maxDisplacement(i)), ' mm']) % Maximum displacement in direction of applied force
    end
    laminate(k).deflection = maxDisplacement;

    figure(3)
    plot(radtodeg(abs(stacks(:,1,k))),maxDisplacement,'k')
    xlabel('Ply angle [deg]')
    ylabel('Deflection [mm]')
    legendEntries(1,k)={['Case ',num2str(k)]};
    hold on
  end

  legend(legendEntries{1,:});