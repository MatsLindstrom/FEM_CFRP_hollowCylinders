% main.m is a Matlab script generating:
%   - The mesh of a hollow cylindrical coordinate system from coordGenerator.m,
%   - Homogenized material properties using effCompMatrix.m, where
%     cylCompMat.m transforms the local coordinate system to its global equivalent
%   - A FEM simulation using FEMcylinder.m, where
%     - stima.m assembles the element stiffness matrix,
%     - f.m describes volume forces such as gravitational forces,
%     - g.m describes surface (traction) forces,
%     - u_d.m describes Dirichlet conditions
%   - Results visualized using show.m

% A stack of different ply configurations can be assembled, simulated and visualized by running stackGenerator.m, stackSim.m and plotLaminateStacking.m 

close all

r_o = 10e-3; % Outer radius of cylinder
dr = 1.5e-3; % Thickness of cylinder
h = 460e-3; % Height of cylinder
n = 20; % Number of equally spaced points
force = 1151; % Applied force in N
r_i = r_o - dr;

meshDim = [r_o,r_i,dr,h,n];

display('Generating mesh...')
[coordinates,elements,neumann,dirichlet] = coordGenerator(meshDim);

nAngles = 2;
nPlies = 9;
scenario = 1; % 1,2 or 3

E = [135,8,8];
nu= [0.27,0.27,0.49];
G = [3.8,3.8,2.7];

phi = zeros(1,9);
S = effCompMatrix(phi,meshDim,E,nu,G);
u = FEMcylinder(S,coordinates,elements,neumann,dirichlet,meshDim,force);
display(['Maximum displacement: ', num2str(max(u(2:3:end))*10^3), ' mm']) % Maximum displacement in direction of applied force

figure(1)
show(dirichlet,neumann,coordinates.*10^3);
magn=1.0;
defCart = [(magn*u(1:3:size(u,1))+coordinates(:,1)).*10^3, ...
            (magn*u(2:3:size(u,1))+coordinates(:,2)).*10^3, ...
            flipud((magn*u(3:3:size(u,1))+coordinates(:,3))).*10^3];
figure(2)
show(dirichlet,neumann,defCart); 

stacks = stackGenerator(nAngles,nPlies,scenario);
laminate = stackSim(stacks,E,nu,G,meshDim,coordinates,elements,dirichlet,neumann,force);
