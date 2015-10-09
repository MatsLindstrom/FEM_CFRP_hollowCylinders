This Matlab program homogenizes and solves for the material response of cylindrical hollow composites consisting of helically oriented fiber reinfoced plies under pure bending loads.

%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%

The program was built as a part of a Master thesis at Chalmers University of Technology, Department of Applied Physics, 2015.

Title of thesis: "Improved material performance in floorball sticks - Implementation of the Finite Element Method for helically fiber-oriented laminae in cylindrical hollow composites"
Author: Mats Lindström, mats.at.lindstrom@gmail.com

Remark: This program is based on a supplement to the paper
"Matlab-Implementation of the Finite Element Method in Elasticity" by J. Alberty, C. Carstensen, S. A. Funken, and R. Klose.

%%%%%%%%%%%%%%%%%%%%%%%%%% Files %%%%%%%%%%%%%%%%%%%%%%%%%%

main.m is a Matlab script generating:
  - The mesh of a hollow cylindrical coordinate system from coordGenerator.m,
  - Homogenized material properties using effCompMatrix.m, where
    cylCompMat.m transforms the local coordinate system to its global equivalent
  - A FEM simulation using FEMcylinder.m, where
    - stima.m assembles the element stiffness matrix,
    - f.m describes volume forces such as gravitational forces,
    - g.m describes surface (traction) forces,
    - u_d.m describes Dirichlet conditions
  - Results visualized using show.m

A stack of different ply configurations can be assembled, simulated and visualized by running stackGenerator.m, stackSim.m and plotLaminateStacking.m 


%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%

Functions for single configurations

  % Coordinates and triangulation of hollow cylinder
  [coordinates,elements,neumann,dirichlet] = coordGenerator(meshDim)

  % Effective material compliance matrix of homogenized material 
  S = effCompMatrix(phi,meshDim,E,nu,G)
  
  % Local ply coordinate system to global
  S = cylCompMat(phi, E, nu, G)
  
  % Displacement via FEM
  u = FEMcylinder(S,coordinates,elements,neumann,dirichlet,meshDim,force)

  % Element stiffness matrix assembly
  stima = stima(vertices,S)

  % Volume forces
  volforce = f(x)

  % Surface forces
  sforce = g(x,normal,dA,meshDim,force)

  % Dirichlet conditions
  [W,M] = u_d(x)

  % Visualize triangulation
  show(dirichlet,neumann,coordinates)

Special functions for multiple configurations

  % Set of stacks consisting of different angles, number of plies, etc.
  stacks = stackGenerator(nAngles,nPlies,scenario)

  % Run FEM for different stacks
  laminate = stackSim(stacks,E,nu,G,meshDim,coordinates,elements,dirichlet,neumann,force)

  % Plot obtained results from stack simulations 
  SCRIPT: plotLaminateStacking.m
