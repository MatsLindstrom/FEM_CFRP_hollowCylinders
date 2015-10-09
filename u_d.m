function [W,M] = u_d(x)
  % Generate Dirichlet boundary conditions

  M = zeros(3*size(x,1),3);
  W = zeros(3*size(x,1),1);
  M(1:3:end,1) = 1;
  M(2:3:end,2) = 1;
  M(3:3:end,3) = 1;
