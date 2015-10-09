function stima = stima(vertices,S)
  % Assemble element stiffness matrix
  % Input:
    % Coordinates of element nodes: vertices
    % Material compliance matrix: S
  % Output: Element stiffness matrix stima

  PhiGrad = [1,1,1,1;vertices']\[zeros(1,3);eye(3)];
  R = zeros(6,12);
  R([1,4,5],1:3:10) = PhiGrad';
  R([4,2,6],2:3:11) = PhiGrad';
  R([5,6,3],3:3:12) = PhiGrad';
  Q = inv(S);
  stima = det([1,1,1,1;vertices'])/6*R'*Q*R;
