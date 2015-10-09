function u = FEMcylinder(S,coordinates,elements,neumann,dirichlet,meshDim,force)
  % Implementation of the Finite Element Method
  % Input:
    % Effective compliance matrix of material: S
    % Coordinates of element nodes: coordinates
    % Element numbering: elements
    % Elements with Neumann conditions: neumann
    % Elements with Dirichlet conditions: dirichlet
    % Mesh dimensions: meshdim = [r_o,r_i,dr,h,n]
      % where meshdim(1:4) is in meters and n is the number of equally spaced nodes
    % Surface force applied to an individual element in Newtons: force
  % Output: displacement vector u

  A = sparse(3*size(coordinates,1),3*size(coordinates,1));
  b = zeros(3*size(coordinates,1),1);

  display('Assembling global stiffness matrix...')

  %Assembly
  for j = 1:size(elements,1)
    I = 3*elements(j,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0];
    A(I,I) = A(I,I) +stima(coordinates(elements(j,:),:),S);
  end

  % Volumeforces - commented out if neglected
  % for j = 1:size(elements,1)
  %   I = 3*elements(j,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0]; 
  %   fs = f(sum(coordinates(elements(j,:),:))/4)';
  %   b(I) = b(I) +det([1,1,1,1;coordinates(elements(j,:),:)'])*[fs;fs;fs;fs]/24;
  % end

  display('Applying boundary conditions...')

  %Neumann conditions
  if ~isempty(neumann)
    for j = 1:size(neumann,1)
      % Normal vector of jth surface triangle
      n = cross( coordinates(neumann(j,2),:)-coordinates(neumann(j,1),:), ...
    coordinates(neumann(j,3),:)-coordinates(neumann(j,1),:));
      % Area of triangle
      dA = norm(n);
      center=sum(coordinates(neumann(j,:),:))/3;
      % Compute surface forces (gm in force per unit area)
      gm = g(center,n/dA,dA,meshDim,force)';
      % Construct right hand side
      I = 3*neumann(j,[1,1,1,2,2,2,3,3,3])-[2,1,0,2,1,0,2,1,0];
      b(I) = b(I) +dA*[gm;gm;gm]/6;   
    end
  end

  %Dirichlet conditions
  dirichletnodes = unique(dirichlet);
  [W,M] = u_d(coordinates(dirichletnodes,:));
  B = sparse(size(W,1),3*size(coordinates,1));
  for k = 0:2
    for l = 0:2
      B(1+l:3:size(M,1),3*dirichletnodes-2+k) = diag(M(1+l:3:size(M,1),1+k));
    end
  end
  mask=find(sum(abs(B)'));
  A = [A,B(mask,:)';B(mask,:),sparse(length(mask),length(mask))];
  b = [b;W(mask,:)];

  %Calculation of the solution
  x = A\b;
  u = x(1:3*size(coordinates,1));
