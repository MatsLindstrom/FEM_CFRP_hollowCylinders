function sforce = g(x,normal,dA,meshDim,force)
  % Application of surface forces
  % Input:
    % Cartesian coordinates of considered element: x = [x,y,z]
    % Surface normal: normal = [x,y,z]
    % Area of surface element: dA
    % Mesh dimensions: meshdim = [r_o,r_i,dr,h,n]
      % where meshdim(1:4) is in meters and n is the number of equally spaced nodes
    % Surface force applied to an individual element in Newtons: force
  % Ouput: Applied surface force

  r_o = meshDim(1);
  dr = meshDim(3);
  h = meshDim(4);
  n = meshDim(5);

  sforce=zeros(size(x,1),3);

  if (x(3)>(h-h/(2*(n+1))) && normal(2)<-0.95 && x(2)<-(r_o-dr/4) && x(1)>0)
  	f=force*10^(-9);
  	fdA=f/dA;
  	display(['Applied force: ',num2str(force),' N'])
  	sforce(1,2)=fdA;
  end
