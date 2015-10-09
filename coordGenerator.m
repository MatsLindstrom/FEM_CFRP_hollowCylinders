function [coordinates,elements,neumann,dirichlet] = coordGenerator(meshDim)
% Produces coordinates and triangulation for a hollow cylinder
% Input: scalar vector meshdim = [r_o,r_i,dr,h,n]
% Output:
%   coordinates = 3*(n+1)*n X 3 matrix of cartesian coordinates
%     - each row number represents a node number
%     - n is the number of nodes in each direction
%   elements = N X 4 matrix defining node numbers of N tetrahedrons
%   neumann = N X 3 matrix defining node numbers of N Neumann surface elements
%   dirichlet N X 3 matrix defining node numbers of N Dirichlet surface elements

r_o = meshDim(1); r_i = meshDim(2); dr = meshDim(3); h = meshDim(4);
n = meshDim(5); % Number of equally spaced nodes in circumferential and z direction

% Generate vector of radial coordinates
r_c = r_i + dr/2;
r = [r_o,r_i,r_c];
r = repelem(r,(n+1)*n);

% Generate vector of circumferential coordinates
theta = 2*pi/n:2*pi/n:2*pi;
theta = repmat(theta,[1,n+1]);
theta = repmat(theta,[1,3]);

% Generate z-axis coordinates
z = 0:h/n:h;
z = repelem(z,n);
z = repmat(z,[1,3]);

% Assemble polar coordinates of cylinder
polarCoords = zeros(3*(n+1)*n,3);
polarCoords(:,1) = r';
polarCoords(:,2) = theta';
polarCoords(:,3) = z';
save('polarCoords.mat','polarCoords')

% Transform to cartesian cartesian coordinates
[x,y,z] = pol2cart(polarCoords(:,2),polarCoords(:,1),polarCoords(:,3));
coordinates = [x,y,z];
save('coordinates.mat','coordinates')

% Generate triangulation (tetrahedrons)
tri = delaunayTriangulation(coordinates);

% Find elements not in cylindrical shell
hullCc = circumcenter(tri);
toRemove = [];
normHullCc = sqrt(hullCc(:,1).^2 + hullCc(:,2).^2);
toRemove = find(normHullCc<r_i);

% Remove elements not in shell from connectivity list
connectivityList = tri.ConnectivityList;
connectivityList(toRemove,:) = [];
elements = connectivityList;
save('elements.mat','elements')

% Find surface elements
hullCon = convexHull(tri);
hullTri = triangulation(hullCon,tri.Points);
hullNormals = faceNormal(hullTri);
hullCc = circumcenter(hullTri);
hullCc = sqrt(hullCc(:,1).^2+hullCc(:,2).^2);

% Find elements on which Dirichlet conditions are applied
bottomNormal = [0 0 -1];
dirichletIndex = find(ismember(hullNormals,bottomNormal,'rows') & hullCc>r_i);
dirichlet = hullCon(dirichletIndex,:);
save('dirichlet.mat','dirichlet')

% Find elements on which Neumann conditions are applied
neumannIndex = find((~ismember(hullNormals,bottomNormal,'rows')) & hullCc>r_i);

allR = [cylRadius(coordinates(elements(:,1),:)),cylRadius(coordinates(elements(:,2),:)),cylRadius(coordinates(elements(:,3),:)),cylRadius(coordinates(elements(:,4),:))];

condition = r_i+dr/4;
[row,col] = find((allR(:,1)>condition & allR(:,2)>condition) |...
    (allR(:,1)>condition & allR(:,3)>condition) |...
    (allR(:,1)>condition & allR(:,4)>condition) |...
    (allR(:,2)>condition & allR(:,3)>condition) |...
    (allR(:,2)>condition & allR(:,4)>condition) |...
    (allR(:,3)>condition & allR(:,4)>condition));
innerTetras = elements;
innerTetras(row,:) = [];

tempRad = [cylRadius(coordinates(innerTetras(:,1),:)),cylRadius(coordinates(innerTetras(:,2),:)),cylRadius(coordinates(innerTetras(:,3),:)),cylRadius(coordinates(innerTetras(:,4),:))];
[row1, col1] = find(tempRad(:,1)>condition);
[row2, col2] = find(tempRad(:,2)>condition);
[row3, col3] = find(tempRad(:,3)>condition);
[row4, col4] = find(tempRad(:,4)>condition);
tempTetras=num2cell(innerTetras,2);

for i = 1:length(row1)
    tempTetras{row1(i)}(1) = [];
end
for i = 1:length(row2)
    tempTetras{row2(i)}(2) = [];
end
for i = 1:length(row3)
    tempTetras{row3(i)}(3) = [];
end
for i = 1:length(row4)
    tempTetras{row4(i)}(4) = [];
end

tempTri = cell2mat(tempTetras);
neumann = [hullCon(neumannIndex,:);tempTri];
save('neumann.mat','neumann')
