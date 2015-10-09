function out = cylRadius(coords)
% Returns the radius of a cylinder given 3-dimensional cartesian coordinates

out = sqrt(coords(:,1).^2 + coords(:,2).^2);
