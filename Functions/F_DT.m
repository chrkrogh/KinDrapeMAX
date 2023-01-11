function z = F_DT(DT,Z,x,y)
% This function returns the (linearly) interpolated z value for the
% Delaunay triangulated mesh. Inputs are:
% - DT: delaunayTriangulation object
% - Z: array with z coordinates at the x and y coordinates used for the
%   triangulation.
% - x,y: query coordinate

if nargin == 3
    % x and y coordinates in two column matrix
    P = x;
elseif nargin == 4
    % separate x and y coordinates
    P = [x(:) y(:)];
end

% Find triangles that enclose query points. The function returns the
% triangle ID and the barycentric coordinates
[ti,bc] = pointLocation(DT,P);

% Find the values of Z at the vertices of each enclosing triangle.
if any(isnan(ti))
    triVals(~isnan(ti),1:3) = Z(DT(ti(~isnan(ti)),:));
    triVals(isnan(ti),1:3) = NaN(length(find(isnan(ti))),3);
else
    triVals = Z(DT(ti,:));
end

%Calculate the sum of the weighted values using the dot product
z = dot(bc',triVals')';
end