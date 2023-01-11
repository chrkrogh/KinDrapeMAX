function [X,Y,Z,z,F,DT,MoldEdge] = Flat(Inp)
% Flat mold 3 x 10 m

x_max = 3;
y_max = 10;

F = @(x,y) zeros(size(x));

[X,Y] = meshgrid(0:0.1:x_max,0:0.1:y_max);
Z = F(X,Y);
DT = delaunayTriangulation(X(:),Y(:));
z = Z(:);

% Define the left, right and bottom mold edge
MoldEdge{1} = [X(:,1) Y(:,1) Z(:,1)];
MoldEdge{2} = [X(:,end) Y(:,end) Z(:,end)];
MoldEdge{3} = [X(1,:)' Y(1,:)' Z(1,:)'];
end