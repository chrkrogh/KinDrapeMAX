function [X,Y,Z,z,F,DT,MoldEdge] = Hemisphere___sph(Inp)
% Hemisphere defined in spherical coordinates
% This gives some badly scaled triangles in the mesh...

[Theta,Phi] = meshgrid(linspace(0,2*pi -(2*pi/100),100),...
    linspace(1e-4,pi/2-pi/20,50));
X = 1*cos(Theta).*sin(Phi);
Y = 1*sin(Theta).*sin(Phi);
Z = 1*cos(Phi);

F = scatteredInterpolant(X(:),Y(:),Z(:),'linear','linear');
DT = delaunayTriangulation(X(:),Y(:));

z = Z(:);

MoldEdge = [];
end