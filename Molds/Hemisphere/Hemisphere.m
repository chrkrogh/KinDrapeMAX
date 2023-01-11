function [X,Y,Z,z,F,DT,MoldEdge] = Hemisphere(Inp)
% Hemisphere defined in cartesian coordinates
[X,Y] = meshgrid(linspace(-1,1,100));
MaskIdx = X.^2 + Y.^2 > 0.99.^2;
X(MaskIdx) = NaN;
Y(MaskIdx) = NaN;
Z = sqrt(1 - X.^2 - Y.^2);

F = scatteredInterpolant(X(~isnan(X)),Y(~isnan(Y)),Z(~isnan(Z)),...
    'linear','linear');
DT = delaunayTriangulation(X(~isnan(X)),Y(~isnan(Y)));

z = Z(~isnan(Z));

% Net boundary (hack the system so that net boundary makes a circle with 
% radius R_NB but still follows the requirements to the left, right 
% (ascending y-coordinates and bottom (ascending x-coordinates) mold edges
R_NB = 0.95;

Theta1 = linspace(180,90,100)';
Theta2 = linspace(0,90,100)';
Theta3 = linspace(180,360,100)';

% "Left" mold edge
MoldEdge{1} = [R_NB*cosd(Theta1) R_NB*sind(Theta1)];
MoldEdge{1}(:,3) = F(MoldEdge{1}(:,1),MoldEdge{1}(:,2));

% "Right" mold edge
MoldEdge{2} = [R_NB*cosd(Theta2) R_NB*sind(Theta2)];
MoldEdge{2}(:,3) = F(MoldEdge{2}(:,1),MoldEdge{2}(:,2));

% "Bottom" mold edge
MoldEdge{3} = [R_NB*cosd(Theta3) R_NB*sind(Theta3)];
MoldEdge{3}(:,3) = F(MoldEdge{3}(:,1),MoldEdge{3}(:,2));
end