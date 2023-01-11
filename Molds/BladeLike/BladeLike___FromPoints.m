function [X,Y,Z,z,F,DT,MoldEdge] = BladeLike___FromPoints(Inp)
% Blade-like mold defined from points and interpolation. This initially was
% used to define the mold geometry but the meshing in Abaqus yields a nicer
% mesh with more regular triangles

YLimits = [-0.2 4.2];
XLimits_root = [-0.7 0.7];
XLimits_aero = [-0.7 0.7];

YNetBdInset = [0.2 0.2];
XNetBdInset = [0.05 0.05];

nDiv = 100;
nDiv2 = 100;
Res = 0.010;

t = linspace(205,335,nDiv);

x_root = -cosd(t) * (sqrt(2)/2) / cosd(t(1));
z_root = sind(t);
y_root = YLimits(1)*ones(1,nDiv);

Pt = [-sqrt(2)/2-0.02 -sqrt(2)/2, -0.35, 0.00,   0.25,  0.45, 0.57  sqrt(2)/2 sqrt(2)/2+0.04; ...
    -sqrt(2)/2 -sqrt(2)/2, -1.00, -0.95, -0.80, -0.75,  -0.80   -0.95 -0.95];

curve = cscvn(Pt);

Pt2 = fnval(curve,linspace(curve.breaks(1),curve.breaks(end),nDiv));

x_aero = Pt2(1,:);
z_aero = Pt2(2,:);
y_aero = YLimits(2)*ones(1,nDiv);

x = zeros(nDiv,nDiv2);
y = zeros(nDiv,nDiv2);
z = zeros(nDiv,nDiv2);
for i = 1:nDiv
    x(i,1:nDiv2) = linspace(x_root(i),x_aero(i),nDiv2);
    y(i,1:nDiv2) = linspace(y_root(i),y_aero(i),nDiv2);
    z(i,1:nDiv2) = linspace(z_root(i),z_aero(i),nDiv2);
end

F = scatteredInterpolant(x(:),y(:),z(:),'linear','linear');
DT = delaunayTriangulation(x(:),y(:));

%[X,Y] = meshgrid(...
%    min([XLimits_root(1) XLimits_aero(1)])+XNetBdInset(1):Res:...
%    max([XLimits_root(2) XLimits_aero(2)])-XNetBdInset(2) , ...
%    YLimits(1)+YNetBdInset(1):Res:YLimits(2)-YNetBdInset(2));

[X,Y] = meshgrid(min(x(:)):Res:max(x(:)) , min(y(:)):Res:max(y(:)));

% Define the net boundary in xy-vertices
LowerRight = [0.66 0.0];
MiddleRight = [0.66 1.0];
UpperRight = [0.70 4.0];
UpperLeft = [-0.70 4.0];
MiddleLeft = [-0.66 1.0];
LowerLeft = [-0.66 0.0];
NetBd = [LowerRight ; UpperRight ; ...
    UpperLeft ; LowerLeft];

[In,On] = inpolygon(X,Y,NetBd(:,1),NetBd(:,2));
X(~(In | On)) = NaN;
Y(~(In | On)) = NaN;
Z = F(X,Y);

% Define the left, right and bottom mold edge
%MoldEdge{1} = [X(:,1) Y(:,1) Z(:,1)];
%MoldEdge{2} = [X(:,end) Y(:,end) Z(:,end)];
%MoldEdge{3} = [X(1,:)' Y(1,:)' Z(1,:)'];

% Define the left, right and bottom mold edge
x_left = linspace(LowerLeft(1),UpperLeft(1),50);
y_left = linspace(LowerLeft(2),UpperLeft(2),50);

x_right = linspace(LowerRight(1),UpperRight(1),50);
y_right = linspace(LowerRight(2),UpperRight(2),50);

x_bottom = linspace(LowerLeft(1),LowerRight(1),25);
y_bottom = linspace(LowerLeft(2),LowerRight(2),25);

MoldEdge{1} = [x_left' y_left' F(x_left',y_left')];
MoldEdge{2} = [x_right' y_right' F(x_right',y_right')];
MoldEdge{3} = [x_bottom' y_bottom' F(x_bottom',y_bottom')];

%%Create stl?
%%First create a 3D triangulation with the same connectivity as the
%%2D Delaunay triangulation
%faces  = delaunay(X(~isnan(X)),Y(~isnan(Y)));
%TR = triangulation(faces,[X(~isnan(X)) Y(~isnan(Y)) Z(~isnan(Z))]);
%%Write
%stlwrite(TR,'BladeLikeMoldSTL_ext.stl','text')