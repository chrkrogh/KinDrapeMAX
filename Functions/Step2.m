function Dra = Step2(Dra,Grid,d,OrgNode,F)
% Step 2: Constrained nodes

% Convergence tolerance
ConTol = 1e-5;

Node = Dra.Node;
P = NaN(prod(Grid-1),4,5);
Dir1 = [1 0 -1 0 ; 0 1 0 -1]';
Dir2 = [0 -1 0 1; 1 0 -1 0]';

% Number of constrained nodes in each quadrant
nConNode = Dra.nIniNode([1 2 ; 3 2 ; 3 4 ; 1 4]);
for i = 1:4
    for j = OrgNode(1) + (0:nConNode(i,1)-1)*(Dir1(i,1)+Dir2(i,1))
        for k = OrgNode(2) + (0:nConNode(i,2)-1)*(Dir1(i,2)+Dir2(i,2))
            % Def. single idx and solver input. Call circle intersec. fun.
            [Idx,CellNo] = LinearCellIdx(Grid,j,k,Dir1,Dir2,i);
            Node(Idx(3,:)) = MoldCircIntersecFun(Node(Idx),F,d,[],ConTol);
            % Put current cell coord. and shear in P array for patch plot
            P(CellNo,1:4,1:3) = Node(Idx);
            [~,~,P(CellNo,1:4,4)] = ShearFun([],Node(Idx),F,0,i);
        end
    end
end

% Compute the fiber orientation deviations
P(:,1:4,5) = FibOriFun(Node);
%P(:,1:4,5) = FibOriFun_legacy(Node);

Dra.P = P;
Dra.Node = Node;
end

function IntersecPt = MoldCircIntersecFun(Vert,F,d,Ang,ConTol)
% Location of new point based on intersection of circle and mold surface
% C, R: center and radius. Vec1, Vec2: perpend. unit vec. spanning circle
if all(isnan(Vert(2,:))) % Step 1
    % Circle is constructed along initial angle, centered in Vert1
    C = Vert(1,:);
    R = d;
    Vec1 = -[cosd(90+Ang) sind(90+Ang) 0];
    Vec2 = [0 0 1];
else % Step 3
    % Circle is intersection of two spheres, centered in Vert 2 and Vert 4
    C = (Vert(2,:) + Vert(4,:))/2;
    R = sqrt(d^2 - norm(C-Vert(2,:))^2);
    Vec1 = (Vert(1,:)-C)/norm(Vert(1,:)-C);
    CircleAxis = (C-Vert(2,:))/norm(C-Vert(2,:));
    Vec2 = cross(Vec1,CircleAxis)/norm(cross(Vec1,CircleAxis));
end
% Find the intersection between the circle and the surface using bisection
Theta = [pi/2 2*pi - pi/2];
IntersecPt = NaN(1,3);
for i = 1:1e3
    % Compute middle point
    Theta_mid = (Theta(1)+Theta(2))/2;
    % Circle pt. in 3D based on center, radius, 2 perp. vectors and ang.
    CircCoor = C + R*cos(Theta_mid).*Vec1 + R*sin(Theta_mid).*Vec2;
    % Added robustness: Sometimes when the mold is bad (e.g. due to
    % extrapolation) R can become complex.
    if ~isreal(CircCoor)
        break
    end
    CircCoor_mold = F(CircCoor(1),CircCoor(2));
    FunVal_mid = CircCoor(3) - CircCoor_mold;
    % Stop or adjust interval based on function value at midpoint
    if abs(FunVal_mid) < ConTol
        IntersecPt = [CircCoor(1) CircCoor(2) CircCoor_mold];
        break
    elseif FunVal_mid > 0
        Theta(1) = Theta_mid;
    else
        Theta(2) = Theta_mid;
    end
end
end

function [Obj, Vert, Shear] = ShearFun(Dx,Vert,F,PreShear,No)
% Calculate unknown node(s) and shear ang. Output: shear ang. sum (obj. in
% step 2), new vertex coord. and vector of shear ang. (for P array)
if length(Dx) == 4 % Step 2
    Vert(3,:) = [Vert(2,1:2)+Dx(1:2) F(Vert(2,1)+Dx(1),Vert(2,2)+Dx(2))];
    Vert(4,:) = [Vert(1,1:2)+Dx(3:4) F(Vert(1,1)+Dx(3),Vert(1,2)+Dx(4))];
end
% Calculate shear angles using edge vectors u and v. Calculate sum
u = Vert([2 3 4 1],:)' - Vert';
v = Vert([4 1 2 3],:)' - Vert';

% Calculate dot product, cross product and norm of cross product
DotProd = sum(u.*v,1);
CrossProd = [u(2,:).*v(3,:) - u(3,:).*v(2,:) ; ...
             u(3,:).*v(1,:) - u(1,:).*v(3,:) ; ...
             u(1,:).*v(2,:) - u(2,:).*v(1,:)];
NormCP = sqrt(CrossProd(1,:).^2 + CrossProd(2,:).^2 + CrossProd(3,:).^2);

% Calculate the shear angles
Shear = (atan2d(NormCP,DotProd)-90).*[1 -1 1 -1]*(-1)^No;
Obj =  abs(sum(Shear - PreShear));
end

function FibOriDev = FibOriFun(Node)
% This function takes the Node array and computes the fiber orientations,
% i.e. the angle between the warp fibers (longitudinal cell edges) and the
% nominal fiber angle projected to each element plane

% The nominal fiber direction, i.e. y-axis
NomDir = [0 1 0];

% Compute the warp fiber angle vector based on the left longitudinal edge 
% of each cell. That is, take the difference along the 2nd dimension
% (columns) from rows 1 to end-1. Each slice/page will then contain a
% vector component
WarpVec = diff(Node(1:end-1,:,1:3),1,2);

% Compute a vector along the weft direction (for cross product calculation 
% later) That is, take the difference along the 1st dimension
% (rows) from columns 1 to end-1. Each slice/page will then contain a
% vector component
WeftVec = diff(Node(:,1:end-1,1:3),1,1);

% Compute the element normal, i.e. by crossing the WarpVec and WeftVec
% arrays along the 3rd dimension
ElemNormal = cross(WeftVec,WarpVec,3);

% Normalize the normal vector by division with the 2-norm (again computed
% along the 3rd dimension)
ElemNormal = ElemNormal ./ vecnorm(ElemNormal,2,3);

%Reshape NomDir to be the same shape as WarpVec
NomVec = reshape(NomDir,1,1,3) .* ones(size(WarpVec));

% Compute the projection of the NomDir onto the element plane
Proj = NomVec - dot(ElemNormal,NomVec,3).*ElemNormal;

% Compute the signed angle between the projection and the WarpVec
AngDev = atan2d(dot(cross(Proj,WarpVec,3),ElemNormal,3),...
    dot(Proj,WarpVec,3));

% The number of cells
nCells = prod(size(Node(:,:,1))-1);

% Reshape the AngDev array so it fits in P (with rows 1:nCells). This also 
% includes repeating the computed cell value 4 times, i.e. so that the 4 
% cell nodes have the same value
FibOriDev = repmat(reshape(AngDev,nCells,1),1,4);
end


function FibOriDev = FibOriFun_legacy(Node)
% This is the legacy version of the FibOriFun. It either computes the fiber
% orientations by computing the angles between the warp fibers 
% (longitudinal cell edges) either by projection to the XY-plane or with 
% 3D angles. None of these approaches fit with the Abaqus convention...

nCells = prod(size(Node(:,:,1))-1);

XYProjection = 0;

% The nominal fiber direction, i.e. y-axis
NomDir = [0 1 0];

% Compute the warp fiber angle vector based on the left longitudinal edge 
% of each cell. That is, take the difference along the 2nd dimension
% (columns) from rows 1 to end-1. Each slice/page will then contain a
% vector component
WarpVec = diff(Node(1:end-1,:,1:3),1,2);
if XYProjection == 1
    WarpVec(:,:,3) = 0.0;
end

% Compute the warp fiber deviations for each course. Reshape NomVec to be
% the same shape as WarpVec
NomVec = reshape(NomDir,1,1,3) .* ones(size(WarpVec));
% Compute the angle projected to the XY-plane
AngDev = atan2d(vecnorm(cross(WarpVec,NomVec),2,3),...
    dot(WarpVec,NomVec,3));

% Reshape the AngDev array so it fits in P (with rows 1:nCells). This also 
% includes repeating the computed cell value 4 times, i.e. so that the 4 
% cell nodes have the same value
FibOriDev = repmat(reshape(AngDev,nCells,1),1,4);
end