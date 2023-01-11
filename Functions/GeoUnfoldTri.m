function [GeoLine, MissingLength] = ...
    GeoUnfoldTri(xy0,Ang,d,nNode,DT,Z,F,IncPlt,Set,i)
% Geodesic line by unfolding triangles to common plane

MissingLength = 0;
% Compute the target length
L_target = nNode * d + 0.5*d;

% Check if starting point is inside the mesh
if isnan(F_DT(DT,Z,xy0))
    % Use extrapolation information
    GeoLine(1,1:3) = [xy0(1) xy0(2) F(xy0(1),xy0(2))];
    IncMove = d;
    MaxIter = ceil(L_target/IncMove) + 1;
    %FirstTriPt = [];
    for Ctr = 2:MaxIter
        % Go IncMove along initial angle
        % Compute new point
        xy_new = xy0(1:2) + (Ctr-1)*IncMove*[cosd(Ang) sind(Ang)];
        xy_new(3) = F(xy_new(1),xy_new(2));
        GeoLine(Ctr,1:3) = xy_new;
        if ~isnan(F_DT(DT,Z,xy_new(1:2)))
            L_path = (Ctr-1)*IncMove;
            FirstTriPt = xy_new(1:2);
            % Set StartPointOnEdge to true. Trick to check all tri edges 
            % later in the code
            StartPointOnEdge = 1;
            EdgePtCollision = 0;
            break
        elseif isnan(F_DT(DT,Z,xy_new(1:2))) && Ctr == MaxIter
            return
        end
    end
elseif isfield(Set,'ForceGeoLineIniAngle') && Set.ForceGeoLineIniAngle
    % Compute first point on GeoLine 'd' away from xy0 along Ang
    % Compute z coordinate of xy0
    xy0(3) = F_DT(DT,Z,xy0);

    % First compute trial point for xy1
    xy1 = xy0(1:2) + d*[cosd(Ang) sind(Ang)];
    xy1(3) = F_DT(DT,Z,xy1);

    % Calculate slope and correct
    if ~any(isnan(xy1))
        SlopeAng = atand((xy1(3)-xy0(3))/norm(xy1(1:2)-xy0(1:2)));
        xy1 = xy0(1:2) + d*cosd(SlopeAng)*[cosd(Ang) sind(Ang)];
        xy1(3) = F_DT(DT,Z,xy1);
    end

    GeoLine(1,1:3) = xy0;
    GeoLine(2,1:3) = xy1;

    % Initialize path length
    L_path = d;

    % Set point for first triangle
    FirstTriPt = xy1(1:2);

    % Initialize counter
    Ctr = 2;
else
    % Small increment backwards along Ang
    dx = cosd(Ang-180)*1e-5;
    dy = sind(Ang-180)*1e-5;

    % Initial geodesic path from starting point and angle
    GeoLine(2,1:3) = [xy0 F_DT(DT,Z,xy0)];
    GeoLine(1,1:3) = [xy0 + [dx dy] F_DT(DT,Z,xy0 + [dx dy])];
    
    % Initialize path length
    L_path = 0.0;

    % Set point for first triangle
    FirstTriPt = xy0;

    % Initialize counter
    Ctr = 2;
end


% Initialize the start triangle as the new triangle
NewTriID = pointLocation(DT,FirstTriPt);
if isempty(NewTriID) || isnan(NewTriID)
    MissingLength = L_target;
    if Set.MeshOrCurveWarning
        fprintf(2,['Geodesic line reached the end of the mesh'...
            ' (generator #%d)\n\n'],i)
    end
    return
end
VertIdx = DT.ConnectivityList(NewTriID,:);

% Calculate the current vertex coordinates
VertCoor(1:3,1:2) = DT.Points(VertIdx,:);
VertCoor(1:3,3) = F_DT(DT,Z,VertCoor(1:3,1),VertCoor(1:3,2));

% Calculate unit normal of start triangle and store as previous normal
% vector
EdgeVec_prev(1,1:3) = VertCoor(2,1:3) - VertCoor(1,1:3);
EdgeVec_prev(2,1:3) = VertCoor(3,1:3) - VertCoor(2,1:3);
NormVec = cross(EdgeVec_prev(1,1:3),EdgeVec_prev(2,1:3));

% Check that the normal points outward and if not flip it
Centroid = mean(VertCoor,1);
PointInDirOfNormal = Centroid + NormVec;
if PointInDirOfNormal(3) < Centroid(3)
    NormVec = cross(-EdgeVec_prev(1,1:3),EdgeVec_prev(2,1:3));
end

% Initialize shared edge vector (values do not matter because first
% triangle is not rotated)
NewSharedEdgeVec = [0 0 1];

% Shareed vertex idx. Just pick vertex 1 and 2 in 1st iteration
NewSharedVertIdx = VertIdx(1:2);

% Different incremental plotting (debugging) possibilities
if IncPlt
    N = neighbors(DT,NewTriID);
    figure
    hold on
    % 3D
    %trimesh(DT.ConnectivityList,DT.Points(:,1),...
    %    DT.Points(:,2),Z(:),'LineWidth',0.1,'FaceColor','none')
    trimesh(DT.ConnectivityList([NewTriID N(~isnan(N))],:),DT.Points(:,1),...
        DT.Points(:,2),Z(:),'edgecolor','k')
    scatter3(xy0(1),xy0(2),F_DT(DT,Z,xy0(1:2)),'kx')
    % 2D projection
%     triplot(DT.ConnectivityList([NewTriID N],:),DT.Points(:,1),...
%         DT.Points(:,2),'b--');
%     scatter(xy0(1),xy0(2),'kx')
%     scatter(VertCoor(:,1),VertCoor(:,2),'ro')
    % axis equal
    xlabel('x'); ylabel('y'); zlabel('z')
    view(3)
end

VertOrder = [1 2; 2 3 ; 3 1];
% Loop until the path has a length defined by the number og nodes and d
% plus a small amount extra
while L_path <= L_target
    % Updates for the current iteration
    VertCoor_prev = VertCoor;
    SharedVertIdx = NewSharedVertIdx;
    SharedEdgeVec = NewSharedEdgeVec;
    NormVec_prev = NormVec;
    TriID = NewTriID;
    Ctr = Ctr + 1;
    
    % Reset
    NewTriID = [];
    IntersecPt = [];
    NewSharedVertIdx = [];
    NewSharedEdgeVec = [];
    
    % Update vertex information for current triangle
    VertIdx = DT.ConnectivityList(TriID,:);
    VertCoor(1:3,1:2) = DT.Points(VertIdx,:);
    VertCoor(1:3,3) = F_DT(DT,Z,VertCoor(1:3,1:2));
    
    % Calculate edge vectors for new triangle
    EdgeVec = diff(VertCoor([1 2 3 1],:),1);
    % Get the edge lengths
    ElemSize2D = vecnorm(EdgeVec(1:3,1:2),2,2);
    
    % Unfold new triangle to same plane as previous triangle. The first two
    % vertices are the shared vertices
    UnfoldVertCoor(1:2,1:2) = DT.Points(SharedVertIdx,1:2);
    UnfoldVertCoor(1:2,3) = F_DT(DT,Z,UnfoldVertCoor(1:2,1:2));
    
    % The last vertex of the unfolded triangle must be rotated 
    % Calculate the angle of rotation (angle between normal vectors)
    % First, calculate the current normal vector as the cross product of 
    % two edge vectors
    NormVec = cross(EdgeVec(1,1:3),EdgeVec(2,1:3));
    NormVec = NormVec/norm(NormVec);
    
    % Check that the normal points outward and if not flip it
    Centroid = mean(VertCoor,1);
    PointInDirOfNormal = Centroid + NormVec;
    if PointInDirOfNormal(3) < Centroid(3)
        NormVec = cross(-EdgeVec(1,1:3),EdgeVec(2,1:3));
    end
    
    % Calculate the fold angle 
    FoldAng = -atan2d(norm(cross(NormVec_prev,NormVec)),...
        dot(NormVec_prev,NormVec));
    
    if abs(FoldAng) > 45
        %keyboard
        MissingLength = L_target - L_path;
        if Set.MeshOrCurveWarning
           fprintf(2,['Geodesic line algorithm encountered large angle difference between\n'...
               'two mesh triangles (> 45 deg) ' ...
               'with %g %% remaining. Aborting.\n\n'], ...
               100*MissingLength/L_target) 
        end 
        break
    end
    
    % Construct unit vector to rotate about, i.e. shared edge
    RotAxis = SharedEdgeVec/norm(SharedEdgeVec);
    
    % Identify the "free" (non-shared) vertex of new triangle 
    FreeVertexIdx = MY_setdiff(VertIdx,SharedVertIdx);
    FreeVertCoor(1:2) = DT.Points(FreeVertexIdx,:);
    FreeVertCoor(3) = F_DT(DT,Z,FreeVertCoor(1:2));
    % Create the vector to rotate as the free vertex coordinate relative to 
    % the mid point on the shared edge
    RotPoint = mean(UnfoldVertCoor(1:2,1:3),1);
    VecToVert = FreeVertCoor - RotPoint;
    
    % Check both ways of rotation (+/- FoldAng), i.e. which produces the
    % most co-planar result
    RotatedVec(1:2,1:3) = 0.0;
    NormC(1:2) = 0.0;
    FoldAngVec = FoldAng*[1 -1];
    for k = 1:2
        % Rotate vector using Rodrigues' rotation formula
        RotatedVec(k,:) = VecToVert*cosd(FoldAngVec(k)) + ...
            cross(RotAxis,VecToVert)*sind(FoldAngVec(k)) + ...
            RotAxis*dot(RotAxis,VecToVert)*(1-cosd(FoldAngVec(k)));
        
        TempPoint = RotatedVec(k,:) + RotPoint;
        
        % Check how co-planar the point is with the previous triangle
        NormC(k) = FourPtCoplanarCheck(VertCoor_prev(1,:),...
            VertCoor_prev(2,:),VertCoor_prev(3,:),TempPoint);
    end
    
    % Pick the most co-planar point
    if NormC(1) < NormC(2)
        UnfoldVertCoor(3,1:3) = RotatedVec(1,:) + RotPoint;
        FoldAng = FoldAngVec(1);
    else
        UnfoldVertCoor(3,1:3) = RotatedVec(2,:) + RotPoint;
        FoldAng = FoldAngVec(2);
    end
    
    % Check which edge of the unfolded triangle intersects with the current 
    % direction of the calculated geodesic path
    if Ctr == 3 || StartPointOnEdge == 1 || EdgePtCollision == 1 
        % In first iteration, check intersections with all edges
        EdgesToLoop = 1:3;
    else
        % In subsequent iterations, only check 2 and 3, i.e. non-shared
        EdgesToLoop = 2:3;
    end
     
    % Construct 2D vector from previous GeoLine point using prev. direction
    q = GeoLine(Ctr-1,1:2);
    s = (GeoLine(Ctr-1,1:2) - GeoLine(Ctr-2,1:2));
    s = s/norm(s) * max(ElemSize2D)*2;
    
    t_iter(1:3) = 0.0;
    u_iter(1:3) = 0.0;
    Cross_rs_iter(1:3) = 0.0;
    StartPointOnEdge = 0;
    EdgePtCollision = 0;
    
    % Loop over edges
    for j = EdgesToLoop
        % Get points of current edge
        EdgePt = UnfoldVertCoor(VertOrder(j,:),:);
        
        % Check intersection with geodesic path
        % Construct 2D vector from 1st vertex to 2nd as p + r
        p = EdgePt(1,1:2);
        r = EdgePt(2,1:2) - EdgePt(1,1:2);
        
        t = Cross2D((q-p),s)/Cross2D(r,s);
        u = Cross2D((p-q),r)/Cross2D(s,r);
        
        % Store values of t, u and Cross_rs for examination
        t_iter(j) = t;
        u_iter(j) = u;
        Cross_rs_iter(j) = Cross2D(r,s);
        
        % The lines intersect if the following three conditions are met
        if Cross2D(r,s) ~= 0 && all([t u] >= 0) && all([t u] <= 1)
            % Calculate intersection point (3D)
            IntersecPt = EdgePt(1,:)+t*(EdgePt(2,:) - EdgePt(1,:));
            
            if norm(IntersecPt - GeoLine(2,:)) == 0 && Ctr == 3
                % The starting point is located on an edge
            
                StartPointOnEdge = 1;

            elseif any(vecnorm(IntersecPt - EdgePt,2,2) == 0)
                % The intersection point is an edge point / vertex
                
                EdgePtCollision = 1;

            else
                
                % Find new triangle that shares the intersecting edge
                AllVertIdx = [SharedVertIdx FreeVertexIdx];
                NewSharedVertIdx = AllVertIdx(VertOrder(j,:));
                NewSharedEdgeVec = EdgePt(2,1:3) - EdgePt(1,1:3);
                
                % Find the new triangle ID by the shared vertex indices
                V = edgeAttachments(DT,NewSharedVertIdx);
                NewTriID = MY_setdiff(V{1,1},TriID);
                
            end
            
            break
        end 
        
        if abs(Cross2D(r,s)) < sqrt(eps) && ~(all([t u] >= 0) && all([t u] <= 1))
            % Co-linear?
            
            %keyboard
            
        end
        
    end
    
    if isempty(NewSharedVertIdx) || EdgePtCollision == 1 || ...
            StartPointOnEdge == 1
        % Adjust point
        Adj = sqrt(eps);
        if EdgePtCollision == 1 || StartPointOnEdge == 1
            PtTemp_ref = IntersecPt(1:2);
            RefPt_a = IntersecPt(1:2);
            RefPt_b = GeoLine(Ctr-1,1:2);
        else
            PtTemp_ref = GeoLine(Ctr-1,1:2);
            RefPt_a = GeoLine(Ctr-1,1:2);
            RefPt_b = GeoLine(Ctr-2,1:2);
        end
        
        ForwardVec = (GeoLine(Ctr-1,1:2) - GeoLine(Ctr-2,1:2));
        ForwardVec = ForwardVec/norm(ForwardVec);
        
        if EdgePtCollision == 1 || StartPointOnEdge == 1
            PerpVec = [0 0]; 
        else
            PerpVec = [-ForwardVec(2) ForwardVec(1)];
        end
        
        % Maybe start with just moving forward and if that does not
        % help also move perpendicular?
        NewPtTemp = PtTemp_ref + ForwardVec*Adj + PerpVec*Adj;
        
        NewTriID = pointLocation(DT,NewPtTemp);
        
        if isempty(NewTriID) || isnan(NewTriID)    
            MissingLength = L_target - L_path;
            if Set.MeshOrCurveWarning
                fprintf(2,'Geodesic line reached the end of the mesh (generator #%d)\n\n',i)
            end
            break
        end
        
        % Check if new triangle is a neighbor with previous
        NeighborID = neighbors(DT,NewTriID);
        
        if ismember(TriID,NeighborID)
            
            VertIdx_new = DT.ConnectivityList(NewTriID,:);
            
            NewSharedVertIdx = intersect(VertIdx,VertIdx_new);
      
        else
            % Check which two vertices have the shortest projected
            % distance to old intersection point
            VertIdx_new = DT.ConnectivityList(NewTriID,:);
            proj_dist(1:3) = 0.0;
            for k = 1:3
                TestVertex = DT.Points(VertIdx_new(k),:);
                a_vec = TestVertex - RefPt_a;
                b_vec = NewPtTemp - RefPt_b;
                
                proj = b_vec * dot(a_vec,b_vec)/dot(b_vec,b_vec);
                proj_dist(k) = norm(proj);
            end
            
            [~,TwoMinIdx] = mink(proj_dist,2);
            
            NewSharedVertIdx = VertIdx_new(TwoMinIdx);
         end

        TestEdgePt(1:2,1:2) = DT.Points(NewSharedVertIdx,:);
        TestEdgePt(1:2,3) = F_DT(DT,Z,TestEdgePt(1:2,1),TestEdgePt(1:2,2));
        NewSharedEdgeVec = TestEdgePt(2,1:3) - TestEdgePt(1,1:3);
        
        IntersecPt(1:2) = NewPtTemp;
        IntersecPt(3) = F_DT(DT,Z,IntersecPt(1),IntersecPt(2));
        
        % Adjust previous point on geodesic path to ensure the direction is
        % calculated correctly
        GeoLine(Ctr-1,:) = GeoLine(Ctr-1,:) + [PerpVec*Adj 0];
        
        FoldAng = 0;
        
        if IncPlt
            scatter3(IntersecPt(1),IntersecPt(2),IntersecPt(3),300,'bx')
        end
        
        % Print warning message
        if StartPointOnEdge == 1 && Set.MeshOrCurveWarning
            fprintf('Corrected geodesic path due to starting point being on an edge\n\n')
        elseif EdgePtCollision == 1 && Set.MeshOrCurveWarning
            fprintf('Corrected geodesic path due collision with a vertex \n\n')
        elseif Set.MeshOrCurveWarning
            fprintf('Corrected geodesic path. Unable to identify neighbor triangle\n\n')
        end
        
    end
    
    % Rotate the intersection point back (fold) with negative angle
    VecToPt = IntersecPt - GeoLine(Ctr-1,:);
    RotatedVecBack = VecToPt*cosd(-FoldAng) + cross(RotAxis,VecToPt)*...
        sind(-FoldAng) + RotAxis*dot(RotAxis,VecToPt)*(1-cosd(-FoldAng));
    
    FoldedIntersecPt = RotatedVecBack + GeoLine(Ctr-1,:);
    % Store as the current point on the geodesic line
    GeoLine(Ctr,:) = FoldedIntersecPt;
    
    % Update path length
    L_path = L_path + norm(GeoLine(Ctr,:)-GeoLine(Ctr-1,:));
       
    if IncPlt
        plot3(UnfoldVertCoor([1 2 3 1],1),UnfoldVertCoor([1 2 3 1],2) , ...
            UnfoldVertCoor([1 2 3 1],3),'m--')
        scatter3(IntersecPt(1), IntersecPt(2), IntersecPt(3),'rd','filled')
        scatter3(FoldedIntersecPt(1), FoldedIntersecPt(2), ...
            FoldedIntersecPt(3),'gs','filled')
        trimesh(DT.ConnectivityList(TriID,:),DT.Points(:,1),...
            DT.Points(:,2),Z(:),'edgecolor','k')
        %quiver3(Centroid(1),Centroid(2),Centroid(3),NormVec(1),...
        %    NormVec(2),NormVec(3),5)
        plot3(GeoLine(Ctr-2 : Ctr-1,1),GeoLine(Ctr-2 : Ctr-1,2),...
            GeoLine(Ctr-2 : Ctr-1,3),'g--')
        axis('tight')
        %keyboard
    end
    
    % Check if current point is the same as the previous
    if all(GeoLine(Ctr,:)==GeoLine(Ctr-1,:))
        keyboard
    end
   
    % Check if the angle to the new point is not too small (-> error)
    Vec1 = GeoLine(Ctr,:) - GeoLine(Ctr-1,:);
    Vec2 = GeoLine(Ctr-2,:) - GeoLine(Ctr-1,:);
    Vec1(3) = 0;
    Vec2(3) = 0;
    GeoLineAng = atan2d(vecnorm(cross(Vec1,Vec2),2,2),dot(Vec1,Vec2,2));
    GeoAngTol = 100;
    if GeoLineAng < GeoAngTol
        GeoLine(Ctr,:) = [];
        MissingLength = L_target - L_path;
        if Set.MeshOrCurveWarning
            fprintf(2,['Geodesic line algorithm encountered angle in XY-plane\n' ...
                'between two points on the line of less than %.4g deg\n' ...
                'with %.3g %% remaining. Aborting... (generator #%d)\n\n'],...
                GeoAngTol,100*MissingLength/L_target,i)
        end
        %keyboard
        break
    end

    % Check if no new triangle is found
    if isempty(NewTriID)
        MissingLength = L_target - L_path;
        if Set.MeshOrCurveWarning
            fprintf(2,'Geodesic line reached the end of the mesh (generator #%d)\n\n',i)
        end
        %keyboard
        break
    end
end
end

%% Aux. functions
function result = Cross2D(a,b)
if length(a) > 2 || length(b) > 2
    error('Vectors must have length 2')
end
result = a(1)*b(2) - a(2)*b(1);
end

function NormC = FourPtCoplanarCheck(P1,P2,P3,P4)
%% Method 1: cross product
C = cross (cross(P2-P1, P3-P1), cross(P2-P1, P4-P1));

NormC = norm(C);
end

function Z = MY_setdiff(X,Y)
check = false(1, max(max(X), max(Y)));
check(X) = true;
check(Y) = false;
Z = X(check(X));
end