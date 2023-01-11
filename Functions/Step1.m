function Dra = Step1(d,Grid,Ang,Org,OrgNode,PreShear,F,DT,z,DraContr,...
    SteerPts,SteerPtsRef,Mold,Plt,Set)
% Step 1: Generators (geodesic from unfolded tri or steering curve)

% Define propagation directions in the Node array
Dir1 = [1 0 -1 0 ; 0 1 0 -1]';
Dir2 = [0 -1 0 1; 1 0 -1 0]';

% Initialize Node
Node = NaN([Grid 3]);

TargetY = (OrgNode(2)-1) * d;% + BottomCourseOffset;

% Convert the steering points to absolute values
[SteerPtsAbs,OrgAdj] = AbsSteeringPts(SteerPts,DraContr,SteerPtsRef,...
    Mold,Org,F,d,TargetY);

% First place origin node
Idx0 = LinearCellIdx(Grid,OrgNode(1),OrgNode(2),Dir1,Dir2,1);
Node(Idx0(1,:)) = [OrgAdj(1), OrgAdj(2), F(OrgAdj(1), OrgAdj(2))];

% Calculate the number of nodes in each generator
nIniNode = [Grid-OrgNode  OrgNode-1];

% For geodesic lines: find transverse angles in XY-plane that produces 
% zero shear on surface at the origin.
% The longitudinal angle is determined from the steering curve or 'Ang'
% from the input file in the case of a geodesic longitudinal generator
if strcmpi(DraContr{1},'Geo') && strcmpi(DraContr{2},'Geo')
    LongAng = Ang+90;
elseif strcmpi(DraContr{1},'Steer') && strcmpi(DraContr{2},'Geo')
    % Place a single node on the steering curve 'd' away from the origin
    LongNode1 = PutNodesOnCurve(OrgAdj,SteerPtsAbs{1},F,d,1,2,Plt,Set);
    LongNode1(3) = [];
    % Find the angle of that node
    LongVec = LongNode1 - OrgAdj;
    LongAng = atan2d(1*LongVec(2) - 0*LongVec(1), ...
        1*LongVec(1) + 0*LongVec(2));   
end

if any(strcmpi(DraContr{1},'Geo') | strcmpi(DraContr{2},'Geo'))
    % Find up to two transverse angles that produces zero shear on surface
    TransAng(1:2) = NaN;
    for i = find(nIniNode([1 3])>0)
        TransAng(i) = AngOnSurf(LongAng,Node(Idx0(1,:)),d,F,PreShear,i,Set);
    end
    AngVec = [TransAng(1) LongAng TransAng(2) LongAng+180];
end

% Loop over each generator (1st: transverse, 2nd: longitudinal, 3rd:
% transverse, 4th: longitudinal)
GenCurve = cell(1,4);
for i = find(nIniNode>0)
    
    % Reset the generator nodes created
    GenNode = [];
    
    % Use geodesic line or steering curve depending on settings
    if ((i == 1 || i == 3) && strcmpi(DraContr{2},'Geo')) || ...
            ((i == 2 || i == 4) && strcmpi(DraContr{1},'Geo'))
        % Geodesic    
        [GenCurve{i}, MissingLength] = GeoUnfoldTri(OrgAdj,AngVec(i),d,...
            nIniNode(i),DT,z,F,Plt.GeoInc,Set,i);

        if isnan(GenCurve{i}(1,3))
            GenCurve{i}(1,3) = F(GenCurve{i}(1,1),GenCurve{i}(1,2));
        end
        
        % Extrapolate geodesic line if the end of the mesh is reached
        if MissingLength > 0 && ~isempty(GenCurve{i})
            % Correct the missing length for the slope (using two points)
            x1 = GenCurve{i}(end-1,1);
            x2 = GenCurve{i}(end,1);
            y1 = GenCurve{i}(end-1,2);
            y2 = GenCurve{i}(end,2);
            % Use the interpolating object to get z in case the GeoLine 
            % algorithm crashed
            z1 = F(x1,y1);
            z2 = F(x2,y2);
            Angle = abs(atan2d(z2-z1,norm([x2-x1 y2-y1])));
            MissingLength_corr = MissingLength*cosd(Angle);
            % Extrapolate linearly upon early return from GeoUnfoldTri
            if i == 1 || 3
                % Use x as independent variable
                ExtrapX = GenCurve{i}(end,1) + ...
                    linspace(0.1,MissingLength_corr*2,5)' * ...
                    sign(GenCurve{i}(end,1)-GenCurve{i}(end-1,1));
                ExtrapY = interp1(GenCurve{i}(:,1),GenCurve{i}(:,2),...
                ExtrapX,'linear','extrap');
                ExtrapZ = F(ExtrapX,ExtrapY);
            elseif i == 2 || 4
                % Use y as independent variable
                ExtrapY = GenCurve{i}(end,2) + ...
                    linspace(0.1,MissingLength_corr*2,5)' * ...
                    sign(GenCurve{i}(end,2)-GenCurve{i}(end-1,2));
                ExtrapX = interp1(GenCurve{i}(:,2),GenCurve{i}(:,1),...
                    ExtrapX,'linear','extrap');
                ExtrapZ = F(ExtrapX,ExtrapY);
            end
            % Append to the previously computed points
            GenCurve{i} = [GenCurve{i} ; [ExtrapX ExtrapY ExtrapZ]];

            if Set.MeshOrCurveWarning
                fprintf(2,['Extrapolated geo. line for generator #%d by'...
                    ' %.3g × d \n\n'],i,MissingLength/d);
            end

        end

        [GenNode, ~] = PutNodesOnCurve(OrgAdj,...
            GenCurve{i}(:,1:2)',F,d,nIniNode(i),i,Plt,Set);
        
    elseif (i == 2 && strcmpi(DraContr{1},'Steer')) || ...
            (i == 1 && strcmpi(DraContr{2},'Steer') && nIniNode(i) > 0) || ...
            (i == 3 && strcmpi(DraContr{2},'Steer') && nIniNode(i) > 0) || ...
            (i == 4 && strcmpi(DraContr{1},'Steer') && nIniNode(i) > 0)
        % Compute which steering points to use
        if i == 1
            SteerPts_curr = SteerPtsAbs{2};
        elseif i == 3
            % Flip the transverse steerings points if they are not
            % already in descending x coordinates
            if issorted(SteerPtsAbs{2}(1,:),'descend')
                SteerPts_curr = SteerPtsAbs{2};
            else
                SteerPts_curr = fliplr(SteerPtsAbs{2});
            end
        elseif i == 2
            SteerPts_curr = SteerPtsAbs{1};
        elseif i == 4
            % Flip the longitudinal steerings points if they are not 
            % already in descending y coordinates
            if issorted(SteerPtsAbs{1}(2,:),'descend')
                SteerPts_curr = SteerPtsAbs{1};
            else
                SteerPts_curr = fliplr(SteerPtsAbs{1});
            end
        end
        % Steering curve
        [GenNode, GenCurve{i}] = ...
            PutNodesOnCurve(OrgAdj,SteerPts_curr,F,d,nIniNode(i),i,Plt,Set);
    end
    
    % Check distances between nodes
    NodeDistDiff = vecnorm(diff([Node(Idx0(1,:)) ; GenNode],1),2,2) - d;
    if any(abs(NodeDistDiff) > 5e-5)
        keyboard
    end
    
    % Store generator nodes in Node array
    NodePos = OrgNode' + (1:nIniNode(i)).*Dir1(i,:)';
    IdxGen = LinearGenIdx(Grid,NodePos);
    Node(IdxGen(:)) = GenNode(:);
    
    % If necessary, adjust the nIniNode array if all nodes could not be
    % placed due to e.g. reaching the end of the mold mesh
    if any(isnan(GenNode))
        
        [NotNaNRows,~] = find(~isnan(GenNode(:,1)));
        
        nIniNode(i) = max(NotNaNRows);
        
        if Set.MeshOrCurveWarning
            fprintf(2,['Cropped generator #%d to fit the mold mesh'...
                '/steering curve \n\n'],i);
        end
    end
end

Dra.Node = Node;
Dra.GenCurve = GenCurve;
Dra.nIniNode = nIniNode;
Dra.SteerPtsAbs = SteerPtsAbs;
Dra.OrgAdj = OrgAdj;
end

function XYAng = AngOnSurf(InpAng,Node1,d,F,PreShear,TransAngNo,Set)

% Make the step smaller in case it could extend outside the boundary?
if isfield(Set,'ForceGeoLineIniAngle') && Set.ForceGeoLineIniAngle
    d_mod = d;
else
    d_mod = d/100;
end

Node2(1:2) = Node1(1:2) + d_mod*[cosd(InpAng) sind(InpAng)];
Node2(3) = F(Node2(1),Node2(2));

% Find angle using bisection
% Set search interval based on transverse angle no, i.e. generator 1 or 3
if TransAngNo == 1
    Theta = InpAng + [-180 0];
    ObjSignControl = 1;
elseif TransAngNo == 2
    Theta = InpAng + [1 179];
    ObjSignControl = -1;
end

XYAng = NaN;
for i = 1:1e3
    % Compute middle point
    Theta_mid = mean(Theta);
    FunVal_mid = AngFun(Theta_mid,Node1,Node2,d_mod,F) - PreShear;
    FunVal_mid = FunVal_mid * ObjSignControl;
    % Stop or adjust interval based on function value at midpoint
    if abs(FunVal_mid) < 1e-5
        XYAng = Theta_mid;
        break
    elseif FunVal_mid > 0
        Theta(2) = Theta_mid;
    else
        Theta(1) = Theta_mid;
    end
end

if isnan(XYAng)
    fprintf(2,'Could not determine angle on surface\n\n')
    keyboard
    return
end

    function AngDiff = AngFun(x,Node1,Node2,d_mod,F)
        % Compute vector to new node based on input angle (1st iteration)
        Node4(1:2) = Node1(1:2) + d_mod*[cosd(x) sind(x)];
        Node4(3) = F(Node4(1),Node4(2));
        
        % Calculate the mold slope angle and correct new node
        SlopeAng = atand((Node4(3)-Node1(3))/norm(Node4(1:2)-Node1(1:2)));
        Node4(1:2) = Node1(1:2) + d_mod*cosd(SlopeAng)*[cosd(x) sind(x)];
        Node4(3) = F(Node4(1),Node4(2));

        % Set up edge vectors in cell
        u = (Node2 - Node1)';
        v = (Node4 - Node1)';

        % Calculate the angle difference
        AngDiff = 90 - atan2d(vecnorm(cross(u,v),2,1),dot(u,v));
    end
end

function Idx = LinearGenIdx(Grid,RowCol)
Idx = RowCol(1,:)' + (RowCol(2,:)'-1)*Grid(1) + ((1:3)-1)*Grid(1)*Grid(2);
end