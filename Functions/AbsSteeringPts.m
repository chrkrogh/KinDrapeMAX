function [SteerPtsAbs,OrgAdj] = AbsSteeringPts(SteerPts,DraContr,...
    SteerPtsRef,Mold,Org,F,d,TargetY)
% This function returns the steering points in absolute, global coordinates
% based on the definition relative to a mold edge (preform division) or
% another curve, e.g. an edge of a previously draped course.
%
% The relative steering pts are defined in xy-coordinates and for the
% longitudinal points a mold slope compensation is carried out using a 
% secant angle.
%
% If the steering points are already defined in the absolute, global
% coordinates, they are just passed on.
% For converting the steering points, the idea is to create a smooth spline
% fit based on the offset distances, that can be added to the reference
% curve.
%
% The function also computes an adjusted origin:
% In the case where two steering curves are defined as generators, the
% origin is adjusted to be their intersection.
% In the case where a longitudinal steering curve and a transverse geodesic
% line is specified as the generators, the adjusted origin is calculated
% based on the y-coordinate (2nd component) of Org

% Initialize
OrgAdj = Org;
SteerPtsAbs = {[],[]};

% Longitudinal direction
if strcmpi(DraContr{1},'Steer')
    if strcmpi(SteerPtsRef{1},'Abs')
        % The steerpts are already in absolute coordinates
        SteerPtsAbs{1} = SteerPts{1};
    elseif any(strcmpi(SteerPtsRef{1},...
            {'Left','Right','Left-CourseEdge','Right-CourseEdge'}))
        % Find the right reference curve to use
        if strcmpi(SteerPtsRef{1},'Left')
            No = 1;
            RefCurve = Mold.MoldEdge{No};
        elseif strcmpi(SteerPtsRef{1},'Right')
            No = 2;
            RefCurve = Mold.MoldEdge{No};
        elseif strcmpi(SteerPtsRef{1},'Left-CourseEdge')
            No = 1;
            RefCurve = Mold.RightEdge;
        elseif strcmpi(SteerPtsRef{1},'Right-CourseEdge')
            No = 2;
            RefCurve = Mold.LeftEdge;
        end
        
        % First correct the input SteerPts such that the mold arc length is
        % taken into account. That is, the points are defined in the
        % XY-plane but the offset distance must be along the mold surface.
        % With small offsets (less than d), the slope should be 
        % approximately constant. Otherwise, use the true arc length.

        % Locate (project) the relative steering points on the reference 
        % curve and store in the variable Pt1. The y-coordinates are the 
        % same as the y-coordinates of the input relative steer pts and the 
        % x-coordinates are found from linear interpolation.
        Pt1(1:size(SteerPts{1},2),2) = SteerPts{1}(2,:);
        % Make a linear fit of the reference curve and use that to locate
        % the x-coordinates of Pt1.
        % First make sure that all the points are in ascending order (this
        % check was introduced due to an incident where a course used as
        % the ref. was curving a lot)
        if ~issorted(RefCurve(:,2))
            % Find the points that are not ascending and remove until end
            RmIdx = find(diff(RefCurve(:,2)) < 0) + 1;
            RefCurve(RmIdx:end,:) = [];
        end
        RefCurveFit = griddedInterpolant(RefCurve(:,2),RefCurve(:,1),...
            'linear','linear');
        Pt1(:,1) = RefCurveFit(SteerPts{1}(2,:));
        % If the offset distance is less than d, the offset distance in the
        % XY-plane can be found from cosine to the secant angle
        if max(abs(SteerPts{1}(1,:))) <= d
            % Query the z-coordinates of Pt1
            Pt1(:,3) = F(Pt1(:,1),Pt1(:,2));
            % Create a 2nd point set: Pt1 plus the offset in the
            % x-direction
            Pt2 = Pt1;
            Pt2(:,1) = Pt2(:,1) + SteerPts{1}(1,:)';
            % Again, query the z-coordinates
            Pt2(:,3) = F(Pt2(:,1),Pt2(:,2));
            % Calculate the secant angle of the mold between Pt1 and Pt2
            Angle = abs(atan2d(abs(Pt2(:,3)-Pt1(:,3)),...
                abs(Pt2(:,1)-Pt1(:,1))));
            % Adjust the input relative steer pts such that their offset
            % distance in the xy-plane equals the desired distance along 
            % the mold surface
            SteerPts{1}(1,:) = SteerPts{1}(1,:) .* cosd(Angle)';
        else
            % If the offset distance is larger than d, use the mold arc
            % length. Create a 100 pt vector emanating from each rel. 
            % steer pt on the ref. curve with a length equal to the offset 
            % distance. Find the two indices where the distance is closest 
            % to the desired offset distance and interpolate linearly in
            % between.
            SteerPts_upd(1,1:size(SteerPts{1},2)) = 0;
            for ONo = 1:size(SteerPts{1},2)
                % The current rel. steer pt
                TargetSteerPt = SteerPts{1}(:,ONo);
                % Create the vector
                TrialVecX = Pt1(ONo,1) + linspace(0,TargetSteerPt(1),100);
                TrialVecY = TargetSteerPt(2) * ones(1,100);
                TrialVecZ = F(TrialVecX,TrialVecY);
                % Compute the cumulative distance of the vector
                CumDist = [0 cumsum(sqrt(diff(TrialVecX).^2 + ...
                    diff(TrialVecY).^2 + ...
                    diff(TrialVecZ).^2))];
                % Find the two distances clostst to the target offset dist.
                [~,idx] = mink(abs(CumDist-abs(TargetSteerPt(1))),2);
                Idx1 = min(idx);
                Idx2 = max(idx);
                % Find the weighting between the two indices
                Frac = (abs(TargetSteerPt(1))-CumDist(Idx1))/...
                    (CumDist(Idx2)-CumDist(Idx1));
                if isnan(Frac)
                    % Most likely due to division with 0, i.e. the values
                    % at the two idx are identical
                    Frac = 1;
                end
                % Use the weighting and corresponding indices to get the
                % interpolated x-coordinate
                SteerPts_upd(1,ONo) = Frac*TrialVecX(max(idx)) + ...
                    (1-Frac)*TrialVecX(min(idx)) - Pt1(ONo,1);
            end
            SteerPts{1}(1,:) = SteerPts_upd;  
        end

        % Now use the (mold slope compensated) relative steerpts and the
        % reference curve to calculate the absolute steer pts
        SteerPtsAbs{1} = OffsetRefPts(SteerPts{1},RefCurve,No);
    end
    
    % Adjust the origin if using a longitudinal steering curve in
    % combination with a transverse geodesic line
    if strcmpi(DraContr{2},'Geo')
        
        % Calculate the course offsetting from the bottom, i.e. root
        % Only for Set.Mode = 'analysis' because it is handled differently
        % with optimization (see the function LayerObj).
        if ~isempty(TargetY)
            % Auto-move origin so that the course starts at TargetY
            % Get z coordinate of AbsSteerPts
            AbsSteerPtsZ = F(SteerPtsAbs{1}(1,:),SteerPtsAbs{1}(2,:));
            % Calculate cumulative distance
            CumDist = [0 cumsum(sqrt(diff(SteerPtsAbs{1}(1,:)).^2 + ...
                diff(SteerPtsAbs{1}(2,:)).^2 + diff(AbsSteerPtsZ).^2))] + ...
                SteerPtsAbs{1}(2,1);
            % Find the two points clostst to the target y
            [~,idx] = mink(abs(CumDist-TargetY),2);
            Idx1 = min(idx);
            Idx2 = max(idx);
            Frac2 = (TargetY-CumDist(Idx1))/(CumDist(Idx2)-CumDist(Idx1));

            Org(2) = Frac2*SteerPtsAbs{1}(2,max(idx)) + ...
                (1-Frac2)*SteerPtsAbs{1}(2,min(idx));
        end

        % Use the y-coordinate (2nd component) of Org as the target
        OrgAdj = AdjOrgSteerPts(SteerPtsAbs,Org);
    end
end

% Transverse direction
if strcmpi(DraContr{2},'Steer')
    if strcmpi(SteerPtsRef{2},'Bottom')
        No = 3;
        RefCurve = Mold.MoldEdge{No};
        SteerPtsAbs{2} = OffsetRefPts(SteerPts{2},RefCurve,No);
    elseif strcmpi(SteerPtsRef{2},'Abs')
        SteerPtsAbs{2} = SteerPts{2};
    end
end

% Adjust origin to be intersection if two steering curves are used
if strcmpi(DraContr{1},'Steer') && strcmpi(DraContr{2},'Steer')
    
    % Find the intersection point and the index of the point on each curve,
    % that is closest to the intersection point
    [IntersecPt,LongIdx,TransIdx]...
        = IntersecOfTwoCurves(SteerPtsAbs{1}',SteerPtsAbs{2}',0,0);
    
    % Take the adjusted origin point as the intersection point
    OrgAdj = IntersecPt;
    
    % Adjust the closest points of the two steering curves to be the
    % intersection point
    SteerPtsAbs{1}(:,LongIdx) = IntersecPt;
    SteerPtsAbs{2}(:,TransIdx) = IntersecPt;
end
end

function AbsSteerPts = OffsetRefPts(SteerPtsSel,RefCurve,No)
% Calculate points to add to the right, left or bottom edge points
% No is 1, 2, or 3 i.e. left, right or bottom

% Left and right edges: offset x coordines. Bottom edge: offset y coord.
% This is organized with Coor1 and Coor2
if No == 1 || No == 2
    % Longitudinal
    Coor1 = 2;
    Coor2 = 1;
else
    % Transverse
    Coor1 = 1;
    Coor2 = 2;
end

if size(SteerPtsSel,2) == 1
    % Add a constant offset
    %PtsToAdd = [SteerPtsSel(1,1) 0]';
    PtsToAdd(1:2,1) = 0;
    PtsToAdd(Coor2,1) = SteerPtsSel(Coor2,1);
else
    Method = 'spline';
    if strcmpi(Method,'linear')
        % Loop over all steering points
        nPtsToAdd = size(RefCurve,1);
        PtsToAdd = NaN(2,nPtsToAdd);
        MinIdx_m1 = 1;
        % Linear interpolation
        for k = 2:size(SteerPtsSel,2)
            % Longitudinal: find nearest y coordinate in edge points
            % Transverse: find nearest x coordinate in edge points
            [~,MinIdx] = min(abs(SteerPtsSel(Coor1,k)-RefCurve(:,Coor1)));
            % Longitudinal: Interpolate linearly between the points in x
            % Transverse: Interpolate linearly between the points in y
            PtsToAdd(Coor2,MinIdx_m1:MinIdx) = linspace(...
                SteerPtsSel(Coor2,k-1),SteerPtsSel(Coor2,k),...
                MinIdx-MinIdx_m1+1);
            % Longitudinal: Set the y coordinates equal to zero
            % Transverse: Set the x coordinates equal to zero
            PtsToAdd(Coor1,MinIdx_m1:MinIdx) = 0;
            % Store previous index
            MinIdx_m1 = MinIdx;
        end
        % Insert the last value in the remaining (possible) empty
        % entries
        PtsToAdd(:,MinIdx:end) = PtsToAdd(:,MinIdx).*...
            ones(2,nPtsToAdd-MinIdx+1);
    elseif strcmpi(Method,'spline')
        % Spline interpolation
        PtsToAdd = zeros(2,size(RefCurve,1));
        % Notice that the x and y coordinates are interchanged such that for
        % each x coordinate in the spline fit, there will be a unique y
        % coordinate
        % The query points
        x_query = RefCurve(:,Coor1)';
        % Create a spline fit based on the input (relative) steering points
        % and evaluate it at the reference points (griddedInterpolant is
        % faster than e.g. spline)
        Spl = griddedInterpolant(SteerPtsSel(Coor1,:),...
            SteerPtsSel(Coor2,:),'spline','spline');
        PtsToAdd(Coor2,:) = Spl(x_query);
    end
end
% Add points to the ref pts.
AbsSteerPts = RefCurve(:,1:2)' + PtsToAdd;

% Make sure all rel. steer pts are covered, i.e. extend AbsSteerPts
if (No == 1 || No == 2) && AbsSteerPts(Coor1,end) < SteerPtsSel(Coor1,end)
    % Create vector between the end of AbsSteerPts and the input (relative)
    % steering pts
    x_query_new = linspace(AbsSteerPts(Coor1,end),SteerPtsSel(Coor1,end),10);
    if strcmpi(Method,'linear')
        % Linear / constant extrapolation
        AbsSteerPts = [AbsSteerPts  [ones(1,10)*AbsSteerPts(Coor2,end); ...
            x_query_new]];
    elseif strcmpi(Method,'spline')
        % Spline extrapolation
        % Extrapolate the reference points until the end of the defined
        % (relative) steering points
        ExtrapFit = griddedInterpolant(RefCurve(:,Coor1)',...
            RefCurve(:,Coor2)','pchip','pchip');
        ExtrapRefPts = ExtrapFit(x_query_new);
        % Append the AbsSteerPts using the extrapolated points and the
        % previous spline fit
        AbsSteerPts = [AbsSteerPts(:,1:end-1) ...
            [ExtrapRefPts + Spl(x_query_new) ; x_query_new] ];
    end
end
end

function OrgAdj = AdjOrgSteerPts(SteerPtsAbs,Org)
% Adjust the origin when using a single steering curve and a geodesic line
% as the generators. The adjusted origin point must be located on the curve
% of absolute longitudinal steering points.
% The target y coordinate of the adjusted origin point is the y-coordinate
% (2nd component) of Org.

% Calculate the target y coordinate
TargetYCoor = Org(2);
% Find the index of the abs. steer pts. y coord. closest to the target y
[~,MinIdx] = min(abs(SteerPtsAbs{1}(2,:)-TargetYCoor));
% Interpolate linearly to determine the x coordinate on the absolute
% steering curve that matches the input y coordinate
NoOfPts = size(SteerPtsAbs{1}(2,:),2);
if SteerPtsAbs{1}(2,MinIdx) < TargetYCoor
    Vec = SteerPtsAbs{1}(:,min(MinIdx+1,NoOfPts)) - ...
        SteerPtsAbs{1}(:,MinIdx);
    RemainingYDist = TargetYCoor - SteerPtsAbs{1}(2,MinIdx);
else
    Vec = SteerPtsAbs{1}(:,max(1,MinIdx-1)) - ...
        SteerPtsAbs{1}(:,MinIdx);
    RemainingYDist = SteerPtsAbs{1}(2,MinIdx) - TargetYCoor;
end
% Create unit vector
if ~all(Vec==0)
    Vec = Vec / norm(Vec);
end
VecScale = abs(RemainingYDist);
OrgAdj = SteerPtsAbs{1}(:,MinIdx)' + Vec'*VecScale;
end