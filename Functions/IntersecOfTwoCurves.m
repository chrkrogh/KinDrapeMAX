function [IntersecPt,MinIdxC1,MinIdxC2] = ...
    IntersecOfTwoCurves(Curve1,Curve2,Cycl1,Cycl2)
% This function finds the intersection of two curves in 2D defined by
% points. Curve1 and Curve2 must be of dimension nPoints x 2
% Cycl1 and Cycl2 control if the curves have cyclic behavior, e.g. if the
% last index plus 1 equals the first index (true) or it means out of bounds
% (false).

if nargin == 2
    Cycl1 = false;
    Cycl2 = false;
end

IntersecPt = [NaN NaN];

% Find the two closest points on the two curves
% Create 3D array with difference between all points. The rows and cols
% are respectively curve 1 and curve 2 points and the
% two pages have x,y coordinates
if size(Curve1,1) == 2
   % When curve 1 has two points, e.g. a a cell edge, use the mean point to
   % add robustness
   %MeanC1 = mean(Curve1,1);
   MeanC1 = sum(Curve1,1)/size(Curve1,1);
   
   PtDiff = cat(3,MeanC1(:,1)-Curve2(:,1)',...
        MeanC1(:,2)-Curve2(:,2)');
else
   PtDiff = cat(3,Curve1(:,1)-Curve2(:,1)',...
        Curve1(:,2)-Curve2(:,2)');
end
% Calculate point distances as the 2-norm along the 3rd dimension
PtDist = vecnorm(PtDiff,2,3);
% Find the minimum point distance
[~,MinIdx] = min(PtDist(:));
[MinIdxC1,MinIdxC2] = ind2sub(size(PtDist),MinIdx);

% Test intersections in a neighborhood of the closest point of the two 
% curves
Neighborhood = [-2 -1 0 1];

Ctr = 0;
nPts = length(Neighborhood)^2;
IntersecPt_temp(1:nPts,1:2) = NaN;
t(1:nPts) = NaN;
u(1:nPts) = NaN;

for i = Neighborhood
    for j = Neighborhood
        
        Ctr = Ctr + 1;
        
        NewIdx1 = TestNewIdx(Curve1,MinIdxC1,i,Cycl1);
        AdjIdx1 = TestNewIdx(Curve1,MinIdxC1,i+1,Cycl1);
        
        NewIdx2 = TestNewIdx(Curve2,MinIdxC2,j,Cycl2);
        AdjIdx2 = TestNewIdx(Curve2,MinIdxC2,j+1,Cycl2);
        
        if isempty(NewIdx1) || isempty(AdjIdx1) || isempty(NewIdx2) ...
                || isempty(AdjIdx2)
            % Skip the iteration
            continue
        end
        
        Pt11 = Curve1(NewIdx1,:);
        Pt12 = Curve1(AdjIdx1,:);
        Pt21 = Curve2(NewIdx2,:);
        Pt22 = Curve2(AdjIdx2,:);
        
        % Compute the intersection point of the lines defined by Pt11,Pt12,
        % Pt21, and Pt22
        [IntersecPt_temp(Ctr,:),t(Ctr),u(Ctr)] = ...
            IntersecOfTwo2DLines(Pt11,Pt12,Pt21,Pt22);
        
        % If 0<= t,u <= 1, then the lines intersect inside the interval
        % defined by the points
        if all([t(Ctr),u(Ctr)] >= 0) && all([t(Ctr),u(Ctr)] <= 1)
            IntersecPt = IntersecPt_temp(Ctr,:);
            break
        end
    end
end

if all(isnan(IntersecPt))
    % If no intersections inside the interval are found, then take the
    % intersection point with the smallest distance outside the interval
    [~,BestIdx] = min(sqrt(t.^2 + u.^2),[],'omitnan');
    
    IntersecPt = IntersecPt_temp(BestIdx,:);
end

end

function NewIdx = TestNewIdx(Curve,BaseIdx,Dir,Cycl)
% Test if a new index is out of bounds (Cycl = false) and set it to empty
% or adjust it in a cyclic behavior (Cycl = true)

NewIdx = BaseIdx + Dir;
if NewIdx < 1
    if Cycl
        % Count backwards from last index
        NewIdx = size(Curve,1) + NewIdx;
    else
        % Skip this iteration
        NewIdx = [];
    end
elseif NewIdx > size(Curve,1)
    if Cycl
        % Count from first index
        NewIdx = NewIdx - size(Curve,1);
    else
        % Skip this iteration
        NewIdx = [];
    end
end
end