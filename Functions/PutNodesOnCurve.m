function [CurveNode, GenCurve] = PutNodesOnCurve(Org,CurvePts,F,d,...
    nNode,i,Plt,Set)
% This function takes a curve defined by points (CurvePts) and a starting
% point, Org, and populates it with nNode nodes spaced with a distance d. F
% is the interpolated surface object and i is the generator number used for
% determing the right progression direction

% Max number of iterations and convergence tolerance
MaxIter = 1e3;
ConTol = 1e-5;

% Coordinate for progression measurement
if i == 1 || i == 3
    PrgCoor = 1;
elseif i == 2 || i == 4
    PrgCoor = 2;
end

% Sign due to direction
if i == 1 || i == 2
    Dir = 1;
elseif i == 3 || i == 4
    Dir = -1;
end

% Create parametric cubic spline fit
if any(isnan(CurvePts(:)))
    CurvePts(:,isnan(CurvePts(1,:))) = [];
end
SplnFt = cscvn(CurvePts);

% Evaluate curve for plotting
if Plt.GenCurve
    GenCurve = fnval(SplnFt,linspace(SplnFt.breaks(1),...
        SplnFt.breaks(end),50))';
    GenCurve(:,3) = F(GenCurve(:,1),GenCurve(:,2));
else
    GenCurve = [];
end
% Retrieve breaks and coefficients from the spline
SplBrk = SplnFt.breaks.';
cf_x = SplnFt.coefs(1:2:end,:);
cf_y = SplnFt.coefs(2:2:end,:);    

% Populate spline curve with nodes using bisection
% Initialize and place first node
CurveNode(1:nNode+1,1:3) = NaN;
CurveNode(1,1:3) = [Org, F(Org(1),Org(2))]; 

% Locate the last spline break point before the first curve point
BrkIdx = find(Dir.*sign(CurvePts(PrgCoor,:)-CurveNode(1,PrgCoor))...
    <=0,1,'last');
if isempty(BrkIdx)
    BrkIdx = 1;
end
%[~,min_idx] = min(sqrt((CurvePts(1,:)-CurveNode(1,1)).^2 + ...
%    (CurvePts(2,:)-CurveNode(1,2)).^2));
L_prev = SplBrk(BrkIdx);

% Loop over number of nodes
EndOfCurveFlag = false;
for NodeIter = 2:nNode + 1
    % Find upper interval
    L_upper = L_prev;
    % Define the base point
    BasePt = CurveNode(NodeIter-1,:);
    % Locate an upper point in an interval that brackets the next node
    for j = 1:MaxIter
        % Go a little further on the curve
        % From MATLAB cscvn doc: the parameter value of the spline is 
        % accumulated square root of chord length
        L_upper = L_upper + sqrt(d);

        % Evaluate the spline at L_upper
        % Using ppval (faster than fnval)
        %TrialPt(1:2) = ppval(SplnFt,L_upper);
        
        % But this can be improved:
        % First find the segment, i.e. between what breaks is L_upper
        %inds = discretize(L_upper, [-inf; SplBrk(2:end-1); +inf]);
        inds = find(L_upper >= [-inf; SplBrk(2:end-1); +inf],1,'last');
        % Shift to local coordinates
        x_shf = L_upper - SplBrk(inds);
        % Compute the independent variable values, i.e. powers of x
        zero  = ones(size(x_shf));
        one   = x_shf;
        two   = one .* x_shf;
        three = two .* x_shf;
        % Compute the polynomial value for both sets of coefficients
        TrialPt(1) = sum( [three two one zero] .* cf_x(inds,:), 2);
        TrialPt(2) = sum( [three two one zero] .* cf_y(inds,:), 2);

        TrialPt(3) = F(TrialPt(1),TrialPt(2));
        
        FunVal_upper = norm(TrialPt-BasePt)...
            *Dir*sign(TrialPt(PrgCoor)-BasePt(PrgCoor)) - d;
        
        if FunVal_upper > 0
            break
        elseif isnan(FunVal_upper) || j == MaxIter
            L_upper = SplnFt.breaks(end);
        end
        
    end
    
    % Set interval on spline parameter
    L = [L_prev L_upper];
    
    if norm(SplnFt.breaks(end) - L_prev) < 0.5*sqrt(d)
       EndOfCurveFlag = true;
    end
    
    % Do a max of 1e3 iterations to locate the node
    for j = 1:MaxIter
        
       L_mid = sum(L)/2;
       
       % Evaluate the spline at L_upper
       % Using ppval (faster than fnval)
       %TrialPt(1:2) = ppval(SplnFt,L_mid);
       
       % But this can be improved (from StackOverflow):
       % First find the segment, i.e. between what breaks is L_mid
       %inds = discretize(L_mid, [-inf; SplBrk(2:end-1); +inf]);
       inds = find(L_mid >= [-inf; SplBrk(2:end-1); +inf],1,'last');
       % Shift to local coordinates
       x_shf = L_mid - SplBrk(inds);
       % Compute the independent variable values, i.e. powers of x
       zero  = ones(size(x_shf));
       one   = x_shf;
       two   = one .* x_shf;
       three = two .* x_shf;
       % Compute the polynomial value for both sets of coefficients
       TrialPt(1) = sum( [three two one zero] .* cf_x(inds,:), 2);
       TrialPt(2) = sum( [three two one zero] .* cf_y(inds,:), 2);
       
       TrialPt(3) = F(TrialPt(1),TrialPt(2));
       
       FunVal_mid = norm(TrialPt-BasePt)...
           *Dir*sign(TrialPt(PrgCoor)-BasePt(PrgCoor)) - d;
       
       if isnan(FunVal_mid)
           if Set.MeshOrCurveWarning
               fprintf(2,['Current point on generator curve #%d is outside ' ...
                   'of mold def. \n\n'],i)
           end
           % Remove the origin point and return what is computed
           CurveNode(1,:) = [];
           return
       elseif abs(FunVal_mid) < ConTol
           % Within tolerance: take current point as the node location
           CurveNode(NodeIter,:) = TrialPt;
           L_prev = L_mid;
           break
       elseif FunVal_mid > 0
           % Adjust the interval
           L(2) = L_mid;
       else
           L(1) = L_mid;
       end
    end

    if j == MaxIter && ~EndOfCurveFlag
        if Set.MeshOrCurveWarning
            fprintf(2,'Could not locate node on generator #%d in %d iterations \n\n',...
                i,MaxIter)
        end
        CurveNode(1,:) = [];
        return
    end
    
    nMissingNodes = (nNode + 1)-NodeIter;
    if EndOfCurveFlag && nMissingNodes > 0 && ~Set.ExtrapSteerCurves
        if Set.MeshOrCurveWarning
            fprintf(2,['Reached the end of the defined ' ...
                'curve for generator #%d. \n'...
                'Stopped placing generator nodes with %d remaining \n\n'],...
                i, nMissingNodes)
        end
        CurveNode(1,:) = [];
        return
    end
    
end

CurveNode(1,:) = [];

end

% function [Val,TrialPt] = FunToMin(alpha,SplnFt,BasePt,PrgCoor,d,F,Dir)
% 
% TrialPt(:,1:2) = ppval(SplnFt,alpha)';
% TrialPt(:,3) = F(TrialPt(:,1),TrialPt(:,2));
% 
% Val = abs(vecnorm(TrialPt-BasePt,2,2)...
%     .*Dir.*sign(TrialPt(:,PrgCoor)-BasePt(:,PrgCoor)) - d);
% 
% end

% % Populate spline curve with nodes using a golden section algorithm /
% fsolve 
% This was found to be slower...
% % Initialize and place first node
% CurveNode(1:nNode+1,1:3) = NaN;
% CurveNode(1,1:3) = [Org, F(Org(1),Org(2))];
% 
% delta = 2*d;
% epsilonGS = 1e-5;
% nMaxGSIter = 1e4;
% 
% Opts = optimoptions(@fsolve,'Display','off','Algorithm','levenberg-marquardt');
% 
% for NodeIter = 2:nNode + 1
%     
%     BasePt = CurveNode(NodeIter-1,:);
%     
%     alpha_sol = GoldenSection(...
%        @(alpha)FunToMin(alpha,SplnFt,BasePt,PrgCoor,d,F,Dir), ...
%        delta, epsilonGS, nMaxGSIter);
%     
%     alpha_sol = fsolve(...
%         @(alpha)FunToMin(alpha,SplnFt,BasePt,PrgCoor,d,F,Dir),d,Opts);
%     
%     [~,TrialPt] = FunToMin(alpha_sol,SplnFt,BasePt,PrgCoor,d,F,Dir);
%     
%     CurveNode(NodeIter,:) = TrialPt;
%     
% end