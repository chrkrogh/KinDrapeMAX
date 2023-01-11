function [ObjFun, OptRes, Dra_L] = LayerObj(x,Mold,Plt,Inp,Set,Opt,F,DT,z)
% This function processes the design variables in the layer optimization to
% model input, calls the layer, processes the output and returns the 
% objective for the optimizer

% For sequential layer optimization
if isfield(Opt,'SeqLayer') && Opt.SeqLayer
    % Subtract the previous layer index from the vector of design variables 
    % to get the right indices for the current layer
    Opt.DesVarIdx = structfun(@(x) x-Opt.PrevIdx ,Opt.DesVarIdx ,'UN' ,0);
end

%% Set the input based on design variables

% Set the course width
for i = 1:length(Opt.DesVarIdx.CourseWidth)
    CourseNo = i;
    DesVarNo = Opt.DesVarIdx.CourseWidth(i);
    Inp.Grid_L{CourseNo}(1) = Inp.BndAndCon.AvailableCourseWidths(x(DesVarNo));
end

% Set the longitudinal OrgNode's
for i = 1:length(Opt.DesVarIdx.OrgNodeLong)
    CourseNo = i;
    DesVarNo = Opt.DesVarIdx.OrgNodeLong(i);
    Inp.OrgNode_L{CourseNo}(2) = x(DesVarNo);
    % Convert the 1st component from pct to index.
    [OrgNode_L_temp,~] = ...
        OrgNodeFromPct([],Inp.OrgNodePct_L{CourseNo},Inp.Grid_L{CourseNo});
    Inp.OrgNode_L{CourseNo}(1) = OrgNode_L_temp(1);
    % Move the origin offset from the bottom accordingly
    if strcmpi(Inp.DraContr_L{2},'Steer')
        % Offset the transverse (root) curve
        Inp.SteerPts_L{2,CourseNo}(2,:) = ...
            Inp.SteerPts_L{2,CourseNo}(2,:) + ...
            (Inp.OrgNode_L{CourseNo}(2)-1) * Inp.d_L{CourseNo};
    elseif strcmpi(Inp.DraContr_L{2},'Geo')
        % Move the 2nd entry in Org
        Inp.Org_L{CourseNo}(2) = (Inp.OrgNode_L{CourseNo}(2)-1) * ...
            Inp.d_L{CourseNo};
    end
end

% Set the first course transverse offsets
CourseNo = 1;
if strcmpi(Inp.LayerPropagation,'right-to-left') && ...
        Inp.OrgNodePct_L{CourseNo}(1) == 0
    % If the 1st course has the steering curve along the left edge and the
    % layer propagates from right to left, subtract the course width from 
    % the offset
    CourseWidthOffset = (Inp.Grid_L{CourseNo}(1)-1)*Inp.d_L{CourseNo};
else
    CourseWidthOffset = 0.0;
end
if ~isempty(Opt.DesVarIdx.FirstTransOffset)
    if isa(Inp.DesVar.FirstTransOffset,'logical')
        SteerPtsVec = 1:size(Inp.SteerPts_L{1,CourseNo},2);
    elseif isa(Inp.DesVar.FirstTransOffset,'numeric') && ...
            length(Inp.DesVar.FirstTransOffset) > 1
        SteerPtsVec = Inp.DesVar.FirstTransOffset;
    end
    % Loop over the number of design variables per course
    for j = 1:size(Opt.DesVarIdx.FirstTransOffset,1)
        DesVarNo = Opt.DesVarIdx.FirstTransOffset(j,1);
        OffsetNo = SteerPtsVec(j);
        Inp.SteerPts_L{1,CourseNo}(1,OffsetNo) = x(DesVarNo) - ...
            CourseWidthOffset;
    end
end

% Set the transverse course offsets (gaps/overlaps)
% Loop over the number of courses
for i = 1:size(Opt.DesVarIdx.TransOffset,2)
    CourseNo = i+1;
    if isa(Inp.DesVar.TransOffset,'logical')
        SteerPtsVec = 1:size(Inp.SteerPts_L{1,CourseNo},2);
    elseif isa(Inp.DesVar.TransOffset,'numeric') && ...
            length(Inp.DesVar.TransOffset) > 1
        SteerPtsVec = Inp.DesVar.TransOffset;
    end
    % Loop over the number of design variables per course
    for j = 1:size(Opt.DesVarIdx.TransOffset,1)
        DesVarNo = Opt.DesVarIdx.TransOffset(j,i);
        OffsetNo = SteerPtsVec(j);
        Inp.SteerPts_L{1,CourseNo}(1,OffsetNo) = x(DesVarNo);
    end
end

% Compute the layer
[Dra_L, nCourses_trim] = ComputeLayer(Mold,Plt,Inp,Set,F,DT,z);

%% Gather data for optimization

% Process the shear angles and fiber angle deviations from each course
% Initialize
AngDevCourse = cell(1,nCourses_trim);
ShearCourse = cell(1,nCourses_trim);
for k = 1:nCourses_trim
    % Compute the mean of the shear angles in each cell for the course
    ShearCourse{k} = mean(Dra_L(k).P(:,:,4),2)'; 
    % Compute the mean of the fiber angle dev. in each cell for the course
    AngDevCourse{k} = mean(Dra_L(k).P(:,:,5),2)'; 
end
AllAngDev = [AngDevCourse{:}];
AllAngDev(isnan(AllAngDev)) = [];
AllShearAngles = [ShearCourse{:}];
AllShearAngles(isnan(AllShearAngles)) = [];
TotalTrimArea = sum([Dra_L(:).TrimArea]);

%% Objective function

% Shear angles
if Inp.Obj.Shear && Inp.Obj.p_Shear == inf
    % Max
    Shear_obj = max(abs(AllShearAngles),[],'omitnan');
elseif Inp.Obj.Shear && Inp.Obj.p_Shear < inf
    % p-norm
    p = Inp.Obj.p_Shear;
    Shear_obj = sum(abs(AllShearAngles).^p,'omitnan')^(1/p);
else
    Shear_obj = 0;
end 

% Angle deviations
TargetFiberAng = Inp.Obj.TargetFiberAng;
if Inp.Obj.AngDev && Inp.Obj.p_AngDev == inf
    % Max
    AngDev_obj = max(abs(AllAngDev-TargetFiberAng),[],'omitnan');
elseif Inp.Obj.AngDev && Inp.Obj.p_AngDev < inf
    % p-norm
    p = Inp.Obj.p_AngDev;
    AngDev_obj = sum(abs(AllAngDev-TargetFiberAng).^p,'omitnan')^(1/p);
else
   AngDev_obj = 0; 
end

% Trim area
if Inp.Obj.TrimArea
    % Scale trim area so it can be added to the other criteria
    ShearAng_eq = Inp.Obj.ShearAng_eq;
    TrimArea_eq = Inp.Obj.TrimArea_eq;
    TrimArea_obj = ShearAng_eq * TotalTrimArea/TrimArea_eq;
else
    TrimArea_obj = 0.0;
end

% Total objective function
ObjFun = Shear_obj + AngDev_obj + TrimArea_obj;

%% Store in struct for output
OptRes.TotalTrimArea = TotalTrimArea;
OptRes.AllShearAngles = AllShearAngles;
OptRes.AllAngDev = AllAngDev;
OptRes.nCourses_trim = nCourses_trim;
% Objective function parts
OptRes.Shear_obj = Shear_obj;
OptRes.AngDev_obj = AngDev_obj;
OptRes.TrimArea_obj = TrimArea_obj;
% Design related data
if ~isempty(x) && ~isempty(Opt.DesVarIdx.CourseWidth)
    OptRes.CourseWidths = (Inp.BndAndCon.AvailableCourseWidths(x(...
        Opt.DesVarIdx.CourseWidth(1:nCourses_trim)))-1).*...
        [Inp.d_L{1:nCourses_trim}];
else
    OptRes.CourseWidths = [];
end

if ~isempty(x) && ~isempty(Opt.DesVarIdx.OrgNodeLong)
    OptRes.OrgNodeLong = x(Opt.DesVarIdx.OrgNodeLong(1:nCourses_trim));
else
    OptRes.OrgNodeLong = [];
end

if ~isempty(x) && ~isempty(Opt.DesVarIdx.TransOffset)
    OptRes.TransOffsets = x(Opt.DesVarIdx.TransOffset(:,1:nCourses_trim-1)');
else
    OptRes.TransOffsets = [];
end

if ~isempty(x) && ~isempty(Opt.DesVarIdx.FirstTransOffset)
    OptRes.FirstTransOffsets = x(Opt.DesVarIdx.FirstTransOffset');
else
    OptRes.FirstTransOffsets = [];
end
end